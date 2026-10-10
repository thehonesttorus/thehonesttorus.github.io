"""Exact great-circle leaf compiler, layer-swept and batched (same mathematics as whest/leaf.py, Prop. 2 / Thm. 3).

The circle x(t) = B (cos t, sin t)^T is cut into arcs on which every layer is linear in (cos t, sin t). Instead of
walking events in t-order, sweep the layers: given the arcs and the post-activation coefficients H_{l-1} (n, A, 2) on
each arc, one matrix-matrix product gives Z_l = W_l H_{l-1} on every arc; each neuron's sinusoid a cos t + b sin t has
its roots at psi +- pi/2 (psi = atan2(b, a)); roots inside an arc's open interval are layer-l walls. The refined arcs
take their gate masks from the sign at their midpoints (no accumulated state, so no drift and no tie bookkeeping:
a near-tie only creates a sub-arc of negligible length). H_l = mask * Z_l(parent arc). Cost per layer 4 n^2 A FLOPs.

Returns the leaf mean of every layer's activations, g (L, n) = (1/2 pi) * integral over the circle, and the number
of walls per layer.
"""
import numpy as np

TWO_PI = 2.0 * np.pi


def compile_leaf_b(W, B, dtype=np.float64):
    L = len(W); n_in = B.shape[0]
    bounds = np.array([0.0, TWO_PI])
    H = np.asarray(B, dtype=dtype).reshape(n_in, 1, 2)                     # (n, A, 2), A = 1 arc
    g = np.zeros((L, W[0].shape[0])); walls = np.zeros(L, dtype=np.int64)
    for l in range(L):
        A = H.shape[1]; n = W[l].shape[0]
        Z = (W[l] @ H.reshape(H.shape[0], 2 * A)).reshape(n, A, 2)
        a, b = Z[..., 0], Z[..., 1]
        psi = np.arctan2(b, a)
        lo, hi = bounds[:-1][None, :], bounds[1:][None, :]
        new = []
        for r in (np.mod(psi + 0.5 * np.pi, TWO_PI), np.mod(psi - 0.5 * np.pi, TWO_PI)):
            inside = (r > lo) & (r < hi) & ((a != 0) | (b != 0))
            new.append(r[inside])
        new = np.concatenate(new); walls[l] = new.size
        bounds = np.sort(np.concatenate([bounds, new]))
        tl, tr = bounds[:-1], bounds[1:]; mid = 0.5 * (tl + tr)
        parent = np.searchsorted(bounds_old := lo[0], mid, side="right") - 1      # arcs of the previous partition
        Zp = Z[:, parent, :]                                                       # (n, A', 2)
        mask = (Zp[..., 0] * np.cos(mid)[None, :] + Zp[..., 1] * np.sin(mid)[None, :]) > 0
        H = Zp * mask[..., None]
        ws, wc = np.sin(tr) - np.sin(tl), np.cos(tl) - np.cos(tr)                 # integrals of cos and sin per arc
        g[l] = (H[..., 0] @ ws + H[..., 1] @ wc) / TWO_PI
    return dict(g=g, walls=walls, arcs=bounds.size - 1)
