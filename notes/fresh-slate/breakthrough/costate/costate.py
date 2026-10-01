"""Minimal co-state stream: first-order Heisenberg-Duhamel correction evaluated in pull-back (source-atom) form.

No n^3 tensors: every source kappa_3 tensor is kept as its structured trilinear form (star diagrams +
exact coincident-index corrections) and evaluated on transported directions U_{s->k}, so the cost is a
handful of n^3 products per live (source s, target k) pair. Runs directly at n = 1024.

Closure switches (the co-state candidates of REPORT.md):
  A      : max age carried exactly (None = all). age = k - s - 1.
  old    : what replaces sources older than A: 'drop' | 'slice' (n^2 (2,1)-slice chain: annealed double edge)
  rank   : None or a function age -> d; transported factors of that age are projected on the top-d
           right-singular subspace of U_{s->k} (Tucker truncation of old content).
  coinc  : include the coincident-index correction atoms (Delta, delta).
"""
import sys, os, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../designs/heisenberg"))
from hd import Gauss, inject  # validated Gaussian quantities and Stein injection (heisenberg stream)
from math import factorial


def source_atoms(G, K=8):
    """Structured third cumulant of a = ReLU(z), z ~ N(m, C) (centred), as a trilinear form:
    T(x,y,z) = sum_c a2_c [x_c (Ry)_c (Rz)_c + perms]            (star, all index values)
             + sum_{p!=r} Dl_pr [x_p y_p z_r + x_p z_p y_r + y_p z_p x_r]   (exact (2,1) coincidences)
             + sum_p dl_p x_p y_p z_p                                (exact triple diagonal)."""
    al, rho = G.alpha, G.rho
    a1, a2 = al[1], al[2]
    R = rho * a1[None, :]
    P = np.zeros_like(rho)
    rk = np.ones_like(rho)
    for k in range(1, K + 1):
        rk = rk * rho
        P += G.g2[k][:, None] * al[k][None, :] * rk / factorial(k)
    star_ppr = 2 * (a1 * a2)[:, None] * a1[None, :] * rho + (a1 ** 2)[:, None] * a2[None, :] * rho ** 2
    Dl = P - star_ppr
    np.fill_diagonal(Dl, 0.0)
    dl = G.k3_a - 3 * a1 ** 2 * a2
    return a2, R, Dl, dl


def pair_slice(a2, U, V, X, dl, coinc=True):
    """S(p,q) = T(u_p, u_p, u_q) for the source with factors U (=U_{s->k}), V = R U, X = Dl U."""
    S = 2 * (U * V).T @ (a2[:, None] * V) + (V * V).T @ (a2[:, None] * U)
    if coinc:
        S += (U * U).T @ (X + dl[:, None] * U) + 2 * (U * X).T @ U
    return S


def predict(Ws, A=None, old="drop", rank=None, coinc=True, K=8, record=None, dtype=np.float64):
    Ws = np.asarray(Ws, dtype=np.float64)
    L, n, _ = Ws.shape
    m = np.zeros(n); C = Ws[0].T @ Ws[0]
    live = []          # list of dict(s, a2, dl, U, V, X)
    Sold = None        # n^2 slice chain for sources older than A (old == 'slice')
    out = []
    nprod = 0
    for l in range(L):
        G = Gauss(m, C, K=K)
        S = np.zeros((n, n))
        for src in live:
            S += pair_slice(src["a2"], src["U"], src["V"], src["X"], src["dl"], coinc)
            nprod += 4 if coinc else 2
        if Sold is not None:
            S += Sold
        if record is not None:
            record.append(dict(layer=l, S=S.copy()))
        if live or Sold is not None:
            D = np.diag(S).copy()
            dEa, dC = inject(G, D, S)
        else:
            dEa, dC = np.zeros(n), np.zeros((n, n))
        Ea = G.Ea + dEa
        out.append(Ea)
        if l + 1 == L:
            break
        W = Ws[l + 1]
        g = G.Phi
        M = g[:, None] * W                        # one linear-response step: diag(Phi_l) W_{l+1}
        # old content handled by an n^2 closure: the (2,1)-slice chain with diagonal (annealed) double edge
        if old == "slice":
            newold = None
            if Sold is not None:
                newold = (M * M).T @ Sold @ M; nprod += 2
            for src in [s for s in live if A is not None and (l - s["s"]) > A]:
                Sa = pair_slice(src["a2"], src["U"], src["V"], src["X"], src["dl"], coinc)
                t = (M * M).T @ Sa @ M; nprod += 2
                newold = t if newold is None else newold + t
            Sold = newold
        # advance live sources; retire those older than A
        nl = []
        for src in live:
            age_next = (l + 1) - src["s"] - 1
            if A is not None and age_next > A:
                continue
            for key in ("U", "V", "X"):
                if key == "X" and not coinc:
                    continue
                src[key] = src[key] @ M; nprod += 1
            if rank is not None:
                d = rank(age_next)
                if d is not None and d < n:
                    # project transported factors on the top-d right singular subspace of U_{s->k+1}
                    _, _, Vt = np.linalg.svd(src["U"], full_matrices=False)
                    P = Vt[:d].T
                    for key in ("U", "V", "X"):
                        if key in src and src[key] is not None:
                            src[key] = (src[key] @ P) @ P.T
            nl.append(src)
        live = nl
        # new source at layer l (feeds z_{l+1})
        a2, R, Dl, dl = source_atoms(G, K)
        src = dict(s=l, a2=a2, dl=dl, U=W.copy(), V=R @ W, X=(Dl @ W) if coinc else None)
        nprod += 2 if coinc else 1
        if A is None or A >= 0:
            live.append(src)
        m = Ea @ W
        C = W.T @ (G.cov_a() + dC) @ W
        nprod += 2
    if record is not None:
        record.append(dict(nprod=nprod))
    return np.array(out)
