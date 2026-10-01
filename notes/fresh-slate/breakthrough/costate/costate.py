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
from hd import Gauss, inject, inject_full2  # validated Gaussian quantities and Stein injection (heisenberg stream)
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
    S = 0.0 if V is None else 2 * (U * V).T @ (a2[:, None] * V) + (V * V).T @ (a2[:, None] * U)
    if coinc:
        S += (U * U).T @ (X + dl[:, None] * U) + 2 * (U * X).T @ U
    return S


def Kh(C, h):
    """(2,1) slice of the scale-field mode z = (1+zeta) x, Cov(zeta, x) = h:  2 C_pp h_q + 4 h_p C_pq."""
    return 2 * np.diag(C)[:, None] * h[None, :] + 4 * h[:, None] * C


def Kh_adj(C, X):
    return 2 * np.diag(C) @ X + 4 * np.sum(C * X, axis=1)


def fit_h(C, S, iters=40):
    """Least-squares h for S ~ Kh(C, h), conjugate gradient on the normal equations (O(n^2) per iteration)."""
    b = Kh_adj(C, S); h = np.zeros_like(b); r = b.copy(); p = r.copy(); rr = r @ r
    for _ in range(iters):
        Ap = Kh_adj(C, Kh(C, p)); a = rr / (p @ Ap); h += a * p; r -= a * Ap
        rn = r @ r
        if rn < 1e-24 * (b @ b): break
        p = r + (rn / rr) * p; rr = rn
    return h


def filt(Sold, how, r=None):
    """Oracle filters on the exact old slice: what part of it does the readout need?"""
    if how == "diag":
        return np.diag(np.diag(Sold))
    if how == "off":
        return Sold - np.diag(np.diag(Sold))
    if how == "rank":
        Uu, s, Vt = np.linalg.svd(Sold)
        return (Uu[:, :r] * s[:r]) @ Vt[:r] + np.diag(np.diag(Sold) - np.einsum("ij,j,ji->i", Uu[:, :r], s[:r], Vt[:r]))
    if callable(how):
        return how(Sold)
    raise ValueError(how)


def predict(Ws, A=None, old="drop", rank=None, coinc=True, K=8, record=None, Ap=0, oldfilter=None, r=None, law=False, keep_old=False, k4=None):
    """old: 'drop' | 'slice' (diagonal double edge only) | 'pool' (renewal: aged-out content is projected on its
    (2,1) slice at the current layer and re-emitted as a secondary form, transported exactly for Ap more layers
    (Ap=None: never re-projected); Ap=0 is the full coincident-support slice chain)."""
    Ws = np.asarray(Ws, dtype=np.float64)
    L, n, _ = Ws.shape
    m = np.zeros(n); C = Ws[0].T @ Ws[0]
    live = []          # forms: dict(s, kind, a2, dl, U, V, X)
    Sold = None
    gacc = 0.0         # old == 'gsm': accumulated scale-mode amplitude (conserved scalar)
    out = []
    nprod = 0
    for l in range(L):
        v = 0.5 * gacc if law else 0.0           # scale variance: slice amplitude gamma = 2 Var(t)
        if law and v != 0:
            # model constraint Var z_p = (1+v) C_x,pp + v m_p^2 >= v m_p^2 for every p (and |v| < 1)
            vmax = 0.95 * np.min(np.diag(C) / np.maximum(m * m, 1e-300))
            if v > vmax or v < -0.5:
                v = min(max(v, -0.5), vmax)
                if record is not None:
                    record.append(dict(layer=l, clamp=True))
            C = (C - v * np.outer(m, m)) / (1 + v)
        G = Gauss(m, C, K=K)
        S = np.zeros((n, n))
        sl = []
        Sfilt = np.zeros((n, n)) if oldfilter is not None else None
        for f in live:
            Sf = pair_slice(f["a2"], f["U"], f["V"], f["X"], f["dl"], coinc or f["kind"] == "pool")
            sl.append(Sf)
            if oldfilter is not None and A is not None and (l - f["s"] - 1) > A:
                Sfilt += Sf
            else:
                S += Sf
            nprod += (2 if f["V"] is not None else 0) + (2 if (coinc or f["kind"] == "pool") else 0)
        if Sold is not None:
            S += Sold
        Kg = 2 * G.m[:, None] * G.C + (G.m[None, :] * np.diag(G.C)[:, None])
        if old in ("gsm", "gsmslice") and gacc != 0.0 and not law:
            S += gacc * Kg
        if record is not None and Sfilt is not None and np.any(Sfilt):
            record.append(dict(layer=l, gsm_oracle=float(np.sum(Kg * Sfilt) / np.sum(Kg * Kg)), gsm_pred=float(gacc),
                               Sold=Sfilt.copy() if keep_old else None, Kg=Kg.copy() if keep_old else None,
                               vecs=dict(one=np.ones(n), m=G.m.copy(), s=G.s.copy(), s2=G.s ** 2, Ea=G.Ea.copy(),
                                         Phi=G.Phi.copy(), p0=G.Ed[0].copy()) if keep_old else None))
        if Sfilt is not None and live and np.any(Sfilt):
            if oldfilter in ("hfit", "hfitd"):
                h = fit_h(G.C, Sfilt)
                Kf = Kh(G.C, h)
                if oldfilter == "hfitd":
                    Kf = Kf - np.diag(np.diag(Kf)) + np.diag(np.diag(Sfilt))
                S += Kf
            elif isinstance(oldfilter, str) and oldfilter.startswith("gsm"):
                # Gaussian scale mixture z = (1+d) x: slice of 2 Var(d) (m (x) C)_sym, one scalar per layer fitted (oracle)
                Kg = 2 * G.m[:, None] * G.C + (G.m[None, :] * np.diag(G.C)[:, None])
                if oldfilter == "gsm_off":
                    Ko = Kg - np.diag(np.diag(Kg)); So = Sfilt - np.diag(np.diag(Sfilt))
                    gam = np.sum(Ko * So) / np.sum(Ko * Ko)
                    S += gam * Ko + np.diag(np.diag(Sfilt))
                else:
                    gam = np.sum(Kg * Sfilt) / np.sum(Kg * Kg)
                    S += gam * Kg
                if record is not None:
                    record.append(dict(layer=l, gsm_gamma=gam))
            else:
                S += filt(Sfilt, oldfilter, r)
        if record is not None:
            record.append(dict(layer=l, S=S.copy()))
        if live or Sold is not None or (gacc != 0.0 and not law):
            D = np.diag(S).copy()
            if k4 is not None and l in k4:          # oracle hook: (K4d, K22, K31) of z_l, second-order injection
                dEa, dC = inject_full2(G, D, S, *k4[l])
            else:
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
        if old in ("slice", "gsmslice"):
            newold = None
            if Sold is not None:
                newold = (M * M).T @ Sold @ M; nprod += 2
            for f, Sf in zip(live, sl):
                if A is not None and (l - f["s"]) > A:
                    if old == "gsmslice":
                        gk = float(np.sum(Kg * Sf) / np.sum(Kg * Kg)); gacc += gk
                        Sf = Sf - gk * Kg
                    t = (M * M).T @ Sf @ M; nprod += 2
                    newold = t if newold is None else newold + t
            Sold = newold
        nl = []
        pool = None
        for f, Sf in zip(live, sl):
            age_next = (l + 1) - f["s"] - 1
            lim = (None if oldfilter is not None else A) if f["kind"] == "src" else Ap
            if lim is not None and age_next > lim:
                if old == "pool":
                    pool = Sf if pool is None else pool + Sf
                if old == "gsm":
                    gacc += float(np.sum(Kg * Sf) / np.sum(Kg * Kg))
                continue
            for key in ("U", "V", "X"):
                if f[key] is None:
                    continue
                f[key] = f[key] @ M; nprod += 1
            if rank is not None:
                d = rank(age_next)
                if d is not None and d < n:
                    _, _, Vt = np.linalg.svd(f["U"], full_matrices=False)
                    P = Vt[:d].T
                    for key in ("U", "V", "X"):
                        if f[key] is not None:
                            f[key] = (f[key] @ P) @ P.T
            nl.append(f)
        live = nl
        if pool is not None:
            # secondary form at layer l: coincident-support tensor with slice `pool`, seen from z_{l+1} through M
            Dl = pool.copy(); dl = np.diag(pool).copy(); np.fill_diagonal(Dl, 0.0)
            live.append(dict(s=l, kind="pool", a2=None, dl=dl, U=M.copy(), V=None, X=Dl @ M)); nprod += 1
        a2, R, Dl, dl = source_atoms(G, K)
        if A is None or A >= 0:
            live.append(dict(s=l, kind="src", a2=a2, dl=dl, U=W.copy(), V=R @ W, X=(Dl @ W) if coinc else None))
            nprod += 2 if coinc else 1
        m = Ea @ W
        Ca = G.cov_a() + dC
        if law and v != 0:
            Ca = (1 + v) * Ca + v * np.outer(Ea, Ea)
        C = W.T @ Ca @ W
        nprod += 2
    if record is not None:
        record.append(dict(nprod=nprod))
    return np.array(out)
