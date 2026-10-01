"""First-order Heisenberg-Duhamel estimator in its n x n source-layer pull-back form (DESIGN.md section 2-3).
No n^3 tensors: for every live (source l, target k) pair the transported direction matrix
U = W_{l+1} Phi_{l+1} W_{l+2} ... Phi_{k-1} W_k is kept, and the diagonal kappa_3 and the (2,1) slice at layer k
are evaluated on its columns from the Gaussian source data of layer l (diag + star diagrams + exact coincident planes).
Equals hd(..., diagrams="star") to rounding."""
import numpy as np
from hd import Gauss, inject


def source_data(G, planes=True):
    a1, a2 = G.alpha[1], G.alpha[2]
    rho = G.rho
    n = len(a1)
    # exact coincident plane P[p,r] = E[ap~^2 ar~] (Hermite series in rho_pr), as in Gauss.source_k3
    P = np.zeros((n, n)); rk = np.ones_like(rho); f = 1.0
    for k in range(1, G.K + 1):
        rk = rk * rho; f *= k
        P += G.g2[k][:, None] * G.alpha[k][None, :] * rk / f
    star_ppr = 2 * (a2 * a1)[:, None] * a1[None, :] * rho + (a1 ** 2)[:, None] * a2[None, :] * rho ** 2
    D2 = P - star_ppr
    np.fill_diagonal(D2, 0.0)
    if not planes:
        D2 = None
    D1 = G.k3_a - 3 * a2 * a1 ** 2
    return dict(rho=rho, a1=a1, a2=a2, D2=D2, D1=D1)


def rebirth_data(G, D, S):
    """Exact first-order response of the coincident patterns of kappa_3(z) through ReLU, minus the mean-gate
    response that the transported directions already apply: (2,1) entries Phi_q[Phi_p - E(a_p) phi_p/s_p] - Phi_p^2 Phi_q,
    diagonal Phi + t phi e2/2 - e1 phi - e1^2 t phi - Phi^3 (zero-correlation forms; judge theory.md section 2)."""
    Phi, ph, t, s = G.Phi, G.ph, G.t, G.s
    e1 = G.Ea / s
    e2 = (G.var_a + G.Ea ** 2) / s ** 2
    c2 = Phi - G.Ea * ph / s
    c3 = Phi + 0.5 * t * ph * e2 - e1 * ph - e1 ** 2 * t * ph
    D2 = ((c2 - Phi ** 2)[:, None] * Phi[None, :]) * S
    np.fill_diagonal(D2, 0.0)
    D1 = (c3 - Phi ** 3) * D
    return dict(rho=None, D2=D2, D1=D1)


def eval_DS(src, U):
    """D[p] = k3[U_p,U_p,U_p], S[p,q] = k3[U_p,U_p,U_q] for the Gaussian-generated k3 of the source."""
    U2 = U * U
    if src['rho'] is not None:
        a1, a2 = src['a1'][:, None], src['a2'][:, None]
        R = src['rho'] @ (a1 * U)
        S = 2 * ((a2 * U * R).T @ R) + (a2 * R * R).T @ U      # star
    else:
        S = np.zeros((U.shape[1], U.shape[1]))
    if src['D2'] is not None:
        D2U = src['D2'] @ U
        S += U2.T @ D2U + 2 * ((U * D2U).T @ U)                # coincident planes
    S += (src['D1'][:, None] * U2).T @ U                       # exact diagonal
    return np.diag(S).copy(), S


def hd_pull(Ws, A=None, planes=True, rebirth=False, trace=None):
    L, n, _ = Ws.shape
    m = np.zeros(n); C = Ws[0].T @ Ws[0]
    live = []          # list of [src, U]
    out = []
    for l in range(L):
        G = Gauss(m, C, K=8)
        if live:
            D = np.zeros(n); S = np.zeros((n, n))
            for src, U in live:
                d, s = eval_DS(src, U); D += d; S += s
            dEa, dC = inject(G, D, S)
            if trace is not None:
                trace.append((l, D, S))
        else:
            dEa, dC = np.zeros(n), np.zeros((n, n))
        Ea = G.Ea + dEa; out.append(Ea)
        if l + 1 < L:
            W = Ws[l + 1]
            gW = G.Phi[:, None] * W
            live = [[src, U @ gW] for src, U in live]
            if A is not None:
                live = live[-A:] if A > 0 else []
            if rebirth and l > 0:
                live.append([rebirth_data(G, D, S), W.copy()])
            live.append([source_data(G, planes), W.copy()])
            Ca = G.cov_a() + dC
            m = Ea @ W; C = W.T @ Ca @ W
    return np.array(out)
