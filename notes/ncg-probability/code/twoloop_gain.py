# Two-loop gain injection: first variation of the coherent pair-sum gain functional under the third (and fourth)
# cumulants of the pre-activation, evaluated on the true Monte Carlo pair cumulants of official network NET.
#   Gamma_l  = (1/n(n-1)) sum_{k!=l} [Cov(z_k^2, z_l^2) - Gaussian(mu,S)] / (q_k q_l)   (pair-average gain at z_l)
#   coherent restriction with z = W h:  Gamma = [sum_{a!=c} C4_ac (u_a u_c - U2_ac) + sum_a k4_a (u_a^2 - v_a)
#                                                + 4 sum_{a,c} C3_ac (u_a t_c - V_ac)] / (n(n-1))
#   one loop: C3, C4 of h from Gaussian y (Mehler);  two loop: first variation under kappa3(y) (Edgeworth, exact
#   Gaussian integration by parts: dE[F] = 1/6 sum kappa_ijk E[d_ijk F]);  kappa4(y) variation = transport test.
import numpy as np, sys, time
from math import factorial
from scipy.special import ndtr
sys.path.insert(0, "../num12")
from closure import relu_coeffs, relu2_coeffs
from gac import gac, inject
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
K = 14; KH = K + 5
D = np.load(f"mc_cum_off{net}.npz"); L, n = D["s1y"].shape
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64)
phi = lambda x: np.exp(-0.5 * x * x) / np.sqrt(2 * np.pi)
def mehler(P, Q, R, al, be):
    out = np.zeros((n, n)); Rj = np.ones((n, n))
    for j in range(K + 1):
        if j > 0: Rj = Rj * R
        out += np.outer(P[j + al], Q[j + be]) * (Rj / factorial(j))
    return out
def cum_pair(m, e2, e11, e21, e12, e22, e3, e4):
    """pair cumulants of h from raw moments (diagonals of e11,e21,e12,e22 must be e2,e3,e3,e4)"""
    C3 = e21 - e2[:, None] * m[None, :] - 2 * m[:, None] * e11 + 2 * (m * m)[:, None] * m[None, :]
    v = e2 - m * m; cov = e11 - np.outer(m, m)
    X2Y2 = (e22 - 2 * m[None, :] * e21 - 2 * m[:, None] * e12 + 4 * np.outer(m, m) * e11
            + (m * m)[None, :] * e2[:, None] + (m * m)[:, None] * e2[None, :] - 3 * np.outer(m * m, m * m))
    C4 = X2Y2 - np.outer(v, v) - 2 * cov * cov
    k4 = np.diag(C4).copy(); np.fill_diagonal(C4, 0.0)
    return C3, C4, k4
def gamma_coh(C3, C4off, k4, W, q, muz):
    P = W * W; Pq = P / q[:, None]
    u = Pq.sum(0); v = (Pq * Pq).sum(0); t = W.T @ (muz / q)
    U2 = Pq.T @ Pq; V = Pq.T @ (W * (muz / q)[:, None])
    pair = np.sum(C4off * (np.outer(u, u) - U2)); diag = np.sum(k4 * (u * u - v)); mean = 4 * np.sum(C3 * (np.outer(u, t) - V))
    return (pair + diag + mean) / (n * (n - 1)), (pair / (n * (n - 1)), diag / (n * (n - 1)), mean / (n * (n - 1)))
def gamma_meas(m):
    mu = D["s1y"][m]; E2 = D["S2y"][m]; S = E2 - np.outer(mu, mu); q = np.diag(E2)
    X = (D["M22y"][m] - np.outer(q, q) - 2 * S * S - 4 * np.outer(mu, mu) * S)
    Xn = X / np.outer(q, q); np.fill_diagonal(Xn, 0.0); Xq = X.copy(); np.fill_diagonal(Xq, 0.0)
    return Xn.sum() / (n * (n - 1)), Xq.sum() / q.sum() ** 2
g = gac(list(Wc)); gam_gac = np.array([d["gamma"] for d in g])
Gm = np.array([gamma_meas(m) for m in range(L)])
print(f"net {net}, T = {int(D['T'])}: measured pair-average gain Gamma_l at z_l (plain / q-weighted) vs GAC gamma_l")
rows = []
t0 = time.time()
for l in range(L - 1):
    mu = D["s1y"][l]; E2 = D["S2y"][l]; S = E2 - np.outer(mu, mu); e2y = np.diag(E2); sig = np.sqrt(np.diag(S)); R = S / np.outer(sig, sig)
    a = mu / sig; Pa = ndtr(a); fa = phi(a)
    W = Wc[l + 1]; q = np.diag(D["S2y"][l + 1]); muz = D["s1y"][l + 1]
    # ---- true cumulants of h (coherent restriction test)
    m = D["s1h"][l]; e2 = np.diag(D["S2h"][l]); e11 = D["S2h"][l]; e21 = D["M21h"][l]; e12 = e21.T.copy(); e22 = D["M22h"][l]; e3 = D["s3h"][l]; e4 = D["s4h"][l]
    C3t, C4t, k4t = cum_pair(m, e2, e11, e21, e12, e22, e3, e4)
    G_true, parts_true = gamma_coh(C3t, C4t, k4t, W, q, muz)
    # ---- Gaussian reference (one loop) on the true (mu, S)
    A = relu_coeffs(mu, sig, KH); B = relu2_coeffs(mu, sig, KH)
    mG = A[0]; e2G = B[0]
    e3G = sig**3 * ((a**3 + 3 * a) * Pa + (a * a + 2) * fa); e4G = sig**4 * ((a**4 + 6 * a * a + 3) * Pa + (a**3 + 5 * a) * fa)
    def raw(al, be, dA=None, dB=None):
        pass
    e11G = mehler(A, A, R, 0, 0); e21G = mehler(B, A, R, 0, 0); e12G = e21G.T.copy(); e22G = mehler(B, B, R, 0, 0)
    for M_, dg in ((e11G, e2G), (e21G, e3G), (e12G, e3G), (e22G, e4G)): np.fill_diagonal(M_, dg)
    C3G, C4G, k4G = cum_pair(mG, e2G, e11G, e21G, e12G, e22G, e3G, e4G)
    f1, parts1 = gamma_coh(C3G, C4G, k4G, W, q, muz)
    # GAC's own inject on the same data (K=14, no U2/V corrections) for reference
    f1_gac = inject(W, muz, D["S2y"][l + 1] - np.outer(muz, muz), (mu, sig, R))
    # ---- kappa3 of y: total and tree (minus the gain-induced scale-mixture part)
    k3y = D["s3y"][l] - 3 * mu * e2y + 2 * mu**3
    K21y = D["M21y"][l] - mu[None, :] * e2y[:, None] - 2 * mu[:, None] * E2 + 2 * (mu * mu)[:, None] * mu[None, :]
    gl = Gm[l, 0]
    k3g = 1.5 * gl * mu * sig**2; K21g = gl * (mu[:, None] * S + 0.5 * (sig**2)[:, None] * mu[None, :])
    def two_loop(k3, K21):
        k3s = k3 / sig**3; K21s = K21 / np.outer(sig**2, sig); K12s = K21s.T
        dm = k3s * A[3] / 6; de2 = k3s * B[3] / 6; de3 = k3s * sig**3 * Pa; de4 = 4 * k3s * sig**4 * (a * Pa + fa)
        def dE(P, Q):
            return (k3s[:, None] * mehler(P, Q, R, 3, 0) + 3 * K21s * mehler(P, Q, R, 2, 1)
                    + 3 * K12s * mehler(P, Q, R, 1, 2) + k3s[None, :] * mehler(P, Q, R, 0, 3)) / 6
        d11 = dE(A, A); d21 = dE(B, A); d12 = dE(A, B); d22 = dE(B, B)
        for M_, dg in ((d11, de2), (d21, de3), (d12, de3), (d22, de4)): np.fill_diagonal(M_, dg)
        C3p, C4p, k4p = cum_pair(mG + dm, e2G + de2, e11G + d11, e21G + d21, e12G + d12, e22G + d22, e3G + de3, e4G + de4)
        return gamma_coh(C3p - C3G, C4p - C4G, k4p - k4G, W, q, muz), np.abs(d12 - d21.T).max() / (np.abs(d21).max() + 1e-300)
    (f2_tot, p2_tot), asym = two_loop(k3y, K21y)
    (f2_tree, p2_tree), _ = two_loop(k3y - k3g, K21y - K21g)
    (f2_gain, _), _ = two_loop(k3g, K21g)
    # ---- kappa4 of y first variation (transport of the carried gain + non-gain kappa4)
    X4 = D["s4y"][l] - 4 * mu * D["s3y"][l] + 6 * mu * mu * e2y - 3 * mu**4
    X31 = (D["M31y"][l] - mu[None, :] * D["s3y"][l][:, None] - 3 * mu[:, None] * D["M21y"][l] + 3 * np.outer(mu, mu) * e2y[:, None]
           + 3 * (mu * mu)[:, None] * E2 - 3 * np.outer(mu**3, mu))
    X22 = (D["M22y"][l] - 2 * mu[None, :] * D["M21y"][l] - 2 * mu[:, None] * D["M21y"][l].T + 4 * np.outer(mu, mu) * E2
           + (mu * mu)[None, :] * e2y[:, None] + (mu * mu)[:, None] * e2y[None, :] - 3 * np.outer(mu * mu, mu * mu))
    k4y = X4 - 3 * sig**4; K31y = X31 - 3 * (sig**2)[:, None] * S; K22y = X22 - np.outer(sig**2, sig**2) - 2 * S * S
    k4s = k4y / sig**4; K31s = K31y / np.outer(sig**3, sig); K22s = K22y / np.outer(sig**2, sig**2); K13s = K31s.T
    dm = k4s * A[4] / 24; de2 = k4s * B[4] / 24; de3 = k4s * sig**3 * fa / 4; de4 = k4s * sig**4 * Pa
    def dE4(P, Q):
        return (k4s[:, None] * mehler(P, Q, R, 4, 0) + 4 * K31s * mehler(P, Q, R, 3, 1) + 6 * K22s * mehler(P, Q, R, 2, 2)
                + 4 * K13s * mehler(P, Q, R, 1, 3) + k4s[None, :] * mehler(P, Q, R, 0, 4)) / 24
    d11 = dE4(A, A); d21 = dE4(B, A); d12 = dE4(A, B); d22 = dE4(B, B)
    for M_, dg in ((d11, de2), (d21, de3), (d12, de3), (d22, de4)): np.fill_diagonal(M_, dg)
    C3p, C4p, k4p = cum_pair(mG + dm, e2G + de2, e11G + d11, e21G + d21, e12G + d12, e22G + d22, e3G + de3, e4G + de4)
    f3, p3 = gamma_coh(C3p - C3G, C4p - C4G, k4p - k4G, W, q, muz)
    rows.append(dict(l=l, G=Gm[l, 0], Gq=Gm[l, 1], Gn=Gm[l + 1, 0], Gtrue=G_true, f1=f1, f1gac=f1_gac, f2=f2_tot, f2tree=f2_tree, f2gain=f2_gain,
                     f3=f3, p1=parts1, p2=p2_tot, ptrue=parts_true, asym=asym))
    print(f"  layer {l:2d} done ({time.time()-t0:.0f}s)", flush=True)
print()
print("A. one-step ledger at z_{l+1}: measured Gamma_{l+1} | coherent sum from TRUE h cumulants | one loop f1 (Gaussian y, true mu,S) [GAC inject] | kappa3 variation f2 (total / tree part / gain-induced part) | kappa4 variation f3 | Gamma_l")
print("  l  | Gamma_{l+1}  | coh(true h) | f1        [inject]   |  f2 total   f2 tree    f2 gain-ind |  f3        | Gamma_l    | f1+f2+f3 | (1+G_l)(1+f1)-1")
for r in rows:
    print(f"  {r['l']:2d} | {r['Gn']:.5f}     | {r['Gtrue']:.5f}     | {r['f1']:.5f}  [{r['f1gac']:.5f}] | {r['f2']:+.5f}  {r['f2tree']:+.5f}  {r['f2gain']:+.5f} | {r['f3']:+.5f}   | {r['G']:.5f}   | {r['f1']+r['f2']+r['f3']:.5f}  | {(1+r['G'])*(1+r['f1'])-1:.5f}")
print()
print("B. critical accumulation vs measured: gamma^hat_l = sum_{j<l} (f1_j [+ f2_j]) against Gamma_l (plain pair average) and GAC's gamma_l")
c1 = np.concatenate([[0], np.cumsum([r["f1"] for r in rows])]); c12 = np.concatenate([[0], np.cumsum([r["f1"] + r["f2tree"] for r in rows])]); c12t = np.concatenate([[0], np.cumsum([r["f1"] + r["f2"] for r in rows])])
print("  l  | Gamma_l meas | q-weighted | GAC gamma_l | sum f1   | sum f1+f2(tree) | sum f1+f2(total) | sum f1+f2 / Gamma")
for l in range(L):
    print(f"  {l:2d} | {Gm[l,0]:.5f}      | {Gm[l,1]:.5f}    | {gam_gac[l]:.5f}     | {c1[l]:.5f}  | {c12[l]:.5f}         | {c12t[l]:.5f}          | {c12[l]/max(Gm[l,0],1e-12):.3f}")
print()
print("C. parts of the one-loop and two-loop terms (pair / diagonal / mean-coupled), and the pair-moment symmetry check of the variation")
for r in rows:
    print(f"  {r['l']:2d} | f1: {r['p1'][0]:+.5f} {r['p1'][1]:+.5f} {r['p1'][2]:+.5f} | f2: {r['p2'][0]:+.5f} {r['p2'][1]:+.5f} {r['p2'][2]:+.5f} | true-h: {r['ptrue'][0]:+.5f} {r['ptrue'][1]:+.5f} {r['ptrue'][2]:+.5f} | asym {r['asym']:.1e}")
np.savez(f"twoloop_off{net}.npz", rows=np.array([(r['l'], r['G'], r['Gq'], r['Gn'], r['Gtrue'], r['f1'], r['f1gac'], r['f2'], r['f2tree'], r['f2gain'], r['f3']) for r in rows]), gam_gac=gam_gac, Gm=Gm)
