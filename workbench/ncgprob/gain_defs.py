# Which "gain" is which: the norm-based E2 functional, the plain and q-weighted pair averages, their per-pair structure,
# the Gaussian part with the true marginal state vs GAC's conditional state, GAC's injection on its own state vs on the
# true state, and the last-layer closure mean error (true moments) against the Edgeworth prediction from the diagonal cumulants.
import numpy as np, sys
from scipy.special import ndtr
sys.path.insert(0, "../num12")
from gac import inject, EG
from ledger import gstep
from closure import relu_coeffs
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
D = np.load(f"mc_cum_off{net}.npz"); L, n = D["s1y"].shape
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64); mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
phi = lambda x: np.exp(-0.5 * x * x) / np.sqrt(2 * np.pi)
# --- GAC with the per-layer conditional state kept
def gac_states(Ws):
    out = []; gam = 0.0; g1 = 1.0
    for l, W in enumerate(Ws):
        if l == 0: mu = np.zeros(n); S = W @ W.T; f = 0.0
        else:
            nu = W @ mbar; Q = W @ Sbar @ W.T; mu0 = nu / g1; S0 = Q - np.outer(mu0, mu0)
            f = inject(W, mu0, S0, prev); gam = gam + f; g1 = EG(gam, "gamma"); mu = nu / g1; S = Q - np.outer(mu, mu)
        M, C = gstep(mu, S); mbar = g1 * M; Sbar = C + np.outer(M, M)
        sig = np.sqrt(np.diag(S)); prev = (mu, sig, S / np.outer(sig, sig))
        out.append(dict(mu=mu, S=S, gam=gam, g1=g1, f=f, Q=(S + np.outer(mu, mu)), nu=g1 * mu))
    return out
G = gac_states(list(Wc))
print(f"net {net}: definitions of the gain at z_l (MC T={int(D['T'])})")
print(" l | norm: Var|z|^2/E^2  Gauss(true mu,S) Gauss(GAC cond) | excess(true) excess(GACcond) | pair avg plain  q-weighted  diag part | GAC gamma_l  f_l(GAC state)  f_l(true marginal)  f_l(true, cond-corrected)")
for l in range(L):
    mu = D["s1y"][l]; E2 = D["S2y"][l]; S = E2 - np.outer(mu, mu); q = np.diag(E2); E = q.sum()
    tot = D["nrm"][l, 3] / D["nrm"][l, 2] ** 2 - 1
    gauss_true = (2 * np.sum(S * S) + 4 * mu @ S @ mu) / E ** 2
    mu0, S0 = G[l]["mu"], G[l]["S"]; E0 = np.trace(S0) + mu0 @ mu0
    gauss_gac = (2 * np.sum(S0 * S0) + 4 * mu0 @ S0 @ mu0) / E0 ** 2
    X = D["M22y"][l] - np.outer(q, q) - 2 * S * S - 4 * np.outer(mu, mu) * S
    dg = np.diag(X).sum() / E ** 2
    Xo = X.copy(); np.fill_diagonal(Xo, 0.0); Xn = Xo / np.outer(q, q)
    plain = Xn.sum() / (n * (n - 1)); qw = Xo.sum() / E ** 2
    if l < L - 1:
        W = Wc[l + 1]; muz = D["s1y"][l + 1]; Sz = D["S2y"][l + 1] - np.outer(muz, muz); sig = np.sqrt(np.diag(S)); R = S / np.outer(sig, sig)
        f_true = inject(W, muz, Sz, (mu, sig, R))
        # conditional-corrected: remove the gain-induced rank-one part gamma/4 mu mu^T from S (and from Sz), gamma = plain pair average
        gl = max(plain, 0.0); Sc = S - (gl / 4) * np.outer(mu, mu); sigc = np.sqrt(np.diag(Sc)); Rc = Sc / np.outer(sigc, sigc)
        f_cond = inject(W, muz, Sz - (plain / 4) * np.outer(muz, muz), (mu, sigc, Rc))
        fg = G[l + 1]["f"]
    else: f_true = f_cond = fg = float("nan")
    print(f" {l:2d} | {tot:.5f}  {gauss_true:.5f}  {gauss_gac:.5f} | {tot-gauss_true:.5f}  {tot-gauss_gac:.5f} | {plain:.5f}  {qw:.5f}  {dg:.5f} | {G[l]['gam']:.5f}  {fg:.5f}  {f_true:.5f}  {f_cond:.5f}")
# --- per-pair structure at the last layer and at layer 8
for l in (8, L - 1):
    mu = D["s1y"][l]; E2 = D["S2y"][l]; S = E2 - np.outer(mu, mu); q = np.diag(E2); sig = np.sqrt(np.diag(S)); a = mu / sig
    X = D["M22y"][l] - np.outer(q, q) - 2 * S * S - 4 * np.outer(mu, mu) * S; Xn = X / np.outer(q, q); np.fill_diagonal(Xn, np.nan)
    row = np.nanmean(Xn, 1)
    print(f"\nlayer {l}: per-neuron gain seen by k (row mean of X_kl/(q_k q_l)): mean {row.mean():.5f}, sd {row.std():.5f}; correlations with a_k {np.corrcoef(row, a)[0,1]:+.3f}, mu_k {np.corrcoef(row, mu)[0,1]:+.3f}, q_k {np.corrcoef(row, q)[0,1]:+.3f}, 1/q_k {np.corrcoef(row, 1/q)[0,1]:+.3f}, sigma_k {np.corrcoef(row, sig)[0,1]:+.3f}")
    qq = np.outer(q, q); iu = np.triu_indices(n, 1); xv = Xn[iu]; qv = qq[iu]; mv = np.outer(mu, mu)[iu]; ov = np.outer(a, a)[iu]
    edges = np.quantile(qv, np.linspace(0, 1, 6))
    print("   X_kl by quintile of q_k q_l:", " ".join(f"{np.mean(xv[(qv >= edges[i]) & (qv < edges[i+1] + (i == 4))]):+.5f}" for i in range(5)))
    edges = np.quantile(ov, np.linspace(0, 1, 6))
    print("   X_kl by quintile of a_k a_l:", " ".join(f"{np.mean(xv[(ov >= edges[i]) & (ov < edges[i+1] + (i == 4))]):+.5f}" for i in range(5)))
    # fit X_kl = g0 + g1 * mu_k mu_l /(q_k q_l) + g2 * (sigma_k^2 sigma_l^2)/(q_k q_l): scale-mixture with separate mean and variance channels
    A_ = np.stack([np.ones_like(xv), mv / qv, (np.outer(sig**2, sig**2)[iu]) / qv], 1); coef, *_ = np.linalg.lstsq(A_, xv, rcond=None)
    pred = A_ @ coef; print(f"   fit X = g0 + g1 mu mu/qq + g2 s2 s2/qq: {coef}, R^2 {1 - np.var(xv - pred)/np.var(xv):.3f} (noise floor per pair ~ {3/np.sqrt(D['T']):.4f})")
# --- last-layer closure mean error from the TRUE pre-activation moments vs the Edgeworth prediction from the diagonal cumulants
l = L - 1; mu = D["s1y"][l]; E2 = D["S2y"][l]; S = E2 - np.outer(mu, mu); sig = np.sqrt(np.diag(S)); a = mu / sig
mG = sig * (a * ndtr(a) + phi(a)); err = mG - mt[l]
k3 = D["s3y"][l] - 3 * mu * np.diag(E2) + 2 * mu**3; X4 = D["s4y"][l] - 4 * mu * D["s3y"][l] + 6 * mu * mu * np.diag(E2) - 3 * mu**4; k4 = X4 - 3 * sig**4
A = relu_coeffs(mu, sig, 8)
e3 = -(k3 / 6) * A[3] / sig**3; e4 = -(k4 / 24) * A[4] / sig**4; e33 = -(k3 * k3 / 72) * A[6] / sig**6
scale = lambda v: (mt[l] @ v) / (mt[l] @ mt[l])
print(f"\nlast layer: closure mean from TRUE (mu,sigma) minus truth: MSE {np.mean(err**2):.3e}, scale component c = {scale(err):+.6f} (8c = {8*scale(err):+.5f});")
print(f"   Edgeworth from diagonal cumulants: kappa3 term MSE-explained {1-np.mean((err-e3)**2)/np.mean(err**2):.3f}, +kappa4 {1-np.mean((err-e3-e4)**2)/np.mean(err**2):.3f}, +kappa3^2 {1-np.mean((err-e3-e4-e33)**2)/np.mean(err**2):.3f}; scale components 8c: k3 {8*scale(e3):+.5f}, k4 {8*scale(e4):+.5f}, k3^2 {8*scale(e33):+.5f}")
print(f"   GAC: gamma_L {G[l]['gam']:.5f} -> 8(1-E[G]) = {8*(1-G[l]['g1']):.5f}; GAC final MSE {np.mean((G[l]['g1']*gstep(G[l]['mu'], G[l]['S'])[0]-mt[l])**2):.3e}")
# per-neuron kappa3, kappa4 of z_L: coherent (common) parts vs scale-mixture forms
print(f"   kappa4/q^2: mean {np.mean(k4/np.diag(E2)**2):+.5f}; scale-mixture 3gamma(1+2a^2)/(1+a^2)^2... fit kappa4 = g*(3 s^4 + 6 mu^2 s^2): g = {np.sum(k4*(3*sig**4+6*mu**2*sig**2))/np.sum((3*sig**4+6*mu**2*sig**2)**2):.5f}")
print(f"   kappa3: fit kappa3 = g*1.5 mu s^2: g = {np.sum(k3*1.5*mu*sig**2)/np.sum((1.5*mu*sig**2)**2):.5f}; corr(kappa3, mu s^2) {np.corrcoef(k3, mu*sig**2)[0,1]:+.3f}; mean kappa3/(mu s^2) {np.mean(k3/(mu*sig**2)):+.5f}")
