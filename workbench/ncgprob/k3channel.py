# The mean-coupled third-cumulant channel: kappa3(z_k) = 1.5 g mu_k sigma_k^2.  Measured g_l per layer (and the pair slice
# D21 = g (mu_a S_ac + 1/2 sigma_a^2 mu_c)), GAC's gamma_l, the chain's own D3 coefficient, and the one-step generation ledger
# of the coherent kappa3 of z_{l+1} from the pair cumulants of h_l: Gaussian (one loop) + kappa3(y) variation (transport of
# the channel) + kappa4(y) variation.
import numpy as np, sys, pickle, time
sys.path.insert(0, "../num12")
from closure import relu_coeffs, relu2_coeffs
from gac import gac
from pairvar import Layer, y_cumulants, h_cumulants
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
D = np.load(f"mc_cum_off{net}.npz"); L, n = D["s1y"].shape
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64)
g = gac(list(Wc)); gam_gac = np.array([d["gamma"] for d in g])
try:
    dump = pickle.load(open(f"../k3work/dump_oracleall_off{net}.pkl", "rb")); chain = {d["layer"]: d for d in dump}
except Exception: chain = {}
def fit(x, y): return (x @ y) / (x @ x)
def coh_k3(W, k3h, C3):
    """coherent third cumulant of z_k = sum_a W_ka h_a: diagonal + pair patterns (a=b != c)"""
    C3o = C3.copy(); np.fill_diagonal(C3o, 0.0)
    return (W**3) @ k3h + 3 * (((W * W) @ C3o) * W).sum(1)
print(f"net {net}: the kappa3 mean channel per layer")
print(" l | measured g_l (fit kappa3 = 1.5 g mu s^2)  corr | D21 fit g1 (mu_a S_ac), g2 (s_a^2 mu_c)/0.5 | GAC gamma_l | chain D3 fit g, corr | kappa3 share of (1.5 mu s^2)-fit variance")
gm = np.zeros(L)
for l in range(L):
    Y = y_cumulants(D, l); mu, S, k3, K21 = Y["mu"], Y["S"], Y["k3"], Y["K21"]; s2 = np.diag(S)
    x = 1.5 * mu * s2; gm[l] = fit(x, k3); c = np.corrcoef(x, k3)[0, 1]
    iu = np.triu_indices(n, 1); f1 = (mu[:, None] * S)[iu]; f2 = (s2[:, None] * mu[None, :])[iu]; yv = K21[iu]
    A_ = np.stack([f1, f2], 1); coef, *_ = np.linalg.lstsq(A_, yv, rcond=None); r2 = 1 - np.var(yv - A_ @ coef) / np.var(yv)
    ch = ""
    if l in chain:
        d = chain[l]; xc = 1.5 * d["mu"] * d["var"]; ch = f"{fit(xc, d['D3']):.5f}, {np.corrcoef(xc, d['D3'])[0,1]:+.3f}"
    print(f" {l:2d} | {gm[l]:.5f}  {c:+.3f} | {coef[0]:.5f}  {coef[1]/0.5:.5f}  (R^2 {r2:.3f}) | {gam_gac[l]:.5f} | {ch:>16s} | {1-np.var(k3-gm[l]*x)/np.var(k3):.3f}")
print("\none-step generation ledger of the coherent kappa3 of z_{l+1} (fit coefficient against 1.5 mu_z s_z^2):")
print(" l | g_{l+1} measured | from TRUE h cumulants | Gaussian y (one loop) | kappa3(y) variation: total / tree / gain-ind (g_l) | kappa4(y) variation | sum | g_l")
t0 = time.time(); led = []
for l in range(L - 1):
    Y = y_cumulants(D, l); mu, S = Y["mu"], Y["S"]; sig = np.sqrt(np.diag(S))
    A = relu_coeffs(mu, sig, 19); B = relu2_coeffs(mu, sig, 19); lay = Layer(mu, S, A, B)
    W = Wc[l + 1]; muz = D["s1y"][l + 1]; s2z = np.diag(D["S2y"][l + 1] - np.outer(muz, muz)); xz = 1.5 * muz * s2z
    C3t, C4t, k4t = h_cumulants(D, l)
    k3t = np.diag(C3t).copy()
    g_true = fit(xz, coh_k3(W, k3t, C3t)); g_1 = fit(xz, coh_k3(W, np.diag(lay.C3G).copy(), lay.C3G))
    dC3, _, _ = lay.var3(Y["k3"], Y["K21"]); g_3 = fit(xz, coh_k3(W, np.diag(dC3).copy(), dC3))
    gl = gm[l]; k3g = 1.5 * gl * mu * sig**2; K21g = gl * (mu[:, None] * S + 0.5 * (sig**2)[:, None] * mu[None, :])
    dC3g, _, _ = lay.var3(k3g, K21g); g_3g = fit(xz, coh_k3(W, np.diag(dC3g).copy(), dC3g))
    dC3t, _, _ = lay.var3(Y["k3"] - k3g, Y["K21"] - K21g); g_3t = fit(xz, coh_k3(W, np.diag(dC3t).copy(), dC3t))
    dC34, _, _ = lay.var4(Y["k4"], Y["K31"], Y["K22"]); g_4 = fit(xz, coh_k3(W, np.diag(dC34).copy(), dC34))
    led.append((l, gm[l + 1], g_true, g_1, g_3, g_3t, g_3g, g_4, gl))
    print(f" {l:2d} | {gm[l+1]:.5f}          | {g_true:.5f}               | {g_1:.5f}               | {g_3:+.5f} / {g_3t:+.5f} / {g_3g:+.5f}   | {g_4:+.5f}            | {g_1+g_3+g_4:.5f} | {gl:.5f}   ({time.time()-t0:.0f}s)", flush=True)
np.savez(f"k3channel_off{net}.npz", gm=gm, gam_gac=gam_gac, led=np.array(led))
print("\naccumulation: g_l measured vs sum of one-loop increments (critical) vs GAC gamma_l")
c1 = np.concatenate([[0], np.cumsum([r[3] for r in led])])
for l in range(L): print(f" {l:2d} | g_l {gm[l]:.5f} | sum one-loop {c1[l]:.5f} | GAC {gam_gac[l]:.5f} | ratio g/GAC {gm[l]/max(gam_gac[l],1e-12):.3f}")
