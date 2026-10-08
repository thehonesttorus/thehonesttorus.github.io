# Derived gain-aware closure: the gain is the coefficient g of the mean-coupled third cumulant kappa3(z_k) = 1.5 g mu_k sigma_k^2,
# generated and transported layer by layer on the closure's own conditional state:
#     g_{l+1} = rho_l g_l + phi_l,
# phi_l = coefficient of the coherent kappa3 of z_{l+1} created from the Gaussian conditional state (one loop),
# rho_l = coefficient created by the first variation under the scale-mixture kappa3 / kappa4 slices of unit gain (the channel's
#         own transport through relu and the next weights).  No fitted constant.  Output E[G](g_L) x closure mean (Gamma law).
import numpy as np, sys, time
sys.path.insert(0, "../num12")
from closure import relu_coeffs, relu2_coeffs, closure
from gac import gac, EG
from ledger import gstep
from pairvar import Layer, sm_slices
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
mode = sys.argv[2] if len(sys.argv) > 2 else "full"   # full: rho and phi derived; crit: rho = 1 (critical) with derived phi
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64); mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
L, n = Wc.shape[0], Wc.shape[1]
try: gm = np.load(f"k3channel_off{net}.npz")["gm"]
except Exception: gm = np.full(L, np.nan)
def fit(x, y): return (x @ y) / (x @ x)
def coh_k3(W, k3h, C3):
    C3o = C3.copy(); np.fill_diagonal(C3o, 0.0)
    return (W**3) @ k3h + 3 * (((W * W) @ C3o) * W).sum(1)
t0 = time.time(); gam = 0.0; g1 = 1.0; out = []; log = []
for l, W in enumerate(Wc):
    if l == 0: mu = np.zeros(n); S = W @ W.T
    else:
        nu = W @ mbar; Q = W @ Sbar @ W.T
        gam = rho * gam + phi if mode == "full" else gam + phi
        g1 = EG(gam, "gamma"); mu = nu / g1; S = Q - np.outer(mu, mu)
    M, C = gstep(mu, S); mbar = g1 * M; Sbar = C + np.outer(M, M)
    out.append(mbar)
    if l < L - 1:
        sig = np.sqrt(np.diag(S)); A = relu_coeffs(mu, sig, 19); B = relu2_coeffs(mu, sig, 19); lay = Layer(mu, S, A, B)
        Wn = Wc[l + 1]; nu_n = Wn @ mbar; Qn = Wn @ Sbar @ Wn.T; muz = nu_n / g1; s2z = np.diag(Qn) - muz**2; xz = 1.5 * muz * s2z
        phi = fit(xz, coh_k3(Wn, np.diag(lay.C3G).copy(), lay.C3G))
        sl = sm_slices(mu, S, 1.0)
        d3a, _, _ = lay.var3(sl["k3"], sl["K21"]); d3b, _, _ = lay.var4(sl["k4"], sl["K31"], sl["K22"])
        rho = fit(xz, coh_k3(Wn, np.diag(d3a + d3b).copy(), d3a + d3b))
        rho3 = fit(xz, coh_k3(Wn, np.diag(d3a).copy(), d3a))
        log.append((l, gam, phi, rho, rho3))
        print(f"  layer {l:2d}: g_l {gam:.5f} (measured {gm[l]:.5f})  phi {phi:.5f}  rho {rho:.4f} (kappa3 part {rho3:.4f})   [{time.time()-t0:.0f}s]", flush=True)
mse = lambda a, b: np.mean((a - b)**2); scl = lambda v: 8 * (mt[-1] @ (v - mt[-1])) / (mt[-1] @ mt[-1])
oc = closure(list(Wc)); gg = gac(list(Wc)); c_or = (oc[-1]["m"] @ (oc[-1]["m"] - mt[-1])) / (oc[-1]["m"] @ oc[-1]["m"])
print(f"\nnet {net} ({mode}): final g_L {gam:.5f} (measured {gm[-1]:.5f}), E[G] {g1:.6f}")
print(f"  final MSE: closure {mse(oc[-1]['m'], mt[-1]):.3e} | GAC {mse(gg[-1]['m'], mt[-1]):.3e} (gamma_L {gg[-1]['gamma']:.5f}) | derived {mse(out[-1], mt[-1]):.3e} | oracle scale {mse(oc[-1]['m']*(1-c_or), mt[-1]):.3e}")
print(f"  8 x residual scale: closure {scl(oc[-1]['m']):+.5f} | GAC {scl(gg[-1]['m']):+.5f} | derived {scl(out[-1]):+.5f}")
print("  per layer MSE (closure / GAC / derived): " + " ".join(f"{mse(oc[l]['m'], mt[l]):.1e}/{mse(gg[l]['m'], mt[l]):.1e}/{mse(out[l], mt[l]):.1e}" for l in range(L)))
np.savez(f"gac2_off{net}_{mode}.npz", log=np.array(log), out=np.array(out))
