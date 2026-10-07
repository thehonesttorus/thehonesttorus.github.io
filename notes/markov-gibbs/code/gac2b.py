# Two-channel derived gain on the closure's own conditional state (gac2 with the kappa4 channel added):
#     (g3, g4)_{l+1} = M_l (g3, g4)_l + (phi3, phi4)_l,
# phi3, phi4 = one-loop coherent kappa3 / kappa4 gains of z_{l+1} from the Gaussian conditional state (pair slices),
# M_l = unit-gain scale-mixture transport (pairvar.sm_slices through Layer.var3 / var4), optionally with rows rescaled
# to sum 1 (positive homogeneity makes the mixture an exact invariant family).  The output mean is E[G](g) x closure
# mean, with g = g3 or g4 (mode).  No fitted constant.
import numpy as np, sys, time
sys.path.insert(0, "../num12")
from closure import relu_coeffs, relu2_coeffs, closure
from gac import gac, EG
from ledger import gstep
from pairvar import Layer, sm_slices
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
feed = sys.argv[2] if len(sys.argv) > 2 else "g3"     # which gain scales the output and the state: g3 | g4
hom = (sys.argv[3] == "hom") if len(sys.argv) > 3 else False
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64); mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
L, n = Wc.shape[0], Wc.shape[1]
try: Z = np.load(f"gac3_off{net}.npz"); gm3, gm4 = Z["gm3"], Z["gm4"]
except Exception: gm3 = gm4 = np.full(L, np.nan)
def fit(x, y): return (x @ y) / (x @ x)
def coh_k3(W, k3h, C3):
    C3o = C3.copy(); np.fill_diagonal(C3o, 0.0)
    return (W**3) @ k3h + 3 * (((W * W) @ C3o) * W).sum(1)
def coh_k4(W, k4h, C4):
    C4o = C4.copy(); np.fill_diagonal(C4o, 0.0); P = W * W
    return (W**4) @ k4h + 3 * ((P @ C4o) * P).sum(1)
t0 = time.time(); g3 = g4 = 0.0; g1 = 1.0; out = []; log = []
for l, W in enumerate(Wc):
    if l == 0: mu = np.zeros(n); S = W @ W.T
    else:
        nu = W @ mbar; Q = W @ Sbar @ W.T
        g3, g4 = r33 * g3 + r34 * g4 + p3, r43 * g3 + r44 * g4 + p4
        g1 = EG(g3 if feed == "g3" else g4, "gamma"); mu = nu / g1; S = Q - np.outer(mu, mu)
    M, C = gstep(mu, S); mbar = g1 * M; Sbar = C + np.outer(M, M)
    out.append(mbar)
    if l < L - 1:
        sig = np.sqrt(np.diag(S)); A = relu_coeffs(mu, sig, 19); B = relu2_coeffs(mu, sig, 19); lay = Layer(mu, S, A, B)
        Wn = Wc[l + 1]; nu_n = Wn @ mbar; Qn = Wn @ Sbar @ Wn.T; muz = nu_n / g1; s2z = np.diag(Qn) - muz**2; x3 = 1.5 * muz * s2z; x4 = 3 * s2z * s2z
        p3 = fit(x3, coh_k3(Wn, np.diag(lay.C3G).copy(), lay.C3G)); p4 = fit(x4, coh_k4(Wn, lay.k4G, lay.C4G))
        sl = sm_slices(mu, S, 1.0)
        dC3, dC4, dk4 = lay.var3(sl["k3"], sl["K21"]); r33 = fit(x3, coh_k3(Wn, np.diag(dC3).copy(), dC3)); r43 = fit(x4, coh_k4(Wn, dk4, dC4))
        dC3, dC4, dk4 = lay.var4(sl["k4"], sl["K31"], sl["K22"]); r34 = fit(x3, coh_k3(Wn, np.diag(dC3).copy(), dC3)); r44 = fit(x4, coh_k4(Wn, dk4, dC4))
        if hom:
            s3, s4 = r33 + r34, r43 + r44; r33, r34, r43, r44 = r33 / s3, r34 / s3, r43 / s4, r44 / s4
        log.append((l, g3, g4, p3, p4, r33, r34, r43, r44))
        print(f"  layer {l:2d}: g3 {g3:.5f} g4 {g4:.5f} ratio {g4/g3 if g3 > 0 else float('nan'):.3f} (measured {gm3[l]:.5f} {gm4[l]:.5f} {gm4[l]/gm3[l] if gm3[l] > 1e-4 else float('nan'):.3f})  phi {p3:.5f} {p4:.5f}  M [{r33:.3f} {r34:.3f}; {r43:.3f} {r44:.3f}]   [{time.time()-t0:.0f}s]", flush=True)
mse = lambda a, b: np.mean((a - b)**2); scl = lambda v: 8 * (mt[-1] @ (v - mt[-1])) / (mt[-1] @ mt[-1])
oc = closure(list(Wc)); gg = gac(list(Wc))
print(f"\nnet {net} (feed {feed}, {'homogeneity-corrected' if hom else 'measured M'}): final g3 {g3:.5f} g4 {g4:.5f} ratio {g4/g3:.3f} (measured {gm3[-1]:.5f} {gm4[-1]:.5f} {gm4[-1]/gm3[-1]:.3f}), E[G] {g1:.6f}")
print(f"  final MSE: closure {mse(oc[-1]['m'], mt[-1]):.3e} | GAC {mse(gg[-1]['m'], mt[-1]):.3e} (gamma_L {gg[-1]['gamma']:.5f}) | derived {mse(out[-1], mt[-1]):.3e}")
print(f"  8 x residual scale: closure {scl(oc[-1]['m']):+.5f} | GAC {scl(gg[-1]['m']):+.5f} | derived {scl(out[-1]):+.5f}")
np.savez(f"gac2b_off{net}_{feed}_{'hom' if hom else 'meas'}.npz", log=np.array(log), out=np.array(out))
