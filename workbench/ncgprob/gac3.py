# Two-gain transport ledger: the coherent third cumulant (kappa3(z_k) = 1.5 g3 mu_k s_k^2) and the coherent fourth
# cumulant diagonal (kappa4(z_k) = 3 g4 s_k^4) of the pre-activation z_{l+1} = W h_l, h_l = relu(y_l), generated per
# layer from the TRUE Monte Carlo marginal state (mu_l, S_l) of y_l by
#   one loop   (Gaussian y_l; the fresh non-Gaussianity of relu),
#   kappa3(y_l) variation (exact first-order Gaussian integration by parts, pairvar.Layer.var3),
#   kappa4(y_l) variation (pairvar.Layer.var4),
# each once with the TRUE slices of y_l (ledger check against the measured g3, g4 of z_{l+1}) and once with the
# scale-mixture slices of unit gain (the 2x2 transport matrix M_l of a two-component gain state).  Then the closed
# recursion (g3, g4)_{l+1} = M_l (g3, g4)_l + (phi3, phi4)_l on the true state, against the measured profile.
# Coherent restrictions used:  kappa3(z_k) ~ sum_a W^3 k3_a + 3 sum_{a!=c} W_ka^2 W_kc K21_ac  (pairs (a,a,c));
#                              kappa4(z_k) ~ sum_a W^4 k4_a + 3 sum_{a!=c} W_ka^2 W_kc^2 K22_ac (pairs (a,a,c,c)).
# The dropped slices of h (K111, K31, K211, K1111) are tested once by the ledger column 'from TRUE h'.
import numpy as np, sys, time
sys.path.insert(0, "../num12")
from closure import relu_coeffs, relu2_coeffs
from pairvar import Layer, y_cumulants, h_cumulants, sm_slices
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
D = np.load(f"mc_cum_off{net}.npz"); L, n = D["s1y"].shape
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64)
def fit(x, y): return (x @ y) / (x @ x)
def coh_k3(W, k3h, C3):
    C3o = C3.copy(); np.fill_diagonal(C3o, 0.0)
    return (W**3) @ k3h + 3 * (((W * W) @ C3o) * W).sum(1)
def coh_k4(W, k4h, C4):
    C4o = C4.copy(); np.fill_diagonal(C4o, 0.0); P = W * W
    return (W**4) @ k4h + 3 * ((P @ C4o) * P).sum(1)
# measured gains of z_l (fit coefficients) and the off-diagonal share of the pre-activation variance
gm3 = np.zeros(L); gm4 = np.zeros(L); r2_4 = np.zeros(L); offshare = np.zeros(L)
for l in range(L):
    Y = y_cumulants(D, l); mu, S = Y["mu"], Y["S"]; s2 = np.diag(S)
    x3 = 1.5 * mu * s2; x4 = 3 * s2 * s2
    gm3[l] = fit(x3, Y["k3"]); gm4[l] = fit(x4, Y["k4"]); r2_4[l] = 1 - np.var(Y["k4"] - gm4[l] * x4) / np.var(Y["k4"])
    if l > 0:
        Sh = D["S2h"][l - 1] - np.outer(D["s1h"][l - 1], D["s1h"][l - 1]); W = Wc[l]
        sd = (W * W) @ np.diag(Sh); offshare[l] = 1 - np.mean(sd / s2)
print(f"net {net}, T = {int(D['T'])}: measured gains of z_l: g3 (fit kappa3 = 1.5 g mu s^2), g4 (fit kappa4 = 3 g s^4, R^2), ratio g4/g3, off-diagonal share of Var z (1 - mean W^2 diag(S_h) / s^2)")
for l in range(L): print(f"  {l:2d} | g3 {gm3[l]:+.5f} | g4 {gm4[l]:+.5f} (R^2 {r2_4[l]:.2f}) | g4/g3 {gm4[l]/gm3[l] if abs(gm3[l])>1e-4 else float('nan'):.3f} | off share {offshare[l]:.3f}")
print("\none-step ledger at z_{l+1} (fit coefficients):  [k3 channel]  measured | TRUE h (pair slices) | one loop | var under TRUE kappa3(y) | TRUE kappa4(y) | sum  ||  [k4 channel] same")
print("  l  | g3 meas  TRUEh    1loop    d3       d4       sum    || g4 meas  TRUEh    1loop    d3       d4       sum    || unit-gain sm transport: r33 r34 | r43 r44   (phi3 phi4)")
t0 = time.time(); rows = []
for l in range(L - 1):
    Y = y_cumulants(D, l); mu, S = Y["mu"], Y["S"]; sig = np.sqrt(np.diag(S))
    A = relu_coeffs(mu, sig, 19); B = relu2_coeffs(mu, sig, 19); lay = Layer(mu, S, A, B)
    W = Wc[l + 1]; Yn = y_cumulants(D, l + 1); x3 = 1.5 * Yn["mu"] * np.diag(Yn["S"]); x4 = 3 * np.diag(Yn["S"])**2
    C3t, C4t, k4t = h_cumulants(D, l)
    t3 = fit(x3, coh_k3(W, np.diag(C3t).copy(), C3t)); t4 = fit(x4, coh_k4(W, k4t, C4t))
    p3 = fit(x3, coh_k3(W, np.diag(lay.C3G).copy(), lay.C3G)); p4 = fit(x4, coh_k4(W, lay.k4G, lay.C4G))
    dC3, dC4, dk4 = lay.var3(Y["k3"], Y["K21"]); d33 = fit(x3, coh_k3(W, np.diag(dC3).copy(), dC3)); d43 = fit(x4, coh_k4(W, dk4, dC4))
    dC3, dC4, dk4 = lay.var4(Y["k4"], Y["K31"], Y["K22"]); d34 = fit(x3, coh_k3(W, np.diag(dC3).copy(), dC3)); d44 = fit(x4, coh_k4(W, dk4, dC4))
    sl = sm_slices(mu, S, 1.0)
    dC3, dC4, dk4 = lay.var3(sl["k3"], sl["K21"]); r33 = fit(x3, coh_k3(W, np.diag(dC3).copy(), dC3)); r43 = fit(x4, coh_k4(W, dk4, dC4))
    dC3, dC4, dk4 = lay.var4(sl["k4"], sl["K31"], sl["K22"]); r34 = fit(x3, coh_k3(W, np.diag(dC3).copy(), dC3)); r44 = fit(x4, coh_k4(W, dk4, dC4))
    rows.append((l, gm3[l + 1], t3, p3, d33, d34, gm4[l + 1], t4, p4, d43, d44, r33, r34, r43, r44))
    print(f"  {l:2d} | {gm3[l+1]:.5f} {t3:.5f} {p3:.5f} {d33:+.5f} {d34:+.5f} {p3+d33+d34:.5f} || {gm4[l+1]:.5f} {t4:.5f} {p4:.5f} {d43:+.5f} {d44:+.5f} {p4+d43+d44:.5f} || {r33:.3f} {r34:.3f} | {r43:.3f} {r44:.3f}   ({p3:.5f} {p4:.5f})   [{time.time()-t0:.0f}s]", flush=True)
rows = np.array(rows)
print("\nclosed two-gain recursion on the true state: (g3,g4)_{l+1} = M_l (g3,g4)_l + (phi3,phi4)_l, M_l = [[r33,r34],[r43,r44]]; also the one-gain (homogeneous) recursion g_{l+1} = (r33+r34) g_l + phi3 and the critical sums")
g3 = g4 = 0.0; gh = 0.0; c3 = c4 = 0.0
print("  l  | measured g3  g4  ratio | derived g3  g4  ratio | one-gain g | critical sum phi3  phi4  ratio")
print(f"   0 | {gm3[0]:.5f} {gm4[0]:.5f}   -   | 0.00000 0.00000   -   | 0.00000 | 0.00000 0.00000   -")
for r in rows:
    l, p3, p4 = int(r[0]), r[3], r[8]; r33, r34, r43, r44 = r[11], r[12], r[13], r[14]
    g3, g4 = r33 * g3 + r34 * g4 + p3, r43 * g3 + r44 * g4 + p4; gh = (r33 + r34) * gh + p3; c3 += p3; c4 += p4
    print(f"  {l+1:2d} | {gm3[l+1]:.5f} {gm4[l+1]:.5f} {gm4[l+1]/gm3[l+1]:.3f} | {g3:.5f} {g4:.5f} {g4/g3:.3f} | {gh:.5f} | {c3:.5f} {c4:.5f} {c4/c3:.3f}")
np.savez(f"gac3_off{net}.npz", rows=rows, gm3=gm3, gm4=gm4, offshare=offshare)
