# Does the quartic weight dependence of the kappa4 pair class carry per-neuron content the chain's mean-field
# row-sum transport drops?  For each source layer s (target t = s+1), with W = W_t:
#   exact pair class      Q   = W^o4 K4 + 3 diag(W^o2 K22 W^o2^T)            (K22 zero-diagonal post-activation (2,2) slice)
#   chain's mean field    Qmf = 2 W^o2 g_prev,  g_prev = cA (K4 + K22 1) + cI (sum K4 + sum K22)   (est_v29 regen core)
#   quenched part         Q - Qmf
# evaluated on the Monte Carlo truth (isolates the weighting from K22 errors) and on the chain's own slices, and compared
# with the chain's per-neuron kappa4 residual (true total kappa4(z_t) - chain g4row_t).
import sys, pickle, numpy as np
from scipy.stats import norm
net = int(sys.argv[1]); srcs = [int(v) for v in sys.argv[2].split(",")]
Z = np.load(f"k4mc_off{net}.npz"); T = float(Z["T"]); n = 1024
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64)
D = {d["layer"]: d for d in pickle.load(open(f"v29k4_off{net}.pkl", "rb"))["dumps"]}
cA, cI = 6.0 / (n + 4.0), -3.0 / ((n + 2.0) * (n + 4.0))
def mf(W2, K4, K22):
    g = cA * (K4 + K22.sum(1)) + cI * (K4.sum() + K22.sum())
    return 2.0 * W2 @ g
def exact(W2, W4, K4, K22):
    return W4 @ K4 + 3.0 * np.einsum("ij,ij->i", W2 @ K22, W2)
def corr(a, b): return float(np.corrcoef(a, b)[0, 1])
def rms(a): return float(np.sqrt(np.mean(a * a)))
print(f"net {net}, Monte Carlo T = {int(T)} per half (two halves)")
for s in srcs:
    t = s + 1; W = Wc[t]; W2 = W * W; W4 = W2 * W2
    res = []
    for hf in (0, 1, None):
        get = (lambda k: Z[f"h{hf}_s{s}_{k}"] / T) if hf is not None else (lambda k: 0.5 * (Z[f"h0_s{s}_{k}"] + Z[f"h1_s{s}_{k}"]) / T)
        x1 = get("s1"); C = get("C") - np.outer(x1, x1); c = np.diag(C).copy()
        M22 = get("M22"); K22 = M22 - np.outer(c, c) - 2 * C * C; np.fill_diagonal(K22, 0.0)
        K4 = get("s4") - 3 * c * c
        z1, z2, z3, z4 = (get(k) for k in ("z1", "z2", "z3", "z4"))
        k4tot = z4 - 4 * z3 * z1 - 3 * z2 * z2 + 12 * z2 * z1 * z1 - 6 * z1 ** 4
        E4, E13, E22, E112 = (get(k) for k in ("E4", "E13", "E22", "E112"))
        R4 = E4; R22 = 3 * (E22 - E4)
        Dg = W2 @ c; Q = np.einsum("ij,ij->i", W2 @ (C * C), W2); G4 = 3 * W4 @ (c * c); G22 = 3 * (Dg * Dg - W4 @ (c * c)) + 6 * (Q - W4 @ (c * c))
        pair_cls = (R4 - G4) + (R22 - G22)
        Qx = exact(W2, W4, K4, K22); Qm = mf(W2, K4, K22)
        res.append(dict(pair_cls=pair_cls, Qx=Qx, Qm=Qm, dq=Qx - Qm, k4tot=k4tot, K22=K22, K4=K4))
    h0, h1, r = res
    noise = lambda key: rms(h0[key] - h1[key]) / 2.0          # rms noise of the two-half mean
    dc, dt = D[s], D[t]
    K22c = dc["K22"].copy(); np.fill_diagonal(K22c, 0.0); K4c = dc["K4v"]
    Qxc, Qmc = exact(W2, W4, K4c, K22c), mf(W2, K4c, K22c); dqc = Qxc - Qmc
    g4 = dt["g4row"]; lam_part = g4 - Qmc; svar = dt["var"] - W2 @ dc["K2v"]
    resid = r["k4tot"] - g4
    sig = np.sqrt(dt["var"]); al = dt["mu"] / sig; lev = (al * al - 1) * norm.pdf(al) / sig ** 3 / 24.0   # d(mean)/d(kappa4)
    print(f"\n source {s} -> target {t}")
    print(f"  checks: pair class (power sums) vs exact contraction of true K22: rel diff {rms(r['pair_cls'] - r['Qx']) / rms(r['Qx']):.2e};"
          f"  chain lam part vs lam*(var - W2 var_prev): corr {corr(lam_part, svar):.4f}")
    print(f"  truth: pair class mean {r['Qx'].mean():+.5f} rms {rms(r['Qx']):.5f} (noise {noise('Qx'):.5f}); mean-field mean {r['Qm'].mean():+.5f}"
          f"; quenched part rms {rms(r['dq']):.5f} (noise {noise('dq'):.5f}), corr with pair class {corr(r['dq'], r['Qx']):+.3f}")
    print(f"  truth: total kappa4 mean {r['k4tot'].mean():+.5f} rms {rms(r['k4tot']):.5f} (noise {noise('k4tot'):.5f})")
    print(f"  chain: pair class (own slices) mean {Qxc.mean():+.5f}, mean-field {Qmc.mean():+.5f}, quenched rms {rms(dqc):.5f}; "
          f"corr(chain quenched, true quenched) {corr(dqc, r['dq']):+.3f}; corr(chain exact, true pair) {corr(Qxc, r['Qx']):+.3f} vs mean-field {corr(Qmc, r['Qx']):+.3f}")
    print(f"  chain g4row: mean {g4.mean():+.5f}; residual true - g4row: mean {resid.mean():+.5f} rms {rms(resid):.5f} "
          f"(noise {noise('k4tot'):.5f}); corr(residual, chain quenched) {corr(resid, dqc):+.3f}, corr(residual, true quenched) {corr(resid, r['dq']):+.3f}")
    for tag, add in (("chain quenched", dqc), ("true quenched", r["dq"])):
        b = float(add @ (resid - resid.mean()) / (add @ add)); left = resid - b * add
        print(f"    adding {tag}: LS coefficient {b:+.3f}; residual rms {rms(resid - resid.mean()):.5f} -> {rms(left - left.mean()):.5f}; at coefficient 1: "
              f"{rms(resid - add - (resid - add).mean()):.5f}")
    print(f"  mean leverage: rms d(mean) of the chain quenched part {rms(lev * dqc):.2e}, of the true quenched part {rms(lev * r['dq']):.2e}, "
          f"of the chain's kappa4 residual {rms(lev * (resid - resid.mean())):.2e} (noise {rms(lev) * noise('k4tot'):.2e})")
    kk = r["K22"][np.triu_indices(n, 1)]; kc = K22c[np.triu_indices(n, 1)]
    print(f"  K22 slice, chain vs truth: corr {corr(kc, kk):+.3f}, rel err {rms(kc - kk) / rms(kk):.3f}; row sums corr {corr(K22c.sum(1), r['K22'].sum(1)):+.3f}")
