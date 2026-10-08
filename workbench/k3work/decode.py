# Decode the chain's one-step residual through the KNOWN fresh weights.
# r_l(i) = m*_l(i) - chain_l(i) with the true mean entering layer l. Fresh-weight lemma:
# r_l(i) ~ (phi_i / 2 sigma_i) * w_i^T dC_{l-1} w_i + (skew channel), so a state-level hypothesis M for
# dC_{l-1} predicts the per-neuron signature s_M(i) = (phi_i/2sigma_i) w_i^T M w_i.
import pickle, sys, numpy as np
from scipy.stats import norm
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
D = {d["layer"]: d for d in pickle.load(open(f"dump_oracleall_off{net}.pkl", "rb"))}
pred = np.load(f"pred_oracleall_off{net}.npy"); mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64)
rng = np.random.default_rng(0)
def quad(W, M):  # diag(W M W^T)
    return np.einsum("ij,ij->i", W @ M, W)
def offd(M):
    M = M.copy(); np.fill_diagonal(M, 0.0); return M
print("layer  rms(r)   R2: trace  mm    lin   mehl2  gainK4  d21   diagv  skewD3  rand | joint  joint-noise-floor")
for l in sorted(D):
    if l - 1 not in D or D[l - 1]["C_pre"] is None or l == 15 and False:
        continue
    p = D[l - 1]; c = D[l]
    W = Wcol[l]
    r = mt[l] - pred[l]
    # post-activation pieces at l-1 (chain state, Gaussian Wick functions)
    s = np.sqrt(p["var"]); a = p["mu"] / s; ph = norm.pdf(a); Ph = norm.cdf(a)
    m = s * (a * Ph + ph); F1 = Ph; F2 = ph / s
    Coff = offd(p["C_pre"]); v_post = s**2 * ((a**2 + 1) * Ph + a * ph) - m**2
    gam = 0.0184  # gain variance near depth (GAC); only the shape matters in the regression
    H = {
        "trace": np.eye(1024),
        "mm": np.outer(m, m),
        "lin": np.outer(F1, F1) * Coff,
        "mehl2": 0.5 * np.outer(F2, F2) * Coff**2,
        "gainK4": offd(0.25 * gam * np.outer(F2 * s**2, F2 * s**2)),
        "d21": offd(np.outer(F1, F2) * p["D21"].T + np.outer(F2, F1) * p["D21"]) * 0.5 if p["D21"] is not None else None,
        "diagv": np.diag(v_post),
    }
    R = rng.standard_normal((1024, 1024)); H["rand"] = offd(R + R.T) * 1e-3
    # sensitivity at layer l
    sl = np.sqrt(c["var"]); al = c["mu"] / sl; sens = norm.pdf(al) / (2 * sl)
    F3l = -al * norm.pdf(al) / sl**2
    feats = {}
    for k, M in H.items():
        if M is None: continue
        feats[k] = sens * quad(W, M)
    feats["skewD3"] = F3l / 6.0 * c["D3"]
    def r2(X):
        X = np.column_stack([np.ones(1024)] + X)
        beta, *_ = np.linalg.lstsq(X, r, rcond=None)
        return 1 - np.var(r - X @ beta) / np.var(r)
    out = {k: r2([f]) for k, f in feats.items()}
    joint = r2([feats[k] for k in feats if k != "rand"])
    noise = 0.0
    print(f"{l:4d}  {np.sqrt(np.mean(r**2)):.2e}  " + "  ".join(f"{out.get(k, float('nan')):.3f}" for k in ["trace","mm","lin","mehl2","gainK4","d21","diagv","skewD3","rand"]) + f" | {joint:.3f}")
