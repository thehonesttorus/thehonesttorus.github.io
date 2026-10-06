# Truth-noise floor of the one-step residual r_l = m*_l - chain_l(true mean in): truth noise at l plus the
# injected noise W (m*_{l-1} noise) gated by Phi.  Var(h) from the chain's Gaussian post variance.
import pickle, numpy as np
from scipy.stats import norm
D = {d["layer"]: d for d in pickle.load(open("dump_oracleall_off0.pkl", "rb"))}
pred = np.load("pred_oracleall_off0.npy"); mt = np.load("../official/truth_off0.npz")["m"].astype(np.float64)
Wcol = np.load("../official/W_off0.npy").astype(np.float64)
N = 1e9
def vpost(d):
    s = np.sqrt(d["var"]); a = d["mu"] / s; ph, Ph = norm.pdf(a), norm.cdf(a)
    m = s * (a * Ph + ph); return s**2 * ((a**2 + 1) * Ph + a * ph) - m**2, Ph
for l in sorted(D):
    if l - 1 not in D: continue
    v_l, Ph_l = vpost(D[l]); v_p, _ = vpost(D[l - 1])
    noise = v_l / N + Ph_l**2 * ((Wcol[l]**2) @ v_p) / N
    mse = np.mean((mt[l] - pred[l])**2)
    print(f"layer {l:2d}: local MSE {mse:.2e}  truth-noise floor {noise.mean():.2e}  signal {mse - noise.mean():.2e}  ({noise.mean()/mse:.0%} noise)")
