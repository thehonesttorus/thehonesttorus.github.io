# One-step closure test of the kappa_4 diagonal (note XXXI).   python k4closure.py NET MCPREFIX [TRUTH]
# (TRUTH: the Monte Carlo file judged against, default full; h1 when the oracle run used h0.)
# Reads chain_off{NET}_o.npz: a chain run with every pre-activation statistic replaced by Monte Carlo truth at every
# layer (oracle D3+D21+G4+WK4M+K31+VAR+COFF), dumping per layer the post-activation (2,2) slice K22 and kappa_4 diagonal
# K4v and variance K2v that the chain's Wick program computes from those true inputs, and g4rowown: the chain's own
# kappa_4 diagonal of the next pre-activation before the oracle overwrote it. So every prediction below is a one-step
# map from true inputs at layer l to kappa_4(z_(l+1)), compared with Monte Carlo truth:
#   chain: the shipped regenerated core (mean-field transport, fitted lambda table, rank-4 quenched pair correction);
#   Q:     the exact pair class 3 (w_i^2)^T K w_i^2 - 2 diag(K)^T w_i^4, K = K22 + diag(K4v) (note XXVIII section 4);
#   Q + D: plus the scale mixture's dropped classes 6 g s_diag^2 s_off^2 + 3 g s_off^4 (note XXI section 5), with
#          s_diag^2 = (W o W) K2v, s_off^2 = var(z_(l+1)) - s_diag^2 and g the mixture gain of layer l read from the
#          true kappa_4 diagonal (g4: least squares of k4 = 3 g var^2) or from the true kappa_3 diagonal (g3: of
#          k3 = 1.5 g mu var).
# Reported: relative error ||pred - MC|| / ||MC||, correlation, mean ratio, and the Monte Carlo noise of the target.
import sys, numpy as np
net, pre = int(sys.argv[1]), sys.argv[2]; tt = sys.argv[3] if len(sys.argv) > 3 else "full"
ch = np.load(f"chain_off{net}_o.npz"); F = np.load(f"{pre}_{tt}.npz"); H0 = np.load(f"{pre}_h0.npz"); H1 = np.load(f"{pre}_h1.npz")
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64)
print(f"net {net}: truth {tt}, {int(F['n'])} samples; one-step predictions of kappa_4(z_(l+1)) from true inputs at layer l")


def stat(p, t):
    return f"rel {np.linalg.norm(p - t) / np.linalg.norm(t):.3f} corr {np.corrcoef(p, t)[0, 1]:.4f} ratio {p.mean() / t.mean():.3f}"


for l in range(1, 15):
    if f"K22_{l}" not in ch.files or f"g4rowown_{l + 1}" not in ch.files:
        continue
    K = ch[f"K22_{l}"].astype(np.float64); K = 0.5 * (K + K.T); np.fill_diagonal(K, ch[f"K4v_{l}"].astype(np.float64))
    W = Wcol[l + 1]; WW = W * W
    Q = 3.0 * np.einsum("ia,ia->i", WW @ K, WW) - 2.0 * (WW * WW) @ np.diag(K)
    t = F["k4"][l + 1].astype(np.float64); nz = np.linalg.norm(H0["k4"][l + 1] - H1["k4"][l + 1]) / 2 / np.linalg.norm(t)
    v1 = F["var"][l + 1].astype(np.float64)
    sd = WW @ ch[f"K2v_{l}"].astype(np.float64); so = v1 - sd
    v0, m0 = F["var"][l].astype(np.float64), F["mu"][l].astype(np.float64)
    g4 = float(np.sum(F["k4"][l] * v0 ** 2) / np.sum(3.0 * v0 ** 4))
    b3 = 1.5 * m0 * v0; g3 = float(np.sum(F["k3"][l] * b3) / np.sum(b3 * b3))
    D4 = 6 * g4 * sd * so + 3 * g4 * so * so
    D3 = 6 * g3 * sd * so + 3 * g3 * so * so
    own = ch[f"g4rowown_{l + 1}"].astype(np.float64)
    # least-squares amplitude of the dropped-class shape on the residual t - Q
    c = float(np.dot(t - Q, D4) / np.dot(D4, D4)) if np.dot(D4, D4) > 0 else float("nan")
    print(f"layer {l:2d}->{l + 1:2d} (MC noise {nz:.3f}; g4 {g4:.4f} g3 {g3:.4f}; mean s_off^2/var {np.mean(so / v1):+.3f}):\n"
          f"    chain {stat(own, t)} | Q {stat(Q, t)}\n"
          f"    Q+D(g4) {stat(Q + D4, t)} | Q+D(g3) {stat(Q + D3, t)} | LS amplitude of D(g4) on MC - Q: {c:.2f}", flush=True)
