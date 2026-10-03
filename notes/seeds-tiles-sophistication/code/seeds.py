# Input as seed: iid directions vs antithetic pairs vs randomly rotated cross-polytope frames (3-designs).
# Same number of forward passes M = 2nK. Prediction (infinite width): MSE ratio MC/frames = P(Z>0)/(2 P(Z even>=4)).
import numpy as np, sys
from scipy.special import gammaln
spec = np.load("spec.npy")
n, L = int(sys.argv[1]), int(sys.argv[2])
nets, K, R, Ktruth = 3, 4, 24, int(sys.argv[3]) if len(sys.argv) > 3 else 300
cn = np.sqrt(2) * np.exp(gammaln((n + 1) / 2) - gammaln(n / 2))      # E|X| for X ~ N(0, I_n)
rng = np.random.default_rng(L * 1000 + n)
def fwd(D, Ws):                                   # D: (M, n) unit directions; returns final post-ReLU outputs
    H = D
    for W in Ws: H = np.maximum(H @ W.T, 0)
    return H
def haar(n):
    Q, Rm = np.linalg.qr(rng.standard_normal((n, n))); return Q * np.sign(np.diag(Rm))
def est_mc(Ws, M):
    D = rng.standard_normal((M, n)); D /= np.linalg.norm(D, axis=1, keepdims=True)
    return cn * fwd(D, Ws).mean(0)
def est_anti(Ws, M):
    D = rng.standard_normal((M // 2, n)); D /= np.linalg.norm(D, axis=1, keepdims=True)
    return cn * fwd(np.vstack([D, -D]), Ws).mean(0)
def est_frames(Ws, K):
    out = 0
    for _ in range(K):
        O = haar(n); out = out + fwd(np.vstack([O.T, -O.T]), Ws).mean(0)   # rows of O.T are the columns O e_i
    return cn * out / K
b = spec[L - 1]
pred_anti = (1 - b[0]) / (2 * b[2::2].sum()); pred_frame = (1 - b[0]) / (2 * b[4::2].sum())
res = {"mc": [], "anti": [], "frames": []}
for k in range(nets):
    Ws = [rng.standard_normal((n, n)) * np.sqrt(2.0 / n) for _ in range(L)]
    truth = est_frames(Ws, Ktruth)
    for name, f in [("mc", lambda: est_mc(Ws, 2 * n * K)), ("anti", lambda: est_anti(Ws, 2 * n * K)), ("frames", lambda: est_frames(Ws, K))]:
        res[name].append(np.mean([np.mean((f() - truth)**2) for _ in range(R)]))
mc, an, fr = (np.mean(res[k]) for k in ("mc", "anti", "frames"))
print(f"n={n} L={L}: MSE  MC={mc:.3e}  antithetic={an:.3e}  frames={fr:.3e} | gain anti={mc/an:.2f} (pred {pred_anti:.2f})  frames={mc/fr:.2f} (pred {pred_frame:.2f})"
      f"  [truth noise ~ frames/{Ktruth//K}]")
