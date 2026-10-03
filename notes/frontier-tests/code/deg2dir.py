# Where does the degree-2 energy live across neurons? Error of a frame average vs an iid average:
# degree-2 energy per neuron-direction = (Var_iid - Var_frame) * n. Split into the mean direction m^ and its complement,
# and the top principal directions of the output covariance.
import numpy as np
from deg2 import net, fwd
def run(Wf, n, reps):
    rng = np.random.default_rng(3); fr, iid = [], []
    for _ in range(reps):
        Q = np.linalg.qr(rng.standard_normal((n, n)))[0].astype(np.float32); fr.append(fwd(Wf, Q).mean(0))
        T = rng.standard_normal((n, n)).astype(np.float32); T /= np.linalg.norm(T, axis=1, keepdims=True); iid.append(fwd(Wf, T).mean(0))
    fr, iid = np.array(fr), np.array(iid); m = iid.mean(0); mh = m/np.linalg.norm(m)
    # covariance of the estimator errors; degree-2 part = Cov_iid - Cov_frame
    Ci, Cf = np.cov(iid.T), np.cov(fr.T); D2 = Ci - Cf
    tot = np.trace(D2); along = mh @ D2 @ mh
    ev, V = np.linalg.eigh(Ci); top = V[:, -8:]
    return tot/np.trace(Ci), along/tot, np.trace(top.T @ D2 @ top)/tot, (mh @ Ci @ mh)/np.trace(Ci)
for (n, L, src) in [(1024, 16, "official"), (256, 32, "he"), (256, 16, "he")]:
    Wf = [np.ascontiguousarray(W.T) for W in np.load("../official/W_off0.npy")] if src == "official" else net(n, L, 7)
    v2, al, t8, cm = run(Wf, n, 64 if n == 1024 else 400)
    print(f"n={n:5d} L={L:2d}: degree-2 share {v2:.3f}; of the degree-2 energy: along mean direction {al:.3f}, in top-8 output PCs {t8:.3f}  (mean direction carries {cm:.3f} of total variance)", flush=True)
