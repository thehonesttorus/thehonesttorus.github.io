"""Spectrum of the gated transport cocycle T_{l<-l'} = W_l D(Phi_l) ... W_{l'+1} D(Phi_{l'+1}) on an official network, and the
energy of transported kappa3 leg blocks captured by its top-k left singular subspace.
  python scripts/diag_cocycle.py DATA NET [layers=4,8,12,15] [ages=1,2,3,4,6,8] [ks=16,32,64,128,256]"""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scipy.special import ndtr
from whest.kprop3 import kprop3_chain
D, net = sys.argv[1], int(sys.argv[2]); kw = dict(layers="4,8,12,15", ages="1,2,3,4,6,8", ks="16,32,64,128,256")
for a in sys.argv[3:]:
    k, v = a.split("="); kw[k] = v
layers = [int(x) for x in kw["layers"].split(",")]; ages = [int(x) for x in kw["ages"].split(",")]; ks = [int(x) for x in kw["ks"].split(",")]
W = np.load(f"{D}/W_off{net}.npy").astype(np.float64); L, n, _ = W.shape
rec = {}; t0 = time.time(); kprop3_chain(W, record=rec); print(f"chain run {time.time()-t0:.0f}s", flush=True)
Phi = [ndtr(rec[l]["m"] / np.sqrt(rec[l]["var"])) for l in range(L)]        # gates of layer l (pre-activation z_l = W[l] h_{l-1})
rng = np.random.default_rng(0)
print("layer l, age a: singular values of T_{l<-l-a} (ratios s_k/s_1 at k = 1,2,4,8,16,32,64,128,256) | fraction of ||T X||_F^2 in top-k left subspace for a random X (n x 2n), k in", ks)
for l in layers:
    for a in ages:
        if l - a < 0: continue
        T = np.eye(n)
        for j in range(l - a + 1, l + 1):           # apply W_j D(Phi_j) for j = l-a+1 .. l  (gates of layer j act on h_{j-1}... use Phi of layer j-1 as the gate before W_j)
            T = W[j] @ (Phi[j - 1][:, None] * T)
        U, s, Vt = np.linalg.svd(T, full_matrices=False)
        X = rng.standard_normal((n, 2 * n)) / np.sqrt(n); TX = T @ X; tot = np.sum(TX ** 2)
        fr = [np.sum((U[:, :k].T @ TX) ** 2) / tot for k in ks]
        idx = [1, 2, 4, 8, 16, 32, 64, 128, 256]
        print(f"l={l:2d} a={a}: s_k/s_1 " + " ".join(f"{s[i-1]/s[0]:.3f}" for i in idx) + "  | top-k energy " + " ".join(f"{f:.3f}" for f in fr) + f"  | ||T||_F^2/s_1^2 = {np.sum(s**2)/s[0]**2:.1f}", flush=True)
