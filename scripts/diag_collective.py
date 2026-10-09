"""Stability of the collective mode under the gated transport: for sources born at layer l' (top eigenvector v of the
pre-activation covariance Soff, p = Phi o v), how parallel is T_{l<-l'} p to the collective direction at layer l
(the pre-activation mean m_l, and the top eigenvector v_l)?  Also the rank of the family {T_{l<-l'} p_{l'}}_{l'<l}
and the spectral energy of Soff in its top 1/4/16/64 eigenpairs.
  python scripts/diag_collective.py DATA NET"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.kprop3c import kprop3c_chain
from whest.kprop3 import wick, _zero_diag
D, net = sys.argv[1], int(sys.argv[2])
W = np.load(f"{D}/W_off{net}.npy").astype(np.float64); L, n, _ = W.shape
rec = {}; kprop3c_chain(W, dict(window=99, k=8), record=rec)
Phi = {l: wick(rec[l]["m"], rec[l]["var"], 1, 1) for l in rec}
S = {0: W[0] @ W[0].T}
for l in range(1, L): S[l] = W[l] @ rec[l - 1]["Ch"] @ W[l].T
v, lam, p, q = {}, {}, {}, {}
np.set_printoptions(precision=3, suppress=True, linewidth=220)
print("layer: top-1/4/16/64 energy fraction of Soff; lam_1/lam_2; cos(v, m)")
for l in range(L):
    Soff = _zero_diag(S[l]); e, V = np.linalg.eigh(Soff); idx = np.argsort(-np.abs(e)); e, V = e[idx], V[:, idx]
    sg = np.sign(V[:, 0] @ rec[l]["m"]) if np.linalg.norm(rec[l]["m"]) > 1e-12 else np.sign(V[:, 0].sum()); sg = sg if sg != 0 else 1.0
    v[l], lam[l] = V[:, 0] * sg, e[0]; p[l] = Phi[l] * v[l]; q[l] = wick(rec[l]["m"], rec[l]["var"], 2, 1) * v[l]
    cap = [np.sum(e[:k] ** 2) / np.sum(e ** 2) for k in (1, 4, 16, 64)]
    m = rec[l]["m"]; print(f"  {l:2d}: {cap[0]:.3f} {cap[1]:.3f} {cap[2]:.3f} {cap[3]:.3f} | {e[0]/abs(e[1]):.2f} | {m @ v[l] / np.linalg.norm(m):+.4f}")
def transport(x, lp, l):
    # T_{l<-l'} x = W_l D(Phi_{l-1}) ... D(Phi_{l'+1}) W_{l'+1} x   (legs are gated by Phi after each readout)
    y = W[lp + 1] @ x
    for j in range(lp + 2, l + 1): y = W[j] @ (Phi[j - 1] * y)
    return y
print("\nreadout layer l | cos(T p_{l'}, m_l) for l' = 0..l-2 | singular values of the normalised family (top 4) | cos(T q_{l'}, m_l) mean")
for l in (4, 6, 8, 10, 12, 15):
    m = rec[l]["m"]; mh = m / np.linalg.norm(m); cols = []; cq = []
    cs = []
    for lp in range(0, l - 1):
        y = transport(p[lp], lp, l); yh = y / np.linalg.norm(y); cols.append(yh); cs.append(yh @ mh)
        z = transport(q[lp], lp, l); cq.append(z @ mh / np.linalg.norm(z))
    P = np.stack(cols, axis=1); sv = np.linalg.svd(P, compute_uv=False)
    print(f"  {l:2d} | {np.array(cs)} | {sv[:4] / sv[0]} | {np.mean(cq):+.3f}")
    print(f"       cos(T p_{{l'}}, v_l): {np.array([c @ v[l] for c in cols])}")
