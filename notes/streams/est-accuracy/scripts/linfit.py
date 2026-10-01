"""Linear-response combination of levers.  For lever k run at scale c_k, D_k = (p_k - p_base)/c_k.  To first order
MSE(x) = MSE_0 + 2 sum_k x_k <D_k, e> + sum_kl x_k x_l <D_k, D_l>   (e = p_base - truth; per MLP, final layer)
Minimised over x (ridge-free); reported in-sample and leave-one-MLP-out, with the predicted gain.

  python linfit.py --res RESDIR --base TAG  lever1:scale1 lever2:scale2 ...
"""
import argparse, json, os
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--res", required=True)
ap.add_argument("--base", required=True)
ap.add_argument("--truth", default="/tmp/claude-0/ea/data/dev6/npy/truth.npy")
ap.add_argument("--layer", type=int, default=-1)
ap.add_argument("levers", nargs="+")
a = ap.parse_args()
T = np.load(a.truth)


def preds(tag):
    z = np.load(os.path.join(a.res, "preds", f"{tag}.npz"))
    return z["preds"][:, a.layer], list(z["mlps"])


pb, mb = preds(a.base)
e = pb - T[mb, a.layer]
D = []
names = []
for lv in a.levers:
    tag, sc = lv.split(":")
    p, m = preds(tag)
    assert m == mb
    D.append((p - pb) / float(sc))
    names.append(tag)
D = np.stack(D)  # (K, M, n)
K, M, n = D.shape
b = np.einsum("kmi,mi->km", D, e) / n      # per MLP
G = np.einsum("kmi,lmi->klm", D, D) / n
mse0 = (e ** 2).mean(1)


def solve(ms):
    return np.linalg.solve(G[:, :, ms].sum(-1), -b[:, ms].sum(-1))


x = solve(list(range(M)))
pred = mse0 + 2 * x @ b + np.einsum("k,klm,l->m", x, G, x)
print("in-sample x:", dict(zip(names, np.round(x, 3))))
print("predicted mean MSE %.4e -> %.4e (delta %.3e)" % (mse0.mean(), pred.mean(), pred.mean() - mse0.mean()))
loo = []
for m in range(M):
    xm = solve([j for j in range(M) if j != m])
    loo.append(2 * xm @ b[:, m] + xm @ G[:, :, m] @ xm)
    print(f"  LOO mlp {mb[m]}: x={np.round(xm, 3)}  delta {loo[-1]:+.3e}")
print("LOO mean delta %.3e  (s.e. %.1e)" % (np.mean(loo), np.std(loo, ddof=1) / np.sqrt(M)))
for k in range(K):
    xs = -b[k].sum() / G[k, k].sum()
    print(f"single {names[k]}: x*={xs:.3f}  gain {(2 * xs * b[k] + xs * xs * G[k, k]).mean():+.3e}  per-MLP <D,e>/|D|^2: {np.round(-b[k] / G[k, k], 2)}")
