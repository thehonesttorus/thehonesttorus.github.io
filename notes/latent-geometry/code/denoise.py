# Is the chain's (2,1) third-cumulant slice error high-rank, so that projecting the read-weighted slice onto the
# collective subspace (top-r eigenvectors of the gated covariance) removes error rather than signal?
#   python denoise.py NET   (needs mc2_off{NET}_{full,h0,h1}.npz and chain_off{NET}_fk_base.npz in the working dir)
# Per layer, in the pair program's read metric R = diag(c) D21 diag(Phi), c = phi(alpha)/sigma (truth's alpha):
#   chain error |R_chain - R_true|^2 / |R_true|^2, the same after R_chain -> R_chain U U^T, the truth's own share
#   outside span(U), and the Monte Carlo noise floor |R_h0 - R_h1|^2 / 4 / |R_true|^2.
import sys, numpy as np
from math import sqrt, pi, erf

net = int(sys.argv[1])
F = {h: np.load(f"mc2_off{net}_{h}.npz") for h in ("full", "h0", "h1")}
ch = np.load(f"chain_off{net}_fk_base.npz")
print("mc2 keys:", sorted(F["full"].files)[:40])
print("chain keys:", sorted(ch.files)[:40])
phi = lambda a: np.exp(-0.5 * a * a) / sqrt(2 * pi)
Phi = np.vectorize(lambda a: 0.5 * (1 + erf(a / sqrt(2))))


def get(f, names, l):
    for nm in names:
        for key in (f"{nm}_{l}", f"{nm}{l}"):
            if key in f.files:
                return f[key].astype(np.float64)
    if any(nm in f.files for nm in names):
        nm = next(nm for nm in names if nm in f.files)
        a = f[nm]
        if a.ndim >= 1 and a.shape[0] > l:
            return a[l].astype(np.float64)
    return None


for l in range(1, 15):
    Dt = get(F["full"], ("D21", "d21", "K21z"), l)
    D0 = get(F["h0"], ("D21", "d21", "K21z"), l); D1 = get(F["h1"], ("D21", "d21", "K21z"), l)
    mu = get(F["full"], ("mu", "m1", "mean"), l); var = get(F["full"], ("var", "v", "variance"), l)
    C = get(F["full"], ("C", "cov", "C_pre", "Cz"), l)
    Dc = get(ch, ("D21",), l)
    if Dt is None or Dc is None or mu is None or var is None:
        print(f"layer {l}: missing ({Dt is None}, {Dc is None}, {mu is None}, {var is None})"); continue
    s = np.sqrt(np.maximum(var, 1e-12)); al = mu / s; c = phi(al) / s; P = Phi(al)
    rd = lambda D: c[:, None] * D * P[None, :]
    Rt, Rc = rd(Dt), rd(Dc); nt = np.sum(Rt * Rt)
    if C is None:
        Coff = get(ch, ("C_off",), l); C = Coff - np.diag(np.diag(Coff)) + np.diag(var)
    ev, U = np.linalg.eigh(P[:, None] * C * P[None, :]); U = U[:, ::-1]
    line = f"layer {l:2d}  chain err {np.sum((Rc - Rt) ** 2) / nt:.4f}"
    if D0 is not None and D1 is not None:
        line += f"  MC floor {np.sum((rd(D0) - rd(D1)) ** 2) / 4 / nt:.4f}"
    for r in (64, 128, 256, 512):
        Ur = U[:, :r]; Pc = (Rc @ Ur) @ Ur.T; Pt = (Rt @ Ur) @ Ur.T
        line += (f"  r{r}: proj err {np.sum((Pc - Rt) ** 2) / nt:.4f} (truth outside {np.sum((Rt - Pt) ** 2) / nt:.4f},"
                 f" err in-span {np.sum((Pc - Pt) ** 2) / nt:.4f})")
    print(line, flush=True)
