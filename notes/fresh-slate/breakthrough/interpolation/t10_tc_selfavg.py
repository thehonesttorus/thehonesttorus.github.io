"""T10: self-averaging readout correction on top of TC (and of the closure): fit e = truth - pred on final-layer
features (alpha, s) over TRAINING networks; evaluate on bench w1024_d16."""
import glob, numpy as np
from common import bench
import tc
from t3_selfavg import basis
def feats(W, **kw):
    p, a, s = tc.predict(W, ret_feats=True, **kw); return p[-1], a, s
files = sorted(glob.glob("results/train/w1024_d16_s*.npz"))
for kw, lab in [(dict(chi=True), "TC+chi"), (dict(inject=False, chi=True), "closure+chi")]:
    X, Y = [], []
    for f in files:
        W = bench.weights_from_seed(int(f.split("_s")[-1][:-4]), 1024, 16); T = np.load(f)["means"][-1]
        p, a, s = feats(W, **kw); X.append(np.stack([a, s], 1)); Y.append(T - p)
    X = np.concatenate(X); Y = np.concatenate(Y)
    S = bench.load_set("w1024_d16"); out = {d: [] for d in (0, 1, 2)}; base = []
    for d in (0, 1, 2):
        c, *_ = np.linalg.lstsq(basis(X, d), Y, rcond=None)
        for i in range(len(S["seeds"])):
            W = bench.weights(S, i); T = S["means"][i][-1]; nz = S["noise"][i]
            p, a, s = feats(W, **kw)
            if d == 0: base.append(((p - T) ** 2).mean() - nz)
            q = p + basis(np.stack([a, s], 1), d) @ c
            out[d].append(((q - T) ** 2).mean() - nz)
    print(lab, f"ntrain {len(files)}: raw {np.mean(base):.3e} ->", " ".join(f"deg{d} {np.mean(v):.3e}" for d, v in out.items()), flush=True)
