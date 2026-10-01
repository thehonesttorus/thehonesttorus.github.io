"""Test T2c (B-pseudorandomness): the a-constant (norm-process) channel of the old content through JMR's mean term.
For a fresh layer, E_W ||x diag(P) W||^2 = sum_i x_i^2 P_i^2 c_i (c_i = ||W_i.||^2) ~ ||x||^2 <P^2 c>: the squared
graph's O(n)-average. So the per-row scalars ||Y_r||^2, Y_r.Z_r, ||Z_r||^2, Z_r.T_r of a source can be frozen when it
turns old and multiplied by scalar gains gamma_u = <P_u^2 c^(u+1)> afterwards; the a-constant part of its D21 is then a
VECTOR transported linearly (O(n^2) per layer once merged over sources). Here the frozen-scalar formula replaces the
exact a-constant part (the fluctuation is kept exact) to measure the approximation in isolation.
usage: python t2_annealed.py MLP [AGE]"""
import sys, json
import numpy as np
sys.path.insert(0, '/home/user/thehonesttorus.github.io/notes/fresh-slate/breakthrough/region')
sys.path.insert(0, '/home/user/thehonesttorus.github.io/notes/fresh-slate/bench')
import bench
import fc_hooked as F
from t2_split import contrib

mlp = int(sys.argv[1]); AGE = int(sys.argv[2]) if len(sys.argv) > 2 else 2
DROPF = len(sys.argv) > 3 and sys.argv[3] == 'dropfluct'
S = bench.load_set('w1024_d16')
W = bench.weights(S, mlp); T = S['means'][mlp]
W64 = W.astype(np.float64); colsq = [(w * w).sum(1) for w in W64]   # c_i of W_l (rows i: input index)
F.GATES = []
stash = {}; errs = {}; cumgain = {}

def legs(src):
    Zf = src["Z"].astype(np.float64)
    Y = src["SP"].astype(np.float64) @ Zf
    Tl = src["Delta"].astype(np.float64) @ Zf if "Delta" in src else None
    return Y, Zf, Tl

def hook(t, sources):
    P, l = F.GATES[-1]                      # iteration l = t - 1 transported with diag(P_l) W_{l+1}
    g = float(np.mean(P * P * colsq[l + 1]))
    for src in sources:
        if "q0" in src:
            src["gain"] *= g
    n = sources[0]["Z"].shape[1]
    Dold = np.zeros((n, n)); Dm_ann = np.zeros(n)
    for src in sources:
        age = t - src["s"]
        if age <= AGE:
            continue
        Y, Zf, Tl = legs(src)
        if "q0" not in src:                 # just turned old: freeze exact row scalars
            src["q0"] = dict(yy=(Y * Y).sum(1), yz=(Y * Zf).sum(1),
                             zz=(Zf * Zf).sum(1), zt=(Zf * Tl).sum(1) if Tl is not None else None)
            src["gain"] = 1.0
        q = src["q0"]; gn = src["gain"]; w = src["w2"].astype(np.float64)
        v = ((w * q["yy"]) @ Zf + 2 * ((w * q["yz"]) @ Y)) * gn / n
        if Tl is not None:
            v += ((q["zz"] @ Tl) + 2 * (q["zt"] @ Zf)) * gn / n
        Dm_ann += v
        Dold += contrib(src)
    if Dold.any():
        Dm = Dold.mean(0)
        stash[t] = (Dm, Dm_ann, Dold)
        errs[t] = float(np.linalg.norm(Dm_ann - Dm) / np.linalg.norm(Dm))

def post(t, D):
    if t not in stash:
        return D
    Dm, Dm_ann, Dold = stash[t]
    if DROPF:
        return D - Dold + Dm_ann[None, :]
    return D + (Dm_ann - Dm)[None, :]

F.HOOK = hook; F.POST = post
p = F.run(W, slices=2, k4mf=True)
raw = float(((p[-1] - T[-1]) ** 2).mean() - S['noise'][mlp])
r = dict(mlp=mlp, variant="mean_annealed" + ("+drop_fluct" if DROPF else ""), age=AGE, raw=raw, rel_err_mean_channel={int(k): round(v, 4) for k, v in errs.items()})
print(json.dumps(r), flush=True)
with open('t2_ablate_results.jsonl', 'a') as f: f.write(json.dumps(r) + "\n")
