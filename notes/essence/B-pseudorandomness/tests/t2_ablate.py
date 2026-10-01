"""Test T2b (B-pseudorandomness): which part of the old (age > AGE) third-order content do the final means need?
Inside FC (slices=2, k4mf=True) on bench w1024_d16, at every target layer the old sources' D21 contribution
D_old is split into its a-constant part 1 (1^T D_old / n) (the O(n)-average over the last fresh layer; equivalently
Cov(norm process, z_b)) and the remainder (Wishart fluctuation). Variants remove one part from the D21 the chain
uses.  usage: python t2_ablate.py MLP VARIANT [AGE]   VARIANT in base, drop_fluct, drop_mean, drop_old"""
import sys, json
import numpy as np
sys.path.insert(0, '/home/user/thehonesttorus.github.io/notes/fresh-slate/breakthrough/region')
sys.path.insert(0, '/home/user/thehonesttorus.github.io/notes/fresh-slate/bench')
import bench
import fc_hooked as F
from t2_split import contrib

mlp = int(sys.argv[1]); var = sys.argv[2]; AGE = int(sys.argv[3]) if len(sys.argv) > 3 else 2
stash = {}
def hook(t, sources):
    old = [src for src in sources if t - src["s"] > AGE]
    stash[t] = sum(contrib(s) for s in old) if old else None
def post(t, D):
    Do = stash.get(t)
    if Do is None or var == "base":
        return D
    Dm = np.broadcast_to(Do.mean(0, keepdims=True), Do.shape)
    if var == "drop_fluct": return D - (Do - Dm)
    if var == "drop_mean": return D - Dm
    if var == "drop_old": return D - Do
    raise ValueError(var)
F.HOOK = hook; F.POST = post
S = bench.load_set('w1024_d16')
W = bench.weights(S, mlp); T = S['means'][mlp]
p = F.run(W, slices=2, k4mf=True)
raw = float(((p[-1] - T[-1]) ** 2).mean() - S['noise'][mlp])
r = dict(mlp=mlp, variant=var, age=AGE, raw=raw)
print(json.dumps(r), flush=True)
with open('t2_ablate_results.jsonl', 'a') as f: f.write(json.dumps(r) + "\n")
