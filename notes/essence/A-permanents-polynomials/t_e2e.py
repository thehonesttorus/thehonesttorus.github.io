"""End-to-end raw MSE (bench w1024_d16) of FC with old sources (age >= A) carried only in the Bethe sector
(rank-one mean-field of the squared legs, oracle from the exact legs), vs dropping them (window) and vs FC.
usage: python t_e2e.py MLPS A"""
import sys, json
import numpy as np
sys.path.insert(0, "../../fresh-slate/bench")
import bench, fc_hook as F
S = bench.load_set("w1024_d16"); A = int(sys.argv[2])
for i in [int(x) for x in sys.argv[1].split(",")]:
    W = bench.weights(S, i); T = S["means"][i]
    F.MFAGE = A; p = F.run(W, slices=2, k4mf=True)
    F.MFAGE = None; q = F.run(W, slices=2, k4mf=True, window=A - 1)
    r = dict(mlp=i, A=A, bethe_raw=float(((p[-1] - T[-1]) ** 2).mean() - S["noise"][i]),
             drop_raw=float(((q[-1] - T[-1]) ** 2).mean() - S["noise"][i]))
    print(json.dumps(r), flush=True)
