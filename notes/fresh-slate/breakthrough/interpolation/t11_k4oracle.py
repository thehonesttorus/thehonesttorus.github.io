"""T11: TC + scalar kappa4 channel with ORACLE excess norm variance X_l = Var|a~_l|^2 - 2||C_l||^2 (MC Var from t1,
C from the TC chain itself), MLP 0 at n=1024."""
import json, numpy as np
from common import bench
import tc
S = bench.load_set("w1024_d16"); W = bench.weights(S, 0); T = S["means"][0]; nz = S["noise"][0]; n = 1024
vt = np.array(json.load(open("results/t1_w1024_d16_0_N65536.json"))["var_tau"])[1:]   # relative Var tau of a_l
p0, Cs = tc.predict(W, chi=True, ret_C=True)
X = np.array([vt[l] * np.trace(C) ** 2 - 2 * (C ** 2).sum() for l, C in enumerate(Cs)])
print("excess/total:", np.round(X / (vt * np.array([np.trace(C) ** 2 for C in Cs])), 3))
for scale in (0.5, 1.0):
    p = tc.predict(W, chi=True, X4=scale * X)
    print(f"scale {scale}: TC {((p0[-1]-T[-1])**2).mean()-nz:.3e} -> TC+k4 {((p[-1]-T[-1])**2).mean()-nz:.3e}  bias {np.mean(p[-1]-T[-1]):+.2e} (TC bias {np.mean(p0[-1]-T[-1]):+.2e})")
