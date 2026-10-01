import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "bench")); sys.path.insert(0, HERE)
import bench, fbt
S = bench.load_set("w64_d16")
orig_readout = fbt.readout
for mode in ["gauss", "mem"]:
    res = {"free": [], "oracle_mu": []}
    for i in range(len(S["seeds"])):
        W = bench.weights(S, i); T = S["means"][i]
        res["free"].append(((fbt.predict(W, mode=mode) - T) ** 2).mean(1))
        # oracle: replace each layer's mean (fed forward) by the truth, keep closure covariance
        layer = [0]
        def ro(mu, var, k3, k4, _T=T):
            v = orig_readout(mu, var, k3, k4); l = layer[0]; layer[0] += 1
            ro.pred.append(v.copy()); return _T[l].copy()
        ro.pred = []
        fbt.readout = ro
        fbt.predict(W, mode=mode)
        fbt.readout = orig_readout
        res["oracle_mu"].append(((np.array(ro.pred) - T) ** 2).mean(1))
    for k, v in res.items():
        v = np.mean(v, 0)
        print(mode, k, " ".join(f"{x:.1e}" for x in v[[0, 1, 3, 7, 11, 15]]))
