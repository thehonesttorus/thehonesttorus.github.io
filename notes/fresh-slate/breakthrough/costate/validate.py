"""Check the pull-back implementation against the full-tensor HD (star diagrams + exact coincidences) at n = 64."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../bench")); sys.path.insert(0, os.path.join(HERE, "../../designs/heisenberg"))
import bench, hd, costate
S = bench.load_set("w64_d16")
for i in range(2):
    W = bench.weights(S, i).astype(np.float64)
    a = hd.hd(W, A=None, diagrams="star", K=8)
    b = costate.predict(W)
    t = S["means"][i] if "means" in S else None
    print(i, "max|full - pullback| =", np.abs(a - b).max(), " final rms corr =", np.sqrt(np.mean((a[-1] - hd.closure(W)[-1]) ** 2)))
