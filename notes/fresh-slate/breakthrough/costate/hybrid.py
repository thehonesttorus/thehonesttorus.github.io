"""Computed co-state + true joint kappa_4 (MC atlas, rank-filtered) injected on top: how close to the oracle?"""
import sys, os, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../bench")); sys.path.insert(0, HERE)
import bench, costate
from oracle2nd import passes
SP = "/tmp/claude-0/-home-user-thehonesttorus-github-io/73934192-bb11-53f5-adb1-9561b059184f/scratchpad"


def filt(X, r, drop_top=False):
    d = np.diag(X).copy(); O = X - np.diag(d)
    U, s, Vt = np.linalg.svd(O)
    lo = 1 if drop_top else 0
    hi = len(s) if r is None else r
    return (U[:, lo:hi] * s[lo:hi]) @ Vt[lo:hi] + np.diag(d)


name, N, i = sys.argv[1], int(float(sys.argv[2])), int(sys.argv[3])
S = bench.load_set(name); Ws = bench.weights(S, i).astype(np.float64)
f = os.path.join(SP, f"atlas_{name}_{i}_{N}.npz")
if os.path.exists(f):
    st = dict(np.load(f))
else:
    p1 = passes(Ws, N, 1000 + i); st = passes(Ws, N, 2000 + i, mean=p1["m"]); np.savez(f, **st)
L = Ws.shape[0]
def k4(r, drop_top):
    out = {}
    for l in range(1, L):
        c2 = st["c2"][l]; v = np.diag(c2)
        K4d = st["d4"][l] - 3 * v * v; K22 = st["s22"][l] - np.outer(v, v) - 2 * c2 * c2; K31 = st["s31"][l] - 3 * v[:, None] * c2
        out[l] = (K4d, filt(K22, r, drop_top), filt(K31, r, drop_top))
    return out
base = dict(A=3, old="gsmslice", law=True, coinc=False, K=21)
outf = os.path.join(HERE, "results_live", f"hybrid_{name}.jsonl")
for lab, kw in [("A3gsl_nc", dict()), ("+k4 full", dict(k4=k4(None, False))), ("+k4 rank8", dict(k4=k4(8, False))),
                ("+k4 full, top mode dropped", dict(k4=k4(None, True))), ("+k4 rank2-8", dict(k4=k4(8, True))),
                ("A3 (no scale law) +k4 full", dict(k4=k4(None, False), _nolaw=True))]:
    kk = dict(base); kk.update({k: v for k, v in kw.items() if not k.startswith("_")})
    if kw.get("_nolaw"): kk.update(law=False, old="slice")
    pred = costate.predict(Ws, **kk)
    raw = float(((pred[-1] - S["means"][i][-1]) ** 2).mean() - S["noise"][i])
    open(outf, "a").write(json.dumps(dict(set=name, mlp=i, variant=lab, raw=raw)) + "\n")
    print(name, i, f"{lab:28s} raw {raw:.3e}", flush=True)
