"""n = 1024 confirmation of the response-matched merge (ALS fitter), one network."""
import json, sys, time
import numpy as np
import rmem, bench
mlp = int(sys.argv[1]) if len(sys.argv) > 1 else 0
S = bench.load_set("w1024_d16")
t_ = time.time()
W, st = rmem.closure_states(bench.weights(S, mlp))
t0, age, tmax = 8, 3, 14
Y0, Z0, w, sid = rmem.legs_at_cut(W, st, t0, t0 - age)
Qs = rmem.future_props(W, st, t0, tmax)
print(f"closure+legs {time.time()-t_:.0f}s, atoms {Y0.shape[0]}", flush=True)
ranks = []
for k, t in enumerate(range(t0, tmax + 1)):
    if rmem.KOFF[t] == 0: continue
    M = rmem.readout(Y0 @ Qs[k], Z0 @ Qs[k], w, st[t]["e"], st[t]["P"])
    sv = np.linalg.svd(M, compute_uv=False); en = np.cumsum(sv ** 2) / np.sum(sv ** 2)
    ranks.append((t, int(np.searchsorted(en, 0.9) + 1), int(np.searchsorted(en, 0.99) + 1)))
print("per-target r90/r99:", ranks, flush=True)
out = dict(mlp=mlp, atoms=int(Y0.shape[0]), ranks=ranks, fits=[])
for R in [128, 256]:
    for fw in [2, 99]:
        t_ = time.time()
        r = rmem.x3_als(W, st, t0, tmax, rmem.KOFF, Y0, Z0, w, Qs, R, sweeps=2, cg_iters=20, gd_iters=25, fitwin=fw)
        r["sec"] = round(time.time() - t_, 1)
        h = r["hist"]
        print(f"R={R} fitwin={fw}: " + " | ".join(f"{x['stage']}{x['sweep']} fit {x['fit']:.4f} all {x['all']:.4f} p{x['passes']}" for x in h), f"[{r['sec']}s]", flush=True)
        out["fits"].append(r)
        json.dump(out, open(f"results/w1024_mlp{mlp}_als.json", "w"), indent=1)
