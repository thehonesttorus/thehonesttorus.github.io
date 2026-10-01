"""Where the coherence of the old atoms lives: Perron direction vs bulk; cross-source Gram (n = 1024, one MLP)."""
import sys, os, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '../../fresh-slate/bench'))
import bench, fcs
S = bench.load_set('w1024_d16'); i = int(sys.argv[1]) if len(sys.argv) > 1 else 0
W = bench.weights(S, i); d = []
fcs.run(W, slices=2, k4mf=True, diag=d, diag_fracs=(0.5,), diag_store={6, 10, 14})
out = []
for r in d:
    if "Dold" not in r: continue
    Do = r["Dold"].astype(np.float64); D = r["D"].astype(np.float64); v = r["perron"]
    Ds = [x.astype(np.float64) for x in r["Ds"]]; Es = np.array(r["Es"])
    G = np.array([[np.sum(a * b) for b in Ds] for a in Ds]); nrm = np.sqrt(np.diag(G))
    cos = G / np.outer(nrm, nrm)
    Pb = Do @ np.outer(v, v); Pa = np.outer(v, v) @ Do
    res = Do - Pb - Pa + np.outer(v, v) @ Do @ np.outer(v, v)
    o = dict(t=r["l"], nsrc=len(Ds), gamma_old=float((Do**2).sum() / r["E"]),
             sum_src_over_total=float(np.trace(G) / (Do**2).sum()),
             mean_offdiag_cos=float((cos.sum() - len(Ds)) / max(1, len(Ds) * (len(Ds) - 1))),
             cos_adjacent=[float(cos[k, k + 1]) for k in range(len(Ds) - 1)],
             frac_b_perron=float((Pb**2).sum() / (Do**2).sum()), frac_a_perron=float((Pa**2).sum() / (Do**2).sum()),
             frac_resid=float((res**2).sum() / (Do**2).sum()), old_over_D=float(np.sqrt((Do**2).sum() / (D**2).sum())),
             cos_old_young=float(np.sum(Do * (D - Do)) / np.sqrt((Do**2).sum() * ((D - Do)**2).sum())),
             gamma_src=[float(G[k, k] / Es[k]) for k in range(len(Ds))])
    out.append(o); print(json.dumps(o), flush=True)
json.dump(out, open(os.path.join(HERE, f'results/coh_m{i}.json'), 'w'), indent=1)
