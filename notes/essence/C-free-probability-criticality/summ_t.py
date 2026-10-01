"""Summarise results/t1_<set>.jsonl: raw (mean ± se), products, ratio to A3gsl_nc, adjusted at Strassen L3 (0.683 u/product, B = 1024 u)."""
import sys, json, numpy as np, collections
d = collections.defaultdict(dict); P = collections.defaultdict(dict)
for l in open(sys.argv[1]):
    r = json.loads(l); d[r['variant']][r['mlp']] = r['raw']; P[r['variant']][r['mlp']] = r['nprod']
base = d['A3gsl_nc']; N = len(base)
print(f"{'variant':16s} {'raw mean':>10s} {'se':>9s} {'products':>8s} {'vs A3gsl_nc':>12s} {'range':>11s} {'adj(L3)':>9s} wins")
for v in sorted(d, key=lambda v: np.mean(list(d[v].values()))):
    x = d[v]
    if len(x) < N: continue
    a = np.array([x[i] for i in sorted(x)]); b = np.array([base[i] for i in sorted(x)]); p = np.mean(list(P[v].values()))
    print(f"{v:16s} {a.mean():10.3e} {a.std(ddof=1)/np.sqrt(len(a)):9.2e} {p:8.0f} {np.mean(a/b):12.2f} {min(a/b):5.2f}-{max(a/b):4.2f} {a.mean()*max(0.1, p*0.683/1024):9.2e} {np.sum(a<b)}/{len(a)}")
