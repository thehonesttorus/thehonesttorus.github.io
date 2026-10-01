"""Tables from fcdil_instr_<mlp>.json: dilation-sector share of old D21 by layer; age spectrum of total and remainder."""
import sys, os, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
mlps = [int(x) for x in sys.argv[1].split(",")]
I = {i: json.load(open(os.path.join(HERE, "results_live", f"fcdil_instr_{i}.json"))) for i in mlps}
print("### share of ||D21_old||_F^2 in the dilation sector, by target layer t (mean over MLPs", mlps, ")")
print("| t | age>2: K-template | (s2,mu) rank-1 off-diag | top-1 off-diag | top-8 off-diag | age>4: K-template | (s2,mu) | top-1 |")
print("|---|---|---|---|---|---|---|---|")
T = sorted({r["t"] for r in I[mlps[0]]})
for t in T:
    rs = [next(r for r in I[i] if r["t"] == t) for i in mlps]
    def g(key, f): 
        vals = [r[key][f] for r in rs if key in r]
        return f"{np.mean(vals):.2f}" if vals else "—"
    print(f"| {t} | {g('old2','share_K')} | {g('old2','share_s2mu_off')} | {g('old2','share_top1_off')} | {g('old2','share_top8_off')} | {g('old4','share_K')} | {g('old4','share_s2mu_off')} | {g('old4','share_top1_off')} |")
print("\n### age spectrum (mean over target layers t >= 8 and MLPs): ||D_s||^2, remainder after the dilation projection, dilation amplitude")
print("| age | ||D_age||_F^2 | remainder ||D - a K||^2 | remainder share | dilation amplitude a (×1e3) |")
print("|---|---|---|---|---|")
agg = {}
for i in mlps:
    for r in I[i]:
        if r["t"] < 8: continue
        for a, n2, rm, ka in zip(r["ages"], r["norm2"], r["rem2"], r["kamp"]):
            agg.setdefault(a, []).append((n2, rm, ka))
prev = None
for a in sorted(agg):
    x = np.array(agg[a]); n2, rm, ka = x[:, 0].mean(), x[:, 1].mean(), x[:, 2].mean()
    ratio = "" if prev is None else f" (rem ratio to previous age {rm / prev:.2f})"
    print(f"| {a} | {n2:.3e} | {rm:.3e} | {rm / n2:.2f} | {1e3 * ka:+.2f} |{ratio}")
    prev = rm
