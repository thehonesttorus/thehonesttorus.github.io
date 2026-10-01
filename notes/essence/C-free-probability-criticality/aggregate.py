"""Aggregate the age-spectrum laws over networks: python aggregate.py results/agespec_w1024_d16_*.json"""
import sys, json, numpy as np
R = [json.load(open(f)) for f in sys.argv[1:]]
print(f"{len(R)} networks")
L = len(R[0]["layers"])
# 1. odds law: scale share of total D21 at layer k vs number of live ages
print("layer  ages  scale-share(mean)  odds/ages (mean, min, max)   scale coh/inc / ages   res coh/inc   res cos mean")
for k in range(1, L):
    sh, od, ci, rc, cm = [], [], [], [], []
    for d in R:
        x = d["layers"][k]
        s = x["Escale_sum"] / (x["Escale_sum"] + x["Eres_sum"]); sh.append(s); od.append(s / (1 - s) / k)
        ci.append(x["Escale_sum"] / x["Escale_inc"] / k); rc.append(x["Eres_sum"] / x["Eres_inc"])
        Gr = np.array(x["resgram"]); dg = np.sqrt(np.diag(Gr)); C = Gr / np.outer(dg, dg)
        cm.append(C[~np.eye(len(dg), dtype=bool)].mean() if k > 1 else 0.0)
    print(f"{k:4d} {k:5d}   {np.mean(sh):.3f}             {np.mean(od):.3f} ({min(od):.3f}, {max(od):.3f})        {np.mean(ci):.3f}            {np.mean(rc):.3f}        {np.mean(cm):+.3f}")
# 2. per-step transfer rates by age
print("\nage step   residual/g^3   scale/g^3   gamma ratio    (median over all sources and networks)")
for a in range(0, L - 2):
    rr, ss, gg = [], [], []
    for d in R:
        g = {x["layer"]: x.get("g") for x in d["layers"]}
        P = {(r["s"], r["k"]): r for r in d["rows"]}
        for s in range(L):
            k = s + 1 + a
            if (s, k) in P and (s, k + 1) in P:
                p, q = P[(s, k)], P[(s, k + 1)]
                rr.append(q["E"] * (1 - q["share"]) / (p["E"] * (1 - p["share"])) / g[k] ** 3)
                ss.append(q["E"] * q["share"] / (p["E"] * p["share"]) / g[k] ** 3)
                gg.append(q["gamma"] / p["gamma"])
    if rr:
        print(f"{a:2d}->{a+1:2d}      {np.median(rr):.3f}         {np.median(ss):.3f}       {np.median(gg):.3f}   (n={len(rr)})")
# 3. free-probability law for the propagators
print("\nage   PR/(n/(2(age+1)))   tr/(2 prod g)")
n = None
for a in range(0, L - 1):
    pr, tr = [], []
    for d in R:
        g = {x["layer"]: x.get("g") for x in d["layers"]}
        for r in d["rows"]:
            if r["age"] == a:
                nn = r["PR"]  # PR
                pr.append(r["PR"]); tr.append(r["tr"] / (2 * np.prod([g[l] for l in range(r["s"] + 1, r["k"])])))
    nw = 1024 if np.mean(pr) * (a + 1) > 400 else 128
    print(f"{a:3d}   {np.mean(pr) * 2 * (a + 1) / nw:.3f}               {np.mean(tr):.4f}")
# 4. free-sector tail at the final kink layer: residual energy in ages > A as a share of all residual and of total
print("\nfinal layer k = L-1: share of residual (free-sector) energy and of total D21 energy in ages > A")
for A in (0, 1, 2, 3, 5, 7):
    fr, ft = [], []
    for d in R:
        rows = [r for r in d["rows"] if r["k"] == L - 1]
        Er = np.array([r["E"] * (1 - r["share"]) for r in rows]); ages = np.array([r["age"] for r in rows])
        x = d["layers"][L - 1]; tot = x["Escale_sum"] + x["Eres_sum"]
        fr.append(Er[ages > A].sum() / Er.sum()); ft.append(Er[ages > A].sum() / tot)
    print(f"A={A}: residual tail {np.mean(fr):.3f} of free-sector energy, {np.mean(ft):.3f} of total D21 energy")
