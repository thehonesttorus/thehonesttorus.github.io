import sys, os, json, numpy as np, collections
HERE = os.path.dirname(os.path.abspath(__file__))
D = "results_live" if os.path.isdir(os.path.join(HERE, "results_live")) else "results"
for fn in sorted(os.listdir(os.path.join(HERE, D))):
    if not fn.endswith(".jsonl") or fn.startswith(("oracle2nd_", "hybrid_")): continue
    rows = [json.loads(l) for l in open(os.path.join(HERE, D, fn))]
    by = collections.OrderedDict()
    for r in rows: by.setdefault(r["variant"], {})[r["mlp"]] = r
    print(f"## {fn[:-6]}")
    print("| variant | MLPs | raw final (mean ± s.e.) | gain vs gauss | n^3 products |")
    print("|---|---|---|---|---|")
    g = by.get("gauss", {})
    for v, d in by.items():
        raws = np.array([d[i]["raw"] for i in sorted(d)])
        common = [i for i in d if i in g]
        gain = np.mean([g[i]["raw"] for i in common]) / np.mean([d[i]["raw"] for i in common]) if common else float("nan")
        se = raws.std(ddof=1) / np.sqrt(len(raws)) if len(raws) > 1 else float("nan")
        print(f"| {v} | {len(raws)} | {raws.mean():.3e} ± {se:.1e} | {gain:.2f} | {d[sorted(d)[0]]['nprod']} |")
