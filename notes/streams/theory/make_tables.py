"""Markdown tables for REPORT.md from the closure2 pair-run pickles.

    python make_tables.py results/width128_closure2_lr_mlp0.txt.pkl results/width128_closure2_lr_mlp1.txt.pkl
"""
import pickle
import sys

import numpy as np

COLS = [("noise", "noise"), ("wick", "Wick"), ("cl1_oracle", "1st, oracle coef (B3=1)"), ("cl1", "1st, B3=0.5"),
        ("cl1+B7", "+B7"), ("cl2", "+2nd order (fixed)"), ("cl2+k56", "+kappa5/6"), ("LR1", "leaf-resummed 1st"),
        ("LR2", "leaf-resummed 2nd"), ("LR2+k56", "LR2+kappa5/6"), ("LRT1", "leaf+T resummed 1st"), ("LRT2", "leaf+T resummed 2nd"), ("fit1", "fitted 1st (oracle fit)"),
        ("fit2", "refit 1st, 2nd fixed"), ("fitLR", "refit LR basis"), ("fitLRT", "refit LRT basis"), ("fitall", "all terms free")]
FIRST = ["B0", "B1", "B2", "B3", "B4", "B5", "B6"]
LEG = [1.0, 3.0, 3.0, 0.5, 1.0, 1.5, 1.5]


def load(paths):
    rows = {}
    for p in paths:
        for l, s, avg in pickle.load(open(p, "rb")):
            rows[l] = (s, avg)
    return [(l, *rows[l]) for l in sorted(rows)]


def eps_table(rows, title):
    have = [c for c in COLS if c[0] in rows[0][2]]
    out = [f"**{title}**", "", "| layer | " + " | ".join(h for _, h in have) + " |", "|---" * (len(have) + 1) + "|"]
    for l, s, avg in rows:
        out.append(f"| {l} | " + " | ".join(f"{100 * avg[k]:.2f}" for k, _ in have) + " |")
    # layer groups
    for lo, hi in ((1, 3), (4, 9), (10, 14)):
        sel = [avg for l, s, avg in rows if lo <= l <= hi]
        if sel:
            out.append(f"| {lo}-{hi} mean | " + " | ".join(f"{100 * np.mean([a[k] for a in sel]):.2f}" for k, _ in have) + " |")
    return "\n".join(out)


def coef_table(rows, title, key):
    out = [f"**{title}**", "", "| layer | " + " | ".join(f"{t} ({c:g})" for t, c in zip(FIRST[1:], LEG[1:])) + " |", "|---" * 7 + "|"]
    for l, s, avg in rows:
        c = np.mean([s[t][key] for t in ("AB", "BA")], axis=0)
        out.append(f"| {l} | " + " | ".join(f"{x:.2f}" for x in c[1:]) + " |")
    return "\n".join(out)


def lr_coef_table(rows, title):
    keys = ["LEAF", "TWOLEAF", "B6", "B7", "G3t"]
    theory = [6.0, -1.0, 1.5, 3.0, 1.0]
    out = [f"**{title}**", "", "| layer | " + " | ".join(f"{k} ({t:g})" for k, t in zip(keys, theory)) + " |", "|---" * 6 + "|"]
    for l, s, avg in rows:
        if "coefLR" not in s["AB"]:
            continue
        c = [np.mean([s[t]["coefLR"][k] for t in ("AB", "BA")]) for k in keys]
        out.append(f"| {l} | " + " | ".join(f"{x:.2f}" for x in c) + " |")
    return "\n".join(out)


def lrt_coef_table(rows, title):
    keys = ["LEAF", "TWOLEAF", "TCL", "B6", "G3t"]
    theory = [6.0, -1.0, 1.0, 1.5, 1.0]
    out = [f"**{title}**", "", "| layer | " + " | ".join(f"{k} ({t:g})" for k, t in zip(keys, theory)) + " |", "|---" * 6 + "|"]
    for l, s, avg in rows:
        if "coefLRT" not in s["AB"]:
            continue
        c = [np.mean([s[t]["coefLRT"][k] for t in ("AB", "BA")]) for k in keys]
        out.append(f"| {l} | " + " | ".join(f"{x:.2f}" for x in c) + " |")
    return "\n".join(out)


def size_table(rows, title, terms):
    out = [f"**{title}**", "", "| layer | " + " | ".join(terms) + " |", "|---" * (len(terms) + 1) + "|"]
    for l, s, avg in rows:
        out.append(f"| {l} | " + " | ".join(f"{100 * s['AB']['size'].get(t, float('nan')):.2f}" for t in terms) + " |")
    return "\n".join(out)


def scales_table(rows, title):
    keys = list(rows[0][1]["AB"]["scales"].keys())
    out = [f"**{title}**", "", "| layer | " + " | ".join(keys) + " |", "|---" * (len(keys) + 1) + "|"]
    for l, s, avg in rows:
        out.append(f"| {l} | " + " | ".join(f"{s['AB']['scales'][k]:.3f}" for k in keys) + " |")
    return "\n".join(out)


if __name__ == "__main__":
    groups = [g.split(",") for g in sys.argv[1:]]
    for g in groups:
        p = ",".join(g)
        rows = load(g)
        print(eps_table(rows, p)); print()
        print(coef_table(rows, "fit1 coefficients", "coef1")); print()
        print(coef_table(rows, "fit2 coefficients", "coef2")); print()
        print(lr_coef_table(rows, "LR refit coefficients")); print()
        print(lrt_coef_table(rows, "LRT refit coefficients")); print()
        print(scales_table(rows, "scales")); print()
