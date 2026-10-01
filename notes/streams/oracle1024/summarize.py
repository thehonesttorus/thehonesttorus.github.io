"""Compact ladder tables from the pair json files written by `stream_oracle.py ladder A B --out X.json`.

    python summarize.py results/width1024_mlp770000_pair.json [...]
prints, per layer, the noise-free replica cross-product eps (sqrt of eps^2, negative if eps^2 < 0) with the
jackknife standard error of eps^2 propagated to eps, for the ladder columns; plus the D21 noise of one replica,
the replica noise of the closure's model error and the D21-space fitted coefficients (replica A)."""
import json, sys
import numpy as np

COLS = [("memless", "slices"), ("wick", "Wick"), ("herm4", "herm rho^4"), ("herm3", "herm rho^3"), ("closure_noK4", "clos-noK4"),
        ("closure", "closure"), ("closure_reg", "clos-uC"), ("fit", "fit"), ("fit_reg", "fit-uC"),
        ("closure_noB0", "clos-noB0"), ("closure_reg_noB0", "clos-uC-noB0"), ("closure_projB0", "clos-projB0"), ("closure_reg_projB0", "clos-uC-projB0")]


def eps(e2, se):
    e = np.sign(e2) * np.sqrt(abs(e2))
    de = se / (2 * max(np.sqrt(abs(e2)), 1e-12))
    return e, min(de, np.sqrt(se))


def table(path):
    d = json.load(open(path))
    rows = d["pair"]
    cols = [(c, h) for c, h in COLS if c + "_cp2" in rows[0]]
    out = [f"## {path}", "", "| l | D21 noise | closure δ-noise | B0 share | B0 proj R² | " + " | ".join(h for _, h in cols) + " | fit coef B0..B6 |",
           "|---" * (len(cols) + 6) + "|"]
    for r in rows:
        cells = []
        for c, _ in cols:
            e, de = eps(r[c + "_cp2"], r[c + "_cp2_se"])
            cells.append(f"{100 * e:.2f} ± {100 * de:.2f}")
        extra = (f"{r.get('B0_share_cp', r['B0_share']):.3f} | {r.get('projR2cp', r['projR2']):.3f} | ") if "B0_share" in r else "– | – | "
        out.append(f"| {r['l']} | {r['noise']:.3f} | {r['closure_dn']:.3f} | " + extra + " | ".join(cells) + " | "
                   + " ".join(f"{c:+.2f}" for c in r["coefA"]) + " |")
    return "\n".join(out)


if __name__ == "__main__":
    print("eps in % of ||D21(l+1)||, replica cross-product estimate ± jackknife s.e.\n")
    for p in sys.argv[1:]:
        print(table(p)); print()
