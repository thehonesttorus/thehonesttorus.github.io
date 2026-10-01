"""Tables from agespec results: per-age energy, dilation share, conservation of gamma, free-law check, coherence."""
import sys, json, numpy as np
d = json.load(open(sys.argv[1]))
rows, layers = d["rows"], d["layers"]
n_layers = len(layers)
g = {x["layer"]: x.get("g") for x in layers}
print("layer  Phi2   g=2E[Phi^2]  t_rms")
for x in layers:
    print(f"{x['layer']:3d}  {x['Phi2']:.3f}  {x.get('g', float('nan')):.3f}  {x['t_rms']:.2f}")

print("\nfree law check: PR*age'/n and tr vs prod g (age' = number of factors = age+1)")
by = {}
for r in rows:
    by.setdefault(r["age"], []).append(r)
n = None
for a in sorted(by):
    rs = by[a]
    prodg = [np.prod([g[l] for l in range(r["s"] + 1, r["k"])]) if r["k"] > r["s"] + 1 else 1.0 for r in rs]
    tr = np.array([r["tr"] for r in rs]); PR = np.array([r["PR"] for r in rs])
    print(f"age {a:2d}  tr/prod(g) {np.mean(tr/np.array(prodg)):.3f}  PR {np.mean(PR):7.1f}")

print("\nper (s,k): E (pair D21 energy), dilation share, gamma")
E = np.full((n_layers, n_layers), np.nan); SH = E.copy(); GA = E.copy()
for r in rows:
    E[r["s"], r["k"]] = r["E"]; SH[r["s"], r["k"]] = r["share"]; GA[r["s"], r["k"]] = r["gamma"]
print("share (rows s, cols k):")
for s in range(n_layers - 1):
    print(f"s={s:2d} " + " ".join("  .  " if np.isnan(SH[s, k]) else f"{SH[s,k]:.2f} " for k in range(n_layers)))
print("gamma ratio gamma(s,k+1)/gamma(s,k) (conservation test; 1 = conserved charge):")
for s in range(n_layers - 2):
    print(f"s={s:2d} " + " ".join("  .  " if (np.isnan(GA[s, k]) or np.isnan(GA[s, k+1]) or GA[s, k] == 0) else f"{GA[s,k+1]/GA[s,k]:5.2f}" for k in range(n_layers - 1)))
print("E ratio E(s,k+1)/E(s,k) vs g_k^3:")
for s in range(n_layers - 2):
    print(f"s={s:2d} " + " ".join("  .  " if (np.isnan(E[s, k]) or np.isnan(E[s, k+1])) else f"{E[s,k+1]/E[s,k]/g[k]**3:5.2f}" for k in range(n_layers - 1)))

print("\nper layer: coherent vs incoherent sums (scale part, residual part) and cross-age residual cosines")
for x in layers:
    if "resgram" not in x: continue
    Gr = np.array(x["resgram"]); dg = np.sqrt(np.diag(Gr)); Cc = Gr / np.outer(dg, dg)
    off = Cc[~np.eye(len(dg), dtype=bool)]
    print(f"L{x['layer']:2d} ages {len(dg):2d} | scale coh/inc {x['Escale_sum']/max(x['Escale_inc'],1e-300):6.2f} | res coh/inc {x['Eres_sum']/x['Eres_inc']:5.2f} | res cos mean {off.mean() if off.size else 0:+.3f} max|.| {np.abs(off).max() if off.size else 0:.3f} | scale share of total {x['Escale_sum']/(x['Escale_sum']+x['Eres_sum']):.2f}")

print("\nper-step ratios by age (median over sources): residual energy / g^3, scale energy / g^3, gamma ratio")
Eres = E * (1 - SH); Esc = E * SH
for a in range(0, n_layers - 2):
    rr, ss, gg = [], [], []
    for s in range(n_layers):
        k = s + 1 + a
        if k + 1 < n_layers and not np.isnan(E[s, k]) and not np.isnan(E[s, k + 1]):
            rr.append(Eres[s, k + 1] / Eres[s, k] / g[k] ** 3); ss.append(Esc[s, k + 1] / Esc[s, k] / g[k] ** 3)
            gg.append(GA[s, k + 1] / GA[s, k])
    if rr:
        print(f"age {a:2d}->{a+1:2d}  res/g^3 {np.median(rr):5.2f}  scale/g^3 {np.median(ss):5.2f}  gamma ratio {np.median(gg):5.2f}  (n={len(rr)})")
