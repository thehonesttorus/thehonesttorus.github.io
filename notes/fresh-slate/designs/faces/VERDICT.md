# Faces stream — verdict (1 Oct 2026)

**Principle.** Faces, barycentres and conditional expectations. Without biases every linear region is a cone. On a face the arrow (ReLU then $W$) is linear, and barycentres of past events transport exactly through it.

## Numbers (bench w1024_d16, 6 MLPs, truth noise 3.6e-8 subtracted)

| estimator | raw at 1024 ± s.e. | cost (units) | adjusted |
|---|---|---|---|
| Gaussian closure, face-measure Mehler (reference) | 4.10e-6 ± 0.33e-6 | ≈ 40 | 4.1e-7 |
| one-step facet tree + (2,1) slice ('mem21') | 2.59e-6 ± 0.14e-6 | ≈ 100 | 2.6e-7 |
| facet births, renormalised legs, depth window 6 ('win6') | 4.45e-7 ± 0.33e-7 | ≈ 410 | 1.8e-7 |
| facet births, renormalised legs, all depths ('w16') | **3.24e-7 ± 0.31e-7** | ≈ 630 (≈ 360 with Strassen) | 2.0e-7 (≈ 1.1e-7) |

The gain over Gaussian closure grows with width: ×2.2 at 64, ×7 at 128, ×10 at 256, ×12 at 512, ×12.7 at 1024. Not competitive: about 70× above the bar on adjusted score, and 10× behind the region stream's FC v2 (3.03e-8 raw).

## What the faces principle contributed that the chain lacked

1. **Non-Gaussianity is carried by facets (codimension 1), not by faces (codimension 0).** The gate field misses about 99 % of the skew, and it is exactly symmetric at p = 1/2 (T1). By E2, $\mathbb E f=\mathbb E\Delta f$ and barycentres are facet integrals. The vertex coefficients of any closure are face measures of increasing codimension: mass, facet density, and its derivatives.
2. **An exact telescoping (E5)**, $\tilde z_l = xP_{x\to l}+\sum_{s<l}\nu_sP_{s\to l}$, with births $\nu_s$ at the facets and face-averaged arrows $P=WD_\beta\cdots W$. It turns "old content" into an explicit object and gives an exact age attribution of every cumulant.
3. **Measured age profile.** The final-layer κ3 is born roughly uniformly over all 15 earlier layers plus the input. The most recent facets contribute about nothing net, so memoryless gluing is impossible in principle. This is the structural reason the full-history lines need all ages.
4. **Ceilings, each with its reason.**
   - Facet births on input-chaos legs carry only about 40 % of the variance at depth, independently of width.
   - Renormalised legs at the birth layer are needed; that is v3.
   - Transport by gate averages $\beta$ alone stalls at 3–4e-7. This matches the judge's diagnosis (mean-gate transport of coincident index patterns).
5. **Exact identities, all checked numerically** (`check_exact.py`):
   - the barycentre–facet identity;
   - Stein–Euler, $\mathbb Ef=\mathbb E\Delta f$;
   - the barycentric decomposition, and $\mathrm{Cov}(a)$ in face data;
   - radial factorisation, $a=R\,b(\theta)$;
   - degree-0 Gram transport, $G_{l+1}=W^\top DGDW$.

## Falsified, charged to the dictionary

- gate-field face state (v0);
- hub-neuron gluing: right shape, 2–3× too large at depth; not the radial mode;
- face-wise-constant Gram decoupling: 5–10× worse than Gaussian;
- births on input legs only (v1–v2);
- naive (2,2)/(3,1) pair terms;
- radial κ4 correction: null;
- curvature-passage term in κ3 only: ≤ 3 %.

The missing object is the coherence between facets of different depths (how a facet of layer $s$ conditions the joint law of layer $l$). Wiener-chaos slice transport carries it; mean-gate face transport does not.

Details, logs and code are in `DESIGN.md` §§7–9, `fbt.py`, `t1*`–`t5*`, and `results/`.
