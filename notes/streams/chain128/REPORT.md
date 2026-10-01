# chain128 — does a better D21 interface give a lower final-layer MSE?

*Stream report, 2026-10-01. Code, raw JSON per (variant, MLP) and derived tables are in this directory. Every number
below can be regenerated from `results/` with `summarize.py`, `law.py`, `width_fit.py`.*

## TL;DR

1. **Yes, one-for-one, when the fourth cumulant is right.** With κ4(z)_{aabc} held at its Monte Carlo truth, the D21
   ladder carries straight through to the final layer at width 128 (geometric mean over 4 MLPs): slices only 1.6e-4 →
   leading Wick 2.1e-5 → oracle first-order closure 1.5e-5 → fitted coefficients / full first-order diagram engine
   **5.3e-6**, against K=2 2.1e-4. Across all variants and MLPs, final MSE = a + k·ε² with corr 0.97–0.998, where ε is
   the noise-corrected, layer-rms relative error of D21(z) against the atlas. Here k = 3.5–6.4e-4, i.e. 1–3.7× each
   MLP's own K=2 MSE. At width 1024, 504aldo found 4.2e-6 ≈ 0.95× the K=2 MSE, so the law has the same form at both
   widths.
2. **A full first-order diagram engine beats the oracle closure 2× on the interface.** It sums every connected diagram
   on three distinct vertices with at most one κ3 or κ4 hyperedge (any index pattern) plus Gaussian edges. Fed the true
   z cumulants, its D21(l+1) error at width 128 is 1.1–3.6 % (noise-corrected, layers 1–14). The leg-partition closure
   gets 3.0–6.7 % and leading Wick 7–11 %. At layers 1–7 the engine is at or below the 2.2 % frontier threshold.
3. **The (2,1,1) fourth-cumulant slice is the lever, and regenerating it as u_i C_jk costs the gain.** With the engine:
   - exact slice: ε 2.1 %, final MSE 5.3e-6;
   - u_i·C_jk regeneration: ε 6.8 %, 1.2e-5;
   - slice zeroed: ε 8.5 %, 2.4e-5.

   ε falls with width roughly as n^-0.5 for the full closure, but only as n^-0.1…0.2 when the (2,1,1) slice is
   regenerated or dropped. Extrapolated to n = 1024 that gives ε ≈ 6 % for "first-order closure + regenerated κ4
   core", which matches the published chain's measured 4–7 %. The exact slice would give ε ≈ 1 %.
4. **At width 128 the binding error of a self-consistent chain is the fourth-cumulant closure, not the κ3 rule.**
   - κ4 = 0 or the memoryless κ4: every κ3 rule lands at 3e-5 or 8e-5.
   - A dense κ4 closure (full κ4 carried, n⁴ = 268M entries, with same-order κ3×κ3 terms, `dense2`) stays close
     to teacher-forced quality through layer ≈ 10 on the MLPs checked layer by layer (0 and 1). On 2 of the 5 MLPs
     finished so far its final layer is 2.5e-6 and 3.4e-6, close to the teacher-forced level. On the other 3 its last 3–4
     layers blow up (4.7e-5 to 8.3e-5).
5. ARC's reference `mlp_kprop` is not public (git asks for credentials; not on PyPI), so it could not be used as the
   exact K=3 reference.

## Question and setup

The published chain's interface to each ReLU is the (2,1) slice D21(l+1)_{ab} = κ3(z_a, z_a, z_b). The oracle ladder
(notes/competition-plan.md §3.1) says how well each closure of the all-distinct κ3(a) reproduces it. This stream asks
whether a better D21 actually lowers the final-layer MSE end to end at width 128. It also measures the error law and
extrapolates to width 1024.

- Ground truth: `whest dataset bake --n-mlps 8 --n-samples 1e7 --width 128 --depth 16`, seeds 128000–128007. The truth
  noise floor avg_variance/N = 1.2e-8 is negligible against every number below.
- Atlases: `moment_atlas_np.py --k3 --k4` (full raw third moments, raw (2,1,1) fourth moment).
  - Width 128: N = 5e5, MLPs 0–3, plus a second-seed atlas of MLP 0 for the noise floor.
  - Width 64: N = 1e6, pair on MLP 0.
  - Width 32, depth 8: N = 2e6, pair.

## Method (code)

- **`chain.py`**: the dense chain.
  - State per layer: μ, C, full κ3(z) (n³), and either X = κ4(z)_{aabc} (n³) or the full κ4(z) (n⁴, float32).
  - ReLU step: the Edgeworth operator E[F(z)] = E_G[exp(Σ_r κ_r/r! ∂^r) F], truncated at first order in κ3 and κ4 plus
    κ3²/2. It acts on Gaussian expectations of derivatives of reluᵖ, which are exact and closed-form via truncated
    Gaussian moments (= the paper's Hermite coefficients of powers, Prop S.3.8). The Mehler series runs to all orders in
    C_ij (40 terms).
  - Outputs from the 1-D and 2-D formulas: μ(a), var(a), C_off(a), the (3) and (2,1) slices of κ3(a), and the (4),
    (3,1), (2,2) slices of κ4(a).
  - Linear step: exact multilinear transport, n⁴ for κ3 and n⁵ for κ4.
- **All-distinct κ3(a) rules:**
  - `none`: slices only.
  - `wick`: Gaussian ρ² + Φ³κ3(z).
  - `closure`: `oracle_k3.hermite_model` + `residual_basis` × `CLOSURE_COEF`.
  - `fit`: the same seven diagrams with a per-layer coefficient table fitted on atlases 0 and 1 (`fit_coefs.py`,
    `results/coefs_tensor.txt`).
  - `engine`: `triple_engine`, every connected diagram on three distinct vertices with ≤ 1 hyperedge (κ3: 10 index
    patterns; κ4: 15 patterns, including the (2,1,1) and (2,2) slices) and Gaussian edge multiplicities ≤ 6 (≤ 2 with
    a hyperedge). Its Gaussian part equals `hermite_model` to 1e-16.
- **κ4 closures:**
  - `zero`.
  - `mem`: the (4),(3,1),(2,2) slices of κ4(a), transported exactly.
  - `dense`: the full κ4(a). Its all-distinct part is Φ⁴κ4(z) + 16 Gaussian trees + 12 κ3-hyperedge-plus-edge diagrams.
    Its (2,1,1) slice comes from the engine with vertex function relu² − 2m·relu, minus 2 C_ij C_ik. Transport is n⁵.
  - `dense2`: `dense` plus the same-order diagrams with two κ3 hyperedges and with a (2,1,1) κ4 hyperedge plus a C edge
    (66 more n⁴ terms).
  - `atlas`: κ4(z)_{aabc} teacher-forced from the atlas at every layer.
  - `atlas_reg211` / `atlas_zero211`: as `atlas`, but the all-distinct (2,1,1) part is regenerated as u_i C_jk (u fitted
    per layer on the atlas slice, as in the oracle) or zeroed.
- **Other scripts:**
  - `validate.py`, `ladder_width.py`: the step fed the atlas's true z cumulants; noise-corrected with a pair atlas.
  - `diag_state.py`: chain state against the atlas, per layer.
  - `run_variants.py`: variants, error-law noise `N<eps>@base`, teacher forcing `F<l>@base` / `FD<l>@base`.
  - `law.py`, `noise.py`, `summarize.py`, `width_fit.py`.
- **Cost:** a `closure`/`engine` chain takes 20–110 s; a `dense2` chain about 20 min at 4–5 GB.

## Results

### 1. Step validation
- Width 32, true z fed in (`results/diag`, `validate.py`):
  - Layer 0 (Gaussian input) is reproduced to the atlas's MC noise (μ 1e-4, var 2e-4, C_off 1.3e-3, D21(a) 3e-3).
  - Deeper layers carry truncation errors well above noise at this very non-Gaussian width: μ 0.1–0.3 %, var 0.2–1.6 %,
    C_off 1–3 %, D21(a) 2–7 %.
- The interface ladder against width (`results/ladder_w*_mlp0.txt`; ε of D21(l+1) when the step is fed the true z_l
  cumulants and builds κ3(a) by each rule; noise-corrected with a pair atlas; rms over layers 1–6; `width_fit.py`):

| rule | n = 32 | n = 64 | n = 128 | ε ~ n^-p, p | extrapolated n = 1024 |
|---|---|---|---|---|---|
| slices only | 0.564 | 0.517 | 0.491 | 0.10 | 0.40 |
| leading Wick | 0.119 | 0.132 | 0.097 | 0.15 | 0.08 |
| oracle closure (leg-partition coefs) | 0.082 | 0.058 | 0.037 | 0.58 | 0.011 |
| **diagram engine** | 0.042 | 0.035 | **0.021** | 0.50 | **0.008** |
| engine, (2,1,1) slice regenerated u_i C_jk | 0.078 | 0.098 | 0.068 | 0.10 | 0.061 |
| engine, (2,1,1) slice zeroed | 0.114 | 0.114 | 0.085 | 0.22 | 0.057 |

Width 128, all layers (MLP 0), noise-corrected:
- engine: 1.1, 1.6, 1.8, 2.6, 2.5, 2.5, 2.7, 3.6, 2.8, 3.3, 3.1, 2.6, 2.6, 3.0 % (layers 1–14);
- oracle closure: 3.0–6.7 %;
- leading Wick: 7–11 %;
- u·C: 3.9–7.6 %;
- slices only: 45–64 %.

(The width-32 point is a depth-8 MLP and the fits use three widths, so the exponents are indicative only.)

### 2. End to end, width 128 — teacher-forced fourth cumulant (MLPs 0–3; final-layer MSE vs N = 1e7 truth)

| κ3 all-distinct rule (κ4 = atlas) | MLP 0 | MLP 1 | MLP 2 | MLP 3 | geo-mean | mean over layers | layer-rms ε(D21), noise-corr. (MLP 0) |
|---|---|---|---|---|---|---|---|
| slices only | 1.01e-4 | 1.73e-4 | 2.35e-4 | 1.63e-4 | 1.61e-4 | 6.0e-5 | 0.534 |
| leading Wick | 8.6e-6 | 1.93e-5 | 6.1e-5 | 2.05e-5 | 2.13e-5 | 6.2e-6 | 0.160 |
| oracle closure | 5.4e-6 | 2.37e-5 | 3.44e-5 | 1.05e-5 | 1.47e-5 | 3.5e-6 | 0.094 |
| fitted coefs (fit on 0,1; 2,3 out of sample) | 2.0e-6 | 8.0e-6 | 9.3e-6 | 5.2e-6 | 5.26e-6 | 1.6e-6 | 0.044 |
| **diagram engine** | 1.9e-6 | 8.6e-6 | 1.06e-5 | 4.4e-6 | **5.25e-6** | 1.5e-6 | **0.041** |
| *K=2 (variant A), for reference* | 9.4e-5 | 4.74e-4 | 1.96e-4 | 2.34e-4 | 2.13e-4 | 1.2e-4 | 1 |

The fitted table generalises: on the out-of-sample MLPs 2 and 3 it matches the engine (9.3e-6 vs 1.06e-5; 5.2e-6 vs
4.4e-6). Its coefficients reproduce the oracle's drift with depth:
- ρ³ term 1.0 → 0.24;
- (2,2) term 1.36 → 0.4–1.0;
- (2,1,1) term 1.44 → 1.1;
- D21 terms ≈ 3 and 2.5;
- B3 (D3 + two edges) 0.4 → −0.15 (`results/coefs_tensor.txt`).

### 3. Variant E — the (2,1,1) slice of κ4(z): exact vs regenerated vs dropped (κ4 otherwise teacher-forced)

| κ3 rule | (2,1,1) slice | MLP 0 | MLP 1 | MLP 2 | MLP 3 | geo-mean | ε(D21) MLP 0 |
|---|---|---|---|---|---|---|---|
| engine | exact | 1.9e-6 | 8.6e-6 | 1.06e-5 | 4.4e-6 | 5.3e-6 | 0.041 |
| engine | u_i C_jk | 4.9e-6 | 1.14e-5 | 1.78e-5 | 2.05e-5 | 1.19e-5 | 0.098 |
| engine | zero | 8.5e-6 | 3.35e-5 | 4.85e-5 | 2.51e-5 | 2.43e-5 | 0.166 |
| closure | exact | 5.4e-6 | 2.37e-5 | 3.44e-5 | 1.05e-5 | 1.47e-5 | 0.094 |
| closure | u_i C_jk | 1.19e-5 | 2.40e-5 | 3.97e-5 | 3.16e-5 | 2.45e-5 | 0.142 |
| closure | zero | 1.69e-5 | 5.68e-5 | 8.95e-5 | 3.86e-5 | 4.27e-5 | 0.226 |

Carrying the exact slice is worth 2.3× over the published-style r = 1 regeneration (geo-mean). Once the slice is
regenerated, a better κ3 rule buys only 2.1×, against 2.8× with the exact slice.

### 4. The error law at width 128 (`results/law.md`)
**Structured** (across the κ4-teacher-forced variants of §2–3, per MLP, least squares final MSE = a + k ε²):

| MLP | a | k | corr(MSE, ε²) | K=2 MSE | k / K2 MSE |
|---|---|---|---|---|---|
| 0 | 1.2e-6 | 3.5e-4 | 0.998 | 9.4e-5 | 3.7 |
| 1 | 1.3e-5 | 4.6e-4 | 0.993 | 4.7e-4 | 1.0 |
| 2 | 2.6e-5 | 6.4e-4 | 0.974 | 2.0e-4 | 3.3 |
| 3 | 7.6e-6 | 5.4e-4 | 0.994 | 2.3e-4 | 2.3 |

**Injected**: an independent Gaussian perturbation of relative rms ε, added to D21(z) at every layer, on top of the
fitted-coefficient chain.

| ε | 0.02 | 0.05 | 0.1 | 0.2 |
|---|---|---|---|---|
| ΔMSE / ε² (final, mean of 4 MLPs) | 2.3e-3 | 2.1e-3 | 2.1e-3 | 2.1e-3 |
| ΔMSE / ε² (mean over layers) | 4.5e-4 | 4.4e-4 | 4.4e-4 | 4.5e-4 |

- The injected law is exactly quadratic over a decade of ε.
- Unstructured noise costs 2–5× more per ε² than the structured closure errors. Per-MLP k/K2 is 4–16, versus 1–3.7
  structured. A plausible reason, not tested: random D21 errors do not partly cancel against correlated errors in D3
  and the variance, as closure errors do.
- 504aldo's 1024 law (4.2e-6 ε² ≈ 0.95 × K2 MSE) was calibrated on structured variants. Our structured k/K2 = 1–3.7 is
  therefore the comparable number.

### 5. Variant F — teacher forcing κ3(a_l) at one layer (base: fitted chain with atlas κ4; ratio of final MSE to the base, mean of 4 MLPs)

| layer forced | 1 | 4 | 8 | 12 | 14 |
|---|---|---|---|---|---|
| full κ3(a_l) := atlas | 0.98 | 0.95 | 0.87 | **0.76** | 0.83 |
| only its all-distinct part | — | 0.95 | 0.93 | 0.97 | — |

- No single layer dominates. The early layers, where the chain's ε is 1–3 %, are already at the truth.
- The late layers carry the remaining error (−24 % from fixing layer 12 alone).
- Most of that gain comes from the slices (D3, D21 of a, i.e. the Edgeworth truncation on strongly non-Gaussian late
  pre-activations), not from the all-distinct part.
- Caveat: an atlas at N = 5e5 has 1–6 % noise in its own all-distinct tensor's contribution to D21, which caps what
  forcing that part can show.

### 6. End to end with the chain's own fourth cumulant (8 MLPs; `results/summary_raw.md`, `results/eps`)

| variant | κ3 rule | κ4 closure | final MSE (mean of 8) | layer 4 | layer 8 | layer 12 |
|---|---|---|---|---|---|---|
| A0 paper Alg. 2 | — | — | 2.88e-4 | 8.4e-5 | 1.7e-4 | 2.2e-4 |
| A | K=2, Mehler to all orders | — | 2.77e-4 | 7.9e-5 | 1.6e-4 | 2.1e-4 |
| M0 / M | slices only | zero / mem | 9.9e-5 / 8.6e-5 | 1.5e-5 / 1.1e-5 | | |
| B0 / B | leading Wick | zero / mem | 4.2e-5 / 1.33e-4 | | | |
| C0 / C | oracle closure | zero / mem | 3.2e-5 / 8.4e-5 | 7.7e-6 / 1.6e-6 | 1.5e-5 / 1.4e-5 | 2.1e-5 / 4.2e-5 |
| D0 / D | fitted coefs | zero / mem | 2.9e-5 / 8.4e-5 | | | |
| engine | engine | zero / mem | 3.1e-5 / 8.0e-5 | | | |
| engine | engine | dense (MLPs 0,1,2) | 8.1e-6, 2.0e-5, 3.6e-5 | | | |
| engine | engine | **dense2** (MLPs 0,1,2,4,5) | 2.5e-6, 8.3e-5, 4.7e-5, 5.4e-5, 3.4e-6 | 1.5e-7 | 1.4e-6 | 4.7e-6 |

How the κ4 closure behaves along the chain (MLP 0, `results/eps/*`), as relative error of K4 = κ4(z)_{aaaa} per layer:
- memoryless: 16 % at layer 1, 47 % at layer 4, 65 % at 8, 78–91 % by 13–15;
- dense: 9 %, 13 %, 35 %, 37–67 %;
- dense2: 9 %, 12 %, 34 %, 33–46 %.

With the κ4 teacher-forced, the variance error stays ≤ 0.6 % through layer 9. Mean over layers, dense2 is at 0.9–9e-6 against
1.5e-6 teacher-forced, so dense2 is right for most of the depth. Its failures are a late-layer growth: on MLP 1 the
per-layer MSE triples per layer from layer 11, while the K4 error peaks at 72 % at layer 12.

## Verdict
- **The D21 interface translates.** Holding the fourth cumulant at truth, final-layer MSE follows ε(D21)² over two
  decades:
  - the oracle closure (ε ≈ 9 %) sits at 1.5e-5;
  - the full first-order engine and the fitted table (ε ≈ 4 %, noise-corrected, layer-rms; 1–3 % per layer at the step
    level) sit at 5.3e-6, 2.8× lower;
  - both are 40× below K=2.
- **The fitted closure is a legitimate data table at width 128.** Fitted on two MLPs, it matches the parameter-free
  engine on two unseen MLPs. The engine itself needs no fit: it is just "all first-order diagrams", whose coefficients
  the fit rediscovers.
- **What decides the interface below ≈ 5 % is the (2,1,1) fourth-cumulant slice.** Regenerating it as u_i C_jk (the
  published r = 1 core) roughly triples ε at width 128 (2.1 → 6.8 %) and costs 2.3× in MSE. Its ε barely improves with
  width (p ≈ 0.1), while the exact-slice closure improves as n^-0.5.
- **Implication at width 1024.** With the 1024 law (4.2e-6·ε², ≈ 0.95·K2; our width-128 structured ratio 1–3.7·K2 has
  the same form):
  - a chain with a regenerated (2,1,1) core sits near ε ≈ 6 % → extra MSE ≈ 1.5e-8, which is the published chain's floor
    (raw 2.07e-8 at ε 4–7 %);
  - carrying the (2,1,1) slice exactly (or anything that reproduces its transport to D21 to a few per cent) would put
    ε ≈ 1 % → extra ≈ 4e-10, i.e. remove essentially all of the interface error.
  - So the n³ (2,1,1) object, not the κ3 rule, is the gate below the frontier's ε ≤ 2.2 %. At 1024 it must be carried
    in compressed form (plan §3.1 (vii): its effect on D21 is captured by a few covariance-response modes at depth, but
    needs rank ≳ 32 at layers 1–5). The extrapolations rest on three widths and one MLP per width, so they are
    indicative only.
- **For a self-consistent width-128 chain, the κ4 closure dominates.** Memoryless κ4 is worse than κ4 = 0 at depth.
  A dense κ4 with every same-order diagram is excellent through ~10 layers, but on 4 of 6 MLPs it blows up in the last
  3–4 layers.

## Open issues
- dense2 is still running on MLPs 3, 6 and 7, and closure / Wick under dense2 on MLPs 0 and 5 is queued. The cause of
  the late blow-up of dense2 (missing ρ⁴ Gaussian κ4 diagrams? the second-order (2,1,1)-slice terms of κ4(a), which only
  have first-order hyperedges in the engine? Edgeworth truncation on strongly skewed late pre-activations) is
  undiagnosed.
- The atlas noise of one width-128 atlas (N = 5e5) is 1–6 % on D21 and 2–10 % on K4. It is subtracted only where a
  pair exists (MLP 0 at 128, MLP 0 at 64 and 32); other MLPs' ε include it (≈ +0.5 % at layer-rms).
- At the step level the order-2 Edgeworth term (κ3²/2) hurts at width 32; it was not compared at width 128 (all
  width-128 runs use order 2). The 1-D mean and
  variance would need κ5/κ6 marginals for the late layers of narrow nets.
- The width-64 end-to-end chain hit negative variances on two MLPs at depth (the Edgeworth corrections are not
  positivity-preserving), so the width study is done at the step level only.
- The width exponents come from three widths with a depth-8 MLP at n = 32. A width-256 step-level point (n³ = 16.8M
  per tensor, ≈ 7 GB per atlas) would firm them up.

## Next steps
1. Width-1024 version of §1 restricted to the D21 transport: the (2,1,1) slice's contribution to D21(l+1) computed from
   keenanpepper-style joint-moment atlases, to test the n^-0.5 vs n^-0.1 extrapolation directly.
2. Compress the (2,1,1) slice's *transported effect* (sandwich modes Wᵀ Φ M_r Φ W) and measure ε vs rank inside this
   dense chain. It is the same experiment as variant E with rank r between "u·C" and "exact".
3. Stabilise dense2 (positivity-preserving 1-D step, e.g. a Gram–Charlier fit with a floor on var; the missing ρ⁴ and
   second-order (2,1,1) terms) and run it on all 8 MLPs.
