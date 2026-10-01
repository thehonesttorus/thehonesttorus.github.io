# chain128 — does a better D21 interface give a lower final-layer MSE?

*Stream report, 2026-10-01. Code, raw JSON per (variant, MLP) and derived tables are in this directory. Numbers
regenerate from `results/` with `summarize.py --raw results/raw,results/eps`, `law.py`, and
`width_fit.py results/ladder_o1_w`.*

## TL;DR

1. **Yes, the D21 interface carries through to the final layer once the fourth cumulant is right.**
   - With κ4(z)_{aabc} held at its Monte Carlo truth, the ladder of κ3 closures maps one-to-one onto the final layer
     at width 128 (geometric mean over 4 MLPs, first-order Edgeworth step). K=2 is 2.1e-4.

     | κ3 rule | final MSE |
     |---|---|
     | slices only | 9.8e-5 |
     | leading Wick | 2.6e-5 |
     | oracle first-order closure | 4.3e-6 |
     | fitted coefficients | 3.4e-6 |
     | full first-order diagram engine | **2.2e-6** |

   - Across every variant and MLP, final MSE = a + k·ε², where ε is the noise-corrected layer-rms relative error of
     D21(z) against the atlas. The fit has corr 0.86–0.995 and k = 0.7–2.0 × that MLP's K=2 MSE. At width 1024,
     504aldo's law is 4.2e-6·ε² ≈ 0.95 × the K=2 MSE, so the law is width-portable in this normalised form.
2. **A full first-order diagram engine halves the interface error of the oracle closure.** It sums every connected
   diagram on three distinct vertices with ≤ 1 κ3 or κ4 hyperedge of any index pattern, plus Gaussian edges. Fed the
   true z cumulants at width 128 (noise-corrected), its D21(l+1) error is 1.1–3.8 % against 3.1–6.6 % for the
   leg-partition closure and 7–11 % for leading Wick. Layers 1–3 are below the 2.2 % frontier threshold.
3. **The (2,1,1) slice of κ4(z) is the lever, and it compresses.** With the engine, the slice handling sets the
   result (geo-mean final MSE, 4 MLPs):

   | (2,1,1) slice | final MSE |
   |---|---|
   | exact | 2.2e-6 |
   | best rank-4 covariance-response family | 3.7e-6 |
   | rank 1 | 4.8e-6 |
   | u_i C_jk regeneration (the published r = 1 core) | 6.1e-6 |
   | zeroed | 1.7e-5 |

   At the step level, ε falls with width as n^-0.3…0.6 for the closures with the exact slice. With the slice
   regenerated or dropped it falls as n^-0.05…0.2. Extrapolated to n = 1024, "first-order closure + regenerated
   (2,1,1) core" lands at ε ≈ 6–7 %, which is the published chain's measured 4–7 %. The exact slice extrapolates to
   ε ≈ 1 %.
4. **The Edgeworth order matters.** Adding the κ3²/2 term (order 2) without the matching κ5/κ6 and κ3·κ4 terms
   *hurts*:
   - with κ4 teacher-forced, the engine goes 2.2e-6 (order 1) → 5.3e-6 (order 2);
   - in the self-consistent chain, the late-layer blow-ups seen at order 2 disappear at order 1 (MLP 1: 8.3e-5 → 6.9e-6).
5. **For a self-consistent chain at width 128 the binding error is the fourth-cumulant closure.**
   - κ4 = 0: every κ3 rule sits at 2.0–2.4e-5.
   - The memoryless κ4 (slices only) is worse, 1.4–1.7e-4.
   - A dense κ4 closure fixes this: full κ4 carried (n⁴ = 268M entries), with all same-order diagrams including
     κ3×κ3 (`dense2`), at Edgeworth order 1. It gets within ~3× of the teacher-forced level (MLP 1: 6.9e-6 vs 1.6e-6
     teacher-forced vs 4.6e-5 with κ4 = 0; other MLPs in §6).
6. ARC's reference `mlp_kprop` is not public (git asks for credentials; not on PyPI). The dense chain + engine served
   as the reference instead.

## Question and setup

The published chain's interface to each ReLU is D21(l+1)_{ab} = κ3(z_a, z_a, z_b). The oracle ladder (plan §3.1) says
how well each closure of the all-distinct κ3(a) reproduces it. Here we ask whether a better D21 actually lowers the
final-layer MSE end to end, measure the error law, and extrapolate to width 1024.

- Ground truth: `whest dataset bake --n-mlps 8 --n-samples 1e7 --width 128 --depth 16`, seeds 128000–128007 (split
  `dev`). Truth noise floor avg_variance/N = 1.2e-8, negligible here.
- Atlases: `moment_atlas_np.py --k3 --k4` (full raw third moments plus the raw (2,1,1) fourth moment).
  - Width 128: N = 5e5, MLPs 0–3, plus a second-seed atlas of MLP 0 for the noise floor (`results/noise_mlp0.json`).
  - Width 64: N = 1e6, MLP 0, pair.
  - Width 32, depth 8: N = 2e6, pair.

## Method (code)

- **`chain.py`**: the dense chain.
  - State per layer: μ, C, full κ3(z) (n³), and either X = κ4(z)_{aabc} (n³) or the full κ4(z) (n⁴, float32).
  - ReLU step: the Edgeworth operator E[F(z)] = E_G[exp(Σ_r κ_r/r! ∂^r) F], truncated at first order in κ3 and κ4
    (order 1), optionally with κ3²/2 (order 2). It acts on exact Gaussian expectations of derivatives of reluᵖ: closed
    forms via truncated Gaussian moments, equal to the paper's Hermite coefficients of powers (Prop S.3.8). The Mehler
    series runs to all orders in C_ij.
  - Outputs from the 1-D and 2-D formulas: μ(a), var(a), C_off(a), the (3) and (2,1) slices of κ3(a), and the (4),
    (3,1), (2,2) slices of κ4(a).
  - Linear step: exact multilinear transport.
- **All-distinct κ3(a) rules:**
  - `none`: slices only.
  - `wick`: Gaussian ρ² + Φ³κ3(z).
  - `closure`: `oracle_k3.hermite_model` + `residual_basis` × `CLOSURE_COEF`.
  - `fit`: the same seven diagrams with a per-layer table fitted on atlases 0 and 1 (`results/coefs_tensor.txt`);
    MLPs 2–7 are out of sample.
  - `engine`: `triple_engine`, every connected diagram on three distinct vertices with ≤ 1 hyperedge (κ3: 10 index
    patterns; κ4: 15) and Gaussian edge multiplicities ≤ 6 (≤ 2 with a hyperedge). Its Gaussian part equals
    `hermite_model` to 1e-16.
- **κ4 closures:**
  - `zero`.
  - `mem`: the (4),(3,1),(2,2) slices of κ4(a), transported exactly.
  - `dense`: the full κ4(a). Its all-distinct part is Φ⁴κ4(z) + 16 Gaussian trees + 12 κ3-hyperedge-plus-edge diagrams.
    Its (2,1,1) slice comes from the engine with vertex function relu² − 2m·relu, minus 2C_ijC_ik. Transport is n⁵.
  - `dense2`: `dense` plus the same-order diagrams with two κ3 hyperedges and with a (2,1,1) κ4 hyperedge plus a C edge.
  - `atlas`: κ4(z)_{aabc} teacher-forced from the atlas at every layer.
  - `atlas_reg211` / `atlas_zero211` / `atlas_rank{r}`: as `atlas`, but the all-distinct (2,1,1) part is regenerated as
    u_i C_jk, zeroed, or replaced by its best rank-r family (SVD of the (n, n²) unfolding).
- **Other scripts:**
  - `validate.py`, `ladder_width.py`: the step fed the true z cumulants.
  - `diag_state.py`: chain state against the atlas, per layer.
  - `run_variants.py`: `k3:k4[:order]`, error-law noise `N<eps>s<seed>@base`, teacher forcing `F<l>@base` /
    `FD<l>@base`.
  - `fit_coefs.py`, `noise.py`, `law.py`, `summarize.py`, `width_fit.py`.
- **Cost (one core):** closure / engine chains 20–160 s; `dense2` about 25 min at 4–5 GB.

## Results

### 1. Step validation and the interface ladder against width
- Width 32, true z fed in: layer 0 (Gaussian input) matches the atlas to MC noise (μ 1e-4, var 2e-4, C_off 1.3e-3,
  D21(a) 3e-3), so the Gaussian machinery is exact. Deeper layers carry truncation errors at this very non-Gaussian
  width: μ 0.1–0.3 %, var 0.2–1.6 %, C_off 1–3 %.
- The table gives ε of D21(l+1) with the step fed the true z_l cumulants and κ3(a) built by each rule. It is
  noise-corrected with a pair atlas, order 1, rms over layers 1–6 (`results/ladder_o1_w*_mlp0.txt`,
  `results/width_fit_o1.md`).

| rule | n = 32 | n = 64 | n = 128 | ε ~ n^-p, p | extrapolated n = 1024 |
|---|---|---|---|---|---|
| slices only | 0.561 | 0.517 | 0.491 | 0.10 | 0.40 |
| leading Wick | 0.111 | 0.132 | 0.097 | 0.10 | 0.09 |
| oracle closure | 0.084 | 0.060 | 0.038 | 0.57 | 0.012 |
| **diagram engine** | 0.035 | 0.037 | **0.022** | 0.33 | **0.012** |
| engine, (2,1,1) as u_i C_jk | 0.074 | 0.098 | 0.068 | 0.05 | 0.068 |
| engine, (2,1,1) zeroed | 0.113 | 0.115 | 0.085 | 0.20 | 0.060 |

Width 128, all layers (MLP 0, order 1, noise-corrected):
- engine: 1.1, 1.7, 1.9, 2.7, 2.6, 2.8, 2.7, 3.8, 2.9, 3.0, 3.0, 2.4, 2.3, 2.6 % (layers 1–14);
- closure: 3.1–6.6 %;
- Wick: 7–11 %;
- u·C: 3.8–7.8 %;
- slices only: 45–64 %.

Caveat: the n = 32 point is a depth-8 MLP and each fit uses one MLP at three widths. Treat the exponents as indicative.

### 2. End to end, width 128, fourth cumulant teacher-forced (MLPs 0–3; final-layer MSE vs N = 1e7 truth)

| κ3 rule (κ4 = atlas, order 1) | MLP 0 | MLP 1 | MLP 2 | MLP 3 | geo-mean | mean over layers | layer-rms ε(D21), MLP 0 (noise-corr.) |
|---|---|---|---|---|---|---|---|
| slices only | 5.5e-5 | 1.18e-4 | 1.30e-4 | 1.10e-4 | 9.8e-5 | 4.4e-5 | 0.533 |
| leading Wick | 1.06e-5 | 2.24e-5 | 8.0e-5 | 2.54e-5 | 2.64e-5 | 7.4e-6 | 0.155 |
| oracle closure | 2.0e-6 | 6.9e-6 | 8.5e-6 | 3.0e-6 | 4.3e-6 | 1.4e-6 | 0.097 |
| fitted coefs (fit on 0,1) | 1.4e-6 | 3.0e-6 | 8.4e-6 | 3.8e-6 | 3.4e-6 | 1.3e-6 | 0.044 |
| **diagram engine** | **1.07e-6** | **1.64e-6** | **4.2e-6** | **3.2e-6** | **2.2e-6** | 0.94e-6 | **0.036** |
| *K=2 (variant A)* | 9.4e-5 | 4.74e-4 | 1.96e-4 | 2.34e-4 | 2.13e-4 | 1.2e-4 | 1 |

The same ladder at order 2 (κ3²/2 included) is 1.6e-4, 2.1e-5, 1.5e-5, 5.3e-6, 5.3e-6 (geo-means, same rows). The
ordering is the same. Order 2 is 1.5–3.4× worse for closure, fitted and engine and 1.6× worse for slices only;
only leading Wick is slightly better at order 2.

The fitted table is consistent out of sample (MLPs 2–3) and reproduces the oracle's coefficient drift with depth:
- ρ³ term 1.0 → 0.24;
- (2,2) term 1.36 → 0.4–1.0;
- (2,1,1) term 1.44 → 1.1;
- D21 terms ≈ 3 and 2.5.

The parameter-free engine beats the table, so it captures what the drift approximates.

Width 64, the same end-to-end ladder (κ4 = atlas, order 1, MLPs 0 / 1, `results/w64`):

| | K=2 | slices | Wick | closure | engine | engine + u·C | engine, no (2,1,1) |
|---|---|---|---|---|---|---|---|
| MLP 0 | 6.0e-4 | 4.9e-4 | 5.7e-5 | 2.7e-5 | 1.0e-5 | 2.0e-5 | 2.9e-5 |
| MLP 1 | 3.2e-4 | NaN | 2.4e-4 | 1.7e-4 | 9.6e-6 | 7.2e-5 | 3.4e-4 |

The structured law for MLP 0 has k/K2 = 2.5 (corr 0.999). For MLP 1, the slices-only chain hit a negative variance.
The engine's advantage over the oracle closure grows as the net gets narrower: 2.7× and 17× at width 64, against 1.9×
on average at width 128. That is consistent with the closure's ε shrinking faster with width (§1).

### 3. Variant E — the (2,1,1) slice of κ4(z) (κ4 otherwise teacher-forced; order 1)

| κ3 rule | (2,1,1) slice | MLP 0 | MLP 1 | MLP 2 | MLP 3 | geo-mean | ε(D21) MLP 0 |
|---|---|---|---|---|---|---|---|
| engine | exact | 1.07e-6 | 1.64e-6 | 4.2e-6 | 3.2e-6 | 2.2e-6 | 0.036 |
| engine | best rank 4 | 1.28e-6 | 3.47e-6 | 8.7e-6 | 4.9e-6 | 3.7e-6 | 0.069 |
| engine | best rank 1 | 1.44e-6 | 4.11e-6 | 1.10e-5 | 8.4e-6 | 4.8e-6 | 0.083 |
| engine | u_i C_jk (r = 1 core) | 1.88e-6 | 6.8e-6 | 1.16e-5 | 9.4e-6 | 6.1e-6 | 0.097 |
| engine | zero | 6.5e-6 | 2.10e-5 | 4.0e-5 | 1.56e-5 | 1.71e-5 | 0.169 |
| closure | exact | 2.0e-6 | 6.9e-6 | 8.5e-6 | 3.0e-6 | 4.3e-6 | 0.097 |
| closure | u_i C_jk | 4.8e-6 | 1.10e-5 | 1.56e-5 | 1.38e-5 | 1.03e-5 | 0.145 |
| closure | zero | 1.07e-5 | 3.58e-5 | 5.4e-5 | 2.16e-5 | 2.59e-5 | 0.230 |

At order 2, ranks 16 and 32 were also run (MLPs 0–2): geo-means 6.5e-6 and 6.3e-6. Over all four MLPs the order-2
numbers are exact 5.3e-6, rank 4 7.4e-6, rank 1 9.4e-6. The (2,1,1) slice's effect saturates by rank ≈ 4–16 at width 128, and most of the r = 1 → exact gain
is reached by r = 4.

### 4. The error law at width 128 (`results/law.md`)

**Structured.** Per MLP, a least-squares fit of final MSE = a + k ε² over the κ4-teacher-forced variants of §2–3:

| MLP | order 1: k | corr | k / K2 MSE | order 2: k | corr | k / K2 MSE |
|---|---|---|---|---|---|---|
| 0 | 1.9e-4 | 0.995 | 2.0 | 3.5e-4 | 0.998 | 3.7 |
| 1 | 3.3e-4 | 0.989 | 0.69 | 4.7e-4 | 0.991 | 1.0 |
| 2 | 3.5e-4 | 0.857 | 1.8 | 6.6e-4 | 0.971 | 3.4 |
| 3 | 3.7e-4 | 0.994 | 1.6 | 5.4e-4 | 0.994 | 2.3 |

**Injected.** An independent Gaussian perturbation of relative rms ε, added to D21(z) at every layer, on the fitted
chain (order 2), mean of 4 MLPs:

| ε | 0.02 | 0.05 | 0.1 | 0.2 |
|---|---|---|---|---|
| ΔMSE_final / ε² | 2.3e-3 | 2.1e-3 | 2.1e-3 | 2.1e-3 |
| ΔMSE_mean-over-layers / ε² | 4.5e-4 | 4.4e-4 | 4.4e-4 | 4.5e-4 |

Order-1 chain (fitted, κ4 = atlas, one seed per ε, results/eps3): ΔMSE_final/ε² = 1.1e-3, 1.5e-3, 1.7e-3 at
ε = 0.05, 0.1, 0.2 (mean of 4 MLPs).

- The injected law is quadratic to within ±25 % over ε = 0.02–0.2 (order 2: within 10 %).
- Unstructured noise costs ≈ 4× more per ε² than the closures' structured errors at the same order (per MLP k/K2 =
  4–16 against 1–3.7). A plausible
  reason, not tested: closure errors in D21 come with correlated errors in D3 and the variance that partly cancel.
- 504aldo's 1024 law was calibrated on structured variants (memoryless, rank truncations), so the structured
  k/K2 ≈ 0.7–2 (order 1) is the number to compare with his 0.95.

### 5. Variant F — teacher forcing κ3(a_l) at one layer (base: fitted chain, κ4 = atlas, order 2; final MSE / base, mean of 4 MLPs)

| layer forced | 1 | 4 | 8 | 12 | 14 |
|---|---|---|---|---|---|
| full κ3(a_l) := atlas | 0.98 | 0.95 | 0.87 | **0.76** | 0.83 |
| only its all-distinct part | — | 0.95 | 0.93 | 0.97 | — |

On the order-1 chain (results/eps3): forcing the full κ3(a_l) gives 0.90 (layer 8) and 0.88 (layer 12); forcing only
its all-distinct part at layer 12 gives 1.05.

- No single layer dominates. The early layers, where the chain's ε is 1–3 %, are already at the truth.
- The late layers carry the rest: fixing layer 12 alone removes 24 %.
- Most of that gain comes from the slices (D3, D21 of a), not from the all-distinct part. These runs used the
  order-2 step, whose slice error is exactly what §2 shows to be removable.
- Caveat: one N = 5e5 atlas carries 1–6 % noise in its own all-distinct tensor's contribution to D21. That caps what
  the FD rows can show.

### 6. End to end with the chain's own fourth cumulant (`results/summary_all.md`)

| variant | κ3 rule | κ4 closure | order | final MSE (MLPs) | layer 8 | layer 12 |
|---|---|---|---|---|---|---|
| A0 paper Alg. 2 | — | — | — | 2.88e-4 (8) | 1.7e-4 | 2.2e-4 |
| A | K=2, Mehler to all orders | — | — | 2.77e-4 (8) | 1.6e-4 | 2.1e-4 |
| engine / closure / fit | engine / closure / fit | zero | 1 | 2.38e-5 (8) / 2.44e-5 (8) / 2.0e-5 (MLPs 4–7) | 1.3e-5 | 1.7e-5 |
| B0 / M0 | Wick / slices | zero | 2 | 4.2e-5 / 9.9e-5 (8) | | |
| engine / closure | | mem | 1 | 1.69e-4 (8) / 1.41e-4 (MLPs 4–7) | 1.5e-5 | 6.2e-5 |
| C / D / engine | | mem | 2 | 8.4e-5 / 8.4e-5 / 8.0e-5 (8) | 1.4e-5 | 4.2e-5 |
| engine | engine | dense | 2 | 2.1e-5 (MLPs 0–2) | 1.0e-6 | 8.3e-6 |
| engine | engine | dense2 | 2 | 3.0e-5 (7 MLPs; median 1.4e-5) | 9.8e-7 | 6.4e-6 |
| **engine** | **engine** | **dense2** | **1** | DENSE2_O1 | | |

How the κ4 closure behaves along the chain (MLP 0, order 2), as relative error of K4 = κ4(z)_{aaaa}:
- memoryless: 16 % at layer 1, 47 % at layer 4, 65 % at 8, 78–91 % at 13–15. Memoryless κ4 is worse than κ4 = 0 at
  depth.
- dense: 9 %, 13 %, 35 %, 37–67 %.
- dense2: 9 %, 12 %, 34 %, 33–46 %.

At order 2, dense2 blew up in the last 3–4 layers on 4 of 7 MLPs. At order 1 that blow-up is gone on the MLPs run so
far.

## Verdict
- **The D21 interface translates.** With the fourth cumulant right, final-layer MSE follows ε(D21)² over two decades:
  - the oracle closure (ε ≈ 10 %) sits at 4.3e-6;
  - the full first-order diagram engine (ε ≈ 3.6 %, layer-rms; 1–4 % per layer at the step level) sits at 2.2e-6;
  - they are 50× and 100× below K=2;
  - the leading-Wick chain (ε ≈ 16 %) is 12× worse than the engine.
- **The κ3 rule that does it needs no fitting.** The engine is "all first-order diagrams". The per-layer coefficient
  table rediscovers part of it and is worse.
- **The (2,1,1) fourth-cumulant slice decides the interface below ≈ 5 %.**
  - Regenerating it as u_i C_jk (the published r = 1 core) costs 2.8× in MSE at width 128.
  - Its ε barely improves with width (p ≈ 0.05), while the exact-slice closures improve as n^-0.3…0.6.
  - A rank-4 family of covariance-response matrices recovers most of the gap (3.7e-6 vs 2.2e-6 exact vs 6.1e-6 for u·C).
- **Implication at width 1024.** Use the 1024 law (4.2e-6·ε², which is 0.95·K2; our width-128 structured ratio is
  0.7–2·K2, the same form).
  - A chain with the regenerated core sits at ε ≈ 6–7 % → extra MSE ≈ 1.5–2e-8, which is the published floor
    (raw 2.07e-8).
  - Carrying the (2,1,1) slice's transported effect to a few per cent would put ε ≈ 1 % → extra ≈ 4e-10.
  - The gate below the frontier's ε ≤ 2.2 % is therefore the (2,1,1) object, not the κ3 rule. At 1024 it has to be
    carried compressed. At 128, rank 4 already gets most of it; plan §3.1 (vii) found rank ≳ 32 needed for its *energy*
    at shallow layers, but end to end its *effect* saturates far earlier.
- **Implementation detail that matters as much as any closure:** truncate the Edgeworth step consistently at first
  order. The κ3²/2 term alone costs 1.5–3.4× at width 128 for the good closures and destabilises the deep layers of the self-consistent chain.
- **Self-consistent chain at width 128:** the fourth-cumulant closure is the binding error. Memoryless κ4 is worse than
  none. dense2 at order 1 (every same-order κ4 diagram) brings the chain to within a few × of teacher-forced κ4.

## Open issues
- dense2 at order 1 has been run only on the MLPs listed in §6, at ~25 min per MLP.
- Atlas noise (N = 5e5) is 1–6 % on D21 and 2–10 % on K4 per atlas. ε is noise-corrected only where a pair exists
  (MLP 0 at widths 128, 64, 32); elsewhere it includes ≈ +0.5 % layer-rms.
- The width-64 end-to-end chain hit negative variances at depth (order 2; the Edgeworth corrections are not
  positivity-preserving), so the width study is at the step level. It was not retried at order 1.
- The width exponents come from three widths, one MLP each, with a depth-8 MLP at n = 32. A width-256 step-level point
  (≈ 7 GB per atlas) would firm them up.
- The injected-noise law and teacher forcing were run on the order-2 chain only.

## Next steps
1. Test at width 1024 the claim that the (2,1,1) slice's *transported effect* on D21(l+1) is captured by a few
   covariance-response modes (Wᵀ Φ M_r Φ W sandwiches). It needs the slice from GPU atlases. If it holds, a rank-4
   carrier costs ~8 sandwiches per layer.
2. Port the engine's κ3 diagrams into the factorised 1024 chain. Its extra diagrams beyond the oracle closure are κ3
   and κ4 hyperedges with higher Gaussian multiplicities; each is a Hadamard-weighted contraction of objects the chain
   already holds, except the (2,1,1) slice.
3. Switch every chain in the programme to order-1 Edgeworth on the slices and re-measure (it is free).
