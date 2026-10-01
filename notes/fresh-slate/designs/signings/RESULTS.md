# signings stream — results log

Machine: 4 shared cores. Truth: own Monte Carlo bake (`truth.py`, weights from `whestbench.generation.sample_mlp`),
float32 forward, float64 accumulation; raw MSE has the per-neuron truth noise var_j/N subtracted.

## T0 — exactness of one copula step (layer 1 is an exact Gaussian copula)

Per-sector check at width 12 (`t0_sectors.py`, 4e6 MC samples of a_1, exact joint cumulants vs the Mehler
generators of `copula.py`):

| sector of κ_k(z_2j) | generator rms | MC rms | diff rms | note |
|---|---|---|---|---|
| κ3 {2,1} | 0.2923 | 0.2925 | 0.0013 | exact (MC noise) |
| κ3 {1,1,1} paths | 0.0821 | 0.0772 | 0.0105 | MC includes triangles |
| κ4 {2,2} | 0.1639 | 0.1646 | 0.0036 | exact |
| κ4 {3,1} | 0.3821 | 0.3825 | 0.0039 | exact |
| κ4 {2,1,1} paths | 0.2301 | 0.2427 | 0.1848 | MC includes triangles |

Entry-level check of κ(a_a,a_a,a_b,a_c) (`t0_211.py`): the full three-vertex Mehler sum *including* the
triangle multigraphs reproduces MC (−0.0511 vs −0.0515, 0.0976 vs 0.0971, ...), the tree-only sum does not
(−0.0305): at width 12 (|ρ| ≈ 0.3) loops are not negligible; the generators are correct.

T1 (loop share vs width, `t0_exact.py`, κ3/κ4 of z_2 vs 2e7-sample MC, tree sectors only):

| n | κ3 rms err / rms κ3 | κ4 rms err / rms κ4 |
|---|---|---|
| 12 | 0.0099 / 0.94 | 0.184 / 1.38 |
| 24 | 0.0085 / 0.54 | 0.209 / 0.80 |
| 48 | 0.0043 / 0.29 | 0.050 / 0.45 |
| 96 | 0.0017 / 0.15 | 0.013 / 0.20 |

MC noise (2e7 samples) ≈ 1.4e-3 on κ3 and ≈ 4e-3 on κ4. At n = 96 the tree-only κ3 is at the noise floor, and the
κ4 loop residual falls ≈ n^{-2} from 48 to 96. **T1 passes: loop (cycle) diagrams with fresh-weight legs are
negligible at the widths that matter, so the determinantal/signed resummation has nothing to do in the step
sums** (DESIGN §1 Step A).

## Layer-by-layer diagnosis at width 64 (`diag_mc.py`, 4e6-sample MC statistics of z_2..z_8, MLP w64 seed 0)

rms errors, estimate vs MC (value after "/" is the rms of the quantity itself):

| layer | v0 copula: var rel | v0 offcov | v0 κ3 | v1 (+ carried (2,1) slice D): var rel | v1 offcov | v1 κ3 | v1 D off-diag | mean(a) MSE v0 / v1 |
|---|---|---|---|---|---|---|---|---|
| z_2 | 7.6e-4 | 6.8e-4 / 0.22 | 3.7e-3 / 0.25 | 7.6e-4 | 6.8e-4 | 3.7e-3 | 2.1e-3 / 0.11 | 3.2e-7 / 3.2e-7 |
| z_3 | 3.3e-2 | 2.7e-2 | 0.14 / 0.34 | 5.9e-3 | 5.4e-3 | 0.083 | 0.049 / 0.14 | 3.7e-5 / 7.2e-6 |
| z_4 | 7.1e-2 | 4.7e-2 | 0.26 / 0.33 | 2.1e-2 | 1.3e-2 | 0.13 | 0.067 / 0.14 | 1.2e-4 / 2.9e-5 |
| z_6 | 8.5e-2 | 4.5e-2 | 0.22 / 0.28 | 3.6e-2 | 1.6e-2 | 0.10 | 0.057 / 0.12 | 3.5e-4 / 8.3e-5 |
| z_8 | — | — | — | 5.7e-2 | 2.3e-2 | 0.11 / 0.32 | 0.059 / 0.13 | — / 1.3e-4 |

Reading: one step out of an exact copula is essentially exact (z_2; the carried slice D matches MC to its noise
once the coincident-index tree terms are included — they were missing in the first v1 and cost 13 %). From
z_3 on, the state is not a Gaussian copula and the error is the copula defect, as DESIGN §4 predicted.

Sector split of κ3(z_3) (`dbg_L3.py`, MC sector values from 3e6 samples of a_2):

| sector | MC rms | v1 error rms | v2 (+ diagonal-birth CP sources) error |
|---|---|---|---|
| {3} | 0.122 | 0.0006 | 0.0006 |
| {2,1} (K21 of a_2 is 10 % off) | 0.217 | 0.048 | 0.048 |
| {1,1,1} all-distinct | 0.097 | 0.071 | 0.056 |

The all-distinct sector is ≈ 70 % missed by the copula's tree paths: it is content born at earlier ReLUs (one
cumulant hyperedge) carried by linear gate chains. Carrying only the diagonal births (v2) recovers ≈ 20 % of it;
the rest is path-born (star-shaped) content of the earlier layer and second-order hyperedge terms. v2 with all
ages is unstable at width 64 (the carried slice drifts, the hyperedge correction makes Cov(a) indefinite on
2/8 MLPs).

## Stage Q, width 64 (shared bench `w64_d16`, 8 MLPs, N = 1e7, truth noise 6.6e-9)

| estimator | raw final MSE | ± s.e. | all-layer MSE |
|---|---|---|---|
| Gaussian closure (bench baseline) | 5.11e-4 | 8.5e-5 | 5.10e-4 |
| v0 copula | 3.76e-4 | 5.8e-5 | 3.25e-4 |
| v1 copula + carried (2,1) slice + hyperedge correction | 1.75e-4 | 2.9e-5 | 1.11e-4 |
| v1 with single-edge paths only (p = q = 1) | 2.65e-4 | 5.8e-5 | 1.34e-4 |
| v2 (v1 + CP sources, all ages) | unstable (NaN on 2/8) | | |

## Stage Q, widths 64–512, depth 16 (final-layer raw MSE, truth noise subtracted; mean ± s.e. over MLPs)

w64: shared bench `w64_d16` (8 MLPs, N = 1e7). w128/256/512: own bakes of `sample_mlp(width, 16, seed=s)` weights
(`truth.py`; N = 8e6 / 3e6 / 1.5e6; truth noise ≤ 8e-8, subtracted per MLP).

| estimator | 64 (8) | 128 (4) | 256 (4) | 512 (3) | fit slope p (raw ∝ n^-p) | raw at 1024 (fit, 1σ) |
|---|---|---|---|---|---|---|
| Gaussian closure (bench `baseline_gauss`) | 5.1e-4 ± 0.9e-4 | 3.2e-4 ± 0.8e-4 | 6.1e-5 ± 0.9e-5 | 2.1e-5 | 1.62 ± 0.23 | 7.3e-6 (4.8e-6 – 1.1e-5) |
| v0 copula (marginals + latent R) | 3.8e-4 ± 0.6e-4 | 2.2e-4 ± 0.5e-4 | 5.1e-5 ± 0.9e-5 | 1.3e-5 (1) | 1.66 ± 0.21 | 4.9e-6 (3.3e-6 – 7.2e-6) |
| **v1** copula + carried (2,1) slice + hyperedge terms | 1.75e-4 ± 0.3e-4 | 6.8e-5 ± 1.7e-5 | 1.41e-5 ± 0.2e-5 | 5.1e-6 ± 0.7e-6 | **1.76 ± 0.12** | **1.4e-6 (1.1e-6 – 1.8e-6)** |
| v1, single-edge paths (p = q = 1; the costed form) | 2.6e-4 | — | 1.43e-5 | 5.1e-6 | | same as v1 at n ≥ 256 |

Calibration caveat: the brief quotes ≈ 4e-5 for Gaussian closure at 1024, while my extrapolation of the bench
Gaussian baseline gives 7e-6. If extrapolation from ≤ 512 is optimistic by the same factor (≈ 5×), v1 at 1024 is
≈ 7e-6. v1 is consistently 0.22–0.34 × the Gaussian baseline at n ≥ 128, which gives a second projection,
0.25 × (7e-6 … 4e-5) = 2e-6 … 1e-5.

Layer profile (v1, width 512): L2 ≈ 0 (at noise), L4 4.5e-7, L8 2.2e-6, L16 5.1e-6: error is born from layer 3 on and
accumulates.

## Teacher forcing at width 128 (`ablate.py`, MC statistics of every z_l from 3e6 samples injected at layers ≤ 15,
final step 15 → 16 computed by v1; MLP w128 seed 0)

| injected at every layer ≤ 15 | final raw MSE |
|---|---|
| nothing (v1) | 1.23e-4 |
| D (true (2,1) slice) | 1.53e-4 |
| true covariance | 1.06e-4 |
| true marginals (m, var, κ3, κ4) | 3.6e-5 |
| marginals + D | 3.1e-5 |
| marginals + covariance + D | 3.7e-5 |

**A single step out of an exactly specified state (marginals, covariance, (2,1) slice) has error 3.6e-5 at width
128**, about 1/3 of v1's accumulated error. Decomposition of that step (`ablate2.py`, `margtest.py`):
- the marginal transport is not the problem: a 4-cumulant (Fleishman) marginal fed with exact cumulants reproduces
  E a to 5e-4 rms and κ3(a) to 0.5 % at depth 15 (skew 0.48);
- the next-layer variance is off by 0.5–0.9 % (relative), almost all of it from the **off-diagonal Cov(a)**: rms entry error
  9.6e-4 on 0.11. The first-order hyperedge term already removes ≈ 80 % of the copula's bivariate defect
  (5.3e-3 → 9.6e-4). The rest is second-order bivariate structure (pairwise fourth-order slices κ(z_a,z_a,z_b,z_b),
  κ(z_a,z_a,z_a,z_b), and Δ² terms), which v1 neither carries nor generates;
- κ3(z') is off by 15–30 % in the same step, split between the {2,1} sector (K21 of a, 10 % off) and the all-distinct sector.

## Teacher forcing at width 256 (same protocol, 1.5e6-sample statistics; MLP w256 seed 0)

| injected at every layer ≤ 15 | final raw MSE |
|---|---|
| nothing (v1) | 2.05e-5 |
| variance only | 1.68e-5 |
| variance + mean | 6.2e-6 |
| standardized shape (κ3, κ4) only | 9.5e-6 |
| marginals (m, var, κ3, κ4) | 4.3e-6 |
| marginals + D | 2.7e-6 |
| marginals + covariance + D | 2.55e-6 |

**One step out of an exact state: 3.6e-5 (n = 128) → 2.55e-6 (n = 256)**, a factor 14, ≈ n^{-3.8} (one MLP per width,
so the slope is uncertain to roughly ±1). Full v1 falls only as n^{-1.8}. At large width v1's error is therefore
**accumulated drift of the carried per-neuron statistics** (means, variances and shapes all contribute; none dominates),
not the closure of a single step.
