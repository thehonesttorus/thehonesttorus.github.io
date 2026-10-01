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
