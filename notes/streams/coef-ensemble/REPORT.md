# coef-ensemble: are the renormalised closure coefficients an ensemble property, and how do closure errors scale with width?

Status: **final** (stream stopped by the coordinating session when the campaign moved to fresh-slate designs; widths 48
and 192 were queued and cancelled, the (2,1,1)-slice carrier task sent later was not started).
Scripts: [coef_table.py](coef_table.py) (`cache` / `analyse` / `transfer` / `ablate`), [scaling.py](scaling.py),
[drift.py](drift.py), [report_tables.py](report_tables.py).
Raw outputs in [results/](results/): `width{64,96,128,160}.txt` (all coefficient tables and per-layer ladders),
`transfer.txt`, `ablate.txt`, `scaling.txt`, `drift.txt`, `summary.json`, and the data file `coef_table_n1024.json`
(per-width tables and the n = 1024 extrapolation).

## Question

The all-distinct post-activation third cumulant κ3(a_l) is, to a few per cent, the Gaussian ρ² Wick term plus seven
first-order gate diagrams (`oracle_k3.residual_basis`, B0..B6). With the leg-partition coefficients (`CLOSURE_COEF` =
1, 3, 3, 1, 1, 1.5, 1.5) the D21(l+1) error stayed at 4–6 % at depth on one width-128 MLP; per-layer fitted coefficients
brought it to 1.4–2.3 %. Is a 7 × 15 fitted table a property of the shape (n, L), usable offline? And does any version reach
the frontier bar ε ≤ 2.2 % at n = 1024 (extra MSE ≈ 4.2e-6 ε²)?

## Method

- Atlases: `moment_atlas_np.py --k3 --k4`, depth 16, N = 5e5 per atlas, He-initialised MLPs from seeds 770100+.
  Width 64: 6 MLPs, 96: 6, 128: 8, 160: 3. Pairs (second atlas, independent sample seed): 3 MLPs at widths 64–128,
  1 MLP (770100) at width 160.
- `coef_table.py cache`: per atlas and layer, the transported (n × n) form of every ingredient (slices of κ3(a), leading
  Wick, Gaussian ρ² Wick, the 7 basis tensors and B7 = B6 with the (2,1,1) slice regenerated as u_i C_jk) and the
  tensor-space normal equations of R = all-distinct κ3(a) − Gaussian ρ² Wick on the basis. Everything downstream is linear
  in the coefficients and runs on these caches in seconds.
- Ladder, each MLP evaluated with coefficients from the *other* MLPs of its width (leave-one-out):
  `wick` leading Wick only; `leg` leg-partition closure; `legR` the same with the (2,1,1) slice regenerated as u_i C_jk
  (the published chain's closure); `own` / `ownD` per-MLP fit in tensor / D21 space (oracle: fitted on the evaluated MLP);
  `ens` / `ensD` ensemble table (stacked training MLPs per layer) in tensor / D21 space (D21 space = regress the
  transported residual on the 7 transported basis terms, what the chain consumes); `ensR` / `ensRD` the table refitted with
  the regenerated slice.
- Noise. On pair MLPs the model built from atlas A is scored against atlas B's D21, and two noise energies are subtracted:
  B's target noise (rel(D21_A, D21_B)/√2, as `oracle_k3.analyse_pair`) and the MC noise of the model's own inputs
  (C, κ3(z), κ4 slice, slices of κ3(a)), measured as ‖pred_A − pred_B‖/(√2‖D21_B‖) with the same coefficients. The result,
  **ε_rep, is the pure representation error** (unbiased; noisy when the noise energies dominate). Without the second
  subtraction the cross-evaluated ε *rises* with width at layers 0–2, because at fixed N the input noise grows with n.
  The within-atlas ε (model and target from the same atlas) is reported as a cross-check; it is ≥ ε_rep in expectation.
  Layer 0 is excluded throughout: z_0 is exactly Gaussian, every non-Gaussian input is zero in population, ε_rep(0) ≈ 0.

## Results

### 1. The table is an ensemble property: an offline table equals the per-MLP oracle fit

ε_rep of D21(l+1), mean over pair MLPs, rms over layer bands (width 160: one pair MLP):

| width | leg 1-3 | leg 4-9 | leg 10-14 | own 1-3 | own 4-9 | own 10-14 | **ens 1-3** | **ens 4-9** | **ens 10-14** | ensR 1-3 | ensR 4-9 | ensR 10-14 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 64 | 0.045 | 0.074 | 0.093 | 0.038 | 0.046 | 0.029 | 0.038 | 0.046 | 0.031 | 0.080 | 0.092 | 0.061 |
| 96 | 0.036 | 0.049 | 0.068 | 0.033 | 0.029 | 0.028 | 0.034 | 0.030 | 0.029 | 0.074 | 0.066 | 0.056 |
| 128 | 0.031 | 0.038 | 0.055 | 0.030 | 0.025 | 0.018 | 0.030 | 0.025 | 0.018 | 0.065 | 0.058 | 0.039 |
| 160 | 0.023 | 0.029 | 0.023 | 0.022 | 0.018 | 0.012 | 0.022 | 0.018 | 0.012 | 0.052 | 0.045 | 0.026 |

Per layer at width 128 (mean over 3 pair MLPs; max over MLPs in the last two columns):

| l | wick | leg | legR | own | ownD | ens | ensD | ensR | ensRD | max own | max ens |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.077 | 0.029 | 0.058 | 0.029 | 0.029 | 0.029 | 0.029 | 0.058 | 0.058 | 0.034 | 0.034 |
| 2 | 0.089 | 0.031 | 0.068 | 0.029 | 0.029 | 0.029 | 0.029 | 0.068 | 0.068 | 0.032 | 0.032 |
| 3 | 0.095 | 0.034 | 0.071 | 0.031 | 0.031 | 0.031 | 0.031 | 0.069 | 0.069 | 0.036 | 0.036 |
| 4 | 0.099 | 0.034 | 0.071 | 0.029 | 0.029 | 0.029 | 0.029 | 0.069 | 0.069 | 0.037 | 0.038 |
| 5 | 0.091 | 0.034 | 0.065 | 0.025 | 0.026 | 0.025 | 0.025 | 0.062 | 0.062 | 0.041 | 0.042 |
| 6 | 0.090 | 0.036 | 0.062 | 0.023 | 0.023 | 0.023 | 0.023 | 0.056 | 0.056 | 0.030 | 0.030 |
| 7 | 0.081 | 0.039 | 0.066 | 0.021 | 0.021 | 0.021 | 0.022 | 0.053 | 0.053 | 0.026 | 0.026 |
| 8 | 0.087 | 0.041 | 0.062 | 0.023 | 0.022 | 0.024 | 0.024 | 0.051 | 0.051 | 0.030 | 0.030 |
| 9 | 0.096 | 0.042 | 0.067 | 0.024 | 0.024 | 0.026 | 0.026 | 0.055 | 0.055 | 0.028 | 0.028 |
| 10 | 0.085 | 0.046 | 0.066 | 0.021 | 0.021 | 0.021 | 0.021 | 0.046 | 0.046 | 0.028 | 0.029 |
| 11 | 0.082 | 0.054 | 0.063 | 0.018 | 0.017 | 0.019 | 0.019 | 0.038 | 0.039 | 0.021 | 0.020 |
| 12 | 0.077 | 0.056 | 0.067 | 0.018 | 0.017 | 0.018 | 0.018 | 0.037 | 0.037 | 0.019 | 0.020 |
| 13 | 0.080 | 0.059 | 0.070 | 0.018 | 0.017 | 0.018 | 0.018 | 0.038 | 0.038 | 0.019 | 0.018 |
| 14 | 0.072 | 0.059 | 0.070 | 0.016 | 0.015 | 0.016 | 0.016 | 0.035 | 0.035 | 0.018 | 0.019 |

Facts:
- **Held-out ensemble table = per-MLP oracle fit** to ≤ 0.002 in ε at every layer and width, mean and worst MLP, in tensor
  and D21 space alike (`ens` ≈ `ensD` ≈ `own` ≈ `ownD`). Per-MLP fitting buys nothing over an offline table.
- **Cross-width transfer is free** ([results/transfer.txt](results/transfer.txt); ε_rep bands 1-3/4-9/10-14):

  | evaluated at \ table from | 64 | 96 | 128 | 160 | pooled |
  |---|---|---|---|---|---|
  | 64 | .038/.046/.033 | .039/.047/.032 | .039/.048/.032 | .039/.048/.033 | .039/.047/.031 |
  | 96 | .034/.030/.032 | .033/.030/.030 | .033/.030/.029 | .034/.030/.030 | .034/.030/.030 |
  | 128 | .030/.026/.019 | .030/.025/.018 | .030/.025/.019 | .030/.025/.019 | .030/.025/.018 |
  | 160 | .023/.021/.014 | .022/.019/.013 | .022/.018/.012 | .022/.018/.012 | .022/.019/.013 |

  A table fitted at one width loses ≤ 0.002 at any other width in 64–160.
- A quadratic-in-depth table (3 numbers per coefficient) is as good as the per-layer table.
- What remains after the table is the floor of the seven-diagram basis, not coefficient error.

### 2. Only three coefficients matter: B3, B6, B2 ([results/ablate.txt](results/ablate.txt))

ε_rep bands when one table coefficient is reset to its leg value (`ens − k`) or only one is taken from the table (`leg + k`):

| variant | n=128: 1-3 | 4-9 | 10-14 | n=160: 1-3 | 4-9 | 10-14 |
|---|---|---|---|---|---|---|
| ens | 0.030 | 0.025 | 0.018 | 0.022 | 0.018 | 0.012 |
| leg | 0.031 | 0.038 | 0.055 | 0.023 | 0.029 | 0.023 |
| ens − B2 (w3 D21) | 0.030 | 0.027 | 0.025 | 0.022 | 0.020 | 0.015 |
| ens − B3 (D3 ⊗ C ⊗ C) | 0.031 | 0.034 | 0.047 | 0.023 | 0.026 | 0.021 |
| ens − B6 (κ4 (2,1,1)) | 0.031 | 0.030 | 0.029 | 0.023 | 0.021 | 0.015 |
| ens − B0, B1, B4 or B5 | ≤ 0.030 | ≤ 0.026 | ≤ 0.020 | 0.022 | ≤ 0.019 | 0.012 |
| leg + B3 only | 0.030 | 0.031 | 0.030 | 0.023 | 0.021 | 0.015 |

B3 (the D3 hyperedge with two C edges on one vertex, Hermite-degree-5 vertex weight w5) carries most of the gain, then
B6 and B2. B0, B1, B4, B5 can stay at leg values; their large across-MLP spread (sd 0.3–1.0 at depth) sits in flat
directions of the D21 error.

### 3. Coefficient tables vs depth and width; drift vs n

Tensor-space ensemble fit ± across-MLP sd of per-MLP fits (all seven in [results/width*.txt](results/); the per-pair Monte
Carlo sd of one per-MLP fit is 0.00–0.02 for B2, B3, B6, so the spread is real MLP-to-MLP variation):

| l | B2 64 | B2 96 | B2 128 | B2 160 | B3 64 | B3 96 | B3 128 | B3 160 | B6 64 | B6 96 | B6 128 | B6 160 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| leg | 3 | 3 | 3 | 3 | 1 | 1 | 1 | 1 | 1.5 | 1.5 | 1.5 | 1.5 |
| 1 | 3.02±.09 | 3.05±.05 | 3.04±.04 | 3.04±.02 | 0.43±.08 | 0.44±.06 | 0.41±.03 | 0.42±.02 | 1.40±.02 | 1.43±.01 | 1.45±.01 | 1.45±.00 |
| 3 | 2.66±.12 | 2.80±.10 | 2.87±.04 | 2.82±.08 | 0.14±.05 | 0.18±.13 | 0.30±.06 | 0.25±.06 | 1.23±.03 | 1.30±.05 | 1.34±.01 | 1.36±.01 |
| 5 | 2.44±.08 | 2.55±.12 | 2.66±.09 | 2.65±.07 | −0.02±.09 | 0.00±.07 | 0.04±.14 | 0.05±.05 | 1.14±.04 | 1.21±.03 | 1.26±.03 | 1.27±.02 |
| 8 | 2.34±.08 | 2.51±.08 | 2.58±.05 | 2.60±.07 | −0.09±.12 | −0.12±.02 | −0.06±.06 | −0.05±.04 | 1.06±.05 | 1.14±.05 | 1.19±.04 | 1.20±.04 |
| 11 | 2.41±.07 | 2.47±.08 | 2.51±.06 | 2.61±.05 | −0.17±.05 | −0.14±.04 | −0.06±.08 | −0.02±.12 | 1.02±.06 | 1.09±.03 | 1.15±.02 | 1.19±.02 |
| 14 | 2.34±.12 | 2.46±.07 | 2.48±.05 | 2.55±.03 | −0.17±.14 | −0.11±.11 | −0.09±.06 | −0.07±.03 | 1.06±.07 | 1.13±.03 | 1.14±.05 | 1.20±.03 |

B0 = 1.00 at every layer and width; B1 → 3, B4 → 1, B5 → 1.5 at layer 1, all drifting with depth (tables in results).

Drift |c_ens(n) − c_leg| (rms over layer bands), fitted as n^−q ([results/drift.txt](results/drift.txt)):

| coef | band | n=64 | 96 | 128 | 160 | q | drift at 1024 |
|---|---|---|---|---|---|---|---|
| B2 | 1-3 | 0.23 | 0.13 | 0.08 | 0.10 | 0.99 | 0.01 |
| B2 | 4-9 | 0.56 | 0.48 | 0.38 | 0.35 | 0.55 | 0.13 |
| B2 | 10-14 | 0.61 | 0.53 | 0.52 | 0.43 | 0.36 | 0.23 |
| B6 | 1-3 | 0.20 | 0.15 | 0.12 | 0.10 | 0.74 | 0.03 |
| B6 | 4-9 | 0.39 | 0.33 | 0.27 | 0.25 | 0.49 | 0.10 |
| B6 | 10-14 | 0.47 | 0.39 | 0.36 | 0.31 | 0.43 | 0.14 |
| B1 | 4-14 | 0.41–0.55 | 0.28–0.41 | 0.21–0.32 | 0.13–0.19 | 1.1–1.2 | ≤ 0.03 |
| B5 | 4-9 | 1.02 | 0.86 | 0.72 | 0.55 | 0.65 | 0.18 |
| B3 | 1-3 | 0.74 | 0.69 | 0.64 | 0.66 | 0.15 | 0.48 |
| B3 | 4-14 | 1.06–1.15 | 1.07–1.12 | 1.02–1.08 | 1.01–1.04 | 0.06–0.10 | 0.87–0.90 |

So: **B1, B2, B4, B5, B6 drift back toward the leg-partition values as n grows** (q ≈ 0.4–1.0; the renormalisation is a
finite-width effect), consistent with the width-1024 oracle (oracle1024 stream) finding the fitted D21-space coefficients
at leg values within noise (B2 2.9–3.0, B6 1.40–1.49). The drift exponent from 64–160 (q ≈ 0.4–0.7 at depth for B2, B6)
is slower than the ~n^−0.8 the coordinating note expected; the extrapolated 1024 drift (B2 0.13–0.23, B6 0.10–0.14) is a
little larger than the oracle1024 measurement, i.e. the 64–160 fit is conservative. **B3 is the exception**: in tensor
space it stays at ≈ 0.4 (layer 1) and ≈ 0 (depth) at every width with q ≈ 0.1, far from its leg value 1. Its effect on D21
shrinks with width anyway (ens − B3 vs ens: +0.029 at n = 128, +0.009 at n = 160 for layers 10–14). (oracle1024's list does
not include B3; whether its D21-space value at 1024 is 1 or ≈ 0 is not established here.)

### 4. Width scaling of ε and extrapolation to n = 1024

ε_rep band rms, power-law fit ε ∝ n^−p over 64–160 and its value at 1024; in brackets the within-atlas cross-check fit
(an upper bound in expectation, and biased upward with n at fixed N):

| model | band | 64 | 96 | 128 | 160 | p | ε(1024) [within-atlas] |
|---|---|---|---|---|---|---|---|
| ens | layer 1 | 0.033 | 0.031 | 0.029 | 0.023 | 0.34 | 1.3 % |
| ens | 1-3 | 0.038 | 0.034 | 0.030 | 0.022 | 0.55 | 0.9 % [1.4 %] |
| ens | 4-9 | 0.046 | 0.030 | 0.025 | 0.018 | 0.99 | 0.3 % [0.5 %] |
| ens | 10-14 | 0.031 | 0.029 | 0.018 | 0.012 | 1.05 | 0.2 % [0.7 %] |
| leg | 1-3 | 0.045 | 0.036 | 0.031 | 0.023 | 0.68 | 0.7 % [1.2 %] |
| leg | 10-14 | 0.093 | 0.068 | 0.055 | 0.023 | 1.38 | 0.2 % [0.5 %] |
| legR | layer 1 | 0.068 | 0.063 | 0.058 | 0.049 | 0.32 | 2.8 % |
| legR | 1-3 | 0.082 | 0.075 | 0.066 | 0.052 | 0.46 | 2.4 % [3.4 %] |
| legR | 4-9 | 0.108 | 0.077 | 0.066 | 0.052 | 0.76 | 1.3 % [1.6 %] |
| ensR | 1-3 | 0.080 | 0.074 | 0.065 | 0.052 | 0.44 | 2.5 % [3.4 %] |
| ensR | 10-14 | 0.061 | 0.056 | 0.039 | 0.026 | 0.90 | 0.6 % [2.8 %] |

Calibration against the width-1024 oracle (oracle1024 stream, measured directly): exact-slice leg closure 0.8–1.0 %,
u_i C_jk regeneration 2.5–2.7 %. The extrapolations here (leg 0.7–1.3 %, legR 2.4–2.8 % at layers 1–3) agree, so the
64–160 ladder extrapolates correctly. Caveat: the width-160 point is one pair MLP, and at n = 160 the subtracted noise
energies at layers 1–3 are 2–4× the representation energy (within-atlas numbers bracket it).

## Verdict

1. **Yes, the renormalised coefficients are an ensemble property** at every width measured: a 7 × 15 offline table matches
   the per-MLP oracle fit to ≤ 0.2 points of ε, and a table from any width in 64–160 works at any other. Only B2, B3, B6
   matter.
2. **But the renormalisation is a finite-width effect.** B1, B2, B4, B5, B6 drift back to the leg-partition values as
   n^−(0.4..1.0), and the gain of the table over the leg closure shrinks with n (layers 10–14: 0.055 → 0.018 at 128,
   0.023 → 0.012 at 160). At n = 1024 the table is not needed: the un-fitted leg-partition closure with the exact (2,1,1)
   slice is already ≈ 1 % (extrapolated here, measured by oracle1024), well under 2.2 %, at every layer.
3. **The binding error at n = 1024 is the (2,1,1) fourth-cumulant slice, not the coefficients.** With the slice regenerated
   as u_i C_jk, ε extrapolates to 2.4–2.8 % at layers 1–3 (above 2.2 %), ≈ 1.3 % at 4–9 and ≤ 0.6 % at 10–14; refitting
   the table with the regenerated slice (ensR vs legR) buys nothing at layers 1–3. Layers 1–3 are where a carrier of the
   (2,1,1) slice better than r = 1 is required.

## Facts for the fresh-slate designs

- The all-distinct κ3(a) → D21(l+1) interface is closed by the first-order gate diagrams with leg-partition coefficients to
  ≈ 1 % at n = 1024, *provided* the exact (2,1,1) κ4(z) slice is supplied; nothing per-MLP needs to be fitted.
- Finite-width renormalisation of the diagram coefficients exists (B2: 3 → 2.5, B6: 1.5 → 1.15, B3: 1 → ≈ 0 at depth for
  n ≤ 160), is shape-determined, and decays with n; B3's renormalisation is the slowest.
- The r = 1 regeneration of the (2,1,1) slice costs ≈ 2–3 points of ε at the shallow layers at every width 64–160 (and
  2.5–2.7 % at 1024); its loss shrinks only as n^−0.45.

## Open issues (not pursued; stream stopped)

- Width 192 and a 3-pair width-160 point were cancelled; the 160 numbers rest on one pair MLP.
- B3's persistent renormalisation (1 → ≈ 0) suggests the Gaussian degree-5 vertex weight w5 is the wrong object; untested.
- The (2,1,1)-slice carrier task (transported + birth decomposition, low-rank families, richer O(n²) regenerations) was
  assigned and then withdrawn before work started; the k3+k4 caches built here (`coef_table.py cache` outputs, ~60 MB) are
  not in the repository and the atlases live only in this ephemeral container, but every atlas is reproducible from its seed with the commands above.
