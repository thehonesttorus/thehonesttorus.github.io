# coef-ensemble: are the renormalised closure coefficients an ensemble property, and how do closure errors scale with width?

Status: **draft with widths 64, 96, 128** (width 160, 192 and 48 atlases are building; numbers below will be extended).
Scripts: [coef_table.py](coef_table.py) (cache / analyse / transfer / ablate), [scaling.py](scaling.py), [report_tables.py](report_tables.py).
Raw outputs: [results/](results/) (`width*.txt` full tables, `transfer.txt`, `ablate.txt`, `scaling.txt`, `summary.json`, and the data file `coef_table_n1024.json`).

## Question

The all-distinct post-activation third cumulant κ3(a_l) is, to a few per cent, the sum of seven first-order gate diagrams
(`oracle_k3.residual_basis`, B0..B6) on top of the Gaussian ρ² Wick term. With the leg-partition coefficients
(`CLOSURE_COEF` = 1, 3, 3, 1, 1, 1.5, 1.5) the D21(l+1) error stays at 4–6 % at depth; fitting the seven coefficients per
layer brought it to 1.4–2.3 % at layers 10–14 on one width-128 MLP. If the fitted coefficients are a property of the shape
(n, L) and not of the MLP, a 7 × 15 table fitted offline is a legitimate data file for an estimator. And the frontier needs
ε ≤ 2.2 % at n = 1024 (extra MSE ≈ 4.2e-6 ε²).

## Method

- Atlases: `moment_atlas_np.py --k3 --k4`, depth 16, N = 5e5 per atlas, He-initialised MLPs from seeds 770100+ (one MLP per
  process). Width 64: 6 MLPs, width 96: 6, width 128: 8; at each width 3 MLPs have a second atlas with an independent
  sample seed (the pair, for Monte Carlo noise). (Widths 48, 160, 192 queued.)
- `coef_table.py cache`: per atlas and layer, the transported (n × n) form of every ingredient (slices of κ3(a), leading Wick,
  Gaussian ρ² Wick, the 7 basis tensors and B7 = B6 with the (2,1,1) slice regenerated as u_i C_jk) plus the tensor-space
  normal equations of the residual R = all-distinct κ3(a) − Gaussian ρ² Wick on the basis. Everything downstream is linear in
  the coefficients, so all tables and leave-one-out evaluations run on these caches in seconds.
- Ladder evaluated on each MLP with coefficients from the *other* MLPs of the width (leave-one-out):
  (a) `wick` leading Wick only; (b) `leg` leg-partition closure; (b′) `legR` the same with the regenerated (2,1,1) slice
  (the published chain's closure); (c) `own` per-MLP tensor-space fit (an oracle: fitted on the evaluated MLP), (c′) `ownD`
  the same in D21 space; (d) `ens` ensemble table, tensor space (stacked residuals of the training MLPs per layer); (d′)
  `ensD` ensemble table in D21 space (regress the transported residual on the 7 transported basis terms, i.e. what the chain
  consumes); (e) `ensR` / `ensRD` ensemble table refitted with the regenerated (2,1,1) slice.
- Noise: on the pair MLPs the model built from atlas A is evaluated against atlas B's D21, and two noise energies are
  subtracted: B's target noise (ε_noise = rel(D21_A, D21_B)/√2, as in `oracle_k3.analyse_pair`) and the Monte Carlo noise of
  the model's own inputs (C, κ3(z), the κ4 slice, the slices of κ3(a)), measured as ‖pred_A − pred_B‖/(√2‖D21_B‖) with the
  same coefficients. The result, ε_rep, is the pure representation error. (Without the second subtraction the
  cross-evaluated ε *grows* with width at layers 0–2, because at fixed N the input noise grows with n; a chain does not see
  that noise, so ε_rep is the relevant number.) Layer 0 is excluded from the readings: z_0 is exactly Gaussian, every
  non-Gaussian input is zero in population, and ε_rep(0) ≈ 0.

## Results

### 1. The table is an ensemble property: leave-one-out costs nothing

ε_rep of D21(l+1), mean over the pair MLPs, rms over layer bands:

| width | leg 1-3 | leg 4-9 | leg 10-14 | own 1-3 | own 4-9 | own 10-14 | **ens 1-3** | **ens 4-9** | **ens 10-14** | ensR 1-3 | ensR 4-9 | ensR 10-14 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 64 | 0.045 | 0.074 | 0.093 | 0.038 | 0.046 | 0.029 | 0.038 | 0.046 | 0.031 | 0.080 | 0.092 | 0.061 |
| 96 | 0.036 | 0.049 | 0.068 | 0.033 | 0.029 | 0.028 | 0.034 | 0.030 | 0.029 | 0.074 | 0.066 | 0.056 |
| 128 | 0.031 | 0.038 | 0.055 | 0.030 | 0.025 | 0.018 | 0.030 | 0.025 | 0.018 | 0.065 | 0.058 | 0.039 |

Per layer at width 128 (mean over pair MLPs; [max over MLPs] for own / ens):

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

Readings.
- **The held-out ensemble table equals the per-MLP oracle fit** to ≤ 0.002 at every layer and width (mean and max), in
  tensor space and in D21 space alike (`ens` ≈ `ensD` ≈ `own` ≈ `ownD`). A per-MLP fit buys nothing over an offline table.
- Fitting in D21 space (what the chain consumes) gives the same coefficients' *effect* as fitting in tensor space; the
  tensor-space fit is the cheaper, better-conditioned one and is used for the tables below.
- **The table is also nearly width-independent**: the table fitted at width 64 evaluated on width-128 MLPs loses ≤ 0.001
  ([results/transfer.txt](results/transfer.txt)):

  | evaluated at \ table from | 64 | 96 | 128 | pooled |
  |---|---|---|---|---|
  | 64 (bands 1-3/4-9/10-14) | .038/.046/.033 | .039/.047/.032 | .039/.048/.032 | .039/.047/.031 |
  | 96 | .034/.030/.032 | .033/.030/.030 | .033/.030/.029 | .034/.030/.030 |
  | 128 | .030/.026/.019 | .030/.025/.018 | .030/.025/.019 | .030/.025/.018 |

- A quadratic-in-depth smoothing of the table (3 numbers per coefficient) is as good as the per-layer table.

### 2. Which coefficients matter ([results/ablate.txt](results/ablate.txt))

ε_rep bands at width 128 when one coefficient of the table is reset to its leg-partition value (`ens−k`), or when only one
coefficient is taken from the table and the rest left at leg values (`leg+k`):

| variant | 1-3 | 4-9 | 10-14 |
|---|---|---|---|
| ens | 0.030 | 0.025 | 0.018 |
| leg | 0.031 | 0.038 | 0.055 |
| ens − B2 (w3 D21) | 0.030 | 0.027 | 0.025 |
| ens − B3 (D3 ⊗ C ⊗ C) | 0.031 | 0.034 | 0.047 |
| ens − B6 (κ4 (2,1,1)) | 0.031 | 0.030 | 0.029 |
| ens − B0, B1, B4 or B5 | ≤ 0.030 | ≤ 0.026 | ≤ 0.020 |
| leg + B3 only | 0.030 | 0.031 | 0.030 |

The gain of the table over the leg-partition closure is carried by three numbers per layer: **B3** (the D3 hyperedge with
two C edges on the same vertex, Hermite-degree-5 vertex weight w5) does most of it, then **B6** (the (2,1,1) κ4 hyperedge)
and **B2** (the D21 hyperedge with the C edge on the doubled vertex, weight w3). B0, B1, B4, B5 can stay at their
leg-partition values. The large across-MLP spread of the per-MLP B1, B4, B5 fits (sd up to 0.3–1.0 at depth) is therefore
harmless: they are flat directions of the D21 error.

### 3. Coefficient tables vs depth and width (tensor-space ensemble fit ± across-MLP sd of per-MLP fits)

The three coefficients that matter (full tables for all seven in [results/width*.txt](results/), and per width in
[results/coef_table_n1024.json](results/coef_table_n1024.json)):

| l | B2 n=64 | B2 n=96 | B2 n=128 | B3 n=64 | B3 n=96 | B3 n=128 | B6 n=64 | B6 n=96 | B6 n=128 |
|---|---|---|---|---|---|---|---|---|---|
| leg value | 3 | 3 | 3 | 1 | 1 | 1 | 1.5 | 1.5 | 1.5 |
| 1 | 3.02 ± .09 | 3.05 ± .05 | 3.04 ± .04 | 0.43 ± .08 | 0.44 ± .06 | 0.41 ± .03 | 1.40 ± .02 | 1.43 ± .01 | 1.45 ± .01 |
| 3 | 2.66 ± .12 | 2.80 ± .10 | 2.87 ± .04 | 0.14 ± .05 | 0.18 ± .13 | 0.30 ± .06 | 1.23 ± .03 | 1.30 ± .05 | 1.34 ± .01 |
| 5 | 2.44 ± .08 | 2.55 ± .12 | 2.66 ± .09 | −0.02 ± .09 | 0.00 ± .07 | 0.04 ± .14 | 1.14 ± .04 | 1.21 ± .03 | 1.26 ± .03 |
| 8 | 2.34 ± .08 | 2.51 ± .08 | 2.58 ± .05 | −0.09 ± .12 | −0.12 ± .02 | −0.06 ± .06 | 1.06 ± .05 | 1.14 ± .05 | 1.19 ± .04 |
| 11 | 2.41 ± .07 | 2.47 ± .08 | 2.51 ± .06 | −0.17 ± .05 | −0.14 ± .04 | −0.06 ± .08 | 1.02 ± .06 | 1.09 ± .03 | 1.15 ± .02 |
| 14 | 2.34 ± .12 | 2.46 ± .07 | 2.48 ± .05 | −0.17 ± .14 | −0.11 ± .11 | −0.09 ± .06 | 1.06 ± .07 | 1.13 ± .03 | 1.14 ± .05 |

The pair-to-pair Monte Carlo sd of a single per-MLP fit is 0.00–0.02 for these three, so the across-MLP spread above is real
MLP-to-MLP variation (the flat directions B1, B5 have sd 0.3–1.0 at depth, shrinking from width 64 to 128). B0 = 1.00 at every layer and width.

Drift with width: B2 and B6 move monotonically back toward their leg-partition values as n grows; a fit c(n) = c∞ + a/√n
over 64–128 gives c∞ ≈ 2.7–3.4 for B2 and 1.4–1.6 for B6 at layers 3–14, i.e. **the renormalisation of B2 and B6 is a
finite-width O(ρ) effect that vanishes as n → ∞**. B3 is different: it sits at ≈ 0.4 at layer 1 and ≈ −0.1 at depth at
every width, with only a weak upward drift (c∞ ≈ 0–0.6). The leg-partition value 1 of B3 is not approached; the D3 ⊗ C ⊗ C
diagram with its Gaussian degree-5 vertex weight is strongly over-counted (the w5 vertex weight is the Gaussian closure of
E[relu⁽⁵⁾(z)], a delta-function derivative that is ill-represented by the Gaussian at a non-Gaussian z — a plausible reason,
not tested). Values at n = 1024 from the c∞ + a/√n fit: [results/scaling.txt](results/scaling.txt) and the JSON file.

### 4. Width scaling of ε and extrapolation to n = 1024 (preliminary, widths 64–128)

Fits per layer band of the mean ε_rep: power law ε ∝ n^−p, and a conservative "floor" model ε = a + b/√n
(meaningless when it extrapolates below 0, i.e. when the data fall faster than n^−1/2; marked —):

| model | band | n=64 | n=96 | n=128 | p | ε(1024) power | ε(1024) floor |
|---|---|---|---|---|---|---|---|
| ens | 1-3 | 0.038 | 0.034 | 0.030 | 0.37 | 0.014 | 0.017 |
| ens | 4-9 | 0.046 | 0.030 | 0.025 | 0.91 | 0.004 | — |
| ens | 10-14 | 0.031 | 0.029 | 0.018 | 0.72 | 0.005 | 0.003 |
| leg | 1-3 | 0.045 | 0.036 | 0.031 | 0.53 | 0.010 | 0.010 |
| leg | 4-9 | 0.074 | 0.049 | 0.038 | 0.95 | 0.005 | — |
| leg | 10-14 | 0.093 | 0.068 | 0.055 | 0.76 | 0.011 | — |
| legR | 1-3 | 0.082 | 0.075 | 0.066 | 0.30 | 0.036 | 0.044 |
| legR | 4-9 | 0.108 | 0.077 | 0.066 | 0.72 | 0.014 | — |
| legR | 10-14 | 0.111 | 0.084 | 0.067 | 0.71 | 0.015 | 0.000 |
| ensR | 1-3 | 0.080 | 0.074 | 0.065 | 0.28 | 0.037 | 0.045 |
| ensR | 4-9 | 0.092 | 0.066 | 0.058 | 0.69 | 0.014 | 0.001 |
| ensR | 10-14 | 0.061 | 0.056 | 0.039 | 0.61 | 0.012 | 0.011 |

Per layer, the worst extrapolated `ens` layer is layer 1 (0.033 / 0.031 / 0.029 at 64/96/128, p = 0.19: 2.0 % power law,
2.3 % floor model), then layers 2–3 (1.2–1.6 %).

## Verdict (preliminary, to be re-checked with widths 160 and 192)

1. **Yes, the renormalised coefficients are an ensemble property.** A 7 × 15 table fitted offline reproduces the per-MLP
   oracle fit to ≤ 0.2 points of ε at every layer, and it even transfers across widths at no measurable cost. Only three
   coefficients per layer matter (B2, B3, B6); the others can stay at their leg-partition values.
2. **What remains after the table is the floor of the seven-diagram basis, not coefficient error.** The question "does the
   table reach 2.2 % at 1024" is therefore the question of how the basis floor scales: from 64 to 128 it falls as
   n^−0.4 at layers 1–3 and n^−0.7..−0.9 deeper.
3. **Extrapolated to n = 1024** (three widths, 2× range, so ±50 % on the extrapolated ε): the ensemble table with the true
   (2,1,1) κ4 slice reaches ≈ 1.4–1.7 % at layers 1–3 and ≤ 0.5 % deeper, i.e. **below 2.2 % at every layer except
   marginally layer 1 (2.0–2.3 %)**. With the slice regenerated as u_i C_jk (`ensR`, the only variant a chain can afford
   today) it does **not**: ≈ 3.7–4.5 % at layers 1–3 (≈ 1.2–1.4 % deeper). The shallow layers 1–3 are where the regenerated
   slice fails and where a better carrier of the (2,1,1) slice is required; the table cannot compensate (refitting the
   table with the regenerated slice, ensR vs legR, gains nothing at layers 1–3).
4. Caveat on the extrapolation: the leg-partition closure itself is extrapolated to ≈ 1 % at 1024 (it falls fast at
   depth), so at n = 1024 the table's gain over the un-renormalised closure is expected to be small; the renormalisation is a
   finite-width effect. The published chain's closure (legR) extrapolates to 3.6–4.4 % at layers 1–3, consistent with the
   4–7 % inferred for it at n = 1024 from its raw error, which is the one external calibration point available.

## Open issues / next steps

- Widths 160 (3 MLPs + 1 pair), 192 (2 MLPs + 1 pair) and 48 (6 + 3) running; the extrapolation above rests on 64–128.
- The direct test is the width-1024 oracle on the GPU atlas (keenanpepper's N = 1e9 joint moments give D21, K22; the (2,1,1)
  slice needs a dedicated job): evaluate `legR`, `ensR` at layers 1–3 against the 2.2 % line.
- B3's strong renormalisation (1 → ≈ 0) at every width suggests the Gaussian degree-5 vertex weight is the wrong object;
  replacing w5 by a measured gate statistic is a candidate for a parameter-free closure.
- Layer-1 is the binding layer for `ens`; the residual after the table there is not yet characterised (hub structure,
  missing second-order diagrams).
