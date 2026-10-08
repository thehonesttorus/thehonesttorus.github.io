# The leaderboard system: renormalized closures, judged in adjusted MSE

Working note XXXIX. It assembles the best estimator built so far and then asks one theoretical question of it. Every
number here is reported as the competition scores it: adjusted MSE = raw MSE x max(0.1, C/B), with B = 2^41, measured
in the scored regime (a warm-up predict on another network, then the measured predict, `code/run_v29w.py`).

## 0. What this note does

1. **Phase A: the assembled system.** The V35 chain of note XXIX (section 12; adjusted 4.1273e-9 on all 100 networks)
   with the output-metric counterterms of `../ray-compiler` section 8 (the 103 per-layer amplitudes, `V47_CAL`). Those
   amplitudes were fitted on the V33 chain; this note measures their transfer to V35 and refits them on V35.
2. **Phase B: the renormalized comparison of closures** (section 2). Every derived term added to the chain so far has
   made the free-running output worse while improving its own statistic: GC1 at the theorem's weight (+2.3%), the
   G D core (+9%), the derived kappa4 shapes (+9% to +50%), the physical (2,2) slice (+5.5%). Note XXXVI section 3k
   read the fitted lam core as a compensating counterterm. If that reading is right, those tests counted the same
   physics twice. The fair test renormalizes both sides.

## 1. Phase A: predictions (stated before the runs)

- **Transfer.** V35 computes what V33 computes except the alpha <= -2.5 drop (neutral) and the L-2 join (-0.43%).
  The V33-fitted amplitudes should transfer: on the held-out networks 50-99, V35 + cal33 against V35, raw **-3.0% to
  -3.8%** (V33 measured -3.76%), C/B unchanged (about 0.180), so adjusted by the same percentage.
- **Refit.** Amplitudes refitted on V35's own errors (networks 0-49) with the stored V33 responses: the predicted
  held-out change within 0.5 points of the transfer's measured change.

## 2. The renormalized comparison

**Counterterms.** A closure omits part of the dynamics; its fitted elements absorb part of what is omitted. That is
what a counterterm does, and the chain's fitted elements behave like counterterms: the lam core, METRIC_C = 2 and the
feed weights were fitted to the chain's own output, not to the slices they write (notes XXI, XXXVI section 3k).

**Double counting.** When a derived term T, a piece of the omitted dynamics, is added while the counterterms keep
the values fitted without it, the physics T carries is counted twice: once explicitly and once in the counterterm.
A clean free-running loss is then expected even when T is right. The ledger of note XXXVI section 3h shows exactly
this signature for GC1 (y2): the repair shrinks its own routes (-5%) and loses cross-route opposition (+15%).

**The fair comparison** renormalizes both sides by one matching condition. The counterterm family is the 103
per-layer amplitudes of the statistics the Wick stage reads (var, coff, D3, D21, g4, K22, K31). They are fitted in
the output metric on the training networks 0-49, for the base and for base + T alike. The two are compared on the
held-out networks 50-99.

**First-order implementation** (`code/refit.py`). The base's linear responses R_j to the 103 amplitudes are stored
(10,300 runs, `../ray-compiler` section 8). To first order, out(X + cal) = out_X + R b. So one free-running run of
X per network suffices: b_X = argmin sum_train |e_X + R b|^2, with the ridge chosen inside the training set.
The instrument has predicted two free-running changes to within their noise (-0.28% predicted and measured; -3.07%
predicted, -3.76% measured).

Two further outputs:
- **The free amplitude a of T.** The direction out_X - out_base is appended to the family and fitted with it; a = 1 is
  the derived weight.
- **A synthetic check.** A candidate whose whole effect lies inside the counterterm span gets a = 0 and no
  renormalized gain (`code/refit.py`, checked on synthetic data before any run).

**Why the flex theorem sorts the candidates.** Note XXXVIII proved that the pair state fixes the omitted (2,1,1) and
(1,1,1,1) classes only modulo an exact kernel F. Every closure of those classes from pair statistics is therefore a
convention, a point of F: the lam core, the G D core KD, the derived scale-mixture shapes. A convention carries no
information that the pair state lacks. A law term carries information: a term of the derived next-order kappa3 law,
such as GC1, or a part of the carried state that the representation truncates, such as the birth M-block residual.
The renormalized comparison should separate the two classes.

## 3. Phase B: candidates and predictions (stated before the runs)

All candidates run on V35, networks 0-99, scored regime, outputs saved.

| tag | switch | class | earlier free-running result |
|---|---|---|---|
| y2 | `V40_FB_SY=2` | law: GC1 at the theorem's weight | +2.3% +- 0.5 (32 nets) |
| x2y2 | `V40_FB_SX=2 V40_FB_SY=2` | law: GC1 + GC2 | one-step equal to y2 |
| rres16 | `V17_R_RES=16` | state: birth M-block residual rank 4 -> 16 | -1.38% +- 0.41, +2.8% FLOPs |
| kd | `V39_KD=1` | convention: G D core in three slices | +9% (nets 0-1), +8% FLOPs |
| pmetric | `V44_PMETRIC=1` | convention: physical-metric lam shapes | +0.28% +- 0.18 |
| wk4m | `V33_WK4M=1` | convention: scale-mixture (2,2) slice | +4.8% +- 0.7 (16 nets), +5.5% +- 0.5 (32) |
| wk4m3 | `V33_WK4M=3` | convention: geometric-mean (2,2) slice | -0.61% +- 0.42 (16 nets) |
| k4d3 | `V31_K4D=3` | convention: derived dropped-class diagonal (gain from kappa4) | not run on V33+ |
| k4d4 | `V31_K4D=4` | the same, gain from kappa3 | not run on V33+ |

**Predictions.** Delta is the renormalized held-out change: refitted X against refitted V35, paired.
- **P1 (law term survives).** For y2, Delta < 0 by more than 2 standard errors, and a lies in [0.5, 1.5].
- **P2.** For x2y2, Delta equals y2's within one standard error.
- **P3 (conventions do not).** For kd, pmetric, wk4m, wk4m3, k4d3 and k4d4, Delta >= -0.5% each.
- **P4 (control).** For rres16, Delta lies within 0.5 points of its own free-running change. Accuracy that comes
  from carrying more of the state is neither created nor destroyed by renormalization.

If P1 fails, the reading of section 2 is wrong for the kappa3 law: y2's loss is not double counting against the
amplitude family, and the opposition it removes lies outside every per-layer amplitude.

**Amendment (before any held-out analysis).** A pilot inside the training set showed that each candidate's own
ridge choice is noisy. It was fitted on networks 0-24, judged on 25-49, and the base refit recovered only -0.20% at
that size. With the ridge chosen per candidate, a candidate compared at ridge 0 against a base at ridge 0.1 measures
the regularization, not the term. The comparison is therefore made in one counterterm family with one
regularization: every candidate, every free-amplitude fit and the joint fit use the base's ridge. The candidate's own
choice is printed alongside. The held-out networks 50-99 had not been analysed when this was fixed.

**Decision (adjusted MSE).** T is adopted only if the measured free-running V35 + T + cal_T beats V35 + cal_35 on the
held-out networks in adjusted MSE, by more than two standard errors, cost included.

**The D21 component audit, asked of the estimator itself.** A parallel line of work (the capture-audit and
depth-profile checkpoints supplied with this round) measured the physical D21 table. It splits into additive parts
(row and column effects; captured 0.998 at rank 2) and spectrally flat zero-row-sum parts, which respond at every
depth. Its gating contract asks what the flat parts are worth at the chain's own compression point.

In V35 the full D21 table is used in the covariance programs. Only the feedback into the birth factors (V18) is
compressed, to rank 2, which keeps the additive part. So the question has a direct estimator form: the feedback at
rank 16 and at full rank (exact), against rank 2.

| tag | switch | prediction |
|---|---|---|
| fb16 | `V18_R_FB=16` | raw -1% to -3%; C/B +10% to +15% (the thin legs carry 2r columns per source); adjusted worse |
| fbx | `V18_R_FB=1024` | raw -1.5% to -4%; adjusted much worse |

The raw range is the precision law of note XXXI applied to the one-step legs-error change at rank 16 (section 3e of
note XXXVI). If the exact feedback's raw gain is below its cost share, no response-aware rank-r scheme can pay: the
exact feedback bounds every such scheme. At rank 2 the attached proposal reduces to the exact additive part in place
of the top-2 range finder, which changes the feedback by at most the non-additive energy in the top two directions.

**The attached proposal at rank 2** (`V49_FB_ADD=1`, added after the audit of section 5, before any run of it). At the
production rank the proposal reduces to the exact additive part of D21, from its row and column sums in O(n^2), in place
of the top-2 range finder. The audit says the best rank 2 already keeps 99.8-100% of the additive part, plus 1-34% of
the flat part, which the exact version drops. Prediction, V35 + fbadd on networks 0-99: raw within +-0.5% of V35; C/B
about -0.1% (no range finder); renormalized change within +-0.5%.

**Strassen leaf minimum 8 for the whole chain** (task of note XXIX section 10). `V26_STRASSEN_MIN=8` on networks
0-15: raw unchanged to rounding; C/B about -2.4%; it is adopted if the measured local residual of the measured call
stays within 10% of V35's and VmHWM within 1 GB.

## 4. Phase A, measured

Scored regime, V35 on all 100 networks (`outputs/pa_rows.json`): raw 2.2899e-8, C/B 0.1802, adjusted **4.1273e-9**. This
reproduces note XXIX section 12 to the last digit.

| held-out networks 50-99 | raw | C/B | adjusted |
|---|---|---|---|
| V35 | 2.2670e-8 | 0.1803 | 4.0882e-9 |
| V35 + cal33 (the V33-fitted amplitudes) | 2.1807e-8 | 0.1804 | **3.9333e-9** |
| change | -3.80% (per net -3.77 +- 0.37, better on 46/50) | +0.0% | **-3.79%** |

- **The transfer prediction holds.** The measured -3.80% is at the edge of the predicted -3.0% to -3.8%.
- **The refit misses its band.** Fitted on V35's own errors (networks 0-49, `outputs/refit_v35.txt`), the refit
  predicts -3.16% +- 0.37 held-out, with coefficients that correlate 0.9992 with cal33's. That is 0.64 points short
  of the measured transfer, outside the predicted 0.5-point band.
- **The bias is systematic.** On V33 the first-order prediction was -3.07% against -3.76% measured. In both cases the
  second order is favourable by about 0.6-0.7 points.
- **For submission.** The amplitudes are refitted on all 100 networks; the held-out figure above is the estimate of
  what that buys.

## 5. The five-component audits of the chain's own slices (`code/d21audit.py`, `outputs/*audit_off{0,1}.txt`)

The attached capture audit split the physical D21 table of one cut into additive parts (S0, S1, A1) and flat
zero-row-sum parts (S2, A2). Here the same exact orthogonal split is applied to three tables:
- the chain's own pre-activation slices at every layer (free-running dumps, networks 0 and 1);
- the Monte Carlo truth (1.6e7 inputs);
- the chain's error.

Each error is noise-corrected with the two independent halves: the noise energy of the full sample is
|comp(h0 - h1)|^2 / 4. The noise is 1-2% of the error at depth and 10-40% at layers 2-5.

**D21 = kappa(z_i, z_i, z_c).**
- **The physical structure holds in the chain's tables.** At every layer, on both networks:
  - the truth is additive plus flat, with S1 = A1 and S2 = A2 to three digits;
  - the flat share falls with depth, from 0.69 at layer 1 to 0.18 at layer 14;
  - the flat parts are spectrally flat. Their rank-16 capture runs from 0.10 at layer 1 to 0.55 at layer 14, with the
    S2 and A2 curves identical (the attached checkpoint measured 0.197 at its layer-5 cut; here 0.235).
- **The chain's error is flat.** Noise-corrected, the error energy is 0.0005 of the truth's at layer 2, rising to
  0.0050-0.0054 at layer 14 (7% rms). 90-93% of it is in S2 + A2.
- **The additive part is accurate.** Measured against each sector's own energy:
  - the additive part (row and column effects) is right to about 2% rms;
  - the flat part is off by about 16% rms at depth.
- **What the flat part is.** It is the incoherent quenched sector that the product-gate transport defect of note XV
  and the closure round feeds. Sums over rows or columns average an incoherent defect away, so the additive part
  survives it.
- **The compression.** The production feedback compresses D21 to rank 2. The best rank 2 keeps 99.8-100% of the
  additive part and 1-34% of the flat part.

**K31 = kappa(z_a, z_a, z_a, z_c)** (the chain's lam C_off against the truth's (3,1) slice transposed, as the V37
oracle uses it).
- **The truth is flat, with an antisymmetric part.** S2 is 0.71-0.88 of its energy and A2 is 0.12-0.29; the additive
  share is 0.002. The symmetric part becomes low-rank with depth: rank 16 holds 0.24 of S2 at layer 2 and 0.85 at
  layer 14.
- **The chain's slice is flat and symmetric.** Its S2 share is 0.998, so the antisymmetric part is missing entirely.
- **The error is large.** Noise-corrected, 0.53-0.75 of the truth's energy: about 74% rms at every depth.
- **The calibration was shrinking a noisy estimate.** The calibration's K31 amplitudes at depth (0.50-0.84) are the
  Wiener shrinkage this error implies: 1/(1 + 0.55) = 0.65.

**K22 = kappa(z_a, z_a, z_b, z_b)** (the chain's (dG_a + dG_b)/3 against the truth).
- **The truth is mostly a constant.** It is S0 0.86-0.97, S1 0.01-0.10 and S2 0.02-0.05.
- **The chain's shape is additive by construction.** It misses the flat part entirely.
- **The error grows with depth.** It is 0.02 of the truth's energy at layer 2 and 0.09-0.11 at layer 14 (32% rms).
- **The output barely reads it.** Its route share in the first-entry ledger of note XXXVI is 0.016 / -0.007.

**What K31 resembles** (`code/k31struct.py`, `code/k31rows.py`). Explained shares of the noise-corrected K31
energy, from the truth's own state:

| candidate | layer 2 | layer 8 | layer 14 |
|---|---|---|---|
| C_off, one coefficient (the chain's shape at the optimal amplitude) | 0.49 / 0.45 | 0.51 / 0.48 | 0.58 / 0.52 |
| C_off, one coefficient per row | 0.50 / 0.45 | 0.52 / 0.50 | 0.62 / 0.56 |
| C var_c (the physical-metric shape), one coefficient | 0.50 / 0.45 | 0.51 / 0.49 | 0.63 / 0.56 |
| C^2 | 0.34 / 0.31 | 0.44 / 0.42 | 0.61 / 0.57 |
| row-scaled flat D21, including the gain mode's 3 var_a / mu_a | 0.02 / 0.02 | 0.09 / 0.08 | 0.23 / 0.19 |
| C and D21, two coefficients per row | 0.50 | 0.56 / 0.54 | 0.69 / 0.63 |
| the same with both coefficients cubic in alpha (8 global numbers) | 0.50 | 0.55 / 0.53 | 0.67 / 0.61 |

(network 0 / network 1.)
- **The relation to D21 is weak.** The D21 row coefficient tracks the unit's own state: correlation +0.85 to +0.97
  with alpha and with its skewness, layers 4-15. It adds only 1-7 points to the covariance shape.
- **No function of the pair state closes K31.** About 40% of its energy is not a function of the pair state at all,
  as the flex theorem (note XXXVIII) requires: the (3,1) slice of the next layer is fed by the omitted classes.

**Which basis concentrates K31** (`code/k31basis.py`, `outputs/k31basis_off{0,1}.txt`; the open question of note
XXXVIII section 9). Captured share of the noise-corrected K31 energy, two-sided projection, at K = 16 (net 0 / net 1):

| basis | layer 2 | layer 8 | layer 15 |
|---|---|---|---|
| K31's own SVD (oracle) | 0.29 / 0.29 | 0.58 / 0.60 | 0.85 / 0.84 |
| top covariance eigenvectors (test S's choice) | 0.12 / 0.11 | 0.37 / 0.37 | 0.66 / 0.61 |
| D21's right singular vectors on both sides | 0.11 / 0.11 | 0.39 / 0.39 | 0.68 / 0.64 |

- **No carried basis comes close.** The best one gains 0.02-0.03 on the covariance eigenvectors at depth. K31's own
  basis would change the cost argument of note XXXVIII, but it is not computable from the carried state.
- **The verdict stands.** The carried-symbol state remains closed for cost.

**Reading.**
- **The sectors fail differently.** In the kappa3 sector the chain transports the coherent parts to about 2% and
  the incoherent flat parts to about 16%. In the kappa4 sector its slices are closures, 32-74% off, and limited by
  the flex, not by transport.
- **One wall under both.** The product-gate correction attaches the correlation matrix to a pair of legs. The
  kappa3 -> kappa4 feed that drives the omitted classes attaches the covariance to one leg. Either way every CP
  component becomes rank n, so the exact term costs n^4 per source-layer (note XV, closure round section 5).
- **For the attached proposal.** Exact additive storage plus a response-aware rank for the flat residual addresses a
  compression loss. The audit puts the chain's flat-sector loss in transport. Phase B's feedback-rank runs price the
  compression side directly.

## 6. Phase B, measured: renormalization repairs the conventions, not the law term (`outputs/refit_phaseB_heldout.txt`)

Train on networks 0-49, judge on 50-99, one ridge for all (0, the base's choice). The base refit predicts -3.16% +- 0.37
held-out. Changes are of raw MSE, paired per network.

| candidate | class | free-running, held-out | renormalized change | better on | free amplitude a |
|---|---|---|---|---|---|
| y2 | law (GC1) | +2.15% +- 0.46 | **+2.40% +- 0.39** | 7/50 | 0.20 |
| x2y2 | law (GC1 + GC2) | +3.09% +- 0.47 | +3.30% +- 0.39 | 4/50 | 0.13 |
| rres16 | state (control) | -1.04% +- 0.26 | -1.37% +- 0.23 | 41/50 | 0.91 |
| kd | convention | +8.40% +- 0.68 | +5.09% +- 0.41 | 2/50 | 0.11 |
| pmetric | convention | +5.76% +- 0.65 | +5.71% +- 0.60 | 3/50 | 0.23 |
| wk4m | convention | +5.19% +- 0.37 | **-0.02% +- 0.10** | 25/50 | 0.50 |
| wk4m3 | convention | -1.09% +- 0.23 | -0.23% +- 0.08 | 35/50 | 1.19 |
| k4d3 | convention | -0.10% +- 0.26 | -0.41% +- 0.21 | 31/50 | 0.47 |
| k4d4 | convention | +2.51% +- 0.51 | +1.25% +- 0.42 | 16/50 | 0.28 |

**Verdicts on the predictions of section 3.**
- **P1 fails.** GC1's renormalized change, +2.40% +- 0.39, equals its free-running change, and its free amplitude is
  0.20. The loss is not double counting against the counterterm family.
- **P2 fails.** x2y2 is 0.9 points worse than y2: doubling the Xt (GC2) feedback costs on top of GC1.
- **P3 holds for all six conventions.** None gains more than 0.5% once renormalized.
- **P4 holds.** The control's renormalized change, -1.37%, is within 0.33 points of its free-running change, with
  a = 0.91.

**The sorting is the reverse of section 2's hypothesis for the law term.**
- **Double counting is real, for the conventions.** Renormalization repairs exactly what the conventions break:
  - wk4m's +5.19% free-running loss vanishes (-0.02% +- 0.10);
  - kd's shrinks by 40% and k4d4's by half;
  - k4d3 turns slightly favourable.

  That is the counterterm reading of note XXXVI section 3k, now measured: the lam core and the other amplitudes had
  absorbed what these shapes add.
- **No convention gains.** None of them gains anything beyond the counterterms, as the flex theorem requires: a
  closure of the omitted classes from pair statistics carries no information the pair state lacks.
- **The law term's loss survives.** It is the signed cancellation of note XXXVI section 3h, and section 5 gives it a
  location. The kappa3 errors y2 shrinks sit in the incoherent (flat, per-neuron) sector. So do the kappa4 errors
  they had offset: K31 74% off, g4 30% per neuron. A per-layer amplitude cannot represent a cancellation between
  incoherent components.
- **A paired repair is needed.** A law-term repair of the kappa3 sector needs the kappa4 sector's incoherent content
  carried at the same order: the flex coordinate, or the n^4 gate-covariance terms of section 5.

**The joint free fit (a diagnostic, not adopted).** All nine directions are appended to the counterterms at once:
-2.69% +- 0.29 held-out, better on 47/50. The amplitudes:

| rres16 | wk4m3 | wk4m | y2 | x2y2 | k4d4 | k4d3 | kd | pmetric |
|---|---|---|---|---|---|---|---|---|
| +0.82 | +1.69 | -0.70 | +0.72 | -0.60 | +0.38 | -0.27 | +0.19 | +0.17 |

- **What the amplitudes say.** The output metric wants the geometric-mean (2,2) shape pushed past itself, the rank-16
  residual (which costs +3.2% C/B), and the Xt (GC2) feedback at 0.4 of its present weight, against the theorem's 2.
- **Why it is not adopted.** These are fitted amplitudes of derived directions, the same status as the 103
  counterterms, and none was a pre-registered candidate. Adopting them would mean switching from derivation to
  calibration.

**Adoption candidates under the rule of section 3.** wk4m3 and k4d3 cost nothing and pass the screen. wk4m3's
geometric-mean slice sqrt(g4_i g4_j)/3 is the gain mode's var_i var_j term exactly, and the gain mode is 93% of the
truth's kappa4 at depth. Both, and the pair, are measured free-running on the held-out networks next (section 7).
