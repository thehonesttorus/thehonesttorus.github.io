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

**What it found** (details in sections 4-9):
- **The assembled system.** V35 plus refitted counterterms is the measured best: adjusted 3.93e-9 on the held-out
  networks, 3.8% below V35 at unchanged cost (section 8; the closing 100-network figure is in section 9).
- **Renormalization sorts the derived terms as the flex theorem says.** It repairs every closure convention's
  free-running loss and lets none of them gain. It does not rescue the law term GC1 (section 6).
- **The audits localize the error by sector.** The chain's D21 error is the incoherent flat sector of the table,
  right to 2% in its additive part and off by 16% in its flat part. Its kappa4 slices are flex-limited closures,
  K31 74% off. No carried basis concentrates K31 (section 5).
- **The exact D21 feedback is the round's one large accuracy lever.** It gives -8.6% raw, three times the
  prediction. Rank 64 already holds 73% of that. Its price sits on the same frontier as the tiers: about one percent
  of raw per percent of cost (section 7).
- **Its cost was the implementation's, and the cost audit removes it (section 11).**
  - **The fold.** To first order the exact feedback is a displacement of the newborn's arm, a -> a + Yt/(2 w2) + Xt/6.
    It costs n^2 per birth and needs no new leg.
  - **The constraint.** Because the arm also carries the birth M block's separable part, the M residual must
    absorb the change. Without that correction the fold is +140%.
  - **The new leaderboard system** (V35 + pair + fold + counterterms, section 11h): adjusted 3.4923e-9 on all 100
    networks, against 3.8815e-9 for section 9's. The protocol's held-out comparison is -9.97% +- 0.58, better on
    50/50, at 1.4% less cost.

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

## 7. The exact D21 feedback: the prediction fails by a factor of three (`outputs/fbx_cold_16nets/`)

`V18_R_FB=1024` (the feedback range finder at full rank, so Qf Bf = D21 exactly) failed the scored harness: it overruns
the 2^41 budget on the warm-up predict. So it was measured in the cold harness (`code/oracle_one.py`, budget 2^44),
paired against V35 in the same harness, networks 0-15:

| | raw | per net | better on | FLOPs (cold) |
|---|---|---|---|---|
| V35 | 2.3520e-8 | | | 0.198 B |
| exact feedback | 2.1501e-8 | **-8.59% +- 0.70** | 16/16 | 1.080 B |

- **Outside the predicted range.** The prediction was -1.5% to -4%: the precision law applied to the one-step legs
  error at rank 16. The measured gain is three times larger, so the prediction is wrong in kind, not in precision.
  The rank-16 legs-error change could not see it. Rank 16 holds almost none of the flat sector (section 5: 1-34%,
  energy-ordered), and the feedback's value lives there.
- **The attached claim holds, at a different point.** The flat parts respond, as the attached checkpoints said, but
  the place they matter is the feedback into the births (GC1, GC2 and the Gamma x Gamma case, note XXXVI section 3f),
  not a compression of the table.
- **It is the incoherent half of the next-order law.** At rank 2 the feedback carries only the additive (coherent)
  part of D21, and doubling that part (y2) is adverse even after renormalization (section 6). At full rank the
  feedback adds the flat (incoherent) part of the same Gamma x C terms, which is where section 5 located the chain's
  kappa3 transport error.
- **The cost is the implementation's.** The rank-1024 thin legs are dense, transported for every source at every
  layer and never confined: 5.5 times the bill. An efficient form must carry the feedback in the main legs, where
  the tiers compress it.

**Pre-registered next (16 networks, cold harness, against V35 in the same harness):**
- **P5 (theorem weight).** GC1 at the theorem's weight on the exact feedback (`V18_R_FB=1024 V40_FB_SY=2`) improves on
  the exact feedback at half weight by 1 to 5 points. At rank 2 the same doubling was adverse.
- **P6.** Adding GC2 at the theorem's weight (`V40_FB_SX=2` as well) is within 1 point of P5's run.
- **P7 (rank).** The gain grows with the feedback rank roughly as the captured flat energy: rank 64 gives 40-50% of
  the exact feedback's gain, rank 256 gives 75-90%.

**The feedback at the production point and at rank 16** (scored regime, 100 networks, paired against V35;
`outputs/pb_rows.json`, `outputs/fbadd_rows.json`):

| | raw | C/B | adjusted | worse on |
|---|---|---|---|---|
| fbadd (exact additive part in place of the top-2 range finder) | +0.70% +- 0.06 | -0.05% | **+0.64%** | 86/100 |
| fb16 (range finder at rank 16) | -2.92% +- 0.10 | +6.44% | **+3.32%** | 100/100 |
| exact (rank 1024, cold harness, 16 networks) | -8.59% +- 0.70 | x5.5 | | 16/16 |

- **The attached proposal at rank 2 is adverse.** The top-2 range finder keeps some flat content (section 5: 1-34%),
  and dropping it costs 0.70%. That is just outside the predicted +-0.5%, in the direction the exact feedback
  explains.
- **The gain grows steadily with rank:** -0.7% for the additive part only, -2.9% at rank 16, -8.6% exact. This is
  the signature of the flat sector, which no rank truncation concentrates.
- **The cost grows with rank too.** The thin feedback legs are carried, unconfined, by every source: +6.4% C/B at
  rank 16 and x5.5 at full rank. Every rank tried loses in adjusted MSE.
- **What the exact feedback is.** It is the GC1/GC2 transport of every existing source, collected into the newborn:
  the arm of centre u is row u of D21. Its information is the whole n x n table at each birth. An efficient exact
  form therefore needs one more dense leg per young source: about double the young tier, roughly +57% of the bill
  against -8.6% raw.

**Strassen leaf minimum 8** (`V26_STRASSEN_MIN=8`, networks 0-15). On V35 it changes nothing: raw and C/B are identical
to the last digit, and the residual is within noise (-1.2% +- 0.5). V35's own leaf minima (`V35_SB_MN`, `V35_CPRE_MN`,
`V35_JN_MN`) already took the saving the V33 measurement had found. Not adopted, and nothing is lost.

**P5-P7, measured** (cold harness, networks 0-15, paired against V35 in the same harness; `outputs/fbx2_cold_16nets/`):

| feedback | raw | better on | FLOPs |
|---|---|---|---|
| rank 64 | -6.26% +- 0.40 | 16/16 | +26% |
| rank 256 | -8.34% +- 0.62 | 16/16 | +108% |
| exact (rank 1024) | -8.59% +- 0.70 | 16/16 | +446% |
| exact, GC1 at the theorem's weight | -7.85% +- 1.37; against exact: +0.76 +- 1.03 | 8/16 vs exact | +446% |
| exact, GC1 and GC2 at the theorem's weight | -7.27% +- 1.47; against exact: +1.38 +- 1.07 | 6/16 vs exact | +446% |

- **P5 fails.** The theorem's weight does not help even with the flat part present. The half weight the hub has
  always carried stays.
- **P6 holds**, weakly: GC2 adds 0.6 points of loss, within 1 point.
- **P7 fails in the favourable direction.** Rank 64 already holds 73% of the exact feedback's gain (40-50% was
  predicted) and rank 256 holds 97%. The gain concentrates faster than the flat energy does: its value is in the
  upper part of the flat spectrum.
- **The adjusted verdict.** At the margin the feedback trades raw for cost at about one for one:
  - rank 16: -2.9% raw for +6.4% C/B;
  - rank 64: -6.3% raw for +26% FLOPs;
  - rank 256: -8.3% raw for +108% FLOPs.

  The thin legs are carried by every source, unconfined, and are re-formed densely at every layer: about 8 n^2 r
  per source-layer.
- **What a winning version would need.**
  - The raw/cost ratio of the feedback must beat the frontier on which the tiers sit (elasticity 1, note XXX
    section 7).
  - The cheapest variants (rank 64 at ages 1-2, rank 32 at ages 3-4, none once confined) estimate at about -4.5% raw
    for +6% C/B: still about neutral.
  - The feedback is the one measured accuracy lever of this round. Its cost structure, not its information, is now
    the obstacle.

## 8. Adoption tests (free-running, held-out networks 50-99, each system with its own counterterms refitted on 0-49)

| system | raw | C/B | adjusted | against V35 + counterterms |
|---|---|---|---|---|
| V35 (no counterterms) | 2.2670e-8 | 0.1803 | 4.0882e-9 | |
| V35 + counterterms (refit) | 2.1799e-8 | 0.1804 | **3.9317e-9** | (-3.83% against V35; predicted -3.16%) |
| + geometric-mean (2,2) slice (wk4m3) | 2.1761e-8 | 0.1804 | 3.9255e-9 | -0.14% +- 0.08, better on 31/50 |
| + derived kappa4 diagonal (k4d3) | 2.1713e-8 | 0.1804 | 3.9169e-9 | -0.39% +- 0.22, better on 31/50 |

- **They match the screen.** Both agree with their renormalized predictions (-0.23 and -0.41).
- **Neither passes the rule.** Each is about 1.8 standard errors, under the two of section 3, so neither is adopted
  alone. The pair is tested next.
- **The first-order refit is pessimistic again.** It predicted -3.16% and measured -3.83%, the third time the
  second order has been favourable by 0.6-0.7 points.
- **The pair** (wk4m3 + k4d3, run free on networks 0-99, `outputs/refit_pair_heldout.txt`). Free-running it is -0.24%
  +- 0.28 held-out. Renormalized it is -0.88% +- 0.22, better on 35/50, slightly better than the two singles summed
  (-0.64%), with free amplitude 0.52.
- **The pair's test and the final runs.** The pair is measured free-running on 50-99 with its own refitted
  counterterms. In the same batch both possible final systems are baked and run on all 100 networks: V35 +
  counterterms, and V35 + pair + counterterms.


## 9. The leaderboard system (superseded by section 11h)

**Adopted** (each part under the rule of section 3, measured free-running on held-out networks). V35, plus:
- the geometric-mean (2,2) slice, `V33_WK4M=3`: wk4m = sqrt(g4_i g4_j)/3, the gain mode's var_i var_j term;
- the derived dropped-class kappa4 diagonal, `V31_K4D=3`: 6 g s_diag^2 s_off^2 + 3 g s_off^4, with the mixture gain
  read from the chain's own kappa4 diagonal, in place of the lam s_off^2 stand-in;
- the 103 output-metric counterterms, refitted on all 100 networks for this configuration.

The file is `code/estimator_final_pair.py` (it was `estimator_final.py` until section 11h replaced it), baked from
`est_v29.py` with `code/bake.py` and `code/estimator_final_pair_config.txt`. It runs with no environment, and a baked
file reproduces its environment-driven
run to the last digit (network 50: raw 1.6962e-8, C/B 0.1819, on two different machines).

| scored regime, all 100 networks (`outputs/pd_rows.json`) | raw | C/B | adjusted | against V35 |
|---|---|---|---|---|
| V35 (note XXIX) | 2.2899e-8 | 0.1802 | 4.1273e-9 | |
| V35 + counterterms | 2.1628e-8 | 0.1803 | 3.8986e-9 | -5.54% +- 0.27 |
| **V35 + pair + counterterms (adopted)** | **2.1526e-8** | **0.1803** | **3.8815e-9** | **-5.95% +- 0.34** |

- **The 100-network figures are in-sample for the counterterms.**
- **The held-out estimate is the protocol's.** Counterterms fitted on networks 0-49 and judged free-running on
  50-99 give adjusted 3.8984e-9 against V35's 4.0882e-9: -4.64% (per net -4.59 +- 0.44).
  - The counterterms alone are -3.83%.
  - The pair adds -0.85% +- 0.23 (better on 33/50), against a predicted -0.88%.
- **Residual and memory are V35's.** The local residual of the measured call is 0.61 s mean and 0.68 s max on a
  loaded fleet machine, against V35's 0.59 / 0.68 in the same harness. The grader counts only our own Python
  between ops (note XXIX section 8). VmHWM is 10.36 GB in the in-process harness.

**What it is not.** It is not the breakthrough this round set out to find.
- The leaders sit at adjusted 1.1-2.1e-9; this system is at 3.9e-9.
- The round located why, sector by sector (sections 5-7):
  - the kappa3 error is incoherent flat-sector transport;
  - the kappa4 slices are flex-limited closures;
  - the one large lever found, the exact D21 feedback (-8.6% raw), is information the chain can use but cannot yet
    afford.

## 10. What follows

1. **The feedback's carrier is the design problem.**
   - The information is real: the full D21 table at each birth, worth -8.6% raw.
   - The present representation re-forms rank-r thin legs densely for every source at every layer, about 8 n^2 r
     multiply-adds per source-layer. At the margin it buys raw at about one for one with cost:
     - rank 16: 0.45% raw per 1% cost;
     - rank 64: 0.24;
     - young-only at rank 64: about 0.5 (estimated).
   - To beat the frontier, a carrier must hold the feedback's upper flat spectrum at well under n^2 r per
     source-layer, for instance without re-forming dense legs for the Hadamard reads. That is the one lever this round
     found with information to spare.
2. **The kappa4 sector needs carried information.** K31 is 74% off and no function of the pair state, and no basis
   computable from the carried state, closes it. The flex theorem now has its quantitative face: about 40% of the
   (3,1) slice is history, not state.
3. **Law-term repairs must be paired at the incoherent level.** GC1 at the theorem's weight fails before
   renormalization, after it, and even with the flat part present (P5). The coherent half-weight the hub carries is
   the one the output tolerates. Any further law term has to come with the kappa4-side term that cancels it at the
   same order.

## 11. The cost audit, and the feedback folded into the arm (stated before the runs)

The question of this section: was any development of this round implemented naively against flopscope's cost rules,
and what does the theory allow at no cost? The rules that matter here:
- matmul (m,k)(k,n) bills 2mkn - mn;
- every written element bills at least 1, and views are free;
- float16 bills like float32;
- transcendentals carry weight 16, and gathers and where weight 4.

At V35's operating point one source-layer of young transport or young hub is about 2n^3. One n^2 elementwise op is
1/2048 of that.

**11a. The four developments, line by line.**
- **The pair (`V33_WK4M=3`, `V31_K4D=3`).** O(n^2) per layer, with two redundancies, now removed bit-identically:
  - the WK4M block formed the C_off^2 class (three n^2 ops) and then discarded it at WK4M = 3;
  - the K4D = 3 diagonal recomputed WW var_prev and WW g_prev (two n^2 matvecs), which the adaptive rule's t_v and
    t_g already hold.

  Together that is about 7 n^2 per layer, roughly 1e-4 B per predict.
- **The counterterms (`V47_CAL`).** Per-layer scalar rescales of the Wick stage's reads: vectors, plus at most four
  n x n arrays per layer. About 4 n^2 per layer; negligible.
- **fbadd (`V49_FB_ADD`).** Row and column sums, O(n^2). Negligible, and adverse anyway (section 7).
- **The exact feedback (`V18_R_FB=1024`).** Naive in kind, not in constant. The rank-r thin legs [F1 | F2] are:
  - transported for every source at every layer, unconfined;
  - re-formed densely (Xt = F1 R1^T, Yt = F2 R2^T) for the Hadamard reads of the hub;
  - contracted again through the thin right factors.

  That is about 8 n^2 r per source-layer, x5.5 the bill at r = n. Section 7 concluded that an efficient exact form
  needs one more dense leg per young source. That is wrong at first order, as 11b shows.

**11b. The fold.** The newborn's hub is B1 = Sym(X1 x P x Y1) with P = I, X1 = 3a + Xt and Y1 = a d(w2) + Yt, where
Xt = 1.5 d(w2) D21 and Yt = 0.5 d(w1) D21^T d(w3). Per centre c it is
3 w2_c Sym(a_c x e_c x a_c) + Sym(a_c x e_c x (3 Yt_c + w2_c Xt_c)) + Sym(Xt_c x e_c x Yt_c). To first order in D21
this is the star with the arm a_c moved to a_c + d_c, where

    d_c = Yt_c / (2 w2_c) + Xt_c / 6     (exact D21, full rank, O(n^2) at birth).

The folded arm rides in the A leg, so the transport, the hub, the joins and both old tiers carry it as they carry any
arm. No new leg and no thin columns are needed, and the rank-2 range finder and its thin legs go.
`V52_FB_FOLD=1` implements it.

**The A leg's second duty.** The A leg also carries the separable part of the birth M block, 3 A d(e) (S_sep = d(e) a^T
is the leading (2,1) Wick term). Folding the arm therefore also moves M by 3 d d(e), at first order. The exact (2,1)
slice read off the Wick table (the ('d21T', (1,2), (2,1), 1/2) and ('d21', (2,2), (1,1), 1/2) rows of pk21, with
pk11's centring) holds two full-rank first-order D21 terms. Until now only the rank-4 residual carried them:

    R1[c,i] = e_c w2_i D21[i,c] / 2,      R2[c,i] = (Phi_c (1 - Phi_c) - mu_c w2_c) Phi_i D21[c,i].

The old sources' gate transport w1^2 D21 w1 lacks both; R2 is about a third of it at alpha = 0. The arm's M-block share
e_c d_c has exactly these two shapes. Subtracting it from the residual before compression (mode 1) keeps M exact up to
the compression, and it shrinks the residual:
- R1 is left at (1/2 - FB_SX/4) of itself, so the arm carries it exactly at the theorem weight FB_SX = 2;
- R2's coefficient becomes g_c - (FB_SY/4) e_c w3_c / w2_c, which is smaller for |alpha| > 0 (for example
  -0.128 -> -0.042 at alpha = 1 and 0.114 -> 0.079 at alpha = -1, at FB_SY = 1), and unchanged at alpha = 0.

So the fold gives the star the exact first-order feedback and gives M part of its exact (2,1) slice, both for free.
What differs from the thin-leg exact feedback:
- **Second order.** The star's Gamma x Gamma terms are (1/2) Xt Yt + (3/(4 w2)) Yt Yt + (w2/12) Xt Xt, against Xt Yt.
- **The K4 -> K3 feed.** Its X3 = A c1 picks up d c1. That is a feed x D21 cross term, about 1e-3 of R2, and M absorbs
  it at birth.

Cost: seven n^2 ops per birth.

**11c. Predictions** (cold harness, networks 0-15, paired against V35 in the same harness, as in section 7):
- **P8 (the fold, `V52_FB_FOLD=1`, weights 1 and 1).** Raw -5% to -10%, better on at least 14 of 16 networks.
  FLOPs -0.5% to -1.5% against V35: the rank-2 legs and their range finder go.
- **P9 (ablation, `V52_FB_FOLD=2`: no residual correction).** Within 1 point of P8. Rank 4 keeps little of a flat
  term either way, so the arm's M-block share does its work with or without the correction.
- **P10 (`V40_FB_SX=2` under the fold).** Better than P8 by 0 to 2 points. The arm then carries R1 exactly in M, and
  that outweighs the star's GC2 overshoot, which cost 0.6 points in the thin-leg exact run (P6).
- **P11 (diagnostics at the exact residual, `V17_R_RES=1024`; not candidates).**
  - V35 at the exact residual: -3% to -8%. R_RES 64 gave -2.2% (note XXXVI section 3l), and the gain grows with rank
    as a flat term's does.
  - The fold at the exact residual beats both the fold and V35 at the exact residual.

**Then, in adjusted MSE** (scored regime, all 100 networks; the protocol of sections 3 and 8). The better of P8 and P10
goes on the adopted system (V35 + pair), with its counterterms refitted on networks 0-49 and judged free-running on
50-99.
- **P12.** Held-out adjusted 3% to 8% below the adopted system's 3.8984e-9.
- **The counterterms cannot mimic it.** They are per-layer amplitudes, and the fold's information is the flat
  sector of D21.

Adoption follows the rule of section 3: held-out, more than two standard errors, cost included.

**11d. P8-P11, measured** (cold harness, networks 0-15, paired against V35 rerun in the same batch;
`outputs/fold_cold_16nets/`). The V35 rerun reproduces section 7's V35 to 0.04% (cross-hardware BLAS rounding).

| variant | raw | better on | FLOPs | prediction |
|---|---|---|---|---|
| fold (`V52_FB_FOLD=1`) | **-5.76% +- 1.20** | 14/16 | **-1.25%** | P8: -5% to -10%, >= 14/16, FLOPs -0.5% to -1.5%: **holds** |
| fold, GC2 at the theorem weight (`V40_FB_SX=2`) | -2.35% +- 1.69 | 11/16 | -1.25% | P10: 0-2 points better than P8: **fails** (3.4 worse) |
| fold, no residual correction (`V52_FB_FOLD=2`) | +139.9% +- 8.4 | 0/16 | -1.26% | P9: within 1 point of P8: **fails** |
| V35, exact residual (`V17_R_RES=1024`) | -3.34% +- 1.18 | 12/16 | +270% | P11: -3% to -8%: holds |
| fold, exact residual | **-11.51% +- 1.21** | 16/16 | +269% | P11: beats both: holds |

- **The fold carries the feedback at no cost.** At the production rank it holds 67% of the thin-leg exact
  feedback's -8.59%, and it removes the rank-2 legs' 1.3% of the bill.
- **At the exact residual it beats the thin-leg exact feedback**: -11.5% against -8.6%. That is the exact
  first-order star plus the exact (2,1) slice.
- **The arm's second duty is large and first order (P9).** Without the correction, M carries 3 d d(e) unopposed and
  the output is 2.4 times worse. P9's reasoning (rank 4 keeps little of a flat term) missed the low-rank part of
  d(e) d^T. D21's additive part u 1^T + 1 v^T (section 5: the dominant part, captured 0.998 at rank 2) makes
  d_add rank <= 4, and the correction removes exactly that part through the residual's top directions.
- **The cost of that.** At R_RES = 4 the correction's low-rank part then competes with the residual's own top four
  directions. That is the likely reason the fold gains 5.7 points more at the exact residual (fold: -5.8 -> -11.5)
  than V35 does (-3.3).
- **GC2 at the theorem weight still loses (P10).** The arm then carries R1 exactly, but the star's overshoot costs
  more. The half weight stays, as in P5 and P6.

**P13 (crowding; stated before its runs).** If the correction's rank-<=4 additive part crowds the residual's own
directions, four more residual columns help the fold more than they help V35. Measured as cold, networks 0-15, paired
against the same batch's V35:
- V35 + `V17_R_RES=8` gains Delta_V;
- fold + `V17_R_RES=8` gains Delta_F over the fold;
- prediction: Delta_F - Delta_V <= -0.5 points;
- R_RES = 16 is run for both as the trend.

Each column costs about 0.46 units (0.26% of the bill: Z transport, MP, PPL, the t1 einsum and the birth range
finder, from the profile of 11e). So R_RES = 8 must buy more than about 1% raw to pay.

**P13, measured** (same harness and pairing; `outputs/fold_crowd_16nets/`):

| variant | raw against V35 | better on | FLOPs |
|---|---|---|---|
| V35 + R_RES 8 | -0.98% +- 0.38 | 12/16 | +0.98% |
| fold + R_RES 8 | **-8.17% +- 1.09** | 16/16 | -0.27% |
| V35 + R_RES 16 | -1.44% +- 0.52 | 12/16 | +2.96% |
| fold + R_RES 16 | -9.24% +- 1.15 | 16/16 | +1.70% |

- **P13 holds.** Four more columns gain 2.4 points on the fold (Delta_F) against 1.0 on V35 (Delta_V), so
  Delta_F - Delta_V = -1.4.
- **The cost matches 11d's estimate.** Each column costs 0.25% of FLOPs.
- **The measured scored-regime effect so far.** The fold on the adopted system at R_RES 4, free-running in the
  scored regime, no counterterms, 64 networks so far: raw -4.71% +- 0.59, C/B -1.37%, adjusted **-6.02% +- 0.58**
  (better on 60/64). The full count is in 11f.

**11e. The additive part through its exact legs (`V52_FB_FOLD=3`; stated before its runs).** P9 and P13 locate the
fold's one cost in the additive part u 1^T + 1 v^T of D21. Through the arm's M-block share it puts a rank-<=4 term into
the birth residual, which then competes with the residual's own directions. Mode 3 removes it at the source:
- the additive part rides the exact rank-2 feedback legs of V49 (row and column sums, all orders, no M share);
- only the flat part D21 - u 1^T - 1 v^T is folded (first order, with the residual correction);
- the flat part's rows at the dropped (saturated) neurons are zeroed as D21's are.

Cost: the fold's -1.25% plus the rank-2 legs' +1.3%, so about V35's bill.
- **P14 (cold, networks 0-15, paired against V35).**
  - Mode 3 at R_RES 4: raw within 1.5 points of fold + R_RES 8 (-8.2%), better on at least 15/16, FLOPs within
    -0.1% to +0.3% of V35.
  - Mode 3 + R_RES 8 gains over mode 3 what R_RES 8 gains on V35 (-1.0), within 0.7 points: no crowding left.
- **P15 (scored, all 100 networks, on the adopted system, no counterterms first).**
  - The better of fold + R_RES 8 and mode 3, chosen on networks 0-49 free-running adjusted, then refitted on 0-49
    and judged free-running on 50-99, beats the pre-registered P12 candidate (the fold at R_RES 4, refitted the same
    way) in held-out adjusted MSE.
  - Expected margin: 1-2.5%, the cold raw advantage less its cost difference.

**P14, measured** (cold, networks 0-15, paired against V35 in the same batch; `outputs/fold_hybrid_16nets/`):

| variant | raw against V35 | better on | FLOPs |
|---|---|---|---|
| mode 3 (`V52_FB_FOLD=3`) | **-11.36% +- 0.71** | 16/16 | -0.01% |
| mode 3 + R_RES 8 | -11.73% +- 0.83 | 16/16 | +0.97% |

- **P14's raw range fails favourably.**
  - Against fold + R_RES 8: -3.38 +- 0.59 points (15/16), more than twice the predicted 1.5.
  - Against the fold at R_RES 4: -5.84 +- 0.61 (16/16).
  - Mode 3 at rank 4 equals the fold at the exact residual (-11.5%), at V35's bill instead of x3.7.
- **No crowding is left.** Four more columns give -0.43 +- 0.23 against V35's -0.98: within the predicted 0.7.
- **The cost prediction holds.**
- **Reading: rank 8 is not enough.** The additive correction (rank <= 4 in exact arithmetic) is not captured
  cleanly by a one-pass range finder at rank 8 competing with the residual's own spectrum. Carrying it on its own
  exact legs is what closes the gap.

**11f. Measured in adjusted MSE** (scored regime, networks 0-99, on the adopted system V35 + pair; `outputs/pe_rows.json`,
`outputs/pf2_rows.json`).

Free-running, no counterterms, paired against the same code's V35 + pair (adjusted 4.1302e-9):

| system | raw | C/B | adjusted | per net | better on |
|---|---|---|---|---|---|
| + fold, R_RES 4 (mode 1) | -4.87% | -1.37% | 3.8752e-9 | -5.96 +- 0.47 | 92/100 |
| + fold, R_RES 8 | -7.42% | -0.30% | 3.8125e-9 | -7.53 +- 0.44 | 95/100 |
| + mode 3 | **-10.62%** | -0.03% | **3.6910e-9** | **-10.53 +- 0.26** | **100/100** |

**Renormalized** (counterterms refitted on networks 0-49 with the stored responses, judged on 50-99;
`outputs/refit_fold1_heldout.txt`, `outputs/refit_fold3_heldout.txt`). The renormalized change against refitted V35:
- adopted pair: -0.88%;
- mode 1: -9.42% +- 0.67 (49/50), free amplitude 0.95;
- fold at R_RES 8: -10.06% +- 0.64, amplitude 1.00;
- mode 3: -10.91% +- 0.41 (50/50), amplitude 1.36.

The counterterms do not absorb the fold: its amplitude stays at its derived value, or above it.

**The decisive test** (free-running, held-out networks 50-99, each system with its own counterterms fitted on 0-49;
`outputs/hold_rows.json`):

| system | raw | C/B | adjusted | against the adopted system |
|---|---|---|---|---|
| adopted: V35 + pair + counterterms (rerun on this code) | 2.1609e-8 | 0.1804 | 3.8977e-9 | (section 9: 3.8984e-9) |
| **+ fold, mode 1** | **1.9726e-8** | **0.1779** | **3.5090e-9** | **-9.97% +- 0.58, better on 50/50** |
| + mode 3 | 1.9484e-8 | 0.1804 | 3.5138e-9 | -9.85% +- 0.35, better on 50/50 |

- **P12 holds, beyond its range.** Held-out adjusted is -10.0% against the adopted system (predicted -3% to -8%), and
  the counterterms cannot mimic the fold. That is twenty standard errors: adopted.
- **P15 fails.** Mode 3 against mode 1, both renormalized: raw -1.09 +- 0.41, but C/B +1.38%, so adjusted
  +0.27 +- 0.42. A tie.
  - Mode 3's free-running lead (4.4 points) is mostly inside the counterterm span. The crowding damage sits in the
    low-rank, coherent part of the (2,1) slice, which per-layer amplitudes can repair.
  - The rank-2 legs mode 3 needs cost 1.38%, which the counterterms cannot repair.
  - Mode 1 is adopted: simpler, cheaper, and no thin legs.

**11g. The cost audit, measured.** Profile of the measured predict (predict #3, `code/prof_chain.py`, networks 0 and
1, `outputs/fold_profile/`), in units of 2n^3:

| component | adopted (section 9) | + fold |
|---|---|---|
| young D21 hub and its thin terms (`_dslices`) | 71.7 | 69.8 |
| young transport (the W [A, P] family) | 53.5 | 53.5 |
| joins (range finders, rotations, Grams) | 23.4 | 23.4 |
| old-tier formings (both tiers) | 20.5 | 20.5 |
| C_pre = W C W^T | 6.6 | 6.6 |
| Wick stage | 1.3 | 1.3 |
| residual legs Z (transport) | 0.7 | 0.7 |
| feedback legs Zf (transport) | 0.36 | 0 |
| kappa4 quenched rank | 0.9 | 0.9 |
| everything else | 1.6 | 1.3 |
| **total** | **180.6** | **178.0** |

By operation: matmul 82%, and the block adds, subtracts and copies of the Strassen recursion 12.7%. Every
elementwise multiply in the chain totals 1.9 units, the QRs 1.8, the saturation gather 0.3.
- **The rank-2 feedback cost 2.4 units (1.3%) for -0.7% of raw** (section 7's fbadd row). The fold that replaces it
  costs about 0.05 units.
- **The two redundancies of 11a cost 1e-4 B and are removed.** The identity check on networks 50-51 reproduces raw to
  the last digit and moves C/B by -0.0001.
- **The Strassen bookkeeping is at its floor under these rules.**
  - The upper levels use classical Strassen: 18 block adds and 4 copies per level.
  - Winograd would need 15 adds but 6 copies for the batched layout, saving one write in 22, about 1 unit.
  - The leaf (16) is the smallest at which a further level pays (block side > 7.5); V35's per-family minima
    already took that.
- **What is left is the tiers' n^3 structure.** The young tier (transport + hub) is 69% of the bill, about 2 units per
  young source-layer. Only a change of tiers moves it, and its price is on the frontier (elasticity about 1). The fold
  is the first lever of this round that moves the frontier itself.

**11h. The leaderboard system (supersedes section 9).** V35, plus the pair (`V33_WK4M=3`, `V31_K4D=3`), plus the fold
(`V52_FB_FOLD=1`), plus its 103 counterterms refitted on all 100 networks. The file is `code/estimator_final.py`, baked
from `est_v29.py` with `code/estimator_final_config.txt`. Section 9's file is kept as `code/estimator_final_pair.py`.

| scored regime, all 100 networks | raw | C/B | adjusted | against section 9's final |
|---|---|---|---|---|
| section 9 final (V35 + pair + counterterms) | 2.1526e-8 | 0.1803 | 3.8815e-9 | |
| **V35 + pair + fold + counterterms (adopted)** | **1.9643e-8** | **0.1778** | **3.4923e-9** | **-10.03% (per net -9.91 +- 0.39, 100/100)** |
| V35 + pair + mode 3 + counterterms (alternative) | 1.9364e-8 | 0.1802 | 3.4901e-9 | -10.08%; against the adopted: +0.08 +- 0.26 |

- **The figures are in-sample for the counterterms.** The protocol's held-out estimate is 11f's 3.5090e-9 against
  3.8977e-9 (-9.97%).
- **Residual and memory.** The measured call's local residual is 0.570 s mean and 0.625 s max, below section 9's 0.61
  and 0.68, because the thin legs are gone. VmHWM is 10.33 GB, unchanged.
- **Against V35 alone** (note XXIX, 4.1273e-9) the system is -15.4% in adjusted MSE.

**11i. What the theory unlocked, and what it says next.**
1. **The exact feedback was never a carrier problem.** Section 10 asked for a carrier that holds the feedback's
   flat spectrum at well under n^2 r per source-layer.
   - The answer is n^2 per birth and nothing after: at first order the feedback is a displacement of the arm, and
     the arm is already carried.
   - Section 7's estimate of one more dense leg per young source (+57% of the bill) is withdrawn.
2. **The arm's double duty is first order and decides everything.**
   - The arm carries the star and the separable part of the M block. Without the residual correction the fold is
     +140%; with it, -5.8% cold and -6% scored.
   - The correction's additive part competes with the rank-4 residual's directions (P13).
   - Removing that competition (mode 3) or renormalizing (the counterterms) closes the rest. The two routes are
     worth the same once cost is counted.
3. **Two full-rank first-order terms of the exact (2,1) slice live only in the birth residual**: R1 = e w2 D21^T / 2
   and R2 = (Phi(1-Phi) - mu w2) Phi D21.
   - The exact residual is worth -3.3% on V35 (P11) and costs +270%.
   - The fold carries part of R1 and R2 for free through the arm. Mode 3 at rank 4 equals the fold at the exact
     residual, which says the remaining residual content is worth little once the arm carries its share.
4. **Next.**
   - The free amplitude 1.36 of mode 3 says the output wants more of the hybrid's direction than first order gives.
     The candidate source is the Gamma x Gamma terms, where the fold and the thin legs differ.
   - Making mode 3 pay needs its rank-2 additive legs confined with their sources at AGE_OLD. Today they are carried
     unconfined to the last layer: 1.38% of the bill, against the 1.1% raw it keeps after renormalization.
