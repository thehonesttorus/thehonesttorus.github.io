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

**Decision (adjusted MSE).** T is adopted only if the measured free-running V35 + T + cal_T beats V35 + cal_35 on the
held-out networks in adjusted MSE, by more than two standard errors, cost included.
