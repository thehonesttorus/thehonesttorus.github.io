# H: rational memory and response-preserving relabelling (coordinator experiments, 1 Oct 2026, 22:00–23:00 UTC)

**Question (from the user's input of 22:00 UTC).** Is the old third-order memory "rational" in the sense of the Connes–Duchamp–Reutenauer construction: does a small space of continuation responses determine its future? And can its histories be relabelled so that exactly the questions the future asks are preserved, rather than the histories themselves (the Kraus-equivalence principle)?

**Model.** `rmem.py`. The Wick part of the first-order memory (coordinator note 1 §2) on the exact Gaussian closure (`region/gclose.py`). A birth at layer s, neuron r has legs Y_s(t) = S_s D_s Z_s(t) and Z_s(t) = W_{s+1} D_{s+1} ⋯ W_t, and contributes w2 [(y∘y)zᵀ + 2(y∘z)yᵀ] to D21(z_t). All old atoms at a cut t₀ share the future propagator Q_{t₀→t}. The readout that matters is the covariance correction E∘D21 with E_ab ≈ (p_a/s_a)P_b, priced by region's K_off(t). The cut is t₀ = 8; "old" means age ≥ 3 at the cut (6 sources, 6n atoms); future targets are t = 8 … 14. Gradient checked to 9 digits.

## Results (width 256, MLP 0 unless stated)

**1. The linear continuation (Hankel) space is large.** Eigen-spectrum of the priced response Gram of the old atoms:

| width | atoms | rank for 90 % / 99 % of response energy | actual memory trajectory: rank for 90 % / 99 % |
|---|---|---|---|
| 256 | 1,536 | 379 / 912 | 173 / 690 |
| 512 | 3,072 | 893 / 2,008 | 251 / 1,391 |

**2. But at each future target the readout is nearly low rank.** Rank of the old memory's priced readout matrix for 90 % / 99 % of its energy, targets 8 → 14:

| width | 90 % | 99 % |
|---|---|---|
| 256 | 10, 8, 7, 6, 5, 5, 4 | 34, 28, 24, 22, 19, 18, 16 |
| 512 | 14, 10, 7, 5, 4, 3, 2 | 61, 49, 39, 32, 29, 25, 23 |

**3. Response-matched merged histories reproduce the whole future.** Fit R new rank-one histories at the cut, transported exactly by the true future propagators, to the old memory's priced future readouts. Residual energy relative to the old memory's priced readout (Adam, 800 iterations), against keeping the R strongest original histories:

| R | n/32 | n/16 | n/8 | n/4 |
|---|---|---|---|---|
| response-matched merge | 14.1 % | 7.7 % | 3.4 % | **1.1 %** |
| pruning to the R strongest | 91 % | 84 % | 84 % | 79 % |

The histories go from 6n = 1,536 to n/4 = 64, a 24× reduction. For comparison, region's tensor-norm CP merge at R = n was 6.7× worse in raw MSE (region REPORT §8). The merged histories lie outside the linear span of the old ones (result 1): the relabelling is nonlinear, new "Kraus operators" rather than combinations of the old ones.

**4. A short lookahead mostly suffices** (R = n/4). Fit only targets t₀ … t₀ + w, evaluate on all future targets:

| window w | residual on the fit window | residual over the whole future | residual beyond the window |
|---|---|---|---|
| 0 (the cut only) | 0.2 % | 13 % | 21 % |
| 1 | 0.7 % | 4.7 % | 12.5 % |
| 2 | 0.9 % | 2.3 % | 7.6 % |
| whole future | 1.1 % | 1.1 % | — |

**5. A cheaper fitter.** For fixed first legs the readout is linear in the third legs, so they come from one least-squares solve (conjugate gradient on the normal equations). Starting from the R strongest histories' first legs (R = n/4), with work counted in passes (one pass transports and reads out R histories at every fitted target):

| fitter | passes | residual over the whole future (whole-future fit / window-2 fit) |
|---|---|---|
| one least-squares solve for the third legs | 21 | 6.2 % / 6.6 % |
| 3 alternating sweeps (LS for third legs, 30 Adam steps on first legs) | 174 | 1.4 % / 2.5 % |
| Adam, 800 iterations | 800 | 1.1 % / 2.3 % |

**6. What the residual means for the score.** Old content (age ≥ 3) holds 42–56 % of the priced D21 energy at targets 8–14 (young 22–35 %, positive cross terms 20–27 %). So 1.1–1.4 % residual on the old part is ≈ 0.6–0.7 % of the total, an amplitude error ≈ 8 %. Region's tolerance is 5–10 % uniform for the whole error budget (N2). Team B's constant-in-a channel (80 % of old energy, carried at ≈ 0 cost) could be removed first, leaving only the traceless remainder to merge.

## Reading

- **The user's principle is quantitatively right here.** Preserving the continuation questions (the priced future readouts) instead of the history labels compresses 6n histories to n/4 with ≈ 1 % residual energy, where history-preserving compressions (pruning; CP in the tensor norm) fail.
- **It is not linear rationality.** The linear continuation space (Hankel rank) is ≈ 1.5–1.7n and grows with width. What is small is the per-target readout (rank 16–61 at 99 %) and the nonlinear (rank-one history) fit to it.
- **Cost, honestly.** Carrying R = n/4 merged histories costs ≈ 7·R·n² per layer ≈ 1.75 u at n = 1024, against ≈ 6n histories otherwise. The open problem is the fit: at n = 1024, R = 256 and a 3-target window, one pass is ≈ 2 u, so 21–174 passes per merge is 40–350 u. Warm-started incremental merges (fold the one newly aged source into the existing merged set each layer, a few CG steps) are the obvious next test.

## Still running when written

Width 512 fits (full window and window 2), and MLPs 1–2 at width 256. Results are appended below when they finish.

## Appended 22:50 UTC: width 1024, and the causal incremental merge

**7. Width 1024 (MLP 0, cut 8, ages ≥ 3: 6,144 histories), alternating fitter, 2 sweeps (113 passes).** Per-target readout rank for 99 % of the energy, targets 8 → 14: 116, 94, 78, 69, 58, 50, 43 (90 %: 21 … 2), so it grows roughly in proportion to width at the cut and more slowly later. Residual over the whole future:

| merged histories | window: cut + 2 layers | whole-future fit |
|---|---|---|
| n/8 = 128 | 6.3 % | 4.2 % |
| n/4 = 256 | 3.8 % | 1.9 % |

So the residual at fixed R/n is similar to width 256 (slightly higher with this shorter fit), and the effect holds at the competition width.

**8. The causal incremental merge (`incr.py`).** Young sources (ages 1 … AY) are exact. Each layer, the source that turns AY + 1 is folded into a merged set of R histories by an alternating fit warm-started from the previous merged set, using only the next three targets (window 2). Evaluated at every target against the exact old content. Width 256, MLP 0; the price-weighted total error is relative to the total priced D21 energy (young + old):

| exact young ages | R | sweeps | total priced error (energy) | amplitude |
|---|---|---|---|---|
| 1–2 | n/4 | 2 | 1.84 % | 13.5 % |
| 1–2 | n/4 | 4 | 1.58 % | 12.6 % |
| 1–2 | n/2 | 2 | 1.02 % | 10.1 % |
| 1–3 | n/4 | 2 | 0.88 % | 9.4 % |
| 1–3 | n/2 | 2 | 0.47 % | 6.9 % |
| 1–4 | n/4 | 2 | 0.43 % | 6.6 % |
| 1–4 | n/2 | 3 | **0.18 %** | **4.3 %** |

- **Errors do not accumulate.** With ages 1–2 exact and R = n/4, the residual on the old content stays at 3.6–6.7 % of its energy at every target from 3 to 14; it does not grow with the number of merges.
- **Inside the tolerance.** Ages 1–4 exact plus n/2 merged histories puts the total priced error at 4.3 % in amplitude, inside region's 5–10 % tolerance (N2).
- **The walls that remain.** (i) The exact young tier (ages 1–4) is the dominant cost (≈ 300 u dense at n = 1024 by team D's accounting); the merged old tier itself costs ≈ 7·R·n² per layer ≈ 3.5 u at R = n/2. (ii) The fit is still ≈ 100 passes per merge (≈ 1,000 passes per network), far too many; a frame-restricted fit that compares readouts only on their low-rank part (rank 43–116 at 99 % per target at n = 1024) should cut a pass by ≈ 50–100×, and the warm start should cut the number of passes. Untested.
- **Running:** the same causal merge at width 1024 (ages 1–4 exact, R = n/2, 3 sweeps), `results/incr_w1024_mlp0_AY4_R0.5.log`.

**9. A cheaper fitter (swap and solve).** No gradient steps: at each merge, replace the weakest fraction of the merged histories by the new source's strongest histories, then re-solve all third legs by one conjugate-gradient least-squares solve. Width 256, MLP 0, window 2:

| exact young ages | R | swap | CG steps | passes per network | total priced error (energy) | amplitude |
|---|---|---|---|---|---|---|
| 1–4 | n/2 | 1/8 | 15 | 232 | 0.79 % | 8.9 % |
| 1–4 | n/2 | 1/4 | 15 | 232 | 0.70 % | 8.4 % |
| 1–4 | n/2 | 1/4 | 30 | 412 | 0.45 % | 6.7 % |
| 1–3 | n/2 | 1/4 | 15 | 248 | 1.31 % | 11.4 % |
| 1–4 | n/2 | ALS, 3 sweeps (for comparison) | | 1,240 | 0.18 % | 4.3 % |

About 5× fewer passes buys 2.5–4× more error energy: the accuracy-for-work curve is the object to improve, not one setting.

## Appended 23:00 UTC: the causal merge at width 1024, and what it can buy

**10. Width 1024, causal merge (MLP 0; ages 1–4 exact, R = n/2 = 512, window 2, 3 alternating sweeps; 1,240 passes, 478 s).** Price-weighted total error **0.69 %** of the D21 energy (**8.3 %** amplitude). At width 256 the same settings gave 0.18 % (4.3 %).

| target t | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|---|---|
| old-content residual, n = 1024 | 1.65 % | 1.26 % | 0.99 % | 0.99 % | 1.71 % | 2.34 % | 3.07 % | 3.06 % | 3.71 % | 4.26 % |
| old-content residual, n = 256 | 0.55 % | 0.63 % | 0.54 % | 0.68 % | 0.77 % | 1.04 % | 1.16 % | 1.17 % | 1.25 % | 1.14 % |
| old share of D21, n = 1024 | 0.02 | 0.06 | 0.10 | 0.17 | 0.23 | 0.30 | 0.35 | 0.38 | 0.39 | 0.42 |

- At width 1024 the residual on the old content grows slowly with the number of merges (≈ 1 % → 4.3 %). At width 256 it is flat. The old share is also larger at 1024 (0.42 against 0.31 at t = 14).
- Two readings, not yet separated:
  1. the fit is less converged at width 1024 with the same iteration budget (15 CG steps on a 512 × 1024 least squares);
  2. the rank needed grows faster than n.
  The whole-future ALS test (§7) gave similar residuals at fixed R/n at the two widths, which favours (1).
- **Edge of tolerance.** The result sits at the edge of region's 5–10 % amplitude tolerance, not comfortably inside it.

**11. What the merged-history carrier can buy in the competition (cost arithmetic).**

Cost model. Total D21 cost ≈ (number of targets, 15) × (resolution per target, in units of n atoms) × (products per unit of resolution) × (Strassen factor, ≈ 0.66 for FC's mix).

- FC's constant is 7 products per unit (6 at age 1). It is 4 without the slice legs.
- For the score bar's price (≤ 0.1 B ≈ 100 u), resolution × constant must stay ≤ ≈ 10 per target.

Comparison with the age-graded frames at equal resolution, width 1024:

| carrier | resolution per deep target | measured loss |
|---|---|---|
| this one: ages 1–4 exact + R = n/2 | 4.5 n | 8.3 % amplitude on priced D21 (score level not measured) |
| team G's uniform law k(a) = c·n/a, ages 1–2 exact, c = 1.5 | ≈ 4.7 n | +16–35 % on FC raw |
| same, c = 2 | ≈ 5.6 n | +1–4 % on FC raw |
| team D's frames, k = 2n/a | ≈ 5.6 n | lossless |

The merge packs everything older than four layers into n/2 atoms, where c = 2 spends ≈ 2.5 n on ages ≥ 5. But it needs a fit that the frames do not, and its score-level loss is unmeasured. It is a different point on the same curve, not a new curve.

**The binding cost is the young pairs, not the memory.**
- Ages 1–2 exact are 29 (source, target) pairs. At FC's constant (6–7 products) that is ≈ 190 u dense, ≈ 125 u wall-feasible.
- Adding the covariance arrow (30 u dense, ≈ 20 u wall-feasible) gives ≈ 145 u ≈ 0.14 B, *before any older content*, even if memory compression were free.
- At FC's raw (3.0e-8) that floors the adjusted score at ≈ 4.2e-9: below the public chain's 5.4e-9, but ≈ 2.6× the leaders' 1.6e-9.
- For comparison, the public chain's young pairs cost 2.17 u each wall-feasible (≈ 4 dense products), against FC's ≈ 4.3.

**Verdict on the rational-memory direction.**
- **What it establishes.**
  - The old memory has a small continuation space at every target: per-target readout rank 43–116 at 99 % at n = 1024.
  - That space can be carried causally by response-matched merged histories without error accumulation at width 256, and with slow growth at width 1024.
- **What it does not establish.** The competition wall has moved: to the per-pair constant of the exact young tier, and to FC's raw accuracy (2–3× the leaders' raw). The merged-history carrier attacks neither.
- **What would move the score.** A representation of ages 1–2 that does not materialise the Hadamard squares (Y∘Y, Y∘Z, Z∘Z, Z∘T) per source per target, or a cut of the per-pair constant (7 → 4 by dropping or sketching the slice legs where their score-level weight is small). Both are measurable inside FC on the six networks.

## Appended 23:00 UTC: can the per-pair constant be cut by keeping slice legs only at young ages?

**12. Score-level ablation: slice legs only up to age K** (`fc_age.py` = region's `fc.py` plus `slice_age`, and a per-pair product counter; `run_sliceage.py`; width 1024; raw with truth noise subtracted).

- **Counter.** Of FC's 825 D21 products, the Wick legs take 465 and the slice legs 360. So slice legs are 44 % of the D21 bill.
- **Sanity check.** Full FC reproduces region's MLP 0 score (3.24e-8).

| slice legs kept for ages ≤ K | slice products | MLP 0 | MLP 1 | MLP 2 | worse than full FC by |
|---|---|---|---|---|---|
| all (full FC) | 360 | 3.24e-8 | 1.81e-8 | 3.03e-8 | — |
| 4 | 162 | 5.04e-8 | 6.29e-8 | 6.93e-8 | 1.6× / 3.5× / 2.3× |
| 2 | 87 | 1.15e-7 | 1.37e-7 | 1.17e-7 | 3.6× / 7.6× / 3.9× |
| 1 | 45 | 1.68e-7 | 2.19e-7 | 1.74e-7 | 5.2× / 12× / 5.7× |

**Negative: the slice memory is as long as the Wick memory.**
- Slice legs older than four layers still carry a factor of 1.6–3.5 in raw MSE. The constant of 7 products per pair cannot be cut by age-truncating the slice legs.
- Row pruning of slice atoms also fails (region: 50 % pruning → 2.5× worse).
- What remains is to treat slice legs like Wick legs, in age-graded frames: team G's c·n/a law already projects both.
