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
