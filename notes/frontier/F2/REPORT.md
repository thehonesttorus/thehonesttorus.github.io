# F2: a finite-state cohort for the deep old sources of V29

Child session F2, 2 Oct 2026. Base: `../estimator_v29r3.py` (504aldo V29, MIT; not edited). Everything here lives in the copy `estimator_f2.py`, which keeps the V29 op stream bit for bit when no `F2_*` flag is set.

## Verdict

1. **The cohort works but does not pay at depth 16.**
   - The cohort is one transported basis Q3 (n × k) plus two static k³ Tucker cores. Its readout is algebraically exact.
   - Its joint span really is fixed: at k = 224, six deep sources are held at 99.98 % of their weighted leg energy, and raw matches the base.
   - The content needs k ≈ 200. Each mover's own leg Gram has a 99 % rank of 122 to 155 (ages 9 to 14).
   - Absorbing a source into a k³ core costs about 3·n·k³, and reading the core costs about n·k³ per layer. At k ≈ 200 that is 5 to 10 u per layer. The whole per-source deep tier it replaces costs only 22 u over the run.
   - At k = 64 the cohort is cheap (−12.5 u, −4.8 % C/B). It keeps about 95 % of the deep content's value, but raw gets worse: +26 % on MLP 0 and 5.6× on MLP 1. No (age, k) point beats the base on adjusted score.
2. **Measured value of the deep content** (dropping it entirely, MLPs 0 and 1):
   - ages ≥ 9: raw 7× worse;
   - ages ≥ 12: 2× worse;
   - ages ≥ 14: +8 %.

   The old memory is long. It cannot be truncated by age, and k = 64 is too small for it.
3. **By-product, the post-transport basis (`F2_POSTT=1`).**
   - What it does: the tier-1 basis is chosen by the range finder on W·G·Ω, i.e. on the legs after the transport, so it is optimal for this layer's readout. In V29 it is chosen on the legs before the transport and then multiplied by W.
   - Tier 2's U is then also chosen in the post-transport metric, at no extra cost.
   - Cost: one extra (n,n)×(n,r) product per join, +1.4 % C/B.
   - At rank 384, raw improves on 5 of 6 dev MLPs (−1.6 % mean), so the adjusted score is flat (−0.15 %).
   - The gain grows when the basis is rank-limited: at R_OLD = 320 it is −5 to −16 % raw per MLP versus the same rank without it. It is an enabler for rank cuts, not a win on its own. Paired on MLPs 0 to 2, POSTT with R_OLD2 = 192 is +0.2 % and POSTT with R_OLD = 320 is +3.3 % (see the table).

## 1. Design

Every term of `_dslices` has the same shape per source: a sum over the hub column j of a product of three legs. In D21, two of the legs sit at the row i and one at the row c. In D3, all three sit at the row i. The legs are:

- A and P;
- the M leg M = P·s + 3A·e + Z Lᵀ;
- the feedback legs Xt = F1 R1ᵀ and Yt = F2 R2ᵀ;
- the row-scaling Yk = P y.

If every leg of a cohort source is written as Q3 × (k × n static coordinates), its whole contribution becomes

```
D21_ic = Σ_{s,t,u} Q3_is Q3_it Q3_cu G_stu        D3_i = Σ_{s,t,u} Q3_is Q3_it Q3_iu H_stu
```

G and H are static k³ cores. They are additive over sources and rotate with the basis: G ← G ×₁T ×₂T ×₃T.

**Build per mover.** The mover is grouped by its third leg:

| third leg | pair (Khatri–Rao) sum |
|---|---|
| A | (a, 2w2·p) + (p, e·p + yt) + (xt, w2/3·p) |
| P | (a + xt/3, a·w2 + yt) + (p, s/3·p) + (m, 2/3·p) |
| Z Lᵀ, Yt, Xt | (p, p/3); (a + xt/3, p); (a, w2/3·p) + (p, yt/3). These are low rank and pass through L, R2ᵀ, R1ᵀ (cheap). |
| feed | ρ = (c1·a + c2·p) pᵀ, then G += (2/3) ψ⊗ρ + (1/3) ρ⊗ψ and H += ρ⊗ψ, with ψ = the Y3 column of ζ |
| D3 (third leg P) | (3a + xt, a·w2 + yt) + (m, p) |

**Readout.** Only the part of G that is symmetric in (s,t) contributes, so the core is packed to P = k(k+1)/2 pair rows. Then:

- KQ = Q3[:, iu] ⊙ Q3[:, ju], an (n × P) matrix;
- Y = KQ @ [G_p | H_p];
- D21 += Y₁ Q3ᵀ and D3 += rowsum(Y₂ ⊙ Q3).

**Cohort state and move.**

- The cohort holds slots [0:kc] (ages > `F2_AGE3`). Tier 2 becomes [kc:kb], with a slab offset `fa2_off` and S2 −= the mover's Gram.
- Once per layer, the oldest tier-2 source moves in, using the previous layer's QU and the not-yet-transported Z/Zf stacks:
  1. The new basis comes from a range finder on the weighted leg Gram: mover plus old cohort, 2 passes.
  2. The old coordinates and cores are projected with T = Q3nᵀ(w1 ⊙ Q3).
  3. The mover's coordinates are formed: fa = Q3nᵀ(w1 ⊙ QU) FA2, likewise fp, ζ = Q3nᵀ(w1 ⊙ Z), φ = Q3nᵀ(w1 ⊙ Zf).
  4. Finally Q3 = W Q3n.
- Between moves Q3 ← W diag(w1) Q3.
- The consumers of the old legs are covered: the D21 hub, the D3 rowsums, the feedback thin legs, the feed pair, the PPL Zᵀ thin term and the M leg. The use-side κ4 terms do not read the legs.

Flags: `F2_AGE3` (0 = off), `F2_K3`, `F2_QPASS3`, `F2_CHECK`, `F2_EXACT`, `F2_DIAG`, `F2_DROP`, `F2_POSTT`, `F2_RAISE` (re-raise instead of the fallback).

## 2. Exactness checks (MLP 0)

- **Algebra (`F2_CHECK=1`).** At every layer the cohort's sources are also run through the ordinary per-source `_dslices`, with dense legs formed from the same coordinates. Max relative error at layers 9 to 15, with up to 6 sources: D3 2e-7 to 7e-7 and D21 8e-7 to 1.5e-6, i.e. float32 rounding.
- **No truncation (`F2_EXACT=1`, one mover at age 14, k = 288).** The basis is the exact span [QU | Z | Zf].
  - End to end, raw is 1.9489e-8 against the base's 1.9422e-8 (+0.35 %).
  - The residual comes from the tier-2 basis being re-fitted without the mover (U changes for the remaining tier-2 sources), not from the core.
- **Full cohort, untruncated in practice (ages ≥ 9, k = 224).** Raw 1.934e-8 against the base's 1.942e-8. C/B is not meaningful there: the k = 224 core and the diagnostics cost about 0.6 B.

## 3. Rank and cost diagnostics (`F2_DIAG=1`, MLP 0, ages ≥ 9)

| layer | captured, mover (k = 64) | old cohort (k = 64) | mover (k = 224) | mover's 90 % / 99 % rank |
|---|---|---|---|---|
| 9 | 0.851 | – | 1.0000 | 75 / 155 |
| 10 | 0.876 | 0.976 | 0.9998 | 68 / 144 |
| 11 | 0.898 | 0.986 | 0.9998 | 61 / 135 |
| 12 | 0.913 | 0.993 | 0.9998 | 56 / 129 |
| 13 | 0.923 | 0.995 | 0.9998 | 52 / 126 |
| 14 | 0.930 | 0.997 | 0.9998 | 49 / 122 |

The rank per mover falls with age, but slowly: from 155 to 122 between ages 9 and 14, about age^-0.4 rather than n/age. Six movers together fit in 224 dimensions at 99.98 %. So the span is fixed but not small.

**Cost law.** Units are n³-scale (1 u = 2n³).

| item | flops |
|---|---|
| absorb one mover (3 dense k³ groups) | ≈ 6 n k³ |
| packed readout per layer | ≈ 2 n k³ |
| core rotation per move | ≈ 12 k⁴ |
| per-source tier 2, measured | ≈ 0.65 u per source-layer |

The cohort for ages ≥ 9 can save at most the deep tier's total cost. That total was measured by dropping the tier: 22.4 u, 8.6 % of the bill.

- Absorbing one source pays back only after ≈ 6 n k³ / (1.3 n³) = 4.6 k³/n² of its own source-layers: about 9 at k = 128 and 49 at k = 224.
- At depth 16 the deepest source stays 6 layers.
- In steady state, with one move and one readout per layer (≈ 8 n k³), the cohort breaks even with N deep sources only when N ≈ 6 k³/n². That is about 12 sources at k = 128 and about 66 at k = 224. At depth 16 there are at most 7.
- Measured at k = 64: the cohort costs ≈ 1.4 u per layer (absorption, readout, Q3 transport). That is 9.9 u in total for 22.4 u of deep tier removed.

So a dense or packed Tucker cohort is profitable only for much deeper nets, or if the deep content were representable at k ≲ 64. It is not. A CP-fitted core (option b) is worse still: one ALS sweep on a k³ core costs ≈ 6 k³ R, which is 8 u at k = 224, R = 256.

## 4. Results

The raw scores use the truth noise subtracted. Steady C/B is the mean over the MLPs that are not first in their process (1, 2, 4, 5); MLPs 0 and 3 run staged Strassen L4. Δ is the paired change against the base on the same MLPs.

| config | MLPs | raw per MLP (e-8) | mean raw | steady C/B | adjusted | Δadj vs base (paired) |
|---|---|---|---|---|---|---|
| base V29r3 | 0,1 | 1.942 / 0.384 | 1.1631e-08 | 0.2540 | 2.954e-09 | +0.0% |
| base V29r3 | 0,1,2 | 1.942 / 0.384 / 2.074 | 1.4666e-08 | 0.2533 | 3.714e-09 | +0.0% |
| base V29r3 | 0,1,2,3,4,5 | 1.942 / 0.384 / 2.074 / 1.194 / 2.531 / 2.420 | 1.7576e-08 | 0.2533 | 4.451e-09 | +0.0% |
| cohort a8 k64 (ages>=9) | 0,1 | 2.451 / 2.154 | 2.3025e-08 | 0.2418 | 5.568e-09 | +88.5% |
| cohort a8 k96 | 0,1 | 2.005 / 1.149 | 1.5772e-08 | 0.2610 | 4.116e-09 | +39.4% |
| cohort a8 k128 | 0,1 | 2.133 / 0.732 | 1.4327e-08 | 0.2979 | 4.267e-09 | +44.5% |
| cohort a10 k64 (ages>=11) | 0,1 | 2.167 / 0.763 | 1.4651e-08 | 0.2500 | 3.663e-09 | +24.0% |
| cohort a10 k96 | 0,1 | 2.016 / 0.492 | 1.2542e-08 | 0.2628 | 3.296e-09 | +11.6% |
| cohort a11 k64 (ages>=12) | 0,1 | 1.999 / 0.466 | 1.2325e-08 | 0.2524 | 3.111e-09 | +5.3% |
| drop ages>=9 (zero state) | 0,1 | 13.980 / 11.840 | 1.2910e-07 | 0.2321 | 2.996e-08 | +914.4% |
| drop ages>=12 | 0,1 | 4.038 / 1.601 | 2.8194e-08 | 0.2475 | 6.979e-09 | +136.3% |
| drop ages>=14 | 0,1 | 2.102 / 0.460 | 1.2809e-08 | 0.2527 | 3.237e-09 | +9.6% |
| POSTT | 0,1 | 1.948 / 0.315 | 1.1317e-08 | 0.2576 | 2.916e-09 | -1.3% |
| POSTT | 0,1,2 | 1.948 / 0.315 / 2.062 | 1.4419e-08 | 0.2569 | 3.705e-09 | -0.3% |
| POSTT | 0,1,2,3,4,5 | 1.948 / 0.315 / 2.062 / 1.136 / 2.521 / 2.396 | 1.7298e-08 | 0.2569 | 4.444e-09 | -0.2% |
| base R_OLD=320 | 0,1 | 2.127 / 0.684 | 1.4056e-08 | 0.2376 | 3.339e-09 | +13.1% |
| base R_OLD=320 | 0,1,2 | 2.127 / 0.684 / 2.309 | 1.7067e-08 | 0.2369 | 4.043e-09 | +8.8% |
| POSTT R_OLD=320 | 0,1 | 2.014 / 0.577 | 1.2952e-08 | 0.2406 | 3.117e-09 | +5.5% |
| POSTT R_OLD=320 | 0,1,2 | 2.014 / 0.577 / 2.207 | 1.5993e-08 | 0.2399 | 3.837e-09 | +3.3% |
| base R_OLD2=192 | 0,1 | 1.900 / 0.453 | 1.1762e-08 | 0.2498 | 2.938e-09 | -0.5% |
| base R_OLD2=192 | 0,1,2 | 1.900 / 0.453 / 2.170 | 1.5074e-08 | 0.2491 | 3.755e-09 | +1.1% |
| POSTT R_OLD2=192 | 0,1 | 1.982 / 0.392 | 1.1868e-08 | 0.2534 | 3.008e-09 | +1.8% |
| POSTT R_OLD2=192 | 0,1,2 | 1.982 / 0.392 / 2.043 | 1.4724e-08 | 0.2527 | 3.721e-09 | +0.2% |

The 6-MLP base reproduces the coordinator's number: mean raw 1.7576e-8 at C/B 0.2533, adjusted 4.451e-9.

## 5. What this says about the theory

- **"Rational memory" holds in the weak form.** The deep content occupies a fixed subspace whose dimension does not grow with the number of sources: six sources in 224 dimensions at 99.98 %. A finite-state cohort (Q3, G, H) is therefore exact in principle.
- **It does not hold in the strong form.** The fixed dimension is about 200, not 43 to 116, in the leg metric that the readout needs. At k = 64 the cohort keeps only 85 to 93 % of each mover.
- **The Hadamard wall moves; it does not disappear.** The per-source cost of 8 n² r turns into a per-cohort cost of about 2 n k³ for readout plus 6 n k³ per absorbed source. Both are cubic in the cohort rank. A cohort pays only where k³ ≪ (number of deep sources) · n · r, which needs depth ≫ 16 at these ranks.
- **Basis quality matters more than basis rank (POSTT).** At fixed rank, choosing the basis for the transported legs is worth −1.6 % raw at r = 384 and −5 to −16 % at r = 320. The old tier's error is set largely by the per-layer re-projection (Tq), which V29 does in the pre-transport metric. A multi-layer lookahead (W_{l+1}W_l), possible because all weights are known, is the natural next step. Combined with re-sweeping R_OLD and R_OLD2 under POSTT, it is the coordinator's knob territory.

## Files

- `estimator_f2.py`: V29r3 plus the cohort tier, the POSTT basis and the probes (MIT, `LICENSE-504aldo-MIT`).
- `base_012.json` and `results/*.json|log`: per-MLP results. `results/table.txt` is the table above.
- `sweep*.sh`: the exact commands.
