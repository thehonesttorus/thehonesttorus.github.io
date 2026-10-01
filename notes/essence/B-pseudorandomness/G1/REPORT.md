# G1: pseudorandomness beyond the gap (team B, generalisation round)

*Target from coordinator note 3, §4 G1. Status: final for this round (2 Oct 2026). Files: `toy.py` and
`toy_results.json` for the toy models; `g1_blocks.py` and `g1_blocks_mlp0_c{1,2}.json` for the n = 1024
measurement inside FC (via `../tests/fc_hooked.py`). Labels: THEOREM (proof given or standard, cited), MEASURED,
CONJECTURE.*

## 0. Result in five lines

1. **Where the λ term goes critical.** Jeronimo–Mittal–Roy's (JMR) λ term at walk length Δ is the operator
   W^Δ − Π.
   - On the d-dimensional torus its ε-rank is Θ(N (ln(1/ε)/Δ)^{d/2}). It is **Dixmier-critical, ∝ 1/Δ, exactly at
     d = 2**, the recurrence threshold, where Σ_Δ (W^Δ − Π) is the Green operator, which lies in L^{1,∞} but not in L¹
     (Weyl).
   - Measured: rank·Δ = 12,000 ± 6 % for Δ = 4…512 at N = 4096, against a prediction of 12,000 (Theorem A).
2. **Free products are universally "d = 2".** For a product of a free layers, the squared-singular law's variance is
   additive under ⊠. So the participation ratio is n/(1 + Σ_i v_i) ∝ n/a, at any gating (Theorem B1).
   - Measured: Ginibre PR·(a+1)/n = 0.96–1.00. Our network's r₉₀·a/n = 0.58–0.63 at ages 4–14, flat, and
     r₉₉·a/n = 1.19–1.38.
3. **A frame frozen at a block head keeps its relative error through the block** (in expectation), whenever the
   continuation is left-orthogonally invariant (Theorem B2).
   - Measured: Ginibre leg error constant to 4 digits across a block.
   - In the network at n = 1024 the leg error is constant to ±1 % through each dyadic block, and the D21 readout error
     *falls* within each block, because the signal concentrates in the Perron direction.
4. **The blockwise estimator.** It uses one frame per dyadic block of walk length, of rank r_j = c·n/2^j, so only
   ⌈log₂ L⌉ frames. Its total resolution is within a constant factor of the per-length optimum.
   - Truncation at Δ₀ fails: on the torus its error stays ≈ 1 until Δ₀ ≈ N ln(1/ε); in the network, region N6.
   - Measured at n = 1024 (MLP 0, 4 sources): **c = 2 gives a per-source D21 error of 0.8–2.6 %** across ages 4–14,
     and c = 1 gives 7–17 %.
5. **Consistency with team G's G10.** Compression works along the amenable shift axis (walk length or age), with
   neuron-space frames *refreshed at dyadic times*. No uniform neuron-space grading exists.
   - On the torus the walk group is amenable, so a uniform vertex grading also exists: the Fourier basis with Weyl
     exponent d.
   - On free layers only the blockwise one does.

## 1. The commutative toy: gapless walks on tori (Theorem A)

**Setting.**
- W = ½(I + A/(2d)) is the lazy simple random walk on the torus ℤ_m^d, with N = m^d.
- Its eigenvalues are λ_θ = 1 − (1/d) Σ_i sin²(θ_i/2), for θ ∈ (2π/m) ℤ_m^d.
- The gap is Θ(1/m²) → 0, so the walk is gapless as m → ∞.
- JMR's λ term at length Δ is W^Δ − Π.

**Theorem A.**
1. *Rank law.* For 1 ≪ Δ ≪ m² and fixed ε,

   r_ε(Δ) := #{θ ≠ 0 : λ_θ^Δ > ε} = (1 + o(1)) · V_d (m/2π)^d · (4d ln(1/ε)/Δ)^{d/2}.

   - Proof: λ_θ^Δ > ε ⇔ |θ|² < 4d ln(1/ε)/Δ · (1 + O(|θ|²)), using sin²(x/2) = x²/4 − O(x⁴). Then count lattice
     points in the ball, whose radius in index space is m R/(2π) ≫ 1.
   - So r_ε(Δ) ∝ Δ^{-d/2}. For d = 2 it is (2/π) N ln(1/ε)/Δ: **order exactly one, Dixmier-critical.** For d = 1 it
     decays more slowly (resolution ∝ √T); for d ≥ 3 it is summable.
2. *Green operator.* G = Σ_{Δ ≥ 0} (W^Δ − Π) = (I − W)⁺ has singular values 1/(1 − λ_θ). The number of them above s
   is Θ(N s^{-d/2}), so G ∈ L^{1,∞} \ L¹ iff d = 2. This is Weyl's law, and Connes' statement that Δ^{-d/2} on a
   d-dimensional space is the prototypical Dixmier-class operator.

   **The Dixmier-critical memory is the Green function at the recurrence threshold**, Σ_Δ p_Δ(0, 0) ~ (1/π) ln T.
3. *Blockwise estimator for product functions.*
   - The test is E_RW[∏_{i=1}^k f_i(x_{t_i})] = N⁻¹ ⟨1, M_{f₁} W^{Δ₁} M_{f₂} ⋯ W^{Δ_{k−1}} M_{f_k} 1⟩.
   - For block j = ⌊log₂ Δ⌋, let Q_j be the span of the modes with λ_θ^{2^j} > ε, and replace each W^{Δ_i} by
     Π + Q_j Λ^{Δ_i} Q_j*.
   - Since λ^Δ ≤ λ^{2^j} on the block, every replacement has operator-norm error ≤ ε. Telescoping with contractions
     (‖M_f‖ ≤ 1) gives total error ≤ (k − 1) ε.
   - Frames: ⌈log₂ T⌉ + 1, where T ≈ (m²/π²) d ln(1/ε) suffices, so O(log N).
   - Resolution used: Σ_Δ r_ε(2^{⌊log₂Δ⌋}) ≤ 2^{d/2} Σ_Δ r_ε(Δ), within a constant of optimal.
   - **Truncation** (W^Δ → Π for Δ > Δ₀) has error ≥ λ_min-gap^{Δ₀} = (1 − Θ(1/m²))^{Δ₀}, which is ≈ 1 unless
     Δ₀ = Ω(m² ln(1/ε)) = Ω(N ln(1/ε)) at d = 2. **Truncation needs all ≈ N lengths; blocks need O(log N).**

**Measured** (`toy.py`, ε = 0.01; N = 4096 for each d).

| Δ | 4 | 8 | 16 | 32 | 64 | 128 | 256 | 512 |
|---|---|---|---|---|---|---|---|---|
| d = 1: rank·√Δ | 5076 | 5329 | 5456 | 5521 | 5552 | 5566 | 5568 | 5567 |
| d = 2: rank·Δ | 12376 | 12736 | 12032 | 12288 | 11776 | 12288 | 11264 | 10240 |
| d = 3: rank·Δ^{3/2} | 26800 | 37700 | 29800 | 30800 | 28700 | 26100 | 24600 | – |

- The predicted constants are 5,596 for d = 1 and 12,000 for d = 2. At d = 2, Δ = 512 is near the finite-size edge
  (Δ ~ m²/8).

## 2. The noncommutative toy: free products (Theorems B1, B2)

**Theorem B1 (free Weyl law; standard free probability).**
- Let X_a = G₁ ⋯ G_a, with the G_i asymptotically free and each |G_i|² of unit mean and variance v_i in the n → ∞
  limit. For example, G_i = D_i W_i with W_i Gaussian and D_i a gate diagonal.
- For free a, b of unit mean, φ(abab) = φ(a²) + φ(b²) − 1. So **the variance is additive under ⊠**:
  Var(|X_a|²) = Σ_i v_i.
- Hence PR(X_a) = n/(1 + Σ_i v_i) = Θ(n/a), and Σ_{a ≤ L} PR(X_a) ≈ (n/v̄) ln L: **Dixmier-critical along the
  product length, with no geometry.**
- Ginibre: v = 1 and PR = n/(a+1), which is the Fuss–Catalan second moment.
- For energy ranks at fixed δ, the large-a triangular law (Newman's Lyapunov spectrum; Tucci's a-th-root limit)
  suggests r_{1−δ} ≈ n ln(1/δ)/(a+1). This is a CONJECTURE at the level of constants: the a-th-root limit does not
  control O(1/a) deviations.

**Theorem B2 (frozen frames survive invariant continuations).**
- Let X ∈ ℝ^{p×n}, let Q be any n × r orthonormal frame, and let B be random with B ~ OB for every O ∈ O(n). Then

  E ‖X(I − QQᵀ)B‖²_F / E ‖XB‖²_F = ‖X(I − QQᵀ)‖²_F / ‖X‖²_F.

- Proof: E[BBᵀ] commutes with every O, so it equals βI.
- For a gated continuation B = D B′ with B′ invariant, the ratio is ‖X(I − QQᵀ)D‖² / ‖XD‖². Choosing Q from XD
  (gate-weighted) restores the identity.
- Concentration around the expectation is Hanson–Wright, with relative fluctuation O(n^{-1/2}) per layer.
- **So the cost of a block is set at its head**: one SVD there, then transport of an r × n frame at n² r per step.

**Measured** (`toy.py`, n = 1024).

| quantity | a = 1 | 2 | 4 | 8 | 12 | 16 |
|---|---|---|---|---|---|---|
| Ginibre PR·(a+1)/n | 1.00 | 1.00 | 0.99 | 1.00 | 0.96 | 0.99 |
| Ginibre r₉₀·a/n | 0.51 | 0.75 | 0.97 | 1.14 | 1.21 | 1.23 |
| gated (φ ~ U(0, 1)) r₉₀·a/n | 0.40 | 0.52 | 0.61 | 0.69 | 0.72 | 0.72 |

Frozen frame on a Ginibre product: relative leg error at each step of the block.
- head a₀ = 2, r = 512: 0.1186, 0.1186;
- head a₀ = 4, r = 256: 0.248–0.250;
- head a₀ = 8, r = 128: 0.326–0.328.

## 3. The measurement on our network (n = 1024, bench w1024_d16, MLP 0, inside FC)

`g1_blocks.py` follows the legs L_s(t) = [Y; Z; T] of the sources born at s = 1, 2, 3, 5.

**(i) The rank law.** For source 1 at age a:

| a | 1 | 2 | 4 | 6 | 8 | 10 | 12 | 14 |
|---|---|---|---|---|---|---|---|---|
| r₉₀ | 428 | 258 | 148 | 104 | 79 | 64 | 54 | 46 |
| r₉₀·a/n | 0.42 | 0.50 | 0.58 | 0.61 | 0.62 | 0.62 | 0.63 | 0.63 |
| r₉₉·a/n | 0.72 | 0.99 | 1.19 | 1.27 | 1.32 | 1.35 | 1.38 | 1.38 |

- Sources 2, 3 and 5 agree within 5 %.
- **Dixmier-critical (∝ 1/a) from age 4 on**, matching the gated free model (0.69–0.72). Team C's PR ≈ n/(2(a+1))
  is the participation-ratio form of the same law.

**(ii) The blockwise frozen frame.**
- Blocks are ages {1}, {2, 3}, {4…7}, {8…15}. At each head the frame is refrozen from the exact legs, with rank
  r = c·n/a₀. Afterwards only the coefficients and the r × n frame are carried.
- Errors are relative Frobenius, for the source's own D21 readout and for its legs.

| c | source | ages 2–3 | ages 4–7 | ages 8–14 |
|---|---|---|---|---|
| 1 | s = 1 | D21 11.3 → 10.8 %, leg 9.7–9.8 % | D21 15.9 → 12.2 %, leg 14.2–14.4 % | D21 13.4 → 8.3 %, leg 16.6–16.8 % |
| 1 | s = 5 | 7.6 → 7.2 % | 13.6 → 13.3 % | 17.0 → 16.5 % |
| 2 | s = 1 | exact (r = n) | D21 1.7 → 1.3 %, leg 1.6 % | D21 2.4 → 1.5 %, leg 3.1 % |
| 2 | s = 2, 3, 5 | exact | 0.8–1.4 % | 1.5–2.6 % |

**Readings.**
- Theorem B2 holds in the quenched network: the leg error is constant to the third digit inside every block.
- The D21 readout error falls inside each block, because transport concentrates the signal in the Perron (dilation)
  direction, which the head frame always contains.
- **c = 2 is close to lossless** at ≤ 2.6 % per source. That agrees with team D's 2n/a law and with region §9
  (Tucker n/4 lossless for ages 5–8).
- The per-source transport cost is ⌈log₂ L⌉ blocks × 2^j layers × n² · 2n/2^j ≈ 2n³ log₂ L, against n³ L exact.
  At L = 16 that is 8 n³ against 16 n³ per source.
  - The saving grows as L/(2 log₂ L). It is modest at L = 16.
  - Merging the sources of one block into a shared head frame (the odometer of coordinator note 2 §3) is what
    removes the factor of L. Region §9 measured that version for one bin.

**Caveats.**
- Refreezing at each head here starts from the exact legs. A causal carrier refreezes the *approximate* legs, so the
  head errors add across the ≈ log₂ L blocks. If they are independent (fresh-weight lemma) the total is
  √(Σ_j ε_j²) ≈ 3–4 % at c = 2. Not run end to end.
- The errors are per-source readout errors; the raw score with every source blocked was not run.
- One network.

## 4. Against team G's G10

- G10 says that with ≥ 2 independent fresh layers (a quantum expander), no grading of neuron space with bounded
  commutators is finitely summable.
- Both toys agree with it:
  - **Torus (amenable ℤ^d).** A uniform grading exists, the Fourier/Laplacian one, finitely summable with exponent
    d. The frames are nested and never refreshed (Q_{j+1} ⊂ Q_j).
  - **Free layers.** No uniform grading exists. The Dixmier law appears only along the **amenable shift axis**
    (product length, i.e. age), and the neuron-space frames are valid **only blockwise**: refreshed at dyadic heads,
    each one inherited unchanged for one block by Theorem B2.
- **The clean statement (CONJECTURE, supported by A, B1, B2 and the measurement).** On a gapless layered computation,
  blockwise compression is possible exactly along the amenable (shift) direction, with block ranks given by the
  ⊠-variance law. Along expanding (free) directions no time-uniform compression exists.
- The ε-rank exponent along the shift axis is d/2 for commutative diffusions and 1 for free products.
- **He-ReLU therefore sits at "d_eff = 2"** for the free reason (additive ⊠-variance), not the geometric one.

## 5. Honest status

| kind | items |
|---|---|
| THEOREM | A (elementary: lattice-point counting, plus Weyl for the Green operator); B1 (additivity of the free variance under ⊠, standard); B2 (one line) |
| MEASURED | the toys; the network rank law and the blockwise frozen frames, on one network and 4 sources, at the per-source readout level |
| NOT DONE | the end-to-end raw score of a causal blockwise carrier; the energy-rank constant of B1 (conjecture) |
| CONJECTURE | the amenable-direction compression principle (§4) |

**The main limitation.** The product-function estimator in Theorem A is for commuting (normal) walks. For free
products, B2 gives the expectation statement only. A full error bound for the k-point readout needs B2 plus
concentration, uniformly over the blocks.
