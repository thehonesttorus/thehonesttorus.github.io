# oracle1024: the K = 3 closure ladder at the real shape (width 1024, depth 16)

Status: **in progress** — method verified; preliminary width-1024 ladder at N = 32k (one MLP); production runs (2 MLPs × 2 replicas × N = 3.5e6) running.

## Question

notes/competition-plan.md §3.1 measured, at widths 32 and 128, how well compressed models of the post-activation
third cumulant κ3(a_l) reproduce the next layer's (2,1) slice D21(l+1)_ab = κ3(z_a, z_a, z_b), the whole interface
of a K = 3 chain to the next ReLU (error law: extra MSE ≈ 4.2e-6 ε²; frontier needs ε ≤ 2.2 %). Does the ladder
(slices only → leading Wick → Gaussian ρ³/ρ⁴ → leg-partition closure → fitted closure; (2,1,1) κ4 slice exact vs
regenerated u_i C_jk) hold at n = 1024, and do the leg-partition and fitted closures reach ε ≤ 2.2 %?

## Method: O(n³) per layer, no n³ tensor

Code: [`stream_oracle.py`](stream_oracle.py).

1. **Every model term is a "triangle form"** X_ijk = f_i g_j h_k A_ij B_jk D_ik (vectors: gates Φ, Gaussian weights
   w2, w3, w5, Hermite coefficients c_p, the regenerated u; matrices: C, C_off, ρ^e, D21(z), K22) or a sum of
   permutations of one (sym3 = 6 relabelled triangle forms). Transport to D21(l+1):
   T(X)_ab = Σ W_ia W_ja W_kb X_ijk. When one of A, B, D is absent (every closure diagram B1, B2, B3, B5, the
   leading Wick terms, the Gaussian ρ² terms and the regenerated κ4 term) it factorises into two n×n GEMM
   products and one Hadamard, O(n³), e.g. A absent: T = [(WᵀfD) ∘ (WᵀgB)] (h∘W). The only true triangles are
   the Gaussian ρ³/ρ⁴ Hermite terms (p, q, r ≥ 2): O(n⁴), computed as one fp32 GEMM per output column a
   (≈ 10 s per term per layer at n = 1024).
2. **All-distinct projection, exactly**: T(all_distinct X) = T(X) − [S_A + S_B + S_C − 2 S_diag], where
   S_A = (W∘W)ᵀ X_iik W, S_B = [(W ∘ (X_kii W))ᵀ W], S_C = [(W ∘ (X_iki W))ᵀ W], S_diag = (W∘W)ᵀ diag(X_iii) W —
   the three two-index restrictions of any triangle form are n×n matrices read off its factors (inclusion–exclusion
   over p = q, q = r, p = r; the pairwise intersections are all p = q = r).
3. **The two n³ objects enter only through sample contractions**, with y = z_l − μ_l (exact sample mean; pass 1 and
   pass 2 regenerate the same samples, so every central moment equals the dense oracle's):
   - Φ³κ3(z) (B0, also inside the Wick column): T = E[x_a² x_b], x = (Φ∘y) W_{l+1}; restrictions from D21(z).
   - (2,1,1) κ4 slice (B6): F_ijk = E[y_i² y_j y_k] − var_i C_jk − 2 C_ij C_ik. sym3 with the doubled index in an
     a-slot or the b-slot gives T(B6) = ⅔ T_ad(Y0) + ⅓ T_ad(Y2) with T(Y0) = E[s_a x_a x_b] − T(Gauss0) and
     T(Y2) = E[x_a² s_b] − T(Gauss2), s = (w2∘y²) W_{l+1}; Gaussian parts are triangle forms; restrictions are
     F_iik = E[y_i³ y_k] − 3 var_i C_ik, F_kii = E[y_k² y_i²] − var var − 2C², F_iii.
   - regenerated slice u_i C_off,jk: u_i ‖C_off‖² = E[y_i² q] − var_i ⟨C, C_off⟩ − 2 (C C_off C)_ii − 2 Σ_k F_iik C_off,ik,
     q = yᵀ C_off y per sample.
   - slices of κ3(a): D3(a), D21(a) = E[ã_i² ã_j] (n×n); target D21(l+1) = E[y_a² y_b] at layer l+1.
   Per sample and layer: 13 GEMM-equivalents of m×n×n (2 in pass 1, 11 in pass 2).
4. **Fitted closure in D21 space**: least squares of D21(l+1) − T(slices) − T(Wick) on the seven transported basis
   matrices T(B0..B6) (and on B0..B5 + T(B6reg) for the regenerated variant).
5. **Noise**: two replicas (independent sample seeds) per MLP. Reported: (i) the oracle_k3 `--pair` convention
   (model from A against target of B, target noise subtracted) and (ii) the replica cross-product
   ε² = ⟨M_A − D_A, M_B − D_B⟩ / ⟨D_A, D_B⟩, which removes model-input noise (transported Φ³κ3(z) and κ4
   contractions are as noisy as the target) as well as target noise.

Convention difference to oracle_k3 (one column only): the Hermite columns' Φ³κ3(z) term uses the empirical gate,
oracle_k3 uses Φ(α) there; the effect is ≤ 0.002 in ε at width 32 (`verify` prints both).

## Verification against the dense oracle

Width 32, depth 5, seed 770000, N = 1e5, same samples (`moment_atlas_np.py --k3 --k4 --sample-seed 11`):
every transported term (slices, Wick, B0–B6, H4/H6/H8, B6reg, u) agrees with the dense einsum to 1e-7 – 1e-5
relative (fp32 accumulation); the ladder columns agree to 4 digits (memless 0.2460/0.2460, wick 0.0259/0.0259,
closure 0.0147/0.0147, closure_reg 0.0225/0.0225 at layer 0; same at layers 1–3). The D21-space fit is ≤ 0.002
better than the oracle's tensor-space fit, as it must be (it minimises the reported quantity).

## Results

### Noise at width 1024, and what it means for the estimate

One replica's D21(l+1) is far noisier at width 1024 than at 128: at N = 32 768 its relative noise is
0.88 (layer 1), 0.49 (5), 0.28 (10), 0.20 (15), i.e. ε_noise ≈ 0.2 √(32768/N) at depth. Extrapolating 1/√N,
ε_noise ≤ 1 % needs N ≈ 1.4e7 per replica at layer 14 and ≈ 4.5e7 at layer 7 (≈ 2.5e8 at layer 1): 15–45 hours of
this 4-core box per MLP pair, out of reach here. Two facts make it unnecessary for the *model error*:
(i) model and target are built from the same samples (x = (Φ∘y)W is the linearisation of z_{l+1}), so the replica
noise of δ = M − D is 2.5–5× below that of D (closure column: 0.04–0.2 at N = 32k vs 0.2–0.9);
(ii) the replica cross-product ⟨δ_A, δ_B⟩/⟨D_A, D_B⟩ is noise-free in expectation, and its jackknife error is small
at depth: at N = 32k the closure's ε² is resolved to ±0.05e-4 (ε ≈ 0.9 ± 0.03 %) at layers 6–14, while layers 0–2
are unresolved. The production N = 3.5e6 cuts all noise by 10.3×.

### Preliminary ladder, width 1024, MLP seed 770000, N = 32 768 × 2 replicas

ε of D21(l+1) in %, replica cross-product (± jackknife s.e.); full table in
[`results/summary_prelim.md`](results/summary_prelim.md), raw output in `results/width1024_mlp770000_pair_N32k.txt`.

| layers | slices only | leading Wick | + Gaussian ρ³ | closure without κ4 | **leg-partition closure** | closure, κ4 regenerated u_i C_jk | **fitted closure** | fitted, regenerated κ4 |
|---|---|---|---|---|---|---|---|---|
| 3–5 | 46–47 | 3.3–4.1 | 3.4–4.1 | 2.4–3.2 | 0.0–0.9 (± 0.3) | 1.9–2.9 | 0–0.8 | 1.9–2.9 |
| 6–9 | 45–46 | 4.0–4.2 | 4.0–4.2 | 2.9–3.1 | 0.8–1.0 | 2.5–2.7 | 0.8–0.9 | 2.5–2.7 |
| 10–14 | 44–46 | 4.3–4.5 | 4.3–4.5 | 3.0–3.2 | **0.84–0.95** | 2.5–2.7 | **0.7–0.85** | 2.5–2.7 |

D21-space fitted coefficients (B0..B6, theory 0, 3, 3, 1, 1, 1.5, 1.5): at every layer 1–14
B1 2.90–3.13, B2 2.89–2.99, B3 0.17–0.46, B4 0.8–1.0, B5 1.25–1.50, B6 1.40–1.49 — **no depth drift**.

### Width trend (same two seed MLPs, streaming pair ladders, layers 7–14; N = 5e5 at 128, 1e6 at 256)

| width | leading Wick | closure without κ4 | leg-partition closure | closure, regenerated κ4 | fitted closure | fit coef B6 at layer 14 |
|---|---|---|---|---|---|---|
| 128 | 7.0–10.3 % | 9.4–12.2 % | 3.7–7.1 % | 5.9–8.3 % | 1.7–3.0 % | 1.05–1.15 |
| 256 | 6.6–8.4 % | 5.7–7.4 % | 2.0–2.9 % | 4.3–5.9 % | 1.55–2.04 % | 1.21–1.27 |
| 1024 (prelim., 1 MLP) | 4.0–4.5 % | 2.85–3.17 % | 0.81–0.98 % | 2.46–2.70 % | 0.72–0.89 % | 1.42 |

(the width-128 streaming numbers reproduce the dense oracle tables of competition-plan §3.1 on other MLPs:
closure 3–7 %, fit 1.4–2.7 %, coefficient drift with depth.)

Readings so far. (1) The leg-partition closure's error falls roughly as n^(−0.8) (deep layers 4.7 → 2.6 → 0.9 %),
much faster than the leading Wick term (≈ n^(−0.25)); at width 1024 the *un-fitted* closure is already below the
2.2 % frontier bar at layers ≥ 3. (2) The coefficient drift seen at width 128 (renormalisation by second-order
diagrams) is a finite-width effect: at 1024 the fitted coefficients equal the leg-partition ones within noise at
every depth, and the fit buys only ~0.1 point over the derived closure. (3) The exact (2,1,1) κ4 slice is
essential: without it 3 %, with the r = 1 regeneration u_i C_jk 2.5–2.7 % (above the bar), exact 0.9 %. The
regeneration's relative capture of the κ4 term gets *worse* with width: of the gap between "no κ4" and "exact κ4" it closes ≈ 55 % at width 128, ≈ 33 % at 256 and ≈ 20 % at 1024.
