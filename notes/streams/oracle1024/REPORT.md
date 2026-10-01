# oracle1024: the K = 3 closure ladder at the real shape (width 1024, depth 16)

Status: **done** (2026-10-01). Width 1024, depth 16, 2 MLPs (seeds 770000, 770001) × 2 independent replicas × N = 3.5e6;
widths 128 and 256 on the same two seeds for the width trend. Stream closed by the coordinator after this run.

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

All ε below are relative rms errors of D21(l+1) in %, estimated noise-free by the replica cross-product
ε² = ⟨M_A − D_A, M_B − D_B⟩ / ⟨D_A, D_B⟩. The jackknife standard error over 32 output-column blocks is
≤ 0.02 points at width 1024, layers ≥ 1. Per-layer tables with errors:
[`results/summary_width1024.md`](results/summary_width1024.md) (production) and
[`results/summary_prelim.md`](results/summary_prelim.md) (widths 128 and 256, and the N = 32k pilot). Raw
printouts are `results/width*_pair.txt`; the full rows, including within-sample and oracle-style `_xc` columns,
are in `results/width*_pair.json`.

### Verdict

At the real shape (n = 1024, L = 16), on both MLPs and at every layer 0–14:

1. **The leg-partition closure with the exact (2,1,1) κ4 slice is below the 2.2 % bar everywhere**: 0.74–1.07 %,
   with the maximum at layer 1. Its error falls about as n^(−0.9) with width (deep layers 128: 3.7–7.1 %,
   256: 2.0–2.9 %, 1024: 0.74–0.91 %).
2. **The fitted closure is also below the bar everywhere**: 0.66–1.07 %. At 1024 it gains only 0.0–0.15 points over
   the derived coefficients.
3. **With the (2,1,1) slice regenerated as u_i C_jk (the published chain's r = 1 core), both closures miss the bar**
   at almost every layer: 2.14–2.69 %. A few layers of MLP 770001 (1, 10, 12) are 2.14–2.15 %, which is marginal;
   the fitted version is 2.10–2.64 %. The regeneration closes only about 20 % of the gap between "no κ4 term"
   (2.4–3.2 %) and "exact κ4" (0.7–1.1 %). At width 128 it closes about 55 %, and at 256 about 33 %.
4. **The transported old content T(Φ³κ3(z)) is about 40 % of D21(l+1) at every width and depth, and it is not
   absorbable.** Dropping it costs 37–42 %. Replacing it by its held-out least-squares projection on everything born
   at the layer (slices, Wick, Gaussian ρ², ρ³, B1–B6) still leaves 35–40 %; the held-out R² is ≤ 0.15 at 1024,
   0.02–0.46 at 256 and 0.04–0.69 at 128. The share does not shrink with width, and absorbability falls with width.

### Ladder by width and depth (both MLPs pooled; ranges over layers and MLPs)

| width | layers | D21 noise, 1 replica | closure δ-noise | slices only | leading Wick | + Gaussian ρ³, ρ⁴ | closure without κ4 | **leg-partition closure** | closure, κ4 = u_i C_jk | **fitted closure** | fitted, κ4 = u_i C_jk | closure without old content | old content → projection on births | old-content share | held-out R² of the projection |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 128 | 1–3 | 0.021–0.044 | 0.007–0.017 | 44.4–50.9 | 7.5–9.6 | 7.5–9.7 | 6.4–8.7 | 2.5–3.5 | 5.5–7.5 | 2.2–3.4 | 5.4–7.3 | 36–41 | 31–40 | 0.36–0.41 | 0.04–0.27 |
| 128 | 4–6 | 0.015–0.022 | 0.005–0.007 | 47.7–54.1 | 9.0–11.0 | 9.0–11.0 | 8.0–10.8 | 2.9–4.4 | 6.6–8.7 | 2.6–3.5 | 6.2–8.1 | 40–47 | 34–39 | 0.40–0.46 | 0.25–0.32 |
| 128 | 7–14 | 0.007–0.018 | 0.002–0.005 | 51.6–61.1 | 7.0–10.3 | 6.7–10.4 | 8.0–12.2 | 3.7–7.1 | 5.9–8.3 | 1.7–3.0 | 3.4–7.6 | 44–62 | 30–42 | 0.44–0.59 | 0.34–0.69 |
| 256 | 1–3 | 0.031–0.058 | 0.010–0.024 | 45.8–53.2 | 5.5–7.7 | 5.5–7.7 | 4.7–6.5 | 1.8–2.4 | 4.1–5.7 | 1.7–2.3 | 4.1–5.7 | 37–43 | 35–41 | 0.37–0.43 | 0.02–0.16 |
| 256 | 4–6 | 0.019–0.028 | 0.005–0.009 | 43.9–50.0 | 6.2–7.8 | 6.2–7.8 | 5.5–6.5 | 1.7–2.4 | 4.2–5.7 | 1.5–2.2 | 3.9–5.6 | 37–43 | 33–38 | 0.37–0.43 | 0.14–0.25 |
| 256 | 7–14 | 0.010–0.018 | 0.002–0.005 | 44.7–55.4 | 6.6–8.4 | 6.5–8.4 | 5.7–7.4 | 2.0–2.9 | 4.3–5.9 | 1.6–2.0 | 3.8–5.3 | 40–54 | 35–44 | 0.40–0.53 | 0.24–0.46 |
| **1024** | 1–3 | 0.066–0.119 | 0.022–0.049 | 46.9–49.2 | 2.9–3.6 | 2.9–3.6 | 2.4–2.9 | **0.87–1.07** | 2.15–2.51 | **0.86–1.07** | 2.14–2.50 | 38–40 | 37–40 | 0.38–0.40 | 0.00–0.03 |
| **1024** | 4–6 | 0.038–0.055 | 0.010–0.017 | 44.4–47.8 | 3.6–4.0 | 3.6–4.0 | 2.8–3.1 | **0.80–0.92** | 2.40–2.65 | **0.78–0.90** | 2.39–2.64 | 38–39 | 36–38 | 0.38–0.39 | 0.05–0.08 |
| **1024** | 7–14 | 0.018–0.035 | 0.004–0.009 | 40.5–46.0 | 3.8–4.4 | 3.8–4.4 | 2.7–3.2 | **0.74–0.91** | 2.14–2.69 | **0.66–0.83** | 2.10–2.64 | 37–42 | 35–39 | 0.37–0.42 | 0.08–0.15 |

Layer 0 (Gaussian input) is omitted. There, every model except slices-only (28 %) is exact to within noise
(|ε| ≤ 0.5 %). Width 128 and 256 use N = 5e5 and 1e6 per replica. The width-128 streaming numbers reproduce the
dense oracle_k3 tables of competition-plan §3.1 on different MLPs: closure 3–7 %, fit 1.4–3 %, the same depth drift.

Approximate scaling of the deep-layer ε with width, from midpoints at 128, 256 and 1024:
- leg-partition closure: n^(−0.9);
- closure without κ4: n^(−0.6);
- regenerated-κ4 closure and fitted closure: n^(−0.5);
- leading Wick: n^(−0.35);
- slices only: constant at 45–60 %.

### Fitted coefficients at width 1024 (D21-space fit, both MLPs, layers 1–14)

| basis term | leg-partition value | width 1024 | width 128 (same seeds) |
|---|---|---|---|
| B1 D21(z) hyperedge + edge (w2, w2) | 3 | 2.98–3.05 | 2.2–3.1 |
| B2 D21(z) hyperedge + edge (w3) | 3 | 2.89–3.01 | 2.45–3.09 |
| B3 D3(z) hyperedge + two edges (w5) | 1 | 0.15–0.40 | −0.2–0.43 |
| B4 Gaussian ρ³ | 1 | 0.82–1.00 | 0.28–0.94 |
| B5 K22 hyperedge + edge | 1.5 | 1.23–1.55 | 0.25–1.44 |
| B6 (2,1,1) κ4 hyperedge | 1.5 | 1.36–1.49 | 1.05–1.45 |

The depth drift that §3.1 read as renormalisation by second-order diagrams is mostly a finite-width effect. At 1024
the coefficients stay within 0.15 of the leg-partition values at all depths, except B3, which stays at 0.2–0.4
at every width (an open discrepancy in a small term). The mild residual drift is B5 1.55 → 1.26 and
B6 1.49 → 1.36 from layer 1 to 14. The coefficients agree between the two MLPs to ±0.05 layer by layer.

### Noise and N

At N = 3.5e6 one replica's D21(l+1) has relative noise 0.12 at layer 1, 0.05 at layer 4 and 0.018–0.035 at
layers 7–14. This scales as 1/√N: at the N = 32k pilot it was 0.88, 0.49 and 0.2–0.4. The 1 % target for D21
itself would need N ≈ 1.1e7 (layer 14) to 4e7 (layer 7) per replica. That is 30–100 h on this 4-core box
(≈ 290 samples/s per single-thread process, 13 GEMM-equivalents per sample and layer), so it was not reached.
It is not needed for the ladder. The model error δ = M − D shares most of its sampling noise with the target,
so the replica noise of δ for the closure is 0.4–0.9 % at layers ≥ 6 (1.0–4.9 % at layers 1–5). The
cross-product estimate removes what remains: the N = 32k pilot already gave closure ε = 0.81–0.98 % at layers 6–14,
matching the production 0.74–0.91 %.

## Open issues

- The exact (2,1,1) κ4 column is an oracle. A chain cannot form κ4(z)_iijk, and its transport here is a sample
  contraction; in a chain it is an O(n⁴) object. What this stream establishes is the target: a carrier of the κ4
  slice's *transport* must close ≥ 70 % of the 2.4–3.2-point no-κ4 gap to reach 2.2 %, and ≥ 90 % to reach 1 %.
  The r = 1 regeneration closes about 20 %.
- The old content (≈ 40 % of D21, held-out R² ≤ 0.15 on the births at 1024) has to be carried by any chain; no
  memoryless closure absorbs it.
- B3's fitted coefficient (0.2–0.4 against a counted 1) is unexplained; the term is small.
- There are only two MLPs, and both are seed-regenerated He-init (the same distribution as the competition's).
- The jackknife over output-column blocks assumes the noise is roughly independent across blocks. The agreement
  between the N = 32k pilot and the production run (≈ 100× more samples) is the empirical check.
- The Hermite columns' Φ³κ3(z) term uses the empirical gate (oracle_k3 uses Φ(α)). This changes ε by ≤ 0.002 at
  width 32, and the Gaussian ρ³, ρ⁴ terms are invisible at 1024 anyway (≤ 0.02 points).

## Files

- [`stream_oracle.py`](stream_oracle.py): streaming statistics (`stats`), dense verification (`verify`), and
  single or pair ladder (`ladder`, with `--gauss-cache` to share the O(n⁴) Gaussian terms between replicas).
- [`summarize.py`](summarize.py): per-layer markdown tables from the pair JSON files.
- [`run_widths.sh`](run_widths.sh) and [`run_prod_ladders.sh`](run_prod_ladders.sh): the small-width and
  production drivers. Production stats: `stream_oracle.py stats --seed {770000,770001} --width 1024
  --n-samples 3500000 --chunk 2048 --sample-seed {1001,1002 | 2001,2002}`. The stats files are 0.5 GB each and
  are not committed.
- `results/verify_width32.txt` and `results/verify_width128_layers0-3.txt`: term-by-term agreement with oracle_k3.
