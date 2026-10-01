# Tropical stream — results log

*Updated as runs land. Every number comes from a script in this folder. The seed and width are stated with each result.*

## R-P1: gate temperature (probe_freeze.py; width 256, depth 16, 40k samples, seed 0)

| layer | 1 | 2 | 4 | 8 | 12 | 16 |
|---|---|---|---|---|---|---|
| ρ_l (mean share of 2nd moment) | 0.32 | 0.49 | 0.67 | 0.85 | 0.86 | 0.88 |
| ρ_l, infinite-width correlation map | 0.32 | 0.49 | 0.68 | 0.83 | 0.90 | 0.93 |
| hot (\|μ/s\| < 1) | 1.00 | 0.86 | 0.58 | 0.30 | 0.26 | 0.24 |
| frozen (\|μ/s\| > 3) | 0.00 | 0.00 | 0.01 | 0.12 | 0.27 | 0.27 |

Freezing saturates: a quarter of the gates remain thermal at every depth. **The dominant-cone (zero-temperature) expansion has no small parameter at L = 16.**

## R-P5: Newton-polytope sign cancellation (p5_newton_cancellation.py, exact)

In the tropical form a_l = max(P_l, Q_l) − Q_l, the Gaussian mean widths of the Newton polytopes are exactly additive under Minkowski sums. Their size relative to the answer E a_l (median over neurons):

| layer | 2 | 4 | 8 | 12 | 16 |
|---|---|---|---|---|---|
| n = 256 | 12 | 5.5e3 | 6.2e8 | 9.3e13 | 6.4e18 |
| n = 1024 | 18 | 2.4e4 | 4.1e10 | 7.0e16 | **1.2e23** |

The growth is ≈ E Σ_i |W_ij| = √(4n/π) ≈ 36 per layer at n = 1024. **Signs do not die at zero temperature here; they are the whole answer.** Any estimator built on the P/Q (support-function, mean-width) decomposition must resolve a 10²³-fold cancellation. Dead.

## R-E2: the tropical-curvature split is non-perturbative (e2_wall_split.py)

The identity I2 splits E a_{l,j} into β (own-wall birth, E[δ(z)|∇z|²]) and T (inherited curvature through the gate). Compared with their factorised forms β_G = sφ(μ/s) and T_G = μΦ(μ/s), using the *true* μ and s, the rms over neurons is:

Width 64, depth 8, seed 0, N = 1e5:

| layer | rms m | m − m_G | β − β_G | T − T_G | corr |
|---|---|---|---|---|---|
| 1 | 0.56 | 4.6e-4 | 1.1e-2 (kernel noise) | 1.1e-2 | −1.00 |
| 2 | 0.74 | 6.9e-3 | 0.19 | 0.19 | −1.00 |
| 4 | 0.84 | 9.0e-3 | 0.36 | 0.37 | −1.00 |
| 8 | 0.69 | 1.0e-2 | 0.72 | 0.72 | −1.00 |

The birth and transport each miss by O(1), as large as the answer, and the two misses cancel exactly. The reason is that the tropical multiplicity |∇z|² is uncentred: E|∇z|² ≈ μ² + s² rather than s². The mean spike, which carries ρ ≈ 0.9 of the second moment at depth, sits inside the wall weights.

The equivalent spherical statement: Δ_S F = −(n−1)F + (walls) for F = f|_S. So the identity says that the total wall mass equals (n−1) times the mean, but locally the walls carry the mean and the fluctuation entangled.

**Consequence:** the "first-order tropical correction" (TCT-1) is O(1), not O(n^{-1/2}). The zeroth order equals the Gaussian closure only after replacing the tropical variable |∇z|² by the moment variable s², i.e. after leaving the tropical picture.

Width 64, depth 16, seed 0, N = 2e5. The split worsens with depth: the tropical multiplicity outgrows the variance.

| layer | 1 | 2 | 4 | 8 | 12 | 16 |
|---|---|---|---|---|---|---|
| rms m | 0.56 | 0.74 | 0.84 | 0.69 | 0.39 | 0.38 |
| m − m_G (true μ, s) | 3e-4 | 6.9e-3 | 8.9e-3 | 1.0e-2 | 6.5e-3 | 7.6e-3 |
| β − β_G | 8e-3 | 0.18 | 0.37 | 0.72 | 0.76 | **1.02** |
| median E\|∇z\|²/s² | 1.00 | 1.48 | 2.39 | 4.98 | 10.4 | **16.7** |
| median E\|∇z\|²/(μ²+s²) | 1.00 | 1.31 | 1.27 | 2.65 | 3.51 | 6.40 |

The input-gradient norm of a deep pre-activation is 17× its variance. Even after the mean is added back it is 6× the second moment. The Jacobian picks up "fan roughness": the tropical hypersurface becomes dense and steep with depth, but its walls cancel in the mean. **The tropical variables are the wrong coordinates for this problem at depth.** Width 128: pending.

## R-E0: exactness of I2 at toy scale

- Analytic check: for degree-1 homogeneous a, E[Δa] = E[a(|x|² − n)] (Gaussian shift identity) = E a, using E r³ = (n+1) E r for χ_n.
- `e0_curvature_identity.py` (line integration): agrees within noise at width 8 / depth 3 and width 12 / depth 4. Its RHS weights 1/|v·∇z| have infinite variance, so its z-scores are unreliable.
- `e0b_curvature_kernel.py` (kernel δ, finite variance). The kernel bias is **O(h), not O(h²)**: H jumps wherever another wall crosses. `e0c_bandwidth_scan.py`, width 16, depth 2: the relative bias is 2.4 %, 1.4 %, 0.85 %, 0.39 %, −0.1 % at h = 0.08, 0.04, 0.02, 0.01, 0.005, so it vanishes linearly. With linear extrapolation over h ∈ {0.01, 0.02}:
  - width 16, depth 3, N = 2e6: every layer within noise (max |z| 1.78 / 1.81 / 1.24, 16 neurons per layer). Layer 1 matches the closed form to 9.5e-4.
  - width 12, depth 4: every layer within noise except one neuron at layer 3 (z = 7). At width 12 there are "dead cones", where every input of a neuron is off and z ≡ 0 on an open set. The per-wall decomposition is ambiguous on their boundaries; this has negligible measure at width ≥ 16.

**I2 holds.**

## Stage Q (shared bench, eval_q)

TCT-0 (= Gaussian field + exact one/two-wall cone lift, float64, `tct0.py`, `stageq_tct0.py`):

| set | MLPs | raw final MSE | ± s.e. | all-layer MSE | bench gauss baseline raw |
|---|---|---|---|---|---|
| w64_d16 | 8 | 4.46e-4 | 7.6e-5 | 4.49e-4 | 5.11e-4 |
| w128_d16 | 8 | 2.82e-4 | 6.2e-5 | 1.44e-4 | 2.89e-4 |

- Width fit raw ∝ n^{-0.66}, giving raw(1024) ≈ 7e-5. The bench Gaussian fit (n^{-0.82}) gives 5.2e-5.
- Cost: 19 u, 0.019 B, multiplier 0.1, so **projected adjusted ≈ 5–7e-6** at 1024.
- The exact bivariate (two-wall) lift is slightly better than the bench's linearised cross-covariance, but the difference is inside the s.e.

This is the zeroth order of every tropical realisation (§2 of DESIGN.md).
