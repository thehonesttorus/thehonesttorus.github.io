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

(widths 64 and 128 at depth 16: pending)

## R-E0: exactness of I2 at toy scale (e0_curvature_identity.py)

Width 8, depth 3, M = 2,000 lines: LHS and RHS agree within MC error at every layer (max z-score 2.3 over 24 neurons). Layer 1 matches the closed form |w|/√(2π).

(M = 60,000 and width 12 / depth 4: pending)

## Stage Q

(pending: bakes running)
