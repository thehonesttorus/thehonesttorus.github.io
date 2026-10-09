# F3. Response-directed disintegration along collective modes

Frame F3 of the localization round. Question: decompose the input X = U Y + X_perp along a k-dimensional subspace U
chosen by theory, compute leaf integrals E[F | Y = y] with a chain on N(Uy, I - UU^T) (or a cheap response update of
one chain), integrate over the leaf space by deterministic quadrature, and say honestly whether any regime beats raw
1.55e-8 at <= 0.2 B or <= 0.1 B.

**Verdict: negative for Phase 2, with a clean reason and one decisive experiment that could overturn it.**

- The theory is exact and short. A k-dimensional linear localization changes any Gaussian-law estimator by the
  path integral of its heat defect over the U-block (Theorem 1). In the second-chaos sector, conditioning on U opens
  every closed walk that visits U once into a path of the leaf, and moves walks visiting U twice or more into the
  transverse quadrature (Theorem 3).
- **The collective subspace is real and stays small at width 1024.** The chain's own first-chaos map J_L
  (the mean-gated Jacobian product) has participation ratio 26 at n = 1024, L = 16. Its top direction carries 14% of
  the second spectral moment, 46% of the fourth and 74% of the sixth (section 3). Higher walk moments, the gain's
  quadratic form among them, concentrate on a handful of outlier directions.
- **Localization along those directions removes a large share of a closure's error at small width.** With exact
  Gauss-Hermite leaves along the top-k response directions, k = 1:
  - 0.58-0.78 of the Gaussian closure's final MSE stays, at n = 32-256 and L = 8-16;
  - 0.58-0.76 of the tree-level (kappa3/kappa4) closure's MSE stays, settling near 0.72 from n = 128 on;
  - random directions leave 0.97-1.00.
- **It never pays, because a leaf is a whole chain.**
  - Leaves cannot share the n^3 work: a leaf model that does not recompute the gates reproduces the root exactly
    (Lemma 7).
  - Any rule that sees the localization's effect needs N >= k + 1 leaves (Theorem 8).
  - On the cost frontier of elasticity 1 that the chain sits on, an N-leaf localization beats the base in adjusted MSE
    only if its MSE ratio is rho < 1/N, at any budget split (Theorem 9).
  - Measured with the minimal rules, rho N is 1.05-1.55 for k = 1, 1.2-2.7 for k = 2 and 1.5-6.1 for k = 4, at every
    width tested and always above 1. At n = 1024 the production chain's error is the incoherent bulk, so rho is
    expected to be 0.70-1.00, against the 0.5 it would need.
- **By-products of independent interest:**
  - the exact bias trickle-down (Duhamel) identity, a probe that locates a chain's bias in input directions;
  - the homogeneity form of the isotropic defect: the midpoint estimator (G + Lap_m G)/2 cuts the Gaussian closure's
    MSE 4-50x at n = 64-128, with correlation 0.83-0.99 between the per-neuron bias and the defect trace. It is not
    affordable at n = 1024 as it stands, so it is handed to frame F1;
  - the walk-opening theorem, relevant to anyone attacking the n^4 closed-walk (triangle) wall.

Code: `code/` (all checks below are reproducible from it); printed outputs: `out_*.txt`.

---

## 1. Objects

### 1.1 Estimators on Gaussian laws

An estimator is a functional G(m, Sigma) of the input law N(m, Sigma) with values in R^n (one entry per output neuron).
The chain, the Gaussian closure (GC), the gain-aware closure (GAC) and tree-level propagation (TLP) are all of this
form. Assumptions, all satisfied by GC, GAC and TLP. The production chain is expected to satisfy them, since all its
fitted constants are dimensionless amplitudes; (A3) is checked in step 0 of section 7:
- (A1) G is C^2 in m and C^1 in Sigma on {Sigma > 0};
- (A2) G is exact on point masses, G(x, 0) = F(x), and continuous as Sigma -> 0. With zero covariance every
  cumulant and closure term vanishes and the chain is the forward pass;
- (A3) G is positively 1-homogeneous, G(lambda m, lambda^2 Sigma) = lambda G(m, Sigma).

The exact functional T(m, Sigma) = E F(m + Sigma^{1/2} Z) satisfies (A1)-(A3) and the heat equation
d_Sigma T = (1/2) Hess_m T. Define the **heat defect** of G:

    D_G(m, Sigma) = d_Sigma G - (1/2) Hess_m G        (a symmetric n x n matrix per output neuron;  D_T = 0).

### 1.2 Linear localization

For U (n x k, orthonormal), P = U U^T and strength c in [0, 1], the leaves are N(U a, I - c P) with
a ~ N(0, c I_k). This is the exact posterior of X given the observation sqrt(tau) U^T X + noise, with
c = tau / (1 + tau); c = 1 is full conditioning on Y = U^T X. The localized estimator is

    G_U^c = E_a G(U a, I - c P).

For T it is constant in c (disintegration).

### 1.3 Leaf space, transverse measure and the NCG reading

The equivalence relation R_U = {(x, x') : U^T x = U^T x'} is a proper smooth groupoid, and it is a fibration.
- Its C*-algebra is C_0(R^k) (x) K(L^2(R^{n-k})), Morita equivalent to C_0(R^k): type I, with a manifold leaf space.
- The transverse measure is gamma_k, the law of Y.
- The "central part of the localization" (stage-8 Theorem A) is the conditional expectation onto the centre
  L^infinity(R^k), that is, E[ . | Y].

F3 is therefore the **type-I corner** of the noncommutative-disintegration picture. That is exactly why it is
computable: a manifold leaf space supports ordinary Gauss-Hermite quadrature. The same fact bounds it, in two ways:
- a smooth k-dimensional transversal can only resolve what is a function of k coordinates;
- the chain's residual at width 1024 is a bulk object, which note XII found carried by no finite pattern transversal
  (the "pathological quotient").

The NCG content here is honest but thin. The noncommutative structure (non-type-I leaf spaces, isotropy, modular
cocycles) enters only when the leaves are defined by the network's gates (frame F2) and not by linear coordinates.

### 1.4 Gaussian-leaf rigidity, and why the radial mode is useless

**Lemma 1.** If N(0, I) = integral of N(0, Sigma_w) P(dw), then Sigma_w = I almost surely.

*Proof.* For every xi, the fourth cumulant of xi . X under the mixture is 3 Var_w(xi^T Sigma_w xi) >= 0. For N(0, I)
it is 0, so xi^T Sigma_w xi is almost surely constant, equal to |xi|^2. Taking a countable dense set of xi gives
Sigma_w = I. QED.

So the quadratic collective modes cannot be disintegrated with Gaussian leaves: the radius, and the gain
G_l ~ |x_l| / |x|, which is to leading order the quadratic form x^T J_l^T J_l x. Only mean-shift (linear) localization
keeps the chain's input class.

The radial mode separately gives nothing. By (A3), conditioning on |X| = r gives r times one spherical problem, so
every homogeneous estimator already "knows" it (Monte Carlo measured x0.993, frontier-tests). The gain's linear
shadow, the top response directions, is what F3 can reach (section 3).

---

## 2. Exact identities

### Theorem 1 (Duhamel / bias trickle-down along a k-dimensional localization). Proved; checked.

Under (A1),

    d/dc G_U^c = - E_a tr( D_G(U a, I - c P) P ),       a ~ N(0, c I_k).

Hence:
- (a) the localized bias is b_U = b - integral_0^1 E_a tr(D_G(Ua, I - cP) P) dc;
- (b) for isotropic localization (P = I) under (A2), the root bias is exactly

      b = G(0, I) - E F = integral_0^1 E_{a ~ N(0, cI)} tr D_G(a, (1-c) I) dc
        = E integral_0^inf tr( D_G(m_t, Sigma_t) Sigma_t^2 ) dt        (Eldan's clock: c = t/(1+t), Sigma_t = I/(1+t)).

*Proof.* Write a = sqrt(c) xi with xi ~ N(0, I_k). Then
d/dc G(U sqrt(c) xi, I - cP) = grad_m G . U xi / (2 sqrt c) - tr(d_Sigma G P).
By Stein's lemma in xi, E[grad_m G . U xi] = sqrt(c) E tr(U^T Hess_m G U), so the first term averages to
(1/2) E tr(Hess_m G P). For (b), integrate from 0 to 1, use G(x, 0) = F(x) (A2) at c = 1, and change the variable.

Integrability near c = 1 follows from (A3): D_G(a, sI) = s^{-1/2} D_G(a / sqrt s, I), an integrable singularity when
D_G(., I) is bounded. QED.

This is the bias analogue of the covariance trickle-down V_t - E V_{t+1} = C^T R C of the brief. In both, a
non-martingale (the posterior covariance; the estimator evaluated on the posterior) drifts by a quadratic form in
the localization direction. For the bias, the form is the heat defect instead of the identity.

*Check* (`code/f3_defect.py`, `out_defect_n64L8s1.txt`): d/dc (b . G_U^c) at c = 0, numerically -1.0048e-2, against
-tr(D_b P_U) = -1.0050e-2.

### Theorem 2 (homogeneity form of the isotropic defect). Proved; checked.

Under (A3):
- m . grad_m G + 2 tr(Sigma d_Sigma G) = G;
- at the root, tr D_G(0, I) = (G - Lap_m G) / 2;
- D_G(0, sI) = s^{-1/2} D_G(0, I).

So the frozen-root evaluation of Theorem 1(b) gives b ~ 2 tr D_G(0, I) = G - Lap_m G. The frozen-corrected
estimator is Lap_m G(0, I), the chain-smoothed form of the exact identity E F = E Lap F (Euler plus Stein).

*Euler check:* 2 tr(d_Sigma g) = 0.9531059 against g = 0.9531059.

**Proposition 2b (the midpoint is exact for scale errors).** If G(m, Sigma) = T(m, (1 + eps) Sigma), then to first
order in eps:
- b = (eps/2) T(0, I);
- tr D_G(0, I) = eps tr d_Sigma T(0, I) = (eps/2) T;

so b = tr D_G(0, I), and the estimator G - tr D_G = (G + Lap_m G)/2 is exact. *Proof:* D_T = 0 and Euler.

The path factor beta = |b|^2 / (b . tr D) is therefore 1 for a pure scale error, against the frozen value 2. The
Gaussian closure's error is gain-dominated, a scale-mixture error, which predicts beta ~ 1.

*Measured* (`code/f3_euler.py`, `out_euler.txt`). MSE ratio to the estimator's own root MSE. "mid" is
(G + Lap_m G)/2 and "lap" is Lap_m G.

| net | GC: mid / lap | GC: beta, corr(b, trD) | GAC: mid / lap | GAC: beta, corr(b, trD) |
|---|---|---|---|---|
| n64 L8 s0 | 0.052 / 1.24 | 0.99, 0.971 | 0.232 / 2.21 | 0.93, 0.923 |
| n64 L8 s1 | 0.085 / 1.85 | 0.89, 0.971 | 0.285 / 2.93 | 0.83, 0.930 |
| n64 L8 s2 | 0.252 / 3.60 | 0.72, 0.967 | 0.296 / 2.26 | 0.98, 0.875 |
| n64 L8 s3 | 0.105 / 1.71 | 0.93, 0.930 | 0.703 / 3.55 | 1.07, 0.744 |
| n64 L16 s0 | 0.161 / 0.86 | 1.24, 0.917 | 0.266 / 2.40 | 0.92, 0.905 |
| n64 L16 s1 | 0.021 / 1.05 | 1.01, 0.986 | 0.072 / 1.16 | 1.03, 0.954 |
| n64 L16 s2 | 0.266 / 1.78 | 1.08, 0.830 | 0.223 / 0.94 | 1.31, 0.892 |
| n128 L8 s0 | 0.134 / 1.44 | 1.03, 0.928 | 0.631 / 5.06 | 0.72, 0.910 |

**Reading.**
- The per-neuron heat-defect trace predicts the per-neuron bias of both closures at correlation 0.74-0.99. The path
  factor beta sits at 0.7-1.3, and the frozen value 2 is wrong.
- The midpoint removes 73-98% of the Gaussian closure's MSE and 30-93% of GAC's. (GAC's NaN entries in section 4 are its
  own numerics: a negative variance on a few leaves.)
- This is isotropic (k = n) localization at first order. It belongs to frame F1 and is passed on there. Its cost at
  n = 1024 is the trace of the chain's input-mean Hessian:
  - by finite differences, 2n chain evaluations;
  - by Hutchinson probes, about 4e6 probes for 1e-4 per-neuron accuracy (per-probe std
    sqrt(2) |Hess_m G_i|_F ~ 0.2);
  - by a forward-Laplacian pass, the tangent of C_l with respect to the input mean, an n^3 object per layer.
- None is affordable as it stands. The scalar b . Lap_m G is cheap by Hutchinson when Hess(b . G) has a high
  participation ratio, and that is the F1 probe proposed in section 7.

### Theorem 3 (walk opening under localization). Proved; checked.

In the second-chaos (generalized chi-square) sector the chain carries, z = L . x + (x^T H x - tr H)/2. Write the
blocks of H in the decomposition U (+) U-perp as H_uu, H_up, H_pp, and L_p = P_perp L. Conditionally on y = U^T x, z
is a generalized chi-square in x_perp with:
- first-chaos vector L_p + H_pu y;
- second chaos H_pp.

Hence

    E_y kappa3(z | y) = 3 L_p^T H_pp L_p + 3 tr(H_up H_pp H_pu) + tr H_pp^3,
    E_y kappa4(z | y) = 12 (L_p^T H_pp^2 L_p + tr(H_up H_pp^2 H_pu)) + 3 tr H_pp^4,

and kappa3(z) = E_y kappa3(z|y) + 3 Cov(E[z|y], Var(z|y)) + kappa3(E[z|y]) is exact.

Split the closed 3-walks by their number of U-visits:

    tr H^3 = tr H_uu^3 + 3 tr(H_uu H_up H_pu) + 3 tr(H_up H_pp H_pu) + tr H_pp^3,

so that:
- walks visiting U **once** become **open walks (paths)** of the leaf, the class the chain computes cheaply (D3, the
  Schur hub);
- walks visiting U **twice or more** move into the **transverse quadrature** (exact);
- walks **avoiding** U stay closed: the n^4 wall shrinks to tr H_pp^r.

This is the brief's "conditioning is not freezing" identity in operational form. The leaf's second chaos is the frozen
compression H_pp, but the excursion term tr(H_up H_pu) = tr(q H (1-q) H q) with q = P_U reappears as the leaf's
first-chaos energy, and the leaf chain must carry it.

*Check* (`code/f3_walk.py`, `out_walk.txt`; n = 12, k = 3):
- E_y kappa3(z|y): formula 2.03944, Monte Carlo 2.03853;
- E_y kappa4(z|y): 38.0162 against 37.9898;
- the total-cumulance decomposition: 1.8610098661 = kappa3, to 10 digits;
- the visit classes sum to tr H^3 exactly.

*Consequence.* A k-dimensional localization removes the closed-walk class in proportion to the U-content of the
walk moments. Section 3 shows that the higher walk moments are concentrated on very few directions.

*Measured on the network's exact second-chaos kernels.* H_i = E[F_i (X X^T - I)] (Stein, Monte Carlo 4e6 samples,
n = 32, L = 8; `code/f3_hwalk.py`, `out_hwalk.txt`). The tables give the share of sum_i tr H_i^r on walks **avoiding**
U.

| r | top response directions, k = 1 / 2 / 4 / 8 | random directions, k = 1 / 2 / 4 / 8 |
|---|---|---|
| 2 | 0.84 / 0.70 / 0.45 / 0.23 | 0.92 / 0.89 / 0.80 / 0.58 |
| 3 | 0.77 / 0.61 / 0.31 / 0.12 | 0.87 / 0.85 / 0.76 / 0.47 |
| 4 | 0.71 / 0.54 / 0.21 / 0.07 | 0.83 / 0.79 / 0.70 / 0.37 |

Split of tr H^3 by U-visits (inside U, two visits, one visit opened into leaf paths, avoiding):

| | inside U | two visits | one visit | avoiding |
|---|---|---|---|---|
| k = 1 | 0.07 | 0.10 | 0.06 | 0.77 |
| k = 4 | 0.27 | 0.25 | 0.18 | 0.31 |

Higher closed walks concentrate faster on the response directions, as the spectral moments of section 3 predict.
At n = 64, L = 8 (8e6 samples; the split halves agree to 0.1% on sum_i tr H_i^3), the share of walks avoiding U is:

| r | top response directions, k = 1 / 2 / 4 / 8 | random directions, k = 1 / 2 / 4 / 8 |
|---|---|---|
| 3 | 0.79 / 0.57 / 0.44 / 0.26 | 0.96 / 0.89 / 0.80 / 0.70 |
| 4 | 0.70 / 0.45 / 0.33 / 0.17 | 0.94 / 0.85 / 0.72 / 0.61 |

The one- and two-direction shares hold from n = 32 to n = 64. The tail spreads with width, as the shallow births'
isotropic legs gain weight.

### Proposition 4 (variance capacity). Proved.

- Cov(E[F | Y]) >= J P_U J^T (Loewner order), with J = E grad F. This holds because Cov(F, U^T X) = E[grad F] U
  exactly (Stein), and the regression residual is orthogonal.
- E Var(F_i | Y) <= E |P_perp grad F_i|^2, by the Gaussian Poincare inequality on each fibre (the bound of Zahm et al.,
  gradient-based dimension reduction, SIAM J. Sci. Comput. 2020; the active-subspace literature).

These bound the **variance** a localization removes. The NNGP chaos weights at L = 16 put 19% of sigma^2 in the first
chaos (`code/f3_nngp.py`: 0.189, 0.153, 0.103, 0.080, ...; tail past order 8 is 0.30). Even the top-64 response
subspace, holding 96% of the first chaos, conditions away only about 18% of sigma^2 plus its higher-chaos share. The
leaves are far from localized. The chain is not wrong about variance, though. The relevant capacity is the bias
capacity.

### Proposition 5 (first-order bias capacity, Ky Fan). Proved (first order); checked.

Let D_b = sum_i b_i D_i at the root. To first order in the strength,

    d MSE / dc |_0 = -(2/n) tr(D_b P_U)   and   tr(D_b P_U) <= sum_{j <= k} lambda_j^+(D_b).

The optimal first-order k-subspace is the top-k eigenspace of D_b.

*Measured* (n = 64, L = 8, s = 1, `out_defect_n64L8s1.txt`):

| | k = 1 | k = 2 | k = 4 | k = 8 |
|---|---|---|---|---|
| share of tr D_b in D_b's top-k eigenvectors | 0.36 | 0.60 | 0.91 | 1.31 |
| share of tr D_b in the response subspace | 0.25 | 0.32 | 0.44 | 0.61 |

The defect operator is concentrated, and partly in the response subspace.

Full-strength localization, final MSE ratio (the first-order prediction in brackets):

| | k = 1 | k = 2 | k = 4 |
|---|---|---|---|
| response directions | 0.60 [0.43] | 0.50 [0.28] | 0.35 [~0] |
| D-oracle (D_b's top-k eigenvectors) | 0.51 [0.18] | 0.27 [<0] | 0.10 [<0] |

First order is optimistic. The defect decays along the path, the same effect that makes beta ~ 1 in Theorem 2.

The response share of the bias-weighted defect, `code/f3_share.py`, `out_share.txt`:

| net | s(1) | s(2) | s(4) | s(8) | s(16) | s(32) | random direction share x n |
|---|---|---|---|---|---|---|---|
| n32 L8 s0, s1 | 0.18, 0.28 | 0.36, 0.47 | 0.37, 0.54 | 0.59, 0.70 | 0.89, 0.92 | 1.00 | 0.8, 1.3 |
| n64 L8 s0, s1 | 0.23, 0.25 | 0.36, 0.32 | 0.43, 0.44 | 0.57, 0.61 | 0.83, 0.82 | 1.00, 0.94 | 0.9, 1.2 |
| n64 L16 s0, s1 | 0.08, 0.31 | 0.32, 0.34 | 0.38, 0.49 | 0.56, 0.66 | 0.81, 0.81 | 0.99, 0.95 | 1.1, 1.4 |
| n128 L8 s0, s1 | 0.17, 0.13 | 0.16, 0.22 | 0.27, 0.40 | 0.50, 0.50 | 0.63, 0.69 | 0.80, 0.86 | 0.9, 0.4 |
| n128 L16 s0, s1 | 0.41, 0.32 | 0.58, 0.40 | 0.75, 0.56 | 0.80, 0.69 | 0.81, 0.87 | 0.91, 0.95 | 0.6, -0.2 |

Here s(k) = sum_{j <= k} D_b(u_j, u_j) / tr D_b over the top response directions u_j.
- The defect has a collective part (the outlier response directions carry 13-41% of it in one direction).
- It also has an isotropic bulk: random directions carry about tr D_b / n each.
- Over doubling of n the bulk share grows: s(8) falls from 0.57-0.70 to 0.50 at L = 8.

---

## 3. The collective subspace at width 1024

**Theorem-sketch 6 (Lyapunov effective dimension).** J_L, the chain's own output first-chaos map, is a product of L
mean-gated matrices diag(Phi_l) W_l. In the free-probability limit:
- the normalized J_L^T J_L is the free multiplicative convolution of L factors of mean 1 and variance v, so its
  second moment is m_2 = 1 + v L (the S-transform is S = 1 - v z + ..., and S^L gives variance L v);
- its participation ratio is PR = n / (1 + v L).

Status: sketch. The factors are not exactly free or identically distributed, because the gates depend on the means,
and v is fitted.

*Measured* (`code/f3_pr.py`, `code/f3_pr4.py`, the chain's own mean-gated Jacobians; `out_pr.txt`, `out_pr4.txt`):

| n, L | PR(J_L) | v fitted | 2nd moment: top 1 / 4 / 16 / 64 | 4th moment: top 1 / 4 / 16 | 6th moment: top 1 / 4 |
|---|---|---|---|---|---|
| 64, 8 | 4.4 | 1.7 | 0.40 / 0.81 / 1.00 / 1.00 | 0.68 / 0.97 / 1.00 | 0.83 / 0.99 |
| 128, 8 | 7.0 | 2.2 | 0.30 / 0.64 / 0.96 / 1.00 | 0.61 / 0.91 / 1.00 | 0.80 / 0.98 |
| 256, 8 | 14.0 | 2.2 | 0.17 / 0.44 / 0.84 / 1.00 | 0.41 / 0.76 / 0.98 | 0.61 / 0.91 |
| 64, 16 | 2.0 | 1.9 | 0.66 / 0.97 / 1.00 / 1.00 | 0.85 / 1.00 / 1.00 | 0.90 / 1.00 |
| **1024, 16** | **25.9** | **2.4** | **0.14 / 0.31 / 0.64 / 0.96** | **0.46 / 0.69 / 0.92** | **0.74 / 0.90** |

- **It agrees with the repo.** The per-layer PR at n = 1024 is 512, 287, 188, 140, 111 (layer 5), ..., 46 (layer 11),
  ..., 26. Frontier-tests measured 120 at layer 5 and 50 at layer 11 for single births, and 95% of the energy in the
  top 64 at age 15.
- **Higher moments concentrate on the outliers.** At width 1024 one direction carries 46% of the fourth moment, the
  moment that governs the gain's quadratic form x^T J^T J x and the 4-walks. Four directions carry 69%, and 90% of the
  sixth moment. In this respect n = 1024, L = 16 resembles n = 256, L = 8. That is the regime where k = 1 localization
  removes 31% of the Gaussian closure's MSE (section 4).
- **Consequence.** The collective part of a closure's error, and the U-visiting walk classes of Theorem 3, do not
  dilute with width as fast as the participation ratio suggests. **Capacity is not what kills F3 at width 1024. Cost
  is** (section 5).

---

## 4. Numerics: what a k-dimensional localization removes

Exact leaves (the estimator run on N(Uy, I - UU^T)), product Gauss-Hermite with 7, 5 and 3 nodes per dimension for
k = 1, 2, 4.
- Convergence: GH-9 against GH-7 and GH-7 against GH-5 agree to 0.5% at k <= 2; GH-3^4 against GH-4^4 differ by
  1-6% at k = 4.
- U is the top-k right singular vectors of the chain's own first-chaos map (resp) or Haar random (rand).
- The ratio is localized MSE over the same estimator's root MSE, final layer, against Monte Carlo truth (2e6
  samples; noise <= 3e-3 of the smallest MSE).
- Files `out_loc*.txt`, `out_tlp*.txt`.

| n, L (nets) | root MSE: GC (TLP) | GC resp, k = 1 / 2 / 4 | GC rand, k = 1 / 2 / 4 | GAC resp, k = 1 / 2 / 4 | TLP resp, k = 1 / 2 |
|---|---|---|---|---|---|
| 32, 8 (4) | 1.0-2.1e-3 | 0.58 / 0.45 / 0.26 | 0.97 / 0.88 / 0.78 | 0.90 / 0.51 / 0.17 (one net; the others give NaN) | |
| 64, 8 (4) | 5.5-7.4e-4 (0.75-1.4e-4) | 0.58 / 0.44 / 0.32 | 0.97 / 0.93 / 0.83 | 0.60 / 0.44 / 0.31 | 0.58 / 0.38 (k = 4: 0.27, two nets); full memory: 0.62 / 0.37 |
| 128, 8 (3) | 1.1-1.4e-4 (3.1e-5) | 0.78 / 0.68 / 0.50 | 0.98 / 0.96 / 0.91 | 0.98 / 0.75 / 0.53 | 0.72 / 0.52 |
| 256, 8 (2) | 3.5-5.4e-5 | 0.69 / 0.62 / 0.54 | 1.00 / 0.99 / 0.97 | | 0.72 / - |
| 64, 16 (3) | 3.9-6.4e-4 | 0.74 / 0.58 / 0.37 | 0.97 / 0.93 / 0.83 | 0.87 / 0.51 / 0.31 (two nets) | 0.73 / 0.47 |
| 128, 16 (2) | 1.3-2.4e-4 | 0.59 / 0.46 / 0.33 | 0.99 / 0.97 / 0.92 | | |

**Reading.**
1. **Response directions are the right ones; random ones do nothing.** A random direction removes about 1/n of the
   error, the isotropic share of tr D_b. The MC Stein Jacobian, the chain's own Jacobian and the stacked layer
   Jacobians give the same subspace to within the spread (n = 64: 0.53, 0.58, 0.54 at k = 1, `out_loc_n64L8.txt`).
   The top-W0 directions are in between (0.82).
2. **For the Gaussian closure the effect does not fade with width at fixed depth as the participation ratio would
   suggest** (n = 128 to 256: 0.78 to 0.69 at k = 1). The closure's error is gain-dominated, and it tracks the
   fourth-moment concentration of section 3.
3. **It survives better closures.**
   - TLP, four to seven times more accurate than GC, loses the same share as GC or a larger one (n = 64, 128).
   - With TLP's full source memory (rmax = L - 1, MSE 1.7-3x lower again) the localizable share at n = 64 is
     unchanged (0.62 / 0.37, `out_tlp_rmax.txt`).
   - TLP's removed share at k = 1 is 0.42 (n = 64), 0.29 (128) and 0.28 (256): it settles near 0.3 from n = 128 on,
     tracking the fourth-moment concentration as GC does.
   - With full source memory at n = 128, where the MSE is halved again, rho is 0.76 / 0.68 at k = 1 and 0.46 / 0.43 at
     k = 2 (two nets).
   - A tree-level closure at n = 1024, L = 16, which has the fourth-moment profile of n = 256, L = 8, would lose about
     25-30% of its MSE to a k = 1 localization.
   - GAC, which has removed the radial/gain part, is the one estimator where k = 1 is null at n = 128 (0.98). Its
     k = 2-4 shares remain.
   - The localized part is therefore not only the gain. Theorem 3 says what else it is: the U-visiting closed walks,
     opened into paths that a tree-level closure handles.
4. **The leaf-space dimension does what a smooth transversal can do and no more.** Even k = 4 with the D-oracle at
   n = 64 leaves 0.10 (Proposition 5), and the isotropic bulk of D_b is untouched by any small k.

**What the production chain's residual is.** At n = 1024 the production chain is 1/200 of GC (2.2e-8 against 4.06e-6,
note XII). Its residual is the incoherent flat-sector kappa3 transport error and the flex-limited kappa4 slices
(notes XXXIX, XLI). Decoding it through the fresh weights finds no low-dimensional structure (frontier-tests). That
statement concerns activation space at each layer. In input space Theorem 1 locates it through the chain's own heat
defect, which has not been measured for the production chain. Hence the decisive experiment of section 7.

The prior for the production chain is **rho_1 in [0.70, 1.00], central value 0.85.** The small-width evidence pulls
toward 0.72: GC and TLP, with or without full source memory, sit at 0.68-0.78 from n = 128 to 256. GAC pulls toward 1
at k = 1. A structural argument says the production chain should sit above the tree-level value.
- The production chain carries the young sources densely. These are deep births, and their legs live in the
  collective subspace (PR(L_b) ~ n/(1 + v b): 46 at b = 11).
- Its measured error sits in the old tier and in the flat sector of the kappa3 transport. Old sources are shallow
  births, whose legs are spread over 100-500 input directions (PR 512, 287, 188, 140, 111 at b = 1-5).
- That is the part of D_b a k <= 2 localization does not reach (the isotropic bulk of Proposition 5).

---

## 5. Cost

**Lemma 7 (no free leaves). Proved.** Suppose the leaf states of a chain are affine in y with y-independent
covariances: means m_l + A_l y, covariances C_l - A_l A_l^T, with gates and all nonlinear coefficients frozen at the
root. Then the localized Gaussian readout equals the root readout exactly.

*Proof.* E_y E relu(m + A y + (C - A A^T)^{1/2} Z) = E relu(N(m, C)). QED.

So the leaves' value lies entirely in recomputing the nonlinear steps (gates, Mehler coefficients, source births) per
leaf. After one layer the leaf covariance differs from the root by the Hadamard terms
(A_k(y) A_k(y)^T - A_k A_k^T) o R^k/k!. These are generically full rank, so the next W^T (.) W costs a full n^3 per
leaf. The kappa3 young tier (about 70% of the bill, note XXXIX section 11g) must likewise be re-run per leaf. A
second-order tangent of the chain along u costs at least the primal plus two tangent sweeps. **A leaf costs a chain.**

**Theorem 8 (cost floor). Proved.** Expand the leaf function y -> G(Uy, I - P) in y. Its degree-2 part carries
(1/2) tr(Hess_m G P), which is what compensates the covariance downdate -tr(d_Sigma G P). For the exact functional the
two cancel identically (Theorem 1 with D_T = 0).

A transverse rule that does not integrate quadratics in y exactly is therefore inconsistent at leading order: it is
wrong by O(tr(Hess_m T P)) even when G = T. For example, the one-node rule y = 0 freezes the leaf, and the frozen leaf
has the reduced variance without the compensating spread of means.

Degree-2 exactness on gamma_k needs N >= k + 1 nodes, since the moment matrix of the degree-<=1 monomials must be
nonsingular. The bound is attained by the regular simplex rule; Stroud's 2k-point rule gives degree 3. With Lemma 7:
**cost(k-localization) >= (k + 1) x cost(leaf chain).**

**Theorem 9 (costed bound on an elasticity-1 frontier). Proved under the stated hypothesis.** Hypothesis: near the
operating point the chain's raw MSE scales as r(C) = r_0 C_0 / C. This is the measured "elasticity about 1" of the
tiers (note XXX section 7, note XXXIX section 7). Let an N-leaf localization use leaves of cost c each and achieve
raw rho r(c). Then for every c,

    adjusted(localized) / adjusted(base) = rho r_0 (C_0/c) max(0.1 B, N c) / (r_0 max(0.1 B, C_0)) >= rho N

when C_0 >= 0.1 B. Equality holds whenever N c >= 0.1 B. **An N-leaf localization beats the base in adjusted MSE iff
rho < 1/N**, however the budget is split between leaves. In particular:
- **the "free" zone below 0.1 B does not help.** Leaves at 0.1 B / N each are N times worse at the margin, which
  cancels the free cost exactly;
- below the young tier's cost floor the elasticity exceeds 1 (four youngest sources only: 7.3e-7, frontier-tests),
  which makes small leaves strictly worse.

**The minimal rules, measured** (`code/f3_minrule.py`, `out_minrule.txt`; GC, response directions):

| net | k = 1: 2-point (GH-7) | k = 2: simplex-3, Stroud-4 (GH-5) | k = 4: simplex-5, Stroud-8 (GH-3^4) |
|---|---|---|---|
| n64 L8 s0 | 0.694 (0.618) | 0.741, 0.684 (0.465) | 0.655, 0.524 (0.339) |
| n64 L8 s1 | 0.632 (0.602) | 0.911, 0.617 (0.495) | 0.828, 0.563 (0.347) |
| n64 L8 s2 | 0.527 (0.495) | 0.403, 0.390 (0.344) | 0.660, 0.400 (0.268) |
| n64 L8 s3 | 0.631 (0.607) | 0.421, 0.472 (0.459) | 0.300, 0.448 (0.317) |
| n128 L8 s0 | 0.773 (0.756) | 0.807, 0.886 (0.763) | 1.002, 0.763 (0.593) |

rho N with the minimal rule is:
- 1.05-1.55 for k = 1, N = 2;
- 1.2-2.7 for k = 2, N = 3 (1.6-3.5 with Stroud-4);
- 1.5-5.0 for k = 4 with simplex-5, 3.2-6.1 with Stroud-8.

**It is above 1 on every network at every width tested, including widths where localization is most effective.**
Even the exact product rules (N = q^k) fail: rho N >= 7 x 0.5 at k = 1. The costed bound holds with margin.

**The answer to the frame's question.**
- **At <= 0.2 B:** no regime. One production chain is 0.203 B, so any localization is >= 0.406 B and needs
  rho_1 < 0.5. Predicted rho_1 is 0.70-1.00; the best case observed at any width is 0.40-0.49 on single small nets,
  and the typical case is 0.56-0.78.
- **At <= 0.1 B:** no regime. By Theorem 9, leaves at 0.05 B need the same rho < 1/2, and the chain at 0.05 B is below
  the young-tier floor.

---

## 6. The design (the best F3 system), costed

**System "Lyapunov leaves" (LL-k).**
1. **Collective direction.** u_1 is the top right singular vector of the chain's own first-chaos map J_L. Five power
   iterations on J_L^T J_L reuse the transported mean-gate maps the chain already carries: 2 x 5 x 16 matvecs,
   about 0.005 units.
2. **Leaves.** For k = 1, two nodes y = +-1 (degree 3). Run the production chain from mu = y u_1,
   C = I - u_1 u_1^T: the first layer's pre-activation law is N(y u_1 W_0, W_0^T W_0 - v v^T) with v = W_0^T u_1, and
   is exact. The rest of the chain is unchanged.
3. **Combine.** The estimate is (Chat(+u_1) + Chat(-u_1))/2. The counterterms are refitted on the leaf-averaged
   output, since the per-layer amplitudes were fitted at the root.
4. **Cost.** 2 x 203 units + 0.01 = 0.406 B. For k = 2 (simplex-3): 0.61 B. For k = 4 (Stroud-8): 1.62 B.

**Predictions (official networks, production chain, free-running):**
- LL-1: raw 1.09-1.55e-8 (rho_1 = 0.70-1.00, central 1.32e-8), C/B 0.406, adjusted **4.4-6.3e-9 against the base
  3.14e-9: worse by 1.4-2.0x**;
- LL-1 with GH-3 nodes (0, +-sqrt 3), 0.61 B: raw within 3% of LL-1, adjusted 2.4-3.0x worse;
- the same localization on the Gaussian closure at n = 1024: rho_1 = 0.60-0.80 (it tracks n = 256, L = 8: 0.69,
  through the fourth-moment concentration of section 3). This calibrates the instrument: if the GC value falls outside
  this band, the small-width extrapolation is wrong.

**Self-verdict: negative** for the scored system. F3's ingredients are not useless:
- the response subspace is the right object;
- the walk-opening theorem says what it does;
- the Duhamel identity measures it.

But the leaf is a chain, and the frontier is elastic.

---

## 7. The decisive minimal experiment

**F3-D (AWS, about 0.5 instance-hour).** Inputs: official networks 0-7 (`W_off{k}.npy`, `truth_off{k}.npz`) and the
production chain `estimator_final_v56.py`, free-running, counterterms on and off.

0. **Homogeneity check.** G(2 m, 4 Sigma) = 2 G(m, Sigma) on one network, root law: assumption (A3). Theorems 1 and
   8-9 do not need it; only Theorem 2 does.
1. **Patch** (lead). An environment hook `F3_U=<path to u.npy>`, `F3_Y=<float>` sets the input law. In `predict`,
   replace `mu = zeros(n)`, `C = eye(n)` by `mu = y u`, `C = eye(n) - u u^T`. Check that the li = 0 `C_pre` branch
   takes the general `einsum("ij,ia,jb->ab", C, w32, w32)` path when C is not the identity (est_v29 lines about
   1164-1178). Also check that the identity path bills the same.
2. **Directions.** u_1, u_2 are the top right singular vectors of the chain's J_L (dump the transported mean-gate map,
   or form J_L in numpy offline: unbilled, diagnostic), and u_r is a Haar random control.
3. **Runs per network (14 chains):**
   - the root;
   - u_1 at y in {0, +1, -1, +sqrt 3, -sqrt 3} (y = 0 is the frozen leaf N(0, I - u u^T), not the root);
   - u_2 at y = +-1;
   - u_r at y = +-1;
   - the Gaussian closure (`code/f3_lib.gc_chain`, float64) at the root and at u_1, y = +-1, as the instrument
     calibration;
   - optional F1 add-on: 8 Hutchinson probe pairs G(+h z_j), G(-h z_j) for b . Lap_m G (16 more chains).
4. **Outputs per network:**
   - final-layer MSE of the root, the 2-point leaf average, the GH-3 average (weights 2/3 on the y = 0 leaf and 1/6 on
     each +-sqrt 3 leaf) and the k = 2 cross average;
   - rho_1, rho_2 and rho_rand;
   - per-layer MSE;
   - the leaf second difference delta_1 = (Chat(+u) + Chat(-u))/2 - Chat(0) and its alignment
     (b . delta_1)/|b|^2 with b = Chat(0) - truth, the measured U-block of Theorem 1;
   - the GC instrument rho_1(GC);
   - optional: beta_hat = |b|^2 / (b . (G - Lap_hat G)/2) for frame F1.
5. **Compute.** 14 x 8 = 112 production-chain runs at a few minutes each on one core: a few core-hours, or about
   10 minutes on the 48-slot fleet.

**Decision rule** (pre-registered here):
- **rho_1 (2-point, mean over 8 networks) < 0.5, better on >= 7/8:** F3 is alive and contradicts this note's
  prediction. Run LL-1 in the scored regime with refitted counterterms; adopt by note XXXIX's rule. Then test k = 2.
- **0.5 <= rho_1 <= 0.9:** the chain's residual has a real collective part (record (b . delta_1)/|b|^2 per layer),
  but F3 does not pay (Theorem 9). Hand the U-block to frames that can carry it inside one chain. A per-layer
  counterterm along delta_1 is cheap to test, since delta_1 needs two chains and would have to be learned across
  networks.
- **rho_1 > 0.9:** close F3 for Phase 2. The prediction, central 0.85, falls in the middle band.
- **Instrument check.** If rho_1(GC) falls outside 0.60-0.80, the small-width extrapolation is wrong, and the verdict
  rests only on the production-chain number.

---

## 8. What this frame establishes, and what it hands on

1. **The leaf space of a linear localization is a manifold, the type-I corner of noncommutative disintegration.**
   Deterministic transverse quadrature is ordinary Gauss-Hermite. Nothing noncommutative is needed, and nothing
   noncommutative helps.
2. **The bias trickle-down (Theorem 1) is the right accounting.** A localization removes exactly the path-integrated
   U-block of the estimator's heat defect. Homogeneity turns the isotropic version into the Euler identity
   E F = E Lap F and the midpoint estimator (G + Lap_m G)/2, which removes 73-98% of the Gaussian closure's MSE at
   small width (Theorem 2). For F1: the open problem is a sub-2n-chain evaluation of Lap_m of a chain.
3. **The walk-opening theorem (Theorem 3).** Conditioning on U turns once-visiting closed walks into leaf paths and
   multiply-visiting ones into quadrature. The closed walks concentrate on the few outlier response directions:
   - fourth-moment share 0.46 in one direction at width 1024;
   - the top-4 response directions touch 69% of tr H^3 at n = 32, and 56% at n = 64.

   So the closed-walk (triangle) wall has a low-dimensional input-space part. The costs inside one chain:
   - the inside-U piece, tr (U^T H_i U)^3, needs only the legs' U-coordinates: k^2 n^2 per source-layer;
   - the once- and twice-visiting pieces need H_i U: k n^3 per source-layer, the price of k extra legs.

   Whether the cheap inside-U piece is worth anything to the production chain depends on how much of the chain's error
   is that omission. It is a cheap diagnostic for whoever attacks the triangle class next (not F3's design).
4. **The cost theorem (Theorems 8 and 9)** is general. On an elasticity-1 frontier, splitting the budget into N leaves
   of any kind (localization nodes, mixture components, pattern windows) pays only if it cuts raw MSE N-fold. This also
   explains note XII's finding that mixtures were "priced out of the contest" at M = 4.

## 9. Files

- Development and checks: `code/f3_lib.py` (nets, Monte Carlo truth with Stein Jacobians, Gaussian closure on
  general input, Gauss-Hermite and Smolyak grids, localization), `f3_gac.py` (GAC on general input; repo code
  imported read-only), `f3_tlp.py` (tree-level propagation on general input, plus its localization experiment),
  `f3_truth.py`.
- Experiments: `f3_loc.py`, `f3_loc2.py` (localization experiments), `f3_defect.py` (heat defect, Duhamel checks,
  D-oracle), `f3_euler.py` (Theorem 2 estimators), `f3_share.py` (defect shares), `f3_walk.py` (Theorem 3),
  `f3_nngp.py` (chaos weights), `f3_pr.py`, `f3_pr4.py` (spectral concentration), `f3_minrule.py` (minimal rules),
  `f3_hwalk.py` (closed-walk shares on exact second-chaos kernels), `f3_tlp_r.py` (TLP with full source memory).
- Outputs: `out_*.txt` beside this file; data in `data/` (tiny networks and their truth).

---

## Referee report

Adversarial referee, frame F3. Scripts and outputs: `ref_f3/` (`chk.py`, `chk2.py`, `quaderr.py`, `qe_*.txt`).

**Bottom line.** The negative verdict survives. The exact identities are
correct but mostly classical. The cost theorem T9 is right in its algebra, but it rests on a hypothesis stronger than
anything measured. Its slogan "the free zone below 0.1 B does not help" holds only under a frontier-shape condition
that nobody has checked. The prediction band for rho_1 is too optimistic and too narrow. The pre-registered
experiment has three specification errors. None of this rescues F3 for Phase 2.

### R1. Independent checks (what is right)

- **T1 (Duhamel).** Checked on a smooth G that is deliberately not a heat solution (n = 4, k = 1, GH-60 transverse
  rule, finite-difference D_G; `chk2.py`): d/dc G_U^c = -3.461700 / -2.349641 / -1.669763 against
  -E tr(D_G P) = -3.461699 / -2.349641 / -1.669761 at c = 0.05 / 0.3 / 0.7. Correct. (A first attempt with an
  oscillatory test function disagreed at c = 0.7; that came from under-resolved quadrature, not from the theorem.)
- **T3 (walk opening).** Exact total-cumulance decomposition at n = 9, k = 3 (`chk.py`): kappa3 = 6.114106105485 on
  both sides. The cross term 3 Cov(E[z|y], Var(z|y)) = 3 (2 L_u^T H_up L_p + tr(H_uu H_up H_pu)), checked by Monte
  Carlo (3.04 against 3.01, 2e6 samples).
- **T2 / P2b.** The Euler identity, tr D_G(0,I) = (G - Lap G)/2, the scaling D_G(0,sI) = s^{-1/2} D_G(0,I), and the
  first-order exactness of the midpoint for G = T(m,(1+eps)Sigma) all re-derive cleanly.
- **P4, P5, L1, T8 and the free-probability second moment of S6** (var(mu ⊠ nu) = var mu + var nu for mean-1
  factors, from S'(0) additivity) are correct.
- **Table entries** checked against `out_loc2_*`, `out_tlp_*`, `out_minrule.txt` and `out_euler.txt`. The rho N
  ranges (1.05-1.55, 1.2-2.7, 1.5-5.0 / 3.2-6.1) are reproduced from `out_minrule.txt`.

### R2. Mathematical gaps and overclaims

1. **T9's hypothesis is not what was measured.** The repo's "elasticity 1" (note XXX section 7.4) is a *local*
   stationarity of raw x C at the operating point. Measured there: dlog MSE / dlog r = -0.30 against
   dlog C / dlog r = +0.33. T9 instead assumes r(C) = r_0 C_0 / C globally, down to 0.05 B.
   - The bound adjusted(loc)/adjusted(base) >= rho N actually needs only r(c) c >= r_0 C_0 for every leaf cost c,
     that is, the base configuration minimizes raw x C over its family. Under that condition the proof goes through
     in both zones. I re-derived it.
   - Without that condition the theorem fails. If the chain's raw were flat on [0.05 B, 0.1 B], two leaves at 0.05 B
     would give ratio = rho, not 2 rho.
   - With a local power law r ~ C^{-e}, the condition reads rho < N^{-e}. For e = 0.5 and N = 2, the bar is
     rho < 0.71, which lies **inside** the author's predicted band.
   - So "the free zone does not help" depends on r(0.05 B) / r(0.1 B) for the production chain, and this has not
     been measured. The leaders, at raw 1.1-1.5e-8 for 0.11-0.15 B, show that chain frontiers below 0.2 B are not
     pinned by our chain's local slope.
   - T9 should be restated with the r(c) c >= r_0 C_0 hypothesis, and the frontier scan should be part of the
     decisive experiment (R5).
2. **The 2-point rule has its own quadrature error on the exact functional. It matters at small width, and it is
   probably negligible at n = 1024.** (New measurement, `quaderr.py`.)
   - Definition: Q = (T(u, I-P) + T(-u, I-P))/2 - E_y T(yu, I-P), with u the top right singular vector of the Monte
     Carlo Stein Jacobian. Computed with common random numbers, a GH-24 reference and a debiased mean square.
   - Algebra: Q = -2 c_4 + 16 c_6 - 132 c_8 + ..., where c_k are the Hermite coefficients of F along u
     (He_4(1) = -2, He_6(1) = 16).
   - Results (relative = Q^2 / sigma^2):

     | n, L, seed | 2-point: Q^2 (relative) | GH-3: Q^2 (relative) | GH-5 |
     |---|---|---|---|
     | 32, 8, 0 | 6.9e-6 (2.5e-5) | 9.2e-8 (3.4e-7) | 0 within noise |
     | 64, 8, 0 | 1.9e-6 (1.3e-5) | 3.2e-8 (2.1e-7) | 0 |
     | 64, 16, 0 | 2.0e-7 (4.5e-6) | 1e-8 (2e-7, at the noise level) | 0 |
     | 128, 8, 0 / 1 | 4.6e-8 / 2.0e-7 (2.1e-7 / 1.2e-6) | 0 / 1.4e-9 | 0 |
     | 128, 16, 0 / 1 | 1.3e-7 / 3.3e-7 (6.7e-7 / 3.4e-6) | 0.8e-9 / 4.4e-9 | 0 |
     | 256, 8, 0 / 1 | 1.5e-8 / 1.7e-8 (6.3e-8 / 1.3e-7) | 0 | 0 |
     | 256, 16, 0 | 1.75e-8 (1.3e-7) | 0 | 0 |
     | 512, 16, 0 (N = 64000) | 3.5e-10 +- about 1e-10 (3.9e-9) | 0 | 0 |

   - Reference level: the production chain's own relative MSE is 1.55e-8 / 0.074 = 2.1e-7.
   - At n = 256 the 2-point rule's error **on the exact function** is 0.3-0.6 of that level.
   - At n = 512 it has fallen to about 2% of it. The fall from n = 256 to 512 is steep (one seed), and it extrapolates
     to Q^2 <~ 1e-10 at n = 1024, under 1% of the base MSE.
   - Consequences:
     - The small-width 2-point rho values (`out_minrule.txt`, n = 64-128) are barely affected. For GC and TLP the
       floor is 0.1-4% of their MSE.
     - At n = 1024 the floor most likely does not decide anything.
     - The note's statement "GH-3 raw within 3% of LL-1" is plausible at n = 1024 for this reason, but it was not
       derived. E1 below is a cheap direct check, and I downgrade it from a kill test to a sanity test.
   - This item does **not** strengthen the negative verdict at n = 1024. I report it because the first runs at
     n <= 256 looked damning, and they are not representative.
3. **The prediction band for rho_1 ignores the evidence that a gain-aware estimator gains nothing or loses.**
   - The production chain carries the gain mode explicitly: the wk4m3 slice sqrt(g4_i g4_j)/3 is "the gain mode's
     var_i var_j term exactly" (leaderboard note, section 6). Its closest small-n proxy is GAC, not GC or TLP.
   - GAC's k = 1 ratios per net: 1.549 (n32 s0), 1.347 (n64 L16 s2), 1.011 and 1.114 (n128 L8 s0, s1), mean 0.98 at
     n = 128. **Localization made GAC worse on 4 of 14 net instances.**
   - The band [0.70, 1.00] with centre 0.85 should be about [0.85, 1.3] with centre near 0.95-1.0. That is before the
     quadrature term of R2.2, which adds to it.
   - Misreport: the section-4 table says GAC n = 32, k = 1 is "0.90 (one net; the others give NaN)". In fact the
     0.904 is the mean over all four nets, one of them 1.549. The NaN occurs only at k = 4.
4. **"A leaf costs a chain" is proved only for fully frozen leaves.**
   - Lemma 7 is a near-tautology. With every coefficient frozen, the chain is linear in the state, so averaging
     commutes with it.
   - Partially recomputed leaves are not covered. Examples: re-gating the mean and diagonal per leaf, and carrying
     the leaf covariance as the root plus Hadamard-rank-structured corrections.
   - The "generically full rank, so n^3" step is a heuristic, not a theorem. It is a plausible one, since the
     corrections (a b^T ∘ R^k) need R^k (diag(b) W) at n^3 per layer.
   - The status should read "proved (frozen); sketch (general)".
   - Input-space localization is also not the only option. A branch at layer l* costs (L - l*)/L of a chain per
     leaf. It is not considered, although the old-tier residual is born shallow, so this probably does not save F3.
5. **T1(b) integrability** needs D_G(., I) bounded as |m| -> infinity, that is, near the point-mass limit. This is
   an extra hypothesis, not implied by (A1)-(A3).
6. **S6 "outliers".** The concentrated top singular directions are the edge of a heavy-tailed Fuss-Catalan bulk,
   not separated outliers. At fixed L the top direction's share of every moment goes to 0 as n -> infinity.
   "Does not dilute with width" is a finite-n statement. It is true at n = 1024 by direct measurement, which is all
   that matters here.
7. **The summary puts "at width 1024" next to the closed-walk shares (56-69%).** Those shares were measured only at
   n = 32-64.

### R3. Cost-accounting errors (minor; none changes the verdict)

- "2 x 203 units" should read 2 x 208 units = 416 units = 0.406 B. Here 0.203 B = 208 units; the note mixes
  thousandths of B with units.
- Power iteration: 2 x 5 x 16 = 160 dense n x n matvecs = 160 x 2n^2 FLOPs = **0.156 units**, not 0.005. It is still
  negligible (1.5e-4 B).
- The leaf's first layer takes the general einsum path (about 2 extra units) unless the rank-1 downdate
  W^T W - v v^T is coded directly (n^2). This is negligible, but it should be stated, because the experiment's
  patch forces the general path.

### R4. Specification errors in F3-D

1. **The "k = 2 cross average" cannot be computed from the listed runs.** Leaves u_1 at +-1 and u_2 at +-1, each with
   C = I - u u^T, are not a 2-D localization. That needs C = I - U U^T with U = [u_1, u_2] and nodes in R^2 (simplex-3:
   three more chains per net).
2. **The delta_1 diagnostic uses the wrong bias.** b = Chat(0) - truth with Chat(0) the *frozen* leaf
   N(0, I - u u^T) is dominated by the leaf's missing variance along u, so b . delta_1 / |b|^2 is about 1 by Theorem 8,
   whatever the chain does. Use the root bias b_root and Delta = (Chat(+u) + Chat(-u))/2 - Chat(root). Report
   rho = 1 + 2 b_root . Delta / |b_root|^2 + |Delta|^2 / |b_root|^2 term by term.
3. **The input patch has hidden-assumption risk.** The note XLI chain grades the input chaos against the *standard*
   Gaussian (kappa3 sources as second chaos with legs built from W columns). A leaf input N(yu, I - uu^T) needs the
   legs re-whitened (Sigma^{1/2} W), not only the C_pre einsum. Required null tests:
   - (a) y = 0 with C = I reproduces the root to the last digit;
   - (b) a weak localization, C = I - eps u u^T at nodes +-sqrt(eps): the change must be linear in eps, with slope
     -tr(D_G P) estimable by a 3-point fit;
   - (c) per-leaf MSE against a per-leaf truth on one net. The leaf truth can be generated cheaply with common random
     numbers, as in `quaderr.py`.

### R5. Corrected minimal decisive experiment (runs before F3-D, cheaper than it)

- **E0 (frontier shape; decides T9).** Production chain on nets 0-7 at C/B ≈ 0.05, 0.1 and 0.2, using the existing
  tier/rank switches; 24 chains. If r(0.05 B) x 0.05 >= r(0.2 B) x 0.2 and r(0.1 B) x 0.1 >= r(0.2 B) x 0.2, T9 holds
  in its corrected form and F3 needs rho_1 < 0.5. If r(0.05 B) / r(0.1 B) < 2, record the local exponent e; the bar
  becomes rho_1 < 2^{-e}.
- **E1 (quadrature floor, no chain).** For nets 0-7 and u_1, measure Q_2pt and Q_GH3 for the *true* network:
  - method: common-random-number MC, 26 forward passes per sample, about 1e5 samples, so about 2.6e6 passes,
    about 9e13 FLOPs per net (minutes on the fleet);
  - sanity rule (it was a kill rule before n = 512 was measured): if mean Q_2pt^2 >= 0.25 x 1.55e-8, LL-1 with 2
    points is dead regardless of the chain, and only GH-3
    or higher remains (needs rho < 1/3).
- **E2.** F3-D with the R4 fixes. The decision rule otherwise stands, with the bar set by E0.

My prediction:
- E0: elasticity at or above 1 below 0.1 B (the young-tier floor);
- E1: Q_2pt^2 <~ 1e-10 (under 1% of base);
- E2: rho_1 ≈ 0.9-1.1 for both the 2-point and the GH-3 rule.

F3 closes.

### R6. Novelty

- T1 is the classical Gaussian-interpolation (heat-semigroup Duhamel) formula, written as a defect integral. It is
  the same calculation as Itô's formula for G(m_t, Sigma_t) along Eldan's localization, or the smart-path
  interpolation of Chatterjee and Talagrand. Calling it "bias trickle-down" is apt but not new.
- T3 is the law of total cumulance for Gaussian quadratic forms. P4 is Stein plus the gradient-based dimension
  reduction of Zahm et al. P5 is Ky Fan. T8 is the classical lower bound on cubature nodes (Mysovskikh / Stroud).
- The NCG reading is honestly labelled as thin (a type-I groupoid), and I agree.
- **Genuinely new here:**
  - the homogeneity midpoint (G + Lap_m G)/2, with the measured path factor beta ≈ 1, which removes 73-98% of GC's
    error at n = 64-128;
  - the empirical finding that the response subspace carries a collective share of closure bias, while random
    directions carry only tr D_b / n;
  - the corrected pricing law of R2.1 for any N-way budget split.

### R7. Relevance at the 1e-8 level

None. Every F3 design costs at least 0.41 B. The best plausible rho_1 for the gain-aware production chain is about
0.9-1.0, against a bar of 0.5 (or 0.71 at best if E0 finds e = 0.5). The midpoint by-product
needs Lap_m of the chain, at 2n chains or an n^4 forward Laplacian. A cheap approximate Laplacian is useless, because
tr D = (G - Lap G)/2 is a difference of O(1) quantities that must be accurate to about 1e-5 per neuron.

### Verdict

**Survives as a negative result**, with the corrections above:
- T9 restated under the hypothesis r(c) c >= r_0 C_0;
- the prediction band widened to [0.85, 1.3];
- the 2-point rule's quadrature floor measured (small at n = 1024);
- the three experiment-specification fixes.

The strongest idea that survives is the corrected T9 pricing law, combined with the Duhamel defect identity used as
an offline diagnostic. Together they say: splitting a total C* >= 0.1 B into N leaves pays only if rho < r(C*) / r(C*/N),
which is at most 1/N on a raw x C-optimal family. Below n = 512 a transverse rule needs degree >= 5 exactness (GH-3) on the true function. At n = 1024 the
degree-3 rule's floor appears to be under 1e-10.

Priority for Phase 2: **1/10**.
