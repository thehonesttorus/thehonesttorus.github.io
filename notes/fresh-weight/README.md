# XLIV. One interpolation under three papers, and the fresh-weight theorem

Status: in progress (9 October 2026). Sections 1-4 are the theory. Section 5 is the pre-registration of the decisive
test; it was committed before the test's data.

## 1. The request

The user sent three papers and a set of ideas from another chat:
- the Sahasrabudhe survey on probabilistic combinatorics at exponentially small scales (arXiv 2512.15077);
- Sah-Sahasrabudhe-Sawhney on the Spielman-Teng conjecture (2405.20308);
- Wang-Lau-Zhou on derandomizing matrix concentration through free probability (2601.08111);
- the other chat's ideas: a KMS half-tilt theorem, sech factors read as Chernoff affinities, Littlewood-Offord theory as
  the spectral theory of the modular Hamiltonian, and Klartag's identity read as Ito's formula for the log-det barrier.

The brief is theory first, with the Phase 2 adjusted MSE as the only judge, and no patching of the moment chain.
The user named the bridges to develop:
- codegree, the intersection of balls, read as the overlap of projections;
- approximate conditional expectations;
- Klartag's growing ellipsoid, and random lattices being locally Poisson;
- local random-walk analysis;
- time-evolving structure, the modular group of a faithful non-tracial state, and Hamiltonians on evolving graphs;
- the Littlewood-Offord concentration function as the object for activations.

The user also warned that means may drop out of some other object, and that the kappa3/kappa4 frame may be the wrong
one.

## 2. What each source does (read in full here; claims checked where marked)

### 2a. The survey (2512.15077)

- **The nibble, controlled by codegrees.** Sphere packing and spherical codes get a log d gain as follows:
  - sample a Poisson process in a box and join points whose balls meet;
  - forget the geometry and keep only the degree Delta and the codegree Delta_2;
  - find an independent set by a random greedy "nibble".

  In high dimension two balls meet in a fraction e^(-|x-y|^2/2) of their volume, so Delta_2 << Delta. The nibble's
  martingale increments are bounded by codegrees,
  sum_y |I ∩ N(y)|^2 <= sum_(y,z in N(x)) |N(y) ∩ N(z)|,
  so the process follows its mean-field trajectory, as on a tree.
- **Klartag's ellipsoid.** Klartag runs a random walk A_t on the space of ellipsoids {x : <x, A_t x> <= 1}:
  - a lattice point that hits the boundary sticks and imposes the linear constraint <y, A_t y> = 1;
  - the walk stops after about d^2/2 contacts;
  - a random lattice looks locally like a Poisson process of intensity 1, so the final volume is about d^2.
- **The concentration function.** rho_eps(v) = max_b P(|<X, v> - b| < eps). The survey covers:
  - Halasz's Fourier bound;
  - the inverse theorems (a large rho forces arithmetic structure in v);
  - the least common denominator (LCD);
  - Tikhomirov's typical Littlewood-Offord theorem: for random v, E_v rho_eps(v) = Theta(eps);
  - approximate negative correlation of small-ball events;
  - the rank-splitting proof for random symmetric matrices.

### 2b. Spielman-Teng (2405.20308)

The theorem is P(sigma_n(M) <= eps n^(-1/2)) = (1 + o(1)) eps + e^(-Omega(n)). It is proved in four moves:
1. **Secular reduction.** A rank-one update, Lemma 3.2, turns the global event into a one-dimensional small-ball event
   |<v, X>| <= eps chi(X). Here v spans ker M*, and chi^2(X) = sum <v_i, X>^2 / sigma_i^2 is a correction from the
   other singular directions.
2. **Truncation.** chi is replaced by its sqrt(log n) most critical terms.
3. **Gaussian replacement at scale o(eps).** X is replaced by a Gaussian vector even though eps can be exponentially
   small. This is Fourier analysis on a smooth bump function:
   - low frequencies are handled by Lindeberg's third-moment bound, using that all relevant vectors are delocalized;
   - high frequencies are killed by the LCD of v.
4. **Rescaling.** For Gaussian Z the density of W_0 = <v, Z> is flat at 0, so P(|W_0| <= eps f) = (eps/eps_0)
   P(|W_0| <= eps_0 f). This moves the problem up to a scale eps_0 = n^(-c), where Tao-Vu universality applies.

The lesson for us: a quantity that depends only on a density at a threshold is universal down to tiny scales, provided
the relevant directions are delocalized and arithmetically unstructured.

### 2c. Free-probability derandomization (2601.08111)

- **The free model.** The moments of X_free = sum_i A_i (x) s_i (free semicirculars) are computable by a non-crossing
  recursion. The corresponding classical moments are not.
- **The interpolation.** For any x, put A_t(x) = A_0 (x) 1 + A(x) (x) 1 + sqrt(1 - t) Xbar_free and a potential
  Phi(t, x), for example a trace moment of A_t(x). A Brownian path x_t ~ N(0, t I) runs from 0 to g.
  - By Ito's formula, E Phi(1, x_1) - Phi(0, 0) = int E[d Phi].
  - The drift is the classical second-order increment in x minus the free increment, which comes from splitting
    sqrt(1 - t) Xbar_free into free pieces. The two differ only by crossing terms, bounded by the matrix alignment
    sigma^(1/2) nu^(1/2) ("intrinsic freeness").
- **Derandomization.** The step only uses second-order statistics of dx_t, so pairwise-independent increments suffice.
  A deterministic walk then keeps Phi from increasing. This is the method of conditional expectations, with the free
  model in the role of the conditional expectation. Linear constraints are handled by Lovett-Meka sticky walks.

### 2d. The other chat's ideas (checked here)

- **KMS half-tilt: correct, and exact.** Take a faithful state with vector xi, a self-adjoint x, and mu_x, the spectral
  measure of log Delta in x xi.
  - S x xi = x xi gives Delta^(1/2) x xi = J x xi.
  - J log Delta J = -log Delta then gives e^s mu_x(ds) = mu_x(-ds).
  - So e^(s/2) mu_x is symmetric, it is the spectral measure in Delta^(1/4) x xi, and its mass is
    <x xi, Delta^(1/2) x xi>.

  For a product state on (x) M_2 with the global flip X = (x) sigma^x, the matrix units |b-bar><b| carry eigenvalue
  sum_i +-nu_i, with nu_i = log(p_i / (1 - p_i)), and half-tilted weight prod_i sqrt(p_i (1 - p_i)). That weight is the
  same for every b. So the normalized measure is the Rademacher law of sum eps_i nu_i, with mass
  prod 2 sqrt(p_i (1 - p_i)) = prod sech(nu_i / 2). (Re-derived by hand above.)
- **Klartag's Lemma 3.3 as Ito for -log det.** For a martingale dA,
  d(-log det A) = -tr(A^-1 dA) + (1/2) tr((A^-1 dA)^2),
  so the volume grows by the Ito term alone. That term is the trace of the free-face projection in the metric
  tr(A^-1 X A^-1 Y). (This agrees with the survey's description.)
- **Out of scope here.** The Coxeter-Lehmer items (the Hecke-Kazhdan formula, the Lehmer gap of rigidity windows, the
  growth-rate inequality) are mathematics with no route to the estimator. They are not pursued in this note.

## 3. The common skeleton, and what it says about the network

### 3a. Proposition 1: the free interpolation is our localization

Let e(m, tau) be any estimator run on the input law N(m, tau I) and exact on point masses, e(x, 0) = F(x) (the
production chain is one; note XLII section 3). Put Phi(t, x) = e(x, (1 - t) I) and let x_t be Brownian motion from 0.
Given x_t the input has law N(x_t, (1 - t) I), so Phi(t, x_t) is an approximate conditional expectation of F(g). Ito's
formula gives

    E F(g) - e(0, I) = int_0^1 E[(d_t + (1/2) Lap_x) Phi(t, x_t)] dt = - int_0^1 E[D(x_t, 1 - t)] dt,
    D = d_tau e - (1/2) Lap_m e,

which is note XLII's identity. In Wang-Lau-Zhou the free model plays e's part, and the drift they bound is this heat
defect. So the heat-defect programme is the free method of conditional expectations, with the chain as the free model.

The dictionary:
- **Non-crossing diagrams** are the chain's tree diagrams: births with open legs, transported by mean gates.
- **Crossings** are its loops: the closed walks of note XLII's Theorem A5 and the gate-covariance loops of note XLIII.
- **Intrinsic freeness** is the power counting that suppresses loops in the bulk, where correlations are O(n^(-1/2)).
  The collective modes break that counting (note XLIII 4c: the residual is about 10x the naive bulk order).

Pairwise independence (the drift needs only second-order statistics of the increment) is the statement that D depends
only on e's Hessian.

### 3b. Proposition 2: the means are a gated transport of threshold sources (exact)

For every unit, E relu(h) = E[h 1(h > 0)] = E h P(h > 0) + Cov(h, 1(h > 0)). Write:
- Gamma_l = diag P(h_l > 0), the true gate probabilities;
- s_l = Cov(h_l, 1(h_l > 0)), unit by unit.

Then

    m_l = Gamma_l W_l m_(l-1) + s_l,   m_L = sum_(k=1..L) Gamma_L W_L ... Gamma_(k+1) W_(k+1) s_k   (m_0 = 0).

For h a functional of the Gaussian input, the Nourdin-Peccati identity gives s = p_h(0) tau_h(0): the density at the
threshold (the concentration function at scale zero) times the Stein kernel there. In note XLII's notation this is
also E[x_L] = E[Lap x_L], the kink sum. So the means are the solution of a linear gated transport equation driven by
threshold-local sources. A chain's errors enter as

    delta m_L = sum_k T_(L<-k) (delta s_k + delta Gamma_k W_k m_(k-1)) + second order.

For a Gaussian unit, d E relu / d v = phi(alpha) / (2 s) = p(0) / 2. So a variance error delta v becomes a mean error
(1/2) (threshold density) delta v. The concentration function at the threshold is the lever arm of every second-order
error.

### 3c. Theorem 3: the fresh-weight theorem

**Statement.** Let A be any symmetric matrix that depends only on W_1, ..., W_l: for instance a chain's error in the
post-activation covariance of layer l, or any correction a chain omits there. Let W = W_(l+1) have iid N(0, sigma^2)
entries (He: sigma^2 = 2/n; measured on network 0: n var = 2.005, excess kurtosis -0.003). Put E = W A W^T. Then,
exactly:
- (i) the diagonal entries E_aa = w_a^T A w_a are independent across units given A, with
  E E_aa = sigma^2 tr A and Var E_aa = 2 sigma^4 ||A||_F^2;
- (ii) off the diagonal, E E_ab = 0, E E_ab^2 = sigma^4 ||A||_F^2, and E E_ab E_ac = 0 for b != c;
- (iii) consequently Var_a(E_aa) / mean_(a != b)(E_ab^2) = 2, whatever the rank or spectrum of A.

*Proof.* These are the Gaussian quadratic-form identities Var(w^T A w) = 2 sigma^4 tr A^2 and
E(w^T A w')^2 = sigma^2 E w^T A^2 w = sigma^4 tr A^2 for independent rows. Independence of A from W_(l+1) holds because
A is built from the first l layers.

**Reading 1: the variance injection is fresh-weight noise.** Any estimator that forms the next pre-activation
covariance by the exact linear map makes per-neuron variance errors delta v_(l+1,a) = w_a^T A_l w_a. These are:
- a common shift sigma^2 tr A_l;
- plus independent per-neuron noise of rms sigma^2 sqrt(2) ||A_l||_F.

Only the Frobenius norm of the previous layer's error matters; its shape does not.

**Reading 2: the nibble.** E_aa is a "degree" and E_ab a "codegree" of the error operator A. Their second moments are
tr A^2 = sum_cd A_cd^2, the weight of closed 2-walks, exactly as the nibble's martingale increments are bounded by
codegrees. This is the precise form of the user's codegree / overlap-of-projections bridge: the injected noise at a
layer is controlled by the overlap of the error operator with itself.

**Reading 3: typical Littlewood-Offord and Spielman-Teng.** E_aa - sigma^2 tr A is a decoupled quadratic
Littlewood-Offord sum with random, unstructured coefficients, so it is close to Gaussian once ||A||_op << ||A||_F.
Its mean sigma^2 tr A is the universal part: it is fixed by the coincident-site (two-site) data of A, as Gaussian
replacement fixes the small-ball law. Its fluctuation is the part tied to the specific direction w_a.

**Reading 4: freeness.** W_(l+1) is free from the past, so the predictable part of W A W^T is its trace (note XIII
(ii)). Theorem 3 measures the unpredictable part exactly.

### 3d. Corollaries

1. **No two-site lunch.** Let Omega_l be a class a chain omits at layer l, for example the all-distinct gate
   covariance of note XLIII 3d. Its effect on the next layer's variances is sigma^2 tr Omega_l, which is coincident-site
   data, plus a fresh-weight chaos of rms sigma^2 sqrt(2) ||Omega_l||_F. That chaos has zero mean given the past. It is
   not determined by any rotation-invariant statistic of Omega_l, because Gaussian W_(l+1) is rotation invariant. So
   no closure built from carried slices can recover it: it has to be computed as Omega_l contracted with the actual
   rows. This is why note XLIII 4e's six two-site completions explain only 14-29% of the D21 error.
2. **The accounting.** Injections made through different fresh weights are uncorrelated. This is the measured fact of
   note XLII section 2 ("injections at different layers are uncorrelated"), and here it is a theorem.

   The final squared error is therefore a sum over layers. Each term is a transport factor times
   2 sigma^4 ||A_l||_F^2 (1/2 p(0))^2, plus the mean-shift and mean-error terms of Proposition 2. The design target of
   any second-moment method is this Frobenius budget, layer by layer, not slice-wise accuracy.
3. **Why truth-free detectors go blind late.** A detector built from the chain's own state, such as the heat defect,
   sees Omega_l only through what the chain computes. The injected noise is a chaos of Omega_l against fresh rows,
   which the chain never forms. This fits note XLII's X* = 0.11 at the final layer. It is an explanation, not a
   theorem: a detector that re-runs the chain on perturbed inputs could in principle encode more.

### 3e. Where the other chat's modular ideas land

- **The one-parameter group.** The natural one-parameter group of this problem is the dilation group: it commutes
  with relu and preserves the quasi-free family (note XIII (iii)). Its "modular Hamiltonian" is the log-gain. The
  radial disintegration is the KMS-like decomposition along it, and note XIII E4 measured that it removes the
  coherent, direction-averaged non-Gaussianity of generic directions.
- **The half-tilt theorem.** It applies to the gate pattern under a product law. There the half-tilted spectral
  measure of the global gate flip is the Rademacher law of sum_c eps_c nu_c, with nu_c the gate log-odds.
- **Where it stops.** Proposition 2 shows that the means need the gate probabilities and the threshold sources, not
  flip affinities. So the sech / Chernoff structure has no foothold in the mean problem. This is recorded as a
  negative scope statement.

## 4. What it says about the system

By Theorem 3 the late floor is fresh-weight noise from whatever the chain omits. Its size is set by Frobenius norms.
Removing it requires computing the omitted class against the actual rows, at a cost proportional to the class's rank.
That leaves exactly two structural routes that could beat the moment chain's elasticity-1 frontier (note XLII F3,
Theorem 9). There is also a fallback.

**Route R: dilation fibres (radial KMS).** Suppose the late law is a gain mixture of quasi-free fibres, h = sqrt(G) y
with y ~ N(mu~, S~) independent of G. Then every two-site slice is rank-structured in (mu, C) with three gain moments
g_k = E G^(k/2). For the third cumulant, writing v~_a = S~_aa:

    k(a,a,b) = (g3 - g1 g2) [(v~_a + mu~_a^2) mu~_b + 2 S~_ab mu~_a] - 2 g1 (g2 - g1^2) mu~_a^2 mu~_b.

This shape matches both features of the D21 error found in note XLIII 4d:
- a C_ab x (function of a) part, like the regression form it correlated with at +0.57;
- a column part (v_a + mu_a^2) mu_b, which is a coupling of neuron b to the layer's energy (20% of the error lives in
  column means).

None of the three forms was among the six completions tried in note XLIII 4e. If they carry the error, the fix costs
O(n^2) plus three scalars per layer, and the three-site transport they stand in for is the dilation mode in disguise.

**Route C: an exact collective block.** Suppose the post-activation covariance error A_l is concentrated in the top-k
collective modes. Then the fresh-weight noise is sum_j dlambda_j (w_a . u_j)^2 plus cross terms. A k-dimensional
latent computed exactly would remove it:
- particles in the latent cost O(P n k) a layer;
- the bulk stays Gaussian at the cost of the covariance transport, which is 3% of the bill;
- unlike F3's leaves, the n^3 work is shared.

**The fallback: allocate cost by the Frobenius budget.** If neither route holds, the late floor is bulk fresh-weight
noise from three-site structure, and it costs one leg transport per unit of rank. The remaining lever on the adjusted
score is then to spend FLOPs where the transport factor times the Frobenius contribution per FLOP is largest. Layers
0-8 inject only 13% of the final squared error (note XLII section 2).

## 5. Pre-registration (committed before the data)

**Data** (cached, networks 0 and 1):
- the production chain's pre-activation state per layer (`chaindump2_{0,1}.npz`: var, C_off, D21, D3);
- Monte Carlo truth at 1.6e7 inputs in two independent halves (`mc2_off{0,1}_{full,h0,h1}.npz`: mu, var, cov, D21).

Every second moment of an error is estimated noise-free, as the product of the chain-minus-truth differences against
the two halves. Script: `code/fw_test.py`.

**T1: the fresh-weight theorem on the chain's own error.**
- *Quantity.* E_s = C_chain(s) - C_true(s) (pre-activation covariance, layer s) and
  R_s = Var_a(E_aa) / mean_(a != b)(E_ab^2), at layers 2-14. E_s = W_s A_(s-1) W_s^T with A built from earlier layers,
  so the theorem predicts R = 2 exactly. It can fail only if the chain's pre-activation stage adds non-transport
  terms, such as the per-unit calibrations or look-ahead counterterms.
- *Prediction.* R in [1.7, 2.3] at >= 10 of the 13 layers on both networks; central 1.95.

**T2: route R.**
- *Quantity.* Noise-free regression of the D21 error dD = D21_true - D21_chain (off the diagonal, layers 6-14) on
  F1 = (v_a + mu_a^2) mu_b, F2 = C_ab mu_a and F3 = mu_a^2 mu_b, built from the true state.
  - The metric is weighted by omega_ab = p_a(0) P(h_b > 0), the leading Mehler weight with which D21 enters the
    post-activation covariance.
  - The report is the share of the weighted noise-free energy explained, with all three fitted per layer in sample.
- *Prediction.* 0.18 (0.10-0.30) at layers 9-14. The births' fresh-weight structure should dominate, so the dilation
  shape should not.

**T3: route C.**
- *Quantity.* Let P_k project on the top-k eigenvectors of C_true(s). The collective share of the error is
  c_k = 1 - ||P_k^perp E_s P_k^perp||_F^2 / ||E_s||_F^2 (off-diagonal, noise-free), for k = 8, 32, 128 at layers
  12-14. A random subspace gives 1 - (1 - k/n)^2: 0.016, 0.061 and 0.234.
- *Prediction.* c_32 = 0.35 (0.20-0.50). For context the same share is also reported for the true non-Gaussian part
  of the covariance, C_true(s) - W_s KG(true state at s - 1) W_s^T. No prediction is registered for it.

**Decision rule.**
- **T1.** If it holds, the Frobenius budget is adopted as the error accounting, and section 3d's corollaries hold for
  the production chain. If it fails, record the layers and look for non-transport terms in the pre-activation stage
  before using the corollaries.
- **Route R.** A T2 share >= 0.50 at layers 9-14 on both networks means building the dilation-fibred readout of the
  (2,1) slice: three gain moments per layer, O(n^2). It is then cold-screened on 16 networks against production. A
  share < 0.30 closes route R.
- **Route C.** c_32 >= 0.60 at layers 12-14 on both networks means building the collective latent block, then a cold
  screen. A value < 0.40 closes route C.
- **Both closed.** The next step is the Frobenius budget map: transport factor times ||A_l||_F^2 contribution against
  FLOPs, by layer and component. It is used to re-allocate the bill, and that re-allocation is pre-registered
  separately.

## 6. Results (networks 0 and 1; `outputs/fw_test_net{0,1}.txt`)

**T1: the fresh-weight theorem holds.** R_s = Var_a(E_aa) / mean E_ab^2 (noise-free):

| layer | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| net 0 | 1.06 | 3.00 | 1.06 | 2.29 | 1.99 | 2.09 | 2.05 | 2.13 | 1.90 | 2.25 | 2.06 | 2.01 | 2.04 |
| net 1 | 4.69 | 2.62 | 1.84 | 2.10 | 1.92 | 1.90 | 2.19 | 2.02 | 1.90 | 2.12 | 1.90 | 2.08 | 1.94 |

- R falls in [1.7, 2.3] at 10 of 13 layers on network 0 and 11 of 13 on network 1, so it passes as registered.
- From layer 5 on R = 2.03 +- 0.11, against the theorem's 2.
- The failures are layers 2-4, where the error is tiny (rms dv/v 5e-5 to 2.4e-4) and non-transport terms dominate.
- The common shift sigma^2 tr A is negative at every layer from 3 on: -0.29 to -0.66 of the per-neuron rms. The chain
  underestimates variances uniformly by about half of its scatter.

So from layer 5 on, the production chain's per-neuron variance error is the theorem's object: fresh-weight noise set
by ||A||_F, plus a common shift.

**T2: route R is neither built nor closed.** Share of the omega-weighted D21 error explained by the three
dilation-fibre forms, cross-fitted:

| layer | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|---|
| net 0 | 0.24 | 0.27 | 0.36 | 0.34 | 0.33 | 0.38 | 0.39 | 0.29 | 0.29 |
| net 1 | 0.29 | 0.29 | 0.36 | 0.26 | 0.32 | 0.33 | 0.34 | 0.24 | 0.29 |

- Almost all of it is F1 = (v_a + mu_a^2) mu_b alone (0.20-0.34), the coupling of neuron b to the layer's energy that
  note XLIII 4d saw as "column means". F2 explains nothing and F3 adds about 0.1.
- The prediction (0.18, 0.10-0.30) is exceeded. The build threshold (0.50) is not reached, and the closing threshold
  (< 0.30) is not met at most layers. Recorded; no build.

**T3: the error is collective, so route C is triggered.** c_k is the collective share of the off-diagonal covariance
error; the random baselines are 0.016 / 0.061 / 0.234.

| layer | net 0 error c_8 / c_32 / c_128 | net 1 error | true non-Gaussian part c_8 / c_32 / c_128 (net 0; net 1) | ||g|| / ||E|| |
|---|---|---|---|---|
| 12 | 0.64 / 0.83 / 0.97 | 0.68 / 0.85 / 0.97 | 0.87 / 0.94 / 0.98; 0.86 / 0.93 / 0.98 | 7.9; 6.9 |
| 13 | 0.70 / 0.86 / 0.98 | 0.66 / 0.85 / 0.98 | 0.88 / 0.94 / 0.99; 0.86 / 0.93 / 0.98 | 6.7; 6.5 |
| 14 | 0.69 / 0.87 / 0.98 | 0.72 / 0.88 / 0.98 | 0.88 / 0.95 / 0.99; 0.88 / 0.94 / 0.99 | 6.3; 6.3 |

- c_32 = 0.83-0.88 on both networks, far above the prediction (0.35, 0.20-0.50) and above the 0.60 build threshold.
- The true non-Gaussian part of the covariance is even more collective: 0.93-0.95 in 32 modes and 0.86-0.88 in 8.
- The chain removes all but about 1/7 of it in Frobenius norm.

**What it means.** Together with T1:
- 85% of the energy of the next layer's per-neuron variance error is sum_(jk) M_jk (w_a . u_j)(w_a . u_k), where u_j
  are the top-32 eigenvectors of the true pre-activation covariance and M is a 32 x 32 symmetric block;
- only the k x k block M is unknown;
- the bulk carries 12-17%.

This is the opposite of note XLII E0's reading that the final error is "the incoherent bulk". E0 measured the final
mean error. In the metric that drives the injection, the Frobenius norm, the error is collective. So the late floor is
a low-dimensional object, which the fresh-weight theorem makes exactly readable: what has to be computed better is the
covariance restricted to about 32 modes.

## 7. What the collective error is (exploratory, not pre-registered)

### 7a. One mode and one number (`code/fw_explore.py`, `outputs/fw_explore_net{0,1}.txt`)

In the eigenbasis of the true pre-activation covariance at layers 8-14:

| quantity | network 0 | network 1 |
|---|---|---|
| c_1: share of the off-diagonal error energy in the top mode alone | 0.33 (layer 8) -> 0.58 (layer 13) | 0.35 -> 0.59 |
| c_8 | 0.45 -> 0.70 | 0.51 -> 0.72 |
| eigenvalue (diagonal) share of the top-8 error block | 0.76-0.94 | 0.78-0.94 |
| top-mode relative eigenvalue error dlam_1 / lam_1 | -4.4e-3 to -6.7e-3 | -4.8e-3 to -8.3e-3 |
| other top-8 modes, dlam_j / lam_j | +-2e-3, mixed signs | the same |
| lam_1 / mean eigenvalue | 35 -> 88 | 35 -> 84 |
| propagated share of the error energy (Gaussian closure of the chain's state vs the true state) | 0.56 -> 0.80 | 0.61 -> 0.82 |

The off-diagonal block excludes the variance diagonal, so dlam_1 here is the off-diagonal contribution.

- **The largest structured component of the late second-moment error is one number per layer.** It is a 0.5-0.8%
  underestimate of the top eigenvalue of the pre-activation covariance, negative at every layer on both networks and
  growing with depth.
- **The top mode is a large spike** (35-88 times the mean eigenvalue).
- **It is mostly inherited**: the propagated share rises to 0.8 by layer 14, so it accumulates.
- **What it explains.** By Theorem 3, a rank-one error dlam u u^T in the post-activation covariance enters the next
  layer's variances as dlam (w_a . u)^2. That has a common part (part of the negative shift in T1) and a part tied to
  the known per-neuron feature (w_a . u)^2.

### 7b. The collective error is the dilation mode (`code/fw_explore2.py`, `code/tau_check.py`)

**The top mode is the mean direction.** |u_1 . mu_hat| = 0.93-0.98 at layers 8-14 on both networks. The one-parameter
group behind it is the dilation group of section 3e: the top mode is the gain, the mode its "modular Hamiltonian"
log G moves.

**Its error budget**, relative to lam_1 (layers 8-14, both networks):
- total -0.45% to -0.83%;
- inherited (through the Gaussian closure of the chain's own earlier state) -0.37% to -0.79%;
- injected per layer -0.02% to -0.14%.

**The injection by first-order Edgeworth class** (the chain's slices at its own state minus the true slices, projected
on r = W_s^T u_1):
- D21 is negative at every layer: -0.06% to -0.17%;
- K31 is positive at every layer and partly compensates: +0.03% to +0.06%;
- K22, kappa3 and kappa4 are below 0.01%;
- the class sum reproduces the measured injection to within a factor of about 2.

The (2,1) slice's projection on the dilation mode is T2's F1 column shape (the coupling to the layer's energy).
Route R's partial signal and route C's collective signal are therefore one object: the chain's D21 under-feeds the
dilation mode's variance.

**Size against the multiplicative heuristic.** The relative variance of the mode, tau = lam_1 / (u . mu)^2, is
stationary from layer 6 on (0.0087-0.0100 on network 0, 0.0077-0.0094 on network 1). A gain that composes
multiplicatively through independent layers, but is carried additively, would leave a relative deficit of about
tau/2. The measured deficit is 0.7-1.4 times tau/2 on network 0 and 0.8-2.1 times on network 1 (layers 6-14). The ratio
grows with depth. The heuristic gets the order and the sign. It is not a law: with stationary tau the gain is
mean-reverting, not a product of independent layer gains. So no parameter-free correction follows from it yet.

**Where this leaves the design.**
- The fresh-weight theorem reduces the late floor to Frobenius noise of the covariance error.
- That error is collective, and its largest single piece is one number per layer, the dilation mode's eigenvalue,
  biased by -0.5% to -0.8%.
- Its effect on the next layer is that number times the known per-neuron feature (w_a . u)^2; on the means it is that
  again times p_a(0) / 2 (Proposition 2).

What is still missing is an exact law for the dilation sector: how log G's increments compose through a relu layer
when the gain is mean-reverting. That is the next theory item. Its value is bounded by a dilation-mode oracle (the true
eigenvalue imposed along the true mean direction at every layer). Both are to be pre-registered before they are run.

### 7c. Why it collapses onto that mode: the dilation mode is critical

From the same outputs, take the transport factor of the top mode's relative error:

    a_s = (propagated part at layer s) / (total at layer s - 1).

| layer | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|
| net 0 | 0.80 | 0.90 | 0.92 | 0.97 | 1.02 | 0.84 |
| net 1 | 0.83 | 0.93 | 1.02 | 1.12 | 0.92 | 1.14 |

**The mode is marginal.** a is about 0.95 on average, so the relative error is carried nearly undamped. The total at
layer 14 is 8.6 times (network 0) and 16 times (network 1) the mean per-layer injection (-7.0e-4 and -5.2e-4 of
lam_1). Every other direction contracts: an injection at layer 3 loses about 9x in squared norm by layer 15 (note
XLII section 2).

This is note XIII (i) made quantitative. He scaling makes the gated transport a unital channel, and the gain is its
critical zero mode. So the chain's late floor is not where the per-layer injection is largest. It is where injections
are integrated without damping.

**The design rule this gives.** Allocate accuracy by inverse spectral gap. An injection error into a mode with
transport factor a costs about 1/(1 - a) layers of accumulation, so the critical mode needs about ten times the
per-layer accuracy of the rest. It is one direction per layer, so exactness there is cheap in principle:
- O(n^2) for the projections;
- O(n^3) per layer for one extra collective leg (section 7b's loop through the mode).

The dilation-sector item of 7b is therefore an exact (conservation-respecting) treatment of one critical mode, not a
collective block of rank 32.

## 8. Stage 9o (other chat) read against these measurements

The user supplied "Ellipsoids, codes and contacts" (stage 9o: Klartag's process, centroid bodies, billiard rigidity).
Each claim was checked by hand.

**Correct as stated.**
- **(I) The ellipsoid averaging principle.** E_N(0,Sigma) y = ((n + d)/n) E r^d <y>_(E_Sigma). This is polar
  coordinates: Gaussian and uniform-ellipsoid laws share the angular part.
- **(II) The centroid-body recursion.** E a_(l+1,k) = h_(Gamma+(mu_l))(w_k) is the definition of E relu<w, x>. The
  Edgeworth coefficient -K3(w,w,w) mu phi / (6 sigma^3) is right: d^3/dmu^3 E relu = -mu phi / sigma^3.
- **(III) The Ito-Price drift.** It is (1/8) E <grad^4 y, Pi_F>. For 1-homogeneous y, E Lap^2 y = (d - 2) d E y = -E y,
  so the unconstrained drift is -(1/8) E y.
- **(IV) The code-moment identity.** For layer 2 with isotropic, zero-mean input,
  Var z_(2,k) = sum_r b_r ||M_r(nu_k)||^2.
- **(V) Exact closure on an open set forces oddness.** On a ray, z = p|xi| + q xi has
  kappa_3 = p (p^2 kappa_3(|xi|) + 3 q^2 sqrt(2/pi)), which vanishes only at p = 0.
- **The token ladder.** (1 - r) f_r = (r + 1)(r + 2) tr f_(r+2) follows from the Euler operator in Fock form,
  x . grad = N + sum_i a_i^2. Its rung r = 0 is E y = 2 tr f_2 = E Lap y, the kink sum.

**Two corrections from the measurements.**
1. **The "design defect" pairing is one channel, not the leading correction.** <K3, sum_k c_k w_k^(x)3> is the direct
   per-neuron mean channel, and the chain already computes it exactly (its D3). The first-order error also flows
   through the next layer's covariance. That channel pairs K3 with sum_(a,b) c_ab w_a (x) w_a (x) w_b, of CP rank up
   to n^2, which is the Khatri-Rao wall. Here it dominates: kappa(a,a,b) carries 70-80% of the variance injection and
   the per-neuron kappa3 1-2% (note XLIII 4d); the late D3 oracle is worth -2.5% / -13%.
2. **The ellipsoid / radial factor is the trivial part of the gain.** The input radius contributes tau = 1/(2n) =
   4.9e-4 to the mean direction's relative variance. The measured critical mode has tau = 0.008-0.010, about 18 times
   larger: it is the network's own angular amplification, which no input-side disintegration separates.

**What is usable.**
- **The ladder is an exact, truth-free set of constraints between chaos traces.** Rung 1 (tr f_3 = 0) and rung 0
  (E y = 2 tr f_2) can be checked on the chain's chaos-graded pieces: the legs (second chaos, built from first-chaos
  vectors only) and V56's third-chaos star.
- **The ladder measures what the legs omit.** Rung 0 makes "chain mean minus 2 tr H_chain" exactly the trace of the
  second chaos the legs omit, unit by unit. That is a free measurement of the omitted birth kernels (the gradient
  fluctuation at each kink). Whether it predicts the chain's error is untested.
- **The frozen-contact Klartag drift (III ii) is the rung-2 defect** of note XLII's ladder, with the first-layer terms
  removed. It is consistent with Theorem A5, and it has no cheap evaluation.

## 9. Stage 13 (other chat's note, 21 pages) tested on our data

The user attached the stage-13 note (`Localization in the commutant, stage 13`) and asked for real engagement with its
novel claims. I read all 21 pages. The note makes four claims that can be tested on data we hold:

1. **Eldan's cumulant drift is the Schur hub** (Prop. 4.1-4.3, Conjecture 4.5): the hub is the between-tier fourth
   cumulant, 3 Var(proj_field(z^2)). Verified by hand here for the diagonal (Cov(z^2, x) = 2HL, so 3 Var(proj) = 12|HL|^2,
   which is note XLI's path class). The (2,2) and (3,1) pairings need the three-site third cumulant contracted once
   against C^(-1) D_i, which is one hub-sized contraction per source-layer: about +40% of the bill. Not pursued as a build.
2. **The dilation charge is conserved and the chain's pieces are mutually inconsistent along it** (Prop. 5.1, Thm 5.4,
   Cor. 5.5). Checked by hand: for z = (1 + delta)(mu + y), Var delta = eps, the first-order cumulants are
   - C: eps (Ct + mu mu^T);
   - kappa3(a,b,c): 2 eps (mu_a Ct_bc + mu_b Ct_ac + mu_c Ct_ab);
   - kappa4(a,a,a,a) = 12 eps Ct_aa^2; (2,2): eps (4 Ct_aa Ct_bb + 8 Ct_ab^2); (3,1): 12 eps Ct_aa Ct_ab.

   The post-ReLU mean E relu(z) = E(1+delta) E relu(g) does not depend on eps at all. A Gaussian closure sees a spurious
   O(eps) through the variance, and the Edgeworth terms cancel it exactly when they carry the mixture's own skewness
   6 eps r and kurtosis 12 eps (the note's cancellation theorem; I re-derived it, and the pieces sum to (eps/2)(1 + r^2) -
   eps r^2 + (eps/2)(r^2 - 1) = 0). It also reproduces the production chain's fitted scale-mixture kappa4 (K4SM:
   k4 = 3 g var^2, K22 = g var var^T, K31 = 3 g d(var) C_off with g = 4 eps), so our fitted "pair" is the kappa4 piece
   of this package with a free gain.
3. **The gain saturates, it does not mean-revert** (Prop. 5.2, Remark 5.3). This contradicts the reading I wrote in
   section 7b ("mean-reverting"). It is falsifiable (P1 below), and I will correct 7b if it holds.
4. **Truth-free exact identities** (ladder rungs 0 and 1, memory budget): not tested here.

**Why claim 2 matters for us.** Earlier rounds saw three facts that look unrelated:
- replacing the chain's kappa4 diagonal by the Monte Carlo value at layers 12-15 makes the last step worse
  (+16% / +34%, note XLIII section 4a; 4.18e-9 -> 7.29e-9 in `frontier-tests`);
- the output-fitted counterterms push D21 down by 1-2% where the slice-level truth wants it 2-4% up (note XLIII 4e);
- a gain-conditioned chain double-counts the radial zero mode (`arrow-mixture`).

All three are what a first-order cancellation among (variance, kappa3, kappa4) along one direction would produce: fix
one piece and the compensation breaks. The tests below measure that directly.

### 9a. Pre-registration (committed before the data)

**Script `code/dil_ledger.py`, networks 0-1, layers 3-14.** Part (A) fits each slice's dilation charge eps_S (the
least-squares coefficient on the package form), for the truth and for the chain's own slices. Part (B) takes the
one-step first-order Edgeworth sensitivities at the true state, dm = phi/(2S) dvar - R phi/(6 S^2) dk3 + (R^2 - 1)
phi/(24 S^3) dk4 + Phi dmu, and evaluates them with the chain's error in each slice (noise-free through the two
independent halves).

Predictions:
- **A1 (truth).** eps_k3 is 0.003-0.008 at layers 8-14 (the note's 0.71 V/4); eps_k4 / eps_k3 is 0.7-1.0 (the note's
  kurtosis-to-skewness slope ratio 1.6-1.8, halved); D21, K22 and K31 give charges within a factor of 2 of eps_k3.
- **A2 (chain).** The chain's Rayleigh charge is within 1.5% of the truth's. Its eps_k3 is within 10%. Its eps_k4 is more
  than 20% off (the memoryless regeneration).
- **B1 (compensation).** At 8 or more of the 12 layers on both networks, corr(var-term, k4-term) < -0.3 and
  rms(sum of the three terms) < 0.8 x their quadrature sum.
- **B2.** The three terms individually have rms within a factor of 3 of the actual one-step mean error (the leak is
  not negligible).

Decision rule: if B1 fails (ratio >= 0.9 and correlations within +-0.15), the compensation story is closed and the
consistency route stops. If B1 holds, the next step is an oracle that imposes the package on the chain's own state, with
and without the counterterms.

**P1 (the note's falsifiable prediction), script `code/lognorm.py`, network 0, N = 40000 Gaussian inputs, float32.** With
Delta_l = log(|h_l|^2 / |h_(l-1)|^2) across inputs:
- the correlation of successive increments lies within +-0.15 at all layers from 4 on, and the AR(1) slope of Delta_(l+1)
  on log|h_l|^2 lies within +-0.1 (mean reversion would make it clearly negative);
- Var Delta_l is 0.8-1.6 times the leading law 4 theta_(l-1)^2 / n at layers 6-15 (theta from the arc-cosine recursion);
- V_l / 4 (V_l = Var_x log|h_l|^2) is 0.0065-0.0090 at layers 8, 12, 15.

### 9b. Results of the log-norm test and the ledger (read after the pre-registration above)

**P1 holds on network 0** (`outputs/lognorm_net0.txt`, N = 40000):
- successive log-norm increments: correlation within +-0.14 at every layer, and no sign preference (mean +0.01);
- AR(1) slope of the next increment on the current log-norm: |b| < 0.025 at every layer (mean reversion would be clearly negative);
- Var Delta_l is 0.84-1.19 times the leading law 4 theta_(l-1)^2 / n at layers 6-16;
- V_l / 4 = 0.0060, 0.0077, 0.0083 at layers 8, 12, 15 (registered 0.0065-0.0090; layer 8 is 8% below the interval's floor), and
  flat from layer 13 (0.0081, 0.0082, 0.0083, 0.0081).

**So the gain saturates; it does not mean-revert. My reading in section 7b was wrong** and is superseded by this section. The
stationary relative variance seen in section 7 is a saturating sum of uncorrelated increments (summable because two inputs
at angle theta disagree on a fraction theta/pi of gates), not a restoring force.

**Ledger, part (A)** (`outputs/dil_ledger_net{0,1}.txt`; charges x 1e-3, layer 14, truth | chain):

| | Rayleigh | k3 | D21 | K22 | K31 | k4 |
|---|---|---|---|---|---|---|
| net 0 | 8.56 \| 8.49 | 5.53 \| 5.27 | 5.48 \| 3.81 | 5.24 \| 4.30 | 5.53 \| 2.35 | 5.27 \| 4.47 |
| net 1 | 7.34 \| 7.27 | 5.01 \| 4.75 | 5.01 \| 3.48 | 5.05 \| 4.16 | 5.12 \| 2.30 | 5.01 \| 4.22 |

- **The truth is a pure common-scale mixture at deep layers.** At layers 12-14 the five slices give one charge, eps = 5.0-5.5e-3,
  to within 1-4% (net 1, layer 14: 5.01, 5.01, 5.05, 5.12, 5.01). The agreement is the cancellation structure of the stage-13
  note (and of our own gain-mode null tuple of note XXXVI) showing up in five independent Monte Carlo slices. At shallow
  layers kappa4's charge is lower than kappa3's (ratio 0.65 at layer 3, 0.85 at 10, 0.95 at 14), as in the stage-13 note.
- **The mixture explains 65-68% of the gain variance**: V_14 / 4 = 0.0082 against eps = 0.0054 (the note's 0.71 and 0.57-0.63).
- **The chain's slices are mutually inconsistent along the dilation**, and strongly so. At layer 14 it carries, as a fraction
  of the truth's charge: Rayleigh 99% / 99%, k3 95% / 95%, k4 85% / 84%, K22 82% / 82%, D21 70% / 69%, K31 42% / 45%
  (nets 0 / 1). The D21 and K31 deficits grow with depth (D21 0.89 at layer 5 -> 0.70 at layer 14; K31 0.66 -> 0.42).

**Ledger, part (B)** (one-step mean-error terms at the true state, rms x 1e-5, layer 14, nets 0 / 1): variance term 4.9 / 5.2,
k3 term 3.9 / 4.2, k4 term 2.7 / 2.7, propagated-mean term 10.2 / 9.1, actual one-step error 11.3 / 11.9. The sum of the three
slice terms is 5.9 / 6.3 (0.86 / 0.88 of their quadrature sum) with pairwise correlations between -0.17 and -0.09.
**B1 fails as registered**: ratios 0.80-0.96, correlations within +-0.25 (the criterion was < -0.3 at 8 of 12 layers and a ratio
< 0.8). The slice errors compensate each other only weakly at one step. **B2 holds**: each slice term is 2.3-4.3 times smaller
than the actual one-step error, so the three slice errors together carry about 27% of its energy; the rest is the inherited
mean error. Consequence: whatever the dilation inconsistency does, it does not do it through compensation among the three
terms at the same layer. What it does through the transport (which the one-step ledger cannot see) is what the oracle below
measures. This agrees with note XXXVI section 3k(3): the per-unit mean error implied by the inconsistency is not aligned with
the chain's one-step query error (cosines -0.08 to +0.05 there).

**What is already known and what is new.** The cancellation theorem (Cor. 5.5) is our note XXXVI null tuple: the exact
direction (t (Sigma + mu mu^T), 6 t mu Sigma, 12 t Sigma Sigma), invisible to every future mean, with the earlier amplitude fits
(kappa4 at 83-85%, (3,1) at 54-57% of the truth). Our earlier tests changed single amplitudes at a time, always at the chain's
unchanged mean amplitude, and each made the output worse (3 g var^2 +53%, x1.19 +58%; note XXXVI 3k). **Not tested before: all
slices' gain amplitudes set jointly to the truth's, including the D21 deficit (30%) that the earlier templates did not see**
(they fitted t3 on D21 at 99%). That joint correction is O(n^2) per layer.

### 9c. Pre-registration of the gain-package oracle (committed before the runs)

Research switch `V62_GPK` in the research copy of the production estimator (`workbench/k3work/estimator_final_v56.py`, off by
default): after the counterterms, at layers 3-14, each chosen slice is shifted along its package form so that its charge equals
the target from `code/gpk_targets.py` (truth); the pre-injection D3 / D21 stay as what the legs carry (the V43 convention), so the
injected content flows into the newborn source. Networks 0 and 1, counterterms on, float32, baseline raw 1.4771e-8 (network 0).

Runs (raw final-layer MSE against the dataset truth, change from the same network's baseline):
- G1: all six slices (variance Rayleigh, k3, D21, K22, K31, k4);
- G2: D21 + K31 only;
- G3: k4 diagonal + K22 only;
- G4: the variance (Rayleigh) only.

Predictions: G1 -10% to -30% (central -18%) on both networks; G2 -8% to -25% (central -15%); G3 between -3% and +10%; G4 between
-5% and +15%. The compensation story predicts that partial corrections (G3, G4) do not help while the joint correction (G1)
does.

Decision rule: G1 <= -10% on both networks -> build the truth-free version (every charge set to the chain's own k3 charge at
deep layers, shapes below layer 10 from the derived ratio profile) and run a cold 16-network screen; -3% to -10% -> record the
ladder and examine which slices matter; > -3% or worse -> the gain-shaped amplitude is not the lever and this line closes.

### 9d. Result of the gain-package oracle (network 0; `outputs/gpk_ladder_net0.txt`): the prediction fails by a wide margin

Baseline raw 1.4771e-8 (cold, local, float32). Every variant that moves the gain-shaped part of D21, K22, K31 or k4 toward the
truth's charge makes the output worse:

| variant (layers 3-14, counterterms on) | raw | change |
|---|---|---|
| G1 all six slices, readout only | 6.239e-8 | +322% |
| G1 all six slices, persisted into the newborn source | 1.433e-7 | +870% |
| G2 D21 + K31, readout only | 5.594e-8 | +279% |
| G2 D21 + K31, persisted | 1.824e-7 | +1135% |
| D21 alone | 5.762e-8 | +290% |
| K31 alone | 1.612e-8 | +9.1% |
| G3 k4 + K22 | 2.314e-8 | +57% |
| G4 variance (Rayleigh) only | 1.4757e-8 | -0.1% |

**Against the registration**: G1 predicted -10% to -30%, observed +322%; G2 predicted -8% to -25%, observed +279%; G3
predicted between -3% and +10%, observed +57%; G4 predicted between -5% and +15%, observed -0.1% (inside). **Decision rule:
G1 worse than -3% -> the gain-shaped amplitude is not the lever and this line closes.** (Network 1 was not run: the
network-0 outcome is far outside any plausible network-to-network spread.)

**Why** (one measurement, `outputs/err_along_mean.txt`): the chain's final-layer error vector has only 0.4-3.7% of its energy
along the mean direction (layers 3-15, both networks), and its neuron-averaged signed error is 2-12% of its rms. The error that
remains is incoherent per-neuron scatter. The coherent (mean-direction) sector is exactly the one the output-fitted counterterms
(the V47 table: per-layer rescales of var, D3, D21, g4, k22, k31, coff) have already tuned out, and the output is extremely
sensitive to it: moving D21's gain amplitude alone by +30% multiplies the raw MSE by 3.9. The chain's deficits in the five slices
(D21 70%, K31 42%, K22 82%, k4 85% of the truth's charge) are therefore not an error to remove; together with its other errors
they are the compensated optimum that the counterterm fit found. Setting them to the truth's values destroys the compensation.
The variance-only correction (G4), which has no compensation partner, is neutral because the chain's own Rayleigh quotient is
already within 1% of the truth's.

**What this settles**
- The stage-13 build item 1 as formulated ("overwrite the dilation components of C, D3, D21 and k4 by the package at O(n^2)")
  fails its pre-registered test. It is the same family as our earlier tests of single amplitudes (note XXXVI 3k, all worse),
  now also for the joint correction.
- The mixture structure itself is real and sharper than we had it (section 9b: five slices, one charge, to 1-4% at deep layers;
  the 65-68% mixture share of V_l / 4). It explains why the counterterm table has the shape it has. It does not provide a lever
  for the score: the chain's visible error is not in the dilation sector.
- Note XLIV section 7's finding (the top-mode eigenvalue 0.5-0.8% low, inherited, transport factor 0.95) stays as measured, but
  its reading as "where the late error lives" is withdrawn. It is where the covariance error has the most energy; the output
  error, which is read through per-neuron variances and means, is not coherent along it.

### 9e. Pre-registration: is the collective block of the covariance error the lever? (committed before the runs)

Research switch `V63_TOPK` ("K:FILE:LAYERS:MODE", default off): at layers 1-14, after the pre-activation covariance is formed,
replace the chain's covariance by the truth's in the subspace of the true top-K eigenvectors U of C_T:
- MODE `block`: C <- C + (D P + P D - P D P), D = C_T - C_chain, P = U U^T (everything that touches the subspace; the bulk-bulk
  error stays);
- MODE `eig`: C <- C + U diag(lambda_j^T - u_j^T C u_j) U^T (the true top-K eigenvalues along the true eigenvectors only).

Monte Carlo covariance (1.6e7 samples) from `mc2_off{N}_full.npz`; its noise projected on K dimensions is about sqrt(K/n) of the
full noise, so the oracle is nearly noise-free. Network 0, counterterms on, baseline 1.4771e-8.

Predictions (raw MSE change): block K = 1: -3% (-8% to +3%; the top-mode oracle G4 was neutral); K = 8: -12% (-25% to -3%);
K = 32: -25% (-50% to -8%); K = 128: -40% (-65% to -15%); eig K = 32: -10% (-25% to +3%).

Decision rule: block K = 32 at or below -20% on network 0 (then confirm on network 1): the collective sector is a lever, and the
next step is a cheaper way to compute that block (a K-dimensional exact latent for the top modes). Between -20% and -5%: the
bulk matters as much; record the ladder. Above -5%: the collective covariance error is not the lever either, and the output
error is in per-neuron structure that no low-dimensional treatment reaches.

### 9f. Results of the top-K covariance oracle (`outputs/topk_net{0,1}.txt`)

Raw MSE, counterterms on, baselines 1.4771e-8 (network 0) and 1.5733e-8 (network 1):

| oracle | network 0 | network 1 |
|---|---|---|
| block K = 1 | 1.518e-8 (+2.8%) | |
| block K = 8 | 1.413e-8 (-4.3%) | |
| block K = 32 | 1.235e-8 (-16.4%) | 1.382e-8 (-12.2%) |
| block K = 128 | 1.007e-8 (-31.8%) | 1.316e-8 (-16.3%) |
| eigenvalues only, K = 32 | 1.604e-8 (+8.6%) | |

**Against the registration**: all block values fall inside the registered intervals (K = 1: -3% [-8, +3]; K = 8: -12% [-25, -3];
K = 32: -25% [-50, -8]; K = 128: -40% [-65, -15]; eigenvalues only: -10% [-25, +3], observed +8.6%, outside), at the weak end.
**Decision band** (K = 32 between -20% and -5%): the bulk matters as much; the collective sector is not the lever. The value
grows gradually with K (1: +3%, 8: -4%, 32: -12 to -16%, 128: -16 to -32%) instead of concentrating in a few modes, and
correcting the eigenvalues alone makes things worse (+8.6%): the same partial-correction failure as the gain oracle. For
reference the full variance oracle at every layer is -56% / -48% (note XLIII 4a). The collective-latent route (a K-dimensional
exact treatment of the top modes, my Route C) is closed as a main lever: a perfect oracle for the top 32 modes is worth 12-16%
and a real implementation would be worth less and cost more.

### 9g. How much of the truth's third and fourth cumulant is the mixture package? (`code/gain_share.py`, `outputs/gain_share_net{0,1}.txt`)

cos^2 between each slice and its package form, truth | chain (percent), network 0, layer 3 -> 14 (network 1 within 1-3 points):

| slice | truth, layer 3 -> 14 | chain, layer 14 | share of the chain's slice-error energy that is package-shaped, layer 14 |
|---|---|---|---|
| k3 diagonal | 78 -> 93 | 93 | 50% |
| D21 | 61 -> 86 | 59 | 32% |
| K22 | 98 -> 94 | 100 | 32% |
| K31 | 44 -> 63 | 96 | 35% |
| k4 diagonal | 99 -> 95 | 96 | 48% |

- **At deep layers 86-94% of the energy of the true k3, D21, K22 and k4 slices is the dilation package** (K31: 44-63%; the chain's
  own K31 is 96-99% package-shaped by construction, so the non-package half of the truth's (3,1) slice, the "history", is not
  represented at all). The non-package remainder is 6-14% of the energy of those slices.
- The chain's slice errors are 32-50% package-shaped in energy, yet correcting that component is harmful (section 9d): the
  Frobenius metric is the wrong weighting. The output reads these slices through coherent combinations that the counterterms have
  balanced.

**What this suggests, and what is untested.** If a design supplied the package exactly at O(n^2) and carried only the remainder
in the sources, the carried content would be about a third of the current amplitude (sqrt of 6-14%). With the same relative
accuracy the absolute error from the sources' transport would fall about 3x in amplitude, or equivalently the sources could be
carried at lower rank or precision for the same absolute error. This is a control-variate argument and it is not tested: the
oracle of section 9d tested the wrong thing for it (truth amplitudes against counterterms fitted for the old state). The test that
matches the claim is next.

### 9h. Pre-registration: does supplying the package analytically make the chain robust to lower rank? (committed before the runs)

The package amplitudes are pinned to the chain's own full-rank effective charges (`code/pin_targets.py`, network 0), so at full
rank the pin is a no-op, and applied to a chain run at reduced ranks it re-supplies, at each layer's readout and without
persistence, the package content that the compression loses. Runs on network 0 (raw MSE, cold FLOPs from flopscope):
- P0: full rank (R_OLD 320, R_OLD2 192) with the pin: a no-op, must reproduce 1.4771e-8 within 0.5%;
- P1: reduced ranks (R_OLD 160, R_OLD2 96), no pin;
- P2: reduced ranks with the pin.

Predictions: P0 within +-0.5%. P1: raw +15% to +60%, FLOPs -6% to -15% (the 128 cliff of note XXIV makes 96 risky; if P1 is above
+150%, repeat at (192, 128)). P2: removes at least a quarter of P1's degradation (P2 - baseline <= 0.75 (P1 - baseline)); 35%
confidence.

Decision rule: P2 removes at least 40% of P1's degradation and P2's adjusted MSE (raw x C/B) beats the baseline's: the package acts
as a control variate for the sources; next build a truth-free version (charges from the chain's own k3 charge and a universal
ratio profile) and run a cold 16-network screen. Removes 15-40%: record; the lever is real but small. Below 15%: close.

### 9i. Correction: the oracle of 9d and the chain deficits of 9b were measured with an implementation error (found by the pre-registered no-op check)

The pinned run P0 (section 9h) was registered as a no-op and was not: raw 1.7542e-8 against 1.4771e-8 (+18.8%) with charge
corrections of about 1e-6. Bisecting (`outputs/sens3_net0.txt`): the hook path with no modification reproduces the baseline
exactly; pinning only D3 gives +18.8% for correction scales -1, 0.2, 1 and 2 alike (a step, not a response), and the effect
grows with depth (layers 3-7: +0.07%, 8-11: +6.5%, 12-14: +10.4%). Cause: the hook multiplied D3 by the saturation mask
(neurons with alpha <= -2.5 dropped, 2% of rows at layer 3 and about 21% at layer 14). The production chain does not mask D3 (only
`V35_SAT_FULL` does, and it is off); it masks the rows of D21. Zeroing D3 of the saturated neurons is the "full drop" variant that
note XXIX found harmful.

A second, deeper error in the same measurement: **the charges of 9b and of the oracle targets were fitted over all neurons**. The
saturated neurons have the largest |mu|, so they carry a disproportionate share of every package form, and the chain's D21 has
zero rows there by design. Restricted to the neurons the output reads (alpha > -2.5, `outputs/dil_ledger_act_net{0,1}.txt`):

| layer 14, charge x 1e-3, truth \| chain | Rayleigh | k3 | D21 | K22 | K31 | k4 |
|---|---|---|---|---|---|---|
| net 0, all neurons (9b) | 8.56 \| 8.49 | 5.53 \| 5.27 | 5.48 \| 3.81 | 5.24 \| 4.30 | 5.53 \| 2.35 | 5.27 \| 4.47 |
| net 0, active set | 8.61 \| 8.54 | 5.56 \| 5.30 | 5.67 \| 5.45 | 5.09 \| 4.34 | 4.91 \| 2.46 | 5.14 \| 4.49 |
| net 1, active set | 7.50 \| 7.43 | 5.08 \| 4.81 | 5.19 \| 4.99 | 4.85 \| 4.21 | 4.33 \| 2.41 | 4.81 \| 4.23 |

**Retracted**: "the chain's D21 carries 70% of the charge" (it carries 96% on the active set; the 30% deficit was the saturation
drop, and it tracked the dropped fraction 2% -> 23% with depth), and "the five slices give one charge to 1-4%" (on the active set
the truth's third-order slices are 5.6-5.7 and its fourth-order slices 4.9-5.1 at layer 14, a ratio about 0.9, i.e. the stage-13
note's kurtosis-to-skewness ratio of 1.7-1.8 and not a pure mixture). **Stands**: the chain's Rayleigh charge is within 1%,
k3 within 5%, D21 within 4%, K22 and k4 low by 13-15%, and K31 low by about 50% (genuine: its regenerated shape is 96-99% package-shaped
and carries none of the truth's non-package half). The oracle results of 9d (G1 +322% etc.) are **withdrawn**: they were produced
with the wrong mask and overcorrected targets and say nothing about the question. The top-K covariance oracle of 9f does not touch these
slices and stands. 9g (cos^2 of each slice with its package form, all neurons) must be recomputed on the active set before it is
quoted. The registered ladder of 9c is rerun with the corrected hook (no D3 mask, charges and corrections on the active set only,
pin from a recorded baseline); the P0 check is now a no-op (1.4779e-8 against 1.4771e-8, +0.05%).
