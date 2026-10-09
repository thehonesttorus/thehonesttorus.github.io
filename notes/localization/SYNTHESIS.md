# XLII synthesis. Localization in a commutative home: error theory, four walls, what survives

Lead synthesis, 9 October 2026.

**Sources.**
- Frames F1-F6 and their referee reports. Referee corrections are applied throughout and take precedence.
- `DATA_FINDINGS.md` (100 official networks).
- README sections 1-4, including the pre-registered heat-defect experiment.
- The four Elicit surveys.
- The stage-8 introduction.

**Status labels.** **proved** (a proof in the cited frame, re-derived by its referee), **checked** (verified
numerically), **measured**, **conjecture**.

**Conventions.**
- 1 unit = 2n^3 FLOPs; B = 1024 units.
- adjusted = raw x max(0.1, C/B). So the score is raw x C on [0.1 B, B], and cost is free below 0.1 B.
- Baseline (V56, scored on 100 networks): raw 1.55e-8 at C/B 0.203 (208 units), adjusted 3.14e-9.
- Leaders: 1.1-2.1e-9 adjusted at 0.11-0.15 B, i.e. raw about 1.0-1.4e-8. Matching the best of them at our cost needs
  raw 5.4e-9 (a 2.9x cut). At 0.1 B it needs raw 1.1e-8 at half our cost.

## 1. The natural home

### 1.1 The heat flow on the ray quotient

**Setting.**
- u(m, Sigma) = E F(m + Sigma^{1/2} Z), and the truth is u(0, I).
- An estimator E(m, Sigma) is the chain run on the input law N(m, Sigma).
  - It is *point-mass exact* if E(x, 0) = F(x). Every cumulant chain is.
  - It is *homogeneous* if E(cm, c^2 Sigma) = c E. Every chain built from Gaussian relu moments and Edgeworth terms
    is.
- The heat operator is H = d_Sigma - (1/2) Hess_m. The truth satisfies Hu = 0 (Price). The *defect* is D = HE.

**Theorem 1 (the error is the integrated heat defect; Dynkin). Proved.**
- **Claim.** For any Gaussian linear-tilt localization (m_t, Sigma_t) of N(0, I), with dm = Sigma C dW and
  dSigma/dt = -Sigma C C^T Sigma (Eldan: C = I), and any C^{2,1}, point-mass-exact E:

      truth - E(0, I) = -E int_0^inf < D(m_t, Sigma_t), Sigma_t C_t C_t^T Sigma_t > dt.

  An estimator that is defect-free and point-mass exact equals u.
- **Where it appears.** F1 Thm 1. Found independently in F3 Thm 1, F4 Prop 5.3, F5 Thm 4.2 and the Elicit survey.
  The F1 referee checked it end to end at n = 6, L = 3.
- **What it is.** Dynkin plus Price. The merge built on it in section 3 is Talay-Tubaro / Zadunaisky defect
  correction.
- **ARC reading.** Iterated estimation becomes the objective statement "E(m_t, Sigma_t) is a martingale", and its
  integrated violation is the error. F1's claim that no efficient estimator can satisfy it is demoted to a remark.

**Theorem 2 (ray quotient and Green weight). Proved (F1 Cor 2); checked.**
- **Setting.** For a homogeneous E put e(y) = E(y, I) and delta(y) = e - y.grad e - Lap e.
- **Claims.**
  - tr D(m, sI) = delta(m/sqrt s)/(2 sqrt s).
  - The error is a pairing with a Green measure:

        truth - e(0) = -(1/2) int_0^inf (1+tau)^{-3/2} E_{y~N(0,tau I)} delta(y) dtau = -<G, delta>.

    G has mass 1 and tail G(tau > T) = (1+T)^{-1/2}.
  - In log time T = log(1+tau), y_T is the repulsive OU process dy = (y/2) dT + dB. The truth solves
    (A - 1/2) v = 0, with A = (1/2)(Lap + y.grad) and the network as boundary data at infinity.
- **The base point.** delta(0) = e - Lap_m e = 2 tr D(0, I). The scale direction is empty: tr d_Sigma E = e/2.
- **The detector form.** If the localization profile b(s) (the error left at posterior variance s) is ~ err s^p
  near s = 1, then **err = -delta(0)/(2p)**.
- **Measured** (Gaussian closure, GC, n = 64):
  - half of the error is made at tau <= 1, under 10% beyond tau = 9;
  - the profile factorizes per neuron to 2-8%;
  - p = 1.1-1.2 for GC and 2.3-2.8 for the dense kappa_3 chain (K3).
- **Reading.** "The measure on the space of leaves passes to the leaf space" is exact here:
  - the leaves are scaling orbits;
  - the estimand is leafwise linear;
  - the transverse measure is the occupation (Green) measure of the reduced flow.

  It is commutative.

**Theorem 3 (Euler-Stein ladder). Proved; checked by Monte Carlo.**
- **Claim.** With c_k = E D^k F: tr_{(k+1,k+2)} c_{k+2} = (1-k) c_k.
  - Rung 0: E F = E Lap F (mu = tr H). Checked: slope 0.983, correlation 1.0000.
  - Rung 1: the third chaos is traceless. Checked at its noise level.
  - Rung 2: E Lap^2 F = -E F. Checked: slope 1.08.
- **Reading.**
  - These are truth-free equalities, so they carry no Gaussian slack (section 1.3).
  - The chain's rung-0 violation at the base point is its defect. Lap e is a second estimator of the truth, the
    Euler-Stein companion.
  - The merge e + a delta = (1+a) e - a Lap e is the oblique projection onto rung 0.
  - Rung 0 already explained why the dressed (Schur) hub beat the bare one.
  - Higher rungs need about n^2 chain runs.

**Theorem 4 (all-orders trickle-down). Proved.**
- H kappa_k = (1/2) sum_{A u B = [k]} grad kappa_|A| (x) grad kappa_|B|, and the posterior cumulants obey a matching
  SDE.
- k = 2 is the Anari-Koehler-Vuong (AKV) trickle-down.
- It accounts for variance; bias is Theorem 1's job.

**Theorem 5 (exact layer telescoping). Proved (F1 Prop 8); checked to 4 digits.**
- **Claim.** delta = 2 sum_l S_l . tr d_l.
  - d_l is layer l's local defect relative to the hierarchy.
  - S_l is the chain's linearization from layer l to the output, applied by a forward pass through all channels
    (referee R3.1).
- **Measured split** (GC, n = 48, L = 8):
  - the mean and variance readouts cancel;
  - the defect sits in the covariance map's cross term (the first-chaos-projected (2,1) slice, n^3 per source-layer)
    and its diag-type term;
  - the n^4 pair Gram carries <= 3%.

**Proposition 6 (first-order visibility).**
- **Claim.** -tr D of a Gaussian readout is the Edgeworth correction of the infinitesimal-leaf mixture
  (kappa3_loc = 3<grad mu, grad v>, kappa4_loc = 3|grad v|^2).
- **Class weights** in the second-chaos model, each relative to the class's own bias:

  | class | weight |
  |---|---|
  | kappa_4 path class | 1 |
  | kappa_3 path class | 2 |
  | the gain's kappa_3 half | 1 |
  | the gain's kappa_4 half | 0 |
  | the readout's closed walks | 0 |

- **Status.** Proved for a single readout (F5 Prop 5.4, with its referee's correction). Open for composed chains: F1
  argues that closed-walk *omissions* in carried tables are seen through d_Sigma.

**Corollary 7 (visibility bound; new here). Proved, given the class model.**
- **Claim.** For uncorrelated class components with energies E_c and defect weights w_c:

      X* = (sum w_c E_c)^2 / (sum E_c . sum w_c^2 E_c) <= (visible share).

  An equal mix of weight-1 and weight-2 classes loses only 10%.
- **Applied to the production chain.**
  - The oracle attribution (note XXXI) puts about 40% of the MSE in the kappa_3 readouts (their open and closed
    parts), 15-30% in the kappa_4 diagonal, and about 30% elsewhere: closed walks, pair gates, truncation,
    calibration.
  - That gives a visible share of about 0.35-0.7, so X* is about 0.3-0.65. This matches the registered 0.45
    [0.25, 0.70].
- **A falsifiable consequence.** X* above 0.7 would show that closed-walk omissions are visible in a composed chain.

### 1.2 The groupoids, and why they are type I or abelian

| groupoid | type | transverse measure | content |
|---|---|---|---|
| scaling R_+ on Gaussian laws and on inputs (rays) | proper free action with a smooth quotient (slice Sigma = I; S^{n-1}): type I | Green measure G; sigma on S^{n-1} | radial Rao-Blackwellization is exact and removes Var(r) E\|F(u)\|^2 = 1/(2n) = 0.69% of sigma^2 |
| rotations O(n) on S^{n-1} | Gelfand pair (O(n), O(n-1)): the commutant is abelian and equals the centre | nu_F = sum a_m delta_m, with kink tail a_k ~ 0.127 (L +- K'(0)) k^{-5/2} | every isotropic (W-oblivious) scheme acts on the MSE by a nonnegative multiplier |
| gate-pattern relation R_l | (+)_s B(L^2(R_s)) (x) 1: type I. Its centre is a commutative AF tower of feasible prefixes | the law of the pattern | E F_i is a facet sum (Stein); tropically, E F_i = (V_1(P) - V_1(Q))/sqrt(2 pi). Exact, but no algorithm |
| linear localization R_U | Morita equivalent to C_0(R^k): type I | gamma_k | ordinary Gauss-Hermite quadrature |

**Why everything is type I or abelian.** Every quantity the chain carries is an expectation, under one classical law,
of commuting variables: pre-activations, gates and Gaussian observations. The GNS representation is therefore a
multiplication algebra (F4 Prop 1.2). The only non-type-I candidate, the tail relation of gate patterns as
L -> inf, is irrelevant at L = 16.

**What each stage-8 theorem becomes.**

| stage-8 theorem | here |
|---|---|
| A: localization in the commutant | the Cartan part is the classical input localizations, and the centre is the scale / pattern / rotation disintegration. "Refines the central one" holds only at terminal time |
| B: modular tilt, KMS response | all Cov^f coincide, the Petz inequalities are equalities, and the Wigner-Yanase skew information is 0. What remains is AKV: gate pinning is an exact linear tilt with Sigma - E Sigma_after = (1/N) Sigma D^{-1} Sigma, and input localization is the Wiener-chaos decomposition |
| C: localization capacity | the law of total variance (with Stein and Poincaré for k directions; a factor-analytic diagonal-plus-rank-K capture for gates) |
| D: Kikuchi-Toeplitz, Berezin-Lieb | every two-sided partition-function sandwich pins the means together (F5 Thm 2.1). The spherical lower symbol is Jensen to 1e-3 (beta_2(l) ~ 2l/n). The gate Kikuchi hierarchy is the Mehler expansion, and its level-2 projection G is read linearly |
| E: pinning, Oppenheim | the spectral-independence constant is 3.4-10.9 and does not grow with n. It controls Glauber mixing, not closure bias (the Elicit pinning survey concurs) |

### 1.3 Where noncommutativity exists, and the Gaussian-slack no-go

1. **The arrow (Schur) algebra of a layer.**
   - **The objects.** M_n as functions on pairs of neurons. The state is the unit data d (variances, or a kappa_4
     diagonal); the modular operator is Delta e_ij = (d_i/d_j) e_ij. The ReLU positive-diagonal gauge acts by
     coboundaries.
   - **Theorem (F4 Thm 2.1; proved).** A closure of an arrow statistic from unit statistics is gauge covariant iff it
     is a weighted geometric (KMS) mean.
   - **What it explains.** The adopted geometric slice, V33_WK4M=3.
   - **Why it is not a lever.** The arithmetic-geometric gap is about cv^2/4 of a slice whose route share is 0.016,
     and covariant counterterm shapes measured worse (PMETRIC +5.7%).
2. **Compressed multiplication operators.**
   - **The structure.** The chain solves a truncated moment problem. Its operators P_k z_a P_k on degree-<=k
     polynomials do not commute; they are the Berezin-Toeplitz structure of stage-8 Theorem D. A closure is a choice
     of commuting (flat) extension (F4 referee).
   - **In chaos form.** The second-chaos tuple (H_i) in (M_n, tr) has vector-state moments L^T H H L, which are open
     walks reached at n^3 by the Schur hub, and trace-state moments tr(H_i H_j H_k), which are closed walks at n^4.
   - **Conclusion.** The noncommutative distribution of the second chaos is the cost wall.
3. **The Gaussian-slack no-go.**
   - **The geometry.** The layer law sits about 2 var^2 inside the moment cone (exact Gaussian slack
     V - Dt^T C^{-1} Dt = 2 C o C). A closure must supply non-Gaussian quantities of size r var^2, with
     r = kappa_4/(2 var^2).
   - **The numbers.** Median r is 0.61 / 0.39 / 0.16 / 0.085 at n = 32-256 (it halves per doubling), and about 0.03
     at n = 1024.
   - **The consequence.** Every certified constraint (moment positivity, flat extension, Petz / Kubo-Ando convexity,
     SDP relaxations) is inactive unless the closure error exceeds **about 1/r ~ 30x** the closed quantity (3-30x;
     median about 15). The chain's kappa_4 error is 0.3x.
   - **Status.** Proved for the diagonal Schur-complement constraint; the truth violates it on 0 neurons. Conjectured
     for the full degree-2 matrix. False in near-null directions, which the output does not read.

## 2. The four walls

**W1. Pinning.**
- **The class.** The first class that pinning or Kikuchi adds beyond the chain's frozen product-gate transport is the
  pair-gate (joint-gate triangle) class.
- **Its value.** At most about 1/3 of the MSE (about 6.5e-9; one network, one layer). Its annealed mean is zero, so
  no counterterm absorbs it. It therefore pays only if carried for **Delta C < C/2 ~ 104 units**.
- **Its cost.** The cheapest CP carrier is **about 40 units per separating direction**, and the factor-analytic K_90
  is predicted at **about 0.04n = 40-45** at layer 10 (referee extrapolation from n <= 512). So only K <= 2 is
  affordable, capturing 3-6% of the MSE.
- **Corrected ratios.** Adjusted x1.53 at K = 4 and x3.07 at K = 16. Exact carrying costs n^4 per source-layer.
- **Neighbours.** Kink-graded mean field fails at O(1). The input radius, the only gate-neutral localization, is
  worth 0.69% and is already used.

**W2. Collective localization: the pricing law.**
- **Statement.** A k-dimensional localization needs N >= k+1 leaves, and a leaf costs a chain (proved for frozen
  leaves). Splitting a total C* >= 0.1 B into N leaves pays only if

      rho < r(C*)/r(C*/N) <= 1/N

  on a raw x C-optimal family (F3 Thm 9, as corrected). Under a local power law r ~ C^{-e} the bar is
  rho < N^{-e}.
- **Measured.** rho N = 1.05-1.55 (k = 1), 1.2-2.7 (k = 2), 1.5-6.1 (k = 4). It is above 1 on every network at
  n <= 256.
- **The production chain.** It carries the gain mode, so rho_1 is predicted in [0.85, 1.3]. That puts LL-1
  (0.406 B) at 5.3-8.2e-9 adjusted.
- **Scope.** The same law prices mixtures and pattern windows.

**W3. Sampling hybrids.**
- **Statement.** For any positive-weight N-node rule with nodes independent of W, **N x MSE >= 0.56 sigma^2 at
  0.1 B and 0.54 sigma^2 at B**. The bound is a certified Delsarte-Yudin LP, and frames plus antipodes come within
  7-11% of it.
- **Consequence.** An unbiased add-on improves the chain by at most b^2/(b^2 + v): **0.24% at 0.1 B** and 2.5% at B,
  while C/B rises to 0.303 or 1.203. Sampling alone floors at 6.1-6.4e-7.
- **The open edges stay closed.**
  - The signed-weight floor (0.27 sigma^2) is a conjecture.
  - W-dependent control variates explain only about 26% of sigma^2 against the 99.86% needed.
  - The error is co-localized with the variance: the oracle directional hybrid gains 0.2-0.7%.

**W4. Certified and convex-order bounds.**
- **The mean is invisible.** Any two-sided MGF sandwich pins the means, so Ising <= boson <= Ising and every
  Berezin-Lieb sandwich bound fluctuations only.
- **Convex order does not propagate.** H_i is indefinite (52% of its eigenvalues are negative), and matched closures
  are convex-incomparable with the truth.
- **Moment sandwiches** (Markov-Krein, 4-8 moments) are **50-1500x** the target width. Closing them needs about 78%
  of the variance certified Gaussian.
- **Gaussian slack** leaves every certified per-neuron constraint inactive by about 30x.
- **Measured in passing.** The gain sign rule is exhausted (truth < chain on 0.47-0.55 of neurons; |corr| < 0.03
  with the gain shapes), so no gain-template counterterm remains.

## 3. Leverage that survives, ranked by expected adjusted gain

| rank | item | added cost | adjusted (from 3.14e-9) | P(pays) | gate |
|---|---|---|---|---|---|
| 1 | **A.** heat-defect merge, late in-place form | +20-50 units (estimated) | 2.1-2.8e-9 at X* kappa_a = 0.45-0.6 | ~0.2 | the X* run plus the defect-injection re-analysis |
| 2 | **B.** V57_CHOL | -7.7 to -8.2 units | 3.02e-9 (-3.6 to -3.9%), raw unchanged | ~0.85 | E-B at b = 128, with the residual measured on grader terms |
| 3 | **C.** frontier scan | none billed (24 offline chains) | gains iff r(0.1 B)/r(0.2 B) < 2.03 | ~0.2 | the scan itself |
| 4 | **D.** small items | ~0 | <= 1-3% | low | listed in D |

Expected values: A about -6%, B about -3%, C about -2%. Run B and C first, since they are cheap; A's measurement goes
in the same fleet session.

### A. The heat-defect merge

**The component.**
- e_DR = e + a delta_hat, with a = -1/(2p): the OLS blend (1+a) e - a Lap e of the chain with its Euler-Stein
  companion.
- a must be refit on the calibrated chain. In the referee's MF test, one counterterm moved the slope from about
  -0.85 to about -0.69, while the explained share held at 0.74-0.92.

**Small-n evidence** (two independent replications).

| chain | detector explains | slope a | held-out merge removes |
|---|---|---|---|
| GC (n = 48-256) | 59-98% | -0.43 pooled | 76% |
| K3, ten times more accurate | 74-88% | -0.18 to -0.22 | 81% |

The parameter-free midpoint (e + Lap e)/2 cuts GC 8.7x, part of that a common shift that production lacks
(de-meaned 0.04-0.46).

**The pending n = 1024 measurement** (README section 4).
- **Setup.** Networks 0-3, counterterms on, float64, saturation drop off. K = 128 directions with 4 evaluations each:
  about 2,050 runs, 1-2 instance-hours.
- **Registered predictions.** X* = 0.45 [0.25, 0.70], slope -0.15 to -0.25, rho 3-6.
- **Registered rule.**
  - X* >= 0.5 with a held-out slope stable within 20%: derive the late-layer local form, and build it only if it
    captures >= 70% of delta for <= +70% of the bill.
  - X* < 0.3: the identity stands as explanation only.

**Additions from this synthesis** (they sharpen the registered rule; they do not change it).
1. **Break-even.** The adjusted ratio is (1 - s X* kappa_a)(C'/C):
   - s is the MSE share the implemented defect covers;
   - kappa_a is its analytic capture;
   - C'/C is the cost factor.

   For the full-depth form (s = 1, the referee's truncated transport at C'/C = 1.5-1.9), the merge needs
   X* kappa_a >= 0.33-0.47. At the registered +70% cap it needs 0.41, so the cap is loose: decide on a predicted
   ratio <= 0.9.
2. **A-late** (cross-reading F1 Thm 5 with DATA_FINDINGS).
   - **Where the error is made.** It is a sum of uncorrelated per-layer injections under contractive mean-gate
     transport (about 9x loss between layers 3 and 15). Layers 12-15 inject 11 / 15 / 18 / 21% of it, 65% together;
     layers 9-15 inject 88%.
   - **The design.**
     - Truncate the defect to the local defects of the last 3-4 layers: about 1/4 of F1's +70-140 units, plus the
       first-chaos legs at about 1 unit per layer.
     - Inject it in place: mean_l <- mean_l + a_l dinj_l. The chain's own transport then carries the correction as
       it carries the error, which removes the referee's separate forward pass (+35-400 units).
   - **Estimated cost.** +20-50 units (C'/C = 1.10-1.25); up to about +85 if the late source-layers are young-dense.
   - **Predictions:**

     | X* kappa_a | A-late raw | A-late adjusted | A-full raw | A-full adjusted |
     |---|---|---|---|---|
     | 0.30 | 1.25e-8 | 2.8-3.2e-9 | 1.09e-8 | 3.3-4.2e-9 |
     | 0.45 | 1.10e-8 | 2.4-2.8e-9 | 8.5e-9 | 2.6-3.3e-9 |
     | 0.60 | 9.5e-9 | 2.1-2.4e-9 | 6.2e-9 | 1.9-2.4e-9 |
     | 0.80 | 7.4e-9 | 1.7-1.9e-9 | 3.1e-9 | 0.9-1.2e-9 |

   - **Choice.** The crossover is at X* kappa_a of about 0.55-0.65. Corollary 7 favours the low side, so A-late is
     the likely build. Only A-full at X* kappa_a >= 0.8 goes below the leaders' front.
3. **The defect-injection test** (no new compute; it gates A-late).
   - **Inputs.** The pending run's per-layer delta_hat_l and the chain dumps.
   - **Per layer.** inj_l = e_l - diag(Phi(alpha_l)) W_l e_{l-1}, and dinj_l is the same with delta_hat in place of
     e. Compute the noise-corrected X*_inj(l) and slope a_l at layers 9-15.
   - **Build A-late iff** X*_inj >= 0.4 at layers 12-15, with a_l stable within 20% (fitted on networks 0-1, applied
     to 2-3).
   - **kappa_a** comes from the exact-defect run on network 0 (2049 chains, as F1 specifies).
   - **Caveat.** dinj_l also carries the one-step image of upstream variance and covariance defects. If those
     dominate, A-late's cost moves toward A-full.
4. **Freeze, do not drop.**
   - **The problem.** The n = 256 pilot found the defect h-dependent from layer 13 on, exactly where 54% of the
     error is injected.
   - **The remedy.** If the n = 1024 pilot fails there, freeze every mask, count, basis, clip and max at the base
     law (F1 R4.1). The frozen chain equals production at the base point and is analytic. Add 8 unfrozen control
     directions, and subtract the truth-free stencil bias (h^2/8) e/(n+2).
   - **Otherwise** only layers 1-12 (35% of the error) are gated.

**Risks.**
1. The visibility ceiling of about 0.6-0.7.
2. kappa_a is unknown. The production chain's local defects sit one level up (kappa_4 and kappa_5 with an input leg),
   and their factoring through first-chaos legs is unshown.
3. a must be refit after calibration.
4. A positive gate is followed by a long build.

### B. V57_CHOL

- **What.** Compute the hub diag(Y M^{-1} Y^T) (M symmetric positive definite) as the column norms of L^{-1} Y^T,
  using Cholesky plus a blocked forward substitution in place of the LU solve.
- **Cost.** 0.706 units per layer (b = 64) or 0.745 (b = 128), against 1.334. Over 13 layers that is -8.2 or -7.7
  units, taking C/B to 0.195.
- **Effect.** Raw is unchanged (float32 relative error 5-8e-7). Adjusted becomes 3.02e-9 (-3.6 to -3.9%).
- **Risk (binding): grader-counted residual time.** At b = 64 the loop issues about 133 small ops per layer
  (180-200 µs each), about +0.1-0.25 s against a 0.4 s gate; b = 128 halves it. A Cholesky breakdown is the minor
  risk: a NaN guard catches it, with eps = 1e-2 or LU as the fallback.
- **Gate (F4 E-B).** Cold on networks 0-15, paired; then scored on 100. Adopt iff adjusted improves by >= 3% with
  raw within +-0.5% and the residual under 0.4 s on grader-equivalent hardware.
- This is not a theory win.

### C. Frontier scan

- **What.** The production chain on networks 0-7 at C/B of about 0.05, 0.1 and 0.2 (24 chains). Use the existing
  switches: V21_R_OLD and V24_R_OLD2, the age boundary, and the young ranks.
- **Why.** W2 and every "spend less" idea rest on r(C) below the operating point, which has never been measured.
  - The measured elasticities are only local: -0.91 and -0.93, score-neutral at the margin.
  - The leaders' raw of 1.0-1.4e-8 at 0.11-0.15 B shows that frontiers below 0.2 B are not pinned by our slope.
- **Decision.**
  - 0.1 B pays iff r(0.1 B)/r(0.203 B) < 2.03; 0.05 B pays only under the same ratio, since cost is free below
    0.1 B.
  - With a local exponent e, the adjusted change at 0.1 B is x 2.03^{e-1}: -6% at e = 0.91, -30% at e = 0.5.
  - Adopt iff paired adjusted improves by >= 3% on >= 6 of 8 networks; record e either way.
- **Prior.** e >= 1 below 0.1 B (the young-tier floor), so likely neutral. Finding out costs minutes.

### D. Small items

- **D1. Visibility prediction** (Corollary 7). Reads the X* outcome as a statement about which omitted classes are
  visible. Free.
- **D2. The inside-U triangle** (F3 Thm 3 with F2 Prop 6).
  - **Object and cost.** tr(U^T H_i U)^3 at k^2 n^2 per source-layer, about 0 units.
  - **Value.** At n = 32 the inside-U share is 0.07 (k = 1) and 0.27 (k = 4); unmeasured at n = 1024. The value is
    at most about 1-3% of raw, and its annealed mean is zero.
  - **Risk.** It overlaps the carried gain mode.
  - **Gate.** Build only if the share on official sources is >= 0.15.
- **D3. V58_BARYGATE** (F4).
  - **What.** The first vertex Phi(mu/sigma_eff), O(n) per layer.
  - **Prediction.** Raw -3% to +2%.
  - **Gate.** A cold 16-network screen. It closes the k = 0 branch.

## 4. Bottom line

**What the programme gave.**
1. **An exact error theory, entirely commutative.** The error is the Green-integrated heat defect of the
   ray-quotient flow. It is accompanied by truth-free ladder equalities, an exact layer telescoping, and per-class
   visibility weights.
   - The NCG vocabulary maps one-to-one onto classical objects:

     | NCG term | classical object |
     |---|---|
     | transverse measure | Green / occupation measure |
     | Cartan subalgebra | classical localizations |
     | centre | scale / pattern / rotation quotient |
     | KMS mean | geometric Schur mean |
     | Toeplitz operator | compressed multiplication operator |

   - The stage-8 lens adds one real thing: a precise address for the cost wall (the noncommutative distribution of
     the second chaos), whose constraints are about 30x slack.
2. **Four no-go walls with numbers.**
   - pinning: < 104 units affordable against >= 40 per direction;
   - collective localization: rho < 1/N;
   - sampling: 0.54-0.56 sigma^2, so 0.24%;
   - certified bounds: mean invisibility, 50-1500x, about 30x.
3. **One positive phenomenon.** At small n, the base-point heat defect is a truth-free, per-neuron error detector
   with an essentially universal coefficient. Its survival at n = 1024 is the one open measurement.
4. **One engineering item**, V57_CHOL.

**What a genuinely new Phase 2 system must look like.**
1. **One deterministic, W-reading chain.** No samples (W3). No leaves or mixtures (W2): localization lives inside
   the chain as tangents, not leaves. No certificates (W4).
2. **Omitted classes are carried as open walks at n^3, or left alone.** The pair-gate / closed-walk class has no
   carrier under C/2 (W1). No localization makes the gates act exactly on the carried object.
3. **Gains must target the incoherent per-neuron error in layers >= 9.** That is where 88% of it is injected, and
   the coherent part is gone (signed mean error < 1e-5 per layer).
4. **The one localization-native design is the self-correcting chain.**
   - It carries the first-chaos tangents of its late tables, evaluates its own heat defect, and cancels it in place:
     A-late, i.e. classical defect correction applied to a moment closure.
   - With CHOL it lands at **2.0-2.7e-9** if X* kappa_a = 0.45-0.6. That is inside the leaders' band, not below its
     front.
5. **The leaders' route looks like cost, not accuracy.** On these walls it can come only from cheaper constants in
   the open-walk young tier (69% of our bill): CHOL-type linear algebra, and whatever e the frontier scan finds.
   Localization does not supply it.

**Next steps, in order.** E-B; then the X* run with the freeze fallback and the defect-injection re-analysis; then
the 24-chain scan. If X* < 0.3, localization closes as an explanation.
