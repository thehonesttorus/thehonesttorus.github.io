# F4 - Modular tilts, KMS/Petz covariances and signed closure inequalities for the gated covariance recursion

Frame F4 of the localization round. Question: is the gating step (relu of a near-Gaussian vector, then W) a modular
tilt or a conditional expectation in a matrix algebra, so that the stage-8 Theorems B (noncommutative trickle-down,
one-sided Petz inequalities, Wigner-Yanase defect) and C (Bures capacity) act on the chain? Targets: (a) signed
bounds on closure errors from operator concavity; (b) whether a geometric or harmonic mean of covariances is the
right object somewhere; (c) whether the trickle-down identity for a bank of response observables gives an exact
accounting of what the truncation drops.

Code: `f4code/` (all runs local, small He networks, row convention z_l = x_(l-1) W_(l-1)); data and outputs: `f4data/`.

## 0. Verdict

**The frame is negative for accuracy, with three structural results and one cost item.**

1. **The gating step is neither a tilt nor a fixed conditional expectation, and in every representation built from
   the network's law the algebra is commutative** (section 1). There Theorem B is the classical
   Anari-Koehler-Vuong identity (its Petz inequalities become equalities, the Wigner-Yanase defect is 0) and Theorem C
   is the law of total variance. The post-activation law is not absolutely continuous with respect to the
   pre-activation law (atoms on the orthant faces), so the gate is a hard pinning, outside the same-algebra modular
   transport. Both theorems
   were checked numerically and hold (section 6.1); they just have no content here.
2. **Noncommutativity lives in one place: the arrow (Schur) algebra of a layer** (section 2). This is M_n as
   functions on pairs of neurons, with unit data (variances, kappa4 diagonals) as the state. The ReLU positive-diagonal
   gauge acts on it by modular-cocycle coboundaries.
   - **Theorem 2.1:** among all ways of closing an arrow statistic from unit statistics, only the weighted geometric
     (KMS) mean is gauge covariant.
   - It explains the adopted geometric (2,2) slice and the K4D = 3 diagonal, and on small nets the geometric mean wins
     28/28 layer-network cases (section 6.3).
   - It is not a lever. At width 1024 the arithmetic-geometric gap is about cv^2/4 of a slice that the output barely
     reads. The renormalized chain also rejected a gauge-covariant K31 shape (PMETRIC, +5.7%), so covariance is
     necessary for exactness but does not by itself improve a fitted counterterm.
3. **(b), the right mean of the covariance, has an exact answer for radial mixtures** (section 3). A k-homogeneous
   readout sees the k-th power mean of the scale.
   - The mean gate (k = 1) sees the Bures-Wasserstein barycenter scale (E sqrt G)^2.
   - The covariance (k = 2) sees the arithmetic E G.
   - The W2 barycenter gives a certified one-sided bound on E relu for every Gaussian mixture (Theorem 3.2), with
     equality on scale families.
   - The chain already carries this through its kappa3/kappa4 Edgeworth terms (the gain to 2%, note XIII E10). On small
     nets the barycentric correction explains 7-56% of the per-neuron mean-gate error, while the 4th-order Edgeworth
     term with exact cumulants explains 96-99% at n = 64. The gate formula is not where the error is; the cumulants
     are.
   - Gates and the arcsine kernel (k = 0) are invariant under radial mixing, so the bare vertex Phi(mu/sigma) of the
     legs and births is radially biased.
     - The barycentric vertex explains 30-68% of the per-neuron vertex defect on small nets (section 6.6).
     - It gives one free candidate, `V58_BARYGATE`: a bounded resummation of the unstable V53 tadpole. It is
       predicted near-neutral and listed to close the branch.
4. **(a), signed closure inequalities: every certified one carries the Gaussian slack, which is about 1/r times the
   quantity being closed** (section 4). Here r = kappa4/(2 var^2), and r falls as 1/n.
   - Measured r is 0.61, 0.39, 0.16, 0.085 at n = 32, 64, 128, 256 (layer 8). That extrapolates to about 0.03 at
     n = 1024, in line with the official-network kurtosis 0.066 (r = 0.033, note XIII).
   - The certified per-neuron constraint therefore activates only when a kappa4 closure error exceeds about 3-30
     times kappa4 itself (median about 15). The chain's error is 0.3 times kappa4.
   - Sharp inequalities exist only inside a model. The second-chaos law gives kappa4 kappa2 >= (4/3) kappa3^2
     (Theorem 4.2, sharp, proved). True laws break it on 0-12% of neurons and satisfy it with a factor of 2-11 to
     spare. Inside the model it is a statement about the model, not a certificate.
5. **(c), trickle-down for a response bank: exact accounting, but of variance, not of the chain's bias** (section 5).
   - Under Gaussian localization of the input it is the Wiener-chaos decomposition (Theorem 5.1).
   - Under localization of the layer field with noise covariance C, its t = 0 rate for the bank {C_i} is exactly the
     adopted Schur hub of note XLI (Theorem 5.2). What the hub omits is revealed only at t > 0, through nonlinear
     regression on the field, which is note XL's n^4 wall again.
   - The mean's bias obeys a different exact identity (Proposition 5.3): it is the integrated heat-equation
     (Euler-Stein) defect of the chain along the Gaussian bridge. It is exact and not computable at this budget.
6. **Estimator.** No frame-derived accuracy change survives.
   - One cost change comes from the frame's re-reading of the hub as a Schur complement: compute it with a Cholesky
     factor and a blocked forward substitution instead of an LU solve.
     - flopscope bills it at 0.706 units per layer instead of 1.334, measured at n = 1024 (section 6.5). That is
       -8.2 units in all, so C/B 0.2028 -> 0.1948 and adjusted about -3.9% at unchanged raw (section 7).
     - The float32 error is about 7e-7 on a realistic M with cond 5.8e4.
   - It is an engineering consequence, not a theory win, and it carries a float32-conditioning risk, which is stated.

## 1. What has to hold for Theorems B and C to bite, and why it does not here

**Stage-8 statements used** (from the supplied intro). State rho on M_d, observables A_j with means m_j. For a
symmetric operator-monotone f, Cov^f_rho(A_i, A_j) = Re <A_i, m_f(L_rho, R_rho) A_j>_HS - m_i m_j, with the
arithmetic mean giving SJ (Jordan), the geometric K (KMS) and the harmonic H (Bures).
- **Theorem B.** For rho' = rho^(1/2)(1 + <A - m, Z>)rho^(1/2), E Z = 0, E ZZ^T = C:
  - m' - m = K Z;
  - E SJ(rho') = SJ(rho) - KCK (exact);
  - E Cov^f(rho') <= Cov^f(rho) - KCK for every f.
- **Theorem C.** Every localization has Cov(m_inf) <= H.

**Proposition 1.1 (commutative collapse).** If rho commutes with every A_j:
- Cov^f = SJ for every f (so K = H = SJ, and the Wigner-Yanase and Fisher defects SJ - K and 4(SJ - H) vanish);
- the modular tilt is the linear tilt rho(1 + <A - m, Z>);
- the inequalities of Theorem B are equalities (the AKV identity);
- Theorem C reads Cov(m_inf) <= Cov, which is the law of total variance.

*Proof.* In a joint eigenbasis A_kl = 0 for k != l, and m_f(lambda, lambda) = lambda, so (eq. covf of the intro) every
Cov^f equals sum_k lambda_k A_kk B_kk - m_A m_B. For the tilt, rho^(1/2) X rho^(1/2) = rho X when X is diagonal. With
all Cov^f equal and E Cov(rho') = Cov(rho) - KCK exact for the Jordan form, the inequalities hold with equality. ∎

**Proposition 1.2 (the network's law generates a commutative algebra).** Every quantity the chain carries is an
expectation, under the classical law of the input, of products of the commuting random variables z_l, y_l and gates
1[z > 0]:
- moments and cumulants of every layer;
- conditional expectations given gates;
- conditional expectations given Gaussian observations of the input or of a layer.

The GNS representation of that law on the algebra they generate is a multiplication algebra, and Proposition 1.1
applies. A non-commutative structure therefore has to come from an algebra the network's variables do not generate.
∎

**Is the gating step a tilt, or a conditional expectation?**
- **Not a tilt.** A tilt is absolutely continuous: d nu'/d nu = h. The law of y = relu(z) gives positive mass
  Pr(z_i <= 0) to the face {y_i = 0}, while the pre-activation law has a density and gives every face mass 0. So
  law(y) is not absolutely continuous with respect to law(z), and no density h exists. In localization language the gate is a hard pinning: the support changes. The
  stage-8 and brief's same-algebra KMS transport explicitly excludes this case.
- **A state-dependent corner, not a channel.** For a fixed realization the map X -> Q X Q, Q = diag(1[z > 0]), is a
  conditional expectation of M_n onto the corner Q M_n Q (idempotent, positive, Q M_n Q-bimodular). The
  post-activation second moment is E[Q(z) z z^T Q(z)], in which the corner is chosen by the vector it compresses.
  - For an A independent of the gate, E[QAQ] = P o A (P = E qq^T) is a Schur-multiplier channel.
  - Operator Jensen holds for it: E[QAQ] - E[Q] A E[Q] = E[(Q - EQ) A (Q - EQ)] >= 0, the "conditioning is not
    freezing" identity of the brief.
  - In the network the dependence makes the update the nonlinear Mehler (Schur-power) map
    Cov(y) = sum_(k>=1) (c_k c_k^T / k!) o R^(o k) o (sigma sigma^T), with c_k the Hermite coefficients of relu at a_i.
- **Proposition 1.3.** Cov(y) >= d(Phi) C d(Phi) (Loewner order) for Gaussian z. Every Mehler term is PSD by Schur's
  product theorem, and the k = 1 term is d(Phi) C d(Phi). So the product-gate (first-chaos) covariance is a certified
  lower bound. It is useless as a bound because the chain carries the full Mehler series (order 2 explicit, orders
  >= 3 at 1e-6 of order 2: frontier tests).

## 2. Where non-commutativity genuinely lives: the arrow algebra of a layer and the ReLU gauge

**The arrow algebra.**
- **Arrows.** The pair groupoid of the n neurons of a layer, with arrows (i, j). Its algebra is M_n, whose elements
  X_ij are pair (arrow) statistics: C_ij, K22_ij = k(a,a,b,b), K31_ab = k(a,a,a,b), D21.
- **State.** A positive unit measure d (the variances, or a kappa4 diagonal) defines omega_d(X) = sum_i d_i X_ii /
  sum d. Its modular operator acts on arrows by Delta e_ij = (d_i / d_j) e_ij, the Radon-Nikodym cocycle of the
  groupoid.
- **Kubo-Ando means.** They act as Schur multipliers X_ij -> m_f(d_i, d_j) X_ij. The KMS (geometric) one is
  X -> d^(1/2) X d^(1/2).
- **Covariance.** The covariance is the KMS sandwich of the correlation, C = D^(1/2) R D^(1/2) with D = diag(C).
  - The Gaussian relu update is Cov(y) = D^(1/2) K(a, R) D^(1/2), with K dimensionless (Mehler in R and the gauge
    invariants a_i = mu_i / sigma_i).

**The gauge.**
- **Action.** Lambda = diag(lambda_i) > 0 acts by z -> Lambda z, implemented on the network by W_(l-1) -> W_(l-1)
  Lambda and W_l -> Lambda^-1 W_l. This leaves F invariant, because relu(Lambda z) = Lambda relu(z).
- **Effect on the modular data.** It moves d -> Lambda^2 d. The modular cocycle changes by the coboundary
  (lambda_i / lambda_j)^2, which is exactly the brief's D' = D h(r) / h(s) with h = lambda^2.
- **Degrees.** A statistic has multidegree e if it scales as prod lambda_i^(e_i):
  - var_i has degree 2e_i, and C_ij has e_i + e_j;
  - K22_ab has 2e_a + 2e_b, and K31_ab has 3e_a + e_b;
  - kappa4_i has 4e_i;
  - R and a are invariant.

**Lemma 2.0.** The exact one-step map (relu, then W) commutes with the gauge, so every exactly computed statistic
transforms with its degree.

**Theorem 2.1 (the gauge selects the weighted geometric mean).**
- **Setting.** X is an arrow statistic of bidegree (p, q), u a unit statistic of degree r > 0, K a gauge-invariant
  arrow kernel and M: (0, inf)^2 -> R continuous.
- **Claim.** The closure X^_ij = M(u_i, u_j) K_ij is covariant under every positive diagonal gauge if and only if
  M(u, v) = M(1,1) u^(p/r) v^(q/r).
- **Symmetric case.** For a symmetric slice (p = q) closed from a unit statistic of degree r = 2p, the geometric mean
  sqrt(uv) is the only covariant Kubo-Ando mean. The arithmetic, harmonic, logarithmic and every other symmetric
  operator mean are excluded unless u is constant across neurons.

*Proof.* Covariance requires M(lambda^r u, mu^r v) = lambda^p mu^q M(u, v) for all lambda, mu > 0 (on a pair with
K_ij != 0). Put u = v = 1, lambda = s^(1/r) and mu = t^(1/r): then M(s, t) = s^(p/r) t^(q/r) M(1,1). The converse is
immediate. A Kubo-Ando mean with M(s,t) = sqrt(st) times a constant is the geometric mean, since m(1,1) = 1. ∎

**Corollary 2.2 (gauge spread is a truth-free error bound).** For a non-covariant closure and two gauges Lambda,
Lambda', map both estimates of the same physical slice to one frame. The truth is gauge invariant, so the larger of
the two errors is at least half their difference, entrywise. ∎

**Corollary 2.3 (size of the lever).** For a symmetric slice, arithmetic - geometric = (sqrt(u_a) - sqrt(u_b))^2 / 2,
whose mean relative size over pairs is cv(u)^2 / 4 + O(cv^3). Pre-activation unit data at width 1024 are concentrated:
the overlap spread of note XIII E8 is sqrt(2/PR), 3-18%. So the covariant and non-covariant closures differ by a few
percent of a slice whose route share in the output is 0.016 (note XXXIX section 5).

**Audit of the chain's closures** (est_v29, regen block lines 1610-1740, with the measured outcomes of notes XXXIX):

| closure | statistic (degree) | form | covariant | measured |
|---|---|---|---|---|
| default (2,2) slice | K22 (2,2) | (dG_a + dG_b) METRIC_C / 6, arithmetic | no | base |
| `V33_WK4M=3` | K22 | sqrt(g4_a g4_b) / 3, geometric | yes | adopted: -1.09% free-running, -0.23% renormalized |
| `V33_WK4M=1` | K22 | geometric + 2 s_a s_b C_ab^2 | yes | +5.19% free-running, -0.02% renormalized |
| lam core, diagonal | kappa4_i (4) | METRIC_C diag(W^T G W), G = diag(g) + lam C_off: dG_i = WW g_prev + lam s_off^2 has degree 2, and the constant METRIC_C = 2 stands in for the quenched metric \|w_i\|^2 | no | base |
| `V31_K4D=3` | kappa4_i | METRIC_C (WW g_prev + k4corr) + 6 g s_d^2 s_o^2 + 3 g s_o^4 | partly: the added classes have degree 4 | adopted (pair) |
| not run | kappa4_i | \|w_i\|^2 dG_i in place of 2 dG_i: the z-covariant completion of the isotropic core, free (n^2) | yes | predicted <= 0.5% raw: \|w_i\|^2 / 2 has rms spread sqrt(2/n) = 4.4% against a 30% kappa4 error; not proposed |
| lam core, K31 | K31 (3,1) | 0.5 METRIC_C lam C_off (degree (1,1)) | no | base |
| `V44_PMETRIC=1` | K31 | the same times var_c / mean(var) (stored transposed, degree (3,1)) | yes | +5.76% free-running, +5.71% renormalized: rejected |

Two covariant replacements helped and two hurt. Gauge covariance is a property of the exact map. The lam core is a
counterterm for omitted classes (note XXXIX section 6), and the output metric prefers its non-covariant shape. On the
true one-step statistics, by contrast, the covariant forms do explain more:
- the official networks' audit: C var_c 0.50-0.63 against C_off 0.49-0.58 of K31 (note XXXIX section 5);
- here, on small nets: section 6.3.

**The operator-valued sources are not states.** In the chaos picture of note XLI, neuron i carries the second-chaos
kernel H_i = sum_m P_im w2_m l_m l_m^T on the input space.
- It is a signed combination of positive rank-one births: w2 >= 0, but the mean-Jacobian leg P is signed. So H_i is
  indefinite and not a state. Lieb and Kubo-Ando concavity, which need positivity, do not act on it.
- The tuple (H_1, ..., H_n) in (M_n, tr) is the one genuinely non-commutative object in the problem. Its moments split
  in two:
  - trace-state moments tr(H_i H_j H_k) and tr H_i^4 are closed walks;
  - vector-state moments L_i^T H_j H_k L_l are open walks.
- Note XL priced the first at n^4 and reached the second through the Schur hub. So the non-commutative distribution of
  the second chaos is exactly the cost wall, and no Petz-type covariance reaches it more cheaply.

**Verdict on (b) in the Schur sense.** The geometric (KMS) mean is the right pairwise mean wherever an arrow statistic
is closed from unit statistics. The chain already uses it where the output reads it, and the remaining
non-covariant closures are counterterms whose shape the output metric chose.

## 3. (b) in the mixture sense: which mean of the covariance a readout sees

**Theorem 3.1 (power-mean hierarchy for radial mixtures).**
- **Setting.** z = sqrt(G) u, with u ~ N(mu~, C~) independent of G >= 0, and f positively homogeneous of degree k.
- **Claim.** E f(z) = E[G^(k/2)] E f(u) = E f(s_k u), where s_k = (E G^(k/2))^(1/k) is the k-th power mean of the
  scale.
- **k = 1 (means).** s_1 = E sqrt(G). The effective Gaussian N(s_1 mu~, s_1^2 C~) is the 2-Wasserstein barycenter
  of the components N(sqrt(g) mu~, g C~) weighted by the law of G. Its covariance is their Bures-Wasserstein mean,
  the fixed point of S = E (S^(1/2) g C~ S^(1/2))^(1/2), solved by S = s_1^2 C~.
- **k = 0 (gates).** For gates 1[z_i > 0], sign correlations and the arcsine kernel
  Pr(z_a > 0, z_b > 0) = 1/4 + arcsin(R_ab)/(2 pi) at zero mean, the law is invariant under radial mixing:
  E f(z) = E f(u). These readouts see the fibre's a~_i and correlations R~, which the barycentric Gaussian has exactly.
- **k = 2 (second moments).** s_2 = (E G)^(1/2): the arithmetic scale.
- **The moment-matched Gaussian is exact for none of k = 0, 1, 2.** It is
  N(s_1 mu~, s_2^2 C~ + (s_2^2 - s_1^2) mu~ mu~^T): the scale fluctuation becomes a rank-one Gaussian spike along the
  mean. This is the "spike is the gain" of note XIII E10.
- **Its bias at k = 0.** The spike biases a_i = mu_i / sigma_i and the correlations, so the bare first vertex
  Phi(mu/sigma) is not the radially invariant Pr(z > 0). This is exactly the "gauge inconsistency" note XL section 2
  named for the bare gate w1 = Phi of the legs and births.

*Proof.* f(sqrt(g) u) = g^(k/2) f(u), then independence. For the barycenter, the W2 barycenter of Gaussians has the
averaged mean, and its covariance solves the Bures-Wasserstein fixed point; S = s^2 C~ gives s = E sqrt(G). ∎

**Theorem 3.2 (certified one-sided barycenter bound).** For any finite Gaussian mixture law of a scalar,
sum p_k N(mu_k, sigma_k^2), with mu = sum p_k mu_k:
- G(mu, sum p_k sigma_k) <= E relu(z), where G(mu, s) = s phi(mu/s) + mu Phi(mu/s);
- equality holds if and only if the (mu_k, sigma_k) are collinear (a scale family);
- G(mu, sum p_k sigma_k) <= G(mu, sigma_mm) for the moment-matched sigma_mm.

*Proof.* G(mu, s) = E relu(mu + s xi) is jointly convex in (mu, s) and 1-homogeneous (it is the perspective of
phi(a) + a Phi(a)). Jensen gives sum p_k G(mu_k, sigma_k) >= G(sum p mu, sum p sigma), with equality if and only if the
points lie on a ray. G is increasing in s, and sigma_mm^2 = sum p (sigma_k^2 + mu_k^2) - mu^2 >= (sum p sigma_k)^2. ∎

Checked on 20,000 random mixtures: no violation, and equality to 1.8e-15 on scale families (section 6.2).

**Corollary 3.3 (sign of the moment-matched closure).**
- **Radial mixtures.** On every radial mixture the moment-matched Gaussian over-estimates the mean gate:
  - E relu(z) = G(mu, sigma_eff), with sigma_eff^2 = s_1^2 sigma~^2 = (s_1^2 sigma^2 - (s_2^2 - s_1^2) mu^2) / s_2^2,
    which is <= sigma^2;
  - to first order in g = Var G (E G = 1), this is sigma_eff^2 = sigma^2 - (g/4)(sigma^2 + mu^2), so
    delta E relu = -(g/8) sigma phi(a)(1 + a^2). That is exactly note XIII's identity.
- **In general** there is no sign. A location mixture is under-estimated: N(+-0.8, 0.36) gives 0.4254 against the
  closure's 0.3989. A scale mixture is over-estimated: 0.3853 against 0.3989.

**What it means for the chain.**
- The mean gate must see sigma_eff and the covariance must see sigma: two different "means" of one covariance.
- The chain does this through the Edgeworth kappa3/kappa4 terms. Note XIII E10 showed the first-order radial
  correction equals the kappa3 + kappa4 Edgeworth terms with kappa3 = 1.5 g mu sigma^2 and kappa4 = 3 g sigma^4.
  The chain's D3 carries g to 2%.
- The all-orders barycentric form differs from it at O(g^2), about 5e-5 relative and coherent, so the per-layer
  counterterms already span it.
- **Small nets (test E, section 6.3).** Against the per-neuron Gaussian-closure error of the mean gate (with true
  mu and sigma):
  - the 4th-order Edgeworth term (plus kappa3^2) with the true cumulants explains 96-99% at n = 64 and 85-99% at
    n = 32, with slope 0.93-1.04;
  - the barycentric correction explains 7-56% of the per-neuron variance but 69-128% of the mean.

  The barycentric correction is the coherent (radial) part, and the rest is quenched non-radial non-Gaussianity.
  **The mean-gate formula is not the error; the cumulants that feed it are.**

## 4. (a) Signed closure inequalities: certified versus sharp, and why they do not meet

**Theorem 4.1 (Gaussian slack).** Let b = (z, z o z) and M = Cov(b) = [[C, Dt], [Dt^T, V]] >= 0, with
Dt_ab = Cov(z_a, z_b^2) = D21_ba + 2 mu_b C_ab and V = Cov(z o z).
- **The constraint.** The certified constraint on the carried statistics is the Schur complement
  V >= Dt^T C^-1 Dt; together with the one-variable Hankel conditions it covers every positivity constraint of degree
  <= 2. On the diagonal, kappa4_b >= LB_b := [Dt^T C^-1 Dt]_bb - 2 var_b^2 - 4 mu_b kappa3_b - 4 mu_b^2 var_b.
- **(i)** For the Gaussian law with the same (mu, C): Dt = 2 C d(mu) and V = 2 C o C + 4 d(mu) C d(mu), so
  V - Dt^T C^-1 Dt = 2 C o C exactly. The diagonal slack is 2 var_b^2.
- **(ii)** For the true law, let gap_b = Var(z_b^2) - [Dt^T C^-1 Dt]_bb and r_b = kappa4_b / (2 var_b^2).
  - Regressing on z_b alone gives the upper bound gap_b <= 2 var_b^2 + kappa4_b - kappa3_b^2 / var_b.
  - The measured lower bound is gap_b >= 0.9 var_b^2 at every width and layer tested (section 6.3, column min
    gap/var^2).
  - So gap_b is of order 2 var_b^2, and gap_b / kappa4_b is about 1/r_b. A kappa4 closure that holds the other carried
    statistics at their true values violates the constraint only if its error exceeds gap_b.
- **The same slack covers the kappa3 slices.** An error delta D21 moves [Dt^T C^-1 Dt]_bb by about 4 mu_b delta
  D21_bb against a slack of 2 var_b^2. The constraint therefore needs a skewness error of order 1, against an actual
  skewness of about 1.5 g a, near 0.03.
- **(iii)** The Kubo-Ando/Petz inequalities of Theorem B add nothing. They are identities in the commutative case
  (Proposition 1.1).
- **(iv)** Operator Jensen for the gate is Proposition 1.3, whose slack (the Mehler tail) the chain carries.

*Proof.* (i) is the Gaussian fourth-moment identity. The upper bound in (ii) is
Var(z_b^2) - Cov(z_b, z_b^2)^2 / var_b = kappa4 + 2 var^2 - kappa3^2/var, since adding regressors only shrinks the
residual. The lower bound in (ii) is measured, not proved. ∎

**Measured** (sections 6.3 and 6.4). r falls as 1/n. At layer 8 (L = 8, seed 0), median r_b is 0.609, 0.388, 0.162
and 0.085 at n = 32, 64, 128 and 256. The median certified slack gap/kappa4 is 1.6, 2.4, 5.8 and 11.0. The truth
never violates the bound (0/8 layers x 4 widths), as a certificate must not.
- **At n = 1024.** r is about 0.03 at depth: note XIII measured excess kurtosis 0.066-0.086 for linear functionals of
  the last layer, so r = 0.033-0.043, and the scale-mixture law gives r_b = 1.5 g (1 + 2 a_b^2) with g = 0.022.
  - The certified constraint activates only for kappa4 errors 3x (|a| = 2) to 30x (a = 0) larger than kappa4 itself;
    the median is about 15x.
  - The chain's per-neuron kappa4 error is 30% (note XLI section 3).
  - So the constraint is never active.
- **Directional caveat (measured).** In near-null eigen-directions of C the ratio r reaches 10^2-10^5 at every width
  tried. These directions are rarely active units: tiny variance, sparse and heavy-tailed.
  - Matrix-level certified constraints can be active there.
  - Those directions carry negligible mean-channel weight. The chain drops saturated rows (a <= -2.5) exactly
    (note XXIX), and the K22 slice they would constrain has output route share 0.016.

**Theorem 4.2 (the second-chaos inequality, sharp).**
- **Setting.** Every generalized chi-square z = L.x + (x^T H x - tr H)/2 (H symmetric, any L).
- **Claim.** kappa4 kappa2 >= (4/3) kappa3^2. The constant is sharp: it is approached when |L| / ||H|| -> inf with L an
  eigenvector of H, and is not attained.
- **Path form** (always, by Cauchy-Schwarz): 12 |H L|^2 >= 12 (L^T H L)^2 / |L|^2 = (4/3) D3_path^2 / |L|^2, with
  D3_path = 3 L^T H L. It is the dominant form when |L| >> ||H||.

*Proof.* In an eigenbasis of H (eigenvalues h_k, components l_k of L):
- kappa2 = sum (h_k^2/2 + l_k^2), kappa3 = sum (h_k^3 + 3 h_k l_k^2) and kappa4 = sum (3 h_k^4 + 12 h_k^2 l_k^2)
  (note XLI's formula);
- write kappa3 = sum [(h_k / sqrt 2)(sqrt 2 h_k^2) + l_k (3 h_k l_k)];
- Cauchy-Schwarz gives kappa3^2 <= kappa2 sum (2 h_k^4 + 9 h_k^2 l_k^2) <= kappa2 (3/4) kappa4;
- for sharpness, the one-mode ratio is (1 + 3s)^2 / (3 (1/2 + s)(1 + 4s)) with s = l^2/h^2, which increases to 3/4. ∎

Checked numerically: the maximum of kappa3^2/(kappa2 kappa4) over 20,000 random instances is 0.750000.

**Measured on true laws** (n = 64, section 6.3).
- The inequality holds on 88-100% of neurons at layers >= 2 (n = 32: 78-100%). It is not a certificate: the third
  chaos (the c3 star, note XL) and higher cumulants break it.
- It is loose: the median kappa4 kappa2 / kappa3^2 is 2.3-11.
- As a "Rayleigh floor" on the chain's path class, 12 |H_i L_i|^2 >= (4/3) D3_i^2 / var_i is free (D3 and var are
  carried). Its ratio to the Schur hub measures how far L_i is from an eigenvector of H_i, nothing more.

**Proposition 4.3 (why signed bounds cannot reach a quenched residual).**
- **The projection.** If the truth lies in a certified convex set K_i, projecting the estimate onto K_i never
  increases the error. It improves only on the set A = {i : estimate outside K_i}, by at most sum over A of
  dist(estimate, K_i)^2.
- **The residual after counterterms.** The per-layer counterterms remove the coherent part of any sign-definite error.
  What remains is quenched, with near-zero mean per layer (notes XIII E7, XXXIX section 5).
- **Consequence.** A signed bound therefore pays only if it is active and its activity correlates with that quenched
  residual. By Theorem 4.1, A is empty for every carried per-neuron statistic at width 1024. ∎

**Verdict on (a).** Operator concavity (Lieb, Kubo-Ando) yields only identities in this commutative problem. Moment
positivity yields certified bounds that are inactive by a factor of 1/r, about 30. Sharp bounds exist only inside a
chaos truncation. The frame provides no sign-definite or two-sided closure bound that can move a kappa4 or kappa3
closure at width 1024.

## 5. (c) The trickle-down identity for a response bank

**Theorem 5.1 (input localization is the chaos decomposition).**
- **Setting.** Y_t = t X + B_t with X ~ N(0, I_n), and M_t(F) = E[F(X) | Y_t] for a bank F = (F_i) in L^2. The
  posterior is X | Y_t ~ N(Y_t/(1+t), I/(1+t)).
- **Claim.** Cov(M_t(F_i), M_t(F_j)) = sum_(k>=1) (t/(1+t))^k k! <f_k^i, f_k^j>. The trickle-down (AKV, or stage-8
  B with commuting observables) rate is
  d/dt Cov(M_t) = (1+t)^-2 E[J_t J_t^T], with J_t = E[grad F | Y_t].
- **Consequences.** The identity accounts for Cov(F) chaos by chaos. The target E F is the zeroth chaos, invariant
  along the localization (M_t is a martingale), so the chain's bias is not a term of it.

*Proof.* m_t = Y_t/(1+t) is Gaussian with Var(m_t) = Cov(m_t, X) = t/(1+t) per coordinate, so its correlation
with X is rho_t = sqrt(t/(1+t)). M_t = E[F(m_t + (1+t)^(-1/2) xi)] is the Ornstein-Uhlenbeck (Mehler) image of F at
rho_t, so its covariance is the Mehler series in rho_t^(2k). The rate follows from the innovation
form dM_t = Sigma_t E_t[grad F] dW with Sigma_t = I/(1+t) (Gaussian integration by parts). ∎

Checked for relu in one dimension at t = 0.1, 1 and 5: the direct and chaos sums agree to 2e-3 relative. The residue
is the chaos sum truncated at k = 60, which converges slowly for relu.

**Theorem 5.2 (the Schur hub is the t = 0 rate of a field localization).**
- **Setting.** Observe the layer field: Y_t = t z + C^(1/2) B_t, where z = z_l has any law nu and C = Cov(z). The
  posterior is a linear tilt of nu (AKV form), and for any bank a in L^2 Kushner-Stratonovich gives:
  - dE_t[a] = Cov_t(a, z)^T C^(-1/2) dW;
  - d<M(a_alpha), M(a_beta)>_t = Cov_t(a_alpha, z)^T C^-1 Cov_t(z, a_beta) dt;
  - Var(a) = E int_0^inf Cov_t(a, z)^T C^-1 Cov_t(a, z) dt.
- **Claim.** For note XLI's cross term C_i (first chaos x second chaos of z_i), Cov(C_i, z_j) = 2 Y_ij, so the adopted
  path class 12 [Y C^-1 Y^T]_ii is exactly 3 x (the t = 0 rate for C_i).
- **What the hub omits.** It is the part of Var(C_i) revealed only at t > 0. That part needs the posterior
  cross-covariances Cov_t(C_i, z), whose leading t-correction is a regression on the quadratic field {z_a z_b}: n^2
  regressors, the n^4 wall of notes XL and XLI.
- **The 4-cycle.** The 3 tr H^4 term sits in the bank element q_i^2 (q_i the second chaos of z_i). Already its t = 0
  rate Cov(q_i^2, z) is a closed walk, tr(H_i^2 H_j). So cheapness is a property of the bank element (open versus
  closed walk), not of t.

*Proof.* The filtering equation for a linear Gaussian observation of z with noise covariance C. Martingale
convergence gives Var(a) = E <M(a)>_inf, since the posterior concentrates at z. The identification with Y is note
XLI's Cov(C_i, z_j) = 2 Y_ij. ∎

For a Gaussian prior, Sigma_t = C/(1+t) and Cov_t(z_i^2, z) = 2 m_(t,i) Sigma_t e_i. The rate is
4 m_(t,i)^2 var/(1+t)^2, with E m_(t,i)^2 = mu^2 + var t/(1+t). Its integral is
4 mu^2 var + 4 var^2 int_0^inf t (1+t)^-3 dt = 4 mu^2 var + 2 var^2 = Var(z_i^2), as the theorem requires.

**Proposition 5.3 (exact bias identity: the integrated heat-equation defect).**
- **Setting.** Let Phi^(m, tau) be any estimator of E F(m + sqrt(tau) xi) (for the chain: the chain run on the input
  law N(m, tau I)), with Phi^(m, 0) = F(m) (exact on point masses) and enough regularity. Then
  E F - Phi^(0, 1) = int_0^1 E_(m ~ N(0, (1-tau) I)) [ (1/2) Delta_m Phi^ - d_tau Phi^ ](m, tau) dtau.
- **Homogeneous form.** With the chain's own homogeneity, Phi^(m, tau) = sqrt(tau) Phi^(m/sqrt(tau), 1), the
  integrand is the Euler-Stein defect (1/2)[Delta Phi^ - (Phi^ - m.grad Phi^)/tau].
  - For the exact map it vanishes identically (heat equation; E F = E Delta F for 1-homogeneous F).
  - At (m, tau) = (0, 1) it is the chain's mean minus the trace of its input Hessian, the Euler identity mu_i = tr H_i
    that note XLI measured at correlation 0.73-0.95 with the sources' H.
- **Proof.** Put h(tau) = E_(m ~ N(0,(1-tau)I)) Phi^(m, tau). Then h(0) = E F and h(1) = Phi^(0,1), and
  h'(tau) = E[d_tau Phi^] - (1/2) E[Delta Phi^]. ∎
- **Reading.** The identity is the stochastic-localization Dynkin formula specialised to the input's Gaussian bridge:
  the drift of the chain along the localization, the companion of the quadratic-variation identity 5.1.
- **Cost.** It is exact, and it is the only exact accounting of the mean's bias the frame supplies. Evaluating it
  needs the chain's Laplacian in n input directions along the bridge, about n chain runs per bridge point. And
  Hutchinson probes cannot resolve a defect that is a 1e-4 difference of two O(1) terms.

**Verdict on (c).**
- **Variance.** Yes, the trickle-down identity for a response bank gives an exact accounting: the chaos decomposition
  for the input, the hub plus a t > 0 remainder for the field.
- **Mean.** The identity is silent on the mean. The mean's bias has its own exact accounting (5.3), but not at the
  budget's price.

## 6. Numerical record

### 6.1 Theorems B and C (`f4code/f4_modular_checks.py`; d = 6, three observables, five random faithful states)

| check | worst over trials |
|---|---|
| K computed via eigen-sum vs tr(rho^(1/2) A rho^(1/2) B) directly | 1.8e-15 |
| m' - m - KZ | 1.8e-15 |
| E SJ(rho') - (SJ - KCK) (exact) | 5.3e-15 |
| min eig[(Cov^f - KCK) - E Cov^f(rho')], f = geo / har / log | >= 6.7e-3 / 4.5e-3 / 5.0e-3 (positive: the inequality holds strictly) |
| min eig(H - Cov(m_inf)) over 200 random pure-state decompositions per state | >= 0.61 (holds; random decompositions are far from Yu's optimum) |
| SJ >= K >= H (min eigenvalues of the differences) | > 0.16, > 0.12 |
| commutative case: max over f of |Cov^f - SJ| | 0.0 |

The positivity of 1 + <A - m, Z> was enforced (minimum eigenvalue 0.31-0.70). The stage-8 claims B and C hold as
stated. Their non-commutative content is the strict gaps, and those are 0 in the commutative case.

### 6.2 Theorems 3.2, 4.2, 5.1 (`f4code/f4_theorem_checks.py`)

```
Thm 4.2: max k3^2/(k2 k4) over 20000 random generalized chi-squares = 0.750000 (bound 0.75)
   MC cumulants k2 5.6294 k3 0.1767 k4 126.1853 vs formula 5.6210 0.2396 125.5917   (2e6 samples; k3 within ~1 s.e.)
Thm 3.2: violations 0/20000; max |truth - barycenter| on scale families 1.8e-15; min gap otherwise -4.4e-16
   moment-matched N(0,1) 0.3989: location mixture truth 0.4254 (closure under), scale mixture truth 0.3853 (closure over)
Thm 5.1: t=0.1: Var M_t direct 0.023391  chaos sum 0.023381
Thm 5.1: t=1.0: Var M_t direct 0.145491  chaos sum 0.145195
Thm 5.1: t=5.0: Var M_t direct 0.267777  chaos sum 0.267311
```

### 6.3 Exact-sample cumulant slices on small He networks (`f4code/f4_cum.py`, `f4code/f4_analyze.py`)

**Setup.**
- n = 32 and 64, L = 8, seeds 0 and 1, 4e6 inputs.
- Two passes over the same stream: exact means first, then central moments in float64. A half-sample split gives
  the K22 noise.
- Full outputs: `f4data/analyze_n32.txt` and `f4data/analyze_n64.txt`.

**Predictions** (made after a 16-wide smoke test, before these runs):
- geometric beats arithmetic for K22 at every layer >= 2;
- the covariant K31 shape beats C_ab, and the wrong-index shape is worst;
- the certified bound is never violated by the truth;
- Edgeworth with true cumulants explains >= 80% of the mean-gate error, and the barycentric correction explains the
  mean but not the per-neuron spread.

n = 64, seed 0 (seed 1 alike; layer 1 is Gaussian and omitted):

| layer | cv(var) | cv(k4) | median r | median gap/k4 | min gap/var^2 | frac k4 k2 >= (4/3) k3^2 | median ratio | K22 1-R^2: ari / geo / har / g var var | K31 1-R^2: 3 var_a C / C / var_b C |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 0.20 | 0.41 | 0.081 | 12.0 | 1.83 | 1.00 | 9.8 | 0.723 / 0.703 / 0.716 / 0.610 | 0.441 / 0.439 / 0.457 |
| 3 | 0.29 | 0.69 | 0.135 | 7.0 | 1.76 | 1.00 | 6.8 | 0.678 / 0.644 / 0.697 / 0.675 | 0.502 / 0.509 / 0.563 |
| 4 | 0.31 | 0.70 | 0.159 | 5.8 | 1.60 | 0.98 | 4.7 | 0.607 / 0.599 / 0.649 / 0.616 | 0.524 / 0.525 / 0.585 |
| 5 | 0.43 | 1.04 | 0.225 | 4.2 | 1.33 | 1.00 | 4.8 | 0.490 / 0.379 / 0.416 / 0.441 | 0.311 / 0.399 / 0.456 |
| 6 | 0.45 | 0.98 | 0.321 | 2.9 | 0.90 | 0.98 | 3.1 | 0.498 / 0.411 / 0.444 / 0.525 | 0.339 / 0.397 / 0.489 |
| 7 | 0.57 | 1.30 | 0.376 | 2.6 | 1.04 | 1.00 | 3.2 | 0.501 / 0.297 / 0.292 / 0.390 | 0.304 / 0.365 / 0.475 |
| 8 | 0.51 | 1.39 | 0.386 | 2.3 | 1.21 | 0.97 | 2.8 | 0.452 / 0.289 / 0.356 / 0.375 | 0.217 / 0.380 / 0.461 |

K22 noise (half-sample) is <= 0.002 of the slice's variance at layers >= 2.
- **K22.** The geometric mean is at least as good as the arithmetic on 28/28 (layer, network) cases at n = 32 and 64.
  The margin grows with cv(k4): at layer 8 it is 0.289 against 0.452 (seed 0) and 0.312 against 0.495 (seed 1). The
  harmonic mean is within a few points of the geometric. A single Schur mean explains only 25-75% of K22.
- **K31.** The covariant shape is best or tied on 13/14 cases at n = 64, and the wrong-index shape is worst on 14/14.
- **The chaos-2 inequality** fails on up to 12% of neurons (n = 64, seed 1, layer 8: 0.88) and up to 22% at n = 32.

Mean gate (test E). e_G = m_true - G(mu, sigma) per neuron; the Edgeworth term includes kappa3^2/72; the barycentric
correction uses sigma_eff^2 = sigma^2 - (g/4)(sigma^2 + mu^2), with g fitted from kappa3 ~ 1.5 g mu sigma^2:

| n, seed | Edgeworth: slope / variance explained, layers 2-8 | barycentric: variance explained / mean explained |
|---|---|---|
| 64, 0 | 0.94-1.04 / 0.973-0.993 | 0.10-0.40 / 0.81-1.28 |
| 64, 1 | 0.93-1.04 / 0.962-0.994 | 0.07-0.47 / 0.69-1.28 |
| 32, 0 | 0.85-1.08 / 0.880-0.970 | 0.10-0.50 / 0.69-1.35 |
| 32, 1 | 0.87-1.07 / 0.845-0.990 | 0.04-0.57 / 0.65-1.11 |

All predictions hold.

### 6.4 Width scaling of the Gaussian slack (`f4code/f4_light.py`; L = 8, seed 0, 1e6 inputs)

| layer | median r: n = 32 / 64 / 128 / 256 | median gap/k4: n = 32 / 64 / 128 / 256 | r of top eigen-direction (n = 256) | max r over eigen-directions (n = 256) |
|---|---|---|---|---|
| 2 | 0.135 / 0.080 / 0.040 / 0.020 | 6.9 / 12.2 / 24.2 / 49.7 | 0.020 | 0.030 |
| 4 | 0.332 / 0.161 / 0.093 / 0.046 | 2.8 / 5.8 / 10.3 / 21.4 | 0.074 | 6.9e3 |
| 6 | 0.505 / 0.317 / 0.142 / 0.066 | 1.8 / 2.9 / 6.8 / 14.6 | 0.100 | 1.3e4 |
| 8 | 0.609 / 0.388 / 0.162 / 0.085 | 1.6 / 2.4 / 5.8 / 11.0 | 0.106 | 1.3e4 |

- r halves per doubling of n from 64 on, and the certified slack doubles.
- The truth violates the certified bound on 0 neurons at every width and layer.
- Extrapolated to n = 1024 at depth 15: r of about 0.03 and slack of about 30, consistent with the official networks'
  kurtosis (note XIII E4).
- Near-null directions are wildly non-Gaussian: rarely active units.

### 6.5 The hub solve, priced and accuracy-checked (`f4code/f4_chol_test.py`, `f4data/chol_test.txt`)

M = C + eps mean(var) I from the Gaussian-closure covariance of a width-1024, depth-16 He network; Y = B C or
Y = B C^(1/2) with B Gaussian / n (the second weights the bottom directions, as note XLI says the hub's quadratic
form does). Units of 2n^3; relative rms error of diag(Y M^-1 Y^T) in float32 against float64:

```
layer 7 eps 0.001 cond 1.9e+04 Y=B C      | LU 1.334 units rel.err 3.8e-08 | chol 0.706 units rel.err 6.6e-07 finite True
layer 7 eps 0.001 cond 1.9e+04 Y=B C^1/2  | LU 1.334 units rel.err 1.3e-07 | chol 0.706 units rel.err 5.1e-07 finite True
layer 7 eps 0.01 cond 1.9e+03 Y=B C      | LU 1.334 units rel.err 3.8e-08 | chol 0.706 units rel.err 6.3e-07 finite True
layer 7 eps 0.01 cond 1.9e+03 Y=B C^1/2  | LU 1.334 units rel.err 5.0e-08 | chol 0.706 units rel.err 5.4e-07 finite True
layer 15 eps 0.001 cond 5.8e+04 Y=B C      | LU 1.334 units rel.err 4.1e-08 | chol 0.706 units rel.err 7.4e-07 finite True
layer 15 eps 0.001 cond 5.8e+04 Y=B C^1/2  | LU 1.334 units rel.err 2.2e-07 | chol 0.706 units rel.err 6.1e-07 finite True
layer 15 eps 0.01 cond 5.8e+03 Y=B C      | LU 1.334 units rel.err 4.1e-08 | chol 0.706 units rel.err 7.6e-07 finite True
layer 15 eps 0.01 cond 5.8e+03 Y=B C^1/2  | LU 1.334 units rel.err 7.0e-08 | chol 0.706 units rel.err 6.2e-07 finite True
```

### 6.6 The first vertex at k = 0 (`f4code/f4_gate.py`, `f4data/gate_n*.txt`; 2e6 inputs)

e = Pr(z_i > 0) - Phi(mu_i / sigma_i) per neuron. The table gives the fraction of its across-neuron variance
explained by the barycentric vertex Phi(mu / sigma_eff) (g fitted from kappa3 ~ 1.5 g mu var) and by the
Edgeworth-dressed vertex with the true kappa3, kappa4 and kappa3^2 terms. The "mean explained" ratios are omitted
because mean(e) is near 0 and the ratios are unstable.

| layer | rms e (n = 64, s0 / s1) | barycentric, s0 / s1 | Edgeworth, s0 / s1 |
|---|---|---|---|
| 2 | 6.2e-3 / 6.5e-3 | 0.00 / -0.05 | 0.956 / 0.966 |
| 3 | 8.7e-3 / 6.8e-3 | 0.44 / 0.30 | 0.941 / 0.931 |
| 4 | 8.9e-3 / 9.3e-3 | 0.55 / 0.67 | 0.908 / 0.915 |
| 5 | 1.1e-2 / 1.1e-2 | 0.38 / 0.63 | 0.836 / 0.902 |
| 6 | 1.5e-2 / 1.3e-2 | 0.49 / 0.50 | 0.840 / 0.764 |
| 7 | 1.6e-2 / 1.3e-2 | 0.46 / 0.68 | 0.587 / 0.853 |
| 8 | 1.6e-2 / 1.5e-2 | 0.57 / 0.65 | 0.657 / 0.877 |

n = 256, seed 0:

| layer | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| rms e | 1.4e-3 | 2.0e-3 | 2.6e-3 | 3.3e-3 | 3.3e-3 | 3.8e-3 | 4.0e-3 |
| barycentric share | -0.08 | 0.38 | 0.54 | 0.56 | 0.58 | 0.58 | 0.68 |
| Edgeworth share | 0.985 | 0.988 | 0.986 | 0.985 | 0.984 | 0.971 | 0.975 |

- **Scale of the defect.** The vertex defect falls about 4x from n = 64 to 256, as 1/n.
- **Radial share.** It is 55-68% at depth and holds or rises with width. That is consistent with note XIII's finding
  at n = 1024 that conditioning on the radius removes the non-Gaussianity of generic directions.
- **Edgeworth with exact cumulants.** It explains about 98% at n = 256. Whether the chain's own (inexact) cumulants
  do as well is what the unstable V53 run probed.

## 7. Estimator consequences

**No accuracy change follows from the frame.** Each candidate fails for the reason given:

| candidate | reason |
|---|---|
| covariant (KMS-mean) replacements of the remaining closures | the lam core's non-covariant shapes are output-metric counterterms; PMETRIC was measured +5.7% |
| barycentric mean gate (k = 1) | inside the Edgeworth terms the chain carries; the remainder is O(g^2) and coherent |
| barycentric first vertex (k = 0) | not excluded: free and stable, but predicted near-neutral; see `V58_BARYGATE` below |
| projection onto certified constraints | inactive by a factor of about 30 at width 1024 (Theorem 4.1) |
| the chaos-2 Rayleigh floor | not certified, and loose |
| Petz or Bures covariances in place of the arithmetic one | identical in the commutative problem |

**The one cost change: the hub as a Schur complement, computed by Cholesky (`V57_CHOL`, proposed).**
- **Current cost.** V56 computes diag(Y M^-1 Y^T), with M = C + eps mean(var) I, by `fnp.linalg.solve(M, Y^T)`
  (est_v29 line 1827). flopscope bills 2n^3/3 + 2n^3, which is 1.333 units per layer, 13 layers, 17.3 units.
- **The structure.** M is symmetric positive definite and only the quadratic form is needed. So
  diag(Y M^-1 Y^T) = squared column norms of L^-1 Y^T, with M = L L^T.
- **Cost of the replacement.**
  - `fnp.linalg.cholesky` is billed n^3/3: 0.167 units.
  - The forward substitution has no triangular primitive in flopscope, so it is done blocked with 64 x 64 blocks:
    16 small `inv` of diagonal blocks (0.004 units), 16 block products (0.06 units) and the strictly-lower updates,
    n^3 (1 - 1/16) (0.47 units).
  - Arrays are immutable, so the solved blocks are concatenated as they are produced. That is about 7.5 n^2 written
    elements, negligible.
  - The column norms are n^2.
  - Total about 0.70 units: saving 0.63 units x 13 = 8.2 units.
- **Measured on a realistic M** (`f4code/f4_chol_test.py`, section 6.5). M is the Gaussian-closure covariance of a
  width-1024 He network at layers 7 and 15, with the regularization of V56 (cond(M) = 1.9e4-5.8e4 at eps = 1e-3).
  - flopscope bills 0.706 units against 1.334 for the LU route.
  - The float32 relative error of the hub vector against a float64 reference is 5-8e-7 (LU: 4e-8-2e-7). The hub is
    a 30%-level kappa4 correction, so both are far below anything that matters.
  - No breakdown at cond 5.8e4.
- **Prediction** (scored regime, 100 networks, existing counterterms not refitted):
  - raw within +-0.3% per network;
  - C/B 0.2028 -> 0.1948;
  - adjusted 3.137e-9 -> 3.01e-9 (-3.9% +- 0.3);
  - about 80 more small ops per layer-solve (about 1000 per predict), so +<=0.03 s of grader-counted residual.
- **Risk.** At eps = 1e-3, cond(M) of the chain's own (non-Gaussian) covariance may exceed the Gaussian closure's
  5.8e4. The float32 Cholesky sufficient condition 1/(n eps_32), about 1.6e4, is already exceeded without failure, but
  a breakdown is possible on some network. It appears as a NaN on diag(L), checkable in n operations. The fallback per
  layer is eps = 1e-2 (-0.46 points raw, note XLI section 10) or the LU solve.
- **This is not a theory win.** It follows from reading the hub as a Schur complement/parallel sum, whose native
  factorization is Cholesky, and it is the only change this frame produced.

**One frame-derived accuracy candidate, cheap but expected near-neutral: the barycentric first vertex
(`V58_BARYGATE`).**
- **The idea.** Theorem 3.1 at k = 0 says the first vertex of every leg, birth and replicated slice should be the
  radially invariant Pr(z_i > 0) = Phi(a~_i), not the bare Phi(mu_i / sigma_i).
  - In the radial model, Phi(a~_i) = Phi(mu_i / sigma_eff,i), with sigma_eff^2 = sigma^2 - (g/4)(sigma^2 + mu^2).
  - g is read from the chain's own D3 by the K4SM fit, g = sum(v D3) / sum(v^2) with v = 1.5 mu var.
- **Why it differs from V53.** This is the bounded, all-orders-in-g resummation of the scale-mixture part of note XL's
  tadpole-dressed vertex (V53). V53 was neutral to adverse on 13 of 16 networks and unstable on 3 because of the
  truncated Edgeworth tails at large |a|. The barycentric form is a Phi of a shifted argument: positive, bounded, with
  no tail blow-up.
- **Small nets** (`f4code/f4_gate.py`; n = 64 with two seeds and n = 256, layers 3-8; section 6.6). It explains
  30-68% of the per-neuron variance of Pr(z > 0) - Phi(a), and 54-68% at n = 256 for layers 4-8. The Edgeworth-dressed vertex with the true cumulants explains 59-97% at n = 64
  (worst at depth, where the Edgeworth series of an indicator degrades) and about 98% at n = 256.
  - At the vertex (k = 0), then, the radial part is a large share of the defect, unlike at the mean gate (k = 1,
    7-56%).
- **Cost.** O(n) per layer: zero units.
- **Prediction** (cold, networks 0-15, paired):
  - raw change in [-3%, +2%];
  - no network worse than +5%, that is, no V53-type instability;
  - renormalized change in [-1%, +0.5%].
- **Decision.** Adopt only if the held-out renormalized change beats the base by more than two standard errors. This
  is listed to close the k = 0 branch, not because it is expected to win: the chain's per-layer counterterms already
  span much of a coherent vertex shift.

## 8. Minimal decisive experiments

**E-A: closes (a) at width 1024.**
- **Inputs.** Official networks 0 and 1, layers 2-15:
  - the existing 1.6e7-sample truth slices of notes XXXIX and XLI (mu, C, kappa3, kappa4 diagonal, D21). Where a
    dump lacks C or D21, regenerate them with one Monte Carlo pass of the `chaos2_diag.py` kind;
  - the chain's dumped g4row, D3, var (the `V21_NO_CONFINE=1` dumps of note XLI section 5).
- **Outputs, per layer.**
  - median r_b, and median and 5th-percentile gap_b / kappa4_b;
  - the activity of the certified bound: the fraction of active neurons with g4row_b < LB_b;
  - the activity of the chaos-2 floor F_b = (4/3) D3_b^2 / var_b: the fraction with g4row_b < F_b;
  - the change of the kappa4-diagonal residual's signal energy (truth - g4row, split-half noise removed) when g4row
    is clipped at F_b.
- **Compute.** One n x n solve per layer per network: under a minute on one instance.
- **Predictions.**
  - median r between 0.02 and 0.06 at layers 8-15;
  - median gap/kappa4 between 8 and 40;
  - certified activity 0.000 at every layer;
  - chaos-2 floor activity < 3%, and the clip changes the residual signal by < 1%.
- **Decision.**
  - Certified activity 0 and |clip effect| < 2%: the signed-inequality route is closed at width 1024.
  - The clip removes >= 5% of the residual signal on both networks: the cold 16-network run of a `V57_K4FLOOR` switch
    (O(n) per layer) follows.

**E-C (optional, closes the k = 0 branch): `V58_BARYGATE`.**
- **Implementation.** Where `V53_TADPOLE` dresses the first vertex of legs, births and replicated slices, use
  Phi(mu / sigma_eff) instead. sigma_eff^2 = max(sigma^2 - (g/4)(sigma^2 + mu^2), sigma^2 / 4), with g from the K4SM
  fit of the chain's own D3. Cost: O(n) per layer.
- **Runs.** Cold, networks 0-15, paired against the same batch's base; then renormalized on 0-99 by note XXXIX's
  protocol if the cold result is in range.
- **Predictions and decision.** As in section 7.

**E-B: the cost change.**
- **Implementation.** `V57_CHOL=1`: Cholesky plus 64-block forward substitution at est_v29 line 1827, with a NaN guard
  and the per-layer eps fallback.
- **Runs.** First cold, networks 0-15, paired against the same code with `V57_CHOL=0`; then scored on all 100
  networks with the existing counterterms.
- **Predictions.** As in section 7; the number of guard fallbacks is predicted to be 0 in 1300 layer solves.
- **Decision.** Adopt if adjusted improves by >= 3% on 100 networks, with raw within +-0.5% and the residual under the
  grader's 0.4 s gate.

## 9. What the frame says about a new system design

The natural home of KMS/Petz structure in this problem is small and precise:
- **The home.** The arrow algebra of each layer, with the ReLU gauge acting by modular-cocycle coboundaries.
- **Its content.** Compute in gauge-invariant coordinates (correlations, a_i), and close arrow statistics with
  (weighted) geometric means. The chain already does both where the output reads them.
- **Everything else.** The stage-8 machinery, applied to this problem, is commutative. Trickle-down is exact chaos
  bookkeeping. Petz and Bures covariances coincide with the classical one. The capacity bound is total variance.
- **Why it must be so.** The input law is a classical Gaussian: the hbar = 0 point of the quasi-free calculus, where
  every Petz covariance of a quantized observable coincides with the classical one, the differences being O(hbar^2).
  - Non-commutativity could re-enter only by a deliberate deformation (hbar > 0), with Berezin-Lieb brackets in the
    manner of stage-8 Theorem D.
  - Each side of such a bracket is a quantum expectation of the whole network, no cheaper than the classical one.
  - For the mean (phi = identity) the bracket collapses to identities.
- **Why signed inequalities cannot help.** The classes the chain still lacks are closed walks: the joint-gate
  triangle, the 4-cycle, the gate covariance. In the localization reading they are trace-state moments of the
  second-chaos tuple, reached either as t > 0 innovations of open-walk banks or as t = 0 rates of closed-walk banks.
  Either way their price is still n^4. Certified inequalities that might bound them are inactive by the 1/r of about 30 Gaussian slack,
  which grows with n.

So F4 does not support a new Phase 2 design. It supports three things:
- keeping the chaos-graded chain;
- the free cost item of section 7 (about -4% adjusted);
- treating any proposal that relies on Petz or KMS distinctions, or on operator-concavity bounds, as vocabulary unless
  it first names a non-commutative algebra that the network's variables do not generate.

## Referee report

Referee checks are in `ref_f4/` (`chk.py`, `bill.py`, `ovh.py`, `ovh2.py`, `blk.py`). Bottom line: **the report
survives as a correct negative result.** I found no fatal mathematical error. There are several overclaims and one
real cost-accounting omission (residual Python time). The single deliverable, V57_CHOL, is confirmed, but it does not
come from the frame.

### R1. Re-checked and correct
- **Thm 4.2.** The cumulant formula kappa_k = k![tr H^k/(2k) + L^T H^(k-2) L/2] holds for a non-diagonal 3x3 H:
  4e6-sample MC gives k2/k3/k4 = 6.346/12.196/93.54 against the formula's 6.340/12.204/93.51.
  - The adversarial maximum of kappa3^2/(kappa2 kappa4) over 20,000 instances is 0.7500000 (spread scales, d up to 5).
    So the constant 4/3 is right and sharp.
  - The proof is a single Cauchy-Schwarz. With L = 0 it gives the known pure-second-chaos form
    kappa3^2 <= (2/3) kappa2 kappa4. The constant 4/3 for the mixed first+second chaos is elementary.
  - The sharp limit is |L| >> ||H||, the near-Gaussian regime, where the inequality says nothing. This is worth
    stating.
- **Thm 4.1 (i).** The Gaussian slack V - Dt^T C^-1 Dt = 2 C o C holds to 2e-15.
  - The LB_b formula, Var(z^2) = kappa4 + 2var^2 + 4mu kappa3 + 4mu^2 var, is correct.
  - The upper bound gap <= kappa4 + 2var^2 - kappa3^2/var is correct.
- **Thm 3.2 and Cor. 3.3.** The location-mixture number 0.42544 against 0.39894 is reproduced.
  - The perspective-convexity proof is right: G = s g(mu/s) with g'' = phi > 0, so G is strictly convex across rays.
  - "Collinear" should read "on a common ray through the origin" (equal mu_k/sigma_k).
- **Thm 3.1.** It is correct, and it is one line: homogeneity plus independence. The W2-barycenter fixed point
  S = s1^2 C~ checks out. The barycenter reading is decoration: the identity E f(z) = E f(s_1 u) holds for every
  1-homogeneous f, which is stronger than any barycenter statement.
- **Thm 5.1.** The rate (1+t)^-2 E[J_t J_t^T] equals d/dt of the Mehler series term by term. I checked this
  analytically: E|J_t|^2 = sum_k k rho^(2(k-1)) k! |f_k|^2. It is the classical Ornstein-Uhlenbeck / heat-flow
  covariance identity (Houdre-Perez-Abreu, Chatterjee), not new.
- **Prop 5.3.** Correct: the fundamental theorem of calculus along the Gaussian bridge. The homogeneous form
  d_tau Phi = (Phi - m.grad Phi)/(2 tau) also checks.
- **Thm 2.1.** Correct (a Cauchy-type functional equation).
  - The 28/28 count of "geometric at least as good as arithmetic" is confirmed in `f4data/analyze_n*.txt`.
- **Cost.** I re-billed in flopscope 0.12.1 at n = 1024.
  - LU route: 1.3343 units. Cholesky + 64-block substitution: 0.7065 units. Float32 relative error 5.7e-7 at cond 4e3.
  - Saving 0.628 x 13 = 8.2 units, so C/B 0.2028 -> 0.1948 and adjusted x0.9606. The arithmetic is right.
  - The floor is n^3/3 + n^3 = 0.667 units; b = 32 reaches 0.683.

### R2. Errors and overclaims
1. **"The gate is a hard pinning, outside the same-algebra modular transport" is a category error.**
   - Absolute continuity of law(relu z) with respect to law(z) is irrelevant. That law is a pushforward, not a
     localization.
   - The localization a gate performs is conditioning on the sign pattern Q. Its density 1_Q/Pr(Q) with respect to
     law(z) exists: it is an indicator tilt, a coordinate pinning in the AKV sense, inside the same commutative
     algebra. Then E relu(z) = sum_Q Pr(Q) E[Q z | Q].
   - What it violates is faithfulness (the density vanishes on a set), not absolute continuity.
   - The conclusion (commutative, so Theorem B is AKV) is unaffected. The stated reason is wrong.
2. **"Covers every positivity constraint of degree <= 2" (Thm 4.1) is false as stated.**
   - The vector (z, z o z) contains only the diagonal quadratic monomials. The full degree-2 moment matrix on
     {1, z_a, z_a z_b} also constrains the K22 and K31 off-diagonal slices.
   - The Gaussian-slack argument very likely extends: the Isserlis slack is of order lambda_a lambda_b in eigen-pairs,
     and the corrections are of order r times that. But the extension is not proved, and it fails in the near-null
     directions the author already flags.
3. **Thm 5.2 is not new.**
   - Note XLI section 7 already states P4_i = 12[Y C^-1 Y^T]_ii = 3 Var(proj_(span z_l) C_i). The t = 0 innovation rate
     of a linear-Gaussian observation is, tautologically, the linear-regression variance. The theorem is a relabelling.
   - Its gloss "what the hub omits is revealed only at t > 0" is wrong in two ways:
     - **(a)** C_i is a function of the input x, not of z_l. So E∫rate = Var(E[C_i|z_l]), not Var(C_i).
     - **(b)** The hub's target is 12|H_i L_i|^2 = 3 Var(J_1^(x) C_i), the input-first-chaos part. The hub's dressed
       projection onto span z_l is neither a lower bound on that target nor nested with it: z_l has higher-chaos
       content, and note XLI's dressed hub beats the first-chaos one for exactly that reason.
   - So no positive "omitted remainder" exists in the sense claimed. The theorem also conflates a rate at t = 0 with
     an integral over t > 0.
4. **"NC could re-enter only by a deliberate hbar-deformation" (section 9) is overstated.**
   - The chain is a truncated moment problem. Its natural operator object is the tuple of compressed multiplication
     operators P_k z_a P_k on polynomials of degree <= k in L^2(law): multivariate Jacobi/Toeplitz operators. They do
     not commute, and they are exactly the Berezin-Toeplitz structure of stage-8 Theorem D.
   - A closure is a choice of commuting extension (Curto-Fialkow flat extension).
   - The author's Thm 4.1 is the degree-2 Hankel piece of this structure. The correct statement is: the
     non-commutativity exists, but every positivity or flatness constraint it supplies inherits the Gaussian slack
     about 1/r, so it is inactive at n = 1024.
5. **Smaller points.**
   - "The harmonic mean is within a few points of the geometric" is not always true: n = 32, seed 1, layer 8 gives
     0.244 against 0.385.
   - The covariant product shape g var_a var_b beats the geometric on 5/28 cases, and the non-covariant harmonic beats
     it on 2/28. So the data support "geometric beats arithmetic", not "covariance selects the geometric".
   - Section 6.5 says note XLI puts the hub's quadratic form in the bottom directions. Note XLI section 5 says the
     opposite: Y lies in C's top (outlier) subspace. The float32 risk is therefore smaller than the author's test
     suggests.

### R3. Cost-accounting omission: residual Python time
The author counts about 80 extra ops per layer and about 1000 per predict, for +0.03 s.
- **Op count.** The 64-block loop issues about 8 ops per block (two slices, a matmul, a subtract, a slice, an inv, a
  matmul, a concatenate) plus about 5 fixed. That is about 133 per layer, about 1700 per predict.
- **Measured overhead.** Here, a flopscope op on an 8x8 array costs 180-200 µs. The whole 16-block hub at n = 128,
  which is pure dispatch overhead, costs 68 ms per call, about 0.9 s over 13 layers on this machine.
- **Grader terms.** Note XXIX's local:grader ratio is about 0.3, which gives about +0.1-0.25 s against a 0.4 s gate.
  Note XLI section 10 already says the margin "must be checked on the grader's terms".
- **This is the binding risk for V57, not float32 conditioning.**
- **Fix.** Use b = 128: 0.745 units, about 69 ops per layer, saving 7.7 units, adjusted about -3.6%. Or a recursive
  2x2 split with 128-leaves. Run E-B with the residual measured on grader-equivalent hardware.
- **Wall time.** On this sandbox the Cholesky route was also slower than LU (7.6 s against 5.2 s) despite half the
  FLOPs, because the BLAS efficiency of small blocks is poor. Only billed FLOPs score, but the wall time should be
  watched.

### R4. Relevance at the 1e-8 level
- Nothing in F4 moves the raw final-layer MSE.
- V57_CHOL is a near-certain -3.6 to -3.9% on adjusted, if the residual gate holds. It is textbook linear algebra
  (an SPD quadratic form by Cholesky), already on note XLI's to-do list ("the solve is the next thing to make
  cheaper"). It is not a consequence of the frame.
- V58_BARYGATE: the vertex defect is rms 4e-3 at n = 256, falling as 1/n to about 1e-3 at 1024, and only 55-68% of it
  is radial. The per-layer counterterms already span coherent vertex shifts. I agree with "near-neutral", and I would
  rank E-C below E-B.
- E-A is cheap but its outcome is essentially determined by the r ~ 1/n argument. It is only worth running as a
  by-product of E-B.

### R5. Strongest surviving idea (corrected)
**The Gaussian-slack no-go**, stated correctly:
- **Claim.** At width n the layer law sits inside the moment cone at a distance of order var^2 from its boundary, in
  every non-degenerate direction. The non-Gaussian quantities a closure must supply are of order r var^2, with
  r = kappa4/(2var^2) ~ 1/n (measured 0.61/0.39/0.16/0.085 at n = 32-256, and about 0.03 at n = 1024 from the
  official kurtosis).
- **Consequence.** Every certified constraint is inactive unless the closure error exceeds about 1/r ~ 30 times the
  quantity being closed. That covers moment positivity, Hankel or flat-extension constraints, operator-convexity
  (Petz, Kubo-Ando) bounds, and SDP relaxations.
- **Status.** Proved for the diagonal Schur-complement constraint, conjectured for the full degree-2 moment matrix,
  and false in near-null directions of C, which the output does not read.
- **Value.** It closes the whole "signed or certified closure" branch for Phase 2 at once.
- **Engineering item.** V57_CHOL with b = 128 and a measured residual.

Verdict: it survives as a negative result. Priority 3/10, essentially all of it the cost item.
