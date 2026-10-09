# F6. Localization couplings for unbiased hybrid estimators: the central reduction and the kink floor

Frame F6 of the localization round. The question was whether a deterministic part (the chain, or a chain on conditional
laws) can be combined with a little sampling through a localization coupling, so that the sampling sees only a tiny
residual variance.

## 0. Verdict

The answer is no, and the result is a clean impossibility bound, sharp to about 10%. It also says why the deterministic
chains win.

1. **The natural home is abelian (section 2).** The input law is O(n)-invariant and F is positively homogeneous.
   - The radial leaves integrate exactly (Theorem 1).
   - The rotation part is a Gelfand pair (O(n), O(n-1)). Its commutant is abelian, so "Cartan = centre" holds in the
     sense of stage-8 Theorem A.
   - So every isotropic localization scheme is a *central* localization. The relevant schemes are Eldan's Gaussian
     channel, the OU/Mehler flow, isotropic noise, Haar-randomized quadrature designs and telescoping over
     localization times.
   - Each one acts on the estimator's MSE through a single measure: the Schoenberg-Godement spectral measure
     nu_F = sum_m a_m delta_m of F's rotation-averaged autocorrelation. This is the transverse measure of the problem
     (Theorem 2).
   - Noncommutativity never enters.
2. **nu_F has a universal heavy tail (section 3).** The tail comes from the ReLU kink. For the He-initialized network
   (infinite-width kernel), a_k ~ (sqrt2 / (2 pi^{3/2})) (L + (-1)^k K'_{L-1}(0)) k^{-5/2} = 0.127 (16 +- 0.015)
   k^{-5/2}. This is Theorem 3 (Bietti-Bach's depth-linear amplitude; the constant is ours and checked numerically).
   - The model reproduces the three official-network measurements that exist:
     - sigma^2: 0.0705 against 0.074;
     - odd-degree share: 54.6% against 54.5%, from the antithetic x0.91;
     - degree >= 4 even share: 30.2% against 30%, from the frame quadrature x0.60.
3. **Black-box floor (Theorems 5, 6).** No randomized quadrature with N forward passes can beat the following bounds
   (sigma^2 = 0.074, any N between 0.1 B and B):
   - signed weights: N * MSE >= 0.27 sigma^2;
   - positive weights: N * MSE >= 0.54 sigma^2, by a certified Delsarte-Yudin LP bound whose dual certificate is a
     cubic.
   - Orthonormal frames + antipodes (the leaders' Phase 1 family) achieve 0.60 sigma^2, within 7-11% of optimal.
   - No spherical 4-design exists below 525,824 nodes (8 B).
   - Adjusted floor: >= 3.1-3.4e-7 (signed) or >= 6.1-6.4e-7 (positive weights) for every budget fraction >= 0.1.
     That is 100x / 200x above our 3.14e-9 and 280x / 560x above the leaders.
4. **Every hybrid class examined inherits the floor (section 5)**, each by its own short argument:
   - chaos control variates with exactly known mean: degree ~750 would be needed;
   - response-directed / low-rank localization: capped by the first-chaos share, 19%;
   - Rao-Blackwellization along Eldan leaves: dominated (homogeneity transfers the hardness);
   - multilevel telescoping over localization times: a Minkowski no-gain theorem, MSE x cost >= c_min sigma^2;
   - Rao-Blackwellization along linear regions: leaves of angular radius 2e-4 carry 5e-7 of the variance;
   - internal-activation control variates: 4.6% sigma^2 kink residual;
   - symmetry-generated zero-variance directions: degrees <= 3 and the odd degrees only.
5. **The value of sampling at the frontier is <= 0.25% of raw at 0.1 B and <= 2.5% at B (Theorem 12; 0.46% and 4.8%
   with signed weights).**
   - Directional / per-neuron pooling would change this only if the chain's error lived in low-variance output
     directions.
   - On tiny networks it does not. The chain's error and the output variance are co-localized: corr(log e_i^2, log v_i)
     = 0.68-0.92. The oracle directional hybrid gains 0.2-0.7% at the official ratio (Proposition 13).

**System design implied:** no sampling component. The minimal decisive experiment (section 8, E1) is an offline audit:
16 networks, about 30 CPU-minutes. It either closes the frame on official networks or reopens it with a number.

---

## 1. What a hybrid has to achieve

**Notation.**
- F(x) = x_L(x) in R^n, with n = 1024 and L = 16.
- x ~ gamma_n = N(0, I_n).
- One forward pass costs 2^25 FLOPs = 1/64 unit, so B = 65536 passes and 0.1 B = 6554 passes.
- The chain's per-neuron bias energy is b^2 = raw = 1.55e-8; the output variance is sigma^2 = 0.074.

**The hybrid class.** Every estimator considered here has the form

    m_hat = D + sum_j w_j (F(x_j) - G(x_j)),      E[G] known (= D, or with a known bias),

with D deterministic, G a surrogate evaluable per node and the nodes x_j random. The class includes:
- plain MC (G = 0, D = 0);
- control variates;
- designs (correlated nodes);
- Rao-Blackwellized estimators (F replaced by a conditional expectation);
- multilevel telescoping (a sum of such terms, with G the previous level).

**What it must achieve.** If its variance is v_eff/N, it beats or matches the chain only when v_eff/N <~ b^2, that is

    v_eff <= N b^2 = 1.0e-4 = 1.4e-3 sigma^2          (N = 6554; raw 1e-8 needs 6.5e-5).

Equivalently, the coupling has to remove 99.86% of the per-sample variance. The rest of this note proves that no
localization coupling comes within two orders of magnitude of that.

## 2. The natural home: the input's localization structure is central

### 2.1 Radial leaves

**Theorem 1 (radial Rao-Blackwellization; proved).**
- Write x = r u with r = |x| and u = x/|x|. Under gamma_n, r ~ chi_n is independent of u ~ sigma (the uniform measure on
  the sphere).
- Positive homogeneity gives F(x) = r F(u), and hence E F = E[r] * integral F d sigma, with
  E r = sqrt2 Gamma((n+1)/2)/Gamma(n/2).
- Take any linear rule sum_j w_j F(x_j) whose radii are independent of its directions. Replacing each node by
  E[r] F(u_j) is the conditional expectation given the directions. It is unbiased and never increases the MSE.
- The variance it removes is Var(r) (E_u F)^2 / E[r^2] ~ m^2 / (2n), which is 0.64% of sigma^2 in the kernel model.
  The frontier tests (`notes/frontier-tests`) measured x0.993.

*Proof.* The rule's conditional expectation given (u_j) is sum_j w_j E[r] F(u_j), by linearity of F in r and
independence. Conditional expectation reduces variance (Rao-Blackwell). QED.

**The NCG reading.** This is the one place where "a measure on the space of leaves passes to the leaf space" holds
literally.
- The leaves are the rays, and F is linear along each leaf.
- The leafwise integral is exact (E[r] times the value at the transversal), and the transverse measure is sigma.
- From here on everything lives on S^{n-1}.

### 2.2 Rotations: a Gelfand pair, so the commutant is the centre

O(n) acts on L^2(S^{n-1}, sigma) by (pi(g)f)(u) = f(g^{-1}u). Write M = pi(O(n))''.
- L^2(S^{n-1}) = Ind_{O(n-1)}^{O(n)} 1 decomposes multiplicity-free into the harmonic spaces H_m, m = 0, 1, 2, ...
  This is because (O(n), O(n-1)) is a Gelfand pair.
- Hence M' is abelian. It is the algebra of zonal convolutions, which by Funk-Hecke is the algebra of multipliers on the
  degree m, so M' ~ l^inf(N).
- Since M' is abelian, M' is contained in M'' = M, and so Z(M) = M ∩ M' = M'.
- **The commutant is its own centre: Cartan = centre.**

**The state and its central measure.** The state relevant for estimating the mean of F is the rotation-averaged
autocovariance

    h_F(t) = integral_{O(n)} Cov_u( F(u), F(g u) ) dg restricted to <u, gu> = t,     h_F(t) = sum_{m>=1} a_m P_m(t),

with:
- P_m the normalized Gegenbauer polynomial in dimension n, P_m(1) = 1;
- a_m >= 0 the energy of F in H_m. This is Schoenberg's theorem, or Godement's Bochner theorem for Gelfand pairs. The
  per-neuron energies are averaged over the n outputs.

The vector state omega_F(pi(g)) = <F, pi(g) F> has GNS commutant l^inf({m : a_m > 0}). By stage-8 Theorem A, its
localization processes are exactly the positive martingales of functions of m, and its central (transverse) measure is

    nu_F = sum_m a_m delta_m        (the spectral measure of h_F).

**Theorem 2 (central reduction; proved).** Consider any isotropic scheme:
- (i) the posterior mean of an O(n)-equivariant observation channel (Eldan's Gaussian channel Y_t = tX + B_t, isotropic
  noise, the OU/Mehler flow);
- (ii) a Haar-randomized N-node linear rule;
- (iii) a telescoping sum of (i)-(ii).

Its variance is a linear functional of nu_F,

    Var = sum_{m>=1} a_m * mu_scheme(m),      mu_scheme(m) >= 0,

with a scheme-specific multiplier:
- (i) |multiplier|^2 of the channel on H_m;
- (ii) Q_m(X, w) = sum_{j,j'} w_j w_j' P_m(<x_j, x_j'>) (Theorem 4).

*Proof.*
- (i) An equivariant Markov operator commutes with pi, so it lies in M'. It therefore acts on H_m by a scalar (Schur's
  lemma, multiplicity one).
- (ii) Haar invariance gives Cov(F(Q x_j), F(Q x_j')) = h_F(<x_j, x_j'>). Expand h_F in P_m.
- (iii) Linearity. QED.

For the Gaussian channel this is exact in Hermite chaos.
- X | Y_t ~ N(Y_t/(1+t), I/(1+t)), so X = sqrt(lam) G + sqrt(1-lam) Z with lam = t/(1+t) and G, Z independent N(0, I).
  This checks the brief's formula.
- Hence M_t = E[F | Y_t] = (P_tau F)(G), with e^{-tau} = sqrt(lam).
- Therefore Var(M_t) = sum_k c_k lam^k, where c_k is the Hermite-chaos energy. Here c_k = a_k + O(1/n), and
  c_k = b_k + O(1/n) below.

**What this means.** For isotropic localization couplings the noncommutative structure is not merely unhelpful; it is
absent: M' is abelian. Every question in this frame is a linear program over the one measure nu_F.
- The schemes not covered by Theorem 2 are the quenched, non-isotropic ones: response-directed tilts, coordinate
  pinning and gate pinning.
- They are still N-node linear rules, so Theorems 5 and 6 bound them. Section 5 treats them individually.

**Checking the brief's claims.** All three are correct for this problem, and all reduce to classical statements because
the algebra is abelian:
- Cov(F_i, u.X | F_t) = u^T Sigma_t E[grad F_i | F_t] is Gaussian integration by parts.
- The trickle-down identity V_t - E V_{t+1} = C^* R C is classical here.
- "Conditioning is not freezing" is illustrated in section 5.5.

## 3. The spectrum nu_F of a deep He ReLU MLP

**The kernel model.** In the infinite-width (NNGP) model, the ensemble-averaged autocorrelation of one output on inputs
of norm sqrt(n) is

    K_L(rho) = f^{oL}(rho),   f(rho) = ( sqrt(1-rho^2) + (pi - arccos rho) rho ) / pi,   f(1) = 1,

so that E x_L^2 = 1 and E (E_x F_i)^2 = K_L(0). Its Taylor coefficients b_k are the ensemble Hermite energies, and
a_m = sum_j b_{m+2j} alpha^{(n)}_{m+2j,j} via the Gegenbauer connection coefficients (computed exactly in
`code/f6_spectrum.py`; a_m = b_m to 3 digits for m <= 8 at n = 1024).

**Theorem 3 (kink law; proved, with Delta-analyticity cited from Chen-Xu 2021).** Let c = 2 sqrt2/(3 pi) = 0.3001.

    K_L(1-eps) = 1 - eps + L c eps^{3/2} + o(eps^{3/2}),      K_L(-1+eps) = K_{L-1}(0) + c K'_{L-1}(0) eps^{3/2} + o(eps^{3/2}),

and therefore

    b_k = (3/(4 sqrt pi)) c (L + (-1)^k K'_{L-1}(0)) k^{-5/2} (1 + o(1)) = 0.1270 (L + (-1)^k K'_{L-1}(0)) k^{-5/2}.

For L = 16, K'_{15}(0) = 0.0152, and the tail is T_K = sum_{k>K} b_k ~ 0.0847 L K^{-3/2} = 1.355 K^{-3/2}.

*Proof.*
- f(1-eps) = 1 - eps + c eps^{3/2} + O(eps^{5/2}): expand arccos(1-eps) = sqrt(2 eps)(1 + eps/12 + ...) and
  sqrt(1-rho^2).
- The fixed point rho = 1 has f'(1) = 1, so the eps^{3/2} terms add along the composition: A_l = A_{l-1} + c.
- At -1, f(-1+eps) = f(1-eps) - (1-eps) = c eps^{3/2} + ... (because f(rho) - rho is even), so the outer composition
  contributes K'_{L-1}(0).
- The Flajolet-Odlyzko transfer gives [rho^k] (1-rho)^{3/2} ~ k^{-5/2}/Gamma(-3/2), with Gamma(-3/2) = 4 sqrt(pi)/3.
- Delta-analyticity of the composed arc-cosine kernel is Chen-Xu (2021), as used by Bietti-Bach (arXiv 2009.14397,
  Cor. 2, where the amplitude "grows linearly with L"). QED.

The ratio b_k / asymptote, computed: 0.25 (k = 10), 0.55 (30), 0.81 (100), 0.93 (300), 0.978 (1000), 0.998 (10^4).

**The numbers** (`data/f6_spec_L16.txt`), against the independent official-network measurements of `notes/frontier-tests`:

| quantity | kernel model (n = 1024, L = 16) | official networks |
|---|---|---|
| sigma^2 = 1 - K_L(0) | 0.0705 | 0.074 |
| E m_i^2 = K_L(0) | 0.929 | 0.944 (0.58^2 + 0.78^2) |
| radial share | 0.64% | 0.7% (x0.993) |
| odd-degree share | 54.6% | 54.5% (antithetic x0.91 per pass => 2 a_even = 0.91) |
| a_1, a_2, a_3 | 18.9%, 15.3%, 10.3% | |
| even degree >= 4 share | 30.2% | 30% (frames+antipodes x0.60 = 2 a_even>=4) |
| frames+antipodes N * Var | 0.0424 (x 0.074/0.0705 = 0.0445) | 0.0442 (7.2e-6 x 6144) |

The model reproduces every official measurement to 1-5%. I therefore use it for all floors, scaled by
0.074/0.0705 = 1.05 when quoting absolute numbers.

At finite width the shares are robust; the norm and the mid-order profile are not (`data/f6_check_*.txt`):
- the odd/even split is within 0.01 at n = 64-128;
- the absolute sigma^2, m^2 and V(lam) drift at L/n >= 0.06;
- n = 1024 (L/n = 0.016) matches, as the table shows.

**Proposition (localization profile and Poincaré deficit; proved, given Theorem 3).** Along Eldan's localization the
removed variance is V(lam) = sum c_k lam^k. Near the leaves (lam = 1 - eps), the unresolved variance is

    E Var(F | Y) = sum c_k (1 - (1-eps)^k) = eps K'(1) - L c eps^{3/2} + O(eps^2).

The first term is the Gaussian Poincaré bound. Equality would hold iff F were linear on the conditional laws. The
eps^{3/2} deficit is the kink tail.
- Computed: D(eps)/sigma^2 = 0.0122, 0.0914, 0.435 at eps = 1e-3, 1e-2, 0.1. The linear Poincaré term would give
  0.014, 0.14, 1.4.
- V(lam)/sigma^2 = 2.1%, 5.9%, 15.4%, 33.3%, 56.5%, 90.9% at lam = 0.1, 0.25, 0.5, 0.75, 0.9, 0.99.
- Variance is resolved only near the leaves. This is the quantitative content of "the kink makes localization slow".

## 4. The black-box floor

**Theorem 4 (exact design variance; proved).** Let a rule estimate integral F d sigma by sum_j w_j F(Q u_j), with
sum_j w_j = 1 and Q Haar-random. Then

    Var = sum_{j,j'} w_j w_j' h_F(<u_j, u_j'>) = sum_{m>=1} a_m Q_m,    Q_m = sum_{j,j'} w_j w_j' P_m(<u_j,u_j'>) = (1/dim H_m) sum_i |sum_j w_j Y_{m,i}(u_j)|^2 >= 0.

- For a fixed (non-random) rule, the MSE averaged over the network ensemble equals this. The reason is that the law of
  W_0 is rotation invariant, so F∘Q has the same law as F.
- Because the formula is linear in a_m, the average over the 100 official networks is the formula with the averaged
  spectrum.

**Lemma (no near-design for degrees >= 4 in budget; proved).**
- For even m >= 4, P_m >= -beta_m on [-1, 1], with beta_4 = 6/((n+4)(n-1)) = 5.70e-6 exactly. The minimum is at
  t^2 = 3/(n+4), from P_4 = ((n+2)(n+4) t^4 - 6(n+2) t^2 + 3)/((n-1)(n+1)).
- Numerically beta_6 = 9.5e-8, beta_8 = 2.7e-9 and beta_10 = 1.1e-10 (`f6_spectrum.py`).
- For any weights with sum w = 1:

      Q_m >= |w|^2 - beta_m (|w|_1^2 - |w|^2) >= |w|^2 (1 - (N-1) beta_m) >= (1 - (N-1) beta_m)/N.

**Theorem 5 (signed-weight floor; proved for any F given its spectrum; numbers from the validated model).**

    N * Var >= sum_{m even >= 4} a_m (1 - (N-1) beta_m)_+ .

Terms are dropped where beta_m is not computed; each dropped term is >= 0, so the bound stays valid.
- At n = 1024 and N = 6554 (0.1 B): N Var >= 0.0210 (0.298 sigma^2); scaled to sigma^2 = 0.074, 0.0220.
- At N = 65536 (B): 0.0191 (0.27 sigma^2); scaled, 0.0200.
- At most 3.7% (0.1 B) or 37% (B) of the degree-4 energy can be cancelled, and nothing of degree >= 6
  (N beta_6 = 6e-4).
- Consistent with Delsarte-Goethals-Seidel: an exact 4-design needs N >= C(n+1, 2) + n = 525,824 nodes (8 B).

**Theorem 6 (positive-weight floor: a certified Delsarte-Yudin LP; proved, certificate checked numerically).**
- For w >= 0, the pair measure mu = sum w_j w_j' delta_{<u_j,u_j'>} is a probability measure on [-1, 1] with
  mu({1}) >= 1/N and integral P_k d mu >= 0. So Var >= min integral h_F d mu over that LP.
- Dual certificate: any g = y_0 + sum_{k>=1} y_k P_k with y_k >= 0 and g <= h_F on [-1, 1] gives

      N * Var >= h_F(1) + (N-1) y_0 - sum_k y_k.

- `code/f6_lpdual.py` solves the LP (HiGHS, Delsarte constraints k <= 40, h_F truncated at degree 200, which is valid
  because the dropped terms are positive definite). It then checks g <= h_F on 2e6 points (violation <= 1.5e-8) and
  subtracts the violation.

| N | certified N * Var floor | / sigma^2 | certificate g | frames + antipodes |
|---|---|---|---|---|
| 2048 (= 2n) | 0.04214 | 0.597 | y1 P1 + y2 P2 + y5 P5 + y0 | 0.0424 (optimal: the LP support is {0, -1}, the cross-polytope) |
| 6554 (0.1 B) | 0.03964 | 0.562 | y1 = 0.01334, y2 = 0.01063, y3 = 0.00556, y0 = -1.6e-7 | 0.0424 (+7%) |
| 65536 (B) | 0.03818 | 0.541 | y1 = 0.01334, y2 = 0.01078, y3 = 0.00710, y0 = -1.2e-8 | 0.0424 (+11%) |

The certificate is a cubic below the network's correlation function.
- Degrees 1 and 2 are cancelled at full weight (y_1 = a_1, y_2 = a_2).
- Degree 3 is cancelled only partially (y_3 < a_3 = 0.0073): cancelling it needs near-antipodal pairs, and they double
  the even degrees.
- **The leaders' Phase 1 family is LP-optimal at N = 2n and within 7-11% of optimal at every budget.**

**Remark (information-based complexity; cited theorem in the GP model, conjecture beyond).**
- For a Gaussian prior on F, adaption and nonlinear algorithms do not help for linear functionals. The optimal method is
  the linear kriging rule, which may have signed weights (Traub-Wasilkowski-Woźniakowski, *Information-Based
  Complexity*, 1988, ch. 6).
- In the NNGP limit, Theorem 5 therefore bounds *every* black-box algorithm that uses N evaluations of F: adaptive
  designs, learned surrogates, Bayesian quadrature.
- For the actual (near-Gaussian, quenched) networks this is a conjecture.

**Corollary (adjusted floor).** Spending a fraction phi >= 0.1 of B on N = 65536 phi samples gives
adjusted = raw x phi >= (N Var)/65536:
- any linear rule: >= 3.1-3.4e-7;
- positive weights: >= 6.1-6.4e-7;
- frames + antipodes: 6.8e-7.

The deterministic system is at 3.14e-9 and the leaders at 1.1-2.1e-9. The empirical "the variance would have to fall
about 573x" of `notes/frontier-tests` is now a theorem: designs can reduce it by at most 1/0.54 = 1.85x with positive weights, or 3.7x with
signed weights.

## 5. Hybrids: each coupling inherits the floor

### 5.1 Control variates with exactly known mean (chaos truncations)

The functions whose means are known exactly, whatever the chain's accuracy, are:
- the mean-zero chaos polynomials G = sum_{1<=k<=K} <T_k, h_k(x)>, with E G = 0 for any tensors T_k (the chain only
  supplies good T_k, e.g. L = E grad F or the second-chaos sources H_i);
- transforms of F under symmetries of gamma_n, which are designs (Theorem 4).

**Theorem 7 (proved, given the spectrum).** By orthogonality of the Wiener chaoses, Var(F - G) >= T_K = sum_{k>K} c_k
for any degree-<=K control variate.
- With the best design on top, the effective residual is 2 sum_{k even > K} c_k ~ T_K.
- Reaching v_eff <= 6.5e-5 (raw 1e-8 at 0.1 B) needs **K >= 748** (T_K = 4.2e-5 at K = 1000, 2.5e-4 at K = 300).
- The chain-computable CVs are degree 1 (L x, n^2 per sample) and degree 2. The degree-2 CV is the sources, evaluated as
  sum_m P_im w2_m (l_m . x)^2 at about 2 passes per sample.
- Both lie inside degrees <= 3, which frames + antipodes already cancel for free. **Chain control variates add exactly
  nothing to the design.**

### 5.2 Response-directed and low-rank localization

**Proposition 8 (ensemble: proved via rotation invariance; quenched first-chaos bound: proved).**
- Observing a rank-r projection Pi x explains, in the ensemble,

      Var E[F | Pi x] = sum_k c_k * dim Sym^k(R^r)/dim Sym^k(R^n) ~ K_L(r/n) - K_L(0),

  that is 0.02% (r = 1), 1.2% (r = 64), 5.9% (r = 256) and 15% (r = 512) of sigma^2.
- For a quenched, response-directed choice, the linear (first-step) part is |Pi L_i|^2 <= |L_i|^2, the first-chaos
  share **18.9%**. That share is exactly the initial resolution rate in the innovation identity
  Cov(F_i, u.X) = u^T E grad F_i.
- `notes/frontier-tests` measured R^2 = 0.16 for linear regression on the top-256 layer-1 frame directions. That is first chaos only:
  85% of the 18.9%.
- So the brief's response-directed disintegration resolves at most about 19% of the variance in its linear step.
  Everything else sits in chaos degree >= 2, which an r-dimensional observation sees only as (r/n)^k.

### 5.3 Rao-Blackwellization along Eldan leaves

**Proposition 9 (variance part proved; dominance is a cost argument under the stated cost model).**
- Replacing F(X) by M_t(Y) = E[F | Y_t] (with lam = t/(1+t)) and placing the Y-nodes in a frames + antipodes design
  gives variance 2 sum_{k even >= 4} c_k lam^k / N.
- This is small for lam <= 0.25 (V4(0.25) = 2.7e-5). But M_t is the original problem again: by homogeneity,
  M_t(Y) = sqrt(1-lam) E[F(mu + Z)] with mu = sqrt(lam/(1-lam)) G.
  - At lam <= 0.25, 92% of first-layer units still sit within one standard deviation of their kink. The conditional law
    is no easier for a cumulant chain.
  - The chain's cost does not depend on the input mean, so each node costs a chain run.
- At >= 32 units per covariance-carrying chain run, 0.1 B buys 3 runs, which forces V(lam) <= 3e-8, i.e.
  lam <= 2e-6. That is the original law.
- At <= 0.05 units per run (a mean-field chain, about 3 passes), the conditional bias is O(sigma), because the
  inter-neuron covariance is about 0.93 of the second moment.
- **Dominated by the deterministic chain at the root** in both cases.

### 5.4 Multilevel telescoping over localization times

**Theorem 10 (Minkowski no-gain; proved).**
- Take levels D_0 (deterministic), D_1, ..., D_J = F (any approximations, e.g. chains on conditional laws at times
  t_1 < ... < t_J), coupled along one path. The level increments Delta_j = D_{j+1} - D_j are sampled with N_j independent
  paths at cost c_j each.
- Then Var x Cost >= (sum_j sqrt(V_j c_j))^2.
- Pathwise, sum_j Delta_j = F(X) - D_0, so sigma^2 = Var(sum_j Delta_j) <= (sum_j sqrt V_j)^2 by Minkowski.
- Hence

      MSE x Cost >= c_min sigma^2,       c_min = min_j c_j.

- Every level of a localization-time ladder needs a chain run per path, so c_min >= c_chain >= c_F. **No ladder beats
  plain MC.**
- In general only levels cheaper than a forward pass can help. Their variance is then that of F - (cheap surrogate),
  which leads back to 5.1, 5.2 and 5.6.

### 5.5 Rao-Blackwellization along linear regions (finest gate pinning)

**Proposition 11 (estimate).** F is exactly linear on each activation cone, so E[F | pattern] = F(E[X | pattern]). The
finest pinning makes F deterministic-linear. Its leaves are tiny:
- Each pre-activation has density 1/sqrt(4 pi) at 0, marginally over units and inputs, at every layer.
- A perturbation of angle delta moves it by E|dz| = 1.13 delta.
- So a random direction meets 16 n (0.282)(1.13) = 5215 delta gate flips. The linear regions have angular radius about
  1.9e-4.
- The within-leaf variance is about K'(1) delta^2 = **5e-7 sigma^2**. The leaf space (the transverse measure of the gate
  partition) carries all the variance and the leaves carry none.
- This is "conditioning is not freezing" in quantitative form: pinning the gates conditions F completely, but it removes
  nothing that a sampler pays for.

### 5.6 Internal-activation control variates

The forward pass yields every hidden activation for free, so G = Phi(a_i) z_{L,i} (the last pre-activation, with slope
Phi(a), a = mu/s) costs nothing.
- In the kernel model, a_i ~ N(0, K_15(0)/(1 - K_15(0))) = N(0, 3.46^2) and s^2 = 0.154 (`data/f6_lastcv.txt`).
- The kink residual is E_a[s^2 (Var relu - Phi^2)] = 3.2e-3 = **4.6% sigma^2**, which is 32x above v_eff.
- The bias is Phi(a)(b_{L-1} W)_i, i.e. the chain's error one layer earlier.
- Linearizing more layers only adds kink residual.

### 5.7 Symmetries and zero-variance directions

The symmetries of gamma_n are O(n); homogeneity adds the dilations as an equivariance of F.
- **Dilations:** exact, 0.64% (Theorem 1).
- **Antipodal map -I** (central in O(n)): acts by (-1)^m and cancels the odd degrees (54.6%) at the price of doubling the
  even degrees. Sign flips are not symmetries of ReLU, so this antithetic is the only exact odd cancellation.
- **Hyperoctahedral orbits** (frames): a 3-design, so degrees 1, 2 and 3 are cancelled.
- **Any finite orbit:** it kills degree m iff it has no invariant harmonic of degree <= m (Sobolev). Degree 4 needs
  525,824 nodes.
- **Within budget,** the symmetry-generated zero-variance subspace is exactly {degree <= 3} + {odd}: 69.9% of sigma^2,
  net x0.60 after the antipodal doubling. No other zero-variance direction exists, since a generic W_0 has trivial
  stabilizer.

### 5.8 The value of sampling at the frontier, and pooled (directional) hybrids

**Theorem 12 (combination bound; proved).**
- Let D have bias energy b^2 and U be independent and unbiased with variance v. The best combination D + c(U - D) has
  MSE b^2 v/(b^2 + v): a relative gain of b^2/(b^2 + v).
- With the certified positive-weight floor v = 0.0416/N, the gain is **0.24% at 0.1 B and 2.5% at B**.
- Signed weights give 0.46% and 4.8%.
- Our C/B would rise from 0.203 to 0.303 or 1.203, so the adjusted score gets worse in every case.

**Pooling across neurons.** The noise covariance of the sample mean is Sigma_out/N. The chain's error e is an n-vector.
Consider the directional (Wiener) combination m_hat = D + A(MC - D) with A = N eps^2 (N eps^2 I + Sigma_out)^{-1}.
- It uses the samples only in eigendirections with lambda_j < N eps^2, which is 1.4e-3 sigma^2 at 0.1 B.
- Its cost is one solve with the chain's own C_L (0.17 units) plus the passes.
- The output covariance is very ill-conditioned (`data/f6_covspec_*.txt`, n = 128-512, L = 16):
  - the top 64 directions hold 89-99.5% of the trace;
  - 9-11% of the directions are exactly zero (dead units);
  - another 5-11% of the live directions lie below the 0.1 B threshold.
- **If the chain's error were isotropic,** the Wiener hybrid would gain 12-13% (live subspace) at 0.1 B and 26-28% at B.

**Proposition 13 (co-localization; measured on tiny nets, mechanism sketched, conjecture for official nets).** The
error is not isotropic: it lives where the variance lives.

*Measured.* Deterministic part: the Gaussian-closure chain on the n = 64 nets (L = 8, 16) with 2e6-sample truth. The
official ratio N eps^2/sigma^2 was imposed (`data/f6_coloc.txt`).
- corr(log e_i^2, log v_i) over live units = 0.68-0.92.
- The 20% lowest-variance live units carry 0.1-4% of the mean error energy (and 0.01-1% of the variance).
- At the 0.1 B ratio, the computable Wiener hybrid is **adverse** (-2% to -6%), against a flat-error prediction of
  +12% to +36%.
- The *oracle* per-direction hybrid (an upper bound on any linear directional scheme) gains **0.2-0.7%** at 0.1 B and
  2-6% at B.

*Mechanism.* For unit i with standardized gate a:
- the readout's sensitivities to upstream errors are Phi(a) (mean) and phi(a)/(2s) (variance);
- its variance is s^2 v(a), with v(a) ~ phi(a)/|a|^3 as a -> -infinity;
- so e_i^2/v_i ~ phi(a)|a|^3 -> 0 on the low-variance (near-dead) units, and the ratio is largest on units near the
  kink, where the variance is large.

Errors and noise are produced at the same kinks and carried by the same linear maps.

*Test-time calibration as a pooled hybrid.* Fitting p per-network amplitudes from the samples costs p v/n per neuron
(OLS), i.e. 6.6e-9 p at 0.1 B, already 0.43 p of our raw. A GLS fit is far more precise but sees only the low-variance
directions, where the error is co-localized-small. Same verdict.

## 6. The estimator this frame implies, with cost

The design is **the null hybrid: no sampling component** in the Phase 2 system. Spend the budget on deterministic,
W-using structure. The chain escapes the floor because it computes the mean directly from W, through Stein / chaos
identities, instead of averaging F over the tail of nu_F.

For completeness, here is the best sampling component and its price, for a system with budget to spare below 0.1 B:

| component | cost (units of 2n^3) | effect |
|---|---|---|
| radial factor E[r] (exact) | 0 | always on; x0.994 variance |
| frames + antipodes, k frames | 32 k (2048 passes per frame) | N Var = 0.0445 (0.60 sigma^2), 1.07-1.11x the certified floor |
| Wiener combiner D + A(MC - D) | 0.17 (Cholesky of C_L + N eps^2 I) | predicted <= 0.5% raw at 0.1 B (Theorem 12 / Proposition 13) |
| chain control variates (degree 1, 2) | 1/16 pass + 2 passes per node | 0 on top of frames (Theorem 7) |

The leaders sit at C/B 0.11-0.15, so they have no free budget either. **Prediction: no entry can gain more than 0.5% of
raw from any sampling component at <= 0.1 B.**

## 7. Predictions (quantitative, falsifiable)

| ID | prediction |
|---|---|
| P-F6-1 | **Official-network localization profile** (4 nets): V(lam)/sigma^2 = 0.154, 0.333, 0.565, 0.909 at lam = 0.5, 0.75, 0.9, 0.99. D(eps)/sigma^2 = 0.0122, 0.0914, 0.435 at eps = 1e-3, 1e-2, 0.1. Parameter-free ratios: D(1e-2)/D(1e-3) = 7.47, D(0.1)/D(1e-2) = 4.76 (10 and 10 for a smooth F). All within 15%. |
| P-F6-2 | Any randomized quadrature with N in [6554, 65536] passes on official nets gives N * MSE >= 0.0416 with positive weights and >= 0.0200 with any weights. Frames + antipodes give 0.0445 +- 3%. |
| P-F6-3 | **Co-localization on official nets** (E1): corr(log e_i^2, log v_i) >= 0.6 for the adopted system; oracle directional hybrid gain <= 2% at N = 6554 and <= 6% at N = 65536; computable Wiener hybrid gain <= 0. |
| P-F6-4 | No hybrid of the classes in section 5 reaches raw < 1e-7 at <= 0.1 B. |

## 8. Minimal decisive experiments (for the lead on AWS)

### E1: co-localization audit (decides the frame; offline, about 30 CPU-minutes)

**Inputs.** Networks 0-15 (W_off{k}.npy, truth_off{k}.npz['m'][-1]), plus the adopted system's final-layer output D_k,
from the stored scored rows or a rerun of `estimator_final_v56.py`.

**Script outline** (`f6_coloc_official.py`):
1. Draw x with T = 2^17 (float32, batched 8192), forward pass, accumulate S = sum F F^T and sum F. Cost: 2^17 passes =
   2048 units = 4.4e12 FLOP per net, about 2 minutes numpy on 32 cores. The chain's C_L may be dumped instead.
2. Form Sigma_out and its live mask (diag > 1e-9 sigma^2), then the eigendecomposition of the live block.
3. Form e = D - m* and e_j = Q^T e.

**Outputs per net.**
- sigma^2;
- corr(log e_i^2, log v_i);
- the share of |e|^2 in eigendirections with lambda_j < tau, for tau/sigma^2 in {3.4e-4, 1.37e-3, 1.37e-2};
- the oracle directional gain 1 - sum_j e_j^2 (lambda_j/N)/(e_j^2 + lambda_j/N) / |e|^2;
- the Wiener gain with eps^2 fitted on nets 0-7 and judged on 8-15, at N in {1638, 6554, 65536}.

The truth noise (1e9 samples, 7e-11) and the covariance noise (n/T = 0.8%) are negligible.

**Decision rule.**
- If the mean oracle gain at N = 6554 is < 5%, F6 is closed: no sampling component, P-F6-3 confirmed.
- If the held-out Wiener gain at N = 6554 exceeds 33% raw, the break-even for our C/B 0.203 -> 0.303, implement the
  combiner (0.17 units + 0.1 B) and run it in the scored regime.
- In between, report the number. A sampling component is then worth it only to an entry below C/B 0.1 with free budget.

### E2: spectral-model confirmation at mid and high order (optional; about 20 minutes on one machine)

E2 tests the kink law, which every floor depends on, beyond the two points (antithetic, frames) already measured.

**Inputs.** Networks 0-3; T = 2^16 pairs (x, x') per net.

**Measurements.** F(x) and F(rho x + sqrt(1-rho^2) x') for rho in {0.5, 0.75, 0.9, 0.99, 1-1e-2, 1-1e-3} (and 0.1 as a
redundant check).
- Record V(rho) = Cov(F(x), F(x_rho)) and D(eps) = E|F(x) - F(x_{1-eps})|^2/2, per neuron then averaged.
- Paired differences make D(eps) low-noise.

**Cost.** 4 x 2^16 x 7 passes = 1.8e6 passes = 6e13 FLOP.

**Decision.** If P-F6-1 holds within 15%, the kink law and every floor in sections 4-5 hold on official nets as stated.
If D(eps) is closer to linear than predicted, F is smoother than the kernel model and Theorems 5-7 must be recomputed
with the measured spectrum (the theorems hold for any spectrum; only the numbers change).

## 9. Checks run (local, all in `loc/code`, outputs in `loc/data`)

**`f6_spectrum.py 16` → `f6_spec_L16.txt` (18 s).** FFT of f^{o16} on the unit circle, 2^18 points:
- sum b_k = 1.000000 and b_k >= 0;
- b_0 = 0.92946 (direct iteration 0.92946);
- the table of section 3;
- the asymptotic ratios;
- T_K and K* = 748;
- V(lam), the rank-r shares, the Gegenbauer a_m (equal to b_m to 3 digits) and beta_k (beta_4 equal to the exact
  formula).

**`f6_lp.py`, `f6_lpdual.py` → `f6_lp.txt`, `f6_lpdual.txt` (6 s).** The primal LP and the certified dual of
Theorem 6. The violation is 1.5e-8 at N = 2048 and <= 2.5e-12 otherwise.

**`f6_check.py` → `f6_check_n64L16.txt`, `f6_check_n128L16.txt`, `f6_check_n64L8.txt`.** Quenched random nets:
- even share 0.449 / 0.453 / 0.450 / 0.432, against kernel 0.453 / 0.453 / 0.454 / 0.422;
- the frame variance (±, 2n nodes) matches the truth-referenced value (1.756e-4 against 1.753e-4), which checks the
  design identity;
- the absolute sigma^2 and V(lam) deviate at small width as noted.

**`f6_covspec.py` → `f6_covspec_n{128,256,512}.txt`.** The output-covariance spectra, live-subspace thresholds and
flat-error Wiener gains of section 5.8.

**`f6_coloc.py` → `f6_coloc.txt`.** Proposition 13.

**`f6_lastcv.py` → `f6_lastcv.txt`.** Section 5.6. It also checks E_a[s^2 v(a)] = 0.0705 = sigma^2, so the per-neuron
Gaussian picture is consistent with the kernel.

## 10. What this says beyond F6

- **The frame is closed at the level of a theorem.** The obstruction is not the coupling but the measure: nu_F has a
  k^{-5/2} tail with amplitude 0.127 L, created by the ReLU kink at every layer and summed along depth (Theorem 3).
  Every isotropic localization, design or telescoping acts on nu_F through a nonnegative multiplier, so the tail is paid
  in full (Theorems 2, 4-6, 10).
- **The NCG frame is the right home, and it is degenerate.**
  - The commutant of the input symmetry is abelian, so Cartan = centre and the transverse measure is all of the
    localization (stage-8 Theorem A in its degenerate case).
  - It delivers a sharp floor (an LP over the Gelfand pair's spherical functions), not an algorithmic unlock.
  - Where noncommutativity would bite (quenched, non-isotropic schemes), the general N-node bounds and the
    co-localization measurement close the door.
- **Why deterministic chains win.** A W-reading estimator never integrates over the transverse spectrum. It computes
  the mean through Stein / chaos identities layer by layer: E F = E[Delta F] = tr H from homogeneity (note XLI),
  through the sources.
  - Its error is set by how well it carries non-Gaussianity: the same eps^{3/2} kink content, but deterministically.
  - The useful lesson for the other frames: the kink amplitude L c is the one number that measures how far the network
    is from a smooth Gaussian functional. Any deterministic closure works against exactly that budget.

## Referee report

Adversarial referee, frame F6. My checks are in `loc/ref_f6/` (`chk_radial.py`, `chk_signed.py`, `chk_cv1.py`,
`kern_cv.py`). Their outputs are quoted below.

### Summary verdict

The note's central quantitative result survives. The sharp part is a certified Delsarte-Yudin LP floor for
W-oblivious (or Haar-randomized), positive-weight N-node rules. Combined with the value-of-sampling bound
b^2/(b^2+v), it shows that adding sampling to the chain is worth at most about 0.25% raw at 0.1 B.

The surrounding claims need correction:
- The framing ("natural home", "clean impossibility for hybrids") is substantially overstated.
- One proof (Theorem 5, signed weights) is invalid as written.
- One premise (Theorem 7: which control variates have known means) is false.
- Theorem 1 has a minor formula error.

The practical conclusion is "no sampling component". That is what the system already does, so the note changes
nothing at the 1e-8 level. It is a closure argument, not an improvement.

### Mathematical errors

**1. Theorem 1, the removed variance is misstated.**
- The removed variance is Var(r) E|F(u)|^2, not Var(r)(E_u F)^2/E r^2 ~ m^2/(2n).
- Derivation: Var(rF(u)) - Var(E[r] F(u)) = E r^2 E F^2 - (E r)^2 E F^2 = Var(r) E[F(u)^2].
- With E|x_L|^2 = 1 per neuron this is 1/(2n), i.e. 0.69% of sigma^2, not 0.64%. That is in better agreement with
  the measured x0.993.
- Tiny-net check (4e5 samples):

  | net | measured removal | Var(r) E\|F(u)\|^2 | author's formula |
  |---|---|---|---|
  | n = 32, L = 8 | 5.97% | 6.01% | 4.53% |
  | n = 64, L = 16 | 4.03% | 4.00% | 3.24% |

- The statement "never increases MSE" is fine.

**2. Theorem 5 (signed weights): the proof is invalid.**
- The lemma's first inequality, Q_m >= |w|^2 - beta_m(|w|_1^2 - |w|^2), bounds every off-diagonal term
  w_j w_j' P_m(t) below by -beta_m |w_j w_j'|. That holds only when w_j w_j' >= 0. For a mixed-sign pair one needs
  P_m <= beta_m, which is false.
- Counterexample: two coincident nodes with w = (2, -1). Then Q_m = 1, while the right-hand side is 5 - 4 beta_m.
- The final inequality, Q_m >= (1 - (N-1) beta_m)/N, also does not follow from "off-diagonal Gram entries >= -beta_m"
  alone. Three coplanar unit feature vectors with positive mutual inner products have 0 in their affine hull, so the
  abstract minimum is 0.
- Whether it holds for actual degree-m zonal features is open.
- Numerics (`chk_signed.py`, n = 6, N = 8, 30 L-BFGS restarts over node sets, optimal signed kriging weights):
  - single degree: min Q_4 = 0.033 against the claimed 0.020, so no violation was found;
  - full spectrum a_m = m^{-5/2}, m <= 20: achieved N Var = 0.189 against the claimed floor 0.013.
- So the signed-weight bound is plausible but unproved. It should be relabelled a conjecture (numerically supported at
  tiny n).
- The adjusted "floor >= 3.1e-7 (any weights)" is likewise not a theorem.
- The 525,824-node statement is fine: the p^2-vanishing argument for exact degree-4 cubature does not use the sign of
  the weights.
- Theorem 6 (positive weights) is correct: I re-derived the dual inequality, and it requires g(1) <= h(1), which
  g <= h gives. It is the bound that should carry the note.

**3. Theorem 3, a typo.** "f(rho) - rho is even" should read "f(rho) - rho/2 is even", since f(-rho) = f(rho) - rho.
- The consequence used, f(-1+eps) = f(1-eps) - (1-eps) = c eps^{3/2}, is correct.
- I re-derived c = 2 sqrt2/(3 pi), the additivity A_L = L c, the -1 singularity amplitude c K'_{L-1}(0) and the
  prefactor 3/(4 sqrt pi). All are correct.

**4. Theorem 7, the premise is false.**
- The note claims that "the functions whose means are known exactly ... are the mean-zero chaos polynomials" and
  symmetry transforms. That is false.
- Every first-layer ridge feature relu(x . w_j) has the exact mean |w_j|/sqrt(2 pi) and infinite chaos degree: it
  carries the layer-1 kink. More generally, any function of a few projections with a closed-form low-dimensional
  Gaussian integral qualifies.
- So "a known-mean CV needs degree 748" is a non sequitur. The orthogonality statement Var(F-G) >= T_K is true only
  for polynomial G.
- Checks:
  - Tiny net, n = 64, L = 16, out-of-sample linear regression (`chk_cv1.py`): the residual share is 0.636 with x
    alone and 0.428 with (x, x_1 - E x_1). The known-mean layer-1 features explain about 21% more of the variance than
    first chaos does.
  - Kernel model, n = 1024 (`kern_cv.py`): a linear functional of x_1 explains 25.7% of sigma^2, against 18.9% for
    x_0. For x_l with the chain-supplied mean, it explains 31.7% (l = 2), 62.5% (l = 8) and 95.4% (l = 15). The
    l = 15 figure reproduces the note's 4.6% residual.
- The conclusion survives, because the best known-mean residual is still about 74% of sigma^2, far above the 0.14%
  needed. The argument has to be replaced, though: "known mean" is not "finite chaos".

**5. Scope error in sections 2.2 and 4.**
- Theorem 4's identification of fixed-rule MSE with the Haar formula requires the node set to be independent of W
  (law of W_0 rotation invariant *given* the design).
- The sentence "they [quenched, non-isotropic schemes] are still N-node linear rules, so Theorems 5 and 6 bound them"
  is wrong. Response-directed tilts, stratification on first-layer gates, W-adapted importance sampling and any node
  placement computed from W are not covered.
- The IBC remark does not repair this. In this problem the algorithm reads W, so F is known, not a sample from a GP
  prior. "Black-box" is a restriction on the information model, and it is exactly the restriction the deterministic
  chain violates.

**6. Proposition 8 is stronger than its proof.**
- The (r/n)^k law is an ensemble statement.
- For a quenched, response-directed choice of directions, the second chaos need not be spread out. Note XLI writes
  H_i as a sum of low-rank legs, and a projection onto the legs' span could capture a large part of the degree-2
  energy of a given neuron.
- Only the first-chaos cap (18.9%) is proved. The rest should be marked as a sketch.

**7. Theorem 10 is near-tautological.**
- It is correct once one assumes that c_min >= the cost of a forward pass and that D_0 is a constant.
- The interesting case is the one it excludes: D_0 a cheap random surrogate with exactly known mean. That reduces to
  item 4.

**8. Minor points.**
- The certificate at N = 2048 has degree 5 (y_5 > 0, y_2 < a_2), not 3; "the certificate is cubic" holds at
  N = 6554 and 65536.
- D(1e-2)/sigma^2 = 0.0914 is the exact FFT value. The two-term kink asymptotic gives 0.074. That is fine, but it shows
  that P-F6-1's ratios rest on the exact spectrum, not on Theorem 3's leading term.

### Overclaims

- **"The natural home is abelian / Cartan = centre."** This is true for the Haar-averaged (W-blind) sector, where it
  is an immediate consequence of multiplicity-freeness. It is not a structural fact about the estimation problem,
  whose relevant data is W, not the input symmetry. It relabels W-oblivious Monte Carlo, which was already known to be
  about 300x off the frontier (pure MC adjusted 1.1e-6), and says nothing about the deterministic chain.
- **"A clean impossibility bound for localization-coupled hybrids, sharp to about 10%."** Rigorous and sharp only for
  positive-weight, W-oblivious designs. The hybrid classes in section 5 are closed case by case, at varying levels of
  rigor: Propositions 8, 11 and 13 are sketches or measurements, Theorem 7's premise is false, and Theorem 10's cost
  premise is assumed.
- **Proposition 13 (co-localization).** Measured with a Gaussian-closure chain on n = 64 nets, not the calibrated
  system. Calibrated counterterms could move the residual error into other directions, so the claim really is open
  until E1.
- **"The constant is ours."** The k^{-5/2} decay and the linear-in-L amplitude are Bietti-Bach (2009.14397) and
  Chen-Xu (2021). Only the explicit prefactor and the (-1)^k K'_{L-1}(0) parity correction are added.

### Cost accounting

- Frames: 2048 passes = 32 units; 0.1 B = 3.2 frames. Correct.
- Floors (N Var/65536) and break-even (33% raw for C/B 0.203 -> 0.303): correct.
- Wiener combiner: about n^3/3 = 0.17 units. Correct if flopscope bills Cholesky at n^3/3; at worst under 1 unit.
- Omitted but negligible: Haar frame generation (QR, about 0.67 units per frame).
- E1: 4.4e12 FLOP per network for 2^17 passes, correct; "2 minutes on 32 cores" is conservative.
- E1 also needs D_k from the adopted system, i.e. a rerun or stored outputs. That is cheap.
- E2: 6e13 FLOP. Correct.

### Relevance at the 1e-8 level

- None directly. The recommendation is the status quo: no sampling.
- The real value is negative: it saves effort on sampling hybrids, with a certified positive-weight number.
- E1 is cheap and well posed, and its decision rule is correct. I would run it only as a sanity audit of dead and
  near-dead units, where any chain error would be fixable for free by one sample.

### Novelty

- Delsarte-Yudin LP bounds: classical (Delsarte-Goethals-Seidel 1977; Yudin 1993; Cohn-Kumar 2007).
- Randomized-design variance as kernel energy: classical (Hickernell's discrepancy; Brauchart-Saff-Sloan-Womersley,
  QMC designs).
- Kink spectrum: Bietti-Bach and Chen-Xu.
- New here: assembling these for this problem, with an exact certificate and validation against three official-network
  measurements. That is competent applied work, not new theory.

### Strongest surviving idea (corrected)

1. For any positive-weight N-node rule whose nodes do not depend on W (or are Haar-randomized),
   N * MSE >= h_F(1) + (N-1) y_0 - sum y_k. With the kernel spectrum this is 0.56 sigma^2 at 0.1 B and 0.54 sigma^2
   at B. The certificate is cubic at those N, and frames + antipodes are within 7-11% of it.
2. So an independent unbiased add-on improves the chain by at most b^2/(b^2+v), about 0.24% raw at 0.1 B.
3. Corrections:
   - the signed-weight version is a conjecture;
   - the radial share is Var(r) E|F(u)|^2 = 0.69%;
   - known-mean control variates include first-layer ridge features (and lower-dimensional Gaussian integrals), which
     raise the explainable share from 18.9% to about 26% in the kernel model. That is still nowhere near the 99.86%
     needed.
