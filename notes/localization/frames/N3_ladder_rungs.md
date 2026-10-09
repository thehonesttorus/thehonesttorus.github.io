# N3. The Euler-Stein ladder as operator-valued moment constraints: a defect hierarchy that sees closed walks

Frame N3 of the NCG round. Status labels:
- **proved**: a complete proof is given here;
- **sketch**: the argument is given, the routine details are not;
- **measured**: numerical, small n, outputs quoted;
- **conjecture**.

Code is in `loc/n3/`, raw outputs in `loc/n3/*.txt`. The scripts are:
- `n3_lib.py`: the exact second-chaos model and the leaf-mixture CGF;
- `a1_check_cgf.py`, `a1c.py`: the CGF against Monte Carlo;
- `a2_profiles.py`: class profiles;
- `a3_gain.py`: the gain model;
- `a4_weights.py`: weight tables;
- `a5_k2check.py`: the closed-form rung-2 defect;
- `b1_bilap.py`, `b2_eval.py`: two rungs on the Gaussian-closure chain;
- `b3_k3bilap.py`: two rungs on the kappa_3 chain;
- `b4_eval_k3.py`: evaluation of the kappa_3-chain runs;
- `c1_closedwalk.py`: closed walks and their noncommutative surrogates;
- `c2_kappa_truth.py`, `c2b_kappa_truth.py`: closed walks against Monte Carlo cumulants;
- `c3_blocks.py`: affordable restrictions of the surrogates;
- `a6_checks.py`: the operator-valued localization expansion and Proposition 1.2;
- `d1_visibility.py`: multi-rung visibility.

Conventions are those of F1:
- units of 2n^3 FLOPs, B = 1024 units;
- s is the posterior variance and u = 1 - s the localization depth;
- e(y) = E(y, I);
- delta = e - y.grad e - Lap e;
- b(s) is the error left after localizing to variance s, with b(1) = err;
- a_k = L^T H^k L are the open walks (vector-state moments), t_k = tr H^k the closed walks (trace-state moments).

## 0. Summary

**The question.** Rung 0 of the Euler-Stein ladder (E Lap F = E F), i.e. first-order localization, is blind to the
closed walks tr H^k and to the gain's kappa_4 half. Does the next rung see them, at what weights, at what cost, and
where does genuine noncommutativity enter?

1. **Theorem A (the leaf mixture is exactly solvable; rung-walk visibility; proved and checked to 4 digits).**
   - *The leaf mixture.* For a second-chaos pre-activation read by the Gaussian closure, the localized estimator at
     depth u is the readout of an explicit Gaussian leaf mixture. Its CGF is a resolvent/log-det deformation of the
     truth's:

         K_u(theta) = theta mu0 + sum_r u^r Omega_r(theta),
         Omega_r = (theta^2/2) <L, p_r(theta H) L> + tr q_r(theta H),

     with p_r and q_r explicit polynomials of degree 2r.
   - *Rung r* (the r-th derivative of the localization profile at s = 1) sees exactly the walks with at most 2r
     letters H. In NCG terms these are the moments of degree <= 2r of the second-chaos operator in both states, the
     vector state omega_L = <L, . L> and the trace tr.
   - *First-visibility rule.* A walk with k letters is first seen at rung ceil(k/2). Rung 1 sees a_1 and a_2 only;
     rung 2 adds t_3, t_4, a_3, a_4 and the products of rung-1 classes.
   - *Profiles.* Every class has an explicit polynomial profile phi_c(s) = 1 - rho_c(1 - s), and its rung weights are
     w_r(c) = (-1)^{r+1} r! [u^r] rho_c. Examples:

     | class | rho_c(u) | w_1 | w_2 |
     |---|---|---|---|
     | kappa_3 path a_1 | 2u - u^2 | 2 | 2 |
     | kappa_4 path a_2 | u + u^2 - u^3 | 1 | -2 |
     | triangle t_3 | 3u^2 - 2u^3 | 0 | -6 |
     | 4-cycle t_4 | 2u^2 - u^4 | 0 | -4 |

   - *The gain.* At leading order in 1/k, its kappa_3 half has the a_2 profile and its kappa_4 half the 4-cycle
     profile. The measured weights at k = 50 are (0.98, -2.01) and (0.00, -3.99). So rung 2 sees the gain's kappa_4
     half at weight -4.
   - *Operator-valued localization.* Under every Gaussian localization, including an operator-valued
     (anisotropic, non-commuting) leaf covariance U, the first-order term contains no trace-state moment. **First-order
     localization is blind to closed walks for every localization scheme.** The second-order closed content is
     (theta^3/2) tr(U H U H^2) + (theta^4/4) tr(U H^2 U H^2).
2. **Theorem B (Gamma-calculus telescoping; where noncommutativity enters; proved, checked).**
   - *Rung r involves the iterated carre-du-champ forms Gamma_k for k <= r* (Bakry-Emery for the flat heat operator).
     Gamma_1 is the tangent Gram: open walks, vector state, hub class n^3. Gamma_2 is the Hilbert-Schmidt pairing of
     Hessians: the trace state of the second-chaos tuple.
   - *The rung-2 chain rule* is written out for an arbitrary layer map. For a Gaussian readout the closed walks enter
     **only** through (1/2) D^2 Phi[Gamma_2(X, X)] = 2 g_mv tr H^3 + 2 g_vv tr H^4, and all other terms are open walks.
   - *No free lunch.* The closed-walk part of rung 2 is exactly the chain's own closed walks. A two-rung defect is
     never a cheaper route to closed walks than computing them; it only certifies their weight.
3. **Theorem C (combining rungs; proved).**
   - Every second-chaos class profile has a double zero at s = 0.
   - The R-rung Hermite-Richardson estimator err_R = sum_r beta_r b^(r)(1), with beta_r = (-1)^{r+1}(R+1-r)/((R+1) r!),
     is exact on every class whose profile has degree <= R + 1. Its homogeneous forms are:
     - R = 1: T = (3e + Lap e)/4;
     - R = 2: T = (15e + 10 Lap e + Lap^2 e)/24, exact on a_1, a_2, a_3, t_3 and capturing 2/3 of t_4 and of the gain's
       kappa_4 half;
     - R = 3: T = (105e + 105 Lap e + 21 Lap^2 e + Lap^3 e)/192.
4. **Measured on small networks (Monte Carlo truth).**
   - *The Gaussian-closure chain, 8 networks, n = 16-64.*
     - Two rungs, held out (leave one network out, global coefficients), leave a residual MSE ratio with geometric
       mean 0.042, against 0.134 for one rung: **3.2x better, on 8/8 networks**.
     - The parameter-free Hermite T_2 reaches 0.085, against 0.135 for F5's parameter-free midpoint.
     - The demeaned explained share rises from 0.54-0.99 to 0.79-1.00.
   - *Caution.* That chain omits both open classes, so part of the gain is separating a_1 from a_2.
   - *The kappa_3 chain* (production-like, 4 networks, n = 32). Two rungs held out reach 0.110 against 0.247 for one
     rung (2.2x better, on 4/4). The demeaned explained share rises from 0.67 to 0.87. The parameter-free Hermite T_2
     reaches 0.144, beating the fitted one-rung merge. The pooled weights (0.577, -0.112) are close to Hermite's
     (0.667, -0.167).
5. **Closed walks as noncommutative moments: Toeplitz/Berezin deterministic equivalents (section 5; measured).**
   - *The structure.* The chain's second chaos is a Toeplitz operator on the birth set: H_i = sum_m d_im |l_m><l_m|,
     with symbol d_i and coherent states the birth legs. This is exactly the form T_f = sum_x f(x)|k_x><k_x| of stage-8
     Theorem D.
   - *The surrogate.* Closed walks are its traces. Keeping the non-crossing (operator-valued semicircular, amalgamated
     over the birth diagonal) overlap diagrams gives O(n M^2) surrogates. Its R^2 against the exact per-neuron tr H^3
     is 0.98-0.997, and against tr H^4 0.93-0.99 (n = 64-256, layers 3-11).
   - *The slope.* It is depth dependent and tends to 1 with width (layer 6: 1.29, 1.18, 1.10 at n = 64, 128, 256), so
     it is calibrated truth-free against exact closed walks offline.
   - *The classical-symbol term alone* (the upper-symbol moment sum_m sigma_im^3) costs about nothing and has
     R^2 0.82-0.93.
   - *Old births carry the closed walks of deep neurons.* A young-tier-only surrogate fails.
6. **Component and prediction.**
   - *The defect merge does not pay.* As a runtime component, the two-rung defect costs +56-118 units in its late
     form and is break-even at best at the production class mix.
   - *What carries the value: C-T.* Rung 2's value is its closed-walk content, and that content is better added
     directly. The overlap closed-walk readout C-T adds the Berezin pair (overlap) diagrams of the triangle and the
     4-cycle to the kappa_3 and kappa_4 diagonals at weight 1, as V56 added the path class.
   - *What C-T leaves out.* The self-loops (the classical symbol) are the births' own skewness, which the chain
     already carries, so they are not added.
   - *Cost.* About 1 unit per young source-layer, i.e. +16 units at layers 12-15.
   - *Prediction.* For the production chain a two-rung defect explains X*_12 = 0.54-0.63 of the per-neuron error,
     against X*_1 = 0.23-0.42 for one rung, about 2x. Rung 2 restricted to open walks adds only 0.05.
   - *C-T prediction.* Raw -5% [-1%, -13%] and adjusted -3% [+6%, -10%]. This is a modest component.
   - *The decisive experiment* is offline and truth-backed (section 7). It measures the overlap closed walks'
     share f_cw of the production residual (prior 0.12 [0.03, 0.30]) and the fidelity of the surrogate. The same
     measurement fixes the two-rung increment, X*_12 - X*_1 = f_cw^total + 0.03-0.06, which at n = 1024 cannot be
     measured directly (it would take 2n^2 chain runs).

## 1. Rungs, operator-valued constraints and the defect hierarchy

### 1.1 The ladder as an operator-valued moment constraint

**Proposition 1.1 (proved).** Let F: R^n -> R be 1-homogeneous with polynomial growth, X ~ N(0, I). The moment
functional phi_F(p) := E[p(X) F(X)] on polynomials satisfies

    phi_F( (|x|^2 - n - 1) p - x.grad p ) = 0      for every polynomial p.                                    (1.1)

On homogeneous p of degree k, (1.1) is the k-th tensor rung tr_{(k+1,k+2)} c_{k+2} = (1 - k) c_k (F1 Theorem 3), with
c_k = E D^k F.

*Proof.* Gaussian integration by parts gives E[g d_j F] = E[(x_j g - d_j g) F]. With g = x_j p, summed over j:
E[p x.grad F] = E[(|x|^2 p - n p - x.grad p) F]. Euler's identity x.grad F = F gives (1.1). For p = x^alpha with
|alpha| = k, expand in Hermite polynomials: |x|^2 x^alpha - (n + 1 + k) x^alpha pairs with the chaos coefficients as
tr c_{k+2} + (k - 1) c_k. QED.

- **Multi-output and products.** The same holds for each output, and for products G = F_a F_b ... with n + 1
  replaced by n + deg G. These product identities are implied by the single-output ones (homogeneity of F gives
  homogeneity of every product), so they carry no new truth-free information. In particular they do not determine
  the triangle tr(H_a H_b H_c) from carried data. They relate it to identity-traced higher chaos, which the single
  rungs already fix.
- **On a chain's chaos representation** (mu_i, L_i, H_i, the third-chaos star, ...) the rungs are consistency
  conditions:
  - rung 0: mu_i = tr H_i (note XLI: correlation 0.73-0.95 for the bare sources);
  - rung 1: the third-chaos kernel is traceless;
  - rung 2: tr_{34} c_4 = -H_i, so the fourth chaos of every neuron has a fixed partial trace. A pure second-chaos
    truncation violates rung 2 maximally.

  These are representation defects, not estimator defects. What follows is about the latter.

### 1.2 The defect hierarchy

Let K := d_s - (1/2) Lap_m act on functions of (m, s), with the input law N(m, sI). Put
g(s) := E_{m ~ N(0,(1-s)I)} E(m, sI) and b(s) := truth - g(s). F1 Theorem 1 gives g'(s) = E_m[K E(m, s)]; iterating
(the m-average is the heat semigroup e^{(1-s)Lap/2}):

    b^(r)(1) = - Delta_r,     Delta_r := (K^r E)(0, 1).                                                       (1.2)

These are the **rung-r defects**: truth-free, per neuron, and zero for the exact functional.

**Proposition 1.2 (homogeneous form; proved).** For a homogeneous estimator, with Lambda_K := Lap^K e(0),

    Delta_r = sum_{K=0}^{r} Lambda_K nabla^K p_r(0) / (2^K K!),     p_r(i) := prod_{j=0}^{r-1} (1/2 - j - i),

where nabla is the forward difference in i. Explicitly:

    Delta_1 = (e - Lap e)/2 = delta(0)/2,
    Delta_2 = (Lap^2 e + 2 Lap e - e)/4 = -(delta + Lap delta)(0)/4,
    Delta_3 = -(Lap^3 e + 9 Lap^2 e + 9 Lap e - 3e)/8.

Each vanishes on the Euler-Stein sequence Lap^K v(0) = (1, 1, -1, 3, -15, ...) v(0).

*Proof.*
- On E(m, sI) = s^{1/2} e(m/s^{1/2}), K maps s^a f(m/sqrt s) to s^{a-1}[(a - A) f](m/sqrt s), with A = (1/2)(Lap + y.grad)
  the repulsive OU generator. Hence Delta_r = [prod_{j<r} (1/2 - j - A) e](0).
- The conjugation 2A = e^{Lap/2} (y.grad) e^{-Lap/2} holds because [Lap, y.grad] = 2 Lap. So
  p(A) = e^{Lap/2} p(y.grad/2) e^{-Lap/2}.
- Evaluate at 0 on the Taylor components: e^{-Lap/2} maps the degree-2i part of Lap-powers through (-1/2)^i/i!, and
  e^{Lap/2} at 0 keeps (1/2)^j Lap^j/j!. Collecting the coefficient of Lambda_K gives
  sum_{j+i=K} p(j) 2^{-K} (-1)^i/(j! i!) = nabla^K p(0)/(2^K K!).
- The second form of Delta_2 uses Lap delta(0) = -Lambda_1 - Lambda_2. QED.

- **Tensor defects.** D^k delta(0) = (1 - k) c_k^e - tr c_{k+2}^e is the violation of the k-th tensor rung by the
  chain's own chaos coefficients c_k^e = D^k e(0).
- **Which tensor rungs the mean error sees.**
  - The isotropic Green integral (F1 Corollary 2) sees only the full traces Lap^j delta(0) of the even ones.
    Delta_2 is the first that contains the traced rung-2 tensor defect.
  - The odd and partially traced rungs enter only through anisotropic localization (section 2, (v)) and through the
    telescoping (section 3), where grad(K X), the gradient of the rung-1 defect of an inner table, is needed.

## 2. Theorem A: the exact leaf mixture and rung-walk visibility

**Setting.**
- The model is z = mu0 + L.x + (1/2)(x^T H x - tr H), x ~ N(0, I_n). This is note XLI's second-chaos
  (generalized chi-square) pre-activation, with cumulants kappa_k = k![t_k/(2k) + a_{k-2}/2].
- E(m, s) = g(mu(m, s), v(m, s)) is the Gaussian-closure readout with g(mu, v) = E relu(N(mu, v)). It is fed the exact
  posterior mean and variance:
  - mu(m, s) = mu0 + L.m + (1/2)(m^T H m - u tr H);
  - v(m, s) = s |L + H m|^2 + (s^2/2) tr H^2.
- The leaf mixture Y_u is the mixture over m ~ N(0, uI) of N(mu(m, s), v(m, s)). Then g(s) = E relu(Y_u): b(s) is the
  readout of a law, not just a number.

**Theorem A.**
1. **(Closed form; proved.)** With x := theta h_j over the eigenvalues h_j of H and l_j := (Q^T L)_j:

       K_u(theta) = theta mu0 + sum_j [ -u x/2 + s^2 x^2/4 - (1/2) log(1 - u x (1 + s x)) ]
                    + sum_j l_j^2 [ theta^2 s/2 + (u theta^2/2)(1 + s x)^2 / (1 - u x (1 + s x)) ].

   Equivalently, in operator form:
   - closed part: -(1/2) log det(I - u theta H (I + s theta H)) - (u theta/2) tr H + (s^2 theta^2/4) tr H^2;
   - open part: (theta^2/2) <L, [s + u (I + s theta H)(I - u theta H (I + s theta H))^{-1} (I + s theta H)] L>.

   At u = 1 this is the truth's CGF; at u = 0 it is the Gaussian closure's.
2. **(Rung expansion; proved.)**

       K_u = theta mu0 + theta^2 V/2 + sum_{r>=1} u^r Omega_r,
       Omega_r(theta) = (theta^2/2) <L, p_r(theta H) L> + tr q_r(theta H),

   with:
   - p_1(x) = 2x + x^2;
   - p_2(x) = -x + x^2 + 3x^3 + x^4;
   - q_1 = 0;
   - q_2(x) = x^3/2 + x^4/4;
   - q_3(x) = -x^3/3 + x^5/2 + x^6/6.

   In general deg p_r = deg q_r = 2r.
3. **(Visibility; proved.)** b^(r)(1) is a polynomial in {a_k, t_k : k <= 2r} times derivatives of g at the base
   law. A walk with k letters H (open: a_k, cumulant order k + 2; closed: t_k, cumulant order k) first appears at rung
   ceil(k/2). A product class prod_i c_i first appears at rung sum_i ceil(k_i/2).
4. **(Profiles and weights; proved.)** For each walk class c (a monomial in the cumulants times the matching Edgeworth
   derivative), b(s) contains err_c phi_c(s) with phi_c(s) = 1 - rho_c(1 - s), where:
   - closed walk t_k: rho_k^c(u) = k sum_j (u^j/j) C(j, k - j) s^{k-j};
   - open walk a_k: rho_k^o(u) = sum_j u^{j+1} C(j + 2, k - j) s^{k-j};
   - product class: rho = product of the factors' rho.

   The rung-r weight is w_r(c) = phi_c^(r)(1) = (-1)^{r+1} r! [u^r] rho_c.
5. **(Operator-valued localization; proved.)** Let the leaves be N(m, I - U), m ~ N(0, U), with 0 <= U <= I an
   arbitrary matrix (any Chen-Eldan control, accumulated). The O(U) term of the CGF is
   theta^3 <L, H U L> + (theta^4/2) <L, H U H L>: open walks with one U insertion, **and no trace-state term**. The
   closed part of the O(U^2) term is

       (theta^3/2) tr(U H U H^2) + (theta^4/4) tr(U H^2 U H^2).

6. **(Several outputs; proved.)** For a tuple (z_a) with (L_a, H_a), the joint leaf-mixture CGF is item 1 with
   theta H -> Theta := sum_a theta_a H_a and theta L -> sum_a theta_a L_a.
   - Rung weights therefore depend only on word length.
   - Noncommutative words tr(H_{a1} ... H_{ak}) enter through tr Theta^k, so tr(H_a H_b H_a H_b) and
     tr(H_a^2 H_b^2) appear at rung 2 with their symmetrized-word coefficients.

*Proof.*
1. Given m, the leaf law of z is mu(m, s) + sqrt(s) w.Z + (s/2)(Z^T H Z - tr H), with w = L + H m. That gives mu and
   v. The CGF of the mixture is E_m exp(theta mu + theta^2 v/2). This is the Gaussian integral of exp(b.m + m^T M m/2)
   with:
   - b = theta (I + s theta H) L;
   - M = theta H + s theta^2 H^2;
   - a constant theta mu0 - theta u tr H/2 + theta^2 s|L|^2/2 + theta^2 s^2 tr H^2/4.

   For m ~ N(0, uI) this equals det(I - uM)^{-1/2} exp((u/2) b^T (I - uM)^{-1} b) times the constant's exponential.
   Diagonalize.
2. Expand in u at fixed theta, using s = 1 - u.
   - *Closed part.* -(1/2) log(1 - u(x + x^2) + u^2 x^2) - ux/2 + (1 - u)^2 x^2/4. The u^1 coefficient is
     (x + x^2)/2 - x/2 - x^2/2 = 0. The u^2 coefficient is (x + x^2)^2/4 - x^2/2 + x^2/4 = x^3/2 + x^4/4. The u^3
     coefficient follows the same way.
   - *Open part.* s + u(1 + sx)^2 sum_j u^j x^j (1 + sx)^j. Each power of u contributes at most two powers of x per
     m-contraction.
3. b(s) = truth - E relu(Y_u) = [exp(sum_k kappa_k(z) D^k/k!) - exp(sum_k kappa_k(Y_u) D^k/k!)] g at (mu0, V), with
   D = d/dmu (mean and variance are equal at every u). The r-th u-derivative at 0 involves only Omega_1, ..., Omega_r.
4. Read the coefficient of each walk from item 1. The closed walk t_k is the coefficient of x^k in
   (1/2) sum_j u^j x^j (1 + sx)^j / j. The open walk a_k is the coefficient of x^k in sum_j u^{j+1} x^j (1 + sx)^{j+2}.
   Products follow because exp is multiplicative in the cumulants.
5. With covariance U the same Gaussian integral has:
   - b = theta (I + theta H S) L and M = theta H + theta^2 H S H, with S = I - U;
   - leaf variance (L + Hm)^T S (L + Hm) + (1/2) tr(SHSH);
   - leaf mean shift -tr(UH)/2.

   Expand -(1/2) log det(I - U^{1/2} M U^{1/2}) = (1/2) tr(UM) + (1/4) tr((UM)^2) + ....
   - *First order.* (1/2) tr(U(theta H + theta^2 H^2)) - (theta/2) tr(UH) - (theta^2/2) tr(U H^2) = 0.
   - *Second order, theta^2 tr(UHUH).* The coefficient is 1/4 - 1/2 + 1/4 = 0. The three contributions come from
     tr(SHSH)/4, from the first-order S inside M, and from tr((UM)^2)/4.
   - *Second order, remainder.* (theta^3/2) tr(UHUH^2) + (theta^4/4) tr(U H^2 U H^2).
6. The joint closure's CGF of a leaf is theta.mu + theta^T V theta/2, with
   V_ab = (L_a + H_a m)^T S (L_b + H_b m) + (1/2) tr(S H_a S H_b). Its quadratic form in m is Theta + Theta S Theta,
   so the integral is item 1 with theta H -> Theta. QED.

### 2.1 Class weights per rung (`a4_weights.py`)

"Hermite R" is the capture of the R-rung Hermite-Richardson estimator of Theorem C (1 = exact).

| class | rho_c(u) coefficients (u^0, u^1, ...) | w_1 | w_2 | w_3 | Hermite R = 1 | R = 2 | R = 3 |
|---|---|---|---|---|---|---|---|
| a_1 (kappa_3 path, D3 open) | 0, 2, -1 | +2 | +2 | 0 | 1.00 | 1.00 | 1.00 |
| a_2 (kappa_4 path, P4) | 0, 1, 1, -1 | +1 | -2 | -6 | 0.50 | 1.00 | 1.00 |
| a_3 (kappa_5 path) | 0, 0, 3, -2 | 0 | -6 | -12 | 0 | 1.00 | 1.00 |
| a_4 (kappa_6 path) | 0, 0, 1, 3, -4, 1 | 0 | -2 | +18 | 0 | 0.33 | 1.25 |
| t_3 (triangle) | 0, 0, 3, -2 | 0 | **-6** | -12 | 0 | 1.00 | 1.00 |
| t_4 (4-cycle) | 0, 0, 2, 0, -1 | 0 | **-4** | 0 | 0 | 0.67 | 1.00 |
| t_5 | 0, 0, 0, 5, -5, 1 | 0 | 0 | +30 | 0 | 0 | 1.25 |
| t_6 | 0, 0, 0, 2, 3, -6, 2 | 0 | 0 | +12 | 0 | 0 | 0.50 |
| a_1^2 (kappa_3^2) | 0, 0, 4, -4, 1 | 0 | -8 | -24 | 0 | 1.33 | 1.00 |
| a_1 a_2 (kappa_3 kappa_4) | 0, 0, 2, 1, -3, 1 | 0 | -4 | +6 | 0 | 0.67 | 1.25 |
| a_2^2 (kappa_4^2) | 0, 0, 1, 2, -1, -2, 1 | 0 | -2 | +12 | 0 | 0.33 | 1.00 |
| gain, kappa_3 half (= a_2, measured) | | +0.98 | -2.01 | | 0.5 | 1.00 | |
| gain, kappa_4 half (= t_4, measured) | | 0.00 | -3.99 | | 0 | 0.67 | |

- **Rung 1 reproduces the corrected F1/F5 class accounting:**
  - kappa_4 path at weight 1;
  - kappa_3 path at 2;
  - gain kappa_3 half at 1;
  - closed walks and the gain kappa_4 half at 0.

  It is exact at all orders, not only at leading order. The rung-1 defect of one Gaussian readout is exactly
  b'(1) = a_1 g''' + (a_2/2) g'''' (`a5_k2check.py`: 6 digits at finite H).
- **Rung 2 sees the closed walks**, at large negative weights: -6 for the triangle and -4 for the 4-cycle and for the
  gain's kappa_4 half.
- **Closed form at finite H** (`a5_k2check.py`; checked against the exact profile to 1e-4):

      b''(1) = g'''(a_1 - t_3) - g''''(a_2 + t_4/2) - 3 a_3 g^(5) - (a_1^2 + a_4) g^(6) - a_1 a_2 g^(7) - (a_2^2/4) g^(8).

### 2.2 Verification (`a1c.py`, `a2_profiles.py`, `a5_k2check.py`)

**The CGF.** Leaf-mixture readout by characteristic-function inversion of item 1, against Monte Carlo over leaves
(2e7 draws; n = 5, random L and H):
- u = 0.3: CF 0.511520, MC 0.511577 +- 5.8e-5;
- u = 0.7: CF 0.494969, MC 0.494964 +- 9.7e-5;
- truth: CF 0.488713, MC 0.488794 +- 1.2e-4.

**The profiles.** Each class is isolated by structure and by parity in a scale epsilon of H, with Richardson in
epsilon. "Measured" is b(s)/b(1) from the exact CF profile; "predicted" is phi_c(s).

| s | 0.10 | 0.25 | 0.40 | 0.55 | 0.70 | 0.85 | 1 | err: measured / Edgeworth |
|---|---|---|---|---|---|---|---|---|
| a_1 measured | 0.0100 | 0.0625 | 0.1600 | 0.3025 | 0.4900 | 0.7225 | 1 | -1.53710e-2 / -1.53710e-2 |
| a_1 predicted | 0.0100 | 0.0625 | 0.1600 | 0.3025 | 0.4900 | 0.7225 | 1 | |
| a_2 measured | 0.0190 | 0.1094 | 0.2560 | 0.4386 | 0.6370 | 0.8309 | 1 | -1.06260e-1 / -1.06260e-1 |
| a_2 predicted | 0.0190 | 0.1094 | 0.2560 | 0.4386 | 0.6370 | 0.8309 | 1 | |
| t_3 measured | 0.0280 | 0.1562 | 0.3520 | 0.5747 | 0.7840 | 0.9392 | 1 | -1.16877e-2 / -1.16895e-2 |
| t_3 predicted | 0.0280 | 0.1562 | 0.3520 | 0.5748 | 0.7840 | 0.9392 | 1 | |
| t_4 measured | 0.0361 | 0.1912 | 0.4093 | 0.6357 | 0.8279 | 0.9554 | 1 | -1.0736e-1 / -1.0760e-1 |
| t_4 predicted | 0.0361 | 0.1914 | 0.4096 | 0.6360 | 0.8281 | 0.9555 | 1 | |
| a_3 measured | 0.0280 | 0.1563 | 0.3522 | 0.5749 | 0.7841 | 0.9393 | 1 | -3.8691e-1 / -3.8668e-1 |
| a_3 predicted | 0.0280 | 0.1562 | 0.3520 | 0.5748 | 0.7840 | 0.9392 | 1 | |

The residual differences (t_4, a_3: 3e-4) are the next order in epsilon.

**Operator-valued localization and Proposition 1.2** (`a6_checks.py`).
- *Theorem A(5).* Generic positive U, not commuting with H, at n = 5, with the exact Gaussian-integral CGF and a
  polynomial fit in the scale of U, at theta = 0.3 and 0.7:
  - the O(U) closed coefficient is 1e-13 (prediction 0);
  - the O(U) open coefficient is +4.985714e-3 and +6.799000e-2, equal to theta^3 L^T H U L + (theta^4/2) L^T H U H L to
    all printed digits;
  - the O(U^2) closed coefficient is -2.308541e-3 and -2.391606e-2, equal to
    (theta^3/2) tr(UHUH^2) + (theta^4/4) tr(UH^2UH^2) to all printed digits.
- *Proposition 1.2* (random polynomial, 1-D). prod_{j<r}(1/2 - j - A) e(0), the difference formula and the explicit
  forms agree to 10 digits for r = 1, 2, 3.

### 2.3 The gain (`a3_gain.py`)

**The model** (F5 referee R1).
- x = (y in R^k, g), z = A(y)(mu + sigma g), A = |y|/sqrt k, k = 50.
- The Gaussian closure is fed the exact noncentral-chi moments.
- The profile b(s) is computed by 80 x 80 Gauss-Laguerre/Hermite quadrature over the leaves (|m_y|^2, m_g).

| alpha = mu/sigma | err | b(s)/b(1) at s = 0.1, 0.25, 0.4, 0.55, 0.7, 0.85 | w_1 | w_2 |
|---|---|---|---|---|
| 0 (kappa_4 half only) | -1.990e-3 | 0.0367 0.1940 0.4132 0.6390 0.8296 0.9558 | -0.000 | -3.990 |
| 1 (kappa_3 half only) | -2.420e-3 | 0.0193 0.1110 0.2594 0.4432 0.6414 0.8336 | +0.980 | -2.010 |
| 0.5 (mixed) | -2.195e-3 | 0.0298 0.1609 0.3518 0.5607 0.7541 0.9067 | +0.396 | -3.171 |

- The kappa_4 half follows the 4-cycle profile (1 - u^2)^2 (0.0361, 0.1914, ...).
- The kappa_3 half follows the kappa_4-path profile s^2(2 - s) (0.0190, 0.1094, ...).
- The differences are O(1/k).

**Proposition 2.1 (sketch, leading order in 1/k).**
- Within a leaf N(m_y, sI), Var(A) = (2s - s^2)/(2k) + O(k^{-2}), a fraction 1 - u^2 of the total. This is the
  radial (isotropic second-chaos) mode split by localization.
- By the law of total cumulance, the kappa_4 class of b is E[kappa_4 | leaf] + 4 Cov(mean_leaf, kappa_3 | leaf). That
  is 12 sigma^4 nu [s^2 + 2s(1 - s)](1 - u^2) = 12 sigma^4 nu (1 - u^2)^2.
- The kappa_3 class is E[kappa_3 | leaf] = 6 mu sigma^2 nu s (1 - u^2) = 6 mu sigma^2 nu s^2 (2 - s).

**Reading.** The gain is a second-chaos object in disguise: an open walk of length 2 (kappa_3 half) and a 4-cycle
(kappa_4 half) of the radial operator. Rung 2 sees its kappa_4 half at weight -4.

### 2.4 What is noncommutative here, and what is not

- **One neuron.** For one neuron the relevant algebra is generated by the single self-adjoint H. Its two states are:
  - the vector state omega_L, with spectral measure sum_j l_j^2 delta_{h_j}, giving the open walks;
  - the trace, with sum_j delta_{h_j}, giving the closed walks.

  Theorem A says the localization filtration of this moment problem is by degree: rung r reveals the moments of
  degree <= 2r of both measures. That is commutative.
- **Genuine noncommutativity** enters in three places:
  1. **Several outputs (item 6).** The covariance and cross-cumulant tables see words in non-commuting H_a. For
     example kappa(z_a, z_a, z_b, z_b) contains both tr(H_a^2 H_b^2) and tr(H_a H_b H_a H_b). Rung 2 sees them through
     tr Theta^4, with length-only weights, but the chain must compute each word.
  2. **Anisotropic localization (item 5).** A leaf covariance U that does not commute with H. The collective
     localization of F3 is the case U = P_k (top-k collective modes). Its rung-2 closed content tr(P H P H^2) consists
     of closed walks that pass through the collective subspace twice.
  3. **The birth-space Toeplitz structure (section 5).** In the chain H_i = Lambda D_i Lambda^T. The closed walk is
     tr((D_i G)^k) for the non-commuting pair formed by the diagonal symbol algebra (D_i) and the Gram of the birth
     legs (G = Lambda^T Lambda).

     This is where an NC tool buys something. Operator-valued free probability with amalgamation over the birth
     diagonal gives the closed walks at n M^2 instead of n^2 M per output (section 5).
- **The Berezin-Toeplitz reading of the leaf mixture.** The localized estimator g(s) averages the Gaussian closure
  over coherent states N(m, sI) (anti-Wick/Toeplitz quantization of the closure). The truth is the Wick (exact) side.
  Theorem A's u-expansion is the semiclassical expansion of the gap between the two, with the localization depth u as
  the expansion parameter. Each order of u is one more Wick contraction of the leaf means, and each contraction closes
  at most two more letters of a walk.

## 3. Theorem B: the Gamma-calculus, the rung-2 telescoping, and where the trace state enters

### 3.1 Iterated carre du champ

For smooth f, g of (m, s) put Gamma_0(f, g) = fg and, for k >= 1,

    Gamma_k(f, g) := sum_{j1..jk} (d_{j1} ... d_{jk} f)(d_{j1} ... d_{jk} g) = <grad^k f, grad^k g>_HS.

**Lemma 3.1 (proved).**

    K(fg) = f Kg + g Kf - Gamma_1(f, g),
    K Gamma_k(f, g) = Gamma_k(Kf, g) + Gamma_k(f, Kg) - Gamma_{k+1}(f, g).

*Proof.* Lap(fg) = f Lap g + g Lap f + 2 grad f . grad g. Apply this to each product d_J f d_J g, and use that K commutes
with d_J (constant coefficients). QED.

- **Bakry-Emery.** With L = Lap/2, the Bakry-Emery carre du champ is Gamma(f) = |grad f|^2/2 and its iterate is
  Gamma_2(f) = |grad^2 f|^2_HS/4 (flat Bochner). So Gamma_1 and Gamma_2 here are the Bakry-Emery forms up to the
  factors 2 and 4.
- **The noncommutative version.** The Cipriani-Sauvageot carre du champ on an operator algebra has the same algebra.
  Here it acts on the chaos algebra, where Gamma_2 becomes a trace pairing (Corollary 3.4).

**Proposition 3.2 (grading of rung r; proved by induction on r with Lemma 3.1).** K^r(Phi o X) is a finite sum of
multilinear terms D^j Phi[A_1, ..., A_j]. Each A_i is K^a X or a nested form Gamma_k(K^a X, ...). Each term has total
weight r, where a K counts 1 and a Gamma_k counts k.
- The only term with a Gamma_r is (-1)^r (1/2) D^2 Phi[Gamma_r(X, X)].
- Hence rung r is the first rung that contains the k-th derivative pairing Gamma_k for k = r.

**Proposition 3.3 (rung-2 chain rule; proved).** For any smooth layer map Phi and table X:

    K^2(Phi o X) = DPhi[K^2 X] + D^2 Phi[KX, KX] - 2 D^2 Phi[Gamma_1(X, KX)] - D^3 Phi[KX, Gamma_1(X, X)]
                   + (1/2) D^2 Phi[Gamma_2(X, X)] + (1/2) D^3 Phi[Gamma_1(X, Gamma_1(X, X))]
                   + (1/4) D^4 Phi[Gamma_1(X, X), Gamma_1(X, X)].                                            (3.1)

The contractions are:
- D^2 Phi[Gamma_1(X, KX)] = sum_ab Phi_ab Gamma_1(X_a, K X_b);
- D^3 Phi[KX, Gamma_1] = sum_abc Phi_abc K X_a Gamma_1(X_b, X_c);
- D^3 Phi[Gamma_1(X, Gamma_1)] = sum_abc Phi_abc Gamma_1(X_a, Gamma_1(X_b, X_c)).

*Proof.* Start from K(Phi o X) = DPhi[KX] - (1/2) D^2 Phi[Gamma_1(X, X)]. Apply K to both terms with Lemma 3.1. The
Gamma_2 term comes from K Gamma_1(X, X) = 2 Gamma_1(X, KX) - Gamma_2(X, X). QED.

**Check.** Take Phi = g(mu, v), a Gaussian readout fed exact moments of the second-chaos model, at the base point.
- *Inputs.* KX = (0, |L|^2), K^2 X = (0, -tr H^2), Gamma_1(X, X) = [[a_0, 2a_1], [2a_1, 4a_2]], and
  Gamma_2(X, X) = [[t_2, 2t_3], [2t_3, 4t_4]].
- *Result.* With g_v = g''/2, (3.1) reproduces term by term the closed form of section 2.1:
  g'''(t_3 - a_1) + g''''(a_2 + t_4/2) + 3a_3 g^(5) + (a_1^2 + a_4) g^(6) + a_1 a_2 g^(7) + (a_2^2/4) g^(8). That form
  is checked against the exact profile to 1e-4.
- *Where the closed walks come from.* They come **only** from (1/2) D^2 Phi[Gamma_2] = (1/2) g_mm t_2 + 2 g_mv t_3 + 2 g_vv t_4.
  The t_2 piece cancels against DPhi[K^2 X] = -g_v t_2 because the variance is exact.

**Corollary 3.4 (Gamma_k in the chaos representation; proved).** At the base point:
- if grad mu_a = L_a and grad^2 mu_a = H_a;
- if v_a and C_ab are the second-chaos variance and covariance, so that grad^2 v_a = 2 H_a^2 and
  grad^2 C_ab = H_a H_b + H_b H_a;

then:

    Gamma_1(mu_a, mu_b) = L_a.L_b,  Gamma_1(mu_a, v_b) = 2 L_a^T H_b L_b,  Gamma_1(v_a, v_b) = 4 L_a^T H_a H_b L_b   (vector states),
    Gamma_2(mu_a, mu_b) = tr(H_a H_b),  Gamma_2(mu_a, v_b) = 2 tr(H_a H_b^2),  Gamma_2(v_a, v_b) = 4 tr(H_a^2 H_b^2),
    Gamma_2(mu_c, C_ab) = 2 tr(H_c H_(a H_b)),  Gamma_2(C_ab, C_cd) = tr((H_a H_b + H_b H_a)(H_c H_d + H_d H_c))   (trace state).

- **Gamma_1 is the vector state.** Gamma_1 evaluates vector-state moments, open walks: the hub class, n^3 per
  source-layer.
- **Gamma_2 is the trace state.** Gamma_2 evaluates the trace state of the second-chaos tuple:
  - for the per-neuron readout, the closed walks tr H_a^3, tr H_a^4;
  - for the covariance map, the genuinely noncommutative words tr(H_c H_a H_b) (the joint-gate triangle among three
    neurons) and tr(H_a H_b H_a H_b) against tr(H_a^2 H_b^2).

**Proposition 3.5 (no free lunch; proved).** In the rung-2 defect of a chain built from Gaussian-type layer maps, the
whole trace-state content is sum_l (linear transport) . (1/2) D^2 Phi_l[Gamma_2(X_l, X_l)].
1. An approximation of Delta_2 whose trace-state part has relative error eps is an approximation of the chain's own
   closed walks, Gamma_2 of its tables, to relative error eps, and conversely.
2. Rung 2 therefore cannot see the closed walks more cheaply than the closed walks can be computed. What it adds is
   their *weight in the error*: -6 for the triangle and -4 for the 4-cycle. The R = 2 Hermite combination returns them
   at weights 1 and 2/3.
3. When a surrogate of the closed walks is available, adding it to the readout at its Edgeworth weight (1 for the
   triangle, 1 for the 4-cycle) uses the same information without the other rung-2 terms. Section 6 shows this is
   strictly cheaper.

### 3.2 The cost of each rung via the telescoping

The telescoping is F1 Proposition 8 at rung 1 and (3.1) iterated over layers at rung 2. By R3.1 of F1's referee it
runs as a forward (tangent-linear) pass through all channels the later layers read. The A-late in-place variant
(SYNTHESIS A) lets the chain's own transport carry it.

| rung | new objects | NC type | cost at n = 1024 (units) |
|---|---|---|---|
| 1 | Gamma_1(X, X): grad mu = L (first chaos), grad v (D21 with an input leg), <grad mu, grad C> (cross), <grad C, grad C> (pair, V2) | vector state (open walks); pair = 2-word trace, measured <= 3% (GC, F1) | hub class, n^3 per source-layer; F1 + referee: +70-140 full depth, A-late +20-50 |
| 2, open part | Gamma_1(X, KX): the gradient of the rung-1 defect (the vector rung-1 tensor); Gamma_1(X, Gamma_1(X, X)): open walks a_3, a_4 (L^T H^3 L, \|H^2 L\|^2); products D^2 Phi[KX, KX], D^3, D^4 (elementwise) | vector state | 2 more hub families per source-layer (Y_2 = L_j . H_i v_i) and one Schur solve per layer: about +8-12 per layer |
| 2, readout closed part | Gamma_2(mu_a, v_a) = 2 tr H_a^3, Gamma_2(v_a, v_a) = 4 tr H_a^4 | trace state, single word per neuron | exact: K_4 contraction, n^3 M per layer (about 8000 units with M = 16n; offline only); NC surrogate (section 5): 2-5 per layer |
| 2, covariance-map closed part | Gamma_2(mu_c, C_ab), Gamma_2(C_ab, C_ab): tr(H_a^2 H_b), tr(H_a H_b H_a H_b), tr(H_a^2 H_b^2) for all n^2 pairs | genuinely noncommutative words | n^2 pairs x a trace each: n^5 exact; surrogate n^2 M^2; **out of budget in any form** |
| r | Gamma_k for k <= r: traces of words of length <= 2r | | grows with the word length |

**Reading.**
- The heat-defect ladder reproduces note XL's cost classification exactly. Rung 1 lives in the vector state (open
  walks, n^3). Rung 2 is the first rung that needs the trace state. The trace state is the K_4 wall per neuron and
  n^5 for the covariance map.
- "Rung 2 sees the closed walks" is true. "Rung 2 is a way to compute the closed walks" is false (Proposition 3.5).
- Its covariance-map part is unaffordable even with surrogates. So a runtime rung-2 defect must drop the
  covariance-map closed words, and its trace-state content is then only the per-neuron closed walks.

## 4. Theorem C: combining rungs

### 4.1 The double zero

**Lemma 4.1 (proved).** In the second-chaos model every class profile satisfies phi_c(0) = phi_c'(0) = 0.

*Proof.* Let s be the leaf variance and decompose z = mu(m, s) + zeta, where the leaf fluctuation zeta has
conditional mean 0 and variance v(m, s) = O(s). Its conditional cumulants of order >= 3 are O(s^2): kappa_3(zeta | m) =
3s^2 w^T H w + s^3 t_3, and the rest similarly.

By the law of total cumulance, kappa_k(z) - kappa_k(Y_u) is a sum of joint cumulants involving at least one
conditional cumulant of zeta of order >= 3, hence O(s^2). Every Edgeworth class of b = truth - E relu(Y_u) is a
polynomial in these differences, so it is O(s^2). Equivalently rho_c(1) = 1 and rho_c'(1) = 0, which the explicit rho
of Theorem A satisfy. QED.

- **Caveat for composed ReLU chains.** Their leaves also contain kinks inside the leaf. F1's kink-count heuristic
  suggests b ~ s^{v/2} near s = 0 for a diagram with v kink vertices. The double zero is therefore a property of the
  smooth per-readout model, tested empirically on networks below.

### 4.2 Hermite-Richardson

**Theorem C (proved).** Fix R >= 1. Suppose b(s) = s^2 q(s) with deg q <= R - 1. Then

    err = b(1) = sum_{r=1}^{R} beta_r^(R) b^(r)(1),     beta_r^(R) = (-1)^{r+1} (R + 1 - r) / ((R + 1) r!).

The estimator T_R := e + sum_r beta_r^(R) b^(r)(1) = e - sum_r beta_r^(R) Delta_r is exact on every class whose
profile has degree <= R + 1 (all of them, by Lemma 4.1).
- *Homogeneous forms* (Proposition 1.2):

      T_1 = (3e + Lap e)/4,
      T_2 = (15e + 10 Lap e + Lap^2 e)/24,
      T_3 = (105e + 105 Lap e + 21 Lap^2 e + Lap^3 e)/192.

  All are consistent: they return v(0) on the Euler-Stein sequence.
- *Captures* (section 2.1 table):
  - R = 1 is exact on a_1 only (a_2 at 1/2);
  - R = 2 is exact on a_1, a_2, a_3 and t_3, and gives 2/3 of t_4 and of the gain's kappa_4 half;
  - R = 3 is exact on all classes with at most four letters.

*Proof.* Write q = sum_{j<R} c_j (s - 1)^j and s^2 = 1 + 2(s - 1) + (s - 1)^2. Then
b^(r)(1) = r! (c_r + 2c_{r-1} + c_{r-2}) for r = 1..R, with c_j := 0 outside [0, R - 1]. This is a triangular system
in (c_0, ..., c_{R-1}). Its solution for c_0 = q(1) = b(1) is linear in b^(r)(1), with the stated coefficients.

Check, with gamma_r := r! beta_r = (-1)^{r+1}(R + 1 - r)/(R + 1). Inserting b^(r)(1) = r!(c_r + 2c_{r-1} + c_{r-2})
into sum_r beta_r b^(r)(1):
- the coefficient of c_0 is 2 gamma_1 + gamma_2 = (2R - (R - 1))/(R + 1) = 1;
- the coefficient of c_j for 1 <= j <= R - 1 is gamma_j + 2 gamma_{j+1} + gamma_{j+2} (with gamma_{R+1} = 0), and the
  bracket (R+1-j) - 2(R-j) + (R-1-j) vanishes.

The homogeneous forms follow by inserting Delta_1, Delta_2, Delta_3. QED.

- **The limit R -> infinity.** beta_r -> (-1)^{r+1}/r!, the Taylor expansion of b from s = 1 to s = 0. The finite-R
  weights (R + 1 - r)/(R + 1) are its Cesaro damping, which the double zero makes exact for polynomial profiles.
- **Relation to F5's T_2.** F5's T_2 = e - Delta_1 + Delta_2/2 assumes b = alpha s + beta s^2. That violates the double
  zero and is exact on a_1 only. Captures: a_2 -> 2, t_3 -> 3, t_4 -> 2.

### 4.3 Multi-rung OLS and the visibility bound

**Corollary 4.2 (proved).** Model the per-neuron error as a sum of uncorrelated class components e_c with energies E_c
and rung weights w_r(c). The best linear combination of the first R rungs explains

    X*_R = r^T S^{-1} r / sum_c E_c,     r_k = sum_c w_k(c) E_c,  S_kj = sum_c w_k(c) w_j(c) E_c.

For R = 1 this is SYNTHESIS Corollary 7.

*Proof.* This is the OLS projection of sum_c e_c onto the span of the R defect combinations. QED.

**Applied to class mixes for the production residual** (`d1_visibility.py`). "Invisible" covers truncation,
compression, calibration and classes with five or more letters.

| scenario (energy shares) | visible share | X*_1 (one rung) | X*_12 (exact rung 2) | X*_12, rung 2 without closed walks | fixed Hermite T_2 |
|---|---|---|---|---|---|
| closed-heavy (a_1 .15, a_2 .10, t_3 .25, t_4 .10, a_3/a_4 .05) | 0.65 | 0.23 | 0.63 | 0.29 | 0.63 |
| central (a_1 .20, a_2 .10, t_3 .15, t_4 .08, a_3/a_4 .05, gain-4 .02) | 0.60 | 0.28 | 0.58 | 0.34 | 0.58 |
| open-heavy (a_1 .30, a_2 .15, t_3 .05, t_4 .03, a_3/a_4 .02) | 0.55 | 0.42 | 0.54 | 0.47 | 0.54 |
| XLII Corollary 7 model (a_1 .25, a_2 .10, t_3 .15, t_4 .10) | 0.60 | 0.33 | 0.59 | 0.35 | 0.59 |

**Reading.**
- Two exact rungs reach essentially the whole visible share, and the parameter-free Hermite T_2 does as well as OLS
  under this model.
- The increment over one rung is +0.12 to +0.40, about **2x** at the central mix.
- Rung 2 restricted to open walks adds only +0.03 to +0.06. **The value of the second rung is its closed-walk
  content.**

### 4.4 Measured: two rungs on small networks (`b1_bilap.py`, `b2_eval.py`, `b3_k3bilap.py`)

**Method.**
- Lap e and Lap^2 e are taken by finite differences in the input mean: 2n + 1 runs at h = 0.02 for Lap, and
  2n^2 + 4n + 1 runs at h = 0.06 for Lap^2 (the 4th-difference pair stencil, F5 `bilaplacian`). Runs are float64.
- Truth is F5's Monte Carlo (noise/MSE < 1e-3 everywhere).
- Then Delta_1 = (e - Lap e)/2 and Delta_2 = (Lap^2 e + 2 Lap e - e)/4.
- The estimators are:
  - T1 = e - Delta_1 (F5's midpoint, p = 1 in variance);
  - H1 = e - Delta_1/2 (Hermite R = 1);
  - H2 = e - (4 Delta_1 - Delta_2)/6 (Hermite R = 2, parameter-free);
  - T2 = F5's e - Delta_1 + Delta_2/2;
  - LOO1 and LOO2 = one- and two-rung OLS with global coefficients fitted on the other networks (held out);
  - oracle1 and oracle2 = per-network fits.
- X1 and X12 are the demeaned explained shares.

**The Gaussian-closure chain** (MSE ratios to GC):

| net | GC MSE | T1 | H1 | **H2** | T2 | LOO1 | **LOO2** | oracle1 | oracle2 | X1 | X12 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| n16 L8 s1 | 1.08e-2 | 0.097 | 0.287 | 0.112 | 0.032 | 0.115 | 0.042 | 0.092 | 0.015 | 0.86 | 0.98 |
| n16 L8 s2 | 3.88e-4 | 0.834 | 0.630 | 0.192 | 0.728 | 0.688 | 0.185 | 0.624 | 0.073 | 0.67 | 0.93 |
| n16 L8 s3 | 2.05e-2 | 0.153 | 0.107 | 0.023 | 0.042 | 0.083 | 0.007 | 0.012 | 0.001 | 0.99 | 1.00 |
| n32 L16 s1 | 1.61e-4 | 0.357 | 0.571 | 0.407 | 0.275 | 0.414 | 0.310 | 0.330 | 0.141 | 0.54 | 0.79 |
| n32 L8 s1 | 1.08e-3 | 0.201 | 0.355 | 0.122 | 0.121 | 0.202 | 0.042 | 0.192 | 0.026 | 0.80 | 0.97 |
| n32 L8 s2 | 2.10e-3 | 0.038 | 0.260 | 0.062 | 0.120 | 0.066 | 0.013 | 0.037 | 0.011 | 0.92 | 0.98 |
| n64 L16 s3 | 1.56e-3 | 0.037 | 0.238 | 0.062 | 0.143 | 0.052 | 0.027 | 0.032 | 0.024 | 0.95 | 0.96 |
| n64 L8 s1 | 5.45e-4 | 0.091 | 0.205 | 0.030 | 0.245 | 0.055 | 0.035 | 0.054 | 0.011 | 0.94 | 0.98 |
| **geometric mean** | | 0.135 | 0.29 | **0.085** | 0.12 | 0.134 | **0.042** | | | | |

- **Two rungs held out are 3.2x better than one rung held out** (geometric mean 0.042 against 0.134), and better on 8/8
  networks. At n = 64 the factor is 1.6-1.9, because one rung already explains 0.94-0.95 there.
- **The parameter-free Hermite T_2** (15e + 10 Lap e + Lap^2 e)/24 beats the parameter-free midpoint T1 by 1.6x in
  geometric mean. It wins on 5/8 networks and on all three where T1 is weak (n16 s2: 0.19 against 0.83).
- **The pooled two-rung OLS weights** are (beta_1, beta_2) = (0.825, -0.297), stable between 6 and 8 networks (0.838,
  -0.324). Hermite's are (0.667, -0.167).
- **Reading through Theorem A.**
  - Captures under the OLS weights: a_1 -> 1.06, a_2 -> 1.42, t_3 -> 1.78, t_4 -> 1.19. Under Hermite's all are 1,
    and t_4 is 0.67.
  - So the fitted weights over-correct every class with a negative rung-2 weight, relative to the single-readout
    model. The composed chain's error therefore carries rung-2 weights smaller in magnitude than a single readout's.
    One global pair (beta_1, beta_2) cannot say by how much per class.
  - This is the composition (transport-gate) correction, which enters through D^2 Phi[KX, KX] and Gamma_1(X, KX) in
    (3.1).
- **Caution.** The Gaussian-closure chain omits both open classes, which have weights (2, 2) and (1, -2). Two rungs
  separate those exactly, while one rung cannot. Part of the 3.2x is therefore open-class separation that the
  production chain does not need, since it carries a_1 (D3) and a_2 (V56's P4).

**The kappa_3 chain K3** (F1's dense third-cumulant chain; it carries kappa_3, births and triangles at birth, and
omits kappa_4 and the feedback).

Four distinct networks (n = 32, L = 8, data-file truth). F5's s1 and s2 are the same networks with another truth
sample, so they are not counted twice. n = 16 is unusable: K3's Hermite series overflow on dead neurons and give NaN.

| net | K3 MSE | T1 | H1 | **H2** | LOO1 | **LOO2** | oracle1 | oracle2 | X1 | X12 |
|---|---|---|---|---|---|---|---|---|---|---|
| n32 L8 s0 | 4.75e-4 | 1.712 | 0.163 | 0.046 | 0.138 | 0.042 | 0.137 | 0.026 | 0.81 | 0.97 |
| n32 L8 s1 | 1.43e-4 | 1.639 | 0.212 | 0.170 | 0.189 | 0.108 | 0.189 | 0.099 | 0.78 | 0.89 |
| n32 L8 s2 | 2.61e-4 | 1.704 | 0.477 | 0.226 | 0.444 | 0.191 | 0.442 | 0.185 | 0.41 | 0.75 |
| n32 L8 s3 | 1.23e-4 | 1.625 | 0.348 | 0.245 | 0.323 | 0.172 | 0.323 | 0.121 | 0.68 | 0.88 |
| **geometric mean / mean** | | 1.67 | 0.275 | **0.144** | 0.247 | **0.110** | | | 0.67 | 0.87 |

- **One rung.** The pooled coefficient on Delta_1 is -0.42, i.e. e - 0.21 delta: F1's K3 slope (-0.18 to -0.22) at
  n = 32. The parameter-free midpoint T1 fails (MSE ratio 1.6-1.7), because K3 needs about half of it. Hermite R = 1
  (the s^2 profile) is right to within the scatter: 0.275 against the OLS 0.247.
- **Two rungs held out: 2.2x lower residual than one rung** (0.110 against 0.247; better on 4/4). The demeaned
  explained share rises from 0.67 to 0.87, and the unexplained incoherent part falls from 0.33 to 0.13.
- **The parameter-free Hermite T_2 (0.144) beats the held-out one-rung fit (0.247)** with no fitted parameter.
- **The pooled two-rung weights (0.577, -0.112) are close to Hermite's (0.667, -0.167).** K3 carries kappa_3, so its
  open classes are not the issue. Its residual is the kappa_4 sector, with the path a_2 (1, -2) and the 4-cycle
  t_4 (0, -4), plus feedback. The second rung's gain here is what Theorem A predicts: the 4-cycle (and products), which
  rung 1 cannot see.

## 5. Closed walks as noncommutative moments: Toeplitz/Berezin deterministic equivalents

### 5.1 The chain's second chaos is a Toeplitz operator on the birth set

In the chain's source representation (note XLI),

    H_i = sum_m d_im |l_m><l_m| = Lambda D_i Lambda^T,     d_im = P_im w2_m,

where m runs over births (neurons of earlier layers), l_m is the birth's first-chaos leg, and P_im its mean-gate
transport to output i.

This is exactly stage-8 Theorem D's form T_f = sum_x f(x) |k_x><k_x|:
- a Toeplitz (anti-Wick) quantization of the per-output **symbol** d_i on the birth set;
- the **coherent states** are the legs l_m;
- the overlap Gram is G = Lambda^T Lambda.

The objects of the theory then read:
- **Open walks** are vector-state moments a_k = <L_i, T^k L_i>.
- **Closed walks** are quantum traces t_k = tr T_{d_i}^k = tr((D_i G)^k): sums over closed birth walks
  m_1 -> ... -> m_k of prod_j d_{m_j} G_{m_j m_{j+1}}.
- **The noncommutative pair.** D_i lies in the commutative diagonal algebra D of birth functions and G does not; the
  closed walk is a mixed moment of this non-commuting pair.
- **The cost.** For all outputs at once, the exact trace costs n^2 M per output (n^3 M per layer), since H_i must be
  formed. That is the K_4 wall in this coordinate system.

### 5.2 The semiclassical (overlap) expansion and its non-crossing part

Write G = G_0 + X with G_0 = diag(G) and X the off-diagonal overlaps. Expand tr((D(G_0 + X))^k) by the number and the
pattern of X factors.
- The X-free term is the **classical (upper-symbol) moment** sum_m (G_mm d_m)^k.
- The two-X terms are the first quantum correction. With the Berezin overlap O_mm' = G_mm'^2/(G_mm G_m'm'), the
  normalized coherent-state kernel, they involve the Berezin transform of the symbol,
  (B sigma)_m = sigma_m + sum_{m' != m} O_mm' sigma_m' with sigma_m = G_mm d_m.

**The surrogate.** Keep exactly the non-crossing pairings of X factors, contracted with the realized variance profile
X_mm'^2. This is the operator-valued semicircular (D-free) approximation, with amalgamation over the birth diagonal.

    t3_hat = sum_m sigma_m^3 + 3 sum_m sigma_m^2 (O sigma)_m,
    t4_hat = sum_m sigma_m^4 + 4 sum_m sigma_m^3 (O sigma)_m + 2 sum_{m != m'} sigma_m^2 sigma_m'^2 O_mm'
             + 2 sum_m sigma_m^2 ((O sigma)_m)^2,

where (O sigma)_m = sum_{m' != m} O_mm' sigma_m'.
- *What it drops.* The overlap cycles X_12 X_23 X_31 (and crossing pairings). These are the genuinely
  noncommutative K_4 contractions.
- *Status.* The expansion is an identity (proved by expanding the product). That the dropped part is small, or
  correlated with the kept part, is measured below.
- *Cost.* One product (G o G) D_i^T for all outputs: M^2 n, done blockwise by source. Within-source blocks cost
  n^3 per source-layer, about 0.5 unit.
- **The collective alternative.** Restrict to a k-dimensional collective subspace U of the input space: t_k inside U
  is tr((U^T H_i U)^k). This is Theorem A(5)'s operator-valued localization at U = P_k, where the rung-2 closed content
  is tr(P H P H^2). Cost: n k^2 M per layer.

### 5.3 Measured (`c1_closedwalk.py`, `c1_n{64,128,256}.txt`)

**Setup.**
- He networks, L = 12.
- The chain-style bare second chaos: first-chaos legs by mean-gate transport, w2 = phi(alpha)/sigma from the Gaussian
  closure, all earlier layers as births.
- Exact t_3 and t_4 for every output neuron, against the surrogates.
- "R^2" is the share of sum_i t_i^2 explained by c x surrogate (one global slope c per layer); "dm" is demeaned.
- "U(H2)" is the top-k eigenspace of sum_i H_i^2; "frac" is the projection coefficient (the share of the magnitude
  captured).

| n | layer | PR(sum H^2)/n | rms t_3 / 3a_1 | rms 3t_4 / 12a_2 | t_3 classical: R^2 (slope) | **t_3 D-semicircular: R^2 (slope)** | t_4 classical: R^2 | **t_4 D-semicircular: R^2 (slope)** | t_3 / t_4 inside U(H2), k = 16: R^2 (frac) |
|---|---|---|---|---|---|---|---|---|---|
| 64 | 6 | 0.207 | 0.129 | 0.084 | 0.853 (9.7) | **0.983 (1.29)** | 0.865 | **0.951 (1.84)** | 0.999 (0.92) / 0.998 (0.90) |
| 64 | 11 | 0.122 | 0.184 | 0.133 | 0.923 (17.9) | **0.991 (1.69)** | 0.949 | **0.982 (3.13)** | 1.000 (0.97) / 1.000 (0.98) |
| 128 | 3 | 0.396 | 0.086 | 0.059 | 0.885 (4.1) | **0.996 (1.02)** | 0.967 | **0.992 (1.07)** | 0.976 (0.42) / 0.938 (0.26) |
| 128 | 6 | 0.278 | 0.093 | 0.058 | 0.818 (8.5) | **0.987 (1.18)** | 0.830 | **0.932 (1.46)** | 0.991 (0.67) / 0.958 (0.59) |
| 128 | 11 | 0.173 | 0.125 | 0.075 | 0.934 (14.9) | **0.990 (1.39)** | 0.846 | **0.953 (2.24)** | 0.997 (0.81) / 0.994 (0.81) |
| 256 | 6 | 0.313 | 0.080 | 0.049 | 0.893 (8.1) | **0.997 (1.10)** | 0.940 | **0.975 (1.25)** | 0.971 (0.36) / 0.875 (0.25) |
| 256 | 11 | 0.214 | 0.102 | 0.057 | 0.920 (12.1) | **0.989 (1.25)** | 0.896 | **0.935 (1.84)** | 0.984 (0.58) / 0.955 (0.57) |

**What the measurements say.**
1. **The cheap NC surrogate tracks the exact per-neuron closed walks.** R^2 is 0.98-0.997 for the triangle and
   0.93-0.99 for the 4-cycle at every size and depth tested. The closed walks change sign across neurons (mean/rms
   0.01-0.10), so this is not a common scale.
2. **The dropped overlap cycles are correlated with the kept diagrams, and they shrink with width.**
   - The slope (exact over surrogate) at layer 6 is 1.29, 1.18, 1.10 for the triangle and 1.84, 1.46, 1.25 for the
     4-cycle (n = 64, 128, 256). The excess roughly halves per doubling, which extrapolates to 1.02-1.06 at n = 1024.
   - At depth the excess is larger, because the legs align with the collective modes. At layer 11 the slopes are
     1.69, 1.39, 1.25 for the triangle and 3.13, 2.24, 1.84 for the 4-cycle (n = 64, 128, 256). The excess shrinks by
     about 0.6 per doubling, which extrapolates to about 1.1 and 1.35 at n = 1024.
   - A per-layer slope is calibrated **truth-free**, offline, against exact closed walks: no ground truth enters.
3. **Inside-U at fixed k loses its share as n grows** (k = 16: frac 0.92 -> 0.67 -> 0.36 at layer 6). This is F1's
   partial-localization law again: the collective dimension of sum_i H_i^2 is about 0.12-0.31 n. Inside-U is not the
   scalable route; the D-semicircular one is.
4. **Affordable restrictions** (`c3_blocks.py`, n = 128). A chain can afford some restrictions and not others:

   | restriction | layer 6: R^2 (slope) | layer 11: R^2 (slope) |
   |---|---|---|
   | classical symbol, \|l\|^2 metric | 0.818 (8.5) | 0.934 (14.9) |
   | classical symbol, full-variance (dressed) metric | 0.823 (3.7) | 0.884 (5.3) |
   | pairs within a birth layer only | 0.947 (4.2) | 0.971 (6.6) |
   | all pairs | 0.987 (1.18) | 0.990 (1.39) |
   | births of the youngest 1 / 2 / 4 layers only: exact triangle | 0.12 / 0.42 / 0.78 | 0.17 / 0.13 / 0.51 |

   - **Old births build the closed walks of deep neurons.** A young-tier-only closed walk misses most of it. Any
     runtime surrogate must include the old tier: per-birth symbols for young sources, a small cubic tensor of the
     rank-r legs for old ones (section 6.2).
   - **Cross-layer pairs mostly set the scale.** Within-layer pairs keep R^2 high, and the per-layer slope absorbs the
     rest.
5. **Size.** In the bare representation the closed walks are 8-18% of the open-walk readout terms in amplitude,
   decreasing slowly with n (about n^{-0.3} at layer 6) and increasing with depth. Their relevance to the production
   residual is not established by these numbers. Note XLI found the dressed path class 4-8x larger than the bare one
   at depth, and the 4-cycle at 24-30% of the dressed path term at n = 1024. Section 5.4 tests the bare closed walks
   against true cumulants.

### 5.4 Are the closed walks real? (`c2_kappa_truth.py`)

Per-neuron true kappa_3 and kappa_4 of the pre-activations come from Monte Carlo (n = 64, L = 12, 2e7 inputs in two
independent halves, noise subtracted). They are compared with the bare second-chaos walk sums:
- kappa_3: open 3a_1 against open plus closed 3a_1 + t_3;
- kappa_4: 12a_2 against 12a_2 + 3t_4.

Residual = unexplained share of the true per-neuron cumulant signal (noise-corrected), with one free slope per
regressor; "both" fits the open and closed walks jointly (`c2_n64.txt`, `c2b_n64.txt`).

| layer | kappa_3: open 3a_1 (slope) | closed t_3 (slope) | both (open coef, closed coef) | dressed classical sum sigma^3 | D-semicircular | kappa_4: open 12a_2 | closed 3t_4 | both |
|---|---|---|---|---|---|---|---|---|
| 3 | 0.154 (2.6) | 0.170 (27) | 0.130 (1.5, 12) | 0.310 | 0.204 | 0.072 | 0.169 | 0.071 |
| 6 | 0.119 (3.8) | **0.048** (31) | 0.048 (-0.3, 33) | 0.118 | 0.046 | 0.293 | **0.208** | 0.194 |
| 9 | 0.114 (5.5) | **0.059** (35) | 0.058 (0.8, 30) | 0.222 | 0.053 | 0.114 | **0.075** | 0.068 |
| 11 | 0.041 (7.6) | **0.035** (42) | 0.035 (1.2, 35) | 0.056 | 0.041 | 0.056 | **0.050** | 0.050 |

The bare walks at weight 1 leave 0.44-0.78 of the signal. They are 3-8x (open) and 27-42x (closed) too small in
scale, because the bare legs miss the dressing (higher chaos in the legs; note XLI's 4-8x for the path class).

**The closed-walk shape is a real feature of the true per-neuron cumulants.**
- At layers 6-11 the exact closed walk alone is a *better* predictor of the true kappa_3 and kappa_4 than the open
  walk alone.
- Once it is in, the open walk adds nothing at layers 6 and 11.
- It carries information beyond the self-loops: the classical (self-loop) symbol alone is worse (0.06-0.31), and the
  pair-corrected surrogate recovers the full closed walk's power (0.041-0.053 at layers 6-11).

**The limit of this test.** The regressors are collinear through a shared per-neuron scale; every free-slope fit
explains 83-96%. So this test cannot apportion the residual of a chain that already carries these classes. That is
what the production-dump regression of section 7 (item 4) does.

## 6. Estimator components, costs and predictions

Costs are in units of 2n^3 at n = 1024. The base is V56 (208 units, C/B 0.203, raw 1.55e-8, adjusted 3.14e-9). The
late share s_late = 0.65 is the share of final MSE injected at layers 12-15 (DATA_FINDINGS).

### 6.1 The runtime two-rung defect merge (A2-late): costed, not recommended

**Definition.** e_A2 = e + beta_1 Delta1_hat + beta_2 Delta2_hat.
- Delta_r_hat comes from the telescoping (F1 Proposition 8 at r = 1, (3.1) at r = 2), truncated to the last four
  layers' local defects and injected in place (SYNTHESIS A-late).
- The covariance-map closed words are dropped (they are unaffordable, section 3.2).
- The readout closed walks come from the C-T overlap surrogate.
- beta is fitted on training networks (truth-fitted, two numbers per injection layer).

| item | units |
|---|---|
| rung 1: local defects plus first-chaos legs (A-late) | 20-50 |
| rung 2 open part: Gamma_1(X, Gamma_1) and Gamma_1(X, KX) hub families plus one solve, 4 layers | 30-50 |
| rung 2 readout closed part (C-T overlap surrogate), 4 layers | 4-16 |
| D^2-D^4 products, in-place transport | about 2 |
| **total** | **56-118 (C'/C = 1.27-1.57)** |

**Value.**
- The exact two-rung visibility is X*_12 = 0.58 (central mix, Corollary 4.2).
- Truncation keeps s_late = 0.65 of it.
- The analytic capture kappa_a is unknown (prior 0.6-0.8).
- About 1/3 of the closed content is lost with the covariance-map words (prior).
- So raw x (1 - 0.65 x 0.58 x 0.7 x 0.85), about x 0.78, and adjusted x 0.99-1.22.
- **Verdict: break-even at best.** The second rung's value is its closed-walk content (section 4.3), and Proposition
  3.5 says that content is better bought directly.

### 6.2 The proposed component: C-T, the overlap closed-walk (Toeplitz-Berezin) readout

**Which part of the closed walk is omitted.**
- The closed walk splits by the Toeplitz expansion of section 5.2: tr T_{d_i}^3 = (self-loops) + (overlap diagrams).
- *Self-loops* are the classical symbol sum_m sigma_im^3. They are each birth's own chi-square skewness, cube-
  transported. The production chain already carries them. Its sources start from the exact Gaussian kappa_3 of the
  relu at birth, which contains that skewness, and transport it on the A, A, P legs. With A_im ~ P_im w1 var_m, the
  leg product A^2 P is the self-loop's P^3 var^2 structure.
- *The overlap diagrams* are what is omitted:
  - pairs: 3 sum_m sigma_m^2 (O sigma)_m;
  - cycles X_12 X_23 X_31 between **distinct** births.

  The overlap triangle with two current-layer births and one carried source is note XLI's joint-gate triangle.
- C-T therefore inserts only the overlap part. Inserting the classical term would double-count the chain's own birth
  kappa_3. The decisive experiment checks this split (section 7, item 4).

**Definition.** Insert the overlap closed walks where V56 inserted the path class, at layers 3-15 (or 12-15) and at
weight 1. Weight 1 is their Edgeworth weight, which rung 2 certifies (Proposition 3.5). For each pre-activation i:

    kappa3_diag_i += c3_l . 3 sum_m sigma_im^2 (O sigma_i)_m,
    kappa4_diag_i += 3 c4_l [ 4 sum_m sigma_im^3 (O sigma_i)_m + 2 sum_{m != m'} sigma_im^2 sigma_im'^2 O_mm'
                              + 2 sum_m sigma_im^2 ((O sigma_i)_m)^2 ],

with:
- sigma_im = P_im w2_m var_m, the dressed symbol (the full birth variance makes rung 0 hold; F1 section 4, note XLI);
- O_mm' the squared dressed correlation between births:
  - within a source it is R_{b,mm'}^2 at the birth layer, which the chain has when the source is born; store R_b o R_b
    per young source;
  - across sources it would come from the Schur Gram At^T C^{-1} At (not costed into the base version);
- c_l the per-layer remainder factors (overlap cycles and cross-source pairs);
- the chain's own closure and transport then carry the terms, as they carry D3 and P4.

**Calibration (truth-free).** c3_l and c4_l are per-layer ratios of exact to surrogate overlap closed walks. They are
computed offline on training networks from the exact dressed walks tr((D_i G)^k) minus their self-loops,
G = At^T C^{-1} At. No ground truth enters.

**Cost.**

| piece | cost |
|---|---|
| (R_b o R_b) Sigma_s^T and (R_b o R_b)(Sigma_s^2)^T per young source-layer | 2n^3, about 1 unit |
| elementwise sums | about 0 |
| with about 4 young sources | about 4 units per layer |
| **C-T at layers 12-15** | **+16 units** (C/B 0.203 -> 0.219) |
| C-T at layers 3-15 | +52 units |
| cross-source and old-tier pairs through the Schur Gram | n^3 per source pair; left to the slopes, extended only if item 4 of section 7 shows they carry the signal |
| memory: R_b o R_b, n x n float32 per young source | 4 MB each |

**Small-n evidence** (bare metric, sections 5.3-5.4).
- The pair-corrected (D-semicircular) surrogate tracks the full per-neuron triangle at R^2 0.98-0.997.
- Within-birth-layer pairs alone: R^2 0.95-0.97, but with slopes 4-7 (cross-layer pairs set the scale).
- Young births alone fail (R^2 0.12-0.78).
- Against Monte Carlo cumulants (section 5.4), the closed walk's per-neuron shape predicts the true kappa_3 and
  kappa_4 at least as well as the bare open walk does, and better at layers 6-11. The self-loop alone does worse.
  The overlap diagrams carry real shape information.
- So the small-n evidence supports the surrogate algebra and shows that closed walks are real. It does not measure
  the overlap part's share of a carried chain's residual; only the experiment can.

**Prediction (conjecture; it hinges on f_cw, the overlap closed-walk share of the production residual, which
section 7 measures).**
- *The prior on f_cw* is 0.12 [0.03, 0.30]. It rests on:
  - W1's bound: the joint-gate (overlap) triangle class is at most about 1/3 of the MSE;
  - note XLI's 4-cycle at 24-30% of the dressed path term;
  - the remark that the remaining kappa_3 defect is closed walks.
- *Expected change.* With within-source fidelity about 0.8 after calibration, and about 0.8 of the overlap error
  removable at weight 1:
  - C-T at layers 12-15 (s_late 0.65): raw -5% [-1%, -13%] at +16 units, so adjusted -3% [+6%, -10%];
  - C-T at layers 3-15: raw -8% [-2%, -20%] at +52 units, so adjusted +17% to -10%. Not preferred unless f_cw >= 0.25.
- *Reading.* This is a modest component, priced honestly. The theory's contribution is to say exactly which object to
  add (the overlap Toeplitz diagrams, not the self-loops) and at which weight (1), and that no defect machinery is
  needed.

**Risks.**
1. f_cw may be at the low end: the chain's D21 hub may carry part of the overlap triangle implicitly.
2. The cross-source overlaps, outside the within-source surrogate, may carry most of the signal. c3 shows that
   cross-layer pairs set the scale. Then the cost rises toward n^3 per source pair.
3. The remainder factors c_l may vary across networks.

### 6.3 How much more of the error does a two-rung defect explain than one rung? (the requested prediction)

1. **Theory plus class model** (Theorem A, Corollary 4.2).
   - For the production residual, X*_1 = 0.23-0.42 and X*_12 = 0.54-0.63, an increment of +0.12 to +0.40.
   - At the central mix that is 0.28 -> 0.58, i.e. **2.1x**.
   - At least 85% of the increment is closed walks. Rung 2 restricted to open walks adds only +0.03 to +0.06.
2. **Measured on small networks.**
   - Gaussian-closure chain (8 networks): held-out residual 0.134 -> 0.042, 3.2x lower; demeaned explained share
     0.84 -> 0.95 on average. That chain's gain is inflated by open-class separation.
   - kappa_3 chain (4 networks, n = 32), the production-like case. Held-out residual 0.247 -> 0.110 (2.2x lower);
     demeaned explained share 0.67 -> 0.87; parameter-free Hermite T_2 at 0.144.
3. **Identity used by the experiment.** By Proposition 3.5 and Corollary 4.2,

       X*_12 - X*_1 = (0.85-1.0) f_cw^total + (0.03 to 0.06).

   - Here f_cw^total is the share of the residual in the closed-walk class: the omitted overlap diagrams plus any
     error in the carried self-loops, since both have closed-walk profiles.
   - Measuring it (section 7, item 4, truth-backed and offline) therefore measures the two-rung increment at n = 1024
     without the 2n^2 + 4n + 1 = 2.1e6 chain runs per network that Lap^2 would need.
4. **Answer.** A two-rung defect explains about **twice** the one-rung share: +0.15 to +0.35 of the per-neuron error
   at the central prior. Both small-n chains agree in direction and size: 2.2x (K3) and 3.2x (Gaussian closure) lower
   held-out residual. Almost all of the increment is closed-walk content.

## 7. The minimal decisive experiment (E-N3; offline on AWS, unbilled)

**Question.** What share f_cw of the production residual sits in the closed-walk class, the class that rung 1 cannot
see and rung 2 sees at weights -6 and -4? How is it split between omitted overlap diagrams and errors in the carried
self-loops? How faithful is the within-source Berezin surrogate?
- This decides C-T.
- By Proposition 3.5 and Corollary 4.2 it also fixes the two-rung increment, X*_12 - X*_1 = (0.85-1.0) f_cw^total +
  0.03-0.06, without the 2.1e6 chain runs per network that Lap^2 would need at n = 1024.

**Inputs.**
- Official networks 0-7 (`W_off{k}`).
- MC caches `mccache_{k}.npz`: per-layer, per-neuron power sums give the true kappa_3 and kappa_4 of every
  pre-activation (2^20 inputs; the kappa_3 noise is about 4e-3 sigma^3 per neuron). Use split halves for the noise if
  the cache has them, else a bootstrap over inputs of one recomputed layer.
- Chain dumps from `est_v29` / `estimator_final_v56` with counterterms on, with a dump patch at layers 3-15:
  - (a) the per-neuron kappa_3 and kappa_4 diagonals the readouts use (D3 input, g4row);
  - (b) per young source: P_s, w2_s, var_s at birth, R_b and the arms At_s;
  - (c) per old source: the rank-r legs;
  - (d) C_l;
  - (e) the per-layer injected mean error inj_l (dumps plus truth means).

**Computation**, per network, at layers l in {6, 9, 11, 12, 13, 14, 15}:
1. **Exact dressed closed walks, split.**
   - Form K_i = C_l^{-1/2} (sum_s At_s D_{s,i} At_s^T) C_l^{-1/2}. That is n^2 M per neuron, about 5n^4 per layer
     with M about 5n.
   - Compute t3_i = tr K_i^3 and t4_i = tr K_i^4.
   - Split them into self-loops (the classical symbol in the same metric) and overlap = total - self.
   - Split the overlap further into within-source and cross-source parts.
2. **Bare walks**, from first-chaos legs, for the dressed/bare ratio.
3. **The within-source Berezin surrogate** (section 6.2): R^2 against the exact overlap, and the per-layer factors c_l.
4. **Residual regressions.**
   - Regress Delta kappa_3 = kappa_3_true - kappa_3_chain on the overlap triangle and on the self-loop triangle
     separately. Do the same for Delta kappa_4 and the 4-cycle parts.
   - Report: noise-corrected removed signal variance, fitted slopes, and the removal at slope 1. The prediction is a
     slope near 1 on the overlap part and a small coefficient on the self-loop part, because the chain carries it.
5. **End-to-end oracle.** The production chain, free-running, on networks 0-7, with its kappa_3 and kappa_4 diagonals
   augmented at weight 1. A file-fed hook supplies:
   - (i) the exact overlap closed walks, at layers 3-15 and at 12-15;
   - (ii) the within-source surrogate times c_l, at 12-15.

   Report raw MSE and cost against the base. f_cw (overlap) is the raw reduction from (i) at layers 3-15.

**Compute.** Items 1-4: 8 networks x 7 layers x about 3 min on a 32-core instance, about 3 instance-hours. Item 5:
24 chain runs, under 1 instance-hour. Total under 5 instance-hours.

**Decision rule.**
- **Build** if the oracle (i) gives raw <= -8% at layers 12-15 (f_cw >= 0.12 effective), the surrogate (ii) gives
  >= 60% of the oracle's reduction, and c_l is stable within +-20% between networks 0-3 and 4-7:
  - implement V58_CT (within-source Berezin overlap readout, layers 12-15);
  - cold-screen it on 16 networks;
  - adopt by the note XXXIX protocol if adjusted MSE improves >= 3% paired, on >= 12/16;
  - then score it on 100 networks with refitted counterterms.
- **Extend to cross-source pairs** (Schur Gram, about n^3 per source pair) only if item 1 shows that the cross-source
  overlap carries more than half of the oracle signal and the oracle at 12-15 is <= -12%.
- **Drop C-T** if the oracle (i) at layers 3-15 is > -4% (f_cw < 0.04). Then:
  - the closed walks are not where the production residual lives;
  - the two-rung increment over one rung is <= 0.1;
  - beyond rung 1 the ladder survives as explanation only.

**Registered predictions.**
- *Item 4.* The overlap triangle and 4-cycle remove 10-30% of the chain's kappa_3 and kappa_4 residual signal at
  layers >= 9, with slope 0.7-1.4. The self-loop parts remove < 5%.
- *Item 5.* The oracle (i) gives raw -12% [-3%, -30%] at layers 3-15 and -8% [-2%, -20%] at 12-15. The surrogate (ii)
  gives 50-80% of the latter.
- *Ratios.* Dressed/bare overlap amplitude 2-6. Within-source share of the overlap 0.3-0.6, with c_l 1.5-3.
- *Two-rung increment.* X*_12 - X*_1 = f_cw^total + 0.03-0.06, i.e. 0.15-0.35, about 2x the one-rung share.

**A local follow-up** (optional; F1's truth): the two-rung K3 test at n = 48 (`b3_k3bilap.py`, 4801 runs per network,
about 1.5 core-hours each), to see whether the 2.2x of section 4.4 holds with width.

## 8. What the frame gives, and what it does not

**Gained (proved, and checked).**
- *Theorem A: the exact solvability of the leaf mixture.*
  - The localization filtration of the second-chaos moment problem is by degree, in both the vector state and the
    trace state.
  - Rung r sees exactly the walks with <= 2r letters, and each class has an explicit profile and rung weight.
  - The gain's two halves are identified with the a_2 and t_4 profiles.
  - **Every Gaussian localization, including a non-commuting operator-valued U, is first-order blind to the trace
    state.**
- *Theorem B (Gamma-calculus).* Rung r is graded by the iterated carre du champ Gamma_k, k <= r. Gamma_2, the
  Hilbert-Schmidt pairing of Hessians, is where the trace state of the chaos tuple enters. Rung 2 sees closed walks
  only through it, so it is never a cheaper route to them (no free lunch).
- *Theorem C (Hermite-Richardson).* Parameter-free multi-rung estimators with explicit weights, exact on all classes of
  bounded degree.
- *Measured.* Two rungs beat one, held out:
  - 2.2x on the kappa_3 chain (the production-like case);
  - 3.2x on the Gaussian-closure chain.

  The parameter-free Hermite T_2 beats every parameter-free and every fitted one-rung merge on the kappa_3 chain.
- *The Toeplitz/Berezin reading* of the chain's second chaos: H_i is the Toeplitz quantization of a birth symbol. Its
  closed walks have cheap noncommutative (operator-valued free, D-semicircular) deterministic equivalents with
  R^2 0.98-0.997, whose remainder (the overlap cycles) shrinks with width.

**Where genuine noncommutativity enters, precisely.**
1. *The vector state against the trace state of the second-chaos operator.* Commutative for one neuron; the split
   between rungs 1 and 2 is still exactly this distinction.
2. *Closed walks as mixed moments of the non-commuting pair* (the commutative birth-diagonal symbol algebra and the
   overlap Gram). This is what operator-valued free probability over the birth diagonal approximates:
   - Shlyakhtenko's D-semicircular elements with a variance profile;
   - deterministic equivalents of Gram matrices with a variance profile (Hachem-Loubaton-Najim).

   The approximation buys the closed walks at n M^2 (or n M) instead of n^2 M per output.
   - *The Toeplitz expansion also separates two parts.* The self-loops (the classical upper symbol) are the births'
     own skewness, which the chain already carries. The overlap diagrams (the first quantum correction, i.e. the
     Berezin transform of the symbol) are what it omits. The joint-gate triangle is an overlap diagram.
   - This split was found during the work (section 6.2). It changed the component: insert the overlaps, not the
     classical term.
3. *Covariance-map words* tr(H_a H_b H_a H_b) and tr(H_c H_a H_b). They are genuinely noncommutative and seen at rung 2,
   but unaffordable at n^2 pairs.
4. *Operator-valued localization U.* Its rung-2 closed content tr(U H U H^2) is the collective-subspace closed walk.
   That route loses its share as n grows, which is F1's partial-localization law.

**Not gained.**
- A cheaper defect, or a cheaper second rung.
- Truth-free merge weights: beta must be fitted with truth. Only the surrogate slopes are truth-free.
- Anything for the covariance-map words.
- The n = 1024 relevance of the closed walks, which is section 7's experiment.

**Novelty.**
- *Classical ingredients:* Gaussian log-det integrals (Imhof), Edgeworth expansions, the Bakry-Emery Gamma_2,
  Hermite/Richardson extrapolation (Talay-Tubaro), Berezin-Toeplitz semiclassics, operator-valued free probability
  (Speicher, Shlyakhtenko), deterministic equivalents.
- *New to our knowledge:*
  - the exact leaf-mixture CGF and the rung-walk visibility theorem for Gaussian-closure readouts;
  - the operator-valued first-order blindness to the trace state;
  - the Hermite-Richardson weights beta_r = (-1)^{r+1}(R+1-r)/((R+1) r!) for heat-defect hierarchies;
  - the Toeplitz reading of the chaos tuple, with a measured D-free closed-walk surrogate.
- *Literature check.* The Elicit RG/Hopf survey (`elicit/rg_hopf`) found no Connes-Kreimer counterterm prescription
  for finite cumulant closures, and states that "identities do not by themselves make the error computable". That is
  Proposition 3.5 in their language. The nc_heat session produced no content.

## Referee report

Adversarial referee, NCG round. Checks are in `loc/ref_n3/`:
- `r1_sym.py`: sympy check of Theorem A(2)-(4), Proposition 1.2 and Theorem C;
- `r2_num.py`, `r2_out.txt`: Theorem A(5) at n = 5;
- `r3_toep.py`, `r3_out.txt`: Toeplitz-overlap surrogate fidelity on synthetic Gram models;
- `r4_hybrid.py`, `r4_out.txt`: the spike-exact hybrid;
- `r5_arith.py`, `r5_out.txt`: visibility and adjusted-cost arithmetic.

**Verdict.** The mathematics is correct. The NC content is real but thin. The proposed component C-T does not survive
the cost arithmetic as written: its central prediction is a loss in adjusted MSE, not -3%.

### R1. What checks out (re-derived independently)

- **Theorem A(1)-(4).**
  - The series of the closed and open leaf-mixture CGF reproduces q_1 = 0, q_2 = x^3/2 + x^4/4 and
    q_3 = x^6/6 + x^5/2 - x^3/3, as well as p_1 and p_2.
  - Every class profile rho in the table is reproduced exactly, with weights t_3 (0, -6, -12), t_4 (0, -4, 0),
    t_5 (0, 0, 30), t_6 (0, 0, 12), a_1 (2, 2, 0), a_2 (1, -2, -6), a_3 (0, -6, -12) and a_4 (0, -2, 18).
  - The sign conventions agree: b^(r)(1) = -Delta_r, w_r = (-1)^{r+1} r! [u^r] rho, and b''(1) matches the Edgeworth
    weights.
- **Proposition 1.2.** Delta_1, Delta_2 and Delta_3 are reproduced from the forward-difference formula. The
  derivation of A^2 e(0) = (Lap^2 e + 2 Lap e)/4 is also correct.
- **Theorem C.** T_1, T_2 and T_3 are reproduced: T_3 = (105, 105, 21, 1)/192. All three return v(0) on the
  Euler-Stein sequence. Theorem C is two-point Hermite interpolation: zeros of order 2 at s = 0 and derivatives at
  s = 1. The weights are new in this context but elementary.
- **Theorem A(5).**
  - The O(U) closed coefficient is about 1e-12.
  - The O(U) open coefficient matches theta^3 L^T H U L + (theta^4/2) L^T H U H L to 10 digits.
  - With L = 0, the O(U^2) coefficient matches (theta^3/2) tr(UHUH^2) + (theta^4/4) tr(UH^2UH^2) to 5-6 digits.
- **Proposition 3.3 and Corollary 3.4.** The Gamma_2 contractions (2 g_mv t_3 + 2 g_vv t_4 = g''' t_3 + g'''' t_4/2)
  are consistent with b''(1).
- **Corollary 4.2.** X*_1 / X*_12 = 0.278 / 0.581 (central), 0.229 / 0.630 (closed-heavy) and 0.417 / 0.543
  (open-heavy) are reproduced.

### R2. Errors

**E1 (cost arithmetic; decisive for the component). The adjusted prediction for C-T has the wrong sign.**
- Adjusted = raw x C/B. At +16 units the bill goes from 208 to 224, a factor of 1.077.
- Raw -5% therefore gives adjusted x 0.95 x 1.077 = **x 1.023, i.e. +2.3%, not -3%.**
- The stated interval [+6%, -10%] is also wrong. The correct interval is [+6.6%, -6.3%].
- Break-even needs raw <= -7.1% at +16 units. The central prior (-5%) loses.
- C-T at layers 3-15 (raw -8% at +52 units) gives x 1.15.

**E2 (cost under-count).**
- The table bills "(R_b o R_b) Sigma_s^T and (R_b o R_b)(Sigma_s^2)^T per young source-layer" as "2n^3, about
  1 unit". These are two n x n x n products, i.e. **2 units**. Only the kappa_3 half needs just the first product.
- The production code confines sources older than AGE_OLD = 4 (estimator_final_v56.py:919). That leaves up to 5
  dense sources per layer, not "about 4".
- The kappa_3-plus-kappa_4 component at layers 12-15 is therefore +32-40 units, not +16. Break-even is then raw
  <= -13% to -16%, the optimistic end of the author's own interval. At the central raw -5% the result is
  adjusted x 1.10-1.12.

**E3 (the build rule contradicts the cost).**
- The rule "oracle <= -8% at 12-15 and surrogate >= 60% of it" admits a surrogate at raw -4.8%. At +16 units that is
  adjusted +2.5%, and at the corrected +32-40 units +10%.
- The rule must be stated in adjusted terms: surrogate raw <= -(Delta C/(C + Delta C)), i.e. <= -7% at +16 units and
  <= -13% at +32 units.

**E4 (internal inconsistency in the design).**
- Section 5.3 item 4 measures that young births alone fail: the exact triangle restricted to the youngest 1/2/4 birth
  layers has R^2 0.12/0.42/0.78 at layer 6 and 0.17/0.13/0.51 at layer 11. It concludes that "any runtime surrogate
  must include the old tier".
- C-T (section 6.2) nevertheless uses only within-young-source pairs and leaves old-tier and cross-source pairs "to
  the slopes c_l".
- A per-layer scalar cannot repair a shape mismatch of R^2 0.13-0.51. The design rests on the case the author's own
  data reject.

### R3. Overclaims

**O1. The surrogate's R^2 0.98-0.997 is not evidence of NC surrogate fidelity.**
- It is R^2 of the *total* t_3, self-loops included, and the self-loops alone already give 0.82-0.93. R^2 is
  dominated by a shared per-neuron scale.
- In `r3_toep.py`, with symbols whose signs are coherent (75% one sign, as for mean-gate transport), even total ~ self
  reaches R^2 0.975-0.994 at a wrong slope of 12.
- The quantity C-T inserts is the *overlap* part (total - self). Its fidelity against the within-source surrogate
  (the inserted object) was never measured. The informative statistic is the slope, and the slopes are far from 1:
  within-layer pairs have slope 4.2-6.6 in section 5.3's own table.

**O2. "The slope tends to 1 with width" is a property of the small-n regime, not a law.**
- In the author's data PR(sum H^2)/n *rises* with n (0.12, 0.17, 0.21 at layer 11), so the collective fraction falls.
- At n = 1024 the production covariance is strongly collective: PR/n is 0.087 at layer 10 and 0.052 at layer 15, and
  the top 128 modes carry 89% (DATA_FINDINGS).
- In a spiked Gram model (8 collective modes carrying about 70% of each leg's energy; `r3_toep.py`), the dropped
  overlap cycles *dominate* the pairs (rms cycle/pair = 1.6, 3.2, 6.1 at n = 64, 128, 256). The D-semicircular slope
  *grows* with n (2.6, 4.0, 7.0) while R^2 stays 0.97-0.99.
- Low-rank collective structure is exactly what makes odd overlap moments coherent. Free (D-semicircular) asymptotics
  is the wrong model for the late layers, where 2/3 of the MSE is injected. The extrapolations 1.02-1.06 and 1.1-1.35
  at n = 1024 are unsupported.

**O3. "Weight 1 is certified by rung 2" adds nothing.**
- Weight 1 is just the cumulant formula kappa_3 = 3a_1 + t_3, kappa_4 = 12a_2 + 3t_4. Note XLI already states it
  (chaos-grading README lines 63, 74, 524). It already names the omitted piece as "the joint-gate triangle".
- The ladder certifies nothing that the Edgeworth expansion did not already give. The author half-concedes this ("no
  defect machinery is needed"). The summary still frames C-T as a product of the ladder.

**O4. "Weight 1" against the truth is contradicted by section 5.4.**
- The bare closed walks predict true kappa_3 and kappa_4 only at slope 27-42.
- The predicted dressing factor of 2-6 leaves a 5-20x scale gap. So either the truth's closed-walk-shaped content is
  not the chain-model's walk at weight 1, or the shared-scale collinearity (which the author concedes) makes the test
  uninformative.
- In either case the "truth-free" calibration of c_l calibrates to a model quantity, the Gaussian-equivalent walks.
  The weight on the truth is not established.

**O5. Proposition 3.5(2), "rung 2 cannot see the closed walks more cheaply than they can be computed", is labelled
proved but is not a theorem.**
- Item 1 is a correct structural identity.
- Item 2 is a complexity claim with no lower-bound argument. Finite-difference or probe estimates of Lap^2 e never
  form tr H^3, and randomized trace estimation of tr H_i^3 costs about 3p n^2 M per layer against n^3 M exact.
- The conclusion is probably right in practice, because the probes are too noisy (F1). It should be labelled a
  heuristic.

**O6. The "identity" X*_12 - X*_1 = (0.85-1.0) f_cw^total + 0.03-0.06 is a model statement.**
- It rests on single-readout weights and uncorrelated classes.
- The author's own Gaussian-closure-chain fit shows that composed chains carry rung-2 weights of smaller magnitude
  than a single readout (OLS over-corrects every negative-weight class: captures 1.42-1.78).
- So measuring f_cw does not "measure the two-rung increment at n = 1024".
- Corollary 4.2's X* values are also fixed by invented class mixes. The production a_1 and a_2 classes are carried by
  the chain (D3, P4), and their residual (dressing and scale errors) need not have the omitted-class weights (2, 2)
  and (1, -2).

**O7. Hidden assumption: "the chain already carries the self-loops".**
- The argument "A^2 P is the self-loop's P^3 var^2 structure" conflates the c_1^2 c_2 Hermite structure (the open walk
  a_1) with the c_2^3 structure (the self-loop sigma^3).
- Note XLI line 63 says D3 = 3 L^T H L and that tr H^3 is omitted, which by its plain reading includes the self-loops.
- If the self-loops are not carried, C-T omits the larger half of the closed walk: the self term dominates the
  synthetic overlap only when the signs are incoherent. The experiment's item 4 tests this, but the design should not
  presuppose the answer.

**O8. Metric mismatch.**
- Within-source O = R_b o R_b is the birth-layer correlation. The calibration target G = At^T C_l^-1 At is the
  current layer's Schur view of the births, which is exact only for the first chaos (note XLI).
- c_l therefore absorbs a metric change as well as the remainder diagrams, and its stability across networks is less
  likely.

### R4. Noncommutative content: what is genuine

1. **Theorem A(5)** is the one genuinely new and useful statement.
   - The full matrix-valued heat defect D = d_Sigma E - (1/2) Hess_m E at the base point, not just its trace,
     contains no trace-state (closed-walk) term for a Gaussian-closure readout.
   - Hence no anisotropic, collective-mode or operator-valued first-order localization (the F3/N2 routes with
     U = P_k) can see closed walks.
   - It is a clean negative result that prunes a branch. It is first-order in U and holds for deterministic U only;
     adaptive Chen-Eldan controls are not covered by "every localization scheme".
2. **Vector state against trace state.** Correct and illuminating. For one neuron it is two commutative spectral
   measures of one self-adjoint H, so it is a grading, not noncommutativity.
3. **The Toeplitz/D-semicircular surrogate.**
   - The computation is "keep the terms of tr((D(G_0 + X))^k) with at most two off-diagonal Gram factors (plus
     X^2 X^2)".
   - The operator-valued semicircular language is accurate in that odd X-moments vanish, but it changes no
     computation beyond that truncation.
   - Small caveat: the t_4 surrogate double-counts the X_ab^4 term, which appears in both the a = c and the b = d
     patterns. It needs a -sum_{a != b} sigma_a^2 sigma_b^2 O_ab^2 correction.
4. **Covariance-map words** are genuinely noncommutative but buy nothing: the author agrees they are unaffordable.

**Net.** One genuine NC-flavoured theorem (A5, negative). One correct but mostly relabelled surrogate. No computation
is changed in a way that helps the production chain.

### R5. Strongest surviving idea, corrected

The real question is a deterministic equivalent for per-neuron closed walks in a *spiked* Gram. The correct
noncommutative tool is not plain D-semicircular. It is a finite-rank (outlier) part treated exactly plus the free
(D-semicircular) bulk:

    t3_hybrid = tr((A D_i A^T)^3) + [semi_3(G) - semi_3(A^T A)],     A = P_k^T Lambda (top-k collective basis of the legs).

- **Synthetic check** (`r4_hybrid.py`, same spiked model): the slope against the exact t_3 is 0.987, 0.991, 0.994 at
  n = 64, 128, 256, with R^2 1.000. Plain semicircular gives 2.6, 4.0, 7.0.
- **Cost.** The k x k compressions for all outputs need about k^2 M n multiply-adds, using symmetry k(k+1)/2.
  - At k = 32 and M = 5n that is about 2.6 units per layer.
  - It reuses the collective basis that DATA_FINDINGS shows carries 89% of the late variance.
- Together these give a closed-walk estimator whose remainder does not grow with collectivity.

Even so, C-T's value is capped by the kappa_3 and kappa_4 diagonal oracles.

### R6. The decisive experiment, corrected (cheaper and sharper)

- **Step 0 (missing, about 1 instance-hour).**
  - Run the existing oracle hook V37_ORACLE = D3 with V37_ORACLE_LAYERS = 12,13,14,15 and MC-truth kappa_3, plus the
    analogous g4row oracle, on networks 0-7.
  - This is the ceiling for *any* diagonal-cumulant augmentation, C-T included. Note XXXI (D3 + D21 at every layer,
    about 40%) suggests that the D3-only ceiling at 12-15 is perhaps 10-20%.
  - **Kill C-T** if the D3 + g4row oracle at 12-15 gives raw > -15%. Closed walks are a fraction of that ceiling,
    and break-even at the corrected cost is -13% to -16%.
- **Item 1 must report what C-T inserts.** That is overlap-only fidelity (slope and demeaned R^2 of exact overlap
  against within-young-source surrogate), plus the young/old and within/cross split. Total-t_3 R^2 is not enough.
- **Add the hybrid (R5)** at k = 32 and 64 as a surrogate arm.
- **State the decision rule in adjusted terms** with the corrected bill (E1-E3).

### Summary of the verdict

- **Theorems A, C and Proposition 1.2: correct and independently verified.**
- **Theorem B:** the identity part is correct; the "no free lunch" complexity claim is a heuristic.
- **Measured small-n gains of two rungs:** plausible, but irrelevant to runtime, which is offline-only at n = 1024.
- **C-T:** mis-costed by about 2-2.5x. Its central prediction is an adjusted *loss* of +2% to +12%, not -3%. It rests on
  the young-only case the author's own data reject, and on an R^2 that does not measure what is inserted.
- **NC content:** A5 is genuine and useful as a pruning result. The Toeplitz surrogate is a relabelled truncation, and
  it is the wrong free model in the collective regime.
- **Priority for Phase 2:** low. Run Step 0 first; most likely it closes C-T.
