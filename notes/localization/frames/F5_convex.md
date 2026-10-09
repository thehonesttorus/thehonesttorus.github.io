# F5. Convex order, Berezin-Lieb sandwiches, and the localization-Richardson midpoint

Frame F5 of the localization round. Question: what do convex-order relations (martingale couplings, Berezin-Lieb
sandwiches, the Ising <= boson <= Ising chain of stage-8 Theorem D) say about E[x_L] for the Phase-2 MLP, and what
estimator do they imply. Everything below is either proved here (with the proof), marked as a sketch or a conjecture,
or measured by a script in `loc/f5/` or `loc/` (named at each result). Small-width checks only; heavy work is specified in
section 7 for the AWS fleet.

## 0. Verdict and summary

**Verdict: marginal.**
- Convex-order sandwiches themselves are closed for Phase 2. They are proved, and measured to be 10^2-10^3 too wide or
  identically Jensen.
- One construction the frame forced out is new and strong: Richardson extrapolation in the localization parameter,
  done deterministically from the estimator's own derivatives.
  - It removes 90-97% of a Gaussian closure's MSE on 5 of 6 networks at widths 64-128 (51% on the sixth), with no
    fitted parameter. At n = 128 the derived coefficient is confirmed to 3-5%.
  - It cannot be bought inside the budget for the production chain. Its exact form costs what the chain's missing
    closed-walk classes cost, n^4, and its cheap part is the weak part.
  - Its value for Phase 2 is as an offline exactness auditor for the chain (section 6.4), decided by experiments E1a
    and E1b (section 7).

**What the frame closes (all proved; numbers measured):**
1. *The mean is invisible to martingale couplings.* Any two-sided partition-function sandwich, including the
   Ising <= boson <= Ising chain of Theorem D, pins the means together. Convex order reaches E F only through a single
   relu applied to an ordered pre-activation (Thm 2.1).
2. *One-layer sandwiches do not propagate.* E F_i is not convex-order monotone in the input law for depth >= 2,
   because H_i = E grad^2 F_i is indefinite. Measured: 52% of its eigenvalues are negative, and a convex-order dilation
   lowers E F_i by 0.019 (Thm 2.3).
3. *Matched closures are never comparable* with the truth (Prop 2.4).
4. *Sharpest moment sandwiches* (Markov-Krein LP) are 50-1500x the target width with 4-8 exact moments. They close only
   if ~78% of the variance is a certified independent Gaussian component (Thm 2.5).
5. *Spherical Berezin-Lieb*, the exact transplant of Theorem D to S^(n-1), where homogeneous networks live. The Berezin
   multipliers are beta_2(l) ~ 2l/n, with beta_2(1) = 2/(n+2) exactly and beta_odd = 0. At affordable levels the lower
   symbol is Jensen to 1e-3 of the variance (Thm 2.6).
6. *Random-leaf localization and Richardson die on variance.* The floor is V_2/K with V_2 = 0.0108 at depth 16, i.e.
   >= 1e6 chain runs (Thm 5.1). Deterministic subspace localization buys ~0.3% per direction at n = 1024.

**What has a definite sign:** the gain mode (Thm 3.1). relu commutes with a scale mixture, so every matched Gaussian
closure overestimates under it; its kappa_3 + kappa_4 template correction is negative at every alpha. Measured: the
Gaussian closure overestimates the post-activation mean on 84-100% of neurons at every layer >= 2 (n = 16 and 64). The
production chain carries the gain, so its residual should have no sign (prediction P-F5-3).

**What localization gives positively:**
- **The bias of any Gaussian-input estimator is its drift along stochastic localization** (Thm 4.2, Dynkin): the
  localization-integrated failure of the heat equation. The exact functional is the Doob martingale.
- Expanding the localization curve at t = 0 gives the **Euler-Stein tower**:
  - T_1 = (C + Delta_a C)/2 is the midpoint of the sandwich formed by the chain and its Euler-Stein companion;
  - T_2 = (3C + 6 Delta C + Delta^2 C)/8.
- The defect calculus (Thm 5.3, Prop 5.4) shows that first-order localization sees exactly the open walks of the second
  chaos, with class-dependent exponents (p = 2 for the kappa_4 path class, p = 4 for the kappa_3 path class). It
  rederives note XLI's path class and its weight with no regression.
- Measured on 12 networks (n = 16-128): T_1 cuts the closure MSE by a geometric mean of 8.7x. At n = 128 the cut is
  11-37x, with correlation 0.95-0.98 between drift and error and fitted coefficient -0.95/-0.97 against the derived -1.

**Cost at n = 1024.** Exact T_1 by finite differences takes 2n+1 chain runs: offline only. The structural drift needs
the covariance-tangent Grams, which are n^4. The local-tangent truncation (~+145 units) keeps the shape but delivers
only 60-87% of the closure MSE cut, against 90-97% exact.

**Decisive experiments:**
- E1a, a Gaussian-closure drift at n = 1024 (~6 core-hours per network). It predicts that T_1 removes >= 90% and reaches
  raw ~5e-8.
- E1b, the production chain's drift (2049 float64 runs per network). Is the production residual drift-visible?


---

## 1. The dictionary

**Setting.** x_0 ~ gamma = N(0, I_n); z_l = x_(l-1) W_(l-1), x_l = relu(z_l), l = 1..L; F(x_0) = x_L. Every z_(l,j) and
x_(l,j) is a positively 1-homogeneous, piecewise-linear function of x_0. Target m* = E F. Convex order: X <=_cx Y iff
E phi(X) <= E phi(Y) for every convex phi; by Strassen, iff a coupling with E[Y | X] = X exists (a martingale coupling).

**(D1) Gaussian stochastic localization of the input is the Mehler semigroup.** With observation Y_t = t X + B_t, the
posterior is N(m_t, s_t^2 I), m_t = Y_t/(1+t), s_t^2 = 1/(1+t). Write X = m_t + s_t g with m_t ~ N(0, (1 - s_t^2) I),
g independent. For f in L^2(gamma),

    E[f(X) | Y_t] = (P_tau f)(X~),   X~ = m_t / sqrt(1 - s_t^2) ~ N(0, I),   e^(-tau) = sqrt(1 - s_t^2),

where P_tau is the Ornstein-Uhlenbeck (Mehler) semigroup. *Proof:* X = e^(-tau) X~ + sqrt(1 - e^(-2 tau)) g. P_tau multiplies
the q-th Wiener chaos by e^(-q tau), so the localization martingale M_t = E[F_i | Y_t] has

    Var(M_t) = sum_(q>=1) (1 - s_t^2)^q V_q,                                                        (1.1)

V_q the chaos-q variance of F_i. The *lower symbol* (Berezin transform) of the network at "time" tau is P_tau F; the
*upper symbol* is F itself; P_tau F <=_cx F in law, increasingly in t.

**(D2) The readout is a support function.** For a random vector x, h_x(w) := E relu(w . x) = sup{ w . E[h(x) x] :
0 <= h <= 1 }, attained at h = 1{w . x > 0}. So each layer's mean vector is the support function of the convex body
K(x_(l-1)) = { E[h x_(l-1)] } evaluated at the columns of W_(l-1): mu_l = h_(x_(l-1))(W_(l-1)). K(xA) = K(x)A, and
h_x(w) - h_x(-w) = w . E x. With an affine term, { E relu(u_0 + u . X) } is the support function of Koshevoy-Mosler's
lift zonoid, and inclusion of lift zonoids is the linear convex order. A bias-free layer only sees u_0 = 0.

**(D3) Chaos and spheres.** By homogeneity F(x) = |x| f(theta), with |x| independent of theta = x/|x|. The degree-q
spherical harmonics of f and the q-th Wiener chaos of F are the same object up to radial factors. The natural home of a
Berezin-Toeplitz quantization of the MLP is therefore the sphere S^(n-1) (section 2.5), the continuum analogue of the
hypercube in Theorem D.

---

## 2. What convex order cannot do

### 2.1 The mean is invisible to martingale couplings (and to every partition-function sandwich)

**Theorem 2.1.** (i) For any localization (M_t) of F_i, E M_t = E F_i, and Z_s(beta) := E e^(beta M_s) <= Z_t(beta) for
s <= t and every real beta. (ii) If two moment generating functions satisfy Z_1 <= Z_2 on a neighbourhood of 0, then
Z_1'(0) = Z_2'(0).
Consequently every two-sided sandwich Z_low <= Z_F <= Z_up that is valid for beta of both signs pins the means together:
E_low = E F = E_up. This holds in particular for every sandwich realised by two localizations of one state, which is how
Theorem D's chain Z_Ising(beta_2 t) <= (2^n/d_l) Tr e^(t T_f) <= Z_Ising(t) arises. Such sandwiches bound fluctuations
(variance, log-MGF), never the mean.
*Proof.* (i) e^(beta x) is convex. (ii) h = Z_2 - Z_1 >= 0 with h(0) = 0 has an interior minimum at 0, so h'(0) = 0.

So the mean enters a convex-order argument only through a non-affine function applied to a convex-ordered argument. In a
bias-free ReLU network that is one relu at one layer: E relu(z) with z ordered.

**The Ising <= boson <= Ising question, answered.** The output's partition function Z_i(beta) = E e^(beta F_i) does have
the sandwich (section 2.5, with phi = exp): E e^(beta B f_i) <= (1/d_l) Tr e^(beta T_(f_i)) <= E e^(beta f_i) on the
sphere. All three share their first derivative at 0 by Theorem 2.1, so the chain bounds Var F_i and higher cumulants,
never E F_i. A one-sided bound on E relu(z) from an MGF (relu(z) <= e^(lambda z - 1)/lambda) has an O(sigma) gap. No
partition-function chain bounds E F usefully.

### 2.2 One-layer sandwiches exist; they do not propagate

**Proposition 2.2 (one layer).** If X <=_cx Y in R^n (or only in the linear convex order), then E relu(w . X) <=
E relu(w . Y) for every w simultaneously (relu(w . x) is convex in x). In zonoid terms, K(X) is contained in K(Y).

**Theorem 2.3 (no propagation through signed layers).** For depth >= 2, law(x_0) -> E F_i is not monotone in convex
order. For Sigma_eps = I + eps u u^T (so N(0, Sigma_eps) >=_cx N(0, I)),

    E F_i(N(0, Sigma_eps)) - E F_i(N(0, I)) = (eps/2) u^T H_i u + O(eps^2),   H_i = E grad^2 F_i = E[F_i (x x^T - I)],

which is negative for u along a negative eigenvector of H_i. For one layer H_i = p_i(0) w_i w_i^T >= 0. For depth >= 2,
H_i is indefinite. The mechanism: linear maps preserve convex order; relu maps convex order to the increasing convex
order (means rise); a negative weight destroys the increasing convex order.
*Proof.* Price's theorem / the Gaussian heat equation, d/d(Sigma) E f(Sigma^(1/2) g) = (1/2) E grad^2 f.
*Measured* (`f5/exp4_identities.py` (b); n = 16, L = 4, 2e6 samples):
- the fraction of negative eigenvalues of H_i is 0.52 on average (min 0.31);
- dilating the input by eps = 0.5 along the most negative eigenvector of H_1 lowers E F_1 from 0.5440 by 0.0186
  (second-order prediction 0.0209), although the dilated input is larger in convex order.

### 2.3 No moment-matched closure is comparable with the truth

**Proposition 2.4.** If X <=_cx Y and E X^2 = E Y^2 < inf, then X =_d Y.
*Proof.* With the martingale coupling Y = X + D, E[D | X] = 0, we get E Y^2 = E X^2 + E D^2, so D = 0.

Every closure that matches the pre-activation mean and variance is therefore either exact or convex-incomparable with
the truth. That covers the Gaussian closure, every Edgeworth-corrected law, and every localization mixture
Q_t = int N(m_t(Y), v_t(Y)) dP(Y), which matches by total variance. The sign of a closure's readout error is set by its
omitted cumulants through the Edgeworth coefficients

    E relu^((m+2))(N(mu, sigma^2)) = (-1)^m sigma^(-1-m) He_m(alpha) phi(alpha),   alpha = mu/sigma,          (2.1)

which change sign with alpha (He_1 at 0, He_2 at +-1). Section 3 isolates the one combination with a definite sign.

### 2.4 The sharpest moment sandwiches are 10^2-10^3 times too wide

**Theorem 2.5 (Markov-Krein extremal problem for the ReLU readout).** For a moment vector in the interior of the moment
cone of order 2k, the set of values of E relu(z) over laws with those moments is a closed interval [L_k, U_k]. Its
endpoints are attained by laws with at most 2k+1 atoms (Richter-Rogosinski / Caratheodory; Karlin-Studden). The same
holds for laws of the form z = M + sqrt(v0) xi with xi ~ N(0,1) independent of M, a certified Gaussian component of
variance v0, with relu replaced by its heat-smoothed version g_v0(m) = sqrt(v0) G(m/sqrt(v0)).
*Proof.* Linear programming over measures with finitely many moment constraints; the extreme points of the feasible set
have at most (#constraints) atoms.

*Measured* (`exp2_mk.py`, LP on a 4801-point grid). Widths are in units of the Phase-2 per-neuron target
rms 1.2e-4. Laws are n=1024-like: sigma = 0.4, Gaussian scale (gain) mixtures.

| law | Gaussian-closure error | width, 4 moments | 6 moments | 8 moments |
|---|---|---|---|---|
| gain mixture Var log G = 0.01, alpha = 0 | +2.0e-4 | 699x | 516x | 373x |
| same, alpha = 1 | +2.4e-4 | 418x | 298x | 221x |
| same, alpha = 2 | +1.4e-4 | 84x | 53x | 51x |
| Var log G = 0.04, alpha = 1 | +9.7e-4 | 420x | 310x | 231x |

- On the final layer of an n = 16, L = 8 network with its exact Monte Carlo moments, the median widths are 0.17
  (2 moments), 0.074 (4) and 0.058 (6), and the truth always lies inside.
- The sandwich is looser than the uncertified Gaussian guess by 30-1000x. The chain's precision does not come from
  moment information; it comes from regularity.

**Regularity is what closes the sandwich.** With a certified independent Gaussian component (two-point gain law,
alpha = 1, at most 0.79 of the variance certifiable):

| certified fraction f = v0/Var z | 0 | 0.25 | 0.5 | 0.7 | 0.78 |
|---|---|---|---|---|---|
| width, 4 moments | 420x | 141x | 47x | 16.7x | 10.4x |
| width, 6 moments | 308x | 58x | 11.7x | 2.6x | 1.2x |

Empirically W_k(f) ~ W_k(1/2) ((1-f)/f)^(k-1). In the language of quantum optics this is the statement that the
law's Glauber-Sudarshan P-function must stay a positive measure after removing most of the variance (a large
"nonclassical depth" margin). The chain can certify neither the exact cumulants (its bulk kappa_3 is about 10% off) nor
the divisibility. **Certified sandwiches are not a route to 1e-8.**

### 2.5 The spherical Berezin-Lieb sandwich: Theorem D transplanted, and why it is Jensen at feasible levels

**Theorem 2.6.** Let H_l be the degree-l spherical harmonics on S^(n-1), d_l = dim H_l, Z_l(theta, eta) the zonal
reproducing kernel (Z_l(theta, theta) = d_l by the addition theorem), k_theta = Z_l(theta, .). For f in C(S^(n-1)) let
T_f = P_l M_f P_l = int f(theta) |k_theta><k_theta| dsigma(theta) and Bf(eta) = <k_eta, T_f k_eta>/d_l. Then:
1. B is the Funk-Hecke operator with kernel d_l P_l(t)^2 (P_l the Gegenbauer polynomial normalised by P_l(1) = 1), a
   Markov operator. Its eigenvalue on H_q is beta_q(l) = d_l int P_l^2 P_q dmu_n, which vanishes for odd q and for q > 2l.
2. For every convex phi: int phi(Bf) dsigma <= (1/d_l) Tr phi(T_f) <= int phi(f) dsigma.
3. The spectral variance of T_f is (1/d_l) Tr T_f^2 - ((1/d_l) Tr T_f)^2 = sum_(q even, 2 <= q <= 2l) beta_q(l) ||f_q||^2.
4. For homogeneous networks and positively homogeneous phi (relu):
   E relu(z_(L,i)) = E|x| int relu(f_i) >= E|x| (1/d_l) Tr relu(T_(f_i)) >= E|x| int relu(B f_i).

*Proof.*
- (1) The Funk-Hecke theorem, plus P_l^2 even and orthogonality of P_q to polynomials of lower degree.
- (2) Lieb-Berezin (Lieb 1973; Simon 1980). For the upper bound, take an eigenbasis v_k of T_f: lambda_k = int f p_k,
  with p_k = |<k_theta, v_k>|^2 a probability density because int |k_theta><k_theta| = P_l. Jensen in each p_k and
  sum_k p_k = Z_l(theta, theta) = d_l give it. For the lower bound, apply Jensen to the spectral measure of T_f in the
  unit vector k_eta / sqrt(d_l), then integrate in eta.
- (3) (1/d_l) Tr T_f^2 = int int f(theta) f(eta) Z_l(theta, eta)^2 / d_l = <f, Bf>.
- (4) |x| is independent of theta.

This is the exact analogue of Theorem D: the hypercube's beta_q = C(q,q/2) C(n-q, l-q/2)/C(n,l) becomes the Gegenbauer
linearisation coefficient.

*Measured* (`f5/berezin_sphere.py`, exact Gauss-Jacobi quadrature). beta_2(1) = 2/(n+2) exactly at n = 8, 64, 1024,
and beta_2(l) ~ 2l/n for l << n:

| n = 1024 | d_l | beta_2(l) | beta_4(l) | Toeplitz spectral variance / Var f (depth-16 chaos profile) |
|---|---|---|---|---|
| l = 1 | 1,024 | 0.00195 | 0 | 3.0e-4 |
| l = 2 | 524,799 | 0.00389 | 1.1e-5 | 6.0e-4 |
| l = 3 | 1.8e8 | 0.00581 | 3.4e-5 | 9.0e-4 |

Two consequences:
- The odd chaos, including the first (19% of the output variance at depth 16), is invisible to every Berezin transform.
- The even chaos is seen with weight about 2l/n.

At any l with d_l affordable, the Berezin-Lieb lower bound equals Jensen, relu(E z), to within 1e-3 of the variance.
The sandwich closes only in the Szego regime l ~ n, where d_l ~ 2^n. The Kikuchi refutation regime l = n^delta works for
spectral norms, which need only the top of the spectrum; an expectation at relative precision 1e-4 needs the whole
spectral law. The small-n check in `f5/exp4ac_out.txt` (n = 8, L = 2, kinked neurons) confirms the ordering and that
the spectral variance is a few percent of Var f.

Kinked neurons of an n = 8, L = 2 network (P(f > 0) between 0.33 and 0.56), 2e6 sphere samples:

| neuron | int relu(Bf) | (1/d_1) Tr relu(T) | (1/d_2) Tr relu(T) | truth int relu(f) | Jensen | spectral var / Var f (l = 1, 2) |
|---|---|---|---|---|---|---|
| 6 | 0.1086 | 0.1114 | 0.1148 | 0.2710 | 0.1086 | 0.027, 0.039 |
| 1 | 0.1355 | 0.1423 | 0.1464 | 0.2541 | 0.1350 | 0.055, 0.081 |
| 5 | 0.0633 | 0.0790 | 0.0834 | 0.1960 | 0.0617 | 0.048, 0.070 |

The sandwich holds and is useless: the Toeplitz value moves off Jensen by 2-10% of the gap to the truth. The Monte
Carlo beta_2(1) is 0.2006, against an exact 0.2000.

---

## 3. What does have a definite sign: the gain-mode rule

**Theorem 3.1 (relu commutes with the gain).** Let z = A w with A >= 0 independent of w. Then:
1. Equivariance: E relu(A w) = E[A] E relu(w), exactly, for every law of w.
2. If w ~ N(mu_w, sigma_w^2) and N_z denotes the Gaussian with the mean and variance of z, then

       E relu(N_z) - E relu(z) = g(mu, v + delta) - g(mu, v) > 0,
       mu = E[A] mu_w,  v = E[A]^2 sigma_w^2,  delta = Var(A) (mu_w^2 + sigma_w^2),

   where g(mu, v) = E relu(N(mu, v)). So the Gaussian closure overestimates at every alpha, with gap
   ~ (delta/2) phi(alpha)/sigma.
3. Cumulant form: with E A = 1 and Var A = nu, the gain template is kappa_3 = 6 mu sigma^2 nu + O(nu^(3/2)) and
   kappa_4 = 12 sigma^4 nu + O(nu^(3/2)). Its Edgeworth correction, from (2.1), is

       (kappa_3/6) E relu''' + (kappa_4/24) E relu'''' = -(nu/2)(mu^2 + sigma^2) phi(alpha)/sigma < 0,

   negative at every alpha. The kappa_3 and kappa_4 terms separately change sign (at alpha = 0 and |alpha| = 1); their
   sum does not.

*Proof.*
- (1) relu(A w) = A relu(w) since A >= 0.
- (2) By homogeneity E relu(z) = E[A] g(mu_w, sigma_w^2) = g(E[A] mu_w, E[A]^2 sigma_w^2), while
  Var z = E[A]^2 sigma_w^2 + Var(A)(mu_w^2 + sigma_w^2). g is strictly increasing in its variance argument, with
  dg/dv = p(0)/2 > 0.
- (3) Expand E z^k = E[A^k] E w^k to first order in the fluctuation of A. The mu^4 and mu^2 sigma^2 coefficients of
  kappa_4 cancel exactly, leaving 12 sigma^4 nu. Then insert (2.1): t_3 + t_4 = nu phi sigma (-alpha^2 + (alpha^2 - 1)/2).

The input radius is an exact instance: z_l(x) = R z_l(sqrt(n) theta), R = |x|/sqrt(n) independent of theta,
Var R ~ 1/(2n). The network's own norm fluctuations |x_l|^2/n are approximate instances, and by note XXXVI's nullalign
they dominate the cumulants at depth: the truth's kappa_4 is the gain mode, with t3 = t4 from layer 11 on.

*Measured* (`exp1_closure_sign.py`). The local readout error is e_l = E relu(z_l) - sigma G(mu/sigma), computed with the
true mean and variance of z_l, i.e. the Gaussian closure's one-step error at true inputs.

n = 16, L = 8, 2e7 samples (`exp1_n16.txt`); "frac<0" counts active neurons:

| layer | GC chain MSE | local e: rms | frac(e<0) | radial part e_rad: rms | frac(e_rad<0) | Edgeworth R^2 (kappa_3 + kappa_4) | corr(e, e_rad) |
|---|---|---|---|---|---|---|---|
| 1 | 8.3e-9 | 2.0e-5 | 0.44 (noise; z_1 is exactly Gaussian) | 8.0e-3 | 1.00 | | 0.07 |
| 2 | 5.6e-4 | 2.4e-2 | 0.94 | 6.8e-3 | 1.00 | 0.975 | 0.85 |
| 3 | 1.4e-3 | 3.8e-2 | 0.94 | 7.4e-3 | 1.00 | 0.986 | 0.79 |
| 4 | 1.9e-3 | 3.2e-2 | 1.00 | 5.6e-3 | 1.00 | 0.979 | 0.70 |
| 5 | 3.1e-3 | 2.3e-2 | 1.00 | 4.1e-3 | 1.00 | 0.944 | 0.53 |
| 6 | 5.2e-3 | 2.6e-2 | 1.00 | 4.8e-3 | 1.00 | 0.972 | 0.67 |
| 7 | 8.6e-3 | 2.3e-2 | 1.00 | 3.8e-3 | 1.00 | 0.970 | 0.49 |
| 8 | 1.1e-2 | 2.5e-2 | 1.00 | 4.2e-3 | 1.00 | 0.851 | 0.61 |

At layer 1 the radial and spherical parts cancel exactly: z_1 is Gaussian, so the spherical law's negative kurtosis
offsets the radial mixture.

n = 64, L = 16, 4e6 samples (`exp1_n64.txt`):

| layers | local e: rms | frac(e<0) | radial part: rms (frac<0) | Edgeworth R^2 (kappa_3, kappa_4 only / with k5, k6, k3^2, k3k4) |
|---|---|---|---|---|
| 2-4 | 6.1e-3 to 8.2e-3 | 0.94-1.00 | 1.5e-3 to 2.0e-3 (1.00) | 0.995-0.998 / 0.998-1.000 |
| 5-9 | 7.8e-3 to 9.0e-3 | 0.95-0.98 | 1.1e-3 to 1.4e-3 (1.00) | 0.982-0.992 / 0.975-0.998 |
| 10-13 | 8.7e-3 to 1.1e-2 | 0.84-0.97 | 0.9e-3 to 1.2e-3 (1.00) | 0.950-0.977 / 0.963-0.988 |
| 14-16 | 4.3e-3 to 6.9e-3 | 0.89-0.95 | 3.4e-4 to 6.8e-4 (1.00) | 0.932-0.967 / 0.79-0.92 |

At the output, truth < GC on 86% of neurons, and mean error / rms error = -0.53. The sign rule survives width 64 and
depth 16. The radial (input-norm) part is exactly sign-definite but small. The sign comes from the network's own gain
fluctuations, which dominate the cumulants at depth.

Summary of both runs:
- The radial part predicted by Theorem 3.1 is negative on 100% of active neurons at every layer, as proved.
- The total local error has the same sign from layer 2 on: on 94-100% of neurons at n = 16, and on 84-98% at n = 64.
- The kappa_3 + kappa_4 Edgeworth terms explain 85-99% (n = 16) and 93-99.8% (n = 64) of it (R^2).
- So in the gain-dominated regime the Gaussian closure overestimates essentially every post-activation mean.
- Propagated through signed weights, the final-layer error keeps a bias: truth < GC on 81% (n = 16) and 86% (n = 64)
  of neurons, with mean/rms -0.58 and -0.53.

**Prediction for the production chain (testable on existing dumps, section 7, E0).** The chain carries the gain
template, with 98% of the kappa_3 and 83-85% of the kappa_4 gain amplitude (note XXXVI). The missing 15% of the kappa_4
amplitude predicts a residual 0.075 nu sigma (alpha^2 - 1) phi(alpha) per neuron: of order 1e-5, i.e. MSE ~ 1e-10, two
orders below the chain's 1.5e-8. Hence the production residual should have no definite sign:
- the fraction of neurons with truth < chain at each layer >= 7 lies in [0.42, 0.58];
- |corr(residual, (1 + alpha^2) phi(alpha) sigma)| < 0.15.

This agrees with the frontier finding that the residual "has no gain or scale component". The sign rule explains why
the gain must be carried, not that anything is left.

---

## 4. Exact identities from localization

**Theorem 4.1 (kink, or Tanaka, decomposition of the mean).** For a bias-free ReLU network with Gaussian input, if every
kink surface {z_(l,j) = 0} is crossed transversally (no pre-activation is identically zero, or proportional to a single
upstream unit, on an open set),

    E F_i = sum_(l=1..L) sum_j E[ delta(z_(l,j)) |grad z_(l,j)|^2 dF_i/dy_(l,j) ]
          = sum_(l,j) p_(l,j)(0) E[ |grad z_(l,j)|^2 dF_i/dy_(l,j) | z_(l,j) = 0 ].

*Proof.* Euler gives F_i = x . grad F_i a.e., and Gaussian integration by parts gives E[x . grad F_i] = E[Delta F_i]
(distributionally, F_i Lipschitz with BV gradient). Across {z_(l,j) = 0}, grad F_i jumps by
(dF_i/dy_(l,j)) grad z_(l,j), so Delta F_i carries the single-layer density (dF_i/dy_(l,j)) |grad z_(l,j)| dS. The
coarea formula gives dS = |grad z| delta(z) dx. Off the kink set F_i is linear, so Delta F_i = 0 there.

Equivalently, by Tanaka's formula along any localization, E relu(z) = relu(E z) + (1/2) E L^0_inf: the mean is the
expected local time at the kinks. That is the only place convex order touches E F.

*Measured* (`f5/exp4_identities.py` (a), smoothed delta delta_h).
- The transversality hypothesis is not decorative. At n = 3 and n = 6 the identity fails by up to 0.07: cones where all
  upstream units are dead, or only one is active, have probability about (n+1) 2^(-n) and produce spurious half-kinks.
- At n = 24, L = 3 (contamination ~1e-6) it holds once the O(h) smoothing bias is removed. The bias is O(h), not O(h^2),
  because E[dF_i/dy | z = u] has a corner at u = 0.

n = 24, L = 3, 1.5e6 samples:
- smoothed-delta estimates at h = 0.2, 0.1, 0.05 approach E F linearly in h;
- the linear extrapolation 2K(0.05) - K(0.1) matches E F to rms 0.0049 over 24 neurons (rms E F = 0.78, i.e. 0.6%),
  within the kink estimator's own noise;
- at n = 6 the same code misses by up to 0.07, from the non-transversal cones.

**Theorem 4.2 (the Dynkin formula for estimators).** Let Psi(a, Sigma) be any estimator defined on Gaussian inputs
N(a, Sigma): the chain, the GC chain, anything. Assume it is C^2 in a, C^1 in Sigma, continuous at Sigma = 0 with
Psi(a, 0) = F(a) (exact on point masses, as every forward-moment chain is). Define its heat defect

    D_Psi = d_Sigma Psi - (1/2) grad_a^2 Psi        (an n x n matrix; zero for the exact functional Psi*(a,Sigma) = E F(a + Sigma^(1/2) g)).

Then along Gaussian stochastic localization of N(0, I), with m_t the posterior mean and Sigma_t = I/(1+t),

    bias(Psi) := Psi(0, I) - E F = int_0^inf E[ D_Psi(m_t, Sigma_t) : Sigma_t^2 ] dt.                               (4.1)

*Proof.* The innovation form of the posterior mean is dm_t = Sigma_t dW~_t, so d<m>_t = Sigma_t^2 dt and
dSigma_t = -Sigma_t^2 dt. Ito's formula gives dPsi(m_t, Sigma_t) = (martingale) - D_Psi : Sigma_t^2 dt. Take
expectations and let t -> inf: Psi(m_t, Sigma_t) -> Psi(X, 0) = F(X), and E F(X) = E F.

The bias of every estimator of this kind is the localization-integrated failure of the heat equation. The chain is
exact iff its defect vanishes on the localization orbit.

---

## 5. The localization-Richardson midpoint

### 5.1 The localization curve and its s^2 law

Define the localization curve of an estimator Psi:

    beta(s) := E_(m ~ N(0, (1 - s^2) I)) [ Psi(m, s^2 I) ] - E F,   s in (0, 1].

beta(1) = bias(Psi) and beta(0+) = 0. *Measured* for the Gaussian-closure chain (GC), n = 16, L = 8
(`exp3_localize.py` (a); 24000 antithetic leaf pairs per point; truth from 2e7 samples):

| t | s^2 | MSE ratio to t = 0 | rms-bias ratio | exponent p in beta ~ s^p | per-leaf-pair variance |
|---|---|---|---|---|---|
| 0.25 | 0.800 | 0.639 | 0.799 | 2.01 | 2.6e-3 |
| 1 | 0.500 | 0.231 | 0.481 | 2.11 | 1.7e-2 |
| 4 | 0.200 | 0.022 | 0.149 | 2.37 | 5.3e-2 |
| 16 | 0.059 | 6.8e-4 | 0.026 | 2.57 | |

The GC bias is linear in s^2 near s = 1 (p = 2.01) and steepens slowly. Geometric reading:
- localization preserves the input radius on average (|m|^2 + n s^2 ~ n), so a leaf is a spherical cap of angular size
  ~ s;
- the leaf's non-Gaussianity comes from the kinks inside the cap, and its readout effect scales like the cap's area.

**Richardson in s^2 works for the bias.** Linear extrapolation from two points of the curve:

| points (t) | {0, 0.25} | {0.25, 1} | {1, 4} | {0, 0.25, 1} (quadratic) |
|---|---|---|---|---|
| MSE ratio to GC | 0.073 | 0.032 | 0.011 | 0.022 |

### 5.2 Random leaves kill it: the variance floor

**Theorem 5.1 (Richardson variance floor).** Estimate beta's curve point by K antithetic leaf pairs,
E^(s) = (1/K) sum_k (Psi(m_k, s^2 I) + Psi(-m_k, s^2 I))/2. Suppose Psi is accurate per leaf, in the sense that the
leaf error is small against leaf fluctuations. Then

    Var E^(s) >= (1/K) sum_(q even >= 2) (1 - s^2)^q V_q,

by (1.1) and antithetic symmetry, which keeps exactly the even chaos. The Richardson combination
R = (E^(s) - s^p Psi(0, I)) / (1 - s^p), unbiased under beta ~ s^p, therefore has

    Var R >= [(1 - s^2)/(1 - s^p)]^2 V_2/K >= min(1, 4/p^2) V_2/K     for every s.

At p = 2 the floor is V_2/K exactly, independent of the localization depth.

The depth-16 chaos profile at infinite width (`chaos_profile.py`) comes from E[F(x)F(x')] = rho^(o16)(r), with
rho(c) = (sqrt(1-c^2) + (pi - arccos c) c)/pi, and Taylor coefficients by a Cauchy integral:

| depth | V_0 | sum V_q (q >= 1) | V_1 | V_2 | V_3 | V_4 | share even q >= 2 | share odd q >= 3 | share q > 20 |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.318 | 0.682 | 0.500 | 0.159 | 0 | 0.013 | 0.27 | 0 | 0.001 |
| 8 | 0.834 | 0.165 | 0.050 | 0.036 | 0.020 | 0.014 | 0.42 | 0.27 | 0.04 |
| 16 | 0.930 | 0.0705 | 0.0133 | 0.0108 | 0.0073 | 0.0056 | 0.45 | 0.36 | 0.12 |

The total 0.0705 matches the measured sigma^2 ~ 0.074. With V_2 = 0.0108, raw 1e-8 needs K >= 1.1e6 leaf pairs, each a
chain run: about 10^6 times the budget. Even the bare curve point E^(s) at t = 0.01 already needs about 100 pairs for
outer noise 1e-8, while the bias moved by 1%. Random-leaf localization at Phase-2 precision is dead.

### 5.3 Deterministic leaves along a few directions

Localize exactly along a k-dimensional subspace U (leaf N(U c, I - U U^T)) with a tensor Gauss-Hermite rule in c, 12
nodes per direction. The outer integral is then exact to quadrature error. `exp3_localize.py` (b), GC, n = 16, L = 8;
response-directed U = top eigenvectors of sum_i L_i L_i^T with L_i = E[F_i x]:

| k (k/n) | response-directed MSE ratio | random U (6 draws) |
|---|---|---|
| 1 (0.06) | 0.731 | 0.888 +- 0.125 |
| 2 (0.12) | 0.469 | 0.739 +- 0.060 |
| 3 (0.19) | 0.280 | 0.739 +- 0.081 |

Per direction this removes about 1.8 k/n of the MSE at random and about 4 k/n response-directed. At n = 1024 that is
0.2-0.4% per direction, for 5-12 chain runs per direction. Dead on cost.

### 5.4 The derivative form: a deterministic Richardson midpoint

The random-leaf obstruction disappears if the curve is extrapolated from its derivative at s = 1, which is
deterministic.

**Theorem 5.2 (localization-Richardson from the heat defect).** Let Psi be positively homogeneous of degree 1
(Psi(c a, c^2 Sigma) = c Psi(a, Sigma)) and satisfy the hypotheses of Theorem 4.2. Then:
1. d beta / d(s^2) at s = 1 equals tr D_Psi(0, I) = Psi(0, I)/2 - Delta_a Psi(0, I)/2.
2. If beta(s) = b s^p (the Richardson ansatz), then

       E F = Psi(0, I) - (2/p) tr D_Psi(0, I) = (1 - 1/p) Psi(0, I) + (1/p) Delta_a Psi(0, I).                    (5.1)

3. At p = 2 this is the midpoint E F = (Psi + Delta_a Psi)/2, and the errors of Psi and of its Euler-Stein companion
   Delta_a Psi are exactly opposite: err(Delta_a Psi) = -err(Psi). The two estimators bracket the truth, and (5.1) is
   the midpoint of the sandwich.

*Proof.*
- (1) By (4.1), d beta/dt at t = 0 is -tr D_Psi(0, I), and ds^2/dt = -1 there. Homogeneity in Sigma at a = 0 gives
  Psi(0, lambda I) = sqrt(lambda) Psi(0, I), hence tr d_Sigma Psi(0, I) = Psi(0, I)/2.
- (2) beta = b s^p gives d beta/d(s^2) at s = 1 equal to (p/2) beta(1), so beta(1) = (2/p) tr D.
- (3) The exact functional satisfies Delta_a Psi* (0, I) = E Delta F = E F (Euler-Stein), so
  err(Delta_a Psi) = Delta_a eps, where eps = Psi - Psi*. The ansatz at p = 2 is equivalent to Delta_a eps = -eps.

**What the midpoint computes (Theorem 5.3, the defect calculus).** D is a second-order operator, so for a composite
h(S(a, Sigma)) of the chain's state,

    D[h o S] = sum_i h_i D[S_i] - (1/2) sum_ij h_ij <grad_a S_i, grad_a S_j>_Sigma.                               (5.2)

Take one Gaussian-closure readout g(mu, v) of a pre-activation whose mean mu and second moment both satisfy the heat
equation (exact Gaussian expectations of the input), so that D[mu] = 0 and D[v] = grad mu grad mu^T. Using g_mumu = 2 g_v,

    -tr D[g] = g_(mu v) <grad mu, grad v>_Sigma + (1/2) g_(vv) |grad v|^2_Sigma
             = (kappa3_loc/6) E relu''' + (kappa4_loc/24) E relu'''',
    kappa3_loc = 3 <grad mu, grad v>_Sigma,   kappa4_loc = 3 |grad v|^2_Sigma.

These are precisely the third and fourth cumulants of the infinitesimal-leaf mixture: the leaves' means and variances
co-vary (skewness) and their variances vary (kurtosis). For the generalized chi-square (second-chaos) model of note XLI
they are path terms: grad mu = L, grad v = 2 H L. So the drift of the GC chain contains, with no fitted parameter, the
two second-chaos corrections the production chain builds by hand: D3's path part and V56's P4 path class. The midpoint
weights both by 1; the correct, class-resolved weights are 1/2 and 1 (Proposition 5.4). Deeper in the chain, (5.2)
composes these corrections through every later layer automatically.

*Proof.* Direct from (5.2) with h = g. The cumulant identities: a Gaussian mixture with component mean mu + grad mu . m
and variance v + grad v . m, m ~ N(0, dt Sigma), has kappa_3 = 3 dt <grad mu, grad v> and kappa_4 = 3 dt |grad v|^2 to
first order in dt. Then insert (2.1).

**Proposition 5.4 (what first-order localization sees, and the class-dependent exponent).** Take the second-chaos
model of note XLI, z = mu_0 + L . x + (1/2)(x^T H x - tr H), read by the Gaussian closure. At (0, I):
grad mu = L and grad v = 2 H L. So

    -tr D = 2 g_(mu v) L^T H L + 2 g_(vv) |H L|^2   <->   Edgeworth with kappa3_loc = 6 L^T H L,  kappa4_loc = 12 |H L|^2,

while the true cumulants are kappa_3 = 3 L^T H L + tr H^3 and kappa_4 = 12 |H L|^2 + 3 tr H^4.
1. First-order localization sees only open walks (L^T H^k L, k <= 2). Closed walks tr H^k are invisible to beta'(1).
   They enter only through higher derivatives of the leaf statistics (grad^2 v = 2 H^2, ...), at higher orders of the
   tower (sketch).
2. The Richardson exponent is class-dependent:
   - for the kappa_4 path class (V56's P4), slope = bias, so p = 2 and the midpoint is exact;
   - for the kappa_3 path class (D3's path part), slope = 2 x bias, so p = 4 and the midpoint over-corrects twofold.
3. The class-resolved estimator C - (1/2)[g_(mu v) part of tr D] - [g_(vv) part of tr D] is path-exact.
4. Under the ansatz that each class keeps its exponent on all of u in [0, 1], a mixture of a p = 2 and a p = 4 class
   makes beta a quadratic in u = s^2, and the second-order tower below is exact for it without knowing the mixture.

*Proof.* Theorem 5.3 with grad_a of the leaf mean and variance at a = 0, then (2.1) with E relu''' = 2 g_(mu v),
E relu'''' = 4 g_(vv) (from g_(mu mu) = 2 g_v).

**The Euler-Stein tower.** The exact functional satisfies E[Delta^j F] = (3 - 2j) E[Delta^(j-1) F]. This holds because
Delta^(j-1) F is homogeneous of degree 3 - 2j and E[Delta G] = E[x . grad G] = deg(G) E G. So
Delta_a^j Psi* (0, I) = (1, -1, 3, -15, ...) x E F: every Delta_a^j Psi is a consistent estimator of +-E F.
Expanding the localization curve at u = 1,

    beta(u) + E F = sum_k (1-u)^k u^(1/2-k) Delta_a^k Psi(0, I) / (2^k k!),
    beta'(1) = Psi/2 - Delta Psi/2,   beta''(1) = -Psi/4 + Delta Psi/2 + Delta^2 Psi/4.

Richardson with beta a polynomial of degree k in u (beta(0) = 0) gives deterministic estimators:

    T_1 = Psi - beta'(1) = (Psi + Delta_a Psi)/2,
    T_2 = Psi - beta'(1) + beta''(1)/2 = (3 Psi + 6 Delta_a Psi + Delta_a^2 Psi)/8.

T_2 is also the level-2 iterate: the midpoint's homogeneous extension Psi_1(a, Sigma) = (Psi + a . grad Psi +
tr(Sigma grad^2 Psi))/2 corrected at p = 4. With exact inputs (T, T, -T), both return T. An Aitken variant estimates p
from the pooled ratio beta''/beta' = p/2 - 1.

### 5.5 Measured

**(a) The midpoint across widths, depths and seeds.** GC chain; Delta_a C by central differences (2n+1 runs, h = 0.02,
float64); truth from 2e6-1e7 samples (truth noise <= 1e-7, below every MSE shown). Sources: `f5/exp5_midpoint.py`,
`f5/loo_kappa.py`.

| network | GC MSE | T_1 = (C + Delta C)/2: MSE ratio | oracle kappa | LOO global kappa: ratio | corr(T - C, -tr D) | corr(err_C, err_Delta) |
|---|---|---|---|---|---|---|
| n=16 L=8 s1 | 1.08e-2 | 0.097 | -0.93 | 0.092 | 0.93 | -0.76 |
| n=16 L=8 s2 | 3.9e-4 | 0.834 | -0.57 | 0.778 | 0.82 | -0.56 |
| n=16 L=8 s3 | 2.05e-2 | 0.153 | -0.73 | 0.097 | 0.99 | -0.99 |
| n=32 L=8 s1 | 1.08e-3 | 0.201 | -0.91 | 0.193 | 0.89 | -0.69 |
| n=32 L=8 s2 | 2.10e-3 | 0.038 | -0.96 | 0.039 | 0.96 | -0.87 |
| n=32 L=16 s1 | 1.61e-4 | 0.357 | -1.25 | 0.376 | 0.74 | -0.05 |
| n=64 L=8 s1 | 5.45e-4 | 0.091 | -0.83 | 0.070 | 0.97 | -0.92 |
| n=64 L=16 s3 | 1.56e-3 | 0.037 | -0.93 | 0.032 | 0.98 | -0.92 |
| n=64 L=16 s4 | 1.33e-4 | 0.487 | -1.30 | 0.505 | 0.75 | +0.16 |
| n=64 L=16 s5 | 5.48e-4 | 0.041 | -0.93 | 0.035 | 0.98 | -0.92 |
| **n=128 L=16 s1** | 1.24e-4 | **0.090** | -0.95 | 0.088 | 0.95 | -0.84 |
| **n=128 L=16 s2** | 3.16e-4 | **0.027** | -0.97 | 0.028 | 0.98 | -0.93 |

Across the 12 networks:
- The geometric-mean MSE ratio is 0.115 for T_1 (no parameter) and 0.106 with a leave-one-out global kappa. The LOO
  kappa is stable at -0.92 to -0.94.
- The two weak cases (n16 s2, n64 s4) have small GC errors and fitted kappa away from -1. Their residual is dominated by
  classes that are not p = 2.
- **Width trend.** At n = 128 the correlation is 0.95-0.98 and the oracle kappa is -0.95 and -0.97, converging to the
  derived -1. The midpoint removes 91-97% of the closure's MSE.
- **Absolute residual after the midpoint**, depth 16: median GC MSE goes from 5.5e-4 (n = 64) to ~2e-4 (n = 128), a
  factor ~2.8 per doubling. T_1's MSE goes from 5.8e-5 to ~1e-5, a factor ~6 per doubling. Three more doublings
  extrapolate, with an uncertainty of about one decade, to GC ~1e-5 and T_1 ~5e-8 at n = 1024, within ~3x of the
  production chain. This is only an extrapolation; E1a measures it.

**(b) The tower and the exponent** (`f5/exp5c_tower.py`; the bi-Laplacian by 4th-order differences, stable over
h = 0.04-0.1):

| network | T_1 | T_2 = (3C + 6 Delta C + Delta^2 C)/8 | Aitken (p_hat) |
|---|---|---|---|
| n=16 s1 | 0.097 | 0.032 | 0.133 (1.77) |
| n=16 s2 | 0.834 | 0.728 | 0.672 (5.46) |
| n=16 s3 | 0.153 | 0.042 | 0.019 (2.54) |
| n=32 s1 | 0.201 | 0.121 | 0.193 (2.26) |
| n=32 s2 | 0.038 | 0.120 | 0.948 (1.05) |

- The level-2 iterate equals T_2. On n = 16 s1 its own drift predicts the midpoint's error with correlation 0.84, and
  the oracle level-2 factor is -0.44 (p ~ 4.5).
- T_2 helps where the error mixes p = 2 and p = 4 classes and hurts where it does not (n32 s2). Aitken is unstable.
- **T_1, or T_kappa with a global kappa, is the robust member of the family.**

**(c) Which tangents carry the drift** (`f5/exp6_truncated.py`). Laplacians of truncated variants of the GC chain, all
equal to C at a = 0: V1 freezes the pre-activation correlations, so only means and variances respond (local tangents);
V2 freezes the whole covariance, so only means respond.

| network | full T_1 | V1: T_1 | V1: oracle kappa -> ratio | V2: oracle kappa -> ratio | corr(Delta C_V1, Delta C) |
|---|---|---|---|---|---|
| n=16 s1 | 0.097 | 0.491 | -1.61 -> 0.405 | -0.51 -> 0.111 | 0.992 |
| n=32 s2 | 0.038 | 0.249 | -1.59 -> 0.129 | -0.31 -> 0.045 | 0.997 |
| n=64 L16 s3 | 0.037 | 0.305 | -1.73 -> 0.154 | -0.22 -> 0.145 | 1.000 |
| n=64 L16 s5 | 0.041 | 0.345 | -1.68 -> 0.218 | -0.17 -> 0.279 | 0.999 |

- The local tangents get the drift's shape, but they miss about 40% of its magnitude (kappa ~ -1.65 is needed, stable
  across networks) and half of its value.
- The rest sits in the correlation tangents, the covariance-tangent Grams <dC_ab, dC_cd>. Section 6.2 shows these are
  the closed (K4) contractions. **The part of the drift that is cheap is the weak part.**

---

## 6. The estimator, its cost, and its prediction

### 6.1 The estimator family

For any Gaussian-input estimator Psi that is homogeneous and exact on point masses (the production chain qualifies):
- **T_1 = (Psi + Delta_a Psi)/2**, the localization-Richardson midpoint;
- **T_kappa = Psi + kappa tr D_Psi**, with one global kappa fitted like a counterterm (T_1 is kappa = -1);
- **T_cr**, class-resolved: weight 1/2 on the g_(mu v) terms of the drift and 1 on the g_(vv) terms (Proposition 5.4).
  This needs the drift decomposed by (5.2), not just its trace;
- **T_2 = (3 Psi + 6 Delta_a Psi + Delta_a^2 Psi)/8**, the second-order tower, with no free parameter.

Here Delta_a Psi = tr(Sigma_1 grad^2_(mu_1) Psi) is the Laplacian of the estimator in the input mean, i.e. in the first
pre-activation mean mu_1 = a W_0 with metric Sigma_1 = W_0^T W_0.

### 6.2 Cost at n = 1024, in units of 2n^3

| route | T_1 | T_2 | note |
|---|---|---|---|
| finite differences, production chain (~200 units/run) | (2n+1) x 200 ~ 4.1e5 units = 400 B | 2n^2 runs ~ 4e8 units | offline only; float64 required (float32 roundoff summed over 1024 directions ~3e-3/neuron at h = 0.05) |
| finite differences, GC chain (~32 units/run: the covariance sandwich, 2 units/layer) | ~6.6e4 units = 64 B | ~6.6e7 units | offline diagnostic (E1a) |
| structural forward drift (5.2), GC chain, local tangents only (V1) | ~+145 units (~9 units/layer: tangents 5, Grams 3, variance drift 1) | | measured: with kappa ~ -1.65 it removes 60-87% of the GC MSE, vs 90-97% for the exact drift (5.5c) |
| structural drift, GC chain, exact | n^4 per layer = 512 units/layer | | the covariance-tangent Grams <dC_ab, dC_cd> are closed (K4) contractions |
| structural drift, production chain | >= n^4 per source-layer | | the tangents of the A, P, M legs are n x n x n per source |
| Hutchinson | variance >= 2 \|grad^2 Psi\|_F^2 ~ 4 V_2 ~ 0.04 per probe: 4e6 probes | | with the chain's own H as a control variate (10% Frobenius residual): >= 1e4 probes x ~3 chains |
| dressed drift (Schur-hub Grams through the current field, as V56) | conjecturally ~+2.5 units/layer ~ +40 units | | speculative; section 6.4 |

The defect calculus does not dodge note XL's cost classification; it reproduces it.
- For the GC chain, the exact drift contains the V2-wall objects (covariance-tangent Grams).
- Its open-walk part (mean and variance tangents) costs n^3 per layer.
- Its closed-walk part costs n^4 per layer.

### 6.3 Predictions

**(P-F5-1) Gaussian closure at n = 1024 (E1a).** Extrapolated from the width trend (5.5a), per network:
- corr(T - C, -tr D) >= 0.95 and the oracle kappa in [-1.05, -0.90];
- T_1 removes >= 90% of the GC chain's MSE;
- T_kappa with kappa from network 0 transfers to network 1 within 5 points;
- in absolute terms, GC raw ~1e-5 (3e-6 to 3e-5) and T_1 raw ~5e-8 (1e-8 to 3e-7);
- the local-tangent drift (V1) needs kappa ~ -1.6 and removes 60-85%.

Why: the GC error is dominated more and more by the p = 2 classes (gain mode, kappa_4 path) as the width grows, and
those are exactly what first-order localization sees.

**(P-F5-2) Production chain at n = 1024, free-running (E1b).**
- corr(T - C_prod, -tr D_prod) in [0.25, 0.65] and kappa_hat in [-1.0, -0.4].
- T_kappa with kappa fitted on net 0 lowers net 1's raw MSE by 10-35%.
- Why (heuristic): the production residual is about 10% incoherent error in bulk kappa_3 transport plus joint-gate
  (Mehler pair) effects. Both should be visible to first-order localization, since averaging product gates over
  infinitesimal leaves produces the joint-gate correlation. The closed-walk part (triangle, 4-cycle) is drift-invisible
  and should hold at least a third of the residual.

**(P-F5-3) Sign audit of the production residual (E0).**
- At every layer >= 7, the fraction of neurons with truth < chain lies in [0.42, 0.58].
- |corr(residual, (1 + alpha^2) phi(alpha) sigma)| < 0.15.

### 6.4 What F5 says about system design

1. **Exactness criterion.** A moment chain is exact iff it is a martingale along stochastic localization of its input:
   it satisfies the heat equation on the posterior orbit. Its bias is its total drift (4.1). This is checkable without
   ground truth, network by network.
2. **The drift is a compiler.** By (5.2) the drift decomposes into local sources at each layer, propagated by the
   chain's own linear response.
   - First-order localization produces only open walks, with their derived weights. The class-resolved factors are 1
     for the kappa_4 path class and 1/2 for the kappa_3 path class.
   - Closed walks need the second-order tower.
   - So the drift, computed once offline on training networks, says which omitted classes carry the residual and at
     what weight, with no regression. For the kappa_4 path class it reproduces V56's derived weight 1. Whether the
     drift also accounts for V56's third-chaos star, at its weight 1, is not established here: third-chaos content
     enters (5.2) only through compositions.
3. **In-budget use is conditional.** The exact drift of the production chain is out of budget at any order. A dressed
   drift (Grams through the current field's covariance via the Schur identity, as V56's hub) is the only candidate for
   an in-budget form. Its accuracy is the same open question V56 answered for P4. Build it only if E1b shows a
   drift-visible residual.
4. **What F5 closes.** Certified sandwiches, localization mixtures with random leaves, Berezin-Lieb lower symbols,
   partition-function chains and subspace localization are out at Phase-2 precision. Each is quantified above by a
   factor of 10^2 or more.


---

## 7. The decisive experiments (AWS)

**E0 (free; existing dumps).** Inputs: the production chain's per-layer post-activation means (free-running) and
truth 'm' per layer on networks 0-15.
- Outputs, per layer >= 7: frac(truth < chain), corr(residual, (1 + alpha^2) phi(alpha) sigma), and corr(residual,
  (alpha^2 - 1) phi(alpha) sigma). Here (alpha, sigma) come from the chain's pre-activation state.
- Decision: P-F5-3 holds -> no sign-definite (gain-mode) component is left. If it fails, one gain-template counterterm
  (amplitude of the kappa_4 template 12 sigma^4 nu) is the fix, at zero cost.

**E1a (Gaussian-closure drift at n = 1024; ~6 core-hours per network).**
- Inputs: W_off{0,1}, truth_off{0,1}['m'].
- Script outline: a standalone float64 numpy GC chain. Port `loc/f5lib.py:gc_chain`, which needs the bivariate ReLU
  moments via Owen's T: ~2e6 owens_t calls and two 1024^3 matmuls per layer, ~10 s per run. Run it at
  a in {0, +-h e_k : k = 1..1024}, h = 0.02 (a enters as mu_1 = a W_0; Sigma_1 = W_0^T W_0 unchanged).
- Outputs: C, Delta_a C = sum_k (C(h e_k) + C(-h e_k) - 2C(0))/h^2, T_1, kappa_hat, corr(T - C, -tr D), MSE ratios.
  Add an h-check on 64 directions at h in {0.01, 0.02, 0.04}.
- Also run the frozen-correlation variant V1 (`f5/exp6_truncated.py`). Report corr(err_(T_1-GC), err_production) and the
  inverse-MSE-weighted blend of T_1-GC with the production chain: two parameter-light estimators with different error
  physics.
- Decision:
  - corr >= 0.7 and T_1 ratio <= 0.5 on both networks -> the drift principle scales to full width; run E1b.
  - corr < 0.4 -> it does not scale; close F5 apart from the sign rule.

**E1b (production chain drift at n = 1024; ~10-20 core-hours per network).**
- Patch `est_v29.py`:
  - line 1046: mu = a, read from an .npy named by an env var, instead of zeros;
  - run in float64 (rebind the module's f32 to float64);
  - freeze every data-dependent discrete decision at its a = 0 value: saturation masks, rank truncations, the old-tier
    bases, the age gate. This keeps the map smooth in a;
  - counterterms off (free-running).
- Run the 2n+1 = 2049 inputs above on networks 0 and 1. Truth as in E1a.
- Outputs: tr D_prod per neuron; corr(T - C_prod, -tr D_prod); kappa_hat on each network; the cross-applied
  T_kappa (kappa from net 0 on net 1, and back); T_1.
- Decision:
  - corr >= 0.4 and the cross-applied T_kappa lowers raw by >= 20% -> the residual is drift-visible. Then decompose the
    drift by layer and class with (5.2), which needs the 2n+1 runs' per-layer states, and cost its dominant open-walk
    terms in a dressed form. Adopt only if one term delivers >= 10% raw for <= +20 units.
  - corr < 0.25 -> the production residual is drift-invisible (closed walks, the second tower). F5 closes, and the
    remaining error is the n^4 class note XLI priced.
- Compute: 2 x 2049 runs x ~10-20 s (float64, single core) ~ 12-24 core-hours. Minutes on the fleet.


---

## 8. Claims ledger

**Novelty.** The classical ingredients are used as such:
- Strassen's theorem, Koshevoy-Mosler lift zonoids, Markov-Krein / Karlin-Studden extremal moment problems;
- Lieb-Berezin inequalities, Funk-Hecke, Tanaka's formula, Ito/Dynkin;
- Price's theorem, the Mehler semigroup, Glauber-Sudarshan nonclassical depth.

The nearest relative of the localization-Richardson construction is zero-noise extrapolation in quantum error
mitigation, which is Richardson in a noise strength from amplified-noise runs. Here the "noise" is the leaf variance,
and its derivatives come from the estimator itself, so no extra noisy runs are needed. A targeted alphaXiv search
(9 October 2026) found no prior statement of the following:
- the estimator-drift identity (4.1) for moment-closure estimators;
- the Euler-Stein tower estimators;
- the class-dependent localization exponents;
- the spherical form of Theorem D.

Priority is not claimed.

| claim | status | checked |
|---|---|---|
| D1: localization = Mehler semigroup; Var M_t = sum (1-s^2)^q V_q | proved | chaos profile vs measured sigma^2 (0.0705 vs 0.074) |
| 2.1: mean invisible to martingale couplings and to any two-sided MGF sandwich; no Ising <= boson <= Ising bound on E F | proved | (elementary) |
| 2.2 / 2.3: one-layer sandwich = zonoid inclusion; no propagation through signed layers (H_i indefinite) | proved | H_i 52% negative eigenvalues; dilation lowers E F by 0.019 (pred. 0.021) |
| 2.4: matched closures are cx-incomparable | proved | |
| 2.5: Markov-Krein sandwich widths 50-1500x target; certified Gaussian fraction ~0.78 needed with 6 moments | proved (LP) + measured | exp2 |
| 2.6: spherical Berezin-Lieb (Theorem D transplanted), beta_2(1) = 2/(n+2), beta_odd = 0, spectral variance = sum beta_q \|\|f_q\|\|^2 | proved | exact quadrature; MC beta_2(1) = 0.2006 vs 0.2000 |
| 3.1: gain-mode sign rule (relu commutes with the gain; GC overestimates; cumulant template sign-definite) | proved | GC local error < 0 on 94-100% of neurons (n=16); radial part 100% |
| 4.1: kink (Tanaka) decomposition, under transversality | proved | n=24: 0.6% relative after O(h) extrapolation; fails at n <= 6 as predicted by non-transversality |
| 4.2: bias = localization-integrated heat defect (estimator drift) | proved (Ito/Dynkin), assuming Psi is C^{2,1}, the martingale term integrable, and E Psi(m_t, Sigma_t) -> E F (dominated convergence) | localization curve E(t) -> T (exp3a) |
| 5.1: random-leaf Richardson variance floor >= min(1, 4/p^2) V_2/K | proved (given per-leaf accuracy) | V_2 = 0.0108 at depth 16 |
| 5.2: derivative Richardson T = Psi - (2/p) tr D; midpoint at p = 2; sandwich err(Delta Psi) = -err(Psi) | proved (under the ansatz) | corr(err_C, err_Delta) between -0.56 and -0.99 on 10 of 12 nets (-0.05 and +0.16 on the two weak ones) |
| 5.3: drift of one GC readout = Edgeworth of the infinitesimal-leaf mixture (kappa3 = 3<grad mu, grad v>, kappa4 = 3\|grad v\|^2) | proved | |
| 5.4: first-order localization sees open walks only; p = 2 for the kappa_4 path, p = 4 for the kappa_3 path; T_2 exact for quadratic curves | proved for the second-chaos model; sketch in general | oracle kappa varies -0.57 to -0.96 across nets, as a class mixture predicts |
| Euler-Stein tower E[Delta^j F] = (3-2j) E[Delta^(j-1) F] | proved | |
| s^2 law of the GC localization curve near s = 1 | measured (n = 16: p = 2.01) | conjecture in general |
| midpoint / tower gains on GC | measured, small n (table 5.5) | width trend: E1a |
| production residual drift-visible | conjecture (P-F5-2) | E1b |
| dressed (hub) in-budget drift | conjecture | after E1b |


---

## Appendix: scripts and outputs (all under `scratchpad/loc/`)

| file | content |
|---|---|
| `f5lib.py` | He nets, Monte Carlo moments, the exact Gaussian-closure chain (bivariate ReLU moments via Owen's T, checked against MC) |
| `exp1_closure_sign.py`, `exp1_n16.txt`, `exp1_n64.txt` | local readout errors, signs, Edgeworth and radial decomposition |
| `exp2_mk.py`, `exp2_out.txt` | Markov-Krein LP sandwiches, certified Gaussian component |
| `exp3_localize.py`, `exp3_out.txt` | localization curve, Richardson, subspace localization, Euler-Stein companion (n = 16) |
| `chaos_profile.py`, `chaos_profile.npy` | V_q of the depth-L output at infinite width |
| `f5/exp4_identities.py`, `f5/exp4_out.txt`, `f5/exp4ac_out.txt` | kink decomposition, H_i indefiniteness and dilation counterexample, spherical Berezin-Lieb |
| `f5/berezin_sphere.py` | exact spherical Berezin multipliers beta_q(l) |
| `f5/exp5_midpoint.py`, `f5/exp5_out.txt`, `f5/exp5d_out.txt`, `f5/exp5_n*_L*_s*.npz` | the midpoint across widths, depths and seeds. The "level-2" lines in exp5_out.txt used a wrong off-origin extension and are superseded by exp5b/exp5c |
| `f5/exp5b_level2.py`, `f5/exp5b_out.txt` | level-2 Richardson with the correct homogeneous extension (h-stability of the bi-Laplacian) |
| `f5/exp5c_tower.py`, `f5/exp5c_out.txt` | the Euler-Stein tower T_1, T_2 and the Aitken exponent |
| `f5/loo_kappa.py` | leave-one-out global kappa over all exp5 networks |
| `f5/exp6_truncated.py`, `f5/exp6_out.txt` | drift carried by local tangents (V1) vs mean tangents (V2) vs full |


---

## Referee report

Adversarial referee for F5. Every check below can be rerun from the scripts in `loc/ref_f5/` (`chk1.py`, `chk_gain.py`,
`chk_mid.py`, `chk_coh.py`, `t_owens.py`). Their outputs are quoted.

### R0. Verdict

- **The negative results are correct.** They are mostly elementary, and they should stay closed: Theorems 2.1, 2.3, 2.4,
  2.5, 2.6 and 5.1.
- **I checked the positive algebra and it is right:**
  - the Dynkin identity (4.1);
  - dbeta/du = tr D, and tr D = Psi/2 - Delta Psi/2 under homogeneity;
  - the Euler-Stein identity;
  - the tower coefficients beta'(1) and beta''(1), and T_2;
  - the defect of one GC readout, kappa3_loc = 3<grad mu, grad v> and kappa4_loc = 3|grad v|^2;
  - the second-chaos cumulants of Prop 5.4;
  - the gain template of Theorem 3.1(3).
- **The midpoint numbers reproduce** with independent Monte Carlo truth and are stable in h.
- **The frame does not yield a Phase-2 estimator.** It has no in-budget component, and nothing in it changes final-layer MSE at
  the 1e-8 level inside B. The author says this too.
- **The extrapolation behind P-F5-1 is too optimistic.** Its stated mechanism, "the gain mode is a p = 2 class", is false in
  the clean model (R1).
- **The E1a decision rule is not tied to Phase-2 relevance** (R3).
- **Verdict: does not survive as a Phase-2 system design. Survives as correct theory plus an offline diagnostic.**

### R1. The gain mode is not a p = 2 class: its kappa_4 half is invisible to first-order localization

The prediction (section 6.3) rests on "the GC error is dominated more and more by the p = 2 classes (gain mode, kappa_4
path)". No derivation is given for the gain. I tested the exact independent-gain model of Theorem 3.1:
- x = (y in R^k, g); z = A(y)(mu + sigma g) with A = |y|/sqrt k and k = 50, so nu ~ 0.01;
- Psi(a, sI) = g(mu_z, v_z) with the exact Gaussian moments of z (noncentral chi);
- -tr D = Psi_s - Delta_a Psi/2, computed by finite differences.

This uses dbeta/du = tr D, which does not need homogeneity. (`ref_f5/chk_gain.py`)

```
mu=+0.0: truth-GC -1.990e-03 | template k3 part -0.000e+00 k4 part -1.985e-03 | -trD -1.790e-08 | T1-truth +1.990e-03
mu=+0.5: truth-GC -2.195e-03 | template k3 part -8.757e-04 k4 part -1.314e-03 | -trD -8.684e-04 | T1-truth +1.327e-03
mu=+1.0: truth-GC -2.419e-03 | template k3 part -2.407e-03 k4 part +0.000e+00 | -trD -2.372e-03 | T1-truth +4.757e-05
mu=-0.7: truth-GC -2.322e-03 | template k3 part -1.522e-03 k4 part -7.922e-04 | -trD -1.506e-03 | T1-truth +8.162e-04
```

**What the numbers show.** The drift reproduces the kappa_3 half of the gain template exactly. It misses the kappa_4 half
(12 sigma^4 nu) entirely.

**Analytic reason.** At (0, I):
- grad v = 2 nu sigma mu e_g, so |grad v|^2 = O(nu^2);
- the kappa_4 template is O(nu);
- it is a trace (closed) object, Var A, and not a path through L.

**Consequences.**
- At alpha = 0 the midpoint removes 0% of the gain error.
- The gain is therefore a mixed class: kappa_3 half at p = 2, kappa_4 half invisible at first order.
- Note XII (arrow-mixture) measured that the gain is about 62% of the GC MSE at n = 1024, where GAC removes 62%. So the
  n = 1024 GC error contains a large component that the single-readout theory says T_1 cannot see.
- At n <= 128 the full chain's covariance-tangent terms evidently recover much of it. Whether they still do at n = 1024
  is exactly what E1a measures. The width-trend prediction has no mechanism behind it.

**A related inconsistency.** Prop 5.4 says "closed walks are invisible at first order". Section 5.5c/6.2 says the exact GC
drift's missing 40% sits in "closed (K4)" covariance-tangent Grams. Both cannot be the same notion of closed.
- Prop 5.4 is a statement about one readout with exact input moments.
- The chain's composed drift does contain closed-type contractions.

The claim "first-order localization sees exactly the open walks" should be restricted to the single-readout second-chaos
model.

### R2. The n = 1024 baseline is already known, and the arithmetic does not reach production

**The baseline.** The Gaussian closure on official network 0 at n = 1024, depth 16, is measured: 4.06e-6 (notes XII and
ncg-probability, 1e9-sample truth). The note predicts "GC ~1e-5 (3e-6 to 3e-5)" from 2.8x per width doubling. The repo's
own numbers give about 3.9x per doubling from n = 256 (6.0e-5) to n = 1024.

**The cut T_1 needs.**
- To match production (1.55e-8), T_1 must remove 99.6% of the GC MSE, a 265x cut.
- The largest cut measured is 37x (n = 128 s2), and the geometric mean at n = 128 is about 20x.
- If the cut stays where it was measured, T_1 lands at 1.1e-7 to 3.7e-7: 7-24x worse than production.
- "T_1 ~ 5e-8" needs the cut to keep growing about 2x per doubling for 3 more doublings. That rate comes from a single
  doubling (64 -> 128) on 3 and 2 seeds, whose GC MSEs scatter by 10x.

**The proposed blend is worth almost nothing.** It mixes T_1-GC with production by inverse MSE. Even with uncorrelated
errors, a 4e-7 estimator improves a 1.5e-8 one by at most about 4%, and T_1-GC costs about 64 B.

### R3. The E1a decision rule does not test what matters

"T_1 ratio <= 0.5 and corr >= 0.7 -> run E1b" is a pass at a raw of about 2e-6, two decades from relevance. The
corrected rules:
- **(a) Relevance.** Report T_1 raw against 1.55e-8. Declare the drift a candidate production ingredient only if T_1-GC
  reaches <= 5e-8 on both networks.
- **(b) Mechanism (the R1 test).** Regress T_1's residual on the gain kappa_4 shape (nu sigma/2)(alpha^2 - 1) phi(alpha).
  A significant coefficient near 1 means the gain's kappa_4 half is invisible at full width.
- **(c) Coherent part.** Report the ratio on de-meaned errors, plus an oracle 2-parameter affine fit (T ~ c0 + c1 C) as a
  null. On the existing networks (`chk_coh.py`):
  - the common shift is 5-52% of GC MSE;
  - the affine oracle reaches ratio 0.02-0.94;
  - T_1 on de-meaned errors gives 0.04-0.46.

  So T_1's per-neuron gain is genuine, but part of the headline is the trivial common shift. The production chain has no
  common shift: |mean error| < 1e-5 per layer (DATA_FINDINGS).

**A numerical trap in E1a.** At a = 0 every layer-1 pre-activation mean is exactly 0, so `f5lib.Phi2` substitutes
|h| = 1e-10. In the Laplacian, the base-point error is multiplied by 2n/h^2 = 5.1e6.

Measured (`chk1.py`): rms dC(0) between the 1e-10 and 1e-12 substitutions is 4.4e-12 at n = 64 and 1.6e-12 at n = 128.
The effect falls with n, but the safe extrapolation to n = 1024 is a 1e-5 to 2e-5 rms Laplacian bias. That is an MSE of
up to about 4e-10, which is not negligible against a 5e-8 target.

Fixes:
- use the exact Phi2(0, 0, rho) = 1/4 + asin(rho)/(2 pi) at zero means, or eps <= 1e-14;
- loop over the 2049 inputs. The f5lib batching allocates an (N, n, n) array, about 16 GB.

Timing on this machine: one GC layer at n = 1024 takes 0.85 s (Owen's T 0.28 s per 2e6 calls), so 2049 runs take about
8 core-hours per network. That matches the author's estimate.

### R4. E1b validity: is the production chain's dependence on a honest?

The midpoint needs Psi(a, Sigma) to be the chain's genuine estimate for N(a, Sigma), and needs it to be jointly
homogeneous: Psi(c a, c^2 Sigma) = c Psi. Three parts of the production chain may assume a centred isotropic input:
- the gain/radial initialisation (chi^2_n radial variance 1/(2n));
- the layer-1 source births;
- calibrated or diagonal kappa_4 tables.

If any of these is hard-coded, Delta_a Psi carries an artefact rather than the drift. Before the 2049-run job, E1b
should run two checks:
- (i) homogeneity: Psi(c a, c^2 I) = c Psi(a, I) to 1e-12 in float64;
- (ii) shifted inputs: Psi(a, I) against Monte Carlo truth at |a| ~ 1 in 4 directions. The chain's error there should be
  of the same order as at a = 0.

The rebind of f32 to float64 also has to be complete (est_v29 is large). One float32 leftover would destroy the second
differences: float32 roundoff divided by h^2, summed over 1024 directions, is about 3e-3.

### R5. Smaller points

- **Thm 2.1.** It is true but nearly tautological: Tr T_f / d_l = int f, because phi affine gives equality. For Theorem D's
  chain at temperatures beta_2 t and t, the shared first derivative forces beta_2 E f = E f, so the chain is only
  two-sided valid when the mean is zero. That is the case for Theorem D's centred Ising forms, which says the same thing
  more simply.
- **Thm 5.2(3) "bracketing".** This is the ansatz, not a theorem. Measured, the companion has the opposite sign on only
  32-70% of neurons, and corr(err_C, err_Delta) is +0.16 and -0.05 on two networks. "Sandwich" overstates it.
- **Thm 5.1.** The floor holds for i.i.d. leaves only. Structured designs (control variates using the chain's own
  Hessian) change the constant, not the conclusion. That part is fine.
- **Cost of the exact structural drift.** Propagating n covariance tangents through W^T X W costs 4n^4 FLOPs, about 2048
  units per layer, not 512. This is immaterial (both are far out of budget), but the table should be corrected.
- **The "dressed drift at +40 units" in 6.2/6.4** has no derivation or measurement. It should be removed from the cost
  table or labelled as a placeholder.
- **The local-tangent variant V1 (+145 units).** It is costed for the GC chain. A production-chain V1 drift would need
  tangents of the A, P, M legs and is not costed.
- **The V2 mean-tangent drift.** It costs about 1-2 units per layer and reaches oracle ratios of 0.045-0.28, comparable to
  V1. But its oracle kappa drifts toward 0 with size (-0.51, -0.31, -0.22, -0.17), so no transferable coefficient exists.
  This should be said, since it is the only variant that is cheap enough.

### R6. Novelty

- The ingredients are classical, as the ledger says:
  - Strassen, lift zonoids, Markov-Krein;
  - Lieb-Berezin (the Funk-Hecke Berezin transform on spheres is standard harmonic analysis);
  - Price, Tanaka;
  - the Gaussian interpolation (smart-path) form of (4.1).
- The gain sign rule is Jensen (E A <= sqrt(E A^2)), and the repo's GAC already uses it.
- What looks new is the specific construction: the heat-defect midpoint T_1 = (Psi + Delta_a Psi)/2, and the tower, for
  moment-closure chains. An alphaXiv search on Richardson extrapolation, Gaussian smoothing, moment closure and heat
  equation found nothing equivalent. Calling it "Richardson in the smoothing variance, with derivatives from the
  estimator" is an accurate description. Its relation to zero-noise extrapolation is fair.

### R7. What survives, corrected

**The Euler-Stein companion.** For any jointly homogeneous Gaussian-input moment chain Psi:
- Delta_a Psi(0, I) is a second consistent estimator of E F; the exact functional gives E[F(|x|^2 - n)] = E F, checked
  numerically (`chk1.py`, max deviation 1.3e-3, within Monte Carlo noise);
- its defect tr D = Psi/2 - Delta_a Psi/2 is a ground-truth-free measure of the chain's local failure of the heat
  equation.

**Measured value on the GC chain at n <= 128.** On 12 networks, -tr D predicts the per-neuron error with corr 0.74-0.98
and kappa about -0.93. I reproduced n16 s1, n32 s2 and n32 L16 s1 to 3 digits with fresh truth, at h = 0.02 and 0.05.

**The correct class accounting.** First-order localization weights:
- the kappa_4 path at 1;
- the kappa_3 path at 2;
- the gain's kappa_3 half at 1;
- the gain's kappa_4 half and the readout's own closed walks at 0.

**What it is for in Phase 2.** It is an offline diagnostic, not an estimator: about 400 B per network for production.
Run it only if R4's validity checks pass. Even then, oracle attribution against truth (note XXXI) already exists, so the
incremental information is modest.

**Highest-value action: E0.** It is free, and its sign and gain-shape audit can confirm or refute "no gain component left"
in the production residual.
