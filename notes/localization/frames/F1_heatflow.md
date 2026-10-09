# F1. The heat-flow defect: localization as the exact error functional of a deterministic estimator

Frame F1 of the localization round. Status labels: **proved** (complete proof here), **sketch** (argument given,
details routine but not written), **measured** (numerical, small n, outputs quoted), **conjecture**.
Code: `loc/f1code/` (`f1_lib.py` Gaussian-closure and mean-field chains; `f1_k3.py` dense third-cumulant chain;
`f1_defect.py`, `f1_defectK.py`, `f1_k3defect.py` defect vs error; `f1_profile.py`, `f1_profile_k3.py` localization
profiles; `f1_dirs.py`, `f1_cv.py`, `f1_k3dirs.py` directional probes; `f1_split.py`, `f1_split2.py` layer
decomposition; `f1_lono.py` held-out merge; `f1_ladder.py` Euler-Stein ladder). Raw outputs `loc/f1_*.txt`, `loc/f1out/`. Tiny networks and Monte Carlo truth are the F3 files
`loc/data/net_n{64,128}_L{8,16}_s*.npz` (2e6 samples) plus my own `loc/f1data/` (n = 48, 3e6; n = 256, 1.2e6).

## 0. Summary

1. **Exact error representation (proved).** For any estimator E(m, Sigma) of u(m, Sigma) = E F(m + Sigma^{1/2} Z)
   that is exact on point masses, the error at the challenge law is minus the integral of its heat-equation defect
   D = d_Sigma E - (1/2) Hess_m E along any Gaussian localization path (Eldan, Chen-Eldan control, or a discrete
   tilt). With homogeneity the integral descends to the leaf space of the scaling foliation: the localization becomes
   a *repulsive* Ornstein-Uhlenbeck flow dy = (y/2) dT + dB, the true functional is its 1/2-eigenfunction, and the
   error is the integral of the reduced defect delta = e - y.grad e - Lap e against an explicit Green (transverse)
   measure of total mass one, density (1/2)(1+tau)^{-3/2} over the Gaussian scales N(0, tau I).
2. **Structure (proved).** The heat equation for the mean lifts, by Hopf-Cole, to an exact *cumulant heat hierarchy*
   (d_Sigma - Hess/2) kappa_k = (1/2) sum_{A|B} grad kappa_A (x) grad kappa_B, whose expected-posterior form is the
   all-orders generalization of the Anari-Koehler-Vuong covariance trickle-down. With homogeneity it yields an exact
   *Euler-Stein ladder* of truth-free constraints, tr_2 c_{k+2} = (1 - k) c_k for the chaos coefficients
   c_k = E D^k F (rung 0 is note XLI's mu = tr H; rung 1: the third chaos is traceless). A finite-moment closure is
   consistent with the law of iterated expectations iff it is affine in moments; so for this problem ARC's iterated
   estimation plus point-mass exactness *forces* the truth, and the error of every deterministic chain *is* its
   integrated violation of iterated estimation.
3. **The measured surprise: the defect at the base point is an almost perfect per-neuron error detector.**
   - *Gaussian closure (dense covariance).* err_i ~ -0.43 delta_i(0) explains 88-96% of the per-neuron error energy
     (77-98% on n = 64, L = 16), with one slope across n = 64, 128, 256 and L = 8, 16 (-0.34 to -0.52, pooled -0.43).
   - *The localization profile.* b(s), the error left after localizing to variance s, behaves as
     b_i(s) ~ err_i s^p with p = 1.1-1.2 at s near 1, and the slope is -1/(2p), as the Duhamel formula requires.
   - *The dense third-cumulant chain K3* (10x lower MSE than Gaussian closure, its error at the next order). Here
     delta explains 74-88% with one slope, -0.18 to -0.22 (p = 2.3-2.8): a better chain makes its error closer to the
     base law. Held out, with one pooled coefficient, the merge e - 0.20 delta cuts K3's MSE by 81% on average (6
     networks), and GC's by 76% (e - 0.48 delta, 4 networks).
   - *Reading (conjecture).* The localization acts on annealed kink statistics while the error amplitude is quenched,
     so the profile factorizes and delta is a quenched detector with an annealed coefficient.
4. **The cost results (Proposition 8 proved, Proposition 9 sketched) and the scaling measurements decide the design.** delta(0) needs the chain's Laplacian
   in the input mean.
   - *Random probing gets worse with n.* One directional defect costs four chain evaluations (0.4-0.8 B for the
     production chain, already over budget). Its single-direction relative noise is 1.5-1.9 at n = 64, 2.1-3.5 at
     n = 128 and 3.1 at n = 256.
   - *Partial localization loses its share.* The share of tr D captured by the top-k response directions falls with n
     (top-1: 0.28, 0.17, 0.14; top-8: 0.55, 0.42, 0.30). This prices every partial-localization quadrature (frame F3)
     at n = 1024: it removes at most the captured share, and it predicted F3's measured n = 64 gain (0.27 against
     0.28).
   - *The analytic defect telescopes exactly over layers* (Proposition 8, verified to four digits). For the Gaussian
     closure it is carried at depth by the off-diagonal covariance map's cross term (the first-chaos-projected (2,1)
     slice) and its diag-type term. The per-neuron mean and variance readouts cancel, and the traced pair Gram (the
     V2 n^4 wall) contributes at most 3%. The cross term factors through the chain's sources with one first-chaos
     leg: n^3 per source-layer, the cost class of a young-hub family.
   - So the defect does not make the next order cheaper, but it does not hit the wall either. It selects the
     localization-visible combination of the next order and weights it by -1/(2p).
5. **Verdict and component.**
   - *Probing is out.* No runtime component works by probing: one directional defect costs four chain evaluations.
   - *Proposed component.* The *defect-Richardson merge* e - c delta, c = 1/(2p) calibrated once per chain class,
     with delta computed analytically by the layer telescoping: first-chaos legs on the sources and one hub-type
     product per source-layer, an estimated +35-70% of the bill.
   - *Prediction.* If the production chain behaves like the small-n kappa_3 chain (X* = 0.6-0.8 of its error
     detected), adjusted MSE goes from 3.14e-9 to 0.85-2.1e-9.
   - *Gate.* One offline experiment decides it (section 10: 1028 float64 chain evaluations, under an hour on one
     large instance), followed by an exact-defect localization on one network.
   - *Main risk.* The production chain's local defects sit one level up (tangents of its kappa_3 sources and kappa_4
     closures), and their cost in the source representation is not established here.
   - *If the gate fails (X* < 0.3).* F1 closes as explanation. It explains why the counterterms cannot reach the
     incoherent error, why partial localization dies with n, and why closures cannot be fixed by consistency more
     cheaply than by carrying the next order.

## 1. Setting and notation

- Input X ~ N(0, I_n); F: R^n -> R^n the bias-free ReLU MLP (row convention x_{l+1} = relu(x_l W_l)); F is
  positively 1-homogeneous and Lipschitz.
- u(m, Sigma) := E F(m + Sigma^{1/2} Z), Z ~ N(0, I). The truth is u(0, I).
- An *estimator* is any map E: R^n x S_+^n -> R^n (a chain run on the input law N(m, Sigma)); e(y) := E(y, I).
- E is *exact on point masses* if E(x, 0) = F(x); every cumulant chain is (all cumulants vanish, every Gaussian or
  Edgeworth readout reduces to relu of the mean). E is *homogeneous* if E(cm, c^2 Sigma) = c E(m, Sigma) for c > 0;
  every chain built from Gaussian relu moments, Mehler series and Edgeworth corrections is (all closure functions
  depend on mu/sigma and scale with sigma).
- The heat operator H := d_Sigma - (1/2) Hess_m (matrix valued; <A, H E> = sum_jk A_jk (dE/dSigma_jk -
  (1/2) d^2E/dm_j dm_k)). u solves H u = 0. The *defect* of E is D := H E.
- Eldan's localization: Y_t = tX + B_t; X | F_t ~ N(m_t, Sigma_t) with Sigma_t = I/(1+t), m_t = Y_t/(1+t),
  dm_t = Sigma_t dW_t. Write s = 1/(1+t) (posterior variance) and tau = t = (1-s)/s.

## 2. The exact representation

**Theorem 1 (Duhamel formula along a localization; proved).** Let (m_t, Sigma_t) solve dm_t = Sigma_t C_t dW_t,
dSigma_t/dt = -Sigma_t C_t C_t^T Sigma_t from (0, I), with C_t adapted, locally bounded, and Sigma_t -> 0 a.s.
(isotropic Eldan: C = I; any Chen-Eldan linear-tilt control). Then N(m_t, Sigma_t) is the conditional law of
X = lim m_t ~ N(0, I) given F_t. Let E be C^{2,1} in (m, Sigma) on Sigma > 0 with derivatives of polynomial growth.
Then for every T,

    E[ E(m_T, Sigma_T) ] - E(0, I) = - E int_0^T < D(m_t, Sigma_t), Q_t > dt,     Q_t := Sigma_t C_t C_t^T Sigma_t,

and if E is exact on point masses and E(m_T, Sigma_T) -> F(X) in L^1,

    truth - E(0, I) = - E int_0^infty < D(m_t, Sigma_t), Q_t > dt.                                   (2.1)

*Proof.* Ito: dE = grad_m E . dm + <d_Sigma E, dSigma> + (1/2) <Hess_m E, d<m>>, with d<m> = Q_t dt = -dSigma. So
dE = (martingale) - <d_Sigma E - (1/2) Hess_m E, Q_t> dt; the stochastic integral is a true martingale under
polynomial growth because m_t has Gaussian tails. Take expectations, then T -> infinity with the boundary condition.
The conditional-law statement is the Gaussian posterior computation (Kalman filter with observation dY = C^T X dt +
dW). QED.

- **Exactness.** D = 0 on the path and the boundary condition give E(0, I) = truth; conversely H u = 0, so u is the
  unique estimator (with polynomial growth) exact on point masses that is defect-free everywhere.
- **Deterministic form (isotropic path).** With g(s) := E_{m ~ N(0,(1-s)I)} E(m, sI) one has g(1) = E(0, I),
  g(0+) = truth and g'(s) = E[tr D(m, sI)] (the m-average contributes -(1/2) E Lap E by the heat semigroup), so
  truth - E(0, I) = - int_0^1 E_{m~N(0,(1-s)I)} tr D(m, sI) ds. The *localization profile*
  b(s) := truth - g(s) = E_m[u(m, sI) - E(m, sI)] is the error left after localizing to variance s; b(1) = err,
  b(0) = 0, b'(s) = -E tr D(m, sI), so b'(1) = -tr D(0, I) (= -delta(0)/2 for homogeneous E, section 3). If the
  profile is a power, b(s) = err s^p, then **err = -delta(0)/(2p)**: the base-point defect determines the error up to
  the profile exponent.
- **Discrete linear tilt (Anari-Koehler-Vuong).** For any measure-valued martingale nu_{k+1} = (1 + <x - m_k, Z_{k+1}>)
  nu_k ending at point masses and any functional E of laws exact on point masses, truth - E(nu_0) = - sum_k
  E[ E(nu_k) - E[E(nu_{k+1}) | F_k] ]: the error is the sum of the one-step iterated-expectation gaps. (2.1) is the
  continuum limit; the Gaussian path is the one on which the chain can actually be evaluated, since the chain only
  accepts Gaussian inputs.
- **Regularity in practice.** The chains of this repo are real-analytic in (m, Sigma) on Sigma > 0 (Phi, phi, Hermite
  polynomials of mu/sigma); the only issue is s -> 0 at neurons with mean exactly 0, a null set. The boundary
  condition was checked numerically: b(0.03) = 0 within Monte Carlo noise on all seven n = 64 networks
  (`f1_profile_*.txt`, |b(0.03)| <= 0.02 |err| against noise 0.01-0.02 |err|).

## 3. Homogeneity: the leaf-space (ray) reduction

The scaling group R_+ acts on Gaussian laws by (m, Sigma) -> (c m, c^2 Sigma). Its orbits are the leaves; the slice
Sigma = I is a transversal. u is leafwise linear (u(c.) = c u), so it is determined by v(y) := u(y, I).

**Corollary 2 (proved).** For homogeneous E, put e(y) = E(y, I) and the *reduced defect*
delta(y) := e(y) - y.grad e(y) - Lap e(y). Then:
1. tr D(m, sI) = (2 sqrt s)^{-1} delta(m / sqrt s).
2. truth - e(0) = -(1/2) int_0^infty (1+tau)^{-3/2} E_{y ~ N(0, tau I)} delta(y) dtau, and the weight
   (1/2)(1+tau)^{-3/2} has total mass 1.
3. In logarithmic localization time T = log(1+t), y_T := Y_t / sqrt(1+t) is the repulsive Ornstein-Uhlenbeck process
   dy = (y/2) dT + dB_T, with generator A = (1/2)(Lap + y.grad). The truth satisfies (A - 1/2) v = 0, and
   e^{-T/2} v(y_T) is the martingale u(m_t, Sigma_t). For any homogeneous estimator
   truth - e(0) = E int_0^infty e^{-T/2} (A - 1/2) e(y_T) dT = -(1/2) E int_0^infty e^{-T/2} delta(y_T) dT.

*Proof.* (1) E(m, sI) = sqrt s e(m/sqrt s); differentiate in s and twice in m. (2) Theorem 1 with C = I:
Q_t = (1+t)^{-2} I, y = m_t sqrt(1+t) ~ N(0, tI). (3) dy = (1+t)^{1/2} dm + (1/2) y dt/(1+t) and d<y> = dt/(1+t) = dT;
apply Ito to e^{-T/2} e(y_T). At y = 0, (A - 1/2) v = 0 reads Lap v(0) = v(0). QED.

- **Reading.** The localization descends to the leaf space as a diffusion; the truth is a harmonic (eigen) object for
  the transverse dynamics, with the network as boundary data at infinity (v(y) ~ F(y) as |y| -> infinity, where the
  input noise is relatively negligible). The *Green measure* G(dy) = (1/2) int (1+tau)^{-3/2} N(0, tau I)(dy) dtau
  is the transverse measure; err = -<G, delta>. In n dimensions |y|^2 / n concentrates at tau, so G is effectively a
  measure on the signal-to-noise ratio tau of the input.
- **Heavy tail.** P(tau > T) = (1+T)^{-1/2} under G. A Taylor expansion of delta around y = 0 therefore diverges term
  by term (E_{N(0,tau)} delta = sum_k (tau/2)^k Lap^k delta(0)/k! against moments int tau^k (1+tau)^{-3/2} = infinity
  for every k >= 1). Whether the error is "made" near the base law is a property of the estimator, measured in section 8.
- **The scale direction carries no information.** d_s E(0, sI) = e(0)/2 identically, by homogeneity. All defect
  information sits in the mean-Laplacian Lap e.

## 4. The Euler-Stein ladder (truth-free constraints)

**Theorem 3 (proved).** Let c_k := E[D^k F(X)] (the k-th chaos coefficient times k!, a symmetric k-tensor, per
output). For 1-homogeneous F and Gaussian X,

    tr_{(k+1,k+2)} c_{k+2} = (1 - k) c_k,     k = 0, 1, 2, ...

*Proof.* v solves Lap v + y.grad v - v = 0 on R^n (Corollary 2(3)). Apply D^k at y = 0, using D^k(y.grad v)(0) =
k D^k v(0); D^k v(0) = c_k and D^k Lap v(0) = tr c_{k+2}. QED.

- Rung 0: E[Lap F] = E F, i.e. mu_i = tr H_i (note XLI's Euler identity, now one rung of an infinite ladder).
- Rung 1: sum_j E[d_j d_j d_k F] = 0 for every k: the third chaos kernel is traceless.
- Rung 2: tr c_4 = -c_2; fully traced, E[Lap^2 F] = -E F.
- **Measured** (`f1_ladder_out.txt`, n = 64, L = 8, 2e6 samples, Hermite-weight Monte Carlo):
  rung 0: slope 0.983, correlation 1.0000 across the 64 outputs (deviation consistent with the common-mode noise of
  the |X|^2 weight); rung 1: rms of E[F X_k (|X|^2 - n - 2)] = 6.9e-3, equal to its noise level (expected 6.2e-3),
  against 2.75e-2 for the first chaos E[F X_k]; rung 2: slope 1.08, correlation 0.998 (noisy weight).
- **Use.** Every rung is a constraint any carried chaos representation must satisfy, at no Monte Carlo cost. The
  chain's bare sources satisfy rung 0 only to correlation 0.73-0.95 (note XLI). Dressing a birth's traced second chaos
  with the full variance (w2 var instead of w2 |l|^2) makes rung 0 hold at Gaussian level, because the chain's mean
  recursion unrolled is sum_b P_{b->L}(w2_b var_b). That is the ladder's reading of why the dressed (Schur) hub beat the
  bare one.

## 5. The cumulant heat hierarchy and the posterior cumulant SDE

**Theorem 4 (cumulant heat hierarchy; proved).** Let kappa_k(m, Sigma) be the joint cumulant tensors of the outputs
F(m + Sigma^{1/2} Z). Then

    H kappa_k = (1/2) sum_{A u B = [k], A, B nonempty} grad_m kappa_{|A|}(A) (x)_sym grad_m kappa_{|B|}(B),   (5.1)

the outer product being over the two input indices of H. *Proof.* For fixed lambda, w = E exp(<lambda, F>) solves
H w = 0; K = log w then satisfies d_Sigma K = (1/2)(Hess K + grad K grad K^T) (Hopf-Cole). Expand in lambda. QED.

- k = 1 is H mu = 0. k = 2: H C = grad mu grad mu^T (symmetrized). Traced at the base with homogeneity (tr d_Sigma
  kappa_k = (k/2) kappa_k), (5.1) gives exact *radial identities*: mu = Lap mu (rung 0), and
  C - L L^T = (1/2) Lap C with L = grad mu = E[grad F] the first chaos: the higher-chaos part of the output covariance is
  half its mean-Laplacian.
- **Tilt identity (proved).** A mean shift of a Gaussian is an exponential tilt, so grad_m kappa_k(F_{i1..ik}) =
  kappa_{k+1}(F_{i1}, .., F_{ik}, X). Hence (5.1) is a statement about mixed (output, input) cumulants.

**Theorem 5 (posterior cumulant SDE; proved).** Let nu be any law on R^d with exponential moments, localized by the
linear tilt d nu_t(z) = nu_t(z) <z - mu_t, C_t dW_t>. Its cumulant tensors satisfy

    d kappa_k = kappa_{k+1}[ . , C dW ] - (1/2) sum_{A u B = [k]} < kappa_{|A|+1}[A, .], C C^T kappa_{|B|+1}[B, .] > dt.

*Proof.* The posterior moment generating function M_t(lambda) = E_t e^{<lambda, z>} is a martingale with
dM = M (grad K(lambda) - grad K(0)) . C dW; apply Ito to log M and read off the coefficient of lambda^k / k!. QED.

- k = 2 gives E[dSigma] = -Sigma C C^T Sigma dt: the AKV covariance trickle-down (Sigma - E Sigma_next = Sigma R Sigma).
  Theorem 5 is its all-orders version: the expected loss of the k-th cumulant is a sum of pair contractions of two
  cumulants whose orders add to k + 2, across the observation channel. In the token language of note XL every drift term is a
  two-token contraction whose internal line is the observation.
- Under input localization the same formula holds for the law of any layer's pre-activations with the observation
  legs replaced by input legs (tilt identity). This is the bridge between the localization and the chain's tables:
  the localization moves a layer's cumulant of order k through that layer's cumulant of order k+1 with one input leg.

## 6. Iterated-expectation consistency, closures, and ARC's law of iterated expectations

**Proposition 6 (consistent closures are affine in moments; proved).** Let an estimator read a finite vector of
linear statistics M(nu) = (int phi_1 d nu, .., int phi_N d nu) of a law and output Psi(M(nu)). It satisfies the
iterated-expectation identity E(nu) = sum_k p_k E(nu_k) for every finite disintegration nu = sum_k p_k nu_k inside a
convex class N iff Psi is affine on M(N). *Proof.* M is affine under mixtures; the identity is Jensen's equality for
Psi on the convex hull. QED.

- **Consequence for this problem.** A chain is a function of the first two moments of the input law. If it satisfied
  iterated estimation along every Gaussian disintegration and were exact on point masses, then D = 0 and E = u by
  Theorem 1; and if it were a function of finitely many moments of a layer law it would make that layer's readout
  affine in those moments, i.e. replace relu by a polynomial. So **no efficient deterministic estimator of this
  problem satisfies iterated estimation; its error is exactly its integrated iterated-estimation gap (2.1).**
- **Gaussian closure is consistent along Gaussian directions.** The Gaussian family is closed under Gaussian
  localization and the Gaussian readout is exact on it, so the Gaussian closure's defect comes only from the
  non-Gaussian part of the localization tangent (section 7). This is the precise sense in which "conditioning is not
  freezing" (the brief's q A^2 q - (q A q)^2 = q A (1 - q) A q): a closure freezes the posterior's higher cumulants,
  the localization does not, and the frozen excursion is the defect.
- **ARC (Christiano, Hilton, Lincoln, Neyman, Xu, arXiv 2410.01290).** Their iterated estimation G(G(Y|Pi)|Pi') =
  G(Y|Pi') is stated subjectively (the outer G estimates its own output) and they turn to objective accuracy over a
  distribution because the subjective form is unwieldy. Here the arguments are input localizations (Pi' <-> F_s,
  Pi <-> F_t), G(Y|Pi) is the chain run on the posterior law, and the outer G can be replaced by the *exact* average
  over the known Gaussian observation law: iterated estimation becomes the objective statement "E(m_t, Sigma_t) is a
  martingale", which is 1-accuracy over the observation distribution and needs no truth. Two consequences:
  1. ARC's desideratum is unattainable for any efficient estimator here (Proposition 6), and its quantitative
     violation is the error, exactly (Theorem 1).
  2. ARC's principle of unpredictable errors is violated in a computable way: delta is a truth-free quantity that
     predicts the error (section 8). Their remedy for a predictor of the error is an OLS merge (their Proposition
     3.4); here the merge is e - c delta with c = 1/(2p) set by the localization profile. After the merge the residual
     is orthogonal to delta (self-accuracy with respect to the defect).

## 7. The chain's defect in chaos-graded form

**Proposition 7 (local defect of a readout; proved for the stated readouts).** Let a neuron's next-layer mean be read
by a closure h(kappa_1, .., kappa_K) of its own marginal cumulants (Gaussian: K = 2; Edgeworth: K = 3, 4). Write
R_k for the right side of (5.1) evaluated on the chain's tables and D_k := H kappa_k^chain - R_k for the chain's
upstream defect in the k-th cumulant. By the chain rule for H,

    tr H (h o kappa) = [ sum_k d_k h . tr R_k - (1/2) sum_{k,k'} d_k d_k' h < grad kappa_k, grad kappa_k' > ]  +  sum_k d_k h . tr D_k .

The bracket is the readout's *local* defect; it vanishes identically for the exact readout (which is linear in the
law), and the last sum transports the upstream defects. For the Gaussian readout h_G(mu, v) = E relu, with
w_k = E f^(k):

    tr H h_G = -(w3/2) <grad mu, grad v> - (w4/8) |grad v|^2,

since h_v = w2/2 cancels the (1/2) h_mumu |grad mu|^2 term (the Gaussian closure is consistent along Gaussian
directions). *Proof.* Chain rule for H on a composition plus (5.1) for k = 1, 2; h_mu v = w3/2 and h_vv = w4/4 follow
from the Gaussian heat equation d_v = (1/2) d_mu^2. QED.

- **Chaos grading.** grad mu_i = L_i (first chaos) and grad v_i = kappa_3(z_i, z_i, X) = 2 H_i L_i + (third-chaos
  terms). So the Gaussian readout's defect is -w3 L^T H L - (w4/2) |H L|^2 + ...: the open walk of length two (the
  open part of D3) and of length three (the path class P4 of note XLI). Closed walks (tr H^3, the joint-gate triangle;
  tr H^4) enter the defect only through the d_Sigma part of H acting on their omission, never through the tangent
  Gram, because d_m tr H^k = 0 at the second-chaos level: the heat operator links them to Laplacians of open walks,
  Lap_m (L^T H^{k-2} L) = 2 tr H^k (for the chi-square model, proved by direct differentiation).
- **Which truncated cumulants feed it, at which order.** If a chain carries kappa_{<= K} and reads them with
  Edgeworth closures, its local defect is fed by kappa_{K+1} with one input leg (the "input-leg slices"
  kappa_{K+1}(z^K, X)) and by the closed walks of order K+1: relative order n^{-(K-1)/2}. For the Gaussian closure
  (K = 2) the feed is kappa_3(z, z, X): the D21 table with one leg moved to the input. For a kappa_3 chain (K = 3) it
  is kappa_4(z, z, z, X) (the (3,1) slice with an input leg, i.e. the K31 sector that note XXXIX found 74% off and
  flex-limited) plus the triangle and the omitted GC feedback.
- **What this says about the counterterms.** The defect is per-neuron and quenched; a per-layer amplitude can only
  rescale its coherent part. The repo's measured ceiling of the counterterms (-3.8% held-out) is consistent with a small coherent share.

## 8. Measurements (small n; truth by Monte Carlo)

Estimators: GC (Gaussian closure, dense covariance, Mehler order 14, exact diagonal), MF (mean field, diagonal
covariance), K3 (dense third-cumulant chain of `f1_k3.py`: first-order Edgeworth readouts, Mehler covariance with all
first-order kappa_3 terms, Gaussian kappa_3 of relus by the trivariate Hermite expansion to total order 5 with exact
coincident slices, first-order gate transport of kappa_3; it omits every kappa_4 effect and the GC1/GC2 feedback).
Lap e(0) by central differences along the n axes, h = 0.05 (h = 0.1 agrees to 3e-4-7e-4 of rms delta).

### 8a. The localization profile b(s) (GC; `f1_profile_n64L{8,16}_gc.txt`)

b(s) = E_m[u(m, sI) - E(m, sI)], m ~ N(0, (1-s)I), 400 outer m (antithetic) x 300 inner samples; beta(s) =
<b(s), err>/|err|^2 (fraction of the error still present after localizing to variance s), resid = |b - beta err|/|b|.

| net | beta(0.1) | beta(0.2) | beta(0.35) | beta(0.5) | beta(0.7) | beta(0.85) | resid at 0.85 |
|---|---|---|---|---|---|---|---|
| n64 L8 s0 | 0.040 | 0.107 | 0.239 | 0.453 | 0.671 | 0.840 | 0.056 |
| n64 L8 s1 | 0.012 | 0.061 | 0.195 | 0.369 | 0.673 | 0.819 | 0.036 |
| n64 L8 s2 | 0.015 | 0.085 | 0.170 | 0.339 | 0.619 | 0.837 | 0.075 |
| n64 L8 s3 | 0.051 | 0.081 | 0.254 | 0.406 | 0.639 | 0.834 | 0.063 |
| n64 L16 s0 | 0.094 | 0.187 | 0.357 | 0.538 | 0.733 | 0.873 | 0.071 |
| n64 L16 s1 | 0.035 | 0.110 | 0.257 | 0.436 | 0.693 | 0.848 | 0.022 |
| n64 L16 s2 | 0.089 | 0.182 | 0.299 | 0.508 | 0.713 | 0.889 | 0.042 |

- **The error is made near the base law.** Half of it is made between s = 0.5 and s = 1 (tau in [0, 1]); less than
  10% below s = 0.1 (tau > 9). The heavy tail of the Green measure is irrelevant: the defect decays fast in tau.
- **The profile factorizes.** b_i(s) ~ err_i phi(s) to 2-8% (resid) for s >= 0.7, and phi(s) ~ s^p with local
  exponent p = 1.07-1.23 (L = 8) and 0.72-1.02 (L = 16) at s = 0.85.
- **Consistency with the Duhamel identity.** b'(1) = -tr D(0, I) = -delta(0)/2. From 8b, <-delta/2, err>/|err|^2 =
  1.01 on n64 L8 s0; the profile's finite difference (1 - beta(0.85))/0.15 = 1.07. The two independent computations
  (Laplacians of the chain vs. Monte Carlo over localized laws) agree.

### 8b. The base-point defect predicts the per-neuron error (`f1_defect_n64*.txt`, `f1_dirs_*.txt`)

Fit err ~ a delta (through the origin); "expl" = 1 - |err - a delta|^2/|err|^2; "demeaned" removes the common mode
from both (the incoherent part only).

| estimator | nets | MSE(err) | slope a | expl | corr | demeaned expl |
|---|---|---|---|---|---|---|
| GC | n64 L8 s0-3 | 5.5-7.4e-4 | -0.47, -0.42, -0.34, -0.43 | 0.95, 0.95, 0.95, 0.92 | -0.97, -0.97, -0.97, -0.93 | 0.94, 0.94, 0.93, 0.87 |
| GC | n64 L16 s0-2 | 3.9-6.4e-4 | -0.52, -0.49, -0.41 | 0.84, 0.98, 0.77 | -0.92, -0.99, -0.83 | 0.84, 0.97, 0.69 |
| GC | n128 L8 s0-2 | 1.1-1.4e-4 | -0.45, -0.44, -0.42 | 0.88, 0.92, 0.96 | -0.93, -0.96, -0.97 | |
| GC | n256 L8 s0 | 4.8e-5 | -0.41 | 0.95 | -0.97 | |
| GC | n48 L8 s0-3 | 3.7-6.0e-4 | -0.57, -0.53, -0.48, -0.36 | 0.87, 0.87, 0.89, 0.59 | -0.91, -0.93, -0.95, -0.80 | |
| MF | n64 L8 s0-3 | 1.6-5.2e-3 | -1.08, -1.05, -0.82, -0.96 | 0.85, 0.84, 0.91, 0.88 | -0.92, -0.91, -0.95, -0.92 | 0.85, 0.83, 0.91, 0.85 |
| MF | n64 L16 s0-2 | 1.3-7.0e-3 | -1.47, -1.54, -1.06 | 0.81, 0.75, 0.82 | -0.87, -0.76, -0.92 | 0.77, 0.58, 0.84 |
| GC, Mehler order 1-4 | n64 (7 nets) | | identical to order 14 within 0.02 | | | |

- **One coefficient per estimator class.** GC: -0.41 to -0.47 on 9 of 11 networks, pooled about -0.43, stable from
  n = 64 to 256. MF: about -1. With one global coefficient instead of the per-network fit the explained share drops
  by at most 0.07 (it is |err - c delta|^2 = |err - a delta|^2 + (a - c)^2 |delta|^2).
- **The merge is a two-estimator blend.** truth ~ e + a delta = (1 + a) e - a Lap e, i.e. 0.57 e + 0.43 Lap e for
  GC. The chain's own mean-Laplacian Lap e is the chain's estimate of E[Lap F] = E F (rung 0 of the ladder), an
  estimator of the same truth whose error is about -1.33 err: anti-correlated with the chain's, which is why the blend
  works.
- **Why the coefficient differs by class (conjecture with a heuristic).** -1/(2p), where p is the localization
  exponent of the leading omitted diagram: localizing to variance s keeps a fraction ~ s^{1/2} of neurons within one
  noise width of their kink, so a diagram with v kink vertices (factors E f^(k), k >= 2) scales as s^{v/2}. MF omits
  the off-diagonal covariance read by one kink (v = 1, p = 1/2, slope -1); GC omits kappa_3 born at one kink and read
  at another (v = 2, p = 1, slope -1/2). Measured: -1.0 and -0.43.
- **Why the profile factorizes (conjecture).** The localization changes the input's signal-to-noise ratio tau, which
  acts on annealed, self-averaging statistics (kink densities per layer), while the per-neuron amplitude of the error
  is quenched (it sits in the downstream weights, which localization does not touch). Hence b_i(s) ~ err_i phi(s) up
  to O(n^{-1/2}), and delta_i is a quenched error detector with an annealed (universal) coefficient. The measured
  residual of the factorization (2-8% at s = 0.85) is consistent with this.

### 8c. Can the defect be probed? (`f1_dirs_*.txt`, `f1_profile_*.txt`)

Directional defect d(v) = <v v^T, d_Sigma E> - (1/2) v^T Hess_m E v (four chain evaluations, all outputs at once);
Hutchinson estimate n E_v d(v) of tr D; "rel. noise" = single-direction rms error / rms tr D. Share = projection
coefficient of sum_{k<=K} d(u_k) onto tr D for the top-K left singular vectors u_k of the true output Jacobian.

| n (L = 8, GC) | rel. noise per random direction | J energy top-1 / top-8 | tr D share top-1 / top-2 / top-4 / top-8 |
|---|---|---|---|
| 64 (4 nets) | 1.53-1.92 | 0.45-0.54 / 0.95-0.97 | 0.24-0.31 / 0.33-0.40 / 0.45-0.47 / 0.49-0.64 |
| 128 (3 nets) | 2.14-3.49 | 0.32-0.36 / 0.84-0.86 | 0.16-0.18 / 0.21-0.25 / 0.28-0.39 / 0.41-0.43 |
| 256 (1 net) | 3.12 | 0.18 / 0.62 | 0.14 / 0.17 / 0.23 / 0.30 |

- **The covariance-direction derivative is an indispensable control variate.** Probing only the mean Hessian
  (Lap e = n E_v v^T Hess v, two chain evaluations per direction) and using the exact homogeneity value
  tr d_Sigma E = e/2 gives single-direction noise 7.1 rms tr D on n64 L8 s0, against 1.81 for the directional defect
  d(v), whose two parts correlate 0.93 across directions (`f1_cv.py`). For a chain whose defect is 1e-4 of its output
  (the production chain) instead of 5e-2 (GC at n = 64), the Hessian-only probe is useless and the control variate
  carries the whole measurement.
- **Random probing gets worse with n.** |D_i|_F / tr D_i ~ rel.noise / sqrt 2 grows from about 1.2 to 2.2: the defect
  matrix is indefinite and its trace is a small residue of cancelling parts. Extrapolated to n = 1024: 4-6 per
  direction, so a useful per-neuron estimate (noise below 0.3 rms) needs K >= 200 directions, about 600-800 chain
  evaluations (3-4 per direction).
- **Response-directed localization is priced by the defect share.** Localizing a subspace U first and integrating it
  exactly (Gauss-Hermite over U, frame F3) replaces the U-part of (2.1) by an exact quadrature, so it removes at most
  the share of tr D in U. At n = 64 the top direction holds 0.28 and F3's k = 1 localization measured an MSE ratio of
  0.53 (amplitude 0.73, removal 0.27): the prediction matches. The share falls roughly like n^{-1/2} (0.28, 0.17,
  0.14 for top-1; 0.55, 0.42, 0.30 for top-8). Extrapolated to n = 1024: top-1 about 0.07, top-8 about 0.15. On the
  official networks the layer-1 frame regression of `../../notes/frontier-tests` already measured R^2 0.007 for the
  top 8 coordinates. **Prediction: partial localization along k <= 8 directions removes less than 25% of the MSE at
  n = 1024, at a cost of at least q^k chain evaluations.**

### 8d. The third-cumulant chain K3 (`f1_k3def_*.txt`, `f1out/k3def_*.npz`)

K3 against truth: MSE 7.96e-5 and 4.11e-5 on n64 L8 s0, s1 (GC: 7.38e-4, 5.53e-4; ratio 0.11 and 0.07). With
kappa_3 switched off K3 reproduces GC to the last digit.

| K3 | MSE(err) | slope a | expl | corr | demeaned expl | implied p = -1/(2a) | rms delta / rms err |
|---|---|---|---|---|---|---|---|
| n64 L8 s0 | 7.96e-5 | -0.204 | 0.884 | -0.919 | 0.844 | 2.45 | 4.6 |
| n64 L8 s1 | 4.11e-5 | -0.200 | 0.800 | -0.845 | 0.713 | 2.50 | 4.5 |
| n64 L8 s2 | 4.40e-5 | -0.184 | 0.795 | -0.910 | 0.829 | 2.72 | 4.9 |
| n64 L8 s3 | 5.94e-5 | -0.217 | 0.830 | -0.912 | 0.831 | 2.30 | 4.2 |
| n48 L8 s0 | 7.61e-5 | -0.212 | 0.739 | -0.846 | | 2.35 | |
| n48 L8 s1 | 9.64e-5 | -0.180 | 0.856 | -0.906 | | 2.77 | |

**Held out, with one coefficient per class** (`f1_lono.py`, `f1_lono.txt`). The slope is fitted on all other
networks and the merge e + a delta is applied to the held-out one:
- K3: leave-one-out slopes -0.195 to -0.204; held-out MSE ratio 0.27, 0.16, 0.12, 0.20, 0.21, 0.18 (mean 0.19, an 81%
  reduction).
- GC (n48, 4 nets): slopes -0.46 to -0.53; ratios 0.16, 0.14, 0.11, 0.54 (mean 0.24).

- **The defect remains an error detector one order up.** For the kappa_3 chain, whose error is ten times smaller
  than GC's and sits at the next order (kappa_4 effects, GC feedback), delta explains 74-88% of the per-neuron error
  energy (71-84% of its incoherent part), again with one slope, now -0.18 to -0.22. With a single pooled coefficient
  the held-out merge removes 73-88% of the MSE.
- **The profile steepens with the order, in the order the kink count predicts.** Implied p = 2.3-2.8 (GC 1.1,
  MF 0.5): a better chain makes its error closer to the base law, and its defect is larger relative to its error
  (rms delta / rms err = 4.2-4.9, about 2p, against about 2.1 for GC). The localization profile measured directly on n48 L8 s0
  confirms it:
  beta(0.5) = 0.32, beta(0.75) = 0.62, beta(0.9) = 0.80 (local exponent 2.1 at s = 0.9; noise 0.15 |err|), against
  beta(0.75) = 0.73-0.77 and beta(0.9) = 0.87-0.93 for GC on the same four n = 48 networks (`f1_profile_gc_n48.txt`,
  `f1_profile_k3_n48.txt`; on n48 s1 the K3 profile was stopped after beta(0.5) = 0.12).
- **Single-direction noise of tr D for K3: 1.86** (n48 L8 s0, 24 directions), the same as GC at n = 64: the defect
  of a better chain is no more concentrated.

### 8e. Which local terms carry the defect (GC; `f1_split.py`, `f1_split2.py`, `f1_split*_n48.txt`)

**Exact layer decomposition (proved in section 9, verified here).** For a chain of layer maps Phi_l, with linear
transport between them, the defect telescopes exactly:

    delta = 2 sum_l S_l . d_l,

where d_l is the *local* defect of layer l's map relative to the hierarchy (5.1), and S_l is the chain's sensitivity
of the output to that layer's post-activation tables. For GC, d_l has three parts:
- the mean readout E relu(z_a);
- the variance readout (diagonal of the post-covariance);
- the off-diagonal covariance map G_ab(mu_a, mu_b, v_a, v_b, C_ab). This one is further split by the tangent Grams
  it uses:
  - "diag-type": <grad mu_a, grad mu_b>, <grad v_a, grad v_b>, <grad mu_a, grad v_b>, which need only the n x n
    per-neuron tangents;
  - "cross": <grad mu, grad C_ab>, <grad v, grad C_ab>, the first-chaos-projected (2,1) slice;
  - "pair": sum_j (d_j C_ab)^2, the traced pair Gram of the gate-covariance (V2) class.

Tangents by central differences, sensitivities by finite differences of the chain restarted at each layer.

| n48 network | reconstruction corr / slope | mean readout | variance readout | cov. map: diag-type | cross | pair (V2) |
|---|---|---|---|---|---|---|
| L = 2, s0 | 1.0000 / 1.000 | +1.000 | 0 | 0 | 0 | 0 |
| L = 3, s0 | 1.0000 / 1.000 | +0.900 | -0.120 | +0.051 | +0.155 | +0.013 |
| L = 5, s0 | 1.0000 / 1.000 | +0.945 | -0.389 | +0.159 | +0.280 | +0.004 |
| L = 8, s0 | 1.0000 / 1.000 | +0.840 | -0.869 | +0.384 | +0.671 | -0.026 |
| L = 8, s1 | 1.0000 / 1.000 | +1.045 | -0.943 | +0.379 | +0.491 | +0.028 |
| L = 8, s2 | 1.0000 / 1.000 | +1.025 | -0.720 | +0.285 | +0.421 | -0.012 |

(Projection coefficients of each part onto delta; they sum to the reconstruction slope.)

- **The decomposition is exact.** It holds to four digits at every depth, which verifies Proposition 7 and the
  telescoping.
- **The V2 wall carries nothing.** The traced pair Gram contributes -0.03 to +0.03 at L = 8.
- **What carries a deep chain's defect.** At depth the mean and variance readouts cancel (+0.84 to +1.05 against
  -0.72 to -0.94). This is why the per-neuron readout terms alone explain only 0-17% (first split, `f1_split_n48.txt`).
  What remains is carried by the off-diagonal covariance map's non-wall terms (+0.71 to +1.06): the cross term
  (projected (2,1) slice) and the diag-type term.
- **The cross term is cheaper than the wall.** By the tilt identity it factors through the chain's sources:
  sum_j d_j mu_a d_j C_ab = sum over births of A_am A_bm (l_m . L_a). The first-chaos leg l_m . L_a is the
  first-chaos cross-covariance between the birth layer and layer l, transported by the mean gate like the P leg. So
  the cross table is a D21-hub-type product with one more leg: n^3 per source-layer, not n^4.
## 9. The cost theorem

**Proposition 8 (exact telescoping; proved).** Let the chain be a composition of nonlinear layer maps Phi_l on
cumulant tables, interleaved with the linear transports z = W^T y. For each table X define its defect relative to the
hierarchy, D(X) := H X - R(X), where R is the right side of (5.1). Then:
1. Linear transports commute with H and R, so D(W^T X W) = W^T D(X) W.
2. For a nonlinear map, D(Phi(X)) = DPhi . D(X) + d(Phi; X), where the local defect is
   d(Phi; X) = DPhi . R(X) - (1/2) D^2 Phi[grad X, grad X] - R(Phi(X)). It vanishes for the exact (law-linear) map.
3. Hence delta = 2 tr D(e) = 2 sum_l S_l . tr d_l, with S_l the chain's linearization from layer l to the output.
*Proof.* The chain rule for H: H(Phi o X) = DPhi . H X - (1/2) D^2 Phi[grad X, grad X]. Substitute H X = R(X) + D(X)
and iterate. Item 1 holds because R is built from products of first derivatives, which transform covariantly under
linear maps. QED. Verified numerically to four digits (section 8e).

**Proposition 9 (cost; sketch).**
- *Generic forward mode.* The "Laplacian state" Lambda_l := sum_j d_j^2 S_l of the chain's state propagates
  linearly, with a source term D^2 Phi_l contracted with the traced tangent Gram sum_j d_j S_l (x) d_j S_l. Exactly
  this costs n chain-equivalents in forward mode, or 2n + 1 chains by finite differences.
- *The factored route.* Proposition 8 needs only the traced Grams in the blocks each layer's Hessian touches. By the
  tilt identity the covariance tangent factors along the chain's own sources: d C_l / d m = kappa_3(z_l, z_l, X) is a
  sum over births of CP terms with the source's two legs and the birth's first-chaos input leg l_m.
  - Traced against a per-neuron tangent (the cross term), it becomes a hub product with the first-chaos leg
    l_m . L_a: n^3 per source-layer.
  - Traced against itself (the pair Gram), it becomes the V2 wall, n^4. It is measured negligible for GC.
  - The per-neuron tangents grad mu_l = L_l (first chaos, mean-gate transport, about 1 unit per layer) and
    grad v_l (the D21 diagonal with an input leg, one hub product per source-layer) complete the list.
- *Cost of an order-K chain's analytic defect.* About one to two extra young-hub families: the young D21 hub is
  about 70 units of the present bill, so +35% to +70%. That is the cost class of carrying the next order, not
  cheaper, but also not the n^4 wall.
- *The diagonal readout terms.* For GC these are the same objects as the production chain's D3 readout (its open
  part) and the V56 path class, with different weights. The defect has w3 L^T H L and (w4/2)|H L|^2; the cumulant
  expansion has (w3/2) L^T H L + (w3/6) tr H^3 and (w4/2)|H L|^2 + (w4/8) tr H^4. A single profile exponent cannot
  reproduce both weights, so -c delta is not the next cumulant order: it is the localization-visible combination of
  it.
- *Caveat for the production chain.* Its defect sits one level up: the tangents of its kappa_3 sources and kappa_4
  closures (kappa_4 and kappa_5 with an input leg). Whether these factor as cheaply as GC's (through first-chaos legs
  on the existing sources) is not established here. That is the main risk of the component in section 10.

## 10. Design: the defect-Richardson merge, its cost, and the decisive experiment

**Component (F1-DR).** e_DR = e - c delta_hat, with c = 1/(2p) fixed once offline per chain class (GC 0.43-0.48;
K3 0.20) and delta_hat the chain's heat-flow defect, computed analytically by Proposition 8: the local defects of every
layer map, transported by the chain's linearization (one adjoint pass of the mean channel). It is the OLS merge of the
chain with its own Euler-Stein estimate Lap e. Costs are in units of 2n^3 = 2^31 FLOPs; the chain itself is 100-200
units.

| route to delta | cost | predicted value at n = 1024 | status |
|---|---|---|---|
| exact finite differences (2n + 1 chains) | 2-4e5 units | the full X* below | offline only |
| random directional defects d(v), covariance-direction control variate (4 chains each) | 400-800 units per direction | per-neuron noise rho/sqrt(K), rho ~ 4-6 | infeasible at any K |
| partial localization along top-k response directions (frame F3), q^k chains | q^k x 100-200 units | removes <= the tr D share: about 0.07 (k = 1) to 0.15 (k = 8) | negative |
| analytic, per-neuron readout terms only | about +16-20 units | about 0: the mean and variance readouts cancel at depth (8e) | negative |
| analytic, full local defect without the V2 pair Gram (readouts + diag-type + cross) | first-chaos legs on the young sources plus one hub-type product per source-layer: about +70-140 units (+35-70%) | reproduces delta to within the V2 share (<= 3% for GC) | the candidate |

**Costed prediction.**
- *Small-n analogue.* For a K3-like chain, the held-out merge with one pooled coefficient removes 73-88% of the MSE
  (section 8d). If the production chain behaves alike, with X* the explained share and kappa_a the share of delta its
  analytic form captures (GC: 0.97-1.03), raw falls by a factor 1 - X* kappa_a: 0.2-0.4 for X* = 0.6-0.8.
- *Cost.* +35% to +70% of the bill (C/B 0.203 -> 0.27-0.35), or less if the young-tier products are shared with the
  existing hub families.
- *Adjusted.* 3.14e-9 x (0.2 to 0.4) x (1.35 to 1.7) = 0.85e-9 to 2.1e-9. Against the present system that is
  -33% to -73%, and the better end is below the leaders (1.1-2.1e-9).
- *Break-even.* The merge loses if X* kappa_a < 0.26-0.41 (the cost factor). Below X* = 0.3 it is not worth deriving.
- *The main risk.* Proposition 9's caveat: the production chain's local defects sit one level up (tangents of its
  kappa_3 sources and kappa_4 closures). If they do not factor through first-chaos legs on the existing sources, the
  cost rises toward the next-order wall and the component dies.

**The decisive experiment (offline, unbilled; the lead runs it on AWS).**
- *Inputs.* The production estimator (`estimator_final_v56.py`, unchanged except the input patch), official networks
  0-3, `truth_off{k}.npz` (per-layer means).
- *Patch.* The input layer must accept a general Gaussian, mu_0 = m W_0 and C_0 = W_0^T Sigma W_0; the later layers
  already handle nonzero means. Run in float64. A directional defect is about 1e-7 of the output, below the float32
  rounding of a 16-layer chain.
- *Script outline.* For each network and r = 1..64 (v_r independent uniform unit vectors, h = 0.5), evaluate
  E(+-h v_r, I), E(0, I +- h v_r v_r^T) and E(0, I). With P_r = v_r v_r^T the directional defect is
  d_r = [E(0, I + h P_r) - E(0, I - h P_r)]/(2h) - [E(h v_r, I) + E(-h v_r, I) - 2 E(0, I)]/(2 h^2),
  tr D_hat = (n/64) sum_r d_r and delta_hat = 2 tr D_hat. Repeat 8 directions at h = 0.25 as a Richardson check in h
  (the truncation is O(h^2 c_4[v,v,v,v]) ~ h^2 mu / n^2). Dump per-layer means at every evaluation, so that the defect
  of every layer's output can be fitted against that layer's truth.
- *Outputs.* Per neuron and per layer: err = truth - e and delta_hat; slope a, explained share, correlation. Also:
  - the across-direction variance, which gives the n = 1024 single-direction noise rho directly;
  - the noise-corrected explained share X*, using corr_true ~ corr_meas sqrt(1 + rho^2/64);
  - the held-out slope (fit on networks 0-1, apply to 2-3).
- *Compute.* 257 chain evaluations per network, 1028 in total, each a 0.2 B-FLOP predict in float64 (about 2x the
  float32 wall time). That is 10-20 core-hours: one 96-core instance for under an hour.
- *Decision rule.*
  - X* >= 0.5 with a slope stable within +-20% across networks: the production chain's error is
    localization-local. Then run the exact defect on network 0 (2n + 1 = 2049 float64 chains with per-layer dumps,
    about 50 core-hours, under an hour on one 96-core instance). Its per-layer defects delta_l = e_l - Lap e_l, layer
    by layer, locate the layers and tables whose local defects (Proposition 8) carry delta.
    - Derive those local defects in the source representation (first-chaos legs on the young and old sources; the
      tangents of the kappa_3 readouts and of the kappa_4 closures' inputs).
    - Validate the analytic delta against the exact one on network 0.
    - Build the merge (V57) only if the analytic form captures >= 70% of delta for <= +70% of the bill. Screen it as
      usual: cold on 16 networks, then scored with counterterms, in adjusted MSE.
  - X* < 0.3: F1 closes as theory, with no runtime component.
  - In between: record the per-layer profile exponents and stop.
- *Prediction for the experiment* (from K3, the nearest small-n analogue): X* = 0.6-0.85, slope -0.15 to -0.25
  (p = 2-3.5, steeper than K3 because the production chain carries more of the kappa_4 sector), rho = 3-6.

## 11. What the frame gives, and what it does not

- **The natural home, for this frame.** The space of Gaussian laws is foliated by the scaling group, and Eldan's flow
  descends to the repulsive OU diffusion on the leaf space.
  - The truth is that diffusion's 1/2-eigenfunction, with the network as boundary data at infinity.
  - Every estimator's error is the pairing of its eigen-defect with the flow's Green (transverse) measure.
  - "The measure on the space of leaves passes to the leaf space" because homogeneity makes the estimand leafwise
    linear.
- **The operator form.** Let P_t be the localization semigroup on functionals of laws, (P_t Psi)(nu) = E Psi(nu_t),
  and G its generator.
  - The truth is G-harmonic, and Theorem 1 is Dynkin's formula: err = -(potential of G applied to G E)(nu_0).
  - Theorem 5 is the carre-du-champ identity for the log-moment functional: the drift of each cumulant is
    -(1/2) Gamma(lower, higher).
  - The brief's operator-valued trickle-down Gamma_s = E_s Gamma_t + Gamma_s(E_t, E_t) is the same identity for
    nested conditional expectations.
  - All of this is a commutative instance of transverse integration. The noncommutative content enters only through
    the algebra the defect acts on: the Wick (chaos) algebra, on which localization acts by a Hopf-Cole (log-state)
    martingale, and the defect's terms sort by the token moves of note XL (V0 readouts, V1 one-leg moves, V2 pair
    walls).
- **Gained.**
  - An exact error functional (Theorem 1) and its leaf-space reduction (Corollary 2).
  - An infinite family of truth-free constraints (Theorem 3).
  - The all-orders trickle-down (Theorems 4-5).
  - A no-go for consistent closures, with its link to ARC (Proposition 6).
  - The chaos form of the defect and its feeds (Proposition 7).
  - A measured, nearly perfect per-neuron error detector with a class-universal coefficient (GC, K3).
  - A pricing of partial localization that matches F3 at n = 64 and predicts its failure at n = 1024.
  - An exact layer telescoping of the defect (Proposition 8). At depth it puts a Gaussian-closure chain's defect in
    the projected-(2,1) cross terms of the covariance map, a hub-class cost, and not in the V2 wall.
- **Not gained.** A cheaper next order. Calibration by consistency is possible in principle (delta is truth-free),
  but it costs either about n chain evaluations or the next order itself (hub-class for GC; unestablished one level
  up). The counterterms'
  truth-fitted amplitudes remain the right tool for the coherent part. The defect shows the incoherent part is
  beyond any per-layer amplitude, and says where it lives.

## Referee report

Adversarial referee for frame F1. Checks were run with my own code, independent of `f1code/`, in `loc/ref_f1/`:
`mf_check.py` (a mean-field chain written from scratch, Corollary 2 and the detector), `radial_check.py` (Theorem 4,
k = 2) and `mf_cal.py` (the detector under truth calibration). I also read the production estimator
`workbench/k3work/estimator_final_v56.py` for the smoothness and input-layer assumptions that the experiment relies
on.

### R1. What is correct (verified)

- **Theorem 1 and Corollary 2.** Correct as stated for C^{2,1} estimators. I re-derived the Ito sign, the
  (2 sqrt s)^{-1} delta(m/sqrt s) reduction, the change of variables to tau and the mass-1 weight. I checked
  Corollary 2 end to end on a tiny net (n = 6, L = 3; MF chain; truth from 4e6 MC; 24-node Gauss-Legendre in s;
  400 antithetic m per node):
  - err = [-0.02338 -0.05292 -0.11597 -0.06774 -0.05304 -0.07828];
  - -int = [-0.02299 -0.05288 -0.11580 -0.06841 -0.05282 -0.07765];
  - the two agree within the MC noise of the s-integral.
- **Theorem 3 (Euler-Stein ladder).** Correct. In one dimension it is the classical recursion of the relu Hermite
  coefficients: c_2 = c_0 = phi(0), c_4 = -c_2, c_6 = -3 c_4.
- **Theorem 4.** Correct. The Hopf-Cole combinatorics hold with ordered pairs (A, B). I checked the radial identity
  C - L L^T = (1/2) Lap C by common-random-number MC (n = 5, L = 2, h = 0.4). Entries agree to 1-3%, which is the
  O(h^2) bias; for example 0.0508 / 0.0498, 0.0581 / 0.0565, 0.0179 / 0.0175.
- **Theorem 5.** Correct. Re-derived via d log M_t; for k = 2 it is the AKV / Eldan covariance drift.
- **Proposition 7 and Proposition 8.** Correct as identities. They are the chain rule for H with exact tangents;
  I checked the Gaussian-readout bracket -(w3/2)<grad mu, grad v> - (w4/8)|grad v|^2 by hand.
- **The detector, replicated independently.** My MF chain gives:

  | net | slope | explained share |
  |---|---|---|
  | n = 64, L = 8, seed 1 | -1.05 | 0.84 |
  | n = 64, L = 8, seed 2 | -0.83 | 0.91 |
  | n = 128, L = 8, seed 3 | -0.81 | 0.92 |
  | n = 6, L = 3 | -0.53 | 0.95 |

  This is consistent with the author's MF rows. The phenomenon is real at small n.

### R2. Mathematical errors and overclaims

1. **Proposition 6 is misapplied, so its stated consequence is overclaimed.**
   - The "only if" (consistent => affine) needs consistency along every finite disintegration inside a *convex*
     class. The chain is defined only on Gaussian laws, which are not a convex class: mixtures of Gaussians are not
     Gaussian. The defect-free functional on Gaussians is u itself, which is not affine in (m, Sigma).
   - For layer laws, the disintegrations induced by input localization are a thin subfamily, so "the readout must be
     affine, i.e. relu replaced by a polynomial" does not follow.
   - What survives is the uniqueness statement in Theorem 1 (Doob/Dynkin): a defect-free estimator exact on point
     masses equals u.
   - "No efficient deterministic estimator of this problem satisfies iterated estimation" is a complexity claim:
     it says u is not efficiently computable. Nothing here proves it. Relabel it as a remark.
2. **The partial-localization bound "removes at most the tr D share" is not proved.**
   - Exact quadrature over U replaces the U-part of the whole path integral (2.1), not just its base-point slope.
     The removed amount can exceed the base-point share.
   - It is a first-order heuristic, matched on a single data point (n = 64, 0.27 vs 0.28).
   - Even at face value the arithmetic is off: an amplitude removal of 0.15 is an MSE removal of
     1 - 0.85^2 = 0.28, not "< 25%".
   - The n^{-1/2} extrapolation rests on one n = 256 network.
3. **The kink-count heuristic for p fails quantitatively at K3.**
   - p = v/2 predicts MF 0.5 and GC 1. The leading K3 omission (kappa_4 / three kinks) would give 1.5, against the
     measured 2.3-2.8.
   - Only the ordering is supported, so the per-class coefficient cannot be predicted and must be fitted with truth.
4. **The coefficient is not a "class constant" under calibration.** The production chain carries truth-fitted
   counterterms (V47_CAL, about 76 per-layer amplitudes). In `mf_cal.py` I added one counterterm-like scalar (gamma on
   the variance channel), fitted it on 2 of 4 networks (gamma = 0.95), and kept everything else fixed (n = 64, L = 8):
   - The explained share survives: 0.74-0.92.
   - The slope moves from -0.76..-0.93 to -0.63..-0.75.
   - So c must be refitted on the calibrated production chain, and its spread grows. The offline experiment does fit
     it, which is correct, but the "fixed once per chain class from small n" language in the estimator summary is
     wrong.
5. **"kappa_a ~ 0.97-1.03" does not support the proposed analytic delta.** The four-digit verification of
   Proposition 8 used exact (finite-difference) tangents and exact finite-difference sensitivities. The costed design
   uses three approximations, and none of them has a measured fidelity:
   - mean-gate transport of the first-chaos legs. This is not the chain's tangent, which also contains
     grad v / grad C feedback and Stein-corrected gates;
   - rank-one (mean-gate) Hadamard multipliers in the tangent transport;
   - the tilt identity applied to *chain* tables. grad_m C^chain is not kappa_3^chain(z, z, X); GC carries no
     kappa_3 at all.

   The kappa_a of the cheap analytic delta is therefore unknown, not about 1. The V2 share (<= 3%) was measured at a
   single size (n = 48, L = 8).
6. **Novelty is over-framed.**
   - Theorem 1 is Dynkin's formula, with the Gaussian heat equation in (m, Sigma) being Price's theorem (1958).
   - "Error = integrated generator defect along the flow, then cancel the leading term by extrapolation" is exactly
     the Talay-Tubaro global-error expansion for weak SDE schemes (1990), and Talay-Tubaro Richardson.
   - The merge e - c delta is classical defect correction (Zadunaisky; Stetter 1978), or a residual-based
     a posteriori estimate in the dual-weighted-residual sense.
   - Theorem 5 is the standard log-MGF computation for Eldan's localization.
   - The "transverse measure" is the commutative occupation (Green) measure. The noncommutative frame contributes
     vocabulary, as section 11 concedes.
   - What is new is the empirical finding: per neuron, quenched, at an essentially universal coefficient for this
     problem. So is the identity "merge = (1 + a) e - a Lap_m e", a blend of the chain with its own rung-0 estimate.

### R3. Cost-accounting errors

1. **Transport mode.**
   - "One adjoint pass of the mean channel gives S_l" is wrong in kind. An adjoint pass gives the sensitivity of
     *one* scalar functional of the output; delta is a 1024-vector. The correct route is a single *forward*
     (tangent-linear) pass that injects tr d_l at every layer and transports the accumulated perturbation, through
     *all* channels the later layers read (mean, dense covariance, and the source tables and their births).
   - Cost of that pass:
     - mean plus dense covariance channels only: about 2-3 units per layer, so +35-50 units;
     - full state: a tangent-linear version of the chain, about 1-2 chain equivalents (+200-400 units).
   - The truncated version's fidelity is unmeasured. The author omits this line entirely.
2. **The local-defect cost is an analogy, not a count.**
   - Cross terms "n^3 per source-layer" over about L^2/2 = 136 (birth, layer) pairs is about 136-270 units if
     young-dense.
   - The +70-140 estimate assumes that most pairs are confined or low-rank, as in the present tiers. Plausible but
     unshown.
   - The production chain's level-up terms (tangents of kappa_3 sources and kappa_4 closures) are uncosted, as the
     author says.
3. **Corrected range.**
   - C/B = 0.203 x (1.5-1.9) = 0.31-0.39 with truncated transport, or 0.203 x (2.3-3.6) = 0.47-0.74 with full
     transport.
   - Adjusted = 3.14e-9 x (1 - X* kappa_a) x cost factor.
   - At X* kappa_a = 0.5: 2.4-3.0e-9 (truncated) or 3.6-5.6e-9 (full). The full-transport case is worse than
     today.
   - At 0.75: 1.2-1.5e-9 or 1.8-2.8e-9.
   - The author's best case (0.85e-9) needs X* kappa_a = 0.8 at +35%, the most favourable corner on both axes.
   - Break-even is X* kappa_a >= 0.33-0.47 (truncated) or 0.57-0.72 (full), not 0.26-0.41.

### R4. The decisive experiment is broken as specified (fixable)

1. **The production chain is not smooth in (m, Sigma).** `estimator_final_v56.py` makes discrete, data-dependent
   choices from mu/sqrt(var) at every layer:
   - the saturation row drop at alpha <= -2.5, with the kept count rounded to multiples of 64
     (SAT_ROUND = 64; argsort-ranked sat_perm / sat_mask);
   - the birth mask (BIRTH_SAT, ranked top-nb; DECK_BIRTH and DECK_ROWS thresholds when set);
   - clip(rr/ref, 0.5, 2) and max(., 0) adaptive gains;
   - rank-k randomized eigen / QR subspace selections.

   A mean shift of h = 0.5 moves each alpha by about 0.01. That flips of order a few rows per evaluation over
   16 x 1024 units, and a 64-row block whenever a count crosses a multiple of 64. The signal is a second difference
   of about h^2 d(v) ~ 1e-7 of the output. Any such jump of 1e-7 or more destroys it, and Theorem 1 does not even
   apply to a discontinuous E, whose distributional Laplacian has jump-layer terms.

   **Fix:** run the experiment on the *mask-frozen* chain. Compute every mask, permutation, rounded count, selected
   subspace basis and clipped or maxed gain once at the base law and hold them fixed for all 257 evaluations. Its
   base value equals the production e exactly, it is analytic, and it is what a runtime analytic delta would
   differentiate anyway.

   **Add a control:** repeat 8 directions unfrozen. If they disagree, that confirms the issue.
2. **The input patch is not a one-liner.** Layer 0 hardcodes C = I: `C_pre = w32^T w32`, and `var = sum w32^2`
   under trim. The f32 casts (`w32`, `.astype(f32)`) are pervasive, so "run in float64" is a code-wide dtype
   patch. It should be validated by checking that the float64 chain reproduces the float32 base output to about
   1e-6 before any differencing.
3. **Stencil bias.**
   - For the exact u the specified stencil is not defect-free. The mean second difference carries
     +(h^2/24) D^4[v^4], while the covariance difference carries only O(h^2 D^6).
   - Hutchinson-averaged and with rung 2 (Lap^2 v = -v), the bias is tr D_hat(u) ~ +(h^2/8) e/(n+2): at h = 0.5,
     about 3% of the expected tr D, and proportional to e. That is small, and the h = 0.25 check catches it, but it
     should be subtracted analytically (it is truth-free).
4. **The noise budget is thin.** With rho = 4-6 and 64 directions the per-neuron noise is 0.5-0.75 rms of tr D. The
   X* correction formula then carries a factor of 1.25-1.56. Near the 0.3 and 0.5 decision thresholds this matters.
   Use 128 directions (1540 evaluations per network, still about 1-2 instance-hours), or report X* with a bootstrap
   interval over directions.
5. **Decision-rule addition.** The exact-defect run on network 0 should also report delta's per-layer decomposition.
   DATA_FINDINGS puts about 2/3 of the final MSE in injections from layers 12-15. If delta is likewise
   late-layer-local, the runtime delta can be truncated to the last 3-4 layers' local defects, cutting the
   local-defect cost by about 4x. The first-chaos legs still run from the input, at about 1 unit per layer.

### R5. Relevance at the 1e-8 level

- **Mechanism.** delta(0) = -err - Lap eps(0), with eps = e - v. "err ~ -delta/(2p)" is therefore the statement
  Lap eps_i(0) ~ (2p - 1) err_i with one ratio: the chain's error, as a function of an input mean shift, decays
  radially and self-similarly.
- **Small-n evidence.** This holds where one omitted diagram class dominates: MF, GC, K3, all n <= 256, and K3 only
  at n = 48-64 and L = 8.
- **The production error is not of that kind.**
  - It is a superposition of uncorrelated per-layer injections (DATA_FINDINGS: shares 3/10/23/11/15/18/21% by layer
    block).
  - Its oracle attribution is about 40% kappa_3 readouts, 15-30% kappa_4, and about 30% other.
  - It contains compression and truncation errors (rank-r legs, saturation drops) and truth-calibrated
    counterterms.
- **Visibility.** Truncation-type and calibration-type errors need not be defect-visible: an error component with
  (1 - y.grad - Lap) eps = 0 at the base is invisible. Different components can have different p, which lowers the
  single-slope share.
- **My prior.** X* for the production chain is 0.3-0.6, centre about 0.45, not 0.6-0.85. At that X*, with
  realistic costs, the merge is roughly break-even in adjusted MSE.
- **Upside.** It can still change the final-layer MSE at the 1e-8 level, since raw is 1.55e-8 and it acts on the
  full per-neuron error. That needs X* kappa_a >= 0.6 together with the truncated transport.

### R6. Verdict

- **The theory is sound but mostly classical.** Dynkin / Price / Talay-Tubaro / defect correction / Hermite
  recursion / Eldan's moment SDE.
- **Proposition 6's consequence and the partial-localization bound are overclaimed.**
- **The real contribution is empirical and independently replicated here.** At small n the base-point heat defect
  delta = e - Lap_m e, using the homogeneity reduction tr d_Sigma E = e/2, is a truth-free per-neuron error predictor
  with a nearly universal coefficient. The coefficient survives calibration only after a refit.
- **The cost section has a mode error and unmeasured approximation fidelities.** Corrected, the plausible outcome
  spans a 2x gain to a 1.5x loss.
- **The experiment needs a mask-frozen float64 chain to be meaningful.** With that correction it is cheap
  (unbilled, about 1-2 instance-hours) and genuinely decisive for X*. Even a positive X* leaves a long, unproven
  build: level-up local defects, approximate tangents, and forward transport.
- **Survives as a gated hypothesis.** Priority is moderate.
