# N2. The matrix-valued defect at a late layer: field localization, two projections, and what a late cut can see

Frame N2 of the NCG round on the heat-defect principle. Status labels:
- **proved**: a complete proof is given here;
- **sketch**: the argument is given and the details are routine;
- **checked**: the identity was verified numerically;
- **measured**: a small-n number with its output file;
- **conjecture**.

Code is in `loc/n2code/` and raw outputs are in `loc/n2out/`. The networks are the F1/F3 He ReLU MLPs (`loc/data`,
`loc/f1data`), with Monte Carlo truth. Every run was single-threaded under `nice`, and none takes more than minutes.

Notation follows F1:
- H = d_Sigma - (1/2) Hess_m is the heat operator, and D = H E is the defect.
- Row convention: z_l = x_(l-1) W_l is the pre-activation at loop index l, and x_l = relu(z_l).
- The output is x_(L-1).
- One unit is 2n^3 FLOPs.

---

## 0. Summary

**Question.** The frame asks four things:
- localize a late layer's Gaussian surrogate N(mu_l, C_l) instead of the input;
- find which pairing of the matrix defect D_l with the anisotropic covariance C_l governs the error, and where
  non-commutativity enters;
- find whether a few collective modes suffice;
- relate the result to F1's full delta, and cost a runtime component.

**Answer in one paragraph.**
- The right late object is the **field (Fisher-metric) trace** tr(D_l C_l) = tr(C^{1/2} D C^{1/2}) of the
  **Gaussian-part (S-frozen)** defect.
- At small n it is an almost perfect, truth-free detector of the error **its suffix creates**:
  - it explains 77-97% of it (median 0.91), with a near-universal slope;
  - it needs only k ~ 4-8 eigen-directions of C_l, because the defect diagonal is flat across modes.
- The non-commuting part of D (about 90% of its Frobenius mass) is provably invisible to it, and is irrelevant.
- The genuine non-commutativity is elsewhere: in the relative position of the input's first chaos and the layer's
  field span (a Halmos pair). That explains why the input-level delta is nearly blind to late-created error, and it
  derives note XLI's bare-to-dressed rule.
- **But a late cut cannot see error made before it** (Proposition 11), and "injected late" does not mean "created
  late". At n = 48-64 the late chains inject 15-79% of their final error in the last half of the layers (last 4 of 8 for K3: 67-79%), yet the
  late field defect detects:
  - for the Gaussian closure, 0-18% of the total error (suffix share 0-8% for 2-4 layer suffixes);
  - for the kappa_3 chain, 0-14% on three of four networks and 48-64% on the fourth (3-4 layer suffixes). Only at
    5 of 8 layers does it reach 44-69%.
- The runtime component is therefore predicted **negative-to-marginal**:
  - the cheap analytic version (two layers, +10-20 units) sees almost nothing;
  - the rolling in-place version (one-step defects at the last 4 layers) removes 12% of K3's MSE held out, but at
    +40-80 production units it loses in adjusted MSE;
  - the useful multi-layer version is infeasible by probing and only sketched analytically.
- One cheap offline measurement on the production chain (section 9) decides it.

**What was established.**

1. **What a late cut localizes (proved).** The chain's state at a cut is s_l = (mu_l, C_l, S_l), where S_l holds the
   kappa_3 sources, the kappa_4 machinery, the masks and the calibration.
   - There are three different localizations of a late layer (section 1):
     - the input localization pushed forward;
     - the field localization of the true layer law (observe z_l through noise of covariance C_l);
     - the **Gaussian-part (S-frozen)** localization.
   - Only the last is truth-free and evaluable at the cut. Its target is heat-exact (H commutes with every
     constant-coefficient cumulant shift).
   - So **every Edgeworth readout with frozen cumulants has zero defect** (Lemma 1). The S-frozen defect lives
     entirely in the suffix's propagation maps.
   - The production readout's non-Edgeworth parts have a closed-form, free "calibration defect": the online
     correction features and the V30 gain term.
2. **The anisotropic identity and the layer delta (proved; checked to 1e-5).**
   - Theorem 2 is Duhamel with a matrix control from (mu, C). In whitened form the rate is tr(D~ R K~ K~^T R), with
     D~ = C^{1/2} D C^{1/2}.
   - Corollary 3 covers the field control K = C^{-1/2} (Sigma_t = C/(1+t)) with homogeneity:

         delta_l := e - mu.grad e - Delta_C e = 2 tr(D_l C_l),   Delta_C = tr(C Hess) the Fisher-metric Laplacian.

   - The leaf-space flow is the repulsive OU started at y_0 = C^{-1/2} mu, with Green weight (1/2)(1+tau)^{-3/2} over
     N(sqrt(1+tau) y_0, tau I).
3. **Which pairing (proved, measured).**
   - **The embedding.** D~ is the geometric (KMS) embedding, the unique operator-mean embedding covariant under the
     ReLU diagonal gauge (Proposition 4). On the commutant of C, though, every Kubo-Ando embedding gives the same
     number, and **spectral localizations see only E_C(D)** (Theorem 5).
   - **The detector.** For GC at n = 64, L = 8 and 16 (5 networks, 18 cut/suffix pairs, 2-6 layer suffixes):
     - the field trace explains **0.77-0.97 (median 0.91)** of the per-neuron suffix error, with slope **-0.545 +- 0.02**
       (err ~ -0.27 delta_l, p_l ~ 1.8);
     - it is the best pairing in 17 of 18 rows and within 0.005 of the best in the other;
     - the alternatives: the geometric mean C # C_1 gives 0.77-0.95, the first-chaos pairing tr(D C_1) 0.67-0.94,
       tr(C D C) 0.63-0.92, tr(D C^{1/2}) 0.61-0.89, and the isotropic tr(D) 0.11-0.59.
   - **Few modes.** d_jj is flat across C's modes, so the top-4 and top-8 modes give 0.72-0.96 and 0.77-0.95.
   - **The non-commuting part.** 86-93% of |D~|_F is off-diagonal and carries none of the power.
4. **Two projections, the keyhole, and dressing (proved, measured).**
   - Input localization sees the layer through P_{H_1}, field localization through P_{Z_l}. These form a Halmos pair
     in generic position.
   - Gamma = C^{-1/2} C_1 C^{-1/2} has median eigenvalue about 0 and maximum 0.4-0.7 at depth. **The input base point
     sees a deep layer through a keyhole of a few collective directions**, and reveals the rest in chaos order along a
     non-commuting path, with Q_l(t) = (1+t)^{-2} sum_k k rho_t^{2(k-1)} C_k (Theorem 7).
   - Replacing P_{H_1} by P_Z in every tangent Gram **is** note XLI's skeleton (dressed) resummation (Theorem 9).
   - With the true kappa_3 at the cut, the dressed generator explains 93-96% of the GC suffix error from the exact
     state and 52-83% of the total, against 46-68% for the bare input-tangent part on the same cuts.
5. **How the late defect relates to F1's delta (proved, measured).**
   - The cut decomposition is tr D_in = (prefix defects transported) + Lambda_l, with
     Lambda_l = <D_l^G, C_1> + (cross + pair) (Theorem 10).
   - For GC, Lambda_l explains 27-94% of the total (growing with suffix length), almost all of it through the cross
     and pair terms: the kappa_3 with an input leg, i.e. accumulated non-Gaussianity. The Gaussian block is small.
   - The S-frozen field defect is the *dressed Gaussian block*. It and F1's delta are near-orthogonal detectors: F1's
     delta explains 1-36% of the suffix error.
6. **Analytic two-layer field defect (proved; checked to 6e-4).**
   - tr(D C_l) = -(E f'''/2) <grad mu', grad v'>_{C_l} - (E f''''/8) |grad v'|^2_{C_l}.
   - It costs a few dense products, because **the mean tangent of a ReLU layer is row-local**: the full Jacobian costs
     one directional derivative.
   - So "k directions instead of n" is the right economy for **probing** (k = 32-64 predicted at n = 1024 from the
     measured spectrum) and the wrong one for analytic tangents.
7. **Rolling, in-place variant (measured; section 7).**
   - Each late layer's one-step field defect is injected into its own readout and carried by the chain.
   - On K3 (n = 48, 4 networks, held-out slope) it gives MSE ratios of 0.99 (last layer), 0.93 (last 2) and **0.88
     (last 4)**, with a per-network spread from -35% to +9%.
   - The slope is stable at -0.6 to -0.9.
   - At the estimated production cost (+10-20 units per layer) this is **adjusted x(1.03-1.22)**. It does not pay
     unless the per-layer cost falls to the GC class (<= 5 units).
8. **Verdict: negative for a runtime component, positive for theory.**
   - The late frame gives an exact, well-conditioned, cheaply truncatable detector of *suffix-created* error, but the
     error is created early.
   - The non-commutativity that matters is the input/field two-projection geometry. The KMS/Kubo-Ando structure
     collapses on the commutant.
   - Registered prediction for the production chain:
     - X_late = 0.08 [0, 0.35] at cut 12 and 0.03 [0, 0.1] at cut 14;
     - top-32/64 capture >= 0.9;
     - d_j flat.
   - The decision experiment (section 9) costs about 15 minutes on one 96-core instance.

---

## 1. What is localized and what is held fixed

### 1.1 The state at a cut and the suffix map

At the pre-activation cut l the production chain (`estimator_final_v56.py`) carries s_l = (mu_l, C_l, S_l):
- mu_l is the pre-activation mean;
- C_l is the dense pre-activation covariance C_pre. At the trimmed last layer only its diagonal is formed;
- S_l holds everything else:
  - the kappa_3 sources (young dense legs A_s, P_s; the old tier-1 basis Q with factors FA, FP; the tier-2 sub-basis U;
    the thin legs Z, Zf and the M legs through S_sep + Rres);
  - the readouts D3 and D21 formed from them;
  - the kappa_4 machinery (the transported diagonal dG, the adaptive lam_l, the V56 Schur-hub path class and star,
    the V38 mixture gains, the V31/V33 pair terms);
  - the masks, ranks and bases (sat_mask, birth masks, randomized bases);
  - the calibration: 103 V47_CAL amplitudes and the 13-feature online mean correction beta_l.

The suffix map is E_l(m, Sigma; S) := the chain from layer l to the output, started from (m, Sigma, S).

### 1.2 Three localizations of a late layer

| localization | what moves at the base point | target that is a martingale | evaluable from the chain? |
|---|---|---|---|
| (I) input Eldan, pushed to the cut | every table: d mu_l = L^T dW (QV C_1 = L^T L); d C_l = (grad_m C_l)[dW] (the kappa_3 with an input leg); d S_l = (grad_m S_l)[dW]; drifts = the chain's own tr H s_l | E F(X) | yes, but only through the prefix's tangents (n chain-equivalents, or forward mode) |
| (F) field localization of the **true** layer law (observe z_l with noise covariance C_l) | posterior cumulants by the cumulant SDE (F1 Thm 5) with control C^{-1}: QV of mu is C; d<mu_a, C_bc> = kappa3_abc; d<C_ab, C_cd> = kappa3_ab. C^{-1} kappa3_cd.; d kappa_3 needs kappa_4 with a field leg, and so on | E G_l(z_l) under the true law | no: it needs the true next-order cumulants |
| (G) **Gaussian-part (S-frozen)** localization of N(mu_l, C_l) | only (m, Sigma) move, by a Chen-Eldan control; S_l is held fixed | u^S(m, Sigma) below | **yes**, at the cut, from the chain alone |

**What (G) holds fixed, precisely.**
- Every table in S_l at layer l is frozen: legs, factors, bases, masks, counts, the transported kappa_4 diagonal,
  lam's ratio, and the calibration constants.
- Every suffix map after l runs normally on the varied (m, Sigma). The births at l and later, the gates of the
  transports, the covariance maps and the hub solves all respond.
- For finite differences, the masks and bases of the suffix are frozen at their base values as well (F1 referee R4.1).
  This is also what an analytic tangent differentiates.

**Lemma 1 (proved).**
- **(a) The Gaussian-part target is heat-exact.**
  - If the layer law is z = g + xi with g ~ N(m, Sigma) independent of a fixed centred xi, then
    u^S(m, Sigma) := E G(g + xi) satisfies H u^S = 0.
  - Formally, for any frozen cumulant tables kappa_{>=3}, the Edgeworth-shifted target
    u^S = exp(sum_k kappa_k . (-d_m)^k / k!) u satisfies H u^S = 0.
- **(b) Edgeworth readouts are heat-exact.** If a readout is h = P_S(d_mu) h_G(mu, v), with P_S a polynomial whose
  coefficients depend only on the frozen tables and h_G(mu, v) = E relu(mu + sqrt(v) Z), then
  (d_v - (1/2) d_mu^2) h = 0.

*Proof.*
- (a) Convolution with N(0, Sigma) is the heat semigroup in m, so E G(g + xi) is the heat flow of the function
  m -> E G(m + xi). The formal version holds because H has constant coefficients in m and the shift operator does not
  depend on Sigma, so they commute.
- (b) h_G solves the heat equation (Price), and P_S(d_mu) commutes with H.

QED.

**Consequences.**
- The S-frozen defect of a one-layer suffix is zero for every exact Edgeworth readout (GC's mean, K3's
  A_0 + A_3 tau/6).
- **A late S-frozen defect is generated entirely by the suffix's propagation maps:** the covariance map, the
  births, the gate transport and the closures.
- For the production readout, only the non-Edgeworth parts survive in the one-layer defect. These are closed-form
  per-neuron functions, free to evaluate (n per layer):
  - **the online mean correction** delta = beta . feats. The heat operator acts on its features as follows:
    - H(1) = H(mu) = H(Phi(alpha)) = 0;
    - H(var) = 1;
    - H(sigma) = 1/(2 sigma);
    - H(alpha) = -alpha/(2 var);
    - H(|alpha|) = -|alpha|/(2 var) for alpha != 0;
    - H(phi(alpha)) = phi/(2 var);
    - H(sigma phi) = phi/sigma;
    - the frozen features D3, |D21| and k4 give 0.
  - **the V30 gain term** -(g/8) mu (Phi - alpha phi), whose defect is (g/4) alpha^2 phi / sigma.

  These formulas were checked by central differences (one-liner in the session; all agree to 1e-6). This
  "calibration defect" is a self-consistency diagnostic of the truth-fitted correction. Nothing here says it predicts
  the residual.

---

## 2. The anisotropic defect identity at a cut

**Theorem 2 (Duhamel at a cut, matrix control; proved).**
- **Setting.**
  - Fix S. Let E(m, Sigma) := E_l(m, Sigma; S) be C^{2,1} with derivatives of polynomial growth.
  - Let (m_t, Sigma_t) solve dm = Sigma K_t dW and dSigma = -Sigma K_t K_t^T Sigma dt from (mu, C), with K_t adapted.
    This is the Kalman-Bucy posterior of N(mu, C) given dY = K^T g dt + dW.
- **Claim.** For every T:

      E[E(m_T, Sigma_T)] - E(mu, C) = -E int_0^T < D(m_t, Sigma_t), Sigma_t K_t K_t^T Sigma_t > dt,

  and the same holds with u^S in place of E and a zero integrand.
- **The profile identity.** Hence (u^S - E)(mu, C) = E[(u^S - E)(m_T, Sigma_T)] - E int_0^T <D, Q_t> dt. With the
  profile b(s) := E[(u^S - E)(m_T, Sigma_T)] at posterior scale s, the base rate is b'(1) = -<D, Q_0>, as in F1.
- **Whitened form.** Put K~ := C^{1/2} K and R_t := C^{-1/2} Sigma_t C^{-1/2} = (I + int_0^t K~ K~^T)^{-1}. The
  integrand is tr(D~_t R_t K~ K~^T R_t), where D~_t := C^{1/2} D(m_t, Sigma_t) C^{1/2}.

*Proof.*
- Ito, exactly as F1 Theorem 1: d<m> = Sigma K K^T Sigma dt = -dSigma, so dE = mart - <H E, Q> dt.
- For u^S, Lemma 1(a) gives H u^S = 0.
- For the whitened form, Sigma_t^{-1} = C^{-1} + int K K^T, and C^{1/2} Sigma_t^{-1} C^{1/2} = I + int K~ K~^T.

QED.

**When the identity reaches the truth.**
- For S = 0 (the Gaussian closure, or any chain whose state at the cut is Gaussian), the localization can be run to
  Sigma -> 0. There E(x, 0) = F_suffix(x), and the identity gives the suffix error exactly.
- For S != 0, a point mass with non-zero kappa_3 is not a law. So only the profile form is available; the detector
  uses only b'(1) and the shape of b near s = 1.

**Corollary 3 (field control, homogeneity, and the layer delta; proved; checked).**
- **Setting.** Take the field control K = C^{-1/2}, so Sigma_t = C/(1+t) and the base rate is tr(D C) = tr D~.
- **Claim 1 (the layer delta).** If the suffix is jointly homogeneous, E(cm, c^2 Sigma; c^deg S) = c E, then

      2 tr(D(mu, C) C) = E - mu.grad_mu E - Delta_C E - Lambda_S E,      Delta_C := sum_ab C_ab d_a d_b,

  where Lambda_S = sum_k k S_k . d_{S_k} is the Euler operator of the frozen tables.
- **Claim 2 (S = 0).** With e(y) := E(C^{1/2} y, C) and y_0 := C^{-1/2} mu:
  - delta_l := e - y.grad e - Lap e at y_0 equals 2 tr(D C);
  - (u - E)(mu, C) = -(1/2) int_0^inf (1+tau)^{-3/2} E_{y ~ N(sqrt(1+tau) y_0, tau I)} delta_l(y) dtau.
- **Reading.** The leaf-space flow is F1's repulsive Ornstein-Uhlenbeck process dy = (y/2) dT + dB, started at
  y_0 = C^{-1/2} mu instead of at 0.

*Proof.*
- The Euler relation is the derivative of homogeneity at c = 1.
- E(m, sC) = sqrt(s) e(C^{-1/2} m / sqrt(s)), so tr(D(m, sC) C) = (2 sqrt s)^{-1} delta_l(C^{-1/2} m / sqrt(s)).
- On the field path the whitened increment is C^{-1/2} dm = dW/(1+t). Then y_t := C^{-1/2} m_t sqrt(1+t) satisfies
  dy = (y/2) dT + dB_T with T = log(1+t), and y_tau ~ N(sqrt(1+tau) y_0, tau I).
- Integrating the rate (1+t)^{-2} tr(D(m_t, C/(1+t)) C) gives the weight (1/2)(1+tau)^{-3/2}, as in F1 Corollary 2.

QED.

- **Checked.** On every GC cut of E1 (`n2out/e1_*.txt`), 2 tr(D C) computed from n directional defects and
  e - mu.grad e - Delta_C e computed from mean derivatives agree to 2e-6-1e-5 relative (correlation 1.00000).
- **Reading.** Delta_C is the Laplace-Beltrami operator of the flat Fisher-Rao metric C^{-1} on the mean fibre. The
  field localization is Brownian motion in that metric, and delta_l is F1's reduced defect written in the layer's own
  Fisher-normal coordinates. It is a commutative heat flow in a non-Euclidean metric. It is not a noncommutative
  object.
- **The scale direction is not empty when S != 0.** At the input, tr d_Sigma E = e/2 was free. At a cut with frozen
  higher tables, <d_Sigma E, C> carries Lambda_S E: the response to rescaling kappa_3 and kappa_4 against the
  covariance. For a non-homogeneous production chain (constant and degree-0 features in beta, clips), <d_Sigma E, C>
  must be computed directly: it is one directional derivative along C itself.

---

## 3. Which pairing governs: KMS embedding, the commutant, and the optimal control

**Proposition 4 (invariants of a matrix defect relative to a covariance; proved).**
1. **Every base rate depends on D only through the whitened defect D~ = C^{1/2} D C^{1/2}.** Its spectrum equals
   spec(D C) = spec(C D). The three are similar through X -> C^{+-1/2} X C^{-+1/2}, a modular-type conjugation.
2. **The Kubo-Ando embeddings collapse on the commutant.** For every symmetric operator mean m_f and every Q in the
   commutant {C}':

       < Q, m_f(L_C, R_C)(D) >_HS = tr(Q C D).

   The geometric (KMS) embedding m_geo(L_C, R_C) D = C^{1/2} D C^{1/2}, the arithmetic (CD + DC)/2 and the harmonic
   one all agree on spectral controls. They differ only on the off-diagonal blocks of D in C's eigenbasis.
3. **Gauge covariance selects the geometric embedding.**
   - Under a linear change of layer coordinates z -> A z, the pair transforms as D -> A^{-T} D A^{-1} and
     C -> A C A^T.
   - The geometric embedding then transforms by an orthogonal conjugation: D~ -> O D~ O^T, with O from the polar
     decomposition A C^{1/2} = (A C A^T)^{1/2} O.
   - Among the embeddings m_f(L_C, R_C), only the geometric one is covariant under the ReLU diagonal gauge
     A = diag(a) > 0.
   - This is F4 Theorem 2.1 transported to the defect.

*Proof.*
1. A similarity preserves the spectrum.
2. In the eigenbasis of C, m_f(L_C, R_C)(D)_jk = m_f(lambda_j, lambda_k) D_jk. A Q in {C}' is block-diagonal on the
   eigenspaces, and m_f(lambda, lambda) = lambda.
3. The polar identity gives (A C A^T)^{1/2} = A C^{1/2} O^T, which is symmetric; substitute it. For uniqueness:
   - with C and A diagonal, covariance requires m_f(a_j^2 lambda_j, a_k^2 lambda_k) = a_j a_k m_f(lambda_j, lambda_k)
     for all a > 0;
   - so m_f(x, y) = sqrt(xy), the F4 argument.

QED.

**Theorem 5 (spectral localizations see only the commutant part of the defect; proved).**
- **Setting.** Suppose K_t K_t^T is in {C}' for all t: for instance a function of C, or any control isotropic within
  the eigenspaces of C.
- **Claim 1.** Then Sigma_t is in {C}' for all t, and the Duhamel integrand of Theorem 2 is

      < E_C(D(m_t, Sigma_t)), Q_t >,    E_C(X) := sum_lambda P_lambda X P_lambda,

  the trace-preserving conditional expectation onto the commutant.
- **Claim 2.** The error, which every control integrates to the same value, is therefore a functional of the
  C-diagonal part of the defect field along any spectral path.
- **Claim 3.** The non-commutative part D - E_C(D) is invisible to every spectral localization. Along a non-spectral
  path it does enter the rate, and since the integral is path-independent its contribution must be compensated
  elsewhere along that path. It is a path gauge, not error signal.

*Proof.*
- {C}' is a *-algebra closed under inversion. So Sigma_t^{-1} = C^{-1} + int K K^T is in {C}', hence Sigma_t is,
  hence Q_t = Sigma_t K K^T Sigma_t is.
- For Q in {C}', tr(D Q) = tr(E_C(D) Q).

QED.

**Computational corollary.**
- A spectral detector needs only the eigen-directional defects d_jj = u_j^T D u_j: n numbers per output, four chain
  evaluations each.
- Truncation to the top k modes leaves the tail sum_{j>k} q(lambda_j) d_jj.
- If d_jj is flat with incoherent sign across modes, the explained share of the truncated detector tracks the
  lambda^2-share of the top k modes. Section 8 uses this for the prediction at n = 1024.

**Proposition 6 (optimal detection control; proved for one output, measured for the family).**
- **One output.** Among controls with unit information rate (1/2) tr(K~ K~^T) = 1/2, the base rate |tr(D~ K~ K~^T)| is
  maximal for rank-one K~ along the top-|eigenvalue| eigenvector of D~, with value ||D~||_op. That eigenvector is a C
  mode only if D~ commutes with C; at n = 64 the per-neuron optima are not C modes (the off-diagonal share is 0.9).
- **A shared control** is what a runtime detector can use. Among spectral controls the field control q = 1 is
  distinguished: it is the only one whose path stays on the scaling ray of the base law (Sigma_t proportional to C).
  So by Corollary 3 the whole Green integral reduces to the one-parameter OU of the ray. That is the structural reason
  a single per-neuron profile, and hence a single slope, can hold.
  - Every other spectral control (q(lambda) = lambda^{p-1}, p != 1) changes the shape of the law along the path. Modes
    then localize at different speeds, and the per-neuron profiles mix scales.

**Measured** (`n2out/e1_n64L8_a.txt`, `e1_n64L16.txt`; GC suffix; truth for the Gaussian surrogate N(mu_l, C_l) by
2^22 antithetic Monte Carlo samples, standard error <= 1.5e-4 against a suffix error of rms 6e-4 to 8e-3). Each entry
is the explained share of the per-neuron suffix error by a one-slope fit through the origin; the demeaned shares are
within 0.03 of these. The field column also gives its correlation with the *total* error: 0.006-0.55 explained, see
`e1_*.txt`.

| net, cut (suffix layers) | suffix share of total MSE | field tr(D C) (slope) | geo C#C_1 | input tr(D C_1), chain / true | tr(C D C) | tr(D C^1/2) | iso tr(D) | top-4 / top-8 / top-16 | F1 input delta vs suffix error |
|---|---|---|---|---|---|---|---|---|---|
| n64 L8 s0, 6 (2) | 0.01 | **0.933** (-0.55) | 0.931 | 0.919 / 0.911 | 0.835 | 0.791 | 0.326 | 0.82 / 0.92 / 0.94 | 0.02 |
| n64 L8 s0, 5 (3) | 0.05 | **0.948** (-0.51) | 0.947 | 0.939 / 0.928 | 0.859 | 0.889 | 0.586 | 0.85 / 0.92 / 0.95 | 0.24 |
| n64 L8 s0, 4 (4) | 0.08 | **0.900** (-0.59) | 0.905 | 0.895 / 0.866 | 0.795 | 0.744 | 0.458 | 0.82 / 0.90 / 0.90 | 0.16 |
| n64 L8 s1, 6 (2) | 0.01 | **0.965** (-0.53) | 0.951 | 0.929 / 0.934 | 0.826 | 0.875 | 0.452 | 0.90 / 0.95 / 0.97 | 0.04 |
| n64 L8 s1, 5 (3) | 0.01 | **0.900** (-0.52) | 0.867 | 0.833 / 0.836 | 0.667 | 0.781 | 0.278 | 0.78 / 0.86 / 0.89 | 0.08 |
| n64 L16 s0, 14 (2) | 0.01 | **0.949** (-0.55) | 0.886 | 0.826 / 0.859 | 0.850 | 0.739 | 0.201 | 0.93 / 0.94 / 0.95 | 0.02 |
| n64 L16 s0, 13 (3) | 0.02 | **0.882** (-0.52) | 0.820 | 0.758 / 0.778 | 0.850 | 0.772 | 0.227 | 0.84 / 0.88 / 0.89 | 0.02 |
| n64 L16 s0, 12 (4) | 0.03 | **0.884** (-0.56) | 0.841 | 0.795 / 0.828 | 0.822 | 0.790 | 0.393 | 0.87 / 0.87 / 0.88 | 0.11 |
| n64 L16 s0, 10 (6) | 0.07 | **0.851** (-0.55) | 0.814 | 0.783 / 0.738 | 0.783 | 0.766 | 0.382 | 0.79 / 0.83 / 0.85 | 0.13 |
| n64 L16 s1, 14 (2) | 0.00 | **0.916** (-0.58) | 0.843 | 0.736 / 0.668 | 0.794 | 0.774 | 0.113 | 0.84 / 0.92 / 0.92 | 0.01 |
| n64 L16 s1, 13 (3) | 0.00 | **0.848** (-0.56) | 0.800 | 0.728 / 0.723 | 0.628 | 0.722 | 0.176 | 0.83 / 0.84 / 0.85 | 0.05 |
| n64 L16 s1, 12 (4) | 0.01 | **0.859** (-0.52) | 0.825 | 0.777 / 0.783 | 0.661 | 0.763 | 0.324 | 0.80 / 0.86 / 0.86 | 0.08 |
| n64 L16 s1, 10 (6) | 0.03 | **0.919** (-0.54) | 0.876 | 0.821 / 0.781 | 0.722 | 0.860 | 0.566 | 0.88 / 0.91 / 0.92 | 0.14 |
| n64 L16 s2, 14 (2) | 0.01 | **0.951** (-0.53) | 0.942 | 0.937 / 0.883 | 0.919 | 0.842 | 0.165 | 0.96 / 0.95 / 0.95 | 0.12 |
| n64 L16 s2, 13 (3) | 0.03 | **0.910** (-0.52) | 0.889 | 0.869 / 0.846 | 0.858 | 0.699 | 0.148 | 0.90 / 0.91 / 0.92 | 0.05 |
| n64 L16 s2, 12 (4) | 0.06 | **0.835** (-0.57) | 0.828 | 0.824 / 0.794 | 0.723 | 0.610 | 0.120 | 0.82 / 0.83 / 0.84 | 0.01 |
| n64 L16 s2, 10 (6) | 0.18 | **0.773** (-0.55) | 0.772 | 0.758 / 0.756 | 0.763 | 0.709 | 0.418 | 0.72 / 0.77 / 0.77 | 0.06 |
| n64 L8 s1, 4 (4) | 0.05 | **0.940** (-0.54) | 0.929 | 0.908 / 0.926 | 0.802 | 0.887 | 0.570 | 0.80 / 0.93 / 0.94 | 0.36 |

Diagnostics on the same cuts:
- **The diagonal is flat.** The rms of d_jj over the C-mode blocks (top 4, 4-16, 16-32, 32-64) is constant within
  +-25%. Example (n64 L8 s0, cut 5): 6.4, 6.9, 6.8 and 6.4e-4.
- **The whitened defect is non-commutative.** The median off-diagonal share of |D~|_F in C's eigenbasis is
  **0.86-0.93**. The median |[D~, Gamma]| / (2 |D~| |Gamma|) is 0.23-0.27.
- **Probing noise.** The median |D~|_F / |tr D~| is 1.7-2.8. A random C-shaped probe v ~ N(0, C) therefore has a
  relative noise of sqrt(2) times that, 2.4-4.0, while the eigen-directional evaluation has none.

**Reading.**
- **Pairing.** The pairing that governs the suffix error is the field trace tr(D C) = tr D~, the KMS/whitened trace.
  The geometric mean of the two natural covariances, C # C_1 (Kubo-Ando), is second. The first-chaos pairing is third.
  The isotropic one fails.
- **The slope is universal:** -0.545 +- 0.02 (sd) over 18 (net, cut) pairs, so err_suffix ~ -0.27 delta_l.
- **Few modes suffice.** At n = 64, 8 of the 64 modes reproduce the detector, because d_jj is flat and C is
  concentrated.
- **The non-commutative part is irrelevant to detection.** 90% of the defect's Frobenius mass is off-diagonal and
  carries no detector power, as Theorem 5 says it cannot at the base point of a spectral path.

---

## 4. Two projections: why the input sees a deep layer through a keyhole

**Theorem 7 (the input localization pushed to a cut; proved).**
- **Claim 1 (revelation profile).** Under the input Eldan localization (posterior X | F_t ~ N(Y_t/(1+t), I/(1+t))),
  with rho_t^2 = t/(1+t) and C_k the k-th Wiener-chaos part of C_l (C_l = sum_k C_k):

      Cov(E[z_l | F_t]) = sum_{k>=1} rho_t^{2k} C_k,      E Cov(z_l | F_t) = C_l - sum_k rho_t^{2k} C_k.

- **Claim 2 (induced control).** In the Gaussian-surrogate reading, the input localization therefore acts on the layer
  as a time-inhomogeneous Chen-Eldan localization with mean quadratic-variation rate

      Q_l(t) = (1+t)^{-2} sum_k k rho_t^{2(k-1)} C_k,    Q_l(0) = C_1 = L^T L,    int_0^inf Q_l dt = C_l.

- **Claim 3 (non-commuting path).** If [C_1, C_l] != 0, the path Sigma_l(t) = C_l - sum_k rho_t^{2k} C_k is a
  non-commuting family. The layer's posterior covariance rotates from the first-chaos frame (t ~ 0) to the
  full-covariance frame (t -> inf).

*Proof.* E[f | F_t] is the Mehler (Ornstein-Uhlenbeck) image of f at correlation rho_t, so its covariance is the chaos
sum (F4 Theorem 5.1). Differentiate in t, using d rho^2/dt = (1+t)^{-2}. If [Sigma_l(t), Sigma_l(t')] = 0 for all t and t', then the
coefficient of rho_{t'}^2 at t = 0 gives [C_l, C_1] = 0. QED.

**Proposition 8 (Halmos pair and canonical correlations; proved).**
- **The two projections.** Let H_1 be the first chaos of L^2(gamma_n) and Z_l = span{z_{l,a} - mu_a}.
- **Claim.** In the orthonormal basis w = C^{-1/2}(z_l - mu) of Z_l, the compression P_Z P_{H_1} P_Z |_Z has matrix

      Gamma = C^{-1/2} C_1 C^{-1/2},   0 <= Gamma <= I.

  Its eigenvalues are cos^2 theta_j of the principal angles between Z_l and H_1, the squared canonical correlations
  between X and z_l.
- **Halmos structure.** The pair (P_{H_1}, P_Z) is the direct sum of a commuting part (theta in {0, pi/2}) and a
  generic part. On the generic part both projections act on 2 x 2 blocks as diag(1, 0) and
  [[c^2, cs], [cs, s^2]], which do not commute.
- **The two base rates.** On the Gaussian block, the input and field base rates are tr(D~ Gamma) and tr D~. They agree
  iff Gamma = I on the support of D~, that is, iff z_l lies in the first chaos (a linear network).

*Proof.* With P_Z f = Cov(f, z) C^{-1} (z - mu) and P_{H_1} f = Cov(f, X) X, the compression has
<w_a, P_{H_1} w_b> = (C^{-1/2} Cov(z, X) Cov(X, z) C^{-1/2})_ab = Gamma_ab. Halmos' two-subspace theorem gives the
decomposition. QED.

**Measured** (E1 and `n2code` one-liners; Gamma uses the chain's first chaos):
- **The keyhole.** Gamma's eigenvalues have **maximum 0.39-0.71 and median 0.000-0.03**. In whitened coordinates
  almost every direction of a deep layer is orthogonal to the first chaos, and the input base point sees the layer only
  through a handful of collective directions.
- **Where the first chaos sits.** The first chaos holds 17-46% of tr C_l (C_1 by the chain) or 30-52% (by Monte Carlo).
  It lies inside the top-8 block of C: 89-99% of tr C_1.
- **The induced path rotates.** |[C_1, C]| / (2 |C_1| |C|) = 0.05-0.10. The off-diagonal share of C_1 in C's
  eigenbasis is 0.20-0.37. So the input-induced path is non-commuting but only mildly.

**Theorem 9 (dressing: field tangents are the skeleton of input tangents; proved).**
- **Claim.** For a functional of finitely many moments of the layer law, the base-point Dynkin rate of the field
  localization (F) is obtained from that of the input localization (I) by replacing every tangent Gram
  <grad_x f, grad_x g> = Cov(f, X) . Cov(X, g) by Cov(f, z) C^{-1} Cov(z, g). In particular:

      <grad mu_a, grad mu_b> = (C_1)_ab                          ->  C_ab
      <grad mu_a, grad C_bc> = kappa3(P_1 z_a, z_b, z_c)          ->  kappa3_abc
      <grad C_ab, grad C_cd> = <P_1(z_a z_b), P_1(z_c z_d)>        ->  kappa3_ab. C^{-1} kappa3_cd.

- **The GC readout.** For a Gaussian readout h_G, F1 Proposition 7's bracket -(w3/2)<grad mu, grad v> -
  (w4/8)|grad v|^2 becomes (w3/2) kappa3_iii + (w4/8) D21_i C^{-1} D21_i^T:
  - the first is the first-order Edgeworth term;
  - the second is the Schur form of the (2,1) slice, the dressed "infinitesimal-leaf" kappa_4 of F5.

*Proof.* P_{H_1} f = Cov(f, X) X and P_Z f = Cov(f, z) C^{-1} (z - mu). The (I) and (F) quadratic variations at t = 0
are the Gram matrices of these projections (Kushner-Stratonovich for a linear Gaussian observation), and
Cov(z_a z_b, z_c) = kappa3_abc + mu_a C_bc + mu_b C_ac. QED.

**Reading.**
- This is note XLI's bare-to-dressed map, now derived: the field localization is the input localization with every
  line resummed into a full propagator.
- It is also why the bare path class was 4-8x too small at depth, where Gamma ~ 0 for most directions.
- Caution from note XLI. The symmetrized-slice Schur form 3 diag(D21 C^{-1} D21^T), which is the dressed pair term
  here, was useless as a kappa_4-diagonal predictor. Only the source-resolved hub Y carried the kappa_4 error. The
  dressed pair is a defect term, not a kappa_4 estimate.

**Measured** (E3, `n2out/e3_n64L8s0.txt`; true kappa_3 at the cut from 2^21 Monte Carlo inputs; GC suffix).
- **What is computed.** The field generator of the true law at the cut is

      A E = -tr(D C) + sum_a d_mu_a <d_C E, kappa3[a]> + (1/2) sum_j d^2 E[V_j, V_j],   V_j = kappa3 x_3 C^{-1/2} e_j.

  It is exact, because the GC suffix reads only (mu, C) of the posterior.
- **The target.** err_A = truth - GC_suffix(true mu*, C*) is the suffix error from the exact state; it includes the
  layer's non-Gaussianity.

| cut (suffix) | err_A share of total MSE / corr | dressed A E explains err_A / total | Gaussian block alone (err_A / total) | dressed cross / pair (err_A) | F1's bare suffix part Lambda_l explains total (E1b, same net) |
|---|---|---|---|---|---|
| 5 (3) | 0.40 / 0.78 | **0.934 / 0.661** | 0.639 / 0.418 | 0.617 / 0.313 | 0.519 |
| 4 (4) | 0.67 / 0.89 | **0.928 / 0.830** | 0.544 / 0.423 | 0.658 / 0.072 | 0.680 |
| 6 (2) | 0.44 / 0.65 | **0.957 / 0.522** | 0.100 / 0.072 | 0.857 / 0.026 | 0.456 |

- **Reading.** Dressing beats bare by 0.07-0.15 of the total at the same cut. But it needs the **true** next-order
  cumulant at the cut. Fed with the chain's own (absent) kappa_3, GC's dressed generator collapses to the Gaussian block.
  For a kappa_3-carrying chain the dressed generator needs kappa_4 with a field leg, which is the next-order wall.

---

## 5. How the late defect relates to the full delta (the cut decomposition)

**Theorem 10 (cut decomposition of the input defect; proved).**
- **Setting.** Let s_l(m, Sigma) be the prefix map. For each carried table X put D(X) := tr H X - R(X), where R is
  the right side of the cumulant hierarchy (F1 (5.1)): R(mu) = 0, R(C) = sum_j d_j mu d_j mu^T = C_1, and
  R(kappa3) = 3 sym sum_j d_j mu (x) d_j C.
- **Claim.**

      tr D_in = d_mu E_suf . D(mu_l) + <d_C E_suf, D(C_l)> + <d_S E_suf, D(S_l)>  +  Lambda_l,
      Lambda_l = < D_l^G , C_1 >  +  X_l,

  where:
  - D_l^G = d_C E_suf - (1/2) Hess_mu E_suf is the S-frozen defect at the cut;
  - X_l collects the cross and pair terms of the covariance and source tangents (the kappa_3 with an input leg).
- **Reading.** The first three terms are the prefix's defects transported by the suffix. Lambda_l is the suffix-local
  part of F1's telescoping (F1 Proposition 8), with input tangents: the "A-late" object of the synthesis.

*Proof.* Chain rule, H(E_suf o s_l) = DE_suf . H s_l - (1/2) D^2 E_suf[grad s_l, grad s_l], then H s_l = R(s_l) +
D(s_l). QED.

**Measured, GC** (`n2out/e1b_n64L8.txt`; four networks n64 L8; tr D_in explains 0.92-0.95 of the total error).
Columns are explained shares of the **total** error.

| cut (suffix) | upstream (transported prefix defects) | Lambda_l (suffix-local, input tangents) | of which Gaussian block <D_l^G, C_1> | of which cross + pair X_l |
|---|---|---|---|---|
| 6 (2) | 0.73-0.91 | 0.27-0.49 | 0.05-0.13 | 0.23-0.45 |
| 5 (3) | 0.74-0.87 | 0.50-0.65 | 0.17-0.34 | 0.33-0.63 |
| 4 (4) | 0.60-0.84 | 0.68-0.82 | 0.27-0.55 | 0.49-0.72 |
| 2 (6) | 0.32-0.84 | 0.82-0.94 | 0.43-0.88 | 0.35-0.81 |

**Measured, K3** (`n2out/e2b_k3cut_n48L8.txt`, `e2_k3_*.txt`; n = 48, L = 8; tr D_in explains 0.74).

| net, cut (suffix) | upstream | Lambda_l (input tangents) | **S-frozen field tr(D_l C_l)** | joint (Lambda_l, field) |
|---|---|---|---|---|
| s0, 5 (3) | 0.19 | 0.32 | **0.48** | 0.49 |
| s0, 4 (4) | 0.13 | 0.31 | **0.64** | 0.64 |
| s1, 6 (2) | | | **0.00** | |
| s1, 5 (3) | 0.60 | 0.05 | **0.04** | 0.14 |
| s1, 4 (4) | 0.52 | 0.18 | **0.00** | 0.18 |
| s2, 5 (3) | | | **0.09** (top-16) | |
| s2, 4 (4) | | | **0.07** (top-16) | |
| s2, 3 (5) | | | **0.44** (top-16) | |
| s3, 5 (3) | | | **0.08** (top-16) | |
| s3, 4 (4) | | | **0.14** (top-16) | |
| s3, 3 (5) | | | **0.69** (top-16) | |

On s1 the full input delta explains 0.86 of the total; on s0, 0.74.

**Proposition 11 (blindness of a late cut; proved, trivial but decisive).**
- **Claim.** Every defect evaluated at or after the cut from the chain's state, whether S-frozen, field, dressed with
  the chain's own tables, or rolling, is a function of (s_l, suffix) only. Two chains that agree on the suffix and on
  s_l have identical late defects, whatever their prefix errors.
- **Consequence.** A late defect can see error made before the cut only through correlation, or through the next-order
  information the suffix's own tangents encode via the tilt identity at the cut.

**Reading of sections 3-5 together.**
1. **The late field defect is an almost perfect detector of the error its suffix creates, and that error is
   usually a small part of the total.**
   - For GC it is 0-8% of the total (2-4 layer suffixes) and 3-18% at 6 layers.
   - For K3 the S-frozen defect detects 0-14% of the total at 3-4 layer suffixes on three of four networks, and
     48-64% on s0.
   - At a 5-of-8 suffix it detects 44-69%: the creation horizon is about 5 layers.
2. **The difference is the type of omission.**
   - GC's omission (kappa_3) is born at every layer and accumulates: old sources keep being read. So error created
     early and read late dominates. Only input tangents see it: F1's Lambda_l through its cross and pair terms, the
     kappa_3 with an input leg.
   - K3's omission (kappa_4) is born locally but also transported. Whether a late cut sees it depends on the
     network: one of four has a late-created error.
3. **Injected late is not created late** (`n2out/e5_injection.txt`, the DATA_FINDINGS decomposition
   e_j = Phi(alpha_j) o (e_(j-1) W_j) + inj_j, with each inj_j pushed to the output by the chain's mean gates).

   | chain | share of the final error injected in the last half of the layers | in the last 2 layers | late field defect detects (total) |
   |---|---|---|---|
   | GC n64 L8 (4 nets) | 0.25-0.48 | 0.12-0.35 | 0.02-0.55 (3-4 layer suffix) |
   | GC n64 L16 (3 nets) | 0.15-0.74 | 0.04-0.24 | 0.006-0.19 |
   | K3 n48 L8 (4 nets) | 0.67-0.79 | 0.40-0.57 | 0.00-0.64, median about 0.08 (3-4 layer suffix); 0.44-0.69 at 5 layers |

   - The production chain's "two thirds in layers 12-15" (DATA_FINDINGS) is the same phenomenon. An injection at
     layer j is the readout error caused by table errors at j, most of which were created by the maps of earlier
     layers.
   - A localization at the cut sees creation, not injection.
4. **The three detectors are complementary.**
   - The input delta is dominated by the upstream terms.
   - The bare suffix part Lambda_l is dominated by cross and pair terms.
   - The S-frozen field defect is the dressed Gaussian block.

---

## 6. The analytic two-layer field defect (row-locality of mean tangents)

**Theorem 12 (two-readout suffix; proved; checked).**
- **Setting.** Take a cut l whose suffix is: readout at l, transport by W = W_(l+1), Gaussian (or frozen-Edgeworth)
  readout at l+1.
- **Claim.**

      tr(D_suf C_l) = -(E f'''(z'_i)/2) g_i - (E f''''(z'_i)/8) q_i,
      g_i = grad mu'_i^T C_l grad v'_i,    q_i = grad v'_i^T C_l grad v'_i,

  where grad = d/d mu_l, mu' = m W and v' = diag(W^T C_y W). The tangents are:

      grad mu'_{ci}  = W_ci Phi(alpha_c)
      grad v'_{ci}   = 2 W_ci [ sum_{k>=1} (A_{k+1,c}/sigma_c) (M_k (A_k o W))_{ci} / k! ] + 2 W_ci^2 A_{0,c} (1 - Phi_c),

  with M_k = R^{o k} (diagonal zeroed) and A_k = sigma^k E f^(k) the Hermite coefficients of the readout at l.
- **Reading.**
  - g_i = Cov(z'_i, z_l) C_l^{-1} kappa3(z_l, z'_i, z'_i): the kappa_3 born at layer l, projected on the layer-l
    field.
  - q_i is the Schur form of the (2,1) cumulant between layer l+1 and the layer-l field.
  - Together these are F1 Proposition 7's bracket with the input metric replaced by the field metric C_l. It is the
    dressing of Theorem 9 one layer back.

*Proof.*
- Write E = h(mu', v'). Then

      tr(D C) = h_mu [<d_C mu', C> - (1/2) Delta_C mu'] + h_v [<d_C v', C> - (1/2) Delta_C v']
                - (1/2)[h_mumu |grad mu'|_C^2 + 2 h_muv <grad mu', grad v'>_C + h_vv |grad v'|_C^2].

- The first bracket is W^T applied to per-neuron heat defects of E relu, which are zero.
- The second bracket is diag(W^T [(d_C - (1/2) Delta_C) C_y] W). For the exact Gaussian covariance map the hierarchy
  (F1 Thm 4, k = 2) gives Phi_a C_ab Phi_b, so the bracket is |grad mu'|_C^2.
- Since h_v = h_mumu/2 (Lemma 1), it cancels the first term of the last bracket.
- With h_muv = E f'''/2 and h_vv = E f''''/4, the formula follows. The tangent formulas come from
  d_mu A_k = A_{k+1}/sigma and d_mu Var relu = 2 A_0 (1 - Phi).

QED.

**Checked** (`n2out/t_analytic2.txt`). Against the central-difference field defect, three nets give correlation
1.000000, with relative differences 6.6e-5, 9.8e-5 and 6.3e-4.
- The g-term carries the norm (share 1.01-1.06); the q-term is 0.16-0.21.
- **Mehler order in the tangent.** Truncating the tangent at order Kd = 1, 2, 3, 4 gives relative errors of 4.7-6.6%,
  0.6-0.9%, 0.16-0.24% and 0.05-0.11%.
- **Metric truncated to the top-r modes of C_l.** Correlation with the full metric is 0.93-0.99 (r = 4), 0.99 (r = 8),
  1.00 (r = 16). The explained share of the suffix error is 0.82-0.93 (r = 4) and 0.92-0.94 (r = 8), against
  0.93-0.95 for the full metric.

**The cost lesson (proved).**
- d/d mu_c of any table produced by a ReLU layer touches only the coefficients of neuron c: row and column c.
- So the full mean-Jacobian of a next-layer table costs **one dense product per coefficient family**, for all n
  directions at once. It is the same price as one directional derivative.
- The collective subspace then only shrinks the final metric contraction, from n^3 to n^2 k.
- **"k directions instead of n" is the right economy for probing** (finite differences, section 3: k ~ 8 of 64, about
  32-64 of 1024 predicted in section 8). It is the wrong economy for analytic tangents, where it saves nothing.

**What changes for the production chain** (sketch). Its readout at l+1 = 15 reads (mu', v', D3', g4'), and D3' depends
on (mu_14, C_14) through the young births at 14 and the gate transport of every older source.
- **The GC g-term is not the production defect.** It is exactly the kappa_3 born at l and read at l+1, which the
  production chain already carries through the newborn source. Using the GC formula alone would count the chain's own
  correction twice.
- **What remains is one level up:**
  - h_{mu,D3} [<grad mu', grad D3'>_C - R(D3')];
  - h_{v,D3} <grad v', grad D3'>_C;
  - the kappa_4-closure terms;
  - the non-Edgeworth readout parts (Lemma 1).
- **Its tangents are still row-local:**
  - grad D3' is one product per fused dslice family (the gate derivative w1'_c on row c of every transported leg) plus
    the newborn's row-local derivative;
  - grad v' adds one product per covariance-term table (C_off, D21, the (2,2) programs).
- **Estimate:** 10-20 units. **Value:** see section 5. A two-layer suffix sees almost nothing of a K3-type error
  (s1: 0.00), so this cheap version is not a candidate on its own.

---

## 7. The rolling field defect, injected in place (small-n analogue of the component)

**Definition.**
- For each late layer c + 1, take the one-step map (mu_c, C_c; S_c frozen) -> m_(c+1). This is readout c, transport,
  readout c+1 (post-activation means at c+1).
- Its field defect r_(c+1) = tr(D C_c) is a vector over the neurons of layer c+1.
- Inject it in place: m_(c+1) <- m_(c+1) + a r_(c+1). The chain's own transport then carries the correction to the
  output.
- This is the "A-late in place" design of the synthesis, with the input tangents replaced by the field tangents of
  each layer. It is the telescoping re-metrized layer by layer with each layer's own Fisher metric C_c.

**Measured** (`n2out/e4_rolling_k3_n48L8.txt`; K3 chain, n = 48, L = 8, four networks; r computed by central
differences in the top-16 of 48 eigen-directions of C_c; the slope a is fitted leave-one-network-out):
- the per-layer correlation between r_(c+1) and the chain's per-layer mean error is weak: +0.05 to -0.66 over four
  networks and layers 3-7, with median about -0.2;
- the injection results are in the table below.

| injected at layers | s0 | s1 | s2 | s3 | mean | geometric mean | fitted a* (leave-one-out) |
|---|---|---|---|---|---|---|---|
| last 1 (7) | 0.883 | 1.206 | 0.933 | 0.956 | 0.995 | 0.988 | -0.5 to -1.3 |
| last 2 (6-7) | 0.827 | 1.066 | 0.969 | 0.876 | 0.935 | 0.930 | -0.6 to -0.8 |
| last 4 (4-7) | 0.654 | 1.085 | 1.010 | 0.786 | **0.884** | **0.866** | -0.6 to -0.9 |
| layers 3-7 | 0.706 | 0.937 | 1.126 | 0.753 | 0.880 | 0.865 | -0.6 to -0.9 |

The entries are held-out final-layer MSE ratios: the slope a is fitted on the other three networks, on a grid of
step 0.1.

**Reading.**
- The one-step field defect, injected in place and carried by the chain, removes **12-13% of the K3 MSE on average**
  when applied over the last half of the network.
- The spread is wide: -35%, -21%, +1% and +9%.
- The slope is stable at about -0.6 to -0.9, so a truth-free class coefficient is plausible.
- A single layer is worth nothing (0.99).
- This is the best runtime-shaped result of the frame. It is still modest, because each one-step defect sees only the
  kappa_4 born at c and read at c+1, not the transported part.

**Pricing at n = 1024.**
- The production one-step defect is estimated at 10-20 units per layer (section 6), so the last 4 layers cost
  40-80 units: C'/C = 1.19-1.38.
- Break-even then needs an MSE ratio of 0.84-0.72. The K3 analogue (0.87-0.88, and over half the depth rather than a
  quarter) falls short: **adjusted x(1.03-1.22), negative.**
- It pays only if the per-layer cost can be brought to <= 5 units, the GC formula's cost class. That needs the D3 and
  g4 tangents of the production layer to be dropped or folded into existing families; none of this is shown.

---

## 8. Estimator components, costs, and the prediction

Costs are in units of 2n^3 at n = 1024. The present bill is 208 units: C/B 0.203, raw 1.55e-8, adjusted 3.14e-9.
A late production layer costs about 13-15 units.

| component | what it computes | cost | detected share (small-n analogue) | status |
|---|---|---|---|---|
| L0. production calibration defect (Lemma 1) | H(beta . feats) + V30 term, closed form | ~0 (n per layer) | none established | diagnostic only |
| L1. analytic two-layer field defect at cut 14 (Thm 12, production version) | readout-15 local defect in the C_14 metric | +10-20 | GC: 0-1% of total; K3: 0.00 (s1) | **negative** |
| L2. rolling one-step field defect, last 4 layers, in place | sum of L1-type defects at 12-15 | +40-80 | K3: held-out MSE x0.88 (spread 0.65-1.09) | **marginal-negative** (adjusted x1.03-1.22); pays only at <= 5 units/layer |
| L3. multi-layer late field defect, cut 12 (4-layer suffix), by finite differences in the top-k modes | tr(D_12 C_12), k = 8-32 | k x 4 x ~50 = 1600-6400 | K3: 0.48-0.64 (s0), 0.04 (s1); GC: 0.02-0.55 | **infeasible at runtime** |
| L3'. L3 by analytic tangents ("tangent sources") | see below | +20-40 (sketch, unvalidated) | as L3 if exact | **gated** |

**L3' (sketch).** The multi-layer field defect from cut l needs the C_l-metric second variation of the suffix, through
3-4 layers. Three structural facts make it affordable in principle:
1. **Row-locality at the cut.** The mean tangent of every table born at l is row-local: d/d mu_c touches row and
   column c.
2. **Rank-2 transport of covariance tangents at first Mehler order.**
   - The first Mehler term is A_1 A_1^T o R. A row-local tangent of C_y is diag(e_c) X + X diag(e_c).
   - Hadamard multiplication by a rank-one kernel keeps it rank-2 per direction. So W^T (.) W keeps the covariance
     tangent at rank 2 per direction, and all n directions cost n^3 per layer.
   - The higher Mehler orders break the rank. In the two-layer check they are a 5-7% correction (Kd = 1 against full,
     section 6).
3. **Tangent sources.** The mu_l-derivative of the young births at l is row-local and is transported linearly, like
   the sources themselves. Carrying "tangent legs" for the 2-3 newest births through the last layers costs about
   2 legs x 3 sources x 3 layers = 18 transport products, plus about 6 hub-type contractions with the C_l metric.

The risks:
- the gate dependence of the transport on mu (a second-order effect);
- the higher Mehler orders;
- the hierarchy residues of the production's non-Gaussian covariance terms.

None of these is measured. The leading-order tangents of the D21/D3 dslices are the same objects as V56's hub
families, so part of the work may share their products.

**Prediction.**
- **Top-k capture at n = 1024 (from flat d_jj and the measured spectrum).**
  - The layer-15 spectrum is fitted to the official values: PR/n = 0.052, top-16 share 0.42, top-128 share 0.89
    (lambda_j proportional to j^{-0.55} e^{-j/96}).
  - That gives lambda^2-shares of 0.71 (k = 8), 0.83 (16), 0.92 (32) and 0.98 (64).
  - At n = 64, the ratio (explained by top-k) / (explained by full) tracked the lambda^2-share: 0.87-0.93 at
    lambda-share 0.55, 0.95-1.0 at 0.75-0.8.
  - **Predicted: k = 32 recovers about 90% of the full field detector, k = 64 about 97%.**
- **Why the error's enrichment in collective modes does not bound this.**
  - The measured enrichment (top-16/64/128 hold 8-15% / 27-43% / 46-63% of the pre-activation mean error) is a
    statement about where the **error vector** lies.
  - The truncation here is on the **localization** side. Each of the k directions produces a defect for every output
    neuron, so a top-k field detector is not confined to correcting inside the top-k subspace.
  - The two facts are consistent: the error lives where C is large, and the field metric weights directions by
    lambda. But the enrichment is not a ceiling.
- **Visibility at n = 1024.** This is the open quantity. The production chain's omissions are a mixture:
  - local: the kappa_4 closure, the joint-gate triangle born from two current births;
  - accumulated: the old-source tier-2 compression, transported table errors.
  
  Its late-suffix visibility should be like the small-n chains, where injection is late and creation is early.
  Predicted one-slope explained share of the total final error:
  - **X_late = 0.08 [0, 0.35]** at cut 12 (4-layer suffix);
  - **0.03 [0, 0.10]** at cut 14 (2-layer);
  - slope on tr(D C) between -0.3 and -0.8.

  The wide upper band reflects the network-to-network spread at small n: for K3, 0.04 against 0.48 at the same cut.
- **Adjusted outcome.** The adjusted ratio is (1 - X)(C'/C).
  - L3' at +30 units (C'/C = 1.14) pays iff X_late >= 0.13, and gives -20% adjusted only at X_late = 0.30.
  - L1 at +15 units pays iff X >= 0.07, which is above its predicted band.
  - **Expected value: slightly negative.** P(X_late >= 0.13) is about 0.2, with a right tail if the production
    error is late-created on most networks.

---

## 9. The minimal decisive experiment (offline; the lead runs it on AWS)

**Purpose.**
- Measure X_late, the share of the production chain's final error visible to a late Gaussian-part field defect. This
  decides L1/L2/L3' outright.
- Settle "local or accumulated" for the production error, which bounds every late-layer design including the
  synthesis's A-late.

**Inputs.**
- `estimator_final_v56.py`, patched as follows:
  - restartable at a pre-activation cut l from (mu, C, S_l);
  - S_l frozen, together with every mask, count, basis, rank selection and clip, at its base-run value (F1 referee
    R4.1);
  - run in float64.
- Official networks 0-7, `truth_off{k}.npz` (per-layer truth means) and `chaindump_{k}.npz` (mu_l, C_l).

**Validation first.**
- The float64 restart at each cut with the unperturbed (mu_l, C_l) must reproduce the float32 production output to
  <= 1e-6 relative.
- A no-op perturbation (h -> 0) must be smooth: compare h and h/2 on 8 directions.

**Runs.**
- For each network and cut l in {12, 13, 14}:
  1. Compute the top 64 eigenvectors u_j of C_l.
  2. For j = 1..64, evaluate the suffix at (mu +- h_m u_j, C) and (mu, C +- h_s u_j u_j^T):
     - h_m = 0.02 sqrt(mean var_l) along the unit vector u_j;
     - h_s = 2e-3 mean var_l;
     - 8 directions are repeated at h/2 (Richardson check) and 8 directions unfrozen (mask sensitivity).
  3. Form d_j = [E(C + h_s P_j) - E(C - h_s P_j)]/(2 h_s) - [E(mu + h_m u_j) + E(mu - h_m u_j) - 2E]/(2 h_m^2), and
     the partial traces T_k = sum_{j<=k} lambda_j d_j for k = 8, 16, 32, 64.
- **Rolling variant (L2).** The same with the one-step map (mu_c, C_c; S_c) -> m_(c+1) for c = 11..14, with 32
  directions each.

**Compute.**
- Per network: 3 cuts x (4 x 64 + 1 + 16) = 819 suffix evaluations.
- A late full layer costs about 15 units and the trimmed last layer about 5. So a suffix from 12, 13 or 14 costs about
  0.24, 0.17 or 0.10 of a 208-unit chain, about 140 chain-equivalents in all.
- The rolling variant adds 4 x 129 evaluations of a one-step map (about 0.07 of a chain each), about 36.
- That is about 175 float64 chain-equivalents per network, 1400 for 8 networks.
- At about 20 s per float64 predict, this is about 8 core-hours: under 15 minutes on a 96-core instance, including
  restarts.

**Outputs.** Per network and cut:
- the explained share of err_final by T_k (one slope through the origin): in-sample, and held out with the slope
  fitted on the other 7 networks;
- the slope;
- the T_k-versus-k capture curve, and the flatness of rms d_j across mode blocks;
- the h/2 and unfrozen controls.
- For the rolling variant, the per-layer correlation of r_(c+1) with the DATA_FINDINGS injection inj_(c+1).
- Also corr(T_64, delta_hat_in) against the pending F1 X* run on the same networks: complementarity.

**Decision rule.**

| outcome | action |
|---|---|
| held-out X_late(cut 12 or 13, k <= 64) >= 0.25 and slope stable within +-25% across networks | derive and validate L3' (analytic tangent sources) on network 0 against the finite-difference T_64. Build it only if it reproduces >= 80% of T_64's explained share at <= +40 units. Then cold-screen 16 networks and score. |
| 0.10 <= X_late < 0.25 | build nothing now. Record X_late and its per-network spread as the late-visibility measurement for every future late-layer design (A-late included). |
| X_late < 0.10 | the late frame closes as explanation. The production error is accumulated and only input tangents (F1) can see it. |
| independently: rolling X(L2) >= 0.15 held out | cost the production L1 per layer (+10-20 units) and screen L2 |

**Registered predictions** (section 8):
- X_late = 0.08 [0, 0.35] at cut 12, 0.03 [0, 0.10] at cut 14, with a large network-to-network spread;
- capture T_32/T_64 >= 0.9;
- d_j flat within +-30% across mode blocks;
- slope between -0.3 and -0.8;
- corr(T_64, delta_hat_in) < 0.5 (complementary detectors).

---

## 10. What the non-commutativity is, and what it buys (ledger)

| where | object | genuinely non-commutative? | what it buys |
|---|---|---|---|
| pairing D with C | D~ = C^{1/2} D C^{1/2} (KMS / geometric embedding), spec(DC) | the pair (D, C) does not commute: the off-diagonal share of D~ is 0.86-0.93 | **the canonical embedding (unique gauge-covariant one, Prop 4).** But the error detector is a trace, and **spectral controls see only E_C(D)** (Thm 5). The non-commuting 90% is path gauge: invisible to the field detector, carrying none of its power. It matters only for probing variance (|D~|_F/tr D~ ~ 2) and per-neuron optimal directions. |
| Kubo-Ando / Bures / KMS means of covariances | C # C_1 as a control | yes ([C_1, C] = 0.05-0.10) | the geometric mean of the input-induced and field controls is the second-best pairing (0.77-0.95). There is no evidence any mean beats the field trace. |
| input vs field localization | Halmos pair (P_{H_1}, P_{Z_l}); Gamma = C^{-1/2} C_1 C^{-1/2} | **yes, the real one.** Two projections in generic position; at depth Gamma's median eigenvalue is ~0 and its maximum 0.4-0.7 | **the keyhole theorem** (Thm 7, Prop 8): the input base point sees a deep layer only through a few collective directions, and reveals the rest in chaos order along a non-commuting path. This explains why F1's input delta is nearly blind to suffix-created error (1-36%), and why a layer's own field trace is the better late detector (0.77-0.97 against 0.67-0.94). |
| dressing | P_{H_1} -> P_Z in every tangent Gram | consequence of the above | **derives note XLI's bare-to-dressed rule** (skeleton resummation = field localization). With the true next-order cumulant, dressed beats bare by 0.07-0.15 of the total error (E3). |
| Wick/chaos algebra | pair Grams <P_1(z_a z_b), P_1(z_c z_d)> = L^T H_a H_c L | yes (F4, synthesis) | unchanged: the pair (V2) terms are small in GC's suffix (dressed pair 0.01-0.31 of the total), and the n^4 wall is untouched. |

**Bottom line of the ledger.**
- Non-commutativity enters the late-layer defect in exactly one load-bearing place: the relative position of the first
  chaos and the layer's field span. That fixes which localization sees what.
- The KMS/modular machinery (Theorems B and C of stage 8) reduces, on the commutant of C, to a single trace. The
  non-commutative remainder of the defect matrix is invisible to every spectral localization.
- The frame's positive content:
  - the layer delta (Cor 3);
  - the field pairing with its universal slope and flat diagonal;
  - the analytic two-layer formula with row-locality;
  - the dressing derivation;
  - the cut decomposition with its blindness proposition.
- The frame's negative content: a late cut cannot see error made before it. Whether that leaves enough for a runtime
  component is the one measurement of section 9.

---

## 11. Files

| file | content |
|---|---|
| `n2code/n2_lib.py` | restartable GC suffix, Monte Carlo suffix truth, directional and full defect matrices |
| `n2code/n2_exp1.py`, `n2out/e1_n64L8_a.txt`, `e1_n64L16.txt` | pairings, top-k, non-commutativity diagnostics (section 3) |
| `n2code/n2_exp1b.py`, `n2out/e1b_n64L8.txt` | GC cut decomposition (Thm 10) |
| `n2code/n2_dressed.py`, `n2out/e3_n64L8s0.txt` | dressed field generator with the true kappa_3 (Thm 9) |
| `n2code/n2_analytic2.py`, `n2out/t_analytic2.txt` | analytic two-layer field defect (Thm 12) |
| `n2code/n2_k3.py`, `n2code/n2_k3cut.py`, `n2code/n2_k3lean.py`; `n2out/e2_*.txt`, `e2b_*.txt`, `e2c_*.txt` | K3 chain: S-frozen defect, cut decomposition |
| `n2code/n2_rolling.py`, `n2out/e4_rolling_k3_n48L8.txt` | rolling in-place injection (section 7) |

---

## Referee report

Adversarial referee for frame N2. My checks are independent of `n2code/` and live in `loc/ref_n2/`:
- `chk_small.py`: the Lemma 1 closed forms, the Edgeworth sign, and the section 8 spectrum fit;
- `chk_green.py`: Corollary 3 end to end with y_0 != 0, on my own GC two-step suffix, against Monte Carlo truth;
- `chk_rank.py`: the rank claim of Conjecture 13.

I also read the rolling-injection raw output (`n2out/e4_rolling_k3_n48L8.txt`) and the synthesis's A-late design.

### R1. What is correct (verified)

1. **Lemma 1 closed forms.** These are correct. Central differences reproduce the V30 defect (g/4) alpha^2 phi/sigma
   to 7 digits at three points (1.437e-2, 3.367e-2, 1.051e-1), and H(sigma phi) = phi/sigma likewise. I re-derived
   H(Phi(alpha)) = 0 and H(phi(alpha)) = phi/(2 var) by hand.
2. **Theorem 2 and its whitened form.** Correct. Sigma^{-1} = C^{-1} + int K K^T, R = C^{-1/2} Sigma C^{-1/2}, and
   the trace identity check out.
3. **Corollary 3 (S = 0), including the shifted start.** Correct.
   - Setup: n = 3, a random dense C, mu != 0, Mehler-12 GC readout, transport, then a Gaussian readout. Truth comes
     from 2^22 antithetic samples; the Green integral uses 24 Gauss-Legendre nodes in r = (1+tau)^{-1/2} with
     2e4 samples per node.
   - u - E = [ 0.006441, -0.038578, -0.068388] (se <= 1.4e-4).
   - -(1/2) int (1+tau)^{-3/2} E_{N(sqrt(1+tau) y_0, tau I)} delta_l = [0.006427, -0.038571, -0.068331].
   - 2 tr(D C) from the covariance and mean differences equals delta_l(y_0) to 1e-7 relative.
   - The F1 referee checked y_0 = 0 only; the shifted start is now checked too.
4. **The remaining identities are correct**, and all are short:
   - Proposition 4(1, 2);
   - Theorem 5 claims 1-2;
   - Theorem 7 (the Mehler image of the Eldan posterior);
   - Proposition 8 (Gamma is the canonical-correlation matrix; 0 <= Gamma <= I by Bessel);
   - Theorem 9 (the Kushner-Stratonovich innovation Grams). I re-derived the dressed GC bracket
     +(w3/2) kappa3_iii + (w4/8) D21_i C^{-1} D21_i^T as the field generator, signs included;
   - Theorem 10 (the chain rule for H on E_suf o s_l);
   - Theorem 12, whose tangent formula I re-derived, including d Var relu/d mu = 2 A_0 (1 - Phi).
5. **The section 8 spectrum fit reproduces.** lambda_j proportional to j^{-0.55} e^{-j/96} gives PR/n = 0.055,
   top-16 0.42, top-128 0.90, and lambda^2-shares 0.71 / 0.83 / 0.92 / 0.98 at k = 8 / 16 / 32 / 64.
6. **The author's own negative verdict on runtime components is right.** Section R3 shows it is, if anything,
   too generous.

### R2. Mathematical errors and overclaims

1. **Lemma 1(a) has the wrong sign.** The Edgeworth shift acting on u(m, Sigma) is exp(sum_k kappa_k (+d_m)^k / k!),
   not (-d_m)^k: the minus sign belongs to the density, where d_z = -d_m.
   - Check: xi two-point skewed, G = x^3. E(m + xi)^3 = 3.125 = u + 3m kappa_2 + kappa_3, while the (-d) form gives
     0.125.
   - This is harmless for heat-exactness (H commutes either way). But any production use of the formula would flip
     every odd-order term.
2. **The headline "layer delta" holds only for S = 0.** The summary and the structured statement give
   delta_l = e - mu.grad e - Delta_C e = 2 tr(D C). For the production chain (S != 0, frozen) the identity is
   2 tr(D C) = delta_l - Lambda_S E, and the chain also carries non-homogeneous parts (degree-0 and constant
   features in beta, clips).
   - The body says this, but the cheap "radial" form must not be used at production: it would mis-state the defect
     by Lambda_S E, which is of the same order as the signal.
   - The section 9 stencil measures <d_Sigma E, C> directly, so the experiment itself is fine.
3. **For S != 0 there is no error theorem, only an Ito identity.** A point mass carrying kappa_3 is not a law, so the
   boundary term needed to turn Theorem 2 into "error = integrated defect" does not exist. The K3 and production
   detectors are therefore purely empirical (b'(1) with a fitted p).
   - The near-universal slope -0.545 +- 0.02 and the 77-97% are S = 0 (GC) results, and GC is the case the
     production chain is not.
   - The relevant numbers are the K3 ones: 0-64% of the total, median about 0.08.
4. **Theorem 5 claim 3 and the "measured" invisibility are tautologies, not findings.**
   - tr(D C) = sum_j lambda_j d_jj is an identity. "The C-diagonal part alone gives the detector's full power" is
     therefore true of every matrix D and carries no evidence.
   - "Path gauge" is an overstatement. The off-diagonal part of D at the base point is genuine information about E;
     only one particular functional of it, the trace with a commuting Q, ignores it.
   - The non-commutativity figure is not large against a baseline. A random symmetric 64 x 64 matrix has an
     off-diagonal Frobenius share of 0.98 (`ref_n2`, GOE, 3 seeds). D~'s 0.86-0.93 means its diagonal holds
     14-26% of |D~|_F^2, against about 3% at random. D~ is substantially aligned with C's eigenbasis, which is the
     opposite of the reading given.
5. **Proposition 4(3) changes no number in the paper.**
   - For every operator mean, tr m_f(L_C, R_C) D = tr(C D).
   - The geometric "KMS embedding" C^{1/2} D C^{1/2} is the defect written in whitened coordinates
     y = C^{-1/2}(z - mu).
   - Its uniqueness under the diagonal gauge is correct but idle: no computed quantity depends on the choice of
     embedding.
6. **Two smaller overstatements.**
   - **Proposition 8 "agree iff".** The input and field rates agree *for all D~* iff Gamma = I on supp D~. For a given
     indefinite D~, tr(D~(I - Gamma)) = 0 can hold otherwise.
   - **Theorem 9 "is the first-order Edgeworth term".** It is the same tensor slice, but the generator coefficient is
     w3/2, against w3/6 for the Edgeworth correction.
7. **Theorem 7 claim 2 holds only in the expected sense.** The pushed-forward posterior of z_l is non-Gaussian, and
   its mean quadratic variation is random. Q_l(t) is the *expected* rate, and "acts as a Chen-Eldan localization of the
   surrogate" is a statement about expectations. This should be said, because the non-commuting-path claim is about
   that expected path only. It is also mild: |[C_1, C]| relative 0.05-0.10.
8. **Conjecture 13 is false as stated** (`chk_rank.py`; n = 40, first-Mehler-order chain, d/d mu_c by central
   differences).

   | after | numerical rank of the covariance tangent | sv_3 / sv_1 |
   |---|---|---|
   | 1 transport | 2 | 1e-10 |
   | 2 transports | 40 of 40 | 0.30 |
   | 3 transports | 40 of 40 | 0.36 |

   - **The reason.** After one transport the mean and variance tangents are full vectors. The gain tangent
     delta(A_1 A_1^T) o R' is then the Hadamard product of a rank-2 kernel with a full-rank correlation R', which is
     full rank.
   - Row-locality, and the "all n directions for n^3 per layer" economy, hold only at the first layer of the suffix.
   - So the claim in section 6, that "k directions instead of n is the wrong economy for analytic tangents", is right
     only for the two-layer suffix. For any longer suffix the k-direction economy is the only one, analytic or not.
9. **The rolling result's "stable slope" is overstated.** In `e4_rolling_k3_n48L8.txt` the per-network optimal slopes
   for "last 4" are -1.1, -0.3, -0.4 and -1.1, and for "last 1" they reach +0.1. The quoted -0.6 to -0.9 is the
   leave-one-out average of three networks.
   - Two of the four networks are not improved (1.085, 1.010).
   - The mean ratio 0.884 has a standard error of about 0.1 from n = 4.
   - The per-layer correlation of the local defect with the layer error is weak (median about -0.2, r^2 about 0.04).
   - This is not evidence for a truth-free class constant.
10. **The A-late extension of the decision rule is invalid.** Section 9 records X_late as "the late-visibility
    measurement for every future late-layer design (A-late included)", and "X_late < 0.10: only input tangents can
    see it". Both are wrong.
    - A-late (SYNTHESIS A.2) evaluates late local defects with input (first-chaos) tangents of the late tables. These
      are not functions of (s_l, suffix), so Proposition 11 does not cover them.
    - The author's own measurements separate the two:
      - GC: Lambda_l (input tangents, mostly X_l) explains 0.27-0.94 of the total, while the S-frozen field defect
        explains 0-0.18;
      - K3: s1 0.05-0.18 against 0.00-0.04.
    - So the S-frozen X_late bounds only S-frozen (field) designs: L1, L2 and L3. It says nothing about A-late, which
      stays gated by the pending F1 X* run and the defect-injection test.
11. **Novelty: nothing here is new mathematics.**
    - Field localization is Chen-Eldan with a constant control.
    - The "KMS embedding" is whitening.
    - The Halmos pair is Jordan/Hotelling principal angles (CCA).
    - Dressing is the Kushner-Stratonovich innovation Gram.
    - Theorem 12 is F1 Proposition 7 in the C_l metric.
    - The new content is empirical, and it is the S = 0 replay of F1's finding at a cut, in whitened coordinates.

### R3. Cost-accounting errors

1. **L3' (+20-40 units) has no support** once Conjecture 13 fails.
   - **What it needs.** The C_l-metric field defect of a suffix of 3 or more layers is
     <d_Sigma E, C> - (1/2) sum_j lambda_j d^2_{u_j} E. The first term is one directional derivative along C. The
     second, by a forward Laplacian (Li et al. 2023), needs k first-order tangents carried through every suffix layer
     plus one accumulated second-order object.
   - **Cost per tangent per full layer.** The covariance channel alone is a dense W^T X W, 2 units. The production
     source and readout tables add an estimated 1-3 units.
   - **From cut 12** (three full layers plus the trimmed one), that is about 10-16 units per tangent. At k = 16-32
     the total is **about +160-500 units**, not +20-40. A finite-difference version (L3) is about the same, so
     forward mode buys only a constant factor.
   - **Break-even.** This needs X_late >= 0.42 (+150), 0.49 (+200) or 0.59 (+300), against a prediction of 0.08.
   - **Corrected rule.** Build nothing late-field-shaped unless held-out X_late >= 0.45 at k <= 32. Under the
     author's own prior that outcome has probability under 0.05.
2. **L2 cost (+10-20 units per layer) is unsupported but plausible.** The one-step defect is genuinely row-local, so
   the Theorem 12 economy applies. But the production version (D3' and g4' tangents, hierarchy residues) is a sketch.
   Even at the optimistic 5 units per layer, the 4-network ratio of 0.88 gives an adjusted ratio of 0.965: a 3.5%
   gain, inside the uncertainty of a 4-network estimate.
3. **The experiment's evaluation count is low.**
   - "4 x 64 + 1 + 16" omits that each h/2 or unfrozen direction needs 4 evaluations. The true count is
     256 + 1 + 32 + 32 = 321 per cut, about 165 + 36 = 200 chain-equivalents per network, and about 1600 in total
     (about 9 core-hours).
   - This is still trivially cheap and does not change the conclusion.

### R4. The decisive experiment: fixable defects

1. **The known late kink is not addressed.**
   - The smoothness pilot (task brief) found that from layer 9-10 on, a kink of about 1e-7, suspected to come from the
     nested tier-2 compression of sources older than 8 layers, dominates second differences. Cuts 12-14 sit squarely
     inside it, and section 9 never mentions it.
   - Take h_m = 0.02 sqrt(mean var) along a unit vector, with mean var of order 1. A 1e-7 jump then enters d_j as
     about 1e-7/(2 h_m^2), roughly 1e-4. That is two orders above the per-direction defect scale measured at the
     input (rms about 1e-6).
   - The h/2 check would detect this but not cure it.
   - **Fix (a).** Before production runs, toggle the tier-2 compression (exact old sources on a single network) or
     freeze its basis, selected rank and any data-dependent max or clip. Show h, h/2, h/4 agreement at cuts 12-14.
   - **Fix (b).** Use whitened steps h_j = 0.1-0.2 sqrt(lambda_j) in each mode. The field trace is exactly a sum of
     whitened second differences, the frozen chain is analytic, and the O(h^2) bias is removed by (h, h/2)
     Richardson. A residual jump's effect then falls by (h_new/h_old)^2, a factor of about 50-500 for the top modes.
   - **Fix (c).** Mean-direction second differences with mask freezing may still hit a frozen-basis convention. The
     frozen-basis tangent differs from the derivative of a moving top-r projection, so one convention has to be
     declared and kept.
2. **The tail beyond the top 64 modes is never measured.**
   - T_64 is used as "full", and T_k/T_64 cannot reveal mass in modes 65-1024, which are 94% of the modes. Flatness
     of d_jj over 64 modes at n = 64 does not imply flatness over 960 tail modes at n = 1024. If d_jj ~ 1/lambda_j in
     the tail (an isotropic whitened defect), the tail dominates the trace.
   - **Add:** 16 Hutchinson probes in the complement, v = P_perp C^{1/2} g. That is 64 more evaluations per cut, and
     it bounds the tail share with a bootstrap interval.
3. **The creation horizon is not sampled.** At small n the visible share jumps only once the suffix covers about 5 of
   8 layers (0.44-0.69). Cuts 12-14 (suffixes of 2-4 layers) cannot see such a horizon.
   - **Add cut 10** (a 6-layer suffix, about 0.38 of a chain per evaluation, +80 chain-equivalents per network).
   - Report X_late(l) as a curve. Without it the "< 0.10 closes the frame" branch can misfire.
4. **The scope of the decision rule needs correcting:** see R2.10 (A-late) and R3.1 (the 0.45 threshold).

### R5. Relevance at the 1e-8 / adjusted level

- **What the detector can see.** The S-frozen late field defect sees only error its suffix *creates*, by
  construction (Proposition 11). It cannot see errors in the frozen tables S_l.
- **Where the production error comes from.** The oracle attribution puts about 40% of the production MSE in the
  kappa_3 readouts (D3, D21) and 15-30% in the kappa_4 diagonal, mostly created by transport and compression across
  many layers. The DATA_FINDINGS injection profile is about readout, not creation, as the author correctly argues.
- **Expected effect: nil.** Raw final-layer MSE would change by at most X_late x 1.55e-8, about 1e-9 at the predicted
  X_late. At the corrected costs every late-field component (L1-L3') raises adjusted MSE.
- **Value of the experiment.** It is diagnostic only: it settles whether production error is suffix-created. It does
  *not* gate A-late (R2.10). Since its engineering (a restartable, mask-frozen float64 chain) is shared with the
  pending F1 X* run, it is worth running as a 15-minute add-on with the R4 fixes, and not worth anything on its own.

### R6. Is the noncommutative content genuine?

| claimed NC object | what it is | does it change a computation? |
|---|---|---|
| D~ = C^{1/2} D C^{1/2} "KMS / geometric embedding" | the defect in whitened coordinates | no: every operator-mean embedding gives the same trace, and no reported number depends on the choice |
| E_C(D), "spectral paths see only the commutant" | tr(D Q) = tr(E_C(D) Q) for [Q, C] = 0 | no: a tautology for any trace detector |
| Kubo-Ando C # C_1 as a control | a non-commuting mean of two covariances | it was tested and lost (0.77-0.95, against 0.77-0.97 for the field trace), so it buys nothing |
| Halmos pair (P_{H_1}, P_Z) | canonical correlations of X and z_l (Gamma). The algebra the two projections generate is type I: M_2 over the angle spectrum plus a commuting part | only through the choice of metric, C against C_1, which is ordinary linear algebra. The "keyhole" (median eigenvalue of Gamma about 0) is a correct and useful diagnostic, but it is CCA |
| Theorem 7 non-commuting induced path | the expected posterior covariance path rotates from the C_1 frame to the C frame | no: mild (0.05-0.10), and nothing is computed from it |
| dressing (Theorem 9) | the field innovation Gram replacing the input Gram | it is a derivation of note XLI's rule (commutative filtering), not an NC effect, and it needs the true next-order cumulant to beat bare |

**Verdict on NC.** The frame's own ledger concedes most of this. Its remaining claim, that the Halmos pair is "the
real one", is correct only in the trivial sense that two generic projections do not commute. Nothing in the
estimator or the error theory uses more than the eigenvalues of Gamma. **The NC content is relabelling.**

### R7. Strongest surviving idea (corrected)

The field-metric (whitened) trace of the S-frozen late defect, tr(D_l C_l) = <d_Sigma E, C> - (1/2) sum_j lambda_j
d^2_{u_j} E, has three uses:
- it is the exact heat-defect of a *suffix* (proved and checked for S = 0, with shifted Green measure from
  y_0 = C^{-1/2} mu);
- with the cut decomposition (Theorem 10) it gives a clean three-way split of the input defect: transported prefix
  defects, input-tangent suffix terms, and the dressed Gaussian block;
- measured on the production chain at cuts 10-14, with whitened Richardson steps, a complement-tail probe and the
  kink fixed, it yields the curve X_late(l), the share of error the last layers *create*.

That curve is a cheap, decisive diagnostic for the whole class of late S-frozen components (L1-L3'). It is not a
component. Every late-field component is negative at corrected costs unless X_late >= 0.45. It does not bound A-late.

### R8. Verdict

- **Mathematics.** Correct apart from the Lemma 1(a) sign, the S = 0-only headline identity and the false
  Conjecture 13. Several statements presented as theorems or measurements are tautologies (Theorem 5 claim 3, the
  diagonal-carries-all-power measurement, Proposition 4(3)).
- **Overclaims.**
  - the rolling slope stability;
  - "path gauge";
  - extending the decision rule to A-late;
  - "k-directions is the wrong economy for analytic tangents" beyond two layers;
  - the L3' cost.
- **NC content.** Relabelling: whitening, CCA and Kushner-Stratonovich.
- **Relevance.** Nil as a component. Diagnostic only, as an add-on to the F1 X* run.
- **Survives as:** theory plus a corrected optional diagnostic. **Priority: low (2/10).**
