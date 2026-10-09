# F2. Activation-pattern pinning and the transverse measure of linear regions

Frame F2 of the localization workflow (9 October 2026). Code and raw outputs: `loc/f2/` (all runs local, tiny or
medium He networks; nothing here touched the official networks). Status labels: **proved** (complete proof here),
**checked** (numerically verified, output quoted), **sketch**, **conjecture**.

## 0. Verdict

The network's own localization is the filtration of its gates. It has a clean and complete exact theory, and that
theory ends at a quantified wall rather than at a new estimator.

1. **Structure (proved).** The patterns form the *centre* of the commutant of the pattern groupoid. Classical input
   localizations (Eldan's process, linear tilts) live in the right Cartan subalgebra and all refine the pattern
   localization. The leaf space is finite, so this is the type-I corner of stage-8 Theorem A.
2. **The mean lives on the codimension-one skeleton (proved, checked to 1e-12).** E[F] is the transverse integral,
   over the facets between adjacent linear regions, of the jump of the leafwise Jacobian (Theorem 1). In tropical
   form it is a difference of first intrinsic volumes (Gaussian mean widths) of two Newton polytopes (Proposition 2).
   The transverse measure of the linear regions is the external-angle measure of those polytopes.
3. **Gate pinning (proved, checked).**
   - It is an exact linear tilt with an exact trickle-down identity (Proposition 3).
   - The gate covariance is exactly the trickle-down of a common-mode localization: it is the covariance that pinning a
     latent removes (Proposition 4).
   - Conditioning is not freezing: E[D A D] = P A P + G o A, a Schur-positive excursion term (Proposition 5).
4. **Mean field and its hierarchy (proved, with one numerical demonstration).**
   - The mean-field (product-gate) law is the chain's closure: the Gaussian closure conditions at order 2 (it uses
     exact pair gate laws) and the kappa3 chain freezes at order 3 (product gates on transported sources).
   - The pinning cluster expansion is the Mehler expansion. Its first omitted level is the pair-gate class, the V2
     move: about a third of the chain's MSE by the closure-round measurement.
   - The kink (facet) grading is not law-intrinsic. Its mean-field truncation is off by O(1) (rms error 0.35-0.9
     against 0.03-0.04 for the Gaussian closure), so it is not a usable grading.
5. **The main quantitative result: the separating dimension (proved, plus a conjecture with numerical support).**
   - The cost of carrying the pair-gate class in any representation that keeps CP legs is set by K_f: the smallest
     latent dimension that makes a layer's gates conditionally independent, to a fraction f of the class in the
     readout metric. The fresh-weight proposition (Proposition 6) proves that the captured fraction is a weighted
     Frobenius mass. This explains note XX's rank-K table.
   - K_f is extensive: K_90 grows linearly in n (n = 128 to 1024, section 4.3), while lambda_max of the gate
     correlation stays at 3-10.
   - K_90/n is n-independent over n = 128-1024: in the gate-covariance metric, 0.11 at layer 10 and 0.08 at layer 14.
     At n = 1024, K_90 is 55-114 latent directions at layer 10 (correlation metric and gate-covariance metric), each
     costing one extra leg transport per source-layer. Exact carrying is a K_4 contraction at n^4 per source-layer.
   - A fresh random n = 1024 He network reproduces the official network's rank-K masses (note XX) to within 0.013,
     so the separating spectrum is annealed and predictable even though the correction is quenched.
6. **System consequence: a clear negative.** The frame yields no deterministic component that pays under the score.
   It does three things instead:
   - it closes the "pin the gates" design direction with a stated reason;
   - it gives a falsifiable scaling law;
   - it specifies one offline experiment that fixes the value of the wall. That value decides whether a
     representation without the product-gate defect is the thing to look for.

## 1. Setting

- **The network.**
  - Input x = x_0 ~ N(0, I_n).
  - Layers z_l = x_l W_l and x_{l+1} = relu(z_l), for l = 0..L-1. The row convention is used, and F = x_L.
  - Gates s_l = 1{z_l > 0} in {0,1}^n, with D_l = diag(s_l).
- **The leafwise map.** On the set where the pattern s = (s_0..s_{L-1}) is fixed, F(x) = x J_s with
  J_s = W_0 D_0 W_1 D_1 ... W_{L-1} D_{L-1}. So F is 1-homogeneous, continuous and piecewise linear.
- **Downstream gains.** G_{l+1} = d x_L / d x_{l+1} = W_{l+1} D_{l+1} ... W_{L-1} D_{L-1}, with G_L = I. It is
  n x n, indexed [j, i] for post-activation j of layer l+1 and output i.
- **Regions.** R_s = {x : pattern(x) = s} is an open polyhedral cone, because every constraint z_{l,j}(x) > 0 or < 0
  is linear on the cone fixed by the earlier gates.
- **Generic position (G).** No layer is entirely inactive on an open set. When (G) fails, x_l = 0 on a cone and every
  downstream z vanishes identically there. Under Gaussian W, (G) fails only on the event that some layer is dead on a
  cone of positive measure. That event is common at n = 8, where it broke the first checks; it occurs for 1.5e-6
  of samples at n = 16 (n1_kink16), and its probability is astronomically small at n = 1024.

## 2. The leaf space and its transverse measure

### 2.1 The pattern groupoid and the commutant (Theorem A, specialised; proved)

Let R_l = {(x, y) : s_{<=l}(x) = s_{<=l}(y)}. This is a measured equivalence relation on (R^n, gamma). Its classes are
the cones R_s, finitely many and of positive measure almost surely. The Haar system is gamma restricted to the
classes.

- **The algebra and its commutant.** The groupoid von Neumann algebra acts on L^2(R_l) = (+)_s L^2(R_s) (x) L^2(R_s)
  as L(R_l) = (+)_s B(L^2(R_s)) (x) 1. Its commutant is (+)_s 1 (x) B(L^2(R_s)).
  - The right Cartan subalgebra of the commutant is (+)_s 1 (x) L^\infty(R_s), which is L^\infty(R^n, gamma).
  - The centre is l^\infty(patterns up to l), which is L^\infty(Pat_l) with Pat_l = sigma(s_0, ..., s_l).
- **What stage-8 Theorem A says here.**
  - (i) Classical localizations of the input law are filtrations in the right Cartan subalgebra. Eldan's process
    (X | Y_t ~ N(Y_t/(1+t), I/(1+t))) and every linear-tilt scheme are of this kind.
  - (ii) The pattern filtration Z_0 ⊂ Z_1 ⊂ ... ⊂ Z_{L-1} is the central localization: a finite tower of
    finite-dimensional commutative algebras (a commutative AF tower). Its Bratteli diagram is the tree of feasible
    pattern prefixes, with multiplicity one and at most 2^n children per vertex. For neuron-by-neuron revelation it
    has nL levels and is binary, with infeasible children pruned (an empty cone).
  - (iii) Every localization that ends at point masses refines the central one, because a point determines its
    pattern.
  - (iv) The transverse measure is the law of the pattern, mu(s) = gamma(R_s). Since F is linear on each leaf,
    E[F] = sum_s mu(s) F(b_s), with b_s = E[x | R_s].
- **Honesty about the noncommutative content.** The leaf space is finite and the algebra is type I. Connes'
  machinery reduces to finite sums, and nothing here is a bad quotient.
  - The only genuinely noncommutative object would appear as depth goes to infinity: the tail relation of gate
    patterns, "agree from some layer on", is a hyperfinite (AF) relation with a nontrivial transverse measure. It is
    irrelevant at L = 16.
  - What does carry content is the codimension-one structure of section 2.2.

### 2.2 The codimension-one skeleton carries the mean

**Theorem 1 (facet / kink representation of the mean).** *Proved; checked to 1.6e-12.*

Let F be a bias-free ReLU MLP and x ~ N(0, I_n). Write Sigma for the set of facets: the (n-1)-dimensional cells of the
region complex. For a facet sigma between R_s and R_s', let nu_sigma be its unit normal pointing into R_s, and let
[grad F_i]_sigma = J_s e_i - J_s' e_i be the jump of the gradient. Then

    E[F_i] = sum_{sigma in Sigma} <[grad F_i]_sigma, nu_sigma> * gamma_{n-1}(sigma),
    gamma_{n-1}(sigma) = int_sigma phi_n dH^{n-1} = (2 pi)^{-1/2} * pi_sigma.          (1)

Here pi_sigma is the probability that a standard Gaussian vector of the hyperplane nu_sigma^perp lies in the facet
cone: a Gaussian solid angle.

Under (G), every facet lies in the zero set of exactly one neuron (l, j), and only gate s_{l,j} flips across it. The
neuron-wise form is then

    E[F_i] = sum_{l=0}^{L-1} sum_{j} E[ delta(z_{l,j}) |grad_x z_{l,j}|^2 G_{l+1}[j, i] ],               (2)

where delta(z)|grad z|^2 dgamma denotes the measure |grad z| phi_n dH^{n-1} restricted to {z_{l,j} = 0}.
Equivalently, E[F] = E[Delta F], with Delta F the distributional Laplacian supported on the kink set.

*Proof.*

1. **One cone.** F is 1-homogeneous, so F_i(x) = <x, grad F_i(x)> almost everywhere (Euler), and grad F_i = J_s e_i
   is constant on the open cone R_s. The Gaussian divergence theorem on a polyhedral cone gives
   E[x 1_{R_s}] = int_{R_s} x phi_n = -int_{R_s} grad phi_n = int_{partial R_s} phi_n nu_in dH^{n-1}. Faces of
   dimension n-2 or less carry no H^{n-1} mass.
2. **Sum over cones.** E[F_i] = sum_s <J_s e_i, E[x 1_{R_s}]>. Each facet is counted twice, with opposite inward
   normals, which gives (1).
3. **The neuron-wise form.** Under (G), a facet lies in {z_{l,j} = 0} for one neuron. On it,
   x_{l+1,j} = relu(z_{l,j}) = 0, so every downstream pre-activation is continuous across it, and generically none of
   them vanishes there. The upstream gates are constant across it. The gradient therefore jumps by
   grad z_{l,j} G_{l+1}[j, i], and nu = grad z/|grad z|, so <jump, nu> = |grad z_{l,j}| G_{l+1}[j, i].
4. **Both orientations.** The sign works out the same way for both orientations of the crossing.

*Checks.*

| check | result |
|---|---|
| exact, input dimension 2, width 16, depth 4, two seeds (`n1b_kink2d.py`; E[F] = sqrt(pi/2) times the circle average of f, kinks bisected to machine precision) | 156 and 136 kink rays; max \|E[F] - (1)\| = 1.6e-12 and 2.5e-12 |
| width 8, where a layer is entirely inactive on an arc, so (G) fails | (2) fails by 7e-2. The facet form (1) still holds, but the boundary of the dead sector is a facet shared by many neurons, and (2) mis-assigns it |
| Monte Carlo, n = 16, L = 3, 4e6 samples, delta estimated by windows at eps = 0.08, 0.04, 0.02 sd (`n1_kink16.py`) | window values 0.389, 0.379, 0.372 against E[F] = 0.368 (output average), converging. E[Y \| z] has a slope jump at z = 0, because downstream gates move on one side only, so the window bias is first order in eps |

**Relation to the literature.**
- Hanin and Rolnick ("Complexity of linear regions in deep networks", 2019) use the same co-area structure, with
  weight |grad z| instead of |grad z|^2 G, to count the (n-1)-volume of the kink set.
- Özkan and Hirsch (arXiv 2609.32695, 2026) give a conditional Kac-Rice formula for that measure.
- Formula (2) is the version weighted by the Gaussian and the jump, and that weight is what yields the mean.

### 2.3 The tropical home: the mean is a mean-width difference

**Proposition 2.** *Proved; the classical inputs are cited.*

1. **The recursion.** Every coordinate of every layer is a difference of support functions,
   x_{l,i} = h_{P_{l,i}} - h_{Q_{l,i}}, where h_K(x) = max_{v in K} <v, x> and P, Q are polytopes. The recursion:
   - base: x_{0,k} = h_{{e_k}} - h_{{0}};
   - linear step: z_{l,i} = h_{A_i} - h_{B_i}, with A_i = sum_a (W+_{ai} P_a + W-_{ai} Q_a) and
     B_i = sum_a (W+_{ai} Q_a + W-_{ai} P_a) (Minkowski sums, W+- the positive and negative parts);
   - ReLU step: x_{l+1,i} = h_{conv(A_i ∪ B_i)} - h_{B_i}.
2. **The mean.** E[h_K(x)] = V_1(K)/sqrt(2 pi) (Sudakov-Tsirelson), so

       E[F_i] = (V_1(P_{L,i}) - V_1(Q_{L,i})) / sqrt(2 pi).

3. **Which steps are exact.** V_1 is Minkowski-additive and positively homogeneous, so the linear layer is exactly
   linear in the pair (V_1(P), V_1(Q)). The only non-additive step is the convex hull of a union, which is the ReLU.
4. **Leaves and facets.** The linear regions of F_i refine the normal fans of P_{L,i} and Q_{L,i}, and facets of the
   fan correspond to edges of the polytopes. Formula (1) restricted to h_P is the edge formula
   V_1(P) = sum_{edges e} |e| gamma_ext(e, P): the jump |v - v'| = |e|, and pi_sigma is the external angle of the
   edge.

*Proof.*
- relu(h_A - h_B) = max(h_A, h_B) - h_B = h_{conv(A ∪ B)} - h_B; c h_K = h_{cK} for c >= 0; h_K + h_K' = h_{K+K'}.
  This is the bias-free case of Zhang, Naitzat and Lim, "Tropical geometry of deep neural networks" (2018).
- The mean-width identity is Tsirelson's theorem. Check on a segment: E relu(<w, x>) = |w|/sqrt(2 pi) =
  V_1([0, w])/sqrt(2 pi).

**Reading.** The transverse measure of the linear regions is the external-angle measure of the Newton polytopes, and
the mean is the Gaussian mean-width difference. This is the natural home for the frame's "a measure on the space of
leaves passes to the leaf space": the leaves are vertices, and the integral is carried by the edges.

The Newton polytopes have astronomically many edges, so this is a home and not an algorithm. It does fix one
structural fact used below: everything linear is free (V_1-additive), and all difficulty sits in conv(A ∪ B), that is,
in the joint law of the gates with the field.

## 3. Pinning the gate law

### 3.1 Linear-tilt form and exact trickle-down

**Proposition 3.** *Proved.*

1. **The tilt.** Revealing one gate s (with p = P(s = 1) in (0,1)) at the value sigma in {0,1} is an exact linear
   tilt:

       d nu_sigma / d nu = 1{s = sigma}/P(s = sigma) = 1 + (s - p) Z,   Z = (sigma - p)/(p(1 - p)),
       E[Z] = 0,   E[Z^2] = 1/(p(1 - p)) =: R.

2. **Bank trickle-down.** For any bank of square-integrable observables a = (a_1..a_m), with V = Cov(a) and
   C = Cov(s, a) (1 x m):

       V - E[ Cov(a | s) ] = C^T R C = Cov(a, s) Cov(s, a) / (p(1 - p)).

3. **Random-order pinning.** For a uniformly random gate g of a layer (N gates, D = diag(p_g(1 - p_g))), the
   bank form applied to the gates themselves gives

       Sigma - E[Sigma_after] = (1/N) Sigma D^{-1} Sigma.

   This is the AKV equation with a hypercube in place of a simplicial complex. As a pure complex: a pattern is a
   facet of the boundary of the n-cross-polytope (choose +e_a or -e_a per coordinate), and pinning is the Gelfand-
   Tsetlin localization of stage-8 Theorem E on its face-poset AF algebra.

*Proof.* For binary s, E[a | s] is affine in s, so Cov(E[a|s]) = Cov(a, s)Cov(s, a)/Var(s). This is the law of total
covariance. The noncommutative version of the brief (V_t - E V_{t+1} = C* R C) reduces to this on the diagonal; no
noncommutativity is involved.

**Tilt of the next layer by one pinned gate.** Pinning gate (l, j) moves the law of x_{l+1} by
Cov(x_{l+1}, s_{l,j}) Z. Here Cov(relu(z_{l,k}), s_{l,j}) is a bivariate Gaussian integral under the closure.
- **The tilt matrix.** T_l = Cov(x_{l+1}, s_l) is n x n, costs O(n^2) per layer, and reaches the pre-activations of
  layer l+1 as T_l W_{l+1} (one n^3 product).
- **What it is in the chain.** It is a V1 move with a localized centre: the chain's star and fold already carry
  that class.

### 3.2 The gate covariance is the trickle-down of a common mode

**Proposition 4.** *Proved; checked to six digits.*

- **The factor law.** Under a Gaussian closure z ~ N(mu, D + B B^T), with D diagonal and B of size n x K, write
  z = mu + B X + D^{1/2} eps with X ~ N(0, I_K). Given X the gates are independent, with
  p_a(X) = Phi((mu_a + B_a X)/sqrt(D_a)).
- **The identity.** By the law of total covariance,

      G_ab := Cov(s_a, s_b) = Cov_X(p_a(X), p_b(X))   (a != b),        G_aa = E_X[p_a(1-p_a)] + Var_X p_a.

  The off-diagonal gate covariance is exactly the covariance that a localization of X removes. In AKV form: under
  Eldan's process on X, the trickle-down term C^T R C integrates to Cov(E[s | X]).
- **First order.** To first order in the correlation, G_ab = phi(alpha_a) phi(alpha_b) rho_ab + O(rho^2). This is
  the first Mehler term,

      P(z_a > 0, z_b > 0) = Phi(alpha_a)Phi(alpha_b) + sum_{k>=1} rho^k/k! He_{k-1}(alpha_a) He_{k-1}(alpha_b) phi phi.

Check (`n3_gatecov.py`, n = 6, K = 2, 60^2-node Gauss-Hermite in X):

| (a,b) | exact bivariate orthant | Cov_X(p_a, p_b) | phi phi rho | rho |
|---|---|---|---|---|
| (0,1) | +0.006300 | +0.006300 | +0.006523 | +0.094 |
| (1,3) | -0.031810 | -0.031810 | -0.031725 | -0.246 |
| (4,5) | -0.019333 | -0.019333 | -0.019293 | -0.124 |

**Corollary (capacity).** Suppose a K-dimensional linear-Gaussian latent renders a layer's gates conditionally
independent. Then the linearized off-diagonal gate covariance D_phi rho D_phi has rank at most K. So a K-dimensional
separating localization captures at most the top-K (weighted) Frobenius mass of G_off, at first order. This is the
classical, commutative form of the stage-8 capacity theorem (Theorem C) for this observable bank.

### 3.3 Conditioning is not freezing

**Proposition 5.** *Proved.*

For a random diagonal gate D = diag(s) and any fixed matrix A,

    E[D A D] = P A P + G o A,     P = diag(p),  G = Cov(s)  (G_aa = p_a(1 - p_a)).

If A is positive semidefinite, G o A is positive semidefinite (Schur product theorem). This is the gate analogue of
q A^2 q - (q A q)^2 = q A (1 - q) A q >= 0: the excursion term that freezing drops.

*Proof.* (E[D A D])_ab = A_ab E[s_a s_b], and s_a^2 = s_a.

**What each estimator freezes.** Order by order in the transported cumulants:

| order | object transported through the gates | Gaussian closure | kappa3 chain | exact |
|---|---|---|---|---|
| 1 | mean, first chaos | per neuron: exact | exact (Stein), plus the fold | — |
| 2 | covariance | **conditioned**: the bivariate law, arcsine kernel, exact G | conditioned | — |
| 3 | bulk kappa3 legs | absent | **frozen**: Phi_a Phi_b Phi_c T_abc | E[s_a s_b s_c] T_abc |

The frozen-minus-conditioned difference at order 3 is sum_pairs G_ab Phi_c T_abc + kappa3(s) T. Its first term is
the V2 move (the gate covariance of notes XVIII-XX, the joint-gate triangle tr H^3 of note XLI).

**Conditioning can also be freezing.** Conditioning on the magnitudes |z_l| determines the signs up to a
zero-temperature dense Ising law: note XX measured couplings of 10^2-10^3 and Dobrushin sums near 1000. A
localization finer than the gates (magnitudes) freezes. The pattern localization itself is not frozen: the gate law of
a layer is weakly correlated, with lambda_max of the gate correlation 3.4-10 (section 4.3).

### 3.4 Local-to-global (Oppenheim/AKV) statements, and why they do not bear on bias

- **The measured constant.** The spectral-independence constant of a layer's gate law is
  eta_l = lambda_max(diag(G)^{-1/2} G diag(G)^{-1/2}). Measured (`n4_sepdim.py`): 3.4-3.6 at layer 1, rising to 6-10
  at depth (up to 10.9 at n = 1024), with no growth in n over n = 128-1024.
- **What it gives.** Via trickle-down on the cross-polytope complex, provided the links satisfy similar bounds (not
  checked), this gives polynomial mixing of gate Glauber dynamics with an exponent O(eta).
- **Why it does not bear on bias.** Every such statement controls the decay of covariance under pinning, that is,
  variance and mixing. None controls the bias of a deterministic closure.
- **The exact version.** The mean-field error of a k-gate product is the telescoping sum of one-pin tilts
  (P(s_b = 1 | s_a = 1) - p_b = G_ab/p_a). It is O(rho) per pair and sign-indefinite, and section 4 measures it
  directly.

**Conclusion.** The Oppenheim machinery has no estimator content at the 1e-7 sigma^2 level.

## 4. Mean field, the pinning hierarchy, and the separating dimension

### 4.1 Mean field and the cluster expansion

The mean-field gate law is the product law ⊗_a Bern(Phi(alpha_a)) given the closure's marginals, applied to every
transported object of order 3 or more.

**Cluster expansion (proved).** By the moment-cumulant formula,

    E[prod_{a in S} s_a] = sum over set partitions pi of S of prod_{blocks} kappa(s_block).

Under the Gaussian closure, kappa(s_a, s_b) = G_ab = phi phi rho + O(rho^2) and kappa_3(s) = O(rho^2). So truncation
at pair clusters is the first-order Mehler correction. A Kikuchi or cluster-variation hierarchy over gate subsets of
size at most k reproduces the Mehler hierarchy to order k-1 in rho, because a connected k-cluster needs at least
k-1 correlation lines. Note XIII E3 measured the next order at 0.08-0.10% of the gate.

In the Kikuchi-Toeplitz language of stage-8 Theorem D, the degree-2 Walsh projection of the gate law is exactly G.
The Berezin-Lieb sandwich bounds convex functionals of such projections, but the transported readouts are linear in
G, so the sandwich says nothing about them.

**The first omitted level.** It is the pair-gate class. Its magnitude was measured once on official network 0, from
the live legs of layer 10 (closure round, section 5):
- the correction is 1.09% of the product-gate D3 readout, incoherent (correlation 0.09), with mean 0.04 of its rms;
- it injects 1.09e-5 rms into the post-activation mean at layer 11;
- with persistence 0.85 per layer, about 6.5e-9 of the then 2.2e-8, roughly a third of the MSE.

### 4.2 The kink grading is not a projection (a negative)

Formula (2) suggests a grading of E[F] by kink layer, and a mean-field version that factorizes each Palm term:
- MF1: E[delta |grad z|^2] E[G];
- MF2: p(0) E|grad z|^2 E[G].

Given the kink weights, MF2 costs O(L n^2) through one backward sweep. The weights themselves need the Jacobian Gram
E|grad z|^2, about one n^3 product per layer. These truncations fail badly (`n1_kink16.py`, n = 16, L = 3; output
average E[F] = 0.368):

| estimator | rms error over outputs |
|---|---|
| Gaussian closure (law-based) | 0.037 |
| kink sum, Palm-exact (window Monte Carlo) | <= 0.015 (window bias) |
| MF1: downstream gain independent of the pinned kink | 0.348 |
| MF2: fully factorized Palm | 0.607 |

At n = 8: 0.033, 0.010, 0.405 and 0.882.

The per-layer exact terms (-0.066, -0.097, +0.528) against MF1's (+0.081, +0.077, +0.528) show the mechanism.
- Pinning a hidden neuron at its kink switches off its own downstream paths, so the Palm-conditioned downstream gain
  differs from the unconditioned one at O(1) once summed over the layer.
- The facet terms are not functions of the law of the pre-activations; they depend on how z is written as a function
  of x (relu(h) - relu(-h) = h has kinks and no mean). So a product closure of the Palm measure is not a projection.
- At n = 1024 the fully factorized version is note XLI's Euler identity mu = tr H, measured at correlation 0.73-0.95
  with the chain's mean. That is the same failure.

**This grading is discarded.**

### 4.3 The separating dimension

**Definition.** For layer l, a separating localization of dimension K is a K-dimensional Gaussian latent X, a
subalgebra of the right Cartan subalgebra generated by K linear functionals of the closure's field, such that the
gates are conditionally independent given X. Under a Gaussian closure z ~ N(mu, C) this means
C - Cov(z, X) Cov(X)^{-1} Cov(X, z) is diagonal, that is, C = D + B B^T with rank(B) = K.

K_f is the smallest K whose best separating localization captures a fraction f of the pair-gate correction in the
readout metric (Proposition 6).

**Proposition 6 (fresh-weight capture).** *Proved; checked.*

- **The objects.** Let M and M' be 3-tensors supported on all-distinct index triples, and let w ~ N(0, (2/n) I) be
  independent of them: the next layer's weight column.
- **The moments.** Then E[M[w,w,w]] = 0, and

      E[M[w,w,w] M'[w,w,w]] = 6 (2/n)^3 <Sym M, Sym M'>_F.

- **The readout.** The pair-gate correction to the next layer's kappa3 diagonal is
  D3corr_i(rho) = M_rho[w_i, w_i, w_i], with M_rho = Sym(rho_ab phi_a phi_b Phi_c T_abc) linear in rho. Over output
  neurons i, the regression slope of D3corr(rho_K) on D3corr(rho) converges to

      <Sym M_{rho_K}, Sym M_rho>_F / ||Sym M_rho||_F^2 =: weighted mass of rho_K,

  with weights omega_ab = sum_c (phi_a phi_b Phi_c T_abc)^2 (up to symmetrization). In general the squared
  correlation is slope x <M_K, M>/||M_K||^2; it equals the slope when rho_K is the omega-orthogonal truncation.

*Proof.* By Isserlis, for all-distinct triples (a,b,c) and (a',b',c'), E[w_a w_b w_c w_a' w_b' w_c'] is nonzero only
when the two index sets coincide, and then equals (2/n)^3 times the number of matchings. The w_i are i.i.d. and
independent of M, so the empirical slope over i converges to the population ratio.

**Consequence.** The pair-gate correction has no annealed part: its mean over the fresh weights is zero. Its capture
by a rank-K truncation is a weighted Frobenius mass.

Checks:
- **Local (`n3_gatecov.py`, n = 512).** Layer 9 of a He network, a 48-term source with covariance-row arms and a
  dense centre leg, 192 fresh neurons:

  | K | unweighted mass | weighted mass | slope | correlation |
  |---|---|---|---|---|
  | 1 | 0.302 | 0.019 | 0.013 | 0.13 |
  | 4 | 0.503 | 0.409 | 0.422 | 0.71 |
  | 16 | 0.820 | 0.814 | 0.834 | 0.94 |
  | 64 | 0.933 | 0.926 | 0.929 | 0.99 |

  The top mode (the mean-direction spike) holds 30% of rho's mass but 2% of the readout-weighted mass. Gates are
  saturated along the mean direction, where phi is small, so the gain spike is not the separator. This matches note
  XIII E3's "top coherent mode carries 0.1-0.3%".
- **The official networks (note XX, n = 1024, layer 10, real sources).** Unweighted masses 0.19 / 0.35 / 0.67 / 0.91
  gave slopes 0.22 / 0.40 / 0.74 / 0.96 and correlations of about sqrt(mass). Proposition 6 explains that table.

**Conjecture 7 (extensivity of the separating dimension).** *Numerical support plus a free-probability sketch.*

For He ReLU networks at fixed depth l, K_f(n, l) = kappa_f(l) n (1 + o(1)), with kappa_f decreasing in l. The
measured K_90 (`n4_sepdim.py`, Monte Carlo with 32768 inputs, 16384 at n >= 512; the noise mass n^2/N is at most 1%
of the signal):

| layer | quantity | n = 128 | n = 256 | n = 512 | n = 1024 | K_90/n |
|---|---|---|---|---|---|---|
| 6 | K_90(rho_off) | 22 | 39 | 71 | 157 | 0.14-0.17 |
| 6 | K_90(G_off) | 23 | 41 | 83 | 168 | 0.16-0.18 |
| 10 | K_90(rho_off) | 6 | 14 | 24 | 55 | 0.047-0.055 |
| 10 | K_90(G_off) | 15 | 27 | 54 | 114 | 0.105-0.117 |
| 14 | K_90(rho_off) | 5 | 7 | 15 | 36 | 0.027-0.039 |
| 14 | K_90(G_off) | 11 | 22 | 43 | 80 | 0.078-0.086 |
| 14 | lambda_max of the gate correlation | 9.7 | 6.0 | 8.0 | 10.9 | |
| 14 | PR(C) | 5.4 | 10.2 | 22.6 | 65.4 | |

The ratio K_90/n is n-independent to within about 10-20% over an eightfold range of widths. That is the content of
the conjecture.

**The n = 1024 run reproduces the official network.** This is a fresh random He network, not an official one. At
layer 10 its rank-K Frobenius masses of rho_off are 0.199 / 0.350 / 0.669 / 0.913 at K = 1 / 4 / 16 / 64. Note XX
measured 0.186 / 0.353 / 0.671 / 0.912 on official network 0. So the separating spectrum is a self-averaging
(annealed) property of He networks, although the correction it governs is quenched. It can be predicted without the
official weights.

P-F2-1 below (the n = 1024 values: 95-120 at layer 10, 80-95 at layer 14, lambda_max below 12) was written before the
n = 1024 run finished. Measured: 114, 80 and 10.9. It holds.

- **The sketch.** The bulk of the depth-l covariance is that of a product of l gated He matrices. By free
  multiplicative convolution its normalized spectrum is n-independent: Fuss-Catalan for ungated Gaussian products.
  So every spectral fraction, K_f/n included, converges. Check (`n6_fusscatalan.py`, n = 512, ungated products):
  PR = 255, 169, 100, 71, 56, 46, 34 at l = 1, 2, 4, 6, 8, 10, 14, against n/(l+1) = 256, 171, 102, 73, 57, 47, 34.
- **The gates.** They add a spike and shorten the effective depth. The He-network PR is about n/(l+1) at moderate
  depth (n = 256: 45.6, 33.9, 26.0, 21.8 at l = 4, 6, 8, 10, against 51, 37, 28, 23), and is lowered further by the
  spike at depth.
- **What this rules out.** The separator is an extensive collective object, not a small hidden variable. Note XX's
  "the omitted memory is small compared with 1024" was the right number (about 64) at n = 1024. As a structure it
  is about 6-10% of the width, and it grows with the width.

### 4.4 The cost theorem

**Theorem 8 (carrying the pair-gate class).** (a) is proved; (b) is a sketch; (c) is arithmetic.

**(a) Exact.** Per source, the D3 readout of the corrected transport is

    sum_{a,b,r,i} W_ai W_bi R'_ab Y_ar Y_br q_ri.

Its index graph is K_4 on {a, b, r, i}, of treewidth 3. Every pairwise contraction order therefore costs n^4 per
source-layer, and fast matrix multiplication gives n^{1+omega} (about n^3.8 with Strassen blocks).
- At n = 1024, one source-layer costs 1024 units at n^4 (about 600 under Strassen).
- The chain has about 120 source-layers.

**(b) Through a separating localization.** Each latent direction k turns the leg pair into a rank-one modification,
(g_k o A_r) (x) (g_k o A_r). That is one extra leg transport per source-layer, with no sharing between sources,
because the A and P legs differ across sources.
- The only shared savings: arms in the covariance's top-64 subspace (96%, frontier tests, the Lyapunov section) let
  the (A, A) pairing share the products W diag(u_k) U across sources, at about 4 units per source-layer for K = 64.
- The (A, P) pairings have dense P legs and cannot share.
- The cheapest CP carrier therefore costs about K_f transports per source-layer.

**(c) Against the score.**
- One transport over the depth is about 0.2 F, roughly 40 units (note XX).
- A component pays only if it removes more than Delta C/(C + Delta C) of the MSE.
- With the class at about 1/3 of the MSE, capturing slope s removes about (1 - (1-s)^2)/3:

| K | slope (note XX) | MSE removed | extra cost | adjusted change |
|---|---|---|---|---|
| 4 | 0.40 | 21% | +160 units (C/B 0.20 -> 0.36) | x1.42 (worse) |
| 16 | 0.74 | 31% | +640 units (C/B -> 0.83) | x2.8 (worse) |
| 64 | 0.96 | 33% | +2560 units (C/B 2.7, over the budget) | not admissible |
| exact | 1 | 33% | about 1.2e5 units | not admissible |

**No K pays, and K_f grows with n.**

## 5. Estimator designs this frame implies, priced

Costs are in units of 2n^3 (budget 1024). Raw is final-layer MSE at n = 1024.

| | design | cost | expected raw | status |
|---|---|---|---|---|
| E0 | mean-field pattern law = Gaussian closure (pair-exact gates) | about 2 units per layer, 32 in total (C/B 0.03) | 4.1e-6 (note XIII E6) | known; adjusted 4.1e-7 |
| E1 | kink-graded mean field (MF2): Sigma_l k_l^T E[G_{l+1}], one backward sweep, plus a gradient Gram per layer | about L n^3, 16 units | O(1) relative error (section 4.2) | **negative** |
| E2 | the kappa3 chain: freeze at order 3 | about 208 units (C/B 0.203) | 1.55e-8 | existing system |
| E3 | chain plus the pair-gate class through a K-dim separating localization | +40 K units | 1.55e-8 x (1 - (1-(1-s_K)^2)/3) | **does not pay** for any K (section 4.4) |
| E4 | chain plus the exact pair-gate class (K_4 contraction) | about 1e5 units | about 1.0e-8 (prediction) | inadmissible |
| E5 | gate-pattern Monte Carlo for the exact-gate transport: sample s from the closure's orthant law, so the transport keeps CP rank | 1 transport per sample | relative noise O(1) per sample on a 1% effect, so about 1e4 samples | **negative** |
| E6 | Hutchinson probes in the latent: unbiased for the K-mode sum | M transports per source-layer | relative error about sqrt(2/(M k_eff)) on an incoherent per-neuron quantity | no better than E3 |

**Conclusion.** The frame produces no competitive deterministic component. Its positive content is the theorems of
sections 2-4 and the decisive measurement below.

## 6. What this says about a Phase-2 design

1. **Where the boundary sits.** The product-gate defect is the boundary between what CP-leg representations can
   carry cheaply (V0 and V1 moves) and what they cannot (V2). The pinning frame identifies it as the first
   conditioned-versus-frozen excursion term, G o T. Proposition 6 with Conjecture 7 prices its cheapest CP carrier
   at about kappa_90 n transports per source-layer: extensive, 0.05-0.11 n at layers 10-14.
2. **Requirement (i): a separating localization.** A representation in which the gate acts exactly on the carried
   object (the conditional-modular note's question) must include a separating localization of dimension at least
   K_f.
3. **Requirement (ii): or carry no third-order object through the gates.** For example: Gaussian mixtures whose
   components are pushed exactly, with the third order generated only by mixing. That fails by the CP-rank argument:
   K components give kappa3 of rank at most O(K), while the truth needs a flat, extensive tensor (notes XV and XL,
   strong minimality).
4. **The input radius is exactly gate-neutral, and small.** By 0-homogeneity, the input-radius sigma-algebra
   sigma(|x_0|) is independent of the whole pattern sigma-algebra Pat_{L-1}: gates are functions of x/|x|, which is
   independent of |x| under N(0, I).
   - So the input-radial localization commutes with the entire gate algebra: an exact commuting square.
   - It carries only the input part of the gain: Var(|x_0|^2)/E^2 = 2/n, about 0.002, against a layer-norm variance
     of 0.03 at depth (note XIII E2/E10). The network-generated gain is angular and couples to the gates.
   - So the one exactly gate-neutral localization is the one that does not matter, and the separating common mode is
     angular and extensive.
5. **The structural conclusion.** No localization in the network's own gate filtration makes the gates exactly
   independent at sub-extensive cost, apart from the input radius, which is gate-neutral and already exploited
   through homogeneity. The remaining defect is a property of the He-critical spectrum (Fuss-Catalan bulk), not of a
   hidden variable. A system that beats the leaders by accuracy at no more than 0.2 B has two options:
   - (a) remove the defect at the source, with an expansion whose reference makes the third-order transport exact
     (this frame found no such reference among localizations);
   - (b) accept it and win on cost, the leaders' apparent route.

## 7. Predictions

- **P-F2-1 (scaling, local; held).** K_90(G_off) at fixed layer doubles with n. Measured 15 / 27 / 54 at layer 10
  and 11 / 22 / 43 at layer 14 for n = 128 / 256 / 512. The prediction for n = 1024, stated before that run
  finished, was 95-120 at layer 10 and 80-95 at layer 14, with lambda_max of the gate correlation below 12 at every
  layer. Measured: 114, 80 and at most 10.9.
- **P-F2-1b (official networks).** On official networks 0-3, the rank-K masses of rho_off at layers 6, 10 and 14
  match those of fresh He networks of the same width to within 0.02, and K_90(G_off) is 105-125 at layer 10.
- **P-F2-2 (capture law, checked at n = 512).** On real sources, the slope of the rank-K truncated pair-gate readout
  equals the readout-weighted Frobenius mass to within 0.05, for K in {4, 16, 64}. Note XX's n = 1024 slopes exceed
  the unweighted masses by 0.03-0.07 (0.40 against 0.35, 0.74 against 0.67). The proposition attributes the excess
  to the readout weights; the check is to compute the weighted mass on those same dumped legs (a cheap rerun of
  `gatecov_rank.py`).
- **P-F2-3 (the value of the wall: the decisive one).** On the official networks, exact one-step injection of the
  first-order pair-gate correction into the D3 and D21 readouts, at every layer from 1 to 14, with the sources
  themselves unmodified:
  - raw changes by **-12% +- 5%** (networks 0-3, the free-running adopted system);
  - the per-layer one-step residual falls by 25-40% at layers 4-14;
  - the correction's rms is 0.9-1.3% of the product-gate D3 at every layer.

  **Reasoning.** The class is about 6.5e-9 of MSE (closure round). The current raw is 1.55e-8, and V56 addresses a
  different class. One-step injection removes only the k = 0 term of the persistence sum: a fraction
  1/sum_k 0.85^{2k}, which is 0.28 of the class, so 0.28 x 6.5e-9 = 1.8e-9, about -12%.
- **P-F2-4 (persistence).** The pair-gate correction born at layer l, carried as an extra source (through its K = 256
  latent, slope at least 0.94), is read at the D3 of layer l + k with an amplitude ratio of 0.80-0.92 per layer. The
  full class is then worth -35% to -45% raw.

## 8. The minimal decisive experiment (for the AWS fleet)

**Purpose.** Fix the value of the wall, P-F2-3 and P-F2-4. This decides whether any exact-gate representation is
worth the team's search.

**Script outline** (`f2_v2inject.py`, a hook in `est_v29.py` behind `V57_V2INJ=1`):

1. **Run.** The adopted system, free-running, without counterterms, on networks 0-3 (W_off{k}.npy, truth_off{k}.npz),
   paired against itself with the hook off.
2. **The correction.** At each layer l = 1..14, after the sources' legs are transported to the pre-activations of
   layer l, compute for every live source the first-order pair-gate correction of its gated kappa3:
   - Delta kappa3(x_{l+1})_abc = sum over leg pairs of G_ab Phi_c T_abc, with G = D_phi R D_phi from the chain's own
     covariance;
   - optionally with the exact bivariate orthant G, as a control.
   Old-tier sources are expanded to dense legs.
3. **The readouts.**
   - Delta D3_i = sum_abc W_ai W_bi W_ci Delta kappa3_abc;
   - Delta D21_ij = sum_abc W_ai W_bi W_cj Delta kappa3_abc;
   both at layer l+1. Per term r these are the contractions diag(W^T d(phi A_r) R d(phi A_r) W) and
   W^T d(phi A_r) R d(phi P_r) W, plus the q-weighted sums. Add them to the chain's readouts there.
4. **Persistence (network 0 only).** At l = 6 and 10, form the correction's K = 256 latent CP form (the top weighted
   eigenvectors of R) and carry it through 4 further layers with the chain's own leg rules. Record its D3 rms at each
   age.

**Outputs.**
- Per-layer and final raw for both arms.
- The rms of the correction relative to D3 per layer.
- The age profile of the persistence run.

**Compute.**
- Per source-layer: about 4 n^4 multiply-adds, about 9e12 FLOP, about 5-10 s at 1-2 TFLOP/s float32.
- About 120 source-layers, so 15-20 minutes per network on one large CPU instance (c7i.48xlarge class); four
  networks in parallel take about 30 minutes wall-clock.
- Persistence: about 2 x 256 x 4 extra leg transports, about 5 minutes.
- Memory: legs plus one n x n R, streamed per term, within the existing 10 GB footprint.
- No FLOP accounting: the hook is an offline diagnostic.

**Decision rule.**
- **Raw at or below -8% on at least 3 of 4 networks, and persistence at least 0.8.** The class is at least 30% of the
  MSE in full. Combined with Theorem 8 and Conjecture 7, no CP carrier pays. The design search should move to
  representations without the order-3 product-gate defect (section 6, option a) and stop adding closure terms. F2
  hands the problem on, with its value stated.
- **Raw at or above -4%.** The class is at most 15% of the MSE. The leaders' gap lies elsewhere (the kappa4 sector,
  or cost). Drop every joint-gate idea; F2 is closed.
- **In between.** Extend to networks 0-7 before deciding.

## 9. Numerical record (all in `loc/f2/`)

- `n1b_kink2d.py` (outputs `n1b_kink2d_s7.out`, `n1b_kink2d_s11.out`): Theorem 1 exact to 1.6e-12 and 2.5e-12.
  `n1b_kink2d.out` is the n = 8 run that failed (G), with a dead arc.
- `n1_kink.py`, `n1_kink16.py` (`.out`): Monte Carlo of (2) with window extrapolation, and the MF1/MF2 failures.
- `n3_gatecov.py` (`.out`): Proposition 4 to six digits; Proposition 6 at n = 512.
- `n4_sepdim.py` (`n4_a.out` for n = 128 and 256, `n4_b.out` for n = 512 and 1024): K_90, PR and lambda_max. At
  n = 1024, layer 10: PR 90.6; rho_off masses 0.199 / 0.350 / 0.669 / 0.913 at K = 1 / 4 / 16 / 64; K_90(rho_off) = 55;
  K_90(G_off) = 114; lambda_max = 9.0.
- `n6_fusscatalan.py` (`n6_fc512.out`): PR = n/(l+1) for Gaussian products.

## 10. Literature used

- Anari, Koehler and Vuong, "Trickle-down in localization schemes" (arXiv 2407.16104): the trickle-down equation.
- Chen and Eldan (2022): localization schemes.
- Eldan (2013): stochastic localization.
- Oppenheim (2018): local-to-global spectral expansion.
- Zhang, Naitzat and Lim (2018): ReLU networks as tropical rational maps.
- Tsirelson (1985) and Sudakov: E sup_K <g, x> = V_1(K)/sqrt(2 pi).
- Hanin and Rolnick (2019); Özkan and Hirsch (arXiv 2609.32695): co-area and Kac-Rice for the kink set.
- Plackett (1954) and Mehler: orthant derivatives and expansion.
- Fuss-Catalan spectra of Gaussian products: Banica, Belinschi, Capitaine and Collins; Neuschel.
- Repo notes XIII (ncg-probability), XVIII-XX (closure round, hidden separator), conditional-modular, XL, XLI.

## Referee report

Adversarial referee, 9 October 2026. My checks are in `loc/f2/ref/` (`r1_tropical.py`, `r2_kink.py`,
`r3_freediag.py`, `r3b_psd.py`, with outputs `*.out`). All were run locally on fresh He networks; none used an
official network.

### R0. Verdict in one paragraph

The negative verdict survives, and it is stronger than the author's own arithmetic makes it. The exact identities
are correct; I re-derived Theorem 1 and Proposition 2 independently, to 1e-10 and 2e-13. The headline "main result",
however, is built on a mis-defined quantity. K_90 in the gate-covariance metric (114 at n = 1024) is measured with
the diagonal forced to zero. A separating latent leaves the diagonal free, and with the diagonal free the same
matrices need 2.5-2.7x fewer directions. Even so, the separating dimension plays no part in the score decision, which
holds at any K by a simple break-even bound (R4). The cost table also overstates what each K removes, by about 1.6x.
The decisive experiment is sensible, but its baseline does not match its prediction. The noncommutative structure
(Prop. 1) is, as the author admits, trivial at finite depth. The tropical identity is a pleasant corollary of known
results. The frame offers no route to the 1e-8 level.

### R1. Identities re-checked

| claim | my check | result |
|---|---|---|
| Theorem 1, neuron form (2) | `r2_kink.py`: input dimension 2, widths 16/16/12, depth 4, 3 seeds; my own root-finding of kink rays and my own upstream Jacobian and downstream G at each ray | max \|E F - sum\| = 6e-11, 5e-11, 1.0e-10 (82-162 rays). The residual is quadrature error. **Holds.** |
| Prop. 2, E F_i = (V_1(P) - V_1(Q))/sqrt(2 pi) | `r1_tropical.py`: n = 2, widths 5/5/4, polygons built by explicit Minkowski sums and conv(A ∪ B), V_1 = perimeter/2 | max diff 2.3e-13; pointwise h_P - h_Q = F to 6e-13. **Holds.** |
| Prop. 3 tilt, trickle-down, random-order pinning | algebra | Correct: law of total covariance with binary s. |
| Prop. 4 identity, first-order Mehler | algebra; the author's 6-digit table | Correct. |
| Prop. 6, E M[w,w,w]M'[w,w,w] = 6(2/n)^3 <Sym M, Sym M'> | Isserlis | Correct for all-distinct supports. |
| Input radius independent of Pat_{L-1} | 0-homogeneity of the gates; \|x\| independent of x/\|x\| | Correct; 2/n is right. |

Theorem 1 deserves a lower billing. It is Gaussian integration by parts, E F = E[<x, grad F>] = E[Delta F], applied
to a piecewise-linear F. The content is the bookkeeping on the facets, and that bookkeeping is correct. "The mean
lives on the codimension-one skeleton" is Stein's identity, not a new structural fact.

### R2. Mathematical errors and mis-statements

1. **The capacity corollary (Prop. 4) is wrong as stated, and the separating dimension is measured with the wrong
   truncation.**
   - A K-dimensional Gaussian separating latent means rho = D' + B B^T. So D_phi rho D_phi is
     *diagonal + rank K*, not "rank at most K".
   - G_off agrees off the diagonal with a rank-K matrix, and the diagonal of that matrix is free: the pair-gate
     correction lives on all-distinct triples, so the diagonal never enters it.
   - The best K-latent capture is therefore the free-diagonal (Frisch / factor-analysis) problem,
     min over rank-K PSD L of ||offdiag(G - L)||_F. This is not the top-K SVD of the zero-diagonal G_off.
   - The corollary's "captures *at most* the top-K Frobenius mass of G_off" is false: the free-diagonal fit captures
     strictly more. Measured on the author's own seeds and Monte Carlo sizes (16384 inputs), with the PSD constraint
     making no difference:

   | n | layer | metric | zero-diagonal K_90 (author) | free-diagonal K_90 | mass at K = 16, zero / free |
   |---|---|---|---|---|---|
   | 128 | 10 | G_off | 15 | 6 | 0.912 / 0.988 |
   | 256 | 10 | G_off | 27 | 10 | 0.859 / 0.961 |
   | 512 | 10 | G_off | 54 | 21 | 0.791 / 0.853 |
   | 128 | 14 | G_off | 11 | 5 | 0.930 / 0.992 |
   | 256 | 14 | G_off | 22 | 8 | 0.873 / 0.977 |
| 512 | 14 | G_off | 43 | 17 | 0.820 / 0.898 |
   | 256 | 10 | rho_off | 14 | 11 | 0.913 / 0.959 |
   | 512 | 10 | rho_off | 24 | 19 | 0.855 / 0.876 |

   - In the gate metric the author's headline K_90 is inflated about 2.5-2.7x. The phi-weighting makes the removed
     diagonal very non-constant, and that spreads mass into the tail.
   - The extrapolation "K_90(G_off) = 114 at n = 1024, layer 10" should read about 40-45. That figure assumes the
     free-diagonal K_90/n of about 0.04, which is still roughly linear: 6, 10, 21 over n = 128-512.
   - Extensivity (Conjecture 7) may survive with kappa_90 of about 0.04 rather than 0.11. P-F2-1 "held" only for
     the mis-specified quantity.
   - Note XX's rank-K slopes (0.40 / 0.74 at K = 4 / 16) also used zero-diagonal truncation. They are therefore
     lower bounds on what a K-latent captures.

2. **The cost table (Theorem 8c) uses the wrong removal formula.**
   - "Capturing slope s removes (1 - (1-s)^2)/3" would hold only if the truncated readout were s times the exact
     one. It is not: its correlation is about sqrt(s).
   - For a projection-like truncation, the residual is |c - c_K|^2 = (1 - s)|c|^2, so the MSE removed is about s/3.
   - Note XX itself applied its formula with s = 1 - (residual rms), not with the slope. Its residual rms is 77% at
     K = 4 and 50% at K = 16, which gives 13.6% and 25% of the MSE removed.
   - F2 instead plugs the slope into the formula and gets 21% and 31%. It overstates the gain about 1.6x at K = 4.
   - Corrected adjusted ratios: K = 4 gives 0.864 x 0.36/0.203 = **x1.53**; K = 16 gives 0.75 x 0.83/0.203 = **x3.07**.
     This strengthens the negative.

3. **"Every complete localization refines the pattern localization" (Prop. 1(iii)) is true only at terminal time.**
   - At finite t, Eldan's sigma-algebra F_t does not contain Pat_l, and the two filtrations are not nested.
   - "Refines" here means only that a point determines its pattern. It is not a statement about the processes.
   - The rest of Prop. 1 is correct but empty: finitely many classes of positive measure, a type-I algebra, finite
     sums. The author says so.

4. **Theorem 8(a) is n^4 only when a source's legs have n terms.**
   - For a source with R_s CP terms, the exact D3 correction costs R_s units per source-layer: one product
     G (diag(a_r) W) per term.
   - The K_4 / treewidth-3 statement is correct for dense r. "Exact = 1.2e5 units" implicitly sets R_s of about n
     for every source.
   - The verdict is unaffected, since even sum_s R_s over all source-layers is far above budget. The statement should
     nonetheless carry its R-dependence.

5. **P-F2-2 says "checked on real source legs at n = 512".** `n3_gatecov.py` uses a synthetic 48-term source
   (covariance-row arms plus a dense centre). The check on real legs is the note-XX n = 1024 table, and there the
   weighted masses were never computed. The capture law is verified on a synthetic source only.

### R3. Hidden assumptions

- **The value of the class (about 1/3 of the MSE, 6.5e-9) rests on one measurement.** That is one network and one
  layer (closure round, layer 10, 48 neurons), scaled by an assumed persistence of 0.85 from note XVII. Every
  quantitative system statement in F2 inherits this number.
- **The baseline of P-F2-3 does not match its prediction.**
  - The -12% is 1.8e-9 / 1.55e-8, with 1.55e-8 being the adopted system *with* counterterms.
  - The experiment's arms run *without* counterterms, where raw is larger (the closure-round figure was 2.2e-8).
    The same absolute removal then reads about -8%, which falls exactly on the decision threshold.
  - The prediction should be stated in absolute MSE (-1.8e-9 +- 0.8e-9), or recomputed on the no-counterterm
    baseline.
- **The D21 part was never in the 6.5e-9.** Adding the D21 injection measures a larger class than the one predicted.
- **The decomposition assumes the class's error is uncorrelated with the rest of the error.** The age-0 removal
  fraction 1 - 0.85^2 = 0.28 is the same for coherent and incoherent ages, so that part is robust. The reduction
  equals |e_0|^2, however, only if <E - e_0, e_0> = 0; the cross term can be of either sign. With four networks this
  is within the ±5% band only if the cross term is small.

### R4. Relevance: the decision does not need Conjecture 7

There is a two-line bound that settles the design question without the separating dimension.
- The pair-gate class is at most 1/3 of the MSE.
- An addition of cost Delta C pays only if the fraction of MSE it removes exceeds Delta C/(C + Delta C).
- So nothing that targets this class can pay unless Delta C < C/2, about 104 units, whatever its representation.
- At about 40 units per latent direction that allows K <= 2. At K = 1-2 even the free-diagonal capture is about
  0.2-0.35 in the gate metric, and the readout-weighted capture is lower, since the mean-direction spike holds 2%.
- The best case removes about 3-6% of the MSE for +20-40% flops.

The extensivity of K_f matters only to n-scaling, and n is fixed at 1024. "Main result" overstates the role of
Conjecture 7 in the decision. The real constraint is break-even against C/2, and that holds for every K.

The same bound kills any sampling route (E5, E6): a 1%-of-D3 incoherent correction with zero annealed mean needs a
relative noise far below 1% per neuron per layer.

At the 1e-8 level the class is real, at about 6.5e-9 of the 1.55e-8. But nothing in F2 can recover it inside the
budget, and both branches of the decision rule end in "build nothing from F2". The experiment is value-of-information
for the team's direction (it closes the joint-gate line or sends the work to defect-free representations). It is not
a route to a score.

### R5. Novelty

- **Theorem 1.** Gaussian integration by parts plus facet bookkeeping. It is the jump-weighted analogue of the
  Hanin-Rolnick co-area count. Correct and tidy; low novelty.
- **Proposition 2.**
  - Zhang-Naitzat-Lim (2018) give the bias-free tropical recursion. Tsirelson gives V_1 = sqrt(2 pi) E h_K.
  - I found no published statement that E F_i = (V_1(P) - V_1(Q))/sqrt(2 pi) for a ReLU network (Exa search; the
    tropical-NN literature is about region counts and decision boundaries). It may be new as a statement.
  - It is a one-line corollary of those two results and has no algorithmic content (the polytopes have astronomically
    many edges).
  - Its one structural lesson, that linear layers are exactly V_1-additive and all difficulty sits in conv(A ∪ B),
    is the standard "the nonlinearity is the gate" observation in new clothes.
- **Props. 3-5 and the cluster/Mehler hierarchy.** Standard (total covariance, Schur product, moment-cumulant,
  tetrachoric/Mehler). They relabel the repo's V1/V2 classes. AKV, Gelfand-Tsetlin and Theorem E are names attached
  to one-step identities; none of them adds a computation.
- **Prop. 6.** A small, genuinely useful lemma: zero annealed mean, slope = weighted mass. It explains note XX's
  table. Note XX had already measured the slope-mass tracking.
- **Section 3.4.** The author is correct that the spectral-independence and AKV machinery bounds variance and mixing,
  not closure bias, so it does not bear on a deterministic estimator. That is worth having said.

### R6. Cost accounting

- The unit conventions are right: 1 unit = 2n^3, and n^4 MACs per source-layer = 1024 units.
- The experiment's compute (about 9e12 FLOP per source-layer, about 1e15 per network) is consistent.
- Two corrections are needed:
  - the R_s-dependence of the exact cost (R2.4);
  - the removal arithmetic (R2.2).
- E1's "16 units" for the kink mean field is plausible, and irrelevant since that design fails at O(1).

### R7. What to change in the text

1. Restate the corollary of Prop. 4 with a free diagonal: "G_off = offdiag(L), rank L <= K". Redefine K_f as the
   free-diagonal (factor-analytic) quantity, and remeasure K_90 at n = 1024 that way. I predict about 40-45 at
   layer 10 and about 30-35 at layer 14 in the gate metric (free-diagonal K_90 is 5 / 8 / 17 at layer 14 for
   n = 128 / 256 / 512).
2. Replace (1 - (1-s)^2)/3 with 1 - (residual)^2, then divide by 3. Use note XX's residual column.
3. State P-F2-3 in absolute MSE and on the experiment's own (no-counterterm) baseline. Measure the D3-only and the
   D21 injections separately, so that the 6.5e-9 is actually tested.
4. Demote Conjecture 7 from "main result" to a scaling remark. Lead the system section with the C/2 break-even bound.
5. Relabel the P-F2-2 n = 512 check as synthetic-source.

### R8. Strongest salvageable idea

A *value-of-information bound* for the whole pair-gate (V2) class. Its D3 readout has zero annealed mean over the
fresh weights (Prop. 6, proved), so no calibrated counterterm can absorb it. Its value is at most about 1/3 of the
current MSE. Any carrier must therefore cost less than C/2, about 100 units in all, to pay.

The free-diagonal (factor-analytic) separating dimension, about 0.04 n per layer rather than 0.11 n, sets how close a
latent carrier can come to that bound. At 40 units per direction it cannot come close.

The offline one-step injection experiment, re-baselined, is the right one-time measurement to close the joint-gate
direction. Its useful outcome is a redirect, not an estimator.
