# The gated transport of cumulants at the independent reference

Working note XXXIV. A theorem and its closed forms, answering obligation 1 of note XXXIII section 6 (the response
cores at alpha != 0) and sharpening the closure statement of notes XXXII and XXXIII. Derived by Gaussian
integration by parts and checked against direct quadrature in `code/check_gated_transport.py` (16 closed forms and
three vanishing cases to rounding, including the cumulant theory's N(0, I) witness). No network experiment.

## 1. Setting

The reference is the product law z_a ~ N(mu_a, s_a^2), independent, with y_a = relu(z_a), alpha = mu/s, Phi = Phi(alpha),
phi = phi(alpha), m = E y = mu Phi + s phi, v = Var y. A first-order tangent of the law of z is a set of cumulant entries;
an entry on index support S with multiplicities k_a >= 1 (sum k_a = r) has Hermite score

    omega * prod_(a in S) He_(k_a)(q_a),   q_a = (z_a - mu_a)/s_a^2,   He_1 = q, He_2 = q^2 - 1/s^2, He_3 = q^3 - 3q/s^2,

with omega the entry times its number of orderings over r!: C_ab for a covariance entry, Gamma_ab/2 for a third-cumulant
entry kappa(z_a, z_a, z_b) = D21_ab, T_abc for an all-distinct third-cumulant entry, K_ab/4, B_ab/6 and V/2 for
fourth-cumulant entries of the (2,2), (3,1) and (2,1,1) patterns, U for an all-distinct one. (C_off is treated as a
tangent, so every statement is first order in the off-diagonal covariance as well.)

The univariate gate coefficients are c(p, k) = E[(d/dz)^k (y - m)^p] with distributional derivatives (Price's theorem for
nonsmooth tests):

| | k = 1 | k = 2 | k = 3 |
|---|---|---|---|
| p = 1 | Phi | phi/s | -alpha phi/s^2 |
| p = 2 | 2m(1 - Phi) | 2 Phi - 2 m phi/s | |
| p = 3 | 3(v - m^2 (1 - Phi)) | 6m(1 - Phi) + 3 m^2 phi/s | 6 Phi - 6 m phi/s - 3 m^2 alpha phi/s^2 |

The terms in phi/s are the gate facets: the boundary density of z at the hinge, entering through delta(z).

## 2. The theorem

**Support preservation.** The first-order change of a joint cumulant of y with index support S' vanishes unless the
tangent entry has support exactly S'. A coordinate outside the entry's support carries the constant score and stays
independent of everything the tangent moves, so every cumulant involving it keeps the value zero to first order; a
coordinate inside the entry's support but outside S' contributes a factor E[He_k] = 0. In particular a unary source (the
diagonal of the third or fourth cumulant, the variance) moves no mixed cumulant, and a pair source moves no cumulant on
three or more neurons.

**Gated transport.** For S' = S, the change is omega times a product of univariate coefficients. For the fourth-order
classes:

    (3,1)    kappa(y_a, y_a, y_a, y_b):   omega [c_a(3, k_a) - 3 v_a c_a(1, k_a)] c_b(1, k_b)
    (2,2)    kappa(y_a, y_a, y_b, y_b):   omega c_a(2, k_a) c_b(2, k_b)
    (2,1,1)  kappa(y_a, y_a, y_b, y_c):   omega c_a(2, k_a) c_b(1, k_b) c_c(1, k_c)
    (1,1,1,1) kappa(y_a, y_b, y_c, y_d):  omega prod c(1, k) = omega Phi_a Phi_b Phi_c Phi_d   (only k = 1 there)

Proof. The score factorises over coordinates and the reference is a product law, so every raw moment's first variation
factorises; Stein's identity E[He_k(q) f(z)] = E[f^(k)(z)] converts each factor to a gate coefficient. In the partition
expansion of the cumulant, every product term containing a moment that involves a coordinate of S with nonzero k but not
the full product vanishes at first order (it either has a zero-mean score factor or a reference cross-moment
E[(y_a - m_a)(y_b - m_b)] = 0); only the full term survives, together with, for the (3,1) pattern, the term
-3 E[(y_a - m_a)^2] E[(y_a - m_a)(y_b - m_b)], whose second factor varies. That gives the bracket. QED.

The transport is support-preserving but not multiplicity-preserving: a y entry of pattern (2,1,1) on {a, b, c} receives the
z entries on {a, b, c} of patterns (2,1,1), (1,2,1) and (1,1,2), each with its own coefficient product.

The cumulant theory's witness is the special case alpha = 0, s = 1, Gamma_12 = t: c(1) = 1/sqrt(2 pi) gives
kappa(y1, y1, y1, y2) = t(3c/8 + 3c^3/2) and kappa(y2, y2, y2, y1) = t(3c/8 - 3c^3/2). The formulas above extend it to every
alpha, every scale, and every pattern of order four.

## 3. Consequences for the chain

**The pair-supported sector is closed under the gate.** The y classes (4), (2,2), (3,1) respond only to pair-supported
z entries: the off-diagonal covariance, D3, D21, and the z fourth-cumulant diagonal, (2,2) and (3,1) slices, all of which
the chain holds. The pair program computes (4) and (2,2). It never computes the y (3,1) slice, which is therefore a pure
representation omission (note XXXII section 3), now with a closed form. In the chain's variables (Gamma = D21,
K^z = wk4m, B^z = wk431^T, C = C_off), to first order around the independent reference:

    B^y_ab = 3 (v_a - m_a^2)(1 - Phi_a) Phi_b C_ab
           + (Gamma_ab / 2) [6 m_a (1 - Phi_a) + 3 (m_a^2 - v_a) phi_a/s_a] Phi_b
           + (Gamma_ba / 2) 3 (v_a - m_a^2)(1 - Phi_a) phi_b/s_b
           + (K^z_ab / 4) [c_a(3,2) - 3 v_a phi_a/s_a] phi_b/s_b
           + (B^z_ab / 6) [c_a(3,3) + 3 v_a alpha_a phi_a/s_a^2] Phi_b
           + (B^z_ba / 6) 3 (v_a - m_a^2)(1 - Phi_a) (-alpha_b phi_b/s_b^2).

Every term is a per-neuron coefficient times an n x n matrix the chain already has: O(n^2) elementwise. The first term is
the one-loop (3,1) class, of the scale mixture's shape (a per-neuron factor times C_ab); the next two are the
third-chaos facet response of the D21 source, which note XXXIII put at leading order; the last three are the gated
transport of the z fourth-cumulant slices.

**The non-closure is the linear layer's.** Since the gate preserves supports, the pair-restricted state would be exactly
closed if the next layer read only pair-supported y entries. It does not: z' = W y mixes supports, and the next pair
slices and diagonal receive the y (2,1,1) and (1,1,1,1) classes through contractions such as
6 sum W_ia^2 W_ib W_ic kappa(y_a, y_a, y_b, y_c). By the theorem those classes have exactly three first-order sources:

1. the all-distinct third-cumulant entries, kappa(y_a, y_a, y_b, y_c) gaining 2 m_a (1 - Phi_a) Phi_b Phi_c T_abc. These live
   in the chain's legs, not in its pair state;
2. the z (2,1,1) and (1,1,1,1) classes themselves, gated: memory of the omitted classes;
3. at second order, products of pair tangents (C_ab C_ac and their mixed terms with D21): the scale mixture's dropped
   classes, 6 g sigma_diag^2 sigma_off^2 on the diagonal (note XXI), belong here.

This is the precise form of E_l for the chain's fourth-cumulant state: the triple and quadruple classes, fed by the legs,
by their own memory and by second-order products of the covariance. Note XXXII's F2 (the defect is regenerated at each
step) says the memory term, source 2, carries little; sources 1 and 3 are what the fitted lam C_off term stands in for.

**What the readout of the triple class costs.** If the triple entries factorise as (a per-neuron vector at a) times (a
matrix in b, c), as the scale mixture's and the leg-fed ones do (2 m_a (1 - Phi_a) times the gated leg contraction), their
contribution to the next diagonal is 6 (H u)_i (w_i^T M w_i)_off: O(n^2) once M is a matrix whose transport is already
computed (the covariance), or a readout of the transported legs with one leg replaced by W^2 (u o alpha_j), the same shape
as the chain's D3 readout. The gated transport does not preserve that factorised form exactly (the multiplicity
patterns on one support are re-weighted), so whether a few product cores carry the triple class to the required
precision is open; SSC's success suggests they can, and the theorem says what they must reproduce.

## 4. What this changes in the plan

- The y (3,1) output is a representation repair with a closed form at O(n^2): no assumption, no new observable.
- The transport of the full pair-supported y tensor (diagonal, (2,2), (3,1)) to the next layer's slices is then the
  coherent tensor of note XXXII section 2, now including the (3,1) class it lacked.
- The triple and quadruple classes are the law-closure frontier. Their first-order sources are named: the legs'
  all-distinct entries and second-order covariance products. Both have explicit gate coefficients.
- The decisive kernel test (note XXXII section 6, step 1) now has a sharper form. From the one-step dump, evaluate the
  exact transport of D(d) + S22(K) + S31(B^y), with B^y from the formula above, in all three slices, and compare with Monte
  Carlo truth. What remains is the triple and quadruple classes plus the Wick truncation and the correlated-reference
  corrections. That remainder decides whether product cores for the triple class are needed, or whether the
  pair-supported repair suffices.

## Code

- `code/check_gated_transport.py`: direct quadrature of the first-order change of each post-activation cumulant under
  each tangent (scores integrated against powers of relu, cumulants by the partition formula, no Stein identity),
  against the closed forms; support-preservation zeros; the N(0, I) witness.
