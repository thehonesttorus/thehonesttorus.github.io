# Activation modules read against the programme

Working note XVI. Two texts by George Jeffreys (with Siu-Cheong Lau): the paper "Noncommutative geometry of
computational models and uniformization for framed quiver varieties" (arXiv 2201.05900) and the Boston University
dissertation "Noncommutative geometry for computing machines" (2022), which compiles that paper with the
quantum-finite-automata paper and the Kaehler/toric paper. The dissertation was read in full. This note records what
the construction is, what the dissertation adds beyond the paper, which pieces are the same objects as ours under
other names, and what the construction leaves out. Nothing here changes the chain; section 5 says what we keep.

## 1. The construction in one paragraph

A computing machine is a framed quiver representation: vertices carry spaces V_i, arrows carry linear maps w_a,
and every vertex has a framing e^(i): C^{n_i} -> V_i. The algebra A of paths is augmented by formal symbols
varsigma_j for the nonlinear activations, giving a near-ring A{varsigma} (a ring without left distributivity:
x(y + z) is not xy + xz). Elements are activation trees: paths with activation nodes. A Karoubi-de Rham differential
on the near-ring satisfies d^2 = 0; the differential of a 0-form is the sum over nodes of chain-rule paths, which is
backpropagation written without choosing dimensions (paper Prop 2.36, dissertation Prop 3.4.14). The near-ring forms
descend to GL(V)-invariant forms on the representation space (Theorem 2.40 / 3.4.18). The moduli of stable framed
representations carries a canonical GL-equivariant metric on each universal bundle,
H_i = (rho_i rho_i^*)^{-1}, rho_i = (w_gamma e^{(t(gamma))})_{paths gamma into i}, with the recursion
H_i^{-1} = e e^* + sum_{h(a)=i} w_a H_{t(a)}^{-1} w_a^* (Lemma 3.18 / 4.2.14, Prop 4.1.8). Activations are applied
to H-whitened coordinates on the framing. Changing the sign of the framing's quadratic form uniformizes the moduli
into compact, Euclidean (ordinary weight space, Remark 3.26 / 4.2.22) and non-compact (space-like Grassmannian,
Theorem 3.17 / 4.2.13) types.

## 2. What the dissertation adds beyond the paper

- **Motivation from the measurement problem (Chapter 2, section 3.3).** A quantum finite automaton has unitary
  state changes and one nonlinear operation, the probabilistic projection sigma_0 onto a basis state with Born
  probabilities. The near-ring is proposed as the algebraic home for linear evolution plus measurement, with a
  pointer to Pascual Jordan's abandoned near-ring attempts. Section 3.3 writes the automaton as a unitary activation
  module: the measurement is the activation, iterated reads give a function from initial to final states.
- **The metric from the moment map (section 4.1).** The moduli is the symplectic quotient of the moment map
  mu_i = e e^* - sum_{t(a)=i} w_a^* w_a + sum_{h(a)=i} w_a w_a^* at level -I. Theorem 4.1.7 proves H_i is positive
  definite on the level set by induction over a topological order of the vertices; Theorem 4.1.14 shows it is the
  metric of the iterated Grassmann bundle (each vertex is a Grassmannian of quotients of framing plus incoming fibres
  over the moduli of the sub-quiver of earlier vertices). The residual symmetry of the metric is the product of
  unitary groups of the framings.
- **Activations from toric geometry (Chapter 5).** The sigmoid is the moment map of P^1 on its open orbit; the
  Guillemin-Abreu symplectic coordinates give, for P^d, the map z -> z / sqrt(1 + |z|^2), a U(d)-equivariant
  symplectomorphism of C^d onto the unit ball (Prop 5.1.4, Lemma 5.1.6). Any projective toric variety gives
  sigma(r) = sum_i e^{2(u_i, r)} u_i / sum_j e^{2(u_j, r)} with u_i the polytope's vertices; softplus is a toric
  Kaehler potential on C (Example 5.1.8). The family version over the moduli is v -> v / sqrt(1 + H_i(v, v))
  (Prop 5.1.9), equivariant under the residual unitary symmetry, and it is also the symplectomorphism of C^n onto the
  hyperbolic ball (Prop 5.1.13). Example 5.2 writes the three-vertex network explicitly: H_1 = I,
  H_2 = (I + b b^* + W_1 W_1^*)^{-1}, H_3 = (I + W_2 W_2^* + W_2 b b^* W_2^* + W_2 W_1 W_1^* W_2^*)^{-1}, and notes
  that for small weights the machine function reduces to the flat W_2 sigma(W_1 v + b).
- **Universal approximation (section 5.3).** Theorem 5.3.2: the flat machine with the P^{d_2} activation is
  L^2-dense on compacta. The proof is Cybenko's rescaling trick in toric language: sigma_t(x) = sigma(t x) converges
  as t -> infinity to a step function constant on the cones of the dual fan (Lemma 5.3.3, the tropical limit), a
  centered simplicial web is the preimage of the fan of P^d under a linear map (Theorem 5.3.8), so step functions on
  webs are reachable and are dense. Remark 5.3.7 notes the webs have no integral structure and no balancing
  condition because the weights are real.

## 3. The same objects under other names

- **Near-ring = the chain's term algebra.** Activation trees with nodes D^{(p)} varsigma are the chain's Hermite
  jets of the gate; the generator of the term tables in the ncg-probability note (composition and substitution
  Hopf algebras, Faa di Bruno) is the near-ring's multiplication table. What the near-ring lacks is a state: an
  expectation functional over the random variables that the activation nodes act on. In the automaton chapter the
  state is explicit (the probability space Omega with the Born rule) and then dropped when the near-ring is formed;
  from then on nothing is averaged. The chain is the computation of near-ring elements under the Gaussian state of
  He-initialized weights and inputs; the gate's statistics (Phi, phi, the pair coefficients) are that state.
- **Moment-map level set = annealed weights.** For He initialization with no bias, mu_l = W_{l+1}^* W_{l+1} -
  W_l W_l^* has expectation zero and fluctuations of order n^{-1/2}. The dissertation's symplectic reduction fixes the
  level exactly and so works in the annealed approximation; E10 showed the gain lives in the quenched fluctuations
  (the mean-coupled third cumulant). The construction therefore describes the part of the problem the leaders and we
  already have, and is silent on the part that sets the score.
- **Metric recursion = linear covariance propagation.** H_i^{-1} = e e^* + sum w_a H_{t(a)}^{-1} w_a^* is the
  covariance of the vertex signal when all framings are fed white noise, with no gates. The chain's covariance step
  is this recursion with the gate's pair coefficients inserted; the dissertation never inserts them because the
  activation is applied after the metric is formed, never inside it.
- **Iterated Grassmann bundle = birth and transport of sources.** Theorem 4.1.14's induction pulls back the
  metric of vertex j from the moduli of the sub-quiver of vertices up to j; later vertices act only by pullback.
  That is the statement that a source is fixed at its birth layer and afterwards only transported (note XV, the
  propagator nesting P_{b'}(l) = P_b(l) P_{b'}(b)). The shared-basis tier of v29 is a numerical range finder for the
  fibres of this bundle.
- **Tropical limit = the gate.** ReLU is the tropical limit of softplus, softplus is the toric potential of C, so
  the challenge's activation sits exactly at the degenerate point t -> infinity of the dissertation's family, where
  the Kaehler potential becomes piecewise linear and the symplectic form concentrates on the walls. The Gaussian
  measure regularizes this point: Phi(alpha) is the Gaussian average of the step, phi(alpha) the average of the
  wall's delta function. The dissertation uses the limit only to prove approximation; our whole apparatus is the
  measure at that limit.
- **Toric activations are radial compactifications.** v / sqrt(1 + H(v, v)) is bounded and contracts along the
  mean direction. Under such an activation the gain mode of E10 (top eigenvector = mean direction, growing with
  depth) would be squashed; the challenge's homogeneous gate is the opposite regime, where the mean direction is the
  unstable direction of the Lyapunov transport.

## 4. What the construction leaves out

- **No state, no expectation.** As above. Without a functional on the near-ring, nothing can be computed: no
  kernel, no Mehler formula, no cumulant, no gain. The GL-invariants of the representation space are cyclic words
  in weights and framings; their evaluation under a Gaussian state is free probability (fresh-weight asymptotic
  freeness of the layers), which is absent.
- **Inputs are never integrated.** The framing e is a fixed map. Gaussian inputs make e e^* a covariance; the
  no-bias challenge makes the first moment-map equation read "input covariance equals identity". The data enters the
  dissertation only through the training loss, never through the geometry.
- **Homogeneity is unused.** The toric family is parametrized by the moment-map level, which is a scale. ReLU has
  no scale, and that symmetry is what gives the exact arc-cosine pair coefficients and the dilation-group description
  of the depth dynamics (note XIV). The dissertation's geometry degenerates at exactly the scale-free point.
- **The approximation theorem does not use the geometry.** The approximants have weights t L with t -> infinity,
  where H_2 = (I + b b^* + W_1 W_1^*)^{-1} -> 0. The metric and symplectic structure play no role in the proof,
  the norm is L^2 rather than uniform, and no rate is given. The universal approximation of the canonical (curved)
  machine function is not proved.
- **Missed connections.** The moment-map level set is the balanced manifold of deep linear networks: the quantities
  W_{l+1}^* W_{l+1} - W_l W_l^* are conserved by gradient flow of a GL-invariant loss (Kempf-Ness, as in the
  implicit-bias literature on deep linear networks), so the dissertation's moduli is where those networks already
  live, and the implicit-bias results are dynamics on its fibres. The tropical hypersurface arrangement of a ReLU
  layer (tropical geometry of neural networks) is the dissertation's Lemma 5.3.3 read as geometry of the network
  rather than as a proof device. The information metric of a network on data (Fisher, natural gradient) occupies
  the slot of the tangent metric H_T of Example 5.2, defined here from weights alone.
- **Oriented cycles.** Remark 4.1.10: the recursion becomes an infinite sum for recurrent quivers and is only
  treated where it converges. The gated version of that sum is the Lyapunov transport, which is where we have the
  convergence theory (note XIV).

## 5. What we keep

- The correct algebraic statement of our term tables: near-ring elements (activation trees) evaluated under a
  state. Writing the chain's generator as a near-ring makes explicit that the only choices are the state and the
  truncation, which is the view that produced E10 and E11.
- The identification of our regime: the tropical point of the toric family, with the Gaussian measure as the
  regularization. It explains why nothing in the symplectic construction can be transferred as a computation: the
  structure degenerates exactly where the challenge lives.
- The confirmation, from the Grassmann-bundle induction, that the age structure (fix at birth, transport by pullback)
  is exact and that the only lossy steps are the range finders at the joins (note XV, section 3).
- A negative result worth stating: the moment-map / whitened coordinates are not the right basis for the old tier.
  The loss is final-layer MSE in raw coordinates, so raw energy, not H-whitened energy, ranks directions; the
  residual unitary symmetry of H is broken by the loss.
