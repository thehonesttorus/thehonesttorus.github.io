# Flexes of the pair state: the omitted classes as hidden coordinates

Working note XXXVIII. This note reads seven papers for their mechanisms and turns the reading into one theorem about the chain, one pre-registered test, and a short list of the links that are only relabellings. Five of the papers come from one research line; the other two are the coding papers of the previous round. The standard is that of notes XVII-XIX: a link counts when it changes a number or a decision.

| arXiv | short name in this note | authors | title |
|---|---|---|---|
| 2510.13777 | subspace designs | Brakensiek, Chen, Dhar, Zhang | *From random to explicit via subspace designs, with applications to local properties and matroids* |
| 2605.02211 | Hamiltonian sparsification | Basu, Brakensiek, Putterman | *Many Hamiltonians are sparsifiable* |
| 2606.09728 | quantum-cut sparsifiers | Basu, Brakensiek, Kothari, Putterman | *Quantum cut sparsifiers* |
| 2609.36761 | covers | Brakensiek, Guruswami, Putterman | *The reach of Abelian covers in hypergraphs* |
| 2603.28700 | multiway cut | Brakensiek, Huang, Potechin, Zwick | *Improved approximation algorithms for multiway cut by large mixtures of new and old rounding schemes* |
| 2608.09347 | rate bounds | Alrabiah, Guruswami | *Binary code rate bounds via classical-quantum channels* (the pretty good criterion and syndrome duality) |
| 2608.00265 | low-degree tests | Alrabiah, Arunachalam, Grewal, Wright | *No low-degree tests for quantum states* |

## 0. What this note establishes

1. **The latent model is independence with amalgamation over a collective algebra.**
   - Conditioning on the collective coordinates t is a projection onto the subalgebra ℬ of functions of t.
   - When the units are conditionally independent given t, every off-diagonal fourth-cumulant class is an evaluation of per-unit symbols on the collective space (Proposition 1). Exact enumeration checks it to 1e-14.
2. **The pair slices are a multiplicity code of those symbols.**
   - The (3,1) row of a unit is a vector symbol g_a.
   - Its (2,1,1) block is a quadratic symbol N_a.
   - The (1,1,1,1) class is one global quartic T4.
   - The (2,2) slice is the symmetrised quadratic symbol plus the covariance of conditional variances.
3. **Flex theorem (Theorem 2).**
   - The diagonal, (2,2) and (3,1) slices determine the omitted (2,1,1) and (1,1,1,1) classes only modulo an exact kernel F of dimension C(K+3, 4) + q(q-1)/2, with q = K(K+1)/2.
   - F is spanned by a free quartic and an antisymmetric rotation of the quadratic symbols. The omitted classes move along F.
4. **Consequences.**
   - Every closure that computes the omitted classes from pair statistics picks a point of F by convention. That includes the lam core, the scale-mode core D and the G D core KD.
   - The truth's point of F is set by the history of the law.
   - The compensation between the fourth-cumulant slices is one missing coordinate seen three times.
   - There are two kinds of invisibility. The null code of note XXXVI is invisible to the readout and can be quotiented out. The flex is invisible to the present state and has to be carried.
5. **The minimal extension and its price.** The state would carry T4 and the per-unit symbols N_a. Their transport into the next layer's slices costs about 1.5-4 units (2n³) per layer at K = 16-32 (section 5).
6. **Test S**, pre-registered in section 6: does the truth's omitted-class transport equal the transport of its own collective symbols? Results are in section 8.
   - The amplitude comes out right with nothing fitted: the coefficient is 0.99-1.08 from 8 -> 9 on, and within 13% everywhere.
   - At depth the share captured reaches 0.59-0.79 of the (3,1) slice at K = 64. This does not reject ℬ-independence.
   - Three registered predictions fail: P1 (narrowly at 8 -> 9, clearly at 2 -> 3), the second half of P2 (convergence in K is slow), and P3 (nothing beyond the scale mode at K = 32).
   - At affordable K the hidden coordinate is the scale mode that the chain already found output-adverse. The carried-symbol extension is closed for cost (section 9).

## 1. The mechanisms, as used here

**Subspace designs (2510.13777).**
- *Definition.* Subspaces H_1, ..., H_n of F^k form an (ℓ, A) design if every ℓ-dimensional W meets them with total dimension at most A. For a folded code, H_i is the kernel of the map from messages to the i-th symbol. The design property says that no low-dimensional space of messages can hide from many symbols.
- *First engine: local equivalence.* The Levi-Mosheiff-Shagrithaya potential is Φ(V, U, R) = -n dim U + Σ_i dim(V_i ∩ U) + Rn dim U. Its threshold over coordinate subspaces is the rate threshold of any local property for random linear codes. Subspace-design codes reach the same threshold through the same intersection inequalities, which is the local equivalence of random and explicit codes.
- *Second engine: a Schubert-calculus lower bound.* Over algebraically closed fields, each condition dim(H_i ∩ W) >= 1 cuts a Schubert variety of codimension s - d + 1 in Gr(d, k). Hence some W can always be made to meet floor(d(k-d)/(s-d+1)) of the H_i.
- *Matroids.* The paper ties maximal recoverability of tensor codes to the bipartite-rigidity and matrix-completion matroids. Both are questions of identifying vertex unknowns from edge measurements.

**Hamiltonian sparsification (2605.02211).**
- A set of PSD terms is *non-redundant* if each term has a private witness: a state annihilated by every other term but not by it.
- The largest non-redundant set is a lower bound on the size of any sparsifier.
- The connectivity parameter N(α) gives the upper bound. It is the size above which every set of terms contains one that is α-dominated by the rest.
- The interaction of kernels is studied through automorphisms of joint ground states.

**Quantum-cut sparsifiers (2606.09728).**
- The quantum-cut Hamiltonian restricted to Hamming weight ℓ is the Laplacian of the level-ℓ Kikuchi (token) graph.
- The harmonic decomposition of functions on a slice (the Johnson-scheme isotypic spaces) is invariant for every graph, and leverage scores are bounded on each invariant block.
- The Alon-Kozma operator inequality compares any expander with the complete graph. It is built on the octopus inequality of Caputo, Liggett and Richthammer behind Aldous' spectral-gap theorem.

**Covers (2609.36761).** Covers are dependencies among hyperedges, at three strengths:
- *even covers*, mod 2;
- *Abelian covers*, integer combinations that produce a hyperedge;
- *Catalan covers*, stacks reducible by adjacent cancellation, which is how non-Abelian groups behave.

Catalan and Abelian covers coincide in arity 3 (by Hurewicz) and separate in arity 4 (by nilpotent groups).

**Multiway cut (2603.28700).**
- The approximation ratio of a mixture of rounding schemes is linear in the mixture weights.
- Hundreds of schemes, each drawn with its own parameter law, are combined by a computationally found LP solution.
- The result is certified with interval arithmetic.

**The pretty good criterion and syndrome duality (2608.09347).**
- *Criterion.* Suppose posterior sampling decodes a classical-quantum channel with bit error below δ. Posterior sampling here is the pretty good measurement, which is also the Petz recovery map. Then every code of distance δ has rate at most the channel's capacity.
- *Syndrome duality.* For a linear code C, the Holevo information of the pure-state channel equals the entropy of a Bernoulli-q string modulo the dual code C^⊥, with q = 1/2 - sqrt(p(1-p)). A transversal Hadamard exchanges the code with its dual and the noise parameter with its complement.
- *Dual covering.* The dual covering argument (Friedman-Tillich, Navon-Samorodnitsky) reads the distance of C as the expansion of a Cayley graph on the syndrome space.

**No low-degree tests (2608.00265).** Suppose the dual code can be decoded from random errors. Then codeword phase states cannot be told apart from states far from every codeword with few copies: low-degree tests do not exist.

## 2. The collective algebra as a projection

**Setting.**
- e = y - E y is a layer's post-activation fluctuation.
- U holds the top-K eigenvectors of Cov(y), the collective subspace of note XXXVII, and t = U^T e are the collective coordinates.
- ℬ is the algebra of functions of t, and E_ℬ = E[. | t] is the conditional expectation onto it.

**Operator-algebra reading.**
- E_ℬ is a projection of norm one onto a subalgebra (Tomiyama).
- The gate D = diag(1[z > 0]) is a projection-valued random variable.
- E_ℬ D is the conditional gate probability, a function on ℬ (diag(Φ((μ + U_z t) / s)) when the pre-activation residual is Gaussian with variance s²).
- The gates do not belong to ℬ. Everything below is about what that leaves.

**Definition (ℬ-independence).** The units are independent with amalgamation over ℬ when the residuals ε_a = e_a - E[e_a | t] are conditionally independent given t. This is the commutative case of operator-valued independence. The conditional cumulants are then the ℬ-valued cumulants, and Brillinger's law of total cumulance is the operator-valued moment-cumulant formula.

**Proposition 1 (evaluation identities).** Assume ℬ-independence with conditional means E[e_a | t] = u_a . t, where u_a is row a of U, and let v_a(t) = Var(ε_a | t). Then, for distinct units a, b, c, d:

    kappa(e_a, e_b, e_c, e_d) = T4[u_a, u_b, u_c, u_d],                         T4  = kappa_4(t),
    kappa(e_a, e_a, e_c, e_d) = u_c^T N_a u_d,                                     N_a = kappa(e_a, e_a, t, t),   K x K,
    kappa(e_a, e_a, e_a, e_b) = g_a . u_b,                                          g_a = kappa(e_a, e_a, e_a, t), K,
    kappa(e_a, e_a, e_b, e_b) = u_b^T N_a u_b + u_a^T N_b u_a - T4[u_a^2, u_b^2] + Cov(v_a(t), v_b(t)),
    kappa_4(e_a)              = 4 g_a . u_a - 6 u_a^T N_a u_a + 3 T4[u_a^4] + 3 Var v_a(t) + E kappa_4(eps_a | t).

*Proof.* Write e_c = u_c . t + ε_c and expand multilinearly. Take any cumulant that contains ε_c (c ≠ a) together with e_a or another unit. Given t, ε_c is independent of the rest and has conditional mean zero, so by total cumulance the cumulant vanishes. The (2,2) and diagonal formulas follow from the law of total cumulance with the conditional variance, skewness and fourth cumulant of ε. ∎

*Check.* `code/check_bindep.py` enumerates a 5-unit model with:
- a 4-atom collective law;
- skewed, heteroscedastic 3-point residuals.

All five identities hold to relative error 1e-15 to 1.5e-14 (`outputs/check_bindep.txt`).

**Reading.** Each unit carries *symbols on the collective space*:
- a vector g_a, whose evaluations at the other units' loadings give its (3,1) row;
- a quadratic form N_a, whose evaluations at pairs of loadings give its (2,1,1) block.

The (1,1,1,1) class is one global quartic T4 evaluated at quadruples of loadings. Values and derivatives of polynomial symbols at evaluation points are what a multiplicity code records (Kopparty-Saraf-Yekhanin; Guruswami-Kopparty). Here the evaluation points are the units' loadings u_a.

## 3. The flex theorem

The chain carries the *pair state*: the fourth-cumulant diagonal, the (2,2) slice and the (3,1) slice. Every closure in the programme reads it. Under Proposition 1 the pair state is a linear image of the symbols (T4, 𝐍, 𝐆, V, X), where:
- 𝐍 = (N_a) and 𝐆 = (g_a);
- V_ab = Cov(v_a, v_b);
- X_a = E κ4(ε_a | t).

**Theorem 2.** Fix the loadings and V. Suppose the rank-one forms u_a u_a^T span Sym²(R^K), the loadings u_b (b ≠ a) span R^K for every a, and n(n+1)/2 + q(q-1)/2 > n(q+1), with q = K(K+1)/2. For n = 1024 the last condition holds comfortably at K <= 32 and fails beyond K = 43. Then the kernel of the map (T4, 𝐍, 𝐆, X) -> (diagonal, (2,2), (3,1)) is exactly

    F = { (dT4,  dN_a = dT4[u_a, u_a, ., .] / 2 + mat(A vec(u_a u_a^T)),  dG = 0,  dX = 0) :
          dT4 in Sym^4(R^K),  A antisymmetric on Sym^2(R^K) },
    dim F = C(K+3, 4) + q(q-1)/2,

and the omitted classes move along it:

    d(2,1,1)_(a; c,d) = dT4[u_a, u_a, u_c, u_d] / 2 + vec(u_c u_d^T)^T A vec(u_a u_a^T),      d(1,1,1,1)_(abcd) = dT4[u_a, u_b, u_c, u_d].

Here vec is a Frobenius-isometric coordinate map on Sym²(R^K).

*Proof.*
1. *The (3,1) slice.* dg_a . u_b = 0 for every b ≠ a forces dg_a = 0.
2. *The (2,2) slice.* Let 𝐐 be the n x q matrix with rows vec(u_a u_a^T). The (2,2) slice is the off-diagonal part of sym(𝐍 𝐐^T), minus the T4 term, plus V. Its diagonal part enters only the fourth-cumulant diagonal, where dX absorbs it. A kernel element therefore satisfies one condition: sym(d𝐍 𝐐^T) minus the T4 term is diagonal.
3. *Solving that condition.*
   - dN_a = dT4[u_a, u_a] / 2 is a particular solution.
   - The solutions with sym(d𝐍 𝐐^T) = 0 are d𝐍 = 𝐐A with A antisymmetric. To see this: d𝐍 𝐐^T is then antisymmetric, so its column space lies in col(𝐐), and 𝐐 has full column rank.
   - Under the dimension condition, no nonzero d𝐍 makes sym(d𝐍 𝐐^T) diagonal without making it zero. For generic loadings this is the statement that an n-dimensional space of diagonals meets the image of d𝐍 -> sym(d𝐍 𝐐^T) only in zero.
4. *The diagonal.* -6 u_a^T dN_a u_a + 3 dT4[u_a^4] + dX_a = -3 dT4[u_a^4] + 3 dT4[u_a^4] - 6 vec(u_a u_a^T)^T A vec(u_a u_a^T) + dX_a = dX_a. So dX = 0. ∎

If V is not held fixed, so that the covariance of conditional variances is a free symmetric matrix, any d𝐍 can be absorbed into V. The (2,1,1) class is then not identified at all.

*Checks.*
- `code/check_bindep.py`: along a random element of each family, the pair slices change by at most 2e-14 and the omitted classes by O(1).
- `code/check_flex_rank.py`: builds the map on a basis and computes its nullity. It equals C(K+3, 4) + q(q-1)/2 exactly once the dimension condition holds:

| K | n | nullity | predicted |
|---|---|---|---|
| 2 | 10 | 8 | 8 |
| 3 | 14 | 30 | 30 |
| 4 | 24 | 80 | 80 |

  Below the condition the kernel is larger, as the proof says (K = 3, n = 9: 33; K = 4, n = 12: 89). See `outputs/check_flex_rank.txt`.

**Sizes of F.**

| K | free quartic, C(K+3, 4) | antisymmetric rotations, q(q-1)/2 |
|---|---|---|
| 16 | 3,876 | 9,180 |
| 32 | 52,360 | 139,128 |

**Three readings.**
- *Coding.* The pair state is a syndrome of the omitted classes, and F is the code whose cosets the syndrome cannot resolve.
- *Subspace designs.* The kernels H_a of the per-unit maps (symbols -> unit a's pair data) share a common intersection of dimension dim F. The slice code is therefore far from a design: a whole subspace of messages hides from every symbol. The Schubert bound of 2510.13777 says that some low-dimensional W always meets many symbol kernels. Here the situation is worse and exact.
- *Rigidity.* This is the reading of the paper's matroid section. Edge (a, b) of the (2,2) slice measures the vertex unknowns N_a, N_b, each at the other vertex's position. The antisymmetric family is a flex that exists for every configuration, the analogue of a trivial motion.

## 4. What the flex explains

**Every pair-slice closure picks a point of F by convention.**
- The lam core is one fitted scalar times C_off. The scale-mode core D of note XXXVI is the radial tangent, with an amplitude fitted per layer or taken from the trace sector. The G D core is KD.
- All three compute the omitted classes from pair-supported statistics.
- By Theorem 2, which already holds inside the most favourable model class, that information cannot determine those classes. Whatever the convention, its error is the truth's flex coordinate minus the convention's.
- The truth's flex coordinate is fixed by the history of the law, not by its present pair state.

Theorem 2 does not say in which direction a convention errs. It does say why note XXXVI's two findings can hold together:
- the (3,1) slice of z' is 80-97% the transport of classes that the pair state omits;
- a closure that captures most of that transport in L2 can still leave the part the readout uses untouched. D explains 62-81% of it, yet KD lowers the slice L2 error by 30% and raises the output error by 9%.

**The compensation is one coordinate seen three times.**
- The transport R of the omitted classes enters the next layer's diagonal, (2,2) and (3,1) slices together (note XXXVI section 3). A missing flex coordinate therefore makes the errors of g4, K22 and K31 correlated.
- Replacing one of them by truth, while the other two keep the flex misfit, breaks that correlation. This is the oracle finding that one true slice alone makes the output worse (true K22 alone: +13% to +25%). It is also the bundle requirement of note XXXVII section 6.2.
- In the language of 2605.02211, the three slices are redundant terms of the readout's energy, because their errors share a witness. A slice is worth correcting only together with the others that see the same coordinate.

**Two kinds of invisibility.**
- *Invisible to the readout.* The quotient programme of note XXXVI (sections 3j-3l) found the null code: state perturbations (G/12, μ.G/4, Σ.G) that no future mean can see, by homogeneity. This invisibility can be quotiented out at no cost.
- *Invisible to the present state.* That is the flex. The downstream layers see it, so it has to be carried.

Syndrome duality frames the pair:
- the null code is the kernel of the readout, and its dual is the adjoint, the suffix kernel of note ray-compiler;
- the flex is the code modulo which the pair state determines the omitted classes.

Note XXXVI concluded that "gauge covariance is not where the remaining error is". Theorem 2 says where an irreducible part of it is.

## 5. The minimal extension

**What would be carried.** Under ℬ-independence the omitted classes become explicit if the state carries:
- T4, C(K+3, 4) numbers;
- the per-unit symbols N_a^off, n q numbers.

Here N_a^off is the collective projection U^T Ψ_a U of unit a's (2,1,1) block Ψ_a, and T4^off is the collective projection of the all-distinct class. Both are exact definitions; they do not assume ℬ-independence, which is what predicts that the projections reproduce the classes.

**The model tensor and its transport.**
- Set Ã = W U and M_a = N_a^off - T4^off[u_a, u_a, ., .].
- The model tensor is Ŷ = T4^off[U^4] + Σ over the 6 slot pairs of δ_pair (U M_a U^T). It equals the model on the omitted classes.
- Its transport W#[Ŷ] has the slices

      diag_i   = T4[Ã_i^4] + 6 Σ_a W_ia^2 Ã_i^T M_a Ã_i,
      (2,2)_ij = T4[Ã_i^2, Ã_j^2] + Σ_a [ W_ia^2 Ã_j^T M_a Ã_j + W_ja^2 Ã_i^T M_a Ã_i + 4 W_ia W_ja Ã_i^T M_a Ã_j ],
      (3,1)_ij = T4[Ã_i^3, Ã_j] + 3 Σ_a [ W_ia^2 Ã_i^T M_a Ã_j + W_ia W_ja Ã_i^T M_a Ã_i ].

- The transport of Ŷ's own pair slices (note XXXVI's T_pair program) is subtracted from these.
- In a chain the subtraction folds into the transport the chain already does: chain' = T_pair(y slices - Ŷ_pair) + W#[Ŷ].

**Cost.**

| step | cost |
|---|---|
| Q_ia = Ã_i^T M_a Ã_i and Σ_a W_ia^2 M_a | O(n^2 K^2) |
| (W o Q) W^T | one n^3 product |
| T4[Ã_i, Ã_i, ., .] | O(n K^4) |

That is about 1.5 units per layer at K = 16 and about 4 at K = 32, before Strassen. The cross term 4 Σ_a W_ia W_ja Ã_i^T M_a Ã_j of the (2,2) slice costs n^3 K over the full slice. Its sum over a is incoherent (random signs). Test S measures its share before a chain would pay for it.

**What the theory does not settle.**
- *Representability.* Does the truth satisfy ℬ-independence well enough that its omitted classes are their own collective projections at an affordable K? Test S measures this.
- *Transport of the symbols.* Two maps generate the symbols at layer l + 1 from those at layer l:
  - the gate, which acts per unit on the conditional moments on ℬ;
  - W, which is exact and linear from y to z'. Since t' = U'^T W U t + U'^T W ε, κ4(t') is A^(⊗4) T4 plus mixed terms that are again symbols.

  The law is only outlined here. Deriving it is the next obligation if test S passes.

  One piece is closed. At a Gaussian collective state the gate's symbol is rank one along the unit's own loading:

      N_a = c_a (Λ u_a)(Λ u_a)^T,      c_a = E[(f^2)''] - 2 E[f']^2 - 2 E[f] E[f''] + E[g''],

  with f and g the conditional ReLU mean and variance as functions of p = u_a . t. It follows from Stein's identity applied twice and is checked by quadrature to 2e-4, the finite-difference error (`code/check_symbol_gate.py`).
  - Its (2,1,1) evaluation is c_a C_ac C_ad on collective covariances: the one-loop (tree) content.
  - Everything beyond it comes from T3, T4 and the transported heteroscedasticity.

## 6. Test S (pre-registered; written before any run)

**Question.** Is the transport R of the truth's omitted classes equal to the transport of their own collective projections, R_pred(K), with nothing fitted?

**Data.** Network 1.
- *R, per transition l -> l + 1.* Computed from mc4: 1.6e7 inputs, R = t - T_pair(y). R is exactly the transport of the y (2,1,1) and (1,1,1,1) classes (note XXXVI section 3).
- *The symbols.* Computed from a new, independent Monte Carlo (`code/mcsym.py`): 8e6 inputs in 16 chunks, two halves.
  - Taken at the y-level of layers 2, 5, 8, 11 and 14.
  - U is the top-64 eigenvectors of the true Cov(y), and t = U^T e.
  - The run accumulates the raw moments from which N_a^off and T4^off are formed exactly, with every index coincidence removed by inclusion-exclusion (`code/symtest.py`).
  - K = 8, 16 and 32 use the leading coordinates of the same t.

**What is computed.**
- R_pred(K) and its split into the (2,1,1) and (1,1,1,1) parts.
- For each slice:
  - the explained share 1 - |R - R_pred|^2 / |R|^2, at coefficient 1 and noise-corrected with the halves of both Monte Carlo runs;
  - the best-fit coefficient and the correlation.
- The same for candidate D at its fitted coefficient, and for D and R_pred fitted jointly.
- The share of the incoherent W_ia W_ja terms.
- A width-7 network check against the brute-force tensor comes first: at K = n, R_pred must equal R exactly on the same samples.

**Predictions.**
- **P1 (representability).** At K = 32, R_pred explains at least 0.5 of R's energy in the (3,1) and diagonal slices at transitions 8 -> 9, 11 -> 12 and 14 -> 15. At 2 -> 3 and 5 -> 6 it explains at least 0.3. The reasons:
  - R is a W²-weighted sum over units, which averages the units' non-collective parts down by about sqrt(n);
  - the y-level law is collective at depth (note XXXVII: the collective part carries 0.91-0.94 of the third-cumulant slices at layer 15).
- **P2 (convergence).** The explained share rises with K at every transition (8 < 16 < 32 < 64). From 8 -> 9 on, Δ(16 -> 32) <= Δ(8 -> 16).
- **P3 (beyond the scale mode).** From 8 -> 9 on, fitting R_pred(32) jointly with D raises the explained share of the (3,1) slice by at least 0.1 over D alone: the symbols carry what D misses.
- **P4 (failure criterion).** Suppose that at K = 64 R_pred explains less than 0.3 of the (3,1) slice from 8 -> 9 on. Then ℬ-independence is rejected for the omitted classes: their content lies in conditional dependence beyond the collective algebra, and a carried-symbol state is not the lever.

## 7. The other links, placed

- **Kikuchi graphs and invariant blocks (2606.09728).**
  - The level-ℓ Kikuchi graph is the ℓ-particle sector of a pair Hamiltonian.
  - The counterpart here: the ℓ-th cumulant lives on Sym^ℓ of the unit space; the chain's slices are classes of index coincidence; under ℬ-independence Sym^ℓ(R^K), the collective "token space", carries the all-distinct class.
  - Block-wise matrix concentration is a sampling tool, and sampling is closed at the precision needed (note ray-compiler, the Funk-Hecke ceiling).
  - Relabelling.
- **The octopus inequality and Aldous' theorem.**
  - Their counterpart is that the slowest k-leg mode is a power of the slowest one-leg mode.
  - That is the concentration of accumulated content in the top covariance eigenspace (note XXXVII section 5).
  - Relabelling.
- **Covers (2609.36761).**
  - Around the independent reference, the cumulant expansion is a sum over connected hypergraphs of tangent entries, with vertex weights Π_v c_v(1, deg v). Its first diagrams are the legs, the facet term, GC1, GC2, the triangle and the double edge of note XXXVI.
  - At a zero-mean unit c(1, k) = 0 for odd k >= 3, so nonlinear vertices need even degree. That even-cover condition is the deck parity of note XXXVII section 1.
  - Catalan against Abelian, local against commutative cancellation, is the planar against non-planar split of the annealed weight expansion. The chain is quenched and sums all pairings.
  - The one use: the diagram grading is the completeness rule of note XXXVI section 3h (a repair must include every diagram of its order). Theorem 2 sharpens it for the fourth-order sector: the missing diagrams are the ones that carry the flex coordinate.
- **Non-redundancy (2605.02211).**
  - It is used in section 4.
  - Sparsification itself is of no use here. Keeping a weighted subset of terms that preserves every energy to 1 ± ε needs about n / ε² terms, which is 1e7 at ε = 1e-2.
- **Multiway-cut mixtures (2603.28700).**
  - The output-metric calibration of note ray-compiler section 8 is a mixture of closure directions with weights fitted on training networks (validated at -3.76%).
  - The paper's lesson is that a large mixture with a rigorous factor-revealing analysis can beat a single elegant scheme. Our factor-revealing analysis is the first-entry ledger.
  - No new number.
- **The pretty good criterion (2608.09347).**
  - The pretty good measurement is the Petz recovery map of note conditional-modular.
  - Posterior sampling as a converse tool has no counterpart here: we estimate a quantity, we do not bound a rate.
- **Syndrome duality, Gaussian form.**
  - It is used in section 4.
  - The exchange of p and 1/2 - sqrt(p(1-p)) is the rotation that exchanges a correlation ρ with sqrt(1 - ρ²).
  - The arc-cosine kernel carries both: relu = z/2 + |z|/2, and E|u||v| is proportional to sqrt(1-ρ²) + ρ arcsin ρ.
  - It changes no computation here.
- **No low-degree tests (2608.00265).**
  - When the dual code is decodable from noise, low-degree statistics cannot distinguish codeword states.
  - Theorem 2 is the same kind of statement: the pair state, a low-order statistic, cannot distinguish laws that differ along F.
- **Local equivalence and potential functions (2510.13777).**
  - The threshold of a local property is decided by the densest coordinate subspace. That is large-width power counting, already used throughout the programme.
  - Relabelling.
- **Covering radius and the loading cloud.**
  - Theorem 2's hypotheses and the conditioning of the (2,2) system ask the loading cloud {u_a} for a frame property: the rank-one forms u_a u_a^T must span Sym²(R^K) robustly, an approximate spherical 2-design.
  - The covering radius of the normalised cloud bounds that conditioning from one side. With n = 1024 points in R^32 it is not the binding quantity; the frame potential is.
  - No number changed.
- **The gate-geometry synthesis (supplied PDF, previous round).**
  - Its arithmetic is correct where it can be checked. Its central claim, that a full-rank gated kernel can be applied at linear cost through interval geometry, is true (its Theorem 9.1 and Proposition 9.2).
  - Every route it builds on that claim is a conditionally integrated sampler: source planes, finite rotation groups, circle tiles. They meet the precision ceiling of note ray-compiler.
  - Its premise that cheap transport need not be low rank does not help the fourth-order flex. The data and Theorem 2 place the readout-relevant hidden content in the collective (low-rank) space and in the bulk, and section 8 measures how slowly the collective part converges.
  - The synthesis says itself that no leaderboard gain follows.

## 8. Results (network 1; `outputs/symtest_off1.txt`, `outputs/symread_off1.txt`)

**Data.**
- Symbols: 1.6e7 new inputs in 16 chunks (`code/mcsym.py`, seed `[4713, 1, chunk]`), centred on the mc4 means. Every mean correction is exact, and the bookkeeping was checked on a width-7 network (`outputs/check_symtest_tiny.txt`):
  - the R identity holds to 1e-12;
  - the symbols match the brute-force projections to 1e-13 at K = 4 and 7;
  - R_pred = R at K = n to 1e-14.
- Noise: R's noise share is 0.18-0.24 at 2 -> 3, 0.03-0.06 at 5 -> 6 and below 0.02 deeper. R_pred's is at most 0.005.

**Explained share of R** (noise-corrected, R_pred at coefficient 1; D at its fitted coefficient; slices diagonal | (2,2) | (3,1)):

| transition | D (fitted coef) | K = 8 | K = 16 | K = 32 | K = 64 |
|---|---|---|---|---|---|
| 2 -> 3 | .385 .222 .276 (2.7) | .072 .059 .074 | .152 .104 .120 | .245 .151 .193 | .366 .275 .309 |
| 5 -> 6 | .553 .371 .456 (3.9) | .267 .166 .214 | .347 .218 .285 | .460 .305 .380 | .591 .423 .499 |
| 8 -> 9 | .615 .458 .523 (4.4) | .375 .282 .321 | .454 .334 .385 | .551 .415 .478 | .668 .527 .592 |
| 11 -> 12 | .761 .598 .638 (4.6) | .643 .460 .487 | .704 .517 .549 | .761 .583 .625 | .816 .671 .709 |
| 14 -> 15 | .817 .680 .721 (4.7) | .704 .567 .612 | .743 .605 .655 | .816 .675 .719 | .874 .758 .791 |

- R_pred's best-fit coefficient is 0.92-1.13 everywhere, and 0.99-1.08 at K >= 16 from 8 -> 9 on.
- The (1,1,1,1) part's norm is 0.09-0.41 of the (2,1,1) part's at K = 32 and 0.17-0.50 at K = 64, rising with depth.
- The incoherent W_ia W_ja term is 5-7% of the predicted (2,2) slice.

**Verdicts on the registered predictions.**
- **P1 fails.**
  - At K = 32 the (3,1) share is 0.478 at 8 -> 9, against the registered 0.5. At 2 -> 3 it is 0.193, and the diagonal 0.245, against the registered 0.3.
  - The other registered thresholds hold: 5 -> 6 at .460 / .380; 11 -> 12 at .761 / .625; 14 -> 15 at .816 / .719.
- **P2 fails in its second half.**
  - The share rises with K at every transition and in every slice.
  - It does not converge: Δ(16 -> 32) exceeds Δ(8 -> 16) at 8 -> 9 (diagonal .097 against .079; (3,1) .093 against .064), at 11 -> 12 ((3,1) .076 against .062) and at 14 -> 15 (diagonal .073 against .039; (3,1) .064 against .043).
  - Each doubling of K adds 0.05-0.12, so the content is spread over many collective modes.
- **P3 fails.**
  - Fitted jointly with D at K = 32, R_pred adds .009 / .017 / .018 to D alone in the (3,1) slice from 8 -> 9 on, against the registered 0.1.
  - At K = 64 it adds .073 / .074 / .072.
- **P4 is not triggered.** At K = 64 the (3,1) share is .592 / .709 / .791 from 8 -> 9 on. ℬ-independence is not rejected for the collective part of the omitted classes.

**Readout-weighted shares (post hoc, `code/symread.py`).** Each slice is weighted by the Hermite coefficients through which the next layer reads it: the mean map's c4 on the diagonal, c(1,2)c(1,2)/4 and c(1,3)Φ/6 in the covariance program. Coefficients are refitted in this metric.

| transition | D: diagonal \| (3,1) (coef) | R_pred K = 32 | R_pred K = 64 | K = 64 jointly with D |
|---|---|---|---|---|
| 8 -> 9 | .564 \| .421 (4.5 \| 4.9) | .516 \| .370 | .635 \| .492 | .643 \| .510 |
| 11 -> 12 | .696 \| .397 (4.7 \| 5.3) | .680 \| .390 | .745 \| .513 | .755 \| .524 |
| 14 -> 15 | .807 \| .391 (5.2 \| 5.7) | .814 \| .408 | .871 \| .546 | .874 \| .555 |

R_pred's coefficient in this metric is 0.96-1.11 at K >= 32 (0.96-1.09 from 8 -> 9 on).

**What the numbers say** (derived from them, nothing fitted to the output).

1. **The collective symbols carry the omitted-class transport with the right amplitude.**
   - Their coefficient is 1 within 9% from 8 -> 9 on, in L2 and in the readout metric, and within 13% everywhere.
   - D needs an amplitude fitted per layer, 2.7 -> 4.7 in Var(r²) units, which note XXXVI could not derive at depth: the trace-sector estimate g_conn falls while the fitted value rises.
   - D's best amplitude also shifts by 10-20% between the two metrics (14 -> 15 (3,1): 4.71 in L2, 5.66 in the readout metric), the mark of a shape that errs in the direction the next layer reads.
   - The symbols measure the generated-scale amplitude that the scale-mode core lacked. Proposition 1's structure holds for the collective part of the truth.
2. **At affordable resolution the hidden coordinate is the scale mode.**
   - At K <= 32 the symbols span what D spans. Fitted jointly they add at most 0.02 to D in L2 and at most 0.05 in the readout-weighted (3,1) slice.
   - The flex coordinate that a K <= 32 collective state can hold is therefore, for the readout, the generated-scale amplitude. Carried at its fitted value, that amplitude made the chain worse (KD, +9%).
3. **Beyond K = 32 the symbols carry what D misses, slowly.**
   - At K = 64 the readout-weighted (3,1) share is .49-.55 against D's .39-.42 at depth.
   - Half of the readout-relevant (3,1) content is still missing at K = 64.
   - Shallow transitions are mostly non-collective: 0.31 at K = 64 at 2 -> 3.

## 9. What follows

**Kept.**
- *The flex theorem.* The pair state determines the omitted classes only modulo F. So every closure the programme has tried for them, and any pair-slice closure, is a convention, and the hidden coordinate has to be carried to be right.
- *The two invisibilities.*
  - The null code is invisible to the readout and is quotiented at no cost.
  - The flex is invisible to the state and is carried at a cost.
- *The structural measurement.* The omitted-class transport is collective with the right amplitude, and at affordable K it is the generated-scale mode.

**Closed for cost.** This is the same verdict as note XXVIII section 3. The carried-symbol state:
- at K <= 32, equals the scale-mode core that the chain already found output-adverse;
- at K = 64, would cost about 16 units per layer for the T4 contraction (n K^4) plus about 4 for the symbols, roughly 0.3 B over the network. The gain it could buy is bounded by the K31-sector oracles: true slice -6% / -22%; column projection on the top 64 modes -3.2% / -13.9% (note XXXVII section 6.2).

The symbol transport law is therefore not derived further.

**Open.**
- Whether a basis other than the top variance eigenvectors of Cov(y) concentrates the omitted classes faster. A basis computable from the carried state would change the cost argument.
- Whether network 0 repeats these numbers. Earlier structural measurements agreed between the two networks to the second digit (notes XXXVI, XXXVII). The decision above does not depend on it.
