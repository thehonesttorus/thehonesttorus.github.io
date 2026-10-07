# The hidden separator, tested

Working note XX. Two syntheses were pasted into the session. The first argues that a projection can destroy visible
Markovianity without destroying a small representation of the underlying process (a hidden sign makes four visible
signs independent; eliminating it creates interactions of every order), connects that to recent Gibbs results on
conditional mutual information after coarse observation (Chen-Rouze, Yang, Bakshi-Liu-Moitra-Tang,
Zhang-Gopalakrishnan), gives a consistency condition for proposed local conditional rules with a Hodge repair, a
Dobrushin-type bound that propagates local rule errors to a readout, a decaying-memory contract for a progressive
generator, and the exact Ising law of the layer-1 signs given their magnitudes; it ends with a "conditional
specification compiler" whose research target is a small, weight-accessible boundary memory that makes the
completion rules local. The second reads two duality papers (Bauer-Cvetko-Vah on skew Boolean algebras with
intersections; Bezhanishvili-Bezhanishvili-Gabelaia-Kurz on bitopological duality for distributive lattices and
Heyting algebras) as the algebra of such a compiler: sections as (region, rule) pairs, rectangular bands as the
bracket, Heyting implication as refinement-stable certification, and minimization through dual observables.

The standard is the one of notes XVII-XIX: a claim counts when it changes a number or a decision. Two tests were run
for this note (`code/`); the rest is measured in the earlier notes and is cited by number.

What the tests say:

- **The hidden separator of the remaining error exists, and it has dimension about 64.** The gate-covariance term of
  note XVIII section 5, the one incoherent transport defect that is measured and sits at the n^4 wall, is a linear
  functional of the pre-activation correlation matrix. Replacing that matrix by its rank-K truncation captures the
  term with slope 0.22 at K = 1, 0.40 at K = 4, 0.74 at K = 16 and 0.96 at K = 64, tracking the Frobenius mass of
  the truncation (0.19, 0.35, 0.67, 0.91). So the omitted memory is the common-mode subspace of the pre-activations,
  the effective rank 60-70 that the Gaussian part of the norm fluctuation already implied. This revises note XVIII's
  "no low-rank route": there is one, at rank 64. It costs K transports per source and layer, and the score's linear
  charge for flops makes every K a loss (section 2).
- **The exact layer-1 Ising model is a zero-temperature dense spin glass.** Its couplings at typical magnitudes have
  median 120-280 and rms 230-670 on the two networks (the precision matrix of a square He matrix has condition
  number 7e5-3e6), and the Dobrushin row sums are 1000-1023 against the uniqueness threshold 1. Every Gibbs result
  the synthesis cites is a high-temperature, bounded-degree result; the network sits at the opposite corner on both
  counts, and the deeper layers saturate the gates along the mean direction, which is also low temperature.
- **The accessible small separators are already in the chain, and the error is orthogonal to them.** A small hidden
  variable corrects a mean through the curvature of the output in that variable, which is a coherent correction
  along a few fixed directions. The chain's error is white: 0.1% coherent, 0.3% along the output mean, 0.8% along
  the top 16 covariance eigenvectors, R^2 0.04 against every derived local feature (note XVIII section 4). The two
  separators that the fibre measurements found (the scale, and the mean-direction projection) are what the chain's
  mean-coupled third-cumulant channel carries.
- **The duality reading is the algebra of what the chain already does**, with one genuine content (the rectangular
  band is the bracket's identity set, which note XIV used) and one dead end for this problem (dual minimisation of
  a full-vector readout minimises nothing; the version that works is the Lyapunov compression of the sources, done
  and measured).

## 1. The separator hypothesis against the error anatomy

The hidden-sign example is exact and the point is right: high-order visible correlations do not imply a high-order
generative description. The question for the chain is whether its remaining error is the shadow of a small
omitted variable. There is a measurement that answers it without naming the variable.

If the true law were a mixture over a hidden variable h of laws the chain computes exactly given h, the chain's
output would be m(h~) for a reference value and the truth would be the average of m(h), so the error would be
(1/2) m''(h~) Var h at leading order: a fixed vector (or, for a d-dimensional h, a vector in a d(d+1)/2-dimensional
span) times a scalar. A small separator makes a coherent error. The scale is the worked case: m(G) = G m by
homogeneity, the error is -(g/8) m along the output mean, and the chain carries it. The measured error of the best
chain (note XVIII section 4) has no such direction: its mean is 0.1% of the MSE, its component along the output mean
0.3%, along the top 16 eigenvectors of the last covariance 0.8%, and no function of the carried last-layer state
predicts it (R^2 0.04 jointly, 0.6% per feature). The fibre measurements of note XVIII's record (E4) found exactly
two accessible hidden scalars: conditioning on the norm removes the excess kurtosis of random directions entirely
(0.066 to -0.001), and conditioning on the mean-direction projection removes most of it (to 0.008); the output
law is, to fourth order in random directions, a two-parameter mixture of Gaussians. Both parameters are in the
chain: the mean-coupled third cumulant 1.5 g mu sigma^2 is their coupling, and the error is orthogonal to it by
construction at the fitted optimum of its amplitude (an amplitude at its optimum has zero projection of the
residual on its derivative direction).

So for every candidate separator that is accessible from the network's own structure the chain's error is
orthogonal. A separator whose curvature direction is uncorrelated with every structural direction is not excluded
by this, but it is not small in any sense that can be found or used: the paste's own condition, "weight-accessible",
is the one it fails.

## 2. The rank-K test on the term where the error lives

Note XVIII section 5 computed the first-order gate-covariance correction to the all-distinct third-cumulant
transport exactly on 48 output neurons of layer 11 from the live source legs of layer 10,

    D3corr_i = 3 sum_abc W_ia W_ib W_ic rho_ab phi_a phi_b Phi_c T_abc,

found it to be 1.09% of the product-gate value, incoherent, about a third of the error, and n^4 per source-layer to
carry, and wrote that there is no low-rank route because rho is the bulk correlation. The paste's hypothesis is
precisely that such a term is the shadow of an omitted separator. The test (`code/gatecov_rank.py`,
`outputs/gatecov_rank_off0_l10.txt`) replaces rho by its rank-K truncation (top-K eigen-directions of the
off-diagonal correlation, diagonal removed) and recomputes the same contraction on the same neurons:

| K | Frobenius mass of rho captured | correlation with the exact term | slope | residual rms / exact |
|---|---|---|---|---|
| 1 | 0.186 | 0.44 | 0.22 | 90% |
| 4 | 0.353 | 0.64 | 0.40 | 77% |
| 16 | 0.671 | 0.87 | 0.74 | 50% |
| 64 | 0.912 | 0.986 | 0.96 | 17% |
| 256 | 0.946 | 0.994 | 0.96 | 11% |

The off-diagonal correlation at layer 10 has eigenvalues 43, 26, 23, 22, ... (16th 13, 64th 3.2, 256th 1.0), and the
correction follows its spectrum with no concentration beyond it: the term is the shadow of a separator of dimension
about 64, the common-mode subspace of the pre-activations, which is the effective rank n_eff = 60-70 that
2 tr S^2 / tr(S)^2 = 0.028 gives at depth. The paste is right that the omitted memory is small compared with 1024.
Note XVIII's "no low-rank route" is withdrawn: the route is rank 64, and the sentence that stands is that it is not
cheap.

**Its cost under the score.** For one direction u with eigenvalue lambda, the correction is 3 lambda T(x_i, x_i,
y_i) with x_i = W_i o phi o u and y_i = W_i o Phi, and for all output neurons at once it is one n^2 x rank product
per source, the same cost as the product-gate transport the chain already does. So K directions cost K transports
per source-layer. The score is MSE x F / B with the chain at F / B = 0.26, so a change is worth it only when the
relative MSE gain exceeds the relative flop increase. One transport over the depth is about a fifth of F. The term
is about a third of the MSE in quadrature (note XVIII), so capturing a fraction s of its rms removes about
(1 - (1 - s)^2) / 3 of the MSE: K = 1 removes 6% of the MSE for about 20% more flops, K = 4 removes 12% for 80%, K =
16 removes 22% for 320%. No K pays. The chain's own compression knobs were tuned by search to the same marginal
rate, which is why they stopped where they did (note XVII: 10% raw per 64 ranks of the shared basis at the edge).

**What the separator is, and why it is medium.** The common-mode subspace is the image of the mean direction and
the few expanding directions of the depth dynamics in the pre-activation correlation; its dimension is set by
the He-critical spectrum (one expanding direction, a marginal bulk), not by a small number of hidden symbols. A
genuinely small separator would need the correlation to be rank 1-4, and it is not: the top direction carries 19%
of its mass. This is the quantitative form of the warning the paste draws from Zhang-Gopalakrishnan: compressing the
hidden state (the chain's product gates assume independent gates given the carried moments) changes the visible
local dependencies, and the change is carried by a subspace of the size of the correlation's effective rank.

## 3. The exact Ising model is at zero temperature

The paste records that for an invertible first layer, Z = W_1 x is Gaussian with covariance Sigma = W_1 W_1^T, and
conditional on the magnitudes r the signs follow an Ising law with J_ij(r) = -r_i (Sigma^-1)_ij r_j, and it notes
that the couplings depend on the precision matrix. The numbers (`code/ising_l1.py`, `outputs/ising_l1.txt`):

| | network 0 | network 1 |
|---|---|---|
| eigenvalues of Sigma, min / max | 1.1e-5 / 8.06 | 3.0e-6 / 7.92 |
| condition number | 7.4e5 | 2.7e6 |
| couplings at typical magnitudes, median / rms / 99th percentile | 120 / 232 / 744 | 281 / 670 / 2345 |
| Dobrushin row sums sum_j tanh |J_ij| , median (max) | 1020 (1023) | 1022 (1023) |
| the same at one random input, median | 1000 | 1008 |

A square He-Gaussian matrix has a Marchenko-Pastur spectrum with edge at zero, so the precision matrix is dominated
by its smallest eigenvalues and the couplings are of order 10^2-10^3: every sign is coupled to essentially every
other sign at a strength that makes the conditional law of each sign a step function of the rest. The Dobrushin
condition needs a row sum below 1; the measured row sums are 1000. This is a zero-temperature, dense, mixed-sign
spin glass, and nothing about it is sampled or marginalised cheaply; what the network does with it (relu(Z) =
r o (1 + sigma) / 2) is the Gaussian computation the chain performs, in which the magnitudes and signs are not
separated. Every Gibbs result cited in the paste's section 3 assumes bounded interaction degree and high
temperature, or finite range on a lattice; the network has degree 1023 and coupling 10^2 at layer 1, and at depth
the gates saturate along the mean direction, which is also low temperature in the pattern coordinates. Their
hypotheses are the negation of the problem, and the paste's own caution ("small covariance entries do not establish
a high-temperature Dobrushin condition") is confirmed by three orders of magnitude.

## 4. The consistency condition, the influence bound and the allocation rule

Three pieces of the paste are about an estimator that has local conditional rules and a sampler. The chain has
neither: it is a moment closure, its "rules" are the integration-by-parts programs, and it never samples. So the
square condition for conditional odds and its Hodge repair have nothing to check, and the Dobrushin comparison
|nu f - mu f| <= delta(f)^T (I - C)^-1 b bounds a stationary law the chain does not construct.

The allocation principle behind the bound is nonetheless the principle the chain was tuned by. Its backward
importance weight w = (I - C)^-T delta(f) is, in the chain, the response of the final mean to content born at
layer b, measured directly in note XVII (rho^age with rho = 0.80-0.92 on the third-cumulant channel, 1 on the mean
direction), and the error tail of a rank-r compression of a source is algebraic in r (note XVII). With cost linear
in rank and an algebraic tail, the paste's allocation b_i = (p a_i / lambda w_i)^(1 / (p + 1)) gives a rank falling
with age as w^(1 / (p + 1)); for p = 1 and seven layers of age that is 0.80^3.5 = 0.46 to 0.92^3.5 = 0.75 of the
newborn rank. The ladder the search found is 224 / 384 = 0.58. That is consistency, not a derivation, and the
search already found the ladder; the rule would save the search, not the flops.

## 5. The decaying-memory contract

The bound D(P || Q) <= (1/2) sum_t epsilon_t^2 for a short-memory decoder with per-step log-likelihood errors
epsilon_t, and the conclusion that exponentially decaying memory needs logarithmic retained history, has one
hypothesis the network violates: note XIX section 1 measured the memory of the depth dynamics to be geometric
(0.80-0.92 per layer) on every mode but the scale, whose transport is exactly 1 by homogeneity. On the scale the
variation sequence is constant, the sum of epsilon_t^2 does not converge with retained history, and the decoder must
carry the mode as a state. The chain does (the gain, then the coherent sources). On the bulk the retained history
is all sixteen ages, which is already the logarithm in the trivial sense; the contract's content for us is the
bookkeeping of per-layer conditional-law errors, which the chain has in moment form: the 1.1% gate-covariance
defect per layer accumulating in quadrature to the output error (note XVIII section 5). The bounded-readout caveat
the paste raises (Gaussian outputs are unbounded; integrate the radial part first) is the homogeneity reduction of
note XIX section 2.

## 6. The duality reading

**The rectangular band is the bracket.** (l, r) * (l', r') = (l, r') satisfies x * x = x and x * (y * z) = x * z =
(x * y) * z, which are the identities of the Smale bracket [x, [y, z]] = [x, z] = [[x, y], z]. This is correct, and it
is the content note XIV section 5 took from the Smale paper: the bracket joins compatible pasts and futures, and the
measured structure (how the joins are weighted) is a separate ingredient. The paste draws the same line ("skew
structure: what can be joined; measured structure: how it must be weighted"). The chain's state is the measured
structure, and the only place the joins failed to be a product was the gate covariance of section 2.

**Sections forget the rule; the chain's state is the rule.** A section is a region of addresses with a rule chosen
on it, and the Boolean quotient keeps the region and forgets the rule. In the chain the region is the gate pattern
and the rule is the carried Gaussian state with its cumulant corrections; note XVI's pullback transport is the
statement that the rule, not the region, is what is carried. The paste's warning that passing to the Boolean
quotient is not mean-preserving is the reason the chain was never a symbolic estimator (note XVIII section 1).

**Heyting certification is the stratification bound.** The refinement-stable certificate "U(q) - L(q) <= epsilon" for
a readout on the completions compatible with a partial address is an interval bound for a Lipschitz readout on a
1024-dimensional Gaussian input. The interval closes only when the address is essentially complete, which is the
exponent 1 + 1 / 512 of note XVIII section 1: cell quadrature that costs more than Monte Carlo. The two-bit example
in the paste (f = max(x - y, 0) decided by (0, ?) or (?, 1)) is exact and does not scale: the output of the network
is a mean over 2^1024 cells and no partial address of accessible size pins it.

**Dual minimisation of a full-vector readout minimises nothing.** The observable space V = span{T_w^T o_j} is the
smallest space containing the readouts and closed under the transposed transitions. The readouts here are all 1024
output coordinates, and the gated propagator W diag(Phi) is full rank, so V is the whole space at every layer and
the reduced state is the state. The version of the idea that applies is to compress the sources by what the
transport does to them, which is the Lyapunov compression of task 6 (note XV, note XVII): the shared basis of rank
384 and the ladder are its measured answer, and the residual rank 4 is its edge. "Minimisation via duality" names
that computation; it does not shrink it.

**The completion compiler.** The paste's architecture (partial rules on certified regions, skew operations for
restriction and priority, exact intersections for agreement, dual minimisation, then pseudorandomisation of what
remains) is a disciplined description of an adaptive estimator with routing. Its costs, which the paste lists
(discovering validity regions, testing predicates, routing, integrating the selected rule), are the costs note
XVIII section 6 priced for the refilling calculus: any cheap retained description leaves 99% of the output variance
unexplained (E9), the frame of 256 directions keeps 16%, antithetic 9%, radial 0.7%, and the budget is short by
573x of what conditional sampling needs. A compiler does not change those integrals.

## 7. What transfers

| from the pastes | what it is here | number or decision |
|---|---|---|
| a projection can hide a small separator | the gate-covariance term is the shadow of the common-mode subspace | rank 64 captures 96% of the term; no K pays under the score; XVIII's "no low-rank route" withdrawn as stated, kept as priced |
| the exact Ising law of the layer-1 signs | a zero-temperature dense spin glass | couplings 10^2-10^3, Dobrushin row sums 1000 |
| high-temperature CMI decay, local recovery | hypotheses are the negation of the problem | degree 1023, coupling 10^2 at layer 1; saturated gates at depth |
| consistency of local odds, Hodge repair | no local rules to check | - |
| influence bound and allocation | the measured response rho^age and the rank ladder | 0.46-0.75 predicted ratio for age 7, 0.58 found by search |
| decaying-memory contract | memory is non-summable on the scale; the chain carries it | row sums 1 (note XIX) |
| rectangular band = bracket | the identity set of the Smale bracket; the measured structure is separate | as note XIV section 5 |
| sections and the Boolean quotient | gate region and carried state; the quotient is not mean-preserving | as note XVIII section 1 |
| Heyting certification | the stratification bound | exponent 1 + 1 / 512 |
| dual minimisation | the Lyapunov compression of the sources, already done | ranks 384 / 224 / 4 |
| completion compiler | the refilling calculus with routing costs | E9: 99% unexplained; 573x short |

The one thing the pastes change is a sentence: the remaining error is not at an n^4 wall with no low-rank route;
it is the shadow of a 64-dimensional separator whose retention costs 64 transports per source-layer, and the
score's price of a flop is what forbids it. That is a sharper statement of the same wall, and it names the design
problem of task 2 precisely: a representation in which the gate acts on the carried object through the common-mode
subspace at less than one transport per direction.

## 8. Sources named in the pastes

Chen-Rouze arXiv:2504.02208; Yang arXiv:2609.38007; Bakshi-Liu-Moitra-Tang arXiv:2510.08542; Zhang-Gopalakrishnan
arXiv:2502.13210; Renault arXiv:math/0211369; Fernandez-Gallo-Maillard arXiv:1106.4188; Chazottes-Ugalde
arXiv:math/0211457; Fawzi-Renner arXiv:1410.0664; Dobrushin-type comparison arXiv:1308.4117; Alev-Rao
arXiv:2405.08927; Gibbsian representations arXiv:2001.03880; Bauer-Cvetko-Vah, Stone duality for skew Boolean
algebras with intersections; Bezhanishvili-Bezhanishvili-Gabelaia-Kurz, Bitopological duality for distributive
lattices and Heyting algebras; Bauer-Cvetko-Vah-Gehrke-van Gool-Kudryavtseva, A non-commutative Priestley duality;
Cvetko-Vah, On skew Heyting algebras; Bezhanishvili-Kupke-Panangaden, Minimization via duality. The pastes'
accompanying files (THEORY.md, the zip of finite checks) were not available in this session; their finite checks
are taken as stated and none of them concerns the challenge networks.
