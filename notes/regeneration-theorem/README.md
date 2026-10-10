# The regeneration theorem for the fourth-order pair slices

Working note XLIX. Theory first: two lemmas, one theorem and its cost, then one pre-registered check of the theorem's
prediction against the Monte Carlo truth (committed before the run). No parameter is fitted anywhere.

## 0. Why this object, from facts already established

- **F1 (note XXXIII).** On the system of note XXXII, making the five local slices true at every layer (D3, D21, the
  kappa4 diagonal, the (2,2) and (3,1) slices), with the chain's own mean, variance and covariance, takes the final
  error to 1.8e-9 on network 1 and to the noise floor on network 0.
- **F2 (note XXXIII).** Fed true inputs, the chain's next-layer fourth-cumulant slices are as wrong as when it runs free:
  diagonal 14-33%, (2,2) 16-32%, (3,1) 74-90%. The defect is made afresh by each layer's map.
- **F5 (note XXXIII).** True variance is worth nothing once the slices are true: the covariance program is exact given
  them.
- **Production** builds the (3,1) slice as lambda C_off (the dilation shape only). The truth's (3,1) slice is 45-53%
  dilation-shaped (note XLIV 9k). Single added classes (the hub, the transported single-site class) each explained under
  2% (note XLV 3a).

So the map that regenerates the fourth-order pair slices from the carried state is the one object between the forced
architecture and its measured ceiling. This note derives it completely at leading order.

## 1. The cluster lemma (exact)

Let z be a random vector whose cross-site dependence is a Gaussian coupling C_ac (a != c) plus cross-site cumulant
hyperedges (kappa3, kappa4 of z with at least two distinct sites), and whose site marginals are arbitrary. Let
h_a = relu(z_a). For linear readouts z'_i = sum_a W_ia h_a, put theta_a = sum_k t_k W_(i_k a).

**Lemma 1 (linked clusters).** With K_a(theta; u) = log E exp(theta h_a(z_a + u)) the site cumulant generating
function under a location shift u,

    log E exp(sum_a theta_a h_a) = [ exp(sum_(a<c) C_ac d_(u_a) d_(u_c) + hyperedge operators) exp(sum_a K_a(theta_a; u_a)) ] at u = 0,

and its logarithm is the sum over connected graphs. Concretely, every joint cumulant of the readouts is a sum over
labelled structures:
- a partition of the legs into vertices (several vertices may sit at one site);
- an assignment of vertices to sites (distinct site blocks are distinct sites);
- a connected multigraph of C-edges between vertices at different sites, plus hyperedges.

The weight of a structure is

    prod_vertices [ d^(deg v) kappa_(|v|)(h_site) * prod_(legs in v) W_(row, site) ]  x  prod_edges C  x  prod 1/m!,

with m! for m parallel C-edges and for m incidences of one hyperedge on one vertex.

*Proof.* Write the Gaussian coupling as the heat operator in the site shifts. Expand exp(sum K_a); each factor K_a is a
vertex at site a. The linked-cluster theorem keeps the connected graphs. The vertex factor is the d-th shift
derivative of K_a, whose Taylor coefficient in theta is the joint cumulant of the vertex's legs. A hyperedge
kappa_(a..a b..) with site multiplicities (r, s) enters the operator with coefficient 1 / (r! s!) after counting index
positions. Its r derivatives distribute over the vertices at site a in r! / prod m_v! ways. Labelled legs absorb the
multinomials of theta^p / p!. ∎

Checks by hand:
- kappa(h_a, h_c) gives Mehler's series.
- kappa(h_a, h_a, h_c, h_c) at order C^2 gives (1/2) C^2 d2V_a d2V_c + C^2 (Phi_a^2 d2V_c + Phi_c^2 d2V_a). This agrees with
  the direct moment computation. The two split-vertex terms are exactly what a site-blocked rule would miss.

Here V = Var h, dV and d2V are its shift derivatives, and Phi = P(z > 0), J2 = p(0) and J3 = -p'(0) are the wall jets of
the true marginal. All site factors are shift derivatives of site cumulants of the true marginal; the shift derivative
of E relu^p is p E relu^(p-1).

## 2. Coherence power counting (which structures are leading in the readout metric)

The next layer reads a slice through fresh rows. By the fresh-row theorem (note XLVIII, Theorem 3), the relevant size
of a slice term X_ij is E|X_ij|^2 over the fresh W, at fixed upstream state. Rows of W are i.i.d. N(0, sigma^2) with
sigma^2 = 2/n. Distinct off-diagonal C entries are uncorrelated, of size n^(-1/2).

**Lemma 2.** For a structure with s sites and e C-edges, E|X_ij|^2 scales as n^(f - e - 4), where f counts the free site
indices of the two copies after the forced identifications:
- A site carrying a j-leg is identified across the two copies.
- So is a site carrying an odd number of i-legs.
- So are both ends of an edge of odd multiplicity.
- A site with an even number of i-legs, no j-leg and only even-multiplicity edges is free (counted twice).

*Proof.* Gaussian pairing. A forced cross-copy pairing identifies indices. A free site self-pairs its W_ia^2 to sigma^2,
and its doubled edges to C^2 > 0, so its index sums coherently in each copy. ∎

The maximum of f - e for the (3,1) slice is 1. Every tree of single edges attains it, for every leg pattern. So do the
structures in which a pair of i-legs sits at one site joined to the rest only by a doubled edge (the "doubly-occupied"
structures). Everything else is O(n^(-1/2)) smaller in amplitude.

**Consequence (why single repairs failed).** The leading (3,1) content is a sum of about twenty structures, all of the
same order. Distinct structures are nearly orthogonal (they pair different index patterns), so any one of them alone
explains a few percent at most. This is what note XLV measured. Only the complete sum can reproduce the slice.

## 3. The theorem and its cost

**Theorem (leading-order regeneration).** For each fourth-order pair slice of z' = W h (diagonal, (2,2), (3,1)), the
sum of the structures of Lemma 1 that are leading by Lemma 2 equals the slice up to relative O(n^(-1/2)) in the readout
metric. The site factors are evaluated under the true site marginals, the C-edges with the carried covariance, and the
hyperedges with the carried cross-site slices.

Two classes are not determined by the carried state:
- the transported (2,1,1) and (1,1,1,1) classes of the previous layer's kappa4(z);
- the cross-site kappa3(z) with three distinct sites.

Their size is a statement about memory. Note XXXIII §1 called it small outside the gain mode; this was measured only in
the gain projection.

**Cost.** Each structure is an einsum over at most four site indices and the two rows. All reduce to products of n x n
matrices: M = W Phi C, Q = W^2 dV C, R = (W o M) J2 C, and the doubled-edge forms W^2 (f C o C). The structures that end
in the same right factor (W^T or M^T) are summed before the last product. That gives about 8 n^3 per layer for all three
slices together, roughly 4 units of 2 n^3, against production's 227 units in total.

## 4. Pre-registration (committed before the run)

**Test.** `code/regen4.py`:
- Network 0, transitions 6->7, 10->11, 13->14.
- Inputs: the true layer-l statistics (mu, var, k3, k4, cov, D21, K22, K31) and W_(l+1).
- Structures are enumerated automatically from Lemma 1 (all leg partitions, site assignments, C-multigraphs with at most
  4 edges). They are kept when f - e >= 1 by Lemma 2, plus all structures with one carried hyperedge and at most one
  C-edge.
- Target: the layer-(l+1) truth on the active set (alpha > -2.5), off-diagonal, using the readout weights:
  - (3,1): J3_i J1_j;
  - (2,2): J2_i J2_j;
  - diagonal: unweighted.
- R^2 is noise-free from the two halves. The best global scale is reported. Nothing is fitted.

**Predictions.**
- P1. (3,1): weighted R^2 >= 0.80 at all three transitions, with scale in [0.85, 1.15]. Production's dilation shape is
  45-53%.
- P2. (2,2) and the diagonal: R^2 >= 0.85.
- P3. No single structure carries more than 50% of the (3,1) slice's energy.

**Decision.**
- If P1 holds (and then on network 1), the regeneration map is closed in carried objects, and the next system is a
  clean, fit-free chain built on it.
- If the (3,1) R^2 is below 0.6, the leading list is incomplete. The residual is then examined against the two
  undetermined classes of §3.

## 4a. A caveat recorded before the results (two expansion parameters, not one)

Lemma 2 counts powers of n with C entries of size n^(-1/2). On network 0 the true off-diagonal correlations are larger:
- rms 0.076 at layer 6 and 0.118 at layer 13, against n^(-1/2) = 0.031;
- a collective outlier eigenvalue of the off-diagonal correlation matrix of 22 (layer 6) and 72 (layer 13);
- eigenvalues at -1, from dead neurons.

With an entry scale c, a structure with f free sites and e edges has E|X|^2 ~ sigma^8 n^f c^(2e). So two parameters
order the structures:
- 1/n per lost free site (coherence);
- c^2 (about 0.06 here, in units of the site variance) per extra edge that adds no site (a loop or multi-edge).

The leading list of §4 is maximal in f - e. Within it, structures with more tree edges dominate, as (n c^2)^e.
The first omitted group is the one with the same free sites and one extra loop edge (e = 4, f - e = 0: 93 structures
for (3,1)). It is suppressed only by c^2 per structure, so as a group it may carry tens of percent of the energy.

The second run (EXPO_MIN = 0, network 0, (3,1)) measures that group directly. A leading-only shortfall of up to about
20% is therefore not a failure of Lemma 1. It is the c^2 tail, and the decision of §4 is read with this group included.

## 4b. Engine check on an exact toy (before the network results)

The one-step run first failed (outputs/regen4_net0_L13_coincidence_bug.txt: (3,1) R^2 -8%). The cause was an
implementation error, not the theory:
- the einsums summed over all site tuples;
- non-adjacent site coincidences therefore reproduced the doubly-occupied structures, which are leading and are also
  enumerated on their own, so they were counted twice.

The engine now sums over distinct site indices exactly, by Moebius inversion over coincidence patterns.

Two checks of the corrected engine:
- **Covariance and variance of z' (two legs).** The enumeration reproduces the layer-14 truth from layer-13 inputs with
  R^2 100.000% and scale 1.0001 (outputs/validate_low_net0_L13.txt). This checks the site factors, the edges and the
  weights.
- **Fourth cumulants on an exact Gaussian toy (code/toy_check.py).** Setup: n = 6 sites, correlations about 0.1, all
  structures up to 4 edges, against a Monte Carlo of 4e8 samples. Every entry agrees to the Monte Carlo noise
  (outputs/toy_check_seed1.txt):

| entry | engine | Monte Carlo |
|---|---|---|
| K31[0,1] | -0.007823 | -0.007825 |
| K31[2,0] | +0.017719 | +0.017717 |
| K22[0,1] | +0.021897 | +0.021901 |
| k4[1] | +0.055611 | +0.055583 |

So Lemma 1 and its implementation are verified exactly. The network run tests only two things: the truncation of
Lemma 2, and the classes of §3 that the state does not carry.

## 5. The same theorem as a TAP / Plefka expansion (the compact form for production)

Write C = D + C_o, with D the site variances. The layer's readout generating function is

    G(theta) = log [ exp( (1/2) d_u^T C_o d_u ) exp( sum_a K_a(theta_a; u_a) ) ]_(u=0).

Its connected expansion splits by loop number.
- **Tree level is a saddle point.** G_tree(theta) = stat_u [ sum_a K_a(theta_a; u_a) - (1/2) u^T C_o^(-1) u ], with
  u = C_o grad_u K. Site labels are summed independently, so trees that revisit a site, the split-vertex structures,
  are automatically included.
- **The first loop correction is -(1/2) log det(I - C_o H(theta))**, with H = diag(d_u^2 K_a(theta_a; u_a)). Its
  two-cycle term (1/4) tr((C_o H)^2) is exactly the doubled-edge (doubly-occupied) family. The longer cycles are the
  c^2-suppressed loop tail of §4a.

So the leading regeneration is the tree level plus the two-cycle of the TAP free energy of one layer. The cumulants of
the readouts are the theta-derivatives along s W_i + t W_j. Expanding the saddle to third order in theta gives the
structure sums as a few matrix products per layer. The production form needs no enumeration.

## 6. A second prediction, written before the results: the memory of the four-site class

Section 3 left the transported (2,1,1) and (1,1,1,1) classes of kappa4(z^l) undetermined. Power counting sizes them.

**Size.**
- The previous layer's leading trees create distinct-site entries kappa4(z^l)_abcd of the same size as its (3,1)
  entries: E X^2 ~ sigma^8 n^4 c^6 J^8.
- Transport by T = W Phi on all four legs multiplies the energy by (sigma^2 n mean Phi^2)^4 ~ 0.9^4 ~ 0.66 per layer of
  age. This uses sigma^2 n mean Phi^2 = 0.9, which is note XLVII's one-step gain.
- The transported four-site class therefore enters K31' at the same order as the newborn structures, discounted by
  0.66 per age. Births from different layers are nearly orthogonal (stage 15, Theorem 3.2).

**Consequence.** Summed over ages, the memory carries about 0.66 / (1 - 0.66) ~ 1.9 times the newborn energy. A
one-step regeneration from the carried slices alone would then reach only about 1 / 2.9 ~ 35% of the (3,1) slice.

**Prediction P4.** Leading-only (3,1) R^2 lies in 0.2-0.5, and the loop tail of §4a adds less than 0.1. If so, P1 fails
for a reason the theory states: the regeneration must be fed the four-site memory. In the chaos state that memory is
not a new tensor:
- the (1,1,1,1) and (2,1,1) classes of kappa4(z^l) are, at leading order, pairs of second-chaos sources contracted
  through the birth Gram (the path or hub class);
- and third-chaos sources (the stars), each a vector pair per neuron per layer, the same storage as the kappa3 sources.

So the design object is the chaos-state formula for the fourth-order slices, read from the carried sources. The
one-step slice map is not.

If the (3,1) R^2 is instead at or above 0.8, the memory is small, and the one-step regeneration of §3 is the design.

## 7. Result (network 0; outputs/fast31_net0_L*.txt, outputs/regen4_net0_L13.txt)

**The engine was rebuilt for speed.** The first fleet runs evaluated each structure as one multi-operand einsum. Its
intermediate steps ran in numpy's scalar loop, at about 1 GFLOP/s and on one thread, so each slice took 25 minutes or
more. code/regen4_fast.py compiles every structure into BLAS matrix products by leaf elimination:
- it agrees with the einsum engine to 7e-15 on all 4806 structure terms (code/equiv_fast.py);
- one slice now runs in 17 s;
- two kinds of term are handled apart. The (i,j)(i,j) structures of the (2,2) slice are n^4 as full matrices and are
  evaluated on 8000-20000 sampled pairs. The 12 cyclic coincidence terms per slice (path ends merged into a triangle)
  have f - e = 0 by Lemma 2 and are skipped.

(3,1) slice, leading structures from the true layer-l state, readout-weighted, noise-free R^2:

| transition | all leading | scale | C-edges only | dilation shape alone (fitted on the target) |
|---|---|---|---|---|
| 0 -> 1 | **100.2%** | 1.006 | 100.2% | 48.9% |
| 1 -> 2 | 56.3% | 1.11 | 36.0% | 49.6% |
| 2 -> 3 | 42.2% | 1.03 | 14.2% | 51.1% |
| 4 -> 5 | 29.3% | 0.95 | 4.4% | 51.9% |
| 6 -> 7 | 19.5% | 0.85 | -1.8% | 46.5% |
| 10 -> 11 | 8.6% | 0.70 | -4.9% | 42.9% |
| 13 -> 14 | 1.2% | 0.53 | -8.3% | 40.0% |

**What this settles.**
- **The theorem is exact where its hypotheses hold.** At 0 -> 1 the layer is exactly Gaussian (z^0 = W_0 x), so no
  cross-site cumulant beyond C exists. There the leading structures of Lemmas 1-2 reproduce the truth's (3,1) slice to
  the Monte Carlo noise, at scale 1.006, on the real network with its real correlations. The truncation of Lemma 2 is
  therefore adequate, and the c^2 loop tail of §4a is negligible.
- **P1 fails at every deeper transition, and it fails as P4 predicted, only more strongly.** The shortfall grows
  monotonically with depth, from 44% at layer 1 to 99% at layer 13. That is the profile of accumulated memory:
  - the cross-site third cumulants of z^l with three distinct sites;
  - the (2,1,1) and (1,1,1,1) classes of kappa4(z^l).

  The carried slices do not contain these, and every layer adds to them. P4 estimated memory at 1.9 times the newborn
  energy (one-step R^2 about 35%). The measured deep-layer values, 1-9%, mean the memory is about 10-100 times the
  newborn energy there.
- **The C-only part becomes anti-correlated at depth.** The newborn Gaussian-closure content is a small piece pointing
  against the transported history. This is the regime where any one-step or Gaussian-closure regeneration, including
  production's lambda C_off and every single-class repair of notes XLIII-XLV, must fail. It explains F2 (74-90% error
  with true inputs) from first principles.

**Consequence for the design.** The fourth-order pair slices cannot be regenerated from the layer's slices. They must
be read from a carried history, exactly as the third-order slices are read from the kappa3 sources. In the chaos
state that history consists of objects of the same kind as the kappa3 sources (§6):
- second-chaos pairs contracted through the birth Gram (path and hub);
- third-chaos sources (stars), primary (born with the wall's third jet) and secondary (the wall's second jet times the
  incoming second chaos);
- second-chaos loops.

The next theory item is the chaos-state formula for all three fourth-order slices, with its cost. Its check is the same
ladder: it must hold at every depth, not only at 0 -> 1.

## 8. Pre-registration: the memory ladder (committed before the run)

**Claim (multilinearity).** Each leading structure is multilinear in the rows of W, one row factor per leg. So the
transport of layer m's newborn fourth cumulant into the slice of layer l+1, by first jets on every leg, is the same
structure sum with W_(m+1) replaced by the transported matrix

    Wt_m = T_(l+1) T_l ... T_(m+2) W_(m+1),     T_k = W_k diag(Phi^(k-1)).

Site factors, C-edges and slice hyperedges are those of layer m.

**Test (code/regen_mem.py).**
- Ladder: sum_(m <= l) newborn_m(Wt_m), with nothing fitted.
- Network 0, (3,1) slice, transitions 2->3, 6->7 and 13->14. Same metric as §7.
- Reported with and without a best-scaled dilation shape. The dilation's own memory is the gain mode, which the
  ladder may or may not reproduce.

**Prediction P5.** The ladder reaches R^2 >= 0.6 at 6->7 and 13->14 (one-step: 19.5% and 1.2%), with scale in
[0.7, 1.3]. If it falls short, the residual memory is the part this transport omits:
- the dressing of coincident legs at intermediate layers;
- the kappa3 -> kappa4 conversion through three-site kappa3 hyperedges.
