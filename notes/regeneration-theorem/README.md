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
