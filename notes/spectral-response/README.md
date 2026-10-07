# The spectral response and the Ward identity, tested on the ledger

Working note XXII. A 46-page theory note was uploaded from a separate session (`ncg_cumulant_theory.pdf`, checkpoint
NCG-20261007-A, written against note XIX without access to the code, the ledger arrays or the networks). It derives:
an exact scalar identity for the three omitted fourth-cumulant classes and the complete influence polynomial of
their first variation (Theorem 2); the exact response of that quantity to a pure Gaussian third-cumulant tangent,
in score form and in gate-facet form (Proposition 3, Theorem 4, Corollary 5); two analytic witnesses; a Gaussian
Bott-Dirac spectral triple of index one in which the response is a matrix element of the third Hermite sector, the
sector the index cancels (section 3); the exact scale Ward identity with reference drift (Theorem 4.4) and the
statement that the scale tangent of a Gaussian reference has no cumulant component above four (Corollary 4.3); the
compression-defect and closure criteria (sections 5-7); the adjoint error identity (section 8); and a conditional
factor-Gaussian evaluation of the response at O(n^2) per outer node (section 9). It ends with what a theorem about
the deep profile would still need (section 10.3).

Two of its results are testable on our ledger, and both were tested here. The rest is placed against the numbers
of notes XIX-XXI.

- **Its two-neuron witness is exact.** The closed form Delta_T(a, b) = (3 t c / 2) a b [a^2 (1 + 4 c^2) + b^2 (1 - 4
  c^2)] was recomputed by factorised adaptive quadrature with recentering (`code/twoneuron_exact.py`): relative
  difference 2e-16 at three rows, and the two per-entry variations 3c/8 +- 3c^3/2 to ten digits. The influence
  polynomial and its gate-facet reduction are right where they can be checked.
- **The Ward identity closes the ledger's rows.** Note XIX reported that the unit-gain transport of a scale mixture,
  computed from the third- and fourth-cumulant slice variations at fixed (mu, S), sums to 0.93 in the third-cumulant
  row and 0.83 in the fourth at depth, and guessed the gap was the fourth-order truncation. The theory note says
  the scale tangent has three more parts, the reference drift (mu -> mu - mu / 8, S -> S + mu mu^T / 4 per unit
  gain) and nothing above order four, and that the correct check is drift + both channels = the Ward derivative.
  Measured (`code/ward.py`, `outputs/ward_off0.txt`): the drift adds +0.04 to +0.07 to the third-cumulant row and
  -0.015 to +0.005 to the fourth; with the dropped output classes of note XXI added to the fourth row, the sums are
  1.04 and 1.05 at layer 2, 1.01 and 1.02 at layers 5-7, 1.00 and 0.98 at layer 8, and 0.96-0.99 and 0.96-0.98 at
  layers 10-15. Both rows sum to 1 within 3-4% at every depth. The guess in note XIX is withdrawn.
- **The exact omitted response is, on these networks, at the noise floor.** The theory note's Delta_T is the
  first variation of the three omitted classes under the incoming pair-slice third-cumulant source. Note XXI
  measured those classes directly and attributed them: at targets 8, 11, 15 the omitted classes are 0.00078,
  0.00121, 0.00337 in g_4 units, the mixture's own transport accounts for 0.00083, 0.00158, 0.00309, the one loop
  for -0.0002 to -0.0001, and what is left for a third-cumulant response is +0.0001, -0.0003, +0.0004, within the
  noise of the (1+1+1+1) class. So Delta r_43 through the omitted classes is at most 2% of g_4 here. The
  third-to-fourth feed that matters (0.0023 against 0.0004 at layer 14, note XXI section 5) is in the retained
  classes, and in the theory note's own decomposition (C.2) it is leakage, not coefficient error: the true
  third-cumulant slices have a component outside the mixture-shaped retained direction.
- **The spectral placement is the Mehler order.** The response lives in the third Hermite sector of the
  Ornstein-Uhlenbeck number operator; in the chain that sector is the Hermite-index-3 programs of the term table
  (the (3,0), (2,1), (1,2), (0,3) Mehler terms of `pairvar.Layer.var3`), which note XVIII section 1 already
  identified as the spectral decomposition of the heat operator. The index cancels the sector: a clean proof that
  the Connes route gives the frame and not the amplitude, as note XVII concluded from the absence of a mobility gap.

## 1. The omitted classes, in the theory note's form and in ours

The theory note's Proposition 1 writes the sum of the three omitted classes as K_miss = E[S_1^4 - 3 S_2^2 + 2 S_4]
- 3 (V^2 - D^2) + 6 R with S_r the power sums of u_i = w_i (H_i - E H_i), V = 1^T C 1, D = tr C and R the ordered
sum of squared off-diagonal covariances. That is the contracted form of note XXI's class decomposition: (4) and
(2,2) are what the retained power sums S_4 and S_2^2 see, 3 D^2 + 6 R is their Gaussian pairing, and the rest is
the three omitted classes. Note XXI's `k4classes.py` computes the same thing class by class with the five power-sum
polynomials, on both networks, and its Gaussian pairings per class are the inclusion-exclusion split of 3 V^2 -
3 D^2 - 6 R. The two are the same identity; ours resolves it into its three parts, which is what showed the (2+1+1)
class to be the whole of it.

Theorem 2's influence polynomial F(u) = P(u) - b^T u - 6 V S_1^2 + 6 D S_2 + 12 u^T C_off u, with the recentering
and covariance terms, is what our first-variation programs compute implicitly for the pair classes: `var3` and
`var4` differentiate the pair moments and then re-form the cumulants (`cum_pair`), which is the same bookkeeping.
For the omitted classes we never formed the influence; we measured the classes on the true law and on the Gaussian
reference and took the difference. Where the two routes can be compared, the two-neuron witness, they agree to
machine precision (`outputs/twoneuron_exact.txt`):

| quantity | quadrature | closed form |
|---|---|---|
| delta cum(H1, H1, H1, H2) | 0.2448438091 | 3c/8 + 3c^3/2 = 0.2448438091 |
| delta cum(H1, H2, H2, H2) | 0.0543629012 | 3c/8 - 3c^3/2 = 0.0543629012 |
| Delta_T at w = (1, 1) | 1.1968268412 | 3c (1 + 4c^2) + ... = 1.1968268412 |

## 2. The Ward identity on the ledger

**What was missing from note XIX's accounting.** The ledger's transport of a unit-gain scale mixture applied the
mixture's third- and fourth-cumulant slices (`pairvar.sm_slices`) to the Gaussian reference at the true (mu, S) and
read the coherent third and fourth cumulants of the next pre-activation. Corollary 4.3 of the theory note says the
scale tangent of a Gaussian reference is exactly four things: a mean shift -mu / 8, a covariance shift mu mu^T / 4,
the third-cumulant direction and the fourth-cumulant direction, with every cumulant above four having zero first
variation. The ledger had the last two. Theorem 4.4 says the sum of all four, applied to the output's fourth
cumulant, equals the scale derivative of that cumulant, which for a positively homogeneous output is exact. The
drift was computed here by central differences of the Gaussian one loop in the direction (-mu / 8, mu mu^T / 4)
(`code/ward.py`, step 0.05 per unit gain), projected on the same patterns as the ledger, and added to the raw rows;
the fourth row also receives the mixture's dropped output classes from note XXI. Network 0:

| l+1 | r33 + r34 (pair) | drift | third-cumulant row | r43 + r44 (pair) | drift | dropped output classes | fourth-cumulant row |
|---|---|---|---|---|---|---|---|
| 2 | 1.005 | +0.036 | 1.041 | 1.041 | +0.004 | +0.005 | 1.050 |
| 4 | 0.960 | +0.062 | 1.022 | 1.000 | +0.004 | +0.025 | 1.029 |
| 6 | 0.937 | +0.071 | 1.008 | 0.978 | +0.003 | +0.036 | 1.017 |
| 8 | 0.944 | +0.053 | 0.997 | 0.912 | -0.000 | +0.065 | 0.977 |
| 10 | 0.926 | +0.059 | 0.985 | 0.910 | -0.006 | +0.086 | 0.990 |
| 12 | 0.928 | +0.048 | 0.976 | 0.865 | -0.009 | +0.109 | 0.965 |
| 14 | 0.916 | +0.047 | 0.963 | 0.855 | -0.015 | +0.141 | 0.982 |
| 15 | 0.933 | +0.038 | 0.971 | 0.834 | -0.015 | +0.144 | 0.963 |

- **The third-cumulant row.** The drift is the mixture's inflation of the covariance along the mean direction seen
  as a change of the Gaussian reference; it carries 4-7% of the unit transport, and with it the row is 1.00 +- 0.04
  at every depth. The 2-3% below 1 at layers 10-15 is the dropped (1+1+1) output class (note XIX: the true-h pair
  restriction is 2% low there). This is the same fact note XXI section 5 found from the other side: the one loop on
  the marginal state re-counts the drift, so a recursion on the marginal state must not add it twice, and one on
  the conditional state must.
- **The fourth-cumulant row.** The drift is nil (the covariance shift along mu mu^T enters the fourth cumulant only
  at second order in the gain), and the whole gap is the dropped output classes, which note XXI derived in closed
  form; with them the row is 1.00 +- 0.04. Homogeneity holds on the ledger to that precision.
- **What the O(g) terms of the Ward derivative are.** Theorem 4.4 gives d kappa_4 / dt = 3 v^2 + 3 m kappa_3 +
  kappa_4 and the normalised-gain derivatives 1 + (1/2 - alpha^2 / 4) g_3 and 1 + 1.5 g_3 alpha^2 + g_4 (1 -
  alpha^2 / 2). These are corrections to the increment of gain added by a fresh scale, of relative size 10-40% of
  that increment on the sigma^4-weighted neurons; they are not corrections to the transport of the gain already
  present, which homogeneity fixes at exactly 1. At depth the increments are 0.001 per layer, so the second-order
  terms are 1e-4 per layer and the first-order ledger is the right object. The theory note states the distinction
  ("the full scale tangent must be distinguished from the ledger's pair-restricted incoming tangent"); the numbers
  above are its test.

## 3. The response, placed

**The exact omitted response is small here.** Corollary 5 gives Delta_T for the ledger's incoming pair source
(tau_i on the diagonal, Gamma_ij on the (i,i,j) entries) as bulk terms with products of three gate indicators plus
single gate facets with the boundary density p_i(0), with no double-facet terms, and Corollary F.4 evaluates it in
closed form at the independent standard reference: diagonal sources contribute nothing and the pair source gives
(3c / 2) sum_{i != j} Gamma_ij w_i w_j [w_i^2 (1 + 4 c^2) + w_j^2 (1 - 4 c^2)], a sum over odd powers of a random-sign
weight row, hence incoherent unless Gamma has structure. For the mixture's Gamma_ij = g (mu_i C_ij + sigma_i^2 mu_j
/ 2) the coherent part is sum_i w_i^3 sigma_i^2 times mu_z, which is O(n^-1/2). Note XXI's attribution measured
the same thing without the formula: the omitted classes on the true law minus the mixture's transport minus the
one loop leave +0.0001, -0.0003, +0.0004 at targets 8, 11, 15, the noise of the (1+1+1+1) class. So on these
networks Delta r_43 through the omitted classes is bounded by 2% of g_4 and is not where the deep profile comes
from. The theory note is explicit that it does not claim otherwise ("a correct one-step formula [must not be]
misreported as a complete derivation of the measured depth profile").

**Where the feed is, in the theory note's language.** Equation (C.2) splits a missed response into coefficient error
on the retained space and leakage into its complement, and says a probe inside the retained space tests only the
first. Note XXI's probe was the true third-cumulant slices of y, and they feed the next fourth cumulant 2-6 times
more than the mixture-shaped slices at the fitted gain, while feeding the next third cumulant exactly as the mixture
does. In (C.2) that is the second term: the physical perturbation has a component outside the two-gain retained
space (the two patterns of the (2,1) slice carry different amplitudes at depth, 0.0164 and 0.0234, equal for a
mixture), and it is observed by the fourth-cumulant query and not by the third. The theory note's caution that a
late failure of the rollout "alone proves neither E != 0 nor that no two-coordinate closure exists" is right as
stated; the identification of the leaking direction in note XXI is what turns the failure into a witness of E != 0
for that probe. Its Theorem 3 then says the minimal augmentation is rank E, and the measured witness is one
direction (the second (2,1) pattern), so a three-coordinate state (two third-cumulant amplitudes and one fourth)
is the next closure to test, if one wants a closed recursion rather than the chain, which carries the full slice.

**The factor-Gaussian evaluation and the budget.** Section 9 evaluates Delta_T at O(n r + n^2) per outer node per
row when Sigma = diag + B B^T exactly with r factors. Our pre-activation covariance is not of that form exactly (the
off-diagonal correlations have rms 0.07-0.12 across all pairs), and approximately it needs r = 64 for 91% of the
off-diagonal Frobenius mass (note XX). With an r-dimensional outer integral of a few hundred nodes and n rows, the
cost is hundreds of n^3 per layer, the same wall note XX priced for the rank-K gate covariance (K transports per
source-layer) and the same reason: the common-mode subspace has dimension 64, not 4. As an offline derivation
tool at a few rows it is fine, and the quantity it would compute is the one note XXI already bounded at 2% of g_4.

**The Connes construction.** The Gaussian Bott-Dirac operator D = sum (eps_j a_j + iota_j a_j^*) with D^2 = N_b +
N_f, compact resolvent, index one, and the response as <S_T, Pi_3 I_f>, the pairing of the source score (a third
Hermite polynomial, N_b-eigenvalue 3) with the third-sector component of the influence: this is the statement, with
a proof, that the response is a Mehler-order-3 matrix element of the Ornstein-Uhlenbeck semigroup, and that an
index, which cancels every positive-energy sector, cannot carry it. In the chain the third sector is the
Hermite-index-3 Mehler terms of the pair programs, and the response to a pair source is what `var3` computes for
the pair classes. The gate regularity result (the gate indicator is in Dom(N_b^{s/2}) only for s < 1/2, ReLU for
s < 3/2) is the spectral form of note XVII's "the Hermite coefficients of the gate decay algebraically": the reason
the Mehler truncation at order 2 costs 1e-4 of the covariance (note XVIII) and not less. The thermal state
sigma_s^{omega_beta} and the remark that the Mehler semigroup is not a modular automorphism group correct note
XIX section 4's shorthand "the modular flow is the Ornstein-Uhlenbeck semigroup": the OU semigroup is the heat
operator of the Gaussian Dirac complex; the modular flow of the thermal state is a different, unitary object. Note
XIX's point (that the KMS/groupoid characterisation adds nothing computable) stands; its wording is now wrong by
one identification and is left as recorded, with this correction.

**The adjoint identity and the chain.** Section 8's identity c^*(x_L - x^_L) = sum v_{l+1}^T delta_l, errors
weighted by their remaining readout sensitivity, is the mechanism behind two of our measurements: the derived
two-channel gain on the conditional state is 27% low in g_3 and gives the best final mean (note XIX section 2.4),
and the chain's fitted fourth-cumulant table beats every physically more accurate replacement (notes XVIII, XXI
section 2.5). The theory note draws the conclusion we drew: replacing a compensating coefficient by an isolated
correct one need not improve the final prediction until the other terms it compensated are accounted for. The
proposed chain experiment of note XXI section 5 is framed that way: it replaces the table by all three of the
pieces it stands for, not by one.

## 4. What transfers

| from the theory note | what it is here | number or decision |
|---|---|---|
| contracted identity for the omitted classes (Prop. 1) | note XXI's class decomposition in summed form | same identity; ours resolves it into (3+1), (2+1+1), (1+1+1+1) |
| influence polynomial and gate-facet response (Thm 2, Cor. 5) | the first-variation programs, extended to the omitted classes | two-neuron witness verified to 2e-16 |
| scale tangent has no cumulants above four (Cor. 4.3) | note XIX's "truncation" guess withdrawn | drift 4-7% on the third row, 0 on the fourth |
| Ward identity with drift (Thm 4.4) | rows of the ledger sum to 1 with drift and dropped classes | 1.00 +- 0.04 at every depth, both rows |
| exact omitted response Delta_T | bounded by the measured remainder of note XXI | at most 2% of g_4 on these networks |
| coefficient error against leakage (C.2) | the non-mixture shape of the third-cumulant slices is leakage | witness: 0.0164 against 0.0234 at layer 15; rank-1 augmentation to test |
| factor-Gaussian evaluation at O(n^2) per node (Thm F.3) | needs r = 64 here; hundreds of n^3 per layer | the same wall as note XX, same cause |
| Bott-Dirac triple, response in the third sector, index cancels it | the Mehler-order-3 programs; the frame, not the amplitude | relabelling with a proof |
| gate regularity s < 1/2 | algebraic Hermite decay of the gate (note XVII) | the Mehler truncation floor of note XVIII |
| modular flow is not the Mehler semigroup | note XIX section 4 corrected | - |
| adjoint error identity (sec. 8) | why a wrong gain can give the best mean, and why the table beats its replacements | notes XIX 2.4, XVIII, XXI 2.5 |
| what a deep-profile theorem needs (10.3) | exact closure fails by one measured direction; the chain carries the full slice | the three-coordinate closure is the test; the chain experiment of XXI is the decision |

## 5. A flop lever proposed from the theory, tested and closed

In the discussion that followed this note a lever was proposed: the coherent part of every source is now known in
closed form (three numbers per layer with the mixture and skew shapes), so if that part occupies many of the ranks
the legs carry, stripping it before compression would lower the rank ladder and save transport flops at equal
error. The estimate offered was 70-100 of the 384 shared ranks, i.e. 20-25% of the flops. The test
(`code/coherent_strip.py`, `outputs/coherent_strip_off0_l10.txt`) forms the (2,1) slice of each source's
third-cumulant tensor from the live legs of layer 10 (uncompressed, correlation 0.994 with the chain's own carried
D21), projects it on the two coherent patterns mu_a S_ac and sigma_a^2 mu_c, and compares the singular spectra of
the slice and of its residual:

| object | coherent share (R^2) | ranks for 50 / 90 / 99% of the Frobenius mass | after stripping the coherent part |
|---|---|---|---|
| all sources, (2,1) slice from the legs | 0.436 | 1 / 21 / 190 | 1 / 48 / 251 |
| chain's carried D21 | 0.442 | 1 / 17 / 220 | 1 / 47 / 308 |
| per source, ages 10 to 1 | 0.09-0.31 | 1 / 22-107 / 91-335 | 1-9 / 34-116 / 107-345 |

The coherent part is the top singular direction and nothing else: it carries half the slice's mass in rank one (the
mean direction), and removing it leaves a flatter spectrum that needs more ranks, not fewer, for any given share
of the mass. Carrying the coherent sector analytically would save about one rank per source. The lever is closed;
the estimate of 20-25% is withdrawn. This is the same fact as note XVII's "the mean direction has a flat spectrum
in the shared basis" read correctly: the ranks are spent on the incoherent, per-neuron content, which is the content
the output needs, and the coherent content costs nothing to carry. What the theory of notes XIX-XXII derived in
closed form is, in the legs, one direction.

Network 1's Ward accounting (`outputs/ward_off1.txt`) matches network 0: with the drift the third-cumulant row is
1.04 at layer 2 falling to 0.99 at layer 11 and 0.92-0.99 at layers 12-15, and the fourth-cumulant row with the
dropped output classes is 1.05 at layer 2 and 0.98-1.00 at layers 6-8.
