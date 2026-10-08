# The deck symmetry and the latent geometry of the cumulant slices

Working note XXXVII. A reading of translation-surface dynamics (branched covers and their deck groups, the
Kontsevich-Zorich cocycle and the deviation of ergodic averages, the wind-tree diffusion rate as a Lyapunov exponent
of one isotypic piece, Masur-type recurrence criteria, the five-distance theorem, genericity on curves) and of the
algorithmic theory of random geometric graphs, turned into falsifiable statements about the chain and tested on Monte
Carlo truth (networks 0 and 1, 1.6e7 inputs, `mc2`) and, where a statement implied a rule, in the chain (32 networks).

What the reading predicts, in the chain's terms:

1. **A deck group.** relu(z) = z + relu(-z) is the Moreau decomposition onto the orthant and its polar cone: two sheets
   over each unit, exchanged by alpha -> -alpha. The hinge coefficients are deck-invariant, so whatever a saturated unit
   contributes through its hinge is the same on both sides (section 1).
2. **A cocycle with a tautological direction and a deviation subspace.** In the KZ cocycle the top exponent belongs to
   the tautological plane, which is computable exactly; deviations of ergodic sums live in a finite-dimensional
   unstable subspace. Here the tautological direction is the gain (scale) mode, exact by homogeneity; the prediction
   is that what the chain misses lies in the top covariance eigenspace, not in random directions (sections 2, 3, 5, 6).
3. **Latent positions.** In a random geometric graph the k-point statistics are functions of latent positions and the
   elementary signals are signed motif counts. Here higher cumulants should be multilinear forms of the units' collective
   loadings plus the gain, and the fresh part should be the tree motifs (stars and paths) of the inter-layer graph
   (sections 4-6).

What the tests say:

- **The deck symmetry is exact and gives no saving beyond the off side.** The symmetric drop rule costs +59% to +156%
  in the chain: the certificate covers hinge content, and an on-saturated row is almost all linear content.
- **The gain is nearly all of the fourth-cumulant diagonal and (2,2) slice (94-98.5% of their energy), and the chain's
  error on them is exactly the rest.** The chain has the gain amplitude right (0-3% of its error is gain-shaped); what it
  lacks is the non-gain part. The (3,1) slice is only 22-65% gain (22-49% in the read metric), and the chain is 73-94% wrong on it; the Monte Carlo
  oracle values that slice at -21.9% (network 0) and -6.0% (network 1) of the final MSE.
- **The fresh part is the tree motifs, and it is small at depth.** The one-step star and path births (formulas verified
  to three digits against Monte Carlo) add 7-12% to the explained share of the true (3,1) slice at layers 2-5 and 0.4-5%
  at layers 9-15. The non-gain content at depth is accumulated, not fresh.
- **The accumulated non-gain content is collective.** On its single index, 76-85% of the non-gain (3,1) energy at layers
  10-15 lies in the top-64 covariance eigenspace (6% for a random 64-dimensional subspace), 94-98% in the top 256; the
  (2,1) slice likewise.
- **Units are not conditionally independent given a 64-dimensional collective state.** The Monte Carlo split of the
  true (3,1) slice into collective, residual and cross parts (section 6) shows that the collective fourth cumulant
  explains the non-gain content only from 6% (layer 3) to 55% (layer 15); the rest is a collective-residual cross term
  whose column index is collective but whose row dependence is unit specific (one shared matrix gain explains 4-10%).
  A pure collective-mode propagation is therefore not closed at any rank we can afford.
- **The output-metric calibration of note ray-compiler section 8 holds out of sample:** -3.76% raw on the 50 held-out
  networks (46/50 better), against -3.07% predicted.

## 1. The deck symmetry

**Statements** (`code/deck_check.py`, `outputs/deck_check.txt`; Gauss-Hermite quadrature at the Gaussian reference):

- *Hinge parity.* With c(1, k) = E[d^k relu(z) / dz^k] at z ~ N(mu, s^2), alpha = mu / s: c(1, k)(alpha) =
  (-1)^k c(1, k)(-alpha) for k >= 2, since every derivative of order >= 2 of relu(z) and of its polar part relu(-z)
  agree; the transmission c(1, 1) = Phi(alpha) is not deck-invariant (Phi(alpha) + Phi(-alpha) = 1). Exact to the
  printed precision; finite differences of the closed form E relu = mu Phi + s phi reproduce c(1, 2) and c(1, 3).
- *Polar cumulants.* The cumulants of order >= 2 of an on-saturated unit's polar part relu(-z) equal those of the
  mirrored off-saturated unit relu(z') at -alpha (variance 1.2e-3, kappa_3 1.0e-3, kappa_4 1.1e-3 at alpha = 2.5, s = 1).
- *Birth split.* For (z_j, z_k) jointly Gaussian, kappa(y_j, y_j, y_k) = kappa(z_j, z_j, y_k) + 2 kappa(z_j, v_j, y_k) +
  kappa(v_j, v_j, y_k), v = relu(-z) (residual 1e-16). The first term has its nonlinearity at k; the j-hub terms carry
  only v_j, so they vanish like alpha_j's tail on both sides: 1.5% and 3.0% of the total at alpha_j = +3.0 and +2.5.

**Measurements** (`code/deckspec.py`, `outputs/deckspec_off{0,1}.txt`, chain dumps). The alpha distribution is symmetric:
at |alpha| >= 2.5 about 4% of units per side at layer 4, 10-12% at layers 7-8, 16% at layer 10, 20-22% at layer 14.
Rows with |alpha| >= 2.5 carry 1e-5 to 3.8e-4 of the read-weighted energy of the (2,1) slice.

**The chain test** (`V48_DECK_ROWS`, `V48_DECK_BIRTH` in `est_v29.py`, 32 networks, paired against the V35 candidate,
mean raw 2.278e-8; `outputs/deck_rule_32nets.txt`):

| at both sides, \|alpha\| >= 2.5 | raw change | better |
|---|---|---|
| drop the (2,1) slice's rows | +58.6% +- 3.0 | 0/32 |
| drop the birth columns | +64.0% +- 2.8 | 0/32 |
| both | +156% +- 5.6 | 0/32 |

The certificate covers hinge content only. The chain's (2,1) rows and its leg columns also carry linearly transmitted
content weighted by Phi (the pass-through of the sources and the separable part of the A leg), which is the whole row
of an on-saturated unit. The deck map keeps |c(1, k)| for k >= 2 and sends Phi to 1 - Phi, so only the off side is free,
and notes XXIX and XXXI already take it. The symmetry is exact and closed as an algorithmic lever; it is kept as a
check (any proposed hinge-level rule must be deck-symmetric).

## 2. The (2,1) slice is collective, and so is the chain's error on it

(`outputs/deckspec_off{0,1}.txt`, `outputs/denoise_off{0,1}.txt`.) In the read metric R = diag(c(1,2)) D21 diag(Phi), the
singular energy at rank 32 / 64 / 128 / 256 / 512 is 0.970 / 0.986 / 0.996 / 0.9996 / 1.0 at layer 14 (0.53 / 0.65 /
0.80 / 0.93 / 0.99 at layer 1). The column (linear) index aligns with the gated covariance diag(Phi) C diag(Phi): its
top-256 eigenvectors capture 99.3-99.7% of R at layers 8-14. The squared index does not (29-71%), as a latent model
predicts: it enters through v_i (x) v_i, whose span has dimension K(K + 1)/2.

The chain's error on R (network 0) is 0.0037-0.0039 of the truth's energy at layers 4-14 (Monte Carlo floor 1e-4 to
3e-4), and 66% (layer 2) to 87% (layer 14) of it lies inside the top-256 span. Projecting the chain's R on that span
never lowers the error and is neutral at rank 512. The error is not high-rank noise a projection could remove; it is a
wrong collective component.

## 3. The gain share of every true slice

`code/latent_slices.py` fits, per layer and slice, the scale-mixture shape (z = G z~, Var G^2 = g, first order in g):
k3 = 1.5 g mu var, D21_ij = (g/2)(2 mu_i C_ij + mu_j var_i), k4 = 3 g var^2, K22_ij = g (var_i var_j + 2 C_ij^2),
K31_ij = 3 g var_i C_ij. "Read" weights rows and columns by the Hermite coefficients the pair program applies
(c(1,2) x Phi for D21, c(1,3) x Phi for K31, c(1,4) for k4). Network 0 (`outputs/latent_slices_off{0,1}.txt`; network 1
agrees to the second digit):

| layer | g | k4 gain share (read) | K22 | K31 plain / read | D21 read | k3 read | chain read error: g4 / K22 / K31 / D21 / D3 |
|---|---|---|---|---|---|---|---|
| 2 | 0.006 | 0.985 | 0.981 | 0.37 / 0.37 | 0.56 | 0.73 | 0.12 / 0.14 / 0.85 / 0.059 / 0.043 |
| 5 | 0.010 | 0.982 | 0.981 | 0.49 / 0.49 | 0.72 | 0.70 | 0.13 / 0.14 / 0.73 / 0.050 / 0.047 |
| 9 | 0.016 | 0.970 | 0.967 | 0.50 / 0.45 | 0.81 | 0.68 | 0.16 / 0.18 / 0.74 / 0.052 / 0.057 |
| 12 | 0.020 | 0.969 | 0.957 | 0.60 / 0.42 | 0.85 | 0.65 | 0.17 / 0.20 / 0.76 / 0.054 / 0.068 |
| 14 | 0.021 | 0.959 | 0.944 | 0.63 / 0.40 | 0.87 | 0.64 | 0.19 / 0.23 / 0.78 / 0.057 / 0.078 |

- **The fourth-cumulant diagonal and (2,2) slice are the gain to 94-98.5%, at every depth.** The non-gain remainder is
  1.5-6% of the energy, 12-24% of the norm.
- **The chain's errors on those slices are that remainder.** Its read error on g4 (12-20%) and K22 (13-23%) equals the
  norm of the truth's non-gain part, and 0-3% of its error is gain-shaped: the chain carries the gain at the right
  amplitude and nothing else.
- **The (3,1) slice is where the gain is weakest and the chain is worst.** 22-65% gain-shaped (22-49% in the read metric,
  which also wants about 30% less amplitude: g_read 0.016 against g_plain 0.022 at depth), and the chain is 73-94%
  wrong, about what a gain-only slice would be (sqrt(1 - 0.4) = 0.77). The production slice is lam C_off, which also
  lacks the var_c factor of the gain shape.
- **What the slices are worth** (noise-extrapolated oracles of note XXXI; base raw 2.2525e-8 / 2.1266e-8 on networks
  0 / 1): the true (3,1) slice alone gives -21.9% / -6.0%; the true (2,2) slice alone makes things worse, +24.7% / +12.7%
  (the statistics' errors compensate each other); the true diagonal, (2,2) and (3,1) slices together -41% / -33%; the
  true (2,1) and third-cumulant diagonals and g4 together -72% / -56%; all five slices -95.9% / -91.4% (replicate on
  network 0: -98%). Given true statistics the chain's readout is therefore wrong by only about 2-9% of the MSE; the
  other 90% or more is the error of the carried statistics themselves.
  The calibration of section 7 independently shrinks the (3,1) slice to 0.50-0.68 at layers 12-14.

## 4. The fresh part: tree motifs of the inter-layer graph

For jointly Gaussian d (covariance C) and y_a = sum_k (c_ak / k!) He_k(d_a), the joint cumulant of four linear reads
of y is a sum over connected multigraphs on the reads without loops, each tree (three covariance factors) with
coefficient 1 (`code/births.py`):

    star at read t:  sum_a w_t(a) c_a3 prod_(s != t) M_(a,s);   path s1-s2-s3-s4:  sum_ab w_s2(a) c_a2 M_(a,s1) C_ab w_s3(b) c_b2 M_(b,s4)

with M = C diag(c_1) W^T, the covariance between a unit of layer l and a linearised unit of layer l + 1. For the
(3,1) slice and the diagonal of the next pre-activation:

    K31_tree = 3 S1 + S2 + 6 P1 + 6 P2,   g4_tree = 4 diag(S2) + 12 rowsum(V o X),
    S1 = (W o c3 o (M^T)^2) M,  S2 = ((M^T)^3 o c3) W^T,  X = W o c2 o M^T,  V = X C,
    P1 = (V o W o c2) M,  P2 = (V o M^T)(W o c2)^T          (five n^3 products; the diagonal needs one).

These are the signed star and path counts of the bipartite graph between layers l and l + 1, with edge weights M.
**Checked** against Monte Carlo with polynomial units of known Hermite coefficients (`code/births_check*.py`,
`outputs/births_check.txt`, 4e7 samples): with c3 = 0 the exact cumulant is paths plus 4-cycles, matched to three
digits on the diagonal and the (3,1) slice; the stars are matched to three digits by their linear part in c3 (from c3
and 2 c3; the raw odd part also carries the c3^3 diagrams).

**Against truth** (`outputs/births_off{0,1}.txt`): evaluated at the true mean and covariance of layer l and fitted
jointly with the gain shape to the true slices of layer l + 1,

| step | K31 read: gain alone / gain + trees | tree amplitude | g4 read: gain / + trees | trees' share of the chain's K31 error |
|---|---|---|---|---|
| 1 -> 2 | 0.37 / 0.49 | 0.74 | 0.985 / 0.988 | 0.28 |
| 3 -> 4 | 0.48 / 0.56 | 0.75 | 0.986 / 0.989 | 0.18 |
| 6 -> 7 | 0.45 / 0.52 | 0.82 | 0.976 / 0.980 | 0.10 |
| 9 -> 10 | 0.43 / 0.48 | 0.82 | 0.969 / 0.974 | 0.05 |
| 13 -> 14 | 0.40 / 0.42 | 0.68 | 0.959 / 0.962 | 0.02 |

(network 0; network 1 agrees). The trees are real (amplitudes 0.7-0.8, correlation with the gain residual +0.36 at
shallow depth) and small where the error is: the one-step content is 2-12% of the (3,1) slice's energy and falls with
depth while the non-gain part grows. The deviation the chain misses is a transported, accumulated quantity, as the
cocycle reading says; a one-step correction cannot reach it. At five n^3 products per layer (about 0.07 B before
Strassen) the trees are also not worth their cost as a stand-alone repair.

## 5. The accumulated non-gain content is collective

A latent model, z - mu = U t + (an independent rest) with U the top-K eigenvectors of C, makes every off-diagonal slice
linear in U_j on its single index (K31_ij = T4[U_i, U_i, U_i, U_j], D21_ij = T3[U_i, U_i, U_j]), so the plain slice's
column space lies in span(U), and the read-weighted slice's in span(diag(Phi) U). `code/latent_span.py`
(`outputs/latent_span_off{0,1}.txt`) measures the share, for the truth and for its gain residual, against K / n for a
random subspace (0.016, 0.062, 0.25):

| layer | K31 residual, plain: K = 16 / 64 / 256 | read | D21 residual, plain | read |
|---|---|---|---|---|
| 4 | 0.23 / 0.54 / 0.88 | 0.24 / 0.57 / 0.91 | 0.23 / 0.57 / 0.90 | 0.24 / 0.59 / 0.93 |
| 8 | 0.41 / 0.71 / 0.92 | 0.43 / 0.73 / 0.95 | 0.36 / 0.70 / 0.93 | 0.37 / 0.72 / 0.95 |
| 12 | 0.56 / 0.81 / 0.95 | 0.58 / 0.82 / 0.97 | 0.45 / 0.78 / 0.96 | 0.46 / 0.80 / 0.97 |
| 15 | 0.61 / 0.84 / 0.97 | 0.62 / 0.85 / 0.98 | 0.56 / 0.83 / 0.97 | 0.57 / 0.84 / 0.98 |

(network 0; network 1 within 0.02.) The concentration grows with depth, the signature of an unstable subspace
collecting the transported content. (An earlier version of this measurement mixed the read weights with the ungated
basis and reported 25-48% at K = 64; that combination rotates the column space out of span(U) and is not a test of
anything.) The column test is necessary, not sufficient: the latent model also fixes the dependence on the repeated
index (cubic in U_i for K31). Section 6 tests the whole model.

## 6. The latent split by Monte Carlo

`code/mclatent.py` draws 8e6 inputs per network (8 chunks of 1e6, `default_rng([4711, net, chunk])`) and, at layers 3, 6,
9, 11, 13, 15, splits each sample's fluctuation d = z_l - mu_l into its collective part c = U U^T d (U = top-64
eigenvectors of the true covariance) and the residual r = d - c, accumulating the raw-moment sums of d, c and r.
`code/latent_an.py` forms their slices. With one gain amplitude g fitted per slice on d, the non-gain part
R_d = K31_d - g S_d splits exactly into R_c, R_r and a cross remainder (c and r are uncorrelated, so the gain shapes
add: S_d = S_c + S_r + 3(var_c C_r + var_r C_c)). A latent model (z - mu = U t + an independent Gaussian rest) predicts
R_d = R_c. Network 0 (`outputs/latent_split_off{0,1}.txt`; network 1 within 0.05 on every explained share):

| layer | non-gain share of K31 energy (MC noise) | R_d explained by collective T4 | by residual | cross share of R_d |
|---|---|---|---|---|
| 3 | 0.62 (0.26) | +0.06 | +0.11 | 0.83 |
| 6 | 0.53 (0.08) | +0.19 | +0.03 | 0.78 |
| 9 | 0.51 (0.03) | +0.31 | +0.01 | 0.68 |
| 11 | 0.47 (0.02) | +0.40 | +0.01 | 0.59 |
| 13 | 0.38 (0.01) | +0.50 | +0.003 | 0.50 |
| 15 | 0.35 (0.006) | +0.55 | +0.002 | 0.45 |

- **The collective model becomes right with depth and never completes.** The collective fourth cumulant
  T4[U_i, U_i, U_i, U_j] explains 6% of the non-gain (3,1) content at layer 3 and 55% at layer 15; the residual's own
  fourth cumulant is negligible at depth. The remaining 45-83% is the CROSS part, the dependence between the collective
  state and the idiosyncratic residual (conditional heteroscedasticity). So units are not conditionally independent given
  a 64-dimensional latent state, and the dependence enters at fourth order.
- **The cross part is collective in its column index and unit specific in its rows.** `code/latent_mgain.py` fits the
  model r | t ~ N(0, diag(s^2)(1 + h(t))) with one modulation h shared by all units, which predicts
  X = diag(s^2) U M' U^T for the cross part X with one KxK matrix M'. It explains 4-10% of X at K = 64 (held out on the
  other half of the samples: 0.4-8.6%), whereas the best row-wise collective matrix (X U U^T) captures 48% (layer 3) to
  66% (layer 15). The same M' predicts the cross part of the fourth-cumulant diagonal with correlation +0.60-0.74. So the
  slope of unit i's residual variance in the collective coordinates is unit specific: it is the t-derivative of
  Cov(z'_i | t) = sum_ab W_ia W_ib Cov(y_a, y_b | t), a third cumulant of (y_a, y_b, t_k), of which the chain carries only
  the a = b part (D21 contracted with the collective eigenvectors).
- **Third cumulant and (2,1) slice.** The collective part alone explains the third-cumulant diagonal to 0.55 (layer 3)
  - 0.94 (layer 15) and the (2,1) slice to 0.45-0.91 (plain), so for them a collective model is nearly complete at depth.
- **What this rules out.** A collective-mode cumulant propagation with a conditionally independent rest closes the
  fourth-order sector only to the extent of the first column of the table. Together with the precision law (about 1e-5
  relative energy needed) it is not a replacement state at any affordable K; its content is the diagnosis that the
  missing fourth-order information is the t-derivative of the conditional covariance.

### 6.1 Pre-registered prediction for the projected-slice oracles

The oracle runs of note XXXI replace a statistic by truth at every layer. `code/k31proj.py` builds the inputs in which
the true K31 is replaced by (colK) its projection K31 U_K U_K^T on the top-K covariance eigenvectors of its single
(column) index, K = 16, 64, and by (gain) its scale-mixture fit 3 g var_i C_ij. They were run (26 chain runs, networks 0
and 1, with the two Monte Carlo halves for the noise extrapolation) before this section was written and not looked at;
the prediction below is derived from the latent theory and the measurements above. Write phi_X = (MSE_base - MSE_X) /
(MSE_base - MSE_trueK31) for the fraction of the true-K31 oracle gain that replacement X retains. If the chain's
sensitivity to an error in the slice were isotropic in the read metric, then

    phi_X = 1 - E_out(X) / E_chain,

with E_out(X) the read energy of the true slice that X lacks (its relative energy outside the span for colK: 1 - f_K,
f_K from `outputs/latent_span_off*.txt`) and E_chain = 0.56-0.72 the relative error energy of the chain's own slice
(read error 0.75-0.85). With f_16 = 0.5-0.8, f_64 = 0.75-0.92, f_256 = 0.93-0.97 over the layers that carry the
sensitivity (10-15), this gives

| replacement | predicted phi (central; range) |
|---|---|
| colK, K = 256 (not run) | 0.9; 0.8-1.0 |
| colK, K = 64 | 0.75; 0.5-0.95 |
| colK, K = 16 | 0.35; 0.0-0.7 |
| gain fit | 0.2; 0.0-0.5 |

The robust prediction is the ordering phi_gain < phi_16 < phi_64 <= 1. If the readout weights the collective directions
more than isotropically (the sums over j are coherent for collective columns and incoherent otherwise), every value
moves up; if phi_64 < 0.5 the idiosyncratic out-of-span part of K31 matters much more than its energy share, and the
column-collective picture is refuted for K31. The gain fit retains little because the chain's own slice (lam C_off) is
already a gain shape with a fitted amplitude: replacing it by a better-fitted gain shape removes only the difference of
two similar errors. Network 1 has a true-K31 gain of only 6%, so its phi values are noise dominated and are reported
but not used for the decision.

### 6.2 Result: the prediction fails in part

Noise-free values (2 MSE(full) - mean MSE(halves); `outputs/k31_projection_oracles.txt`), K31 replaced at every layer:

| network (base raw) | true K31 | colK, K = 64 | colK, K = 16 | gain fit |
|---|---|---|---|---|
| 0 (2.2525e-8) | -22.5% | -13.9% | **+8.8%** | +2.0% |
| 1 (2.1266e-8) | -7.2% | -3.2% | -3.5% | +5.8% |

Retained fractions phi (network 0, where the true-K31 gain is large): colK64 +0.62, colK16 **-0.39**, gain **-0.09**;
network 1 (true-K31 gain only 7%, so noise dominated): +0.44, +0.48, -0.80.

Against section 6.1: colK64 falls inside the registered range (0.62 in 0.5-0.95, below the central 0.75); colK16 falls
outside it (-0.39 against 0.0-0.7) and the gain fit marginally outside (-0.09 against 0.0-0.5); the registered
ordering phi_gain < phi_16 < phi_64 FAILS on network 0 (observed phi_16 < phi_gain < phi_64). The premise of the
prediction, an isotropic response of the output to errors of the slice in the read-metric energy and no interaction with
the other statistics, is refuted.

What the failure says (derived from these runs and the oracle table of section 3, not fitted):

1. **Truncating K31 to its top 16 collective modes is worse than the chain's own slice** (+8.8%), although the dropped
   part carries only 20-70% of the read energy (section 5) against the chain's own error energy of 56-72%. The output's
   sensitivity to the dropped part is therefore several times larger than its energy share: the modes 17-64 of the
   column index are read.
2. **The statistics' errors compensate each other, so one slice cannot be improved in isolation.** The same table
   gives +24.7% / +12.7% for a TRUE K22 alone, and here a truth-fitted gain-shaped K31 (3 g var_i C_ij with g fitted on
   the slice) is worse than the chain's slot, which is lam C_off without the var factor at a fitted amplitude. The chain
   sits at an error-cancelling configuration; the all-true oracle (-91% to -96%) is a statement about the whole bundle.
3. **Minimal collective rank for K31.** Whatever supplies the (3,1) slice must carry its column index at K >= 64
   (colK64 retains 62% of the gain on network 0); K = 16 does harm.

Decision: a component acting on K31 (a latent model, a response correction, tree births) must be judged inside a
coherent bundle with g4, K22 and D21 against the all-true direction, never against the K31 oracle alone. The
isotropic energy bookkeeping used in 6.1 is withdrawn; sensitivity-weighted (adjoint) bookkeeping is required.

## 7. The calibration holds out of sample

The 103-direction output-metric projection of note ray-compiler section 8 (fitted on networks 0-49: train -6.04%,
predicted held-out -3.07% +- 0.38) was run free on networks 50-99 (`V47_CAL` with the fitted coefficients, production
configuration, paired; `outputs/calval_50nets.txt`): mean raw 2.2692e-8 -> 2.1838e-8, **-3.76%** (per network -3.73%
+- 0.37), better on 46/50, at unchanged cost. The free run does slightly better than the linear prediction.

## 8. What follows

Sections 6 and 6.2 left one question open: can the omitted fourth-order content be carried as a collective state? Note XXXVIII (`notes/pair-flex`) answers it in two parts.
- **It has to be carried.** Under conditional independence given the collective coordinates, the pair slices determine the omitted (2,1,1) and (1,1,1,1) classes only modulo an exact kernel of dimension C(K+3, 4) + q(q-1)/2 (its Theorem 2). Any closure computed from pair statistics is therefore a convention.
- **Carrying it does not pay at affordable K.**
  - Their transport is reproduced by their own collective projections with nothing fitted (coefficient 1 within 9% from layer 8 on).
  - Convergence in K is slow: 0.59-0.79 of the (3,1) slice at K = 64.
  - At K <= 32 the content equals the scale mode, which section 4's core and note XXXVI's KD already tested in the chain.

The collective-state extension is closed for cost there.
