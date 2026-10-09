# Chaos grading: the network's Wiener chaos as the chain's bookkeeping, and the Schur hub

Working note XLI. This round asked what the algorithmic papers supplied with it say about the chain when their essences
are carried over, not copied. The papers were:
- the multiway-cut rounding mixtures;
- the subspace designs;
- Hamiltonian and quantum cut sparsification;
- the token-graph and Kikuchi-graph spectra and the Quantum Max-Cut algorithms;
- the one-dimensional Gibbs-state preparation;
- spectrum estimation;
- Hamiltonian learning at any temperature;
- the folded Reed-Solomon and multiplicity codes;
- the shallow-light Steiner trees;
- the sigma-inverse mean curvature flow;
- f-connectedness.

The essence that turns out to bite is shared by the free-fermion and Bogoliubov papers:
- quadratic fluctuations around a mean field form a solvable sector;
- every observable of that sector is a walk sum, a matrix function of one operator.

In the chain that sector is the second Wiener chaos of the input. This note proposes the chaos grading as a replacement
for one inherited component, the bookkeeping of cumulants by index class (diagonal, slices, bulk). It derives from the
grading a construction, the Schur hub, that moves one omitted kappa4 class from n^3 per source-layer to one solve per
layer. The predictions in section 4 are stated before their runs.

## 1. The grading

**The expansion.** Every pre-activation is a function of the Gaussian input x, so it has a Wiener chaos expansion
z_(l,i) = mu + I_1(L_i) + I_2(H_i / 2) + I_3(K_i) + ... with symmetric kernels on the input space. The pieces of the
chain read as follows in this grading:
- **First chaos.** L_l = E[grad_x z_l]. It is transported exactly by the mean gate: L_(l+1) = W_(l+1) d(Phi_l) L_l.
  This is Stein's lemma, so no pair correction ever enters the first chaos.
- **Second chaos.** The kappa_3 sources are the second chaos. A birth at neuron m of layer b contributes the kernel
  w2_m l_m l_m^T, with l_m the first-chaos vector of z_(b,m) and w2 = E f''. It is carried to layer l by the source's
  P leg (the mean Jacobian): H_(l,i) = sum_b sum_m P_im w2_m l_m l_m^T.
  - The chain's arm is the first-chaos cross-covariance, A_im + P_im w1_m var_m = L_(l,i) . l_m. The diagonal term
    w1 var is carried by the M block's separable part.
  - The M block's exact Gaussian coefficient e = 2 E[f] (1 - Phi) equals the second-chaos value 2 w2 w1 var at
    alpha = 0. Elsewhere it also holds the single site's higher chaos.
  - So the star and the M block together are, to that accuracy, the second-chaos star with the full arm
    At = A + P d(w1 var).
- **Third chaos.** It comes from two places:
  - the He_3 births (c3 = E f''', the level-4 star of note XL);
  - the chaos product of a vertex's slope fluctuation with the incoming second chaos, I_1(L_a) I_2(H_a) -> I_3. This is
    the kappa3 x C class of note XL. Its first-chaos remainder 2 w2 H_a L_a is the folded feedback.

**Why the grading is the consistent one.**
- The chaos components are orthogonal, so truncating at a chaos order is a projection.
- The index-class grading mixes chaos orders: every slice and every bulk class collects pieces of several. Truncating an
  index class is therefore not a projection, and note XL measured the consequence: the level-4 star alone (a third-chaos
  piece with the diagonal of its arm missing) was adverse, +2.7%.
- In the chaos grading the C^3 order of the kappa4 bulk splits into two consistent sub-systems:
  - the second chaos squared (the path class P4, including pairs of births at different layers);
  - the third chaos (the star with the full arm, plus the product term).

**The second-chaos sector is solvable.** Truncated at the second chaos, a pre-activation is a generalized chi-square.
Its cumulant generating function is the Gaussian integral of a quadratic exponent:

    K_i(t) = -(1/2) log det(I - t H_i) - (t/2) tr H_i + (t^2/2) L_i^T (I - t H_i)^(-1) L_i,
    kappa_k(z_i) = k! [ tr H_i^k / (2k) + L_i^T H_i^(k-2) L_i / 2 ],

so that
- kappa_3 = 3 L^T H L + tr H^3 (the chain's D3, plus the triangle);
- kappa_4 = 12 |H L|^2 + 3 tr H^4;
- every cumulant is a walk sum on one operator: open walks L^T H^k L (paths) and closed walks tr H^k (cycles).

This is the token-walk picture of note XL made exact, and it matches the papers.
- It is the free-fermion statement of the token-graph papers: the many-token sector of a quadratic Hamiltonian is fixed
  by one-body data through a determinant.
- It is the Bogoliubov statement of the dense Quantum Max-Cut algorithm (Jiang-Kwon-Wright): fluctuations around the
  mean field, which here is the first chaos, are a quadratic sector whose observables are matrix functions of one flat
  interaction.

The joint-gate correction (V2) is a closed walk: two births at the current layer and one carried source make tr H^3,
the triangle. Its contraction graph is K_4, so it stays a wall at n^4, as note XL found.

## 2. The Schur hub

**The path class is a regression.** The path part of the kappa4 diagonal is 12 |v_i|^2, with
v_i = H_i L_i = sum_s sum_m X_(s,im) l_m and X_s = P_s o At_s o w2_s. Directly that is two costs:
- 12 sum_(s,s') X_s G_(ss') X_(s')^T, where G is the Gram matrix of the birth vectors;
- one n^3 product per pair of sources, or per source-layer for the within-source part (note XL section 8).

Three identities turn it into a regression on quantities the chain already holds:
- Y_ij = L_j . v_i = sum_s [X_s At_s^T]_ij is the first of the two D21 hub products, once the hub reads the full arm At
  (same cost; the hub's two right factors become At and P).
- While L_l is invertible, G_(ss') = At_s^T C_l^(-1) At_(s'): the birth Gram is the current layer's view of the births,
  by the Schur identity Cov(z_b) = Cov(z_b, z_l) Cov(z_l)^(-1) Cov(z_l, z_b) of the first chaos.
- Hence

      |v_i|^2 = [ Y C_l^(-1) Y^T ]_ii    for every pair of sources at once.

In probabilistic terms this is the variance of z_i^2 explained by linear regression on the layer's field: the
positive-semidefinite completion of [[C, D21_i], [D21_i^T, Var z_i^2]], applied to the second-chaos part of D21.

**Where the papers enter.**
- **Shallow-light Steiner trees.** The Schur hub is a Steiner point. It replaces the S^2 direct source-pair connections
  by one shared object, C_l^(-1), through which every pair connects. Its depth, the amplification of errors through
  the hub, is the condition number of C_l. That is large: a square Wishart matrix has smallest eigenvalue about
  1/n^2 of its mean.
  - An incoherent relative error delta in Y is amplified by about n delta^2 relative to the signal, which is order one
    at delta = 5%. So the hub must be regularized, (C + eps)^(-1).
  - This drops the part of v_i in the near-null directions of L_l: about (2/pi) sqrt(eps / var) of its energy for a
    Marchenko-Pastur spectrum.
  - The trade between eps (the light weight) and amplification (the depth) is the shallow-light trade.
- **f-connectedness and strong minimality.** The identity needs L_l to see every input direction. No small separator
  cuts the births off from the current layer, and the conditioning measures how close the network comes to one.
- **Spectrum estimation (bucketing plus moment matching).** The kappa4 diagonal is the second local moment
  L^T H^2 L of the operator whose first moment L^T H L the chain reads as D3. The hub estimates that second moment from
  carried readouts. It does not reconstruct H, much as power sums are estimated without the eigenvalues.

## 3. What the grading predicts about the chain's kappa4 error

The chain's kappa4 diagonal is a closure:
- a transported diagonal content (W o W) g;
- a fitted lambda C_off term;
- a rank-4 quenched pair class.

Its per-neuron error is 30% rms at depth, and the g4 oracle is worth 15-30% of the MSE (note XXXI). The path class of
the second chaos is a law term: per-neuron, quenched information the chain does not carry. The flex theorem (note
XXXVIII) does not forbid it, because it is computed from the sources' chaos structure (which neuron each piece of H
belongs to), not from the pair state. The same kappa_3 tensor with a different H gives a different kappa_4. That
non-uniqueness is the flex kernel, and the sources resolve it.

## 4. Predictions for the diagnostic (stated before the run)

`workbench/k3work/chaos2_diag.py`, networks 0 and 1. The setup:
- the adopted fold system with every source dense (`V21_NO_CONFINE=1`, so all legs are available);
- legs dumped at layers 3, 5, ..., 15;
- Monte Carlo truth from 1.6e7 inputs, with split halves for the noise.

| | quantity | prediction |
|---|---|---|
| D1 | the second-chaos star with the full arm, 2 Y + Q, against the chain's D21 | relative error <= 0.3 (the identification of section 1); chain arms against the first-chaos arms L_l L_b^T: median relative error <= 0.2 |
| D2 | the Schur hub against the exact all-pairs path class (v_i from the first-chaos maps) | at eps = 1e-2 mean(var): correlation >= 0.9, mean within 15%; at eps = 0: unstable (correlation < 0.9 or mean off by more than 50%) |
| D3 | the exact path class, demeaned, against the chain's per-neuron kappa4 residual (truth - g4row) | at layers >= 7 on both networks: at least 20% of the residual's signal variance removed, slope in [0.3, 1.5] |
| D4 | cross-source pairs | mean of the within-source part / all pairs <= 0.8 at layers >= 9 |
| D5 | the third-chaos star with the full arm against the V55 star (arm A) | removes more of the residual's signal |

**Decision.**
- If D3 holds, the hub goes into the chain:
  - `V56_P4`: the D21 hub reads the full arm, and g4row gains the regularized path class with its mean taken out of
    the closure's lambda;
  - it is tested in the cold harness on networks 0-15 and then in the scored regime with refitted counterterms, in
    adjusted MSE.
- Cost: one regularized solve per layer, about 1-2 units each, against a g4 oracle of 15-30% of the MSE.
- If D3 fails, the second chaos does not hold the chain's kappa4 error, and the third chaos (D5) is the next place to
  look.

## 5. The diagnostic, measured (`outputs/diag_v1_off{0,1}.txt`, `outputs/diag_v2_off{0,1}.txt`)

The setup, as stated in section 4: the fold system with every source dense (raw 2.262e-8 and 2.047e-8 on networks 0
and 1), and Monte Carlo truth from 1.6e7 inputs. The chain's kappa4-diagonal residual (truth - g4row) is about 30% of
the truth's rms at depth, and it is almost all signal: the split-half noise is 5e-5 to 2e-4 against a residual signal
of 7e-4 to 1.6e-3.

| | prediction | measured | |
|---|---|---|---|
| D1 | 2Y + Q within 0.3 of the chain's D21; chain arms within 0.2 (median) of the first-chaos arms L_l L_b^T | D21: 0.22 (layer 3) to 0.52 (layer 14). Arms: 0.30 (layer 3) rising to 1.38 (layer 15) | **fails** |
| D2 | the Schur hub (eps = 1e-2) correlates >= 0.9 with the exact all-pairs class, mean within 15%; unstable at eps = 0 | correlation 0.95, but means 6.6-9x apart. Stable at eps = 0 (same correlation, mean 10% higher) | **fails as stated** |
| D3 | the exact path class alone removes >= 20% of the residual's signal at layers >= 7, slope in [0.3, 1.5] | 6-41% (network 0), 10-35% (network 1); slopes 0.36-5.6 | **fails** |
| D4 | within-source / all pairs <= 0.8 at layers >= 9 | 0.89-1.14 in the chain's metric | **fails** |
| D5 | the full-arm star removes more than the V55 star | alone 0-13% against 0-7%. Jointly it is the essential partner (below) | holds, weakly alone |

**What failed, and why.**
- The chain's arms are not first-chaos cross-covariances. They are the full ones: the transport moves the whole
  covariance by the mean gate.
- The higher-chaos share of the cross-covariance between a birth layer and the current layer grows with their
  separation, until it exceeds the first-chaos part. The "exact" path class with first-chaos arms is therefore
  4-8x too small at depth, and its slopes (up to 5.6) say so.
- The Schur hub uses the chain's own full covariance and arms, a self-consistent metric. It does not suffer from this.
- It is also well conditioned in practice. The hub product Y lies in the covariance's top (outlier) subspace, where
  the arms are born, so the solve never sees the near-null directions that section 2 worried about.

**What holds: the two consistent pieces of the grading carry most of the chain's kappa4 error.** The residual was
regressed on the demeaned candidates; "removed" is the share of the residual's signal variance.
- "fit h0 -> h1" means fitted on one Monte Carlo half and scored on the other.
- The slopes are the fitted weights of the two pieces.

| layer | network 0: removed (fit h0 -> h1) | slopes: path, star | network 1: removed (fit h0 -> h1) | slopes: path, star |
|---|---|---|---|---|
| 3 | 0.64 (0.62) | 1.07, 0.88 | 0.66 (0.69) | 1.11, 0.86 |
| 5 | 0.70 (0.72) | 1.13, 0.77 | 0.67 (0.67) | 1.11, 0.74 |
| 7 | 0.70 (0.70) | 1.12, 0.74 | 0.68 (0.69) | 1.09, 0.68 |
| 9 | 0.69 (0.69) | 1.12, 0.69 | 0.67 (0.68) | 1.12, 0.66 |
| 11 | 0.72 (0.73) | 1.15, 0.69 | 0.74 (0.73) | 1.20, 0.68 |
| 13 | 0.81 (0.83) | 1.16, 0.68 | 0.71 (0.69) | 1.15, 0.66 |
| 14 | 0.80 (0.80) | 1.16, 0.69 | 0.75 (0.74) | 1.19, 0.71 |

The two pieces are the Schur-hub path class 12 diag(Y (C + eps)^-1 Y^T), with the hub's own exact-e left factor,
and the third-chaos star with the full arm, 4 sum c3 P At^3. Together they remove 64-81% of the error's signal at
every layer on both networks.

**The weights are nearly constant.**
- The path class sits at 1.07-1.20 at every layer of both networks. That is the derived weight plus the 4-cycle: the
  4-cycle 3 tr H^4, computed for eight neurons, is 24-30% of the path term and correlated with it.
- The star sits at 0.66-0.88. Its partner in the third chaos, the product term of the kappa3 x C class, is left out.
- Each piece alone does much less: the path class 15-51%, the star 0-13%. They are complementary halves of the error,
  as the grading says the second and third chaos are.

**The control confirms the flex kernel.**
- Replacing the physical hub product Y by the symmetrized slice (3 diag(D21 (C + eps)^-1 D21^T), or its diagonal-metric
  form) removes 0-6% at layers 5-14, and 0-14% jointly with the star.
- So the same kappa_3 numbers carry nothing about the kappa_4 diagonal. Only the assignment of each second-chaos piece
  to the neuron that carries it (the sources' structure, section 3) does.

**Also measured.**
- **The star.** The star's four arm-power terms (A^3 P, A^2 P^2 t, A P^3 t^2, P^4 t^3) have unstable separate
  weights, while the full-arm star is stable. The V55 star (A^3 P alone) is an inconsistent truncation, as section 1
  said.
- **The Euler identity mu_i = tr H_i** of homogeneity (E[F] = E[Delta F]), with the sources' H, correlates 0.73-0.95
  with the chain's mean.
- **The diagonal-metric hub** (12 sum_j Y_ij^2 / var_j, n^2 given Y) alone removes about what the Schur hub removes
  (16-50%). But it is 30-60x too large, because Y lies in the covariance's top subspace, where C^-1 is far below
  diag(C)^-1. Jointly with the star it reaches 36-66%, against 64-81% for the Schur metric.

**Value and cost.**
- The g4 oracle is worth 15-30% of the MSE (note XXXI). Removing about 70% of the oracle's error variance is therefore
  worth roughly 10-20% of raw MSE, if the effect is near linear.
- The price:
  - the hub product Y, which is the D21 hub split so that its first half reads the full arm (no new n^3 work);
  - one regularized solve per layer (1.33 units at flopscope's pricing, 2n^3/3 + 2n^3 for an n x n right-hand side);
  - the star's n^2 per source-layer.
- The trimmed last layer forms no hub, so its kappa4 diagonal needs separate treatment.
- Section 6 states the chain test.

## 6. The chain test (stated before its runs)

**`V56_P4` in `workbench/k3work/est_v29.py`.**
- **The hub.** The young D21 hub runs as its two contractions in two families on half the slots each. The full arm
  At = A + P d(t), with t = w1 var recorded at birth, is formed in place in the A slots and restored after:

      sum_s LA_s At_s^T  +  sum_s (LP_s - LA_s o t_s) P_s^T  =  D21 hub (unchanged up to rounding).

  The first family is 2 Y. Each family runs in the buffers the fused one owns, at the same products.
- **The path class.** 12 diag(Y (C + 0.01 mean(var))^-1 Y^T), one solve per layer (1.33 units).
- **The star.** 4 sum c3 P At^3 over the dense slots (n^2 per source-layer).
- **How they enter.** Both are added to the closure's kappa4 diagonal with their active-neuron means removed (weights
  A, B), after the closure's own (2,2) slice and mixture gains have read it, as V55 did.
- **Modes.**
  - Mode 1: the correction enters everything.
  - Mode 2: the closure's transported memory (K4v -> g_prev, k4q) is kept free of it, through its leading image
    w1^4 d. The path class is recomputed in full from the sources at every layer, so its carried part must not also
    ride the transported diagonal.
- **Coverage.** Only dense (young) sources are covered, so the screen runs with every source dense (`V21_NO_CONFINE=1`)
  for all variants and for the base.

**Predictions.** Cold harness, networks 0-15. The base is the fold system with `V21_NO_CONFINE=1`, rerun in the same
batch, and every variant is paired against it.

| tag | switches | prediction |
|---|---|---|
| m2fit | `V56_P4=2 V56_A=1.15 V56_B=0.7` (the diagnostic's weights) | **P22:** raw -8% to -20%, better on >= 13/16 |
| m2der | `V56_P4=2 V56_A=1 V56_B=1` (derived weights) | P23: worse than m2fit by 0-6 points |
| m1fit | `V56_P4=1 V56_A=1.15 V56_B=0.7` | P24: worse than m2fit by >= 2 points (the carried path content counted twice) |
| m2path | `V56_P4=2 V56_A=1.15 V56_B=0` (path class only) | P25: raw -3% to -12% |
| m2star | `V56_P4=2 V56_A=0 V56_B=0.7` (star only) | P26: within +-3% (the star alone carries 0-13% of the error signal) |

FLOPs: + one solve per layer, about 18-19 units over 14 layers. The hub split adds no products.

**Decision.**
- If P22 holds, the production version follows:
  - the old tier's half of Y (the factor-space family split the same way, plus one extra r-lift) and its star (from
    the formed legs);
  - a cheaper solve (Galerkin in the covariance's top subspace, where Y lies, section 5).
- Then the scored regime with refitted counterterms, in adjusted MSE (note XXXIX protocol).

**Amendment after a smoke test on network 0, before the screen** (`V56_LMIN`, `V56_LAST`; the pre-registered
variants run as specified).
- m2fit as specified gives +3.5% raw on network 0 (2.342e-8 against 2.262e-8).
- Per layer, though, its MSE falls sharply from layer 3 to layer 14: -18% at layer 3, -27% at layer 11, -14% at
  layer 14.
- The loss comes from two ends:
  - **Layers 1-2.** They get worse (+36%, +15%). The diagnostic never examined them, and at layer 1 the closure's
    lambda term was fitted to the bulk born at layer 0, the same content.
  - **The trimmed last layer.** It forms no hub, so only the star acts there. The diagnostic's star alone has a
    negative slope at layer 15 (it is the path class's partner, not a correction by itself). The last step then adds
    7.2e-9 against the base's 3.8e-9.
- Cost on network 0: +6.9% FLOPs (+18.4 units, the 14 solves).

| tag | switches (all with V21_NO_CONFINE=1) | prediction |
|---|---|---|
| a1 | m2fit + `V56_LMIN=3 V56_LAST=0` (layers 3-14) | A1: raw -6% to -15%, better on >= 13/16 |
| a2 | m2fit + `V56_LMIN=3 V56_LAST=2` (layers 3-15, the last layer's hub family and covariance formed for it) | A2: better than a1 by 2-6 points |
| a3 | m2der + `V56_LMIN=3 V56_LAST=2` | A3: within +-4 points of a2 |
| a4 | m2path + `V56_LMIN=3 V56_LAST=2` | A4: raw -2% to -10% |

P22 as specified is expected to fail.

## 7. Why the dressed hub is the right object, and what each paper contributed

**The path class is a regression, and the hub dresses its lines.** Write C_i = (L_i . x) q_i for the cross term
between a neuron's first and second chaos.
- Its first-chaos projection is 2 v_i . x, so the chaos-2 path class is 12 |v_i|^2 = 3 Var(J_1 C_i).
- Its covariance with the field is Cov(C_i, z_j) = 2 Y_ij. The hub therefore computes three times the variance of
  C_i's linear projection on the whole field z_l, through the full covariance:

      P4_i = 12 [Y C^-1 Y^T]_ii = 3 Var( proj_(span z_l) C_i ).

- The "exact" version uses the first chaos only: bare propagators on the internal line (the birth Gram L_b L_b'^T)
  and on the external lines (the arms L_l L_b^T).
- The hub uses the full covariance and the chain's full arms. This is the skeleton expansion of the same diagram,
  every line a full propagator, with the self-energy insertions resummed.
- At depth the higher-chaos share of the cross-covariances exceeds the first-chaos share (section 5, D1). The bare
  diagram is then 4-8x too small and the dressed one carries the error at weight about 1.15.
- What is left over is accounted for by the two omitted pieces:
  - the 4-cycle (24-30% of the path term in the bare version, correlated with it) lifts the path weight from 1;
  - the third chaos's product term (the kappa3 x C class) lowers the star's weight from 1 to about 0.7.

**The papers, essence by essence.**
- **Free fermions and the token graphs (APS, BBKL, AGM).** The second chaos is the quadratic sector. Its observables
  are walk sums on one operator, so the cost classes of note XL become walk shapes:
  - open walks (paths) through the current layer's hub are cheap;
  - closed walks (cycles: the joint-gate triangle tr H^3, the 4-cycle) keep their K4 contraction graph and stay at
    n^4.
- **The dense bipartite Quantum Max-Cut algorithm (JKW).** Fluctuations around the mean field form a quadratic
  (Bogoliubov) sector driven by the flat part of the interaction; the additive part acts trivially. Here the mean field
  is the first chaos, the quadratic sector the second chaos.
  - The measured value sits exactly in the flat, quenched assignment of second-chaos pieces to neurons: the
    symmetrized-D21 control carries nothing.
  - This is note XL's finding that value lies in quenched content, now located inside one diagram.
- **Optimal spectrum estimation (bucketing and moment matching).** The kappa4 diagonal is the second local moment
  L^T H^2 L of the operator whose first moment L^T H L the chain reads as D3. The hub estimates it from carried rows,
  without the operator.
- **The one-dimensional Gibbs states (positive splits).** The Schur complement is the canonical positive split:
  Cov(z_b) = (the part the current layer explains) + (a positive remainder). The hub keeps the explained part; the
  regularization discards the remainder.
- **Shallow-light Steiner trees.** The hub is a Steiner point that connects the S^2 source pairs through one shared
  object. Its depth, the conditioning, turned out shallow, because Y lies in the covariance's outlier subspace.
- **Subspace designs, f-connectedness, strong minimality.** The identity needs every birth to be visible from the
  current layer, with no small separator; the measured stability at eps = 0 says the network supplies that.
- **Hamiltonian learning at any temperature.** The hub is an identification of a hidden quantity, |H_i L_i|^2, from
  local data (a D21-type row and the covariance) by an exact regression identity: learning from marginals. The flat
  approximations explain the earlier tadpole failure (truncated Edgeworth vertices are not flat outside the band); that
  lesson was not used here.
- **Hamiltonian and cut sparsification.** The useful form was the restricted one: preserve one functional family (the
  path class of the diagonal) instead of the operator, through one hub instead of S^2 terms.
- **The multiway-cut rounding mixtures.** The essence that transfers is the one notes XXXIX and XL had measured:
  mixtures of law terms help, mixtures of conventions do not. The two chaos pieces are law terms, and their fitted
  weights are near derived values (1 plus the 4-cycle share; 1 minus the product share), not free.
- **Multiplicity codes.** The chaos kernels are expected derivatives (Stein): a multiplicity representation of each
  neuron at the Gaussian point. Homogeneity relates its 0th and 2nd "symbols": mu_i = tr H_i (E F = E Delta F), measured
  at correlation 0.73-0.95 with the sources' H.
- **The sigma-inverse mean curvature flow.** It supplied only the pattern of a sharp inequality with a rigid case:
  Var(z_i^2) >= D21_i^T C^-1 D21_i, with equality at the chaos-2 model. The sharp form that worked is the
  source-resolved one, not the D21 bound.

## 8. The cold screen, measured (`outputs/screen56_cold_16nets/`)

Networks 0-15, every variant paired against the fold system with every source dense (`V21_NO_CONFINE=1`), rerun in
the same batch (mean raw 2.118e-8, 0.2585 B).

| tag | switches | raw against the base | better on | FLOPs | prediction |
|---|---|---|---|---|---|
| m2fit | P4=2, weights (1.15, 0.7), all layers, star alone at the last | +6.46% +- 3.28 | 5/16 | +7.0% | P22, -8% to -20%: **fails** |
| m2der | P4=2, weights (1, 1), same | +20.58% +- 3.47 | 1/16 | | P23: **fails** |
| m1fit | P4=1, weights (1.15, 0.7) | +6.54% +- 3.28 | 5/16 | | P24, worse than m2fit by >= 2: **fails** (equal: the memory channel is negligible) |
| m2path | the path class alone at 1.15 | +95.0% +- 8.2 | 0/16 | | P25, -3% to -12%: **fails** |
| m2star | the star alone at 0.7 | +79.4% +- 3.9 | 0/16 | | P26, within +-3%: **fails** |
| a1 | weights (1.15, 0.7), layers 3-14 | -11.28% +- 2.27 | 14/16 | +6.0% | A1, -6% to -15%: holds |
| a2 | weights (1.15, 0.7), layers 3-15 (last layer full) | -12.15% +- 2.41 | 14/16 | +10.1% | A2, better than a1 by 2-6: **fails** (0.9) |
| **a3** | **weights (1, 1), layers 3-15 (last layer full)** | **-20.87% +- 1.12** | **16/16** | **+10.1%** | A3, within +-4 of a2: **fails** (8.7 points better) |
| a4 | the path class alone, layers 3-15 | +130.8% +- 10.2 | 0/16 | | A4: **fails** |

Per-layer MSE of a3 against the base: -16% at layer 3, -23% to -26% at layers 4-8, and -21% to -22% from layer 9 to
the output.

**What the screen says.**
- **The two chaos pieces are one correction.** Each alone is strongly adverse (+79% to +131%). Together, at the derived
  weights (1, 1), they remove a fifth of the MSE on every network. This is the grading's claim made in the chain: the
  second-chaos path class and the third-chaos star are complementary halves of the kappa4 diagonal's quenched error, and
  neither is a truncation of anything by itself.
- **The derived weights beat the one-step fit by 8.7 points.** The one-step regression's (1.15, 0.7) is attenuated
  by the candidates' own errors. Free running prefers the weights the theory gives.
- **Layers 1-2 must be left out, and the last layer must have both pieces.** With the star alone there it costs more
  than it saves.
- **The closure's memory channel is negligible** (m1fit = m2fit). The post-activation kappa4 diagonal is dominated by
  the Gaussian nonlinearity's own fourth cumulant.

**The production version** (`V56_OLD=1`: the old tier's factor families split the same way plus one r-lift, and the
star over every slot, since the old tier's legs are formed dense). On network 0 against the production fold system:
- full solve: 1.789e-8 against 2.320e-8, **-22.9% raw**, FLOPs +12.8% (0.2215 B against 0.1964 B), adjusted about -13%;
- the Galerkin solve in Y's top row space (`V56_GAL=64`): +5.8%. Y is not concentrated in 64 dimensions, so the
  full solve stays.

## 9. The production screen (stated before its runs)

Networks 0-15, paired against the production fold system (B35 + `V52_FB_FOLD=1`) rerun in the same batch. All
variants use `V56_P4=2 V56_A=1 V56_B=1` with the full solve.

| tag | further switches | prediction |
|---|---|---|
| p3 | `V56_LMIN=3 V56_LAST=2` | **P27:** raw -15% to -25%, better on >= 15/16; FLOPs +12% to +14%; adjusted -7% to -14% |
| p3L0 | `V56_LMIN=3 V56_LAST=0` | P28: raw worse than p3 by 2-8 points; FLOPs +8% to +10% |
| p3m2 | `V56_LMIN=2 V56_LAST=2` | P29: within +-3 points of p3 |
| p3e | `V56_LMIN=3 V56_LAST=2 V56_EPS=0.001` | P30: within +-2 points of p3 |

If P27 holds, the winner goes to the scored regime: all 100 networks, counterterms refitted on 0-49, held-out 50-99,
adjusted MSE (note XXXIX protocol).

**The scored regime (stated before its runs).** This is the adopted system: V35 + pair (`V33_WK4M=3 V31_K4D=3`) +
fold, in the scored harness (`run_v29w.py`: warm-up predict, then the measured one). It runs on all 100 networks,
free-running without counterterms, plus `V56_P4=2 V56_A=1 V56_B=1 V56_LMIN=3 V56_LAST=2`. Each variant is paired
against the adopted system's own free-running scored outputs (note XXXIX `pe`); ten of those are rerun in the batch to
check that the code with V56 off is bit-identical.

| tag | | prediction |
|---|---|---|
| v56 | eps = 0.01 | **P31:** raw -15% to -25%, C/B +11% to +14%, adjusted -5% to -14%, better on >= 95/100 |
| v56e | eps = 0.001 | P32: raw within +-2 points of v56 |

The protocol then follows note XXXIX: counterterms refitted on networks 0-49 with the stored responses, then judged
free-running on 50-99 against the adopted system with its own counterterms (3.5090e-9).
- **P33:** the held-out adjusted MSE is at least 5% below 3.5090e-9.

## 10. The scored regime, measured (`outputs/scored_100nets/`)

The adopted system (V35 + pair + fold) plus `V56_P4=2 V56_A=1 V56_B=1 V56_LMIN=3 V56_LAST=2`, run in the scored
harness on all 100 networks, free-running without counterterms. It is paired against the adopted system's own
free-running scored runs (note XXXIX `pe`).
- The ten reruns of the adopted system in this batch reproduce those runs bit for bit on 10/10 networks: with V56
  off, the code is the adopted system.
- Sixteen runs were killed by the out-of-memory killer (48 slots of about 10 GB on one 384 GB instance) and rerun at
  20 slots.

| | raw | C/B | adjusted | per net (adjusted) | better on (raw / adjusted) |
|---|---|---|---|---|---|
| adopted, free-running | 2.1800e-8 | 0.17777 | 3.8752e-9 | | |
| **+ V56, eps = 0.01** | **1.6686e-8 (-23.46%)** | 0.20279 (+14.08%) | **3.3841e-9 (-12.67%)** | -12.34 +- 0.53 | 100/100 / 97/100 |
| + V56, eps = 0.001 | 1.6579e-8 (-23.95%) | 0.20279 (+14.08%) | 3.3625e-9 (-13.23%) | -12.90 +- 0.53 | 100/100 / 97/100 |

- **P31.** Raw, adjusted and the "better on" count hold. C/B is +14.08% against the predicted +11% to +14%,
  marginally outside.
- **P32 holds** (0.5 points).
- **The free-running V56 system is already below the adopted system with its counterterms**: 3.38e-9 against
  3.4923e-9 in-sample and 3.5090e-9 held-out.

**The cost** (profile of network 0, `profile_*_net0.txt`): +25.8 units of 2n^3.
- The 13 solves are 17.3 units (exactly flopscope's 2n^3/3 + 2n^3 each).
- The last layer's families, its covariance sandwich and the old tier's lifts are 7.2 units.
- The elementwise work is 1.4 units.
- The solves are two thirds of the price. The Galerkin shortcut failed (section 8), so the solve is the next thing
  to make cheaper.

**Residual time and memory.**
- The measured call's local residual rises from 0.634 s to 0.781 s mean (max 0.899 s). The grader counts only the
  Python, and measured 0.17 s for a system this harness put near 0.57 s (note XXIX), so the margin to 0.4 s must be
  checked on the grader's terms before submission.
- The memory high-water mark falls from 10.30 GB to 9.89 GB.
