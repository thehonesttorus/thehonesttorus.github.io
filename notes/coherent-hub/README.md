# The pair sector is not closed under the weights: the coherent doubly-occupied hub

Working note XLV. It reads three documents from the other line of work (stage 15 "deformed products, quasi-free lifts
and the error cocycle"; "Wick geometry, covariance holonomy, and the closure of Kikuchi dynamics"; "Branching, angular
geometry, and joint dynamics") against the measured anatomy of our estimator's error. It derives one concrete
consequence, a term that the regeneration of the fourth-cumulant pair slices has never contained, and pre-registers
a test before running it.

## 0. Where the error is, from measurements already made

These facts constrain any design, theory-native or not.

| fact | source |
|---|---|
| With every per-layer statistic the chain's local maps read made true (D3, D21, kappa4 diagonal, (2,2), (3,1)), the error falls 91-98% | note XXXIII F1 |
| Fed true inputs, the chain's next-layer fourth-cumulant slices are as wrong as when it runs free: diagonal 14-33%, (2,2) 16-32%, (3,1) 74-90%. The defect is made afresh by each layer's map | note XXXIII F2 |
| Removing the compression of the carried third-cumulant sources (exact dense legs) is worth about 8% | V21 lean reference, estimator header |
| The final error is incoherent: 0.4-3.7% of its energy along the mean direction | note XLIV 9d |
| Moving any slice's dilation amplitude toward the truth raises the MSE (2-56%) | note XLIV 9j |
| The truth's (3,1) slice is only 45-53% dilation-package shaped at deep layers; the chain's own is 96-97% package shaped | note XLIV 9k |
| Production builds the (3,1) slice as lambda C_off and the (2,2) slice as a rank-one var-class outer product | `estimator_final_v56.py`, regeneration block |

So the binding constraint is the per-layer regeneration map of the pair slices, not compression, not amplitudes
along the collective mode, and not the representation class.

## 1. What the three documents say, in this problem's terms

**Stage 15** (verified there on network 0):
- The error of any estimator is a cocycle: births at walls transported by the first jet.
- Gaussian rows make distinct Hermite multi-indices orthogonal; the classical form of Pauli-path orthogonality.
- Given the computed state, the anisotropic part of a birth has zero conditional mean, so only computation, not
  calibration, reaches it.

**Holonomy note.** Two statements carry over exactly.
- *Hard-core normalizer theorem.* An orthogonal substitution preserves the square-free quadratic sector only if it is a
  signed permutation. A generic rotation sends z_1 z_2 to cos(2 theta) z_1 z_2 + (sin(2 theta)/2)(z_1^2 - z_2^2): the
  doubly-occupied terms leak in.
  - The chain's pair slices D21, K22, K31 are exactly the hard-core (Kikuchi) sector in the neuron basis.
  - The weight matrix is the generic substitution between consecutive neuron bases.
  - So the pair sector at layer l+1 receives the doubly-occupied classes of layer l: kappa(h_a, h_a, .) for the third
    cumulant and kappa(h_a, h_a, h_c, h_d) for the fourth.
- *The invariant ladder.* The rotation-invariant operators on Gaussian Wick space are su(1,1):
  - K_- is the heat operator;
  - K_0 is the dilation (Euler) operator;
  - K_+ is the radial pair creation.
  This is Howe duality for (O(n), sl(2,R)). Two consequences:
  - The annealed (O(n)-averaged) dynamics acts on each harmonic sector only through these three operators.
  - Everything reachable without instance-specific computation lies in this sector: the dilation charge of notes XLIII
    and XLIV, the counterterms, and the isotropic births of stage 15's Corollary 3.3. This is why every amplitude
    correction along it has failed.

  Two further statements do not carry over at our precision:
  - *The all-degree sparsifier.* Its sample count is q ~ W log n / (beta n eps^2). The variance error the chain must
    beat is about 1e-3 relative, so eps ~ 1e-3 gives q >> n^2.
  - *The spherical sampling dynamics.* It is a Monte Carlo estimator. Its variance per evaluation is the output
    variance, 0.079, against a target of 1e-8.

**Branching note.**
- The product defect Gamma(A, B) is exactly what a closure that propagates separate means discards.
- Multiplicity is the information behind a label: a pair statistic K31_ij aggregates content arriving by several
  histories that transform differently at the next layer. One of these histories is the doubly-occupied site.

## 2. Coherence power counting, and the term it selects

Take the fresh row w_i of W_(l+1) (entries of variance 2/n), and the next-layer pre-activations z' = W h. Every slice
of z' is a contraction of a post-activation cumulant of h with rows. Under the row randomness:
- A contraction in which some upstream index a meets the row an odd number of times has mean zero. Its fluctuation
  is down by n^(-1/2): it is incoherent.
- A contraction in which index a meets row i twice carries W_ia^2 > 0 and adds coherently.

For the (3,1) slice K31'_ij = kappa(z'_i, z'_i, z'_i, z'_j):
- the diagonal, (3,1) and (2,2) classes of h, and the all-distinct stars, are incoherent (order n^(-3/2));
- the doubly-occupied class kappa(h_a, h_a, h_c, h_d) with a on two of the three i-slots is coherent and of order
  1/n, the order of K31 itself.

The pair state never holds this class: it is an n^3 object.

**Its leading form is exact and cheap.** Write α = mu/S, and let Φ and φ be the normal cdf and density at α. Apply
the wall-jet theorem, i.e. Stein's lemma twice on the Hermite expansion of relu at the true pre-activation law:

    kappa(h_a, h_a, h_c, h_d) = G_a  Phi_c C_ac  Phi_d C_ad,       G = 2 [ Phi (1 - Phi) - alpha phi Phi - phi^2 ]
    kappa(h_a, h_a, h_c)      = Q_a  Phi_c C_ac,                   Q = 2 m (1 - Phi),  m = E relu(z_a)

Contract them with the rows. P = W diag(Phi) C_off is Cov(z', z) in first jet. Then:

    H31_ij = 3 sum_a W_ia^2 G_a P_ia P_ja        (one n^3 product)
    H22_ij = sum_a G_a (W_ia^2 P_ja^2 + W_ja^2 P_ia^2)
    H4_i   = 6 sum_a W_ia^2 G_a P_ia^2
    HD21_ij = sum_a W_ia^2 Q_a P_ja,   HD3_i = 3 sum_a W_ia^2 Q_a P_ia

The third-cumulant pair births (HD21, HD3) are what the chain's D21 births already carry. The fourth-cumulant ones
(H31, H22, H4) have no counterpart in the production regeneration. Two remarks on what H31 is:
- It is not package-shaped. It is the cross-covariance squared, weighted by the instance's W_ia^2, against the
  package's proportionality to C'.
- It is quenched. The annealed replacement W_ia^2 -> 2/n keeps the shape (6/n) sum_a G_a P_ia P_ja. The difference is
  the fresh-weight fluctuation: Var(W_ia^2) / E(W_ia^2)^2 = 2, the identity of note XLIV T1. Section 3 measures both
  forms, so the test also measures directly what the estimator gains by coupling its structure to the given weights.

## 3. Pre-registration: one-step regeneration against the 1.6e7-input truth (committed before the data)

**Setup.**
- Networks 0 and 1, target layers s = 4, 6, 8, 10, 12, 13, 14, 15.
- The H terms are built from the true layer s-1 state (mu, C from the full Monte Carlo) and W_s.
- They are compared with the true slices at layer s and with the chain's own slices (chaindump2, s <= 14).
- All metrics use the neurons with true alpha > -2.5 (active by active for matrices, off-diagonal).
- Energies of truth and of chain error are noise-free (products across the two independent halves).
- Two metrics:
  - entrywise;
  - the slice's own readout into the next layer's variance, diag W_(s+1) T(X) W_(s+1)^T. Here T is the exact
    bivariate-Edgeworth map of `var_ladder.py`, taken at the true layer s state; this is the metric the output feels.
- Script: `code/hub_test.py`.

**Predictions** (layers s >= 8, both networks):
1. In the joint fit truth K31 ~ c_p f31(package) + c_h H31, the coefficient c_h lies in [0.5, 1.5].
2. Regressing the chain's (3,1) error (truth minus chain) on H31 gives a coefficient in [0.5, 1.5] and explains at least
   25% of that error's energy in the variance-readout metric.
3. H22 explains at least 10% of the chain's (2,2) error energy, and H4 at least 10% of its diagonal error energy, at
   coefficients in [0.5, 1.5].
4. Control: HD21 explains at most 10% of the chain's D21 error, because the chain carries these births.
5. In predictions 1-3 the quenched form explains more than the annealed form.

Prior, stated in advance: H31 has a collective part along the dilation direction, through the top eigenvector of C.
That part is partly degenerate with the package, which could leave little for H31 to explain beyond it. I give
prediction 2 about 45%.

**Decision.**
- If prediction 2 holds on both networks: implement H31 (and H22/H4 where 3 holds) as a parameter-free regeneration
  term, one n^3 product per layer, behind a switch. Screen it cold on networks 0-7 with counterterms on and off.
- If prediction 2 fails: the coherent doubly-occupied class is not the missing (3,1) content, and the note records
  what the residual is instead.

## 4. Pre-registration of the build gate: all slices true at once, late layers only (committed before the runs)

F1 (note XXXIII) measured the all-slice, all-layer oracle on an older system. The build decision needs the same on
the production estimator, and its split by depth.

**Setup.**
- Production estimator with counterterms. V37_ORACLE=D3,D21,G4,WK4M,K31 (all five readout slices true together),
  with the consistent wiring V41_ORC_CONSIST=1: the legs keep carrying the chain's own content.
- Network 0, Monte Carlo full and halves. Noise-extrapolated a = 2 MSE(full) - mean(MSE(h0), MSE(h1)), against
  1.4771e-8.
- (A) layers 10-14; (B) every layer.

**Expected.** (B) -80% to -95%, i.e. F1 replicates on the production system. (A) -40% to -70%.

**Gate.**
- If (A) is at least -40%, the target of the next system is a consistent recomputation of the late slices. It is
  priced against the 2x that the adjusted score allows (C/B from 0.2 to at most 0.4 for a 2x raw gain, break even).
- If (A) is under -25% while (B) is large, the slices must be right at every depth, and the cost question is global.
