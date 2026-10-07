# Cumulant geometry and controlled memory, placed

Working note XXIV. A 62-page continuation of the theory notes was uploaded (`ncg_geometry_and_memory.pdf`,
checkpoint NCG-20261007-B). It derives: the conditional mechanism of the omitted fourth-cumulant classes in a
Gaussian factor reference (law of total cumulance: a pair source perturbs conditional covariances, whose connected
factor correlations with other neurons make the omitted classes); their orders under weak common factors ((3,1) at
O(1), (2,1,1) at O(eps^2), (1,1,1,1) at O(eps^4)); a quadratic-residual decomposition kappa_4(Y) = t' C^-1 t +
E R^2 - 2 v^2 with the envelope identity d E R^2 = E[h R^2]; a relative Fredholm determinant whose fourth
coefficient is the fourth cumulant on the probability vacuum; a Witten deformation in which the Gaussian Dirac index
stays one while the harmonic state moves by h/2; a sharp modular cross-ratio bound on maximal correlation, r <=
tanh(Delta/4); a balanced-truncation reduction for a nonautonomous linear system with the certificate sum_l
sigma_{r_l+1}(O_l W_l^{1/2}); a symmetry-channel split; and counterexamples showing that homogeneity plus a spectral
gap, adjacent one-step kernels, or unary smoothing do not imply short memory (sections 12-16). It states plainly that
none of it is evaluated on the two networks.

## What it confirms

**The class structure we measured has its mechanism.** Note XXI decomposed the fourth cumulant of the next
pre-activation by index class on the Monte Carlo law: (3,1) zero, (2,1,1) the whole omitted part and equal to the
scale mixture's 6 g s_diag^2 s_off^2, (1,1,1,1) at 3 g s_off^4 and within noise. The weak-factor orders of section 8
are exactly this ordering once eps^2 is read as the common-mode share of the variance: the (2,1,1) class is first
order in the off-diagonal variance s_off^2, the (1,1,1,1) class second order, and the (3,1) class, though O(1) per
entry, enters the coherent projection with odd powers of the weights and averages out. The law-of-total-cumulance
derivation is the conditional explanation of a measured fact.

**The memory counterexamples are our record.** Theorem 14.2 (a deterministic layer has score-map norm exactly one),
section 15 (homogeneity plus a gap does not give short memory) and Proposition CL.3 (a persistent global mode hidden
from adjacent kernels) are the formal versions of note XIX section 1: the response of the output to content born at
layer b decays geometrically on the third-cumulant bulk (0.80-0.92 per layer) and not at all on the scale and mean
direction, which is why the gain had to be carried as a global variable. The symmetry split of section 12 is that
global variable.

## What it changes

Nothing in the chain, for the reason given in notes XXII and XXIII: the quantity whose mechanism and evaluation the
programme refines, the response of the omitted classes to the incoming third cumulant, is at the noise floor on these
networks (at most 2% of g4, note XXI section 5), and the fourth-cumulant feed that matters is in the retained classes,
where the chain's first-order programs are already exact. The quadratic residual E[h R^2], the relative determinant
and the moving harmonic state are three exact descriptions of that same small response.

**One candidate lever, with an honest prior.** The balanced-truncation theorem of section 11 chooses the reduced
coordinates by the singular values of O_l W_l^{1/2}, reachability and observability together. The chain's
compression of its third-cumulant sources (the shared basis of rank 384 and the nested tier of 224) uses a range
finder on the transported content, which is reachability only. Weighting by observability, the sensitivity of the
final output means to each direction of carried content, could keep fewer ranks at equal error. The prior is weak:
the output readout is all 1024 means with equal weight, the bulk transport is marginal (He criticality), and the
error anatomy of note XVIII found the output error white across neurons and uniform in alpha, which is what an
isotropic observability Gramian produces. Where observability is not isotropic is the last two layers, through the
saturated gates (38% of output neurons have alpha < -1 and their means barely respond to third-cumulant content), and
the chain already trims the last layer. A test would need the chain's linearised output sensitivity to its carried
content at an intermediate layer, which is a dump of the legs plus a backward pass through the chain's term programs
on the Azure VM; it is recorded here as the one item of this note worth a measurement, behind the bounded levers
already queued.

## Measurement: observability at the gateway (V32, Azure VM, `code/v32_join_post.diff`, `outputs/`)

**The prior above was wrong in one respect, and the correction is exact.** Isotropy of the observability Gramian
holds on average over the weights, not for the given weights. For one square He matrix, W^T W has the
Marchenko-Pastur spectrum of ratio one, spread over [0, 4] times its mean, so about a fifth of all directions are
damped more than tenfold by a single layer. And the chain truncates in exactly the wrong coordinates: the join at
layer l re-projects the joiner's legs and the old core onto rank 384 in the post-gate rows of layer l - 1, and every
later use reads the legs only after W_l (Qc = W_l Qn, legs A_s = Qc FA_s). The error that matters is
W_l (I - Pi) X, so the one-step observability Gramian of the join is W_l^T W_l, with no approximation. Two
implementations, both behind flags in `est_v29.py`:

- `V32_JOIN_POST=1`: the range finder and the truncation in post-W coordinates (Qn' = range of W G W^T, factors
  Qn'^T W X, Qc = Qn'). Three extra n x n x r products per join.
- `V32_JOIN_POST=2`: keep the pre-W subspace, orthonormalise W Qn and project W X onto it orthogonally (factors
  (W^T Qt)^T X, rotation (W^T Qt)^T Qp). One extra product and a QR. For a fixed subspace the orthogonal projection
  in the metric where the error is read is optimal (Pythagoras), so this variant never loses in that metric.
- `V32_JOIN_SMM=1`: the join's n x n x r products (range finder, factor projections, basis transport) through the
  Strassen family the rest of the chain already uses (they were plain matmuls). Pure arithmetic.

Sixteen official networks, paired against the shipped configuration (best100), all variants in parallel on the VM:

| variant | raw | C/B | adjusted |
|---|---|---|---|
| post-W range finder + truncation, rank 384 | -0.76% +- 0.43 (10/16) | 0.2672 | +3.34% |
| post-W range finder + truncation, rank 320 | +4.21% +- 0.71 | 0.2459 | -0.14% |
| shipped join, rank 320 (control) | +8.40% +- 0.57 | 0.2370 | +0.12% |
| post-W projection only, rank 384 | -1.13% +- 0.39 (12/16) | 0.2627 | +1.22% |
| post-W projection only, rank 320 | +3.59% +- 0.58 | 0.2418 | -2.38% |
| Strassen join products only, rank 384 | -0.06% +- 0.10 | 0.2468 | -3.88% |
| Strassen join + post-W projection, rank 384 / 352 / 320 / 288 | -1.13% / +0.91% / +3.61% / +9.42% | 0.2515 / 0.2418 / 0.2326 / 0.2238 | -3.09% / -4.91% / -6.08% / -4.57% |

Three readings. The whole effect is the projection metric, not the subspace: keeping the pre-W subspace and
projecting after W beats the full post-W range finder at every rank, at a third of the cost (one pass of a
warm-started finder already sits on the forward Lyapunov subspace; tilting its sketch by W^T W only moves it off).
At rank 384 the projection removes 1.13% of the 1.4% that the whole shared basis costs (note XVII), and at rank 320
it halves the truncation loss (+8.4% to +3.6%), which moves the rank optimum down. My cost estimate for the
first variant ("about 1% of the budget") was 1% of B, i.e. 4% of the chain's bill; the projection-only variant fixes
that.

**Validation on all 100 networks** (Strassen join products, post-W projection, shared basis 320, nested 192;
`outputs/v32_100nets.txt`): raw +4.30% +- 0.31, C/B 0.2326, adjusted **-5.45% +- 0.31**, mean adjusted
**5.432e-9** against 5.742e-9 for the shipped chain. Ranks 336 and 304 on nets 0-15 give -5.88% and -5.34% against
-6.08% for 320 on the same nets, so 320 is the optimum.

**The renormalised metric (`V32_JOIN_POST=3`, `outputs/v32_renorm_metric.txt`).** The post-W metric is still flat
across the rows of layer l, but a row's content reaches every later layer through its own gate (the legs are scaled
by w1 = Phi(alpha) before the next transport), and its local reads (D3, D21 at layer l) are what is left. So the
metric in which the join's error is read is diag(om) with om_i = Phi(alpha_i)^2 + c: Mallat's renormalisation of
the wavelet coefficients (note XXV), with the weight set by observability instead of variance. alpha is exact for
the mean (mu is formed before the join) and uses the diagonal transport of the previous post-activation variance
for sigma, which is enough for a weight. Same cost as the plain projection (a row scaling and a weighted QR). With
the factor rotations and the D21 lift also through the Strassen family (`V32_ROT_SMM=1`, raw unchanged, 1.3% of
the bill saved):

| variant (rank 320 / nested 192 unless stated) | nets | raw | C/B | adjusted |
|---|---|---|---|---|
| post-W projection + Strassen rotations | 16 | +3.61% | 0.2295 | -7.34% |
| renormalised metric, c = 0.1 | 16 | +2.20% +- 0.49 | 0.2296 | -8.56% |
| renormalised metric, c = 0.5 / 0.03 | 16 | +2.73% / +2.73% | 0.2296 | -8.08% / -8.08% |
| renormalised metric, c = 0.1, rank 288 / 352 | 16 | +7.47% / +0.12% | 0.2213 / 0.2383 | -7.32% / -7.02% |
| **renormalised metric, c = 0.1** | **100** | **+2.98% +- 0.27** | **0.2296** | **-7.86%; mean adjusted 5.2925e-9** |

This is the new best configuration: `V29_WARM_JOIN=1 V29_WARM_FB=1 V17_R_RES=4 V26_STRASSEN=6 V26_STRASSEN_MIN=16
V32_JOIN_SMM=1 V32_ROT_SMM=1 V32_JOIN_POST=3 V32_JP_C=0.1 V21_R_OLD=320 V24_R_OLD2=192`, mean adjusted 5.2925e-9
on the 100 official networks against 5.742e-9 shipped. The renormalisation is worth 1.4 points of raw at zero
cost on top of the flat post-W projection, so the observability content of the join is not exhausted by the
gateway W; the gate is the second factor of the same Gramian.

## What transfers

| from the note | what it is here | decision |
|---|---|---|
| conditional mechanism of the omitted classes, weak-factor orders | explains note XXI's measured class structure | confirmation |
| quadratic residual, envelope identity E[h R^2] | an exact description of the omitted response | the response is at the noise floor here |
| relative determinant, Witten moving harmonic state | further exact descriptions of the same displacement | relabelling with proofs |
| modular cross-ratio bound on maximal correlation | a sufficient condition for memory decay on a conditional Markov family | the network is deterministic; the decay we have is measured, not certified |
| memory counterexamples (norm one, nilpotent delay, hidden global mode) | the formal record of note XIX section 1 | confirmation |
| balanced truncation with observability | the join's projection taken in the post-W metric, where the legs are read | measured: post-W projection -5.45%, gate-renormalised metric + Strassen rotations -7.86% adjusted on 100 networks (new best, 5.2925e-9) |
