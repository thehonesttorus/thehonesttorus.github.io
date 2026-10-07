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

## What transfers

| from the note | what it is here | decision |
|---|---|---|
| conditional mechanism of the omitted classes, weak-factor orders | explains note XXI's measured class structure | confirmation |
| quadratic residual, envelope identity E[h R^2] | an exact description of the omitted response | the response is at the noise floor here |
| relative determinant, Witten moving harmonic state | further exact descriptions of the same displacement | relabelling with proofs |
| modular cross-ratio bound on maximal correlation | a sufficient condition for memory decay on a conditional Markov family | the network is deterministic; the decay we have is measured, not certified |
| memory counterexamples (norm one, nilpotent delay, hidden global mode) | the formal record of note XIX section 1 | confirmation |
| balanced truncation with observability | a possible rank saving in the source compression | candidate; weak prior; measurable on the VM |
