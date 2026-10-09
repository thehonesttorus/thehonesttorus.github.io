[?] Deterministic estimation of the expected outputs of deep neural networks under Gaussian input, without sampling: cumulant propagation; analytic moment propagation / assumed-density filtering / moment matching for Bayesian neural networks; deterministic variational inference; Edgeworth and Gram-Charlier corrections; finite-width 1/n corrections to NNGP kernel recursions (Yaida; Roberts-Yaida-Hanin 'Principles of Deep Learning Theory'; Dyer-Gur-Ari Feynman diagrams; Hanin); Gaussian-mixture and Rao-Blackwellized/conditional moment propagation; unscented and sigma-point transforms in high dimension; and the Alignment Research Center program on heuristic estimators ('Formalizing the presumption of independence' with cumulant propagation, 'Towards a law of iterated expectations for heuristic estimators', mechanistic estimation, estimating tail risk / low-probability estimation). Which methods reach errors far below Monte Carlo at equal FLOPs for width-1024 depth-16 ReLU MLPs, and what limits them? Hypothesize how localization (conditioning on a few shared input directions or collective modes such as the input radius or a gain mode, revealing activation patterns layer by layer, iterated-expectation consistency) could give a strictly better cost-accuracy frontier.

[?] Finding evidence on deterministic neural-network moment estimation

[?] Working

[?] Working

[?] Working

[?] Working

[?] Indexed 38 sources on neural-network moment estimation

[?] Working

[?] Working

[?] Working

[?] Reviewed indexed sources on deterministic neural-network moments

[?] Reviewed a source on neural-network tail-risk estimation

[?] Working

[?] Reviewed sources on mechanistic estimation for wide random networks

[?] Working

[?] Working

[?] Reviewed wide-MLP moment estimation and covariance propagation sources

[?] Finding evidence on deterministic neural-network moment estimation

[?] Separating demonstrated FLOP advantages from finite-width extrapolations

[?] Indexed 38 sources and retrieved material on neural-network tail risk, mechanistic estimation for wide random MLPs, and estimating their expected outputs without sampling.

[?] Working

[?] Working

[?] Working

[?] Indexed 44 sources; one paper’s full text was unavailable

[?] Working

[?] Reviewed indexed sources and searched the wide-network study for limitations

[?] Working

[?] Working

[?] Working

[?] Reviewed source index and extracted neural-network efficiency passages

[?] Working

[?] Reviewed variational inference source and covariance propagation summaries

[?] Separating demonstrated FLOP advantages from finite-width extrapolations

[?] Synthesizing the cost–accuracy frontier and a conditional-localization hypothesis

[?] Indexed 44 sources and examined citation records; full text for “Non-Gaussian Processes and Neural Networks at Finite Widths” could not be retrieved.

[?] Working

[?] qShANRCA:48-49,175-206,211-221,228-230,245-258,789-792,815-816

[?] qShANRCA:48-49,175-206,211-221,228-230,245-258,789-792,815-816

[?] Working

[?] Resolved evidence on efficient expected-output estimation for wide random MLPs

[?] Working

[?] Prepared methods section comparing evidence and estimands

[?] **Factorized, augmented cumulant propagation is the strongest demonstrated candidate—but the evidence does not establish a winner for width-1024, depth-16 ReLU MLPs.** ARC’s 2026 work directly addresses the fixed-network Gaussian-input expectation and demonstrates substantial equal-FLOP advantages at smaller architectures. Its experiments cover widths up to 256 and hidden-layer counts up to 12; depth is also the main weakness of its accuracy scaling. Consequently, “far below Monte Carlo” at your requested architecture is a plausible extrapolation, not an established result. <citations>qShANRCA:48-49,202-218,228-230,789-792</citations>

The distinction between **randomly initialized** and **trained** networks is decisive. The demonstrated advantage concerns particular networks drawn from an initialization ensemble, with estimation error then averaged over initializations. It is not a guarantee for arbitrary trained weights. <citations>qShANRCA:102-105,410-415</citations>

## 1. What actually beats sampling—and what does not yet qualify

**Cumulant propagation with the right contractions and factorization.** ARC propagates low-order activation cumulants, uses Hermite expansions to handle nonlinearities, and contracts the resulting terms through the actual weight matrices. The important improvements are not simply “keep skewness and kurtosis”:

- **Power cumulants** handle repeated neuron indices correctly, including exact Gaussian marginal variances at covariance order.
- **Trace information beyond the nominal order** retains collective fluctuations that ordinary truncation would discard.
- **Factorized higher-order tensors** reduce computational cost without changing the algorithm’s output; augmentation retains additional trace-containing components. <citations>qShANRCA:124-160,175-207</citations>

At width 256, the factorized variants report approximately **1,000-fold advantages for two hidden layers and 10–100-fold advantages for four hidden layers** in their MSE–FLOP comparisons. But at eight hidden layers, the fourth-order variants begin to underperform sampling. This is strong evidence for an advantageous region of the frontier, not uniform superiority. <citations>qShANRCA:254-258</citations>

The theorem is narrower than the headline: it assumes **polynomial activations, fixed depth, independent Gaussian initialization, and width tending to infinity**. In the appropriate tolerance regime, factorization gives runtime proportional to \(n/\varepsilon^2\), versus \(n^2/\varepsilon^2\) for sampling. ReLU performance is supported empirically rather than by that theorem as stated. <citations>qShANRCA:98-99,211-223</citations>

**Analytic Gaussian moment propagation / ADF / moment matching.** These are useful baselines, especially with full covariance rather than diagonal variance. Analytic covariance formulas can make the nonlinear Gaussian-to-moment step exact or arbitrarily precise. They do **not** make an entire deep network’s propagation exact: the next layer generally receives a non-Gaussian distribution that has been replaced by a Gaussian. The retrieved covariance work explicitly distinguishes the exact Gaussian-input covariance calculation from the Gaussian assumptions used in propagation. <citations>9auoQGRg:16-19,57-60</citations>

My assessment is that this family can give excellent cheap estimates, but its unresolved closure bias prevents treating it as a route to arbitrarily small error at fixed width and depth. Full covariance also changes the economics: dense covariance propagation through a square affine layer has cubic cost, versus quadratic cost for one forward pass. The covariance paper reports this same cubic-versus-sample-count distinction for DVI and Monte Carlo VI. <citations>9auoQGRg:115</citations>

**Deterministic variational inference for BNNs.** DVI removes Monte Carlo gradient noise by approximating neural-network moments deterministically and has demonstrated predictive benefits. That is not the same benchmark as accurately integrating a particular fixed-weight MLP over Gaussian inputs. A method can improve posterior optimization while retaining substantial moment-closure or posterior-approximation error. The retrieved DVI paper establishes deterministic moment approximation and regression performance, not a matched-cost victory on your architecture. <citations>fpUhKhxQ:3-7</citations>

**Edgeworth / Gram–Charlier corrections.** These are better viewed as expansion machinery than as standalone competitors. The useful question is which cumulants and contractions are retained, how the nonlinear step is evaluated, and whether the representation remains cheap. ARC’s algorithm supplies these missing implementation and error-analysis ingredients through cumulants, Hermite coefficients, and diagram sums. <citations>qShANRCA:124-160,169-207</citations>

As a mathematical caution, a truncated density correction need not remain positive, and small central-distribution error does not imply small relative tail-probability error. For ReLU, the kink also makes a naive smooth-function Taylor argument inappropriate. Expanding observables in a Gaussian-orthogonal basis is often preferable to constructing and then integrating a possibly invalid approximate density—but still requires control of omitted terms.

**Finite-width NNGP corrections and Feynman-diagram theory.** The ensemble/instance distinction matters more than the expansion order. Yaida’s finite-width construction produces weakly non-Gaussian **priors** by integrating over random weights; Dyer–Gur-Ari’s diagrammatic work analyzes wide-network asymptotics and training corrections. These are valuable tools, but an ensemble kernel recursion is not itself an estimator of \(\mathbb E_X[f_W(X)]\) for the supplied weights. <citations>wMwwWCwA:1-4,rLl5jZYg:1-4,qShANRCA:444-445</citations>

For example, averaging a zero-mean final weight matrix over the ensemble can annihilate a mean that is nonzero for an individual realization. Adding an ensemble finite-width correction does not automatically recover that realization-specific mean. The theory becomes operational here when translated into **weight-dependent contractions**, as in quenched cumulant propagation.

**Gaussian mixtures and conditional propagation.** Mixtures are a genuine alternative to a single Gaussian closure. Retrieved work propagates a split-and-merge Gaussian mixture, uses a Wasserstein splitting criterion, and reports convergence guarantees and accurate output-density estimation. However, it does not establish the requested equal-FLOP comparison. <citations>LPQfXbcQ:1-6</citations>

The structural limit is component growth: approximating many intersecting activation boundaries can require many components, while merging them can erase precisely the dependencies needed downstream. A mixture over a few shared latent modes looks much more promising than indiscriminate splitting in the full activation space.

**Unscented and sigma-point transforms.** I would treat these as deterministic quadrature baselines, not established winners here. Their weakness follows from their construction: matching a few polynomial moments does not ensure that a limited collection of points resolves many ReLU boundaries. Tensor quadrature becomes expensive in high dimension; sparse rules help only when the integrand has favorable effective dimension or interaction structure. The promising role is **sigma points over selected collective coordinates**, with analytic integration of the residual—not sigma points over every input coordinate.

## 2. Why width 1024 and depth 16 remain uncertain

ARC conjectures an error scale of the form
\[
\mathrm{MSE}_K\lesssim c_K(L/n)^K,
\]
with empirical depth scaling roughly consistent with it. Its factorized runtime bound has a leading dependence \(L^2n^K\) for higher orders, rather than the basic algorithm’s \(Ln^{K+1}\). These are order-dependent statements with unknown constants, not numerical error forecasts. <citations>qShANRCA:219-221,789-792</citations>

For your architecture, \(L/n\) is small, which favors a width expansion. But depth can increase both approximation error and the cost of maintaining the factorized representation. Increasing cumulant order may therefore buy less accuracy per FLOP than at shallow depth. **I would benchmark the lower-order and augmented variants before assuming that the highest affordable order wins.**

There are three further limits:

1. **Fixed-width accuracy:** increasing the order is not known to give a sampling-competitive frontier as tolerance tends to zero at fixed width. ARC explicitly leaves this regime open. <citations>qShANRCA:412</citations>
2. **Training-induced structure:** trained networks can amplify selected higher-order dependencies that an initialization-based expansion regards as negligible. ARC identifies tracking those deviations as an essential unresolved extension. <citations>KvXnqCgQ:55-62</citations>
3. **FLOPs versus implementation:** the paper adjusts FLOP counts to remove redundant symmetric-tensor operations without implementing corresponding optimized kernels. The authors also report that their methods often lose in wall-clock time. Equal-FLOP superiority should not be read as current GPU superiority. <citations>qShANRCA:249-251,KvXnqCgQ:76-77</citations>

## 3. Localization: a plausible strictly better frontier

My hypothesis is **conditional, output-aware cumulant propagation with persistent shared coordinates**. The gain would come from representing structured dependence compactly rather than raising the global tensor order.

### Condition on shared directions, then integrate the residual analytically

For isotropic Gaussian input, choose orthonormal columns \(U\) and write
\[
X=Uz+\xi,\qquad z\sim\mathcal N(0,I),\quad
\xi\sim\mathcal N(0,I-UU^\top),
\]
with \(z\) and \(\xi\) independent. Then
\[
\mathbb E[f(X)]=\mathbb E_z[m(z)],\qquad
m(z)=\mathbb E_\xi[f(Uz+\xi)].
\]

Propagate conditional moments for the residual and integrate \(z\) using deterministic quadrature. Candidate directions should capture **shared gate changes and their effect on the requested output**, not merely explain input or activation variance. Preserve the same coordinates across layers; re-Gaussianizing them independently at each layer would discard the intended benefit.

Why could this help? Even if coordinates are nearly Gaussian conditional on a common mode, mixing over that mode produces higher-order dependence. For example, common conditional variance fluctuations generate fourth-order cross-cumulants. A low-dimensional latent representation can encode these collectively rather than storing a large tensor.

The qualification is important: conditioning improves the *representation opportunity*, not automatically the estimator. Bad directions, inaccurate conditional closure, or expensive quadrature can all make the frontier worse.

### Separate input radius—but do not confuse radial and angular error

There is an exact simplification for the bias-free ReLU architecture. Write \(X=RS\), with independent Gaussian radius \(R\) and uniform spherical direction \(S\). Positive homogeneity gives
\[
f(RS)=R f(S),\qquad
\mathbb E[f(X)]=\mathbb E[R]\,\mathbb E[f(S)].
\]

This is a mathematical identity, not a proposed approximation. Radius can therefore be integrated exactly before estimating the angular expectation. It prevents the propagation algorithm from repeatedly approximating the same radial fluctuation, but **does not solve the angular gate problem**. Biases break this simple factorization. A deeper-layer gain mode is potentially more valuable when it changes gate probabilities or carries angular dependence, rather than merely reproducing input radius.

### Reveal selected activation patterns—not all patterns

Conditioning on a complete gate pattern makes a ReLU network affine inside its activation region. Exact integration then becomes a Gaussian polyhedral-probability and truncated-moment problem. That shifts the difficulty rather than eliminating it: the number of regions and their constraints can be enormous.

A practical localization scheme would instead branch on gates or gate groups with large **output-relevant approximation error**, retain their conditioning constraints, and leave the remaining computation under a conditional cumulant closure. Stable gates can be treated as locally affine. Merely forcing a gate on or off without updating the conditional input distribution would not constitute correct conditioning.

### Enforce iterated-expectation consistency

For a branch event \(A\), compare the parent estimate with
\[
\widehat p_A\widehat\mu_A+
(1-\widehat p_A)\widehat\mu_{A^c}.
\]

A discrepancy is a useful refinement signal: the estimator can already predict that its own answer changes under a tractable decomposition. But agreement is **not an error certificate**; parent and children can share a bias.

This is an implementable analogue of ARC’s unpredictable-errors motivation, not a solved consequence of its heuristic-estimator theory. The iterated-expectations work identifies barriers and leaves the appropriate formalization open. Separately, ARC’s tail-risk discussion explicitly warns that layer-by-layer projection can lose shared computational dependencies and create correlated errors. <citations>PCUQgjwQ:21-36,riRJOfag:238-255</citations>

### What would establish strict improvement?

For deterministic quadrature weights \(w_j\), conditional closure error \(e(z_j)\), and quadrature error \(e_Q\), the resulting error decomposes as
\[
\widehat\mu-\mu=e_Q+\sum_j w_j e(z_j).
\]

A strictly better frontier requires conditioning to reduce the *integrated closure error* enough to repay direction discovery, conditional propagation, and quadrature. Signed-error cancellation can help, but cannot be assumed.

The decisive experiment would compare this localized scheme with the best factorized/augmented cumulant variants, full-covariance closure, and strong integration baselines on identical networks and budgets. It should include ablations for radius alone, shared directions, gain modes, selective gate branching, and consistency refinement. If localization only helps after expensive uncounted mode discovery, it has not improved the frontier.

**My strongest bet is persistent shared-mode conditioning plus selective trace/cumulant corrections—not a global Gaussian mixture or a blanket increase in cumulant order.** This directly targets the dependence that depth repeatedly transports and Gaussian projection repeatedly forgets. For randomly initialized networks it may improve constants or depth scaling; for trained networks it offers a way to capture structured deviations that the initialization expansion misses. Neither improvement is yet demonstrated for your architecture.

Finally, mean accuracy and tail accuracy must remain separate objectives. ARC reports promising low-probability estimates, but its earlier tail-risk experiments had not beaten importance sampling. Tail claims need comparison against importance sampling and event-specific validation, not ordinary Monte Carlo alone or accurate central moments. <citations>KvXnqCgQ:28-30,riRJOfag:68</citations>

[?] **Factorized, augmented cumulant propagation is the strongest demonstrated candidate—but the evidence does not establish a winner for width-1024, depth-16 ReLU MLPs.** ARC’s 2026 work directly addresses the fixed-network Gaussian-input expectation and demonstrates substantial equal-FLOP advantages at smaller architectures. Its experiments cover widths up to 256 and hidden-layer counts up to 12; depth is also the main weakness of its accuracy scaling. Consequently, “far below Monte Carlo” at your requested architecture is a plausible extrapolation, not an established result. <citations>qShANRCA:48-49,202-218,228-230,789-792</citations>

The distinction between **randomly initialized** and **trained** networks is decisive. The demonstrated advantage concerns particular networks drawn from an initialization ensemble, with estimation error then averaged over initializations. It is not a guarantee for arbitrary trained weights. <citations>qShANRCA:102-105,410-415</citations>

## 1. What actually beats sampling—and what does not yet qualify

**Cumulant propagation with the right contractions and factorization.** ARC propagates low-order activation cumulants, uses Hermite expansions to handle nonlinearities, and contracts the resulting terms through the actual weight matrices. The important improvements are not simply “keep skewness and kurtosis”:

- **Power cumulants** handle repeated neuron indices correctly, including exact Gaussian marginal variances at covariance order.
- **Trace information beyond the nominal order** retains collective fluctuations that ordinary truncation would discard.
- **Factorized higher-order tensors** reduce computational cost without changing the algorithm’s output; augmentation retains additional trace-containing components. <citations>qShANRCA:124-160,175-207</citations>

At width 256, the factorized variants report approximately **1,000-fold advantages for two hidden layers and 10–100-fold advantages for four hidden layers** in their MSE–FLOP comparisons. But at eight hidden layers, the fourth-order variants begin to underperform sampling. This is strong evidence for an advantageous region of the frontier, not uniform superiority. <citations>qShANRCA:254-258</citations>

The theorem is narrower than the headline: it assumes **polynomial activations, fixed depth, independent Gaussian initialization, and width tending to infinity**. In the appropriate tolerance regime, factorization gives runtime proportional to \(n/\varepsilon^2\), versus \(n^2/\varepsilon^2\) for sampling. ReLU performance is supported empirically rather than by that theorem as stated. <citations>qShANRCA:98-99,211-223</citations>

**Analytic Gaussian moment propagation / ADF / moment matching.** These are useful baselines, especially with full covariance rather than diagonal variance. Analytic covariance formulas can make the nonlinear Gaussian-to-moment step exact or arbitrarily precise. They do **not** make an entire deep network’s propagation exact: the next layer generally receives a non-Gaussian distribution that has been replaced by a Gaussian. The retrieved covariance work explicitly distinguishes the exact Gaussian-input covariance calculation from the Gaussian assumptions used in propagation. <citations>9auoQGRg:16-19,57-60</citations>

My assessment is that this family can give excellent cheap estimates, but its unresolved closure bias prevents treating it as a route to arbitrarily small error at fixed width and depth. Full covariance also changes the economics: dense covariance propagation through a square affine layer has cubic cost, versus quadratic cost for one forward pass. The covariance paper reports this same cubic-versus-sample-count distinction for DVI and Monte Carlo VI. <citations>9auoQGRg:115</citations>

**Deterministic variational inference for BNNs.** DVI removes Monte Carlo gradient noise by approximating neural-network moments deterministically and has demonstrated predictive benefits. That is not the same benchmark as accurately integrating a particular fixed-weight MLP over Gaussian inputs. A method can improve posterior optimization while retaining substantial moment-closure or posterior-approximation error. The retrieved DVI paper establishes deterministic moment approximation and regression performance, not a matched-cost victory on your architecture. <citations>fpUhKhxQ:3-7</citations>

**Edgeworth / Gram–Charlier corrections.** These are better viewed as expansion machinery than as standalone competitors. The useful question is which cumulants and contractions are retained, how the nonlinear step is evaluated, and whether the representation remains cheap. ARC’s algorithm supplies these missing implementation and error-analysis ingredients through cumulants, Hermite coefficients, and diagram sums. <citations>qShANRCA:124-160,169-207</citations>

As a mathematical caution, a truncated density correction need not remain positive, and small central-distribution error does not imply small relative tail-probability error. For ReLU, the kink also makes a naive smooth-function Taylor argument inappropriate. Expanding observables in a Gaussian-orthogonal basis is often preferable to constructing and then integrating a possibly invalid approximate density—but still requires control of omitted terms.

**Finite-width NNGP corrections and Feynman-diagram theory.** The ensemble/instance distinction matters more than the expansion order. Yaida’s finite-width construction produces weakly non-Gaussian **priors** by integrating over random weights; Dyer–Gur-Ari’s diagrammatic work analyzes wide-network asymptotics and training corrections. These are valuable tools, but an ensemble kernel recursion is not itself an estimator of \(\mathbb E_X[f_W(X)]\) for the supplied weights. <citations>wMwwWCwA:1-4,rLl5jZYg:1-4,qShANRCA:444-445</citations>

For example, averaging a zero-mean final weight matrix over the ensemble can annihilate a mean that is nonzero for an individual realization. Adding an ensemble finite-width correction does not automatically recover that realization-specific mean. The theory becomes operational here when translated into **weight-dependent contractions**, as in quenched cumulant propagation.

**Gaussian mixtures and conditional propagation.** Mixtures are a genuine alternative to a single Gaussian closure. Retrieved work propagates a split-and-merge Gaussian mixture, uses a Wasserstein splitting criterion, and reports convergence guarantees and accurate output-density estimation. However, it does not establish the requested equal-FLOP comparison. <citations>LPQfXbcQ:1-6</citations>

The structural limit is component growth: approximating many intersecting activation boundaries can require many components, while merging them can erase precisely the dependencies needed downstream. A mixture over a few shared latent modes looks much more promising than indiscriminate splitting in the full activation space.

**Unscented and sigma-point transforms.** I would treat these as deterministic quadrature baselines, not established winners here. Their weakness follows from their construction: matching a few polynomial moments does not ensure that a limited collection of points resolves many ReLU boundaries. Tensor quadrature becomes expensive in high dimension; sparse rules help only when the integrand has favorable effective dimension or interaction structure. The promising role is **sigma points over selected collective coordinates**, with analytic integration of the residual—not sigma points over every input coordinate.

## 2. Why width 1024 and depth 16 remain uncertain

ARC conjectures an error scale of the form
\[
\mathrm{MSE}_K\lesssim c_K(L/n)^K,
\]
with empirical depth scaling roughly consistent with it. Its factorized runtime bound has a leading dependence \(L^2n^K\) for higher orders, rather than the basic algorithm’s \(Ln^{K+1}\). These are order-dependent statements with unknown constants, not numerical error forecasts. <citations>qShANRCA:219-221,789-792</citations>

For your architecture, \(L/n\) is small, which favors a width expansion. But depth can increase both approximation error and the cost of maintaining the factorized representation. Increasing cumulant order may therefore buy less accuracy per FLOP than at shallow depth. **I would benchmark the lower-order and augmented variants before assuming that the highest affordable order wins.**

There are three further limits:

1. **Fixed-width accuracy:** increasing the order is not known to give a sampling-competitive frontier as tolerance tends to zero at fixed width. ARC explicitly leaves this regime open. <citations>qShANRCA:412</citations>
2. **Training-induced structure:** trained networks can amplify selected higher-order dependencies that an initialization-based expansion regards as negligible. ARC identifies tracking those deviations as an essential unresolved extension. <citations>KvXnqCgQ:55-62</citations>
3. **FLOPs versus implementation:** the paper adjusts FLOP counts to remove redundant symmetric-tensor operations without implementing corresponding optimized kernels. The authors also report that their methods often lose in wall-clock time. Equal-FLOP superiority should not be read as current GPU superiority. <citations>qShANRCA:249-251,KvXnqCgQ:76-77</citations>

## 3. Localization: a plausible strictly better frontier

My hypothesis is **conditional, output-aware cumulant propagation with persistent shared coordinates**. The gain would come from representing structured dependence compactly rather than raising the global tensor order.

### Condition on shared directions, then integrate the residual analytically

For isotropic Gaussian input, choose orthonormal columns \(U\) and write
\[
X=Uz+\xi,\qquad z\sim\mathcal N(0,I),\quad
\xi\sim\mathcal N(0,I-UU^\top),
\]
with \(z\) and \(\xi\) independent. Then
\[
\mathbb E[f(X)]=\mathbb E_z[m(z)],\qquad
m(z)=\mathbb E_\xi[f(Uz+\xi)].
\]

Propagate conditional moments for the residual and integrate \(z\) using deterministic quadrature. Candidate directions should capture **shared gate changes and their effect on the requested output**, not merely explain input or activation variance. Preserve the same coordinates across layers; re-Gaussianizing them independently at each layer would discard the intended benefit.

Why could this help? Even if coordinates are nearly Gaussian conditional on a common mode, mixing over that mode produces higher-order dependence. For example, common conditional variance fluctuations generate fourth-order cross-cumulants. A low-dimensional latent representation can encode these collectively rather than storing a large tensor.

The qualification is important: conditioning improves the *representation opportunity*, not automatically the estimator. Bad directions, inaccurate conditional closure, or expensive quadrature can all make the frontier worse.

### Separate input radius—but do not confuse radial and angular error

There is an exact simplification for the bias-free ReLU architecture. Write \(X=RS\), with independent Gaussian radius \(R\) and uniform spherical direction \(S\). Positive homogeneity gives
\[
f(RS)=R f(S),\qquad
\mathbb E[f(X)]=\mathbb E[R]\,\mathbb E[f(S)].
\]

This is a mathematical identity, not a proposed approximation. Radius can therefore be integrated exactly before estimating the angular expectation. It prevents the propagation algorithm from repeatedly approximating the same radial fluctuation, but **does not solve the angular gate problem**. Biases break this simple factorization. A deeper-layer gain mode is potentially more valuable when it changes gate probabilities or carries angular dependence, rather than merely reproducing input radius.

### Reveal selected activation patterns—not all patterns

Conditioning on a complete gate pattern makes a ReLU network affine inside its activation region. Exact integration then becomes a Gaussian polyhedral-probability and truncated-moment problem. That shifts the difficulty rather than eliminating it: the number of regions and their constraints can be enormous.

A practical localization scheme would instead branch on gates or gate groups with large **output-relevant approximation error**, retain their conditioning constraints, and leave the remaining computation under a conditional cumulant closure. Stable gates can be treated as locally affine. Merely forcing a gate on or off without updating the conditional input distribution would not constitute correct conditioning.

### Enforce iterated-expectation consistency

For a branch event \(A\), compare the parent estimate with
\[
\widehat p_A\widehat\mu_A+
(1-\widehat p_A)\widehat\mu_{A^c}.
\]

A discrepancy is a useful refinement signal: the estimator can already predict that its own answer changes under a tractable decomposition. But agreement is **not an error certificate**; parent and children can share a bias.

This is an implementable analogue of ARC’s unpredictable-errors motivation, not a solved consequence of its heuristic-estimator theory. The iterated-expectations work identifies barriers and leaves the appropriate formalization open. Separately, ARC’s tail-risk discussion explicitly warns that layer-by-layer projection can lose shared computational dependencies and create correlated errors. <citations>PCUQgjwQ:21-36,riRJOfag:238-255</citations>

### What would establish strict improvement?

For deterministic quadrature weights \(w_j\), conditional closure error \(e(z_j)\), and quadrature error \(e_Q\), the resulting error decomposes as
\[
\widehat\mu-\mu=e_Q+\sum_j w_j e(z_j).
\]

A strictly better frontier requires conditioning to reduce the *integrated closure error* enough to repay direction discovery, conditional propagation, and quadrature. Signed-error cancellation can help, but cannot be assumed.

The decisive experiment would compare this localized scheme with the best factorized/augmented cumulant variants, full-covariance closure, and strong integration baselines on identical networks and budgets. It should include ablations for radius alone, shared directions, gain modes, selective gate branching, and consistency refinement. If localization only helps after expensive uncounted mode discovery, it has not improved the frontier.

**My strongest bet is persistent shared-mode conditioning plus selective trace/cumulant corrections—not a global Gaussian mixture or a blanket increase in cumulant order.** This directly targets the dependence that depth repeatedly transports and Gaussian projection repeatedly forgets. For randomly initialized networks it may improve constants or depth scaling; for trained networks it offers a way to capture structured deviations that the initialization expansion misses. Neither improvement is yet demonstrated for your architecture.

Finally, mean accuracy and tail accuracy must remain separate objectives. ARC reports promising low-probability estimates, but its earlier tail-risk experiments had not beaten importance sampling. Tail claims need comparison against importance sampling and event-specific validation, not ordinary Monte Carlo alone or accurate central moments. <citations>KvXnqCgQ:28-30,riRJOfag:68</citations>

[?] Evidence on cumulant propagation, analytic covariance propagation, deterministic Bayesian inference, finite-width theory, Gaussian mixtures, and ARC heuristic estimation was compared by estimand and benchmark applicability. Localization proposals were developed as mathematical hypotheses rather than reported experimental results.

[?] Synthesizing the cost–accuracy frontier and a conditional-localization hypothesis

[?] {"createdAt": "2026-10-09T13:26:02.997Z", "kind": "session_idle"}

[?] Reviewed excerpts from “Estimating the Expected Output of Wide Random MLPs More Efficiently Than Sampling” and staged a methods comparison; the delivered message began but was truncated.
