# Localization for deterministic approximation: results, error identities, and ReLU conjectures

# Localization for deterministic approximation: established results and Gaussian–ReLU proposals

## 1. What is established—and what is not

**There are localization-based deterministic approximation bounds, especially for mean-field free energy. There is not, in the literature located here, a general efficient deterministic integration algorithm obtained merely by replacing the posterior mean in stochastic localization with AMP/TAP.** Eldan’s entropy-efficient decomposition is a particularly direct positive answer: it converts control of the entropy lost during localization and the covariance remaining in the components into an explicit mean-field log-partition-function error bound. The algorithmic localization of El Alaoui–Montanari–Sellke, in contrast, is a randomized sampler with deterministic approximate-mean computations inside it. <citations>0,1,2,3</citations>

For the Gaussian-input network in the question, the cleanest proposal is **not to approximate the input posterior**: that posterior is already exactly Gaussian. Instead, evaluate the given moment-propagation rule on conditional Gaussian inputs, and use its defect in the backward heat equation as a signed bias correction. The exact error representation below is a derivation, not a published network-specific localization theorem. Efficient deterministic evaluation of that correction and a useful certificate remain conjectural.

### Scope and evidence qualification

This is a broad, non-exhaustive map of the retrieved literature, not a certified bibliography of every known use. “Localization” below means Eldan/Chen–Eldan measure-valued localization, or an explicitly identified equivalent Gaussian-channel/Föllmer/Polchinski construction. Ordinary analytic interpolation, spatial correlation decay, and polymer expansions are separated from that meaning. Some sources were accessible only through abstracts, and some PDF equations were fragmented by extraction; only recoverable formulas are stated quantitatively. The newer survey is Shi–Tian–Zhang, *Perspectives on Stochastic Localization*, arXiv:2510.04460; the trickle-down paper is Anari–Koehler–Vuong, arXiv:2407.16104. <citations>4,5</citations>

The network dimensions and the approximately thousand-fold accuracy advantage are supplied premises, not verified empirical findings. All network formulas below concern a **scalar output**, or one fixed linear functional of a vector output. Weights, biases, output normalization, the particular closure, the accuracy metric, and the reference used to establish “truth” have not been supplied.

## 2. Published uses and adjacent deterministic methods

### 2.1 Direct positive result: localization gives mean-field free-energy error bounds

**Eldan, *Taming correlations through entropy-efficient measure decompositions with applications to mean-field approximation*, arXiv:1811.11530.** For a measure \(\mu\) with finite covariance and well-defined entropy relative to a background measure, Theorem 2 constructs a mixture \(\mu=\mathbb E_\theta\mu_\theta\), for any positive-definite matrix \(L\), satisfying the following inequalities. <citations>6</citations>

$$
H(\mu)-\mathbb E H(\mu_\theta)
\le \log\det(I+\operatorname{Cov}(\mu)L)
\le \operatorname{Tr}(\operatorname{Cov}(\mu)L),
$$

$$
\mathbb E\operatorname{Cov}(\mu_\theta)\preceq L^{-1},\qquad
\mathbb E[\operatorname{Cov}(\mu_\theta)L\operatorname{Cov}(\mu_\theta)]
\preceq\operatorname{Cov}(\mu).
$$

These are directional entropy–covariance trade-offs, not a construction of a cheaply enumerable mixture. <citations>7</citations>

For the product-reference quadratic Gibbs model

$$
f(x)=x^\top Jx+h^\top x,\quad
\mu(dx)=Z^{-1}e^{f(x)}\nu^{\otimes n}(dx),\quad
\mathcal F_{\rm MF}=\sup_{\xi\text{ product}}\{\mathbb E_\xi f+H_{\nu^{\otimes n}}(\xi)\},
$$

Theorem 5 gives the explicit deterministic inequality

$$
\boxed{0\le\log Z-\mathcal F_{\rm MF}
\le3\log\det(I+\operatorname{Cov}(\mu)|J|)},\qquad |J|=(J^2)^{1/2}.
$$

The convention here has **no factor one-half in the quadratic Hamiltonian**. This is a genuine approximation-error theorem obtained from localization decomposition. It is not a guarantee that global optimization of the mean-field objective is easy. The covariance on the right is also generally unknown. <citations>8,9,1</citations>

If each single-spin support has diameter at most \(D\), a computable corollary is

$$
\log Z-\mathcal F_{\rm MF}
\le3\operatorname{rank}(J)\log\bigl(1+D^2n\|J\|_{\rm op}\bigr).
$$

Thus low rank can give a logarithmic rather than polynomial dimension dependence. This result covers the stated discrete and continuous bounded-spin product-reference settings; it is not a theorem about a deep network’s hidden activations. <citations>10,11</citations>

**Eldan, *Gaussian-width gradient complexity, reverse log-Sobolev inequalities and nonlinear large deviations*, arXiv:1612.04346.** For the uniform reference \(u\) on \(\{-1,1\}^n\), Corollary 2 gives the following bound. <citations>12</citations>

$$
0\le \log\mathbb E_u e^f
-\sup_{\xi\text{ product}}\{\mathbb E_\xi f-D(\xi\|u)\}
\le64\operatorname{Lip}(f)^{2/3}\mathcal D(f)^{1/3}n^{2/3},
$$

where \(\partial_i f\) is half the difference between the two coordinate-flipped values, \(\operatorname{Lip}(f)=\max_{i,x}|\partial_i f(x)|\), and \(\mathcal D(f)\) is the Gaussian width of the discrete-gradient image augmented by zero. The proof framework uses entropy-optimal stochastic control. This is another direct deterministic free-energy approximation bound, not an FPTAS for arbitrary \(f\). <citations>13,12,14</citations>

**Related structural results.** Eldan–Gross, *Decomposition of mean-field Gibbs distributions into product measures*, places most mixture mass near solutions of \(m\approx\tanh(\nabla f(m))\), under low gradient complexity and suitable Lipschitz assumptions. This connects localization to mean-field stationary points but does not supply deterministic mixture weights or a general error bound for every observable. El Alaoui–Montanari’s *An Information-Theoretic View of Stochastic Localization*, arXiv:2109.00709, gives an elementary Gaussian-channel proof of an entropy–covariance decomposition. <citations>15,16,17</citations> For \(Y=\sqrt\tau X+Q^{1/2}Z\), \(\tau\sim\mathrm{Unif}[1,2]\), \(Q\succ0\), it proves

$$
\mathbb E\Sigma_\theta\preceq Q,\quad
I(X;\theta)\le\tfrac12\log\det(I+2Q^{-1}\Sigma),\quad
\mathbb E[\Sigma_\theta Q^{-1}\Sigma_\theta]\preceq\Sigma.
$$

Neither a small information budget nor small component covariances automatically makes the outer mixture integral computationally tractable. <citations>15,16,17</citations>

### 2.2 AMP/TAP posterior means: deterministic inner oracle, randomized outer algorithm

**El Alaoui–Montanari–Sellke, *Sampling from the Sherrington–Kirkpatrick Gibbs measure via algorithmic stochastic localization*.** The retrieved FOCS version proves: for any fixed \(\varepsilon>0\), \(\beta_0<1/2\), algorithm parameters can be chosen independently of \(N\) so that, for \(\beta\le\beta_0\), with probability \(1-o_N(1)\) over the GOE disorder, the following guarantee holds. <citations>3</citations>

$$
W_{2,N}(\mu_A^{\rm alg},\mu_A)\le\varepsilon,
\qquad \text{runtime }O(N^2),
$$

where \(W_{2,N}^2=\inf\mathbb E\|X-Y\|^2/N\). The algorithm explicitly draws Gaussian increments and performs randomized final rounding. <citations>18,2,3</citations>

Proposition V.6 gives asymptotic accuracy of the AMP posterior-mean oracle along the localization path for \(\beta<1\); the sampler’s stated theorem has the stricter \(\beta<1/2\) condition. These are distinct claims. The method reduces **sampling to approximate mean computation**, not mean computation to deterministic integration of the path. <citations>19,20,21</citations>

The random-linear-model analogue, *Sampling from the Random Linear Model via Stochastic Localization Up to the AMP Threshold*, likewise states convergence of a sample distribution in smoothed KL below its specified AMP noise threshold. <citations>22</citations>

**A published TAP-oracle correction along localization:** *Sampling from Spherical Spin Glasses in Total Variation via Algorithmic Stochastic Localization* explicitly constructs a correction to the AMP-selected TAP fixed point to improve the tilted-mean estimator. It obtains a polynomial-time sampler with vanishing total-variation error for mixtures satisfying \(\xi''(s)<(1-s)^{-2}\) for all \(s\in[0,1)\). This is directly relevant to correcting an approximate moment oracle, but the final algorithm remains a sampler; the retrieved abstract does not provide the correction's formula. It should not be confused with a covariance-conservation-based certificate for arbitrary observables. <citations>23</citations>

**A genuinely deterministic TAP result, but not a localization quadrature theorem:** *Mean-field variational inference with the TAP free energy: Geometric and statistical properties in linear models* proves asymptotic consistency of certain TAP minima for posterior moments and evidence in Gaussian-design Bayesian linear regression. Under its compactly supported product prior, proportional-dimensional regime, and unique nondegenerate replica-symmetric minimizer assumptions, Theorem 3.2 gives coordinatewise mean-square convergence for Lipschitz posterior test functions. Local optimization uses AMP initialization followed by natural-gradient descent in the computationally easy regime. These assumptions do not transfer to arbitrary fixed ReLU weights. <citations>24,25,26,27</citations>

### 2.3 Heat flow, Polchinski evolution, Föllmer drift, and Schrödinger bridges

Shi–Tian–Zhang explicitly identifies stochastic localization with a Gaussian posterior channel and, after time changes, diffusion, Polchinski semigroup, and Schrödinger-bridge constructions. Its Polchinski perspective supplies exact expectation evolution and backward/forward semigroup equations. These furnish a natural starting point for Duhamel error representations; an exact representation is not itself an efficient deterministic algorithm. <citations>28,29,30</citations>

The dynamic bridge perspective supplies the drift-energy identity

$$
D(P\|W)=\tfrac12\mathbb E_P\int_0^1\|b_s\|^2ds
$$

under the usual absolute-continuity/integrability conditions, with data-processing bounds for marginals. Föllmer drift is the point-start Wiener-reference special case. Deterministically solving the associated heat, Hamilton–Jacobi, or Fokker–Planck equations would evaluate expectations, but a general high-dimensional solver with a localization-derived complexity guarantee was not located. <citations>31</citations>

Adjacent approximation results include *Weak approximation of Schrödinger–Föllmer diffusion*: the retrieved abstract states weak convergence of time-discretized terminal distributions under mild regularity, but not an explicit rate usable here. Föllmer-process stability results for log-Sobolev, Talagrand, and Shannon–Stam inequalities give distributional approximation information under substantial covariance/log-concavity/spectral-gap hypotheses; they are not neural-network closure certificates. <citations>32,33,34</citations>

**General Duhamel weak-error formula, derived consequence.** Suppose the true and approximate diffusion generators are \(L_s\) and \(\widehat L_s\), with the same initial law, and \(\widehat v\) solves \(\partial_s\widehat v+\widehat L_s\widehat v=0\), \(\widehat v(T)=F\). Under conditions justifying Itô and integration,

$$
\mathbb E_P F(X_T)-\mathbb E_{\widehat P}F(\widehat X_T)
=\int_0^T\mathbb E_P[(L_s-\widehat L_s)\widehat v(s,X_s)]ds.
$$

For equal diffusion matrices, the integrand is \((b_s-\widehat b_s)\cdot\nabla\widehat v\); differing diffusion covariances add \(\tfrac12\operatorname{Tr}[(a_s-\widehat a_s)\nabla^2\widehat v]\). This is the precise error representation behind a drift-oracle correction. Evaluating its true-law expectation and the approximate backward solution still requires work; declaring the oracle deterministic does not make those tasks tractable. The survey's semigroup equations supply the corresponding localization evolution framework. <citations>29</citations>

### 2.4 Barvinok, correlation decay, and cluster expansions: deterministic, but not Eldan localization

These are relevant analogies and possible computational tools, **not established equivalents of a random measure-valued localization scheme** in the sources retrieved.

**Barvinok/Taylor interpolation.** A zero-free complex neighborhood permits approximation of \(\log Z(z)\) by a truncated Taylor series, provided coefficients can be computed efficiently. As a concrete theorem, Patel–Regts, *Deterministic Polynomial-Time Approximation Algorithms for Partition Functions and Graph Polynomials*, Theorem 1.4 gives a deterministic multiplicative \(\varepsilon\)-approximation in \((|V|/\varepsilon)^{O(1)}\) for fixed maximum degree \(\Delta\), fixed spin alphabet size, and interaction entries satisfying \(|A_{ij}-1|\le0.34/\Delta\). Multiplicative accuracy is defined using log-magnitude and phase errors. <citations>35,36,37</citations>

**Correlation decay.** The overlap is real: *Barvinok’s interpolation method meets Weitz’s correlation decay approach* combines these methods and gives deterministic approximate counting of sink-free orientations on graphs of minimum degree at least three, with factor \(e^\varepsilon\) and runtime \(O(n(m/\varepsilon)^7)\). It is not a trickle-down error-correction theorem. In a correlation-decay computation, the actual approximation mechanism is truncating a conditional recursion with a proven boundary-influence tail—not integrating a Brownian localization trajectory. <citations>38,39</citations>

**Cluster/contour expansions.** Helmuth–Perkins–Regts, *Algorithmic Pirogov–Sinai theory*, combines contour representations with truncated Taylor counting. A stated consequence is an FPTAS for the hard-core partition function at sufficiently high fugacity on lattice subsets with suitable boundary conditions. This illustrates that deterministic expansion methods are not restricted to a naive high-temperature expansion; they can operate around a controlled phase. No conversion of a localization covariance identity into such an expansion for a dense deep ReLU network was located. <citations>40</citations>

**Marginals from partition functions require their own error analysis.** For a source parameter \(h\), \(\partial_h\log Z\) is an expectation. A uniform additive free-energy error \(\epsilon\) does not justify differentiating an arbitrary approximation. Derived consequence: if \(|G''|\le M\) on \([-a,a]\) and the approximate values at \(\pm a\) each have error at most \(\epsilon\), then the central-difference marginal estimate has error at most \(\epsilon/a+Ma/2\). Thus free-energy accuracy alone generally loses precision when converted into moments. The partition-function derivative relation is explicitly noted by Eldan. <citations>41</citations>

## 3. What trickle-down and conservation can—and cannot—certify

### 3.1 Exact identity and the covariance convention

For the discrete affine update

$$
\frac{d\nu_{k+1}}{d\nu_k}(x)=1+\langle x-a_k,Z_{k+1}\rangle,
\quad \mathbb E[Z_{k+1}\mid\mathcal F_k]=0,
$$

assume the density is nonnegative. With \(\Sigma_k=\operatorname{Cov}(\nu_k)\) and \(C_k=\mathbb E[Z_{k+1}Z_{k+1}^\top\mid\mathcal F_k]\), Anari–Koehler–Vuong’s equation (8) is

$$
\boxed{\Sigma_k-\mathbb E[\Sigma_{k+1}\mid\mathcal F_k]=\Sigma_kC_k\Sigma_k.}
$$

This exact discrete identity concerns **affine density tilts**, not an arbitrary finite exponential-tilt step. Nonzero affine tilts cannot stay nonnegative on an everywhere-positive Gaussian density with unbounded support; the Gaussian setting instead uses the continuous process. <citations>42,43</citations>

For

$$
d\nu_t(x)=\nu_t(x)\langle x-a_t,B_t\,dW_t\rangle,
\qquad Q_t=B_tB_t^\top,
$$

the corresponding continuous identity is

$$
\Sigma_0=\mathbb E\Sigma_T+
\int_0^T\mathbb E[\Sigma_tQ_t\Sigma_t]dt.
$$

The paper uses a symmetric driving matrix and writes its square: confusing that amplitude with the increment covariance causes a missing-square error. <citations>44</citations>

For \(g_t=\mathbb E_{\nu_t}F\), let \(v_t=\operatorname{Cov}_{\nu_t}(X,F)\). Derived from the measure-valued SDE,

$$
dg_t=v_t^\top B_t\,dW_t,\qquad
\operatorname{Var}_{\nu_0}(F)
=\mathbb E\operatorname{Var}_{\nu_T}(F)+
\int_0^T\mathbb E[v_t^\top Q_tv_t]dt.
$$

This gives the variance of a stopped conditional-expectation estimator. It does **not** identify the bias of a deterministic closure for \(g_t\).

### 3.2 Published covariance bounds can enter a free-energy certificate

Anari–Koehler–Vuong Theorem 3: for an Ising interaction \(J\succeq0\) with constant diagonal \(\alpha\), let \(\eta=\alpha/\|J\|_{\rm op}\). Then

$$
\|\operatorname{Cov}(\mu)\|_{\rm op}\le q_\eta(\|J\|_{\rm op}),
\quad q_\eta(z)=r(\eta z)+\int_0^zq_\eta(y)^2dy,
$$

$$
r(t)=\mathbb E_{\varepsilon=\pm1,G\sim N(0,1)}
[1-\tanh^2(t\varepsilon+\sqrt tG)].
$$

The approximate tensorization-of-entropy constant is at most \(\exp(\int_0^{\|J\|_{\rm op}}q_\eta(z)dz)\), when the bound is finite. <citations>45</citations>

**Derived cross-paper consequence, not claimed as a theorem of either paper:** if a compatible covariance estimate \(\operatorname{Cov}(\mu)\preceq K\) is available, Eldan’s theorem gives

$$
\log Z-\mathcal F_{\rm MF}
\le3\log\det(I+|J|^{1/2}K|J|^{1/2}).
$$

Taking \(K=\kappa I\) gives \(3\log\det(I+\kappa|J|)\). To use the trickle-down theorem for \(\kappa\), reconcile the two papers’ Hamiltonian conventions and any diagonal shifts first. Diagonal shifts leave the discrete Ising law unchanged but alter the free-energy normalization and spectral bound. This supplies a genuine way for trickle-down covariance control to strengthen a mean-field certificate. <citations>1,46</citations>

### 3.3 Entropy is useful only when it controls the approximation actually used

The entropy-efficient mean-field result above is a positive answer to the question about conservation producing error bounds. Chen–Eldan-style “entropy conservation” also supports modified log-Sobolev and mixing estimates; that use is not automatically an approximation-error result. The survey’s appendix explicitly follows this route. <citations>47</citations>

Here are two useful **derived** bridges to expectation error:

* If \(P,Q\) are diffusion path laws with identity diffusion coefficient, the same initial law, and appropriate Girsanov conditions, then \(D(P\|Q)=\tfrac12\mathbb E_P\int\|b-\widehat b\|^2dt\). Endpoint KL is no larger. For bounded \(F\), this gives \(|\mathbb E_PF-\mathbb E_QF|\le\operatorname{osc}(F)\sqrt{D(P\|Q)/2}\). These require an integrated drift-error estimate; a covariance conservation residual alone is not enough. The underlying drift-energy principle is in the survey. <citations>48</citations>
* For a Gaussian reference \(G=N(m,S)\), \(S\succ0\), and a Euclidean \(L_F\)-Lipschitz function, the Gaussian transport-entropy inequality implies \(|\mathbb E_PF-\mathbb E_GF|\le L_F\sqrt{2\|S\|_{\rm op}D(P\|G)}\). This is a useful conditional Gaussian-closure certificate only if KL can itself be bounded. Singular ReLU hidden-state laws can have infinite KL to a nonsingular Gaussian, making this particular route vacuous. Föllmer stability results should not be interpreted as removing that obstruction. <citations>49,50</citations>

A simple mathematical counterexample explains the limitation: a Rademacher variable and a standard Gaussian have identical mean and variance, but their expected positive parts are respectively \(1/2\) and \(1/\sqrt{2\pi}\). Consequently, matching means/covariances—even exactly—cannot certify a general ReLU expectation. This is a direct calculation, not an empirical comparison.

## 4. Exact localization error representation for the specified network

All formulas in this section are **derived here** from Gaussian conditioning and Itô’s formula; they are not attributed to a network-specific paper. Their localization foundation is the Gaussian-channel representation in the survey. <citations>30</citations>

### 4.1 Input localization is explicit

For \(X\sim N(0,I_d)\), isotropic localization has

$$
\nu_t=N\left(m_t,(1+t)^{-1}I_d\right),\quad
m_t=\frac{c_t}{1+t},\quad dm_t=(1+t)^{-1}dW_t.
$$

Set \(s=t/(1+t)\). In distribution, the mean process becomes standard Brownian motion \(B_s\), and

$$
X\mid\mathcal F_s\sim N(B_s,(1-s)I_d),\qquad X=B_1.
$$

Thus

$$
u(s,m)=\mathbb E[F(m+\sqrt{1-s}Z)],\quad
\partial_su+\tfrac12\Delta_mu=0,\quad u(1,m)=F(m).
$$

For the network in the question, take \(d=1024\); width and depth enter through \(F\), not through a new posterior law.

### 4.2 Signed correction for any conditional moment-propagation rule

Let \(\widehat u(s,m)\) be the existing deterministic moment-propagation estimate applied to input \(N(m,(1-s)I)\), with \(A=\widehat u(0,0)\). Assume sufficient regularity and integrability; otherwise work on \([0,T]\), \(T<1\), and retain the endpoint term. Define

$$
R(s,m)=\partial_s\widehat u(s,m)+\tfrac12\Delta_m\widehat u(s,m).
$$

If \(\widehat u(1,m)=F(m)\) and the residual is integrable, Itô gives the exact identity

$$
\boxed{\mathbb E F(X)-A
=\int_0^1\mathbb E_{M\sim N(0,sI)}R(s,M)\,ds.}
$$

This is a **signed weak-error/Duhamel correction**, not a variance bound. A single random localization path would estimate it stochastically; evaluating the outer Gaussian integrals by certified deterministic methods would make the correction deterministic.

Without terminal exactness, the full identity is

$$
\mathbb EF-A
=\mathbb E[F(B_1)-\widehat u(T,B_T)]
+\int_0^T\mathbb E R(s,B_s)ds.
$$

A deterministic approximation \(C_T\) to the residual integral, with error at most \(\epsilon_Q\), yields

$$
|\mathbb EF-(A+C_T)|\le\epsilon_Q+\epsilon_T,
\quad \epsilon_T\ge\mathbb E|u(T,B_T)-\widehat u(T,B_T)|.
$$

These bounds are exact but only become informative when \(\epsilon_Q,\epsilon_T\) are computable. The integration problem has been relocated to a potentially simpler residual, not eliminated.

### 4.3 ReLU kinks cannot be discarded

A ReLU network is continuous and piecewise affine. Its classical Hessian vanishes inside activation cells, but its distributional curvature lives on cell boundaries. Therefore “the Hessian is zero almost everywhere, so the correction vanishes” is wrong. One must differentiate the conditional Gaussian-smoothed estimator for \(s<1\), account for boundary terms, or retain a terminal remainder. Neither width nor depth alone ensures that the residual is small, low rank, or cheaply integrable.

For a finite network, a valid global Lipschitz bound is the product of layer operator norms, including the output map. It may be much too large for a useful certificate. A finite-norm Gaussian-input network has finite first and second output moments, so the Brownian and variance identities are not obstructed by unbounded ReLU outputs.

## 5. Concrete conjectures and conditional certificates

The following are **proposed research directions**, not existing guarantees. Each states the unresolved assumption rather than hiding it in “localization.”

### Conjecture A: low-dimensional residual correction

There exists a modest rank \(r\ll d\), an orthonormal \(U\in\mathbb R^{d\times r}\), and a computable residual model \(\widetilde R(s,U^\top m)\) such that

$$
\int_0^T\mathbb E|R(s,B_s)-\widetilde R(s,U^\top B_s)|ds\le\epsilon_{\rm proj}.
$$

If deterministic quadrature of the resulting \((r+1)\)-dimensional integral has error \(\epsilon_Q\), and the terminal defect is bounded by \(\epsilon_T\), then **the following consequence is a theorem**, not a conjecture:

$$
|\mathbb EF-(A+C)|\le
\epsilon_{\rm proj}+\epsilon_Q+\epsilon_T.
$$

The conjecture is the existence and economical certification of this low-rank residual for the particular network. Candidate directions should come from the residual’s sensitivity or curvature, not just from the output-gradient covariance: small gradient energy does not by itself certify the residual’s approximation error.

### Conjecture B: a signed second-order closure correction dominates higher orders

Build a hierarchy \(\widehat u^{(1)},\widehat u^{(2)},\ldots\), where the first rule is the existing closure and the next includes selected pairwise gate/preactivation information. Let \(R^{(j)}=(\partial_s+\Delta/2)\widehat u^{(j)}\). The testable conjecture is that, for this network and a specified cost budget,

$$
\int_0^T\mathbb E|R^{(2)}(s,B_s)|ds+\epsilon_T^{(2)}
\ll
\int_0^T\mathbb E|R^{(1)}(s,B_s)|ds+\epsilon_T^{(1)}.
$$

A sharper version asks that the signed first residual correction predicts the actual baseline bias, while the next residual supplies a much smaller remainder. The covariance-budget identity is a useful consistency check for each closure, but is not the missing bound on higher-order gate correlations. Nor does the sign or form of an Onsager/TAP correction transfer automatically from a spin glass to this network.

### Conjecture C: terminal activation-cell certificates make truncated localization practical

For a mean \(m\) away from activation boundaries, the current activation pattern defines a polyhedral cell and an affine extension \(\ell_m\) with \(\ell_m(m)=F(m)\). Let \(\sigma=\sqrt{1-T}\), and let \(p_{\rm exit}(m)\) bound the probability that \(m+\sigma Z\) exits this cell. If \(L_F\) bounds the Lipschitz norm of both the network and its cellwise affine pieces, then the following **derived certificate** holds:

$$
|u(T,m)-F(m)|
\le2L_F\sigma\sqrt{d\,p_{\rm exit}(m)}.
$$

Reason: the difference from the affine extension is zero inside the cell and at most \(2L_F\sigma\|Z\|\) outside; the affine extension has Gaussian mean \(F(m)\), and Cauchy–Schwarz completes the bound.

A computable union bound uses the cell’s defining affine gate inequalities. If a gate’s signed margin at \(m\) is \(b_j(m)>0\), with normal \(a_j\), then

$$
p_{\rm exit}(m)\le
\min\left\{1,\sum_j\Phi\left(-\frac{b_j(m)}{\sigma\|a_j\|}\right)\right\}.
$$

These normals must be obtained from the **complete fixed-pattern input-space cell**, including downstream gates—not by treating hidden preactivations as independent Gaussians. Degenerate/zero-margin gates need separate handling. For the actual terminal approximation,

$$
|u(T,m)-\widehat u(T,m)|
\le |F(m)-\widehat u(T,m)|+2L_F\sigma\sqrt{d\,p_{\rm exit}(m)}.
$$

The conjecture is that adaptive deterministic partitioning of high-probability mean space makes the integral of this certificate small at acceptable cost. Many near-zero gate margins or a loose product-of-norms bound may defeat it.

### Conjecture D: certify the residual in a weighted norm, not by path conservation alone

If a deterministic procedure proves

$$
\int_0^T\mathbb E_{N(0,sI)}|R(s,M)|ds\le\epsilon_R,
$$

and bounds the terminal defect by \(\epsilon_T\), then \(|\mathbb EF-A|\le\epsilon_R+\epsilon_T\). This is already a mathematical certificate. The conjecture is that Gaussian-weighted interval bounds, low-rank representations, or sparse expansions of **the residual** can prove a bound close to the observed small error without prohibitive dimensional dependence. Unweighted suprema over all input space and global Lipschitz constants are likely to be too pessimistic. Signed corrections can be much smaller than this absolute-residual bound; accuracy and certification are separate goals.

## 6. How to interpret and test the reported accuracy advantage

Treat the reported advantage as an observation to validate, not a consequence of localization. If “a thousand times more accurate” means a thousand-fold smaller RMSE than independent Monte Carlo at matched cost, then at that cost the baseline error would be roughly \(\operatorname{sd}(F)/(1000\sqrt N)\), where \(N\) is the number of complete forward evaluations the Monte Carlo budget purchases. This is a conditional scaling calculation, not a performance measurement.

Before testing the conjectures:

1. Specify absolute versus relative error, scalar versus vector norm, and fixed-instance versus across-network error. Establish the reference independently of the approximation being judged; a Monte Carlo reference whose uncertainty exceeds the deterministic error cannot resolve the claimed advantage.
2. Include the full deterministic cost: covariance propagation, derivatives of the closure, trace/Hessian calculations, residual integration, and certification. Randomized trace estimation or randomized cubature makes the resulting procedure hybrid, not fully deterministic.
3. Test affine networks and a single affine–ReLU layer, where the conditional Gaussian expectation is explicit. In these cases an exact conditional oracle must have zero heat residual. Then test smaller deep networks with independently resolved expectations.
4. Measure \(R\), its signed integral, its absolute integral, and terminal defects separately. A closure can be very accurate by cancellation while giving a weak absolute-error certificate.
5. Compare against deterministic alternatives and structured randomized integration, not only raw Monte Carlo. A large raw-Monte-Carlo gap is not itself evidence for a localization mechanism.

**Recommendation:** retain the existing estimate as the baseline; first derive and evaluate its conditional-Gaussian heat residual. Investigate low-rank residual quadrature and activation-cell terminal bounds before building an AMP/Föllmer sampler. The published entropy-efficient mean-field theorems establish that localization can yield approximation certificates, but a certificate for this network needs additional control of the observable-specific residual or gate-boundary error—not just covariance conservation.
