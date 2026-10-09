# Stochastic localization, exact RG, and defect-derived counterterms

# Stochastic localization, exact RG, and defect-derived counterterms

The localization–Polchinski–heat-semigroup connection is an exact identity after specifying the Gaussian reference measure, covariance schedule, and direction of time. The Föllmer process is the associated heat-semigroup Doob transform. Shi–Tian–Zhang give an explicit particle and posterior identification; Bauerschmidt–Bodineau–Dagallier develop the covariance-dependent version and use the resulting flow for functional inequalities. <citations>0,1,2,3,4,5,6</citations>

A truncation has a residual in the full equation, and that residual yields exact error identities through Duhamel’s formula or Itô’s formula. These identities do **not** by themselves make the error computable: the residual, the propagator, or its expectation may require precisely the information the closure discarded. Hopf algebras provide established composition, subtraction, and modified-equation machinery, but the retrieved literature does not establish a general Connes–Kreimer counterterm prescription for finite cumulant closures of deep-network expectations.

**Scope and status.** This is a primary-source map of the principal precise links, not an exhaustive theorem census. Sections 1 and 3 describe published structures; Section 2 gives explicit mathematical derivations and connects them to published approximation methods; Section 4 proposes a new estimator design. Formulas called “derived” are deductions here, not results attributed to a named paper. Finite-dimensional formulas avoid the additional regularization and domain issues of infinite-dimensional field theory.

## 1. The exact dictionary

### 1.1 Gaussian convolution and the Polchinski potential equation

Let \(C_s\) be a differentiable, increasing positive-semidefinite covariance, with \(C_0=0\), and write

$$
L_s=\tfrac12\dot C_s:\nabla^2,\qquad
Z_s(x)=\mathbb E[e^{-V_0(x+\xi_s)}],\quad \xi_s\sim N(0,C_s),\qquad V_s=-\log Z_s.
$$

Subject to integrability and differentiation under the integral,

$$
\partial_sZ_s=L_s Z_s,\qquad
\partial_sV_s=\tfrac12\dot C_s:\nabla^2V_s
-\tfrac12\nabla V_s^\top\dot C_s\nabla V_s.
$$

The first equation is linear heat evolution; the second is its logarithmic, nonlinear Polchinski form. This is the Gaussian-convolution construction in Bauerschmidt–Bodineau–Dagallier, Section 3.2. It is important not to call the potential itself a solution of the linear heat equation. <citations>4</citations>

In Shi–Tian–Zhang’s increasing particle-time convention, the *remaining* covariance decreases:

$$
Z_\tau=P_{1-\tau}e^{-V_1},\qquad
\partial_\tau Z_\tau=-\tfrac12\Delta Z_\tau,\qquad
\partial_\tau V_\tau=-\tfrac12\Delta V_\tau+\tfrac12|\nabla V_\tau|^2.
$$

These are their equations (16)–(17). The sign difference is a time-direction difference, not a different RG equation. <citations>2,3</citations>

“Wilson exact RG” is broader than this particular representation. Cutoff actions, derivative expansions, Legendre-transformed effective actions, and rescaled dimensionless flows need not have literally this PDE: transformations can introduce scaling terms or change the functional equation. Aoki et al. explicitly distinguish Wilsonian and cutoff 1PI effective actions and discuss the Legendre relation to Polchinski’s equation. <citations>7,8,9</citations>

### 1.2 The Polchinski Markov semigroup is a Doob transform

For \(h_\tau=P_{1-\tau}f>0\), define, for \(\sigma\leq\tau\),

$$
Q_{\sigma,\tau}g(x)
=\frac{P_{\tau-\sigma}(h_\tau g)(x)}{h_\sigma(x)}.
$$

**Derived directly from heat composition:** \(Q_{\sigma,\tau}1=1\),
\(Q_{\sigma,\tau}Q_{\tau,\rho}=Q_{\sigma,\rho}\), and its generator is

$$
\mathcal A_\tau g=\tfrac12\Delta g+\nabla\log h_\tau\cdot\nabla g
=\tfrac12\Delta g-\nabla V_\tau\cdot\nabla g.
$$

This is Shi–Tian–Zhang’s Polchinski generator (14); their Kolmogorov identities are equation (15). Bauerschmidt–Bodineau–Dagallier give the covariance-dependent two-parameter semigroup in Proposition 3.5. A two-parameter inhomogeneous Markov semigroup is not automatically a one-parameter subgroup. <citations>1,10,11</citations>

### 1.3 Gaussian-channel stochastic localization

For a prior \(\mu\), let \(X\sim\mu\) and observe

$$
Y_t=tX+B_t.
$$

The posterior is

$$
\mu_{t,y}(dx)=\frac{e^{y\cdot x-t|x|^2/2}\mu(dx)}{Z(t,y)},\qquad
Z(t,y)=\int e^{y\cdot x-t|x|^2/2}\mu(dx).
$$

In the observation filtration,

$$
dY_t=m_tdt+dW_t,\qquad m_t=\mathbb E[X\mid Y_t].
$$

This is the Bayesian channel/innovation description in Shi–Tian–Zhang, Section 3. Here \(W\) is the innovation Brownian motion, not necessarily the original noise \(B\). <citations>12,13,14</citations>

**Derived posterior identities:** for an integrable test function \(G\),

$$
d\mu_t(G)=\operatorname{Cov}_{\mu_t}(G,X)\cdot dW_t,
\qquad dm_t=\operatorname{Cov}_{\mu_t}(X)\,dW_t.
$$

Thus conditional expectations are martingales, while the changing posterior carries the information. Moreover,

$$
\partial_tZ=-\tfrac12\Delta_yZ,\qquad
\partial_tK=-\tfrac12(\Delta_yK+|\nabla_yK|^2),\quad K=\log Z.
$$

The negative logarithm \(-K\) therefore obeys the same backward-covariance Polchinski form. These equations follow simply by differentiating the posterior normalizer; they require finite tilted moments on the relevant domain.

For a non-isotropic observation \(dY=A_tXdt+A_t^{1/2}dB\), the precision increment is \(A_tdt\), the posterior tilt replaces \(tI\) by \(\int_0^tA_sds\), and the same formulas use the corresponding covariance contractions. Adapted schedules require their predictability and a well-posed observation model.

### 1.4 The explicit localization–Polchinski particle identification

Shi–Tian–Zhang define particles \(v_\tau\) and fluctuation measures by

$$
dv_\tau=\frac{m_\tau-v_\tau}{1-\tau}d\tau+dW_\tau,
\qquad
\pi_\tau^v(dx)\propto
\exp\!\left(\frac{v\cdot x}{1-\tau}-\frac{\tau|x|^2}{2(1-\tau)}\right)\pi_0(dx).
$$

Their Theorem 4 identifies these with localization by

$$
t=\frac{\tau}{1-\tau},\qquad c_t=\frac{v_\tau}{1-\tau},\qquad
\pi_t=\pi_\tau^{v_\tau}.
$$

This identifies posterior measures and processes under deterministic time/space change—not merely an analogy between PDEs. The Brownian motions must also be transformed under the time change. <citations>15,16,0</citations>

### 1.5 Föllmer drift and its terminal-law interpretation

Let the target \(\mu\) have density \(f=d\mu/d\gamma\) relative to \(\gamma=N(0,I)\). With appropriate existence and entropy assumptions,

$$
dU_\tau=\nabla\log P_{1-\tau}f(U_\tau)d\tau+dW_\tau,
\quad U_0=0,\quad U_1\sim\mu.
$$

This is the Föllmer process. Its energy-minimizing drift satisfies

$$
H(\mu\mid\gamma)=\tfrac12\mathbb E\int_0^1
|\nabla\log P_{1-\tau}f(U_\tau)|^2d\tau.
$$

The Brownian transport map paper gives the heat-semigroup drift, minimum-energy property, and entropy identity in Section 2, and identifies the conditional terminal law with stochastic localization in Lemma 4.1. Bauerschmidt–Bodineau–Dagallier explain the same identification with their covariance/time conventions in Section 5. <citations>17,18,19,20,21,6</citations>

**Derived bridge formula:** Brownian motion conditioned on terminal value \(X\) has \(U_\tau\mid X\sim N(\tau X,\tau(1-\tau)I)\). Consequently

$$
\operatorname{Law}(X\mid U_\tau=u)\propto
\exp\!\left(\frac{u\cdot X}{1-\tau}-\frac{\tau|X|^2}{2(1-\tau)}\right)\mu(dX),
\quad
\nabla\log P_{1-\tau}f(u)=\frac{\mathbb E[X\mid U_\tau=u]-u}{1-\tau}.
$$

This gives the same particle SDE as above. The Gaussian-reference density ratio matters: using a Lebesgue density in place of \(d\mu/d\gamma\) without compensation changes the drift.

### 1.6 Why the functional-inequality results matter—and what they do not give

Bauerschmidt–Bodineau–Dagallier’s multiscale Bakry–Émery criterion controls the log-Sobolev constant through integrated Hessian information along the exact renormalized potentials. Their account emphasizes that sufficient improvement of nonconvexity along the flow yields the inequality, and that the proof iterates entropy decomposition across scales. <citations>22,23</citations>

This supplies a possible stability tool for approximation analysis, but it is **not** already a theorem that a small finite-coupling RG residual gives a small observable error. Such a theorem needs a specified norm, endpoint conditions, and stability of the relevant propagator.

## 2. Truncation defects: exact identities, bounds, and corrections

### 2.1 What a closure omits

Write the full evolution as \(\partial_sV=\mathcal R_s(V)\). A finite ansatz \(\widetilde V(g_1,\ldots,g_M)\) gives a residual

$$
d_s=\partial_s\widetilde V-\mathcal R_s(\widetilde V).
$$

A projected flow may enforce \(\Pi d_s=0\), while \((I-\Pi)d_s\neq0\). This is a **definition**, not a small-error claim. If the manifold is invariant under the exact equation, a finite ansatz need not have any residual at all.

For example, with scalar covariance speed \(q\), \(\widetilde V=a x^2+b x^4\) generates an \(x^6\) term \(-8qb^2x^6\) in \(\mathcal R(\widetilde V)\). An ordinary quartic projection discards it and hence has an \(x^6\) residual. This computation illustrates the nonlinear obstruction; it does not prove that the discarded term is the dominant observable error.

Published ERG studies support both caution and exceptions. Morris finds polynomial truncations that initially approach the untruncated answer but then cease converging and generate spurious solutions. Aoki et al. describe coordinate-dependent projections and improved expansion about the potential minimum, including special large-\(N\) “perfect coordinates.” These are not universal residual-error theorems. <citations>24,25,26,27</citations>

A moment or cumulant closure similarly replaces unretained quantities in the hierarchy by functions of retained ones. Its reduced equations can be satisfied exactly while the reconstructed full distribution or generating function violates the full PDE. Some closures do not reconstruct a positive distribution at all; a full-PDE residual must then be defined on the formal generating series or another explicitly chosen reconstruction, not on a nonexistent probability density.

For the scalar Gaussian-channel normalizer, the hierarchy can be made completely explicit. Define \(\kappa_n(t,y)=\partial_y^nK(t,y)\); these are the tilted posterior cumulants. Differentiating the backward heat equation gives the **derived exact hierarchy**

$$
\partial_t\kappa_n=-\tfrac12\kappa_{n+2}
-\tfrac12\sum_{k=0}^n\binom nk\kappa_{k+1}\kappa_{n-k+1}.
$$

This is a derivative at fixed observation coordinate \(y\), not the Itô evolution along \(Y_t\). A closure replacing \(\kappa_{M+1},\kappa_{M+2}\) changes the right-hand sides at the top retained orders. Reconstructing \(\widetilde K\) and substituting it into \(\partial_t\widetilde K+\tfrac12(\partial_y^2\widetilde K+(\partial_y\widetilde K)^2)\) exposes a full generating-function defect, including terms outside the retained span. Without a consistent reconstruction, a residual in retained equations alone cannot certify the full expectation error.

### 2.2 Linear heat defect: Duhamel identity

**Derived identity.** Let \(\partial_sZ=L_sZ\), and let a sufficiently regular approximation have

$$
r_s=\partial_s\widetilde Z-L_s\widetilde Z.
$$

For the exact heat propagator \(H_{s,t}\),

$$
Z_t-\widetilde Z_t
=H_{0,t}(Z_0-\widetilde Z_0)-\int_0^tH_{s,t}r_s\,ds.
$$

In a norm in which the propagator contracts,

$$
\|Z_t-\widetilde Z_t\|
\leq\|Z_0-\widetilde Z_0\|+\int_0^t\|r_s\|ds.
$$

The contraction claim applies to forward heat evolution, not unrestricted backward heat inversion. Polynomial potentials and unbounded observables often require weighted spaces rather than a global sup norm.

For \(\widetilde Z=e^{-\widetilde V}\), the two residuals are related by

$$
r_s=-\widetilde Z_s d_s.
$$

An additive potential residual and a weight residual therefore have different error scales. Recovering \(V\) from \(Z\) requires a lower bound on \(Z\); recovering normalized expectations additionally requires control of normalization.

### 2.3 Nonlinear potential defect and dual weighting

**Derived exact difference equation.** For \(e=V-\widetilde V\) in forward covariance time,

$$
\partial_se=L_se-
\tfrac12\big[\dot C_s(\nabla V+\nabla\widetilde V)\big]\cdot\nabla e-d_s.
$$

Thus the error is propagated by a diffusion–transport operator depending on the exact and approximate solutions. It is not generally the bare heat convolution of \(d\). A stability estimate yields a bound of the form

$$
\|e_t\|\leq K(t,0)\|e_0\|+\int_0^tK(t,s)\|d_s\|ds.
$$

For a goal functional, the residual should instead be paired with the appropriate adjoint sensitivity. A residual large in an irrelevant direction can have little effect on the goal; a small residual in an unstable direction can matter greatly. These are the mechanisms behind residual/dual-weighted error analysis, not special privileges of RG terminology.

### 2.4 Observable defect along a diffusion

**Derived identity.** If \(dX_s=b_s(X_s)ds+\sigma_s(X_s)dW_s\), generator \(\mathcal A_s\), and \(\widetilde u(T,\cdot)=F\), define

$$
D_s=(\partial_s+\mathcal A_s)\widetilde u(s,\cdot).
$$

Itô’s formula gives

$$
\mathbb E[F(X_T)]-\widetilde u(0,X_0)
=\mathbb E\int_0^TD_s(X_s)ds.
$$

A terminal mismatch contributes \(\mathbb E[F(X_T)-\widetilde u(T,X_T)]\). Sufficient integrability to take expectations is essential. If the diffusion is stochastic localization, this is the localization-generator residual—not automatically a bare heat residual in the observation coordinates.

With the true backward solution \(u\), \((\partial_s+\mathcal A_s)u=0\). The identity explains the weak-error approach in Bally–Talay: errors are expanded through a backward parabolic equation, and the leading weak-error coefficient is an integrated quantity along the exact diffusion. Their theorem imposes smooth-coefficient and uniform Hörmander-type hypotheses. <citations>28,29</citations>

### 2.5 Four correction methods that must not be conflated

**Defect correction.** For a linear equation, solve

$$
(\partial_s-L_s)c=-r_s,\qquad c(0)=0,
$$

and replace \(\widetilde Z\) by \(\widetilde Z+c\). For a nonlinear equation, solve the linearized residual equation

$$
(\partial_s-D\mathcal R_s[\widetilde V])c=-d_s
$$

with compatible boundary/normalization conditions. This is a derived Newton-type construction. If the correction is restricted to the same inadequate ansatz, unresolved residual components can remain.

**Richardson/Romberg.** If a genuine asymptotic expansion holds,
\(A_h=J+c_ph^p+O(h^{p+1})\), then

$$
A_h^{\mathrm R}=\frac{2^pA_{h/2}-A_h}{2^p-1}
$$

removes the leading coefficient without knowing \(J\). Bally–Talay explicitly describe weak-error expansions and extrapolation across Euler step sizes. There is no corresponding guarantee merely from comparing successive cumulant orders: a cumulant order is not a step size with an established power law. <citations>30,31</citations>

**RG improvement of singular perturbation series.** Chen–Goldenfeld–Oono absorb secular terms into running amplitudes, introduce an arbitrary intermediate scale, and demand independence of that scale. Their examples connect the resulting RG equations to amplitude equations and multiple scales. This can remove growing perturbative errors; it does not supply a general posterior residual bound for arbitrary finite closures. <citations>32,33,34,35,36,37</citations>

**Renormalization conditions.** Krajewski–Martinetti prescribe relevant/marginal parameters at a low-energy scale rather than divergent boundary data at a high-energy scale. Their paper also organizes Wilsonian flows and ODE numerical expansions with Hopf algebras. This is the most direct retrieved bridge between RG boundary prescriptions and numerical composition. But physical renormalization conditions usually contain prescribed parameter values; they are not intrinsically truth-free. An estimator needs conditions whose values follow from its known input, architecture, or exact identities. <citations>38,39,40,41</citations>

## 3. Hopf algebras: the established bridges and the missing theorem

### 3.1 Feynman graphs, Birkhoff factorization, and counterterms

In the Connes–Kreimer graph Hopf algebra, multiplication is disjoint union and the coproduct schematically is

$$
\Delta\Gamma=\Gamma\otimes1+1\otimes\Gamma+
\sum_{\gamma}\gamma\otimes\Gamma/\gamma,
$$

where the sum is over admissible superficially divergent subgraphs, including appropriate disjoint unions and external-structure labels. It organizes nested/subdivergent contributions. <citations>42,43</citations>

For a regularized character \(\phi\) valued in Laurent series, choose a pole-part projection \(R\). With convolution \(*\),

$$
\phi=\phi_-^{-1}*\phi_+,\qquad
\bar\phi(\Gamma)=\phi(\Gamma)+\sum_\gamma\phi_-(\gamma)\phi(\Gamma/\gamma),
$$
$$
\phi_-(\Gamma)=-R\bar\phi(\Gamma),\qquad
\phi_+(\Gamma)=(I-R)\bar\phi(\Gamma).
$$

These are the recursive Birkhoff/BPHZ subtraction formulas. For general target algebras, multiplicativity requires suitable algebraic hypotheses, commonly a Rota–Baxter splitting—not an arbitrary linear error projection. <citations>44,45</citations>

### 3.2 The one-parameter RG statement

In the dimensional-regularization/minimal-subtraction setting, with grading automorphism \(\theta\) and counterterm character \(\gamma_-\), Connes–Marcolli write

$$
\mathrm{rg}_t=\lim_{z\to0}\gamma_-(z)\,
\theta_{tz}(\gamma_-(z)^{-1}).
$$

Under the locality/scale-independence structure of that setup, this is a one-parameter subgroup generated by the beta function. Their universal formulation has a canonical homomorphism \(\mathbb G_a\to U\), whose image represents the RG under a graded representation of the universal group. These are perturbative algebraic statements, not an identification of the universal group with the localization Markov semigroup. <citations>46,47,48,49</citations>

Three objects should remain distinct: Gaussian convolution is generally an analytic semigroup; Polchinski particles have a time-inhomogeneous Markov evolution; the algebraic RG is a subgroup in a character group. Reparametrization or formal inversion does not erase these distinctions.

### 3.3 Trees and numerical analysis: a genuine application

Lundervold–Munthe-Kaas explain B-series order conditions by matching rooted-tree coefficients with the exact flow. The Butcher group is the character group of the rooted-tree Connes–Kreimer Hopf algebra; backward error analysis obtains a modified vector field by a formal logarithm, expressed through the Eulerian idempotent. The rooted-tree algebra and the Feynman-graph algebra are related constructions, not interchangeable lists of generators. <citations>50,51,52</citations>

Krajewski–Martinetti go further toward the requested link: their Wilsonian analysis uses ordered Feynman diagrams, while their nonlinear differential-equation expansions have the B-series convolution composition law. They discuss renormalization conditions and effective scale-dependent couplings within the same framework. This is a concrete bridge, though not a moment-closure error theorem. <citations>39,53,38,40,41</citations>

Bronasco–Laurent provide a stochastic numerical counterpart: exotic aromatic series, Hopf-algebraic composition, postprocessing, and explicit modified-vector-field formulas for ergodic SDEs. Their backward-error theorems concern weak/invariant-measure accuracy under their assumptions, not finite-network Gaussian expectations. This nevertheless shows that correction coefficients can be derived algebraically from a method and a generator rather than fitted to exact output values. <citations>54,55,56,57</citations>

The Talay–Tubaro-type weak expansions and these tree expansions are compatible approaches to numerical error organization. The retrieved evidence supports both separately; it does not establish that every such expansion has a QFT-style Birkhoff decomposition. <citations>30,29,56</citations>

### 3.4 Moments, cumulants, and Wick renormalization

Ebrahimi-Fard–Patras–Tapia–Zambotti give a coalgebra/convolution treatment of classical cumulants and moments and describe Wick polynomials as a Hopf-algebra deformation induced by multivariate moments. The full text was not recovered in this search, so that attribution is limited to the abstract. <citations>58</citations>

**Independent formal construction:** on the polynomial Hopf algebra with primitive variable generators, let the unital functional \(M\) assign joint moments. Then \(K=\log_*M\) assigns joint cumulants and \(M=\exp_*K\) reconstructs moments. For example, evaluating the convolution exponential on a monomial sums over set partitions of its variable occurrences. The formal series terminate degree by degree. Crucially, the moment functional is generally **not a character** on this ordinary polynomial algebra: expectation does not multiply across dependent variables. A cumulant convolution algebra and the character group used for B-series therefore require an explicit bridge, not a silent identification.

What is missing is a theorem that combines these cumulant structures with a finite closure, a specified subtraction projection, and a residual-controlled observable error bound. Neither moment–cumulant inversion nor Wick ordering automatically supplies that theorem. No direct application of QFT Birkhoff counterterms to deep-ReLU cumulant propagation was identified in the retrieved sources; this is a search outcome, not a proof of absence.

## 4. A truth-free counterterm hypothesis for Gaussian ReLU expectations

### 4.1 An exact heat-defect identity in variables the estimator already knows

Let \(F:\mathbb R^d\to\mathbb R\) be a fixed finite feed-forward ReLU network output and let \(X\sim N(m_0,\Sigma_0)\), initially with \(\Sigma_0\succ0\). Define

$$
G(m,\Sigma)=\mathbb E[F(m+\Sigma^{1/2}Z)],\qquad Z\sim N(0,I).
$$

**Derived identity:** for a symmetric covariance direction \(H\),

$$
D_\Sigma G[H]=\tfrac12H:\nabla_m^2G,
\qquad G(m,0)=F(m).
$$

The derivative uses symmetric-matrix directions, avoiding off-diagonal coordinate conventions. For positive-definite covariance, Gaussian smoothing justifies differentiation for a finite ReLU network’s at-most-linear growth. At zero covariance, use limits or weak derivatives; ordinary classical second derivatives of \(F\) omit its kink measures.

Let \(A(m,\Sigma)\) be the deterministic cumulant-propagation approximation. Define its heat defect tensor by

$$
\mathcal D_A[H]=D_\Sigma A[H]-\tfrac12H:\nabla_m^2A.
$$

Under the Gaussian observation channel, the posterior covariance and mean obey

$$
\Sigma_t=(\Sigma_0^{-1}+tI)^{-1},\qquad
 dm_t=\Sigma_t dW_t,\qquad \dot\Sigma_t=-\Sigma_t^2.
$$

Itô’s formula now gives the exact identity

$$
J-A(m_0,\Sigma_0)
=-\int_0^\infty\mathbb E\big[\mathcal D_A(m_t,\Sigma_t)[\Sigma_t^2]\big]dt,
\qquad J=\mathbb E[F(X)],
$$

provided \(A(m,0)=F(m)\), endpoint convergence is uniformly integrable, and the defect integral exists. Otherwise add
\(\mathbb E[F(X)-A(X,0)]\). This is the proposed estimator’s integrated heat-defect statement in a precise convention. It is a derivation here, not a theorem from the surveyed papers.

A useful finite-horizon variant is to choose a deterministic decreasing covariance path \(\Sigma(s)\), from \(\Sigma_0\) to zero, and let \(dm_s=(-\dot\Sigma(s))^{1/2}dW_s\). Then \(m_T\sim N(m_0,\Sigma_0)\) and

$$
J-A(m_0,\Sigma_0)
=\int_0^T\mathbb E\big[\mathcal D_A(m_s,\Sigma(s))[\dot\Sigma(s)]\big]ds.
$$

This finite-horizon construction is a Gaussian martingale realization; it need not be the identical finite-time Gaussian-channel filtration. It avoids an infinite-time endpoint without claiming a new localization equivalence.

### 4.2 The computability obstruction

The integral remains an expectation over random posterior means. Its exact evaluation is not made cheap by being called a defect. For Gaussian input, however,

$$
m_t\sim N(m_0,\Sigma_0-\Sigma_t).
$$

Consequently deterministic quadrature or another controlled Gaussian approximation can evaluate the defect expectation. If that secondary computation has certified total error \(\eta\), then adding its estimate of the signed integral corrects \(A\) up to \(\eta\), plus any endpoint error. That is a genuinely truth-free correction: it uses \(F\), the input distribution, derivatives of \(A\), and a controlled residual integration—not labeled exact network expectations.

It is useful only if the residual is easier to integrate than the original network. Possible reasons include cancellations, concentration near a small number of activation boundaries, low-rank covariance dependence, or a smaller effective dimension. None follows automatically from the identity.

### 4.3 Conditions that can determine counterterms without output fitting

**Proposed design.** Choose correction functions \(B_\alpha\) and set

$$
A^{\rm ren}=A+\sum_\alpha c_\alpha B_\alpha.
$$

Determine their coefficients from exact structural conditions and residual cancellation:

1. **Zero-noise matching:** \(A^{\rm ren}(m,0)=F(m)\). Counterterms must vanish there, or supply an explicitly known boundary mismatch.
2. **Heat consistency:** cancel selected contractions or asymptotic coefficients of \(\mathcal D_{A^{\rm ren}}\), not just a mismatch in retained cumulant ODEs.
3. **Gaussian composition:** enforce, to a chosen order,
   \(A^{\rm ren}(m,\Sigma_1+\Sigma_2)=P_{\Sigma_1}[A^{\rm ren}(\cdot,\Sigma_2)](m)\).
   Exact enforcement everywhere would recover the original heat problem; finite tests only control selected components.
4. **Path consistency:** different admissible covariance decompositions should yield the same endpoint estimate to the target order. This plays the role of independence from an arbitrary renormalization scale.
5. **Exactly solvable sectors:** require exactness for affine networks and single Gaussian ReLU units, and preserve known input symmetries. These supply analytic benchmark conditions, not empirical fits.
6. **Normalization and affine observables:** if correcting a reconstructed law rather than one output, preserve total mass and the prescribed moments that truly are known. Hidden-layer moments must not be called known merely because the closure predicts them.

For a unit with preactivation \(a\sim N(\mu,\sigma^2)\), an exact condition available without simulation is the derived Gaussian integral

$$
\mathbb E[a_+]=\sigma\varphi(\mu/\sigma)+\mu\Phi(\mu/\sigma).
$$

This helps select a kink-aware basis. It does not make all hidden preactivations Gaussian.

For small corrections, a projected system can take the form

$$
\sum_\alpha c_\alpha\,\ell_\beta(\mathcal D_{B_\alpha})
=-\ell_\beta(\mathcal D_A),
$$

supplemented by boundary/normalization constraints. The tests \(\ell_\beta\) may be analytic coefficient functionals or deterministically computed weighted residual integrals. If \(A\) is nonlinear in the counterterm coefficients, use the corresponding Jacobian and iterate. Rank deficiency means the conditions do not determine unique counterterms; instability means that increasing correction order can worsen the estimator.

### 4.4 Where Hopf algebra could enter—and what must be proved

**Hypothesis, not an established application.** Build a graded combinatorial algebra whose generators represent connected cumulant contributions and heat-defect insertions. Let the coproduct encode how a connected contribution contains smaller correction substructures. Let a character encode composition of propagation/correction steps. Derive corrections recursively from a selected “bad-part” projection and a normalization prescription.

Three gates are essential:

- Prove that the coproduct actually models the network/closure’s substitution and composition operations. A ReLU computation graph is not automatically a Feynman graph Hopf algebra.
- Specify the grading and expansion parameter. Depth, covariance size, cumulant order, and integration step size are not interchangeable small parameters.
- Prove algebraic compatibility of the subtraction projection and then establish analytic control of the remaining heat defect. A formal Birkhoff decomposition alone does not prove a convergent or accurate estimator.

The immediately defensible starting point is B-series-style order matching or a modified-equation correction for the cumulant propagation algorithm, followed by the observable heat-defect bound. A graph-Hopf presentation is valuable if it resolves nested correction combinatorics; it is not required for the first useful truth-free correction.

### 4.5 ReLU-specific failure modes and the attainable guarantee

A Taylor expansion using ordinary derivatives inside activation regions misses the distributional Hessian supported on kink boundaries. Gaussian smoothing converts those terms into finite boundary contributions. When \(m\) lies on an activation boundary, leading corrections can scale as \(\sqrt\varepsilon\) under covariance \(\varepsilon\Sigma\), rather than integer powers of \(\varepsilon\). Away from boundaries, crossing contributions can instead be exponentially small. Thus an assumed integer-power Richardson expansion needs justification or smoothing; it cannot be imported mechanically.

Finite moment information also does not generally identify the expectation of a piecewise-linear observable under an arbitrary hidden-layer law. Structural conditions can select a consistent closure without proving that it approximates the true distribution. The guarantee must come from an independently controlled residual norm or integral, not from satisfying a few renormalization conditions.

The feasible claim is conditional and precise: **counterterms determined from boundary data, heat-consistency conditions, and controlled residual integration can reduce or bound the Gaussian expectation error without fitting exact output values.** The stronger claim—universal, inexpensive, exact counterterms from finitely many cumulants—is not supported. The central research question is whether the network’s defect has a lower-complexity representation than its full Gaussian integral.
