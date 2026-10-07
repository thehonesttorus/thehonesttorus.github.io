# Correlated boundary-response compilation
## A source–observable calculus, minimum-energy Gaussian transport, and exact conditional wall integration

**Date:** 7 October 2026  
**Status:** research development with proofs and offline finite checks. The identities below are derived here from the supplied manuscript and standard Gaussian/operator tools. No claim of publication priority is made. No current competition implementation, challenge weights, private ground truth, or official FLOP meter was used.

## 0. What is new relative to the supplied manuscript

The uploaded *The missing third-to-fourth cumulant feed* identifies the omitted fourth-order multiplicities `(31),(211),(1111)`, derives their scalar influence, gives a Gaussian score and gate-facet representation, and supplies a small-factor conditional algorithm. It explicitly does not establish small common-factor rank or a cheap deep closure for the supplied networks.

The preceding review supplied an all-row formula at an independent Gaussian reference. This note develops a different result:

> For an arbitrary positive-definite correlated Gaussian reference and a supplied third-order pair source, a fixed coherent sum of all requested missing fourth-order responses has an exact expectation representation costing O(mn+n^2) per evaluation after O(mn^2+n^3) preparation. The integrand does not require unknown third-order baseline moments. One Gaussian coordinate can be integrated analytically through all its gate crossings in the same asymptotic per-evaluation complexity.

Here `m` is the number of downstream query rows, not an activation mean. This is **one fixed coherent linear combination of m responses**, not a free computation of m independently usable responses. At m=n the evaluation is quadratic; constructing the supplied covariance, source, and compiled matrices is separately charged.

A second result turns the Gaussian score into a quadratic **minimum-energy transport velocity**. This realizes the source as the tangent of an ordinary positive pushforward law and replaces source multiplication by a first derivative of the readout. Its energy minimality is a Hilbert-space statement, not a proof of optimal computational cost or minimum estimator variance.

The third result is a Ward change of integrand. It preserves the pairing, not the score as a universal probability tangent. It can be added as an exactly zero-mean control without silently replacing the ledger's pair source by a full scale source.

None of these results asserts a dimension-independent integration-node count, an autonomous two-gain closure, or a numerical prediction of the reported deep-layer deficit.

## 1. Fixed coherent target and notation

Let H in R^n have enough moments for the derivatives below, let

\[
m=\mathbb EH,\qquad x=H-m,\qquad C=\mathbb E[xx^\top].
\]

Let W in R^{m_out x n} have rows w_r^T. Let omega in R^{m_out} be a **fixed** coherent readout vector, possibly signed. Define B=W^{odot 2}, C_off=C-diag(diag C),

\[
V_r=w_r^\top Cw_r,\qquad D_r=\sum_i w_{ri}^2 C_{ii},
\]

and compile

\[
M=B^\top\operatorname{diag}(\omega)B,\quad c_4=\operatorname{diag}M,
\quad N=M\odot C_{\rm off}.
\]

M,N are symmetric. Neither needs to be positive when omega is signed. N is not being declared a Dirichlet operator or a covariance.

For a sample x put y=Wx and s=B(x^{odot 2}). Set

\[
\mathcal P(x)=\omega^\top(y^{\odot4}-3s^{\odot2})
               +2c_4^\top x^{\odot4}.
\]

### Proposition 1: exact pooled omitted cumulant

The coherent sum of the omitted `(31),(211),(1111)` fourth-cumulant sectors is

\[
\mathcal K_{\rm miss}
=\mathbb E\mathcal P(x)-3\omega^\top(V^{\odot2}-D^{\odot2})
 +6\sum_{ij}M_{ij}(C_{\rm off})_{ij}^2.
\tag{1}
\]

**Proof.** For one row, the retained `(4),(22)` fourth moments are `3 S2^2 - 2 S4`. Subtract them from S1^4. The complete covariance-pairing subtraction is 3V_r^2, and the retained part is 3D_r^2+6 sum_{i!=j} w_ri^2 w_rj^2 C_ij^2. Summing with weights omega yields (1); M collects the latter pair contractions. All ordered-sum conventions are retained. QED.

This identity is basis-dependent in exactly the same sense as the specified multiplicity classes: W and the sectors refer to the original neuron coordinates.

## 2. Compiling the moving-mean correction removes a baseline oracle

Let L be the derivative of a probability law, with L1=0, and let d=L(H) be the derivative of its mean. Let all unadorned means, covariances and coefficients be evaluated at the reference law.

Define

\[
\mathcal J(x)=\mathcal P(x)-6(\omega\odot V)^\top y^{\odot2}
 +6(\omega\odot D)^\top s+12x^\top Nx.
\tag{2}
\]

### Theorem 2: coherent influence without baseline third moments

If L(f)=E[s_src f], then

\[
\boxed{\dot{\mathcal K}_{\rm miss}
 =\mathbb E\left[s_{\rm src}\mathcal J(x)-D_d\mathcal P(x)\right].}
\tag{3}
\]

Moreover,

\[
\boxed{
D_d\mathcal P
=4(\omega\odot y^{\odot3})^\top Wd
 -12(\omega\odot s)^\top B(x\odot d)
 +8(c_4\odot d)^\top x^{\odot3}.}
\tag{4}
\]

**Proof.** Differentiation of the reference centering gives

\[
\frac d{dt}\mathbb E_t\mathcal P(H-m_t)\big|_0
=L\mathcal P-\mathbb E[\nabla\mathcal P]^\top d.
\]

Since Ex=0, covariance differentiation is `dot C = L(xx^T)`. Therefore

- `dot V_r = L(y_r^2)`;
- `dot D_r = L(s_r)`;
- the derivative of the last term of (1) is `12 L(x^T N x)`.

These terms give L J. Finally E[grad P]^T d=E[D_d P]. Direct differentiation of P gives (4). QED.

The source report stores b=E grad P as baseline third-order contractions. Equation (3) integrates b^T d in the same pass instead. No population third moment has been approximated or set to zero. The required population mean m, covariance C, and mean derivative d remain.

If omega, W, the gain denominator, or the fitted output projection changes with the perturbation, its derivative contributes separately. In particular, for g4=kappa4/(3v^2),

\[
\dot g_4=\dot\kappa_4/(3v^2)-2g_4\dot v/v.
\]

Those contributions cannot be hidden inside a fixed-omega claim.

## 3. Exact Gaussian pair source at arbitrary covariance

Let Z~N(mu,Sigma), Sigma positive definite, and H=relu(Z). Write Q=Sigma^{-1}, q=Q(Z-mu), sigma_i^2=Sigma_ii, alpha_i=mu_i/sigma_i.

The supplied source is specified by tau_i=T_iii and Gamma_ij=T_iij for i!=j, with Gamma_ii=0 and no fully distinct source entries. Gamma need not be symmetric. Its three symmetric tensor permutations have the same coefficient.

Compile r_i=sum_j Gamma_ij Q_ij. The source score is

\[
\boxed{
 s_T(Z)=\frac16\tau^\top(q^{\odot3}-3\operatorname{diag}Q\odot q)
 +\frac12(q^{\odot2}-\operatorname{diag}Q)^\top\Gamma q-q^\top r.}
\tag{5}
\]

This is the Gaussian Hermite contraction, with the factor 1/2 accounting for three permutations of each `(i,i,j)` source entry.

Its post-ReLU mean derivative is explicit:

\[
\boxed{d_i=-\frac{\tau_i\alpha_i\phi(\alpha_i)}{6\sigma_i^2}.}
\tag{6}
\]

**Proof.** The i-th marginal readout depends only on z_i, so all mixed mean derivatives vanish. Its third mean derivative is `-alpha_i phi(alpha_i)/sigma_i^2`; multiply by tau_i/6. Equivalently integrate the normal derivative of the marginal density at zero. QED.

The score represents a real law tangent. For example,

\[
p_\epsilon=p_0\frac{(1+\epsilon s_T/2)^2}{1+\epsilon^2\mathbb E s_T^2/4}
\]

is positive and normalized, and differentiates to p_0 s_T at zero. Its first cumulant variation is T in order three and zero in every other order. This is an infinitesimal statement, not an exact model of a finite deep-network non-Gaussian law.

### Preparation and evaluation cost

Given the population C and source data:

- construct M and V with two matrix products, O(m_out n^2);
- factor/invert Sigma, O(n^3), unless already available;
- build N, r and other vectors, O(n^2+mn);
- evaluate (3) per Gaussian sample with a fixed number of matrix-vector products, O(m_out n+n^2).

For m_out=n, the latter is quadratic even at fully correlated Sigma. Baseline C can be computed from bivariate Gaussian ReLU moments to prescribed numerical precision; that preparation and its error remain separate. No small-factor decomposition is assumed.

## 4. A source can be realized as a minimum-energy input transport

This section gives a different implementation of the same derivative, rather than a new source assumption.

For the Gaussian density p define the weighted Ornstein-Uhlenbeck number operator

\[
\mathsf N f=(z-\mu)\cdot\nabla f-\operatorname{tr}(\Sigma\nabla^2f).
\]

Gaussian integration by parts gives

\[
\mathbb E[f\mathsf Ng]=\mathbb E[(\nabla f)^\top\Sigma\nabla g].
\]

The third score satisfies N s_T=3s_T. Define

\[
\boxed{v_T(z)=\frac13\Sigma\nabla s_T(z).}
\tag{7}
\]

### Theorem 3: transport, source duality, and minimal energy

For suitable polynomial-growth tests f, including the local ReLU polynomials above,

\[
\boxed{\mathbb E[s_Tf]=\mathbb E[v_T\cdot\nabla f].}
\tag{8}
\]

The pushforward random variable

\[
Z_\epsilon=Z+\epsilon v_T(Z)
\]

is an ordinary probability model for every epsilon and has precisely the source tangent s_T at zero. Among square-integrable vector fields v inducing that same tangent, v_T minimizes

\[
\mathcal A(v)=\mathbb E[v^\top\Sigma^{-1}v],
\]

and

\[
\boxed{\mathcal A(v_T)=\tfrac13\mathbb E s_T^2.}
\tag{9}
\]

**Proof.** Equation (8) is the Dirichlet identity applied to N s_T=3s_T. Differentiating f(Z+epsilon v_T) yields (8); for locally Lipschitz polynomial-growth ReLU tests, the nonsmooth hyperplanes have Gaussian measure zero and polynomial bounds dominate difference quotients. The map need not be globally invertible to define its pushforward.

If v is another realization, w=v-v_T satisfies E[w dot grad f]=0 on the test domain. By Sobolev approximation use f=s_T. Then E[w dot grad s_T]=0, so v_T and w are orthogonal in the kinetic inner product. Thus A(v)=A(v_T)+A(w). Finally A(v_T)=E[grad s_T^T Sigma grad s_T]/9=E[s_T^2]/3. QED.

For a general symmetric whitened source tensor T_tilde, E s_T^2=||T_tilde||_F^2/6, hence the minimum energy is ||T_tilde||_F^2/18. This says nothing by itself about evaluation FLOPs or optimal Monte Carlo variance.

For the pair source, the quadratic velocity has a cheap expression. Put u=q^{odot2}-diag Q. Then

\[
\boxed{
 v_T=\frac13\left[\tfrac12\tau\odot u
       +q\odot(\Gamma q)+\tfrac12\Gamma^\top u-r\right].}
\tag{10}
\]

**Proof.** Differentiate (5) with respect to q and use Sigma grad_z=grad_q. QED.

### A genuine integration gauge, not a topological invariance claim

For a vector field v define the Gaussian divergence

\[
\delta_\gamma v=(z-\mu)^\top\Sigma^{-1}v-\nabla\cdot v.
\]

For the velocity (7), delta_gamma v_T=s_T. Its product rule gives

\[
\boxed{s_T f-v_T\cdot\nabla f=\delta_\gamma(fv_T).}
\]

The Gaussian expectation of the right side is zero under the integrability conditions used above. Thus score multiplication and source transport differ by a genuine mean-zero divergence. This is the precise, state-dependent sense in which an integration representative may be changed. It is not a claim that arbitrary cyclic coboundaries, K-theory homotopies, or changes to the weights preserve the requested real-valued mean.

This gives the alternative response integrand

\[
\boxed{
 \dot{\mathcal K}_{\rm miss}
 =\mathbb E\left[(v_T\odot\mathbf1_{Z>0})^\top\nabla_x\mathcal J(x)
                       -D_d\mathcal P(x)\right].}
\tag{11}
\]

Only first derivatives of ReLU appear. No facet-delta products are evaluated. The gradient is

\[
\nabla\mathcal J
=4W^\top[\omega\odot(y^{\odot3}-3V\odot y)]
-12x\odot B^\top[\omega\odot(s-D)]
+8c_4\odot x^{\odot3}+24Nx.
\tag{12}
\]

Again this costs O(mn+n^2) per node. A finite nonzero epsilon does introduce higher-order changes to the law; this construction is not an uncorrected replacement of the competition input distribution.

### The operator-geometric meaning

Whitening turns the Gaussian Dirichlet form into d* d. Equation (7) solves the source-divergence equation in the gradient subspace. The remainder freedom consists of divergence-free currents orthogonal to this gradient. Thus the source/readout pairing can be evaluated either as multiplication by a score or as a flux paired with the differential of the observable. This uses the analytic, measured part of the operator calculus, not an integer index.

## 5. A full-scale Ward change of integrand

This section concerns the *full* third component of the manuscript's regular scale tangent:

\[
T^{\rm sc}_{ijk}=\tfrac12(\mu_i\Sigma_{jk}+\mu_j\Sigma_{ik}+\mu_k\Sigma_{ij}).
\]

It generally differs from the matched pair source. The associated tangent operator is

\[
\mathcal L_3 f=\tfrac14\mathbb E[D_\mu\Delta_\Sigma f].
\]

Let

\[
\ell=\mu^\top Q(Z-\mu),\quad a^2=\mu^\top Q\mu,
\quad\mathscr E f=Z\cdot\nabla f.
\]

### Theorem 4: Ward reduction to a first derivative

\[
\boxed{
\mathcal L_3 f=
\tfrac14\mathbb E\left[\ell(\mathscr E f-f)-(\ell^2-a^2)f\right].}
\tag{13}
\]

**Proof.** Gaussian integration by parts first gives

\[
E D_\mu\Delta_\Sigma f
=E\Delta_\Sigma D_\mu f
=E[(Z-\mu)\cdot\nabla D_\mu f].
\]

Use `[mathscr E,D_mu]=-D_mu` to write the final expression as

\[
E[D_\mu(\mathscr E f-f)-D_\mu^2f].
\]

The first and second Gaussian scores in direction mu are ell and ell^2-a^2. Substitute them. Smooth approximation or distributional Gaussian integration extends the equality to polynomial-growth ReLU polynomials. QED.

For a homogeneous degree-k f, this simplifies to

\[
\mathcal L_3f=\tfrac14 E[((k-1)\ell-\ell^2+a^2)f].
\tag{14}
\]

The influence J is not homogeneous because it contains frozen means and covariance coefficients. Use its actual Euler derivative, not `4J`.

The full third score in whitened coordinates xi and a=L^{-1}mu is

\[
s_3^{\rm sc}=\tfrac14(a\cdot\xi)(\|\xi\|^2-n-2).
\]

The minimum-energy physical velocity of that source is

\[
v_3^{\rm sc}
=\tfrac1{12}\left[\mu(\|\xi\|^2-n-2)+2(a\cdot\xi)(Z-\mu)\right].
\]

This is an exact structured transport for that tangent; no such form is asserted for an arbitrary source.

### A safe way to use the Ward identity on a different source

Let

\[
C_W=\tfrac14[\ell(\mathscr E\mathcal J-\mathcal J)
 -(\ell^2-a^2)\mathcal J]-s_3^{\rm sc}\mathcal J.
\]

Then E C_W=0. For a matched pair-source estimator Z_T from (3), the random variable

\[
Z_T+\beta C_W
\]

has the same expectation for any fixed beta. The variance-minimizing beta, when Var C_W>0, is `-Cov(Z_T,C_W)/Var(C_W)`. Learning beta on the same production points and retrospectively modifying their estimates does not inherit this simple unbiasedness argument. Independent pilots or predictable updates are required.

Neither the Ward nor the transport representation universally reduces variance. They are exact alternative representatives of a fixed pairing; their cost and variance must be compared.

## 6. Exact line completion with arbitrary correlated covariance

Let Sigma=L L^T, choose a unit vector a in whitened coordinates, and decompose

\[
\xi=\eta+at,\quad\eta\sim N(0,I-aa^\top),\quad t\sim N(0,1),\quad\eta\perp t.
\]

Thus

\[
Z=z_0+bt,\quad z_0=\mu+L\eta,\quad b=La.
\]

Conditional on eta, every gate changes at most once, at t_i=-z0_i/b_i when b_i!=0. On each resulting interval,

\[
x(t)=A+Bt.
\]

Here A,B are interval-specific vectors, not the earlier matrix B=W^2. To avoid confusion the implementation calls them `a,b` within its interval state.

The source score s_T(z0+bt) is a cubic polynomial in t. J is quartic, and D_d P is cubic. Therefore the response integrand on each interval is a polynomial of degree at most seven.

### Theorem 5: analytic integration along the omitted coordinate

For interval endpoints l,u, define

\[
I_k(l,u)=\int_l^u t^k\phi(t)\,dt.
\]

Then

\[
I_0=\Phi(u)-\Phi(l),\qquad I_1=\phi(l)-\phi(u),
\]

\[
I_k=l^{k-1}\phi(l)-u^{k-1}\phi(u)+(k-1)I_{k-2}.
\tag{15}
\]

If c_{jk} are the integrand coefficients on interval j,

\[
\boxed{\mathbb E[Z_T\mid\eta]=\sum_j\sum_{k=0}^7c_{jk}I_k(t_j,t_{j+1}).}
\tag{16}
\]

**Proof.** Gaussian conditioning gives the independent standard t. The local activation is affine on each interval, so the polynomial-degree assertion follows. Equation (15) is integration by parts using phi'=-t phi. Summation proves (16). QED.

This is an exact integral identity with elementary arithmetic and univariate normal special functions. A float64 implementation is not an exact-arithmetic or interval-certified numerical answer.

The theorem covers the local one-ReLU response in this note. Restricting an entire deep network to a line can create many more than n breakpoints.

## 7. Why the sweep stays quadratic rather than recomputing every interval

At a crossing, only one coordinate changes its affine coefficients. If the crossing coordinate is j,

\[
\delta A_j=\operatorname{sign}(b_j)z_{0,j},\qquad
\delta B_j=|b_j|.
\]

The centered activation is continuous at the crossing even though its two coefficients change.

Maintain the following interval polynomial data:

- WA and WB;
- W^2 A^2, 2W^2(A odot B), W^2 B^2;
- W^2(A odot d), W^2(B odot d);
- NA, NB and the three coefficients of x(t)^T N x(t);
- five coefficients of c4 dot x(t)^4;
- four coefficients of (c4 odot d) dot x(t)^3.

Changing coordinate j updates each vector by a scaled stored column of W, W^2, or N. The work per crossing is O(m+n). Constructing its degree-seven integrand costs O(m) and integrating it costs constant work. There are at most n crossings.

Thus

\[
\boxed{C_{\rm line}=O(n^2+mn+n\log n).}
\tag{17}
\]

This includes generic dense source-line preparation; fixed direction-dependent quantities can be cached.

### An array-only compiler removes the width-dependent Python loop

Sort the crossing indices pi. All WA/WB and W^2 updates are prefix sums of scaled columns in that order. The apparent sequential updates of the quadratic form can also be compiled.

Let N_pi be N reordered by pi and L_pi its strict lower triangle. If deltaA,deltaB are the crossing increments, the j-th row just before its update sees

\[
(NA_\mathrm{base})_{\pi_j}+(L_\pi\delta A)_j,
\]

and similarly for B. Hence the quadratic-coefficient increments are available from two triangular matrix-vector products and prefix sums. The remaining fourth/cubic monomial coefficients use prefix sums of scalar increments.

The supplied `vectorized.py` implements this. Its Python loops are bounded by the polynomial degree, not by n. Its peak workspace is O(mn+n^2), larger than the scalar sweep's workspace but still without a three- or four-index tensor. Actual gathering, triangular masking, prefix operations, memory traffic, and precision must be metered before using it in a competition.

## 8. A concrete measured groupoid and its noncommuting operations

The line disintegration has a literal completion relation. Two Gaussian descriptions (eta,t) and (eta,u) are related when they have the same eta. On each fiber, composition of kernels is

\[
(K_1*K_2)(\eta;t,u)=\int K_1(\eta;t,s)K_2(\eta;s,u)\,d\gamma_1(s).
\]

The constant kernel 1 is the conditional-integration projection P_a. It satisfies P_a^2=P_a. Multiplication by a readout generally does not commute with it. For bounded f,

\[
P_a M_f P_a=M_{P_af}P_a.
\]

These are concrete operator operations on the disintegrated probability space. This is a continuous-fiber pair relation with a conditional Gaussian measure, not a claim that the network has already been realized as an étale Penrose groupoid.

Let F_a=2P_a-I. On the constant cyclic vector, for square-integrable f,

\[
[F_a,M_f]1=-2(I-P_a)f,
\]

hence

\[
\boxed{
\operatorname{Var}(f)-\operatorname{Var}(P_af)
=\|(I-P_a)f\|_2^2
=\tfrac14\|[F_a,M_f]1\|_2^2.}
\tag{18}
\]

For unbounded polynomial-growth f these expressions are understood on the displayed vector and appropriate dense domains, or through bounded truncations. No bounded-commutator spectral triple, compactness, or index theorem is inferred.

This is the operative NCG connection: keep the projection representing forgetting/integration in the algebra, and measure its noncommutation with the requested observable. The line-fill formula makes that projection computable in this case.

## 9. Predictable adaptive directions and resource allocation

For every fixed a, (16) preserves the target expectation and reduces variance. If a_t is chosen using weights and samples observed strictly before a fresh outer draw, then the conditional line estimator Y_t satisfies

\[
\mathbb E[Y_t\mid\mathcal F_{t-1}]=\dot{\mathcal K}_{\rm miss}.
\]

For a predetermined number N of draws,

\[
\mathbb E\left|\frac1N\sum_tY_t-\dot{\mathcal K}_{\rm miss}\right|^2
=\frac1{N^2}\sum_t\mathbb E\operatorname{Var}(Y_t\mid\mathcal F_{t-1}).
\tag{19}
\]

**Proof.** Apply the tower property for unbiasedness; centered increments are martingale differences, making cross terms vanish. QED.

Choosing a direction using the very current sample being corrected requires a different conditional-law argument. Data-dependent stopping and self-normalized ratios also require separate treatment.

### A small rigorous direction-selection certificate

Let b_f=E[xi f(xi)] for the fixed scalar response integrand in whitened coordinates. For a unit a,

\[
\boxed{\|f-P_af\|_2^2\geq(a^\top b_f)^2.}
\]

**Proof.** The variable a^T xi is a unit-norm mean-zero function of the integrated coordinate, independent of eta. It is orthogonal to every function of eta. Therefore its inner product with f-P_af is a^T b_f; apply Cauchy--Schwarz. QED.

When b_f is known and nonzero, a=b_f/||b_f|| certifies removal of at least ||b_f||^2 variance. It need not be the globally best direction and gives no guarantee when that projection vanishes. Estimating b_f requires pilot data and an uncertainty allowance; using a from previous batches preserves (19). This provides a mathematically justified first adaptive signal instead of assuming a weight PCA direction is optimal.

An independent sampler, an expander schedule of previously selected directions, or a deterministic schedule can select which conditional expectation to apply. None of these choices by itself creates a useful gap on the readout. In particular, random single-coordinate heat-bath updates on an n-dimensional product Gaussian have a 1/n gap on first-chaos functions; a universal dimension-free improvement cannot be asserted.

For independent outer samples, a line estimator with variance V_line and cost c_line beats a point estimator with variance V_point and cost c_point at matched sampling work only if

\[
V_{\rm line}/V_{\rm point}<c_{\rm point}/c_{\rm line},
\]

with preparation costs also included. Rao–Blackwellization alone proves only the variance inequality, not this cost inequality.

## 10. What has been removed and what has not

Removed for the stated local coherent response:

1. A need to assume diagonal-plus-small-rank preactivation covariance.
2. A per-output-row dense pair loop at every integration node.
3. A need for baseline third moments solely to account for moving centering.
4. A need to materialize missing fourth-cumulant entries.
5. A need to sample one chosen Gaussian coordinate or explicitly handle delta distributions.
6. A need for width-dependent Python loops in the array-only line implementation.

Not removed:

1. Preparation of the actual local covariance and source tangent.
2. Any cost of assembling the pair source from the chain's separately retained histories.
3. Correct identification of the physical incoming tangent and coherent output projection.
4. Error of replacing a true deep non-Gaussian law by a Gaussian reference plus a tangent.
5. The remaining outer integral and its node count at the required accuracy.
6. Finite-amplitude or autonomous depth propagation error.
7. Numerical conditioning and signed cancellation, especially at tiny residuals.
8. The cost of n independently reusable corrections when only a single coherent scalar was compiled.

The compiled construction is a local evaluator, not an established replacement for the existing whole-network estimator.

## 11. Executed checks and evidence scope

`results.json` records all inputs used in the small Gaussian test through their definitions in `checks.py`.

### Arbitrary finite laws

Eighteen positive finite-law tests checked (3) against a positive-law finite difference and, at small n, an independently constructed fourth-cumulant tensor. Maximum scaled derivative discrepancy was 7.18e-11; tensor identity discrepancy was 2.44e-15.

### Correlated Gaussian line identities

Fourteen constructed full-rank covariance tests checked the cubic source restriction, incremental/rebuilt interval sweeps, and independent scalar quadrature. Discrepancies were at floating-point rounding scale.

### Actual one-layer Gaussian tangent, not a challenge network

The principal example had

\[
\mu=(0.45,-0.30),\quad
\Sigma=\begin{pmatrix}1.30&0.43\\0.43&0.80\end{pmatrix},
\]

four signed query rows and fixed coherent weights, diagonal source `(0.15,-0.09)`, and pair source `Gamma_12=0.12`, `Gamma_21=-0.07`.

Direct score integration, minimum-energy transport integration, and conditional line integration all gave approximately

\[
0.0440303793847996.
\]

An independent first-variation finite difference gave 0.0440303793700254 at its recorded step. The full-scale source and its Ward representative both gave approximately 0.0167328331378553. These are different source tangents and intentionally different numbers.

Quadrature-estimated integrand variances were:

- score representative: 0.4482780308;
- transport representative: 0.1047978677;
- score with an oracle Ward control: 0.2495012970;
- analytically integrated line: 0.0141673470.

The line reduced variance by a factor 31.64 in this two-dimensional example. This is not a sample-cost speedup, not a high-dimensional variance theorem, and not a leaderboard measurement. The oracle control coefficient is a diagnostic, not an independently learned production coefficient.

Gaussian quadrature was split at physical gates and truncated at ten marginal standard deviations for its two-dimensional reference calculations. Line integrals use exact moment formulas evaluated in float64. The implementation does not produce outward-rounded quadrature or tail certificates.

### Array-only realization

The prefix implementation agreed with the incremental sweep through n=m=256 on constructed inputs. Those cases use arbitrary symmetric baseline matrices for algebra checks, not a claim that every such matrix is a true rectified-Gaussian covariance. The timing values are local single-call diagnostics, not official FLOP counts or grader residual time.

### Gaussian transport identity

Polynomial test pairings through several mixed degrees and the kinetic-energy identity agreed to rounding. The minimum transport energy was 0.0065863415181926, equal to E[s_T^2]/3. Finite differences through the actual pushforward `Z+epsilon v_T(Z)` also approached the same response.

## 12. Relation to the wider research programme

The geometric operation here is a true conditional filling: eta specifies an address and t the omitted fiber. Along that fiber, the supplied weights explicitly determine where local rules change. The mean of the fiber is computed, not guessed. The algebra of multiplication and conditional integration records which nonlinear queries commute with forgetting and what they lose when they do not.

The full tangent can also be represented by a minimum-energy Gaussian transport. This is the literal input-transformation direction the earlier programme sought, but at a specified infinitesimal source order. It is not a claim that a quadratic deformation samples the true deep law at finite amplitude.

The strongest remaining theorem would bound the number and cost of conditional fills required for the actual weight-dependent response, including accumulated baseline and reference errors. It may involve adapted directions, retained boundary contexts, short operator memory, or a different integrable family. The current work gives executable local primitives and exact pairings on which such a theorem could act.

It would be incorrect to infer from these formulas that the network already has a Penrose lattice, a uniform expander gap, a small hidden Gibbs separator, or a topological formula for its mean. It would also be incorrect to infer a tensor lower bound merely from the number of formal missing entries: the all-row pooled calculation explicitly avoids enumerating them.

## Sources and provenance

1. **Supplied manuscript:** *The missing third-to-fourth cumulant feed*, checkpoint NCG-20261007-A, 7 October 2026, uploaded as `ncg_cumulant_theory (1)(1).pdf`. Sections 2, 4, 6–9 provide the omitted-sector definition, Gaussian tangent, Ward identity, and separation between exact response and efficient closure. Its definitions are preserved here. The coherent compilation, minimum-energy realization, Ward representative, and batched line evaluator are developed in this note.
2. **Supplied preceding review:** `codex_manuscript_review/REVIEW.md`, independent-reference all-row response and the explicit correlation/source-shape correction gap.
3. **Alain Connes**, *Noncommutative Geometry* (1994), Chapters I and IV and Chapter V section 10. Operator algebras of measured spaces, quantized differentials, and projections implementing conditional expectations supply the conceptual context. No theorem there is asserted to give the current numerical coefficient or its computational complexity. Author-hosted text: https://alainconnes.org/wp-content/uploads/book94bigpdf.pdf
4. **Felix Voigtlaender**, *A general version of Price's theorem*, arXiv:1710.03576v2. Rigorous covariance differentiation of Gaussian expectations, including nonsmooth/distributional tests. The mean-direction/Dirichlet identities used here are proved directly. https://arxiv.org/abs/1710.03576
5. **Houman Owhadi and Clint Scovel**, *Conditioning Gaussian measure on Hilbert space*, arXiv:1506.04208v2. Gaussian conditional measures and covariance shorting. Our one-dimensional orthogonal Gaussian disintegration is its elementary finite setting, proved directly. https://arxiv.org/abs/1506.04208

6. **Ivan Nourdin and Giovanni Peccati**, *Stein's method on Wiener chaos*, arXiv:0712.2940v5. Established background for integration by parts connecting Gaussian chaos, the Ornstein--Uhlenbeck operator, and response calculations; no neural compression theorem is attributed to it. https://arxiv.org/abs/0712.2940

Known ingredients such as Gaussian integration by parts, the Dirichlet variational principle, Rao–Blackwellization, and polynomial Gaussian quadrature are not claimed as new discoveries. The contribution developed here is their target-specific arrangement, the closed coherent formulas, and the two implementations with explicit scope and costs.
