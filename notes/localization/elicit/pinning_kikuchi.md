# Pinning, Kikuchi lifts, and deterministic expectations of a quenched ReLU network

# Pinning, Kikuchi lifts, and deterministic expectations of a quenched ReLU network

## Conclusions

There are rigorous bridges between **pinning, weighted simplicial complexes, spectral independence, and down–up walks**; between **Kikuchi free-energy Hessians and subset-indexed spectral algorithms**; and between **graph Kikuchi lifts, token graphs, and exclusion processes**. Those bridges do not, by themselves, give a deterministic numerical approximation to a Gaussian-input neural-network expectation. Their established conclusions are primarily spectral contraction, sampling, detection, recovery, or refutation—not integration. <citations>0,1,2,3,4</citations>

The closest established deterministic approximation framework is **pinning/correlation rounding plus convex moment/entropy hierarchies**, rather than an unmodified Kikuchi matrix. Jain–Koehler–Risteski obtain quantitative mean-field free-energy error bounds and deterministic hierarchy algorithms. Tree-reweighted and hypertree variational constructions give certified one-sided partition-function bounds. Correlation decay and zero-free interpolation give deterministic counting schemes under different structural hypotheses. None of these hypotheses has been established for the gate pushforward law of the network specified here. <citations>5,6,7,8,9,10</citations>

For this network, a viable proposed lift must retain **gate–input or gate–amplitude moments**, not only gate probabilities. It must also distinguish an architecture layer from a simplicial face-size level. A plausible research program combines conditional Gaussian/spherical integration, mixtures of conditional mean-field approximations, and nested convex bounds on amplitude-weighted gate moments. Trickle-down could help certify covariance or residual bounds within this program, but a new observable-specific theorem is needed to turn those bounds into output error.

**Scope and status.** The setting is exactly the supplied one: input dimension and hidden widths 1024, 16 fixed bias-free layers, one quenched He-initialized realization, and randomness only in the Gaussian input. No weights were supplied, so no numerical output vector, realized influence spectrum, or instance-specific error certificate is computed. Literature results below are distinguished from deductions made here and from conjectural transfers. Several recent results were available only as abstracts; where full-text displays were lost in retrieval, only recoverable statements are reported. The survey identifies no established end-to-end deterministic network-expectation theorem combining these tools; this is a search finding, not a proof that no such work exists.

## 1. What pinning and trickle-down actually establish

### Weighted complexes and links

A probability law on assignments can be represented by a multipartite weighted complex. Give coordinate/value pairs their own vertices, represent each supported complete assignment by a facet, weight it by its probability, and take downward closure. A face specifies a partial assignment, and its weighted link represents the conditional law. The corresponding top-face down–up walk is a coordinate heat-bath/Glauber walk. This construction is explicit in Anari–Liu–Oveis Gharan. <citations>11</citations>

For gates, use vertices **(gate index, gate value)**, not just active neurons. This retains zeros, incompatible literals, and support restrictions. Face-size levels count how many gates have been specified; they do not count network depth. Conditioning on a positive-probability face is an exact discrete localization. This last application is a construction for the present problem, not a claim that the complex is an expander.

### Oppenheim’s one-step descent

In the random-walk eigenvalue convention, a useful precise form is:

* the weighted local walk \(P\) is connected;
* every vertex-link walk satisfies \(\lambda_2(P_v)\le a\), with \(0\le a\le 1/2\);
* then

$$
\lambda_2(P)\le \frac{a}{1-a}.
$$

This is Oppenheim’s trickle-down theorem as restated in the matrix trickle-down paper. Connectivity is essential; this is a bound on the second-largest eigenvalue, not automatically on the largest absolute nontrivial eigenvalue. <citations>12</citations>

**Deduction.** Iterating the scalar map gives \(a\mapsto a/(1-ta)\) after \(t\) steps, while the intermediate bounds remain in the applicable range. This shows why a small bound at one high-codimension level need not remain informative after many descents. One cannot repeatedly descend an order-one constant across the entire gate complex and expect a useful uniform bound.

Oppenheim’s original work establishes descent of link spectral gaps and uses these to control higher Laplacians. Its assumptions include connectivity of the relevant links and compatible weighted structures. <citations>13,14</citations>

### Kaufman–Oppenheim: structure of high-order walks

For a pure locally expanding complex, Kaufman–Oppenheim decompose cochains into components associated with different levels. For the upper walk on \(k\)-dimensional faces, the characteristic shrinkage factors are as follows. <citations>15</citations>

$$
\frac{k+1-j}{k+2}+O(\lambda),\qquad 0\le j\le k,
$$

with distinct one-sided and two-sided versions. Two-sided expansion gives approximate eigenspaces and clustering of the spectrum around these values; the one-sided theory gives mixing/contraction without requiring two-sided expansion. The stated formula is the paper’s informal decomposition description, not a uniform constant-free error estimate for arbitrary \(k\). <citations>15,16,17</citations>

This supports a multilevel decomposition of an observable once a suitable weighted complex and expansion hypotheses are available. It does **not** say that replacing a dependent assignment law by a product law has small error, or that truncating an arbitrary observable’s interaction degree is accurate. <citations>16,17</citations>

### Spectral independence

A convenient binary-law convention uses the off-diagonal influence

$$
\Psi_{ij}=\Pr(S_j=1\mid S_i=1)-\Pr(S_j=1\mid S_i=0),\qquad \Psi_{ii}=0,
$$

omitting deterministic coordinates. This equals \(\operatorname{Cov}(S_i,S_j)/\operatorname{Var}(S_i)\), so it is similar to a symmetric normalized covariance matrix with its diagonal removed. These identities follow directly from binary conditional expectation.

Spectral independence requires influence-eigenvalue bounds **throughout the family of conditional laws**, not just a small unconditioned covariance matrix. If a law on \(N\) binary coordinates is \((\eta_0,\ldots,\eta_{N-2})\)-spectrally independent, its assignment complex has the following local expansion parameters. <citations>18</citations>

$$
\left(\frac{\eta_0}{N-1},\frac{\eta_1}{N-2},\ldots,\frac{\eta_{N-2}}1\right).
$$

This is the explicit spectral-independence/local-expansion bridge. Anari–Liu–Oveis Gharan use it with local-to-global walk theorems to prove mixing, including for the hardcore model below uniqueness. Their counting consequence is randomized; they separately identify Weitz’s deterministic counting algorithm as prior work. <citations>19,20,21,22</citations>

For the proposed network lift, the difficult task is not defining this matrix. It is bounding it after feasible pinnings, dealing with coordinates whose conditional variance vanishes, and obtaining conditional probabilities without already solving the integration problem.

### Chen–Eldan localization and entropic independence

Chen–Eldan associate localization martingales of measures with Markov chains. Spectral and entropic independence can be recovered through martingale arguments, without requiring the high-dimensional-expander formalism. Their principal guarantees concern mixing, KL contraction, and functional inequalities. <citations>23</citations>

Anari–Koehler–Vuong formulate trickle-down for linear-tilt localization and apply it to spin-glass and Ising sampling. This is a direct conceptual bridge between the two localization languages, but not a deterministic integration algorithm. <citations>24</citations>

For a law \(\mu\) on \(k\)-element sets, an entropic-independence statement has the form

$$
D(\nu D_{k\to1}\Vert\mu D_{k\to1})
\le \frac{C}{k}D(\nu\Vert\mu),
$$

where the down operator samples a uniformly chosen element. The retrieved paper describes exactly this entropy-contraction notion and proves that spectral independence under arbitrary external fields implies entropic independence. Its applications are modified log-Sobolev inequalities and fast down–up sampling. The external-field hypothesis is stronger than checking only ordinary pinnings. <citations>25</citations>

A recent sparse-localization result weakens the requirement to a restricted family of pinnings fixing at most a fraction \(c\) of the coordinates, with an order-\(c^{-1}\) loss, and obtains entropic stability/independence. This is relevant to avoiding impossible all-pinning requirements, but the available abstract does not provide an output-integration theorem or exact constants. <citations>26</citations>

**Key distinction:** entropy contraction toward a target law is not the same as a bound on that law’s KL distance to a product law. Likewise, a spectral gap says a dynamics converges; it does not make its exponentially large transition matrix or its conditional probabilities cheap to evaluate deterministically.

## 2. Kikuchi matrices, token graphs, and refutation

### The spectral Kikuchi lift

For an even-uniform hypergraph with edge weights \(a_e\), a basic level-\(r\) lift is indexed by \(r\)-subsets:

$$
K_{S,T}=a_{S\triangle T}\quad\text{when }S\triangle T\text{ is a hyperedge},
$$

and zero otherwise. In particular, a \(2q\)-edge exchanges \(q\) occupied vertices with \(q\) unoccupied vertices. The unweighted adjacency definition is explicit in the recent Johnson-sparsification work; weights are the natural linear extension. <citations>27</citations>

Wein–El Alaoui–Moore motivate symmetric-difference matrices through the Hessian of a Kikuchi free energy. An order-\(r\) matrix approximately appears as a block of a level-\(R\) Hessian when \(R\ge r+p/2\) for even tensor order \(p\). The resulting spectral detection and refutation algorithms have rigorous analyses, but the free-energy/Hessian derivation itself is explicitly described as heuristic. These two statuses must not be conflated. <citations>2,28,29,30,31</citations>

Their even-arity random-XOR application gives polynomial-time strong refutation above the \(N^{k/2}\) scale and a continuum of subexponential trade-offs. Strong refutation computes a bound valid for every assignment of the realized instance; its usefulness is high-probability over the random instance. It is not estimation of the uniform assignment average or of a partition function. <citations>32,33,34</citations>

### Guruswami–Kothari–Manohar and Hsieh–Kothari–Mohanty

Guruswami–Kothari–Manohar give an \(N^{O(r)}\)-time refutation trade-off for smoothed Boolean CSPs at a constraint count of order

$$
\widetilde O(N)(N/r)^{k/2-1},
$$

with the stated success probability over the input smoothing. Their work connects semi-random XOR refutation to short even covers in worst-case hypergraphs and resolves the corresponding hypergraph Moore-bound conjecture up to logarithmic losses. <citations>35</citations>

Hsieh–Kothari–Mohanty use a reweighted Kikuchi matrix and an edge-deletion step to simplify the Moore-bound argument, reducing the loss to a single logarithmic factor for uniformity greater than two and recovering the classical graph bound without loss. Their ideas also simplify the smoothed-CSP refutation trade-off. This is a concrete example where **normalizing/reweighting a lift matters**, rather than merely taking a higher tensor power. <citations>36</citations>

A 2026 abstract by Schmidhuber–Hastings claims a normalized Kikuchi hierarchy attaining the sharp trade-off without logarithmic loss for every arity, together with detection, recovery, refutation, and degree-controlled SoS certificates. This update is included as an abstract-level result, not a full-text-verified theorem or a numerical integration guarantee. <citations>37</citations>

### Exact token/exclusion connection

For an ordinary graph, the level-\(r\) Kikuchi graph is precisely the token graph: vertices are \(r\)-subsets; adjacent states differ by moving one token across an edge into an empty vertex. <citations>38,3</citations>

With nonnegative symmetric edge rates, **minus the token-graph Laplacian** is the continuous-time exclusion generator. The diagonal holding-rate term matters: adjacency alone is not the generator. Caputo–Liggett–Richthammer prove equality of interchange and single-particle random-walk spectral gaps; their exclusion-process projection and spectral inclusions yield the same gap for every nontrivial particle-number sector. This is a Laplacian/rate-normalization statement, not equality of adjacency spectra or of discrete-time normalized-walk gaps. <citations>39,40,41</citations>

“Symmetric power” needs a convention: token states forbid collisions, whereas the ordinary bosonic symmetric tensor power includes repeated occupations. Fermionic exterior powers introduce signs. A token adjacency is therefore not automatically either of those operators. This follows directly from their different state spaces and matrix entries.

A particularly relevant recent bridge is Kothari’s result that random even-uniform hypergraph Kikuchi lifts spectrally approximate the complete-hypergraph lifts. Its proof uses blocks of Johnson eigenspaces and band-locality of single-edge operators. This is a genuine subset-lift/spectral-geometry connection, not a trickle-down theorem for arbitrary probability measures. <citations>42,43,44</citations>

### Why these operators do not directly transfer to gates

The neural-network matrices are directed between layers and have signed entries. The number of active gates fluctuates. A gate can change as the input direction moves without another gate switching off. Thus the literal gate process is neither symmetric exclusion nor a fixed-particle-number token walk.

A useful lift can still be built on **subsets of gate labels**, but these subsets should represent retained interactions or test functions, not actual conserved particles. A signed Kikuchi matrix can control an algebraic interaction norm; to become a Markov generator it needs a separate positivity/reversibility construction. Taking absolute values of weights supplies positivity only by changing the problem and potentially destroying cancellations.

## 3. Where certified deterministic approximation is actually known

### Correlation decay: deterministic marginals and counting

Weitz proves a deterministic approximation scheme for the hardcore partition function on graphs of fixed maximum degree \(\Delta\), for

$$
\lambda<\lambda_c(\Delta)
=\frac{(\Delta-1)^{\Delta-1}}{(\Delta-2)^\Delta},\qquad \Delta\ge3.
$$

The method is based on reducing marginal computation to a tree recursion and truncating using correlation decay. For fixed parameters strictly below threshold, the scheme is fully polynomial in input size and requested accuracy. <citations>9,45,20</citations>

**Error mechanism, stated as a general deduction.** If a ratio recursion has a proved root-error bound \(Ce^{-ct}\) at depth \(t\), deterministic enumeration of the truncated recursion gives an absolute marginal error of that size. Partition functions are recovered by a self-reduction product of conditional probabilities. If each required probability is bounded below by \(b>0\), additive errors \(\delta\le b/2\) give log-product error at most \(2N\delta/b\). Thus one needs accuracy of order \(\epsilon/N\), not merely a constant-quality marginal. The constants and contraction metric are model-specific.

Correlation-decay counting also works in some hypergraph regimes where standard strong spatial mixing fails. A retrieved result gives deterministic FPTASs for minimum hyperedge size three/maximum degree six, and for sufficiently large minimum edge size \(k\ge1.66\Delta\). This is a useful warning against treating conventional strong spatial mixing as necessary, but it still requires a specialized recursion analysis. <citations>46</citations>

### Zero-free interpolation: deterministic partition functions

Barvinok-style interpolation approximates \(\log Z\) by a truncated Taylor series along a zero-free domain. Algorithmic usefulness requires both zero-freeness and efficient computation of the retained coefficients. Bounded-degree local enumeration can supply the latter. The retrieved zero-free partition-function work gives deterministic schemes for bounded-degree classical and quantum local systems under its zero-free hypotheses. <citations>47,48,49,50</citations>

A concrete quantitative example is the paper combining interpolation and correlation-decay recursions: for a graph with minimum degree at least three, it approximates the number of sink-free orientations within a factor \(e^\epsilon\), deterministically, in time \(O(N(M/\epsilon)^7)\), where \(N\) and \(M\) are vertex and edge counts. <citations>51</citations>

**Elementary remainder bound.** If \(Z(z)=Z(0)\prod_{i=1}^d(1-z/\zeta_i)\), all roots satisfy \(|\zeta_i|\ge R>|z|\), and \(q=|z|/R<1\), truncating \(\log Z\) after degree \(t\) has error at most

$$
\frac{d q^{t+1}}{(t+1)(1-q)}.
$$

This follows by summing the logarithm’s geometric-series tail. It gives an explicit deterministic error once roots are excluded and coefficients are computable. It supplies neither fact for a gate-generated model automatically.

### Pinning, correlation rounding, and mean-field free energy

Jain–Koehler–Risteski establish, for Ising models, the mean-field gap

$$
0\le \log Z-F_{\mathrm{MF}}
\le O\!\left((N\|J\|_F)^{2/3}\right),
$$

in their interaction normalization. In particular, \(\|J\|_F=o(\sqrt N)\) implies an \(o(N)\) free-energy error. Their convex-hierarchy/correlation-rounding algorithms are deterministic before the paper’s separate random-subsampling section. <citations>5,52</citations>

The operative idea is to condition on a small coordinate set and replace the remaining law by a conditional product approximation. The paper explicitly relates this to statistical-physics pinning. The low-order moments needed for conditioning and rounding are supplied by a convex relaxation, rather than assumed known. <citations>53</citations>

This is the closest known structural template to the proposed network program. However, an \(o(N)\) error in \(\log Z\) may still mean an exponentially large relative error in \(Z\), and says nothing immediate about an individual output coordinate. Moreover, dense He weights do not identify a gate Ising interaction matrix \(J\), let alone verify its Frobenius-norm regime.

### Certified variational bounds versus ordinary Kikuchi estimates

For a discrete model, the exact variational identity is

$$
\log Z=\sup_q\{\mathbb E_q H+\mathcal H(q)\}.
$$

Any explicit tractable family of genuine probability laws gives a lower bound. Tree-reweighted variational inference gives upper bounds from convex combinations of tractable structures. Wainwright–Jaakkola–Willsky prove the convex tree bound and describe nested hypertree extensions yielding progressively tighter bounds, related to convexified Kikuchi free energies. Optimizing over hypertrees brings its own hardness. <citations>7,8</citations>

Ordinary Kikuchi uses locally consistent region beliefs and an inclusion–exclusion entropy. A global optimum of that surrogate is not necessarily close to the true answer; adding regions is not a universal monotonic-improvement theorem. The graph-cover characterization explicitly explains why Bethe and Kikuchi partition estimates can differ from the exact partition function. Concavity results for reweighted Kikuchi concern optimization of the surrogate, not its approximation error. <citations>54,55</citations>

Consequently, “systematically improvable with error bars” should mean **nested valid lower and upper bounds**, not just a sequence of more elaborate generalized-belief-propagation fixed points.

### A deterministic continuous-integration precedent

Gamarnik–Smedira’s *Integrating High-Dimensional Functions Deterministically* gives quasi-polynomial deterministic integration for products of weakly varying differentiable local potentials on a rectangular domain. It discretizes the integral and uses correlation decay. Its bounded-degree, bounded-gradient, near-constant-potential conditions are restrictive; hard gate constraints and the dense deep network do not directly satisfy them. <citations>56,57,58,59</citations>

This is a positive deterministic-integration result in a nearby language, not a ready-made algorithm for the specified expectation.

### Turning log-partition bounds into expectation bounds

**Deduction useful across these methods.** For a bounded observable \(h\), set

$$
F(t)=\log\mathbb E_\mu e^{t h};\qquad F'(0)=\mathbb E_\mu h.
$$

If certified bounds \(L(t)\le F(t)\le U(t)\) are available, convexity gives, for \(t>0\),

$$
\frac{L(0)-U(-t)}t
\le \mathbb E h
\le \frac{U(t)-L(0)}t.
$$

Here \(F(0)=0\) is exactly known. Small free-energy gaps are useful only relative to the finite-difference step; an additive error \(\epsilon\) typically contributes \(\epsilon/t\). For \(0\le h\le B\), \(F''(t)\le B^2/4\), supplying an explicit curvature remainder. This is a rigorous route from variational partition bounds to expectations, unlike differentiating an uncontrolled Kikuchi optimum.

## 4. Quantization and coherent states: an analytic connection, not a shortcut

Berezin quantization maps classical symbols to operators; covariant/lower symbols map operators back to coherent-state expectations. The retrieved work on operator convolutions provides a rigorous framework for these representations and Berezin–Lieb inequalities. <citations>60,61,62</citations>

A precise finite-dimensional version can be derived as follows. Let normalized states \(|z\rangle\) resolve the identity,

$$
\int |z\rangle\langle z|\,d\nu(z)=I,
$$

and suppose a self-adjoint operator has a real upper symbol,

$$
A=\int a_{\rm up}(z)|z\rangle\langle z|\,d\nu(z),
\qquad a_{\rm low}(z)=\langle z|A|z\rangle.
$$

Scalar Jensen inequalities in the eigenbasis give, for a convex function \(\phi\),

$$
\int\phi(a_{\rm low})\,d\nu
\le \operatorname{Tr}\phi(A)
\le \int\phi(a_{\rm up})\,d\nu.
$$

Choosing \(\phi(t)=e^{-\beta t}\) sandwiches an **operator partition function** between classical-symbol integrals. If the two symbols differ uniformly by at most \(\delta\), their bounding log integrals differ by at most \(|\beta|\delta\). These are analytic deterministic bounds; numerical certification still requires evaluating the integrals and bounding the symbol discrepancy. The argument above states its own normalization and assumptions because the retrieved full-text displays did not preserve all formulas. The generalized Berezin–Lieb result is supported by the source. <citations>63</citations>

For the network, no canonical coherent-state resolution, quantization parameter, or symbol estimate has been identified. Gaussian coherent-state phase space, a real input sphere, a Boolean gate cube, and a fixed-cardinality subset space are different objects. Choosing an operator whose trace equals the desired expectation is possible in many artificial ways; it helps only if its symbols are simpler to integrate and the approximation gap is computable.

A potentially useful analogy is that conditional-expectation or feature-compression operators smooth an observable, while an amplitude-aware moment lift retains increasingly rich features. But treating that compression as “Toeplitz quantization” does not supply a semiclassical error bound. For indicator gates, uniform symbol-smoothing errors are especially problematic at boundaries; the continuous final ReLU output is a better target for approximation than the discontinuous gates themselves.

## 5. Exact structure of the specified network

The results in this section are direct deductions from the supplied network definition; they do not rely on an annealed random-weight approximation.

### Radial reduction

Write \(X=RU\), where \(U\) is uniform on \(S^{n-1}\), \(R\) is independent with the chi distribution, and define \(f(u)=x_L(u)\). Positive homogeneity gives

$$
\mathbb E x_L(X)=c_n\mathbb E_U f(U),
\qquad
c_n=\sqrt2\,\frac{\Gamma((n+1)/2)}{\Gamma(n/2)}.
$$

This removes the unbounded radial variable exactly. It does not remove the high-dimensional angular integral.

Let

$$
K=\prod_{\ell=0}^{L-1}\|W_\ell\|_{\rm op}.
$$

Every output coordinate satisfies \(0\le f_j(u)\le K\), and the vector map is \(K\)-Lipschitz. This worst-case constant is computable from the actual weights, but can be too large to provide practical error bars. Sharper restricted-Jacobian bounds are desirable.

### A full gate pattern defines a cone; a partial pattern generally does not

For a complete gate assignment \(s\), put \(D_\ell(s)=\operatorname{diag}(s_\ell)\) and

$$
A_s=W_0D_0(s)W_1D_1(s)\cdots W_{L-1}D_{L-1}(s).
$$

On its feasible region \(C_s\),

$$
x_L(X)=XA_s.
$$

All preactivation inequalities become linear in \(X\) once the preceding gate assignments are fixed; hence \(C_s\) is a polyhedral cone, with strict/non-strict inequalities matching the gate convention. Some patterns are infeasible or have zero Gaussian probability. Boundary conventions matter when an entire previous layer is zero, because downstream preactivations can then vanish on a positive-probability set.

A first-layer pin is a half-space. A later-layer pin with preceding gates unpinned is generally a **union of such cone pieces**, potentially nonconvex. Therefore continuous log-concave localization or truncated-Gaussian formulas for a single convex polyhedron cannot be applied to arbitrary deep partial pinnings without further decomposition.

The exact target is

$$
\mathbb E x_L(X)
=\sum_s \left(\mathbb E[X\mathbf1_{C_s}]
ight)A_s.
$$

Cone probabilities alone do not provide the vector first moments appearing here.

### Path expansion identifies the missing moments

For output \(j\), each layered path contributes a product of weights, one input coordinate, and one gate per layer:

$$
\mathbb E x_{L,j}
=\sum_{i_0,\ldots,i_{L-1}}
\left(\prod_{\ell=0}^{L-1}W_{\ell,i_\ell i_{\ell+1}}\right)
\mathbb E\!\left[X_{i_0}\prod_{\ell=0}^{L-1}s_{\ell,i_{\ell+1}}\right],
\qquad i_L=j.
$$

Thus a gate-only law is insufficient. Moreover, replacing gates by independent Bernoulli variables **independent of \(X\)** yields zero for the path expansion, since \(\mathbb E X=0\). Except for a degenerate zero network, that cannot be the desired nonnegative activation mean. “Independent-gate mean field” needs a more careful definition.

### Exact first-layer checks

For a nonzero first-layer column \(w_j\),

$$
\mathbb E\operatorname{ReLU}(Xw_j)=\frac{\|w_j\|}{\sqrt{2\pi}},
\qquad
\mathbb E[X\mathbf1_{Xw_j>0}]
=\frac{w_j^\top}{\|w_j\|\sqrt{2\pi}}.
$$

For two columns with correlation \(\rho\),

$$
\Pr(s_i=s_j=1)=\frac14+\frac{\arcsin\rho}{2\pi},
$$

and, writing \(\theta=\arccos\rho\),

$$
\mathbb E[x_{1,i}x_{1,j}]
=\frac{\|w_i\|\|w_j\|}{2\pi}
\left(\sin\theta+(\pi-\theta)\cos\theta\right).
$$

These are exact Gaussian two-dimensional integrations, derivable by rotating to the span of the columns. They supply deterministic initial moments and tests for any proposed implementation, without implying Gaussianity of deeper preactivations.

### Quenched means invalidate a naive zero-mean closure

For later layers,

$$
\mathbb E_X[x_\ell W_\ell]=(\mathbb E_Xx_\ell)W_\ell,
$$

which generally is not zero for the fixed matrix. Conditional on the preceding network, if a *new* He column were averaged over its initialization, this mean would have variance \(2\|\mathbb E_Xx_\ell\|^2/n\). That identity is a diagnostic across possible weight draws, not permission to average the already-fixed weights.

In particular, nonzero activation means need not produce vanishing preactivation means as width grows. Deep gates are not automatically fair coins in the quenched problem. A Gaussian moment closure should propagate the realized mean as well as covariance. If a scalar surrogate preactivation is \(N(m,v)\), its rectified mean is

$$
\sqrt v\,\varphi(m/\sqrt v)+m\Phi(m/\sqrt v),
$$

with the deterministic \(v=0\) limit. This is an exact formula for the surrogate, not an error bound for replacing the true preactivation law by that surrogate.

## 6. A defensible hierarchy: proposals and conditional guarantees

### Proposal A: pinning plus mixtures of conditional approximations

For a chosen gate set \(T\), retain a mixture over its feasible assignments:

$$
\mathbb E f(U)=\sum_a p_T(a)\,\mathbb E[f(U)\mid S_T=a].
$$

Approximate **within each branch**, retaining conditional amplitudes, rather than replacing the full law by one product distribution. This mirrors correlation rounding’s conditional-product strategy, while acknowledging that branch laws can be unions of cones. The analogy to established pinning/correlation rounding is structural, not a transferred theorem. <citations>53,64</citations>

A simple rigorous certificate follows from the mixture identity. If branch weights have total \(\ell_1\) error at most \(\delta\), and each coordinate’s conditional-mean error is at most \(\epsilon_a\), then

$$
|\mu_j-\widehat\mu_j|
\le c_n\left(K\delta+\sum_a p_T(a)\epsilon_a\right).
$$

A discarded set of branches with true total mass at most \(\tau\) contributes at most \(c_nK\tau\) per coordinate. These are derived bookkeeping bounds. They require certified branch masses and conditional means, not sample estimates presented as deterministic certificates.

**Hypothesis:** adaptive pinning can expose a modest number of branches where amplitude-aware dependencies are small enough for useful conditional mean-field bounds. This must be proved or certified on the realized matrices; it is not implied by initialization alone.

### Proposal B: amplitude-aware region and moment lifts

Maintain, for selected gate subsets \(A\), quantities such as

$$
p_A=\mathbb E\prod_{g\in A}s_g,
\qquad
v_{A,i}=\mathbb E\left[U_i\prod_{g\in A}s_g\right].
$$

For overlap tests and spectral compression, use Gram entries

$$
M_{(A,i),(B,j)}
=\mathbb E\left[U_iU_j\prod_{g\in A\cup B}s_g\right].
$$

This matrix is positive semidefinite by construction. With uncentered \(\{\pm1\}\)-valued gate characters, products instead reduce through symmetric difference, giving an algebraic resemblance to Kikuchi matrices. But the law-dependent weights and amplitude factors are different from random-XOR clause weights.

The path identity shows that mixed moments involving one input coordinate and at most \(L\) gates suffice algebraically for exact output means. Enumerating them all is prohibitive. Gate-only moments through degree \(L\) do not suffice, and a moment PSD constraint alone does not enforce the network’s geometric gate semantics.

A certified relaxation should include:

* consistency and normalization of region probabilities;
* deterministic gate support constraints from the realized network;
* positive moment/localizing matrices;
* spherical constraints such as \(\sum_iU_i^2=1\);
* certified intervals for any moments supplied by deterministic integration;
* a bounded-output observable or a uniformly certified surrogate for it.

Optimize the output linear functional over a nested outer feasible family to obtain upper and lower bounds. Adding valid constraints shrinks the interval. At the full exact model this is exact by construction, but neither practical convergence rates nor polynomial complexity are established here. An ordinary Kikuchi entropy replacement does not deliver this nesting automatically.

### Proposal C: output-specific approximation rather than global independence

Let \(h_j(U,S)\) be a retained-feature approximation to \(f_j(U)\). A certified uniform residual

$$
\|f_j-h_j\|_\infty\le\epsilon_j
$$

implies output error at most \(c_n\epsilon_j\) once \(\mathbb Eh_j\) is integrated exactly. More generally, a certified \(L^2\) residual gives the same form with its \(L^2\) norm, by Cauchy–Schwarz. Moment intervals supply an additional explicit coefficient-weighted integration error.

This is where a weighted subset lift could help: estimate the high-order residual relevant to the output rather than demand that the entire gate law be close to a product measure. Signed path coefficients make absolute-sum bounds extremely pessimistic, so a PSD/Gram or operator-norm certificate is preferable when available.

**A necessary warning, proved by example:** pairwise independence alone does not control high-degree observables. Under a uniform even-parity law on sufficiently many bits, low-order marginals can agree with the independent law while the full parity expectation differs maximally. This example does not refute all-pinning spectral independence; it refutes the proposed shortcut from small unconditioned pair correlations to accurate path-product estimates.

### A conditional error theorem for Gaussian/independent-coordinate closure

A useful derived target for certification is a transport-error recurrence. Let \(\mu_\ell\) be the true activation law and \(\nu_\ell\) a deterministically represented surrogate. Form the linear pushforward \(\widetilde\nu_{\ell+1}=(z\mapsto zW_\ell)_\#\nu_\ell\), replace it by a tractable Gaussian or independent-coordinate law \(G_{\ell+1}\), and set \(\nu_{\ell+1}=\operatorname{ReLU}_\#G_{\ell+1}\). If the replacement has a certified Wasserstein-1 error \(\rho_\ell\), using Euclidean transport cost, then

$$
e_{\ell+1}\le \|W_\ell\|_{\rm op}e_\ell+\rho_\ell,
\qquad e_\ell=W_1(\mu_\ell,\nu_\ell).
$$

Starting from the exact Gaussian input gives

$$
\|\mathbb E x_L-\mathbb E_{\nu_L}x\|_2
\le\sum_{\ell=0}^{L-1}\rho_\ell
\prod_{r=\ell+1}^{L-1}\|W_r\|_{\rm op}.
$$

The proof uses transport contraction under a linear map, the 1-Lipschitz property of ReLU, the triangle inequality, and a coupling bound for vector means. Matching means and covariances does not certify \(\rho_\ell\); nor does a gate influence bound alone. Pinning and region moments could be used to bound the replacement residual branch by branch. This identifies a precise missing theorem and the downstream amplification that any claimed mean-field guarantee must address.

### What trickle-down would need to add

A credible theorem would have four parts:

1. Define conditional gate/amplitude link operators and verify their connectivity or handle support components.
2. Prove or certify appropriately normalized conditional covariance bounds for the relevant pinnings.
3. Use matrix/localization trickle-down to propagate those bounds to the retained hierarchy levels.
4. Prove an **observable-specific approximation theorem** converting the resulting bounds into a computable mean-field or truncation error for \(f_j\).

Parts 1–3 resemble established spectral-independence and matrix trickle-down machinery. Part 4 is the missing network theorem. The matrix trickle-down coloring work illustrates induction on conditioned combinatorial models, not this amplitude-weighted observable transfer. <citations>65,66,67</citations>

Entropy-based certification offers another route. If a genuine approximate joint law \(q\) on angular inputs and gates satisfies a certified KL bound \(D(q\Vert\mu)\le\eta\), then boundedness and Pinsker’s inequality give

$$
|\mathbb E_q f_j-\mathbb E_\mu f_j|
\le K\sqrt{\eta/2}.
$$

Multiply by \(c_n\) for the Gaussian target. But a discrete gate KL bound cannot be substituted for this joint bound unless the conditional amplitude law is also controlled. Independently randomized gates may give positive probability to impossible input/gate pairs, making this KL direction infinite. Functional inequalities alone do not supply the needed \(\eta\).

### Hypotheses worth testing, not conclusions

**H1 — Conditional amplitude decorrelation.** After a judicious small family of gate pinnings, the output-relevant mixed moments have a rapidly decaying interaction residual.

**H2 — Weighted rather than scalar descent.** Influence bounds weighted by downstream sensitivity are substantially sharper than uniform all-gate bounds, allowing useful descent despite redundant gates, near-deterministic branches, and dense connectivity.

**H3 — Certifiable hierarchy gap.** A nested amplitude-aware region/moment relaxation has an empirically small and eventually provably shrinking upper–lower gap for typical fixed He realizations.

**H4 — Quantization only after a concrete construction.** A feature-compression/coherent-state formulation becomes useful if its upper/lower symbols and resolution measure can be specified and their integration gap bounded. Until then it is a guiding analogy, not a certification method.

Each hypothesis is falsifiable. A result holding with high probability over initialization would still need to state its failure probability and distinguish a typical-realization theorem from an a posteriori certificate for the supplied realization.

## 7. An unconditional deterministic baseline with an explicit error bound

A deterministic integration guarantee exists for every fixed finite network, without expansion or independence. It is generally astronomically expensive.

**Construction and proof.** Partition \([-T,T]^n\) into boxes of side at most \(h\). At each center \(z_b\), evaluate the network. Weight that vector by the exact Gaussian box probability, a product of one-dimensional Gaussian-CDF differences:

$$
Q=\sum_b\Pr(X\in b)x_L(z_b).
$$

The Lipschitz bound gives interior Euclidean error at most \(Kh\sqrt n/2\). Outside the cube, \(\|x_L(X)\|\le K\|X\|\). Cauchy–Schwarz and a Gaussian union-tail bound give

$$
\|\mathbb E x_L(X)-Q\|_2
\le \frac{Kh\sqrt n}{2}
+K\sqrt{2n^2e^{-T^2/2}}.
$$

Certified CDF and arithmetic intervals add their own computable numerical-error term. Choosing \(T\) and \(h\) makes this an arbitrary-accuracy deterministic vector approximation, with roughly \((2T/h)^n\) evaluations. This is a derived baseline, not a claim of feasibility at the specified dimension.

An angular partition can similarly exploit the exact radial reduction, provided its cell probabilities are certified. The purpose of pinning and a Kikuchi/region hierarchy would be to replace exponential global cubature by a much smaller certified calculation—not to establish that deterministic integration is possible at all.

## 8. Recommended mathematical and computational program

Start with exact first-layer means, pair moments, and Gaussian gate correlations. Propagate a **quenched mean-and-covariance Gaussian closure** as a clearly labeled deterministic heuristic. Never impose fair deep gates or independence from the input.

Next, develop small pinned-region or mixed-moment calculations with certified integration intervals. Select regions by downstream sensitivity and by influence diagnostics, while treating diagnostics as evidence rather than proof. Compare these against the heuristic and against a deterministic coarse integration certificate on reduced-dimensional or small-width test networks.

Build nested output upper/lower bounds before attempting a universal trickle-down theorem. The valuable deliverable is a shrinking interval for each output coordinate. A useful theorem would explain why those intervals shrink cheaply for typical fixed networks, and what instance-specific quantities certify this behavior.

The most defensible synthesis is therefore:

**Gaussian/spherical geometry supplies the measure and amplitude moments; pinning supplies conditional branches; correlation rounding supplies the conditional-product template; convex region/moment hierarchies supply monotone certificates; trickle-down may supply covariance control; Kikuchi/token lifts supply an interaction-indexed operator language; Berezin–Lieb supplies a possible analytic sandwich only after a genuine symbol construction.**

The existing literature establishes important pieces of this chain, but not the chain itself for the specified network. A low-level deterministic method with nonvacuous, realized-weight error bars remains a research target—not a consequence of current mixing or random-XOR refutation theorems.
