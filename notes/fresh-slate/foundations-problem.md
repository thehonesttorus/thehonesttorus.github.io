# Foundations: the problem from first principles

### What any estimator of the per-neuron means of a deep He-initialised ReLU network must know, how accurately, at what cost, and how this scales

*Fresh-slate foundations note, Part 1 (2026-10-01). Written for the six design streams. It derives the problem's exact identities and symmetries, the information an estimator needs about each layer's law and to what accuracy, the cost floors at n = 1024, and the scalings in width and depth. It makes no reference to the internals of any existing estimator. Inputs used: [BRIEF.md](BRIEF.md) §§1 and 3 (task, metering, neutral facts), the measured prices of [../streams/costmodel/REPORT.md](../streams/costmodel/REPORT.md), and the starter-kit reference numbers quoted in [../competition-phase2.md](../competition-phase2.md). All other numbers were computed for this note (Monte Carlo on He networks of widths 16–256 and on one width-1024 network of the shared bench, depth 16; methods and code in Appendix A).*

Sibling note: [foundations-unlocks.md](foundations-unlocks.md) (Part 2, the theoretical unlocks) derives several of the same identities independently (its unlocks 1 and 3–5); where both notes state one, the agreement is recorded inline, as the programme's rule C3 asks.

Labels. **Exact**: an identity proved here (proof or proof sketch given). **Measured**: a number from a run made for this note. **Fact**: taken from the cited source. **Estimate**: an extrapolation, with its basis stated. Citations of textbook results from memory are marked *(from memory)*.

---

## 0. Summary

Five statements carry the note; each is derived in the section named and checked numerically.

1. **The answer is a pairing of the input law with a piecewise-linear observable, and every sensitivity is a face statistic** (Sections 2.1, 2.6, 3.1). Exactly, $\mathbb E[a_L]=\sum_\sigma P(C_\sigma)\,b_\sigma M_\sigma$ (linear in the Gaussian barycentres $b_\sigma$ of the activation cones) $=\sum_{\text{walls}}p_{z}(0)\,\mathbb E[\|\nabla z\|^2\cdot\text{downstream gain}\mid z=0]$ (Euler + Stein; checked to 1.2 % at n = 64). The exact first-order sensitivities of a neuron's mean to its own cumulants are $\mathbb E[\mathrm{relu}^{(k)}(z)]/k!$: the gate probability, half the density at the wall, and the slope and curvature of that density. The first layer is an exactly known Gaussian ($W_1$ enters only through $W_1^\top W_1$), and positive homogeneity makes the overall scale a neutral mode: a relative scale error at any layer reaches the output with gain exactly 1.

2. **The readout needs each final neuron's $(m,s^2,\lambda_3,\lambda_4)$ and nothing more** (3.1–3.2). Given the exact quenched $(m_i,s_i^2)$, the Gaussian readout errs by $4.2\times10^{-4}$ rms at n = 1024 (MSE $1.7\times10^{-7}$, 17 times the bar); with $\lambda_3,\lambda_4$ the residual is ≤ $1.7\times10^{-5}$ (MSE ≤ $3\times10^{-10}$). Tolerances for MSE 2.5·10⁻⁹ from each source: the mean to $5.6\times10^{-5}$ relative, the variance to 0.14 %, the $\lambda_3$ term to ≈ 12 % and the $\lambda_4$ term to ≈ 40 % (measured at n = 1024). The non-Gaussianity grows linearly in depth, $\mathrm{rms}\,\lambda_3\approx(5$–$7)\,l/n$ and $\lambda_4\approx(3.6$–$6)\,l/n$ (0.10 and 0.056 at the last layer at n = 1024), and is to a large part a common scale mixture ($\lambda_{3,i}\approx b_lt_i$ with $R^2$ up to 0.89, $\lambda_4$ independent of $t$), the imprint of the norm process that positive homogeneity factors out.

3. **One layer back, errors enter through Isserlis pairings: trace channels are not averaged, and every index pattern counts** (3.3–3.4). With $W_{l+1}$ independent of the state and of its errors, $\kappa_k(z_{l+1,i})=K_k(a_l)[w_i^{\otimes k}]$: a bias in the average per-neuron variance passes undiminished (so it must be ≤ 0.13 %), while spread errors are averaged by the $n$ inputs ($\Sigma_{L-1}$ to ≈ 0.6 % in Frobenius norm). Every index-coincidence pattern of $K_3$ and $K_4$ carries leading-order weight at every width measured, so the all-distinct parts ($n^3$, $n^4$ entries) can be neither dropped nor stored; the neuron-averaged fourth cumulant is one scalar, the thin-shell excess $3\sigma^4[\mathrm{Var}\|X\|^2-2\|\Sigma\|_F^2]$ (checked to 0.1–2 % at n = 1024). Through the gate, $\Sigma_l$ needs the pairwise $(2,1)$, $(3,1)$, $(2,2)$ cumulant slices of $z_l$ (per-neuron cumulants do not help it at all): predicting $\Sigma_l$ from the exact pairwise Gaussian part alone induces next-layer errors of ≈ $(3$–$8)\times10^{-4}$ at n = 1024, adding the third-order slices ≈ $(5$–$7)\times10^{-5}$, adding the fourth-order slices ≈ $2\times10^{-5}$, against a per-layer budget of ≈ 4·10⁻⁵.

4. **In depth, scale errors never decay, coherent second-order errors hardly decay, random mean errors do** (3.5). The dilation mode has gain 1; the expected Laplacian of the pulled-back readout (the response to a pure covariance error) is $O(1)$ at every depth, so a coherent relative variance error must stay below ≈ 3·10⁻⁴ at shallow and ≈ 10⁻³ at deep layers, and below ≈ 10⁻⁴ if it has the same sign at all layers (n = 64); random mean errors are damped by the derivative of the arc-cosine correlation map, $\prod_kf'(\rho_k)$ = 0.013 from the input to 0.87 from layer 15 (the measured gains exceed it by at most a factor 2.4), and accumulate to ≈ ×2.4, so every layer's means are needed to ≈ 4·10⁻⁵ rms. The bar sits at $\approx0.4\,(L/n)^2$: every first-order ($1/n$) effect must be computed to ≈ 10–40 %.

5. **Second order is affordable and unavoidable; everything beyond it must be generated, not stored** (Section 4). At the 0.1 B floor a layer has 6.4 units (≈ 6 dense or 11 Strassen products) and ≈ 1,100 metered calls. Exact per-neuron quadratic forms with a dense covariance (0.56–1 unit per layer) are a floor: no compression of $\Sigma$ that ignores $W_{l+1}$ reaches 0.1 % below rank 0.47 n (last layer) to 0.93 n (second layer), measured at n = 1024. Carrying one $n\times n$ object costs 1–2 units per layer; one contraction pass of an explicit third-order tensor costs ≥ 1024 units ($=B$); plain sampling needs $7.5\times10^6$ samples (114 B) for raw $10^{-8}$. A design therefore gets the covariance level for 16–32 % of the floor budget and has about four or five $n\times n$-equivalents per layer for all higher-order content, which must be produced from $O(n^2)$ objects at $O(n^3)$ cost.

---

## 1. The object, and the accuracy the score asks for

**Network.** $x\sim\gamma_n=\mathcal N(0,I_n)$, $n=1024$. Row-vector convention (the estimator receives $W_l$ in the `x @ W` layout):
$$z_1=xW_1,\qquad a_l=\mathrm{relu}(z_l),\qquad z_{l+1}=a_lW_{l+1},\qquad l=1,\dots,L-1,\quad L=16,$$
with $W_l$ i.i.d. $\mathcal N(0,2/n)$ entries, independent across $l$, no biases. The scored vector is $\mu_L=\mathbb E[a_L]\in\mathbb R^n$ (all layers are returned; only the last is scored).

**Arrows.** In the programme's language an arrow is a linear layer with the ReLU immediately before it: $A_{l+1}(\zeta)=\mathrm{relu}(\zeta)\,W_{l+1}$, so $z_{l+1}=A_{l+1}(z_l)$ and $z_L=A_L\circ\cdots\circ A_2(xW_1)$. The first map $x\mapsto xW_1$ is an arrow whose gate is trivial (the input is not a post-activation), and the readout $\zeta\mapsto\mathrm{relu}(\zeta)$ is the gate half of a final arrow with identity linear part. Throughout, $p_l$ is the law of $z_l$ and $P_l$ the law of $a_l$; $\mu_l=\mathbb E a_l$, $\Sigma_l=\mathrm{Cov}(a_l)$, $K_k(a_l)$ the $k$-th joint cumulant tensor of $a_l$.

**Faces.** The activation pattern of layer $l$ is $\sigma_l(x)=\{j: z_{l,j}(x)>0\}$; a history $\sigma=(\sigma_1,\dots,\sigma_L)$ is constant on an open polyhedral cone $C_\sigma$ (Section 2.1). Walls are the codimension-one pieces $\{z_{l,j}=0\}$.

**The accuracy target.** Raw MSE $\approx10^{-8}$ (and below) is a root-mean-square per-neuron error of $10^{-4}$ ($7\times10^{-5}$ for $5\times10^{-9}$). The final-layer means have root-mean-square $\approx0.95$ (zero-prediction MSE 0.9095 on the public mini split, **Fact**, competition-phase2.md) and per-neuron standard deviation $\approx0.27$ (mean per-neuron variance 0.0748, **Fact**). So the target is a relative accuracy of $\approx10^{-4}$ on $O(1)$ quantities, and $\approx4\times10^{-4}$ of the per-neuron spread. For orientation (**Fact**, BRIEF §1 and bench/RESULTS.md, six width-1024 networks): a Gaussian covariance closure scores raw $4.3\times10^{-6}$ and plain sampling at the 0.1 floor $1.1\times10^{-5}$; the target is ≈ 400 times below the former.

---

## 2. Exact identities and symmetries (part a)

### 2.1 Positive homogeneity: the sphere reduction, faces as cones, the barycentric and wall formulas

Every $z_l$ and $a_l$ is positively homogeneous of degree one: $a_l(tx)=t\,a_l(x)$ for $t\ge0$.

**(i) Sphere reduction (Exact).** Write $x=ru$, $r=\|x\|\sim\chi_n$, $u$ uniform on $S^{n-1}$, independent. For every statistic homogeneous of degree $k$,
$$\mathbb E\big[a_{l_1,i_1}\cdots a_{l_k,i_k}\big]=\mathbb E[r^k]\int_{S^{n-1}}a_{l_1,i_1}(u)\cdots a_{l_k,i_k}(u)\,d\sigma(u),\qquad \mathbb E[a_L]=\mathbb E[r]\int_{S^{n-1}}a_L\,d\sigma,$$
with $\mathbb E[r]=\sqrt2\,\Gamma(\tfrac{n+1}2)/\Gamma(\tfrac n2)=31.992188$ at $n=1024$ ($\approx\sqrt{n-1/2}$) and $\mathbb E[r^2]=n$. The joint law of all layers is a scale mixture, $r\times$(the network on the sphere), with $r$ independent. Two consequences:
- the radial fluctuation puts the rank-one term $\frac{\mathrm{Var}\,r}{(\mathbb Er)^2}\mu_l\mu_l^{\top}=4.884\times10^{-4}\,\mu_l\mu_l^{\top}\ (\approx\mu_l\mu_l^\top/2n)$ into every $\Sigma_l$ — the purely radial part of the spike along the mean direction (its Rayleigh quotient along $\hat\mu_l$ is $\|\mu_l\|^2/2n$, about half the per-neuron mean square, against an average eigenvalue equal to the per-neuron variance; at depth the ratio is $\approx 6$ from the radius alone, Section 5);
- the problem is a problem about the uniform measure on the sphere and the piecewise-linear map on it; Gaussian and spherical statements differ only by moments of $r$.

**(ii) Faces are cones; the Jacobian is face-constant (Exact).** On $C_\sigma$ every pre-activation is a fixed linear functional of $x$, and $a_L(x)=x\,M_\sigma$ with $M_\sigma=W_1D_1W_2D_2\cdots W_LD_L$, $D_l=\mathrm{diag}(1_{\sigma_l})$. The Jacobian $J(x)=\partial a_L/\partial x=M_{\sigma(x)}$ is homogeneous of degree zero, and Euler's identity reads $a_L(x)=x\,J(x)$.

**(iii) Barycentric formula (Exact).** Summing over faces,
$$\mathbb E[a_L]=\sum_\sigma\mathbb E[x\,1_{C_\sigma}]\,M_\sigma=\sum_\sigma P(C_\sigma)\,b_\sigma\,M_\sigma,\qquad b_\sigma=\mathbb E[x\mid x\in C_\sigma]=\frac{\mathbb E[r]}{P(C_\sigma)}\int_{C_\sigma\cap S^{n-1}}u\,d\sigma(u).$$
The scored vector is linear in the face barycentres; the arrows act on a face linearly; all of the nonlinearity is in which faces carry mass and where their barycentres sit. (The Gaussian barycentre of a cone is $\mathbb E r$ times its spherical barycentre.) The same statement holds at every layer for the orthant faces of $z_l$: $\mathbb E[a_l]=\sum_{S\subseteq[n]}P(\sigma_l=S)\,\mathbb E[z_l\mid\sigma_l=S]\,D_S$. (Also unlock 1 of foundations-unlocks.md.) This is the concrete content, in this problem, of "objects are barycentres": the mean of the next layer is the mass-weighted sum of linear images of face barycentres. (It is a realisation, not a definition for the abstract theory.)

**(iv) Wall formula (Exact).** Euler gives $a_{L,i}=x\cdot\nabla a_{L,i}$, and Gaussian integration by parts (Stein's identity, valid for the piecewise-constant gradient in the distributional sense) gives $\mathbb E[x\cdot\nabla f]=\mathbb E[\Delta f]$. The distributional Laplacian of a stack of ReLU layers is a sum of single-wall terms, so
$$\boxed{\;\mathbb E[a_{L,i}]=\sum_{l=1}^{L}\sum_{j=1}^{n}\mathbb E\Big[\delta(z_{l,j})\,\|\nabla_x z_{l,j}\|^2\,\frac{\partial a_{L,i}}{\partial a_{l,j}}\Big]=\sum_{l,j}p_{z_{l,j}}(0)\;\mathbb E\Big[\|\nabla_xz_{l,j}\|^2\frac{\partial a_{L,i}}{\partial a_{l,j}}\,\Big|\,z_{l,j}=0\Big].\;}$$
Each wall contributes the density of its pre-activation at the threshold times the conditional mean, on the wall, of its squared normal speed times the downstream gain. For $L=1$ this is $\mathbb E\,\mathrm{relu}(\langle w,x\rangle)=\|w\|^2\varphi(0)/\|w\|=\|w\|/\sqrt{2\pi}$. **Measured** (n = 64, L = 16, 6·10⁵ samples, symmetric second differences of $a_L$ under isotropic input noise): $\mathbb E[\Delta_xa_L]$ equals $\mathbb E[a_L]$ neuron by neuron to 1.2 % rms (Monte Carlo error ≈ 2 %), correlation 0.9998. The mean is carried entirely by the walls (codimension-one faces), every layer's walls contributing. (Found independently as unlock 4 of foundations-unlocks.md, with a cavity reading of the downstream gain.)

**(v) Dilations and Euler at every cut (Exact).** Let $F_{l\to L}$ be the map $a_l\mapsto a_L$ and $J_{l\to L}=\partial a_L/\partial a_l$. Euler at the cut $l$ gives $a_L=a_lJ_{l\to L}$ pointwise, hence
$$\mathbb E[a_L]=\mu_l\,\mathbb E[J_{l\to L}]+\mathbb E\big[(a_l-\mu_l)J_{l\to L}\big]\qquad\text{for every }l,$$
a mean-transported part and a fluctuation part. **Measured** (n = 64): the mean-transported part is 21 % (rms, relative to $\|\mu_L\|$) of the answer from $a_1$, 46 % from $a_2$, 68 % from $a_4$, 89 % from $a_7$ and 99 % from $a_{15}$. On laws, every arrow commutes with the dilations $D_c$ ($a\mapsto ca$), so the dilation generator is an eigenvector with eigenvalue exactly 1 of every linearised arrow: **a relative scale error $\varepsilon$ in the state at any layer produces exactly the relative error $\varepsilon$ in every final mean, to first order, and is never damped.** With $\|\mu_L\|_{\rm rms}\approx0.95$ the scale component of the state must be right to $\approx10^{-4}$ at every layer.

### 2.2 Rotation invariance of the input: the first layer is exactly known

$\gamma_n$ is $O(n)$-invariant, so the joint law of $(z_1,\dots,z_L)$ depends on $W_1$ only through $C_1=W_1^{\top}W_1$, and $z_1\sim\mathcal N(0,C_1)$ exactly (**Exact**). $C_1$ is Wishart: diagonal $2(1+O(n^{-1/2}))$, correlations $\rho_{jk}\approx\mathcal N(0,1/n)$, spectrum Marchenko–Pastur with ratio one on $[0,8]$, smallest eigenvalues $O(n^{-2})$ (hard edge) *(from memory)*: the first pre-activation is a full-rank but very anisotropic Gaussian. Everything about layer 1 that involves at most three neurons is elementary (**Exact**): $P(z_{1,j}>0)=\tfrac12$; $\mathbb Ea_{1,j}=\sqrt{C_{jj}/2\pi}$; $\mathrm{Var}\,a_{1,j}=C_{jj}(\tfrac12-\tfrac1{2\pi})$; $\mathbb E[a_{1,j}a_{1,k}]=\frac{\sqrt{C_{jj}C_{kk}}}{2\pi}\big(\sqrt{1-\rho^2}+(\pi-\arccos\rho)\rho\big)$ (arc-cosine kernel); $P(z_{1,j}>0,z_{1,k}>0)=\tfrac14+\tfrac{\arcsin\rho}{2\pi}$; triple orthant probabilities $\tfrac18+\tfrac1{4\pi}\sum\arcsin\rho$. Four-neuron orthant probabilities are not elementary. Parity: $x\mapsto-x$ sends $z_1\mapsto-z_1$; $a_1=\tfrac12z_1+\tfrac12|z_1|$ splits into an odd (Gaussian, linear in $x$) and an even part, and every mixed moment of odd total degree in $z_1$ vanishes. No rotation invariance survives the first ReLU: the gate commutes only with permutations and positive diagonal scalings.

### 2.3 Independence of the layers: quenched and annealed

**Conditional structure (Exact).** $W_{l+1}$ is independent of $(W_1,\dots,W_l)$, hence of $P_l$. Conditionally on $P_l$, the $n$ marginal laws of $z_{l+1}$ are the push-forwards of $P_l$ under $n$ independent Gaussian linear functionals $a\mapsto a\cdot w_i$, $w_i\sim\mathcal N(0,\tfrac2nI)$: **the per-neuron laws of a layer are $n$ independent random one-dimensional projections of the previous layer's law.** All their cumulants are explicit multilinear forms:
$$\mathbb E z_{l+1,i}=\mu_l\cdot w_i,\qquad \kappa_k(z_{l+1,i})=K_k(a_l)[w_i,\dots,w_i]\quad(k\ge2),$$
and joint cumulants of several next-layer neurons are the mixed contractions.

**Quenched versus annealed.** The target $\mathbb E_x[a_L\mid W]$ is quenched: an average over the input only, for the given weights. Annealing over $W_{l+1}$ alone makes all neurons of layer $l+1$ identical in law; annealing over all weights gives one number per $(n,L)$ (the infinite-width "mean-field" values). What concentrates as $n\to\infty$ (self-averages): normalised traces $\tfrac1n\mathrm{tr}\,\Sigma_l$, $\tfrac1n\|\mu_l\|^2$, the spectral distribution of $\Sigma_l$, the empirical distribution over neurons of standardized biases and gate probabilities, one scalar per layer such as the thin-shell excess of Section 3.3. What does not: everything indexed by a neuron. $m_i=\mu_l\cdot w_i$ is an $O(1)$ random variable; $s_i^2=w_i^\top\Sigma_lw_i$ fluctuates around $\tfrac2n\mathrm{tr}\,\Sigma_l$ by the relative amount $\sqrt{2\|\Sigma_l\|_F^2}/\mathrm{tr}\,\Sigma_l=\sqrt{2/r_{\rm eff}}$, $r_{\rm eff}=(\mathrm{tr}\Sigma)^2/\|\Sigma\|_F^2$. **Measured**: replacing the quenched $s_i^2$ by the annealed $\tfrac2n\mathrm{tr}\Sigma$ in an otherwise exact Gaussian readout costs $6.5\times10^{-3}$ rms per neuron at n = 1024, against $4.2\times10^{-4}$ with the quenched $s_i^2$; the relative spread of $s_i^2$ around its annealed value is 7 % at $a_2$ rising to 21 % at $a_{15}$ there. (The same quenched/annealed split, with the annealed twirl $\mathbb E_W[W^\top SW]=\tfrac2n\mathrm{tr}(S)I$, is unlock 5 of foundations-unlocks.md.) Ensemble knowledge can only supply self-averaging quantities; $(m_i,s_i^2)$ and the per-neuron higher cumulants must be computed per network.

### 2.4 Symmetries

*Quenched (exact for the given network).*
- Hidden-unit relabelling: $(W_l,W_{l+1})\mapsto(W_lP,P^\top W_{l+1})$ with $P$ a permutation; $a_l\mapsto a_lP$. Any estimator should be equivariant.
- Positive rescaling between layers: $(W_l,W_{l+1})\mapsto(W_lD,D^{-1}W_{l+1})$, $D$ positive diagonal: $a_l\mapsto a_lD$, later layers unchanged; and separate positive homogeneity $W_l\mapsto cW_l\Rightarrow a_L\mapsto ca_L$ ($c>0$).
- Input rotations: $W_1\mapsto QW_1$, $Q\in O(n)$, leaves the law of every layer unchanged (Section 2.2): row permutations and row sign flips of $W_1$ are symmetries. Column sign flips of any $W_l$ are not ($\mathrm{relu}(-z)\ne\mathrm{relu}(z)$).
- The ReLU identity $\mathrm{relu}(z)=\tfrac12z+\tfrac12|z|$ gives $\mathbb E[a_l]=\tfrac12\mu_{l-1}W_l+\tfrac12\mathbb E|z_l|$: half of every mean propagates exactly linearly; the whole difficulty is the even functional $\mathbb E|z_l|$.

*Annealed (the weight ensemble).* Each $W_l$'s law is bi-orthogonally invariant and invariant under row and column sign flips, independently across $l$. Ensemble averages of statistics odd in any single $W_l$ vanish; e.g. the third cumulant of a pre-activation has zero ensemble mean (its sign is largely that of the neuron's standardized bias, Section 3.3), while the fourth has a non-zero, coherent ensemble mean.

### 2.5 Always-on and always-off neurons

**Exact.** Neuron $(l,j)$ is always-on iff $z_{l,j}\ge0$ on the whole sphere (then $a_{l,j}=z_{l,j}$, a linear function of $x$), always-off iff $z_{l,j}\le0$ (then $a_{l,j}\equiv0$). At layer 1 neither can happen: $P(z_{1,j}>0)=\tfrac12$ exactly. For $l\ge2$, $z_{l,j}=a_{l-1}\cdot w_j$ with $a_{l-1}$ in the positive orthant; exact always-on/off has negligible probability for Gaussian weights, but *effective* determinism is common at depth. With $t_j=m_j/s_j$ the standardized bias, the infinite-width theory gives $t\sim\mathcal N(0,\tau_l^2)$, $\tau_l^2=\rho_l/(1-\rho_l)$, where $\rho_l$ is the correlation of $z_l$ between two independent inputs, $\rho_1=0$, $\rho_{l+1}=f(\rho_l)$, $f(\rho)=\big(\sqrt{1-\rho^2}+(\pi-\arccos\rho)\rho\big)/\pi$ (**Derivation**; the arc-cosine recursion is standard *(from memory)*). Fractions of neurons with $|t|>3$ (gate-flip probability $<1.3\times10^{-3}$) and $|t|>4.75$ ($<10^{-6}$):

| layer $l$ | 2 | 4 | 8 | 12 | 16 |
|---|---|---|---|---|---|
| $\tau_l$ (infinite width) | 0.68 | 1.24 | 2.06 | 2.78 | 3.46 |
| $\vert t\vert >3$ | 0 | 1.5 % | 14.5 % | 28 % | 39 % |
| $\vert t\vert >4.75$ | 0 | 0 | 2.1 % | 8.8 % | 17 % |
| measured $\mathrm{std}(t)$, n = 1024 (one network) | 0.69 | 1.24 | 2.10 | 3.02 | 3.70 |
| measured $\vert t\vert >3$ / $\vert t\vert >4.75$, n = 1024 | 0 / 0 | 0.9 % / 0 | 15.5 % / 2.3 % | 33 % / 13 % | 43 % / 22 % |
| measured $\mathrm{std}(t)$, n = 128 | 0.64 | 1.08 | 1.54 | 1.79 | 2.05 |

At n = 1024 the measured spread of $t$ is the infinite-width one to within 7 %, and 3.5 %, 15 % and 26 % of the neurons of $z_8$, $z_{12}$, $z_{16}$ did not change gate once in $1.2\times10^6$ samples. For such neurons $\mathbb Ea=m$ (on) or $0$ (off) up to $O(s\varphi(t)/t^2)$, i.e. below $10^{-7}$ for $|t|>4.75$: the network restricted to them is an exactly known linear (or zero) map. The fraction grows with depth because the mean direction absorbs the energy: the share of $\mathbb E\|a_l\|^2$ in $\|\mu_l\|^2$ is $f(\rho_l)$ = 0.32 at $a_1$, 0.83 at $a_8$, 0.93 at $a_{16}$ (infinite width; the public mini split has 0.9095/(0.9095+0.0748) = 0.92 at $a_{16}$).

### 2.6 Arrows acting on laws and on observables

- **Schrödinger picture.** $p_{l+1}=(A_{l+1})_\#p_l$: linear in the law, affine on the simplex of laws, mapping Dirac masses to Dirac masses. No finite-dimensional family is invariant (Gaussians are not); the exact law is a mixture over faces of linear images of the input Gaussian restricted to cones, $p_l=\sum_\sigma (x\mapsto xM^{(l)}_\sigma)_\#(\gamma_n|_{C_\sigma})$.
- **Heisenberg picture.** $g\mapsto g\circ A_{l+1}$: linear and multiplicative on functions (a composition, or Koopman, operator). The readout pulled back to layer $l$ is $F_{l\to L,i}=\mathrm{relu}_i\circ A_L\circ\cdots\circ A_{l+1}$, continuous and piecewise linear with exponentially many pieces.
- **Duality at every cut (Exact).** $\mu_{L,i}=\langle p_l,F_{l\to L,i}\rangle$ for every $l$: a design may push laws forward, pull observables back, or meet at a cut.
- **Barycentres (Exact).** For a law $p$ and an arrow $A$ that is affine on each face $F$: $\mathrm{bary}(A_\#p)=\sum_F p(F)\,A_F\big(\mathrm{bary}\,p(\cdot\mid F)\big)$. The Jensen gap $\mathbb E\,\mathrm{relu}(z)-\mathrm{relu}(\mathbb Ez)=\tfrac12(\mathbb E|z|-|\mathbb Ez|)\ge0$ is the entire nonlinearity of one arrow on one neuron.
- **Sensitivity (Exact).** The first-order change of $\mu_{L,i}$ under a perturbation $\delta p_l$ is $\langle\delta p_l,F_{l\to L,i}\rangle$. In cumulant coordinates, for any base law, $\partial\mu_{L,i}/\partial\kappa_\alpha=\mathbb E[\partial^\alpha F_{l\to L,i}]/\alpha!$ (Section 3.1). Every derivative of the ReLU beyond the first is a wall distribution, so **every sensitivity in this problem is a face statistic**: gate probabilities, densities of pre-activations at their walls and their derivatives, joint densities at corners where walls meet, multiplied by gated downstream Jacobians.

---

## 3. The information requirement (part b)

The question is which quantities about each layer's law an estimator must know, and how accurately, for a raw final-layer MSE of $10^{-8}$ or below. It is answered backwards from the readout: first the final neuron's own law (3.1–3.2), then the previous layer's joint law through the linear half of the arrow (3.3), then through the gate half (3.4), then through depth (3.5). Section 3.6 collects the answer.

### 3.1 The readout and its exact sensitivities

$\mu_{L,i}=\mathbb E\,\mathrm{relu}(z_{L,i})$ is a linear functional of the marginal law of $z_{L,i}$.

**Exact sensitivities.** For any law $p$ of a scalar $z$ with finite moments, the direction that changes the $k$-th cumulant alone is $\delta p=(-\partial_z)^kp/k!$ (its characteristic function is $(iu)^k\hat p(u)/k!$), and integration by parts gives, for any $g$,
$$\frac{\partial\,\mathbb E[g(z)]}{\partial\kappa_k}=\frac{\mathbb E[g^{(k)}(z)]}{k!}.$$
The same holds for joint cumulants with a multi-index ($\partial\mathbb E g/\partial\kappa_\alpha=\mathbb E[\partial^\alpha g]/\alpha!$). For $g=\mathrm{relu}$, whose derivatives are the gate and then the wall distribution and its derivatives:
$$\frac{\partial\mu}{\partial\kappa_1}=P(z>0),\qquad \frac{\partial\mu}{\partial\kappa_2}=\tfrac12p_z(0),\qquad\frac{\partial\mu}{\partial\kappa_3}=-\tfrac16p_z'(0),\qquad\frac{\partial\mu}{\partial\kappa_4}=\tfrac1{24}p_z''(0).$$
**The readout is sensitive to the mean through the gate probability, to the variance through the density of $z$ at the wall, and to the third and fourth cumulants through the slope and curvature of that density at the wall** (**Exact**). An equivalent exact form of the readout itself is $\mathbb E\,\mathrm{relu}(z)=\mathbb E[z]\,P(z>0)+p_z(0)\,\tau_z(0)$ with the Stein kernel $\tau_z$ (unlock 3(c) of foundations-unlocks.md); the two agree, since expanding $\tau_z(0)$ in cumulants reproduces the series below. With $z\sim\mathcal N(m,s^2)$, $t=m/s$: $\Phi(t)$, $\varphi(t)/2s$, $-t\varphi(t)/6s^2$, $(t^2-1)\varphi(t)/24s^3$.

**The full series (Exact).** With $\hat z=(z-m)/s$ and Hermite moments $h_k=\mathbb E[He_k(\hat z)]$,
$$\mathbb E\,\mathrm{relu}(z)=s\big[t\Phi(t)+\varphi(t)\big]+s\varphi(t)\sum_{k\ge3}\frac{(-1)^kHe_{k-2}(t)}{k!}\,h_k,$$
because $\mathbb E[\mathrm{relu}(m+sG)He_k(G)]=s(-1)^kHe_{k-2}(t)\varphi(t)$ for $k\ge2$. In standardized cumulants $\lambda_k=\kappa_k/s^k$: $h_3=\lambda_3$, $h_4=\lambda_4$, $h_5=\lambda_5$, $h_6=\lambda_6+10\lambda_3^2$, $h_7\ni35\lambda_3\lambda_4$, $h_8\ni35\lambda_4^2$. Since $\lambda_3,\lambda_4=O(1/n)$ and $\lambda_5,\lambda_6=O(1/n^2)$ (Section 5), the *first-order* readout keeps $\lambda_3,\lambda_4$:
$$\mu\approx s\big[t\Phi+\varphi\big]-\tfrac{s\,t\varphi}{6}\lambda_3+\tfrac{s(t^2-1)\varphi}{24}\lambda_4 .$$
Analogous series give $P(z>0)=\Phi(t)+\varphi(t)\big[\tfrac{He_2(t)}6\lambda_3-\tfrac{He_3(t)}{24}\lambda_4+\dots\big]$ and $p_z(0)=\tfrac{\varphi(t)}{s}\big[1-\tfrac{He_3(t)}6\lambda_3+\tfrac{He_4(t)}{24}\lambda_4+\dots\big]$.

**Check on Monte Carlo (n = 64, depth 16, 2×10⁶ samples, 32 neurons at each of four layers, two networks).** Each sensitivity was measured as the response of $\mathbb E\,\mathrm{relu}(z)$ to a perturbation of the law with known cumulants, averaged analytically sample by sample (no extra noise): a shift (gate probability), an independent Gaussian increment (wall density), an antithetic exponential increment (slope), a Laplace-minus-Gaussian increment of equal variance (curvature). Relative rms error of the Gaussian (G) and first-order (E1) predictions:

| layer | median $\vert \lambda_3\vert $, $\lambda_4$ | $P(z>0)$: G / E1 | $p(0)$: G / E1 | $p'(0)$: G / E1 | $p''(0)$: G / E1 |
|---|---|---|---|---|---|
| $z_2$ | 0.10–0.12, 0.13–0.14 | 1.0–1.1 % / 0.1 % | 3.1–3.3 % / 0.5–0.6 % | 10–11 % / 6 % | 10–14 % / 11 % |
| $z_4$ | 0.14–0.24, 0.28–0.33 | 1.5 % / 0.2 % | 5.1–5.4 % / 1.2–1.6 % | 13–16 % / 5–6 % | 20–24 % / 10–14 % |
| $z_8$ | 0.45–0.59, 0.57–0.66 | 2.3 % / 0.5–0.6 % | 10 % / 2.5–3.2 % | 20–23 % / 9–12 % | 31–43 % / 14–24 % |
| $z_{16}$ | 0.66–0.84, 1.1 | 1.8–3.3 % / 0.7–1.4 % | 19–21 % / 6–8 % | 35–36 % / 25–28 % | 53–73 % / 26–46 % |

The shift response reproduces the sample gate probability to $3\times10^{-5}$ (relative), as it must. The Gaussian sensitivities are wrong by an amount proportional to the non-Gaussianity of the neuron's own law, and the first-order correction removes most of it while $\lambda_{3,4}\lesssim0.3$. At $n=1024$, where $\lambda_{3,4}$ are 7–20 times smaller (Section 5), the Gaussian values of the gate probability and the wall density are accurate to ≈ $1\times10^{-3}$ absolute at the last layer (Monte Carlo noise of the check ≈ $4\times10^{-4}$) and better than ≈ 0.3 % of $1/s$ (the noise level of the check).

**Readout truncation (widths 32–1024, depth 16; networks per width in brackets).** Root-mean-square over neurons of the error of each readout *given the exact per-neuron cumulants* (same samples, so Monte Carlo noise cancels); the last two columns are the rms sizes of the $\lambda_3$ and $\lambda_4$ terms themselves:

| width (networks) | layer | Gaussian, exact $(m,s)$ | $\times n$ | + $\lambda_3,\lambda_4$ (first order) | $\times n^2$ | second order | $\lambda_3$ term | $\lambda_4$ term |
|---|---|---|---|---|---|---|---|---|
| 32 (4) | $z_{2}$ | 1.3e-02 | 0.40 | 1.1e-03 | 1.1 | 6.8e-04 | 1.0e-02 | 3.8e-03 |
| 32 (4) | $z_{8}$ | 1.4e-02 | 0.46 | 2.9e-03 | 3.0 | 8.4e-03 | 1.4e-02 | 8.8e-03 |
| 32 (4) | $z_{16}$ | 6.6e-03 | 0.21 | 2.1e-03 | 2.2 | 2.8e-03 | 9.5e-03 | 6.2e-03 |
| 64 (2) | $z_{2}$ | 6.6e-03 | 0.42 | 3.5e-04 | 1.4 | 8.3e-05 | 5.7e-03 | 2.1e-03 |
| 64 (2) | $z_{8}$ | 6.6e-03 | 0.42 | 8.2e-04 | 3.4 | 4.1e-04 | 7.4e-03 | 3.0e-03 |
| 64 (2) | $z_{16}$ | 5.3e-03 | 0.34 | 1.1e-03 | 4.5 | 1.0e-03 | 5.8e-03 | 3.8e-03 |
| 128 (1) | $z_{2}$ | 2.8e-03 | 0.36 | 6.6e-05 | 1.1 | 3.0e-05 | 2.4e-03 | 1.0e-03 |
| 128 (1) | $z_{8}$ | 4.2e-03 | 0.54 | 3.9e-04 | 6.4 | 8.6e-05 | 4.3e-03 | 1.8e-03 |
| 128 (1) | $z_{16}$ | 2.3e-03 | 0.29 | 4.2e-04 | 6.8 | 2.1e-04 | 2.4e-03 | 1.5e-03 |
| 256 (1) | $z_{2}$ | 1.6e-03 | 0.40 | 3.6e-05 | 2.4 | 2.1e-05 | 1.4e-03 | 5.0e-04 |
| 256 (1) | $z_{8}$ | 2.0e-03 | 0.50 | 8.1e-05 | 5.3 | 1.5e-05 | 2.0e-03 | 8.4e-04 |
| 256 (1) | $z_{16}$ | 1.5e-03 | 0.38 | 1.2e-04 | 7.7 | 2.9e-05 | 1.6e-03 | 7.5e-04 |
| 1024 (1) | $z_{2}$ | 4.2e-04 | 0.43 | ≤ 7.3e-05 | (floor) | ≤ 5.3e-05 | 3.6e-04 | 1.5e-04 |
| 1024 (1) | $z_{8}$ | 5.1e-04 | 0.52 | ≤ 2.8e-05 | (floor) | ≤ 2.0e-05 | 5.1e-04 | 1.8e-04 |
| 1024 (1) | $z_{16}$ | 4.2e-04 | 0.43 | ≤ 1.7e-05 | (floor) | ≤ 1.1e-05 | 4.2e-04 | 1.3e-04 |

Widths 32–256: networks with weights `default_rng(1000+seed)`, $10^7$ samples ($8\times10^6$ at 256); width 1024: bench network 7301001, $1.2\times10^6$ samples, where the first- and second-order columns sit at the sampling floor of the check and bound the residual from above. Final-layer MSE of the three readouts: 5.5·10⁻⁵ / 5.8·10⁻⁶ / 1.0·10⁻⁵ (n = 32), 3.0·10⁻⁵ / 1.4·10⁻⁶ / 1.2·10⁻⁶ (64), 5.1·10⁻⁶ / 1.7·10⁻⁷ / 4.6·10⁻⁸ (128), 2.2·10⁻⁶ / 1.4·10⁻⁸ / 8.6·10⁻¹⁰ (256), **1.7·10⁻⁷ / ≤ 2.9·10⁻¹⁰ / ≤ 1.2·10⁻¹⁰ (1024)**. The Gaussian readout's error is $\approx(0.2$–$0.57)/n$ at every width; the first-order readout's is $O(n^{-2})$ (×$n^2$ = 1–8 up to n = 256). At $n\le64$, where $\lambda_{3,4}\gtrsim0.3$ at depth, the second-order truncation is not uniformly better than the first: the expansion is asymptotic.

### 3.2 What the final layer requires

At the final layer the standardized biases are spread, $t\sim\mathcal N(0,\tau_{16}^2)$ with $\tau_{16}\approx3.5$ (infinite width; 3.70 measured at n = 1024), and $s\approx0.39$ (0.42 measured). Root-mean-square sensitivities over this population (**Derivation** from the Gaussian values, which are accurate to the stated level at n = 1024):

| quantity of $z_{L,i}$ | rms sensitivity of $\mu_{L,i}$ | size at n = 1024 | tolerance, rms, for MSE ≤ 2.5·10⁻⁹ from this source alone | relative |
|---|---|---|---|---|
| mean $m_i$ | $\langle\Phi(t)^2\rangle^{1/2}=0.66$ | rms 1.35 | $7.6\times10^{-5}$ | $5.6\times10^{-5}$ |
| variance $s_i^2$ | $\langle(\varphi/2s)^2\rangle^{1/2}=0.23$ | ≈ 0.15 | $2.2\times10^{-4}$ | 0.14 % |
| $\lambda_{3,i}$ | $\langle(st\varphi/6)^2\rangle^{1/2}=0.0081$ | rms 0.10 (measured, n = 1024); its readout term $4.2\times10^{-4}$ rms | $6\times10^{-3}$ | ≈ 12 % of the $\lambda_3$ term |
| $\lambda_{4,i}$ | $\langle(s(t^2-1)\varphi/24)^2\rangle^{1/2}=0.0025$ | 0.056, coherent (measured, n = 1024); its readout term $1.3\times10^{-4}$ rms | $2\times10^{-2}$ | ≈ 40 % of the $\lambda_4$ term |
| $\lambda_5,\lambda_6,\lambda_3^2,\dots$ | — | $O(n^{-2})$ | — | negligible: the first-order readout's own bias is ≤ $1.7\times10^{-5}$ rms |

Three conclusions. (i) The mean must be essentially exact; it is: $m_i=\mu_{L-1}\cdot w_i$, so this is a requirement on $\mu_{L-1}$ (Section 3.3). (ii) The variance must be right to 0.14 % per neuron. (iii) The neuron's own third and fourth cumulants are *necessary* (ignoring them costs $1.7\times10^{-7}$ in MSE) and *sufficient* at first order, with loose relative tolerances (≈ 12 % of the $\lambda_3$ term and ≈ 40 % of the $\lambda_4$ term). Measured at n = 1024 the rms sensitivities are 0.68 (mean) and 0.21 (variance), within 10 % of the values above. Raw $5\times10^{-9}$ tightens every tolerance by $\sqrt2$.

### 3.3 One step back through the linear half: the joint law of $a_{L-1}$

**Exact transfer and its error channels.** By Section 2.3, $m_i=\mu\cdot w_i$, $s_i^2=w_i^\top\Sigma w_i$, $\kappa_3(z_i)=K_3[w_i^{\otimes3}]$, $\kappa_4(z_i)=K_4[w_i^{\otimes4}]$, where $(\mu,\Sigma,K_3,K_4)$ describe $a_{L-1}$ and the $w_i$ are independent of them **and of any error an estimator makes in them** (as long as the estimate of layer $L-1$ is formed before $W_L$ is used). For an error tensor $\delta T$ of order $k$, Isserlis' theorem with $w\sim\mathcal N(0,\sigma^2I)$, $\sigma^2=2/n$, gives the mean and spread of the induced per-neuron errors (**Exact**):

| order | mean over neurons (coherent channel) | variance over neurons |
|---|---|---|
| 1 | 0 | $\sigma^2\Vert \delta\mu\Vert ^2$ |
| 2 | $\sigma^2\,\mathrm{tr}\,\delta\Sigma$ | $2\sigma^4\Vert \delta\Sigma\Vert _F^2$ |
| 3 | 0 | $\sigma^6\big[6\Vert \delta K_3\Vert _F^2+9\Vert v\Vert ^2\big]$, $v_c=\sum_j\delta K_{3,jjc}$ |
| 4 | $3\sigma^4\sum_{jk}\delta K_{4,jjkk}$ | $\sigma^8\big[24\Vert \delta K_4\Vert _F^2+72\Vert \mathrm{tr}_{12}\delta K_4\Vert _F^2\big]$ |

*Coherent channels* (traces) are not averaged down by the $n$ inputs: a bias of $\varepsilon$ in every per-neuron variance of $a_{L-1}$ shifts every $s_i^2$ by $2\varepsilon$. *Incoherent channels* are averaged: an error of rms $\varepsilon$ spread over all $n^2$ entries of $\Sigma$ moves $s_i^2$ by only $2\sqrt2\,\varepsilon$. Partial traces (the vector $v$ for $K_3$, the matrix $\mathrm{tr}_{12}$ for $K_4$) sit in between: an error that is coherent along the repeated index is amplified by up to $n$ relative to the Frobenius channel. **Check with the actual weights** (n = 128, layers $a_4$ and $a_{15}$): an error of +1 % on every variance shifts every next-layer $s_i^2$ by $5.60\times10^{-3}$ (predicted $\sigma^2\mathrm{tr}\,\delta\Sigma=5.65\times10^{-3}$) with a spread across neurons 6.5 times smaller; a random symmetric error with the same Frobenius norm produces no shift and the same spread ($1.0\times10^{-3}$, predicted $0.95\times10^{-3}$). The ratio of the two channels grows as $\sqrt{n/2}$, to ≈ 22 at $n=1024$.

**Tolerances at n = 1024** (from Section 3.2, each source alone at MSE 2.5·10⁻⁹, **Derivation**):
- $\mu_{L-1}$: rms per-neuron error ≤ $5.4\times10^{-5}$ (random direction; directions aligned with $W_L$'s structure do not exist since $W_L$ is independent).
- $\Sigma_{L-1}$, coherent: $\tfrac1n\mathrm{tr}\,\delta\Sigma\le1.1\times10^{-4}$, i.e. the average per-neuron variance of $a_{L-1}$ to ≈ 0.13 % (its value is 0.091 at n = 1024).
- $\Sigma_{L-1}$, incoherent: $\|\delta\Sigma\|_F\le0.078$, i.e. rms entry error $\le7.8\times10^{-5}$, ≈ 0.6 % of $\|\Sigma_{L-1}\|_F$.
- $K_3(a_{L-1})$: $\big[6\|\delta K_3\|_F^2+9\|v\|^2\big]^{1/2}\le4.8$, against 118 for $K_3$ itself (measured at n = 1024 as $\mathrm{rms}_i\,\kappa_3(z_{L,i})/\sigma^3$): ≈ 4 % in this contraction norm for an error spread like a random tensor, ≈ 12 % for an error proportional to the true third cumulant (which concentrates on neurons with large $|t|$, where the sensitivity is small).
- $K_4(a_{L-1})$: the coherent double trace $\sum_{jk}K_{4,jjkk}$ to ±48 against its value 168 (measured at n = 1024), i.e. ≈ 28 %, and the per-neuron part to a similar relative accuracy.

**Every index pattern counts (Derivation, Measured).** Split $K_k$ by the coincidence pattern of its indices — for $k=3$: $(3)$, $(2,1)$, $(1,1,1)$; for $k=4$: $(4),(3,1),(2,2),(2,1,1),(1,1,1,1)$. A pattern with $p$ distinct indices has $\Theta(n^p)$ entries. For a law whose pairwise correlations are $O(n^{-1/2})$ (layer 1 exactly; deeper layers in the bulk), a connected joint cumulant among $p$ distinct neurons is a tree of $p-1$ correlation edges, of size $O(n^{-(p-1)/2})$, so **every pattern carries Frobenius mass $\Theta(n)$** and contributes to $\kappa_k(z_i)$ at the same order $O(n^{(1-k)/2})$. Measured at widths 24–64 (Appendix A.4; noise-corrected by independent replicas):

$n\times$ rms over neurons of each pattern's contribution to $\lambda_3(z_{l+1,i})$ (contributions divided by the layer's mean $s^3$):

| layer | n = 24: (3) / (2,1) / (1,1,1) / total | n = 32: (3) / (2,1) / (1,1,1) / total | n = 64: (3) / (2,1) / (1,1,1) / total |
|---|---|---|---|
| $a_{1}$ | 7.6 / 7.9 / 2.1 / 11.1 | 8.4 / 5.9 / 2.0 / 11.3 | 6.6 / 5.7 / 2.2 / 11.4 | 
| $a_{2}$ | 11.2 / 14.0 / 3.1 / 23.6 | 6.7 / 15.9 / 6.0 / 21.5 | 7.1 / 11.7 / 7.0 / 20.3 | 
| $a_{4}$ | 10.2 / 32.2 / 18.4 / 25.2 | 16.0 / 31.9 / 22.7 / 42.4 | 6.8 / 29.4 / 13.8 / 40.0 | 
| $a_{8}$ | 18.8 / 48.5 / 23.6 / 55.1 | 13.9 / 59.1 / 20.3 / 83.6 | 6.5 / 54.5 / 22.6 / 67.7 | 
| $a_{15}$ | 11.1 / 44.7 / 23.6 / 45.9 | 13.5 / 88.3 / 35.1 / 87.8 | 10.7 / 74.9 / 35.6 / 86.5 |

$n\times$ rms of each pattern's contribution to $\lambda_4(z_{l+1,i})$ (widths 24 and 32, where the full $K_4$ fits):

| layer | n = 24: (4) / (3,1) / (2,2) / (2,1,1) / (1,1,1,1) / total | n = 32: (4) / (3,1) / (2,2) / (2,1,1) / (1,1,1,1) / total |
|---|---|---|
| $a_{1}$ | 12 / 7 / 12 / 5 / 1 / 20 | 12 / 7 / 6 / 2 / 1 / 14 | 
| $a_{4}$ | 15 / 30 / 67 / 45 / 22 / 45 | 25 / 54 / 53 / 48 / 32 / 104 | 
| $a_{15}$ | 18 / 63 / 78 / 66 / 33 / 99 | 18 / 61 / 155 / 229 / 108 / 236 |

At $a_1$ (the ReLU of an exact Gaussian) the diagonal pattern is the largest single part of $\lambda_3$, the off-diagonal ones together are as large, and $n\times$rms is width-independent (total 11.1, 11.3, 11.4 at n = 24, 32, 64), as the counting predicts. From $a_4$ on, the $(2,1)$ pattern is the largest part of $\lambda_3$ and the diagonal the smallest (≤ 15 % of the total at $a_8$–$a_{15}$, n = 64); for $\lambda_4$ the $(2,2)$, $(2,1,1)$ and $(1,1,1,1)$ patterns exceed the diagonal by factors 3–10 at depth. No pattern can be dropped, and the all-distinct parts ($n^3$ and $n^4$ entries) can be neither dropped nor stored (Section 4).

**The coherent fourth cumulant is one scalar per layer (Exact).** For centred $X=a_l-\mu_l$,
$$\mathbb E_w\big[\kappa_4(z_{l+1,i})\big]=3\sigma^4\sum_{jk}K_{4,jjkk}=3\sigma^4\Big[\mathrm{Var}\,\|X\|^2-2\|\Sigma_l\|_F^2\Big]:$$
the neuron-averaged fourth cumulant of the next pre-activation is the *thin-shell excess* of the current layer — how much more the squared norm of the centred activation fluctuates than it would for a Gaussian with the same covariance. **Measured**: at n = 1024 the identity gives the neuron-averaged $\kappa_4$ to 0.1–2.2 % at all five layers checked ($z_3,z_5,z_9,z_{13},z_{16}$); at n = 128 it holds to within 1–9 % at 13 of 15 layers and 15–21 % at the other two, the residual being the finite sample of $n$ columns.

**The per-neuron non-Gaussianity is largely a common scale mixture (Derivation, Measured).** Positive homogeneity writes $z_{l+1,i}=\|a_l\|\,(\hat a_l\cdot w_i)$. If the direction factor is close to Gaussian (a random projection of a high-dimensional direction) and nearly independent of the norm, every pre-activation of the layer is a Gaussian scale mixture with one common mixing variable $\rho=\|a_l\|/\mathbb E\|a_l\|$; with $\varepsilon^2=\mathrm{Var}\,\rho$ the law of total cumulance gives, to leading order, $\lambda_{3,i}\approx6\varepsilon^2t_i$ and $\lambda_{4,i}\approx12\varepsilon^2$, independent of $t_i$. Measured at n = 128: $\lambda_3$ is a linear function of $t$ with $R^2$ = 0.46 ($z_2$), 0.76 ($z_4$), 0.86 ($z_8$), 0.79 ($z_{16}$), with slope 0.9–0.5 times $6\varepsilon^2$ computed from the norm statistics (shallow to deep); $\lambda_4$ has no $t^2$ dependence (coefficient ≤ 0.03 against an intercept of 0.07–0.76) and its mean is 0.46–0.90 times $12\varepsilon^2$. At n = 1024 (bench network): $R^2$ = 0.52, 0.78, 0.87, 0.89 at $z_2,z_4,z_8,z_{16}$, slope 1.00, 0.78, 0.75, 0.67 times $6\varepsilon^2$, $\lambda_4$ without $t^2$ dependence (coefficient ≤ $3\times10^{-4}$) and 0.47–0.72 times $12\varepsilon^2$. In the readout the common-scale part carries 30–46 % of the $\lambda_3$ term and 58–80 % of the $\lambda_4$ term at n = 1024 (34–54 % and 55–72 % at n = 128); the rest is specific to the neuron. This is the per-neuron face of the mean-direction spike of $\Sigma_l$ (Section 5), and the reason why the third cumulant, whose ensemble mean is zero, is nevertheless far from random given the neuron's own bias.

### 3.4 One step back through the gate half: the pairwise law of $z_{L-1}$

$\Sigma_{L-1}=\mathrm{Cov}(\mathrm{relu}(z_{L-1}))$ is a pairwise functional of the law of $z_{L-1}$. By the multivariate form of 3.1 its exact sensitivities are joint face statistics (**Exact**): for $j\ne k$,
$$\frac{\partial\,\mathbb E[a_ja_k]}{\partial C_{jk}}=P(z_j>0,z_k>0),\quad\frac{\partial\,\mathbb E[a_ja_k]}{\partial\kappa(z_j,z_j,z_k)}=\tfrac12\,\mathbb E[\delta(z_j)1_{z_k>0}],\quad\frac{\partial\,\mathbb E[a_ja_k]}{\partial\kappa(z_j,z_j,z_k,z_k)}=\tfrac14\,p_{jk}(0,0),\quad\frac{\partial\,\mathbb E[a_ja_k]}{\partial\kappa(z_j,z_j,z_j,z_k)}=\tfrac16\,\mathbb E[\delta'(z_j)1_{z_k>0}]$$
(the first is Price's theorem in the Gaussian case *(from memory)*): the joint gate probability, the density on one wall times the other gate, the density at the corner where two walls meet. So the second-order structure of $a$ needs, from $z$: the pairwise means and covariances (Gaussian part), the $(2,1)$ slices $\kappa(z_j,z_j,z_k)$, and the $(3,1)$ and $(2,2)$ slices of the fourth cumulant.

**Information test (Measured).** How much of the pairwise law of $z_l$ is needed to know $\Sigma_l$? Predict $\Sigma_l$ from the *exact* pairwise statistics of $z_l$ truncated at cumulant order 2 (G: means and covariances, bivariate Gaussian integrals), order 3 (E3: plus the marginal and $(2,1)$ third cumulants, entering through the face-statistic sensitivities above) or order 4 (E4: plus the $(4)$, $(3,1)$, $(2,2)$ fourth cumulants); feed the error of each prediction through the next layer's quadratic forms and readout sensitivity. This is a one-step oracle measurement of information content, not an estimator (widths 16–128, noise-corrected by replicas). Root-mean-square per-neuron error induced at layer $l+1$ by the error of $\Sigma_l$ alone, averaged over layers $l=2,\dots,15$ (layer 1 is exact for every closure, since $z_1$ is Gaussian):

| width | Gaussian (G) | + marginal $\kappa_3$ only | + $(2,1)$ slices (E3) | + marginal $\kappa_3,\kappa_4$ only | + $(3,1)$, $(2,2)$ slices (E4) | off-diagonal relative error of $\Sigma_l$: G / E3 / E4 |
|---|---|---|---|---|---|---|
| 16 | 1.6e-02 | 2.3e-02 | 1.9e-02 | 2.4e-02 | 4.8e-03 | 30 % / 13.6 % / 3.4 % |
| 24 | 5.8e-03 | 6.7e-03 | 2.6e-03 | 5.6e-03 | 7.3e-04 | 18 % / 5.3 % / 2.0 % |
| 32 | 4.1e-03 | 4.8e-03 | 2.4e-03 | 3.8e-03 | 4.8e-04 | 15 % / 5.6 % / 1.4 % |
| 64 | 3.6e-03 | 4.0e-03 | 1.1e-03 | 3.6e-03 | 3.5e-04 | 15 % / 3.5 % / 1.2 % |
| 128 | 2.0e-03 | 2.2e-03 | 5.3e-04 | 2.2e-03 | 1.4e-04 | 16 % / 2.8 % / 0.8 % |
| fitted exponent, widths 32–128 | -0.51 | -0.55 | -1.10 | -0.41 | -0.91 |  |

Readings. (i) The per-neuron (marginal) cumulants do not improve the covariance at all: the second column is no better than the first, nor the fourth than the third. What matters are the *pairwise* slices. (ii) The order-2 (Gaussian) prediction's relative error on the off-diagonal covariance does not fall with width (15–16 % from n = 32 to 128); its induced error falls only as $\approx n^{-1/2}$ over the measured range. (iii) The $(2,1)$ slices cut it 4–8-fold, and the fourth-order slices another 4-fold, with errors falling roughly as $1/n$. (iv) Against a per-layer budget of ≈ $4\times10^{-5}$ (each layer's mean errors reach the output with summed gain ≈ 5, Section 3.5): at $n=1024$ the order-2 prediction induces ≈ $(3$–$8)\times10^{-4}$ (6–18× over), the third-order slices ≈ $(5$–$7)\times10^{-5}$ (at the budget), the fourth-order slices ≈ $2\times10^{-5}$ (within it). (v) Consistency: a Gaussian closure injects at every layer this error plus the Gaussian readout's own (≈ $4\times10^{-4}$, Section 3.1); accumulated with the gains of Section 3.5 (factor ≈ 2.4) that is ≈ $(1$–$2)\times10^{-3}$ rms at the output, and the bench measures $2.1\times10^{-3}$ rms (raw $4.3\times10^{-6}$) for a Gaussian covariance closure at $n=1024$. The orders agree.

### 3.5 Through depth

Errors made at layer $l$ reach the output through the pulled-back observable $F_{l\to L}$ (Section 2.6). Three channels behave differently.

- **Scale (dilation) channel: gain exactly 1** (Section 2.1 v).
- **Mean channel.** A mean error $\delta\mu_l$ changes $\mu_L$ by $\delta\mu_l\,\mathbb E[J_{l\to L}]$, the input-averaged gated propagator. If the gates were independent across layers, $\mathbb E[J]$ would be a product of $W_k\,\mathrm{diag}(\pi_k)$, and a random error would lose power by $2\langle\pi_k^2\rangle$ per layer. In the infinite-width theory $2\langle\pi^2\rangle=2P(z(x)>0,z(x')>0)=\tfrac12+\tfrac{\arcsin\rho_k}{\pi}=f'(\rho_k)$ for independent inputs $x,x'$ (**Derivation**): the damping of mean errors is the derivative of the correlation map, 0.5 at layer 1 rising to 0.87 at layer 16; the cumulative power gain from $a_l$ to the output is 0.013 from the input, 0.027 from $a_1$, 0.21 from $a_7$, 0.87 from $a_{15}$. **Measured** (exact reverse-mode averages): $\|\mathbb E J_{l\to L}\|_F^2/n$ = 0.034 / 0.0092 from the input, 0.074 / 0.018 from $a_1$, 0.21 / 0.050 from $a_4$, 0.37 / 0.12 from $a_7$, 0.50 / 0.24 from $a_{10}$, 0.89 / 0.84 from $a_{15}$ (n = 64 / 128), against the gate-independent products 0.014 / 0.0046, 0.029 / 0.0092, 0.13 / 0.040, 0.34 / 0.10, 0.46 / 0.24, 0.89 / 0.83 built from the measured gate probabilities. Gates of different layers are correlated, which raises the gain at shallow layers by up to ×2.4 (n = 64) and ×2 (n = 128); at depth the product is exact. Summed over source layers 1–15 the gains are 6.2 (n = 64), 3.7 (n = 128) and 5.0 (infinite width), so independent spread mean errors of rms $\varepsilon$ made at every layer reach the output as ≈ $(2.2$–$2.7)\,\varepsilon$. The averaged propagator is anisotropic: its top singular value is 1.2–3.2 (n = 64) and 0.9–2.8 (n = 128), of the order of the spectral norm of one weight matrix ($2\sqrt2$), and its top output direction is the final mean direction (overlap 0.6–0.99). Errors independent of the downstream weights see the Frobenius gain above; errors aligned with the top input directions would be amplified rather than damped.
- **Second-order (trace) channel.** Adding independent isotropic noise of variance $\varepsilon^2$ to $a_l$ (a pure covariance error $\delta\Sigma_l=\varepsilon^2I$) changes $\mu_{L,i}$ by $\tfrac{\varepsilon^2}2\mathbb E[\Delta_{a_l}F_{l\to L,i}]$, a sum over all downstream walls. At $l=0$ this is $\tfrac{\varepsilon^2}2\mu_{L,i}$ exactly (wall formula); at the last cut it is $\tfrac{\varepsilon^2}2p_{z_{L,i}}(0)\|w_i\|^2$, the variance sensitivity of 3.1 (checked at n = 64: rms 1.403 and mean 0.859 from the measured wall densities, against 1.401 and 0.859 from the noise injection; at n = 1024 the same formula gives rms 0.84, mean 0.49). **Measured** (n = 64): $\mathbb E[\Delta_{a_l}F]$ has rms 1.1–1.7 and mean +0.47 to +0.93 per output neuron at every source layer $l=1,\dots,15$, with no decay in depth. A coherent relative error $\varepsilon_v$ in the per-neuron variances of $a_l$ therefore moves the final means by $\approx\tfrac12\varepsilon_v\bar v_l\,\mathbb E[\Delta F]$, where $\bar v_l$ is the mean per-neuron variance of $a_l$ (0.68 at $l=1$ falling to ≈ 0.08 at $l=15$): measured at n = 64, a coherent relative error in the per-neuron variances of a single layer, using the whole budget alone, must stay below $2.8\times10^{-4}$ at $a_1$, $(4$–$5)\times10^{-4}$ at $a_2$–$a_4$, $(6$–$8)\times10^{-4}$ at $a_5$–$a_{13}$ and $1.2\times10^{-3}$ at $a_{15}$. Because the mean of $\mathbb E[\Delta F]$ is positive at every layer, errors of one sign at several layers add: a relative bias $\varepsilon$ at all of layers 2–15 moves the final means by $1.24\,\varepsilon$ on average, so such a bias must be below $8\times10^{-5}$. Spread (zero-mean, $W$-independent) variance errors are averaged by the next layer's quadratic forms and are $\sqrt n$ times less harmful.

### 3.6 What must be known about each layer's law, and how accurately

For raw MSE $10^{-8}$ at $n=1024$, each line using the whole budget alone (halve the tolerances, roughly, when all lines share it; divide by $\sqrt2$ for $5\times10^{-9}$). "Coherent" means the same sign across neurons; "spread" means independent of the next layer's weights.

| about the law of layer $l$ | needed because | tolerance | where |
|---|---|---|---|
| the overall scale (dilation component of the state) | propagated with gain exactly 1 | $10^{-4}$ relative, at every layer | 2.1 (v) |
| per-neuron means $\mu_{l,i}$ | the next layer's $m$ is $\mu_l W_{l+1}$ exactly; mean errors reach the output with summed gain ≈ 5 | ≈ $4\times10^{-5}$ rms per neuron at every layer (spread errors) | 3.3, 3.5 |
| per-neuron pre-activation variances $s_{l,i}^2$ | each layer's readout of $\mu_{l,i}$ | 0.14 % per neuron (spread); coherent part below ≈ $3\times10^{-4}$ (shallow) to $10^{-3}$ (deep) for one layer alone, below $8\times10^{-5}$ if of one sign at all layers | 3.2, 3.5 |
| per-neuron $\lambda_{3},\lambda_4$ of $z_{l,i}$ | each layer's readout | the $\lambda_3$ term to ≈ 12 %, the $\lambda_4$ term to ≈ 40 % of their sensitivity-weighted size (n = 1024); a common-scale (norm-process) part carries 30–45 % of the first and 60–80 % of the second | 3.2, 3.3 |
| the covariance $\Sigma_l$ off the diagonal | the next layer's quadratic forms $w_i^\top\Sigma_lw_i$ | ≈ 0.6 % in Frobenius norm (spread); the trace to 0.13 % | 3.3 |
| the pairwise cumulant slices of $z_l$: $(2,1)$; $(3,1)$, $(2,2)$ | $\Sigma_l=\mathrm{Cov}(\mathrm{relu}(z_l))$ through the gate | the Gaussian part alone is not enough; with the third- and fourth-order slices the induced error is within the per-layer budget (per layer at n = 1024: Gaussian part only ≈ 3–8·10⁻⁴, with third-order slices ≈ 5–7·10⁻⁵, with fourth-order slices ≈ 2·10⁻⁵) | 3.4 |
| every index pattern of $K_3(a_l)$, $K_4(a_l)$ | the next layer's $\lambda_3,\lambda_4$ and pairwise slices | contraction norms to ≈ 4–12 % ($K_3$) and ≈ 30 % (coherent part of $K_4$); the all-distinct parts at leading order | 3.3 |

The requirements chain: per-neuron fourth-order information at layer $l+1$ needs joint fourth-order information at layer $l$, and pairwise slices at layer $l$ need the full tensors at layer $l-1$. The hierarchy does not close by truncation in the order of the cumulants, because every index pattern counts; it closes, if at all, through structure (the $1/n$ expansion, the common-scale mixture, the tree-like connected cumulants) — which is exactly what a design has to supply.

---

## 4. Cost floors at n = 1024 (part c)

Prices are the flopscope 0.12.1 measurements of [../streams/costmodel/REPORT.md](../streams/costmodel/REPORT.md) (**Fact**). 1 unit = one dense $(1024\times1024)@(1024\times1024)$ product $=2^{31}$ FLOPs; $B=1024$ units; the 0.1 floor is 102.4 units, **6.4 units per layer** at $L=16$.

| operation, per layer | FLOPs | units | remark |
|---|---|---|---|
| touch $W_l$ once (one vector through the layer) | $2n^2=2^{21}$ | 0.001 | mean propagation for all 16 layers: 0.016 u |
| one sample through the whole network | $L(2n^2-n)\approx2^{25}$ | 0.0156 per sample | $B$ = 65,536 samples; the floor = 6,554 |
| $r$ vectors through one layer, $(n,n)@(n,r)$ | $2n^2r$ | $r/n$ | 6.4 u/layer buys ≈ 6,500 vectors per layer |
| one elementwise pass over an $n\times n$ array | $n^2$ | 0.0005 | $\approx$ 13,000 passes per layer at the floor |
| one transcendental pass over $n\times n$ (exp/log 16, norm.cdf 48 f32 / 96 f64 billed per element) | $16$–$96\,n^2$ | 0.008–0.047 | pairwise nonlinear statistics are cheap |
| per-neuron quadratic forms $s_i^2=w_i^\top\Sigma w_i$, dense $\Sigma$ | $\approx2n^3$ | 1.0 dense, 0.556 Strassen L5 | unavoidable at the required accuracy (below) |
| congruence $W^\top\Sigma W$ (propagate an $n\times n$ symmetric object) | $\approx2$–$4n^3$ | 1.03 (Strassen sym3) – 2.0 dense | "one $n\times n$ object carried" |
| one elementwise pass over an $n\times n\times n$ array | $n^3=2^{30}$ | 0.5 per FLOP/element | a dense f32 $n^3$ array is 4 GiB, the grader's single-array cap |
| contract a dense third-order tensor with all $w_i^{\otimes3}$ | $\ge2n^4$ | $\ge1024$ | $=B$: impossible |
| transport a dense third-order tensor by $W$ in all legs | $\approx6n^4$ | $\approx3072$ | $3B$ per layer: impossible |

Other binding limits (**Fact**, costmodel and starter kit): 0.4 s of Python-side residual time at ≈ 0.022 ms per metered call on the grader, i.e. **≤ ≈ 18,000 metered calls per network, ≈ 1,100 per layer**; a Strassen-L5 product family costs ≈ 113 calls; 8 GB of memory (an $n\times n$ float32 array is 4 MB); 120 s of wall time on 2 vCPUs.

**The least work any method must do.**

1. *Reading the weights* is free in FLOPs (0.016 u for all layers) but every per-network (quenched) quantity depends on them.
2. *Per-neuron variances.* The readout needs $s_i^2$ at the final layer to ≈ 0.14 % (Section 3.2), and every layer's $s_i^2$ to comparable accuracy because their errors propagate (Section 3.5). The annealed value $\tfrac2n\mathrm{tr}\,\Sigma$ misses the quenched $s_i^2$ by the relative amount $\sqrt{2/r_{\rm eff}}$ (7 % at $a_2$ rising to 21 % at $a_{15}$, measured at n = 1024), two orders above the tolerance, and truncating $\Sigma$ to rank $k$ leaves a per-neuron error of relative size $\approx\sqrt2\,\|\Sigma-\Sigma_k\|_F/\mathrm{tr}\,\Sigma$, which must be $\lesssim10^{-3}$, i.e. the residual must hold $\lesssim$ 0.5 % of the Frobenius norm. No compression of $\Sigma$ that is blind to $W_{l+1}$ achieves this unless its spectrum decays fast, which it does not (at n = 1024 a 0.1 % truncation needs rank ≥ 0.93, 0.81, 0.68, 0.55, 0.47 n at $a_2,a_4,a_8,a_{12},a_{15}$, and 1 % needs 0.45, 0.30, 0.19, 0.15, 0.12 n). So **the quadratic forms with a dense $n\times n$ second-order object are the floor: ≈ 0.56–1 unit per layer, 9–16 units for the network (9–16 % of the floor budget).**
3. *Second-order state.* If the state carries $\Sigma_l$ (or any $n\times n$ symmetric object) it pays ≈ 1–2 units per layer to transport it, ≈ 16–32 units in all; the gate step on such an object is elementwise and costs 0.01–0.1 units.
4. *Beyond second order.* Explicit third- and fourth-order joint objects are excluded by factors of 100–500 (rows 8–10 of the table). Since every index pattern of $K_3$ and $K_4$ carries leading-order information (Section 3.3), higher-order joint content must be carried implicitly: in a form whose contractions with $w_i^{\otimes k}$ cost $O(n^2)$ per neuron, i.e. $O(n^3)$ per layer, generated from objects of size $O(n^2)$.
5. *Sampling as a floor comparison.* Plain sampling has per-neuron error variance 0.0748/N (**Fact**); raw $10^{-8}$ needs $N=7.5\times10^6$ samples = 117,000 units = 114 B. At the floor ($N=6{,}554$) its raw MSE is $1.1\times10^{-5}$, 1,100 times the target. A sampling component used as a correction must therefore remove ≥ 99.9 % of the per-sample variance of $a_{L,i}(x)$, i.e. come with a surrogate that matches the network function pointwise to ≈ 3 % rms in $L^2(\gamma_n)$.

**What this implies within $B=1024$.** At the floor each layer has 6.4 units ≈ 6 dense products ≈ 11 Strassen products, and ≈ 1,100 metered calls. A covariance-level computation (carrying one $n\times n$ second-order object by the congruence $W^\top\Sigma W$, whose diagonal is the set of quadratic forms) uses 1–2 units of that. Whatever a design carries beyond second order must fit in the remaining ≈ 4.4–5.4 units per layer: **about four or five carried $n\times n$ objects, or about 4,500–5,500 vectors pushed through $W_l$, per layer.** Above the floor the multiplier grows linearly: doubling the cost from 0.1 B to 0.2 B must halve the raw MSE to break even, so the floor regime is the design point.

---

## 5. How the difficulty scales with width and depth (part d)

The infinite-width theory supplies the $O(1)$ skeleton; the per-network corrections come in powers of $n^{-1/2}$; the accuracy target sits at a definite place in that hierarchy.

| order | quantities | depth dependence | measured |
|---|---|---|---|
| $O(1)$, self-averaging | per-neuron second moment $\mathbb Ez^2\approx2$ (He-critical); cross-input correlation $\rho_l$ ($\rho_1=0$, $\rho_{l+1}=f(\rho_l)$); share of the energy in the mean $f(\rho_l)$; distribution of $t=m/s$, $\mathcal N(0,\rho_l/(1-\rho_l))$; gains $f'(\rho_l)$ | $\rho$: 0, 0.32, 0.49, … 0.92 at $z_{16}$; $1-\rho_l\approx c/l^2$ | std of $t$ at n = 1024: 0.69, 1.24, 2.10, 3.02, 3.70 ($z_2,z_4,z_8,z_{12},z_{16}$) against 0.68, 1.24, 2.06, 2.78, 3.46 at infinite width; only 0.64, 1.08, 1.54, 1.79, 2.05 at n = 128 |
| $O(1)$, quenched | $m_i$, $t_i$, $\pi_i$, $\mu_{l,i}$, $p_{z_{l,i}}(0)$ | spread of $t$ grows: $\tau_l$ = 0.68 → 3.46 | — |
| $O(n^{-1/2})$ | correlations between distinct neurons; relative fluctuation of $s_i^2$ around $\tfrac2n\mathrm{tr}\Sigma$, $\sqrt{2/r_{\rm eff}}$; radius $\sqrt{\mathrm{Var}\,r}/\mathbb Er=1/\sqrt{2n}$ | $r_{\rm eff}/n$ falls with depth | $\sqrt n\times$ rms correlation 0.74 ($a_1$), 0.87 ($a_2$), 1.1–1.3 ($a_4$), 1.4 ($a_8$) at n = 32–128; $r_{\rm eff}/n$ = 0.42 → 0.044 from $a_2$ to $a_{15}$ at n = 1024, as at n = 32–256; $s_i^2$ spreads 7–21 % around its annealed value at n = 1024 |
| $O(n^{-1})$ | per-neuron $\lambda_3$ (zero ensemble mean, sign largely that of $t$) and $\lambda_4$ (coherent, positive; per-neuron spread about half its mean at n = 128); error of the Gaussian readout given exact $(m,s)$; entries of connected 3-point cumulants; radial spike $\mu\mu^\top/2n$ | grows with depth, roughly like $l/n$ | at n = 1024: $n\,\mathrm{rms}\,\lambda_3/l$ = 5.0 → 6.9 and $n\lambda_4/l$ = 5.0 → 3.6 from $z_2$ to $z_{16}$ (4.2–5.5 and 4.6–6.9 at n = 128); $n\times$(Gaussian-readout error) = 0.43–0.57 at n = 1024, 0.2–0.57 at n = 32–256 |
| $O(n^{-3/2})$ | entries of connected 4-point cumulants among distinct neurons; the unstructured (Frobenius-channel) part of each $\lambda_4$ | | |
| $O(n^{-2})$ | $\lambda_5,\lambda_6,\lambda_3^2,\lambda_3\lambda_4,\lambda_4^2$; error of the first-order readout | | $n^2\times$(first-order readout error) = 1–8 at n = 32–256, growing with depth; at n = 1024 below the sampling floor ($\le1.7\times10^{-5}$ at $z_{16}$) |

**Where the target sits (Measured).** The rms per-neuron error target $10^{-4}$ at $n=1024$ is $0.1/n$, or $\approx0.4\,(L/n)^2$. The error of the Gaussian readout given the exact quenched $(m_i,s_i)$ — a first-order-in-$1/n$ quantity, $(0.2$–$0.57)/n$ at every width — is $4.2\times10^{-4}$ rms on the width-1024 bench network (MSE $1.7\times10^{-7}$, 17 times the target), and the first-order readout's residual is ≤ $1.7\times10^{-5}$ (MSE ≤ $3\times10^{-10}$). So **the bar requires every first-order ($1/n$) effect to be computed to roughly 10–40 % relative accuracy, while second-order ($n^{-2}$) effects can be dropped at the readout** — but not necessarily in the propagation, where they accumulate over 15 layers with the gains of Section 3.5 and where coherent errors are not averaged.

**Depth.** Three things grow with depth: the non-Gaussianity of each neuron ($\mathrm{rms}\,\lambda_3\approx(5$–$7)\,l/n$ and $\lambda_4\approx(3.6$–$6)\,l/n$, measured at n = 128–1024 while $l/n\lesssim1/8$; 0.10 and 0.056 at the last layer at n = 1024), the anisotropy of the second-order structure (the mean-direction spike holds 0.6 % of $\mathrm{tr}\,\Sigma$ at $a_2$ and 11 % at $a_{15}$ at n = 1024, its eigenvalue 107 times the mean one there, ≈ 20 times what the input radius alone puts in it; $r_{\rm eff}/n$ falls accordingly), and the fraction of nearly deterministic gates (Section 2.5). The damping of mean errors weakens with depth ($f'(\rho_l)\to1$). The expansion parameter of the finite-width corrections is $L/n=1/64$ at the competition shape; the smoke-test shape (width 256, depth 32) has $L/n=1/8$, eight times larger: a design that leans on the asymptotics must still not fail there (the smoke test is not scored, but a failure is catastrophic).

---

## 6. Degrees of freedom any design must choose (part e)

Stated abstractly, without presupposing a representation. Each is a genuine choice: the problem's mathematics constrains it (the constraint is given) but does not fix it.

1. **The picture and the cut.** What travels: laws forward (push-forward by arrows), observables backward (pull-back of the readout), or both, paired at one or more cuts (Section 2.6; $\mu_{L,i}=\langle p_l,F_{l\to L,i}\rangle$ for every $l$). *Constraint:* the quenched information enters only through the $W_l$, linearly in each arrow's linear half; whichever object meets $W_l$ pays $n^2$ per vector or $n^3$ per $n\times n$ object.
2. **The state space at a layer.** Which coordinates on a law are kept: cumulants or moments, face statistics (masses and barycentres of faces, gate probabilities, wall densities), values of transforms, conditional or graphical structure, samples, or something else; at which index resolution (per neuron, per pair, higher); and in which representation (explicit, factored, generated on demand). *Constraints:* it must deliver, at the readout layer, each neuron's $(m_i,s_i^2,\lambda_{3,i},\lambda_{4,i})$ to the tolerances of Section 3.2; it must deliver the second-order structure of every layer to the tolerances of Sections 3.3–3.5 (≈ 0.1 % coherent, ≈ 0.6 % in Frobenius norm); it must contain leading-order information from every index pattern of the third and fourth joint cumulants without storing them (Section 3.3); and it must fit ≈ 6 units per layer (Section 4).
3. **The arrow on states, in two halves.** (a) *The gate half* — the map from the state of $z_l$ to the state of $a_l=\mathrm{relu}(z_l)$ — is where every approximation of the problem lives: it is not closed on any finite family, and its exact sensitivities are face statistics (Sections 2.6, 3.1, 3.4). (b) *The linear half* — the transport by $W_{l+1}$ — is exact and multilinear, and sets the cost: $O(n^3)$ per carried $n\times n$ object, $O(n^2)$ per carried vector.
4. **The readout.** The map from the final state to $\mathbb E\,\mathrm{relu}(z_{L,i})$. Its exact derivatives are the gate probability and the wall density and its derivatives (Section 3.1); at $n=1024$ the per-neuron first non-Gaussian order is necessary and sufficient (Section 3.2).
5. **Error control.** How the design keeps (a) the scale (dilation) component of its state error below $10^{-4}$ at every layer, since it is propagated with gain exactly one; (b) its coherent (trace-type, sign-definite) errors below the tolerances that the averaging over $n$ inputs does not relax; (c) which quantities it computes per network and which it takes from the ensemble (only self-averaging quantities can be ensemble constants, Section 2.3); (d) how it uses the exact identities as checks or constraints (dilation and Euler at every cut, the wall formula, the exactly known first layer, the thin-shell identity).
6. **Allocation over depth.** Errors injected at layer $l$ reach the output through gains that are small for random mean errors at shallow layers but equal to one for scale errors and of order one for coherent second-order errors at every depth (Section 3.5). The budget per layer and the accuracy per layer are therefore a design choice with a definite optimum, not a uniform split.

---

## 7. Messages to the design streams

Facts and identities from this note that each stream may find load-bearing. They are offered as constraints and questions, not as designs.

- **Faces and barycentres.** The answer is exactly linear in the face barycentres (2.1 iii) and exactly a sum of wall terms (2.1 iv); every sensitivity is a face statistic (2.6, 3.1, 3.4): gate probabilities, densities at walls and their slopes, densities at corners where two walls meet. A state built from face data has the right dual variables by construction; the question is which face data close under an arrow at $O(n^3)$.
- **Bethe and cavity on pseudorandom geometry.** Conditionally on a layer's law, each next-layer neuron sees an independent random projection (2.3). Connected joint cumulants among distinct neurons are tree-like, of size $n^{-(p-1)/2}$, and every index pattern matters at leading order (3.3): this is the regime where tree factorisations are natural and where storing the tree's leaves is not possible. The coherent part of the fourth cumulant is one global scalar per layer (the thin-shell excess), which a tree-local scheme must still get right.
- **Matchings, hafnians, signings.** The error-transfer table of 3.3 is Isserlis' theorem: errors reach the next layer through pairings, and the pairings that close inside one copy of the tensor (traces) are the coherent, unaveraged channels. Joint moments of ReLU outputs of a Gaussian layer expand, through the Hermite expansion of the ReLU, into sums over pairings of Hermite legs between neurons (the diagram formula); at layer 1 everything with at most three neurons is elementary and four-neuron orthant probabilities are not (2.2).
- **Tropical skeleton and temperature.** At depth a large fraction of gates is effectively deterministic (39 % with $|t|>3$ and 17 % with $|t|>4.75$ at layer 16 in the infinite-width limit; 43 % and 22 % measured at n = 1024), and on them the network is exactly linear; the Jensen gap $\tfrac12(\mathbb E|z|-|\mathbb Ez|)$ is the only nonlinearity, and it is controlled by the wall density (2.5, 3.1).
- **Heisenberg picture and Dirichlet forms.** $\mu_{L,i}=\langle p_l,F_{l\to L,i}\rangle$ at every cut; first-order sensitivities are the cumulant-weighted derivatives $\mathbb E[\partial^\alpha F]/\alpha!$; the second-order (trace) channel is the expected Laplacian of the pulled-back readout, $O(1)$ at every depth; at the input it equals the answer itself (Euler–Stein). The averaged propagator $\mathbb E[J_{l\to L}]$ is anisotropic, its top singular value of the order of a weight matrix's spectral norm.
- **Markov networks, CMI, exchange relations.** Layers are conditionally independent given the previous law (2.3); the gates are not independent across layers (the averaged propagator exceeds the gate-independent product by up to ×2.4 at n = 64 and ×2 at n = 128, 3.5); the coherent quantities (scale, traces, thin-shell excess) are global and must be carried exactly, whatever local structure a design exploits.

## 8. What this note does not settle

- Width extrapolations from small widths can be badly pre-asymptotic: the Gaussian closure's error falls as $n^{-0.82}$ between widths 64 and 128 but as $\approx n^{-2}$ between 128 and 1024 (BRIEF §4). The per-neuron statements of this note (readout errors, cumulant sizes, the bias distribution, the second-order structure) were therefore checked directly on one width-1024 bench network; the joint-structure ladders (3.3–3.4) rest on widths 16–128 with one or two networks per width and are quoted with a range of exponents.
- The depth responses (3.5) were measured at widths 32–128; their $n=1024$ values are estimated from the infinite-width theory where it applies (mean channel) and assumed $O(1)$ where it does not (trace channel). At the last cut the trace response is 0.84 rms at n = 1024 against 1.40 at n = 64, so the n = 64 tolerances of the trace channel are conservative by up to ≈ 1.7×.
- The one-step information test (3.4) isolates the error of the gate half at one layer given exact inputs; how such errors compound over 15 layers depends on the design, and is governed by the gains of 3.5.
- Nothing here ranks designs. The requirements and floors are the same for every design; which state and which arrow meet them is the streams' question.

## Appendix A. Numerical methods (reproducible from the code below)

All runs: pure numpy (float64 unless stated), `OPENBLAS_NUM_THREADS=1`, depth 16, $W_l$ drawn as `default_rng(1000+seed).standard_normal((n,n))*sqrt(2/n)` for $l=1..16$ in order, inputs from separate generators. "Replicas" are independent input streams for the same network; squared signals are estimated by cross-products of replicas, which removes the Monte Carlo noise bias.

- **A.1 Per-neuron laws** (Sections 2.3, 2.5, 3.1–3.2, 5): widths 32 (4 networks), 64 (2), 128 (1), each $10^7$ samples; width 256 (1 network, $8\times10^6$ samples); width 1024: bench network `w1024_d16` seed 7301001, weights regenerated bit-identically to the whest bake (float32 forward pass, float64 accumulation), $1.2\times10^6$ samples, its Monte Carlo means agreeing with the bench truth within the combined noise. Per layer and neuron: power sums of $z-c$ to order 6 ($c$ a pilot mean), $\sum\mathrm{relu}(z)$, $\sum1_{z>0}$, a window count for $p_z(0)$; for $a_l$ the covariance and the squared-norm moments. Two replicas by alternating chunks.
- **A.2 Readout sensitivities** (Section 3.1): width 64, 2 networks, $2\times10^6$ samples, 32 neurons at $z_2,z_4,z_8,z_{16}$; perturbations with known cumulants averaged analytically per sample.
- **A.3 Depth responses** (Sections 2.1, 3.5): averaged Jacobians $\mathbb E[J_{l\to L}]$ for all $l$ by reverse mode (width 64: $2\times10^5$ samples; width 128: $10^5$); second-order responses by antithetic isotropic noise with common random numbers (width 64: $6\times10^5$ samples; width 32: $10^6$).
- **A.4 Joint structure** (Sections 3.3–3.4): two passes over identical samples (exact centring), $4\times10^6$ samples per replica, two replicas per network; widths 16, 24, 32 (full $K_3$ and $K_4$ of $a_l$ at five layers), 64 (full $K_3$), 128 (pairwise slices only, $3\times10^6$ per replica).

<details><summary>Core code (numpy; about 120 lines)</summary>

```python
import numpy as np
from scipy.special import ndtr
phi = lambda u: np.exp(-0.5*u*u)/np.sqrt(2*np.pi)
R   = lambda u: u*ndtr(u) + phi(u)                 # E relu(u + G)
def net(n, seed, L=16):
    g = np.random.default_rng(1000+seed)
    return [g.standard_normal((n, n))*np.sqrt(2.0/n) for _ in range(L)]

# A.1: per-layer accumulation (one chunk x of shape (B, n)); S[l,k] += sum((z-c)^k), k<=6
def accumulate(W, x, c, S, Rs, Ps):
    a = x
    for l, Wl in enumerate(W):
        z = a @ Wl; y = z - c[l]; p = np.ones_like(y)
        for k in range(1, 7): p *= y; S[l, k] += p.sum(0)
        a = np.maximum(z, 0); Rs[l] += a.sum(0); Ps[l] += (z > 0).sum(0)

# Section 3.1: first- and second-order readouts from standardized cumulants
def readout(s, t, l3, l4, l5=0, l6=0, order=1):
    ak = lambda k: s*(-1)**k*np.polynomial.hermite_e.hermeval(t, np.eye(k-1)[k-2])*phi(t)
    r = s*(t*ndtr(t) + phi(t))
    if order >= 1: r += l3*ak(3)/6 + l4*ak(4)/24
    if order >= 2: r += l5*ak(5)/120 + (l6/720 + l3**2/72)*ak(6) + l3*l4*ak(7)/144 + l4**2*ak(8)/1152
    return r

# A.2: responses with known cumulants, averaged analytically per sample Z (array of z values), scale s
def responses(Z, s):
    e = 0.02*s; shift = (np.maximum(Z+e, 0) - np.maximum(Z-e, 0))/(2*e)              # -> P(z>0)
    v = (0.05*s)**2; gauss = (np.sqrt(v)*R(Z/np.sqrt(v)) - np.maximum(Z, 0))/(v/2)     # -> p(0)
    e3 = 0.15*s                                                                         # Y = e3(E-1), E ~ Exp(1)
    up = np.where(Z >= e3, Z, e3*np.exp(np.minimum(Z-e3, 0)/e3))                         # E relu(Z+Y)
    dn = np.where(Z+e3 > 0, Z + e3*np.exp(-np.maximum(Z+e3, 0)/e3), 0.0)                # E relu(Z-Y)
    skew = ((up - dn)/2)/(-(e3**3)/3)                                                   # -> p'(0)
    b = 0.15*s                                                                          # Laplace(b) vs N(0, 2b^2)
    kurt = (np.maximum(Z, 0) + 0.5*b*np.exp(-np.abs(Z)/b) - np.sqrt(2)*b*R(Z/(np.sqrt(2)*b)))/(b**4/2)  # -> p''(0)
    return [r.mean(0) for r in (shift, gauss, skew, kurt)]

# A.3: averaged Jacobians E[d a_L / d a_l] for all l at once (reverse mode), chunk x of shape (B, n)
def avg_jacobians(W, x, EJ):
    a = x; gates = []
    for Wl in W:
        z = a @ Wl; g = (z > 0)*1.0; a = z*g; gates.append(g)
    Jb = gates[-1][:, None, :]*np.eye(len(W[0]))[None]          # d a_L / d z_L
    for l in range(len(W)-1, -1, -1):
        JA = np.einsum('jk,bki->bji', W[l], Jb)                  # d a_L / d a_l  (l = 0: input)
        EJ[l] += JA.sum(0)
        if l > 0: Jb = gates[l-1][:, :, None]*JA

# A.4: bivariate Gaussian pieces of the information test of Section 3.4 (predictions of Cov(relu(z)), vectorised over pairs j<k)
def T00(mj, mk, Cjj, Ckk, Cjk, ny=400):                          # E[relu(zj) relu(zk)]
    yq, wq = np.polynomial.legendre.leggauss(ny); sj = np.sqrt(Cjj); lo = np.maximum(-mj/sj, -12.0)
    y = 0.5*(12-lo)[:, None]*yq + 0.5*(12+lo)[:, None]; w = 0.5*(12-lo)[:, None]*wq
    beta = Cjk/Cjj; v = Ckk - Cjk**2/Cjj; zeta = mj[:, None] + sj[:, None]*y
    u = (mk[:, None] + beta[:, None]*(zeta - mj[:, None]))/np.sqrt(v)[:, None]
    return (zeta*phi(y)*np.sqrt(v)[:, None]*R(u)*w).sum(1)
def Tdelta(p, h, mj, mk, Cjj, Ckk, Cjk):     # E[delta^(p)(zj) h(zk)] = (-1)^p d^p/dzeta^p [p_j(zeta) E(h(zk)|zj=zeta)] at 0
    sj = np.sqrt(Cjj); beta = Cjk/Cjj; sv = np.sqrt(Ckk - Cjk**2/Cjj); hz = 1e-3*sj
    def Psi(z):
        u = (mk + beta*(z - mj))/sv
        H = {'relu': sv*R(u), 'gate': ndtr(u), 'delta': phi(u)/sv}[h]
        return phi((z - mj)/sj)/sj*H
    f = [Psi(k*hz) for k in (-2, -1, 0, 1, 2)]
    if p == 0: return f[2]
    if p == 1: return -(f[0] - 8*f[1] + 8*f[3] - f[4])/(12*hz)
    return (-f[0] + 16*f[1] - 30*f[2] + 16*f[3] - f[4])/(12*hz**2)
# E[a_j a_k] (E4) = T00 + sum over |alpha| in {3,4} of kappa_alpha/alpha! * E[d^alpha (relu_j relu_k)], e.g.
#   + D21[j,k]/2 * Tdelta(0,'gate',j,k) + K22[j,k]/4 * Tdelta(0,'delta',j,k) + K31[j,k]/6 * Tdelta(1,'gate',j,k)
#   + kappa3(j)/6 * Tdelta(1,'relu',j,k) + kappa4(j)/24 * Tdelta(2,'relu',j,k)   (and the j<->k terms)

# A.4: index-pattern masks for a third-order tensor (fourth order analogously)
def patterns3(n):
    i, j, k = np.indices((n, n, n)); eq = (i == j)*1 + (j == k) + (i == k)
    return {'(3)': eq == 3, '(2,1)': eq == 1, '(1,1,1)': eq == 0}
# contribution of a pattern to kappa_3(z_i):  np.einsum('ijk,ia,ja,ka->a', K3*mask, W, W, W)
```
</details>

