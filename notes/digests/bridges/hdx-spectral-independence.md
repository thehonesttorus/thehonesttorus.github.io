# Bridge digest: high-dimensional expanders and spectral independence

### A Gibbs distribution is a weighted simplicial complex whose links are its pinnings; local spectral expansion of that complex is spectral independence; local-to-global theorems turn it into mixing. Where the "approximate local Markov property" sits in this picture.

*Bridge digest for prong 3, 2026-10-01. Written from the primary sources listed in §8, read in this session (full text or targeted page retrieval through alphaXiv; ECCC for Dinur–Kaufman). Labels follow programme rule P1.2: **Fact** (a source's own statement, with its theorem number), **Derivation** (proved here; the finite identities are checked numerically by [`check_hdx_si.py`](check_hdx_si.py), §9), **Interpretation** (my reading, not in the sources), **Conjecture**, **Guard**. Inputs used from the other notes, as vocabulary only: Note 1 §§4–6 ([conditional arrow algebra](../../conditional-arrow-algebra.md): Lüders conditioning, KMS/Gibbs states on histories, the Jenčová–Petz sufficiency criterion), Note 2 §1.1 ([simplicial complex](../../simplicial-complex-as-decomposition.md): the face algebra $D_0(K)\cong C(K)$, states are probability measures on faces), prong 3 §§2, 4, 5 ([local-to-global unlocks](../../local-to-global-unlocks.md), unlocks U2 and U6), prong 2 §3.2 ([MLP bridge](../../mlp-bridge.md): memory measured as a conditional mutual information). Nothing below refers to a particular realization.*

---

## 0. The answer in brief

1. **The bridge is a theorem.** For any probability law $\mu$ on $\{0,1\}^n$ (or $[q]^n$), Anari–Liu–Oveis Gharan build a pure, $n$-partite, weighted simplicial complex $X_\mu$ whose faces are partial configurations and whose **links are exactly the conditional laws (pinnings) of $\mu$**. The second eigenvalue of the link walk at a face $\tau$ equals $\lambda_{\max}(\Psi_{\mu^\tau})/(n-|\tau|-1)$, where $\Psi$ is the pairwise influence matrix (ALO Thm 1.11 and Thm 3.1, an *identity*, not an inequality). Hence "spectral independence of $\mu$ and all its pinnings" and "local spectral expansion of $X_\mu$" are the same statement (ALO Thm 1.5), and the Glauber dynamics is the top down-up walk of $X_\mu$. The local-to-global theorems of Kaufman–Mass, Dinur–Kaufman, Kaufman–Oppenheim and Alev–Lau then bound its gap (ALO Thm 1.3); the entropy versions (Chen–Liu–Vigoda; Anari–Jain–Koehler–Pham–Vuong) give optimal $O(n\log n)$ mixing. §§1–3.
2. **The common proof architecture is a chain rule along a filtration of conditional expectations** ("pin one more coordinate"): Garland's identities, the entropy decomposition, the variance martingale of Chen–Eldan and the trickle-down equation of Anari–Koehler–Vuong are the same law of total variance (or entropy, or covariance) applied level by level, with the link controlling the loss at each level (§2.8).
3. **What plays the "approximate local Markov property".** Not spectral independence. In the classical HDX encoding the exact Markov property of a Markov random field is *built in*: it is the statement that the link at a face depends only on the face's restriction to a separating boundary, and that the link at a separating face is a join (product) of the links of the components. It is used in three places (heredity of the class under links, the tree recursion that certifies spectral independence, and the "shattering" step of the optimal-mixing proof) and *not at all* in the local-to-global theorems, which hold for arbitrary $\mu$. Spectral independence plays the role of **strong spatial mixing / uniform clustering** (ALO call their bound "very similar to aggregate strong spatial mixing"). The classical scheme "(exact Markov) + (strong spatial mixing) ⇒ quasi-linear mixing", quoted by Chen–Rouzé, becomes in HDX language "(links are local and split at separators) + (every link expands) ⇒ global gap". For non-commuting quantum Gibbs states there are no pinnings, hence no links; the **approximate local Markov property (Chen–Rouzé Thm III.1; Yang, arXiv:2609.38007, Thm II.1) is what reconstructs an approximately local "link"**, and the time-averaged detailed-balanced Lindbladian with single-Pauli jumps on $A$ is precisely the quantum counterpart of *resampling $A$ from its link* (the classical block heat-bath $E_A$, which is one restriction–co-restriction step of the down-up walk). §4.
4. **Where noncommutative geometry enters, as far as these sources go:** at the level of operator algebras, not of spectral triples or $K$-theory. Exact Markov ⇔ a commuting square of conditional expectations (Derivation, §4.5); approximate Markov ⇔ a quantitative defect of commuting squares (Bardet–Capel–Rouzé's approximate tensorization) or a sufficiency defect (Petz recovery, Fawzi–Renner); this is the same circle of theorems as Note 1 §6 (Jenčová–Petz). §4.5.
5. **For the programme:** unlock U6 of the prong-3 note is a theorem, not a conjecture, once a state on the face algebra $D_0(K)$ is read as a law on $2^V$ (§6.2); the measured "memory" of prong 2 is exactly the relative-entropy distance of the history law from the Gibbs/KMS family of Note 1 (Derivation, §6.4); the Markov property is *not* needed by the local-to-global machinery, so the non-Markov history law is no obstruction to the HDX route (§6.4).

A factual note on the user's framing: arXiv:2609.38007 (T. H. Yang, "Improved estimate of local Markovianity for quantum Gibbs states", 29 Sep 2026) states in its AI usage statement that "the results presented in this paper were found by GPT-6 Astra Pro"; it does not attribute them to Claude.

---

## 1. Definitions

### 1.1 Weighted pure complexes, marginals, links

**Fact** (Alev–Lau [AL] §2.2, following Dikstein–Dinur–Filmus–Harsha; Kaufman–Oppenheim [KO] §2). A simplicial complex $X$ on a ground set $U$ is a downward-closed family of subsets (faces); $X(j)$ are the faces of dimension $j$, i.e. of size $j+1$, with $X(-1)=\{\emptyset\}$; $X$ is **pure of dimension $d$** if every face lies in a face of $X(d)$. A **weighted** pure complex $(X,\Pi)$ carries a probability $\Pi_d$ on $X(d)$ and the induced marginals
$$
\Pi_j(\alpha)=\frac1{j+2}\sum_{\beta\in X(j+1),\ \beta\supset\alpha}\Pi_{j+1}(\beta),
$$
equivalently: draw $\beta\sim\Pi_d$ and a uniformly random $j$-face of $\beta$. (KO and Oppenheim use an unnormalised weight $m$ with $m(\tau)=\sum_{\sigma\supset\tau,\ \sigma\in X(k+1)}m(\sigma)$, "balanced" weights; Oppenheim's homogeneous weight is $m(\tau)=(n-k)!\,\#\{\sigma\in X(n):\tau\subseteq\sigma\}$.)

**Fact** ([AL] §2.2). The **link** of $\alpha\in X(j)$ is $X_\alpha=\{\beta\setminus\alpha:\alpha\subseteq\beta\in X\}$, of dimension $d-j-1$, with the conditional weights
$$
\Pi^\alpha_l(\tau)=\Pr_{\beta\sim\Pi_{j+1+l}}[\beta=\alpha\cup\tau\mid\beta\supset\alpha]=\frac{\Pi_{j+l+1}(\alpha\cup\tau)}{\binom{|\alpha\cup\tau|}{|\alpha|}\Pi_j(\alpha)} .
$$
Its 1-skeleton $G_\alpha=(X_\alpha(0),X_\alpha(1),\Pi^\alpha_1)$ is a weighted graph. A $d$-dimensional complex is **$(d+1)$-partite** if $U=U_0\sqcup\dots\sqcup U_d$ with every top face meeting every part exactly once.

### 1.2 Local walks and local spectral expansion

**Fact** ([AL] §2.3). The local walk at $\alpha$ is $M_\alpha=D_\alpha^{-1}A_\alpha$ on $G_\alpha$, $M_\alpha(x,y)=\Pi^\alpha_1(\{x,y\})/(2\Pi^\alpha_0(x))$; it is self-adjoint in $\langle\cdot,\cdot\rangle_{\Pi^\alpha_0}$, has top eigenvalue $1$ and stationary law $\Pi^\alpha_0$. Put
$$
\gamma_j:=\max_{\alpha\in X(j)}\lambda_2(M_\alpha),\qquad j=-1,\dots,d-2 .
$$
$X$ is a **(one-sided) $\gamma$-local spectral expander** if $\gamma_j\le\gamma$ for all $-1\le j\le d-2$ (Def. 1.1 of [AL], after Oppenheim and [KO]). **Two-sided** ([KO] Def. 1.2; Dinur–Kaufman's "$\lambda$-HD expander", [DK] Def. 1.4): additionally the smallest eigenvalue $\nu_\tau\ge-\gamma$.

**Fact** ([KO] §1.5, Cor. 4.7). Oppenheim's own definition only constrains the top links ($\dim\alpha=d-2$) and requires all links of dimension $>0$ to be connected. The two definitions agree up to the change $\lambda'=\lambda/(1+(n-1)\lambda)$: if $\mu_{n-1}\le\lambda/(1+(n-1)\lambda)$ on the top links (and links are connected) then $X$ is a one-sided $\lambda$-local spectral expander.

**Fact** (Oppenheim [O14] Def. 1.1, phrased with Laplacians $\Delta^+_{\tau,0}=I-M_\tau$). For $\lambda>\frac{n-1}n$: $X$ has **$\lambda$-local spectral expansion** if $X$ and all its links of dimension $>0$ are connected and every 1-dimensional link has Laplacian spectral gap $\ge\lambda$.

**Guard (partiteness).** In a $(d+1)$-partite complex the link graphs are partite; [DK] §1.2 records that the LSV Ramanujan complexes have link eigenvalue $-1/d$, so only one-sided expansion is available; ALO Claim 3.2 and Remark 3.4 record the general fact for partite complexes: $P_\emptyset$ has the eigenvalue $-\frac1{d-1}$ with multiplicity at least $d-1$, the part indicators $\mathbb 1_{U_i}$ being eigenvectors of $P_\emptyset-\frac d{d-1}\mathbb 1\pi^\top$ with that eigenvalue. Every complex $X_\mu$ of §3 is partite, so only one-sided theorems apply to it.

### 1.3 Up and down operators; higher-order walks

**Fact** ([AL] §2.4). For $f\in\mathbb R^{X(j)}$, $g\in\mathbb R^{X(j+1)}$:
$$
[U_jf](\beta)=\frac1{j+2}\sum_{x\in\beta}f(\beta\setminus x),\qquad
[D_{j+1}g](\alpha)=\sum_{\beta\supset\alpha,\ \beta\in X(j+1)}\frac{\Pi_{j+1}(\beta)}{(j+2)\Pi_j(\alpha)}\,g(\beta),
$$
adjoint to each other: $\langle g,U_jf\rangle_{\Pi_{j+1}}=\langle D_{j+1}g,f\rangle_{\Pi_j}$. The **down-up** walk $P^\triangledown_j=U_{j-1}D_j$ and the **up-down** walk $P^\triangle_j=D_{j+1}U_j$ are positive semi-definite, self-adjoint in $\Pi_j$, share their non-zero spectrum, so $\lambda_2(P^\triangledown_{j+1})=\lambda_2(P^\triangle_j)$; the **non-lazy up-down walk** is $P^\wedge_j=\frac{j+2}{j+1}\big(P^\triangle_j-\frac1{j+2}I\big)$. On the bipartite incidence graph $X(j)$–$X(j+1)$ with edge weights $\frac1{j+2}\Pi_{j+1}(\beta)$, the square of the walk matrix is $\mathrm{diag}(P^\triangledown_{j+1},P^\triangle_j)$. Longer walks: $P^\triangle_{a,b}=D_{a+1}\cdots D_b\,U_{b-1}\cdots U_a$.

(Notation warning: [AL] name $U_j$ and $D_{j+1}$ by their action on functions; as moves of the walker, $U_j$ goes *down* and $D_{j+1}$ goes *up*, [AL] Remark 2.7.)

**Fact** (Chen–Liu–Vigoda [CLV] §1.2, §5.1.1, sizes instead of dimensions). For $0\le r<s\le n$ the **order-$(s,r)$ down-up walk** $P^\vee_{s,r}$ on $X(s)$ removes $s-r$ uniformly random elements and re-adds $s-r$ elements from the link distribution; $P^\vee_{s,r}=(P^\downarrow_s\cdots P^\downarrow_{r+1})(P^\uparrow_r\cdots P^\uparrow_{s-1})$, a product of two mutually adjoint operators (Fact A.6). For a law on $[q]^V$ the order-$(n,n-1)$ walk is the Glauber dynamics; the order-$(n,n-\ell)$ walk is the **heat-bath block dynamics** that resamples a uniformly random set of $\ell$ sites.

### 1.4 Cochains and Laplacians (for Garland)

**Fact** ([O14] §§3–5). With inner products on $C^k(X,\mathbb R)$ weighted by $m$, the differential $d_k$, its adjoint $\delta_k$, the upper and lower Laplacians $\Delta^+_k=\delta_{k+1}d_k$, $\Delta^-_k=d_{k-1}\delta_k$; $\tilde H^k(X,\mathbb R)\cong\ker\Delta^+_k\cap\ker\Delta^-_k$. For a face $\tau$ the localisation $\varphi_\tau$ of a cochain is its restriction to the link $X_\tau$. In the 1-dimensional link of a $(n-2)$-face with the homogeneous weight, $\Delta^+_{\tau,0}$ is the usual normalised graph Laplacian.

---

## 2. Local-to-global theorems

### 2.1 Garland's method

**Fact** ([O14] Lemma 5.4, "local to global"). Let $X$ be pure $n$-dimensional, weighted, all links of dimension $>0$ connected, $0\le k\le n-1$. If $\bigcup_{\tau\in\Sigma(k-1)}\mathrm{Spec}(\Delta^+_{\tau,0})\setminus\{0\}\subseteq[\lambda,\kappa]$, then for every $\varphi\in C^k(X,\mathbb R)$
$$
(k+1)\|\varphi\|^2\Big(\kappa-\frac k{k+1}\Big)-\kappa\|\delta\varphi\|^2\ \ge\ \|d\varphi\|^2\ \ge\ (k+1)\|\varphi\|^2\Big(\lambda-\frac k{k+1}\Big)-\lambda\|\delta\varphi\|^2 .
$$
*Proof architecture* (as written there): split $\varphi_\tau$ into its constant part $\Delta^-_{\tau,0}\varphi_\tau$ and the rest; apply the link gap $\lambda\|(\varphi_\tau)^1\|^2\le\|d_\tau\varphi_\tau\|^2\le\kappa\|(\varphi_\tau)^1\|^2$ in every link; sum over $\tau\in\Sigma(k-1)$ using the localisation identities $\sum_\tau\|\varphi_\tau\|^2\propto\|\varphi\|^2$, $\sum_\tau\|d_\tau\varphi_\tau\|^2\propto\|d\varphi\|^2$ and $\|\Delta^-_{\tau,0}\varphi_\tau\|^2=\|\delta_{\tau,0}\varphi_\tau\|^2$, which sum to $\|\delta\varphi\|^2$.

**Fact** ([O14] Cor. 5.6; this is Garland's vanishing theorem in Oppenheim's normalisation). If moreover $\kappa\ge\lambda>\frac k{k+1}$ then $\tilde H^k(X,\mathbb R)=0$, $C^k=\ker\Delta^+_k\oplus\ker\Delta^-_k$ and $\mathrm{Spec}(\Delta^+_k)\setminus\{0\}\subseteq[(k+1)\lambda-k,\ (k+1)\kappa-k]$. (On $\ker\delta$ the lower bound reads $\|d\varphi\|^2\ge((k+1)\lambda-k)\|\varphi\|^2>0$: no harmonic cochains.) With $k=d-1$ this is the bound $\lambda^{(d-1)}(X)\ge1+d\varepsilon-d$ quoted in the prong-3 note from Lubotzky's Thm 2.5.

**Fact** ([O14] Thm 5.8, "very local to very global"). With $f(x)=2-\frac1x$ and $f^j$ its iterates: if the top links ($\Sigma(n-2)$) have non-zero spectrum in $[\lambda,\kappa]$ with $\lambda>\frac{n-1}n$, then for all $0\le k\le n-1$: $\tilde H^k(X,\mathbb R)=0$ and $\mathrm{Spec}(\Delta^+_k)\setminus\{0\}\subseteq[(k+1)f^{n-1-k}(\lambda)-k,\ (k+1)f^{n-1-k}(\kappa)-k]$.

### 2.2 Trickle-down

**Fact** ([O14] Lemma 5.1, "descent in links"). If $\bigcup_{\sigma\in\Sigma(k)}\mathrm{Spec}(\Delta^+_{\sigma,0})\setminus\{0\}\subseteq[\lambda,\kappa]$ then $\bigcup_{\tau\in\Sigma(k-1)}\mathrm{Spec}(\Delta^+_{\tau,0})\setminus\{0\}\subseteq[2-\frac1\lambda,\ 2-\frac1\kappa]$.

**Fact** (the same, in walk eigenvalues: [AL] Thm 1.3 and Thm 2.5; [KO] Lemma 4.5). If $\lambda_2(G_\beta)\le\gamma\le\frac12$ for every $k$-face and $G_\alpha$ is connected for every $(k-1)$-face, then $\lambda_2(G_\alpha)\le\frac\gamma{1-\gamma}$ for every $(k-1)$-face; [KO] add the two-sided half $\nu_k\ge\nu_{k+1}/(1-\nu_{k+1})$. Iterated ([AL] Cor. 1.4, 2.6; [KO] Cor. 4.6): if all links are connected,
$$
\gamma_j\le\frac{\gamma_{d-2}}{1-(d-2-j)\gamma_{d-2}} .
$$
Two reasons this is the step that makes the local certificate checkable ([AL] §1): the top links have edge weights in $\{0,1\}$ under the uniform weight, and in the regime $\gamma=O(1/d^2)$ the descent is lossless.

**Fact** (trickle-down as a covariance identity: Anari–Koehler–Vuong [AKV] Eq. (8), (12), Thm 101). For the pinning filtration of a law $\nu_0$ on $\binom{[n]}k$ ($\nu_t$ = posterior after observing $t$ uniformly random elements of $S\sim\nu_0$), the law of total covariance gives
$$
\mathrm{cov}(\nu_t)=\mathbb E[\mathrm{cov}(\nu_{t+1})\mid\mathcal F_t]+\frac1{k-t}\,\mathrm{cov}(\nu_t)N_t^{-1}\mathrm{cov}(\nu_t),
$$
$N_t$ the diagonal matrix of link marginals. Oppenheim's theorem in this language (Thm 101): if a.s. $\Sigma_1\preceq C\,N_1$ and the $1\leftrightarrow k$ up-down walk has $\lambda_2<1$, then $\Sigma_0\preceq\frac{C(k-2)}{k-1-C}N_0$. Since $C$-spectral independence of a $k$-homogeneous law is $\mathrm{cov}(\nu)\preceq Ck\,\Pi$ ($\Pi=\mathrm{diag}(\nu(i)/k)$; [AKV] Lemma 98, from ALO), trickle-down *is* the backward induction of spectral independence through the levels of the filtration.

### 2.3 Kaufman–Mass

**Fact** (Kaufman–Mass [KM], Def. 1.3, 1.4, 1.6, 1.7; Thms 2.1, 3.1, 3.2). The walk on dimension $i$ moves from an $i$-face to a random $(i+1)$-face containing it and then to a *different* $i$-face inside it. $X$ is an **$\alpha$-skeleton expander** if for every face $\sigma$ (including $\emptyset$) and every $S\subseteq X_\sigma(0)$, $\|E(S)\|_\sigma\le\|S\|_\sigma(\|S\|_\sigma+\alpha)$ (link 1-skeletons look random up to $\alpha$). $X$ is an **$\epsilon$-colorful expander** if every $i$-cochain $W$ with $0<\|W\|\le\frac12$ has $\|\psi(W)\|/\|W\|\ge\epsilon$, where $\psi(W)$ are the $(i+1)$-faces that meet $W$ but are not covered by it. Thm 2.1: if $\alpha<(\sqrt[2^d]2-1)/\sqrt2$ then $X$ is $\epsilon$-colorful with $\epsilon=\big((\sqrt[2^d]2-1-\sqrt2\alpha)/(2\sqrt2 d)\big)^d$. Thm 3.1: an $\epsilon$-colorful expander has all high-order walks $\mu$-rapidly mixing with $\mu=1-\frac{\epsilon^2}{2(d+1)^2}$ (through conductance $\Phi(G_i)\ge\epsilon/(i+2)$, Lemma 3.4, and Cheeger). Thm 4.1: $q$-thick Ramanujan complexes with $q>q_0(d)$ qualify. (The radicals are my reading of the extracted typography "$2^d\sqrt2$" in the PDF text.) [KM] §1.4 explain why Garland alone does not give this: Garland bounds the spectrum orthogonal to coboundaries, not orthogonal to the constants.

### 2.4 Dinur–Kaufman

**Fact** (Dinur–Kaufman [DK], ECCC TR17-089). Def. 1.1: a $d$-dimensional complex is a **$c$-agreement expander** if its underlying graph is connected and some test distribution $D$ on pairs of faces has $\Upsilon(X,D)=\inf_f\mathrm{disagree}_D(f)/\mathrm{dist}(f,\mathrm{Global})>c$; i.e. $\mathrm{agree}_D(f)>1-\varepsilon$ implies a global $g$ with $\Pr_s[f_s=g|_s]>1-\varepsilon/c$. Def. 1.4: **$\lambda$-HD expander** = every link of dimension $\ge1$ has 1-skeleton with all non-trivial normalised eigenvalues in $[-\lambda,\lambda]$ (two-sided). Thm 1.7: on a $d$-dimensional $\lambda$-HD expander the one-up walk $A_{k,k+1}$ on $k$-faces has $\lambda_2\le1-\frac1{k+1}+O(k\lambda)$ (the complete complex has $1-\frac1{k+1}-o_n(1)$). Thm 1.8: $\lambda_2(A_{k,k+t})\le\frac{k+1}{t+k+1}+O(tk\lambda)$. Thm 1.9 ("double sampler"): three-layer incidence graphs $[n]$–$V\subset\binom{[n]}k$–$W\subset\binom{[n]}d$ with $|V|+|W|+|E|=O(n)$ and $\lambda(G(U,V))^2\le\frac1k+\gamma$, $\lambda(G(V,W))^2\le\frac kd+\gamma$. Thm 5.1 (the form of Thm 1.6 proved): if $X$ is a $d$-dimensional $\lambda$-HD expander with $k^2<d$ and $\lambda<1/d$, then for the test "choose $r\in X(2k)$, then $s_1,s_2\subset r$ independently", $\mathrm{agree}(f)>1-\varepsilon$ implies $\Pr_s[f_s\equiv g|_s]>1-O(\varepsilon)$ with $g$ the pointwise majority. Lemma 1.5: explicit bounded-degree $\lambda$-HD expanders exist for every $\lambda,d$ (low-dimensional skeleta of LSV complexes).

*Proof architecture* ([DK] §1.3). (i) **Decreasing differences**: study the variance of a function pushed down through all dimensions simultaneously; the differences between successive levels are themselves $\lambda$-approximately decreasing (the authors credit a discussion with R. Eldan "regarding martingales"); this yields Thm 1.7. (ii) Write $A_{k,k+t}=B^\dagger B$ with $B$ a product of one-step incidence operators and telescope the singular values (Thm 1.8). (iii) Reduce agreement on $X$ to agreement on the complete complexes inside each top face (Dinur–Steurer), using the double-sampling property twice.

**Fact** ([DK] §1.6, their words). "Agreement expansion is a kind of approximate cohomology with local coefficients": exact agreement on overlaps is the sheaf condition, a global section; agreement testing studies "approximate sections".

### 2.5 Kaufman–Oppenheim: beyond the spectral gap

**Fact** ([KO] Thm 1.4 = Thm 5.4). On a pure $n$-dimensional one-sided $\lambda$-local spectral expander, for every $k\le n-1$ and $\varphi\in C^k_0$ (orthogonal to constants), $\|M^+_k\varphi\|\le\big(\frac{k+1}{k+2}+\frac{k+1}2\lambda\big)\|\varphi\|$. Thm 1.3/5.2 (one-sided decomposition): there are $j$-cochains $\varphi_j$, $0\le j\le k$, with $\|\varphi\|^2=\sum_j\|\varphi_j\|^2$ and $\langle M^+_k\varphi,\varphi\rangle\le\sum_j\big(\frac{k+1-j}{k+2}+f(k,j)\lambda\big)\|\varphi_j\|^2$; the identity behind it is $(d_k)^*d_k=(k+2)M^+_k$ together with a Garland-type formula for $\|d_k\varphi\|^2$ in terms of link operators. Thm 1.5/5.9 (two-sided): with $U^j_k=V^j_k\cap(V^{j-1}_k)^\perp$, $V^j_k=d_{j\nearrow k}(C^j_0)$ (functions on $k$-faces induced from $j$-faces), $C^k_0=\bigoplus_jU^j_k$ orthogonally, and for small $\lambda$, $\mathrm{Spec}(M^+_k)\subseteq\{1\}\cup\bigcup_j[\frac{k+1-j}{k+2}\pm\frac{\sqrt{k+1}}{k+2}\varepsilon_k]$. The abstract names the obstruction: the gap of high-order walks "is inherently small, due to natural obstructions (called coboundaries)": the slow eigenvectors are functions lifted from lower dimensions.

### 2.6 Alev–Lau: the product formula

**Fact** ([AL] Thm 1.5 = Thm 3.1). For a pure $d$-dimensional weighted complex and $0\le k\le d$,
$$
\lambda_2(P^\triangledown_k)=\lambda_2(P^\triangle_{k-1})\le1-\frac1{k+1}\prod_{j=-1}^{k-2}(1-\gamma_j).
$$
Thm 3.2: $P^\triangle_k$ has at most $|X(r)|$ eigenvalues above $1-\frac1{k+2}\prod_{j=r}^{k-1}(1-\gamma_j)$. Cor. 3.4: if $\gamma_{k-2}\le\frac1{k+1}$ and links are connected, $\lambda_2(P^\triangledown_k)\le1-\frac1{(k+1)^2}$. Prop. 3.3 (tightness): $\lambda_2(P^\triangledown_d)\ge1-\frac2{d+1}$ when $2(d+1)\le|X(0)|$, so an $O(1/k)$-local spectral expander has optimal gap up to constants. Cor. 3.5: $\lambda_2(P^\triangle_{a,b})\le(1+\gamma)^{b-a}\frac{a+1}{b+1}$, improving [DK]'s $\frac{a+1}{b+1}+O(a(b-a)\gamma)$.

*Proof architecture* ([AL] §3.1). **Garland identities** (Lemma 3.7, from [KO] and DDFH): for $f\in\mathbb R^{X(j)}$,
$$
\langle f,f\rangle_{\Pi_j}=\mathbb E_{\alpha\sim\Pi_{j-1}}\|f_\alpha\|^2_{\Pi^\alpha_0},\quad
\langle f,P^\triangledown_jf\rangle_{\Pi_j}=\mathbb E_{\alpha}\|J_\alpha f_\alpha\|^2_{\Pi^\alpha_0},\quad
\langle f,P^\wedge_jf\rangle_{\Pi_j}=\mathbb E_{\alpha}\langle f_\alpha,M_\alpha f_\alpha\rangle_{\Pi^\alpha_0},
$$
$f_\alpha$ the restriction to $\{\alpha\cup x\}$, $J_\alpha$ the projection on constants. Since $M_\alpha-J_\alpha\preceq\lambda_2(M_\alpha)I$ on the complement of constants, **Lemma 3.6**: $P^\wedge_k-P^\triangledown_k\preceq\gamma_{k-1}(I-P^\triangledown_k)$, i.e. $P^\wedge_{j+1}\preceq\gamma_jI+(1-\gamma_j)P^\triangledown_{j+1}$. Induction on $k$ using that $P^\triangledown$ and $P^\triangle$ share non-zero spectrum, and $P^\triangle_{j+1}=\frac{j+2}{j+3}P^\wedge_{j+1}+\frac1{j+3}I$, gives the product.

### 2.7 Variance and entropy contraction

**Fact** ([CLV] Def. 5.3, Thm 5.4; Guo–Mousa obtained it independently, [CLV] Remark 5.5; the variance version is [CLV] Thm A.9, also Kaufman–Mass). $(X,w)$ satisfies **$(\alpha_0,\dots,\alpha_{n-2})$-local entropy contraction** if for every global $f^{(n)}\ge0$, every $k$ and every $\tau\in X(k)$, $\mathrm{Ent}_{\pi_{\tau,2}}(f^{(2)}_\tau)\ge(1+\alpha_k)\mathrm{Ent}_{\pi_{\tau,1}}(f^{(1)}_\tau)$. Then with $\Gamma_0=1$, $\Gamma_i=\prod_{j<i}\alpha_j$, $(X,w)$ has order-$(k,n)$ **global entropy contraction** $\mathrm{Ent}_{\pi_k}(f^{(k)})\le(1-\kappa)\mathrm{Ent}_{\pi_n}(f^{(n)})$ with
$$
\kappa=\frac{\sum_{i=k}^{n-1}\Gamma_i}{\sum_{i=0}^{n-1}\Gamma_i}.
$$
Fact A.8: local variance contraction with $\alpha_k$ is local spectral expansion with $\zeta_k=\frac{1-\alpha_k}{1+\alpha_k}$. Thm 5.6: if $(X,w)$ is $(b_0,\dots,b_{n-1})$-**marginally bounded** ($\pi_{\tau,1}(i)\ge b_k$ for $\tau\in X(k)$) and a $(\zeta_0,\dots,\zeta_{n-2})$-local spectral expander, then local entropy contraction holds with
$$
\alpha_k=\max\Big\{1-\frac{4\zeta_k}{b_k^2(n-k)^2},\ \frac{1-\zeta_k}{4+2\log\frac1{2b_kb_{k+1}}}\Big\}.
$$
Fact 5.2: global entropy contraction at rate $\kappa$ gives, for $P^\vee_{s,r}$ and $P^\wedge_{r,s}$: Poincaré constant $\kappa$, MLSI constant $\kappa$, entropy decay rate $\kappa$, $T_{\mathrm{mix}}\le\lceil\frac1\kappa(\log\log\frac1{\pi^*_s}+\log\frac1{2\varepsilon^2})\rceil$, and Gaussian concentration $\Pr[|f-\mu f|\ge a]\le2e^{-\kappa a^2/2c^2}$ for $c$-Lipschitz $f$.

### 2.8 The proof architecture, stated once

All the local-to-global proofs above run on one identity. Four forms of it:

| form | identity | source |
|---|---|---|
| Garland | global quadratic forms of $f$ on $X(j)$ are $\Pi_{j-1}$-averages of link quadratic forms of the restrictions $f_\alpha$ | [AL] Lemma 3.7; [O14] Lemma 5.4 |
| entropy chain rule | $\mathrm{Ent}_{\pi_k}(f^{(k)})=\mathrm{Ent}_{\pi_j}(f^{(j)})+\sum_{\tau\in X(j)}\pi_j(\tau)\mathrm{Ent}_{\pi^{\tau,k-j}}(f^{(k-j)}_\tau)$ | [CLV] Lemma 5.1 (from Cryan–Guo–Mousa) |
| variance martingale | $\mathrm{gap}(P)=\inf_\varphi\mathbb E[\mathrm{Var}_{\nu_\tau}\varphi]/\mathrm{Var}_\nu\varphi$ for the chain $P_{x\to y}=\mathbb E[\nu_\tau(x)\nu_\tau(y)/\nu(x)]$ attached to a localization process; for linear tilts $\mathbb E[\mathrm{Var}_{\nu_{t+1}}\varphi\mid\nu_t]-\mathrm{Var}_{\nu_t}\varphi=-\langle v_t,C_tv_t\rangle\ge-\Vert C_t^{1/2}\mathrm{Cov}(\nu_t)C_t^{1/2}\Vert \,\mathrm{Var}_{\nu_t}\varphi$ | Chen–Eldan [CE] Prop. 19, Claim 22, Eq. (14)–(15) |
| total covariance | $\mathrm{cov}(\nu_t)=\mathbb E[\mathrm{cov}(\nu_{t+1})\mid\mathcal F_t]+\mathrm{cov}(\nu_t)\mathrm{cov}(Z_{t+1}\mid\mathcal F_t)\mathrm{cov}(\nu_t)$ | [AKV] Eq. (8) |

**Interpretation.** Each is the law of total variance (entropy, covariance) along the filtration "pin one more coordinate", i.e. along a tower of conditional expectations $\mathbb E[\cdot\mid\mathcal F_0]\leftarrow\mathbb E[\cdot\mid\mathcal F_1]\leftarrow\cdots$. The local hypothesis says that *one* step of the tower loses at most a fraction of the variance or entropy, the fraction being a link quantity (a second eigenvalue, an influence-matrix eigenvalue, an entropic-independence constant); the product of the per-level retentions is the global constant. Nothing global about $X$ is used beyond the tower. This is the precise form of the prong-3 slogan "restriction composed with co-restriction" (local-to-global unlocks §4, row 2), and §4.5 below is its operator-algebraic form.

---

## 3. A Gibbs distribution as a weighted simplicial complex

### 3.1 The encoding

**Fact** (Anari–Liu–Oveis Gharan [ALO] §1.1). For $\mu$ on $2^{[n]}$, $X_\mu$ is the pure $n$-dimensional (size convention: top faces have $n$ elements) $n$-partite complex on the ground set $\{1,\bar1,\dots,n,\bar n\}$ with parts $U_i=\{i,\bar i\}$: for each $S\in\mathrm{supp}\,\mu$ the top face $\sigma_S$ contains $i$ for $i\in S$ and $\bar i$ for $i\notin S$, with weight $w(\sigma_S)=\mu(S)$; lower faces by downward closure, $w(\tau)=\sum_{\sigma\supset\tau}w(\sigma)$. Its top down-up walk $P^\vee_n$ ("remove a uniformly random element, add one back with probability $\propto w$") **is exactly the Glauber dynamics** for $\mu$.

**Fact** ([CLV] §1.2, §5.1). For $\mu$ on $[q]^V$, $|V|=n$: ground set $\{(v,i)\}$, faces are partial configurations $(U,\tau)$, $\tau\in\Omega_U$; $\pi_k(U,\tau)=\binom nk^{-1}\mu(\sigma_U=\tau)$; $\pi_n=\mu$. The order-$(n,n-1)$ down-up walk is the Glauber dynamics. The same object appears as the **homogenization** $\mu^{\hom}$ of a law on a product space ([AJKPV-U] Def. 28): $\sigma\mapsto\{(\sigma_1,1),\dots,(\sigma_n,n)\}\in\binom\Omega n$.

### 3.2 Links are pinnings

**Fact** ([ALO] proof of Thm 1.5). "Conditioning on an element $i$ being 'in' corresponds exactly to taking the link of $X_\mu$ w.r.t. $i$. Similarly, conditioning on $i$ being 'out' corresponds exactly to taking the link w.r.t. $\bar i$." Generally the link of the face $(U,\tau)$ is the complex of the pinned law $\mu^\tau=\mu(\cdot\mid\sigma_U=\tau)$ on $V\setminus U$.

**Fact** ([CLV] proof of Claim 1.16). The local walk at a face has stationary law $\pi_{(U,\tau),1}((v,i))=\frac1{n-k}\mu(\sigma_v=i\mid\sigma_U=\tau)$: the vertex marginals of a link are the single-site marginals of the pinned law, divided by the number of free sites.

### 3.3 The influence matrix and spectral independence

**Fact** ([ALO] Def. 1.1, 1.2). $\Psi_\mu(i,j)=\Pr[j\mid i]-\Pr[j\mid\bar i]$ for $i\ne j$, $\Psi_\mu(i,i)=0$. $\mu$ is **$\eta$-spectrally independent** if $\lambda_{\max}(\Psi_\mu)\le\eta$, and **$(\eta_0,\dots,\eta_{n-2})$-spectrally independent** if every pinning of $i$ coordinates is $\eta_i$-spectrally independent. Always $\eta_i\le n-i-1$; product laws are $(0,\dots,0)$-SI; "if $\mu$ is $d$-homogeneous and all measures obtainable from $\mu$ by conditioning are negatively correlated, then $\mu$ is $(1,1,\dots,1)$-spectrally independent"; the law giving mass $\frac12$ to $\{1,\dots,\frac n2\}$ and to its complement has $\lambda_{\max}=n-1$. A row- or column-sum bound suffices: $\lambda_{\max}\le\max_i\sum_j|\Psi(i,j)|$.

**Fact** ([CLV] Def. 1.10, 1.11, multi-spin). $\Psi^\tau_\mu((u,i),(v,j))=\mu(\sigma_v=j\mid\sigma_u=i,\sigma_\Lambda=\tau)-\mu(\sigma_v=j\mid\sigma_\Lambda=\tau)$; $\eta$-SI means $\lambda_1(\Psi^\tau_\mu)\le\eta$ for every pinning. (Feng et al. use a total-variation version that implies it; [CLV] §1.1.)

**Fact** (equivalent forms). (a) [CE] Fact 23: for laws on $\{\pm1\}^n$ with $\Psi(\nu)_{ij}=\mathbb E[X_i\mid X_j=1]-\mathbb E[X_i\mid X_j=-1]$, $\Psi=\mathrm{Cov}(\nu)\,\mathrm{diag}(\mathrm{Cov}\,\nu)^{-1}$ and $\|\mathrm{Cor}(\nu)\|_{op}=\rho(\Psi(\nu))$, $\mathrm{Cor}=$ the correlation matrix; equivalently $\eta$-SI is $\mathrm{Cov}(\nu)\preceq(\eta+1)\mathrm{diag}(\mathrm{Cov}\,\nu)$ ([CE] Remark 37). (b) [AJKPV-U] Lemma 21: for $\mu$ on $\binom{[n]}k$ and $P=U_{1\to k}D_{k\to1}$, $\lambda_2(P)=\lambda_{\max}(\Psi^{cor}_\mu)/k$. (c) [AJKPV-U] Fact 20: $\eta$-SI iff the Hessian of $\log g_\mu(z^{1/\eta})$ at $\vec1$ is negative semi-definite, $g_\mu$ the generating polynomial (local concavity at one point).

### 3.4 The theorem: spectral independence is local spectral expansion

**Fact** ([ALO] Thm 1.11, Thm 3.1, Thm 1.5). For every law $\mu$ on $2^{[n]}$ the eigenvalues of $\Psi_\mu$ are real, and the spectrum of the 1-skeleton walk $P_\emptyset$ of $X_\mu$ is, as a multiset,
$$
\mathrm{spec}(P_\emptyset)=\mathrm{spec}\big(\tfrac1{n-1}\Psi_\mu\big)\ \cup\ \{-\tfrac1{n-1}\}^{\times(n-1)}\ \cup\ \{1\},
\qquad\text{so}\qquad
\lambda_2(P_\emptyset)=\tfrac1{n-1}\lambda_{\max}(\Psi_\mu)
$$
($\lambda_{\max}(\Psi_\mu)\ge0$ because the spectrum is real and $\mathrm{tr}\,\Psi_\mu=0$, so the trivial $-\frac1{n-1}$ never wins). Applied to every pinning: an $(\eta_0,\dots,\eta_{n-2})$-SI law gives a $\big(\frac{\eta_0}{n-1},\frac{\eta_1}{n-2},\dots,\frac{\eta_{n-2}}1\big)$-local spectral expander $X_\mu$. Multi-spin: [CLV] Claim 1.18 (= Thm 8 of Chen–Galanis–Štefankovič–Vigoda), $\zeta_k=\frac\eta{n-k-1}$.

*Proof architecture* ([ALO] §3, Claims 3.2–3.3, Remark 3.4). (i) $P_\emptyset(i,j)=\frac1{n-1}\Pr[j\mid i]$ etc., reversible w.r.t. $\pi(i)=\frac1n\Pr[i]$, $\pi(\bar i)=\frac1n\Pr[\bar i]$. (ii) **Partiteness produces trivial eigenvectors**: the part indicators $\mathbb 1_i=e_i+e_{\bar i}$ satisfy $Q_\emptyset\mathbb 1_i=-\frac1{n-1}\mathbb 1_i$ for $Q_\emptyset=P_\emptyset-\frac n{n-1}\mathbb 1\pi^\top$; "a generalization of the fact that the transition matrix of a bipartite graph always also has eigenvalue $-1$". (iii) Deflate them: $M_\emptyset=P_\emptyset-\frac n{n-1}\mathbb 1\pi^\top+\frac n{n-1}\sum_i\mathbb 1_i(\pi_i)^\top$ has the spectrum of $P_\emptyset$ with the $n$ trivial eigenvalues replaced by $0$. (iv) $M_\emptyset=\begin{pmatrix}A&-A\\B&-B\end{pmatrix}$ with $A(i,j)=\frac1{n-1}(\Pr[j\mid i]-\Pr[j])$, $B(i,j)=\frac1{n-1}(\Pr[j\mid\bar i]-\Pr[j])$, and $\frac1{n-1}\Psi_\mu=A-B$; $\det(xI-M_\emptyset)=x^n\det(xI-(A-B))$. The identity is checked numerically in §9 (C1).

**Interpretation.** The theorem has nothing to do with graphs or Gibbs measures: it is a statement about *any* law on a product of two-point sets and its conditionals. "Gibbs distribution" enters only when one asks how to *certify* the hypothesis (§3.9, §4.2).

### 3.5 What it yields

**Fact** ([ALO] Thm 1.3, Thm 1.6 = [AL], Remark 1.7). For an $(\eta_0,\dots,\eta_{n-2})$-SI law the Glauber dynamics has spectral gap
$$
\ge\ \frac1n\prod_{i=0}^{n-2}\Big(1-\frac{\eta_i}{n-i-1}\Big);
$$
if $(X,w)$ is $(\frac\alpha{d-1},\frac\alpha{d-2},\dots,\frac\alpha1)$-local spectral expander then $\lambda_2(P^\vee_d)\le1-\frac1{d^{1+\alpha}}$. (ALO state Alev–Lau with size indexing, $\lambda_2(P^\vee_d)\le1-\frac1d\prod_{k=0}^{d-2}(1-\alpha_k)$.) Mixing from any $\tau$: $t_\tau(\varepsilon)\le\frac1{1-\lambda^*}\log\frac1{\varepsilon\pi(\tau)}$ (Diaconis–Stroock, [ALO] Thm 2.3).

**Fact** (hardcore model, [ALO] Thm 1.8, Thm 1.13, Remark 1.10, Lemma 1.12). For $G$ of max degree $\le\Delta$ and $\lambda=(1-\delta)\lambda_c(\Delta)$, $\lambda_c(\Delta)=\frac{(\Delta-1)^{\Delta-1}}{(\Delta-2)^\Delta}$, the hardcore law is $(\eta_0,\dots)$-SI with $\eta_i\le\min\{C(\delta),\frac\lambda{1+\lambda}(n-i-1)\}$, $C(\delta)\le\exp(O(1/\delta))$; Glauber mixes in $O\big(((1+\lambda)n)^{1+C(\delta)}\log\frac1{\varepsilon\mu(\tau)}\big)$, with no dependence on $\Delta$; FPRAS for $Z_G(\lambda)$ up to the uniqueness threshold (Cor. 1.9), complementing Sly's hardness above it. The bound is on the total influence *on* a vertex, $\sum_u|\Psi(u,v)|\le C(\delta)$ (Thm 1.13); Chen–Liu–Vigoda later bounded the total influence *of* a vertex by $O(1/\delta)$ ([ALO] Remark 1.14).

The constant in §3.5 is sharp in shape but not in value: §9 (C5) computes, for an Ising chain on 6 sites, Glauber gap $0.066$ against the ALO bound $0.036$.

### 3.6 Optimal mixing: Chen–Liu–Vigoda

**Fact** ([CLV] Thm 1.12, Remark 1.13). If $G$ has max degree $\le\Delta$ and $\mu$ is a totally-connected Gibbs distribution of a spin system on $G$ that is **$b$-marginally bounded** and **$\eta$-SI**, then Glauber satisfies the MLSI with constant $\frac1{C_1n}$, $C_1=(\Delta/b)^{O(\eta/b^2+1)}$; explicitly for $n\ge\frac{24\Delta}{b^2}(\frac{4\eta}{b^2}+1)$, $C_1=\frac{18\log(1/b)}{b^4}\big(\frac{24\Delta}{b^2}\big)^{4\eta/b^2+1}$ and $T_{\mathrm{mix}}\le\lceil C_1n(\log n+\log\log\frac1b+\log\frac1{2\varepsilon^2})\rceil$. Consequences (their Thms 1.1–1.7): $O(n\log n)$ for antiferromagnetic 2-spin systems up to uniqueness, colourings of triangle-free graphs for $q\ge(\alpha^*+\delta)\Delta$, $\alpha^*\approx1.763$, matchings in $O(m\log n)$ (monomer–dimer is $\min\{2\lambda\Delta,2\sqrt{1+\lambda\Delta}\}$-SI, Thm 2.10).

*Proof architecture* ([CLV] §2). Four lemmas.
1. **Approximate tensorization** (Def. 2.1): $\mathrm{Ent}(f)\le C_1\sum_v\mu[\mathrm{Ent}_v(f)]$; **$\ell$-uniform block factorization** (Def. 2.2): $\frac\ell n\mathrm{Ent}(f)\le C\binom n\ell^{-1}\sum_{|S|=\ell}\mu[\mathrm{Ent}_S(f)]$.
2. Lemma 2.7: $\ell$-uniform block factorization with $C$ $\iff$ order-$(n-\ell,n)$ global entropy contraction with $\kappa$, $C\kappa=\ell/n$ (a direct consequence of the chain rule Lemma 5.1).
3. Lemma 2.5: $b$-marginally bounded and $\eta$-SI imply $\lceil\theta n\rceil$-uniform block factorization with $C=(2/\theta)^{4\eta/b^2+1}$, for $n\ge\frac2\theta(\frac{4\eta}{b^2}+1)$ (from Thm 5.4 + Thm 5.6, with $\hat\Gamma_k$ telescoping).
4. Lemma 2.3 (**shattering**): for $\theta\le\frac{b^2}{12\Delta}$, block factorization implies approximate tensorization with $C_1=\frac{18\log(1/b)}{b^4}C$. "The intuition … is that for $\ell$ as large as $\theta n$, if one picks a uniformly random subset $S$ … $G[S]$ … is disconnected into many small connected components … Since the conditional Gibbs distribution $\mu^\tau_S$ is a product distribution of each connected component, we can use entropy factorization for product distributions" (Lemma 4.1: $\mathrm{Ent}^\tau_S(f)\le\sum_{U\in\mathcal C(S)}\mu^\tau_S[\mathrm{Ent}_U(f)]$; Lemma 4.3: $\Pr(|S_v|=k)\le\frac\ell n(2e\Delta\theta)^{k-1}$).

Step 4 is the only place in this chain where the **Markov property** (conditional independence across a separator) and bounded degree are used (§4.2).

### 3.7 Entropic independence and fractional log-concavity

**Fact** (Anari–Jain–Koehler–Pham–Vuong [EI-I] Def. 2, Thm 4, Thm 5). $\mu$ on $\binom{[n]}k$ is **$(1/\alpha)$-entropically independent** if $D_{KL}(\nu D_{k\to1}\|\mu D_{k\to1})\le\frac1{\alpha k}D_{KL}(\nu\|\mu)$ for all $\nu$. Thm 4: with $p=\mu D_{k\to1}$, $(1/\alpha)$-EI $\iff$ $g_\mu(z_1^\alpha,\dots,z_n^\alpha)^{1/k\alpha}\le\sum_ip_iz_i$ on $\mathbb R^n_{\ge0}$ (the transformed generating polynomial lies below its tangent at $\vec1$); $\lambda*\mu$ is $(1/\alpha)$-EI for every external field $\lambda$ $\iff$ $\mu$ is **$\alpha$-fractionally log-concave** ($\log g_\mu(z^\alpha)$ concave). Thm 5: $\alpha$-FLC (or $(1/\alpha)$-EI of all links) gives $D_{KL}(\nu D_{k\to\ell}\|\mu D_{k\to\ell})\le(1-\kappa)D_{KL}(\nu\|\mu)$ for $\ell\le k-\lceil1/\alpha\rceil$, hence an MLSI $\Omega(k^{-1/\alpha})$ for the $k\leftrightarrow\ell$ walk. Picture of the three notions (their Fig. 1): SI = local convexity of a level set of $g_\mu(z^\alpha)$ at $\vec1$; EI = the level set lies above its tangent at $\vec1$; FLC = the level set is convex.

**Guard** ([EI-I] §1). Spectral independence alone does not give entropy contraction: the uniform law on the edges of a constant-degree expander graph has $2\leftrightarrow1$ walk = lazy random walk on the graph, with $\Omega(1)$ gap but mixing time $\simeq\log n$, while constant-factor entropy contraction would give $\log\log n$. This refutes Liu's conjecture that $O(1)$-SI implies MLSI $\Omega(1/k)$. The external-field hypothesis replaces the bounded-degree and Markov hypotheses of [CLV]: "Our framework makes no assumptions on marginals of $\mu$ or the degrees of the underlying graphical model".

**Fact** (Ising application, [EI-I] abstract and §1.3). Glauber on an Ising model whose interaction matrix has spectrum in an interval of length $<1$ mixes in $O(n\log n)$ (needle decomposition to rank one, then non-uniform FLC, then $D_{n\to1}$ entropy contraction, then induction).

**Fact** (restricted version, [EI-II] Def. 18, Prop. 19, Prop. 39, Thm 37). $(\eta,\epsilon)$-**spectral domination**: $\lambda_{\max}(\Psi^{cor}_{\lambda*\mu})\le\eta$ for all fields $\lambda\in(0,1+\epsilon]^n$; it implies FLC of $\mu^{\hom}$ on a cone and, for $\nu$ that are $C$-bounded w.r.t. $\mu$, $D_{KL}(\nu^{\hom}D_{n\to1}\|\mu^{\hom}D_{n\to1})\le\frac{\eta'}nD_{KL}(\nu^{\hom}\|\mu^{\hom})$, $\eta'=\max\{2\eta,\sqrt{\log C/\log(1+\epsilon)}\}$; consequence: sampling hardcore and Ising up to tree uniqueness in $\tilde O_\delta(n)$ time independent of $\Delta$.

### 3.8 Field dynamics, localization schemes, universality

**Fact** (Chen–Feng–Yin–Zhang [CFYZ] §2.1, Lemmas 2.3–2.4). The **field dynamics** $P^{FD}_\theta$ for $\mu$ on $\{\pm1\}^V$: put each $v$ with $\sigma_v=-1$ into $S$, each $v$ with $\sigma_v=+1$ into $S$ with probability $\theta$; resample $\sigma_S$ from $\pi^{\sigma_{V\setminus S}}_S$ with $\pi=\mu^{(\theta)}$, $\mu^{(\theta)}(\sigma)\propto\mu(\sigma)\theta^{\|\sigma\|_+}$. It is reversible w.r.t. $\mu$. Mixing lemma: complete $\eta$-SI gives $\lambda^{FD}_{gap}(\mu,\theta)\ge(\theta/2)^{2\eta+7}$. Comparison lemma: $\lambda^{GD}_{gap}(\mu)\ge\lambda^{FD}_{gap}(\mu,\theta)\cdot\lambda^{GD}_{\text{min-gap}}(\mu^{(\theta)})$, where the min-gap is over all pinnings. Result: optimal $\Omega(1/n)$ Glauber gap for antiferromagnetic 2-spin systems in the uniqueness regime, for *every* max degree.

**Fact** (Chen–Eldan [CE] Defs. 3, 5, 6; §2.2; Thm 24). A **localization process** is a measure-valued martingale $(\nu_t)$ with $\nu_t(A)\to\{0,1\}$; a **localization scheme** assigns one to each $\nu$; the attached chain is $P_{x\to y}=\mathbb E[\nu_\tau(x)\nu_\tau(y)/\nu(x)]$, reversible w.r.t. $\nu$. Coordinate-by-coordinate localization (pin coordinates in random order) with $\tau=n-1$ gives Glauber; with $\tau=n-\ell$ the $\ell$-block Glauber; subset-simplicial-complex localization gives the down-up walk; the **negative-fields localization** (increasing negative external fields with compensating pinnings to $+1$) gives the field dynamics; Gaussian channel (stochastic) localization gives the restricted Gaussian dynamics. Thm 24 (reformulation of [ALO] Thm 1.3): if $\rho(\Psi(R^u\nu))\le\eta_{|u|_1}$ for all pinnings $u$, the $k$-Glauber dynamics has gap $\ge\prod_{i=0}^{n-k-1}(1-\frac{\eta_i}{n-i})$; "In our proof, we completely bypassed the need to use the notion of high-dimensional expanders or the up-down walk." Thm 46/47 (annealing): concatenating two localizations, a variance (entropy) retention $\varepsilon$ in the first stage and a gap $\delta$ of the second stage's chain give gap $\varepsilon\delta$.

**Fact** (Anari–Koehler–Vuong [AKV], Thm 1, Thm 2, Eq. (10)). The trickle-down equation for stochastic localization with driving matrix $C_t$: $\mathrm{cov}(\nu_t)=\mathbb E[\mathrm{cov}(\nu_T)\mid\mathcal F_t]+\int_t^T\mathbb E[\mathrm{cov}(\nu_s)C_s^2\mathrm{cov}(\nu_s)\mid\mathcal F_t]ds$. Used to pass the convexity barrier $\beta=0.25$ for the SK model ($O(n\log n)$ Glauber up to $\beta\approx0.295$), and to sample antiferromagnetic Ising on $d$-regular **expander** graphs with $\max\{|\lambda_2(A)|,|\lambda_n(A)|\}\le\lambda$ for all $\beta<\frac1{2\lambda}$ in $\tilde O(nd)$ (polarized walk), i.e. up to $\Theta(1/\sqrt d)$ on random regular graphs instead of the worst-case $\Theta(1/d)$.

**Fact** (universality, [AJKPV-U] Thm 1 = Thm 30, Thm 32, Thm 43). If the $k\leftrightarrow k-1$ down-up walk of $\mu$ on $\binom{[n]}k$ has gap $\frac1{Ck}$, then $\mu$ is $C$-SI. Proof (four lines, their §3.1): apply the Poincaré inequality to $f(S)=\sum_{i\in S}v_i$, use the law of total variance over $S_{k-1}$, and $\mathrm{Var}\le\mathbb E[\cdot^2]$ on the single removed element: $\mathrm{Var}[\sum_{i\in S}v_i]\le C\,\mathbb E[\sum_{i\in S}v_i^2]$. So an $O(n)$ relaxation time of Glauber is *equivalent*, up to the losses of local-to-global, to $O(1)$-SI under all pinnings; and an $O(n)$ relaxation time under all external fields is equivalent to FLC. Thm 32: an inductive trickle-down for the $k\leftrightarrow k-1$ gap, $C''=C\frac{k-1-C}{k-2C}$. Thm 43: a contractive coupling ($\kappa\le1-\epsilon/n$, any metric) implies $1/\epsilon$-SI.

### 3.9 How spectral independence is certified

**Fact** (summarised from [AJKPV-U] §1, [ALO] §§1.6, 4–5, [CLV] §2.3). Known routes: (i) correlation decay on Weitz's self-avoiding-walk tree, decoupling the influence of a set into single-vertex influences, and the potential method on the tree recursion $R_r=\lambda\prod_u\frac1{1+R_u}$ ([ALO] Thm 1.13; the decay they obtain is "very similar to the notion of aggregate strong spatial mixing" of Mossel–Sly and Blanca–Chen–Vigoda); (ii) Dobrushin uniqueness / contractive coupling (Liu; Blanca–Caputo–Chen–Parisi–Štefankovič–Vigoda; [AJKPV-U] Thm 43); (iii) zero-freeness of the partition function (Alimohammadi–Anari–Shiragur–Vuong; Chen–Liu–Vigoda "stability"); (iv) universality from a known $O(n)$ relaxation time. [ALO] §4 adds: on **amenable** graphs (subexponential ball growth, e.g. $\mathbb Z^d$) strong spatial mixing gives $O(1)$-SI "in a black-box fashion", recovering the classical equivalence of strong spatial mixing and optimal mixing on lattices (Dyer–Sinclair–Vigoda–Weitz; Weitz), but this "exclude[s] most graphs such as expanders".

**Interpretation (two different expanders).** In this theory *the complex $X_\mu$ is the expander one wants*, while *expansion of the interaction graph $G$ is an obstacle* to certifying it: on a non-amenable $G$, pointwise correlation decay must beat the exponential growth of spheres, which is exactly the tree-uniqueness condition $|f'_d(\hat R_d)|<1$ up to degree $\Delta$. [AKV] Thm 2 shows the other side: spectral expansion of $G$ *helps* sampling antiferromagnetic Ising through a different chain. These are different objects and should not be conflated in the transfer (§6.5).

---

## 4. Where the Markov property sits, and what plays the "approximate local Markov property"

### 4.1 The exact statements

**Fact** (Hammersley–Clifford, in the form stated by Yang [Y] §I and Kato–Brandão [KB] §I.A; the 1971 manuscript itself was not re-read here). Positive Markov random fields on a graph are exactly the Gibbs distributions with clique-local interactions; on a 1D chain, positive Markov chains are exactly $p(x)=Z^{-1}\exp(-\sum_ih_i(x_i,x_{i+1}))$.

**Fact** ([KB] §I.A, citing Leifer–Poulin and Brown–Poulin). A full-rank quantum state on a chain is a quantum Markov chain ($I(A_1\dots A_{i-1}:A_{i+1}\dots A_n\mid A_i)=0$ for all $i$) iff it is $Z^{-1}\exp(-\sum_ih_{i,i+1})$ with **commuting** terms; for networks this equivalence holds only for triangle-free graphs ([KB] footnote 5).

**Fact** (CMI and recovery). $I(A:C\mid B)_\rho=S(\rho_{AB})+S(\rho_{BC})-S(\rho_B)-S(\rho_{ABC})$; $I=0$ iff exact recovery $\rho_{ABC}=(\mathcal R_{B\to AB}\otimes\mathrm{id}_C)(\rho_{BC})$ (Hayden–Jozsa–Petz–Winter, [Y] ref. 3); Fawzi–Renner: $I(A:C\mid B)\ge-2\log F(\rho_{ABC},\Delta_{B\to BC}(\rho_{AB}))$ for some channel ([KB] Eq. (3)).

**Fact** (Kato–Brandão [KB] Thm 1, Thm 4; "a variant of the Hammersley–Clifford theorem for 1D quantum systems"). An $\varepsilon$-approximate Markov chain on an open 1D chain of $m$ blocks has $S(\rho\|e^{-H}/Z)\le\varepsilon m$ for some $H=\sum_ih_{A_i,A_{i+1}}$ (proof: maximum-entropy state with the same nearest-neighbour marginals, iterated strong subadditivity, Pythagorean identity). Conversely, 1D Gibbs states of short-range $H$ admit a recovery map with error $e^{-q(\beta)\sqrt{d(A,C)}}$.

**Fact** (Chen–Rouzé [CR] Thm III.1, §IV.C, Remark III.1.1; Fig. 2). For $H$ of interaction degree $\le d$ and a region $A$, let $\mathcal R_{A,t}=\frac1t\int_0^te^{s\mathcal L_A}ds$, $\mathcal L_A=\sum_{a\in P^1_A}\mathcal L_a$, where $P^1_A=\{X_i,Y_i,Z_i\}_{i\in A}$ and $\mathcal L_a$ is the exactly detailed-balanced Lindbladian with jump $A_a$ and Metropolis weight. Then, with $\beta_0=1/4d$,
$$
\big\|\mathcal R_{A,t}[\mathrm{Tr}_A\rho_\beta\otimes\tau_A]-\rho_\beta\big\|_1\le|A|^22^{2|A|}\cdot\begin{cases}r(\beta,d)\,t^{-128\beta_0^4/(\beta^3(\beta+5\beta_0))}&\beta>4\beta_0\\ r'(\beta,d)\,t^{-2\beta_0/(\beta+5\beta_0)}&\beta\le4\beta_0\end{cases}\ \le\ r\,e^{\mu|A|}t^{-\lambda},
$$
$\tau_A$ maximally mixed. Time-averaging gives the gap-free bound $\mathcal E_A(\mathcal R^\dagger_{A,t}[X])\le\frac2t$ on the Dirichlet form (Cor. VII.1), and Lieb–Robinson makes the map quasi-local, radius $\sim\mathrm{Poly}(\beta)(|A|+\log\frac1\epsilon)$. Consequences: CMI decay $C_\beta|A||C|e^{C_\beta\min\{|A|,|C|\}-c_\beta r}$ ([Y] Table I), and quasi-local Gibbs preparation under uniform clustering. Their stated classical template: "(Exact Markov) + (Strong spatial mixing) ⟹ MCMC mixes in quasi-linear time" (Martinelli).

**Fact** (Yang [Y] Thm II.1, §I). For finite-range $H$ and a full partition $\Lambda=A\sqcup B\sqcup C$, $I(A:C\mid B)_{\rho_{\Lambda,\beta}}\le C_\beta e^{C_\beta g_{A|A^c}-c_\beta r}$, $g_{A|A^c}=\sum_{X\cap A\ne\emptyset\ne X\cap A^c}\|\Phi(X)\|_\infty\lesssim|\partial_eA|$: boundary instead of volume in the exponent. The proof is *static*: "we compare the Gibbs state with the reference state obtained by removing interactions across the cut $A\,|\,BC$. The reference state is a product across this cut, and the CMI is bounded by the decrease in relative entropy between the two states when $C$ is traced out."

**Derivation** (the inequality behind Yang's reduction). For any $\sigma=\sigma_A\otimes\sigma_{BC}$,
$$
D(\rho_{ABC}\|\sigma)-D(\rho_{AB}\|\sigma_{AB})=I(A:C\mid B)_\rho+\big[D(\rho_{BC}\|\sigma_{BC})-D(\rho_B\|\sigma_B)\big]\ \ge\ I(A:C\mid B)_\rho .
$$
*Proof.* Expand both relative entropies as $-S(\cdot)-\mathrm{Tr}(\cdot\log\sigma)$; the $\sigma_A$ terms cancel, the entropy terms give $S(\rho_{AB})-S(\rho_{ABC})$, and adding and subtracting $S(\rho_{BC})-S(\rho_B)$ gives the CMI plus the bracket, which is $\ge0$ by monotonicity under $\mathrm{Tr}_C$. $\square$ So the CMI is bounded by a **loss of relative entropy under a coarse-graining** — a sufficiency defect in the sense of Note 1 §6.

**Fact** (Bakshi–Liu–Moitra–Tang [BLMT] Def. 5.1, Lemma 5.2, Thm 2.1, Thm 2.3). A quantum Dobrushin influence matrix $D^{(\Phi)}_{ij}=\max\|\Phi_i(\rho-\sigma)\|_{W_1}$ over $j$-neighbours ($\mathrm{tr}_j(\rho-\sigma)=0$), in the quantum Wasserstein-1 norm of De Palma–Marvian–Trevisan–Lloyd; $\|D^{(\Phi)}\|_{1\to1}\le1-\gamma/n$ implies a unique fixed point and mixing in $\frac n\gamma\log\frac n\varepsilon$ steps (quantum path coupling). For a KMS-detailed-balanced Lindbladian with jumps $\sigma^{1/4}P\sigma^{-1/4}$, $P$ single-qubit Paulis: rapid mixing in $O(\log(n/\varepsilon))$ time for $\beta<1/(10^4K^3b^2d)$, and $I_\sigma(A:C\mid B)=O(a^a|A||C|)e^{-\mathrm{dist}(A,C)/\zeta}$ for $\beta<\beta_c$ (global Markov at high temperature, via a recovery map $e^{\mathcal L_{CMI}t}$ supported near $C$).

### 4.2 The pipeline, ingredient by ingredient

| ingredient of the classical HDX/SI theory | uses the Markov property? | its role |
|---|---|---|
| encoding $X_\mu$; links = pinnings ([ALO] §1.1) | no; defined for every $\mu$ | the complex of conditionals |
| SI ⇔ local spectral expansion ([ALO] Thm 1.11, 3.1) | no | an identity for any law on $\{0,1\}^n$ |
| local-to-global ([AL]; [KO]; [CLV] Thm 5.4; [EI-I] Thm 5; [CE] Thm 24) | no | chain rule along the pinning filtration |
| heredity: every link of a spin system on $G$ is a spin system on $G[V\setminus\Lambda]$ with modified fields | **yes** (conditional law depends only on the boundary) | lets one argument certify SI at every link |
| certifying SI by the SAW tree and the tree recursion ([ALO] Thm 1.13) | **yes** (Weitz's tree; conditional independence of subtrees given the root) | correlation decay ⇒ bounded total influence |
| shattering, block factorization ⇒ approximate tensorization ([CLV] Lemma 2.3, 4.1, 4.3) | **yes**, plus bounded degree | optimal $O(n\log n)$ |
| entropic independence under external fields ([EI-I]) | no ("no assumptions on … the underlying graphical model") | replaces shattering by a hypothesis on all tilts |

### 4.3 The identification

**Interpretation.** Three statements, each anchored in the facts above.

(a) **Exact Markov property = a structural property of the links.** For a Markov random field on $G$, the link of $X_\mu$ at a face $(U,\tau)$, restricted to a region $A\subseteq V\setminus U$, depends on $\tau$ only through $\tau|_{\partial A}$ when $U\supseteq\partial A$; and if $U$ separates $V\setminus U$ into components $V_1,\dots,V_m$, the link at $(U,\tau)$ is the **join** $X_{\mu^\tau|V_1}*\cdots*X_{\mu^\tau|V_m}$ (product law ⇒ join of complexes), whose influence matrix is block-diagonal. This is what [CLV] Lemma 4.1 exploits and what makes the class hereditary.

(b) **Spectral independence = aggregate strong spatial mixing = uniform clustering, not Markov.** $\Psi_{\mu^\tau}(u,v)$ compares the links at $\tau\cup\{u\}$ and $\tau\cup\{\bar u\}$ on the *unpinned* site $v$: the influence travels through unpinned territory. The Markov property says nothing about it (it is the trivial statement that the influence vanishes when everything in between is pinned).

(c) **Martinelli's template is the HDX template.** "(Exact Markov) + (SSM) ⇒ quasi-linear mixing" reads, in HDX language, "(links are local and split as joins at separators) + (every link is a uniform local spectral expander) ⇒ global gap / MLSI". The quantum programme of Chen–Rouzé, Yang and Bakshi et al. keeps this template and has to *prove* the first half, because for non-commuting $H$ there is no pinning: conditioning a quantum state on a region is not an operation on states. **The approximate local Markov property is the quantum surrogate of ingredient (a)**: it asserts that the "link" of the Gibbs state at the complement of $A$, realised as a recovery channel, exists and is quasi-local around $A$, with error $e^{c|A|-r/\xi}$ (Chen–Rouzé) or $e^{c|\partial A|-r/\xi}$ (Yang). Uniform clustering (Chen–Rouzé's other hypothesis) is the surrogate of (b), and the quantum Dobrushin condition of [BLMT] is the quantum surrogate of route (ii) of §3.9 — a condition that, classically, *implies* spectral independence ([AJKPV-U] Thm 43).

### 4.4 The time-averaged single-Pauli-jump Lindbladian as a quantum link-resampling map

**Derivation** (classical counterpart; checked in §9, C4). Let $\mu$ be a positive law on $\Omega^V$, $A\subseteq V$, and $\mathcal L_A=\sum_{a\in A}(E_a-I)$ the heat-bath generator on $A$ ($E_a$ resamples site $a$ from $\mu(\cdot\mid\sigma_{V\setminus a})$). Assume the $A$-restricted chain is irreducible for every boundary condition. Then:
1. $\mathcal L_A$ is self-adjoint in $L^2(\mu)$, its fixed points are the functions of $\sigma_{A^c}$, and $e^{t\mathcal L_A}\to E_A:=\mu(\cdot\mid\sigma_{A^c})$; by the mean ergodic theorem also $\mathcal R_{A,t}:=\frac1t\int_0^te^{s\mathcal L_A}ds\to E_A$.
2. **Exact recovery.** For $\nu=\mu_{A^c}\otimes\mathrm{unif}_A$, $\nu E_A=\mu$: resampling $A$ from its conditional erases whatever was there.
3. **Locality ⇔ Markov.** $E_A$ acts on $\sigma_A$ through $\mu(\sigma_A\mid\sigma_{A^c})$; for a Markov random field this equals $\mu(\sigma_A\mid\sigma_{\partial A})$, so the recovery map is supported on $A\cup\partial A$ with zero error.
4. **Gap-free Dirichlet bound.** By spectral calculus, with $\phi_t(x)=\frac{e^{tx}-1}{tx}$, $\mathcal E(\mathcal R_{A,t}f)=\langle\phi_t(\mathcal L_A)f,-\mathcal L_A\phi_t(\mathcal L_A)f\rangle\le\sup_{y\ge0}\frac{(1-e^{-ty})^2}{t^2y}\|f\|^2\le\frac{0.41}t\|f\|^2$.
5. **HDX reading.** $E_A$ is the composite "restrict to the face $\sigma_{A^c}$, co-restrict by sampling from the link": one step of the order-$(n,n-|A|)$ down-up walk with the removed block fixed to be $A$. Averaging $E_A$ over uniformly random $A$ of size $\ell$ gives $P^\vee_{n,n-\ell}$, the block dynamics of [CLV].

**Interpretation.** Chen–Rouzé's $\mathcal R_{A,t}$ is item 1–4 for a non-commuting Gibbs state: the detailed-balanced generator with single-site Pauli jumps on $A$ is the quantum heat-bath on $A$ (the single-qubit Paulis generate all operators on a site, so its fixed points are what does not "act on $A$", dressed by the state); its long-time or time-averaged limit is the quantum $E_A$, i.e. **the quantum link of the face $A^c$**, and their Thm III.1 is item 2 with an error; item 3 fails exactly (the dressing $\sigma^{\pm1/4}P\sigma^{\mp1/4}$ is only quasi-local), which is why a shield of width $r$ replaces the one-layer separator $\partial A$; item 4 is their Cor. VII.1. Two further points of contact:

- **Fact** (Bardet–Capel–Rouzé [BCR] Thm 1, §2.4). The conditional expectation onto the fixed-point algebra of a Davies generator coincides with the one generated by a Petz recovery map: the quantum "$E_A$" is simultaneously a thermal dynamics limit and a recovery map. Chen–Rouzé's map is a time-averaged version for non-commuting $H$.
- **Interpretation.** Chen–Rouzé remark that a polynomial **local gap** of $\mathcal L_A$ would improve their bound. Classically the gap of $\mathcal L_A$, minimised over boundary conditions, is $\lambda^{GD}_{\text{min-gap}}$ of [CFYZ] — the Glauber gap *inside the links* — and inside each link it is controlled by local-to-global from the spectral independence of that link and its pinnings. So the quantum "local gap condition" is the counterpart of "local spectral expansion of all links", one level up.

**Interpretation (what is *not* deep).** The time average is a device: $\frac1t\int_0^te^{s\mathcal L}ds$ converges to the fixed-point projection with a $1/t$ Dirichlet bound and no spectral gap. Nothing in the HDX theory corresponds to it beyond the Cesàro mean of a semigroup; the structural content is in the *limit object* $E_A$ and in its locality.

### 4.5 Operator-algebraic form: commuting squares, conditional expectations, sufficiency

**Derivation** (classical; checked in §9, C3). Let $\mathcal F_{\le\ell},\mathcal F_{\ge\ell},\mathcal F_\ell$ be the σ-algebras of a process up to, from, and at time $\ell$ (or, for a field, of $A^c$, $B^c$ and $(A\cup B)^c$). Then
$$
\mathbb E[\,\cdot\mid\mathcal F_{\le\ell}]\ \mathbb E[\,\cdot\mid\mathcal F_{\ge\ell}]=\mathbb E[\,\cdot\mid\mathcal F_\ell]
\iff \text{past}\perp\text{future}\mid\text{present}.
$$
*Proof.* For $g$ measurable w.r.t. $\mathcal F_{\ge\ell}$, $\mathbb E[g\mid\mathcal F_{\le\ell}]=\mathbb E[g\mid\mathcal F_\ell]$ for all bounded $g$ is the definition of the conditional independence; applying both sides to $f$ and using the tower property gives the operator identity, and conversely. $\square$ For a field: $E_AE_B=E_{A\cup B}$ iff $\sigma_{A\setminus B}\perp\sigma_{B\setminus A}\mid\sigma_{(A\cup B)^c}$; note that $A\cap B$ is *unpinned*, so this is a clustering statement across the buffer $A\cap B$, not a Markov statement.

**Fact** ([BCR] Eq. (1.2), (1.3), Def. 2, Thm 2, Thm 3). For von Neumann algebras $\mathcal M\subset\mathcal N_1,\mathcal N_2\subset\mathcal N$ with conditional expectations forming a **commuting square** ($E_1E_2=E_2E_1=E_{\mathcal M}$): $D(\rho\|E_{\mathcal M*}\rho)\le D(\rho\|E_{1*}\rho)+D(\rho\|E_{2*}\rho)$, of which strong subadditivity is the case of partial traces. Without commuting squares, **approximate tensorization** $\mathrm{AT}(c,d)$: $D(\rho\|E_{\mathcal M*}\rho)\le c(D(\rho\|E_{1*}\rho)+D(\rho\|E_{2*}\rho))+d$, with $c=\frac1{1-c_1}$, $c_1=\max_i\|E^{(i)}_1E^{(i)}_2-E^{(i)}_{\mathcal M}:L^1(\tau_i)\to L^\infty\|<1$; the additive $d$ "is exclusively due to the non-commutativity of the underlying algebras" (it vanishes for classical Hamiltonians embedded in quantum systems, Thm 3). Classically (Cesi; Dai Pra–Paganoni–Posta) this "quasi-factorization" is "the main ingredient in modern proofs of modified logarithmic Sobolev inequalities" for lattice spin systems.

**Interpretation.**
- The HDX local-to-global theorems are chain rules along a tower of *commuting* conditional expectations (pinning is commutative); [CLV]'s block factorization is the many-block version of $\mathrm{AT}(c,0)$. [BCR] is the two-step non-commutative version. **No $k$-step "quantum trickle-down" or quantum local-to-global theorem with link-level hypotheses was found in this session's searches** (alphaXiv discovery, Exa; §8); the nearest objects are [BCR], the double-localization gap assembly of Hayakawa et al. ([HSL] Thm 1: interface gap $a$, local gaps $g_{loc}$ on overlapping regions, cover/overlap/Schur conditions ⇒ $K\succeq(1-\sqrt{1-\eta})\min\{a,\tau\}I$), and [BLMT]'s quantum Dobrushin.
- Exact Markov at a level ⇔ a commuting square; approximate Markov ⇔ either a small commuting-square defect (sup-norm, $c_1$ of [BCR]) or a small sufficiency defect (averaged, CMI; Petz/Fawzi–Renner). These are the operator-algebra statements underlying both sides of the bridge, and they are the same family as Note 1 §6: the Jenčová–Petz theorem characterises sufficiency of a subalgebra by invariance of relative entropy under restriction, which is the zero-defect case of the Derivation in §4.1.

### 4.6 Verdict on the three intuitions (tested, not assumed)

| intuition | verdict | grounds |
|---|---|---|
| Gibbs states / MRFs share an essence with what drives expander theory | **supported at theorem level** | [ALO] Thm 1.5/1.11: SI of a law = local spectral expansion of its complex of conditionals; Glauber = down-up walk; [AJKPV-U] Thm 30: optimal Glauber relaxation ⇒ SI. The shared essence is the local-to-global chain rule along conditioning (§2.8), and "expander" refers to the links of $X_\mu$, not to the interaction graph (§3.9) |
| the time-averaged detailed-balanced single-Pauli-jump Lindbladian on $A$ hides a deep analogue | **supported, with a precise dictionary** | it is the quantum version of the link-resampling map $E_A$ (restriction ∘ co-restriction), i.e. of one block step of the down-up walk; what it adds is a quantitative locality of the quantum link (§4.4). The time averaging itself is a technical device |
| noncommutative geometry ties these together | **partially supported** | the ties found are operator-algebraic: conditional expectations, commuting squares, Petz recovery, sufficiency, KMS detailed balance (§4.5). No spectral triple, $K$-theory or cyclic cohomology appears in any of the bridges read. The missing theorem is a non-commutative local-to-global (tower) inequality with link-level hypotheses |

---

## 5. Known theorems connecting the areas

Each row is a theorem (or an equivalence) linking two of: Markov random fields / Gibbs measures (MRF), expanders and HDX, local-to-global (LtG), Markov semigroups / functional inequalities (MS), noncommutative / operator algebras (NC). "Read" = read in this session.

| # | link | statement | source | status |
|---|---|---|---|---|
| 1 | MRF ↔ HDX | links of $X_\mu$ = pinnings; $\lambda_2(P_\emptyset)=\lambda_{\max}(\Psi_\mu)/(n-1)$; SI ⇔ local spectral expansion | [ALO] Thm 1.5, 1.11, 3.1 | read |
| 2 | HDX ↔ MS | Glauber = top down-up walk of $X_\mu$; block dynamics = $P^\vee_{n,n-\ell}$ | [ALO] §1.1; [CLV] §1.2 | read |
| 3 | LtG | Garland: link gaps $>\frac k{k+1}$ ⇒ $\tilde H^k=0$ and global Laplacian bounds | [O14] Lemma 5.4, Cor. 5.6, Thm 5.8 | read |
| 4 | LtG | trickle-down $\gamma\mapsto\gamma/(1-\gamma)$ | [O14] Lemma 5.1; [AL] Thm 1.3; [KO] Lemma 4.5 | read |
| 5 | LtG ↔ MS | product formula for $\lambda_2$ of down-up walks | [AL] Thm 1.5, 3.2 | read |
| 6 | LtG ↔ MS | local entropy contraction ⇒ global, $\kappa=\sum_{i\ge k}\Gamma_i/\sum_i\Gamma_i$; marginals + local spectral expansion ⇒ local entropy contraction | [CLV] Thm 5.4, 5.6 | read |
| 7 | MRF ↔ MS | SI + marginal bounds + bounded degree ⇒ MLSI $\Omega(1/n)$ | [CLV] Thm 1.12 | read |
| 8 | MRF ↔ MS | SI under all external fields ⇔ FLC ⇒ entropic independence ⇒ MLSI | [EI-I] Thm 4, 5 | read |
| 9 | MS ↔ MRF | Glauber gap $\frac1{Cn}$ ⇒ $C$-SI (universality) | [AJKPV-U] Thm 1, 30 | read |
| 10 | MS ↔ MRF | contractive coupling / Dobrushin ⇒ SI | [AJKPV-U] Thm 43 (also Liu; Blanca et al.) | read (statement) |
| 11 | MRF ↔ MS | correlation decay (tree uniqueness) ⇒ SI; on amenable graphs SSM ⇒ SI; SSM ⇔ optimal mixing on lattices | [ALO] Thm 1.13, §1.3, §4 | read |
| 12 | LtG ↔ MS | localization schemes: gap $=\inf\mathbb E\mathrm{Var}_{\nu_\tau}/\mathrm{Var}$; SI = approximate variance conservation; annealing | [CE] Prop. 19, Thm 24, 46 | read |
| 13 | LtG | trickle-down = law of total covariance in linear-tilt localizations | [AKV] Eq. (8)–(12), Thm 101 | read |
| 14 | HDX ↔ (testing / sheaves) | HD expanders ⇒ agreement expanders, linear-size double samplers; "approximate cohomology with local coefficients" | [DK] Thm 1.6–1.9, §1.6 | read |
| 15 | HDX ↔ MS | skeleton expansion ⇒ colorful expansion ⇒ rapid mixing of all high-order walks | [KM] Thm 2.1, 3.1, 3.2 | read |
| 16 | HDX | coboundaries are the obstruction; decomposition of cochains by dimension | [KO] Thm 1.3–1.5 | read |
| 17 | MRF ↔ NC | quantum Hammersley–Clifford: full-rank quantum Markov chains = commuting 1D Gibbs states; approximate Markov ⇒ $S(\rho\Vert e^{-H}/Z)\le\varepsilon m$ | [KB] §I.A, Thm 1 | read |
| 18 | MRF ↔ NC ↔ MS | quantum Gibbs states are locally Markov at all $\beta$; recovery map = time-averaged detailed-balanced Lindbladian with single-Pauli jumps on $A$ | [CR] Thm III.1 | read |
| 19 | MRF ↔ NC | CMI $\le C_\beta e^{C_\beta\lvert\partial_eA\rvert-c_\beta r}$, static proof by a cut Gibbs state and relative-entropy loss | [Y] Thm II.1 | read |
| 20 | MS ↔ NC | quantum Dobrushin (Wasserstein influence matrix) ⇒ rapid mixing and global CMI decay at high temperature | [BLMT] Lemma 5.2, Thm 2.1, 2.3 | read |
| 21 | MS ↔ NC | commuting squares ⇒ SSA-type tensorization; approximate tensorization for non-commuting conditional expectations; Davies conditional expectation = Petz conditional expectation | [BCR] Thm 1, 2, 3 | read |
| 22 | MS ↔ NC | abstract gap assembly from an interface gap and local gaps (double localization) | [HSL] Thm 1 | read (statement) |
| 23 | NC ↔ MRF | CMI $=0$ ⇔ exact (Petz) recovery; Fawzi–Renner approximate recovery | HJPW; Fawzi–Renner, as cited in [Y], [KB] | cited, not read |
| 24 | NC (Note 1) | sufficiency of a subalgebra ⇔ invariance of relative entropies ⇔ Connes cocycles in the subalgebra | Jenčová–Petz, Thm 1 (read for Note 1) | read earlier |

---

## 6. Relation to the prong-3 note

### 6.1 Corrections and additions to [local-to-global-unlocks.md](../../local-to-global-unlocks.md)

- **Identifier.** The note cites "Oppenheim [1407.8517] … *Part I: Descent of spectral gaps*, Discrete Comput. Geom. 2018". arXiv:1407.8517 (v4) is titled *Local spectral expansion approach to high dimensional expanders* and contains the descent lemma as Lemma 5.1 and Garland's method as Lemma 5.4 / Cor. 5.6; the paper titled *… part I: Descent of spectral gaps* is arXiv:1709.04431 (alphaXiv listing). The statements quoted in the note are correct; the identifier for "Part I" should be 1709.04431 (programme rule P3.3).
- **§2.4, ALO.** The note's last line asked for the influence-matrix definition to be re-read before U6 is made precise; it is now §3.3, and the theorem (§3.4) is an identity, valid for any law, with no Gibbs or graph hypothesis.
- **§2.5, Dinur–Kaufman.** Add Thm 1.7–1.9 (one-up and $t$-up walks, double samplers) and the proof form actually given (Thm 5.1: $k^2<d$, $\lambda<1/d$, majority decoding); the method of decreasing differences is the variance chain rule of §2.8 in a two-sided form.
- **§4, row "transport operator".** The four identities of §2.8 make the row precise: the transport operator is a level of a tower of conditional expectations, and the "invariant controlled by local data" is the per-level retention of variance or entropy.
- **§4, "the obstruction is always a coboundary".** Confirmed by [KO]'s abstract and decomposition theorems; add that in the probabilistic setting the slow directions of the down-up walk are the functions *lifted from lower levels of the pinning tower* ($V^j_k=d_{j\nearrow k}C^j_0$), which for $X_\mu$ are functions of few coordinates.

### 6.2 U6 is a theorem once the dictionary is fixed

**Derivation** (pure translation; no new mathematics). Note 2 §1.1: a state $\omega$ on the face algebra $D_0(K)\cong C(K)$ is a probability measure $\mu_\omega$ on the faces of $K$, i.e. on a downward-closed subset of $2^V$; the vertex projections $e_u$ are the indicators "$u\in$ face". Lüders conditioning by $e_u$ or $1-e_u$ (Note 1 §4.2) is pinning $u$ in or out:
$$
\omega^{e_u}(a)=\frac{\omega(e_uae_u)}{\omega(e_u)},\qquad
\Psi_\omega(u,v)=\omega^{e_u}(e_v)-\omega^{1-e_u}(e_v),
$$
which is ALO's $\Psi_{\mu_\omega}$. Hence, with $n=|V|$: if $\omega$ and all its iterated Lüders conditionings on vertex projections satisfy $\lambda_{\max}(\Psi)\le\eta_i$ ($i$ = number of conditionings), then the walk "pick a vertex $u$ uniformly; re-decide $e_u$ by Lüders-conditioning $\omega$ on the current status of all other vertices" has spectral gap $\ge\frac1n\prod_{i=0}^{n-2}(1-\frac{\eta_i}{n-i-1})$ ([ALO] Thm 1.3), and the complex $X_{\mu_\omega}$ (on $\{u,\bar u\}_{u\in V}$) is a $(\frac{\eta_i}{n-i-1})_i$-local spectral expander ([ALO] Thm 1.5). Irreducibility holds when the support is downward closed (every face connects to $\emptyset$). **So U6's implication is a Fact; the only conjectural content left in U6 is whether the states that arise satisfy SI**, which is a question for prong 2 to measure and prong 1 to explain.

Two complexes must be kept apart: $K$ itself (faces = included-vertex patterns; links = "in"-pinnings only; relevant when $\omega$ is supported on the facets of a pure $K$, the homogeneous case of U2, where $C$-SI ⇔ $\lambda_2(U_{1\to k}D_{k\to1})\le C/k$, [AKV] Lemma 98) and $X_{\mu_\omega}$ (faces = partial in/out decisions; $|V|$-partite; only one-sided expansion possible, §1.2 Guard).

### 6.3 U2 in the light of §3

The negative-correlation case is explicit in [ALO] §1: a $d$-homogeneous law all of whose conditionals are negatively correlated is $(1,\dots,1)$-SI. Strong log-concavity (matroids) gives $0$-local spectral expansion of all links; fractional log-concavity is the graded relaxation that still gives MLSI ([EI-I] Thm 5). For U2 the hereditary property the prong-3 note asks for ("an exchange property for faces") is now replaceable by a weaker, checkable one: **$\alpha$-FLC of the facet-generating polynomial**, equivalently SI of all tilts, equivalently ([AJKPV-U] §1.1) $O(k)$ relaxation of the down-up walk under all external fields.

### 6.4 Candidate unlocks added by this digest

**U7 (memory = distance to the Gibbs/KMS family).** *Derivation* (classical, checked in §9, C2). For any law $P$ on histories $(\sigma_0,\dots,\sigma_L)$ and its Markov projection $Q(\sigma)=P(\sigma_0)\prod_iP(\sigma_{i+1}\mid\sigma_i)$,
$$
D(P\|Q)=\sum_{i=1}^{L-1}I(\sigma_{i+1};\sigma_{<i}\mid\sigma_i),\qquad
\min_{R\ \text{Markov}}D(P\|R)=D(P\|Q),
$$
the second by the Pythagorean identity for the log-linear family of Markov chains (whose sufficient statistics are the consecutive pair marginals, which $P$ and $Q$ share). By Note 1 §5.2–5.3 the Gibbs/KMS states on complete histories, with free cocycle $F$ and origin weights, are exactly the Markov chains on the quiver (Doob transform). Hence **the "memory" measured by prong 2 (E2, with the whole past) sums to the exact relative-entropy distance of the realized history law from the KMS family of Note 1**, and lag-1 values are lower bounds on the per-layer terms. This is the classical shadow of Kato–Brandão Thm 1 (approximate Markov ⇒ close to Gibbs, error $\varepsilon m$), and it gives prong 2's 25–35 % figure a precise meaning in prong 1's vocabulary. *What must be true:* nothing beyond the definitions; it is an identity. *What it unlocks:* a principled test for dictionary v2 (D2 of the MLP bridge): the minimal sufficient sub-frame is the one that drives this sum to zero.

**U8 (Markov at a layer = commuting square in the history algebra).** *Derivation* (§4.5) + *Conjecture.* In Note 1's tower $D_0\subset D\subset C^*(\Lambda)$, with the conditional expectations of the history algebra $D$ onto "functions of the history up to layer $\ell$", "from layer $\ell$ on" and "of the face at layer $\ell$" (all inside the commutative $D$), the commuting-square identity holds iff the law on histories is Markov at $\ell$ (proved). *Conjecture:* for a state on the full $C^*(\Lambda)$ with coherences across re-routings (Note 1 §4.1; not a KMS state, whose restriction to $D$ determines it), the right non-commutative statement is an approximate commuting square in the sense of [BCR] for state-preserving conditional expectations onto the past, future and present subalgebras, with vanishing defect equivalent to a quantum Markov-chain structure along the layers, hence, by the quantum Hammersley–Clifford theorem of §4.1, to a Gibbs form with commuting layer-to-layer terms. *Caution:* this is a statement about **one** state. Sufficiency of the face for the whole Gibbs family (Note 1 §6, "$F$ is a coboundary") is a different statement: every KMS state of Note 1 is Markov, whatever $F$ is, so its classical commuting square holds for every $F$. *What must be true:* that state-preserving conditional expectations onto the face, past and future subalgebras exist (for the commutative $D$ they always do; in $C^*(\Lambda)$ this needs the subalgebras to be invariant under the state's modular flow (Takesaki's criterion, cited from memory); for KMS states Note 1 §5.2 identifies that flow with weighted depth, a gauge action under which the face, past and future algebras are invariant, but for a state with coherences it has to be checked) and compose like [BCR]'s $E_1,E_2,E_{\mathcal M}$; this is a prong-1 question about the relative position of $D_0$, the past algebra and the future algebra.

**U9 (local-to-global does not need the Markov property).** *Interpretation with a theorem behind it.* Every local-to-global theorem in §2 and [ALO] Thm 1.3, [CE] Thm 24, [EI-I] Thm 5 hold for arbitrary laws. Therefore the non-Markov history law of prong 2 is no obstruction to the HDX route: the relevant complex is $X_\mu$ for $\mu$ the law of the face process (or of vertex activity, §6.2), and the relevant certificate is SI of its pinnings. What *is* lost without a Markov property is (i) heredity (pinnings need not be "of the same kind") and (ii) shattering (no exact block factorization); [EI-I] replaces both by a hypothesis on all external fields, which in the framework's language are **tilts of the state by elements of the face algebra**, $\omega\mapsto\omega(e^{h}\,\cdot\,)/\omega(e^h)$ with $h=\sum_u h_ue_u\in D_0$ (Interpretation).

### 6.5 Guards specific to this bridge

- **Pinning is commutative.** Links of $X_\mu$ exist because conditioning on a commutative frame is a state-level operation. For a non-commutative algebra the replacement is a conditional expectation onto a fixed-point or Petz subalgebra, and towers of those do not commute ([BCR]'s additive $d$). Any transfer of §2 beyond the commutative frame $C(V)$ needs a non-commutative chain rule, which is not available in the literature found (§4.5).
- **Partite ⇒ one-sided.** $X_\mu$ is $n$-partite; only one-sided theorems ([KO] Thm 1.4, [AL], [CLV], [EI-I]) apply. This matches the prong-3 guard on layered structures.
- **SI ≠ MLSI.** Constant SI gives only an $n^{-(1+\eta)}$-type gap through [ALO]/[AL], and no entropy contraction without extra hypotheses (marginals + bounded degree + Markov for [CLV]; all external fields for [EI-I]); the expander-edge example of [EI-I] shows the extra hypothesis is necessary.
- **Two expanders.** Expansion of the complex of conditionals is the goal; expansion of the interaction graph is, for worst-case correlation decay, an obstacle (§3.9). A statement in the framework that "the quiver is an expander" is not a statement about $X_\mu$.
- **Constants.** The quantum Markov bounds carry $e^{\mu|A|}$ ([CR]) or $e^{C|\partial A|}$ ([Y]) prefactors and polynomial-in-$t$ rates; classical exact Markov has none. Any use for regions $A$ comparable to the system size is outside their scope.
- **Drag test (rule C2).** Everything in §§3–4 is stated for a law on a product space or a state on a commutative frame, or for von Neumann algebras with conditional expectations; it would be stated identically for a tiling's patch frequencies or a Bratteli diagram's path measure, so it passes the drag test. The specific spin models (hardcore, Ising) are examples, not definitions.

---

## 7. Messages to the other prongs

**To prong 1 (theory).** (i) Define, for a state $\omega$ on $D_0(K)$, the influence matrix $\Psi_\omega$ by Lüders conditioning (§6.2) and record ALO's identity as a Fact of the framework; this turns U6 into a statement about states. (ii) Compute the commuting-square defect of the triple (past algebra, future algebra, face algebra at layer $\ell$) inside $D$ and inside $C^*(\Lambda)$ for KMS states; prove or refute that it vanishes iff $F$ is a coboundary (U8). (iii) External fields on the frame are tilts by $e^h$, $h\in D_0$; FLC of a state is then "SI of all tilts", a condition on the face algebra alone (U9).

**To prong 2 (bridge).** (i) Report memory with the whole past, layer by layer; their sum is the exact KL distance to the KMS family (U7), a calibrated quantity (the Markov surrogate gives $0$). (ii) Measure $\lambda_{\max}(\Psi)$ for the law of unit activity (or of faces) and for a sample of its pinnings; estimate $\eta_i$ against $i$. A growth $\eta_i\approx c$ with $c$ independent of $|V|$ is the HDX certificate; growth with $|V|$ is the obstruction. Calibrate with a product surrogate ($\Psi=0$) and a two-cluster surrogate ($\lambda_{\max}=n-1$). (iii) For the trivial-nerve finding (E6), note that HDX certificates are about the *weights* on the complex (links as conditional laws), not about the homology of the support: a nearly complete complex can still have large $\eta$.

**To prong 3 (unlocks).** U6 is a Fact modulo SI; U7 is an identity; U8 and U9 are new. The NCG-specific open problem distilled from this bridge: **a non-commutative local-to-global inequality** — a tower of conditional expectations $E_0\leftarrow E_1\leftarrow\cdots\leftarrow E_k$ on a finite-dimensional von Neumann algebra, with per-level hypotheses on "links" (relative positions of consecutive levels), implying contraction of $D(\rho\|E_{0*}\rho)$ — whose two-level case is [BCR] Thm 2 and whose commutative case is [CLV] Thm 5.4.

---

## 8. Sources read in this session

- [ALO] N. Anari, K. Liu, S. Oveis Gharan, *Spectral independence in high-dimensional expanders and applications to the hardcore model*, [2001.00303](https://arxiv.org/abs/2001.00303) (FOCS 2020): full text §§1–3 (Defs. 1.1–1.4, Thms 1.3, 1.5, 1.6, 1.8, 1.11, 1.13, 3.1, Claims 3.2–3.3, Remark 3.4, Lemma 1.12, Remarks 1.7, 1.10, 1.14, §1.4, §2), §4 opening.
- [CLV] Z. Chen, K. Liu, E. Vigoda, *Optimal mixing of Glauber dynamics: entropy factorization via high-dimensional expansion*, [2011.02075](https://arxiv.org/abs/2011.02075) (STOC 2021): Defs. 1.10, 1.11, 1.15, 1.17, 2.1, 2.2, 2.6, 5.3, A.5, A.7; Thms 1.12, 1.19, 2.9, 2.10, 5.4, 5.6, A.9; Lemmas 2.3, 2.5, 2.7, 2.8, 4.1–4.4, 5.1; Claims 1.16, 1.18; Facts 5.2, A.6, A.8; §7.
- [AL] V. L. Alev, L. C. Lau, *Improved analysis of higher order random walks and applications*, [2001.02827](https://arxiv.org/abs/2001.02827): §§1–3 (Def. 1.1, Thms 1.2, 1.3, 1.5, 2.5, 3.1, 3.2, Cors. 1.4, 1.6, 1.11, 1.12, 2.6, 3.4, 3.5, Prop. 3.3, Lemmas 3.6, 3.7, §2.2–2.4).
- [KO] T. Kaufman, I. Oppenheim, *High order random walks: beyond spectral gap*, [1707.02799](https://arxiv.org/abs/1707.02799): abstract, §1, §2, §4.2 (Lemma 4.5, Cors. 4.6–4.7), §5 (Thms 5.2, 5.4, 5.6, 5.9, 5.10, Cors. 5.3, 5.8, 5.11).
- [O14] I. Oppenheim, *Local spectral expansion approach to high dimensional expanders*, [1407.8517](https://arxiv.org/abs/1407.8517): Def. 1.1, §2, §4, §5 (Lemmas 5.1, 5.4, 5.15, Cors. 5.6, 5.17, Thm 5.8, Prop. 5.11), §7 (Thm 7.12, Cor. 7.13), §8 (Thm 8.9, 8.12, Cors. 8.10, 8.13). (The *Part I* paper is [1709.04431](https://arxiv.org/abs/1709.04431), listed but not opened.)
- [KM] T. Kaufman, D. Mass, *High dimensional random walks and colorful expansion*, [1604.02947](https://arxiv.org/abs/1604.02947): §1, Thms 2.1, 3.1–3.3, 4.1, Lemma 4.5, §1.4.
- [DK] I. Dinur, T. Kaufman, *High dimensional expanders imply agreement expanders*, ECCC [TR17-089](https://eccc.weizmann.ac.il/report/2017/089/) (FOCS 2017): §1 (Defs. 1.1, 1.2, 1.4, Thms 1.3, 1.6–1.9, Lemma 1.5, §1.6), §2.3, §3 (Lemma 3.1), §5 (Thm 5.1, Claims 5.2–5.3), §7.
- [EI-I] N. Anari, V. Jain, F. Koehler, H. T. Pham, T.-D. Vuong, *Entropic independence I*, [2106.04105](https://arxiv.org/abs/2106.04105): §1 (Defs. 1–3, Thms 4, 5, Fig. 1), §1.3, §3 (proof of Thm 4), Prop. 31, App. A (Example 38), App. B.
- [EI-II] Same authors, *Entropic independence II*, [2111.03247](https://arxiv.org/abs/2111.03247): abstract, §1, Defs. 16–18, 22, 24, 27–30, Prop. 19, 39, Thm 37, 43, Cor. 42.
- [CE] Y. Chen, R. Eldan, *Localization schemes*, [2203.04163](https://arxiv.org/abs/2203.04163): §§1.1, 2 (Defs. 3, 5, 6, 10, 11, Prop. 18, 19, examples), 3.1 (Claim 22, Eq. (13)–(15), Fact 23, Thm 24, Remark 25), 3.2.2 (Prop. 35, Lemma 36), 3.3 (Thm 42), 4 (Thms 46, 47).
- [AKV] N. Anari, F. Koehler, T.-D. Vuong, *Trickle-down in localization schemes and applications*, [2407.16104](https://arxiv.org/abs/2407.16104): §1, §3 (Eq. (7)–(12), Thm 53), §4.1 (Thm 54), App. A.1 (Lemmas 98–100, Thm 101).
- [AJKPV-U] N. Anari, V. Jain, F. Koehler, H. T. Pham, T.-D. Vuong, *Universality of spectral independence with applications to fast mixing in spin glasses*, [2307.10466](https://arxiv.org/abs/2307.10466): §1, §2 (Def. 18, 22, 27, 28, Fact 20, Lemma 21, Thm 23), §3 (Thms 30, 32), §4.1 (Thms 33, 36), §5.1 (Thm 43).
- [CFYZ] X. Chen, W. Feng, Y. Yin, X. Zhang, *Rapid mixing of Glauber dynamics via spectral independence for all degrees*, [2105.15005](https://arxiv.org/abs/2105.15005): §2 (field dynamics, Prop. 2.2, Lemmas 2.3, 2.4, 2.7, 2.8, Thm 2.5), §3.2, §4 (Lemma 4.1), §5.1.
- [CR] C.-F. Chen, C. Rouzé, *Quantum Gibbs states are locally Markovian*, [2504.02208](https://arxiv.org/abs/2504.02208): abstract, §I.A–B, §II, §III (Thm III.1, Remark III.1.1), §IV.C, §V, §VII.A (Lemmas VII.2, VII.3), §VIII.
- [Y] T. H. Yang, *Improved estimate of local Markovianity for quantum Gibbs states*, [2609.38007](https://arxiv.org/abs/2609.38007): abstract and AI usage statement, §I, Table I, §II (Thm II.1 and setting), parts of §V and §VI, §VII, references. (Read for the role of the approximate Markov property only; its proofs are not digested here.)
- [KB] K. Kato, F. G. S. L. Brandão, *Quantum approximate Markov chains are thermal*, [1609.06636](https://arxiv.org/abs/1609.06636): §I, §II (Thms 1–4, Cors. 5–6), §III.A–D.
- [BLMT] A. Bakshi, A. Liu, A. Moitra, E. Tang, *A Dobrushin condition for quantum Markov chains*, [2510.08542](https://arxiv.org/abs/2510.08542): §1, §2.1, §5 (Def. 5.1, Lemmas 5.2, 5.5, 5.8, Thms 2.1, 2.2, 5.7, Remark 5.9), §6 opening (Thm 2.3).
- [BCR] I. Bardet, Á. Capel, C. Rouzé, *Approximate tensorization of the relative entropy for noncommuting conditional expectations*, [2001.07981](https://arxiv.org/abs/2001.07981) (Ann. Henri Poincaré 2022): §1, §2.4, Def. 2, Thms 1–3, Prop. 4, §4.1, §5, §6.
- [HSL] R. Hayakawa, A. Southwell, C. M. G. Leditto, K.-C. Chen, M.-H. Hsieh, *Double localization for quantum Gibbs sampler gaps*, [2609.39802](https://arxiv.org/abs/2609.39802): abstract, §1, §2.1 (Thm 1), Table 1 (statements only).
- Listings consulted without opening the paper: Štefankovič–Vigoda et al. monograph [2307.13826](https://arxiv.org/abs/2307.13826) (abstract), Blanca et al. [2103.07459](https://arxiv.org/abs/2103.07459) (abstract), Bergamaschi [2606.26090](https://arxiv.org/abs/2606.26090) (title, abstract).

**Not retrieved in this digest.** Sibling files in this directory written in the same session cover what this one leaves out: `arxiv-2609.38007.md` (the paper), `chat-2609.38007-retrieval-status.md` (the shared chat) and `expanders.md` (expander graphs); they were not used for any statement here. Not retrieved for this digest: the ChatGPT share link given by the user for 2609.38007, the Semantic Scholar page of the Hammersley–Clifford manuscript (its statement is quoted from [Y] and [KB]), T. Tao's 245B notes, and the local file `C:\Users\User\Downloads\expander_survey.pdf` (a Windows path, not reachable from this Linux session). None was needed for the statements above; the user's own reading of the chat and the expander notes may add context not reflected here. The search for a non-commutative (quantum) version of spectral independence or of a $k$-level local-to-global theorem (alphaXiv discovery, two queries; Exa, one query) returned none; this is reported as "not found", not as "does not exist".

---

## 9. Numerical checks

Script: [`check_hdx_si.py`](check_hdx_si.py), numpy, seed 0, under a minute. Output of the run used for this note:

| id | statement | computed |
|---|---|---|
| C1 | [ALO] Thm 3.1: $\mathrm{spec}(P_\emptyset)=\mathrm{spec}(\Psi/(n-1))\cup\{-\frac1{n-1}\}^{n-1}\cup\{1\}$, three random full-support laws on $\{0,1\}^5$ | max deviation $6.7\cdot10^{-16}$ |
| C2 | U7: $D(P\Vert Q_{\text{Markov}})=\sum_iI(\sigma_{i+1};\sigma_{<i}\mid\sigma_i)$, random law on $\{0,1,2\}^4$ | $0.3554112364$ = $0.3554112364$ |
| C3 | §4.5: commuting-square defect $\Vert E_{\le1}E_{\ge1}f-E_1f\Vert _\infty$ | generic law $0.20$; its Markov projection $1.1\cdot10^{-16}$ |
| C4 | §4.4: heat-bath on $A=\{2,3\}$ of a 6-site Ising chain recovers $\mu$ from $\mu_{A^c}\otimes\mathrm{unif}_A$; limit depends on $\sigma$ only through sites $1,4$; Dirichlet form of the time average $\le c/t$ | $\Vert \cdot\Vert _1=3.6\cdot10^{-14}$; locality true; $t\cdot\mathcal E=0.196,0.174,0.046,0.014$ at $t=1,4,16,64$ (below the gap-free $0.41$) |
| C5 | §3.5: ALO lower bound vs true Glauber gap, Ising chain $n=6$, $\beta=0.4$ | gap $0.0660$, bound $0.0358$; $\eta_i=0.930,0.864,0.748,0.590,0.363$ |

These check algebra and identities, not claims about any realization.
