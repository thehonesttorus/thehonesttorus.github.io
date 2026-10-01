# Local-to-global unlocks

### What high-dimensional expanders and tiling algebras actually buy, and what the barycentre-arrow framework would need in order to buy the same

*Prong 3 of the [research programme](research-program.md). Written at the abstract level: nothing below refers to a particular realization. Inputs from the other prongs used here: the objects of Note 1 ([conditional arrow algebra](conditional-arrow-algebra.md)) and Note 2 ([simplicial complex as decomposition](simplicial-complex-as-decomposition.md)), as vocabulary only. Status labels: **Fact** (cited, with the statement as the source gives it), **Reduction** (the complexity gain, stated as global object / local certificate / explicit bound / cost), **Conjecture** (what would have to be true in the framework), **Guard** (what does not transfer).*

---

## 0. The hunt

The question for this prong is not "is there a random walk" but "what is the exact mechanism by which a local condition, checked on small pieces, controls a global quantity with an explicit error, and what work does that make unnecessary". Three literatures do this with theorems:

1. **expander graphs**: one number (the spectral gap) turns global averages into local steps;
2. **high-dimensional expanders (HDX)**: the number becomes local itself (spectra of links), and the walk becomes a walk on faces whose steps are restriction and co-restriction;
3. **tiling algebras**: the invariants of an uncountable, aperiodic hull are computed from a finite substitution matrix and a direct limit.

The pattern that survives all three is in §4; the candidate unlocks for the framework are in §5; what does not transfer is in §6.

---

## 1. Expander graphs: the prototype reduction

**Fact.** For a reversible chain with transition operator $P$, stationary law $\pi$ and second eigenvalue $\lambda_2<1$, $\|P^tf-\pi(f)\|_{L^2(\pi)}\le\lambda_2^{\,t}\|f-\pi(f)\|_{L^2(\pi)}$. The number of steps needed to estimate a global average $\pi(f)$ to accuracy $\varepsilon$ is $O(\log(1/\varepsilon)/(1-\lambda_2))$ and does not depend on the size of the state space (Lubotzky [1712.02526](https://arxiv.org/abs/1712.02526), §1, for the equivalent definitions; the expander mixing lemma for the "could be anywhere" statement: edge counts between any two sets are what a random graph would have, up to $\lambda_2\sqrt{|S||T|}$).

**Reduction 1 (global average from local steps).** Global object: $\pi(f)$, a sum over the whole space. Local certificate: $\lambda_2$, one number. Explicit bound: $\lambda_2^{\,t}$. Cost: $t$ evaluations of $P$, each local. The walk stays in a small region per step, but because the stationary law is unique and the gap is uniform, the region it ends in is representative: *globality is not where the walk is, it is that it could be anywhere with the right law*.

The weakness, and the reason HDX exists: for a graph, $\lambda_2$ is itself a global quantity. Certifying it means looking at the whole graph.

---

## 2. High-dimensional expanders: the certificate becomes local

Setting: a pure $d$-dimensional simplicial complex $X$ (every face is in a $d$-face). The **link** of a face $\alpha$ is $X_\alpha=\{\beta\setminus\alpha:\alpha\subset\beta\in X\}$; $G_\alpha$ is the weighted graph underlying the link (vertices $=$ faces $\alpha\cup\{v\}$, edges $=$ faces $\alpha\cup\{v,w\}$, weighted by the number of $d$-faces containing them). $X$ is a **$\gamma$-local spectral expander** if $\lambda_2(G_\alpha)\le\gamma$ for every face $\alpha$ of dimension $\le d-2$ (one-sided: only the second largest eigenvalue is controlled; two-sided: also the smallest). The **down-up walk** $P^\triangledown_k$ on $k$-faces removes a uniformly random vertex and adds a uniformly random one back (restriction, then co-restriction); the up-down walk $P^\triangle_{k-1}$ is the other order.

### 2.1 Garland: global cohomology from local spectra

**Fact** (Garland 1973, as summarised in Lubotzky [1712.02526](https://arxiv.org/abs/1712.02526) Thm 2.5). If $\dim X=d$ and for every face $F$ of dimension $d-2$ the link $\ell k_X(F)$ (a graph) has normalised spectral gap $\lambda^{(0)}(\ell k_X(F))\ge\varepsilon$, then $\lambda^{(d-1)}(X)\ge1+d\varepsilon-d$. For $\varepsilon>1-1/d$ this is positive, so $H^{d-1}(X,\mathbb R)=0$; applying it to skeleta gives all $H^j$, $1\le j\le d-1$. Garland proved Serre's conjecture (vanishing of real cohomology of cocompact lattices in $p$-adic groups of rank $\ge2$) this way, using that the vertex links of Bruhat–Tits buildings are spherical buildings over $\mathbb F_q$, whose spectral gaps are known and tend to $1$ as $q\to\infty$ (Lubotzky §2.3).

**Reduction 2 (topology from links).** Global object: $H^j(X,\mathbb R)$. Local certificate: the spectral gaps of the codimension-2 links, which are *graphs*, finite and unweighted. Explicit bound: $1+d\varepsilon-d$. Cost: one eigenvalue computation per small link, instead of the rank of a global coboundary matrix. The mechanism is an inequality between the global Laplacian and a sum of link Laplacians (Garland's "$p$-adic curvature").

### 2.2 Oppenheim: trickling down, so only the top links need checking

**Fact** (Oppenheim [1407.8517](https://arxiv.org/abs/1407.8517), as stated in Alev–Lau [2001.02827](https://arxiv.org/abs/2001.02827) Thm 1.3 and Cor 1.4). Let $X$ be pure $d$-dimensional. If $\lambda_2(G_\beta)\le\gamma\le\tfrac12$ for every face $\beta$ of dimension $k$ and $G_\alpha$ is connected for every face $\alpha$ of dimension $k-1$, then $\lambda_2(G_\alpha)\le\gamma/(1-\gamma)$ for every face $\alpha$ of dimension $k-1$. Iterating: if $\lambda_2(G_\beta)\le\gamma\le1/d$ for all $(d-2)$-faces and all links are connected, then every $k$-face $\alpha$ has $\lambda_2(G_\alpha)\le\gamma/(1-(d-2-k)\gamma)$.

**Reduction 3 (one level of links certifies all levels).** The local spectral expansion of *every* link, exponentially many of them with implicit weights, is certified by the links of the top codimension, whose edge weights are $0$ or $1$, plus connectivity. Alev–Lau note the reduction is "basically lossless" in the regime $\gamma=O(1/d^2)$. This is the step that makes the local certificate *checkable*.

### 2.3 Kaufman–Oppenheim and Alev–Lau: walks on faces mix, with the gap written out

**Fact** (Kaufman–Oppenheim [1707.02799](https://arxiv.org/abs/1707.02799), Thm 1.4, one-sided version). On a pure $n$-dimensional one-sided $\lambda$-local spectral expander, for a $k$-cochain $\varphi$ orthogonal to the constants, $\|M_k^+\varphi\|\le\big((1-\tfrac1{k+2})+\tfrac{k+1}2\lambda\big)\|\varphi\|$. Their abstract states the obstruction exactly: the spectral gap of high-order walks is *inherently small* ($\approx1/(k+2)$) because of coboundaries lifted from lower dimensions, so "beyond spectral gap" means decomposing a $k$-cochain into parts coming from dimensions $j\le k$, each shrunk by $\tfrac{k+1-j}{k+2}$ plus an error in $\lambda$. Negative eigenvalues of the links do not matter: one-sided expansion suffices, which is what Ramanujan complexes and matroid complexes provide.

**Fact** (Alev–Lau [2001.02827](https://arxiv.org/abs/2001.02827), Thm 1.5 and Cor 1.6). With $\gamma_j:=\max\{\lambda_2(G_\alpha):\dim\alpha=j\}$, for $0\le k\le d$,
$$
\lambda_2(P^\triangledown_k)=\lambda_2(P^\triangle_{k-1})\le1-\frac1{k+1}\prod_{j=-1}^{k-2}(1-\gamma_j).
$$
In particular if $\gamma_{k-2}\le1/(k+1)$ and the links are connected up to dimension $k-2$, then $\lambda_2(P^\triangledown_k)\le1-1/(k+1)^2$. The earlier Kaufman–Oppenheim bound was $1-\tfrac1{k+1}+\tfrac{k\gamma}2$, which needs $\gamma<2/k^2$ to say anything.

**Reduction 4 (sampling a facet from local gaps).** Global object: the uniform (or weighted) law on the $k$-faces. Local certificate: the $\gamma_j$. Explicit bound: the product formula; the mixing time is $O\big((k+1)^2\log(1/\varepsilon\pi_{\min})\big)$ under the corollary. Cost: single-vertex exchanges. The down-up walk is the composite *restriction $\circ$ co-restriction*; nothing in the proof uses a global coordinate.

### 2.4 Anari–Liu–Oveis Gharan–Vinzant: counting from sampling, with the local condition inherited for free

**Fact** ([1811.01816](https://arxiv.org/abs/1811.01816), Thm 1.1). If $\mu$ is a $d$-homogeneous strongly log-concave distribution on subsets of $[n]$ then its down-up chain $P_\mu$ has, for each $0\le k\le d-1$, at most $|X(k)|$ eigenvalues above $1-\tfrac{k+1}d$; hence spectral gap $\ge1/d$ and mixing time from $\tau$ at most $d\log\frac1{\varepsilon\mu(\tau)}$. The uniform law on the bases of any matroid is strongly log-concave (their Section 5, self-contained; also Huh–Wang via Hodge theory, per Alev–Lau footnote 2), equivalently the matroid complex is a $0$-local spectral expander, and the reason is **hereditary**: every link of a matroid complex is a matroid complex (a contraction), so the local condition propagates by itself and trickling down applies with $\gamma=0$. Corollaries: an FPRAS for counting bases of a matroid given by an independence oracle; the Mihail–Vazirani conjecture (Thm 1.5: the bases-exchange graph has expansion $\ge1$).

**Reduction 5 (enumeration replaced by a walk whose steps are the structure's own restrictions).** Global object: $|\mathcal B(M)|$, exponential to enumerate, hard to count exactly. Local certificate: a hereditary property (strong log-concavity of the generating polynomial, equivalently log-concavity of every contraction's quadratic form). Explicit bound: gap $1/d$. Cost: polynomially many single-element exchanges, plus the Jerrum–Valiant–Vazirani counting-to-sampling reduction. This is the sharpest example of the theme: *the proof never looks at the whole complex because the condition is closed under taking links*.

**Fact** (Anari–Liu–Oveis Gharan [2001.00303](https://arxiv.org/abs/2001.00303), abstract). A distribution $\mu$ is *spectrally independent* if an associated correlation (influence) matrix has bounded largest eigenvalue for $\mu$ and all of its conditional distributions; spectral independence implies that the associated simplicial complex is a local spectral expander, hence (via Kaufman–Mass, Dinur–Kaufman, Kaufman–Oppenheim, Alev–Lau) that the Glauber dynamics mixes rapidly. Application: Glauber dynamics for the hardcore model mixes in polynomial time up to the uniqueness threshold, improving Weitz's quasi-polynomial correlation-decay algorithm. The local certificate here is correlation decay, a property of pairs, checked under every pinning.

### 2.5 Dinur–Kaufman: global consistency from linearly many local views

**Fact** (Dinur–Kaufman, *High dimensional expanders imply agreement expanders*, FOCS 2017, pp. 974–985, [doi:10.1109/FOCS.2017.94](https://doi.org/10.1109/FOCS.2017.94); ECCC [TR17-089](https://eccc.weizmann.ac.il/report/2017/089/)). A local assignment on a $d$-dimensional complex $X$ gives each $d$-face $s$ a function $f_s:s\to\{0,1\}$. It is *global* if some $g:X(0)\to\{0,1\}$ restricts to all $f_s$. $X$ is a $c$-agreement expander if, for the test distribution $D$ on pairs of faces, $\mathrm{agree}_D(f)\ge1-\varepsilon$ implies there is a global $g$ with $\Pr_s[f_s=g|_s]\ge1-\varepsilon/c$ (Def. 1.1). Theorem 1.6: for $\lambda$ small enough in terms of $d$, a $\lambda$-HD expander of dimension $d^2$ has $k$-skeleta ($k\le d$) that are $c$-agreement expanders, with $c$ an absolute constant; Theorem 1.3: explicit bounded-degree agreement expanders exist (from the Ramanujan complexes of Lubotzky–Samuels–Vishne), with $O_d(n)$ faces where the complete complex, the only previously known case, has $\approx n^{d+1}$.

**Reduction 6 (gluing from sparse overlaps).** Global object: whether a family of local functions comes from one global function, and if not, how far it is. Local certificate: pairwise agreement on overlapping faces, tested on a bounded-degree family. Explicit bound: $\varepsilon/c$. Cost: $O(n)$ local views instead of $n^{d+1}$. This is the local-to-global statement behind PCPs, and it is a statement about *restrictions*: $f_s$ is the restriction of a hypothetical $g$, and agreement of restrictions on overlaps is the sheaf condition, tested approximately and sparsely.

---

## 3. Tilings: invariants of an uncountable hull from a finite matrix

**Fact** (Kellendonk [cond-mat/9403065](https://arxiv.org/abs/cond-mat/9403065); Bellissard's gap-labelling, cited there). To a tiling one attaches the groupoid of translations on its hull and the C\*-algebra $\mathcal A_T$; a local Schrödinger operator on the tiling has, by Shubin's formula, integrated density of states equal to a trace per unit volume on spectral projections, and *if $E$ lies in a gap then* $\mathrm{IDS}(E)\in\mathrm{tr}_*\big(K_0(\mathcal A_T)\big)\cap[0,1]$ (his eq. (7), "part of the abstract gap labelling theorem of Bellissard"). For substitution tilings the self-similarity gives an AF subalgebra whose $K_0$ is the direct limit of the substitution matrices and whose trace is the Perron–Frobenius eigenvector: the module of *pattern frequencies*. For products of one-dimensional tilings this frequency module already exhausts the gap labels.

**Fact** (Bellissard–Benedetti–Gambaudo [math/0109062](https://arxiv.org/abs/math/0109062), abstract). For a large class of tilings including Penrose and the icosahedral ones, the continuous hull is a minimal lamination with flat leaves and Cantor transversal, and it is the projective limit of a sequence of branched, oriented, flat compact manifolds (finite telescopic approximations); the positive invariant measures form a convex cone canonically sitting in the projective limit of the top homology groups of these branched manifolds; a gap-labelling theorem follows. (Anderson–Putnam, ETDS 1998, is the substitution-tiling version: the hull as an inverse limit of a finite CW complex under the substitution map, with cohomology the direct limit.)

**Fact** (Connes, *Noncommutative Geometry* 1994, Ch. II §3, read for Note 1). The Penrose algebra is AF, built from the Bratteli diagram of the substitution; $K_0\cong\mathbb Z^2$ with the order determined by the golden ratio; the unique trace is the frequency measure, under which every finite past of a given symbol is equally likely.

**Reduction 7 (an aperiodic infinite object reduced to its substitution matrix).** Global object: the spectrum's gap structure of an operator on an infinite aperiodic system; the invariant measures of the hull; its $K$-theory. Local certificate: the finite list of patch types and the finite matrix saying how each patch decomposes under one substitution step. Explicit bound: none needed, the direct limit is exact; what is approximate is the finite telescopic stage, and successive stages are nested so the error is monotone. Cost: a Perron–Frobenius eigenvector and a direct limit of finitely generated groups.

The "could be anywhere" principle appears here as **unique ergodicity**: pattern frequencies exist uniformly over the tiling, so any patch sampled anywhere is representative, and the uniqueness of the trace on the AF algebra is the algebraic form of that statement. The expander gap is the quantitative version of the same thing (rate at which a local sample becomes representative); uniqueness of the trace is the qualitative version (that it does at all).

---

## 4. The common mechanism

Across the seven reductions the same four ingredients appear, and nothing else is used.

| ingredient | expanders / HDX | tilings | in the framework's vocabulary |
|---|---|---|---|
| **a hereditary class**: restricting an object gives an object of the same kind | links of a complex are complexes; links of matroid complexes are matroid complexes | a patch of a tiling is a tiling patch; a substitution step is again a substitution | faces of the layer complex, sub-faces, restrictions to a sub-frame (Note 2 §2); the link of a face |
| **a transport operator built from restriction and its adjoint** | down-up walk $=$ restrict then co-restrict; Glauber dynamics; the Laplacian $\delta^*\delta$ | substitution / de-substitution; the tail-equivalence relation | the arrows and their adjoints; Lüders conditioning and its dual (Note 2 §3); the derived walk and the transfer operator (Note 1 §5.3) |
| **a spectral or $K$-theoretic invariant of the composite controlled by local data, with a formula** | Garland $1+d\varepsilon-d$; Oppenheim $\gamma/(1-\gamma)$; Alev–Lau $1-\frac1{k+1}\prod(1-\gamma_j)$; ALOV $1/d$; DK $\varepsilon/c$ | $K_0=\varinjlim$ of substitution matrices; trace $=$ PF eigenvector; IDS at gaps $\in\mathrm{tr}_*K_0$ | the KMS/Gibbs states and the gap of the transfer operator; the cocycle class of $F$ (Note 1 §6); the dimension group of the endpoint AF algebra (Note 1 §3.4) |
| **homogeneity**: a unique stationary law / unique trace, so a local sample is globally representative | expander mixing lemma; unique stationary law | unique ergodicity; unique trace | centrality of the state $\iff$ $F$ a coboundary (Note 1 §6); uniqueness of KMS states |

The profound reduction is always the second row acting on the first: *because the class is hereditary, the operator that moves one step is also the operator that passes to a smaller object of the same kind, so a bound on small objects is a bound on one step, and iterating one step is all a global computation ever needs*.

Two further observations that are not in the sources but follow from reading them side by side.

- **The obstruction is always a coboundary.** Kaufman–Oppenheim's "inherently small gap" comes from cochains lifted from lower dimensions; Garland's bound is on the Laplacian orthogonal to coboundaries; the gap labels are classes in $K_0$, i.e. projections modulo equivalence; and in Note 1 the failure of sufficiency is exactly a non-trivial cohomology class of the cocycle $F$. The local-to-global theorems are theorems about *what is left after quotienting by coboundaries*.
- **One-sidedness matches layering.** A layered (partite) complex has links that are partite, hence have an eigenvalue $-1/d$ (Dinur–Kaufman's remark about the LSV complexes); only one-sided local spectral expansion can hold. Kaufman–Oppenheim and ALOV show the one-sided hypothesis is enough for optimal mixing. So the HDX theory that applies to layered structures is precisely the one-sided theory, and it loses nothing.

---

## 5. Candidate unlocks for the framework (conjectures, with what must be true)

Each is stated for the abstract objects of Notes 1–2: a layered family of complexes $K_\ell$, the quiver $\Lambda$ of arrows between faces, the history complex and the coherence complex, a cocycle $F$, states and their restrictions. None mentions a realization.

**U1 (local gap $\Rightarrow$ global state from finite-depth statistics).** *Conjecture.* If the layer complexes $K_\ell$ and the history complex are one-sided $\gamma$-local spectral expanders (checked, by trickling down, on codimension-2 links only), then the transfer operator of the derived walk (Note 1 §5.3) has a spectral gap bounded below by an Alev–Lau-type product, and the KMS/Gibbs state on complete histories is determined up to $\lambda^t$ by the statistics of histories of length $t$. *What must be true:* that the derived walk is a composite of restriction and co-restriction in a pure complex (it is, inside a layer: exclusion and inclusion of a vertex; across layers it is an arrow, and the history complex is what makes consecutive layers faces of one complex, Note 2 §6). *Work for prong 1:* define the links of the history complex (links of a chain in an order complex are joins of interval complexes) and compute what local spectral expansion of a join means in terms of the factors.

**U2 (hereditary face structure $\Rightarrow$ enumeration replaced by conditioning steps).** *Conjecture.* If a layer complex $K_\ell$ is a matroid complex, or more generally its facet-generating polynomial is strongly log-concave, then the uniform law on its facets, equivalently the maximally coherent states $\omega_\sigma$ of the facet algebra (Note 2 §3), is sampled by the chain whose steps are *Lüders exclusion of a vertex followed by its adjoint inclusion*, with spectral gap $\ge1/\dim$. *Why it would be an unlock:* counting or integrating over the faces of a layer is replaced by polynomially many conditioning operations, which are the framework's own arrows, and the proof uses nothing but the hereditary property. *What must be true:* an exchange property for faces, i.e. that from two facets one can always move one vertex from the larger to the smaller side. Whether the complexes that arise as nerves of restrictions (Note 2 §1) have this property is a question about restrictions, not about coordinates.

**U3 (coboundary expansion $\Rightarrow$ sufficiency certified locally).** Note 1 §6 says the barycentre is sufficient iff $F$ is a coboundary. Coboundary expansion (Linial–Meshulam, Gromov; Lubotzky §3.2, Def. 3.3) is the inequality $\operatorname{dist}(F,B^1)\le\|\delta F\|/h^1(X)$ for 1-cochains: the global distance from the coboundaries is bounded by a *local* quantity, the coboundary $\delta F$ evaluated on 2-cells, over an expansion constant. *Conjecture.* For the 2-complex whose 2-cells are the commuting squares of the layered quiver (two arrows out of $\sigma$ and two arrows into $\tau''$ forming a square through $\tau,\tau'$), $\delta F$ on a square is the alternating sum $F(\sigma\to\tau)+F(\tau\to\tau'')-F(\sigma\to\tau')-F(\tau'\to\tau'')$, and if that square complex is an $h$-coboundary expander then the failure of sufficiency of the barycentre, a global quantity, is at most $\max_{\text{squares}}|\delta F|/h$. *Why it would be an unlock:* insufficiency would be measured on four arrows at a time. *Work for prong 1:* identify the square complex of the quiver with the 2-skeleton of the path-space groupoid's classifying space, where the class of $F$ lives, and check whether its coboundary expansion is forced by the layered structure.

**U4 (agreement expansion $\Rightarrow$ a state from linearly many restrictions).** *Conjecture.* If the coherence complex $K(C(V)\subset A)$ of Note 2 §1.3 is a $c$-agreement expander, then a state on $A$ is determined, up to $\varepsilon/c$, by its restrictions to a bounded-degree family of faces that agree pairwise on overlaps up to $\varepsilon$; in particular reconstruction needs $O(|V|)$ restrictions rather than all faces. *What must be true:* an agreement test for *states* rather than $\{0,1\}$-valued functions. Abramsky–Brandenburger's sheaf formulation of contextuality (Note 2 §5) is the exact obstruction when the gluing fails; the Dinur–Kaufman theorem is the statement that on an expander the obstruction is detectable sparsely. The precise analogue to look for is a "noncommutative agreement expander": the test distribution on pairs of overlapping faces, and the fidelity of the two restricted states on the overlap.

**U5 (the history algebra as a finite telescopic approximation).** *Conjecture.* When the layered structure is repeated (the same incidence pattern at every depth), the history algebra of Note 1 is an inductive limit of the finite algebras $\bigoplus_\sigma M_{N_\ell(\sigma)}$ with the incidence matrix of the quiver as the connecting map, exactly as the Penrose algebra is the limit of its substitution matrices. Then: $K_0$ is the dimension group of the incidence matrix; the central (tracial) states are the Perron–Frobenius data; uniqueness of the trace is the "could be anywhere" property of the derived walk; and the class of the cocycle $F$ pairs with $K_0$ the way the integrated density of states pairs with the trace, giving a *gap labelling of the sufficiency defect*: a countable set of values that the defect can take, computed from the finite incidence matrix. *What must be true:* repetition of the incidence pattern (or at least eventual periodicity of the Bratteli diagram), which is the one assumption that makes a direct limit computable. Note 1 §3.4 already identifies the Bratteli diagram; what is missing is the stationarity.

**U6 (correlation decay on the frame $\Rightarrow$ mixing of the vertex walk).** *Conjecture.* The spectral-independence criterion transfers verbatim to the frame $C(V)$ of a layer: if the influence matrix of the law on faces (how conditioning on vertex $u$ being in or out of the face changes the marginal of $v$) has bounded largest eigenvalue for the law and all its pinnings, then the single-vertex walk on faces mixes in time polynomial in $|V|$. *Why it matters:* pinning is exactly restriction to a sub-frame (Note 2 §2), so the hypothesis is a statement about how restrictions act on restrictions, which is the question the framework was built to ask first.

---

## 6. Guards: what does not transfer, and the discipline for not being dragged

- **Purity and dimension.** HDX theorems are for pure complexes; a layer complex of restrictions need not be pure (facets of different dimension, which is the "dimension can change by more than one" feature). Either work with the pure part of each dimension, or use the weighted theory (Alev–Lau's $G_\alpha$ are weighted graphs and tolerate non-uniform facets). This is a question for prong 1, not a reason to change the objects.
- **Acyclicity.** The quiver of arrows is acyclic and layered; walks on it are not reversible chains. The HDX machinery applies to the *within-layer* restriction structure and to the history complex; the cross-layer walk is the Doob transform of Note 1 §5.3 and its "mixing" is convergence of the backward partition function, not of a reversible chain. Keep the two separate.
- **Bounded degree is the point, not an obstacle.** The only reason Dinur–Kaufman's $O(n)$ faces is remarkable is that the complete complex has $n^{d+1}$. A framework whose face complexes are small has nothing to gain from U4; the gain is exactly for frames with many vertices and few relevant faces.
- **One-sided only.** Because of partiteness, two-sided local expansion is impossible for layered structures; any transfer must go through the one-sided theorems (Kaufman–Oppenheim Thm 1.4, ALOV, Alev–Lau), never through the two-sided decomposition theorem.
- **No coordinates were used.** Every reduction above is stated for abstract complexes, groupoids and algebras. This is why prong 3 can speak to prong 1 directly. It is also the test for whether a proposed transfer is being dragged toward a realization: if a statement needs a carrier space, a metric or an embedding to be made, it belongs to prong 2 as a hypothesis about a dictionary, not here.

---

## 7. Messages to the other prongs

**To prong 1 (theory).** Three definitions are needed before any of U1–U6 can be proved or refuted: (i) the links of the history complex and the meaning of one-sided local spectral expansion for a join of interval complexes; (ii) the square 2-complex of the layered quiver and its identification with the cochain complex in which $[F]$ lives, with a coboundary-expansion constant; (iii) the stationarity (repetition) hypothesis under which the history algebra is a direct limit with one incidence matrix. The suggested order is (ii), since Note 1 §6 already defines the quantity it would control.

**To prong 2 (bridge).** The quantities that prong 3 would like measured on any realization, because they are the certificates in the reductions above and are coordinate-free: the second eigenvalues $\gamma_j$ of the links of each layer complex (codimension 2 first); the alternating sums of $F$ around squares of the quiver; the number of faces versus the number of faces needed for agreement expansion; whether consecutive layer complexes repeat their incidence pattern. A negative result on any of these is a statement about the dictionary used to read the realization (programme rule P2.1), not about the reductions, which are theorems.

---

## Sources read for this note

- A. Lubotzky, *High dimensional expanders*, [1712.02526](https://arxiv.org/abs/1712.02526) (ICM 2018 survey): §1 expander definitions; §2.2 Garland's method, Thm 2.5; §2.3 Bruhat–Tits buildings and Ramanujan complexes; §3.2 coboundary expansion, Def. 3.3, Thm 3.5.
- I. Oppenheim, *Local spectral expansion approach to high dimensional expanders Part I: Descent of spectral gaps*, [1407.8517](https://arxiv.org/abs/1407.8517), Discrete Comput. Geom. 2018 (the trickling-down theorem, quoted here in the form given by Alev–Lau).
- T. Kaufman, I. Oppenheim, *High order random walks: beyond spectral gap*, [1707.02799](https://arxiv.org/abs/1707.02799) (abstract, Thms 1.3–1.5); *High dimensional expanders and coset geometries*, [1710.05304](https://arxiv.org/abs/1710.05304) (STOC 2018 as *Construction of new local spectral high dimensional expanders*): elementary, self-contained bounded-degree local spectral expanders from coset geometries of elementary-matrix groups, acting simply transitively on the top faces of a pure, partite clique complex.
- V. L. Alev, L. C. Lau, *Improved analysis of higher order random walks and applications*, [2001.02827](https://arxiv.org/abs/2001.02827), STOC 2020: Thms 1.2, 1.3, 1.5, Cors 1.4, 1.6, §1.2 applications.
- N. Anari, K. Liu, S. Oveis Gharan, C. Vinzant, *Log-concave polynomials II: high-dimensional walks and an FPRAS for counting bases of a matroid*, [1811.01816](https://arxiv.org/abs/1811.01816), STOC 2019: Thm 1.1, Thm 1.5, §1.4, §5.
- N. Anari, K. Liu, S. Oveis Gharan, *Spectral independence in high-dimensional expanders and applications to the hardcore model*, [2001.00303](https://arxiv.org/abs/2001.00303), FOCS 2020 (abstract read; the influence-matrix definition should be re-read in full before U6 is made precise).
- T. Kaufman, D. Mass, *High dimensional random walks and colorful expansion*, [1604.02947](https://arxiv.org/abs/1604.02947), ITCS 2017 (Thms 1.5, 1.8–1.10; the combinatorial "colorful expansion" route to rapid mixing).
- I. Dinur, T. Kaufman, *High dimensional expanders imply agreement expanders*, FOCS 2017, pp. 974–985, [doi:10.1109/FOCS.2017.94](https://doi.org/10.1109/FOCS.2017.94); ECCC [TR17-089](https://eccc.weizmann.ac.il/report/2017/089/): Def. 1.1, Lemma 1.5, Thms 1.3 and 1.6. (Not arXiv 1703.00766, which is a different paper.)
- H. Garland, *p-adic curvature and the cohomology of discrete subgroups of p-adic groups*, Ann. of Math. 97 (1973), via Lubotzky.
- J. Kellendonk, *Noncommutative geometry of tilings and gap labelling*, [cond-mat/9403065](https://arxiv.org/abs/cond-mat/9403065), Rev. Math. Phys. 1995: §1.2 (eq. (7)), §2, and the AF algebra of a substitution.
- J. Bellissard, R. Benedetti, J.-M. Gambaudo, *Spaces of tilings, finite telescopic approximations and gap-labelling*, [math/0109062](https://arxiv.org/abs/math/0109062), Comm. Math. Phys. 2006 (abstract).
- J. E. Anderson, I. F. Putnam, *Topological invariants for substitution tilings and their associated C\*-algebras*, Ergodic Theory Dynam. Systems 18 (1998).
- A. Connes, *Noncommutative Geometry*, Ch. II §3 (Penrose tilings), as read for Note 1.
