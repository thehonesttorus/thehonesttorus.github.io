# The conditional arrow algebra

### A noncommutative-geometric home for abstract barycentre arrows, their projections, and their states

*Research note, 2026-10-01. Continues the "barycentre arrows" thread: objects are abstract barycentres of abstract simplices; an arrow is a linear layer together with the ReLU immediately before it; functions on arrows, not the weights themselves, form the algebra.*

---

## 0. What this note does

The standing conjecture is:

> preceding ReLU information + linear layer $\;\rightsquigarrow\;$ an arrow between barycentres of abstract simplices, carrying projection data (which vertices are included, hence the face dimension, which may jump arbitrarily).

The question asked was: go deep into why noncommutative geometry (NCG) exists, and find the natural home for the objects (barycentres), the arrows, the functions on arrows with their convolution product (a C\*-algebra or perhaps only an operator system), the projections, and the states, together with the way projections and states interact.

The answer developed here, in one sentence: **the natural home is the groupoid/semigroupoid C\*-algebra of the layered quiver of faces, viewed through Connes's measure-theoretic layer of NCG; there the arrows are partial isometries between face projections, states are measures on histories with coherences, conditioning by an arrow is the Murray–von Neumann transport of states, the network weights enter only as a real cocycle on arrows, the states adapted to that cocycle are KMS (Gibbs) states whose modular flow is layer depth, and "is the current barycentre a sufficient statistic?" becomes "is the cocycle a coboundary?".** Finite resolution (depth cut-offs, overlap tolerance between faces) produces operator systems, not algebras, exactly in the sense of Connes–van Suijlekom.

Everything below is labelled as one of:

- **Fact** — a theorem from the literature, with the source.
- **Derivation** — proved here from those facts (short proofs given); the finite-dimensional identities were also checked numerically (Section 9).
- **Expectation** — a stated conjecture or a modelling choice not yet fixed.

No coordinates are chosen for barycentres anywhere. No Euclidean carrier is assumed. The weights are never identified with elements of the arrow algebra.

---

## 1. Why noncommutative geometry exists, in Connes's own terms

This section is a faithful précis of Connes, *Noncommutative Geometry* (1994), Introduction and Chapter II §3, and of his later surveys ([Year 2000](https://arxiv.org/abs/math/0011193), [spectral standpoint](https://arxiv.org/abs/1910.10407)). The point of rehearsing it is that each motivation below reappears, essentially verbatim, in the barycentre-arrow problem.

### 1.1 Spaces that classical tools cannot see

Gelfand duality identifies a compact space $X$ with the commutative C\*-algebra $C(X)$. Connes's starting observation is that there is an abundant class of spaces for which this duality collapses: the space of Penrose tilings, the leaf space of a foliation, the dual of a discrete group, the orbit space of a group action. They are all *quotients* $Y/\!\sim$ of a good space by an equivalence relation whose classes are dense. The quotient topology is then trivial, measure theory is trivial, and (Connes, Chapter I Appendix C) even the effective cardinality is wrong: "any explicitly constructed map from such a set to the real line fails to be injective".

There are two ways to quotient (Connes, *Year 2000*, §IV):

1. keep only functions constant on classes, $A=\{f: f(a)=f(b)\}$ — this is the classical quotient and, for dense classes, gives $\mathbb C$;
2. keep both points and let them "speak to each other": adjoin to functions on $\{a,b\}$ the identification of $a$ with $b$, obtaining $M_2(\mathbb C)$.

The second operation applied to an equivalence relation, or better a **groupoid**, gives its convolution algebra. "The noncommutativity of the algebra betrays the identification of points of the ordinary space, and when the quotient happens to exist as an ordinary space, the algebra is Morita equivalent to the commutative algebra of functions on the quotient." Connes insists that groupoids are not a luxury: Heisenberg replaced the frequency *group* of a classical system by the *groupoid* of quantum transitions, and imitating the group convolution algebra for that groupoid *is* matrix multiplication.

This is exactly the situation of a layered network of faces: two histories (sequences of faces, one per layer) that currently sit at the same face are "identified" by the barycentre label, but they are not equal. Option 1 (functions of the current face only) discards the distinction. Option 2 keeps it as a groupoid of re-routings and gives a noncommutative algebra. The whole of Section 3 is option 2.

### 1.2 The hierarchy: measure theory first

Connes organises NCG by the classical hierarchy measure theory → topology → differential → metric, and stresses that the new phenomena appear already at the coarsest level:

- **Measure theory = von Neumann algebras.** Type I is classical measure theory up to multiplicity. Type II exhibits *continuous dimension*: subspaces are classified by a real number, with $\dim(E\wedge F)+\dim(E\vee F)=\dim E+\dim F$ still holding (von Neumann's continuous geometry). Type III is reduced to type II by the modular theory. And the headline: **"Noncommutative measure spaces evolve with time!"** — a von Neumann algebra $M$ has a canonical homomorphism $\delta:\mathbb R\to \mathrm{Out}(M)$, the class of the modular flow $\sigma^\varphi_t$ of any faithful normal state $\varphi$, independent of $\varphi$ (Connes's thesis, via Tomita–Takesaki and the Radon–Nikodym cocycle theorem). The relation between the KMS condition of statistical mechanics and the modular operator "remains one of the deepest points of contact between physics and pure mathematics".
- **Continuous dimension as density.** For a measured foliation, the Murray–von Neumann dimension of a projection is the transverse measure of a transversal, i.e. the *density* of a discrete set on the generic leaf: "the continuous dimensions bear the same relation to ordinary dimensions as densities bear to the counting of finite sets". For Penrose tilings it is the frequency of a motif.
- **Topology = C\*-algebras and K-theory.** Projections (self-adjoint idempotents) are the noncommutative characteristic functions; finite projective modules are the vector bundles (Serre–Swan); $K_0$ is built from equivalence classes of projections, exactly as in Murray–von Neumann; for AF algebras the ordered $K_0$ (dimension group) is a complete invariant (Bratteli–Elliott).
- **Metric = spectral triple.** The distance is a supremum over functions rather than an infimum over paths, $d(\varphi,\psi)=\sup\{|\varphi(a)-\psi(a)|:\ \|[D,a]\|\le 1\}$, and therefore makes sense for discrete, non-connected spaces and *between states*. Connes's two-point example: $A=\mathbb C\oplus\mathbb C$, $D=\begin{pmatrix}0&\mu\\ \mu&0\end{pmatrix}$ gives distance $\mu^{-1}$.

### 1.3 The Penrose tiling algebra, because it is a history space

Connes's opening example (Chapter II §3) is built from sequences: a tiling corresponds to a sequence $z=(z_n)$ of 0's and 1's with $z_n=1\Rightarrow z_{n+1}=0$; two sequences give the same tiling iff they agree from some index on. So the space of tilings is $X=K/\mathcal R$ with $K$ a Cantor set of admissible sequences and $\mathcal R$ *eventual agreement*. Every class is dense, the quotient topology is trivial.

The algebra $A$ is the convolution algebra of $\mathcal R$: matrices $(a_{z,z'})$ indexed by pairs of eventually-agreeing sequences. Connes then exhibits it as an AF algebra: on the finite set $K_n$ of admissible words of length $n+1$ take the relation $\bar{\mathcal R}_n=\{(z,z'): z_n=z'_n\}$, "same current symbol"; functions on $\bar{\mathcal R}_n$ form $A_n=M_{k_n}(\mathbb C)\oplus M_{k'_n}(\mathbb C)$, where $k_n,k'_n$ count words ending in 0, resp. 1; the inclusions $A_n\subset A_{n+1}$ are given by $k_{n+1}=k_n+k'_n$, $k'_{n+1}=k_n$, i.e. by the matrix $\begin{pmatrix}1&1\\1&0\end{pmatrix}$. The invariants: $K_0(A)=\mathbb Z^2$ ordered by the half-plane of the golden ratio; a unique trace $\tau(a)=\int_K a_{z,z}\,d\mu(z)$ whose values on projections are $(\mathbb Z+\tfrac{1+\sqrt5}{2}\mathbb Z)\cap\mathbb R_+$; weak closure the hyperfinite II$_1$ factor.

Read with the present problem in mind: $K_n$ is "the complex at layer $n$" (two faces), a word is a *history of faces*, $A_n$ is the algebra of operators that mix histories agreeing in their current face, the block sizes $k_n,k'_n$ are *numbers of histories arriving at a face*, and the unique trace is the unique measure on histories that is invariant under re-routing of the past. This is a Bratteli diagram, and Section 3.4 says that a layered network of complexes is one too.

### 1.4 What operator systems add (Connes–van Suijlekom)

Two recent papers extend the basic construction to the case where the relation is **not transitive** ([Spectral truncations](https://arxiv.org/abs/2004.14115), [Tolerance relations](https://arxiv.org/abs/2111.02903)). Two sources of non-transitivity are named: (i) a spectral cut-off $P$ turns the algebra $A$ into $PAP$, which is a self-adjoint subspace of $B(PH)$ but not an algebra; (ii) a *tolerance relation* (reflexive, symmetric, not transitive) such as $d(x,y)<\varepsilon$, or "two vertices of a simplicial complex are joined by an edge" — they explicitly mention simplicial complexes without the Kan property. In both cases one gets an **operator system** $E$: matrix-ordered, with a positive cone, states, pure states, and a distance formula, but no product. Two invariants measure the distance to an algebra: the C\*-envelope $C^*_{\mathrm{env}}(E)$ and the **propagation number** $\mathrm{prop}(E)$, the least $n$ such that products of $\le n$ elements of $E$ span a C\*-algebra. Facts used below:

- For a reflexive symmetric relation $R$ on a finite set $X$, $E(R)=\mathrm{span}\{e_{xy}:(x,y)\in R\}\subset M_{|X|}$; the $0$–$1$ matrix of $R$ is positive semidefinite **iff** $R$ is an equivalence relation, in which case $E(R)$ is a subalgebra (Lemma 3.8, Cor. 3.9).
- If the graph of $R$ is connected, $C^*_{\mathrm{env}}(E(R))=M_{|X|}$ and $\mathrm{prop}(E(R))=\mathrm{diam}(R)$ (Prop. 3.7).
- A vector state restricted to $E(R)$ is pure iff $R$ restricted to the support of the vector generates the full equivalence relation on that support (Prop. 3.11).
- Band matrices of width $N$ in $M_p$ have propagation number $\lceil p/N\rceil$ (Example 4.9); for $d(x,y)<\varepsilon$ on a path metric space of diameter $\delta$, $\mathrm{prop}=\lceil\delta/\varepsilon\rceil$ (Theorem 4.11).

So "C\*-algebra or operator system" has a precise answer: an algebra exactly when the composition relation is an equivalence relation (or a category closed under composition); an operator system whenever composition is cut off by resolution, depth, or tolerance. Both carry states and projections-as-effects.

### 1.5 Where barycentres already live, coordinate-free

Choquet theory does not need a metric. **Fact** (Kadison 1951; Alfsen 1971; see the [nc Choquet survey](https://arxiv.org/abs/2412.09455)): every compact convex set $K$ is affinely homeomorphic to the state space of the order-unit space $A(K)$ of continuous affine functions on it; a compact convex set is the state space of a *commutative* unital C\*-algebra iff it is a Bauer simplex (Choquet simplex with closed extreme boundary); every point of $K$ is the barycentre (resultant) $\int x\,d\mu(x)$ of a boundary measure, uniquely iff $K$ is a simplex. **Fact** (Kennedy–Shamovich, [nc Choquet simplices](https://arxiv.org/abs/1911.01023)): the noncommutative version replaces $K$ by the nc state space of an operator system, and nc Bauer simplices are exactly nc state spaces of C\*-algebras.

Consequence for the objects: an abstract simplex with vertex set $V$ is the state space $S(C(V))=\Delta(V)$; a face $\sigma\subset V$ is the face $\Delta(\sigma)=S(C(\sigma))$, the states supported on the hereditary subalgebra $p_\sigma C(V)$; "the barycentre of $\sigma$" is a distinguished interior point $b_\sigma\in\Delta(\sigma)$ (the resultant of a measure on the vertices of $\sigma$, e.g. the uniform one). Nothing metric was used; the locally convex structure is all. The spectral distance of §1.2 is an *optional* extra layer. The objects of the arrow category will be the *labels* $\sigma$; the barycentre $b_\sigma$ is what a state restricted to the face algebra sees (Section 4.1).

---

## 2. The objects and the arrows, stated abstractly

### 2.1 Data

For each layer $\ell=0,\dots,L$ there is an abstract simplicial complex $K_\ell$ on the vertex set $V_\ell$ (the units of the layer); its faces $\sigma\in K_\ell$ are the admissible "included-vertex" patterns. Put $E^0=\bigsqcup_\ell K_\ell$ (objects: one per face, equivalently one per abstract barycentre) and
$$
E^1=\bigsqcup_{\ell<L}\Lambda_\ell,\qquad \Lambda_\ell\subset K_\ell\times K_{\ell+1}\times(\text{multiplicity index}),
$$
an arrow $\gamma:\sigma\to\tau$ meaning: *the ReLU that produced the face $\sigma$ at layer $\ell$, followed by the linear layer, lands in the face $\tau$ at layer $\ell+1$.* Write $s(\gamma)=\sigma$, $r(\gamma)=\tau$. Each object carries the dimension label $d(\sigma)=|\sigma|$; an arrow carries the pair $(d(\sigma),d(\tau))$ with no constraint, which is the "projection information allowing jumps between non-adjacent dimensions". The quiver $\Lambda=(E^0,E^1,r,s)$ is finite, layered and acyclic; its path category is the semigroupoid of composable arrow words.

**Expectation (the one open modelling decision).** The weights determine *which* arrows exist (which $\tau$ are attainable from $\sigma$) and possibly multiplicities; how to read $\Lambda_\ell$ off $W_{\ell+1}$ is not fixed here. Everything below is conditional on a choice of $\Lambda$, and nothing below depends on coordinates.

### 2.2 Three different jobs the weights can do

The weights are not elements of the arrow algebra. In the standard constructions they can enter in exactly three distinct ways, each a known theory:

| role of the weights | object | theory |
|---|---|---|
| decide which arrows exist (and multiplicities) | the quiver $\Lambda$ | graph / semigroupoid C\*-algebras (Cuntz–Krieger; [Exel](https://arxiv.org/abs/math/0703182)) |
| a real number $F(\gamma)$ on each arrow ("action", "energy") | a cocycle $c_F$ on the path groupoid; a dynamics $\alpha^F$ | generalised gauge actions, KMS states ([aHLRS](https://arxiv.org/abs/1205.2194); [Christensen–Thomsen](https://arxiv.org/abs/1505.04751); [Neshveyev](https://arxiv.org/abs/1106.5912)) |
| a measure $\lambda_\sigma$ on the arrows leaving $\sigma$ | a C\*-correspondence over $C(E^0)$ | topological quivers ([Muhly–Tomforde](https://arxiv.org/abs/math/0312109)); transfer operators and interactions ([Exel](https://arxiv.org/abs/math/0012084), [Exel](https://arxiv.org/abs/math/0409267)) |

Sections 3–6 use the first two; the third is the Heisenberg/Schrödinger packaging of "questions go backward, states go forward" (Section 5.3).

### 2.3 Exel's semigroupoid conditions, checked against $\Lambda$

Exel's construction (Definition 2.1 and §19 of [math/0703182](https://arxiv.org/abs/math/0703182)) takes a *categorical* semigroupoid (a small category with identities removed), requires every element to be monic, intersecting pairs to have least common multiples, and the absence of *springs* (elements $f$ with $\Lambda^f=\emptyset$, i.e. nothing composable on their domain side). For the layered quiver: composition is unique, so monicity holds (equal composites have equal factors), and two arrows with a common extension are nested, so their least common multiple is the longer one; but every arrow leaving layer $0$ is a spring, since layer-$0$ faces receive nothing. Exel's Theorem 18.4 (semigroupoid algebra $\cong C^*$ of the groupoid of germs on the tight spectrum) therefore does not apply verbatim; the same algebra is reached by the graph-algebra route, which handles sources directly (an Huef–Laca–Raeburn–Sims, Kajiwara–Watatani). Exel's identification of the tight spectrum with the boundary path space is what survives: **the units of the groupoid are histories with a full past**, not faces.

---

## 3. Functions on arrows: the algebra and its three resolutions

Conventions follow Raeburn's book and [aHLRS](https://arxiv.org/abs/1205.2194): a path $\mu=\mu_1\mu_2\cdots\mu_n$ has $s(\mu_i)=r(\mu_{i+1})$, so $\mu_n$ is traversed first, $s(\mu)=s(\mu_n)$ is where the history starts and $r(\mu)=r(\mu_1)$ where it currently ends; $|\mu|=n$; vertices are paths of length $0$; $E^*$ is the set of finite paths.

### 3.1 The Toeplitz–Cuntz–Krieger family on histories

**Fact.** A Toeplitz–Cuntz–Krieger $\Lambda$-family consists of mutually orthogonal projections $\{p_\sigma\}_{\sigma\in E^0}$ and partial isometries $\{s_\gamma\}_{\gamma\in E^1}$ with
$$
s_\gamma^*s_\gamma=p_{s(\gamma)},\qquad \sum_{\gamma\in G}s_\gamma s_\gamma^*\le p_\tau\ \ \text{for finite } G\subset\{\gamma: r(\gamma)=\tau\}.
$$
The universal such family generates the Toeplitz algebra $\mathcal T(\Lambda)$; imposing in addition the Cuntz–Krieger relation $p_\tau=\sum_{r(\gamma)=\tau}s_\gamma s_\gamma^*$ at every face $\tau$ that receives at least one arrow gives $C^*(\Lambda)$. Both are spanned by $s_\mu s_\nu^*$ with $s(\mu)=s(\nu)$, with the product rule
$$
(s_\mu s_\nu^*)(s_\alpha s_\beta^*)=\begin{cases}s_{\mu\alpha'}s_\beta^* & \alpha=\nu\alpha'\\ s_\mu s_{\beta\nu'}^* & \nu=\alpha\nu'\\ 0&\text{otherwise.}\end{cases}
$$
A faithful representation lives on $\ell^2(E^*)$: $Q_\sigma h_\mu=[r(\mu)=\sigma]\,h_\mu$ and $T_\gamma h_\mu=h_{\gamma\mu}$ if $s(\gamma)=r(\mu)$, else $0$ (aHLRS §3).

Read as dynamics: $h_\mu$ is a *history*; $p_\sigma$ asks "does the history currently end at the face $\sigma$?"; $s_\gamma$ *extends the history by the arrow $\gamma$*. The arrow is now an operator between face projections: its initial projection is the source face, its range projection $s_\gamma s_\gamma^*$ ("arrived via $\gamma$") is a sub-projection of the target face. This is the precise algebraic form of "an arrow carrying its projection information".

**Derivation (structure for a finite layered quiver).** Let $N(\sigma)=\#\{\mu\in E^*: s(\mu)=\sigma\}$ be the number of forward paths starting at $\sigma$, the trivial one included. Then
$$
\mathcal T(\Lambda)\;\cong\;\bigoplus_{\sigma\in E^0}M_{N(\sigma)}(\mathbb C),\qquad
C^*(\Lambda)\;\cong\;\bigoplus_{\sigma_0\in K_0}M_{N(\sigma_0)}(\mathbb C).
$$
*Proof.* $\ell^2(E^*)=\bigoplus_\sigma\ell^2(\{\mu: s(\mu)=\sigma\})$ and every $s_\mu s_\nu^*$ preserves the summands. In the summand of $\sigma$, the element $p_\sigma-\sum_{r(\gamma)=\sigma}s_\gamma s_\gamma^*$ is the rank-one projection onto the trivial path $h_\sigma$ (numerically confirmed in §9), hence $s_\mu(p_\sigma-\sum s_\gamma s_\gamma^*)s_\nu^*=e_{\mu\nu}$ for all $\mu,\nu$ starting at $\sigma$, so each summand is full. The Cuntz–Krieger relation kills precisely the block of every face that receives an arrow, leaving the blocks of the layer-$0$ faces. $\square$

So $C^*(\Lambda)$ is "operators that mix histories with the same origin"; the groupoid behind it is $\{(\mu\lambda,\ |\mu|-|\nu|,\ \nu\lambda)\}$, **re-routings of the recent portion of a history over a shared past** $\lambda$. Its unit space is the boundary path space: histories whose start lies in layer $0$ (Exel's tight spectrum, Farthing–Muhly–Yeend's boundary paths). Vertex projections $p_\tau$ at layer $\ell$ project onto the $h_\tau=\#\{\mu: s(\mu)\in K_0,\ r(\mu)=\tau\}$ histories arriving at $\tau$. Reversing all arrows gives the mirror algebra $\bigoplus_{\tau\in K_L}M_{N^-(\tau)}$, "operators that mix histories with the same output face"; that is the Christensen–Thomsen convention.

### 3.2 Three nested resolutions

Inside $C^*(\Lambda)$ sit two commutative subalgebras:

- the **face algebra** $D_0=\mathrm{span}\{p_\sigma\}\cong C(E^0)$ — this is the barycentre level: a state restricted to $D_0$ is a point of the simplex $\bigoplus_\ell\Delta(K_\ell)$, barycentric coordinates over faces;
- the **history algebra** $D=\mathrm{span}\{s_\mu s_\mu^*\}\cong C(\partial\Lambda)$, functions on histories with full past; its projections form the Boolean algebra of cylinder sets, Exel's semilattice of idempotents.

$D_0\subset D\subset C^*(\Lambda)$ is the refinement hierarchy: face $\to$ history $\to$ history with coherences across re-routings. Everything in the earlier transcript about "support projections" lives in $D$; everything about "the final mask identifies far fewer cases than the history" is the inclusion $D_0\subsetneq D$.

### 3.3 Convolution and the source/target commutator

**Fact.** For an étale groupoid $G$ with unit space $X$, $(f*g)(\gamma)=\sum_{\alpha\beta=\gamma}f(\alpha)g(\beta)$, $f^*(\gamma)=\overline{f(\gamma^{-1})}$, and for $b\in C(X)$: $(b*f)(\gamma)=b(r\gamma)f(\gamma)$, $(f*b)(\gamma)=f(\gamma)b(s\gamma)$. Hence
$$
[b,f](\gamma)=\big(b(r\gamma)-b(s\gamma)\big)f(\gamma).
$$
This is the "before/after the arrow" distinction from the earlier session, now with $X$ = histories. In the generator picture it reads $p_\tau s_\gamma=s_\gamma=s_\gamma p_\sigma$ for $\gamma:\sigma\to\tau$, and $b\,s_\gamma-s_\gamma\,b=(b(\tau)-b(\sigma))s_\gamma$ for $b\in D_0$.

### 3.4 The Bratteli diagram and the two dimensions

Adjoin a root $v_0$ with one arrow to every layer-$0$ face. Then $(E^0,E^1)$ is a finite Bratteli diagram with incidence matrices $F_\ell(\tau,\sigma)=\#\{\gamma:\sigma\to\tau\}$ and heights $h^{(\ell+1)}=F_\ell h^{(\ell)}$ (Bezuglyi–Karpel, [1503.03360](https://arxiv.org/abs/1503.03360) §2.2). The associated AF filtration
$$
A_0\subset A_1\subset\cdots\subset A_L,\qquad A_\ell=\bigoplus_{\tau\in K_\ell}M_{h_\tau}(\mathbb C),
$$
is the sequence of "same current face" algebras, exactly Connes's $A_n$ for Penrose tilings (§1.3); the inclusion $A_\ell\subset A_{\ell+1}$ is the incidence matrix. Its ordered $K_0$ is the dimension group $\varinjlim(\mathbb Z^{K_\ell},F_\ell)$; with infinitely many layers (or a recurrent network, which is a stationary diagram) the space of complete histories modulo eventual agreement is a Penrose-type bad quotient and $A=\varinjlim A_\ell$ is a genuinely noncommutative AF algebra.

**Two dimensions, not one.** NCG attaches to a face $\tau$ the Murray–von Neumann dimension $\mathrm{Tr}(p_\tau)=h_\tau$, the *number of histories reaching it* (a density in the infinite case). This is not the simplex dimension $d(\tau)=|\tau|$. Both are available: $d$ is a function in $D_0$, $h$ is the trace of a projection. They must not be conflated; the earlier numbers "1636 cells versus 141 masks" are, in this language, $h$ versus the number of faces, and that experiment measured $h$ on a slice of concrete activation geometry, not on the abstract quiver.

### 3.5 Change of dimension is a coboundary

The dimension observable $d=\sum_{u\in V}e_u\in D_0$, $e_u=\sum_{\sigma\ni u}p_\sigma$ (the vertex projections, commuting but not orthogonal), has along an arrow the increment $d(r\gamma)-d(s\gamma)=(\delta d)(\gamma)$. The "projection information capturing the change of dimension" is therefore the coboundary of the function $d$ on faces. By Section 6 a coboundary never produces history dependence; history dependence comes only from the part of the weight cocycle that is *not* a coboundary. The zero of a ReLU at vertex $u$ is, at the face level, the projection $1-e_u$: states killed by $e_u$ never visit a face containing $u$.

---

## 4. States, projections, and conditioning along an arrow

### 4.1 What a state is, at each resolution

A state $\varphi$ on $C^*(\Lambda)=\bigoplus_{\sigma_0}M_{N(\sigma_0)}$ is a convex combination of block states; the block of $\sigma_0$ is a full matrix algebra, so its states are density matrices on histories from $\sigma_0$. Restricted to $D$: a probability on histories. Restricted to $D_0$: barycentric coordinates over faces, layer by layer. Pure states are vector states inside one block; the decomposition of a mixed state into pure ones is non-unique exactly because the block is noncommutative (the "two orthonormal pairs" remark of the earlier session, now placed). Tracial states form a Choquet simplex with one extreme point per block (per input face; per output face in the mirror algebra).

**Fact (Cauchy–Schwarz).** $|\varphi(s_\mu s_\nu^*)|^2\le\varphi(s_\mu s_\mu^*)\varphi(s_\nu s_\nu^*)$: a history of zero mass has zero coherence with every other history. The "zero kills the row and column" observation is this inequality in the arrow algebra.

### 4.2 Conditioning by a partial isometry

For a partial isometry $v$ and a state $\varphi$ with $\varphi(v^*v)>0$,
$$
\varphi^{v}(a)=\frac{\varphi(v^*av)}{\varphi(v^*v)}
$$
is a state supported under $vv^*$. This is the *conditional arrow algebra* operation and it has two readings on the same formula:

- $v=s_\gamma$: $\varphi^{s_\gamma}$ is $\varphi$ **transported forward along $\gamma$**, defined whenever $\varphi(p_{s(\gamma)})>0$, and supported under "arrived via $\gamma$". Composition: $\varphi^{s_\mu}=(\varphi^{s_{\mu_n}})^{s_{\mu_{n-1}}\cdots}$.
- $v=s_\gamma^*$: $\varphi^{s_\gamma^*}$ is $\varphi$ **conditioned on the event "the last step was $\gamma$" and moved back** to the source face; this is Bayes' rule in the algebra.

**Derivation (law of total probability).** For a face $\tau$ not in layer $0$, the Cuntz–Krieger relation gives
$$
\varphi(p_\tau)=\sum_{r(\gamma)=\tau}\varphi(s_\gamma s_\gamma^*),\qquad
\varphi^{s_\gamma}(s_\gamma s_\gamma^*)=1,
$$
so the mass at a face is the sum over last arrows, and forward transport along a *fixed* arrow lands with certainty — the choice *between* arrows is not in $\varphi^{s_\gamma}$ but in the numbers $\varphi(s_\gamma s_\gamma^*)$, i.e. in $\varphi|_D$. No stationary walk is assumed anywhere.

### 4.3 Murray–von Neumann comparison is the "dimension-changing arrow"

$s_\gamma^*s_\gamma=p_\sigma$ and $s_\gamma s_\gamma^*\le p_\tau$ say $p_\sigma\precsim p_\tau$: the source face projection is Murray–von Neumann subequivalent to the target one, and $p_\tau$ is *partitioned* by its incoming arrows. Under any tracial state $\tau_{\mathrm{tr}}$ this gives $\tau_{\mathrm{tr}}(p_\tau)=\sum_{r(\gamma)=\tau}\tau_{\mathrm{tr}}(p_{s(\gamma)})$, which is the height recursion $h^{(\ell+1)}=F_\ell h^{(\ell)}$. The dimension *grows* along arrows in this sense even when $|\tau|<|\sigma|$: NCG dimension counts histories, not vertices (§3.4).

---

## 5. Where the weights go: a cocycle, a dynamics, and the derived walk

### 5.1 Generalised gauge action

Give each arrow a real number $F(\gamma)$ (the "action" of the step; one natural candidate is a log-likelihood or a log-gain extracted from the restricted linear map, but this note does not fix it). **Fact** (Christensen–Thomsen §2; aHLRS §1): the universal property gives a one-parameter group $\alpha^F_t(p_\sigma)=p_\sigma$, $\alpha^F_t(s_\gamma)=e^{itF(\gamma)}s_\gamma$, hence $\alpha^F_t(s_\mu s_\nu^*)=e^{it(F(\mu)-F(\nu))}s_\mu s_\nu^*$ with $F(\mu)=\sum_i F(\mu_i)$; on the groupoid it is the cocycle $c_F(\mu\lambda,\,k,\,\nu\lambda)=F(\mu)-F(\nu)$. The plain gauge action is $F\equiv1$: time is layer depth.

### 5.2 KMS states are Gibbs measures on histories

**Fact** (Renault; Neshveyev Thm 1.3): KMS$_\beta$ states of $C^*(G)$ for the dynamics of a cocycle $c$ correspond to probability measures on the unit space that are quasi-invariant with Radon–Nikodym cocycle $e^{-\beta c}$ (plus traces on isotropy, trivial here since the quiver is acyclic), composed with the conditional expectation onto $C(G^{(0)})$. **Fact** (Christensen–Thomsen Lemma 2.1, Carlsen–Larsen): on the graph algebra these are $\omega_\psi(S_\mu S_\nu^*)=\delta_{\mu\nu}e^{-\beta F(\mu)}\psi_{r(\mu)}$ with $\psi$ a normalised almost-harmonic vector for the weighted matrix $A(\beta)_{vw}=\sum_{e: v\to w}e^{-\beta F(e)}$ (their conventions reverse arrows; for the plain gauge action aHLRS Thm 3.1 and Cor. 6.1, with $\rho(A)=0$ for an acyclic graph so every $\beta\in\mathbb R$ occurs).

**Derivation (complete-history corner).** Cut $C^*(\Lambda)$ down by the projection $P_L=\sum_{\tau\in K_L}p_\tau$ onto histories that have reached the output layer. The corner is $\bigoplus_{\sigma_0}M_{N_L(\sigma_0)}$ with matrix units $e_{\mu\nu}$ indexed by complete histories, and $\alpha^F_t=\mathrm{Ad}(u_t)$ with $u_t=\mathrm{diag}(e^{itF(\mu)})$. On a full matrix algebra the KMS$_\beta$ state of $\mathrm{Ad}(e^{itH})$ is unique and equals $e^{-\beta H}/\mathrm{Tr}\,e^{-\beta H}$. Hence the KMS$_\beta$ states of the corner are exactly
$$
\mathbb P_\beta(\mu)=\frac{w_{s(\mu)}}{Z_{s(\mu)}(\beta)}\,e^{-\beta F(\mu)},\qquad Z_{\sigma_0}(\beta)=\sum_{\mu\ \text{complete},\ s(\mu)=\sigma_0}e^{-\beta F(\mu)},
$$
Gibbs measures on complete histories with energy the accumulated arrow action and a free weight $w$ on the origin (on the output face in the mirror algebra; on nothing in the pair-groupoid algebra of all histories). The KMS condition was verified numerically to $5\cdot10^{-17}$ (§9).

Two readings of Connes's "noncommutative measure spaces evolve with time":

- the modular flow of $\mathbb P_\beta$ is $\sigma^{\mathbb P_\beta}_t=\alpha^F_{-\beta t}$: **the canonical time of the state is (weighted) layer depth**;
- Connes–Rovelli's thermal-time hypothesis ([gr-qc/9406019](https://arxiv.org/abs/gr-qc/9406019)) read backwards: here the time flow is given (the layers), and the states for which that flow *is* their modular flow are precisely the Gibbs measures above. A state on histories that is not of this form has a modular flow that is not depth; the Connes cocycle $[D\varphi:D\mathbb P_\beta]_t$ measures the discrepancy.

### 5.3 The walk is derived, not assumed

Write $Z_\sigma=\sum_{\nu:\ \sigma\to\text{output}}e^{-\beta F(\nu)}$ for the backward partition function. Then $\mathbb P_\beta$ is the Markov chain on faces with kernel
$$
\mathbb P(\gamma\mid\text{at }\sigma)=e^{-\beta F(\gamma)}\frac{Z_{r(\gamma)}}{Z_\sigma},
$$
the Doob $h$-transform of $e^{-\beta F}$. The backward quantity $Z$ is computed by the **transfer operator** $(\mathcal L f)(\sigma)=\sum_{\gamma: s(\gamma)=\sigma}e^{-\beta F(\gamma)}f(r\gamma)$ applied to $f\equiv1$ on the output layer: a *question* about the output propagated backward through the layers; its dual pushes *states* forward. This is Exel's endomorphism/transfer-operator pair ([math/0012084](https://arxiv.org/abs/math/0012084)), symmetrised in his *interactions* $(\mathcal V,\mathcal H)$ with $\mathcal V\mathcal H\mathcal V=\mathcal V$, $\mathcal H\mathcal V\mathcal H=\mathcal H$ ([math/0409267](https://arxiv.org/abs/math/0409267)), whose introduction says it exactly: the endomorphism "accounts for the future", the transfer operator "for the past", and when both directions are uncertain one needs the symmetric pair. The C\*-correspondence of Muhly–Tomforde packages the same thing with the measures $\lambda_\sigma$ on outgoing arrows as the inner product. So "states go forward, questions go backward, a walk emerges" is not a hope; it is the Heisenberg–Schrödinger duality of a completely positive map, and the walk is the normalised transfer operator.

The input distribution (Gaussian or otherwise) enters as a state on $D_0$ restricted to layer $0$, i.e. the origin weights $w_{\sigma_0}$, not as a KMS condition. Whether the *actual* data law on histories is close to some $\mathbb P_\beta$ is an empirical question this note does not touch.

---

## 6. Sufficiency: when is the barycentre enough?

**Fact** (Petz; Jenčová–Petz [math-ph/0412093](https://arxiv.org/abs/math-ph/0412093), Thm 1): for a family $\{\varphi_\theta\}$ of normal states dominated by a faithful $\omega$ on a von Neumann algebra $M$ and a subalgebra $M_0\subset M$, the following are equivalent: $M_0$ is sufficient (a coarse-graining into $M_0$ preserves every $\varphi_\theta$); the Connes cocycles $[D\varphi_\theta:D\omega]_t$ lie in $M_0$ for all $t$; the $\omega$-preserving generalised conditional expectation onto $M_0$ leaves every $\varphi_\theta$ invariant; the relative entropies $S(\varphi_\theta,\omega)$ are unchanged by restriction. In the commutative case this is the classical factorisation criterion (their Prop. 3), and for exponential families $\varphi_\theta=[\omega^{\sum\theta_ia_i}]$ sufficiency of $M_0$ is equivalent to $\sigma^\omega_t(a_i)\in M_0$ (Thm 6). The minimal sufficient subalgebra is the one generated by the cocycles.

Apply this to the Gibbs family $\{\mathbb P_\beta\}_{\beta\in\mathbb R}$ of §5.2, which is the exponential family generated by the single observable $F\in D$ (the accumulated action), with the subalgebra $D_0$ of face functions.

**Derivation.** Fix the origin weights. For complete histories $\mu,\mu'$ with the same endpoint, the factorisation $d\mathbb P_\beta/d\mathbb P_0=g_\beta(r(\mu))\,h(\mu)$ for all $\beta$ forces $e^{-\beta(F(\mu)-F(\mu'))}=h(\mu)/h(\mu')$ for all $\beta$, hence $F(\mu)=F(\mu')$. Therefore:

- the **output face** is a sufficient statistic for the family iff $F(\mu)$ depends only on $r(\mu)$ on complete-past histories, i.e. $F=\delta U$ is a coboundary, $F(\gamma)=U(r\gamma)-U(s\gamma)$, with the potential $U$ constant on the input layer;
- the pair **(input face, output face)** is sufficient iff $F$ is a coboundary with arbitrary $U$;
- otherwise the minimal sufficient statistic is $F(\mu)$ itself, the accumulated action, and the face is *not* sufficient: the history matters.

Three further equivalences make the statement structural rather than statistical. When $F=\delta U$: $\mathbb P_\beta(\mu)\propto e^{-\beta U(r\mu)}$, so all histories arriving at a face are equally likely — this is exactly **tail-invariance** of the measure (Bezuglyi–Karpel §5.2: an $\mathcal R$-invariant measure gives equal mass to any two finite paths with the same endpoint), which is the condition for a **trace** on the endpoint AF algebra $A_L$ (central measure, Vershik–Kerov), which in turn is the case where the modular flow of the state on $A_L$ is **trivial**. Summarising:

$$
\text{face sufficient}\iff F\sim 0\text{ in }H^1(\Lambda)\iff \mathbb P_\beta\text{ central}\iff \mathbb P_\beta\text{ tracial on }A_L\iff\text{modular flow trivial.}
$$

Connes's unique trace on the Penrose algebra is the central measure: under it every past leading to the current symbol is equally likely, which is the history-independence that sufficiency of the face expresses. A non-coboundary cocycle is the algebraic signature of history dependence, and "noncommutative measure spaces evolve with time" is the statement that the barycentre is not sufficient. Note also §3.5: the dimension change $\delta d$ is a coboundary by construction, so face-dimension bookkeeping alone can never make the barycentre insufficient; only the weight-derived part of $F$ can.

The numerical check (§9) confirms all three cases: random $F$ — neither statistic sufficient; $F=\delta U$ with $U$ arbitrary — the pair sufficient, the face alone not; $F=\delta U$ with $U=0$ on layer $0$ — the face sufficient.

---

## 7. Finite resolution: where the operator systems are

### 7.1 Depth truncation

Let $E_k=\mathrm{span}\{s_\mu s_\nu^*:\ |\mu|,|\nu|\le k,\ s(\mu)=s(\nu)\}$, functions on re-routings that only touch the last $k$ layers, taken inside the Toeplitz algebra $\mathcal T(\Lambda)$, where the elements $s_\mu s_\nu^*$ form a basis (their number is $\sum_\sigma N(\sigma)^2=\dim\mathcal T(\Lambda)$). $E_k$ contains the unit (the sum of the face projections) and is self-adjoint, so it is an operator system; $E_m=\mathcal T(\Lambda)$ iff $m\ge L$.

**Derivation.** In $\mathcal T(\Lambda)$: $C^*_{\mathrm{env}}(E_k)=\mathcal T(\Lambda)$ and $\mathrm{prop}(E_k)=2\lceil L/k\rceil-1$; in the Cuntz–Krieger quotient the propagation number of the image of $E_k$ is at most this, since products pass to the quotient. *Proof.* $\mathcal T(\Lambda)$ is a finite direct sum of simple blocks, so a boundary ideal is a sum of blocks; the quotient by any block sends that block's vertex projections, which lie in $E_k$, to $0$, so it is not completely isometric on $E_k$; hence the Šilov ideal is $0$ and the envelope is $\mathcal T(\Lambda)$. For the propagation number, products of basis elements are basis elements or $0$, so it suffices to track the pair of lengths $(|\mu|,|\nu|)$ of a product $s_\mu s_\nu^*$ as one more factor $s_\alpha s_\beta^*\in E_k$ is multiplied on the right. By the product rule either $|\alpha|\le|\nu|$, and the pair becomes $(|\mu|,\ |\nu|-|\alpha|+|\beta|)$, or $|\alpha|>|\nu|$, and it becomes $(|\mu|+|\alpha|-|\nu|,\ |\beta|)$. So a factor can lengthen $\mu$ by at most $k$ only at the price of resetting $|\nu|$ to at most $k$, and can lengthen $\nu$ by at most $k$ only while leaving $\mu$ alone. To reach $|\mu|=|\nu|=L$ one must therefore first build $\mu$ ($\lceil L/k\rceil$ factors, the last of which may leave $|\nu|\le k$) and then build $\nu$ from at most $k$ to $L$ ($\lceil L/k\rceil-1$ further factors), and this is achieved by the word $s_{\alpha_1}\cdots s_{\alpha_{m-1}}(s_{\alpha_m}s_{\beta_m}^*)s_{\beta_{m+1}}^*\cdots$. $\square$ For $k=L$ this is $1$; for $k=1$ it is $2L-1$; the formula is confirmed combinatorially in §9. It is the depth analogue of the band-matrix count $\lceil p/N\rceil$ of Connes–van Suijlekom, with the factor $2$ because a re-routing has a past side and a future side that must be built one after the other: the number of compositions needed before a $k$-layer description closes into an algebra counts the $k$-layer blocks of a history and of its re-routing.

### 7.2 Overlap tolerance between faces

On a single complex $K_\ell$ the relation $\sigma\,R_k\,\tau\iff|\sigma\cap\tau|\ge k$ (or adjacency in the face poset) is reflexive and symmetric but not transitive: a **tolerance relation** on faces. **Facts** (§1.4): $E(R_k)\subset M_{|K_\ell|}$ is an operator system; it is an algebra iff $R_k$ is an equivalence relation; if the overlap graph is connected its C\*-envelope is the full matrix algebra and its propagation number is the diameter of the overlap graph; its pure states are the vector states whose support is $R_k$-connected. So "we only resolve faces up to sharing $k$ vertices" is a bona fide Connes–van Suijlekom finite-resolution space, with the diameter of the overlap graph as its intrinsic resolution scale. (The same object is the "non-commutative graph" of quantum graph theory.)

### 7.3 Projections become effects

Under a spectral or depth cut-off $P$, the face projections $p_\sigma$ become $E_\sigma=Pp_\sigma P$ with $0\le E_\sigma\le1$, $\sum E_\sigma=P$, and the defect $E_\sigma-E_\sigma^2=Pp_\sigma(1-P)p_\sigma P\ge0$; two sharp commuting face tests acquire the commutator $[E_\sigma,E_\tau]=P(p_\tau(1-P)p_\sigma-p_\sigma(1-P)p_\tau)P$. These are the identities from the earlier session; their home is the operator system $PC^*(\Lambda)P$, and the state space of that system is a noncommutative convex set in the Kennedy–Shamovich sense, not a simplex.

---

## 8. The optional metric layer

A finite quiver carries natural spectral triples, e.g. $D=\sum_\gamma\ell(\gamma)^{-1}(s_\gamma+s_\gamma^*)$ on $\ell^2$ of histories, for which $[D,b]$ is the weighted difference of $b$ across arrows (§3.3), so $\|[D,b]\|\le1$ is a Lipschitz constraint and Connes's formula gives a distance between barycentric states. **Fact** (Iochum–Krajewski–Martinetti, [hep-th/9912217](https://arxiv.org/abs/hep-th/9912217)): on finite commutative spaces this distance is finite iff the points are connected, is bounded by the path length, is in general *not* the path metric (explicit three-point formula; the four-point case is not solvable by radicals), and every finite metric can be realised by some triple satisfying all the axioms. Connes–van Suijlekom extend the distance to operator-system triples $(E,H,D)$ (Definition 3.1 of the truncation paper). This layer is optional: Sections 2–7 use no metric, which answers the worry that Choquet theory smuggles in a distance.

---

## 9. Numerical sanity checks

The finite-dimensional claims were checked on a random layered quiver with face counts $(3,4,3,2)$ and random arrow multiplicities in $\{0,1,2\}$ (script: `notes/checks/graph_algebra_checks.py`, seed 0, `--full` for the rank computations):

| claim | computed | predicted |
|---|---|---|
| Toeplitz–Cuntz–Krieger relations on $\ell^2(E^*)$ | hold | — |
| finite paths $\lvert E^*\rvert$ | 214 | — |
| $\dim\mathcal T(\Lambda)$ | 9390 | $\sum_\sigma N(\sigma)^2=9390$ |
| $\dim C^*(\Lambda)$ | 8548 | $\sum_{\sigma_0\in K_0}N(\sigma_0)^2=8548$ |
| $p_\tau-\sum_{r(\gamma)=\tau}s_\gamma s_\gamma^*$ | projection onto the trivial path at $\tau$ | same |
| complete histories by endpoint / origin | 47, 48 / 40, 12, 43 | heights $h_\tau$, $N_L(\sigma_0)$ |
| KMS$_\beta$ condition for the origin-block Gibbs measure, random $F$, $\beta=0.7$ | max violation $5.6\cdot10^{-17}$ | 0 |
| face sufficient / (origin, face) sufficient, random $F$ | no / no | no / no |
| same, $F=\delta U$, $U$ arbitrary | no / yes | no / yes |
| same, $F=\delta U$, $U=0$ on layer 0 | yes / yes | yes / yes |
| $\mathrm{prop}(E_k)$ in $\mathcal T(\Lambda)$, $L=3$, $k=1,2,3$ | 5, 3, 1 | $2\lceil L/k\rceil-1=5,3,1$ |

These are checks of algebra, not experiments about any network.

---

## 10. What is established, what is not, and what to do next

Established here (given a quiver $\Lambda$ and a cocycle $F$): the structure of $\mathcal T(\Lambda)$ and $C^*(\Lambda)$; the identification of face / history / coherent resolutions; the conditioning calculus; the Gibbs form of the KMS states on complete histories and the identification of their modular flow with weighted depth; the derived Doob walk; the coboundary criterion for sufficiency of the barycentre and its equivalence with centrality, traciality and trivial modular flow; the propagation number of depth truncations.

Not established: how $\Lambda_\ell$ and $F$ are to be read off the actual weights (the one modelling decision, §2.1); whether any real network's law on histories is near a Gibbs state; any claim about the activation geometry of a particular network. For a finite network the NCG *topology* is trivial (each block is Morita equivalent to a point; $C^*(\Lambda)$ is Morita equivalent to $C(K_0)$), so the content is in the *measure theory*: dimensions of face projections, the cocycle class of $F$, states and their modular flows — precisely the level at which Connes says the new phenomena live.

Next steps, in order of leverage:

1. Fix a rule $W_{\ell+1}\mapsto(\Lambda_\ell,F|_{\Lambda_\ell})$ and compute, for a real network's quiver, the heights $h_\tau$ and the dimension group; compare $h$ with the face counts.
2. Test the coboundary criterion: compute the cohomology class of $F$ (the function $F(\mu)$ on histories to a fixed face) — its spread across histories is the exact amount by which the barycentre fails to be sufficient.
3. Compute $\mathrm{prop}(E_k)$ against observed depth of history dependence, and the overlap-tolerance propagation numbers per layer.
4. Only then put a Dirac operator on the quiver and compare the Connes distance between barycentric states with transport along the derived walk.

---

## Sources read for this note

- A. Connes, *Noncommutative Geometry*, Academic Press 1994: Introduction §§1, 2, 5; Chapter II §3 (Penrose tilings). [PDF](https://alainconnes.org/wp-content/uploads/book94bigpdf.pdf)
- A. Connes, *Noncommutative geometry, Year 2000*, [math/0011193](https://arxiv.org/abs/math/0011193); *A short survey of noncommutative geometry*, [hep-th/0003006](https://arxiv.org/abs/hep-th/0003006); *Noncommutative geometry, the spectral standpoint*, [1910.10407](https://arxiv.org/abs/1910.10407).
- A. Connes, W. van Suijlekom, *Spectral truncations in noncommutative geometry and operator systems*, [2004.14115](https://arxiv.org/abs/2004.14115); *Tolerance relations and operator systems*, [2111.02903](https://arxiv.org/abs/2111.02903).
- A. Connes, C. Rovelli, *Von Neumann algebra automorphisms and time-thermodynamics relation*, [gr-qc/9406019](https://arxiv.org/abs/gr-qc/9406019).
- R. Exel, *Inverse semigroups and combinatorial C\*-algebras*, [math/0703182](https://arxiv.org/abs/math/0703182); *A new look at the crossed-product of a C\*-algebra by an endomorphism*, [math/0012084](https://arxiv.org/abs/math/0012084); *Interactions*, [math/0409267](https://arxiv.org/abs/math/0409267).
- P. Muhly, M. Tomforde, *Topological quivers*, [math/0312109](https://arxiv.org/abs/math/0312109).
- A. an Huef, M. Laca, I. Raeburn, A. Sims, *KMS states on the C\*-algebras of finite graphs*, [1205.2194](https://arxiv.org/abs/1205.2194).
- J. Christensen, K. Thomsen, *Finite digraphs and KMS states*, [1505.04751](https://arxiv.org/abs/1505.04751).
- S. Neshveyev, *KMS states on the C\*-algebras of non-principal groupoids*, [1106.5912](https://arxiv.org/abs/1106.5912).
- S. Bezuglyi, O. Karpel, *Bratteli diagrams: structure, measures, dynamics*, [1503.03360](https://arxiv.org/abs/1503.03360).
- A. Jenčová, D. Petz, *Sufficiency in quantum statistical inference*, [math-ph/0412093](https://arxiv.org/abs/math-ph/0412093).
- B. Iochum, T. Krajewski, P. Martinetti, *Distances in finite spaces from noncommutative geometry*, [hep-th/9912217](https://arxiv.org/abs/hep-th/9912217).
- M. Kennedy, E. Shamovich, *Noncommutative Choquet simplices*, [1911.01023](https://arxiv.org/abs/1911.01023); K. Davidson, M. Kennedy, *Noncommutative Choquet theory: a survey*, [2412.09455](https://arxiv.org/abs/2412.09455) (for Kadison's representation theorem and Bauer's theorem).
- Also consulted earlier in this thread: L. Helmer, B. Solel, *Weighted Cuntz–Krieger algebras*, [2108.05601](https://arxiv.org/abs/2108.05601) (weights as operators on the Fock correspondence — a fourth possible role for the weights, not used above).
