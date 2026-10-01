# Hammersley–Clifford: Markov = Gibbs, its Möbius-inversion proof, why positivity is needed, and its approximate and quantum versions

### Bridges digest: the classical theorem lemma by lemma; what fails without positivity (support shape, intersection, local cohomology); quantum conditional independence as Petz sufficiency; quantum Hammersley–Clifford and where it fails; recovery maps; "approximate quantum Markov networks are thermal"; decay of conditional mutual information. Then a test of which of this touches expanders, local-to-global and NCG.

*Prong-3 reading note (see [research-program.md](../../research-program.md)), 2026-10-01. Inputs used from the other prongs, as vocabulary only: Note 1 ([conditional-arrow-algebra.md](../../conditional-arrow-algebra.md): KMS/Gibbs states on histories §5.2, the coboundary criterion §6, Jenčová–Petz sufficiency), [local-to-global-unlocks.md](../../local-to-global-unlocks.md) (the four-ingredient table §4 and U3), [mlp-bridge.md](../../mlp-bridge.md) §3.2 (25–35 % memory of the face process). Sibling digests in this folder that this file builds on and does not repeat: [arxiv-2609.38007.md](arxiv-2609.38007.md) (Yang; the Chen–Rouzé Lindbladian in full), [expanders.md](expanders.md), [hdx-spectral-independence.md](hdx-spectral-independence.md). Numerical checks: [check_hammersley_clifford.py](check_hammersley_clifford.py) (§15).*

**Labels.** **Source** = stated in a text retrieved in this session (quotations in "…" are verbatim from the extraction; displayed formulas are my LaTeX transcription of it). **Mine** = my own derivation, reading or name for something. Bridges carry **THEOREM** (proved by a source, or by an elementary argument given here and checked numerically), **KNOWN-LINK** (a link stated in a source, or a standard fact; standard facts not re-read in this session are marked *from memory*), **ANALOGY**, **SPECULATION**.

**Sigla.** [HC] Hammersley–Clifford 1971 (not retrieved, §0). [Cl90] P. Clifford, "Markov random fields in statistics", in G. Grimmett, D. Welsh (eds.), *Disorder in Physical Systems* (OUP 1990) 19–32; its §2 restates the main results and the method of [HC]. [Be] Besag, JRSS B 36 (1974) 192–236, with the printed discussion, which includes written contributions by Hammersley and by Clifford. [Gr] Grimmett, Bull. LMS 5 (1973) 81–84. [Pr] Preston, Adv. Appl. Prob. 5 (1973). [Mo] Moussouris, J. Stat. Phys. 10 (1974). [Po] D. Pollard, Yale handout "Hammersley–Clifford theorem for Markov random fields" (2004). [Sh] D. Shah, MIT 6.438 Lecture 3 (2014). [Yu] H. Yu, "Markov Random Fields", Helsinki lecture slides (2010). [Du] Duke STA 345 notes "Introduction to Graphical Models" (2010). [GL] Gandolfi–Lenarda, MEMOCS 4 (2016) 407–. [CM] Chandgotia–Meyerovitch, arXiv:1305.0808v3. [Ch] Chandgotia, arXiv:1406.1849v2. [Pe] Peters, arXiv:1403.0408 (Exa highlights only). [HJPW] Hayden–Jozsa–Petz–Winter, quant-ph/0304007. [LP] Leifer–Poulin, arXiv:0708.1337. [BP] Brown–Poulin, arXiv:1206.0755. [LZ] Lauritzen–Zwiernik, arXiv:2605.19453v1. [FR] Fawzi–Renner, arXiv:1410.0664. [JRSWW] Junge–Renner–Sutter–Wilde–Winter, arXiv:1509.07127. [ILW] Ibinson–Linden–Winter, quant-ph/0611057. [KB] Kato–Brandão, arXiv:1609.06636. [KKB] Kuwahara–Kato–Brandão, arXiv:1910.09425. [Kuw] Kuwahara, arXiv:2407.05835v3. [CR] Chen–Rouzé, arXiv:2504.02208. [BLMT] Bakshi–Liu–Moitra–Tang, arXiv:2510.08542. [C26] Chen, "Note on strong quantum Markov properties", arXiv:2605.02877. [Y] Yang, arXiv:2609.38007 (through the sibling digest).

---

## 0. Retrieval status

**The manuscript itself was not retrieved.** Every route failed, and failed again when the review re-tried them independently:
- `www.statslab.cam.ac.uk/~grg/books/hammfest/hamm-cliff.pdf` (the copy linked from Grimmett's page for *Disorder in Physical Systems*, which lists it as "Bonus paper: Markov fields on finite graphs and lattices, by J. M. Hammersley and P. Clifford, unpublished manuscript, 1971"): `curl` → proxy CONNECT 403; WebFetch → `EGRESS_BLOCKED`; alphaXiv → "Failed to fetch"; Exa → `CRAWL_LIVECRAWL_TIMEOUT`. The Wayback Machine copy: `curl` → 403; Exa → an empty page.
- The Oxford Research Archive record (`ora.ox.ac.uk/objects/uuid:4ea849da-1511-4578-bb88-6a8d02f457a6`): `curl` → 403. Exa returns the catalogue page only. It reads: "This is a PDF file obtained by scanning the original unpublished 1971 typescript that has been lost. The PDF file cannot be read as text because the scan is of low precision and poorly aligned. The paper is frequently cited but no copy is currently available." It also says the file is "Available personally from" Clifford.
- The Semantic Scholar page given in the brief: the alphaXiv reader returned a schema error.

**What stands in for it (added in review).** Two texts by the authors themselves were retrieved:
- **[Cl90] §2 (pp. 21–26).** Full text via Exa, from a PDF of the whole volume (`lib.ysu.am/disciplines_bk/74370a8b08cd84bd0536f26501c88858.pdf`). Clifford writes that the 1971 paper "was never published and only a few copies were distributed … there is, perhaps inevitably, some confusion about the exact contents. The method of proof in the unpublished paper is constructive and the operator techniques used are unusual. For these reasons it seems appropriate to take this opportunity to state the main results and to describe the methods by which they were obtained. This is done in Section 2." §2.0 and §3.8 below are taken from it.
- **Hammersley's and Clifford's written contributions to the discussion of [Be]** (Clifford p. 228; Hammersley pp. 230–231). Read with the alphaXiv reader on a second, text-bearing scan (`cise.ufl.edu/~anand/fa11/Besag_Spatial_interaction.pdf`). They give the authors' own account of why [HC] stayed unpublished (§2.4).

So **nothing below is quoted from [HC] itself**. The content of [HC] is quoted from [Cl90] §2, its authors' restatement, and from [Be] §3, which says "Our definitions will closely follow those of Hammersley and Clifford". Later formulations come from [CM] Theorems 3.1–3.2, [BP] Theorem 1, [LP] Theorem 4.1, [GL] Theorem 3.1 and the abstract of [Gr].

**Classical sources.**
- [Be]: the JSTOR scan hosted at `stat.cmu.edu`. The alphaXiv reader returned 46 empty pages (image PDF). Exa returned an OCR of pp. 192–195 (Summary, §1, §2 up to eq. (2.3)) and stopped there. *Added in review:* a second scan (`cise.ufl.edu`) gave the alphaXiv reader the text of pp. 192–201, 203 and the discussion, pp. 228–235. Most displayed equations are missing from that extraction. Besag's §3 statement and proof are now quoted from it (§3.6). The displayed expansion of $Q(x)$ is still taken from the [Du] notes, which say they follow it.
- [Gr]: abstract fragment and reference list only (Exa library page). The full text (LMS/Wiley; Grimmett's site blocked) was not read.
- [Pr]: abstract only (Cambridge Core page via Exa).
- [Mo]: abstract only (ADS via Exa). The counterexample is taken from [GL] Example 3.3 (full text of pp. 407–417 via Exa) and from [Yu], which cites Lauritzen (1996).
- Spitzer 1971, Averintsev 1970, Sherman 1973, Dobrushin 1968, Lauritzen's book (1996), Geiger–Meek–Sturmfels 2006: not retrieved; they are cited only as other sources cite them. (Clifford's 1990 survey [Cl90] was retrieved in review.)
- [Po], [Sh], [Yu], [Du]: full or near-full text via Exa.
- [CM]: full text via alphaXiv, saved and grepped.
- [Ch]: pp. 1–2, 4–5, 9–12, 14–15, 21, 26 via alphaXiv.

**Quantum sources.**
- Full text: [BP], [HJPW], [ILW], [C26].
- Selected pages:
  - [LP]: pp. 1, 3, 5, 8, 10, 13, 15–16, 18–23, 48–49.
  - [FR]: pp. 1–4, 15–19, 22, 29–31.
  - [JRSWW]: pp. 1–12, 16–20, 24.
  - [KB]: pp. 1–11, 13, 17, 24, 29, 31.
  - [KKB]: pp. 1–4, 6, 10–11, 17, 21, 23.
  - [Kuw]: pp. 1–3, 18, 24, 85–86, 89, 91, 94.
  - [CR]: pp. 1–2, 4–5, 9–10, 13.
  - [BLMT]: pp. 1–5, 8–9, 12, 41, 43, 45–46.
  - [LZ]: pp. 1–3, 6–9, 13–14, 18, 22, 26–29.
- Exa highlights only: [Pe]; Lauritzen–Zwiernik's companion "Bayesian networks of density operators" (arXiv:2607.27876, abstract and a few theorem statements); a 2017 thesis by M. Gasse (one historical sentence, §2.4).
- Not retrieved, known only as cited:
  - Accardi–Frigerio (1983), via [HJPW] §V.
  - Petz (1986, 1988), via [HJPW] Theorem 3 and [JRSWW].
  - Poulin–Hastings (arXiv:1012.2050), via [BP] Theorem 5.
  - Bravyi–Vyalyi, via [BP].
  - Chen–Kastoryano–Gilyén (arXiv:2311.09207), via [C26] §II.
  - Bergamaschi–Chen–Vazirani (arXiv:2510.08538), via [C26] Theorem II.1.
  - Kato–Kuwahara (arXiv:2504.02235).
  - Sutter–Fawzi–Renner (arXiv:1504.07251).
  - Jouneghani et al. (2014), cited by [KKB] as a second proof of the commuting case.

**Corrections to the brief.**
1. The 2609.38007 paper is attributed in its own AI statement to "GPT-6 Astra Pro", not to Claude. See [arxiv-2609.38007.md](arxiv-2609.38007.md) §0, which also records that the ChatGPT share page could not be fetched; it was not retried here.
2. `C:\Users\User\Downloads\expander_survey.pdf` is on the user's machine and was not available.

---

## 1. The classical objects

### 1.1 Setting (Source: [Po] §1, [Yu], [GL] §2)

$G=(V,E)$ is a finite simple undirected graph with $|V|=n$, and each $v\in V$ carries a finite state set $S_v$. $\Omega=\prod_v S_v$, and $P$ is a probability on $\Omega$. Notation:
- $x_A$ is the restriction to $A\subseteq V$.
- $\partial A=\mathrm{bd}(A)$ is the set of neighbours of $A$ outside $A$, and $\mathrm{cl}(A)=A\cup\partial A$.
- A **complete set** (often loosely "clique") is a set of pairwise adjacent vertices.
- $S$ **separates** $A$ from $B$ if every path from $A$ to $B$ meets $S$ ([Yu]; [BP] says "$B$ shields $A$ from $C$").

### 1.2 Positivity: three strengths

- **(Pos-B)** [Be] §2, verbatim: "if $x_1,\dots,x_n$ can individually occur at the sites $1,\dots,n$, respectively, then they can occur together. Formally, if $P(x_i)>0$ for each $i$, then $P(x_1,\dots,x_n)>0$. This is called the positivity condition by Hammersley and Clifford (1971) and will be assumed throughout the present paper. It is usually satisfied in practice." [Yu] glosses it as "if $y_j$'s can occur singly they can occur together."
- **(Pos-S)** strict positivity, $P(\omega)>0$ for all $\omega\in\Omega$ ([Po] Def. ⟨1⟩(i); [GL] marks it by a superscript "$>$"). This is the form in [HC] itself. [Cl90] §2.1 builds "the positivity condition $P(\chi)>0,\ \forall\chi\in C$" into the definition of "Markovian", and Hammersley (discussion of [Be], p. 230): "we assumed a positivity condition, namely that no probability should be zero."
- **(Pos-safe)** [CM] §3.1, verbatim: "A topological Markov field $X\subset A^V$ is said to have a safe symbol if there exists an element $\star\in A$ such that for all $x\in X$ and $A\subset V$, $y$ given by $y_n=x_n$ for $n\in A$, $\star$ for $n\in A^c$ is also an element of $X$." [CM] §1: "The key assumption in the proof of the Hammersley–Clifford Theorem is the existence of a so-called 'safe symbol'." [Ch] §1 calls it "a positivity assumption on the MRF given by the presence of a safe symbol in the support, also referred to as the vacuum state."

*Mine.* (Pos-B) says the support is the product of the single-site supports. Shrinking each $S_v$ to that support turns it into (Pos-S). (Pos-S) gives every state the role of a safe symbol. (Pos-safe) is strictly weaker: it asks only that the support be closed under "switch any set of sites to the vacuum". §3 shows this closure is all the proof uses.

### 1.3 Four Markov properties and factorisation (Source: [Yu] slides 21–22, after Lauritzen 1996)

- **(F)** $P$ factorises: $p(x)=\prod_{C}\phi_C(x_C)$ over complete sets $C$, with $\phi_C\ge0$.
- **(G)** Global: for disjoint $A,B,S$, if $S$ separates $A$ from $B$ then $X_A\perp X_B\mid X_S$.
- **(L)** Local: for all $v$, $X_v\perp X_{V\setminus\mathrm{cl}(v)}\mid X_{\mathrm{bd}(v)}$.
- **(P)** Pairwise: for non-adjacent $v,v'$, $X_v\perp X_{v'}\mid X_{V\setminus\{v,v'\}}$.

Source [Yu]: "Fact: (F) ⇒ (G) ⇒ (L) ⇒ (P) … (All inclusions are strict generally.)" "The version of Hammersley–Clifford theorem shown earlier establishes (L) ⇔ (F) under the positivity condition … There is another version for (P) ⇔ (F) under the same condition (see Graphical Models by Lauritzen, 1996)."

Variants found in the sources:
- [LP] §4.1 defines a Markov network by the local property for **all subsets** $U$: $H(U:V\setminus(U\cup n(U))\mid n(U))=0$.
- [GL] §2 defines "pair-Markov" as (P). Their "global-Markov" says: for disjoint non-neighbouring $A,B$, $A\perp B\mid \Lambda\setminus(A\cup B)$.
- [CM] §2.1 defines a Markov random field by the same condition for finite separated sets. On infinite graphs it adds a "global" version for infinite sets.

**Graphoid axioms** (Source: [LP] (85)–(89)):
- symmetry;
- decomposition $I(U,W\cup Y|X)\Rightarrow I(U,W|X)$;
- weak union $I(U,W\cup Y|X)\Rightarrow I(U,W|X\cup Y)$;
- contraction $I(U,W|X)\wedge I(U,Y|X\cup W)\Rightarrow I(U,W\cup Y|X)$;
- for a *positive* graphoid, also **intersection** $I(U,W|X\cup Y)\wedge I(U,Y|W\cup X)\Rightarrow I(U,W\cup Y|X)$.

[LP] Theorem 4.6, quoting Lauritzen: "The undirected graph dependency model is equivalent to the dependency model obtained by setting $I(U,V-(U\cup n(U))|n(U))$ for all $U\subseteq V$ … and demanding closure under the positive graphoid axioms." [LP]: "although its closure under the positive graphoid axioms is equivalent to the Global Markov Property, this is not the case for a graphoid that doesn't satisfy intersection." [LP]: "if $P(V)$ is positive for all possible valuations of the variables then the associated dependency model is actually a positive graphoid."

### 1.4 Gibbs distributions and potentials (Source: [Po] Def. ⟨2⟩, [BP] Thm 1, [Yu])

$P$ is **Gibbs for $G$** if $P(x)=\frac1Z\exp\big(\sum_{C}\phi_C(x_C)\big)$ with real $\phi_C$ on complete sets ([BP] (7)). [Yu]: "$\{\phi_C\}$ is called a potential. Note: the potential … is not unique." Under (Pos-S), Gibbs is the same as (F) with positive factors.

---

## 2. The theorem

### 2.0 What [HC] proved, in the authors' formulation (Source: [Cl90] §2; added in review)

Notation of [Cl90] §2.1, verbatim where quoted:
- $G=(Z,E)$ is finite, with $Z=\{z_1,\dots,z_n\}$. $X+Y$ is the union, $X-Y$ the difference, and "$\partial Y=\{x:(x,y)\in E,\ x\notin Y,\ y\in Y\}$". A set is a clique "if and only if it is a singleton or if every member of $Y$ is a neighbour of every other member of $Y$".
- Each site has a finite colour set $C_i$ with more than one element, and "every set contains a colour which we can agree to call black". "Let $\chi^Y$ denote the colouring obtained from $\chi$ by changing the colours on the sites in $Y$ to black." $\chi_Y$ is the partial colouring on $Y$. (The extraction loses sub- and superscript placement; here blackening is a superscript and restriction a subscript.)
- "A set $Y$ is said to be light relative to $\chi$ if no site in $Y$ is black under the colouring $\chi$." $L_\chi$ is the set of cliques that are light relative to $\chi$.
- $P$ satisfies "the positivity condition $P(\chi)>0,\ \forall\chi\in C$".
- **Condition $M(X)$.** "We say that $P$ is Markovian for the set $X$ if and only if it satisfies the positivity condition and $P(\chi)/P(\chi_{Z-X})=P(\chi_{X+\partial X})/P(\chi_{\partial X})$." "If we postulate $M(z)$ for all singleton sets $z\in Z$, we say $P$ is locally Markovian. If we postulate $M(X)$ for all $X\subseteq Z$ we say it is globally Markovian."

**Theorem 1 of [HC]** (as [Cl90] states it). "Global and local Markov properties are equivalent."

**Theorem 2 of [HC].** "$P$ is Markovian if and only if it can be written in the form
$$P(\chi)/P(\chi^Z)=\exp\sum_{Y\in L_\chi}Q(\chi_Y),$$
where $Q$ is an arbitrary real-valued function of light colourings on cliques. Furthermore, if $P$ is Markovian then the associated function $Q$ is given by
$$Q(\chi_Y)=\sum_{X\subseteq Y}(-1)^{|X|}\log P(\chi^{(Z-Y)+X}),\qquad\forall\,Y\in L_\chi."$$

Hammersley's one-line summary (discussion of [Be], p. 230): "Essentially our theorem states that the probabilities associated with a Markov field must satisfy certain algebraic identities."

*Remarks (Mine unless marked).*
- Theorem 2's $Q$ is Lemma A's $\Phi_Y$ (§3.1) with black as the vacuum. $\chi^{(Z-Y)+X}$ leaves exactly $Y-X$ unblackened, and $(-1)^{|X|}=(-1)^{|Y\setminus(Y-X)|}$. "Light" is Lemma A(2). So the explicit Möbius formula is already in [HC], not only in [Gr]. (Mind the notational clash: in §3, $x^B$ *keeps* $B$ and puts the vacuum elsewhere, while [HC]'s $\chi^Y$ *blackens* $Y$.)
- [Cl90]'s restatement of [HC] has no pairwise property (P). Their "global" property, $M(X)$ for all $X$, is equivalent to separation (G) with no positivity. (G) ⇒ $M(X)$ because $\partial X$ separates $X$ from $Z-(X+\partial X)$. Conversely, if $S$ separates $A$ from $B$, let $X$ be the union of the components of $Z-S$ that meet $A$. Then $\partial X\subseteq S$, and $M(X)$ followed by weak union and decomposition gives $A\perp B\mid S$.
- In [HC] the graph is not given in advance (Source: [Be] p. 196, "Our definitions will closely follow those of Hammersley and Clifford"): "site $j(\ne i)$ is said to be a neighbour of site $i$ if and only if the functional form of $P(x_i\mid x_1,\dots,x_{i-1},x_{i+1},\dots,x_n)$ is dependent upon the variable $x_j$". The graph is the minimal graph of the single-site conditionals.
- *Checked in review (scratch script, not part of [check_hammersley_clifford.py](check_hammersley_clifford.py)).* Setting: $C_4$, three colours, colour 0 black. For one random positive Gibbs law and one generic positive law, $M(X)$ holds iff $\beta_X\log P=\log P$ (§3.8), for each of the 15 non-empty $X$. For the Gibbs law, $\max|Q|$ over non-cliques is $8.9\cdot10^{-15}$, and Theorem 2 reconstructs $\log P$ to $4.4\cdot10^{-15}$.

### 2.1 Finite graphs

**Theorem (Hammersley–Clifford; finite form).** Let $G$ be finite, $\Omega$ finite, and $P$ strictly positive. Then (P), (L), (G) and (F) are equivalent. Equivalently, $P$ is Markov for $G$ iff $\log P=\sum_{C\ \text{complete}}\phi_C(x_C)+\text{const}$.

This is the modern form, as the secondary sources below state it. [HC]'s own statement (§2.0) is (L) ⇔ (G) ⇔ (F) with an explicit formula for the potential; [Cl90]'s restatement contains no pairwise version; [Yu] cites Lauritzen (1996) for (P) ⇔ (F).

Sources, in their own words:
- [GL] Theorem 3.1: "Given a graph $\mathcal G=(\Lambda,\mathcal B)$, a random field $P$ is $\mathcal B$-pair-Markov$^>$ if and only if it is $\mathcal B$-$\emptyset$-Gibbs for some potential $\phi$."
- [BP] Theorem 1: "Let $G=(V,E)$ be a graph and $P(V)$ be a positive probability distribution … The pair $(P(V),G)$ is a positive Markov network if and only if … $P(V)=\frac1Ze^{H(V)}$ where $H(V)=\sum_{Q\in\mathcal C}h_Q(Q)$."
- [LP] Theorem 4.1 gives the same with "$\psi(C)$ a positive function".
- [Pr] abstract: "It is shown that the set of Markov random fields and Gibbs states with nearest neighbour potentials are the same for any finite graph. The set of Markov random fields is also shown to be the same as the equilibrium states of time-reversible birth/death processes with nearest neighbour interactions defined on the graph."

### 2.2 Uniqueness: the canonical potential

The potential is not unique, but once normalised at a vacuum it is. See §3.3 for the argument. [GL] §4, at the end of the proof of their Theorem 3.1 (just before Lemma 4.2): "The interaction thus identified is unique, except for the value of $J_{\sigma(0)}$".

### 2.3 Infinite graphs, specifications, and Markov cocycles (Source: [CM] §3)

For infinite locally finite $G$ the measure is not the primary object: the specification is. [CM]: "It may happen that two distinct Markov random fields have the same specification, as in the case of the 2-dimensional Ising model in low temperature." The objects:

- A **topological Markov field** $X\subset A^V$ is closed and such that "for all finite $F\subset V$ and $x,y\in X$ satisfying $x|_{\partial F}=y|_{\partial F}$, there exist $z\in X$" equal to $x$ on $F$ and $y$ off $F$. The support of a Markov random field is one ([CM] §2.1). On a finite graph, "$X\subset A^V$ is the support of an MRF if and only if it is a topological Markov field" ([Ch] p. 5).
- The **homoclinic relation** is $\Delta_X=\{(x,y)\in X\times X: x_n=y_n$ for all but finitely many $n\}$.
- A **$\Delta_X$-cocycle** is $M:\Delta_X\to\mathbb R$ with $M(x,z)=M(x,y)+M(y,z)$. It is a **Markov cocycle** if, whenever $x|_{F^c}=y|_{F^c}$, "the value $M(x,y)$ is determined by $x|_{F\cup\partial F}$ and $y|_{F\cup\partial F}$."
- "There is a clear bijection between Markov cocycles and Markov specifications on $X$": $M(x,y)=\log\Theta_{F,y|_{\partial F}}(y|_F)-\log\Theta_{F,x|_{\partial F}}(x|_F)$.
- [CM] §1: "Following Petersen and Schmidt we utilize the formalism of cocycles for the homoclinic equivalence relation … these are logarithms of the Radon–Nikodym derivatives with respect to the homoclinic relation, that is, logarithms of the ratio of probabilities of configurations which differ at only finitely many sites."
- The **Gibbs cocycle** of a nearest-neighbour interaction $\phi$ is $M_\phi(x,y)=\sum_{W}\phi(y|_W)-\phi(x|_W)$. [CM]: "A Borel probability measure $\mu$ is a Gibbs state with nearest neighbor interaction $\phi$ if and only if its $\Delta_X$-Radon–Nikodym cocycle is $e^{M_\phi}$". Also: "a Markov random field $\mu$ is adapted to $X$ if and only if the measure $\mu$ is non-singular with respect to $\Delta_X$."
- Notation: $\mathbf M_X$ is the vector space of Markov cocycles and $\mathbf G_X\subset\mathbf M_X$ the subspace of nearest-neighbour Gibbs cocycles. On a subshift, $\mathbf M^\sigma_X$ and $\mathbf G^\sigma_X$ are the shift-invariant versions.

**Theorem 3.1 of [CM] (weak version).** "Let $X$ be a topological Markov field with a safe symbol. Then: (1) Any Markov random field with $\mathrm{supp}(\mu)=X$ is a Gibbs state for a nearest neighbor interaction. (2) Further if $X$ is a subshift, any shift-invariant Markov random field with $\mathrm{supp}(\mu)=X$ is a Gibbs state for a shift-invariant nearest neighbor interaction." Part (2) "is not a part of the original formulation".

**Theorem 3.2 of [CM] (strong version).** "An inspection of the original proof … gives": (1) $\mathbf M_X=\mathbf G_X$; (2) $\mathbf M^\sigma_X=\mathbf G^\sigma_X$.

**Proposition 3.3 of [CM].** With a safe symbol, every MRF adapted to $X$ has full support $X$.

[CM] also note that "for a topological Markov field $X$ defined over a finite graph $G$, $\mathbf M_X$ is finite dimensional; the problem of determining which Markov cocycles are Gibbs amounts to solving a finite (but possibly large) system of linear equations. The resulting equations are essentially the 'balanced conditions' mentioned in [Moussouris]."

### 2.4 History, as the retrieved sources tell it

- Grimmett's abstract (Source): "Averintsev [1] and Spitzer [2] proved that the class of Markov fields is identical to the class of Gibbs ensembles when the domain is a finite subset of the cubic lattice and each site may be in either of two given states. Hammersley and Clifford [3] proved the same result for the more general case when the domain is the set of sites …" (truncated in the extraction). The reference list includes Rota, "On the foundations of combinatorial theory I. Theory of Möbius functions" (1964).
- Wikipedia (via Exa): "The relationship between Markov and Gibbs random fields was initiated by Roland Dobrushin and Frank Spitzer … Simpler proofs using the inclusion–exclusion principle were given independently by Geoffrey Grimmett, Preston and Sherman in 1973, with a further proof by Julian Besag in 1974."
- [Be] §2: "the celebrated Hammersley–Clifford theorem which, sadly, has remained unpublished by its authors."
- [Du]: "Besag (1974) is usually credited with the first published proof … but an alternate approach based on the Möbius Inversion Theorem was taken independently by Grimmett (1973)."
- **The authors' own account** (added in review; Hammersley, written discussion of [Be], pp. 230–231). "In proving this result, we assumed a positivity condition, namely that no probability should be zero." They disagreed with Besag that relaxing it would be "of little practical significance": "in many of the most important practical applications to statistical mechanics, the physical system is subject to constraints which prevent the system from assuming certain forbidden states". They expected a limiting argument to work ("the Principle of the Irrelevance of Algebraic Inequalities"). "So much for hand-waving mathematics. On the other hand, wriggle as we might, we were unable to convert this reasoning into a watertight argument. The very good reason for our failure was the unexpected discovery by a graduate student, Mr John Moussouris, of a counter-example!" Moussouris's thesis "supersedes the Hammersley–Clifford paper (hence our decision to leave our paper unpublished)". And: "by and large, it can now be said that the positivity condition is 'conquered'."
- Clifford, in the same discussion (p. 228): "Whatever the historical reasons for not publishing in 1971 the paper has clearly been superseded by the work of others". Besag's reply (p. 235): "I spent longer than I should care to admit, a couple of years ago, trying to overcome the condition without any form of success."
- [Cl90] §1. The work was done at Berkeley in summer 1971. Hammersley's lectures there included "Spitzer's (1971) characterisation of two-state MRFs on a square lattice. This characterisation had been obtained independently by Averintsev (1970). Hammersley and I were able to generalise the results to arbitrary graphs and lattices, and to identify the central importance of the clique functions, as terms in the potential of a generalised Gibbs distribution." Hammersley sent the paper to Besag, "who had already obtained partial results for rectangular lattices (Besag 1972). Besag then wrote to Hammersley with a much simpler, analytical proof of the general result". Also: "A simple derivation is also possible using the factorisation theorem of Brook (1964)."
- [Be] §3: the theorem "superseded the comparatively pedestrian results which had been obtained for 'nearest-neighbour' systems on the k-dimensional finite cubic lattice (Spitzer, 1971; Besag, 1972a). However, the original method of proof is circuitous and requires the development of an operational calculus (the 'blackening algebra')." And: "The positivity condition remains as yet unconquered".
- A 2017 thesis (M. Gasse; secondary) says that Hammersley and Clifford "postponed their publication in hope of relaxing" positivity. The authors' account above confirms the attempt. The reason it gives for leaving the paper unpublished is that Moussouris's counterexample superseded it.
- [Mo] abstract: "We give background for these two types, review proofs that they are in fact identical for systems with nonzero probabilities, and explore the new behavior that arises with constraints." Its keywords include "inversion formula for potentials", "barriers and wells", "strongly Markovian systems".

---

## 3. The Möbius-inversion proof, in full

Sources for the argument: [Po] Lemma ⟨4⟩ and Theorem ⟨7⟩ (general finite alphabets); [Sh] (binary, "adapted from [Grimmett]"); [GL] Lemma 4.1 (spin-product basis). *Corrected in review:* the operator form is not the digest's invention. It is [HC]'s own "blackening algebra" ([Cl90] §2.2; §3.8 below). In their notation $\varepsilon_i$ is the pure blackening operator $B_{\{i\}}$ and $\Delta_i=1-B_{\{i\}}$. Only the symbols $\varepsilon_i,\Delta_i$ and the remarks marked Mine are the digest's.

**Notation.** Fix a vacuum $o\in\Omega$. For $x\in\Omega$ and $B\subseteq V$ let $x^B$ agree with $x$ on $B$ and with $o$ off $B$. Put $g=\log P$. For $i\in V$ define on functions $f:\Omega\to\mathbb R$
$$\varepsilon_i f(x)=f(x\text{ with }x_i:=o_i),\qquad \Delta_i=\mathrm{id}-\varepsilon_i .$$
The $2n$ operators commute. $\varepsilon_i$ and $\Delta_i$ are complementary idempotents ($\varepsilon_i\Delta_i=0$). And $\sum_{A\subseteq V}\prod_{i\in A}\Delta_i\prod_{i\notin A}\varepsilon_i=\prod_i(\Delta_i+\varepsilon_i)=\mathrm{id}$ (Mine).

### 3.1 Lemma A (Möbius inversion on the Boolean lattice at a vacuum)

For $A\subseteq V$ define
$$\Phi_A(x)=\sum_{B\subseteq A}(-1)^{|A\setminus B|}\,g(x^B)\;=\;\Big(\prod_{i\in A}\Delta_i\prod_{i\notin A}\varepsilon_i\Big)g\,(x).$$
Then:
1. $\Phi_A$ depends on $x$ only through $x_A$; $\Phi_\emptyset=g(o)$.
2. If $A\ne\emptyset$ and $x_i=o_i$ for some $i\in A$, then $\Phi_A(x)=0$.
3. $g(x)=\sum_{A\subseteq V}\Phi_A(x)$; more generally $g(x^A)=\sum_{B\subseteq A}\Phi_B(x)$.

*Proof (Source: [Po]).*
1. Each $g(x^B)$ with $B\subseteq A$ ignores $x_j$ for $j\notin A$.
2. Pair each $B\not\ni i$ with $\tilde B=B\cup\{i\}$. Then $g(x^B)=g(x^{\tilde B})$ because $x_i=o_i$, and the two terms carry opposite signs.
3. The coefficient of $g(x^B)$ in $\sum_{A\subseteq T}\Phi_A$ is $\sum_{E\subseteq T\setminus B}(-1)^{|E|}$, which equals $1$ if $B=T$ and $0$ otherwise.

$\square$ In Rota's language this is the Möbius function $\mu(B,A)=(-1)^{|A\setminus B|}$ of the Boolean lattice, which [Gr] cites. In the operator form it is the expansion of the identity into the $2^n$ commuting idempotents of the Boolean algebra generated by $\{\varepsilon_i\}$. This is [HC]'s own route: [Cl90] (2.1) expands $\beta=\prod_z(B_z+B^*_z)$ in exactly this way, and the last step of their proof of Theorem 2 uses $\prod_{z\in Y}(1-B_z)=\sum_{X\subseteq Y}(-1)^{|X|}B_X$.

### 3.2 Lemma B (vanishing off complete sets; the key identity)

Let $P$ satisfy (P) and let every configuration $x^B$ ($x\in\mathrm{supp}\,P$, $B\subseteq V$) lie in $\mathrm{supp}\,P$. If $A$ contains non-adjacent $i\ne j$, then $\Phi_A\equiv0$.

*Proof (Source: [Po] Thm ⟨7⟩; [Sh] eqs. (3)–(9)).* Group the subsets of $A$ into blocks $\{B,\,B\cup i,\,B\cup j,\,B\cup ij\}$ with $B\subseteq A\setminus\{i,j\}$:
$$\Phi_A(x)=\sum_{B\subseteq A\setminus\{i,j\}}(-1)^{|A\setminus B|}\,L_B(x),\qquad L_B(x)=\log\frac{P(x^{B\cup ij})\,P(x^{B})}{P(x^{B\cup i})\,P(x^{B\cup j})}.$$
Let $z$ be the common value of these four configurations off $\{i,j\}$. Then
$$\frac{P(x^{B\cup ij})}{P(x^{B\cup j})}=\frac{P(X_i=x_i\mid X_j=x_j,X_{\rm rest}=z)}{P(X_i=o_i\mid X_j=x_j,X_{\rm rest}=z)},\qquad \frac{P(x^{B\cup i})}{P(x^{B})}=\frac{P(X_i=x_i\mid X_j=o_j,X_{\rm rest}=z)}{P(X_i=o_i\mid X_j=o_j,X_{\rm rest}=z)} .$$
(P) says the conditional law of $X_i$ given $(X_j,X_{\rm rest})$ does not depend on $X_j$, so the two ratios agree and $L_B=0$. [Sh] writes the same step as $a_0a_{ij}=a_ia_j$. $\square$

**The key estimate is an identity.** The log cross-ratio $L_B$ vanishes exactly. In operator form, $L_B(x)=\big(\Delta_i\Delta_j g\big)(x^{B\cup ij})$ (Mine).

### 3.3 Theorem and uniqueness

**Theorem.** Under the hypotheses of Lemma B,
$$P(x)=\exp\Big(\sum_{C\ \text{complete}}\Phi_C(x_C)\Big),$$
with $\Phi_C$ vanishing whenever a coordinate of $C$ is at the vacuum. The proof is Lemma A(3) plus Lemma B.

*Uniqueness (Mine; [GL] Lemma 4.2 for their basis).* Suppose $g=\sum_A\Psi_A$ with $\Psi_A$ depending on $x_A$ and vanishing at vacuum coordinates. Then $\Psi_D(x^B)=\Psi_D(x)[D\subseteq B]$. Therefore
$$\sum_{B\subseteq A}(-1)^{|A\setminus B|}g(x^B)=\sum_D\Psi_D(x)\sum_{D\subseteq B\subseteq A}(-1)^{|A\setminus B|}=\Psi_A(x),$$
so $\Psi=\Phi$. The normalised ("canonical", term *from memory*) potential is unique.

**Where each hypothesis is used.**

| hypothesis | where | what it buys |
|---|---|---|
| finiteness of $\log P$ at every $x^B$ | Lemma A | $\Phi_A$ is defined |
| $P(X_j=\cdot, X_{\rm rest}=z)>0$ at the contexts $z$ of $x^B$ | Lemma B | conditional probabilities at the mixed configurations exist |
| $P(X_i=o_i\mid\cdot)>0$ | Lemma B | division by the vacuum probability |
| (P) for non-adjacent pairs, **only** at contexts of the form $x^B$ | Lemma B | cross-ratio $=1$ |
| nothing else: no (L), no (G), no graph structure beyond the non-edges, no infinite-volume input | — | — |

All three positivity uses are exactly the safe-symbol closure (Pos-safe). Strict positivity is more than the proof needs. **Mine:** the hypothesis is a *hereditary* closure property of the support: it must contain the whole Boolean lattice $\{x^B\}_{B\subseteq V}$ below each of its points. This is the "hereditary class" row of [local-to-global-unlocks.md](../../local-to-global-unlocks.md) §4 in its most elementary form.

### 3.4 The converse (F) ⇒ (G), with no positivity

(Source: [Yu] slide 7 for the local version; [BP] proof of Thm 3 for the "enlarge $A$" step. Assembly mine.)
1. Suppose $P=\prod_C\phi_C$ with $\phi_C\ge0$, and $S$ separates $A$ from $B$.
2. Enlarge $A$ to $A'$, the union of $A$ with every component of $V\setminus S$ that misses $B$. Set $B'=V\setminus(A'\cup S)$.
3. No complete set meets both $A'$ and $B'$ (it would contain an edge across $S$). So $P(x)=f(x_{A'\cup S})\,h(x_{B'\cup S})$.
4. The factorisation criterion gives $X_{A'}\perp X_{B'}\mid X_S$, and decomposition gives $X_A\perp X_B\mid X_S$.

Zero factors are allowed throughout, so **Gibbs ⇒ Markov is exact and positivity-free.** For local (L): [Yu] shows $p(y_j\mid y_{-j})\propto\exp(-\sum_{C\ni j}\phi_C(y_C))$ depends only on $y_{\mathrm{cl}(j)}$.

### 3.5 Brook–Besag identity: conditionals determine the joint law

(Source: [Be] §2 (2.2) as transcribed by [Yu].) Under positivity, for any $y,y'$:
$$\frac{P(y)}{P(y')}=\prod_{j=1}^n\frac{p(y_j\mid y_1,\dots,y_{j-1},y'_{j+1},\dots,y'_n)}{p(y'_j\mid y_1,\dots,y_{j-1},y'_{j+1},\dots,y'_n)} .$$
[Yu]: "Under the positivity condition, $P$ is uniquely determined by its full conditional distributions … Note: $P$ here is general, not necessarily a Markov random field." [Be] §2 adds that the conditional formulation "is subject to some unobvious and highly restrictive consistency conditions … (Brook, 1964)."

*Mine.* The right-hand side multiplies single-site conditional ratios along the path $y'=z^0\to z^1\to\dots\to z^n=y$, with $z^k=(y_1..y_k,y'_{k+1}..y'_n)$. Each step is a "pivot move" in [CM]'s sense. Positivity puts every $z^k$ in the support. This is [CM] (3.2), $M(x,y)=\sum_iM(x^{(i)},x^{(i+1)})$ along a chain of pivots, and the consistency conditions are path-independence: the cocycle identity.

### 3.6 Two other proofs in the sources

- **Besag's expansion** (Source: [Du], "which we follow"). With $Q(x)=\log\{p(x)/p(0)\}$,
  $$Q(x)=\sum_ix_iG_i(x_i)+\sum_{i<j}x_ix_jG_{ij}(x_i,x_j)+\dots+x_1\cdots x_nG_{1\dots n}.$$
  Theorem ([Be] p. 198, verbatim; read in review): "for any $1\le i<j<\dots<s\le n$, the function $G_{i,j,\dots,s}$ in (3.3) may be non-null if and only if the sites $i,j,\dots,s$ form a clique. Subject to this restriction, the $G$-functions may be chosen arbitrarily."
  Proof ([Be], same page): $Q(x)-Q(\mathbf x_i)$ "can only depend upon $x_i$ itself and the values at sites which are neighbours of site $i$". If $\ell$ is not a neighbour of $1$, "Putting $x_i=0$ for $i\ne1$ or $\ell$, we immediately see that $G_{1,\ell}(x_1,x_\ell)=0$". "Similarly, by other suitable choices of $x$, it is easily seen successively that all 3-, 4-, …, $n$-variable $G$-functions involving both $x_1$ and $x_\ell$ must be null."
  Besag's set-up follows [HC]: finitely many values per site, and "the value zero is available at each site", which "ensures that, under the positivity condition, an entire realization of zeros is possible" (zero is the safe symbol). His corollary: "In the Hammersley–Clifford terminology, the local and global Markovian properties are equivalent." [GL] remark that Besag's proof "works only for the binary case and has some problematic steps", which [GL] §3 softens to "with some unclear steps in the proofs".
- **[GL] Lemma 4.1 (spin products).** With states $\Omega_x\subset\mathbb R$ and the Vandermonde matrix $V(x)=(r^s)$,
  $$P(\omega)=\frac1Z e^{\sum_{\sigma}J_\sigma\omega^\sigma},\qquad J_\sigma=\sum_\omega V^{-1}_{\sigma,\omega}\log P(\omega).$$
  Pairwise CI kills $J_\sigma$ whenever $\sigma_x\sigma_y\ne0$ for a non-edge, using $\sum_{\omega_x}V^{-1}_{\sigma_x,\omega_x}=\delta_{\sigma_x,0}$. This is the same inversion in the monomial basis instead of the vacuum-indicator basis (Mine).

### 3.7 The essence: a discrete integrability theorem (Mine)

Under positivity, (P) for the non-edge $\{i,j\}$ is equivalent to
$$\Delta_i\Delta_j\log P\equiv0 .$$
Indeed, CI is equivalent to the cross-ratio identity $P(a,b,z)P(a',b',z)=P(a,b',z)P(a',b,z)$ for all $a,a',b,b',z$, and under positivity the case $a'=o_i,b'=o_j$ implies the rest. Since $\Phi_A=\prod_{A^c}\varepsilon\,\prod_A\Delta\,g$ and the $\Delta$'s commute, $\Delta_i\Delta_jg=0$ kills every $\Phi_A$ with $A\supseteq\{i,j\}$. So:

> **Hammersley–Clifford is a discrete Poincaré lemma / separation-of-variables theorem:** a function on a product of finite sets whose mixed second differences vanish across every non-edge of $G$ is a sum of functions of complete sets.

The classical analogue for smooth functions: $\partial_i\partial_jg=0$ implies $g=g_1(x_{\hat i})+g_2(x_{\hat j})$. Probability enters only to translate "conditional independence" into "vanishing mixed difference of the log-density". Positivity enters only to make $\log P$ and the vacuum substitutions exist. **No rate, gap or expansion constant appears anywhere: the theorem is exact and combinatorial.**

### 3.8 The original proof: the blackening algebra (Source: [Cl90] §§2.2–2.3; added in review)

[Be] called [HC]'s proof "circuitous" because it "requires the development of an operational calculus (the 'blackening algebra')". [Cl90] reproduces it. Let $R$ range over real functions on the set $C$ of colourings.
1. **Pure operators.** "$B_YR(\chi)=R(\chi^Y)$." Then "$B_XB_Y=B_YB_X=B_{X+Y}$ so that pure operators commute", and "Every pure operator is a projector". Their finite linear combinations ("mixed blackening operators") form a commutative algebra with unit $1=B_\emptyset$.
2. **Lemma 1.** "If $X\subseteq Y$ then $(1-B_X)B_Y=0$."
3. **The Markov projectors.**
   - $\beta_X=B_X+B_{Z-(X+\partial X)}-B_{Z-\partial X}=1-(1-B_X)(1-B_{Z-(X+\partial X)})$.
   - $B^*_z=B_{Z-(z+\partial z)}(1-B_z)$, so that $\beta_z=B_z+B^*_z$.
   - $\beta=\prod_z\beta_z=\sum_{Y}B_{Z-Y}B^*_Y$ (their (2.1)).
4. **Lemma 2.** "If $Y\ne\emptyset$ and $Y$ is not a light clique relative to $\chi$, then $B^*_YR(\chi)=0$." A non-clique contains non-neighbours $x,y$, and Lemma 1 kills $B_{Z-(y+\partial y)}(1-B_x)$. A black site $z$ is killed by $1-B_z$.
5. **Lemma 3.** The invariant set $I(\beta)=\{R:\beta R=R\}$ "consists of those functions $R\in\mathcal R$ which have the representation $R(\chi)=S(\chi^Z)+\sum_{X\in L_\chi}S(\chi^{Z-X})$ for some $S\in\mathcal R$."
6. **Lemma 4.** "If $X\subseteq Z$, then $I(\beta)\subseteq I(\beta_X)$." The only graph input: "if $Y$ is a clique it cannot be partly in $X$ and partly in $Z-(X+\partial X)$."
7. **Lemma 5.** "The invariant set $I(\beta)$ is given by $\cap_{z\in Z}I(\beta_z)$."
8. **The Markov condition is linear in $\log P$.** $M(X)$ is equivalent to "$P(\chi)/P(\chi^X)=P(\chi^{Z-(X+\partial X)})/P(\chi^{Z-\partial X})$" (their (2.8)). With "$R(\chi)=\log P(\chi)$" this reads "$\beta_XR(\chi)=R(\chi)$": "Condition $M(X)$ is therefore equivalent to $R\in I(\beta_X)$."
9. **The proofs.**
   - Theorem 1: "if $P$ is locally Markovian then $R\in\cap_{z\in Z}I(\beta_z)=I(\beta)\subseteq I(\beta_X)$ by Lemma 4 and hence $P$ is globally Markovian."
   - Theorem 2: Lemma 3 with $Q(\chi_X):=S(\chi^{Z-X})$. Then $Q(\chi_Y)=\prod_{z\in Y}(1-B_z)R(\chi^{Z-Y})$, expanded by $\prod_{z\in Y}(1-B_z)=\sum_{X\subseteq Y}(-1)^{|X|}B_X$.

*Reading (Mine).*
- **One commutative algebra of idempotents.** The whole theorem is linear algebra there. The Markov property at $X$ says that $\log P$ is fixed by the idempotent $\beta_X$. Local ⇒ global is "a product of commuting idempotents is the idempotent onto the intersection of their ranges" (Lemma 5), plus one combinatorial fact (Lemma 4). No estimate, rate or gap enters (cf. §3.7, §13.1).
- **Why the condition is linear.** The $B_Y$ are composition operators $R\mapsto R\circ b_Y$ along the idempotent maps $b_Y:\chi\mapsto\chi^Y$. They are therefore unital $*$-endomorphisms of the commutative algebra $\mathbb R^C$, and so they commute with $\log$: $B_Y\log P=\log B_YP$. That is what makes "$P$ is Markov" a *linear* condition on $\log P$.
- **Where positivity enters.** Only so that $\log P$ and the ratios in (2.8) exist at the blackened colourings. That is the safe-symbol closure of §3.3, with black as $\star$.
- *Checked in review (scratch, $C_4$, $q=3$, not in the script).*
  - The four $\beta_z$ are idempotent and commute.
  - $\beta$ is idempotent, and $\operatorname{rank}\beta=\dim\bigcap_zI(\beta_z)=25=1+\sum_C(q-1)^{|C|}$: one parameter for the all-black value plus $(q-1)^{|C|}$ light colourings per clique $C$.
  - $\beta_X\beta=\beta$ for all 15 non-empty $X$ (Lemma 4).

*Quantum shadows (Sources [LP] App. B, [BP] Theorem 2; reading Mine).*
- **[LP] keeps the blackening algebra.** [LP]'s proof of their Theorem 4.7 is this inversion with "blacken" replaced by "compress to a pure product vacuum $|\alpha\rangle$". Their (175)–(176) set $J_U=\langle\alpha|_{V-U}H_V|\alpha\rangle_{V-U}\otimes I_{V-U}$ and $K_U=\sum_{W\subseteq U}(-1)^{|U-W|}J_W$. Their Lemma B.2, $\langle\alpha|_uK_U|\alpha\rangle_u=0$, is "light". The local Markov condition kills $K_U$ for non-cliques. [BP] take the trace as the vacuum instead (§7.2).
- **What survives and what fails.** Both quantum vacua give commuting, idempotent, unital CP maps, so the Möbius bookkeeping survives intact. Multiplicativity does not survive: a compression or a partial trace is not a $*$-homomorphism, so it does not commute with $\log$.
  - Markov ⇒ clique-local $\log\rho$ still holds ([LP] Thm 4.7, [BP] Thm 2). Its proof uses only linearity.
  - The classical step "Markov ⇔ $\log P$ is a fixed point of $\beta_X$" has no quantum counterpart that is linear in $\log\rho$ alone. Classically, (2.8) turns a statement about marginals into one about values at blackened colourings. Quantumly, the Markov condition stays linear only in the logs of the marginals (Ruskai's identity, §6.1), and those are not linear in $\log\rho$. That is why the converse direction fails (§7.5).
- *Checked in review (scratch, not in the script).*
  - [LP]'s pure-vacuum inversion of $\log\rho$ for the [BP] wheel state (random product vacuum) gives $\max\|K_U\|=4.9\cdot10^{-14}$ over non-cliques, and reconstructs $\log\rho$ to $9\cdot10^{-15}$.
  - A single-qubit vacuum compression $E$ has $\|E(\sigma^x\sigma^z)-E(\sigma^x)E(\sigma^z)\|=0.56$, so it is not multiplicative.

---

## 4. Why positivity is needed

### 4.1 What breaks

Without positivity the chain (F) ⇒ (G) ⇒ (L) ⇒ (P) still holds (§3.4, [Yu]). Each reverse arrow can fail, for a different reason.

### 4.2 (L) ⇏ (G): the intersection axiom fails

**Example (Source: [Yu] slide 17).** On the path $1-2-3-4-5$ take "$X_1=X_2$, $X_3\perp X_2$, $X_4=X_2+X_3$, $X_5=X_4$". Then every $X_j$ is independent of the rest given its neighbours, yet "$X_1,\dots,X_5$ is not a Markov chain". The local property holds and the global one fails.

**Mechanism (Source: [LP], [Pe]).** The graph closure of (L) is (G) only under intersection (§1.3), and classical CI satisfies intersection under positivity. [Pe] gives the sharp version, for densities continuous in the conditioning variables: "A necessary and sufficient condition for the intersection property is that all path-connected components of the support of the density are equivalent, that is, they can be connected by axis-parallel lines." [Pe] also states "a weaker condition than the intersection property still holds" (the "weak intersection property").

**Mine (made quantitative in §13, B3).** Suppose $A\perp B\mid CD$ and $A\perp D\mid BC$, and write $h(b,c,d)=P(A\in\cdot\mid b,c,d)$. Then $h$ is invariant under moving $b$ with $(c,d)$ fixed, and under moving $d$ with $(b,c)$ fixed, *within the support*. Intersection asks that $h$ depend on $c$ alone. That holds iff, for each $c$, the graph on $\mathrm{supp}(B,C,D)$ whose moves change $b$ or $d$ is connected. This is irreducibility of the **two-block Gibbs sampler** that alternately resamples $B$ given $(C,D)$ and $D$ given $(B,C)$. Positivity makes the sampler irreducible. A support like $\{X_A=X_B=X_D\}$ makes it reducible (§15, H8: $c_F=1$).

### 4.3 (G) ⇏ (F): Moussouris

**Example 3.3 of [GL], verbatim data.** "Take $\Lambda=\{1,2,3,4\}$, $\mathcal B_2=\{(1,2),(2,3),(3,4),(4,1)\}$ and let $P_M$ be the uniform distribution on
$$\Omega_0=\{(0,0,0,0),(1,0,0,0),(1,1,0,0),(1,1,1,0),(1,1,1,1),(0,1,1,1),(0,0,1,1),(0,0,0,1)\}."$$
[GL]: $\{1\}\perp\{3\}\mid\{2,4\}$ and $\{2\}\perp\{4\}\mid\{1,3\}$, "so that $P_M$ is $\mathcal B_2$-pair-Markov, and is (on this graph) also $\mathcal B_2$-global-Markov." Yet it "cannot be $\mathcal B_2$-$\mathcal B_2$-Gibbs".

Contradiction, as [Yu] gives it:
- $p(0,0,0,0)=\tfrac18$ forces $\phi_{12}(0,0)\phi_{23}(0,0)\phi_{34}(0,0)\phi_{41}(0,0)\neq0$.
- $p(0,0,1,0)=0$ and $p(0,0,1,1)=\tfrac18$ then force $\phi_{23}(0,1)\ne0$ and $\phi_{34}(1,0)=0$.
- This contradicts $p(1,1,1,0)=\tfrac18=\phi_{12}(1,1)\phi_{23}(1,1)\phi_{34}(1,0)\phi_{41}(0,1)$.

The 2017 thesis (Exa) puts the reason in one line: "each combination of {A,B}, {B,C}, {C,D} or {D,A} has a positive probability, while some joint combinations … have zero probabilities. The only way to obtain a zero probability … would be to set one of the clique potentials … to zero, which would immediately result in a zero probability for the corresponding pairwise combination."

*Checked (§15, H3).*
- Every law on $\Omega_0$, whatever its weights, has both pairwise CMIs equal to $0$. Each of the 8 cross-ratio relations of [GL] (5) reads $0=0$ on this support.
- Every edge pattern occurs in $\Omega_0$, so no edge factor can vanish.
- The single-flip graph of $\Omega_0$ is an 8-cycle (all degrees 2) containing **no commuting square**. It is the staircase $0000\to1000\to1100\to1110\to1111\to0111\to0011\to0001\to0000$, which flips each coordinate twice.

### 4.4 Gandolfi–Lenarda's sharpening: two obstructions (Source, with a structural reading that is Mine)

**[GL]'s reading of Moussouris.** "The probability in Moussouris example is $\mathcal B_2$-$\Lambda$-Gibbs, actually even $\emptyset$-$\Lambda$-Gibbs". It is uniform Bernoulli with one global hard-core constraint. Their Lemma 5.1: "a probability $P$ is $\mathcal B$-$\Lambda$-Gibbs if and only if it has a $\mathcal B$-global-Markov strictly positive extension." Their Example 6.2: any strictly positive law on $\Omega_0$ has such an extension.

**[GL] Lemma 5.2 (a worse example).** The support is
$$\Omega_{00}=\{(1,1,1,1),(0,1,1,1),(1,0,1,1),(0,0,1,1),(1,1,0,0),(1,0,0,0),(0,1,0,0),(0,0,0,0)\},$$
with $P^*(1,1,1,1)=\tfrac29$ and $\tfrac19$ elsewhere. It is global-Markov on $C_4$ but has **no** strictly positive Markov extension. Proof by a chain of the cross-ratio relations:
- (I) gives $\hat P(0,1,0,1)/\hat P(1,1,0,1)=\tfrac12$;
- (VII), (II) and (VIII) give $\hat P(1,1,0,1)=\hat P(1,0,0,1)=\hat P(0,0,0,1)=\hat P(0,1,0,1)$;
- these contradict each other.

[GL] Example 6.3 says $P^*$ is $\mathcal B_3$-$\mathcal B_3$-Gibbs and lists, among others, "$J_{\{1,2,3\}}=\log 2$", $J_{\{3,4,1\}}=\log\frac23$ and $J_{\{4,1\}}=\log\frac32$.

*Review note: the source is inconsistent here.* [GL] define $\mathcal B_3=\mathcal B_2\cup\{\{2,4\}\}$ in their Example 6.1. That graph's triangles are $\{1,2,4\}$ and $\{2,3,4\}$. The sets $\{1,2,3\}$ and $\{3,4,1\}$ listed in 6.3 are the triangles of $\mathcal B_2\cup\{\{1,3\}\}$. (Mine:) on $\Omega_{00}$ we have $x_3=x_4$, so $x_1x_2x_3=x_1x_2x_4$, and either chord carries $P^*$'s 3-body term.

*Checked (§15, H3).* Treat log-ratios on the support as data and the 8 cross-ratio equations as constraints on all 16 log-probabilities. The least-squares residual is $0$ for Moussouris' uniform law (an extension exists) and $0.231$ for $P^*$ (none exists).

**Two obstructions (Mine; elementary).** (F) with $\phi_C\ge0$ holds iff both:
- **(a) support shape:** $\mathrm{supp}\,P=\{x:x_C\in L_C\ \forall C\}$ for some sets $L_C$ of allowed clique patterns. The support must be a nearest-neighbour constraint space, in [Ch]'s term "n.n.constraint space".
- **(b) local cohomology:** on the support, $\log P=\sum_C\varphi_C(x_C)$ with $\varphi_C$ defined on $L_C$.

(⇒ take $L_C=\{\phi_C>0\}$; ⇐ put $e^{\varphi_C}$ on $L_C$ and $0$ off it.)
- Moussouris' $P_M$ fails (a) only: $\log P_M$ is constant, but $\Omega_0$ is not cut out by edge patterns.
- $P^*$ fails (b) only. $\Omega_{00}=\{x_3=x_4\}$ *is* an edge constraint, but $\log P^*$ carries the 3-body monomial $x_1x_2x_3$.
- The support constraint $x_3=x_4$ **contracts the edge $\{3,4\}$**. That turns $C_4$ into a triangle on $\{1,2,(34)\}$, which carries a 3-body term that is legitimate *on the contracted graph*. Pairwise CI on $C_4$ cannot see it, because Lemma B would need configurations with $x_3\ne x_4$.

Obstruction (b) is exactly $\mathbf M_X\ne\mathbf G_X$ in [CM]'s language. [Ch] p. 2: "the difference of their dimensions measures the extent to which $X$ satisfies the conclusion of the Hammersley–Clifford theorem." On $\Omega_{00}$, a short count (Mine) gives:
- Markov cocycles: all functions on the 8 points modulo constants, dimension 7, because every law on $\Omega_{00}$ is pairwise Markov.
- Gibbs cocycles: edge monomials $\{1,x_1,x_2,s,x_1x_2,x_2s,sx_1\}$ with $s=x_3=x_4$, dimension 6 modulo constants.

So $\mathbf M/\mathbf G$ is spanned by $x_1x_2s$.

### 4.5 Positive results without positivity (Source)

- **Chordal graphs.** [GL] Theorem 3.2, after Lauritzen: "If the graph $\mathcal G=(\Lambda,\mathcal B)$ is chordal, then a random field $P$ is $\mathcal B$-global-Markov if and only if it is $\mathcal B$-$\mathcal B$-Gibbs for some potential $\phi$." [Ch]: "a bipartite graph is decomposable if and only if it is a forest." Paths and trees are chordal, so on them (G) ⇔ (F) with zeros allowed. But (L) ⇏ (G) still fails on the path (§4.2).
- **Safe symbol.** [CM] Theorems 3.1–3.2 hold on any countable locally finite graph.
- **Pivot property.** [CM] §3.2: "for all $(x,y)\in\Delta_X$ … there exists a finite sequence of points … such that each $(x^{(i)},x^{(i+1)})$ is a pivot move". Examples: safe symbol; the 3-coloured chessboard; $r$-colourings of $\mathbb Z^d$ with $r\ge 2d+2$ (Prop. 3.4); homomorphisms to dismantlable graphs.
  - Prop. 3.5: pivot property on a shift-invariant topological Markov field $\Rightarrow\dim\mathbf M^\sigma_X\le|B_{\{0\}\cup\partial\{0\}}(X)|^2$.
- **The $X_r$ family** (3-coloured chessboard at $r=3$): $x_n-x_m=\pm1\bmod r$ on edges of $\mathbb Z^d$, $d\ge2$, $r\ne1,2,4$.
  - Prop. 5.2: "$\dim\mathbf M^\sigma_{X_r}=r$", with basis $M_i(x,y)=\sum_nN_i(\hat x_n,\hat y_n)$ through height-function lifts $\hat x$.
  - Prop. 5.3: "$\dim\mathbf G^\sigma_{X_r}=r-1$" (the coefficients must satisfy $\sum\alpha_i=0$).
  - Cor. 5.4: every shift-invariant Markov cocycle is $M_0+\alpha\hat M$ with $M_0$ Gibbs and $\hat M(x,y)=\sum_n(\hat y_n-\hat x_n)$, "the volume between height functions".
  - Prop. 5.5: nevertheless $\mathbf M_{X_r}=\mathbf G_{X_r}$, so every Markov cocycle is Gibbs for some *non*-shift-invariant interaction.
  - Cor. 5.6: $\mathbf G^\sigma_X\ne\mathbf G_X\cap\mathbf M^\sigma_X$.
  - Thm 6.1: "Any shift-invariant Markov random field adapted to $X_r$ is a Gibbs state for some shift-invariant nearest neighbor interaction." The reason is Prop. 6.2: an invariant MRF whose cocycle has $\alpha\ne0$ must be **frozen**. The intuition [CM] give: "$\sum\alpha_i>0$ indicates an inclination to raise the height function. However $\sigma$-invariance implies the existence of a well defined 'slope' … Unless this slope is extremal … this will lead to a contradiction."
  - Prop. 7.1: adapted shift-invariant MRFs on $X_r$ are fully supported or frozen.
- **No pivot property.** [CM] §9 constructs a $\mathbb Z^2$ subshift with "an uncountable family of linearly independent shift-invariant Markov cocycles", hence "a shift-invariant Markov random field which is not Gibbs for any shift-invariant finite range interaction". The Markov specification "cannot be 'given by a finite number of parameters'".
- **Folding.** [Ch] Theorems 4.1–4.2: on a bipartite graph, strong config-folds and config-unfolds of Hammersley–Clifford spaces are Hammersley–Clifford. Cor. 4.8: $\mathbf M^{Gr}_X/\mathbf G^{Gr}_X\cong\mathbf M^{Gr}_{X^a}/\mathbf G^{Gr}_{X^a}$ under a strong config-fold. Consequences: $\mathrm{Hom}(G,H)$ is Hammersley–Clifford for dismantlable $H$, and $\mathrm{Hom}(G,C_4)$ for bipartite $G$.
- **Dimension one.** [CM] §1: "When the underlying graph is a Cayley graph of $\mathbb Z$ the Markov random field is shift-invariant and $A$ is finite, the conclusions of the Hammersley–Clifford Theorem hold without any further assumptions … in that setting any Markov random field is a stationary Markov chain. Even when the underlying graph is a Cayley graph of $\mathbb Z$, this conclusion can fail for countable $A$ … or if we drop the assumption of shift-invariance [Dobrushin 1968]."

---

## 5. The converse direction, dynamics, and approximate classical versions

### 5.1 Gibbs ⇒ Markov is exact; Markov belongs to the specification

§3.4 gives the converse for every graph, every potential and every zero pattern on the same graph. It is a property of the conditional structure, not of uniqueness or mixing:
- [CM]: two MRFs can share a specification (low-temperature 2D Ising).
- [Y] (via [arxiv-2609.38007.md](arxiv-2609.38007.md) §10): "a classical Gibbs distribution retains its Markov property even at a thermal phase transition."

### 5.2 Markov fields are the equilibria of reversible local dynamics (Source: [Pr] abstract)

"The set of Markov random fields is also shown to be the same as the equilibrium states of time-reversible birth/death processes with nearest neighbour interactions defined on the graph." *Mine:* this is the classical ancestor of two quantum statements:
- "approximate detailed balance ⇒ local Markov" ([C26] Theorem II.1, from Bergamaschi–Chen–Vazirani);
- "recovery by a detailed-balanced Lindbladian" ([CR]).

Every MRF is stationary for *some* reversible local dynamics, but that dynamics' spectral gap is unconstrained.

### 5.3 Approximate classical Hammersley–Clifford

**(i) The exact identity (Source: [ILW] §II, eqs. (11)–(21)).** "for any joint distribution $P_{XYZ}$ of three random variables …
$$I(X:Z|Y)=\min\{D(P\|Q): Q\text{ Markov}\}."$$
The minimiser keeps $P_Y$, $P_{X|Y}$ and $P_{Z|Y}$. The extraction of their (12) reads "$P_{Z|X}$", but the derivation (18) sets $Q_{Z|Y}=P_{Z|Y}$, so (12) is a typo. For chains, the sibling check C2 of [hdx-spectral-independence.md](hdx-spectral-independence.md) verifies $D(P\|Q_{\rm Markov})=\sum_iI(X_{i+1};X_{<i}\mid X_i)$.

**(ii) The maximum-entropy route (Source: [KB] Theorem 1; [LZ] for chordal graphs).**
- [KB]'s proof (§9.1 below) uses only strong subadditivity and maximum entropy, so it applies verbatim to commuting (classical) states: an $\varepsilon$-approximate Markov chain is within $m\varepsilon$ in relative entropy of a nearest-neighbour Gibbs law. [KB] note "Theorem 1 is not restricted to full-rank states."
- For chordal graphs the maxent completion is the junction-tree formula, and the defect is [LZ]'s "global information" $g\mathrm I(\mathcal G)_\rho$. "In the two-clique case this is the conditional mutual information." "For classical positive distributions on a chordal graph, every pairwise consistent family of clique marginals has a Markov extension, given by the junction-tree formula, and this extension is also the maximum-entropy completion" ([LZ] Remark 3.15).

**(iii) A pointwise approximate Hammersley–Clifford (Mine; checked, §15 H2).** Let $P>0$ and suppose $|\Delta_i\Delta_j\log P(x)|\le\delta$ for every non-edge and every $x$. By Lemma B's block expansion,
$$\|\Phi_A\|_\infty\le2^{|A|-2}\,\delta\quad(A\text{ not complete}).$$
With $P_G=\exp(\sum_{C}\Phi_C)/Z_G$ over complete $C$, put $R=\sum_{A\text{ not complete}}\Phi_A$. Then $Z_G=\sum_xP(x)e^{-R(x)}$ gives $|\log Z_G|\le\|R\|_\infty$, and so
$$\|\log P-\log P_G\|_\infty\le2\sum_{A\ \text{not complete}}2^{|A|-2}\delta\le\tfrac12\,3^n\delta,\qquad D(P\|P_G)\le\tfrac12\,3^n\delta .$$
- The per-coefficient constant $2^{|A|-2}$ is attained (H2 ratio $1.000$).
- The global constant is exponential in $n$. That is the price of a sup-norm, Boolean-lattice argument.
- The entropic routes (i)–(ii) cost $O(\#\text{separators})$ on chains and chordal graphs. No retrieved source gives an entropic approximate Hammersley–Clifford for non-chordal classical graphs (Mine: not found, not claimed absent).

### 5.4 Heredity: pinning is free, marginalising fills in (Derivation)

**Pinning.** If $P=\prod_C\phi_C$ on $G$ and $P(x_S=s)>0$, then
$$P(\cdot\mid x_S=s)\propto\prod_C\phi_C(x_{C\setminus S},s_{C\cap S}),$$
and each $C\setminus S$ is complete in the induced graph $G[V\setminus S]$. So **pinnings of Markov fields are Markov fields on induced subgraphs, with external fields**: exact, positivity-free.
- This is precisely the hereditary property used by spectral independence. Links of the Gibbs complex are pinnings ([hdx-spectral-independence.md](hdx-spectral-independence.md) §3.2).
- Hammersley–Clifford guarantees the pinned law is again in the class.

**Marginalising.** Summing out $x_S$ produces a factor on $\partial S$, i.e. the neighbours of $S$ are filled in to a clique (Mine; standard variable elimination, *from memory*). [LP] Lemma 4.9 is the quantum version, stated exactly: tracing out $u$ from a quantum Markov network gives one on $G''$, "the graph obtained by adding … an edge between every pair of distinct neighbors of $u$".

---

## 6. Quantum conditional independence

### 6.1 Equality in strong subadditivity: Petz sufficiency, structure, recovery (Source: [HJPW])

$I(A:C|B)_\rho=S(AB)+S(BC)-S(ABC)-S(B)\ge0$ (Lieb–Ruskai). [HJPW] (6) writes it as a relative-entropy difference under the partial trace:
$$I(A:C|B)=S(\rho_{ABC}\|\rho_A\otimes\rho_{BC})-S(\rho_{AB}\|\rho_A\otimes\rho_B).$$

**Theorem 3 of [HJPW] (Petz).** "$S(\rho\|\sigma)=S(T\rho\|T\sigma)$ if and only if there exists a quantum operation $\hat T$ such that $\hat TT\rho=\rho$, $\hat TT\sigma=\sigma$." On the support of $T\sigma$:
$$\hat T\alpha=\sigma^{1/2}T^*\big((T\sigma)^{-1/2}\alpha(T\sigma)^{-1/2}\big)\sigma^{1/2}.$$
Specialised (their (10)–(11)), equality in SSA iff $\rho_{ABC}=(\mathrm{id}\otimes\hat R)\rho_{AB}$, with the Petz map of $\rho_{BC}$ for $\mathrm{Tr}_C$ ([FR] (15)):
$$X_B\mapsto\rho_{BC}^{1/2}\big(\rho_B^{-1/2}X_B\rho_B^{-1/2}\otimes\mathrm{id}_C\big)\rho_{BC}^{1/2}.$$

**Theorem 6 of [HJPW] (structure).** "A state $\rho_{ABC}$ … satisfies strong subadditivity … with equality if and only if there is a decomposition of system $B$ as $\mathcal H_B=\bigoplus_j\mathcal H_{b^L_j}\otimes\mathcal H_{b^R_j}$ … such that
$$\rho_{ABC}=\bigoplus_jq_j\,\rho_{Ab^L_j}\otimes\rho_{b^R_jC}."$$

*Proof architecture ([HJPW] §V and App. A).*
1. From the Markov condition, $\varphi=\mathrm{Tr}_C\circ\hat R$ leaves invariant every state $\mu$ obtained from $\rho_{AB}$ by conditioning on effects $0\le M\le1$ on $A$.
2. **Koashi–Imoto (Thm 10)** applied to that invariant family gives the decomposition. Its algebraic proof:
   - **Lemma 12**: for a channel $F$, the Cesàro mean $P^*=\lim_N\frac1N\sum_{n\le N}(F^*)^n$ "is a conditional expectation onto the $*$-subalgebra $\mathcal A_F=\{X:F^*(X)=X\}=\{B_i,B_i^*\}'$", the commutant of the Kraus operators;
   - $\mathcal A_0=\bigcap_F\mathcal A_F$ is realised by one $F_0$ (an average);
   - **Lemma 13**: every finite-dimensional $*$-subalgebra is $\bigoplus_j B(\mathcal H_{b^L_j})\otimes1$, and every unital CP projection onto it is $X\mapsto\bigoplus_j\mathrm{Tr}_{b^R_j}(\Pi_jX\Pi_j(1\otimes\omega_j))\otimes1$.
3. The Stinespring unitary of $\hat R$ has the form $\bigoplus_j1\otimes U_j$ (their (15)), which yields the product form.

Further statements in [HJPW]:
- **Cor. 7:** if $I(A:C|B)=0$ then $\rho_{AC}$ is separable. "Conversely, for each separable state $\rho_{AC}$ there exists an extension $\rho_{ABC}$ such that $I(A:C|B)=0$."
- **Interpretation:** "there is information in the $B$ system which can be obtained by a non-demolition measurement, conditioned upon which the quantum state factorises."
- **Accardi–Frigerio** (as [HJPW] state them): a state is Markovian if for every $m$ there is a unital CP map $T_{m,m+1}:\mathcal A_{m+1}\to\mathcal A_m$ "which leaves the state … invariant and the subalgebra $\mathcal A_{m-1}$ fixed … quasi-conditional expectation; its dual is the quantum analogue of the Markov kernel."
- **Open question** at the end of [HJPW]: "if a state almost satisfies strong subadditivity, does it mean that its structure is close in some sense to the form of theorem 6?" (Answered in §8.)

**Two algebraic restatements.**
- [BP] (3)–(4): $I=0\iff\rho_{ABC}=\Lambda_{AB}\Lambda_{BC}$ with $[\Lambda_{AB},\Lambda_{BC}]=0$. For $\rho>0$, $\rho=e^{H_{AB}+H_{BC}}$ with $[H_{AB},H_{BC}]=0$.
- Ruskai's logarithmic form ([LP] (28)): $\log\rho_{ABC}+\log\rho_B=\log\rho_{AB}+\log\rho_{BC}$.

### 6.2 Graphoid axioms, and intersection for full-rank states

[LP] Theorem 4.5: "The quantum dependency model is a graphoid." Symmetry, decomposition, weak union and contraction all follow from SSA and the chain rule. [LP]: "The analogous quantum property would be to require that $\rho_V$ is a strictly positive operator … but we have not been able to prove that this property implies intersection … unlike in the classical case, we cannot conclude that the global Markov property holds for positive quantum Markov networks."

**Resolved: Proposition 2.5 of [LZ] (2026).** "For strictly positive density operators, quantum conditional independence satisfies intersection: (Q5) $A\perp\!\!\!\perp_QB\mid C\cup D$ and $A\perp\!\!\!\perp_QD\mid(B\cup C)\Rightarrow A\perp\!\!\!\perp_Q(B\cup D)\mid C$." "The intersection property ensures that different variants of Markov properties (pairwise, local, and global) are equivalent."

*Proof architecture ([LZ] App. A).*
- **A.1–A.4 (variational representation).** On the weighted operator space $\langle X,Y\rangle_\rho=\mathrm{Tr}(X^*Y\rho)$,
  $$D(\rho\|\tau)-D(\rho_A\|\tau_A)=\int_0^\infty\Big[\max_Zf^{\rho,\tau}_t(Z)-\max_Xf^{\rho,\tau}_t(J_AX)\Big]dt,$$
  with $J_AX=X\otimes I_B$. The global maximiser is unique: $X_t=(tI+\Delta_{\tau,\rho})^{-1}I$ for the relative modular operator $\Delta_{\tau,\rho}$.
- **Prop. A.9 (equality criterion).** Equality holds iff for every $t>0$ there is $Y_t$ with $X_t=J_A(Y_t)$: "the global maximizer already belongs to the subspace".
- **Prop. A.10.** Equality iff $\rho=\tau^{1/2}(\tau_A^{-1/2}\rho_A\tau_A^{-1/2}\otimes I_B)\tau^{1/2}$ (the Petz formula). The proof expands the resolvent in powers of $t^{-1}$ and uses a polynomial $p$ with $p(s)=s^{-1/2}$ on the spectra.
- **A.6 (intersection).**
  1. Set $\tau=\rho_A\otimes\rho_{BCD}$. Then $D(\rho\|\tau)-D(\rho_{ACD}\|\tau_{ACD})=I(A:B|CD)$ and $D(\rho\|\tau)-D(\rho_{ABC}\|\tau_{ABC})=I(A:D|BC)$.
  2. Both vanish, so by A.9, $X_t\in(L(\mathcal H_{ACD})\otimes I_B)\cap(L(\mathcal H_{ABC})\otimes I_D)$.
  3. "This intersection is $L(\mathcal H_A\otimes\mathcal H_C)\otimes I_B\otimes I_D$" (a basis-expansion argument).
  4. A.9 once more gives $I(A:BD|C)=0$.

*Mine.* The classical intersection axiom needs a *connected support* (§4.2). The quantum one needs a *faithful state*, and then reduces to the fact that the tensor subalgebras $L(ACD)\otimes1$ and $L(ABC)\otimes1$ intersect in $L(AC)\otimes1$. That is a commuting-square fact about tensor products (Popa's term, *from memory*). Faithfulness is the quantum form of "the vacuum substitutions stay inside the support".

### 6.3 The quantum marginal problem for Markov completions (Source: [LZ])

- **Definition.** $T(\mathcal R)=\exp\{\log\rho_{A\cup C}+\log\rho_{B\cup C}-\log\rho_C\}=\rho_{A\cup C}\odot\rho_{B\cup C}\odot\rho_C^{-1}$. Here $C$ is the separator, and $\odot$ is the "exp of sum of logs" product that [LP] also use.
- **Theorem 3.1.**
  - $\mathrm{Tr}\,T(\mathcal R)\le1$, by "Lieb's three-matrix inequality" and partial-trace cancellation down to $\sum_i\int_0^\infty\lambda_i^2(t+\lambda_i)^{-2}dt=1$.
  - $\mathrm{Tr}\,T=1\iff T\in M(\mathcal R)\iff$ a Markov completion exists.
  - "When these conditions hold, $T(\mathcal R)$ is the unique Markov completion."
- **Lemma 3.2.** $D(\omega\|T)+1-\mathrm{Tr}\,T=I(A:B|C)_\omega+\Delta_{\mathcal R}(\omega)$, with $\Delta_{\mathcal R}\ge0$. Their $D$ for an unnormalised second argument includes "$-1+\mathrm{Tr}(T)$". For $\omega=\rho$ this is $\mathrm{Tr}\,\rho(\log\rho-\log T)=I$; checked to $3\cdot10^{-15}$, §15 H6.
- **The operator $K$.** The two one-sided reconstructions agree iff $K=\rho_{A\cup C}^{1/2}\rho_C^{-1/2}\rho_{B\cup C}^{1/2}$ is normal, and then $T=KK^*=K^*K$.
- **Theorem 3.10 (chordal).** The same trace criterion holds; the completion is unique and is the maximum-entropy element.
- **Remark 3.15.** "the maximum-entropy completion may exist and have log-linear form without being quantum Markov."
- **Example 4.3 / Lemma 4.6.** Three-qubit Pauli marginals that are consistent and strictly feasible when $\varepsilon^2+\delta^2<1$ "admit no Markov completion when $\varepsilon\delta\ne0$."
- **Remark 4.8.** Markov states can have non-commuting overlapping marginals, e.g. $\rho=\sigma_1\otimes\sigma_{23}$.

---

## 7. Quantum Hammersley–Clifford

### 7.1 Definitions

- **[LP] §4.4:** $(G,\rho_V)$ is a quantum Markov network if the quantum graphoid satisfies the *local* Markov property on $G$. It is positive if $\rho_V$ has full rank.
- **[BP] §II.B:** a quantum Markov network requires "$I(A:C|B)=0$ for all disjoint subsets $A,B,C\subset V$ such that $B$ shields $A$ from $C$" (the global property). "the condition for regions $A$, $B$, and $C$ that span only a subset of the vertices follow[s] from strong subadditivity."
- **[KKB]:** "$\rho_{V_0}$ is the quantum Markov network on $V_0$" if $I(A:C|B)=0$ for $d_{A,C}>0$.

By [LZ] Prop. 2.5, all exact definitions coincide for full-rank states. This is **not** the open "pairwise/local/global" question of [CR] and [Kuw]. That question concerns the *approximate*, size-dependent versions: [CR] "In the quantum case, these three are not known to be equal, even allowing for approximations"; [Kuw] "the conditions for such equivalence remain unclear". The sources do not conflict.

### 7.2 Markov ⇒ Gibbs of a clique-local (not necessarily commuting) Hamiltonian

**[LP] Theorem 4.7.** "If $(G,\rho_V)$ is a positive quantum Markov network then there exist positive operators $\sigma_C$ acting on the cliques of $G$ … such that $\rho_V=\bigodot_{C\in\mathcal C}\sigma_C$." The proof, "very similar to a standard proof for the classical case", is in their App. B.

**[BP] Theorem 2.** "If the pair $(\rho,G)$ is a positive quantum Markov network, then the state $\rho$ can be expressed as $\rho=e^H$ where $H=\sum_{Q\in\mathcal C}h_Q$ … on the particles located in cliques $Q$." Their proof, in full:
1. Decompose $H=\log\rho=\sum_{X\subseteq V}K_X$ into cumulants with "$\mathrm{Tr}_YK_X=0$ for any $Y\subseteq X$" ($Y\ne\emptyset$). These are Hilbert–Schmidt orthogonal, and $\mathrm{Tr}(HK_X)=\mathrm{Tr}K_X^2$.
2. If $X$ is not a clique, pick non-adjacent $a,c\in X$. $B=V-a-c$ shields $a$ from $c$, so $H=H_{aB}+H_{Bc}$ (by (4)).
3. Then $\mathrm{Tr}(HK_X)=\mathrm{Tr}(H_{aB}[\mathrm{Tr}_cK_X])+\mathrm{Tr}(H_{Bc}[\mathrm{Tr}_aK_X])=0$.
4. Hence $K_X=0$. $\square$

*Mine (THEOREM, elementary; checked §15 H5).* With $E_i=\tau_i\otimes\mathrm{Tr}_i$ (normalised partial trace on site $i$, tensored back with the identity), the cumulant is
$$K_X=\prod_{i\in X}(\mathrm{id}-E_i)\prod_{i\notin X}E_i\,(\log\rho).$$
This is **Lemma A with the vacuum evaluation $\varepsilon_i$ replaced by the tracial conditional expectation $E_i$**. [BP]'s step 3 is Lemma B: pairwise Markov with full rank gives $\log\rho=H_{aB}+H_{Bc}$, hence $(\mathrm{id}-E_a)(\mathrm{id}-E_c)\log\rho=0$, the quantum mixed difference. A matrix algebra has no characters to evaluate at, so the vacuum must be a state. *Corrected in review:* the trace is not the only choice. [BP] use the normalised trace. [LP]'s own proof of their Theorem 4.7 (App. B, Lemmas B.1–B.2) uses a pure product vacuum $|\alpha\rangle$, that is, [HC]'s blackening with "blacken" read as "compress to $|\alpha\rangle$" (§3.8). With a product reference measure the same construction is the Hoeffding/Efron–Stein decomposition (*from memory*; cf. [arxiv-2609.38007.md](arxiv-2609.38007.md) B5). Only the **pairwise** condition and **full rank** are used, as classically.

### 7.3 Commuting clique Hamiltonian ⇒ Markov

**[BP] Theorem 3.** "Let … $H=\sum_{Q\in\mathcal C}h_Q$, $[h_Q,h_{Q'}]=0$ … Then $(\rho,G)$ is a positive quantum Markov network for $\rho=\frac1Ze^H$." Proof:
1. For full partitions ($B$ shields $A$ from $C=V-A-B$), "no clique $Q$ can overlap simultaneously with region $A$ and $C$". Assign each clique to $AB$ or $BC$; the two sums commute, and (4) gives $I=0$.
2. For general shielded triples, expand $A,C$ to $AA',CC'$ filling $V$ with $B$ still shielding, and use $I(AA':CC'|B)\ge I(A:C|B)$ (SSA).
3. "the above result may be extended to non-positive density matrices that are projectors onto eigenspaces of arbitrary commuting Hamiltonians, including for example the projector onto the code space of any local stabilizer code."

### 7.4 Markov ⇒ *commuting* Hamiltonian: chains, trees, triangle-free graphs; fails with triangles

- **Chains (Source: [BP] §II.A).** Iterating (4) site by site gives $\rho=e^H$, $H=\sum_kh_{k,k+1}$ with $[h_{k,k+1},h_{k',k'+1}]=0$: "a complete equivalence between positive quantum Markov chains and local commuting Hamiltonians in 1D." [KB] §I.A states the same.
- **Trees (Source: [BP] Theorem 5, Poulin–Hastings).** "Let $G$ be a tree … $(\rho,G)$ is a positive quantum Markov network if and only if … $H=\sum_{Q}h_Q$, $[h_Q,h_{Q'}]=0$."
- **Triangle-free graphs (Source: [BP] Theorem 4).** "If the pair $(\rho,G)$ is a positive quantum Markov network and $G$ contains only two-body cliques, then … $H=\sum_{Q}h_Q$, $[h_Q,h_{Q'}]=0$."
  - *Lemma 1 (genuine operators):* "Let $H_{AB}$ be a genuine operator on $AB$ and $H_{BC}$ be a genuine operator on $BC$. Then $[H_{AB},H_{BC}]$ is either 0 or a genuine operator on $AB'C$ where $B'\subseteq B$ and $B'\ne\emptyset$." Proof: operator-Schmidt decompositions; "the commutator $[G^j_B,R^k_B]$, if non-zero, cannot be proportional to the identity for finite dimensions."
  - *Step 1:* all two-body cumulants commute. Expand $[H_{AB},H_{BC}]=0$ in cumulants; the terms $[K_{ab},K_{bc}]$ have pairwise distinct supports, so each vanishes by Lemma 1; vary the shields.
  - *Step 2:* split each one-body cumulant $K_u=h_u+\sum_{v\in N(u)}G^v_u$. Take $A\subset N(u)$, $B=N(A)$, $C=V-A-B$; derive $[K^A_u,K^C_u]=[K_{au},K^C_u]=[K^A_u,K_{uc}]=0$ (their (30)); apply Theorem 5 to the star $G(u)$.
  - *Where triangle-freeness is used:* (24) "cannot contain two-body cumulants $K_{aa'}$ because this would create a triangle", and in (26) "the two-body cumulants $K_{bb'}$ … cannot include node $u$."
- **Failure with triangles (Source: [BP] §III.C, Fig. 2).** On the 5-vertex wheel (rim $1,2,3,4$, centre $5$), take $H=h_5+h_C+h_4+h_B$ with
  $$h_5=\sigma^z_1\sigma^z_2\sigma^y_5,\quad h_C=\sigma^z_2\sigma^z_3\sigma^x_5,\quad h_4=\sigma^z_3\sigma^z_4\sigma^y_5,\quad h_B=\sigma^z_4\sigma^z_1\sigma^x_5 .$$
  - The only shields are $(1;\{2,4,5\};3)$ and $(2;\{1,3,5\};4)$. For each, $H_{AB}$ and $H_{BC}$ commute ($h_5+h_B$ vs $h_C+h_4$; $h_5+h_C$ vs $h_B+h_4$), so the Gibbs state is a positive quantum Markov network.
  - "Nevertheless, $\rho$ is not the Gibbs state of a local commuting Hamiltonian as the individual terms in the Hamiltonian do not mutually commute."
  - *Checked (§15 H4):* both CMIs below $1.3\cdot10^{-15}$ at $\beta\in\{0.3,1,2.5\}$; $\|[h_5,h_C]\|=11.3$; a generic triangle Hamiltonian on the same wheel gives $I(1:3|245)=0.116$.
- **Shield-commuting Hamiltonians (Source: [BP] Def. 1).** $H$ is shield commuting "when it can be written as $H=H_{AB}+H_{BC}$ with $[H_{AB},H_{BC}]=0$ whenever $B$ shields $A$ from $C$." By (4), positive quantum Markov networks are exactly the Gibbs states of shield-commuting Hamiltonians. Two further remarks from [BP]:
  - "there is no obvious efficient way to test it";
  - the lattice version of the example is fixed by coarse-graining ("no 'large-scale' obstructions to commutation").
- **Open question (Source: [BP] §IV).** "whether there exist Hamiltonians that satisfy this global commutation property … but that cannot be transformed into a local commuting Hamiltonian under a suitable renormalization procedure." If so, they "would reveal a new phase of matter that exhibits quantum non-locality without long range entanglement."
- **Splitting (Source: [BP]).** "for Hamiltonians containing only two-body commuting interactions, the Hilbert space of each vertex must split into a direct sum of factor spaces as in Eq. 2" [Bravyi–Vyalyi]. "commutation of genuine operators that act non-trivially on more than two subsystems does not imply that the local Hilbert spaces will split."

*Mine.* Classically the factors of a factorisation commute automatically, so the only obstructions are support-topological (§4.4). Quantumly, each shield $B$ comes with its own Koashi–Imoto splitting $\bigoplus_jB(\mathcal H_{b^L_j})\otimes1$ (§6.1). A commuting clique Hamiltonian exists when the splittings attached to different shields are compatible. Triangles let overlapping shields carry incompatible splittings. **The quantum obstruction is algebraic (incompatible conditional-expectation structures), not topological.**

### 7.5 Non-commuting local Gibbs ⇏ Markov

- **[LP] Example 4.8.** For the antiferromagnetic Heisenberg chain on three qubits, $\rho(\beta)=e^{-\beta H}/Z$ "has the form eq. (92), but for any finite $\beta$ it has a non-zero mutual information between A and C conditioned on B" (Fig. 2).
- *Checked (§15 H5):*
  - $I(A:C|B)=0.082,\ 0.315,\ 0.412,\ 0.415$ bits at $\beta=0.5,1,2,4$;
  - yet the cumulants of $\log\rho$ on $\{A,C\}$ and $\{A,B,C\}$ vanish to $10^{-14}$;
  - the diagonal (Ising) chain has $I=0$.
- [Kuw] (S.29) context: "when the Hamiltonian is non-commuting, the quantum Markov property (S.29) breaks down in the exact sense."
- [C26]: "For quantum Gibbs states, however, an exact analog fails [Kuw24, BLMT25, KK25]."

*Mine: where exactly the converse breaks.* Classically three statements are equivalent for $P>0$:
- (i) $\Delta_a\Delta_c\log P=0$ (the joint log-density splits);
- (ii) $\log P_{ABC}-\log P_{AB}-\log P_{BC}+\log P_B=0$ (the log-marginals split);
- (iii) $I(A:C|B)=0$.

Quantumly, (ii) ⇔ (iii) (Ruskai/Petz, full rank) and (iii) ⇒ (i) ([BP]), but **(i) ⇏ (iii)**: the Heisenberg chain satisfies (i) and not (iii). The gap is that a partial trace of $e^{H_{AB}+H_{BC}}$ is not $e^{\text{(local)}}$ unless the terms commute. The reduced state's "effective Hamiltonian" ([KKB] (10)–(11): $\tilde H_L=H_L+\Phi_L$; [Kuw]: "entanglement Hamiltonian", "Hamiltonian of mean force") is only quasi-local. The CMI is the expectation of the **mixed second difference of the log-marginals** over the square $B\subset AB,BC\subset ABC$:
$$I(A:C|B)=-\mathrm{Tr}\,\rho\,\big[\ell(AB)+\ell(BC)-\ell(ABC)-\ell(B)\big],\qquad\ell(L)=\log\rho_L\otimes1_{L^c},$$
which is the operator [KKB] expand, $\tilde H(A:C|B)$ (their (A20)–(A21)). So:
- the Hammersley–Clifford Möbius inversion acts on the *joint* $\log P$ along the lattice of vacuum substitutions;
- the CMI acts on the *marginal* logs along the lattice of regions.

Classically both vanish together; quantumly the second is the one that matters, and only the first is controlled by locality of $H$.

### 7.6 Without full rank; marginals; fill-in

- **Code states and topological order (Source: [BP] §I).** The projector onto a stabilizer code space is a quantum Markov network. "In a system with topological order, the conditional mutual information would be zero for any topologically trivial region A, however, there must exist a non-trivial region, for instance when A is a ribbon wrapping around a torus, such that the mutual information … conditioned on the boundary of A, is equal to the topological entanglement entropy … The state obtained from the uniform mixture of all ground states … can form a Markov network however."
- **GHZ (Source: [KB] §II).** The GHZ state is "locally" Markov (each one-site-traced state is an exact quantum Markov chain), but "$I(A:C|B)=1$ for any tripartition ABC of the whole system where B shields A from C."
- **Cluster state (Source: [Kuw] (S.621)).** For the 1D cluster Hamiltonian with periodic boundary, the reduced state on the odd sublattice has "an infinite correlation length for the CMI". [KKB] add that the cluster state "is globally a Markov network, but not for particular selections of $V_0$."
- **Partitions (Source: [Kuw]).** "the condition $\Lambda=A\sqcup B\sqcup C$ is crucial. If $A\sqcup B\sqcup C\subset\Lambda$, even for commuting Hamiltonians, there exists a counterexample to (7) at low temperatures."
- **Fill-in (Mine).** These are all fill-in phenomena (§5.4; [LP] Lemma 4.9). Tracing out a site of an exact quantum Markov network keeps it Markov, but only on the graph where that site's neighbours are joined. Repeated tracing cascades the fill-in until the graph is nearly complete. A marginal of a Markov network is therefore Markov only for a much denser graph. Classically this is variable elimination; quantumly it is the same lemma.

---

## 8. Recovery maps and approximate quantum Markov chains

### 8.1 Fawzi–Renner (Source: [FR])

**Theorem 5.1.** "For any density operator $\rho_{ABC}$ on $A\otimes B\otimes C$, where $A$, $B$, and $C$ are separable Hilbert spaces, there exists a trace-preserving completely positive map $\mathcal T_{B\to BC}$ … such that
$$2^{-\frac12I(A:C|B)_\rho}\le F\big(\rho_{ABC},(\mathcal I_A\otimes\mathcal T_{B\to BC})(\rho_{AB})\big)."$$
In finite dimensions the map is a rotated Petz map, $X_B\mapsto V_{BC}\rho_{BC}^{1/2}(\rho_B^{-1/2}U_BX_BU_B^\dagger\rho_B^{-1/2}\otimes\mathrm{id}_C)\rho_{BC}^{1/2}V^\dagger_{BC}$.

**Consequences ([FR] p. 2–3).**
- (5): $\inf_\sigma D_{1/2}(\rho\|\sigma)\le I$.
- (6): $\frac1{\ln2}\Delta(\rho_{ABC},\sigma_{ABC})^2\le I(A:C|B)$, with $\Delta$ the trace distance.
- Converse (10): $I(A:C|B)\le7\log_2(\dim A)\sqrt{\Delta}$ for $\Delta\le\frac1{11}$; "a term proportional to the logarithm of the dimension of A is necessary in general."
- Classical case (11)–(14): with $B$ classical, the recovered state is a Markov chain and "$D(\rho_{ABC}\|\sigma_{ABC})=I(A:C|B)_\rho$".
- Remark 5.3: the recovery map can be chosen to reproduce the $B$ and $C$ marginals exactly.
- "the reconstructed state … is not necessarily a Markov chain. (Note that this is a major difference to the classical case …)"

**Proof architecture.**
1. Typicality bounds on the relative entropy via one-shot entropies ($D_H^\epsilon$, $D^\epsilon_{\max}$), App. A and §2.
2. A de Finetti reduction (§3), used for the fidelity of permutation-invariant states (§4).
3. Main step (§5): typical-eigenspace projectors $\Pi_{B^n}$, $\Pi_{B^nC^n}$ of $\rho_B^{\otimes n},\rho_{BC}^{\otimes n}$ (Lemma 2.5). The relative-entropy bound (108) is turned into fidelity, then split over polynomially many eigenvalue classes (Lemma B.7). Lemma 4.2 produces unitaries $U_B$, $V_{BC}$ at $\mathrm{poly}(n)$ cost.
4. Take the $n$-th root, using $H(BC)-H(B)-H(C|AB)=I(A:C|B)$.

### 8.2 Universal, explicit recovery: modular-time-averaged Petz map (Source: [JRSWW])

**Theorem 2.1.** For $\sigma\ge0$, $\rho$ with $\mathrm{supp}\,\rho\subseteq\mathrm{supp}\,\sigma$, and any channel $\mathcal N$,
$$D(\rho\|\sigma)\ge D(\mathcal N(\rho)\|\mathcal N(\sigma))-2\int_{\mathbb R}dt\,\beta_0(t)\log F\big(\rho,(\mathcal R^{t/2}_{\sigma,\mathcal N}\circ\mathcal N)(\rho)\big),$$
where
$$\mathcal R^t_{\sigma,\mathcal N}(X)=\sigma^{-it}\mathcal P_{\sigma,\mathcal N}\big(\mathcal N(\sigma)^{it}X\mathcal N(\sigma)^{-it}\big)\sigma^{it},\qquad\beta_0(t)=\frac\pi2\big(\cosh(\pi t)+1\big)^{-1}.$$

**Remark 2.2.** By concavity, one map $\mathcal R_{\sigma,\mathcal N}=\int\beta_0(t)\mathcal R^{t/2}_{\sigma,\mathcal N}dt$ works. Remark 2.4 lists its properties: universality ("does not depend on $\rho$"), perfect reconstruction of $\sigma$, normalisation, stabilisation.

**Cor. 4.1.** $I(A:C|B)_\rho\ge-2\log F(\rho_{ABC},\mathcal R_{B\to BC}(\rho_{AB}))$, with $\mathcal R$ depending only on $\rho_{BC}$ ($\sigma=\mathrm{id}_A\otimes\rho_{BC}$, $\mathcal N=\mathrm{Tr}_C$).

**Proof.**
- **Lemma 3.1:** the Rényi difference $\tilde\Delta_\alpha\to D(\rho\|\sigma)-D(\mathcal N\rho\|\mathcal N\sigma)$ as $\alpha\to1$, and $\tilde\Delta_{1/2}=-2\log F(\rho,\mathcal P\circ\mathcal N(\rho))$.
- **Lemma 3.2:** Hirschman's strengthening of the Hadamard three-line theorem, applied to $G(z)=([\mathcal N\rho]^{z/2}[\mathcal N\sigma]^{-z/2}\otimes\mathrm{id}_E)U\sigma^{z/2}\rho^{1/2}$ on the strip.
- **The kernel:** $\beta_\theta(t)=\frac{\sin\pi\theta}{2\theta(\cosh\pi t+\cos\pi\theta)}\to\beta_0$ as $\theta\searrow0$.
- **Checked (§15 H7):** on random 3-qubit states, $I-(-2\log_2F)\ge0.145$ bits, with $\mathrm{Tr}\,\mathcal R(\rho_{AB})=1.0000$.

*Mine (KNOWN-LINK, NCG).* $X\mapsto\sigma^{it}X\sigma^{-it}$ is the modular automorphism group of the faithful state $\sigma$ on a matrix algebra (sign conventions *from memory*). So $\mathcal R^t$ is the Petz map intertwined with the **modular flows** of $\sigma$ and $\mathcal N(\sigma)$, and the universal recovery map is a **modular-time average** of the Petz map against the strip kernel $\beta_0(t)=\frac\pi4\mathrm{sech}^2(\pi t/2)$. Compare the KMS-strip kernels $1/|\sinh\pi t|$ in [Y] and $1/\cosh(2\pi t/\beta)$ in [CR] ([arxiv-2609.38007.md](arxiv-2609.38007.md) B14).

### 8.3 Small CMI does not mean close to an exact quantum Markov chain (Source: [ILW])

- **Theorem 4:** $\Delta(\rho):=\inf\{D(\rho\|\mu):\mu\text{ Markov}\}\ge I(A:C|B)_\rho$. **Theorem 2:** the infimum is attained, with $\dim\hat B\le d_B^4$.
- **Example 6:** $|\psi(x)\rangle$ on three qubits with $\Delta(\rho)/I\ge-\ln2/\ln x+O(1)\to\infty$ as $x\to0$.
- **Example 7:** the purification of the maximally mixed symmetric state has $I<1$ bit but $\Delta\ge\log d$, showing "that a log-dimensional factor is also necessary."
- Conclusion: "the characterisation of quantum Markov chains in terms of vanishing quantum conditional mutual information is not robust."

[FR] resolve this by measuring closeness through *recovery*, with no Markov target.

---

## 9. "Approximate quantum Markov networks are thermal", and decay of CMI

### 9.1 Kato–Brandão (Source: [KB])

**Definition.** $\rho_{A_1\dots A_m}$ is a quantum $\varepsilon$-approximate Markov chain if $I(A_1\dots A_{i-1}:A_{i+1}\dots A_m\mid A_i)\le\varepsilon$ for all $i$.

**Theorem 1.** "Let $\rho_{A_1\dots A_m}$ be a quantum $\varepsilon$-approximate Markov chain on a 1D open chain … Then there exists a short-range Hamiltonian $H=\sum_{i=1}^{m-1}h_{A_i,A_{i+1}}$ … such that $S(\rho\|e^{-H}/Z)\le\varepsilon m$." "Note that Theorem 1 is not restricted to full-rank states."

*Proof in full.*
1. Let $\sigma$ be the maximum-entropy state with $\sigma_{A_iA_{i+1}}=\rho_{A_iA_{i+1}}$. For any $\omega$ in the local Gibbs family $\mathcal E(\mathcal X)$, the Pythagorean theorem gives $S(\rho\|\omega)=S(\rho\|\sigma)+S(\sigma\|\omega)$ (Weis). With $\omega$ maximally mixed, $S(\rho\|\sigma)=S(\sigma)-S(\rho)$.
2. SSA, iterated: $S(\sigma)\le\sum_{i=1}^{m-2}[S(A_iA_{i+1})-S(A_{i+1})]+S(A_{m-1}A_m)$. Every term is a two-site entropy, so it can be evaluated on $\rho$.
3. Approximate Markov, telescoped: that sum is $\le S(\rho)+(m-1)\varepsilon$.
4. $\sigma$ lies in the closure of $\mathcal E(\mathcal X)$, so some local Gibbs $\omega$ has $S(\sigma\|\omega)\le\varepsilon$.
5. Pythagoras: $S(\rho\|\omega)\le m\varepsilon$. $\square$

Only SSA and maximum entropy are used. No positivity, no commutation, no locality estimate.

**Theorem 2 (closed chain).** The same bound with three-site terms, under either:
- (i) $I(A_i:A\setminus A_{i-1}A_iA_{i+1})\le\varepsilon$, or
- (ii) a uniform Markov property ("Tr$_{A_i}(\rho)$ is a quantum $\varepsilon$-approximate Markov chain").

**Theorem 3.** If the state with any one block traced out is $\varepsilon$-approximately Markov, then for $K=\Theta(N)$
$$\min_{\mu\in\mathcal E^K_{nn}}S(\rho_X\|\mu)=I(A:C|B)_\rho+\epsilon(N,\delta),\qquad|\epsilon|\le cN^{5/2}\delta^{1/16},\ \delta=8\sqrt\varepsilon+2^{-N}.$$
"the value of the conditional mutual information for the whole system approximately represents the distance from the set of local Gibbs states."

**Theorem 4 (converse).** "Let $H=\sum_ih_i$ be a short-range 1D Hamiltonian with $\|h_i\|\le1$ … For an inverse temperature $\beta>0$ and any partition $ABC$ with $d(A,C)\ge l_0$, there exists a CPTP-map $\Lambda_{B\to BC}$ … $\|\rho^H_{ABC}-\Lambda_{B\to BC}(\rho^H_{AB})\|_1\le e^{-q(\beta)\sqrt{d(A,C)}}$, where $q(\beta)=ce^{-c'\beta}$ if the correlation length of $\rho$ is $\xi=e^{O(\beta)}$."

**Corollary 5.** As transcribed from the extraction:
$$I(A:C|B)\le6\Big(d(A,C)+\frac{8(1+\frac{q(\beta)}2\sqrt{d(A,C)})}{q(\beta)^2}\Big)e^{-\frac{q(\beta)}2\sqrt{d(A,C)}} .$$

**Corollary 6.** Depth-two circuits with gates on $O(e^{O(\beta)}\log^2(n/\varepsilon))$ sites prepare 1D Gibbs states.

*Proof architecture of Thm 4.*
- Quantum belief propagation: $\frac{d}{ds}e^{-\beta H(s)}=-\frac\beta2\{e^{-\beta H(s)},\Phi^{H(s)}_\beta(V)\}$, their (38). It makes the effect of a local change of $H$ quasi-local.
- Lemma 9: a CP, trace-non-increasing approximate recovery with error $C_1(\beta)e^{-q_1(\beta)l}$.
- A repeat-until-success channel (their (60)): split $B$ into blocks $B_i$ ($|B_i|=2l$) and buffers $\bar B_i$ ($|\bar B_i|=l$), so that $d(A,C)=\tfrac32l^2-l$. Failed branches are traced out together with their buffer.
- 1D exponential decay of correlations (Araki) bounds the effect of failures. This is where $\xi=e^{O(\beta)}$, and hence $q(\beta)$, enters. The $\sqrt d$ comes from $l\sim\sqrt d$.
- [KB] §V: in 2D "the success probability decays too rapidly, and therefore our strategy does not work."

**Conjecture 1 of [KB].** "$I(A:C|B)_\rho\le Ce^{-cd(A,C)}$" in $D$ dimensions. On extending Theorem 1 to higher dimensions: "we may need additional conditions for general graphs. Although we do not know any counter-example, we also could not find any partial result."

### 9.2 High temperature, cluster expansion (Source: [KKB])

- **Setting.** $H=\sum_{|X|\le k}h_X$ with $\sum_{X\ni v,\mathrm{diam}X\ge R}\|h_X\|\le f(R)$, $f(1)\le1$.
- **Theorem 1.** For finite range $r$ and $\beta<\beta_c:=1/(8e^3k)$, "the Gibbs state $\rho$ is an approximate Markov network on an arbitrary subset $V_0\subseteq V$" with
  $$I_\rho(A:C|B)\le e\,\min(|\partial A_r|,|\partial C_r|)\frac{(\beta/\beta_c)^{d_{A,C}/r}}{1-\beta/\beta_c}.$$
- **Theorem 2.** Quasi-locality of the effective Hamiltonian: $\|\Phi_L-\Phi_{\partial L_l}\|\le\frac{e^{4\beta}(\beta/\beta_c)^{l/r}}{1-\beta/\beta_c}|\partial L_r|$.
- **Theorem 3 (power law $f(R)=R^{-\alpha}$).** For $\beta<\beta_c/11$, $d_{A,C}\ge2\alpha$: $I\le\beta\min(|A|,|C|)C_\beta/d^\alpha_{A,C}$, with $C_\beta=\frac{11e^{1/k}/\beta_c}{1-11\beta/\beta_c}$.
- **Method.** A "generalized cluster expansion" of $\tilde H_{\vec a}(A:C|B)=\log\tilde\rho^{AB}_{\vec a}+\log\tilde\rho^{BC}_{\vec a}-\log\tilde\rho^{ABC}_{\vec a}-\log\tilde\rho^B_{\vec a}$ in the coupling parameters. Prop. 6: only connected clusters linking $A$ and $C$ survive.
- **Caveat (Source: [CR] footnote 1).** "the proof requires expansions of operator-valued partial trace functionals, whose correctness remains unclear." Kato–Kuwahara (abstract, sibling digest) point to "intrinsic divergence problems rooted in the Baker–Campbell–Hausdorff formula".

### 9.3 Any temperature: pairwise in $D$ dimensions, global in 1D (Source: [Kuw])

- **Theorem (Supp. Thms 4–5).** $G_I(R)=D_{AC}\,e^{-c_1R/(\beta^{D+1}\log R)+c_2\log(\beta|AC|)}$; in 1D, $G_I(R)=e^{-c_3R/\beta+c_4\beta\log(\beta R)}$.
- "For $D=1$, the result does not depend on the subsystem sizes … We thus conclude that the global Markov property holds in one-dimensional quantum Gibbs states at arbitrary temperatures."
- "for $D\ge2$ … the current bound is insufficient to prove the global Markov property due to the growth of the coefficient $D_{AC}$, which increases as $e^{\Omega(|A|+|C|)}$." Hence only the pairwise property.
- **Correlation length.** The main text says $O(\beta)$ in 1D and $O(\beta^{D+1})$ for $D\ge2$. The Remark after Supp. Thm 5 says "correlation length of $\tilde O(\beta^2)$" in 1D, the $e^{\Theta(\beta)\log R}$ prefactor moving the onset. Both are as stated; I have not reconciled them.
- "much smaller than that of the bipartite correlation function, which can be infinitely large at critical points in high dimensions or at least exponentially increases with $\beta$ as $e^{O(\beta)}$ in one dimension."
- **Method.** Effective ("entanglement") Hamiltonians on subsystems; Lemma 16 / Cor. 17 bound $I$ via a regularised partial-trace operator $P_{L,\tau}=e^{-\tau Q_L}$, with $\tau=R/\Theta(\beta^{D+1}\log R)$; belief-propagation operators in 1D.

### 9.4 Local Markov at any temperature, bounded degree (Source: [CR]; full treatment in [arxiv-2609.38007.md](arxiv-2609.38007.md) §8)

- **Table I of [CR].**
  - Pairwise: $\exp(c|AC|-\mathrm{dist}/\xi)$ [Kuw24].
  - Local: $|A||C|\exp(c\min(|A|,|C|)-\mathrm{dist}/\xi)$ (their Thm III.1).
  - Global: $0$ "if $\mathrm{dist}(A,C)\ge1$" for commuting or classical Hamiltonians; $\exp(-\mathrm{dist}/\xi)$ in 1D [KB19, Kuw24].
- **Recovery map.** $\mathcal R_{A,t}=\frac1t\int_0^t\exp\big(s\sum_{a\in\mathcal P^1_A}\mathcal L_a\big)ds$, with all single-site Paulis on $A$. "In fact, the recovery map is a time-averaged detailed-balanced Lindbladian based on single-Pauli jumps on A."
- **Corollary III.2.** For $\mathrm{dist}(A,C)\ge4e^2\beta d$: $I\lesssim r'|A||C|\exp(\mu'\min(|A|,|C|)-\lambda'\mathrm{dist}(A,C))$.
- **Key steps.** $\mathcal E_A(\mathcal R^\dagger_{A,t}[X])\le2/t$, "unconditional polynomial decay … under time-averaging … independent of the spectral gap". Then $\|\mathcal R^\dagger_{A,t}[X]-(\mathcal R^\dagger_{A,t}[X])_{-A}\|\le e^{\mu|A|}\mathcal E_A^\nu\le re^{\mu|A|}/t^\nu$.
- **Remark III.1.2.** "The exponential dependence on |A| is hard to remove unconditionally … If the present argument can be combined with a faster mixing time or spectral gap analysis, one might be able to improve the exponential dependence on |A|, hence establishing the global Markov property."
- **[Y]** replaces $e^{c|A|}$ by $e^{Cg_{A|A^c}}$ (cut strength) with a static, modular-cocycle proof; see the sibling digest.

### 9.5 Global Markov at high temperature via a quantum Dobrushin condition (Source: [BLMT])

- **Definition 5.1.** For $\Phi=\Phi_1+\dots+\Phi_n$, the influence matrix is $D^{(\Phi)}_{ij}=\max_{\rho,\sigma\ j\text{-neighbours}}\|\Phi_i(\rho-\sigma)\|_{W_1}$, using the quantum Wasserstein-1 norm of De Palma et al. The condition is $\|D^{(\Phi)}\|_{1\to1}\le1-\gamma/n$.
- **Lemma 5.2.** This gives a unique fixed point and $\|\Phi^\tau(\rho)-\sigma\|_{\rm tr}\le\varepsilon\|\rho-\sigma\|_{\rm tr}$ for $\tau\ge\frac n\gamma\log\frac n\varepsilon$ (Banach, with $\frac12\|X\|_{\rm tr}\le\|X\|_{W_1}\le\frac n2\|X\|_{\rm tr}$).
- **Theorem 2.1.** For $\beta<1/(10000K^3b^2d)$, Lindbladian evolution for $O(\log(n/\varepsilon))$ time mixes. "quantum analog of Weitz's result."
- **Theorem 2.3 ("globally Markovian").** For $\beta<1/(Ce^{16a}bK^2(d+1))$, $I_\sigma(A:C|B)=O(a^a|A|\cdot|C|)\exp(-\mathrm{dist}(A,C)/\zeta)$.
- **Method.**
  - Fact 6.1: $I\le O(|A|)\|\sigma-\rho\|^{1/2}_{\rm tr}$ for the recovered $\rho$.
  - Recovery runs the Lindbladian restricted to $C_\Delta$: $\|e^{\mathcal L_{\rm CMI}t}(\sigma_{AB}\otimes\rho_C)-\sigma_{ABC}\|\le e^{-t/2}|C|+e^{2t}10^{-\Delta}$, with $t=\mathrm{dist}/4$, $\Delta=\mathrm{dist}/8$.
  - "disagreements cannot teleport across the system" (a light cone in transport plans).
- Classically, "strong temporal mixing … implies strong spatial mixing [Wei04]"; "Weitz showed that strong spatial mixing and strong temporal mixing are equivalent."

### 9.6 Metastability and the strong Markov property (Source: [C26])

- **Theorem II.1 ([BCV25]).** Approximate detailed balance for all single-site Paulis on $A$ gives $\|\sigma-\mathcal R_{A,t}[\mathcal N_A[\sigma]]\|_1\le e^{\mu|A|}t^{-\lambda}+c|A|t\epsilon_{ADB}^{1/2}$.
- **Strong local Markov (Def. I.4).** $\|\mathcal M_{AB}[K\sigma K^\dagger]-\sigma\,\mathrm{Tr}[K\sigma K^\dagger]\|_1\le\epsilon$ for all $\|K\|\le1$ on $A$: recovery of every post-selected branch.
- **Lemma I.1 / III.1.** Clustering between $A$ and $C$ plus approximate detailed balance on $AB$ ⇒ strong Markov, with error $\epsilon_{AC}+e^{\mu|AB|}(t^{-1}+\sqrt{\epsilon_{ADB}})^\lambda+c|AB|t\sqrt{\epsilon_{ADB}}$.
- **Lemma I.2 / III.2.** Strong Markov (with any recovery map) ⇒ $4\epsilon$-clustering.
- **Consequences.**
  - Repeatability (Lemma IV.1): single-copy tomography with error $2e^{-2r\tau^2}+r\epsilon$.
  - "locally very close or locally very far" (Lemma IV.4).
  - Local extremality (Lemma IV.5): $\|\sigma^A_1-\sigma^A_2\|_1\le2\sqrt{2\epsilon/(p_1p_2)}$, "reminiscent of quantum error-correcting codes".

### 9.7 Constants at a glance

| result | geometry, temperature | bound on $I(A:C\vert B)$ | size dependence |
|---|---|---|---|
| [KB] Cor. 5 | 1D, any $\beta$ | $\sim d\,e^{-q(\beta)\sqrt d/2}$, $q=ce^{-c'\beta}$ | none (global) |
| [Kuw] | 1D, any $\beta$ | $e^{-c_3R/\beta+c_4\beta\log(\beta R)}$ | none (global) |
| [Kuw] | $D$-dim, any $\beta$ | $D_{AC}e^{-c_1R/(\beta^{D+1}\log R)+c_2\log(\beta\vert AC\vert)}$ | $e^{\Theta(\vert A\vert+\vert C\vert)}$ (pairwise) |
| [CR] | bounded degree, any $\beta$ | $\vert A\vert\vert C\vert e^{\mu'\min(\vert A\vert,\vert C\vert)-\lambda'\mathrm{dist}}$ | exp. in the smaller region (local) |
| [Y] | bounded degree, any $\beta$ | $C_\beta e^{C_\beta g_{A\vert A^c}-c_\beta r}$ | exp. in cut strength |
| [KKB] | any graph, $\beta<1/(8e^3k)$ | $e\min(\vert\partial A_r\vert,\vert\partial C_r\vert)(\beta/\beta_c)^{d/r}/(1-\beta/\beta_c)$ | boundary (caveat §9.2) |
| [BLMT] | lattice, $\beta<1/(Ce^{16a}bK^2(d+1))$ | $O(a^a\vert A\vert\vert C\vert)e^{-\mathrm{dist}/\zeta}$ | polynomial (global) |
| commuting / classical | any graph, any $\beta$ | $0$ beyond the interaction range | — |

---

## 10. Ledger: which implications hold exactly, approximately, or fail

| # | statement | status | source |
|---|---|---|---|
| 1 | classical, positive: (P)⇔(L)⇔(G)⇔(F), canonical potential unique | **exact** | [HC] via [Cl90] Thms 1–2 (§2.0), [GL] Thm 3.1, [Yu], §3 |
| 2 | classical, any: (F)⇒(G)⇒(L)⇒(P) | **exact** | [Yu], §3.4 |
| 3 | classical, no positivity: (L)⇒(G) | **fails** (intersection) | [Yu] 5-variable chain; [Pe] |
| 4 | classical, no positivity: (G)⇒(F) | **fails**: support shape or local cohomology | [Mo] via [GL] Ex. 3.3; [GL] Lemma 5.2 |
| 5 | classical, chordal graph: (G)⇔(F) with same-graph zeros | **exact** | [GL] Thm 3.2 (Lauritzen) |
| 6 | classical, safe symbol, any countable graph: Markov cocycles = Gibbs cocycles | **exact** | [CM] Thms 3.1–3.2 |
| 7 | classical, $X_r$ ($\mathbb Z^d$ colourings): $\mathbf M=\mathbf G$ but $\dim\mathbf M^\sigma=r>r-1=\dim\mathbf G^\sigma$; invariant MRFs nonetheless Gibbs | **mixed** | [CM] Props 5.2–5.5, Thm 6.1 |
| 8 | classical, $\mathbb Z^2$ without pivot property: invariant MRF not Gibbs for any invariant finite-range interaction | **fails** | [CM] §9 |
| 9 | classical, $\mathbb Z$, shift-invariant, finite alphabet: every MRF a Markov chain and Gibbs | **exact**; fails without invariance or for countable alphabet | [CM] §1 |
| 10 | classical, approximate: $I=\min_{\rm Markov}D(P\Vert Q)$; chain maxent bound $\le m\varepsilon$; pointwise $\Vert\Phi_A\Vert\le2^{\vert A\vert-2}\delta$ | **exact identity / approximate** | [ILW]; [KB] Thm 1; §5.3 (Mine) |
| 11 | quantum, full rank: pairwise ⇔ local ⇔ global (exact) | **exact** (2026) | [LZ] Prop. 2.5 |
| 12 | quantum, full rank: Markov network ⇒ $\log\rho$ clique-local | **exact** | [LP] Thm 4.7; [BP] Thm 2 |
| 13 | quantum: commuting clique $H$ ⇒ Markov network (also commuting projectors) | **exact** | [BP] Thm 3 |
| 14 | quantum, full rank: Markov ⇒ commuting clique $H$ | **exact** on chains, trees, triangle-free graphs; **fails** with triangles | [BP] §II.A, Thms 4–5, Fig. 2 |
| 15 | quantum: Markov ⇔ Gibbs of a shield-commuting $H$ | **exact** | [BP] Def. 1 with (4) |
| 16 | quantum: local (non-commuting) Gibbs ⇒ exact Markov | **fails** | [LP] Ex. 4.8; [Kuw]; [C26] |
| 17 | quantum: local Gibbs ⇒ approximately Markov | 1D global any $\beta$; $D$-dim pairwise ([Kuw]) and local ([CR], [Y]) any $\beta$; global at high $T$ ([BLMT]; [KKB] caveat); **global at low $T$, $D\ge2$: open** | §9 |
| 18 | quantum: approximately Markov ⇒ approximately local Gibbs | 1D open chain **holds** ($\le m\varepsilon$); closed chain under extra assumptions; $D\ge2$ **open** | [KB] Thms 1–3 |
| 19 | quantum: small CMI ⇒ close to an exact quantum Markov chain | **fails** (dimension factors) | [ILW] |
| 20 | quantum: small CMI ⇔ good recovery | **holds**: $2^{-I/2}\le F$; converse with $\log\dim A$ | [FR]; [JRSWW] |
| 21 | quantum: strong (post-selected) Markov ⇔ clustering (given metastability) | **holds**, both directions | [C26] |
| 22 | quantum: Markov completion of overlapping marginals | exists iff $\mathrm{Tr}\,T(\mathcal R)=1$ (two cliques; chordal); maxent completion need not be Markov | [LZ] Thms 3.1, 3.10, Rem. 3.15 |

---

## 11. Proof architectures side by side: where each hypothesis is used

| proof | hypothesis | used in | key step |
|---|---|---|---|
| [HC] original (§3.8; [Cl90]) | positivity | $\log P$ and the blackened ratios (2.8) exist | $M(X)\iff\beta_X\log P=\log P$; $I(\beta)=\bigcap_zI(\beta_z)\subseteq I(\beta_X)$ for commuting idempotents |
| classical HC (§3) | safe symbol / positivity | Lemmas A–B: all vacuum substitutions in the support | cross-ratio identity $L_B=0$ (exact) |
| | pairwise CI only | Lemma B | — |
| [BP] Thm 2 | full rank | $\log\rho$ exists; (4) | partial-trace cumulants; $\mathrm{Tr}(HK_X)=\mathrm{Tr}K_X^2=0$ |
| [BP] Thm 4 | triangle-free | (24), (26) | Lemma 1 (genuine-operator commutators) |
| [LZ] intersection | strict positivity | unique maximiser $X_t=(t+\Delta_{\tau,\rho})^{-1}I$ | intersection of tensor subalgebras |
| [KB] Thm 1 | none beyond finite dimension | SSA telescoping; maxent duality | Pythagorean theorem for $\mathcal E(\mathcal X)$ |
| [KB] Thm 4 | 1D; finite range | QBP; Araki clustering; repeat-until-success | $e^{-q\sqrt d}$ from $d\sim l^2$ |
| [FR] | none | typicality, de Finetti | fidelity of permutation-invariant states; $\mathrm{poly}(n)^{1/n}\to1$ |
| [JRSWW] | none | Rényi difference at $\alpha=\frac12$ and $\alpha\to1$ | Hirschman's three-line theorem |
| [KKB] | $\beta<\beta_c$ | absolute convergence of the cluster expansion | connected clusters linking $A$, $C$ |
| [Kuw] | finite range; any $\beta$ | effective-Hamiltonian locality; PTP operator | $\tau$ tuned to $R/(\beta^{D+1}\log R)$ |
| [CR] | bounded degree; any $\beta$ | Lieb–Robinson; KMS detailed balance | $\mathcal E_A\le2/t$ (time averaging, gap-free); $e^{\mu\vert A\vert}$ from the Pauli expansion |
| [BLMT] | $\beta<\beta_c$ | quantum Dobrushin (contraction in $W_1$) | path coupling; light cone of disagreements |

---

## 12. What is genuinely new

**In the literature after 1971 (Source).**
1. **Specifications as cocycles** ([CM]): Markov specifications are locally determined Radon–Nikodym cocycles on the homoclinic relation. Hammersley–Clifford becomes $\mathbf M_X=\mathbf G_X$, and the failure has a dimension ($\mathbf M/\mathbf G$).
   - First explicit computations: $X_r$ has a one-dimensional invariant non-Gibbs part given by the height cocycle, which only frozen invariant measures can carry.
   - [Ch] shows $\mathbf M/\mathbf G$ is invariant under folding.
2. **The two readings of Moussouris** ([GL]): the uniform example is Gibbs with global constraints. The genuine failure is "no strictly positive Markov extension", a chain of cross-ratio relations that closes inconsistently.
3. **Quantum conditional independence** ([HJPW]): it is Petz sufficiency, with the direct-sum/tensor structure from Koashi–Imoto.
   - [BP] reduce quantum Hammersley–Clifford to commutation, and show triangle-free graphs suffice and triangles break it.
   - [LZ] (2026) close the intersection problem and characterise Markov completions by $\mathrm{Tr}\,T(\mathcal R)=1$.
4. **Robust versions:**
   - in distance, they fail ([ILW]);
   - in recovery, they hold ([FR]), with a universal, explicit, modular-averaged Petz map ([JRSWW]);
   - approximately Markov 1D states are approximately Gibbs ([KB]).
5. **CMI decay** (§9), now proved:
   - globally in 1D and at high temperature;
   - pairwise, locally and with boundary dependence in $D$ dimensions at every temperature.

   [C26]: the strong (post-selected) Markov property is equivalent to clustering.

**In this digest (Mine; labelled where stated).**
1. **The integrability form** (§3.7): Hammersley–Clifford ⇔ "vanishing mixed differences across non-edges ⇒ clique decomposition". The vacuum Möbius inversion, the spin-product inversion and the quantum cumulants are one construction with three choices of conditional expectation: evaluation, product measure, trace. [LP]'s pure-vector vacuum is a fourth. *Review note:* the operator form itself is [HC]'s 1971 blackening algebra (§3.8), so it is not new here. What is added is the reading of the pairwise property as $\Delta_i\Delta_j\log P=0$.
2. **Two classical obstructions** (§4.4): support shape and local cohomology. Moussouris is purely the first and [GL]'s $P^*$ purely the second. $P^*$'s 3-body term is explained by edge contraction ($x_3=x_4$). Checked numerically.
3. **Where the quantum converse fails** (§7.5): joint-log locality versus marginal-log locality. CMI is the mixed second difference of the log-marginals, which is what [KKB] expand.
4. **Classical intersection as connectivity** (§4.2): positivity ⇒ intersection is irreducibility of the two-block Gibbs sampler. Made quantitative in §13 B3 with constant $1/(1-c_F)$; checked numerically.
5. **Pinning is free, marginalising fills in** (§5.4, §7.6). This is the classical source of the hereditary class used by spectral independence. The quantum analogue of pinning is the strong Markov property, which costs clustering (§13 B5).

---

## 13. Testing the user's intuition (all Mine unless a source is named)

### 13.1 Verdict

The intuition is **right about one shared skeleton and wrong about the engine.**

- **Wrong about the engine.** Hammersley–Clifford contains no expansion, no spectral gap and no rate.
  - It is an exact, local, combinatorial integrability theorem (§3.7).
  - Its one hypothesis, positivity, is topological: the support must contain the Boolean lattices of vacuum substitutions. That makes the support's move graph connected (intersection) and free of holonomy (factorisation).
  - Markov structure holds whatever the mixing behaviour, including where mixing fails: at phase transitions, and on expanders at low temperature ([Y], [CM]; [arxiv-2609.38007.md](arxiv-2609.38007.md) B1–B2; [expanders.md](expanders.md) §9.7).
- **Right about the skeleton.** Both theories are written in the algebra of **commuting local conditional expectations** and their Möbius (Efron–Stein) decomposition:
  - Hammersley–Clifford fixes the *support* of the decomposition of $\log P$ (complete sets only);
  - expansion and mixing bound the *size* of what local conditional expectations fail to remove (contraction of products such as $E_{\neg B}E_{\neg D}$).

  The two meet exactly where a qualitative statement must become quantitative:
  - approximate pairwise ⇒ approximate global Markov (B3);
  - local ⇒ global quantum Markov, which the sources prove by mixing or contraction ([CR] App. B; [BLMT]);
  - strong Markov ⇔ clustering ([C26]).
- **Noncommutative geometry** enters in three places, all at the level of operator algebras and modular theory, not spectral triples:
  - Gibbs measures are KMS states of groupoid cocycles, and Markov specifications are locally determined cocycles (B6);
  - quantum Markov = Petz sufficiency; the universal recovery is a modular-time average; the 2026 intersection theorem is a relative-modular-operator argument (B7);
  - the quantum obstruction is algebraic: incompatible Koashi–Imoto splittings (§7.4).

  *Added in review:* the original 1971 proof is itself operator-algebraic. [HC] prove local ⇒ global inside a commutative algebra of idempotent composition operators, and the Markov property is a fixed-point equation for one of them (§3.8, B11).

  No retrieved Hammersley–Clifford result uses a Dirac operator or a spectral distance (B10).

### 13.2 Bridges

**B1 — THEOREM (Source [Pr]).** *Markov fields are exactly the stationary laws of reversible nearest-neighbour dynamics.* Every MRF on a finite graph is the equilibrium of a time-reversible birth/death process with nearest-neighbour interactions, and conversely. The *existence* of a local reversible dynamics is equivalent to Markov; its *gap* is the expander-type quantity and is not constrained by Hammersley–Clifford. Quantum descendants: [CR] (recovery by a KMS-detailed-balanced Lindbladian) and [C26] Thm II.1 (approximate detailed balance ⇒ local Markov).

**B2 — KNOWN-LINK (Sources [HJPW] Lemma 12, [CR] §IV.C).** *The same Cesàro average produces the recovery in the exact and approximate cases.*
- [HJPW]: the mean-ergodic average $\lim_N\frac1N\sum(F^*)^n$ of the recovery channel is the conditional expectation onto its fixed-point algebra, whose structure *is* the Markov decomposition.
- [CR]: the time average $\frac1t\int_0^te^{s\mathcal L}ds$ of a detailed-balanced Lindbladian gives a recovery map with Dirichlet form $\le2/t$, "independent of the spectral gap".

Both use averaging, not a gap. A gap upgrades the polynomial rate $t^{-\lambda}$ to exponential ([CR] App. B). That is the expander contribution, and it is where $e^{|A|}$ would go.

**B3 — THEOREM (elementary, Mine; checked §15 H8).** *Quantitative intersection is controlled by the Friedrichs angle, i.e. the spectral gap of a two-block Gibbs sampler.* Let $h(b,c,d)=P(A=1\mid b,c,d)$ in $L^2(P_{BCD})$. Write $E_{\neg B}=E[\cdot\mid C,D]$ and $E_{\neg D}=E[\cdot\mid B,C]$ (each integrates out one block), $T=E_{\neg B}E_{\neg D}$, and $c_F=\|T|_{L^2(C)^\perp}\|$. Then
$$\|h-E[h\mid C]\|\le\frac{\|h-E_{\neg B}h\|+\|h-E_{\neg D}h\|}{1-c_F}.$$
*Proof:*
1. $(I-T)h=(h-E_{\neg B}h)+E_{\neg B}(h-E_{\neg D}h)$.
2. $g=h-E[h|C]$ is orthogonal to $L^2(C)$ and $(I-T)g=(I-T)h$.
3. $\|(I-T)g\|\ge(1-c_F)\|g\|$.

So "$A\perp B|CD$ and $A\perp D|BC$ approximately" ⇒ "$A\perp BD|C$ approximately", at the cost $1/(1-c_F)$. Here $c_F$ is the cosine of the angle between $L^2(C,D)$ and $L^2(B,C)$, a correlation/expansion constant of the two-block sampler.
- Positivity makes $c_F<1$ (irreducible sampler; exact intersection).
- The support $\{X_A=X_B=X_D\}$ gives $c_F=1$ with both hypotheses exact and $\|h-E[h|C]\|=0.5$ (H8).
- Over 200 random positive laws the inequality's ratio is at most $0.669$.

This is the precise point where Hammersley–Clifford's positivity and expander-type spectral quantities are the *same* hypothesis at two resolutions: irreducibility versus a quantitative gap. In entropy form the analogue is approximate tensorization ([arxiv-2609.38007.md](arxiv-2609.38007.md) B6). It matches [CR]'s route to the global Markov property ("uniform local gap"), and the [expanders.md](expanders.md) §11.2 angle inequality.

**B4 — ANALOGY → SPECULATION.** *Positivity as simple connectivity of a square complex; Hammersley–Clifford as vanishing of a local first cohomology.*
- Build a cube complex $K(X)$ on the support $X$:
  - vertices: the configurations;
  - edges: pivot moves (single-site changes inside $X$);
  - 2-cells: commuting squares $x,x^i,x^j,x^{ij}\in X$ for distinct sites.
- A specification is a 1-cochain $c$ on pivot edges whose value at site $i$ depends on $x_{i\cup\partial i}$ (locality).
- Consistency ([Be] §2) is closedness of $c$ around every cycle of the pivot graph.
- (P) is closedness around every square at a *non-edge*, which reads $\Delta_i\Delta_j\log P=0$ (§3.7).

Under positivity $K(X)$ is the full product complex, the cycle space is generated by squares, and Hammersley–Clifford says: **local closed cochains are coboundaries of local (clique) potentials**.
- Moussouris' support is an 8-cycle with **no** squares (checked), so cycles exist that no square fills.
- [CM]'s height cocycle on $X_r$ is a genuine non-local class, arising from the lift of $\mathbb Z_r$-valued maps to $\mathbb Z$.

This is the same species of object as prong 3's U3 ([local-to-global-unlocks.md](../../local-to-global-unlocks.md) §5): the square 2-complex of the layered quiver, with $\delta F$ evaluated on squares.
- *What must be true:* a precise definition of the local cochain complex on $K(X)$ for which $\mathbf M_X/\mathbf G_X\cong H^1_{\rm loc}(K(X))$. I have not proved this identification.
- If it holds, a *coboundary expansion* constant of $K(X)$ would be exactly a quantitative (approximate) Hammersley–Clifford for supports without a safe symbol.

**B5 — KNOWN-LINK + ANALOGY.** *Pinning heredity is free classically; quantumly it is the strong Markov property and costs clustering.*
- Classically, pinnings of Markov fields are Markov fields on induced subgraphs (§5.4). That is the hereditary class on which trickle-down and spectral independence run ([hdx-spectral-independence.md](hdx-spectral-independence.md) §3.2).
- Quantumly there is no pinning; [Y]: "a marginal with an arbitrary unobserved exterior is not silently identified with a local Gibbs state".
- The closest operational analogue is post-selection on a local measurement. Recovery for *every* branch is the strong Markov property, which [C26] prove equivalent to clustering.

So the quantum hereditary class is not free: it costs a decorrelation (mixing-type) hypothesis. This is the strongest precise sense in which "expander-type" input is *necessary* in quantum Hammersley–Clifford-type theory.

**B6 — KNOWN-LINK (NCG; Sources [CM] §3 and Note 1 §5.2).** *Gibbs measures are KMS states of locally determined groupoid cocycles.*
- [CM]: MRFs adapted to $X$ are the measures non-singular for $\Delta_X$, and $\mu$ is Gibbs for $\phi$ iff its Radon–Nikodym cocycle is $e^{M_\phi}$.
- Note 1 (Renault; Neshveyev Thm 1.3): KMS$_\beta$ states of $C^*(G)$ for a cocycle $c$ are quasi-invariant measures with Radon–Nikodym cocycle $e^{-\beta c}$ (plus isotropy traces, none for an equivalence relation).
- So, with $\Delta_X$ given its inductive-limit étale topology (*from memory*; AF for full shifts), **Gibbs measures for $\phi$ are KMS states of $C^*(\Delta_X)$ for the dynamics of $M_\phi$**, and Markov specifications are the locally determined cocycles.

Hammersley–Clifford asks whether a locally determined cocycle is a *formal coboundary of a local potential*: $M_\phi=\delta H$ with $H=\sum_W\phi(x|_W)$ an infinite sum, not a function. This sits next to Note 1 §6's criterion. "Face sufficient ⇔ $F$ a genuine coboundary ⇔ the measure is central/tracial" is the trivial-cocycle case: $M\equiv0$ gives tail-invariant measures, the traces. The two questions are two quotients of the same cocycle space:
- **genuine coboundaries:** sufficiency, centrality;
- **local-potential coboundaries:** Gibbsianity, finite memory.

**B7 — THEOREM-level links (Sources [HJPW], [JRSWW], [LZ], [Y]).** *Quantum Markov theory is modular theory.*
- Quantum CI is Petz sufficiency (§6.1).
- Universal recovery is the Petz map averaged over modular time with the strip kernel $\beta_0$ (§8.2).
- [LZ]'s intersection proof and [Y]'s CMI bound use the *same* variational/resolvent representation of a relative-entropy difference through the relative modular operator. [LZ] need the maximiser $X_t=(t+\Delta_{\tau,\rho})^{-1}I$ to lie *exactly* in a subalgebra (intersection); [Y] bounds its *leakage* (sibling digest B7).
- Combined with B3, an approximate quantum intersection theorem would need the Friedrichs angle between $L(ACD)\otimes1$ and $L(ABC)\otimes1$ *in the $\rho$-weighted (KMS) inner product*. That angle is $<1$ for product $\rho$ and is the natural "local gap" otherwise (SPECULATION).

**B8 — ANALOGY.** *Two time averages.*
- The classical heat-bath update of a block $A$ (forget $A$, resample from $P(\cdot\mid x_{\partial A})$) is an *exact* recovery map, and is the Markov property.
- [CR]'s map has the same shape: forget $A$ by the Pauli twirl (the tracial vacuum of §7.2), then re-thermalise $A$ by a Cesàro-averaged detailed-balanced Lindbladian.
- [JRSWW]'s map averages the Petz map over *modular* time.

Dissipative time (CR) versus modular time (JRSWW): both average a one-parameter group or semigroup tied to the state against a strip-type kernel. Whether one can be obtained from the other is open; I found no source relating them.

**B9 — KNOWN-LINK (Source [Kuw]; [expanders.md](expanders.md) §9.4 and §9.7 for the expander side).** *The CMI correlation length is polynomial in $\beta$ while the bipartite correlation length can diverge.* Conditional decorrelation (Markov) is controlled by locality alone. Unconditional decorrelation (clustering) is the expansion/mixing property and fails at criticality. This is the cleanest one-sentence separation of the two theories in the sources.

**B10 — SPECULATION (flagged as such).** *A spectral-triple reading.* §8.2's modular averaging and §7.2's tracial vacuum suggest a Dirichlet form $\mathcal E(X)=\sum_a\|[A_a,X]\|_\rho^2$ ([arxiv-2609.38007.md](arxiv-2609.38007.md) B9) as a Lipschitz seminorm, with the Markov property as "the seminorm restricted to $A$ is controlled by its restriction to $\partial A$". No retrieved Hammersley–Clifford result has this form. It is recorded only as the place where a Connes-distance formulation would have to start.

**B11 — THEOREM (Source [Cl90] Lemmas 4–5) + ANALOGY (Mine; added in review).** *[HC]'s local ⇒ global is an intersection-of-fixed-spaces theorem for commuting idempotents. The [CR] recovery map has the same shape with non-commuting generators, and that is where a gap enters.*
- **[HC] (§3.8).** $M(z)\iff\log P\in I(\beta_z)$. The $\beta_z$ commute, so $\beta=\prod_z\beta_z$ is the idempotent onto $\bigcap_zI(\beta_z)$, and one finite product reaches that intersection exactly. Lemma 4 then gives every $M(X)$.
- **[CR] (§9.4).** $\mathcal L_A=\sum_{a\in\mathcal P^1_A}\mathcal L_a$, each $\mathcal L_a$ KMS-detailed-balanced.
  - In the KMS inner product each $\mathcal L_a^\dagger$ is self-adjoint and $\le0$. So $\ker\mathcal L_A^\dagger=\bigcap_a\ker\mathcal L_a^\dagger$ (elementary: $\langle X,-\mathcal L_A^\dagger X\rangle=\sum_a\langle X,-\mathcal L_a^\dagger X\rangle$, a sum of non-negative terms).
  - The time average $\mathcal R^\dagger_{A,t}=\frac1t\int_0^te^{s\mathcal L_A^\dagger}ds$ tends to the projection onto that intersection (mean ergodic theorem, *from memory*).
  - [CR]'s bound $\mathcal E_A(\mathcal R^\dagger_{A,t}[X])\le2/t$ is a gap-free rate for this convergence.
- **Same shape, two differences.** In both, the global object is the intersection of single-site fixed spaces.
  - (i) [HC]'s idempotents act on $\log P$; [CR]'s generators act on observables and states.
  - (ii) [HC]'s commute, so one product lands exactly on the intersection. [CR]'s do not, so the intersection is reached only asymptotically: polynomially by averaging, and exponentially only with a local gap ([CR] Remark III.1.2; App. B).
  - [CR] App. B draws the same line at the level of Hamiltonians. Their local-gap condition "unconditionally holds for $(H,\{A_a\})$ being a commuting Hamiltonian and local jumps $\mathcal P^1_A$ on a region $A$ with a local gap independent of the global system size". But "For general noncommutative Hamiltonians (with local jumps), we do not know of any a priori bound on the local gap, even assuming high temperature." Their Cor. B.2: a uniform local gap $c|A|^{-c'}$ upgrades the local Markov bound to a global one, $\mathrm{Poly}(|A|,|C|)\exp(-\mathrm{dist}(A,C)/\xi)$.
- **Where expansion enters.** For two non-commuting orthogonal projections, the speed of alternating projections onto the intersection is set by the Friedrichs angle (von Neumann–Halperin, *from memory*). B3 is that statement for classical intersection. So the most precise meeting point found between Hammersley–Clifford and expander-type quantities is this. *Expansion is the speed at which non-commuting local projections reach the intersection that, in the commuting world, a single product of [HC]'s idempotents reaches exactly.*
- **The two halves of [CR]'s map.** Its forgetful step (replace $A$ by $\tau_A$) is a blackening of $A$ with the tracial vacuum. Its re-thermalising step corresponds classically to the heat-bath resampling of $X$, which is exactly the operation that $M(X)$ makes $\partial X$-local.

---

## 14. Messages to the other prongs (programme rule C1)

**To prong 1 (theory).**
1. **Adopt [CM]'s Markov-cocycle formalism on the path groupoid** (tail/homoclinic relation of the Bratteli path space, Note 1 §1.3, §3.4):
   - a *Markov cocycle of range $k$*: $M(\mu,\nu)$ for histories differing in finitely many layers is determined by the window $F\cup\partial_kF$;
   - a *Gibbs cocycle of range $k$*: comes from a potential on depth windows of length $k+1$;
   - and the quotient $\mathbf M/\mathbf G$.

   This answers mlp-bridge's question (a) in definitions, not in values. With positivity, depth is a path graph and Hammersley–Clifford with blocks holds, so a law with finite Markov length $k$ is a range-$k$ Gibbs law. Without positivity the path is still chordal, so (G) ⇔ (F) with zeros allowed ([GL] Thm 3.2), but (L) ⇒ (G) can fail (§4.2).

   It passes the drag test: the definitions are stated for any étale equivalence relation with a locality structure.
2. **Two different coboundary questions (B6).** Note 1 §6's "is $F$ a coboundary" (sufficiency) and Hammersley–Clifford's "is $M$ a local-potential coboundary" (finite memory) are distinct quotients of one cocycle space. Record both.
3. **The intersection theorem (§6.2) as a commuting-square statement.** For Note 1's algebras, check whether the face algebras of different layers inside the history algebra form commuting squares in the KMS inner product. If they do, "pairwise sufficiency ⇒ global sufficiency" holds. If not, B3 says the loss is $1/(1-c_F)$.
4. **Directed version.** The quiver is a layered DAG. Lauritzen–Zwiernik's companion "Bayesian networks of density operators" (arXiv:2607.27876; abstract via Exa only) states that for positive definite operators the ordered, local and global *directed* quantum Markov properties are equivalent. It also states that kernel-based ("extrinsic") constructions are order-independent "exactly for transitive DAGs". That is the precise result to read before putting coherences on histories.

**To prong 2 (bridge).**
1. **Memory as fill-in.** The measured 25–35 % memory of the face process ([mlp-bridge.md](../../mlp-bridge.md) §3.2) is what Hammersley–Clifford theory predicts for a *marginal*, or a coarse-graining, of a Markov system (§5.4, §7.6). The face sequence is a function of the network's hidden-state chain, and functions of Markov chains are Markov only under lumpability (Kemeny–Snell, *from memory*). Rule P2.1 reading: D2 (minimal sufficient sub-frame) *is* a search for a lumpable coarse-graining.
2. **Positivity fails in the realization.** The quiver has many forbidden transitions, so measure the *global* CMI given the full past as well as lag-$j$ pairwise quantities. Without positivity, pairwise or local zeros do not imply global ones (§4.2).
3. **Calibration.** The surrogate (P2.4) must respect the support (same zero transitions). Otherwise a support-shape effect of type §4.4(a) can masquerade as memory.

**To prong 3 (unlocks).** Add Hammersley–Clifford as a reduction in the §1–3 format.

| item | entry |
|---|---|
| global object | the joint law on $n$ variables, $\prod\vert S_v\vert$ numbers |
| local certificate | vanishing log cross-ratios across non-edges (pairwise CI) plus a safe symbol |
| explicit bound | exact equality $\log P=\sum_{C}\Phi_C$; approximately, $\Vert\Phi_A\Vert\le2^{\vert A\vert-2}\delta$ pointwise (§5.3), or $\le m\varepsilon$ in relative entropy on chains ([KB]) |
| cost | $\sum_C\prod_{v\in C}\vert S_v\vert$ parameters |

Guards:
- no rate and no gap: Hammersley–Clifford is the *hereditary-class* ingredient, not the *spectral* one;
- hereditary under pinning, not under marginalisation or quantum truncation;
- the quantum version needs commutation (classically free) and, for heredity, clustering (B5).

B3 is the bridge to the spectral ingredient: intersection degrades by the inverse block-sampler gap.

---

## 15. Numerical checks

Script: [check_hammersley_clifford.py](check_hammersley_clifford.py) (pure numpy). Output of `python3 check_hammersley_clifford.py`:

| id | statement checked | result |
|---|---|---|
| H1 | Lemma A/B: positive Gibbs law on $C_4$ ($q=3$) | $\max\vert\Phi_A\vert$ off complete sets $=6.2\cdot10^{-15}$; reconstruction error $5.3\cdot10^{-15}$; generic law: $4.70$ |
| H2 | §5.3 pointwise approximate Hammersley–Clifford | $\max\vert\Phi_A\vert/(2^{\vert A\vert-2}\max\vert\log\text{cross-ratio}\vert)=1.000$ (bound tight) |
| H3 | Moussouris support | CMIs of random laws on it $=0$; every edge pattern occurs; flip graph an 8-cycle with 0 commuting squares |
| H3 | [GL] Lemma 5.2 | $P^*$ pairwise CMIs $=0$; positive-Markov-extension residual $0$ (uniform Moussouris) vs $0.231$ ($P^*$) |
| H4 | [BP] wheel example | CMIs $\le1.3\cdot10^{-15}$; $\Vert[h_5,h_C]\Vert=11.3$; generic triangle $H$: $I(1{:}3\vert245)=0.116$ |
| H5 | [LP] Ex. 4.8 | $I(A{:}C\vert B)=0.082,0.315,0.412,0.415$ bits at $\beta=0.5,1,2,4$; cumulants $K_{AC},K_{ABC}$ of $\log\rho$ $\le4\cdot10^{-14}$; Ising chain $I=0$ |
| H6 | [LZ] Thm 3.1, Lemma 3.2 | $\max(\mathrm{Tr}\,T-1)=-1.3\cdot10^{-3}$; $\mathrm{Tr}\rho(\log\rho-\log T)=I$ to $3\cdot10^{-15}$; Markov state: $\mathrm{Tr}\,T=1$, $T=\rho$ to $10^{-15}$ |
| H7 | [JRSWW] Cor. 4.1 | $\min[I-(-2\log_2F)]=0.145$ bits over 5 random 3-qubit states; recovered trace $1.0000$ |
| H8 | §13 B3 quantitative intersection | ratio $\le0.669$ over 200 random positive laws; degenerate law $X_A=X_B=X_D$: $c_F=1.000$, hypotheses exact, $\Vert h-E[h\vert C]\Vert=0.500$ |

These are checks of algebra and of retrieved theorems on small instances, not experiments.

*Review checks (scratch scripts, not added to the script above; reported in §2.0 and §3.8):*
- [Cl90]'s $M(X)$ is equivalent to $\beta_X\log P=\log P$ for each $X$ on $C_4$, $q=3$.
- $\operatorname{rank}\beta=\dim\bigcap_zI(\beta_z)=25$, and $\beta_X\beta=\beta$.
- Theorem 2's $Q$ vanishes off cliques ($8.9\cdot10^{-15}$).
- [LP]'s pure-vacuum inversion on the [BP] wheel gives non-clique $K_U$ of size $4.9\cdot10^{-14}$.

---

## 16. Sources read for this digest

**Classical.**
- J. Besag, *Spatial interaction and the statistical analysis of lattice systems*, JRSS B 36 (1974) 192–236: pp. 192–195 (OCR via Exa of the JSTOR scan at stat.cmu.edu). In review: pp. 192–201, 203 and the discussion, pp. 228–235 (alphaXiv reader on the cise.ufl.edu scan; most equations not extracted).
- P. Clifford, *Markov random fields in statistics*, in G. Grimmett, D. Welsh (eds.), *Disorder in Physical Systems* (OUP 1990) 19–32: §§1–2 in full (added in review; Exa, PDF of the volume at lib.ysu.am).
- G. Grimmett, *A theorem about random fields*, Bull. LMS 5 (1973) 81–84: abstract and references (Exa).
- C. J. Preston, *Generalized Gibbs states and Markov random fields*, Adv. Appl. Prob. 5 (1973) 242–261: abstract (Cambridge Core via Exa).
- J. Moussouris, *Gibbs and Markov random systems with constraints*, J. Stat. Phys. 10 (1974) 11–33: abstract and keywords (ADS via Exa).
- D. Pollard, Yale Stat 251 handout (2004), full. D. Shah, MIT 6.438 Lecture 3 (2014), proof section. H. Yu, Helsinki slides (2010), full. Duke STA 345 notes (2010), §4.1.
- A. Gandolfi, P. Lenarda, *A note on Gibbs and Markov random fields with constraints and their moments*, MEMOCS 4 (2016): pp. 407–417.
- N. Chandgotia, T. Meyerovitch, *Markov random fields, Markov cocycles and the 3-colored chessboard*, [1305.0808](https://arxiv.org/abs/1305.0808) v3: full text.
- N. Chandgotia, *Generalisation of the Hammersley–Clifford theorem on bipartite graphs*, [1406.1849](https://arxiv.org/abs/1406.1849) v2: pp. 1–2, 4–5, 9–12, 14–15, 21, 26.
- J. Peters, *On the intersection property of conditional independence…*, [1403.0408](https://arxiv.org/abs/1403.0408): Exa highlights only.
- Wikipedia, "Hammersley–Clifford theorem": history paragraph (Exa). ORA catalogue record (Exa). M. Gasse, thesis (2017): one sentence (Exa).

**Quantum.**
- P. Hayden, R. Jozsa, D. Petz, A. Winter, [quant-ph/0304007](https://arxiv.org/abs/quant-ph/0304007): full.
- M. Leifer, D. Poulin, *Quantum graphical models and belief propagation*, [0708.1337](https://arxiv.org/abs/0708.1337): pp. 1, 3, 5, 8, 10, 13, 15–16, 18–23, 48–49.
- W. Brown, D. Poulin, *Quantum Markov networks and commuting Hamiltonians*, [1206.0755](https://arxiv.org/abs/1206.0755): full.
- S. Lauritzen, P. Zwiernik, *The (Markov) marginal problem for density operators*, [2605.19453](https://arxiv.org/abs/2605.19453) v1: pp. 1–3, 6–9, 13–14, 18, 22, 26–29. Their *Bayesian networks of density operators*, [2607.27876](https://arxiv.org/abs/2607.27876): Exa highlights only.
- O. Fawzi, R. Renner, [1410.0664](https://arxiv.org/abs/1410.0664): pp. 1–4, 15–19, 22, 29–31.
- M. Junge, R. Renner, D. Sutter, M. Wilde, A. Winter, [1509.07127](https://arxiv.org/abs/1509.07127): pp. 1–12, 16–20, 24.
- B. Ibinson, N. Linden, A. Winter, [quant-ph/0611057](https://arxiv.org/abs/quant-ph/0611057): full.
- K. Kato, F. Brandão, [1609.06636](https://arxiv.org/abs/1609.06636) v3: pp. 1–11, 13, 17, 24, 29, 31.
- T. Kuwahara, K. Kato, F. Brandão, [1910.09425](https://arxiv.org/abs/1910.09425) v2: pp. 1–4, 6, 10–11, 17, 21, 23.
- T. Kuwahara, [2407.05835](https://arxiv.org/abs/2407.05835) v3: pp. 1–3, 18, 24, 85–86, 89, 91, 94.
- C.-F. Chen, C. Rouzé, [2504.02208](https://arxiv.org/abs/2504.02208) v1: pp. 1–2, 4–5, 9–10, 13. More in the sibling digest.
- A. Bakshi, A. Liu, A. Moitra, E. Tang, [2510.08542](https://arxiv.org/abs/2510.08542) v1: pp. 1–5, 8–9, 12, 41, 43, 45–46.
- C.-F. Chen, *Note on strong quantum Markov properties*, [2605.02877](https://arxiv.org/abs/2605.02877) v1: full.
- T. H. Yang, [2609.38007](https://arxiv.org/abs/2609.38007): through [arxiv-2609.38007.md](arxiv-2609.38007.md).

**Found by search, not read beyond a snippet:**
- *Abstract Markov random fields* ([2407.02134](https://arxiv.org/abs/2407.02134)): "Markov random fields are known to be fully characterized by properties of their information diagrams";
- *Lattice supported distributions and graphical models* ([2411.03139](https://arxiv.org/abs/2411.03139)): "We prove a generalization of the Ham[mersley–Clifford theorem]";
- *Belavkin–Staszewski quantum Markov chains* ([2501.09708](https://arxiv.org/abs/2501.09708));
- *A structural theory of quantum metastability* ([2510.08538](https://arxiv.org/abs/2510.08538)).

**Not retrieved** (§0): the Hammersley–Clifford manuscript (its content is known through [Cl90] §2); Besag's displayed equations in §3; Grimmett's full text; Spitzer; Averintsev; Sherman; Lauritzen 1996; Geiger et al. 2006; Petz 1986/1988; Accardi–Frigerio 1983; Poulin–Hastings 2011; Chen–Kastoryano–Gilyén 2023.

**Statements marked *from memory*** were not re-read here and should be checked before use (rule P3.3): the modular-group sign convention; Popa's commuting squares; the étale topology on homoclinic relations; Kemeny–Snell lumpability; Efron–Stein as Möbius inversion for a product measure; variable-elimination fill-in.
