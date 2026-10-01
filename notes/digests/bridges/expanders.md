# Expanders, and the bridges to Gibbs measures, Markov semigroups and NCG

### Digest for the bridges stream: what expander theory actually proves, and which of its connections to Markov random fields, Glauber dynamics, quantum Gibbs samplers and noncommutative geometry are theorems

*Bridges digest, 2026-10-01. Written for the [research programme](../../research-program.md). Read first: [research-program.md](../../research-program.md), [local-to-global-unlocks.md](../../local-to-global-unlocks.md) (prong 3, whose §1 already uses the expander gap as "Reduction 1"), [conditional-arrow-algebra.md](../../conditional-arrow-algebra.md) (Note 1: KMS/Gibbs states on histories, the coboundary criterion, Jenčová–Petz sufficiency), [mlp-bridge.md](../../mlp-bridge.md) (prong 2, dictionary v1). Standing constraint respected: nothing here reads objects as activation vectors; every transfer to the framework is stated for abstract faces, arrows, algebras and states.*

*Sibling digests in this folder (written in parallel, not read for this file):* [hdx-spectral-independence.md](hdx-spectral-independence.md) covers §9.5 in depth (Gibbs distributions as weighted simplicial complexes whose links are pinnings); [arxiv-2609.38007.md](arxiv-2609.38007.md) and [chat-2609.38007-retrieval-status.md](chat-2609.38007-retrieval-status.md) cover Yang's paper and the Chen–Rouzé Lindbladian of §9.9. Where they and this file overlap, disagreements should be resolved against the sources, not by majority (rule C4).

**Labels.** **Fact** = a statement taken from a source that was opened in this session, with its numbering. **Fact (memory)** = a standard statement that was *not* re-read here; to be checked before use (programme rule P3.3). **Derivation** = a short argument made here. **Interpretation** = my reading. **Conjecture** = what would have to be true in the framework. **Guard** = what does not transfer.

---

## 0. Sources, retrieval, and the verdict in brief

**Primary sources read in full or in the relevant sections.**

| source | route | what was read |
|---|---|---|
| T. Tao, *254B Notes 1: Basic theory of expander graphs* (2 Dec 2011), [blog](https://terrytao.wordpress.com/2011/12/02/245b-notes-1-basic-theory-of-expander-graphs/) | Exa fetch, full text | §§1–4: Defs 1–4, Lemma 3, Exs 6–13, Prop 14 with proof, Rem 17, Exs 18–24, 27, Props 28–29, Cor 30, Rems 36–38 |
| S. Hoory, N. Linial, A. Wigderson, *Expander graphs and their applications*, Bull. AMS 43(4) (2006) 439–561 | Exa fetch of `cs.huji.ac.il/~nati/PAPERS/expander_survey.pdf`, full text (355 k characters) | §§2, 3, 4.4–4.6, 5, 6–7 (statements), 9, 10 (Def 10.3), 11.1 (Kazhdan constant), 12, 13.3 |
| E. Mossel, A. Sly, *Exact thresholds for Ising–Gibbs samplers on general graphs*, Ann. Probab. 41 (2013), [0903.2906](https://arxiv.org/abs/0903.2906) | alphaXiv | Thms 1–3, Lemmas 1–3, 13, §§1.4, 2.2–2.3, 5.4 |
| M. Dyer, L. A. Goldberg, M. Jerrum, *Matrix norms and rapid mixing for spin systems*, AAP 19 (2009), [math/0702744](https://arxiv.org/abs/math/0702744) | alphaXiv | §1 (Lemmas 1–3), §2.2, §3.1 (Lemmas 16–17, Cors 18–19, 25), §3.2 |
| R. Eldan, F. Koehler, O. Zeitouni, *A spectral condition for spectral gap*, [2007.08200](https://arxiv.org/abs/2007.08200) | alphaXiv | Thm 1 with proof, Def 5, Thm 6 (Wu), Lemmas 7–9, Thms 10–12, §5 |
| Z. Chen, K. Liu, E. Vigoda, *Optimal mixing of Glauber dynamics: entropy factorization via high-dimensional expansion*, [2011.02075](https://arxiv.org/abs/2011.02075) (STOC 2021) | alphaXiv | Defs 1.9–1.11, Thm 1.12, Rem 1.13, Claim 1.18, Thm 1.19, Def 2.1, Thm 2.9, Thms 6.3–6.4 |
| N. Berger, C. Kenyon, E. Mossel, Y. Peres, *Glauber dynamics on trees and hyperbolic graphs*, [math/0308284](https://arxiv.org/abs/math/0308284) | alphaXiv | Props 1.1–1.3, Thms 1.4–1.5, Lemma 2.1, §§3–4, Problem 4 |
| M. J. Kastoryano, F. G. S. L. Brandão, *Quantum Gibbs samplers: the commuting case*, CMP, [1409.3435](https://arxiv.org/abs/1409.3435) | alphaXiv | Thms 1, 3 (informal), Defs 14–15, 21, Prop 20 with proof, Thms 23, 28, 30–31 |
| C.-F. Chen, C. Rouzé, *Quantum Gibbs states are locally Markovian*, [2504.02208](https://arxiv.org/abs/2504.02208) | alphaXiv | abstract, Fig. 2 caption, eq. (1.1), Thm III.1, §IV.C, Lemmas VII.2–VII.3, Cors III.1–III.2 |
| T. H. Yang, *Improved estimate of local Markovianity for quantum Gibbs states*, [2609.38007](https://arxiv.org/abs/2609.38007) (29 Sep 2026) | alphaXiv | abstract, §I, Table I, eq. (1.2), §VII |
| R. Willett, G. Yu, *Higher index theory for certain expanders and Gromov monster groups I/II*, [1012.4150](https://arxiv.org/abs/1012.4150), [1012.4151](https://arxiv.org/abs/1012.4151) | Exa (arXiv HTML and author PDF) | abstracts, §1 (Cor 1.7), Ex 5.3 (basic Kazhdan projection), Cor 6.2 |
| M. B. Hastings, *Random unitaries give quantum expanders*, PRA 76 (2007), [0706.0556](https://arxiv.org/abs/0706.0556) | Exa (ar5iv) | definition, eq. (3), main result, eq. (12) |
| A. Marcus, D. Spielman, N. Srivastava, *Interlacing families I*, Ann. Math. 182 (2015), [1304.4132](https://arxiv.org/abs/1304.4132) | Exa | abstract, §2, Thms 5.3, 5.5, 5.6 |

**Could not be retrieved.**
- The user's local file `C:\Users\User\Downloads\expander_survey.pdf` is not reachable from this container. Hoory–Linial–Wigderson is used as the stand-in. *Observation, unverified:* the HLW survey is served from Linial's homepage under exactly the file name `expander_survey.pdf`, so the local file is plausibly the same document; if it is a different survey, the section numbers below refer to HLW, not to it.
- `arxiv.org` and `pages.uoregon.edu` refuse direct connections through the proxy (403 on CONNECT); arXiv papers were read through alphaXiv or Exa instead. The Levin–Peres–Wilmer book was therefore not opened; its results are quoted only as other papers state them.
- Not re-read, quoted from memory and marked so: the Witsenhausen/Hirschfeld–Gebelein–Rényi maximal-correlation identity, the Kesten–McKay density, the Jerrum–Sinclair conductance inequality, Popa's commuting squares, the precise form of Hammersley–Clifford, the Stroock–Zegarlinski / Martinelli–Olivieri equivalences (quoted only as Chen–Rouzé and Mossel–Sly report them).

**Verdict on the user's intuition, in brief (Interpretation; details in §§8–11).**
1. *Gibbs states / Markov random fields share an essence with what drives expanders.* **Partly yes, in a precise form.** The Markov property of a Gibbs measure is a *static, exact, hereditary* statement (conditionals of a Gibbs measure on a subregion are Gibbs measures there, and depend only on the separator); it holds at every temperature, including at phase transitions (Yang, §I). It is not expansion. What expander theory supplies is the *other* half: a uniform decorrelation constant for a local operator. Every theorem that connects the two (Dobrushin, path coupling, Hayes, Mossel–Sly, spectral independence, Kastoryano–Brandão) has the shape *hereditary conditional structure + uniform decorrelation of conditionals ⇒ global spectral gap*, which is exactly the shape of high-dimensional-expander local-to-global theorems. The algebra that proves the zig-zag theorem (split a function into its conditional mean and its fluctuation; bound the cross term by an angle) is the same algebra that proves Kastoryano–Brandão's Prop 20.
2. *Expansion of the interaction graph helps Gibbs mixing.* **Only at high temperature, and only through spectral norms.** The Ramanujan bound of a random regular interaction graph sets the regime $\beta<1/(4\sqrt{d-1})$ for fast mixing of a diluted spin glass (Eldan–Koehler–Zeitouni §5), better by $\sqrt d$ than Dobrushin. At low temperature expansion is the *enemy*: random regular graphs give $\exp(\Omega(n))$ mixing beyond the tree threshold and are the hardness gadgets. Mossel–Sly say outright that block-dynamics and strong-spatial-mixing methods fail on expanders because boundary is proportional to volume.
3. *The "time-averaged detailed-balanced Lindbladian based on single-Pauli jumps on A".* The phrase is the caption of Fig. 2 of **Chen–Rouzé, arXiv:2504.02208**, not of 2609.38007; Yang's 2609.38007 is explicitly a *static* proof that does not use that Lindbladian. Its mechanism is deliberately **gap-free** (the Dirichlet form of the time average is $\le 2/t$ unconditionally). The expander-type quantity is what Chen–Rouzé *avoid*: a local spectral gap of the $A$-restricted Lindbladian, which they say would remove their $e^{\mu|A|}$ prefactor. So the bridge is real but sits one step to the side of where it was sensed: the Lindbladian is the noncommutative heat-bath resampling of $A$, i.e. a dynamically produced conditional expectation, and the expander question it raises is the *trickling-down* question for the gaps of these local generators (§9.9, §11).
4. *NCG ties these together.* There are real theorems here, not only analogies: property (T)/(τ) produces expanders and the Kazhdan constant is equivalent to the gap (HLW §11); expanders produce non-compact ghost projections in Roe algebras and hence counterexamples to the coarse Baum–Connes conjecture (Higson; Higson–Lafforgue–Skandalis; Willett–Yu); detailed balance for quantum Gibbs samplers is defined through the KMS inner product, i.e. through modular theory; quantum expanders are channels with a gap and obey a quantum Alon–Boppana bound. The *single* essence across them is my interpretation (§8, §10): **a spectral gap is what makes a global projection a norm-limit of local operators, while leaving it invisible to every local observer.**

---

## 1. Definitions: three equivalent views and one group-theoretic one

Conventions. HLW: $G=(V,E)$ undirected, $d$-regular, $n=|V|$, loops and multi-edges allowed; $E(S,T)$ is the set of *directed* edges from $S$ to $T$; adjacency eigenvalues $d=\lambda_1\ge\lambda_2\ge\dots\ge\lambda_n$; $\lambda(G)=\max(|\lambda_2|,|\lambda_n|)$; $\hat A=A/d$. An **$(n,d)$-graph** is $d$-regular on $n$ vertices; an **$(n,d,\alpha)$-graph** additionally has $\lambda(G)\le\alpha d$. Tao: undirected, loop-free, multiplicity-free in Notes 1 (loops allowed later for Cayley graphs).

### 1.1 Combinatorial

**Fact** (HLW Def 2.1–2.2). Edge boundary $\partial S=E(S,\bar S)$; **expansion ratio** $h(G)=\min_{|S|\le n/2}|\partial S|/|S|$. A sequence of $d$-regular graphs of increasing size is a **family of expanders** if $h(G_i)\ge\varepsilon>0$ for all $i$. Vertex expansion and the expansion of sets of a given size are the two other variants (HLW §4.6, §10). Explicitness (Def 2.3): *mildly explicit* (the $j$-th graph in time poly$(j)$), *very explicit* (the $k$-th neighbour of $v$ in time poly$(\log n)$).

**Fact** (HLW Def 10.3). A bipartite graph with left degree $D$ is a **$(K_{\max},\epsilon)$-lossless expander** if every set of $K\le K_{\max}$ left vertices has at least $(1-\epsilon)DK$ neighbours; equivalently most neighbours of a small set are *unique* neighbours.

### 1.2 Spectral

**Fact** (Tao Def 4, Rem 5). For $k$-regular $G$, $\lambda_1=k\ge\lambda_n\ge-k$ (Lemma 3). $G$ is a **one-sided $\varepsilon$-expander** if $\lambda_2\le(1-\varepsilon)k$, **two-sided** if also $\lambda_n\ge-(1-\varepsilon)k$; equivalently the normalised Laplacian $\Delta=1-\frac1kA$ has a spectral gap $\varepsilon$. Qualitatively (Ex 6): $\lambda_2=k$ iff disconnected; $\lambda_n=-k$ iff a bipartite component exists. So a single graph is "an expander for some $\varepsilon$" iff connected; the notion only has content quantitatively or for families.

**Fact** (Tao Ex 13, Poincaré form). $G$ is a one-sided $\varepsilon$-expander iff $\|\nabla f\|_{\ell^2}^2\ge2k\varepsilon\|f\|_{\ell^2}^2$ for all mean-zero $f$, where $|\nabla f(v)|^2=\sum_{w\sim v}|f(w)-f(v)|^2$. This is the bridge sentence to Markov semigroups: *a spectral gap is a Poincaré inequality for the Dirichlet form of the walk.*

**Fact** (Tao Ex 7). $\sum\lambda_i=0$, $\sum\lambda_i^2=nk$, hence $\max(|\lambda_2|,|\lambda_n|)\ge\sqrt k-o(1)$ (HLW Claim 2.8: $\lambda^2\ge d\frac{n-d}{n-1}$).

**Fact** (Tao Exs 9–12). The cycle $\mathbb Z/n$ has eigenvalues $2\cos(2\pi j/n)$ and is not an expander family; $K_n$ has $\lambda_2=\dots=\lambda_n=-1$; $K_{n/2,n/2}$ has $\lambda_n=-n/2$, $\lambda_2=\dots=\lambda_{n-1}=0$ (one-sided $1$-expander, not two-sided); for the complement graph $\lambda_i(G^c)=-1-\lambda_{n+2-i}(G)$.

### 1.3 Random-walk mixing

**Fact** (HLW Thm 3.2, Thm 3.3, Lemma 3.4). For an $(n,d,\alpha)$-graph and any probability vector $p$,
$$
\|\hat A p-u\|_2\le\alpha\|p-u\|_2,\qquad \|\hat A^tp-u\|_2\le\alpha^t,\qquad \|\hat A^tp-u\|_1\le\sqrt n\,\alpha^t,
$$
$u$ the uniform vector. Proof architecture: $u$ is invariant and $p-u\perp u$ lies in the span of the non-trivial eigenvectors.

**Fact** (Tao Ex 24, Ex 27). For $k$-regular families and any $\alpha>1/2$: two-sided expander family $\iff$ $\|\mu^{(i)}-1/|V_n|\|_{\ell^2}\le|V_n|^{-\alpha}$ for all $i\ge C\log|V_n|$ and all starting vertices. One-sided $\iff$ the same for the *lazy* walk. Rem 25: the uniform law is a "very strong attractor" for all initial laws.

**Fact** (HLW §3.1.2). With $H_2(p)=-2\log\|p\|_2$ and $p=u+f$, $\mu=\|f\|/\|p\|$: $H_2(\hat Ap)\ge H_2(p)-\log(1-(1-\alpha^2)\mu^2)$; Rényi-2 entropy strictly increases unless $p$ is uniform, faster for smaller $\alpha$. For Shannon entropy the per-step increase is governed by the **log-Sobolev constant**, which HLW place "beyond the scope of this survey". This is the point where expander theory hands over to the theory of Markov semigroups (§9.5: modified log-Sobolev).

### 1.4 Group-theoretic: the Kazhdan constant

**Fact** (HLW Def 11.18 and the inequality after it). For a group $H$ and $S\subset H$,
$$
K(H,S)=\min_{v\perp\mathbf 1}\max_{h\in S}\frac{\|r(h)v-v\|^2}{\|v\|^2}
=\min_{\rho\ne1\ \text{irred}}\ \min_{v}\max_{h\in S}\frac{\|\rho(h)v-v\|^2}{\|v\|^2},
\qquad \frac{K(H,S)}{2|S|}<g(H,S)<\frac{K(H,S)}{2},
$$
$g(H,S)$ the spectral gap of the normalised Cayley graph. Claim 11.19: if each element of $\tilde S$ is a product of at most $m$ elements of $S$, then $K(H,S)\ge K(H,\tilde S)/m$ (not true for the gap). Margulis's first explicit expanders came from Kazhdan's property (T) of $SL_3(\mathbb Z)$; Lubotzky's from Selberg's $3/16$ theorem ("property $(\tau)$") for $SL_2(\mathbb Z)$ quotients (HLW §11 intro and §2.2). This is the first NCG bridge (§9.11): the gap is the *uniform isolation of the trivial representation*.

---

## 2. Cheeger's inequality and its converse

**Fact** (HLW Thm 2.4 = Thm 4.11; Dodziuk 1984, Alon–Milman 1985, Alon 1986). For finite connected $d$-regular $G$,
$$
\frac{d-\lambda_2}{2}\ \le\ h(G)\ \le\ \sqrt{2d(d-\lambda_2)} .
$$
Tao's normalisation (Rem 17, eq. (5)): if $\lambda_2=(1-\varepsilon)k$ then $\frac{\varepsilon}{2}k\le h(G)\le\sqrt{2\varepsilon}\,k$.

**Fact** (tightness, HLW §4.5). Lower bound tight for the hypercube $Q_d$: $h=1$ (on a subcube) and gap $d-\lambda=2$. Upper bound tight for the cycle: $h(C_n)=\Theta(1/n)$, gap $\Theta(1/n^2)$.

**Fact** (continuous original, HLW Thm 4.9, Cheeger 1970). For a compact Riemannian manifold with Cheeger constant $h(M)=\inf_A\mu_{n-1}(\partial A)/\min(\mu_n(A),\mu_n(M\setminus A))$, the first positive Laplace eigenvalue satisfies $\lambda\ge h^2/4$ (Buser gives the converse direction).

**Proof architecture.**
- *Gap ⇒ expansion* (easy; Tao Prop 14 (i)⇒(ii), HLW §4.5.1). Test the Rayleigh quotient on $1_F-\frac{|F|}{|V|}1$ (Tao) or $f=|\bar S|1_S-|S|1_{\bar S}$ (HLW). Tao obtains $|\partial F|\ge\varepsilon k|F|/4$.
- *Expansion ⇒ gap* (hard; Tao Prop 14 (ii)⇒(i)). Two reductions and one integral. (1) Mean-zero $f$ is split $f=f_+-f_-$; $\langle Af_+,f_-\rangle\ge0$ lets one treat each sign separately; one of them lives on at most half the vertices, and a Markov-inequality truncation handles the other. (2) For non-negative $f$ supported on $\le|V|/2$ vertices use the **"wedding cake" (co-area) decomposition** $f=\int_0^\infty1_{F_t}dt$, $F_t=\{|f|>t\}$, $\|f\|^2=\int2t|F_t|dt$, and bound $\langle A1_{F_s},1_{F_t}\rangle$ by $k|F_t|$ or by $(k-c)|F_s|$ (the latter from $|\partial F_s|\ge c|F_s|$), splitting at $s=(1-\varepsilon)t$. The co-area step is what turns a statement about *sets* into a statement about *all functions*; it is the same step that proves the Riemannian inequality.
- *Reversible Markov chains.* HLW note that the most useful generalisation is to reversible chains, where edge expansion becomes **conductance** (Jerrum–Sinclair 1989) and has "a huge impact on the analysis of convergence of Monte Carlo algorithms". **Fact (memory):** for a reversible chain with conductance $\Phi$, $\Phi^2/2\le1-\lambda_2\le2\Phi$. This is the form used to prove *slow* mixing of Gibbs samplers by exhibiting a bottleneck in configuration space (§9.7).

**Fact** (converse of the expander mixing lemma, HLW Lemma 2.6, Bilu–Linial). If $\big||E(S,T)|-d|S||T|/n\big|\le\rho\sqrt{|S||T|}$ for all disjoint $S,T$, then $\lambda\le O(\rho(1+\log(d/\rho)))$, and this is tight. So discrepancy and spectral gap are equivalent up to a logarithm, a tighter link than Cheeger's.

**Fact** (expansion of small sets).
- Tanner (HLW Thm 4.15): an $(n,d,\alpha)$-graph has vertex expansion $\Psi_V(G,\rho n)\ge1/(\rho(1-\alpha^2)+\alpha^2)$. Proof: compare $\|\hat A1_S\|^2\le\rho n(\rho+\alpha^2(1-\rho))$ with the Cauchy–Schwarz lower bound $\rho n|S|/|\Gamma(S)|$.
- Kahale: near-Ramanujan graphs give vertex expansion $\approx d/2$ for small linear sets, and $d/2$ is essentially the limit of spectral methods (HLW §4.6.1).
- Random graphs (HLW Thm 4.16): for fixed $d\ge3$ and every $\delta>0$ there is $\varepsilon>0$ with, for almost every $(n,d)$-graph, $\Psi_E(G,\varepsilon n)\ge d-2-\delta$, and for almost every $d$-regular bipartite graph $\Psi_V(G,\varepsilon n)\ge d-1-\delta$; these are best possible. So **random graphs beat the eigenvalue bound** on small sets; lossless expanders (HLW §10) are the explicit version.

---

## 3. The expander mixing lemma and its walk version

**Fact** (HLW Lemma 2.5; Tao Ex 18). For a $d$-regular graph with $\lambda=\lambda(G)$ and all $S,T\subseteq V$,
$$
\Big|\,|E(S,T)|-\frac{d|S||T|}{n}\Big|\ \le\ \lambda\sqrt{|S||T|},
$$
with Tao's refinement $\sqrt{|S|(1-|S|/n)\,|T|(1-|T|/n)}$ on the right. Proof: expand $1_S,1_T$ in an orthonormal eigenbasis; the $v_1=\mathbf1/\sqrt n$ component gives $d|S||T|/n$; the rest is $\sum_{i\ge2}\lambda_i\alpha_i\beta_i$, bounded by Cauchy–Schwarz. **Tao Ex 19** gives the one-sided variant: one-sided expander family $\iff$ $|E(F_1,F_2)|\le(k-c)\sqrt{|F_1||F_2|}$ for sets of size $\le|V|/2$.

**Consequences** (HLW §2.4; Tao Exs 20–23). In an $(n,d,\alpha)$-graph: every independent set has size $\le\alpha n$; $\chi(G)\ge1/\alpha$; diameter $O(\log n)$ (balls grow by a factor $1+\epsilon$ with $\epsilon=\tfrac12-\alpha$ while below $n/2$). Removing $m$ edges leaves a component of size $\ge n-Cm$ (Tao Ex 21). Lipschitz functions concentrate: $|\{|f-M|\ge\lambda K\}|\le Cne^{-c\lambda}$ (Tao Ex 23).

**Reading** (HLW §3.2). The lemma compares two experiments: pick $(i,j)$ uniformly from $V^2$, or pick $i$ uniformly and $j$ a random neighbour. Success probabilities for "$i\in S,j\in T$" differ by at most $\alpha$: $\big|\frac{|S||T|}{n^2}-\frac{|E(S,T)|}{dn}\big|\le\alpha$. *The edge measure is $\alpha$-close to the product measure.*

**Fact** (walks resemble independent samples; HLW Thm 3.6, Ajtai–Komlós–Szemerédi 1987, Alon–Feige–Wigderson–Zuckerman 1995). For $B\subset V$ of density $\beta$ and a walk $X_0,\dots,X_t$ from a uniform start,
$$
\Pr[\forall i:\ X_i\in B]\le(\beta+\alpha)^t .
$$
Architecture: $\Pr=\|(P\hat A)^tPu\|_1$ with $P$ the coordinate projection on $B$ (Lemma 3.7), and the one-line key lemma
$$
\|P\hat AP\|_{2\to2}\le\beta+\alpha\qquad\text{(Lemma 3.8)},
$$
proved by writing $Pv=u+z$, $z\perp u$: the constant part costs $\beta$ (Cauchy–Schwarz on a set of density $\beta$), the fluctuation costs $\alpha$. Refinements: if $\beta>6\alpha$, $\beta(\beta+2\alpha)^t\ge\Pr\ge\beta(\beta-2\alpha)^t$ (Thm 3.9); for any $K\subset\{0..t\}$, $\Pr[X_i\in B\ \forall i\in K]\le(\beta+\alpha)^{|K|-1}$ (Thm 3.10); for time-varying sets of densities $\beta_i$, $\Pr\le\prod_{i<t}(\sqrt{\beta_i\beta_{i+1}}+\alpha)$ (Thm 3.11). Gillman (1998) gives the Chernoff bound for expander walks (cited in HLW §3.2).

**Interpretation.** Lemma 3.8 is the cleanest local-to-global statement in the subject: *conditioning on a local event (projecting onto $B$) does not destroy decorrelation; it only adds the density of the event.* The same structure ($\|PTP\|\le$ "mass" + "gap") reappears in Kastoryano–Brandão's Prop 20 (§9.8) with conditional expectations in place of coordinate projections.

---

## 4. Alon–Boppana, the tree, and Ramanujan graphs

**Fact** (spectrum of the tree; HLW §5.1, partial proof after Friedman). The $d$-regular tree $T_d$ has $\mathrm{spec}(A_{T_d})=[-2\sqrt{d-1},2\sqrt{d-1}]$. Architecture: $\lambda\in\mathrm{spec}$ iff $\delta_v\notin\mathrm{Range}(\lambda-A)$; reduce to spherical $f$, $f(u)=x_i$ at distance $i$, with recurrence $\lambda x_0=dx_1+1$, $\lambda x_i=x_{i-1}+(d-1)x_{i+1}$; roots $\rho_{1,2}$ of $\lambda\rho=1+(d-1)\rho^2$. For $|\lambda|<2\sqrt{d-1}$ the roots have modulus $1/\sqrt{d-1}$, so $|x_i|=\Theta((d-1)^{-i/2})$ against $\Theta((d-1)^i)$ vertices at distance $i$: $\|f\|_2=\infty$, so $\lambda$ is in the spectrum. For $|\lambda|>2\sqrt{d-1}$ an $\ell^2$ solution exists. **The edge of the spectrum is the point where the squared decay of an eigenfunction exactly balances the branching of spheres.** (This $\ell^2$ balance reappears as the Kesten–Stigum threshold in §9.6.)

**Fact** (Alon–Boppana; HLW Thm 2.7, Thm 5.3 [Nilli 1991, Friedman 1993], Cor 5.4). There is a constant $c$ ($\approx2\pi^2$ from the second proof) such that every $(n,d)$-graph of diameter $\Delta$ has
$$
\lambda_2(G)\ \ge\ 2\sqrt{d-1}\,\big(1-c/\Delta^2\big),\qquad\text{hence}\qquad \lambda_2\ge2\sqrt{d-1}\,(1-O(1/\log^2n)).
$$
Two proofs. (I) Moment method: $\mathrm{tr}A^{2k}\ge n\,t_{2k}$, $t_{2k}$ = closed walks of length $2k$ at a vertex of $T_d$ (every regular graph has at least as many as its universal cover); applied to $f=\delta_s-\delta_t$ at distance $\Delta$. (II) Rayleigh quotient of a lifted eigenfunction of a truncated tree (Claim 5.5, Lemma 5.6, Cor 5.7: $\lambda_2(A)\ge\lambda_2(T_k)$).

**Fact** (Serre; HLW Thm 5.8). For every $d$ and $\varepsilon>0$ there is $c(\varepsilon,d)>0$ such that every $(n,d)$-graph has at least $cn$ eigenvalues $>2\sqrt{d-1}-\varepsilon$. (Cioabă's proof via $\mathrm{tr}(A+dI)^k$.)

**Fact** (universal-cover form; HLW Thm 6.6, Greenberg–Lubotzky). If $\{G_i\}$ share a universal cover $T$, then $\lambda(G_i)\ge\rho(T)-o(1)$. HLW Def 6.7: $G$ is **Ramanujan** if $\lambda(G)\le\rho(\hat G)$ (spectral radius of the universal cover). In the regular case, Def 5.11: $\lambda(G)\le2\sqrt{d-1}$.

**Fact** (constructions).
- Lubotzky–Phillips–Sarnak 1988, Margulis 1988, Morgenstern 1994 (HLW Thm 5.12): for every prime $p$ and $k\ge1$ there are infinitely many $d$-regular Ramanujan graphs with $d=p^k+1$. The LPS graph $X^{p,q}$ ($p,q\equiv1\bmod4$ distinct primes) is a Cayley graph of $PGL(2,q)$ with the $p+1$ generators coming from Jacobi's four-square representations of $p$; it is $(p+1)$-regular of size $\Theta(q^3)$; the optimal bound uses Deligne's proof of the Ramanujan conjecture over finite fields (HLW §11 intro).
- Marcus–Spielman–Srivastava (Ann. Math. 2015, Thms 5.3, 5.5, 5.6): **every graph $G$ has a signing $s$ with all eigenvalues of $A_s$ at most $\rho(T)$**, $T$ its universal cover; hence (via 2-lifts, Bilu–Linial Lemma 3.1) infinite families of $d$-regular *bipartite* Ramanujan graphs for every $d\ge3$, and $(c,d)$-biregular ones with non-trivial eigenvalues $\le\sqrt{c-1}+\sqrt{d-1}$. Architecture: the expected characteristic polynomial of a random signing is the **matching polynomial** $\mu_G$ (Godsil–Gutman); $\mu_G$ is real-rooted (Heilmann–Lieb) and its roots are bounded by $\rho(T)$ because $\mu_G$ divides the matching polynomial of Godsil's **path tree** (the tree of self-avoiding walks); the signings form an **interlacing family**, so some member has largest root at most that of the expectation. Non-bipartite Ramanujan graphs of all degrees: open as of HLW (Conj 5.13); MSS settle the bipartite case. Later work on the edge statistics of random regular graphs is reported to imply the non-bipartite case for all degrees; not retrieved here, check before citing.
- Bilu–Linial (HLW Thm 6.12): mildly explicit $(n,d,\alpha)$-graphs with $\alpha d=O(\sqrt{d\log^3d})$, by 2-lifts (Conjectures 6.8–6.10: every ($d$-regular Ramanujan) graph has a 2-lift with new eigenvalues in $[-2\sqrt{d-1},2\sqrt{d-1}]$).

**Interpretation (one shared object).** Godsil's path tree in MSS is the same tree of self-avoiding walks that Weitz uses to compute Gibbs marginals on a graph from a tree (§9.4, §9.5; Chen–Liu–Vigoda Thm 6.3 use it for the monomer–dimer model "extending known results on the univariate matching polynomial"). *The universal cover / SAW tree controls both the extremal spectrum of a graph and the extremal decay of Gibbs correlations on it.* This is a genuine common object, not an analogy: the matching polynomial *is* the partition function of the monomer–dimer Gibbs measure.

---

## 5. Random regular graphs and random matrices

**Fact** (existence by counting; Tao §4). For even $k=2l$ large, take $l$ uniform permutations $\pi_1..\pi_l$ of $[n]$ and join $v\sim\pi_i(v)$. Prop 28: the result is $2l$-regular with probability $\ge c_k-o(1)$ (Bonferroni inequalities; the failure count is asymptotically Poisson, Ex 35; McKay's "swapping" alternative, Rem 36). Prop 29: there is $\varepsilon(k)>0$ with $\Pr[2l\text{-regular and not a one-sided }\varepsilon\text{-expander}]=o(1)$, by a union bound over $F\subset F'$ with $\Pr[\pi_i(F)\subset F']\le((r+r')/n)^r$. Rem 38: the permutation model is *contiguous* to the uniform model, so $1-o(1)$ of all $k$-regular graphs ($k\ge3$) are $\varepsilon_k$-expanders. Ex 37: two-sided as well.

**Fact** (bulk; HLW §7.1).
- Wigner (Thm 7.1): for a real symmetric matrix with i.i.d. entries of variance $\sigma^2$ and finite moments, the empirical eigenvalue law of $A_n/\sqrt n$ converges to the semicircle on $[-2\sigma,2\sigma]$.
- McKay (Thm 7.2): for $d$-regular graphs with $C_k(G_n)=o(|V(G_n)|)$ cycles of each fixed length $k\ge3$ (true a.s. for random regular graphs), the empirical spectral law converges to the Kesten–McKay law. **Fact (memory):** its density is $\frac{d\sqrt{4(d-1)-x^2}}{2\pi(d^2-x^2)}$ on $|x|\le2\sqrt{d-1}$; it tends to a semicircle after rescaling as $d\to\infty$ (stated in HLW).

**Fact** (edges; HLW §7.2–7.3).
- Füredi–Komlós, Vu (Thm 7.4): zero-mean entries bounded by $K$, variance $\sigma^2$: with probability $1-o(1)$, all $|\lambda_i|<2\sigma\sqrt n+O(n^{1/3}\log n)$.
- Broder–Shamir (Thm 7.5): almost every $2d$-regular graph (permutation model) has $\lambda=O(d^{3/4})$, by the trace method.
- Friedman lifts (Thm 7.9): for almost all high lifts of $G$, new eigenvalues are $\le\sqrt{\lambda_1\rho}+o(1)$.
- **Friedman** (Thm 7.10): for every $\varepsilon>0$, $\Pr[\lambda(G)\le2\sqrt{d-1}+\varepsilon]=1-o(1)$ for a random $(n,d)$-graph. A new proof and extension to random lifts: Bordenave, arXiv 1502.04482 (cited by EKZ). Explicit near-Ramanujan graphs of every degree: Mohanty–O'Donnell–Paredes, STOC 2020 (cited by EKZ).
- Experiments (HLW table, 1000 random 4-regular graphs, $n$ up to 400 000): $\Pr[\lambda<2\sqrt3]\approx0.62$–$0.68$, median and mean approach $2\sqrt3$ from below, standard deviation $\to0$; Novikoff conjectures Tracy–Widom fluctuations and a limiting probability strictly between 0 and 1 of being Ramanujan. *(Not retrieved here: later work claiming to prove this; to be checked before citing.)*

**Proof architecture: the trace (moment) method.** Estimate $\mathrm{tr}A^{2k}$ combinatorially as a count of closed walks; subtract $\lambda_1^{2k}$; for $k$ large the remainder is dominated by $\lambda_2^{2k}$. Large $k$ helps the dominance but makes the combinatorics harder; Friedman's proof (very long) and Bordenave's (non-backtracking operator) are refinements. The same method gives Alon–Boppana (I) and Hastings' quantum version (§9.10).

---

## 6. Zig-zag and replacement products

**Fact** (definition; HLW §9.3, Def 9.3). $G$ an $(n,m,\alpha)$-graph with edges at each vertex numbered $e_v^1..e_v^m$; $H$ an $(m,d,\beta)$-graph on $[m]$. Replace each $v$ by a **cloud** $\{v\}\times[m]$. The **replacement product** $G\,ⓡ\,H$ has the edges of $G$ between clouds plus a copy of $H$ in each cloud. The **zig-zag product** $G\,ⓩ\,H$ on $V(G)\times[m]$ joins $(v,i)$ to $(u,j)$ iff there are $k,l$ with $(i,k),(l,j)\in E(H)$ and $e_v^k=e_u^l$: a step is *random step in the cloud, deterministic step between clouds, random step in the new cloud*. It is an $(nm,d^2)$-graph.

**Fact** (Zig-Zag Theorem; HLW Thm 9.1, Reingold–Vadhan–Wigderson 2002). $G\,ⓩ\,H$ is an $(nm,d^2,\varphi(\alpha,\beta))$-graph with
$$
\alpha,\beta<1\Rightarrow\varphi<1,\qquad \varphi\le\alpha+\beta,\qquad \varphi\le1-(1-\beta^2)\frac{1-\alpha}{2},
$$
and $\varphi\ge\max\{\alpha,\beta\}$ (footnote 16: zig-zag cannot improve expansion). Martin–Randall proved the third bound for the replacement product "in a different context (analysis of expansion in Markov chains via decomposition, which may be viewed as reversing the replacement product)"; the second bound need not hold for it.

**Proof architecture** (HLW, weak form $\varphi\le\alpha+\beta+\beta^2$). Transition matrix $Z=\tilde BP\tilde B$, $\tilde B=\hat B\otimes I_n$ (cloud steps), $P$ the involutive permutation matrix of the inter-cloud step. Split $f\perp\mathbf1$ as
$$
f=f^{\parallel}+f^{\perp},\qquad f^\parallel=\text{cloud average},\quad f^\perp\ \text{sums to 0 on every cloud}.
$$
Then $\tilde Bf^\parallel=f^\parallel$, $\|\tilde Bf^\perp\|\le\beta\|f^\perp\|$, and $\langle f^\parallel,Pf^\parallel\rangle=\langle g,\hat A_Gg\rangle\le\alpha\|f^\parallel\|^2$ with $g(v)=\sqrt m f^\parallel(v,i)$. Hence $|\langle f,Zf\rangle|\le\alpha\|f^\parallel\|^2+2\beta\|f^\parallel\|\|f^\perp\|+\beta^2\|f^\perp\|^2$, whose maximum over $\|f^\parallel\|^2+\|f^\perp\|^2=1$ is the top eigenvalue of $\begin{pmatrix}\alpha&\beta\\\beta&\beta^2\end{pmatrix}$, at most $\alpha+\beta+\beta^2$.

**Entropy reading** (HLW §9.4). If the conditional law on clouds is far from uniform, the $H$-step raises entropy; if it is near uniform, the middle step is a permutation (keeps total entropy) whose $G$-marginal is a genuine $G$-step (raises the $G$-marginal's entropy), so the $H$-marginal's entropy must drop, and the third step is in the good case. "The key is that Step 2 is simultaneously a permutation ... and an operation whose $G$-marginal is simply a random step on $G$."

**Fact** (explicit family; HLW §9.2, Prop 9.2). With $H$ a $(d^4,d,1/4)$-graph (found by exhaustive search), $G_1=H^2$, $G_{n+1}=(G_n)^2\,ⓩ\,H$: every $G_n$ is a $(d^{4n},d^2,1/2)$-graph. Tensoring makes it very explicit.

**Fact** (Cayley graphs; HLW §11.2, Alon–Lubotzky–Wigderson 2001). Under conditions on the groups and generators, $C(A,S_A)\,ⓩ\,C(B,S_B)$ is a Cayley graph of the semidirect product $A\rtimes B$; iterated semidirect products give Cayley expanders (Cor 11.28, Rozenman–Shalev–Wigderson). HLW §11.4: expansion is not a group property (it depends on generators).

**Fact** (SL = L; HLW §9.5, Reingold 2005). Make the input $D$-regular with $D=d^{16}$; with $H$ a $(d^{16},d,1/2)$-graph set $G_1=G$, $G_{i+1}=(G_i\,ⓩ\,H)^8$. From the third zig-zag bound, $1-\mu_i\ge\frac38(1-\lambda_i)$ and $\lambda_{i+1}=\mu_i^8\le\max(\lambda_i^2,\tfrac12)$, so after $k=O(\log n)$ rounds $G_k$ is an $(nd^{16k},d^{16},3/4)$-graph (Claim 9.4) whose neighbourhoods are computable in logspace (Claim 9.5). Since a connected $D$-regular graph already has $\alpha<1-\Omega(1/n^2)$, $O(\log n)$ doublings of the gap reach a constant; expanders have logarithmic diameter, so $s$–$t$ connectivity is in L.

**Interpretation.** The zig-zag proof is a *two-level* local-to-global theorem: global gap from (i) the gap of the coarse graph acting on cloud averages and (ii) the gap inside clouds acting on fluctuations, with the cross term controlled by $\beta$. This is the decomposition of a function into its conditional expectation onto a coarse algebra and the orthogonal fluctuation, i.e. a **two-projection (angle) estimate**. It is structurally the same estimate as Garland's method and trickling down (local-to-global note §2) and as Kastoryano–Brandão's Prop 20 (§9.8).

---

## 7. Applications: derandomisation, hardness, codes, embeddings

**Error reduction** (HLW §3.3.1). An RP algorithm with one-sided error $\beta$ using $k$ random bits: walk $t$ steps on an explicit $(2^k,d,\alpha)$-graph with $\alpha<\beta$ and take the conjunction; error $\le(\beta+\alpha)^t$ with $m+O(t)$ random bits (Thm 3.6). BPP (error $\le1/10$): majority over a walk; with $\alpha+\beta\le1/8$, $\Pr[\text{fail}]\le2^t(\beta+\alpha)^{(t-1)/2}=O(2^{-t/2})$ (Thm 3.10 and a union bound). HLW's table: $t$ independent repetitions give $2^{-t}$ with $tm$ bits; a walk on an $(n,d,1/40)$-graph gives $2^{-t/2}$ with $m+O(t)$ bits.

**Hardness of approximating clique** (HLW §3.3.2, Thm 3.13; [ALM98], proof following [AFWZ95]). The Berman–Schnitger randomised graph-product reduction is derandomised by taking all $t$-tuples that are length-$(t-1)$ walks on an $(n,d,\alpha)$-expander: $\omega(H')\le(\delta_2+2\alpha)^tm$ if $\omega(G)\le\delta_2n$ and $\omega(H')\ge(\delta_1-2\alpha)^tm$ if $\omega(G)\ge\delta_1n$ (Claims 3.15–3.16, via Thm 3.9). Approximating $\omega$ within $n^\epsilon$ is NP-hard.

**Codes** (HLW §12). For a $k$-left-regular bipartite graph $G$ ($n$ left "variables", right "checks"), $C(G)$ is the binary code with parity-check matrix the biadjacency matrix, and $L(G,d)$ the minimal left expansion of sets of size $\le d$ (definitions reconstructed from HLW's proofs; the formula text was garbled in the extraction).
- Sipser–Spielman (Thm 12.8): if $L(G,d)>k/2$ then $\mathrm{dist}(C(G))\ge d$. Proof: expansion $>k/2$ forces a **unique neighbour** of every set of size $\le d$, which witnesses a violated check.
- Sipser–Spielman (Thm 12.9): if $L(G,d)>\tfrac34k$ and $y$ is within $d/2$ of a codeword $x$, repeatedly flipping a variable with more unsatisfied than satisfied checks returns $x$ after a linear number of iterations. Proof: $|U|+|S|>\tfrac34k|A|$ and $|U|+2|S|\le k|A|$ give $|U|>\tfrac12k|A|$, so a flippable variable exists while errors remain, $|U|$ decreases, and $|A_i|$ never reaches $d$. The parallel version converges in $O(\log n)$ phases with lossless expanders [CRVW02].
- Gilbert–Varshamov (Thm 12.2) is the random benchmark.

**Metric embedding** (HLW §13). Linial–London–Rabinovich (Thm 13.8): for $k\ge3$, $\epsilon>0$, an $(n,k)$-graph with $\lambda_2\le k-\epsilon$ has $c_2(G)=\Omega(\log n)$ (least distortion into $\ell^2$), matched by Bourgain's upper bound. Thm 13.9 (the Poincaré inequality for embeddings): for $f:V\to\mathbb R^n$,
$$
\mathbb E_{(u,v)\in V\times V}\|f(u)-f(v)\|^2\ \le\ \frac{k}{k-\lambda_2}\ \mathbb E_{(u,v)\in E}\|f(u)-f(v)\|^2 .
$$
*Expanders are the metric spaces that Hilbert space cannot see at bounded distortion.* This is the hinge to coarse geometry and NCG (§9.11).

---

## 8. The mechanism in one paragraph: why spectral gap = local-to-global decorrelation

*Interpretation, assembled from the Facts above.* A step of the walk is local (it moves mass only along edges), yet a gap $1-\alpha$ says that after one step *every* function's deviation from its global mean has shrunk by the factor $\alpha$, uniformly over all functions orthogonal to the constants. Read probabilistically, $\alpha$ is the largest correlation any observable of the present can have with any observable of the next step: **Fact (memory)**, the Hirschfeld–Gebelein–Rényi maximal correlation of $(X_0,X_1)$ for a stationary reversible chain is the second largest singular value of the transition operator, here $\lambda(G)/d$ (Witsenhausen 1975). The expander mixing lemma is exactly this statement for indicator functions: the edge measure is $\alpha$-close to the product measure (§3). Decorrelation compounds ($\alpha^t$ after $t$ steps, Thm 3.3), and it survives local conditioning at an additive price ($\|P_B\hat AP_B\|\le\beta+\alpha$, Lemma 3.8), which is why a dependent walk samples like independent points. Cheeger says decorrelation can fail only through a **bottleneck**, a region whose boundary is small relative to its volume, i.e. a region that *remembers*. Alon–Boppana says decorrelation cannot beat the **universal cover**: locally every $d$-regular graph is a tree, on which eigenfunctions decay like $(d-1)^{-r/2}$ against $(d-1)^r$ growth; Ramanujan graphs decorrelate as fast as their own local geometry permits. So "local-to-global" here means: **a uniform bound on one-step correlations, which can be certified locally when the structure is hereditary (links, clouds, conditionals), propagates to all scales because composition multiplies the bound.** In the Gibbs world the same three statements reappear with "graph" replaced by "configuration space with single-site moves": gap $\Leftrightarrow$ no bottleneck (conductance), one-step decorrelation controlled by *influences* (Dobrushin is an $\ell^1$ contraction, Hayes and spectral independence are $\ell^2$ contractions), and the universal cover (the self-avoiding-walk tree) fixes the thresholds. What expander theory does **not** supply is the Markov property itself: conditional independence given a separator is static and exact (Hammersley–Clifford), and holds at every temperature. The Markov property is the *hereditary* ingredient that lets the decorrelation constant be checked on small pieces; the gap is the *global* conclusion.

---

## 9. The bridges that are theorems

Setting for §§9.1–9.7 (classical). A **spin system** on a finite graph $G=(V,E)$ with spins $[q]$: Gibbs distribution $\mu(\sigma)\propto\prod_{uv\in E}A(\sigma_u,\sigma_v)\prod_vh(\sigma_v)$ (Chen–Liu–Vigoda §1.1). Ising: $\mu(\sigma)\propto\exp(\sum_{uv}\beta_{uv}\sigma_u\sigma_v+\sum_vh_v\sigma_v)$ (Mossel–Sly Def 1), or $\nu_0\propto\exp(\tfrac12\langle x,Jx\rangle+\langle h,x\rangle)$ on $\{\pm1\}^n$ (EKZ eq. (1)). **Glauber dynamics** (Gibbs sampler, heat bath): pick $v$ uniformly, resample $\sigma_v$ from $\mu(\cdot\mid\sigma_{V\setminus v})$ (Mossel–Sly Def 2). It is reversible with respect to $\mu$. **Markov property**: the conditional law of $\sigma_v$ depends only on the neighbours; more generally $\mu(\sigma_A\mid\sigma_{A^c})$ depends only on $\sigma_{\partial A}$. **Fact (memory)**, Hammersley–Clifford: a strictly positive law is Markov with respect to $G$ iff it is a Gibbs law with clique potentials (Yang §I: "The Hammersley–Clifford theorem expresses this connection between positive Markov random fields and Gibbs distributions with clique-local interactions").

### 9.1 Spectral gap, Dirichlet form, relaxation and mixing

**Fact** (Mossel–Sly eqs. (5)–(6); EKZ eq. (4), Thm 10 = LPW Thm 20.6). For a reversible chain with stationary $\pi$, relaxation time $\tau=1/\mathrm{gap}$ has the variational form
$$
\tau=\sup_f\frac{2\sum_\sigma\pi(\sigma)f(\sigma)^2}{\sum_{\sigma\ne\tau}Q(\sigma,\tau)(f(\sigma)-f(\tau))^2},\quad \textstyle\sum\pi f=0,\qquad Q(\sigma,\tau)=\pi(\sigma)P(\sigma\to\tau),
$$
and $\tau\le\tau_{\rm mix}\le\tau\big(1+\tfrac12\log(\min_\sigma\pi(\sigma))^{-1}\big)$. For Glauber dynamics in continuous time the Dirichlet form is $\mathcal E_\nu(\varphi,\varphi)=\mathbb E_\nu\sum_i(\mathbb E_\nu[\varphi\mid X_{\sim i}]-\varphi)^2$ and a gap $\gamma$ is the Poincaré inequality $\mathrm{Var}(\varphi)\le\gamma^{-1}\mathcal E(\varphi,\varphi)$; then $\max_x\|P_t(x,\cdot)-\pi\|_{TV}\le\epsilon$ for $t\ge\gamma^{-1}\log(1/(\epsilon\min\pi))$. Continuous-time relaxation is $n$ times faster than discrete (Mossel–Sly §1.2.2; BKMP Def 1.1).

**Interpretation.** This is Tao's Ex 13 verbatim with the graph replaced by the configuration graph (Hamming graph restricted to the support) weighted by $Q$: *a Gibbs sampler is a random walk on a weighted graph whose vertices are configurations, and its gap is that graph's expansion.* Everything in §§1–3 applies to it; the difficulty is that this graph has $q^n$ vertices, so its expansion must be certified from the local structure of $\mu$.

### 9.2 Dobrushin uniqueness, path coupling and matrix norms

**Fact** (definitions; Dyer–Goldberg–Jerrum §2.2). For single-site updates $P^{[j]}$ preserving $\pi$, let $\mu_j(x,\cdot)$ be the law of the new spin at $j$. The **influence** of site $i$ on site $j$ is $\hat\varrho_{ij}=\max_{(x,y)\in S_i}d_{TV}(\mu_j(x,\cdot),\mu_j(y,\cdot))$, $S_i$ = pairs differing only at $i$. A **dependency matrix** is any $R\ge(\hat\varrho_{ij})$ entrywise. Random update $P^\dagger=\frac1n\sum_jP^{[j]}$; systematic scan $P^\to=\prod_jP^{[j]}$.

**Fact** (the conditions; DGJ §1). **Dobrushin**: $\|R\|_1<1$ (max column sum: total influence *on* a site); Dobrushin (1968/1970) showed this implies **uniqueness of the infinite-volume Gibbs measure**. **Dobrushin–Shlosman**: $\|R\|_\infty<1$ (max row sum); also implies uniqueness. **Hayes** (FOCS 2006): $\|R\|_2<1$ suffices for rapid mixing; for symmetric $R$ this is $\lambda(R)<1$.

**Fact** (DGJ Lemma 1, Lemma 17, Cor 18). If $\|R\|\le\mu<1$ in **any** matrix norm, random-update Glauber mixes in $\hat\tau_r(\epsilon)\sim n(1-\mu)^{-1}\ln((1-\mu)^{-1}J_n/\epsilon)$, $J_n$ the norm of the all-ones matrix; for $\|\cdot\|_1,\|\cdot\|_\infty,\|\cdot\|_p$: $\tau_r(\epsilon)\le n(1-\mu)^{-1}\ln(n/\epsilon)$. Systematic scan (Lemma 2): $\hat\tau_s(\epsilon)\sim(1-\mu)^{-1}\ln((1-\mu)^{-1}J_n/\epsilon)$ scans; Dobrushin $\Rightarrow\tau_s\le(1-\mu)^{-1}\ln(n/\epsilon)$ (Cor 25). For symmetric $R$ with zero diagonal and $\lambda(R)=\lambda<1$ (Lemma 3): $\hat\tau_s(\epsilon)\sim(1-\tfrac12\lambda)(1-\lambda)^{-1}\ln((1-\lambda)^{-1}n/\epsilon)$.

**Proof architecture** (DGJ §3.1, §3.2). **Path coupling** (Bubley–Dyer 1997): put weights $\delta_i$ on sites, use the path metric $d_\delta(x,y)=\sum_i\delta_i1[x_i\ne y_i]$, couple only pairs differing at one site (same site choice, maximal coupling of the new spin), and get $\beta_{t+1}\le\beta_tR^\dagger$ for the row vector of disagreement probabilities, $R^\dagger=\frac{n-1}nI+\frac1nR$, $\|R^\dagger\|\le1-(1-\mu)/n$ (Lemma 16); iterate and use the coupling lemma. **Dobrushin uniqueness**: the same estimate acting on the right, $\delta(P^\dagger f)\le R^\dagger\delta(f)$ for the vector of site-oscillations $\delta_i(f)$ (Lemma 30). DGJ: "the path-coupling approach is essentially equivalent to Dobrushin uniqueness"; row vectors vs column vectors are the two dual pictures (states forward, observables backward; compare Note 1 §5.3).

**Derivation (Ising influences).** For the heat-bath Ising update, $\mathbb P(\sigma_j=+\mid\cdot)=\frac{e^a}{e^a+e^{-a}}$ with $a=h_j+\sum_i\beta_{ij}\sigma_i$; flipping $\sigma_i$ shifts $a$ by $2\beta_{ij}$, and $\max_a\frac12|\tanh a-\tanh(a-2\beta)|=\tanh\beta$. So $\varrho_{ij}\le\tanh|\beta_{ij}|$, i.e. $R\le\tanh(\beta)A_G$ entrywise for uniform coupling $\beta$. Consequences: Dobrushin $\Leftarrow\Delta\tanh\beta<1$ ($\Delta$ = max degree); Hayes $\Leftarrow\lambda_{\max}(A_G)\tanh\beta<1$ (Perron–Frobenius monotonicity). **The largest adjacency eigenvalue of the interaction graph, not its degree, is the Hayes certificate**; for $d$-regular graphs they coincide, for irregular, planar or tree-like graphs Hayes gains (DGJ §1: Hayes applies a new estimate of the largest eigenvalue of planar graphs).

**Fact** (ℓ²-Dobrushin Poincaré inequality; EKZ Def 5 and Thm 6 = Wu, Ann. Probab. 2006, Thm 2.1). With the influence matrix $A$ (zero diagonal, $A_{ij}=\max\|P[X_i\in\cdot\mid x_{\sim i}]-P[X_i\in\cdot\mid x'_{\sim i}]\|_{TV}$ over $x,x'$ differing at $j$), $(1-\|A\|_{OP})\mathrm{Var}(\varphi)\le\mathcal E_\nu(\varphi,\varphi)$.

### 9.3 Spectral conditions on the interaction matrix; the Ramanujan regime

**Fact** (Bauerschmidt–Bodineau, JFA 2019, as quoted by EKZ). For $J\succeq0$ with $\|J\|_{OP}<1$ and any law $\rho$ on $\{\pm1\}^n$: $D(\rho\|\nu_0)\lesssim\frac1{1-\|J\|_{OP}}\sum_i\mathbb E_\nu|\partial_i\sqrt{d\rho/d\nu_0}|^2$ (an LSI with the hypercube gradient, which EKZ show can be $e^{\Theta(\beta\sqrt n)}$ weaker than the Glauber Dirichlet form on SK).

**Fact** (EKZ Thm 1). For $\nu_0$ with $0\preceq J\prec I$ and every $\varphi$,
$$
(1-\|J\|_{OP})\,\mathrm{Var}_{\nu_0}(\varphi)\ \le\ \mathcal E_{\nu_0}(\varphi,\varphi)
$$
(Glauber Dirichlet form). Since $\langle x,Jx\rangle$ changes only by a constant under $J\mapsto J+cI$ on $\{\pm1\}^n$, the hypothesis is effectively $\lambda_{\max}(J)-\lambda_{\min}(J)<1$ (EKZ: "we can clearly assume without loss of generality that $J$ is symmetric and positive definite"). Thm 11: continuous-time mixing for $t\ge\frac1{1-\|J\|}\big((1+2\|J\|)n+2|h|_1+\log\frac1\epsilon\big)$; discrete $O\big(\frac{n^2+|h|_1n+n\log(1/\epsilon)}{1-\|J\|}\big)$.

**Proof architecture** (EKZ §2). **Stochastic localization / needle decomposition**: an SDE $dF_t(x)=F_t(x)\langle C_t(x-a_t),dW_t\rangle$, $dJ_t/dt=-C_t^2$, drives $\nu_0$ through a martingale of measures $\nu_t$ with quadratic part $J_t=J-\int_0^tC_s^2ds$ until rank $\le1$ (stopping time $T\le\frac12\mathrm{Tr}J$), keeping $\int\varphi\,d\nu_t$ nearly constant ($C_t$ almost kills the direction $V_t$). Two facts close it: rank-one Ising measures have influence matrix norm $\le|u|^2\le\|J\|$ (Lemma 7, via $\tanh$ 1-Lipschitz and Perron–Frobenius), hence Wu's Poincaré inequality (Lemma 8); and the Glauber Dirichlet form is a **supermartingale** along the localization (Lemma 9, by convexity of $|\cdot|^2$ in the conductance formula $\mathcal E=\sum_{x\sim y}\frac{\nu(x)\nu(y)}{\nu(x)+\nu(y)}(\varphi(x)-\varphi(y))^2$). Thm 12: every such $\nu_0$ is a mixture of rank-one models $w_{u,v}$ with $|u|\le\|J\|$, preserving $\int\varphi$ and not increasing conductances.

**Fact (the Ramanujan bridge; EKZ §5).** SK model ($J_{ij}\sim N(0,\beta^2/n)$): spectrum of $J$ in $[-2\beta-\epsilon,2\beta+\epsilon]$ a.a.s., so Poincaré and polynomial mixing for all $\beta<1/4$, whereas Dobrushin needs $\beta=O(1/\sqrt n)$. **Diluted SK on a random $d$-regular graph** with Rademacher couplings $\pm\beta$: "it follows from a version of Friedman's Theorem that $\|J\|_{OP}\le\beta(2\sqrt{d-1}+\epsilon)$ a.a.s.", giving the Poincaré inequality and fast mixing for all $\beta<\frac1{4\sqrt{d-1}}$, "whereas the model is only in Dobrushin's uniqueness regime for $\beta=O(1/d)$ — note that up to constants the latter bound is tight for general Ising models on arbitrary $d$-regular graphs" (citing Galanis–Štefankovič–Vigoda and Sly–Sun).

**Interpretation.** This is the sharpest theorem-level answer to "expansion of the interaction graph and Gibbs mixing": for signed (frustrated) interactions, *the Ramanujan property of the interaction graph enlarges the provable high-temperature regime by a factor $\sqrt d$*, because the relevant certificate is a spectral norm and Alon–Boppana/Friedman fix the spectral norm of a $d$-regular signed adjacency at $2\sqrt{d-1}$ rather than $d$. For ferromagnetic couplings nothing is gained: **Derivation**, for $J=\beta A_G$ on a $d$-regular graph, $\lambda_{\max}(J)-\lambda_{\min}(J)=\beta(d-\lambda_{\min}(A_G))\ge\beta d$, so the condition is no better than Dobrushin's $\Delta\tanh\beta<1$ up to constants.

### 9.4 Exact thresholds on general graphs: the universal cover decides

**Fact** (Mossel–Sly Thm 1). For $d\ge2$ and $\beta>0$ with
$$
(d-1)\tanh\beta<1,
$$
there are $0<\lambda^*(d,\beta)$, $C(d,\beta)<\infty$ such that on *any* graph of maximum degree $d$ on $n$ vertices the discrete-time mixing time of the Gibbs sampler for the ferromagnetic Ising model with all couplings $\le\beta$ and **arbitrary external fields** is $\le Cn\log n$, and the continuous-time gap is $\ge\lambda^*$, also on infinite graphs of maximum degree $d$. **Thm 2**: on $G(n,d/n)$ with $d\tanh\beta<1$, w.h.p. $n^{1+c/\log\log n}\le\tau_{\rm mix}\le n^{1+C/\log\log n}$. Both are **tight**: on random $d$-regular graphs (resp. $G(n,d/n)$) with no field, mixing is w.h.p. $\exp(\Omega(n))$ when $(d-1)\tanh\beta>1$ (resp. $d\tanh\beta>1$).

**Fact** (meaning of the threshold; Mossel–Sly §1). $(d-1)\tanh\beta\le1$ is the uniqueness threshold of the ferromagnetic Ising model on the infinite $d$-regular tree (Lyons); Weitz showed spatial mixing for $(d-1)\tanh\beta\le1$ on any graph of maximum degree $d$; $d\tanh\beta<1$ is uniqueness on the Poisson($d$) Galton–Watson tree. Mossel–Sly: "our results are the first rigorous results establishing exact thresholds for dynamics on random graphs in terms of spatial thresholds on trees", and **"since most graphs of bounded degree are expanders, the strong spatial mixing technique does not apply to them"**; block dynamics "cannot be extended to the nonamenable setting since the bounds rely crucially on the small boundary-to-volume ratio".

**Proof architecture** (Mossel–Sly §§1.3–2, 5). A general theorem from three local conditions on balls $B(v,R)$: volume $\mathrm{Vol}(R,X)$, local mixing $\mathrm{LM}(R,T)$, spatial mixing $\mathrm{SM}(R)$: $\sum_{u\in S(v,R)}a_u\le\frac14$ with $a_u$ the effect on $\sigma_v$ of the boundary spin at $u$. Monotone coupling from all-$+$ and all-$-$; **censoring** (Peres–Winkler: for monotone systems, omitting updates only slows coupling) restricts updates to $B(v,R)$; the triangle inequality splits the disagreement at $v$ into local mixing error and boundary effect, giving $\max_u\Pr[X^+_{S+t}(u)\ne X^-_{S+t}(u)]\le\frac12\max_u\Pr[X^+_S(u)\ne X^-_S(u)]$, hence a constant gap. SM via **Weitz's self-avoiding-walk tree** (Lemma 13: marginals on $G$ equal marginals on $T_{\rm saw}(G,v)$ with fixed leaves) and tree recursions: $\sum_ua_u\le\sum_{\ell\ge R}d(d-1)^{\ell-1}\tanh^\ell\beta=\frac{d(d-1)^{R-1}\tanh^R\beta}{1-(d-1)\tanh\beta}\le\frac14$ (Lemma 3). LM via cut-width (§9.6).

**Interpretation.** The tree of self-avoiding walks is a universal cover with boundary conditions. *Alon–Boppana says the spectrum of a regular graph cannot beat its universal cover; Mossel–Sly say the Gibbs sampler on a bounded-degree graph mixes exactly as well as the uniqueness of its universal-cover tree allows.* The common mechanism is branching $(d-1)$ against per-edge decay ($\tanh\beta$ for Gibbs influence, $1/\sqrt{d-1}$ per edge for $\ell^2$ eigenfunctions).

### 9.5 Spectral independence: correlation spectrum = link spectrum

**Fact** (definitions; Chen–Liu–Vigoda Defs 1.9–1.11). $\mu$ on $[q]^V$ is **$b$-marginally bounded** if $\mu(\sigma_v=i\mid\sigma_\Lambda=\tau)\ge b$ for every pinning and feasible $i$. Influence matrix under pinning $\tau$: $\Psi^\tau_\mu((u,i),(v,j))=\mu(\sigma_v=j\mid\sigma_u=i,\sigma_\Lambda=\tau)-\mu(\sigma_v=j\mid\sigma_\Lambda=\tau)$ for $u\ne v$, zero on the diagonal blocks; its eigenvalues are real. $\mu$ is **$\eta$-spectrally independent** if $\lambda_1(\Psi^\tau_\mu)\le\eta$ for **every** pinning $\tau$ (Anari–Liu–Oveis Gharan 2020). Assumption 1.8: $\mu$ is **totally connected** (every pinned support is connected under single-site changes).

**Fact** (CLV Claim 1.18, from Chen–Galanis–Štefankovič–Vigoda Thm 8). $\eta$-spectral independence implies that the weighted simplicial complex of partial configurations is a $(\zeta_0,\dots,\zeta_{n-2})$-**local spectral expander** with $\zeta_k=\eta/(n-k-1)$. *The influence matrix under pinning $\tau$ is, up to normalisation, the non-trivial part of the random walk on the link of $\tau$.*

**Fact** (CLV Thm 1.12, Rem 1.13). If $G$ has maximum degree $\le\Delta$ and $\mu$ is totally connected, $b$-marginally bounded and $\eta$-spectrally independent, Glauber dynamics satisfies the **modified log-Sobolev inequality** with constant $1/(C_1n)$, $C_1=(\Delta/b)^{O(\eta/b^2+1)}$, and $T_{\rm mix}(\epsilon)=(\Delta/b)^{O(\eta/b^2+1)}\,O(n\log(n/\epsilon))$. Explicitly, for $n\ge\frac{24\Delta}{b^2}(\frac{4\eta}{b^2}+1)$: $C_1=\frac{18\log(1/b)}{b^4}\big(\frac{24\Delta}{b^2}\big)^{4\eta/b^2+1}$ and $T_{\rm mix}\le\lceil C_1n(\log n+\log\log\frac1b+\log\frac1{2\epsilon^2})\rceil$.

**Fact** (CLV Thm 1.19, the local-to-global step). For a pure $n$-dimensional weighted simplicial complex that is $(b_0..b_{n-1})$-marginally bounded with $(\zeta_0..\zeta_{n-2})$-local spectral expansion, the order-$(s,r)$ down-up and order-$(r,s)$ up-down walks satisfy MLSI with $\kappa=\sum_{k=r}^{s-1}\Gamma_k/\sum_{k=0}^{s-1}\Gamma_k$, $\Gamma_0=1$, $\Gamma_k=\prod_{j<k}\alpha_j$, $\alpha_k=\max\{1-\frac{4\zeta_k}{b_k^2(s-k)^2},\frac{1-\zeta_k}{4+2\log(1/(2b_kb_{k+1}))}\}$; $T_{\rm mix}\le\lceil\kappa^{-1}(\log\log\frac1{\pi^*_s}+\log\frac1{2\epsilon^2})\rceil$. Proof architecture (CLV §2): **approximate tensorization of entropy** $\mathrm{Ent}(f)\le C_1\sum_v\mu[\mathrm{Ent}_v(f)]$ (Def 2.1) is deduced from **uniform block factorization** into blocks of linear size $\ell=\theta n$ (Lemma 2.3), which is proved by the Alev–Lau local-to-global scheme adapted to entropy (Thm 2.9).

**Fact** (applications; CLV §1). $O(n\log n)$ mixing on bounded-degree graphs for: the hard-core model for $\lambda<\lambda_c(\Delta)=\frac{(\Delta-1)^{\Delta-1}}{(\Delta-2)^\Delta}$ (the uniqueness threshold on the $\Delta$-regular tree); all antiferromagnetic 2-spin systems that are up-to-$\Delta$ unique with gap $\delta$; $q$-colourings of triangle-free graphs for $q\ge(\alpha^*+\delta)\Delta$, $\alpha^*\approx1.763$; matchings in $O(m\log n)$ with $\eta=\min\{2\lambda\Delta,2\sqrt{1+\lambda\Delta}\}$ (Thms 2.10, 6.1–6.4, the edge-influence bound proved on the SAW tree). Complexity (CLV §1, citing Weitz and Sly): Weitz's correlation-decay FPTAS for $\lambda<(1-\delta)\lambda_c(\Delta)$; for $\lambda>\lambda_c$ no FPRAS unless NP = RP (Sly 2010; Sly–Sun; Galanis–Štefankovič–Vigoda).

**Interpretation.** This is the theorem the user's intuition was reaching for, and it is already in the programme's local-to-global note as U6. *Decay of correlations under every pinning (a Gibbs/MRF statement) is literally the local spectral expansion of the complex of pinnings (an HDX statement), and the HDX local-to-global theorem turns it into a global log-Sobolev inequality.* The MRF structure is used through hereditarity: a pinned Gibbs measure is a Gibbs measure of the same kind on the remaining sites, so "every link" means "every pinning".

### 9.6 Trees, hyperbolic graphs, cut-width; spectral gap implies decay of correlations

**Fact** (BKMP Thm 1.4; continuous time; $\epsilon=(1+e^{2\beta})^{-1}$, so $1-2\epsilon=\tanh\beta$). Ising on the $b$-ary tree $T_r$ with $n_r$ vertices, free boundary:
1. at all temperatures $\tau_2=n_r^{O(\log(1/\epsilon))}$, and $\lim_r\log\tau_2/\log n_r$ exists;
2. if $\tanh\beta>1/\sqrt b$: $\tau_2=\Theta\big(n_r^{\log_b(b\tanh^2\beta)}\big)$ (unbounded), $\Theta(\log n_r)$ at equality; $\tau_2=n_r^{\Omega(\log(1/\epsilon))}$ as $\epsilon\to0$;
3. if $\tanh\beta<1/\sqrt b$: $\tau_2=O(1)$, for every external field.

The tree's Gibbs regimes: uniqueness for $\tanh\beta<1/b$; for $1/b<\tanh\beta<1/\sqrt b$ infinitely many Gibbs measures yet bounded relaxation time; reconstruction (typical boundaries bias the root) for $\tanh\beta>1/\sqrt b$. **The gap of Glauber dynamics is governed by the Kesten–Stigum (census-reconstruction) threshold $b\tanh^2\beta=1$, not by uniqueness $b\tanh\beta=1$.** Lower bounds come from low-conductance cuts of configuration space (global majority of boundary spins; recursive majority for very low temperature, §3); upper bounds from block dynamics plus path coupling with the weighted Hamming metric $\sum_v\theta^{|v|}1(\sigma_v\ne\eta_v)$, $\theta=1/\sqrt b$ (§4).

**Fact** (cut-width; BKMP Def 1.4, Prop 1.1, Lemma 2.1, Props 1.2–1.3). Cut-width $\xi(G)$: least $\max_k|E(\{v_1..v_k\},\{v_{k+1}..v_n\})|$ over orderings. Ising: $\tau_2\le n\,e^{(4\xi(G)+2\Delta)\beta}$ (canonical paths along the ordering); colourings with $q\ge\Delta+2$: $\tau_2\le(\Delta+1)n(q-1)^{\xi(G)+1}$. $\xi(T_r^{(b)})<(b-1)r+1$. Planar graphs with Cheeger constant $\ge c>0$, degree $\le\Delta$ and no separating cycles in balls: $\xi(G_r)\le C\log n_r$, so polynomial mixing at every temperature on balls of hyperbolic tilings; yet (Prop 1.3) at low enough temperature such graphs have long-range correlations $\mathbb E[\sigma_u\sigma_v]\ge\delta$. **Polynomial mixing coexists with long-range order on non-amenable planar graphs.** Problem 4 asks for which graphs a lower bound $\tau_2\ge e^{c\xi(G)}$ holds at low temperature: "Such a lower bound is known to hold for boxes in a Euclidean lattice, our results imply its validity for regular trees, and we can also verify it for expander graphs." (Expanders have linear cut-width, so this gives $e^{\Omega(n)}$ at low temperature.)

**Fact** (BKMP Thm 1.5; general bounded-degree graphs). If $\tau_2(G_r)=O(1)$ then for every fixed finite $A$ there is $c_A>0$ with $\mathrm{Cov}[f,g]\le e^{-c_Ar}\sqrt{\mathrm{Var}f\,\mathrm{Var}g}$ whenever $f$ depends on $\sigma_A$ and $g$ on the spins at distance $r$; equivalently $I[\sigma_A;\sigma_r]\le e^{-c'_Ar}$. Proof: disagreement percolation (van den Berg; Zegarlinski) and first-passage percolation with a Peierls bound. It holds even with multiple Gibbs measures. Problem 3 asks for the converse (fails in some lattices with plus boundary, per Martinelli).

**Fact** (gap ⇒ bounded covariance; EKZ Rem 14). A Poincaré inequality with constant $1-\|J\|$ applied to linear functions gives $\mathrm{Var}\langle w,X\rangle\le|w|^2/(1-\|J\|)$, i.e. $\|\Sigma\|_{OP}=O(1)$ for the covariance matrix.

**Interpretation (the $\ell^2$ threshold is one threshold).** The Kesten–Stigum condition $b\lambda_2^2<1$ for the edge channel ($\lambda_2=\tanh\beta$ is the second eigenvalue of the binary symmetric channel with flip probability $\epsilon$) and the Alon–Boppana edge $\lambda=2\sqrt{d-1}$ (where $(d-1)|\rho|^2=1$) are both the point where branching balances the *square* of the per-edge decay, i.e. where an $\ell^2$ (second-moment) quantity stops being summable over a tree. Uniqueness ($b\lambda_2<1$) is the $\ell^1$ (worst-case) threshold. Dobrushin is $\ell^1$; Hayes and spectral independence are $\ell^2$. *This explains why spectral (ℓ²) certificates reach further than Dobrushin's (ℓ¹) and why the gap on the tree tracks Kesten–Stigum.*

### 9.7 Expansion is the enemy at low temperature

**Fact** (Mossel–Sly §1, citing their refs [4, 7, 26], not opened here; CLV §1 citing Sly 2010 and follow-ups; EKZ §5 citing Galanis–Štefankovič–Vigoda, CPC 2016, and Sly–Sun, FOCS 2012). On random $d$-regular graphs the Ising Glauber dynamics without field mixes in $\exp(\Omega(n))$ w.h.p. when $(d-1)\tanh\beta>1$; for the hard-core and antiferromagnetic Ising models, mixing on almost all random $d$-regular *bipartite* graphs is exponential beyond the tree uniqueness threshold, and approximate counting is NP-hard there (unless NP = RP). BKMP Problem 4: the low-temperature lower bound $e^{c\xi(G)}$ is verified for expanders.

**Interpretation.** The mechanism is the one Cheeger's inequality names, applied in configuration space: in an expander every balanced cut of the *interaction* graph has $\Omega(n)$ edges, so every configuration interpolating between the two ordered phases pays energy $\Omega(\beta n)$, the magnetisation (or sublattice-occupation) cut of configuration space has conductance $e^{-\Omega(n)}$, and the conductance bound (§2) gives an exponentially small gap. Random bipartite regular graphs are the hardness gadgets *because* they are expanders: expansion forces phase coexistence to be global. So the interaction graph's expansion and the configuration-space walk's expansion are **anti-correlated at low temperature** and **positively related at high temperature only through spectral norms** (§9.3).

### 9.8 Amenable lattices, and the quantum commuting case: gap ⇔ strong clustering

**Fact (as reported).** On $\mathbb Z^d$ (amenable), strong spatial mixing implies a uniform spectral gap and log-Sobolev inequality for Glauber dynamics; Chen–Rouzé §I.B summarise Martinelli's Saint-Flour notes as "(Exact Markov) + (Strong spatial mixing) ⟹ MCMC mixes in quasi-linear time", strong mixing being $\sup_\tau\|\rho^{B,\tau}_A-\rho^{B,\tau^x}_A\|_{TV}\le Ce^{-\mathrm{dist}(A,x)/\xi_s}$. Mossel–Sly: this is "only known for amenable graphs and for a strong form of spatial mixing". **Fact (memory):** Stroock–Zegarlinski 1992 and Martinelli–Olivieri 1994 for the equivalences on $\mathbb Z^d$.

**Fact** (Kastoryano–Brandão, CMP 2016). Setting: $\Lambda\subset\mathbb Z^d$, $r$-local bounded **commuting** potential, Gibbs state $\rho=e^{-\beta H}/\mathrm{tr}$. Conditional expectations $\mathbb E_A$ (two kinds: the Davies/Liouvillian projector $\mathbb E^{\mathcal L}_A$ and the "minimal" static one $\mathbb E^\rho_A$; $\mathbb E_A(f)$ is supported on $A^c$). Conditional covariance $\mathrm{Cov}_A(f,g)=|\langle f-\mathbb E_Af,g-\mathbb E_Ag\rangle_\rho|$.
- Def 14 (**weak clustering**): $\mathrm{Cov}(f,g)\le c\|f\|_{2,\rho}\|g\|_{2,\rho}e^{-d(\Sigma_f,\Sigma_g)/\xi}$.
- Def 15 (**strong clustering**): for $A\cap B\ne\emptyset$, $\mathrm{Cov}_{A\cup B}(\mathbb E_A f,\mathbb E_B f)\le c\|f\|^2_{2,\rho}e^{-d(B\setminus A,A\setminus B)/\xi}$. It "incorporates the strong mixing (or complete analyticity) condition for classical systems".
- **Prop 20** (general, no Gibbs or locality assumption used): if $\mathrm{Cov}_{A\cup B}(\mathbb E_Af,\mathbb E_Bf)\le\epsilon\,\mathrm{Var}_{A\cup B}(f)$ for all $f$, then
$$
\mathrm{Var}_{A\cup B}(f)\le(1-2\epsilon)^{-1}\big(\mathrm{Var}_A(f)+\mathrm{Var}_B(f)\big).
$$
Proof: expand $0\le\|(\mathrm{id}-\mathbb E_{A\cup B})(\mathrm{id}-\mathbb E_A-\mathbb E_B)f\|^2$ and use contractivity.
- **Thm 23** (§VI.A, strong clustering ⇒ gap): on $\Lambda=(\mathbb Z/l)^d$, $l\ge l_0$, strong clustering w.r.t. $\mathbb E^{\mathcal L}$ implies the Davies generator has a gap independent of $|\Lambda|$. Architecture: the conditional gap $\lambda_\Lambda(A)=\inf_f\langle f,-\mathcal L_A f\rangle_\rho/\mathrm{Var}_A(f)$ is unchanged when $A$ is doubled to $A\cup B$ with overlap of side $\ge\sqrt L$ (Prop 20 with $\epsilon=ce^{-\sqrt L/\xi}$), iterated over scales; Lemma 22 localises the conditional gap to $A_\partial$. Thm 26 (§VI.B, gap ⇒ strong clustering): the introduction attributes one direction to a scale-doubling argument "reminiscent of the analogous classical result" and the other to a map onto frustration-free Hamiltonians and the detectability lemma; since §VI.A visibly uses the scale-doubling, the detectability lemma belongs to §VI.B.
- **Thm 28, Prop 29**: in 1D weak ⇔ strong clustering; all 1D commuting Gibbs samplers are gapped (Araki's clustering, MPS transfer-operator gap). **Thms 30–31**: above a size-independent temperature $T_c(r,d)$ both samplers are gapped (Knabe-type argument).

**Interpretation.** Prop 20 is the zig-zag estimate (§6) and AKS Lemma 3.8 (§3) in conditional-expectation language: the global variance is controlled by local variances provided the two conditional expectations are at a definite **angle** (their "cross-covariance" is a fraction $\epsilon<1/2$ of the variance). Strong clustering is the statement that this angle defect decays with the distance between $A\setminus B$ and $B\setminus A$. **Fact (memory):** in subfactor theory two conditional expectations with $\mathbb E_A\mathbb E_B=\mathbb E_B\mathbb E_A=\mathbb E_{A\cap B}$ form a *commuting square* (Popa). *Strong clustering reads as an approximate commuting-square condition with defect $e^{-d/\xi}$, and the gap theorem as "approximate commuting squares at every scale ⇒ global gap".* That is my interpretation; it is the operator-algebraic form of "Markov + decorrelation ⇒ gap".

### 9.9 Noncommuting quantum Gibbs states: the single-Pauli-jump Lindbladian

**Fact** (Chen–Rouzé 2504.02208, Thm III.1). Hamiltonian $H=\sum_\gamma H_\gamma$ of interaction degree $\le d$, Gibbs state $\rho_\beta$, region $A$, $\rho_{\beta,-A}=\mathrm{Tr}_A[\rho_\beta]\otimes\tau_A$ ($\tau_A$ maximally mixed). Let
$$
\mathcal R_{A,t}[\cdot]=\frac1t\int_0^t\exp(s\mathcal L_A)[\cdot]\,ds,\qquad \mathcal L_A=\sum_{a\in P^1_A}\mathcal L_a,\quad P^1_A=\{X_i,Y_i,Z_i\}_{i\in A},
$$
each $\mathcal L_a$ the exactly detailed-balanced (KMS) Lindbladian with the single-qubit Pauli jump $A_a$ and Metropolis weight. With $\beta_0=1/(4d)$,
$$
\|\mathcal R_{A,t}[\rho_{\beta,-A}]-\rho_\beta\|_1\le|A|^2\,2^{2|A|}\times\begin{cases}r(\beta,d)\,t^{-\frac{128\beta_0^4}{\beta^3(\beta+5\beta_0)}}&\beta>4\beta_0,\\ r'(\beta,d)\,t^{-\frac{2\beta_0}{\beta+5\beta_0}}&\beta\le4\beta_0,\end{cases}
$$
so $\le re^{\mu|A|}t^{-\lambda}$ with $r,\mu>0$, $0<\lambda<1$ depending on $\beta,d$. Fig. 2's caption: "the recovery map is a time-averaged detailed-balanced Lindbladian based on single-Pauli jumps on A."

**Proof architecture** (Chen–Rouzé §IV, §VII). (i) **Gap-free decay**: for the time average, the Dirichlet form obeys $\mathcal E_A(\mathcal R^\dagger_{A,t}[X])\le2/t$ unconditionally (Cor VII.1); "time-averaging provides a different mechanism to obtain a small Dirichlet form that is independent of the spectral gap"; "such a property is generally false for the Lindblad evolution itself without time-averaging, because $\mathcal L^\dagger_A$ may have arbitrarily small eigenvalues". (ii) The Dirichlet form controls commutators with the jumps in the KMS-weighted norm (a commutator-square formula; a Hölder-like inequality; regularised operator Fourier transforms at low temperature), so a slowly changing operator nearly commutes with every $A_a$, hence (after reducing high-weight Paulis to single Paulis, Cor VIII.1) is nearly independent of $A$. (iii) **Quasi-locality**: by Lieb–Robinson, $\|\mathcal L^\dagger_{A,\ell}-\mathcal L^\dagger_A\|_{\infty\to\infty}\lesssim|A|(e^{-c'\ell/(d\beta)}+2^{-\ell})$ for $\ell\ge4e^2\beta d$ (Lemma VII.3), and truncation errors in $\mathcal R$ grow linearly in $t$ (Lemma VII.2). (iv) Choosing $t=e^{(\mu|A|+m\ell)/\lambda+1}$ gives a quasi-local recovery map (Cor III.1) and $I(A{:}C|B)\lesssim\log(\dim C)\sqrt\Delta$ (Cor III.2). Discussion §V: "the bound on CMI grows exponentially with the size $|A|$, which comes from the possibility of an exponentially long mixing time"; a local gap of the $A$-Lindbladians decaying polynomially in $|A|$ would improve it.

**Derivation (why the time average is gap-free).** For any self-adjoint generator $\mathcal L\le0$ on a Hilbert space (here the KMS inner product, by detailed balance), with spectral variable $x\ge0$ for $-\mathcal L$: $\mathcal E(\mathcal R_tX)=\int\big(\frac{1-e^{-tx}}{tx}\big)^2x\,d\mu_X(x)\le\frac{c}{t}\|X\|^2$ with $c=\sup_{u>0}(1-e^{-u})^2/u\approx0.407$ (attained at $u\approx1.26$; computed numerically here). This is the von Neumann mean ergodic theorem with a rate in Dirichlet form; it uses only self-adjointness, i.e. **detailed balance**, which is why KMS detailed balance is essential. With a gap $\gamma$ one would instead get $\|e^{t\mathcal L}X-\mathbb EX\|\le e^{-\gamma t}\|X\|$: the expander route.

**Fact** (Yang 2609.38007, Thm II.1 and Table I). For finite-range interactions of bounded degree at every $\beta>0$: $I(A{:}C|B)_\rho\le C_\beta e^{C_\beta g_{A|A^c}-c_\beta r}$, $g_{A|A^c}$ the interaction strength across the cut ($\lesssim|\partial_eA|$ on $\mathbb Z^D$), improving Chen–Rouzé's $C_\beta|A||C|e^{C_\beta\min\{|A|,|C|\}-c_\beta r}$. **"Our proof is static, in contrast to the dynamical approach of Chen and Rouzé"**: compare the Gibbs state with the cut Gibbs state (boundary interactions removed), reduce CMI to a relative-entropy loss $\delta_{AB}(\rho,\sigma)$, and control it through a resolvent formula for the relative modular operator $\Delta=L_\sigma R_{\rho^{-1}}$ and moment/tail bounds on $K=\log\Delta$ (Lemmas III.2, V.3–V.4). Table I also records: Brown–Poulin (commuting Gibbs states are exactly Markov), Kato–Brandão ($C_\beta e^{-c_\beta\sqrt r}$ in 1D), Kuwahara, Kato–Kuwahara (CMI decay under rapid mixing or uniform clustering), **Bakshi–Liu–Moitra–Tang ("A Dobrushin condition for quantum Markov chains", arXiv 2510.08542): $C_\beta|A||C|e^{-c_\beta r}$ at high temperature via a quantum Dobrushin condition and rapid mixing**, and Rosa-Ruiz–Scandi–Capel–Alhambra (fixed points of rapidly mixing local Lindbladians). Classical baseline (Yang §I): "a classical Gibbs distribution retains its Markov property even at a thermal phase transition"; vanishing CMI ⇔ exact recovery (Petz; Hayden–Jozsa–Petz–Winter), and Fawzi–Renner make it quantitative.

**Interpretation.** Classically, running heat-bath Glauber dynamics only on $A$ with the outside frozen converges to $\mu(\cdot\mid\sigma_{A^c})$, which by the Markov property depends only on $\sigma_{\partial A}$: the dynamics on $A$ *is* an implementation of the conditional expectation onto the outside, and its locality is exact. In the noncommuting case there is no explicit local conditional expectation (Kastoryano–Brandão's framework fails, their §IX), so Chen–Rouzé *manufacture* one by time-averaging a local detailed-balanced generator, with locality supplied by Lieb–Robinson and convergence supplied by the gap-free mean ergodic theorem. **The expander analogue sits in the missing ingredient**: a uniform local gap for $\mathcal L_A$, which is a *trickling-up* question (single-site generators trivially gapped ⇒ gap of $\mathcal L_A$ for large $A$?). Classically the answer is "yes, under Dobrushin / spectral independence / strong clustering"; quantumly Bakshi–Liu–Moitra–Tang give it at high temperature. Yang's static route shows that for the *Markov* statement the dynamics can be bypassed altogether, with the cut size $|\partial A|$ (an isoperimetric quantity) in the exponent: the boundary-to-volume ratio, i.e. amenability, re-enters exactly as in §9.4.

### 9.10 Quantum expanders

**Fact** (Hastings 2007, after Ben-Aroya–Schwartz–Ta-Shma and Hastings's earlier paper). A **quantum expander** is a CP trace-preserving unital map $\mathcal E(M)=\sum_{s=1}^DA(s)^\dagger MA(s)$ on $N\times N$ matrices with $N$ large, $D$ small, eigenvalue $1$ (eigenvector $\mathbb 1/\sqrt N$) and $|\lambda_a|\le1-\delta$ for $a>1$. For $A(s)=U(s)/\sqrt D$ with $U(1..D/2)$ Haar-random and $U(s+D/2)=U(s)^\dagger$ (Hermitian case), $|\lambda_2|\to\lambda_H=\frac{2\sqrt{D-1}}{D}$ in probability as $N\to\infty$: "the same as the recently proven tight bound [Friedman] in the classical case". **Quantum Alon–Boppana** (eq. (12)), valid for every such map: $|\lambda_2|\ge\lambda_H(1-O(\ln\ln N/\ln N))$. Proof: trace method with Schwinger–Dyson equations for averages over $U(N)$. Pisier [Pis14] gave "a different point of view on Hastings' quantum expander result" in operator-space terms (as described in arXiv 1811.08847; Pisier's paper itself was not opened, and its arXiv number is not given here). Random isometry channels are (generalised) quantum expanders for Kraus rank $k\ge169$ (arXiv 1811.08847, *On the spectral gap of random quantum channels*, abstract and §7 seen; authors not checked here).

**Interpretation.** A quantum expander is a "noncommutative Cayley graph": the state space $M_N$ replaces $\ell^2(V)$, conjugation by unitaries replaces edges, and the gap is again isolation of the trivial (maximally mixed) eigenvector. Gibbs-sampler Lindbladians are *non-unital, detailed-balanced* channels whose fixed point is $\rho_\beta$ rather than $\mathbb 1/N$; their gap is the quantum analogue of the Glauber gap, and the KMS inner product replaces $\ell^2(\pi)$.

### 9.11 Noncommutative geometry: property (T), Roe algebras, ghosts

**Fact** (HLW §11, §2.2). Kazhdan's property (T) of $SL_3(\mathbb Z)$ gives a uniform gap for the Cayley graphs of all its finite quotients with a fixed generating set (Margulis 1973, existential); property $(\tau)$ via Selberg's $3/16$ for $SL_2$ quotients (Lubotzky); Ramanujan bounds via Deligne. $K(H,S)/(2|S|)<g(H,S)<K(H,S)/2$ (§1.4).

**Fact** (Willett–Yu I, Ex 5.3 and Cor 6.2; abstracts of I and II). For a space of graphs $X=\sqcup G_n$ that is an expander, the graph Laplacian $\Delta\in\mathbb C[X]$ has spectrum in $\{0\}\cup[c,2]$, so the spectral projection at $0$ — the **basic Kazhdan projection** $p=\prod_np^{(n)}$, $p^{(n)}$ the projection onto constants on $G_n$ — lies in the Roe algebra $C^*(X)$ ("the limit exists in the norm topology using the 'spectral gap' of $\Delta$"). Its matrix entries are $1/|G_n|\to0$, so it is a **ghost**: "non-compact ghost operators have a definite global existence (as non-compact), while simultaneously being 'locally almost invisible'". If moreover girth $\to\infty$: the coarse Baum–Connes assembly map for $X$ is **injective but not surjective** ($[p]\notin$ image); the maximal coarse assembly map is an isomorphism (paper II), and "geometric property (T)" obstructs it. Gromov monster groups (coarsely containing such expanders) have Baum–Connes assembly maps with certain coefficients injective but not surjective, with the maximal version an isomorphism. Earlier: Higson (non-surjectivity for Margulis-type expanders); Higson–Lafforgue–Skandalis (for any expander, either coarse Baum–Connes fails surjectivity or Baum–Connes with coefficients for an associated groupoid fails injectivity). Yu (2000): coarse Baum–Connes holds for bounded-geometry spaces that coarsely embed in Hilbert space (as stated in arXiv 2511.22438's introduction); expanders do not coarsely embed (HLW Thm 13.8 at the level of distortion).

**Interpretation (the NCG essence of a spectral gap).** In the Roe algebra, "local" means *finite propagation*. A spectral gap is exactly what lets the global projection onto the invariant (constant) functions be written as a norm limit $f(\Delta)$ of polynomials in the Laplacian, i.e. of finite-propagation operators, with $f(0)=1$, $f|_{[c,2]}=0$. So **the gap makes the global average a limit of local steps** (Reduction 1 of the local-to-global note, now as an operator-algebra statement), while the ghost property says **no finite window sees it**. The coarse Baum–Connes conjecture predicts that $K$-theory of $C^*(X)$ is assembled from local data; expanders produce a class that is assembled by local operators but not from local $K$-homology. This is the precise sense in which expanders are "local-to-global at the level of analysis but not at the level of topology". It is my reading, but every ingredient is a theorem.

---

## 10. One table: the common ingredients across four theories

*Interpretation, aligned with the "common mechanism" table of [local-to-global-unlocks.md §4](../../local-to-global-unlocks.md).*

| ingredient | expanders | Gibbs / Glauber | quantum Gibbs samplers | NCG (Roe / group algebras) |
|---|---|---|---|---|
| hereditary class | subgraphs; clouds (zig-zag); links (HDX) | pinned Gibbs measure is Gibbs on the rest (DLR consistency, Markov property) | restricted Gibbs states $\rho_A$; local Lindbladians $\mathcal L_A$ | finite-propagation operators; quotients $\Gamma/\Gamma_n$ |
| restriction operator and its adjoint | cloud average $f\mapsto f^\parallel$; coordinate projection $P_B$ | conditional expectation $\mathbb E[\cdot\mid\sigma_{A^c}]$; heat-bath block update | $\mathbb E_A$ (commuting case); $\mathcal R_{A,t}$ (noncommuting) | spectral functions $f(\Delta)$; Kazhdan projection |
| local certificate | gap of $H$ and of $G$; link gaps; $\Vert P\hat AP\Vert \le\beta+\alpha$ | Dobrushin $\Vert R\Vert _1$; Hayes $\Vert R\Vert _2$; spectral independence $\lambda_1(\Psi^\tau)$; strong spatial mixing | strong clustering; local gap; quantum Dobrushin | Kazhdan constant $K(H,S)$; property (T)/(τ) |
| global conclusion | gap, $O(\log n)$ mixing, EML | $O(n\log n)$ mixing, uniqueness, decay of correlations | size-independent gap; quasi-local recovery; CMI decay | $p\in C^*(X)$; coarse assembly not surjective |
| obstruction | bottleneck (Cheeger); universal cover (Alon–Boppana) | phase coexistence on expanders; tree uniqueness / Kesten–Stigum thresholds | exponential local mixing ($e^{\mu\vert A\vert }$) | ghost classes invisible to local $K$-homology |
| algebra of the proof | $f=f^\parallel+f^\perp$; $2\times2$ angle matrix | path coupling ⇔ Dobrushin (row/column duality) | Prop 20: $\mathrm{Var}_{A\cup B}\le(1-2\epsilon)^{-1}(\mathrm{Var}_A+\mathrm{Var}_B)$ | functional calculus across a spectral gap |

The row "algebra of the proof" is the strongest common thread found: **a family of conditional expectations indexed by a hereditary family of pieces, whose pairwise angle defects (cross-covariance over variance, as in Kastoryano–Brandão's Prop 20) stay uniformly below $\tfrac12$ at a definite scale, has a global gap.** In the framework's vocabulary (Note 1 §6; Note 2 §3), conditional expectations onto face algebras and the Jenčová–Petz generalised conditional expectations are already the central objects.

---

## 11. What this gives the programme

Stated for the abstract objects of Notes 1–2 (layer complexes $K_\ell$ on frames $V_\ell$, the quiver $\Lambda$ of arrows, the cocycle $F$, the KMS/Gibbs states $\mathbb P_\beta$ on histories, conditional expectations onto face and history algebras). Drag test (rule C2): each item would be stated identically for a Bratteli diagram or a substitution tiling; where it would not, it is marked as a prong-2 hypothesis.

### 11.1 Derivation: the history law is a one-dimensional Markov random field, and the sufficiency defect obeys a Poincaré inequality on the re-routing graph

Write a complete history forward as $\mu=(\sigma_0\xrightarrow{\gamma_1}\sigma_1\xrightarrow{\gamma_2}\cdots\xrightarrow{\gamma_n}\sigma_n)$ (Note 1 writes paths right to left; nothing below depends on the convention), $F(\mu)=\sum_iF(\gamma_i)$. Note 1 §5.2: $\mathbb P_\beta(\mu)\propto w_{\sigma_0}e^{-\beta F(\mu)}$.

1. **$\mathbb P_\beta$ is a nearest-neighbour Gibbs measure on the path graph of layers**, with spins = faces (and arrows, if multiplicities are present), pair potential $\beta F$ on consecutive layers, hard constraints given by the quiver, and boundary field $\log w$ at layer $0$. Its Markov property is Hammersley–Clifford in dimension one; this is why Note 1 can say that every Gibbs law on histories is Markov, and why mlp-bridge §3.2 reads measured memory as "outside the Gibbs class".
2. **Squares.** A *square* $\square$ is a pair of length-2 paths with common ends, $\sigma\xrightarrow{\gamma}\tau\xrightarrow{\gamma'}\rho$ and $\sigma\xrightarrow{\eta}\tau'\xrightarrow{\eta'}\rho$ ($\tau=\tau'$ allowed with parallel arrows); $\delta F(\square)=F(\gamma)+F(\gamma')-F(\eta)-F(\eta')$, the alternating sum of U3 in the local-to-global note.
3. **Re-routing graph.** For fixed endpoints $(\sigma_0,\rho)$ let $\mathcal H$ be the set of length-$n$ histories from $\sigma_0$ to $\rho$, and join $\mu\sim\mu'$ if they differ by one **square flip** at some position $i\in\{1..n-1\}$ (only $(\gamma_i,\sigma_i,\gamma_{i+1})$ changes). This is exactly the configuration graph of Glauber dynamics for the chain of item 1 with both ends pinned.
4. **Heat-bath re-routing chain.** $Q$: pick $i$ uniformly, resample $(\gamma_i,\sigma_i,\gamma_{i+1})$ among length-2 paths from $\sigma_{i-1}$ to $\sigma_{i+1}$ with weights $e^{-\beta(F(\gamma_i)+F(\gamma_{i+1}))}$. It is reversible for $\pi=\mathbb P_\beta(\cdot\mid\sigma_0,\rho)\propto e^{-\beta F}$ (the origin weight cancels).
5. **Bound.** For $\mu\sim\mu'$ across a square, $F(\mu)-F(\mu')=\pm\delta F(\square)$. The Poincaré inequality for $Q$ (§9.1) applied to $g=F$ gives
$$
\mathrm{Var}_{\mathbb P_\beta}\big(F(\mu)\,\big|\,\sigma_0,\rho\big)\ \le\ \frac{1}{2\,\mathrm{gap}(Q_{\sigma_0,\rho})}\ \mathbb E_{\pi\otimes Q}\big[\delta F(\square)^2\big]\ \le\ \frac{\max_\square\delta F(\square)^2}{2\,\mathrm{gap}(Q_{\sigma_0,\rho})}.
$$
By Note 1 §6 the left side is exactly the amount by which the pair (input face, output face) fails to be sufficient ("the spread of $F(\mu)$ across histories to a fixed face"). So **the global sufficiency defect is bounded by the local square defects divided by the spectral gap of the re-routing graph.** If the re-routing graph is disconnected the gap is $0$ and no bound follows: a cocycle with vanishing square sums that still varies on a fibre can exist only when that fibre's re-routing graph is disconnected. If $\delta F\equiv0$ and the gap is positive, $F$ is constant on the fibre, which is exactly pair-sufficiency by Note 1 §6's factorisation argument. *Caution for prong 1 (Derivation):* constancy of $F$ on $(s,r)$-fibres is **weaker** than $F$ being a coboundary. At depth 1 with layers $\{a,b\}\to\{c,d\}$ fully joined and no multiplicities, every fibre is a single arrow, so every $F$ makes the pair sufficient, while coboundaries $U(r\gamma)-U(s\gamma)$ satisfy $F(ac)-F(ad)-F(bc)+F(bd)=0$ and span only 3 of the 4 dimensions. Note 1 §6's second bullet ("the pair is sufficient iff $F$ is a coboundary with arbitrary $U$") therefore needs an extra hypothesis; pair-sufficiency is the vanishing of $\delta F$ on all squares that occur inside complete histories, and a coboundary additionally needs $\Phi(\sigma_0,\tau)-\Phi(\sigma_0',\tau)$ to be independent of $\tau$, where $F(\mu)=\Phi(s\mu,r\mu)$. The first bullet (output face alone, histories with full past ending at every layer) is unaffected.
6. **Converse (Cheeger reading).** A re-routing graph with a bottleneck admits cocycles with tiny square defects and a large sufficiency defect (constant on each side of the cut, different across it). The local certificate controls the global defect iff the re-routing graph has no bottleneck: this is the expander moral, transplanted.

**Reading.** This is U3 of the local-to-global note in its first rigorous, degree-1 form, with one correction: what local square sums control is the distance of $F$ from the cochains that are constant on fibres (pair-sufficiency), not from the coboundaries; the controlling constant is the spectral gap of the square-flip walk, a Poincaré (degree-1) shadow of the coboundary-expansion constant U3 asked for. The gap itself can be certified by §9.2 when constraints are soft: with influences $\varrho$ of each neighbouring layer on the resampled block, Dobrushin on the path graph ($2\varrho<1$) gives the path-coupling contraction $1-(1-2\varrho)/(n-1)$ per step (DGJ Lemma 16), and a contraction of this kind bounds the gap below by $(1-2\varrho)/(n-1)$ (**Fact (memory)**: contraction in a metric implies a spectral gap, Chen 1998 / LPW Thm 13.1; Mossel–Sly use the same step via LPW Cor 12.6). With hard constraints (most face pairs not joined by arrows) influences can be $1$, and connectivity of the re-routing graph (total connectivity, CLV Assumption 1.8) is the first thing to check.

### 11.2 Derivation: the angle inequality applies to the framework's conditional expectations

Kastoryano–Brandão's Prop 20 uses "only very general properties of the conditional expectations; in particular we have not assumed that $\rho$ is a Gibbs state, nor that $\mathbb E$ has any local structure". Hence for any faithful state $\omega$ and any two pieces $S,T$ (sub-frames of a layer, or windows of depth in the history algebra) with $\omega$-preserving conditional expectations $\mathbb E_S,\mathbb E_T,\mathbb E_{S\cup T}$ satisfying the relations used in that proof (to be checked against KB §III before use): if $\mathrm{Cov}_{S\cup T}(\mathbb E_Sf,\mathbb E_Tf)\le\epsilon\,\mathrm{Var}_{S\cup T}(f)$ with $\epsilon<\frac12$, then $\mathrm{Var}_{S\cup T}(f)\le(1-2\epsilon)^{-1}(\mathrm{Var}_S(f)+\mathrm{Var}_T(f))$.

**Conjecture E1 (angle certificate ⇒ gap for sub-frame resampling).** If, for a state on the faces of a layer complex, the angle defect $\epsilon(S,T)$ of overlapping sub-frames decays with a distance between $S\setminus T$ and $T\setminus S$ (in the 1-skeleton of $K_\ell$, or in the overlap cosine $G$ of Note 2), then the dynamics that resamples a random sub-frame from its conditional law has a gap independent of the width $|V_\ell|$. *What must be true:* a doubling family of pieces (as in KB Thm 23) inside the complex; on a non-amenable complex the boundary of a piece is proportional to its volume and the doubling argument fails (Mossel–Sly), so the HDX route (links instead of boxes, §9.5) is the one to try. *Guard:* in a layer of fixed width there is no limit; "independent of width" needs a family.

### 11.3 Derivation and Conjecture: spectral independence for laws on faces

**Derivation.** A law on $\{0,1\}^{V}$ supported on the faces of a simplicial complex $K$ is **totally connected** in the sense of CLV Assumption 1.8: pinning a set $S$ in and a set $O$ out leaves $\{\sigma\in K:S\subseteq\sigma,\ \sigma\cap O=\emptyset\}$, which is closed under removing a free vertex (because $K$ is down-closed), so every element is joined to $S$ by single-vertex removals. Single-vertex toggling is the framework's Lüders exclusion/inclusion of a vertex (Note 2 §3).

**Conjecture E2 (U6 made precise).** If the law on faces of a layer is $\eta$-spectrally independent under all pinnings (pinnings = restrictions to sub-frames with prescribed included and excluded vertices) then single-vertex toggling mixes in time polynomial in the width; with $b$-bounded marginals and a bounded-degree dependency structure, in $O(n\log n)$ (CLV Thm 1.12). *What must be true:* (i) a lower bound $b$ on conditional inclusion probabilities, which fails if some vertices are almost never included; (ii) a bounded-degree interaction graph, which is not automatic for an arbitrary law on faces, so only the polynomial (Alev–Lau/ALO) form transfers without it. *Status:* the implication is a theorem once the hypotheses hold; the conjecture is that natural states of the framework satisfy them.

### 11.4 Conjecture (prong-2 hypothesis): a Ramanujan-type certificate for pairwise face laws

If the law on faces of a layer is read as a pairwise exponential family on $\{\pm1\}^{V_\ell}$ with interaction matrix $J$ (a dictionary choice: a maximum-entropy pairwise model of the face law), EKZ Thm 1 gives a Poincaré inequality for single-vertex toggling with constant $1-(\lambda_{\max}(J)-\lambda_{\min}(J))$, and the Hayes form needs only $\lambda_{\max}$ of the influence matrix. *Drag test:* "a state on a frame of commuting projections is a pairwise exponential family" is not a statement about tilings or Bratteli diagrams, so this stays in prong 2.

### 11.5 Speculation: a KMS-detailed-balanced square-flip Lindbladian on the history algebra

Note 1 §5.2 supplies what Chen–Rouzé need and the classical theory does not: a faithful **KMS state** on a noncommutative algebra (the complete-history corner $\bigoplus_{\sigma_0}M_{N_L(\sigma_0)}$) whose modular flow is weighted depth. A KMS-detailed-balanced Lindbladian is defined from a KMS state and a set of jump operators. The natural "single-Pauli jumps on $A$" are the **square flips inside a depth window $A$**, realised as partial isometries between matrix units of histories that differ by one flip (Interpretation). Then:
- the fixed-point algebra of the window generator should be the algebra of observables invariant under re-routing inside $A$, i.e. the observables that see the history only outside $A$ and through the window's end faces;
- the time-averaged map $\mathcal R_{A,t}$ is a gap-free, quasi-local (in depth) recovery map of a state from its restriction outside $A$, by the derivation in §9.9 (only self-adjointness in the KMS inner product is used);
- for diagonal states (measures on histories) this collapses to classical heat-bath re-routing (§11.1); the content is for **states with coherences between histories** (Note 1 §4.1), where a noncommutative approximate Markov property along depth (CMI decaying with window width) is the graded replacement for the Markov property that the realised history law lacks (mlp-bridge §3.2).

*What must be true:* a notion of locality in depth strong enough for a Lieb–Robinson-type bound (finite depth range of the dynamics $\alpha^F$, which holds since $\alpha^F_t(s_\mu s_\nu^*)=e^{it(F(\mu)-F(\nu))}s_\mu s_\nu^*$ multiplies by phases); and a family (stationary quiver, U5) so that "decay with window width" has meaning. Nothing here is proved.

### 11.6 Speculation: the sufficiency projection as a Kazhdan/ghost projection

For a stationary layered quiver (U5), the conditional expectation onto re-routing-invariant functions of histories (the "sufficiency projection") is the spectral projection at $0$ of the re-routing Laplacian of §11.1. If the re-routing graphs form an **expander family** as depth grows, then by the Willett–Yu mechanism (§9.11) this projection is a norm limit of finite-propagation (finite-window) operators and, since fibres grow, a ghost: present globally, invisible to every finite window. Its $K$-theory class would be a candidate for the "gap label of the sufficiency defect" that U5 asks for. *What must be true:* stationarity (a direct-limit structure), connectivity and uniform gap of the re-routing graphs, and a coarse structure on the history space in which flips have bounded propagation. Labelled speculation.

### 11.7 Guards

- **Families, not graphs.** Expansion is a statement about families (Tao Ex 6). A finite network with fixed widths and depth has no expander property; the size parameter must be width (number of units, or of histories) or depth under stationarity.
- **Low temperature.** Expansion of the interaction structure creates bottlenecks at large $\beta$ (§9.7). A gap theorem for re-routing at large $\beta$ in $\mathbb P_\beta\propto e^{-\beta F}$ should not be expected when $F$ has competing near-minimal histories separated by high-action flips.
- **Markov ≠ expansion.** The Markov property (exact, static) does not imply decorrelation, and decorrelation does not imply Markov; do not import "MRF ⇒ fast mixing".
- **Partite structure.** Layered and bipartite graphs have $\lambda_n=-d$; use lazy walks (Tao Ex 27) and one-sided bounds, as the local-to-global note already requires for HDX.
- **Hard constraints and marginals.** Dobrushin, Hayes and spectral-independence theorems need bounded marginals and soft or totally connected constraints; quivers with sparse arrows violate the first by design.
- **Commuting vs noncommuting.** Kastoryano–Brandão's framework is for commuting potentials and fails otherwise (their §IX); Chen–Rouzé need bounded interaction degree and Lieb–Robinson; the framework's history algebra is a direct sum of full matrix algebras, so any "locality" must be imported from the quiver, not from a lattice.

---

## 12. Explicit list of known theorems connecting the areas

| link | theorem | source (as read) | status |
|---|---|---|---|
| MRF/Gibbs ↔ expanders | Dobrushin $\Vert R\Vert _1<1$ ⇒ unique Gibbs measure and $O(n\log n)$ Glauber mixing; Dobrushin–Shlosman $\Vert R\Vert _\infty<1$; any matrix norm $<1$ suffices | DGJ Lemmas 1–3, Cors 18, 25 | Fact |
| MRF/Gibbs ↔ expanders | Hayes: $\Vert R\Vert _2<1$ ⇒ rapid mixing; for Ising $\lambda_{\max}(A_G)\tanh\beta<1$ | DGJ §1; derivation §9.2 | Fact + Derivation |
| MRF/Gibbs ↔ expanders | EKZ: $(1-\Vert J\Vert )\mathrm{Var}\le\mathcal E$; Ramanujan regime $\beta<1/(4\sqrt{d-1})$ for diluted SK on random $d$-regular graphs via Friedman | EKZ Thm 1, §5 | Fact |
| MRF/Gibbs ↔ expanders | Mossel–Sly: $(d-1)\tanh\beta<1$ ⇒ $Cn\log n$ on every graph of max degree $d$; exponential beyond on random regular graphs | Mossel–Sly Thms 1–2 | Fact |
| MRF/Gibbs ↔ HDX | $\eta$-spectral independence ⇒ local spectral expansion $\zeta_k=\eta/(n-k-1)$ ⇒ MLSI, $O(n\log n)$ with bounded degree and marginals | CLV Claim 1.18, Thms 1.12, 1.19 | Fact |
| MRF/Gibbs ↔ trees | BKMP: Glauber gap on $b$-ary tree bounded iff $b\tanh^2\beta<1$; polynomial at all temperatures; cut-width bound $ne^{(4\xi+2\Delta)\beta}$ | BKMP Thm 1.4, Prop 1.1 | Fact |
| Markov semigroups ↔ Gibbs | gap $O(1)$ ⇒ exponential decay of point-to-set correlations on bounded-degree graphs | BKMP Thm 1.5 | Fact |
| Markov semigroups ↔ Gibbs (quantum, commuting) | gap of Davies/heat-bath generator size-independent ⇔ strong clustering; 1D always gapped; high temperature gapped | KB Thms 23, 26, 28, 30–31 | Fact |
| Markov semigroups ↔ Gibbs (quantum, noncommuting) | time-averaged single-Pauli-jump KMS Lindbladian is a quasi-local recovery map; CMI decays exponentially with shielding distance (prefactor $e^{\mu\vert A\vert }$) | Chen–Rouzé Thm III.1, Cors III.1–III.2 | Fact |
| Gibbs ↔ Markov (static) | CMI $\le C_\beta e^{C_\beta g_{A\vert A^c}-c_\beta r}$ by comparison with the cut Gibbs state | Yang Thm II.1 | Fact |
| Markov semigroups ↔ expanders | gap ⇔ Poincaré inequality for the walk's Dirichlet form | Tao Ex 13; EKZ eq. (4) | Fact |
| Markov semigroups ↔ expanders | Cheeger / conductance: $\frac{d-\lambda_2}2\le h\le\sqrt{2d(d-\lambda_2)}$; Jerrum–Sinclair for reversible chains | HLW Thm 4.11; JS (memory for the constants) | Fact / Fact (memory) |
| Markov semigroups ↔ expanders | entropy increase per step ⇔ log-Sobolev constant | HLW §3.1.2 (pointer only) | Fact (pointer) |
| expanders ↔ local-to-global | zig-zag theorem: $\varphi\le\alpha+\beta$, $\varphi\le1-(1-\beta^2)(1-\alpha)/2$; SL = L | HLW Thm 9.1, §9.5 | Fact |
| expanders ↔ local-to-global | Garland, trickling down, Kaufman–Oppenheim, Alev–Lau, ALOV, Dinur–Kaufman | local-to-global note §2 | Fact (there) |
| expanders ↔ random matrices | Wigner, Füredi–Komlós, McKay, Friedman ($2\sqrt{d-1}+\epsilon$), Bordenave | HLW §7; EKZ refs | Fact |
| expanders ↔ spectral extremality | Alon–Boppana, Serre, Greenberg–Lubotzky; LPS/Margulis/Morgenstern; MSS (signings, interlacing, matching polynomial, path tree) | HLW §5–6; MSS | Fact |
| expanders ↔ Gibbs (shared object) | Godsil's path tree (MSS) = Weitz's self-avoiding-walk tree (Mossel–Sly Lemma 13, CLV Thm 6.3) | MSS §2–5; Mossel–Sly §5.1; CLV §6 | Fact (object); Interpretation (meaning) |
| expanders ↔ NCG | property (T)/(τ) ⇒ expanders; $K/(2\vert S\vert )<g<K/2$ | HLW §11 | Fact |
| expanders ↔ NCG | basic Kazhdan projection is a non-compact ghost in $C^*(X)$; coarse Baum–Connes injective, not surjective, for expanders of large girth; Gromov monsters | Willett–Yu I–II; Higson; HLS | Fact |
| expanders ↔ metric geometry ↔ NCG | $c_2(G)=\Omega(\log n)$ for expanders (LLR); Bourgain $O(\log n)$; Yu: coarse embeddability ⇒ coarse BC | HLW Thms 13.1, 13.8; 2511.22438 intro | Fact |
| expanders ↔ quantum channels | random unitary channels: $\vert \lambda_2\vert \to2\sqrt{D-1}/D$; quantum Alon–Boppana | Hastings eqs. (3), (12) | Fact |
| Gibbs ↔ NCG | KMS states of groupoid algebras = quasi-invariant measures with Radon–Nikodym cocycle $e^{-\beta c}$; on the history corner, Gibbs measures with modular flow = weighted depth | Note 1 §5.2 (Neshveyev; Christensen–Thomsen) | Fact (there) |
| Gibbs ↔ NCG | Petz recovery: vanishing CMI ⇔ exact recovery; Fawzi–Renner quantitative | Yang §I | Fact (as cited) |
| Gibbs ↔ NCG | commuting Hamiltonians: Gibbs states exactly Markov (Brown–Poulin) | Yang Table I | Fact (as cited) |
| Gibbs ↔ subfactors | strong clustering as approximate commuting square | §9.8 | Interpretation (commuting squares from memory) |
| framework | sufficiency defect $\le\max\delta F^2/(2\,\mathrm{gap})$ on the re-routing graph | §11.1 | Derivation |

---

## 13. Messages to the other prongs

**To prong 1 (theory).** (o) A correction to check (rule C4): Note 1 §6's equivalence "pair (input face, output face) sufficient iff $F$ is a coboundary with arbitrary $U$" fails at depth 1 (counterexample in §11.1, item 5); the correct local form of pair-sufficiency is vanishing square sums $\delta F$ on the squares occurring in complete histories, which is what §11.1 controls quantitatively. Note 1 §9 checked only the forward direction (coboundary ⇒ pair sufficient), which is true. (i) §11.1 gives U3 a rigorous first form: define the square-flip (re-routing) graph of a layered quiver with fixed endpoints, its $\mathbb P_\beta$-weighting, and its gap; the open question is whether this gap is controlled by link data (HDX-style trickling down along depth) and how it behaves under the stationarity of U5. (ii) Kastoryano–Brandão's angle inequality (§11.2) applies to the framework's conditional expectations; the definition needed is the angle defect of two overlapping faces under a state, and its relation to the overlap cosine $G$ of Note 2. (iii) Note 1's KMS structure is exactly the input of a KMS-detailed-balanced Lindbladian; defining one with square-flip jumps (§11.5) would give a noncommutative, graded version of the Markov property for states with coherences between histories.

**To prong 2 (bridge).** Quantities computable from dictionary v1 data, each targeting a dichotomy: (a) the square defects $\delta F(\square)$ on all squares of the realised quiver, and the spectral gap of the heat-bath re-routing chain on histories with fixed endpoints; then check the Poincaré bound of §11.1 against the measured endpoint $R^2$ and centrality defect (E3). A large defect with small square defects must come with a small re-routing gap (a bottleneck); if not, the estimator or the dictionary is wrong (rule P2.4). (b) Conditional mutual information between future and past *given the whole window of width $w$*, as a function of $w$ (E2 currently conditions only on the present face): exponential decay in $w$ is the approximate Markov property; non-decay means an effective potential of range at least the depth. (c) For the law on faces of a layer: lower bounds on conditional inclusion probabilities and the top eigenvalue of the influence matrix under a sample of pinnings (spectral independence). Negative outcomes are charged first to the dictionary (rule P2.1).

**To prong 3 (unlocks).** U3 as written bounds the distance of $F$ from the coboundaries; by the depth-1 counterexample of §11.1 the quantity that local square sums control is the distance from fibre-constant cochains (pair-sufficiency), and §11.1 gives that control with the re-routing gap as constant. §9.5 promotes U6 from conjecture to "theorem once the hypotheses hold" (CLV Thm 1.12; total connectivity is automatic for down-closed supports, §11.3), leaving the hypotheses as the work. §9.8 and §10 add the conditional-expectation (angle) form of the common mechanism, which is the form closest to the framework's objects. §9.11 adds an NCG row: a spectral gap puts the global projection into the algebra of local operators while keeping it a ghost.

---

## 14. Sources read for this digest

- T. Tao, *254B, Notes 1: Basic theory of expander graphs*, What's new, 2 December 2011, <https://terrytao.wordpress.com/2011/12/02/245b-notes-1-basic-theory-of-expander-graphs/> (full text).
- S. Hoory, N. Linial, A. Wigderson, *Expander graphs and their applications*, Bull. Amer. Math. Soc. 43 (2006) 439–561, <https://www.cs.huji.ac.il/~nati/PAPERS/expander_survey.pdf> (full text; stand-in for the user's local `expander_survey.pdf`).
- E. Mossel, A. Sly, *Exact thresholds for Ising–Gibbs samplers on general graphs*, Ann. Probab. 41 (2013) 294–328, [0903.2906](https://arxiv.org/abs/0903.2906).
- M. Dyer, L. A. Goldberg, M. Jerrum, *Matrix norms and rapid mixing for spin systems*, Ann. Appl. Probab. 19 (2009) 71–107, [math/0702744](https://arxiv.org/abs/math/0702744).
- R. Eldan, F. Koehler, O. Zeitouni, *A spectral condition for spectral gap: fast mixing in high-temperature Ising models*, [2007.08200](https://arxiv.org/abs/2007.08200).
- Z. Chen, K. Liu, E. Vigoda, *Optimal mixing of Glauber dynamics: entropy factorization via high-dimensional expansion*, STOC 2021, [2011.02075](https://arxiv.org/abs/2011.02075).
- N. Berger, C. Kenyon, E. Mossel, Y. Peres, *Glauber dynamics on trees and hyperbolic graphs*, [math/0308284](https://arxiv.org/abs/math/0308284).
- M. J. Kastoryano, F. G. S. L. Brandão, *Quantum Gibbs samplers: the commuting case*, Commun. Math. Phys. (journal details not checked), [1409.3435](https://arxiv.org/abs/1409.3435).
- C.-F. Chen, C. Rouzé, *Quantum Gibbs states are locally Markovian*, [2504.02208](https://arxiv.org/abs/2504.02208).
- T. H. Yang, *Improved estimate of local Markovianity for quantum Gibbs states*, [2609.38007](https://arxiv.org/abs/2609.38007).
- R. Willett, G. Yu, *Higher index theory for certain expanders and Gromov monster groups I, II*, [1012.4150](https://arxiv.org/abs/1012.4150), [1012.4151](https://arxiv.org/abs/1012.4151) (abstracts, §1, Ex 5.3, Cor 6.2).
- M. B. Hastings, *Random unitaries give quantum expanders*, Phys. Rev. A 76 (2007) 032315, [0706.0556](https://arxiv.org/abs/0706.0556).
- A. Marcus, D. Spielman, N. Srivastava, *Interlacing families I: bipartite Ramanujan graphs of all degrees*, Ann. of Math. 182(1) (2015), [1304.4132](https://arxiv.org/abs/1304.4132).
- Seen only as abstracts or introductions: *On the spectral gap of random quantum channels*, [1811.08847](https://arxiv.org/abs/1811.08847); *K-theory of ghostly ideals for ℓ^p-coarsely embeddable spaces*, [2511.22438](https://arxiv.org/abs/2511.22438) (for the statement of Yu's theorem and the role of ghosts).
- Cited through the papers above, not opened: Dobrushin (1968, 1970); Dobrushin–Shlosman; Bubley–Dyer (path coupling); Hayes (FOCS 2006); Wu (Ann. Probab. 2006); Bauerschmidt–Bodineau (JFA 2019); Weitz (STOC 2006); Sly (2010); Sly–Sun (FOCS 2012); Galanis–Štefankovič–Vigoda (CPC 2016); Anari–Liu–Oveis Gharan (FOCS 2020); Martinelli (Saint-Flour 1997); Bakshi–Liu–Moitra–Tang ([2510.08542](https://arxiv.org/abs/2510.08542)); Brown–Poulin ([1206.0755](https://arxiv.org/abs/1206.0755)); Kato–Brandão ([1609.06636](https://arxiv.org/abs/1609.06636)); Higson; Higson–Lafforgue–Skandalis (GAFA 2002); Bordenave ([1502.04482](https://arxiv.org/abs/1502.04482)); Mohanty–O'Donnell–Paredes (STOC 2020).
