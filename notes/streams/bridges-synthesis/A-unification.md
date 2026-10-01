# A — The unifying mathematics

### One stage (a faithful state and its "forget-$A$" subalgebras), two defects (the angle between two of them; the leakage of a Connes cocycle out of one of them), three mechanisms (gap, Cesàro mean, modular resolvent). A dictionary across Markov random fields, expanders, high-dimensional expanders, detailed-balance quantum Markov semigroups and noncommutative geometry, with every identification labelled and the user's intuition tested.

*Bridges-synthesis stream, angle A, 2026-10-01. Prong 3 of the [research programme](../../research-program.md), with messages to prongs 1 and 2 (§11). Inputs, all read in full for this note: the bridge digests [arxiv-2609.38007.md](../../digests/bridges/arxiv-2609.38007.md), [arxiv-2504.02208.md](../../digests/bridges/arxiv-2504.02208.md), [hammersley-clifford.md](../../digests/bridges/hammersley-clifford.md), [expanders.md](../../digests/bridges/expanders.md), [hdx-spectral-independence.md](../../digests/bridges/hdx-spectral-independence.md), [nc-dirichlet-lindblad.md](../../digests/bridges/nc-dirichlet-lindblad.md), [transfer-spectrum-measurement.md](../../digests/bridges/transfer-spectrum-measurement.md), [chat-2609.38007-retrieval-status.md](../../digests/bridges/chat-2609.38007-retrieval-status.md); the programme notes [research-program.md](../../research-program.md), [conditional-arrow-algebra.md](../../conditional-arrow-algebra.md) (Note 1), [simplicial-complex-as-decomposition.md](../../simplicial-complex-as-decomposition.md) (Note 2), [local-to-global-unlocks.md](../../local-to-global-unlocks.md), [mlp-bridge.md](../../mlp-bridge.md); for competition context [competition-plan.md](../../competition-plan.md) §3.1, §6b, §7 and the stream reports under [streams/](../). Nothing here changes the objects of Notes 1–2. Every definition is stated for an inclusion of finite-dimensional algebras with a faithful state, so it reads the same for a spin lattice, a Bratteli diagram or a tiling (drag test C2).*

**Labels.** Every claim carries exactly one label.
- **THEOREM**: a published theorem, cited to the digest section where it was read, or an elementary statement proved in full here.
- **KNOWN-LINK**: a published connection between two areas, cited.
- **DERIVED**: proved here, with the proof included; where marked, also checked numerically (checks C1–C11, §12 and the Appendix).
- **ANALOGY**: a structural similarity, with the exact point where it breaks.
- **CONJECTURE**: a precise statement, with what must be true.
- **SPECULATION**.

Numbers from the checks are evidence attached to a labelled claim, not claims of their own.

**Sigla.** [D-Y] Yang digest (arXiv:2609.38007). [D-CR] Chen–Rouzé digest (arXiv:2504.02208). [D-HC] Hammersley–Clifford digest. [D-EXP] expanders digest (Tao 254B Notes 1 and Hoory–Linial–Wigderson, the latter very probably the user's `expander_survey.pdf`, [D-EXP] §0). [D-HDX] high-dimensional expanders / spectral independence digest. [D-NC] noncommutative Dirichlet forms digest. [D-TS] transfer-spectrum measurement. [D-CHAT] chat retrieval status. Inside those: [CR] Chen–Rouzé, [Y] Yang, [CKG] Chen–Kastoryano–Gilyén, [KB] Kastoryano–Brandão, [BCR] Bardet–Capel–Rouzé, [ALO] Anari–Liu–Oveis Gharan, [CLV] Chen–Liu–Vigoda, [LMRRW] Liu–Mohanty–Raghavendra–Rajaraman–Wu, [VW] Vernooij–Wirth, [JRSWW] Junge–Renner–Sutter–Wilde–Winter, [C26] Chen 2605.02877, [BLMT] Bakshi–Liu–Moitra–Tang, [DLL] Ding–Li–Lin, [HJPW] Hayden–Jozsa–Petz–Winter, [CNNR] Carey–Neshveyev–Nest–Rennie, [Pet] Peterson (the last four through [D-NC] and [D-HC]). Here [KB] is Kastoryano–Brandão; [D-Y] and [D-HC] use the same siglum for Kato–Brandão.

**Retrieval.** No new source was opened for this note; everything cited is cited through the digests, which quote their sources. The ChatGPT conversation behind [Y] was not retrievable ([D-CHAT] §1), so nothing here describes it. Standard facts used *from memory* (rule P3.3) are marked so where they occur: the alternating-projection rate (Aronszajn, Kayalar–Weinert), the identity "Friedrichs cosine $=\lVert PQ-P_{U\cap V}\rVert$" (Deutsch), Takesaki's criterion, the Pimsner–Popa constant and Popa's commuting squares as names, Petz's equality condition, Witsenhausen's maximal-correlation identity, the Diaconis–Saloff-Coste comparison theorem, Chebyshev acceleration. Each is either used only as a name or is also checked numerically.

---

## 0. The answer in brief

1. **There is one stage.** Fix a faithful state $\varphi$ (density $\rho$) on a finite-dimensional algebra $M$, and the family of **forget-$A$ subalgebras** $N_A\subset M$: the observables that do not see the region $A$. All five areas are statements about where the subspaces $L^2(N_A)$ sit inside the state's standard $L^2$-space (the KMS embedding $X\mapsto\rho^{1/4}X\rho^{1/4}$; in the commutative case simply $L^2(\mu)$). §1.
2. **On that stage there are two independent defects, not one.** The user's intuition merges them.
   - **Defect I, the angle** $c(A,B)=\lVert E_AE_B-E_{A\cup B}\rVert$ between two forget-algebras. One functional gives four familiar things, depending on the pair it is evaluated on. It is **zero on separated pairs exactly when the field is Markov**; it **decays across a buffer** under strong spatial mixing; it **equals $\lambda(G)/d$ on the two ends of an edge** (expander mixing); and on adjacent binary spins it is the **largest conditional correlation, at most the geometric mean of the two Dobrushin influences**. Spectral independence is a frame bound for the single-site versions, and the local-to-global step is the sharp two-projection inequality $(1-c)\,\mathrm{Var}_{A\cup B}\le\mathrm{Var}_A+\mathrm{Var}_B$. All DERIVED, §3.
   - **Defect II, the cocycle leakage.** For quantum states the commuting-square form of the Markov property is unavailable. The state-preserving conditional expectations it needs do not exist, by Takesaki, even for commuting Gibbs states embedded in matrix algebras. What survives is: **$I(A{:}C|B)_\rho=0$ iff the Connes cocycle $\rho_{BC}^{it}\rho^{-it}$ stays inside $M_{AB}\otimes1$ for all $t$**. Quantitatively, $I(A{:}C|B)\le(\tfrac1\alpha+3)\,d_A^{2\alpha/(1+\alpha)}q^{2\alpha/(1+\alpha)}$, with $q$ the KMS-strip-weighted leakage and $d_A^2$ the Jones index of the inclusion. DERIVED, §4. Classically defect II vanishes exactly when defect I vanishes on the separated pair (both say "Markov"), although the two numbers differ. No counterpart of defect II appears in the expander literature digested here ([D-EXP]).
3. **The three mechanisms share one proof shape. They are not one theorem.** Each approximates the projection onto an $L^2(N)$ by a function of a positive operator canonically attached to the state, and then needs a light cone.
   - Gap: powers of averaged conditional expectations.
   - Chen–Rouzé: a Cesàro mean of a KMS-symmetric generator, made to work without a gap by **energy duality**: the defect to be recovered is orthogonal to $L^2(N_A)$ and is paired against a vector of small Dirichlet energy.
   - Yang: the resolvent of the relative modular operator, expanded in modular time with the KMS-strip kernel $1/\lvert\sinh\pi t\rvert$.

   The shared shape is DERIVED instance by instance. The claim that it is "one mechanism" is an ANALOGY, and it breaks at the light cone and at the size factor. §5.
4. **The time-averaged detailed-balanced single-Pauli Lindbladian, in this language.** It is the Cesàro mean of the heat flow of the KMS-twisted derivation of the inclusion $N_A=M_{A^c}\subset M$: a gap-free, quasi-local, exactly stationary surrogate for a $\rho$-preserving conditional expectation onto $N_A$ that does not exist. At $\beta=0$ it is literally the Cayley Laplacian of the Pauli group of $A$ with its single-site generators: the twirl is its Kazhdan projection, the Laplacian gap is $4$ for every $\lvert A\rvert$ (tensorization), and the Kazhdan constant is $8/(3\lvert A\rvert)$, so this is not an expander family (DERIVED). **Its genuinely new idea:** *recovery is not mixing*. What the proof needs is coercivity of the dressed single-site derivation **relative to the algebra $N_A$**, not convergence to its exact kernel, and [CR] prove a Hölder form of that relative coercivity unconditionally, at every temperature. §6.
5. **A correction to the "expander half" of [CR] (DERIVED, with checks).** For a non-commuting Hamiltonian whose fixed-point algebra $F_A$ is trivial, the spectral gap of $L_A$ above its exact kernel is exponentially small in the distance from $A$ to the farthest site. This is the quantity in [CR]'s Definition B.1 as the digest records it. In the checks it falls about 16-fold per site: $1.8\cdot10^{-2}$, $1.1\cdot10^{-3}$, $7.2\cdot10^{-5}$. The Poincaré constant **relative to $N_A$**, $\kappa(A)$, stays between $0.23$ and $0.54$ ($n\le5$; it decreases in $n$ by shrinking steps). Recovery obeys $\lVert\mathcal R_{A,t}[\rho_{-A}]-\rho\rVert_1\le\sqrt{0.41/t}\,\chi_{\rm KMS}(\rho_{-A}\Vert\rho)\,\kappa(A)^{-1/2}$. So the expander-type certificate that would turn local into global Markov is a **relative** inequality. In this gap-free route the $\chi^2$ factor still carries an exponential in $|A|$: it grows by a factor of about 2–3.5 per site of $A$ (checks C10–C11), while the starting relative entropy $D(\rho_{-A}\Vert\rho)=O(\beta|A|)$ is linear (DERIVED; check C11). Two routes around it are an **entropic** relative inequality (a relative modified log-Sobolev inequality, which starts from the linear relative entropy) and a relative inequality with an exponential rate in $t$; both must also control the drift of the $A^c$-marginal (C-2 (ii)). CONJECTURES, §5.4.
6. **NCG ties the areas at two layers.** The measure-theoretic layer is modular theory (Petz sufficiency, the Connes cocycle, Takesaki, the KMS embedding, the Jones index). The first-order layer is the Dirichlet form as a squared derivation, property (T) as innerness of derivations, and the Roe-algebra statement that a gap makes a global projection a norm limit of local operators. At the measure-theoretic layer, and for Dirichlet forms as squared derivations, the ties are literal, theorem by theorem. Property (T) and the Roe-algebra statement are theorems about groups and expanders; carried over to the $\beta>0$ Markov problem they are ANALOGIES (§5.7, §7). Spectral triples and Connes/Carlen–Maas metrics are the $L^\infty$ and transport faces of the same derivation. In finite dimension they carry no Markov content, and no Carlen–Maas $W_2$ branch is available for the KMS-only samplers (it is constructed only under GNS symmetry, [D-NC] L8). §7.
7. **Verdict on the intuition** (§8). Several things are literally theorems: the [ALO] identity, the four regimes of $c$, the $\beta=0$ Kazhdan picture, and [CR]'s conditional "local gap ⇒ global Markov" (T13). The hypothesis of the last one, read literally, fails for non-commuting $H$ (D-7). Its relative replacement is where the expander essence would enter, and that version is a CONJECTURE (C-1, C-2), not a theorem. Markov and expansion are **two halves of one local-to-global mechanism**: heredity (exact, combinatorial) plus uniform decorrelation (spectral). The intuition is **false** where it says that the Markov property is an expansion property. Markov holds at every temperature, including where expansion fails; expansion of the interaction graph *destroys* Yang's boundary gain; and the single-Pauli generator is not a quantum expander.

---

## 1. The stage

### 1.1 Standard form, modular objects

- **THEOREM** ([D-NC] §2.1, §2.4). For a faithful density $\rho$ on $M=M_n(\mathbb C)$ the KMS inner product $\langle X,Y\rangle_\rho=\mathrm{Tr}(X^\dagger\rho^{1/2}Y\rho^{1/2})$ is the pull-back of Hilbert–Schmidt under the symmetric embedding $X\mapsto\rho^{1/4}X\rho^{1/4}$. KMS-symmetric quantum Markov semigroups correspond one-to-one to conservative completely Dirichlet forms on this space (Cipriani; Goldstein–Lindsay; finite-dimensional form in Carlen–Maas). In the commutative case $M=L^\infty(\Omega)$, $\rho=\mu$, this is $L^2(\mu)$ and the reversible Markov chains.
- **THEOREM** ([D-Y] §3.4). The modular group is $\sigma_t(X)=\rho^{it}X\rho^{-it}$. The relative modular operator of a second faithful state $\sigma$ is $\Delta(Y)=\sigma Y\rho^{-1}$ on Hilbert–Schmidt space, and the Connes cocycle is $(D\sigma{:}D\rho)_t=\sigma^{it}\rho^{-it}$, with $\Delta^{it}\rho^{1/2}=\sigma^{it}\rho^{-it}\rho^{1/2}$.

### 1.2 The forget-$A$ algebras and three kinds of conditional expectation

Sites $V$; a region $A\subset V$. The forget-$A$ algebra $N_A$:

| area | $M$, state | $N_A$ |
|---|---|---|
| Markov random field / Gibbs | $L^\infty(\Omega^V)$, $\mu$ | functions of $\sigma_{V\setminus A}$ |
| expander graph | $L^\infty(V\times V)$, uniform law of a directed edge $(X_0,X_1)$ | sites $\{0,1\}$: $N_{\{0\}}=$ functions of $X_1$ |
| HDX / spectral independence | $L^\infty$ of the top faces of $X_\mu$, weight $\mu$ | links: conditioning on a face $(U,\tau)$ is compression by $1_{[\sigma_U=\tau]}$ ([D-HDX] §3.2) |
| detailed-balance QMS | $\bigotimes_{x}M_{d_x}$, Gibbs $\rho_\beta$ | $1_A\otimes M_{A^c}$ |
| NCG | von Neumann algebra $M$, faithful normal $\varphi$ | an inclusion $N\subset M$; its relative commutant $N'\cap M$ is the "fibre" |

Three kinds of conditional expectation onto $N$:
- (a) **$\varphi$-preserving.** **THEOREM** (Takesaki, *from memory*; used in [D-CR] §6.1 and [D-NC] §4.4): it exists iff $N$ is $\sigma^\varphi$-invariant. Classically it always exists (the modular group of a commutative algebra is trivial).
- (b) **Trace-preserving.** It always exists. For $N_A$ it is the twirl $E^\tau_{N_A}(X)=d_A^{-2}\sum_{S}SXS^\dagger$ over a unitary error basis of $A$ (Paulis for qubits; [D-CR] §6.2).
- (c) **Generalised (Petz / Accardi–Cecchini).** It always exists and in general is not idempotent ([D-HC] §6.1).

**Lattice law.** $N_A\cap N_B=N_{A\cup B}$. **DERIVED** (proof below).
- Quantum tensor case: the law always holds.
- Classical case: it holds for a strictly positive $\mu$, and in general it is exactly the *intersection axiom* of conditional independence.
- *Proof.*
  - Quantum: $(1_A\otimes M_{A^c})\cap(1_B\otimes M_{B^c})=1_{A\cup B}\otimes M_{(A\cup B)^c}$ by expansion in a product basis ([D-HC] §6.2 quotes the same fact from Lauritzen–Zwiernik).
  - Classical: a function of $(\sigma_B,\sigma_R)$ that equals $\mu$-a.e. a function of $(\sigma_A,\sigma_R)$ ($R=V\setminus(A\cup B)$, $A,B$ disjoint) depends on $\sigma_R$ alone as soon as every pinning of $\sigma_R$ leaves a connected two-block move graph on the support. [D-HC] §4.2 identifies this connectivity with the intersection axiom; positivity gives it. $\square$

### 1.3 The two defects, defined

- **Angle.** Let $E_1,E_2$ be $\varphi$-symmetric conditional expectations, i.e. orthogonal projections in $L^2(\varphi)$, onto $N_1,N_2$. Set
$$c(N_1,N_2):=\lVert E_1E_2-E_{N_1\cap N_2}\rVert_{L^2(\varphi)} .$$
  This is the cosine of the Friedrichs angle between $L^2(N_1)$ and $L^2(N_2)$ (identity *from memory*, Deutsch; consistent with check C2). Equivalently, it is the largest correlation between an element of $L^2(N_1)$ and an element of $L^2(N_2)$, both orthogonal to $L^2(N_1\cap N_2)$. Write $c(A,B):=c(N_A,N_B)$.
- **Cocycle leakage.** For a faithful reference $\sigma$ and a subalgebra $N$ with trace-preserving $E^\tau_N$:
$$\eta_N(t)=\lVert u_t-E^\tau_N(u_t)\rVert_\infty,\qquad u_t=(D\sigma{:}D\rho)_t,\qquad q_N=\tfrac12\int_{\mathbb R}\frac{\eta_N(t)}{\lvert\sinh\pi t\rvert}\,dt .$$
  This is [Y]'s modular leakage ([D-Y] §3.4); the general-$N$ form is in [D-CHAT] §7.1.
- **Modular core.** $N^\sigma$ is the largest $\sigma^\varphi$-invariant subalgebra of $N$ ([D-CR] §6.1). It decides whether an exact recovery can be *local* (§4.1).

### 1.4 Dirichlet forms are squared derivations; kernels are commutants

- **THEOREM** ([D-NC] §3.0 L3; [VW] Thms 2.4–2.5). A KMS-symmetric generator on $M_n$ has $\langle X,-\mathcal L X\rangle_\rho=\sum_j\lVert[V_j,X]\rVert_\rho^2$ with $\{V_j\}=\{V_j^*\}$. On the standard form it is $\lVert\delta X\rVert^2$ for a derivation twisted by $\sigma_{\mp i/4}$.
- For the [CKG]/[CR] sampler the $V$'s are explicit: Gaussian-filtered single Paulis dressed in modular time, $\mathcal E(X)=\sum_a\iint g(t)h(\omega)\lVert[\hat A^a(\omega,t),X]\rVert^2_\rho$, $g(t)=1/(\beta\cosh(2\pi t/\beta))$. **THEOREM** ([D-CR] §3.5, [CR] Lemmas X.1–X.3).
- **THEOREM** ([VW] Thm 2.5 with KMS self-adjointness; in Lindblad form $\ker\mathcal L=\{G,L_j,L_j^\dagger\}'$, [DLL] Lemma 3 via [D-NC] §2.2; [CR] (4.1)). $\ker\mathcal L=\{V_j\}'$. For jumps generating $M_A$ this gives $\ker\mathcal L_A^\dagger\subseteq N_A$.

---

## 2. The dictionary

### 2.1 Objects

The table is a dictionary of objects, not a list of claims; each identification in it is labelled where it is stated or proved (§1, §3–§7). One entry is weaker than the rest: "HC as vanishing of local $H^1$" is SPECULATION ([D-HC] B4 has not proved the identification).

| object | Markov random fields / Gibbs | expanders | HDX / spectral independence | detailed-balance QMS | NCG |
|---|---|---|---|---|---|
| state space $L^2(\varphi)$ | $L^2(\mu)$ | $L^2$(uniform), or $L^2$(edge law) | $L^2(\Pi_d)$ on top faces of $X_\mu$ | KMS space of $\rho_\beta$ | standard form $L^2(M,\varphi)$ |
| forget-$A$ algebra $N_A$ | functions of $\sigma_{A^c}$ | functions of one endpoint | link algebras (pinnings) | $1_A\otimes M_{A^c}$ | inclusion $N\subset M$ |
| restrict–co-restrict operator | heat-bath block update $E_A$ | one walk step $\hat A$ ($=E_{\{0\}}E_{\{1\}}$ read on functions of one endpoint) | down-up walk; Glauber $=\frac1n\sum_vE_v$ | $\mathcal R_{A,t}\circ(\tau_A\otimes\mathrm{Tr}_A)$ ([D-CR] §7.3) | (generalised) conditional expectation |
| Dirichlet form | $\sum_v\mathbb E[\mathrm{Var}_v f]$ | $\sum_{x\sim y}\lvert f(x)-f(y)\rvert^2$ | link Laplacians (Garland) | $\sum_a\iint gh\lVert[\hat A^a(\omega,t),X]\rVert_\rho^2$ | $\lVert\delta X\rVert^2$, twisted derivation |
| Markov property | $c=0$ on separated pairs (D-1); clique potentials (HC) | (two-site system; vacuous) | links split as joins at separators ([D-HDX] §4.3) | cocycle stays in $M_{AB}$ (D-4); shield-commuting $H$ ([D-HC] §7.4) | Petz sufficiency; commuting square (classical) |
| decorrelation certificate | Dobrushin $\lVert R\rVert<1$; strong spatial mixing | $\lambda/d$; Cheeger $h$ | spectral independence $\eta$ (D-3) | strong clustering; local gap; relative $\kappa(A)$ (§5.4); quantum Dobrushin | Kazhdan constant; $c(N_1,N_2)$ |
| local-to-global step | Martinelli template; Dobrushin contraction | EML, AKS walks, zig-zag | Garland, trickle-down, Alev–Lau product formula | [KB] Thm 23; [BCR] approximate tensorization; [CR] Cor B.2 | alternating projections; two-projection inequality (D-2) |
| time average | locally stationary laws ([LMRRW]) | lazy walk | — | Cesàro $\mathcal R_{A,t}$; Davies secular; Gaussian window | mean ergodic projection; modular-time averages ([Y], [JRSWW]) |
| size factor | $q^{\lvert A\rvert}$ configurations | degree $d$ | $\lvert X(k)\rvert$ | Jones index $d_A^2$ (D-5); [CR]'s $2^{2\lvert A\rvert}$, equal in value for qubits but from word costs ([D-Y] §8.8); cut weight $g_{A\vert A^c}$ ([Y]) | Jones / Pimsner–Popa index |
| obstruction | bottleneck; phase coexistence | Cheeger bottleneck; Alon–Boppana | coboundaries lifted from lower dimension | trivial $F_A$ (modular core collapses); $e^{\mu\lvert A\rvert}$ | non-modular-invariant $N$; ghost classes |
| metric / transport | Hamming $W_1$, Ornstein $\bar d$ | graph metric; distortion $\Omega(\log n)$ | — | De Palma–Marvian–Trevisan–Lloyd $W_1$; Carlen–Maas $W_2$ (GNS only) | Connes distance, Rieffel Lip-norm |
| degree-1 (cohomological) layer | HC as vanishing of local $H^1$ ([D-HC] §3.7, B4) | coboundary expansion | Garland vanishing; agreement expansion | triangles: incompatible Koashi–Imoto splittings ([D-HC] §7.4) | Connes cocycle class; sufficiency iff coboundary (Note 1 §6) |

### 2.2 Theorems that realise the rows

| # | statement | areas | label | where |
|---|---|---|---|---|
| T1 | positive law: Markov ⇔ Gibbs with clique potentials; the proof is a Möbius inversion over commuting idempotents | MRF | THEOREM | [D-HC] §2–3, §3.8 |
| T2 | spectral independence of $\mu$ under all pinnings = local spectral expansion of $X_\mu$, with $\lambda_2(P_\emptyset)=\lambda_{\max}(\Psi_\mu)/(n-1)$ | MRF ↔ HDX | THEOREM | [D-HDX] §3.4 ([ALO] Thm 3.1) |
| T3 | local spectral expansion ⇒ gap of the down-up walk (product formula); with marginal bounds and bounded degree ⇒ MLSI and $O(n\log n)$ | HDX ↔ MS | THEOREM | [D-HDX] §2.6, §3.6 |
| T4 | $(1-\lVert J\rVert)\mathrm{Var}\le\mathcal E$; the Ramanujan regime $\beta<1/(4\sqrt{d-1})$ for diluted SK | MRF ↔ expanders | THEOREM | [D-EXP] §9.3 (Eldan–Koehler–Zeitouni) |
| T5 | approximate commuting square ⇒ $\mathrm{Var}_{A\cup B}\le(1-2\epsilon)^{-1}(\mathrm{Var}_A+\mathrm{Var}_B)$; strong clustering ⇒ size-independent gap, and conversely gaps on every subregion with a projective expectation (Davies) ⇒ strong clustering (commuting $H$) | QMS | THEOREM | [D-EXP] §9.8 ([KB] Prop 20, Thms 23, 26) |
| T6 | approximate tensorization with constant $(1-c_1)^{-1}$, $c_1$ an $L^1\to L^\infty$ clustering constant, plus an additive term for non-commuting algebras | QMS ↔ NCG | THEOREM | [D-HDX] §4.5 ([BCR]) |
| T7 | time-averaged single-Pauli KMS Lindbladian recovers $\rho$ from $\rho_{-A}$ with error $re^{\mu\lvert A\rvert}t^{-\lambda}$, quasi-locally; hence local Markov | QMS | THEOREM | [D-CR] §2 ([CR] Thm III.1, Cors III.1–2) |
| T8 | $I(A{:}C\vert B)\le C_\beta e^{C_\beta g_{A\vert A^c}-c_\beta r}$ by modular-cocycle leakage | QMS ↔ NCG | THEOREM | [D-Y] §4–5 ([Y] Thm II.1) |
| T9 | CMI $=\delta_{AB}(\rho,\rho_{-A})$, the Petz sufficiency defect of $M_{AB}\otimes1$ for $\{\rho,\rho_{-A}\}$ | QMS ↔ NCG | THEOREM (elementary, proved there) | [D-CR] §6.3 |
| T10 | KMS-symmetric QMS ↔ Dirichlet forms ↔ twisted derivations | QMS ↔ NCG | THEOREM | [D-NC] §3 L2–L3 |
| T11 | (T) ⇔ closable derivations inner (finite factors); subexponential spectral growth ⇒ amenable | expanders ↔ NCG | THEOREM | [D-NC] §5.4 |
| T12 | expander: the Kazhdan projection is a norm limit of finite-propagation operators and is a ghost; for expanders of girth $\to\infty$ coarse Baum–Connes is not surjective | expanders ↔ NCG | THEOREM | [D-EXP] §9.11 (Willett–Yu) |
| T13 | uniform local gap ⇒ global Markov | QMS | THEOREM (conditional; stated "by the same reasoning") | [D-CR] §2, App. B (Cor B.2) |
| T14 | strong (post-selected) Markov ⇔ clustering, given approximate detailed balance | QMS | THEOREM | [D-HC] §9.6 ([C26]) |
| T15 | universal recovery = Petz map averaged over modular time, density $\frac\pi2(\cosh\pi t+1)^{-1}$ | NCG | THEOREM | [D-HC] §8.2 ([JRSWW] Thm 2.1) |
| D-1 to D-8 | the propositions proved below, with their corollaries | all | DERIVED | §3–6 |

---

## 3. Defect I: the angle

### 3.1 One functional, four regimes

**D-1 (DERIVED; checks C1).** Let $\mu$ be a strictly positive law on $\Omega^V$, $E_A$ the conditional expectation that integrates out $A$, and $c(A,B)=\lVert E_AE_B-E_{A\cup B}\rVert_{L^2(\mu)}$.

- (a) **Markov.** For disjoint $A,B$: $c(A,B)=0$ iff $\sigma_A\perp\sigma_B\mid\sigma_R$, with $R=V\setminus(A\cup B)$. For non-adjacent pairs on a graph, $c=0$ for all such pairs is the pairwise Markov property, hence (T1) the Gibbs property.
- (b) **Decorrelation across a buffer.** If $A\cup B=V$, then $E_{A\cup B}$ is the mean and $c(A,B)$ is the maximal correlation between $\sigma_{V\setminus B}$ and $\sigma_{V\setminus A}$, two blocks separated by the buffer $A\cap B$, which is integrated out.
- (c) **Expander mixing.** For the uniform law of a directed edge $(X_0,X_1)$ of a $d$-regular graph $G$: $c(\{0\},\{1\})=\lambda(G)/d$ with $\lambda(G)=\max(\lvert\lambda_2\rvert,\lvert\lambda_n\rvert)$. The expander mixing lemma is this inequality evaluated on indicators. (The edge law is not strictly positive. The formula holds for every $d$-regular $G$; it is also the Friedrichs cosine of §1.3 exactly when the lattice law $N_{\{0\}}\cap N_{\{1\}}=\mathbb C$ holds, i.e. when $G$ is connected and non-bipartite. For bipartite $G$ the left side is $1=\lambda(G)/d$, while the cosine relative to the two-dimensional intersection is $\lvert\lambda_2\rvert/d$; on the 8-cycle, $1$ against $0.707$.)
- (d) **Dobrushin.** For adjacent binary spins $i,j$: $c(\{i\},\{j\})=\max_r\lvert\mathrm{corr}(\sigma_i,\sigma_j\mid\sigma_{\rm rest}=r)\rvert=\max_r\sqrt{\Psi_r(i,j)\Psi_r(j,i)}\le\sqrt{R_{ij}R_{ji}}$, the geometric mean of the two Dobrushin influences. For Ising this is $\le\tanh\lvert J_{ij}\rvert$.

*Proof.*
- (a) $N_A=L^2(\sigma_{B\cup R})$, $N_B=L^2(\sigma_{A\cup R})$, and by the lattice law (§1.2) $N_A\cap N_B=L^2(\sigma_R)$.
  - Two $\sigma$-algebras form a commuting square over their intersection iff $\mathbb E[g\mid\sigma_{B\cup R}]=\mathbb E[g\mid\sigma_R]$ for all $g\in L^2(\sigma_{A\cup R})$, i.e. iff $\sigma_A\perp\sigma_B\mid\sigma_R$ ([D-HDX] §4.5: the derivation there, and that digest's own check C3).
  - $c=0$ is the commuting square.
- (b) The two ranges are $L^2(\sigma_{V\setminus A})$ and $L^2(\sigma_{V\setminus B})$, the intersection is the constants, and the norm of $E_AE_B-E_V$ ($E_V$ the mean) is the largest correlation between the two ranges.
- (c) $E_{\{1\}}h=\mathbb E[h\mid X_0]=\tilde h(X_0)$ and $E_{\{0\}}\tilde h(X_0)=\mathbb E[\tilde h(X_0)\mid X_1]=(\hat A\tilde h)(X_1)$, with $\hat A$ the normalised adjacency matrix. So $E_{\{0\}}E_{\{1\}}-E_{\{0,1\}}$ ($E_{\{0,1\}}$ the mean) has norm $\lVert\hat A-J\rVert_{L^2(\mathrm{unif})}=\lambda(G)/d$ (Witsenhausen, *from memory*; checked). Taking $f=1_S$, $g=1_T$ in the correlation bound gives $\lvert\,\lvert E(S,T)\rvert-d\lvert S\rvert\lvert T\rvert/n\rvert\le\lambda\sqrt{\lvert S\rvert\lvert T\rvert}$, as in [D-EXP] §3.
- (d)
  - The operator is block-diagonal over the value $r$ of the other spins.
  - For a binary pair the maximal correlation is $\lvert\mathrm{corr}\rvert$.
  - $\mathrm{corr}^2=\frac{\mathrm{Cov}}{\mathrm{Var}_i}\frac{\mathrm{Cov}}{\mathrm{Var}_j}=\Psi_r(i,j)\Psi_r(j,i)$ with [ALO]'s $\Psi(i,j)=P(j\mid i)-P(j\mid\bar i)$ ([D-HDX] §3.3).
  - The Dobrushin influence of $j$ on $i$ is $R_{ji}=\max_r\lvert\Psi_r(j,i)\rvert$ ([D-EXP] §9.2).
  - The Ising bound $R\le\tanh\lvert J\rvert$ is [D-EXP] §9.2's derivation. $\square$

*Check C1* (Ising chain, $n=9$, $\beta=0.6$, $h=0.25$):
- separated pairs: $c(\{1\},\{3\})=2.8\cdot10^{-16}$ and $c(\{0,1\},\{4,5\})=1.9\cdot10^{-16}$;
- adjacent pair: $c(\{1\},\{2\})=0.465$;
- across a buffer of width $w=1,\dots,5$: $c=0.257$, $0.130$, $0.065$, $0.033$, $0.017$. The ratio per site, $0.504$, matches the transfer-matrix ratio $\lambda_2/\lambda_1=0.5004$;
- random 4-regular multigraph on 60 vertices: $c=0.860929=\lambda/d$;
- adjacent random-Ising pairs: $c=\max_r\lvert\mathrm{corr}\rvert$ exactly, and $0.186,0.901,0.354\le\tanh\lvert J\rvert=0.243,0.922,0.388$.

**Reading.** The Markov property and expansion are **the zero set and the uniform smallness of one functional, evaluated on different pairs**. Markov concerns separated pairs, where the buffer is pinned. Decorrelation concerns overlapping pairs, where the buffer is integrated out, and in the expander case the pair of endpoints. This is the precise sense in which they share an essence. The independence of the two halves is equally visible (check C1 at $\beta=2$, $h=0$): $c(\{1\},\{3\})=2\cdot10^{-16}$ (exact Markov), while the buffer angles are $0.93$, $0.86$, $0.80$ at $w=1,3,5$ (weak decorrelation).

### 3.2 The local-to-global step is the sharp two-projection inequality

**D-2 (DERIVED; check C2).** Let $P,Q$ be orthogonal projections onto closed subspaces $U,V$ of a Hilbert space, $R$ the projection onto $W=U\cap V$, and $c$ the Friedrichs cosine. For every $g\perp W$,
$$(1-c)\lVert g\rVert^2\le\lVert g-Pg\rVert^2+\lVert g-Qg\rVert^2 ,$$
and the constant is attained whenever $U\ne V$: $\lambda_{\min}\big((2-P-Q)\vert_{W^\perp}\big)=1-c$.

*Proof.*
- Since $R\le P$, $RPg=Rg=0$, so $Pg\in U\ominus W$; likewise $Qg\in V\ominus W$. Hence $\lvert\langle Pg,Qg\rangle\rvert\le c\,a\,b$ with $a=\lVert Pg\rVert$, $b=\lVert Qg\rVert$.
- Then $a^2+b^2=\langle(P+Q)g,g\rangle\le\lVert(P+Q)g\rVert\lVert g\rVert$ and $\lVert(P+Q)g\rVert^2\le a^2+b^2+2cab\le(1+c)(a^2+b^2)$. So $a^2+b^2\le(1+c)\lVert g\rVert^2$, and $\lVert g-Pg\rVert^2+\lVert g-Qg\rVert^2=2\lVert g\rVert^2-a^2-b^2\ge(1-c)\lVert g\rVert^2$.
- Sharpness: take unit principal vectors $u\in U\ominus W$, $v\in V\ominus W$ with $\langle u,v\rangle=c$, $Pv=cu$, $Qu=cv$ (they exist in finite dimension when $c>0$; if $c=0$, a unit vector in whichever of $U\ominus W$, $V\ominus W$ is non-zero already gives equality), and set $g=u+v$. Then $\lVert g\rVert^2=2+2c$ and the right side is $2(1-c^2)=(1-c)\lVert g\rVert^2$. $\square$

*Check C2:* over 300 random pairs, $\min[\lambda_{\min}(2-P-Q\vert_{W^\perp})-(1-c)]=-1.3\cdot10^{-15}$; equality, as claimed.

**Corollary (DERIVED).** With $\mathrm{Var}_A(f):=\lVert f-E_Af\rVert^2$ and the lattice law,
$$\mathrm{Var}_{A\cup B}(f)\le\frac{\mathrm{Var}_A(f)+\mathrm{Var}_B(f)}{1-c(A,B)}.$$
This has the shape of [KB] Prop 20 (T5), which states $(1-2\epsilon)^{-1}$ under a covariance hypothesis, and of [BCR]'s approximate tensorization (T6) in variance form. *Proof:* apply D-2 to $g=f-E_{A\cup B}f$. $\square$

**THEOREM** (Kayalar–Weinert, *from memory*; check C2 reproduces it to $2\cdot10^{-15}$): $\lVert(PQ)^k-R\rVert=c^{2k-1}$. Alternating the two restriction–co-restriction operators reaches the intersection at a rate set by the angle alone.

### 3.3 Spectral independence is a frame bound

**D-3 (DERIVED; check C3).** For a law $\mu$ on $\{0,1\}^n$ let $P_i$ be the projection of $L^2(\mu)$ onto $\mathbb C(\sigma_i-\mathbb E\sigma_i)$. Then
$$\lambda_{\max}(\Psi_\mu)=\Big\lVert\sum_iP_i\Big\rVert-1 .$$
So $\eta$-spectral independence says that the $n$ single-site fluctuation lines are almost orthogonal in the frame sense, $\sum_iP_i\le1+\eta$. Under every pinning it is the hypothesis of T2–T3.

*Proof.*
- For binary spins $\Psi(i,j)=\mathrm{Cov}_{ij}/\mathrm{Var}_i$ ($i\ne j$), so $\Psi=\mathrm{diag}(\mathrm{Var})^{-1}\mathrm{Cov}-I$, which is similar to $\mathrm{Cor}-I$.
- With $v_i$ the normalised centred spins, $\sum_iP_i=\sum_iv_iv_i^*$ has the same non-zero spectrum as the Gram matrix $(\langle v_i,v_j\rangle)=\mathrm{Cor}$.
- Hence $\lambda_{\max}(\Psi)=\lambda_{\max}(\mathrm{Cor})-1=\lVert\sum P_i\rVert-1$ (cf. Chen–Eldan Fact 23 and the normalisation guard in [D-HDX] §3.3). $\square$

*Check C3:* four random laws on $\{0,1\}^5$, agreement to $10^{-12}$.

### 3.4 "Hereditary conditional structure + uniform decorrelation ⇒ global gap"

**THEOREM instances.**
- [ALO]/Alev–Lau (T2–T3): the frame bound of the single-site fluctuation lines (D-3), under every pinning ⇒ gap of the average of the single-site expectations, the Glauber walk.
- [KB] Thm 23 (T5): angle bounds for overlapping boxes, decaying with the overlap, at every scale ⇒ size-independent gap.
- [BCR] (T6): an entropy analogue of D-2, with an $L^1\to L^\infty$ clustering constant in place of the $L^2$ cosine and an additive term for non-commuting algebras ([D-HDX] §4.5).

**ANALOGY** (the zig-zag theorem and Dobrushin's condition). The zig-zag bound is a $2\times2$ angle computation between cloud averages and the inter-cloud permutation ([D-EXP] §6). Dobrushin is an $\ell^\infty$-oscillation contraction ([D-EXP] §9.2). Both have the shape "pairwise decorrelation, iterated". *Break:* neither is an $L^2$ statement about conditional expectations onto a lattice of subalgebras. Zig-zag works with an unrelated permutation, and Dobrushin works in the oscillation seminorm, not in $L^2$.

**Where heredity enters, and where it breaks quantumly.**
- *Classical:* pinning a Markov random field gives a Markov random field on the induced subgraph ([D-HC] §5.4, THEOREM, elementary). So the links of $X_\mu$ are of the same kind, and one local argument certifies every link ([D-HDX] §4.2 table).
- The lattice law that $E_{A\cup B}$ needs is the intersection axiom, which positivity implies. In finite dimension $c(A,B)=\lVert E_AE_B-E_{A\cup B}\rVert<1$ iff the lattice law holds (the Friedrichs cosine relative to the true intersection is always $<1$). So **Hammersley–Clifford's positivity secures the qualitative statement "the relevant angle is $<1$", and the expander gap is its quantitative version.** This is [D-HC] B3 (THEOREM, elementary, proved there), restated on this stage.
- *Quantum:* there is no pinning. Truncation is not conditioning, and marginals of a non-commuting Gibbs state are not local Gibbs states ([D-Y] B4); the closest operational analogue of pinning is post-selection on a local measurement ([D-HC] B5). The operational substitute, recovery of every post-selected branch (the strong Markov property), is *equivalent to clustering* (T14). **So quantumly the heredity half itself costs a decorrelation hypothesis: the two halves are no longer independent.** THEOREM ([C26] via [D-HC] §9.6); the reading is [D-HC] B5.

---

## 4. Defect II: the cocycle

### 4.1 Takesaki's obstruction removes the commuting-square definition

**THEOREM** (Takesaki's criterion, *from memory*), with an elementary computation. In the classical Ising ring of [D-CR] §6.1 ($H=\sum_iZ_iZ_{i+1}+0.3\sum_iZ_i$ on 5 sites, $A=\{0\}$) embedded in $M_{2^5}$, $\sigma_t(1_A\otimes X_1)=X_1\,e^{2i\beta t(Z_0Z_1+Z_1Z_2+0.3Z_1)}$ depends on $Z_0$, so the forget-site algebra $N_A$ (dimension 256) is not modular-invariant; [D-CR] §6.1 finds numerically that its modular core $N_A^\sigma$ has dimension 96. So **even a commuting Gibbs state has no $\rho$-preserving conditional expectation onto $M_{A^c}$**. The $\rho$-preserving conditional expectations that §3 used (classically, the $E_A$) are therefore unavailable in matrix algebras. KMS-orthogonal projections onto $L^2(N_A)$ still exist, but they are not conditional expectations, and the commuting-square form of the Markov property has no direct quantum meaning (C-4 in §10 asks what survives).

What survives is the sufficiency form.

### 4.2 Quantum Markov is cocycle localisation

**D-4 (DERIVED).** For faithful $\rho_{ABC}$:
$$I(A{:}C|B)_\rho=0\iff\rho_{BC}^{it}\,\rho^{-it}\in M_{AB}\otimes1_C\quad\text{for all }t\in\mathbb R$$
(here $\rho_{BC}$ stands for $1_A\otimes\rho_{BC}$).

*Proof.*
- By T9, $I(A{:}C|B)=D(\rho\Vert\rho_{-A})-D(\rho_{AB}\Vert(\rho_{-A})_{AB})$, the loss of relative entropy of the pair $\{\rho,\rho_{-A}\}$ under restriction to $N=M_{AB}\otimes1$.
- By the Jenčová–Petz theorem (Note 1 §6, Fact; Petz's equality condition), with $\rho$ faithful the loss vanishes iff $N$ is sufficient for the pair, iff $(D\rho_{-A}{:}D\rho)_t\in N$ for all $t$.
- Finally $(D\rho_{-A}{:}D\rho)_t=(\tau_A\otimes\rho_{BC})^{it}\rho^{-it}=d_A^{-it}(1_A\otimes\rho_{BC}^{it})\rho^{-it}$, and the scalar phase does not affect membership. $\square$

*Classical reading.* With $\rho=\mathrm{diag}\,p$, the cocycle is $p(a\mid b,c)^{-it}$. It lies in $L^\infty(a,b)$ iff $p(a\mid bc)$ does not depend on $c$, which is the Markov property. *Check C4:*
- a Koashi–Imoto Markov state $\rho_{Ab_L}\otimes\rho_{b_RC}$ has leakage $q=1.1\cdot10^{-14}$;
- a commuting Ising chain has $q=3\cdot10^{-16}$;
- a random state has $q=0.68$.

**Why the quantum converse of Hammersley–Clifford fails, in one line** (DERIVED from D-4 and [D-HC] §7.5). For $H=H_{AB}+H_{BC}$ the cocycle is generated by $\log\rho-\log\rho_{BC}$.
- If the terms commute, $\log\rho_{BC}=-\beta H_{BC}+(\text{function on }B)$, and the cocycle stays in $M_{AB}$.
- If they do not commute, $\log\rho_{BC}$ is no longer $-\beta H_{BC}$ plus an operator on $B$: it is the effective Hamiltonian, or Hamiltonian of mean force, of [D-HC] §7.5, whose quasi-locality the sources establish only under extra hypotheses ([D-HC] §9.2: high temperature, with a caveat on the expansion). Whenever the CMI is non-zero (e.g. the Heisenberg chain of [D-HC] §7.5) the cocycle leaves $M_{AB}$ by D-4, and the Markov property is at best approximate.
- Hammersley–Clifford's Möbius inversion acts on the *joint* $\log\rho$, through the commuting reference expectations. The Markov property is about the *marginal* $\log\rho_{BC}$. Classically the two coincide; quantumly only the first is controlled by locality of $H$.

### 4.3 Quantitative form; the Jones index appears

**D-5 (DERIVED; check C4).** For faithful $\rho_{ABC}$ and $0<\alpha\le1$,
$$I(A{:}C|B)_\rho\le\Big(\frac1\alpha+3\Big)\,d_A^{\frac{2\alpha}{1+\alpha}}\;q^{\frac{2\alpha}{1+\alpha}},\qquad q=\tfrac12\int\frac{\big\lVert u_t-E^\tau_{AB}(u_t)\big\rVert_\infty}{\lvert\sinh\pi t\rvert}dt,\quad u_t=\rho_{BC}^{it}\rho^{-it}.$$
In particular $I\le4\,d_A\,q$ (take $\alpha=1$).

*Proof.*
- [Y] Cor III.5 (THEOREM, [D-Y] §5.5) holds for any faithful pair $(\rho,\sigma)$ and $R=AB$: $\delta_{AB}(\rho,\sigma)\le(\frac1\alpha+3)[\mathrm{Tr}\rho^{1+\alpha}\sigma^{-\alpha}]^{1/(1+\alpha)}q_{AB}^{2\alpha/(1+\alpha)}$. Take $\sigma=\rho_{-A}$ and use T9.
- **Pimsner–Popa bound** ([D-CR] §6.2): $\rho_{-A}=d_A^{-2}\sum_SS\rho S^\dagger\ge d_A^{-2}\rho$ (the term $S=1$), so $\rho\le d_A^2\rho_{-A}$.
- Since $x\mapsto x^{-\alpha}$ is operator monotone decreasing for $\alpha\in(0,1]$, $\rho_{-A}^{-\alpha}\le d_A^{2\alpha}\rho^{-\alpha}$, and so $\mathrm{Tr}\rho^{1+\alpha}\rho_{-A}^{-\alpha}=\mathrm{Tr}\,\rho^{\frac{1+\alpha}2}\rho_{-A}^{-\alpha}\rho^{\frac{1+\alpha}2}\le d_A^{2\alpha}$. $\square$

*Check C4* (all $M_\alpha\le d_A^{2\alpha}$; the bound at $\alpha=1$):

| state | CMI | $q$ | bound |
|---|---|---|---|
| random 3-qubit | 0.342 | 0.683 | 3.59 |
| transverse-field Ising, 4 qubits, $\beta=1$ | $1.79\cdot10^{-4}$ | $9.9\cdot10^{-3}$ | $5.4\cdot10^{-2}$ |
| same, $\beta=2$ | $5.1\cdot10^{-3}$ | $7.7\cdot10^{-2}$ | 0.54 |

In the two Gibbs cases CMI$/q^2=1.84$ and $0.85$. That is consistent with [Y]'s middle-range estimate $j(s)\le(1+s)q^2/s^2$ being the dominant term, so the bound's linear power of $q$ is loose there.

**Reading** (DERIVED comparison). D-5 is the "erase-$A$" version of [Y]'s theorem.
- The size factor is $\sqrt{\text{Jones index}}=d_A$. [D-CR] §6.2 reads [CR]'s $2^{2\lvert A\rvert}$ as this index. [D-Y] §8.8, which re-read [CR]'s proof, traces it instead to per-string Gibbs-conjugation and Leibniz costs, the twirl being a normalised average. The two agree in value for qubits ($4^{\lvert A\rvert}=d_A^2$); only D-5's factor is literally an index.
- The cocycle $\rho_{BC}^{it}\rho^{-it}$ involves the modular flow of a *marginal*, $\log\rho_{BC}$, which Lieb–Robinson bounds for $H$ do not control (§4.2).
- [Y] keeps the same inequality and changes the reference to the cut state $\gamma_A\otimes\gamma_{BC}$ (T8). The size factor becomes $\mathrm{Tr}\rho^{1+\alpha}\sigma^{-\alpha}\le e^{3\alpha\beta g_{A|A^c}}$ ([D-Y] §5.6, Lemma IV.2). The cocycle becomes the interaction-picture unitary generated by the cut $V$, which Lieb–Robinson controls (Lemma IV.1).
- **The trade is index for interface.** The exponents of the two prefactors are in the ratio $3\beta g_{A|A^c}/(2\lvert A\rvert\log2)$ (qubits), a Cheeger ratio up to a constant. It tends to $0$ along Følner sequences (large boxes in $\mathbb Z^D$) and is bounded below on expanders (THEOREM, elementary, [D-Y] B1).

### 4.4 The two defects side by side

| | defect I: angle $c(N_A,N_B)$ | defect II: cocycle leakage $q$ |
|---|---|---|
| concerns | two subalgebras, relative position | one subalgebra and a pair of states |
| zero means | commuting square: Markov on separated pairs; independence across a buffer | sufficiency: exact Petz recovery; CMI $=0$ |
| quantitative tool | two-projection inequality (D-2); alternating projections | [Y] Lemma III.4 / Cor III.5 (D-5) |
| classical case | carries both Markov and decorrelation | vanishes exactly when defect I vanishes on separated pairs |
| quantum case | needs $\rho$-preserving conditional expectations, which do not exist (Takesaki); only KMS-orthogonal projections survive | survives; this is approximate Markov |
| expander theory | yes: EML, Kazhdan, spectral independence | no counterpart |

The rows restate D-1, D-2, D-4 and D-5 (DERIVED) and §4.1. The "no counterpart" entry is a reading of [D-EXP], not a theorem.

---

## 5. Three mechanisms, one proof shape

### 5.1 The shape

| mechanism | positive operator attached to the state | function approximating a projection | what is measured | light cone | size factor |
|---|---|---|---|---|---|
| gap (H+D⇒G) | $I-$ average of conditional expectations; down-up walk | $x^k$ or $e^{-tx}$ across a gap | angles / frame bound on every link | walk steps | none when heredity is exact |
| [CR] | $-\mathcal L_A$, the KMS Dirichlet generator of the dressed single-site derivation | Cesàro $\phi_t(x)=\frac{1-e^{-tx}}{tx}$, no gap | Dirichlet energy $\le c_*/t$, paired against a defect orthogonal to $L^2(N_A)$ (D-6) | Lieb–Robinson for $e^{s\mathcal L_A}$, error linear in $t$ | word costs $2^{2\lvert A\rvert}$ ([D-Y] §8.8; equal in value to the index $d_A^2$ for qubits) |
| [Y] | relative modular operator $\Delta_{\sigma,\rho}$, $K=\log\Delta$ | resolvent $(\Delta+s)^{-1}$, written as a $1/\lvert\sinh\pi t\rvert$-weighted average of $\Delta^{it}$ | squared distance of the resolvent vector from $L^2(N_R)\rho^{1/2}$ (Carlen–Vershynina), bounded by cocycle leakage | Lieb–Robinson in modular time against $e^{-\pi\lvert t\rvert}$ | Rényi moment $e^{O(\beta g)}$ |

Each row is a THEOREM in its source (T3/T5, T7, T8). **The claim that the three rows are "one mechanism" is an ANALOGY.** All three estimate a distance to $L^2(N)$ in the standard form of the state, by functional calculus of a canonical positive operator plus a light cone. *It breaks* in two places:
- The light cones are of three different kinds: combinatorial steps, Lindblad time and modular time.
- The size factors are of three different kinds: none, a word cost on the region (equal in value to an index) and an interface. No single inequality specialises to all three.

### 5.2 The Chen–Rouzé row made exact: energy duality

**D-6 (DERIVED; checks C5, C9).** Let $X\ge0$ be self-adjoint on a finite-dimensional Hilbert space, $P_0$ the projection onto $\ker X$, $\phi_t(x)=(1-e^{-tx})/(tx)$ (with $\phi_t(0)=1$), and $c_*=\sup_{u>0}(1-e^{-u})^2/u=0.4073$ (the Cesàro constant; not the angle $c$ of §3).
- (i) $\langle\phi_t(X)O,X\phi_t(X)O\rangle\le\frac{c_*}t\lVert O\rVert^2$.
- (ii) For $v\perp\ker X$: $\lvert\langle\phi_t(X)O,v\rangle\rvert\le\min\big(\sqrt{c_*/t}\,\lVert X^{-1/2}v\rVert,\ t^{-1}\lVert X^{-1}v\rVert\big)\lVert O\rVert$.

*Proof.*
- (i) Spectral theorem: $x\phi_t(x)^2=t^{-1}(1-e^{-u})^2/u$ with $u=tx$.
- (ii) $\langle\phi_t(X)O,v\rangle=\langle X^{1/2}\phi_t(X)O,X^{-1/2}v\rangle$, then Cauchy–Schwarz and (i). The second bound uses $\langle X\phi_t(X)O,X^{-1}v\rangle$ and $\sup_{x>0}x\phi_t(x)=t^{-1}\sup_u(1-e^{-u})=t^{-1}$. $\square$

**Corollary (DERIVED) — gap-free recovery.** Let $\mathcal L_A$ be KMS-symmetric with $\mathcal L_A[\rho]=0$ and $F_A=\ker\mathcal L_A^\dagger\subseteq N_A$ (THEOREM for single-Pauli jumps, §1.4). Set $W=\rho^{-1/2}(\rho_{-A}-\rho)\rho^{-1/2}$. Then:
- **(a)** $W$ is KMS-orthogonal to all of $N_A$, hence to $F_A$;
- **(b)** $\lVert\mathcal R_{A,t}[\rho_{-A}]-\rho\rVert_1\le\min\big(\sqrt{c_*/t}\,\lVert W\rVert_{H^{-1}},\ t^{-1}\lVert W\rVert_{H^{-2}}\big)$, with $\lVert W\rVert^2_{H^{-1}}=\langle W,(-\mathcal L_A^\dagger)^{+}W\rangle_\rho$ and $\lVert W\rVert_{H^{-2}}=\lVert(-\mathcal L_A^\dagger)^+W\rVert_\rho$;
- **(c)** with the **relative Poincaré constant**
$$\kappa(A):=\inf_X\frac{\langle X,-\mathcal L_A^\dagger X\rangle_\rho}{\lVert X-\Pi_AX\rVert_\rho^2}\qquad(\Pi_A=\text{KMS-orthogonal projection onto }L^2(N_A)),$$
$$\lVert\mathcal R_{A,t}[\rho_{-A}]-\rho\rVert_1\le\sqrt{c_*/t}\;\chi_{\rm KMS}(\rho_{-A}\Vert\rho)\;\kappa(A)^{-1/2},\qquad\chi_{\rm KMS}:=\lVert W\rVert_\rho .$$

*Proof.*
- (a) For $Z=1_A\otimes z$: $\langle Z,W\rangle_\rho=\mathrm{Tr}[Z^\dagger(\rho_{-A}-\rho)]=0$, because $\rho_{-A}$ and $\rho$ have the same marginal on $A^c$.
- (b), (c): set-up.
  - Duality and stationarity give $\lVert\mathcal R_{A,t}[\rho_{-A}]-\rho\rVert_1=\sup_{\lVert X\rVert\le1}\lvert\mathrm{Tr}[(\mathcal R^\dagger_{A,t}X)(\rho_{-A}-\rho)]\rvert$. The supremum may be taken over Hermitian $X$.
  - Write $Y=\mathcal R^\dagger_{A,t}X=\phi_t(-\mathcal L_A^\dagger)X$. Then $\mathrm{Tr}[Y(\rho_{-A}-\rho)]=\langle Y,W\rangle_\rho$, and by (a) this equals $\langle Y-\Pi_AY,W\rangle_\rho$.
- (b) Apply D-6(ii) with $X\mapsto-\mathcal L_A^\dagger$ and $O=X$, using $\lVert X\rVert_\rho\le\lVert X\rVert$ ([CR] Lemma II.1).
- (c) Bound by $\lVert Y-\Pi_AY\rVert_\rho\lVert W\rVert_\rho\le\kappa^{-1/2}\langle Y,-\mathcal L^\dagger Y\rangle_\rho^{1/2}\lVert W\rVert_\rho$, then D-6(i). Note $\lVert W\rVert_{H^{-1}}\le\kappa^{-1/2}\lVert W\rVert_\rho$, because $W\in\mathrm{Ran}(1-\Pi_A)$. $\square$

*Checks.*
- **C5** (Davies sampler with single-site Pauli jumps, 5-qubit rings, $\beta=1.5$). The $H^{-2}$ bound holds at $t=1,10,10^2,10^3$ with ratios actual/bound of $0.07$–$0.19$, and both bounds hold throughout.
  - Classical ring, $\lvert A\rvert=1$: $\lVert W\rVert_{H^{-2}}=4.39$; errors $0.313$, $0.036$, $3.6\cdot10^{-3}$, $3.6\cdot10^{-4}$.
  - Transverse-field ring: $\lVert W\rVert_{H^{-2}}=9.44$; errors $0.668$, $0.173$, $0.0175$, $1.75\cdot10^{-3}$.
- **C9** (DLL/CKG-type Gaussian-filter sampler, open transverse-field chain, $\beta=1.3$, $t=100$). The (c) bound holds with errors $0.009$–$0.018$ against bounds $0.16$–$0.44$.

**Reading** (DERIVED, from the corollary). The Cesàro mean makes the Dirichlet energy small for free (i). Everything hard about recovery is the **inverse-Dirichlet size of the defect $W$** in (b), or the relative coercivity $\kappa(A)$ times $\chi_{\rm KMS}$ in (c). The defect is exactly orthogonal to $N_A$; that is the content of exact recovery at $t=\infty$ ([D-CR] §6.1).

[CR]'s proof bounds a Hölder surrogate of (c) unconditionally. It moves the commutators onto $Y$ through the twirl, pairs against $\rho$ itself (so no $\chi^2$ is paid), and compares string commutators with single-site Dirichlet energy (Lemmas IX.5, VIII.1, X.4). The price is an exponent $\lambda\ll1$ and the factor $2^{2\lvert A\rvert}$ ([D-CR] §3.0, §3.8). D-6 shows what an exponent-one coercivity would buy. It also shows the best gap-free rate: $t^{-1}$, against $W$'s $H^{-2}$ norm.

### 5.3 The literal "local gap" collapses; in the checks the relative one does not

**D-7 (DERIVED; checks C8, C10).** For any operator $X$,
$$\mathrm{gap}(\mathcal L_A)\cdot\mathrm{dist}_\rho(X,F_A)\le\lVert\mathcal L_A^\dagger X\rVert_\rho ,$$
where the gap is the smallest non-zero eigenvalue of $-\mathcal L_A^\dagger$ in the KMS geometry. Combine this with $\mathcal L^\dagger_{A,\ell}X=0$ for $X$ supported outside the $\ell$-patch (elementary: the truncated sampler's jumps and coherent term are supported in the patch) and with [CR] Lemma VII.3 (THEOREM, [D-CR] §3.9): $\lVert\mathcal L^\dagger_{A,\ell}-\mathcal L^\dagger_A\rVert_{\infty\to\infty}\lesssim\lvert A\rvert(e^{-c'\ell/(d\beta)}+2^{-\ell})$ for $\ell\ge4e^2\beta d$. Hence
$$\mathrm{gap}(\mathcal L_A)\ \lesssim\ \lvert A\rvert\big(e^{-c'\ell/(d\beta)}+2^{-\ell}\big)\,\frac{\lVert X\rVert}{\mathrm{dist}_\rho(X,F_A)}\qquad\text{for every }X\text{ supported at distance }\ge\ell\text{ from }A .$$

*Proof.*
- Let $P_0$ be the KMS projection onto $F_A$ and $Y=X-P_0X\perp F_A$.
- Then $\mathrm{gap}\,\lVert Y\rVert^2\le\langle Y,-\mathcal L^\dagger Y\rangle=\langle Y,-\mathcal L^\dagger X\rangle\le\lVert Y\rVert\lVert\mathcal L^\dagger X\rVert_\rho$.
- Finally $\lVert\cdot\rVert_\rho\le\lVert\cdot\rVert_\infty$. $\square$

**Consequence (DERIVED).** In the commuting case $F_A$ contains every far observable ([D-NC] K3: $F_A=1_A\otimes D_{\partial A}\otimes M_{\rm far}$), so the distance is $0$ and the bound is empty. When $F_A$ is trivial, as in the non-commuting examples of [D-NC] K3 and C8 (the ring of [D-CR] §6.1 also keeps the reflection fixing $A$), $\mathrm{dist}_\rho(X,F_A)=\lVert X-\rho(X)\rVert_\rho=O(1)$ for a far single-site $X$. **The gap above the exact kernel is then exponentially small in the distance from $A$ to the farthest site, with only the prefactor $\lvert A\rvert$ of Lemma VII.3.**

*Checks C8 and C10* (Gaussian-filter KMS sampler, open transverse-field chain, $\beta=1.3$):

| | gap above exact kernel | $\kappa(A)$ relative to $N_A$ | $\chi_{\rm KMS}(\rho_{-A}\Vert\rho)$ |
|---|---|---|---|
| commuting, $\lvert A\rvert=1$, $n=3,4,5$ | $0.352$, $0.352$, $0.352$ ($\dim F_A=8,32,128$) | $0.644$ (all $n$) | 2.0–2.2 |
| non-commuting, $\lvert A\rvert=1$, $n=3,4,5$ | $1.76\cdot10^{-2}$, $1.11\cdot10^{-3}$, $7.2\cdot10^{-5}$ ($\dim F_A=1$) | $0.415$, $0.263$, $0.227$ | 3.1–3.3 |
| non-commuting, $n=5$, $\lvert A\rvert=1,2,3$ | $7.2\cdot10^{-5}$, $1.0\cdot10^{-3}$, $1.8\cdot10^{-2}$ | $0.227$, $0.366$, $0.543$ | 3.25, 8.92, 30.8 |
| commuting, $n=5$, $\lvert A\rvert=1,2,3$ | $0.352$, $0.191$, $0.141$ | $0.644$, $0.191$, $0.141$ | 2.15, 4.03, 12.1 |

Readings of the table:
- In the non-commuting rows the gap depends only on the distance $\ell$ from $A$ to the far end: $\ell=2,3,4$ give $\approx1.8\cdot10^{-2}$, $1.1\cdot10^{-3}$, $7\cdot10^{-5}$, i.e. $\approx e^{-2.75\ell}$.
- The Rayleigh quotient of the far-end $Z$ with $F_A$ projected out, an upper bound on the gap by the min–max principle, is $7.6\cdot10^{-2}$, $8.2\cdot10^{-3}$, $7.6\cdot10^{-4}$.
- The recovery error at $t=100$ stays at $0.017$–$0.018$ for $n=3,4,5$ (C9), although the gap falls 250-fold. **Recovery is not mixing.**

**Consequence for [CR] App. B (DERIVED, on the digest's reading).** [D-CR] §2 records Def. B.1 as $-\lambda\langle X,\mathcal L^\dagger X\rangle_\rho\le\langle X,\mathcal L^{\dagger2}X\rangle_\rho$ for all $X$, uniformly over restricted Gibbs states. For a KMS-self-adjoint $\mathcal L^\dagger\le0$ this is exactly $\lambda\le\mathrm{gap}$ above the exact kernel (diagonalise: $\lambda x\le x^2$ on the spectrum).
- By D-7 the hypothesis cannot hold uniformly in the size of the restricted region for a non-commuting $H$ whose single-region fixed-point algebra is trivial.
- This *explains* [CR]'s remark that "we do not know of any a priori bound on the local gap, even assuming high temperature" for non-commuting $H$. It does not contradict any theorem: Cor B.2 is a correct implication with a hypothesis that, read literally, is not available there.
- [CR]'s App. B text itself was not re-read for this note, only the digest's transcription; if B.1 is meant only on patches of radius $O(d\beta)$, this consequence does not apply.

### 5.4 Where the expander essence would have to enter: two conjectures

**CONJECTURE C-1 (relative local gap).** For bounded-degree $H$ and fixed $\beta$ there are $c_\beta>0$ and $p\ge0$ such that the single-Pauli KMS sampler satisfies $\kappa(A)\ge c_\beta\lvert A\rvert^{-p}$, uniformly in the system size and in the shape of the region.
- *What must be true:* the dressed single-site derivation controls, linearly, the KMS distance of any operator from the algebra $1_A\otimes M_{A^c}$. That is an exponent-one version of [CR] Lemma X.4 combined with the twirl.
- *Evidence:* C9–C10, for $n\le5$ only. For the non-commuting chain, $\kappa$ decreases with $n$ by shrinking steps ($0.415$, $0.263$, $0.227$ for $n=3,4,5$, consistent with a positive limit) and does not decay with $\lvert A\rvert$ ($0.227$, $0.366$, $0.543$). For the commuting chain it falls with $\lvert A\rvert$ ($0.644$, $0.191$, $0.141$ for $\lvert A\rvert=1,2,3$; three points fix no rate) and does not depend on $n$.
- *What it buys* (DERIVED from D-6(c)): recovery error $\le\sqrt{c_*/t}\,\chi_{\rm KMS}\,c_\beta^{-1/2}\lvert A\rvert^{p/2}$, i.e. the rate $t^{-1/2}$ in place of $t^{-\lambda}$, $\lambda\approx1/(2d^4\beta^4)$.
- *What it does not buy:* $\chi_{\rm KMS}(\rho_{-A}\Vert\rho)$ grows by a factor of about 2–3.5 per site of $A$ in C10 (2.15→4.03→12.1 commuting, 3.25→8.92→30.8 non-commuting) and C11. Classically it is $\big(\sum_x\mu_{-A}(x)^2/\mu(x)-1\big)^{1/2}$, exponential in $\lvert A\rvert$ unless the conditional law on $A$ is uniform. So the $\chi^2$ route leaves an exponential prefactor in $\lvert A\rvert$, as an $L^2$ warm start does for Markov chains ([D-NC] §2.8, $\sqrt{1/\sigma_{\min}}$ against $\sqrt{2\log(1/\sigma_{\min})}$). With the gap-free rate $t^{-1/2}$ this prefactor costs a time $t\sim\chi_{\rm KMS}^2$, exponential in $\lvert A\rvert$; with an exponential rate $e^{-\kappa t}$ it would cost only a time $\log\chi_{\rm KMS}/\kappa$, polynomial in $\lvert A\rvert$ when $\kappa\ge c_\beta\lvert A\rvert^{-p}$.

**CONJECTURE C-2 (relative modified log-Sobolev ⇒ global Markov).** Suppose the entropy production of the region's semigroup dominates the relative entropy *conditional on $N_A$*:
$$\mathrm{EP}_{\mathcal L_A}(\omega)\ge\alpha\,D_A(\omega\Vert\rho),\qquad D_A(\omega\Vert\rho)=D(\omega\Vert\rho)-D(\omega_{A^c}\Vert\rho_{A^c}),\qquad\alpha\ge c_\beta\lvert A\rvert^{-p}.$$
Here $D_A$ is the conditional relative entropy that [D-Y] §10 identifies with [Y]'s $\delta_{A^c}$. Suppose this holds uniformly over regions and restricted Gibbs states. Then the semigroup itself, not its Cesàro mean, recovers: $\lVert e^{t\mathcal L_A}[\rho_{-A}]-\rho\rVert_1\le\mathrm{poly}(\lvert A\rvert)\,e^{-\alpha t/2}$, up to the Lieb–Robinson truncation. As in [CR] App. B (Lemmas B.1–B.2), the exponential decay is used directly. With the truncation and continuity steps of [CR] Cors III.1–III.2, this gives $I(A{:}C|B)\le\mathrm{Poly}(\lvert A\rvert,\lvert C\rvert)e^{-\mathrm{dist}(A,C)/\xi}$: the global Markov property.
- *What must be true:*
  - (i) The starting entropy is linear, not exponential. **DERIVED:** $D(\rho_{-A}\Vert\rho)\le2\beta\sum_{\gamma:\,\mathrm{supp}\,h_\gamma\cap A\ne\emptyset}\lVert h_\gamma\rVert=O(\beta d\lvert A\rvert)$.
    - $D(\rho_{-A}\Vert\rho)=\beta(\mathrm{Tr}\rho_{-A}H-\mathrm{Tr}\rho H)-(S(\rho_{-A})-S(\rho))$.
    - The entropy difference is $\ge0$ by subadditivity: $S(\rho)\le\lvert A\rvert\log2+S(\rho_{BC})=S(\rho_{-A})$.
    - Terms of $H$ supported off $A$ have equal expectations, and each other term changes by at most $2\lVert h_\gamma\rVert$.
    - *Check C11* (open transverse-field chain, $n=6$, $\beta=1.3$, $\lvert A\rvert=1,\dots,4$): $D(\rho_{-A}\Vert\rho)=0.71,1.38,2.06,2.75$ (commuting) and $1.04,2.04,3.04,4.05$ (non-commuting), linear and within the bound. Over the same range $\chi_{\rm KMS}=2.0,4.4,10.9,24.0$ and $3.2,9.2,29.2,87.6$, exponential.
    - This is the same gain as $\sqrt{2\log(1/\sigma_{\min})}$ against $\sqrt{1/\sigma_{\min}}$ in [D-NC] §2.8.
  - (ii) The drift of the $A^c$-marginal under $e^{t\mathcal L_A}$ is controlled. The classical heat-bath generator vanishes on $N_A$; $\mathcal L_A^\dagger$ vanishes only on $F_A\subseteq N_A^\sigma$ and is merely quasi-local on the rest of $N_A$.
  - (iii) A Pinsker step from $D$ to trace distance.
- *Status:*
  - The classical analogue is the entropy-factorisation route of [CLV] and [BCR] ([D-HDX] §2.7, §4.5).
  - [D-HDX] §4.5 reports that no $k$-step quantum trickle-down was found.
  - This conjecture is where "expander essence" would make the global quantum Markov property a theorem. What it asks for is a relative, entropic, hereditary functional inequality, not a spectral gap (D-7 rules the literal gap out).

### 5.5 Yang's row in the same shape

- **THEOREM** ([D-Y] §5.3–5.4).
  - $\delta_R(\rho,\sigma)=\int_0^\infty j(s)\,ds$, where $j(s)=\inf_{y\in L^2(N_R)\rho^{1/2}}\lVert(\Delta+s)^{1/2}(y-(\Delta+s)^{-1}\rho^{1/2})\rVert^2$ is a weighted squared distance from the resolvent vector to the subalgebra's subspace (Carlen–Vershynina).
  - $(\Delta+s)^{-1}\rho^{1/2}=B_s\rho^{1/2}$, with $B_s$ a $1/\sinh(\pi t)$-weighted average of the Connes cocycle $u_t$.
  - Hence $j(s)\le(1+s)q^2/s^2$.
- **THEOREM** ([D-Y] §5.6). For the cut pair, $u_t$ is the interaction-picture propagator of $V$. Lieb–Robinson growth $e^{v\beta\lvert t\rvert}$ against the kernel decay $e^{-\pi\lvert t\rvert}$ gives $q\lesssim(1+\beta g)e^{-\pi\mu_{LR}(r-R_0)/(\pi+\beta v)}$.
- In the language of §1.3, [Y] measures **defect II for the pair (state, cut reference): the operator-norm distance of the Connes cocycle from $N_{AB}$**, averaged over modular time. ANALOGY, as a reading of [CR] against [Y]; [D-CHAT] §7.2 makes the same move, as its own interpretation ("the static proof does not do away with the time average; it moves it", from Lindblad time to modular time). *Break:* the Lindblad-time average is a Cesàro mean of a contraction semigroup, the modular-time average a kernel integral of a unitary cocycle, and neither is obtained from the other.

### 5.6 The modular-time kernels have KMS-strip Fourier transforms

**DERIVED** Fourier identities (check C7, 8 digits), with an **ANALOGY** between the three kernels ([D-Y] B14, which labels it so; [D-HC] §8.2):
- [CR]: $\int g(t)e^{-ixt}dt=\frac1{2\cosh(\beta x/4)}$ for $g(t)=\frac1{\beta\cosh(2\pi t/\beta)}$ ([CR] Lemma X.2).
- [JRSWW]: $\int\beta_0(t)e^{-ixt}dt=\frac x{\sinh x}$ for $\beta_0(t)=\frac\pi2(\cosh\pi t+1)^{-1}=\frac\pi4\mathrm{sech}^2(\pi t/2)$.
- [Y]: $\int_0^\infty\frac{\sin vt}{\sinh\pi t}dt=\frac12\tanh\frac v2$ ([Y] App. B.1).

All three are Fourier transforms of functions analytic in a horizontal strip whose width is fixed by the KMS condition. That is the one feature the three "time averages" over modular time share. *Break:* the kernels do different jobs. [Y]'s $1/\lvert\sinh\pi t\rvert$ beats Lieb–Robinson growth; [CR] use $g$ only near $t=0$ (Lemma X.4), and their light cone is cut by the Gaussian filter instead ([D-Y] §6); [JRSWW]'s $\beta_0$ is a probability density over rotated Petz maps. *Proof of the second identity:* $\int\mathrm{sech}^2(at)e^{-ixt}dt=\pi x/(a^2\sinh(\pi x/2a))$ with $a=\pi/2$. $\square$

### 5.7 The gap row in coarse geometry

**KNOWN-LINK** ([D-EXP] §9.11, Willett–Yu). With a spectral gap, the projection onto constants is $f(\Delta)$ for a continuous $f$ equal to $1$ at $0$ and to $0$ on $[\lambda,2]$ ($\lambda$ the gap). It is therefore a norm limit of finite-propagation operators (the basic Kazhdan projection), and a ghost. *From memory* (Chebyshev acceleration): a polynomial of degree $k$ achieves norm error $\le2e^{-k\sqrt{2\lambda}\,(1+o(1))}$, so the optimal propagation is $\sim\lambda^{-1/2}\log(1/\epsilon)$.

**ANALOGY** (to [CR]). D-6 is the gap-free row. The Cesàro mean approximates the projection only in energy norm, never in operator norm without a gap. That suffices for recovery because recovery tests against a single vector $W$, which is orthogonal to $L^2(N_A)$ and has controlled inverse-Dirichlet size. *Break:*
- Roe algebras concern operators on $\ell^2(X)$ over a coarse space, and the projection onto constants.
- Here the operators are superoperators with Lieb–Robinson tails, and the target is a subalgebra $N_A$ whose trace-preserving projection is already local.
- The non-local object is the $\rho$-preserving expectation that does not exist, so there is no ghost phenomenon to transfer.

---

## 6. The time-averaged detailed-balanced Lindbladian with single-Pauli jumps on $A$

### 6.1 Word by word

| word | mathematical content | label | where |
|---|---|---|---|
| detailed-balanced | KMS-symmetric, so a Dirichlet form, so the squared norm of a derivation twisted by $\sigma_{\mp i/4}$; explicitly, Gaussian-filtered single Paulis dressed in modular time with kernel $g$. GNS symmetry needs jumps that are exact Bohr-frequency eigenoperators (Alicki), generically non-local for non-commuting $H$, and [CKG]'s filtered quasi-local construction cannot be GNS-symmetric | THEOREM | §1.4; [D-NC] §3.2, §4.3 ([CKG] App. E) |
| single-Pauli jumps | a generating set of $M_A$, so $\ker\mathcal L_A^\dagger\subseteq\{P^1_A\}'=N_A$: the vertical derivation of the inclusion $N_A\subset M$ | THEOREM | [D-CR] §3.5 (4.1) |
| on $A$ | quasi-local by Lieb–Robinson; truncation error linear in $t$ | THEOREM | [D-CR] §3.9 |
| time-averaged | Cesàro mean $\phi_t(-\mathcal L_A)$; energy $\le c_*/t$ with no gap; limit $=$ KMS projection onto $F_A\subseteq N_A^\sigma$ (modular core). The limit recovers exactly but is non-local for non-commuting $H$ | THEOREM (energy $\le2/t$, [CR] Cor VII.1), DERIVED (the constant $c_*$, D-6(i); the limit) | [D-CR] §3.7, §6.1 |
| classical shadow | block heat bath $=$ resampling $A$ from its link $=$ one block step of the down-up walk | DERIVED | [D-HDX] §4.4 (check C4 there) |

**In one sentence (DERIVED, assembling the rows).** $\mathcal R_{A,t}$ is the Cesàro mean of the heat flow of the KMS-twisted derivation of $N_A=M_{A^c}\subset M$. It is a gap-free, quasi-local, exactly $\rho$-preserving surrogate for a $\rho$-preserving conditional expectation onto $N_A$, which does not exist (§4.1). It succeeds on the one input that matters because the defect $W=\rho^{-1/2}(\rho_{-A}-\rho)\rho^{-1/2}$ is exactly orthogonal to $L^2(N_A)$ and is controlled by Dirichlet energy relative to $N_A$ (D-6).

### 6.2 At $\beta=0$ it is a Cayley Laplacian with its Kazhdan projection

**D-8 (DERIVED; check C6).** Let $\beta=0$, so $\rho\propto1$ and the [CR] generator is, up to a constant, $\mathcal L_A(X)=\sum_{i\in A,P}(P_iXP_i-X)$ ([D-NC] §3.4).
- (i) $\mathcal L_A$ is the generator of the Cayley-graph walk of the conjugation action of the Pauli group $G_A\cong\mathbb Z_2^{2\lvert A\rvert}$ with the $3\lvert A\rvert$ single-site generators.
- (ii) The twirl $E^\tau_{N_A}=4^{-\lvert A\rvert}\sum_S S\cdot S$ is the Haar average, i.e. the Kazhdan projection onto invariant vectors.
- (iii) The gap of $-\mathcal L_A$ (a *sum* over generators) is $4$ and $\mathrm{Var}_A\le\frac18\mathcal E_1$ (Efron–Stein), independent of $\lvert A\rvert$. The Kazhdan constant in the max-over-generators form of HLW Def 11.18 ([D-EXP] §1.4), taken in the regular representation, is $K(G_A,P^1_A)=8/(3\lvert A\rvert)$, and the normalised walk $\frac1{3\lvert A\rvert}\sum_s\mathrm{Ad}\,s$ has gap $4/(3\lvert A\rvert)$.
- (iv) The "string" Dirichlet form $\mathcal E_{\rm str}(O)=4^{-\lvert A\rvert}\sum_S\lVert[S,O]\rVert_2^2$ (the complete Cayley graph, a perfect expander) equals $2\mathrm{Var}_A(O)$. The canonical-path comparison along words gives $\mathcal E_{\rm str}\le\frac{\lvert A\rvert}4\mathcal E_1$.

*Proof.*
- (i), (ii) Conjugation by Paulis acts diagonally on Pauli strings $T$ through the characters $\chi_T(S)=\pm1$ (commute or anticommute).
- (iii)
  - A non-trivial character $\chi_T$ equals $-1$ on at least two single-site generators: at a site where $T_i\ne1$, two of the three Paulis anticommute with $T_i$. So the Laplacian eigenvalue on $T$ is $4w_A(T)\ge4$.
  - Kazhdan constant. For a unit vector $v\perp1$ with weights $p_\chi=\lvert v_\chi\rvert^2$ on characters, $\lVert sv-v\rVert^2=4P_p(\chi(s)=-1)$. Hence $\max_s\ge\frac4{3\lvert A\rvert}\mathbb E_p\#\{s:\chi(s)=-1\}\ge\frac8{3\lvert A\rvert}$, with equality for $p$ uniform on the $3\lvert A\rvert$ single-site strings (each generator $P_i$ anticommutes with exactly two of them).
  - $\sum_{P}\lVert[P_i,T]\rVert_2^2=8\lVert T\rVert_2^2$ if $T_i\ne1$ (two of the three Paulis anticommute), so $\mathcal E_1(T)=8\,w_A(T)\lVert T\rVert_2^2$ while $\mathrm{Var}_A(T)=\lVert T\rVert^2_2$ for $w_A(T)\ge1$.
- (iv) $4^{-\lvert A\rvert}\sum_S\lVert[S,O]\rVert_2^2=2\lVert O\rVert^2-2\langle O,E^\tau O\rangle$.
  - Leibniz: $[S,O]=\sum_k(\prod_{i<k}S_i)[S_k,O](\prod_{i>k}S_i)$, and Hilbert–Schmidt norms are unitarily invariant, so $\lVert[S,O]\rVert_2^2\le w\sum_k\lVert[S_k,O]\rVert_2^2$.
  - For fixed $i$ and fixed $P$, the number of strings with $S_i=P$ is $4^{\lvert A\rvert-1}$, and $w\le\lvert A\rvert$. $\square$

*Check C6:* $\mathcal E_{\rm str}/2=\mathrm{Var}_A$ to $10^{-14}$; $\max\mathrm{Var}_A/(\mathcal E_1/8)=0.47$ on random operators ($\lvert A\rvert=3$).

**ANALOGY** (for $\beta>0$). [CR]'s proof is a **Diaconis–Saloff-Coste comparison** (*from memory*) between the complete-graph Dirichlet form, which the twirl provides, and the single-site Dirichlet form, which the dynamics provides, carried out in the KMS geometry. *It breaks* exactly where the KMS norm stops being invariant under multiplication by Paulis. Each Leibniz prefix and suffix then costs an imaginary-time conjugation $\lVert\rho_{\beta_0}S\rho_{\beta_0}^{-1}\rVert\le2^{w}$ ([CR] Cor IX.1, THEOREM, [D-Y] §8.8), and frequency filtering (Lemma IX.5) is needed to keep that cost finite at low temperature, at the price of a Hölder exponent. **On this reading [CR]'s $2^{2\lvert A\rvert}$ is an upper bound on (a Hölder form of) the comparison constant of the Pauli group's two Cayley graphs in the twisted geometry; in the tracial geometry the constant is $\lvert A\rvert/4$ (D-8).** The gap-free Hölder-type inequality *is* what [CR] Lemma X.4 proves; the comparison-theorem reading is the analogy.

**Guard** (THEOREM, [D-NC] §4.4, §5.3). The $\beta=0$ generator is not a quantum-expander family. Its degree $3\lvert A\rvert$ is unbounded, its Kazhdan constant is $8/(3\lvert A\rvert)\to0$ (D-8), and its $\lvert A\rvert$-independent gap comes from tensorization. Bounded-degree quantum expanders (Hastings; Ben-Aroya–Schwartz–Ta-Shma) are a different object.

### 6.3 What is genuinely new, and what is not

- **Not new: the Cesàro mean.** It is the mean ergodic theorem with a rate in Dirichlet form (D-6(i)). [HJPW]'s Lemma 12 uses the Cesàro mean of a channel to produce the conditional expectation onto its fixed-point algebra ([D-HC] B2), and [LMRRW] use the classical gap-free version $\mathbb E_t\mathcal E(f_t,\log f_t)\le\mathrm{KL}/T$ ([D-Y] B3). KNOWN-LINK.
- **Not new: exact detailed balance with quasi-local jumps.** That is [CKG]. KNOWN-LINK.
- **New, by [CR]'s own account** ([D-CR] §4, THEOREM): a connection between the Dirichlet form, a dynamic quantity, and commutators in the KMS inner product, a static quantity. Concretely, Lemma X.4 ("small Dirichlet form ⇒ small commutator with each jump") holds at every temperature thanks to the Bohr-frequency (Gaussian) regularisation of imaginary-time conjugation (Lemma IX.2: a Gaussian filter beats any $e^{\beta\omega}$).
- **New, the idea, in this note's language** (DERIVED for the structure, via D-6, D-7 and checks C8–C9; the attribution of intent is ANALOGY):
  - Recovery of a region needs coercivity of the region's derivation *relative to the algebra that forgets the region*, not mixing of the region's dynamics to its exact fixed points.
  - [CR] supply that relative coercivity in Hölder form, unconditionally.
  - In expander terms: the classical locally-stationary theory turns small Dirichlet energy into correct conditionals only under a *hereditary* modified log-Sobolev inequality over all pinnings ([LMRRW] Lemma 3.5, [D-CR] §7.1). [CR] show that, for recovering one region, an *unconditional, non-hereditary, single-region* Hölder coercivity suffices, at a price exponential in $\lvert A\rvert$.
  - The exact-kernel gap that a "mixing" proof would need collapses exponentially with distance (D-7), while in the checks the relative constant does not (C9–C10, $n\le5$; C-1 conjectures this in general). That is why the gap-free route is not merely convenient here; it is the right route.
- **New in [Y], the same theorem seen statically** (THEOREM [Y]; reading in §4.3, §5.5): the dynamics is replaced by the modular flow of the pair (state, cut reference), the erase reference by the cut reference (index → interface), and coercivity by the light cone of the Connes cocycle.

---

## 7. What noncommutative geometry contributes

| NCG layer | what it supplies here | label | the exact point of contact or break |
|---|---|---|---|
| modular theory (Tomita–Takesaki, Connes cocycle, KMS) | Markov $=$ cocycle localisation (D-4); quantitative form with index (D-5); [Y]'s leakage; Takesaki's obstruction (§4.1); KMS-strip kernels (§5.6) | THEOREM / DERIVED | literal; the quantum Markov theory *is* modular theory of the inclusion $M_{AB}\otimes1\subset M$ for a pair of states |
| inclusions and subfactor geometry (index, Pimsner–Popa, commuting squares) | $\rho\le d_A^2\rho_{-A}$ gives D-5's $d_A$ (and equals [CR]'s $2^{2\lvert A\rvert}$ in value, though not in origin, §4.3); classical Markov $=$ commuting square (D-1a); expansion $=$ Friedrichs angle (D-1, D-2) | DERIVED; names KNOWN-LINK (*from memory*) | literal in finite dimension. *Break:* the subfactor theory of infinite index or type III, which the names suggest, is not used |
| Dirichlet forms as squared derivations | the [CR] Dirichlet form, explicitly (T10); kernels as commutants; $\kappa(A)$ as relative coercivity of the derivation (§5.4) | THEOREM / DERIVED | literal |
| property (T), Kazhdan constants, amenability | (T) ⇔ derivations inner (T11); at $\beta=0$ the single-Pauli generator is the Cayley Laplacian of the Pauli group, twirl $=$ Kazhdan projection, Laplacian gap 4, Kazhdan constant $8/(3\lvert A\rvert)$ (D-8) | THEOREM / DERIVED | literal at $\beta=0$. ANALOGY at $\beta>0$: a "twisted Kazhdan constant" for the KMS-twisted bimodule; *break:* [Pet] is proved for tracial finite factors and untwisted bimodules |
| coarse geometry (Roe algebras, Kazhdan projections, ghosts) | gap ⇒ global projection is a norm limit of local operators (T12); the gap row of §5.1 | THEOREM (expanders); ANALOGY (Markov) | *break:* §5.7; the quantum object to localise is a non-existent conditional expectation, not a projection onto constants |
| spectral triples, Connes distance, Rieffel Lip-norms | the $L^\infty$ face of the same derivation; Lip-norm iff primitive; radius $\le\sqrt{N/(\lambda\rho_{\min})}$ (gap ⇒ small diameter); at $\beta=0$ the single-Pauli Connes distance is De Palma–Marvian–Trevisan–Lloyd's $W_1$ to within $\tfrac23$–$\tfrac32$ | DERIVED in [D-NC] §3.3–3.4 | construction plus constants. In finite dimension the axioms are automatic and neither Markov proof uses the metric |
| Carlen–Maas $W_2$ | gradient-flow geometry for GNS-symmetric semigroups | THEOREM ([D-NC] L8) | **absent** for the KMS-only [CKG]/[CR] samplers. Ricci or Talagrand statements about them are unsupported |
| $K$-theory, cyclic cohomology, modular spectral triples | none for Markov or expansion in finite dimension; [CNNR]'s modular index might label sufficiency defects for stationary diagrams | SPECULATION ([D-NC] B-NC3) | finite-dimensional NCG topology is trivial (Note 1 §10) |

**Verdict.** NCG ties the five areas together where the ties are *operator-algebraic*: the measure-theoretic layer (modular theory, inclusions) and the first-order layer (derivations and Dirichlet forms; Kazhdan constants and Roe-algebra locality for groups and expanders). Those ties are literal. Carried to the $\beta>0$ Markov problem, the Kazhdan and Roe-algebra rows are ANALOGIES (§5.7, §6.2). At the metric and topological layers the contribution is a reformulation with constants, or speculation. This agrees with [D-HDX] §4.6, [D-HC] §13.1 and [D-NC] §9, reached there by other routes (rule C3).

---

## 8. The user's intuition, tested

**Where the "same essence" is literally a theorem.**
1. The pairwise influence spectrum of a Gibbs measure under all pinnings *is* the link spectrum of a high-dimensional expander (T2), and that is a frame bound on single-site fluctuation spaces (D-3).
2. The Markov property, strong spatial mixing and the expander mixing lemma are values of one functional, $c(A,B)$, and Dobrushin influences bound it on adjacent pairs (D-1). The step from local to global is the sharp two-projection inequality (D-2) and its iterations (T3, T5, T6).
3. At $\beta=0$ the "single-Pauli jumps on $A$" generator is the Cayley Laplacian of the Pauli group of $A$, and the twirl is its Kazhdan projection (D-8). The $\lvert A\rvert$-independent gap is a tensorization (sum) phenomenon, while the Kazhdan constant decays like $1/\lvert A\rvert$.
4. A uniform local gap would turn the local quantum Markov property into the global one (T13, a theorem as stated). For non-commuting $H$ with trivial $F_A$ its hypothesis, read literally, fails (D-7). Read with the relative constant $\kappa(A)$ instead, the statement is no longer a theorem: through D-6(c) it keeps an exponential $\chi_{\rm KMS}$ prefactor (C-1), and the entropic version is CONJECTURE C-2.

**Where it is two halves of one mechanism.** Every local-to-global theorem in these areas pairs:
- **heredity**: the class is closed under passing to links or pinnings, and on separated pairs the angle is exactly zero (Markov / Hammersley–Clifford; exact, combinatorial, no rate);
- **uniform decorrelation**: the angle is bounded away from 1 on overlapping pairs or links (expansion, spectral independence, strong clustering, Dobrushin).

[D-HC] B3 shows the hinge: Hammersley–Clifford's positivity secures "angle $<1$" (the lattice law, §3.4), the gap is "angle $\le c_0<1$". Quantumly the two halves are coupled, since heredity (strong Markov) costs clustering (T14).

**Where it is false.**
1. *"The Markov property is an expansion property."* False. $c=0$ on separated pairs at every temperature, while the buffer angle is $0.93$ (C1, $\beta=2$). Classical Gibbs measures stay Markov at phase transitions ([D-Y] §10). Both all-temperature quantum proofs use no gap.
2. *"Expansion of the interaction graph helps Markov structure."* False for the static bound. Yang's improvement from volume to boundary is a Cheeger-numerator gain and vanishes on expanders ([D-Y] B1). At low temperature expansion makes the bottlenecks that kill mixing ([D-EXP] §9.7).
3. *"The time-averaged single-Pauli Lindbladian is a quantum expander."* False: unbounded degree, a tensorization gap, and Hastings' quantum expanders are another object (§6.2 guard).
4. *"Its essence is mixing."* False: recovery succeeds at $t=100$ while the exact-kernel mixing time is $\gtrsim10^4$ (C9, $n=5$). What is needed is relative coercivity (D-6, D-7).
5. *"Yang's static proof hides an expander argument."* False. It is isoperimetry plus Lieb–Robinson plus modular theory ([D-CHAT] §7.3).

---

## 9. For the programme, and the competition context

### 9.1 Abstract objects (prong 1; stated for inclusions $N\subset M$ with a faithful state)

- **The history corner is the commuting case. The literal local gap is the right notion there** (DERIVED from [D-NC] B-NC1 and D-6).
  - Note 1's KMS states $\mathbb P_\beta$ have modular flow equal to weighted depth, which fixes the face algebra $D_0$. So $D_0$ is modular-invariant and the Takesaki expectation onto it exists: defect II's obstruction is absent.
  - The re-routing Lindbladian's Cesàro limit is that expectation ([D-NC] §6.1), provided at most one endpoint fibre is a single history ([D-NC] §6.1(b)); otherwise the kernel is larger than $D_0$.
  - By D-6, for any state $\omega$ on one origin block, $\lVert\mathcal R_t(\omega)-E_*\omega\rVert_1\le\min\big(\sqrt{c_*/t}\lVert W\rVert_{H^{-1}},\,t^{-1}\lVert W\rVert_{H^{-2}}\big)$ with $W=\mathbb P_\beta^{-1/2}(\omega-E_*\omega)\mathbb P_\beta^{-1/2}$.
  - *Proof.* $E\,\mathcal L^\dagger=0$ gives $E\mathcal R_t^\dagger=E$. So $\mathrm{Tr}[X(\mathcal R_t\omega-E_*\omega)]=\mathrm{Tr}[(\mathcal R_t^\dagger X-EX)(\omega-E_*\omega)]=\langle\mathcal R_t^\dagger X,W\rangle_{\mathbb P_\beta}$, using $\mathrm{Tr}[(EX)(\omega-E_*\omega)]=\mathrm{Tr}[XE_*(\omega-E_*\omega)]=0$. Moreover $W\perp D_0=\ker\mathcal L^\dagger$ by the same identity, so D-6(ii) applies. $\square$
  - With a gap $\lambda$ of the re-routing (or flip) generator, $\lVert W\rVert_{H^{-1}}\le\lambda^{-1/2}\lVert W\rVert_{\rm KMS}$. This is [D-EXP] §11.1's Poincaré bound in dynamical form.
- **Two defects, two measurements** (CONJECTURE-level programme content, stated for abstract faces).
  - **Defect I.** Define the angle of the triple (past algebra, future algebra | face algebra at layer $\ell$) inside the history algebra $D$, i.e. the conditional maximal correlation of past and future given the present face. It is zero iff the law is Markov at $\ell$ (D-1a; U8 of [D-HDX] §6.4). Unlike the CMI it is a sup-type, $L^2$ quantity, and its uniform smallness is the expander-type certificate.
  - **Defect II.** For states with coherences between histories (Note 1 §4.1; [D-Y] B13), define the leakage of $(D\omega{:}D\mathbb P_\beta)_t$ out of $D_0$. [Y]'s Cor III.5 in its subalgebra form ([D-CHAT] §7.1: proof sketch and numerical checks there) bounds the sufficiency defect of $D_0$ for the pair $(\omega,\mathbb P_\beta)$ by this leakage, with the Rényi moment $\mathrm{Tr}\,\omega^{1+\alpha}\mathbb P_\beta^{-\alpha}$ as size factor; an index replaces it only for an erase-type reference, as in D-5. [D-Y] B8 already gives $\eta\le\lvert\beta-\beta'\rvert\lvert t\rvert\max_\tau\mathrm{osc}F$ for diagonal states.
  - Both pass the drag test: they are defined for any inclusion with a faithful state.
- **Arrows.** An arrow is a linear layer together with the ReLU immediately before it (Note 1 §2.1). Nothing above uses activation vectors. Faces enter only as the projections generating $D_0$, and the angle is computed in $L^2$ of the state on histories.

### 9.2 Competition context (prong 2; dictionary-level, rule P2.1)

- **Annealed against quenched** (DERIVED in [D-TS] §2.5). The ensemble-averaged covariance transfer is a conditional expectation onto the scalars (all non-trivial eigenvalues $0$: a perfect expander), and the averaged third-order transfer is $0$. A single realised layer is a one-Kraus channel, with no averaging and therefore no isolated invariant subspace.
- In this note's language the old $\kappa_3$ content lives in what the conditional expectation does not see. The propagator law $\mathrm{PR}\approx n/(2\cdot\text{age})$ has no gap, so mode counts are a fixed fraction of $n$ ([D-TS] §3.4), consistent with the old-content stream's verdict that none of the carriers it tested is cheap enough at $n=1024$ ([old-content/REPORT.md](../old-content/REPORT.md)).
- Per P2.1 that negative is charged to the dictionary (N1: the $\Phi^3$ pass-through; N2: the noAD bookkeeping), not to the theory. The only gapped direction is the rank-one mean mode ([D-TS] B6).
- **The face law** has measured unpinned $\eta_0=1.7$–$8$ ([D-TS] §3.5). By D-3 this is a frame bound $\lVert\sum P_i\rVert-1$ on gate fluctuation lines. Pinned values, the actual hypothesis of T2–T3, need third-order gate statistics that the atlases do not store.
- **Honest scope.** This synthesis gives **no direct lever on the binding competition error**. Per [competition-plan.md](../../competition-plan.md) §6b and the [oracle1024](../oracle1024/REPORT.md) and [chain128](../chain128/REPORT.md) reports, that error is the carrier of the $(2,1,1)$ fourth-cumulant slice. The closest contact is §9.2's first bullet, which predicts the failure of expander-style compression of quenched content and agrees with the measurements. Nothing here should change Lines A–C.

---

## 10. Conjectures and open problems, precisely

1. **C-1 (relative local gap)**, §5.4: $\kappa(A)\ge c_\beta\lvert A\rvert^{-p}$ uniformly. Test: extend C9–C10 to $n=6,7$ on a 2D patch, and to the [CR] Metropolis weight in place of the Gaussian filter.
2. **C-2 (relative MLSI ⇒ global Markov)**, §5.4. The quantum $k$-step trickle-down that [D-HDX] §7 calls for would be its proof strategy: certify relative MLSI on two-site links, then trickle down.
3. **C-3 (sharpness of D-5's power).** For Gibbs states of bounded-degree $H$, $I(A{:}C|B)\le C\,q^2\log(1/q)$ for the erase-reference leakage. Evidence: CMI$/q^2\in\{1.84,0.85\}$ in C4. *What must be true:* the small-resolvent tail of [Y] Lemma III.4 is $O(q^2\log(1/q))$ for the pair $(\rho,\rho_{-A})$, i.e. the negative spectral tail of $\log\Delta_{\rho_{-A},\rho}$ is controlled by $\log d_A^2$.
4. **C-4 (angle form of approximate quantum intersection).** Approximate quantum intersection ([D-HC] B7) holds with constant $(1-c_\rho)^{-1}$, where $c_\rho$ is the Friedrichs cosine of $L^2(M_{ACD}\otimes1)$ and $L^2(M_{ABC}\otimes1)$ in the **KMS** geometry of $\rho$. *What must be true:* D-2 applies verbatim, with "$\mathrm{Var}$" replaced by the squared KMS distance to the subalgebra subspaces (KMS-orthogonal projections exist even though conditional expectations do not), and the CMI can be bounded by those distances. The second step is the open part.
5. **SPECULATION.** The defect-II leakage $q$ for the cut reference plays the role, in the global Markov problem, that the angle plays in the gap problem. A "coarse" global Markov statement might follow from an angle bound between the cocycle orbits of overlapping cuts, rather than from any spectral gap. [Y] §VII's 2D rescaling corollary (7.2) is the only result of this type known here.

---

## 11. Messages to the other prongs (rule C1)

**To prong 1 (theory).**
1. Record the **two defects** of an inclusion with a faithful state: the angle $c(N_1,N_2)$ (§1.3) and the cocycle leakage $q_N$. Record D-1 (Markov $=$ zero angle on separated pairs, classically) and D-4 (Markov $=$ cocycle localisation, quantum). Both are stated for abstract inclusions.
2. Note 1's KMS family is the commuting case. There the Takesaki expectation onto $D_0$ exists, the exact-kernel gap is the right local gap, and D-6 bounds the re-routing recovery by inverse-Dirichlet norms (§9.1).
3. For any future extension with coherent re-routing energies: the literal local gap will collapse (D-7), and the right constant is the Poincaré constant relative to the face algebra, $\kappa$.

**To prong 2 (bridge).**
1. Measure the conditional maximal correlation of (past faces, future faces | present face) per layer. It is the $L^2$ Markov defect (D-1a) and complements E2's CMI. Calibrate with the Markov surrogate (P2.4): it must read $0$ there.
2. Pinned spectral independence of the face law needs third-order gate statistics. By D-3, reporting $\lVert\sum_iP_i\rVert$ under a few pinnings is enough.
3. Negative results go to the dictionary first (P2.1).

**To prong 3 (unlocks).** Add a row to the §4 table of [local-to-global-unlocks.md](../../local-to-global-unlocks.md):

| global object | local certificate | bound | cost |
|---|---|---|---|
| the conditional expectation onto a forget-algebra (or its sufficiency defect) | angle $c$ (expander / heredity side), or cocycle leakage $q$ (Markov side), or relative coercivity $\kappa$ (recovery side) | D-2 (sharp $(1-c)^{-1}$), D-5 ($4d_Aq$), D-6 ($\sqrt{c_*/t}\chi\kappa^{-1/2}$) | one Friedrichs angle; one modular-time integral; one Cesàro run |

Add two guards:
- the literal local gap collapses for non-commuting dynamics (D-7);
- the gap-free $\chi^2$ route keeps an exponential in the region size, so the expander-type input should be entropic (C-2) or give an exponential rate in $t$.

---

## 12. Numerical checks

All numpy, run for this note; scratch scripts (the task allowed one file), with condensed code in the Appendix. C5 imports `check_2504_fixed_point_recovery.py`; C8–C10 import `check_nc_dirichlet_lindblad.py` (both in `notes/digests/bridges/`).

| id | claim | result |
|---|---|---|
| C1 | D-1 (a)–(d) | separated pairs $c\le2.8\cdot10^{-16}$; buffer decay $0.257\to0.017$ at ratio $0.504$ (transfer ratio $0.5004$); 4-regular graph $c=\lambda/d=0.860929$; adjacent pairs $c=\max_r\lvert\mathrm{corr}\rvert\le\tanh\lvert J\rvert$; low $T$: $c=2\cdot10^{-16}$ (separated) against buffer $0.93$; $\lVert(E_AE_B)^{200}-E_{A\cup B}\rVert=6\cdot10^{-14}$ |
| C2 | D-2 sharp; Kayalar–Weinert | $\min[\lambda_{\min}-(1-c)]=-1.3\cdot10^{-15}$ over 300 pairs; $\lvert\lVert(PQ)^3-R\rVert-c^5\rvert\le2.2\cdot10^{-15}$ |
| C3 | D-3 | $\lambda_{\max}(\Psi)=\lVert\sum P_i\rVert-1$ to $10^{-12}$ (4 laws on $\{0,1\}^5$) |
| C4 | D-4, D-5 | Markov states $q\le1.1\cdot10^{-14}$; random state CMI $0.342\le3.59$; TFIM CMI $1.79\cdot10^{-4}\le5.4\cdot10^{-2}$ and $5.1\cdot10^{-3}\le0.54$; $M_\alpha\le d_A^{2\alpha}$ always |
| C5 | D-6 (b) | Davies sampler, 5-qubit rings: actual $\le$ bounds at $t=1..10^3$, ratio to the $H^{-2}$ bound $0.07$–$0.19$ |
| C6 | D-8 | $\mathcal E_{\rm str}/2=\mathrm{Var}_A$ to $10^{-14}$; Efron–Stein ratio $\le0.47$ |
| C7 | §5.6 kernels | three Fourier identities to 8 digits (Yang's to $2\cdot10^{-5}$, grid) |
| C8 | D-7 | non-commuting: gap $1.76\cdot10^{-2}$, $1.11\cdot10^{-3}$, $7.2\cdot10^{-5}$ ($n=3,4,5$), $\dim F_A=1$; commuting: $0.352$ constant, $\dim F_A=4^{n-1}/2$ |
| C9 | D-6 (c), C-1 | $\kappa=0.415,0.263,0.227$ (non-commuting), $0.644$ (commuting); recovery errors $0.009$–$0.018\le$ bounds $0.16$–$0.44$ |
| C10 | D-7, C-1, $\chi$ growth | gap a function of the distance to the far end only; $\kappa$ not decaying in $\lvert A\rvert$ (non-commuting); $\chi_{\rm KMS}$ $\times1.9$–$3.5$ per site of $A$ |
| C11 | C-2 (i) | $D(\rho_{-A}\Vert\rho)$ linear in $\lvert A\rvert$ ($\approx0.7$–$1.0$ per site) and below $2\beta\sum\lVert h_\gamma\rVert$; $\chi_{\rm KMS}$ $\times2.2$–$2.5$ (commuting) and $\times2.9$–$3.2$ (non-commuting) per site ($n=6$, $\lvert A\rvert\le4$) |

These are checks of identities and inequalities on small systems, not experiments. The evidence for C-1 is $n\le5$ only.

---

## Appendix: condensed check code

```python
import numpy as np, itertools
# ---- C1: commuting-square defect c(A,B) = ||E_A E_B - E_{AuB}||_{L2(mu)}, E_A integrates out A
def configs(n): return np.array(list(itertools.product([-1,1], repeat=n)))
def Eop(X, mu, A):
    keep=[i for i in range(X.shape[1]) if i not in A]; M=np.zeros((len(X),len(X))); G={}
    for i,k in enumerate(map(tuple, X[:,keep])): G.setdefault(k,[]).append(i)
    for g in G.values(): g=np.array(g); M[np.ix_(g,g)]=(mu[g]/mu[g].sum())[None,:]
    return M
def c(X, mu, A, B):
    D=np.sqrt(mu); T=Eop(X,mu,A)@Eop(X,mu,B)-Eop(X,mu,sorted(set(A)|set(B)))
    return np.linalg.norm(D[:,None]*T/D[None,:], 2)
n=9; X=configs(n); E=sum(X[:,i]*X[:,i+1] for i in range(n-1))+0.25*X.sum(1); mu=np.exp(0.6*E); mu/=mu.sum()
print(c(X,mu,[1],[3]), [c(X,mu,list(range(0,2+w)),list(range(2,n))) for w in range(1,6)])
# expander: X0,X1 = ends of a uniform directed edge of a d-regular graph P (adjacency, rows sum d):
#   mu2 = P.ravel()/P.sum() on pairs (i,j) with P_ij>0; c(X2,mu2,[0],[1]) == second |eigenvalue| of P/d

# ---- C2: two-projection inequality: lambda_min(2-P-Q on W-perp) == 1 - ||PQ - P_W||  (random subspaces sharing W)

# ---- C3: for mu on {0,1}^n: Psi = diag(Var)^-1 Cov - I ;  max eig(Psi) == ||sum_i v_i v_i^T|| - 1,
#          v_i = sqrt(mu)*(x_i - m_i)/||.||

# ---- C4: CMI vs leakage of u_t = (1_A x rho_BC)^{it} rho^{-it} out of M_AB (x) 1
def mfun(R,f): w,V=np.linalg.eigh((R+R.conj().T)/2); return (V*f(w))@V.conj().T
def ptr(R,dims,keep):
    k=len(dims); L='abcdefghij'; a=list(L[:k]); b=[L[k+i] if i in keep else L[i] for i in range(k)]
    d=int(np.prod([dims[i] for i in keep]))
    return np.einsum(''.join(a)+''.join(b)+'->'+''.join(a[i] for i in keep)+''.join(b[i] for i in keep),
                     R.reshape(dims+dims)).reshape(d,d)
def leakage(rho, dA, dB, dC, h=0.004, T=8):
    dims=[dA,dB,dC]; rBC=ptr(rho,dims,[1,2]); w,V=np.linalg.eigh(rho); wb,Vb=np.linalg.eigh(rBC)
    def eta(t):
        U=np.kron(np.eye(dA),(Vb*np.exp(1j*t*np.log(wb)))@Vb.conj().T)@((V*np.exp(-1j*t*np.log(w)))@V.conj().T)
        return np.linalg.norm(U-np.kron(ptr(U,dims,[0,1])/dC,np.eye(dC)),2)
    ts=np.arange(h/2,T,h); return 0.5*h*sum((eta(t)+eta(-t))/np.sinh(np.pi*t) for t in ts)
# bound (alpha=1): CMI <= 4 * sqrt(Tr rho^2 rho_{-A}^{-1}) * q  <=  4 d_A q

# ---- C8-C10: Gaussian-filter KMS sampler with single-site Paulis on A (DLL form, from check_nc_dirichlet_lindblad.py):
#   S,_,_,_ = dll_lindbladian(H, beta, [site_op(P,a,n) for a in A for P in PAULIS], lambda nu: np.exp(-(beta*nu)**2/8))
#   Wh = kron(rho^{1/4}.T, rho^{1/4});  K = -sym(Wh S Wh^{-1})   (K >= 0, KMS-symmetrised, column-stacking vec)
#   gap   = smallest nonzero eigenvalue of K
#   Pi_A  = orthogonal projector onto Wh*vec(1_A (x) E_jk);  Q = I - Pi_A;  kappa = 1/||Q K^+ Q||
#   W     = Wh^{-1} vec(rho_{-A} - rho);  chi = ||W||;   bound(t) = sqrt(0.4073/t) * chi / sqrt(kappa)
#   C5 uses davies_heis(...) from check_2504_fixed_point_recovery.py with the same symmetrisation (row-major vec).
# ---- C11: D(rho_-A||rho) = Tr rho_-A (log rho_-A - log rho)  vs  2*beta*sum ||h_gamma|| over terms meeting A;
#          chi_KMS^2 = Tr(rho_-A rho^-1/2 rho_-A rho^-1/2) - 1
```

---

## Sources

All through the digests named in the header, at the sections cited in the text; no new source was opened. *From memory*, not reopened (rule P3.3), and each either only named or also checked numerically:
- Kayalar–Weinert (1988) and Aronszajn: the alternating-projection rate (checked, C2);
- Deutsch: Friedrichs cosine $=\lVert PQ-P_W\rVert$ (consistent with C2);
- Takesaki's criterion for state-preserving conditional expectations (as used in [D-CR] §6.1 and [D-NC] §4.4);
- the Pimsner–Popa constant and Popa's commuting squares (names only; the inequalities used are proved);
- Petz's equality condition for monotonicity of relative entropy (via Jenčová–Petz, Note 1 §6);
- Witsenhausen's maximal correlation (checked, C1);
- the Diaconis–Saloff-Coste comparison theorem (name only; D-8 proves the comparison used);
- Chebyshev acceleration (§5.7, name and rate only).
