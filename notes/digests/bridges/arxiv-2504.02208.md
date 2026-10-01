# Chen–Rouzé, *Quantum Gibbs states are locally Markovian* (arXiv:2504.02208): full digest, and its relation to arXiv:2609.38007

### What is proved, how each lemma works, where every hypothesis is spent, and which bridges to expanders, local-to-global theory, NCG and the barycentre-arrow programme are theorems

*Bridges digest (prong 3 input), 2026-10-01. Programme inputs read first: [research-program.md](../../research-program.md) (rules P1.4, P3.1–P3.4, C1–C3), [local-to-global-unlocks.md](../../local-to-global-unlocks.md) (§4 common mechanism, U1–U6, §6 guards), [conditional-arrow-algebra.md](../../conditional-arrow-algebra.md) (Note 1: §3–§6, KMS states on histories, Doob walk, Jenčová–Petz sufficiency), [simplicial-complex-as-decomposition.md](../../simplicial-complex-as-decomposition.md) (Note 2: §2–§6, Lüders exclusion, frame composition), [mlp-bridge.md](../../mlp-bridge.md) (dictionary v1, E2 memory, D1–D5). Sibling digests consulted and cross-checked, not copied: [expanders.md](expanders.md), [chat-2609.38007-retrieval-status.md](chat-2609.38007-retrieval-status.md).*

**Labels.** **Source** = stated in text retrieved in this session, with its number. **Source (memory)** = a standard fact not opened in this session; check before use (P3.3). **Derivation** = proved here (with a numerical check where marked). **Interpretation**, **Analogy**, **Conjecture**, **Guard** as in the programme. Bridge hints carry THEOREM / KNOWN-LINK / ANALOGY / SPECULATION.

---

## 0. Source status and verdict

**Retrieved.**

| item | route | what was read |
|---|---|---|
| C.-F. Chen, C. Rouzé, *Quantum Gibbs states are locally Markovian*, arXiv:2504.02208**v1** (3 Apr 2025), 35 pp. | alphaXiv `get_paper_content` (fullText), 157k characters | everything: §§I–XI, App. A (Lieb–Robinson, dissipative and coherent parts), App. B (local gap), references |
| T. H. Yang, *Improved estimate of local Markovianity for quantum Gibbs states*, arXiv:2609.38007v1 (29 Sep 2026) | alphaXiv `answer_pdf_queries` (abstract, §§I–II, Table I, §VII) and `get_paper_content` fullText (95k characters) | abstract; Thms II.1–II.3; §III (Lemmas III.1–III.4, Cor. III.5) and §IV (Lemmas IV.1–IV.2, proof of Thm II.1) in full; §§V–VI proof skeleton; §VII |
| K. Kato, T. Kuwahara, *Clustering of CMI via quantum belief-propagation channels*, arXiv:2504.02235v2 | alphaXiv `answer_pdf_queries` | abstract, Conjectures 1–3, Def. 1, Assumptions 3–4, Thms 2–3, the proof step that uses Chen–Rouzé (their Lemma 8) |
| K. Liu, S. Mohanty, P. Raghavendra, A. Rajaraman, D. X. Wu, *Locally stationary distributions*, arXiv:2405.20849v4 (cited by Chen–Rouzé as [LMR+24]) | alphaXiv `answer_pdf_queries` | abstract, Def. 1.2, Thm 1.3, Lemmas 3.1–3.2, §1 discussion |

**Not retrieved.** The arXiv listing page (`arxiv.org` is egress-blocked: `EGRESS_BLOCKED`), so I could not see whether 2504.02208 has versions after v1. alphaXiv returned "Failed to fetch PDF for 2504.02208v2", and every mirror found by Exa (alphaXiv, Emergent Mind, ResearchGate, fugumt) shows v1 only. 2609.38007 cites it without a version, and its Table I reproduces the v1 bound. **Everything below is v1.** The ChatGPT conversation behind 2609.38007 was not re-attempted here; the sibling file [chat-2609.38007-retrieval-status.md](chat-2609.38007-retrieval-status.md) logs ten failed routes. **Source:** 2609.38007's AI-usage statement reads "The results presented in this paper were found by GPT-6 Astra Pro. The author verifies and organizes the proof and takes full responsibility for the correctness of the paper." So the model credited is GPT-6 Astra Pro, not Claude.

**Verdict in five lines.**
1. **Source.** For every Hamiltonian of bounded interaction degree and every $\beta$, erasing a region $A$ of the Gibbs state is undone, up to $re^{\mu|A|}t^{-\lambda}$, by running the exact-detailed-balance (KMS) Lindbladian with single-site Pauli jumps on $A$ for a random time in $[0,t]$. The map is quasi-local. Consequences: CMI $\lesssim|A||C|e^{\mu'\min(|A|,|C|)-\lambda'\mathrm{dist}(A,C)}$, and Gibbs preparation under uniform clustering without the uniform-Markov assumption.
2. **No spectral gap is used.** The engine is *time averaging* (Dirichlet form $\le2/t$) plus a nonlinear inequality "small Dirichlet form ⇒ small commutators with the jumps". That inequality is made to work at low temperature by Gaussian frequency filtering of imaginary-time conjugation.
3. **The expander-type certificate appears in exactly one place** (App. B). A uniform inverse-polynomial *local spectral gap* of the region-restricted samplers would upgrade local Markov to **global** Markov and quasi-polynomial preparation to polylog. That is the precise point where the user's intuition (Markov ↔ expansion) becomes a theorem, conditional on the gap.
4. **2609.38007 improves $e^{O(|A|)}$ to $e^{O(\text{cut strength})}=e^{O(|\partial_eA|)}$ by a static proof.** It compares with the cut Gibbs state and bounds a Petz-type sufficiency defect by the modular-time leakage of a Connes cocycle. It imports Chen–Rouzé's Lemma IX.3, gives up the explicit physical recovery map, and extends to quasi-local and power-law interactions.
5. **NCG is the common language, as identities rather than metaphors:** CMI is exactly a Petz sufficiency defect for the pair $(\rho,\rho_{-A})$ (§6.3); the volume prefactor is a Jones/Pimsner–Popa index (§6.2); the $t\to\infty$ recovery map is a modular-invariant conditional expectation (§6.1). These tie into Note 1 §6 directly.

---

## 1. Setting and definitions (Source, §§I–II of 2504.02208)

### 1.1 Hamiltonians and geometry
- $n=|\Lambda|$ qubits; $H=\sum_{\gamma\in\Gamma}H_\gamma$ with $\|H_\gamma\|\le1$ (few-body terms).
- **Interaction graph**: vertices $\Gamma$; $\gamma_1\sim\gamma_2$ iff $\mathrm{Supp}(H_{\gamma_1})\cap\mathrm{Supp}(H_{\gamma_2})\ne\emptyset$ (self-loops allowed). For $A\subset\Lambda$, $A\sim\gamma$ iff $A\cap\mathrm{Supp}(H_\gamma)\neq\emptyset$. **$d$** = maximal degree, a constant independent of $n$.
- $\mathrm{dist}(A,B)=\min\{\ell:\exists\gamma_1,\dots,\gamma_\ell,\ A\sim\gamma_1\sim\cdots\sim\gamma_\ell\sim B\}$.
- **Local patch** $H_\ell:=\sum_{\gamma:\mathrm{dist}(\gamma,A)<\ell-1}H_\gamma$.
- Remark II.2.2: $n$ does not enter the arguments; the authors believe the results formalise for infinite systems.
- $\rho_\beta=e^{-\beta H}/\mathrm{Tr}\,e^{-\beta H}$. $\tau_A$ = maximally mixed on $A$. **The erased state** is $\rho_{\beta,-A}:=\mathrm{Tr}_A[\rho_\beta]\otimes\tau_A$. $\mathsf P_A$ = non-trivial Pauli strings on $A$ ($4^{|A|}-1$ of them); $\mathsf P^1_A=\{X_i,Y_i,Z_i\}_{i\in A}$ ($3|A|$ single-site Paulis).

### 1.2 Markov notions
- $I(A:C|B)_\rho=S(\rho_{AB})+S(\rho_{BC})-S(\rho_B)-S(\rho_{ABC})$. Exact Markov: $I(A:C|B)_\rho=0\iff\exists\mathcal R_{AB}$ with $\mathcal R_{AB}[\rho_{BC}]=\rho$ (HJPW04, BP12).
- Classical Gibbs distributions of local Hamiltonians are Markov at every $\beta>0$. The pairwise, local and global Markov properties are equivalent for positive distributions (Fig. 1), but "In the quantum case, these three are not known to be equal, even allowing for approximations."
- **Conjecture [KKB20]** as quoted: for short-range $H$ on a $D$-dimensional lattice at any $\beta$ and any tripartition $\Lambda=A\sqcup B\sqcup C$, $I(A:C|B)_{\rho_\beta}\le\mathcal D(\mathrm{dist}(A,C))$ with $\mathcal D$ superpolynomially decaying. The three versions: **pairwise** ($|A|,|C|=O(1)$), **local** ($\min(|A|,|C|)=O(1)$), **global** ($|A|,|C|=O(|\Lambda|)$).
- Table I of the paper:

| property | sufficient condition | bound on $I(A:C\mid B)$, $ABC=\Lambda$ |
|---|---|---|
| pairwise | $D$-dim lattice, any $\beta$ [Kuw24] | $\exp(c\vert AC\vert -\mathrm{dist}/\xi)$ |
| **local** | **degree $d$, any $\beta$ (Thm III.1)** | $\vert A\vert \vert C\vert \exp(c\min(\vert A\vert ,\vert C\vert )-\mathrm{dist}/\xi)$ |
| global | commuting or classical, any $\beta$ | $0$ if $\mathrm{dist}\ge1$ |
| global | 1-dim, any $\beta$ [KB19, Kuw24] | $\exp(-\mathrm{dist}/\xi)$ |

### 1.3 The sampler (CKG23), the KMS inner product, detailed balance
- **Lindbladian (2.1)**: for jumps $\{A^a\}$ closed under adjoint, $\|A^a\|\le1$,
$$\mathcal L[\cdot]=-i[B,\cdot]+\sum_a\int\gamma(\omega)\Big(\hat A^a(\omega)\,\cdot\,\hat A^a(\omega)^\dagger-\tfrac12\{\hat A^a(\omega)^\dagger\hat A^a(\omega),\cdot\}\Big)d\omega .$$
- **Operator Fourier transform (2.2)**: $\hat A(\omega)=\frac1{\sqrt{2\pi}}\int e^{iHt}Ae^{-iHt}e^{-i\omega t}f(t)\,dt=\sum_{\nu\in B(H)}A_\nu\hat f(\omega-\nu)$. Here $A_\nu=\sum_{E_2-E_1=\nu}P_{E_2}AP_{E_1}$ (the Bohr components), $\hat f(\omega)=(\sigma\sqrt{2\pi})^{-1/2}e^{-\omega^2/4\sigma^2}$ and $f(t)=e^{-\sigma^2t^2}\sqrt{\sigma\sqrt{2/\pi}}$. The final bounds take $\sigma=1/\beta$.
- **Weights**: Metropolis $\gamma(\omega)=\exp(-\beta\max(\omega+\beta\sigma^2/2,0))$ (2.4); Gaussian $\gamma_G(\omega)=\exp(-(\omega+\omega_\gamma)^2/2\sigma_\gamma^2)$ with $\beta(\sigma_\gamma^2+\sigma^2)=2\omega_\gamma$ (2.5).
- **KMS inner product** $\langle X,Y\rangle_\rho=\mathrm{Tr}[X^\dagger\rho^{1/2}Y\rho^{1/2}]$ and norm $\|X\|_\rho$. Lemma II.1: $\|X\|_\rho\le\|X\|$, $\langle X,Y\rangle_\rho\le\|X\|\|Y\|$.
- **Thm II.1 [CKG23]**: (2.1) is KMS-$\rho_\beta$-detailed balanced, hence $\mathcal L[\rho_\beta]=0$ exactly.
- **Thm II.2 [CKG23, CKBG23]**: $e^{t\mathcal L}$ can be simulated to $\epsilon$ in diamond norm with $\tilde O(|P|t\beta)$ Hamiltonian-simulation time, $\tilde O(1)$ resettable ancillas, $\tilde O(|P|t)$ block-encodings of $\frac1{\sqrt{|P|}}\sum_a|a\rangle\otimes A^a$, and $\tilde O(|P|t)$ other gates. Remark II.2.1: the time average is implemented by sampling $s\sim\mathrm{Unif}[0,t]$ and running to time $s$.

### 1.4 The recovery map
$$\mathcal R_{A,t}[\cdot]:=\frac1t\int_0^t\exp(s\mathcal L_A)[\cdot]\,ds,\qquad \mathcal L_A:=\sum_{a\in\mathsf P^1_A}\mathcal L_a\qquad(3.1)$$
Each $\mathcal L_a$ is (2.1) for the single jump $A^a$ with Metropolis weight. Fig. 2 caption, verbatim: **"In fact, the recovery map is a time-averaged detailed-balanced Lindbladian based on single-Pauli jumps on A."** The recovery radius is $\approx\mathrm{Poly}(\beta)(|A|+\log(1/\epsilon))$ (Fig. 2).

---

## 2. Main results, exactly as stated (Source)

**Theorem III.1 (quasi-local recovery via time-averaged Gibbs sampling).** Let $H$ have interaction degree $\le d$, let $A\subset\Lambda$, and take $\mathcal R_{A,t}$ with jumps $\mathsf P^1_A$, Metropolis weight and $\sigma=1/\beta$, $t>0$. Then, with $\beta_0:=1/4d$,
$$\|\mathcal R_{A,t}[\rho_{\beta,-A}]-\rho_\beta\|_1\le|A|^2\,2^{2|A|}\times\begin{cases}r(\beta,d)\,t^{-\frac{128\beta_0^4}{\beta^3(\beta+5\beta_0)}}&\beta>4\beta_0,\\[2pt] r'(\beta,d)\,t^{-\frac{2\beta_0}{\beta+5\beta_0}}&\beta\le4\beta_0,\end{cases}$$
for explicit $r,r'$ (not written out in the paper). Hence there are $r,\mu>0$ and $0<\lambda<1$ depending only on $\beta,d$ with $\|\mathcal R_{A,t}[\rho_{\beta,-A}]-\rho_\beta\|_1\le re^{\mu|A|}t^{-\lambda}$ (3.2).
- Remark III.1.1: the statement still holds with $\beta_0=(1-\epsilon)/2d$. At $\beta_0=1/2d$, $\|\rho_{\beta_0}A\rho_{\beta_0}^{-1}\|$ may grow with $n$.
- Remark III.1.2: "The exponential dependence on $|A|$ is hard to remove unconditionally using the current Lindbladian approach. Right now, it appears due to the slow, inverse polynomial decay of the Dirichlet form." A faster mixing-time or gap analysis "might" give the global Markov property.

**Corollary III.1 (quasi-locality).** $\mathcal R_{A,t}$ is approximated by a strictly local $\mathcal R_{A,t,\ell}$ (the sampler built from $H_\ell$, supported within distance $\ell$ of $A$): $\|\mathcal R_{A,t,\ell}-\mathcal R_{A,t}\|_{1\to1}\lesssim|A|t(e^{-c'\ell/d\beta}+2^{-\ell})$. With $m=\min(\ln2,c'/d\beta)$ and a time $t_*(\ell)$ (from the proof, $t_*=e^{(\mu|A|+m\ell)/(1+\lambda)}$, which balances the two errors; the extracted statement prints the exponent ambiguously), $\|\rho_\beta-\mathcal R_{A,t_*,\ell}[\rho_{\beta,-A}]\|_1\lesssim r\exp(\mu'|A|-m\lambda\ell/2)$ with $\mu'=\mu+2$.

**Corollary III.2 (quantum Gibbs states are locally Markov).** For $\Lambda=ABC$ with $A$ shielded by $B$ and $\mathrm{dist}(A,C)\ge4e^2\beta d$,
$$I(A:C|B)_{\rho_\beta}\lesssim r'|A||C|\exp\big(\mu'\min(|A|,|C|)-\lambda'\,\mathrm{dist}(A,C)\big),$$
with $r',\mu',\lambda'$ depending only on $\beta,d$. Remark III.1.4: the prefactor is exponential in the smaller region and linear in the larger. The recovery statement itself does not refer to the system size; the loss comes from the $\log(\dim)$ factor in converting trace distance to entropy. The paper also states: "the current argument does not handle the case with more refined partitions ABCD [KKB20]."

**Definition III.1 (uniform clustering).** $(H,\beta)$ is uniformly clustering if for all $A,C\subseteq X\subseteq\Lambda$ with $\mathrm{dist}(A,C)\ge\ell$, $\mathrm{Cov}_{\rho^X_\beta}(A,C)\le\mathrm{Poly}(|A||C|)\|A\|\|C\|e^{-\ell/\xi}$. Here $\rho^X_\beta$ is the Gibbs state of $H_X=\sum_{Z\subseteq X}h_Z$ and $\mathrm{Cov}_\sigma(A,B)=|\mathrm{Tr}\sigma AB-\mathrm{Tr}\sigma A\,\mathrm{Tr}\sigma B|$.

**Theorem III.2 ([BK19, Thm 5], local indistinguishability).** On a $D$-dimensional lattice, uniform clustering implies, for $ABC=X\subset\Lambda$ with $\mathrm{dist}(A,C)=\ell$: $\|\mathrm{Tr}_{BC}[\rho^X]-\mathrm{Tr}_B[\rho^{AB}]\|_1\le e^{c'\beta|\partial_BC|}\big(\mathrm{Poly}(|A|,\ell^D)e^{-\ell/2\xi}+e^{-\ell/c}\big)$.

**Corollary III.3 (uniform clustering ⇒ quasi-local preparation).** For a finite-range $H$ on $[-L,L]^D$ under Def. III.1 there is a channel of dissipative gate complexity $e^{O(\log^D(n/\epsilon))}$ outputting $\rho'$ with $\|\rho'-\rho_\beta\|_1\le\epsilon$. Remark III.2.1: this improves to $\log^{O(1)}(n/\epsilon)$ if the $\mathcal L_{A,\ell}$ have gap at least inverse-polynomial in $|A|$.

**Appendix B (local gap).**
- **Def. B.1**: $\mathcal L$ is $\lambda$-locally-gapped if $-\lambda\langle X,\mathcal L^\dagger X\rangle_\rho\le\langle X,\mathcal L^{\dagger2}X\rangle_\rho$ for all $X$. It is *uniformly* locally gapped if this holds for every restricted Gibbs state $\rho^X$, $X\subset\Lambda$.
- Quoted: "The above unconditionally holds for $(H,\{A^a\})$ being a commuting Hamiltonian and local jumps $\mathsf P^1_A$ … with a local gap independent of the global system size $\lambda\ge f(|A|)>0$. For general noncommutative Hamiltonians (with local jumps), we do not know of any a priori bound on the local gap, even assuming high temperature."
- **Lemma B.2**: locally gapped ⇒ $\mathcal E(e^{\mathcal L^\dagger t}X)\le e^{-2\lambda(t-1)}\|X\|_\rho^2$ for $t>1$.
- **Cor. B.1**: a $c|A|^{-c'}$ uniform local gap gives gate complexity $O(n\,\mathrm{polylog}(n/\epsilon))$.
- **Cor. B.2**: under the same gap, $I(A:C|B)_{\rho_\beta}\le\mathrm{Poly}(|A|,|C|)\exp(-\mathrm{dist}(A,C)/\xi)$, i.e. **global Markov**.

---

## 3. Proof architecture, lemma by lemma (Source unless marked)

### 3.0 The chain in one line
$$\|\rho-\mathcal R[\rho_{-A}]\|_1\ \overset{\text{duality, fixed point}}{=}\ \sup_{\|X\|\le1}\big|\mathrm{Tr}[(\mathcal R^\dagger X-(\mathcal R^\dagger X)_{-A})\rho]\big|\ \overset{\text{twirl}}{\le}\ 2^{-2|A|-1}\sum_{S\in\mathsf P_A}\|[S,[S,\mathcal R^\dagger X]]\|_\rho\ \overset{\text{IX.5, VIII.1}}{\lesssim}\ \sum_{a\in\mathsf P^1_A}\|[A^a,\mathcal R^\dagger X]\|_\rho^{\,\nu}\ \overset{\text{X.4}}{\lesssim}\ \mathcal E_A(\mathcal R^\dagger X)^{\nu'}\ \overset{\text{VII.1}}{\le}\ (2/t)^{\nu'} .$$
In words (§IV of the paper): "the long-time averaging map $\mathcal R_{A,t}$, which keeps updating region A using all supported single qubit Pauli jumps, must mix A very thoroughly, 'conditioned' on everything else."

### 3.1 Step 0–1: duality and the Pauli twirl
- **Duality.** $\|\rho_\beta-\mathcal R_{A,t}[\rho_{\beta,-A}]\|_1=\sup_{\|X\|\le1}|\mathrm{Tr}[(\mathcal R^\dagger_{A,t}[X]-(\mathcal R^\dagger_{A,t}[X])_{-A})\rho_\beta]|$. It uses the fixed point $\mathcal R_{A,t}[\rho_\beta]=\rho_\beta$, which is where exact detailed balance (Thm II.1) is spent. So recoverability means $\mathcal R^\dagger X$ is nearly trivial on $A$.
- **Twirl identity** (§XI.A): $\rho-\rho_{-A}=2^{-2|A|-1}\sum_{S\in\mathsf P_A}[S,[S,\rho]]$. Then $\mathrm{Tr}[B[A,[A,C]]]=\mathrm{Tr}[C[A,[A,B]]]$ and Cauchy–Schwarz in the KMS inner product ($\|I\|_\rho=1$) give $\le2^{-2|A|-1}\sum_S\|[S,[S,\mathcal R^\dagger X]]\|_\rho$. Checked numerically to $10^{-17}$ (§6).

### 3.2 Regularising imaginary-time conjugation (the low-temperature engine, §IX)
The obstacle (§IV.A): naive Hölder, $\|SO\|_\rho\le\|\rho^{1/4}S\rho^{-1/4}\|\,\|O\|_\rho$, involves $\|e^{\beta H}Ae^{-\beta H}\|$, which "may generally diverge at low temperatures ∼ $e^{\Omega(n)}$ [PGPH23]". The paper calls this "a genuinely noncommutative phenomenon absent in commuting or classical Hamiltonians": a local operator can change the energy by a lot with exponentially small amplitude (Fig. 3).
- **Lemma IX.1**: $A=\frac1{\sqrt{2\sigma\sqrt{2\pi}}}\int\hat A(\omega)\,d\omega$.
- **Lemma IX.2 (key identity)**: $e^{\beta H}\hat A(\omega)e^{-\beta H}=e^{\beta\omega+\sigma^2\beta^2}\hat A(\omega+2\sigma^2\beta)$ (completing the square in the Gaussian). Hence $\|e^{\beta H}\hat A(\omega)e^{-\beta H}\|\le\frac{e^{\sigma^2\beta^2}}{\sqrt{\sigma\sqrt{2\pi}}}e^{\beta\omega}$. This depends only on $\omega$, not on $n$: "the bounds are now entirely (quasi)-local".
- **Lemma IX.3 (convergence radius; bounded degree used here)**: for single-site $\|A\|\le1$ and $|\beta|<1/2d$, $\|e^{\beta H}Ae^{-\beta H}\|\le(1-2d|\beta|)^{-1}$. Proof: $\|\mathrm{ad}_H^k(A)\|\le k!(2d)^k$, because at step $k$ at most $kd$ terms can contribute.
- **Cor. IX.1**: weight-$w$ Pauli $S$: $\|e^{\beta H}Se^{-\beta H}\|\le(1-2d|\beta|)^{-w}$. At $\beta_0=1/4d$ this is $2^w\le2^{|A|}$. **This is the first source of $e^{O(|A|)}$.**
- **Cor. IX.2**: $\|\hat A(\omega)\|\le\frac{e^{-\beta\omega+\sigma^2\beta^2}}{\sqrt{\sigma\sqrt{2\pi}}}\|e^{\beta H}Ae^{-\beta H}\|$. The proof "borrows" $e^{\pm\beta H}$ and uses that the OFT commutes with imaginary-time conjugation. So large energy changes carry exponentially small amplitude whenever conjugation at a *small* $\beta_0$ is bounded, which Lemma IX.3 guarantees.

### 3.3 The Hölder-like inequality (Lemma IX.5)
- **Lemma IX.4 (loose, high $T$)**: $\|[A,O]\|_\rho\le(\|\rho^{1/4}A\rho^{-1/4}\|+\|\rho^{-1/4}A\rho^{1/4}\|)\|O\|_\rho$.
- **Lemma IX.5**: for $\|O\|,\|A\|\le1$ and $\beta>4\beta_0>0$,
$$\|AO\|_{\rho_\beta},\|OA\|_{\rho_\beta}\lesssim\Big(\tfrac{e^{\sigma^2\beta'^2}}{\beta'\sigma}+\tfrac{e^{\sigma^2\beta_0^2}}{\beta_0\sigma}\Big)\|O\|_{\rho_\beta}^{4\beta_0/\beta}\big(\|\rho_{\beta_0}A\rho_{\beta_0}^{-1}\|+\|\rho_{\beta_0}^{-1}A\rho_{\beta_0}\|\big),\qquad\beta'=\beta/4-\beta_0 .$$
- Proof: split $A=c^{-1}(\int_{|\omega|\le\Omega}+\int_{|\omega|>\Omega})\hat A(\omega)d\omega$.
  - Low frequencies: Lemma IX.4 plus Cor. IX.2 at $\beta'$, paying $e^{\beta'\Omega}$.
  - High frequencies: drop the weighted norm (Lemma II.1), and Cor. IX.2 at $\beta_0$ gives $e^{-\beta_0\Omega}$ (eq. 9.1).
  - Balance with $e^\Omega=\|O\|_\rho^{-4/\beta}$. Both pieces become $\|O\|_\rho^{4\beta_0/\beta}$.
  - **The exponent loss $4\beta_0/\beta<1$ is the price of low temperature.** As $\beta\to4\beta_0$ it tends to $1$, recovering Lemma IX.4.

### 3.4 Global → local jumps (§VIII)
- Lemma VIII.1 is the Leibniz rule $[\prod_iA_i,O]=\sum_j(\prod_{i<j}A_i)[A_j,O](\prod_{i>j}A_i)$.
- **Cor. VIII.1**: for $S=\prod_{i=1}^wA_i$ and $\beta>1/d$, $\|[S,O]\|_\rho\lesssim2^wr'\sum_j\|[A_j,O]\|_\rho^{16\beta_0^2/\beta^2}$ (Lemma IX.5 used twice). For $\beta\le1/d$, $\lesssim2^w\sum_j\|[A_j,O]\|_\rho$.
- Remark VIII.0.1: the proof goes through without this step, but high-weight Pauli jumps carry an exponential simulation overhead that "cannot be improved by any additional mixing time assumption".

### 3.5 The Dirichlet form is a sum of squared commutators (§X)
- **Lemma X.1 ([RFA24, Lemma C.2])**: if the transition part is $\sum\alpha_{\nu_1\nu_2}A^a_{\nu_1}[\cdot]A^{a\dagger}_{\nu_2}$, put $h_{\nu_1\nu_2}=\alpha_{\nu_1\nu_2}e^{\beta(\nu_1+\nu_2)/4}$. Then $\mathcal E(X,Y):=-\langle X,\mathcal L^\dagger Y\rangle_\rho=\sum_a\sum\bar\alpha_{\nu_1\nu_2}\mathrm{Tr}[\sqrt\rho[A^a_{\nu_1},X]^\dagger\sqrt\rho[A^a_{\nu_2},Y]]$, with $\bar\alpha=h/2\cosh(\beta(\nu_1-\nu_2)/4)$.
- **Lemma X.2 (Gaussian weight)**: $\mathcal E(X,Y)=\sum_a\iint g(t)h_G(\omega)\,\mathrm{Tr}[\sqrt\rho[\hat A^a(\omega,t),X]^\dagger\sqrt\rho[\hat A^a(\omega,t),Y]]\,dt\,d\omega$.
  - Here $\hat A(\omega,t)=e^{iHt}\hat A(\omega)e^{-iHt}$, $g(t)=\frac1{\beta\cosh(2\pi t/\beta)}\ge0$ and $h_G(\omega)=e^{-\beta\omega_\gamma/4}e^{-\omega^2/2\sigma_\gamma^2}\ge0$.
  - Proof: write $1/2\cosh(\beta(\nu_1-\nu_2)/4)=\int g(t)e^{-i(\nu_1-\nu_2)t}dt$, then use Lemma IX.2 to move $e^{\pm\beta H/4}$ through.
- **Lemma X.3 (Metropolis)**: the same with $h(\omega)=e^{-\sigma^2\beta^2/8}e^{-|\omega|\beta/2}$. Proof: write Metropolis as a superposition of Gaussians, $\gamma=\int g_x\gamma_G^x\,dx$ ([CKG23, Prop. II.4]), and use $\int_0^\infty e^{-s^2-a^2/s^2}ds=\frac{\sqrt\pi}2e^{-2|a|}$.
- So $\mathcal E_a(X)=\iint g\,h\,\|[\hat A^a(\omega,t),X]\|_\rho^2$: "an elegant, manifestly PSD quadratic form of commutators". **This is the static–dynamic link.**
- **Consequence (4.1)**: $\mathcal L_a^\dagger X=0\iff\mathcal E_a(X)=0\iff[\hat A^a(\omega,t),X]=0$ a.e. $\Rightarrow[A^a,X]=0$ for all $a\in\mathsf P^1_A$. Hence $\ker\mathcal L_A^\dagger\subseteq\mathbf 1_A\otimes B(\mathcal H_{A^c})$ (§6.1 uses this).

### 3.6 Small Dirichlet form ⇒ small commutators (Lemma X.4)
For $\|A\|,\|O\|\le1$ and any $\beta,\beta_0>0$,
$$\|[A,O]\|_{\rho_\beta}\lesssim d^2|A|\Big(\tfrac{e^{\sigma^2\beta_0^2}}{\beta_0\sigma}+\tfrac{e^{\sigma^2\beta^2/16}}{\sqrt{g(1)}\beta\sigma}\Big)^{\frac{\beta+4\beta_0}{\beta+5\beta_0}}\big(\|\rho_{\beta_0}A\rho_{\beta_0}^{-1}\|+\|\rho_{\beta_0}^{-1}A\rho_{\beta_0}\|\big)^{\frac{\beta}{\beta+5\beta_0}}\mathcal E(O)^{\frac{2\beta_0}{\beta+5\beta_0}} .$$
- **Step 1 (time smearing; bounded degree used).** $t=0$ has measure zero in $\mathcal E$, so average over $|t|\le\epsilon$. The error is $\|\int_{-\epsilon}^{\epsilon}(A-A(t))dt\|\le\frac{\epsilon^3}3\|[H,[H,A]]\|\lesssim\epsilon^3d^2|A|$ (10.1); the first order cancels by symmetry.
- **Step 2 (frequency truncation).**
  - For $|\omega|\le\Omega$, Cauchy–Schwarz against the weight $g(t)h(\omega)$ gives $\lesssim\sqrt\epsilon\,e^{\sigma^2\beta^2/16}e^{\Omega\beta/4}\sqrt{\mathcal E(O)}/(\sqrt{g(\epsilon)}\beta)$ (10.2).
  - Here "the diverging reciprocal $1/h(\omega)=e^{\sigma^2\beta^2/8}e^{|\omega|\beta/2}$ is the reason why we had to truncate".
  - For $|\omega|>\Omega$, use (9.1): $\lesssim\epsilon e^{-\beta_0\Omega}(\dots)$ (10.3).
- **Step 3.** Optimise $\Omega$, then $\epsilon\le1/(d\sqrt{|A|})$.

### 3.7 Time averaging gives a gap-free decay (§VII)
- **Lemma VII.1**: $\|\mathcal L^\dagger\mathcal R_t^\dagger[O]\|\le2/t$ for $\|O\|\le1$, since $\frac1t\mathcal L^\dagger\int_0^te^{s\mathcal L^\dagger}ds=\frac1t(e^{t\mathcal L^\dagger}-I)$.
- **Cor. VII.1**: for KMS-detailed-balanced $\mathcal L$, $\mathcal E(\mathcal R^\dagger_t[O])\le2/t$.
- "Note that such a property is generally false for the Lindblad evolution $e^{\mathcal L^\dagger_At}$ itself without time-averaging, because $\mathcal L^\dagger_A$ may have arbitrarily small eigenvalues. Thus, time-averaging provides a different mechanism to obtain a small Dirichlet form that is independent of the spectral gap."
- The paper notes that the interplay of time averaging, stationarity and Dirichlet form "has also been recently exploited in full glory in the recent analysis of classical slow-mixing Markov chains [LMR+24]" (§7.1 below).

### 3.8 Assembly and exponent bookkeeping (§XI.A)
- **Low temperature** ($\beta>4\beta_0$), applied in this order:
  1. Lemma IX.5 with Cor. IX.1 ($\|\rho_{\beta_0}S\rho_{\beta_0}^{-1}\|\le2^{|A|}$) peels the outer commutator: exponent $4\beta_0/\beta$.
  2. Cor. VIII.1: exponent $16\beta_0^2/\beta^2$.
  3. Lemma X.4: exponent $2\beta_0/(\beta+5\beta_0)$.
  4. Jensen, $\mathbb E|x|^\alpha\le(\mathbb E|x|)^\alpha$, restores the $1/3|A|$ normalisation.
  5. Cor. VII.1.
- Product of exponents: $\frac{4\beta_0}\beta\cdot\frac{16\beta_0^2}{\beta^2}\cdot\frac{2\beta_0}{\beta+5\beta_0}=\frac{128\beta_0^4}{\beta^3(\beta+5\beta_0)}$, as in Thm III.1.
- **High temperature** ($\beta\le4\beta_0$): Lemma IX.4, then Cor. VIII.1 without exponent loss, then Lemma X.4: exponent $2\beta_0/(\beta+5\beta_0)$.
- Collected prefactor: $|A|^2 2^{2|A|}$.

### 3.9 Quasi-locality (§VII.A and App. A)
- **Lemma VII.2 (Duhamel)**: $\|\mathcal R^\dagger_{A,t,\ell}-\mathcal R^\dagger_{A,t}\|_{\infty\to\infty}\lesssim t\|\mathcal L^\dagger_{A,\ell}-\mathcal L^\dagger_A\|_{\infty\to\infty}$.
- **Lemma VII.3**: for $\ell\ge4e^2\beta d$, $\|\mathcal L^\dagger_{A,\ell}-\mathcal L^\dagger_A\|_{\infty\to\infty}\lesssim|A|(e^{-c'\ell/d\beta}+2^{-\ell})$. Remark VII.0.1: the tail falls exponentially in $\ell$ while the error accumulates linearly in $t$, "Thus, the quasi-locality holds for exponential times."
- **Lemma A.1 (Lieb–Robinson; bounded degree via path counting $|A|d(d-1)^{p-1}$)**: $\|e^{iH_\ell t}Ae^{-iH_\ell t}-e^{iHt}Ae^{-iHt}\|\lesssim\|A\||A|\frac{(2d|t|)^\ell}{\ell!}$.
- **Lemma A.2 (dissipative part)**: uses the "purified jump" $V=\int\hat A_H(\omega)\otimes|\omega\rangle d\omega$, which gives $\|\mathcal D-\mathcal D_\ell\|_\diamond\le4\|V-V'\|$. The difference is bounded by the Gaussian time tail of $f$ plus Lieb–Robinson, with $T=\ell/2de^2$. As extracted: $\lesssim\|A\|(e^{-\sigma^2\ell^2/2d^2}+|A|e^{-\ell})$.
- **Cor. A.1 / A.2 (coherent part $B_a$)**:
  - Norm bound independent of $n$: the $1/t$ singularity of the kernel $b_2$ is tamed by symmetrising the integral.
  - $\|B_a-B^\ell_a\|\lesssim e^{-c'\ell/d\beta}+|A|2^{-\ell}$ for $\ell\ge4e^2d\beta$, with $T=\ell/4e^2d\beta-1$.

### 3.10 Corollaries
- **Cor. III.1.** $2\Delta\lesssim re^{\mu|A|}t^{-\lambda}+t|A|(e^{-c'\ell/d\beta}+2^{-\ell})$; choose $t$ to balance the two terms.
- **Cor. III.2.** Take $\ell=\mathrm{dist}(A,C)-1$, so $H_\ell$ lives in $B$ and $\mathcal R_{A,t,\ell}$ acts on $AB$. Then $I(A:C|B)\le\Delta\log\dim C+h_2(\Delta)\lesssim\log(\dim C)\sqrt\Delta$ ([Wil11, Thm 11.10.5]).
  - *Own reading:* the $\min(|A|,|C|)$ comes from applying the same argument with $A$ and $C$ exchanged, using the symmetry of CMI; the linear factor is $\log\dim$ of the other region.
- **Cor. III.3.** Patching after [BK19] on a coloured tiling $A^{h,j}_-\subset A^{h,j}\subset A^{h,j}_+$, $h=1..h_0$ colours, with same-colour patches disjoint.
  - Apply the trace-out-and-recover maps $\mathcal F_{A^{h,j}}=\mathcal R_{A^{h,j},t,\ell}\circ(\tau\otimes\mathrm{Tr}_{A^{h,j}})$ colour by colour.
  - Local indistinguishability (Thm III.2) handles the punched holes; telescoping sums the errors.
  - Side length $\log n$ forces $|A|\sim\log^Dn$, hence $\ell=\Theta(\log^D(n/\epsilon))$ and $t_*=e^{\Theta(\log^D(n/\epsilon))}$, giving quasi-polynomial cost. *(Own reading: the $e^{\mu|A|}$ of Thm III.1 forces $\ell\gtrsim|A|$, and $t_*\sim e^{|A|}$.)*

### 3.11 Where each hypothesis is spent

| hypothesis | spent in |
|---|---|
| exact KMS detailed balance (CKG23) | fixed point in the duality step; $\mathcal E$ as a KMS quadratic form (Cor. VII.1); the explicit forms X.1–X.3 |
| bounded degree $d$, $\Vert H_\gamma\Vert \le1$ | Lemma IX.3 (radius $1/2d$; it sets $\beta_0=1/4d$ and the $2^{\vert A\vert }$); Lemma X.4 step 1 ($d^2\vert A\vert $); Lemma A.1 path counting; $\ell\ge4e^2\beta d$ |
| Gaussian filter width $\sigma$ ($=1/\beta$) | Lemma IX.2 (Gaussian beats exponential); Lemma A.2 time tail |
| Metropolis weight | $h(\omega)\propto e^{-\vert \omega\vert \beta/2}$ (Lemma X.3), whose reciprocal forces the cut-off $\Omega$ |
| jumps $=\mathsf P^1_A$ | the twirl reduces recovery to commutators with $\mathsf P_A$; Cor. VIII.1 reduces to single sites; (4.1) puts the kernel inside $\mathbf 1_A\otimes B(\mathcal H_{A^c})$ |
| time averaging | Lemma VII.1: $2/t$ with no gap |
| full tripartition $ABC=\Lambda$, $A$ shielded | Cor. III.2: the recovery map must live on $AB$ and the state must be the global Gibbs state |
| lattice, finite range, uniform clustering | Cor. III.3 only |
| local gap | App. B only (conditional corollaries) |

### 3.12 Constants traced (Derivation, from the stated formulas)
- With $\beta_0=1/4d$ and $\beta\gg1/d$: $\lambda=\frac{128\beta_0^4}{\beta^3(\beta+5\beta_0)}\approx\frac1{2d^4\beta^4}$ and $m\approx c'/d\beta$.
- So the decay rate in Cor. III.2, $\approx m\lambda/4$, is $\approx c'/(8d^5\beta^5)$: **a CMI correlation length growing like $d^5\beta^5$.**
- The recovery time $t_*$ is exponential in $|A|/\lambda\sim d^4\beta^4|A|$.
- For comparison, 2609.38007's explicit decay length (4.23), $\xi_\beta=\frac{1+\alpha_0}{2\alpha_0}\frac{\pi+\beta v}{\pi\mu_{LR}}$ with $\alpha_0=\min(1,1/4\beta dJ)$, grows like $\frac{2dJv}{\pi\mu_{LR}}\beta^2$. *(Both are upper bounds; neither paper claims optimality.)*

---

## 4. What is genuinely new

**Per the authors (Source, abstract and §V).**
1. Local Markov at **every** $\beta$ for bounded-degree interaction graphs (any dimension, any graph). This settles an open problem of [Kuw24].
2. The recovery map is a **thermalisation dynamics**, explicit and physical, depending mostly on the nearby Hamiltonian.
3. Two tools are announced: a regularisation scheme for imaginary-time-evolved operators at arbitrarily low temperature; and "a connection between the Dirichlet form, a dynamic quantity, and the commutator in the KMS inner product, a static quantity."
4. The uniform-Markov assumption of [BK19] is removed from clustering-based preparation.

**Assessment (Interpretation).**
- The genuinely new *mechanism* is that **no mixing is needed for a Markov statement**. Time averaging alone makes the Dirichlet form $O(1/t)$, and a nonlinear (Hölder-type) inequality converts "small Dirichlet form" into "nearly in the commutant of the jumps". That is a *local Poincaré inequality without a gap*, at the price of a power $\lambda<1$ and a constant $e^{\mu|A|}$.
- The genuinely new *technique* is the Bohr-frequency (OFT) regularisation. Analyticity of imaginary-time conjugation at small $\beta_0<1/2d$ is turned into a bound valid at all $\beta$, using that a Gaussian filter in frequency beats any exponential $e^{\beta\omega}$.
- Imported: Lieb–Robinson, entropic continuity, [BK19] patching, and [RFA24]'s Dirichlet-form expression.
- The losses are visible and the paper names them: $e^{\mu|A|}$ (Remark III.1.2) and $\mathrm{poly}(\beta)$ length scales. ABCD partitions are not handled.

---

## 5. Relation to arXiv:2609.38007 (T. H. Yang)

### 5.1 What 2609.38007 proves (Source)
**Setting.** A full partition $\Lambda=A\sqcup B\sqcup C$, "essential to the cut-Hamiltonian argument". Write $r=d(A,C)$. The cut interaction is $V=\sum_{X\cap A\ne\emptyset\ne X\cap A^c}\Phi(X)$, with **cut strength** $g_{A|A^c}=\sum_{X\text{ crossing}}\|\Phi(X)\|_\infty$.

- **Thm II.1 (finite range).** Assume $\|h_\gamma\|\le J$, $\mathrm{diam}\le R_0$, and each term overlaps at most $d$ terms (counting itself). For every $\beta>0$: $I(A:C|B)\le C_\beta\exp(C_\beta g_{A|A^c}-c_\beta r)$.
  - On $\mathbb Z^D$, $g_{A|A^c}\le C|\partial_eA|$, so $I\le C_\beta e^{C_\beta|\partial_eA|-c_\beta r}$.
  - The abstract says this "improves the exponential volume dependence of the local Markov bounds established by Chen and Rouzé [1]".
- **Thm II.2 (stretched-exponential, unrestricted body size, $D\ge2$).** Assume $\sum_{X\ni x,y}\|\Phi(X)\|\le J_{\mu,\theta}e^{-\mu d(x,y)^\theta}$, $0<\theta\le1$. Then $I\le C_\beta\exp[C_\beta(1+|\partial_eA|)^{\theta/D}-c_\beta(1+r)^{\theta/(D+1-\theta)}]$.
- **Thm II.3 (power law, $\zeta>D$, unrestricted body size).**
  - For $D<\zeta<2D$: $I\le C_\beta(1+|A|)^{b_\zeta+1}[\log(e+r)]^{-b_\zeta}$.
  - An endpoint version holds at $\zeta=2D$.
  - For $\zeta>2D$: $I\le C_\beta(1+|A|)^{b_\zeta+1}(1+r)^{-\nu_\zeta}$, where $b_\zeta=1+2\zeta/D$ and $\nu_\zeta=b_\zeta(\zeta-2D)/(\zeta-D+b_\zeta)$.
- Abstract: "Our proof is static, in contrast to the dynamical approach of Chen and Rouzé."

### 5.2 How (Source, §§III–IV)
- **Lemma III.1 (cut comparison).** Let $\sigma=\gamma_A\otimes\gamma_{BC}$ (Gibbs states of $H_A$ and $H_{BC}$) and $\delta_R(\rho,\sigma)=D(\rho\|\sigma)-D(\rho_R\|\sigma_R)$. Then:
  - $\delta_{AB}(\rho,\sigma)=I(A:C|B)_\rho+D(\rho_{BC}\|\gamma_{BC})-D(\rho_B\|\gamma_B)$, so $I(A:C|B)\le\delta_{AB}(\rho,\sigma)$ by data processing;
  - $|\log Z-\log Z_0|\le\beta g_{A|A^c}$.
- **Lemma III.2.** A resolvent representation $\delta_R=\int_0^\infty j(s)ds$, with $\Delta=L_\sigma R_{\rho^{-1}}$ the relative modular operator.
- **Lemma III.3.** An exact variational formula for $j(s)$ over $\mathrm{Ran}\,W_R$, where the isometry satisfies $W_R^\dagger\Delta W_R=\Delta_R$.
- **Lemma III.4.** $\delta_R\le\langle\rho^{1/2},\log(1+e^{-\ell-K})\rho^{1/2}\rangle+e^{-\ell}+\mathfrak q_R^2(e^\ell+2\ell)$.
  - The **modular leakage** is $\eta_R(t)=\|\sigma^{it}\rho^{-it}-\mathsf E_R(\sigma^{it}\rho^{-it})\|_\infty$, with $\mathsf E_R$ the tracial conditional expectation.
  - The integrated leakage is $\mathfrak q_R=\frac12\int\eta_R(t)/|\sinh\pi t|\,dt$.
- **Cor. III.5.** $\delta_R\le(\frac1\alpha+3)[\mathrm{Tr}\rho^{1+\alpha}\sigma^{-\alpha}]^{1/(1+\alpha)}\mathfrak q_R^{2\alpha/(1+\alpha)}$.
- **Lemma IV.1.** $\mathfrak q_{AB}\le C_0(1+\beta g)\exp[-\lambda_\beta(r-R_0)_+]$ with $\lambda_\beta=\pi\mu_{LR}/(\pi+\beta v)$.
  - The modular cocycle $\sigma^{it}\rho^{-it}=e^{-i\beta tH_0}e^{i\beta t(H_0+V)}$ is the interaction-picture evolution of the cut $V$.
  - Lieb–Robinson keeps it near $AB$ for $|t|\lesssim r/\beta v$; the kernel $1/\sinh\pi t$ suppresses longer modular times.
- **Lemma IV.2.** $\mathrm{Tr}(\rho^{1+\alpha}\sigma^{-\alpha})\le\exp[\alpha\beta g-\frac g{2dJ}\log(1-2\alpha\beta dJ)]$, and in particular $\le e^{3\alpha\beta g}$ for $\alpha\le\alpha_0=\min(1,1/4\beta dJ)$.
  - The proof says verbatim: "We use the nested-commutator counting argument of [1, Lemma IX.3 and its proof], applied here to an original interaction term $h_\gamma$ with the interaction norm bound J restored."
  - This gives $\|e^{-uH_0}Ve^{uH_0}\|\le g/(1-2dJ|u|)$, and Grönwall finishes.

### 5.3 The relation, stated precisely

| | Chen–Rouzé 2504.02208 | Yang 2609.38007 |
|---|---|---|
| object controlled | trace distance $\Vert \mathcal R_{A,t}[\rho_{-A}]-\rho\Vert _1$ (recovery) | $I(A:C\mid B)\le\delta_{AB}(\rho,\gamma_A\otimes\gamma_{BC})$ (relative-entropy loss) |
| logical direction | explicit recovery ⇒ CMI (continuity, loses $\log\dim C$) | CMI ⇒ recovery maps exist "through the recovery theorem" (Fawzi–Renner); none is constructed |
| reference state | the erased state $\rho_{-A}=\tau_A\otimes\rho_{BC}$ | the cut Gibbs state $\gamma_A\otimes\gamma_{BC}$ |
| what carries the decay | Dirichlet-form relaxation in Lindblad time, then Lieb–Robinson for the sampler | Lieb–Robinson for the **modular cocycle** in modular time, kernel $1/\lvert\sinh\pi t\rvert$ |
| imaginary-time input | Lemma IX.3 for Pauli strings on $A$ ⇒ $2^{w}\le2^{\lvert A\rvert}$ (volume) | the **same Lemma IX.3** for the cut $V$ ⇒ $e^{3\alpha\beta g}$ (boundary) |
| low-$T$ device | Gaussian OFT filter, frequency truncation | Rényi moment at small $\alpha\le1/4\beta dJ$ plus resolvent cut-offs |
| prefactor | $\lvert A\rvert\lvert C\rvert e^{\mu'\min(\lvert A\rvert,\lvert C\rvert)}$ | $C_\beta e^{C_\beta g_{A\lvert A^c}}$, independent of $\lvert C\rvert$ and $n$ |
| decay length (traced, §3.12) | $\sim d^5\beta^5$ | $\sim\beta^2$ |
| interactions | bounded degree, any graph | finite range (any metric site set); also stretched-exponential and power law on $\mathbb Z^D$ |
| extras | quasi-polynomial Gibbs preparation (Cor. III.3); global Markov under a local gap (App. B) | coarse-graining in $D=2$ (7.2); polynomial-in-$\lvert A\rvert$ trade-off (7.3) |

**Derivations (own, from the two statements).**
- *Symmetrised 2609 bound.* CMI is symmetric and swapping $A$ with $C$ keeps a full partition. So $I(A:C|B)\le C_\beta\exp(C_\beta\min(g_{A|A^c},g_{C|C^c})-c_\beta r)$. With $g\le dJ|A|$ this implies CR's Cor. III.2 form up to constants. That holds for finite-range $H$, and for CR's bounded-degree setting once the sites carry the metric induced by the interaction graph and the Lieb–Robinson constants of 2609's (4.8) are obtained by path counting (own remark, also made by the sibling digest §7.3). **So, for full partitions, 2609 Thm II.1 implies CR's CMI corollary.** It does not contain CR's Thm III.1, Cor. III.1, Cor. III.3 or App. B, which are statements about an explicit map and its cost.
- *Where volume becomes boundary.* In CR the $e^{O(|A|)}$ enters at Cor. IX.1: conjugating weight-$|A|$ Pauli strings at $\beta_0$. It is then amplified by $t_*\sim e^{\mu|A|/\lambda}$. In 2609 only the operator that *differs* between $\rho$ and the reference is conjugated, namely the cut $V$, whose norm is the cut strength. In other words, 2609 conjugates the *interface*, CR the *region*.
- *Index versus interface divergence (§6.2).* CR's reference $\rho_{-A}$ satisfies $D_{\max}(\rho\|\rho_{-A})\le2|A|\log2$, the log of a Jones index. 2609's reference satisfies $D_{1+\alpha}(\rho\|\gamma_A\otimes\gamma_{BC}):=\frac1\alpha\log\mathrm{Tr}\rho^{1+\alpha}\sigma^{-\alpha}\le3\beta g_{A|A^c}$ (from (4.15)). **Changing the reference state from "erase A" to "cut the bonds" is exactly what turns the volume into the boundary.**
- *What 2609 does not give.* An explicit, physical, quasi-local recovery channel. 2609 gets recovery maps only "through the recovery theorem": Fawzi–Renner, which in 2609's words bounds "the recovery error in terms of the CMI". That yields *existence* of a map on $AB$, not a construction.
  - The converse direction, recovery error ⇒ CMI, is what CR use via [Wil11, Thm 11.10.5]. Kato–Kuwahara's Lemma 1 states a version of it with constant $7\log_2\min(D_A,D_C)$.
  - CR's map is the only one of the two that is a thermalisation dynamics implementable by Thm II.2.

**Convergence check (rule C3).** The sibling [chat-2609.38007-retrieval-status.md](chat-2609.38007-retrieval-status.md) reached two of these statements by a different route, from the 2609 text and its own checks:
- §7.1: Cor. III.5 holds for any subalgebra with a trace-preserving conditional expectation;
- §7.3: the boundary gain is isoperimetric and disappears on expanders.

I re-read the 2609 proof and agree on both. The only tensor-product inputs are the isometry $W_R$ in (3.11) and $\mathsf E_R$ in (3.13), both of which have subalgebra analogues. The cut functional $g$ is the numerator of a Cheeger ratio.

### 5.4 Other work that uses Chen–Rouzé as a black box (Source)
- **Kato–Kuwahara 2504.02235v2** (concurrent; CR call it a "stronger global Markov property, but only at high temperatures").
  - Their Lemma 8 is CR's Cor. III.2 verbatim ("Corollary III.2 in Ref. [67]").
  - It is used inside their belief-propagation-channel proof of Theorem 3: under uniform clustering there is an approximate BP channel with error $e^{\Theta(\beta)-\Theta(1)\kappa_\beta(r/\xi_\beta)^{1/D}}+\Theta(n)e^{-\Theta(r)/\tilde\xi_\beta}$.
  - The local Markov property enters exactly where a small region $|A|\propto\ell^D$ has to be recovered.
- **Rosa-Ruiz–Scandi–Capel–Alhambra** (2606.28054, as described in 2609's text, not opened): extend CR's local Markov property to $k$-local power-law Gibbs states with $\alpha>D$.

---

## 6. Own derivations, with a numerical check

Script: [`check_2504_fixed_point_recovery.py`](check_2504_fixed_point_recovery.py) (numpy). It uses a 5-qubit ring with $A=\{0\}$ and $\beta=1.5$, comparing a classical Ising ring ($ZZ+0.3Z$) with a transverse-field Ising ring ($ZZ+0.9X+0.3Z$). The generator is a Davies generator with Metropolis rates and single-site Pauli jumps on $A$. It is KMS-symmetric, which is all that §6.1 uses; CR's CKG23 sampler is the quasi-local member of the same class.

### 6.1 The $t\to\infty$ limit of CR's map: exact recovery, a conditional expectation, and why the theorem is a finite-$t$ trade-off
**Derivation.** Let $\mathcal L_A$ be KMS-symmetric for a faithful $\rho$, with fixed-point set $\mathcal F_A=\ker\mathcal L_A^\dagger$. Suppose $\mathcal F_A\subseteq\mathcal N_A:=\mathbf 1_A\otimes B(\mathcal H_{A^c})$, which CR's (4.1) gives for jumps $\mathsf P^1_A$. Let $\mathcal P_A=\lim_{t\to\infty}\mathcal R_{A,t}$.
- $\mathcal P_A^\dagger$ is the KMS-orthogonal projection onto $\mathcal F_A$. For any $X$, $\mathcal P_A^\dagger X\in\mathcal N_A$, so $(\mathcal P^\dagger_AX)_{-A}=\mathcal P^\dagger_AX$.
- Since $\rho$ and $\rho_{-A}$ agree on $\mathcal N_A$, $\mathrm{Tr}[X\mathcal P_A(\rho_{-A})]=\mathrm{Tr}[(\mathcal P_A^\dagger X)\rho_{-A}]=\mathrm{Tr}[(\mathcal P_A^\dagger X)\rho]=\mathrm{Tr}[X\rho]$.
- **So $\mathcal P_A[\rho_{-A}]=\rho$ exactly.** It holds for every $\omega$ with $\mathrm{Tr}_A\omega=\mathrm{Tr}_A\rho$.

Further structure (Source (memory): Tomiyama, a norm-one projection onto a C\*-subalgebra is a conditional expectation; Takesaki, a $\rho$-preserving conditional expectation onto $\mathcal M_0$ exists iff $\sigma^\rho_t(\mathcal M_0)=\mathcal M_0$):
- $\mathcal P_A^\dagger$ is the $\rho$-preserving conditional expectation onto the algebra $\mathcal F_A$.
- Hence $\mathcal F_A\subseteq\mathcal N_A^\sigma:=\{X\in\mathcal N_A:\sigma^\rho_t(X)\in\mathcal N_A\ \forall t\}$, the largest modular-invariant subalgebra of $\mathcal N_A$. Here $\sigma^\rho_t=\mathrm{Ad}\,e^{-i\beta tH}$, so $\mathcal N_A^\sigma$ is the largest $\mathrm{ad}_H$-invariant subspace of $\mathcal N_A$.

**Check (all to machine precision unless noted).**

| quantity | classical Ising | transverse-field Ising |
|---|---|---|
| twirl identity (§3.1) | $1.4\cdot10^{-17}$ | $1.4\cdot10^{-17}$ |
| KMS symmetry defect | $1.3\cdot10^{-15}$ | $6.9\cdot10^{-10}$ (frequency-binning tolerance) |
| $\dim\mathcal F_A$ / $\dim\mathcal N_A^\sigma$ / $\dim\mathcal N_A$ | 96 / 96 / 256 | **2 / 2** / 256 |
| $\mathcal F_A\subseteq\mathcal N_A$ | yes ($8\cdot10^{-14}$) | yes ($9\cdot10^{-10}$) |
| $\lVert\mathcal P_A[\rho_{-A}]-\rho\rVert_1$ | $7.7\cdot10^{-16}$ | $1.1\cdot10^{-13}$ |
| $\mathcal P_A^\dagger(X_{\text{far}})$ | $=X_{\text{far}}$ ($1.6\cdot10^{-14}$): **local** | $\ne X_{\text{far}}$ (error 1.0): **non-local** |
| $\lVert\mathcal R_{A,t}[\rho_{-A}]-\rho\rVert_1$ at $t=1,10,10^2,10^3$ | 0.313, 0.036, 0.004, 0.000 | 0.668, 0.173, 0.018, 0.002 |

- **Classical case.** $\mathcal F_A=(M_1\oplus M_2\oplus M_1)_{\{1,4\}}\otimes B(\mathcal H_{\{2,3\}})$, i.e. the operators off $A$ that commute with the boundary field $Z_1+Z_4$ that $A$ feels. Dimension check: $6\cdot16=96$. **Own reading:** this is the operator form of the classical heat bath depending only on $\partial A$.
- **Non-commuting case.** $\mathcal F_A=\mathrm{span}\{\mathbf 1,R\}$, with $R$ the ring reflection fixing site 0 ($\|[H,R]\|=4\cdot10^{-16}$, $\|\mathcal L_A^\dagger R\|=5\cdot10^{-15}$). So $\mathcal P_A$ is, up to a symmetry, the global replacement channel.

**Interpretation.** Recovery from $\rho_{-A}$ is free if locality is ignored (the $t=\infty$ map recovers exactly), so the entire content of Thm III.1 is the *finite-$t$ trade-off*: Dirichlet relaxation ($t$ large) against Lieb–Robinson spreading ($t$ small). For commuting $H$ the infinite-time map is already local, which is the exact Markov property (Brown–Poulin). For non-commuting $H$ the modular flow pushes $\mathcal N_A$ out of itself ($\mathcal N_A^\sigma$ collapses), and only a quantitative truncation survives. The NCG reading: **exact local Markov is modular invariance of the far subalgebra; approximate local Markov is approximate modular invariance for short modular times.** 2609.38007 uses literally that: Lieb–Robinson for $\sigma^{it}\rho^{-it}$, weighted by $1/\sinh\pi t$.

### 6.2 The volume prefactor is a Jones index
**Derivation.** $\rho_{-A}=4^{-|A|}\sum_{S\in\{\mathbf 1\}\cup\mathsf P_A}S\rho S\ge4^{-|A|}\rho$, so $\rho\le4^{|A|}\rho_{-A}$ and $D_{\max}(\rho\|\rho_{-A})\le2|A|\log2$. Checked: the minimum eigenvalue of $4^{|A|}\rho_{-A}-\rho$ is $\ge0$ in both models.

$4^{|A|}=\dim(\mathcal H_A)^2$ is the Jones index of $\mathcal N_A\subset B(\mathcal H)$, and $4^{-|A|}$ is its Pimsner–Popa constant (names: Source (memory)). CR's $2^{2|A|}$ enters through exactly this twirl representation of $\mathsf E_{\mathcal N_A}$ (§3.1). 2609 replaces this by an interface quantity, $D_{1+\alpha}(\rho\|\gamma_A\otimes\gamma_{BC})\le3\beta g$ (§5.3).

### 6.3 CMI is a Petz sufficiency defect for CR's own pair
**Derivation.**
- $D(\rho\|\tau_A\otimes\rho_{BC})=S(\rho_{BC})-S(\rho)+|A|\log2$ and $D(\rho_{AB}\|\tau_A\otimes\rho_B)=S(\rho_B)-S(\rho_{AB})+|A|\log2$. Subtracting,
$$I(A:C|B)_\rho=\delta_{AB}(\rho,\rho_{-A})=D(\rho\|\rho_{-A})-D(\rho_{AB}\|(\rho_{-A})_{AB}).$$
- So **CMI is the loss of distinguishability between the Gibbs state and its $A$-erased version under restriction to the subalgebra $B(\mathcal H_{AB})\otimes\mathbf 1_C$**: the Petz sufficiency defect of that subalgebra for the pair $\{\rho,\rho_{-A}\}$.
- This is 2609's Lemma III.1 with $\sigma=\rho_{-A}$, where the correction $D(\rho_{BC}\|\sigma_{BC})-D(\rho_B\|\sigma_B)$ vanishes.
- Checked: transverse-field Ising, $B=\{1,4\}$, $C=\{2,3\}$: both sides $1.679819\cdot10^{-2}$, difference $1.6\cdot10^{-15}$.
- **Consequence:** CR (recover $\rho$ from $\rho_{-A}$) and 2609 (bound the sufficiency defect with a better-behaved second state) are two treatments of the *same pair of states and the same subalgebra*. This is exactly the form of Note 1 §6 (Jenčová–Petz).

---

## 7. Bridges

### 7.1 Expander graphs
- **THEOREM (CR App. B, conditional).** A **uniform inverse-polynomial local spectral gap** of the region-restricted KMS samplers, holding for all regions and all restricted Gibbs states, implies **global Markov**, $I\le\mathrm{Poly}(|A|,|C|)e^{-\mathrm{dist}/\xi}$ (Cor. B.2), and $O(n\,\mathrm{polylog})$ preparation (Cor. B.1).
  - This is the one place in the paper where an expander-type certificate drives a Markov statement.
  - Without the gap, time averaging yields only local Markov, and the $e^{\mu|A|}$ is attributed to "the possibility of an exponentially long mixing time" (§V).
  - Uniformity over restricted states is the *hereditary* form of the certificate, in the shape of "every link is a spectral expander" (local-to-global-unlocks §2, §4).
- **KNOWN-LINK (CR cite it; text retrieved).** Liu–Mohanty–Raghavendra–Rajaraman–Wu [LMR+24 = 2405.20849] come from the expander and spectral-independence community.
  - Thm 1.3: for a reversible chain started anywhere, a uniformly random time $t\sim[0,T]$ with $T=\frac1{\delta\varepsilon}\log\frac1{\pi_{\min}}$ gives an $\varepsilon$-locally-stationary law with probability $\ge1-\delta$.
  - Lemma 3.1: $\mathbb E_{t\sim[0,T]}\mathcal E(f_t,\log f_t)\le\mathrm{KL}(\nu\|\pi)/T$.
  - Locally stationary laws have "conditional marginals matching" the stationary law (their §1, Lemma 3.5).
  - CR's Lemma VII.1 and Cor. VII.1 are the Heisenberg-picture, quadratic-form, quantum version ($\mathcal E(\mathcal R_t^\dagger O)\le2/t$). The conclusion "the conditional structure on $A$ is right" is CR's recovery of $A$ given the rest.
- **ANALOGY (precise).** Expanders: the gap is the Poincaré constant, $\mathrm{Var}(f)\le\mathrm{gap}^{-1}\mathcal E(f)$, with exponent 1. CR's chain proves a gap-free local Poincaré inequality, $\|X-\mathsf E_{\mathcal N_A}X\|_{\rho}\le2^{-2|A|-1}\sum_{S}\|[S,[S,X]]\|_\rho\lesssim r|A|^22^{2|A|}\mathcal E_A(X)^{\lambda}$ for all $\|X\|\le1$, with $\lambda<1$. (Derivation: Lemmas IX.5, VIII.1 and X.4 hold for any such $X$, not only $\mathcal R^\dagger X$. The first step is the Heisenberg form of the twirl identity, with $\mathsf E_{\mathcal N_A}(X)=4^{-|A|}\sum_{S}SXS$.) The local gap of Def. B.1 is the exponent-1 version, relative to $\mathcal F_A$. "Time average plus Hölder-Poincaré" is the slow-mixing substitute for "iterate plus gap".
- **Derivation (degeneration on expander interaction graphs).** CR's Thm III.1 and Cors III.1–2 use only bounded degree, so they apply on expander interaction graphs. But there:
  - $\mathrm{dist}(A,C)\le\mathrm{diam}=O(\log n)$;
  - the recovery radius $\ell\sim(\mu|A|+\log1/\epsilon)/m\lambda$ covers $\sim|A|d^\ell$ sites, i.e. the whole system once $\ell\gtrsim\log_dn$;
  - 2609's boundary gain vanishes, since $g_{A|A^c}\gtrsim h|A|$ under edge expansion $h$ (with term norms bounded below).

  Expansion is the *enemy* of Markov-type locality (no room for shields) and the *friend* of mixing. This agrees with the sibling's §7.3.

### 7.2 Markov random fields, Hammersley–Clifford, and the heredity half of local-to-global
- **KNOWN-LINK (CR §I.B, quoted).** "(Exact Markov) + (Strong spatial mixing) ⟹ MCMC mixes in quasi-linear time" [Mar99]. CR position Thm III.1 as the quantum *Markov half*, and Def. III.1 as the clustering half.
- **Interpretation.** Read against the programme's §4 table: the Markov property is what makes the class **hereditary under pinning**. Conditioning on the outside of $A$ gives a Gibbs measure of the same kind on $A$, with the boundary as a field (the $\mathcal F_A$ of the classical check is literally the commutant of that boundary field). Spectral independence, strong spatial mixing, or a uniform local gap are the **local certificate**, and the gap is the **global bound**.
- So Hammersley–Clifford and expansion are not one essence. They are the two complementary halves that every local-to-global sampling theorem needs. In the quantum setting, CR supply the first half *approximately*, and App. B shows the second half closing the loop.
- **Guard (Source, Kato–Kuwahara v2 §I).** For $A\cup B\cup C\subsetneq\Lambda$ (marginals with a traced-out exterior), CMI need not decay at low temperature (topological order: Castelnovo–Chamon, Hastings, their [44, 45]). Both CR and 2609 require full partitions. *Marginalisation can create memory.*

### 7.3 High-dimensional expanders: CR's patch map is a quantum down-up step
- **ANALOGY (precise).** $\mathcal F_A=\mathcal R_{A,t}\circ(\tau_A\otimes\mathrm{Tr}_A)$, the paper's "trace-out-and-recovery" map, has two steps:
  - the **down** step $\tau_A\otimes\mathrm{Tr}_A=\mathsf E_{\mathcal N_A}^{\,*}$ forgets $A$ (restriction to the subalgebra off $A$);
  - the **up** step re-creates $A$ by a KMS-symmetric local dynamics run for a random time.

  Its classical limit is the block heat bath. Glauber dynamics is the down-up walk on the complex of partial assignments (local-to-global-unlocks §2.3–2.4, Alev–Lau, ALO20).
- **Source (memory).** The ideal quantum "up" step is the Petz recovery map, the KMS-adjoint of the inclusion $B(\mathcal H_{A^c})\hookrightarrow B(\mathcal H)$ (the Accardi–Cecchini generalised conditional expectation). CR replace it by a physical, quasi-local surrogate.
- **ANALOGY.**
  - The coloured patching of Cor. III.3 (same-colour patches disjoint and updated in parallel; local indistinguishability glues overlaps) has the shape of a *partite systematic-scan* down-up walk, and of agreement-expander gluing (Reduction 6).
  - The quantitative glue is local indistinguishability (Thm III.2), not a spectral bound.

### 7.4 Noncommutative geometry
- **KNOWN-LINK.** The OFT is a Gaussian-smoothed **Arveson spectral decomposition** of $A$ for the time evolution, which is the modular flow of $\rho_\beta$ up to the rescaling $t\mapsto-\beta t$. Lemma IX.2 is the KMS analytic-continuation identity on spectral subspaces ($\rho A_\nu\rho^{-1}=e^{-\beta\nu}A_\nu$), smoothed. Note 1 §1.2 quotes Connes: the KMS–modular relation is "one of the deepest points of contact".
- **KNOWN-LINK (memory: Cipriani–Sauvageot).** KMS-symmetric Markov semigroups are noncommutative Dirichlet forms, and these are squared norms of derivations into Hilbert bimodules.
  - CR's Lemma X.2/X.3 is such a representation, with derivation $X\mapsto([\hat A^a(\omega,t),X])_{a,\omega,t}$ into $L^2(a,\omega,t;\,gh)\otimes$ KMS-$L^2$.
  - It is the noncommutative graph Laplacian whose "edges" are the dressed jumps.
  - The expander gap is the Poincaré constant of the discrete derivation $f\mapsto f(x)-f(y)$; Note 1 §3.3's $[b,f](\gamma)=(b(r\gamma)-b(s\gamma))f(\gamma)$ is the same derivation on the quiver.
- **THEOREM-level identities (§6).** (i) CMI equals the Petz defect of $(\rho,\rho_{-A})$ on $B(\mathcal H_{AB})$. (ii) The volume prefactor is an index. (iii) The infinite-time recovery map is the $\rho$-preserving conditional expectation onto a modular-invariant subalgebra $\mathcal F_A\subseteq\mathcal N_A^\sigma$. 2609 adds the Connes cocycle $\sigma^{it}\rho^{-it}$ as the controlling object (Lemma III.4).
- **SPECULATION (memory: Hastings; Ben-Aroya–Schwartz–Ta-Shma).** "Uniformly locally gapped" (Def. B.1) is a family of *relative quantum expanders*, one per region and hereditary under restriction. A trickling-down theorem for such families (gap for single-site samplers under all restricted states ⇒ gap $c|A|^{-c'}$ for blocks) would be the Oppenheim analogue and would make App. B unconditional. No such theorem is known to me, and CR say no a-priori bound is known.

### 7.5 The barycentre-arrow programme (stated for abstract objects; drag test C2 passed unless marked)
1. **Derivation-level bridge to Note 1 §6.** By §6.3, the Markov defect of any state with respect to a tripartition is a Jenčová–Petz sufficiency defect for the pair (state, state with $A$ erased). For a layered history algebra take $A$ = a depth window, $B$ = its two boundary layers, $C$ = the rest. *Memory across a buffer is a sufficiency defect.* This holds verbatim for Bratteli path spaces and tiling patch algebras.
2. **The programme analogue of "time-averaged detailed-balanced Lindbladian with single-site jumps on A" (Derivation in the commutative resolution, Conjecture beyond).**
   - On complete histories, take the jumps on a window $W$ to be the re-routings that change a history only inside $W$ with the same endpoints: the finite part of the tail-equivalence groupoid of Note 1 §1.3, §3.1.
   - Give them Metropolis rates $\min(1,e^{-\beta\Delta F})$. Detailed balance with respect to $\mathbb P_\beta$ *is* the KMS quasi-invariance of Note 1 §5.2: Radon–Nikodym cocycle $e^{-\beta c_F}$, and $\mathbb P_\beta(\kappa\alpha'\lambda)/\mathbb P_\beta(\kappa\alpha\lambda)=e^{-\beta(F(\alpha')-F(\alpha))}$ by additivity.
   - For additive $F$ the $t=\infty$ map is exact, local Bayes resampling of the window given its boundary faces, as in the classical Ising row of §6.1.
   - For laws with memory the rates depend on the whole past and future. The CR-type statement becomes a **Conjecture**: if that dependence decays (quasi-local Radon–Nikodym cocycle), the time-averaged window dynamics recovers the window with error $e^{O(\text{window boundary})}e^{-w/\xi}$, and memory across a buffer of width $w$ decays at that rate.
   - In the commutative resolution this is elementary, since conditioning on the buffer leaves only the direct cut interactions, which is 2609's mechanism. **CR's machinery is needed only in a noncommutative resolution** (coherent carrier, frame composition with overlap kernel $G$, Note 2 §3), where the "Hamiltonian" on histories has off-diagonal terms.
3. **U1 sharpened (Conjecture).** App. B suggests replacing "local spectral expansion ⇒ transfer-operator gap" by "**uniform** local gap of window dynamics, over all restricted (pinned) laws ⇒ global Markov of the history law (poly dependence on window sizes)". The uniformity over pinnings is what makes it hereditary. U6 (spectral independence on the frame) is the natural certificate.
4. **Index versus interface as a measurable dichotomy (Conjecture, for prong 2).** For a window $W$:
   - compare $D_{\max}$ of the law relative to the $W$-erased law (volume or index scale, $\le\log\#\{\text{window paths}\}$);
   - compare the Rényi divergence to a *cut* reference that decouples $W$ from its exterior;
   - a small interface divergence is the 2609 certificate for boundary-controlled memory.
5. **SPECULATION.** mlp-bridge measured 25–35 % memory in the face process. Two sources fit the CR/2609 picture, besides a non-additive cocycle:
   - (a) the face process is a *measurement* (restriction to the commutative face algebra) of a state whose generator has coherences, so it is non-Markov even when the full state is a KMS state of a finite-range generator;
   - (b) random sub-frames are *marginals*, and marginals of Markov laws are not Markov (§7.2 guard).

   Both are charges against the dictionary (P2.1, N6), not against the theory.

### 7.6 The user's intuition, tested
- *"Gibbs states / approximate local Markov ≈ the essence of expanders."* **Partly supported.** Markov (heredity under conditioning) and expansion (a spectral certificate) are distinct, and they are the two halves of every local-to-global theorem in the programme's §4. CR prove the first unconditionally and show that the second (a uniform local gap) would give global Markov (App. B). The *mechanism* of CR's unconditional theorem is not expansion but its slow-mixing substitute, time averaging, which is the same tool the expander/spectral-independence community uses for locally stationary laws (LMR+24).
- *"Time-averaged detailed-balanced Lindbladian on single-Pauli jumps on A hides a deep analogue."* **Supported in a precise form.** It is a quantum down-up (block heat-bath) step. Its infinite-time limit is the modular-invariant conditional expectation onto the fixed-point algebra off $A$, exact but non-local for non-commuting $H$ (§6.1). The theorem is the quantitative finite-time compromise. The programme analogue is the KMS-detailed-balanced window resampling on histories, item 7.5.2. 2609 shows the dynamics is *dispensable for CMI bounds*: the modular-time average of the Connes cocycle replaces Lindblad time. It is *not dispensable for an explicit recovery channel*.
- *"NCG ties them together."* **Supported at the level of identities** (§6: Petz defect, index, modular invariance, Dirichlet forms as derivations, Arveson spectrum). A deeper unification (a trickling-down theorem for KMS samplers; quantum expanders as the certificate) is **speculation**.

---

## 8. Guards (P3.4)
- **Tensor-product and full-partition structure.** CR and 2609 both need partial traces over sites and a full partition $ABC=\Lambda$. History path spaces have hard constraints ($s(\mu_i)=r(\mu_{i+1})$). Hammersley–Clifford needs positivity, and marginals can carry memory.
- **Bounded degree and Lieb–Robinson.** Both proofs spend bounded degree on (i) the imaginary-time radius $1/2d$ and (ii) path counting for light cones. The programme has no carrier metric (P1.4). A light cone must be defined from quiver path length, and arrows are not reversible (local-to-global-unlocks §6).
- **Constants.** $\lambda\approx1/2d^4\beta^4$ makes CR's map practically non-local at low $T$ (radius $\sim d^5\beta^5$). Do not transfer the *rate*, only the *structure*.
- **The commutative case is easy.** In the programme's commutative (history-algebra) resolution, KMS states are exactly Markov and the CR theorem is trivial. Its tools matter only for noncommutative carriers.
- **v1 only;** constants $r,r'$ are not explicit in the text.

---

## 9. Messages to the other prongs

**To prong 1.**
1. (§6.3) Record "memory across a buffer = Petz sufficiency defect of the pair (state, state with the window erased)" next to the coboundary criterion of Note 1 §6. Together with the sibling's subalgebra form of 2609 Cor. III.5, this gives a *quantitative* coboundary criterion: the defect is controlled by the modular-time leakage of the Connes cocycle out of the boundary subalgebra.
2. (§6.1) Define, for a window $W$ of the history algebra, the relative commutant $\mathcal N_W$ and its modular core $\mathcal N_W^\sigma$. Exact window-Markov holds iff the KMS-symmetric window dynamics has a *local* fixed-point algebra. Ask when $\mathcal N_W^\sigma$ is local for a given cocycle.
3. Is there a trickling-down theorem for uniform local gaps of KMS-symmetric window dynamics? It would make U1 and App. B unconditional.

**To prong 2.**
1. Replace (or complement) E2's lagged $I(\sigma_{\ell+1};\sigma_{\ell-j}\mid\sigma_\ell)$ by **buffer CMI**, $I(\text{future}_{>\ell}:\text{past}_{<\ell-w}\mid\text{buffer}_{[\ell-w,\ell]})$, as a function of $w$. Exact vanishing at $w=R_0$ is the finite-range classical Gibbs signature; exponential decay is the CR/2609 class. Calibrate with the Markov surrogate (P2.4).
2. Charge memory first to marginalisation (random sub-frames, N6) and to measurement in the face basis, before charging it to the cocycle (§7.5.5).

**To prong 3.** Add to the §4 table a reduction: *global object*: recovery of an erased region / CMI; *local certificate*: none (time averaging) for local Markov, a uniform local gap for global Markov, cut strength plus light-cone velocity for the boundary version; *bound*: Thm III.1 / Cor. B.2 / 2609 Thm II.1; *cost*: $e^{O(|A|/\lambda)}$ time, $\mathrm{poly}(\beta)$ radius. Record the expander guard: shields need room, and expansion removes it.

---

## Sources read for this note
- C.-F. Chen, C. Rouzé, *Quantum Gibbs states are locally Markovian*, arXiv:2504.02208v1 (full text: §§I–XI, App. A–B, references).
- T. H. Yang, *Improved estimate of local Markovianity for quantum Gibbs states*, arXiv:2609.38007v1 (abstract, §§I–IV in full, §§V–VI skeleton, §VII, Table I, references).
- K. Kato, T. Kuwahara, *Clustering of conditional mutual information via quantum belief-propagation channels*, arXiv:2504.02235v2 (abstract, §I, Def. 1, Assumptions 3–4, Thms 2–3, Lemma 8 and its use).
- K. Liu, S. Mohanty, P. Raghavendra, A. Rajaraman, D. X. Wu, *Locally stationary distributions: a framework for analyzing slow-mixing Markov chains*, arXiv:2405.20849v4 (abstract, §1.1, Def. 1.2, Thm 1.3, Lemmas 3.1–3.2).
- Sibling digests: [expanders.md](expanders.md) (headings and §0), [chat-2609.38007-retrieval-status.md](chat-2609.38007-retrieval-status.md) (§§7–8).
- From memory, flagged where used: Tomiyama, Takesaki (conditional expectations), Pimsner–Popa and Jones (index), Cipriani–Sauvageot (Dirichlet forms as derivations), Accardi–Cecchini and Petz (recovery as KMS adjoint), quantum expanders (Hastings; Ben-Aroya–Schwartz–Ta-Shma), Martinelli's lectures beyond the sentence quoted by CR.
