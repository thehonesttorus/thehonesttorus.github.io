# Noncommutative Dirichlet forms, detailed-balanced Lindbladians, and the metric of noncommutative geometry

### The chain Lindbladian → Dirichlet form → derivation → Dirac operator → metric, which links are theorems, what the "time-averaged detailed-balanced Lindbladian based on single-Pauli jumps on A" is in this language, and what the chain gives the barycentre-arrow programme

*Prong-3 reading note (see [research-program.md](../../research-program.md)), 2026-10-01. Inputs used from the other prongs:*
- *Note 1 ([conditional-arrow-algebra.md](../../conditional-arrow-algebra.md)): the complete-history corner, the Gibbs/KMS family $\mathbb P_\beta$, the modular flow $\alpha^F_{-\beta t}$, the face algebra $D_0$, the coboundary criterion of §6, and the optional Dirac operator of §8.*
- *[local-to-global-unlocks.md](../../local-to-global-unlocks.md): the four-ingredient table of §4 and the unlocks U1, U3, U4, U5.*
- *The sibling digests in this directory: [arxiv-2609.38007.md](arxiv-2609.38007.md) (§8 gives the Chen–Rouzé construction in full; bridges B5, B7–B10, B14), [expanders.md](expanders.md) (§9.9–9.11, §11.1, §11.5) and [hdx-spectral-independence.md](hdx-spectral-independence.md) (§4.4).*

*Nothing below changes the objects of Note 1. Numerical checks: [check_nc_dirichlet_lindblad.py](check_nc_dirichlet_lindblad.py) (pure numpy; all pass, §10).*

**Labels.** **Source** means stated in a text I retrieved in this session; theorem and equation numbers are the source's own. **Mine** means my derivation or interpretation. Bridges are labelled **THEOREM**, **KNOWN-LINK**, **ANALOGY** or **SPECULATION**. *From memory* marks a statement I did not open in this session.

**Sigla.**
- [CM17] Carlen–Maas, *Gradient flow and entropy inequalities for quantum Markov semigroups with detailed balance*, arXiv:1609.01254, JFA 273 (2017).
- [CM20] Carlen–Maas, *Non-commutative calculus, optimal transport and functional inequalities in dissipative quantum systems*, arXiv:1811.04572, J. Stat. Phys. 178 (2020).
- [CKG] Chen–Kastoryano–Gilyén, *An efficient and exact noncommutative quantum Gibbs sampler*, arXiv:2311.09207.
- [DLL] Ding–Li–Lin, *Efficient quantum Gibbs samplers with Kubo–Martin–Schwinger detailed balance condition*, arXiv:2404.05998.
- [KB] Kastoryano–Brandão, *Quantum Gibbs samplers: the commuting case*, arXiv:1409.3435.
- [CR] Chen–Rouzé, *Quantum Gibbs states are locally Markovian*, arXiv:2504.02208.
- [VW] Vernooij–Wirth, *Derivations and KMS-symmetric quantum Markov semigroups*, arXiv:2303.15949 (CMP 2023).
- [Ver] Vernooij, *On the existence of derivations as square roots of generators of state-symmetric quantum Markov semigroups*, arXiv:2203.12307.
- [Wir] Wirth, *A noncommutative transport metric and symmetric quantum Markov semigroups as gradient flows of the entropy*, arXiv:1808.05419.
- [CS17] Cipriani–Sauvageot, *Amenability and subexponential spectral growth rate of Dirichlet forms on von Neumann algebras*, arXiv:1611.01749, Adv. Math. 322 (2017).
- [CGIS] Cipriani–Guido–Isola–Sauvageot, *Spectral triples for the Sierpinski gasket*, arXiv:1112.6401, JFA 266 (2014).
- [DM] D'Andrea–Martinetti, *A view on optimal transport from noncommutative geometry*, arXiv:0906.1267.
- [Rie] Rieffel, *Metrics on state spaces*, arXiv:math/9906151.
- [RD] Datta–Rouzé, *Concentration of quantum states from quantum functional and transportation cost inequalities*, arXiv:1704.02400.
- [DMTL] De Palma–Marvian–Trevisan–Lloyd, *The quantum Wasserstein distance of order 1*, arXiv:2009.04469.
- [KT] Kastoryano–Temme, *Quantum logarithmic Sobolev inequalities and rapid mixing*, arXiv:1207.3261.
- [BST] Ben-Aroya–Schwartz–Ta-Shma, *An explicit construction of quantum expanders*, arXiv:0709.0911.
- [Har] Harrow, *Quantum expanders from any classical Cayley graph expander*, arXiv:0709.1142.
- [Pet] Peterson, *A 1-cohomology characterization of property (T) in von Neumann algebras*, Pacific J. Math. 243 (2009).
- [CMo] Connes–Moscovici, *Type III and spectral triples*, arXiv:math/0609703.
- [CNNR] Carey–Neshveyev–Nest–Rennie, *Twisted cyclic theory, equivariant KK theory and KMS states*, arXiv:0808.3029.
- [Voi] Voiculescu, *Perturbations of operators and non-commutative condensers*, arXiv:2503.23248.
- [Y] Yang, arXiv:2609.38007, read through the sibling digest.

---

## 0. Retrieval status

**Read in this session.** The papers above were read through alphaXiv page retrieval. Retrieval returns the pages relevant to my queries, so not every paper was read in full; each statement below is quoted from a page that was returned. [Pet] was read as the msp.org PDF text, via Exa. For Cipriani's survey I found only the bibliographic metadata (MaRDI, Politecnico repository); no open copy turned up.

**Read only through other sources.** None of the following was opened. Each statement is taken from the source named after it, and is marked as such where it is used.
- *Cipriani–Sauvageot, Derivations as square roots of Dirichlet forms* (JFA 201, 2003): via [CS17 §2(e)], [Wir Thm 1.9 as described], [VW §1] and [Ver §1].
- Cipriani 1997 (JFA 147) and Goldstein–Lindsay 1995 (Math. Z. 219): via [VW §5] and [CS17 §2].
- Albeverio–Høegh-Krohn 1977 and Davies–Lindsay 1992/93: via [VW §1], [CM17 refs] and [RD refs].
- Fagnola–Umanità 2007/2010: via [DLL §2.2] and [CM20 §2].
- Alicki 1976: via [CM17 Thm 3.1] and [CM20 Thm 2.4].
- Sauvageot 1989/1990: via [Pet §0].
- Hastings, arXiv:0706.0556: via [BST Thm 3], [Har] and [expanders.md](expanders.md) §9.10.
- Pask–Rennie 2006: via [CNNR].
- Rouzé–França–Alhambra [RFA24]: via [CR] and [DLL]. I have not checked its identifier.

**Not retrievable.**
- The ChatGPT share page for arXiv:2609.38007. A parallel run documents eleven failed routes in [chat-2609.38007-retrieval-status.md](chat-2609.38007-retrieval-status.md). I did not retry hosts that are blocked by the proxy. Nothing here describes that conversation.
- `C:\Users\User\Downloads\expander_survey.pdf` is on the user's machine, not in this container.

**What I re-verified from primary text.** The sibling digests are machine-generated, so I re-read the following in [CR] myself:
- Lemmas VII.1–VII.3, VIII.1 and X.1–X.4;
- Theorem II.1, Theorem III.1 and its proof;
- Appendix B.

What §4 says about [CR] rests on that reading.

---

## 1. The answer in brief

**What the phrase means.** "Time-averaged detailed-balanced Lindbladian based on single-Pauli jumps on A" is the caption of [CR] Fig. 2. It names
$$
\mathcal R_{A,t}=\frac1t\int_0^t e^{s\mathcal L_A}\,ds,\qquad \mathcal L_A=\sum_{a\in P^1_A}\mathcal L_a,\qquad P^1_A=\{X_i,Y_i,Z_i\}_{i\in A}.
$$
Each $\mathcal L_a$ is the [CKG] generator with the single coupling $A^a$. The three words correspond to three separate mathematical facts:
- **"detailed-balanced"** means KMS-symmetric, not GNS-symmetric. It is exact for every Hamiltonian, thanks to a coherent correction term.
- **"single-Pauli jumps on A"** means the couplings form a generating set of the local algebra $M_A$, so the commutant of the couplings is $M_{A^c}$.
- **"time-averaged"** means a Cesàro mean of the semigroup. This is distinct from two other time averages in the same literature (§4.2).

**What it does.** It is the explicit recovery map in [CR] Thm III.1: $\|\mathcal R_{A,t}[\rho_{-A}]-\rho_\beta\|_1\le re^{\mu|A|}t^{-\lambda}$. Together with quasi-locality, this gives the decay of conditional mutual information. In Dirichlet-form language (Mine, assembled from [CR]), the mechanism has five steps:
1. KMS symmetry makes the generator a symmetric Dirichlet form.
2. Cesàro averaging makes that form $\le c/t$ with no spectral gap. This is the mean ergodic theorem.
3. The form is the KMS-weighted energy of a derivation whose components are commutators with Gibbs-filtered single Paulis dressed in modular time.
4. A small derivation means approximate commutation with every single Pauli on $A$, hence approximate triviality on $A$. This step loses a Hölder exponent and a factor $2^{O(|A|)}$ for word length.
5. Observables that are trivial on $A$ cannot tell $\rho$ from $\rho_{-A}$.

**The chain is a chain of theorems in finite dimensions, with one branch missing (§3).**

| step | object | what makes it hold |
|---|---|---|
| 1 | KMS-detailed-balanced Lindbladian | theorem |
| 2 | completely Dirichlet form on the standard form | theorem |
| 3 | derivation into a bimodule, twisted by $\sigma_{\mp i/4}$ | theorem ([VW]) |
| 4 | Dirac operator | construction; trivially a spectral triple in finite dimension |
| 5 | Connes metric, a Wasserstein-1 | metric exactly when the semigroup is primitive (elementary) |

The **missing branch** is the Carlen–Maas Wasserstein-2 gradient-flow metric. It is constructed only for GNS-symmetric generators [CM20 §2.4], and the [CKG]/[CR] generators are KMS-symmetric but not GNS-symmetric (check K1).

**Verdict on the three intuitions.**
1. *Gibbs/Markov and expanders share an essence.* This is right at exactly one place. Uniform coercivity of the "vertical" derivations on regions $A$ (a uniform local gap, equivalently a uniform Kazhdan-type constant for single-site generators) turns the proved local Markov property into the open global one ([CR] Cor B.2, Remark III.1.2). The Markov property itself is not an expansion property: [CR] and [Y] prove it with no gap.
2. *The time-averaged single-Pauli Lindbladian hides a deep analogue.* Yes, and it can be made a theorem inside the programme. The analogue is a **re-routing Lindbladian** on the complete-history corner, GNS-detailed-balanced for $\mathbb P_\beta$ with jumps the re-routing matrix units (§6.1, check K7). Its kernel is the face algebra $D_0$. Its Cesàro limit is the $\mathbb P_\beta$-preserving conditional expectation onto $D_0$. It recovers $\mathbb P_{\beta'}$ for all $\beta'$ iff the barycentre is sufficient, and the recovery error is bounded by $2\tanh(|\beta-\beta'|\,\mathrm{osc}F/4)$. This is the Chen–Rouzé construction with every approximation removed: Note 1's modular Hamiltonian is diagonal in the history basis, which is the "commuting case".
3. *NCG ties these together.* Partly right. The shared object is a **first-order differential calculus**, a derivation into a Hilbert bimodule. It carries:
   - the Dirichlet form, as the $L^2$ energy;
   - the Connes metric, through the $L^\infty$ Lipschitz seminorm;
   - the spectral gap, as coercivity off the kernel;
   - property (T), since (T) holds iff all closable derivations are inner [Pet];
   - amenability, through subexponential spectral growth [CS17];
   - the Markov conditional expectations, as kernels.

   The guards are in §6.7. The quarter-twist by the modular group cannot be removed for non-tracial states [Ver]. The Wasserstein-2 branch is missing in the KMS case. Finite-dimensional NCG is topologically trivial, so all content is in constants. The single-Pauli generator at $\beta=0$ is not a quantum expander.

---

## 2. Definitions

### 2.1 Modular data in finite dimensions (Source unless marked)

- **The state and its modular objects.** $\sigma$ is a full-rank density matrix. The modular operator is $\Delta_\sigma(X)=\sigma X\sigma^{-1}$ [CM20 Def 1.1]. The modular group, with its analytic continuation, is $\sigma_z(X)=\sigma^{iz}X\sigma^{-iz}$ [VW §2]. The Gibbs state is $\rho_\beta=e^{-\beta H}/Z$.
- **Bohr components** [CKG §II.A]. They are $A_\nu=\sum_{E_2-E_1=\nu}P_{E_2}AP_{E_1}$, satisfying $e^{iHt}Ae^{-iHt}=\sum_\nu A_\nu e^{i\nu t}$, $[H,A_\nu]=\nu A_\nu$ and $(A_\nu)^\dagger=(A^\dagger)_{-\nu}$. *Mine:* hence $\Delta_{\rho_\beta}A_\nu=e^{-\beta\nu}A_\nu$, so Bohr components are exactly the eigenvectors of the modular operator.
- **Inner products** [CM20 (2.3)–(2.4)]:
  - GNS: $\langle X,Y\rangle_{\rm GNS}=\mathrm{Tr}(\sigma X^*Y)$;
  - KMS: $\langle X,Y\rangle_{\rm KMS}=\mathrm{Tr}(X^*\sigma^{1/2}Y\sigma^{1/2})$;
  - BKM: $\langle X,Y\rangle_{\rm BKM}=\int_0^1\mathrm{Tr}(X^*\sigma^sY\sigma^{1-s})ds$;
  - the general family $\langle A,B\rangle_f=\mathrm{Tr}[A^*f(\Delta_\sigma)B\sigma]$ [CM17 (2.5)].

  Index conventions differ between sources. In [CM17] and [DLL] the GNS case is $s=1$; in [CM20] it is $s=0$.
- **The KMS (symmetric) embedding.** *Mine, one line:* $X\mapsto\sigma^{1/4}X\sigma^{1/4}$ is an isometry from $(M_n,\langle\cdot,\cdot\rangle_{\rm KMS})$ onto Hilbert–Schmidt space. Four literatures use this one embedding (KNOWN-LINK):
  - it is $i_\omega(x)=\Delta_\omega^{1/4}x\xi_\omega$ with $\xi_\omega=\sigma^{1/2}$, the "symmetric embedding" of noncommutative potential theory [CS17 §2];
  - [KB] (79) uses it to map a Gibbs sampler to a frustration-free Hamiltonian with ground state $|\sqrt\rho\rangle$: $\hat{\mathcal L}(f)=\rho^{1/4}\mathcal L(\rho^{-1/4}f\rho^{-1/4})\rho^{1/4}$ is Hermitian (their Table I);
  - [CKG] Def II.1 calls $D(\rho,\mathcal L)=\rho^{-1/4}\mathcal L[\rho^{1/4}\cdot\rho^{1/4}]\rho^{-1/4}$ the "discriminant" and purifies it into "frustration-free parent Hamiltonians";
  - [KT] Def 6 states detailed balance as $\Gamma_\sigma\circ\mathcal L=\mathcal L^*\circ\Gamma_\sigma$. I read $\Gamma_\sigma(X)=\sigma^{1/2}X\sigma^{1/2}$ off their identity $\mathrm{Var}_\sigma(\Gamma_\sigma^{-1}\rho)=\chi^2(\rho,\sigma)$ in the proof of Thm 22; the definition itself was not on the retrieved pages.

### 2.2 Quantum Markov semigroups and detailed balance (Source)

A QMS is a semigroup $P_t=e^{t\mathcal L}$ of unital completely positive maps. Its generator has Lindblad form [DLL (2.3)]: $\mathcal L(X)=i[G,X]+\sum_j(L_j^\dagger XL_j-\tfrac12\{L_j^\dagger L_j,X\})$.

**Definitions.**
- **GNS detailed balance** (σ-DBC): each $P_t$ is self-adjoint for $\langle\cdot,\cdot\rangle_{\rm GNS}$ [CM17 Def 2.10; CM20 Def 2.3].
- **KMS detailed balance:** each $P_t$ is self-adjoint for $\langle\cdot,\cdot\rangle_{\rm KMS}$. Equivalently $\mathcal L^\dagger[\cdot]=\rho^{-1/2}\mathcal L[\rho^{1/2}\cdot\rho^{1/2}]\rho^{-1/2}$ [CKG Def II.1].
- **BKM detailed balance** is defined in the same way [CM20 Def 2.2].

**Facts.**
- [CM17 Thm 2.9; CM20 Lemma 2.1, credited to Alicki]. If a real map is self-adjoint for $\langle\cdot,\cdot\rangle_s$ for some $s\neq\tfrac12$, then it commutes with $\Delta_\sigma$ and is self-adjoint for every $s$, KMS included. Hence GNS implies KMS. The converse fails: [CM17 App. B] builds KMS-symmetric generators on $M_2$ that do not commute with $\Delta_\sigma$.
- [CM20 Thm 2.9]. If an ergodic QMS's forward equation is the gradient flow of $\mathrm{Ent}_\sigma$ for some $C^1$ Riemannian metric, then each $P_t$ is BKM-symmetric. The BKM-symmetric class strictly contains the GNS class (§2.4). "Only when each $P_t$ is self-adjoint with respect to the GNS inner product do we have a construction of such a Riemannian metric."
- [CKG App. E]. Take a transition part $\sum\alpha_{\nu_1\nu_2}A_{\nu_1}(\cdot)A_{\nu_2}^\dagger$. Detailed balance for the $s$-inner product forces $\alpha_{\nu_1,\nu_2}=\alpha_{-\nu_2,-\nu_1}e^{-\beta(1-s)\nu_1-\beta s\nu_2}$. For $s\neq\tfrac12$ and $\beta\neq0$ this forces $\alpha_{\nu_1\nu_2}=0$ whenever $\nu_1\ne\nu_2$ (E5). In their words: "the only existing Lindbladian that satisfies (E5) is the Davies' generator, which requires resolving the level spacing using a (exponentially) long Hamiltonian simulation time."

**Structure theorem, GNS case** (Alicki; [CM17 Thm 3.1], [CM20 Thm 2.4]). $P_t$ satisfies GNS detailed balance iff
$$
\mathcal L=\sum_je^{-\omega_j/2}\mathcal L_j,\qquad \mathcal L_j(A)=V_j^*[A,V_j]+[V_j^*,A]V_j,
$$
with $\{V_j\}=\{V_j^*\}$ and $\Delta_\sigma V_j=e^{-\omega_j}V_j$; equivalently $[V_j,H]=-\omega_jV_j$ with $H=-\log\sigma$. [CM17] adds the normalisations $\tau[V_j^*V_k]=\delta_{jk}$ and $\tau[V_j]=0$. **Converse:** any family $\{V_j\}=\{V_j^*\}$ of $\Delta_\sigma$-eigenvectors gives, through (3.3), a QMS satisfying GNS detailed balance.

**Structure theorem, KMS case** ([DLL] Thm 10, revisiting Fagnola–Umanità [FU07 Thms 7.2–7.3] along Amorim–Carlen [AC21]). $\mathcal L$ satisfies KMS detailed balance iff
$$
\mathcal L(X)=i[G,X]+\sum_j\big(L_j^\dagger XL_j-\tfrac12\{L_j^\dagger L_j,X\}\big),\qquad \Delta_\sigma^{-1/2}L_j=L_j^\dagger,\qquad G=-i\tanh\circ\log(\Delta_\sigma^{1/4})\Big(\tfrac12\sum_jL_j^\dagger L_j\Big).
$$
The condition on the jumps is equivalent to $\Delta^{-1/4}L_j$ being self-adjoint, that is $L_j=\Delta^{1/4}\tilde A_j$ with $\tilde A_j=\tilde A_j^\dagger$ ([DLL] (3.1)). [DLL] Lemma 9 gives a further equivalent: the completely positive part and the $K^\dagger X+XK$ part are separately KMS-self-adjoint.

**Kernel** ([DLL] Lemma 3, citing Wolf Thm 7.2). If a full-rank invariant state exists, then $\ker\mathcal L=\{G,L_j,L_j^\dagger\}'$. Primitivity means this commutant is $\mathbb C1$.

### 2.3 The Gibbs samplers of this literature (Source)

1. **Davies generator** [DLL (1.2); KB §IV.A]. It comes from the weak-coupling limit:
   $$\mathcal L^\dagger_{\rm Davies}(\rho)=-i[H+H_{LS},\rho]+\sum_{a,\omega}\gamma_a(\omega)\big(\hat A^a(\omega)\rho\hat A^a(\omega)^\dagger-\tfrac12\{\hat A^a(\omega)^\dagger\hat A^a(\omega),\rho\}\big),$$
   with $\gamma_a(-\omega)=e^{\beta\omega}\gamma_a(\omega)$ and $\sigma_\beta\hat A^a(\omega)\sigma_\beta^{-1}=e^{-\beta\omega}\hat A^a(\omega)$. The source says:
   - "the Davies semigroup is exactly the class of QMS with GNS-DBC, up to the coherent term";
   - "derived from a secular approximation";
   - implementing it "requires accurately resolving all the Bohr frequencies … impractically long Hamiltonian simulation time" ([DLL] Remark 7).

   For commuting $H$ it is local ([KB] Lemma 11(2)).
2. **Heat-bath generator** [KB (34)]. $\mathcal L^H_A(f)=\sum_{k\in A}(E^\rho_k(f)-f)$ with $E^\rho_k(f)=\mathrm{tr}_k[\eta^\rho_kf\eta^{\rho\dagger}_k]$ and $\eta^\rho_k=(\mathrm{tr}_ke^{-\beta H})^{-1/2}e^{-\beta H/2}$. It is local for commuting potentials ([KB] Lemma 12).
3. **The [CKG] sampler.** It uses the operator Fourier transform with Gaussian filter $f(t)\propto e^{-\sigma_E^2t^2}$, the Metropolis weight $\gamma(\omega)=\exp(-\beta\max(\omega+\tfrac1{2\beta},0))$ at $\sigma_E=1/\beta$, and a coherent term $B$; [arxiv-2609.38007.md](arxiv-2609.38007.md) §8.1 reproduces all of it.
   - [CKG] Thm I.1: exact KMS detailed balance, hence $\mathcal L_\beta[\rho_\beta]=0$.
   - [CKG] Thm I.2: $e^{\mathcal L_\beta t}$ is implementable at a cost of $\tilde O(t\beta)$ total Hamiltonian simulation time.
   - [DLL] generalises it to finitely many jumps and any filter $q$ with $q(-\nu)=\overline{q(\nu)}$: $L_a=\sum_\nu q_a(\nu)e^{-\beta\nu/4}A^a_\nu=\int f_a(t)A^a(t)dt$, with $G$ given by (3.10). It "encompasses the construction of Chen, Kastoryano, and Gilyén as a special instance" (§4).
   - [DLL] Prop 20: quasi-locality at every temperature, $\|L_a-L_a^{(r)}\|\le\|A_a\|\big(O(J\beta/\sqrt r)\big)^r$.

**[KB] main theorem** (Thms 23 and 26; informal statement Thm 1). For commuting local Hamiltonians, the Davies and heat-bath generators have a gap independent of system size iff the Gibbs state satisfies strong clustering (Def 15):
$$\mathrm{Cov}_{A\cup B}(E_A f,E_Bf)\le c\|f\|_{2,\rho}^2e^{-d(B\setminus A,A\setminus B)/\xi}.$$
- This always holds in 1D (§VII) and above a size-independent critical temperature (Thms 30–31).
- The source is explicit about the limitation (§IX): "in general $L^D,L^H,E^\rho$, and $E^L$ all become non-local when $H$ is non-commuting, and very little of the framework can be recovered."

### 2.4 Noncommutative Dirichlet forms (Source)

**Definition** ([CS17] Def 2.1, for a σ-finite von Neumann algebra $N$ with faithful normal state $\omega$ and standard form $(N,L^2(N,\omega),L^2_+,J_\omega)$). A densely defined, nonnegative, lower semicontinuous quadratic form $\mathcal E$ on $L^2$ is:
- **real** if $\mathcal E[J_\omega\xi]=\mathcal E[\xi]$;
- **Dirichlet** if it is real and **Markovian**: $\mathcal E[\xi\wedge\xi_\omega]\le\mathcal E[\xi]$ for real $\xi$, where $\xi\wedge\xi_\omega$ is the Hilbert projection onto $C_\omega=\{\eta=J\eta:\ \xi_\omega-\eta\in L^2_+\}$;
- **completely Dirichlet** if every matrix amplification is Dirichlet.

**Correspondence theorem.**
- Dirichlet forms correspond to Markovian self-adjoint semigroups on $L^2$, and through $i_\omega$ to "modular $\omega$-symmetric" completely positive contractive semigroups on $N$: $\omega(S_t(x)\sigma^\omega_{-i/2}(y))=\omega(\sigma^\omega_{-i/2}(x)S_t(y))$ ([CS17] (2.4), citing Cipriani 1997).
- [VW §5]: "There is a one-to-one correspondence between quantum Dirichlet forms on $L^2(M)$ and KMS-symmetric quantum Markov semigroups on $M$ (see [Cip97, Theorem 4.11], [GL95, Theorem 5.7])." Here "quantum Dirichlet form" means a conservative completely Dirichlet form.
- In finite dimension [CM20 §2.2]: "The class of bilinear forms $\mathcal E$ defined in terms of a self-adjoint QMS … through (2.9) is, by definition, the class of conservative completely Dirichlet forms." Their §3 develops Beurling–Deny theory in finite dimension. The source flags its Theorem 3.8 as new; it "singles out the KMS inner product". I did not retrieve its statement.

**History** (from the sources). The tracial theory was initiated by Albeverio–Høegh-Krohn (1977) and developed by Davies–Lindsay (1992, 1993) [VW §1]. Sauvageot (1989, 1990) "makes a connection between quantum Dirichlet forms and differential calculus" [Pet §0].

### 2.5 Derivations, bimodules, carré du champ (Source)

A **Hilbert bimodule** carries commuting left and right $*$-actions; a **derivation** satisfies $\delta(AB)=A\delta(B)+\delta(A)B$ [Ver Defs 2.1–2.3].

**Theorem (tracial; Cipriani–Sauvageot 2003, as stated in [CS17 §2(e)] and [Ver §1]).**
- The $L^2$ generator of a $\tau$-symmetric QMS on a C\*-algebra satisfies $L^{(2)}=\delta^*\delta$ for a densely defined closable derivation $\delta$ into a Hilbert bimodule. The derivation is "essentially unique".
- Conversely, the closure of $\|\partial a\|^2$ is a Dirichlet form whenever $\partial$ is a closable derivation.

**Theorem (KMS, finite dimension; [VW] Thm 2.4).** Let $\mathcal L$ be the generator of a KMS-symmetric QMS on $M_n$, with [VW]'s sign convention $\mathcal L=\lim(\mathrm{id}-\Phi_t)/t\ge0$. There exist:
- a Hilbert space $\mathcal H$ with a unital $*$-homomorphism $\pi_l$ and a commuting unital $*$-antihomomorphism $\pi_r$;
- an antilinear isometric involution $J$;
- a map $\delta$ with $\delta(A^*)=J\delta(A)$ and $\mathcal H=\mathrm{span}\,\pi_l(A)\delta(B)$,

satisfying
$$
\delta(AB)=\pi_l(\sigma_{-i/4}(A))\delta(B)+\pi_r(\sigma_{i/4}(B))\delta(A),\qquad \langle A,\mathcal L(B)\rangle_\rho=\langle\delta(A),\delta(B)\rangle_{\mathcal H}.
$$
**[VW] Thm 2.5:** there are $V_1,\dots,V_N$ with $\{V_j\}=\{V_j^*\}$ and $\langle A,\mathcal L(B)\rangle_\rho=\sum_j\langle[V_j,A],[V_j,B]\rangle_\rho$.

Beyond finite dimension:
- [VW] Thm 4.2 treats uniformly continuous semigroups on general von Neumann algebras; there $\delta$ is bounded and inner, $\delta(a)=\pi_l(a)\xi_0-\xi_0J\pi_r(a)^*J$.
- [VW] Thm 5.4 treats arbitrary KMS-symmetric semigroups: a closed $\delta$ on $\mathrm{dom}(\mathcal L_2^{1/2})$ with $\mathcal L_2=\delta^*\delta$. Uniqueness is open (Remark 5.6).
- The key tool is the "$\mathcal V$-transform", which solves $\tfrac12(\Delta^{1/4}S\Delta^{-1/4}+\Delta^{-1/4}S\Delta^{1/4})=T$ and preserves positivity.

**No-go theorem** ([Ver], Examples 5.3–5.4). For GNS-symmetric generators on $M_3$, a derivation into a Hilbert bimodule (with $*$-actions) satisfying $L^{(2)}=\delta^*\delta$ need not exist, whether the $L^2$ structure is GNS or KMS. "Our results can be seen as proof that this generalisation to twisted derivations cannot be avoided."

**GNS case, explicit** ([CM20] Prop 2.5 and (2.13)). With $\partial_jA=[V_j,A]$ from Alicki's form,
$$\mathcal LA=-\sum_j\partial^\dagger_{j,\sigma}\partial_jA,\qquad\text{i.e.}\qquad -\langle A,\mathcal LB\rangle_{\rm KMS}=\sum_j\langle\partial_jA,\partial_jB\rangle_{\rm KMS}.$$
In their words, "it is the KMS inner product that is more natural: the Dirichlet form … can be expressed in terms of a 'squared gradient'." Prop 4.11 gives $\ker\mathcal L=\ker\nabla$.

**Carré du champ.** Classically [CGIS (4.2)]: $\Gamma(f,g)(h)=\tfrac12(\mathcal E(f,hg)-\mathcal E(fg,h)+\mathcal E(g,fh))$, "(up to the constant 1/2) the Hochschild co-boundary of the 1-cocycle $\phi(f_0,f_1):=\mathcal E(f_0,f_1)$". Noncommutatively, $\Gamma(a)(x)=\langle x\partial a,\partial a\rangle_{\mathcal H}$ [Wir §2].

### 2.6 Spectral triples, Lip-norms, Connes distance, twisted and modular triples (Source)

- **Spectral triple** [DM §1; CGIS §1]: $(A,\mathcal H,D)$ with $[D,a]$ bounded and $D$ of compact resolvent.
- **Connes distance** [DM (1.1)]: $d_D(\varphi,\psi)=\sup\{|\varphi(a)-\psi(a)|:\|[D,a]\|\le1\}$.
- **Lip-norm** [Rie Def 5.1]: a seminorm $L$ with
  1. $L(a)=0$ iff $a\in\mathbb Re$;
  2. $L$ lower semicontinuous;
  3. $\{L\le1\}$ totally bounded in $A/\mathbb Re$.

  In finite dimension, $L(a)=\|[D,a]\|$ is a Lip-norm iff $[D,a]=0$ implies $a\in\mathbb CI$ [Rie §7].
- **Radius criterion** [Rie Prop 2.2]: $\rho_L\le2r$ on $S(A)$ iff $\|\tilde a\|^\sim\le rL(a)$ for all $a$, where $\|\tilde a\|^\sim=(\max a-\min a)/2$. The radius of the state space is the best constant in an "oscillation ≤ Lipschitz" inequality.
- **Commutative case** [DM Prop 2.1]: on a complete Riemannian spin manifold, $d_D=W_1$ on all states, by Kantorovich duality.
- **Noncommutative case** [DM]: pure states do not form a path metric space. On the other hand, $d_D(\varphi_s,\varphi_t)=|s-t|\,d_D(\varphi_0,\varphi_1)$ along segments of $S(A)$ (1.9).
- **Twisted spectral triples** [CMo Def 3.1]:
  - $Da-\sigma(a)D$ is bounded, with unitarity $\sigma(a^*)=(\sigma^{-1}(a))^*$;
  - under (1PG), $\sigma$ is the value at $t=-i$ of a one-parameter group, namely the modular group in the foliation example (Remark 2.4);
  - $d_\sigma(a)=Da-\sigma(a)D$ is a derivation into the bimodule $a\cdot\omega\cdot b=\sigma(a)\omega b$ (Prop 3.4);
  - the Chern character lands in ordinary cyclic cohomology.
- **Modular spectral triples** [CNNR]. The data are a circle action $\sigma$ on $A$, the fixed-point algebra $F$, the expectation $\Phi(a)=\frac1{2\pi}\int\sigma_t(a)dt$, and a KMS$_\beta$ weight $\phi$.
  - $D$, the generator of $\sigma$ on the C\*-module completion over $F$, defines $[D]\in KK^{\mathbb T}_1(A,F)$.
  - $\phi_D=\mathrm{Tr}_\phi(e^{-\beta D/2}\cdot e^{-\beta D/2})$ is a weight whose modular group is implemented by $D$ and which is a trace on $N^\sigma$.
  - The spectral flow of "modular partial isometries" is a residue of a twisted cyclic cocycle (Thms 1.1–1.4).
  - Examples: on Cuntz $O_n$, $sf(S_\alpha S_\alpha^*D,S_\alpha S_\beta^*DS_\beta S_\alpha^*)=(|\beta|-|\alpha|)n^{-|\alpha|}$; on Araki–Woods factors, $sf=-n(1+e^\beta)^{-n}$.
  - "For modular unitaries $u_v$, $sf_{\phi_D}(D,u_vDu_v^*)$ is just Araki's relative entropy of the two KMS weights $\phi_D$ and $\phi_D\circ\mathrm{Ad}\,u_v$."
  - The graph-algebra case is Pask–Rennie [28 in CNNR], which I did not open.
- **Spectral triples from Dirichlet forms** [CGIS].
  - Circle (Thm 3.8): $D_\alpha=\begin{pmatrix}0&\partial_\alpha\\\partial^*_\alpha&0\end{pmatrix}$ with $\partial^*_\alpha\partial_\alpha=\Delta^\alpha$ gives a spectral triple of dimension $1/\alpha$, and the energy $\mathcal E_\alpha$ is recovered from $D$ by a residue formula.
  - Sierpiński gasket, zeta function (Thm 4.3): $Z_D(s)=4\zeta(\alpha s)/(1-3\cdot2^{-s})$.
  - Sierpiński gasket, metric (Cor 5.3): $f\mapsto\|[D,f]\|$ is a Lip-norm in Rieffel's sense, and for $\beta>\alpha$ the Connes metric is bi-Lipschitz to $\rho_{\rm geo}^\beta$.
  - Sierpiński gasket, energy (Thm 5.5): at the energy dimension $\delta_D=\max\{\alpha^{-1},2-\frac{\log5/3}{\beta\log2}\}$ the residue $\mathrm{Res}\,\mathrm{tr}(|D|^{-s/2}|[D,f]|^2|D|^{-s/2})$ equals a constant times the standard Dirichlet form. The energy dimension $d_E=\log(12/5)/\log2\approx1.26$ is smaller than $d_H=\log3/\log2\approx1.58$.
  - Thm 5.4: the Fredholm module pairs non-trivially with $K_1$.

### 2.7 Transport metrics (Source)

**Carlen–Maas $W_2$** [CM17 Def 7.1, Thms 7.5–7.6].
- The "multiplication by $\rho$" operator is $[\rho]_\omega(A)=\int_0^1e^{\omega(1/2-s)}\rho^sA\rho^{1-s}ds$ ([RD] Lemma 2).
- A path $\rho_t$ with velocity field $V$ satisfies the continuity equation $\dot\rho+\mathrm{div}([\rho]_{\vec\omega}V)=0$, and $W_{2,\mathcal L}$ is the Benamou–Brenier infimum of $\int_0^1\|V\|^2_{\mathcal L,\rho_t}dt$ over such paths.
- Under GNS detailed balance, the forward equation is the gradient flow of $D(\cdot\|\sigma)$.

**Wirth's infinite-dimensional version** [Wir]. It builds $W_2$ from the Cipriani–Sauvageot derivation for tracially symmetric quantum Dirichlet forms. Under the Bakry–Émery-type estimate GE$(K,\infty)$ it proves an EVI$_K$ gradient flow (Thm 6.26), existence of geodesics (Thm 7.7) and $K$-convexity of the entropy (Thm 7.12). Example 4.19 leaves open "how one can generalize the construction of the metric $W$ to the case of infinite-dimensional quantum Markov semigroups satisfying the detailed balance condition for a non-tracial state."

**Datta–Rouzé $W_1$** [RD].
- Lipschitz constant (2.32): $\|f\|_{\rm Lip}=\big(\tfrac1d\sum_jc_j(e^{-\omega_j/2}+e^{\omega_j/2})\|\partial_jf\|_\infty^2\big)^{1/2}$.
- Distance: $W_{1,\mathcal L}(\rho,\sigma)=\sup\{|\mathrm{Tr}f(\rho-\sigma)|:\|f\|_{\rm Lip}\le1\}$, which is a Connes-type distance.
- Lemma 6: $W_{1,\mathcal L}\le\sqrt d\,W_{2,\mathcal L}$.
- Thm 3: TC$_2(c_2)$ implies TC$_1(dc_2)$.
- Thm 4 (quantum Otto–Villani): MLSI$(\alpha_1)$ implies TC$_2(1/\alpha_1)$.
- Thm 6: TC$_2$ implies the Poincaré inequality.
- Thm 8: TC$_1$ implies Gaussian concentration.

**De Palma–Marvian–Trevisan–Lloyd $W_1$** [DMTL].
- Two states are **neighbouring** if they coincide after one qudit is discarded. The norm $\|X\|_{W_1}$ is (8), and $\tfrac12\|X\|_1\le\|X\|_{W_1}\le\tfrac n2\|X\|_1$ (Prop 2).
- Prop 6: on diagonal states it equals the Hamming $W_1$, Ornstein's $\bar d$.
- Dual Lipschitz constant (Prop 8): $\|H\|_L=2\max_i\min_{H^{(i)}}\|H-I_i\otimes H^{(i)}\|_\infty$.
- Prop 15: $\frac{d^2}{d^2-1}\max_i\|H-E_iH\|_\infty\le\|H\|_L\le2\max_i\|H-E_iH\|_\infty$, where $E_i$ replaces the $i$th qudit by the maximally mixed state.
- Thm 2 (quantum Marton): $\|\rho-\sigma\|_{W_1}\le\sqrt{\tfrac n2S(\rho\|\sigma)}$ for product $\sigma$.
- Thm 3: Gaussian concentration.

### 2.8 Gap, log-Sobolev, mixing [KT] (Source)

- **Gap** (Def 10): $\lambda=\min\mathcal E_2(g)/\mathrm{Var}_\sigma(g)$.
- **Log-Sobolev constants** (Def 11): LS$_p$. Thm 15: LS$_2$ is equivalent to hypercontractivity. Thm 16: $\alpha_1\le\lambda$ for primitive reversible generators.
- **Mixing** (Thm 22): $\|\rho_t-\sigma\|_{tr}\le\sqrt{1/\sigma_{\min}}\,e^{-\lambda t}$, and $\le\sqrt{2\log(1/\sigma_{\min})}\,e^{-\alpha_1t}$.
- **Depolarising semigroup:** $\lambda=\gamma$ and $\alpha_1\ge\gamma/2$. Mixing therefore takes $O(\log\log d)$, against $O(\log d)$ from the χ² bound.
- **Quantum expanders.** For primitive, reversible, unital generators of Kraus rank $D$ ((144)):
  $$\frac{2(1-2/d)\lambda}{\log(d-1)}\le\alpha_2\le\log D\,\frac{4+\log\log d}{2\log(3d/4)}.$$
  The source: "for expanders, even though the spectral gap is asymptotically independent of the dimension $d$, the LS$_2$ constant will always decrease logarithmically with the dimension."

### 2.9 Quantum expanders (Source)

- **Definition** [BST Defs 1.1–1.2]. $G=\frac1D\sum_dU_d\cdot U_d^\dagger$ is a $(N,D,\lambda)$ quantum expander if it fixes $\tilde I$ and $\|G(A)\|_2\le\lambda\|A\|_2$ on traceless $A$.
- **Zig-zag** [BST] Thm 1: from a $(N_1,D_1,\lambda_1)$ and a $(D_1,D_2,\lambda_2)$ expander one gets a $(N_1D_1,D_2^2,\lambda_1+\lambda_2+\lambda_2^2)$ expander. Thm 2: explicit $(D^{8t},D^2,\lambda+O(\lambda^2))$ families.
- **Existence** [BST] Thm 3, citing Hastings: a $(D^8,D,4\sqrt{D-1}/D)$ expander exists; footnote: $(1+O(D^{-16/15}\log D))\frac{2\sqrt{D-1}}D$.
- **Harrow's construction** [Har] (3). $\mathcal E(\rho)=\frac1{|\Gamma|}\sum_{g\in\Gamma}r_\lambda(g)\rho r_\lambda(g)^\dagger$ on an irrep satisfies $\lambda_2(\mathcal E)\le\lambda_2(W_\Gamma)$, where $W_\Gamma$ is the Cayley walk. The proof decomposes $V_\lambda\otimes V_\lambda^*$ into irreps.
- Hastings' random-unitary limit $2\sqrt{D-1}/D$ and the quantum Alon–Boppana bound are in [expanders.md](expanders.md) §9.10.

---

## 3. The chain, link by link

### 3.0 Status of each link

| # | link | status | where |
|---|---|---|---|
| L1 | For every $H$ and every self-adjoint set of couplings there is a KMS-detailed-balanced Lindbladian fixing $\rho_\beta$ exactly, quasi-local for local $H$ | THEOREM | [CKG] Thm I.1; [DLL] Thm 10, Prop 20; check K1 |
| L2 | KMS-symmetric QMS correspond to conservative completely Dirichlet forms on the standard form, via $X\mapsto\rho^{1/4}X\rho^{1/4}$ | THEOREM | [Cip97 Thm 4.11], [GL95 Thm 5.7], via [VW §5]; [CM20 §2.2] in finite dimension |
| L3 | The Dirichlet form is $\|\delta\cdot\|^2$ for a derivation into a bimodule, twisted by $\sigma_{\mp i/4}$; in finite dimension, $\sum_j\|[V_j,\cdot]\|^2_{\rm KMS}$ with $\{V_j\}=\{V_j^*\}$ | THEOREM | [VW] Thms 2.4, 2.5, 4.2, 5.4; tracial: Cipriani–Sauvageot 2003; GNS: [CM20] Prop 2.5 |
| L3′ | For the [CKG]/[CR] sampler the $V$'s are explicit: Gibbs-filtered single Paulis dressed in modular time, with the KMS-strip weight $g(t)=1/(\beta\cosh(2\pi t/\beta))$ | THEOREM (source) + check | [CR] Lemmas X.1–X.3 (citing [RFA24] Lemma C.2); my [DLL]-form version §3.2, check K2 |
| L3″ | An untwisted derivation into a Hilbert bimodule can fail to exist for non-tracial states | THEOREM | [Ver] Ex. 5.3 |
| L4 | Derivation to Dirac operator: $D=\begin{pmatrix}0&\partial^*\\\partial&0\end{pmatrix}$, or $D=\bigoplus_j\begin{pmatrix}0&V_j\\V_j^*&0\end{pmatrix}$ for inner derivations | construction; spectral-triple axioms trivial in finite dimension (Mine); infinite-dimensional theorems only in examples | [CGIS] Thm 3.8, Cor 5.3, Thm 5.5 |
| L5 | The Connes distance of that $D$ is a genuine metric on $S(A)$ iff the QMS is primitive | THEOREM (elementary, Mine) | [Rie §7] + [DLL] Lemma 3 / [CM20] Prop 4.11 |
| L6 | A spectral gap $\lambda$ gives $d_D(\rho,\psi)\le\sqrt{N/(2\lambda\rho_{\min})}$ for all states $\psi$ | THEOREM (elementary, Mine) | §3.3; check K8 |
| L7 | Connes distance = $W_1$ | THEOREM classically ([DM] Prop 2.1). Quantum: [RD]'s $W_1$ is Connes-type; at $\beta=0$ the single-Pauli Connes distance and the [DMTL] $W_1$ agree within factors $2/3$ and $3/2$ (Mine, check K6) | §3.4 |
| L8 | Carlen–Maas $W_2$ / gradient-flow branch | THEOREM under GNS detailed balance ([CM17] Thm 7.6; [Wir] in the tracial case). **NOT AVAILABLE** for KMS-only generators: only the necessary condition (BKM) is known [CM20 Thm 2.9], and the [CKG]/[CR] generators are not GNS (K1) | §3.5 |
| L9 | MLSI ⇒ TC$_2$ ⇒ TC$_1$ ⇒ Gaussian concentration; TC$_2$ ⇒ Poincaré; $\alpha_1\le\lambda$ | THEOREM | [RD] Thms 3, 4, 6, 8; [KT] Thm 16 |

### 3.1 L1–L2

[CKG] Thm I.1 and [DLL] Thm 10 produce KMS-symmetric generators for any $H$. Through L2 each is a quantum Dirichlet form on Hilbert–Schmidt space. Its similarity transform $\hat{\mathcal L}=\rho^{1/4}\mathcal L(\rho^{-1/4}\cdot\rho^{-1/4})\rho^{1/4}$ is Hermitian and $\le0$. This is the [KB] Table I correspondence "reversible Liouvillian ↔ frustration-free Hamiltonian, Gibbs state ↔ ground state $|\sqrt\rho\rangle$". It is the only structure used by the gap-free Cesàro argument (§4.3).

### 3.2 L3′: the [CR] Dirichlet form is the energy of an explicit twisted derivation

**Source** ([CR] Lemma X.1, citing [RFA24, Lemma C.2]). Write the transition part as $\sum_{\nu_1,\nu_2}\alpha_{\nu_1\nu_2}A^a_{\nu_1}[\cdot]A^{a\dagger}_{\nu_2}$ and set $h_{\nu_1\nu_2}=\alpha_{\nu_1\nu_2}e^{\beta(\nu_1+\nu_2)/4}$. Then
$$
\mathcal E(X,Y):=-\langle X,\mathcal L^\dagger[Y]\rangle_\rho=\sum_a\sum_{\nu_1,\nu_2}\bar\alpha_{\nu_1\nu_2}\,\mathrm{Tr}\big[\sqrt\rho[A^a_{\nu_1},X]^\dagger\sqrt\rho[A^a_{\nu_2},Y]\big],\qquad \bar\alpha_{\nu_1\nu_2}=\frac{h_{\nu_1\nu_2}}{2\cosh(\beta(\nu_1-\nu_2)/4)} .
$$
**[CR] Lemma X.2** (Gaussian weight):
$$\mathcal E(X,Y)=\sum_a\iint g(t)h^G(\omega)\mathrm{Tr}[\sqrt\rho[\hat A^a(\omega,t),X]^\dagger\sqrt\rho[\hat A^a(\omega,t),Y]]\,dt\,d\omega,$$
with $\hat A(\omega,t)=e^{iHt}\hat A(\omega)e^{-iHt}$ and $g(t)=\frac1{\beta\cosh(2\pi t/\beta)}$. The proof uses $\frac1{2\cosh(\beta x/4)}=\int g(t)e^{-ixt}dt$. **[CR] Lemma X.3** (Metropolis weight) replaces $h^G$ by $h(\omega)=e^{-\sigma^2\beta^2/8}e^{-|\omega|\beta/2}$. In the source's words, this is "an elegant, manifestly PSD quadratic form of commutators."

**Mine, in the [DLL] parametrisation** (checked numerically, K2). Take one self-adjoint coupling $A$ and a real even filter $q$. Then $L=\sum_\nu e^{-\beta\nu/4}q(\nu)A_\nu=\Delta^{1/4}\tilde A$ with $\tilde A=\sum_\nu q(\nu)A_\nu=\tilde A^\dagger$, and
$$
\mathcal E(X)=\sum_{\nu_1,\nu_2}\frac{q(\nu_1)q(\nu_2)}{2\cosh(\beta(\nu_1-\nu_2)/4)}\mathrm{Tr}\big[\sqrt\rho[A_{\nu_1},X]^\dagger\sqrt\rho[A_{\nu_2},X]\big]=\int_{\mathbb R}g(t)\,\big\|[\tilde A(t),X]\big\|_\rho^2\,dt,\qquad \tilde A(t)=e^{iHt}\tilde Ae^{-iHt},\quad \int g=\tfrac12 .
$$
The numerical agreement is $<10^{-13}$ relative error in both forms, over three random noncommuting 3-qubit Hamiltonians with $A=\{X_0,Y_0,Z_0\}$ and $\beta=1.3$.

**The bimodule and the derivation** (Mine; identification only).
- Bimodule: $\mathcal H_\partial=\bigoplus_aL^2(\mathbb R,g\,dt)\otimes L^2_{\rm KMS}(M_n,\rho)$.
- Derivation: $\partial X=\big(t\mapsto[\tilde A_a(t),X]\big)_a$, with $\mathcal E(X)=\|\partial X\|^2$.
- This is the continuum version of [VW] Thm 2.5 with $V_t=\sqrt{g(t)}\tilde A(t)$ self-adjoint.
- Transported to the standard form, $\delta(X)=\rho^{1/4}[V,X]\rho^{1/4}$ satisfies the twisted Leibniz rule $\delta(XY)=\sigma_{-i/4}(X)\delta(Y)+\delta(X)\sigma_{i/4}(Y)$, with $\sigma_{-i/4}(X)=\rho^{1/4}X\rho^{-1/4}$. This is [VW] Thm 2.4(iv). The quarter-twist *is* the KMS embedding of §2.1.

The kernel $g$ is the KMS-strip kernel; [Y] uses its relative $1/|\sinh\pi t|$ (sibling B14).

**What this upgrades.** The sibling [arxiv-2609.38007.md](arxiv-2609.38007.md) B9 ("the Dirichlet form is a squared Connes commutator") was labelled ANALOGY→SPECULATION. With [VW] Thm 2.5 and [CR] Lemma X.2 it is a **THEOREM**: the [CR] Dirichlet form is the KMS-$L^2$ energy of a derivation that is inner in each modular-time slice.

### 3.3 L4–L6: Dirac operator, Lip-norm, radius (Mine)

**The finite spectral triple.** Take inner derivations $\partial_j=[V_j,\cdot]$ with $\{V_j\}=\{V_j^*\}$. Set $\mathcal H=\mathbb C^n\otimes\mathbb C^{2N}$, $\pi(a)=\bigoplus_j(a\oplus a)$ and $D=\bigoplus_j\begin{pmatrix}0&V_j\\V_j^*&0\end{pmatrix}$. Then
$$\|[D,\pi(a)]\|=\max_j\max(\|[V_j,a]\|,\|[V_j^*,a]\|)=\max_j\|[V_j,a]\|.$$
The spectral-triple axioms are automatic in finite dimension. For the [CKG] sampler, use the direct integral over modular time $D=\int^\oplus\tilde A_a(t)\,dt$. Then $\|[D,\pi(a)]\|=\operatorname{ess\,sup}_t\|[\tilde A_a(t),a]\|=\sup_t\|[\tilde A_a,\alpha_{-t}(a)]\|$: the Lipschitz seminorm is the sup over modular time of commutators with the filtered coupling.

**$L^2$ energy against $L^\infty$ Lipschitz.** The Dirichlet form is the $\ell^2$-over-jumps, KMS-$L^2$ energy of the same derivation whose $\ell^\infty$, operator-norm size is $\|[D,a]\|$. Since $\|Y\|_{\rm KMS}\le\|Y\|_\infty$ ([CR] Lemma II.1),
$$\mathcal E(a)\le N\|[D,\pi(a)]\|^2,$$
or $\mathcal E(a)\le\tfrac12\sum_a\sup_t\|[\tilde A_a(t),a]\|^2$ in the continuous case (check K8).

**L5.** $\|[D,a]\|=0$ iff $a\in\{V_j\}'=\ker\mathcal L$, by [CM20] Prop 4.11 or [DLL] Lemma 3. By [Rie §7] the seminorm is therefore a Lip-norm, and $d_D$ is a metric giving the weak-\* topology, **iff the QMS is primitive**. Otherwise $d_D=\infty$ between states that differ on the fixed-point algebra. *The Connes metric is blind exactly to what the dynamics conserves.*

**L6.** Let $a$ be self-adjoint with $\rho(a)=0$. Then
$$\|a\|_\infty\le\rho_{\min}^{-1/2}\|a\|_{\rm KMS}\le(\lambda\rho_{\min})^{-1/2}\mathcal E(a)^{1/2}\le\sqrt{N/(\lambda\rho_{\min})}\,\|[D,a]\|.$$
The first step uses $\mathrm{Tr}(Y^*\rho^{1/2}Y\rho^{1/2})=\sum|Y_{ij}|^2\sqrt{p_ip_j}\ge p_{\min}\|Y\|_2^2$; the second is Poincaré. Hence $d_D(\rho,\psi)\le\sqrt{N/(\lambda\rho_{\min})}$ for every state $\psi$, with $N/2$ in place of $N$ in the continuous [CKG] case. The radius of $S(A)$ in the Connes metric is controlled by the spectral gap, the analogue of "expanders have small diameter". The $\rho_{\min}^{-1/2}$ loss is the same one that appears in the χ² mixing bound [KT Thm 22]. In examples the bound is crude (K8 values $10^3$–$10^4$); it is a qualitative bridge, not a sharp estimate.

**Infinite dimension.** For fractals, [CGIS] closes the loop as theorems: Dirichlet form → derivation → $D$ → (residue) → the same Dirichlet form (Thm 5.5), and the Connes metric reproduces the geodesic metric up to a power (Cor 5.3).

### 3.4 L7: at infinite temperature the single-Pauli Connes distance is the De Palma–Marvian–Trevisan–Lloyd $W_1$

**Mine (elementary; check K6).** For qubits, $E_iH=\tfrac14\sum_{P\in\{I,X,Y,Z\}}P_iHP_i$, so $H-E_iH=\tfrac14\sum_{P\in\{X,Y,Z\}}P_i[P_i,H]$. Write $M(H)=\max_{i,P}\|[P_i,H]\|$. Two bounds follow:
- $\|H-E_iH\|\le\tfrac34\max_P\|[P_i,H]\|$;
- since $[P_i,E_iH]=0$, $\|[P_i,H]\|\le2\|H-E_iH\|$.

So $\tfrac12M(H)\le\max_i\|H-E_iH\|\le\tfrac34M(H)$. With [DMTL] Prop 15 at $d=2$:
$$
\tfrac23M(H)\le\|H\|_L\le\tfrac32M(H),\qquad\text{hence}\qquad \tfrac23\,d_{D_P}\le W_1^{\rm DMTL}\le\tfrac32\,d_{D_P}.
$$
Here $D_P=\bigoplus_{i,P}P_i$ on $\mathbb C^{2^n}\otimes\mathbb C^{3n}$, and $M(H)=\|[D_P,\pi(H)]\|$. On 200 random 3-qubit $H$ the ratio $\max_i\|H-E_iH\|/M(H)$ ranged over $[0.521,0.683]$.

**Why this is the infinite-temperature chain.** At $\beta=0$ the [CR] generator is $\mathcal L_A(X)=\sum_{i\in A,P}(P_iXP_i-X)$. Its derivation is $([P_i,\cdot])$, so its Dirac operator is $D_P$. Hence:
- the Connes distance of the infinite-temperature single-Pauli Lindbladian is the [DMTL] quantum $W_1$, up to the factor $3/2$;
- on diagonal states it is the Hamming $W_1$ (Ornstein's $\bar d$). For diagonal $f$, $\|[X_i,f]\|=\max_x|f(x)-f(x\oplus e_i)|$ and $[Z_i,f]=0$, so $M(f)$ is the Hamming-Lipschitz constant;
- Marton's inequality and Gaussian concentration hold for it ([DMTL] Thms 2–3).

At $\beta>0$ the [CKG] dressing replaces $P_i$ by $\tilde P_i(t)$. I have not seen a $W_1$ built from the dressed Paulis, nor its transport-entropy inequality for Gibbs states, in the sources read.

### 3.5 L8: the Wasserstein-2 branch is missing at $\beta>0$

[CM17] and [CM20] build $W_{2,\mathcal L}$ only for GNS detailed balance. [CM20] Thm 2.9 shows that any gradient-flow structure forces BKM symmetry, and the BKM class strictly contains the GNS class; the gap between the two is open. Check K1 shows that the [DLL]/[CKG] generators used by [CR] do not commute with $\Delta_\rho$ (relative size of $[\mathcal L,\Delta]$ of order $10^3$). By [CM17] Thm 2.9 they are therefore not GNS-symmetric. **No Carlen–Maas metric is known for the time-averaged single-Pauli Lindbladian at $\beta>0$.** Only the $W_1$/Connes branch of the chain exists there. The GNS alternative, Davies, has a $W_2$ but is non-local ([DLL] Remark 7; [CKG] App. E).

---

## 4. The phrase, unpacked, and its role in proving approximate Markov properties

### 4.1 The exact objects (Source: [CR] §II, (1.1), (3.1))

$$\mathcal R_{A,t}[\cdot]=\frac1t\int_0^t\exp(s\mathcal L_A)[\cdot]\,ds,\qquad \mathcal L_A=\sum_{a\in P^1_A}\mathcal L_a .$$
- Each $\mathcal L_a$ is (2.1) of [CKG] with the single jump $A^a$, the Metropolis weight and energy width $\sigma=1/\beta$.
- The comparison state is $\rho_{\beta,-A}=\mathrm{Tr}_A\rho_\beta\otimes\tau_A$, with $\tau_A$ maximally mixed.
- Implementation (Remark II.2.1): sample $s$ uniformly in $[0,t]$ and run $e^{s\mathcal L}$.

### 4.2 "time-averaged": three different time averages (Mine; ingredients Source)

1. **T1, Davies / secular.** $A_\nu=\lim_{T\to\infty}\frac1{2T}\int_{-T}^Te^{-i\nu t}e^{iHt}Ae^{-iHt}dt$. *Mine, elementary from $e^{iHt}Ae^{-iHt}=\sum_\nu A_\nu e^{i\nu t}$.* This is the mean ergodic projection of the unitary group $e^{-i\nu t}\mathrm{Ad}\,e^{iHt}$ onto its fixed space, so it projects onto spectral subspaces of the **modular derivation** $\delta_H=i[H,\cdot]$. It needs times of order (Bohr-frequency spacing)$^{-1}$, which can be exponential ([DLL] Remark 7: the energy–time uncertainty). It yields exact $\Delta$-eigen jumps, hence GNS detailed balance (Alicki), and non-locality.
2. **T2, [CKG] Gaussian window.** $\hat A(\omega)=\frac1{\sqrt{2\pi}}\int e^{iHt}Ae^{-iHt}e^{-i\omega t}f(t)dt$ with $f\propto e^{-t^2/\beta^2}$: T1 smoothed over a time of order $\beta$. It is quasi-local by Lieb–Robinson ([DLL] Prop 20). The jumps are only approximate $\Delta$-eigenoperators, so GNS detailed balance fails. Exact KMS detailed balance is restored by the coherent term ([CKG] Cor II.2: $B_\nu$ with $\tanh(-\beta(\nu_1-\nu_2)/4)$ weights).
3. **T3, [CR] Cesàro mean of the semigroup.** $\mathcal R_t=\frac1t\int_0^te^{s\mathcal L}ds$ is the mean ergodic projection of the **dissipative semigroup** onto $\ker\mathcal L=\ker\partial$, reached at a gap-free rate (§4.3).

**Structure (Mine).** T1 and T3 are both von Neumann mean-ergodic projections. T1 belongs to the *modular* derivation $\delta_H$, the "time" generator. T3 belongs to the *dissipative* derivation $\partial$, the "transport" generator. The detailed-balance conditions are compatibility conditions between the two:
- **GNS:** the components of $\partial$ are $\delta_H$-homogeneous, $[V_j,H]=-\omega_jV_j$ ([CM20] Thm 2.4).
- **KMS:** weaker. Each jump is real for the quarter-twisted involution, $\Delta^{-1/4}L_j=(\Delta^{-1/4}L_j)^\dagger$ ([DLL] Thm 10), and the derivation is twisted by $\sigma_{\mp i/4}$ ([VW]).

[CMo]'s twisted spectral triples twist by the modular group at $t=-i$, and [CNNR]'s modular spectral triples take the generator of the modular circle action as the Dirac operator. On the NCG side, the "time" derivation is therefore already the Dirac operator in type III index theory (§6.3). That the dissipative and modular derivations are two Dirac-type operators with detailed balance as their compatibility condition is my reading, labelled ANALOGY.

### 4.3 "detailed-balanced": why KMS, and what it buys

- **Why KMS and not GNS.** GNS detailed balance with local couplings is impossible for noncommuting $H$ without T1, by [CKG] App. E and Alicki. KMS is the strongest symmetry compatible with quasi-local filtered jumps.
- **What it buys.** A symmetric Dirichlet form, and with it the spectral theorem. For any $X$, with $x$ the spectral variable of $-\mathcal L$ in the KMS geometry,
  $$\mathcal E(\mathcal R_tX)=\int x\Big(\frac{1-e^{-tx}}{tx}\Big)^2d\mu_X(x)\le\frac{c}{t}\|X\|^2_{\rm KMS},\qquad c=\sup_{u>0}\frac{(1-e^{-u})^2}{u}\approx0.4073 .$$
  This derivation is in [expanders.md](expanders.md) §9.9 and is checked here at $t=1,10,10^2,10^3$ (K4). [CR] Cor VII.1 states the operator-norm version, $\mathcal E(\mathcal R^\dagger_t[O])\le2/t$ for $\|O\|\le1$. The source: "such a property is generally false for the Lindblad evolution $e^{\mathcal L^\dagger_At}$ itself without time-averaging, because $\mathcal L_A^\dagger$ may have arbitrarily small eigenvalues. Thus, time-averaging provides a different mechanism to obtain a small Dirichlet form that is independent of the spectral gap."
- **Classical analogue.** "Locally stationary distributions" ([LMRRW] Lemma 3.1; sibling B3).

### 4.4 "single-Pauli": a generating set, the noncommutative hypercube, and what happens at $\beta>0$

**What single-site Paulis provide.**
- Products of single-site Paulis give every Pauli string on $A$, so the commutant of $P^1_A$ is $1_A\otimes M_{A^c}$.
- Passing from single-site to string commutators costs the word length ([CR] Lemma VIII.1 is the Leibniz expansion $[\prod_iA_i,O]=\sum_j\prod_{i<j}A_i[A_j,O]\prod_{i>j}A_i$; Cor VIII.1 bounds it with $2^w$ and a Hölder exponent).
- High-weight jumps would be exponentially costly to simulate ([CR] Remark VIII.0.1).

**At $\beta=0$: the noncommutative hypercube** (Mine; check K5).
- The generator is $\mathcal L_A(X)=\sum_{i\in A,P}(P_iXP_i-X)$. On a Pauli string $S$, $\mathcal L_A(S)=-4\,w_A(S)\,S$, with $w_A$ the weight on $A$.
- Hence the kernel is $1_A\otimes M_{A^c}$ and the gap is $4$, independent of $|A|$.
- The Efron–Stein inequality $\|X-E_AX\|_2^2\le\sum_{i\in A}\|X-E_iX\|_2^2$ holds (max ratio 0.713 in K5).
- Harrow's mechanism [Har], read for the Pauli group: $M_{2^{|A|}}$ splits into the $4^{|A|}$ characters labelled by Pauli strings, and the conjugation walk is a Cayley walk on $\mathbb Z_2^{2|A|}$.
- As a channel, $\frac1{3|A|}\sum P_i\cdot P_i$ has second eigenvalue $1-\frac4{3|A|}$. This is **not** a quantum expander family: the degree $3|A|$ is unbounded and the normalised gap tends to $0$. Its dimension-free continuous-time gap is a *tensorization* phenomenon, not an expansion phenomenon.

**At $\beta>0$: the exact vertical structure is destroyed** (Mine; check K3). The [DLL] single-Pauli-on-$A$ generator ($n=3$, $A=\{0\}$, $\beta=1.3$) has kernel:
- dimension **1** for random noncommuting $H$ (three trials), against $4^{n-1}=16$ at $\beta=0$;
- dimension **8** for a commuting Ising $H$, consisting of operators diagonal on the boundary qubit, tensored with anything on the rest.

The commuting case is the classical Markov property seen as a fixed-point algebra: conditioning on the boundary configuration ([KB] Lemma 12, "locally primitive"). In the noncommuting case no exact "fibre over $A^c$" survives. [CR] therefore never use the kernel. They use approximate commutation with $P^1_A$ at finite $t$, together with quasi-locality (§4.6). This is the dynamical side of the sibling's B10. By Takesaki's theorem (*from memory*), no $\rho$-preserving conditional expectation onto $M_{A^c}$ exists, because $M_{A^c}$ is not modular-invariant.

### 4.5 "on A": locality (Source)

- [CR] Lemma VII.3: $\|\mathcal L^\dagger_{A,\ell}-\mathcal L^\dagger_A\|_{\infty\to\infty}\lesssim|A|(e^{-c'\ell/(d\beta)}+2^{-\ell})$ for $\ell\ge4e^2\beta d$.
- Lemma VII.2: truncation errors grow at most linearly in $t$, "Thus, the quasi-locality holds for exponential times."
- Every $\mathcal L_a$ fixes $\rho_\beta$, so $\mathcal R_{A,t}[\rho_\beta]=\rho_\beta$.

### 4.6 The proof of [CR] Thm III.1 in derivation language (Mine, following [CR] §XI.A)

1. **Twirl.** $\rho-\rho_{-A}=2^{-2|A|-1}\sum_{S\in P_A}[S,[S,\rho]]$. Hence $|\mathrm{Tr}[X\mathcal R(\rho-\rho_{-A})]|\le2^{-2|A|-1}\sum_S\|[S,[S,\mathcal R^\dagger X]]\|_\rho$, using stationarity and KMS Cauchy–Schwarz. In NCG terms the right-hand side is a second-order quantity in the *string* derivations $[S,\cdot]$.
2. **Peeling the outer commutator** (Lemma IX.5 with the Gibbs-conjugation bound $\le2^{|A|}$ at $\beta_0=1/4d$). This reduces to first order.
3. **Strings to single sites** (Cor VIII.1): $2^w$ and a Hölder exponent $16\beta_0^2/\beta^2$.
4. **Commutator from the derivation** (Lemma X.4):
   $$\|[A,O]\|_\rho\lesssim d^2|A|\,(\cdots)\,\mathcal E(O)^{2\beta_0/(\beta+5\beta_0)} .$$
   The proof smears over Heisenberg time $|t|\le\epsilon$, cuts frequencies at $\Omega$ (because $1/h(\omega)=e^{\sigma^2\beta^2/8}e^{|\omega|\beta/2}$ diverges), and applies Cauchy–Schwarz against $g\,h$. In derivation language: **a single component of the derivation, evaluated at modular time $0$, is controlled by the total $L^2$ energy, at a Hölder cost.** That is a regularity statement about the derivation of §3.2.
5. **Cesàro** (Cor VII.1): $\mathcal E_A(\mathcal R^\dagger_{A,t}X)\le2/t$.
6. **Chain:** $|A|^2\,2^{2|A|}\,r\,t^{-\lambda}$ with $\lambda=\frac{128\beta_0^4}{\beta^3(\beta+5\beta_0)}$ at low temperature (exponent as in Thm III.1).

**What the argument is (Mine).** Steps 1–4 prove a **relative Poincaré-type inequality**: the distance of an observable from $1_A\otimes M_{A^c}$, measured through twirl commutators, is bounded by a power of the energy of the single-site derivation $\partial_A$. Note that $1_A\otimes M_{A^c}$ is not $\ker\partial_A$ at $\beta>0$ (K3). At $\beta=0$ this is Efron–Stein with exponent 1 and a polynomial constant (K5). At $\beta>0$, [CR] get exponent $\frac{128\beta_0^4}{\beta^3(\beta+5\beta_0)}\ll1$ and constant $2^{O(|A|)}$.

### 4.7 Role in approximate Markov properties, and where expansion enters (Source + Mine)

- **Local Markov property.** $\mathcal R_{A,t}\circ(\tau_A\otimes\mathrm{Tr}_A)$ is a quasi-local approximate recovery map. With continuity, $I(A{:}C|B)\lesssim\log(\dim C)\sqrt\Delta$ ([CR] Cor III.2): CMI $\lesssim r'|A||C|\exp(\mu'\min(|A|,|C|)-\lambda'\mathrm{dist}(A,C))$.
- **Global Markov property under a gap.** [CR] Def B.1: $\mathcal L$ is $\lambda$-locally gapped if $-\lambda\langle X,\mathcal L^\dagger X\rangle_\rho\le\langle X,\mathcal L^{\dagger2}X\rangle_\rho$, uniformly for every restricted Gibbs state $\rho_X$. Cor B.2: a uniform local gap $c|A|^{-c'}$ implies CMI $\le\mathrm{Poly}(|A|,|C|)e^{-\mathrm{dist}(A,C)/\xi}$. Remark III.1.2: "If the present argument can be combined with a faster mixing time or spectral gap analysis, one might be able to improve the exponential dependence on |A|, hence establishing the global Markov property."
- **Where expansion enters (Mine).** A uniform local gap is uniform coercivity of the vertical derivation $\partial_A$ off its kernel. In [Pet]'s language (§5.4) it is a *uniform Kazhdan constant* for the single-site generating set, acting on KMS-twisted bimodules. This is the one precise place where the expander essence enters the Markov theory. The Markov property itself needs no gap: [CR] use the Cesàro mean, [Y] uses the modular cocycle.

### 4.8 The abstract template (Mine)

Data: an inclusion $N\subset M$ of finite-dimensional von Neumann algebras, a faithful state $\varphi$, and a generating set $\mathcal G$ of the relative commutant $N'\cap M$ (the "fibre"). Build the KMS-detailed-balanced Lindbladian with couplings $\mathcal G$ ([DLL] Thm 10) and Cesàro-average it.
- **If $N$ is $\sigma^\varphi$-invariant** and $\mathcal G$ can be chosen of $\Delta$-eigenvectors: the generator is GNS, its kernel is $\mathcal G'$, and the Cesàro limit is the Takesaki conditional expectation onto $\mathcal G'$, which is $N$ when $(N'\cap M)'\cap M=N$. Then the Jenčová–Petz sufficiency test (Note 1 §6, condition iii) is run **dynamically**: a family is sufficient for $N$ iff each member is fixed by the Schrödinger dual of the limit. §6.1 is this case for the programme, checked numerically.
- **If $N$ is not modular-invariant** (noncommuting Gibbs states): the kernel collapses (K3), no exact expectation exists, and the finite-$t$ Cesàro map is an approximate, quasi-local recovery. That is [CR].

The construction depends only on the inclusion and the state, so it passes the drag test. It reads the same for a Bratteli diagram, a tiling AF algebra or a measurement context.

---

## 5. Spectral gap, log-Sobolev, mixing and quantum expanders: what is proved

### 5.1 Functional inequalities (Source)

- Spectral-gap and LSI mixing bounds: [KT] Thm 22, §2.8.
- $\alpha_1\le\lambda$: [KT] Thm 16.
- Chain MLSI ⇒ TC$_2$ ⇒ TC$_1$ ⇒ Gaussian concentration, and TC$_2$ ⇒ PI ⇒ exponential concentration: [RD] Figure 2, Thms 3, 4, 6, 7, 8.
- [CM17] (abstract): uniform convexity of relative entropy along the $W_2$ geometry (a Ricci bound) gives "new inequalities for the decay of relative entropy". The exact constants were not retrieved.
- Complete MLSI for finite-dimensional GNS-symmetric QMS: [GR22], cited in [VW §1] (not opened).

### 5.2 Gibbs samplers (Source)

- **Commuting Hamiltonians.** Gap ⇔ strong clustering; gapped in 1D and at high temperature ([KB] Thms 23, 26, 30, 31).
- **Noncommuting Hamiltonians.** [RFA24] "estimated the spectral gap of the efficient quantum Gibbs sampler in [CKG23] in the high-temperature regime, by mapping the Lindbladian … to a Hamiltonian … and then analyzing its spectral properties by the stability of gapped Hamiltonians" ([DLL] Remark 21; not opened). Bakshi–Liu–Moitra–Tang: Dobrushin-type high-temperature rapid mixing (sibling digest; abstract only).
- **General temperature, noncommuting.** [CR] App. B: "For general noncommutative Hamiltonians (with local jumps), we do not know of any a priori bound on the local gap, even assuming high temperature" (sibling §8.6, re-verified).

### 5.3 Quantum expanders against Gibbs samplers (Source + Mine)

[BST], [Har] and Hastings give quantum expanders: constant degree $D$, gap $\Omega(1)$ on $M_N$. [KT] (144) shows such maps have $\alpha_2=O(\log D\cdot\log\log N/\log N)$, so their mixing time is $\Omega(\log N)$.

*Mine.* A local Gibbs sampler on $n$ qubits has $3n=O(\log N)$ jump families on $M_N$ with $N=2^n$. Its degree grows with $\log N$, and when it is gapped the gap is a *tensorization* (product-structure) gap, not an expander gap. The expander notion relevant to Gibbs sampling is the high-dimensional one: link spectra, spectral independence and local-to-global, as in [hdx-spectral-independence.md](hdx-spectral-independence.md). It is not the bounded-degree quantum expander. The bounded-degree object would matter only for a *design* question: sampling a *maximally mixed* fibre state with few Kraus operators.

### 5.4 Expansion as a property of derivations: (T), amenability and Haagerup through Dirichlet forms (Source)

- **[Pet] Thm 0.1** (vN-algebra Delorme–Guichardet). For a separable finite factor $N$ the following are equivalent:
  1. $N$ has property (T);
  2. $N$ does not have property Γ, and for any weakly dense $*$-subalgebra $N_0\ni1$ containing a non-Γ set, every densely defined closable derivation on $N_0$ into a Hilbert $N$–$N$ bimodule is inner;
  3. there is a weakly dense $*$-subalgebra $N_0$, countably generated as a vector space, such that every closable derivation into a Hilbert $N$–$N$ bimodule whose domain contains $N_0$ is inner.
- **[Pet] Thm 3.2.** Condition (b), $\exists F,K$ such that $\|\xi_0-\xi\|\le K\max_{x\in F}\|x\xi-\xi x\|$ for some central vector $\xi_0$, is the Kazhdan-pair form of (T). Condition (d): every closable conservative symmetric c.c.n. map is $\|\cdot\|_1$-bounded. The source notes these "are in fact extensions of generators associated to completely Dirichlet forms".
- **[CS17] Thm 3.15.** A Dirichlet form with subexponential spectral growth, i.e. $\mathrm{Tr}\,e^{-tL}<\infty$ for all $t>0$, forces $N$ to be amenable. Remark 3.16: on a non-amenable $N$ every Dirichlet form has exponential spectral growth. Example 3.11(ii): for a negative-definite word length, the spectral growth rate of $\mathcal E_\ell$ is the growth rate of $(\Gamma,S)$.
- **[CS17] §1.** On $L(\Gamma)$, $\mathcal E_\ell[a]=\sum_g\ell(g)|a(g)|^2$ is a Dirichlet form iff $\ell$ is conditionally negative definite, with discrete spectrum iff $\ell$ is proper. Caspers–Skalski: Haagerup property (H) ⇔ existence of a Dirichlet form with discrete spectrum.
- **[Voi] §14** (stated without proof there). For a finitely generated group with generating tuple $\gamma$, the classical $J$-capacity of condensers on the Cayley graph equals the noncommutative quasicentral condenser modulus of $\lambda(\gamma)$: $\mathrm{cap}_J(X_1,X_2)=k_J(\lambda(\gamma);P_{\ell^2(X_1)},P_{\ell^2(X_2)})$. Moreover $\mathrm{cap}_p(G)>0$ or $=0$ is $p$-hyperbolicity or $p$-parabolicity (Yamasaki).

**Reading (Mine).** Every property in the expander family can be phrased through derivations or Dirichlet forms:
- **(T):** derivations are inner, with a uniform Kazhdan constant.
- **non-amenability:** spectral growth is exponential.
- **(H):** discrete spectrum.
- **expansion of one graph:** coercivity of its graph derivation $f\mapsto f(r\gamma)-f(s\gamma)$, which is $[b,f](\gamma)$ in the groupoid algebra (Note 1 §3.3).

This is the operator-algebraic common ground the user senses. The ghost-projection face of it is in [expanders.md](expanders.md) §9.11.

---

## 6. Bridges to the programme

### 6.1 B-NC1, THEOREM (elementary, Mine; check K7): the re-routing Lindbladian on the complete-history corner

**Setting** (Note 1 §5.2). Take one origin block of the complete-history corner, $M_N$ with $N$ the number of complete histories, and the endpoint map $r$. The state is $\mathbb P_\beta=\mathrm{diag}(e^{-\beta F(\mu)})/Z$. The jumps are the re-routing matrix units $V_{\mu\nu}=e_{\mu\nu}$, for $\mu\ne\nu$ with $r(\mu)=r(\nu)$; in Note 1's notation these are $s_\mu s_\nu^*$ restricted to complete histories. Since $\Delta_{\mathbb P_\beta}e_{\mu\nu}=e^{-\beta(F(\mu)-F(\nu))}e_{\mu\nu}$, set $\omega_{\mu\nu}=\beta(F(\mu)-F(\nu))$ and
$$
\mathcal L_\beta(X)=\sum_{\mu\neq\nu,\ r\mu=r\nu}c_{\{\mu\nu\}}\,e^{-\omega_{\mu\nu}/2}\big(e_{\mu\nu}^*[X,e_{\mu\nu}]+[e^*_{\mu\nu},X]e_{\mu\nu}\big),\qquad c_{\{\mu\nu\}}=c_{\{\nu\mu\}}>0 .
$$

**Claims.**
- (a) $\mathcal L_\beta$ satisfies GNS (hence KMS) detailed balance for $\mathbb P_\beta$. This is Alicki's converse, [CM17] Thm 3.1. Check K7: KMS asymmetry $10^{-16}$, $[\mathcal L,\Delta]=4\cdot10^{-16}$.
- (b) $\ker\mathcal L_\beta=\{e_{\mu\nu}\}'=\mathrm{span}\{p_\tau\}$, the endpoint face algebra $D_0$ restricted to the corner. Reason: the fibre matrix units generate $\bigoplus_\tau M_{h_\tau}$, whose commutant in $M_N$ is $\bigoplus_\tau\mathbb C1_\tau$. K7 finds dimension 3 for 3 endpoints.
- (c) $\lim_t\frac1t\int_0^te^{s\mathcal L_\beta}ds=E$, with $E(X)=\sum_\tau\frac{\mathbb P_\beta(p_\tau X)}{\mathbb P_\beta(p_\tau)}p_\tau$ the $\mathbb P_\beta$-preserving (Takesaki) conditional expectation onto $D_0$. K7 error: $10^{-15}$.
- (d) $\mathcal E(X)=\sum c\,\|[e_{\mu\nu},X]\|^2_{\rm KMS}$, by [CM20] Prop 2.5. K7 error: $2\cdot10^{-16}$.
- (e) For $\beta'\ne\beta$ the following are equivalent:
  - $\mathbb P_{\beta'}$ is stationary for $\mathcal L_\beta$;
  - $E_*\mathbb P_{\beta'}=\mathbb P_{\beta'}$;
  - $F$ is constant on every endpoint fibre;
  - with a single origin, the face is a sufficient statistic (Note 1 §6), equivalently $F=\delta U$ with $U$ constant on the input layer.

  K7: for random $F$, $\|\mathcal L^*\mathbb P_{\beta'}\|_1=0.58$; for $F=U\circ r$ it is $10^{-16}$.
- (f) **Quantitative.**
  $$\|E_*\mathbb P_{\beta'}-\mathbb P_{\beta'}\|_1\le2\tanh\big(|\beta-\beta'|\max_\tau\mathrm{osc}_\tau F/4\big).$$
  K7: $0.080\le0.487$.

**Proof of (f).** $E_*\mathbb P_{\beta'}$ keeps the endpoint marginal of $\mathbb P_{\beta'}$ and replaces each fibre law by the fibre Gibbs law at $\beta$. Hence $\|E_*\mathbb P_{\beta'}-\mathbb P_{\beta'}\|_1=\sum_\tau\mathbb P_{\beta'}(p_\tau)\|G^\tau_\beta-G^\tau_{\beta'}\|_1$. On a fibre, $dG^\tau_{\beta'}/dG^\tau_\beta$ takes values in an interval of multiplicative width $e^{\delta}$ with $\delta=|\beta-\beta'|\mathrm{osc}_\tau F$. For such a ratio $r$ with mean $1$, $\mathrm{TV}=\mathbb E(r-1)_+\le\frac{(M-1)(1-m)}{M-m}$, which over $M/m=e^\delta$ is maximised at $\tanh(\delta/4)$. ∎

**Dictionary with [CR]** (the analogy the user sensed, now exact):

| [CR] (noncommuting lattice Gibbs state) | programme (Note 1 history corner) |
|---|---|
| region $A$; forget $A$: $\rho_{-A}=\mathrm{Tr}_A\rho\otimes\tau_A$ | endpoint fibre; forget the history given the face |
| single-Pauli couplings generating $M_A$ | re-routing matrix units generating $\bigoplus_\tau M_{h_\tau}$ (or local flips, §6.5) |
| KMS-detailed balance, made exact by a coherent term | GNS-detailed balance, exact and as local as the jumps, since $H_F$ is diagonal in the history basis (the commuting case; sibling B10) |
| Cesàro mean $\mathcal R_{A,t}$, gap-free | Cesàro mean, converging to the Takesaki expectation onto $D_0$ |
| approximate recovery; CMI decays with distance | exact recovery iff sufficiency; defect $\le2\tanh(\Delta\beta\,\mathrm{osc}F/4)$ |
| $2^{O(|A|)}$ from word length in the Pauli group | diameter of the re-routing (flip) graph in a fibre (§6.5) |
| uniform local gap ⇒ global Markov (Cor B.2) | uniform fibre gap ⇒ Poincaré bound on the sufficiency defect ([expanders.md](expanders.md) §11.1) |

So the time-averaged detailed-balanced Lindbladian is, in the programme, the dynamical implementation of Petz–Jenčová sufficiency. Its Cesàro limit is the conditional expectation onto the barycentre (face) algebra. [CR]'s difficulties all come from the absence of that exact expectation, and in Note 1's KMS family it is present. **Guard:** the realised history laws of [mlp-bridge.md](../../mlp-bridge.md) §3.2 are outside the Gibbs class, because they have memory. The construction applies to the Markov part, i.e. the fitted kernel, and to finite-range extensions (sibling B11).

### 6.2 B-NC2, THEOREM (elementary, Mine; not numerically checked): the arrow-jump Lindbladian and Note 1's Dirac operator

Note 1 §5 gives $\alpha^F_t=\mathrm{Ad}\,e^{itH}$ on the Toeplitz algebra $\mathcal T(\Lambda)=\bigoplus_\sigma M_{N(\sigma)}$ acting on $\ell^2(E^*)$, with $Hh_\lambda=F(\lambda)h_\lambda$. Let $\omega=\bigoplus_\sigma w_\sigma e^{-\beta H}|_\sigma/Z_\sigma$ with all $w_\sigma>0$; it is a faithful KMS$_\beta$ state. Then $\Delta_\omega s_\gamma=e^{-\beta F(\gamma)}s_\gamma$. By Alicki's converse, the jumps $\{s_\gamma,s_\gamma^*\}$ define a GNS-detailed-balanced QMS with $\omega_\gamma=\beta F(\gamma)$ and Dirichlet form $\sum_\gamma c_\gamma(\|[s_\gamma,X]\|_{\rm KMS}^2+\|[s^*_\gamma,X]\|_{\rm KMS}^2)$.
- **Kernel.** Each $s_\gamma$ preserves the origin blocks. Within a block the restrictions generate the whole block, because $p_\sigma|_{\rm block}$ is the projection onto the trivial path. So the kernel is the **centre**, one scalar per origin: this "arrow walk" equilibrates everything except the origin face.
- **Dirac operator.** Its $D$ (§3.3) is the direct-sum form of Note 1 §8's $D=\sum_\gamma\ell(\gamma)^{-1}(s_\gamma+s_\gamma^*)$. Its Connes distance is a metric on states with fixed origin weights, and by L6 its radius is controlled by the gap of the noncommutative quiver Laplacian.

This places Note 1 §8's optional metric layer inside the chain: it is the Connes metric of the arrow-Lindbladian's derivation.

### 6.3 B-NC3, KNOWN-LINK + SPECULATION: modular spectral triples as the NCG home of Note 1's KMS structure

- **KNOWN-LINK.** [CNNR] build index theory from exactly the data of Note 1 §5: a C\*-algebra with a circle action and a KMS state. $D$ is the generator of the action, the fixed-point algebra carries the trace $\phi_D$, and the pairing is spectral flow, with Araki relative entropy in examples. The graph-algebra instance is Pask–Rennie (cited there).
- **Mine.**
  - For Note 1's plain gauge action ($F\equiv1$) the fixed-point algebra is the AF core $\mathrm{span}\{s_\mu s_\nu^*:|\mu|=|\nu|\}$.
  - On the complete-history corner the plain gauge action is trivial, since all complete histories have length $L$. The content is in a weighted $F$.
  - Integer-valued $F$ gives a circle action, to which [CNNR] applies verbatim. Real $F$ gives an $\mathbb R$-action, not covered by [CNNR].
- **SPECULATION.** U5's "gap labelling of the sufficiency defect" might be a [CNNR] modular index pairing. [CNNR] identify $sf_{\phi_D}(D,uDu^*)$ with Araki relative entropy, and Note 1 §6 identifies insufficiency with a non-trivial modular cocycle. What must be true: a stationary Bratteli diagram (U5), an integer-valued action cocycle, and modular partial isometries implementing re-routings.

### 6.4 B-NC4, THEOREM (elementary, Mine) + ANALOGY with (T): sufficiency is innerness of the modular derivation by a face observable

On the corner with one origin, $\delta_F=i[H_F,\cdot]$ is the generator of $\alpha^F$. Every derivation of $M_N$ is inner, so the question is *by what*.

**Claim.** $\delta_F$ is inner by an element of $D_0+\mathbb C1$ iff $F(\mu)=U(r\mu)+c$ on complete histories, i.e. iff the face is sufficient. Moreover
$$\inf_{U,c}\|H_F-U\circ r-c\|_\infty=\tfrac12\max_\tau\mathrm{osc}_\tau F,$$
which is the quantity in sibling B8 and in §6.1(f). Elementary, by the sup-norm approximation of $F$ by functions of the endpoint.

**ANALOGY.** [Pet]'s (T) is "every closable derivation is inner, uniformly". The programme's sufficiency is a *relative* innerness: inner by the subalgebra $D_0$. Its defect is the distance above. The relative notion corresponds to rigidity of inclusions (Popa, cited in [Pet]), a direction not explored here.

**Guard.** $H_F$ as a Dirac operator has $\ker[H_F,\cdot]\supseteq$ all diagonal elements, so its Connes distance is infinite between distinct measures on histories. It is a "time" operator, as in [CNNR], not a metric one. The metric comes from transport derivations (§6.1, §6.2).

### 6.5 B-NC5, SPECULATION: local flips as the analogue of single-site Paulis, and where U3 meets [CR]

Replace all re-routings $e_{\mu\nu}$ by **flips**: $e_{\mu\nu}$ for $\mu,\nu$ differing in one intermediate face. These are the squares of U3 and of [expanders.md](expanders.md) §11.1. They generate the fibre algebra iff the flip graph of each fibre is connected.
- Commutators with long re-routings are then Leibniz-expanded along flip paths, which is the analogue of [CR]'s $2^w$, with $w$ the flip distance.
- The analogue of a uniform local gap is the gap of the $\mathbb P_\beta$-weighted flip Laplacian on the fibres. U3's coboundary-expansion constant would bound the sufficiency defect by square defects $\delta F(\square)$.

What must be true:
1. connectivity of fibre flip graphs (a property of the quiver);
2. a gap that is uniform in depth for a stationary diagram (U5);
3. for the noncommuting extension (sibling B13), a depth-locality bound replacing Lieb–Robinson.

### 6.6 B-NC6, KNOWN-LINK: transport-entropy inequalities on faces (for prong 2)

Restricted to diagonal states, the β=0 single-Pauli chain gives Ornstein's $\bar d$ with Marton's inequality ([DMTL] Thm 2, product reference) and Gaussian concentration (Thm 3). Faces are $\{0,1\}$-patterns on a layer's units. Under any law on faces close to a product law in relative entropy, Hamming-Lipschitz statistics of faces concentrate. One example is the face dimension $d(\sigma)$, which is 1-Lipschitz. This is a classical consequence. It is stated here because its noncommutative extension is the chain of §3.4, and because it gives prong 2 a calibrated test: the measured spread of $d(\sigma)$ against $\sqrt{n\,S(\cdot\|\text{product})/2}$.

### 6.7 Guards

- **KMS twist.** For non-tracial states, untwisted Hilbert-bimodule derivations can fail to exist [Ver]. Every derivation statement at $\beta>0$ carries $\sigma_{\mp i/4}$ or a KMS-weighted norm.
- **Missing Wasserstein-2.** No Carlen–Maas gradient-flow metric is known for KMS-only generators (§3.5). Statements about Ricci curvature or Talagrand TC$_2$ for the [CR] dynamics are unsupported.
- **Finite-dimensional NCG.** Spectral-triple axioms are trivial in finite dimension, so the content is in constants: $\rho_{\min}$, $N$, gaps, Hölder exponents. [CGIS]-type theorems need infinite dimension or a family.
- **Connes metric degeneracy.** It is infinite between states differing on $\ker\partial$ (L5). For the re-routing dynamics, $d_D=\infty$ between states with different face marginals: the metric only sees within-fibre (history) differences. That is the right behaviour for measuring insufficiency, and the wrong one for comparing faces.
- **Single-Pauli at $\beta=0$ is not an expander** (§4.4). Its gap is a tensorization gap. Bounded-degree quantum expanders are a different object (§5.3).
- **Commuting vs noncommuting.** Note 1's KMS family is the commuting case. §6.1 is exact there, and any import of [CR]'s *estimates* (rather than its architecture) is vacuous. The [CR] mechanism matters for the programme only after coherent, off-diagonal "energies" are introduced (sibling B13).
- **Memory.** The realised history law is not Gibbs ([mlp-bridge.md](../../mlp-bridge.md) §3.2). B-NC1 applies to Gibbs laws and their finite-range extensions, not to the raw empirical law.

---

## 7. Explicit list of known theorems connecting the areas

| link | theorem | source (as read) | status |
|---|---|---|---|
| Markov semigroups ↔ NC Dirichlet forms | KMS-symmetric QMS ↔ quantum Dirichlet forms on the standard form | [VW §5] citing Cip97 Thm 4.11, GL95 Thm 5.7; [CS17] Def 2.1, (2.4) | Fact (via secondary) |
| Markov semigroups ↔ NC Dirichlet forms | finite dimension: self-adjoint QMS ⇔ conservative completely Dirichlet forms; Beurling–Deny | [CM20] §2.2, §3 | Fact |
| Dirichlet forms ↔ derivations | tracial: $L=\delta^*\delta$ for an essentially unique closable derivation into a Hilbert bimodule | Cipriani–Sauvageot 2003 via [CS17 §2(e)], [Ver] | Fact (via secondary) |
| Dirichlet forms ↔ derivations | KMS: $\mathcal L=\delta^*\delta$, $\delta$ twisted by $\sigma_{\mp i/4}$; finite dimension: $\sum_j\|[V_j,\cdot]\|^2_{\rm KMS}$ | [VW] Thms 2.4, 2.5, 4.2, 5.4 | Fact |
| Dirichlet forms ↔ derivations | untwisted derivation can fail to exist for non-tracial states | [Ver] Ex. 5.3 | Fact |
| Gibbs samplers ↔ derivations | Dirichlet form of the [CKG] sampler = KMS-strip-weighted squared commutators with dressed jumps | [CR] Lemmas X.1–X.3 ([RFA24] C.2); §3.2 | Fact + Derivation (K2) |
| Gibbs samplers ↔ detailed balance | GNS-DB ⇔ Alicki form with Δ-eigen jumps; KMS-DB ⇔ [DLL] Thm 10 form; $s\neq1/2$ forces Davies (non-local) | [CM17] Thm 3.1; [DLL] Thm 10; [CKG] App. E | Fact |
| Gibbs samplers ↔ Gibbs states | exact KMS-DB quasi-local Lindbladian for any $H$ | [CKG] Thm I.1; [DLL] Prop 20 | Fact |
| Gibbs samplers ↔ Markov property | time-averaged single-Pauli KMS Lindbladian is a quasi-local recovery map; CMI decays exponentially in distance, with prefactor $e^{\mu|A|}$ | [CR] Thm III.1, Cor III.2 | Fact |
| gap ↔ Markov property | uniform local gap ⇒ global Markov property | [CR] Def B.1, Cor B.2 | Fact |
| gap ↔ clustering (commuting) | gap of Davies/heat-bath ⇔ strong clustering; 1D; high temperature | [KB] Thms 23, 26, 30, 31 | Fact |
| gap/LSI ↔ mixing | $\sqrt{1/\sigma_{\min}}e^{-\lambda t}$, $\sqrt{2\log(1/\sigma_{\min})}e^{-\alpha_1t}$; $\alpha_1\le\lambda$; hypercontractivity ⇔ LS$_2$ | [KT] Thms 15, 16, 22 | Fact |
| LSI ↔ transport ↔ concentration | MLSI ⇒ TC$_2$ ⇒ TC$_1$ ⇒ Gaussian concentration; TC$_2$ ⇒ PI | [RD] Thms 3, 4, 6, 8 | Fact |
| Dirichlet forms ↔ W$_2$ | GNS-DB QMS = gradient flow of relative entropy for the Carlen–Maas metric; gradient flow ⇒ BKM | [CM17] Thm 7.6; [CM20] Thm 2.9; [Wir] Thms 6.26, 7.7 (tracial, infinite-dim) | Fact |
| Connes distance ↔ W$_1$ | equal on complete spin manifolds | [DM] Prop 2.1 | Fact |
| Connes distance ↔ Lip-norms | Lip-norm ⇔ Connes metric gives weak-\* topology; radius = Poincaré-type constant | [Rie] Def 5.1, Prop 2.2, §7 | Fact |
| quantum W$_1$ ↔ single-Pauli derivation | $\frac23 d_{D_P}\le W_1^{\rm DMTL}\le\frac32d_{D_P}$ | [DMTL] Prop 15 + §3.4 | Derivation (K6) |
| Dirichlet forms ↔ spectral triples | Sierpiński gasket: triple from a deformed derivation recovers dimension, measure, metric (bi-Lipschitz) and the Dirichlet form (residue); K$_1$ pairing | [CGIS] Thms 3.8, 4.3, 5.4, 5.5, Cor 5.3 | Fact |
| KMS ↔ NCG (type III) | twisted spectral triples, $\sigma$ = modular group at $-i$ | [CMo] Def 3.1, Remark 2.4, Prop 3.4 | Fact |
| KMS ↔ NCG (index) | modular spectral triples from circle actions with KMS weights; spectral flow = twisted residue cocycle; Araki relative entropy | [CNNR] Thms 1.1–1.4, §6 | Fact |
| expanders ↔ derivations | (T) ⇔ closable derivations inner (vN Delorme–Guichardet); Kazhdan-pair form | [Pet] Thms 0.1, 3.2 | Fact |
| expanders ↔ Dirichlet forms | subexponential spectral growth ⇒ amenable; (H) ⇔ discrete spectrum (Caspers–Skalski); cnd functions ↔ Dirichlet forms on $L(\Gamma)$ | [CS17] Thm 3.15, §1, Ex. 3.11 | Fact |
| expanders ↔ quantum channels | zig-zag; explicit constant-degree; random unitaries; Cayley-graph irreps; LS$_2$ of expanders $O(\log D\log\log d/\log d)$ | [BST] Thms 1–3; [Har] (3); [KT] (144) | Fact |
| graph potential theory ↔ NC operator theory | Cayley-graph condenser capacity = quasicentral modulus of $\lambda(\gamma)$ | [Voi] §14 | Fact (stated without proof in a survey) |
| standard form ↔ Gibbs samplers ↔ frustration-free Hamiltonians | the $\rho^{1/4}\cdot\rho^{1/4}$ embedding | [CS17] §2; [KB] (79), Table I; [CKG] Def II.1 | Fact (identification Mine) |
| programme | re-routing Lindbladian: kernel $D_0$, Cesàro = Takesaki expectation, recovery ⇔ sufficiency, $2\tanh$ bound | §6.1 | Derivation (K7) |
| programme | sufficiency ⇔ modular derivation inner by face observable; distance = ½ max osc | §6.4 | Derivation |

---

## 8. Messages to the other prongs (programme rule C1)

**To prong 1 (theory).**
1. **Record §6.1 as a Derivation in Note 1 §6.** A Lindbladian on the history corner whose jumps are re-routings over a fixed face is GNS-detailed-balanced for $\mathbb P_\beta$. Its Cesàro limit is the Takesaki expectation onto $D_0$, and it recovers $\mathbb P_{\beta'}$ exactly iff the barycentre is sufficient, with defect $\le2\tanh(|\beta-\beta'|\,\mathrm{osc}F/4)$. This gives the Jenčová–Petz criterion a dynamical form, and a Dirichlet form $\sum c\|[e_{\mu\nu},\cdot]\|^2$ whose kernel is exactly the face algebra.
2. **Define the flip sub-family (§6.5)** and ask for the gap of its weighted fibre Laplacian. That gap is the programme's "uniform local gap", and the place where U3's coboundary expansion would enter.
3. **Note 1 §8's Dirac operator is a choice of derivation (§6.2).** Its natural source is the arrow-jump Lindbladian, and its Connes metric is finite only on states with equal origin weights (L5). Keep the "time" derivation $\delta_F$, whose innerness by face observables is the sufficiency question (§6.4), separate from the "transport" derivations, which give metrics.
4. **Integer cocycles.** If $F$ is integer-valued (or rescaled to be), [CNNR]'s modular index theory applies to Note 1's data. That is the concrete version of U5's gap labelling to try first.

**To prong 2 (bridge).**
- The re-routing chain of §6.1 is computable from dictionary v1's fitted kernel. Report two quantities:
  - $\max_\tau\mathrm{osc}_\tau F$ together with the tanh bound;
  - the gap of the flip Laplacian on fibres. A large oscillation with a large flip gap would contradict [expanders.md](expanders.md) §11.1's Poincaré bound, and would be charged to the estimator first (P2.4).
- For a Hamming-Lipschitz statistic of faces, compare its empirical spread with Marton's $\sqrt{n\,S/2}$ (§6.6).
- As everywhere, a negative result indicts the dictionary first (P2.1).

**To prong 3 (unlocks).** Add to the §4 table of [local-to-global-unlocks.md](../../local-to-global-unlocks.md) a fifth column, "derivation language":
- **hereditary class** ↔ an inclusion $N\subset M$ and the fibre algebra $N'\cap M$;
- **restriction–co-restriction operator** ↔ the Cesàro mean of the fibre Lindbladian;
- **invariant controlled by local data** ↔ coercivity of the fibre derivation (gap or Kazhdan constant);
- **homogeneity** ↔ modular invariance of $N$, which makes the expectation exact.

Two cautions:
- The quantum Gibbs/Markov reduction of [arxiv-2609.38007.md](arxiv-2609.38007.md) §13 is the case where modular invariance fails. Its cost $e^{|A|}$ is the price of that failure.
- Bounded-degree quantum expanders are not the relevant expander notion for Gibbs samplers (§5.3).

---

## 9. What the user's intuition gets, in one paragraph (Mine)

The essence common to Gibbs/Markov structure, expanders and NCG is a **derivation into a bimodule together with its kernel**. The kernel is the conditioned-upon algebra, the Markov expectation. The derivation's $L^2$ energy is the Dirichlet form and Lindbladian. Its $L^\infty$ size is the Lipschitz seminorm and Connes metric. Its coercivity off the kernel is the gap, and uniformly over all representations it is property (T), the source of expanders. The time-averaged detailed-balanced single-Pauli Lindbladian is the **Cesàro mean of the heat flow of the vertical derivation of the inclusion $M_{A^c}\subset M_\Lambda$**, made KMS-symmetric by the quarter-twist. It is a dynamical approximation to a conditional expectation that does not exist, because $M_{A^c}$ is not modular-invariant. In the programme the corresponding inclusion $D_0\subset$ (history corner) *is* modular-invariant. So the same construction is exact there and becomes a dynamical test of the barycentre's sufficiency (§6.1). Expansion enters only as the uniformity that would upgrade local to global (§4.7): for [CR] it is the open global Markov problem, and for the programme it is the flip-graph gap (§6.5).

---

## 10. Numerical checks ([check_nc_dirichlet_lindblad.py](check_nc_dirichlet_lindblad.py), pure numpy; run output)

| id | statement | result |
|---|---|---|
| K1 | [DLL] Thm 10 sampler with single-Pauli couplings on $A=\{0\}$, $n=3$, $\beta=1.3$, random noncommuting $H$: KMS-symmetric, fixes $\rho_\beta$, $G$ Hermitian, **not** GNS | asymmetry $\le2.4\cdot10^{-15}$, $\|\mathcal L^*\rho\|\le4\cdot10^{-16}$; $[\mathcal L,\Delta]$ of relative size $3\cdot10^3$–$2\cdot10^4$ (3 trials) |
| K2 | Dirichlet form = frequency-domain commutator sum ([CR] X.1) = $\int g\|[\tilde A(t),X]\|^2_\rho dt$ | relative errors $\le4.4\cdot10^{-14}$ (frequency) and $\le9.2\cdot10^{-14}$ (time quadrature) |
| K3 | kernel of the single-Pauli-on-$A$ generator | $\beta=0$: 16 $=4^{n-1}$; $\beta=1.3$ noncommuting: 1 (3 trials); commuting Ising: 8 (classical boundary ⊗ rest) |
| K4 | gap-free Cesàro decay $\mathcal E(\mathcal R_tX)\le0.4073\|X\|^2_{\rm KMS}/t$ | holds at $t=1,10,100,1000$ |
| K5 | $\beta=0$ hypercube: $\mathcal L(S)=-4w_A(S)S$; Efron–Stein | exact on all 64 strings; max Efron–Stein ratio 0.713 |
| K6 | $\tfrac12M\le\max_i\|H-E_iH\|\le\tfrac34M$ | ratio in $[0.521,0.683]$ over 200 random $H$ |
| K7 | re-routing Lindbladian (§6.1), 3 endpoint fibres of sizes 3, 4, 2 | KMS/GNS $\le4\cdot10^{-16}$; dim ker = 3; Cesàro = Takesaki to $10^{-15}$; Dirichlet = commutator sum to $2\cdot10^{-16}$; random $F$: $\|\mathcal L^*\mathbb P_{\beta'}\|_1=0.58$ and recovery defect $0.080\le0.487$; $F=U\circ r$: both $\sim10^{-16}$ |
| K8 | $\mathcal E(a)\le\frac12\sum\sup_t\|[\tilde A(t),a]\|^2$; Poincaré; radius bound | holds; bounds $1.1\cdot10^3$–$5.3\cdot10^3$ (crude, as expected) |

These are checks of finite-dimensional identities and inequalities, not experiments about networks.

---

## 11. Sources read for this digest

- E. A. Carlen, J. Maas, arXiv:1609.01254: §1, §2 (Defs 2.1–2.2, 2.10, Lemma 2.8, Thm 2.9), §3 (Thm 3.1, Remarks 3.2–3.4), §5 (derivations, Lemma 5.9, Thm 5.10), §6 (Bose OU example), §7 (Def 7.1, Thms 7.5–7.6), App. A (end of proof of Thm 3.1), App. B (KMS ⇏ GNS).
- E. A. Carlen, J. Maas, arXiv:1811.04572: §1, §2 (Lemma 2.1, Def 2.2–2.3, Thm 2.4, Prop 2.5, (2.13), Prop 2.7, Thm 2.9), §3 (abstract Beurling–Deny), Props 4.11–4.12, §5.1–5.2.
- C.-F. Chen, M. J. Kastoryano, A. Gilyén, arXiv:2311.09207: abstract, §I (Thms I.1–I.2), §II.A–C (Def II.1, Props II.1–II.2, Cor II.2), App. D, App. E (s-detailed balance).
- Z. Ding, B. Li, L. Lin, arXiv:2404.05998: §1, §2 (Lemmas 1–3, 5, 6, Remark 7, Lemma 8, Lemma 9, Thm 10, Cor 11), §3.1, §3.4 (Prop 20), Remark 21, §4.
- M. J. Kastoryano, F. G. S. L. Brandão, arXiv:1409.3435: §I, §IV (Davies, heat-bath, Lemmas 11–13), §V (Defs 14–15), §VI (Prop 20, Def 21, Table I, (79)), §VIII (Thms 30–31), §IX.
- C.-F. Chen, C. Rouzé, arXiv:2504.02208: §I, §II (Thms II.1–II.2, Lemma II.1), §V, §VII (Lemmas VII.1–VII.3, Cor VII.1), §VIII (Lemma VIII.1, Cor VIII.1), §X (Lemmas X.1–X.4), §XI.A, App. B.
- M. Vernooij, M. Wirth, arXiv:2303.15949: §1, §2 (Thms 2.4–2.5, V-transform), §3 (Prop 3.12, Def 3.13), §4 (Prop 4.1, Thm 4.2), §5 (Thms 5.2, 5.4, Remarks 5.3, 5.5, 5.6).
- M. Vernooij, arXiv:2203.12307: §1–§5 (Thm 3.1, Thm 4.1, Examples 5.2–5.4).
- M. Wirth, arXiv:1808.05419: Introduction, §1 overview, §4 (Defs 4.1–4.12, Examples 4.16–4.19, Prop 4.20), §6.1–6.2 excerpts, §7 summary.
- F. Cipriani, J.-L. Sauvageot, arXiv:1611.01749: §1, §2 (Defs 2.1–2.2, symmetric embedding, (2.4), item (e)), §3 (Def 3.9, Lemma 3.13, Thm 3.15, Remark 3.16, Cor 3.17, Ex. 3.11), §4 (Defs 4.3–4.4, Thm 4.6).
- F. Cipriani, D. Guido, T. Isola, J.-L. Sauvageot, arXiv:1112.6401: §1, §3.3 (Thm 3.8, Prop 3.9), §4 (energy form, Thm 4.3), §5 (Cor 5.3, Thms 5.4–5.5, summary).
- F. D'Andrea, P. Martinetti, arXiv:0906.1267: §1, §2 (Prop 2.1, §2.3–2.4), §4 (noncommutative examples), §5.
- M. A. Rieffel, arXiv:math/9906151: §0–§5 (Prop 2.2, Thms 4.1–4.2, Def 5.1, Thm 5.2), §7, §11.
- N. Datta, C. Rouzé, arXiv:1704.02400: §1, §2.4–2.5 (Lemmas 2–6, Prop 1, (2.32)), §3 (Thms 3–4), §4 (Thms 8–9), §6.
- G. De Palma, M. Marvian, D. Trevisan, S. Lloyd, arXiv:2009.04469: §I–§V (Defs 5–8, Props 2, 6–9), §VII–VIII (Thms 2–3), App. A–B (Props 14–15).
- M. J. Kastoryano, K. Temme, arXiv:1207.3261: §I, §II.A (Defs 6–7, 10), §III (Def 11, Thms 15–16), §IV (Lemma 21, Thm 22), §V ((137)–(144), Cor 27, Lemma 28, Prop 31).
- A. Ben-Aroya, O. Schwartz, A. Ta-Shma, arXiv:0709.0911: full text (Defs 1.1–1.2, Props 3.1–3.2, Thms 1–3, §4).
- A. W. Harrow, arXiv:0709.1142: full text.
- J. Peterson, Pacific J. Math. 243 (2009) 181–199: §0 (Thm 0.1, Cor 0.2), §1.2, §3 (Def 3.1, Thm 3.2 with proof excerpts), via the msp.org PDF text through Exa.
- A. Connes, H. Moscovici, arXiv:math/0609703: full text.
- A. L. Carey, S. Neshveyev, R. Nest, A. Rennie, arXiv:0808.3029: §1 (Thms 1.1–1.4), §2.1, §4.2, §5.3, §6 (examples).
- D.-V. Voiculescu, arXiv:2503.23248: §10–§15 (the condenser dictionary and the Cayley-graph identity).
- Metadata only (MaRDI, Politecnico repository): F. Cipriani, *Noncommutative potential theory: a survey*, J. Geom. Phys. 105 (2016) 25–59.
- Via the sibling digests in this directory (not re-read unless stated in §0): Yang arXiv:2609.38007; Rosa-Ruiz et al.; Liu–Mohanty–Raghavendra–Rajaraman–Wu; Hastings arXiv:0706.0556; Willett–Yu; Lubotzky; Alev–Lau; Anari–Liu–Oveis Gharan.
- *From memory, not opened:* Takesaki's theorem on conditional expectations and modular invariance; the Delorme–Guichardet theorem for groups (the von Neumann-algebra version [Pet] was read); Popa's rigidity of inclusions (named in [Pet]).
