# Team F: the noncommutative local-to-global framework — REPORT (v0)

*1 Oct 2026, v0 (≈ 21:30 UTC). Scope: Garland / Oppenheim / Alev–Lau, spectral and entropic independence, Chen–Eldan localization, commuting squares and quasi-factorization; a precise noncommutative trickle-down statement; the dictionary to the competition object; from ≈ 22:30 the synthesis across teams A–E (§7, to come in the final). Labels: **Theorem** (proved here, elementary, and checked numerically), **Fact** (cited), **Measured** (run here, numbers in `results/`), **Synthesis**, **Conjecture**, **Speculation**.*

*Builds on, does not repeat: `notes/local-to-global-unlocks.md` §4–5, `notes/digests/bridges/hdx-spectral-independence.md` (§2 trickle-down, §4.5 BCR, §7 "no k-step quantum trickle-down found"), `notes/streams/bridges-synthesis/A-unification.md` (D-2 two-projection inequality, D-3 "SI is a frame bound", C-2 the missing quantum k-step trickle-down), `C-competition.md` (R1 doubled leg W∘W = (2/n)11ᵀ + E, R11 no transport gap, "no direct lever").*

## 0. Summary

1. **The spectral trickle-down has no noncommutative obstruction once it is written in the forget lattice** (Theorem F1, §2.1). For any lattice of KMS-symmetric conditional expectations E_A ("forget A"), closed under intersection of ranges, the heat-bath dynamics (1/n)Σ_i E_i has gap ≥ ∏_{m=2}^{n} (1 − (1+η_m)/m), where η_m + 1 is the frame bound ‖Σ_{i∈A}(E_{A∖i} − E_A)‖ at level |A| = m, and η_m ≤ λ_max of the matrix of pairwise commuting-square defects ‖E_{A∖i}E_{A∖j} − E_A‖. Pinning is never used: the operator norm on the whole L² is the supremum over pinnings classically, and N_A-bimodularity replaces disintegration noncommutatively. The classical case is exactly Anari–Liu–Oveis Gharan Thm 1.3 / Alev–Lau. Checked on a non-commuting 4- and 5-qubit Gibbs ring: the bound holds, is exact at β = 0 and within 1.3–1.7× at β = 0.3.
2. **The noncommutative failure is entropic, and it is one step** (§2.2). Classical local-to-global for entropy (Chen–Liu–Vigoda Thm 5.6) converts spectral to entropic contraction *inside a link*, where the densities are single-site marginals with bounded ratio b. With no links, the noncommutative comparison must be made on a whole region, at a cost that grows like e^{c|A|} (the χ_KMS growth already measured in A-unification C10–C11). **Repair hypothesis:** the Pimsner–Popa constant of the one-site inclusion N_A ⊂ N_{A∖i} (a Jones index, O(d²), not exponential) plays the role of b. Conjecture F2 states the resulting iterable approximate tensorization; its two known corners (m = 2: Bardet–Capel–Rouzé; exact commuting lattice: Lemma F2a, proved here) are both consistent with it.
3. **Dictionary, precise (§3).** The fresh-weight average E_W over rotation-invariant Gaussian weights acts on tensors through the **Brauer algebra**: a sum over pairings. The *through-string* pairings are the commuting-square part: they give the Frobenius currency and the incoherent sum over layers. The *caps* are Jones projections with loop value n and Markov trace τ = 1/n: they are partial traces, i.e. the coincident-pattern channel, the column sums of N7, the trace channel and the dilation charge. He initialisation puts the cap sector exactly at eigenvalue 1 (loop n × variance 2/n × gate density ½ = 1), a unipotent sector with no mixing (angle 0). By Proposition F3 an angle-0 sector cannot be quasi-factorized; it must be **pinned** (localized), which is Chen–Eldan / Bauerschmidt–Bodineau's move.
4. **Measured at n = 1024 (§5.2, three networks).** The cap (trace) channel holds 0.66–0.79 of the energy of the old (2,1) slice at layers ≥ 8. With costate's scale mode the union holds 0.82 at layer 15 (random-direction control ≤ 0.005). The remaining through-string bulk (18 % of old energy at A = 1, still 12.5 % beyond age 3 and 3.5 % beyond age 7) forgets only slightly faster than the cap sector. Priced by the ε² law, that bulk costs ≈ 3e-7 at A = 3. So the cap sector is most of the memory, not all of it, and the bulk is the decisive remainder.
5. **Best transfer (§4, T1 + T2).** Carry the cap sector at law level by localization: the scale field h_l plus the collective coordinate g pinned at K quadrature nodes, O(K n²) per layer on top of one Gaussian closure. Treat the through-string bulk as *iterable*: since weight layers form an exact commuting lattice, per-layer compression errors of the bulk add incoherently and do not compound. Kill criteria and costs are in §4.

## 1. Instances, to the level of their proofs

### 1.1 Garland → Oppenheim → Alev–Lau, and its spin-system form (ALO, Chen–Eldan Thm 24)

- **Fact** (Oppenheim, arXiv:1709.04431; Alev–Lau, arXiv:2001.02827 Thm 1.5; ALO, arXiv:2001.00303 Thm 1.3). For a pure weighted complex whose links have second eigenvalues γ_j, λ₂(P^▽_k) ≤ 1 − (1/(k+1)) ∏_{j=−1}^{k−2}(1 − γ_j). For a spin system with η_i-spectral independence under all pinnings of i sites, the Glauber gap is ≥ (1/n) ∏_{i=0}^{n−2}(1 − η_i/(n−i−1)).
- **Proof skeleton** (Garland's method). Two ingredients:
  - (i) *localization*: the walk at level k is the average over (k−1)-faces τ of the walk in the link of τ, ⟨f, Mf⟩ = E_τ⟨f_τ, M_τ f_τ⟩;
  - (ii) *the mean of the localized function is the global walk*: the constant part of f_τ in the link of τ is (Mf)(τ).

  Then ⟨Mf, f⟩ ≤ λ‖f‖² + (1−λ)‖Mf‖², and the spectral theorem gives μ ≤ λ/(1−λ).
- **Where the classical proof uses commutativity.** Only in (i): a link is a *pinning*, a disintegration of L²(μ) over the values of the pinned coordinates. That is a disintegration over the centre of the conditioning algebra L^∞(x_τ), which is the whole algebra because it is abelian. Noncommutatively the conditioning algebra is a factor and has no such disintegration. This is the step the "quantum trickle-down" literature has been unable to reproduce (hdx digest §4.5, §7).
- **What §2.1 shows.** Step (i) is unnecessary. Written in the *forget lattice* (algebras of everything except A, not marginals on A), the classical theorem is a statement about a lattice of orthogonal projections. The supremum over pinnings is the operator norm.

### 1.2 Quasi-factorization: the two-algebra case and why it does not iterate

- **Fact** (Bardet–Capel–Rouzé, arXiv:2001.07981, as recorded in the hdx digest §4.5).
  - If E₁E₂ = E₂E₁ = E_M (a commuting square), then D(ρ‖E_M*ρ) ≤ D(ρ‖E₁*ρ) + D(ρ‖E₂*ρ).
  - Otherwise the inequality holds with constant 1/(1 − 2c₁) and an additive noncommutativity correction d, where c₁ = ‖E₁E₂ − E_M : L¹(σ) → L^∞‖ and d = 0 for classical Hamiltonians.
- **Fact** (A-unification D-2). The L² version is the sharp two-projection inequality (1 − c)‖g − Rg‖² ≤ ‖g − Pg‖² + ‖g − Qg‖², with c the cosine of the Friedrichs angle.
- **Why it does not iterate.** Applied to a chain of n regions by repeated bisection, the constants multiply: ∏(1 − 2c)^{−1} over log n levels, with c bounded away from 0 for adjacent regions. The classical theory iterates *level by level in the number of free sites* (§1.1), not by bisection. The entropic version of that level-by-level step is what is missing (§2.2).
- **Rewriting BCR in level form.** By the chain rule for nested conditional expectations, the two-algebra inequality is equivalent to Σ_{i=1,2} D(E_i*ρ‖E_M*ρ) ≤ (1 + 2c₁) D(ρ‖E_M*ρ). That is the m = 2 case of the entropic frame bound in Conjecture F2, with η^{ent}_2 = 2c₁.

### 1.3 Localization schemes (Chen–Eldan, arXiv:2203.04163)

- **Fact** (hdx digest §2.8). A localization process is a measure-valued martingale ν_t, with ν_0 = μ and ν_∞ a point mass. Its associated Markov chain resamples from ν_τ. The gap is bounded below by the *approximate conservation of variance*, inf_f E Var_{ν_τ}(f) / Var_μ(f). For stochastic localization this ratio is controlled by ∫‖Cov(ν_t)‖_op dt. Coordinate-by-coordinate pinning reproduces ALO.
- **The estimation form, which is the one relevant to us.** E_μ[f] = E[E_{ν_τ} f] exactly. If a closure G is accurate on the localized laws ν_τ (because they are closer to Gaussian), then E_μ f ≈ Σ_k w_k G(ν_{τ_k}). The error is E[G(ν_τ) − E_{ν_τ} f] plus quadrature error. The variance-conservation ratio becomes a *non-Gaussianity-removal ratio*: how much of κ3, κ4 the pinning removes.
- **Fact** (Chen–Eldan Prop 19, 21, 27). gap(P) = inf_φ E[Var_{ν_τ}φ] / Var_ν φ. Under (κ_i)-variance conservation, E[Var_{ν_i}φ | ν_{i−1}] ≥ (1 − κ_i) Var_{ν_{i−1}}φ, one has gap ≥ ∏(1 − κ_i), and likewise for entropy and MLSI.
  - Annealing (Thms 46/47): if a first localization preserves an ε fraction of variance, and the localized measures have gap ≥ δ, the gap is ≥ εδ.
  - Thm 49 (Ising): ρ_LS(Glauber on ν) ≥ ε · exp(−2‖J‖_op ∫₀¹ α(λ) dλ), where α bounds ‖Cov‖ along the interpolation.
- **"Pin the bad directions" (Facts).**
  - Bauerschmidt–Bodineau, arXiv:1712.03676, Thm 1. One Gaussian-convolution ("renormalisation") step splits M⁻¹ = c⁻¹ id + B⁻¹. Conditionally on the field the measure is a product, and the field law is uniformly convex, giving an LSI with constant (2/γ)(1 + 2n‖M‖/(n − ‖M‖)).
  - Eldan–Koehler–Zeitouni, arXiv:2007.08200, Thm 1: (1 − ‖J‖_op) Var φ ≤ E(φ, φ) for 0 ⪯ J ≺ Id. Stochastic localization reduces J to rank-one models.
  - Most relevant to Proposition F3: Boban–Li–Oveis Gharan, arXiv:2609.13138, a **rank-1-perturbed trickle-down** that subtracts the rank-one part from the cross terms (SK, β < 1/2 + 5·10⁻⁵).
  - Mikulincer–Sohn, arXiv:2512.22803: Hubbard–Stratonovich + localization with the rank-one part handled by coordinate localization.

  These are the classical instances of "factorize the commuting part, pin the defect sector". (Constants of the last two are from alphaXiv summaries, not the PDFs.)
- **Quantum analogue (proposal 2 of the input).** The probe found no paper using quantum filtering or continuous weak measurement as a localization scheme for gap or MLSI bounds. With Theorem F1 in hand, the noncommutative statement to aim at is the entropic Chen–Eldan Prop 27 for a measurement-induced martingale of states ρ_t. Its "κ_i" would be the entropic frame constants of Conjecture F2 restricted to the measured observable's algebra. That is a Conjecture, not pursued in v0.

### 1.4 Jones' tower and the Markov trace (the fresh-weight lemma in its native form)

- **Fact** (Jones). In the tower N ⊂ M ⊂ M₁ = ⟨M, e⟩, the Jones projection e satisfies e x e = E_N(x) e for x ∈ M, and the Markov trace satisfies tr(x e) = τ tr(x), with τ = [M:N]^{−1}. So the new projection is independent of the past in the trace.
- In §3.1 this is the fresh-weight lemma literally: the caps of the Brauer algebra of O(n) are Jones projections with loop value δ = n. Their normalised trace property is "the next layer's weights are independent of everything computed before".

## 2. The noncommutative generalisation

### 2.1 Theorem F1 (k-level spectral trickle-down in a forget lattice)

**Setting.**
- H is a Hilbert space and V = {1, …, n} a set of pieces.
- For A ⊆ V, R_A ⊆ H is a closed subspace with R_∅ = H and R_A = ∩_{i∈A} R_{{i}}. Forgetting more gives a smaller range: B ⊆ A ⇒ R_A ⊆ R_B.
- E_A is the orthogonal projection onto R_A.
- **Quantum instance.** H = L²(M, σ) with the KMS inner product; L_i are σ-KMS-detailed-balanced local Lindbladians (Chen–Kastoryano–Gilyén, arXiv:2311.09207); R_i = ker L_i. Then R_A = ker Σ_{i∈A} L_i, and E_A = lim e^{t L_A} is a σ-preserving, KMS-symmetric conditional expectation onto the fixed-point algebra (Frigerio).
- **Classical instance.** R_A = L²(x_{V∖A}), so E_A is the heat bath on A.

**Definitions.**
- For |A| = m ≥ 2 and i ∈ A, let P^A_i := E_{A∖i} − E_A. This is the orthogonal projection onto H^A_i := R_{A∖i} ⊖ R_A, the "one more piece free" subspace.
- The frame constant is 1 + η_m := max_{|A| = m} ‖Σ_{i∈A} P^A_i‖.
- The pairwise defect matrix is C^A_{ij} := ‖E_{A∖i}E_{A∖j} − E_A‖ for i ≠ j (zero diagonal). This is the approximate-commuting-square defect of the square (R_{A∖i}, R_{A∖j}; R_A).

**Theorem F1.**
- (a) η_m ≤ max_{|A|=m} λ_max(C^A).
- (b) For every f ∈ H and every A with |A| = m ≥ 2: avg_{i∈A} ‖f − E_{A∖i} f‖² ≥ (1 − (1+η_m)/m) ‖f − E_A f‖².
- (c) The heat-bath operator P = (1/n)Σ_i E_i satisfies: gap(P) on R_V^⊥ ≥ ∏_{m=2}^{n} (1 − (1+η_m)/m) = (1/n) ∏_{m=2}^{n}(1 − η_m/(m−1)).
- (d) In the quantum instance, if each −L_i ≥ g (1 − E_i) in the KMS order, then gap(Σ_i L_i) ≥ n g ∏_{m=2}^{n}(1 − (1+η_m)/m).

*Proof.*
- (a) P^A_i P^A_j = E_{A∖i}E_{A∖j} − E_A, because E_A E_B = E_B E_A = E_A whenever B ⊆ A. Then ‖Σ P_i‖ = ‖T*T‖ for the synthesis map T: ⊕H_i → H, and T*T is the block matrix with identities on the diagonal and P_iP_j off it. Its norm is at most 1 + λ_max of the matrix of block norms.
- (b) Write V_A(f) := ‖f − E_A f‖². Since R_A ⊆ R_{A∖i}, Pythagoras gives V_A(f) − V_{A∖i}(f) = ‖P^A_i f‖². Summing over i ∈ A: Σ_i (V_A − V_{A∖i}) = ⟨f, Σ_i P^A_i f⟩ ≤ (1 + η_m) ‖(1 − E_A) f‖² = (1 + η_m) V_A, because every H^A_i ⊥ R_A. Rearranging gives (b).
- (c) Let a_m := avg_{|A|=m} V_A(f). A uniform m-set with one uniform element removed is a uniform (m−1)-set, so a_{m−1} ≥ (1 − (1+η_m)/m) a_m. Iterate from a_n = ‖f − E_V f‖² down to a_1 = ⟨f, (1 − P) f⟩.
- (d) Standard comparison. ∎

**Remarks.**
- *Classical case.* Fibrewise over pinnings of x_{V∖A}, H^A_i is the space of centred functions of x_i. So the frame bound is "Var(Σ f_i) ≤ (1 + η) Σ Var f_i for every pinning", which is spectral independence (A-unification D-3; ALO Thm 1.11). The operator norm on L²(μ) equals the supremum of the fibre norms, because Σ P^A_i commutes with multiplication by L^∞(x_{V∖A}).
- *What replaces pinning noncommutatively.* Each E_{A∖i} is an N_A-bimodule map (N_A := R_A as an algebra). So Σ_i P^A_i commutes with left and right multiplication by N_A, and the frame bound holds in every inner product reweighted by N_A, e.g. by the relative density of E_A*ρ. This is the noncommutative substitute for "the inequality holds in every link". It is what an entropic version must exploit (§2.3).
- *What does not transfer for free.* Positivity. The heat-bath operator is a quantum channel only when the E_A are genuine conditional expectations. For non-commuting Gibbs states the KMS-orthogonal projection onto "operators trivial on A" is not one (Takesaki; A-unification §4.1). The theorem itself does not care, which is why §5.1 can test it on those subspaces.
- *Novelty.* The classical content is ALO / Alev–Lau. The literature probe (alphaXiv + web, abstracts and main theorems, Oct 2026) found **no k-level quantum spectral-independence or trickle-down theorem**. The closest are:
  - two-block / martingale gap theorems: arXiv:2505.08991, gap(D) ≥ min_x gap(D_x) · gap(**H**) via ‖P_A P_B − P_{A∪B}‖;
  - Hayakawa–Southwell–Leditto–Chen–Hsieh, arXiv:2609.39802, "double localization" gap assembly;
  - LaRacuente's multiplicative quasi-factorization for J conditional expectations, arXiv:1912.00983: constant → 1 as ‖E₁⋯E_J − E‖_⋄ → 0, O(ln C) in the Pimsner–Popa index C.

  None is level-by-level with a product formula. Absence from a search is not proof, and the statement is elementary, so the likely status is "folklore but unwritten". Its value here is that it locates the noncommutative difficulty exactly (§2.2).

**Measured** (`nc_trickle.py`, `results/nc_trickle.txt`). The setting is a ring of n qubits with a random Heisenberg + transverse/longitudinal-field Hamiltonian (non-commuting), the KMS inner product, and R_A = operators acting trivially on A.

| n | β | η₂, η₃, η₄, η₅ | λ_max(C) at each level | bound ∏ | true gap | ratio |
|---|---|---|---|---|---|---|
| 4 | 0 | 0, 0, 0 | 0, 0, 0 | 0.2500 | 0.2500 | 1.00 |
| 4 | 0.3 | 0.366, 0.602, 0.970 | 0.366, 0.605, 0.983 | 0.0750 | 0.0993 | 1.33 |
| 4 | 1.0 | 0.894, 1.677, 2.652 | 0.894, 1.678, 2.652 | 0.0005 | 0.0078 | 15.7 |
| 5 | 0.3 | 0.362, 0.580, 0.714, 0.630 | 0.362, 0.581, 0.740, 0.785 | 0.0581 | 0.0987 | 1.70 |

The bound holds everywhere. It is tight for products and within 2× at high temperature. Part (a) is nearly an equality, so **pairwise commuting-square defects determine the frame constants**. At β = 1 the level factors approach 0 (η₂ = 0.89), and the bound, though valid, is 16× pessimistic.

### 2.2 Where the noncommutative theory actually fails: the entropic level step

- **The chain rule survives.** For genuine conditional expectations with R_A ⊆ R_{A∖i}, D(ρ‖E_A*ρ) = D(ρ‖E_{A∖i}*ρ) + D(E_{A∖i}*ρ‖E_A*ρ) (Petz). So the entropic analogue of (b) is equivalent to the **entropic frame bound**

  EFB(η): Σ_{i∈A} D(E_{A∖i}*ρ ‖ E_A*ρ) ≤ (1 + η^{ent}_m) D(ρ ‖ E_A*ρ),

  and EFB at every level gives D(ρ‖E_V*ρ) ≤ [∏_m (1 − (1+η^{ent}_m)/m)]^{−1} (1/n) Σ_i D(ρ‖E_i*ρ). That is iterable approximate tensorization, hence MLSI for the heat bath and, by comparison, for Σ L_i.
- **The exact failing step.** Classically, EFB is derived from the spectral frame bound inside a link (Chen–Liu–Vigoda, arXiv:2011.02075, Thm 5.6): after pinning, the objects E_{A∖i}*ρ are single-site marginals, their density ratios are bounded by the marginal bound b, and a reverse Pinsker inequality at ratio b converts χ² back into relative entropy. Noncommutatively there is no link. The comparison must be made between ρ and E_A*ρ on the whole region, where the best ratio bound is the Pimsner–Popa constant of N_A ⊂ M, i.e. d^{−2|A|}: exponential. This is the same exponential that A-unification measured as χ_KMS growth (×2–3.5 per site, C10–C11). It is also the e^{μ|A|} in Chen–Rouzé's Thm III.1 (arXiv:2504.02208).

### 2.3 Repair hypothesis and Conjecture F2

**Lemma F2a (exact commuting lattice; proved here).** If E_S E_T = E_{S∪T} for all S, T ⊆ A (every square commutes), then Σ_{i∈A} D(E_{A∖i}*ρ‖E_A*ρ) ≤ D(ρ‖E_A*ρ), i.e. η^{ent} = 0.

*Proof.* Order A = {i₁, …, i_m} and put B_k = A ∖ {i₁, …, i_k}. The chain rule telescopes: D(ρ‖E_A*ρ) = Σ_k D(E_{B_k}*ρ ‖ E_{B_{k−1}}*ρ). Apply the channel E_{A∖i_k}* to both arguments of the k-th term. By commutativity E_{A∖i_k}E_{B_k} = E_{A∖i_k} and E_{A∖i_k}E_{B_{k−1}} = E_A, so data processing gives D(E_{B_k}*ρ‖E_{B_{k−1}}*ρ) ≥ D(E_{A∖i_k}*ρ‖E_A*ρ). Sum over k. ∎

(For tensor-product σ this is superadditivity of relative entropy with respect to a product reference. The lemma is the commuting-square form of it.)

**Repair hypothesis (H-PP).** For every A and i ∈ A, the one-site inclusion N_A ⊂ N_{A∖i} has Pimsner–Popa constant λ_i, i.e. E_A(X) ≥ λ_i X for X ∈ N_{A∖i}⁺, with λ_i ≥ λ₀ > 0 independent of |A|. For a qudit site with a commuting Hamiltonian this is a local Jones index (Chen–Rouzé digest §6.2: "4^{|A|} is a Jones index"). Dually, E_{A∖i}*ρ ≤ λ_i^{−1} E_A*ρ on N_{A∖i}: a *one-site* ratio bound, the analogue of b-marginal boundedness.

**Conjecture F2 (noncommutative entropic trickle-down).** Assume H-PP and the multiplicative pairwise defects c^A_{ij} := ‖E_{A∖i}E_{A∖j} − E_A : L¹(σ) → L^∞(σ)‖_{cb}. Then EFB holds with η^{ent}_m ≤ K(λ₀) · max_{|A|=m} λ_max(c^A), for a function K of the one-site index only. Consequently, if Σ_j c^A_{ij} ≤ η for all i and A (a noncommutative Dobrushin / spectral-independence condition), the heat-bath dynamics satisfies approximate tensorization with constant n^{O(K η)}, and the Lindbladian Σ_i L_i satisfies an MLSI with constant ≥ g / n^{1+O(Kη)}.

- *Consistency with known cases.*
  - At m = 2 it is BCR (Thm 2: c = 1/(1 − c₁) plus an additive d, which vanishes classically).
  - For an exact commuting lattice it is Lemma F2a (η = 0). For two algebras that case is Gao–Junge–LaRacuente arXiv:1710.10038 Cor 2.3, constant 1; they also show that holding for all ρ forces the commuting square.
  - LaRacuente arXiv:1912.00983 already has a J-fold multiplicative quasi-factorization with an O(ln C) Pimsner–Popa dependence. That is the strongest evidence that H-PP is the right hypothesis. What it lacks is the *level-by-level* form with a pairwise-defect matrix, which is what makes the constant polynomial rather than exponential in the number of pieces.
- *The proof route this suggests.* Run CLV Thm 5.6 with N_A-valued (bimodule) relative entropy in place of the pinned relative entropy, using:
  - the bimodule invariance of the frame bound (§2.1 remark);
  - a reverse Pinsker inequality at the one-site PP ratio, which holds for operator ratios bounded by λ₀^{−1} (operator convexity of t log t);
  - the conditional relative entropy of Capel–Lucia–Pérez-García as the N_A-valued object (Yang digest B6).

  The step to check first is whether the reverse Pinsker inequality can be applied to E_{A∖i}*ρ against E_A*ρ *uniformly in the N_A-component*. That is precisely what pinning gives for free classically.
- *Status.* Conjecture. The noncommutative corrections d of BCR (non-zero for non-commuting Hamiltonians) may force an additive term, which would make the statement a "quasi" tensorization only.

### 2.4 Proposition F3 (a central defect must be pinned, not factorized)

**Proposition F3.** Let P, Q be orthogonal projections with R = P ∧ Q, and let Z be a projection commuting with P and Q such that PQ − R = Z(PQ − R)Z. Then:
- on Ran(1 − Z), P and Q commute and ‖(1 − Z)(f − Rf)‖² ≤ ‖(1 − Z)(f − Pf)‖² + ‖(1 − Z)(f − Qf)‖² exactly;
- on Ran Z the constant is 1/(1 − c_Z), with c_Z = ‖Z(PQ − R)‖.

If c_Z = 1 (an angle-0 sector, i.e. a conserved quantity), no factorization constant exists on Z. The only way to recover a global statement is to **condition on the Z-sector**, pinning its value and averaging over its law. Then each conditional problem lives on Ran(1 − Z), where the square commutes.

*Proof.* Z commutes with P, Q and R, so everything splits into the two blocks. The first block is the commuting-square case; the second is A-unification D-2. ∎

This is the abstract reason that localization (Chen–Eldan, Bauerschmidt–Bodineau) and not factorization is the right tool for conserved or critical sectors. It is the shape used in T1.

## 3. The dictionary, precise (seed 4)

| role | object in the competition | statement | status |
|---|---|---|---|
| pieces | **pairings** of the fresh Gaussian weight legs (Brauer diagrams), not neurons | E_W[W^{⊗2k}] acts on tensors in index space as a sum over perfect matchings of the 2k legs (Isserlis / hafnian) | Fact |
| commuting-square part | **through-strings**: pairings that connect copy-1 legs with copy-2 legs | off coincident output indices only these survive: E_W[(WᵀδW)²_{ab}] = (2/n)²‖δ‖²_F for a ≠ b. This is the Frobenius currency of the fresh-weight lemma (region §1) | Theorem (Isserlis, below) |
| defect of the square | **caps** (self-pairings within one copy) = partial traces | on coincident output indices a = b the extra terms are (2/n)²[(tr δ)² + ‖δ‖²_F]. For κ3: Σ_i δ_{iik}, the column sums. These are the coincident-pattern channel (theory judge §2), the column means of N7, and C-competition R1's Perron mode of W∘W | Theorem + Synthesis |
| Jones projection | a cap e on legs (j, j+1): e² = n·e, so p = e/n is a projection | Markov trace tr(x p) = n^{−2}·tr(x)-type independence = "fresh weights are independent of the past". Index [M:N] = n²: far beyond the "wall at 4" | Synthesis (the O(n) Brauer algebra is the Schur–Weyl dual of the weights' rotation invariance) |
| books close | Theorem F1 at η = 0 on the lattice of weight layers {W_l} (independent ⇒ exact commuting lattice) | errors tensorize over layers. The measured 2–6 % (46 % on MLP 1) is the inter-layer defect, i.e. cap-sector cross terms. The frame-bound form predicts MSE ∈ [1 + λ_min, 1 + λ_max] × Σ_l K(l)‖δ_l‖² with λ from the inter-layer correlation matrix of propagated local errors | Synthesis; test T3 |
| down step | E_W over the newest layer (annealed transport): onto Brauer invariants | the slice chain S_{k+1} ≈ (M∘M)ᵀ S_k M (costate C3) keeps exactly the diagonal of the Khatri–Rao square; the dropped a ≠ b part is through-string bulk of the same order | Fact (costate C3) + reading |
| up step | the quenched layer with fixed W (one Kraus operator, no gap, A-unification §9.2) followed by the ReLU's per-neuron re-birth | re-birth = the exact coincident response (theory §2 table) | Fact |
| angle per sector | caps: cos θ = 1 (eigenvalue exactly 1, Euler relu′(z)z = relu(z), costate C6); through-strings: transported with mean zero, energy retained slowly (R11: Lyapunov gap 0.013–0.018) | the cap sector is unipotent at He-criticality: loop value n × weight variance 2/n × gate density ½ = 1 | Synthesis |
| gates | diagonal projections Γ_l = diag 1(z_l > 0); a gated propagator is an alternating product of projections and Gaussian matrices | two free projections of trace ½: arcsine angle law (foundations F9.1, team C) | Fact |
| OU noise operator | P_ρ = down–up with angle ρ on the Gaussian input; Wiener chaos k contracted by ρ^k | K_l(ρ) = E[a_l(x)a_l(x′)ᵀ] interpolates between global (ρ = 0, m mᵀ) and local (ρ = 1). Its slope 1 at ρ = 1 is the same unipotent wall as the cap sector | Fact + Synthesis |
| dilation charge | the cap sector's invariant: q = ‖ã‖²/n (the O(n) Casimir) and its correlations with the state (the scale field h = Cov(q, z)) | conserved charge of every arrow (costate C6) = the angle-0 sector of Proposition F3 ⇒ must be pinned (localized), never factorized | Synthesis |

**Isserlis check of rows 2–3.** Let W have i.i.d. N(0, 2/n) entries and δ be symmetric. Then E[W_ia W_jb W_ka W_lb] has three pairings:
- (ia, jb)(ka, lb) needs a = b and gives (tr δ)²;
- (ia, ka)(jb, lb) gives ‖δ‖²_F for every (a, b);
- (ia, lb)(jb, ka) needs a = b and gives ‖δ‖²_F.

The through-string (ia, ka)(jb, lb) is the only pairing that survives off the coincident pattern. ∎

**Why this is not the per-neuron drawing.** The pieces are the *pairings of weight legs*, and the decomposition is by diagram type (through-string versus cap), not by neuron, layer or age. The same cap appears at every layer and every age, and its image (partial traces) is a handful of n-vectors per layer whatever the number of neurons or sources. The bulk is defined only as what the caps leave, and it carries no neuron basis at all (it is high-rank and incoherent, §5.2).

## 4. Transfers

### T1. Pin the cap sector (localization on the dilation / trace sector at law level)

- **Roles.**
  - Pieces: the cap-sector coordinates, i.e. the collective scale q_l(x) and the top principal coordinate g_l(x) (cos 0.98 with the mean at depth; interpolation §8).
  - Down: condition on (q, g).
  - Up: average over their law with K Gauss–Hermite nodes.
  - Angle: 0 on the pinned sector, so Proposition F3 says it must be pinned.
  - Global quantity: the final means.
- **Not the per-neuron drawing.** The pinned object is one O(n)-invariant scalar (or a scalar plus its correlation field h ∈ ℝⁿ) per layer, shared by all neurons and all ages.
- **Prediction.**
  - Pinning carries the cap share of old content at law level, to all cumulant orders at once: ≈ 0.80–0.85 of old-slice energy at depth (§5.2), 93 % of per-neuron κ3 variance (theory §2, trace channel) and ≈ 78 % of coherent kurtosis (interpolation §8).
  - It cannot carry the through-string bulk (≈ 15–20 % of old energy). That bulk, dropped, costs ≈ 3.8e-6 × 0.7 × 0.18 ≈ 5e-7 by the ε² law (region N2). So T1 alone predicts raw ≈ (3–5)e-7 at n = 1024, not the bar.
- **Cost.** O(K n²) per layer on top of one Gaussian closure (≈ 2 u per layer): h′ = Wᵀ(Φ ∘ h) is a matrix–vector product; the rank-1 conditional update of C per node is O(n²). Total ≈ 32 u + K·O(L n²) ≈ 0.03 B.
- **Cheapest decisive tests** (both already queued by other streams; F supplies the reason and the kill line):
  - (a) costate §6, scale field at law level: must beat A3gsl_nc's 3.25e-7 at ≤ 300 products;
  - (b) bethe `localized()` (pin the Perron input direction, K-node GH): must cut v4's 4.3e-7 by ≥ 2×.
- **Kill criterion.** If (a) and (b) each fail to reach ≤ 2e-7 raw at n = 1024, the cap sector is not being pinned by these constructions. Then charge the dictionary first: the pinned coordinate must be the O(n) Casimir q and h, not a linear input direction. Bethe's input-direction pinning captures only the first-chaos part of g, while q is second chaos in x.
- **Explains** established fact 4 (Section 3 of the brief: leaders carry memory at about zero cost) only partly: ~80 % of memory is one sector that costs nothing to carry.

### T2. The through-string bulk is iterable: compress it per layer, and errors will not compound

- **Roles.**
  - Pieces: the bulk old content (the old slice minus the cap sector).
  - Down: per-layer re-fit to a merged R-atom representation (CP merge).
  - Up: transport the atoms' legs by Φ, W (3 products per R = n atoms).
  - Angle: the commuting square of the weight lattice is exact (c = 0, independent layers), so by Theorem F1 / Lemma F2a at η = 0, per-layer compression errors tensorize. They add as Σ_l K(l)‖e_l‖² rather than compounding multiplicatively.
- **Prediction.**
  - Tolerance: the bulk is ≈ 0.18 of old energy, and old content is ≈ 0.7 of D21 at depth. A uniform relative merge error e on the bulk costs ≈ 3.8e-6 × 0.7 × 0.18 × e², so the full bar budget (≈ 1e-8) allows e ≈ 12 %, and half of it e ≈ 9 %.
  - This is looser than region's ≈ 7 % on the whole of D21, because the cap sector is removed first and carried exactly by T1.
- **Cost.** 3 u per layer for transport plus the ALS re-fit (≈ 2–4 u per sweep at R = n), i.e. ≈ 8–10 u per layer, inside the leaders' 7–10 u/layer envelope.
- **Cheapest decisive test.** Region §6's experiment, sharpened: inside FC, carry the cap sector of old sources by T1 and CP-merge only the bulk at R ∈ {n/2, n}, on MLPs 0–2 at n = 1024.
- **Kill criterion.** If the bulk at R = n has merge error > 15 % per layer, or raw > 6e-8, then the bulk is not compressible: a generic high-rank random tensor has CP rank ≫ n, and region found each source atom-complete. The memory question then reduces to whether the leaders spend their error budget on the bulk (region §6 (b)).
- **Risk, stated now.** Removing the structured cap sector leaves the *least* compressible part. Our measurement (§5.2: bulk high-rank, slowly decaying in age) makes the negative branch likely.

### T3. Books close as a frame bound over the weight lattice

- **Roles.**
  - Pieces: per-layer local errors δ_l.
  - Down/up: propagation by the exact closure.
  - Angle: the inter-layer correlation matrix ρ_{ll′} of the propagated per-neuron error vectors.
- **Prediction.** By Theorem F1(a)–(b) in its two-sided form, final MSE ∈ [1 + λ_min(ρ), 1 + λ_max(ρ)] × Σ_l K(l)‖δ_l‖². The 46 % over-prediction on MLP 1 must then appear as a negative eigenvalue concentrated on cap-sector (trace) errors at neighbouring layers.
- **Test.** Inject each layer's measured local error separately and propagate. That costs L Gaussian closures per network (≈ 0.5 B in prototype time, minutes), using region's MC atlas.
- **Kill criterion.** If the off-diagonal mass of ρ is not concentrated in the trace channel, the cap/through-string split does not explain the cancellation.
- **Value.** It turns region's K(l) pricing into a certified two-sided bound, and identifies which coherent errors cancel. This is low competition value, but it is the cleanest test of the dictionary.

## 5. Tests run

### 5.1 Theorem F1 on non-commuting quantum Gibbs states

Table in §2.1 (`nc_trickle.py`; runs in seconds; n ≤ 5 because the operator space has dimension 4^n).

### 5.2 How much of the old (2,1) slice is the cap sector, at n = 1024 (`trace_share.py`)

**Method.**
- Inputs: costate's exact first-order co-state on bench `w1024_d16`, with the exact old slice S_old (sources of age > A) recorded at each layer.
- Shares measured, as fractions of the off-diagonal energy ‖S_old‖²_F:
  - **trace[u]**: the projection onto the cap channel {u xᵀ + y uᵀ}, for u ∈ {s, s², 1, m, Φ};
  - **trace[s², m]**: the two-vector channel;
  - **scale mode**: costate's dense template 2 m_p C_pq + m_q C_pp;
  - **union**: the union of the scale mode and trace[s², m].
- Control: the channel for a random u.

**Measured (A = 1; all three MLPs agree within ±0.05).**

| layer | 4 | 6 | 8 | 10 | 12 | 15 |
|---|---|---|---|---|---|---|
| scale mode (MLP 0) | 0.39 | 0.56 | 0.69 | 0.73 | 0.75 | 0.79 |
| trace[s²] (MLP 0) | 0.52 | 0.64 | 0.71 | 0.73 | 0.73 | 0.76 |
| trace[Φ] (MLP 0) | 0.39 | 0.46 | 0.54 | 0.50 | 0.56 | 0.53 |
| trace[s², m] (MLP 0 / 1 / 2) | 0.52 / 0.55 / 0.55 | 0.64 / 0.68 / 0.61 | 0.72 / 0.72 / 0.66 | 0.74 / 0.74 / 0.67 | 0.75 / 0.74 / 0.71 | 0.79 / 0.76 / 0.76 |
| **union scale ∪ trace (MLP 1 / 2)** | 0.56 / 0.56 | 0.71 / 0.64 | 0.75 / 0.69 | 0.78 / 0.72 | 0.79 / 0.75 | **0.82 / 0.82** |
| random-u control | 0.002 | 0.002 | 0.001 | 0.001 | 0.005 | 0.001 |

**Age profile (MLP 0, layer 15).** These are tail energies of the old slice for sources older than A, split into the cap channel trace[s², m] and the bulk (the rest). Cross-age terms are ignored; ages are nearly orthogonal, foundations F8.3.

| A | E_old | E_cap | E_bulk | cap share | union share | bulk as % of the A = 1 old energy |
|---|---|---|---|---|---|---|
| 1 | 5.76 | 4.53 | 1.23 | 0.79 | ≈ 0.82 (MLPs 1, 2) | 21 % |
| 2 | 4.74 | 3.78 | 0.97 | 0.80 | — | 17 % |
| 3 | 3.90 | 3.17 | 0.72 | 0.81 | 0.87 | 12.5 % |
| 5 | 2.46 | 2.05 | 0.41 | 0.83 | 0.88 | 7.1 % |
| 7 | 1.40 | 1.20 | 0.20 | 0.86 | 0.90 | 3.5 % |

**Readings.**
- **The cap channel is real and dominant.** It holds 0.66–0.79 of old energy at layers ≥ 8 on all three networks, against a 0.1–0.5 % control. With costate's dense scale mode the union holds 0.82 at layer 15 (A = 1) and up to 0.90 (A = 7).
- **The cap sector is the scale direction, not the gate direction.** The choice of u among s, s², 1 hardly matters (cos ≥ 0.95 between them), while the u = Φ channel is clearly weaker (0.50–0.56).
- **The bulk forgets, but slowly.** Per two ages its tail energy falls by ≈ 0.57–0.49, against ≈ 0.65–0.59 for the cap tail. So the through-string sector has a small but non-zero angle, and the cap sector a smaller one. This is consistent with C-competition R11 (no transport gap) and with costate's "A = 7 still leaves 20 %".
- **Price of the uncarried bulk** (ε² law, region N2, with old content ≈ 0.7 of D21 at depth). Carrying the cap sector exactly and a window of A ages leaves an uncarried bulk of ≈ 3.8e-6 × 0.7 × f:
  - ≈ 5.6e-7 at A = 1;
  - ≈ 3.3e-7 at A = 3;
  - ≈ 1.9e-7 at A = 5;
  - ≈ 0.9e-7 at A = 7.

  That is Frobenius pricing. Costate's oracles show the readout reads old content through fewer modes, so these are upper estimates. Even so, **pinning the cap sector plus a window does not reach the bar on its own**. The bulk has to be carried by something (T2) or the readout metric has to be shown to suppress it.


## 6. Honest assessment

- **Theorem**:
  - F1 (spectral k-level trickle-down in any forget lattice; elementary, checked numerically);
  - Lemma F2a (entropic, exact commuting lattice);
  - Proposition F3 (central defect);
  - the Isserlis split of the fresh-weight average into through-strings and caps.
- **Measured**: the F1 checks, and the cap share of old content at n = 1024 (one network so far).
- **Synthesis**: the Brauer/Jones reading of the fresh-weight lemma; the identification of caps with the coincident pattern, column sums, trace channel and dilation charge; He-criticality as unit loop eigenvalue of the cap sector.
- **Conjecture**: F2 (noncommutative entropic trickle-down under one-site Pimsner–Popa bounds).
- **Speculation**: that the leaders' memory carrier is "cap sector pinned + bulk merged".
- **Risks.**
  - F1 may be folklore (being checked).
  - F2 may need an additive noncommutativity correction (BCR's d).
  - On the competition side, the measured bulk (≈ 15–20 % of old energy, long-lived, high-rank) is exactly the part that no sector argument carries. If T2's merge fails, the local-to-global lens explains ≈ 80 % of memory and leaves the decisive 20 % to brute force or to error-budget shifting.

## 7. Synthesis across teams A–E

To come in the final (from ≈ 22:30 UTC).

## 8. Sources

- Oppenheim, arXiv:1709.04431; Alev–Lau, arXiv:2001.02827; Anari–Liu–Oveis Gharan, arXiv:2001.00303; Kaufman–Oppenheim, arXiv:1707.02799.
- Chen–Liu–Vigoda, arXiv:2011.02075; Chen–Eldan, arXiv:2203.04163.
- Bardet–Capel–Rouzé, arXiv:2001.07981; Chen–Rouzé, arXiv:2504.02208; Yang, arXiv:2609.38007; Chen–Kastoryano–Gilyén, arXiv:2311.09207.
- Jones, *Index for subfactors* (1983); Pimsner–Popa (1986).
- Repo: the bridge digests and streams listed at the top.

- Literature probe (abstracts and main theorems):
  - Gao–Junge–LaRacuente arXiv:1710.10038 (Cor 2.3) and 1909.01906; LaRacuente arXiv:1912.00983;
  - arXiv:2505.08991; Hayakawa et al. arXiv:2609.39802;
  - Bauerschmidt–Bodineau arXiv:1712.03676; Eldan–Koehler–Zeitouni arXiv:2007.08200; Boban–Li–Oveis Gharan arXiv:2609.13138; Mikulincer–Sohn arXiv:2512.22803; Anari–Koehler–Vuong arXiv:2407.16104.

## Files

- `nc_trickle.py`: Theorem F1 check.
- `trace_share.py`: cap-sector share of old content at n = 1024 (uses costate's code read-only).
- `results/`: outputs.
