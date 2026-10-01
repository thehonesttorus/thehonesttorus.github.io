# Team F: the noncommutative local-to-global framework — REPORT (final)

*1 Oct 2026, v1 final (synthesis included). Scope: Garland / Oppenheim / Alev–Lau, spectral and entropic independence, Chen–Eldan localization, commuting squares and quasi-factorization; a precise noncommutative trickle-down statement; the dictionary to the competition object; the synthesis across teams A–E (§7). Labels: **Theorem** (proved here, elementary, and checked numerically), **Fact** (cited), **Measured** (run here, numbers in `results/`), **Synthesis**, **Conjecture**, **Speculation**.*

*Builds on, does not repeat: `notes/local-to-global-unlocks.md` §4–5, `notes/digests/bridges/hdx-spectral-independence.md` (§2 trickle-down, §4.5 BCR, §7 "no k-step quantum trickle-down found"), `notes/streams/bridges-synthesis/A-unification.md` (D-2 two-projection inequality, D-3 "SI is a frame bound", C-2 the missing quantum k-step trickle-down), `C-competition.md` (R1 doubled leg W∘W = (2/n)11ᵀ + E, R11 no transport gap, "no direct lever").*

## 0. Summary

1. **The spectral trickle-down has no noncommutative obstruction once it is written in the forget lattice** (Theorem F1, §2.1). For any lattice of KMS-symmetric conditional expectations E_A ("forget A"), closed under intersection of ranges, the heat-bath dynamics (1/n)Σ_i E_i has gap ≥ ∏_{m=2}^{n} (1 − (1+η_m)/m), where η_m + 1 is the frame bound ‖Σ_{i∈A}(E_{A∖i} − E_A)‖ at level |A| = m, and η_m ≤ λ_max of the matrix of pairwise commuting-square defects ‖E_{A∖i}E_{A∖j} − E_A‖. Pinning is never used: the operator norm on the whole L² is the supremum over pinnings classically, and N_A-bimodularity replaces disintegration noncommutatively. The classical case is exactly Anari–Liu–Oveis Gharan Thm 1.3 / Alev–Lau. Checked on a non-commuting 4- and 5-qubit Gibbs ring: the bound holds, is exact at β = 0 and within 1.3–1.7× at β = 0.3.
2. **The noncommutative failure is entropic, and it is one step** (§2.2). Classical local-to-global for entropy (Chen–Liu–Vigoda Thm 5.6) converts spectral to entropic contraction *inside a link*, where the densities are single-site marginals with bounded ratio b. With no links, the noncommutative comparison must be made on a whole region, at a cost that grows like e^{c|A|} (the χ_KMS growth already measured in A-unification C10–C11). **Repair hypothesis:** the Pimsner–Popa constant of the one-site inclusion N_A ⊂ N_{A∖i} (a Jones index, O(d²), not exponential) plays the role of b. Conjecture F2 states the resulting iterable approximate tensorization; its two known corners (m = 2: Bardet–Capel–Rouzé; exact commuting lattice: Lemma F2a, proved here) are both consistent with it.
3. **Dictionary, precise (§3).** The fresh-weight average E_W over rotation-invariant Gaussian weights acts on tensors through the **Brauer algebra**: a sum over pairings. The *through-string* pairings are the commuting-square part: they give the Frobenius currency and the incoherent sum over layers. The *caps* are Jones projections with loop value n (normalised projection e/n, Markov trace τ = 1/n², index n²): they are partial traces, i.e. the coincident-pattern channel, the column sums of N7, the trace channel and the dilation charge. He initialisation puts the cap sector exactly at eigenvalue 1 (loop n × variance 2/n × gate density ½ = 1), a unipotent sector with no mixing (angle 0). By Proposition F3 an angle-0 sector cannot be quasi-factorized; it must be **pinned** (localized), which is Chen–Eldan / Bauerschmidt–Bodineau's move.
4. **Measured at n = 1024 (§5.2, three networks).** The cap (trace) channel holds 0.66–0.79 of the energy of the old (2,1) slice at layers ≥ 8. With costate's scale mode the union holds 0.82 at layer 15 (random-direction control ≤ 0.005). The remaining through-string bulk (18 % of old energy at A = 1, still 12.5 % beyond age 3 and 3.5 % beyond age 7) forgets only slightly faster than the cap sector. Priced by the ε² law, that bulk costs ≈ 3e-7 at A = 3. So the cap sector is most of the memory, not all of it, and the bulk is the decisive remainder.
5. **Transfers, measured (§4, §5.4).**
   - Pinning the cap sector (T1) lands at 7.6–7.7e-7 on all three networks at first order (young A = 1 exact). That is 60–70 % of the old-content gain, the same cap that team A's no-go and team D's Perron oracle hit.
   - The bulk is readout-relevant and needs rank ≈ 256 per layer (T2).
   - Neither is a carrier at one layer's cost.
6. **Synthesis (§7).**
   - *The mechanism.* All of A–E's transfers factorize one conditional expectation, the Brauer/Weingarten fresh-weight average E_W, which has exactly two sectors:
     - the cap sector: angle 0, unipotent, must be pinned;
     - the through-string free sector: per-step angle g_l = 2τ(P_x P_{x′}), a two-projection trace, with energy g³ per step and rank n/(2(a+1)); it tensorizes by Theorem F1 at η = 0.
   - *Agreement.* Five teams' numbers agree on this split: cap share 0.75–0.82, free decay 1.02–1.10·g³, rank law ≈ 2n/a.
   - *The theorem shape for memory at one layer's cost.* An approximate commuting square between the carrier and the transport, iterable for free over independent layers. Linear carriers face a conjectured harmonic floor Ω(n³ log L) (Conjecture F4).
   - *The live escape.* Node sufficiency, from Bethe's oracle: memory is needed only through 3 n-vectors per layer. Add law-level cap pinning, plus free-sector cores merged without fit in a shared basis per age bin (Conjecture F5).
   - *Measured (S1, §5.5).* The shared-basis cores form an exact commuting square with transport. Causally they carry **all content of ages ≥ 8 in 64 dimensions at excess ≈ 5e-9 over exact first order, for ≈ 1–2 u per layer**. Deep memory is cheap; the bill is ages 1–7.

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

**Proposition F2b (linear response: the entropic frame constant is never better than the spectral one; proved here).**

*Setting.* The E_A are σ-preserving conditional expectations, so each N_A is modular-invariant and E_A commutes with the modular operator Δ_σ (Takesaki).

*Claim.* Then lim_{t→0} sup_X [Σ_i D(E_{A∖i}*ρ_t‖E_A*ρ_t)] / D(ρ_t‖E_A*ρ_t) = 1 + η_m, with ρ_t = σ + tX. Hence η^{ent}_m ≥ η_m, with equality in linear response.

*Proof.* To second order, D(σ + tX ‖ σ + tY) = (t²/2) ‖X − Y‖²_{BKM,σ}. The dual maps E_A* are orthogonal projections in the BKM metric, because they commute with Δ_σ. So both sides become quadratic forms of the same nested projections, and the ratio's supremum is ‖Σ_i P^A_i‖ in the BKM metric. Σ_i P^A_i commutes with Δ_σ and is self-adjoint in every metric of the form ⟨·, f(Δ_σ)·⟩, so its norm, which equals its spectral radius, is the same as in the KMS metric. ∎

*Reading.* Conjecture F2 is therefore a statement about the *non-perturbative* regime only. Near equilibrium the entropic and spectral local-to-global theories coincide exactly in the noncommutative setting. What H-PP has to control is how far from σ the ratio can drift.

### 2.4 Proposition F3 (a central defect must be pinned, not factorized)

**Proposition F3.** Let P, Q be orthogonal projections with R = P ∧ Q, and let Z be a projection commuting with P and Q such that PQ − R = Z(PQ − R)Z. Then:
- on Ran(1 − Z), P and Q commute and ‖(1 − Z)(f − Rf)‖² ≤ ‖(1 − Z)(f − Pf)‖² + ‖(1 − Z)(f − Qf)‖² exactly;
- on Ran Z the constant is 1/(1 − c_Z), with c_Z = ‖Z(PQ − R)‖.

If c_Z = 1 (an angle-0 sector, i.e. a conserved quantity), no factorization constant exists on Z. The only way to recover a global statement is to **condition on the Z-sector**, pinning its value and averaging over its law. Then each conditional problem lives on Ran(1 − Z), where the square commutes.

*Proof.* Z commutes with P, Q and R, so everything splits into the two blocks. The first block is the commuting-square case; the second is A-unification D-2. ∎

This is the abstract reason that localization (Chen–Eldan, Bauerschmidt–Bodineau) and not factorization is the right tool for conserved or critical sectors. It is the shape used in T1.

### 2.5 G5: Theorem F1 at the wall (Dixmier-critical trickle-down)

Write ε_m := η_m/(m − 1) for the *local excess* at level m (m free pieces). Theorem F1(c) reads n·gap ≥ ∏_{m=2}^{n}(1 − ε_m). Since log(1 − x) ≥ −x − x² for 0 ≤ x ≤ ½:

**Corollary F1-crit (proved).** If ε_m ≤ ½ for all m, then

  n · gap ≥ exp(−Σ_{m≤n} ε_m − Σ_{m≤n} ε_m²).

The global gap is therefore governed by **one number, the divergence of Σ_m ε_m**. This gives three regimes, which are exactly the three regimes of coordinator note 3.

| regime | local excess ε_m | global gap | name |
|---|---|---|---|
| gapped | Σ ε_m < ∞ (e.g. η_m ≲ m^{−δ}, correlations decaying faster than the free region grows) | n·gap ≥ c > 0: optimal Θ(1/n) | geometric / summable |
| **Dixmier-critical** | ε_m ~ η̄/m, i.e. η_m → η̄ (bounded spectral independence, local gaps closing like 1/m) | n·gap = n^{−η̄ + o(1)} | **the exponent of polynomial loss is the Dixmier trace** Tr_ω(diag ε_m) = lim (1/log n) Σ_{m≤n} ε_m |
| sub-critical log | ε_m ~ c/(m log m) | n·gap ≳ (log n)^{−c} | the global gap closes only logarithmically |

**Readings.**
- **Wall ≡ L^{1,∞}.** In the Dixmier-critical row, the operator D_ε = diag(ε_m) lies in the weak trace class L^{1,∞}. The global loss is its Dixmier trace, which ignores every finite set of levels and sees only the scale-invariant tail. This is exactly the structure coordinator note 2 found for the age-graded memory: k(a) = 2n/a, total Σ_a k(a) ≈ 2n log L.
  - So the answer to "does F1 have a Dixmier-critical analogue" is **yes, and it is F1 itself, read through the Dixmier trace**. Bounded spectral independence (the standard ALO hypothesis) *is* the critical case.
  - The uniform-η form gives n^{−sup η}. The level-resolved form gives n^{−Tr_ω(ε)}, with Tr_ω(ε) = the log-average of η_m. That average is never larger than the sup and can be much smaller when η_m is large only on a bounded range of scales.
- **Logarithmic closing** needs η_m ≈ c/log m. That means spectral independence must *improve* slowly with the size of the free region, which is the signature of a system sitting exactly at a threshold where correlation decay is marginal.
- **Sharpness.** At β = 0 (product state) the product formula is exact (§2.1 table, 0.2500 = 0.2500), and Alev–Lau Prop 3.3 gives tightness in the complex setting. So the three regimes are not artefacts of the bound.
- **Noncommutative validity.** Unchanged: the corollary is pure algebra on Theorem F1, which holds for any lattice of KMS-symmetric conditional expectations.

**Which sampling or counting problem it would improve (Speculation, with the hypothesis made explicit).** The level-resolved exponent pays off at a **uniqueness threshold**, where the uniform bound sup_m η_m degrades but its log-average need not. The natural targets are:
- the hardcore model at λ_c(Δ);
- antiferromagnetic two-spin systems at the tree-uniqueness boundary on bounded-degree graphs;
- the quantum analogue, Gibbs samplers at the threshold temperature of a commuting or weakly non-commuting Hamiltonian, through the noncommutative F1.

If at the threshold the pinned-region influence bound behaves like η_m ≤ C + c log m (influence sums truncated by the free-region size), the uniform form gives nothing finite, while F1-crit gives n·gap ≥ n^{−C} e^{−(c/2) log² n}, i.e. quasi-polynomial. If η_m is bounded on average (finite log-mean), F1-crit gives a polynomial bound with the averaged exponent. Annealed counting (Štefankovič–Vempala–Vigoda) inherits whichever bound holds. **What must be checked** before any claim is the level-resolved influence profile η_m at the threshold. Recent threshold-mixing results should be read with this in mind; their constants were not verified here.

**The competition reading.** Our free sector is Dixmier-critical in the age variable (rank 2n/a, team D; energy g³ per step, team C), and the measured carrier pays log L in leg dimension. S1 (§5.5) shows the same phenomenon from the carrier side: the shared rank needed falls faster than 2n/a once the content is old (64 dimensions for all ages ≥ 8). In F1-crit language, the old levels have *summable* excess, and only ages ≲ 8 are critical.

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
- **Status (final).** Measured at Stein level with young A = 1 exact, the cap sector gives 7.6–7.7e-7 on three networks (§5.4). Team A's Bethe-sector FC gives 1.05e-6, and team D's Perron sandwich 1.5e-6. At law level, costate's A3gsl_nc gives 3.25e-7 and team C's trace-adaptive window 2.43e-7. **Alive as a component, dead as a carrier.**

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
- **Status (final).**
  - In the readout, the bulk needs rank ≈ 256 (§5.4).
  - Region's CP merge of all ages > 2 to R = n gives 2.0e-7 at 3,510 u: killed on cost (team E T4).
  - Team D's age multiresolution (each source in 2n/a directions of its own propagator) is lossless as an oracle at O(n³ log L).
  - T2 as a CP merge is dead. As fit-free merging in a shared basis it survives as Conjecture F5, decided by S1 (§5.5).

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


### 5.3 Is the bulk random? Rank profile of the bulk slice (`bulk_rank.py`, MLP 0)

The bulk is the old off-diagonal slice minus the cap channel and minus the scale mode. The table gives top-k singular-value energy shares, k = 1 / 16 / 64 / 256 / 512.

| object | layer 6 | layer 10 | layer 15 |
|---|---|---|---|
| i.i.d. Gaussian n × n control | 0.004 / 0.058 / 0.210 / 0.623 / 0.894 | (same) | (same) |
| whole old slice, A = 1 | 0.64 / 0.77 / 0.91 / 0.993 / 0.999 | 0.74 / 0.87 / 0.96 / 0.996 / 0.999 | 0.79 / 0.92 / 0.98 / 0.997 / 0.999 |
| bulk, A = 1 (0.33 / 0.21 / 0.16 of old) | 0.03 / 0.35 / 0.75 / 0.985 / 0.999 | 0.05 / 0.46 / 0.82 / 0.989 / 0.999 | 0.07 / 0.56 / 0.88 / 0.991 / 0.999 |
| bulk, A = 3 (0.38 / 0.17 / 0.13 of old) | 0.03 / 0.37 / 0.78 / 0.992 / 0.999 | 0.06 / 0.49 / 0.86 / 0.993 / 0.999 | 0.07 / 0.56 / 0.90 / 0.993 / 0.999 |

**Reading.** The bulk is *not* random as a matrix: 64 directions hold 75–90 % (control 21 %), 256 hold ≥ 98.5 % (control 62 %). Its symmetric and antisymmetric parts carry equal energy (0.500), as a generic non-symmetric slice would. So the bulk is high-rank compared with the cap sector (rank ≲ 4 plus one dense template) but low-rank compared with n. That matches team D's age-dependent rank law (coordinator note, correction): content of age a lives in about 2n/a propagator directions.

### 5.4 What the readout needs: oracle filters on the old slice inside costate's first-order co-state (`cap_filter.py`, A = 1, n = 1024)

Each filter keeps the exact diagonal and a part of the old off-diagonal slice, injected at every layer. Raw = final MSE − truth noise.

| filter on old content (age > 1) | MLP 0 | MLP 1 | MLP 2 |
|---|---|---|---|
| none (diagonal only) | 1.36e-6 | 1.59e-6 | 1.46e-6 |
| cap sector (cap channel ∪ scale mode) | 7.74e-7 | 7.62e-7 | 7.64e-7 |
| cap + bulk rank 16 | 6.82e-7 | — | — |
| cap + bulk rank 64 | 5.09e-7 | 5.21e-7 | 4.77e-7 |
| cap + bulk rank 256 | 4.20e-7 | 4.72e-7 | (running) |
| plain rank 64 of the whole off-diagonal (costate's filter) | 4.84e-7 | 4.72e-7 | (running) |
| all (exact first order) | 3.92e-7 | 4.03e-7 | 4.13e-7 |

**Reading.**
- In the readout metric the cap sector recovers 60 % / 70 % / 66 % of the old-content gain on MLPs 0 / 1 / 2 (e.g. (1.36 − 0.774)/(1.36 − 0.392) on MLP 0), against ≈ 80 % of its energy. It lands at 7.6–7.7e-7 on all three networks.
- The bulk is readout-relevant. Rank 64 recovers 69 % of the remaining gain and rank 256 recovers 93 %.
- **The cap channel is not the readout-optimal subspace.** A plain rank-64 SVD (4.84e-7) beats cap + bulk rank 64 (5.09e-7). The readout reads old content through a ≈ 64–256-dimensional subspace that only partly aligns with the O(n)-invariant cap directions.
- Consequence for T1: pinning the cap sector alone, at Stein level, would at best take the A = 1 co-state from 1.36e-6 to ≈ 7.7e-7. Costate's law-level treatment of the same sector (A3gsl_nc, 3.25e-7) does better than Stein level, which is the evidence that *pinning at law level* (Proposition F3), rather than first-order projection, is what pays.

### 5.5 Test S1: one shared basis per dyadic age bin, so that cores merge without fit (`s1_shared.py`, `costate_f.py`, MLP 0, n = 1024)

**Construction.** Inside costate's exact first-order co-state, every source of age ≥ a_min is restricted to a k-dimensional target-space basis that is **shared by all sources in its dyadic age bin** [2^j, 2^{j+1}). The oracle takes the top-k eigenvectors of Σ_{s∈bin} U_sᵀU_s at each layer. For comparison, "own" restricts each source to its own top-k directions (costate's rank hook).

**Why sharing matters (theorem-level reason).**
- Once all members of a bin live in one basis B (n × k), each (2,1)-slice term of each source has the form Σ_{ijl} G_{ijl} B_{pi}B_{pj}B_{ql}. The core G (k³) is a *sum over members with no fit*, and it does not change under transport: transport acts only on B, as B ← orth(Mᵀ B), with the k × k change of basis applied to G.
- So the bin's projection P_bin and the transport form an **exact commuting square**, T_t P_bin(t) = P_bin(t+1) T_t. The only defect is the one-time projection at bin entry (Theorem F1 at η = 0 then sums these entry defects incoherently).
- This is the "fit-free linear re-binning that commutes with transport" that team E §4.9 left as its open hatch.

| variant (ages ≥ a_min restricted) | raw | ratio to exact |
|---|---|---|
| exact first order, all pairs | 3.92e-7 | 1 |
| shared k = 256, a_min = 4 | 3.97e-7 | 1.01 |
| shared k = 128, a_min = 4 | 4.92e-7 | 1.25 |
| own k = 128, a_min = 4 | 4.91e-7 | 1.25 |
| shared k = 64, a_min = 4 | 6.39e-7 | 1.63 |
| **shared k = 64, a_min = 8** | **3.94e-7** | **1.01** |
| shared k = 2n/a, a_min = 4 | 3.93e-7 | 1.00 |
| *causal* cohort k = 64, a_min = 8, one cohort per 8 layers | **3.97e-7** | **1.01** |
| *causal* cohort k = 32, a_min = 8 | 4.27e-7 | 1.09 |
| *causal* cohort k = 256, a_min = 4, cohorts of 4 layers | 4.12e-7 | 1.05 |
| *causal* cohort k = 128, a_min = 4 | 5.51e-7 | 1.40 |
| MLP 1: exact / shared k = 64, a_min = 8 | 4.03e-7 / 4.11e-7 | 1 / 1.02 |
| MLP 1: causal k = 64, a_min = 8 / causal k = 32, a_min = 8 / causal k = 256, a_min = 4 | 4.18e-7 / 4.37e-7 / 4.03e-7 | 1.04 / 1.08 / 1.00 |

*Causal* means: each source is projected once, when it reaches age a_min, onto its own top-k target directions. It joins the open cohort, whose older members are re-projected onto the newest member's basis (a k × k change on the cores). Nothing else is ever fitted, and transport is exact afterwards. Building the basis needs only a rank-k range finder of one n × n factor, O(k n²), not an oracle.

**Readings.**
- **Sharing is free.** One basis per dyadic bin is as good as each source's own basis (4.92e-7 vs 4.91e-7 at k = 128).
- **Deep memory is small.** All content of ages ≥ 8 is read losslessly through one shared 64-dimensional basis (n/16). Ages 4–7 need ≈ 256 (n/4). This decays faster than the per-source law 2n/a (team D), which would give 256 at age 8: sharing plus g³ decay compress further.
- **Cost of the merged-core readout.** It is 2nk³ + 2n²k per bin per layer:
  - k = 64: 0.27 u;
  - k = 128: 2.1 u;
  - k = 256: 16.6 u.

  Basis transport is 2n²k (≤ 0.25 u). Building a core costs ≈ 2nk³ per term once per entering source.

  So bin [8, 16) costs ≈ 0.3–1 u per layer, essentially free. Bin [4, 8) at k = 256 costs more than its 4 exact pairs would. **The bill sits in ages 1–7, not in deep memory.**
- **The causal version holds up.** Deep memory (ages ≥ 8) is carried causally at k = 64 with excess ≈ 5e-9 over exact first order, and at k = 32 with ≈ 3.5e-8. Ages ≥ 4 at k = 256 carry an excess of ≈ 2e-8.
- **Accounting at n = 1024, L = 16.** FC's pairs with age ≥ 8 are 28 of 120, ≈ 196 u of its 840 u. The causal cohort carrier replaces them by ≈ 1–2 u per layer:
  - one core build per entering source, ≈ 6 terms × 2nk³ ≈ 1.6 u at k = 64;
  - the readout, ≈ 0.3 u;
  - the basis transport, ≈ 0.1 u.

  That saves ≈ 20 % of FC at no measurable accuracy cost.
- **The real bill.** Ages 1–7 (≈ 92 pairs, ≈ 640 u) stay at near-full rank. So at this width and depth, **"memory" in the deep sense is cheap. The binding object is the mid-range content of ages 1–7**, whose free sector has not yet decayed (g³ ≈ 0.6 per step at depth).
- **Caveat.** These are first-order (Stein-level) co-state numbers on one network (MLP 1 in `results/s1_w1024_mlp1.jsonl`). The excess that matters is absolute, since FC sits at 3e-8. The a_min = 8, k = 64 excess (≈ 2–5e-9) is below the bar's resolution. The a_min = 4, k = 128 excess (1.6e-7) is not.

## 6. Honest assessment

- **Theorem** (proved here, elementary, checked):
  - F1 (spectral k-level trickle-down in any forget lattice);
  - Lemma F2a (entropic, exact commuting lattice);
  - Proposition F2b (the entropic frame constant equals the spectral one in linear response);
  - Proposition F3 (a central defect must be pinned);
  - the Isserlis split of the fresh-weight average into through-strings and caps;
  - the exact commuting square between a transported shared basis and transport, which is what makes merged cores fit-free (§5.5).
- **Measured** (n = 1024, bench w1024_d16):
  - cap share of old (2,1)-slice energy, 0.82 at depth on three networks;
  - the age profile of cap and bulk;
  - the bulk rank profile;
  - readout filters: cap 7.6–7.7e-7 on three networks, bulk needing rank ≈ 256;
  - S1: shared-basis cores, causal, deep memory at k = 64 with excess ≈ 5e-9 (MLP 0; MLP 1 in `results/`).

  All competition measurements are at first-order (Stein) level inside costate's co-state, not inside FC.
- **Synthesis**:
  - the Brauer/Jones reading of the fresh-weight lemma;
  - the two-sector mechanism uniting A–E;
  - the consistency ledger;
  - the approximate-commuting-square formulation of a carrier.
- **Conjecture**:
  - F2 (noncommutative entropic trickle-down under one-site Pimsner–Popa bounds);
  - F4 (harmonic floor for linear carriers);
  - F5 (node-level iterable tensorization).
- **Speculation**: how the leaders carry ages 1–7.
- **Risks.**
  - F1 is probably folklore. The literature probe found no k-level quantum statement, but that is not proof.
  - F2 may need BCR's additive noncommutativity correction.
  - S1's excess is measured against first order (4e-7 floor). Inside FC (3e-8) the k = 64 deep-memory excess must be re-measured.
  - F4 is supported by every linear carrier measured, but is unproven.

## 7. Synthesis across teams A–E

*Sources: the five final REPORTs in `notes/essence/*/`, `COORDINATOR-NOTE-1.md` with its 21:10 correction, and F's own measurements. Section numbers refer to the teams' reports.*

### 7.1 One noncommutative mechanism behind the transfers

**Statement (Synthesis, built from proved pieces).** Every competition transfer the five teams made is an attempt to factorize one conditional expectation: the **fresh-weight average E_W**, the Weingarten/Brauer average over a rotation-invariant Gaussian layer. It has an exact two-sector structure.

1. **The cap sector.** Brauer caps (Jones projections, loop value n) are partial traces. This is the O(n)-invariant part: dilation, trace channel, coincident patterns, the Bethe/diagonal pairing, the Perron mode.
   - Its transfer eigenvalue is exactly 1: loop n × variance 2/n × gate density ½ = 1. This matches Euler's identity and the costate C6 Jordan block, and team C's Proposition C1, which finds a 2×2 unipotent block.
   - So its angle is 0. By Proposition F3 it cannot be factorized; it must be **pinned** (localized).
2. **The through-string ("free") sector.** Its pricing is tracial and Frobenius:
   - team D Proposition D1;
   - team B §2.2 (Schur orthogonality);
   - region §1.

   Its contents are mutually orthogonal across ages and layers (cross-age cosines ≤ 0.03, team C §4). Its per-step angle is a **two-projection angle**, g_l = 2E[Φ_l²] = 2τ(P_x P_{x′}) = ½ + arcsin(ρ_l)/π: the trace of the two replica gate projections at angle θ_l = arccos ρ_l (team E §2.2, team C). Energy contracts like g_l³, one factor per leg, at 1.02–1.10 g³ (team C), and rank falls like the participation ratio n/(2(a+1)).

   By Theorem F1 with η = 0 on the lattice of independent weight layers, local errors in this sector **tensorize**: they add, they do not compound. That is the books-close identity and region's K(l) price list, reproduced from first principles by team E's wedge calculus to ±0.1 in log.

**Each team's transfers in this language.**

| team | transfer | sector it acts on | role of E_W | measured verdict |
|---|---|---|---|---|
| A | Bethe / diagonal pairing (T1), gauge-group lattice | cap: the (Z₂)ⁿ / O(n) average is the commuting-square lattice | down step = gauge average | the cap carries ≈ 80 % of memory energy but leaves 0.40–0.48 D21 error, ≥ 1.05e-6 end-to-end. **No-go for gauge-invariant carriers** |
| A | Barvinok / rank-k fluctuation (T2), Godsil–Gutman (T3) | free sector, sampled or truncated | sketching the through-strings | killed: rank 256 still 8–10 %; GG variance ≥ 10⁵× the budget |
| B | orthogonal hybrid (T1), Schur orthogonality (§2.2) | free sector: martingale increments over F_l | E_W as Doob conditional expectation | derivation = Theorem F1 at η = 0. Cancellation test not run |
| B | JMR squared graph (T2), INW seed (T5) | free sector (λ_eff ≈ 1, no contraction in the quenched propagator) | ignore-first-step needs a gap that is absent | killed / derived negative |
| C | two-sector trichotomy, trace-adaptive window (T2) | both: parabolic cap + free g³ | the angle g_l sets the window | **only measured gain**: 2.43e-7 (0.75× A3gsl_nc), 6/6 networks |
| D | D1 tracial price; Kyng–Sachdeva (3.1); Perron sandwich (3.3); age multiresolution (3.7) | free sector (D1, 3.1, 3.7); cap (3.3) | E_W = the tracial state (fresh consumer) | 3.1 and 3.3 killed. **3.7 oracle lossless**: k(a) = 2n/a, 3.31e-8 vs FC 3.24e-8, O(n³ log L). Causal version not run |
| E | wedge calculus (T1), age law (T2), dyadic merging (T4) | the price of both sectors from θ_l; free-sector age law | E_W over the fresh W = the down step | T1 positive (K(l) to ±0.1 log); T4 killed on cost (CP merge R = n: 2.0e-7 at 3,510 u) |
| F | pin the cap (T1); iterable bulk (T2) | cap / free | Proposition F3 / Theorem F1 at η = 0 | cap at Stein level 7.6–7.7e-7 (60–70 % of the old gain); bulk needs rank ≈ 256 in the readout |

**Consistency ledger.** Independent measurements of the same sector agree.

| quantity | value | measured by |
|---|---|---|
| cap/dilation share of old (or all-age) D21 energy at the last layer | 0.82 (C, all ages), 0.82 (F, union, A = 1, three networks), ≈ 0.8 (A, Bethe sector), 0.75 (D, Perron second leg at t = 15) | C §4, F §5.2, A §4, D §4.2 |
| cap-only end-to-end | 7.6–7.7e-7 (F, Stein level, young exact at A = 1); 1.05e-6 (A, FC with memory replaced by its Bethe value at A = 3); 1.5e-6 (D, Perron q = 1) | F §5.4, A T1, D 3.3 |
| free-sector decay per step | (1.02–1.10)·g³, g³ ≈ 0.13 → 0.6 with depth (C); tail ratio ≈ 0.7–0.77 per age at layer 15 (F) | C §4, F §5.2 |
| free-sector rank | PR = n/(2(a+1)) (C, E); resolution ≈ 4·PR = 2n/a lossless (D); bulk slice 64 → 75–90 %, 256 → ≥ 98.5 % (F); 121–191 modes at 10 % (C §6.2) | C, D, E, F |
| pricing | K(l) from θ_l to ±0.1 log (E); D1 exact (D); books close 2–6 % (region) | E, D, region |

**What the mechanism explains among the established facts (brief §3).**
- "Old content nearly orthogonal to the present and to other ages" (F8.3) is through-string orthogonality: free-sector contents of different ages are Brauer-orthogonal, and the cap part is coherent across ages. D's mean cosine between old sources, 0.38–0.46, is consistent with coming from the cap (inference, not separately measured).
- "Every age matters" (N6) has two causes: angle 0 in the cap, and a slow free angle (g³ ≈ 0.6 at depth).
- "Ensemble mean zero, purely quenched" (F8.5): the free sector has zero E_W-image by construction, so every annealed or gauge-averaged carrier, i.e. every conditional expectation onto invariants, misses it (A's no-go, E's Conjecture E3).
- "The leaders carry memory at about zero cost" is **not** explained. §7.2 says exactly what would explain it.

### 7.2 The theorem shape that would give a carrier of memory at about one layer's cost

**Requirement.** A carrier is a family of states s_t with an update s_{t+1} = Φ_t(s_t, W_{t+1}) costing O(n³) (≈ 1–4 products), plus a readout of the per-layer quantities the means need (region §4: D21 to 5–10 %, node κ3/κ4 to ≈ 7 % at layers ≥ 5).

By Theorem F1 at η = 0 (exact independence of the weight layers), the global error is the incoherent sum Σ_t K(t)·δ_t², where δ_t is the per-layer **commuting-square defect** between the carrier and the true transport:

  δ_t := ‖R_{t+1}(Φ_t(s_t)) − T_t(R_t(s_t))‖_read,

with R_t the readout map and T_t the exact (first-order chaos) transport. So the theorem to look for is an **approximate commuting square for the carrier**,

  R_{t+1} ∘ Φ_t ≈ T_t ∘ R_t, with defect δ_t ≤ 5–10 % at layers 6–13 (team E §4.8 gives the per-layer tolerance),

*iterable for free* because the weight lattice commutes exactly. This is Theorem F1's structure: local defects are global errors, with no compounding. The whole difficulty is the local defect.

**Two measured constraints on the defect.**
- (i) **Cap.** The cap sector has angle 0, so a carrier must hold it exactly (Proposition F3). Holding it at **law level**, i.e. pinning the dilation/collective coordinate, does better than holding it at first order. Costate's law-level A3gsl_nc at 3.25e-7 beats its exact all-pairs first order at 4.0e-7 with 32 % of the products, and C's trace-adaptive window takes it to 2.43e-7.
- (ii) **Free sector — the harmonic floor (Conjecture F4).** Content of age a needs a resolution of ≈ 2n/a directions (D's oracle). It is orthogonal across ages (C) and quenched with zero E_W-image (F8.5). The conjecture is that any carrier that is a **linear image of the first-order chaos content** (the R_t of a conditional expectation onto a subspace of atoms or legs) has defect ≤ 10 % only if its state dimension per target is ≳ n Σ_{a ≤ a*} 2/a ≈ 2n ln t. That gives cost Ω(n³ log L) per layer, which is D's oracle cost.

  *Status:* a conjecture, supported by every linear carrier measured:
  - shared subspaces: q = n/4 is 6× worse;
  - per-source bases;
  - pruning;
  - Kyng–Sachdeva;
  - CP merge (R = n reaches only 2.0e-7);
  - F's bulk ranks.

**Consequence.** One-layer cost requires a carrier that is **not a linear image** of the first-order content. The local-to-global literature offers exactly three non-linear escape shapes, and the measurements say which are live.

1. **Localization (Chen–Eldan / Bauerschmidt–Bodineau), applied to the cap sector.** Pin a few collective coordinates and run memoryless closures on the conditional laws. This is live and cheap (O(K n²) per layer), but by the measurements it carries ≲ 70 % of the old gain at first order. Its law-level advantage is real (3.25e-7, 2.43e-7) but still ≈ 25× short.
2. **Node sufficiency (Bethe / commuting square of the pair algebra over the node algebra).** Bethe's deciding experiment (designs/bethe DESIGN §14) shows that v4's memoryless pair, mean and readout machinery, given the *true node beliefs* (v, κ3, κ4 per neuron) at layers ≥ 6, reaches ≤ 3.9e-8 (MC noise floor), from 4.3e-7. In NC terms: given the node algebra, the pair algebra and the past form an approximate commuting square, E_past E_pair ≈ E_node, with defect below 4e-8.
   - So **memory is needed only through 3 n-vectors per layer**.
   - The memory problem reduces to producing per-neuron (v, κ3, κ4) at depth to ≈ 7 %. This is a *diagonal* readout of the free sector: node κ3 = Σ_s 1ᵀ diag(w2_s)(Y_s ∘ Y_s ∘ Z_s) at the diagonal. In core form it is a contraction G_{ii′j} B_{ai}B_{ai′}D_{aj} in which **cores of the same age bin, expressed in a shared transported basis, add without any fit**: the "fit-free linear re-binning that commutes with transport" that team E left as its one open hatch.
3. **Inter-layer cancellation** (team E hatch 4, B T1, F T3). The incoherent sum is only an upper bound when the frame matrix ρ_{ll′} has negative eigenvalues. MLP 1's 46 % over-prediction is such a case. The leaders could exploit it, but this cannot be engineered without knowing the sign structure. Low prior.

**The precise theorem shape F proposes (Conjecture F5, node-level iterable tensorization).** For He-ReLU networks at n → ∞ with L/n → 0:
- (a) **Node sufficiency.** The final-layer means depend on the past only through the node beliefs (v_t, κ3_t, κ4_t) ∈ (ℝⁿ)³ of each layer, up to a defect o(1/n²) per layer in the readout metric. This is Bethe's oracle, as a theorem.
- (b) **Sector split of node beliefs.** Node beliefs split as (cap part, computable at law level from the pinned collective coordinates at O(n²) per layer) + (free part). The free part's node-diagonal is a sum over age bins of core contractions in a shared QR-transported basis, merged without fit.
- (c) **Iterability.** By Theorem F1 at η = 0, per-layer defects in (a) and (b) add incoherently.

**Measured status of (b) (§5.5, test S1).**
- **Fit-free merged cores in a shared, transported basis are measured to work, causally.** All content of ages ≥ 8 is carried in one 64-dimensional cohort basis per bin with excess ≈ 5e-9 over exact first order, at ≈ 1–2 u per layer.
- Ages ≥ 4 need k ≈ 256 (excess ≈ 2e-8 causal), where the n k³ core readout (≈ 16 u per bin per layer) is no cheaper than exact pairs.
- **So the theorem shape closes on cost for deep memory and not for ages 1–7.** At n = 1024, L = 16 the remaining bill (≈ 640 u in FC's accounting) is mid-range content whose free sector has not yet decayed. A one-layer-cost estimator therefore needs, in addition, a cheaper representation of ages 1–7. The candidates are:
  - node sufficiency (a), which reduces what must be read to 3 n-vectors per layer but not the transport;
  - the public chain's young tier;
  - a lower-rank readout of bin [4, 8) (HOSVD of its core, untested).

### 7.3 Deciding tests the synthesis implies (cheapest first)

- **S1 (shared basis per age bin; run, §5.5).** It was run end-to-end on the readout rather than at node level.
  - Deep memory (ages ≥ 8) closes at k = 64, causally.
  - Ages ≥ 4 need k ≈ 256.
  - Robustness on MLP 1 is in `results/s1_w1024_mlp1.jsonl`.
  - Next: the same inside FC (exact slices + κ4 mean field) to confirm the excess at the 3e-8 level, and HOSVD truncation of the bin-[4, 8) core.
- **S2 (Bethe node sufficiency with computed cap and true free node beliefs).** In bethe v4, inject the true node beliefs minus their cap part, and supply the cap part from the law-level scale field. This tests whether (a) and (b) compose.
- **S3 (frame bound over layers, F T3).** Inject each layer's local error separately and measure the inter-layer correlation matrix. This tests the cancellation route and the Theorem-F1 reading of books-close.
- **Already decided by others.**
  - Costate's scale field at law level: queued.
  - D's causal age-multiresolution: kill if raw > 5e-8 at k = 2n/a.
  - C's trace-adaptive window: 2.43e-7 measured.

### 7.4 What the user's frame bought, honestly

- **The noncommutative side.**
  - The spectral local-to-global theorem is fully noncommutative (Theorem F1).
  - The entropic one fails at exactly one step, which is repairable by a one-site Jones-index hypothesis (Conjecture F2, consistent with BCR, GJL and LaRacuente). Near equilibrium it holds exactly (Proposition F2b).
  - The fresh-weight lemma *is* Jones' Markov-trace independence in the Brauer algebra of O(n), with index n².
- **The competition side.** The frame reorganizes all of A–E into one two-sector picture whose numbers agree across five teams. It explains why every gauge-invariant or annealed carrier is capped (they are conditional expectations onto the cap) and why local errors do not compound (Theorem F1, η = 0). It turns "carry memory cheaply" into an approximate-commuting-square defect with a measured per-layer tolerance.
- **What it did not buy.** A carrier. No measured construction has both an O(n²) state and a small free-sector defect. The live route is node sufficiency (Bethe's oracle) plus law-level cap pinning plus shared-basis free cores (S1).

### 7.5 Check of coordinator notes 2–3 against F's dictionary

**Ghost reading (note 3 §2).** The fact is correct: for a disjoint union of expanders, ⊕P_k (projections onto constants) is a non-compact ghost in the Roe algebra and obstructs coarse Baum–Connes (Higson; Higson–Lafforgue–Skandalis). Two corrections to the analogy:
- **The HLS ghost exists *because* of the spectral gap.** P = f(Δ) for a function f isolating 0, which is possible only if the gap is uniform. Our memory sits at the wall with no gap (team D F9.2; team B λ_eff ≈ 1). So "a ghost of the same kind" holds only for the ghost *property* (entries → 0, not compact), not for the *mechanism*.
- **The block ranks are not uniformly bounded.** The HLS blocks have rank 1 for every k. Our blocks have rank ≈ n/4 (ages 5–8, region §9) or 64 = n/16 (ages ≥ 8, F §5.5), growing with n. In coarse-geometric terms the memory is "blockwise finite rank" only at fixed n. The correct statement is **uniformly bounded rank per block relative to n**, a ratio, which is a weaker (L^{p,∞}-type) condition than HLS's.

**Hyperfiniteness reading (note 3 §2).** Connes 1989 (*Compact metric spaces, Fredholm modules, and hyperfiniteness*): a finitely summable Fredholm module forces a hypertrace, so reduced group C*-algebras of non-amenable groups admit none. θ-summable ones can exist. Two flags:
- At every finite n the layer algebras are finite-dimensional, hence trivially hyperfinite. The obstruction can only be **asymptotic**: in the n → ∞ limit fresh Gaussian layers generate a free semicircular family (the L(F_L) situation), and the statement must be that the summability degree p(n) of any global carrier diverges with n.
- The age axis is ℤ (amenable). The proposed escape, AF/hyperfinite along age and free within a block, is consistent with F's two-sector picture: the cap sector is the amenable, invariant part, and the free sector is the non-amenable part, frozen blockwise.

  Status: synthesis, plausible. It is not a theorem until the quantitative p(n) statement is formulated.

**Update (team G, Theorem G10; proof checked here).** The finite-n version of the hyperfiniteness reading that F asked for exists.
- *Statement.* If the fresh layers form a λ-quantum expander (Hastings; Pisier), any grading D of neuron space with ‖[D, U_i]‖ ≤ 1 has a counting function growing like (1 + (1 − λ)/4)^m per level step. So no finitely summable grading is uniformly compatible with k ≥ 2 independent layers.
- *Proof check.* F checked its three steps: the expander bound on ⟨P, ΦP⟩, the commutator bound on the far block, and the shell rank. All are correct as written for orthogonal layers.
- *Open part.* The extension to gated Gaussian steps (not unitary) is open, as G states.
- *In F's language.* The free sector is the expander direction, compressible only along the amenable age axis. The cap sector is the invariant (λ = 1) direction.
- *Agreement with S1.* S1's carrier uses bases that are *transported*, B ← orth(MᵀB), i.e. covariant under the time shift. It is not a static grading, so S1 is exactly the kind of carrier G10 allows, and the static shared subspaces that failed (region, 6×) are the kind it forbids.
- *Correction to note 2 (team G, agreed).* [F, X] has rank ≤ 2n, hence is compact. What is width-flat is its value on the state, so F's "not compact" wording above should be read as "not small in the state norm".

**Note 2 (Fredholm).** Consistent with F's dictionary:
- F = 2P − 1 with P the gauge average equals the cap/through-string split (P = E_W onto the Brauer invariants);
- [F, a] ≠ compact is the statement that the free sector has no finite-rank image under the invariant projection;
- "books close = D² over the layer filtration" is Theorem F1 at η = 0 (orthogonal martingale increments).

The proposed Dixmier lower bound Ω(n³ L log L) for covariant linear carriers is F's Conjecture F4 in the same words. **S1's result constrains it.** The log is paid only over ages ≲ 8: beyond that, a shared 64-dimensional basis is lossless at n = 1024 (MLPs 0, 1). So any Dixmier lower bound must be stated over the critical age window, not over all of L.

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
- `bulk_rank.py`: rank profile of the bulk slice.
- `cap_filter.py`: readout filters (cap, cap + rank-k bulk).
- `costate_f.py`: a copy of costate's `predict` with two added hooks, `shared` (oracle basis per dyadic bin) and `causal` (cohort bases); otherwise unchanged.
- `s1_shared.py`: runner for S1.
- `results/`: outputs.
