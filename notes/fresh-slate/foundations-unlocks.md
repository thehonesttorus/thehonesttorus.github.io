# Foundations, part 2: the theoretical unlocks

*Fresh slate, 1 Oct 2026. This note distils the programme's theoretical unlocks into a numbered list. For each unlock it gives the precise statement, its status, and its computational meaning for an estimator of quenched per-neuron means. Contract: [BRIEF.md](BRIEF.md). Part 1, [foundations-input-verified.md](foundations-input-verified.md) (the BRIEF calls it `foundations.md`), verifies the user's 1 Oct input; Appendix A lists the claims of that input used here, with their status, updated against that note's §7. An adversarial check of this note (1 Oct, after the BRIEF's 18:55 UTC correction) is logged at the end.*

---

## 0. How to read this note

**Inputs read in full.**
- The programme: [research-program.md](../research-program.md), [conditional-arrow-algebra.md](../conditional-arrow-algebra.md) (Note 1), [simplicial-complex-as-decomposition.md](../simplicial-complex-as-decomposition.md) (Note 2), [local-to-global-unlocks.md](../local-to-global-unlocks.md) and [mlp-bridge.md](../mlp-bridge.md).
- Every digest in [digests/bridges/](../digests/bridges/), and the two bridge syntheses [A-unification.md](../streams/bridges-synthesis/A-unification.md) and [B-programme.md](../streams/bridges-synthesis/B-programme.md).
- The BRIEF and the 1 Oct input.
- The six v0 design notes under [designs/](designs/), read for their exact identities and first measurements, not for their estimators.

The primary sources opened in this session are listed under Sources, with what was checked in each.

**Status labels.** Each unlock carries one label; sub-claims are labelled separately where they differ.
- THEOREM: published, with source. "(from memory)" marks a standard result not re-opened in this session.
- DERIVED: proved here (proof given or sketched in one line) or in a programme note (cited by section). "Checked" means a numerical check in §12.
- CONJECTURE: a precise statement, not proved, given with its evidence and a test.
- SPECULATION: a direction, not yet precise enough to be false.
- ANALOGY: a correspondence of form between two settings, asserted as a theorem in neither (used in unlock 40(d) and Appendix A).

**Format of an unlock.**
- **Statement**: definitions and hypotheses.
- **Status.**
- **Computational meaning**: what an estimator of quenched per-neuron means can then compute exactly, cheaply, or with a controlled error.
- **Guard**: where the unlock stops applying, when there is such a limit.

**Stream keys (the six design principles).**
- [F] faces and barycentres
- [B] Bethe/cavity on pseudorandom geometry
- [M] matchings, hafnians, signings
- [T] tropical skeleton and temperature
- [H] Heisenberg picture and Dirichlet forms
- [K] Markov networks, CMI, exchange relations

⊕ marks a cross-cutting unlock, one that serves at least three streams essentially. Its keys list every stream it serves. §8 tabulates them.

**The realisation dictionary, the only network-specific input.** The general statements are about a state on an algebra, its conditional expectations, and arrows between them. They meet the competition only through this dictionary:
- **object at layer l**: a face σ ⊆ [n], i.e. a sign pattern of the pre-activation z_l. It comes with its conditional state ω_σ (Lüders conditioning on the face); its barycentre is the first moment of ω_σ.
- **arrow σ → τ**: the diagonal projection D_σ (the ReLU on that face) followed by W_{l+1}, so z_{l+1} = z_l D_{σ_l} W_{l+1}.
- **history** μ = (σ_1, …, σ_L): a path in the layered quiver of faces. Its cone C_μ ⊆ R^n is the set of inputs that realise it.

Negative findings below are charged to this dictionary first (BRIEF rule 3), and each one says so.

**Notation.**
- x ~ γ = N(0, I_n); a_0 = x, z_l = a_{l−1} W_l, a_l = relu(z_l), g_l = 1[z_l > 0].
- The W_l have i.i.d. N(0, 2/n) entries and are independent across l.
- For a scalar random variable Z: μ_Z = E Z, s_Z² = Var Z, and p_Z is its density.
- φ, Φ and Φ̄ = 1 − Φ are the standard normal density and distribution functions.

---

## 1. The stage: exact identities of the object (all cross-cutting)

### 1. Face-history decomposition and the radial factorisation ⊕ [F, T, H, K]

**Statement.**
- For a history μ = (σ_1, …, σ_L), let C_μ = {x : g_l(x) = 1_{σ_l}, l = 1..L} and M_μ^{(l)} = W_1 D_{σ_1} W_2 D_{σ_2} ⋯ D_{σ_{l−1}} W_l.
- Each nonempty C_μ is an open polyhedral cone: on the cone of the history up to l − 1 the constraint on z_l is linear in x. The cones partition R^n up to a null set.
- On C_μ, z_l(x) = x M_μ^{(l)}. Hence, for every layer l and unit j:

  E z_l(j) = Σ_μ γ(C_μ) ⟨b_μ, M_μ^{(l)} e_j⟩,  b_μ = E[x | C_μ].
- Write x = r u, with r = ‖x‖ ~ χ_n independent of u uniform on S^{n−1}. Then γ(C_μ) b_μ = E r · E[u 1_{C_μ}(u)], and γ(C_μ) is the normalised solid angle of C_μ.
- a_l(x)/‖a_l(x)‖_1 lies in the relative interior of the face simplex Δ(σ_l(x)).

**Status.** DERIVED (piecewise linearity and positive homogeneity; the BRIEF states the homogeneity).

**Computational meaning.**
- *Exactly:* the radial law factors out as the single scalar E r = √2 Γ((n+1)/2)/Γ(n/2). Everything else is spherical geometry of cones.
- The estimand is linear in the barycentric vector measure ν(μ) = γ(C_μ) b_μ on histories. Once the history is fixed, no nonlinearity remains.
- The whole difficulty is the number of histories and the geometry of their cones. Layer 1 alone has 2^n faces, since W_1 is invertible. Every later unlock compresses this sum; none enumerates it.

### 2. ReLU is Lüders conditioning; the arrow is linear given its face ⊕ [F, H, K]

**Statement.** Let ω be the law of z_l, p_σ = 1[g_l = 1_σ] the face projections in the commutative face algebra D_0, and ω_σ = ω(p_σ · p_σ)/ω(p_σ). For every test function h:
- E h(a_l) = Σ_σ ω(p_σ) ω_σ(h ∘ D_σ);
- E[z_{l+1} 1_A] = E[a_l 1_A] W_{l+1} for every event A, in particular every face event of the layers ≤ l (z_{l+1} = a_l W_{l+1} holds pointwise).

In Note 2's frame form, where ω_σ is the barycentre vector state of the facet algebra (a coherent superposition over the vertices of σ), excluding a vertex u ∈ σ is Lüders conditioning: ω_σ^{1−e_u} = ω_{σ∖u}. This ω_σ is not the commutative one defined above: conditioned on the atom p_σ, p_σ(1 − e_u) = 0 for u ∈ σ.

**Status.** DERIVED (Note 2 §2; faces design E1).

**Computational meaning.**
- The only nonlinear operation in the network is choosing a face, i.e. conditioning. Between choices everything is a linear push-forward of conditional states.
- *Exactly:* the first moment of z_{l+1} on any event is a linear image of the corresponding barycentre of a_l.
- What is not linear is how a conditional state re-splits at the next layer's walls. Unlock 3 says how little of that is needed.

### 3. Ridge reduction and the exact mean recursion ⊕ [F, B, K, H, T]

**Statement.**
- (a) μ_{l+1} := E z_{l+1} = m_l W_{l+1}, with m_l = E a_l, exactly.
- (b) Let w_k be the k-th column of W_{l+1}. Then m_{l+1}(k) = E relu(⟨a_l, w_k⟩) depends on the law of a_l only through its n one-dimensional push-forwards along the lines w_k, equivalently through its characteristic function on n lines.
- (c) For any scalar Z with a density continuous at 0:

  E relu(Z) = μ_Z P(Z > 0) + p_Z(0) τ_Z(0).

  Here τ_Z is the Stein kernel of Z, defined by τ_Z(t) p_Z(t) = ∫_t^∞ (s − μ_Z) p_Z(s) ds; τ_Z ≡ s_Z² when Z is Gaussian. So the per-unit mean at every layer is the exact drift from (a) times the unit's own face mass, plus its own wall density times the Stein kernel at the wall.
- (d) W_{l+1} is independent of the law of a_l (BRIEF §3), so the n lines are random directions independent of the state.

**Status.** DERIVED: (a) and (b) by linearity of expectation, (c) by one integration by parts, (d) by the independence of layers. The Gaussian closure is the special case P(Z > 0) = Φ(μ/s), p_Z(0) τ_Z(0) = s φ(μ/s). Unlock 37 gives the Malliavin form of τ_Z.

**Computational meaning.**
- *Exactly:* given the exact drift μ_Z, the means need exactly two scalars per unit per layer: the face mass P(z_{l,j} > 0) and the wall term p(0) τ(0).
- With the exact drift, the error of any closure decomposes exactly as μ × (face-mass defect) + (wall defect). A drift error adds (μ − μ̂) × (the closure's face mass).
- The question the state must answer at layer l is n one-dimensional laws along fresh random directions, not the joint law. This is consistent with the EscAI oracle (BRIEF §3), but the oracle says more than (b): it measures that per-neuron cumulants to fourth order at every layer suffice, a statement made in the cumulant dictionary. (b) is exact and singles out no parametrisation of the n laws.
- *With controlled error:* for a fixed state and random lines the projections are Gaussian at leading order, under thin-shell and weak-correlation conditions on the centred state (random-projection central limit theorems: Sudakov; Diaconis–Freedman; from memory). The Gaussian closure is that leading order.
- The Gaussian closure (full covariance propagation) leaves raw ≈ 4.3e-6 at n = 1024, i.e. ≈ 2.1e-3 rms per unit (bench w1024_d16, 6 MLPs; BRIEF §1 as corrected at 18:55 UTC). Its rms error behaves as ≈ 2/n for 128 ≤ n ≤ 1024 (raw 2.9e-4, 5.4e-5, 1.8e-5, 4.3e-6 at n = 128, 256, 512, 1024; bench RESULTS.md), not as n^{−1/2}. That is consistent with the O(n^{−1/2}) quenched per-unit terms living in the variance, which the closure carries. The bar is ≈ 1.3e-4 rms (raw 1.6e-8), ≈ 16× below, so whatever supplies the closure's missing order must be right to about 6 %. (An earlier version used the BRIEF's superseded 4e-5 and said 0.2 n^{−1/2} and 2 %.)

### 4. Barycentres are boundary integrals; the mean is the Gaussian mass of the walls ⊕ [F, T, H, B]

**Statement.**
- (a) For a polyhedral cone C with facets e and outward unit normals ν_e:

  E[x 1_C] = −Σ_e ν_e γ_{n−1}(e),

  where γ_{n−1} is the Gaussian surface measure. This is the divergence theorem with ∇γ = −xγ.
- (b) Let f be positively 1-homogeneous, continuous, and piecewise linear on a fan. Then

  E f = E[x · ∇f] = ⟨γ, Δf⟩ = Σ_walls γ_{n−1}(e) [∂_ν f]_e,

  by Euler's identity and then Gaussian integration by parts. Δf is a signed measure on the codimension-1 skeleton, the tropical hypersurface. (b) is (a) summed over the histories of unlock 1: each interior facet appears twice with opposite normals, leaving the jump of the linear map.
- (c) For a ReLU output, take traces of Hess a_l(j) = δ(z_l(j)) ∇z_l(j)∇z_l(j)ᵀ + 1[z_l(j) > 0] Σ_i W_l(i,j) Hess a_{l−1}(i) and unroll:

  E a_L(k) = Σ_{l ≤ L} Σ_j E[δ(z_l(j)) ‖∇_x z_l(j)‖² ∂a_L(k)/∂a_l(j)].
- (d) Cavity reading. On the wall {z_l(j) = 0} the unit's activation is 0, so there the downstream network coincides with the cavity network in which unit (l, j) is deleted. ∂a_L(k)/∂a_l(j) is that cavity network's linear response to re-inserting the unit.
- (e) On the sphere, Δ_S F = −(n−1)F + (walls) for F = f|_S, so the total wall mass is (n − 1) times the spherical mean.

**Status.** DERIVED. Found independently by the tropical design (I2, checked at widths 8–16, depths 2–4: R-E0) and the faces design (E2), and measured neuron by neuron at n = 64, L = 16 in [foundations-problem.md](foundations-problem.md). The cavity reading (d) is also in the tropical design (§3, its transport term); the history-sum derivation of (b) is added here.

**Computational meaning.**
- *Exactly:* the mean of every neuron is a sum, over upstream units, of wall integrals. Each is (Gaussian density of the unit's pre-activation at 0) × (wall-conditional expectation of slope² × cavity response).
- Interiors of faces contribute nothing; only codimension-1 faces carry the mean, and frozen units contribute little.
- It is the exact meet-in-the-middle split: the wall measure is a forward (state) object, the cavity response a backward (question) object.

**Guard.** The identity is exact, but its wall-by-wall factorisation is not perturbative at He initialisation.
- The weight ‖∇_x z‖² is uncentred and counts the q-th Wiener chaos q times (unlock 37). Its mean exceeds s² by the mean spike and by the high-chaos roughness that depth builds: median E‖∇z‖²/s² = 1.5, 5.0, 16.7 and E‖∇z‖²/(μ² + s²) = 1.3, 2.7, 6.4 at layers 2, 8, 16 (width 64, tropical R-E2). It is not ≈ μ² + s² at depth.
- The own-wall and inherited terms each miss their factorised forms by O(1), and the two misses cancel (tropical R-E2, correlation −1.00).
- The centred form, unlock 37, removes this.

### 5. Fresh randomness: annealed twirl, quenched single Kraus operator, sign sectors ⊕ [B, M, H, K]

**Statement.** Fix the law of a_l and let W = W_{l+1} be independent of it.
- (a) *Annealed.* E_W[Wᵀ S W] = (2/n) Tr(S) I for every symmetric S. The averaged second-order map is twice the trace-preserving conditional expectation onto the scalars (the factor 2 is the He gain): a twirl, with all non-trivial eigenvalues 0. Also E[W^{⊗odd}] = 0.
- (b) *Quenched.* S ↦ Wᵀ S W is completely positive with a single Kraus operator.
- (c) For w ~ N(0, (2/n) I) and symmetric B: Var(wᵀ B w) = 8‖B‖_F²/n². This is a relative fluctuation √(2/PR(B)) around (2/n) Tr B, with PR(B) = (Tr B)²/‖B‖_F².
  - Bulk objects (PR ∝ n) fluctuate per unit at O(n^{−1/2}), and their layer averages at O(1/n).
  - A spike (PR = O(1), e.g. the mean direction) fluctuates per unit at O(1).
- (d) *Sign sectors.* The group Z_2^{n×n} of sign flips of W preserves its law.
  - Averaging a quenched quantity over it is the conditional expectation onto functions of |W|.
  - For a polynomial in W, the average keeps exactly the monomials in which every weight occurs to an even power: the paired (even) sector. It still depends on the quenched magnitudes |W|.
  - The remainder is the sign-odd, quenched sector.

**Status.**
- (a)–(c): THEOREM (Wick/Isserlis; elementary). (a) and (b) are also derived in transfer-spectrum-measurement (finding 5 and §2.5). (c) checked (C6).
- (d): DERIVED, since E_S ∏ S_e^{k_e} = 1 iff every k_e is even. It is the network analogue of Godsil–Gutman (unlock 22); unlike there, the even sector is not a tree object (unlock 22(c)).

**Computational meaning.**
- *Exactly and cheaply:* every annealed prediction is a scalar recursion, and self-averaging bulk quantities (normalised traces, layer averages, spectra) are annealed-computable to O(1/n) (second-order freeness, Mingo–Speicher, from memory).
- *Warning:* the scored quantity is per unit, and there three quenched parts enter.
  - From layer 2 on, each unit's mean carries an O(1) sign-odd term: its alignment ⟨m_{l−1}, w_k⟩ with the mean direction. The exact recursion of unlock 3(a) handles it at no cost.
  - O(n^{−1/2}) sign-odd terms pair unpaired fresh weights with the state's joint structure, e.g. Σ_{i≠j} W_ik W_jk Cov(a_i, a_j) in the unit's variance.
  - O(n^{−1/2}) sign-even terms carry the quenched magnitudes, e.g. Σ_i (W_ik² − 2/n) Var(a_i). They are not smaller: at layer 2, n = 1024, their rms is 0.060 against 0.044 for the sign-odd part, on a variance of 1.36 (C9). Replacing each unit's quenched variance by its annealed value costs 6.5e-3 rms per neuron at n = 1024 ([foundations-problem.md](foundations-problem.md)).
- The sign average removes the two sign-odd parts and keeps the sign-even one; the annealed average removes all three. Neither can be the estimator; they are references against which the quenched sectors are measured.

---

## 2. Faces and barycentres [F]

### 6. Three resolutions and the arrow algebra of the face quiver [F, K, H]

**Statement (Note 1 §2–3; Note 2 §1–3).**
- The layered quiver Λ of faces carries a Toeplitz–Cuntz–Krieger family: partial isometries s_γ, one per arrow, with s_γ* s_γ = p_{source} and Σ_{γ into τ} s_γ s_γ* ≤ p_τ.
- For a finite layered quiver, T(Λ) ≅ ⊕_σ M_{N(σ)} and C*(Λ) ≅ ⊕_{σ_0} M_{N(σ_0)}, where N(σ) counts the forward paths starting at σ, the trivial path included (Note 1 §3.1). The Bratteli heights satisfy h_{l+1} = F_l h_l, with F_l the incidence matrix.
- Three resolutions: D_0 ⊂ D ⊂ C*(Λ).
  - D_0: functions of the current face. This is the Stanley–Reisner algebra C(K) of Note 2.
  - D: functions of the history.
  - C*(Λ): histories with coherences.
- Conditioning on a path v is φ^v(a) = φ(v* a v)/φ(v* v). The Cuntz–Krieger relation is the law of total probability over the last arrow.

**Status.** THEOREM for the graph-algebra facts (as quoted in Note 1) + DERIVED (Note 1, Note 2).

**Computational meaning.**
- Every history sum of unlock 1 is a state on the path algebra, and the resolution is a design variable.
  - A state on D_0 costs one conditional state per face.
  - A state on D costs one per history, which is exponential.
  - A state on C*(Λ) costs more.
- A transfer recursion (states forward, questions backward; Exel) computes any history sum layer by layer. It is exact when the state is Markov at the chosen resolution (unlock 44). Otherwise its error is the memory (unlock 45).

### 7. KMS states are Gibbs measures on histories [F, K, T]

**Statement (Note 1 §5).**
- Take a real cocycle F on arrows, extended additively to histories, and the gauge action α_t(s_γ) = e^{itF(γ)} s_γ.
- The KMS_β states of the path algebra restrict to Gibbs measures P_β(μ) ∝ w(μ) e^{−βF(μ)} on histories. They are realised by a Doob h-transformed walk, and the modular flow is weighted depth.
- Every P_β is Markov at every layer, for every F and β. Conversely, every positive Markov chain on the quiver is KMS for some F (B-programme §4, §10; hdx digest U7).

**Status.** THEOREM (Renault, Exel, Laca–Neshveyev-type results, as quoted in Note 1) + DERIVED (Note 1; B-programme).

**Computational meaning.**
- The family {P_β} interpolates between counting histories (β = 0, the Bratteli heights) and the max-plus ground state (β → ∞; unlock 28).
- When the true history law is (close to) a member, the estimator's state is one Perron vector per layer and one transfer matrix per layer of arrows: exact and cheap.
- How far the law is from the family is measured, not assumed (unlock 45).

### 8. Sufficiency of the barycentre is a cohomology class, certified by squares [F, K]

**Statement (B-programme R1–R3, R8–R10; corrects Note 1 §6).** Let the history law be P_β ∝ w e^{−βF}.
- The pair (input face, output face) is sufficient for {P_β} iff F vanishes on bigons.
- F is a coboundary iff, in addition, the induced function Φ on the reachability graph R ⊆ K_0 × K_L is a coboundary there.
- The obstruction space is O(Λ) = Z_1(Γ)/p(B). It vanishes when R is a forest, or when each component of R has a waist face; any complete layered quiver of depth ≥ 2 has one.
- Which statistic is sufficient depends on how the family is normalised (R3).
- Quantitatively:
  - Var_β(F | ends) ≤ S²/(2 gap), with S the largest square defect (R8);
  - if H^1(X_□) = 0, then dist_2(F, B^1) ≤ ‖δ_□ F‖_2/√μ_1 (R10);
  - μ_1 = n² was measured for complete layered quivers of width n.

**Status.** DERIVED (B-programme §3, §6.4, §7). R1 is a counterexample to Note 1 §6 at depth 1. That μ_1 does not depend on depth is a CONJECTURE.

**Computational meaning.**
- When the square defects of the arrow cocycle are small, the face (barycentre) is a sufficient statistic. An estimator can then carry one number per face instead of one per history.
- The loss is bounded by a local quantity, four arrows at a time, divided by √μ_1 = n (in the ℓ² distance of R10). The Fisher-information loss is exactly E Var(F | ends).

**Guard.** The unlock requires the Gibbs form, i.e. positivity. Dictionary v1 fails it: its fitted F has residual 0.42–0.5 against coboundaries, and its face process has 25–35 % memory (mlp-bridge). The failure is charged to the dictionary.

### 9. Low-codimension face masses are elementary; high-codimension ones are not [F, H, M]

**Statement.** Let (z_1, …, z_k) be a centred Gaussian with correlations ρ_ij.
- k = 2 (Sheppard): P(z_1 > 0, z_2 > 0) = 1/4 + arcsin ρ_12/(2π).
- k = 3: P(z_1, z_2, z_3 > 0) = 1/8 + (arcsin ρ_12 + arcsin ρ_13 + arcsin ρ_23)/(4π).
- k ≥ 4: the orthant mass is a Schläfli-type function with no elementary closed form. With non-zero means, k = 2 needs Owen's T and larger k needs numerical integration (Genz).
- Plackett's reduction: ∂P(all > 0)/∂ρ_ij = φ_2(0, 0; ρ_ij) P(rest > 0 | z_i = z_j = 0), a codimension-2 face mass. So every orthant mass is an integral of lower-dimensional ones along a correlation path.
- Counting:
  - A central arrangement of N hyperplanes in general position in R^d has 2 Σ_{i<d} C(N−1, i) regions (Cover; Wendel; Schläfli). So layer 1 alone (N = d = n) has all 2^n orthants.
  - Zhang–Naitzat–Lim Thm 6.3 bounds the number of linear regions of the whole network by ∏_{l=1}^{L−1} Σ_{i ≤ d} C(n_l, i), for hidden widths n_l ≥ d and a linear output layer.

**Status.** THEOREM (Sheppard; the trivariate formula; Plackett; Cover; from memory, except that the two formulas were checked in C3 and Zhang–Naitzat–Lim was verified).

**Computational meaning.**
- *Exactly and cheaply:* masses and barycentres of faces of codimension ≤ 3 in any Gaussian layer cost O(1) work each: elementary functions when the layer is centred, as layer 1 is (exactly Gaussian), and Owen's T or a one-dimensional quadrature when the means are non-zero. Pair and triple face data of a Gaussian layer therefore cost O(n²) and O(n³) elementwise operations.
- *Never:* enumeration of faces.
- An estimator that needs the face mass of four or more units must get it by gluing (unlocks 17, 45–47) or by integrating along a path (unlock 38).

### 10. Depth windows and tolerance: how many compositions carry history [F, K]

**Statement (Note 1 §7).**
- The operator system E_k generated by depth windows of length k has C*-envelope T(Λ) and propagation number prop(E_k) = 2⌈L/k⌉ − 1 (checked there at L = 3: 5, 3, 1).
- On a single layer complex, the overlap relation |σ ∩ τ| ≥ k is a tolerance relation. It is reflexive and symmetric but not transitive, and its operator system is an algebra iff it is an equivalence.
- Overlapping faces compose with the overlap cosine as weight (Note 2's composition fork).

**Status.** DERIVED (Note 1; the tolerance facts are as quoted there from Connes–van Suijlekom).

**Computational meaning.** A window-based estimator with windows of k layers needs 2⌈L/k⌉ − 1 compositions to represent arbitrary history dependence. Its depth budget should be compared with the measured memory profile, which is experiment E5 of mlp-bridge (still open).

---

## 3. Bethe and cavity on pseudorandom geometry [B]

### 11. Trees are exact; the computation tree carries any sparse graph [B]

**Statement.**
- (a) On a tree, belief propagation (sum-product) computes all marginals exactly, in one sweep each way.
- (b) For 2-spin systems, the marginal at v in G equals the marginal at the root of the self-avoiding-walk tree T_SAW(G, v), with boundary conditions at the leaves that close cycles (Weitz).
- (c) Loopy BP has a unique fixed point and converges in either of two cases:
  - if ρ(A) < 1, where A_{i→j, k→l} = tanh|J_ij| δ_{il} 1[k ∈ ∂i ∖ j] is the weighted non-backtracking operator on directed edges (Mooij–Kappen Cor. 3). It and Dobrushin's condition do not imply each other: in Mooij–Kappen's Table I, Dobrushin holds and Cor. 3 fails in 170 of 50 000 random N = 4 models. Their local-field refinement (Cor. 4) was the strongest in their trials;
  - if the Gibbs measure on the computation tree is unique (Tatikonda–Jordan).

**Status.** THEOREM: (a) classical; (b) Weitz 2006, from memory; (c) Mooij–Kappen verified, with Tatikonda–Jordan as cited there.

**Computational meaning.** On a sparse pseudorandom geometry (large girth, bounded degree), node and edge beliefs computed by the cavity recursion are exact up to the decay of correlations along the computation tree. The certificate is one spectral radius of an O(#edges) operator.

**Guard.** Both certificates use |J| and are ℓ¹ in nature. On a dense layer with couplings of size n^{−1/2} and random signs, ρ(A) ≈ Σ_k tanh|J_ik| ~ n^{1/2}, so both fail. Unlock 12 is the dense replacement.

### 12. Dense pseudorandom couplings: Bethe becomes TAP/AMP, certified in ℓ² [B, H]

**Statement.**
- (a) For a symmetric coupling matrix with i.i.d. entries of variance σ²/n, ‖J‖_op → 2σ (the semicircle edge), while every row ℓ¹ norm is ≈ σ√(2n/π). Spectral (ℓ²) conditions are therefore dimension-free, and Dobrushin (ℓ¹) conditions fail by a factor √n.
- (b) EKZ: an Ising measure with interaction J satisfies (1 − ‖J‖) Var_μ f ≤ E(f, f) for the Glauber Dirichlet form, with the normalisation of the bridge digests.
- (c) For dense weak couplings the cavity equations acquire the Onsager reaction term (TAP). At high temperature they are the fixed point of approximate message passing, whose iterates obey a scalar state evolution (Bolthausen).
- (d) Spectral independence is an ℓ² notion: it is λ_max of the influence matrix (unlock 17).

**Status.** THEOREM: (a) Füredi–Komlós / Bai–Yin, from memory; (b) EKZ as quoted in the bridge digests; (c) from memory. The √n comparison is DERIVED from (a).

**Computational meaning.**
- Across the width of a layer, the right gluing rule is tree-shaped but carries the Onsager correction.
- Its validity is controlled by a signed spectral quantity, never by absolute row sums.
- So the measured spectral independence η ≈ 2–8 of the face law is a quantity of the right (ℓ²) kind, bounded in n so far. At that size it certifies only a bounded variance factor, not convergence or precision (unlock 17). A Dobrushin-type influence sum over the n/2 hot units is O(√n).

### 13. Power counting for a fresh layer: loops closed among the units a fresh layer touches are suppressed; loops inside the state's history are not ⊕ [B, M, K]

**Statement.** Fix the law ω of a_l, and let w = w_k be a fresh column. Organise the deviation of the law of ⟨a_l, w⟩ from its Gaussian approximation by two things: which units of layer l the fresh weights touch, and the dependences of ω that join those units.
- (a) Sign-paired (even) contributions are of the order their counting gives. Each sign-unpaired weight turns a sum over units into a random-sign sum, which costs n^{−1/2} unless the summand is coherent in sign. The mean direction is the one coherent case: one unpaired weight, O(1).
- (b) Contributions whose dependence graph through ω contains a cycle are suppressed by an extra n^{−1/2} per independent cycle, relative to tree-shaped ones.
- (c) Hence, at the needed precision, only tree-shaped (star, path) contractions of ω's dependences enter the per-unit laws at layer l + 1. But tree-shaped contractions of ω's pairwise dependences already contribute at the same order as the single-site terms. Node beliefs alone are wrong at the leading quenched order; edge beliefs are necessary.

**Status.**
- DERIVED. This is Wick counting for the fresh layer, done independently in the signings design (§1, table) and the bethe design (§2, table). It is consistent with the measured n^{−0.8} residual of first-order gate diagrams (BRIEF §3).
- CONJECTURE: (b) holds uniformly in depth for the states the network actually produces. The signings design's test T1 measured it at layer 2 only, one step out of an exact Gaussian copula at n = 12–96: the tree-only κ3 reaches the Monte Carlo floor by n = 96, and the κ4 loop residual falls ≈ n^{−2} from 48 to 96. Depth is untested.

**Computational meaning.**
- *Cheaply:* the readout of layer l + 1 from a state at layer l needs only tree-shaped contractions. These are matrix products of the fresh layer with node and edge data of the state; no loop sums are needed.
- Beyond the mean recursion, quenched information enters through node data weighted by the quenched magnitudes |W| (unlock 5) and through the state's joint structure: its edges, and how they were built through depth.
- That is where loops matter and where unlocks 43–47 apply.

**Guard.** The suppression factor of (b) is the typical pairwise correlation of the state, n^{−1/2} only up to a factor that grows with depth: the rms pre-activation correlation reaches 0.43 at layer 15 at n = 128 (bethe design §8). For a covariance of participation ratio PR the rms correlation is ≈ PR^{−1/2}; with the propagator law PR ≈ n/(2·age) of unlock 18 this would be ≈ 0.15 at n = 1024, layer 15 (an estimate, not a measurement).

### 14. Bethe is the limit of covers; for the network the lift limit is exactly computable ⊕ [B, M, F]

**Statement.**
- (a) Bethe approximations of the permanent:
  - For nonnegative θ, perm_B(θ) = lim sup_{M→∞} perm_{B,M}(θ), where perm_{B,M} is the M-th root of the average permanent over M-covers (Vontobel Thm 39).
  - per_B(A) ≤ per(A) ≤ 2^{n/2} per_B(A) (Gurvits; Anari–Rezaei Thm 4), tight for I_{n/2} ⊗ J_2.
  - On large-girth d-regular bipartite graphs, the per-vertex growth of perfect matchings is Schrijver's ½ ln((d−1)^{d−1}/d^{d−2}) (ACFK Thm 1.5), which is the Bethe value.
- (b) Network form. Let G^{(M)} be a random M-lift of the layered network:
  - each unit becomes M copies;
  - each weighted edge (i → j) becomes a uniformly random perfect matching between the copies, carrying the same weight W_ij;
  - each input coordinate becomes M i.i.d. copies.

  As M → ∞, the depth-L computation tree of a unit has distinct leaves with probability tending to 1. In the limit the parents of every unit are independent, and its pre-activation law is the convolution of its parents' scaled laws: φ_{z_{l,k}}(t) = ∏_i φ_{a_{l−1,i}}(W_{ik} t).

**Status.** (a) THEOREM (Vontobel, Gurvits, Anari–Rezaei, ACFK; all verified). (b) DERIVED, sketch: random lifts are locally tree-like. The rate in M is not established.

**Computational meaning.**
- *Exactly:* the tree (Bethe) reference of the quenched network is a well-defined object. It is computable by n² one-dimensional convolutions per layer (characteristic functions on a grid). It keeps every quenched (signed) weight, and its single-site laws are exactly computable and non-Gaussian. They are the lift's laws, not the network's: from layer 2 on they omit every covariance between parents, which changes each unit's variance at the leading quenched order (unlock 5) and more at depth.
- Measured: node beliefs alone (the lift limit) give raw 2.5e-3 and 3.3e-3 at n = 64 and 128, not improving with width (bethe design §8). The cover limit's loop corrections are the covariance between parents and the old content (signings design §8, item 4).
- The difference between the quenched network and its lift is exactly the loop content.
- (a) also says that for permanent-type sums the Bethe value is accurate only to e^{O(n)}, i.e. to O(1) per site. The Bethe reference is a per-site object, never a precision tool on its own, consistent with unlock 13(c).

### 15. Locality comes from zero-freeness; hard constraints are non-local unless the geometry expands [B, M, T]

**Statement.**
- (a) Heilmann–Lieb: the matching polynomial of a graph of maximum degree d has only real roots, all of modulus ≤ 2√(d−1).
  - The power sums of its roots count closed tree-like walks, i.e. walks in Godsil's path tree (ACFK Remark 3.6).
  - On large-girth sequences the root distribution tends to Kesten–McKay (ACFK Thm 4.1).
- (b) The matching entropy per vertex is estimable: it converges along every Benjamini–Schramm convergent sparse sequence (ACFK Thm 1.2, activity 1). The monomer–dimer free energy at any positive activity λ, ½∫ ln(1 + λx²) dρ_G, follows in the same way from the weak convergence of the matching measures (ACFK Thm 3.5).
  - The number of perfect matchings is not estimable, even for d-regular bipartite graphs (Thm 1.8).
  - In some such graphs one edge lies in all but a c^n fraction of the perfect matchings (Thm 1.7).
- (c) On bipartite δ-expanders every edge has p(e) ≥ (1/d) n^{−2 ln(d−1)/ln(1+δ)} (Thm 1.9). With Gamarnik–Katz, the growth of perfect matchings becomes local.
- (d) General principle (Barvinok; Patel–Regts; from memory). Suppose a partition function has no zeros in a neighbourhood of the path from 0 to the target activity. Then its logarithm is determined to ε by O(log(N/ε)) Taylor coefficients, each a sum over connected local structures.

**Status.** THEOREM: (a)–(c) verified in ACFK; (d) from memory.

**Computational meaning.**
- An estimator may trust a local (bounded-radius) computation when the relevant generating function is zero-free along its interpolation path (a sufficient condition, not a characterisation). The order (equivalently, the radius of the local structures) needed is O(log(N/ε)/log(R/|λ|)), where |λ| is the target activity and R the distance to the nearest zero, after a conformal map of the zero-free region to a disk.
- Hard (zero-temperature) constraints break locality unless the geometry expands.
- *Network reading:* the linear-to-ReLU homotopy relu_γ has its singularities at γ = (1 ± i)/2, which gives the geometric rate 0.71 per order at γ = ½ (tropical design P6). A conformal re-parametrisation of the path that avoids these points is the Barvinok-type repair. It was not tried (SPECULATION).

### 16. Uniqueness, reconstruction, and why expansion hurts separators [B, K]

**Statement.**
- On the d-regular tree (branching number d − 1), Gibbs uniqueness for the zero-field ferromagnetic Ising model holds iff (d−1) tanh β ≤ 1, an ℓ¹ condition.
- Reconstruction (Kesten–Stigum) occurs when (d−1) tanh² β > 1, an ℓ² condition on the second eigenvalue.
- Mossel–Sly: for the ferromagnetic Ising model with arbitrary external fields on any graph of maximum degree d, Glauber dynamics mixes in O(n log n) when (d−1) tanh β < 1.
- Yang's CMI bound has a boundary gain that vanishes on expanders (Yang digest, B1). Separator bounds are therefore weakest exactly where tree bounds are strongest.

**Status.** THEOREM (as quoted in the expanders and hdx digests; Kesten–Stigum and Mossel–Sly also from memory).

**Computational meaning.**
- Between the ℓ² and ℓ¹ thresholds, local tree computations still work on average, through spectral criteria, although worst-case influence sums do not. Dense pseudorandom layers live in this window (unlock 12).
- At low temperature (frozen faces with strong couplings), expansion creates rigidity (unlock 15(b)), and neither tree gluing nor separator gluing is local.

### 17. Spectral independence of the face law: what it gives and what it does not [B, K, F]

**Statement.** Let ω be a law on {0,1}^V (the face law, restricted to the non-deterministic vertices), with D = diag Var(e_u).
- The influence matrix is Ψ(u, v) = P(v | u) − P(v | ū) = Cov_uv/Var_u for u ≠ v, with Ψ(u, u) = 0, i.e. Ψ = D^{−1} Cov − I. Then λ_max(Ψ) ≤ η_0 iff Cov ⪯ (1 + η_0) D.
- Hence every linear observable Σ_u a_u e_u has variance at most (1 + η_0) times its product-state value (B-programme R13; markov design P4, citing Chen–Eldan Remark 37).
- The local-to-global theorems need spectral independence of all pinnings (U6). These are:
  - ALO: λ_2 of the down-up walk from the link spectra;
  - CLV: entropy factorisation and MLSI;
  - trickle-down;
  - Alev–Lau.

  They hold for arbitrary laws, Markov or not (U9).
- Measured:
  - unpinned η_0 ≈ 1.7–5.6 at n = 128 and 2–8 at n = 1024;
  - 472 of 1024 gates are uncertain (0.02 < p < 0.98) at depth;
  - at n = 128, a zero-mean Gaussian surrogate with the same correlations gives η_0 = 15–29, against 1.7–5.6 measured. Part of the gap is only the gate count (the surrogate keeps all 128 gates, the measurement the 53–78 uncertain ones); per gate the measured η_0/m is 2–3× below the surrogate's.

**Status.** DERIVED (R13) + THEOREM (ALO, CLV, as read in the hdx digest) + measured (transfer-spectrum §3.5). That η_0 stays bounded in n is a CONJECTURE: it levels off between n = 512 and 1024 at layers 9–12, and is still rising at layers 13–15.

**Computational meaning.**
- *With a controlled error:* any linear statistic of the gate field has variance within a factor 1 + η_0 ≈ 3–9 of its independent-gate value. This includes the gate-field part of every next-layer pre-activation.
- That certifies near-product structure up to a bounded factor, not to the bar's precision.
- Mixing certificates would need pinned η, and are numerically vacuous at η ≈ 7 (R13).
- Use it to bound corrections and to choose what is carried, not as an estimator.

### 18. No gap in the signed bulk: the free-probability law [B, H, K]

**Statement.**
- Products of gated Gaussian layers J_{s→t} = W_{s+1} Φ_{s+1} ⋯ W_t obey free multiplicative convolution at the competition's shape (d/n = 1/64, the free regime).
- Participation ratios compose as r(a ⊠ b) − 1 = (r(a) − 1) + (r(b) − 1), with PR = n/r. Measured: PR ≈ n/(2·age).
- The number of modes k_ε needed for a fixed accuracy of the old-content transport has k_ε/n tending to a constant as n grows.
- The only width-independent structure is one outlier along the mean direction μ_z: the BBP mechanism, iterated. Its overlap ⟨u_1, μ_z⟩² is 0.85 at n = 1024, age 15.

**Status.** THEOREM (free multiplicative convolution; rank-one multiplicative spikes, Benaych-Georges–Nadakuditi; as read in the transfer-spectrum digest, which labels them KNOWN-LINK) + DERIVED (the PR composition law) + measured (transfer-spectrum-measurement §2.2–2.6, §3.4, B2, B6; Hanin–Nica for the polymer reading). PR ≈ n/(2·age) holds to 4–6 % for ages ≤ 4–6; PR·age/n falls to 0.39 at age 15 for n = 1024. That k_ε/n is constant is a CONJECTURE beyond the generic sources tested.

**Computational meaning.**
- *Negative and decisive, in the dictionary measured:* no fixed-rank compression of old third-order content (the (2,1) slice) transported along the gated propagators works, in propagator, transported-covariance or HOSVD bases. The rank needed is a constant fraction of n, ≈ 0.3 n at the brief's accuracy.
- The one positive (Perron) mode, the mean direction, can be carried separately and exactly at rank one.

### 19. Orientation scrambling: compression error is predicted by the propagator alone [B, H]

**Statement.**
- J_{s→t} = R · W_{s+1}, where W_{s+1} is Gaussian, independent of the source created at layer s, and right-orthogonally invariant.
- Up to the dependence of the downstream gates on W_{s+1}, the orientation of the source relative to J's right singular vectors is therefore Haar-random.
- The expected error of any projection built from J depends on the source only through O(n)-invariant quadratic data.

**Status.** DERIVED (transfer-spectrum §2.4, B5). Verified for k_{2%} at 11–12 of 15 ages. The per-k curves deviate by up to 40 % at ages 6–14 and by up to 1.9× at age 15.

**Computational meaning.**
- *Cheaply:* the error of truncating any carried object to k directions is predicted before it is computed, from the singular values of the gated propagator.
- Carrier design reduces to spectral data of J. Combined with unlock 18, it says how much must be carried.

### 20. Collective coordinates and near-product conditionals [B, K, T]

**Statement.**
- (a) Hubbard–Stratonovich. Let P(s) ∝ exp(h·s + ½ sᵀ V Vᵀ s) on {0,1}^n, with V ∈ R^{n×r}. Then P = ∫ π(g) ⊗_i Bern(sigmoid(h_i + (Vg)_i)) dg exactly, with π(g) ∝ φ_r(g) ∏_i (1 + e^{h_i + (Vg)_i}): a mixture, over an r-dimensional collective variable g, of product measures.
  - If J = V Vᵀ + J', the mixture components have interaction J'.
  - The EKZ condition applies to them when ‖J'‖ < 1.
- (b) Network:
  - At depth one outlier direction (the mean; unlock 18) carries 0.88 of the second moment (ρ_16 at width 256).
  - Correlations among the uncertain gates fall as the mean shift freezes units: at n = 128, η_0 = 1.7–5.6 measured against 15–29 for the zero-mean surrogate, partly a gate-count effect (unlock 17); mean |Cor_ij| among the uncertain gates is still 0.08–0.10 at depth.

**Status.** (a) THEOREM (Hubbard–Stratonovich, exact; the general localisation scheme of Chen–Eldan is read in the hdx digest). (b) measured.

**CONJECTURE.** Conditioning the face law of a deep layer on its top one or two principal coordinates reduces its unpinned spectral independence to O(1) with a small constant.
- Gates are 0-homogeneous, so these coordinates are angular.
- Test: η_0 of the conditional gate law, given the top principal coordinate of the pre-activation field.

**Computational meaning.** If the conjecture holds:
- Per-unit means become an integral, over one or two collective coordinates, of near-product conditional laws. The integral is a quadrature with 10–30 nodes. Each conditional law is computable by tree gluing, with a certified variance factor close to 1.
- That is cheap with a controlled error. It is the one place where a low-dimensional structure is real: the rank-one Perron mode.

---

## 4. Matchings, hafnians, signings [M]

### 21. Gaussian expectations are matching sums; quasi-free states [M, H, B]

**Statement.**
- (a) Isserlis–Wick. For a centred Gaussian vector with covariance R, E[g_{i_1} ⋯ g_{i_{2m}}] = haf(R[i_1..i_{2m}]), the sum over perfect matchings of the index multiset. For a complex Gaussian with covariance A, E ∏_i |g_i|² = per(A).
- (b) Diagram formula. Let F̂_v(d) = E[F_v(g) He_d(g)] be the Hermite profiles. Then E ∏_v F_v(g_{i_v}) is a sum over loopless multigraphs on the vertices, weighted by ∏_v F̂_v(deg v) and ∏_edges ρ^{m_e}/m_e! (a loopless hafnian of the half-edge matrix).
  - Joint dependence beyond the product is the connected part.
  - The two-vertex sector resums by Mehler's kernel (signings design P2).
- (c) Fermionic (Grassmann) Gaussian integrals are Pfaffians or determinants.

**Status.** THEOREM ((a) and (c) classical; (b) the Mehler/diagram formula as stated in the signings design).

**Computational meaning.**
- *Exactly:* every joint expectation of ReLU functionals of a Gaussian or Gaussian-copula layer is a weighted matching sum.
  - Its two-vertex part is closed-form (Mehler).
  - Its quasi-free resummation is one determinant and one resolvent (unlock 23); the degree-≤ 2 part itself (vertex-disjoint paths and cycles) is an α-permanent, not a determinant.
- Bosonic sums (hafnians, permanents) are #P-hard in general (Valiant, from memory).
- An estimator must therefore use the structure of the sum: trees (unlock 13) and quasi-free resummations (unlock 23). It must never evaluate the full hafnian.

### 22. Sign-averaging is the even-sector projection; a quenched network is one signing ⊕ [M, B, T]

**Statement.**
- (a) Godsil–Gutman. For a graph G with adjacency matrix A and a uniformly random signing s, E_s det(x − A_s) = μ(G, x), the matching polynomial. The cycle terms of the determinant expansion carry an odd power of some edge sign and vanish on average; the matchings survive.
- (b) The moments of the matching polynomial's roots count closed tree-like walks, i.e. walks in Godsil's path tree, the tree of self-avoiding paths from a vertex, a finite subtree of the universal cover (ACFK Remark 3.6).
  - A signing is a 2-lift, i.e. a Z/2 gauge field.
  - Interlacing families give a signing whose largest new eigenvalue is ≤ 2√(d−1); for bipartite graphs this yields Ramanujan 2-lifts (Marcus–Spielman–Srivastava, Interlacing Families I; from memory).
- (c) Network form (unlock 5(d)). Averaging any quenched quantity over the signs of a fresh layer gives its even sector: every fresh weight occurs to an even power, and the quenched magnitudes |W| remain. In the network this is not a tree object. Three references must be kept apart:
  - the annealed average E_W, a scalar recursion with no per-unit information;
  - the sign average, which keeps |W| (e.g. Σ_i W_ik² Var a_i) and the state's pair dependences met by paired weights (e.g. Σ_{i≠j} W_ik² W_jk² κ(a_i, a_i, a_j, a_j)), and drops every odd term, including the diagonal Σ_i W_ik³ κ3(a_i);
  - the Bethe lift of unlock 14, which keeps every quenched signed weight and drops the parents' joint structure.

  From layer 2 on, the per-unit mean has two sign-odd parts:
  - an O(1) part, its alignment with the mean direction;
  - an O(n^{−1/2}) part, from unpaired weights joined by the state's joint structure.

  It also has a sign-even quenched part of the same O(n^{−1/2}) order (unlock 5, C9).

**Status.** (a), (b) THEOREM (Godsil–Gutman; Godsil; ACFK; MSS). (c) DERIVED.

**Computational meaning.**
- *Exactly:* the decomposition quenched = (sign average) + (sign-odd sectors) is canonical and computable sector by sector.
- The sign average is cheap and wrong per unit at O(1), because it drops the mean alignment. It is neither the annealed object nor the Bethe lift, and it is not free of quenched information: its even sector carries the |W| fluctuations at O(n^{−1/2}).
- A design built on an annealed or sign-averaged object must therefore carry:
  - the coherent sign-odd sector through the exact mean recursion (unlock 3(a));
  - the incoherent sign-odd sector through the state's edges (unlock 13);
  - the sign-even quenched sector through the quenched magnitudes (each unit's variance).
- A design built on the Bethe lift keeps the signed weights and must add the parents' joint structure (unlock 14).

### 23. The free (determinantal) resummation of an interacting Gaussian sum [M, H, K]

**Statement.**
- (a) If every vertex profile is the exponential of a polynomial of degree ≤ 2, F_v(g) = c_v exp(b_v g + ½ d_v (g² − 1)) (a quasi-free vertex), then E ∏_v F_v(g_v) is a Gaussian integral of the exponential of a quadratic form, summed exactly by det^{−1/2}(I − D R) and a resolvent. Here −½ log det(I − D R) = ½ Σ_k tr((DR)^k)/k sums closed walks that may revisit vertices (free bosons). With Grassmann variables the walks become vertex-disjoint cycles with a sign per cycle, det^{+1} (free fermions).
- (b) The sub-sum of the multigraph expansion of unlock 21(b) over multigraphs of maximum degree ≤ 2 (vertex-disjoint paths and cycles, each vertex used once) is not this determinant. Its 2-regular part is the α-permanent of the zero-diagonal correlation matrix at α = ½, a hard-core loop gas. Example: two vertices with profiles F̂ = (0, 0, 1) give the polynomial ρ²/2, which vanishes at ρ = 0 and nowhere else, while a constant times a Gaussian integral of the exponential of a real quadratic form is either identically zero or never zero (C10).
- (c) After renormalising the vertex weights, the two agree on every diagram that visits each vertex at most once. They differ by walks that revisit a vertex (in the two-vertex example, the multiple edges of multiplicity ≥ 4), which the determinant counts with quasi-free weights.

**Status.** (a) THEOREM (Gaussian and Grassmann integrals). (b)–(c) DERIVED here and checked (C10). They correct the signings design's P3(i), which equated the degree-≤ 2 sector with the determinant. Its estimator never used the determinant (signings design §8), so its results stand.

**Computational meaning.**
- *Exactly and cheaply:* the quasi-free resummation of any sum over pair interactions is one determinant and one resolvent, O(n³). It is not the paths-and-cycles sector: it adds vertex-revisiting walks with quasi-free weights, so as an approximation to that sector it is controlled only to leading order after renormalising the vertex weights.
- Beyond it, vertices of degree ≥ 3 (the interacting hubs) and the hard-core corrections both need separate treatment.
- This is the matching-side form of the separator rule of unlock 47: quasi-free parts are glued by linear algebra, interacting parts by trees. By signings T1, loops among fresh-weight legs are negligible at n ≥ 96, so the network's one-step sums never needed this resummation (unlock 13).

### 24. Determinantal evaluation needs Pfaffian orientations; dense layers have none [M, T]

**Statement.**
- Kasteleyn: a planar graph has an orientation that turns its perfect-matching count into a Pfaffian. Genus g needs 4^g Pfaffians (Galluccio–Loebl; Tesler; Cimasoni–Reshetikhin, verified in foundations-input-verified.md).
- Little: a bipartite graph is Pfaffian iff it contains no even subdivision of K_{3,3} as a central subgraph.
- Robertson–Seymour–Thomas (and McCuaig): a brace is Pfaffian iff it is the Heawood graph or is built from planar braces by repeated 4-sums. This is recognisable in O(n³).
- Lieb–Loss Thm 3.1: on every planar bipartite graph, with any hopping amplitudes, the canonical flux (π through every square face, 0 through every hexagon; in general π for faces of length 4k and 0 for 4k + 2) maximises |det H| for the hopping matrix H carrying those phases, and gives det H = ±D², with D the dimer partition function. This is Kasteleyn's face rule, and with Lieb–Loss's Theorem A.1 it is another proof of Kasteleyn's theorem.
- Lieb (1994): the same flux minimises the half-filled-band energy on the lattices covered by his theorem. For arbitrary planar bipartite graphs and amplitudes it does not (Lieb–Loss §VII(B): a four-box graph whose limit is a ring of eight sites).

**Status.** THEOREM (Little, RST with its O(n³) recognition algorithm, Lieb–Loss Thm 3.1 and §VII(B), and Lieb verified; the genus count, 2^{2g} = 4^g Pfaffians indexed by spin structures, is verified in foundations-input-verified.md via Cimasoni–Reshetikhin).

**Computational meaning (guard).**
- K_{n,n} is not Pfaffian for n ≥ 3. No signing of a dense layer-to-layer graph turns its matching sums into determinants, so exact determinantal evaluation of the network's matching-type sums is unavailable.
- Determinants enter only through quasi-free sectors (unlock 23) or randomised signings (unlock 25).
- The flux picture (signs as a ground-state gauge field) is SPECULATION for the network. Its weight signs are physical, not a gauge, because the ReLU is not odd.

### 25. Randomised determinants: the Godsil–Gutman estimator and algebra-valued signs [M]

**Statement.**
- For a nonnegative n × n matrix A, let B_ij = ε_ij √A_ij with i.i.d. uniform signs. Then E det(B)² = per(A).
- The critical ratio (second moment over squared mean) depends on the sign algebra:
  - ≤ 3^{n/2} for real signs;
  - ≤ 2^{n/2} for complex signs (Karmarkar et al.);
  - ≤ (3/2)^{n/2} for quaternions;
  - (1 + O(2^{−k/2}))^{n/2} for signs in the Clifford algebra with k generators, hence bounded once 2^k ≳ n² (Chien–Rasmussen–Sinclair). No polynomial-time algorithm is known for these noncommutative determinants beyond the quaternion case k = 3 (Moore–Russell's summary, arXiv:0906.1702).
- On random instances the ratio is small: polynomial for random 0/1 matrices (Frieze–Jerrum), and subexponential w.h.p. (Costello–Vu).

**Status.** THEOREM (the critical ratios of Karmarkar et al. and Chien–Rasmussen–Sinclair, re-checked through Moore–Russell's summary; the random-instance results from memory). Unbiasedness checked for n = 3, 4 by exact enumeration (C4).

**Computational meaning.**
- Bosonic quantities (permanents, hafnians) of nonnegative weight matrices have unbiased determinantal estimators. Their variance is controlled by the algebra of the signs and is small on random instances.
- For the network this is SPECULATION. Signed lifts of the quenched layers could serve as control variates: evaluate a cheap estimator on the network and on its signed lifts, and use the difference to isolate loop content (signings design §7).
- No identity yet makes the network's mean a permanent.

### 26. Real stability, interlacing and paving [M, B]

**Statement.**
- (a) Mixed characteristic polynomials are real-rooted and form interlacing families (MSS). Consequences:
  - Kadison–Singer: pure states of the diagonal subalgebra of B(ℓ²) extend uniquely.
  - Weaver's KS_2 form (MSS II Cor 1.5): if Σ_i v_i v_iᵀ = I and ‖v_i‖² ≤ δ, there is a partition into r parts with ‖Σ_{i∈S_j} v_i v_iᵀ‖ ≤ (1/√r + √δ)².
  - Anderson paving (Thm 6.1): every zero-diagonal self-adjoint T is (r, ε)-pavable with r = (6/ε)^4. That is, there are coordinate projections P_1..P_r summing to I with ‖P_i T P_i‖ ≤ ε‖T‖. r ≥ 1/ε² is necessary.
- (b) Gurvits: for real-stable generating polynomials, capacity bounds give per ≥ n!/n^n for doubly stochastic matrices (van der Waerden); capacity is computed by Sinkhorn scaling. The Bethe lower bound per ≥ per_B is also Gurvits's, but his proof goes through Schrijver's permanental inequality; proofs by real stability (Anari–Oveis Gharan) and by 2-lifts (Csikvári) came later (Anari–Rezaei §1.2).
- (c) Real stability of a multi-affine generating polynomial implies negative dependence: it is strongly Rayleigh (Borcea–Brändén–Liggett, from memory).

**Status.** THEOREM ((a) MSS II verified; (b) through Anari–Rezaei; (c) from memory).

**Computational meaning.**
- (a) is SPECULATION for the network. Pave the zero-diagonal matrix T = D^{−1/2} Cov(g) D^{−1/2} − I of the face law. Its top eigenvalue is η_0, and ‖T‖ ≤ max(η_0, 1).
  - Paving splits the units of a layer into r = (6‖T‖/ε)^4 blocks, inside each of which the face law has spectral independence ≤ ε. Random partitions do comparably well on random-like T.
  - The face law would then be near-product inside blocks and coupled only across blocks: a certified block mean field.
  - At n = 1024 the theorem's block count is vacuous: r = (6‖T‖/ε)^4 ≥ 6^4 = 1296 > n even at ε = ‖T‖. Only the random-partition heuristic is usable at this width.

**Guard for (b) and (c).** The network's face law has positive correlations, so it is not strongly Rayleigh. Real-stability tools (capacity, SLC exchange walks, the U2 conjecture of local-to-global-unlocks) do not apply to it directly. They apply to matching-type objects, such as Godsil–Gutman's sign average of a determinant; the network's sign averages are not of that type (unlock 22(c)).

---

## 5. Tropical skeleton and temperature [T]

### 27. ReLU networks are tropical rational maps [T, F]

**Statement (Zhang–Naitzat–Lim).**
- With integer weights (rational weights after scaling), every bias-free ReLU network is a tropical rational map: a difference of two tropical polynomials (Thm 5.2). With real weights it is a tropical rational signomial map (Prop 5.6).
- Every neuron is z = P − Q, with P and Q convex, positively homogeneous and piecewise linear: support functions of polytopes (Newton polytopes).
  - A linear layer acts by Minkowski sums, with W = W⁺ − W⁻.
  - A ReLU acts by a convex hull: relu(P − Q) = max(P, Q) − Q.
- Linear regions are the cones of the common refinement of the normal fans. They correspond to vertices of upper faces.
- The number of regions is at most ∏_{l=1}^{L−1} Σ_{i ≤ d} C(n_l, i), for hidden widths n_l ≥ d and a linear output layer (Thm 6.3).
- ReLU networks with biases are exactly the continuous piecewise-linear functions (Arora et al.). Bias-free ones are exactly the positively homogeneous ones (DERIVED: a max–min lattice representation with linear pieces, and max(a, b) = b + relu(a − b)).

**Status.** THEOREM (verified in Zhang–Naitzat–Lim; the P − Q recursion as in the tropical design §1).

**Computational meaning.**
- The faces of unlock 1 are the cells of a tropical fan, and the walls of unlock 4 are its tropical hypersurface.
- At zero temperature, the per-neuron mean is a difference of Gaussian mean widths, E a = w_G(conv(A ∪ B)) − w_G(B). Unlock 29 shows that this difference cancels catastrophically.

### 28. Temperature: Maslov dequantisation of the ReLU [T, H]

**Statement.**
- max(u, v) = lim_{T→0} T log(e^{u/T} + e^{v/T}), so relu = lim_{T→0} softplus_T.
- Each tropical polynomial becomes a subtraction-free exponential sum Z_T(x) = Σ_v e^{⟨v, x⟩/T} over the vertices of its Newton polytope. Minkowski sums become products and convex hulls become sums.
- For a bounded density p that is twice differentiable at 0: E softplus_T(Z) = E relu(Z) + (π²/6) T² p(0) + O(T⁴) (checked, C11). If p is only continuous at 0, the remainder is o(T²).
- For the network, F(T) = F(0) + c_2 T² + …, where c_2 is a density-weighted Gaussian mass of the tropical hypersurface (tropical design I4).

**Status.** THEOREM (Maslov dequantisation, standard) + DERIVED (tropical design I4).

**Computational meaning.**
- The temperature is a deformation parameter whose first correction is again a wall mass. That makes it a diagnostic: F(T) − F(0) ∝ T² × (hypersurface mass).
- It is not a cheaper route: nothing is cheaper at finite T, and the high-T series is only asymptotic (tropical design P4).

### 29. Signs do not die at zero temperature; they are the whole answer ⊕ [T, M, H]

**Statement.**
- In the P − Q form, the Gaussian mean widths of the Newton polytopes are exactly additive under Minkowski sums and grow by ≈ E Σ_i |W_ij| = √(4n/π) per layer.
- At n = 1024, their median ratio to the answer is 18 at layer 2, 2.4·10⁴ at layer 4 and 1.2·10²³ at layer 16 (tropical design R-P5).
- The signed path sum E a_L = Σ_paths ∏ W · E[x_{i_0} ∏ gates] is a directed polymer in weak disorder. The per-step energy ln|W| has standard deviation ≈ 1.1, against the single-path (Viterbi) threshold √(2 ln n) ≈ 3.7 (tropical design P2). No path dominates, and the sum self-averages.

**Status.** DERIVED + measured (tropical design R-P5, P2; the polymer identity is Hanin–Nica's, as quoted in the transfer-spectrum digest).

**Computational meaning (guard).**
- Any subtraction-free (positive, tropical) representation of the quenched network must resolve an exponential cancellation: 10²³ at the competition shape.
- The max-plus skeleton (dominant cone, dominant path, assignment problem) carries no usable information.
- This is charged to the dictionary "tropical skeleton = dominant cone or path", not to the theory. The theory's tropical content survives in the exact wall identities 4 and 37, which are signed.

### 30. Gate temperature: frozen units are exactly linear arrows ⊕ [T, F, H]

**Statement.**
- For a unit whose pre-activation is N(m, s²), let t = |m|/s. Then E relu(Z) − relu(m) = s[φ(t) − tΦ̄(t)] ≤ s φ(t)/(1 + t²) (by the Mills-ratio bound Φ̄(t) ≥ φ(t) t/(1 + t²)).
- On a sub-frame where every gate is constant with probability 1 − ε, the arrow D_σ W is a fixed linear map up to an event of mass ε, and no new joint structure is created there. All splitting of conditional states happens at hot units, weighted by their wall densities p(0) ≈ φ(t)/s (Price, unlock 38).
- Measured at width 256, depth 16:
  - the hot fraction (|μ/s| < 1) is 1.00, 0.86, 0.58, 0.30, 0.26, 0.24 at layers 1, 2, 4, 8, 12, 16;
  - the frozen fraction (|μ/s| > 3) is 0.27 at depth 16.
- At n = 1024, 472 of 1024 gates are uncertain (0.02 < p < 0.98) at depth.

**Status.** DERIVED (checked, C7) + measured (tropical design R-P1; transfer-spectrum §3.5).

**Computational meaning.**
- *Exactly, at the bar:* a unit can be replaced by its linear arrow only if t ≳ 4, where the per-unit error is ≈ 7e-6·s (7.1e-6·s at t = 4).
  - At t = 3 the error is 3.8e-4·s.
  - At t = 2.054 (gate certainty 0.98) it is 7.3e-3·s, sixty times the bar's rms when s = 1.
- *Cheaply:* the non-linear (face-splitting) work of a layer is confined to its hot sub-frame, about half the layer at depth for n = 1024.
- The tropical skeleton is a cost saving on the frozen half, not an approximation scheme (unlock 29).

### 31. Depth cools, slowly and incompletely [T, B]

**Statement.**
- For the annealed (infinite-width) network, the correlation of activations at two independent inputs evolves by the arc-cosine map c ↦ (√(1−c²) + (π − arccos c) c)/π, starting from c_0 = 0.
- This c_l is also the mean's share of the second moment, ρ_l = ‖E a_l‖²/E‖a_l‖².
- c_l = 0.318, 0.494, 0.681, 0.834, 0.897, 0.929 at l = 1, 2, 4, 8, 12, 16.
- 1 − c_l ~ 9π²/(2l²) only asymptotically: l²(1 − c_l) = 43.7 at l = 2000, against 9π²/2 = 44.4, while 1 − c_16 = 0.071.
- Measured at width 256: ρ_16 = 0.88. Freezing saturates near 25 % hot.

**Status.** DERIVED (Cho–Saul arc-cosine kernel; the iteration checked, C5) + measured (tropical design R-P1).

**Computational meaning.**
- The tropical (zero-temperature) limit is approached polynomially in depth and is not reached at L = 16.
- The mean direction carries most of the second moment. It must be propagated exactly (unlock 3(a)), and every centred quantity must be defined relative to it (unlock 37).
- The rate law tells a design how fast the savings from frozen units grow with depth.

### 32. Crossing from zero to finite temperature needs a structure [T, M, K]

**Statement.**
- The formal tropicalisations of the permanent and the determinant coincide: both become the assignment problem, which is solvable in polynomial time. The zero-temperature limit of T log|det| can fall below the assignment value, or to −∞, when top-weight permutations of opposite sign cancel (foundations-input-verified.md, C15). Signs survive the limit here too (unlock 29).
- The finite-temperature permanent is #P-hard (Valiant, from memory).
- Each known polynomial crossing rests on a structure:
  - *positivity:* the FPRAS for nonnegative permanents by Markov chains (Jerrum–Sinclair–Vigoda, from memory), and the Birkhoff–Hopf contraction of positive transfer operators (unlock 40);
  - *planarity or Pfaffianity:* Kasteleyn (unlock 24);
  - *real stability:* Gurvits (unlock 26);
  - *cluster periodicity:* Chin Thm 1.1. For every Zamolodchikov periodic B-matrix, T_i(t + h_Γ + h_Δ) = T_{σ(i)}(t), with σ an automorphism of order ≤ 2. The proof is tropical (maximal green sequences, tropical c-vectors, the separation formula). It lifts to the full system by the Inoue–Iyama–Keller–Kuniba–Nakanishi theorem that tropical periodicity implies full periodicity.

**Status.** THEOREM (Chin verified; the others as marked).

**Computational meaning (guard).**
- A design that starts from a zero-temperature skeleton must name the structure that carries it to finite temperature.
- For the quenched network none is present:
  - the weights are signed;
  - the layer graphs are dense and non-Pfaffian;
  - the face law is not strongly Rayleigh;
  - no periodicity is known.
- The tropical route is therefore closed, at the level of the dictionary, unless a new structure is exhibited.

### 33. Exchange relations tropicalise: CMI is a function of a y-ratio ⊕ [T, K, M]

**Statement.**
- Take a Gaussian with covariance M on a ⊔ B ⊔ c, with a and c single coordinates (blocks: last bullet). The Desnanot–Jacobi (Dodgson) relation in subtraction-free form reads

  det M[aB] det M[Bc] = det M[aBc] det M[B] + (det M[aB|Bc])².
- Put y = (det M[aB|Bc])²/(det M[aBc] det M[B]). Then I(a : c | B) = ½ log(1 + y).
- Markov (a ⊥ c | B) is y = 0: the degeneration of the exchange relation to one monomial.
- Under the tropical limit (log scale, max-plus), ½ log(1 + y) → ½ max(0, log y).
- For blocks a and c the CMI is ½ log(det M[aB] det M[Bc]/(det M[aBc] det M[B])) = −½ Σ_i log(1 − ρ_i²), over the canonical correlations ρ_i of a and c given B (unlock 46(b); foundations-input-verified.md, Theorem A). A single cross minor no longer suffices.

**Status.**
- DERIVED from Desnanot–Jacobi (checked, C1).
- Identifying y with a cluster y-variable of an octahedron-recurrence seed is the user's Idea A. As an organising structure for Markov networks it is SPECULATION here; foundations-input-verified.md §2 derives part of it (separator changes as cube moves; graphical models as coordinate strata for chordal graphs on ≤ 6 vertices) and leaves its CONJECTURE A4* open.

**Computational meaning.**
- *Exactly:* Gaussian and Gaussian-copula CMIs (unlock 46(d)) are ratios of minors, computable in O(|B|³). Along a nested family of separators they are updated by condensation, one exchange relation per step.
- *Tropically:* the CMI is the positive part of a log-ratio. That is the relevant form in the frozen regime, where determinants are dominated by single terms.
- Whether mutation (changing the separator) has a Markov-network meaning is open.

### 34. Total positivity and faithfulness [T, K]

**Statement (Fallat et al.).**
- An MTP₂ distribution that is also a graphoid is faithful to its pairwise independence graph (Thm 6.1).
- MTP₂ laws are closed under marginalisation and conditioning. They are upward-stable, singleton-transitive, compositional semigraphoids (Thm 5.3).
- A Gaussian is MTP₂ iff its precision matrix is an M-matrix (Karlin–Rinott).
- Coarsening can change the concentration graph (Example 6.2).

**Status.** THEOREM (verified).

**Computational meaning (guard).**
- Under total positivity, Markov structure can be read off pairwise data (faithfulness), so local tests certify global separations.
- The network's layers have signed weights, and its face law has correlations of both signs, so it is not MTP₂.
- Faithfulness therefore cannot be assumed; separations must be measured (unlock 45).

---

## 6. Heisenberg picture and Dirichlet forms [H]

### 35. States forward, questions backward ⊕ [H, F, K]

**Statement.**
- Let K_l be the Koopman operator of arrow l, (K_l g)(z_l) = g(relu(z_l) W_{l+1}), and put g_{l,j} = K_l ⋯ K_{L−1} r_j for the readout r_j = relu(z_{L,j}). Then E_{μ_L} r_j = E_{μ_l} g_{l,j} for every l.
- K_l is linear and positive on observables, although the forward map on laws is nonlinear in any finite parametrisation. In Note 1's language this is Exel's transfer operator: states forward, questions backward.
- At every layer there are exactly n questions, one per output unit. Under the linearised (mean-gate) transport, each is a ridge function along one direction ū_{l→L}(·, j) of layer-l space.

**Status.** THEOREM (duality of push-forward and pull-back; Exel's transfer operator as in Note 1) + DERIVED (the rank-n count).

**Computational meaning.**
- Any split point l gives an exact meet-in-the-middle identity.
- The forward state needs to be accurate only on the span of the n pulled-back questions, and the backward object has fixed size n per layer at every depth.
- Forward compression of old signed content in a fixed basis fails in the dictionary measured (unlock 18). Backward evaluation at full rank n costs only matrix products: the heisenberg design pulls slices back to the source layer with about four n × n matrix products per (source, target) pair, and no intermediate n³ tensor.

### 36. Duhamel telescoping and the bilinear (doubly robust) error ⊕ [H, B, K, F]

**Statement.**
- Let ν_l be any reference chain with ν_1 = μ_1 and ν_{l+1} = Π(T_l ν_l). Here T_l is the exact push-forward by arrow l and Π is any projection onto a tractable family. The heisenberg design takes the Gaussian law with the same mean and covariance; that makes covariance propagation, an existing estimator that BRIEF rule 1 admits only as a baseline, the forward backbone. The identity is equally valid with other Π: the product-of-marginals projection, whose chain is the Bethe lift (unlock 14), a face-law family, or a Gibbs family on histories (unlock 7). Put ρ_l = T_{l−1} ν_{l−1}. Then, exactly,

  E_{μ_L} r_j − E_{ν_L} r_j = Σ_{l=2}^{L} (ρ_l − ν_l)[g_{l,j}].
- For any model ĝ_l of the pulled-back question,

  truth − (reference + Σ_l (ρ_l − ν_l)[ĝ_l]) = Σ_l (ρ_l − ν_l)[g_l − ĝ_l].

**Status.** DERIVED (heisenberg design, Theorems 1–2; a one-line telescoping). It is the commutative analogue of Chen–Rouzé's telescoping against the erased state (A-unification §5).

**Computational meaning.**
- *With a controlled error that is a product:* the error of reference plus correction is bilinear, (local defect of one arrow) × (model error of the question).
- A cheap forward reference and a cheap backward model therefore give an error of the order of the product of their errors. Example, with the Gaussian closure as reference: the total correction is the closure's error, ≈ 2.1e-3 rms per unit at n = 1024, and the bar is ≈ 1.3e-4 (unlock 3). If the defects are computed exactly from the reference, a question model with ≈ 6 % relative error meets the bar. Measured: the first-order heisenberg realisation removes a factor ≈ 10 in raw MSE at n = 1024, i.e. it leaves ≈ 30 % of the correction (CONVERGENCE.md).
- All old content sits inside the exact g_l, so the identity has no memory truncation. Truncations enter only through ĝ.

### 37. Centred wall identity: the Malliavin–Stein form of the per-unit mean ⊕ [H, T, F]

**Statement.** Let Z = z_{l,j}(x), a Lipschitz function of the Gaussian input with a density continuous at 0, and let L be the Ornstein–Uhlenbeck generator. Then

 E relu(Z) = μ_Z P(Z > 0) + E[δ(Z) Γ_Z],
 Γ_Z = ⟨∇Z, −∇L^{−1}(Z − μ_Z)⟩ = ∫_0^1 ⟨∇Z(x), E'[∇Z(u x + √(1−u²) x')]⟩ du,

with x' an independent copy of x.
- E Γ_Z = Var Z exactly.
- E[Γ_Z | Z = t] is the Stein kernel τ_Z(t) of unlock 3(c).
- Equivalently, by the Houdré–Pérez-Abreu covariance identity,

  E[(Z − μ) 1(Z > 0)] = ∫_0^1 E[δ(Z(x)) ⟨∇Z(x), ∇Z(x_u)⟩] du.

  This is a wall integral of the two-replica overlap of input-space gradients at x and a u-correlated copy x_u.
- In Wiener-chaos terms, the centred weight counts every chaos once, E Γ_Z = Σ_q ‖J_q Z‖². The uncentred Euler–Stein weight counts the q-th chaos q times, E‖∇Z‖² = Σ_q q ‖J_q Z‖².

**Status.**
- THEOREM for the ingredients (Nourdin–Peccati Stein-kernel formula; Houdré–Pérez-Abreu covariance identity; Mehler representation of −DL^{−1}; from memory). DERIVED for the application to the ReLU readout.
- Checked on a hot unit (C8):
  - E Γ_Z = 0.40966 against Var Z = 0.40967;
  - wall form 0.2478 against 0.2473, with Monte Carlo error ≈ 7e-4;
  - the uncentred weight's mean is 0.550.
- CONJECTURE: at He initialisation the centred wall correction p_Z(0)(τ_Z(0) − Var Z) is of the order of the non-Gaussianity of Z, O(n^{−1/2}) per unit, so the centred split is perturbative where the uncentred one is not (unlock 4, guard).
  - Support: τ_Z ≡ Var Z for a Gaussian, and the Stein discrepancy E|τ_Z(Z) − Var Z| controls the distance to the Gaussian (Nourdin–Peccati).
  - Test: §11, tropical.
  - Constraint from measurement: the face-mass and wall defects together make the Gaussian closure's error, ≈ 2/n rms per unit for 128 ≤ n ≤ 1024 (unlock 3). Either each defect is O(1/n) at these widths, or their O(n^{−1/2}) parts cancel. The test should report the two defects separately, with their width scaling.

**Computational meaning.**
- *Exactly:* the per-unit mean is drift × face mass + wall density × the wall-conditional mean of a weight whose global mean is Var Z. The weight carries no mean spike.
- At the level of the identity, this repairs the obstruction that made the uncentred split non-perturbative (unlock 4, guard).
- The Gaussian closure is the approximation "τ_Z(0) = Var Z, with Gaussian face mass and wall density".
- The exact correction is a wall–overlap correlation: how the two-replica overlap of input-space gradients differs on the wall from its average. Equivalently, it compares the linear maps M_μ of the histories met by x and x_u. That is a pair-of-histories (replica) object of the face algebra, not a moment.

### 38. Price's theorem: covariance derivatives are face masses one codimension down [H, F, M]

**Statement.**
- For Z ~ N(m, C) and g of moderate growth: ∂E g(Z)/∂C_ij = E[∂_i ∂_j g(Z)] for i ≠ j, ½ E ∂_i² g for i = j, and ∂E g/∂m_i = E ∂_i g.
- Along the interpolation C_t = (1−t) C_0 + t C_1 (the Slepian–Kahane smart path):

  E_{C_1} g − E_{C_0} g = ½ ∫_0^1 Σ_ij (C_1 − C_0)_ij E_{C_t}[∂_i ∂_j g] dt.
- For g = ∏_{i∈S} 1[z_i > 0] or ∏ relu(z_i), ∂_i ∂_j g involves δ(z_i) δ(z_j): a codimension-2 face density (Plackett's reduction, unlock 9).

**Status.** THEOREM (Price 1958; Plackett 1954; Kahane 1986; from memory). The derivative of Sheppard's formula, 1/(2π√(1−ρ²)) = φ_2(0, 0; ρ), is an instance.

**Computational meaning.**
- *Exactly:* the sensitivity of any face mass or ReLU expectation to the joint structure of a Gaussian layer is a lower-dimensional face mass.
- Perturbation and interpolation estimators, which move from a tractable covariance to the true one, cost one codimension at a time.
- This is also the precise sense in which new non-Gaussian structure is born only at walls (unlock 30).

### 39. Dirichlet forms are squared derivations; Lipschitz constants of the questions govern error propagation ⊕ [H, B, K]

**Statement.**
- (a) Completely Dirichlet forms on the KMS-embedded L² space correspond one-to-one to KMS-symmetric Markov semigroups (Cipriani; Goldstein–Lindsay).
  - Noncommutative Dirichlet forms are squares of twisted derivations (Vernooij–Wirth).
  - In the commutative case the Connes distance is a Kantorovich–Rubinstein W_1 distance. The digest has this for a complete Riemannian spin manifold, with the geodesic metric (Dirac operator). For a Dirichlet form with carré du champ Γ, sup{|φ(f) − ψ(f)| : Γ(f) ≤ 1} is W_1 for the intrinsic metric when Γ(f) ≤ 1 characterises the 1-Lipschitz functions, as for strongly local forms (from memory).
- (b) Hence |ω(q) − ω̃(q)| ≤ Lip(q) · W_1(ω, ω̃) for every observable q.
  - For a pulled-back question, Lip(g_{l,j}) is bounded by the norms of the gated propagators from l to L.
  - Measured gains (transfer-spectrum §2.6): a generic direction has mean-square gain 2E[Φ²] ≈ 0.58–0.95 per layer, a contraction (0.5 at layer 0; one measured layer at 1.07). The mean direction is transported with 3.1–7.8 times the bulk gain after 8–15 layers (n ≥ 256).

**Status.** (a) THEOREM (as read in the nc-Dirichlet digest). (b) DERIVED (Kantorovich–Rubinstein duality) + measured.

**Computational meaning.**
- *With a controlled error:* a state error introduced at layer l reaches the output multiplied by the Lipschitz constants of the pulled-back questions.
- Errors in the bulk are damped layer by layer; errors along the mean direction are amplified.
- An error budget should therefore go first to the mean direction, which is rank one and exact by unlock 3(a), and then to the bulk at the layers nearest the output.

### 40. Positivity quantified: a hereditary certificate for forgetting and local gaps ⊕ [H, K, B, F]

**Statement (B-programme §6).**
- (a) If p_i/q_i ∈ [c, c e^D], then TV(p, q) ≤ tanh(D/4), and the bound is sharp.
- (b) For a positive chain with kernels A_l > 0, every Doob transform has Dobrushin coefficient ≤ tanh(Δ(A_l)/4). This covers every pinning and every boundary change. Δ is the projective (Hilbert, cross-ratio) diameter, and the certificate is hereditary.
- (c) The single-site window sampler has absolute spectral gap ≥ (1 − 2 tanh(Δ/4))/w for every window length w and every boundary condition. Block samplers of length b > e^{Δ/2} − 1 have gap ≥ (b − (e^{Δ/2} − 1))/(w + b − 1) (R6).
- (d) The restriction of the KMS state to the first t layers depends on the terminal condition at most through ∏_{l ≥ t} tanh(Δ_l/4) in total variation (R7). This has the shape of U1 but is not U1; the B-programme's critique pass relabelled that reading ANALOGY.

**Status.** DERIVED (B-programme Lemmas 6.1–6.3, Thm 6.4, Cor 6.5; checked there). The ingredients are classical: Birkhoff's contraction coefficient, and Dobrushin's condition run through the Dyer–Goldberg–Jerrum path coupling (B-programme §6).

**Computational meaning.** In the positive sector (face masses, Gibbs laws on histories, the rectified mean direction), an estimator can:
- truncate history after a depth computed from one layer kernel, with a guaranteed total-variation error;
- resample windows with a guaranteed gap;
- condition freely, since conditioning never spoils the certificate.

**Guard.** Positivity only.
- The old-content transport is a signed sum over histories through open gates and lies outside (B-programme §11). No gap is predicted there, and none is observed (unlock 18).
- Dictionary v1's fitted kernels have zero transitions, so Δ = ∞ at a single layer. The onset depth of positivity must be reported.

### 41. Recovery is not mixing: Cesàro means, energy duality and window re-routing [H, K]

**Statement.**
- (a) For a reversible generator 𝓛 with Cesàro mean R_t = t^{−1} ∫_0^t e^{s𝓛} ds:
  - R_t → Π, the conditional expectation onto the jump components;
  - E(R_t f) ≤ 0.41 t^{−1} ‖f‖², with no spectral gap. The constant is sup_u (1 − e^{−u})²/u ≈ 0.4073.
- (b) For window re-routing on histories, reversibility is exactly KMS_β quasi-invariance with cocycle e^{−βc_F}.
  - If the window fibres are flip-connected, Π = E_P[· | outside the window], and every Q that agrees with P outside the window is recovered exactly: QΠ = P.
  - For Gibbs laws the recovery map sees only the two faces adjacent to the window (R5).
- (c) In the noncommutative case recovery obeys ‖R_{A,t}[ρ_{−A}] − ρ‖_1 ≤ √(0.41/t) χ_KMS(ρ_{−A} ‖ ρ) κ(A)^{−1/2}. Here κ(A) is the Poincaré constant relative to N_A. It stays O(1) while the literal local gap collapses exponentially (A-unification §0 item 5; [CR]).

**Status.** DERIVED (A-unification §5; B-programme §5) + THEOREM (Chen–Rouzé, as read in the digest).

**Computational meaning.**
- A conditional expectation that cannot be written down can be approximated by a time-averaged local dynamics. The error decays like t^{−1/2} without a gap, provided the dynamics is coercive relative to the target algebra.
- In the commutative diagonal resolution the recovery is exact and local. This matters only for states with coherences (unlock 49).

### 42. A gap turns a global projection into a short polynomial of local operators [H, B]

**Statement.**
- Let P be a self-adjoint contraction with spectral gap λ, i.e. spectrum in [−1, 1 − λ] ∪ {1}.
  - The projection onto its top eigenspace is approximated in norm to ε by a polynomial of degree O(λ^{−1/2} log(1/ε)) in P (Chebyshev acceleration).
  - More crudely, by P^k with k = O(λ^{−1} log(1/ε)), provided the gap is absolute (spectrum in [−1 + λ, 1 − λ] ∪ {1}); with −1 in the spectrum use ((I + P)/2)^k.
- In coarse geometry: with a gap, the global projection is a norm limit of finite-propagation operators (Roe algebra; Kazhdan projections). For an expander this projection is a non-compact ghost: it lies in the Roe algebra, yet no finite window sees it, and its K-theory class is not assembled from local data (coarse Baum–Connes fails surjectivity; Higson, Higson–Lafforgue–Skandalis, Willett–Yu, as read in the expanders digest).

**Status.** THEOREM (Chebyshev, standard; Roe and Kazhdan as read in the expanders digest; the degree bound from memory).

**Computational meaning.**
- Where a gap exists (the positive sector, the Perron mean mode), extracting the stationary part costs O(λ^{−1/2} log(1/ε)) applications of local operators.
- Where there is no gap (the signed bulk), no polynomial shortcut exists, and the content must be carried (unlocks 18, 35).

---

## 7. Markov networks, CMI and exchange relations [K]

### 43. Hammersley–Clifford is discrete integrability [K, F]

**Statement.**
- A positive law on a product space is Markov for a graph G iff Δ_i Δ_j log P = 0 for every non-adjacent pair (i, j), where Δ_i is a discrete difference in coordinate i.
- Möbius inversion of log P over subsets (the blackening algebra) gives the interaction potentials Φ_A, which vanish off the cliques.
- Approximately: if every non-adjacent mixed difference is ≤ δ, then ‖Φ_A‖ ≤ 2^{|A|−2} δ for every non-clique A (hammersley-clifford digest).
- Pinning (conditioning on a set) preserves the Markov property with the induced graph. Marginalising fills in: it adds edges between the neighbours of the removed set.
- Quantum version: Gibbs states of commuting clique Hamiltonians are quantum Markov networks. The converse (full rank, Markov ⇒ commuting clique Hamiltonian) holds on chains, trees and triangle-free graphs, and fails with triangles (Brown–Poulin; Leifer–Poulin; hammersley-clifford digest).
- Intersection holds quantitatively, with Friedrichs-angle constant 1/(1 − c_F).

**Status.** THEOREM (Hammersley–Clifford; Möbius) + DERIVED (the pointwise approximation and quantitative intersection, in the hammersley-clifford digest and A-unification).

**Computational meaning.**
- *Exactly:* whether a carried state may be represented by local potentials is decided by local mixed differences.
- An approximately Markov state has approximately clique-supported potentials, with an explicit constant.
- Pinning is free and marginalising costs fill-in, so a design should condition (pin) rather than marginalise wherever it can.

### 44. Depth is exactly Markov at full resolution; memory is fill-in ⊕ [K, F, B]

**Statement.**
- z_{l+1} = relu(z_l) W_{l+1} is a deterministic function of z_l, so (z_1, …, z_L) is a (degenerate) Markov chain for every network.
- The face process σ_l = sgn z_l is a function (lumping) of it. It is Markov for every initial law iff the chain is strongly lumpable (Kemeny–Snell); for the given input law the condition is weak lumpability (Rosenblatt; from memory). Otherwise its memory is created by marginalising the coordinates within faces: fill-in (unlock 43).
- Every Gibbs law on histories is Markov at every layer (unlock 7).
- Markov at a layer is equivalent to a commuting square of the past, future and present conditional expectations inside the history algebra D (U8; the commutative case is proved).

**Status.** DERIVED (elementary; U8 as in the hdx digest).

**Computational meaning.**
- Memory is a property of the resolution, not of depth. The full-resolution state, the law of z_l, is a perfect separator, so an estimator that carried it exactly would need no history.
- Every coarser state (faces; a Gaussian field plus sites; windows) has a memory cost, and that cost is exactly computable (unlock 45).
- Measured:
  - dictionary v1's face process has 25–35 % memory (mlp-bridge);
  - about 40 % of the pairwise third-order structure feeding the next layer, the (2,1) slice κ3(z_a, z_a, z_b), is older than one layer (BRIEF §3), a measurement made in the cumulant dictionary.

  Both are charged to the resolution of the dictionary.

### 45. Memory = Σ CMI = distance to the Gibbs family = Σ Petz defects ⊕ [K, H, F]

**Statement.**
- (a) Let P be any law on (σ_0, …, σ_L), with Markov projection Q(σ) = P(σ_0) ∏ P(σ_{i+1} | σ_i). Then D(P ‖ Q) = Σ_i I(σ_{i+1} : σ_{<i} | σ_i). This equals the relative-entropy distance of P from the closure of the KMS family (U7; checked to 10 digits in the hdx digest, C2).
- (b) The following are equivalent (B-programme Prop 4.1):
  - P is Markov at layer l;
  - σ_l is Fisher–Neyman sufficient for the family whose parameter is the past;
  - I(past : future | σ_l) = 0;
  - the (past, σ_l) subalgebra is Petz-sufficient for a pair of states. Its defect is exactly that CMI.

  "Markov" and "sufficient for the Gibbs family" concern two different Petz families and are logically independent (R4).
- (c) Junction-tree form. Take blocks V_1, …, V_k with separators S_t ⊆ V_{<t}. The projection p_JT = ∏_t p(V_t | S_t) satisfies D(p ‖ p_JT) = Σ_t I(V_t : V_{<t} ∖ S_t | S_t) (markov design P3). The Gaussian case is checked in C2.

**Status.** THEOREM/DERIVED (the chain rule; Petz and Jenčová–Petz; as cited).

**Computational meaning.**
- *With a controlled, exactly budgeted error:* the price of carrying only a separator is a sum of CMIs in KL, one per cut. The separator can be one layer's face state, a window, or a quasi-free field plus sites.
- These CMIs are measurable from samples at small width, and computable in closed form for quasi-free states (unlock 46).
- The window length of a depth junction tree is chosen by the measured decay of these CMIs with age.

### 46. Gaussian exchange relations: Dodgson, Koteljanskii, chordal determinants, copula invariance [K, M, T]

**Statement.**
- (a) Desnanot–Jacobi: det M[aBc] det M[B] = det M[aB] det M[Bc] − det M[aB|Bc] det M[Bc|aB], where the subtracted term is a square for symmetric M. Hence I(a : c | B) = −½ log(1 − r²) = ½ log(1 + y) (unlock 33).
- (b) Koteljanskii: for positive definite M and index sets S and T,

  I(S∖T : T∖S | S∩T) = ½ log(det M_S det M_T / (det M_{S∪T} det M_{S∩T})) ≥ 0.
- (c) For a chordal graph with cliques 𝒞 and separators 𝒮 (a junction tree):

  log det M ≤ Σ_C log det M_C − Σ_S log det M_S,

  with equality iff the Gaussian is Markov on the graph (iterated Koteljanskii; equivalently, maximum entropy given the clique marginals). Half the gap, ½(Σ_C log det M_C − Σ_S log det M_S − log det M), equals the KL divergence to the Markov projection and the sum of the cut CMIs.
- (d) CMI is invariant under injective transformations of each coordinate. So for a Gaussian copula (z_i = T_i(g_i), T_i increasing, g Gaussian), the CMI of z is the Gaussian CMI of the latent correlation.

**Status.** THEOREM ((a), (b), (d) classical; (c) maximum entropy) + checked (C1; C2, where half the log-det gap, 1.2277, equals Σ CMI = 1.2277, and a tridiagonal precision gives 0 to 4e-16).

**Computational meaning.**
- *Exactly and cheaply:* for every quasi-free or copula state, every cut CMI is a ratio of determinants. The Gaussian parts of the faces, signings and markov designs' states are of this kind.
- The memory budget of unlock 45 then costs O(s³) per separator of size s.
- The Markov defect can be monitored inside the estimator at no sampling cost.

**Guard.** (d) needs a Gaussian copula. It can hold for pre-activations (exactly at layer 1; approximately where the state is modelled as a copula), never for activations: relu is not injective, so the CMIs of a_l and z_l differ.

### 47. Separator along depth, tree across width: pseudorandomness is the switch ⊕ [K, B, M]

**Statement.**
- (a) Separator gluing (junction trees, Hammersley–Clifford, Yang's CMI bounds) is exact, or controlled by CMI. For a general law it costs exp(separator size) (treewidth). For a quasi-free (Gaussian) law, conditioning on a separator is a Schur complement, O(s³), and its CMI is a determinant ratio (unlock 46).
- (b) Tree gluing (Bethe, cavity, TAP) is approximate, controlled by the decay of correlations along computation trees.
  - It is strong on pseudorandom geometry, and weak at low temperature and on amenable geometry.
  - Expanders are the worst case for separators, since Yang's boundary gain vanishes on them. They are the best case for trees (unlocks 12, 16) under correlation decay: uniqueness for marginals, unfrustrated models for free energies (foundations-input-verified.md §3.2, T3–T4, and §3.3). The competition scores marginals.
- (c) The network supplies both geometries.
  - Depth is a path of exact separators (unlock 44), but each separator is a whole layer: size n, and not quasi-free because of the faces.
  - Width is a dense pseudorandom bipartite geometry through fresh weights, where only tree-shaped contractions survive, as far as measured (unlock 13: layer 2).

**Status.** DERIVED (a synthesis of the cited theorems: the user's Idea B, checked against the programme's results).

**Computational meaning.** The design rule it admits under the quasi-free reading (see the bias guard below):
- Glue along depth through quasi-free separators: Gaussian or Gaussian-copula fields, carried exactly by linear algebra at O(n³) per layer.
- Attach the non-quasi-free content (faces, hubs, unary non-Gaussian potentials) as local decorations, glued across width by tree rules with the Onsager correction.
- Charge the residual to two measurable defects: the cut CMIs along depth (unlock 45) and the cycle content across width (unlock 13).

**Guard (bias).** Read literally, this rule is covariance propagation with local cumulant decorations (hubs, sources): the structure of the existing estimators that BRIEF rule 1 admits only as baselines, and the structure on which all six fresh-slate designs converged (CONVERGENCE.md). The principle does not single out quasi-free separators. That choice is made for linear-algebra convenience and is a dictionary item. Separators of other kinds (face states, unlock 8; Gibbs laws on histories, unlock 7; coherent states, unlock 49) are equally admissible and have not been taken to an estimator.

**Guard.** Old content (BRIEF §3) is a cross-layer joint structure with no low-rank form. Neither rule is local there; the Heisenberg evaluation (unlocks 35–36) is the alternative.

### 48. Two defects on one stage: angles and cocycle leakage [K, H, B]

**Statement (A-unification).** Fix a faithful state and its forget-A algebras N_A.
- (a) The angle c(A, B) = ‖E_A E_B − E_{A∪B}‖ behaves as follows:
  - it is zero on separated pairs iff the field is Markov;
  - it decays across a buffer under strong spatial mixing;
  - it equals λ(G)/d on the two ends of an edge (expander mixing);
  - on adjacent binary spins it is the largest conditional correlation, at most the geometric mean of the two Dobrushin influences.

  The sharp two-projection inequality is (1 − c) Var_{A∪B} ≤ Var_A + Var_B. Alternating projections converge at rate ‖(E_A E_B)^k − E_{A∪B}‖ = c^{2k−1}, where E_{A∪B} projects onto N_A ∩ N_B (Aronszajn; Kayalar–Weinert; from memory).
- (b) For quantum states, I(A : C | B) = 0 iff the Connes cocycle ρ_{BC}^{it} ρ^{−it} stays in M_{AB} ⊗ 1. Quantitatively, I ≤ (1/α + 3) d_A^{2α/(1+α)} q^{2α/(1+α)}, with q the KMS-strip-weighted leakage and d_A² the Jones index.

**Status.** DERIVED (A-unification §3–4) + THEOREM (as marked).

**Computational meaning.**
- *With a controlled error:* gluing two local conditional expectations loses a factor 1/(1 − c) in variance. Examples are two overlapping windows, two sub-frames, or a forward and a backward sweep.
- Iterated local corrections (alternating projections, forward–backward sweeps) converge geometrically at the Friedrichs-angle rate. The number of sweeps an iterative estimator needs is therefore predicted by one angle.

### 49. Where noncommutative Markov tools would carry content: coherent histories [K, H]

**Statement (B-programme Lemma 8.3, Conjectures 8.4–8.5).**
- Across a depth cut at layer a, ℓ²(histories) = ⊕_{σ ∈ K_a} H^<_σ ⊗ H^>_σ. The Gibbs state of an additive cocycle is a classical mixture, over the cut face, of product states: it is exactly Markov.
- Conjecture 8.4: let H = H_F + V, with window-local coherent terms V_X (|X| ≤ R_0, ‖V_X‖ ≤ J, bounded overlap). Then I_ρ(A : C | B) ≤ C_β exp(C_β g − c_β (b − a)) across a buffer of b − a layers.
- Conjecture 8.5: window Lindbladians have a noncommutative uniform local gap under a finite noncommutative projective diameter.

**Status.** DERIVED (8.3) + CONJECTURE (8.4, 8.5).

**Computational meaning.**
- Chen–Rouzé and Yang add something beyond the commutative unlocks only for a design whose state carries coherences between histories, i.e. a non-diagonal state on C*(Λ).
- For such a state, CMI decay across a buffer would bound the error of forgetting everything beyond a window.
- No current design carries such a state. This is the precise place where the user's NCG direction would enter an estimator.

### 50. Agreement and local-to-global reconstruction of states [K, F, B]

**Statement.**
- U4 (CONJECTURE): suppose the coherence complex of a state's restrictions is a c-agreement expander. Then the state is determined, up to ε/c, by its restrictions to a bounded-degree family of faces that agree pairwise on overlaps up to ε.
  - Reconstruction would then need O(|V|) restrictions instead of all faces.
  - Dinur–Kaufman proves this for {0,1}-valued functions; the version for states is open. The sheaf obstruction is Abramsky–Brandenburger contextuality.
- U6 and U9: local-to-global mixing from pinned spectral independence holds for arbitrary laws. External fields are tilts by e^h, with h in the face algebra.

**Status.** CONJECTURE (U4) + interpretation backed by theorems (U6, U9, hdx digest).

**Computational meaning.**
- If U4 is true, a global face law can be glued from O(n) pairwise-consistent local face laws, with an error linear in the local inconsistency. That would certify pairwise face data as a sufficient state.
- Test: at widths 64–256, compare the pairwise consistency of a layer's local face laws with the global error of the glued law.

---

## 8. Cross-cutting unlocks

● = the unlock is used essentially by that stream.

| # | unlock | F | B | M | T | H | K | what it gives every stream it touches |
|---|---|---|---|---|---|---|---|---|
| 1 | face-history decomposition, radial factorisation | ● | | | ● | ● | ● | the estimand as a linear functional of a vector measure on histories |
| 2 | ReLU = Lüders conditioning | ● | | | | ● | ● | the only nonlinearity is conditioning |
| 3 | ridge reduction, exact mean recursion, two scalars per unit | ● | ● | | ● | ● | ● | the exact target: face mass and wall term per unit per layer |
| 4 | barycentres as boundary integrals; mean = wall mass | ● | ● | | ● | ● | | the exact interior/boundary duality (uncentred) |
| 5 | fresh randomness: twirl, one Kraus operator, sign sectors | | ● | ● | | ● | ● | what is annealed-computable, and why the score is not |
| 13 | power counting for a fresh layer | | ● | ● | | | ● | tree contractions across width suffice (measured at layer 2); joint structure of the state is necessary |
| 14 | Bethe = cover limit; network lift limit | ● | ● | ● | | | | an exact, computable tree reference |
| 22 | sign-averaging = even-sector projection | | ● | ● | ● | | | the sign average keeps the weight magnitudes, drops the mean alignment, and is neither the annealed nor the Bethe object |
| 29 | signs are the whole answer at T = 0 | | | ● | ● | ● | | closes the dominant-cone/path route |
| 30 | frozen units are linear arrows | ● | | | ● | ● | | where face-splitting work lives; exact threshold t ≳ 4 |
| 33 | CMI = ½ log(1 + y) | | | ● | ● | | ● | exchange relations, their degeneration and tropical limit |
| 35 | states forward, questions backward | ● | | | | ● | ● | n questions per layer, at full rank |
| 36 | Duhamel telescoping, bilinear error | ● | ● | | | ● | ● | error = product of two small errors |
| 37 | centred wall identity | ● | | | ● | ● | | the per-unit mean without the mean spike in the weights |
| 39 | Lipschitz constants of pulled-back questions | | ● | | | ● | ● | where state errors are amplified (mean direction) |
| 40 | positivity quantified | ● | ● | | | ● | ● | certified forgetting in the positive sector only |
| 44 | depth exactly Markov at full resolution | ● | ● | | | | ● | memory is a property of the resolution |
| 45 | memory = Σ CMI = KL to the Gibbs family | ● | | | | ● | ● | an exact error budget for any separator state |
| 47 | separator along depth, tree across width | | ● | ● | | | ● | the gluing rule for each direction of the network; its quasi-free reading is a dictionary choice (guard) |

---

## 9. The five most important statements

**I. The estimand has an exact skeleton: a linear mean recursion plus two wall/face scalars per unit (unlocks 3, 37, with 1 and 4).**
- E z_{l+1} = E[a_l] W_{l+1} exactly.
- For every unit, E relu(z) = μ P(z > 0) + E[δ(z) Γ_z]. Here Γ_z = ⟨∇z, −∇L^{−1}(z − μ)⟩ is a two-replica overlap of input-space gradients, and its mean is exactly Var z.
- Gaussian closure is the approximation "Gaussian face mass, Gaussian wall density, Γ ≡ Var z". Given the exact drift, its error splits exactly into drift × (face-mass defect) + (wall defect). It is ≈ 2/n rms per unit for 128 ≤ n ≤ 1024, ≈ 16× the bar at n = 1024 (unlock 3).
- Globally, E f = E Δf: the mean is the Gaussian mass of the tropical hypersurface, which by summation by parts is the statement that barycentres of faces are facet integrals.
- The uncentred form of this identity puts the mean spike and the high-chaos roughness of depth into the wall weights and is non-perturbative (tropical R-E2). The centred form does not.

**II. Fresh randomness organises the width: annealed = twirl, sign average = even sector, and the score has three quenched parts (unlocks 5, 13, 22, 12).** W_{l+1} is independent of the law of a_l. Consequently:
- the annealed one-step map is a perfect expander;
- averaging over the fresh layer's signs keeps exactly the even sector, the network analogue of Godsil–Gutman; it keeps the quenched magnitudes |W| and is not a tree object;
- the quenched per-unit information has three parts: an O(1) coherent sign-odd alignment with the mean direction, which the exact recursion carries; an O(n^{−1/2}) incoherent sign-odd pairing of fresh weights with the state's joint structure; and an O(n^{−1/2}) sign-even part carried by the magnitudes, as large as the second (C9).

Loops closed among the units a fresh layer touches are suppressed by the state's typical pairwise correlation (n^{−1/2} up to a factor that grows with depth; measured at layer 2 only). So across width the gluing rule is tree-shaped (Bethe/TAP with the Onsager term, certified only by ℓ² quantities), and everything beyond the mean recursion and the quenched magnitudes lives in the state's joint structure.

**III. Depth is exactly Markov at full resolution; any coarser state pays a computable sum of CMIs (unlocks 44–47).**
- Memory is a property of the resolution. Its price is D(P ‖ P_JT) = Σ_cuts CMI. This equals the KL distance to the Gibbs-on-histories (KMS) family and the sum of the Petz sufficiency defects.
- For quasi-free (Gaussian or Gaussian-copula) separators, conditioning is a Schur complement, and each CMI is a ratio of minors: for single coordinates, I = ½ log(1 + y), with y the ratio of the two monomials of the Desnanot–Jacobi exchange relation, and Markov is the degeneration y = 0; for blocks, the Koteljanskii ratio (unlock 33).
- One design rule follows if separators are taken quasi-free: glue along depth through Gaussian or copula fields, attach the face (non-quasi-free) content as local decorations glued across width by trees, and set the window length by the measured CMI decay. That reading reproduces covariance propagation with cumulant decorations, the existing estimators' structure (unlock 47, guard). The quasi-free choice is a dictionary item, not part of the principle.

**IV. States forward, questions backward, with a bilinear error (unlocks 35, 36, 39).**
- The pairing of the law at layer l with the pulled-back readout is invariant. Duhamel telescoping writes the error of any forward reference exactly as Σ_l (local defect of arrow l)[exact pulled-back question].
- With a model of the question, the remainder is bilinear: defect × model error.
- There are exactly n questions per layer, so backward evaluation works at full rank. Forward compression of old signed content in a fixed basis fails in the dictionary measured (k_ε/n is constant; unlock 18).
- State errors reach the output through the Lipschitz constants of the questions: contracted in the bulk (gain 2E[Φ²] ≈ 0.58–0.95 per layer; one measured layer at 1.07) and amplified along the mean direction.

**V. Positivity is the only certified source of forgetting and compression; the scored content is signed (unlocks 40, 18, 29, 12).**
- The positive sector has explicit, hereditary certificates. This sector is the face masses, Gibbs laws on histories, and the rectified mean direction. Its certificates are TV ≤ tanh(Δ/4), uniform window gaps, forgetting ∏ tanh(Δ_l/4), and one Perron/BBP mode.
- The signed sector has none:
  - the free-probability bulk has no gap (PR ≈ n/(2·age), k_ε/n constant);
  - the path sum is a weak-disorder polymer with no dominant path;
  - the zero-temperature Newton-polytope form cancels by 10²³ at n = 1024.
- No design may rely on decay, low rank, or a tropical skeleton of the signed sector in the dictionaries measured so far (propagator, covariance and HOSVD bases for old cumulant slices; dominant cone or path). There it must be carried at full rank or evaluated backward; a new dictionary would have to be tested afresh (BRIEF rule 3). On the dense pseudorandom geometry only spectral (ℓ²) certificates are usable; Dobrushin-type (ℓ¹) ones fail by √n.

---

## 10. Guards, and negative results charged to the dictionary

Each guard is stated with the dictionary item it is charged to (BRIEF rule 3). None of them refutes a general statement.

- **Uncentred wall split** (unlock 4; tropical R-E2). Charged to the dictionary "walls weighted by input-space slope²". The centred identity (37) is the repair at the level of the identity. Whether its wall correction is O(n^{−1/2}), or O(1/n) as the closure's measured error suggests, is the CONJECTURE stated in unlock 37, tested as in §11.
- **The zero-temperature skeleton** (unlocks 29, 30, 32). The 10²³ cancellation, the weak-disorder polymer, and hot gates (a quarter of them at depth 16) are charged to "dominant cone or dominant path as the state". The tropical content of the theory survives as the signed wall identities and as a cost saving on the frozen sub-frame.
- **Forward compression of signed content** (unlocks 18, 19): no gap, k_ε/n constant, no low-rank form below ≈ 0.3 n. Charged to "carry old content in a fixed small basis". The theory's answer is backward evaluation (unlocks 35–36).
- **Dictionary v1** (unlocks 8, 44): 25–35 % face-process memory, and a residual of 0.42–0.5 against coboundaries. Charged to the resolution of v1. Unlock 45 prices it exactly.
- **Positivity hypotheses** (unlocks 7, 8, 40, 41). The certificates apply to Gibbs laws on histories and to positive kernels. Signed transport is outside them, which is a statement about scope, not a failure.
- **ℓ¹ certificates on dense layers** (unlock 11): Mooij–Kappen, Dobrushin, and BP certificates based on |J| fail by √n. Use ℓ² certificates (unlock 12).
- **Spectral independence at η ≈ 7** (unlock 17) gives a bounded variance factor, not precision. Mixing certificates built on it are numerically vacuous.
- **Determinantal, real-stable and totally positive tools** (unlocks 24, 26, 34). K_{n,n} is not Pfaffian, and the face law is neither strongly Rayleigh nor MTP₂. These tools apply to matching-type objects only; the network's sign averages are not such objects (unlock 22(c)).
- **Sign averages** (unlock 22) are wrong per unit at O(1), and they are not annealed objects: they keep the quenched magnitudes. They are references, never estimators.
- **Determinantal resummation** (unlock 23). The determinant counts vertex-revisiting walks; the hard-core paths-and-cycles sector is an α-permanent. Charged to the dictionary "quasi-free = degree-≤ 2 sector" (signings design P3(i)).
- **Copula invariance** (unlock 46(d)) does not pass through relu: it can hold only for pre-activations, and only where they are (modelled as) a Gaussian copula.
- **Annealed asymptotics** (unlock 31): 9π²/(2l²) overestimates 1 − c_16 by a factor 2.5. Use the iterate.
- **The user's Idea C** (a periodic quantum circuit from quantum cluster mutations). The construction is THEOREM-level (Kashaev–Nakanishi; foundations-input-verified.md §4). No computational meaning for quenched means was found; it is not used.

---

## 11. Messages to the design streams

- **faces.**
  - Carry the per-unit triple (face mass, wall density, Stein kernel at the wall). With the exact mean recursion it gives the readout exactly (unlock 3(c)).
  - Barycentres are facet integrals (unlock 4(a)).
  - In test T1, split the one-step error into μ × ΔP and Δ(p τ): they have different origins (face masses vs wall overlaps).
- **tropical.**
  - Redo R-E2 with the centred split of unlock 37, i.e. P(z > 0) and E[δ(z) Γ_z] against their Gaussian values, instead of I2's uncentred split. The centred weight has mean exactly Var z.
  - Exact linearisation at the bar needs t = |μ|/s ≳ 4 (unlock 30).
  - 1 − c_16 = 0.071, not the asymptotic 0.17 (unlock 31).
- **bethe.**
  - The random-lift limit (unlock 14(b)) is an exact, computable tree reference whose non-Gaussian single-site laws are exactly computable; they are the lift's, not the network's. Use it to measure loop content as quenched minus lift; that difference already contains the covariance between parents, at the leading quenched order.
  - Dense certificates must be ℓ², with the Onsager term (unlock 12).
- **signings.**
  - The sign-sector decomposition (unlocks 5(d), 22) is the theorem behind §1 Step A. Two corrections to DESIGN.md: a block of three legs on one index (W³) has zero mean and is unpaired, so the {3} sector of κ3 is O(n^{−1}), not O(n^{−1/2}) (the bethe design's table agrees; n·rms is constant from n = 256 to 4096); and P3(i) equates the degree-≤ 2 sector with a determinant, which unlock 23 corrects.
  - The Pfaffian guard (unlock 24) closes exact determinantal evaluation on dense layers.
  - Signed-lift control variates (unlock 25) need an identity that makes the network's mean a sign average before they can be used.
- **heisenberg.**
  - The Lipschitz gains of unlock 39 give an error budget: amplified along the mean direction, contracted in the bulk.
  - The centred identity (unlock 37) is the readout question in Gaussian-space form.
  - Theorem 2's bilinearity means ĝ needs ≈ 6 % relative accuracy on the correction (≈ 1.3e-4 absolute on a ≈ 2.1e-3 correction at n = 1024), not the bar's absolute accuracy on the mean. The first-order realisation leaves ≈ 30 % (unlock 36).
- **markov.**
  - For the quasi-free (Gaussian, copula) part of the state, the CMI budget is closed-form (unlock 46). Monitor it inside the estimator.
  - Measure cut CMIs against age to set d (unlock 45).
  - Copula invariance can hold for pre-activations only, and only where they are a Gaussian copula (exactly at layer 1).

---

## 12. Numerical checks

Run in the session scratchpad with numpy (scipy is not installed in the whest environment; normal functions come from math.erf), OPENBLAS_NUM_THREADS=1, seeds fixed. Each check supports the claim it is attached to; none is a claim of its own. C9–C11 were added by the adversarial check; C1, C2, C4, C5 and C8 were re-run then from the script below (C8 at 1e6 samples: E Γ_Z 0.40988 against Var Z 0.40967).

| check | claim | result |
|---|---|---|
| C1 | Dodgson exchange relation; I(a:c\|B) = ½ log(1 + y) (unlocks 33, 46) | three random 6 × 6 SPD matrices: relation residual ≤ 2e-15; CMI equal to 10 digits (0.2632371887, 1.2337629262, 0.1329056055) |
| C2 | chordal (path) junction tree: half the log-det gap = KL to the Markov projection = Σ CMI; zero iff Markov (unlock 46(c)) | generic 6 × 6: 1.227709586 = 1.227709586; tridiagonal precision: −4.4e-16 and 0 |
| C3 | Sheppard and trivariate orthant formulas (unlock 9) | Monte Carlo (8e6 samples) within 1σ in all six cases, e.g. 0.24689 ± 0.00015 vs 0.24702; 0.08717 ± 0.00010 vs 0.08726 |
| C4 | Godsil–Gutman: E_ε det(ε ∘ √A)² = per(A) (unlock 25) | exact enumeration over all sign patterns: n = 3, 0.1946602504 = 0.1946602504; n = 4, 1.1579296180 = 1.1579296180 |
| C5 | arc-cosine iterate from c_0 = 0 (unlock 31) | c_l = 0.318, 0.494, 0.681, 0.834, 0.897, 0.929 (l = 1, 2, 4, 8, 12, 16); 1 − c_16 = 0.0705 against 9π²/512 = 0.1735; l²(1 − c_l) = 43.72 at l = 2000 against 44.41 |
| C6 | Var(wᵀBw) = 8‖B‖_F²/n² (unlock 5(c)) | n = 200, 2e4 samples: 3.899 against 3.947 (Monte Carlo error ≈ 1 %) |
| C7 | E relu(Z) − relu(m) = s[φ(t) − tΦ̄(t)] (unlock 30) | 2e7 samples with a control variate: 8.643e-3 vs 8.643e-3; 7.123e-3 vs 7.138e-3; 1.919e-4 vs 1.911e-4; 3.83e-6 vs 3.57e-6 |
| C8 | centred wall identity (unlock 37), Z = Σ_i v_i relu(⟨x, w_i⟩), d = 3, h = 6, hot unit (μ/s = −0.07) | 4e6 samples: E Γ_Z = 0.40966 vs Var Z = 0.40967 (exact arc-cosine value); smooth form E[(Z−μ) tanh(Z/ε)] = E[tanh'(Z/ε) Γ_Z/ε] at ε = 0.3, 0.1, 0.03: 0.4581/0.4575, 0.4909/0.4898, 0.4954/0.4942; wall form 0.2478 vs 0.2473 (ε = 0.02 kernel); E relu Z = 0.22909 vs μP + wall = 0.22852; uncentred E‖∇Z‖² = 0.550 vs μ² + Var = 0.412 |
| C9 | sign-even and sign-odd quenched parts of a unit's variance at layer 2 (unlocks 5, 22) | exact Cov(a_1) by the arc-cosine kernel, one network per width: rms of Σ_i (W_ik² − 2/n) Var a_i against Σ_{i≠j} W_ik W_jk Cov(a_i, a_j) is 0.117 vs 0.086 at n = 256 and 0.060 vs 0.044 at n = 1024, on an annealed variance of 1.36; through the Gaussian readout, 8.6e-3 vs 6.4e-3 rms per unit at n = 1024 |
| C10 | degree-≤ 2 multigraph sector against det^{−1/2} (unlock 23) | two vertices, F̂ = (0, 0, 1): Gauss–Hermite E[He₂(g_1) He₂(g_2)]/4 = ρ²/2 at ρ = 0.1, 0.3, 0.6 to machine precision. Two vertices with F̂ = (1, 0, d): sector 1 + d²ρ²/2 = 1.045, against the quasi-free det^{−1/2}(I − DR) normalised at ρ = 0, 1.250 (d = 0.5, ρ = 0.6) |
| C11 | softplus temperature expansion (unlock 28) | Z ~ N(0.3, 1), quadrature: (E softplus_T − E relu − (π²/6)T² p(0)) / (E softplus_T − E relu) = −3.9e-2, −1.0e-2, −2.6e-3 at T = 0.2, 0.1, 0.05, i.e. an O(T⁴) remainder |

C8 uses the closed form of the Mehler-smoothed gradient for a two-layer unit. With ∇Z(y) = Σ_i v_i 1[⟨y, w_i⟩ > 0] w_i:

E'[∇Z(u x + √(1−u²) x')] = Σ_i v_i Φ(u ⟨x, w_i⟩ / (√(1−u²) ‖w_i‖)) w_i,

and ∫_0^1 Φ(u s/√(1−u²)) du = ∫_0^{π/2} Φ(s tan θ) cos θ dθ, tabulated on a grid of s.

Compact reproduction script for C1, C2, C4, C5 and C8. Its random instances differ from the table's because the random stream is consumed differently; the identities hold to the same precision.

```python
import numpy as np, itertools, math
rng = np.random.default_rng(1); Phi1 = lambda t: 0.5*math.erfc(-t/math.sqrt(2))
d = lambda M, r, s=None: np.linalg.det(M[np.ix_(r, r if s is None else s)])
# C1  Dodgson and CMI = 1/2 log(1+y)
A = rng.standard_normal((6, 8)); M = A @ A.T; a, B, c = [0], [1, 2, 3, 4], [5]
y = d(M, a+B, B+c)**2 / (d(M, a+B+c)*d(M, B))
print(0.5*np.log(d(M, a+B)*d(M, B+c)/(d(M, B)*d(M, a+B+c))), 0.5*np.log1p(y))
# C2  path junction tree: gap = sum of cut CMIs; zero for a tridiagonal precision
def jt(M):
    k = len(M); ld = lambda r: np.linalg.slogdet(M[np.ix_(r, r)])[1]
    gap = 0.5*(sum(ld([i, i+1]) for i in range(k-1)) - sum(ld([i]) for i in range(1, k-1)) - ld(list(range(k))))
    cmi = sum(0.5*(ld(list(range(t))) + ld([t-1, t]) - ld([t-1]) - ld(list(range(t+1)))) for t in range(2, k))
    return gap, cmi
print(jt(M))
# C4  Godsil-Gutman, exact enumeration
n = 3; A = rng.uniform(0, 1, (n, n)); S = np.sqrt(A)
gg = np.mean([np.linalg.det(np.array(b).reshape(n, n)*S)**2 for b in itertools.product([-1, 1], repeat=n*n)])
per = sum(np.prod([A[i, p[i]] for i in range(n)]) for p in itertools.permutations(range(n))); print(gg, per)
# C5  arc-cosine iterate
f = lambda c: (np.sqrt(1-c*c) + (np.pi - np.arccos(c))*c)/np.pi; c = 0.0
for l in range(16): c = f(c)
print(1 - c, 9*np.pi**2/512)
# C8  centred wall identity for Z = sum_i v_i relu(<x,w_i>)
dd, h = 3, 6; rng = np.random.default_rng(5)
W = rng.standard_normal((dd, h))*np.sqrt(2/dd); v = rng.standard_normal(h)*np.sqrt(2/h) + 0.12
nw = np.linalg.norm(W, axis=0); G = W.T @ W; cs = np.clip(G/np.outer(nw, nw), -1, 1); th = np.arccos(cs)
K = np.outer(nw, nw)/(2*np.pi)*(np.sin(th) + (np.pi - th)*cs); mu = v @ nw/np.sqrt(2*np.pi); var = v @ K @ v - mu**2
xg, wg = np.polynomial.legendre.leggauss(400); tq = (xg+1)*np.pi/4; wq = wg*np.pi/4; sg = np.linspace(-15, 15, 6001)
Ig = np.array([np.sum(wq*np.vectorize(Phi1)(s*np.tan(tq))*np.cos(tq)) for s in sg])
x = rng.standard_normal((4_000_000, dd)); pre = x @ W; g = (pre > 0)*1.0; Z = np.maximum(pre, 0) @ v
gam = np.einsum('si,ik,sk->s', g*v, G, np.interp(pre/nw, sg, Ig)*v); eps = 0.02
print(gam.mean(), var)                                                     # E Gamma = Var Z
print(((Z-mu)*(Z > 0)).mean(), (np.exp(-0.5*(Z/eps)**2)/(eps*np.sqrt(2*np.pi))*gam).mean())  # wall form
```

---

## Appendix A. Claims of the 1 Oct input used here, with their status

| claim of the input | status here | used in |
|---|---|---|
| Heilmann–Lieb real roots; all-matchings growth is local (ACFK) | THEOREM (ACFK Thm 3.3, Thm 1.2; verified) | 15 |
| perfect matchings are not local; one edge in all but a c^n fraction | THEOREM (ACFK Thm 1.7, 1.8; verified) | 15 |
| expansion restores locality (with Gamarnik–Katz); Schrijver's constant at large girth | THEOREM (ACFK Thm 1.9, 1.5; verified) | 14, 15 |
| per_B ≤ per ≤ 2^{n/2} per_B; Bethe = limit over covers | THEOREM (Gurvits; Anari–Rezaei Thm 4; Vontobel Thm 39, as a lim sup; verified) | 14 |
| Kasteleyn (planar); 4^g Pfaffians at genus g | THEOREM (planar case via Lieb–Loss, verified; genus count verified in foundations-input-verified.md via Cimasoni–Reshetikhin) | 24 |
| Pfaffian bipartite graphs (RST, McCuaig) with the Heawood exception | THEOREM (RST verified: the brace characterisation and Little's K_{3,3} criterion) | 24 |
| Kasteleyn's face rule equals Lieb's flux rule | THEOREM for the determinant's modulus on every planar bipartite graph (Lieb–Loss Thm 3.1); for the half-filled energy only on the lattices of Lieb (1994), false in general (Lieb–Loss §VII(B)) | 24 |
| Godsil–Gutman; roots ≤ 2√(d−1); MSS Ramanujan 2-lifts; Kadison–Singer | THEOREM (Godsil–Gutman, KS/paving verified; Ramanujan 2-lifts verified in foundations-input-verified.md, MSS I §5) | 22, 26 |
| "the diagonal is the classical shadow" | ANALOGY: D_0 ⊂ C*(Λ) plays the role of the diagonal masa (unlock 6); paving as a block mean field is SPECULATION (unlock 26) | 6, 26 |
| tropically per = det = assignment | THEOREM as formal tropicalisations; the zero-temperature limit of the determinant falls below the assignment value when top permutations cancel (foundations-input-verified.md, C15) | 32 |
| Chin's proof: tropical c-vectors and T-systems, lifted by positivity and periodicity | THEOREM (Chin 2602.15140 Thm 1.1 verified; the lift is the IIKKN periodicity theorem) | 32 |
| the general permanent is #P-hard; Kasteleyn and JSV cross to finite temperature | THEOREM (from memory); JSV's crossing rests on positivity | 32 |
| cluster variables and T-systems are matching sums; Kuo condensation; urban renewal; Goncharov–Kenyon; Galashin–Pylyavskyy | THEOREM (from memory); context only, no computational use found | 33 |
| MTP₂ distributions are faithful | THEOREM with the graphoid hypothesis (Fallat et al. Thm 6.1, verified) | 34 |
| Idea A: the Dodgson identity, and CMI as a function of its two terms | DERIVED and checked (C1); sharpened to I = ½ log(1 + y) | 33, 46 |
| Idea A: a cluster structure on Gaussian states with Markov networks as boundary strata | SPECULATION here; partly DERIVED in foundations-input-verified.md §2 (cube moves change separators; graphical models as coordinate strata on chordal graphs with ≤ 6 vertices); its CONJECTURE A4* is open, and MTP₂ is not the positive part | 33 |
| Idea B: separator vs tree local-to-global, with pseudorandomness as the switch | DERIVED for the network as unlock 47 (its quasi-free reading is a dictionary choice, unlock 47 guard); the tree half needs correlation decay; the quantum half is a CONJECTURE (B-programme 8.4) | 47, 49 |
| Idea C: a provably periodic quantum circuit | construction THEOREM-level per foundations-input-verified.md §4; not used | — |

---

## Sources

**Opened and checked in this session (the statements used here were read in the source).**
- N. Anari, A. Rezaei, *A tight analysis of Bethe approximation for permanent*, arXiv:1811.02933. Thm 4: per ≤ √2^n bethe. Gurvits's per ≥ bethe as quoted (Thm 3); the definition of bethe(A); tightness at I_{n/2} ⊗ J_2.
- M. Abért, P. Csikvári, P. Frenkel, G. Kun, *Matchings in Benjamini–Schramm convergent graph sequences*, arXiv:1405.3271. Thms 1.2, 1.5, 1.7, 1.8, 1.9, 3.3, Remark 3.6, Thm 4.1.
- P. O. Vontobel, *The Bethe permanent of a non-negative matrix*, arXiv:1107.4196. Thm 39 (the cover limit), Thm 32, Lemma 48.
- L. Zhang, G. Naitzat, L.-H. Lim, *Tropical geometry of deep neural networks*, arXiv:1805.07091. Thm 5.2, Prop 5.6, Thm 6.3.
- A. Marcus, D. Spielman, N. Srivastava, *Interlacing families II*, arXiv:1306.3969. KS_2, Cor 1.5, Thm 6.1 (paving).
- S. Fallat, S. Lauritzen, K. Sadeghi, C. Uhler, N. Wermuth, P. Zwiernik, *Total positivity in Markov structures*, arXiv:1510.01290. Thms 5.3, 6.1, Example 6.2.
- E. H. Lieb, M. Loss, *Fluxes, Laplacians and Kasteleyn's theorem*, Duke Math. J. 71 (1993) 337–363, arXiv:cond-mat/9209031.
- E. H. Lieb, *Flux phase of the half-filled band*, Phys. Rev. Lett. 73 (1994) 2158.
- I. Chin, *Half-periodicity of Zamolodchikov periodic cluster algebras*, arXiv:2602.15140. Thm 1.1 and the structure of its proof.
- S. Chien, L. Rasmussen, A. Sinclair, *Clifford algebras and approximating the permanent*. The Godsil–Gutman estimator and the critical ratios quoted in unlock 25.
- C. Moore, A. Russell, *Approximating the permanent via nonabelian determinants*, arXiv:0906.1702. §1: the critical ratios of Karmarkar et al. and Chien–Rasmussen–Sinclair, and the absence of a polynomial-time algorithm for the Clifford determinants beyond the quaternions (opened by the adversarial check).
- J. M. Mooij, H. J. Kappen, *Sufficient conditions for convergence of the sum-product algorithm*, arXiv:cs/0504030.
- N. Robertson, P. D. Seymour, R. Thomas, *Permanents, Pfaffian orientations, and even directed circuits*, Ann. of Math. 150 (1999) 929–975 (with Little's theorem as stated there).

**Re-opened by the adversarial check (1 Oct).** ACFK (Thms 1.2, 1.5, 1.7–1.9, 3.3, 3.5, Remark 3.6, Thm 4.1), Anari–Rezaei (Thms 3–4, §1.2, the tight example), Vontobel (Def. 38, Thm 39), Zhang–Naitzat–Lim (Thm 5.2, Props 5.5–5.6, Thm 6.3 with its hypotheses), MSS II (Cor 1.5, Thm 6.1, the r ≥ 1/ε² remark), Fallat et al. (Thms 5.3, 6.1, Example 6.2, §4.1), Chin (Thm 1.1, Thm 2.14 = the IIKKN periodicity theorem), Mooij–Kappen (Cor. 3, §V, Table I), RST (main theorem, Little's theorem, Algorithm 9.7 in O(n³)), Lieb–Loss (Thm 3.1, §VII(B), Appendix).

**Through the programme's digests (cited by digest section).** Chen–Rouzé arXiv:2504.02208; Yang arXiv:2609.38007; Hanin–Nica arXiv:1812.05994; ALO, CLV, Chen–Eldan, Alev–Lau, Oppenheim, Dinur–Kaufman (hdx digest); EKZ, Hoory–Linial–Wigderson, Willett–Yu (expanders digest); Cipriani, Goldstein–Lindsay, Vernooij–Wirth, Carlen–Maas (nc-Dirichlet digest); Hammersley–Clifford, Brown–Poulin, Leifer–Poulin (hammersley-clifford digest); Fawzi–Renner and JRSWW (Chen–Rouzé digest).

**From memory (standard; not re-opened).** Isserlis 1918; Sheppard 1899; Plackett 1954; Price 1958; Kahane 1986; Cover 1965; Wendel 1962; Kemeny–Snell 1960; Sudakov 1978; Diaconis–Freedman 1984; Mingo–Speicher 2006; Weitz 2006; Mossel–Sly 2013; Kesten–Stigum 1966; Godsil–Gutman 1981; Godsil 1981; Marcus–Spielman–Srivastava, Interlacing Families I; Barvinok 2016; Patel–Regts 2017; Borcea–Brändén–Liggett 2009; Valiant 1979; Jerrum–Sinclair–Vigoda 2004; Galluccio–Loebl 1999; Tesler 2000; Füredi–Komlós 1981; Bai–Yin 1988; Bolthausen 2014; Aronszajn 1950; Kayalar–Weinert 1988; Hubbard 1959 and Stratonovich 1957; Cho–Saul 2009; Nourdin–Peccati (Stein kernel and Malliavin calculus, 2009–2012); Houdré–Pérez-Abreu 1995; Maslov dequantisation.

---

## Critique log (adversarial check, 1 Oct 2026)

Every claim was re-derived or re-computed; the sources marked as verified were re-opened (Sources). Numbers checked in the scratchpad are C1–C11. Corrected:

1. **Unlock 3, §9 I, unlock 36, §11 heisenberg (scale of the closure error).** The note used the BRIEF's superseded 4e-5 for the Gaussian closure at n = 1024 and wrote "0.2 n^{−1/2} rms, next order right to 2 %, ĝ to 2–3 %". Measured: raw 4.3e-6, ≈ 2.1e-3 rms per unit, scaling ≈ 2/n for 128 ≤ n ≤ 1024. The required relative accuracy is ≈ 6 %.
2. **Unlocks 5, 22, §8 row 22, §9 II (sign sectors).** "The quenched per-unit content is sign-odd" and "sign average = annealed or tree object" were false. The sign average keeps the quenched magnitudes |W|, whose O(n^{−1/2}) contribution to each unit's variance is as large as the sign-odd edge term (C9). The sign average, the annealed average and the Bethe lift are three different references. The averaged covariance map is twice the conditional expectation onto the scalars.
3. **Unlock 23 (determinantal sector).** "The degree-≤ 2 multigraph sector equals det^{−1/2}(I − DR)" was false; it repeated the signings design's P3(i). The determinant sums vertex-revisiting walks; the hard-core sector is an α-permanent (C10).
4. **Unlock 4 guard, §9 I.** E‖∇z‖² is not ≈ μ² + s² at depth: it is 6.4 (μ² + s²) at layer 16 (tropical R-E2). The excess is the chaos over-count of unlock 37, not only the mean spike.
5. **Unlock 46(c), C2.** The log-det gap is twice the KL divergence and the CMI sum; the checked number was already half the gap.
6. **Unlock 6.** N(σ) counts forward paths starting at σ (Note 1 §3.1), not histories ending at σ.
7. **Unlock 2.** The barycentre vector state of Note 2 and the commutative atom-conditioned ω_σ were identified; for the latter, p_σ(1 − e_u) = 0.
8. **Unlock 40(d).** "U1 as a theorem" reversed the B-programme's own relabelling to ANALOGY.
9. **Unlock 42.** P^k needs an absolute gap; ghost projections are the gap mechanism's instance for expanders, not an exception to it.
10. **Unlock 11.** Mooij–Kappen's spectral condition does not dominate Dobrushin's (their Table I).
11. **Unlock 26.** Gurvits's per ≥ per_B goes through Schrijver's inequality, not capacity; MSS paving is vacuous at n = 1024 (r ≥ 1296).
12. **Smaller.** Unlock 8 (√μ_1, not μ_1); unlock 9 (non-centred layers need Owen's T or a quadrature; ZNL hypotheses); unlock 12 ("usable certificate" contradicted unlock 17); unlock 14 (the lift's single-site laws are not the network's); unlock 15 (zero-freeness is sufficient, not "exactly when"; activity 1 vs any activity); unlock 16 (Mossel–Sly is ferromagnetic); unlock 17 (Ψ has zero diagonal; the surrogate's 15–29 is at n = 128); unlock 18 (undefined label KNOWN-LINK; the negative result is stated for its dictionary); unlock 22(b) (the path tree is not the universal cover); unlock 24 (Lieb's energy result is lattice-specific); unlock 25 (Clifford determinants are not known to be computable efficiently); unlock 27 (Arora's theorem needs biases); unlock 28 (O(T⁴) needs smoothness at 0); unlock 30 (7.1e-6 at t = 4; t = 2.054; "sixty times" holds for s = 1); unlock 32 (formal tropicalisation only); unlock 33 (singletons; blocks need canonical correlations); unlock 39 (Connes distance = W_1 as the digest has it; one layer at gain 1.07); unlock 43 (quantum HC fails with triangles); unlock 44 (the 40 % is the (2,1) third-order slice; strong vs weak lumpability); unlock 46 guard (copula invariance needs a copula); unlock 48 ("at most" the geometric mean; E_{A∪B}).

Bias, flagged in place (unlock 3, unlock 36, unlock 47 guard, §9 III and V). The note's design rule, "glue along depth through quasi-free separators and decorate across width", reads as covariance propagation plus cumulant hubs and sources. BRIEF rule 1 admits these existing estimators only as baselines, and all six fresh-slate designs converged on them (CONVERGENCE.md). The Duhamel identity of unlock 36 was quantified only with the Gaussian closure as Π. The EscAI oracle, a cumulant-ladder measurement, was called "the precise content" of an exact identity. Negative results measured in one dictionary (old cumulant slices in propagator bases) were stated for "signed content" in general.
