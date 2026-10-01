# Foundations, part 2: the theoretical unlocks

*Fresh slate, 1 Oct 2026. This note distils the programme's theoretical unlocks into a numbered list. For each unlock it gives the precise statement, its status, and its computational meaning for an estimator of quenched per-neuron means. Contract: [BRIEF.md](BRIEF.md). Part 1 (`foundations.md`, when it appears) verifies the user's 1 Oct input; Appendix A lists the claims of that input used here, with their status.*

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

In Note 2's form, excluding a vertex u is Lüders conditioning: ω_σ^{1−e_u} = ω_{σ∖u}.

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
- *Exactly:* the means need no approximation beyond two scalars per unit per layer: the face mass P(z_{l,j} > 0) and the wall term p(0) τ(0).
- The error of any closure decomposes exactly as μ × (face-mass defect) + (wall defect).
- The question the state must answer at layer l is n one-dimensional laws along fresh random directions, not the joint law. This is the precise content of the EscAI oracle (per-neuron marginal laws suffice; BRIEF §3).
- *With controlled error:* for a fixed state and random lines the projections are Gaussian at leading order, under thin-shell and weak-correlation conditions on the centred state (random-projection central limit theorems: Sudakov; Diaconis–Freedman; from memory). The Gaussian closure is that leading order.
- Gaussian closure leaves ≈ 0.2 n^{−1/2} rms per unit (raw 4e-5 at n = 1024). The bar is ≈ 0.004 n^{−1/2}, so the next order must be right to about 2 %.

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

**Status.** DERIVED. Found independently by the tropical design (I2, checked at width 8, depth 3: R-E0) and the faces design (E2). Parts (d) and the history-sum derivation of (b) are added here.

**Computational meaning.**
- *Exactly:* the mean of every neuron is a sum, over upstream units, of wall integrals. Each is (Gaussian density of the unit's pre-activation at 0) × (wall-conditional expectation of slope² × cavity response).
- Interiors of faces contribute nothing; only codimension-1 faces carry the mean, and frozen units contribute little.
- It is the exact meet-in-the-middle split: the wall measure is a forward (state) object, the cavity response a backward (question) object.

**Guard.** The identity is exact, but its wall-by-wall factorisation is not perturbative at He initialisation.
- The weight ‖∇_x z‖² is uncentred: its mean is ≈ μ² + s², not s². The mean spike sits in the wall weights.
- The own-wall and inherited terms each miss their factorised forms by O(1), and the two misses cancel (tropical R-E2, correlation −1.00).
- The centred form, unlock 37, removes this.

### 5. Fresh randomness: annealed twirl, quenched single Kraus operator, sign sectors ⊕ [B, M, H, K]

**Statement.** Fix the law of a_l and let W = W_{l+1} be independent of it.
- (a) *Annealed.* E_W[Wᵀ S W] = (2/n) Tr(S) I for every symmetric S. The averaged second-order map is the conditional expectation onto the scalars: a twirl, with all non-trivial eigenvalues 0. Also E[W^{⊗odd}] = 0.
- (b) *Quenched.* S ↦ Wᵀ S W is completely positive with a single Kraus operator.
- (c) For w ~ N(0, (2/n) I) and symmetric B: Var(wᵀ B w) = 8‖B‖_F²/n². This is a relative fluctuation √(2/PR(B)) around (2/n) Tr B, with PR(B) = (Tr B)²/‖B‖_F².
  - Bulk objects (PR ∝ n) fluctuate per unit at O(n^{−1/2}), and their layer averages at O(1/n).
  - A spike (PR = O(1), e.g. the mean direction) fluctuates per unit at O(1).
- (d) *Sign sectors.* The group Z_2^{n×n} of sign flips of W preserves its law.
  - Averaging a quenched quantity over it is the conditional expectation onto functions of |W|.
  - For a polynomial in W, the average keeps exactly the monomials in which every weight occurs to an even power: the paired, tree sector.
  - The remainder is the sign-odd, quenched sector.

**Status.**
- (a)–(c): THEOREM (Wick/Isserlis; elementary). (a) and (b) are also derived in transfer-spectrum-measurement (finding 5 and §2.5). (c) checked (C6).
- (d): DERIVED, since E_S ∏ S_e^{k_e} = 1 iff every k_e is even. It is the network form of Godsil–Gutman (unlock 22).

**Computational meaning.**
- *Exactly and cheaply:* every annealed prediction is a scalar recursion, and self-averaging bulk quantities (normalised traces, layer averages, spectra) are annealed-computable to O(1/n) (second-order freeness, Mingo–Speicher, from memory).
- *Warning:* the scored quantity is per unit, and there the quenched content is sign-odd.
  - From layer 2 on, each unit's mean carries an O(1) sign-odd term: its alignment ⟨m_{l−1}, w_k⟩ with the mean direction. The exact recursion of unlock 3(a) handles it at no cost.
  - Beyond it, the O(n^{−1/2}) sign-odd terms pair unpaired fresh weights with the state's joint structure.
- The sign average removes both. The annealed or tree object can therefore never be the estimator; it is the reference against which the quenched sector is measured.

---

## 2. Faces and barycentres [F]

### 6. Three resolutions and the arrow algebra of the face quiver [F, K, H]

**Statement (Note 1 §2–3; Note 2 §1–3).**
- The layered quiver Λ of faces carries a Toeplitz–Cuntz–Krieger family: partial isometries s_γ, one per arrow, with s_γ* s_γ = p_{source} and Σ_{γ into τ} s_γ s_γ* ≤ p_τ.
- For a finite layered quiver, T(Λ) ≅ ⊕_σ M_{N(σ)} and C*(Λ) ≅ ⊕_{σ_0} M_{N(σ_0)}, where N(σ) counts the histories ending at σ. The Bratteli heights satisfy h_{l+1} = F_l h_l, with F_l the incidence matrix.
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
- The loss is bounded by a local quantity, four arrows at a time, divided by μ_1 = n². The Fisher-information loss is exactly E Var(F | ends).

**Guard.** The unlock requires the Gibbs form, i.e. positivity. Dictionary v1 fails it: its fitted F has residual 0.42–0.5 against coboundaries, and its face process has 25–35 % memory (mlp-bridge). The failure is charged to the dictionary.

### 9. Low-codimension face masses are elementary; high-codimension ones are not [F, H, M]

**Statement.** Let (z_1, …, z_k) be a centred Gaussian with correlations ρ_ij.
- k = 2 (Sheppard): P(z_1 > 0, z_2 > 0) = 1/4 + arcsin ρ_12/(2π).
- k = 3: P(z_1, z_2, z_3 > 0) = 1/8 + (arcsin ρ_12 + arcsin ρ_13 + arcsin ρ_23)/(4π).
- k ≥ 4: the orthant mass is a Schläfli-type function with no elementary closed form. With non-zero means, k = 2 needs Owen's T and larger k needs numerical integration (Genz).
- Plackett's reduction: ∂P(all > 0)/∂ρ_ij = φ_2(0, 0; ρ_ij) P(rest > 0 | z_i = z_j = 0), a codimension-2 face mass. So every orthant mass is an integral of lower-dimensional ones along a correlation path.
- Counting:
  - A central arrangement of N hyperplanes in general position in R^d has 2 Σ_{i<d} C(N−1, i) regions (Cover; Wendel; Schläfli). So layer 1 alone (N = d = n) has all 2^n orthants.
  - Zhang–Naitzat–Lim Thm 6.3 bounds the number of linear regions of the whole network by ∏_l Σ_{i ≤ d} C(n_l, i).

**Status.** THEOREM (Sheppard; the trivariate formula; Plackett; Cover; from memory, except that the two formulas were checked in C3 and Zhang–Naitzat–Lim was verified).

**Computational meaning.**
- *Exactly and cheaply:* masses and barycentres of faces of codimension ≤ 3 in any Gaussian layer cost O(1) elementary functions each. This includes layer 1, which is exactly Gaussian. Pair and triple face data of a Gaussian layer therefore cost O(n²) and O(n³) elementwise operations.
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
  - if ρ(A) < 1, where A_{i→j, k→l} = tanh|J_ij| δ_{il} 1[k ∈ ∂i ∖ j] is the weighted non-backtracking operator on directed edges (Mooij–Kappen, which improves on Dobrushin);
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
- So the measured spectral independence η ≈ 2–8 of the face law is a usable certificate. A Dobrushin-type influence sum over the n/2 hot units is O(√n).

### 13. Power counting for a fresh layer: loops of fresh weights are suppressed, loops through the state are not ⊕ [B, M, K]

**Statement.** Fix the law ω of a_l, and let w = w_k be a fresh column. Organise the deviation of the law of ⟨a_l, w⟩ from its Gaussian approximation by two things: which units of layer l the fresh weights touch, and the dependences of ω that join those units.
- (a) Sign-paired (even) contributions are of the order their counting gives. Each sign-unpaired weight turns a sum over units into a random-sign sum, which costs n^{−1/2} unless the summand is coherent in sign. The mean direction is the one coherent case: one unpaired weight, O(1).
- (b) Contributions whose dependence graph through ω contains a cycle are suppressed by an extra n^{−1/2} per independent cycle, relative to tree-shaped ones.
- (c) Hence, at the needed precision, only tree-shaped (star, path) contractions of ω's dependences enter the per-unit laws at layer l + 1. But tree-shaped contractions of ω's pairwise dependences already contribute at the same order as the single-site terms. Node beliefs alone are wrong at the leading quenched order; edge beliefs are necessary.

**Status.**
- DERIVED. This is Wick counting for the fresh layer, done independently in the signings design (§1, table) and the bethe design (§2, table). It is consistent with the measured n^{−0.8} residual of first-order gate diagrams (BRIEF §3).
- CONJECTURE: (b) holds uniformly in depth for the states the network actually produces. The signings design's test T1 measures it.

**Computational meaning.**
- *Cheaply:* the readout of layer l + 1 from a state at layer l needs only tree-shaped contractions. These are matrix products of the fresh layer with node and edge data of the state; no loop sums are needed.
- All quenched information beyond the mean recursion lives in the state's joint structure: its edges, and how they were built through depth.
- That is where loops matter and where unlocks 43–47 apply.

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
- *Exactly:* the tree (Bethe) reference of the quenched network is a well-defined object. It is computable by n² one-dimensional convolutions per layer (characteristic functions on a grid). It has exact non-Gaussian single-site laws, all the quenched weights, and no joint structure.
- The difference between the quenched network and its lift is exactly the loop content.
- (a) also says that for permanent-type sums the Bethe value is accurate only to e^{O(n)}, i.e. to O(1) per site. The Bethe reference is a per-site object, never a precision tool on its own, consistent with unlock 13(c).

### 15. Locality comes from zero-freeness; hard constraints are non-local unless the geometry expands [B, M, T]

**Statement.**
- (a) Heilmann–Lieb: the matching polynomial of a graph of maximum degree d has only real roots, all of modulus ≤ 2√(d−1).
  - The power sums of its roots count closed tree-like walks, i.e. walks in Godsil's path tree (ACFK Remark 3.6).
  - On large-girth sequences the root distribution tends to Kesten–McKay (ACFK Thm 4.1).
- (b) The matching entropy per vertex, i.e. the monomer–dimer free energy at any positive activity, is estimable: it is continuous under Benjamini–Schramm convergence (ACFK Thm 1.2).
  - The number of perfect matchings is not estimable, even for d-regular bipartite graphs (Thm 1.8).
  - In some such graphs one edge lies in all but a c^n fraction of the perfect matchings (Thm 1.7).
- (c) On bipartite δ-expanders every edge has p(e) ≥ (1/d) n^{−2 ln(d−1)/ln(1+δ)} (Thm 1.9). With Gamarnik–Katz, the growth of perfect matchings becomes local.
- (d) General principle (Barvinok; Patel–Regts; from memory). Suppose a partition function has no zeros in a neighbourhood of the path from 0 to the target activity. Then its logarithm is determined to ε by O(log(N/ε)) Taylor coefficients, each a sum over connected local structures.

**Status.** THEOREM: (a)–(c) verified in ACFK; (d) from memory.

**Computational meaning.**
- An estimator may trust a local (bounded-radius) computation exactly when the relevant generating function is zero-free along its interpolation path. The order (equivalently, the radius of the local structures) needed is O(log(N/ε)/log(R/|λ|)), where |λ| is the target activity and R the distance to the nearest zero, after a conformal map of the zero-free region to a disk.
- Hard (zero-temperature) constraints break locality unless the geometry expands.
- *Network reading:* the linear-to-ReLU homotopy relu_γ has its singularities at γ = (1 ± i)/2, which gives the geometric rate 0.71 per order at γ = ½ (tropical design P6). A conformal re-parametrisation of the path that avoids these points is the Barvinok-type repair. It was not tried (SPECULATION).

### 16. Uniqueness, reconstruction, and why expansion hurts separators [B, K]

**Statement.**
- On the d-regular tree (branching number d − 1), Gibbs uniqueness for the zero-field ferromagnetic Ising model holds iff (d−1) tanh β ≤ 1, an ℓ¹ condition.
- Reconstruction (Kesten–Stigum) occurs when (d−1) tanh² β > 1, an ℓ² condition on the second eigenvalue.
- Mossel–Sly: on any graph of maximum degree d, Glauber dynamics mixes in O(n log n) when (d−1) tanh β < 1.
- Yang's CMI bound has a boundary gain that vanishes on expanders (Yang digest, B1). Separator bounds are therefore weakest exactly where tree bounds are strongest.

**Status.** THEOREM (as quoted in the expanders and hdx digests; Kesten–Stigum and Mossel–Sly also from memory).

**Computational meaning.**
- Between the ℓ² and ℓ¹ thresholds, local tree computations still work on average, through spectral criteria, although worst-case influence sums do not. Dense pseudorandom layers live in this window (unlock 12).
- At low temperature (frozen faces with strong couplings), expansion creates rigidity (unlock 15(b)), and neither tree gluing nor separator gluing is local.

### 17. Spectral independence of the face law: what it gives and what it does not [B, K, F]

**Statement.** Let ω be a law on {0,1}^V (the face law, restricted to the non-deterministic vertices), with D = diag Var(e_u).
- The influence matrix is Ψ(u, v) = P(v | u) − P(v | ū) = Cov_uv/Var_u. Then λ_max(Ψ) ≤ η_0 iff Cov ⪯ (1 + η_0) D.
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
  - a Gaussian surrogate with the same correlations gives η_0 = 15–29.

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

**Status.** DERIVED + KNOWN-LINK + measured (transfer-spectrum-measurement §2.2–2.6, §3.4, B2, B6; Hanin–Nica for the polymer reading). That k_ε/n is constant is a CONJECTURE beyond the generic sources tested.

**Computational meaning.**
- *Negative and decisive:* no fixed-rank compression of signed content carried across depth works. The rank needed is a constant fraction of n, ≈ 0.3 n at the brief's accuracy.
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
  - Correlations among the uncertain gates fall as the mean shift freezes units: η_0 = 2–8 measured, against 15–29 for the zero-mean surrogate.

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
  - Its degree-≤ 2 part (paths and cycles) is a determinant (unlock 23).
- Bosonic sums (hafnians, permanents) are #P-hard in general (Valiant, from memory).
- An estimator must therefore use the structure of the sum: trees (unlock 13) and determinantal sectors (unlock 23). It must never evaluate the full hafnian.

### 22. Sign-averaging is the tree projection; a quenched network is one signing ⊕ [M, B, T]

**Statement.**
- (a) Godsil–Gutman. For a graph G with adjacency matrix A and a uniformly random signing s, E_s det(x − A_s) = μ(G, x), the matching polynomial. The cycle terms of the determinant expansion carry an odd power of some edge sign and vanish on average; the matchings survive.
- (b) The moments of the matching polynomial's roots count closed tree-like walks, i.e. walks in Godsil's path tree, a universal-cover object (ACFK Remark 3.6).
  - A signing is a 2-lift, i.e. a Z/2 gauge field.
  - Interlacing families give a signing whose largest new eigenvalue is ≤ 2√(d−1); for bipartite graphs this yields Ramanujan 2-lifts (Marcus–Spielman–Srivastava, Interlacing Families I; from memory).
- (c) Network form (unlock 5(d)). Averaging any quenched quantity over the signs of a fresh layer gives its paired (tree) sector. From layer 2 on, the per-unit mean has two sign-odd parts:
  - an O(1) part, its alignment with the mean direction;
  - an O(n^{−1/2}) part, from unpaired weights joined by the state's joint structure.

**Status.** (a), (b) THEOREM (Godsil–Gutman; Godsil; ACFK; MSS). (c) DERIVED.

**Computational meaning.**
- *Exactly:* the decomposition quenched = (sign average) + (sign-odd sectors) is canonical and computable sector by sector.
- The sign average is the annealed or tree object: cheap, and wrong per unit at O(1). The scored information is entirely sign-odd.
- A design built on a tree, Bethe or annealed object must therefore carry the sign-odd sectors explicitly:
  - the coherent one through the exact mean recursion (unlock 3(a));
  - the incoherent one through the state's edges (unlock 13).

### 23. The free (determinantal) sector of an interacting Gaussian sum [M, H, K]

**Statement.** In the multigraph expansion of unlock 21(b), take the sub-sum over multigraphs of maximum degree ≤ 2: disjoint unions of paths and cycles. It equals a Gaussian integral of the exponential of a quadratic form, summed exactly by det^{−1/2}(I − D R) and a resolvent (free bosons; with fermionic signs, det^{+1}) (signings design P3(i)).

**Status.** THEOREM (Gaussian integral identity, as stated in the signings design).

**Computational meaning.**
- *Exactly and cheaply:* the quasi-free part of any sum over pair interactions is one determinant and one resolvent, O(n³).
- Only vertices of degree ≥ 3, the genuinely interacting hubs, need local (tree) treatment.
- This is the matching-side form of the separator rule of unlock 47: quasi-free parts are glued by linear algebra, interacting parts by trees.

### 24. Determinantal evaluation needs Pfaffian orientations; dense layers have none [M, T]

**Statement.**
- Kasteleyn: a planar graph has an orientation that turns its perfect-matching count into a Pfaffian. Genus g needs 4^g Pfaffians (Galluccio–Loebl; Tesler; from memory).
- Little: a bipartite graph is Pfaffian iff it contains no even subdivision of K_{3,3} as a central subgraph.
- Robertson–Seymour–Thomas (and McCuaig): a brace is Pfaffian iff it is the Heawood graph or is built from planar braces by repeated 4-sums. This is recognisable in O(n³).
- Lieb–Loss: on planar bipartite graphs, the phases that maximise |det| are determined completely, which gives another proof of Kasteleyn.
- Lieb: the optimal flux for the half-filled band is π per square plaquette and 0 per hexagon. This matches Kasteleyn's face rule.

**Status.** THEOREM (Little, RST, Lieb–Loss and Lieb verified; the genus count from memory).

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
  - bounded by a constant for Clifford-algebra signs of growing dimension (Chien–Rasmussen–Sinclair).
- On random instances the ratio is small: polynomial for random 0/1 matrices (Frieze–Jerrum), and subexponential w.h.p. (Costello–Vu).

**Status.** THEOREM (verified through Chien–Rasmussen–Sinclair). Unbiasedness checked for n = 3, 4 by exact enumeration (C4).

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
- (b) Gurvits: for real-stable generating polynomials, capacity bounds give per ≥ n!/n^n for doubly stochastic matrices (van der Waerden), and per ≥ per_B. Capacity is computed by Sinkhorn scaling.
- (c) Real stability of a multi-affine generating polynomial implies negative dependence: it is strongly Rayleigh (Borcea–Brändén–Liggett, from memory).

**Status.** THEOREM ((a) MSS II verified; (b) through Anari–Rezaei; (c) from memory).

**Computational meaning.**
- (a) is SPECULATION for the network. Pave the zero-diagonal matrix T = D^{−1/2} Cov(g) D^{−1/2} − I of the face law. Its top eigenvalue is η_0, and ‖T‖ ≤ max(η_0, 1).
  - Paving splits the units of a layer into r = (6‖T‖/ε)^4 blocks, inside each of which the face law has spectral independence ≤ ε. Random partitions do comparably well on random-like T.
  - The face law would then be near-product inside blocks and coupled only across blocks: a certified block mean field.

**Guard for (b) and (c).** The network's face law has positive correlations, so it is not strongly Rayleigh. Real-stability tools (capacity, SLC exchange walks, the U2 conjecture of local-to-global-unlocks) do not apply to it directly. They apply to sign-averaged (matching) objects.

---

## 5. Tropical skeleton and temperature [T]

### 27. ReLU networks are tropical rational maps [T, F]

**Statement (Zhang–Naitzat–Lim).**
- With integer weights (rational weights after scaling), every bias-free ReLU network is a tropical rational map: a difference of two tropical polynomials (Thm 5.2). With real weights it is a tropical rational signomial map (Prop 5.6).
- Every neuron is z = P − Q, with P and Q convex, positively homogeneous and piecewise linear: support functions of polytopes (Newton polytopes).
  - A linear layer acts by Minkowski sums, with W = W⁺ − W⁻.
  - A ReLU acts by a convex hull: relu(P − Q) = max(P, Q) − Q.
- Linear regions are the cones of the common refinement of the normal fans. They correspond to vertices of upper faces.
- The number of regions is at most ∏_l Σ_{i ≤ d} C(n_l, i) (Thm 6.3).
- ReLU networks are exactly the continuous piecewise-linear functions (Arora et al.).

**Status.** THEOREM (verified in Zhang–Naitzat–Lim; the P − Q recursion as in the tropical design §1).

**Computational meaning.**
- The faces of unlock 1 are the cells of a tropical fan, and the walls of unlock 4 are its tropical hypersurface.
- At zero temperature, the per-neuron mean is a difference of Gaussian mean widths, E a = w_G(conv(A ∪ B)) − w_G(B). Unlock 29 shows that this difference cancels catastrophically.

### 28. Temperature: Maslov dequantisation of the ReLU [T, H]

**Statement.**
- max(u, v) = lim_{T→0} T log(e^{u/T} + e^{v/T}), so relu = lim_{T→0} softplus_T.
- Each tropical polynomial becomes a subtraction-free exponential sum Z_T(x) = Σ_v e^{⟨v, x⟩/T} over the vertices of its Newton polytope. Minkowski sums become products and convex hulls become sums.
- For a density p continuous at 0: E softplus_T(Z) = E relu(Z) + (π²/6) T² p(0) + O(T⁴).
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
- *Exactly, at the bar:* a unit can be replaced by its linear arrow only if t ≳ 4, where the per-unit error is ≤ 7e-6·s.
  - At t = 3 the error is 3.8e-4·s.
  - At t = 2.05 (gate certainty 0.98) it is 7.3e-3·s, sixty times the bar's rms.
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
- The tropicalisations of the permanent and the determinant coincide: both become the assignment problem, which is solvable in polynomial time.
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
- Take a Gaussian with covariance M on a ⊔ B ⊔ c. The Desnanot–Jacobi (Dodgson) relation in subtraction-free form reads

  det M[aB] det M[Bc] = det M[aBc] det M[B] + (det M[aB|Bc])².
- Put y = (det M[aB|Bc])²/(det M[aBc] det M[B]). Then I(a : c | B) = ½ log(1 + y).
- Markov (a ⊥ c | B) is y = 0: the degeneration of the exchange relation to one monomial.
- Under the tropical limit (log scale, max-plus), ½ log(1 + y) → ½ max(0, log y).

**Status.**
- DERIVED from Desnanot–Jacobi (checked, C1).
- Identifying y with a cluster y-variable of an octahedron-recurrence seed is the user's Idea A. As an organising structure for Markov networks it is SPECULATION.

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
- Forward compression of signed content is impossible (unlock 18). Backward evaluation at full rank n costs only matrix products: the heisenberg design pulls slices back to the source layer at 4n³ per (source, target) pair, with no intermediate n³ tensor.

### 36. Duhamel telescoping and the bilinear (doubly robust) error ⊕ [H, B, K, F]

**Statement.**
- Let ν_l be any reference chain with ν_1 = μ_1 and ν_{l+1} = Π(T_l ν_l). Here T_l is the exact push-forward by arrow l and Π is any projection onto a tractable family, for example the Gaussian law with the same mean and covariance (the Gaussian closure, used here only as the reference). Put ρ_l = T_{l−1} ν_{l−1}. Then, exactly,

  E_{μ_L} r_j − E_{ν_L} r_j = Σ_{l=2}^{L} (ρ_l − ν_l)[g_{l,j}].
- For any model ĝ_l of the pulled-back question,

  truth − (reference + Σ_l (ρ_l − ν_l)[ĝ_l]) = Σ_l (ρ_l − ν_l)[g_l − ĝ_l].

**Status.** DERIVED (heisenberg design, Theorems 1–2; a one-line telescoping). It is the commutative analogue of Chen–Rouzé's telescoping against the erased state (A-unification §5).

**Computational meaning.**
- *With a controlled error that is a product:* the error of reference plus correction is bilinear, (local defect of one arrow) × (model error of the question).
- A cheap forward reference and a cheap backward model therefore give an error of the order of the product of their errors. Example: the defects are O(n^{−1/2}) per unit and computed exactly from the reference, so a question model with ≈ 2–3 % relative error meets the bar, which needs the first-order correction to ≈ 2 %.
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
  - In the commutative case, the associated Connes distance between states is the Kantorovich–Rubinstein W_1 distance for the intrinsic metric (as read in the nc-Dirichlet digest).
- (b) Hence |ω(q) − ω̃(q)| ≤ Lip(q) · W_1(ω, ω̃) for every observable q.
  - For a pulled-back question, Lip(g_{l,j}) is bounded by the norms of the gated propagators from l to L.
  - Measured gains (transfer-spectrum §2.6): a generic direction has mean-square gain 2E[Φ²] ≈ 0.58–0.95 per layer, a contraction. The mean direction is transported with 3.1–7.8 times the bulk gain after 8–15 layers (n ≥ 256).

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
- (d) The restriction of the KMS state to the first t layers depends on the terminal condition at most through ∏_{l ≥ t} tanh(Δ_l/4) in total variation (R7; U1 as a theorem).

**Status.** DERIVED (B-programme Lemmas 6.1–6.3, Thm 6.4, Cor 6.5; checked there).

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
  - More crudely, by P^k with k = O(λ^{−1} log(1/ε)).
- In coarse geometry: with a gap, the global projection is a norm limit of finite-propagation operators (Roe algebra; Kazhdan projections). Ghost projections are the exception (Willett–Yu, as read in the expanders digest).

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
- The quantum version holds for commuting Gibbs states (Brown–Poulin; Leifer–Poulin).
- Intersection holds quantitatively, with Friedrichs-angle constant 1/(1 − c_F).

**Status.** THEOREM (Hammersley–Clifford; Möbius) + DERIVED (the pointwise approximation and quantitative intersection, in the hammersley-clifford digest and A-unification).

**Computational meaning.**
- *Exactly:* whether a carried state may be represented by local potentials is decided by local mixed differences.
- An approximately Markov state has approximately clique-supported potentials, with an explicit constant.
- Pinning is free and marginalising costs fill-in, so a design should condition (pin) rather than marginalise wherever it can.

### 44. Depth is exactly Markov at full resolution; memory is fill-in ⊕ [K, F, B]

**Statement.**
- z_{l+1} = relu(z_l) W_{l+1} is a deterministic function of z_l, so (z_1, …, z_L) is a (degenerate) Markov chain for every network.
- The face process σ_l = sgn z_l is a function (lumping) of it. A lumping is Markov iff the lumpability condition holds (Kemeny–Snell; Rosenblatt; from memory). Otherwise its memory is created by marginalising the coordinates within faces: fill-in (unlock 43).
- Every Gibbs law on histories is Markov at every layer (unlock 7).
- Markov at a layer is equivalent to a commuting square of the past, future and present conditional expectations inside the history algebra D (U8; the commutative case is proved).

**Status.** DERIVED (elementary; U8 as in the hdx digest).

**Computational meaning.**
- Memory is a property of the resolution, not of depth. The full-resolution state, the law of z_l, is a perfect separator, so an estimator that carried it exactly would need no history.
- Every coarser state (faces; a Gaussian field plus sites; windows) has a memory cost, and that cost is exactly computable (unlock 45).
- Measured:
  - dictionary v1's face process has 25–35 % memory (mlp-bridge);
  - 40 % of the pairwise joint structure feeding the next layer is older than one layer (BRIEF §3).

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

  with equality iff the Gaussian is Markov on the graph (iterated Koteljanskii; equivalently, maximum entropy given the clique marginals). The gap equals the KL divergence to the Markov projection and the sum of the cut CMIs.
- (d) CMI is invariant under injective transformations of each coordinate. So for a Gaussian copula (z_i = T_i(g_i), T_i increasing, g Gaussian), the CMI of z is the Gaussian CMI of the latent correlation.

**Status.** THEOREM ((a), (b), (d) classical; (c) maximum entropy) + checked (C1; C2, where the gap 1.2277 equals Σ CMI = 1.2277 and a tridiagonal precision gives 0 to 4e-16).

**Computational meaning.**
- *Exactly and cheaply:* for every quasi-free or copula state, every cut CMI is a ratio of determinants. The Gaussian parts of the faces, signings and markov designs' states are of this kind.
- The memory budget of unlock 45 then costs O(s³) per separator of size s.
- The Markov defect can be monitored inside the estimator at no sampling cost.

**Guard.** (d) applies to pre-activations, not to activations. relu is not injective, so the CMIs of a_l and z_l differ.

### 47. Separator along depth, tree across width: pseudorandomness is the switch ⊕ [K, B, M]

**Statement.**
- (a) Separator gluing (junction trees, Hammersley–Clifford, Yang's CMI bounds) is exact, or controlled by CMI. For a general law it costs exp(separator size) (treewidth). For a quasi-free (Gaussian) law, conditioning on a separator is a Schur complement, O(s³), and its CMI is a determinant ratio (unlock 46).
- (b) Tree gluing (Bethe, cavity, TAP) is approximate, controlled by the decay of correlations along computation trees.
  - It is strong on pseudorandom geometry, and weak at low temperature and on amenable geometry.
  - Expanders are the worst case for separators, since Yang's boundary gain vanishes on them. They are the best case for trees (unlocks 12, 16).
- (c) The network supplies both geometries.
  - Depth is a path of exact separators (unlock 44), but each separator is a whole layer: size n, and not quasi-free because of the faces.
  - Width is a dense pseudorandom bipartite geometry through fresh weights, where only tree-shaped contractions survive (unlock 13).

**Status.** DERIVED (a synthesis of the cited theorems: the user's Idea B, checked against the programme's results).

**Computational meaning.** The design rule this produces:
- Glue along depth through quasi-free separators: Gaussian or Gaussian-copula fields, carried exactly by linear algebra at O(n³) per layer.
- Attach the non-quasi-free content (faces, hubs, unary non-Gaussian potentials) as local decorations, glued across width by tree rules with the Onsager correction.
- Charge the residual to two measurable defects: the cut CMIs along depth (unlock 45) and the cycle content across width (unlock 13).

**Guard.** Old content (BRIEF §3) is a cross-layer joint structure with no low-rank form. Neither rule is local there; the Heisenberg evaluation (unlocks 35–36) is the alternative.

### 48. Two defects on one stage: angles and cocycle leakage [K, H, B]

**Statement (A-unification).** Fix a faithful state and its forget-A algebras N_A.
- (a) The angle c(A, B) = ‖E_A E_B − E_{A∪B}‖ behaves as follows:
  - it is zero on separated pairs iff the field is Markov;
  - it decays across a buffer under strong spatial mixing;
  - it equals λ(G)/d on the two ends of an edge (expander mixing);
  - on adjacent binary spins it is the geometric mean of the two Dobrushin influences.

  The sharp two-projection inequality is (1 − c) Var_{A∪B} ≤ Var_A + Var_B. Alternating projections converge at rate ‖(E_A E_B)^k − E_{A∩B}‖ = c^{2k−1} (Aronszajn; Kayalar–Weinert; from memory).
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
| 13 | power counting for a fresh layer | | ● | ● | | | ● | tree contractions across width suffice; joint structure of the state is necessary |
| 14 | Bethe = cover limit; network lift limit | ● | ● | ● | | | | an exact, computable tree reference |
| 22 | sign-averaging = tree projection | | ● | ● | ● | | | the scored information is sign-odd |
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
| 47 | separator along depth, tree across width | | ● | ● | | | ● | the gluing rule for each direction of the network |

---

## 9. The five most important statements

**I. The estimand has an exact skeleton: a linear mean recursion plus two wall/face scalars per unit (unlocks 3, 37, with 1 and 4).**
- E z_{l+1} = E[a_l] W_{l+1} exactly.
- For every unit, E relu(z) = μ P(z > 0) + E[δ(z) Γ_z]. Here Γ_z = ⟨∇z, −∇L^{−1}(z − μ)⟩ is a two-replica overlap of input-space gradients, and its mean is exactly Var z.
- Gaussian closure is the approximation "Gaussian face mass, Gaussian wall density, Γ ≡ Var z". Its error splits exactly into drift × (face-mass defect) + (wall defect).
- Globally, E f = E Δf: the mean is the Gaussian mass of the tropical hypersurface, which by summation by parts is the statement that barycentres of faces are facet integrals.
- The uncentred form of this identity puts the mean spike into the wall weights and is non-perturbative (tropical R-E2). The centred form does not.

**II. Fresh randomness organises the width: annealed = twirl, sign average = tree, score = sign-odd (unlocks 5, 13, 22, 12).** W_{l+1} is independent of the law of a_l. Consequently:
- the annealed one-step map is a perfect expander;
- averaging over the fresh layer's signs is exactly the paired, tree sector (Godsil–Gutman);
- the quenched per-unit information is sign-odd. It has two parts: an O(1) coherent alignment with the mean direction, which the exact recursion carries, and an O(n^{−1/2}) incoherent pairing of fresh weights with the state's joint structure.

Loops formed by fresh weights alone are suppressed by n^{−1/2}. So across width the gluing rule is tree-shaped (Bethe/TAP with the Onsager term, certified only by ℓ² quantities), and everything beyond the mean recursion lives in the state's joint structure.

**III. Depth is exactly Markov at full resolution; any coarser state pays a computable sum of CMIs (unlocks 44–47).**
- Memory is a property of the resolution. Its price is D(P ‖ P_JT) = Σ_cuts CMI. This equals the KL distance to the Gibbs-on-histories (KMS) family and the sum of the Petz sufficiency defects.
- For quasi-free (Gaussian or Gaussian-copula) separators, conditioning is a Schur complement, and each CMI is a ratio of minors: I = ½ log(1 + y), with y the ratio of the two monomials of the Desnanot–Jacobi exchange relation. Markov is the degeneration y = 0.
- Design rule: glue along depth through quasi-free separators, attach the face (non-quasi-free) content as local decorations glued across width by trees, and set the window length by the measured CMI decay.

**IV. States forward, questions backward, with a bilinear error (unlocks 35, 36, 39).**
- The pairing of the law at layer l with the pulled-back readout is invariant. Duhamel telescoping writes the error of any forward reference exactly as Σ_l (local defect of arrow l)[exact pulled-back question].
- With a model of the question, the remainder is bilinear: defect × model error.
- There are exactly n questions per layer, so backward evaluation works at full rank. Forward compression of signed content is impossible (k_ε/n is constant).
- State errors reach the output through the Lipschitz constants of the questions: contracted in the bulk (gain 2E[Φ²] < 1 per layer) and amplified along the mean direction.

**V. Positivity is the only certified source of forgetting and compression; the scored content is signed (unlocks 40, 18, 29, 12).**
- The positive sector has explicit, hereditary certificates. This sector is the face masses, Gibbs laws on histories, and the rectified mean direction. Its certificates are TV ≤ tanh(Δ/4), uniform window gaps, forgetting ∏ tanh(Δ_l/4), and one Perron/BBP mode.
- The signed sector has none:
  - the free-probability bulk has no gap (PR ≈ n/(2·age), k_ε/n constant);
  - the path sum is a weak-disorder polymer with no dominant path;
  - the zero-temperature Newton-polytope form cancels by 10²³ at n = 1024.
- No design may rely on decay, low rank, or a tropical skeleton of the signed sector. It must be carried at full rank or evaluated backward. On the dense pseudorandom geometry only spectral (ℓ²) certificates are usable; Dobrushin-type (ℓ¹) ones fail by √n.

---

## 10. Guards, and negative results charged to the dictionary

Each guard is stated with the dictionary item it is charged to (BRIEF rule 3). None of them refutes a general statement.

- **Uncentred wall split** (unlock 4; tropical R-E2). Charged to the dictionary "walls weighted by input-space slope²". The centred identity (37) is the repair at the level of the identity. Whether its wall correction is O(n^{−1/2}) is the CONJECTURE stated in unlock 37, tested as in §11.
- **The zero-temperature skeleton** (unlocks 29, 30, 32). The 10²³ cancellation, the weak-disorder polymer, and hot gates (a quarter of them at depth 16) are charged to "dominant cone or dominant path as the state". The tropical content of the theory survives as the signed wall identities and as a cost saving on the frozen sub-frame.
- **Forward compression of signed content** (unlocks 18, 19): no gap, k_ε/n constant, no low-rank form below ≈ 0.3 n. Charged to "carry old content in a fixed small basis". The theory's answer is backward evaluation (unlocks 35–36).
- **Dictionary v1** (unlocks 8, 44): 25–35 % face-process memory, and a residual of 0.42–0.5 against coboundaries. Charged to the resolution of v1. Unlock 45 prices it exactly.
- **Positivity hypotheses** (unlocks 7, 8, 40, 41). The certificates apply to Gibbs laws on histories and to positive kernels. Signed transport is outside them, which is a statement about scope, not a failure.
- **ℓ¹ certificates on dense layers** (unlock 11): Mooij–Kappen, Dobrushin, and BP certificates based on |J| fail by √n. Use ℓ² certificates (unlock 12).
- **Spectral independence at η ≈ 7** (unlock 17) gives a bounded variance factor, not precision. Mixing certificates built on it are numerically vacuous.
- **Determinantal, real-stable and totally positive tools** (unlocks 24, 26, 34). K_{n,n} is not Pfaffian, and the face law is neither strongly Rayleigh nor MTP₂. These tools apply to sign-averaged (matching) objects only.
- **Sign averages** (unlock 22) are wrong per unit at O(1). They are references, never estimators.
- **Copula invariance** (unlock 46(d)) does not pass through relu: it holds for pre-activations only.
- **Annealed asymptotics** (unlock 31): 9π²/(2l²) overestimates 1 − c_16 by a factor 2.5. Use the iterate.
- **The user's Idea C** (a periodic quantum circuit from quantum cluster mutations). No computational meaning for quenched means was found; it is not used. SPECULATION.

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
  - The random-lift limit (unlock 14(b)) is an exact, computable tree reference with exact non-Gaussian single-site laws. Use it to measure loop content as quenched minus lift.
  - Dense certificates must be ℓ², with the Onsager term (unlock 12).
- **signings.**
  - The sign-sector decomposition (unlocks 5(d), 22) is the theorem behind §1 Step A.
  - The Pfaffian guard (unlock 24) closes exact determinantal evaluation on dense layers.
  - Signed-lift control variates (unlock 25) need an identity that makes the network's mean a sign average before they can be used.
- **heisenberg.**
  - The Lipschitz gains of unlock 39 give an error budget: amplified along the mean direction, contracted in the bulk.
  - The centred identity (unlock 37) is the readout question in Gaussian-space form.
  - Theorem 2's bilinearity means ĝ needs ≈ 2–3 % relative accuracy, not 2e-4 absolute.
- **markov.**
  - For the quasi-free (Gaussian, copula) part of the state, the CMI budget is closed-form (unlock 46). Monitor it inside the estimator.
  - Measure cut CMIs against age to set d (unlock 45).
  - Copula invariance holds for pre-activations only.

---

## 12. Numerical checks

Run in the session scratchpad with numpy (scipy is not installed in the whest environment; normal functions come from math.erf), OPENBLAS_NUM_THREADS=1, seeds fixed. Each check supports the claim it is attached to; none is a claim of its own.

| check | claim | result |
|---|---|---|
| C1 | Dodgson exchange relation; I(a:c\|B) = ½ log(1 + y) (unlocks 33, 46) | three random 6 × 6 SPD matrices: relation residual ≤ 2e-15; CMI equal to 10 digits (0.2632371887, 1.2337629262, 0.1329056055) |
| C2 | chordal (path) junction tree: log-det gap = KL to the Markov projection = Σ CMI; zero iff Markov (unlock 46(c)) | generic 6 × 6: 1.227709586 = 1.227709586; tridiagonal precision: −4.4e-16 and 0 |
| C3 | Sheppard and trivariate orthant formulas (unlock 9) | Monte Carlo (8e6 samples) within 1σ in all six cases, e.g. 0.24689 ± 0.00015 vs 0.24702; 0.08717 ± 0.00010 vs 0.08726 |
| C4 | Godsil–Gutman: E_ε det(ε ∘ √A)² = per(A) (unlock 25) | exact enumeration over all sign patterns: n = 3, 0.1946602504 = 0.1946602504; n = 4, 1.1579296180 = 1.1579296180 |
| C5 | arc-cosine iterate from c_0 = 0 (unlock 31) | c_l = 0.318, 0.494, 0.681, 0.834, 0.897, 0.929 (l = 1, 2, 4, 8, 12, 16); 1 − c_16 = 0.0705 against 9π²/512 = 0.1735; l²(1 − c_l) = 43.72 at l = 2000 against 44.41 |
| C6 | Var(wᵀBw) = 8‖B‖_F²/n² (unlock 5(c)) | n = 200, 2e4 samples: 3.899 against 3.947 (Monte Carlo error ≈ 1 %) |
| C7 | E relu(Z) − relu(m) = s[φ(t) − tΦ̄(t)] (unlock 30) | 2e7 samples with a control variate: 8.643e-3 vs 8.643e-3; 7.123e-3 vs 7.138e-3; 1.919e-4 vs 1.911e-4; 3.83e-6 vs 3.57e-6 |
| C8 | centred wall identity (unlock 37), Z = Σ_i v_i relu(⟨x, w_i⟩), d = 3, h = 6, hot unit (μ/s = −0.07) | 4e6 samples: E Γ_Z = 0.40966 vs Var Z = 0.40967 (exact arc-cosine value); smooth form E[(Z−μ) tanh(Z/ε)] = E[tanh'(Z/ε) Γ_Z/ε] at ε = 0.3, 0.1, 0.03: 0.4581/0.4575, 0.4909/0.4898, 0.4954/0.4942; wall form 0.2478 vs 0.2473 (ε = 0.02 kernel); E relu Z = 0.22909 vs μP + wall = 0.22852; uncentred E‖∇Z‖² = 0.550 vs μ² + Var = 0.412 |

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
| Kasteleyn (planar); 4^g Pfaffians at genus g | THEOREM (planar case via Lieb–Loss, verified; genus count from memory) | 24 |
| Pfaffian bipartite graphs (RST, McCuaig) with the Heawood exception | THEOREM (RST verified: the brace characterisation and Little's K_{3,3} criterion) | 24 |
| Kasteleyn's face rule equals Lieb's flux rule | THEOREM (Lieb 1994 verified: π per square plaquette, 0 per hexagon) | 24 |
| Godsil–Gutman; roots ≤ 2√(d−1); MSS Ramanujan 2-lifts; Kadison–Singer | THEOREM (Godsil–Gutman, KS/paving verified; Ramanujan 2-lifts from memory) | 22, 26 |
| "the diagonal is the classical shadow" | ANALOGY: D_0 ⊂ C*(Λ) plays the role of the diagonal masa (unlock 6); paving as a block mean field is SPECULATION (unlock 26) | 6, 26 |
| tropically per = det = assignment | THEOREM | 32 |
| Chin's proof: tropical c-vectors and T-systems, lifted by positivity and periodicity | THEOREM (Chin 2602.15140 Thm 1.1 verified; the lift is the IIKKN periodicity theorem) | 32 |
| the general permanent is #P-hard; Kasteleyn and JSV cross to finite temperature | THEOREM (from memory); JSV's crossing rests on positivity | 32 |
| cluster variables and T-systems are matching sums; Kuo condensation; urban renewal; Goncharov–Kenyon; Galashin–Pylyavskyy | THEOREM (from memory); context only, no computational use found | 33 |
| MTP₂ distributions are faithful | THEOREM with the graphoid hypothesis (Fallat et al. Thm 6.1, verified) | 34 |
| Idea A: the Dodgson identity, and CMI as a function of its two terms | DERIVED and checked (C1); sharpened to I = ½ log(1 + y) | 33, 46 |
| Idea A: a cluster structure on Gaussian states with Markov networks as boundary strata | SPECULATION | 33 |
| Idea B: separator vs tree local-to-global, with pseudorandomness as the switch | DERIVED for the network as unlock 47, refined by quasi-free separators; the quantum half is a CONJECTURE (B-programme 8.4) | 47, 49 |
| Idea C: a provably periodic quantum circuit | SPECULATION; not used | — |

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
- J. M. Mooij, H. J. Kappen, *Sufficient conditions for convergence of the sum-product algorithm*, arXiv:cs/0504030.
- N. Robertson, P. D. Seymour, R. Thomas, *Permanents, Pfaffian orientations, and even directed circuits*, Ann. of Math. 150 (1999) 929–975 (with Little's theorem as stated there).

**Through the programme's digests (cited by digest section).** Chen–Rouzé arXiv:2504.02208; Yang arXiv:2609.38007; Hanin–Nica arXiv:1812.05994; ALO, CLV, Chen–Eldan, Alev–Lau, Oppenheim, Dinur–Kaufman (hdx digest); EKZ, Hoory–Linial–Wigderson, Willett–Yu (expanders digest); Cipriani, Goldstein–Lindsay, Vernooij–Wirth, Carlen–Maas (nc-Dirichlet digest); Hammersley–Clifford, Brown–Poulin, Leifer–Poulin (hammersley-clifford digest); Fawzi–Renner and JRSWW (Chen–Rouzé digest).

**From memory (standard; not re-opened).** Isserlis 1918; Sheppard 1899; Plackett 1954; Price 1958; Kahane 1986; Cover 1965; Wendel 1962; Kemeny–Snell 1960; Sudakov 1978; Diaconis–Freedman 1984; Mingo–Speicher 2006; Weitz 2006; Mossel–Sly 2013; Kesten–Stigum 1966; Godsil–Gutman 1981; Godsil 1981; Marcus–Spielman–Srivastava, Interlacing Families I; Barvinok 2016; Patel–Regts 2017; Borcea–Brändén–Liggett 2009; Valiant 1979; Jerrum–Sinclair–Vigoda 2004; Galluccio–Loebl 1999; Tesler 2000; Füredi–Komlós 1981; Bai–Yin 1988; Bolthausen 2014; Aronszajn 1950; Kayalar–Weinert 1988; Hubbard 1959 and Stratonovich 1957; Cho–Saul 2009; Nourdin–Peccati (Stein kernel and Malliavin calculus, 2009–2012); Houdré–Pérez-Abreu 1995; Maslov dequantisation.
