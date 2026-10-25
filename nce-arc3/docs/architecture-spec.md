# Final architecture: the Noncommutative Correspondence Engine (`nce`) for ARC-AGI-3

Synthesis of the four panel designs under the three judgments. Spine: **fock-ladder** (first on two of three rankings: Fock vector of a frame, charge sectors, graded residual, second-quantized goals, goal-splitting, bound-state promotion, gauge reduction). Grafts: **moment-map** (family dimensions as summand counts, animation frames as serial compositions, RESET learned, goal positives on the predicted post-state, the probing loop, the Hankel-rank classification of strata), **two-projection** (stratum refinement under residuals, decision-relevant interleaving, sparse-diff hashing of predictions), **kms-library** (the exact E6-corrected promotion inequality and its corollaries, units in the vertex algebra, Letter/ClauseSet/Frontier/Store interfaces, residual-guided frontier pruning, Zipf cost reading of item 3). Every smuggled prior the judges found is removed in Section 0; every spec bug they found (d = 1 cycles, missing fill rule, missing `close` edge, goal positives on the wrong state, the hard depth cap, unbounded per-action loops, the empty pre-completion goal set) is fixed where it occurs.

Theorem labels T1–T13 are the brief's; E1–E12 and items 1–11 are the amendments'; C1–C11 are the brief's constraints; N1–N17 the negative checklist. Monograph references are to the `.tex` (Definition 1.1, Prop. 1.3, Thm. 1.4, Thm. 1.5, Thm. 2.1, Thm. 3.1, Thm. 3.2, Thm. 4.1, Prop. 9.3, Components 1–5).

---

## 0. Governing rule and status tags

**Governing rule.** Nothing in this document comes from a game. No frame was inspected, no source read, no mechanic named. Every component below is traced to a theorem or amendment item by label; where a choice was needed that the theory does not make, the choice is tagged and justified from the interface facts of brief Part A only (the palette is `0..15`, the lattice is `64×64`, action ids are `0..7`, `levels_completed` and the state flag are observed, `win_levels` is given).

**Status tags**, used on every component: **[T]** derived from a theorem the brief marks sound; **[T\*]** derived by transfer with the amendment's citation; **[D]** derived inside a model (a convention forced by making a [T] statement operational, with the model stated); **[C]** conjecture, used only to motivate a fallthrough, never a primary rule.

**Audit of the judges' smuggled priors and spec bugs.** Each is either removed or re-derived; nothing is kept on the basis of "reasonable".

| # | Finding (design) | Resolution in this spec |
|---|---|---|
| 1 | `Trans_ℓ`, d = 4ℓ (fock-ladder): word-metric prior dressed as a Fock grade | Removed. The translation group enters as two families: `Step` = the generator set S of the gauge group's Dirichlet form, d = 4, and `Move` = the group as a set, d = 16128. `Step^k` is the k-th tensor power of ⊕_{g∈S}A and has 4^k summands (Thm 1.5(iii)), so it outweighs `Move` for k ≤ 6. Small displacements are preferred only as a consequence of KMS order [T]. |
| 2 | `Const(Ψ_0)` hard-codes RESET (fock-ladder, kms-library cost model, two-projection cost model) | RESET is an ordinary action id whose word is learned from its first demonstration (Section 4). A `Const` letter exists (a constant correspondence is a legitimate many-to-one letter whose fan-in is its empirical domain size) but nothing is assumed about which action, if any, it explains. The planner may use a verified RESET word (fock-ladder's good idea, kept). |
| 3 | `enclosure(ξ,η)` imports the lattice boundary (fock-ladder) | Removed. |
| 4 | "adjacency at distance ≤ L" (fock-ladder) | Removed. The only predicates on the relative gauge element are `h = e` (gauge coincidence) and `h ∈ S` (one generator apart), both group-theoretic with no constant. |
| 5 | Newborn-seeking over grade-2 contexts as a curiosity bonus (fock-ladder) | Demoted to the last fallthrough and restated in the form C5 licenses verbatim: when hypotheses agree and no goal is known, choose the action whose predicted outcome changes the *goal observables* in a way not yet seen (C5, third clause). Goal-splitting (Kato on the goal coherent state) runs first. Tagged [D] with the concern recorded. |
| 6 | Cell weights `1/(r_F·|fibre|)` (two-projection): background down-weighting | Removed. All inner products use the flat one-hot metric; Thm 3.1(iii) gives slopes, not mass per stratum. |
| 7 | Goal alignment `1 − |Δt|₁/128`, `1 − m/m₀` (two-projection): distance-to-target shaping | Removed. Pre-completion ranking uses Component 4 alignment only in the form "does any observed transition move the observables this goal word reads" (Section 6.3). |
| 8 | Colour bijections in the gauge group with Haar prior (two-projection) | Removed. Recolouring is a letter `Recolour` with d = 15 (the other palette colours); the gauge group is D4 ⋉ Z² and shrinks only by Section 8.4. |
| 9 | Parry top layer `±e1,±e2` for untried actions (two-projection, kms-library) | Removed. Untried actions are tried once each (an action with no demonstration has every outcome newborn, T2/C8); ignorance-state actions are scored only from partial words or in newborn contexts (Section 5), never from a hand-picked layer. |
| 10 | Nominal d = 2 for α, κ, μ, ω (kms-library) | Removed. Every dimension is a summand count fixed by the interface (palette 16, lattice 4096, D4 8, `win_levels`) or a snapshot of an observed alphabet (Section 3.3). |
| 11 | Nearest-neighbour matching by translation norm (kms-library) | Removed. Matching maximizes KMS weight of the explaining word under the `Step`/`Move`/`Rot`/`Recolour` families (Section 4.3), the same quantity synthesis maximizes. |
| 12 | Painting in HN order as a silent background rule (all four) | Stated explicitly as a rendering convention [D] (Section 4.4): the HN flag is a filtration and the canvas realizes it. The colour of a *vacated* cell is never assumed: it is a fitted parameter of the clause (`Fill`, d = 16), fixed by the first demonstration and verified exactly. |
| 13 | Vacuum stratum untracked, canvas fill, `6@Ω` (moment-map) | Removed. Every stratum's components are grade-1 vectors; no stratum is ambient. |
| 14 | `s = e_c ⊗ (1,x1,x2)` and the "exact" HN closed form (moment-map) | Removed. `s = e_c` (one-hot) is the minimal gauge-covariant nowhere-zero map; the flag is exactly colours by mass (Section 2.1). |
| 15 | "Vacuum-coloured components fully enclosed: keep" (moment-map) | Removed. |
| 16 | ACTION6 reduced to one action per type (moment-map) | Replaced by the `SelPointed` selector (the clicked component is the argument) plus one representative cell per component in the candidate set (C2: cells of one component are unitarily equivalent in M_{n_k}). |
| 17 | "Top goal carries > ½ of the weight" threshold, H = 8 (fock-ladder) | Threshold removed: goal-splitting always maximizes W(s)(W_tot − W(s)), which targets a dominant goal automatically (Section 6.4). Horizons and node counts are compute budgets with anytime fallbacks, listed in Section 9, not behaviour constants. |
| 18 | Coherent set = equal-weight ties only (fock-ladder); first exact word only (moment-map) | Every exact clause set found within the synthesis budget is kept, per action id, ordered by weight (T7); the cap is a compute cap. |
| 19 | Hard depth cap `L_max(T)` counting letters against `r(A_0)` (kms-library, fock-ladder) | Removed as a refusal rule. T10's ceiling is used as the *promotion trigger*: an exact word longer than the ceiling acts (T6 and T7 guard it) but is flagged `beyond_ceiling`, and its subwords enter the maternity ward at once (Section 4.6, 8.2). |
| 20 | Death priced at zero (kms-library); "never plan through GAME_OVER" (moment-map) | Death is priced by replay length under the verified model (Section 5, 7), never forbidden, never free. |
| 21 | Goal positives on the pre-state (two-projection) or on the next level's first frame (kms-library) | Positives are evaluated on the model's predicted post-state of the completing action (Section 6.2). The two-projection contradiction (s\* both kills and keeps) cannot arise. |
| 22 | `G_live` empty before the first completion (moment-map code) | The goal coherent set is initialized from the Parry prior at level start, filtered by "false on every observed frame" (Section 6.2). |
| 23 | d = 1 cycles (fock-ladder `Select`/`Annihilate`/`Const`) so β_c does not exist | The d = 1 subgraph is acyclic by construction and asserted in code (Section 3.2): d = 1 edges are `SelPointed`, `SelAll` (F→C) and `Over` (F→⊤, absorbing); every C→F edge has d ≥ 16; every F→F loop has d ≥ 2. |
| 24 | Missing `close` edge, undefined fill, unspecified Parry cursor (fock-ladder) | A clause ends at C (the components persist) or at F through `Annihilate`; no `close` edge is needed. Fill is item 12. The cursor is specified in Section 4.6. |
| 25 | Unbounded `resynthesize`, 20k-state BFS, 30-iteration Tyler loop (moment-map); `align` needing a QR per candidate (two-projection) | Every loop has a wall-clock budget with carry-over (Section 9); Tyler iteration is dropped (the gauge is already compact, so canonicalization is orbit-representative selection); `align` is replaced by residual-guided pruning (Section 4.6). |

**Negative checklist.** N1 (pixel vector): the representation is the Fock vector of components; pixels are evidence only in the exact test. N2 (frame hash as state): hashes exist only in the ignorance table and the planner's visited set. N3 (pixel-diff novelty): the learning signal is the model-conditional newborn residual, never "changed". N4 (random exploration): the Kato rule; no ε. N5 (LLM): none. N6 (game DSL): every letter is a gauge element, a Fock ladder operator, a colour-fibre map, or a gate on a base-algebra invariant. N7 (magic constants): the only constants are interface facts and compute budgets; β_c is computed. N8 (convolutional bias): none. N9 (RL): none. N10 (discard animation): intermediate frames are serial compositions (Section 4.8). N11 (per-level amnesia): Section 8.5. N12 (blob click at random pixel): one deterministic representative cell per component. N13 (agent-written code): none. N14 (free-text memory): the library is the memory. N15 (grid = full state): hidden state is detected and recorded (Section 4.7). N16 (search over the true environment): search runs only inside the verified model. N17 (baseline-tied knobs): no h, no action cap anywhere.

---

## 1. The system in one page

A frame is a probability measure on the 64×64 lattice with a one-hot colour label. Its Harder–Narasimhan flag (T8, Thm 3.1(iii)) orders colour fibres by mass; inside a fibre, the kernel of the Dirichlet form generated by the translation gauge unitaries (§1.3, Component 1) yields 4-connected components. Components are the grade-1 vectors of a Fock space over the typed base algebra; pairs are grade 2; masses, counts and the state flag are grade 0 (T13, item 9). A transition is an operator on this Fock space. It is decomposed first into charge sectors by the gauge torus (Prop. 1.3: per-colour component-count differences), then, in charge 0, into moves T_{hξ}T_ξ\* matched by maximum Parry–KMS weight (T5), with the gauge element read off by transform-then-match (first moments give the translation, second moments give the D4 element, fibres give the colour; C11, E7, item 2), and finally verified cell for cell (T6).

Learning is closed-loop residual heating (T1) organized by grade (T2): every surviving hypothesis for the action taken predicts the next Fock vector; the newborn residual is computed; hypotheses with nonzero residual die; when all die, new clause sets are enumerated from the typed library in Parry–KMS order at the critical β_c of the current library (T5, T11), fitted on transformed residuals, and accepted only on exact reproduction of every stored transition of that action. All exact clause sets found within the compute budget are kept as the action's coherent set (T7). A lookup table of stored transitions exists as the ignorance state and never acts (T12).

Exploration is the Kato criterion (T3): the next experiment is the (action, datum) on which the coherent set's predicted frames have the largest KMS-expected tangent tax Σ p_h p_h' sin²θ, divided by a cost that prices predicted death by replay length. When a goal plan exists, candidates are restricted to the plans' first actions under each surviving hypothesis (T7: commit only when a decision requires it).

Goals are second-quantized predicates dΓ(p), Γ(p) of words from the same library (C6), ranked by the same prior, killed by every observed non-winning post-state, confirmed on the predicted post-state of each completing action, carried across levels. Before any completion the agent runs Kato on the goal coherent state: it plans to the reachable state that splits the goal weight most evenly, and each arrival either wins or kills a slope class of goal words.

Memory is the library. Composite words used often enough are promoted to letters iff the code length of everything held decreases after β_c is recomputed exactly (T5 with E6, T11); conserved pairs are promoted to grade-1 types (C2); the gauge group shrinks to the commutant of the synthesized words (T8); promotion is forced when synthesis hits the in-context depth ceiling (T10). Nothing is carried across games except the library under the reversed test (C7).

```
obs ─► perceive (T8, §1.3, T13) ─► Ψ_t
Ψ_{t-1}, a, data, frames ─► charge + match (Prop 1.3, T5, C11) ─► GradedTransition(s)   [animation: ∅-steps, T4]
GradedTransition ─► residual test on H[a] (T1, T2, T6) ─► deaths; if H[a] = ∅ ─► synthesize (T5, T11, T12, C11)
Ψ_t, levels, flag ─► goal/death update (C6, T7)
H, G, Ψ_t ─► decide: interleaved Kato (T3, E10, T7) | plan (T4) | goal-split (T3 on goals) | goal-observable newborn (C5)
level end ─► promote / bound states / gauge reduction / recompute β_c (T5+E6, T10, T11, T8, C2)
```

Persistent state (the "weights", C1): `Lib` (letters, dimensions, β_c, Perron vectors, gauge subgroup H, promoted letters with abstracted type slots, bound types, split rules), `H[a]` (coherent set per action id), `Store[a]` and `Store∅[a]` (demonstrations), `Ignorance[a]`, `Goals`, `Deaths`, `Negatives`, `Positives` (pending and confirmed), `Frontier[a]` (resumable Parry cursor).

---

## 2. Base algebra construction from the frame-as-measure (T8)

### 2.1 The frame measure and its HN flag [T]

Read a frame F : X → {0..15}, X = Z_64², as the probability measure ν = 4096⁻¹ Σ_x δ_x with the measurable map s(x) = e_{F(x)} ∈ C¹⁶. This s is the minimal one that is nowhere zero, carries the colour fibre, and is covariant under D4 ⋉ Z² (the gauge unitaries permute cells and fix colours); no affine coordinates are adjoined (audit item 14).

Thm 3.1(ii): a Tyler balanced metric exists iff ν(s⁻¹V) < dim V/16 for every proper subspace V. For V = span{e_c : c ∈ S} this reads Σ_{c∈S} ν_c < |S|/16 for every colour set S. Taking S = {argmax ν_c}: a 16-colour one-hot measure has some ν_c ≥ 1/16, so stability always fails. Thm 3.1(iii) then gives the unique HN flag. Because s is supported on the coordinate axes, the destabilizing subspaces are coordinate subspaces and the slope of a colour block is its mass; the flag is the sequence of coordinate subspaces spanned by colours in strictly decreasing order of mass, colours of equal mass forming one polystable block. **Strata = colour fibres ordered by mass, ties grouped.** No stratum is named; the top stratum merely has the largest slope. The HN *rank* of a colour is its position among the distinct masses.

### 2.2 Components as the kernel of the gauge Dirichlet form [T]

§1.3 attaches to (ℓ^∞(X), counting measure) and a symmetric Markov semigroup the Cipriani–Sauvageot correspondence with derivation ∂ and carré du champ Γ. The only semigroup the theory hands us on the lattice is the one generated by the gauge group's translation generators (Component 1): Δ = Σ_{g∈S}(I − U_g), S = {±e₁, ±e₂} the minimal D4-invariant generating set of Z². Then ∂f = (U_g f − f)_{g∈S} and Γ(f,h) is supported on lattice-adjacent cell pairs. Restricted to a colour fibre, ker ∂ is the space of locally constant functions, whose minimal projections are the **4-connected components** of the fibre. 4-connectivity is the generator set, not a choice; diagonal adjacency is a length-2 word in S and appears only as a composite. Each component indicator 1_ξ is a frame vector of the fibre's base algebra (item 9).

### 2.3 Canonical form modulo the gauge group [T]/[D]

The gauge group is G = H ⋉ Z² with H ≤ D4, initially H = D4 (Component 1; reduced only by Section 8.4). Thm 3.1(ii) says whitening reduces an invertible gauge to a compact one; our gauge is already compact, so there is nothing to whiten: canonicalization is **orbit-representative selection** [D: the model is "the gauge is the discrete group Component 1 names"]. The canonical form of a component is the lexicographically minimal image of its bounding-box-cropped mask under the |H| ≤ 8 elements of H; the gauge element of the component is (g, t) with g ∈ H the minimizer and t ∈ Z² the bounding-box corner of the observed mask. The second-moment tensor Σ_ξ = n⁻¹ Σ_{x∈ξ}(x − x̄)(x − x̄)ᵀ is D4-covariant and is used as a bucket key before the exact lexmin (it is the moment map at G = I; nothing more is claimed for it). The **type** of a component is τ(ξ) = (colour, shape key) where the shape key is the canonical mask; the shape alphabet Ξ is the set of keys seen so far in the game.

### 2.4 Grade-2 observables from the same structures [T]

For a pair (ξ, η): `contact(ξ,η) := Γ(1_ξ, 1_η) ≠ 0` (shared lattice edges, §1.3); `overlap(ξ,η) := ⟨1_ξ, 1_η⟩ ≠ 0` as cell sets after a predicted move (the frame Gram of Thm 3.2); the relative gauge element `h(ξ,η) := (t_η − t_ξ, g_ξ⁻¹g_η) ∈ Z² ⋊ H` with the two predicates `h = e` (coincidence) and `h ∈ S` (one generator apart); type equality, colour equality, shape equality modulo colour. Nothing else is a relation (audit items 3, 4).

### 2.5 Strata earn their status (C2): split rules and bound types [D]

Two mechanisms, both triggered only by stored transitions, change what perception emits. **Downward** (two-projection): if a verified-or-candidate word's residual changes a strict subset of a component's cells on some transition, the component is re-filtered by ∂ restricted to changed and unchanged cells, and the split is stored as a rule `shape_key → partition` applied at perception for the rest of the game (the HN filtration of the transition residual). **Upward** (fock-ladder, Section 8.3): a pair whose relative gauge element and contact are conserved under every stored transition in which either member moves is promoted to one grade-1 type; perception thereafter emits the union as one component with its own canonical form. Multi-colour things earn existence by stability under the observed correspondences, never by a prior.

### 2.6 The Fock vector of a frame and the algorithm

```
Component : colour:int, cells:int32[n] (flat indices, sorted), n:int, bbox:(r0,c0,r1,c1),
            g:int (0..7), t:(int,int), key:bytes, type_id:int, rank:int, moments:(f,f,f)
FockVector: frame:int8[64,64], mass:f32[16], count_c:int16[16], count_t:dict[type_id,int],
            flag:str, level:int, comps:list[Component], cell2comp:int32[4096],
            contact:csr_bool[n_comp,n_comp], hash:bytes

def perceive(frame, level, flag, lib):                     # T8, §1.3, T13, C2
    mass = bincount(frame.ravel(), minlength=16) / 4096
    rank = dense_rank_desc(mass)                           # HN flag: ties share a rank (Thm 3.1(iii))
    comps = []
    for c in colours_present(mass):
        lab, k = label4(frame == c)                        # ker ∂ on the fibre; 4-stencil = generator set S
        for i in range(1, k+1):
            cells = flatnonzero(lab == i)
            for piece in lib.split_rules.apply(cells):     # §2.5 downward refinement
                key, g, t = canonical(piece, lib.H)        # lexmin over H ≤ D4 of the cropped mask
                comps.append(Component(c, piece, len(piece), bbox(piece), g, t, key,
                                       lib.type_id(c, key), rank[c], moments(piece)))
    comps = lib.bound_types.bind(comps)                    # §2.5 upward: conserved pairs → one vector
    contact = contact_matrix(comps, frame.shape)           # Γ(1_i,1_j) ≠ 0 via the 4-stencil
    return FockVector(frame, mass, counts_by_colour(comps), counts_by_type(comps), flag, level,
                      comps, cell2comp(comps), contact, sha1(frame))

def canonical(cells, H):
    M = crop(mask_of(cells))
    imgs = [(apply_d4(M, h).tobytes(), h) for h in H]
    key, g = min(imgs)
    return key, g, bbox_corner_after(cells, g)
```

**Complexity.** `label4` on ≤ 16 binary planes of 4096 cells: O(16·4096) ≈ 1–2 ms with `scipy.ndimage.label` (or a numpy union–find); canonical forms O(Σ_ξ 8·n_ξ) ≤ O(8·4096) ≈ 1–3 ms with a cache keyed by the moment bucket; contacts by four shifted comparisons of `cell2comp`: O(4·4096). Total 3–8 ms per frame; frames with several hundred components stay under 15 ms because the per-component work is bounded by total cell count. There is no per-component iteration (audit item 25).

### 2.7 Hankel rank of a stratum in the action variable [T] (moment-map graft)

For each type τ, the tuple of verified words (w_{τ,a})_a over action ids has Hankel rank 0 (all identity: static), 1 (same word for all a: autonomous), or > 1 (action-dependent) by Thm 1.5(i) in Fliess's form (E11 restricts (i) to the Hadamard/Cauchy clauses that are used here). This classification is computed from H, not assumed, and is read by the planner (static types are frozen in the abstract state) and by goal ranking (a `Fix[τ]` goal word requires rank 0 under the verified words).

---

## 3. The graded library: atoms, dimensions, typing, Bratteli paths (T4, T13)

### 3.1 Vertices and paths [T]

V = {F (the frame; grade 0, vacuum sector), C (one selected component; grade 1), P (a selected pair; grade 2), ⊤ (terminal: the GAME_OVER flag; absorbing)}. Promotion adds *letters*, never vertices (type slots are abstracted, Section 8.1); bound states add *types* to Ξ, never vertices. Keeping |V| = 4 keeps β_c's Perron data meaningful and stable within a level (judge 3's objection to type-vertices). An explanation of an action is a set of Bratteli paths from F (T13, item 9): a clause ends at C (the selected components persist, transformed) or at F (annihilated). Two words compose by ⊗_B when the intermediate vertex matches (a word ending at C continues with a C→C letter).

### 3.2 Units live in the vertex algebra, not in the edge set [T]

A specific gauge unitary (one D4 element, one translation, the identity) has statistical dimension 1 (Thm 1.5(iii)). Thm 2.1(iv) requires d_e > 1 for β_c to exist: a cycle of d = 1 edges would pin r(A_β) ≥ 1 for all β. Component 1 already places the gauge unitaries inside A_0 = ⊕ M_{n_i} ⋊ D4. So units are **parameters fitted from data** (Section 4.3), never searched or weighed. What is searched and weighed is the *family* a parameter ranges over (Section 3.3), and an invertible verified word is adjoined to the gauge group of its type, never to E (kms-library corollary (ii), Section 8.1). The code asserts at every library change that the d = 1 subgraph is acyclic and that r(A_0) > 1 (Section 3.4).

### 3.3 Letters and their statistical dimensions [T]

Convention (Thm 2.1(i), Thm 1.5(iii)): a letter is a *family* of correspondences, the direct sum over its parameter set; its dimension is the number of summands (for unitary summands) or the sum of their fan-ins. The data selects the summand, as the Cuntz isometries select a summand in Thm 2.1. Dimensions are fixed by interface facts (palette 16, lattice 4096, |D4| = 8, `win_levels`) or by a **snapshot** |Ξ_0| := max(2, |Ξ|) of the shape alphabet taken at level start and at every promotion (T11: β_c is recomputed when the library changes; an observed new shape extends a family's domain but does not change its charged d until the next recompute). This is the [D] convention that keeps the Parry order stable within a level.

| letter | edge | d | derivation |
|---|---|---|---|
| `SelPointed` | F→C | 1 | the ACTION6 datum is a cell; it lies in one component (C2 lifts the datum to the base algebra) |
| `SelAll` | F→C | 1 | second quantization Γ: act on every component |
| `SelColour[c]` | F→C | 16 | selection by colour fibre (palette) |
| `SelRank[r]` | F→C | 16 | selection by HN rank (at most 16 distinct masses) |
| `SelType[τ]` | F→C | \|Ξ_0\| | selection by canonical type (T8) |
| `Step[g]` | C→C | 4 | ⊕_{g∈S} U_g: the gauge group by its generators (Component 1, §1.3) |
| `Move[v]` | C→C | 16128 | ⊕_{v∈[−63,63]²∖0} U_v: the translation group as a set |
| `Rot[g]` | C→C | 7 | ⊕_{g∈D4∖e} U_g |
| `Recolour[c']` | C→C | 15 | colour-fibre relabelling to another palette colour (C3; fan-in 1 or 2 per member, 15 members) |
| `Fill[c]` | C→C | 16 | the colour taken by cells a moved component vacates: T_{hξ}T_ξ\* annihilates ξ before creating hξ, and the vacated region must be assigned a fibre (palette) |
| `Copy` | C→C | 2 | ρ ⊕ ρ: creation of a second copy of an existing frame vector (Thm 1.5(iii)); followed by units for placement |
| `Annihilate[c]` | C→F | 16 | T_ξ\* with the vacated region's fibre as parameter (palette) |
| `Create[τ,g,t]` | F→C | \|Ξ_0\|·8·4096 | T_ξ from the vacuum: shape, orientation, position |
| `Bind[R,τ']` | C→P | 6·(\|Ξ_0\|+1) | adjoin a second argument standing in relation R ∈ {contact, overlap, h = e, h ∈ S, same colour, same shape} to the first; τ' a type or "any" |
| `Gate` | P→C | 2 | apply the body ⊕ identity according to the pair predicate: a 2-summand direct sum (T4 additivity, C3 "conditionals gated by the newborn residual") |
| `Const[Ψ\*]` | F→F | fan-in snapshot ≥ 2 | a constant correspondence; d = number of distinct pre-states it has been observed to map to Ψ\* |
| `LevelGate[l]` | F→F | `win_levels` | grade-0 gate on `levels_completed` |
| `Over` | F→⊤ | 1 | transition into the GAME_OVER flag; ⊤ is absorbing so no d = 1 cycle |

No letter is a game concept: each is a gauge element family, a Fock ladder operator, a colour-fibre map, or a gate on a base-algebra invariant (C3). The identity is the empty word (the vacuum of the word Fock space) and is not an edge. d = 1 edges: `SelPointed`, `SelAll`, `Over`; the subgraph they span (F→C, F→⊤) is acyclic. Every cycle through F passes an `Annihilate` (d = 16), a `Const` (d ≥ 2) or a `LevelGate` (d ≥ 6).

A **word** is a path with each family member fixed by data. A **clause** is `selector (F→C) · [Bind·Gate]? · body (C→C letters and fitted units) · terminator ∈ {stay, Annihilate[c]}`. A **creation clause** is `Create[τ,g,t]` or `SelType[τ]·Copy·units`. A **hypothesis** (ClauseSet) for an action id is a parallel composition of clauses acting on orthogonal summands of A_0 plus an optional `Over` and an optional `LevelGate` wrapper; its KMS weight is the product of its clauses' weights (the KMS state factorizes over orthogonal gauge slices, Prop. 1.3) and its dimension is the sum of its clauses' products (T4). Components selected by no clause are acted on by the identity.

### 3.4 The Parry–KMS prior at β_c [T]

A_β(u,v) = Σ_{e:u→v} d_e^{−β} on the 4×4 vertex matrix. r(A_β) is strictly decreasing (E6), r(A_0) > 1 (there are loops of multiplicity > 1), r(A_β) → 0 because the d = 1 subgraph is acyclic; so β_c with r(A_{β_c}) = 1 exists uniquely (Thm 2.1(iv)). The study's E3 correction is adopted: the cylinder weight v_{r(w)}Π d^{−β_c} is the state; the *conditional* kernel

    P(w | s) = (v_{r(w)} / v_{s(w)}) · Π_j d_{e_j}^{−β_c}

is the program prior given the start type, which is what a search from the frame type F uses (Prop. 9.3(i)). Item 3 read operationally (kms-library): at β_c the cylinder weights are Zipf with exponent 1, so the expected number of candidates examined before w in weight order is ≈ 1/P(w|s). **Search cost is inverse KMS weight**; the enumerator minimizes it by being best-first in P(w|s), and the promotion criterion (Section 8) is literally a code-length comparison in this currency.

```
def critical_state(lib):                                   # T5, T11; |V| = 4
    assert acyclic(lib.edges_with(d=1)) and perron_root(lib.A(0.0)) > 1
    lo, hi = 0.0, 64.0
    while hi - lo > 1e-9:
        mid = (lo+hi)/2
        lo, hi = (mid, hi) if perron_root(lib.A(mid)) > 1 else (lo, mid)
    beta_c = (lo+hi)/2; A = lib.A(beta_c)
    v = right_perron(A); u = left_perron(A); u /= u @ v; v /= v.sum()
    pi = {e: u[e.src]*v[e.dst]*e.d**(-beta_c) for e in lib.edges}     # stationary edge measure
    Lam = sum(pi[e]*log(e.d) for e in lib.edges)                        # Λ = −λ'(β_c), mean log-dimension
    b_fam = max(number of letter families leaving v for v in lib.V)    # branching ratio for the T10 ceiling
    return beta_c, u, v, pi, Lam, b_fam

def log_P(word, lib):                                      # conditional Parry kernel (E3 reading)
    return log(lib.v[word.dst]) - log(lib.v[word.src]) - lib.beta_c * sum(log(e.d) for e in word.letters)
```

Cost: 40 bisection steps of a 4×4 eigenvalue: microseconds. Recomputed at level start (snapshot refresh) and at every promotion, demotion or gauge reduction (T11), never otherwise.

### 3.5 The Bratteli-path data structure [T]

```
Letter    = namedtuple("Letter", "name src dst d apply fit abstract")
   # apply(obj, params, ctx) -> transformed obj(s) | None (annihilated) | TERMINAL
   # fit(pairs, ctx) -> one params shared across all pairs, or None (C11: moments, then exact)
   # abstract: for promoted letters, the defining word with type slots as variables (§8.1)
Word      = tuple[(letter_id, params)]            # a Bratteli path; vertices implied by src/dst chaining
Clause    = (selector: Word F→C, bind: Word C→P→C | None, body: Word C→C, term: Word C→F | None)
ClauseSet = (clauses: list[Clause], creations: list[Word F→C], over: bool, level_gate: int | None)
Hypothesis= (cs: ClauseSet, logw: float, n_verified: int, beyond_ceiling: bool)
Library   = (letters, edges, V, H ≤ D4, beta_c, u, v, pi, Lam, b_fam, Xi_snapshot,
             split_rules, bound_types, promoted: list[Letter], gauge_units: dict[type -> set])
```

---

## 4. The closed-loop solver

### 4.1 Stored transitions [T]

`Store[a]` holds every transition of action id a as `(Ψ_t, data, Ψ_{t+1}, flag_{t+1}, level_t)`, with Ψ_{t+1} the settled (last) frame. `Store∅[a]` holds the animation sub-steps (Section 4.8). Nothing unexplained is ever dropped (kms-library store policy); explained transitions may be summarized once a verified word reproduces them, but the frames are kept for the exact test. A transition is a demonstration pair of the environment's correspondence (T1) and is reused by every later synthesis, every level (Section 8.5).

### 4.2 Charge sectors first (Prop. 1.3) [T]

E_0 = ⊕_c E_c decomposes by colour fibre, so the torus T¹⁶ acts on the Fock space with multi-weight = component count per colour. The charge of a transition is n = count_{t+1} − count_t ∈ Z¹⁶; within colour c, creations minus annihilations equals n_c. Charge is the **grade-0 residual** and is computed before any matching: a hypothesis whose predicted charge differs from the observed charge is dead without rendering anything (R_0 ≠ ∅), and synthesis knows from n how many creations and annihilations its clause set must contain. Mass mismatch per colour and the flag mismatch complete R_0.

### 4.3 Matching: charge-0 moves by maximum KMS weight, transform-then-match [T]

A move is T_{hξ}T_ξ\* (eq. 1.3). The matching kernel is the Pimsner inner product ⟨ξ,η⟩_A, realized as cell overlap after transport; the assignment maximizes the KMS weight of the explaining word (T5), which is the same objective synthesis uses. Within a bucket of equal shape key (colour may differ: `Recolour`), the cost of pairing ξ ∈ Ψ_t with η ∈ Ψ_{t+1} is the negative log-weight of the cheapest body word mapping ξ to η:

    cost(ξ,η) = β_c·[ 1[v≠0]·min(|v|₁·log 4, log 16128) + 1[g_ξ≠g_η]·log 7 + 1[c_ξ≠c_η]·log 15 + 1[v≠0 or c_ξ≠c_η]·log 16 ]

(`Step^k` versus `Move`, `Rot`, `Recolour`, and the forced `Fill`). Unmatched components cost `Annihilate` (log 16) or `Create` (log(|Ξ_0|·8·4096)); a `Copy` of an existing type costs log 2 plus placement. Transform-then-match (C11, E7, item 2): v is read from first moments (bounding-box corners, linear in position), g from second moments (the squared variable: enumerate the 8 D4 images of Σ_ξ and keep those equal to Σ_η, then confirm on the canonical mask), the colour map from the fibres; the exact cell check comes last. Permutations of identical components are kept as alternatives *only* when they yield different fitted words (if every component of type τ moves by the same v, every permutation is the same word), which removes the combinatorial explosion two-projection admitted.

```
def match(Ψ0, Ψ1, lib):                                    # Prop 1.3, eq. 1.3, T5, C11
    charge = Ψ1.count_c - Ψ0.count_c
    pairs, left0, left1 = [], [], []
    for key in shape_keys(Ψ0) | shape_keys(Ψ1):
        A, B = Ψ0.by_key[key], Ψ1.by_key[key]
        C = full((len(A)+len(B), len(A)+len(B)), 0.0)      # square with dummy rows/cols for unmatched
        for i, ξ in enumerate(A):
            for j, η in enumerate(B): C[i,j] = pair_cost(ξ, η, lib)
            C[i, len(B):] = annihilate_cost(lib)
        C[len(A):, :len(B)] = create_cost(lib)
        rows, cols = linear_sum_assignment(C)
        for i, j in zip(rows, cols):
            if i < len(A) and j < len(B): pairs.append((A[i], B[j], gauge_element(A[i], B[j])))
            elif i < len(A): left0.append(A[i])
            elif j < len(B): left1.append(B[j])
    assert per_colour(len(left1) - len(left0)) == charge
    alts = distinct_word_permutations(pairs)              # keep only permutations that change the fitted word
    return GradedTransition(charge, pairs, created=left1, annihilated=left0, alts=alts,
                            flag=Ψ1.flag, pointed=Ψ0.cell2comp[data] if data else None)
```

Cost: a handful of assignments on ≤ 50×50 blocks, 1–3 ms.

### 4.4 Rendering a predicted Fock vector to a frame [D]

The exact test compares frames cell for cell, so a prediction must be a frame (the three designs that left this undefined are repaired here). Model: the HN flag is a filtration V_1 ⊂ … ⊂ V_r (Thm 3.1(iii)); the canvas realizes the filtration by painting strata in flag order, rarest last (on top). The colour of a cell a component vacates is **not** assumed to be the top stratum: it is the clause's `Fill[c]` parameter (d = 16), fitted from the first demonstration in which that clause vacates cells and verified exactly thereafter. Before any such demonstration, the unfixed `Fill` branch is a direct sum over the palette and contributes maximal disagreement to the Kato score (Section 5), so the agent learns the fill rule from one action rather than assuming it.

```
def render(Ψ, effects):                                    # effects: per-component fate from predict()
    canvas = Ψ.frame.copy()
    for comp, fate in effects:                             # 1. vacate
        if fate.kind in ('move', 'annihilate', 'recolour', 'reshape'):
            canvas.flat[comp.cells] = fate.fill            # Fill[c] parameter (fitted, d = 16)
    survivors = [(fate.new_comp) for comp, fate in effects if fate.kind != 'annihilate'] + effects.created
    for comp in sorted(survivors, key=lambda k: -k.mass_of_colour):   # HN flag order: rarest painted last
        canvas.flat[comp.cells] = comp.colour
    return canvas

def predict(hyp, Ψ, data, lib):                            # a ClauseSet applied to a Fock vector
    effects = []
    for cl in hyp.cs.clauses:
        for comp in select(cl.selector, Ψ, data):          # SelPointed / SelAll / SelColour / SelRank / SelType
            if cl.bind and not gate_open(cl.bind, comp, Ψ): continue      # pair predicate; else identity
            effects.append((comp, apply_body(cl.body, cl.term, comp, Ψ, lib)))
    effects.created = [apply_creation(w, Ψ, lib) for w in hyp.cs.creations]
    flag = 'GAME_OVER' if hyp.cs.over else 'NOT_FINISHED'
    return effects, flag
```

Where two predicted components claim the same cell the flag order decides; a word whose truth depends on a different overlap rule dies by T6 and is re-synthesized with a `Bind[overlap]·Gate` clause, so the data's own rule replaces the convention wherever it matters.

### 4.5 Exact verification (T6) [T]

A hypothesis survives a stored transition iff its predicted flag equals the observed flag and `array_equal(render(predict(h, Ψ_t, data)), Ψ_{t+1}.frame)`. The residual Hankel rank on `Store[a]` is zero iff this holds for every stored transition of a, with one parameter vector shared across them (a parameter that must change per transition is a table, not a word: T12). Nothing less is accepted. The graded newborn residual (T2) of a failing hypothesis is recorded for synthesis: R_0 (charge, mass, flag), then R_1 (matched components whose predicted gauge/colour differ; unpredicted creations/annihilations), then R_2 (the part of R_1 jointly explicable by a pair-conditioned correction); the newborn-ratio energy(R_k)/energy(R) (item 11) tells synthesis at which grade to start.

```
def residual_test(a, Ψ0, data, Ψ1, lib):
    for h in list(H[a]):
        effects, flag = predict(h, Ψ0, data, lib)
        if flag != Ψ1.flag or not array_equal(render(Ψ0, effects), Ψ1.frame):
            H[a].remove(h); Dead[a].append((h, graded_residual(effects, flag, Ψ0, Ψ1)))
    Store[a].append((Ψ0, data, Ψ1)); Ignorance[a][(Ψ0.hash, data)] = Ψ1.hash
    if not H[a]: Frontier[a] = Frontier[a] or new_cursor(a, lib); synthesize(a, lib, budget.slice('synth'))
```

### 4.6 Hypothesis elimination and synthesis in Parry–KMS order at β_c (T5, T11) [T]

Synthesis is a best-first enumeration of clause sets keyed by Σ log P(clause | F), resumable across actions (the heap is kept in `Frontier[a]`), organized by grade (T2), with the frontier pruned by the residual (kms-library): a clause is pushed only if its selector can select a component touched by the current unexplained residual, and only if its body's letter type can produce the observed change (a displacement admits `Step`/`Move`, a change of D4 class admits `Rot`, a colour change admits `Recolour`, a disappearance admits `Annihilate`, an appearance admits `Create` or `Copy`). Selector fitting is a *linear* test on invariant vectors (C11): a selector `SelX[p]` fits iff p is constant on the touched components and false on every untouched component of the same body behaviour; this is a pass over components, not a search over frames (the one-sided law of item 2 is why: matching the residual's indicator against invariants costs the number of invariants, not its square).

The grade ladder: find the clause set of grade 1 (selector·body) explaining the largest exact subset of `Store[a]`; the rest is the grade-1 newborn residual; search grade-2 clauses (selector·`Bind`·`Gate`·body) on that residual only; the complete hypothesis is the direct sum. Every exact clause set found within the budget is kept (T7); the per-action coherent set is ordered by weight and capped at K_max = 32 as a compute cap.

**The in-context depth ceiling as a trigger, not a refusal (T10).** L_max(T) = ⌊log(1+T)/log b_fam⌋ + 1 with T = |Store[a]| in component-level demonstration pairs and b_fam the number of letter *families* leaving C (8 initially), never the number of parameter edges. A clause set whose longest clause exceeds L_max is accepted if exact (T6 and T7 guard it) but flagged `beyond_ceiling`; its subwords are the first candidates in the maternity ward (Section 8.2). The ceiling never makes the agent refuse evidence it has; it tells the agent when it is relying on depth the fixed algebra cannot legitimately supply, which is exactly what item 7 says.

```
def synthesize(a, lib, budget):                            # T5, T11, T2, T6, T7, T12, C11
    cur = Frontier[a]; unexplained = cur.unexplained or Store[a]; found = []
    while budget.remaining() and cur.heap:
        logw, cs = cur.pop()                               # best-first in Σ log P(clause | F)
        params = fit_all(cs, unexplained, lib)             # one shared parameter vector; None if inconsistent
        if params is None: continue
        hyp = Hypothesis(cs.with(params), logw, 0, longest_clause(cs) > L_max(len(Store[a]), lib.b_fam))
        explained = [tr for tr in unexplained if exact(hyp, *tr, lib)]
        if len(explained) == len(unexplained):
            found.append(hyp)
            if len(found) >= K_MAX or (found and logw < found[0].logw - log(64)): break
            continue
        if explained and cur.grade_allows_partial():       # Wold: a partial exact word fixes its subset
            cur.accept_partial(hyp, explained); unexplained = [tr for tr in unexplained if tr not in explained]
            cur.restrict(unexplained); cur.next_grade_if_exhausted()
        for cs2 in expand(cs, unexplained, lib):           # residual-guided one-letter extensions
            cur.push(logw + delta_logP(cs2, cs, lib), cs2)
    if found:
        H[a] = assemble_direct_sums(cur.partials, found); Frontier[a] = None
        if any(h.beyond_ceiling for h in H[a]): maternity_ward(a, H[a], lib)        # §8.2
    else:
        cur.unexplained = unexplained; Frontier[a] = cur    # resume next step; until then Ignorance[a] is the only model
        if cur.exhausted(): force_promotion(a, cur.partials, lib)                    # T10, §8.2

def expand(cs, residual, lib):
    out = []
    for cl_idx in open_clauses(cs) + [NEW]:
        touched = residual_components(residual, cs)       # lifted to the base algebra (C11)
        for sel in selectors_fitting(touched, lib):       # linear test on invariant vectors
            for letter in letters_admitted_by(residual_change_type(touched), lib):   # residual-guided
                out.append(cs.extend(cl_idx, sel, letter))
        if cs.grade >= 2:
            for R, tau in bind_options(touched, lib): out.append(cs.gated(cl_idx, R, tau))
    return out
```

Cost: fitting is O(|Store[a]|·#components) per candidate (moments only), exact test O(|Store[a]|·4096) only for candidates whose fit succeeds; at ~50 µs per candidate the default slice of 60 ms examines ~1000 candidates, which under Zipf (item 3) reaches words of weight ≥ 10⁻³ of the top. Anything deeper resumes next step.

### 4.7 Memorizing versus generalizing (T12) [T]

`Ignorance[a]` is the lookup table (Ψ.hash, data) → Ψ'.hash: the compact part K of 0 → K → T_E → O_E → 0 (item 6). It never selects an action (C10). It has two legitimate uses: as a revisit detector in the planner's visited set, and as the **non-determinism detector**: if two stored transitions of a share (pre-state hash, data) but differ in post-state, no word with constant parameters can exist (T6); the action is marked `nondeterministic_on_observables` and its Kato score is zero for the rest of the level unless a grade-0 gate (`LevelGate`, a count predicate) separates the two cases, which synthesis tries first. This is the honest response to hidden state (N15): recorded, not papered over, and it cannot consume a level's actions.

### 4.8 Animation frames as serial compositions (T4) [T] (moment-map graft, repaired)

When `obs.frames` has k > 1 frames, the pair (frames[0] → frames[−1]) is the transition of action a and is what `Store[a]`, verification and planning use (the settled state; judge 3's `plan()` objection). The intermediate pairs (frames[i] → frames[i+1]) are stored in `Store∅[a]` as transitions under the null action ∅: Thm 1.5(iii) says the composite is the serial composition of its steps, so each step is a demonstration of a letter. Synthesis of the ∅-word for a runs in the same cursor with its own budget; if an exact ∅-word exists, (i) its letters' parameters are extra demonstrations for the a-word's fit (a `Step` seen five times in one animation pins `Step^5` versus `Move`), and (ii) the a-word may be *written* as `first_step ∘ ∅^{k−1}`, which is shorter in the enlarged alphabet after promotion. If the ∅-steps are inconsistent across animations of the same action, `Store∅[a]` is marked non-deterministic and ignored; the settled transition is unaffected. Intermediate frames therefore never kill a hypothesis about the settled transition (judge 2's objection) and are never discarded (N10).

```
def observe(self, obs, aid, data):                         # harness protocol; T1, T4
    frames = obs['frames']; Ψ1 = perceive(frames[-1], obs['levels_completed'], obs['state'], self.lib)
    if self.prev is None: self.prev = (Ψ1, obs); return
    Ψ0 = self.prev[0]
    if obs['levels_completed'] > Ψ0.level:                 # completion: a goal demonstration (§6.2)
        self.goals.positive(Ψ0, aid, data); self.on_level_end(Ψ1); self.prev = (Ψ1, obs); return
    self.budget.start_step()
    residual_test(aid, Ψ0, data, Ψ1, self.lib)             # §4.5–4.6 (includes GAME_OVER via the flag and Over)
    if len(frames) > 1:
        steps = [perceive(f, Ψ0.level, 'NOT_FINISHED', self.lib) for f in frames]
        for F0, F1 in zip(steps[:-1], steps[1:]): Store∅[aid].append((F0, None, F1))
        residual_test_null(aid, self.lib)
    self.goals.negative(Ψ1)                                # every non-winning post-state is a certified negative
    if obs['state'] == 'GAME_OVER': self.deaths.positive(Ψ0, aid, data)
    self.prev = (Ψ1, obs)
```

---

## 5. The Kato exploration rule as an algorithm with cost (T3)

**Derivation.** Thm 4.1(i): the minimum-norm interpolant of a target through a hypothesis subspace pays the tangent tax ‖tan Θ y‖², zero iff the subspace contains the target; normalized, the tax fraction is sin²Θ. Thm 4.1(ii) in the E10-corrected form: alignment changes at rate ≤ ‖G‖·‖[P,Q]‖_HS with ‖[P,Q]‖²_HS = ½ Tr sin² 2Θ, attained along the commutator. For two rank-one prediction projections at angle θ (cos θ = agreeing cells / 4096 in the flat one-hot metric; audit item 6), ½ sin² 2θ = 2 sin²θ cos²θ, which equals 2 sin²θ to first order and *decreases* past θ = π/4 (the Kato saddle two-projection patched by hand). The quantity the next observation will actually realize is the tangent tax: whichever hypothesis is true, every other hypothesis pays sin²θ against it and dies (T6). So the score is the KMS-expected tangent tax of the coherent set against itself,

    D(a, data) = Σ_{h<h'} p_h p_h' sin²θ(ŷ_h, ŷ_h')  =  ½(1 − Tr ρ_a²),   ρ_a = Σ_h p_h |ŷ_h⟩⟨ŷ_h|,

the purity deficit of the coherent prediction state (kms-library), monotone in disagreement for all θ, zero iff all predictions coincide, with p_h ∝ KMS weight (T7). Predictions are rendered as sparse diffs against Ψ_t and hashed; identical predictions merge into one class (two-projection), so the pairwise loop runs over distinct frames only.

**Candidates.** Every available simple action; for ACTION6, one representative cell per component (C2: the cells of one component are unitarily equivalent in M_{n_k}; the representative is the canonical first cell of the component's sorted index array, deterministic, N12), capped at 64 by HN order as a compute cap.

**Coherent set per candidate.** (i) An **untried** action (no stored transition this game) is tried once, in `available_actions` order, before any score is computed: with no demonstration every outcome is newborn (T2) and no ε-rule is involved (fock-ladder; judge 2's "every action is a demonstration"). (ii) An action with a nonempty H[a] uses H[a]. (iii) An action in the **ignorance state** uses the predictions of its frontier's accepted partial words (which act, C10; the table does not); if there are none, its score is zero unless the current context is **newborn** for that action: the pair (type multiset of Ψ_t, contact Gram by type, data's component type) has never appeared in `Store[a]`. In a newborn context the unfixed branches of the top-layer families (identity, `Step`, `Rot`, `Recolour`, `Annihilate`, `Create`, with `Fill` unfixed) are direct sums whose pairwise tax is maximal, so D = ½(1 − 1/m) with m the number of applicable families (moment-map's span projection, computed rather than hand-set). (iv) A `nondeterministic_on_observables` action scores zero (Section 4.7). Rules (iii)–(iv) are the fall-through judge 2 asked for, derived from T2 and T6 with no constant N.

**Cost (C8, brief A.3).** c(a, data) = 1 + P_death(a, data) · replay(Ψ_t), where P_death is the KMS-weight fraction of `Deaths` words (and `Over` clauses) predicting GAME_OVER for this candidate, and replay(Ψ_t) is the planner's path length from the level-initial Fock vector to Ψ_t under the verified model, or the number of actions since level start when no such path is found. Dying is cheap when the state is cheap to recover; it is never forbidden (audit item 20).

**Decision-relevant interleaving (T7).** If a goal target exists (Section 6.4), compute the plan to it under each surviving hypothesis combination that the plan touches; if all first actions agree, take that step (D = 0 there is not a reason to experiment); if they disagree, the candidate set is restricted to those first actions, and the disagreement is the experiment (two-projection).

```
def kato_select(Ψ, avail, lib, goal_target):               # T3 (Thm 4.1(i)-(ii), E10), T7, C5, C8
    for a in avail:
        if status[a] == 'untried': return (a, representative_data(a, Ψ))       # try once, no score
    cands = candidates(avail, Ψ)                           # simple actions + one cell per component (≤ 64)
    if goal_target is not None:
        firsts = {plan(Ψ, goal_target, Hc, lib)[0] for Hc in hypothesis_combinations(H, goal_target.actions)}
        if len(firsts) == 1: return firsts.pop()
        cands = [c for c in cands if c in firsts] or cands
    best, bestval = None, 0.0
    for (a, data) in cands:
        if status[a] == 'nondeterministic': continue
        hyps = H[a] or partial_words(Frontier[a])
        if not hyps:
            if not newborn_context(a, Ψ, data): continue
            m = len(applicable_families(Ψ, data, lib)); D = 0.5*(1 - 1/m)
        else:
            classes = {}                                   # hash(sparse diff) -> (weight, diff)
            for h in hyps:
                diff = sparse_diff(render(predict(h, Ψ, data, lib)), Ψ.frame)
                classes.setdefault(hash(diff), [0.0, diff])[0] += exp(h.logw)
            Z = sum(w for w,_ in classes.values()); D = 0.0
            for (wi, di), (wj, dj) in pairs(classes.values()):
                cos = 1 - len(symdiff(di, dj))/4096; D += (wi/Z)*(wj/Z)*(1 - cos*cos)
        val = D / cost(a, data, Ψ, lib)
        if val > bestval: best, bestval = (a, data), val
    return best                                            # None when every candidate has D = 0

def cost(a, data, Ψ, lib):
    p_death = sum(exp(g.logw) for g in Deaths if g.true_on(predict_consensus(a, Ψ, data))) / Z_deaths
    return 1 + p_death * replay_length(Ψ, lib)
```

**Complexity.** ≤ 64 candidates × ≤ 32 hypotheses sparse renders (O(cells of moved components) each) plus ≤ 64 × C(8,2) symmetric differences of index arrays: 10–30 ms typical, ≤ 60 ms worst; the candidate cap and the class merge bound it.

---

## 6. Goal inference as correspondence synthesis (C6)

### 6.1 Goal words and their quantifiers [T]

A goal is a predicate on the Fock space: a word from the same library read as an observable and lifted to the many-component space by **second quantization** (fock-ladder): for a grade-k predicate p, dΓ(p) counts the k-tuples satisfying p and Γ(p) asserts that all do. These are the canonical functorial lifts of a one-particle operator to F(E) and are the only quantifiers admitted; counting goals and universal goals are derived, not written.

| grade | goal family | d |
|---|---|---|
| 0 | `Zero[c]`: count_c = 0 (equivalently mass_c = 0) | 16 |
| 0 | `Eq[c,c']`: count_c = count_{c'} | 120 |
| 1 | `Absent[τ]` / `Exist[τ]`: dΓ(1_τ) = 0 / ≥ 1 | \|Ξ_0\| each |
| 1 | `Count[τ] = n`: dΓ(1_τ) = n, n fitted from positives | \|Ξ_0\|·4096 |
| 1 | `Fix[τ]`: τ has Hankel rank 0 under the verified words (Section 2.7) | \|Ξ_0\| |
| 2 | `Some[R,τ₁,τ₂]`: dΓ(p_R(τ₁,τ₂)) ≥ 1; `All[R,τ₁,τ₂]`: Γ(p_R(τ₁,τ₂)) | 6·\|Ξ_0\|² each |

R ranges over the six relations of Section 2.4 (gauge coincidence is C6's "injective pairing under a gauge element"; `Recolour` goals are `Exist[τ']` with τ' the recoloured type). Goal words are ranked by the same conditional Parry kernel at β_c as transition words (C6: low-dimension words first). Death words are the same families with GAME_OVER transitions as positives (type ⊤).

### 6.2 Demonstrations and the exact test [T]

**Negatives.** The environment evaluates the win after every action, so every observed post-state Ψ_{t+1} with no increment is a certified negative: any goal word true on it dies. **Positives.** The post-state of the completing action is never observed (the next frame is the next level), so g is evaluated on Φ_{a_T}(Ψ_T) predicted by the agreed H[a_T] (moment-map, fock-ladder); if H[a_T] disagrees or is empty, the positive is stored as *pending* and evaluated as soon as the coherent set agrees (the level-2 data usually settles it). Exact test (T6): g survives iff false on every negative and true on every confirmed positive. Because positives live on post-states and negatives on observed states, the two-projection contradiction cannot arise (audit item 21). Before the first completion the goal coherent set is initialized at level start from the Parry prior restricted to words **false on every frame observed so far** (a goal already satisfied cannot be the goal; kms-library), which repairs moment-map's empty `G_live` (audit item 22). Goal words are kept coherent (T7) and carried across levels; a word consistent with every completion so far ranks first (C6, C7).

```
def goal_update(G, Ψ_post_observed=None, positive=None, lib=None):
    if Ψ_post_observed is not None:                        # negative demonstration
        G[:] = [g for g in G if not g.eval(Ψ_post_observed)]
    if positive is not None:                               # (Ψ_T, a_T, data_T)
        if agrees(H[positive.a]):
            Ψ_hat = settle(predict(consensus(H[positive.a]), positive.Ψ, positive.data, lib))
            G[:] = [g for g in G if g.eval(Ψ_hat)]; Positives.append(Ψ_hat)
        else: Pending.append(positive)
    if not G: G[:] = regenerate_from_prior(Negatives, Positives, lib)   # consistent with all demonstrations
```

### 6.3 Pre-completion ranking by Halmos alignment (Component 4) [D]

Before any positive exists, the theory's only guide is the alignment of a goal word with the newborn residual stream (C6, Component 4). Model: the observables a goal word reads form a subspace Q_g of the observable space (masses, counts by colour and type, contact Gram entries, gauge positions by type); the observed transitions span a subspace Q_r of observables that have actually changed. A(g) := dim(Q_g ∩ Q_r)/dim Q_g, the fraction of g's inputs that any action has ever moved (Tr(P_g Q_r) in this model). Ranking is ω(g)·A(g): a goal whose inputs no action moves is deprioritized because nothing the agent can do tests it. No distance-to-target term exists (audit item 7).

### 6.4 Goal-splitting: Kato on the goal coherent state [T]

The goal coherent state is {g, ω_g}. For a reachable state s let W(s) = Σ_{g true at s} ω_g and W_tot = Σ ω_g. Arriving at s has two outcomes (win, or kill every g true at s), and the Kato quantity for a two-outcome split of weights W and W_tot − W is proportional to W(W_tot − W) (T3 applied to the goal projection). The target is the reachable state maximizing the lexicographic key (W(s)(W_tot − W(s)), W(s), −plan length): the first component halves the goal weight per probe when the state is mixed; when one goal word dominates, every state with it true ties at 0 on the first key and wins on the second, which *is* the probing loop (moment-map). On arrival without an increment, every g true at s dies: a whole slope class per probe (fock-ladder, moment-map). With a verified RESET word, states before the probe are cheap to return to.

```
def goal_split_target(Ψ, G, lib, budget):
    reach = plan_reachable(Ψ, lib, budget)                 # dict state -> path, under consensus words
    W_tot = sum(exp(g.logw) for g in G); best = None
    for s, path in reach.items():
        W = sum(exp(g.logw) for g in G if g.eval(s))
        key = (W*(W_tot - W), W, -len(path))
        if W > 0 and (best is None or key > best[0]): best = (key, s, path)
    return best[2] if best else None
```

### 6.5 The last fallthrough (C5, third clause) [D], flagged

If no goal state is reachable within the planning budget, the brief's rule applies verbatim: choose the action whose predicted outcome changes the goal observables (type counts, contact Gram, positions by type) in a way not yet seen in this level, ties broken by Kato disagreement. This is where judge 1's concern (a coverage rule in the object algebra) lives; it is reached only after Kato, interleaving and goal-splitting all return nothing, and every step it takes is still a demonstration. It is tagged [D] with the model "newborn information about the goal correspondence = a change in the goal observables not in the store" and listed as a falsifier in Section 13.

---

## 7. Planning over the verified model (T4, T13)

Planning runs only inside the model (N16). An abstract state is a Fock vector's component records (types, gauge elements, colours; static types frozen by Section 2.7); words are deterministic maps (`predict` without rendering); a path in the plan is a Bratteli path of vertices and letters (T13). The search is uniform-cost best-first with a node budget, death priced (never forbidden), the planner's visited set keyed by the abstract state's hash (the only hash use besides the ignorance table, N2). A verified RESET word (`Const`) is an ordinary edge. ACTION6 branching is one candidate per component, pruned to components whose selected clause predicts a change under the consensus word. When several hypotheses survive on an action the plan touches, the plan is computed under each (Section 5 interleaving). Plans are recomputed every step (exact models are cheap to query); any surprise during execution is a demonstration and returns control to Section 4.

```
def plan(Ψ, target, Hc, lib, node_budget):                 # Hc: one hypothesis per action id (consensus or a combination)
    frontier = [(0.0, 0, Ψ, [])]; seen = {abstract_hash(Ψ)}; expanded = 0; reach = {}
    while frontier and expanded < node_budget:
        cost, _, s, path = heappop(frontier); expanded += 1
        if target(s): return path
        reach[s] = path
        for (a, data) in actions_at(s, Hc, lib):           # simple actions + pruned ACTION6 cells + verified RESET
            s2 = settle(predict(Hc[a], s, data, lib))
            if s2.flag == 'GAME_OVER': c2 = cost + 1 + replay_length(s, lib)
            else: c2 = cost + 1
            h2 = abstract_hash(s2)
            if h2 in seen: continue
            seen.add(h2); heappush(frontier, (c2, expanded, s2, path + [(a, data)]))
    return None                                            # reach is returned by plan_reachable for §6.4
```

Budget: ≤ 1500 expansions per step at ~20 µs each (component-record arithmetic, no rendering) ≈ 30 ms; anytime: with no complete plan, goal-splitting consumes the reachable set already computed.

---

## 8. Library growth (E6, T5, T10) and cross-level carry-over (C7, T9)

### 8.1 The corrected promotion inequality [T] (kms-library, adopted verbatim)

Let the held explanations be H = {Ω_i} (one clause set per action id, every level, plus verified goal and death words), with code length at the critical state ℓ(H) = −Σ_i log P(Ω_i | s_i) = β_c Σ_i log D(Ω_i) − Σ_i log(v_{r_i}/v_{s_i}), D(Ω) the product of the dimensions of all letters in Ω. By item 3, ℓ(H) is the log of the expected enumeration cost of re-deriving H, which is what the next level will have to do with H's words as factors (T9).

A verified word w from type a to b with product dimension D_w = Π_j d_{e_j} and empirical fan-in d_new ≤ D_w (its Wold newborn projection as a single correspondence, Thm 1.4; measured as the many-to-one-ness of the realized composite on its observed domain) is a candidate letter e_new : a → b. Adjoining it perturbs the Perron root by δλ = u_a v_b d_new^{−β_c} =: π_new, and λ'(β) = −Σ_e π_e log d_e = −Λ < 0, so

    δβ_c = π_new / Λ > 0           (E6: β_c rises; compression raises the entropy rate).

Rewriting the n occurrences of w in H by e_new gives D(Ω_i') = D(Ω_i)(d_new/D_w)^{n_i}, and

    ℓ'(H') − ℓ(H) = − n β_c log(D_w / d_new) + δβ_c · L' − Δ_v,

with L' = Σ_i log D(Ω_i') the energy of the rewritten held set and Δ_v the Perron-vector shift. **Promote iff n β_c log(D_w/d_new) > (π_new/Λ) L' − Δ_v**: the bits saved by compression on all n uses must exceed the alphabet tax (the rise in β_c times everything still held). In code the exact form is used: recompute β_c', u', v' by bisection and compare Σ log P directly, so the sign error of Prop. 9.3(ii) cannot recur and no first-order approximation is trusted.

Corollaries that shape the system: (i) no compression, no promotion (d_new = D_w gives a zero left side); (ii) d_new = 1 means w is a bijection on its data: it is adjoined to the gauge group of its type, never to E (Definition 1.1's absorption of vacuum loops; Component 1's crossed product), which is how the gauge group grows by exactly the invertible correspondences the environment has been seen to apply; (iii) the tax scales with what is held, so late in a game only frequent or strongly compressing words are promoted; (iv) demotion across games is the same inequality reversed (C7). κ_lib = 0 is the only constant and is the theory's default. Promoted letters store their defining word with **type slots abstracted** (moment-map), so they instantiate on new types in later levels and games. Candidates: every subword of length ≥ 2 of every held explanation, evaluated at level end and whenever a new clause set is verified (so compression is available within a level).

```
def promotion_delta(lib, held, w, d_new):                  # exact form; the inequality above is its first order
    if d_new == 1: lib.gauge_units[w.src].add(w); return None      # corollary (ii)
    if d_new >= w.product_dim(): return +inf               # corollary (i)
    lib2 = lib.with_letter(Letter(f"P{len(lib.promoted)}", w.src, w.dst, d_new,
                                  apply=compose(w), fit=fit_compose(w), abstract=abstract_types(w)))
    lib2.refresh()                                         # critical_state: β_c', u', v' (T11)
    held2 = [om.rewrite(w, lib2.letters[-1]) for om in held]
    return -sum(log_P(om2, lib2) for om2 in held2) + sum(log_P(om, lib) for om in held)

def maternity_ward(lib, held, store):                      # T5 + E6, T10, T11, C7
    changed = True
    while changed:
        changed = False
        for w, n in sorted(subwords(held, min_len=2).items(), key=lambda kv: -kv[1]):
            d_new = empirical_fanin(w, store)              # Wold-projected many-to-one-ness on observed domain
            if promotion_delta(lib, held, w, d_new) < 0:
                lib.adjoin(w, d_new); lib.refresh(); held[:] = [om.rewrite(w, lib.letters[-1]) for om in held]
                changed = True; break
    return lib, held

def demote_pass(lib, held_new_game):                       # corollary (iv), at game start
    for e in list(lib.promoted):
        if -sum(log_P(om.expand(e), lib.without(e)) for om in held_new_game) < -sum(log_P(om, lib) for om in held_new_game):
            lib.remove(e); lib.refresh()
```

### 8.2 Promotion is forced at the ceiling (T10) [T]

When a verified clause set is `beyond_ceiling` (Section 4.6), or when synthesis exhausts its cursor with residual remaining, the subwords of the beyond-ceiling or partial words enter `maternity_ward` first. Even a marginal promotion lifts the ceiling because depth is counted in the enlarged alphabet. If no candidate passes, the episode's ceiling has been reached for this action: the partial words act (C10) and the table never does.

### 8.3 Bound states become particles (C2, T2, T13) [D]

A pair (ξ, η) whose relative gauge element h(ξ,η) and contact are conserved across every stored transition in which at least one of them moves (grade-2 newborn identically zero on `Store`) is a conserved grade-2 vector. Its promotion creates a new canonical type in Ξ (the union cell set with its own canonical form); perception thereafter emits it as one grade-1 vector, and clauses select it by type. The compressed word is the pair of clauses `SelType[τ_ξ]·body ⊕ SelType[τ_η]·body` → `SelType[τ_B]·body`, judged by the same inequality with d_new = the number of distinct bound types instantiated. Requires at least two transitions with motion; otherwise nothing is promoted.

### 8.4 Gauge reduction (T8) [D]

After each level, compute H' = {σ ∈ H : for every action a, conjugating every clause set in H[a] by σ (transforming its fitted units) yields a clause set in H[a'] for some action a'}; H' is a subgroup. Canonicalization uses H' thereafter: types that H' does not identify split into distinct types, and selections by type sharpen. The translation part is never reduced. This computes the brief's "gauge group = whatever the synthesized correspondences commute with". Stored Fock vectors are re-canonicalized lazily when read.

### 8.5 What persists [T]

Across levels: `Lib`, β_c, promoted letters, bound types, split rules, H', `H[a]` (as the initial coherent sets of level k+1, killed by the first residual test if wrong; kms-library), `Store[a]` (level-k transitions remain demonstrations unless a verified level-(k+1) word contradicts them, in which case that action's older transitions are retired and the agent records that its word is level-dependent, a `LevelGate` whose dimension `win_levels` the prior penalizes), `Goals` (confirmed first), `Deaths`, `Ignorance`. Across games in one run: `Lib` only, under `demote_pass` on the new game's first held set; nothing else is carried (C7, N11 inverted).

---

## 9. Action budget policy under RHAE and the 5× cap

Facts (brief A, A2): S_l = min(1.15, (h/a)²), every action counts including RESETs, later levels weigh more, h is unobservable, the gateway has no known action cap, the official API's 5h cutoff is not known to apply, abandonment is driven by the 9-hour clock, and a game left untouched scores zero.

Policy, all of it derived from those facts and the theory, none of it from h (N17):

1. **No budget knob refers to h or to a cap.** The agent cannot see h; the only way to score is to minimize expected actions per level, which is what Sections 4–7 do by construction: try-once (every action a demonstration), Kato (every experiment maximally eliminating), interleaving (no experiment when a decision does not need one), goal-splitting (each probe halves goal weight), promotion (shorter words next level). If a 5h cutoff exists, it changes nothing the agent can act on.
2. **Never stop acting.** Levels are never abandoned for action count; a slow level keeps the option value of unlocking heavier later levels. Abandonment happens only when the per-game wall clock the adapter is given runs out.
3. **Compute budgets** (the only constants, all compute, all anytime with carry-over): per-action soft target 0.2 s, hard cap 0.5 s (A2: 0.4 core-hours per game at one process per game); synthesis slice 60 ms (resumable cursor), Kato ≤ 60 ms (candidate cap 64, class merge), planning ≤ 1500 expansions, promotion at level end ≤ 50 ms per candidate pass, perception ≤ 15 ms. A 15-minute heartbeat guard forces an action out of the current coherent set if any loop overruns.
4. **Death pricing**, not death avoidance: cost 1 + P_death · replay (Section 5, 7). After GAME_OVER the harness or gateway submits RESET; it is counted and recorded as a demonstration of action 0.
5. **Floor policy under exception**: any uncaught error inside `act` falls back to Kato on whatever coherent sets exist, then to the try-once rule over `available_actions`; the error is logged (never a frame). Every game emits an action every turn.

---

## 10. Module list for the `nce` package

```
nce/
  __init__.py
  budget.py        class Budget: start_step(); slice(name) -> ctx; remaining() -> float; overrun() -> bool
  perceive.py      class Component; class FockVector
                   def perceive(frame: np.ndarray, level: int, flag: str, lib: Library) -> FockVector
                   def canonical(cells: np.ndarray, H: list[int]) -> tuple[bytes, int, tuple[int,int]]
                   def contact_matrix(comps: list[Component], shape) -> scipy.sparse.csr_matrix
                   class SplitRules: apply(cells) -> list[np.ndarray]; learn(comp, changed_cells)
                   class BoundTypes: bind(comps) -> list[Component]; add(pair_type)
  library.py       class Letter(NamedTuple); class Library
                   def initial_library(win_levels: int) -> Library
                   def critical_state(lib) -> tuple[float, np.ndarray, np.ndarray, dict, float, int]
                   def log_P(word: Word, lib: Library) -> float
                   Library.with_letter(letter) -> Library; Library.refresh(); Library.snapshot_alphabet(xi: int)
                   Library.assert_beta_c_exists()
  words.py         class Word; class Clause; class ClauseSet; class Hypothesis
                   def predict(hyp, Ψ, data, lib) -> tuple[Effects, str]
                   def render(Ψ, effects) -> np.ndarray
                   def settle(effects_and_flag, Ψ, lib) -> FockVector       # abstract post-state for planning/goals
                   def fit_all(cs, transitions, lib) -> Params | None       # C11: moments then exact
                   def exact(hyp, Ψ0, data, Ψ1, lib) -> bool
                   def graded_residual(effects, flag, Ψ0, Ψ1) -> Residual   # R0, R1, R2 + newborn ratio
  match.py         class GradedTransition
                   def match(Ψ0: FockVector, Ψ1: FockVector, data, lib) -> GradedTransition
                   def pair_cost(ξ, η, lib) -> float
  store.py         class Store: append(a, tr); transitions(a); null(a); retire_level(a, level); never_drop_unexplained()
                   class Ignorance: put(a, key, post_hash); nondeterministic(a) -> bool
  synth.py         class ParryCursor: push(logw, cs); pop(); restrict(unexplained); accept_partial(hyp, explained);
                                      grade_allows_partial(); next_grade_if_exhausted(); exhausted()
                   def synthesize(a, lib, store, H, frontier, budget) -> None
                   def expand(cs, residual, lib) -> list[ClauseSet]
                   def L_max(T: int, b_fam: int) -> int
  kato.py          def kato_select(Ψ, avail, H, frontier, store, deaths, lib, goal_target, budget) -> tuple[int, dict|None] | None
                   def cost(a, data, Ψ, deaths, lib) -> float
                   def newborn_context(a, Ψ, data, store) -> bool
  goals.py         class GoalWord: eval(Ψ) -> bool; observables() -> set
                   class GoalSet: negative(Ψ); positive(Ψ_T, a, data); resolve_pending(H); regenerate(lib); alignment(g, store)
                   def goal_split_target(Ψ, G, lib, budget) -> list[(int, dict|None)] | None
                   def goal_observable_newborn_step(Ψ, G, H, store, lib) -> tuple[int, dict|None]   # §6.5 fallthrough
  plan.py          def plan(Ψ, target, Hc, lib, node_budget) -> list | None
                   def plan_reachable(Ψ, lib, budget) -> dict[FockVector, list]
                   def replay_length(Ψ, lib) -> int
  promote.py       def empirical_fanin(w, store) -> int
                   def promotion_delta(lib, held, w, d_new) -> float | None
                   def maternity_ward(lib, held, store) -> tuple[Library, list]
                   def demote_pass(lib, held) -> Library
                   def bound_state_candidates(store) -> list[pair_type]
                   def gauge_reduce(lib, H_sets) -> list[int]
  agent.py         class NCEAgent:
                       def __init__(self, game_id: str, win_levels: int, available_actions: list[int])
                       def act(self, obs: dict) -> tuple[int, dict | None]
                       def observe(self, obs: dict, action_id: int, data: dict | None) -> None
                       def on_game_over(self, obs: dict) -> None
                       def on_level_end(self, Ψ_new: FockVector) -> None       # promotion, bound states, gauge, β_c
  kaggle_adapter.py class MyAgent(Agent)        # ARC-AGI-3-Agents template; translates FrameData <-> obs dict
  tests/            test_perceive.py test_library.py test_words.py test_match.py test_synth.py test_kato.py
                    test_goals.py test_plan.py test_promote.py test_agent_protocol.py synth_env.py
```

**Harness protocol** (exactly `harness.py`): `agent = NCEAgent(game_id, win_levels, available_actions)`; each turn `action_id, data = agent.act(obs)` with `obs = dict(frames=[np.int8 64×64, ...], state, levels_completed, win_levels, available_actions, full_reset)`; after the environment step `agent.observe(obs_new, action_id, data)`; on `state in ("GAME_OVER", "NOT_PLAYED")` the harness submits RESET itself and calls `agent.on_game_over(obs)` first; `observe` is then called with `action_id = 0`, which records RESET's demonstration. `data` for ACTION6 is `{"x": int, "y": int}` with x the column and y the row of the representative cell; `None` otherwise.

```
class NCEAgent:
    def act(self, obs):
        self.budget.start_step()
        try:
            Ψ = self.prev[0] if self.prev else perceive(obs['frames'][-1], obs['levels_completed'], obs['state'], self.lib)
            if self.level_start is None: self.level_start = Ψ
            self.goals.resolve_pending(self.H)
            target = goal_split_target(Ψ, self.goals.G, self.lib, self.budget.slice('plan'))
            choice = kato_select(Ψ, obs['available_actions'], self.H, self.frontier, self.store,
                                 self.deaths.G, self.lib, target, self.budget.slice('kato'))
            if choice is None and target: choice = target[0]
            if choice is None: choice = goal_observable_newborn_step(Ψ, self.goals.G, self.H, self.store, self.lib)
        except Exception as e:
            log(e); choice = floor_policy(obs)
        self.last = choice
        return choice[0], ({"x": choice[1][1], "y": choice[1][0]} if choice[0] == 6 else None)

    def on_level_end(self, Ψ_new):
        held = all_held_explanations(self.H, self.goals, self.deaths)
        self.lib, held = maternity_ward(self.lib, held, self.store)
        for b in bound_state_candidates(self.store): try_promote_bound(self.lib, b, held)
        self.lib.H = gauge_reduce(self.lib, self.H); self.lib.snapshot_alphabet(len(self.lib.Xi)); self.lib.refresh()
        self.level_start = Ψ_new; mark_all_actions_untried_only_if_no_H(self)
```

**Kaggle adapter.** `MyAgent(Agent)` from the ARC-AGI-3-Agents template implements `is_done(frames, latest_frame)` (true on `WIN`) and `choose_action(frames, latest_frame) -> GameAction`: it builds `obs` from `latest_frame.frame` (list of 64×64 grids), `latest_frame.state.name`, `.levels_completed`, `.win_levels`, `[a.value for a in latest_frame.available_actions]`, `.full_reset`; calls `nce.observe(obs, prev_action_id, prev_data)` when a previous action exists; returns `GameAction.RESET` when the state is `NOT_PLAYED` or `GAME_OVER` (calling `on_game_over`), else maps `act`'s `(id, data)` to `GameAction.from_id(id)` with `.set_data(data)` for ACTION6. One `NCEAgent` per game; the library object is shared across games in the same process under `demote_pass`. The adapter enforces the per-game wall clock and the 15-minute heartbeat (Section 9). Exact class and method names are checked against the template at packaging time; no other code depends on them.

---

## 11. Test plan

**Unit tests per module**, all on synthetic frames generated by the test suite (random masks, random colours, known gauge elements; never a game frame):

- `perceive`: HN flag equals colours by mass with ties grouped on 10³ random histograms; 4-connectivity (diagonal-only neighbours are distinct components); canonical key invariant under all of D4 × Z² on 10⁴ random shapes and distinct for non-equivalent shapes; split rules and bound types applied; `perceive < 15 ms` on 300-component frames.
- `library`: β_c exists and is unique on the initial library and on 10³ random libraries with acyclic d = 1 subgraphs; the assertion fires on a d = 1 cycle; r(A_β) is decreasing (E6 check); adjoining any edge raises β_c; `log_P` is a probability over words of fixed length from F.
- `words`/`match`: on synthetic transitions with known moves, rotations, recolourings, creations, annihilations and fills, `match` recovers the decomposition exactly and `fit_all` recovers the parameters from one demonstration; charge conservation asserted; `render ∘ predict` reproduces the generated post-frame; permutation ambiguity yields alternatives only when fitted words differ.
- `synth`: synthetic environments whose mechanics are clause sets *sampled from the library's own Parry prior* (grade 1 and grade 2, with gates on contact/coincidence); measure demonstrations-to-verification: grade-1 words from 1 demonstration, gated words from ≤ 4 (one open, one closed gate); the resumable cursor gives identical results split across steps; `beyond_ceiling` flag set correctly; residual-guided pruning never prunes the true word.
- `kato`: with two hypotheses that differ on exactly one action, that action is chosen in 100% of 10³ trials; untried actions are chosen first; non-deterministic actions score zero; the ignorance fall-through never selects the same (action, context) twice; interleaving returns the agreed first step when plans agree.
- `goals`: positives evaluated on the predicted post-state; a synthetic win predicate from each family is recovered after one completion and survives level 2; goal-splitting halves goal weight per probe on a synthetic 3-level environment; the pre-completion filter removes satisfied words.
- `plan`: BFS over a synthetic model finds the shortest path with death priced not forbidden; RESET edge used when shortest.
- `promote`: numerical verification of δβ_c = π_new/Λ against bisection on random libraries; corollaries (i)–(iv) hold; a synthetic level-3 whose word is a composite of level-1 words is solved under the depth trigger only with promotion enabled (the T10 test); demotion removes letters unused in a new synthetic game.
- `agent_protocol`: a fake harness drives `act`/`observe`/`on_game_over` with the exact `obs` dict and asserts per-step wall clock < 0.5 s over 2000 steps.

**Blind scoreboard.** `python harness.py nce.agent:NCEAgent --procs 8 --time-per-game 1200 --out run.json`, reading only `levels`, `level_actions`, `baselines` and `score` from the output. Diagnostics logged per game, all agent-internal and permitted: actions until every available action has a nonempty H[a]; actions in the ignorance fall-through; coherent-set sizes; synthesis steps resumed; promotions and demotions; bound-state and gauge-reduction events; goal words alive per level; per-step time. Triage is by diagnostic class reproduced in synthetic environments, never by inspecting a frame.

**Theory-predicted ablations** on the scoreboard (each toggled by a flag): no promotion (T10 predicts later levels fall); no Kato (random among available actions; T3 predicts more actions to agreement); last-frame only (T4 predicts longer words and slower verification where animations exist); no interleaving (T7 predicts extra experiments); no goal-splitting (probing only). An ablation that does *not* move the scoreboard in the predicted direction is a falsifier (Section 13).

---

## 12. Two-week build plan with milestones and cuts

The calendar in A2 is eight days to the final submission, so the plan is written for fourteen days with the eight-day cut line marked; items after the line are stretch.

| day | milestone | acceptance |
|---|---|---|
| 1 | `perceive`, `canonical`, contacts, `FockVector`; `library` with β_c, assertions, snapshot | perception tests pass; β_c tests pass |
| 2 | `words` (`predict`, `render`, `fit_all`, `exact`), `match`, `store`, `Ignorance` | synthetic transition suites recovered from one demonstration |
| 3 | `synth`: cursor, residual-guided `expand`, grade ladder, partials, `beyond_ceiling`; `agent` skeleton with try-once and the harness protocol | first blind scoreboard run (levels/actions only) |
| 4 | `kato`: purity-deficit score, candidates, cost with replay, newborn-context rule, non-determinism rule; animation ∅-store | Kato tests; second scoreboard run |
| 5 | `goals`: families, lifts, negatives/positives on post-state, pending, pre-completion filter, alignment; `plan`; goal-splitting; interleaving | synthetic 3-level suite; third scoreboard run |
| 6 | `promote`: exact inequality, corollaries, maternity ward at level end and on verification, forced promotion, type-slot abstraction, cross-level carry, demotion | T10 synthetic test passes; fourth run |
| 7 | Budgets everywhere, heartbeat, floor policy, profiling to 0.2 s typical / 0.5 s hard; Kaggle adapter; offline packaging; license | 2000-step protocol test under cap; dry-run notebook |
| 8 | **Cut line.** Freeze, adversarial code review against Section 0's table and N1–N17, first submission | submission made |
| 9–10 | Bound-state promotion (8.3); stratum split rules (2.5 downward); gauge reduction (8.4) | synthetic multi-colour and rotation suites |
| 11 | `LevelGate` retirement logic; `Const` and planning through RESET; death words in planning | synthetic suites |
| 12 | Ablation runs (Section 11); triage by diagnostics | ablation table |
| 13 | Hankel-rank classification (2.7) feeding the planner's frozen types and `Fix[τ]` | tests |
| 14 | Buffer; second submission pick | — |

**What is cut if time runs out, in order:** gauge reduction (types stay D4-canonical; cost: more hypotheses, no wrong ones), bound-state promotion (multi-colour things stay grade-2; cost: longer words), stratum split rules (sub-component mechanics stay unexplained; cost: ignorance on those actions), Hankel-rank freezing (planner slower), `LevelGate` retirement (level-dependent actions re-synthesized per level). Never cut: exact verification, the coherent set, Kato, goal positives on the post-state, the promotion inequality, the budgets, the floor policy.

---

## 13. Honest expected failure modes and what would falsify the design

**Failure modes.**

1. *Hidden state.* Mechanics depending on unobserved quantities admit no zero-residual word. The agent detects it (Section 4.7), stops spending Kato actions on it, and acts on partial words and grade-0 gates; levels whose win depends on such an action are reached only by goal-splitting through other actions, often not at all.
2. *Non-determinism.* Exact verification kills every hypothesis under stochastic transitions; there is no approximate-fit actor by principle (T6, T12). Expected to score near zero on stochastic environments.
3. *Vocabulary of goals.* Goals outside the second-quantized families (order of events, grade ≥ 3 relations, arithmetic) are never confirmed; such levels are completed only by goal-splitting over the expressible words or by the C5 fallthrough, at several multiples of h.
4. *Strata misaligned with causal units.* A multi-colour or diagonal thing is several grade-1 vectors until bound-state promotion (which needs two transitions with motion); until then explanations are longer and disambiguation costs actions. Parts moving differently stay at grade 2.
5. *Fill rule richer than a constant colour.* If vacated cells reveal what was occluded (memory of the canvas), `Fill[c]` with constant c dies; the action falls to ignorance until a `Bind[overlap]·Gate` clause or promotion expresses it. This is the price of not assuming a background.
6. *Click geometry inside a component.* `SelPointed` cannot distinguish cells within one component; if the response depends on the cell, all clause sets for ACTION6 die and the action is in ignorance until split rules (2.5) isolate the responsive cells, which requires the residual to touch a strict subset of the component.
7. *RHAE floor.* Levels with h ≈ 6–10 cannot absorb try-once plus a few Kato steps; near zero on them is accepted.
8. *Budget blindness.* Without h the agent may run goal-splitting long on a level whose goal vocabulary is inadequate; the only brake is the clock.
9. *Promotion stalls.* If the inequality rejects every candidate while the ceiling binds, the episode cannot deepen (T10) and compositional later levels are lost.
10. *Perception caps.* Frames with many hundreds of components stay under budget but the ACTION6 candidate cap (64 by HN order) can hide the right click; coherent-set caps can drop a true hypothesis on textured frames.

**Falsifiers** (each measured by the Section 11 diagnostics, none by a frame).

- If, on most games, the goal word confirmed at the first completion is killed early in level 2, then goal semantics on the predicted post-state or the goal families are wrong (C6 as read here fails).
- If actions-to-agreement (every available action with a nonempty H[a]) typically exceeds ≈ 3·|available_actions|, the premise that a grade-1 word is fixed by one demonstration and verified by the rest is false for this benchmark, and the one-sided law (item 2) does not buy what the design assumes.
- If the no-Kato ablation does not increase actions-to-agreement, T3 as an exploration criterion is doing no work.
- If promotions fire but synthesis time and word length on later levels do not fall, or if later-level completions do not depend on promotion in the ablation, T10's mechanism is not the bottleneck the design says it is.
- If more than half of a typical game's actions are spent in ignorance fall-through or the C5 fallthrough, exact verification is too brittle for these environments and the T6 standard itself is what fails.
- If the last-frame-only ablation scores the same, the serial-composition reading of animations (T4) carries no information here.
- If gauge reduction ever shrinks H to a subgroup under which later-verified words are *not* covariant, the "gauge = commutant of the synthesized words" reading of T8 is unstable and must be dropped.

Each falsifier names the theorem whose operational reading would be refuted; none can be repaired by adding a game-specific rule without violating the governing rule, which is the point of listing them.
