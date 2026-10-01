# Foundations, part 1: the 1 Oct input, verified

*Fresh slate, 1 Oct 2026. This note checks every factual claim of the user's input of 1 Oct, [input-2026-10-01-matchings-bethe-cluster.md](input-2026-10-01-matchings-bethe-cluster.md), against sources opened in this session or against proofs given here. It then develops the input's three ideas A, B and C into precise statements. Contract: [BRIEF.md](BRIEF.md). Part 2, [foundations-unlocks.md](foundations-unlocks.md), uses some of these claims (its Appendix A); §7 lists where this note refines it. The BRIEF and Part 2 call this note `foundations.md`.*

---

## Summary: the five most important statements

1. **The input's sections 1–4 hold up; five corrections matter for design.** Of the 50 factual claims there, 38 are CORRECT, 11 IMPRECISE, 1 UNVERIFIED and none INCORRECT.
   - Permanent and determinant coincide only as *formal* tropicalisations. The zero-temperature limit of a determinant keeps its signs and can cancel completely (C15).
   - Kasteleyn's face rule is the flux that maximises |det| on every planar bipartite graph (Lieb–Loss). It minimises the half-filled ground-state energy only on the lattices covered by Lieb's theorem; in general it does not (Lieb–Loss's counterexample, reproduced in C12).
   - MTP₂ faithfulness needs the graphoid hypothesis. MTP₂ is also a different positivity from the matrix total positivity where cluster algebras began.
   - The Jerrum–Sinclair–Vigoda chain crosses to finite temperature for *every* nonnegative matrix, through the conductance of its state space, not through expansion of the graph.
   - Bethe gets counts right only to a factor e^{O(n)}: per_B ≤ per ≤ 2^{n/2} per_B. Bethe is a per-site (free-energy) object.
2. **Idea A, block form** (THEOREM, proved in §2, checked C1–C4, C18). For vector-valued a and c:
   - det M[aBc] · det M[B] = det M[aB] · det M[Bc] · ∏ᵢ (1 − ρᵢ²), where ρᵢ are the canonical correlations of a and c conditional on B. Hence I(a : c | B) = −½ Σᵢ log(1 − ρᵢ²).
   - The input's cross-minor ratio generalises to ∏ᵢ ρᵢ², which can vanish without Markovianity.
   - Markov ⟺ rank M[aB|Bc] = |B| ⟺ all |a|·|c| scalar exchange relations degenerate at once.
   - The top-k conditional canonical variates are the optimal rank-k Gaussian memory, with residual −½ Σ_{i>k} log(1 − ρᵢ²).
3. **Idea A, structure.**
   - For Gaussians, "a Markov network is the simultaneous degeneration of exchange relations" is a theorem.
   - Kenyon–Pemantle's cluster-type structure on principal and almost-principal minors goes further (hexahedron move = six cluster mutations):
     - a cube move replaces the conditioning sets of three partial correlations ("mutation = change of separator", DERIVED);
     - a Gaussian graphical model is a coordinate stratum of a chart iff some chart's conditioning sets separate every non-edge (DERIVED);
     - by exhaustive search over all charts for n ≤ 6, this covers every chordal graph on ≤ 6 vertices and every connected graph on ≤ 6 vertices except the triangular prism (C21).
   - MTP₂ is not the positive part: it is not an orthant of a chart (C21), and Chepuri–George–Speyer Thm 5.1 rules out a sign match with Grassmannian positivity.
   - The identity does not extend to quasi-free fermions (C19).
   - The general question stays OPEN, now posed precisely (§2.3).
4. **Idea B: the classical half is a set of theorems with explicit hypotheses, which the input states too loosely.**
   - Bethe free energies are exact on locally tree-like graphs for ferromagnetic Ising and for bipartite two-spin systems *without* correlation decay (Dembo–Montanari; Sly–Sun).
   - Marginals and efficient algorithms *do* need uniqueness (Weitz; Sly–Sun). The competition scores marginals.
   - Perfect matchings need expansion or large girth (ACFK; Gamarnik–Katz).
   - Separator bounds become vacuous on expanders for regions of size ≳ log n.
   - The quantum half is stated as Conjecture QB, with a test (§3.4).
   - Box spaces: expander ⟺ spectral gap of the profinite action (proved, §3.3).
5. **Idea C is well defined but has no route to the competition.**
   - It is well defined as automorphisms of the quantum torus and as unitary operators on L²(ℝ^I) built from Faddeev's Φ_b. Each layer is a set of commuting local gates.
   - It returns exactly, up to a phase that is 1 with tropical signs (Kashaev–Nakanishi).
   - Its half period acts as a permutation: Chin's theorem, checked on A_m × A_n T-systems in exact arithmetic (C17).
   - It is not a C*-automorphism of Connes' torus, since 1 + U is not invertible there.
   - Its dilogarithm value is (π²/6)·N₋, an integer multiple of π²/6, not a central charge.

**Recommendation (§5).**
- For the competition, push A: an exact KL and memory instrument that can be run at n = 1024 within hours.
- For the programme, develop B.
- Park C.

**Pointers for the six design streams.**
- [F] Faces: §1.5. A zero-temperature (frozen-face) limit keeps sign cancellations (C15), consistent with Part 2 unlock 29.
- [B] Bethe:
  - §1.1: Bethe is exact only per site; the cover limit is a lim sup.
  - §3: free energy versus marginals; uniqueness is needed for marginals.
- [M] Signings: §1.2–1.4 (C8–C12, C13–C14).
- [T] Tropical: §1.5 and §4 (C15–C17).
- [H] Heisenberg: §3.4 and §4. Idea C is a Heisenberg-picture object.
- [K] Markov/CMI: §2 (Theorems A and A′, C1–C7, C18–C21) and §3.1.

---

## 0. How to read this note

**Verdicts on the input's claims.**
- CORRECT: true as stated, with a source (theorem number) opened in this session or a proof here.
- IMPRECISE: true after a correction, which is given as **Fix**.
- INCORRECT: false; the counterexample or the correct statement is given.
- UNVERIFIED: neither a retrieved source nor a proof here. This includes references to conversations we do not have.

**Labels on statements developed here** (as in Part 2).
- THEOREM: published, source opened.
- DERIVED: proved here. "Checked" means a numerical check in §6.
- CONJECTURE: precise, with evidence and a test.
- OPEN: a precise question.
- SPECULATION: a direction not yet precise.
- "(from memory)": a standard fact not re-opened in this session. No verdict rests on such a fact alone.

**Numbering.**
- The input's claims are I1–I79 in order of appearance. Appendix A lists every claim with its verdict and streams.
- Numerical checks are C1–C21 (§6). They use numpy only, with OPENBLAS_NUM_THREADS=1, fixed seeds and small sizes.

**Stream keys** (as in Part 2): [F] faces, [B] Bethe/cavity, [M] matchings/signings, [T] tropical, [H] Heisenberg, [K] Markov/CMI.

**Notation.**
- M is a symmetric positive definite (PD) matrix. M[X] is the principal submatrix on X, and M[X|Y] has rows X and columns Y, in the order written.
- p_K = det M[K] is a principal minor. a_{ij|K} = det M[iK|jK] is an almost-principal minor, with K listed in the same order after i and after j.
- ρ_{ij|K} is the Gaussian partial correlation. per, haf and Pf are the permanent, hafnian and Pfaffian; pm(G) is the number of perfect matchings of G.

---

## 1. Claim-by-claim verdicts (input sections 1–4)

### 1.1 Matchings, Bethe and locality [B, M]

**I1.** "The monomer–dimer partition function has no phase transition (Heilmann–Lieb: the matching polynomial has only real roots)." **CORRECT.**
- Heilmann–Lieb (CMP 25 (1972) 190–232), abstract: the partition function "has its zeros on the imaginary axis when the dimer activities are nonnegative", so "no monomer-dimer system can have a phase transition … except, possibly, when the monomer density is minimal (x=0)".
- Equivalently, the matching polynomial is real-rooted (ACFK Thm 3.3; MSS I Thms 3.1–3.2).
- The excluded point x = 0 is exactly the perfect-matching regime of I3–I8.

**I2.** "So its growth rate is fixed by local neighbourhood statistics for every bounded-degree graph (ACFK)." **CORRECT**, read as estimability.
- ACFK Thm 1.2: along every Benjamini–Schramm-convergent sequence of graphs of bounded degree, ln M(G_n)/v(G_n) converges. The same holds for the monomer–dimer free energy at every positive activity.
- Uniform form (proof here). Suppose two sequences of graphs of maximum degree d have local statistics converging to the same limit but values ln M/v that stay ε apart. Interleaving them gives a convergent sequence along which ln M/v does not converge, a contradiction. So for every ε there is a radius R(ε, d) such that radius-R statistics determine ln M/v to ±ε, uniformly over such graphs.
- The mechanism is the one the input names. Real-rootedness gives ln M/v = ∫ ½ ln(1 + x²) dρ_G(x) for the matching measure ρ_G, whose moments count closed tree-like walks (ACFK Remark 3.6).
- Check (C9). On large-girth 3-regular graphs ρ_G tends to the Kesten–McKay law, and ∫ ½ ln(1+x²) dKM₃ = 0.58157540 = ½ ln(16/5). This equals the tree (Bethe) value ln(1 + 3r) − (3/2) ln(1 + r²) with r = ½ to all printed digits.

**I3–I5.** "Perfect matchings: not local in general … one edge lies in all but an exponentially small fraction of perfect matchings … a global rigidity no local view can detect." **CORRECT.**
- ACFK Thm 1.7: there are d-regular bipartite graphs with an edge e such that p(e) > 1 − cⁿ, where p(e) is the fraction of perfect matchings containing e.
- ACFK Thm 1.8: ln pm(G)/v(G) is not estimable, even among d-regular bipartite graphs.

**I6–I7.** "On bipartite expanders, no edge probability can be exponentially close to 1, and the growth rate becomes locally determined (via Gamarnik–Katz)." **CORRECT** for d-regular bipartite δ-expanders.
- ACFK Thm 1.9: every edge has p(e) ≥ (1/d) n^{−2 ln(d−1)/ln(1+δ)}. The d edges at a vertex have probabilities summing to 1, so no p(e) is within an exponentially small distance of 1.
- ACFK §1.3: estimability for such graphs "can be deduced from Corollary 1 of Gamarnik–Katz".
- Gamarnik–Katz (math/0702039), Cor 1: 1 ≤ Z(λ,G)/(λⁿ per) ≤ exp(O(n λ⁻¹ log⁻¹(1+α) log Δ)) on α-expanders. The perfect-matching count is therefore a uniformly controlled λ → ∞ limit of estimable monomer–dimer quantities.
- Gamarnik–Katz Thm 2: a deterministic (1+ε)ⁿ approximation in polynomial time for constant Δ and α.

**I8.** "On large-girth d-regular bipartite graphs it equals Schrijver's constant ((d−1)^(d−1)/d^(d−2))^(1/2) per vertex." **CORRECT.**
- ACFK Thm 1.5: ln pm(G)/v(G) → ½ ln((d−1)^{d−1}/d^{d−2}) as the girth tends to infinity.
- Large girth alone suffices; no expansion is assumed.
- Check (C9): ∫ ln|x| dKM_d matches ½ ln((d−1)^{d−1}/d^{d−2}) to 3e-6 for d = 3, 4, within the quadrature error at the log singularity.

**I9–I10.** "That constant is the Bethe permanent. Pretend the graph is a tree and use the exact Markov factorization valid on trees." **CORRECT.**
- Definitions (Vontobel Cor 15; Anari–Rezaei Def 2):
  - per_B(A) = exp(−min_P F_B(P)) over doubly stochastic P;
  - F_B(P) = Σ P log(P/A) − Σ (1−P) log(1−P).
  - F_B is convex (Vontobel Thm 20 and Lemma 21: the Bethe entropy is ½Σ_rows S + ½Σ_cols S with S concave).
- Proof for the biadjacency matrix A of any d-regular bipartite graph:
  - P = A/d is doubly stochastic.
  - On its support it satisfies the stationarity condition log(P(1−P)/A) = αᵢ + β_j.
  - Convexity makes P = A/d the minimiser, so per_B(A) = ((d−1)^{d−1}/d^{d−2})ⁿ for every girth (Vontobel §VII-E carries out the same evaluation).
- Checked to 10 digits on a 3-regular graph with n = 10 and a 4-regular graph with n = 9 (C9).
- Consequences:
  - Schrijver's bound per ≥ ((d−1)^{d−1}/d^{d−2})ⁿ (ACFK Thm 1.6) is Gurvits's per ≥ per_B, specialised to regular graphs.
  - Bethe's free energy is the tree (junction-tree) formula: edge terms minus (deg − 1) times vertex terms, applied as if the graph had no cycles.

**I11.** "For every nonnegative matrix, per_B ≤ per ≤ 2^(n/2) · per_B (Gurvits; Anari–Rezaei)." **CORRECT.**
- Anari–Rezaei Thm 3 (Gurvits): per ≥ bethe. Thm 4: per ≤ √2ⁿ bethe, tight for I_{n/2} ⊗ J₂.
- Check (C9): over 40 random nonnegative matrices, min per/per_B = 1.000000 and max per/(2^{n/2} per_B) = 0.85; I_{n/2} ⊗ J₂ gives per/per_B = √2ⁿ exactly.

**I12.** "The Bethe permanent is also a limit over graph covers (Vontobel)." **IMPRECISE.**
- **Fix:** Vontobel Thm 39 states per_B = lim sup_{M→∞} (average permanent over M-covers)^{1/M}, a lim sup.

**I13.** "On pseudorandom geometry, a global count is computed by a local, tree-shaped Markov factorization." **IMPRECISE.**
- **Fix:** what the tree computes is the per-site growth rate (the free energy density), not the count. Even for a d-regular bipartite graph, per/per_B ranges over [1, 2^{n/2}] (I11); for example 70/17.76 for the 3-regular graph of C9.
- "Without that, global constraints can freeze it" is CORRECT for perfect matchings (I3–I5).

### 1.2 Signs and topology [M]

**I14.** "Kasteleyn: on a planar graph, signing the edges turns the permanent into a determinant." **IMPRECISE.**
- **Fix:** for a *bipartite* planar graph, a signing of the biadjacency matrix gives pm = per = |det|. For a general planar graph, an orientation turns the hafnian into a Pfaffian: pm = |Pf|.
- Lieb–Loss's Thm 3.1 reproves the planar case.

**I15.** "On a genus-g surface you need 4^g determinants, one per spin structure." **CORRECT** (Pfaffians, or determinants in the bipartite case).
- Cimasoni–Reshetikhin (math-ph/0608070): Z = 2^{−g} Σ_ξ (−1)^{Arf(ξ)} Pf(A^{K_ξ}), summed over the 2^{2g} equivalence classes of Kasteleyn orientations. These classes correspond to spin structures, and the Arf invariants give the signs.
- The 4^g count is Kasteleyn's; Galluccio–Loebl and Tesler proved it (as cited there).

**I16.** "The obstruction is the surface's first Z/2-cohomology." **IMPRECISE.**
- **Fix:** H¹(Σ; Z/2) ≅ (Z/2)^{2g} acts freely and transitively on the classes of Kasteleyn orientations (Cimasoni–Reshetikhin Thm 3.2). It *indexes* the 4^g terms; it is not an obstruction class.
- The obstruction to a single Pfaffian is that no one orientation is Kasteleyn-correct for cycles in every homology class.

**I17.** "Robertson–Seymour–Thomas (and McCuaig) characterized the bipartite graphs that admit such signs: planar pieces glued together, plus one exception." **CORRECT.**
- RST (Ann. of Math. 150 (1999); math/9911268), Thm 1.3 = 6.8: a brace has a Pfaffian orientation iff it is isomorphic to the Heawood graph or is obtained from planar braces by repeated 4-sums (trisums).
- General bipartite graphs reduce to braces through 0-, 2- and 4-sums along tight cuts (§7.2). Recognition is polynomial (O(n³)).
- Equivalent problems: Pólya's permanent problem and even directed circuits (Vazirani–Yannakakis, §1.1). Little's K_{3,3} criterion is quoted in §1.2.
- McCuaig's independent proof is from memory.

**I18.** "The exception is the Heawood graph, the incidence graph of the Fano plane, whose nontrivial eigenvalues ±√2 beat the Ramanujan bound." **CORRECT.**
- Check (C8):
  - spectrum ±3 (simple) and ±√2 (multiplicity 6 each), against 2√2 = 2.83;
  - per(N) = 24 perfect matchings, and some signing reaches |det| = 24 (a Pfaffian orientation; RST §6.3 gives one);
  - K_{3,3}: per = 6, and every signing has |det| ≤ 4 (Little).

### 1.3 Signs and physics [M, H]

**I19–I21.** "Determinants are fermions; permanents are bosons. Free-fermion (matchgate) computations are classically easy, while boson sampling is believed hard." **CORRECT.**
- Terhal–DiVincenzo (quant-ph/0108010): non-interacting fermion circuits are classically simulable (Valiant's matchgates).
- Aaronson–Arkhipov (1011.3245): linear-optical amplitudes are permanents. Exact classical sampling would collapse the polynomial hierarchy; approximate hardness rests on stated conjectures.

**I22.** "Kasteleyn's face rule puts sign product −1 around faces of length 4k and +1 around faces of length 4k+2." **CORRECT** (bipartite planar form).
- Check (C11), with random positive weights:
  - On a 4 × 5 grid with one interior edge deleted (square faces plus one hexagon), signs solving the face rule over GF(2) give |det K| = 688.584 = the weighted perfect-matching count.
  - Imposing −1 on the hexagon too gives 642.080.
  - On a honeycomb (brick wall) patch, all signs + already give |det| = the count.

**I23.** "This coincides with Lieb's flux rule for half-filled free fermions." **IMPRECISE.**
- The face rule is the *canonical flux* of Lieb–Loss (cond-mat/9209031): flux (ℓ − 1)π through a face of length 2ℓ.
  - Lieb–Loss Thm 3.1: the canonical flux maximises |det T| = ∏σᵢ(K) on every planar bipartite graph, for all hopping magnitudes. This is another proof of Kasteleyn.
- The half-filled ground-state energy is −Σσᵢ(K), a different functional of the same singular values.
  - Lieb (cond-mat/9410025) proves that π per square plaquette (and 0 per hexagon) minimises it, under periodicity in one direction and reflection symmetry.
  - For arbitrary planar bipartite graphs and magnitudes, the energy statement is false (Lieb–Loss §VII(B)).
- Check (C12):
  - 4 × 4 grid, unit hoppings: E(π) = −13.2528 < min over 3000 random fluxes = −13.0766 < E(0) = −10.9443.
  - Lieb–Loss's four-square example with spokes weakened to 0.2: E(π,π,π,π) = −5.2284 > E(π,π,π,0) = −5.3168.
- **Fix:** "The face rule is the |det|-maximising flux on every planar bipartite graph. It also minimises the half-filled energy on the lattices covered by Lieb's theorem (square: π; hexagonal: 0), but not in general."

**I24.** "The sign that tames the permanent is a ground-state gauge field." **IMPRECISE.**
- **Fix:** it is the gauge field that maximises |det| (the product of single-particle energies). It is the ground-state field only where I23's lattice conditions hold.
- For the network, Part 2 unlock 24 keeps the flux picture as SPECULATION. Weight signs there are physical, not a gauge.

### 1.4 Signs, pseudorandomness and operator algebras [M, B, H]

**I25–I26.** "Godsil–Gutman: average det(x − A_s) over random signings s and you get the matching polynomial. Heilmann–Lieb bound its roots by 2√(d−1)." **CORRECT.**
- MSS I (1304.4132) Thm 3.6 (Godsil–Gutman); Heilmann–Lieb, as quoted in MSS I §3.
- Check (C10), Petersen graph:
  - exact average over all 2^15 signings equals μ(G, x), with m_k = 1, 15, 75, 145, 90, 6, to 2e-14;
  - roots are real, with largest root 2.6314 ≤ 2√2.

**I27–I28.** "A signing defines a double cover — a random Z/2 gauge field. Marcus–Spielman–Srivastava turned this into bipartite Ramanujan graphs of every degree." **CORRECT.**
- The new eigenvalues of a 2-lift are those of A_s (Bilu–Linial, as used in MSS I §5).
- MSS I §5 (Thm 5.5 as read): every graph has a signing with λ_max(A_s) ≤ the largest root of μ_G ≤ 2√(d−1). Iterated 2-lifts of K_{d,d} give bipartite Ramanujan graphs of every degree d ≥ 3.
- Check (C10): on the Petersen graph, min_s λ_max(A_s) = 2.000 ≤ 2.631.

**I29–I30.** "The same interlacing-polynomial method then solved the Kadison–Singer problem. Pure states on the diagonal subalgebra of B(ℓ²) extend uniquely." **CORRECT.**
- MSS II (1306.3969): each pure state on the diagonal von Neumann algebra D ⊂ B(ℓ²) has a unique extension.
- The method uses interlacing families plus a multivariate barrier argument. Thm 1.4 and Cor 1.5 give Weaver's form; Thm 6.1 gives Anderson paving.

**I31.** "That diagonal is exactly the 'classical shadow' from our groupoid discussion." **UNVERIFIED.** That conversation is not available. Part 2 (unlock 6, Appendix A) treats it as an analogy: D₀ ⊂ C*(Λ) plays the role of the diagonal masa.

**I32.** "Gurvits's real-stable-polynomial proof of van der Waerden's permanent bound is the ancestor of this method." **CORRECT.** MSS II credit Gurvits's van der Waerden proof as an inspiration for their barrier method, with real stability in the sense of Borcea–Brändén.

**I33.** "So the strongest bridge between pseudo-randomness, permanents and NCG is already a theorem." **IMPRECISE.**
- **Fix:** the theorems are real, but Kadison–Singer is a statement about a commutative masa in B(ℓ²): operator algebras, not noncommutative geometry proper.
- The bridge that is a theorem runs from interlacing families to Ramanujan lifts and to paving.

### 1.5 Zero temperature [T, M, F]

**I34.** "Tropicalization is the zero-temperature limit of a partition function: log-sum-exp becomes max." **CORRECT** for sums of positive terms.
- Proof: max w ≤ T log Σᵢ e^{wᵢ/T} ≤ max w + T log N.

**I35.** "Tropically, permanent and determinant coincide; both become the assignment problem, solvable in polynomial time." **IMPRECISE.**
- As formal tropicalisations, they coincide: the tropical determinant is the tropical permanent, max_σ Σ A_{iσ(i)}, since tropical arithmetic has no subtraction. The assignment problem is a linear program over the Birkhoff polytope, hence polynomial.
- The zero-temperature *limit* of T log|det e^{A/T}| equals the assignment value only when the maximising permutations do not cancel, for example when the maximiser is unique.
- Check (C15):
  - For a 6 × 6 Gaussian A, T log per and T log|det| both tend to the assignment value 4.153169 (agreeing to 1e-6 at T = 0.01).
  - For A = [[1,2],[2,3]], where an even and an odd permutation tie at 4, T log per → 4 but det e^{A/T} = 0 at every T.
  - For [[1,2],[2,3.001]], T log|det| reaches the assignment value only once T ≲ 10⁻³.
- **Fix:** "Tropically per and det have the same tropicalisation. Their zero-temperature limits agree only if the top-weight permutations do not cancel."

**I36.** "Zero temperature, where signs die." **IMPRECISE**, for the same reason: signs die in the formal tropicalisation but not in the limit of a signed sum (C15). This is Part 2 unlock 29 for the network: there, signs are the whole answer.

**I37.** "Chin's proof works exactly there, with c-vectors and tropical T-systems, then lifts to the full system via positivity and the periodicity theorem." **CORRECT.**
- Chin (2602.15140 v2), Thm 1.1: T_i(t + h_Γ + h_Δ) = T_{σ(i)}(t), with σ of order ≤ 2.
- Proof structure:
  - a maximal green sequence (Prop 3.1);
  - −C is a permutation matrix (Cor 3.5);
  - the separation formula;
  - Inoue–Iyama–Keller–Kuniba–Nakanishi's periodicity theorem (Thm 2.14) together with GHKK sign-coherence;
  - folding and transposition through tropical T-systems, using positivity (Lemma 4.3).
- Check (C17), exact rationals on A_m ⊗ A_n for (m,n) = (1,2), (2,2), (1,3), (2,3), (3,3), (2,4), (3,4):
  - T(i, j, u + h + h') = T(m+1−i, n+1−j, u), with σ the 180° rotation, in every case;
  - σ = identity fails wherever it is testable;
  - the full period 2(h + h') holds.

**I38.** "For general permanents, the lift from zero to finite temperature is #P-hard." **CORRECT** for exact evaluation.
- The 0/1 permanent is #P-complete (Valiant; cited in Vontobel §I-A, RST and Gamarnik–Katz).
- For nonnegative matrices, randomised *approximation* is polynomial (the JSV FPRAS, cited in Vontobel §I-B). The hardness is about exact values or signed matrices.

**I39.** "Planarity (Kasteleyn) and expansion-based sampling (JSV) are the classic ways across." **IMPRECISE.**
- **Fix:** JSV's Markov chain mixes rapidly for *every* nonnegative matrix. Its analysis bounds the conductance of the chain's state space; the bipartite graph need not expand.
- Other crossings:
  - deterministic e^{O(n)} approximation via Bethe/Sinkhorn (Gurvits; Linial–Samorodnitsky–Wigderson; Anari–Rezaei);
  - deterministic (1+ε)ⁿ approximation on expanders (Gamarnik–Katz Thm 2).
- Part 2 unlock 32 lists the structures that carry each crossing.

**I40.** "Cluster algebras are a setting where a tropical-to-full transfer is a theorem." **CORRECT.** The IIKKN periodicity theorem, as used by Chin; and Kashaev–Nakanishi Props 2.4, 3.4: tropical periodicity ⟺ classical periodicity ⟺ quantum periodicity of seeds.

### 1.6 Cluster algebras made of matchings [M, T, K]

**I41–I42.** "Cluster variables … are generating functions of perfect matchings of planar snake graphs (Musiker–Schiffler–Williams). These are weighted permanents, made into determinants by Kasteleyn." **CORRECT.**
- MSW (0906.0748) Thm 4.9: for an ordinary arc γ, x_γ = (1/cross(T,γ)) Σ_P x(P) y(P), summed over the perfect matchings P of the snake graph. The crossing monomial is a Laurent monomial.
- Notched arcs use γ-symmetric matchings (Thms 4.16, 4.20).
- Snake graphs are planar and bipartite, so I14 applies.

**I43.** "Octahedron-recurrence solutions are sums over perfect matchings (Speyer)." **CORRECT.** Speyer (math/0402452): the Laurent expansions have coefficients 1 and their terms are perfect matchings of "crosses-and-wrenches" graphs. Dodgson condensation is the λ = −1 Robbins–Rumsey case.

**I44.** "Kuo condensation, an identity among matching counts with vertices deleted, is an exchange relation." **CORRECT.**
- Kuo (math/0304090) Thm 2.1, for a, b, c, d on one face in cyclic order (a, c of one colour and b, d of the other): M(G) M(G−abcd) = M(G−ab) M(G−cd) + M(G−ad) M(G−bc).
- It has the form of the three-term Plücker relation, the exchange relation of the Grassmannian cluster algebra. That reading is from memory: boundary measurements of planar bipartite graphs are Plücker coordinates of Gr_{≥0}.
- Check (C13): 6153.460582 on both sides, for a weighted 4 × 4 grid.

**I45–I47.** "A mutation on a planar bipartite graph is urban renewal … preserves the dimer partition function up to a factor … domino shuffling … samples random tilings exactly … integrable structure replaces the expansion that Markov-chain sampling needs." **CORRECT.**
- Propp (math/0111034): urban renewal with factor (ac + bd) for a square face with weights a, b, c, d in cyclic order. Generalised domino shuffling samples exactly.
- Goncharov–Kenyon (1107.5588): the spider move is a cluster Poisson mutation.
- Check (C14): M(G) = 312.15298547 = (ac+bd) M(G′), where the new square carries weights c, d, a, b divided by (ac+bd), attached by edges of weight 1.
- The last sentence is interpretive. Domino shuffling is an exact sampler and does not rely on mixing.

**I48.** "Dimer models on torus graphs are cluster integrable systems (Goncharov–Kenyon)." **CORRECT** (Goncharov–Kenyon 1107.5588, main theorem; the octahedron recurrence appears there as an example).

**I49.** "Planar Ising models embed in the positive orthogonal Grassmannian (Galashin–Pylyavskyy — the same pair whose periodicity classification Chin builds on)." **CORRECT.**
- Galashin–Pylyavskyy (1807.03282) Thm 2.3: the closure of the space of boundary correlation matrices of planar Ising networks is homeomorphic to OG_{≥0}(n, 2n), a ball. Its cells are indexed by matchings, and Kramers–Wannier duality is the cyclic shift.
- Chin extends their classification of Zamolodchikov-periodic quivers.

**I50.** "Total positivity, where cluster algebras began, governs Markov structure: MTP₂ distributions are faithful to their concentration graphs (Fallat et al.)." **IMPRECISE.**
- Fallat–Lauritzen–Sadeghi–Uhler–Wermuth–Zwiernik (1510.01290):
  - Thm 6.1 needs the graphoid hypothesis (for example, a positive density), and Ex. 5.4 (p₀₀₀ = p₁₁₁ = ½) is MTP₂ without faithfulness.
  - For Gaussians, MTP₂ ⟺ the precision matrix is an M-matrix (Karlin–Rinott). Every MTP₂ Gaussian is faithful to its concentration graph.
- Check (C5): an M-matrix precision on a 4-cycle with a pendant path gives
  - 240 statements i ⊥ j | K, with no negative partial correlation and no mismatch between conditional independence and graph separation;
  - a signed 4-cycle (non-MTP₂) has Σ₀₂ = 0 although no set separates 0 from 2.
- "Total positivity, where cluster algebras began": cluster algebras began with total positivity of *matrices* (Lusztig; Fomin–Zelevinsky; Chin's introduction).
  - MTP₂ is a lattice-supermodularity of densities. For a Gaussian it constrains the inverse covariance, not the minors of the covariance.
  - The two positivities are different (§2.3).
- **Fix:** "MTP₂ graphoids (in particular all MTP₂ Gaussians) are faithful to their concentration graphs; MTP₂ is related to, but is not, the matrix total positivity of cluster theory."

---

## 2. Idea A: Dodgson, CMI and the exchange-relation picture [K, M, T]

### 2.1 The input's scalar identity (I51–I56)

**I51–I53.** Take a and c to be single indices. Then:

det M · det M[B] = det M[aB] · det M[Bc] − (det M[aB|Bc])²,

I(a : c | B) = −½ log(1 − r²), with r² = (det M[aB|Bc])²/(det M[aB] det M[Bc]),

and Markovianity is the degeneration of the relation to one monomial. **CORRECT** for scalars.
- These are Desnanot–Jacobi and the Gaussian entropy formula. In addition, (−1)^{|B|} r is the partial correlation ρ_{ac|B}.
- Check (C1): residual ≤ 3e-18, and the three CMI expressions agree to 12 digits.
- The same holds in Part 2's subtraction-free form p_{aB} p_{Bc} = p_{aBc} p_B + a²_{ac|B}, i.e. I = ½ log(1 + y) with y = a²_{ac|B}/(p_{aBc} p_B) (Part 2 unlock 33).
- For blocks (|a|, |c| > 1) the scalar formula does not generalise as written; §2.2 gives the correct form.

**I54.** "Koteljanskii's determinant inequality is the statement CMI ≥ 0." **CORRECT** for PD matrices.
- det M_S det M_T ≥ det M_{S∪T} det M_{S∩T} is I(S∖T : T∖S | S∩T) ≥ 0 for N(0, M). This is log-det submodularity (Madiman–Tetali 0901.0044; Cover–Thomas, SIAM J. Matrix Anal. Appl. 9 (1988) 384–392).
- Check (C7), 5000 random pairs (S, T): no violation. The minimum is 0 exactly when S ⊆ T or T ⊆ S, and 4.6e-7 otherwise.

**I55–I56.** "A Gaussian is Markov on a chordal graph iff its determinant factorizes as ∏ cliques / ∏ separators. That is the junction-tree formula the Bethe approximation imposes everywhere." **CORRECT.**
- For a junction tree (C₁, …, C_m; separators S_t) the following identity holds for every PD M (DERIVED, chain rule):

  ½(Σ log det M_C − Σ log det M_S − log det M) = KL(N(0,M) ‖ N(0,M̂)) = Σ_t I(C_t∖S_t : H_{t−1}∖S_t | S_t), with M̂⁻¹ = Σ_C [M_C⁻¹]⁰ − Σ_S [M_S⁻¹]⁰.

  It vanishes iff every cut CMI vanishes, i.e. iff the Gaussian is Markov on the graph.
- Check (C6), cliques {0,1,2}, {1,2,3}, {2,3,4}, {3,5}: the three expressions agree (1.308621985157). They vanish (≤ 7e-16) at the Markov projection, which keeps the clique marginals.
- Bethe imposes this formula for the tree of factors and variables (edge terms minus (deg − 1) vertex terms) on graphs that are not trees.

### 2.2 The block identity (vector-valued a and c)

**Setting.** M is PD on a ⊔ B ⊔ c, with |a| = p, |B| = s, |c| = q.
- The Schur complement S = M[ac] − M[ac, B] M[B]⁻¹ M[B, ac] is the conditional covariance of (x_a, x_c) given x_B, with blocks S_aa, S_ac, S_cc.
- The **conditional canonical correlations** ρ₁ ≥ … ≥ ρ_m ≥ 0, with m = min(p, q), are the singular values of R = S_aa^{−1/2} S_ac S_cc^{−1/2}. They lie in [0, 1).

**Theorem A** (DERIVED; checked C2–C4).
1. det M[aBc] · det M[B] = det M[aB] · det M[Bc] · ∏_{i=1}^{m} (1 − ρᵢ²).
2. For x ~ N(0, M):

   I(x_a : x_c | x_B) = ½ log(det M[aB] det M[Bc] / (det M[aBc] det M[B])) = −½ Σᵢ log(1 − ρᵢ²) = ½ log(1 + y),

   with y = 𝓡/(det M[aBc] det M[B]) and 𝓡 = det M[aB] det M[Bc] − det M[aBc] det M[B] ≥ 0.
3. **The remainder as an exchange relation with many terms.** Write positions 1, …, p + q for a ∪ c (a first) and order rows and columns with B first. Then

   𝓡 = − Σ_{J ≠ a, |J| = p} ε(J) · det M[Ba | BJ] · det M[Bc | BJᶜ], with ε(J) = (−1)^{Σ_{i≤p} i + Σ_{j∈J} j}.

   - For p = q = 1 this is the single square (det M[aB|Bc])², and (1)–(2) reduce to I51–I52.
   - For p, q ≥ 2 the terms have mixed signs (C3: sign pattern + + + − −). So 𝓡 ≥ 0 comes from Fischer's inequality, not term by term.
4. **The input's ratio, generalised.** For p = q, r² = (det M[aB|Bc])²/(det M[aB] det M[Bc]) = ∏ᵢ ρᵢ². Hence:
   - −½ log(1 − r²) ≤ I, with equality iff m = 1 or I = 0 (C3: 0.0011 against the true 0.1400);
   - r = 0 iff *some* ρᵢ = 0, which does not imply Markovianity (C4: det M[aB|Bc] = 1e-16 while I = 0.2231 = −½ log(1 − 0.6²)).
5. **Markov criterion.** x_a ⊥ x_c | x_B ⟺ S_ac = 0 ⟺ rank M[aB|Bc] = |B| ⟺ det M[Bi | Bj] = 0 for all i ∈ a, j ∈ c. That is, all p·q scalar Dodgson relations (aᵢ, B, c_j) degenerate. C4: rank 3 = |B| in the Markov case; rank 4 and 5 otherwise.

**Proof.**
- Sylvester / Schur: for I, J ⊆ a ∪ c with |I| = |J|, det M[B∪I | B∪J] = det M[B] · det S[I|J] (B first in rows and columns).
- Hence det M[aBc] = det M[B] det S, det M[aB] = det M[B] det S_aa and det M[Bc] = det M[B] det S_cc.
- det S = det S_aa det S_cc det(I − RRᵀ) gives (1).
- (2) is the Gaussian entropy h = ½ log((2πe)^k det).
- (3) is the generalised Laplace expansion of det S along the rows of a, rewritten with Sylvester.
- (4): det S_ac/√(det S_aa det S_cc) = det R = ±∏ρᵢ; and 1 − ∏ρᵢ² ≥ ∏(1 − ρᵢ²), with equality iff m = 1 or all ρᵢ = 0.
- (5): eliminating with M[B] gives rank M[aB|Bc] = s + rank S_ac, and det M[Bi|Bj] = det M[B] (S_ac)_{ij}. ∎

**Theorem A′: conditional CCA is the optimal Gaussian memory** (DERIVED; checked C18). For k ≤ m,

min over U ∈ ℝ^{p×k} of I(x_a : x_c | x_B, Uᵀx_a) = −½ Σ_{i>k} log(1 − ρᵢ²),

attained by the top-k conditional canonical variates of a.
- Proof:
  - Chain rule: I(a : c | B) = I(Uᵀa : c | B) + I(a : c | B, Uᵀa).
  - The first term is −½ log det(I − Vᵀ RRᵀ V), with V = S_aa^{1/2} U orthonormalised.
  - By Poincaré separation (Cauchy interlacing for compressions), the eigenvalues of VᵀRRᵀV are dominated term by term by ρ₁², …, ρ_k², with equality at the top singular vectors. ∎
- C18 (p = 4, s = 3, q = 3), residual at k = 1: 0.849584 against the formula 0.849584. The best of 2000 random U gives 0.883.

**Copula extension.** CMI is invariant under injective maps of each coordinate, so Theorems A and A′ hold verbatim for a Gaussian copula, with M the latent correlation matrix (Part 2 unlock 46(d)). Guard: this covers pre-activations z, not activations relu(z), since relu is not injective.

**Fermionic caveat** (checked C19).
- For a quasi-free fermionic state with correlation matrix C, S(X) = Σ h(νᵢ(C_X)) with h the binary entropy. This is not ½ log det.
- At points where the Gaussian (log-det) criterion holds exactly (S_ac = 0, Gaussian CMI ≤ 2e-16), the fermionic CMI is 2.4e-4 to 3.6e-3 (four random 4-mode examples).
- So the determinantal exchange-relation picture is a classical Gaussian (and copula) statement. For quasi-free fermions the input's question needs other coordinates. OPEN.

### 2.3 What "Markov network = degenerate exchange relation" means, precisely

Four statements, from proved to open.

**A1 (one conditional independence; THEOREM).**
- The square trinomial a_{ij|K}² − p_{iK} p_{jK} + p_K p_{ijK} = 0 (Desnanot–Jacobi) holds for every symmetric matrix (Boege–D'Alì–Kahle–Sturmfels, "The geometry of gaussoids", 1710.07175).
- For N(0, M), x_i ⊥ x_j | x_K ⟺ a_{ij|K} = 0 ⟺ ρ_{ij|K} = a_{ij|K}/√(p_{iK} p_{jK}) = 0 (BDKS). Equivalently, the exchange relation p_{iK} p_{jK} = p_K p_{ijK} + a²_{ij|K} degenerates to a binomial.
- Check (C5): every trinomial on a random 5 × 5 PD matrix vanishes to 3e-16, and the partial-correlation formula holds to 6e-16.

**A2 (a Markov network; THEOREM).**
- N(0, M) is Markov on G ⟺ the pairwise square trinomials (i, j, V∖ij) degenerate for every non-edge ij. Pairwise ⟺ global holds for positive densities (Hammersley–Clifford), and for Gaussians it is the zero pattern of M⁻¹.
- The global Markov property is then the degeneration of the trinomials (i, j, K) for *every* K that separates i from j.
- So a Gaussian Markov network is exactly a simultaneous degeneration of exchange relations. Theorem A(5) is the block version.

**A3 (the positive part in the statistical sense; THEOREM).**
- MTP₂ Gaussians are exactly the PD matrices whose almost-principal minors are all ≥ 0: MTP₂ is closed under marginalisation, and ρ_{ij|V∖ij} ≥ 0 for all pairs is the M-matrix condition.
- Under MTP₂, conditional-independence strata are graph separations (faithfulness; Fallat et al. Thm 6.1).
- BDKS Thm 5.6: positive gaussoids are exactly the graphical models, and their realisation spaces are open balls of dimension |E| + n.
- Check (C5).

**A4 (a cluster structure in which mutation changes the separator; partly DERIVED, partly OPEN).**

What exists:
- Kenyon–Pemantle ("Principal minors and rhombus tilings", 1404.1354), on the principal and odd almost-principal minors of a matrix:
  - Thm 4.3: the ideal of their relations is generated by translates of one relation, the hexahedron relation, which is "a composition of six cluster mutations".
  - Charts are rhombus tilings of a 2n-gon, and moves between charts are cube moves. In every chart the entries are Laurent polynomials (Thm 4.4).
  - On symmetric matrices the relation becomes the Kashaev relation, and each face satisfies Dodgson: |F(f)|² = F(a)F(c) + F(b)F(d) (Thm 5.2). Checked C21.
  - PD matrices are exactly the "positive networks": vertex signs σ(v) = (−1)^{⌊|S|/2⌋}, face signs free (Thm 5.7).
- Read statistically (DERIVED here):
  - A face of a chart at vertex S for the pair (i, j) carries ±a_{ij|S}, i.e. the partial correlation ρ_{ij|S} up to a positive factor.
  - Each pair occurs exactly once per chart, so a chart is a choice of one conditioning set S_{ij} per pair.
  - A cube move on [S, S ∪ {i,j,k}] (i < j < k) replaces the faces (ij|S), (jk|S), (ik|S∪j) by (ik|S), (jk|S∪i), (ij|S∪k): **one composite mutation (six cluster mutations) changes the separators of three partial correlations.** This is the input's "mutation would mean changing a separator", made literal. For n = 3 it exchanges the charts {ρ₁₂, ρ₂₃, ρ_{13|2}} and {ρ₁₃, ρ_{23|1}, ρ_{12|3}}.
  - In the standard chart the faces are the contiguous minors a_{i,j|i+1..j−1}, a D-vine.
- **Proposition (DERIVED; checked C21).** N(0, M) is a Markov chain on the path 1 − 2 − … − n iff every face of the standard chart with |i − j| ≥ 2 vanishes. So path Markov networks are a coordinate stratum of a chart, where the face relations degenerate to binomials.
  - Proof: (⇒) the global Markov property. (⇐) by contraction, x_n ⊥ x_{n−2} | x_{n−1} and x_n ⊥ x_{n−3} | x_{n−2}, x_{n−1} give x_n ⊥ x_{n−3}, x_{n−2} | x_{n−1}. Iterating gives x_n ⊥ x_{1..n−2} | x_{n−1}; then induct on n. ∎

Obstructions (DERIVED or THEOREM):
- **MTP₂ is not the positive part.** In the standard chart for n = 3, (ρ₁₂, ρ₂₃, ρ_{13|2}) = (0.9, 0.1, 0.9) is PD with all chart faces positive, yet ρ_{23|1} = −0.869 (C21).
  - So MTP₂ is not an orthant of the chart. Any orthant containing an interior MTP₂ point also contains non-MTP₂ points.
  - Face signs of PD matrices are not preserved by cube moves either: the PD point (0.5, 0.5, −0.5) has all three twisted chart-1 faces positive, but its chart-2 face ρ₁₃ = −0.125 (C21).
  - In KP's identification, PD itself has alternating vertex signs, so the all-positive part of the cluster structure contains no PD matrix.
- **Chepuri–George–Speyer** ("Electrical networks and Lagrangian Grassmannians", 2106.15418; AIHP D 13 (2026) 191–216), §5.2 Thm 5.1. For m ≥ 3, no embedding of symmetric m × m matrices into Gr(m, 2m) (by signed column placement) makes the principal and almost-principal minors appear with their signs as Plücker coordinates. "Positive gaussoid" positivity is therefore not totally-nonnegative-Grassmannian positivity, and cannot be repaired by rearranging columns or signs.
  - Contrast: positroid varieties are cluster varieties whose totally positive part is where all cluster variables are positive (Galashin–Lam, Thm 3.5 and Cor 4.4), and Karpman's TNN Lagrangian Grassmannian (1510.04386) is a different positivity again.
- **In plain minor coordinates the separator change is not subtraction-free.** The edge trinomial p_K a_{ij|kK} = p_{kK} a_{ij|K} − a_{ik|K} a_{jk|K} (C5) gives a_{ij|kK} as a rational function that vanishes at points with all arguments positive. Example: the MTP₂ path i − k − j with K = ∅, where a_{ij} > 0 but a_{ij|k} = 0. A subtraction-free expression is positive on the positive orthant, so no such expression exists.
  - KP's sign twists make the hexahedron relation subtraction-free, at the price of the sign mismatch above.
- Quasi-free fermions do not fit (C19).

**Proposition A4′: which Markov networks are strata** (DERIVED, sketch; computer-checked for n ≤ 6).

*Statement.* Let G be a graph on [n]. The Gaussian graphical model of G is the PD part of a coordinate stratum of a KP chart T iff, for every non-edge ij of G, the chart's conditioning set S_ij(T) separates i from j in G. Relabelling the vertices is allowed, since charts depend on the order of the labels.

*Proof sketch.*
- (⇒) Generic MTP₂ points of the model are faithful (Fallat et al. Thm 6.1; C5). So a face vanishing on the whole model has a separating conditioning set. Counting dimensions forces exactly one vanishing face per non-edge.
- (⇐) Separation gives model ⊆ stratum.
  - In a chart, the vertex values along a path, together with all face values, are free coordinates (KP §5.2). The inverse map is polynomial in the face values and Laurent in principal minors: KP Prop 5.5 for the standard chart, then the Kashaev form (7)–(10) of cube moves, which divides only by vertex values.
  - Hence the stratum is irreducible of dimension n(n+1)/2 − #non-edges = n + |E|, the dimension of the model. Both closures therefore coincide.
  - PD points of that closure satisfy (Σ⁻¹)_ij = 0 on non-edges. ∎

*Exhaustive search* over all 2, 8, 62 and 908 charts (n = 3, …, 6) and all connected graphs up to isomorphism:
- every connected graph on ≤ 5 vertices is realised;
- on 6 vertices, all 58 chordal graphs and 53 of the 54 non-chordal graphs are realised. The exception is the triangular prism C₃ □ K₂: all 60 of its labellings fail.

*Local check (C21).* At random model points of a 6-path, a 6-star, a 4-cycle and a chordal 6-vertex graph:
- the non-edge faces vanish (≤ 6e-20);
- their Jacobian has full rank, so the stratum has dimension n + |E| and coincides locally with the model.

**Answer to the input's open question (Gaussian case).**
- *Stratification half* ("Markov networks are boundary strata where exchange relations collapse"): holds in KP's cluster-type structure for most graphs, including every chordal graph on ≤ 6 vertices. On each such stratum the Dodgson face relations of the vanishing faces degenerate to binomials. It fails for at least one graph, the prism.
- *Mutation half* ("mutation = change of separator"): holds literally (cube moves).
- *Positivity half* ("MTP₂ would be the positive part"): false, both for KP's structure and for Grassmannian signs.

**CONJECTURE A4\*.** Every decomposable (chordal) Gaussian graphical model is a coordinate stratum of a KP chart under a suitable labelling. Evidence: all chordal graphs on ≤ 6 vertices. Test: the same search at n = 7 (rhombus tilings of the 14-gon).

**OPEN.** Characterise the realised graphs (the prism is the smallest exception found), and find the analogue for quasi-free fermions (C19).

### 2.4 Status of Idea A

| statement | status |
|---|---|
| scalar Dodgson; CMI = −½ log(1 − r²) = ½ log(1 + y) | THEOREM (classical), checked C1 |
| block identity, CMI = −½ Σ log(1 − ρᵢ²), many-term remainder, rank criterion (Theorem A) | DERIVED, checked C2–C4 |
| conditional CCA optimal rank-k memory (Theorem A′) | DERIVED, checked C18 |
| Koteljanskii = CMI ≥ 0; chordal gap = KL = Σ CMI | THEOREM / DERIVED, checked C6–C7 |
| Markov network = simultaneous degeneration of square trinomials (A1–A2) | THEOREM |
| MTP₂ = all almost-principal minors ≥ 0; strata = graphical models (A3) | THEOREM (BDKS Thm 5.6; Fallat et al. Thm 6.1) |
| cube move = six mutations = change of three separators; path Markov chains = a coordinate stratum (A4) | DERIVED from Kenyon–Pemantle, checked C21 |
| graphical models as coordinate strata of KP charts: criterion (Prop A4′); all chordal graphs on ≤ 6 vertices realised, the prism not | DERIVED (sketch), exhaustive search n ≤ 6 (C21) |
| every chordal model is a stratum (A4*) | CONJECTURE |
| MTP₂ is the positive part | false for KP's charts (C21) and for Grassmannian signs (CGS Thm 5.1) |
| characterisation of realised graphs; quasi-free fermions | OPEN |

### 2.5 What Idea A gives an estimator [K, B, M]

- **Exact KL accounting.** Any design whose state is, or contains, a Gaussian or Gaussian-copula field knows exactly what it loses across a cut, at O((p + q + s)³) per separator:
  - the price of forgetting beyond a separator is −½ Σ log(1 − ρᵢ²) (Theorem A);
  - the information-optimal rank-k memory is the top-k conditional canonical variates (Theorem A′).
- **Guard.** It measures the covariance (Gaussian or copula) channel only.
  - The depth chain is exactly Markov at full resolution (Part 2 unlock 44). Memory exists only relative to a compressed state, which is where Theorem A′ applies.
  - The old content of BRIEF §3 lives in third cumulants, outside this channel unless the state is a copula of a Gaussian field.
- **First experiment, hours.** At width 256 and on the bench set w1024_d16:
  - estimate from samples the covariance of (z_{l−k}, the retained modes of z_l, z_{l+1});
  - compute the conditional canonical spectrum and k(ε) = min{k : residual ≤ ε};
  - compare the optimal modes with the top singular directions of the gated propagators, which BRIEF §3 reports as the carriers of old content.
- **Expected outcome.** A flat spectrum (k(ε) ≈ 0.3 n) would confirm, in information units, that no low-rank Gaussian memory exists, and would charge the negative result to the dictionary: the separator is a whole layer.

---

## 3. Idea B: two local-to-global principles [B, K, H]

The input's picture:
- *separator-based* gluing screens a region by its boundary and is strong on amenable geometry;
- *tree-based* gluing (Bethe, cavity, belief propagation) is strong on pseudorandom geometry while correlations decay, and breaks down where expansion creates glassiness;
- expanders are the worst case for the first and the best case for the second.

Below, each half is stated as theorems with hypotheses, followed by corrections and the quantum half.

### 3.1 The separator principle

- **S1, Hammersley–Clifford (THEOREM).** Let p be a *strictly positive* law on a finite product space. Then p factorises over the cliques of G ⟺ it satisfies the global Markov property ⟺ the local property ⟺ the pairwise property.
  - Positivity is needed: Moussouris's counterexample (via Gandolfi–Lenarda Ex. 3.3; digest [hammersley-clifford.md](../digests/bridges/hammersley-clifford.md)).
  - Brown–Poulin (1206.0755) Thm 1 states it in this form.
- **S2, exact gluing cost (DERIVED; checked C6).** For any law and any junction tree, KL(p ‖ ∏p_C/∏p_S) = Σ_t I(C_t∖S_t : H_{t−1}∖S_t | S_t).
  - The cost of separator gluing is exactly the sum of cut CMIs.
  - For Gaussians it is the log-det gap (§2.1), costing O(s³) per separator of size s. For general laws it costs exp(s) (treewidth).
- **S3, quantum, commuting (THEOREM).** Brown–Poulin:
  - Thm 3: Gibbs states of commuting Hamiltonians whose terms sit on cliques are quantum Markov networks (I(A:C|B) = 0 when B separates A from C).
  - Thm 4: the converse holds on triangle-free graphs; a triangular-lattice example shows it fails in general.
  - Trees: Thm 5 (Poulin–Hastings).
- **S4, quantum, all temperatures, local (THEOREM).**
  - Yang (2609.38007) Thm II.1: finite-range Hamiltonians, bounded degree, *every* β > 0:

    I(A:C|B) ≤ C_β exp(C_β g_{A|A^c} − c_β r),

    with r = dist(A, C) and a boundary gain g ≤ C|∂_e A| on ℤ^D.
  - Chen–Rouzé (2504.02208) Cor III.2, every β with dist(A, C) ≥ 4e²βd:

    I(A:C|B) ≲ r′ |A||C| exp(μ′ min(|A|,|C|) − λ′ dist(A,C)).

  - 1D, all temperatures: Kato–Brandão, C_β exp(−c_β √r) (as tabulated in Yang, Table I).
  - Global Markovianity (|A|, |C| ∝ volume) at low temperature is open (both papers); Chen–Rouzé Cor B.2 obtains it under a uniform local gap.
- **S5, where separators fail (DERIVED).**
  - On an expander, |∂_e A| ≥ h|A|, so S4 gives I ≤ C exp(C′ h |A| − c r). This is vacuous unless r ≳ |A|.
  - Diameters are O(log n), so separator screening says nothing for |A| ≳ log n. Classically, treewidth is Θ(n) on expanders, so exact separator gluing costs exp(Θ(n)).
  - "Yang's bounds say little on expanders" is CORRECT (I65).

### 3.2 The tree principle

- **T1, trees (THEOREM, from memory: Pearl; Yedidia–Freeman–Weiss).** On a tree factor graph, belief propagation (BP) is exact and the Bethe free energy equals −log Z.
- **T2, covers (THEOREM).**
  - Permanents: per_B = lim sup_M (E per over M-covers)^{1/M} (Vontobel Thm 39).
  - Random M-lifts converge locally to the universal cover, a tree, as M → ∞.
  - Every factor graph: the analogous graph-cover characterisation of the Bethe partition function (Vontobel, "Counting in graph covers", 1012.0065; cited in the permanent paper, not opened).
- **T3, free energies on locally tree-like graphs, with no uniqueness needed (THEOREM).**
  - Dembo–Montanari (0804.4726) Thm 2.4: the ferromagnetic Ising model on uniformly sparse graph sequences converging locally to Galton–Watson trees. (1/n) log Z → the Bethe prediction for *all* β ≥ 0 and B ∈ ℝ. BP converges exponentially fast for B > 0 (Thm 2.6), and local marginals are approximated (Thm 2.7).
  - Sly–Sun (1203.2602) Thm 4: for every two-spin model (antiferromagnetic Ising, hard-core) on random d-regular (bipartite) graphs, the free energy equals the Bethe prediction, also beyond uniqueness.
  - Monomer–dimer at every activity: estimable on every BS-convergent sequence (ACFK Thm 1.2), with the tree value at large girth (C9).
- **T4, marginals and algorithms need uniqueness (THEOREM).**
  - Weitz: below the uniqueness threshold λ_c(Δ), strong spatial mixing makes marginals bounded-radius tree computations and gives an FPTAS (from memory; quoted by Sly–Sun).
  - Sly–Sun Thms 1–2: above λ_c(d) = (d−1)^{d−1}/(d−2)^d, approximating the hard-core partition function on d-regular graphs is NP-hard.
  - Sly–Sun Thm 5: on bipartite random regular graphs the law splits into ± phases, so a vertex marginal is not determined by its neighbourhood.
- **T5, zero temperature (THEOREM).** Perfect matchings are not estimable (ACFK Thm 1.8). They are estimable on d-regular bipartite expanders (ACFK §1.3 with Gamarnik–Katz) and take the tree (Schrijver = Bethe) value at large girth (ACFK Thm 1.5).
- **T6, glassiness (UNVERIFIED here, from memory).** In random constraint satisfaction problems past condensation, the replica-symmetric (Bethe) free energy is wrong and one-step replica symmetry breaking (1RSB) is needed (Krzakala et al.; Ding–Sly–Sun for k-SAT). The cause is frustration on locally tree-like graphs, not expansion as such: ferromagnets on the same expanders obey T3.

### 3.3 Verdicts on Idea B (I58–I65)

**I58.** "Separator-based (Hammersley–Clifford, Yang's CMI bounds) … strong on amenable geometry." **CORRECT** with S1's positivity hypothesis. Yang's bound decays in r only after paying exp(C|∂A|), which is small relative to volume only on amenable geometry.

**I59.** "Tree-based … strong on pseudorandom geometry while correlations decay along the tree. It breaks down where expansion creates glassiness." **IMPRECISE.**
- **Fix:**
  - Free energies are Bethe-exact on locally tree-like graphs for unfrustrated models (ferromagnets, bipartite two-spin systems), *beyond* uniqueness (T3).
  - Decay (uniqueness) is what marginals and worst-case algorithms need (T4).
  - The breakdown comes from frustration (T6), not from expansion itself.
- For the competition this matters: the score is a vector of marginals (per-neuron means), so a tree/cavity design needs a decay hypothesis along its computation trees, not only a correct free energy.

**I60.** "Expanders are the worst case for the first principle and the best case for the second, and the permanent shows both behaviours." **IMPRECISE.**
- Worst for separators: CORRECT (S5).
- Best for trees only for unfrustrated models, or for zero-temperature counting with expansion (T3, T5). For frustrated models, expanders are bad for both principles.
- The permanent shows both: CORRECT. Perfect matchings are rigid without expansion (ACFK Thm 1.7) and local with it (ACFK Thm 1.9 with Gamarnik–Katz).

**I61.** "Substitution towers converge to the plane, which is amenable. Towers of graph covers converge to a tree, which is not." **IMPRECISE.**
- **Fix:** a tower of finite covers G_k → G₀ converges locally to the universal cover (a tree) iff the girth (injectivity radius) tends to ∞. For normal towers, that means the corresponding subgroups of π₁(G₀) intersect trivially; otherwise the limit is an intermediate cover.
- Substitution (inflation) towers of a tiling converge locally to the tiling space, an amenable unimodular limit.

**I62.** "Bethe permanents are limits over covers, and Ramanujan graphs come from 2-lifts." **CORRECT**, with the cover limit a lim sup (Vontobel Thm 39) and the lifts as in MSS I §5.

**I63–I64.** "Towers of finite quotients ('box spaces') are expanders exactly when the profinite action has a spectral gap. These are the same box spaces whose ghost projections break the coarse Baum–Connes conjecture." **CORRECT**, with the hypotheses below.

*Box-space theorem* (proof here; standard).
- *Setting.* Γ is generated by a finite symmetric set S. N₁ ⊇ N₂ ⊇ … are finite-index normal subgroups. K = lim← Γ/N_k carries Haar probability μ, and M = |S|⁻¹ Σ_{s∈S} π(s) is the Markov operator.
- *Statement.* The Cayley graphs Cay(Γ/N_k, S) form an expander family iff Γ acting on L²₀(K, μ) has a spectral gap. Equivalently, Γ has property (τ) with respect to (N_k).
- *Proof.*
  - L²(K, μ) is the closure of ∪_k V_k, where V_k is the space of functions factoring through Γ/N_k, and V_k ≅ ℓ²(Γ/N_k) equivariantly.
  - On V_k ⊖ ℂ, M is the random-walk operator of Cay(Γ/N_k, S).
  - A uniform bound ⟨Mf, f⟩ ≤ (1−ε)‖f‖² on every V_k ⊖ ℂ extends to L²₀ by density. Conversely, each V_k ⊖ ℂ is an invariant subspace of L²₀.
  - Expansion ⟺ a uniform one-sided gap (Cheeger / Alon–Milman). ∎
- *Coarse Baum–Connes, the hypotheses actually proved.*
  - Higson–Lafforgue–Skandalis (GAFA 12 (2002) 330–354; introduction opened): expanders give counterexamples.
  - Willett–Yu (1012.4151): for expanders of large girth the *maximal* coarse assembly map is an isomorphism, while the ghost (Kazhdan-type) projection lies outside the image of the *reduced* one. So it is the reduced conjecture that fails. Also, geometric property (T) of a box space ⟺ property (T) of Γ.
  - The ℓ^p paper ("Expanders are counterexamples to the ℓ^p coarse Baum–Connes conjecture", EMS Press): for box spaces of hyperbolic groups (for example free groups, SL₂(ℤ)) that are expanders, i.e. with (τ) relative to (N_k) and ∩N_k = {e}, the Kazhdan projection is not in the image of the ℓ^p assembly map.
  - The general statement "every expander box space violates coarse Baum–Connes surjectivity" is from memory.

**I65.** "The quantum frontier is a quantum Bethe theory for Gibbs states on expander interaction graphs, for example via the Markov entropy decomposition on random lifts." Stated precisely as Conjecture QB below.

### 3.4 The quantum half

**What is proved.**
- High temperature: for β below the convergence radius of the cluster expansion (≈ 1/(e d J), from memory), local reduced states and free energy densities are analytic and determined by local structure on any bounded-degree graph. Random lifts converge locally to the d-regular tree T_d, so tree-determined limits follow. (Theorem-level; from memory.)
- Commuting Hamiltonians: S3 (exact Markov structure; trees).
- Local Markov at all temperatures: S4.
- Variational side: Poulin–Hastings (1012.2050), the Markov entropy decomposition (MED), gives F_MED ≤ F, a lower bound on the free energy obtained by bounding entropy with conditional entropies over Markov shields. Its dual is a quantum belief propagation.

**CONJECTURE QB (quantum Bethe exactness on random lifts).**
- *Setting.* H₀ is a fixed finite d-regular base graph. h is a fixed local Hamiltonian density: two-qubit edge terms ‖h_e‖ ≤ J and one-site terms. G_N is a uniform random N-lift of H₀, H_N the lifted Hamiltonian and ρ_N = e^{−βH_N}/Z_N.
- *Quantum tree uniqueness threshold.* β_u(d, h) is the supremum of β such that, on the depth-R ball of T_d, the root's reduced Gibbs state becomes independent of the boundary state on the leaves as R → ∞, uniformly over boundary states.
- *Claims.* For β < β_u:
  - (i) (1/N) log Z_N converges in probability to a deterministic f(β);
  - (ii) one- and two-site reduced states of ρ_N converge to those of the tree with the unique boundary condition;
  - (iii) the MED lower bound with radius-1 shields is asymptotically tight: F_MED/N → f.
- *Expectation for β > β_u.* By analogy with T3 vs T4: (i) may survive for unfrustrated h (quantum ferromagnets), but (ii) fails.
- *Evidence.*
  - QB holds at high temperature, which is the cluster-expansion statement above.
  - The commuting (classical) restrictions reduce to T3–T4 and Brown–Poulin.
  - The CMI bounds S4 are vacuous on these graphs (S5), so they neither support nor contradict QB.
- *Test, about a day.*
  - System: transverse-field Ising, h_e = −J ZZ, h_v = −Γ X, on random 3-regular lifts with N = 12–20 (exact diagonalisation or typicality).
  - Measure: (1/N) log Z, ⟨X_v⟩ and ⟨Z_u Z_v⟩ against β, over many lifts. Compare with an N → ∞ extrapolation and with the MED / quantum-cavity value.
  - QB predicts convergence with corrections exponentially small in the girth below β_u, and lift-to-lift scatter of marginals above it.
  - Caveat: at N ≤ 20 the girth is 3–5, so the test is indicative only.

QB concerns Gibbs states on expanders. The network-facing quantum conjecture is a different one: Part 2 unlock 49, B-programme Conj 8.4, on CMI decay across depth windows for coherent histories. That is where the noncommutative direction would enter an estimator.

### 3.5 Network reading (the realisation dictionary)

Part 2 unlock 47 already applies the two principles to the network: separators along depth (exactly Markov at full resolution, unlock 44), trees across width (dense pseudorandom fresh weights, unlock 13). This note adds three precisions:
- Per-neuron means are marginals, so the tree principle needs a decay hypothesis (T4), not only a correct free energy (T3).
  - The operational test is to perturb an earlier layer and measure how a per-unit mean responds, as a function of depth.
  - BRIEF §3 reports bounded spectral independence (η ≈ 2–8 at width 1024): bounded influence, which is weaker than decay.
- The separator principle along depth is exact only for uncompressed separators; Theorem A′ prices compression exactly in the Gaussian channel (§2.5).
- Expansion across width is what makes separators fail there (S5); "no low-rank old content" (BRIEF §3) is its expected signature.

---

## 4. Idea C: the quantum bipartite belt [H, T]

### 4.1 The construction made precise (I66–I67)

*Data.* A bipartite quiver with skew-symmetric exchange matrix B and vertex set I = white ⊔ black. There are no arrows inside a colour class, so b_ij = 0 there.

*Quantum torus.* 𝒯_q is generated by Y_i^{±1} with Y_i Y_j = q^{2b_ij} Y_j Y_i. Use its skew field of fractions, or a completion.

*Quantum mutation at k* (Kashaev–Nakanishi 1104.4630 Prop 3.1; Keller 1102.4148): μ_k = Ad(Ψ_q(Y_k^{ε})) ∘ τ_{k,ε}.
- Ψ_q is the quantum dilogarithm.
- τ_{k,ε} is a monomial (tropical) change of variables.
- The result does not depend on the sign ε.

*Layers.* μ_◦ = ∏_{k white} μ_k, and μ_• likewise.
- Within a colour class the Y_k commute (b = 0), so the factors commute and each layer is well defined.
- Ad Ψ_q(Y_k) moves only Y_k and its neighbours, so each gate is local in the Heisenberg picture.
- "A brickwork circuit" is therefore **CORRECT** as a statement about automorphisms (μ_• μ_◦)^t of 𝒯_q.

*Well-definedness at three levels.*
1. **Formal** (|q| < 1, or q formal): μ_k is an automorphism of the completed quantum torus. Checked (C20) at q = 0.6 and q = 0.3 + 0.5i, to total degree 12–14:
   - the pentagon identity — the A₂ case of quantum periodicity — E(V) E(U) = E(U) E(−UV) E(V) for VU = qUV, with E(x) = Σ xⁿ/(q;q)_n;
   - Schützenberger's E(U) E(V) = E(U + V);
   - all residuals ≤ 2e-14, while the sign variant E(U) E(−VU) E(V) fails by O(1).
2. **Analytic** (|q| = 1, q = e^{iπb²}, b real): Fock–Goncharov's positive representations on L²(ℝ^I), with Y_k = e^{x̂_k}.
   - Here μ_k is implemented by the unitary Φ_b(x̂_k) (Faddeev's non-compact dilogarithm) composed with a metaplectic (Gaussian) unitary for τ.
   - KN Thm 4.6: the operator of a period is a scalar λ with |λ| = 1. Thm 4.7: λ = 1 when the period is written with the tropical sign sequence.
   - So the circuit is a genuine unitary circuit on L²(ℝ^I): non-Gaussian single-quadrature gates and Gaussian couplers.
   - It is *not* a finite-dimensional tensor-product circuit. At roots of unity there are finite-dimensional cyclic representations (from memory), with normalisation issues.
3. **C\*-algebraic** (Connes' torus, see I71): μ_k is *not* an automorphism. For irrational θ the gauge action U ↦ zU makes spec(U) = 𝕋, so 1 + q^{2a−1}U is not invertible in A_θ.

### 4.2 Periodicity (I68–I70)

**I68.** "Periodicity transfers to the quantum Y-system through the same tropical data; this is how Kashaev and Nakanishi extract quantum dilogarithm identities." **CORRECT.**
- KN Props 2.4 and 3.4: a sequence of mutations is a ν-period tropically ⟺ classically ⟺ quantum mechanically. ν-periods are allowed; KN's A₂ example is the 5-step period with ν a transposition and tropical signs ε = (+,+,−,−,−).
- KN Thm 3.5: the product of quantum dilogarithms in tropical form equals 1.
- Keller (1102.4148) corroborates.

**I69.** "So this interacting circuit should be exactly periodic up to a phase, and by Chin its half-period acts as a permutation." **CORRECT** for the period; **CORRECT** classically for the half-period; quantum half-period DERIVED (sketch).
- Period: KN Thm 4.6–4.7 (exact, phase 1 with tropical signs). Zamolodchikov periodicity of the classical system is Keller's theorem (from memory; used by Chin and KN).
- Classical half-period: Chin Thm 1.1, checked in exact arithmetic on A_m ⊗ A_n T-systems with σ the 180° rotation (C17).
- Quantum half-period: Chin's half-period is a σ-period at the tropical level, and KN Prop 3.4 transfers tropical ν-periods to quantum ones. So the half-period circuit is σ (a permutation of qubits/modes) times a phase.
- Both ingredients are proved, but neither paper states the combination; that is why it is labelled DERIVED (sketch).

**I70.** "Its classical shadow is a dilogarithm identity whose value is a central charge." **IMPRECISE.**
- The classical limit is the dilogarithm identity of KN Thm 2.7: Σ over the period of Rogers' L equals (π²/6)·N₋, an integer multiple of π²/6.
- Check (C16), A₂: (6/π²) Σ_{k=1}^{5} L(y_k/(1+y_k)) = 3.000000000000 and (6/π²) Σ L(1/(1+y_k)) = 2.000000000000 for random positive initial data.
- Central charges are per-step values at the *constant* solution. At y = φ (the golden ratio), (6/π²) L(1/(1+φ)) = 0.4000000000, the Lee–Yang effective central charge 2/5 (from memory), and 5 × 2/5 = 2.
- **Fix:** "a dilogarithm identity whose value is an integer multiple of π²/6, equal to (period) × (effective central charge of the constant solution)".

**Entanglement and CMI (the input's open question).**
- At the period the circuit is a phase times the identity (at the half period, times σ), so entanglement entropies and CMIs of any normalisable state return exactly.
- In between, the Gaussian couplers generate entanglement and the Φ_b gates make it non-Gaussian.
- The question is well posed in L²(ℝ^I), with normalisable inputs such as Gaussian wavepackets. It is SPECULATION whether the returning dynamics bounds the intermediate CMI.

### 4.3 The noncommutative torus, Effros–Shen and Penrose (I71–I73)

**I71.** "At irrational q the quantum torus is Connes' noncommutative torus." **IMPRECISE.**
- For q = e^{2πiθ} with θ irrational, the Laurent algebra with *unitary* generators is the dense polynomial core of the rotation algebra A_θ.
- But quantum cluster theory uses a different *-structure: self-adjoint, positive Y_k in Fock–Goncharov's representations.
- Mutations do not act on A_θ (§4.1, level 3).
- **Fix:** "The quantum torus has Connes' torus as one C*-completion (unitary *-structure), on which cluster mutations do not act. They act on its fraction field and in the positive representations."

**I72.** "At the golden angle it embeds (Pimsner–Voiculescu) in the Effros–Shen AF algebra with matrix [[1,1],[1,0]]." **CORRECT.**
- Pimsner–Voiculescu, J. Operator Theory 4 (1980) 201–210 (IMAR scan opened). Rieffel (Pacific J. Math. 93 (1981)) describes it: A_θ "can be embedded in one of the special AF algebras constructed by E. G. Effros and C. L. Shen whose K₀ group is Z + Zα, ordered as a subgroup of the real line".
- The Effros–Shen Bratteli matrices are [[a_k, 1],[1, 0]] for θ = [0; a₁, a₂, …]. For θ = 1/φ all a_k = 1. For the golden angle 1/φ² = [0; 2, 1, 1, …] this holds from the second level on, which does not change the stable isomorphism class.

**I73.** "… which is Morita-equivalent to the Penrose algebra." **CORRECT.**
- Connes' *Noncommutative Geometry* (introduction, opened): the Penrose-tiling C*-algebra is AF with K₀ = ℤ² ordered by the golden-ratio half-plane.
- Kellendonk–Lenz–Savinien (1412.5442) note that its paths are those of the Bratteli diagram with inclusion graph A₄.
- Both dimension groups are order-isomorphic to ℤ + ℤφ ⊂ ℝ. By Elliott's classification of AF algebras they are stably isomorphic, hence strongly Morita equivalent (Brown–Green–Rieffel). Elliott and Brown–Green–Rieffel are from memory.

### 4.4 Assessment

- The construction is well defined (levels 1–2) and exactly periodic up to a phase.
- It is a rare explicit example of an interacting (non-Gaussian) unitary dynamics with guaranteed exact return, which makes it of mathematical interest.
- It does not touch the competition:
  - its state space is L²(ℝ^I) (or a formal torus), not a law on faces;
  - nothing in the network supplies a periodic exchange matrix.
- Part 2 Appendix A ("SPECULATION; not used") stands for the competition. For the programme it is a THEOREM-level example of the Heisenberg-picture objects of stream [H].

---
