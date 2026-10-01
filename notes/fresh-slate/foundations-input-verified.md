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
2. **Idea A, block form** (DERIVED: proved in §2, checked C1–C4 and C18). For vector-valued a and c:
   - det M[aBc] · det M[B] = det M[aB] · det M[Bc] · ∏ᵢ (1 − ρᵢ²), where ρᵢ are the canonical correlations of a and c conditional on B. Hence I(a : c | B) = −½ Σᵢ log(1 − ρᵢ²).
   - The input's cross-minor ratio generalises to ∏ᵢ ρᵢ², which can vanish without Markovianity.
   - Markov ⟺ rank M[aB|Bc] = |B| ⟺ all |a|·|c| scalar exchange relations degenerate at once.
   - The top-k conditional canonical variates are the optimal rank-k Gaussian memory, with residual −½ Σ_{i>k} log(1 − ρᵢ²).
3. **Idea A, structure.**
   - For Gaussians, "a Markov network is the simultaneous degeneration of exchange relations" is a theorem.
   - Kenyon–Pemantle's cluster-type structure on principal and almost-principal minors goes further (hexahedron move = six cluster mutations):
     - a cube move replaces the conditioning sets of three partial correlations ("mutation = change of separator", DERIVED);
     - a Gaussian graphical model is a coordinate stratum of a chart iff some chart's conditioning sets separate every non-edge (DERIVED, sketch);
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
- Uniform form (proof here). If uniformity failed, there would be pairs (G_n, H_n) of graphs of maximum degree d whose radius-n statistics differ by less than 1/n while their values of ln M/v differ by at least ε. The space of local statistics is compact, so a subsequence of G_n converges, and H_n converges to the same limit. Interleaving the two gives a convergent sequence along which ln M/v does not converge, a contradiction. So for every ε there is a radius R(ε, d) such that radius-R statistics determine ln M/v to ±ε, uniformly over such graphs.
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
- What a single Pfaffian gets wrong is the sign of each perfect matching, which depends on the mod-2 homology class of its symmetric difference with a reference matching. The 4^g terms with Arf signs cancel exactly these signs.

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
  - a signed 4-cycle (non-MTP₂) has Σ₀₂ = 0, i.e. 0 ⊥ 2 marginally, although the empty set does not separate 0 from 2: it is unfaithful.
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
   - For p, q ≥ 2 the terms have mixed signs in every instance tried (C3). So 𝓡 ≥ 0 comes from Fischer's inequality, not term by term.
4. **The input's ratio, generalised.** For p = q, r² = (det M[aB|Bc])²/(det M[aB] det M[Bc]) = ∏ᵢ ρᵢ². Hence:
   - −½ log(1 − r²) ≤ I, with equality iff m = 1 or I = 0 (C3: 0.0011 against the true 0.1400);
   - r = 0 iff *some* ρᵢ = 0, which does not imply Markovianity (C4: det M[aB|Bc] ≈ 1e-15 while I = 0.2231 = −½ log(1 − 0.6²)).
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
- C18 (p = 4, s = 3, q = 3):
  - k = 1: residual 0.4452538821, equal to the formula; the best of 4000 random U gives 0.4511;
  - k = 2: residual 0.0715099386, equal to the formula; the best random U gives 0.1078.

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
  - Charts are rhombus tilings of a 2n-gon, and moves between charts are cube moves. In every chart the matrix entries are Laurent polynomials in the chart's variables (KP §1; Thm 4.4 for the standard chart).
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
  - Brown–Poulin (1206.0755) Thm 1 states it in this form (positive Markov network ⟺ Gibbs form e^{Σ_Q h_Q} over cliques).
- **S2, exact gluing cost (DERIVED; checked C6).** For any law and any junction tree, KL(p ‖ ∏p_C/∏p_S) = Σ_t I(C_t∖S_t : H_{t−1}∖S_t | S_t).
  - The cost of separator gluing is exactly the sum of cut CMIs.
  - For Gaussians it is the log-det gap (§2.1), costing O(s³) per separator of size s. For general laws it costs exp(s) (treewidth).
- **S3, quantum, commuting (THEOREM).** Brown–Poulin (1206.0755):
  - Thm 2 (Leifer–Poulin): every positive quantum Markov network is the Gibbs state of a clique-local Hamiltonian, which need not be commuting.
  - Thm 3: Gibbs states of commuting clique-local Hamiltonians are quantum Markov networks (I(A:C|B) = 0 whenever B shields A from C).
  - Thm 4: on graphs whose only cliques are edges (triangle-free), positive quantum Markov networks are exactly these Gibbs states. Fig. 2 gives a positive quantum Markov network on a graph with triangles (it tiles the plane) that is not of this form.
  - Thm 5 (Poulin–Hastings): the same equivalence on trees.
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
   - KN Thm 4.6: the operator of a period is a scalar λ with |λ| = 1, times the operator of the permutation ν for a ν-period. Thm 4.7: λ = 1 when the period is written with the tropical sign sequence.
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
- *On the noncommutative geometry of tilings* (arXiv:1412.5442; Julien–Kellendonk–Savinien, authors from memory) notes that Connes described the Penrose tilings as the paths of the Bratteli diagram with inclusion graph A₄.
- Both dimension groups are order-isomorphic to ℤ + ℤφ ⊂ ℝ. By Elliott's classification of AF algebras they are stably isomorphic, hence strongly Morita equivalent (Brown–Green–Rieffel). Elliott and Brown–Green–Rieffel are from memory.

### 4.4 Assessment

- The construction is well defined (levels 1–2) and exactly periodic up to a phase.
- It is a rare explicit example of an interacting (non-Gaussian) unitary dynamics with guaranteed exact return, which makes it of mathematical interest.
- It does not touch the competition:
  - its state space is L²(ℝ^I) (or a formal torus), not a law on faces;
  - nothing in the network supplies a periodic exchange matrix.
- Part 2 Appendix A ("SPECULATION; not used") stands for the competition. For the programme it is a THEOREM-level example of the Heisenberg-picture objects of stream [H].

---

## 5. Recommendation

**Ranking.**
- For the competition: A > B > C.
- For the programme: B > A > C.

**Push A first, for the competition, as an instrument.**
- It is the only one of the three that gives *exact* computations on the competition's objects today:
  - the KL price of a Gaussian or copula separator (Theorem A);
  - the information-optimal rank-k memory (Theorem A′);
  - a rank test for Markovianity.
- Each costs O(s³) per cut, and the Stage Q experiment of §2.5 takes hours.
- It decides one design question quickly: does a low-rank Gaussian memory across depth exist, and if so which modes carry it?
- Its limit is explicit. It sees the covariance channel only. The Gaussian closure already reaches raw ≈ 4.3e-6 (BRIEF §1), so the error that remains is non-Gaussian. A filters designs and prices their separators; it is not an estimator and cannot reach raw 1e-8 by itself.

**Develop B, for the programme.**
- B is the precise form of the programme's local-to-global question. Both classical halves are theorems with explicit hypotheses (§3.1–3.2).
- It corrects a point that matters for every tree or cavity design: marginals need decay, free energies do not (I59).
- It isolates the quantum half in one testable conjecture (QB), next to the network-facing Part 2 Conj 8.4. That is where the NCG direction would enter.
- It also absorbs Idea A's best structural result: Kenyon–Pemantle charts realise "mutation = change of separator" and give Markov networks as coordinate strata (§2.3). This is the algebra of the separator principle.

**Park C.**
- It is well defined and returns exactly, but it has no route into the competition.
- Its noncommutative-geometry link (Connes' torus) is weaker than the input states (§4).

---

## 6. Numerical checks

**Setting.** Run in the session scratchpad: numpy 2.4.6 and Python 3.11 in the whest environment; scipy, mpmath and sympy are not installed. OPENBLAS_NUM_THREADS=1, fixed seeds, all sizes ≤ 16 × 16 (exhaustive searches up to 2^15 signings and 908 charts).

**Two sources of numbers.**
- The table reports the session runs quoted in §§1–4.
- The two compact scripts below recheck every identity on fresh random instances. Their random streams differ, so instance-dependent numbers differ while the identities hold to the same precision. Their outputs follow each script.

| check | claim (where used) | session result |
|---|---|---|
| C1 | scalar Dodgson; CMI three ways (I51–I52) | 3 random 6 × 6 PD: residual ≤ 3e-18; CMI agree to 12 digits; (−1)^{\|B\|} r = partial correlation to 10 digits |
| C2 | block identity, Theorem A(1)–(2) | (p,q,s) = (2,2,3), (2,3,3), (3,2,2), (3,3,4): \|lhs/rhs − 1\| ≤ 7e-15; CMI = −½ Σ log(1−ρᵢ²) to 12 digits |
| C3 | naive ratio = ∏ρᵢ²; many-term remainder, Theorem A(3)–(4) | r² = ∏ρᵢ² to 1e-16; naive −½ log(1−r²) = 0.0011 against CMI 0.1400; 𝓡 = sum of 5 signed products (0.47598955426 both), signs mixed |
| C4 | rank criterion; vanishing cross minor without Markov, Theorem A(4)–(5) | ρ = (0.6, 0): det M[aB\|Bc] = 1e-15, rank 4 = \|B\|+1, CMI 0.223144; ρ = 0: rank 3, CMI 0; ρ = (0.6, 0.3): rank 5 |
| C5 | gaussoid trinomials; partial correlations; MTP₂ faithfulness; unfaithful signed 4-cycle (A1–A3, I50) | residuals ≤ 6e-16; 240 statements i ⊥ j \| K: 0 negative partial correlations, 0 CI/separation mismatches; signed 4-cycle Σ₀₂ = −1e-17 |
| C6 | chordal gap = KL = Σ CMI; zero at the Markov projection (I55) | 1.308621985157 three ways; projection ≤ 7e-16 |
| C7 | Koteljanskii (I54) | 5000 random pairs: no violation; min 0 for nested pairs, 4.6e-7 otherwise |
| C8 | Heawood spectrum and Pfaffian signing; K₃,₃ (I17–I18) | ±3, ±√2 (×6); per = 24 = max \|det\| over the 2⁸ gauge-fixed signings; K₃,₃: per 6, max \|det\| 4 |
| C9 | Bethe = Schrijver on regular graphs; Gurvits and Anari–Rezaei; tightness; Kesten–McKay constants (I2, I8–I11) | per_B = (4/3)^10 = 17.7577266338 and (27/16)^9 = 110.9670481672 to 10 digits, at P = A/d; 40 random: min per/per_B = 1.000000, max per/(2^{n/2} per_B) = 0.848; I ⊗ J₂: per/per_B = √2ⁿ (n = 2, 4, 6); ∫½ln(1+x²)dKM₃ = 0.58157540 = tree value = ½ ln(16/5); ∫ln\|x\|dKM_{3,4} = Schrijver to 3e-6 |
| C10 | Godsil–Gutman; Heilmann–Lieb; MSS signing (I25–I28) | Petersen graph, all 2^15 signings: \|E_s charpoly − μ\| ≤ 1.8e-14; real roots, max 2.6314 ≤ 2√2; min_s λ_max(A_s) = 2.0000 |
| C11 | Kasteleyn face rule with a hexagonal face; honeycomb (I22) | 4 × 5 grid: \|det\| = count (297.889); with a hexagon: face rule 688.584 = count, −1 on every face 642.080; honeycomb, all +: \|det\| = count (20.138) |
| C12 | \|det\| versus energy; Lieb's π flux; Lieb–Loss counterexample (I23–I24) | 4 × 4: E(π) = −13.2528 < min of 3000 random fluxes −13.0766 < E(0) = −10.9443; \|det K\| 36 at π, 0 at 0; four squares, spokes 0.2: E(π⁴) = −5.2284 > E(π,π,π,0) = −5.3168 (spokes 0.05: −4.9284 vs −5.2322); 2 × 3 ladder: no counterexample in 400 random magnitude sets |
| C13 | Kuo condensation (I44) | weighted 4 × 4 grid: 6153.460582 on both sides (another instance: 3546.312743) |
| C14 | urban renewal (I45) | M(G) = 312.15298547 = (ac+bd) M(G′) (another instance: 235.81896693) |
| C15 | tropical limits (I34–I36) | 6 × 6 Gaussian: T log per and T log\|det\| → 4.153169 (4.153170 both at T = 0.01); tie [[1,2],[2,3]]: per → 4, det ≡ 0; [[1,2],[2,3.001]]: T log\|det\| = 3.540, 3.977, 4.0005, 4.0010 at T = 10⁻¹…10⁻⁴ |
| C16 | A₂ dilogarithm identity; constant solution (I70) | (6/π²) Σ L(y/(1+y)) = 3.000000000000 and Σ L(1/(1+y)) = 2.000000000000 for 3 random seeds; L(½) = π²/12; (6/π²) L(1/(1+φ)) = 0.4000000000 |
| C17 | Chin half-periodicity on A_m ⊗ A_n T-systems (I37, I69) | exact rationals, 7 pairs (m, n): σ = rotation holds in all; identity fails wherever testable; full period 2(h+h′) holds |
| C18 | conditional CCA = optimal memory, Theorem A′ | (p,s,q) = (4,3,3): k = 1 residual 0.4452538821 = formula, best of 4000 random 0.4511; k = 2: 0.0715099386 = formula, best random 0.1078 |
| C19 | fermionic caveat (§2.2) | 4 complex 4-mode examples: log-det CMI ≤ 2e-16, fermionic CMI 2.4e-4 to 3.6e-3 |
| C20 | quantum pentagon; Schützenberger (§4.1) | q = 0.6 and 0.3 + 0.5i, degree ≤ 14: residuals ≤ 2e-14; sign variant off by O(1) |
| C21 | Kenyon–Pemantle charts (A4, Prop A4′) | Markov chain: contiguous faces of lag ≥ 2 ≤ 1e-19; (0.9, 0.1, 0.9) → ρ_{23\|1} = −0.8692; (0.5, 0.5, −0.5): chart-1 twisted faces (0.5, 0.5, 0.375), chart-2 face −0.125; KP face relation ≤ 2e-10 (4 × 4, minors ~10³); chart counts 2, 8, 62, 908; strata: all connected graphs on ≤ 5 vertices; on 6, all 58 chordal and 53/54 non-chordal (prism fails); stratum Jacobian of full rank for path₆, star₆, C₄ and a chordal 6-graph |

**Script A** (Idea A: C1–C7, C18, C19, C21).

```python
# Idea A checks C1-C7, C18, C19, C21 (numpy only)
import numpy as np, itertools
rng = np.random.default_rng(1)
D = lambda M, r, c=None: np.linalg.det(M[np.ix_(r, r if c is None else c)]) if len(r) else 1.0
spd = lambda n: (lambda A: A @ A.T/(n+2))(rng.standard_normal((n, n+2)))
cmi = lambda M, a, B, c: 0.5*np.log(D(M, a+B)*D(M, B+c)/(D(M, B)*D(M, a+B+c)))
def isq(S): w, V = np.linalg.eigh(S); return V @ np.diag(w**-0.5) @ V.T
def canon(M, a, B, c):   # conditional canonical correlations of a and c given B
    I = a+c; S = M[np.ix_(I, I)] - M[np.ix_(I, B)] @ np.linalg.solve(M[np.ix_(B, B)], M[np.ix_(B, I)]); p = len(a)
    return np.linalg.svd(isq(S[:p, :p]) @ S[:p, p:] @ isq(S[p:, p:]), compute_uv=False), S
# C1 scalar Dodgson; CMI = -1/2 log(1-r^2) = 1/2 log(1+y)
M = spd(6); a, B, c = [0], [1, 2, 3, 4], [5]; x = D(M, a+B, B+c)
print("C1", D(M, a+B+c)*D(M, B) - D(M, a+B)*D(M, B+c) + x**2, cmi(M, a, B, c), -0.5*np.log(1-x**2/(D(M, a+B)*D(M, B+c))), 0.5*np.log1p(x**2/(D(M, a+B+c)*D(M, B))))
# C2-C3 block identity; naive ratio = prod rho^2; remainder = signed sum of products of cross minors
M = spd(7); a, B, c = [0, 1], [2, 3, 4], [5, 6]; rho, _ = canon(M, a, B, c)
print("C2", D(M, a+B+c)*D(M, B)/(D(M, a+B)*D(M, B+c)*np.prod(1-rho**2)) - 1, cmi(M, a, B, c), -0.5*np.log(1-rho**2).sum())
print("C3", D(M, a+B, B+c)**2/(D(M, a+B)*D(M, B+c)), np.prod(rho**2), -0.5*np.log(1-np.prod(rho**2)))
ac = a+c; t = [-(-1)**(3+sum(j+1 for j in J))*D(M, B+a, B+[ac[j] for j in J])*D(M, B+c, B+[ac[j] for j in range(4) if j not in J])
              for J in itertools.combinations(range(4), 2) if J != (0, 1)]
print("C3", D(M, a+B)*D(M, B+c) - D(M, a+B+c)*D(M, B), sum(t), np.sign(t))
# C4 det M[aB|Bc] = 0 with CMI > 0; Markov iff rank M[aB|Bc] = |B|
for X in (np.diag([0.6, 0.0]), np.zeros((2, 2))):
    MB = spd(3); C = rng.standard_normal((4, 3)); M = np.zeros((7, 7)); I = [0, 1, 5, 6]
    M[np.ix_(I, I)] = np.block([[np.eye(2), X], [X.T, np.eye(2)]]) + C @ np.linalg.solve(MB, C.T); M[np.ix_(I, B)] = C; M[np.ix_(B, I)] = C.T; M[np.ix_(B, B)] = MB
    print("C4", D(M, a+B, B+c), np.linalg.matrix_rank(M[np.ix_(a+B, B+c)], tol=1e-9), cmi(M, a, B, c))
# C5 gaussoid trinomials; MTP2 faithfulness; an unfaithful signed 4-cycle
M = spd(5); worst = 0
for i, j in itertools.combinations(range(5), 2):
    for r in range(4):
        for K in itertools.combinations([v for v in range(5) if v not in (i, j)], r):
            K = list(K); worst = max(worst, abs(D(M, [i]+K, [j]+K)**2 - D(M, [i]+K)*D(M, [j]+K) + D(M, K)*D(M, [i, j]+K)))
            for k in K:
                L = [v for v in K if v != k]
                worst = max(worst, abs(D(M, L)*D(M, [i, k]+L, [j, k]+L) - D(M, [k]+L)*D(M, [i]+L, [j]+L) + D(M, [i]+L, [k]+L)*D(M, [j]+L, [k]+L)))
E = [(0, 1), (1, 2), (2, 3), (3, 0), (2, 4), (4, 5)]; A = np.zeros((6, 6))
for u, v in E: A[u, v] = A[v, u] = rng.uniform(0.3, 1)
M = np.linalg.inv((np.linalg.eigvalsh(A).max() + 0.5)*np.eye(6) - A); adj = {v: {w for e in E for w in e if v in e and w != v} for v in range(6)}
def sep(i, j, K):
    seen, st = {i}, [i]
    while st:
        for w in adj[st.pop()] - set(K) - seen:
            if w == j: return False
            seen.add(w); st.append(w)
    return True
res = [(D(M, [i]+list(K), [j]+list(K))/np.sqrt(D(M, [i]+list(K))*D(M, [j]+list(K))), sep(i, j, K)) for i, j in itertools.combinations(range(6), 2)
       for r in range(5) for K in itertools.combinations([v for v in range(6) if v not in (i, j)], r)]
print("C5", worst, len(res), sum(pc < -1e-12 for pc, _ in res), sum((abs(pc) < 1e-10) != s_ for pc, s_ in res), np.linalg.inv(np.array([[1, .3, 0, .3], [.3, 1, .3, 0], [0, .3, 1, -.3], [.3, 0, -.3, 1]]))[0, 2])
# C6 chordal gap = KL to the Markov projection = sum of CMIs; C7 Koteljanskii
M = spd(6); Cl = [[0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 5]]; Sp = [[1, 2], [2, 3], [3]]
def gap(M):
    ld = lambda S: np.linalg.slogdet(M[np.ix_(S, S)])[1]; Kh = np.zeros((6, 6)); H = Cl[0]; tot = 0
    for C_ in Cl: Kh[np.ix_(C_, C_)] += np.linalg.inv(M[np.ix_(C_, C_)])
    for S_ in Sp: Kh[np.ix_(S_, S_)] -= np.linalg.inv(M[np.ix_(S_, S_)])
    for C_, S_ in zip(Cl[1:], Sp): Hn = sorted(set(H) | set(C_)); tot += 0.5*(ld(H) + ld(C_) - ld(S_) - ld(Hn)); H = Hn
    KM = Kh @ M; return 0.5*(sum(map(ld, Cl)) - sum(map(ld, Sp)) - ld(list(range(6)))), 0.5*(np.trace(KM) - 6 - np.linalg.slogdet(KM)[1]), tot, Kh
g = gap(M); print("C6", g[:3], gap(np.linalg.inv(g[3]))[:3])
M = spd(7); sets = [sorted(set(rng.choice(7, rng.integers(1, 7), replace=False).tolist())) for _ in range(4000)]
print("C7", min(np.log(D(M, S)*D(M, T)/(D(M, sorted(set(S) | set(T)))*D(M, sorted(set(S) & set(T))))) for S, T in zip(sets[::2], sets[1::2])))
# C18 conditional CCA: the top-k canonical variates are the optimal rank-k memory
M = spd(10); a, B, c = [0, 1, 2, 3], [4, 5, 6], [7, 8, 9]; rho, S = canon(M, a, B, c); U = np.linalg.svd(isq(S[:4, :4]) @ S[:4, 4:] @ isq(S[4:, 4:]))[0]
def resid(V):   # I(a:c | B, V^T a) = I(a:c|B) - I(V^T a : c | B)
    k = V.shape[1]; T = np.zeros((k+6, 10)); T[:k, :4] = V.T; T[k:, 4:] = np.eye(6); Mt = T @ M @ T.T
    return cmi(M, a, B, c) - cmi(Mt, list(range(k)), list(range(k, k+3)), list(range(k+3, k+6)))
print("C18", resid(isq(S[:4, :4]) @ U[:, :1]), -0.5*np.log(1-rho[1:]**2).sum(), min(resid(rng.standard_normal((4, 1))) for _ in range(2000)))
# C19 quasi-free fermions: the Gaussian (log-det) Markov condition does not make the fermionic CMI vanish
h = lambda v: float(np.sum(-v*np.log(v) - (1-v)*np.log(1-v))); Sf = lambda C, X: h(np.linalg.eigvalsh(C[np.ix_(X, X)]))
Q, _ = np.linalg.qr(rng.standard_normal((4, 4)) + 1j*rng.standard_normal((4, 4))); C = Q @ np.diag(rng.uniform(0.1, 0.9, 4)) @ Q.conj().T
v = (C[np.ix_([0], [1, 2])] @ np.linalg.solve(C[np.ix_([1, 2], [1, 2])], C[np.ix_([1, 2], [3])]))[0, 0]; C[0, 3] = v; C[3, 0] = np.conj(v)
print("C19", np.linalg.eigvalsh(C)[[0, -1]], Sf(C, [0, 1, 2]) + Sf(C, [1, 2, 3]) - Sf(C, [1, 2]) - Sf(C, [0, 1, 2, 3]))
# C21 Kenyon-Pemantle charts: path Markov chains are the stratum {contiguous faces of lag >= 2 = 0}; MTP2 is not an orthant
K = np.diag(rng.uniform(1.5, 2.5, 6)) + (lambda e: np.diag(e, 1) + np.diag(e, -1))(rng.uniform(-0.6, 0.6, 5)); M = np.linalg.inv(K)
r12, r23, r13_2 = 0.9, 0.1, 0.9; r13 = r13_2*np.sqrt((1-r12**2)*(1-r23**2)) + r12*r23; P = np.linalg.inv(np.array([[1, r12, r13], [r12, 1, r23], [r13, r23, 1]]))
print("C21", max(abs(D(M, list(range(i, j)), list(range(i+1, j+1)))) for i in range(6) for j in range(i+2, 6)), -P[1, 2]/np.sqrt(P[1, 1]*P[2, 2]))
# C21b Kenyon-Pemantle charts = rhombus tilings (generated by cube moves); a graphical model is a coordinate stratum of a chart
#      iff some chart's conditioning set S_ij separates i and j for every non-edge ij (Prop A4'). Exhaustive for n <= 5; the prism at n = 6.
from collections import deque
def charts(n):
    std = frozenset((i, j, frozenset(range(i+1, j))) for i in range(n) for j in range(i+1, n)); seen = {std}; dq = deque([std])
    while dq:
        T = dq.popleft()
        for i, j, k in itertools.combinations(range(n), 3):
            for S in {f[2] for f in T} - {f[2] for f in T if {i, j, k} & f[2]}:
                A_ = {(i, j, S), (j, k, S), (i, k, S | {j})}; B_ = {(i, k, S), (j, k, S | {i}), (i, j, S | {k})}
                for X, Y in ((A_, B_), (B_, A_)):
                    if X <= T and (U := frozenset((T - X) | Y)) not in seen: seen.add(U); dq.append(U)
    return [{(i, j): S for i, j, S in T} for T in seen]
def realised(n, E, CH):
    adj = {v: {w for e in E for w in e if v in e and w != v} for v in range(n)}
    def sp(i, j, K):
        seen, st = {i}, [i]
        while st:
            for w_ in adj[st.pop()] - K - seen:
                if w_ == j: return False
                seen.add(w_); st.append(w_)
        return True
    return any(all(sp(i, j, T[(i, j)]) for i, j in itertools.combinations(range(n), 2) if (i, j) not in E) for T in CH)
def classes(n):   # connected graphs on n vertices up to isomorphism, as (edge list, realised by some labelling?)
    P = list(itertools.combinations(range(n), 2)); perms = list(itertools.permutations(range(n))); CH = charts(n); out = {}
    for mask in range(1 << len(P)):
        E = [P[k] for k in range(len(P)) if mask >> k & 1]; adj = {v: {w for e in E for w in e if v in e and w != v} for v in range(n)}
        seen, st = {0}, [0]
        while st:
            for w_ in adj[st.pop()] - seen: seen.add(w_); st.append(w_)
        if len(seen) < n: continue
        key_ = min(sum(1 << P.index(tuple(sorted((p[u], p[v])))) for u, v in E) for p in perms)
        out[key_] = out.get(key_, False) or realised(n, set(E), CH)
    return len(CH), len(out), sum(out.values())
print("C21b", [classes(n) for n in (3, 4, 5)])
CH6 = charts(6); prism = [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (0, 3), (1, 4), (2, 5)]
print("C21b", len(CH6), any(realised(6, {tuple(sorted((p[u], p[v]))) for u, v in prism}, CH6) for p in itertools.permutations(range(6))))
```

Output (each line is labelled by its check; C4's second line is the Markov case, C21b lists (charts, connected graphs, realised graphs) for n = 3, 4, 5 and then the prism test at n = 6):

```
C1 1.6263032587282567e-19 0.26323718866065815 0.2632371886606589 0.2632371886606586
C2 -5.551115123125783e-16 0.3257586949072381 0.3257586949072378
C3 0.02419519527738521 0.024195195277385284 0.012246353865999586
C3 0.0002870190631535566 0.0002870190631535563 [ 1.  1.  1. -1. -1.]
C4 -1.1970534667108688e-16 4 0.22314355131421168
C4 1.0168400978048905e-30 3 -4.4464432136239495e-14
C5 1.3322676295501878e-15 240 0 0 -1.2388464238195343e-17
C6 (np.float64(1.0752524561019645), np.float64(1.0752524561019652), np.float64(1.0752524561019643)) (np.float64(0.0), np.float64(1.6023737137301802e-31), np.float64(0.0))
C7 0.0
C18 0.21776083504293137 0.21776083504293106 0.2229244272081447
C19 [0.30919386 0.69551248] 0.0012302363544853812
C21 1.1144416595990837e-18 -0.8691865541560537
C21b [(2, 2, 2), (8, 6, 6), (62, 21, 21)]
C21b 908 False
```

**Script B** (matchings, signs, tropical and cluster checks: C8–C17, C20).

```python
# Matchings, signs, tropical and cluster checks C8-C17, C20 (numpy + fractions only)
import numpy as np, itertools, math, random
from fractions import Fraction
rng = np.random.default_rng(2)
def per(A):   # Ryser
    n = len(A); return (-1)**n*sum((-1)**len(S)*np.prod(A[:, list(S)].sum(1)) for r in range(1, n+1) for S in itertools.combinations(range(n), r)) if n else 1.0
# C8 Heawood: spectrum; a signing with |det| = per; K33 has none
N = np.zeros((7, 7)); [N.__setitem__((i, (i+d) % 7), 1) for i in range(7) for d in (0, 1, 3)]
print(np.round(np.linalg.eigvalsh(np.block([[0*N, N], [N.T, 0*N]])), 6))
def maxdet(Nb):   # all signings up to row/column gauge: fix + on a spanning tree
    m = len(Nb); E = list(zip(*np.nonzero(Nb))); seen = {0}; tree = set(); st = [0]
    while st:
        u = st.pop()
        for (i, j) in E:
            if (i == u and m+j not in seen) or (m+j == u and i not in seen):
                w = m+j if i == u else i; seen.add(w); st.append(w); tree.add((i, j))
    free = [e for e in E if e not in tree]
    return max(abs(np.linalg.det(np.where(np.isin(np.arange(m*m), [i*m+j for (i, j), s in zip(free, sg) if s < 0]).reshape(m, m), -Nb, Nb))) for sg in itertools.product([1, -1], repeat=len(free)))
print(per(N), maxdet(N), per(np.ones((3, 3))), maxdet(np.ones((3, 3))))
# C9 Bethe permanent by mirror descent / Sinkhorn; regular closed form; I (x) J2 tightness
def sk(X, it=50):
    for _ in range(it): X = X/X.sum(1, keepdims=True); X = X/X.sum(0, keepdims=True)
    return X
def bethe(A):
    P = sk(A.copy()); mk = A > 0
    for _ in range(3000): P = sk(np.where(mk, np.sqrt(P*A/np.clip(1-P, 1e-300, None)), 0))
    Pm, Am = P[mk], A[mk]; Q = 1-Pm; return math.exp(-(Pm*np.log(Pm/Am)).sum() + (Q*np.log(np.where(Q > 0, Q, 1))).sum())
while True:
    A = sum(np.eye(10)[rng.permutation(10)] for _ in range(3))
    if A.max() == 1: break
print(bethe(A), (4/3)**10, per(A), bethe(np.kron(np.eye(3), np.ones((2, 2)))), per(np.kron(np.eye(3), np.ones((2, 2)))))
r = []
for _ in range(20):
    n = int(rng.integers(3, 8)); X = rng.uniform(0, 1, (n, n))*(rng.uniform(0, 1, (n, n)) < 0.8)
    if per(X) > 0: b = bethe(X); r.append((per(X)/b, per(X)/b/2**(n/2)))
print(min(x for x, _ in r), max(y for _, y in r))   # per >= per_B and per <= 2^(n/2) per_B
th = (np.arange(400000) + 0.5)*np.pi/400000
for d in (3, 4):   # Kesten-McKay: all matchings (tree value) and perfect matchings (Schrijver) at large girth
    x = 2*np.sqrt(d-1)*np.cos(th); w = d*4*(d-1)*np.sin(th)**2/(2*np.pi*(d*d - x*x))*np.pi/400000; r_ = (-1 + math.sqrt(4*d-3))/(2*(d-1))
    print((w*0.5*np.log1p(x*x)).sum(), math.log(1 + d*r_) - d/2*math.log(1 + r_*r_), (w*np.log(abs(x))).sum(), 0.5*math.log((d-1)**(d-1)/d**(d-2)))
# C10 Godsil-Gutman on the Petersen graph, Heilmann-Lieb, an MSS signing
E = [(i, (i+1) % 5) for i in range(5)] + [(i, i+5) for i in range(5)] + [(5+i, 5+(i+2) % 5) for i in range(5)]
m = [sum(len({v for e in S for v in e}) == 2*k for S in itertools.combinations(E, k)) for k in range(6)]
mu = np.zeros(11); mu[0::2] = [(-1)**k*m[k] for k in range(6)]
def As(sg): A = np.zeros((10, 10)); [A.__setitem__((u, v), s) or A.__setitem__((v, u), s) for (u, v), s in zip(E, sg)]; return A
ev = [np.linalg.eigvalsh(As(sg)) for sg in itertools.product([1, -1], repeat=15)]
print(m, np.abs(np.mean([np.poly(e) for e in ev], 0) - mu).max(), np.roots(mu).real.max(), min(e[-1] for e in ev), 2*math.sqrt(2))
# C11 Kasteleyn face rule on a 4x5 grid with one interior edge removed (a hexagonal face)
m_, n_ = 4, 5; V = [(i, j) for i in range(m_) for j in range(n_)]
E = [((i, j), (i, j+1)) for i in range(m_) for j in range(n_-1)] + [((i, j), (i+1, j)) for i in range(m_-1) for j in range(n_) if (i, j) != (1, 2)]
F = []
for i in range(m_-1):
    cols = [j for j in range(n_) if ((i, j), (i+1, j)) in E]
    F += [[((i, j1), (i+1, j1)), ((i, j2), (i+1, j2))] + [((r, j), (r, j+1)) for r in (i, i+1) for j in range(j1, j2)] for j1, j2 in zip(cols, cols[1:])]
def solve2(F, rule):   # GF(2): sum of sign bits around a face of length 2l = rule(l)
    R = [[int(e in f) for e in E] + [rule(len(f)//2)] for f in F]; piv = []; r = 0
    for col in range(len(E)):
        k = next((i for i in range(r, len(R)) if R[i][col]), None)
        if k is None: continue
        R[r], R[k] = R[k], R[r]; R = [[x ^ y for x, y in zip(Ri, R[r])] if i != r and Ri[col] else Ri for i, Ri in enumerate(R)]; piv.append(col); r += 1
    x = [0]*len(E); [x.__setitem__(col, R[i][-1]) for i, col in enumerate(piv)]; return x
blk = [v for v in V if sum(v) % 2 == 0]; wht = [v for v in V if sum(v) % 2]; w = rng.uniform(0.5, 2, len(E))
def K(x):
    M = np.zeros((10, 10))
    for (u, v), s, wt in zip(E, x, w): b, c = (u, v) if sum(u) % 2 == 0 else (v, u); M[blk.index(b), wht.index(c)] = (-1)**s*wt
    return M
print(per(np.abs(K([0]*len(E)))), abs(np.linalg.det(K(solve2(F, lambda l: (l+1) % 2)))), abs(np.linalg.det(K(solve2(F, lambda l: 1)))))
Vh = [(i, j) for i in range(4) for j in range(6)]; Eh = [((i, j), (i, j+1)) for i in range(4) for j in range(5)] + [((i, j), (i+1, j)) for i in range(3) for j in range(6) if (i+j) % 2 == 0]
bh = [v for v in Vh if sum(v) % 2 == 0]; wh = [v for v in Vh if sum(v) % 2]; Kh = np.zeros((12, 12))
for (u, v) in Eh: b_, c_ = (u, v) if sum(u) % 2 == 0 else (v, u); Kh[bh.index(b_), wh.index(c_)] = rng.uniform(0.5, 2)
print(per(Kh), abs(np.linalg.det(Kh)))   # honeycomb (brick wall): all faces hexagons, all signs + already Kasteleyn
# C12 Lieb-Loss VII(B): four squares, weak spokes at the centre; flux pi in every square is not energy-minimising
def energy(t, phi):   # 3x3 grid, horizontal bonds real, vertical bonds carry Peierls phases; phi = fluxes of the 4 squares
    T = np.zeros((9, 9), complex); k = 0
    for i in range(3):
        for j in range(2): T[3*i+j, 3*i+j+1] = t[k]; k += 1
    for i in range(2):
        th = 0
        for j in range(3): T[3*i+3+j, 3*i+j] = t[k]*np.exp(1j*th); k += 1; th += phi[2*i+j] if j < 2 else 0
    e = np.linalg.eigvalsh(T + T.conj().T); return e[e < 0].sum()
def egrid(m, phi):   # m x m grid, unit hoppings, flux phi[i][j] through square (i,j)
    T = np.zeros((m*m, m*m), complex)
    for i in range(m):
        for j in range(m-1): T[m*i+j, m*i+j+1] = 1
    for i in range(m-1):
        th = 0
        for j in range(m): T[m*i+m+j, m*i+j] = np.exp(1j*th); th += phi[i][j] if j < m-1 else 0
    e = np.linalg.eigvalsh(T + T.conj().T); return e[e < 0].sum()
print(egrid(4, [[np.pi]*3]*3), egrid(4, [[0]*3]*3), min(egrid(4, rng.uniform(0, 2*np.pi, (3, 3)).tolist()) for _ in range(1000)))
t = np.ones(12); t[[2, 3, 7, 10]] = 0.2
print(energy(t, [np.pi]*4), energy(t, [np.pi, np.pi, np.pi, 0]))
# C13 Kuo condensation and C14 urban renewal on a weighted 4x4 grid (perfect matchings via permanents)
def PM(wt, rm=()):
    B_ = sorted({b for b, _ in wt} - set(rm), key=str); W_ = sorted({c for _, c in wt} - set(rm), key=str)
    if len(B_) != len(W_): return 0.0
    A = np.zeros((len(B_), len(W_)))
    for (b, c), x in wt.items():
        if b in B_ and c in W_: A[B_.index(b), W_.index(c)] = x
    return per(A)
col = lambda v: (v[0]+v[1]) % 2 if v[0] != 'n' else 1 - (v[1][0]+v[1][1]) % 2
key = lambda u, v: (u, v) if col(u) == 0 else (v, u)
wt = {key((i, j), (i+di, j+dj)): rng.uniform(0.5, 2) for i in range(4) for j in range(4) for di, dj in ((0, 1), (1, 0)) if i+di < 4 and j+dj < 4}
a, b, c, d = (0, 0), (0, 3), (3, 3), (3, 0)
print(PM(wt)*PM(wt, (a, b, c, d)), PM(wt, (a, b))*PM(wt, (c, d)) + PM(wt, (a, d))*PM(wt, (b, c)))
sq = [(1, 1), (1, 2), (2, 2), (2, 1)]; ws = [wt[key(sq[k], sq[(k+1) % 4])] for k in range(4)]; Dl = ws[0]*ws[2] + ws[1]*ws[3]
wt2 = {k: x for k, x in wt.items() if k not in [key(sq[k_], sq[(k_+1) % 4]) for k_ in range(4)]}
nv = [('n', v) for v in sq]
for k in range(4): wt2[key(sq[k], nv[k])] = 1.0; wt2[key(nv[k], nv[(k+1) % 4])] = ws[(k+2) % 4]/Dl
print(PM(wt), Dl*PM(wt2))
# C15 tropical limit of A=[[1,2],[2,3]]: both permutations weigh 4; T log per -> 4, but det(exp(A/T)) = 0 for every T
w = np.array([1 + 3.0, 2 + 2.0]); s_ = np.array([1, -1])   # weights and signs of the two permutations
for T in (0.1, 0.01): print(w.max() + T*np.log(np.exp((w - w.max())/T).sum()), (s_*np.exp((w - w.max())/T)).sum())
# C16 Rogers dilogarithm over the A2 period
Li2 = lambda x: sum(x**k/k**2 for k in range(1, 400)) if x <= 0.5 else math.pi**2/6 - math.log(x)*math.log(1-x) - Li2(1-x)
L = lambda x: Li2(x) + 0.5*math.log(x)*math.log(1-x); y = [0.7, 3.1]
for k in range(3): y.append((1+y[-1])/y[-2])
print(6/math.pi**2*sum(L(v/(1+v)) for v in y), 6/math.pi**2*sum(L(1/(1+v)) for v in y))
phi = (1+5**0.5)/2; print(6/math.pi**2*L(1/(1+phi)))   # constant solution y = phi: per-step value 2/5
# C17 T-system on A3 x A3: T(i,j,u+h+h') = T(4-i,4-j,u), exact rationals
mm = nn = 3; H = mm+nn+2; rnd = random.Random(0)
T = {(i, j, u): Fraction(rnd.randint(1, 9), rnd.randint(1, 9)) for i in range(1, mm+1) for j in range(1, nn+1) for u in (0, 1) if (i+j+u) % 2 == 0}
for u in range(1, 2*H+2):
    for i in range(1, mm+1):
        for j in range(1, nn+1):
            if (i+j+u+1) % 2 == 0:
                T[(i, j, u+1)] = (math.prod(T[(k, j, u)] for k in (i-1, i+1) if 1 <= k <= mm) + math.prod(T[(i, k, u)] for k in (j-1, j+1) if 1 <= k <= nn))/T[(i, j, u-1)]
print(all(T[(i, j, u+H)] == T[(mm+1-i, nn+1-j, u)] for (i, j, u) in list(T) if u < 2), all(T[(i, j, u+H)] == T[(i, j, u)] for (i, j, u) in list(T) if u < 2))
# C20 pentagon E(V)E(U) = E(U)E(-UV)E(V) when VU = qUV, E(x) = sum x^n/(q;q)_n, truncated at degree 12
DEG, qq = 12, 0.3+0.5j
def mul(X, Y): 
    Z = {}
    for (a1, b1), x in X.items():
        for (c1, d1), y_ in Y.items():
            if a1+b1+c1+d1 <= DEG: Z[(a1+c1, b1+d1)] = Z.get((a1+c1, b1+d1), 0) + x*y_*qq**(b1*c1)
    return Z
def E_(X):
    R, P, poch = {}, {(0, 0): 1.0}, 1.0
    for k in range(DEG+1):
        if k: P = mul(P, X); poch *= 1 - qq**k
        for kk, v in P.items(): R[kk] = R.get(kk, 0) + v/poch
    return R
lhs = mul(E_({(0, 1): 1}), E_({(1, 0): 1})); rhs = mul(mul(E_({(1, 0): 1}), E_({(1, 1): -1})), E_({(0, 1): 1}))
print(max(abs(lhs.get(k, 0) - rhs.get(k, 0)) for k in set(lhs) | set(rhs)))
```

Output, in order:
- C8: the spectrum (printed over two lines), then per and max |det| for the Heawood graph and K₃,₃;
- C9: regular graph [per_B, (4/3)^10, per] and I⊗J₂ [per_B, per]; random bounds [min per/per_B, max per/(2^{n/2} per_B)]; Kesten–McKay for d = 3 and 4 [all-matchings integral, tree value, ln|x| integral, Schrijver];
- C10: [m_k, Godsil–Gutman residual, max root, min λ_max, 2√2];
- C11: hexagon grid [count, face rule, all −1], then honeycomb [count, |det|];
- C12: 4 × 4 grid [π, 0, best random], then Lieb–Loss [π⁴, (π,π,π,0)];
- C13; C14; C15 [T log per, det] at T = 0.1 and 0.01;
- C16: period sums, then the constant solution;
- C17: [rotation, identity];
- C20: the residual.

```
[-3.       -1.414214 -1.414214 -1.414214 -1.414214 -1.414214 -1.414214
  1.414214  1.414214  1.414214  1.414214  1.414214  1.414214  3.      ]
24.0 23.999999999999993 6.0 4.0
17.75772663381251 17.757726633812588 64.0 1.0 8.0
1.372059886216939 0.7376817330562153
0.5815754049028405 0.5815754049028404 0.14384334671649238 0.14384103622589042
0.6613554692676729 0.6613554692676732 0.26162667118420074 0.26162407188227393
[1, 15, 75, 145, 90, 6] 1.791795878336444e-14 2.63144517283292 2.000000000000001 2.8284271247461903
1481.711869063205 1481.7118690630593 1429.996123327475
1.2760626385133946 1.2760626383655738
-13.252758550612269 -10.94427190999916 -12.970749074898766
-5.22842712474619 -5.316753705465365
1176.1907681936414 1176.1907681936377
149.85988508344963 149.85988508344997
4.0693147180559945 0.0
4.0069314718056 0.0
3.0 2.0000000000000004
0.3999999999999999
True False
3.049565105563204e-16
```

Not reproduced in the scripts (session runs only):
- C17 for the six other pairs (m, n);
- the n = 6 exhaustive stratum search of C21 (86 s);
- the Jacobian-rank check of C21;
- C12's 2 × 3 ladder search.

They use the same functions as the scripts.

---

## 7. Where this note refines Part 2

Each item names the place in [foundations-unlocks.md](foundations-unlocks.md) to update.

- **Unlock 24, last bullet, and its Appendix A row** ("Kasteleyn's face rule equals Lieb's flux rule", THEOREM). Refine to: the face rule is the |det|-maximising canonical flux on every planar bipartite graph (Lieb–Loss Thm 3.1). It minimises the half-filled energy on Lieb's lattices, but not in general (Lieb–Loss §VII(B); C12).
- **Unlock 32, first bullet, and its Appendix A row** ("tropically per = det = assignment", THEOREM). Refine to: they are equal as formal tropicalisations. The zero-temperature limit of log|det| falls below the assignment value, or to −∞, when top-weight permutations cancel (C15). This agrees with Part 2's own unlock 29.
- **Appendix A, genus count "from memory"**: now verified (Cimasoni–Reshetikhin: 2^{2g} Pfaffians, Arf signs, an H¹(Σ; Z/2)-torsor).
- **Appendix A, Ramanujan 2-lifts "from memory"**: now verified (MSS I §5) and checked (C10).
- **Unlocks 33, 45, 46.**
  - For vector-valued blocks, the CMI is −½ Σ log(1 − ρᵢ²) (Theorem A), not ½ log(1 + y) with a single cross minor. The memory budget of unlock 45 at a block separator is Theorem A.
  - Theorem A′ adds the optimal rank-k memory.
- **Appendix A, "Idea A: a cluster structure … SPECULATION"**: now partly DERIVED.
  - Kenyon–Pemantle cube moves change separators.
  - Graphical models are coordinate strata under a criterion that covers all chordal graphs on ≤ 6 vertices (Prop A4′).
  - CONJECTURE A4\* remains, and MTP₂ is not the positive part.
- **Unlock 47(b)** ("expanders … are the best case for trees"). Add the hypotheses: unfrustrated models for free energies (T3), and uniqueness for marginals (T4). The competition scores marginals.
- **Appendix A, "Idea C … SPECULATION; not used"**: the construction is THEOREM-level (KN; §4). It is still not used.
- **Consistent, no change needed**: unlock 14(a) (lim sup), unlock 26, unlock 34 (graphoid hypothesis), and Appendix A's JSV row ("JSV's crossing rests on positivity").

---

## Appendix A. Every claim of the input, with its verdict

The verdict counts in the Summary are taken from this table.
- I1–I50 (input §§1–4): 38 CORRECT, 11 IMPRECISE, 1 UNVERIFIED.
- I51–I79 (ideas and "how much to trust this"): 20 CORRECT, 6 IMPRECISE, 1 UNVERIFIED, 2 OPEN/CONJECTURE.

| # | claim (abridged) | verdict | evidence | streams |
|---|---|---|---|---|
| I1 | monomer–dimer: no phase transition (real roots) | CORRECT | Heilmann–Lieb abstract; ACFK Thm 3.3 | B, M |
| I2 | its growth rate is local for bounded degree | CORRECT (estimability, uniform form proved) | ACFK Thm 1.2; C9 | B |
| I3 | perfect matchings not local in general | CORRECT | ACFK Thms 1.7–1.8 | B, M |
| I4 | an edge in all but a cⁿ fraction of perfect matchings | CORRECT | ACFK Thm 1.7 | B, M |
| I5 | a global rigidity no local view detects | CORRECT | ACFK Thm 1.8 | B |
| I6 | expanders: no p(e) exponentially close to 1 | CORRECT (d-regular bipartite δ-expanders) | ACFK Thm 1.9 | B |
| I7 | growth rate local via Gamarnik–Katz | CORRECT | ACFK §1.3; GK Cor 1, Thm 2 | B |
| I8 | Schrijver's constant per vertex at large girth | CORRECT (no expansion needed) | ACFK Thm 1.5; C9 | B, M |
| I9 | that constant is the Bethe permanent | CORRECT | proof (Vontobel convexity); C9 | B, M |
| I10 | Bethe = tree Markov factorisation | CORRECT | Vontobel Cor 15; T1 | B, K |
| I11 | per_B ≤ per ≤ 2^{n/2} per_B | CORRECT | Anari–Rezaei Thms 3–4; C9 | B, M |
| I12 | Bethe permanent = limit over covers | IMPRECISE: lim sup | Vontobel Thm 39 | B |
| I13 | pseudorandom geometry: global count from local tree factorisation | IMPRECISE: growth rate, not count | I11 | B, M |
| I14 | Kasteleyn: signing turns per into det on planar graphs | IMPRECISE: bipartite per → det; general haf → Pf | Lieb–Loss Thm 3.1 | M |
| I15 | genus g: 4^g determinants, one per spin structure | CORRECT (Pfaffians) | Cimasoni–Reshetikhin | M |
| I16 | the obstruction is H¹(Σ; Z/2) | IMPRECISE: H¹ indexes (torsor), not an obstruction | CR Thm 3.2 | M |
| I17 | RST/McCuaig: planar pieces glued, plus one exception | CORRECT | RST Thm 1.3/6.8, §7.2 | M |
| I18 | Heawood, Fano incidence graph, eigenvalues ±√2 | CORRECT | C8 | M |
| I19 | determinants fermions, permanents bosons | CORRECT | TD; AA | M, H |
| I20 | matchgates classically easy | CORRECT | Terhal–DiVincenzo | M, H |
| I21 | boson sampling believed hard | CORRECT | Aaronson–Arkhipov | M |
| I22 | face rule: −1 on 4k-faces, +1 on 4k+2-faces | CORRECT | C11 | M |
| I23 | coincides with Lieb's flux rule | IMPRECISE: \|det\|-maximiser; energy only on Lieb's lattices | Lieb–Loss Thm 3.1, §VII(B); Lieb 1994; C12 | M, H |
| I24 | the sign is a ground-state gauge field | IMPRECISE: \|det\|-maximising field | as I23 | M, H |
| I25 | Godsil–Gutman | CORRECT | MSS I Thm 3.6; C10 | M, B |
| I26 | roots ≤ 2√(d−1) | CORRECT | Heilmann–Lieb via MSS I; C10 | M, B |
| I27 | signing = double cover = Z/2 gauge field | CORRECT | Bilu–Linial via MSS I §5 | M |
| I28 | MSS: bipartite Ramanujan graphs of every degree | CORRECT | MSS I §5; C10 | M, B |
| I29 | the same method solved Kadison–Singer | CORRECT | MSS II | M |
| I30 | pure states on the diagonal extend uniquely | CORRECT | MSS II | M, H |
| I31 | the diagonal is the "classical shadow" | UNVERIFIED (missing context) | Part 2 treats it as an analogy | F, H |
| I32 | Gurvits's van der Waerden proof is the ancestor | CORRECT | MSS II acknowledgement | M |
| I33 | strongest bridge to NCG already a theorem | IMPRECISE: operator algebras, not NCG proper | — | H |
| I34 | tropicalisation: log-sum-exp → max | CORRECT (positive weights) | proof | T |
| I35 | tropically per = det = assignment | IMPRECISE: formally yes, in the limit no | C15 | T, M |
| I36 | at zero temperature signs die | IMPRECISE: formally only | C15; Part 2 unlock 29 | T, F |
| I37 | Chin: tropical proof lifted by positivity and periodicity | CORRECT | Chin Thm 1.1; C17 | T |
| I38 | finite-temperature permanent #P-hard | CORRECT (exact; approximation is FPRAS) | Valiant via Vontobel §I | T, M |
| I39 | Kasteleyn and expansion-based JSV are the crossings | IMPRECISE: JSV works for all nonnegative matrices | Vontobel §I-B; GK | T, M |
| I40 | cluster algebras: tropical-to-full transfer is a theorem | CORRECT | IIKKN; KN Props 2.4, 3.4 | T |
| I41 | cluster variables = snake-graph matchings | CORRECT (ordinary arcs, crossing monomial) | MSW Thm 4.9 | M, T |
| I42 | weighted permanents made determinants by Kasteleyn | CORRECT | planar bipartite | M |
| I43 | octahedron recurrence = sums over perfect matchings | CORRECT | Speyer | M, T |
| I44 | Kuo condensation is an exchange relation | CORRECT (three-term Plücker form) | Kuo Thm 2.1; C13 | M, K |
| I45 | urban renewal preserves the dimer sum up to a factor | CORRECT | Propp; Goncharov–Kenyon; C14 | M |
| I46 | domino shuffling samples exactly | CORRECT | Propp | M |
| I47 | integrable structure replaces expansion | CORRECT (descriptive) | — | M |
| I48 | torus dimers are cluster integrable systems | CORRECT | Goncharov–Kenyon | M, T |
| I49 | planar Ising in OG_{≥0}; same pair as Chin builds on | CORRECT | Galashin–Pylyavskyy Thm 2.3 | M, T |
| I50 | total positivity governs Markov structure; MTP₂ faithful | IMPRECISE: graphoid needed; MTP₂ ≠ matrix TP | Fallat et al. Thm 6.1, Ex 5.4; C5 | K, T |
| I51 | Dodgson identity | CORRECT (scalar a, c); block form = Theorem A | C1–C4 | K |
| I52 | I = −½ log(1 − r²) | CORRECT (scalar); blocks: −½ Σ log(1 − ρᵢ²) | C1–C3 | K |
| I53 | Markov = degeneration to one monomial | CORRECT (scalar); blocks: all \|a\|·\|c\| relations | Theorem A(5); C4 | K |
| I54 | Koteljanskii = CMI ≥ 0 | CORRECT (PD) | Cover–Thomas; C7 | K |
| I55 | chordal Markov ⟺ det = ∏ cliques/∏ separators | CORRECT | proof; C6 | K |
| I56 | the junction-tree formula Bethe imposes everywhere | CORRECT | T1 | B, K |
| I57 | open: cluster structure, Markov strata, mutation = separator change, MTP₂ positive | partly DERIVED (Prop A4′), CONJECTURE A4\*, MTP₂ part false | §2.3; C21 | K, T, M |
| I58 | separator principle strong on amenable geometry | CORRECT (positivity; boundary gain) | S1, S4 | K, B |
| I59 | tree principle strong while correlations decay; glassiness | IMPRECISE: decay needed for marginals, not free energies | T3–T6 | B |
| I60 | expanders worst for separators, best for trees; permanent shows both | IMPRECISE: best only for unfrustrated models | S5, T3, T5 | B, K |
| I61 | substitution towers → plane; covers → tree | IMPRECISE: needs girth → ∞ | §3.3 | B, F |
| I62 | Bethe permanents are limits over covers; Ramanujan from 2-lifts | CORRECT (lim sup) | Vontobel Thm 39; MSS I | B, M |
| I63 | box spaces expander ⟺ profinite spectral gap | CORRECT | proof §3.3 | B, H |
| I64 | the same box spaces break coarse Baum–Connes | CORRECT under stated hypotheses | HLS; Willett–Yu; ℓ^p paper | H |
| I65 | quantum Bethe on random lifts; Yang says little on expanders | CONJECTURE QB; second part CORRECT | §3.4; S5 | B, K, H |
| I66 | quantum mutation = Ad Ψ_q ∘ monomial | CORRECT | KN Prop 3.1; Keller | H, T |
| I67 | layers of commuting local gates (brickwork) | CORRECT (Heisenberg picture) | b = 0 within colours | H |
| I68 | periodicity transfers via tropical data (KN) | CORRECT | KN Props 2.4, 3.4, Thm 3.5 | H, T |
| I69 | exactly periodic up to a phase; half-period a permutation | CORRECT (period; classical half-period); quantum half-period DERIVED | KN Thms 4.6–4.7; Chin; C17 | H, T |
| I70 | classical shadow's value is a central charge | IMPRECISE: (π²/6)·N₋ = period × c_eff | KN Thm 2.7; C16 | T |
| I71 | at irrational q, the quantum torus is Connes' torus | IMPRECISE: different *-structure; mutation not a C*-automorphism | §4.1 | H |
| I72 | golden angle: PV embedding into Effros–Shen [[1,1],[1,0]] | CORRECT | Pimsner–Voiculescu; Rieffel | H |
| I73 | Morita equivalent to the Penrose algebra | CORRECT | Connes; Elliott + BGR (from memory) | H |
| I74 | "Established: sections 1–4" | IMPRECISE: 11 imprecisions (§1) | Appendix A | all |
| I75 | the pentagon generalisation is an immediate corollary of Chin plus Keller | UNVERIFIED (claim not located in either paper) | — | T, H |
| I76 | Idea A identity elementary and correct; cluster structure open | CORRECT (scalar); open part now partly answered | §2 | K |
| I77 | Idea B a synthesis of proven results; quantum half untested | CORRECT | §3 | B, K |
| I78 | Idea C well defined, interest unproven | CORRECT | §4 | H |
| I79 | "local-to-global" names several mechanisms | CORRECT (descriptive) | §3 | all |

The input's source list includes Newman–Vardi, *FPRAS approximation of the matrix permanent in practice*. No claim depends on it, and it was not opened.

---

## Sources

**Opened and checked in this session** (the statement used was read in the source).
- M. Abért, P. Csikvári, P. Frenkel, G. Kun, *Matchings in Benjamini–Schramm convergent graph sequences*, arXiv:1405.3271. Thms 1.2, 1.5–1.9, 3.3, Remark 3.6, §1.3.
- O. J. Heilmann, E. H. Lieb, *Theory of monomer-dimer systems*, CMP 25 (1972) 190–232. Abstract.
- D. Gamarnik, D. Katz, *A deterministic approximation algorithm for computing the permanent of a 0,1 matrix*, arXiv:math/0702039 (JCSS 76 (2010)). Thm 1 (Bayati–Gamarnik–Katz–Nair–Tetali), Thm 2, Cor 1.
- N. Anari, A. Rezaei, *A tight analysis of Bethe approximation for permanent*, arXiv:1811.02933. Def 2, Thms 3–4, tightness.
- P. O. Vontobel, *The Bethe permanent of a non-negative matrix*, arXiv:1107.4196 (IEEE Trans. Inf. Theory, 2013). Cor 15, Thm 20, Lemma 21 (convexity), §VII-E (the d-regular evaluation), Thm 39 (lim sup over covers), §I (#P, JSV).
- D. Cimasoni, N. Reshetikhin, *Dimers on surface graphs and spin structures I*, arXiv:math-ph/0608070. The 2^{2g} Pfaffian formula, Arf invariants, Thm 3.2.
- N. Robertson, P. D. Seymour, R. Thomas, *Permanents, Pfaffian orientations, and even directed circuits*, Ann. of Math. 150 (1999); arXiv:math/9911268. §§1.1–1.3, 6.3, 6.8, 7.2–7.3.
- E. H. Lieb, M. Loss, *Fluxes, Laplacians and Kasteleyn's theorem*, Duke Math. J. 71 (1993); arXiv:cond-mat/9209031. Thm 3.1, §VII(B), App. A.
- E. H. Lieb, *Flux phase of the half-filled band*, PRL 73 (1994); arXiv:cond-mat/9410025. Hypotheses and statement.
- A. Marcus, D. Spielman, N. Srivastava, *Interlacing families I*, arXiv:1304.4132. Thms 3.1–3.2, 3.6, §5. *Interlacing families II*, arXiv:1306.3969. Kadison–Singer, Thm 1.4, Cor 1.5, Thm 6.1, the credit to Gurvits.
- B. M. Terhal, D. P. DiVincenzo, arXiv:quant-ph/0108010. S. Aaronson, A. Arkhipov, arXiv:1011.3245. Abstracts.
- I. Chin, *Half-periodicity of Zamolodchikov periodic cluster algebras*, arXiv:2602.15140 v2. Thm 1.1, Prop 3.1, Cor 3.5, Thm 2.14 (IIKKN), Lemma 4.3.
- G. Musiker, R. Schiffler, L. Williams, arXiv:0906.0748. Thms 4.9, 4.16, 4.20.
- D. Speyer, arXiv:math/0402452. E. Kuo, arXiv:math/0304090, Thm 2.1. J. Propp, *Generalized domino-shuffling*, arXiv:math/0111034.
- A. B. Goncharov, R. Kenyon, *Dimers and cluster integrable systems*, arXiv:1107.5588.
- P. Galashin, P. Pylyavskyy, *Ising model and the positive orthogonal Grassmannian*, arXiv:1807.03282. Thm 2.3.
- S. Fallat, S. Lauritzen, K. Sadeghi, C. Uhler, N. Wermuth, P. Zwiernik, *Total positivity in Markov structures*, arXiv:1510.01290. Thm 6.1, Ex 5.4, Karlin–Rinott.
- M. Madiman, P. Tetali, arXiv:0901.0044 (Koteljanskii as log-det submodularity). T. Cover, J. Thomas, *Determinant inequalities via information theory*, SIAM J. Matrix Anal. Appl. 9 (1988) 384–392 (abstract).
- T. Boege, A. D'Alì, T. Kahle, B. Sturmfels, *The geometry of gaussoids*, arXiv:1710.07175. Square and edge trinomials, CI ⟺ a_{ij|K} = 0, Thm 5.6, non-realisable gaussoids, Sullivant.
- R. Kenyon, R. Pemantle, *Principal minors and rhombus tilings*, arXiv:1404.1354 (full text). Hexahedron relation, Thms 4.3, 4.4, 5.2, 5.7, Prop 5.4–5.5, Kashaev form (7)–(10).
- S. Chepuri, T. George, D. E. Speyer, *Electrical networks and Lagrangian Grassmannians*, arXiv:2106.15418; Ann. Inst. Henri Poincaré D 13 (2026) 191–216. §5.2, Thm 5.1.
- R. Karpman, *Total positivity for the Lagrangian Grassmannian*, arXiv:1510.04386 (abstract). P. Galashin, T. Lam, *Positroid varieties and cluster algebras* (author's PDF): Thm 3.5, Cor 4.4.
- R. M. Kashaev, T. Nakanishi, *Classical and quantum dilogarithm identities*, arXiv:1104.4630. Props 2.4, 3.1, 3.4, Thms 2.7, 3.5, 4.6, 4.7. B. Keller, *On cluster theory and quantum dilogarithm identities*, arXiv:1102.4148.
- M. Pimsner, D. Voiculescu, *Imbedding the irrational rotation C*-algebra into an AF-algebra*, J. Operator Theory 4 (1980) 201–210 (IMAR scan). M. A. Rieffel, *C*-algebras associated with irrational rotations*, Pacific J. Math. 93 (1981).
- A. Connes, *Noncommutative Geometry* (1994), introduction (the Penrose AF algebra). *On the noncommutative geometry of tilings*, arXiv:1412.5442 (A. Julien, J. Kellendonk, J. Savinien; authors from memory).
- N. Higson, V. Lafforgue, G. Skandalis, *Counterexamples to the Baum–Connes conjecture*, GAFA 12 (2002) 330–354 (introduction). R. Willett, G. Yu, *Higher index theory for certain expanders and Gromov monster groups II*, arXiv:1012.4151. *Expanders are counterexamples to the ℓ^p coarse Baum–Connes conjecture* (EMS Press). Mimura–Ozawa–Sako–Suzuki, *Group approximation in Cayley topology and coarse geometry III* (AGT 15 (2015)).
- M. Brown, D. Poulin, *Quantum Markov networks and commuting Hamiltonians*, arXiv:1206.0755. Thms 1–5. D. Poulin, M. B. Hastings, *Markov entropy decomposition*, arXiv:1012.2050.
- A. Dembo, A. Montanari, *Ising models on locally tree-like graphs*, arXiv:0804.4726. Thms 2.4, 2.6, 2.7.
- A. Sly, N. Sun, *Counting in two-spin models on d-regular graphs*, arXiv:1203.2602. Thms 1, 2, 4, 5.
- Through the programme's digests: Yang, arXiv:2609.38007, Thm II.1 ([arxiv-2609.38007.md](../digests/bridges/arxiv-2609.38007.md)); Chen–Rouzé, arXiv:2504.02208, Cor III.2, Cor B.2 ([arxiv-2504.02208.md](../digests/bridges/arxiv-2504.02208.md)); Hammersley–Clifford and the Moussouris counterexample ([hammersley-clifford.md](../digests/bridges/hammersley-clifford.md)).

**From memory (standard; not re-opened).**
- Valiant (1979); Jerrum–Sinclair–Vigoda (2004); Bilu–Linial (2006); McCuaig (2004).
- Galluccio–Loebl and Tesler (cited in CR); Keller's periodicity theorem; Fock–Goncharov positive representations; Faddeev's Φ_b.
- Elliott's classification of AF algebras; Brown–Green–Rieffel; Lubotzky on property (τ); Cheeger and Alon–Milman.
- Weitz (2006); Pearl and Yedidia–Freeman–Weiss; cluster-expansion locality at high temperature; Kato–Brandão (as tabulated by Yang).
- 1RSB and condensation (Krzakala et al.; Ding–Sly–Sun); the Lee–Yang effective central charge 2/5; boundary measurements and Plücker relations for planar bipartite graphs (Postnikov; Lam).
