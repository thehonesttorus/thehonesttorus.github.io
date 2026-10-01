# Input supplied by the user on 1 Oct 2026 (verbatim)

Supplied as context for the fresh-slate design. It appears to be an assistant's reply from another conversation (it refers to "our groupoid discussion", "last time" and "my first-turn framing"). Nothing in it has been verified here yet; `foundations.md` records the verification.

---

The sharper statement concerns the graph itself, and it is the rigorous form of your intuition:

 All matchings: always local. The monomer–dimer partition function has no phase transition (Heilmann–Lieb: the matching polynomial has only real roots). So its growth rate is fixed by local neighbourhood statistics for every bounded-degree graph (Abért–Csikvári–Frenkel–Kun).
 Perfect matchings: not local in general. The same paper builds d-regular bipartite graphs where one edge lies in all but an exponentially small fraction of perfect matchings. That is a global rigidity no local view can detect.
 Expansion removes the rigidity. On bipartite expanders, no edge probability can be exponentially close to 1, and the growth rate becomes locally determined (via Gamarnik–Katz). On large-girth d-regular bipartite graphs it equals Schrijver's constant ((d−1)^(d−1)/d^(d−2))^(1/2) per vertex.
 That constant is the Bethe permanent. Pretend the graph is a tree and use the exact Markov factorization valid on trees. For every nonnegative matrix, per_B ≤ per ≤ 2^(n/2) · per_B (Gurvits; Anari–Rezaei). The Bethe permanent is also a limit over graph covers (Vontobel).
 On pseudorandom geometry, a global count is computed by a local, tree-shaped Markov factorization. Without that, global constraints can freeze it.

 3. The sign factor is where the deep structure hides
    Topology.

 Kasteleyn: on a planar graph, signing the edges turns the permanent into a determinant.
 On a genus-g surface you need 4^g determinants, one per spin structure. The obstruction is the surface's first Z/2-cohomology.
 Robertson–Seymour–Thomas (and McCuaig) characterized the bipartite graphs that admit such signs: planar pieces glued together, plus one exception. The exception is the Heawood graph, the incidence graph of the Fano plane, whose nontrivial eigenvalues ±√2 beat the Ramanujan bound.
 Physics.

 Determinants are fermions; permanents are bosons. Free-fermion (matchgate) computations are classically easy, while boson sampling is believed hard.
 Kasteleyn's face rule puts sign product −1 around faces of length 4k and +1 around faces of length 4k+2. This coincides with Lieb's flux rule for half-filled free fermions. The sign that tames the permanent is a ground-state gauge field.
 Pseudo-randomness and operator algebras.

 Godsil–Gutman: average det(x − A_s) over random signings s and you get the matching polynomial. Heilmann–Lieb bound its roots by 2√(d−1).
 A signing defines a double cover — a random Z/2 gauge field. Marcus–Spielman–Srivastava turned this into bipartite Ramanujan graphs of every degree.
 The same interlacing-polynomial method then solved the Kadison–Singer problem. Pure states on the diagonal subalgebra of B(ℓ²) extend uniquely. That diagonal is exactly the "classical shadow" from our groupoid discussion.
 Gurvits's real-stable-polynomial proof of van der Waerden's permanent bound is the ancestor of this method.
 So the strongest bridge between pseudo-randomness, permanents and NCG is already a theorem.

 Zero temperature, where signs die. Tropicalization is the zero-temperature limit of a partition function: log-sum-exp becomes max. Tropically, permanent and determinant coincide; both become the assignment problem, solvable in polynomial time.

 Chin's proof works exactly there, with c-vectors and tropical T-systems, then lifts to the full system via positivity and the periodicity theorem. For general permanents, the lift from zero to finite temperature is #P-hard. Planarity (Kasteleyn) and expansion-based sampling (JSV) are the classic ways across. Cluster algebras are a setting where a tropical-to-full transfer is a theorem.

 4. Cluster algebras are already made of matchings
    Cluster variables. In surface cluster algebras they are generating functions of perfect matchings of planar "snake graphs" (Musiker–Schiffler–Williams). These are weighted permanents, made into determinants by Kasteleyn.
    T-systems. Octahedron-recurrence solutions are sums over perfect matchings (Speyer). Kuo condensation, an identity among matching counts with vertices deleted, is an exchange relation.
    Mutation as graph surgery. A mutation on a planar bipartite graph is "urban renewal", a local surgery that preserves the dimer partition function up to a factor. Iterated on the Aztec diamond it becomes domino shuffling, which samples random tilings exactly. Here integrable structure replaces the expansion that Markov-chain sampling needs.
    Other bridges.
    Dimer models on torus graphs are cluster integrable systems (Goncharov–Kenyon).
    Planar Ising models embed in the positive orthogonal Grassmannian (Galashin–Pylyavskyy — the same pair whose periodicity classification Chin builds on).
    Total positivity, where cluster algebras began, governs Markov structure: MTP₂ distributions are faithful to their concentration graphs (Fallat et al.).
 5. Three ideas worth developing (speculative from here)
    A. The Markov property is a degenerate exchange relation. Take a Gaussian with covariance M on a ⊔ B ⊔ c. Dodgson condensation (the determinant form of the octahedron recurrence) reads:

 det M · det M[B] = det M[aB] · det M[Bc] − (det M[aB|Bc])²

 Here M[aB|Bc] is the cross minor with rows aB and columns Bc. The conditional mutual information is:

 I(a:c|B) = −½ log(1 − r²),   r² = (det M[aB|Bc])² / (det M[aB] · det M[Bc])

 So CMI is a function of the ratio between the two terms of an exchange relation. Markovianity is exactly the degeneration of that relation to a single monomial.

 Two classical facts fit this picture:

 Koteljanskii's determinant inequality is the statement CMI ≥ 0.
 A Gaussian is Markov on a chordal graph iff its determinant factorizes as ∏ cliques / ∏ separators. That is the junction-tree formula the Bethe approximation imposes everywhere.
 The open question: is there a cluster structure on Gaussian (or quasi-free fermionic) states in which Markov networks are the boundary strata where exchange relations collapse? Mutation would mean changing a separator, and MTP₂ would be the positive part.

 B. Two local-to-global principles, with pseudo-randomness as the switch. This corrects my first-turn framing:

 Separator-based (Hammersley–Clifford, Yang's CMI bounds): screen a region by its boundary. This is strong on amenable geometry such as lattices and tilings.
 Tree-based (Bethe, cavity, belief propagation): treat neighbourhoods as trees and use the exact Markov factorization there. This is strong on pseudorandom geometry (large girth, expanders) while correlations decay along the tree. It breaks down where expansion creates glassiness.
 Expanders are the worst case for the first principle and the best case for the second, and the permanent shows both behaviours.

 Your "sequence of refinement" has a precise twin: towers of covers.

 Substitution towers converge to the plane, which is amenable. Towers of graph covers converge to a tree, which is not.
 Bethe permanents are limits over covers, and Ramanujan graphs come from 2-lifts.
 Towers of finite quotients ("box spaces") are expanders exactly when the profinite action has a spectral gap. These are the same box spaces whose ghost projections break the coarse Baum–Connes conjecture.
 The quantum frontier is a quantum Bethe theory for Gibbs states on expander interaction graphs, for example via the Markov entropy decomposition on random lifts. It would complement Yang's separator bounds, which say little on expanders.

 C. A quantum circuit that provably returns.

 Construction. Quantize the bipartite belt. On the quantum torus underlying a quantum cluster algebra, each mutation is, up to a monomial change of variables, conjugation by a quantum dilogarithm. Then μ◦ and μ• become layers of commuting local gates on white and then black vertices: a brickwork circuit.
 Periodicity. Periodicity transfers to the quantum Y-system through the same tropical data; this is how Kashaev and Nakanishi extract quantum dilogarithm identities. So this interacting circuit should be exactly periodic up to a phase, and by Chin its half-period acts as a permutation. Its classical shadow is a dilogarithm identity whose value is a central charge.
 Open questions. How do entanglement and CMI evolve under a circuit guaranteed to return?
 Link to last time. At irrational q the quantum torus is Connes' noncommutative torus. At the golden angle it embeds (Pimsner–Voiculescu) in the Effros–Shen AF algebra with matrix [[1,1],[1,0]], which is Morita-equivalent to the Penrose algebra from last time.
 How much to trust this
 Established: sections 1–4. The pentagon generalization is an immediate corollary of Chin plus Keller.
 Idea A: the identity is elementary and certainly correct; whether a cluster structure organizes Markov networks is open.
 Idea B: the two-principle picture is my synthesis of proven results; its quantum half is an untested conjecture.
 Idea C: a well-defined construction whose interest is unproven.
 "Local-to-global" names several genuinely different mechanisms: exact recurrences, inequality cascades, tree factorizations and spectral gaps. The bridges above are real exactly where one mechanism provably feeds another.

 Which do you want to push: A (a concrete computation), B (a precise quantum conjecture you could test numerically), or C (the explicit circuit for the smallest case)?

 Sources:

 Chin, Half-periodicity of Zamolodchikov periodic cluster algebras
 Newman & Vardi, FPRAS Approximation of the Matrix Permanent in Practice
 Abért, Csikvári, Frenkel & Kun, Matchings in Benjamini–Schramm convergent graph sequences
 Fallat et al., Total positivity in Markov structures
 Higson, Lafforgue & Skandalis, Counterexamples to the Baum–Connes conjecture
 Connes, Noncommutative Geometry
