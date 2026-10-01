# The research programme: three prongs

### How the abstract theory, the bridge to networks, and the hunt for unlocks are kept separate and made to talk

*Charter for the work in this repository from 2026-10-01 on. It records the three prongs, what each may take from the others, the rules that stop one prong from distorting another, and where each prong's artefacts live.*

---

## 0. Why three prongs

The framework (abstract barycentres as objects, arrows = the ReLU that produced a face together with the following linear map, functions on arrows forming a C\*-algebra or operator system, projections and states, a simplicial complex read as a nerve) has a precise home in noncommutative geometry (Note 1, Note 2). Two risks follow from that success.

- The general theory could be dragged into something specific and niche by whatever a particular network happens to show. The framework is sensed to be more general than any realization; its statements must remain the kind that are true for a tiling, a Bratteli diagram or a measurement context as much as for a network.
- An experiment on a network could be read as a verdict on the theory, when it is a verdict on the way the network was read into the theory. A negative result most likely means the association of an actual network to the abstract objects was naive, not that the theory is heading the wrong way.

The remedy is a division of labour with explicit interfaces.

| prong | object of study | deliverable | artefacts |
|---|---|---|---|
| **1 — theory** | the abstract objects: faces, arrows, the quiver, its algebras at three resolutions, cocycles, states, conditioning, sufficiency, the nerve and coherence complexes | theorems and definitions, labelled Fact / Derivation / Conjecture / Interpretation | [conditional-arrow-algebra.md](conditional-arrow-algebra.md), [simplicial-complex-as-decomposition.md](simplicial-complex-as-decomposition.md), [checks/graph_algebra_checks.py](checks/graph_algebra_checks.py) |
| **2 — bridge** | specific realizations: an actual ReLU network and an input law, read into the abstract objects by an explicit, versioned dictionary | targeted experiments that test the theory's dichotomies on the realization; a record of what each dictionary can and cannot see | [mlp-bridge.md](mlp-bridge.md), [experiments/mlp_bridge.py](experiments/mlp_bridge.py) |
| **3 — unlocks** | the mechanism by which local conditions control global quantities with explicit bounds, in frontier work (high-dimensional expanders, high-order random walks, log-concave sampling, agreement testing, tiling algebras and gap labelling) | the essence of each reduction, and what would have to be true of the abstract objects for the same reduction to hold there, stated as conjectures | [local-to-global-unlocks.md](local-to-global-unlocks.md) |

---

## 1. Rules for prong 1 (theory)

- **P1.1** Definitions are made for the abstract objects only. A definition may be *motivated* by a realization but must be stated so that it applies verbatim to a Bratteli diagram, a substitution tiling or a family of commuting projections. If it cannot be, it is a prong-2 hypothesis about a dictionary.
- **P1.2** Every statement carries a label: Fact (cited, with the source's own statement), Derivation (proved here), Conjecture, Interpretation. Numerical checks accompany derivations where possible, as in `checks/`.
- **P1.3** Prong 1 takes *questions* from prongs 2 and 3, never *answers*. A measured number on a network is a reason to ask whether a construction exists, not evidence that it does.
- **P1.4** The standing modelling constraints hold: objects are abstract barycentres, not activation vectors, centroids, activation-weighted states or coordinates; the arrow includes the ReLU *before* the linear map; projection information is the included-vertex pattern and the face dimension, which may jump by more than one; no carrier space, metric or walk is assumed before the action of arrows on restrictions is understood; weights are not elements of the arrow algebra.

## 2. Rules for prong 2 (bridge)

- **P2.1 Negative-result rule.** A result that contradicts an expectation is charged first to the dictionary. The write-up must name, in advance, the naive points of the dictionary in use, and after the experiment say which of them could have produced the result and what the next dictionary version changes.
- **P2.2** Experiments test *dichotomies* stated by prong 1 (sufficient or not; coboundary or not; groupoid or frame composition; trivial or non-trivial nerve), never "does the theory work". Each experiment names the statement it targets.
- **P2.3** The dictionary is explicit and versioned. Every object of the theory is given its reading in a table with a status column; anything the theory leaves open (Note 1 §2.1) is marked as a choice.
- **P2.4** Every estimator is calibrated: a surrogate with known answer (e.g. a Markov resampling of the histories) is run alongside so that estimator artefacts are not reported as findings.
- **P2.5** Stand-ins and real networks are distinguished. Stand-ins (random or numpy-trained) establish that the pipeline measures what it claims; conclusions about trained networks of interest wait for their weights. The script accepts weight files for that reason.

## 3. Rules for prong 3 (unlocks)

- **P3.1** Each unlock is written as *global object / local certificate / explicit bound / cost*, from the sources' own theorem statements, with the statement quoted or cited by number.
- **P3.2** A transfer to the framework is a Conjecture with an explicit "what must be true". It is stated for the abstract objects and never optimised for a particular realization; the test is whether it would be stated identically for a tiling or a Bratteli diagram.
- **P3.3** Sources are read, not remembered. A citation from memory is marked as such until the source has been opened; a wrong identifier is corrected in the note.
- **P3.4** Guards are recorded: what in the source theory does not transfer (purity, reversibility, two-sidedness, bounded degree) and why.

## 4. Communication between prongs

- **C1** Findings travel as **messages**: each note ends with a "Messages to the other prongs" section, and begins by stating which inputs from other prongs it used. A finding never travels as a redefinition.
- **C2** The drag test. Before anything from prong 2 or 3 is adopted by prong 1, ask: would it be stated the same way for a tiling, a Bratteli diagram, a family of projections? If not, it stays where it came from.
- **C3** Convergence is checked, not assumed. When two prongs reach the same statement (e.g. "the obstruction is a coboundary", reached by Note 1 §6 and by prong 3 §4), it is recorded in both with the route each took, so that agreement is informative rather than circular.
- **C4** Disagreement is not resolved by majority. If prong 2 measures something prong 1 says cannot happen, prong 1 checks its hypotheses (the realization may be outside the class the theorem is about: this is exactly what happened with the non-Markov history law, mlp-bridge.md §3.2) and prong 2 checks its dictionary; neither changes the other's objects.

## 5. Status board (2026-10-01)

| prong | established | open, next |
|---|---|---|
| 1 | the conditional arrow algebra, its three resolutions, conditioning calculus, KMS/Gibbs states, coboundary criterion for sufficiency, propagation numbers; the complex as nerve and coherence complex, faces as non-liftable decompositions, the composition fork, the history complex | definitions requested by prong 3: links of the history complex, the square 2-complex in which $[F]$ lives and its coboundary expansion, stationarity of the Bratteli diagram; questions from prong 2: cocycles with memory, uniqueness of the $G$-weighted composition, the minimal sufficient sub-frame |
| 2 | dictionary v1 and the experiment pipeline, calibrated by a Markov surrogate; on stand-ins: the face process has 25–35 % memory, $F$ is far from a coboundary, arrows are $G$-continuous, nerves are trivial, training is invisible to v1 | dictionary v2: the minimal sufficient sub-frame (D2), weight-only arrows (D1), faces at a tolerance (D3); runs on the actual competition weights |
| 3 | seven reductions written as global / local / bound / cost; the common mechanism (hereditary class, restriction–co-restriction operator, invariant controlled by local data, homogeneity); six candidate unlocks U1–U6 with their "what must be true"; guards | prong 1's three definitions; measurement of link spectra and square sums on realizations under dictionary v2 |
