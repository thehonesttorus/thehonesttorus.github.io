# nce-arc3: the Noncommutative Correspondence Engine for ARC-AGI-3

An ARC-AGI-3 agent derived from the operator-algebraic account of in-context learning in
*Modular Correspondences, Subfactor Indices, and Two-Projection Geometry* (the monograph in
`docs/`-adjacent notes is referred to below by its theorem numbers). It is not a neural network,
not a reinforcement learner, not a search over the true environment, and it contains no game
knowledge: every rule below is traced to a theorem of the monograph, to one of the referee's
amendments, or to an interface fact of the competition harness (palette `0..15`, lattice
`64×64`, action ids `0..7`, the state flag, `levels_completed`, `win_levels`).

The design document is `docs/architecture-spec.md`; the exact promotion inequality with the
referee's correction E6 is derived in `docs/promotion-inequality.md`. This README is the map
from the code to the theory and the record of what the system does and does not yet do.

## The idea in one paragraph

A frame is a finite measure on the lattice with values in the palette; its perception is a
Harder–Narasimhan flag of colour strata (Tyler fixed point of the one-hot measure) and, inside
each stratum, the kernel of the translation-generated Dirichlet form (connected components) with
D4 ⋉ Z² canonical forms. A transition under an action is a *correspondence* between two such
base algebras, and the agent's whole job is to find the shortest exact word for each action in a
typed library of correspondence families (a Cuntz–Pimsner algebra over a finite graph), ordered
by the Parry measure of the library's critical KMS state at β_c. Words are verified exactly on
every stored demonstration, never fitted approximately; hypotheses that survive form a coherent
set; actions are chosen by the Kato two-projection rule (the action whose predicted outcomes are
least pure across the coherent set, priced by the replay cost of death); goals are second-quantized
predicates on the predicted post-state learned from level completions and split by planning
inside the verified model; and the library grows at level ends by the exact code-length test
(adjoining a composite raises β_c, so it is admitted only when the corpus of verified words gets
shorter).

## Code map, with status tags

Tags: **[T]** forced by a theorem the brief marks sound; **[T\*]** by transfer with the
amendment's citation; **[D]** a convention needed to make a [T] statement operational, with the
model stated in the docstring; **[I]** an interface fact.

| Module | What it is | Theory |
|---|---|---|
| `nce/tyler.py` | Moment map, Tyler fixed point, whitening, HN filtration of a measure by slope gaps | T8 (Tyler/HN perception) [T] |
| `nce/algebra.py` | Base algebra of a frame: colour strata, components as the Dirichlet-form kernel, D4⋉Z² canonical forms, frame vectors, split-aware constructor | T8, T13 [T]; 4-connectivity from the generator set of Z² [T]; the top stratum as canvas [D] |
| `nce/match.py` | Decomposition of a transition into gauge moves, recolourings, reshapes, occlusions, creations, annihilations by maximising KMS weight of the explaining word | Prop. 1.3, Thm 1.5(iii), C11 [T]; containment moves and in-place congruent pairs [D] |
| `nce/library.py`, `nce/kms.py` | The typed library (vertices F, C, P, TOP; letters as families with summand-count dimensions), spectral radius, Perron vectors, β_c by bisection, the conditional Parry kernel, exact promotion decision | T5, T11, Thm 2.1, E3, E6 [T] |
| `nce/words.py` | Clauses (selector · gate · body · fill · term), clause sets, creations, death gates, level gates, exact prediction and rendering, word contexts | T4 (semiring compositionality), T6 (exact verification), T9 (level gates) [T]; rendering in HN order and the `Fill`/`REVEAL` member [D] |
| `nce/synth.py` | Parry-ordered synthesis: behaviour classes from the decomposition, selector/gate fitting, admissible-fill intersection, readings in Parry order (identity before Rot), LevelGate retry, death-gate fitting, coherent set of exact survivors | T5, T6, T7, T12 [T]; sub-component (cell) classes [D] |
| `nce/store.py` | Demonstrations, the ignorance table (non-determinism detector), chronological canvas memory | T1, T12 [T]; canvas memory under occlusion [D] |
| `nce/split.py` | Downward refinement of strata under residuals: stable-remainder split rules | Section 2.5 of the spec, two-projection graft [D] |
| `nce/kato.py` | Candidates, purity deficit D = ½(1 − Tr ρ²), the newborn rule relative to each word's own context, death pricing, tie-breaking by the ignorance table | T2, T3 (with E10), T7, C8 [T] |
| `nce/goals.py` | Goal words as second-quantized predicates, negatives from every non-winning post-state, positives from completions, alignment and splitting | T7, Component 4 [T]; the pre-completion prior restricted to false-on-observed [D] |
| `nce/plan.py` | Uniform-cost search inside the verified model, death-aware | Section 7 of the spec; N16 [T] |
| `nce/maternity.py` | Library growth at level ends: composites with realised-summand dimensions, exact E6 test, abbreviation of weighed words | T5, T10, E6, Thm 2.1 (Wold newborn compression) [T] |
| `nce/agent.py` | The in-episode loop: try once, goal-splitting, Kato, fallthrough; observe, verify, re-synthesise, promote | Sections 5–6, 10 of the spec |
| `nce/proxy.py`, `nce/engine.py`, `kaggle/my_agent.py` | Process-per-game proxy with in-process fallback, the Kaggle `Agent` adapter | [I] |
| `tools/harness.py` | Fast engine-level blind scoreboard: levels, action counts, RHAE, count-only diagnostics | [I] |
| `tools/play_local_framework.py` | The deployment-faithful loop through the vendored official framework | [I] |
| `tools/build_kaggle_notebook.py` | Bundles `nce/` into `notebooks/submission.ipynb` | [I] |

### Letter families (dimensions are summand counts, Thm 1.5(iii))

`SelPointed` 2 (component, cell) · `SelAll` 1 · `SelColour` 16 · `SelRank` 16 · `SelType` |Ξ₀| ·
`Step` 4 · `Move` 16128 · `Rot` 7 · `Recolour` 15 · `Fill` 17 (palette + the remembered canvas) ·
`Layer` 2 · `Grow` 64 · `Erode` 4 · `Copy` 2 · `Annihilate` 16 · `Create` |Ξ₀|·8·4096 ·
`Bind` 6(|Ξ₀|+1)+2(16+|Ξ₀|+8) (pair predicates, unary tests on the component or the pointed one, couplings) ·
`Meter` 2·16·4097 (a colour mass against a threshold) · `Gate` 2 · `Const` 2 · `LevelGate` win_levels · `Over` 1.
The d = 1 subgraph is acyclic, so β_c exists (Thm 2.1(iv)); this is asserted at every refresh.

Two further structures complete a word:

- **The memorised part (T12).** A clause set may carry exact (pre-state, datum, post-state) memories for the
  demonstrations its generalising part does not reproduce. It is admitted only when its letters plus
  β_c·log(16·4096) per memorised cell are a shorter code than memorising every demonstration.
- **Coupled clauses (pair vertex P).** A clause gated by `cofire` fires exactly when the clause fitted to
  another behaviour class changed something in the same step (two correspondences tensored into one event).

## The governing rule

No frame of any game was inspected while building this, no game source was read, no mechanic
was catalogued. Public games are used only as a blind scoreboard (levels completed, action
counts) and as a source of agent-internal *counts* (how many hypotheses, how many unexplained
classes, histograms of changed-cell counts). The negative checklist of priors found in prior
agents (N1–N17 in `docs/architecture-spec.md`, Section 0) is enforced: no pixel vectors, no
frame-hash states beyond the ignorance table and planner visited set, no novelty bonuses, no random
exploration, no language models, no game DSL, no magic constants tied to baselines, no discarded
animation frames, no per-level amnesia (earlier levels stay as data behind a LevelGate), no
true-environment search.

## Running

```
cd nce-arc3
python3 -m pytest -q tests/                      # core math, algebra, matcher, proxy, toy mechanics
python3 tools/harness.py nce.agent:NCEAgent --games ls20,vc33 --max-actions 300 --procs 4
python3 tools/play_local_framework.py --games vc33 --max-actions 200     # official framework loop
python3 tools/build_kaggle_notebook.py           # -> notebooks/submission.ipynb
```

The toy tests in `tests/test_mechanics.py` are not games: each isolates one correspondence family
(a gauge move with a lethal wall and its death gate, a cell-level recolour under the pointed
selector, an in-place D4 change found on the second reading, boundary growth, an
action-independent countdown bar that becomes an `Erode` clause in every word, a level-gated
rule, and an exact promotion).

## Status

See `docs/status.md` for the current blind scoreboard and the count-only diagnostics that
drove each change. The honest summary at the time of writing: the engine forms exact words and
plans on every toy mechanic; on the public games the words do not yet close, and the diagnostics
point at two remaining structural gaps that the theory names but the code does not yet implement:

1. **Atoms coarser than a cell and finer than a component** (lattice tilings of a region by an
   observed atom shape, Section 2.5): the pointed selector has only the component and the cell.
2. **Factoring the action-independent flow** (the modular flow σ_t): components that change
   identically under every action are handled today only when the change is an `Erode`/`Grow`
   boundary step; a general autonomous clock word is not yet learned.

Both are extensions of existing families, not new priors.

## License

MIT, see `LICENSE`.
