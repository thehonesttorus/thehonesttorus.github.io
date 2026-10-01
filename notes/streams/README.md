# Parallel streams of the Phase 2 campaign

## Fresh-slate design (from 1 Oct 2026, 17:45 UTC)

On 1 Oct the user ruled that the competition system starts from a fresh slate, driven by the programme's theoretical principles, with no adaptation of existing estimators. The shared contract is [../fresh-slate/BRIEF.md](../fresh-slate/BRIEF.md); the user's input of 1 Oct is [../fresh-slate/input-2026-10-01-matchings-bethe-cluster.md](../fresh-slate/input-2026-10-01-matchings-bethe-cluster.md). Each design stream derives an estimator from one principle, prototypes it, and measures it in two stages (quick small-width runs projected to n = 1024, then the real harness).

| stream | where | driving principle |
|---|---|---|
| `../fresh-slate/designs/faces/` | [session](https://claude.ai/code/session_01C8zycHwyntKNFLnm15bGU6) | faces, Gaussian barycentres and conditional expectations, glued local-to-global |
| `../fresh-slate/designs/bethe/` | [session](https://claude.ai/code/session_01KhWeB4VLg4VycL145eaA45) | tree-based local-to-global on pseudorandom geometry (Bethe, cavity, loop series, covers) |
| `../fresh-slate/designs/signings/` (finished: raw 1.04e-6 at 1024, adjusted ≈ 2.3e-7) | [session](https://claude.ai/code/session_01Q5meUvTjCCiZqgC9ckDpmJ) | matchings, hafnians and signings (determinants, Pfaffians, gauge averages) |
| `../fresh-slate/designs/tropical/` (finished: not competitive) | [session](https://claude.ai/code/session_01Nxk2j9cDDW1yXsPEw5ZzdJ) | tropical skeleton of the network and its lift to finite temperature |
| `../fresh-slate/designs/heisenberg/` | [session](https://claude.ai/code/session_013o8FK6ftrxpJEts6ucFHWM) | observables pulled back through the arrows; Dirichlet forms on Gaussian space |
| `../fresh-slate/designs/markov/` (finished: raw 4.18e-7 at 1024, adjusted 2.2e-7) | [session](https://claude.ai/code/session_0161pRau9Ub5kdDyNkrkxhcp) | Markov networks, conditional mutual information, exchange relations, recovery maps |
| `../fresh-slate/breakthrough/costate/` | [session](https://claude.ai/code/session_013EjyY2ofzhY3rEZoAzxCrd) (from 19:09) | breakthrough round: the minimal co-state — the smallest observable algebra the readout needs, closed under pull-back |
| `../fresh-slate/breakthrough/interpolation/` | [session](https://claude.ai/code/session_01JXmPTixz71xn11Hzx62BYj) (from 19:09) | breakthrough round: interpolation in the disorder — smart paths, Bolthausen conditioning (Onsager memory), covers |
| `../fresh-slate/breakthrough/region/` | [session](https://claude.ai/code/session_01TjtRHMJmmfJ9TaXzrdKnH8) (from 19:09) | breakthrough round: the feasible region — what any ≤ 0.15 B, raw ≤ 1.5e-8 estimator must carry, and a design inside it |
| `../fresh-slate/bench/` | [session](https://claude.ai/code/session_019vnQdfQcg2TrCQcyC1WYW2) (was est-accuracy) | design-agnostic benchmark sets, truth, the Stage Q evaluator, the Stage P harness |
| `../fresh-slate/scaffold/` | [session](https://claude.ai/code/session_013Q5iH6p39ofihtZhCPeTka) (was est-cost) | grader-safe flopscope template, porting kit with measured primitive costs, residual budget |
| `../fresh-slate/CONVERGENCE.md` | coordinating session (19:15) | what the six designs converged on after round one, the gap (30–60×), and the open question for round two |
| `../fresh-slate/foundations*.md` | background workflow, coordinating session | first-principles problem structure, the programme's unlocks, verification of the user input, neutral measured facts |

## Earlier streams (moment-chain era, wrapping up as fact-finding)

The measurements below are kept as facts about the object (summarised neutrally in BRIEF §3 and in foundations-facts.md); their algorithms are not a starting point for the fresh-slate designs. The submission stream continues only as insurance.

| stream | where it runs | question | deliverable |
|---|---|---|---|
| `submission/` | [sibling session](https://claude.ai/code/session_01MJTw1ZuvoKiLkEWpcU1skv) | Is the public V29 / V25 chain a safe fallback submission (FLOPs, residual-time margin against the 0.4 s cap, 120 s wall, 5 s setup, 8 GB, non-suite smoke shapes, packaging)? | validated bundles, go/no-go table, upload instructions |
| `coef-ensemble/` | [sibling session](https://claude.ai/code/session_01VYMnZMdNTUc4ZQEZJ9yZoj) | Coefficient drift vs width (answered at 1024: no drift), then, from 17:00, a cheap carrier for the (2,1,1) fourth-cumulant slice. | coefficient tables, held-out D21 errors, width extrapolation |
| `chain128/` | [sibling session](https://claude.ai/code/session_016cmvEqAAVYWsNUDTZ4wopY) | Does a better D21 interface lower the final-layer MSE of a full chain? Dense K=3 chain at width 128 with closure variants and teacher forcing. | **finished**: REPORT.md (closure ladder maps one-to-one onto final MSE; MSE = a + k·ε²; the (2,1,1) κ4 slice is the lever) |
| `old-content/` (finished) | [sibling session](https://claude.ai/code/session_01PyRVXWTvmQeHszAB16oDry) | What cheap representation carries the transported old-source third-cumulant content? **Verdict: none; every carrier's size grows with width, and absorbing it into renormalised births fails (10–35× the target).** | per-source tracker, carrier × (D21 error, cost) table |
| `oracle1024/` (finished: final noise-free ladder at width 1024) | [sibling session](https://claude.ai/code/session_012m7DoPtgXQv4s6iAqRWKak) | The closure ladder at the real shape (width 1024) by streaming Monte Carlo, without n³ tensors. | stream_oracle.py, the ladder at width 1024 |
| `est-accuracy/` | [sibling session](https://claude.ai/code/session_019vnQdfQcg2TrCQcyC1WYW2) (from 17:04) | Lower the public chain's raw MSE at width 1024 through its fourth-cumulant handling and births. | lever × (raw, C/B, residual, adjusted) table, bundle |
| `est-cost/` | [sibling session](https://claude.ai/code/session_013Q5iH6p39ofihtZhCPeTka) (from 17:04) | Cut the public chain's FLOP bill and flopscope call count at unchanged numbers (residual ≤ 0.2 s). | ledger before/after, bundle |
| `theory/` | background workflow, coordinating session | Which second-order diagrams explain the depth drift of the fitted coefficients? (each claim adversarially verified) | closure2.py, derivations, toy validation |
| `costmodel/` | background workflow, coordinating session | Metered flopscope cost of every operation a closure-based chain needs, and a design calculator. | cost.py, op and design tables |

Also from the coordinating session's support workflow: `../digests/phase2-intel-2026-10-01.md` (new public information since 10 Sep) and `../../infra/azure/CHECKS.md` (offline validation and fixes of the Azure grid scripts).

## Round 3: the local-to-global essence programme (from 1 Oct 20:45 UTC)

Brief: `notes/essence/BRIEF.md`; the user's mandate verbatim: `notes/essence/INPUT-2026-10-01-local-to-global.md`.

| team | directory | session | scope |
|---|---|---|---|
| A | `notes/essence/A-permanents-polynomials/` | [session](https://claude.ai/code/session_01FVAU2TZ37gB5fKRxEXB4S2) | permanents, hafnians, Barvinok interpolation, Bethe permanent, Clifford estimators, Lorentzian polynomials |
| B | `notes/essence/B-pseudorandomness/` | [session](https://claude.ai/code/session_01EAVH9mdGoca8mpKptzZc4u) | expander-walk pseudorandomness, weighted PRGs and Richardson precision amplification, PRGs for polytopes |
| C | `notes/essence/C-free-probability-criticality/` | [session](https://claude.ai/code/session_01SnorBLeaw2toS2DkuZcanK) | free probability, two-projection geometry, the wall at 2 and He-init criticality |
| D | `notes/essence/D-elimination-sparsification/` | [session](https://claude.ai/code/session_01J7YMAiCbTc5HPjyS9a5L2k) | approximate Gaussian elimination, sparsifiers, Kadison–Singer, operator scaling, intrinsic freeness |
| E | `notes/essence/E-tilings-cones/` | [session](https://claude.ai/code/session_01JaEDBoti5YcGCy8LAXqmhf) | tilings, gap labelling, conic intrinsic volumes, random conical tessellations |
| F | `notes/essence/F-nc-local-to-global/` | [session](https://claude.ai/code/session_01Jg65pQarkMSSY1zxCFTtwM) | noncommutative trickle-down, localization schemes, commuting squares; cross-team synthesis |

Experiments running alongside: costate (dilation-sector share of the memory at n = 1024), region (CP-merged old tier; FC under the wall-feasible Strassen mix), bethe (localization estimator conditioned on the collective coordinate). Faces and heisenberg are writing final verdicts and stopping.

### Round 3 status (22:30 UTC)

- **Memory carriers that work** (inside FC at n = 1024; notes/fresh-slate/CONVERGENCE.md): team D's causal frozen frames at k = 2n/age (lossless on 6/6, ≈ 0.53 B wall-feasible; constant ≈ 9 products per unit of resolution); team F's cohort basis for ages ≥ 8 (≈ 1–2 u per layer, 1.01–1.04×); team B's constant-in-a channel (80 % of old energy at ≈ 0 cost). The bill is now ages 1–7.
- **Theory:** team G's Theorem G10 (no finitely summable grading of neuron space is compatible with ≥ 2 fresh layers; compression only along the age axis); team F's G5 (Dixmier-critical trickle-down: the polynomial-loss exponent is a Dixmier trace); team B's G1 (gapless torus: the λ-term rank is Dixmier-critical exactly at Pólya's d = 2; free products are universally "d_eff = 2"); team E's G3 (the carrier holds at depth 32 and its cost grows as L log L; the memory's tiling class is Fibonacci/Penrose).
- **Coordinator experiments** (`notes/essence/H-rational-memory/`): response-matched merged histories reproduce the old memory's future at R = n/4 with ≈ 1 % residual energy (n = 256), where pruning leaves 79 % and tensor-norm CP fails; the fit cost is the open problem.
- **Bench complete** (bench session): truth for w64 … w1024 (w1024_d16: 6 MLPs, N = 2e6, noise ≤ 4.4e-8) and the smoke shape w256_d32.
- **Background workflows** (judge panel, foundations check, bridges-synthesis critic) stopped at ≈ 20:10 UTC; finished parts are in the repo.
