# Parallel streams of the Phase 2 campaign

## Fresh-slate design (from 1 Oct 2026, 17:45 UTC)

On 1 Oct the user ruled that the competition system starts from a fresh slate, driven by the programme's theoretical principles, with no adaptation of existing estimators. The shared contract is [../fresh-slate/BRIEF.md](../fresh-slate/BRIEF.md); the user's input of 1 Oct is [../fresh-slate/input-2026-10-01-matchings-bethe-cluster.md](../fresh-slate/input-2026-10-01-matchings-bethe-cluster.md). Each design stream derives an estimator from one principle, prototypes it, and measures it in two stages (quick small-width runs projected to n = 1024, then the real harness).

| stream | where | driving principle |
|---|---|---|
| `../fresh-slate/designs/faces/` | [session](https://claude.ai/code/session_01C8zycHwyntKNFLnm15bGU6) | faces, Gaussian barycentres and conditional expectations, glued local-to-global |
| `../fresh-slate/designs/bethe/` | [session](https://claude.ai/code/session_01KhWeB4VLg4VycL145eaA45) | tree-based local-to-global on pseudorandom geometry (Bethe, cavity, loop series, covers) |
| `../fresh-slate/designs/signings/` | [session](https://claude.ai/code/session_01Q5meUvTjCCiZqgC9ckDpmJ) | matchings, hafnians and signings (determinants, Pfaffians, gauge averages) |
| `../fresh-slate/designs/tropical/` | [session](https://claude.ai/code/session_01Nxk2j9cDDW1yXsPEw5ZzdJ) | tropical skeleton of the network and its lift to finite temperature |
| `../fresh-slate/designs/heisenberg/` | [session](https://claude.ai/code/session_013o8FK6ftrxpJEts6ucFHWM) | observables pulled back through the arrows; Dirichlet forms on Gaussian space |
| `../fresh-slate/designs/markov/` | [session](https://claude.ai/code/session_0161pRau9Ub5kdDyNkrkxhcp) | Markov networks, conditional mutual information, exchange relations, recovery maps |
| `../fresh-slate/bench/` | [session](https://claude.ai/code/session_019vnQdfQcg2TrCQcyC1WYW2) (was est-accuracy) | design-agnostic benchmark sets, truth, the Stage Q evaluator, the Stage P harness |
| `../fresh-slate/scaffold/` | [session](https://claude.ai/code/session_013Q5iH6p39ofihtZhCPeTka) (was est-cost) | grader-safe flopscope template, porting kit with measured primitive costs, residual budget |
| `../fresh-slate/foundations*.md` | background workflow, coordinating session | first-principles problem structure, the programme's unlocks, verification of the user input, neutral measured facts |

## Earlier streams (moment-chain era, wrapping up as fact-finding)

The measurements below are kept as facts about the object (summarised neutrally in BRIEF §3 and in foundations-facts.md); their algorithms are not a starting point for the fresh-slate designs. The submission stream continues only as insurance.

| stream | where it runs | question | deliverable |
|---|---|---|---|
| `submission/` | [sibling session](https://claude.ai/code/session_01MJTw1ZuvoKiLkEWpcU1skv) | Is the public V29 / V25 chain a safe fallback submission (FLOPs, residual-time margin against the 0.4 s cap, 120 s wall, 5 s setup, 8 GB, non-suite smoke shapes, packaging)? | validated bundles, go/no-go table, upload instructions |
| `coef-ensemble/` | [sibling session](https://claude.ai/code/session_01VYMnZMdNTUc4ZQEZJ9yZoj) | Coefficient drift vs width (answered at 1024: no drift), then, from 17:00, a cheap carrier for the (2,1,1) fourth-cumulant slice. | coefficient tables, held-out D21 errors, width extrapolation |
| `chain128/` | [sibling session](https://claude.ai/code/session_016cmvEqAAVYWsNUDTZ4wopY) | Does a better D21 interface lower the final-layer MSE of a full chain? Dense K=3 chain at width 128 with closure variants and teacher forcing. | chain.py, MSE per variant, measured error law |
| `old-content/` (finished) | [sibling session](https://claude.ai/code/session_01PyRVXWTvmQeHszAB16oDry) | What cheap representation carries the transported old-source third-cumulant content? **Verdict: none; every carrier's size grows with width, and absorbing it into renormalised births fails (10–35× the target).** | per-source tracker, carrier × (D21 error, cost) table |
| `oracle1024/` | [sibling session](https://claude.ai/code/session_012m7DoPtgXQv4s6iAqRWKak) | The closure ladder at the real shape (width 1024) by streaming Monte Carlo, without n³ tensors. | stream_oracle.py, the ladder at width 1024 |
| `est-accuracy/` | [sibling session](https://claude.ai/code/session_019vnQdfQcg2TrCQcyC1WYW2) (from 17:04) | Lower the public chain's raw MSE at width 1024 through its fourth-cumulant handling and births. | lever × (raw, C/B, residual, adjusted) table, bundle |
| `est-cost/` | [sibling session](https://claude.ai/code/session_013Q5iH6p39ofihtZhCPeTka) (from 17:04) | Cut the public chain's FLOP bill and flopscope call count at unchanged numbers (residual ≤ 0.2 s). | ledger before/after, bundle |
| `theory/` | background workflow, coordinating session | Which second-order diagrams explain the depth drift of the fitted coefficients? (each claim adversarially verified) | closure2.py, derivations, toy validation |
| `costmodel/` | background workflow, coordinating session | Metered flopscope cost of every operation a closure-based chain needs, and a design calculator. | cost.py, op and design tables |

Also from the coordinating session's support workflow: `../digests/phase2-intel-2026-10-01.md` (new public information since 10 Sep) and `../../infra/azure/CHECKS.md` (offline validation and fixes of the Azure grid scripts).
