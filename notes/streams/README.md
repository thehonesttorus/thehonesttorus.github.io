# Parallel streams of the Phase 2 campaign (started 1 Oct 2026, 15:13 UTC)

Each stream runs in its own container (a sibling cloud session) or as a background workflow in the coordinating session, writes only under its own directory here, and keeps a `REPORT.md` current. The coordinating session integrates finished results into [../competition-plan.md](../competition-plan.md).

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
