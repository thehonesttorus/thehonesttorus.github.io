# Competition plan: sixteen days, three lines of attack, one grid

*Written 2026-10-01 against [competition-phase2.md](competition-phase2.md) (rules), [competition-landscape.md](competition-landscape.md) (public frontier) and [infra/azure](../infra/azure/README.md) (the grid). Submissions close 17 Oct 23:59 UTC; team freeze 2 Oct 23:59 UTC; write-up 24 Oct. The plan is revised as results come in; the revision log is at the end.*

---

## 0. The target and the arithmetic

| what | number |
|---|---|
| money line (public board, 1 Oct) | adjusted ≤ 1.6–1.8e-9 |
| equivalent raw at the floor (C/B ≤ 0.10) | ≤ 1.6e-8 |
| equivalent raw at C/B = 0.15 | ≤ 1.1e-8 |
| best public method (504aldo V29) | raw 2.13e-8 at 0.253 B → 5.4e-9 |
| same chain without the old-source tier | 0.150 B but raw 1.3–2.9e-7 |
| reference augmented K=3, un-metered | raw 1.80e-8 |
| ground-truth floor | 7.5e-11 |

So the problem is: deliver the accuracy of an augmented K=3 chain (≈1.5e-8) **or better**, at **≤ 0.12 B**, with zero failures on 100 fresh MLPs. Three lines, in order of expected value per day.

---

## 1. Line A: reproduce, then compress (days 1–4)

1. Run 504aldo's V25 and V29 through our runner on `mini` (100) and `full` (1,000) and on fresh seeds; record raw, C/B, residual time, memory. Expected: raw 2.1–2.3e-8, C/B 0.253/0.367, residual 0.17–0.34 s. This calibrates the grid against the public board (504aldo reports local 8-dump numbers 6.5 % above the board).
2. Re-run the oracle tests on **1,000 networks instead of 8** using keenanpepper's N = 1e9 joint moments (true D21, D3, κ₄ diagonal at every layer): the D21 ε² law, the value of each old source by age, the λ law for κ₄. Cheap on the grid, and it tells us where the 1.5–2× to the leaders can come from with a thousand-network statistical power nobody has used yet.
3. Try the compressions the author did **not** try, each measured against the ε ≤ 2 % rule: (a) a shared basis learned offline (fixed, data-independent subspaces for old sources at each age, fitted on the 14k corpus, so that no range finder is paid online); (b) a **learned memory** for the old tier: predict the old sources' D21 contribution from the young tier's state and the Ω-sketches (the author's +0.02 R² was on 8 networks with cheap features; the sketch targets exist for thousands); (c) a response closure for κ₄ beyond one scalar (SSC-style: a few covariance-response modes with a fitted coupling), evaluated on the same 1,000-network oracle surface.

Line A's deliverable by day 4 is a submission at or below 504aldo's score with a measured margin on fresh seeds, plus the oracle tables.

## 2. Line B: trajectory-fitted corrections at scale (days 3–10)

The census's one positive learned result is a per-layer corrector fitted DAgger-style on the chain's own rolled-forward trajectories; yujigon showed the 64 constants of the cheap family absorb accumulated drift. Nobody has fitted such a corrector on top of a K=3 chain with thousands of networks. On the grid:

1. Roll the chain forward on the d8b corpus (14,048 networks with N = 1e8 marginal targets; 2,048 with the teacher bake), recording the live state at every layer.
2. Fit per-layer corrections sequentially (layer ℓ fitted on states already corrected at layers < ℓ), with the state features that are cheap online (μ/σ, variances, D3, the κ₄ diagonal, the λ observable, the shared-basis energies), frozen before any evaluation; validate on `mini`/`full` and on fresh seeds; report the public-split number last.
3. Keep every fitted table small and shape-guarded (the smoke-test trap), and keep the derivation of what each table absorbs.

Expected value: the Phase 1 analogue gave 8.4× over plain kprop; on a chain already at its truncation floor the honest expectation is 1.2–2×, which is the whole gap.

## 3. Line C: the structural line (days 1–16, in the background)

This is the prong-2 realization of the programme in [research-program.md](research-program.md), kept as a separate line so it cannot be dragged by Lines A and B and they cannot be dragged by it. The public forensics say the Phase 1 front was "deterministic, non-moment, angular"; the Euler–Stein facet identity writes E[h] as a sum over activation cones; our framework's objects are exactly those cones (faces) and the arrows between them. Concretely:

1. Measure on the Phase 2 shape what the earlier notes measured on toy networks: the face process (activation patterns per layer) on He-init 1024×16 networks under the Gaussian state, its memory, the coboundary residual of the derived kernel, the overlap-continuity of arrows. The Azure atlas jobs can record pattern statistics alongside moments at no extra cost.
2. Ask the question the field has not: is the final-layer mean better represented by a **gate decomposition** E[h_L] = Σ over gate configurations of the last two layers than by moments? Gate probabilities and gate–activation cross moments (`gate_GG`, `gate_GX` in the community bakes) are the natural state; a closure in that state is a different projection, and the census's open problem O3 (the joint law of near-threshold neurons) is exactly where the moment chains lose accuracy.
3. Anything from Line C that produces a number competitive with Line A is promoted into the submission path; anything that does not stays as a write-up contribution (the algorithmic-contribution prize is discretionary and judged on the write-up).

---

## 4. The grid, as used by the three lines

| job type | pool | what | scale |
|---|---|---|---|
| estimator evaluation on parquet shards | CPU pools, 4 vCPU per task | `whest run` with graded caps, wall cap relaxed | 1,000 MLPs in ≈15 min at 100 tasks |
| chain roll-outs with state recording | CPU pools | the chain in plain numpy (not metered), dumping per-layer state | 14k networks |
| fresh-seed bakes | GPU pool | `whest dataset bake --torch`, N = 1e8 (floor 7.5e-10) or 1e9 | 1,000 networks ≈ 125 GPU-hours at 1e8 |
| moment atlases | GPU pool | marginal moments to order 6, gates, dense pair blocks | as needed |
| fits | one large CPU node or GPU | ridge / small nets on recorded trajectories | hours |

Staging first: the public dataset (77 GB), keenanpepper's full-split joint moments (337 GB), the d8b corpus and teacher bake, the Ω-sketch bakes. All Hugging Face → blob, once, in-region.

---

## 5. Submission discipline

- Every candidate runs on `mini`, on `full`, and on ≥ 200 fresh-seed MLPs with N ≥ 1e8 truth before a submission slot is spent; 10 slots per UTC day, failures count.
- Every candidate passes: `whest validate`, the subprocess runner at the graded caps with `--max-threads 1`, a depth-32 and a width-4/depth-2 shape probe, a cold-start `setup()` timing under 2.5 s, peak memory under 6 GB, residual under 0.25 s on every MLP of `mini`.
- Randomness only from `mlp.seed` / `ctx.seed`; no per-MLP adaptive constants fitted on in-sample statistics; every shipped table documented with what it absorbs.
- The designated final submission is chosen on fresh-seed MSE, not on the public board.

---

## 6. Timeline

| days | milestone |
|---|---|
| 1–2 (Oct 1–2) | grid live (needs network access or a machine with `az`); datasets staged; 504aldo V25/V29 reproduced locally and on `mini`; **team freeze by 2 Oct 23:59 UTC** |
| 3–4 | 1,000-network oracle tables; first submission (reproduction with margin) |
| 5–8 | Line A compressions and Line B fits in parallel on the grid; daily submissions of the best fresh-seed candidate |
| 9–12 | integrate; robustness campaign; Line C measurements written up |
| 13–16 (Oct 13–16) | freeze, designate, keep one fallback submission at a known score |
| 17 Oct | submissions close 23:59 UTC |
| by 24 Oct | algorithmic-contribution write-up |

---

## Revision log

- 2026-10-01: first version, written before any grid job has run; numbers are the field's, not ours.
