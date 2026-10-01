# ARC White-Box Estimation Challenge 2026, Phase 2: the rules, the environment, the score

*A complete study of the round we are competing in, assembled on 2026-10-01 from the Phase 2 announcement, the starter kit (whest-starterkit at commit of 2026-10-01, whestbench 0.16.1, flopscope 0.12.1), the dataset card, the forum, and the public leaderboard. The official Rules page on aicrowd.com is rendered client-side and could not be fetched from this container; where the forum quotes it (sections 5.2, 5.6, 12, the single-submission clause) that is noted. Everything else below is quoted from the authoritative sources named in each section.*

---

## 0. Deadlines first

| Milestone | Date and time (UTC) |
|---|---|
| **Registration and team freeze** | **2 October 2026, 23:59** |
| Phase 2 submissions close | 17 October 2026, 23:59 |
| Algorithmic-contribution write-up deadline | 24 October 2026, 23:59 |

Source: Phase 2 announcement (forum, 23 Aug 2026). The older timeline table on the challenge overview page (Phase 2 "August 1 – September 19", results "October 1") is superseded. The team-freeze deadline is tomorrow relative to the date of this note: whoever will be on the team must be registered and in the team by then.

Phase 2 carries **$100,000** in prizes (the overview page lists first, second and third leaderboard prizes; the exact split is on the Rules page), plus a discretionary **algorithmic-contribution prize** judged from a technical write-up (write-up deadline above), plus **community-contribution prizes** of $500–$5,000 for useful feedback, bug and mispricing reports, and forum contributions.

---

## 1. The task

Given the weights of a randomly initialised ReLU MLP, return the expected post-ReLU activation of every neuron at every layer under standard normal inputs, within a fixed analytical FLOP budget.

| item | Phase 2 value | source |
|---|---|---|
| width $n$ | 1024 | announcement, rounds.md |
| depth $d$ (number of weight matrices, each followed by ReLU) | 16 | announcement |
| parameters | 16.8 M | announcement |
| forward pass | $h_\ell = \max(0,\ h_{\ell-1} W_\ell)$, $h_0 = x$, **no biases, no linear output layer** | dataset card |
| weights | i.i.d. $\mathcal N(0, 2/n)$ (He), float32, shape `[depth, width, width]`, layout `(input, output)` so a row vector multiplies on the left | dataset card, problem-setup.md |
| input | $x \sim \mathcal N(0, I_{1024})$ | dataset card |
| target | $Y_{\ell,j} = \mathbb E_x[h_\ell(x)_j]$ for all $\ell=1..16$, $j=1..1024$; the estimator returns a `(16, 1024)` array | contract |
| ground truth | Monte Carlo with $N=10^9$ inputs per MLP, float64 accumulation, float32 output; baked on GPU (torch) | dataset card |
| ground-truth noise floor | `avg_variance / N` = 0.0748 / 1e9 = **7.5e-11** MSE | ground-truth.md |

Only the **final-layer row** enters the ranking. The other 15 rows are a diagnostic (`all_layers_mse`).

---

## 2. The score

$$
\text{adjusted\_final\_layer\_score}=\frac1M\sum_{m=1}^M \text{MSE}^{\text{final}}_m\cdot\max\!\big(0.1,\ C_m/B\big),\qquad
\text{MSE}^{\text{final}}_m=\frac1n\sum_{i=1}^n\big(\hat Y_{d,i}-Y_{d,i}\big)^2,
$$

with $B = 2^{41} = 2{,}199{,}023{,}255{,}552$ FLOPs and $C_m = F_m$, the FLOPs flopscope metered during `predict()` and **nothing else** (residual wall time is no longer priced; $\lambda = 0$). Lower is better. $M$ = 100 MLPs per graded submission (50 public + 50 holdout, CHANGELOG).

**Consequences that matter for design.**

- **The compute multiplier is floored at 0.1.** Below 10 % budget use ($C_m \le 2.2\times10^{11}$ FLOPs, about 6,550 forward passes), the ranking is a pure MSE ranking. Every bundled example is already at the floor. The forum post of 17 Sep 2026 verified the rule on graded rows: `adjusted = mse × max(0.1, util)`.
- **Above the floor, the trade is linear.** Spending $u > 0.1$ of the budget multiplies MSE by $u$, so extra compute pays only if MSE falls faster than $1/u$. Plain Monte Carlo falls exactly as $1/u$, which is why a plain sampler's adjusted score is flat at about $1.2\times10^{-6}$ from 6,500 passes upward (problem-setup.md). The current leaders sit at utilisation 0.11–0.16, i.e. *just above the floor*: their marginal compute buys accuracy at slightly better than the $1/u$ rate, which is the signature of a method with a sampling component on top of an analytic one.
- **Any failure on any MLP is catastrophic.** A blown cap (FLOP budget, 400 ms residual, 120 s wall), an exception, a wrong shape or a non-finite value zeroes that MLP's prediction **and forces its multiplier to 1.0**. A zero prediction has MSE ≈ 0.9 (the all-zeros baseline scores 0.9095), so one failed MLP in 100 adds ≈ $9\times10^{-3}$ to a mean that the leaders hold at $1.6\times10^{-9}$. One failure in the 100-MLP suite is a seven-orders-of-magnitude loss. **Robustness dominates everything**: every code path must be exercised on the full public split before a submission, and the estimator must degrade gracefully (fall back to a cheaper method) rather than risk a cap.
- **`setup()` is worse than per-MLP failure.** Its 5 s cap is per run, it runs 5–15 times per submission (one per worker process, replaced workers re-run it), and overrunning it once fails the *whole* submission with `SETUP_TIMEOUT`. `setup()` is for loading shipped files, not computing.

**Where the field is (public leaderboard, 1 Oct 2026, 50 rows read).**

| rank | team | adjusted | final-layer MSE | utilisation |
|---|---|---|---|---|
| 1 | J2W | 1.6e-9 | 1.50e-8 | 0.110 |
| 2 | marius_binner | 1.7e-9 | 1.14e-8 | 0.151 |
| 3 | Luna | 1.8e-9 | 1.41e-8 | 0.131 |
| 4 | suliman_tadros | 2.0e-9 | 1.58e-8 | 0.130 |
| 5 | mliston | 2.1e-9 | 1.70e-8 | 0.126 |
| 6–11 | (six teams) | 2.4–2.5e-9 | 1.5–2.0e-8 | 0.13–0.16 |
| ~20 | | 3.1e-9 | 1.7–1.9e-8 | 0.17–0.18 |
| ~48 | | 4.9e-9 | 2.1–2.2e-8 | 0.22–0.24 |

Reference points from the starter kit (mini split, 100 MLPs, N = 1e9 truth): zeros 0.9095; random 0.6856; mean propagation **2.2e-4** (8.66e7 FLOPs); covariance propagation **4.05e-6** (5.17e10 FLOPs, 2.35 % of budget; adjusted 4.1e-7); plain Monte Carlo 1,000 passes 7.8e-5, 10,000 passes 7.6e-6 (15.3 % of budget), adjusted optimum ≈ 1.2e-6 at ~6,500 passes.

So: the leaders' MSE of 1.1–1.5e-8 is **270–370× better than covariance propagation**, **500× better than optimal plain sampling**, and still **150–200× above the ground-truth floor**. The money line is adjusted ≈ 1.6e-9, i.e. MSE ≤ 1.6e-8 at the floor or ≤ 1.1e-8 at 15 % utilisation. The forum's observed "cluster at 2.13e-8" is not the stock covariance baseline (someone tested; 18× off) and is unidentified.

---

## 3. The environment the submission runs in

| resource | limit | source |
|---|---|---|
| execution | isolated, standardised, **CPU-only** container, pinned dependencies, **no network** | overview page |
| grader software | whestbench 0.16.0 + flopscope[server] 0.12.0 (kit pins 0.16.1 / 0.12.1; FLOP counts identical for everything the kit does) | pyproject.toml |
| CPUs | **2 vCPUs pinned to the solution process, 14 to the flopscope backend** | CHANGELOG |
| solution memory | **8 GB** (Phase 1 box was 64 GB) | announcement, CHANGELOG |
| `predict()` wall clock | **120 s per MLP** | announcement |
| residual wall time | **0.4 s per MLP**, hard cap, zeroes the MLP | announcement |
| `setup()` | **5 s** per run, fails the whole submission | CHANGELOG, contract |
| FLOP budget | $2^{41}$ per MLP, enforced analytically; `BudgetExhaustedError` before the op that would cross | contract |
| suite | 100 MLPs (50 public + 50 holdout), one worker serves many MLPs (state is **not** reset between MLPs), several workers per submission | CHANGELOG, contract |
| submissions | **10 per team per UTC day**, shared across members; first week of Phase 2 gave 10 extra for failed ones | announcement |
| seeds | `ctx.seed` (per run) and `mlp.seed` (per MLP) are supplied; randomness must be drawn from them via `fnp.random.default_rng(...)`; seeding from a constant "may cost you prize eligibility" | contract |

**What the code may use (the allowed-code rule, enforced).** Exactly three things: the grader's Python interpreter, the flopscope client API (`import flopscope as flops`, `import flopscope.numpy as fnp`, plus the `whestbench` contract types), and the pure-Python standard library for control flow and bookkeeping. Prohibited, non-exhaustively: vendored numpy/scipy/BLAS; compiled kernels in any form (wheels, `.so`, static binaries, runtime codegen); ctypes/cffi/FFI; asyncio, threads, subprocesses, multiprocessing, `concurrent.futures`; any computation while a flopscope op is in flight (callbacks, hooks, lazy objects); modifying, monkeypatching or intercepting the flopscope client, transport or accounting. **Packing** several independent values into one machine element (bit-packing booleans into integers, etc.) to reduce metered work is prohibited even though every op is metered; lower-precision dtypes and quantisation are not packing by themselves. The residual 400 ms is "plumbing, not a second budget": meaningful arithmetic in Python outside flopscope is a breach whether or not the cap is respected, and unusual residual usage attracts review.

**Data files are allowed.** Weights, lookup tables, calibration constants, any precomputed artifact, shipped alongside `estimator.py`, loaded in `setup()` (pickle-free `.npz` through flopscope; multi-file submissions must be packaged as a folder with `whest package --estimator .`). "Where the Sponsor considers that a bundled file is not genuinely data, the Sponsor may decline to treat it as data and may disqualify the submission." **This is the lever that makes offline compute unlimited**: anything universal (not tied to specific MLPs or seeds) that can be computed offline and shipped as numbers is legitimate.

**Review.** Automated checks on every submission; flagged ones go to agent-assisted validation; a person reviews the rest; review continues after grading, so a listed score can be invalidated later. Rules questions and mispricing reports go to arc-whestbench@aicrowd.com; "If you believe flopscope is mispricing an operation, report it rather than working around the accounting."

**Final evaluation.** The public leaderboard is not the ranking. Before the deadline each team designates **one** submission; after the deadline the grader re-executes it on a **fresh, unseen MLP suite generated from a private seed** not used during the open phases, and prize ranking is decided exclusively from that re-run. Overfitting to specific MLPs or seeds "will be penalized accordingly". Statistically close finishes are resolved first by generating additional private MLPs, then by all-layer MSE and effective compute. (Phase 1 used a "select up to 2 submissions on the Submissions page, team scope" mechanism; expect the same page for Phase 2's single designation.)

---

## 4. The FLOP meter (what things cost at $n = 1024$)

flopscope bills `charged = flop_cost × weight × dtype_rate × complex_factor`, weight tiers {0, 1, 4, 16}, **float64 costs 2× float32** (billing dtype follows NumPy promotion, so one stray float64 operand doubles the whole expression; `fnp.zeros`, `ones`, `eye` default to float64; `flops.stats.norm.*` promotes to float64).

| op | cost | at $n=1024$ |
|---|---|---|
| matmul `(M,K)@(K,N)` | $MN(2K-1)$ | $1024^2\cdot2047 \approx 2.15\times10^9$ for one $n\times n$ product |
| one forward pass of one input (16 layers) | $2dn^2$ | 33,554,432 arithmetic; 33,637,376 metered |
| budget in forward passes | | **65,536** (65,374 metered); 10 % floor = 6,553 passes |
| elementwise `+ − * /`, comparisons, `maximum`, `sqrt`, fills, copies, `astype`, `reshape` | 1 per element | $n^2$ = 1.05e6 per matrix |
| reductions `sum/max/min` | $N-1$; `mean` $N$; `var` $4N$ | |
| gathers, 3-arg `where` | 4 per element | |
| `exp`, `log`, trig, `x ** y` (even `** 2`) | 16 per element | write `x*x` |
| `stats.norm.pdf` / `cdf` | ≈ 54 / 96 per element (float64) | |
| `rng.standard_normal` | 16 per element float32 (32 float64); `uniform` 6 | 1024 inputs: 16k FLOPs each |
| sort/argsort/unique/searchsorted | ≈ $4N\lceil\log_2 N\rceil$ | |
| free | `zeros`, `empty`, views (`.T`, basic slicing), no-copy `asarray`, RNG construction, flopscope dispatch | |
| symmetric arrays | tagged symmetric arrays (`as_symmetric`, `X.T @ X`, Gram/covariance updates) are billed at the cheaper symmetric cost | overview, CHANGELOG |

Reference costs from the kit at this shape: mean propagation 86,639,616 FLOPs (0.004 %); covariance propagation 51,709,240,799 FLOPs (2.351 %), of which 99.66 % is one `einsum` per layer (the $O(n^3)$ covariance update). The full budget is 42× the covariance-propagation cost: at depth 16, a width of ≈ 3,570 would exhaust it with that method. A $K$-th-order cumulant method with a factorised $O(L^2 n^K)$ linear step (companion paper) costs, naively, $16^2\cdot1024^3 \approx 2.7\times10^{11}$ at $K=3$ (12 % of budget) and $2.8\times10^{14}$ at $K=4$ (130× over budget), so **$K = 3$ is the frontier of the exact-cumulant family at this shape and $K=4$ needs structural tricks (low rank, sparsity, sampling of the tensor)**.

The authoritative per-op table is flopscope's cost-model reference; `ctx.op_log` and `ctx.summary_dict()` give the metered breakdown per op with the resolved dtype.

---

## 5. The dataset

`aicrowd/arc-whestbench-public-2026`, revision **`v2-phase2`** (must be pinned; `v1-phase1` has the Phase 1 shape). Two independent splits with disjoint seeds: `mini` (100 MLPs, 7 GB parquet, 7 files) and `full` (1,000 MLPs, 70 GB, 63 files, ≈ 16 MLPs per file); also a `prepared/` Arrow export. Columns: `mlp_id`, `mlp_name`, `mlp_seed` (int64; protocol `whestbench_explicit_per_mlp_seeds` v3.0, the estimator seed is derived locally), `weights` float32 `[16,1024,1024]`, `all_layer_means` float32 `[16,1024]`, `final_means` float32 `[1024]`, `avg_variance` float64, `sampling_budget_breakdown` JSON. Baked 2026-08-13 on RTX 3090s with torch 2.4.1, whestbench 0.16.0, flopscope 0.12.0, bit-exact under the recorded determinism config; re-bake recipe published (`whest dataset bake --mlp-seeds ...`). The graded suite is 100 *other* MLPs (50 public + 50 holdout) and the final private suite is fresh again.

---

## 6. The local harness

`whest validate` (contract smoke test on a 4×2 MLP, no budget enforced), `whest run --runner local|subprocess|docker` (since whestbench 0.16.0 the defaults *are* the graded round: `--flop-budget 2199023255552 --wall-time-limit 120 --setup-timeout 5 --residual-wall-time-limit 0.4`; `--dataset hf://aicrowd/arc-whestbench-public-2026@v2-phase2 --split mini|full` scores against the baked truth; without `--dataset` it generates 10 MLPs and 200,000-sample truth, floor 3.7e-7), `whest package -o submission.tar.gz` (folder packaging for multi-file), `whest login` / `whest submit --watch`. `--format json` gives the machine-readable report (`run_config`, per-MLP `flops_used`, `residual_wall_time_s`, `budget_exhausted`, `time_exhausted`, `residual_wall_time_exhausted`, `error_code`, `mean_score_multiplier`, `mean_compute_utilization`, best/worst MLP scores). The subprocess runner reproduces the grader's transport and catches dirty imports, stdout writes and runaway memory; one worker serves the whole suite.

---

## 7. What this means for our programme

1. **The competition quantity is a state pushed forward through arrows.** $Y_{\ell}$ is the first moment of the law of $h_\ell$, and the ReLU acts exactly where the face structure of the earlier notes lives: $\mathbb E[\max(0,z_j)] = \mathbb E[z_j;\ z_j>0]$ is an expectation conditioned on the vertex $j$ being in the face. The cumulant-propagation methods of the companion paper are the "transport of a state's low-order statistics through one arrow" made explicit; depth is where they degrade ($\text{MSE}\sim c_K (L/n)^K$), and depth is where the earlier notes locate history dependence. That is the honest point of contact, and prong 2's dictionary v2 should be written against *this* realization (random He-initialised 1024×16 networks, Gaussian state), not against trained toy networks.
2. **Three levers only.** (a) accuracy per FLOP of the analytic part (which cumulants, which factorisation, which corrections); (b) a sampling component whose variance is reduced by the analytic part (control variates: the leaders' utilisation of 0.11–0.16 says they do this); (c) offline precomputation shipped as data (universal tables, calibrations fitted on thousands of random MLPs, which is exactly what cloud scale buys: the public `full` split plus any number of freshly baked MLPs).
3. **Robustness is a design constraint, not a QA step.** One failed MLP costs more than everything else combined (§2). The estimator needs a measured FLOP envelope with margin, a residual-time envelope with margin under the 2-vCPU solution process, an 8 GB memory envelope, idempotent `setup()` under 5 s, and a fallback path.
4. **Validation must mirror the private re-run.** Evaluate on fresh seeds (bake our own MLPs with the published recipe) as well as on `mini`/`full`; any tuning on the public MLPs is what the private re-run is designed to punish.

---

## Sources

- Phase 2 announcement, forum, 23 Aug 2026: https://discourse.aicrowd.com/t/phase-2-of-the-arc-white-box-estimation-challenge-is-live/18197
- Scoring mechanics post, forum, 17 Sep 2026: https://discourse.aicrowd.com/t/phase-2-the-compute-budget-is-not-a-lever-and-the-scoring-rule-is-why-adjusted-mse-x-max-0-1-util/18219
- Phase 1 posts on residual time, instrumented share and the budget-exceedance fix (18132, 18129, 18184, 18143, 18122).
- Starter kit: https://github.com/AIcrowd/whest-starterkit (README, CHANGELOG, docs/concepts/{allowed-code, scoring-model, problem-setup, ground-truth}.md, docs/reference/{rounds, estimator-contract, flopscope-primer}.md).
- Harness and meter: https://github.com/AIcrowd/whestbench, https://github.com/AIcrowd/flopscope.
- Dataset card: https://huggingface.co/datasets/aicrowd/arc-whestbench-public-2026 (revision `v2-phase2`).
- Leaderboard: https://www.aicrowd.com/challenges/arc-white-box-estimation-challenge-2026/leaderboards (read 1 Oct 2026).
- Companion paper: W. Wu, V. Lecomte, M. Winer, G. Robinson, J. Hilton, P. Christiano, *Estimating the expected output of wide random MLPs more efficiently than sampling*, arXiv:2605.05179.
- Challenge overview: https://www.aicrowd.com/challenges/arc-white-box-estimation-challenge-2026 ; Rules: https://www.aicrowd.com/challenges/arc-white-box-estimation-challenge-2026/challenge_rules (not fetchable from this container; binding text to be read on the site).
