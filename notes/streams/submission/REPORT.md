# Submission stream: a validated fallback submission (504aldo V29 / V25)

*Stream deliverable for the ARC White-Box Estimation Challenge 2026, Phase 2. Written 1 Oct 2026. Every number below is measured on this stream's box and comes from a JSON under `results/`. `results/summary.md` lists every run.*

## Verdict

| bundle | what it is | go/no-go | why |
|---|---|---|---|
| **`bundles/v29r2`** | V29 + robustness wrapper + view-memo residual cut (2 rounds) | **GO: primary upload** | Same arithmetic as the graded V29 (#330093: LB 5.40e-9, no failed MLP): 0 FLOP difference and bit-identical MSE on all 6 dev MLPs against the upstream file. It has the lowest V29 residual: 0.26 s steady state and 0.33 s on a worker's first MLP in grader-transport emulation, against 0.33–0.39 s for upstream V29. It passes every robustness case. |
| `bundles/v25` | V25 + robustness wrapper | **GO: safe fallback** | 3× margin on the residual cap (0.11–0.13 s emulated, 0.19–0.23 s in-process), 3× on wall (33–39 s), 3.0 GB in-process. MSE equals V29's (1.757e-8 vs 1.755e-8 on dev, noise-subtracted), but C/B is 0.3667, so the expected LB is ≈ 7.8e-9 vs ≈ 5.4e-9. |
| `bundles/v29r` | V29 + wrapper + view memo (round 1 only) | superseded by v29r2 | Steady residual 0.28–0.33 s. |
| `bundles/v29` | V29 + robustness wrapper only | no-go (use v29r2) | Same residual as upstream (0.32–0.39 s emulated). Do not upload the **upstream** V29 or V25 files as they are: both fail the robustness cases (below). |

**Upload order:** `packages/v29r2.tar.gz` first. If its graded run shows any failed MLP (`n_failed_mlps > 0`, most plausibly a 0.4 s residual breach on a slower grader box), upload `packages/v25.tar.gz` next and designate that one.

**Before uploading, a rules question the team must answer** (not a technical one): these bundles are team 504aldo's public, MIT-licensed code with small local changes. Attribution is kept: the LICENSE file and an in-file MIT notice ship in every package. Whether a team may submit another team's public code as its own entry (and designate it for the final private re-run) is for the team to check against the Rules page (§5 / single-submission clause), which could not be fetched from this container.

## Question

Can we have a submission ready to upload the moment aicrowd.com is reachable that is known not to fail on any MLP? The failure modes are the FLOP budget, the 120 s wall, the 0.4 s residual, the 5 s setup, 8 GB of memory, and crashes on off-suite or odd inputs. It should also score near the public 504aldo numbers (V29: 0.2526 B, raw 2.13e-8; V25: 0.3667 B).

## Method

- **Candidates:** `estimators/estimator_v29.py` and `estimator_v25.py` from github.com/504aldo/whest-p2-cumulant-k3 @ `1a4083f` (MIT). Upstream sha256: V29 `86d9ca9b…8e27`, V25 `c0ae6f12…4b20`.
- **Harness:** whestbench 0.16.1 + flopscope 0.12.1 (`/root/whest`), with the graded caps (2^41 FLOPs, 120 s wall, 0.4 s residual, 5 s setup) and `--max-threads 2`. `scripts/measure_run.py` wraps `whest run --format json --profile`. It also samples the worker's peak RSS from outside the process and subtracts the truth noise floor `avg_variance / N` from the final-layer MSE.
- **Grader-transport emulation (`scripts/grader_emul.py`).** The grader runs the solution against **flopscope-client** while a **flopscope-server** backend holds the arrays and does the numpy work (2 vCPU solution / 14 vCPU backend). The emulator installs flopscope-client 0.12.1 in `/root/fcli` and runs flopscope-server 0.12.1 with 3 BLAS threads and the grader's 4 GiB per-array cap (the server's default of 100 MB breaks V29). The estimator runs in a 1-thread client, with setup and weight upload in their own sessions, mirroring `subprocess_worker`. Residual is the client's own `wall − dispatch` (flopscope-client `_decompose_timing`), the quantity the 0.4 s cap applies to on the grader. FLOPs and MSE agree with the in-process harness to the bit, apart from a +16.7M-FLOP client-side billing difference on V29 (0.0008% of C).
- **Dev set:** `whest dataset bake --width 1024 --depth 16 --n-mlps 6 --n-samples 2000000`, seeds 7301001–7301006 (385 MB, not committed; rebuild with `scripts/` + the seeds). Truth noise floor 3.64e-8, larger than the estimators' error, hence the subtraction (standard error of the noise-subtracted 6-MLP mean ≈ 1e-9).
- **Robustness sets:** 1 MLP each at 1024×32, 512×16, 1024×4 and 256×8, plus the 1024×16 MLP with W×10 and with |W| (all-positive); N = 2000, because only failures and caps matter here.
- **Box:** 4 vCPU, 16 GB. "Idle" runs had nothing else running. The robustness sweep of the upstream files ran while a bake used 3 cores ("loaded").

## Results

### 1. Robustness (subprocess runner, RLIMIT_AS 8 GB; `results/robust_contended/`, `results/robust_idle/`)

| set | V29 upstream | V25 upstream | v29 / v29r / v29r2 (patched) | v25 (patched) |
|---|---|---|---|---|
| 1024×32 | **FAIL**: MemoryError (RSS 6.9 GB) | **FAIL**: FLOP budget (C/B 0.9985) | ok: fallback, C/B 0.122, res 0.02–0.03 s, wall 1.6–2.7 s | ok: C/B 0.122, res 0.03 s |
| 512×16 | ok (res 0.375 s loaded; one earlier loaded run **FAILED** at 0.406 s) | ok (res 0.16 s) | ok: fallback, C/B 0.0075, res 0.01 s | ok: C/B 0.0075 |
| 1024×4 | ok | ok | ok: C/B 0.013 | ok |
| 256×8 | ok | ok | ok: C/B 0.0004 | ok |
| 1024×16, W×10 | **FAIL**: SymmetryError (NaN) | **FAIL**: SymmetryError | ok: chain raises, fallback; C/B 0.106, res 0.15 s, wall 26 s | ok: C/B 0.063 |
| 1024×16, \|W\| | **FAIL**: SymmetryError (NaN) | **FAIL**: SymmetryError | ok: fallback; C/B 0.091, res 0.12 s | ok: C/B 0.106 |

The `inf` MSE on |W| is the harness squaring true means of order 1e23 in float32. It isn't an estimator failure (`n_failed = 0`, finite predictions).

**The robustness patch (`scripts/patch_safe.py`, diffs in `patches/`)** renames the upstream `predict` to `_predict_main` and adds a new `predict`:
- any shape other than 1024×16 (e.g. the grader's depth-32 smoke test) goes straight to `_safe_predict`, a float64 covariance propagation (ReLU Gaussian closure, ≈ 0.12 B at 1024×32, a few MB);
- at 1024×16 the upstream chain runs unchanged. It falls back to `_safe_predict` only if the chain raises something other than a flopscope budget/time exhaustion, or returns non-finite values (one `isfinite` + `all` = 32,767 FLOPs, 0.000001% of C).

### 2. Parity: patched bundles vs upstream (client/server emulation, dev set, 6 MLPs)

| pair | FLOPs difference per MLP | final-layer MSE |
|---|---|---|
| v29 vs V29 upstream | +32,767 (the finiteness check) on all 6 | bit-identical on all 6 |
| v25 vs V25 upstream | +32,767 on all 6 | bit-identical on all 6 |
| v29r vs v29 | 0 on all 6 | bit-identical on all 6 |
| v29r2 vs v29 | 0 on all 6 | bit-identical on all 6 |

### 3. Dev set, grader-transport emulation (`results/emul/`, `results/emul_r23/`)

| bundle (repeat) | fails | C/B (first MLP / steady) | residual per MLP (s), MLPs 1–6 | wall max (s) | server RSS | MSE − noise |
|---|---|---|---|---|---|---|
| V29 upstream | 0/6 | 0.2671 / 0.2526 | 0.321, 0.348, 0.354, 0.349, 0.332, 0.329 | 99 | 5.6 GB | 1.7548e-8 |
| v29 (r1) | 0/6 | 0.2671 / 0.2526 | 0.343, 0.392, 0.390, 0.379, 0.373, 0.378 | 107 | 5.6 GB | 1.7548e-8 |
| v29 (r2) | 0/6 | same | 0.321, 0.355, 0.350, 0.329, 0.355, 0.350 | 98 | 5.6 GB | 1.7548e-8 |
| v29r (r1) | 0/6 | same | 0.334, 0.314, 0.301, 0.287, 0.288, 0.277 | 100 | 5.6 GB | 1.7548e-8 |
| v29r (r2) | 0/6 | same | 0.334, 0.319, 0.312, 0.327, 0.288, 0.298 | 102 | 5.6 GB | 1.7548e-8 |
| v29r (r3) | 0/6 | same | 0.351, 0.321, 0.311, 0.280, 0.282, 0.292 | 99 | 5.6 GB | 1.7548e-8 |
| **v29r2 (r1)** | 0/6 | same | **0.325, 0.301, 0.286, 0.257, 0.258, 0.261** | 98 | 5.6 GB | 1.7548e-8 |
| V25 upstream | 0/6 | 0.3667 | 0.112, 0.114, 0.120, 0.124, 0.121, 0.120 | 38 | 1.9 GB | 1.7566e-8 |
| v25 (r1) | 0/6 | 0.3667 | 0.132, 0.122, 0.124, 0.129, 0.114, 0.122 | 39 | 1.9 GB | 1.7566e-8 |
| v25 (r2) | 0/6 | 0.3667 | 0.119, 0.126, 0.123, 0.120, 0.117, 0.115 | 37 | 1.9 GB | 1.7566e-8 |

C/B is data-independent: steady 555,400,514,030 FLOPs = 0.25257 B. The first predict of a worker runs Strassen level 4 (upstream design) at 0.2671 B, and the second 0.2540 B (pool growth is billed). Wall: the emulated server had 3 BLAS threads; the grader's backend has 14. The 98–107 s wall is therefore a pessimistic bound, but it is the thinnest margin after the residual. V29 was graded on the real backend without a time failure.

### 4. Dev set, in-process harness (`results/dev6_idle/`; arrays live in the solution process)

| run | fails | residual per MLP (s) | wall (s) | peak RSS | MSE − noise |
|---|---|---|---|---|---|
| v25, subprocess r1/r2/r3 | 0/6 ×3 | 0.195–0.231 | 33–37 | 3.0 GB | 1.7564e-8 |
| v29, subprocess r1 | **6/6** | 0.463–0.484 | 70–85 | 6.9 GB (2 MLPs killed by RLIMIT_AS → WORKER_EOF) | – |
| v29, local r1 | **6/6** (all residual) | 0.475–0.520 | 77–87 | 6.3 GB | – |

In-process, every basic slice is pure Python and lands in the residual, and every array counts against the 8 GB address-space limit. On the grader the arrays live in the flopscope server and a slice is a counted round trip. The in-process V29 numbers therefore do **not** predict the grader, and V29's graded run confirms that (no failure, 100 MLPs). They do show that V29 can't be validated with the stock local harness on this box. Use the emulator.

### 5. Setup window (`results/setup/`)

Worker spawn + imports + module load + `setup()`, as SubprocessRunner times it (5 spawns each): v25 ≤ 1.04 s, v29 ≤ 1.08 s, v29r ≤ 1.04 s, against the 5 s cap. `setup()` alone in client mode takes 0.10–0.13 s (V29's ~45 warm-up ops are round trips).

### 6. Residual margin and how v29r/v29r2 cut it (task 6)

Profile of one steady V29 predict in client mode (`cProfile`, `scripts/grader_emul.py --profile-out`): 28,246 server round trips, of which 14,644 are basic slices (`__getitem__`). The estimator's own Python is ~0.14 s even under the profiler. The rest of the residual is per-call client bookkeeping outside the timed span (context-manager setup, frame checks, RemoteArray construction and its weakref finaliser). **So residual scales with the call count, and half the calls re-slice the same persistent pooled buffers.**

- **v29r** (`scripts/patch_views.py`): `_Strassen` memoizes the views of its own pooled buffers (batch prefixes, the 7 combo columns, quadrants, swapaxes, `[None]`, `[0]`). Only registered pool views are memoized; any other operand is sliced fresh. Any pool growth clears the memo. Round trips per steady predict fall from 28.2k to 20.2k.
- **v29r2** (`scripts/patch_views2.py`): `_predict_core` / `_dslices` / `_hub2` hand the kernel memoized slices of the permanent `_Pool` buffers (leg slabs, factor slabs, `lap4`, hub, `abbuf` rows), so the kernel's memo also hits for top-level operands. Per-layer arrays (W, Qc, QU) are untouched.
- Effect (client mode, this box): steady residual 0.33–0.39 s (v29) → 0.28–0.33 s (v29r) → 0.26 s (v29r2), and 0.32–0.35 s → 0.33–0.35 s → 0.325 s on a worker's first MLP. FLOPs and outputs are unchanged (section 2). The memo is fully warm only from a worker's 4th MLP: pools grow during MLPs 1–2 (first call at Strassen level 4, then level 5).
- Margin left (worst MLP, this box): v29r2 0.325 s, i.e. 19% under the cap; V25 0.13 s, 69% under.

## Open issues

1. **The grader's client speed is unknown.** Upstream V29 passed on the grader, and v29r2 does strictly less client-side work than upstream (same ops minus ~8k+ slice round trips). Absent grader drift, v29r2's grader residual is therefore below the graded V29's. A slower grader box would hit V29-family bundles first. V25 is the hedge.
2. A worker's first MLP (0.33 s here) is now v29r2's worst. Further cuts would need the pools allocated in `setup()` (moves GBs of server-side allocation into a 5 s window that runs 5–15 times per submission: not worth the setup-timeout risk) or fewer compute ops (changes FLOPs).
3. 120 s wall: 98–107 s in emulation with a 3-thread backend; the grader's 14-thread backend should be well under. Not independently verifiable here.
4. The in-process local harness can't validate V29-family bundles on this box (residual and RLIMIT_AS). Use `scripts/grader_emul.py`.
5. Rules question on submitting third-party MIT code (see Verdict).

## Upload instructions

(See the end of this file, updated with package hashes once `scripts/verify_package.sh` has run.)
