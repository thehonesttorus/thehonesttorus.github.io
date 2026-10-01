# Submission stream: a validated fallback submission (504aldo V29 / V25)

*Stream deliverable for the ARC White-Box Estimation Challenge 2026, Phase 2. Written 1 Oct 2026. Every number below is measured on this stream's box and comes from a JSON under `results/`. `results/summary.md` lists every run (regenerate it with `python3 scripts/summarize.py`).*

## First: why "V29 patched, all dev MLPs failed" (16:00 status) and what fixed it

V29 was never broken. The **stock local harness** was the wrong judge. In-process, flopscope keeps every array in the solution process and charges every basic slice as residual Python. Two things followed. First, on the idle box V29's residual was 0.46–0.52 s on every MLP, over the 0.4 s cap. Second, under the subprocess runner's 8 GB RLIMIT_AS the second predict ran out of memory: one MemoryError fell back, and the next MLP killed the worker (WORKER_EOF). The grader is client/server: arrays live in the flopscope backend, and a slice is a counted round trip. Under grader-transport emulation (`scripts/grader_emul.py`) the same patched V29 passes 6/6 MLPs with residual 0.32–0.40 s.

**The fix for the thin residual margin** was to memoize views of pooled buffers (v29r3): round trips per MLP went from 28,246 to 16,670 (slices from 14,644 to 3,163), and residual from mean 0.383 / max 0.402 s to 0.289 / 0.317 s, with FLOPs and outputs bit-identical. Per call that is ≈ 0.012 → ≈ 0.017 s per 1k calls of client bookkeeping here (0.383 s / 28.2k vs 0.289 s / 16.7k). V29's published 13,121 "ops" counts compute ops only; slices are extra round trips.

## Verdict

| package | what it is | go/no-go | key numbers |
|---|---|---|---|
| **`packages/v29r3.tar.gz`** | upstream V29 + robustness wrapper + residual cut (view memo, 3 rounds) | **GO: primary upload** | Arithmetic is the graded V29's (#330093: LB 5.40e-9, no failed MLP): 0 FLOP difference against v29 and bit-identical MSE on all 6 dev MLPs; upstream + 32,767 FLOPs, the finiteness check. In grader-transport emulation the residual is mean 0.289 s, max 0.317 s, against upstream-equivalent v29's 0.383 / 0.402 s in the same interleaved window. Passes every robustness case. C/B 0.2526 (first MLP of a worker 0.2671). Setup ≤ 1.07 s. |
| `packages/v25.tar.gz` | upstream V25 + robustness wrapper | **GO: safe fallback** | Residual 0.11–0.13 s emulated, 0.19–0.23 s in-process (≥ 2× margin both ways). Wall 33–39 s, in-process RSS 3.0 GB. Same MSE as V29 (1.757e-8 vs 1.755e-8 on dev, noise-subtracted), but C/B 0.3667, so the expected LB is ≈ 7.8e-9 vs ≈ 5.4e-9. |
| `packages/v29r2.tar.gz`, `bundles/v29r` | earlier residual-cut rounds | superseded by v29r3 | v29r2: mean 0.298 s, max 0.350 s (same A/B window). |
| `packages/v29.tar.gz` | upstream V29 + robustness wrapper only | **no-go** | Residual equals upstream's: max 0.402 s, with 3 of 12 MLPs over the cap on this box in the A/B window. |
| upstream `estimator_v29.py` / `estimator_v25.py` as published | | **no-go** | Both fail robustness: V29 runs out of memory at 1024×32, V25 exhausts the FLOP budget at 1024×32, and both raise on scaled or all-positive weights. Upstream V29 also sits at 0.32–0.40 s residual here. |

**Upload order:** `v29r3` first. If its graded report shows any failed MLP (`n_failed_mlps > 0`, most plausibly a residual breach on a grader box slower than this one), upload `v25` next and designate that one.

**A rules question for the team before uploading:** these packages are team 504aldo's public, MIT-licensed code with local changes; the LICENSE file and an in-file MIT notice with attribution ship in every package. Whether a team may submit another team's public code as its own entry, and designate it for the final private re-run, has to be checked against the Rules page (not fetchable from this container).

## Question

Can we have a submission ready to upload the moment aicrowd.com is reachable that is known not to fail on any MLP? The failure modes are the FLOP budget, the 120 s wall, the 0.4 s residual, the 5 s setup, 8 GB of memory, and crashes on the grader's non-suite smoke shapes or odd inputs. It should also score near the public 504aldo numbers (V29: 0.2526 B, raw 2.13e-8, LB 5.40e-9; V25: 0.3667 B).

## Method

- **Candidates:** `estimators/estimator_v29.py` and `estimator_v25.py` from github.com/504aldo/whest-p2-cumulant-k3 @ `1a4083f` (MIT). Upstream sha256: V29 `86d9ca9b28e6…9d8e27`, V25 `c0ae6f12d27d…dd4b20`.
- **Harness:** whestbench 0.16.1 + flopscope 0.12.1 (`/root/whest`), with the graded caps (2^41 FLOPs, 120 s wall, 0.4 s residual, 5 s setup) and `--max-threads 2`. `scripts/measure_run.py` wraps `whest run --format json --profile`. It also samples the worker's peak RSS from outside the process and subtracts the truth noise floor `avg_variance / N`.
- **Grader-transport emulation (`scripts/grader_emul.py`).** The grader runs the solution against **flopscope-client** while a **flopscope-server** backend holds every array and does the numpy work (2 vCPU solution / 14 vCPU backend). The emulator runs flopscope-client 0.12.1 (`/root/fcli`) in a 1-thread client against flopscope-server 0.12.1 with 3 BLAS threads and the grader's 4 GiB per-array cap. Without that cap the server's default of 100 MB breaks V29 on its first large op. Setup, weight upload and fetch run in their own sessions, mirroring `subprocess_worker`. Residual is the client's `wall − dispatch` (flopscope-client `_decompose_timing`), the quantity the 0.4 s cap applies to on the grader. FLOPs and MSE agree with the in-process harness to the bit, apart from a +16.7M-FLOP client-side billing difference on V29 (0.003% of C).
- **Dev set:** `whest dataset bake --width 1024 --depth 16 --n-mlps 6 --n-samples 2000000`, seeds `[7301001 … 7301006]`, 35 min on 3 threads, 385 MB (not committed). Truth noise floor 3.64e-8, larger than the estimators' error, hence the subtraction (standard error of the noise-subtracted 6-MLP mean ≈ 1e-9).
- **Robustness sets:** 1 MLP each at 1024×32, 512×16, 1024×4 and 256×8 (seed 4242), plus the 1024×16 MLP with W×10 and with |W| (all-positive), N = 2000, because only failures and caps matter here. Built by `whest dataset bake` and a parquet edit (`scripts/make_adversarial.py`). `scripts/make_datasets.sh` rebuilds every dataset used here.
- **Box:** 4 vCPU, 16 GB. "Idle" = nothing else running. The box's speed drifted by ~10% over the day, so residual comparisons that matter were run **interleaved in the same window** (section 3b).

## Results

### 1. Robustness (subprocess runner, RLIMIT_AS 8 GB; `results/robust_contended/`, `results/robust_idle/`)

| set | V29 upstream | V25 upstream | v29 / v29r / v29r2 / v29r3 | v25 |
|---|---|---|---|---|
| 1024×32 | **FAIL**: MemoryError (RSS 6.9 GB) | **FAIL**: FLOP budget (C/B 0.9985) | ok: fallback, C/B 0.122, res 0.02–0.03 s, wall 1.6–2.7 s | ok: C/B 0.122, res 0.03 s |
| 512×16 | ok (res 0.375 s loaded; one earlier loaded run **FAILED** at 0.406 s) | ok (res 0.16 s) | ok: fallback, C/B 0.0075, res 0.01 s | ok: C/B 0.0075 |
| 1024×4 | ok | ok | ok: C/B 0.013 | ok |
| 256×8 | ok | ok | ok: C/B 0.0004 | ok |
| 1024×16, W×10 | **FAIL**: SymmetryError (NaN) | **FAIL**: SymmetryError | ok: chain raises → fallback; C/B 0.106, res 0.14–0.15 s, wall 25–26 s | ok: C/B 0.063 |
| **256×32 (grader smoke shape, whestbench #149)** | **FAIL**: residual 0.487 s (in-process, idle) | ok: C/B 0.037, res 0.26 s | ok (v29r3): fallback, C/B 0.0019, res 0.016 s | ok: C/B 0.0019, res 0.021 s |
| 1024×16, \|W\| | **FAIL**: SymmetryError (NaN) | **FAIL**: SymmetryError | ok: fallback; C/B 0.091, res 0.12 s | ok: C/B 0.106 |

The `inf` MSE on |W| is the harness squaring true means of order 1e23 in float32. It isn't an estimator failure (`n_failed = 0`, finite predictions).

**The robustness patch (`scripts/patch_safe.py`, diffs in `patches/`)** renames the upstream `predict` to `_predict_main` and adds a new `predict`:
- any shape other than 1024×16 (the grader's depth-32 smoke test, say) goes straight to `_safe_predict`, a float64 covariance propagation (ReLU Gaussian closure, ≈ 0.12 B at 1024×32, a few MB);
- at 1024×16 the upstream chain runs unchanged. It falls back only if the chain raises something other than a flopscope budget/time exhaustion, or returns non-finite values (one `isfinite` + `all` = 32,767 FLOPs).

### 2. Parity (client/server emulation, dev set, 6 MLPs)

| pair | FLOPs difference per MLP | final-layer MSE |
|---|---|---|
| v29 vs V29 upstream | +32,767 on all 6 | bit-identical on all 6 |
| v25 vs V25 upstream | +32,767 on all 6 | bit-identical on all 6 |
| v29r, v29r2, v29r3 vs v29 | 0 on all 6 | bit-identical on all 6 |
| packaged (extracted tarball) vs bundle: v29, v29r2, v29r3, v25 | 0 (first 2 dev MLPs) | bit-identical; sha256 of the packaged `estimator.py` equals the bundle's |

### 3a. Dev set, grader-transport emulation, all runs (`results/emul/`, `results/emul_r23/`)

| bundle (run) | fails | residual per MLP (s), MLPs 1–6 | wall max (s) | server RSS | MSE − noise |
|---|---|---|---|---|---|
| V29 upstream | 0/6 | 0.321, 0.348, 0.354, 0.349, 0.332, 0.329 | 99 | 5.6 GB | 1.7548e-8 |
| v29 (2 runs) | 0/6 | 0.343–0.392 / 0.321–0.355 | 107 / 98 | 5.6 GB | 1.7548e-8 |
| v29r (3 runs) | 0/6 | 0.277–0.334 / 0.288–0.334 / 0.280–0.351 | ≤ 102 | 5.6 GB | 1.7548e-8 |
| v29r2 (3 runs) | 0/6 | 0.257–0.325 / 0.284–0.334 / 0.285–0.350 | ≤ 105 | 5.6 GB | 1.7548e-8 |
| V25 upstream | 0/6 | 0.112–0.124 | 38 | 1.9 GB | 1.7566e-8 |
| v25 (2 runs) | 0/6 | 0.114–0.132 / 0.115–0.126 | 39 / 37 | 1.9 GB | 1.7566e-8 |

C/B is data-independent: steady state is 555,400,514,030 FLOPs = 0.25257 B. A worker's first predict runs Strassen level 4 (upstream design) at 0.2671 B, and its second 0.2540 B (pool growth is billed).

### 3b. Interleaved A/B, same time window (`results/emul_ab/`; order v29, v29r2, v29r3, then reversed)

| bundle | residual, run 1 (s) | residual, run 2 (reverse order) (s) | mean of 12 | max | over 0.4 s |
|---|---|---|---|---|---|
| v29 (= upstream residual) | 0.372, 0.401, 0.402, 0.371, 0.389, 0.402 | 0.349, 0.366, 0.387, 0.387, 0.398, 0.373 | 0.383 | 0.402 | **3** |
| v29r2 | 0.350, 0.310, 0.287, 0.272, 0.273, 0.285 | 0.329, 0.326, 0.293, 0.279, 0.293, 0.286 | 0.298 | 0.350 | 0 |
| **v29r3** | 0.317, 0.300, 0.295, 0.304, 0.291, 0.266 | 0.276, 0.300, 0.288, 0.273, 0.266, 0.290 | **0.289** | **0.317** | 0 |

Wall per MLP was 98–107 s with the emulated server's 3 BLAS threads, against the 120 s cap. The grader's backend has 14 threads, and V29's arithmetic was graded without a time failure. This is the second-thinnest margin and can't be measured more faithfully here.

### 4. Dev set, in-process harness (`results/dev6_idle/`; every array lives in the solution process)

| run | fails | residual per MLP (s) | wall (s) | peak RSS | MSE − noise |
|---|---|---|---|---|---|
| v25 subprocess, 3 runs | 0/6 ×3 | 0.195–0.231 | 33–37 | 3.0 GB | 1.7564e-8 |
| v29 subprocess | **6/6** | 0.463–0.484 | 70–85 | 6.9 GB (2 MLPs killed by RLIMIT_AS → WORKER_EOF) | – |
| v29 local | **6/6** (all residual) | 0.475–0.520 | 77–87 | 6.3 GB | – |
| v29r2 local | **6/6** (all residual) | 0.414–0.480 | 75–91 | 6.3 GB | – |
| v29r3 local | **6/6** (all residual) | 0.407–0.455 | 73–89 | 6.3 GB | – |

In-process, every basic slice is pure Python and lands in the residual, and every array counts against the 8 GB address-space limit. On the grader the arrays live in the flopscope server and a slice is a counted round trip. **The stock local harness therefore can't validate V29-family bundles on this box** (V29's graded run shows the grader behaves like the emulator, not like this). Use `scripts/grader_emul.py`. The memo also cuts the in-process residual (v29 0.475–0.520 s → v29r3 0.407–0.455 s), but not below the cap: in-process, the ~13.5k compute ops alone cost more client-side Python than in client mode.

### 5. Setup window (`results/setup/`)

Worker spawn + imports + module load + `setup()`, as SubprocessRunner times it (5 spawns each): v25 ≤ 1.04 s, v29 ≤ 1.08 s, v29r ≤ 1.04 s, v29r2 ≤ 1.20 s, v29r3 ≤ 1.07 s, against the 5 s cap. `setup()` alone in client mode takes 0.10–0.13 s.

### 6. Residual margin and how it was cut (task 6)

A cProfile of one V29 predict in client mode (`grader_emul.py --profile-out`) shows 28,246 server round trips, 14,644 of them basic slices. The estimator's own Python is ~0.14 s even under the profiler. The rest of the residual is per-call client bookkeeping outside the timed dispatch span: context-manager setup, frame checks, RemoteArray construction and its weakref finaliser. **Residual therefore scales with the number of flopscope calls, and half of V29's calls re-slice the same persistent pooled buffers.** Three rounds, each the same ops in the same order (bit-identical, section 2):

| round | change | round trips / slices per predict (MLP 2) | A/B residual mean / max |
|---|---|---|---|
| (v29) | none | 28,246 / 14,644 | 0.383 / 0.402 s |
| 1: v29r (`scripts/patch_views.py`) | `_Strassen` memoizes views of its own pooled buffers: batch prefixes, the 7 combo columns, quadrants, swapaxes, `[None]`, `[0]`. Only registered pool views are memoized; any other operand is sliced fresh. | 20,242 / 6,714 | (not in A/B) |
| 2: v29r2 (`scripts/patch_views2.py`) | `_predict_core` / `_dslices` / `_hub2` pass memoized slices of the permanent `_Pool` buffers (leg and factor slabs, `lap4`, hub, `abbuf` rows), so top-level operands hit the memo too | 19,169 / 5,641 | 0.298 / 0.350 s |
| 3: v29r3 (`scripts/patch_views3.py`) | per-key invalidation: the first allocation of a pool key no longer wipes the memo, and growing key K drops only K's views. Rounds 1–2 wiped everything 133 times in a worker's first predict and 37 times in its second. | 16,670 / 3,163 | **0.289 / 0.317 s** |

Server peak RSS is unchanged (5.56–5.57 GB), so the memo retains no old buffers.

## Grader lessons for any design

1. **Validate under client/server, not only in-process.** The grader runs the solution on flopscope-client against a flopscope-server backend. In-process `whest run` mis-states both residual and memory: here, upstream V29 read 0.47–0.52 s in-process but 0.32–0.40 s emulated, and 6.9 GB in the solution process vs ~2.4 GB client + 5.6 GB server. `scripts/grader_emul.py` emulates it; set the server's `FLOPSCOPE_MAX_ARRAY_BYTES` to the grader's 4 GiB (the default of 100 MB kills large ops).
2. **Residual ≈ (number of flopscope calls) × (client bookkeeping per call)**, about 10–17 µs per call on this box, not the estimator's own arithmetic. Every call counts, free views included: slices, `[None]`, swapaxes, reshape. Budget calls, not FLOPs. Keep views of persistent buffers and reuse them; write into pooled buffers with `out=`; disable `gc` inside `predict` (V29 does).
3. **A worker's first and second predicts are the residual worst case** (pool allocation, cache warm-up, lazy library init). Each submission has 5–15 workers, so each pays this 5–15 times. Warm op signatures in `setup()` (cheap); don't allocate GBs there.
4. **Memory:** in-process runs hit RLIMIT_AS 8 GB long before the grader would (arrays live server-side there). Still keep any single array under 4 GiB.
5. **Setup:** the 5 s window includes interpreter spawn and imports (≈ 1.0 s measured here). `setup()` itself should only load and warm up.
6. **Smoke shapes:** the grader also runs a non-suite MLP (256 wide × 32 deep, HF `aicrowd/whestbench-smoke-mlp`, whestbench #149). One failure there fails the whole submission, and `whest validate` doesn't catch it. Gate the suite-specific path on `(width, depth) == (1024, 16)` and send everything else to a cheap, shape-generic float64 path. Upstream V29 fails residual at 256×32 in-process; V25/V29 run out of FLOP budget or memory at 1024×32.
7. **Fail soft at the suite shape:** wrap the main chain. On any exception other than flopscope budget/time exhaustion, or on a non-finite output, fall back to the cheap path. One zeroed MLP costs ≈ 9e-3 of mean score, i.e. everything.
8. **Client parity traps:** no `x.shape = ...` (flopscope #267, raises on the grader; neither bundle here has it), arrays are immutable on the client (only `out=` writes), pass `str` paths, keep stdout clean (it is the IPC pipe).
9. **Accounting hygiene (16 Sep fair-accounting rule):** V29's discounts are flopscope-validated symmetry tags (`as_symmetric` raises on a non-symmetric input; we saw it fire on NaN weights), the documented aliased-Gram 0.5× on layer 0's `W^T W`, and Strassen-Winograd as flopscope ops (permitted 11 Sep). No packing and no reliance on stale tags was found. A diagonal write voids a tag, so V29 re-tags after `fill_diagonal`.

## Open issues

1. **The grader's client speed is unknown.** Upstream V29 passed on the grader. v29r3 issues the same compute ops with ~11.6k fewer slice round trips per predict, so its grader residual should sit below the graded V29's by roughly the 25% measured here. A grader box ≥ 25% slower than the one that graded #330093 could still breach on the V29 family. V25 is the hedge (≥ 3× margin in emulation).
2. The 120 s wall (98–107 s emulated with a 3-thread backend) can't be measured faithfully here (section 3b).
3. The stock in-process harness fails V29-family bundles on this box (residual, and RLIMIT_AS on the second predict). That is expected (section 4) but means local re-validation needs the emulator.
4. Robustness sets used N = 2000 truth, so only failure flags are meaningful there, not MSE.
5. The rules question on submitting third-party MIT code (see Verdict).

## Upload instructions

Prerequisites: the team is registered and frozen by **2 Oct 2026 23:59 UTC**; submissions close 17 Oct 23:59 UTC; 10 submissions per team per UTC day (failures count).

```bash
# 0. tools (any machine with Python 3.10+; the grader pins whestbench 0.16.0 / flopscope 0.12.0, the kit 0.16.1 / 0.12.1)
python3 -m venv ~/whest && ~/whest/bin/pip install "whestbench==0.16.1" "flopscope==0.12.1"
export PATH=~/whest/bin:$PATH
cd notes/streams/submission

# 1. check the committed archive is the validated one
sha256sum packages/v29r3.tar.gz     # 331047facafc78d5f31a843e6271dc2f5851d26b6fc04b54bda5f624dded53f1
whest validate-package packages/v29r3.tar.gz
#    (or rebuild it: cd bundles/v29r3 && whest package --estimator . --output ../../packages/v29r3.tar.gz --yes;
#     the manifest timestamp changes the tarball hash, the estimator.py sha256 stays 60fa9c12...b589)

# 2. log in and submit (the API key is on the AIcrowd profile page; the challenge slug defaults to
#    arc-white-box-estimation-challenge-2026)
whest login --api-key "$AICROWD_API_KEY"
whest submit packages/v29r3.tar.gz \
    --description "K3 cumulant chain (504aldo V29, MIT) + robustness wrapper + residual cut" --watch --watch-timeout 1800
#    equivalent from the folder (preview first; the dry run was checked: 2 files, 31.8 KB archive):
#      (cd bundles/v29r3 && whest submit --estimator . --yes --dry-run)
#      (cd bundles/v29r3 && whest submit --estimator . --yes --description "..." --watch --watch-timeout 1800)
#    alternative: upload packages/v29r3.tar.gz through the challenge's Submissions page.

# 3. read the graded report: expect n_failed_mlps = 0, mean_compute_utilization ~0.253-0.256,
#    final_layer_mse ~2.1e-8, adjusted ~5.4e-9 (V29's graded #330093: 5.40e-9).
#    If n_failed_mlps > 0: submit the fallback and designate it instead:
sha256sum packages/v25.tar.gz       # 212e34e60f90c282c417b706c404583d067c47cf5374a46805ca53aec760f963
whest submit packages/v25.tar.gz \
    --description "K3 cumulant chain (504aldo V25, MIT) + robustness wrapper" --watch --watch-timeout 1800
#    expected: n_failed_mlps = 0, utilization 0.3667, adjusted ~7.8e-9.

# 4. designate the chosen submission for the final private re-run on the Submissions page (one per team).
```

Package contents (each): `estimator.py` (MIT notice + attribution in the header), `LICENSE` (504aldo, MIT), `manifest.json`. No data files, so whestbench issue #119 (assets dropped by packaging) does not apply.

## Files

- `bundles/{v29r3,v29r2,v29r,v29,v25}/`: estimator.py + LICENSE (upload folders).
- `packages/*.tar.gz`: archives built by `whest package --estimator . --yes`, checked with `whest validate-package`, extracted and re-run.
- `patches/`: unified diffs, upstream → each bundle and round to round.
- `scripts/`: `patch_safe.py`, `patch_views{,2,3}.py` (all asserted textual patches, reproducible from upstream: `patch_safe.py UPSTREAM bundles/v29/estimator.py`, then `patch_views.py` v29→v29r, `patch_views2.py` v29r→v29r2, `patch_views3.py` v29r2→v29r3), `make_datasets.sh`, `make_adversarial.py`, `measure_run.py`, `grader_emul.py`, `robustness.sh`, `validate_dev*.sh`, `emul_dev.sh`, `setup_time.py`, `verify_package.sh`, `summarize.py`.
- `results/`: every run's JSON (with per-MLP FLOPs, timings, MSE), `summary.md`.
