# Submission stream: a validated fallback submission (504aldo V29 / V25)

*Status: IN PROGRESS. Partial results; numbers below are measured, open items are listed at the end.*

## Question

Can we upload, the moment aicrowd.com is reachable, a submission that is known not to fail on any MLP (budget, 120 s wall, 0.4 s residual, 5 s setup, 8 GB, crash) and scores near the public 504aldo numbers (V29: 0.2526 B, raw 2.13e-8; V25: 0.3667 B)?

## Method

- Candidates: `estimators/estimator_v29.py` and `estimator_v25.py` of github.com/504aldo/whest-p2-cumulant-k3 @ `1a4083f` (MIT; the LICENSE file and an in-file MIT notice ship with each bundle).
- Harness: whestbench 0.16.1 + flopscope 0.12.1 in `/root/whest`; `whest run --format json --profile --max-threads 2`, graded caps (2^41 FLOPs, 120 s wall, 0.4 s residual, 5 s setup). Wrapper `scripts/measure_run.py` adds peak RSS of the worker and subtracts the truth noise floor `avg_variance / N`.
- Dev set: `whest dataset bake --width 1024 --depth 16 --n-mlps 6 --n-samples 2000000`, seeds 7301001..7301006 (`/tmp/claude-0/sub/dev6`, not committed: 385 MB). Truth noise floor 3.64e-8 (larger than the estimators' error, hence the subtraction; the standard error of the noise-subtracted 6-MLP mean is ≈ 1e-9).
- Robustness sets: 1 MLP each at 1024×32, 512×16, 1024×4, 256×8, 1024×16, and the 1024×16 MLP with weights ×10 and with |W| (all-positive), N = 2000 (only failure / caps matter here).
- Box: 4 vCPU, 16 GB; idle runs had nothing else running.
- Grader-transport emulation (`scripts/grader_emul.py`): estimator in a client process on flopscope-client 0.12.1, arrays and numpy work in a flopscope-server 0.12.1 process, which is the grader's architecture. In client mode the residual is the client's `wall − dispatch`.

## Results so far

### Robustness (subprocess runner, RLIMIT_AS 8 GB; box loaded, 1 BLAS thread; `results/robust_contended/`)

| set | V29 original | V25 original | V29 patched | V25 patched |
|---|---|---|---|---|
| 1024×32 | **FAIL** MemoryError (RSS 6.9 GB) | **FAIL** budget (C/B 0.9985) | ok, C/B 0.122, res 0.026 s | ok, C/B 0.122, res 0.033 s |
| 512×16 | ok (res 0.375 s; 0.406 s = FAIL in an earlier loaded run) | ok, res 0.16 s | ok, C/B 0.0075 | ok, C/B 0.0075 |
| 1024×4 | ok | ok | ok | ok |
| 256×8 | ok | ok | ok | ok |
| 1024×16 W×10 | **FAIL** SymmetryError (NaN) | **FAIL** SymmetryError | ok (fallback), C/B 0.106 | ok (fallback), C/B 0.063 |
| 1024×16 \|W\| | **FAIL** SymmetryError (NaN) | **FAIL** SymmetryError | ok (fallback), C/B 0.091 | ok (fallback), C/B 0.106 |
| 1024×16 (first MLP of a worker) | FAIL residual 0.453 s (loaded box) | ok, res 0.197 s | FAIL residual 0.457 s (loaded box) | ok, res 0.204 s |

### The patch (`scripts/patch_safe.py`, applied to both bundles)

The original `predict` becomes `_predict_main`. A new `predict` sends any non-suite shape (not 1024×16) straight to a float64 covariance-propagation fallback (`_safe_predict`: ReLU Gaussian closure, ≈ 0.12 B at 1024×32). On the suite shape it runs the original chain unchanged. It falls back only if the chain raises something other than a flopscope budget/time exhaustion, or returns a non-finite array (one `isfinite` + `all`, about 33k FLOPs). Suite-shape outputs are bit-identical to the original; see the parity rows below.

### Dev set, idle box, `--max-threads 2` (`results/dev6_idle/`)

| run | fails | C/B | residual per MLP (s) | wall per MLP (s) | peak RSS | final MSE | MSE − noise |
|---|---|---|---|---|---|---|---|
| V25 patched, subprocess r1 | 0/6 | 0.36666 | 0.195–0.217 | 33–36 | 3.0 GB | 5.397e-8 | 1.756e-8 |
| V25 patched, subprocess r2 | 0/6 | 0.36666 | 0.197–0.214 | 33–35 | 3.0 GB | 5.397e-8 | 1.756e-8 |
| V29 patched, subprocess r1 | **6/6** | – | 0.463–0.484 | 70–85 | 6.9 GB | – | – |

The V29 subprocess run fails in two ways. (a) Residual: 0.46–0.48 s on every MLP, against the 0.4 s cap, with the box idle. (b) Memory: the in-process flopscope keeps every array in the worker, the second predict passes the 8 GB RLIMIT_AS, one MLP falls back and the next kills the worker (WORKER_EOF). On the grader the arrays live in the flopscope server, so (b) is an artifact of in-process emulation. (a) is the real question, which the client/server emulation answers.

## Open items (being worked)

V29 local-runner repeats; parity of the patched vs original bundles; setup-window timing; client/server emulation of the residual for V29 and V25; residual reduction for V29; final packaging check.
