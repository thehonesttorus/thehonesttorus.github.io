# Shared benchmark for the fresh-slate design streams

Stage Q (numpy prototypes, widths 64–512) and Stage P (flopscope at 1024 × 16) evaluation, per `../BRIEF.md` §4.

| file | what |
|---|---|
| `bench.py` | the sets: `list_sets()`, `load_set(name)` (seeds, width, depth, N, avg_variance, truth means (m, L, n) float64, noise = avg_variance/N), `weights(S, i)` regenerates the (L, n, n) float32 weights from the seed, bit-identical to `whest dataset bake` (verified at widths 64 and 1024) |
| `sets/` | one `<name>.json` (seeds, N, avg_variance) + `<name>.npz` (all_layer_means) per set; small files only |
| `eval_q.py` | Stage Q evaluator: `evaluate(predict, sets, units_1024)` → per-MLP final-layer MSE, raw (minus truth noise), all-layer and per-layer MSE, width fit raw ∝ n^-p extrapolated to 1024, projected adjusted = raw × max(0.1, units/1024); calibration baselines `--baseline gauss|mc` |
| `run_p_inproc.py` | Stage P quick runner: a whest `Estimator` in-process on a bench set (weights from seeds, metered under 2^41) |
| `run_p.py`, `eval_p.md` | Stage P grader-faithful runs through `whest run --runner subprocess` with the caps, and the robustness checklist |
| `RESULTS.md` | baseline numbers on every set |

Usage:

```python
import sys; sys.path.insert(0, "notes/fresh-slate/bench")
import eval_q
def predict(W):            # W: (L, n, n) float32, x @ W convention; return (L, n) post-ReLU means
    ...
eval_q.evaluate(predict, ["w64_d16", "w128_d16", "w256_d16"], units_1024=110)
```

CLI: `python eval_q.py --module my.py --func predict --sets w64_d16,w128_d16 --units 110 --json out.json`.

Sets (final; seeds consecutive from the first listed; w512/w1024 truth noise exceeds the best errors: compare paired):

| name | width × depth | MLPs | seeds | N | truth noise |
|---|---|---|---|---|---|
| w1024_d16 | 1024 × 16 | 6 | 7301001… | 2e+06 | ≤ 4.4e-08 |
| w128_d16 | 128 × 16 | 8 | 128001… | 6e+07 | ≤ 4.0e-09 |
| w256_d16 | 256 × 16 | 8 | 256001… | 2e+07 | ≤ 5.7e-09 |
| w256_d32 | 256 × 32 | 2 | 256101… | 1e+07 | ≤ 6.6e-09 |
| w512_d16 | 512 × 16 | 4 | 512001… | 2e+06 | ≤ 4.5e-08 |
| w64_d16 | 64 × 16 | 8 | 64001… | 1e+08 | ≤ 1.1e-09 |
