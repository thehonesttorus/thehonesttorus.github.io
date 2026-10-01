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

Sets (see `python bench.py`; 'q_' sets are the quick first bakes, replaced by higher-N bakes of the SAME seeds as they
finish — the networks never change, only the truth noise falls):

| name | width × depth | MLPs | seeds | N | truth noise |
|---|---|---|---|---|---|
| see `python bench.py` | | | | | |
