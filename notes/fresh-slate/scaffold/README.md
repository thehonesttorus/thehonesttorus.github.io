# Fresh-slate scaffold: grader-safe template, porting kit, residual budget

*Status: first usable version (2026-10-01). Template verified end to end; the cost tables in KIT.md and
RESIDUAL.md are being filled in by `probe_kit.py` / `probe_residual.py`.*

Design-agnostic infrastructure for the fresh-slate design streams (notes/fresh-slate/BRIEF.md). Nothing here
encodes an estimation method: the placeholder `design()` in the template is the Gaussian covariance
closure, a baseline that only exercises the frame.

| file | what |
|---|---|
| `estimator.py` | the template: shape gate, guarded fallback, accountant, pools, float32 discipline, warm-ups; replace `Estimator.design()` |
| `kit.py` | numpy → flopscope porting kit: products, same-object Grams, sandwiches, Gaussian elementwise, linalg, gathers, pools, a lower-op Strassen–Winograd engine with level choice |
| `KIT.md` | measured unit and call prices of the kit's primitives at n = 1024 and batched small sizes |
| `RESIDUAL.md` | the residual-time budget model: seconds per flopscope call (idle box, 2 threads) |
| `verify.sh`, `run_summary.py` | end-to-end check with `whest run` under the graded caps on a baked 1024×16 set, the 256×32 smoke shape and 512×16 |

## Using the template

1. Copy `estimator.py` into `notes/fresh-slate/designs/<key>/`. Implement `design(self, mlp, budget)`: return the
   per-layer means, shape (L, n), on the suite shape 1024 × 16. Everything else is handled by the frame.
2. Develop with `WHEST_STRICT=1` (exceptions and non-finite outputs raise instead of silently taking the safe
   path) and `WHEST_ACCT=1` (per-phase accounting). Wrap blocks in `with ACCT.phase("name"):` and print
   `ACCT.report()` after a predict: calls, units, residual ms and wall ms per phase. Ship with both off.
3. Check: `./verify.sh path/to/estimator.py [n_mlps]`. It bakes (once, under /tmp/claude-0/scaffold-data) a
   1024 × 16 set with N = 20,000 (the truth noise is ~4e-6: this checks caps, failures, cost and residual, not
   accuracy), the 256 × 32 smoke MLP and one 512 × 16 MLP, then runs `whest run --runner subprocess
   --max-threads 2` with the caps (120 s wall, 0.4 s residual, 5 s setup) and prints one line per set.

Verified with the placeholder design (2026-10-01, this box, 4 vCPU, nothing else heavy running):

| set | fails | C/B | residual per MLP | wall |
|---|---|---|---|---|
| 256 × 32 smoke (safe path) | 0/1 | 0.0000 | 0.023 s | 0.2 s |
| 512 × 16 (safe path) | 0/1 | 0.0000 | 0.014 s | 0.1 s |
| 1024 × 16 (design, 31 u) | 0/2 | 0.0303 | 0.027, 0.021 s | 0.7 s |

## The rules the frame enforces, and the traps it avoids

- **Shape gate.** Only (1024, 16) runs `design()`. Every other shape, including the grader's 256 × 32 smoke test,
  goes to `_safe_predict` (float64 diagonal mean field, ~6 n² FLOPs per layer, any width and depth). A 16-row table
  indexed by layer on a depth-32 smoke MLP has already failed a public submission.
- **Guarded fallback.** Any exception except flopscope budget/time exhaustion, and any non-finite or mis-shaped
  output, returns the safe path instead of a zeroed MLP (multiplier 1.0). On the suite shape a fallback is a large
  raw loss, so develop under `WHEST_STRICT=1`.
- **No numpy** in the estimator (the grader image has flopscope and the standard library only). Use `fnp`.
- **No `x.shape = ...`** (flopscope #267: free in-process, keeps a stale symmetry tag, raises on the client).
- **No stale symmetry tags.** Only tags created deliberately and legitimately: `flops.as_symmetric` on a truly
  symmetric array, or the aliased same-object Gram `einsum('ji,jk->ik', X, X)`. Never benefit from tags left by
  shuffles, `choose`, `.shape` (#264/#265/#267); under the 16 Sep fair-accounting text that is disqualifiable,
  even after grading.
- **No tagged 3-operand einsum on indefinite inputs.** `einsum('ij,ia,jb->ab', C_tagged, W, W)` checks the symmetry
  of its float32 output with `allclose(atol=1e-6)` and raises `SymmetryError` on indefinite O(1) inputs.
- **Pools are `fnp.empty` + one `copyto(b, 0.0)`**, never `fnp.zeros/ones/eye` (those return symmetry-tagged
  arrays; writing a plain result into a tagged `out=` buffer raises).
- **Do not return a pooled buffer** from `predict` (the next predict overwrites it): return a copy.
- **float32 by default.** Weights can arrive float64: cast once per layer (`w.astype(f32)`); a float64 operand bills
  the whole op 2×. `flops.stats.norm.*` always bills float64: cast its result back.
- **gc disabled inside predict** (generation-2 collections of 35–110 ms land in the residual).
- **setup() ≤ 5 s in total, every worker spawn.** Use it only for tiny warm-up calls, one per op signature
  (library init otherwise costs 1–30 ms inside the first predict).
- **Residual is per call, not per FLOP.** ~0.02–0.04 ms per flopscope call whatever its size; a fresh (n, n)
  result ~0.1 ms; `out=` into pooled buffers is the cheap form. `reshape`, `transpose`, `swapaxes`, `moveaxis`,
  `diagonal` are logged calls (reshape is billed 1 FLOP per element): build such views once per pooled buffer.
  On the grader (client/server) even a basic slice is a round trip: cache slices of persistent buffers too.
- **Memory.** 8 GB; on the grader arrays live in the flopscope server, in local emulation in the worker.

## Status of the other deliverables

- `kit.py`: written; the Strassen engine is the one measured in notes/streams/est-cost (same products as
  504aldo's V29, 79 calls per L5 family instead of 113).
- `KIT.md`, `RESIDUAL.md`: being measured (next push).
