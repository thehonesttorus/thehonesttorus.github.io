# Stream est-cost: cut the public chain's FLOP bill and flopscope call count at unchanged numbers

*Status: STOPPED at 17:45 UTC on 2026-10-01 on the user's instruction (the competition system is to be designed from a
fresh slate; adapting the public chain is out). Everything below is measured unless marked.*

**Stop note.** Measured before stopping: V29's steady ledger reproduces exactly (260.06 units, 13,121 flopscope calls);
the five Strassen families hold 7,170 calls and the rest is ~400 small elementwise calls per layer. A lower-op Strassen
engine (`sw2.py`: 6-D quadrant grids, two blocks per op, 15 calls per level instead of 22, fused leaf 19 instead of 25)
cuts an L5 family from 113 to 79 calls (hub 127 → 87) at the same price and error. Bundle v29d (engine + the d=2 term
program on 12 matrix parts instead of 52 materialized products + the cumulant conversion as one einsum + dead K4 riders
removed) bills **258.89 units (V29 260.06) in 9,869 calls (V29 13,123, −25 %)** at steady state, with final-layer outputs
within float32 noise of V29 (rms rel 8.0e-7, max 1.0e-6). v29e (C_pre copies paired) runs 9,827 calls. Not done:
idle-box / client-server residual, robustness checks, further elementwise fusions, FLOP levers. Measured exchange rates
(estimated from the op table): keeping C_pre's block family vs going dense ≈ 140 calls per unit; old-tier families
dense ≈ 65 calls per unit; more Strassen on the join products ≈ 160–190 calls per unit. Line-level call attribution of
v29e is in `results/prof_v29e_lines.json`. Note from the submission stream: on the grader every basic slice is also a
client/server round trip, so slice memoization matters as much as op count.

## Question

Can the public 504aldo V29 chain (MIT) be billed fewer FLOPs and, above all, run with far fewer flopscope calls
(its residual, 0.46-0.52 s per MLP in the local 2-thread runner, fails the 0.4 s cap), without changing its outputs
beyond float32 noise and without any accounting trick (fair-accounting rule of 16 Sep: no packing, no stale tags)?

## Method

- `prof.py`: in-process flopscope 0.12.1 meter (whestbench 0.16.1 venv), random He MLP 1024 x 16 (cost is
  data-independent), 2-3 predicts per process (the first is V29's staged L4 call), steady state = last call.
  Units = 2^31 FLOPs, B = 1024 units. Per family / layer: units, calls, and residual attributed to the op that
  follows each wall gap. Saves the predictions of every call for parity.
- Parity: max and rms relative difference of the final-layer output against the patched V29 bundle on the same MLP.
- Residual numbers in this section were taken on a loaded box (a dataset bake was running): only call counts and
  units are comparable across rows here; idle-box residuals follow in the validation section.

## Results so far

### V29 ledger reproduced (steady call, namespace-tagged V29)

260.06 units, 13,121 calls (exactly the published ledger). Calls by family: hub 1,811, young_transport 1,753, birth
1,485, cpre 1,418, old_legs 1,105, shared 1,083, nonlin 990, pk2k 705, fb 533, wick 448, elem 424, regen 313, ...
The five Strassen families hold 7,170 calls (≈ 60 % of the residual gaps); the rest is ~400 small elementwise
calls per layer.

### Lever 1: a lower-op Strassen engine (`sw2.py`, bundled in `bundles/v29c`)

The same Strassen-Winograd products and the same pools as V29, but every level works on 6-D quadrant grids
(b, P, 2, h, 2, q) and does two quadrant blocks per op wherever they form a strided view (column pairs, the
diagonal, the anti-diagonals; a single block broadcasts). 15 calls per level instead of 22; the fused leaf does the
seven products as four matmuls of pairs (19 calls instead of 25). Grids of pooled buffers are built once.

| one family, n = 1024, b = 4 | V29 calls | new calls | units / product V29 → new | rel. error vs f64 |
|---|---|---|---|---|
| plain mm, L5 | 113 | 79 | 0.5441 → 0.5452 | 8.1e-6 → 8.2e-6 |
| hub (sum_k X_k Y_k^T), L5 | 127 | 87 | 0.5465 → 0.5468 | 5.6e-6 → 5.6e-6 |
| plain mm, L3 (n = 256) | 67 | 46 | 0.7070 → 0.7070 | 1.7e-6 → 1.7e-6 |

Whole chain (v29c = patched V29 + engine only), He MLP:

| | units call 1 / steady | calls call 1 / steady | final-layer rel. diff vs V29 (rms / max) |
|---|---|---|---|
| V29 (patched bundle) | 273.47 / 260.06 | 12,907 / 13,123 | – |
| v29c | 274.57 / 261.08 | 11,847 / 11,192 | 7.8e-7 / 1.1e-6 |

(The +1 unit of v29c's second call is the one-time reshape of the L5 pools allocated in that call; to be checked on a
third call.)

## Open items (being worked)

Elementwise fusions (term program on 12 matrix parts instead of 52 materialized products; cumulant conversion as one
einsum; dead K4 riders), C_pre form, FLOP levers and their call prices, idle-box residual, robustness checks, bundle.
