# Residual-time budget model (price a design in calls as well as units)

*Measured 2026-10-01 by `probe_residual.py` on this container (4 vCPU, idle apart from the probe, 2 BLAS threads):
in-process flopscope 0.12.1 (`results/residual_inprocess.json`) and the grader's transport, flopscope-client 0.12.1
against a flopscope-server 0.12.1 process (`results/residual_clientserver.json`). 2,000 calls per kind, best of 3.*

The grader caps the Python-side **residual** (wall − backend − flopscope overhead) at **0.4 s per MLP**; exceeding it
zeroes the MLP. Residual is paid per flopscope call, almost independent of the array size.

## Per-call residual (ms)

| call kind | in-process | client/server | client/server wall per call |
|---|---|---|---|
| vector add, out= | 0.022 | 0.014 | 1.19 |
| vector add, fresh result | 0.024 | 0.015 | 1.27 |
| small matmul (64×64), out= | 0.027 | 0.014 | 1.15 |
| (1024, 1024) add, out= | 0.031 | 0.021 | 2.90 |
| (1024, 1024) add, fresh result | 0.036 | 0.025 | 3.16 |
| add on strided slab views, views prebuilt | 0.024 | 0.015 | 1.38 |
| the same, slicing inside the call (3 slices) | 0.031 | **0.044** | 2.60 |
| a basic slice alone (no op) | 0.001 | **0.0125** | 0.66 |
| reshape (logged view) | 0.020 | 0.013 | 1.04 |
| vector norm.cdf + astype (2 calls) | 0.057 | 0.025 | 2.08 |

## The model

- **In-process (local `whest run`, the stricter case on this box):** residual ≈ 0.025 ms × calls (0.02–0.036 by
  size), + the Python logic between calls, + ~0.1 ms per fresh (n, n) result in a real chain. Whole-chain measurement:
  504aldo's V29 makes 13,121 calls and measured 0.46–0.52 s idle here (0.035–0.04 ms per call including its glue);
  its 9,869-call rewrite (notes/streams/est-cost) is proportionally lower.
- **Client/server (the grader's transport):** residual ≈ 0.014 ms per op + **0.0125 ms per basic slice** (every
  slice is a round trip) + Python logic. Wall ≈ 1–3 ms per call (round-trip latency): 10,000 calls are ≈ 10–30 s of
  the 120 s wall cap before any compute.
- **Budget for a 2× margin (0.2 s):** ≈ 6,000–8,000 calls per MLP in-process with prebuilt views, ≈ 10,000 ops
  client/server if slices are cached; count slices as calls on the grader.
- **Rules of thumb:** write results with `out=` into pooled buffers; build views of persistent buffers once
  (slices, reshapes, transposes, diagonals); batch many same-shape products into one call (no FLOP discount, one
  call); a Strassen family costs ~15 calls per level whatever its batch (KIT.md), so use fewer, fatter families; keep
  n-vector arithmetic in few calls (stack vectors into (k, n) arrays and operate once); disable gc in predict.
- **Setup window.** `setup()` itself took 0.085 s for the template, but interpreter start + imports + load measured
  4.8 s in the client/server emulation on this box (the 5 s setup cap's scope on the grader is the setup call; the
  submission stream measured 2.4–3.5 s windows for V25/V29). Keep setup() to tiny warm-ups.

## Template under both transports (2 MLPs, 1024 × 16, placeholder design, 31 u)

| run | fails | C/B | residual per MLP | wall max |
|---|---|---|---|---|
| `whest run --runner subprocess --max-threads 2` | 0/2 | 0.0303 | 0.027, 0.021 s | 0.7 s |
| client/server emulation (`notes/streams/submission/scripts/grader_emul.py`, 2 server threads) | 0/2 | 0.0303 | 0.013, 0.014 s | 1.3 s |

Emulator setup note: the client venv must hold **flopscope-client only**. Installing whestbench into it pulls the
full flopscope and silently overwrites the client package (everything then runs in-process); fix with
`pip uninstall -y flopscope flopscope-client && pip install --no-deps flopscope-client==0.12.1`.
