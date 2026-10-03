# Code for note IX

Everything is numpy/scipy. Networks are He-initialised ReLU MLPs (no biases) with Gaussian inputs.
`W_n{n}_L{L}_s{seed}.npy` is written by `truth.py`.

| file | what |
|---|---|
| `truth.py n L seed T` | Monte Carlo ground truth (per-layer activation means and variances) |
| `closure.py` | Hermite jets of ReLU and ReLU^2; Gaussian closure with the actual weights (mean + full covariance via Mehler) |
| `edgeworth.py` | one-point tree-level cumulants of `z = W g(y)`: D3, P3, T3, D4, PP4, PD4; Edgeworth mean shift |
| `twopoint.py` | two-point tree-level cumulants K21 = k(y_k,y_k,y_l), K22 = k(y_k,y_k,y_l,y_l) |
| `corrected.py` | `shift_coeffs` (jets a3,a4,a6,b3,b4,b6 used by the Edgeworth step) |
| `corrected3.py` | `tree_closure`: tree-level propagation (TLP) with linear-response transport (`rmax`), Hermite order `J`, optional gap-blob channel G1 |
| `blob.py` | gap-attached blob channel G1 (a 4-point blob of the previous layer absorbed by one gap) |
| `lite.py` | annealed O(n^2)-per-layer corrections (D3, coherent fourth-cumulant terms, annealed responses); no gain at depth |
| `diagsrc.py` | diagonal-source TLP (single-neuron sources only, ~5 n^3 products per layer per step) |
| `residue.py` | **closure + gain residue**: weights-only O(n^2) gain injection (Prop. 7.4), transported with factor `tau`; m <- (1 - gamma/8) m |
| `gain.py` | oracle scale coefficient of the closure error versus measured gain variance gamma/8, layer by layer |
| `inject.py` | per-layer gain injection: formula versus MC from an exactly Gaussian source |
| `tau_scan.py`, `shape_test.py` | residue correction across networks; uniform-scale versus gap-shaped correction |
| `calib_mc.py`, `calib_mc2.py` | calibrating the residue by Monte Carlo instead (noisy), with/without the radial trick |
| `bias_ablate.py`, `compare_diag.py`, `compare_lite.py` | which parts of TLP remove the scale error; cheap variants |
| `run_save.py`, `est1024.py`, `eval1024.py` | saved estimates (widths 256/512/1024) and the width-1024 evaluation |
| `compare5.py n L seed '{cfg}'` | layer-by-layer MSE of closure and TLP variants against ground truth |
| `test_layer2.py` | one-step exactness at layer 2 and term ablation |
| `gauss_source.py` | tree cumulants against MC when the source layer is exactly Gaussian (layers 2, 6) |
| `cycles.py` | triangle (first cycle) versus tree diagrams by depth; row-norm contraction versus 2<P^2> |
| `zstats.py`, `test_blob.py` | pre-activation cumulants by MC; residual of trees + linear response versus G1 |

`outputs/` holds the printed results quoted in the note.

Typical use (run from this directory):

```
python3 truth.py 256 8 0 1.6e7
python3 compare5.py 256 8 0 '{"J1r4": dict(rmax=4, J=1)}'
```
