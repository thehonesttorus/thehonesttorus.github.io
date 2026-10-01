# signings stream — results log

Machine: 4 shared cores. Truth: own Monte Carlo bake (`truth.py`, weights from `whestbench.generation.sample_mlp`),
float32 forward, float64 accumulation; raw MSE has the per-neuron truth noise var_j/N subtracted.

## T0 — exactness of one copula step (layer 1 is an exact Gaussian copula)

Per-sector check at width 12 (`t0_sectors.py`, 4e6 MC samples of a_1, exact joint cumulants vs the Mehler
generators of `copula.py`):

| sector of κ_k(z_2j) | generator rms | MC rms | diff rms | note |
|---|---|---|---|---|
| κ3 {2,1} | 0.2923 | 0.2925 | 0.0013 | exact (MC noise) |
| κ3 {1,1,1} paths | 0.0821 | 0.0772 | 0.0105 | MC includes triangles |
| κ4 {2,2} | 0.1639 | 0.1646 | 0.0036 | exact |
| κ4 {3,1} | 0.3821 | 0.3825 | 0.0039 | exact |
| κ4 {2,1,1} paths | 0.2301 | 0.2427 | 0.1848 | MC includes triangles |

Entry-level check of κ(a_a,a_a,a_b,a_c) (`t0_211.py`): the full three-vertex Mehler sum *including* the
triangle multigraphs reproduces MC (−0.0511 vs −0.0515, 0.0976 vs 0.0971, ...), the tree-only sum does not
(−0.0305): at width 12 (|ρ| ≈ 0.3) loops are not negligible; the generators are correct.

T1 (loop share vs width, `t0_exact.py`, κ3/κ4 of z_2 vs 2e7-sample MC, tree sectors only):

| n | κ3 rms err / rms κ3 | κ4 rms err / rms κ4 |
|---|---|---|
| 12 | 0.0099 / 0.94 | 0.184 / 1.38 |
| 24 | 0.0085 / 0.54 | 0.209 / 0.80 |
| 48 | 0.0043 / 0.29 | 0.050 / 0.45 |
