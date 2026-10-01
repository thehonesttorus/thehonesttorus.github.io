# Baselines on the bench sets (calibration only, not designs)

gauss = Gaussian covariance closure with linearised cross-covariance (`eval_q.baseline_gauss`, float64);
mc = plain Monte Carlo with 6554 samples (the 0.1-floor budget at n = 1024: 0.1 B / (2 L n²)).
raw = final-layer MSE − truth noise; mean over the set's MLPs (± s.e. across MLPs).

| set | MLPs | N | truth noise | gauss raw | mc raw | gauss all-layer |
|---|---|---|---|---|---|---|
| w1024_d16 | 6 | 2e+06 | 3.6e-08 | 4.299e-06 ± 3.6e-07 | 1.086e-05 ± 8.5e-07 | 2.434e-06 |
| w128_d16 | 8 | 1e+07 | 1.3e-08 | 2.890e-04 ± 6.2e-05 | 1.241e-05 ± 3.0e-06 | 1.475e-04 |
| w256_d32 | 2 | 1e+07 | 5.0e-09 | 1.171e-04 ± 7.5e-05 | 6.638e-06 ± 4.4e-06 | 6.533e-05 |
| w512_d16 | 4 | 2e+06 | 3.8e-08 | 1.784e-05 ± 2.0e-06 | 1.014e-05 ± 3.7e-07 | 1.012e-05 |
| w64_d16 | 8 | 1e+07 | 6.6e-09 | 5.106e-04 ± 8.5e-05 | 6.655e-06 ± 7.5e-07 | 5.095e-04 |

gauss width fit (depth-16 sets): raw ∝ n^-1.67, extrapolated raw(1024) = 6.31e-06; measured at 1024: 4.30e-06

Caution: a two-width power-law fit from widths 64–128 is NOT a reliable 1024 projection. For the Gaussian
closure it predicts 5.2e-5 while the measured 1024 value is 4.3e-6 (12x lower; the small widths are in a different
regime). Fit on the widest sets you can afford (256, 512) and check against w1024_d16 when possible.
