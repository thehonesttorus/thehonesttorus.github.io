# Baselines on the bench sets (calibration only, not designs)

gauss = Gaussian covariance closure with linearised cross-covariance (`eval_q.baseline_gauss`, float64);
mc = plain Monte Carlo with 6554 samples (the 0.1-floor budget at n = 1024: 0.1 B / (2 L n²)).
raw = final-layer MSE − truth noise; mean over the set's MLPs (± s.e. across MLPs).

| set | MLPs | N | truth noise | gauss raw | mc raw | gauss all-layer |
|---|---|---|---|---|---|---|
| w64_d16 | 8 | 1e+07 | 6.6e-09 | 5.106e-04 ± 8.5e-05 | 6.655e-06 ± 7.5e-07 | 5.095e-04 |
