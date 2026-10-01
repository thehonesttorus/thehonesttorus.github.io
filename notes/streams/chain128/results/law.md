## (1) injected D21 noise at every layer, base fit:atlas

| eps | dMSE final mlp0 | dMSE final mlp1 | dMSE final mlp2 | dMSE final mlp3 | k = dMSE/eps^2 (mean) | dMSE mean-over-layers / eps^2 |
|---|---|---|---|---|---|---|
| 0.02 | 5.26e-08 | 1.91e-06 | 1.08e-06 | 6.82e-07 | 2.33e-03 | 4.45e-04 |
| 0.05 | 8.05e-07 | 9.69e-06 | 7.74e-06 | 2.95e-06 | 2.12e-03 | 4.37e-04 |
| 0.1 | 3.93e-06 | 3.65e-05 | 3.22e-05 | 1.01e-05 | 2.07e-03 | 4.36e-04 |
| 0.2 | 1.83e-05 | 1.51e-04 | 1.34e-04 | 3.79e-05 | 2.13e-03 | 4.45e-04 |

## (2) teacher forcing (kappa3(a_l) replaced by the atlas's)

| variant | final mlp0 | final mlp1 | final mlp2 | final mlp3 | ratio to base (mean) |
|---|---|---|---|---|---|
| F1@fit:atlas | 1.78e-06 | 8.04e-06 | 9.36e-06 | 5.27e-06 | 0.98 |
| F4@fit:atlas | 1.94e-06 | 8.03e-06 | 8.68e-06 | 4.69e-06 | 0.95 |
| F8@fit:atlas | 1.49e-06 | 7.68e-06 | 9.44e-06 | 3.87e-06 | 0.87 |
| F12@fit:atlas | 1.55e-06 | 5.12e-06 | 7.00e-06 | 4.47e-06 | 0.76 |
| F14@fit:atlas | 1.81e-06 | 5.66e-06 | 7.13e-06 | 4.92e-06 | 0.83 |
| FD4@fit:atlas | 1.95e-06 | 7.97e-06 | 8.73e-06 | 4.59e-06 | 0.95 |
| FD8@fit:atlas | 1.86e-06 | 7.80e-06 | 9.33e-06 | 4.11e-06 | 0.93 |
| FD12@fit:atlas | 2.02e-06 | 8.49e-06 | 8.33e-06 | 4.67e-06 | 0.97 |
| fit:atlas (base) | 1.97e-06 | 8.01e-06 | 9.35e-06 | 5.21e-06 | 1 |

## (3) across variants (teacher-forced kappa4): mean eps(D21) over layers 1-15 vs final MSE

| variant | mlp | mean eps D21 (raw) | mean eps D21 (noise-corr.) | final MSE |
|---|---|---|---|---|
| engine:atlas:1 | 0 | 0.043 | 0.036 | 1.07e-06 |
| engine:atlas_rank4:1 | 0 | 0.073 | 0.069 | 1.28e-06 |
| fit:atlas:1 | 0 | 0.051 | 0.044 | 1.41e-06 |
| engine:atlas_rank1:1 | 0 | 0.087 | 0.083 | 1.44e-06 |
| engine:atlas_reg211:1 | 0 | 0.100 | 0.097 | 1.88e-06 |
| engine:atlas | 0 | 0.047 | 0.041 | 1.89e-06 |
| fit:atlas | 0 | 0.050 | 0.044 | 1.97e-06 |
| closure:atlas:1 | 0 | 0.100 | 0.097 | 2.01e-06 |
| engine:atlas_rank32 | 0 | 0.060 | 0.055 | 2.13e-06 |
| engine:atlas_rank16 | 0 | 0.067 | 0.063 | 2.28e-06 |
| engine:atlas_rank4 | 0 | 0.078 | 0.074 | 2.46e-06 |
| engine:atlas_rank1 | 0 | 0.090 | 0.087 | 3.17e-06 |
| closure:atlas_reg211:1 | 0 | 0.147 | 0.145 | 4.82e-06 |
| engine:atlas_reg211 | 0 | 0.101 | 0.098 | 4.93e-06 |
| closure:atlas | 0 | 0.097 | 0.094 | 5.41e-06 |
| engine:atlas_zero211:1 | 0 | 0.171 | 0.169 | 6.54e-06 |
| engine:atlas_zero211 | 0 | 0.168 | 0.166 | 8.50e-06 |
| wick:atlas | 0 | 0.162 | 0.160 | 8.56e-06 |
| wick:atlas:1 | 0 | 0.157 | 0.155 | 1.06e-05 |
| closure:atlas_zero211:1 | 0 | 0.232 | 0.230 | 1.07e-05 |
| closure:atlas_reg211 | 0 | 0.144 | 0.142 | 1.19e-05 |
| closure:atlas_zero211 | 0 | 0.228 | 0.226 | 1.69e-05 |
| none:atlas:1 | 0 | 0.533 | 0.533 | 5.52e-05 |
| none:atlas | 0 | 0.534 | 0.534 | 1.01e-04 |
| engine:atlas:1 | 1 | 0.044 | nan | 1.64e-06 |
| fit:atlas:1 | 1 | 0.042 | nan | 3.04e-06 |
| engine:atlas_rank4:1 | 1 | 0.067 | nan | 3.47e-06 |
| engine:atlas_rank1:1 | 1 | 0.077 | nan | 4.11e-06 |
| engine:atlas_reg211:1 | 1 | 0.089 | nan | 6.82e-06 |
| closure:atlas:1 | 1 | 0.077 | nan | 6.93e-06 |
| fit:atlas | 1 | 0.053 | nan | 8.01e-06 |
| engine:atlas | 1 | 0.063 | nan | 8.57e-06 |
| engine:atlas_rank4 | 1 | 0.080 | nan | 8.77e-06 |
| engine:atlas_rank16 | 1 | 0.073 | nan | 8.80e-06 |
| engine:atlas_rank32 | 1 | 0.068 | nan | 9.59e-06 |
| engine:atlas_rank1 | 1 | 0.090 | nan | 1.03e-05 |
| closure:atlas_reg211:1 | 1 | 0.146 | nan | 1.10e-05 |
| engine:atlas_reg211 | 1 | 0.088 | nan | 1.14e-05 |
| wick:atlas | 1 | 0.126 | nan | 1.93e-05 |
| engine:atlas_zero211:1 | 1 | 0.197 | nan | 2.10e-05 |
| wick:atlas:1 | 1 | 0.127 | nan | 2.24e-05 |
| closure:atlas | 1 | 0.079 | nan | 2.37e-05 |
| closure:atlas_reg211 | 1 | 0.140 | nan | 2.40e-05 |
| engine:atlas_zero211 | 1 | 0.191 | nan | 3.35e-05 |
| closure:atlas_zero211:1 | 1 | 0.271 | nan | 3.58e-05 |
| closure:atlas_zero211 | 1 | 0.264 | nan | 5.68e-05 |
| none:atlas:1 | 1 | 0.592 | nan | 1.18e-04 |
| none:atlas | 1 | 0.591 | nan | 1.73e-04 |
| engine:atlas:1 | 2 | 0.043 | nan | 4.20e-06 |
| fit:atlas:1 | 2 | 0.044 | nan | 8.40e-06 |
| closure:atlas:1 | 2 | 0.106 | nan | 8.45e-06 |
| engine:atlas_rank4:1 | 2 | 0.067 | nan | 8.73e-06 |
| fit:atlas | 2 | 0.050 | nan | 9.35e-06 |
| engine:atlas | 2 | 0.053 | nan | 1.06e-05 |
| engine:atlas_rank1:1 | 2 | 0.085 | nan | 1.10e-05 |
| engine:atlas_reg211:1 | 2 | 0.093 | nan | 1.16e-05 |
| engine:atlas_rank32 | 2 | 0.061 | nan | 1.22e-05 |
| engine:atlas_rank4 | 2 | 0.073 | nan | 1.34e-05 |
| engine:atlas_rank16 | 2 | 0.065 | nan | 1.37e-05 |
| engine:atlas_rank1 | 2 | 0.088 | nan | 1.54e-05 |
| closure:atlas_reg211:1 | 2 | 0.145 | nan | 1.56e-05 |
| engine:atlas_reg211 | 2 | 0.094 | nan | 1.78e-05 |
| closure:atlas | 2 | 0.106 | nan | 3.44e-05 |
| closure:atlas_reg211 | 2 | 0.143 | nan | 3.97e-05 |
| engine:atlas_zero211:1 | 2 | 0.166 | nan | 4.02e-05 |
| engine:atlas_zero211 | 2 | 0.163 | nan | 4.85e-05 |
| closure:atlas_zero211:1 | 2 | 0.234 | nan | 5.44e-05 |
| wick:atlas | 2 | 0.172 | nan | 6.12e-05 |
| wick:atlas:1 | 2 | 0.160 | nan | 8.04e-05 |
| closure:atlas_zero211 | 2 | 0.230 | nan | 8.95e-05 |
| none:atlas:1 | 2 | 0.580 | nan | 1.30e-04 |
| none:atlas | 2 | 0.580 | nan | 2.35e-04 |
| closure:atlas:1 | 3 | 0.102 | nan | 2.97e-06 |
| engine:atlas:1 | 3 | 0.049 | nan | 3.20e-06 |
| fit:atlas:1 | 3 | 0.058 | nan | 3.77e-06 |
| engine:atlas | 3 | 0.051 | nan | 4.43e-06 |
| engine:atlas_rank4:1 | 3 | 0.080 | nan | 4.88e-06 |
| fit:atlas | 3 | 0.057 | nan | 5.21e-06 |
| engine:atlas_rank32 | 3 | 0.063 | nan | 5.43e-06 |
| engine:atlas_rank16 | 3 | 0.070 | nan | 6.94e-06 |
| engine:atlas_rank1:1 | 3 | 0.096 | nan | 8.41e-06 |
| engine:atlas_reg211:1 | 3 | 0.115 | nan | 9.38e-06 |
| closure:atlas | 3 | 0.097 | nan | 1.05e-05 |
| engine:atlas_rank4 | 3 | 0.081 | nan | 1.05e-05 |
| closure:atlas_reg211:1 | 3 | 0.160 | nan | 1.38e-05 |
| engine:atlas_rank1 | 3 | 0.098 | nan | 1.52e-05 |
| engine:atlas_zero211:1 | 3 | 0.183 | nan | 1.56e-05 |
| engine:atlas_reg211 | 3 | 0.114 | nan | 2.05e-05 |
| wick:atlas | 3 | 0.191 | nan | 2.05e-05 |
| closure:atlas_zero211:1 | 3 | 0.242 | nan | 2.16e-05 |
| engine:atlas_zero211 | 3 | 0.180 | nan | 2.51e-05 |
| wick:atlas:1 | 3 | 0.188 | nan | 2.54e-05 |
| closure:atlas_reg211 | 3 | 0.156 | nan | 3.16e-05 |
| closure:atlas_zero211 | 3 | 0.238 | nan | 3.86e-05 |
| none:atlas:1 | 3 | 0.540 | nan | 1.10e-04 |
| none:atlas | 3 | 0.540 | nan | 1.63e-04 |

## (4) structured law per MLP (Edgeworth order 2): final MSE = a + k eps^2 over the kappa4-teacher-forced variants

| mlp | a | k | corr(MSE, eps^2) | K=2 MSE (variant A) | k / K2 MSE |
|---|---|---|---|---|---|
| 0 | 1.04e-06 | 3.48e-04 | 0.998 | 9.41e-05 | 3.70 |
| 1 | 1.09e-05 | 4.73e-04 | 0.991 | 4.74e-04 | 1.00 |
| 2 | 2.01e-05 | 6.63e-04 | 0.971 | 1.96e-04 | 3.37 |
| 3 | 7.04e-06 | 5.39e-04 | 0.994 | 2.34e-04 | 2.30 |

## (4) structured law per MLP (Edgeworth order 1): final MSE = a + k eps^2 over the kappa4-teacher-forced variants

| mlp | a | k | corr(MSE, eps^2) | K=2 MSE (variant A) | k / K2 MSE |
|---|---|---|---|---|---|
| 0 | 1.08e-06 | 1.91e-04 | 0.995 | 9.41e-05 | 2.03 |
| 1 | 5.64e-06 | 3.26e-04 | 0.989 | 4.74e-04 | 0.69 |
| 2 | 1.79e-05 | 3.51e-04 | 0.857 | 1.96e-04 | 1.79 |
| 3 | 3.61e-06 | 3.67e-04 | 0.994 | 2.34e-04 | 1.57 |

## (5) all recorded variants on the atlas MLPs: final MSE per MLP | geometric mean | mean-over-layers (mean)

| variant | mlp0 | mlp1 | mlp2 | mlp3 | geo-mean final | mean over layers |
|---|---|---|---|---|---|---|
| engine:atlas:1 | 1.07e-06 | 1.64e-06 | 4.20e-06 | 3.20e-06 | 2.20e-06 (4) | 9.38e-07 |
| fit:atlas:1 | 1.41e-06 | 3.04e-06 | 8.40e-06 | 3.77e-06 | 3.42e-06 (4) | 1.30e-06 |
| engine:atlas_rank4:1 | 1.28e-06 | 3.47e-06 | 8.73e-06 | 4.88e-06 | 3.71e-06 (4) | 1.60e-06 |
| closure:atlas:1 | 2.01e-06 | 6.93e-06 | 8.45e-06 | 2.97e-06 | 4.32e-06 (4) | 1.43e-06 |
| fit:atlas | 1.97e-06 | 8.01e-06 | 9.35e-06 | 5.21e-06 | 5.26e-06 (4) | 1.59e-06 |
| engine:atlas_rank1:1 | 1.44e-06 | 4.11e-06 | 1.10e-05 | 8.41e-06 | 4.84e-06 (4) | 1.97e-06 |
| engine:atlas | 1.89e-06 | 8.57e-06 | 1.06e-05 | 4.43e-06 | 5.25e-06 (4) | 1.50e-06 |
| engine:atlas_rank32 | 2.13e-06 | 9.59e-06 | 1.22e-05 | 5.43e-06 | 6.06e-06 (4) | 1.84e-06 |
| engine:atlas_reg211:1 | 1.88e-06 | 6.82e-06 | 1.16e-05 | 9.38e-06 | 6.11e-06 (4) | 2.27e-06 |
| engine:atlas_rank16 | 2.28e-06 | 8.80e-06 | 1.37e-05 | 6.94e-06 | 6.60e-06 (4) | 1.96e-06 |
| engine:atlas_rank4 | 2.46e-06 | 8.77e-06 | 1.34e-05 | 1.05e-05 | 7.43e-06 (4) | 2.19e-06 |
| closure:dense | 1.02e-05 | — | — | — | 1.02e-05 (1) | 2.90e-06 |
| engine:atlas_rank1 | 3.17e-06 | 1.03e-05 | 1.54e-05 | 1.52e-05 | 9.35e-06 (4) | 2.60e-06 |
| closure:atlas_reg211:1 | 4.82e-06 | 1.10e-05 | 1.56e-05 | 1.38e-05 | 1.03e-05 (4) | 2.92e-06 |
| engine:atlas_reg211 | 4.93e-06 | 1.14e-05 | 1.78e-05 | 2.05e-05 | 1.19e-05 (4) | 3.34e-06 |
| engine:dense | 8.12e-06 | 2.03e-05 | — | — | 1.28e-05 (2) | 2.90e-06 |
| closure:atlas | 5.41e-06 | 2.37e-05 | 3.44e-05 | 1.05e-05 | 1.47e-05 (4) | 3.54e-06 |
| engine:atlas_zero211:1 | 6.54e-06 | 2.10e-05 | 4.02e-05 | 1.56e-05 | 1.71e-05 (4) | 5.18e-06 |
| engine:zero:1 | 4.45e-06 | 4.62e-05 | 2.54e-05 | 3.09e-05 | 2.01e-05 (4) | 1.31e-05 |
| closure:atlas_reg211 | 1.19e-05 | 2.40e-05 | 3.97e-05 | 3.16e-05 | 2.45e-05 (4) | 5.51e-06 |
| wick:atlas | 8.56e-06 | 1.93e-05 | 6.12e-05 | 2.05e-05 | 2.13e-05 (4) | 6.18e-06 |
| closure:zero:1 | 4.66e-06 | 4.62e-05 | 2.46e-05 | 3.45e-05 | 2.07e-05 (4) | 1.30e-05 |
| engine:atlas_zero211 | 8.50e-06 | 3.35e-05 | 4.85e-05 | 2.51e-05 | 2.43e-05 (4) | 6.48e-06 |
| closure:atlas_zero211:1 | 1.07e-05 | 3.58e-05 | 5.44e-05 | 2.16e-05 | 2.59e-05 (4) | 6.68e-06 |
| engine:zero | 6.04e-06 | 4.83e-05 | 4.12e-05 | 4.18e-05 | 2.66e-05 (4) | 1.54e-05 |
| wick:atlas:1 | 1.06e-05 | 2.24e-05 | 8.04e-05 | 2.54e-05 | 2.64e-05 (4) | 7.35e-06 |
| closure:zero | 6.47e-06 | 4.87e-05 | 4.07e-05 | 4.71e-05 | 2.79e-05 (4) | 1.54e-05 |
| engine:dense2 | 2.52e-06 | 8.29e-05 | 4.71e-05 | 1.41e-05 | 1.93e-05 (4) | 5.57e-06 |
| closure:atlas_zero211 | 1.69e-05 | 5.68e-05 | 8.95e-05 | 3.86e-05 | 4.27e-05 (4) | 9.90e-06 |
| closure:mem | 7.19e-05 | 9.87e-05 | 1.14e-04 | 1.10e-04 | 9.71e-05 (4) | 2.47e-05 |
| none:atlas:1 | 5.52e-05 | 1.18e-04 | 1.30e-04 | 1.10e-04 | 9.83e-05 (4) | 4.40e-05 |
| none:atlas | 1.01e-04 | 1.73e-04 | 2.35e-04 | 1.63e-04 | 1.61e-04 (4) | 5.96e-05 |
| engine:mem:1 | 9.75e-05 | 1.93e-04 | 2.82e-04 | 2.16e-04 | 1.84e-04 (4) | 3.94e-05 |
| k2:zero | 9.41e-05 | 4.74e-04 | 1.96e-04 | 2.34e-04 | 2.13e-04 (4) | 1.23e-04 |
