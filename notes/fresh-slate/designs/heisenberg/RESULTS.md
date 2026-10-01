# Heisenberg–Duhamel stream — results log

All numbers are final-layer MSE against Monte Carlo truth (truth noise stated; subtract where it matters). Depth 16 unless stated; He-init networks, x @ W convention.

## R0. Toy exactness check of the Duhamel identity (Theorem 1)

`toy_duhamel.py` (n = 8, L = 3, N = 2e7 per expectation): every term (ρ_l − ν_l)[g_l] is estimated by Monte Carlo with the exact network after layer l.

| quantity | value |
|---|---|
| rms(truth − reference) | 6.03e-2 |
| rms(truth − reference − Σ_l terms) | 2.1e-4 (MC standard error of the sum ≈ 1.8e-4) |
| term l = 2 (age 0): rms exact / rms κ₃-Stein pairing / corr | 2.63e-2 / 2.29e-2 / 0.972 |
| term l = 1 (needs one transport step): rms / pairing / corr | 5.12e-2 / 3.38e-2 / 0.875 |

The identity holds to Monte Carlo precision. The first-order Stein pairing captures the age-0 term well and the transported term less well, at a width (8) where ε is O(1).

## R1. First prototype, own Monte Carlo truth (N = 1e7), width 64, 2 random networks

| variant | net 0 final MSE | net 1 final MSE |
|---|---|---|
| reference (Gaussian closure) | 3.82e-3 | 1.15e-3 |
| HD, A = 0 (age-0 sources only) | 2.99e-3 | 9.37e-4 |
| HD, A = 1 | 1.87e-3 | 6.18e-4 |
| HD, A = ∞ (all ages, full κ₃ source series, decoupled-gate transport) | **2.16e-4** | **5.00e-5** |

Gain of 18× and 23× for A = ∞, and almost none for A ≤ 1: the content that matters at the readout is old, as the brief's facts predicted. Width 32 is non-perturbative (gain 4.5× on one net; the other went NaN through an indefinite corrected covariance).

## R2. Oracle diagnostic (`oracle.py`, width 64, N = 1e7): where is the residual?

The true per-layer cumulants of z_l (diagonal κ₃, (2,1) slice, diagonal κ₄) are measured by Monte Carlo and injected with the same Stein formulas.

| variant | net 0 | net 1 |
|---|---|---|
| reference | 3.82e-3 | 1.15e-3 |
| true κ₃ (D, S), first-order injection, own (m, C) chain | 1.71e-4 | 4.83e-4 |
| + true diagonal κ₄ and κ₃² (marginal only) | 1.97e-4 | 1.20e-3 |
| true (m, C) only, no injection | 4.15e-5 | 4.75e-5 |
| **true (m, C) + true κ₃, κ₄ injection** | **2.8e-6** | **4.3e-6** |

Readings. (i) HD's computed sources and transport are as good as the true κ₃ under the same injection (2.2e-4 vs 1.7e-4; 5.0e-5 vs 4.8e-4), so they are not the bottleneck at width 64. (ii) The bottleneck is **covariance drift**: the same injection with the true (m, C) is ≈ 1000× better than the reference. (iii) Second-order marginal terms fix the early layers (layer 2: 3.9e-6 → 1.0e-7) but not the deep ones; what is missing is second order on the **off-diagonal covariance** (joint κ₄ slices κ₄(ppqq), κ₄(pppq) and κ₃·κ₃ products). The readout variance wᵀ Cov(a) w sums n² entries, each with an O(ε²) error, coherently.
