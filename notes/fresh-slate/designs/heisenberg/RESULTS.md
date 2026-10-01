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

## R3. Second-order oracle (`oracle2.py`, `oracle3.py`, width 64, N = 5e6): which second-order pieces matter

True cumulants of z_l injected with the full second-order Stein/Edgeworth formula on the whole covariance (`inject_full2`: κ₃ diag + (2,1) slice; κ₄ diag, κ₄(pppq), κ₄(ppqq); κ₃² terms), own (m, C) chain.

| injection | net 0 | net 1 |
|---|---|---|
| reference | 3.83e-3 | 1.15e-3 |
| first order (true κ₃) | 1.74e-4 | 4.78e-4 |
| **full second order (true κ₃, κ₄ slices)** | **2.25e-5** | **2.35e-5** |
| same, without joint κ₄ (diag κ₄ + κ₃² only) | 2.6e-4 | NaN |
| joint κ₄(ppqq) only | NaN | NaN |
| joint κ₄(pppq) only | 9.8e-5 | 6.4e-5 |
| both joint κ₄ slices, no off-diagonal κ₃² | 2.25e-5 | 3.1e-5 |
| first order with true (m, C) | 3.0e-5 | 1.1e-5 |
| full second order with true (m, C) | 2.8e-6 | 4.2e-6 |

Reading: the joint fourth-cumulant slices κ₄(pppq) and κ₄(ppqq) are necessary at width 64 and must be used *together* (either alone destabilises the covariance); the off-diagonal κ₃² terms are negligible. With them the own-chain injection gains 50–170× over the reference.

## R4. Stage Q, shared bench (first-order HD-∞: all ages, full κ₃ source series, first-order injection)

| set | MLPs | reference raw | HD-∞ raw (± s.e.) | gain |
|---|---|---|---|---|
| w64_d16 | 8 | 4.46e-4 | 9.75e-5 ± 1.8e-5 | 4.6× |
| w128_d16 | 8 | 2.82e-4 | 3.16e-5 ± 5.8e-6 | 8.9× |

Width exponent of HD-∞ first order: raw ∝ n^{-1.6} between 64 and 128 (reference: n^{-0.66} on this pair). Extrapolated with p = 1.6 to 1024: raw ≈ 1.1e-6, far from the 1.5e-8 bar. First order alone is not competitive, consistent with R2–R3.

## R5. Accuracy of the computed (age-0) joint κ₄ slices (`gen_k4_slices`, width 64, layer 7, vs MC with N = 4e6; MC-vs-MC noise 3–4e-4)

| slice | corr | rms true | rms error |
|---|---|---|---|
| κ₄(pppq) | 0.973 | 7.3e-3 | 2.2e-3 |
| κ₄(ppqq) | 0.909 | 5.4e-3 | 3.6e-3 |
| κ₄ diag | 0.969 | 1.5e-2 | 4.7e-3 |

Diagrams used: second-chaos chain (exact pair contractions), degree-3 star, exact single-site; adding the He₂ 4-cycle made no difference. The residual is truncation in ρ: at width 64 the mean |ρ| between pre-activations is 0.15 (layer 3) to 0.32 (layer 15). It falls as n^{-1/2}: 0.054–0.14 at width 512, so the diagram truncation improves with width.

## R6. Width 256 (own networks, seeds 500–502, own MC truth N = 2e6, truth noise 3.6–5.7e-8), first-order HD-∞ (`w256.py`, Kt = 3)

| net | reference | HD-∞ | gain |
|---|---|---|---|
| 0 | 7.11e-5 | 5.19e-6 | 13.7× |
| 1 | 6.71e-5 | 1.03e-5 | 6.5× |
| 2 | 3.04e-5 | 6.82e-6 | 4.5× |
| mean | 5.6e-5 | **7.4e-6 ± 1.2e-6** | 7.6× |

HD-∞ width exponent: p = 2.09 (128 → 256) and 1.86 (64 → 256, mixing bench and own networks). These are 3 networks, not the bench set.

## R7. Cheap realisation of first order (w64 bench, 8 MLPs)

| variant | raw |
|---|---|
| HD-∞, full κ₃ source series | 9.75e-5 |
| HD-∞, **diagonal + star diagrams only** (the n³-cost pull-back of DESIGN §2) | 1.00e-4 |
| + triangle | 1.00e-4 |
| ages A = 7 (full series) | 1.05e-4 |
| ages A = 3 | 1.54e-4 |

The cheap diagrams lose nothing. Truncating at 8 ages costs 8 %; at 4 ages, 57 %.

## R8. Second-order attempts (what failed, and why)

| variant (w64 bench unless stated) | raw | reading |
|---|---|---|
| HD-2: first-order HD-∞ + computed age-0 joint κ₄ slices (`gen_k4_slices`), full second-order injection | 2.38e-4 (one MLP blew up to 1.1e-3) | worse than first order |
| oracle: true κ₃/κ₄ but joint κ₄ replaced by the **exact** (MC) age-0 generated slices, own nets | NaN | age-0 joint κ₄ is the wrong object: it destabilises |
| oracle: all cumulants replaced by exact age-0 ones | 2.3e-3 | old content is essential at both orders |
| HD-2k4: all-age κ₄ tensor with exact single-site κ₄(a) sources, transported with decoupled gates (full n⁴ tensors) | 1.32e-4 | single-site sources are not where joint κ₄ comes from |

Charged to the realisation: the joint κ₄ slices that the readout needs (R3) are dominated by *old*, *multi-site* content. That content is generated by κ₃ passing through later ReLUs (δ-insertion terms, κ₃ ⊗ κ₃ → κ₄) and by multi-site diagrams. Neither is in the transport used here (linear, decoupled gates, first order in the incoming cumulant). Capturing it is the second-order Duhamel term, (ρ − ν)[g − ĝ] with ĝ itself first-order corrected, and it was not built in this session.
