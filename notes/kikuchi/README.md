# Kikuchi algebras and adaptive final means: implementation and first measurements

Source: `kikuchi_algebras_and_adaptive_final_means_1.pdf` (10 October 2026 note).

## Code

- `whest/kikuchi.py` holds every construction of the note:
  - moment interval and rank-two cut relaxation (Sec. 4);
  - bounded-effect restart (Sec. 5) and squared-shell state (Sec. 6);
  - even and odd cyclic hierarchies with the residual certificate, via Lanczos in mpmath (Secs. 7–8);
  - pair overlap, pair up/down chain, three lengths and two-sign heat bath (Sec. 9);
  - conditional brackets, cell intervals and injection (Sec. 10);
  - gate tower (Sec. 3);
  - Kikuchi Laplacian and sparsifier congruence (Sec. 11);
  - Schur refinement with enrichment (Sec. 12);
  - observable-code bound (Sec. 13).
- `scripts/check_kikuchi.py` runs 42 checks, and all pass (`check_kikuchi.txt`). They include:
  - the Sec. 9.1 table L1, L2, U2, U1, R1, R2, reproduced to 1e-10;
  - the Rayleigh closure 1/(2√π);
  - Example 3.2;
  - Theorem 3.1 at k = 0, 1, 2;
  - exact spectral termination;
  - the remaining finite statements.
- `scripts/kik_contract.py` is the single-pass contraction. One float32 forward pass per 4096-input batch accumulates:
  - every layer's power sums up to order 8 (16 at the last layer);
  - E z_+;
  - nested penultimate-gate cells (K = 8, pilot greedy order, eq. 22);
  - single-gate candidate matrices;
  - pair-overlap sums.
  
  Run on AWS: 8 nets × 4.3M inputs, 96 tasks on 352 vCPU, finished in **266 s** with streamed results.
- `scripts/kik_merge.py` builds the per-neuron certificate pipeline (`kik1_certificates.txt`).
- `scripts/kik_ladder.py` gives the shape ladder per layer (`kik1_shape_ladder.txt`).
- `scripts/kik_chain.py` and `scripts/kik_chain_profile.py` run the diagonal-cumulant chain with oracle cumulants (`kik1_chain.txt`).

## What the measurements say (nets 0–7, final layer; RMS over 1024 neurons; target 6e-5)

| certificate / estimate | RMS | comment |
|---|---|---|
| first interval [t_+, (t+s)/2], exact t, q | cert 3.4–4.1e-2 | half-gap Δ ≈ 0.06; most final neurons are strongly one-signed |
| restart effect E_ν g_c | ∈ [0, 0.80] | must be known to ≈1.6e-3 to reach the target |
| cyclic hierarchy, k = 1…4 (moments to E Z^16) | kept cert 2.8e-2 → 7.2e-3 | algebraic convergence (kink of sqrt at 0); R_k/(U_k−L_k) grows with k |
| nested cells on 8 penultimate gates | cert × 0.82 | about 2.5% per gate; the best single gate gives a median 4% (max 26–31%) |
| pair overlap m = EP − κ E_ν 1/max | correction ≈ 9.5 = EP | corr(P,Q) ≈ 0.91; cancellation of two terms of size ~10 |
| Gaussian shape at exact (t, q) | 4.2–4.8e-4 | the final preactivation's shape defect |
| Hermite shape with exact κ3 | 1.3–1.6e-4 | |
| **Hermite shape with exact κ3, κ4** | **1.0–1.15e-5** | 5× below target; same-sample comparison |
| Gaussian closure contractions (A2) | mean 1.5–2.3e-3 | t err 2.2–3.2e-3, q err 3.9–5.6e-3; its interval covers the truth for only 73–81% of neurons |

- **Shape ladder at every layer.** With each layer's own (t, q, κ3, κ4), the four-term Hermite mean map matches E z_+ to
  3.8e-5 at layer 1, falling to 1.1e-5 at layer 15. With κ3 only it is about 1.5e-4; with the Gaussian about 4.5e-4.
- **Accuracy each final contraction needs** (each alone) to reach 6e-5:

  | contraction | needed accuracy |
  |---|---|
  | t | 9e-5 |
  | q | 2.7e-4 |
  | κ3 | 7e-3 (≈7% of its RMS) |
  | κ4 | 2.4e-2 (≈37% of its mean) |

- **Diagonal-cumulant chain** (exact mean recursion t_l = W_l m_(l−1), four-cumulant mean map, oracle per-neuron κ3, κ4):

  | variant | final mean error |
  |---|---|
  | Gaussian closure | 2.0e-3 |
  | oracle κ3, κ4 (closure variances) | 1.7e-3 |
  | oracle κ3, κ4 + oracle per-neuron variances (closure correlations) | 1.2e-4, about the noise injected by the oracle variances |

  With correct variances fed in, the one-step relative variance error the Gaussian pair kernel makes grows from
  7e-4 at layer 1 to 4% at layer 15.

So the note's certificates are valid on the real networks, but at feasible order they are 100–600× too wide to drive the
estimate. The sharp quantitative message concerns which contractions the final means need. They need:

- the four diagonal cumulants (t, var, κ3, κ4) of each layer's pre-activations;
- the exact linear mean recursion between layers.

The remaining obligation is the pre-activation **variance**, which rests on the off-diagonal post-activation covariance.
That is exactly where the Gaussian pair kernel fails. The means themselves are not where it fails.
