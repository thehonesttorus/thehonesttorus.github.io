# Frontier tests

These are measurements on official WhestBench Phase 2 networks (width 1024, depth 16, 1e9-sample truth), taken
while looking for a route to the top of the leaderboard.

**Reference points:**
- The leaderboard top has raw final-layer MSE 1.1-1.5e-8 at 11-15% of the budget.
- The open-source factorized third-cumulant chain (504aldo, MIT; ARC's arXiv:2605.05179) gives raw 2.1-2.3e-8 at
  25-27% here (`run_k3.py`).

**Generator / evaluator designs** (`ge_suite.py`). Each is tested in the strongest form we could build:

| design | effect on Monte Carlo cost at fixed MSE |
|---|---|
| plain Monte Carlo, 10% floor | MSE 1.2e-5 |
| radial split (exact identity) | variance x 0.993 |
| best possible importance proposal q* ~ phi ‖R‖ (exact ceiling for the averaged readout) | variance x 0.986 |
| antithetic pair | x 0.91 per pass (`varstruct.py`) |
| corrected restriction hierarchy (W1 singular frame, exact psi at layer 1, coupled levels, optimal multilevel allocation) | work x 15.4 |
| paired wall-source correction around the top-256 frame | Var(D) = 0.54 Var(F) at 3 passes |
| layer-1 frame regression (top 8 / 256 coordinates) | R^2 0.007 / 0.16 (`cv1024.py`) |
| orthonormal frames + antipodes + exact radius (the leaders' Phase 1 quadrature family) | MSE 7.2e-6 vs 1.2e-5 i.i.d. at 6144 passes, a 1.7x gain (`frames.py`) |

To match the frontier by sampling at the 10% floor, the variance would have to fall by about 573x.

**Old third-cumulant sources versus the gain** (`run_win.py`, `analyze_old.py`, `vfield.py`, `run_fill.py`):

| setting | final MSE |
|---|---|
| full chain | 2.19e-8 |
| 4 youngest sources only | 7.33e-7 |
| 6 youngest sources only | 3.21e-7 |
| 8 youngest sources only | 1.30e-7 |

- The gain's exact third-cumulant form fits the old (2,1) slice with R^2 rising from 0.04 at layer 6 to 0.43 at
  layer 15. A per-neuron gain field reaches 0.50, and it is 0.92 correlated with the mean.
- A parameter-free gain fill for the dropped sources improves the windowed chains by only 1-8%; 2-3x amplitude is worse.
- The full chain's residual has no gain or scale component, and no correction we can compute predicts it: the best
  removes 2.4%.

**Input-chaos windows** (`chaos_rank.py`):
- Each third-cumulant source is the second-chaos (expected-Hessian) term of one birth layer, a quadratic form in the
  input.
- A single deep birth is low-rank in input space: participation ratio 120 at layer 5, down to 50 at layer 11.
- The union of old births is not: participation ratio 350-410.

**Where the chain's error is** (official network 0, the true post-activation mean injected at the input of every
layer; `run_oracle_dump.py`, `noisefloor.py`, `decode.py`, `tdecode.py`, `subdecode.py`, `gating2.py`,
`mehler_probe.py`, `run_mehler.py`, `k4theory.py`, `run_g4.py`, `run_gaink4.py`):

- The one-step residual r_l = m*_l - chain_l(true mean in) is 90-96% signal after the 1e9-sample truth noise
  is subtracted; it grows from 0.9e-9 (layer 4) to 4e-9 (layers 14-15).
- One retention factor closes the budget: sum_l 0.9^(15-l) r_l = 2.16e-8 against the measured 2.19e-8. The
  defects enter a mean channel that loses almost nothing per layer, and the final error is about 70x a
  per-layer increment.
- Decoding r_l through the fresh weights (the next-layer error of neuron i is w_i^T Delta w_i for a covariance
  defect Delta, and the cubic contraction for a kappa_3 defect) finds no low-dimensional structure: the
  gain, the Mehler orders, the kappa_4 shapes, the D21 scale and the common mode explain at most a few percent,
  one step back or transported through all layers, in and out of sample (R^2 0.02-0.04); the top covariance
  eigen-subspaces explain no more than random subspaces of the same dimension.
- The kappa_3-driven terms are large (readout skew 5e-4, variance moved by D21 2-3e-4 per neuron) against
  r ~ 6e-5, so an incoherent 10% error of the bulk third-cumulant tensor accounts for the whole residual.
- One ingredient is identified: the chain gates every all-distinct kappa_3 entry by Phi_a Phi_b Phi_c, while
  the exact first-order gate is the joint orthant probability. The difference is 1.3-1.9% of the gate per
  layer, mean zero, spread over the bulk (its top coherent mode carries 0.1-0.3%), and the first-order Mehler
  pair term phi_a phi_b R_ab Phi_c captures it to 0.1%.
- Completing the Gaussian Mehler series of the chain's off-diagonal covariance to order 5 (it stops at order 2)
  changes the final MSE by 0.05%: the missing orders are 1e-6 of the second-order term.
- The kappa_4 diagonal is not the bottleneck: the parameter-free gain scale-mixture kappa_4 at the last layer
  gives 2.185e-8 against the chain's 2.191e-8 (it matches the fitted lambda); at all layers 2.298e-8; scaling the
  chain's kappa_4 to the Monte Carlo value makes the last step worse (4.18e-9 -> 7.29e-9).

**The chain already contains the gain** (`run_gain.py`, `outputs/gain_conditioned_chain.txt`): see
`../arrow-mixture/README.md`. Running the chain on the gain-conditional state with the output rescaled by E[G_L]
double-counts the radial zero mode (2.19e-8 -> 1.84e-6); the chain's own kappa_3 sources and second-order terms
carry about 70% of E[G_L] - 1 and its kappa_4 channel the remaining 30% (kappa_4 channel off: 2.11e-7).
