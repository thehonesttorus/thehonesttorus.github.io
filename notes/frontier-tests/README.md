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
