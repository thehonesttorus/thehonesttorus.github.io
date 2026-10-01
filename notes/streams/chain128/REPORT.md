# chain128 — does a better D21 interface give a lower final-layer MSE?

*Stream report, started 2026-10-01. Status: in progress (code written, ground truth baking, atlases pending).*

## Question
Build a dense K=3 moment chain at width 128 (full κ3, n³ = 2M entries) and test end to end whether the closures
identified by the oracle ladder (notes/competition-plan.md §3.1) translate into lower final-layer MSE, and measure
the error law MSE vs ε(D21) at width 128.

## Method (code in this directory)
- `chain.py` — dense chain: μ, C, full κ3(z), κ4(z)_{aabc}; ReLU step by the Edgeworth operator (first order in κ3, κ4,
  plus κ3²/2) on exact Gaussian expectations of derivatives of reluᵖ (Mehler series to all orders in C); all-distinct
  κ3(a) from `oracle_k3.hermite_model` / `residual_basis` / `CLOSURE_COEF`; memoryless κ4 closure transported exactly.
- `validate.py` — layer-by-layer check of the ReLU step fed the atlas's true pre-activation cumulants.
- `fit_coefs.py` — per-layer coefficient table (variant D).
- `run_variants.py` — variants A0, A, M, B, C, D, E, F, N (error law) against baked truth.

## Ground truth
`whest dataset bake --n-mlps 8 --n-samples 1e7 --width 128 --depth 16`, seeds 128000–128007 (split `dev`).

## Reference implementation
ARC's `mlp_kprop` (arXiv 2605.05179): `git clone https://github.com/alignment-research-center/mlp_kprop` asks for
credentials (private or nonexistent); not on PyPI. Not available; the dense chain here is the reference.

## Results
(pending)
