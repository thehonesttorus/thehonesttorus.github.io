# The gain-aware closure

Working note XI (results so far; the full write-up is still to come). It replaces the empirical transport constant
tau = 0.95 of notes IX-X with a derivation. It evaluates the result on the official WhestBench Phase 2 networks,
with a flopscope implementation, and tests three proposals from other conversations: a "Ruelle bundle" over sign
chambers, a "dual truncation" over wall codimension, and a "cyclic filtration" probe.

**Exact ledger.** The closure error at the output is exactly a sum of local defects, each carried to the output by
the closure itself:

- closure - truth = sum_l (R_l - R_{l+1}), where R_l is the closure started from the true moments at layer l.
- Telescoping is checked to 1e-19 (`ledger.py`).
- Each local defect is short-memory tree terms plus a persistent gain-mixture term (gamma/8) sigma phi(a)(1 + a^2)
  (`local_compare.py`, `mixture_local.py`).
- The gain behind that mixture term accumulates critically (tau = 1).

**Gain-aware closure (GAC, `gac.py`).** Carry the pre-activations as a scale mixture z = G y, with
y ~ N(mu, S), E G^2 = 1 and Var G^2 = gamma.

- ReLU is homogeneous, so a scale mixture passes through a layer unchanged. The gain therefore accumulates
  critically: gamma_l = gamma_{l-1} + f_l, where f_l is the weights-only O(n^2) injection of note IX.
- Each layer matches the marginal moments (E z = E[G] mu, E zz^T = S + mu mu^T) and propagates the conditional
  Gaussian with the ordinary closure.
- The output is E[G_L] times the closure mean of y_L.
- No fitted constant. The discount that tau = 0.95 imitated comes from the closure's own transport of the
  conditional state.

Final-layer MSE against Monte Carlo truth:

| network | closure | tau = 1 | tau = 0.95 | **GAC** | oracle scale |
|---|---|---|---|---|---|
| 256, seed 0 | 6.29e-5 | 2.12e-5 | 1.40e-5 | **1.33e-5** | 1.32e-5 |
| 256, seed 1 | 7.35e-5 | 3.11e-5 | 2.13e-5 | **1.94e-5** | 2.04e-5 |
| 256, seed 2 | 4.49e-5 | 6.72e-5 | 2.73e-5 | **1.95e-5** | 2.03e-5 |
| 512, seed 0 | 1.61e-5 | 6.55e-6 | 5.93e-6 | **6.64e-6** | 5.58e-6 |
| 1024, seed 0 | 4.14e-6 | 2.09e-6 | 1.50e-6 | **1.44e-6** | 1.48e-6 |
| 1024, seed 1 (new) | 4.43e-6 | 2.46e-6 | 1.57e-6 | **1.76e-6** | 1.57e-6 |
| 1024, seed 2 (new) | 3.30e-6 | 1.78e-6 | 1.27e-6 | **1.35e-6** | 1.26e-6 |

- "Oracle scale" is the best single scalar applied to the closure mean, fitted to the truth.
- GAC changes shape as well as scale, so it can beat the oracle scale (width 256, seeds 1 and 2; layers 2-7
  everywhere).
- At width 1024, GAC still under-corrects slightly: 8 x the residual scale is +0.002 to +0.003.

**Official WhestBench Phase 2 networks** (`code/official/`). These are the 100 networks of the `mini` split of
`aicrowd/arc-whestbench-public-2026@v2-phase2`, with 1e9-sample truth. Final-layer MSE, mean over the networks:

| estimator | mean MSE | flopscope cost (measured) | score = MSE x max(0.1, C/B) |
|---|---|---|---|
| Monte Carlo, full budget | ~1.2e-6 (avg_var / 65,536) | 100% | ~1.2e-6 |
| bundled covariance baseline / closure | 3.85e-6 | 2.35% | 3.9e-7 |
| closure + residue, tau = 0.95 (fitted at width 256) | 1.37e-6 | | |
| **GAC** (no fitted constant) | **1.42e-6** | **2.38%** | **1.4e-7** |
| **GAC + kappa_3 trees** (48 networks so far, 48/48 better than GAC) | **1.12e-6** | **8.25%** | **1.1e-7** |
| oracle scale (fitted to truth) | 1.31e-6 | | |

- tau = 0.95 beats GAC on 74/100 networks, by about 3% on average. GAC's residual scale is systematic:
  8 x scale = +0.0029 +- 0.0014.
- The trees remove that residual and also correct shape, so they beat the oracle scale.
- The trees are the fresh one-point kappa_3 terms of the previous layer (`gac_k3.py`):
  - D3: single walls, about a quarter of the gain;
  - P3 with Hermite orders <= 2, and T3: wall pairs, about three quarters.
- Flopscope implementation (`code/flopscope/gac_flopscope.py`):
  - float32 throughout, with the centred covariance tagged symmetric as in the bundled baseline;
  - the relu^2 Hermite coefficients come for free from B_j = 2 sigma A_{j-1};
  - it reproduces the numpy MSE exactly;
  - wall time is 4.8 s and residual Python time 0.14 s on an idle machine (limits 120 s and 0.4 s).

**Tests of the other proposals.**

- **Ruelle bundle over 2^6 sign chambers** (`bundle_local.py`). Tested at its strongest, with exact chamber
  probabilities and chamber-conditional moments taken from 1e6 samples. The test asks how much of the one-step
  local defect each partition removes (MSE relative to the closure step; random labels give 0.996-1.000):
  - signs of the top-6 right singular vectors of W (the proposal): 0.72-0.96, and the scale part is unchanged (one layer is undefined because a chamber is almost empty);
  - signs of the top-6 covariance eigenvectors of the layer (oracle moments): 0.82 at layer 2, falling to 0.12
    at layer 16;
  - 64 bins of the gain |h|^2 (oracle moments): 0.27-0.44;
  - GAC's analytic one-dimensional gain base, no samples (`gac_local.py`): 0.27-0.50.
  - The transversal matters. The data's depth-outlier subspace is the one worth conditioning on.
- **Diffraction** (`bragg.py`, `diffraction.py`). The law of the deep representation is absolutely continuous:
  there are no peaks up to t = 64.
  - Its characteristic function is a Gaussian broadened by the gain mixture, and GAC predicts it from the weights.
  - At t = 3: measured 0.026, Gaussian 0.011, GAC 0.023.
- **Codimension hierarchy** (`codim.py`). This splits the one-step tree terms into single-wall terms (D3, D4) and
  wall-pair terms (P3, T3, PP4, PD4).
  - Single walls dominate only at layer 2. From layer 4 on, wall pairs dominate. At width 1024, layer 16, pairs
    carry 93% of the gain injection and 87% of the rms mean shift.
  - The pair sum is not sparse. The top n pairs by |rho| carry -3% to +3% of it, and the top 10n carry at most 18%.
  - Its coherent part is the quadratic form u^T F0 u, which GAC evaluates exactly in O(n^2).
  - The same split for the kappa_3 mean corrections on official networks: single walls cut GAC's MSE 5-8%, and wall
    pairs cut it 19-22%.
- **Cyclic-filtration probe** (`probe64.py`). Width 64, depth 8, 10 networks, 1e7-sample truth. MSE relative to the
  closure:
  - wall page (closure + diagonal wall current): 0.69;
  - 16-cylinder page (exact cell moments): 0.81;
  - GAC: 0.28.
  - The page difference d1 correlates only weakly with the closure residual: mean +0.23, range -0.22 to +0.67.
    Each page alone correlates +0.66 / +0.74, because they capture the same part.
  - The non-scale forgetting rate 2<P(1-P)> stays >= 0.15 per layer, a Theta(1) gap. The residual is
    dominated by the exact zero mode, the gain.

Files (`code/`):
- `gac.py` (the estimator), `gac_k3.py` / `gac_tree.py` (GAC + trees), `gac_local.py`
- `flopscope/gac_flopscope.py` (contest implementation), `flopscope/run_flopscope.py`
- `official/` (extract the official parquet shards, evaluate), `probe64.py`
- `ledger.py`, `hmoments.py`, `transport.py`, `channels.py`
- `local_compare.py`, `mixture_local.py`
- `bundle_local.py`, `codim.py`, `bragg.py`, `diffraction.py`
- printed results in `code/outputs/`

These scripts import from `../../trees-in-the-gaps/code`. Run them as:

```
PYTHONPATH=../../trees-in-the-gaps/code python3 gac.py 1024 16 0
```

Weights and truth files come from `truth.py` there (`hmoments.py` for the per-layer second moments).
