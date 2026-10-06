# The arrow mixture: a test of the pattern-window picture

Working note XII (results). It tests the proposal that the law of a deep layer should be carried as a measure on
short *arrows* (a window of activation-pattern labels over a few consecutive layers) with Gaussian fibres, a
mixture of cone-Gaussians indexed by a path, instead of a cumulant hierarchy. The claim to be tested was: "global
kappa_3 is mostly the between-chamber displacement; keep the label and centre, and the Gaussian closure becomes
exact on the complement of the window; the window is set by the border depth, not by the width."

**Implementation (`code/arrow.py`).** The law of each layer is a discrete base measure along k directions times a
shared Gaussian fibre. At every layer the fibre is split along k directions into Gauss-Hermite nodes (q per
direction, product grid, or a total-level sparse grid with the per-direction variance restored), each node is
pushed through the ReLU as a conditional Gaussian, and the fibre is the node-averaged post covariance. Components
born at a layer persist through w-1 further ReLU steps (the path memory, the border depth) and are then merged by
moment matching, their spread returned to the fibre. Second moments are exact at every merge; w = 1 is the
ordinary Gaussian closure (checked: 1e-31). Bases: `eig` (top-k eigenvectors of the fibre), `neuron` (top-k neurons
by single-wall kappa_3, value nodes), `sign` (the literal cone split of a neuron into its two half-spaces,
truncated-Gaussian children), `mean` (the mean direction), `rand` (control). `gain` adds the GAC scale mixture as a
radial coordinate with infinite memory (k = 0 with gain reproduces GAC to 1e-29).

**Prediction from the fresh-weights lemma.** The next layer's weights are independent of everything the closure
computed, so each neuron reads a random projection of the whole state error. A window of k neurons out of n then
carries about k/n of the single-wall injection and (k/n)^2 of the wall-pair injection, and pairs dominate at depth
(93% of the gain injection at width 1024, layer 16, note XI). The radial zero mode (the gain) is the one coordinate
whose share does not shrink with the width.

**Results, final-layer MSE against Monte Carlo truth.**

Width 64, depth 8, four networks, 1e7-sample truth (`code/grid64.py`):

| estimator | MSE | ratio to closure |
|---|---|---|
| Gaussian closure | 6.01e-4 | 1 |
| K3 chain (shape-generic base path, no kappa_4 riders) | 3.7e-4 | 0.62 |
| **GAC** (radial coordinate only) | **1.33e-4** | **0.22** |
| eig k=1, q=7, memory 1 / 2 / 3 / 4 (no gain) | 5.84 / 5.65 / 5.28 / 5.20e-4 | 0.97 / 0.94 / 0.88 / 0.87 |
| mean direction k=1, memory 3 | 5.71e-4 | 0.95 |
| random direction k=1, memory 3 (control) | 5.87e-4 | 0.98 |
| neuron k=1, q=7, memory 3 | 4.94e-4 | 0.82 |
| sign (cone) k=1 / k=3, memory 3 | 7.30e-4 / 6.73e-4 | 1.22 / 1.12 |
| GAC + neuron k=1, memory 2 / 3 | 1.21 / 1.19e-4 | 0.20 |
| GAC + neuron k=2, q=5, memory 2 | 1.72e-4 | 0.29 |
| GAC + sign (cone) k=4, memory 2 | 3.18e-4 | 0.53 |

Width 256, depth 16, three networks (`code/grid_n.py`, `code/diag256.py`):

| estimator | mean of 3 networks | ratio to closure |
|---|---|---|
| Gaussian closure | 6.04e-5 | 1 |
| **GAC** | **1.74e-5** | **0.29** |
| neuron k=1, q=7, memory 1 / 2 / 3, no gain (network 0; closure 6.29e-5) | 6.21 / 6.20 / 6.82e-5 | 0.99 / 0.99 / 1.08 |
| GAC + eig k=4, q=3, level-2 sparse grid, memory 2 | 1.60e-5 | 0.27 |
| GAC + neuron k=4, q=3, level-2 sparse grid, memory 2 | 5.51e-5 | 0.91 |
| GAC + neuron k=8, q=3, level-1 sparse grid, memory 2 | 5.43e-4 | 9.0 |

Width 1024, depth 16, official WhestBench network 0 (1e9-sample truth):

| estimator | MSE |
|---|---|
| Gaussian closure | 4.06e-6 |
| **GAC** | **1.56e-6** |
| GAC + eig k=4, memory 2 | pending |
| K3 chain (reference) | 2.19e-8 |

**Reading.**

- The radial coordinate is the whole story of the coherent part: GAC removes 78% (width 64), 79% (256), 62%
  (1024) of the closure MSE with one scalar per layer and no fitted constant.
- A directional base on top of it is worth 10% at width 64 (one neuron), nothing at width 256, and it turns
  negative as soon as it is carried for more than one layer. More directions make it worse, not better. The
  literal cone split is worse than the closure at every width because two truncated-Gaussian children are a
  two-node quadrature for every other neuron.
- The reason it turns negative is specific. With a shared fibre the mixture carries the between-node skewness
  along the base but drops the coupling between the node position and the node's own post-ReLU covariance.
  That coupling is third-order structure (the (2,1) slice kappa_3(h_i, h_i, h_j)), and the factorised chain
  carries exactly it. Keeping per-node fibres would restore it at a cost of M n^3 per layer, M components:
  at width 1024 that is 1.5 M budget units per layer, so the exact mixture is priced out of the contest at
  M = 4.
- So the claim fails in the direction the lemma predicts: at width 1024 the third-order content of the law is
  bulk. There is no finite pattern window whose conditional Gaussians are exact on the complement; the
  window that would be needed is O(n) neurons (participation ratio 350-410 for the old births, note
  frontier-tests). The pattern transversal is the pathological quotient: no finite set of its coordinates acts
  coherently on the mean, and only its statistics do, which is why statistics are 300x cheaper than addresses
  per unit of accuracy.

**What the hull picture leaves standing.** The deep law has one slow coordinate, radial, with critical
accumulation (ReLU homogeneity), and GAC carries it exactly to all orders for O(n^2) a layer. Everything else
that matters for the mean is bulk third- and fourth-order structure, and the only polynomial-cost
representation of it is the factorised cumulant chain with its birth-time decoration. The leaderboard chain's
memoryless kappa_4 with a fitted lambda per layer is the fitted shadow of the radial zero mode (its diagonal
matches the gain scale-mixture law at the last layer, note frontier-tests). The construction the theory
dictates is therefore the chain run on the gain-conditional state with the output rescaled by E[G_L], no fitted
constant, which also tells the source window what it no longer has to carry.

**The GAC-conditioned chain (`frontier-tests/code/run_gain.py`, official network 0).** Running the K3 chain on the
gain-conditional state (GAC bookkeeping: E z = E[G] mu, E zz^T fixed; output E[G_L] x conditional mean) is
monotonically worse than the chain itself: injection x0 reproduces the chain (2.191e-8), x0.1 gives 4.27e-8, x1
gives 1.84e-6, and with the chain's kappa_4 channel switched off as well 1.98e-6. The chain with its kappa_4
channel off and no conditioning gives 2.11e-7. So the gain is already inside the chain: its sources and
second-order terms carry about 70% of E[G_L] - 1 and the kappa_4 channel the remaining 30%, and the scale
mixture on top double-counts it. The radial zero mode and the low-order marginal cumulants are two descriptions
of one object, and the chain's description is the more complete one at its accuracy (the gain's cumulants
beyond the fourth are O(gamma^2), about 6e-10 in MSE). There is no free win in exchanging the fitted kappa_4 for
GAC; what is left of the chain's error is the transport of bulk third-order structure.

Files (`code/`): `arrow.py` (the estimator), `grid64.py`, `grid_n.py`, `diag256.py`; printed results in
`outputs/`. The scripts import `closure.py` and `gac.py` from `../../gain-aware-closure/code` and the weights and
truth files from `../../trees-in-the-gaps/code/truth.py`.
