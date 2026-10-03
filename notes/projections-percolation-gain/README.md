# Projections, percolation and the gain

Working note X (self-reviewed draft). It tests three sets of ideas against WhestBench Phase 2 (width 1024,
depth 16, He ReLU, Gaussian inputs):

1. the cone-tiling / groupoid estimator from another conversation;
2. the percolation-and-states construction prompted by Maxmud-Suomala, arXiv:2607.01291;
3. a "type-III spikes / KK-pairing" reading from another assistant.

**Rules, verified from the WhestBench source.**
- Score = final-layer MSE against 1e9-sample truth x max(0.1, C/B), with B = 2^41.
- flopscope only: matmul ~2mnk, float64 billed 2x, transcendentals 16x.
- Residual Python time is capped at 0.4 s.
- The bundled covariance baseline costs 5.2e10 FLOPs (2.4%) and reaches MSE 4.17e-6, the same as our closure.

**Exact identities that hold.**
- Projection memory is criticality: phi(e^2) = phi(e) for gate projections is He criticality. Mean-field loses
  2<phi(e) - phi(e)^2> per layer.
- Capacity is variance: the Lyons energy is the second moment of a pruned (Horvitz-Thompson) computation. With
  signed contributions it is divided by the squared cancellation ratio.
- The residue is a pressure gap: under the dilation flow, scale mixtures pass through layers unchanged, and the
  closure error is e^{P(1)} - e^{P(2)/2}, with P(beta) = log E G^beta.

**What fails (measured).**
- Cone mixtures: 3% gain with the proposed W1 frames, at most 13% with the best frame at depth 16, at 64 closures.
  The cylinders carry 3-4% of the log-gain variance.
- Percolation-sketched tree sums: relative error 1.7-36 at 50-1% of the work. The sums are random-signed, so
  compression is impossible.
- "Type-III spikes": the proposal's own test fails. There are 0 spikes above the Marchenko-Pastur edge in 80 He
  matrices, and the layer algebras are abelian (type I).
  - What is true: depth creates outliers in the layer covariance, and 5-45% of the non-scale residual sits on its
    top 4 eigenvectors.

**What survives: measure the gain, not the bias.**
- gamma = Var(G^2) is a coherent average over n^2 neuron pairs, measurable from 1000 samples (1.5% of budget)
  with no fitted constant. Then m <- (1 - gamma/8) m_closure.
- Width 1024: MSE 1.56e-6 (weights-only 1.50e-6, oracle 1.48e-6, full-budget MC 1.03e-6). At the 0.1 multiplier
  that is ~6.6x better than Monte Carlo. One network.

**Where percolation is real.** Depth is a critical branching process (note I). Maxmud-Suomala Theorem A(i) suggests
a criterion for when finite-width gain stops the deep representation from collapsing. This is a conjecture, and
it concerns depth much larger than width.

Files:
- `projections_percolation_gain.pdf` / `.tex`
- `code/`:
  - `cylinders.py` (cone mixture)
  - `perc_sketch.py` (percolation sketch)
  - `gainstruct.py` (log-gain structure, tilt test)
  - `calib_gain.py`, `calib_gain2.py` (residue calibration)
  - `spikes.py` (spike test)
  - `outputs/`

These scripts import `closure.py` and `residue.py` from `../trees-in-the-gaps/code`. Run them with:

```
PYTHONPATH=../../trees-in-the-gaps/code python3 calib_gain2.py 256 16 0 250,1000 20
```

The weight and ground-truth files are produced by `truth.py` there.
