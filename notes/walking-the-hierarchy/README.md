# Walking the hierarchy

Working note VII (self-reviewed draft). It builds the expander analogue proposed in discussion: a Markov
walk on a Cantor set whose transition kernel depends only on the last point of contact x^y of two tree
paths. That vertex also defines the ultrametric d(x,y) = kappa(x^y). Convolution in the groupoid algebra
is then a Markov step.

- Exact spectrum on any tree: lambda_v = 1 - P(contact above v) - psi(v) mu[v] on each Haar detail space,
  together with a converse (which spectra come from a Markov kernel) and a convex program for the optimal walk.
- Explicit walk bounds: exact variance formula and Hoeffding tails.
- No free lunch: walks with nonnegative spectrum (expander walks, heat kernels) never beat independent samples
  per evaluation. Gains need sibling jumps (negative eigenvalues) or stratification over the tail groupoid.
- kappa from the weights: if F is L-Lipschitz for kappa in box dimension D, the error is O(L T^(-1/2-1/D)),
  and this rate is optimal.
- Random ReLU networks (width 256): gains from a frame built from the weights alone fall from 1.30x to 1.06x
  to 1.02x at depths 4, 8, 16; pilot frames give 1.51x, 1.33x, 1.15x. Mean-field reason: the variance sits at
  Hermite degree ~ L^2 (critical Galton-Watson law, note IV).

Files:
- `walking_the_hierarchy.pdf` / `.tex`: the note.
- `code/tree_spectrum.py`: spectrum formula, reversibility, walk variance vs. simulation.
- `code/energy.py`: stratification gains on exact Gaussian Cantor codings (outputs in `code/energy_*_output.txt`).
