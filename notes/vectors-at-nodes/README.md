# Vectors at nodes

Working note V (self-reviewed draft). It starts from the super-expander picture: rows at nodes, random input,
and the question of what graph or tree ReLU builds.

- A bias-free ReLU MLP is random-hyperplane message passing on node vectors.
- Random inputs code each neuron as a point of a zero-entropy Cantor set:
  - split depth is Geometric(angle/pi);
  - three neurons split first according to Gromov products;
  - the number of cylinders is Cover's count C(k,d).
- The mean is an edge functional of the tile graph, E F = sum_e gamma_{d-1}(e) * bend_e. This is the V_1
  edge formula for virtual Newton polytopes, with an exact neuron/layer decomposition.
- Expansion comes for free from the Gaussian input:
  - every tile graph satisfies a sharp L^1 Poincare inequality for every metric target;
  - lumped OU kernels have gap >= 1 - e^{-s}.
- What the network contributes is cancellation in the labels: absolute boundary mass ~ (L-1) sqrt(n)/pi
  against an O(1) mean.

- `vectors_at_nodes.pdf` / `.tex`: the note.
- `code/edge.py`, `code/edge_big.py`: edge/neuron formula versus Monte Carlo (output in `edge_big_output.txt`).
- `code/scan.py`, `code/scan_deep.py`: boundary mass and cancellation index versus width and depth
  (output in `scan_output.txt`).
- `code/ou_gap.py`: lumped Ornstein-Uhlenbeck spectra on activation tiles.
- `code/codes.py`: neuron-code split laws, Cover counts, collision masses.
