# Seeds, tiles and sophistication

Working note II (self-reviewed draft): the input of a random ReLU network as a seed,
depth as critical stripping of input randomness, exact seed identity and lower bound,
and the network's own tile hierarchy.

- `seeds_tiles_sophistication.pdf` / `.tex`: the note.
- `code/`: `seeds.py <n> <L> [Ktruth]` and `seeds32.py` (challenge shape, float32) compare i.i.d.,
  antithetic and rotated cross-polytope seeds; `stripping.py` measures frozen activation bits.
  These read `spec.npy`, produced by `../pseudorandom-arrows-pv/code/spectrum.py`.
