# Trees in the gaps

Working note IX (self-reviewed draft). It asks how the code graph (object 1) and the weight tree (object 2)
should work together to estimate activation means from the weights. The answer follows Connes's treatment
of Cantor sets (*Noncommutative Geometry*, IV.3.epsilon), where the quantized differential lives on the gaps.

- **Gaps are arrows.** The mean is exactly a sum over gaps (bit-flip facets): Gaussian Minkowski content of
  each facet times the jump of the gradient across it. This sum cancels massively (the Hausdorff-versus-harmonic
  phenomenon), so it has to be summed analytically.
- **Every neuron is a node and a gap.** In the Hermite expansion at the neuron's own scale:
  - the degree-1 jet sigma*P(bit on) is the node weight;
  - the jets of degree >= 2, (-1)^j sigma He_{j-2}(a) phi(a), are gap weights (derivatives of the density at the kink);
  - the weights enter only as edges (W and the correlations R).
- **Tree theorem.** Cumulants are sums over connected diagrams, with nodes as leaves, gaps as internal vertices
  and correlations as edges.
  - Tree-shaped diagrams are exactly those computable by leaf-to-root message passing, at O(n^3) per shape.
  - Each cycle costs a further factor max|rho|. Measured: the first cycle is 1% of the trees at every depth to 16.
- **Ideal hierarchy.** Under fresh He weights:
  - Gaussian closure is O(1): the Dixmier part.
  - Trees are O(1/n), with kappa_3 terms random-signed and kappa_4 terms coherent.
  - Each cycle multiplies by max|rho|.
- **Depth.**
  - Non-Gaussianity is forgotten at rate 1 - 2<P(1-P)> per layer, i.e. twice the mean code-bit variance
    (verified to within 4%).
  - ReLU homogeneity makes the **scale mode exactly critical** (eigenvalue 1 at every layer): the pole.
- **The residue is the gain.** At depth, the closure error is mostly one scalar times the means, and that scalar
  is the input-to-input variance of the network's gain: c ~ Var(G^2)/8, within 15% from layer 5 on.
  - The gain injected per layer has an O(n^2) weights-only formula (checked against MC to 4-18%).
  - Its transport is near-critical (tau = 0.95 used empirically).
  - c also has a random, network-dependent part.

Numbers (He, Gaussian inputs, ground truth 1.6e7 / 4e6 / 1e7 samples at widths 256 / 512 / 1024):

- At layer 2, tree-level propagation (TLP) removes the closure error down to the noise floor.
- TLP improves on closure 16x at depth 4, 10x at depth 8 and 8x at depth 16 (width 256). It costs 0.8-1.8 budgets
  as implemented, so it is an accuracy reference.
- Closure + gain residue costs ~3% of the WhestBench budget and cuts closure MSE 1.6-4.5x.
- At width 1024, depth 16 (the WhestBench Phase-2 shape, not used to tune anything):
  - closure + residue reaches MSE 1.50e-6 (oracle scale 1.48e-6; closure 4.1e-6; full-budget MC 1.03e-6);
  - its score under MSE x max(0.1, C/B) is ~7x better than full-budget Monte Carlo (5x with the parameter-free
    tau = 1, 2.5x for plain closure);
  - TLP reaches 8.9e-7 raw (2.3e-7 with the exact residue), but costs ~1.8 budgets as implemented.
  - Caveat: one network at width 1024; the residue has a network-dependent random part (three seeds at 256).

Files:
- `trees_in_the_gaps.pdf` / `.tex`
- `code/` (see `code/README.md`), with the quoted outputs in `code/outputs/`.
