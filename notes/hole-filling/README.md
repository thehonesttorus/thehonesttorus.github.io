# Removing is injective: hole filling in the Penrose hull, and what it says about the chain

Working note XIV. The question: a Penrose tiling with a bounded hole has exactly one legal filling, determined by the
tiling outside the hole. Why, what it costs, and which parts of the estimation programme are of that kind.

## 1. The fact and its mechanism

Take the cut-and-project description: vertices are the projections pi(m) of lattice points m with internal image
pi^perp(m) inside a window W + gamma (for the rhombus tiling, pentagons in a two-dimensional internal space, edges along
lattice directions). Three steps.

- **The outside fixes the parameter.** The internal images of the vertices outside any ball are dense in the window
  and those of the non-vertices are dense in its complement (equidistribution), so the closed window, hence gamma, is
  determined by the complement of any bounded set.
- **The parameter fixes the hole, up to the boundary.** A lattice point of the hole is a vertex iff its internal image
  is inside W + gamma. The only ambiguity is a lattice point whose image lies on the boundary: the singular case.
- **The ambiguity is unbounded.** A window edge through gamma is a segment in a lattice direction, and the lattice
  points whose images lie on it form a dense subset of the segment whose physical projections run along an infinite
  line: a Conway worm. Singular tilings come in pairs (Robinson's positive and negative resolutions) that differ along
  the whole worm, and a worm is exactly what gives a half-plane two completions.

So two Penrose tilings that agree outside a bounded set agree everywhere: T -> T restricted to the complement of B is
injective on the hull for every bounded B. The torus parametrization is one-to-one off the singular set, and on the
singular set its fibres contain no pair at bounded distance.

**The counterexample next door.** For the Fibonacci chain (any Sturmian word) the window is an interval, its boundary
is two points, and a singular parameter puts one lattice point on each. The two resolutions differ in two consecutive
letters, ab against ba: the upper and lower mechanical words of the same slope and intercept differ at exactly two
neighbouring positions when n alpha + rho is an integer. Removing a two-letter patch is not injective on a dense
measure-zero set of the hull. The dichotomy is the codimension of the window boundary: one-dimensional internal space
gives point ambiguities and local flips; a two-dimensional internal space with polytopal windows in lattice directions
gives worms. In dynamical terms the hull is an almost one-to-one extension of its maximal equicontinuous factor (the
torus), and "removing is injective" is the statement that no singular fibre contains two points at bounded distance.

## 2. What filling costs

Injectivity says the outside determines the hole; it does not say from how far.

- A lattice point of the hole whose internal image lies at distance eps from the window boundary is decided only once
  gamma is known to within eps. The outside to radius rho pins gamma to about 1/rho (two physical and two internal
  dimensions), so that point needs the outside to radius about 1/eps.
- A hole of radius R holds about R^2 lattice points; the one nearest the boundary sits at eps ~ 1/R^2. Full
  determination needs the outside to radius ~R^2. The forcing radius of a hole is quadratic in its size.
- Given the outside only to radius rho, the undetermined points are those within 1/rho of the boundary, about R^2/rho
  of them, lying on about R/rho worms. But the consistent fillings are the cells into which the lines through those
  images cut the 1/rho-ball of possible gammas: polynomially many. The conditional entropy of the hole given the
  outside to radius rho is about 2 log(R^2/rho) bits, not R^2/rho bits.
- The filling program is therefore: the parameter to the needed resolution, then cut-and-project (linear in the
  number of tiles). A hole of R^2 tiles is specified by O(log R) bits. Periodic tilings need O(1); a random tiling
  needs R^2. Conway's empires are the converse: a bounded patch forces an unbounded, sparse set of far tiles along its
  Ammann bars, because the bars are straight and the parameter they fix is global.

## 3. The dictionary, with today's numbers

- **The defect of conditional filling.** Filling a hole by conditional expectation E[. | outside] fails the
  Rota-Baxter identity of Connes-Kreimer by exactly the conditional covariance E(ab) - E(a)E(b) (the identity in the
  pasted critique). Removing is injective iff that covariance vanishes for every hole, iff the conditional expectation
  is a character on the hole algebra, iff local counterterms renormalise exactly. The chain is a conditional-expectation
  scheme on its retained state (the Gaussian data and the carried kappa_3 slices); its defect was measured in E7, E8,
  E9 and the frontier tests: not small in any affordable address. It is not a Birkhoff decomposition, and the earlier
  Connes-Kreimer note was wrong to present it as one.
- **The coherent part of the defect is the gain, and it lives in the odd channel.** E10: the one rigid, forced
  structure is the mean-coupled third cumulant kappa_3(z_k) = 1.5 g mu_k sigma_k^2, transported with retention 0.9 and
  fed by a one-loop source; its "empire" (what it forces at infinity) is the critical mean channel. But an old source
  does not force the output only through this channel: dropping the sources older than 8 transports costs 6x in the
  chain and refilling their gain content recovers 1-8% (frontier tests). What survives of an old source is quenched,
  18-23% of the (2,1) energy per layer, confined to the slowly decaying top Lyapunov directions (rank ~3n/8 at age 4),
  and it is exactly the content whose carrier costs n^2 r per leg per layer. That is the worm with a bent trajectory:
  pinned by neither end.
- **The maximal equicontinuous factor of depth is the dilation group.** ReLU homogeneity is the one exact symmetry of
  the layer map; the radial coordinate is the one address whose fibres the outside fills (E4: quasi-free fibres only
  under radial conditioning; the arrow mixture's gain base; GAC). Everything else is mixing: no further eigenvalues,
  no further forced addresses, and E7's residual is a function of no per-neuron state.
- **Worms are walls, but ours bend.** A worm is straight, so its two ends pin it. A wall of layer l >= 2 is a
  piecewise-linear hypersurface bent at every earlier wall it crosses; the bends are the kappa_3 trees. Its forcing
  radius is maximal: the whole weight tensor, which is the content of "quenched". The structures that are straight
  are the first-layer hyperplanes, the dilation direction, and at depth the mean direction, which is the top
  covariance eigenvector (overlap^2 0.94 at layer 14) with two thirds of its variance the gain's rank-one term. These
  are the network's Ammann bars, and they are exactly what the existing estimators already read. E9 measured that
  the first-layer bars force less than 1% of the deep tiles.
- **Entropy of the hole.** The chain's residual at 2e-8 needs bulk third-cumulant data of the size of the hole; the
  program is as long as the hole. The exchange rate of 300x between statistics and addresses is this statement.

## 4. What it buys

A clean statement of the floor, and the identification of the single Penrose-like structure (the odd gain channel and
its derived transport, E10). It does not supply a filling rule the outside can run cheaply, because the hull of the
code process is not an almost one-to-one extension of a low-dimensional factor beyond the radial one.
