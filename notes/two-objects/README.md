# Two objects

Working note VIII (self-reviewed draft). It sets up two objects:
- object 1: a code graph whose nodes are binary activation codes and whose arrows are shifts (the de Bruijn
  quantum graph, with Cuntz algebra O_2);
- object 2: the tree of codes and its Cantor set, with weights from the network.

Builds on Bendikov-Grigor'yan-Pittet-Woess (isotropic Markov semigroups on ultrametric spaces) and Khatua-Roy
(path correspondences of quantum graphs).

- Object 1's path correspondences are the levels of object 2.
- A backtracking walk on the code tree exits with exactly the code law (verified, chi2/dof = 1.02).
- Duality (BGPW/Kigami): the induced process on the Cantor set is isotropic, with jump kernel depending on the
  last point of contact and spectrum {1/G(v,o)} from the walk's Green function.
- Weights -> kappa = mu^(1/s): the spectral zeta function has a simple pole at s with residue s/h (h = entropy rate
  of the code process), and the Dixmier trace is integration against mu, so activation means are Dixmier traces
  (verified numerically).
- Run length T(eps) = 2 eps^-2 sum_v G(v,o) e_v, and a two-phase criterion R^2_l > l/L.
- Measured on width-256 networks: weight-chosen code cells explain at most ~9% of output variance (gain <= 1.10x),
  so no gain over Monte Carlo for random He networks; it would pay with internal randomness or low-entropy codes.

Files: `two_objects.pdf` / `.tex`; `code/tree_walk_exit.py`, `code/zeta_dixmier.py`, `code/code_tree_energy.py`
(+ output).
