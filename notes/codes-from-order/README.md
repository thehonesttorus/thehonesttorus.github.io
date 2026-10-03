# Codes from order

Working note VI (self-reviewed draft). It covers the frontier linking expanders, quantum codes and
pseudorandomness, and uses the programme's NCG as a computational instrument (spectral localizer, Anderson-Putnam/Rauzy
complexes, gap labels, Michon trees, Jewett-Krieger hierarchies) on Li-Boyle quasicrystal codes.

Main results:
- Knill-Laflamme for translation-orbit codes: the error functional is the canonical trace of the hull algebra;
  local indistinguishability = unique ergodicity; recoverability = trivial homoclinic relation.
- Entropy price: k <= N h_nu(l) <= N S(rho_l)/l <= N log2 N / l in 1D (tight: (1-o(1)) N log2 N);
  in d dimensions, k r^d <~ N log2 N.
- Exact type-class count (BEST theorem, Kirchhoff determinant of the Rauzy graph), and a Betti bound
  k <= N (b1 - 1) max nu. Substitution and Sturmian statistics give k l = O(N); Li-Boyle's Fibonacci codes
  sit within a small factor of this.
- Sturmian codes: at most 3 entanglement eigenvalues, each in Z + alpha Z (gap labels). Bounded type is
  equivalent to Bellissard-Julien embeddability and to uniformly conditioned local states.
- Inner (amenable/hyperfinite) structure cannot give linear distance; Ta-Shma's construction gives explicit
  hierarchical quantum expanders of size O(N / eps^(2+o(1))), improving note I.

Files:
- `codes_from_order.pdf` / `.tex`: the note.
- `code/necklace_bound.py`: exhaustive necklace type classes vs. the entropy bound; explicit
  Knill-Laflamme checks.
- `code/lb_fib.py`: Li-Boyle Fibonacci approximant codes vs. the entropy and Betti bounds.
- `code/uniform_type.py`: exact BEST/Kirchhoff counts (tightness).
- `code/sturmian.py`: Fibonacci entanglement spectra, gap labels, S(l), l h(l).
