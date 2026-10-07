# Conditional contraction at the chain's third-cumulant readouts

Working note XXX. A kernel-level test of the conditional-contraction development (the supplied `DEVELOPMENT.md`): an
expensive source-to-readout contraction written as an exact expectation over auxiliary sign variables, its three
conditional fillings, a convex target-preserving gauge whose obstruction is a weighted cycle (Hodge) defect, and the
score criterion (E_0 + V) max(0.1, F/B) < E_0 max(0.1, F_0/B). The development asks for one experiment before any
system change: declare a real expensive contraction and its oracles, compare the fillings that can actually be
implemented, calibrate the gauge, and measure the variance budget. This note does that on the chain's real legs
(`est_v35` with the leg dump; official networks 0 and 1, layers 5, 9, 13; `code/cohere.py`, `outputs/cohere_*.txt`).

## 1. The readouts in the development's form

**D3** (the (3,) slice read at each layer) is a cubic contraction over the internal index a = (term, source k, birth
neuron j):

  D3_i ~ sum_a X1_ia X2_ia X3_ia,  terms (A, A, 3 w2 P) and (M, P, P),  M = P s + 3 A e + Z L^T,

with r = 2kn internal terms (10,240 at layer 5, 26,624 at layer 13). These two terms reproduce the chain's D3 at
correlation 0.998 and carry 84-86% of its energy (the rest is the feedback and feed terms).

**D21** is the matrix readout sum_a LA_ia A_ca + LP_ia P_ca (correlation 0.976-0.990 with the chain's D21).

**Oracles.** The legs are transported linearly, A_(l+1) = W diag(w1) A_l, so packets A eps, P eps and their column-
scaled versions (M eps, P (w2 o theta)) cost 2 n^2 per column per layer: the RAW cubic estimator is implementable.
The filled cuts need paired oracles such as (A o P) eps, and the Hadamard product does not commute with transport,
(WA) o (WP) != W (A o P); lifting it needs T (x) T on n^2-sized objects. Under transport only the raw cut is cheap.
The dense legs themselves are exact quadrature with r basis probes, so a packet scheme pays only if its single-probe
relative variance stays below m delta^2 with m < r.

## 2. Measurements (single-probe variance summed over rows, relative to sum_i d_i^2)

| network, layer | r | kappa = abs(d)/sum abs(t) | raw cubic: ungauged / best per-row gauge (Holder) | filled per-row floor | A or BC / B or AC / C or AB: ungauged -> shared convex gauge | cycle loss shared/floor (A, B, C cut) |
|---|---|---|---|---|---|---|
| 0, 5 | 10,240 | 0.040 | 5.3e6 / 1.0e6 | 273 | 1412 -> 729 / 567 -> 424 / 1033 -> 741 | 2.67 / 1.55 / 2.71 |
| 0, 9 | 18,432 | 0.042 | 1.5e7 / 1.5e6 | 254 | 2512 -> 647 / 860 -> 395 / 1508 -> 755 | 2.55 / 1.55 / 2.97 |
| 0, 13 | 26,624 | 0.050 | 2.6e7 / 1.3e6 | 167 | 2837 -> 403 / 961 -> 259 / 1916 -> 517 | 2.42 / 1.55 / 3.10 |
| 1, 5 | 10,240 | 0.042 | 5.3e6 / 1.1e6 | 285 | 1394 -> 758 / 583 -> 442 / 1081 -> 778 | 2.66 / 1.55 / 2.73 |
| 1, 9 | 18,432 | 0.042 | 1.5e7 / 1.6e6 | 263 | 2471 -> 664 / 863 -> 408 / 1648 -> 787 | 2.52 / 1.55 / 2.99 |
| 1, 13 | 26,624 | 0.046 | 3.2e7 / 1.7e6 | 233 | 3570 -> 558 / 1233 -> 361 / 2481 -> 725 | 2.40 / 1.55 / 3.12 |

D21 (sampled entries): kappa 0.023-0.029, filled per-entry floor 520-830 per probe.

**Where the cancellation lives** (`outputs/cohere_where_*.txt`): within one source, over its birth neurons, kappa is
0.06-0.08; across the 2k term-sources it is 0.51-0.58. The newest source carries 13-32% of sum_k abs(d_k). There is
no head: the 16 / 128 / 1024 largest terms of a row carry 4-5% / 17-22% / 51-65% of its absolute mass.

## 3. Reading

- **The development's geometry is confirmed on real data.** The shared gauge is a convex problem that L-BFGS solves
  in a few hundred iterations, and it removes 1.3x-7x of the filled cuts' variance, more at depth. Its cycle loss
  against the per-row optima (which no common gauge can reach) is a stable 1.55 for the B|AC cut and 2.4-3.1 for the
  other two, the same at every layer and on both networks. B|AC, keeping [A, P] as the retained address, is the
  best cut throughout.
- **The floor is set by sign coherence, which no gauge touches.** A positive gauge preserves every t_ia; the
  per-row floor of a filled cut is (sum_a abs(t_ia))^2 / d_i^2 = 1/kappa_i^2. Against exact quadrature (r basis
  probes, zero variance) a cut with single-probe relative variance V breaks even at relative error
  delta* = sqrt(V / r). For the best implementable-in-principle cut, B|AC with the shared gauge, delta* is 0.20-0.21 at
  layer 5, 0.15 at layer 9 and 0.10-0.12 at layer 13 (networks 0, 1): it would match the dense legs only if D3 could
  carry a 10-20% random error, and only if its paired oracle cost no more than a transported leg column. The raw
  cubic, the only cut transport implements, has delta* of 7 to 10 (Holder floor 1.0-1.7 x 10^6 per probe): about 10^3 times
  more probes than r at any useful precision. D21 is worse than D3.
- **The incoherence is quenched.** Per source, D3_i = <A_i^2 o w2, P_i>: a nonnegative weight against P_i, a
  transported weight leg whose signs are those of this network's W. Its net is 6-8% of its absolute mass (about 200
  effective independent terms), spread over the whole birth layer. That is exactly the fluctuation-scale, weight-
  specific content that annealed constructions lose (notes XXVIII-XXIX), and it is why the chain keeps these legs
  dense: there is no address in which a few terms carry the readout.

**Why no address fixes this, for this readout.** D3 is an odd-order form S_b(p_i, p_i, p_i) of the transported
Jacobian rows p_i (the P legs), with the birth source S_b = sum_j (a_b e_j) (x) (a_b e_j) (x) w2_j e_j. Any
representation S_b = sum_m lambda_m (v_m)^(x)3 has terms lambda_m (v_m . p_i)^3 whose signs follow sign(v_m . p_i),
which for quenched random-sign rows is a coin per term. So kappa_i^2 ~ 1/R_eff for every representation, the filled
floor is ~R_eff per probe, and probes beat exact quadrature in that representation only when delta^2 > R_eff / R,
that is, only when the term magnitudes are very uneven. Here R_eff is about 200 per source of n = 1024 birth neurons
(kappa_within 0.07), which is the delta* above. An even-order readout with sign-definite terms (for example a
diagonal weighted by W o W against a positive kernel) is where this development can pay; D3 and D21 are not.

**Verdict against criterion (30).** At the chain's third-cumulant readouts no implementable conditional filling
undercuts exact quadrature: the raw cut loses by about three orders of magnitude, and the filled cuts, which would
break even only at a 10-20% random error in D3, have no paired-factor oracle under transport ((WA) o (WP) !=
W (A o P)). The two obstructions are measured, not assumed. A representation that carried a paired oracle, or an
even-order readout with coherent terms, is what would reopen the question; section 1's arithmetic test is the
first thing such a representation must pass.
