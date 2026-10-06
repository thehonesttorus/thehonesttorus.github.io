# Every computation of the chain, read through the theory

Working note XV. The current best system is the open-source factorized third-cumulant chain in its v29 form
(504aldo, MIT; 1681 lines; raw 2.29e-8 at 0.267 B on official network 0, adjusted 6.1e-9). This note lists what it
computes per layer, what each object is in the language of the earlier notes (fresh weights, the gain channel of E10,
Lyapunov transport, the Hadamard wall), where the approximations sit, and which of the candidate wins survive a cost
count. Units: one unit = 2n^3 = 2.15 GFLOP; the budget is 1024 units; v29 spends 260.

## 1. The per-layer computation

| step | what is computed | cost (dense count; v29 bills ~0.5-0.6x under Strassen) | reading |
|---|---|---|---|
| linear | mu = W mu; W C rides the transport family; C_pre = (W C) W^T by three block products of the 2x2 partition; var, C_off | 1 + 0.75 u | the Gaussian closure's second moment; exact |
| young transports | A_s <- W (w1 A_s), P_s <- W (w1 P_s) for the <= 4 sources younger than 5 layers, as one batched family | 2 u per source | the propagator B = W diag(Phi) applied to the legs: first-order (product-gate) linear response |
| old tier | Qc <- W diag(w1) Qc; dense legs re-formed from the shared basis, A_s = Qc FA_s, P_s = Qc FP_s (tier 2 through Qc U) | 0.37 + 0.75 u per source | the legs confined to the top-384 forward Lyapunov subspace; exact between joins |
| join (one per layer) | one-pass randomized range finder on the weighted leg Gram of joiner + old (sketch = a weight slice), QR, rotation of all old factors into the new basis | ~1.75 + 2 u per join, 12 joins (~40 u, 15% of the bill) | the only lossy step of the old tier: 12 successive truncations |
| dslices | AP = A*P, PP = P*P, MP, LA, LP (Hadamard products of legs); D3 by triple products; D21 = sum_s LA_s A_s^T + LP_s P_s^T (hub family) + thin terms | 2 u per young source (hub), 0.75 u per old source (through the factors) | the (3,) and (2,1) slices: the identification (i,i) that no basis survives (Khatri-Rao rank r^2) |
| kappa_4 regen | dG = (W*W)(g_prev - lam var_prev) + lam var; g4row = 2 dG; wk4m, wk431 from dG and lam C_off; adaptive lam from mean(dG)/mean(var) | n^2 | memoryless fourth: the annealed (isotropic) transport of the post-ReLU (2,2) content plus a fitted counterterm |
| wick + nonlinear terms | W_all (n x 21 Hermite/sigma table); ~20 Hadamard products of {c_off, d21, wk4m, d3, g4, wk431}; two fused einsums -> pk11, pk21, pk22, pk1..pk4; the moment-to-cumulant table -> K11, K21, K22, K1..K4 | ~150 n^2 (0.1 u) | the Gaussian + Edgeworth closure of the post-ReLU moments to second order in the carried cumulants, Mehler order 2 |
| birth | a_b = w1 * C_off (A leg), P = I; S21 = K21 - rep21; S_sep = e_b C_off w1 rides on A; Rres = S21 - S_sep -> rank 16 (two-pass range finder); D21 feedback -> rank 16 thin legs; K4 -> K3 feed | n^2 r | the fresh third cumulant as the star core Sym(C (x) I (x) C diag(w2)); the slice residual that the star does not hold |

Everything else (K4_vec, the harmonic projection constants cA, cI, the online correction rider which is off) is O(n).

## 2. Structural identities the code does not state

- **One propagator per source.** A leg born as a_b = w1 * C_off and a P leg born as I are transported by the same
  B's, so at every later layer A_s = P_s C_b exactly, with P_s = B_l ... B_{b+1}. The whole source is
  (P_s (x) P_s (x) P_s) applied to a frozen birth core T_b = Sym(C_b (x) I (x) C_b diag(w2)) plus thin terms: the chain
  is the sum over birth layers of propagated star cores.
- **Nesting.** P_{b'}(l) = P_b(l) P_{b'}(b) for b' < b: every older source's legs are the younger source's propagator
  times its own frozen legs. This is why one shared basis works for all old sources (their top left singular
  subspace is the propagator's, which the most recent layers fix), and why it cannot save more than a factor two:
  forming a leg from the basis costs n^2 r, the same as transporting it, and the (i,i) identification of the slices
  needs the formed leg.
- **Cost law.** The bill is the CP rank of the summed cores (n hubs per source) times the number of live sources,
  because the (2,1) slice of a hub term is a Hadamard product of two legs. Horner-style re-association of the nested
  propagators, adjoint (backward) transport of the readout, and a merged Tucker core all land on the same count or
  worse (the core costs n r^3 per layer to read). This is the arithmetic form of "quenched content has no cheap
  carrier" (note XIV).
- **The old tier is exact between joins.** Qc <- B Qc with static factors reproduces B A_s exactly; accuracy is lost
  only at the twelve joins, where joiner + old are re-projected onto rank 384. The join is therefore the place where
  the basis quality matters, and the sketch it starts from is a free choice.

## 3. Where the approximations are, with the theory's verdict

| approximation | size | verdict |
|---|---|---|
| product gate Phi_a Phi_b Phi_c on the all-distinct entries of the transport | 1.5-2.4% of the gate per layer, zero mean (E3); compounds along a source's life | the residual's identified ingredient; its coherent part vanishes for distinct indices (the only surviving pairing is i = j, a = b, which is the slice the pair program already does exactly), so no O(n^2) correction exists; the exact fix is n^4 |
| Mehler order 2 in the slice programs | 0.05% (frontier tests) | done |
| memoryless kappa_4 with a fitted lambda and a transported diagonal | 90% of the full channel; the exact core is -9% raw at 2-4x the cost | E10/E11: the fitted diagonal is a sharp counterterm of the chain's own truncation (3 g sigma^4 double counts); the fourth order is slaved to the third, so memoryless-from-kappa_3 is the physical model |
| rank-16 residual leg and rank-16 feedback | higher rank is worse raw (16/32/64 -> 2.35/2.39/2.41e-8) | the slice residual S21 - S_sep contains the leg transport's error at the slice level; the truncation is regularisation, not a loss of physics |
| twelve join truncations at rank 384 / 224 | +3.4% / +1.2% raw | tested below: the sketch |
| harmonic projection of the (2,2) slice to a diagonal core | the quenched (2,2) content is dropped each layer | the augmented chain carries it at 7 units per source-layer for -9% raw: not affordable |

## 4. Candidate wins, costed

- **Warm-started joins (zero cost).** Replace the weight-slice sketch of the join's range finder by the transported
  previous basis Qp (and the nested sketch by the previous sub-basis U). Lyapunov continuity says the target subspace
  is mostly span(Qp); one pass from there should beat one pass from a random slice. Measured below.
- **A cheaper join sketch.** The Gram products of the finder (4 n^2 r) are 60% of the join; a sketch [Qp | A_j Om_s |
  P_j Om_s] with a 32-column Om_s would cut them to 0.13 u, but the joiner needs ~3n/8 new directions (the rank law),
  and a 64-direction sketch of it fails unless its content lies almost entirely in span(Qp). Depends on the first
  test; worth at most 20 units (8% of the bill).
- **Strassen level 6.** Per product, the next level saves 1/8 of the matmul count and adds 74 n^2 of block
  additions: net -6% of the matmul bill, -5% of the total, if the residual clock allows and float32 holds. Measured
  below.
- **A g-based adaptive lambda.** E10's channel coefficient g (from the chain's own D3, free) instead of
  mean(dG)/mean(var); the current rule is worth 0.6-1.3%, so this is bounded by about 1%.
- **Coherent top-up of D3/D21.** The chain's D3 coefficient is 1.5-2% below the measured channel on both networks;
  the frontier tests bound any gain/scale correction of the residual at 2.4% of the MSE, and the derived recursion
  costs ~0.3 u a layer in its annealed form: neutral.
- **Not wins:** symmetric or dtype pricing (exhausted, F77); a cheaper exact algorithm for the hub (none: the
  contraction is a non-symmetric matmul on legs); dropping or merging hubs (flat importance); a Tucker core of the
  old tier (n r^3 per layer); Cholesky-aliased Gram for C_pre (1.7 u against 0.95 u).

## 5. Results (official networks 0 and 1; `code/est_v29.py` = v29 with the variants behind environment flags,
`code/run_v29.py`; raw = final-layer MSE, adjusted = raw x C/B under the Phase 2 rule; outputs in `outputs/`)

| variant | net 0 raw | C/B | net 0 adjusted | net 1 raw | C/B | net 1 adjusted |
|---|---|---|---|---|---|---|
| v29 as shipped | 2.290e-8 | 0.2671 | 6.116e-9 | 2.086e-8 | 0.2671 | 5.571e-9 |
| warm join (sketch = transported old basis) | 2.222e-8 | 0.2671 | 5.934e-9 | 2.072e-8 | 0.2671 | 5.535e-9 |
| warm join + warm nested sketch | 2.242e-8 | 0.2671 | 5.987e-9 | 2.070e-8 | 0.2671 | 5.527e-9 |
| warm join + warm feedback sketch | 2.183e-8 | 0.2671 | 5.830e-9 | 2.063e-8 | 0.2671 | 5.510e-9 |
| warm join + warm residual sketch | 2.255e-8 | 0.2673 | 6.027e-9 | | | |
| all four warm sketches | 2.221e-8 | 0.2673 | 5.936e-9 | | | |
| two cold passes (v29 + QPASS 2) | 2.259e-8 | 0.2873 | 6.490e-9 | | | |
| two warm passes | 2.233e-8 | 0.2873 | 6.415e-9 | | | |
| Strassen level 6, leaf 16 | 2.293e-8 | 0.2624 | 6.019e-9 | | | |
| **warm join + warm feedback + level 6** | **2.183e-8** | **0.2624** | **5.728e-9** | **2.066e-8** | **0.2624** | **5.421e-9** |

- The warm-started join is strictly better than a second cold pass (which costs 7.6% more and still loses), and a
  second warm pass adds nothing: the join's basis is saturated at zero cost, as the Lyapunov argument predicts (the
  target subspace is the transported previous one plus what the joiner brings).
- The feedback finder benefits from the same warm start (the (2,1) feedback basis of consecutive layers is related
  by transport); the slice-residual finder does not (warm residual is worse than cold): Rres is newly generated
  content and transport error, not a propagated object. That is a measurement of what the residual leg holds.
- Level 6 saves 1.8% of the bill at an unchanged raw (2.293 against 2.290e-8), less than the per-product count
  suggests because the block additions grow; the residual clock was not measured here and must be checked on the
  grader's machine before shipping.
- Combined, zero-cost warm starts plus level 6 take the adjusted score down 6.3% (network 0) and 2.7% (network 1):
  real, derived, small. On the public board (5.40e-9) that is about 5.1-5.2e-9, rank unchanged. The 1.7x to the
  leaders is not in this chain's computations; it is in its representation (section 2).
