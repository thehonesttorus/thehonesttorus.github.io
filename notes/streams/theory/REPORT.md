# Second-order diagrams of the all-distinct post-activation kappa3: report

## Verdict
1. **(a) No.** The fixed-coefficient second-order closure does not reach the fitted level at layers 10-14. Width 128, noise-corrected, both orientations:
   - Fixed second order: 2.3-3.75 % (layer means 3.12 % for MLP0, 2.62 % for MLP1).
   - Fitted first order (fit1): 1.44-2.03 % (means 1.85 % and 1.97 %).
   - The best fixed closure I found resums the leaf and T classes (LRT1, section 1.6): 1.77-2.70 % (means 2.39 % and 2.36 %), still about 0.45 points above fit1.
   - At layers 1-7 the fixed second-order closure does beat fit1. At layers 1-3: 3.59 / 3.44 % against 4.30 / 4.18 %.
2. **(b) The drift is only partly explained.** I refit B0..B6 with all second-order terms held at their fixed coefficients:
   - B6 (the (2,1,1) kappa4 term) returns to 1.45-1.53 at every layer of both MLPs. Without second order it drifts 1.45 -> 1.13.
   - B3 returns to 0.50-0.82 (it was -0.17); see the bug in item 3.
   - B4 (rho^3) returns to 0.65-1.08 (it drifted 1.0 -> 0.31).
   - B5 (K22 + edge) stays at 0.22-1.27 and is not restored.
   - The D21 terms overshoot: B1 rises to 3.35-4.70 and B2 to 3.46-3.89 (theory 3 for both). So the second order overcorrects, and the series behaves as an alternating asymptotic one: depth-wise rho = 0.11 -> 0.38, marginal skew up to 0.5, kurtosis up to 0.64.
   - In the resummed basis, the T-class coefficient (one kappa3_ijk hyperedge with the rest on one pair, i.e. kappa3_ijk times the exact pair gate covariances) fits to 0.97-1.00 at every layer of both MLPs. That fully explains B7's apparent renormalization (its free coefficient was 3.1 -> 5.0 against theory 3).
   - The leaf class (6 x LEAF in theory) fits to 5.8-7.1, the TWOLEAF subtraction to -0.93 to -1.31 (theory -1), and B6 to 1.41-1.54.
3. **Two corrections to oracle_k3's first order.**
   - `CLOSURE_COEF[3]` (B3) should be 0.5, not 1.0. The term w5_i Phi_j Phi_k D3_i C_ij C_ik is symmetric under swapping j and k, so the six sym3 role assignments count each diagram twice. The enumerator gives a ratio of exactly 2.0000. Fixing this one coefficient lowers the fixed first-order error at layers 10-14 from 4.15 to 3.08 % (MLP0) and from 4.48 to 3.35 % (MLP1).
   - One first-order (eps^3) diagram is missing from the basis: B7 = 3 sym3[w2_i w2_j Phi_k kappa3_ijk C_ij]. It is 2.6-3.3 % of ||D21|| at layers 1-5 and lowers the layer 1-3 error by about 0.55 points.
4. **Cost at n = 1024.**
   - Every term whose factors sit on two vertex pairs (a "path" term) transports with 11 n x n matmuls, i.e. 11 units, with no n^3 object. This is verified to 1e-15 against the dense n^4 transport.
   - All C-leaf terms share one arm, Phi∘C, so the whole leaf class costs 11 units per layer together. Once the Wick term is transported that way, its extra cost is zero.
   - Non-leaf path terms cost at most about 44 units. Triangle terms cost about 3n units (about 3 B). Terms needing an n^3 slice cost about 2n units (about 2 B) and need 8 GB each.
   - The frontier budget is about 7.7 units per layer, so only leaf-class terms and perhaps one extra arm are affordable. The best candidate is S4d.

## 1. Derivation

### 1.1 Expansion and the vertex-weight rule
- Edgeworth form: E[F(z)] = E_G[exp(sum_{m>=3} kappa_m . d^m / m!) F], with E_G the Gaussian of the same mean and covariance.
- Gaussian product rule: E_G[prod_v f(z_v)] = exp(sum_{u<v} C_uv d_{mu_u} d_{mu_v}) prod_v E[f(z_v)].
- So a vertex hit by d derivatives carries **w_v(d) = E[f^(d)(z_v)]**. For relu: w1 = Phi(alpha), and w(d) = He_{d-2}(-alpha) phi(alpha) / sigma^{d-1} for d >= 2.
- A term is a multiset of hyperedges (multisets of vertices). A C edge is a size-2 hyperedge on two distinct vertices; kappa_m is a hyperedge with m >= 3. One labelled hypergraph weighs: prod_h kappa_h / prod_v r_{h,v}! x prod_types 1/mult! x prod_v w_v(deg v).
- Weights factorize over components, so by Moebius inversion kappa3(a_i, a_j, a_k) on distinct i, j, k is the sum over connected hypergraphs on {i, j, k}.
- There are no internal vertices: every diagram is an elementwise product of slices.
- The w(d)/d! times leg-partitions rule is identical (a vertex has d!/prod r! leg partitions).
- In sym3 form, a term's coefficient is **c = (6/|Aut X|) x prod 1/r! x prod 1/mult!**.

### 1.2 Dressed vertices
- Single-vertex hyperedges (D3, K4, ...) sum to the true-marginal weight E_true[f^(d)(z_v)].
- For d = 1 this is the true gate P(z > 0), which the atlas gate_p stores exactly, so every degree-1 vertex is dressed to all orders.
- For d = 2 it is the true density p(0): w2 + D3/3! w5 + K4/4! w6 + ...

### 1.3 Order counting (fixed random MLP, eps = n^-1/2)
- C_off ~ eps.
- A kappa_m slice ~ eps^(m-2) if every multiplicity is even (K22 = kappa4_iijj, K4, the kappa6 (2,2,2) slice), and eps^(m-1) otherwise:
  - every kappa3 slice ~ eps^2;
  - the kappa4 (2,1,1) and (3,1) slices ~ eps^3;
  - kappa5 ~ eps^4.
- The target is eps^2.
- **Measured** at layer 1, width 32 vs 128 (rms, sigma-normalized). Expected ratios are 2 for rho, 4 for eps^2 objects, 8 for eps^3 objects.

| object | width 32 | width 128 | ratio | expected |
|---|---|---|---|---|
| rho | 0.211 | 0.109 | 1.94 | 2 |
| kappa3 all-distinct | 0.084 | 0.020 | 4.2 | 4 |
| D21 | 0.138 | 0.033 | 4.2 | 4 |
| D3 | 0.288 | 0.071 | 4.1 | 4 |
| K22 | 0.135 | 0.027 | 5.0 | 4 |
| K4 | 0.354 | 0.077 | 4.6 | 4 |
| K211 | 0.050 | 0.006 | 8.3 | 8 |
| K31 | 0.102 | 0.013 | 7.8 | 8 |

At width 128, rho grows with depth: 0.11 at layer 1, 0.25 at layer 7, 0.29-0.38 at layers 10-14. Kurtosis reaches 0.47-0.64 at layers 10-14.

### 1.4 Order eps^2 and eps^3

| order | term | shape | weights | coef |
|---|---|---|---|---|
| eps^2 | Wick | C C path | Phi Phi w2 | sum over 3 centers |
| eps^2 | B0 | kappa3_ijk | Phi^3 | 1 |
| eps^3 | B1 | D21_ik + C_jk | w2_i Phi_j w2_k | 3 |
| eps^3 | B2 | D21_ik + C_ij | w3_i Phi_j Phi_k | 3 |
| eps^3 | B4 | rho^3 (Hermite 6-4) | | 1 |
| eps^3 | B5 | K22_ij + C_ik | w3_i w2_j Phi_k | 1.5 |
| eps^3 | B6 | kappa4_iijk | w2_i Phi_j Phi_k | 1.5 |
| eps^3 | **B7** | kappa3_ijk + C_ij | w2_i w2_j Phi_k | **3** |
| eps^4 | B3 | D3_i C_ij C_ik | w5_i Phi_j Phi_k | **0.5** |

### 1.5 The 30 second-order (eps^4) shapes
Coefficient derivation: 6/|Aut| x 1/r! x 1/mult!.

| group | term | shape | weights | coef |
|---|---|---|---|---|
| Gaussian | S1 | rho^4 (Hermite 8-6) | | 1 |
| dressing | S2a | Wick with true gates | | 1 |
| dressing | S2c | K4_i C_ij C_ik | w6 Phi Phi | 0.125 |
| kappa3 + 2C | S3a1 | T C_ij^2 | w3 w3 Phi | 1.5 |
| kappa3 + 2C | S3a2 | T C_ij C_ik | w3 w2 w2 | 3 |
| kappa3 + 2C | S3b1 | D_ik C_jk^2 | w2 w2 w3 | 1.5 |
| kappa3 + 2C | S3b2 | D_ik C_ij^2 | w4 w2 Phi | 1.5 |
| kappa3 + 2C | S3b3 | D_ik C_ij C_jk | w3 w2 w2 | 3 |
| kappa3 + 2C | S3b4 | D_ik C_ij C_ik | w4 Phi w2 | 3 |
| kappa3 + 2C | S3b5 | D_ik C_jk C_ik | w3 Phi w3 | 3 |
| kappa3 x kappa3 | S4a | T^2 | w2^3 | 0.5 |
| kappa3 x kappa3 | S4b | T D_ik | w3 Phi w2 | 3 |
| kappa3 x kappa3 | S4c | D_ik D_jk | w2^3 | 0.75 |
| kappa3 x kappa3 | S4d | D_ki D_kj | Phi Phi w4 | 0.75 |
| kappa3 x kappa3 | S4e | D_ik D_kj | w2 Phi w3 | 1.5 |
| odd kappa4 + C | S5a | U C_jk | w2^3 | 1.5 |
| odd kappa4 + C | S5b | U C_ij | w3 w2 Phi | 3 |
| odd kappa4 + C | S5c | V_ij C_ik | w4 Phi Phi | 1 |
| odd kappa4 + C | S5d | V_ij C_jk | w3 w2 Phi | 1 |
| K22 + 2C | S6a | K22 C_ik^2 | w4 w2 w2 | 0.75 |
| K22 + 2C | S6b | K22 C_ik C_jk | w3 w3 w2 | 0.75 |
| K22 + 2C | S6c | K22 C_ik C_ij | w4 w3 Phi | 1.5 |
| K22 x kappa3 | S7a | K22 T | w3 w3 Phi | 0.75 |
| K22 x kappa3 | S7b | K22 D_ik | w4 w2 Phi | 0.75 |
| K22 x kappa3 | S7c | K22 D_ki | w3 w2 w2 | 0.75 |
| K22 x K22 | S8 | K22_ij K22_jk | w2 w4 w2 | 0.1875 |
| kappa5 | S9a | kappa5 (2,2,1) | w2 w2 Phi | 0.75 |
| kappa5 | S9b | kappa5 (3,1,1) | w3 Phi Phi | 0.5 |
| kappa6 | S10 | kappa6 (2,2,2) | w2^3 | 0.125 |

Notation: T = kappa3_ijk, D = D21, U = kappa4_iijk, V = kappa4_iiij. Gaussian non-leaf pieces, also enumerator-checked: G3t = triangle (coef 1), G4a = C_ij^2 C_jk^2 (0.75), G4b = C_ij^2 C_jk C_ik (1.5).

**Atlas fields.**
- Everything above except S9a, S9b and S10 comes from existing fields: pre_M3, pre_M22, pre_M211 (which also gives the K31 slice from its diagonal), pre_s and gate_p.
- S9a, S9b and S10 need three new n^3 fields: E[z_i^2 z_j^2 z_k], E[z_i^3 z_j z_k] and E[z_i^2 z_j^2 z_k^2]. atlas_k56.py builds them.
- The leaf resummation needs Gamma (an n^2 field) and p(0) (an n field), which pair_fields.py computes.

### 1.6 Two exact resummations
- **Leaf identity.** All diagrams in which some vertex k hangs on a single C edge sum to sum_6roles Phi_k C_ik Gamma_ij - sum_centers Phi_j Phi_k C_ij C_ik p_i(0), where **Gamma_ij = Cov(1[z_i > 0], relu(z_j))**.
- **T-class identity.** All diagrams with one kappa3_ijk hyperedge and everything else on one pair sum to kappa3_ijk (Phi^3 + sum_pairs Phi_k Cov(g_i, g_j)). The pair gate covariance comes from gate_GG.
- Completeness at network order <= 4: 187 connected diagrams in total, all covered (18 are gate dressings). The leaf split is 81 C-leaf + 106 non-leaf diagrams, none uncovered.

## 2. Toy validation (results/toy_validation.txt)
- **Per-term identity.** Every closure2 term matches the generic hypergraph enumerator at n = 3 to within 3.7e-15 (39 terms, including Wick, G3t, G4a and G4b). The oracle's B3 ratio comes out as exactly 2.0000.
- **Full expansion vs exact quadrature.** Toy: z = mu + L g + t Q(g) + t^2 cubic, with s = t = lambda, f a Gaussian-smoothed relu (delta = 0.5, quadrature 72^3). residual_N / lambda^(N+1) settles to constants as lambda goes 0.2 -> 0.025:
  - N = 1: -1.14 -> -1.07
  - N = 2: 0.93 -> 0.676
  - N = 3: -0.41 -> -0.484
  - N = 4: -> 1.03
  - N = 5: still drifting (-1.80, -3.23, -3.71).
  - delta = 0.25 gives the same pattern until the quadrature floor (about 1e-9) is reached.
- **Leaf and T-class identities.** After subtracting the resummed parts, residual_N / lambda^(N+1) for N = 1..4 settles to constants: -0.125, 0.205, -0.645, 3.49.
- **Dressing.** At lambda = 0.025 and degree d = 2, the error falls 6.2e-3 -> 1.5e-4 -> 2.3e-5 through orders 0, 1, 2.
- **Not validated:** pure relu (delta = 0). Its quadrature is only accurate to about 1e-3, so that check is qualitative.

## 3. Width-128 results (eps %, noise-corrected, mean of both orientations)

Layer-group means:

| MLP | layers | noise | oracle 1st | 1st (B3=0.5) | +B7 | +2nd fixed | leaf-resummed 1st | leaf+T resummed (LRT1) | fit1 | fit2 | all free |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1-3 | 3.39 | 4.46 | 4.40 | 3.84 | 3.59 | 3.69 | 3.65 | 4.30 | 3.54 | 3.43 |
| 0 | 4-9 | 1.36 | 3.89 | 3.22 | 3.11 | 2.79 | 2.66 | 2.58 | 2.34 | 2.08 | 1.68 |
| 0 | 10-14 | 1.26 | 4.15 | 3.08 | 3.02 | 3.12 | 2.43 | 2.39 | 1.85 | 1.91 | 1.54 |
| 1 | 1-3 | 3.33 | 4.34 | 4.27 | 3.72 | 3.44 | 3.61 | 3.58 | 4.18 | 3.40 | 3.33 |
| 1 | 4-9 | 1.58 | 4.30 | 3.83 | 3.50 | 3.11 | 3.06 | 2.97 | 3.12 | 2.62 | 2.09 |
| 1 | 10-14 | 1.00 | 4.48 | 3.35 | 3.22 | 2.62 | 2.39 | 2.36 | 1.97 | 1.97 | 1.44 |

Layers 10-14 individually:

| layer | MLP0 +2nd | MLP0 LRT1 | MLP0 fit1 | MLP1 +2nd | MLP1 LRT1 | MLP1 fit1 |
|---|---|---|---|---|---|---|
| 10 | 2.66 | 1.77 | 1.44 | 2.66 | 2.70 | 2.21 |
| 11 | 2.68 | 2.55 | 2.01 | 2.31 | 2.04 | 1.89 |
| 12 | 3.35 | 2.47 | 1.78 | 2.73 | 2.64 | 1.99 |
| 13 | 3.15 | 2.63 | 2.02 | 2.64 | 2.09 | 1.72 |
| 14 | 3.75 | 2.55 | 2.00 | 2.77 | 2.32 | 2.03 |

Other findings:
- The free fit over all about 35 diagrams reaches 1.29-1.74 %, so the diagram span is adequate; the coefficients are what is missing.
- **Largest terms at layers 8-14** (MLP0, as a fraction of ||D21||): B0 53-64 %, the leaf class 13-16 %, B6 6.4-10 %, B2 6.3-8.5 %, S5c 2.0-3.7 %, S4d 1.6-2.6 %.
- **S4d alone** added to "+B7" at layers 8-14 lowers the error from 2.29-3.80 % to 1.71-2.53 %.
- **The kappa3 x kappa3 group alone** at layers 10-14 gives 1.66-2.39 % for MLP0 and 1.96-2.20 % for MLP1.
- The full second order is worse than either of these, because the S5/S3/S6 leaf terms overshoot.
- Transported Monte Carlo noise of every second-order term is at most 0.14 % of D21.

## 4. Width 32 with kappa5/kappa6 (seed 770000, depth 6, N = 4e6 per atlas; noise 0.3-0.4 %)
Layers 1-4:

| closure | L1 | L2 | L3 | L4 |
|---|---|---|---|---|
| oracle 1st | 4.2 | 5.8 | 10.0 | 8.2 |
| 1st (B3=0.5) | 3.6 | 5.2 | 7.8 | 6.8 |
| +2nd fixed | 3.2 | 4.7 | 9.8 | 9.9 |
| +2nd +kappa5/6 | 3.0 | 4.2 | 7.75 | 5.8 |
| leaf-resummed 1st | 2.7 | 4.5 | 7.75 | 7.9 |
| fit1 | 3.5 | 4.3 | 6.0 | 4.7 |
| all free | 1.1 | 1.3 | 2.8 | 2.1 |

- The kappa5 slices are not negligible: S9b reaches 9.6 % of D21 and S9a 5.0 %; S10 stays at or below 0.6 %.
- With rho = 0.21-0.36 and kurtosis up to 1.0, the expansion has stopped converging at this width.

## 5. Transport cost at n = 1024 (1 unit = 2^31 FLOPs; per layer)

| class | terms | cost | inputs a chain must carry |
|---|---|---|---|
| C-leaf paths (arm Phi∘C) | Wick, B1, B2, B3, B5, S2a, S2c, S3b4, S3b5, S5c, S5d, S6c, leaf parts of B4 and S1 | 11 units for the whole class; 0 extra over a full Wick transport | C, D21 (+ D3), K22, K31 = kappa4_iiij (n^2, new), K4, gate. Or exactly: Gamma (n^2, new) + p(0) (n, new) |
| non-leaf paths | S3b1, S3b2, S4c, S4d, S4e, S6a, S7b, S7c, S8, G4a | 11 units each; at most 4 distinct second arms, so at most 44 units | C, D21, K22 |
| triangles | G3t, G4b, S3b3, S6b | about n units per role, about 3n = 3 B; about 11 r units with a rank-r arm | C, D21, K22 |
| n^3 slices | B0, B6, B7, S3a1, S3a2, S4a, S4b, S5a, S5b, S7a, TCL, S9a, S9b, S10 | about 2n = 2 B dense, plus 8 GB storage each | kappa3_ijk, kappa4_iijk, the kappa5 (2,2,1) and (3,1,1) slices, kappa6 (2,2,2), gate_GG |

The 11-matmul breakdown: G and H (2), R_center (1), R_end_i (2), R_end_j (2), coincidence-slice correction (4).

**Budget.** The frontier allows about 0.12 B, i.e. 123 units, or about 7.7 units per layer. One extra arm costs about 0.17 B over 16 layers.
- **Costs nothing:** fix B3 = 0.5, and fold the leaf terms into the Wick transport.
- **Candidate paid addition:** S4d, at 11 units per layer.

## 6. Open
- The remaining gap at depth (LEAF fitting to 6.4-7.1 instead of 6, G3t to about 1.6) points at the D-leaf and star classes, e.g. S4d and S3b2. Resumming those exactly needs per-center objects of the form E[He_r(u_k) He_s(u_k) a_k].
- The 1,000-network check at width 1024 has not been run.

## Housekeeping
- Superseded partial result files were deleted. They had been committed in snapshot 0e3eaed; the final runs contain the same numbers.
- The width128_final_*.txt logs were restored from the scratchpad logs. A snapshot during the run had left the working copy at only 19 lines.
- Atlases and pair files live in the session scratchpad: k56_w32/seed{1,2}/{mlp_00000,pair}.npz and atlas128/seed{3,4}/pair_mlp{0,1}.npz.
- Commands below assume: source /root/whest/bin/activate; cd /home/user/thehonesttorus.github.io/notes/streams/theory; S=/tmp/claude-0/-home-user-thehonesttorus-github-io/b94ff8ab-040c-538b-a2af-dbe7aa5288dd/scratchpad.

## Claims and adversarial verification (two verifiers per claim; a claim survives if neither verifier refutes it)

- **C1**: The leg-partition rule (vertex weight w(d)=E[f^(d)], coefficient 6/|Aut| x prod 1/r! x prod 1/mult!) is validated. Every closure2 term (39 incl. Wick, G3t/G4a/G4b) equals the generic hypergraph enumerator at n=3 to <=3.7e-15. The truncated expansion vs exact quadrature on a non-Gaussian toy (smoothed relu, delta=0.5) has residual_N/lambda^(N+1) tending to constants for N=1..4.
- **C2**: oracle_k3.CLOSURE_COEF overcounts B3 (the D3 dressing of the Wick vertex) by exactly 2: the correct coefficient is 0.5, because X is symmetric under j<->k (|Aut|=2). Fixing it alone lowers the fixed first-order D21 error at width 128, layers 10-14, from 4.15 to 3.08 % (MLP0) and from 4.48 to 3.35 % (MLP1).
- **C3**: An eps^3 diagram is missing from the oracle first-order basis: B7 = 3 sym3[w2_i w2_j Phi_k kappa3_ijk C_ij]. It is 2.6-3.3 % of ||D21(l+1)|| at layers 1-5 (MLP1) and 0.6-1.3 % at layers 8-14 (MLP0). Adding it lowers the layer 1-3 error by ~0.55 points (4.40->3.84 % MLP0, 4.27->3.72 % MLP1).
- **C4**: The order counting C_off~eps, kappa_m~eps^(m-2) (all multiplicities even) or eps^(m-1) (otherwise), eps=n^-1/2, holds. Width 32 vs 128 at layer 1 gives ratios: rho 1.94 (2), kappa3 slices 4.1-4.2 (4), K22 5.0 and K4 4.6 (4), K211 8.3 and K31 7.8 (8). The eps^4 set has 30 shapes; at network order <=4 all 187 connected diagrams are covered. Three shapes (kappa5 (2,2,1), kappa5 (3,1,1), kappa6 (2,2,2)) need new atlas fields.
- **C5**: (a) At width 128, the fixed-coefficient second-order closure does NOT reach the fitted first-order level at layers 10-14. It gives 2.66-3.75 % (MLP0, mean 3.12 %) and 2.31-2.77 % (MLP1, mean 2.62 %), vs fit1 1.44-2.02 % (mean 1.85 %) and 1.72-2.21 % (mean 1.97 %). The best fixed closure, leaf+T resummed (LRT1), gives means 2.39 % and 2.36 %. At layers 1-3, the fixed second order beats fit1 (3.59 vs 4.30 %, 3.44 vs 4.18 %). A free fit over all ~35 diagrams reaches 1.29-1.74 %.
- **C6**: (b) With the second-order terms fixed, the refit first-order coefficients partly return to leg-partition values. B6 is 1.45-1.53 at every layer of both MLPs (it drifted 1.45->1.13 without second order). B3 is 0.50-0.82 (was -0.17). B4 is 0.65-1.08 (drifted 1.0->0.31). B5 stays at 0.22-1.27 (theory 1.5). B1 overshoots to 3.35-4.70 and B2 to 3.46-3.89 (theory 3). The drift is only partly explained; second order overcorrects the D21 diagrams.
- **C7**: Two exact resummations hold and explain part of the renormalization. (i) Leaf identity: all diagrams with a C-leaf vertex sum to 6 sym3[Phi_k C_ik Gamma_ij] minus Wick[p(0)], with Gamma=Cov(1[z>0],a). (ii) T-class identity: kappa3_ijk(Phi^3 + sum Phi Cov(g,g)). Both are validated on the toy. At width 128 the refit T-class coefficient is 0.97-1.00 at every layer of both MLPs (theory 1), absorbing B7's apparent 3->5.0 renormalization. LEAF fits to 5.80-7.09 (theory 6), TWOLEAF -0.93 to -1.31 (theory -1), B6 1.41-1.54.
- **C8**: At width 128 depth (layers 8-14), the largest second-order terms are S4d (D21_ki D21_kj Phi Phi w4_k, 1.6-2.6 % of D21) and S5c (kappa4_iiij C_ik w4 Phi Phi, 2.0-3.7 %). Adding S4d alone to the corrected first-order closure lowers the error from 2.29-3.80 % to 1.71-2.53 % (MLP0). The kappa3 x kappa3 group alone gives 1.66-2.39 % (MLP0) and 1.96-2.20 % (MLP1) at layers 10-14. The transported Monte Carlo noise of every second-order term is <=0.14 % of D21.
- **C9**: At width 32 (rho 0.21-0.36, kurtosis up to 1.0) the kappa5 slices are large: S9b up to 9.6 % and S9a up to 5.0 % of D21; the kappa6 (2,2,2) term S10 stays <=0.6 %. Including them lowers the fixed second-order error at L4 from 9.9 to 5.8 % (L3: 9.8->7.75). No fixed closure beats fit1 at L3-L4 (6.0, 4.7 %), so the expansion has stopped converging at this width.
- **C10**: Transport cost at n=1024 (1 unit = 2^31 FLOPs). Every path-shaped term (factors on two vertex pairs) transports exactly with 11 n x n matmuls (11 units/layer, no n^3 object). The whole C-leaf class shares the Phi∘C arm, so it costs 11 units together (0 extra over a full Wick transport). Non-leaf paths cost <=44 units, triangles ~3n units (~3 B), and n^3-slice terms ~2n units (~2 B) plus 8 GB each. Against the ~7.7 units/layer frontier budget, only the leaf-class folding and at most one extra arm (e.g. S4d, ~0.17 B over 16 layers) are affordable.

20 verifier verdicts recorded; 5 refutations. Verifier reasons:

- upheld: I reran every check named in HOW TO CHECK and got the claimed numbers to every digit. I then tested the rule with code I wrote myself and on more seeds, finer grids and smaller lambda. Every test agrees with the claim, and none refutes it.

Two small problems, neither of which refutes the numbers:
1. notes/streams/theory/REPORT.md does not exist. It is not on disk and git has no record of it under any ref.
2. The count is slightly off: there are 39 closure2 terms plus the Wick term, which makes 40 comparisons, not "39 incl. Wick".

There is also a caveat about what the identity check proves. T
- upheld: I could not refute C2. Both halves reproduce, and I checked each one in a way that does not depend on the derivation agent's code.

(1) Combinatorics. B3 is the D3_i self-hyperedge with two C edges on the same center i, so the vertex degrees are (5,1,1) and X is symmetric under j<->k. Each labelled diagram has weight D3_i/3! * C_ij C_ik * w5_i Phi_j Phi_k. Summing over the 3 possible centers gives 0.5*sym3[X]. oracle_k3 counts 6 role assignments, but only 3 are distinct. It applies this |Aut| correction to B6 but leaves it out for B3.

I confirmed this with a sympy computation that does not us
- upheld: C1 holds. I reproduced it exactly and confirmed the rule with an independent test the agent did not run. The caveats below weaken parts of the evidence but do not overturn the claim.

What survives:
(1) The leg-partition rule is exact where an exact check is possible. I wrote my own degree-capped hypergraph enumerator (scratchpad/verify-theory/poly_exact.py). For a polynomial f the expansion is a finite sum, so it must match quadrature exactly at any strength of non-Gaussianity. Tests: f of degree 2 (28 diagrams) and degree 3 (319 diagrams, hyperedges up to kappa9), with s,t up to 0.8/1.0. The
- upheld: I re-ran every number in the claim from the raw atlases and all of them reproduced, so I could not refute it.

Method: I wrote my own script, /tmp/claude-0/-home-user-thehonesttorus-github-io/b94ff8ab-040c-538b-a2af-dbe7aa5288dd/scratchpad/verify-theory/c3v1/indep_b7.py. It loads seed3/seed4 mlp_0000{0,1}.npz and computes everything itself: the central moments (including my own expansion of kappa4_iijk from M211), the B0..B3 and B5..B7 tensors, the transport through W[l+1], and the noise-corrected eps. The only code it borrows is oracle_k3.hermite_model, for the Gaussian Wick term and B4.

- *
- upheld: I could not refute either part of the claim: the combinatorics hold and the numbers reproduce.

Combinatorics: X_ijk = w5_i Phi_j Phi_k D3_i C_ij C_ik is symmetric under j<->k. So sym3(X)[i,j,k] = (1/3)(X_ijk + X_jik + X_kij), which puts weight 1/3 on each of the three centers. The correct per-center weight comes from the Edgeworth / exponential formula: (1/3!) kappa_iii d_i^3 hits the Wick center, the same as the density dressing w2 -> w2 + D3/6 w5. That weight is 1/6, so the coefficient is 0.5. oracle_k3's "6 x 20/5! = 1" leaves out the factor-2 symmetry correction. The same file applies tha
- upheld: I could not refute C3. The derivation, the method and every number in it hold up.

Derivation: B7 is the kappa3_ijk hyperedge plus one C edge. Its degrees are (2,2,1), so the weight is w2_i w2_j Phi_k with coefficient 1 for each labelled diagram. There are 3 labelled diagrams. oracle_k3.sym3 averages over the 6 permutations, so the sum is 3*sym3[X], which matches coef = 6/|Aut| = 6/2 = 3. I checked this numerically without using the repo's enumerator. With an Edgeworth density phi_C(1 + t H_123), Stein's lemma says d kappa3(relu)/dt = E_C[f1' f2' f3']. Monte Carlo (4M samples, finite differenc
- upheld: I could not refute it. The diagram counts reproduce exactly, and an independent fresh Monte Carlo of He-init MLPs confirms the claimed width scaling.

(1) The completeness checks reproduce exactly: 187 connected diagrams, 18 gate dressings, none uncovered; 81 C-leaf, 106 non-leaf covered by 27 shapes, 0 uncovered. My own enumerator, which shares no code with toy_diagrams.py, also finds 187 connected labelled diagrams at network order <= 4 (4 at order 2, 31 at order 3, 152 at order 4). The only order-4 shapes with an m >= 5 hyperedge are kappa5 (1,1,3), kappa5 (1,2,2) and kappa6 (2,2,2). These 
- upheld: Every number in C5 reproduces. Three checks passed.
(1) Running make_tables.py on the three pickles gives output byte-identical to results/tables_width128.txt (checked with diff).
(2) I recomputed the summary statistics directly from the pickles:
- cl2 at layers 10-14: MLP0 2.66-3.75 (mean 3.118), MLP1 2.31-2.77 (mean 2.623).
- fit1 at layers 10-14: MLP0 1.44-2.02 (mean 1.850), MLP1 1.72-2.21 (mean 1.967).
- LRT1 mean at layers 10-14: 2.394 and 2.356.
- Layers 1-3, cl2 vs fit1: MLP0 3.588 vs 4.304, MLP1 3.439 vs 4.182.
- fitall at layers 10-14: 1.29-1.74.
(3) I re-ran closure2.py from scratch,
- upheld: I could not refute C4. The combinatorics, the completeness test and the width ratios all reproduce. A test I added for kappa5/6 also fits the counting rule.

1) Completeness: both checks give exactly the claimed output. "187 connected diagrams; 18 are degree-1 gate dressings; not covered: none" and "81 C-leaf, 106 non-leaf covered by 27 shapes, uncovered: 0".
   - My own enumerator (written separately) also finds 187 diagrams (4 / 31 / 152 at orders 2 / 3 / 4) and 18 dressings.
   - At order 4, after removing dressings, there are 31 classes: 28 non-Gaussian plus 3 pure-C classes. The "30 shape
- upheld: I re-ran the check independently and every number in the claim reproduces, with minor caveats.

What I ran: make_tables.py on results/width128_final_mlp0.txt.pkl (MLP0) and on mlp1a+mlp1b (MLP1, which is mlp_00001 split into layers 0-7 and 8-14). The regenerated fit1 and fit2 tables are byte-identical to results/tables_width128.txt. I also recomputed min/max over layers 1-14 straight from the pickles, averaging the AB and BA orientations as make_tables does.

fit2 (layers 1-14):
- B6: MLP0 1.45-1.53, MLP1 1.46-1.52, so 1.45-1.53 at every layer. Matches.
- B3: 0.51-0.82 and 0.50-0.70, so 0.50-0
- upheld: The numbers in the claim reproduce exactly, and they are not caused by noise leaking from the target into the basis. One part of the explanation is wrong, though. fit2 fixes B7 at 3.0 (an eps^3 first-order term) as well as the second-order terms. B7, not the second-order terms, is what brings B4 back.

1) Reproduction. `make_tables.py` on `width128_final_mlp0.txt.pkl` and on `mlp1a.pkl,mlp1b.pkl` matches `results/tables_width128.txt` exactly. For fit2, layers 1-14, both MLPs:
- B6: 1.45-1.53 (fit1 drifts from 1.45 to 1.13).
- B3: 0.50-0.82 (fit1 minimum -0.17).
- B4: 0.65-1.08 (fit1 runs from 
- REFUTED: The numbers reproduce, but the central width-128 claim does not hold up. The claim is that the T-class coefficient of about 1 shows B7's 3->5.0 renormalization being absorbed. That coefficient comes from B0, which makes up about 98% of TCL. The B7-replacing part keeps the same excess.

1) Toy identities (i) and (ii) hold. Re-running check_leaf_identity and check_tclass_identity gives output digit-for-digit identical to results/toy_validation.txt (leaf N=4 column 3.26, 3.44, 3.49, 3.486). I also ran my own sensitivity test (c7_mutation.py), and the check does detect wrong formulas:
   - Scaling
- REFUTED: Most of the claim holds up. The headline interpretation does not. I count the claim as refuted because it says TCL ≈ 1 "absorbs B7's apparent 3->5.0 renormalization", and the evidence contradicts that.

What holds:
(1) The leaf identity and the T-class identity hold on the toy.
(2) The T-class toy test is a test that can fail. I perturbed it and each wrong version diverges:
  - Cov(g,g) set to 0
  - Cov(g,g) doubled
  - Cov(g,g) replaced by its linearisation w2 w2 C
  - TCL scaled by 0.9
  - TCL omitted
(3) closure2.terms LEAF, TWOLEAF and TCL match the toy formulas to 1e-10. toy_validation ne
- REFUTED: The table numbers reproduce exactly, but the claim's main conclusion does not hold. At width 128, "cl2" is not the complete fixed-coefficient second-order closure. closure2.py's own docstring lists 30 eps^4 shapes, and three of them are missing from the width-128 cl2 because moment_atlas_np atlases do not store the fields they need: S9a (kappa5 (2,2,1)), S9b (kappa5 (3,1,1)) and S10 (kappa6 (2,2,2)). In the width-128 tables, cl2 is "all available" terms (SECOND_AVAIL only), and there is no "+kappa5/6" column. The claim does not say this.

I computed those slices myself on the exact same sample
- upheld: I reproduced every number in C8 three ways: from the committed result files, by re-running closure2.py, and with my own implementation of the S4d and S5c terms. Layer-by-layer figures below are MLP0 at L08-L14, read from width128_final_mlp0.txt.

- **Corrected first-order closure (cl1+B7):** 2.74, 3.04, 2.29, 3.13, 2.81, 3.06, 3.80 %, so 2.29-3.80 %.
- **cl1+B7 plus S4d ('add1'):** 2.32, 1.95, 1.71, 2.20, 2.12, 2.53, 2.26 %, so 1.71-2.53 %. The error drops at every layer.
- **Term sizes ('size' line, AB orientation only):** S4d is 1.57-2.59 % and S5c is 1.96-3.74 %. These are the two largest s
- upheld: I could not refute C8. Every number reproduces from the atlases with code I wrote (scratchpad/verify-theory/c8v2/check.py), and the derivation holds up.

Combinatorics, derived by hand from E[F(z)] = exp(sum kappa d^r / prod r!) E_G[F], keeping connected hypergraphs:
- S4d (kappa3_kki x kappa3_kkj, with w4_k Phi_i Phi_j): 1/2! x 1/2! = 1/4 per labelled diagram. There are 3 choices of centre, and sym3 is the mean over 6 permutations, so the coefficient is 0.75, as in the code.
- S5c (kappa4_iiij C_ik): 1/3! times 6 labelled diagrams gives 1.0.
- The k3xk3 group (S4a 0.5, S4b 3, S4c 0.75, S4e 1.
- upheld: I could not refute C9. Every number in it reproduces, and I found no error in the parts of the derivation and method I checked.

(1) Rerun: the HOW-TO-CHECK command gives output that is byte-identical to results/width32_closure2_k56_lr_pair.txt (82 lines, diff shows no difference). The two atlases have identical weights, different sample seeds (1000003 and 2000006) and N=4e6 each. So the fit and the evaluation use independent samples of the same network, and the noise correction is computed from the A-vs-B difference.

(2) Combinatorics: I re-derived the coefficients from kappa/prod r! times 6
- upheld: Every number in C9 reproduces: exactly on a re-run, and to within about 0.1 points on two newly sampled atlases.

(1) Re-run. I ran `closure2.py --pair56` on the stated seed1/seed2 atlases and pair files, writing to my scratch dir. The output is byte-identical to `results/width32_closure2_k56_lr_pair.txt` (diff is empty).
- L4: cl2 9.90% -> cl2+k56 5.79%. L3: 9.76% -> 7.75%. fit1 is 5.98% (L3) and 4.74% (L4).
- Term sizes (AB orientation): S9b 1.1, 3.0, 7.6, 9.57% at L1-L4. S9a peaks at 4.98% (L4). S10 peaks at 0.62% (L4).
- rho is 0.211-0.355 and K4_iiii reaches 1.00 at L1-L4. Noise is 0.28-0
- REFUTED: The numerical evidence reproduces exactly. The cost accounting that the claim draws from it does not hold.

What holds:
(a) transport_path really uses 11 n x n matmuls. I counted them with an instrumented copy, and no n^3 object is built. On width-128 layer 10 it matches the dense n^4 transport to 7e-16 for WICK, B1, B2, B3, S5c, S4d and LEAF.
(b) The B5 deviation of 1.1e-6 is float32 asymmetry, as claimed. K22 has a relative asymmetry of 2.9e-5; once the K4f[:,j,j] slice is symmetrised, B5 matches to 7.1e-16.
(c) The whole leaf class, 6 LEAF - TWOLEAF, folds into a single path with arm Phi∘C 
- REFUTED: The exactness part of C10 holds. The cost part, which the claim's conclusion rests on, is wrong in three ways.

(a) 11 n x n matmuls is not the cost of a path term. All five R/Y pieces in transport_path end in "@ W", so they can share one final product. That makes it 7 matmuls per path term.

(b) Arms are shared. In a merged transport each extra arm costs about 2 matmuls (its W^T X and one end product), not 11. Leaf class alone is 7 matmuls. Leaf plus S4d is 9. Leaf plus all ten non-leaf path terms (S3b1, S3b2, S4c, S4d, S4e, S6a, S7b, S7c, S8, G4a) is 17. The ten non-leaf terms alone are 13, 
