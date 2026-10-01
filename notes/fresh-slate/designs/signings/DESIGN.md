# Design stream "signings": matchings, hafnians and signings — DESIGN v0

Status: v0 (1 Oct 2026). Principle → derivation → estimator → cost → error model → falsification.
Prototype, Stage Q tables and verdict follow in RESULTS.md as they are measured.

## 0. The principle, stated precisely

**P1 (Wick / Isserlis).** For a centred Gaussian vector g with covariance R,
E[g_{i_1} ⋯ g_{i_{2m}}] = haf(R[i_1..i_{2m}]) — a sum over the perfect matchings of the index multiset.

**P2 (Hermite–Mehler form of P1).** For functions F_v of unit-variance coordinates g_{i_v},
with Hermite profiles F̂_v(d) = E[F_v(g) He_d(g)],

  E[ Π_v F_v(g_{i_v}) ] = Σ_{multigraphs Γ on the vertices, no self-loops}  Π_v F̂_v(deg_Γ v) / ... · Π_{e=(v,w)} ρ_{i_v i_w}^{m_e} / m_e! ,

i.e. the expectation is a *loopless hafnian* of the half-edge matrix (each vertex v carries deg v half-edges;
half-edges are matched across distinct vertices; a matched pair (v,w) is weighted by ρ_{i_v i_w}).
Joint cumulants are the same sum restricted to connected Γ. The two-vertex sector is summed in closed form
by Mehler's kernel: Cov(F(g_a), G(g_b)) = Σ_{d≥1} F̂(d) Ĝ(d) ρ_ab^d / d!.

**P3 (determinants and signings).** Two sectors of P2 are free (determinantal):
(i) the degree-≤2 sector (Γ a disjoint union of paths and cycles) is the Gaussian integral of an exponential
of a quadratic form, summed exactly by det^{-1/2}(I − DR) and a resolvent (free bosons; with signs, det^{+1},
free fermions); (ii) the degree-≤1 sector (matchings) is the matching polynomial, which equals the average of
det(x − A_s) over uniformly random signings s (Godsil–Gutman), the average being exact and the signing a
random Z/2 gauge field (a random 2-lift). Loops (cycles of Γ) are exactly what the signing average or the
determinant has to handle; trees (Bethe sector) factorise locally.

The design question set by the brief: *which of these sums does the problem actually need, in its own
terms, and can a determinant or a signed average evaluate them?*

## 1. The sums the problem needs (derivation)

Layer step: z (pre-activations at layer l), a = ReLU(z), z' = a W with W = W_{l+1} i.i.d. N(0, 2/n),
**independent of the law of z**. The score needs the per-neuron law of z'_L to about fourth order (brief §3).

  κ_k(z'_j) = Σ_{i_1..i_k} W_{i_1 j} ⋯ W_{i_k j} κ(a_{i_1}, …, a_{i_k}).                         (1)

**Step A — matchings of the weight legs.** The k legs of (1) carry fresh, independent Gaussian weights.
Apply P1 to W itself: classify each term of (1) by the set partition π of the k legs induced by equal indices
(the "weight matching"). A block of size ≥2 is *paired* (its weight product has a positive mean: coherent);
a singleton block is *unpaired* (a random sign). A quenched W contributes the realised value, but its order
of magnitude is the Wick count: each unpaired leg costs a random-sign sum (√ instead of linear). With
ρ_{ab} = O(n^{-1/2}) off the diagonal, this gives the order of every sector of (1):

| k | sector (π) | order of its contribution to κ_k(z'_j) | needed? (target ≈ 2 % of κ3, 10 % of κ4) |
|---|---|---|---|
| 2 | {2} diagonal Var a_i | 1 | yes |
| 2 | {1,1} off-diagonal Cov(a_a,a_b) | n^{-1/2} | yes, entrywise to ≈ 3·10^{-3} relative |
| 3 | {3} | n^{-1/2} | yes |
| 3 | {2,1} κ(a_a,a_a,a_b) | n^{-1} | yes (≈ 6·10^{-4} target at n = 1024) |
| 3 | {1,1,1} paths a–b–c | n^{-1} | yes |
| 3 | {1,1,1} triangles (a loop) | n^{-3/2} | no |
| 4 | {4}, {2,2} | n^{-1} | yes |
| 4 | {2,1,1} | n^{-1} (path through the paired vertex) | yes |
| 4 | {3,1}, {1,1,1,1} loops | n^{-3/2} | no |

**Consequence (the first finding of this stream, by derivation).** Because the fresh weights are independent
of the state, *every loop in a diagram whose vertices all carry fresh weight legs is suppressed by an extra
n^{-1/2}*: in (1) only tree diagrams (paths/stars, the Bethe sector) survive at the needed precision.
The determinantal resummation P3(i) and the signing average P3(ii) are therefore *not* needed for the
fresh-weight contraction — the combinatorial sum of the readout is already tree-like and costs matrix
products. Loops matter only *inside the state*: in how the joint law of z itself was built across depth.
Falsification test T1 (§6) checks the n^{-3/2} claim directly.

**Step B — the state that makes P2 applicable.** P2 evaluates every κ(a_{i_1},…) in (1) exactly if z is a
*Gaussian copula*: z_i = T_i(g_i), g ~ N(0, R) with unit diagonal, T_i increasing (the per-neuron marginal
transport). Then a_i = ReLU(T_i(g_i)) = F_i(g_i) is a function of one latent coordinate and all joint
cumulants of a are connected-multigraph (hafnian) sums with vertex profiles F̂_i(d), and with
(F_i²)^(d) for the repeated-index sectors. The state is therefore

  S_l = ( marginal law of each z_i  [m_i, s_i, γ_i, κ_i → T_i],  latent correlation R (n×n) ).

This is not covariance propagation: the marginals are non-Gaussian and the latent R is the matching
(Mehler) inverse of the pre-activation covariance, not the covariance itself; all third/fourth joint cumulants
of a are *generated* by P2 from (T, R) rather than carried.

## 2. The estimator (per-layer operation and readout)

Notation: h_d(i) = E[F_i(g) He_d(g)], q_d(i) = E[F_i(g)² He_d(g)], c_d(i) = E[F_i(g)³ He_d(g)] (d ≤ K, K = 6),
computed per neuron by Gauss–Hermite quadrature (Q = 40–64 nodes). R0 = R with zero diagonal.

1. **Mean**: m'_j = Σ_i W_ij h_0(i).
2. **Covariance (sector {2} + {1,1})**: Cov(a)_ab = Σ_{d=1}^{K} h_d(a)h_d(b) ρ_ab^d/d! (a≠b; Mehler), Var(a_i) exact;
   C' = Wᵀ Cov(a) W.
3. **Third cumulant of z'_j**:
   κ3' = Σ_i W_ij³ κ3(a_i) + 3 Σ_{a≠b} W_aj² W_bj κ21(a,b) + 3 Σ_b W_bj Σ_{p,q≥1} h_{p+q}(b)/(p!q!) X^{(p)}_{bj} X^{(q)}_{bj} (− a=c terms),
   κ21(a,b) = Σ_d [q_d(a) − 2h_0(a)h_d(a)] h_d(b) ρ_ab^d/d!,   X^{(p)} = (R0^{∘p}/p!)… (W∘h_p) (one matmul per p).
4. **Fourth cumulant of z'_j**: Σ_i W⁴ κ4(a_i) + 3 Σ_{a≠b} W_a² W_b² κ22(a,b) + 6 Σ_a W_a² q̃_2(a) (X^{(1)}_a)² + …
   (κ22 from the Mehler sums of the profiles of a, a², centred; the {2,1,1} path through the paired vertex
   reuses X^{(1)}).
5. **New state**: marginal transport T'_j from (m', C'_jj, κ3', κ4') by Cornish–Fisher; latent R'_jk by
   solving Σ_{d=1}^{3} t_d(j) t_d(k) r^d/d! = C'_jk elementwise (Newton, 3 steps), t_d the Hermite profile of T'.
6. **Readout** at every layer: E[a'_j] = E[ReLU(T'_j(g))] by quadrature.

## 3. Cost at n = 1024 (units: 1 unit = one 1024³ matmul = 2^31 FLOPs; float32)

| item per layer | units |
|---|---|
| C' = Wᵀ Cov(a) W (two products; second is a symmetric-output product) | 2.0 (1.5 with the aliased form) |
| X^{(1)}, X^{(2)}, X^{(3)} (matmuls) | 3.0 |
| κ21 contraction (K21 W, then column sums with W∘W) | 1.0 |
| κ22 contraction ((W∘W)ᵀ K22 (W∘W), output diagonal only: n² per column → matmul) | 1.0 |
| Mehler/elementwise (K = 6 powers of R0, Newton inversion, ≈ 60 flops/entry) | 0.03 |
| quadrature (n·Q·K) | < 0.001 |
| **total** | **≈ 7–8 per layer → ≈ 115 units ≈ 0.11 B** |

Trims: the X^{(3)} and κ22 products are candidates to drop if T1/T2 show them small (→ ≈ 5/layer, 80 units, at
the 0.1 floor). Python-side calls: ≈ 40 per layer, ~700 total, far inside the residual cap.

## 4. Error mechanism and predicted scaling

The estimator is exact (up to truncations K, p,q ≤ 3) for one layer step out of a Gaussian copula; layer 1 is an
exact copula (z_1 = xW_1 Gaussian). The error is the *copula defect*: the true law of z_l is not a Gaussian
copula; its joint cumulants beyond those generated by P2 from (T, R) — the content born at earlier ReLUs and
transported linearly (brief §3, "old content") — are dropped every layer. Prediction: the defect in the {2,1}
and {1,1} sectors is O(n^{-1}) per layer with an O(1) prefactor that grows with depth; κ3' and C' off-diagonals
then carry relative errors of fixed order, so the final-mean error is ≈ n^{-1}·φ-scale per neuron, i.e.
raw MSE ∝ n^{-2}·c(L) at fixed depth. If the copula captured ≥ 60 % of the third-order content (the brief's
"40 % old content" says it cannot capture more than ≈ 60 % of the {2,1} slice), raw MSE at 1024 would be
≈ 1e-8 – 1e-7. That is the decisive number.

## 5. What the principle buys and what is lost

- Buys: an exact, closed-form generator of all joint cumulants from (T, R) (Mehler resums all multi-edges);
  a principled selection of the diagrams needed (weight matchings); cost dominated by ≈ 7 matmuls per layer.
- Lost: (a) all non-copula joint structure (old content), (b) loop diagrams (claimed negligible, T1),
  (c) Cornish–Fisher marginal truncation (fourth order).
- Determinants/signings are *not* used in the readout, by the derivation in §1 Step A. Their natural place is
  the state (loops through depth); §7 records the candidate.

## 6. Cheapest falsification tests

- **T0 (exactness, toy):** width 12, depth 2: the copula is exact at layer 1, so κ2..κ4 and E[a_2] from §2
  must match high-precision Monte Carlo to its noise. Failure = derivation/implementation bug.
- **T1 (loops):** at widths 32–256, layer 2, compute the triangle sector exactly (n³ per neuron) and compare with
  the tree sectors: its share must fall as n^{-1/2}.
- **T2 (Stage Q):** widths 64/128/256, depth 16, 4 MLPs each: raw final-layer MSE (truth noise subtracted)
  and its width slope. The design is dead if raw MSE at 256 is above ≈ 1e-6 with a slope shallower than n^{-1}
  (projects above 1e-7 at 1024, i.e. adjusted > 1e-8).

## 7. Where signings would enter (next iteration, if T2 shows the copula defect dominates)

The dropped content is a sum over *paths through depth* that end on loops closed by the next layer's weights.
A 2-lift (random signing of the layer weights) of the network kills exactly the diagrams whose cycles have odd
signed holonomy; averaging the copula estimator over signed lifts versus the quenched network isolates the
loop (old-content) part as a difference of two cheap evaluations — a candidate correction to test once the
size of the copula defect is known.
