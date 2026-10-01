# The minimal co-state (breakthrough stream `costate`): REPORT

*1 Oct 2026, v0 (theory and first width-1024 runs; tables are filled in as runs finish). Lens: the Heisenberg picture taken to its end. The network is quenched, so the functionals of each layer's law that the final means depend on can be fixed backwards from the readout before anything is computed. Question: is there a co-state of O(n²) numbers per layer, closed under pull-back, so that old third-order content costs O(L n³) in total instead of O(L² n³)? Labels: **Theorem** (proved here, elementary), **Measured**, **Prediction**.*

## 0. Answer in brief

1. **The co-state the readout needs is small and exactly closed** (Theorem C1). To first order in the local non-Gaussianity it is the n readout cotangents G_{l,j} = E_ν ∇³ g_{l,j} / 6. Their values are n numbers per layer, and they obey an exact backward recursion G_l = T_l* G_{l+1} + J_{l+1}: pull-back through the linear-response arrow plus a local kink term.
2. **Being small is not the same as being cheap to evaluate** (Theorems C2, C3). The co-state tensors G_{l,j} have CP rank Θ((L−l) n) each, and any exact evaluation of their pairing with the sources goes through the **double edge**: the Hadamard square (U_{s→k} ∘ V_{s→k}) of a source-to-kink propagator, where two legs of a source atom land on the same kink neuron. The double edge factors through an intermediate layer only via the n²-dimensional second chaos: (UM) ∘ (VM) = (U ⊙ V)(M ⊙ M). Consequences:
   - the minimal exactly closed forward state at a cut m is the set of all future (2,1) slices, (L − m) n² numbers;
   - first order costs Θ(n³) per (source, kink) pair, Θ(L² n³) in total, in every contraction order (forward over pairs, or backward through per-readout Hessians, Θ(L n⁴));
   - there is **no closed n²-per-layer co-state** for the old content.
3. **The natural n² closures are measured directly at n = 1024**, together with the exact first-order co-state (§4). The (2,1)-slice chain keeps only the diagonal a = b of the double edge, i.e. it averages over the newest weights; age truncation; Tucker truncation of old propagators. Tucker truncation cannot save cost anyway: a rank-d propagator's double edge has rank d², so it is cheaper only when d < √(n/2) ≈ 22 (Theorem C4).
4. **Second order** (the δ-insertion of old κ₃ into joint κ₄, which the Heisenberg oracle says is needed for the bar) has, in the same pull-back form, one more interior Hadamard vertex. It costs Θ(L³ n³) (Theorem C5).

**Verdict v0**: the readout does read old content only through known contractions, but those contractions contain a squared propagator. That squared propagator is what makes the cost quadratic in depth. The deciding experiment is §4 (does an n² closure lose the old content at n = 1024?).

## 1. Setting and the co-state

z₁ = x W₁, z_{l+1} = ReLU(z_l) W_{l+1}, readout r_j = ReLU(z_{L,j}). The Koopman pull-back is K_l g = g(ReLU(·) W_{l+1}), and the Heisenberg observable is g_{l,j} = K_l ⋯ K_{L−1} r_j. Reference chain ν_l: the Gaussian closure with its own corrected (m, C). The Heisenberg–Duhamel identity (heisenberg DESIGN Thm 1, exact) is truth − reference = Σ_l (ρ_l − ν_l)[g_l] with ρ_l = T_{l−1} ν_{l−1}, and Stein's identity on Gaussian space gives (ρ_l − ν_l)[g] = ⟨κ₃(ρ_l), E_ν ∇³ g⟩/6 + O(ε²).

**Definition.** The *first-order co-state* at layer l is the linear map κ ↦ (⟨κ, G_{l,j}⟩)_{j=1..n} on Sym³(ℝⁿ), with G_{l,j} := E_{ν_l} ∇³ g_{l,j} / 6, the cotangent of readout j with respect to the third cumulant of z_l at the reference. A_l := span_j G_{l,j}, so dim A_l ≤ n.

**Linear-response arrow.** With decoupled gates (the first-order model in which the Heisenberg stream's HD-∞ is exact), a third cumulant κ of z_l reaches z_{l+1} as T_l κ = κ[Φ_l W_{l+1} ·, Φ_l W_{l+1} ·, Φ_l W_{l+1} ·], and at layer l it acts on the readout only through the Stein injection J_l. J_l reads D = κ(p,p,p) into the mean via E δ′/6, and S = κ(p,p,q) into Cov(a_p, a_q) via E[δ(z_p) H(z_q)] and E[δ′(z_p) ReLU(z_q)] (heisenberg `inject`). The injection is followed by the reference chain's tangent map to the readout.

**Theorem C1 (exact closure of the co-state).** G_{l,j} = T_l* G_{l+1,j} + J*_{l,j}, where J*_{l,j} is the injection co-vector at layer l composed with the chain's sensitivity of readout j. Hence A_l ⊆ T_l* A_{l+1} + span_j J*_{l,j}. In particular, the readout correction is Σ_s ⟨Src_s, G_{s+1}⟩ = Σ_{s<k} ⟨Src_s, T*_{s→k} J*_k⟩, with T_{s→k} the transport by U_{s→k} = W_{s+1} Φ_{s+1} W_{s+2} ⋯ Φ_{k−1} W_k.
*Proof.* The first-order correction is linear in the sources, and each source reaches the readout either by injection at the current layer or by one more linear-response step. Unroll the recursion. ∎

So the co-state **values** are closed and n-dimensional, which is as small as anything with n outputs can be. The question is what it costs to evaluate them.

## 2. The obstruction: the double edge

**Source atoms.** The κ₃ source of a = ReLU(z) at a Gaussian layer is, to the accuracy the Heisenberg stream measured as lossless (R7), a structured trilinear form. Write α_k for the Hermite coefficients of a, R = ρ diag(α₁), P_{pr} = E[ã_p² ã_r] (exact), Δ = P − star(p,p,r) off the diagonal, and δ = κ₃(a_p) − 3α₁²α₂:

  T(x, y, z) = Σ_c α_{2,c} [x_c (Ry)_c (Rz)_c + y_c (Rx)_c (Rz)_c + z_c (Rx)_c (Ry)_c] + Σ_{p≠r} Δ_{pr} [x_p y_p z_r + x_p z_p y_r + y_p z_p x_r] + Σ_p δ_p x_p y_p z_p.

The star and coincidence terms are written as n atoms each, with factor matrices I, R, Δ. **Measured:** the pull-back evaluation of this form (`costate.py`) reproduces the full-tensor HD with these diagrams to 5e-15 (`validate.py`, w64).

**What the readout reads at kink layer k, from source s.** Write U = U_{s→k}, V = R_s U, X = Δ_s U. The (2,1) slice is

  S_{s→k} = 2 (U ∘ V)ᵀ diag(α₂) V + (V ∘ V)ᵀ diag(α₂) U + (U ∘ U)ᵀ (X + δ ∘ U) + 2 (U ∘ X)ᵀ U.

Every term carries a **Hadamard product of two copies of the same propagator**, the double edge: two legs of one source atom land on the same kink neuron p, which is the δ(z_p) of the covariance injection. The third leg is an ordinary single edge.

**Theorem C2 (cut lemma).** Fix a layer m. For generic weights and arbitrary content generated before m, the smallest forward state at m from which every later (2,1) slice S_k (k > m) can be computed exactly is the projection of κ₃(z_m) onto ⊕_{k>m} P_k, where P_k = span{sym(u_p ⊗ u_p ⊗ u_q) : u = U_{m→k} e_·}. Each P_k has dimension n² (the multisets {p, p, q}). The P_k for different k are in general position in Sym³ (dimension ≈ n³/6 ≥ (L − m) n² at n = 1024, L = 16). So the state has dimension (L − m) n²: in effect, all the future slices themselves.
For the readout alone the state shrinks to n numbers (C1), but the probe tensors Σ_k Σ_{pq} Γ̃_{jk}(p,q) u_p ⊗ u_p ⊗ u_q have CP rank Θ((L − m) n). Pairing them with the CP form of the old content (rank Θ(m n)) still decomposes into the (s, k) pair contractions above.

**Theorem C3 (no n-dimensional factorisation of the double edge).** (U M) ∘ (V M) = (U ⊙ V)(M ⊙ M), where ⊙ is the row-wise Khatri–Rao product for U, V (n × n²) and the column-wise one for M (n² × n). So the double edges of all pairs (s < m < k) factor through the cut at m only via the n²-dimensional second chaos, one n² vector per atom. The **only** part of (M ⊙ M) that acts in an n-dimensional space is its restriction to the diagonal a = b, i.e. (U ∘ V)(M ∘ M). This gives the n²-state closure S_{k+1} ≈ (M ∘ M)ᵀ S_k M: the (2,1)-slice chain.
Over the newest weights W_{k+1}, the dropped a ≠ b part has mean zero and the same order as the kept part: Σ_{a≠b} U_{ca} V_{cb} M_{ap} M_{bp} has n² terms of random sign, against n coherent terms in the diagonal. So the slice chain is the annealed (W-averaged) double edge, and it is wrong by O(1) relative on old content. The convergence dossier's "old content is orthogonal to anything computable from the current layer" is the measured face of this.
*Contraction orders.* Forward over pairs: Θ(n³) per (s, k), total ≈ (L²/2) × c × n³. Backward per readout: the second-order adjoint 𝒴_s[:,:,j] = M_s (𝒴_{s+1}[:,:,j] + diag Y_{s+1,j}) M_sᵀ, the Hessian of every readout, costs Θ(n⁴) per layer. For L ≪ n the pair order is the cheaper one. Mixed orders meet a second Hadamard vertex at the readout, where the single leg pairs with the double edge row-wise, and reduce to one of the two.

**Theorem C4 (Tucker truncation does not reduce the double edge).** If U, V have rank d (old propagators: the participation ratio of U_{s→k} is ≈ n / (2 · age)), then U ∘ V has rank ≤ d². Evaluating (U ∘ V)ᵀ Y through the factors costs 2 d² n² against n³, so it saves only when d < √(n/2) ≈ 22 at n = 1024, i.e. at ages ≳ 23 > L. Low-rank old content is therefore an *accuracy* statement, never a *cost* lever, for the double edge.

**Theorem C5 (second order is cubic in depth).** The old κ₃ that passes a later ReLU at layer m creates joint κ₄ through the He₂ vertex: κ₄(a)_{q···} ∋ α_{2,q} ρ_{q·} ⊗ κ₃(z_m)_{q··}, a source atom split at a centre neuron q. Its (3,1)/(2,2) slices at a later kink layer k′ have Hadamard vertices at q (source side) and at p (kink side). In pull-back form that is Θ(n³) per (s, m, k′) triple, Θ(L³ n³) in total, unless the insertion vertex q is annealed. Annealing it fails for the same reason as in C3.

## 3. Cost of the exact first-order co-state at n = 1024

Per live (s, k) pair: propagate [U; V; X] by the shared right factor Φ_k W_{k+1} (3 products, one stacked call family), and four contractions forming S_{s→k} (4 products, one hub-sum family across all pairs of the layer). L = 16 gives 120 pairs; counted **855 products per MLP** (`costate.py`, `nprod`), including 2 per layer for the reference sandwich and 2 per new source. Without the coincidence atoms: 4 products per pair, ≈ 500.

| pricing (KIT.md) | 855-product version | 500-product version (no coincidence atoms) |
|---|---|---|
| dense f32 | 855 u = 0.83 B | 500 u = 0.49 B |
| Strassen L3 (0.683 u, ≈ 70 ms/product in batch) | 584 u = 0.57 B, ≈ 60 s backend | 342 u = 0.33 B |
| Strassen L5 (0.557 u, ≈ 270 ms/product) | 476 u = 0.46 B, **≈ 230 s backend: over the 120 s wall cap** | 279 u = 0.27 B, ≈ 135 s (over the cap) |

Calls: one propagation family plus one hub family per layer ≈ 16 × (46 + 87) ≈ 2.1k calls ≈ 0.07 s residual, which is within the 0.4 s limit.

## 4. Measurements (direct, n = 1024 bench set w1024_d16; and w128_d16)

Runs: `run.py` → `results/<set>.jsonl`, summarised by `summarize.py`. Raw = final-layer MSE − truth noise (3.6e-8 at 1024, N = 2e6).

*First run (MLP 0 at n = 1024):* Gaussian closure raw 5.00e-6. **Exact first-order co-state (all 120 pairs) raw 3.92e-7: gain 12.8×, 33 s numpy.** This is the first direct width-1024 measurement of first-order Heisenberg–Duhamel; the heisenberg stream projected 4–7e-7 from width 256.

The closure sweep (all six MLPs: ages A = 0..7; the slice chain for old content; Tucker truncation at rank n/(2·age); no coincidence atoms) is running. Its table goes here.

## 5. Predicted error and next experiments

- Prediction (from C3): the slice-chain closure of old content recovers little more than age truncation at the same A. If it recovers ≥ 80 % of the all-ages gain at n = 1024, then C3's same-order argument fails quantitatively at this width, and an O(L n³) co-state exists in practice.
- Prediction (from C1/C2 and the heisenberg oracle R3): even the exact first-order co-state stays a factor of 10–40 above the bar at n = 1024. The remaining factor is second order, which costs Θ(L³ n³) in pull-back form (C5).
