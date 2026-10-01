# The minimal co-state (breakthrough stream `costate`): REPORT

*1–2 Oct 2026, v1 (theory; direct width-1024 measurements on all six bench networks; w128 oracles). Lens: the Heisenberg picture taken to its end. The network is quenched, so the functionals of each layer's law that the final means depend on can be fixed backwards from the readout before anything is computed. Question: is there a co-state of O(n²) numbers per layer, closed under pull-back, so that old third-order content costs O(L n³) in total instead of O(L² n³)? Labels: **Theorem** (proved here, elementary), **Measured**, **Prediction**.*

## 0. Answer in brief

1. **The co-state the readout needs is small and exactly closed, but not cheap to evaluate** (C1–C3). To first order it is the n readout cotangents, whose values satisfy an exact backward recursion. Evaluating them, though, needs the double edge (U_{s→k} ∘ V_{s→k}): two legs of a source atom meeting at one kink neuron. That factors through an intermediate layer only via the n²-dimensional second chaos, so exact first order costs Θ(n³) per (source, kink) pair, Θ(L² n³) in total, and second order Θ(L³ n³) (C5). No closed n²-per-layer co-state exists for old content in general. Every causal n² closure measured (diagonal slice chain, (2,1)-support renewal, Tucker truncation) loses a large part of the old content at n = 1024 (§4.1).
2. **But the readout reads old content through a much smaller object than the forward state** (§4.2, oracle). With young content (age ≤ 1) exact, the old (2,1) slice is needed only through its diagonal plus a rank-1 spike: 4.6e-5 against 3.1e-5 for all of it, at w128, where rank 8 is enough. At 1024 the rank-1 spike gives ≈ 80 % and the rest is high-rank at Stein level (§4.2b), but a law-level treatment of the spike beats every Stein-level filter. The spike lies along (s², m) and is the **dilation (global scale) mode**. Positive homogeneity makes that mode an **exactly conserved charge** of the true dynamics: Theorem C6, F#(P ⋆ ν) = P ⋆ F#ν. It is also the only exactly closed symmetry sector for generic weights (C7). This is the same object as the bethe stream's "global order parameter", seen here from the readout side in κ₃.
3. **Carrying it costs one scalar per layer.** Content of age > A is projected onto the scale-mixture slice when it retires (an O(n²) dot product). The accumulated scale variance v is used at law level (ReLU moments from the de-scaled covariance, scale re-added: exact by homogeneity). **Measured directly on all six w1024_d16 networks:**

| variant | raw (6 nets, paired) | vs Gaussian | products | units (dense / Strassen L3) |
|---|---|---|---|---|
| Gaussian closure | 4.10e-6 | 1 | 30 | 0.03 B |
| exact first-order co-state (all 120 pairs) | 4.23e-7 ± 0.21e-7 | 9.7× | 855 | 0.83 / 0.57 B |
| A = 1, old dropped | 1.90e-6 | 2.2× | 218 | |
| **A = 1 + scale law (A1gl)** | **6.34e-7 ± 0.71e-7** | 6.5× | **218** | 0.21 / 0.15 B |
| A = 2 + scale law (A2gl) | 4.96e-7 ± 0.38e-7 | 8.3× | 309 | 0.30 / 0.21 B |
| A = 3 + scale law (A3gl) | 4.24e-7 ± 0.29e-7 | 9.7× | 393 | 0.38 / 0.26 B |
| A = 3 + scale law + slice-chain residual (A3gsl) | 3.91e-7 ± 0.22e-7 | 10.5× | 435 | 0.42 / 0.29 B |
| A = 2 + scale law, no coincidence atoms (A2gl_nc) | 4.33e-7 ± 0.25e-7 | 9.5× | 183 | 0.18 / 0.12 B |
| A = 3 + scale law, no coincidence atoms (A3gl_nc) | 3.48e-7 ± 0.20e-7 | 11.8× | 231 | 0.23 / 0.15 B |
| **A3gsl, no coincidence atoms (A3gsl_nc)** | **3.25e-7 ± 0.17e-7** | **12.6×** | **273** | **0.27 / 0.18 B** |

   The scale-mode co-state matches the exact all-pairs first order at **46 % of its products** (A3gl), and A3gsl edges below it (ratio 0.68–1.11 per network, mean 0.93). **Without the coincidence atoms (at 1024 they are a 1/n effect, and with the scale law they hurt) it is better and cheaper: A3gsl_nc beats the exact all-pairs first order on all six networks (ratio 0.69–0.91) at 32 % of its products.** Best adjusted: A3gsl_nc ≈ 3.25e-7 × 0.18 ≈ 5.9e-8 and A2gl_nc ≈ 4.3e-7 × 0.12 ≈ 5.2e-8 with Strassen L3. That is better than the adjusted numbers the convergence dossier reports for the other fresh designs (faces ≈ 1.1e-7 with Strassen, bethe 1.07e-7). At w128 the law-level scale mode even beats all pairs by 1.2×: A3gsl 2.64e-5 against 3.13e-5, which is the mode's own second-order (κ₄-spike) content.
4. **Verdict.** The co-state lens answers the dossier's question (ii) **positively in effect and negatively in principle**. Old third-order content is not O(L n³) as an exact object. Its *readout-relevant* part, however, is one conserved scalar (the dilation charge) plus a residual that a short exact window (A ≤ 3) and the slice chain absorb, which takes it from O(L² n³) to O(A L n³). But the first-order co-state itself sits at raw ≈ 4e-7 at n = 1024, about 40× above the bar's raw ≈ 1e-8. Best adjusted ≈ 5–6e-8 at 0.12–0.18 B (A2gl_nc, A3gsl_nc; raw 3.3–4.3e-7). That is ahead of or level with the other fresh designs (faces raw 3.2e-7, markov raw 2.4e-7, bethe raw 1.04e-6), and ≈ 35× short of 1.6e-9. What is missing is **not old content any more, and not the first-order κ₃ model either**: at w128 the computed co-state equals first-order injection of the *true* κ₃ (2.70e-5 against 2.77e-5, four networks). **The whole remaining gap is joint κ₄ (the second-order Duhamel term): 15× at w128, and its readout-relevant part is concentrated, with rank 8 recovering 8× and rank 32 recovering 11× (§4.6).** The next object for this lens is therefore a few readout-relevant joint-κ₄ modes, identified in closed form as the dilation mode was for κ₃.

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

## 4. Measurements: closures of the old content (direct n = 1024 and w128)

Runs: `run.py` → `results/<set>[_suffix].jsonl`, summarised by `summarize.py`. Raw = final-layer MSE − truth noise (3.6e-8 at 1024, N = 2e6). "Products" = n³ products counted by the prototype. A = largest age carried exactly (age = k − s − 1).

**4.1 The exact first-order co-state and the natural n² closures**

| variant | w1024 raw (3 MLPs) | gain | w128 raw (8 MLPs) | gain | products |
|---|---|---|---|---|---|
| Gaussian closure | 4.82e-6 | 1 | 2.82e-4 | 1 | 60 |
| **all pairs (exact first order)** | **4.03e-7 ± 0.06e-7** | **12.0** | 3.13e-5 | 9.0 | 855 |
| all pairs, no coincidence atoms | 5.08e-7 | 9.5 | 6.09e-5 | 4.6 | 495 |
| A = 0 / 1 / 2 / 3 | 3.03e-6 / 2.09e-6 / 1.44e-6 / 1.05e-6 | 1.6 / 2.3 / 3.3 / 4.6 | 1.92e-4 / 1.24e-4 / 8.6e-5 / 6.1e-5 | 1.5 / 2.3 / 3.3 / 4.6 | 120 / 218 / 309 / 393 |
| A = 5 / 7 | 6.57e-7 / 4.84e-7 | 7.3 / 10.0 | 4.0e-5 / 3.2e-5 | 7.0 / 8.9 | 540 / 659 |
| A = 0 / 1 / 3 + diagonal slice chain for older | 1.57e-6 / 1.17e-6 / 6.87e-7 | 3.1 / 4.1 / 7.0 | 8.0e-5 / 6.0e-5 / 3.7e-5 | 3.5 / 4.7 / 7.6 | 174 / 268 / 435 |
| Tucker truncation of propagators at rank n/(2·age) | 6.97e-7 | 6.9 | 4.2e-5 | 6.7 | 855 (C4: no saving) |
| pool renewal (re-project to (2,1) support, exact after) | — | — | A0: 6.8e-5 … 8.9e-5 | ≤ 4.2 | — |

Readings. (i) **The exact first-order co-state at n = 1024 measures raw 4.0e-7, 12× below the Gaussian closure.** That is the first direct width-1024 number for first-order Heisenberg–Duhamel (projected 4–7e-7 from width 256), and it ties with the best fresh design so far (faces, 3.2e-7). (ii) Old content is essential at 1024 and decays slowly with age: A = 7 still leaves 20 %. (iii) Every causal n² closure loses a large share (C2–C3). Re-projecting onto (2,1) support even once (pool renewal) loses most of it: the all-distinct part matters.

**4.2 What the readout actually reads of the old content (oracle filters, w128).** The old content is carried exactly. Only what is injected at each layer is filtered: young (age ≤ A) exact, old (age > A) filtered.

| injected old slice | w128 raw |
|---|---|
| nothing (A = 1) | 1.24e-4 |
| diagonal only (per-neuron skewness) | 1.06e-4 |
| off-diagonal only | 9.0e-5 |
| exact diagonal + **rank-1** off-diagonal (SVD) | **4.65e-5** |
| exact diagonal + rank 8 / rank 64 | 3.30e-5 / 3.13e-5 |
| A = 0, exact diagonal + rank 8 | 3.27e-5 |
| everything (all pairs) | 3.13e-5 |

**The readout reads the old content almost entirely through a rank-1 spike plus a few modes.** In the readout metric, the old (2,1) slice needs about one mode, not the ≥ 0.3 n modes the Frobenius-metric studies needed (old-content stream). Its anatomy (`spike.py`, w128, A = 1): the off-diagonal old slice is 4–5× the diagonal in norm. Its top singular pair carries 43 % of the energy at layer 4, 69 % at 6, 82 % at 8 and 91–95 % from layer 10 on. The left singular vector aligns with the per-neuron scale s (cos 0.8–0.93), the right with the pre-activation mean m (cos 0.88–0.96).

**4.2b The same filters at n = 1024** (3 networks, paired; all pairs 4.03e-7, A1 2.09e-6 on the same networks):

| injected old slice (age > A, Stein level, oracle) | raw |
|---|---|
| A = 1, diagonal only | 1.47e-6 |
| A = 1, one scale scalar fitted to the exact old slice | 1.08e-6 |
| A = 1, exact diagonal + rank 1 | 7.25e-7 |
| A = 1, + scale-field vector h fitted (4.3), + diagonal | 7.17e-7 |
| A = 1, exact diagonal + rank 8 / rank 32 | 6.66e-7 / 5.47e-7 |
| A = 3, exact diagonal + rank 8 | 4.82e-7 |
| A = 3, one scale scalar fitted | 6.08e-7 |
| *non-oracle* A3gl_nc / A3gsl_nc (law-level scale) | 3.5e-7 / 3.25e-7 (6 nets) |

At 1024, the diagonal plus a rank-1 spike gives about 80 % of the old content's MSE reduction at A = 1. The last 20 %, at Stein level, needs ≳ 32 modes (the old-content stream's growth with n, now in the readout metric). But **every Stein-level filter, even an oracle one, loses to the non-oracle law-level scale mode plus an exact window**. The scale mode's higher cumulants are worth more than all the remaining modes of the old (2,1) slice.

**4.3 The spike is the global scale (dilation) mode.** If z = t·x with a scalar t (E t = 1, Var t = v) independent of x ~ N(m, C), then κ₃(z_p, z_p, z_q) = 2v (2 m_p C_pq + m_q C_pp) + O(v², κ₃(t)). Off the diagonal this is dominated by 2v s_p² m_q, rank one along (s², m), which is exactly the measured spike.

**Theorem C6 (the dilation sector of the co-state is exactly closed).** Every arrow is positively homogeneous, F_l(t z) = t F_l(z) for t > 0. Hence for every law ν and every mixing law P on (0, ∞), F_l#(P ⋆ ν) = P ⋆ (F_l#ν), where P ⋆ ν is the law of t·z with t ~ P independent of z ~ ν. The mixing law is transported **unchanged**, a conserved charge of the exact dynamics. For the homogeneous readout, E_{P⋆ν}[r] = E_P[t] · E_ν[r]. So the co-state has a sector that is exactly closed under pull-back and whose dimension does not grow: the scalar law P. Pathwise, ReLU′(z) z = ReLU(z) (Euler), so the scale direction sym(m ⊗ C) of the cumulant dynamics is an eigen-direction of eigenvalue 1 of the exact dynamics. The decoupled-gate linear response of first-order HD does not have this property (E a ≠ Φ m), so it leaks the mode.
In the programme's terms, the dilations form a one-parameter automorphism group commuting with every arrow, like a gauge action. The scale-mixture states form the family it generates, and its fixed-point (invariant) sector is a sufficient sub-algebra for the scale direction. Per C1–C3, all the rest of the old content costs pairs.

**4.4 Carrying the scale mode (non-oracle).** Content that ages out (age > A) is projected, at the layer where it retires, onto the scale-mixture slice K = 2 m ⊗ C + diag(C) ⊗ m: γ = ⟨S_old, K⟩ / ⟨K, K⟩, an O(n²) dot product of matrices already computed. The amplitudes are accumulated (conserved charge) and used as v = γ/2 at every later layer, in two ways:
- **gp**: first-order Stein injection of γ K;
- **gl** (law level): ReLU moments from the de-scaled covariance C_x = (C − v m mᵀ)/(1 + v), and Cov(a) = (1 + v) Cov_x(a) + v E a E aᵀ. This uses homogeneity exactly, E[ReLU(t x_p)] = E[t] E[ReLU(x_p)], so the mode's higher cumulants (its κ₄ spike, the bethe stream's "global order parameter") come in at no cost.
- **gs / gsl**: the same, plus the residual (S_old − γ K) carried by the diagonal slice chain.

Extra cost over A-truncation: O(n²) per layer for gp/gl, 2 products per layer for the slice chain.

| variant (w128, 8 MLPs) | raw | products |
|---|---|---|
| A = 1 / 2 / 3, old dropped | 1.24e-4 / 8.6e-5 / 6.1e-5 | 218 / 309 / 393 |
| oracle: A = 1 + best single scalar per layer (γ fitted to the exact old slice) | 5.25e-5 | — |
| gp: A = 0 / 1 / 2 / 3 (first-order injection) | 2.7e-4 / 7.3e-5 / 4.6e-5 / 3.7e-5 | 120 / 218 / 309 / 393 |
| gs: A = 1 / 3 | 6.4e-5 / 3.2e-5 | 268 / 435 |
| **gl: A = 1 / 2 / 3 (law level)** | **3.66e-5 / 3.10e-5 / 2.81e-5** | 218 / 309 / 393 |
| **gsl: A = 1 / 3** | **3.04e-5 / 2.64e-5** | 268 / 435 |
| all pairs (exact first order) | 3.13e-5 | 855 |

Projecting at age 0 is wrong: the source has not yet turned into the scale mode, and A0gp is worse than A0. From age 1 on, the projection works. Over a fitted scalar at first order (A1gsm oracle 5.25e-5), the law-level version gains another 1.4×, and gl/gsl beat the exact all-pairs first order at a quarter to half of its products. That is direct evidence that the scale mode carries the second-order content: its κ₄ spike, which is coherent.

**4.5 Width 1024 (the decisive set; six networks, paired, truth noise subtracted).** Table in §0. Further readings:
- A first-order scalar (gp) recovers half of the old content: A1gp 9.4e-7 against A1 1.90e-6. The law-level treatment of the same scalar recovers most of the rest (A1gl 6.3e-7). The difference is the mode's κ₄ (and higher) content, which homogeneity gives exactly.
- Per network, A3gsl / all pairs = 1.11, 0.97, 1.10, 0.68, 0.83, 0.91.
- (Per-layer MSEs in the result files include the per-layer truth noise Var(a_l)/N. That is 3.4e-7 at layer 1, where every method is exact, and it falls with depth to 3.6e-8 at layer 16. Read only the final layer.)

**C6 stated abstractly (for the framework; rule 3).** Let a group G act on every object of a layered system by maps that intertwine every arrow, F ∘ g = g ∘ F. In the Koopman picture this is an automorphism group α_g of each layer algebra with α_g K = K α_g. Then:
- (i) G-mixtures of states are transported with their mixing law unchanged: F#(P ⋆ ν) = P ⋆ F#ν.
- (ii) The conditional expectation onto G-invariant observables commutes with the transfer operator.
- (iii) A readout that is a G-eigen-observable (α_g r = χ(g) r) reads the mixing law only through E_P[χ].
So the co-state contains an exactly closed sector, the G-mixing law, whose size depends only on G and not on depth. Everything outside it is subject to the double-edge obstruction (C3). For a Bratteli diagram or a tiling with a scaling symmetry the statement reads the same. Here G = (ℝ_{>0}, ×) acts by dilation, χ(t) = t, and C7 says it is the whole symmetry group for generic weights. The operational content is a split: carry the invariant sector at law level, for free; spend pairs only on the rest.

**Theorem C7 (uniqueness of the exactly closed symmetry sector, sketch).** Suppose a family of invertible maps g_l intertwines the arrows, F_l ∘ g_l = g_{l+1} ∘ F_l, with F_l = ReLU then W_{l+1}. Commuting with the coordinatewise ReLU (and preserving the orthant structure on which it is linear) forces g_l to be a positive diagonal scaling composed with a permutation. Intertwining a dense W with i.i.d. continuous entries (W g = g′ W) then forces g = c I almost surely. So the dilations are the only continuous symmetry of a generic network, and the scale-mixture law is the only exactly closed sector of the co-state beyond the n readout values themselves. Everything else in the old content must be paid for in pairs, or approximated.

**4.6 What remains: the second-order oracle (w128, `oracle2nd.py`).** True Monte Carlo cumulant slices (two passes, N = 4e6; heisenberg oracle3 atlas code) are injected at every layer into the own (m, C) chain. Table below; `table2nd.py` regenerates it.

| (w128, networks 0, 1, 2, 3) | raw per network | geometric mean |
|---|---|---|
| Gaussian closure | 4.74e-04 / 1.96e-04 / 2.34e-04 / 6.11e-04 | 3.40e-04 |
| computed, exact first order (all pairs) | 4.94e-05 / 3.48e-05 / 4.06e-05 / 3.70e-05 | 4.01e-05 |
| computed, best co-state A3gsl_nc | 3.02e-05 / 2.69e-05 / 2.21e-05 / 2.97e-05 | 2.70e-05 |
| oracle: true κ₃, first-order injection | 1.01e-05 / 3.39e-05 / 2.24e-05 / 7.62e-05 | 2.77e-05 |
| oracle: true κ₃ + joint κ₄ (second order) | 2.49e-06 / 2.72e-06 / 1.63e-06 / 9.42e-07 | 1.80e-06 |
| oracle: κ₄ slices rank 1 off-diagonal + exact diagonal | 5.08e-06 / 9.68e-06 / 1.21e-05 / 3.89e-06 | 6.94e-06 |
| oracle: κ₄ slices rank 8 + diagonal | 3.10e-06 / 4.71e-06 / 3.04e-06 / 2.71e-06 | 3.31e-06 |
| oracle: κ₄ slices rank 32 + diagonal | 3.55e-06 / 3.16e-06 / 2.14e-06 / 1.73e-06 | 2.54e-06 |
| oracle: κ₄ diagonal only | 1.65e-05 / 1.44e-05 / 1.58e-05 / 1.23e-05 | 1.47e-05 |

Reading (4 networks, geometric means):
- (a) **The computed co-state is as good as the true κ₃ at first order.** A3gsl_nc 2.70e-5 against 2.77e-5 for first-order injection of the *true* κ₃, so the first-order model has nothing left to give. Network by network, the true κ₃ alone is erratic (1.0e-5 to 7.6e-5). Without its second-order partner it can destabilise the chain, as heisenberg R2–R3 saw at w64. The computed co-state with the scale law is steadier.
- (b) **Joint κ₄ is the whole remaining gap: 15× at w128** (2.77e-5 → 1.80e-6). **Its readout-relevant part is concentrated**: rank 1 off the diagonal (plus the exact diagonal) recovers 4×, rank 8 recovers 8×, rank 32 recovers 11×. The diagonal alone destabilises the chain (1.5e-5).

Both are second-order Duhamel content. C5 says that carrying them exactly costs Θ(L³ n³), while the concentration measured here says the readout needs only a few modes of joint κ₄, as it did for old κ₃. The obvious first candidate is the scale mode's own κ₄. The law-level scale treatment already contains it, but only for the scale variance that the κ₃ projection detects.

## 5. Cost at n = 1024 and grader notes for the best variants

| variant | products | dense f32 | Strassen L3 (0.683 u, ≈ 70 ms/product batched) | wall (L3) | raw | adjusted (L3) |
|---|---|---|---|---|---|---|
| A1gl | 218 | 0.21 B | 149 u = 0.15 B | ≈ 15 s | 6.3e-7 | ≈ 9.2e-8 |
| A2gl | 309 | 0.30 B | 211 u = 0.21 B | ≈ 22 s | 5.0e-7 | ≈ 1.0e-7 |
| A3gsl | 435 | 0.42 B | 297 u = 0.29 B | ≈ 30 s | 3.9e-7 | ≈ 1.1e-7 |
| A2gl_nc | 183 | 0.18 B | 125 u = 0.12 B | ≈ 13 s | 4.3e-7 | ≈ 5.2e-8 |
| A3gl_nc | 231 | 0.23 B | 158 u = 0.15 B | ≈ 16 s | 3.5e-7 | ≈ 5.4e-8 |
| A3gsl_nc | 273 | 0.27 B | 186 u = 0.18 B | ≈ 19 s | 3.25e-7 | ≈ 5.9e-8 |
| all pairs | 855 | 0.83 B | 584 u = 0.57 B | ≈ 60 s | 4.2e-7 | ≈ 2.4e-7 |

The no-coincidence variants use 4 instead of 7 products per pair (U and V propagation, two star contractions). Calls: per layer, one stacked propagation family ([U; V; X] of all live pairs times Φ_k W_{k+1}) and one hub-sum family for the slices, ≈ 2.1k calls in total. The scale-mode bookkeeping is O(n²) elementwise: a handful of calls per layer. Memory: (A + 1) × 3 n×n per live pair, ≤ 100 MB.

## 6. Verdict and the deciding experiment

**Verdict.** 
- (a) **Closure, exact: disproved** for any co-state of o(L n²) numbers per cut that a causal forward pass could carry (C2). The obstruction is the double edge, which closes only in the second chaos (C3). The readout co-state is exactly closed (C1), but evaluating it costs Θ(L² n³) at first order and Θ(L³ n³) at second order (C5).
- (b) **Closure, readout-weighted: established empirically for the old content.** The readout reads old κ₃ through one conserved scalar, the dilation charge (C6, unique by C7), plus a residual that a window of three ages and the slice chain absorb. This turns O(L² n³) into O(A L n³) with A = 1–3 at no loss at n = 1024, and is measured.
- (c) **First order is exhausted.** The exact first order is ≈ 4e-7 raw at 1024, the scale-mode co-state is 3.25e-7, and at w128 the computed co-state equals the true-κ₃ first order. The bar needs the second-order (joint κ₄) content: 15× at w128, readout-concentrated at rank ≈ 8–32.

**Deciding experiment for what comes next.** Run the oracle hierarchy of §4.6 at w256 and w512 on bench seeds: best computed co-state → true κ₃ at first order → + true joint κ₄ → + κ₄ at rank r. Monte Carlo atlases at 1024 are out of reach. At w128 the measured factors are (a) ≈ 1× for the κ₃ model and (b) ≈ 15× for joint κ₄, with rank 8 recovering 8×. The experiment measures how (b) and its rank concentration scale with n:
- If (b) stays ≳ 10× at w512 and rank ≲ 32 keeps most of it, a second-order co-state is the path to raw ≈ 1e-8 at 1024. It would carry a few readout-relevant joint-κ₄ modes, with the scale law already holding the scale part. The modes must then be identified in closed form (candidates: the scale mode's own κ₄, i.e. the bethe Q spike, and the correlated scale field h of §4.3), and their generation computed by δ-insertion at O(r L n³).
- If (b) falls with n, the first-order co-state's ≈ 3e-7 at 1024 is close to what this lens can reach, and the gap to the bar lies elsewhere.

## 7. Files

- `costate.py`: the estimator. Pull-back source atoms (star + exact coincidences); per-pair slices; closures `old = drop | slice | pool | gsm | gsmslice`; `law=True` for the law-level scale mode; oracle filters (`oldfilter = diag | off | rank | gsm | gsm_off | hfit | hfitd`).
- `validate.py`: pull-back = full-tensor HD to 5e-15.
- `run.py`, `summarize.py`, `rebuild.py`: runs and tables (`results/*.jsonl`; rows of the w1024 scale batch rebuilt from saved predictions after a git incident, exact).
- `oracle2nd.py`, `table2nd.py`: the second-order oracle (true MC κ₃/κ₄ slices, own chain, rank filters) and its table. `hybrid.py`: computed co-state + true κ₄ (rank-filtered) on top. `k4anatomy.py`: top κ₄ modes against scale templates.
- `spike.py`, `anatomy.py`: old-slice anatomy. `renorm.py`: the ensemble-renormalisation test (cosine of the slice-chain and true old corrections 0.6–0.7 at w128; one fitted factor does not close the gap).
