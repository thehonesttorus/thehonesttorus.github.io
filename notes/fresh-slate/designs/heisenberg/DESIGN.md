# Heisenberg–Duhamel estimator (design stream `heisenberg`) — DESIGN v0

*Fresh-slate design stream, 1 Oct 2026. Principle: the Heisenberg picture on Gaussian space (Koopman pull-back of the readout, Ornstein–Uhlenbeck / heat semigroup, Stein–Malliavin integration by parts). This file is v0: principle, derivation, estimator, cost at n = 1024, error mechanism, falsification test. Measurements go in RESULTS.md; status labels: **Theorem** (proved here, elementary), **Approx** (a modelling step with a stated error), **Prediction**, **Measured**.*

## 1. Principle, stated precisely

Write the network on pre-activations: z₁ = xW₁, z_{l+1} = ReLU(z_l) W_{l+1}, readout r_j(z_L) = ReLU(z_{L,j}). Let μ_l be the law of z_l (μ₁ is exactly Gaussian). The Koopman (transfer) operator of arrow l acts on observables, (K_l g)(z_l) = g(ReLU(z_l) W_{l+1}), and the **Heisenberg observable** at layer l is g_{l,j} = K_l K_{l+1} ⋯ K_{L−1} r_j, so E_{μ_L}[r_j] = E_{μ_l}[g_{l,j}] for every l (the Schrödinger and Heisenberg pictures agree; the pairing is invariant).

On Gaussian space the reference measure γ_{m,C} carries the OU/heat semigroup, whose generator's Dirichlet form E|C^{1/2}∇f|² is the commutative prototype of the noncommutative Dirichlet forms of the bridges digests. Three of its structural facts are used in an essential way:

- **(S1) Stein / Malliavin integration by parts.** For a law ρ close to γ = γ_{m,C} (same mean and covariance) and any observable g, (ρ − γ)[g] = Σ_{k≥3} (1/k!) ⟨κ_k(ρ), E_γ[∇^k g]⟩ + (products of cumulants) — the Edgeworth/Stein pairing: the *only* things a non-Gaussian law contributes are its cumulant tensors contracted with **Gaussian-averaged derivatives of the observable**. Equivalently E_γ[∇^k g] = ∇_m^k (P g)(m, C) where (Pg)(m,C) = E_{γ_{m,C}} g is the heat-smoothed observable (smoothing commutes with derivatives; ∂_C = ½∇_m²).
- **(S2) Pull-back is linear.** K_l is linear and positive on observables even though the forward map on (moments of) laws is nonlinear. All non-linearity of the problem sits in the choice of reference measures, never in the operator acting on the observable.
- **(S3) Duhamel (variation of constants) for the semigroup of arrows** — the commutative analogue of Chen–Rouzé's telescoping of a time-averaged Lindbladian against the erased state: the error of any reference chain is a sum of *local* defects, each paired with the *exact* Heisenberg observable at that time.

## 2. Derivation

**Reference chain (Schrödinger, closure).** ν₁ = μ₁; ν_{l+1} = Π(T_l ν_l), where T_l is the exact push-forward by one arrow and Π maps a law to the Gaussian with the same mean and covariance. (For Gaussian ν_l the mean and covariance of T_l ν_l are exact and closed-form — bivariate ReLU moments — so ν is computable; it is the "Gaussian closure", here only the *reference*, not the estimator.)

**Theorem 1 (Heisenberg–Duhamel identity, exact).** For every output j,

  E_{μ_L}[r_j] − E_{ν_L}[r_j] = Σ_{l=2}^{L} (ρ_l − ν_l)[g_{l,j}],  ρ_l := T_{l−1} ν_{l−1}.

*Proof.* χ_l := E_{ν_l}[g_{l,j}]; χ₁ = truth (ν₁ = μ₁), χ_L = reference. χ_{l−1} − χ_l = E_{ν_{l−1}}[K_{l−1} g_l] − E_{ν_l}[g_l] = (T_{l−1}ν_{l−1} − ν_l)[g_l]. Sum. ∎

Each term pairs the **local non-Gaussianity** created by one arrow acting on a Gaussian (ρ_l − Πρ_l: cumulants of order ≥ 3 only, of size ε ~ n^{−1/2} per coordinate) with the **exact** pulled-back observable at that layer. No approximation has been made: all "old content" (structure created at layer l and carried to later layers) lives inside g_{l,j}, which contains every later ReLU.

**Theorem 2 (doubly robust remainder).** If ĝ_l is any model of g_l and Δ̂_l := (ρ_l − ν_l)[ĝ_l], then

  truth − (reference + Σ_l Δ̂_l) = Σ_l (ρ_l − ν_l)[g_l − ĝ_l].

The error is **bilinear**: (local non-Gaussianity, size ε) × (model error of the Heisenberg observable's Gaussian-averaged derivatives). Consequently the model of g_l only needs to be good to relative accuracy δ for the correction to be good to ε·δ. (Immediate from Theorem 1.)

**Using (S1).** Δ_l = (1/6)⟨κ₃(ρ_l), E_{ν_l}∇³g_l⟩ + (1/24)⟨κ₄, E∇⁴g_l⟩ + (1/72)⟨κ₃⊗κ₃, E∇⁶g_l⟩ + O(ε³). κ₃(ρ_l) = W_l^{⊗3} κ₃^{G}(ReLU(z_{l−1})), z_{l−1} ~ ν_{l−1} Gaussian: the **source**, computable from (m, C) of the reference alone by Gaussian (Hermite / Mehler / Price) calculus: with ã_p = Σ_{k≥1} α_{p,k} He_k(z̃_p)/k! and ρ the correlation matrix,

  κ₃^{G}(a)_{pqr} = Σ_{i,j,k≥0} α_{p,i+j} α_{q,i+k} α_{r,j+k} ρ_{pq}^i ρ_{pr}^j ρ_{qr}^k / (i! j! k!)

(a sum of *diagrams*: i = j = 1, k = 0 is the "star" centred at p, i = j = k = 1 the "triangle", …; p = q = r the diagonal).

**What the Heisenberg observable sees (kink expansion).** g_l is piecewise linear; its distributional third derivative is a sum over kinks of later neurons (k, p), k ≥ l: δ′(z_{k,p})·β_{k,p}·(u_{l→k,p})^{⊗3} plus two-kink terms δ(z_{k,p})δ(z_{k′,q})(u_p)^{⊗2}⊗u_q, where u_{l→k,p} = ∂z_{k,p}/∂z_l (the gated propagator column) and β the back-propagated sensitivity. Under ν (Gaussian averaging, (S1)) and decoupled gates (**Approx G**: gate ⟂ gate, Φ = P(z > 0)), ⟨κ₃, E∇³g_l⟩ reduces to:

- **diagonal pulled-back cumulants** D_{l→k}(p) = κ₃(ρ_l)[ū_p, ū_p, ū_p] — the third cumulant that the source at l puts on neuron (k, p), feeding E[a_{k,p}] through E δ′;
- **(2,1) pulled-back cumulants** S_{l→k}(p, q) = κ₃(ρ_l)[ū_p, ū_p, ū_q] — feeding Cov(a_{k,p}, a_{k,q}) through E[δ(z_p) H(z_q)] (two-kink terms);

with ū_{l→k} = Φ_l W_{l+1} Φ_{l+1} ⋯ W_k (mean-gate propagator, Φ = diag P(z > 0) under ν). Then everything downstream of layer k is handled by the reference chain's own tangent (a mean/covariance perturbation of a Gaussian, i.e. ĝ = "A exact transport steps, then the closure"). This is the model ĝ_l of Theorem 2.

**Key computational point (the Heisenberg advantage).** S_{l→k} is evaluated by **pulling the slice back to the source layer**: with Ũ = W_l ū_{l→k} (n × n, directions in a_{l−1} space), the star diagram gives

  S(p, q) = Σ_r α_{r,2} [2 Ũ_{rp} R_{rp} R_{rq} + Ũ_{rq} R_{rp}²],  R = ρ_{l−1} diag(α_{·,1}) Ũ,

plus the diagonal Σ_r κ₃(a_r) Ũ_{rp}² Ũ_{rq}: **4 n³ products per (source, target) pair**, no intermediate n³ tensor, no low-rank assumption on old content (the brief's fact that no representation below ≈ 0.3n modes exists is a statement about forward states; the backward evaluation needs exactly the n directions of the target neurons, at full rank).

## 3. The estimator HD-A

State: reference (m_l, C_l); per source layer l and age a ≤ A the transported direction matrix Ũ_{l→l+a}.

Per layer k (forward sweep; linear superposition lets all sources be summed in one sweep):
1. Gaussian quantities of z_k ~ N(m_k, C_k): s, Φ, Hermite coefficients α_{·,0..4}, Eδ^{(i)}.
2. Pulled-back cumulants at k from every source l with k − l ≤ A: D(p) = Σ_l D_{l→k}(p), S(p,q) = Σ_l S_{l→k}(p,q) (and the diagonal κ₄, local, age 0).
3. Stein injection: E[a_p] += D_p Eδ′(z_p)/6 (+κ₄, κ₃² terms); Cov(a)_{pq} += Edgeworth terms with S, Sᵀ, D (bivariate Gaussian expectations of δ, δ′, H — closed form, O(n²)).
4. Reference step: m_{k+1} = E[a_k] W, C_{k+1} = Wᵀ Cov(a_k) W.
5. Transport directions (pulled back to the source): Ũ_{l→k} = W_l Φ_l W_{l+1} Φ_{l+1} ⋯ Φ_{k−1} W_k (columns = directions in a_{l−1}-space whose pairing with a_{l−1} is z_{k,p} under the linearised transport); update Ũ_{l→k+1} = Ũ_{l→k} Φ_k W_{k+1}: one n³ product per live (l, k) pair.

Readout: E[a_{L,j}] from the corrected Gaussian quantities at layer L (with the layer-L injection).

**Cost at n = 1024 (units = 1024³ products; Strassen L5 price 0.556 from notes/streams/costmodel).**

| item | products per layer | units/layer (dense / L5) |
|---|---|---|
| reference sandwich Wᵀ Cov W + mean | 1 sandwich | 1.03 (sym3) |
| per live (l, k) pair: Ũ update 1, R 1, star contractions 2, diagonal 1 | 5 | 5.0 / 2.8 |
| injection, Gaussian quantities, bivariate CDFs | O(n²) | ≈ 0.1 |

Live pairs per layer = A + 1. Total ≈ 16 × (1.1 + 2.8 (A+1)): **A = 0: 63 u (0.06 B); A = 1: 108 u (0.105 B); A = 3: 197 u (0.19 B).** Memory: (A+1) n×n matrices. Flopscope calls: ≈ 120 per Strassen family → keep families fat (batch all live pairs of a layer into one family).

## 4. Error mechanism and predicted scaling

Theorem 2 gives the error exactly as Σ_l (ρ_l − ν_l)[g_l − ĝ_l]. The model ĝ_l is wrong through:
1. **Age truncation** (sources older than A not transported); the brief's fact: 40 % of the (2,1) structure is older than one layer; old content concentrates on top propagator singular directions with participation ratio ≈ n/(2·age).
2. **Gate decoupling (Approx G)** — joint gate law has spectral independence η ≈ 2–8, so E[H_pH_qH_r] ≠ Φ_pΦ_qΦ_r. Correction: pairwise gate law E[H_pH_q] (arcsine/bivariate CDF, O(n²)) in the transport of S.
3. **Missing diagrams**: triangle (α₂³ρ³) and higher ρ-orders in the source; κ₃ transported through a ReLU with δ-insertions; joint κ₄.
4. **Edgeworth truncation**: O(ε³) per layer.

*Prediction.* With the reference closure at raw MSE M₀ ≈ c₀/n (measured ≈ 4e-5 at n = 1024), and every truncation complete, the remainder is O(ε²) per layer: raw MSE ≈ c₂/n², i.e. a factor ≈ n/c better than closure — at n = 1024 a gain of 10³ is needed to reach ≈ 4e-8, so the truncations must leave ≲ 3 % of the first-order correction in rms.

## 5. Cheapest falsification test

At widths 64/128/256, depth 16, compute three things against baked truth:
- (F1) the reference closure error (M₀);
- (F2) **HD-∞ ceiling**: the same pairing with the *full* source tensor κ₃^{G} (all diagrams, Hermite series), all ages, transported as full n³ tensors with decoupled gates (feasible at n ≤ 256);
- (F3) HD-A with the cheap diagrams (diagonal + star), A = 0, 1, 3.

The design is dead if F2/F1 is not ≲ 1/30 at n = 256 *and* improving with n (the ceiling of the realisation, not of the principle: a failure would be charged first to Approx G and to the κ₃-only Edgeworth order, which F2 does not remove). It is uncompetitive if F3 at the affordable A is not within 2× of F2.

Toy exactness check (width 8–16, depth 2–4): verify Theorem 1 term by term — each (ρ_l − ν_l)[g_l] computed by high-precision Monte Carlo (large N, exact network after layer l) must sum to truth − reference; and check the κ₃-Stein pairing against the exact term.

## 6. Status after Stage Q (v1, same day) — verdict

**What was confirmed.**
- Theorem 1 holds numerically (R0), and the doubly robust structure behaves as predicted. With first-order Stein pairing, full κ₃ sources and all ages (HD-∞), the gain over the reference grows with width: 4.6× (64), 8.9× (128), 7.6× (256, mean of 3; 4.5–13.7×). The HD raw MSE falls as n^{-1.9 to -2.1}, against roughly n^{-1} for the reference, as §4 predicted: O(ε²) remainder against O(ε).
- The source-layer pull-back realisation is sufficient for first order: the diagonal + star diagrams lose nothing against the full Hermite series, and 8 ages keep 92 % of the gain (R7). So first order at n = 1024 costs, with Strassen L5 pricing (kit: 0.557 u per product), about 92 (source, target) pairs × 5 products × 0.557 + 16 × 1.03 (reference) ≈ **270 units ≈ 0.26 B**. A low-rank pull-back of old ages (participation ratio ≈ n/(2·age)) could bring this to ≈ 120–150 u; that has not been measured. Flopscope calls: ≈ 92 Strassen families × 79 calls ≈ 7.3k, at the residual limit. Batching each layer's pairs into one family (16 families ≈ 1.3k calls) is required.
- The oracle chain (R2, R3) locates the remainder exactly: after first order, the error is covariance drift caused by the joint fourth-cumulant slices κ₄(pppq), κ₄(ppqq). With them (true values), the gain over first order is 7–20× at width 64.

**Projection at n = 1024 (first-order HD).** Raw = 7.4e-6 × (256/1024)^p with p = 1.9–2.1: **raw ≈ 4–7e-7 (range including network scatter 2e-7 – 1e-6)**. At 0.26 B: **adjusted ≈ 1.0–1.8e-7**, or 5–9e-8 if the low-rank old-age compression reaches 0.13 B. This is about 30× better than the Gaussian reference (≈ 4e-6 adjusted) and 7–20× better than plain sampling (1.2e-6), and 30–100× short of the 1.6e-9 bar.

**What would make it competitive.** The second-order Duhamel term with *old* joint κ₄: per (source, target) pair, the (3,1) and (2,2) slices of κ₄ that κ₃ becomes when it passes through later ReLUs (δ-insertions) and that κ₃ ⊗ κ₃ generates. The R3 oracle bounds the gain at 7–20× (width 64) over first order, and its relative weight should fall like ε with width. A rough optimistic projection is raw ≈ 3–6e-8 at 1024, adjusted ≈ 1e-8 at ≈ 0.3 B: still above the bar unless the remaining third-order drift is also controlled. Honest verdict: **the Heisenberg–Duhamel principle organises the problem correctly and produces a cheap, provably second-order-accurate first stage, but as realised here it is not competitive. Reaching the bar needs the full second-order (κ₄-transport) layer, at which point the computation coincides with the adjoint form of a second-order cumulant chain.** The principle-specific gains are the exact error accounting, age-resolved source pull-back without forward tensors or low-rank assumptions, and the bilinear remainder. None of them changes the order of the needed object.

**Deciding experiment (as planned in v1).** At width 256, inject the true (MC) κ₄(pppq), κ₄(ppqq) on top of the computed first-order HD (own chain). If the gain over first order is < 5× there, the second-order layer cannot bring the design to the bar, and the stream should stop. If it is ≥ 10×, build the δ-insertion / κ₃⊗κ₃ → κ₄ transport in pull-back form: per pair, the slices are κ₃-source diagrams evaluated on transported direction pairs, at ≈ 3–4 products per pair.

**Deciding experiment, result (R9–R10).** At widths 64 and 128, adding the *true* joint κ₄ slices on top of computed first-order HD makes it **worse** (gain 0.27–0.83). The cause has been located: the transported κ₃ slices carry 22–31 % error from the first transport step on (decoupled gates, no δ-insertions). The source evaluation is exact. So the order of work, if this stream continues, is fixed:
1. exact first-order linear response of κ₃ through a ReLU: joint gate law plus δ-insertions. In pull-back form these are additional diagrams on the transported direction pairs, at ≈ 2 more products per pair. Target: slice error ≲ 3 % at width 64.
2. Only then old joint κ₄ (κ₃ → κ₄ through later ReLUs, κ₃ ⊗ κ₃), whose ceiling is R3's 7–20×.

**Final verdict (this session).** HD first order at n = 1024: projected **raw ≈ 5e-7 (2e-7 – 1e-6), adjusted ≈ 1.3e-7 at ≈ 0.26 B (5e-8 if old ages compress to 0.13 B)**. That is not competitive with 1.6e-9. A competitive version needs steps 1 and 2 above. Both are diagrammatic pull-back extensions of the same structure, and the R3 ceiling (true κ₃ and κ₄ on own chain, 7–20× over first order at width 64) suggests the result would land near 1e-8 adjusted, still short unless the remaining (third-order) covariance drift also scales away fast. The theoretical content that survives is Theorems 1–2 (exact Duhamel accounting and the bilinear remainder) and the source-layer pull-back. The pull-back evaluates old content at full rank with n³ work per (source, target) pair and no forward tensors, which is exactly what the brief's "no low-rank representation below 0.3 n" fact makes hard for forward designs.
