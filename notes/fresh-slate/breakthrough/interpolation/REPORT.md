# Interpolation in the disorder — breakthrough stream REPORT

*Stream `interpolation`, 1 Oct 2026. Lens: the scored quantity F_j(W) = E_x[a_{L,j}] is quenched (one network), but W is
i.i.d. Gaussian, so Gaussian integration by parts **in the weights** (smart paths, Lindeberg, cavity, Bolthausen
conditioning, covers) is available. Labels: **Theorem** (proved here), **Measured** (bench `w1024_d16`, 6 MLPs, truth
noise 3.6e-8 subtracted), **Prediction**. All numbers are raw final-layer MSE at n = 1024 unless stated.*

v0 (≈ 2 h in). Code in this directory; results in `results/`.

## 0. Summary

1. **(a) Smart paths in the weights do not close for a quenched quantity — but they split it.** The interpolation
   W_t = e^{−t}W + √(1−e^{−2t}) W̃ is the Ornstein–Uhlenbeck semigroup in the weights; the exact identity it yields is
   Stroock's chaos formula F(W) = Σ_p (1/p!)⟨E[D^p F], H_p(W)⟩: every kernel is *annealed* (an ensemble average,
   computable in principle by mean-field theory), every quenched bit sits in Hermite polynomials of the actual W, and
   the integrand of the path is as hard as F. What the lens does give is a split of any estimator's error into a
   **self-averaging part** — the disorder-conditional expectation E_W[error | cheap per-neuron quenched statistics], a
   deterministic function that is the same for every network of the ensemble — and a genuinely quenched remainder.
   **Measured: 58 % of the Gaussian closure's final-layer MSE is self-averaging** (a smooth function of each neuron's
   own α_j = m_j/s_j and s_j, fitted on 5 networks, removes 3.39e-6 → 1.41e-6 on the 6th, leave-one-MLP-out). Cost of
   the correction at run time: O(n).
2. **(b) Bolthausen conditioning: no Onsager memory in a feedforward net with independent layers.** The conditional law
   of W_{l+1} given everything computed from W_{≤l} is its prior, so there is no AMP reaction term; the old content of
   the dossier is not an Onsager memory but quenched transport. Conditioning W_{l+1}'s columns on the statistics the
   estimator does reveal (m_j = μ·w_j, v_j = w_jᵀCw_j) gives the Bayes-optimal readout over the disorder given those
   statistics: a Gaussian mixture over per-input **order parameters** (λ(x) = a·μ/|μ|², τ(x) = |a − μ|²/E). Measured:
   Var τ_l grows linearly with depth, ≈ 2.1 × (2/n) per layer (35 × chi at layer 16, s.d. 0.26), layer-to-layer
   correlation 0.97, about twice what the Gaussian closure's own covariance predicts. **But the mixture readout is no
   better than the Gaussian one at exact moments** (rms one-step defect 4.3e-4 Gaussian vs 4.6–5.5e-4 mixture at layer
   16). Part of the reason is a double count in that test (it mixed with the *total* norm variance, whose Gaussian part
   2‖C‖² is not kurtosis); the correct chaos-0 kurtosis uses the excess X only and gives a modest gain (item 7). The
   lasting use of the order parameter is item 5: its covariance with each neuron is the trace channel.
3. **(c) Covers: precisely negative.** A lift that shares the input is function-identical to the base network (each
   copy receives identical pre-activations by induction), so it carries no information; a lift that also lifts the
   input thins every loop through the input with probability ½ per 2-lift, and the tower converges to the tree, whose
   mean is exactly computable (independent neurons; cumulants add) — **measured 2.05e-4, 48× worse than the Gaussian
   closure (4.30e-6)**. The quantities that matter are exactly the loops a cover cuts.
4. **Free exact corrections found on the way.** The radial factor of positive homogeneity (E a_L = E r · E_{S^{n−1}} u_L,
   r ~ χ_n): running the closure on the angular network multiplies all means by E r/√n = 1 − 1/(4n) + …, which lowers
   the closure from 4.30e-6 to **3.53e-6 (−18 %, 6/6 MLPs)**; at exact moments the one-step Gaussian readout's mean
   defect (−2.0e-4 at layer 16) matches the radial factor (−2.2e-4). Gaussian-exact bivariate ReLU covariance (Hermite
   order 2 instead of the linearised one) adds 3.35e-6.

5. **The breakthrough-shaped finding: old content of the per-neuron third cumulant is carried by ONE vector per layer.**
   Expanding κ3(z_{l+1,p}) = κ3(a_l)[w_p^{⊗3}] in Wiener chaos of the column w_p, the chaos-1 part is
   3σ²(Wᵀt_l)_p with the **trace channel** t_l = E[|ã_l|²ã_l] = Cov(|ã_l|², a_l) — the covariance of every neuron with
   the per-input order parameter of (b). **Measured at n = 1024: chaos 1 explains 92–93 % of the variance of the
   per-neuron κ3 from layer 8 on (corr 0.965; MC noise 6 % of the signal).** And t obeys a closed recursion derived
   from first-order Edgeworth plus the chaos-0 contraction of the transport Gram (no fitted coefficients):
   t_{l+1} = src_G + ½σ²(Σγ)·Φ∘(Wᵀt_l) + ½σ²(t_l·Wu)·φ/s, **R² = 0.990–0.998 per layer** against MC (oracle input).
   Cost: O(K n²) per layer on top of the closure.
6. **Trace-channel closure (TC, `tc.py`): raw 1.45e-6 (1.31e-6 with the chi factor) at n = 1024, 6/6 MLPs**, against
   4.10e-6 for the order-2 Gaussian closure and 4.30e-6 for the bench closure; cost ≈ the closure's (≈ 17 units ≈
   0.017 B, multiplier 0.1) → **adjusted ≈ 1.3e-7**: the level of the best fresh design (HD first order, projected
   1–1.8e-7 at 0.26 B) at a fifteenth of its cost, and ≈ 80× above the bar.
7. **What is left is quenched.** A disorder-conditional (self-averaging) correction trained on 23 fresh networks gains
   nothing on top of TC (1.31 → 1.28e-6): the closure's self-averaging error *was* the trace channel. The scalar κ4
   channel (chaos 0 of κ4(z_p): 3σ⁴X, X = excess norm variance, oracle value) adds 1.29 → 1.02e-6 (MLP 0). Injecting
   the Monte Carlo t instead of the propagated one gives 1.51 → 1.37e-6 (MLP 0): the propagation of t is not the
   bottleneck either. The remainder is the chaos ≥ 2 content (the 7 % of κ3 outside the trace channel, the
   off-rank-one part of the (2,1) slice, the κ4 slices through a matrix-valued trace M = Σ_i κ4(a)_{ii··}).

**Verdict v1.** (a), (b), (c) do not give an exact O(Ln³) representation of the quenched mean; each fails for a stated
reason. But the lens's own decomposition — **Wiener chaos in the weights** — identifies which part of the old content is
low-dimensional: the chaos-1 trace of κ3, one n-vector per layer, with an exact-to-1 % O(n²) recursion. That turns
the Gaussian closure into a 3.1–3.3× better estimator at the same cost. It is not a path to 1e-8 on its own: the
remaining error is quenched chaos ≥ 2 content, which is where the expensive designs live. The useful export to the
other streams is (i) the trace channel, which any chain can carry for free (and should subtract before compressing
old content — it is the rank-one "spike" of the (2,1) slice), and (ii) the chi factor.

## 1. The task in the lens

a_0 = x ~ N(0, I_n), z_l = a_{l−1}W_l, a_l = ReLU(z_l), W_l i.i.d. N(0, σ²), σ² = 2/n, layers independent. Target
F_j(W) = E_x[a_{L,j}] for the given W. Per-neuron quenched fluctuations of F_j around anything annealed are O(1)
(the mean direction m_j = μ·w_j alone is O(1)); the bar needs ≈ 1e-4 rms.

## 2. (a) Smart path = OU semigroup in the weights

**Theorem A (interpolation identity).** Let W̃ be an independent copy, W_t = e^{−t}W + √(1−e^{−2t})W̃ and
P_tF(W) = E_{W̃}F(W_t). Then P_0F = F (quenched), P_∞F = E_W F (annealed), ∂_tP_tF = 𝓛P_tF with
𝓛 = Σ_{l,i,j}(σ²∂²_{W_{l,ij}} − W_{l,ij}∂_{W_{l,ij}}), and

  F(W) = E F + Σ_{p≥1} J_pF(W),  J_pF(W) = (1/p!) ⟨E[D^pF], H_p(W)⟩_σ,  P_tF = Σ_p e^{−pt} J_pF.

*Proof.* Mehler's formula and Stroock's formula for the Gaussian measure on the weights (finite-dimensional). ∎

Consequences, stated as facts about this problem:
- Stein (Gaussian IBP) applies only to the averaged copy W̃. The smart-path derivative
  d/dθ E_{W̃}F(√θW + √(1−θ)W̃) = ½E[⟨∇F, W⟩/√θ − σ²ΔF] keeps the quenched term ⟨∇F(W_θ), W⟩, which needs the same
  joint law over x as F itself. No solvable endpoint (infinite width, fresh weights per path, orthogonal or diagonal
  weights, a tree cover) turns this into a cheap integrand: the quenched content is the whole of J_{≥1}F.
- **Chaos 1 of the readout is exact in the closure.** For the last layer, with the exact law of a = a_{L−1} fixed,
  F_j(w) = E_x ReLU(a·w) and the closure G(w) = ψ(μ·w, wᵀCw). E_w∇F = E_x[a·P(a·w > 0)] = ½μ and E_w∇G = ½μ exactly
  (w → −w symmetry), so J_1(F − G) = 0: the closure's readout error starts at chaos 2 (radial/norm structure,
  E_x[aaᵀ/|a|] vs the Gaussian kernel) and chaos 3 (κ3(a)[w^{⊗3}], the quenched third cumulant).
- **The useful split.** For per-neuron quenched statistics f_j (cheap: α_j = m_j/s_j, s_j, …), write the error of a
  reference as e_j = ē(f_j) + ẽ_j, ē(f) = E_W[e_j | f_j = f]. ē is a deterministic function, identical for every
  network of the ensemble (exchangeability of neurons and of networks), so it is a property of (n, L, reference) and
  can be learned once, offline, from simulated networks, or derived from 1/n field theory. ẽ_j is the quenched
  remainder. This is the precise sense in which "interpolation in the disorder" helps a quenched estimate: it moves
  the self-averaging part out of the run-time computation.

**Measured (T3, `t3_selfavg.py`).** Closure (K = 2, chi factor), final-layer error regressed on
{s·α^k, s·φ(α)α^k, α^k}, k ≤ deg, leave-one-MLP-out over the 6 bench networks:

| features | raw MSE (LOO) | per MLP |
|---|---|---|
| none (closure + chi + K2) | 3.39e-6 | — |
| deg 0 (bias, s, sφ) | 2.49e-6 | 1.95–2.99e-6 |
| deg 1 | 1.53e-6 | 1.20–2.03e-6 |
| deg 2 | **1.41e-6** | 1.08–1.91e-6 |
| deg 4 | 1.41e-6 | 1.08–1.92e-6 |

The fit has ≤ 15 coefficients on 5120 neurons; it transfers between networks, i.e. it is self-averaging, as the
theory says. (More training networks are being baked: `gen_train.py`, 24 fresh seeds at 1024, N = 2.6e5.)

## 3. (b) Bolthausen conditioning and AMP

**Theorem B1 (no Onsager term).** Let 𝓕_l = σ(W_1, …, W_l). Any estimator state computed from W_{≤l} is 𝓕_l-
measurable, and W_{l+1} ⟂ 𝓕_l, so Law(W_{l+1} | 𝓕_l) = N(0, σ²)^{⊗n×n}. In AMP the Onsager term arises because the
same matrix is applied to iterates that depend on it; in a feedforward net with independent layers each matrix acts
once, on a law that does not depend on it. Hence the infinite-width state evolution has no memory term, and **the
old content of the dossier is not an Onsager memory**; it is quenched transport by the fixed matrices (every product
W_l Φ_l W_{l+1}⋯ is computed explicitly by any chain; conditioning adds nothing).

**Theorem B2 (MMSE readout given revealed statistics).** Condition column w = w_j on m = μ·w and v = wᵀCw. Given x,
a·w is Gaussian to O(n^{−1/2}) with mean λ(x)m, λ = a·μ/|μ|², and variance
σ²|P_μ^⊥a|² + (ãᵀCã/‖C‖_F²)(v − σ² tr C) (regression of (a·w)² on wᵀCw). So E[F_j | m_j, v_j] is a mixture of
ψ(λm_j, ·) over the joint law of the per-input **order parameters** (λ(x), |ã(x)|², ãᵀCã) — the finite-n
fluctuating order parameters of mean-field theory; at infinite width they are constant and the readout is ψ(m_j, v_j).

**Measured (T1, `t1_norm.py`, MLP 0, N = 65 536).** τ_l = |a_l − μ_l|²/E:

| layer | 0 (|x|²) | 1 | 4 | 8 | 12 | 16 |
|---|---|---|---|---|---|---|
| Var τ in units of 2/n | 1.01 | 3.25 | 8.88 | 15.6 | 23.5 | 35.1 |
| closure's own 2‖C‖²/tr²C (units 2/n) | — | 1.5 | 4.2 | 7.8 | 13.9 | 18.5 |
| corr(τ_l, τ_{l+1}) | 0.56 | 0.79 | 0.91 | 0.94 | 0.96 | — |

The order parameter is a slow, strongly persistent latent (a single scalar carrying "old content"), with about half
of its variance non-Gaussian (Σ_pq κ4(ppqq) excess). One-step readout defect at layer 16 with the exact first two
moments of z (on the empirical law): rms 4.3e-4 (mean −2.0e-4) Gaussian; 5.5e-4 pure scale mixture; 4.6e-4 mixture
plus coherent shift β_j(τ − 1). **The MMSE-over-disorder shape does not beat the Gaussian shape**, because the
quenched projection a·w_j sees the norm fluctuation only through its overlap with the few directions that carry it.
The readout is not where the closure loses: exact first two moments + Gaussian readout ≈ 2e-7 MSE, so ≈ 95 % of the
closure's 4.3e-6 is moment drift through the chain.

## 4. (c) Covers and towers of covers

**Theorem C1.** In a k-lift of the layered graph that does not lift the input (every copy of an input node is x_i),
every copy of neuron (l, j) has the same pre-activation, by induction on l (Σ_i over the lifted neighbours of copy c
of j is Σ_i x_i W_ij for every assignment of copies). The lift computes the base network exactly: no information.

**Theorem C2.** If the input is lifted with independent copies, a 2-lift with signs S_ij keeps a loop through the
input with weight depending on the sign product; on average the covariance of layer-1 pre-activations of the copies is
thinned (off-diagonal pairs with S_ij ≠ S_ik drop out), and as k → ∞ the cover is the tree: independent neurons,
z_j = Σ_i w_ij a_i with independent a_i, exactly computable (cumulants add; κ3(z_j) = Σ w³κ3(a_i) = O(n^{−1/2})).

**Measured (T0, `t0_baselines.py`).** Tree cover (diagonal closure) 2.05e-4; annealed variance (state evolution with
the quenched mean only) 2.61e-4; Gaussian closure 4.30e-6. Covers remove precisely the quenched loops that carry the
answer (48×). The Bethe-on-covers principle does not apply to a network whose computation is a deterministic function
of a shared input.

## 5. The trace channel: chaos 1 of the third cumulant, and the TC estimator

**Theorem T (chaos split of a projected cumulant).** For z = aW with columns w_p i.i.d. N(0, σ²I) independent of the law
of a, and κ = κ3(a), the Hermite (Wick) decomposition of w_iw_jw_k gives, exactly,

  κ3(z_p) = κ[w_p, w_p, w_p] = 3σ² (Wᵀt)_p + κ[:w_p w_p w_p:],  t_c = Σ_i κ_iic = E[|ã|²ã_c] = Cov(|ã|², a_c),

and for the (2,1) slice, a ≠ b: κ3(z_a, z_a, z_b) = σ²(Wᵀt)_b + (chaos 2 in w_a) ⊗ (chaos 1 in w_b). The chaos-1
parts are the only parts that survive averaging over all columns but one; they are carried by the single vector t.
The (2,1) slice's chaos-1 part is rank one (all rows equal Wᵀt): it is the "spike" that the dossier saw emerging with
depth, and it should be removed before any low-rank compression of old content. ∎

**Measured fraction (t5, t9; MLP 0, n = 1024, exact-centred two-pass MC, N = 131 072).**

| z layer | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 |
|---|---|---|---|---|---|---|---|---|
| var(κ3(z_p)) explained by 3σ²(Wᵀt)_p | 0.44 | 0.78 | 0.87 | 0.90 | 0.92 | 0.92 | 0.93 | 0.92 |

**Recursion (derived; `t6_recursion.py`).** With f_a = (ReLU(z_a) − μ_a)², g_c = ReLU(z_c), first-order Edgeworth in
κ3(z) for Cov(f_a, g_c) (pure terms cancel against marginals) and the chaos-0 contraction of the transport Gram
W diag(γ) Wᵀ ≈ σ²(Σγ)I:

  t_{l+1,c} = Σ_a Cov_G(f_a, g_c) + ½σ²(Σ_aγ_a) Φ_c(Wᵀt_l)_c + ½σ²(t_l·W u) φ_c/s_c,
  γ_a = E f_a'' = 2Φ_a − 2μ_aφ_a/s_a,  u_a = E f_a' = 2μ_a(1 − Φ_a),

with the Gaussian source Σ_a Cov_G(f_a, g_c) = Σ_k Σ_{a≠c} E[f_a He_k]E[g_c He_k] ρ_ac^k/k! + E_G[ã_c³] (Mehler,
O(Kn²)). Against MC t_l as input, R² of t_{l+1} = 0.94 (layer 2), 0.98–0.99 (3–5), 0.993–0.998 (6–16); the fitted
coefficient of the transport term exceeds the derived one by 3–6 % per layer. Self-consistently propagated, t has
correlation 0.98–0.995 with the exact-centred MC t at every layer but its scale drifts to 0.71× by depth (the 3–6 %
per layer compounding: an O(ε) term missing from the transport coefficient, not yet identified). Injecting the MC t
instead improves the final MSE by only 9 % (1.29 → 1.18e-6), so this is not the bottleneck.

**Estimator TC (`tc.py`).** Gaussian closure with Hermite-order-2 covariance; per layer y = Wᵀt_{l−1};
κ3(z_p) = 3σ²y_p injected (first-order Edgeworth) into E a = ψ − κ3 αφ/(6s²) and E a² (+κ3 φ/(3s)); the chaos-1 (2,1)
slice σ²y_b injected into Cov(a) as the rank-2 update ½σ²[(φ/s) ⊗ (Φy) + (Φy) ⊗ (φ/s)]; then t_l by the recursion.
Cost over the closure: O(Kn²) per layer (one matvec Wᵀt and the Mehler source), i.e. nothing at n = 1024.

| set | closure (K=2) | TC mean only | **TC** | TC + chi | gain |
|---|---|---|---|---|---|
| w64_d16 (8) | 4.50e-4 | 8.64e-4 | 4.53e-4 | 4.58e-4 | 1.0× |
| w128_d16 (8) | 2.82e-4 | 2.65e-4 | 1.15e-4 | 1.09e-4 | 2.5× |
| w256_d32 (2) | 1.21e-4 | 1.49e-4 | 5.28e-5 | 5.07e-5 | 2.3× |
| **w1024_d16 (6)** | 4.10e-6 | 3.60e-6 | **1.45e-6** | **1.31e-6** | 3.1× |

(The bench's linearised closure is 4.30e-6 at 1024.) The covariance injection is essential (mean-only injection is
worse than nothing at small width). Per-layer MSE, MLP 0: closure grows 4e-7 → 5.0e-6 over 16 layers; TC 3.3e-7 →
1.6e-6 — the trace channel removes most of the depth growth.

## 6. Next (in progress)

1. **Locate TC's remainder** (`t13_hiN.py`, N = 2M): oracle per-neuron z variances vs oracle κ3 in the TC chain.
2. **Is the non-trace 7 % of κ3 young?** (`t14_age0.py`): correlation of the chaos ≥ 2 residual with the age-0
   Gaussian source. If yes, TC + the age-0 star diagram (≈ 3 products per layer) completes the per-neuron κ3.
3. **The matrix trace of κ4.** κ4(z_p,z_p,z_q,z_q) and κ4(z_p,z_p,z_p,z_q) have chaos-0/2 parts carried by the scalar
   X = Σ_ij κ4_iijj (tested as oracle: 1.29 → 1.02e-6) and the matrix M = Σ_i κ4(a)_{ii··} = κ3(τ, a, a) − 2C²,
   which a chain can carry with one extra sandwich per layer (≈ 1 unit).

## Files

`common.py` (helpers), `t0_baselines.py` (closure, chi factor, tree cover, annealed variance), `t1_norm.py`
(order parameters, mixture readouts), `t2_closures.py` (Hermite-order-K Gaussian closure), `t3_selfavg.py`
(self-averaging split), `gen_train.py` (training networks), `tc.py` (trace-channel closure), `t4`–`t14` (tests, see headers). Results: `results/*.json|log`.
