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
   16): the norm fluctuation is anisotropic (it lives in the collapsed top covariance directions, which the closure's C
   carries), not isotropic as annealed conditioning assumes. Negative, precise.
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

**Verdict v0.** The lens does not produce an O(Ln³) exact representation of old content (a, b and c are negative on that
question, with reasons). It does produce a cheap, general, additive tool: **disorder-conditional (self-averaging)
correction of any reference chain, at O(n) run-time cost, trained offline on the known weight ensemble.** On the
Gaussian closure (≈ 17 units ≈ 0.017 B, multiplier 0.1) it gives raw ≈ 1.4e-6 → adjusted ≈ 1.4e-7, which already
matches the best fresh design's projection (Heisenberg–Duhamel, ≈ 1–1.8e-7 adjusted at 0.26 B) at a fifteenth of the
cost; it is still ≈ 100× above the bar on its own. The deciding questions (§5) are whether the self-averaging fraction
grows when the correction is applied *inside* the chain at every layer (drift is the dominant error) and whether it
survives on top of a first-order (κ3) reference.

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

## 5. Next (in progress)

1. **In-chain self-averaging correction (the deciding experiment for this stream).** Apply ē_l(f) at every layer
   (to μ_l, and to the variance/diagonal of C_l), each fitted on training networks with teacher-forced targets, so
   the correction removes drift where it is born instead of at the readout. Prediction: if the self-averaging fraction
   of drift is ≈ 60 % per layer as at the readout, the final MSE falls well below 1.4e-6.
2. **On top of a κ3 reference.** Does the residual of a first-order (old-content) chain keep a large self-averaging
   fraction? If yes, the corrections compose (orthogonal error parts); if no, self-averaging is exhausted by the
   Gaussian level.
3. **Theory of ē.** The chaos-0 parts of the κ4-driven corrections are traces (Σ_pq κ4(ppqq) is exactly the order-
   parameter excess measured in T1), i.e. scalar recursions at O(n²) per layer: a derived, not fitted, ē.

## Files

`common.py` (helpers), `t0_baselines.py` (closure, chi factor, tree cover, annealed variance), `t1_norm.py`
(order parameters, mixture readouts), `t2_closures.py` (Hermite-order-K Gaussian closure), `t3_selfavg.py`
(self-averaging split), `gen_train.py` (training networks). Results: `results/*.json|log`.
