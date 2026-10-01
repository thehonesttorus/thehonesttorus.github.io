# Breakthrough stream "region": the feasible region at n = 1024, and a design inside it

*Status: v2, final for this session (2 Oct 2026, ≈ 02:30 UTC). Part 1 is measured on bench MLPs 0–2 of w1024_d16.
Part 2 is a self-contained first-order Wiener-chaos estimator ("FC"), measured on all 6 bench MLPs at n = 1024 and
paired with the exact Gaussian closure. Every number is at n = 1024, depth 16, on the shared bench networks.*

**Summary.**

1. **The books close.** At n = 1024, the exact Gaussian closure's final MSE is the incoherent sum over layers of its
   local errors times measured transfer coefficients K(l), to 2–6 % (46 % on one MLP, where errors cancel between
   layers).
   - ≈ 70 % of that MSE is the non-Gaussian correction of the off-diagonal covariance.
   - ≈ 30 % is the per-neuron mean readout.
2. **Necessary conditions** (§4). Any estimator at the bar must carry:
   - a dense covariance arrow, accurate to ≈ 0.3 % of ‖C_off‖_F at layers 6–12;
   - D21 at every layer to ≲ 5–10 %, including content of every age (dropping ages > 4 costs 25× the bar);
   - the joint fourth cumulant **only through 2 n-vectors per layer**: the diagonal κ4 and the column means of the
     (2,2) slice. These must be within ≈ 30 %: × 0.7 costs nothing, × 0.3 costs 3×, and omitting the column means
     costs 6×;
   - the diagonal κ3 and κ4 for the readout.

   Not necessary: the n² κ4 slices, the (3,1) slice, second-order Edgeworth on the covariance, coherent directions,
   and early-layer accuracy.
3. **Sufficient construction** (§5). FC carries first-order chaos with two-time cross-covariances, every age, exact
   slices, and a new O(n²)-per-layer mean-field recursion for the κ4 vectors. The fresh-weight lemma makes the
   all-distinct κ4 drop out of every index sum. It scores raw **3.0e-8 ± 0.35e-8 on all 6 bench MLPs** (exact Gaussian
   closure: 4.10e-6, ≈ 135× worse). Its D21 error compounds from 2.5 % to 10 % with depth; the second-order birth
   diagrams are the identified next lever.
4. **Cost verdict.** Dense, FC costs ≈ 840 units (0.82 B), ≈ 0.5 B with Strassen. Adjusted that is ≈ 1.5–2.5e-8, so
   it does not compete. Every cheap form of the old content was measured and fails at n = 1024:
   - a per-source rank-256 row basis: 4× worse;
   - a shared Oseledets subspace q = n/4: 6× worse (q = n/2 is lossless but no cheaper);
   - Wick-atom pruning to n/4: 100× worse;
   - slice-atom pruning to n/2: 2.5× worse.

   Each source is an atom-complete sum: all n atoms carry comparable energy. In a bilinear (matrix-product)
   realisation the all-age D21 then costs ≥ ≈ 4 products per (source, target) pair, i.e. ≥ 480 units (§6). The leaders'
   ≤ 154 units imply they do not compute this object atom by atom. §6 states the remaining open problem precisely.

Files: `gclose.py` (exact bivariate-Gaussian closure), `lr_probe.py` (transfer coefficients), `mc1024.py` +
`mc_analyse.py` + `price.py` (n = 1024 Monte Carlo atlas in two independent halves, noise-corrected norms, pricing),
`fc.py` (the Part 2 estimator), `diag_fc.py`, `check_b0.py`, `check_birth.py`, `check_slices.py` (layer-by-layer
diagnosis against the atlas). Results: `results/` (`lr_mlp*.json`, `mc_*.json`, `price_*.json`, `diag_fc_*.txt`).
Atlases (`data/`, 0.7 GB each) are not committed; they regenerate with `python mc1024.py MLP 262144` (≈ 7 min, 3 threads).

## 0. The bar in this stream's units

- **Accuracy.** Raw ≤ 1.5e-8 means a final-layer per-neuron error ≤ **1.22e-4 rms**. Final-layer activations have
  standard deviation 0.24–0.30, so this is ≈ 4e-4 of a standard deviation.
- **FLOPs.** Cost ≤ 0.11–0.15 B = 113–154 units, i.e. **7–10 units per layer**.
- **Residual time.** ≤ 0.4 s ≈ 6–8k flopscope calls with margin, i.e. **≈ 400–500 calls per layer**.
- **Calibration.** The exact bivariate-Gaussian covariance closure (`gclose.py`, no linearisation) scores raw
  **4.10e-6 ± 0.33e-6** on the 6 bench MLPs; the bench's linearised baseline scores 4.30e-6. Linearising the
  cross-covariance costs only 5 %. The bar is ≈ 270× below the exact Gaussian closure.
- **Truth noise.** The bench truth noise at the final layer is 3.6e-8, above the bar, and is subtracted. At earlier
  layers the per-layer truth noise is var(a_l)/N. It dominates layers 0–6, so all-layer MSEs are not informative there.

## 1. A structural lemma: fresh weights make Frobenius norms the currency

W_{l+1} is independent of the law of a_l and of anything an estimator computes from W_1..W_l. So any error δ in a
layer-l object enters layer l+1 as a Gaussian chaos in the fresh columns w_c of W_{l+1}. The second moments of that
chaos depend only on rotation invariants of δ, taken jointly with the true state.

- **Covariance error.** For a symmetric error δC in C(a_l), δσ²_c = w_cᵀ δC w_c has mean (2/n) tr δC and variance
  2 (2/n)² ‖δC‖_F². The off-diagonal image of δC is Wigner-like, with entry variance (2/n)² ‖δC‖_F².
- **Third-cumulant error.** For an error δκ in a third cumulant, δκ3(z_c) has mean 0 and variance ≈ 6 (2/n)³ ‖δκ‖_F²,
  and ‖δD21(z_{l+1})‖_F² ≈ (16/n) ‖δκ‖_F².

**Consequence:** an error that is incoherent with the true state costs a final-MSE price proportional to its squared
Frobenius norm, with a per-layer coefficient that does not depend on its shape. §2 checks this with coherent probes. The
lemma has a corollary we did not expect (§4, N7): structure that is coherent in *index sums* (row sums of slices) is
transported O(√n) more strongly than incoherent structure. That is where the fourth cumulant enters.

## 2. Transfer coefficients K(l) (measured, `lr_probe.py`)

Method: central differences of the exact Gaussian-closure dynamics around its own trajectory. K is the increase in the
final MSE per unit squared error injected into the post-activation state of layer l:
- K_mean: per (rms δm)²;
- K_diag: per (rms δvar)²;
- K_off: per ‖δC‖_F² for a random symmetric off-diagonal error;
- K_coh: per ‖δC‖_F² for an error shaped like C_off;
- K_top: per ‖δC‖_F² for an error along the top eigencomponent.

MLP 0 (MLPs 1, 2 agree within ±30 %, `results/lr_mlp*.json`):

| l | K_mean | K_diag | K_off | K_coh | K_top |
|---|---|---|---|---|---|
| 0 | 0.023 | 1.0e-4 | 2.6e-8 | 2.3e-8 | 3.4e-7 |
| 1 | 0.041 | 1.5e-3 | 7.4e-8 | 6.4e-8 | 4.8e-7 |
| 2 | 0.063 | 1.9e-4 | 8.6e-8 | 1.2e-7 | 2.5e-7 |
| 3 | 0.109 | 7.9e-4 | 1.8e-7 | 1.5e-7 | 2.0e-7 |
| 4 | 0.122 | 4.8e-4 | 2.3e-7 | 2.0e-7 | 2.4e-7 |
| 5 | 0.171 | 4.9e-4 | 2.7e-7 | 3.2e-7 | 3.9e-7 |
| 6 | 0.236 | 6.6e-4 | 4.4e-7 | 4.6e-7 | 2.1e-7 |
| 7 | 0.367 | 5.3e-4 | 5.4e-7 | 4.7e-7 | 2.2e-7 |
| 8 | 0.379 | 1.4e-3 | 6.0e-7 | 4.9e-7 | 1.4e-7 |
| 9 | 0.386 | 9.0e-4 | 6.2e-7 | 5.6e-7 | 6.4e-8 |
| 10 | 0.509 | 1.7e-3 | 7.0e-7 | 6.8e-7 | 4.3e-8 |
| 11 | 0.693 | 8.7e-4 | 7.6e-7 | 6.1e-7 | 5.1e-8 |
| 12 | 0.704 | 2.2e-3 | 6.3e-7 | 5.1e-7 | 2.3e-8 |
| 13 | 0.665 | 7.7e-4 | 4.7e-7 | 4.4e-7 | 1.3e-8 |
| 14 | 0.924 | 3.1e-4 | 2.9e-7 | 2.4e-7 | 9.1e-9 |
| 15 | 1 | — | — | — | — |

Readings:
- **Mean errors contract by ≈ 0.8 per layer in energy**, so early local mean errors are cheap.
- **Off-diagonal covariance errors cost 2.6e-8 per ‖δC‖_F² at layer 0**, rising to ≈ 7e-7 at layers 9–12 and falling at
  13–14.
- **K_coh ≈ K_off.** An error shaped like C_off does no more harm than a random one; the top eigendirection does less
  harm at depth.
- **K_diag is a single-draw estimate.** It couples to the trace and fluctuates 3× between neighbouring layers.

## 3. Where the Gaussian closure's error lives (MC atlas, N = 2 × 262,144 per MLP, half cross products)

The local errors of the exact Gaussian closure at every layer:
- ΔC_l = C(a_l)^{MC} − G(μ_l, C(z_l)^{MC}), the closure fed the true pre-activation mean and covariance;
- Δm_l and Δvar_l, the same for the per-neuron mean and variance.

Priced with K (MLP 0 per layer; the totals for MLPs 1–2 are below the table):

| l | ‖ΔC_off‖_F² | ‖ΔC‖/‖C_off‖ | price off | rms² Δvar | price diag | rms² Δm (− truth noise) | price mean |
|---|---|---|---|---|---|---|---|
| 1 | 0.71 | 5.3 % | 5.2e-8 | 2.3e-6 | 3.4e-9 | 2.2e-7 | 9.1e-9 |
| 3 | 0.82 | 5.9 % | 1.5e-7 | 2.1e-6 | 1.6e-9 | 4.4e-7 | 4.8e-8 |
| 5 | 0.78 | 6.1 % | 2.1e-7 | 1.5e-6 | 7.3e-10 | 3.4e-7 | 5.8e-8 |
| 7 | 0.74 | 6.5 % | 4.0e-7 | 9.9e-7 | 5.3e-10 | 3.0e-7 | 1.1e-7 |
| 9 | 0.69 | 6.0 % | 4.3e-7 | 7.4e-7 | 6.7e-10 | 3.1e-7 | 1.2e-7 |
| 11 | 0.50 | 5.7 % | 3.8e-7 | 4.7e-7 | 4.1e-10 | 2.4e-7 | 1.7e-7 |
| 13 | 0.51 | 5.1 % | 2.4e-7 | 4.6e-7 | 3.5e-10 | 2.2e-7 | 1.4e-7 |
| 15 | 0.42 | 4.8 % | — | 3.4e-7 | — | 1.7e-7 | 1.7e-7 |

Totals per MLP:

| MLP | price off | price diag | price mean | priced sum | measured Gaussian raw | off after D21 term | off after D21 + κ4(2,2) |
|---|---|---|---|---|---|---|---|
| 0 | 3.78e-6 | 1.3e-8 | 1.50e-6 | 5.29e-6 | 5.00e-6 | 2.25e-8 | 4.4e-9 |
| 1 | 4.49e-6 | 8.1e-9 | 2.45e-6 | 6.95e-6 | 4.75e-6 | 2.77e-8 | 6.0e-9 |
| 2 | 3.24e-6 | 9.1e-9 | 1.55e-6 | 4.80e-6 | 4.72e-6 | 2.16e-8 | 4.3e-9 |

Readings:
- **The books close to 2–6 % on MLPs 0 and 2.** The incoherent superposition over-predicts MLP 1 by 46 %, i.e. local
  errors cancel partly between layers there. This is EscAI's "deferred-mean" linearity, now resolved by channel and layer.
- **Anatomy.**
  - ≈ 70 % of the error is the non-Gaussian correction of the off-diagonal covariance. ΔC is 5–6.5 % of C_off at every
    layer ≥ 1, spread over layers 3–13.
  - ≈ 30 % is the per-neuron mean readout, mostly at layers ≥ 7.
  - < 0.3 % is the per-neuron variance (≈ 1e-8: small, but at the bar's scale).
- **What suffices locally.** Apply first-order bivariate Edgeworth with true slices and factorised Gaussian weights:
  - The D21 term ½ (E[δ(z_a) 1(z_b > 0)] D21_ab + transpose) carries 92–100 % of ‖ΔC‖² at layers 1–15.
  - The κ4 (2,2) term ¼ E f''_a E f''_b κ4(z_a, z_a, z_b, z_b) takes the covariance channel from ≈ 2.5e-8 to ≈ 5e-9.
  - The (3,1) term is worth 3e-10.
- **Structure of ΔC.** It is 98 % incoherent with respect to C_off (squared cosine 0.008–0.023).
- **The readout channel after first-order Edgeworth with true κ3 and κ4 diagonals is not resolved here.**
  - It sits at the truth-noise floor: subtracting var(a)/N per layer is uncertain by ±5e-8 per layer, as layer 0, where
    the truth is exactly Gaussian, shows.
  - EscAI's oracle puts the K3/K4 readout bias at 1.2e-10 at the final layer.

## 4. Necessary conditions for any estimator at ≤ 0.15 B with raw ≤ 1.5e-8

The final MSE ≈ Σ_l [K_off(l) ‖δC_l‖_F² + K_mean(l) rms²(δm_l) + K_diag(l) rms²(δvar_l)], verified above to within the
inter-layer cancellation.

**N1. A dense, essentially exact Gaussian covariance arrow at every layer.**
- **Tolerance.** If the whole budget goes to covariance errors, C(a_l) must be known to
  ‖δC_l‖_F ≲ (1.5e-8 / 16 / K_off(l))^{1/2}. That is ≈ 0.04 at layers 7–12, i.e. **≈ 0.3 % of ‖C_off‖_F**, and
  0.1–0.2 at layers 0–3.
- **Consequence.** Low-rank, sketched or banded representations of C must reach this. EscAI's measured covariance sketch
  (0.70 % relative at layer 14) is 2–3× too coarse.
- **Cost.** The nonlinear entrywise map C(z) → C(a) needs every entry of C(z). The arrow costs 2 units per layer
  (1.1 with Strassen L5), 18–32 units in total.

**N2. The non-Gaussian covariance correction must be carried at every layer from 1 to 14, to ≲ 6 % uniform.**
- **Price.** Omitting it costs 3.2–4.5e-6. The price is quadratic, ≈ 3.8e-6 ε² for a uniform relative error ε. So
  ε ≤ 6.3 % for the whole budget, and **ε ≤ 4.5 % for half of it**.
- **This is the "ε² law" (504aldo's 4.2e-6 ε²), derived.** Its coefficient is the priced energy of ΔC.
- **Weighting.** Weighted by K_off, layers 6–12 carry 60 % of the price and layers 0–2 only 3 %. So ε may be ≈ 3× larger
  at layers 0–2 and 14 than at 6–12.
- **End-to-end check.** FC's D21 error grows from 2.5 % (layer 2) to 10 % (layer 15) and it still scores 2.4e-8.

**N3. Within ΔC, the D21 term is necessary everywhere, the κ4 (2,2) slice is necessary (1.8e-8 if dropped), and the
(3,1) slice is not (3e-10).**

**N4. The per-neuron readout needs the diagonal κ3 and κ4 of z_l at layers ≥ ≈ 5, to ≲ 7 %.** The Gaussian readout
costs 1.5–2.5e-6. Layers 0–3 need almost nothing (K_mean ≤ 0.11).

**N5. The per-neuron variance correction (the κ3 term of var(a)) is marginally necessary** (≈ 1e-8 if dropped).

**N6. Old content is necessary at n = 1024, and dominates D21 at depth** (measured inside FC, §5, MLP 0, κ4 supplied):

| sources carried | raw (MLP 0) | ε(D21) at layer 14 |
|---|---|---|
| all ages | 3.5e-8 | 9.6 % |
| ages ≤ 4 | 8.5e-7 | 75 % |
| ages ≤ 2 | 1.9e-6 | 89 % |

Content older than 4 layers is ≈ 70 % of D21 at layers ≥ 8 (the old-content stream extrapolated 0.5–0.8 from widths
≤ 256; it does not shrink). Dropping it costs 25× the bar. By N2 it must be carried at every age to ≈ 7 % relative at
layers 6–12.

**N7. The joint fourth cumulant is needed only through n-vectors.** FC was run with the atlas's κ4 supplied in
different forms (MLP 0):

| κ4 of z_l supplied as | raw |
|---|---|
| nothing | 3.5e-7 |
| full (2,2), (3,1) slices + diagonal, birth slices only | 2.0e-7 |
| full slices + diagonal, everywhere (births, covariance, readout) | 3.5e-8 |
| **column means of the (2,2), (3,1) slices + diagonal** | **3.8e-8** |
| column means of (2,2) only + diagonal (no (3,1)) | 3.9e-8 |
| column means × 0.7 | 3.6e-8 |

Where κ4 matters most is not the covariance term of N3. It is the **(2,1) slice of the next post-activation third
cumulant**: κ3(a_i, a_i, a_k) gets a first-order correction ½ K22(z)_ik w2_k (Φ_i − m_i w2_i) + (3,1) terms. Its
transport to D21(z_{l+1}) is dominated by the coherent column sums Σ_i (W_ia² ≈ 2/n). Without κ4, FC's transported
birth slices are 20 % wrong; with it they match MC within noise (layers 1–2, `check_birth.py`). The coherent sum is
O(√n) larger than the incoherent part, which is why the n-vector of column means suffices.

**Not necessary (the "obvious" requirements that the measurements remove):**
- **Coherent or special directions in the covariance correction.** ΔC is 98 % incoherent; K_coh ≈ K_off;
  K_top ≤ K_off at depth. No privileged basis helps, and no cheap coherent part can stand in for the whole.
- **Second-order Edgeworth on the covariance at n = 1024.** First order with true slices leaves 4–6e-9 of 3.2–4.5e-6.
  The covariance-drift binding found at widths 64–128 (heisenberg R2–R3) is a small-width effect.
- **The joint fourth cumulant as n² slices.** Only its diagonal and the column means of the (2,2) slice are needed, to
  ≈ 30 % (N7).
- **Exact third-order content.** The means need it only through two objects:
  - the n × n matrix D21 (weighted by E[δ(z_a) 1(z_b > 0)]) to ≲ 5–10 % in Frobenius;
  - the n-vector diag κ3 to ≲ 7 % at late layers.

  All ages contribute to these contractions (N6), so old content is needed, but no more precisely than this.
- **Accuracy at early layers.** Layers 0–2 carry ≤ 3 % of the covariance price and K_mean ≤ 0.06.

## 5. Part 2: the first-order Wiener-chaos estimator (FC)

**Principle.** With independent Gaussian weights and Gaussian input, every activation is a functional of one Gaussian
field (the fresh-weight lemma, §1).
- **Creation.** At each ReLU, non-Gaussianity appears as second Wiener chaos ½ w2 :g²: (w2 = E relu'').
- **Transport.** The first-order term of relu(g + q) is relu'(g) q, whose second-chaos projection is Φ q. Old
  non-Gaussianity therefore travels linearly through the gates.
- **What the third cumulant then needs.** Only two-time objects:
  - the cross-covariance Y_s(l) = Cov(g_s, z_l) = S_s diag(Φ_s) Z_s(l). This is Stein's lemma, exact for Gaussian g;
  - the gated propagator Z_s(l) = W_{s+1} diag(Φ_{s+1}) ⋯ W_l;

  with

  D21(z_l)_ab = Σ_{s<l} Σ_r w2_{s,r} [ Y_ra² Z_rb + 2 Y_ra Z_ra Y_rb ].

  This is the Heisenberg/Duhamel identity written in the chaos basis.
- **Slices.** Repeated indices on the same neuron are O(1), not suppressed by correlations. So each layer adds a
  slice-supported source with atoms (Z_r, Z_r, (ΔZ)_r). Δ = the exact Gaussian per-pair and per-neuron slices of κ3(a)
  (1-D Gauss–Hermite with the conditional ReLU mean), plus their first-order Edgeworth corrections by the current D21,
  diagonal κ4 and K22 column means, minus what the chaos already carries.
- **Fourth order: the mean-field recursion** (`k4mf.py`, new). Per layer:
  - (i) per-neuron κ4(a) by univariate Edgeworth (κ3/6, κ4/24, κ3²/72 terms on truncated-normal derivatives);
  - (ii) column means of the post-activation (2,2) slice from the Gaussian Mehler terms (C, C², Cov²), the K22 and
    D21 corrections;
  - (iii) transport by W_{l+1} using (1/n) Σ_i W_ri W_r'i → (2/n) δ_rr' in coherent sums.

  The all-distinct κ4 is dropped. By the fresh-weight lemma it enters every index sum incoherently. Cost O(n²) per
  layer. Against the atlas it tracks the K22 column means to 6–35 % (ratio 0.70–1.00 by depth) and the diagonal κ4 to
  30–55 %, which is inside the N7 tolerance.
- **The rest.** The covariance arrow is the exact Gaussian kernel plus the first-order bivariate Edgeworth terms of §3
  (D21, K22, K31). The readout and variance are first-order Edgeworth with κ3 and κ4.

**Bench measurement, w1024_d16, all 6 MLPs** (`run_bench.py`, `results/bench_fc_mf_*.json`; raw = final MSE − bench
truth noise):

| MLP | 0 | 1 | 2 | 3 | 4 | 5 | mean ± s.e. |
|---|---|---|---|---|---|---|---|
| exact Gaussian closure | 5.00e-6 | 4.75e-6 | 4.71e-6 | 3.53e-6 | 3.51e-6 | 3.10e-6 | 4.10e-6 ± 0.33e-6 |
| **FC (self-contained)** | **3.24e-8** | **1.81e-8** | **3.03e-8** | **2.28e-8** | **3.89e-8** | **3.96e-8** | **3.03e-8 ± 0.35e-8** |

For reference: the public chain is ≈ 2.1e-8 at 0.25 B; the leaders 1.14–1.5e-8 at 0.11–0.15 B.

**Ablations and oracles (MLP 0 unless stated):**

| variant | raw |
|---|---|
| FC without κ4, second chaos only | 5.2e-7 |
| FC without κ4, exact slices | 3.5e-7 |
| κ4 diagonal-source model only (per-neuron Gaussian κ4 transported, `k4own`) | 2.7e-7 |
| κ4 from the atlas: full slices, used only in birth slices | 2.0e-7 |
| κ4 from the atlas: diagonal only, no (2,2) column means | 2.3e-7 |
| κ4 from the atlas: column means × 0.3 | 1.1e-7 |
| κ4 from the atlas: column means × 0.7 | 3.6e-8 |
| κ4 from the atlas: column means (MLPs 0, 1, 2) | 3.8e-8, 7.6e-9, 2.6e-8 |
| **κ4 from the mean-field recursion (self-contained)** | **3.2e-8** |
| old content: only sources of age ≤ 4 / ≤ 2 | 8.5e-7 / 1.9e-6 |
| per-source row basis k = 256 for ages > 2 | 1.4e-7 |
| shared Oseledets subspace q = 256 / 512 for ages > 2 | 2.3e-7 / 3.6e-8 |
| atom pruning, Wick and slice atoms to 50 % / 25 % (MLPs 0–2) | 1.3e-6 / 3.4e-6 |
| Wick atoms to 25 % only / slice atoms to 50 % only | 3.1e-6 / 8.2e-8 |

**D21 error of FC against the atlas.** The profile is the same on every MLP: 2.4–2.9 % at layer 2, 4.6–4.9 % at layer
5, 6.2–6.5 % at layer 8, 8–8.6 % at layer 11 and 9.6–10 % at layer 14. The growth is the compounding of the omitted
second-order (hyperedge) birth diagrams. oracle1024's leg-partition closure takes the one-step error from ≈ 4 % to
≈ 1 %. With them, N2's ε² law predicts the covariance channel at ≲ 5e-9, i.e. FC near 1e-8.

**Cost at n = 1024 (units; prototype in numpy, 48 s per MLP on 2 threads).**
- There are 120 (source, target) pairs. Each needs:
  - transport of the Y, Z and slice-leg rows: 3 products;
  - the two Wick contractions (Y∘Y)ᵀ diag(w2) Z and (Y∘Z)ᵀ diag(w2) Y: 2 products;
  - the two slice contractions: 2 products.
- That is ≈ 840 units, plus a 30-unit covariance arrow and O(n²) everything else: **≈ 0.85 B dense, ≈ 0.5 B with
  Strassen L5**.
- **Adjusted ≈ 1.5–2.5e-8**: 16× better than the Gaussian closure's 4.1e-7 at the floor, but 10–15× above the bar.
- The 40-node birth-slice quadrature is negligible in FLOPs but ≈ 800 calls per layer. A port needs a closed form or
  ≤ 8 nodes for the residual cap.

## 6. Verdict and the precise open problem

- **What the region says.** The leaders' accuracy (≈ 1–1.5e-8) is reachable with the information listed in §4, and FC
  shows that a principled first-order representation of exactly that information gets within 2× of it (3.0e-8 on the
  6 bench MLPs), with a known next lever, the second-order birth diagrams.
- **Accuracy is not the binding problem at n = 1024. Cost is.** The binding object is the all-age D21 contraction
  D21(t) = Σ_{s<t} Σ_r w2_{s,r} [Y_ra² Z_rb + 2 Y_ra Z_ra Y_rb] (+ slice analogue), needed to ≈ 5–10 % in Frobenius at
  layers 6–12. The measurements establish three properties of it at n = 1024:
  - (i) content of every age is needed (N6);
  - (ii) the legs do not compress into subspaces below n/2;
  - (iii) every atom of every source carries comparable energy, so atoms cannot be pruned.
- **A cost bound under those properties.** For a bilinear realisation (matrix products over the atom index), the
  inner dimension at target t is t·n atoms, so the contractions alone cost ≥ 2 t n³ and the transports ≥ 2 t n³.
  That is ≥ 4 × 120 = 480 units (≥ 270 with Strassen), against the leaders' 113–154 units in total.
- **Consequence for the leaders.** They cannot be computing this object atom by atom. Either:
  - (a) they carry an *aggregated* state whose D21 readout is O(n³) per layer, i.e. a merged representation of the old
    third cumulant with ≲ n atoms re-fitted each layer (CP merging; the old-content stream measured R = n merged atoms
    at 0.7–1.4 % at widths 64–128, and at n = 1024 our tolerance is the looser ≈ 7 %); or
  - (b) their error budget is spent differently: e.g. cruder D21 at layers where K_off is small (0–3 and 14 carry
    ≲ 15 % of the price), compensated by inter-layer cancellation of the kind MLP 1 shows.
- **The concrete next experiment.** Inside FC, replace the old sources (age > 2) by one CP-merged atom set of size
  R ∈ {n/2, n} per layer (ALS in the q-space of the young legs). Measure raw on the 6 bench MLPs and price the merge.
  If R = n gives ≤ 5e-8, an O(L n³) estimator at ≈ 0.15 B with raw ≈ 3e-8 exists, and the hyperedge diagrams are the
  lever to the bar.

## 7. What this stream contributes to the fresh-slate programme

- **A fresh-weight lemma**, with its quantitative consequence, measured: incoherent errors are priced by Frobenius norm
  with layer coefficients K(l). The Gaussian closure's whole error is accounted for channel by channel.
- **The first derivation of the "ε² law" from first principles.** Its coefficient is the priced energy of the
  non-Gaussian covariance correction.
- **The identification of what the joint fourth cumulant contributes:** two n-vectors per layer, closed by a
  memoryless mean-field recursion, again a fresh-weight consequence.
- **A self-contained, principled estimator at 3.0e-8 at n = 1024.** It is the strongest fresh design so far: the
  convergence dossier's best projection was ≈ 4–7e-7. Its cost gap is reduced to one explicitly stated object.
