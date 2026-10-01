# Breakthrough stream "region": the feasible region at n = 1024, and a design inside it

*Status: v0 (1 Oct 2026, ≈ 20:15 UTC). Part 1 (necessary conditions) has measured numbers on one bench MLP (w1024_d16 #0).
MLPs 1 and 2 are being added. Part 2 has a prototype under test (§5). All numbers are at n = 1024, depth 16, on the shared
bench networks, unless stated.*

Files: `gclose.py` (exact bivariate-Gaussian closure), `lr_probe.py` (linear-response transfer coefficients),
`mc1024.py` + `mc_analyse.py` (Monte Carlo atlas at n = 1024 in two independent halves, noise-corrected norms),
`fc.py` + `diag_fc.py` (the Part 2 prototype and its layer-by-layer diagnosis against the atlas). Results: `results/*.json`.
Large atlases (`data/`, 0.7 GB each) are not committed; they regenerate with `python mc1024.py MLP 262144`.

## 0. The bar in this stream's units

- Raw ≤ 1.5e-8 means the final-layer per-neuron error is ≤ **1.22e-4 rms**. Final-layer activations have standard
  deviation ≈ 0.25–0.30 (avg_variance 0.056–0.087 on the bench), so this is ≈ 4e-4 of a standard deviation.
- Cost ≤ 0.11–0.15 B = 113–154 units, i.e. **7–10 units per layer** (1 unit = one dense 1024³ product).
- Residual time ≤ 0.4 s ≈ 6–8k flopscope calls with margin, i.e. **≈ 400–500 calls per layer**.
- Calibration: the exact bivariate-Gaussian covariance closure (`gclose.py`, no linearisation) scores raw
  **4.10e-6 ± 0.33e-6** on the 6 bench MLPs, and the bench's linearised baseline 4.30e-6. So the linearised
  cross-covariance is not where the Gaussian closure loses; the bar is ≈ 270× below the exact Gaussian closure.

## 1. A structural lemma: fresh weights make Frobenius norms the currency

W_{l+1} is independent of the law of a_l and of anything an estimator computes from W_1..W_l. Hence any error δ that an
estimator makes in a layer-l object enters layer l+1 as a Gaussian chaos in the fresh columns w_c of W_{l+1}, whose second
moments depend only on a few rotation invariants of δ (jointly with the true state). For a symmetric error δC in C(a_l):

  δσ²_c = w_cᵀ δC w_c,   E = (2/n) tr δC,   Var = 2 (2/n)² ‖δC‖_F²,

and the off-diagonal image W_{l+1}ᵀ δC W_{l+1} is Wigner-like with entry variance (2/n)² ‖δC‖_F². For an error δκ in a
third cumulant: δκ3(z_c) = ⟨δκ, w_c^{⊗3}⟩ has mean 0 and variance ≈ 6 (2/n)³ ‖δκ‖_F², and ‖δD21(z_{l+1})‖_F² ≈ (16/n) ‖δκ‖_F².
**Consequence:** an error that is not aligned with the true state ("incoherent") costs, to first order, a price
proportional to its squared Frobenius norm, with a per-layer coefficient that does not depend on its shape. Errors that
are aligned (a scale error of C, an error along the top eigenvector) can be priced differently, and we measure both.

## 2. The transfer coefficients K(l) (measured, `lr_probe.py`, MLP 0)

Central differences of the exact Gaussian-closure dynamics around its own trajectory (the transfer is a property of the
propagation; the closure is only the linearisation point). K = increase of the final MSE per unit squared error injected
into the post-activation state of layer l: per (rms δm)² for the means, per (rms δvar)² for the diagonal of C(a_l), per
‖δC‖_F² for off-diagonal errors (random symmetric = "off", the shape of C_off itself = "coh", the top eigencomponent
of C(a_l) = "top").

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

Readings. (i) Mean errors contract by ≈ 0.8 per layer in energy (0.023 at l = 0): early local mean errors are cheap.
(ii) Off-diagonal covariance errors are priced at 2.6e-8 per ‖δC‖_F² at layer 0, rising to ≈ 7e-7 at layers 9–12 and
falling again at 13–14 (the last covariance only reaches the means through one variance). (iii) K_coh ≈ K_off: a
C_off-shaped error is no more harmful than a random one of equal Frobenius norm; the top eigendirection is *less*
harmful at depth. The lemma's "Frobenius is the currency" holds at the level of ±30 %. (iv) K_diag fluctuates 3× between
neighbouring layers (single random draw; the diagonal couples to the trace, i.e. to the mean direction): treat it as
≈ 1e-3 per (rms δvar)².

## 3. Where the Gaussian closure's error lives (MC atlas, N = 2 × 262,144, noise-corrected by half cross products)

At every layer, ΔC_l = C(a_l)^{MC} − G(μ_l, C(z_l)^{MC}) is the *local* error of the exact Gaussian closure fed the true
pre-activation mean and covariance; Δm_l and Δvar_l the same for the per-neuron mean and variance. Priced with K:

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
| **sum over 16 layers** | | | **3.78e-6** | | **1.3e-8** | | **1.50e-6** |

**The books close:** the priced local errors sum to 5.29e-6; the exact Gaussian closure's measured raw MSE on this MLP is
5.00e-6. So, at n = 1024, the final MSE is the incoherent sum over layers of local errors times K(l), to ≈ 6 %
(this is EscAI's "deferred-mean" linearity result, now resolved by channel and layer).

Anatomy of the 4–5e-6: **71 % is the non-Gaussian correction of the off-diagonal post-activation covariance**
(ΔC ≈ 5–6.5 % of C_off at every layer ≥ 1, spread over layers 3–13), **28 % is the per-neuron mean readout**
(mostly layers ≥ 7, K_mean ≥ 0.37), **0.3 % is the per-neuron variance** (1.3e-8 — small, but at the bar's scale).

What suffices locally (true slices, first-order bivariate Edgeworth, factorised Gaussian weights; `mc_analyse.py`):

| local correction applied to ΔC | priced residual of the covariance channel |
|---|---|
| none (Gaussian closure) | 3.78e-6 |
| D21 term ½(E[δ(z_a)1(z_b>0)] D21_ab + transpose) | 2.25e-8 |
| + κ4 (3,1) term | 2.21e-8 |
| + κ4 (2,2) term ¼ E f''_a E f''_b κ4(z_a,z_a,z_b,z_b) | **4.4e-9** |

and for the readout, first-order Edgeworth with the true κ3, κ4 diagonals leaves ≤ 1.6e-7 summed over layers (this
number is at the noise floor of the per-neuron mean truth, ±1e-8 per layer; EscAI's oracle puts the K3/K4 readout bias at
1.2e-10). The D21 term alone carries 92–100 % of ‖ΔC‖² at layers 1–15 (68 % at layer 0, where ΔC is tiny). ΔC is
98 % incoherent with respect to C_off (squared cosine 0.008–0.023).

## 4. Necessary conditions for any estimator at ≤ 0.15 B with raw ≤ 1.5e-8 (v0)

Write the final MSE as Σ_l [K_off(l) ‖δC_l‖_F² + K_mean(l) rms²(δm_l) + K_diag(l) rms²(δvar_l)] (verified to 6 % above).

**N1. A dense, essentially exact Gaussian covariance arrow at every layer.** C(a_l) must be known to
‖δC_l‖_F ≲ (1.5e-8 / 16 / K_off(l))^{1/2} if covariance errors get the whole budget: ≈ 0.04 at layers 7–12, i.e.
**≈ 0.25–0.3 % of ‖C_off‖_F (13–16)**, and ≈ 0.1–0.2 at layers 0–3. Any low-rank, sketched or banded representation of
C must reach this; EscAI's measured covariance sketch (0.70 % relative at layer 14) is 2–3× too coarse. The arrow
W_{l+1}ᵀ C W_{l+1} costs 2 units per layer dense (1.1 with Strassen L5): 18–32 units in total. Not avoidable: the
nonlinear entrywise map C(z) → C(a) needs every entry of C(z).

**N2. The non-Gaussian covariance correction must be carried at every layer from 1 to 14 to ≲ 6 % (uniform).**
Its full omission costs 3.8e-6 and its price is quadratic: ≈ 3.8e-6 ε² for a uniform relative error ε, i.e.
ε ≤ 6.3 % for the whole budget, **ε ≤ 4.5 % for half of it**. This *is* the "ε² law" (504aldo's 4.2e-6 ε²), now derived:
its coefficient is the priced energy of ΔC. Weighted by K_off, layers 6–12 carry 60 % of the price; layers 0–2 only 3 %,
so ε may be ≈ 3× larger at layers 0–2 and 14 than at 6–12.

**N3. Within ΔC, the D21 term is necessary everywhere and the κ4 (2,2) slice is necessary too; the (3,1) slice is
not.** Dropping the (2,2) term costs 1.8e-8 (more than the whole bar), dropping (3,1) costs 3e-10.

**N4. The per-neuron readout needs the diagonal κ3 and κ4 of z_l at layers ≥ ≈ 5 to ≲ 7 %.** Gaussian readout costs
1.5e-6; first-order Edgeworth with true diagonals is at the noise floor. Its price per layer is K_mean(l) × rms²(Δm), so
layers 0–3 need almost nothing.

**N5. The per-neuron variance correction (κ3 term of var(a)) is marginally necessary** (1.3e-8 if dropped).

**Not necessary (the "obvious" requirements that the measurements remove):**

- *Coherent or special directions.* Nothing about ΔC is special: it is 98 % incoherent, and K_coh ≈ K_off, K_top ≤ K_off
  at depth. A representation needs no privileged basis for the correction; conversely, no "cheap coherent part" exists
  to carry instead of the whole.
- *Second-order Edgeworth on the covariance at n = 1024.* First order with true slices leaves 4.4e-9 of 3.8e-6
  (a factor 860). The covariance-drift binding found at widths 64–128 (heisenberg R2–R3) is a small-width effect;
  at n = 1024 the κ4 slices enter only through the (2,2) term, ≈ 2 % of ‖ΔC‖² at layers 1–2 and ≤ 0.1 % at depth.
- *Exact third-order content.* What the means need from κ3 is (a) the n × n matrix ΔC_l (equivalently D21 weighted by
  E[δ(z_a)1(z_b>0)]) to ≲ 5 % in Frobenius, and (b) the n-vector diag κ3(z_l) to ≲ 7 % at late layers. Old content
  is necessary only through these contractions, and only at the K-weighted accuracy above (§6 measures how much of it
  that is).
- *Accuracy at early layers.* Layers 0–2 carry ≤ 3 % of the covariance price and K_mean ≤ 0.06: a crude closure there
  costs nothing measurable.

## 5. Part 2 prototype: the first-order Wiener-chaos chain (FC), in progress

**Principle.** With independent Gaussian weights and Gaussian input, every activation is a functional of one Gaussian
field; non-Gaussianity is created at each ReLU as second Wiener chaos (½ w2 :g²:, w2 = E relu'') and carried to later
layers linearly through the gates (the first-order term of relu(g + q) is relu'(g) q, whose second-chaos projection is
Φ q). The joint third cumulant then needs only **two-time objects**: the cross-covariance Y_s(l) = Cov(g_s, z_l) =
S_s diag(Φ_s) Z_s(l) (Stein's lemma, exact for Gaussian g) and the gated propagator Z_s(l), with

  D21(z_l)_ab = Σ_{s<l} Σ_r w2_{s,r} [ Y_ra² Z_rb + 2 Y_ra Z_ra Y_rb ],

plus a slice-supported source per layer (exact per-neuron and per-pair Edgeworth slices minus what the chaos already
carries). The covariance arrow is the exact Gaussian kernel + the D21 term of §3; the readout is first-order Edgeworth.
This is the Heisenberg/Duhamel identity written in the chaos basis, and the fresh-weight lemma is what makes the
cross-covariance the only object the third cumulant needs. Cost if done densely: 4 products per (source, target) pair →
≈ 480 units (0.47 B); §6 is about cutting this.

**First measurements (MLP 0):**

| variant | raw final MSE | D21 error ε vs MC (layers 2–15) |
|---|---|---|
| exact Gaussian closure | 5.00e-6 | — |
| FC, second chaos only | 5.16e-7 | 9 % (l = 1), 16–24 % |
| FC, + exact Gaussian slices of each birth | 3.93e-7 | — |
| FC, + first-order non-Gaussian slices | 3.51e-7 | 2.4 % (l = 1), 14–25 % |

10–14× better than the Gaussian closure. The chain's D21 is right at layer 1 (2.4 %) and 14–25 % off from layer 2 on.
The oracle1024 stream measured ≈ 4 % for the same class of model (leading Wick + exact slices + gated old transport)
given true inputs, so the 14–25 % is a defect of this prototype, being located now with an atlas that also records the
post-activation slices (`MC_A=1`).

## 6. Next (in progress)

1. Locate the FC D21 defect (post-activation slice atlas); target ε ≈ 4 % → predicted raw ≈ 3.8e-6 × 0.04² + readout
   ≈ 1e-8.
2. Price the old content at n = 1024 inside FC: windowed variants (ages ≤ w only), D21 ε and end-to-end MSE.
3. Cost reduction: per-source row bases from the Oseledets subspaces of the propagators (rank k per old source →
   ≈ 6k/n units per pair instead of 4); measure the k needed for the N2 tolerance.
4. Replicate §§2–3 on MLPs 1–2 (atlases done or running).
