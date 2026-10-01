# Tropical / temperature design stream — DESIGN v1 (v0 + verdict, §6)

*1 Oct 2026. Fresh-slate stream "tropical". Status: v0 (derivation, realisations, predictions, falsification tests). Results go in RESULTS.md as they land.*

## 0. Summary of v0

The principle is stated in §1. It yields two exact identities. These make a clean "zero-temperature skeleton → Gaussian lift" estimator possible, which I call TCT (tropical-curvature transport, §3). However, v0 also argues, and a first measurement supports, that **at He initialisation with n = 1024, L = 16 the network sits deep in the high-temperature phase on both tropical axes**:

- (a) the *gate axis*: a quarter to all of the gates of every layer are thermal (|μ/s| < 1), at every depth;
- (b) the *path axis*: the signed path sum is a directed polymer in weak disorder, β_eff ≈ 1.1 against β_c = √(2 ln n) ≈ 3.7.

So the zero-temperature skeleton carries no information that survives the lift. The lift is the whole job, and its natural small parameter is n^{-1/2}, not the temperature. Every realisation I can derive has the Gaussian closure as its zeroth order (projected raw ≈ 4e-5 at 1024). The leading corrections are wall–gate correlations, which are joint, cross-layer and high-rank: the same bottleneck the brief documents.

The prototype and Stage Q runs test this prediction. They ask whether any tropical-exact ingredient gives a measurable gain over the zeroth order at fixed cost.

## 1. The principle, stated precisely

A bias-free ReLU MLP f = a_L is a **tropical rational map**. Split W = W⁺ − W⁻. By induction every neuron is z = P − Q, with P and Q convex, positively homogeneous and piecewise linear: support functions h_A, h_B of polytopes (Newton polytopes). The recursion is:

- a linear layer: Minkowski sums, (P', Q') = (Σ W⁺P + Σ W⁻Q, Σ W⁺Q + Σ W⁻P);
- a ReLU: a convex hull, relu(P − Q) = max(P, Q) − Q.

The linear regions are the cones of a polyhedral fan, the common refinement of the normal fans. There are no biases, so every cell is a cone through 0.

**Temperature.** Replace max by T·log Σ exp(·/T). Then relu → softplus_T, and each tropical polynomial becomes a subtraction-free (positive) exponential sum Z_T(x) = Σ_v exp(⟨v, x⟩ / T) over the vertices of its Newton polytope:

- a Minkowski sum becomes a product of partition functions;
- a convex hull becomes a sum.

f_T → f as T → 0 ("zero temperature, where signs die"). The Gaussian mean of one neuron is a difference of quenched free energies of two Gaussian-disorder models. At T = 0 it is a difference of expected suprema, i.e. of Gaussian mean widths: E a = w_G(conv(A ∪ B)) − w_G(B).

### Exact identities used (all verified at toy scale, §5 E0)

- **(I1) Homogeneity / Euler.** f(x) = x · ∇f(x), and ∇f is piecewise constant on the cones.
- **(I2) Tropical curvature (Stein + Euler).** E f = E[x · ∇f] = E[Δf]. Here Δf is the distributional Laplacian of a piecewise-linear function: a signed measure on the codimension-1 skeleton of the fan (the tropical hypersurface), with density on each wall equal to the jump of the normal derivative. **The Gaussian mean of a ReLU network equals the Gaussian mass of its tropical hypersurface, weighted by tropical multiplicities.** For a convex part (P or Q) the measure is positive.

  Per layer, with the coarea formula: E[a_{l,j}] = Σ_{k ≤ l} Σ_m E[ δ(z_{k,m}) |∇_x z_{k,m}|² · S_{(l,j) ← (k,m)} ]. Here S is the downstream gated propagator ∂a_{l,j}/∂a_{k,m} (S = 1 for (k,m) = (l,j)).

  Split at the last layer, this is an exact two-term recursion:

  E[a_{l,j}] = β_{l,j} + E[ g_{l,j} · (Δa_{l−1} W_l)_j ]

  - β_{l,j} = E[δ(z_{l,j}) |∇z_{l,j}|²] is the **own-wall birth**;
  - g = 1(z > 0) is the gate;
  - the second term is the inherited curvature seen through the open gate.

  Check at depth 1: E relu(w·x) = |w|² · φ(0) / |w| = |w| / √(2π). ✓
- **(I3) Cone measures.** Gaussian integrals of relu-products of jointly Gaussian variables are Gaussian measures of polyhedral cones:
  - 2 variables: the arc-cosine kernel;
  - 3 variables: spherical-triangle areas (Van Oosterom–Strackee);
  - and so on.

  These are the "exact tropical lift" of a Gaussian field through one layer of walls.
- **(I4) Temperature expansion at T → 0.** softplus_T(z) − relu(z) is even, with ∫ = (π²/6)T². So for any density p continuous at 0: E softplus_T(z) = E relu(z) + (π²/6) T² p(0) + O(T⁴). For the network, F(T) = F(0) + c₂T² + c₄T⁴ + …, with c₂ = (π²/6) × (the same Gaussian mass of the tropical hypersurface as in I2, density-weighted).

## 2. Where can a small parameter come from? (the phase diagram)

A "tropical skeleton → finite-temperature lift" estimator is good when the lift is perturbative. There are three candidate small parameters.

**(P1) Gate temperature τ_{l,j} = s_{l,j}/|μ_{l,j}|** (s, μ are the std and mean of the pre-activation). τ → 0 means the gate is frozen, and the network on its typical inputs is a single cone, i.e. linear. Only then does the zero-temperature skeleton (the dominant cone, the sign pattern of the mean) carry the answer, with an expansion in the number of thermal gates.

Measured (width 256, depth 16, 40k samples, `probe_freeze.py`):

| layer | 1 | 2 | 4 | 8 | 12 | 16 |
|---|---|---|---|---|---|---|
| fraction hot (\|μ/s\| < 1) | 1.00 | 0.86 | 0.58 | 0.30 | 0.26 | 0.24 |
| fraction frozen (\|μ/s\| > 3) | 0.00 | 0.00 | 0.01 | 0.12 | 0.27 | 0.27 |

The mean's share of the second moment, ρ_l, is 0.88 at depth 16; the infinite-width correlation map predicts 0.93. Freezing saturates near 25 % hot. **There is no zero-temperature limit with depth at L = 16.**

**(P2) Path disorder.** E[a_L] = Σ_paths Π W · E[x_{i₀} Π gates] is a signed directed polymer on a complete layered graph (n choices per step, L steps):

- Energy per step: ln|W| has std ≈ 1.1 (the log of a half-normal).
- The single-path (tropical, Viterbi) phase needs β ≳ √(2 ln n) ≈ 3.7.

So the network is in weak disorder: no path dominates, and the sum self-averages (annealed = quenched at leading order, fluctuations O(n^{-1/2})). The max-plus path skeleton is irrelevant.

**(P3) Width, ε = n^{-1/2}.** Each upstream wall enters z_{l,j} with amplitude |W_ij| ~ n^{-1/2}, so the *upstream* fan is microscopic and thermalises (CLT) into a smooth Gaussian-like field. Only the neuron's **own** wall is macroscopic. This is the one small parameter the tropical picture actually provides:

- **zeroth order:** a Gaussian field plus the one macroscopic own wall per neuron, which is the Gaussian closure / arc-cos lift;
- **first order:** one macroscopic upstream wall (wall–gate correlations, O(n^{-1/2}) per neuron);
- and so on.

The per-neuron error target is rms 1.2e-4 ≈ 0.004 · n^{-1/2} at n = 1024. So the first order must be right to about 2 %, and the second to about 10–20 %.

**(P4) Analytic continuation in T (I4).** This would need F(T) at some T > 0 that is cheaper than F(0). Nothing is cheaper at finite T:

- MC variance is T-independent;
- the high-T series (expansion about the linear network) has Gaussian-moment coefficients that grow like (2k)!!, so it is asymptotic only;
- the continuation from large T to T = 0 crosses the singularities of softplus at z = ±iπT, which sit on the integration path.

**Dead as an estimator.** It survives only as I4, an exact statement that F(T) − F(0) ∝ T² × (hypersurface mass): a diagnostic, not a lift.

**(P5) Newton-polytope (P − Q) evaluation.** E a = w_G(conv(A ∪ B)) − w_G(B). The two terms are each O(√n) × the answer, so the evaluation suffers a √n cancellation. The mean width of a convex hull has no closed form; Kubota/Crofton reduce it to random 1-D projections, i.e. Monte Carlo. **Dead.**

**(P6) Linear-to-ReLU homotopy.** Take relu_γ(z) = (1 − γ)z + γ|z|, with γ = ½ for ReLU. Each order in γ adds one macroscopic wall and is exactly computable by cone measures (I3) on the *Gaussian* linear network. The radius of convergence, however, is set by s(γ)² = (1−γ)² + γ², which vanishes at γ = (1 ± i)/2. The geometric ratio at γ = ½ is 0.71/order, so about 27 orders are needed for 1e-4. Resumming the scale (normalising by s(γ) at each γ) makes the zeroth order the Gaussian closure again: this is P3. **Dead as a series; it collapses into P3.**

## 3. The estimator that follows (TCT: tropical-curvature transport)

These are the identity I2 organised by P3 (wall amplitude).

**State at layer l.**
- Per-neuron mean m_l = E a_l (n).
- The thermal field of the pre-activations: mean μ_l and covariance C_l (n × n). This is the CLT limit of the microscopic upstream walls (P3), which is where the field comes from, not an assumption of the method.
- The wall weights β_l (n) and gates Φ_l = P(z_l > 0) (n).

**Per-layer operation (zeroth order).**
1. μ_l = m_{l−1} W_l and C_l = W_lᵀ K_{l−1} W_l. Here K_{l−1} is the second-moment kernel of a_{l−1}, given by the exact cone lift (I3, the arc-cos formula) of the field at l − 1.
2. Births: β_{l,j} = p_{z_{l,j}}(0) · E[|∇z_{l,j}|² | z_{l,j} = 0] ≈ s_{l,j} φ(μ_{l,j}/s_{l,j}). The ≈ uses the Gaussian density at the wall and |∇z|² ≈ s², the Euler/Stein consistency for a linear-field neuron.
3. Transport: E[g ⊙ (Δa_{l−1} W)] ≈ Φ_l ⊙ (m_{l−1} W_l), i.e. gates independent of upstream walls.
4. m_l = β_l + Φ_l ⊙ μ_l. This is algebraically the Gaussian closure sφ(μ/s) + μΦ(μ/s) — **the zeroth order of the tropical expansion is the Gaussian closure, derived rather than assumed.**

**Readout.** m_L.

**First-order lift (the tropical correction).** There are two factorisation errors in I2, both O(n^{-1/2}) per neuron:
- (i) **Birth:** the density at the wall is non-Gaussian, and the gradient norm is correlated with being on the wall;
- (ii) **Transport:** the gate g_{l,j} is correlated with the upstream wall (k, m). This is the cavity effect: on wall (k, m), a_{k,m} = 0 exactly, so downstream sees the cavity network without neuron (k, m), and the conditioning on z_{k,m} = 0 shifts the other neurons of layer k by their covariance.

At one upstream layer (k = l − 1), the transport correction for a jointly Gaussian field is a sum of bivariate cone integrals: (I3) with 2 walls (the target gate and the upstream wall). Only these terms are computable without carrying non-Gaussian joint state. Terms with k < l − 1 need the cross-layer joint law: the "old content" of the brief.

## 4. Cost at n = 1024 (units; costmodel prices)

- Zeroth order (TCT-0 = the Gaussian closure): per layer, one symmetric sandwich Wᵀ K W at 1.03 u (sym3, Strassen L5) plus O(n²) elementwise. Total ≈ 16–19 u ≈ 0.016–0.019 B. The multiplier sits at the floor 0.1.
- First order at one layer of depth (TCT-1): the bivariate-cone transport correction per (target j, upstream m) pair is an n × n elementwise kernel (erf/arccos at about 16–48 FLOPs per element) contracted with W:
  - a few Hadamard-then-contract products, ≈ 1 u each, plus 48n² FLOPs per element type (≈ 0.02 u);
  - per layer ≈ 3–5 u; total ≈ 50–80 u ≈ 0.05–0.08 B, still at the floor.
- Deeper (k < l − 1) corrections: transport of an n × n wall–field cross object per age. That is O(age) n³ per layer without compression, ≈ L²/2 · 1 u ≈ 128 u. A compressed basis is needed for it.

## 5. Error mechanism, predicted scaling, falsification

**Predictions.**
- TCT-0 raw MSE ∝ n^{-1} × const at fixed depth. This is the brief's Gaussian-closure 4e-5 at 1024, so ≈ 1.6e-4 at 256 and ≈ 6e-4 at 64.
- TCT-1 removes the k = l − 1 transport part of the O(n^{-1/2}) error. Using the brief's "~40 % of (2,1) structure is older than one layer", the best case is that TCT-1 removes ~60 % of the first-order error amplitude, i.e. ~85 % of the MSE. That gives raw ≈ 6e-6 at 1024: still 400× above the bar. A competitive design needs the old-content term (k < l − 1) to ~2 %, which is not a tropical object.

**Cheapest falsification tests.**
- E0 (exactness, toy): verify I2 at width 8–16 / depth 2–4 by direct high-precision integration.
- E1 (phase diagram): hot fraction and ρ by layer at widths 64–1024. Done at 256: §2 P1.
- E2 (where the zeroth-order error lives): at width 64–128, decompose truth − TCT-0 per layer into birth and transport errors via I2 (MC with exact gradients). If the transport error is dominated by k = l − 1 walls, TCT-1 is worth building. If it is spread over ages, TCT-1 caps far from the bar.
- E3 (Stage Q): TCT-0 and TCT-1 at widths 64/128/256, depth 16, ≥ 4 MLPs, raw MSE minus truth noise, width-scaling fit, projection to 1024.

## 6. Verdict (v1, 1 Oct evening)

**Measured.** All in RESULTS.md.

1. No zero-temperature small parameter exists at He init, L = 16: a quarter or more of the gates are thermal at every depth (R-P1).
2. The Newton-polytope (P − Q) form has a 10²³-fold sign cancellation at n = 1024, growing ×√(4n/π) per layer (R-P5).
3. The tropical-curvature identity I2 is exact (R-E0), but its natural split is non-perturbative: births and inherited transport each miss by O(1) and cancel (corr −1.00). The multiplicity E|∇z|² reaches 17 s² at depth 16 (R-E2).
4. TCT-0, the only workable realisation, equals the Gaussian closure with an exact two-wall lift: raw 4.5e-4 (w64) and 2.8e-4 (w128) on the bench, projected raw ≈ 5–7e-5 at 1024, 19 u (R-Stage Q).
5. Temperature as a deformation parameter makes the closure 400× more accurate at T ≈ 0.8, but the extrapolation back to T = 0 buys at most 1.6–2.5×, and only optimistically (R-E3).

**Measured at n = 1024** (bench w1024_d16, 6 MLPs):
- TCT-0 is raw 4.10e-6 ± 3.3e-7 at 19 u, i.e. adjusted 4.1e-7. My earlier 64–128 projection of 7e-5 was 17× pessimistic, because small-width scaling is pre-asymptotic.
- Adding the optimistic T-extrapolation gain (1.6–2.5× at 64–128; untested at 1024) gives at best raw ≈ 2–4e-6 and adjusted ≈ 2–4e-7. The cost is several closure passes, still at the 0.1 floor.
- The temperature bridge with a sampled residual (E4) does no better than plain sampling.

**Best tropical estimator: adjusted ≈ 4e-7 (range 2e-7 to 5e-7), about 250× above the bar of 1.6e-9. Not competitive.**

**Charged to what.** Following the brief's rule 3, the failure is charged to the realisation, not to the theory. The tropical dictionary (fan, Newton polytopes, multiplicities, max-plus paths) describes the network exactly, but in *uncentred* coordinates:

- the mean spike (ρ ≈ 0.9 of the second moment at depth) sits inside every tropical quantity;
- the answer is a small difference of large tropical terms: walls vs transport, P vs Q, signed paths;
- the network is in the high-temperature phase on both tropical axes (gates and paths).

The one small parameter the picture yields is n^{-1/2}: upstream walls are microscopic and thermalise, and the neuron's own wall is the only macroscopic one. That collapses onto moment/diagram closures, which carry the joint, cross-layer, high-rank structure that the tropical picture does not compress.

**What would make it competitive.** One of:

- (i) a regime with frozen gates (hot fraction → 0 at the layers that dominate the error), e.g. much deeper or trained networks — not this task;
- (ii) an independent, cheap estimate of the per-neuron wall densities p_{z}(0) and wall-conditional gradients. That is precisely the per-neuron law problem; given it, I2 and I4 become exact readouts;
- (iii) a *centred* tropical calculus, i.e. a fan/Newton-polytope structure for a − (r/E r)·m, in which the cancellations are removed analytically. I know of no such construction.

**Deciding experiment** (to reopen the route): measure the hot fraction and the I2 birth/transport error at n = 1024 on the layers where the closure error is born. If either shrinks with n, faster than the closure error itself, a tropical lift becomes perturbative there. Everything at 64–256 says it does not: R-E2 shows the cancellation is O(1) and grows with depth.

**Reusable exact identities:**
- (I1) Euler: a(x) = x·∇a(x), with ∇a piecewise constant on cones.
- (I2) Tropical curvature: E a = E[Δa] = Σ_{k,m} E[δ(z_{k,m}) |∇z_{k,m}|² S_{←(k,m)}]. Equivalently, the total wall mass equals (n − 1) × the spherical mean. Verified at toy scale with a finite-variance kernel estimator, whose bias is O(h).
- (I3) The exact two-wall (bivariate truncated-normal) lift via Owen's T, in `tct0.py`.
- (I4) The temperature bridge: E softplus_T(z) = E relu(z) + (π²/6) T² p(0) + O(T⁴).
  - Its practical content is that the closure error collapses ~400× by T ≈ 0.8 pre-activation std (w64).
  - Neither extrapolation (≤ 2.5×) nor a sampled residual (none) transfers that accuracy back to T = 0.

**What this stream hands to the others:**

- I2 as an exact consistency check: Σ walls = E a, i.e. (n − 1) × the spherical mean;
- the exact two-wall (bivariate) lift in `tct0.py`, which beats the bench's linearised cross-covariance slightly;
- the T-deformation fact that closure error collapses ~400× by T ≈ 0.8 s. Useful only to a design that can transfer finite-T accuracy without analytic continuation. Positivity along depth (the coordinator's suggestion (b)) has no concrete realisation from this principle: the path sum is signed (R-P5), so there is no positive transfer to exploit.
