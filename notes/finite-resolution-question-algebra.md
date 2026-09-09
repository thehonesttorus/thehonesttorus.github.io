# Inputs as states on a finite-resolution question algebra

*What the 0/1 activation code is, in what sense an input is "not a state on the questions", and why finite resolution makes the first natural and the second exact. Specialised to the ARC White-Box Estimation Challenge (WhestBench).*

Research note, 9 September 2026. Written after reading the challenge specification and Phase 2 rules, the ARC companion paper (Wu, Lecomte, Winer, Robinson, Hilton, Christiano, arXiv:2605.05179), the Phase 1 write-up for submission 327801, and the full ChatGPT transcript (43 user turns, 150 assistant turns). Identities marked (✓) were checked numerically by [`checks/finite_resolution_checks.py`](checks/finite_resolution_checks.py); the recorded output is in [`checks/output.txt`](checks/output.txt). Nothing here is a leaderboard result.

---

## 0. Summary

The position to explore was: *an input is not a state on the questions labelled by the parameter values, but a finite-resolution view of the question algebra might realise it naturally.* The outcome is a split verdict, and the split is the whole content.

1. **The input is a state on the parameter-labelled questions, at every resolution.** In the algebra of functions of the input, the weight row `w` labels the question `ℓ_w(x) = w·x`; the input distribution `γ` is a state; the challenge target `E_γ[h_{L,i}]` is literally `γ` restricted to the operator system spanned by the layer-`L` neuron questions. On a finite-resolution truncation (Wiener chaos of degree `≤ K`, which is the Connes–van Suijlekom spectral truncation of the Gaussian Hodge–Dirac triple) the restriction is *exact*, because the vacuum vector is retained. That half of the conjecture is true and needs no new structure.
2. **What fails is transport, not statehood.** A ReLU layer is a unital `*`-homomorphism on the full question algebra (completely positive, linear). But the pullback of a resolution-`K` question has non-zero chaos components of every even degree (Lemma 1), so no polynomial-size resolution is closed under a layer. The induced map on *restricted* states does not exist; every "propagation of moments" is a choice of extension (a closure), and the extension is what is nonlinear. "The input is not CP" is precisely this: finite-resolution states do not compose through ReLU layers. The ambiguity is exact, `2·dist(a, S)`, and it is a commutator norm with the resolution projection.
3. **The 0/1 code is motivated, but as a derived object.** Positive homogeneity (no biases) forces the only meaningful threshold to be zero, makes every activation region a cone, and makes the network linear on each cell. The code algebra is the smallest finite-dimensional subalgebra that resolves the *sign* of every parameter-labelled question; the neuron questions are the degree-one elements of the module it generates. At resolution 1 the sharp code becomes an effect-valued POVM `{e_c}` on `ℂ ⊕ (input space)`: its diagonal in the vacuum is the code law, its **coherences** between the vacuum and the input directions are exactly the responses `∇_h p_c(0) = E[x 1_c]` of the transcript's turn 178, and the mean activations are vacuum-to-weight matrix elements `⟨ξ_w, e_w Ω⟩` (✓). The code carries the amplitudes, but only jointly and only through coherences that a probabilities-only reading discards. That is the sense in which the code is "a state on the parameters": it is not a state, it is a POVM whose pairing with the vacuum state and with the weights gives everything.
4. **Both leading methods in the challenge are finite-resolution truncations**, of the same kind but at different places: ARC's cumulant propagation truncates per layer in the Hermite/chaos basis and closes with a formal Gaussian-reference completion (non-positive for `K ≥ 3`); the leading Phase 1 sampling designs truncate once at the source in the spherical-harmonic basis (a spherical 5-design is a positive state that agrees with the uniform state on harmonics of degree `≤ 5`). The first compounds closure defects with depth (✓, four orders of magnitude over 16 layers at width 256 for `K = 1`; on 16 real Phase 2 networks `K = 1` finishes about 200× worse than budget-matched sampling), the second does not, which is why sampling variants led Phase 1 at 256×32 and why ARC made Phase 2 wider and shallower.

Section 6 lists what this suggests, as directions rather than results.

---

## 1. The challenge in the framework's vocabulary

Bias-free ReLU MLP, `h_0 = x`, `h_ℓ = (W_ℓ h_{ℓ-1})_+`, input `x ~ γ = N(0, I_n)`, weights i.i.d. `N(0, 2/n)`. Phase 2: `n = 1024`, `L = 16`. Target: the `L × n` matrix `Y_{ℓ,i} = E_γ[h_{ℓ,i}(X)]`; ranked on the last row; budget `2^41 ≈ 2.2·10^12` FLOPs per network, about 65,536 forward passes; score `MSE_final · max(0.1, C/B)` (the overview page still shows a 0.5 floor; the Phase 2 rules use 0.1); all arithmetic through `flopscope`.

| Framework object | In the challenge |
|---|---|
| question algebra `𝒜` | `L^∞(ℝ^n, γ)`, bounded functions of the input; linear forms are affiliated unbounded elements |
| state | `γ` (the input law), or `δ_x` (a single input) |
| parameter value `w` (a weight row) | the signed comparison `ℓ_w = w·x`; its spectral projection `P_w = 1{ℓ_w > 0}` (the bit); its positive part `(ℓ_w)_+ = ℓ_w P_w` (the activation) |
| a layer | the unital `*`-homomorphism `α_W : f ↦ f ∘ F_W`, `F_W(h) = (Wh)_+` (Heisenberg picture), dual to the pushforward of states |
| neuron question at layer `ℓ` | `a^{(ℓ)}_j = h_{ℓ,j}` as a function of `x`, i.e. `α_{W_1} ∘ ⋯ ∘ α_{W_ℓ}` applied to the coordinate `j` |
| the target | `γ` restricted to the operator system `S_ℓ = span{1, a^{(ℓ)}_1, …, a^{(ℓ)}_n}`: `Y_{ℓ,·} = γ|_{S_ℓ}` |

Two consequences of positive homogeneity (no biases) are used throughout: `h_ℓ(tx) = t·h_ℓ(x)` for `t > 0`, so the radius integrates out exactly, `E f(X) = c_n E f(U)` with `U` uniform on the sphere and `c_n = √2 Γ((n+1)/2)/Γ(n/2)`; and every joint activation region is a polyhedral cone through the origin at every depth, on which the network is linear, `h_L(x) = A_c x`.

## 2. Where the transcript stands

The general-first phase (turns 44–92) built the objects this note uses: states and questions paired by probability, operations acting forward on states and backward on questions (turn 47); sufficiency and recovery as what a description may forget (55); identification by adding intertwiners rather than deleting distinctions, the code algebra `C*(P_1, …, P_N)` and its commutant `⊕_c B(H_c)` (59); finite-resolution ambiguity, Rieffel's commutator representation of distance to a subalgebra, and the threshold defect (148, 153); the transport measure and certified refinement (155, 160); the code with its response to a change of state (178); the question that determines its own optimal probe (181); and effects rather than bit vectors (184), where the direct answer to "can the code be a state on the parameters?" was: for a fixed input the code is an *effect on parameter space*, evidence against a prior state on parameters, not a state.

The specialisation to the challenge (from turn 93) produced exact reductions first, radial integration, zonoid support functions, signed boundary flux (97), the weights as constructors of observables and the virtual-polytope compiler (123). It then drifted, in turns 118–143, into estimator tuning around an unproved Gaussian closure. The complaint at turns 139 and 144 was correct and the assistant agreed at 148 and 153: the break was "exact reduction → convenient completion → tune". Nothing in the specialisation itself was wrong; the order was. The theory that should have governed the closures is the finite-resolution theory below, and it happens to be exactly what the last user question asks for.

## 3. In what sense is an input "not a state on the questions"?

### 3.1 Two dual pictures

The pairing `⟨w, x⟩ = w·x` is bilinear, so there are two readings.

**Picture A** (evaluation). Questions are functions of the input, labelled by `w`; the input is the state `δ_x`; the activation is `δ_x((ℓ_w)_+)`, the bit is `δ_x(P_w)`, the Gaussian mean is `γ((ℓ_w)_+)`. All layer maps are `*`-homomorphisms. Nothing is missing.

**Picture B** (identification). Questions are functions of the parameter, labelled by `x`: `λ_x(w) = w·x`, an element affiliated with `L^∞(parameter space)`; the parameter value is the state `δ_w`; the input is an *observable*, and the bit `1{λ_x > 0}` is a projection on parameter space. This is the turn 184 reading: a distribution over parameters is the state, the observed code is evidence. Detector tomography (turn 181 §7) lives here.

The two pictures agree on every number, `(w·x)_+ = δ_x((ℓ_w)_+) = δ_w((λ_x)_+)`, and differ on which side is the algebra. "The input is not a state" is true in Picture B and false in Picture A. For the challenge, Picture A is the relevant one, and the input *is* a state on the parameter-labelled questions. So the obstruction the transcript was circling is not statehood.

### 3.2 Where it actually fails: finite-resolution states do not transport

A state `ω` on the operator system `S_W = span{1, ℓ_{w_1}, …, ℓ_{w_n}}` of one layer's parameter-labelled questions is determined by the vector `(ω(ℓ_{w_i}))_i`: states on the linear question system are mean vectors. The activation `(ℓ_w)_+` is not in `S_W`. So the map one would like, "state on `S_ℓ` ↦ state on `S_{ℓ+1}`", must send means to means, and it does not exist as a function:

**Lemma 1 (no finite polynomial resolution is closed under a ReLU layer).** Let `z ~ N(0,1)`. The probabilists' Hermite coefficients of ReLU, `c_k = E[z_+ He_k(z)]`, are

`c_0 = 1/√(2π)`, `c_1 = 1/2`, `c_k = He_{k-2}(0)/√(2π)` for `k ≥ 2`,

so `c_{2m} = (−1)^{m−1} (2m−3)!! / √(2π) ≠ 0` for every `m ≥ 1` and `c_k = 0` for odd `k ≥ 3` (✓). Hence for any linear form `z = w·x` the pullback `α_W(a) = z_+` of a degree-one question has non-zero Wiener-chaos components of every even degree, and `P_K α_W(a) ≠ α_W(a)` for every `K`. The same holds for Hermite expansions about any reference `N(μ, σ²)`.

*Proof.* `He_k φ = (−1)^k φ^{(k)}`, so `c_k = (−1)^k ∫_0^∞ z φ^{(k)}(z) dz`; two integrations by parts leave `(−1)^k φ^{(k−2)}(0) = He_{k−2}(0)/√(2π)`. ∎

Consequences. The restricted state at resolution `K` (moments of the preactivations up to order `K`, equivalently their cumulants up to order `K`) does not determine the restricted state one layer later. Any "propagation" is a choice of extension of the restricted state from `S` to the larger system `S + α_W(S')`. Positive extensions exist and form a convex set; the worst-case ambiguity of a question `a` is exactly

`sup{|φ(a) − ψ(a)| : φ|_S = ψ|_S} = 2·dist(a, S)`

(transcript turns 148 §4 and 153 §2, order-unit duality), and for the actual restricted state the consistent answers form the interval `[sup{s(b): b ≤ a, b ∈ S}, inf{s(b): b ≥ a, b ∈ S}]`.

This is the precise content of the phrase "the input is not CP so not a state on the questions": **the layer is completely positive and linear on the full algebra and the input is a state on every finite-resolution question system, but the finite-resolution states do not compose through ReLU layers.** The nonlinearity of mean propagation, `μ ↦ E[φ(N(Wμ, σ²))]`, is the nonlinearity of a *completion*, not of the underlying dynamics. Sherman's theorem sharpens the picture: `A_sa` is a lattice (so that `a ↦ a ∨ 0 = a_+` is a lattice operation compatible with sums) only when `A` is commutative; in a noncommutative question algebra `a_+` still exists by functional calculus, but the lattice identities that make the piecewise-linear calculus work (turn 123 §8) are gone.

## 4. The finite-resolution question algebra in which the code becomes natural

### 4.1 The Gaussian spectral triple and its truncations

The input state carries its own differential calculus. Let `d` be exterior differentiation on `L²(γ) ⊗ Λ(ℝ^n)` and `d*_γ` its Gaussian adjoint; the Hodge–Dirac operator `D_γ = d + d*_γ` satisfies, on functions,

`D_γ² = N = −Δ + x·∇`,

the Ornstein–Uhlenbeck number operator, whose spectrum is `ℕ` with eigenspaces the Wiener chaoses (multivariate Hermite polynomials of degree `k`). Its commutator with a multiplication operator is Clifford multiplication by `df`, so

`‖[D_γ, M_f]‖ = ‖∇f‖_∞`,

and Connes' distance formula gives the Kantorovich–Rubinstein distance `W_1` for the Euclidean metric: `d(δ_x, δ_y) = |x − y|`, and `d(γ, N(h, I)) = |h|`. This is the "implicit Dirac metric" of the problem. It is not imposed; it is determined by the input state, and the mean-shift probes of turn 178 move the state along its geodesics.

The Connes–van Suijlekom spectral truncation of this triple is, on scalar functions, `P_K` = projection onto chaos of degree `≤ K`, and the truncated question system is the operator system

`Q_K = P_K L^∞(γ) P_K ⊂ B(P_K L²(γ))`,

of matrix size `C(n+K, K)`. For `K = 1` the retained space is `ℂ ⊕ ℝ^n`: the constant plus the input space itself. Two facts make this the right home for the conjecture.

* **The input state is exact on `Q_K`.** `P_K Ω = Ω` for the vacuum `Ω = 1`, so `γ(f) = ⟨Ω, fΩ⟩ = ⟨Ω, P_K f P_K Ω⟩`. Restricting the questions to finite resolution loses nothing about *their* expectations.
* **Sharp questions become effects.** For a projection `P` (a bit, or a code cell) the compression `e = P_K M_P P_K` satisfies `0 ≤ e ≤ 1` and `e − e² = P_K M_P (1 − P_K) M_P P_K ≥ 0` with `‖e − e²‖ = ‖[P_K, M_P]‖²` (turn 184 §3A). Unsharpness *is* the commutator with the resolution.

### 4.2 The code at resolution 1, explicitly

For a bit `P_w = 1{w·x > 0}` with unit normal `ŵ`, in the orthonormal basis `(1, x_1, …, x_n)` of `P_1 L²(γ)`:

```
e_w = P_1 M_{P_w} P_1 = [ 1/2          ŵᵀ/√(2π) ]
                        [ ŵ/√(2π)     (1/2) I_n  ]
```

(✓). Its spectrum is `{1/2 − 1/√(2π), 1/2 (n−1 times), 1/2 + 1/√(2π)} ≈ {0.101, 0.5, 0.899}` (✓); `e_w + e_{−w} = P_1`; and the unsharpness has norm exactly `1/4 = ‖[P_1, M_{P_w}]‖²` (✓). The finite-resolution bit of a half-space is an unsharp effect with a universal spectrum, independent of `n` and of `‖w‖`.

For a cell `c` of the joint code (indicator `1_c`) the compressed effect is the block matrix

`e_c = [[p_c, m_cᵀ], [m_c, S_c]]`, `p_c = γ(c)`, `m_c = E[x 1_c]`, `S_c = E[x xᵀ 1_c]`,

and `{e_c}` is a POVM on `ℂ ⊕ ℝ^n`: `Σ_c e_c = P_1` (✓). Positivity of `e_c` is the classical second-order moment-matrix condition; the compressed cells do not commute (✓, `max ‖[e_c, e_c']‖ ≈ 0.08` on a 4×4×4 example), which is the finite-resolution origin of noncommutativity in this problem, in the sense of turn 184 §5: individual effects are exact, their products are not.

**Mean activations are coherences.** Put `ξ_w = (0, w) ∈ ℂ ⊕ ℝ^n`. Then

`E[(w·x)_+] = ⟨ξ_w, e_w Ω⟩ = ‖w‖/√(2π)` (✓),

and at depth, with `A_c` the linear map of the network on cell `c`,

`E[h_{L,j}] = Σ_c ⟨ξ_{A_cᵀ e_j}, e_c Ω⟩ = Σ_c e_jᵀ A_c m_c` (✓).

This is turn 178's reconstruction identity `E F = Σ_c A_c ∇_h p_c(0)` in operator form: the response of the code law to a mean shift is the coherence row `m_c` of the resolution-1 effect, and the weights (through `A_c`) are what the coherences are paired with. So the exact statement behind "the 0/1 code as a state on the parameters" is:

> At resolution 1 the code is a POVM `{e_c}` on `ℂ ⊕ (input space)`. The input is the vacuum state `Ω`. The code *law* is the diagonal of the POVM in that state. The *amplitudes* are the off-diagonal (vacuum ↔ input) coherences, paired with the parameter-defined linear actions. A reading of the code that keeps only probabilities discards the amplitudes; a reading that keeps the effects keeps them exactly.

### 4.3 Why zero threshold, why a 0/1 code at all

* Homogeneity removes every threshold but zero: `(w·x)_+ = (w·x)·1{w·x > 0}`, and a threshold family `1{w·x > t}` would only be needed with biases. The layer-cake identity `E z_+ = ∫_0^∞ P(z > t) dt` (turn 123 §1) shows what the other thresholds would carry; here they are replaced by the *linear* factor, which is why the module structure below appears.
* The code algebra `C_ℓ = C*(P^{(k)}_j : k ≤ ℓ)` is the smallest finite-dimensional subalgebra of `L^∞(γ)` resolving the sign of every parameter-labelled question up to layer `ℓ`; its spectrum is the set of realised joint codes (cones). The neuron questions are degree-one elements of the `C_ℓ`-module generated by the coordinates: `a^{(ℓ)}_j = Σ_c 1_c · (A^{(ℓ)}_c x)_j`. Everything nonlinear is in the code algebra; everything linear is in the fibre. In geometric language the network is a section of the trivial bundle over the finite space `Spec(C_L)` with fibre `ℝ^n`.
* The tower `C_1 ⊂ C_2 ⊂ ⋯ ⊂ C_L` of refining conic partitions is the exact, commutative, exponentially large finite-resolution question algebra. The state at level `ℓ` needed for means is `(p_c, m_c)_c`. It transports to level `ℓ+1` only if the full conditional law on each cone is kept (a Gaussian truncated to a polyhedral cone; the sub-cell probabilities are Gaussian cone measures, Kendall 1941, Genz–Bretz), not the mean: this is the "code together with its response" theme of turns 178 and 181 (mean shifts give `m_c`, covariance shifts give `S_c`, the full deformation family gives the conditional law).
* In this tower the conditional expectations `E_{C_ℓ}` commute with each other (tower property), so **depth alone creates no noncommutativity** (turn 87 §3). Noncommutativity enters through the resolution projections `P_K`, which do not commute with the multiplication operators, or through incomparable codes (the same input under perturbed weights). The algebra generated by multiplications and conditional expectations is the Jones basic construction `⊕_c B(L²(cell c))`, the transcript's "piece algebra".

So the 0/1 code is motivated, but as a *derived* object: the sign resolution of the parameter questions forced by homogeneity, whose finite-resolution compression is an effect-valued POVM whose coherences carry amplitudes. What was *not* motivated in the early wall/percolation picture was treating the code law `(p_c)` as the state. The code law alone forgets scale (`(cx)_+` has the same code for every `c > 0`, turn 123 §1); the coherences do not.

## 5. Both leading methods are truncations, and this explains the leaderboard

### 5.1 Cumulant propagation is the resolution-`K` state with a formal completion

ARC's algorithm tracks cumulants up to order `K` of the preactivations layer by layer. Moments up to order `K` and cumulants up to order `K` carry the same information: the state restricted to polynomials of degree `≤ K`, i.e. the chaos truncation `Q_K` relative to the layer's reference Gaussian. The linear step is exact (cumulants are multilinear). The nonlinear step needs the state on `α_W(Q_K)`, which by Lemma 1 has all even chaos degrees; the diagram-summation formula truncated at cumulant order `K` is the completion "cumulants above `K` vanish about the reference Gaussian". For `K ≥ 3` that completion is not a probability measure (Marcinkiewicz's theorem), so it is a formal, non-positive extension: precisely the "invented completion" of turn 148 §5. Its per-layer error is a threshold defect `Φ((q)_+) − (Φ(q))_+`, which integrates to the Schwarz defect `Φ(q²) − Φ(q)²` and has norm `‖[P, π(q)]‖²` (turn 148 §2). ARC prove the defect is `O(1/n^K)` at random weights for fixed depth and conjecture `MSE ≲ c_K (L/n)^K`. The telescoping identity of turn 148 §6 makes the total error a sum of transported per-layer defects with no cancellation guarantee, which is the depth breakdown they report empirically.

Illustration (✓, `K = 1`, width 256, depth 16, one He-initialised network, Monte Carlo reference with `10^6` samples; the Monte Carlo floor is `< 10^{-6}` throughout): variance-normalised MSE of mean propagation by layer

| layer | 1 | 2 | 4 | 8 | 12 | 16 |
|---|---|---|---|---|---|---|
| MSE / var | 8.1e-7 | 2.8e-4 | 1.4e-3 | 5.5e-3 | 8.6e-3 | 1.8e-2 |

Layer 1 is exact (the closure is vacuous there); every later layer adds a defect.

On real Phase 2 networks the same picture holds at scale. Sixteen networks of the public `mini` split (`v2-phase2`, 1024×16) were checked in a remote sandbox by regenerating each network's weights from its `mlp_seed` with the dataset's seed protocol (`SeedSequence(seed).spawn(3)[0]` feeding `default_rng`, He scale `√(2/n)`, float32 cast). The regeneration was verified against the exact layer-1 means `‖w_j‖/√(2π)`, which matched the `10^9`-sample ground truth to the Monte Carlo floor (max deviation about `1e-4` over 1024 neurons, the expected extreme of a `2.6e-5` per-neuron floor). `K = 1` mean propagation against the shipped `all_layer_means` ([`checks/phase2_mean_propagation_check.py`](checks/phase2_mean_propagation_check.py), output in [`checks/phase2_output.txt`](checks/phase2_output.txt)) gives

| layer | 1 | 2 | 4 | 8 | 12 | 16 |
|---|---|---|---|---|---|---|
| MSE / mean² (average of 16 networks) | 2.2e-9 | 8.3e-5 | 2.0e-4 | 2.8e-4 | 2.7e-4 | 3.0e-4 |

At the scored layer the average MSE is `2.7e-4` against an average neuron variance of `7.8e-2` (`MSE/var ≈ 3.5e-3`): roughly 180–270× worse than budget-matched naive Monte Carlo (`var/65536 ≈ 1.2e-6`), and five orders of magnitude above the top public score quoted in the transcript (`3.8e-9`, turn 133). Two things are visible here that the width-256 run does not show: the per-layer defect is smaller at width 1024, as ARC's `1/n^K` scaling predicts, and the accumulated error grows only slowly after about eight layers.

### 5.2 A spherical design is the same truncation at the source

After radial reduction the question algebra is `C(S^{n−1})` with the uniform state, whose Dirac operator has the spherical harmonics as eigenspaces. A spherical `t`-design is a finitely supported *positive* state that agrees with the uniform state on all harmonics of degree `≤ t`: the truncation is applied once, to the source, and the estimator evaluates the whole network at each node. The Phase 1 write-up for submission 327801 uses 129 Kerdock mutually unbiased bases and their antipodes, 66,048 unit inputs forming a spherical 5-design, with moment-repair corrections on top. Its error is the degree->5 harmonic content of the network function on the sphere; depth changes that spectrum but does not compound the error, because no closure is applied between layers. The organisers' Phase 2 note says the best Phase 1 methods "tended to be based on variations of sampling", and that the network was made wider and shallower "to make cumulant-propagation-style methods more competitive": exactly the trade-off between a truncation applied per layer with closure (error accumulates in `L`, shrinks in `n`) and a truncation applied once (error depends on the spectrum of the compiled observable).

### 5.3 The ceiling in one line

For a question `a` and a retained system `S`, the worst-case error of *any* estimator that only sees the state on `S` is `dist(a, S)`. For the polynomial truncation of a single ReLU'd Gaussian linear form of range `R`, `dist` decays like `R/K` (Bernstein's constant for `|x|`); for the design it is the distance to harmonics of degree `≤ t`. The state-specific interval is far tighter at random weights, which is why both methods work at all; the interval, not the point estimate, is what a certified method would report (turn 160).

## 6. What the finite-resolution view suggests

Directions, not results.

1. **Attach the resolution to the state-weighted action, not to the algebra** (turn 174). The truncation subspace at layer `ℓ` should be spanned by what the next layer asks: the `n` preactivation directions `w'_k·h` and their threshold structure, not an isotropic polynomial basis. The minimal resolution that transports *means* through one layer is the set of one-dimensional marginal laws of the next preactivations (`E z_+ = ∫_0^∞ P(z > t) dt`); through two layers it is their joint law. Polynomial resolutions are the wrong basis for a threshold nonlinearity: turn 155 derives an `ε³` boundary-layer discrepancy that survives matching moments of any fixed order (that derivation is the transcript's, checked there symbolically, not re-derived here).
2. **Propagate compressed effects rather than moments.** The resolution-1 effects `e_k = [[p_k, m_kᵀ],[m_k, S_k]]` of the next layer's bits are the objects whose coherences carry amplitudes, and the POVM constraint `Σ_c e_c = P_1` (joint measurability, turn 184 §5) is a consistency condition that Gaussian cell-closures violate. By Lemma 1 any such scheme is again a closure, but one whose defect is the commutator `‖[P, α_W(a)]‖` with an explicit resolution projection, hence measurable rather than guessed.
3. **Certify.** The positive-extension interval can accompany any closure; the last layer admits cheap brackets `b⁻ ≤ F ≤ b⁺` (turn 153 §7). Certificates cost FLOPs and the score rewards accuracy per FLOP, so they pay only if they steer refinement (turn 160's monotone certified refinement is the model).
4. **Where the Dirac operator is.** If the "Dirac metric at finite resolution" the transcript kept reaching for has a precise home in this problem, it is twofold: the Gaussian Hodge–Dirac `D_γ` that the input state supplies (its distance is `W_1`, its spectral truncations are the chaos resolutions), and the grading `2P − 1` of a resolution projection, whose commutators with the parameter-labelled questions are, by Rieffel's theorem, the best-approximation errors. The invariant controlling every closure in the challenge is `‖[P, M_q]‖` for the questions `q` the next layer asks.

## 7. Limits of this note

* No estimator is proposed and no score is claimed. The transcript's best own estimator was about 57× above the then-top public Phase 2 score by its own comparison (turn 133).
* The identity checks use small synthetic networks; the official dataset was not reachable from the environment that ran them. The Phase 2 comparison in §5.1 ran in a separate sandbox on 16 networks of the public `mini` split, using only the shipped ground-truth means and seeds (weights regenerated and verified as described there). It covers `K = 1` only; ARC's `K ≥ 2` variants were not re-implemented.
* The transcript's attached files (sandbox links, Lean corpora, Mac experiment directories) were not accessible; their claims were not verified.
* Only the main text of the ARC paper was read; statements about the method's internals follow it and the published pseudocode.

## Appendix: references used

Connes–van Suijlekom, *Spectral truncations in noncommutative geometry and operator systems* (CMP 2021) and *Tolerance relations and operator systems* (2022); van Suijlekom, *Gromov–Hausdorff convergence of state spaces for spectral truncations* (2021); Rieffel, *Leibniz seminorms and best approximation from C*-subalgebras* (2011); Sherman (1951) on lattice-ordered self-adjoint parts; Marcinkiewicz (1939) on polynomial cumulant generating functions; Kendall (1941) and Genz–Bretz (2009) on Gaussian orthant probabilities and the tetrachoric series; Wu et al., arXiv:2605.05179; the WhestBench dataset card (`aicrowd/arc-whestbench-public-2026`, revision `v2-phase2`); the AIcrowd Phase 2 announcement (23 August 2026); Phase 1 write-up for submission 327801.
