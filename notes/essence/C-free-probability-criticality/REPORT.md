# Team C: free probability, two-projection geometry and the wall at 2 — REPORT

*Essence programme, round 3. v1, 1–2 Oct 2026 (v0 ≈ 22:15 UTC; v1 adds the T1/T2/T4 tests at n = 1024). Labels: **Theorem** (published, with source), **Derived** (proved here, elementary), **Measured** (this directory, numbers below), **Synthesis**, **Conjecture**. Every n = 1024 number is on the six bench networks `w1024_d16` unless stated otherwise.*

## 0. Answer in brief

1. **The process in this domain.** In this domain "forget one piece, re-randomise it, converge" is alternation between two conditional expectations, i.e. two projections P and Q. Halmos's theorem splits the space into 2×2 blocks indexed by an angle θ. In each block alternating projections contract by cos²θ, and alternating reflections rotate by 2θ. The trace of a 2×2 block (A'Campo: λ + λ⁻¹ + 2 = μ²) gives three regimes:
   - elliptic (|tr| < 2): periodic;
   - parabolic (tr = 2): unipotent, a Jordan block, linear growth;
   - hyperbolic (|tr| > 2): exponential, wild, mixing.
2. **Where the network sits.** The wall at 2 is literally the order–chaos boundary χ₁ = σ_w² E[φ′²] = 1 of deep information propagation, and He-initialised ReLU sits on it (χ₁ = 2 · ½ = 1).
   - *Pathwise* (one input), every layer is an isometry in mean: tr(JᵀJ)/n = 1 per layer. This is the parabolic case. Its signatures are linear-in-depth spectral variance (Pennington–Schoenholz–Ganguli: σ²_{JJᵀ} = 2L for Gaussian ReLU), log-normal norms with variance ∝ depth (Hanin–Nica: β = 5 Σ 1/n_l), and polynomial (≈ 9π²/(2l²)) rather than exponential convergence of correlations.
   - *Annealed* (gates averaged over the input, which is what every moment method transports), the per-layer gain on a generic direction is strictly inside the wall. **Derived:** g_l = 2 E_i[Φ_i²] = ½(1 + E_i[(2Φ_i − 1)²]) = χ₁ − 2 E_i[Φ_i(1 − Φ_i)]. So the forgetting rate of a layer is exactly twice the mean Bernoulli variance of its gates, the Jensen gap between the pathwise and the annealed gate.
3. **The decisive measurement (seed 2), done here directly at n = 1024 on all six networks** (`agespec.py`, the exact first-order co-state, all 120 source–kink pairs). Third-order content splits into two sectors with different algebra:
   - **Dilation (unipotent) sector.** The amplitude γ of each pair's (2,1) slice on the scale-mixture template K = 2 m ⊗ C + diag(C) ⊗ m is retained at 0.86 → 0.95 per step and approaches 1 with age (under the decoupled-gate transport; exactly 1 in the true dynamics by costate C6). It adds **coherently** across ages: (Σ_s γ_s)² / Σ_s γ_s² = 0.88–0.94 × the number of ages.
   - **Free (wild, mixing) sector.** The deflated residuals of different ages are mutually orthogonal (mean cross-age cosine ≤ 0.03; coherent/incoherent energy ratio 0.94–1.06), i.e. asymptotically free. Each is transported with the free-probability Frobenius gain: energy ratio per step = **(1.02–1.10) × g_l³**.
   - **No periodic sector.** No elliptic sector exists for generic weights; costate C7 says the dilations are the only continuous symmetry.
   - **The odds law** (new, measured). At every layer k = 1…15, the dilation share s_k of the total old-and-young D21 energy satisfies **s_k / (1 − s_k) = (0.31 ± 0.01) × k**, from 0.30 to 0.33 at each layer, mean over 6 networks. This is the coherent N² sum of a Jordan sector against the incoherent N sum of a free sector.
4. **The free law of the gated propagators holds at n = 1024 to 1 %** (measured):
   - tr(UᵀU)/n = 2 Π g_l to 1e-4 at ages ≤ 5;
   - PR = n/(2(age + 1)) to ≤ 1.1 % at ages 0–5.
   - Beyond age 6 the Perron outlier (the mean direction, F2.2) breaks freeness: PR falls 5 % below at age 9 and 20 % below at age 14, and the trace exceeds the free product by 5 %.
   - The F9.1 mismatch at ages 9–13 is this rank-one spike, i.e. a BBP-type outlier of a spiked free model.
5. **Design consequence (seed 2), priced.**
   - The part that does not mix is exactly the dilation/Perron sector, and in the true dynamics it is a single conserved scalar (one Jordan block). Everything else mixes geometrically at rate g³ ≈ 0.13 (layer 0) to 0.6 (depth) per step, with no spectral gap needed: the forgetting is the gate variance.
   - So: carry the dilation scalar exactly at O(n²), and truncate the free sector at an age set by the cumulative free trace Π g³, not by a fixed age.
   - At the final layer the free sector beyond age A holds 5.6 % (A = 3), 2.9 % (A = 5) and 1.3 % (A = 7) of the total D21 energy (measured).
   - The window costs ≈ 3–4 products per pair (costate prices). Costate's A3gsl_nc (273 products, raw 3.25e-7) is this design with a uniform window.
   - **Measured (v1).** A window set by the cumulative free trace Π (2EΦ²)³ ≥ τ beats the uniform window on all six networks at n = 1024:
     - raw **2.43e-7** at 308 products (τ = 0.05, 0.75× A3gsl_nc), the lowest raw of any first-order design in the campaign;
     - raw 2.68e-7 at 269 products (τ = 0.1, 0.83×);
     - best adjusted ≈ **4.4e-8** (raw 2.91e-7 at 226 products, Strassen L3 price), against costate's best of 5.2–5.9e-8.
   - The other predicted refinement, a leak-corrected dilation amplitude (T1), is killed. The early decay of γ is a transient onto the invariant direction, not a leak.
   - T4 confirms the Perron diagnosis: with the mean gates permuted (traffic-independent), PR follows the free law to 2 % at every age.
6. **What this does not do.** It explains the first-order old content and prices it. It does not reach the bar. Costate §4.6 shows the remaining 15× (at w128) is joint κ₄, and that is again concentrated in the dilation sector. The theory here says why: the dilation sector is the only parabolic (non-mixing) sector, so it is the only one whose higher cumulants accumulate with depth (Var log t ∝ depth).

## 1. Instances

### I1. Two projections: Halmos, the dihedral algebra, alternating projections, and the free pair

- **Statement (Theorem; Halmos 1969; Raeburn–Sinclair, Math. Scand. 65 (1989)).**
  - Two projections P, Q in general position are unitarily equivalent to P = [[1, 0], [0, 0]] and Q = [[c², cs], [cs, s²]] ⊗ (functions of an angle operator 0 < Θ < π/2), with c = cos Θ, s = sin Θ.
  - The universal C*-algebra C*(p, q) = C*(Z₂ * Z₂) (the infinite dihedral group, generated by the reflections 2p − 1 and 2q − 1) is {f ∈ C([0, π/2], M₂) : f(0), f(π/2) diagonal}.
  - Its irreducible representations are the angles θ (2-dimensional) plus four 1-dimensional ones.
- **The process.** The down step is P (forget the complement of ran P); the up step is Q.
  - **Alternating projections:** (PQ)ᵏ → P∧Q at rate cos^{2k} θ_F, where θ_F is the Friedrichs angle (von Neumann 1933; Deutsch's book, Best approximation in inner product spaces, ch. 9).
  - **Alternating reflections:** (2P − 1)(2Q − 1) acts on each block as rotation by 2θ.
  - **Payoff.** Kaczmarz / ART (CT reconstruction), POCS, ADMM, and the Strohmer–Vershynin randomised Kaczmarz with expected rate 1 − κ(A)⁻² (J. Fourier Anal. Appl. 2009). The angle is the condition number.
- **The free pair (Theorem, standard; Nica–Speicher Lectures, free multiplicative convolution of Bernoulli laws).**
  - Take two free projections of trace ½. The spectral law of PQP on ran P is the arcsine law dx/(π√(x(1 − x))) on [0, 1]. Equivalently, **the Halmos angle θ is uniform on [0, π/2]**.
  - Hence the product of the two symmetries, (2P − 1)(2Q − 1), is a Haar unitary in distribution: the eigenvalues e^{±2iθ} are uniform.
- **Derived.** For the free half-projections, alternating projections do not mix geometrically.
  - τ((PQP)ᵏ)/τ(P) = ∫ cos^{2k}θ dθ/(π/2) = C(2k, k)/4ᵏ ~ (πk)^{−1/2}.
  - The angle law has no gap at θ = 0, so the decay is polynomial: a critical (parabolic) behaviour manufactured by freeness alone.
  - This is the two-projection face of "no spectral gap, deep in the free regime" (F9.2).

### I2. Two reflections as a Coxeter element: A'Campo, Zamolodchikov, Smith, Jones, Kruglyak–Roiter

- **Statement (Theorem).** Take a bipartite graph Γ with adjacency eigenvalues μ. The bipartite Coxeter element c = s₊ s₋ (each factor a product of commuting reflections) has eigenvalues λ with λ + λ⁻¹ + 2 = μ² (A'Campo, Invent. Math. 33 (1976)). With μ = 2 cos θ, λ = e^{±2iθ}.
  - **Dynkin (spectral radius < 2, Smith's theorem):** c has finite order h, the Coxeter number. Zamolodchikov periodicity of Y-systems follows (Keller, Ann. Math. 177 (2013), arXiv 1001.1531). There are finitely many indecomposables: Gabriel's theorem, and for Hilbert-space representations, Kruglyak–Roiter (arXiv math/0307163).
  - **Affine (spectral radius = 2):** c has eigenvalue 1 with a Jordan block of size 2. The null root δ is the eigenvector and the defect ∂ the linear functional, so cᵏ x = x + k ∂(x) δ on the relevant lattice (Dlab–Ringel, Mem. AMS 173 (1976)). Preprojective dimension vectors grow linearly and regular tubes are periodic. Representation type: tame.
  - **Wild (spectral radius > 2):** c has a real eigenvalue λ > 1 (Ringel, Math. Ann. 300 (1994)). Dimension vectors grow exponentially. Representation type: wild.
  - **The same wall elsewhere.** Jones: subfactor indices below 4 are 4 cos²(π/n). Kruglyak–Rabanovich–Samoilenko: the set of α for which n projections sum to αI is discrete for n ≤ 4 and contains a band for n ≥ 5 (the formula quoted in the input file is not verified by us).
- **Payoff.** Finite type gives finite classification and exact algorithms. The affine wall is the boundary of tame classification. The computational content is classification complexity, not running time.

### I3. Criticality of deep networks: dynamical isometry and products of random matrices

- **Statement (Theorem).**
  - **Order and chaos (Poole et al., arXiv 1606.05340; Schoenholz et al., arXiv 1611.01232).** The correlation map c ↦ f(c) has slope χ₁ = σ_w² E[φ′(h)²] at c = 1. χ₁ < 1 is ordered, χ₁ > 1 chaotic, and the correlation depth scale ξ_c = −1/log χ₁ diverges at χ₁ = 1.
  - **Pennington–Schoenholz–Ganguli (arXiv 1711.04735; 1802.09979).** For J = Π D_l W_l with active fraction p, σ²_{JJᵀ} = L(μ₂/μ₁² − 1 − s₁) at χ = 1. For ReLU this is **2L** with Gaussian weights and **L** with orthogonal weights: ReLU cannot achieve dynamical isometry, and the spectral variance grows linearly in depth.
  - **Hanin–Nica (arXiv 1812.05994).** (n₀/n_d)‖Mu‖² ≈ exp N(−β/2, β), with β = 5 Σ 1/n_l for He ReLU.
- **Derived.** For He ReLU, f(c) = (√(1 − c²) + (π − arccos c) c)/π. With ε = 1 − c, the map is ε ↦ ε − (2√2/(3π)) ε^{3/2} + …, so 1 − c_l ≈ 9π²/(2 l²). The convergence is a parabolic fixed point with polynomial approach. The subagent also derived this; it is not checked against a published source (Hayou et al. 2019 is believed to state a rate of this form).
- **The three regimes.** Map the A'Campo trichotomy onto information propagation:
  - chaotic ↔ wild (positive Lyapunov exponent, exponential separation, expanders);
  - critical ↔ affine (unipotent, linear growth of variances, polynomial rates);
  - ordered ↔ hyperbolic contraction (|λ| < 1), *not* the elliptic side.

  The elliptic (periodic, Dynkin) side has no counterpart in a real 1-dimensional correlation map. It needs a finite-order symmetry, and a generic ReLU network has none (costate C7).
- **Payoff.** Depth scales for trainability. Free probability gives the full Jacobian spectrum through S-transforms in O(1) work, instead of O(n³) per product.

## 2. The noncommutative statement

**Known (Theorem):** the two-projection C*-algebra (I1); Coxeter functors and the A'Campo correspondence (I2); free multiplicative convolution, S(ab) = S(a)S(b), and second-order freeness (Mingo–Speicher, arXiv math/0405191: traces of products of independent Gaussian, Wishart or Haar matrices have O(1) Gaussian fluctuations with a universal covariance). Male's traffic independence (arXiv 1111.4662) extends asymptotic freeness to permutation-invariant families. Diagonal matrices with i.i.d. entries are traffic-independent of Gaussian matrices (Cor. 1.9).

**The setting.** Layer algebra A_l: functions of z_l; the content space is Sym³ (third cumulants), on which the transfer acts by T_l κ = κ(M_l ·, M_l ·, M_l ·), M_l = diag(Φ_l) W_{l+1}.
- *The down step* is the gate's conditional expectation over the input, E_x[D_x] = diag(Φ), a projection-valued random variable replaced by its mean.
- *The up step* is the fresh weight W_{l+1}. One layer makes the source's orientation Haar relative to the downstream propagator (F9.6).
- *The angle* of the layer is cos² θ_l := g_l = 2 E Φ², the Frobenius gain of one leg.

**Proposition C1 (Derived; exact identities).**
- (a) Pathwise, E_x tr(J_xᵀ J_x)/n = 2 E_x E_i[D_{x,i}] = 1 (χ₁ = 1), so every layer sits on the wall.
- (b) Annealed, g_l = tr(M_lᵀ M_l)/n → 2 E_i[Φ_i²] = ½(1 + E_i[(2Φ_i − 1)²]) = 1 − 2 E_i[Φ_i(1 − Φ_i)] as n → ∞ (the limit uses E W_{ij}² = 2/n).
- (c) Hence the forgetting of a generic leg per layer is twice the mean gate variance, and in Halmos terms the layer's angle is the gate's own two-projection angle. With Φ_i = cos²θ_i, 2Φ_i − 1 = cos 2θ_i. A frozen gate (θ ∈ {0, π/2}) forgets nothing; a fair gate (θ = π/4) forgets half.
- (d) Under the exact dynamics the dilation direction is an eigen-direction of eigenvalue exactly 1 (costate C6: Euler's φ′(z)z = φ(z) makes the annealed step inherit the pathwise criticality on that one direction). With per-layer births the pair (charge, accumulated scale variance) evolves by [[1, 0], [β_l, 1]]: a unipotent Jordan block of size 2, the affine (μ = 2) case of I2.

**Conjecture C2 (the critical trichotomy for quenched third-order content).** Setting: He-initialised bias-free ReLU, width n → ∞, depth L fixed, first-order (decoupled-gate) transport of third-order content. Split the content at each layer as γ K_k ⊕ R, where K_k is the dilation template and R ⊥ K_k. Then:
- (i) **Parabolic sector.** The dilation amplitudes of different ages add coherently. Each source's overlap with K relaxes onto the invariant direction: its retention rate ρ_D(a) ↑ 1 with age a, and the early ρ_D < 1 is the decay of the non-invariant transient (T1 rules out reading it as a leak). On the invariant direction itself the eigenvalue is exactly 1 (C6).
- (ii) **Free sector.** The residuals R_{s→k} of different sources s are asymptotically free of second order: their normalised Gram matrix tends to the identity, with off-diagonals O(n^{−1/2}) plus a rank-one Perron correction. Each one's energy obeys E‖R_{s→k+1}‖² = g_k³ E‖R_{s→k}‖² (1 + O(Perron)).
- (iii) **No elliptic sector** exists for generic weights.
- (iv) **Perron correction.** The only non-free correction is the mean direction. The mean gate Φ_{l+1} depends on W_{l+1} through m_{l+1} = E[a_l] W_{l+1}, a rank-one coupling, so the propagator is a spiked free product: its PR falls below n/(2(age + 1)) once the spike separates (BBP), which happens at age ≈ 6–8 at n = 1024.

Measured support: §4, every clause at n = 1024 to the stated accuracy.

**The precise next theorem (Conjecture C3, second-order freeness of memory).** Under the hypotheses of C2, the old content of age a at layer k, O_{a,k} := Σ_{s = k−a−1} S_{s→k}, has E_W O = 0 (F8.5). Its covariance over the weight ensemble is given by the second-order free cumulants of the gated product, i.e. by the S-transform data of {WᵀW, D_l²} alone.
- Consequently ‖O_{a,k}‖² / ‖birth_s‖² = Π_{l=s+1}^{k−1} g_l³ · (1 + o(1)) in the bulk.
- The readout-relevant norm of the memory is predicted from n numbers per layer (the Φ_l), without computing any tensor.

This is the noncommutative version of "old content is a quenched fluctuation": it lives at the level of second-order freeness, which is exactly where annealing (first-order freeness) sees zero.

## 3. Transfers to our problem

Common dictionary (anti-naivety rule). The pieces are **age sectors of the third-order content** (one per birth layer), split into the dilation template and the free residual. They are tensors in Sym³, not neurons. Per-neuron gates enter only through **one scalar per layer**, g_l, and one vector, the Perron direction. The down step is annealing the gate; the up step is re-randomisation by the fresh weight; the angle is g_l; the global quantity is the D21 slice at the kink layer and, through it, the readout.

### T1. Leak-corrected dilation charge — **killed at w128 (§4.5)**

*Result.* Scaling the retiring amplitudes by 1.33–2.5 hurts: A3gsl_nc goes from 2.56e-5 to 3.4–4.8e-5, and A1gl_nc from 3.5e-5 to 4.4–6.8e-5. The optimum is the uncorrected projection (factor 1–1.25). The reading below was wrong and has been corrected in C2(i). The early decay of γ (0.86 per step) is not a leak of a conserved quantity. It is the transient of the non-invariant part of a star atom's overlap with K. Only a true scale mixture is invariant (C6), and the retention ratio rising to 0.95 with age is the content *relaxing onto* the eigenvalue-1 direction. Projecting at retirement, after the transient, is therefore right, and the self-overlap and κ₄-visible scale variances (costate §4.7) are a different object: the second-order, correlated scale field.

*Original prediction, kept for the record:*

- **Prediction.** Decoupled-gate transport leaks the dilation amplitude at 0.86, 0.87, 0.88, 0.89 per step over ages 0→4 (measured), while the exact dynamics conserves it (C6). The costate carrier projects content when it retires at age A + 1, after the leak, so it under-counts the scale variance by Π ρ_D ≈ 0.58 at A = 3 and 0.75 at A = 1.
  - This explains part of costate §4.7: the κ₃-visible v (≈ 0.02) is below the self-overlap v (0.074) and the κ₄-spike v (0.16).
  - Predicted fix: multiply each retiring amplitude by 1/Π_{a ≤ A} ρ_D(a), from the measured ρ_D table, or project each source at age 1 and conserve the result from then on.
- **Test.** A3gsl_nc and A1gl_nc with the leak correction, w128 first, then 6 × w1024.
- **Cost.** Zero extra products.
- **Kill.** Raw at 1024 not better than 3.25e-7 (A3gsl_nc) at equal products.
- **Status.** Queued for v1.

### T2. Trace-adaptive window for the free sector (seed 2's design, priced) — **confirmed at n = 1024 (§4.5)**

*Result.* Rule: keep pair (s, k) exact while Π_{l=s+1}^{k} (2 E Φ_l²)³ ≥ τ, otherwise retire it into the dilation charge (and, for "gsl", the slice chain). This beats costate's uniform A = 3 window on all six networks at equal or fewer products.
- τ = 0.1: raw 2.68e-7 ± 0.13e-7 at 269 products, 0.83× (0.71–0.93), against 3.25e-7 at 273.
- τ = 0.05: raw **2.43e-7 ± 0.11e-7** at 308 products, 0.75× (0.60–0.86). This is the lowest raw of any first-order design measured in the campaign; the exact all-pairs first order is 4.0e-7 at 855 products.
- Best adjusted: gl with τ = 0.1, raw 2.91e-7 at 226 products, **≈ 4.4e-8 adjusted** at Strassen L3 (0.15 B). Costate's best was 5.2–5.9e-8.

The only input is one scalar per layer, g_l = 2 E Φ_l², which is the free-probability trace law (tr UᵀU/n = 2 Π g to 1e-4). The design gain is theory-driven: the free sector forgets at the gate's Jensen gap, which varies 5× in rate across depth, so a window measured in ages mis-allocates.

*Original prediction:*

- **Prediction.** The free residual of source s at kink k has energy ∝ Π_{l=s+1}^{k−1} 1.03 g_l³ (measured law). Early sources die fast (g ≈ 0.5–0.7, so 0.13–0.35 per step); late sources die slowly (g ≈ 0.8–0.97, so 0.5–0.9 per step).
  - A uniform age window wastes products on early sources and drops late ones too soon.
  - Rule: carry pair (s, k) exactly while Π_{l=s+1}^{k−1} g_l³ ≥ τ, and fold the rest into the dilation charge plus the slice chain.
- **Cost model.** At 3–4 products per live pair, a budget equal to A3's.
- **Test.** At equal products, trace-adaptive window against the uniform A = 3 window, on 6 × w1024.
- **Kill.** No gain at equal products.
- **Explains.** "Every age matters" (region N6) and the measured slow age decay at depth (F8.2): g → 1 as the gates freeze.

### T3. Memory norm from free probability alone (seed 5)

- **Prediction.** The energy of the free-sector memory at layer k is Σ_s ‖birth_s‖² Π g³. The dilation share follows the odds law s/(1 − s) ≈ 0.31 k. This gives an a-priori error budget for any truncation from the n-vectors Φ_l alone, so the size of the binding object can be predicted before any tensor is computed.
- **Measured.** The residual transport is 1.02–1.10 × g³, and the odds constant is 0.30–0.33 at every layer.
- **Kill.** If the readout-weighted (not Frobenius) tail fails to follow the law; costate §4.2 found the readout weighting differs strongly from Frobenius for the scale mode.

### T4. Spiked free product: deflate the Perron direction from the propagators — **diagnostic confirmed at n = 1024 (§4.6)**

*Result.* With the mean gates randomly permuted per layer, PR follows the free law n/(1 + (a + 1) + Σ_l (r(Φ_l²) − 1)) to ≤ 2 % at every age 1–13. With the true gates it falls to 0.77× the law at age 13. With one input's pathwise 0/1 gates it is flat at n/(2(a + 1)) (0.98–1.04).
- So the sub-free tail of F9.1 is exactly the coupling of the annealed gate to W through the mean, a failure of Male's traffic-independence hypothesis by one rank-one direction.
- It is an annealed phenomenon: pathwise gates show no outlier.
- The v0 caveat came from a wrong guess: r(Φ²) is 1.6–2.1 at n = 1024 (the mean gates are nearly frozen), not 1.2.

*Original statement:*

- **Prediction (C2(iv)).** The propagators are a free product plus one outlier. The outlier both breaks the free law (PR −20 % at age 14) and carries the dilation template's m-leg (scale/g³ ≈ 1.2 per step against 1.03 for the bulk).
  - Splitting U = u₁σ₁v₁ᵀ + U_bulk makes the bulk exactly free, so second-order free formulas apply to it.
  - The Perron part is rank one: O(n²) per pair instead of n³ for that leg.
- **Test.** At n = 1024, permute Φ_l independently of W (this kills the coupling); PR should return to the free law at all ages.
- **Kill.** PR still falls below the free law with Φ permuted.

## 4. Tests run (all in this directory)

`agespec.py` reuses costate's validated source atoms and Stein injection (exact first order, no coincidence atoms, all pairs) and records every (s, k) pair. `aggregate.py` → `results/aggregate_w1024.txt`. Run time ≈ 31 s per network at n = 1024.

**4.1 The odds law and coherence (6 × w1024).**

| k (ages live) | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 15 |
|---|---|---|---|---|---|---|---|---|---|
| dilation share of D21 energy | 0.24 | 0.38 | 0.47 | 0.56 | 0.64 | 0.72 | 0.76 | 0.79 | 0.82 |
| odds / k | 0.32 | 0.30 | 0.30 | 0.31 | 0.30 | 0.32 | 0.32 | 0.32 | 0.31 |
| dilation coherent/incoherent ÷ k | 1 | 0.94 | 0.92 | 0.90 | 0.89 | 0.89 | 0.90 | 0.91 | 0.92 |
| residual coherent/incoherent | 1 | 0.97 | 0.96 | 0.96 | 0.94 | 0.94 | 0.95 | 0.99 | 1.06 |
| mean cross-age residual cosine | — | −0.04 | −0.02 | −0.00 | +0.00 | +0.01 | +0.01 | +0.02 | +0.03 |

The w128 control (network 0) is the same picture, noisier: odds grow to 0.93 share at k = 15, and residual coherence rises to 1.8 at depth. At small width the free sector is less free (pre-asymptotic, as the BRIEF warns).

**4.2 Per-step transfer rates by age (6 × w1024, medians over sources).**

| step (age) | 0→1 | 1→2 | 2→3 | 3→4 | 5→6 | 7→8 | 9→10 | 11→12 | 13→14 |
|---|---|---|---|---|---|---|---|---|---|
| free residual energy ÷ g³ | 1.03 | 1.02 | 1.02 | 1.03 | 1.02 | 1.04 | 1.05 | 1.08 | 1.10 |
| dilation energy ÷ g³ | 1.18 | 1.20 | 1.19 | 1.17 | 1.23 | 1.22 | 1.23 | 1.27 | 1.30 |
| dilation amplitude γ retention | 0.86 | 0.87 | 0.88 | 0.89 | 0.92 | 0.92 | 0.93 | 0.94 | 0.95 |

g_l = 2E[Φ²] over layers 0…14 (network 0): 0.50, 0.58, 0.63, 0.71, 0.75, 0.78, 0.76, 0.82, 0.89, 0.83, 0.84, 0.84, 0.97, 0.87, 0.85. The other five networks agree to ±0.05.

**4.3 Free law of the gated propagators (6 × w1024).**

| age | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 | 14 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PR ÷ (n/(2(age + 1))) | 1.00 | 1.01 | 1.01 | 1.01 | 1.00 | 1.00 | 0.99 | 0.97 | 0.93 | 0.87 | 0.80 |
| tr(UᵀU)/n ÷ 2Πg | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.002 | 1.007 | 1.016 | 1.032 | 1.052 |

Resolved in v1 (§4.6). The v0 caveat assumed r(Φ²) ≈ 1.2. The measured r(Φ²) is 1.6–2.1 (the mean gates are nearly frozen), so the free law ≈ n/(2(age + 1)) is predicted at young ages. The old-age deficit disappears when the gates are permuted, which shows it is the Φ–W (Perron) coupling.

**4.4 Free-sector tail at the final kink layer (6 × w1024).** Share of energy in ages > A:

| A | 0 | 1 | 2 | 3 | 5 | 7 |
|---|---|---|---|---|---|---|
| of the free-sector energy | 0.77 | 0.59 | 0.45 | 0.33 | 0.17 | 0.08 |
| of the total D21 energy | 0.13 | 0.10 | 0.075 | 0.056 | 0.029 | 0.013 |

**4.5 Design tests (T1, T2; `costate_c.py` = costate.py plus two switches, `t1.py`, `summ_t.py`).**
- `costate_c.py` adds `leakfac` and `tau` to costate.py; with both at their defaults it reproduces costate exactly (A3gsl_nc 3.247e-7 = costate's 3.25e-7).
- Raw = final-layer MSE − truth noise, paired on the same networks.
- Adjusted = raw × max(0.1, 0.683 u × products / 1024 u), i.e. Strassen L3 at costate's prices. That is a price, not a measured wall time.

| variant (6 × w1024) | raw | products | vs A3gsl_nc (range) | wins | adjusted (L3) |
|---|---|---|---|---|---|
| costate A3gsl_nc (uniform window A = 3) | 3.25e-7 ± 0.17e-7 | 273 | 1 | — | 5.9e-8 |
| **Tgsl_nc, τ = 0.05** | **2.43e-7 ± 0.11e-7** | 308 | 0.75 (0.60–0.86) | 6/6 | 5.0e-8 |
| Tgsl_nc, τ = 0.03 | 2.53e-7 | 334 | 0.78 | 6/6 | 5.6e-8 |
| Tgsl_nc, τ = 0.1 | 2.68e-7 ± 0.13e-7 | 269 | 0.83 (0.71–0.93) | 6/6 | 4.8e-8 |
| Tgsl_nc, τ = 0.15 | 2.83e-7 | 245 | 0.88 | 5/6 | 4.6e-8 |
| **Tgl_nc, τ = 0.1** | 2.91e-7 ± 0.15e-7 | **226** | 0.90 (0.80–1.04) | 5/6 | **4.4e-8** |
| Tgl_nc, τ = 0.2 / 0.3 | 3.9e-7 / 5.8e-7 | 182 / 148 | 1.21 / 1.80 | | 4.7e-8 / 5.8e-8 |

At w128 (8 networks): Tgsl_nc τ = 0.1 wins 8/8 (0.83×), τ = 0.15 wins 7/8 (0.82×). T1, leak factor on the retiring amplitude (w128, A3gsl_nc, geometric means): 0.5 → 2.93e-5, 0.75 → 2.53e-5, **1 → 2.30e-5**, 1.25 → 2.27e-5, 1.7 → 2.93e-5, 2.5 → 4.22e-5. The optimum is at 1–1.25, so T1 is killed.

**4.6 Traffic freeness of the gates (T4, `t4.py`, network 0, w1024, chain from s = 1).** Each entry is PR × 2(a + 1)/n at age a, with the free law (a separate row), computed from the measured r(Φ_l²):

| age | 1 | 2 | 5 | 7 | 9 | 11 | 13 |
|---|---|---|---|---|---|---|---|
| true mean gates | 1.10 | 1.09 | 1.06 | 0.99 | 0.95 | 0.88 | 0.77 |
| permuted mean gates | 1.10 | 1.09 | 1.07 | 1.06 | 1.05 | 1.02 | 1.01 |
| free law with r(Φ²) | 1.10 | 1.10 | 1.07 | — | — | — | 1.03 |
| pathwise gates, one input | 1.00 | 1.00 | 1.01 | 1.00 | 1.03 | 1.03 | 1.04 |

r(Φ_l²) over layers 0…14 is 1.00, 1.62, 1.85, 1.86, 1.91, 1.95, 2.05, 1.95, 1.89, 2.04, 2.07, 2.07, 1.83, 2.05, 2.13: at n = 1024 the mean gates are nearly frozen projections.

## 6. The coordinator's final question: a carrier for the quenched free sector, or a no-go

*Added after coordinator note 1 and the 21:25 cross-team update. That update reports these costate FC dilation numbers, relayed and not re-run here: inside FC, keeping ages ≤ 4 alone gives raw 7.8e-7, adding the dilation sector of older content gives 3.7–4.4e-7, and full FC gives 3.2e-8. So the free sector beyond age 4 is worth ≈ 3.5e-7.*

**6.1 Correction to the premise of note 1 §4 (Measured).** The age spectrum of the content-transfer operator is **not a power law in energy**.
- The free sector loses energy geometrically, at 1.02–1.10 × g_l³ per step (§4.2), where g_l = 2EΦ_l² = 1 − 2E Var(gate). That is ≈ 0.13 per step at layer 0 and ≈ 0.6 per step at depth.
- Only its **rank** follows a power law: PR = n/(2(age + 1)).
- "No spectral gap" (F9.2) does not mean "no mixing". The mixing gap is supplied by the gate variance, not by the spectrum of the propagator. This is the gap between the pathwise wall (χ₁ = 1) and the annealed transfer (g < 1), Proposition C1.
- The only part that does not mix is the dilation sector: its retention rises to 0.95 per step, and its contributions add coherently across ages (§4.1–4.2). Content at the parabolic point is exactly that one sector. Everything else is strictly inside the wall.

**6.2 What a carrier must reproduce at one layer (Measured, `agespec.py` with RANKS=1, networks 0 and 1, n = 1024).** Take the dilation-deflated free residual of ages > 4 and compute the number of singular modes of its D21 slice needed for relative error ε:

| layer | old ages | rank of the sum, ε = 5 / 10 / 20 % | sum of per-age ranks, ε = 10 % | young (≤ 4) sum, ε = 10 % |
|---|---|---|---|---|
| 6 | 1 | 255 / 191 / 126 | 191 | 399 |
| 9 | 4 | 205 / 147 / 90 | 573 | 355 |
| 12 | 7 | 188 / 133 / 80 | 814 | 315 |
| 15 | 10 | 173 / 121 / 72 | 995 | 272 |

(Network 0; network 1 agrees to within ±8 modes.)
- The old free slice needs only **≈ n/8 to n/5 modes at 10 %**, fewer than the age-5 source alone.
- Ranks are **not** additive across ages: the sum needs 121 modes where the per-age ranks add to 995. Freeness would make modes additive at equal energies, but the energies fall geometrically (g³), so the youngest old ages dominate.
- A CP carrier of R atoms yields slices of rank ≤ 3R. The in-layer requirement is therefore only R ≳ 40–60, and **there is no rank no-go at a single layer**.

**6.3 Where the obstruction actually is.** The in-layer readout is cheap. The cost lies in two places:
- **Transport (Theorem, costate C2–C3).** Future slices depend on the content through the second chaos, via the double edge. Any carrier whose state is a function of slices (O(n²) numbers per cut) is O(1)-wrong. Measured: the slice chain and pool renewal lose most of the old content.
- **Merging.** A CP carrier must merge atoms of different ages. The ages sit in mutually free frames (cross-age cosines ≤ 0.03, F9.6), so a *shared* subspace for all ages needs the union of their frames: region's shared Oseledets basis at n/4 is 6× worse. And merging by ALS is itself expensive: region's CP merge to R = n costs 3,510 u and reaches only 2.0e-7.

**The no-go (Conjecture C4, for a class).** Consider carriers that, at every layer, (a) hold the old content as at most R CP atoms in the current layer's coordinates and (b) obtain them by re-fitting a merged approximation of the transported atoms. Assume second-order freeness of different ages (C2(ii)). Then reaching 10 % on the free sector's future slices at every layer needs either R ≥ Σ_{ages} r_{10%}(age, weighted by g³-energy) per merge, or merge work ≥ R² n per layer, with constants that exceed ≈ 7 u per layer at n = 1024.
- The in-layer measurement of 6.2 shows the per-layer requirement itself is mild. What the conjecture puts beyond budget is re-fitting a merged approximation at every layer in a moving frame.
- *Not proved.* What is missing is a lower bound on the merge cost for free low-rank tensors; the time ran out.

**6.4 Escape hatches (what the leaders may be doing instead), in order of how directly 6.1–6.3 point to them.**
1. **Merge once, in a static frame.** Coordinator note §2 writes the memory as a static sum of rank-one birth terms in input space, read through the frame Γ_t. If merging is done there, once per age bin (the note's age multiresolution), the moving-frame re-fit disappears. Two facts make the bins small:
   - 6.2 shows the read-out need is n/8–n/5 modes, dominated by the youngest old bin;
   - 6.1 shows the bin's energy falls by g³ ≈ 0.6 per step, so geometric age bins with resolution set by Π g³ (the T2 rule) are the matched multiresolution, not bins set by PR.
2. **Carry only the trace-weighted young old ages exactly.** T2 at τ = 0.05 already moves the free-sector boundary by energy rather than age. The gain was 0.75×, not 25×, so this alone is not it.
3. **Second-order, not first-order, content.** Costate §4.6: with true κ₃ the first-order model is exhausted at w128, and the remaining 15× is joint κ₄. A leader at raw 1.1–1.5e-8 must carry it, which suggests they treat old content at law level (exact homogeneity, a sampled scale field) rather than as cumulant tensors. The free sector beyond age 4 would then be carried only implicitly.

**Cheapest decisive test for hatch 1.** At n = 1024, take the old free content of ages 5–8 at cut m:
- compress it once in the input-space static frame to R ∈ {n/8, n/4, n/2} atoms (HOSVD of the birth atoms' γ-vectors);
- read it at every later layer through Γ_t;
- report the D21 error of the free sector and the total raw inside FC.

Kill: no R ≤ n/2 reaches 10 % in D21.

## 5. Honest assessment

- **Theorem:** I1–I3 as cited. **Derived:** Proposition C1 (elementary), the polynomial decay C(2k, k)/4ᵏ for the free half-projections, and the parabolic correlation rate.
- **Measured:** §4, first order and decoupled gates, D21 Frobenius metric, n = 1024. **Synthesis:** the identification of the wall at 2 with χ₁ = 1, and the sector map (dilation on the wall, free bulk inside by the gate variance, Perron slightly outside, no elliptic sector). **Conjecture:** C2, C3.
- **Risks.**
  - (a) Everything is in the decoupled-gate first-order model. The dilation retention < 1 is a leak of that model (C6), so the true dynamics may be even more dilation-dominated. That supports the design but means the measured rates are not the truth's.
  - (b) The D21 Frobenius metric is not the readout metric. Costate §4.2 showed the readout weights the dilation spike much more heavily, so the tail prices in §4.4 are conservative for the free sector.
  - (c) The odds constant 0.31 is an empirical regularity. Its value is not derived.
  - (d) None of this addresses joint κ₄, which costate identifies as the remaining 15×. The theory says only that it should again be dominated by the parabolic sector.
  - (e) **Where the design stands.** The theory-driven window is a real but constant-factor gain: 0.75–0.9× raw, about 25 % better adjusted. At ≈ 4.4e-8 adjusted it is still ≈ 27× from the bar (1.6e-9). The first-order model is close to exhausted (costate §4.6), so the next factor must come from second-order content in the parabolic sector.
  - (f) The τ rule uses the leg-wise Frobenius law g³ and ignores readout weighting and the Perron excess (×1.2 per step on the dilation template). A readout-weighted τ is the obvious untested refinement.
