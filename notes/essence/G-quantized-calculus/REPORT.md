# Team G: quantized calculus, Fredholm modules and the age axis — REPORT

*Status: v0 (1 Oct 2026). §1–2 are proofs or exact identities (labelled **thm**/**id**) unless marked **fails** or
**conj**. §3 is measured on the bench (`w256_d16`, `w512_d16`, `w1024_d16`) inside FC (slices = 2, mean-field κ4),
with team D's age-resolution oracle. Code: `fcg.py` (team D's `fcs.py` plus a per-source spectral diagnostic and
arbitrary age profiles), `run_g.py` (T3 sweep), `run_prof.py` (profiles), `analyse_g.py`. Results: `results/`.*

## 0. The answer in brief

1. **(a) and (b) are one statement, and it is exact.** The sign-gauge average of team A, applied to the second-chaos
   operator of a birth, is the conditional expectation E_D of M_n onto its diagonal masa D = ℓ^∞(neurons). With
   F = 2E_D − 1 on L²(M_n, Tr), the quenched memory of an atom is (−½)[F, X]Ω, and the annealed (Bethe) transport is
   the compression E_D ∘ Ad(G) ∘ E_D. **The one-layer cut is a commuting square iff every column of G = diag(P)W has at
   most one non-zero entry**, i.e. iff G is a Bratteli (AF, tail-equivalence) connecting matrix. Gaussian G is at the
   opposite extreme. The Hochschild defect of the compression is exactly ¼ P[F,a][F,b]P (Connes' identity).
2. **Note 2's "not Fredholm" is the wrong diagnosis.** [F, λ(X)] has rank ≤ 2n in an n²-dimensional space, so it is
   asymptotically compact (relative rank 2/n → 0) for every n. What does *not* vanish is its value on the state
   vector, ‖[F, X]Ω‖/‖X‖ ≈ 0.9–1 per atom and 0.44–0.48 for the summed old memory (team A). The obstruction is size on
   the readout subspace, not non-compactness.
3. **(c) holds but is Pythagoras.** The books-close identity is the orthogonality of the spectral projections of any
   Christensen–Ivan operator Σ λ_l(E_l − E_{l−1}) on the weight filtration. The weights K(l) are a gain law, not the
   spectrum of an operator acting on the error.
4. **(d) is exact for degree 1 and fails as stated for degree 2.** The memory H^c_t = E[∇²z_{t,c}] is the
   Gaussian-weighted sum over upstream walls of the quantized differentials [F, ∇z_{t,c}] of the wall module (the
   distributional Hessian of a piecewise-linear function is its jump measure). But every cyclic cocycle of a module
   whose F swaps the two sides of each wall separately is a sum over *single* walls. Pairs of walls at different
   depths are not a cyclic pairing of that module.
5. **(e) holds in the right variables.** The age-resolution operator R = diag(k(a)/n) is in L^{1,∞} and its Dixmier
   trace is c exactly when k(a) = c·n/a. The aging step is the unilateral shift on ℓ²(ages); with the Toeplitz
   polarization, Index(P u P) = −1 is the one birth per layer. The age odometer is the 2-adic adding machine on the
   birth times; a sound version merges blocks with a delay (§2.1).
6. **T3 (measured, widths 256 and 512; 1024 in flight).**
   - The *relative* truncation loss δ(c)/δ_drop is width-universal: ≈ 0.03 at c = 1, ≈ 0.003 at c = 1.5, ≈ 3·10⁻⁴ at
     c = 2, at both n = 256 (8 networks) and n = 512 (4). It decays like e^{−4.6c}.
   - The participation ratio of an age-a propagator is n/(2a) at both widths (a·PR/n = 0.50 ± 0.03 for a = 3–10).
   - But the c needed for "lossless" (loss ≤ 10 % of FC's raw) drifts up with width: c* ≈ 1.24 (n = 256), 1.39 (512),
     ≈ 1.58 (1024, team D's points). FC's own error falls faster with n than the dropped memory does, so the ratio
     to be resolved grows by ≈ e^{0.8} per doubling.
   - Dixmier coefficient at t = L = 16 (young ages 1–2 at full rank): (2 + 1.881c*)/ln 16 = 1.56 / 1.66 / 1.79.
7. **T4.** Conditional theorem (§4.1): under exponential spectral tails with rate ∝ (a + 1)/n, additive age losses
   and a 1/a energy profile, Σ_a k(a) ≥ (n/β) ln(1/ε)(H_{L+1} − 1), sharp. **Converse (§4.2): the Dixmier law is
   forced only by the 1/a energy profile.** With geometric or faster-than-1/a energy decay, the optimal total is O(n)
   in L, and the per-source data (a·k_ε(a)/n falls linearly in a, energy share falls ≈ geometrically) say we are in the
   subcritical regime. **The cost form Ω(n³ L log L) is false asymptotically** (§4.3): shared dyadic frames plus exact
   Tucker cores make cost independent of the number of ages per block. At L = 16, n = 1024 neither escape saves
   much: Σ_{a=3}^{16} 1/a = 1.88.

## 1. The dictionary made rigorous

Notation: n width; G_s = diag(P_s) W_{s+1} the mean-gate step; Z_s(t) = G_s ⋯ G_{t−1}; for a birth (s, r) the leg
y_r = row r of Y_s(t) ∈ ℝⁿ (target index), so the source's D21 contribution is Σ_r w_r [(y_r∘y_r) z_rᵀ + 2 …]
(coordinator note 1, §2). Second chaos H₂ ≅ Sym(n) ⊂ M_n with the Hilbert–Schmidt inner product.

### 1.1 Entries (a)+(b): the gauge polarization is the diagonal masa, and the cut is a failed commuting square

**Definition G1.** D ⊂ M_n the diagonal masa, E_D(X) = diag(X₁₁, …, X_nn) its trace-preserving conditional
expectation, H = L²(M_n, Tr) with cyclic vector Ω = 1, λ(X) left multiplication, F = 2E_D − 1 (a self-adjoint unitary
on H, F² = 1).

**Proposition G1 (id).**
1. E_D(X) = ∫ S X S dS over the sign group S ∈ {±1}ⁿ (diagonal), so team A's (Z₂)ⁿ gauge average acting on second-chaos
   operators is E_D. The same holds for the hyperoctahedral group; the O(n) average is the smaller projection
   X ↦ (tr X/n)·1 (team A's trace channel).
2. X − E_D(X) = −∫ [S, X] S dS: the off-diagonal (quenched) part is an average of commutators with gauge unitaries.
3. [F, λ(X)]Ω = −2(X − E_D X), so ‖[F, λ(X)]Ω‖₂ = 2‖X − E_D X‖₂.
4. rank [F, λ(X)] ≤ 2n (since rank E_D = n on H), so the normalised trace τ = Tr/n² of |[F, λ(X)]|^p is ≤ (2/n)‖X‖^p → 0:
   **the module (H, F) over λ(M_n) is p-summable for every p, uniformly, in the relative (Breuer) sense.**

*Proof.* 1: averaging S_i S_j over signs gives δ_ij. 2: X − SXS = [X, S]S… summed. 3: F(XΩ) − X F Ω = (2E_D X − X) − X.
4: [E_D, λ(X)] = E_D λ(X) − λ(X) E_D has range in ran E_D + λ(X) ran E_D. ∎

*Reading.* For one atom X = y yᵀ, ‖E_D X‖₂/‖X‖₂ = ‖y‖₄²/‖y‖₂² ≈ √(3/n): the quenched part is essentially all of the
atom. Coherent summation over atoms puts weight on the diagonal; team A measured ‖D_old − Bethe‖/‖D_old‖ = 0.44–0.48 at
n = 1024, flat in width. So note 2's sentence "[F, a] is not compact … our system is not a Fredholm module" is
**incorrect as stated**: the commutator is always compact (finite rank, relative rank 2/n). The obstruction is that its
value on the one vector that matters is O(1). The right statement is (G1.3) plus the measured 0.44.

**Proposition G2 (Connes' identity for compressions, id).** For a projection P and F = 2P − 1,
P a P b P − P ab P = ¼ P[F, a][F, b]P. *Proof:* P[P,a][P,b]P = −P a(1−P) b P. ∎ So a compression a ↦ PaP is a
homomorphism exactly when the quantized differentials vanish on ran P, and its Hochschild coboundary is the
degree-2 quantized product.

**Proposition G3 (the cut as a commuting square, thm).** Let Φ_G(X) = Gᵀ X G on M_n (the exact transport of a
second-chaos operator through one mean-gate step; it is how y yᵀ is carried). The annealed (Bethe) transport is
E_D Φ_G E_D, i.e. x ↦ x(G∘G) on diagonals. The square (D ⊂ M_n, Φ_G) commutes, E_D Φ_G = E_D Φ_G E_D, **iff every
column of G has at most one non-zero entry.**

*Proof.* (GᵀXG)_jj = Σ_k G_kj² X_kk + Σ_{k≠l} G_kj G_lj X_kl. The second sum vanishes for all off-diagonal X iff
G_kj G_lj = 0 for all k ≠ l and all j. ∎

- Such G are exactly the connecting matrices of an AF (Bratteli) system of commutative algebras, the maps of
  Christensen–Ivan's and Pearson–Bellissard's settings (a tree: each child has one parent). On a tail-equivalence
  Cantor set the masa is preserved by every step, so the conditional expectations commute with the dynamics.
- For He-Gaussian G the off-diagonal-to-diagonal leak per unit norm is Σ_j Σ_{k≠l} G_kj²G_lj² ≈ 2g²/n per entry
  but it acts on the n² − n off-diagonal entries, which carry ≈ all of an atom's norm. **Our propagators are the
  maximally non-AF case.** This is the precise reason the Cantor-set / AF machinery of note 2 §1 does not apply
  verbatim: it needs the commuting square that Gaussian transport breaks at every cut (this is C2/C3 of costate, now
  as a commuting-square statement).
- Entry (b) of note 2 is G3 with "Hadamard product" read as the product of D: the homomorphic part (G∘G) is the
  compression, the off-diagonal terms are the Hochschild defect, G2.

### 1.2 Entry (c): the books-close identity and Christensen–Ivan

**Proposition G4 (id, weak).** Let F_l = σ(W_1, …, W_l), E_l = E[·|F_l] on L²(Ω_W), Δ_l = E_l − E_{l−1}. For any
real λ_l, D_λ = Σ_l λ_l Δ_l is self-adjoint with spectral projections Δ_l, and for any square-integrable error e,
‖e‖² = Σ_l ‖Δ_l e‖². If one *defines* the local error ε_l by Δ_l e = T_{l→L} ε_l with gain ‖T_{l→L}‖² = K(l), the
books-close identity MSE = Σ_l K(l)‖ε_l‖² is ⟨ε, D_K² ε⟩ with D_K = Σ_l K(l)^{1/2} Δ_l.

**Where it fails as a spectral triple.** Christensen–Ivan need finite-dimensional A_k and λ_k → ∞ so that D has
compact resolvent; here each Δ_l L² is infinite-dimensional and L is fixed. The content of (c) is entirely in the
gain law K(l) (team E) and in locality of ε_l, not in the operator. **Status: holds, but carries no extra information.**

### 1.3 Entry (d): the wall module

**Definition G5.** Let 𝒲 be the set of facets of the activation fan of the network (on ℝ^{n₀} or the sphere). For each
facet f with sides R_f^±, H_𝒲 = ⊕_f ℂ², grading Γ = ⊕ diag(1, −1), F = ⊕ σ_x (swap), and a piecewise-continuous
function a acts by π(a) = ⊕ diag(a(x_f⁺), a(x_f⁻)) (one-sided limits at a point of f). Weighted trace
Tr_γ(T) = Σ_f γ(f) tr T_f with γ(f) the Gaussian surface measure of f.

**Proposition G5 (thm).**
1. [F, π(a)]_f = (a(x_f⁻) − a(x_f⁺)) J with J = [[0, 1], [−1, 0]]: singular values are the jumps. For a continuous
   piecewise-linear z, [F, π(z)] = 0 and [F, π(∇z)]_f = [∇z]_f, the kink.
2. **(Degree 1, the memory.)** For z = z_{t,c}, E_x[∇²z] = Σ_f ∫_f [∇z]_f ⊗ ν_f dγ_f (the distributional Hessian of a
   PL function is its jump measure on walls). The jump of ∇z_{t,c} across the wall of upstream neuron (s, r) is
   ∂z_{t,c}/∂a_{s,r} · ∇z_{s,r} (one-sided), whose mean-gate Gaussian average is Z_s(t)_{rc} γ_{s,r}. So
   H^c_t (coordinator note 1, §2) is the Γ-free degree-1 functional a ↦ Tr_γ(ν ⊗ |[F, π(∇a)]|·sign) of this module,
   summed over the walls of all earlier layers. **The memory is a weighted trace of quantized differentials over
   upstream walls.**
3. **(Degree 2 fails.)** The even character τ(a₀, a₁, a₂) = Tr_γ(Γ a₀[F, a₁][F, a₂]) equals
   −Σ_f γ(f) Δ_f a₀ Δ_f a₁ Δ_f a₂ (all jumps at the same wall). Every cyclic cocycle of a block-diagonal F localises on
   single walls, so "pairs of walls at different depths" is not a pairing of this module. The two-wall object (birth
   wall at depth s, readout wall at depth t) is the product of a degree-1 functional 2 with the readout's own
   δ′(z_{t,c}); it is a cumulant, not a cocycle. To make it a cocycle one needs an F that couples walls, i.e. a
   non-local polarization. **Status: (d) degree 1 thm; degree 2 fails as stated.**

**Summability of the wall module at He-criticality (thm, scaling).** For a fixed readout c, the singular values of
[F, π(∇z_{t,c})] weighted by γ are |Z_s(t)_{rc}|·‖γ_{s,r}‖ over the nL upstream walls (s, r), each O(n^{−1/2}). So
Tr_γ|[F, ∇z]|^p ≍ L·n·n^{−p/2}: **uniformly 2-summable in width** (p = 2 is the threshold), with the depth dependence
carried by the per-age energy Σ_r w2 Z_s(t)_{rc}² (measured in §3.3: falls ≈ geometrically in age).

### 1.4 Entry (e): the age axis is a dimension-1 spectral triple; the shift has index −1

**Definition G6.** H_age = ℓ²({1, 2, …}) ⊗ ℂⁿ, the age-resolution operator R = Σ_a (k(a)/n) e_a e_aᵀ ⊗ 1, and the age
Dirac operator D_age = N (the number operator a ↦ a).

**Proposition G6 (thm).**
1. If k(a) = c·n/a then R = c|D_age|^{−1} ∈ L^{1,∞} \ L¹, Tr_ω(R ⊗ (1/n)) = c, and the total per-target resolution is
   Σ_{a≤L} k(a) = c n H_L ≈ c n ln L. The "Dixmier coefficient" Σ_a k(a)/(n ln L) → c is the Dixmier trace.
2. (D_age, ℓ^∞) has spectral dimension 1. If instead a·k(a)/n → 0, R ∈ L¹ and the total resolution is O(n) in L.
3. The aging step u: e_a ↦ e_{a+1} is the unilateral shift; on ℓ²(ℤ) with the Toeplitz polarization F = sign(N)
   (P = projection on ages ≥ 1), [F, u] has rank one and Index(PuP) = −1. **The global integer is the one new source
   born per layer**; the only non-compact datum is the boundary age. (This is the honest version of the user's "F as
   a shift": the shift is the algebra element, F is the Hardy polarization of the age axis.)

### 1.5 The Fredholm module on the dyadic age tree (Pearson–Bellissard), and the odometer

- **Space.** Birth times s ∈ ℕ embed densely in the 2-adic integers ℤ₂; the time step s ↦ s + 1 is the odometer
  (adding machine), a minimal isometry of ℤ₂ whose orbit relation is tail equivalence of binary digit sequences (up to
  the orbit of −1). This is the exact sense in which "aging is a shift whose orbits are tail classes".
- **Module (Pearson–Bellissard).** Vertices v of the dyadic tree on birth times = blocks B_{j,m} = [m 2^j, (m+1)2^j).
  At each v choose ξ_v^± in its two children; H_PB = ⊕_v ℂ², π(f) = ⊕ diag(f(ξ_v⁺), f(ξ_v⁻)), F = ⊕ σ_x (swap the
  sibling blocks), D_PB = ⊕_v λ_v F_v.
  - With λ_v = 2^j (inverse 2-adic diameter) this is the standard triple of ℤ₂: ζ(s) = Tr|D_PB|^{−s} = Σ_j 2·(2^{−j})^{s}·#{v at level j}
    converges for s > 1 on any finite window per level; spectral dimension 1 (= the 2-adic Hausdorff dimension).
  - With λ_v = 2^j/c (inverse resolution share k_j/n of a block whose members have age ≈ 2^j) the same module carries
    the age-resolution law; its |D|^{−1} is R of G6 restricted to block representatives.
  - [F_v, π(f)] = f(ξ_v⁺) − f(ξ_v⁻) is the jump of the per-age carrier between sibling blocks; for the resolution law
    f(a) = c/a the jump at level j and position a is ≈ c 2^{j−1}/a²: the old siblings are nearly indistinguishable,
    which is why sharing a frame across a dyadic block is plausible (and what team D's odometer test checks).
- **The odometer carrier, made sound.** Bentley–Saxe binary-counter blocks (bit j of t ⇒ a block of 2^j consecutive
  birth times) contain young members whenever the low bits of t vanish (t = 2^J gives one block of ages 1 … 2^J), so
  they cannot be read at a resolution set by age. The sound version merges two sibling level-j blocks into their
  parent only when the parent's youngest member reaches age 2^{j+1}. Then every level-j block has ages in
  [2^j, 3·2^j), at most two blocks live per level, and merges at level j fire with period 2^{j+1}: the carry pattern of
  the odometer, delayed by one block.
- **Merging is exact only in a shared frame (id).** If the members of a block share a target frame U (n × k):
  Y_s = A_s Uᵀ, Z_s = B_s Uᵀ, then the block's D21 is D_ab = Σ U_ai U_ai′ U_bj C_{ii′j} with core
  C = Σ_{s∈block} Σ_r w_{sr} a_{sr} ⊗ a_{sr} ⊗ b_{sr}, summable across members; transport is exact
  (Gᵀ U = U′R ⇒ C ← C ×₁ R ×₂ R ×₃ R); readout costs n k³ + n² k per block. (Sent to team D.)

## 2. Note on what the dictionary buys

- G1–G3 identify the campaign's obstruction with a single, classical failure: the masa of neurons is not preserved by
  Gaussian transport, so no AF / tail-equivalence (commuting-square) structure exists on the neuron side. The only
  place a genuine AF structure exists is the **age axis** (G6, §1.5), which is where the user's Cantor-set intuition
  lands correctly.
- The index statements are integers (one birth per layer); as note 2 §4 says, they do not evaluate quenched means.

## 3. Tests run (T3)

### 3.1 Score-level sweep of the uniform law k(a) = c·n/a (ages ≥ 3; ages 1–2 exact)

δ(c) = mean final-layer (p_c − p_FC)², noise-free; δ_drop = the same for dropping every age > 2.

| width (nets) | FC raw | δ_drop/FC raw | δ(1)/δ_drop | δ(1.5)/δ_drop | δ(2)/δ_drop | loss/FC raw at c = 1 / 1.5 / 2 |
|---|---|---|---|---|---|---|
| 256 (8) | 1.8–4.2e-6 | 4–14 | 0.028–0.062 | 0.002–0.008 | 1–6 e-4 | 0.22–0.39 / 0.020–0.046 / 0.002–0.006 |
| 512 (4) | 2.5–3.7e-7 | 24–26 | 0.021 | 0.002–0.003 | 2–4 e-4 | 0.47–0.51 / 0.058–0.070 / 0.005–0.008 |
| 1024 (team D, m0) | 3.24e-8 | (in flight) | — | — | — | ≈1.5 / ≈0.15 / ≈0.02 (raw differences) |

- **Universal:** δ(c)/δ_drop ≈ A e^{−4.6c} at both widths (per Δc = 0.5 a factor ≈ 10).
- **Not universal:** the dropped memory relative to FC's own error grows ≈ 3–5× per doubling of n, so the c needed for
  fixed *relative-to-FC* accuracy drifts: c*(10 %) ≈ 1.24, 1.39, 1.58 at n = 256, 512, 1024.
- **Dixmier coefficient per target at t = 16:** (2 + 1.881 c*)/ln 16 = 1.56 / 1.66 / 1.79 (with c = 2: 2.08).

### 3.2 Per-source spectral diagnostic (pure diagnostic, exact contributions downstream)

For each source of age a ≥ 3 at each target, the relative error of its D21 contribution when all legs are projected on
the top-k right singular vectors of Z_s(t), on a grid k/n ∈ {1/64, …, 3/4}; k_ε(a) = k at relative error ε.

| n, MLP | a·PR(Z)/n at a = 3 / 6 / 10 / 15 | a·k_{0.01}(a)/n at a = 3 / 6 / 9 / 12 | energy share η_a at a = 3 / 6 / 9 / 12 |
|---|---|---|---|
| 256, m0 | 0.51 / 0.50 / 0.50 / 0.47 | 1.27 / 1.09 / 0.62 / 0.20 | 0.067 / 0.020 / 0.007 / 0.002 |
| 256, m1 | 0.50 / 0.46 / 0.38 / 0.28 | 1.31 / 1.24 / 0.91 / 0.42 | 0.053 / 0.018 / 0.009 / 0.006 |
| 512, m0 | 0.50 / 0.49 / 0.44 / 0.34 | 1.19 / 0.97 / 0.57 / 0.27 | 0.053 / 0.019 / 0.009 / 0.005 |
| 512, m1 | 0.50 / 0.49 / 0.43 / 0.33 | 1.31 / 1.22 / 0.83 / 0.41 | 0.058 / 0.020 / 0.010 / 0.006 |

(Grid floor k/n = 1/64 makes a·k/n ≥ a/64 at the oldest ages.)

- **PR(Z_s(t)) = n/(2a)** at both widths (team C's law with a in place of a + 1).
- **Per-source resolution falls faster than 1/a:** a·k_ε(a)/n decreases roughly linearly in a. The tail exponent of
  the per-source error, e(k) ≈ exp(−E(a) k/n), has E ≈ 7.4, 20, 35 at a = 3, 7, 10 (n = 256, m0): superlinear in a.
- **Energy share decays ≈ geometrically** (×0.65–0.75 per age to a ≈ 9, slower beyond), consistent with team C's
  g³ law for the free sector.

## 4. T4: the Dixmier lower bound

### 4.1 Conditional theorem (thm)

*Setting.* At target t the old memory is a sum of age blocks B_a (a = 3 … L). A *linear-frame carrier* reads block a
through orthogonal projections of rank k_a applied to its target legs (any frames: transported, shared, oracle). Let
τ_a(k) = min over rank-k projections of the relative loss of B_a. By the mode-unfolding argument (projecting more
legs only adds orthogonal error, and Eckart–Young on one leg), τ_a(k) ≥ Σ_{j>k} σ_j²(B_a^{(1)})/‖B_a‖².

*Hypotheses.*
- (H1, spectral tail) τ_a(k) ≥ exp(−β(a + 1)k/n) for 0 ≤ k ≤ n (exponential tail with rate linear in age; the
  Lyapunov / Fuss–Catalan heuristic: σ_j² ≈ g^a(1 − j/n)^a gives tail (1 − k/n)^{a+1} on the propagator).
- (H2, additivity) the losses at different ages add: total relative loss = Σ_a η_a τ_a(k_a), η_a = ‖B_a‖²/Σ‖B‖²
  (age near-orthogonality F8.3 and the fresh-weight lemma).
- (H3, critical energy) η_a = 1/((a + 1)(H_{L+1} − 1)) (energy share ∝ 1/(a + 1)).

*Theorem G7.* Any carrier with total relative loss ≤ ε satisfies Σ_a k_a ≥ (n/β) ln(1/ε) (H_{L+1} − 1), and the bound
is attained by k_a = n ln(1/ε)/(β(a + 1)).

*Proof.* Write S = Σ_a 1/(a + 1) (= H_{L+1} − 1 when every age 1 … L is carried) and b_a = β(a + 1)/n. Lagrange weak
duality: if Σ_a η_a e^{−b_a k_a} ≤ ε then for every μ ≥ 0, Σ k_a ≥ Σ_a min_{k≥0}[k + μ η_a e^{−b_a k}] − μ ε. Under H3,
θ := μ η_a b_a = μβ/(nS) is the same for every a; for θ > 1 each minimum is attained at e^{−b_a k} = 1/θ with value
(1 + ln θ)/b_a. So Σ k_a ≥ (nS/β)(1 + ln θ − θε) for every θ > 1; θ = 1/ε gives (nS/β) ln(1/ε). Attainment: with
k_a = n ln(1/ε)/(β(a + 1)) each age contributes η_a ε, total ε. ∎

### 4.2 Converse: the Dixmier law is forced by H3 and nothing weaker (thm)

Keep H1–H2 and let η_a ≤ C a^{−1−δ} (δ > 0) or η_a ≤ C q^a. The optimal allocation is
k_a = (n/(β(a+1))) ln(μ η_a b_a)₊ (μ fixed by the constraint), which vanishes beyond an age a_max(ε) independent of L
(a_max ≍ (C/ε)^{1/δ} or log_{1/q}(C/ε)). Hence min Σ k_a ≤ (n/β) H_{a_max} ln(μ max_a η_a b_a) = O_ε(n), **bounded in L**. With a flat
profile η_a = 1/L the optimum is ≈ (n/(2β)) ln²(H_L/ε): still o(n ln L). So the n ln L law needs exactly the 1/a energy
profile (or a tail rate that does not grow with age).

**What the data say (§3.2).** η_a falls ≈ geometrically, and the per-source tail rate grows faster than linearly in a.
Both push to the subcritical side of G7/G8: **the measured c·n/a law is an over-resolution of old ages, not a
Dixmier-critical necessity.** The uniform-c law works because it equalises per-age loss, not because it is optimal.
Profiles that resolve old ages less (k = n max(2/a − γ, 1/64), k = c n/a^p with p = 1.5, 2) are being run at
n = 256, 512, 1024 (`run_prof.py`; 22–37 % fewer old atoms than c = 2).

### 4.3 The cost form is false asymptotically (thm, construction)

With shared frames per dyadic block (§1.5) and exact cores, block j costs min(c n³, n k_j³ + n² k_j) per layer with
k_j = c n/2^j, so the total per layer is ≈ c n³ min(log₂ L, log₂(c^{2/3} n^{1/3})) + O(n³): Θ(n³ log n), not
n³ log L, once L ≫ n^{1/3}. Total *resolution counted with multiplicity* is still c n ln L, so the resolution form of
note 2's conjecture can hold while the cost form fails. At n = 1024, L = 16 the crossover age 2√n = 64 exceeds L: no
gain for the competition.

## 5. Honest assessment

- **Theorems:** G1–G3 (exact; the commuting-square characterisation is the cleanest new statement), G5 degree 1,
  G6, G7 (conditional) and its converse, the cost construction. **Failures recorded:** note 2's "not Fredholm"
  (wrong diagnosis), (c) carries no information beyond Pythagoras, (d) at degree 2.
- **Measured:** T3 at n = 256 (8 nets) and 512 (4); n = 1024 in flight. The universal object is the relative loss
  law δ(c)/δ_drop ≈ A e^{−4.6c}; k(a)/n at fixed relative-to-FC accuracy is *not* width-universal (drifts ≈ +0.17 in c
  per doubling).
- **Risk:** the per-source diagnostic measures D21-level loss, not score-level; the score-level profile runs decide.
  Widths ≤ 512 are pre-asymptotic per the brief.
