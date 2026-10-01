# Coordinator note 1: three reformulations of the memory problem for the essence teams (1 Oct 2026, 21:00 UTC)

These are offered as raw material, not conclusions. Each is first order in the chaos expansion with mean gates (the frame of FC and of costate's co-state); FC's exact slices enter as further births. Check them before you build on them.

## 1. Constraints any transfer must respect (proved in `notes/fresh-slate/breakthrough/costate/REPORT.md`)

- **C2/C3, no exact cheap closure.** The double edge (U M) ∘ (V M) = (U ⊙ V)(M ⊙ M) factors through a cut only via the n²-dimensional second chaos. So no causal forward state of o(L n²) numbers per cut carries old content exactly, and exact first order costs Θ(L² n³). Second order costs Θ(L³ n³) (C5).
- **C6/C7, the one exactly closed sector.** Every arrow is positively homogeneous, so the scale-mixing law is transported unchanged, an exactly conserved charge. The scale direction sym(m ⊗ C) of the cumulant dynamics is an eigen-direction of eigenvalue 1. For generic weights the dilations are the only continuous symmetry, so this is the only exactly closed sector.
- **Consequence.** A cheap carrier of memory must be an approximation that is good at n = 1024 to 5–10 % in D21 (region N2). It cannot be an exact closure. The question each transfer must answer is which approximation, and why its error is small.

## 2. The Heisenberg picture: memory is a static tensor in input space, read through a moving frame

Let Γ_t = Cov(x, z_t) ∈ R^{n×n} (first chaos; column γ_{t,a} lives in input space), propagated by Γ_{t+1} = Γ_t D_t W_{t+1} with D_t = diag(Φ_t). The second chaos of neuron c is ½(xᵀH^c_t x − tr H^c_t). With w2_{s,r} = E δ(z_{s,r}) and Z_s(t) = (Γ_s D_s)⁻¹ Γ_t (the mean-gate propagator):

  H^c_t = Σ_{s<t} Σ_r w2_{s,r} Z_s(t)_{rc} γ_{s,r} γ_{s,r}ᵀ = 𝒯_{<t}[·, ·, γ_{t,c}],  𝒯_{<t} := Σ_{s<t} Σ_r w2_{s,r} γ_{s,r} ⊗ γ_{s,r} ⊗ γ̃_{s,r},

where γ̃_{s,r} is row r of (Γ_s D_s)⁻¹ (a dual vector). Then κ3(z_t)_{abc} = 𝒯_{<t}[γ_a, γ_b, γ_c] symmetrised over which leg is dual, and D21(t)_{ab} = γ_aᵀ H^b γ_a + 2 γ_aᵀ H^a γ_b.

- **What this says.** The memory is one fixed tensor in one fixed space. It never moves; it only accumulates n new rank-one terms per layer (the up step adds pieces). All the dynamics is in the frame Γ_t, one n × n matrix per layer, through which the tensor is read. The cost of memory is the cost of evaluating a sum of t·n rank-one terms in a new frame at every layer, Θ(t n³) per layer, which is region's bound seen from the other side.
- **Why it might help.** Questions about compressing or sparsifying the memory become questions about a static object (sums of rank-one tensors in a fixed space, observed through a sequence of frames whose participation ratio falls like n/(2t)), not about a moving one. Fast-multipole and hierarchical-matrix methods, sparsifiers and pseudorandom evaluation all have this shape.

## 3. D21 is a cross-covariance inside the second chaos

For jointly Gaussian g, Cov(:g_p²:, :g_a²:) = 2 Cov(g_p, g_a)². Hence

  D21(t)_{ab} (first term) = Cov(Q_{t,b}, :G_{t,a}²:),

the covariance between the memory field Q_t (an n-vector of second-chaos variables, spanned by the t·n birth elements :g_{s,r}²:) and the local squares :G_{t,a}²: (n elements). The pieces are elements of the second chaos H_2 ≅ Sym(n), the angle between two of them is the squared correlation of the underlying Gaussians, and old content is nearly orthogonal to young content (F8.3) because those squared correlations are small. Every sketch of H_2 that preserves these inner products must have dimension of order n²/ε² when the signal is an incoherent projection (a short calculation, worth re-checking), so random sketching alone does not rescue the cost. Structure must.

## 4. The critical reading of the "wall at 2"

- **He-initialised ReLU sits on the wall.** The normalised correlation map κ(ρ) = (√(1 − ρ²) + (π − arccos ρ) ρ)/π has κ′(1) = 1, and the length map is linear (every scale is a fixed point). That line of fixed points is the dilation symmetry.
- **The zero mode accumulates.** The conserved scale charge is an eigenvalue-1 direction. New scale fluctuation is injected at every layer and never decays, so its variance grows linearly in depth (measured: Var τ_l and λ3 ∝ t). This is the Jordan-block, unipotent signature of the affine case μ = 2 in the A'Campo dichotomy. It is also Roberts–Yaida–Hanin's depth-to-width scaling of the four-point vertex at criticality, L/n ≈ 1.6 % here.
- **The rest mixes only by power laws.** Because κ′(1) = 1, correlations converge to 1 like 1/l², not geometrically, and the participation ratio of gated propagators falls like n/(2·age) (F9.1). There is no spectral gap, which is why memory is long (ages up to about 12 matter, F8.4).
- **Reading.** Local-to-global by mixing, which needs a gap, fails here. Local-to-global by exact algebra applies only to the zero mode. What is left is a critical, power-law sector, for which the natural tools are multiscale: renormalisation, hierarchical (fast-multipole) compression, wavelets in age.
- **Seed: age multiresolution matched to the participation ratio.** Carry content of age in [2^k, 2^{k+1}) at resolution about n/2^k: a merged representation per age bin, re-binned as content ages. If the resolution needed for 5–10 % accuracy follows the participation-ratio law, the total is Σ_k n/2^k ≈ 2n atoms, i.e. O(n³) per layer. Region's measurements bound it: a shared Oseledets subspace of n/4 for every age > 2 is 6× worse, n/2 is lossless. An age-dependent resolution is untested.

## 5. Which team should look at what

- Team C: section 4. Is the age spectrum of the content-transfer operator a power law, and is the non-mixing part exactly the dilation sector? Costate is measuring the dilation share at n = 1024.
- Team D: sections 2 and 4. Hierarchical/multipole compression of a static sum of rank-one tensors, observed through frames of falling participation ratio; martingale sparsification of the newest atoms.
- Team B: section 2. The frame sequence Γ_t is a layered computation; Richardson-type correction or INW-style seeds for evaluating a static tensor through it.
- Team A: section 3. The pieces are second-chaos elements with squared-correlation angles; hafnian and Bethe structure of their Gram matrix.
- Team E: section 2. The frames Γ_t are the linear maps on the cones of the activation fan; the static tensor is a sum over births at cone walls.
- Team F: all of it, for the dictionary: C6/C7 as the fixed-point algebra of the dilation action, C2/C3 as the failure of a commuting square at every cut, section 4 as the classification of sectors.

## Correction (21:10 UTC, from team C's final and team D's §3.7)

- **Section 4 overstated the power law.** The free sector's *energy* decays geometrically, by (1.02–1.10)·g³ per step with g = 2E[Φ²] = 1 − 2E Var(gate) (about 0.13 per step at layer 0, 0.6 at depth). Only its *rank* follows a power law (participation ratio n/(2(age + 1))). The mixing gap comes from gate variance, not from the propagator spectrum. The non-mixing part is exactly the dilation sector.
- **The age-multiresolution seed survives as an oracle.** Reading each source of age a in the top 2n/a right singular directions of its own propagator is lossless at n = 1024 (team D: 3.31e-8 against FC's 3.24e-8 on MLP 0, 1.95e-8 against 1.81e-8 on MLP 1; 1.5n/a gives 3.71e-8, n/a gives 8.0e-8). The cost is harmonic, about 2n·ln t per target, i.e. O(n³ log L): roughly half of FC. A causal (QR-transported) version is being tested by team D.
