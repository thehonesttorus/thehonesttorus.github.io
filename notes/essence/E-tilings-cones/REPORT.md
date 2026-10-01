# Team E: tilings, aperiodic order and conic geometry

*v0, 1 Oct 2026, late evening (UTC). Status labels: **Theorem** (published result, source given), **Derived** (proved here or a short derivation given), **Measured** (number from a script in this directory or a cited stream file), **Synthesis**, **Speculation**. arXiv ids marked "(verify)" were quoted from memory and not opened in this session.*

## 0. Result in one paragraph

The strongest thing this team found is a **conic "wedge calculus"**. It is a parameter-free, annealed model that reproduces three established facts of the campaign at n = 1024 from 2-D conic geometry alone:
- the region stream's measured transfer coefficients K_mean(l) and K_off(l), all 15 layers, to 11 % rms in log, which is within the MLP-to-MLP scatter;
- region's N6 ("dropping ages > 4 costs 25× the bar"): predicted relative D21 error 0.72 when ages ≤ 4 are kept and 0.86 when ages ≤ 2 are kept, against measured 0.75 and 0.89;
- the faces stream's flat age profile (T3).

The mechanism is one conic quantity. The Gaussian measure of the wedge in which one neuron is on for both of two independent inputs is 1/4 + arcsin(ρ_l)/(2π) (Sheppard), where ρ_l is the replica correlation on the arc-cosine orbit. Given W this wedge measure equals the neuron average of Φ_j², so it is the per-leg energy factor of every gated transport. ReLU at He scale is critical: ρ_l → 1 and the wedge factor tends to 1 − 3/l. Memory is therefore a power law, not geometric. Meanwhile the data cone narrows, so the birth strength (facet density on the cone) falls by about 2.5× per early layer. The product of the two is nearly flat in age, and that is why every age matters.

This *explains* the old-content facts. It does not by itself make them cheaper. The model is annealed and the binding object is quenched (F8.5). That is the general verdict for this domain (§2.4): tilings and cones turn local data into global truth when an ergodic or group average (patch frequencies, intrinsic volumes, Perron vectors) is the target. Our target is the quenched remainder after that average, so the transfer that survives is the *pricing* and *budget* side, not a cheaper carrier.

## 1. Instances

### 1.1 Dimers: Kasteleyn, height functions with coupling from the past, and the Cohn–Kenyon–Propp variational principle

- **Exact counting (Theorem).** Kasteleyn (1961) and Temperley–Fisher: for a planar graph, a Pfaffian orientation (every face has an odd number of clockwise edges) makes Pf(K) equal the number of perfect matchings. A #P-complete count becomes P in general, but only on planar graphs.
  - **Local-to-global content.** A *local* sign condition, one per face, makes a *global* signed sum cancel exactly to an unsigned one.
  - **The boundary of the method.** Valiant: counting matchings is #P-hard in general. A genus-g graph needs 4^g Pfaffians (Galluccio–Loebl; Tesler).
- **Sampling (Theorem).** The flip chain on lozenge tilings is the down–up skeleton:
  - **Pieces:** hexagons.
  - **Down step:** remove a hexagon of three rhombi.
  - **Up step:** re-insert it rotated.

  Height functions make the chain monotone, so coupling from the past (Propp–Wilson 1996) gives exact samples from the top and bottom states alone. Wilson (2004, math/0102193 (verify)) proves mixing in O(n³ log n) for lozenge tilings of an n×n region. The certificate is a contraction of the height-function lattice under the monotone coupling, i.e. an eigenfunction (a discrete Laplacian of the heights).
- **The variational principle (Theorem; Cohn–Kenyon–Propp, JAMS 2001, math/0008220 (verify)).**
  - **Statement.** The scaled height function of a uniform random tiling of a large region concentrates on the unique maximiser of ∫ ent(∇h). The local surface tension ent(s) is the entropy per site of the translation-invariant ergodic Gibbs measure of slope s.
  - **Explicit forms.** Kenyon–Okounkov–Sheffield (math-ph/0311005 (verify)) give ent for every bipartite periodic dimer model through the spectral curve P(z, w). Kenyon–Okounkov (math-ph/0507007 (verify)) show the maximiser solves the complex Burgers equation.
  - **Payoff.** An integral over exponentially many tilings becomes a 2-D PDE, with fluctuations O(log N) (a Gaussian free field).
- **Noncommutative form (Theorem; Kenyon, "Spanning forests and the vector bundle Laplacian", Ann. Probab. 2011, arXiv 1001.4028 (verify)).**
  - **The generalisation.** Put an SU(2) connection on the graph. The quaternion determinant of the connection Laplacian is a sum over cycle-rooted spanning forests, weighted by the trace of the holonomy around each cycle. This is matrix-tree and Kasteleyn with matrix-valued (noncommutative) edge weights.
  - **Where it stays exact.** Exactness survives because SU(2) traces of holonomies are class functions, and the quaternion determinant is multiplicative on the relevant subalgebra.
  - **Connection to team A.** This is the same territory as the Clifford and quaternion permanent estimators, on the *easy* side of the Chien–Harsha–Sinclair–Srinivasan dichotomy.

### 1.2 Substitution tilings: Perron–Frobenius frequencies, gap labelling, and Connes' Penrose algebra

- **Frequencies (Theorem).** For a primitive substitution σ with substitution matrix M_σ, the hull is uniquely ergodic. Tile frequencies form the right Perron vector of M_σ, and patch frequencies are the Perron vector of the induced substitution on patches. The convergence rate of frequencies along a sequence of supertiles is set by |λ₂/λ₁|, the angle.
  - **The down–up skeleton.** Desubstitution forgets the finest level. Recognizability (Mossé; Solomyak) guarantees the forgotten level can be reconstructed uniquely.
  - **Payoff.** An ergodic average over an infinite tiling reduces to a d×d eigenproblem.
- **Gap labelling (Theorem).**
  - **Statement.** For a Schrödinger operator on a repetitive tiling with finite local complexity, the integrated density of states at a spectral gap lies in the trace image τ(K₀(C(Ω) ⋊ ℤ^d)). That image is the ℤ-module generated by the patch frequencies, i.e. μ(C(Ξ, ℤ)) over the canonical transversal Ξ.
  - **Status.** Bellissard conjectured it. It was proved independently by Bellissard–Benedetti–Gambaudo, Kaminker–Putnam (math/0205102 (verify)) and Benameur–Oyono-Oyono (math/0112113 (verify)).
  - **Local-to-global content.** A *global* spectral invariant is pinned to a countable module generated by *local* patch frequencies, because K-theory is homotopy-invariant and the trace on K₀ is integral over the transversal.
- **Connes' Penrose space (Theorem).** Penrose tilings modulo isometry form a "bad" quotient. Its C*-algebra is the AF algebra of the Fibonacci Bratteli diagram, with K₀ ≅ ℤ² ordered by the golden ratio, and the unique trace takes values in ℤ + ℤφ⁻¹. Traces on an AF algebra are exactly positive harmonic functions on its Bratteli diagram. For a *stationary* diagram the trace is the Perron vector.
- **Random substitutions (Theorem; Rust–Spindeler and others, 2017– (verify)).** Frequencies are almost surely the Perron vector of the *expected* substitution matrix. This is the annealed statement, and it is the version that transfers (§3.3).

### 1.3 Conic geometry: deletion–restriction, angle sums, intrinsic volumes, Birkhoff contraction

- **Deletion–restriction (Theorem; Zaslavsky 1975).**
  - **Statement.** r(A) = r(A∖H) + r(A^H), and χ_A(t) = χ_{A∖H}(t) − χ_{A^H}(t). The number of regions is |χ_A(−1)|.
  - **Klivans–Swartz (DCG 2011, arXiv 1001.3727 (verify)).** For a central arrangement, Σ over chambers C of the k-th conic intrinsic volume v_k(C) equals the absolute coefficient of t^k in χ_A(t). The Gaussian-measure statistics of all chambers *together* depend only on the intersection lattice.
  - **Valuation form (Derived).** Any simple valuation φ (Gaussian measure, Gaussian first moment) satisfies a deletion–restriction rule: removing H merges C₊ and C₋ into C. For a *piecewise-linear function* f the weights differ on the two sides, and the exact rule is f = f^{del} + f^{res}:
    - f^{del} replaces the kink at H by its face-averaged slope;
    - f^{res} has its second derivative supported on H (a facet term).

    Iterating over all hyperplanes is the faces stream's exact telescoping E5. Iterating twice gives the codimension-2 flats H ∩ H′, which are the "folding" and "curvature passage" terms the faces and markov streams found missing (16–20 % of κ3). The Gaussian form of the codimension-2 step is **Plackett's identity** (Biometrika 1954): ∂P(chamber)/∂ρ_ij = (Gaussian density of the codimension-2 face {z_i = z_j = 0}) × (conditional chamber probability). That is deletion–restriction in the correlation parameter.
- **Kinematic formulae and phase transitions (Theorem).**
  - **ALMT, "Living on the edge", arXiv 1303.6672.** For a convex cone C and a uniformly random rotation Q, P(C ∩ QD ≠ {0}) jumps from ≈ 0 to ≈ 1 when δ(C) + δ(D) crosses the ambient dimension, within a window O(√n). Here δ = E‖Π_C g‖² is the statistical dimension.
  - **McCoy–Tropp (arXiv 1308.5265 (verify)).** The intrinsic volumes concentrate at δ.
  - **Random tessellations.** Cover–Efron and Wendel give exact counts and angle statistics for random central arrangements; Hug–Schneider (arXiv 1508.07668 (verify)) give typical-cone and zero-cone laws.
  - **The elementary case used below (Sheppard).** For a bivariate standard normal with correlation ρ, P(both > 0) = 1/4 + arcsin(ρ)/(2π). This is the Gaussian solid angle of a 2-D wedge.
- **Birkhoff–Hopf (Theorem).**
  - **Statement.** A linear map taking a closed cone into itself, with projective diameter Δ of the image, contracts Hilbert's projective metric by tanh(Δ/4). Nonlinear Perron–Frobenius theory (Lemmens–Nussbaum, CUP 2012) extends this to order-preserving 1-homogeneous maps, which are nonexpansive.
  - **The down–up instance.** Sinkhorn's convergence proof (Franklin–Lorenz 1989): row normalisation then column normalisation is two conditional expectations, and the angle is tanh(Δ/4).
  - **Noncommutative version (Theorem; Reeb–Kastoryano–Wolf, arXiv 1102.5170 (verify)).** Birkhoff contraction on the PSD cone for positive and CP maps. Gurvits's operator scaling is its computational form.

### 1.4 The negative side: local rules that compute

- **Undecidability (Theorem).** Wang tiling is undecidable (Berger 1966). Robinson (1971) gives an aperiodic set whose hierarchical squares can host Turing machines. Cubitt–Pérez-García–Wolf (Nature 2015, arXiv 1502.04573) encode Robinson tilings plus quantum phase estimation into a translation-invariant 2-D Hamiltonian: the spectral gap is undecidable.
- **Synthesis.** Local-to-global fails when local rules can carry a *hierarchical* signal whose effect on the global quantity is not bounded by any local certificate. In the positive instances the global quantity is a *trace* of something homotopy-invariant or a contraction (Perron vectors, K₀ traces, Birkhoff). In CPW it is a ground-state energy *gap*, and no trace controls it.
- **The analogue for us (Synthesis).** A deep ReLU network with fixed weights can carry hierarchical signals in principle. The fresh-weight lemma is what rules that out on average. The quenched remainder is "computation-like" only to O(n^{-1/2}).

### 1.5 The common skeleton across these instances (Synthesis)

| instance | pieces | down step | up step | angle or certificate | global object |
|---|---|---|---|---|---|
| lozenge flips | hexagons | remove 3 rhombi | re-insert rotated | height-lattice contraction (Wilson) | uniform measure; limit shape (CKP) |
| Kasteleyn | edges | — (exact) | — | Pfaffian orientation, one local sign per face | perfect-matching count |
| substitution | tiles of level k | desubstitute | re-substitute | \|λ₂/λ₁\| of M_σ | patch frequencies = Perron vector |
| gap labelling | patches | restrict to transversal | — | K₀-integrality | IDS at gaps ∈ ℤ[frequencies] |
| Sinkhorn / Birkhoff | rows, columns | normalise rows | normalise columns | tanh(Δ/4) | doubly stochastic scaling |
| arrangement | hyperplanes | delete H | restrict to H | Möbius function of the lattice | χ_A, Σ v_k over chambers |

In every positive row the global object is either an exact algebraic identity or an *average* (ergodic, Gaussian or Perron). This decides what can transfer.

## 2. The noncommutative generalisation

### 2.1 Known
- **Tiling algebras.**
  - C(Ω) ⋊ ℤ^d, or the groupoid C*-algebra of the transversal (Kellendonk), with gap labelling as the trace image of K₀ (§1.2).
  - Pattern-equivariant cohomology (Kellendonk–Putnam; Sadun) computes K-theory from local patch data.
  - Penrose's algebra is AF, with traces equal to harmonic functions on the Bratteli diagram.
- **Dimers.** Kenyon's connection Laplacians: quaternion determinants, with holonomy traces as the weights (§1.1).
- **Cones.** Hilbert's metric on the PSD cone and Birkhoff contraction for CP maps (Reeb–Kastoryano–Wolf), and operator Sinkhorn (Gurvits).

### 2.2 The noncommutative object behind the wedge calculus (Derived)
For a layer with gate pattern G_l = diag(g_l(x)), a projection, the transport of a quenched error δ by the fresh layer is the map δ ↦ δ W G_l. Average over the fresh W and over the input pair (x, x′), and the second-moment action on errors becomes the completely positive map
Θ_l(X) = 2 · E_{x,x′}[G_l(x) X G_l(x′)] (restricted to the diagonal algebra for mean errors; the Hadamard/Schur multiplier by E[g_i g_j] for the off-diagonal covariance).
- **Perron value.** Its Perron value on the diagonal is λ_l = 2 τ(G(x)G(x′)) = 1/2 + arcsin(ρ_l)/π. This is the trace of a product of **two projections at angle θ_l = arccos ρ_l**, which is exactly the two-projection (dihedral) algebra of the input mandate. The angle that measures how far "forget this" and "forget that" are from independent is the wedge angle.
- **Criticality.** He-ReLU is the unipotent case (team C's "wall at 2"): θ_l → 0 like 3π/l. The CP maps have no uniform Birkhoff gap. Their product contracts like (s/L)^3 for mean errors and (s/L)^6 for off-diagonal errors, not geometrically.

### 2.3 The precise next statements
- **Proposition E1 (Derived; proof sketch in §4.1).**
  - **Hypotheses.** Bias-free ReLU MLP; i.i.d. N(0, 2/n) weights; x ~ N(0, I_n); first order in the injected error; n → ∞ at fixed L.
  - **Conclusion (means).** For a mean error δm injected at the post-activation of layer s, E_W ‖δm_L‖²/‖δm_s‖² → Π_{l=s+1}^{L} (1/2 + arcsin(ρ_l)/π), with ρ_0 = 0 and ρ_{l+1} = (√(1−ρ_l²) + (π − arccos ρ_l) ρ_l)/π.
  - **Conclusion (off-diagonal covariance).** For an off-diagonal covariance error, K_off(s) = (8/n²) Σ_{k>s} [Π_{s<u<k} λ_u²] q_k K_mean(k), with q_k = E[φ(t)²]/(4σ_k²), t ~ N(0, ρ_k/(1−ρ_k)), σ_k² = 2(1−ρ_k).
- **Conjecture E2 (noncommutative Birkhoff at criticality).**
  - **Setting.** Let Θ_l be the CP maps of §2.2 on M_n(ℂ) with the PSD cone.
  - **Conjecture.** Their Hilbert-metric contraction coefficients satisfy 1 − c(Θ_l) = θ_l/π + O(θ_l²). Hence the composition Θ_L ∘ ⋯ ∘ Θ_{s+1} has projective contraction Π(1 − θ_l/π) ≈ (s/L)^3, uniformly in n.
  - **Consequence.** The power-law memory is a property of the cone geometry of the gate projections, not of any basis. Any carrier of quenched content must keep a polynomially long memory.
- **Conjecture E3 (annealed gap labelling for the network; a rigidity statement).**
  - **Setting.** Let A_L be the algebra generated by the gate projections and the fresh weights. The ensemble trace is the expectation over W and x.
  - **Conjecture.** Every ensemble-averaged quantity built from the moment chain (K coefficients, energies per age, frequencies of neuron types (t_j)) is a polynomial in the wedge values {arcsin ρ_l / π} and Gaussian facet integrals. That is, it lies in an explicit "trace image".
  - **Consequence.** The *quenched* quantities the competition needs have zero component there (F8.5).

  This makes precise why the tiling or K-theory transfer can price old content but cannot compute it.

### 2.4 The verdict for this domain (Synthesis)
- **What the domain's theorems compute.** Averages over a group or an ergodic hull: translations for tilings, rotations for kinematic formulae, the Perron direction for Birkhoff.
- **What our ensemble provides.** It has exactly such a group: rotation invariance of the fresh weights (the fresh-weight lemma). The averages are the annealed (kernel or wedge) quantities, and §4 shows they are now computed essentially exactly.
- **What the competition needs.** The quenched remainder, which has zero ensemble mean.

So the tiling and cone theorems transfer as **pricing, budgeting and structural laws**, not as cheaper carriers of the binding object, unless a carrier is found that exploits the power-law (critical) memory itself (§3.4).

## 3. Transfers to our problem

### 3.1 T1: the wedge calculus as the exact price list (tested, positive)
- **Roles.**
  - Pieces: the two replicas' gate projections G_l(x), G_l(x′).
  - Down step: the conditional expectation over the fresh layer W_{l+1}, which forgets everything but rotation invariants (fresh-weight lemma).
  - Up step: re-gating at layer l+1.
  - Angle: θ_l = arccos ρ_l, the wedge angle.
  - Global quantity: the final-MSE price of any local error.
- **Why it is not the per-neuron drawing.** No neuron or layer state is drawn. The object is the pair of input replicas and the 2-D wedge each fresh neuron cuts, a conic intrinsic quantity of the input space, iterated along the arc-cosine orbit.
- **Prediction and test.** `wedge_transfer.py`, log `wedge_transfer.log`. It reproduces region's K_mean(l) and K_off(l) (MLPs 0–2, n = 1024): mean log-ratio of measured to predicted +0.04 (K_mean) and +0.01 (K_off), with sd 0.11 for both, which is within the MLP-to-MLP scatter. No fitted parameter. The shape of K_off is reproduced too: rise from 2.5e-8 to a peak of 7.6e-7 at l = 11, then a fall to 3.9e-7 at l = 14. It comes from the two-step path off-diagonal → next pre-activation variance → mean.
- **Established facts explained.**
  - Region §2: "mean errors contract by ≈ 0.8 per layer". This is the wedge factor λ_l = 0.60 → 0.87, i.e. a power law, not a constant factor.
  - N4 and the "early layers need almost nothing" item.
- **Use.**
  - K(l) is now analytic at any width and depth (including the 256 × 32 smoke shape), so error budgets can be set *before* building a design.
  - Region's hypothesis (b) for the leaders, "cruder D21 where K_off is small", can be evaluated exactly: layers 0–2 and 14 carry 13 % of the Σ K_off weight.
- **Cost.** O(L) scalars.
- **Kill criterion.** None for the price list; it already matched. For E2: if at depth 32 (256 × 32 bench) the measured K_mean deviates from Π λ by more than 2×, the critical power law is wrong.

### 3.2 T2: the age law of old content (tested against stream numbers, positive)
- **Roles.**
  - Pieces: facets, the bent hyperplanes z_{s,r} = 0 where births happen (codimension 1). Each birth has strength = facet density on the data cone.
  - Down step and angle: the transport of each of the three legs by one wedge factor per layer.
  - Global quantity: the age spectrum of D21 at each target layer.
- **Prediction.** `wedge_ages.py`, log `wedge_ages.log`. Birth strength falls as B_s/B_0 = 1, 0.41, 0.15, 0.065, 0.031, … The data cone narrows: σ_s² = 2(1−ρ_s) shrinks, so fewer facets cut it. The three-leg transport per layer λ³ rises from 0.22 to 0.67 by criticality. The age shares at layer 14 are 0.14, 0.12, 0.11, 0.10, 0.09, … for ages 1, 2, 3, …: nearly flat.
- **Comparison with measurements.**
  - Dropping ages > 4 predicts a D21 relative error of 0.72; region N6 measured 0.75.
  - Dropping ages > 2 predicts 0.86; measured 0.89.
  - Faces T3 (width 64) measured a uniform age profile.
- **Established facts explained.** BRIEF §3 "every age matters", N6 and T3, from two conic numbers per layer.
- **What it predicts beyond those.**
  - **At greater depth** (L = 32, the smoke shape): λ³ → 1 − 9/l and B_s ∝ σ_s^7 ~ s^{-7}. The age profile tilts towards recent ages like (s/t)^9 × s^{-7}·…, so it stays polynomially flat. Old content remains necessary at depth 32. A window of w layers loses a share that is a fixed fraction, about 1 − w/t, not exponentially small.
  - **At greater width**, nothing changes. The law is width-independent at leading order, consistent with "old content does not shrink with n".
- **Kill criterion.** Instrumenting FC at n = 1024 to measure the per-age Frobenius share of D21 (the next test, §4.3): if the measured shares at layer 14 differ from the predicted ones by more than 2× at ages 1–8, T2 is wrong.

### 3.3 T3: substitution / Bratteli structure of the old-content reservoir (structural; negative for cost, positive for pricing)
- **Roles.**
  - Level = layer.
  - The "patch" state is the third-order reservoir T(t) ∈ (ℝⁿ)^{⊗3}, the CP sum of FC's atoms (y_r, y_r, z_r) over all sources.
  - Connecting map: M_t^{⊗3} with M_t = diag(Φ_t) W_{t+1}, the same matrix on all three legs.
  - Births are added at each level: T(t+1) = (T(t) + B_t) M_t^{⊗3}.
  - The readout is a "trace": the partial diagonal T_{aab} (D21) and the diagonal T_{aaa}.
- **Why it is not the per-neuron drawing.** The object is one tensor state per layer, a stationary Bratteli-type system, with the readout as a trace functional. The neurons enter only as the basis in which the trace is taken.
- **The tiling answer.** In a substitution tiling, a stationary Bratteli diagram makes traces Perron vectors: a d-dimensional computation. Here the connecting maps are dense Gaussian (fresh at every level, no recurrence of "patches"). The annealed version is stationary and is computed exactly by T1–T2. The quenched version has no finite "dimension group", because every level reads an independent random projection of an n³-dimensional reservoir (F8.3).
- **What follows.**
  - Information count: L·n² numbers are read from a reservoir of CP rank n·L. This is consistent with region §6: no compression below ≈ n atoms per age.
  - The positive part is pricing. The energies per (source, target) pair are known analytically (T2), so the optimal allocation of rank or accuracy across pairs is computable without Monte Carlo.
- **Kill criterion for "the reservoir has a Perron structure".** In FC at n = 1024, project the old-content D21 on the top-k singular directions of the gated propagator Z_s(t). If k ≤ n/8 does not capture ≥ 90 % for ages ≥ 4, there is no Perron shortcut. Region's per-source rank-256 test (4× worse) already points that way.

### 3.4 T4: a critical-memory carrier (speculation, to test in the final)
- **The idea.** At criticality the transport factor per leg is 1 − θ_l/π with θ_l ≈ 3π/l. A unipotent (Jordan-block) system with power-law memory is the setting where a **multiscale / Horner-type** summation beats explicit (source, target) pairs: memory decaying like (s/t)^p is approximately self-similar in log-depth.
- **Proposal.** Merge sources in dyadic blocks of age (1, 2, 3–4, 5–8, 9–16), as in fast multipole or hierarchical-matrix treatment of power-law kernels. Carry one merged CP atom set per block, re-fitted when blocks merge. Pairs then fall from O(L²) to O(L log L): 120 → ≈ 16 × 4 = 64 at L = 16.
- **What makes it differ from region's failed pruning.** Region pruned atoms *within* a source. This merges *across* sources, with the merge schedule set by the power law.
- **Prediction.** If the energy within a dyadic block is dominated by its near-uniform age profile (T2) and the blocks' legs are close in the Z-propagator geometry (successive sources differ by one wedge factor), a block merge with R = n atoms keeps ≥ 90 % of the block's D21.
- **Cost model.** ≈ 64 pairs × 5 units ≈ 320 u dense, ≈ 180 u with Strassen L5, against FC's 840 / 500.
- **Kill criterion.** In FC at n = 1024, merge ages 5–8 and 9–16 into R = n atoms by ALS. If raw degrades by more than 1.5× (≥ 4.5e-8), T4 is dead.

## 4. Tests run

### 4.1 Derivation behind T1 (Derived)
1. Given W, the gates of two independent inputs are independent, so P(g_j(x) = g_j(x′) = 1 | W) = Φ_j². Averaging over j and over the fresh weights gives the annealed two-replica orthant probability 1/4 + arcsin(ρ)/(2π) (Sheppard).
2. A mean error δm moves the next pre-activation mean by δm W, with ‖δm W‖² ≈ 2‖δm‖² by the fresh-weight lemma. It then moves the post-activation mean by Φ_j δμ_j. Together this gives the factor 2·mean(Φ²) = λ.
3. For an off-diagonal covariance error δC:
   - the off-diagonal image transports with Schur weights E[g_i g_j] ≈ Φ_iΦ_j, giving the factor λ²;
   - the diagonal of W^T δC W has variance (8/n)‖δC‖²_F in total;
   - a variance error moves the mean by φ(t)/(2σ) per unit of σ².
4. Pre-activation second moment = 2 per neuron (He), split as 2ρ (squared mean) and 2(1−ρ) (variance). This gives t_j ~ N(0, ρ/(1−ρ)).
5. The final MSE is per neuron, hence the extra 1/n in K_off.

### 4.2 Numbers (Measured against stream data)
- **K table.** `wedge_transfer.py` / `wedge_transfer.log`: see §3.1.
- **Age law.** `wedge_ages.py` / `wedge_ages.log`: see §3.2.

Both use only region's published per-layer coefficients and N6 numbers. No new network evaluation was needed, and both run in < 1 s.

### 4.3 Next (for the final report)
- **(a)** Instrument FC (`notes/fresh-slate/breakthrough/region/fc.py`, copied into this directory, not edited) at n = 1024 on bench MLP 0 to record the per-(s, t) Frobenius energy of the D21 contributions. This tests T2 directly.
- **(b)** The T4 dyadic-merge test on the same MLP, priced as in §3.4.

## 5. Honest assessment
- **Theorem.** §1 statements (sources given; some arXiv ids to verify).
- **Derived.** Proposition E1 (first order, n → ∞, annealed). Not a rigorous proof: concentration of Φ_j² averages and the neglect of neuron-to-neuron spread of σ_j are assumed.
- **Measured.** E1 reproduces region's K(l) (3 MLPs, 15 layers) to 11 % rms in log, and the age model reproduces N6 to within 4 %. The age model uses one Wick term, a diagonal approximation of S_s and incoherent sums. That it lands within 4 % of N6 may be partly lucky; test (a) decides.
- **Synthesis.** The wedge angle is the two-projection angle of the mandate. Criticality is the unipotent case. Tilings and cones transfer as pricing, not carriers.
- **Speculation.** T4 (dyadic merging across sources) and Conjecture E3.
- **Risks.**
  - T1 and T2 explain facts but do not lower cost.
  - T4 may fail for the reason region's pruning failed: atom completeness may hold across sources too.
