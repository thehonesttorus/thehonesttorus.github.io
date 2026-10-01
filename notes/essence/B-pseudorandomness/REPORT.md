# Team B: pseudorandomness of walks and of layered computation

*Status: v1.1, final for this session (2 Oct 2026); v1.1 adds T3b, the Richardson range no-go. v0 had §§0–5 without the T2 measurements; v1 adds the T2 tests
(§3 T2, §4), closes T6, and adds the second network for T4. Contract: [../BRIEF.md](../BRIEF.md); mandate: [../INPUT-2026-10-01-local-to-global.md](../INPUT-2026-10-01-local-to-global.md).
Tests: [tests/](tests/). Every number quoted from other streams carries its source; everything else is marked
THEOREM (checked in the primary source), DERIVED (derived here, sketch given), MEASURED (run here), SYNTHESIS or
CONJECTURE.*

## 0. Summary

1. **The essence in pseudorandomness is the hybrid argument**, and every instance below is a refinement of it:
   swap one piece of the true process for one piece of the idealised one, bound the change, and sum over pieces.
   - Lindeberg / invariance principles swap one coordinate.
   - INW/Nisan swap one half of a branching program's randomness.
   - Jeronimo–Mittal–Roy (JMR) "ignore the first step" deletes the first vertex of a walk.
   - Richardson / weighted PRGs correct an approximate product by inserting one exact step.

   What decides the payoff is the **angle** in each instance:
   - the spectral gap λ of the walk;
   - the third moment in Lindeberg;
   - the width w of the program;
   - ε₀ = ‖I − M L‖ for Richardson;
   - the smallest irrep dimension D for quasirandom groups.
2. **Our problem has the hybrid argument in an unusually strong form, and none of the gaps.**
   - Fresh weights make the hybrid increments orthogonal (Schur orthogonality, i.e. the fresh-weight lemma), so the
     errors of the layers add in squares. That is why the region stream's books close (region §3).
   - The analogue of JMR's λ is the fluctuation of the "squared graph" P D Pᵀ of the gated transport. That is a
     square Wishart matrix: its fluctuation is O(1) in operator norm and equal to its mean in Frobenius norm, so
     **λ_eff ≈ 1**.
   - The analogue of INW's width is the Frobenius rank of the old content, which is **∝ n** (F9.3).
   - The analogue of Richardson's ε₀ for any *projective* compression is not even defined. Projections are fixed
     points of the error-reduction map (§3.3).

   The pseudorandomness theorems therefore give *explanations* of the established facts (N6, N7, F9, costate C3) and
   clean negatives for the incoherent half of the memory. For the coherent half they give a working carrier (item 6).
3. **The one positive structure: D-quasirandomness.** The fresh layer's O(n) symmetry splits a third-order tensor
   into its trace (vector) part, of dimension n, and its traceless part, of dimension ≈ n³/6. Readouts with coincident
   indices see the trace part with a gain Θ(n) over its Frobenius share.
   - This is exactly JMR's "fraction of the trivial irrep" factor, and it explains why κ4 is needed only through
     n-vectors (region N7).
   - It predicts the same for every coherent channel.
   - It cannot reduce the incoherent half of D21, which the final means do see (region §3: ΔC is 98 % incoherent).
4. **The JMR squared graph appears literally in the costate** (§3.2, DERIVED at mean-field gates).
   - The all-age effect of third-order content on output c is a sum over sources only.
   - The sum over targets collapses into a backward Lyapunov sandwich 𝓜ᶜ_s = P_s [diag(e∘vᶜ)_{s+1} + 𝓜ᶜ_{s+1}] P_sᵀ.
     This is "integrate out the vertex shared by the two copies", i.e. the squared graph.
   - Its price is one n × n state per output, so forward costs L² n³ (source–target pairs) and backward costs L n⁴
     (per output). This is costate C3 seen from the pseudorandomness side.
5. **Test (MEASURED, n = 1024, bench MLP 0):** a Z₂-structured design (antithetic pairs) fools every odd function
   exactly, the analogue of JMR's Theorem 1.7 for Cayley expanders.
   - Its variance gain for the sampled D21 is **0.98–1.04× at layers 3–16** (N = 32,768 forward passes).
   - The reason: the even part of the centred pre-activations already holds **η = 33 % (layer 3) to 45 % (layer 16)**
     of their variance, so parity designs remove nothing at depth.
   - The sampled noise reproduces F6.10 (0.29 at layer 10, 0.19 at 15).
   - MLP 1 agrees: gain 0.98–1.04×, η = 0.33–0.45.
6. **JMR's mean term carries the coherent half of the memory at O(n²) per layer** (MEASURED inside FC, n = 1024,
   bench MLP 0).
   - **The split.** Split the old (age > 2) part of D21 into its part constant in the repeated index a, i.e. the
     O(n)-average over the last fresh layer, which equals Cov(norm process, z_b), and the Wishart fluctuation.
     - The constant part holds **79–81 % of the old D21 energy** at layers 8–15.
     - Each part is needed: dropping it costs **6.8e-7** (constant part) and **1.0e-6** (fluctuation), against FC's
       3.24e-8. The two penalties add up to the whole-memory penalty of 1.77e-6.
   - **The carrier.** JMR's mean term (the squared graph averaged over the fresh layer) freezes each source's row
     scalars when it turns old and multiplies them by one scalar gain per layer, ⟨P² c⟩. That makes the constant part
     a **single n-vector transported by P_t: O(n²) per layer, ≈ 0.02 units in total**.
   - **The result.**
     - With the constant part from this recursion and the fluctuation exact, FC scores **3.67e-8** (+13 %).
       On MLP 1 it scores 1.84e-8 against 1.81e-8 (+2 %).
     - The vector's own error is 2 % at layers 4–7 and grows to 9.5 % at layer 15.
     - Keeping only ages ≤ a exact and the constant channel beyond them gives 1.0e-6, 4.0e-7 and 1.8e-7 at
       a = 2, 4, 6. That is 1.8–2.1× better than dropping old content (region: 8.5e-7 at a = 4), but still
       6–30× above FC.
   - **The verdict.** The incoherent (Wishart) half of the memory, ≈ 20 % of its D21 energy, is the irreducible
     part, and the fresh-weight lemma prices it in full. The coherent half is solved.

## 1. Instances (depth over breadth)

### I1. Jeronimo–Mittal–Roy: the generalised "ignore the first step" (arXiv 2507.14445), THEOREM

**Statement (Thm 3.1, read in the source).**
- Setting:
  - X is a λ-spectral expander, and x⃗ is a random walk on it.
  - S = {i₁ < … < i_k} is a set of indices with gaps Δ_j = i_{j+1} − i_j.
  - f_j : X → M_{d_j}(ℂ) satisfy E_X f_j = 0 and ‖f_j(x)‖_op ≤ 1.
- Conclusion:

  ‖E_{x⃗ ∼ RW}[f₁(x_{i₁}) ⊗ ⋯ ⊗ f_k(x_{i_k})]‖_op ≤ Σ_{I ∈ 𝓘_k} λ^{Σ_{i∈I} Δ_i},

  where 𝓘_k is the set of index sets I with {1, k−1} ⊆ I ⊆ [k−1] that contain at least one of every two
  consecutive indices.
- Corollary 3.2: for non-trivial unitary irreps of a group G labelling the vertices, every f_j = ρ_j ∘ φ satisfies
  the hypotheses.
- Fooling results:
  - symmetric functions over any alphabet are fooled at O(|Σ| λ);
  - group products 1{x₁⋯x_k = t} at (2λ)^{k/2} (Thm 4.10);
  - for D-quasirandom G on "pseudo-Cayley" graphs, symmetric class functions are fooled at O(√|G| λ / D)
    (Thm 1.8 / Prop 5.5);
  - the bound carries the extra factor ⟨χ_triv, χ₁⋯χ_k⟩, the fractional multiplicity of the trivial irrep in
    ρ₁ ⊗ ⋯ ⊗ ρ_k, which is ≤ 1/D² for k = 2;
  - lower bound Ω(λ) for every group (Thm 1.9).

**The pieces, the down and up steps, the angle.**
- Two quantities are carried:
  - ε_j = ‖E N_j‖, the mean of the "future" N_j(y) = E[g_j | x_{i_j} = y];
  - ζ_j² = E‖N_j(y)‖², its second moment.
- **Down step** (bounding ε_j): split the walk at x_{i_{j+1}} and apply the expander mixing lemma (EML) on the
  graph X^{Δ_j}. This gives ε_j ≤ λ^{Δ_j} ζ_{j+1}.
- **The "ignore the first step" move** (bounding ζ_j):
  - E_y ‖N_j(y) v‖² is the overlap of *two independent continuations from the same vertex y*.
  - The factor h_j(y) costs at most ‖f_j‖²_op ≤ 1 (unitarity), and then y is integrated out.
  - Two walks from a common vertex are one walk of length 2Δ_j between their endpoints, i.e. on the **squared graph**.
  - EML on it gives ζ_j ≤ ε_{j+1} + λ^{Δ_j} ζ_{j+1}.
- **Angle:** λ^{Δ} is the cosine between "functions of x_i" and "functions of x_{i+Δ}" after removing the mean.
- **Global:** the recursion closes into λ^{η(S)}.

**Computational payoff.**
- Expander walks are PRGs with seed log|V| + k log d instead of k log|V| for symmetric, product and word tests,
  which is used in code distance amplification and in derandomised Chernoff bounds.
- The quasirandom gain is the dimension factor 1/D: noncommutativity spreads mass over large irreps that a
  low-complexity test cannot see.

### I2. Precision amplification for products of transition matrices (AKMPSV, arXiv 1912.04524), THEOREM

**Statement.**
- Lift a length-k walk to L = I − P_k ⊗ W, where P_k is the shift on the path.
- L⁻¹ is block lower-triangular with blocks W^{ℓ−s}, so every power of W is a block of L⁻¹.
- Given any B with ‖I − B L‖_F ≤ α < 1, P_m = Σ_{i≤m} (I − B L)^i B satisfies ‖I − P_m L‖_F ≤ α^{m+1}
  (Lemma 6.2, read in the source).
- So low-precision approximations of the products (error 1/poly) become ε-precision at O(log(1/ε)) extra terms,
  i.e. space O(log N · log log(1/ε)) (Thm 1.3, §8).
- The approximate inverse comes from an LU / Schur complement that eliminates every other layer of a cycle-lifted
  graph. This turns C_{2^k} ⊗ W into C_{2^{k−1}} ⊗ W², i.e. integrate out one layer and get the squared walk.
- That squared walk is in turn replaced by the Rozenman–Vadhan derandomised square, which is a unit-circle
  approximation (Thm 5.9).

**The pieces, the down and up steps, the angle.**
- **Pieces:** the layers of the lifted graph.
- **Down step:** Schur-complement out the odd layers (exact marginalisation).
- **Up step:** substitute a sparse approximate square, then correct with one exact application of L per
  Richardson round.
- **Angle:** α = ‖I − B L‖, measured in a norm F adapted to the Laplacian. This is the only place where expansion
  enters.
- **Global:** inverse-polynomial-precision walk probabilities in Õ(log N) space for Eulerian graphs.

**Weighted PRG relatives** (cited from memory, to be checked if used):
- Braverman–Cohen–Garg, arXiv 1711.08911: signed weights;
- Cohen–Doron–Renard–Sberlo–Ta-Shma (CCC'21): error reduction for weighted PRGs;
- Pyne–Vadhan (CCC'21): pseudodistributions that beat all PRGs;
- Hoza (RANDOM'21): better pseudodistributions.

They all use the first-order error-reduction identity: with Ã[i, j] ≈ A_i ⋯ A_j,

A[1, k] ≈ Σ_i Ã[1, i−1] A_i Ã[i+1, k] − Σ_i Ã[1, i] Ã[i+1, k],

which has error O(ε²). Each term uses approximate segments and **one exact step**. That is the "forget one piece,
reconstruct it exactly" move, iterated.

### I3. Lindeberg replacement, invariance and fooling polytopes, THEOREM (from the standard sources)

**Mossel–O'Donnell–Oleszkiewicz** (arXiv math/0503503).
- Statement: for a multilinear polynomial Q of degree d with low influences τ and a test ψ with bounded third
  derivative, |E ψ(Q(X)) − E ψ(Q(G))| ≤ O(d · 9^d · τ^{1/2}) ‖ψ‴‖_∞.
- **Pieces** are coordinates.
- **Down / up:** replace X_i by G_i (Lindeberg).
- **Angle:** the per-swap change is a third-order Taylor term, weighted by hypercontractive moment bounds.
- **Global:** the sum over swaps.

**Harsha–Klivans–Meka** (arXiv 0912.4884) and **O'Donnell–Servedio–Tan** (arXiv 1808.04035) fool intersections
of k halfspaces, i.e. polytopes.
- **Bentkus' mollifier** smooths the indicator of a polytope with derivatives polylog(k).
- **Nazarov's bound** puts the Gaussian surface area of any k-facet polytope at O(√log k). It converts smoothing
  error to probability error at only a √log k cost.
- **Seed length** polylog(k)/ε^{O(1)}, and in O'Donnell–Servedio–Tan poly(log k, 1/ε) · log n.

The point for us is that **the dependence on the number of facets is logarithmic**, because the mollifier is
smooth uniformly in the direction of the facets.

### I4. INW / Nisan, THEOREM (standard; Impagliazzo–Nisan–Wigderson STOC'94, Nisan Combinatorica 1992)

- **Statement:** a read-once branching program of width w and length L is fooled to error ε by recursively
  recycling the second half's seed through an expander edge. The seed length is O(log L · log(wL/ε)).
- **Pieces:** the halves.
- **Down step:** condition on the midpoint state, which is only log w bits.
- **Up step:** regenerate the second half's randomness from an expander neighbour.
- **Angle:** the expander's λ against the w-state midpoint distribution.
- **Global:** the recursion over log L levels.
- **Payoff:** BPL ⊆ SC (Nisan), and Saks–Zhou BPL ⊆ L^{3/2}.

The load-bearing quantity is the **width**: the second half is nearly independent of the first *given a short
summary*.

### I5. Quasirandom groups (Gowers, arXiv 0710.3877; Bourgain–Gamburd; Sarnak–Xue), THEOREM (standard)

- **Statement (Gowers):** if every non-trivial irrep of G has dimension ≥ D, then for f, g of mean zero,
  ‖f ∗ g‖₂ ≤ √(|G|/D) ‖f‖₂ ‖g‖₂. Product sets of density > D^{-1/3} cover G.
- **Mechanism:** Schur orthogonality, E_g ρ(g)_{ij} conj(ρ(g)_{kl}) = δ_ik δ_jl / d_ρ. A large irrep spreads
  the mass of every fixed vector.
- **Sarnak–Xue / Bourgain–Gamburd:** large multiplicities of irreps plus few small eigenvalues give a spectral gap
  (expanders from SL₂(𝔽_p)).
- For us, **the fresh Gaussian weight is a Haar-like O(n) element with radial part**. Its Weingarten / Schur
  orthogonality *is* the fresh-weight lemma (§2.2).

## 2. The noncommutative generalisation

### 2.1 What is known

- **Matrix-valued (operator-valued) test functions:** JMR Thm 3.1 above is already noncommutative in the test
  functions. Matrix expander Chernoff (Garg–Lee–Song–Srivastava, arXiv 1704.03864) is its large-deviation sibling.
- **Quantum expanders** (Hastings, arXiv 0706.0556; Ben-Aroya–Schwartz–Ta-Shma, arXiv 0709.0911) replace the walk
  by a unital channel with a gap relative to the completely depolarising one.
- **Finitely correlated states** (Fannes–Nachtergaele–Werner 1992) are the noncommutative Markov chains. Their
  two-point correlations decay as the gap of the transfer channel to the power of the distance.
- **Weingarten calculus** for Haar O(n) and U(n), and asymptotic freeness of Haar-rotated matrices (Voiculescu).
  This is the noncommutative form of "fresh randomness decorrelates" (team C's domain).

### 2.2 The fresh-weight lemma is Schur orthogonality (DERIVED, standard)

Let W have an O(n)-bi-invariant law with E W W ᵀ = (2/n) n I. Decompose Sym³(ℝⁿ) = V_(3) ⊕ V_(1):
- the traceless part, of dimension C(n+2, 3) − n;
- the trace part sym(δ ⊗ v), which is an n-vector v.

For the coincident-index readout D21(T)_ab = (T ×³ W)_aab:
- E_W[D21] = 0 (odd);
- E_W ‖D21(T)‖²_F ≈ (16/n) ‖T_(3)‖²_F + Θ(1) ‖T_(1)‖²_F.

The trace part is transported coherently, through the W_ia² → 2/n diagonal, with gain Θ(n) over its Frobenius share.
- The first term is the region stream's "‖δD21‖² ≈ (16/n) ‖δκ‖²".
- The second is the "coherent index sums travel O(√n) more strongly" of region N7.

**This is JMR's D-quasirandom factor.** A test sees a large irrep only through its trivial-irrep content, at 1/D.
Here the coherent readouts (column means, traces) see V_(3) at 1/n and V_(1) at O(1).

### 2.3 The precise next statement (CONJECTURE, with proof sketch)

**Noncommutative ignore-first-step theorem (a quantum-Markov-chain version of JMR Thm 3.1).**

Setting:
- 𝒜 is a von Neumann algebra with a faithful normal state φ, and 𝓔 : 𝒜 → 𝒜 is a φ-preserving unital completely
  positive (UCP) transfer map.
- λ := ‖𝓔 − φ(·)1‖ on L²(𝒜, φ), in the KMS (or GNS) inner product. Assume 𝓔 is KMS-symmetric.
- For gaps Δ₁, …, Δ_{k−1} and f_j ∈ 𝒜 ⊗ M_{d_j} with (φ ⊗ id) f_j = 0 and ‖f_j‖ ≤ 1, define the multi-time
  correlation

  ω_S(f) := (φ ⊗ id)(f₁ · 𝓔^{Δ₁}(f₂ · 𝓔^{Δ₂}(f₃ ⋯ 𝓔^{Δ_{k−1}}(f_k)⋯))).

Then

‖ω_S(f)‖_op ≤ Σ_{I ∈ 𝓘_k} λ^{Σ_{i∈I} Δ_i}.

Proof sketch:
- **ε step:** |(φ ⊗ id)(f_j 𝓔^Δ(N))| = |⟨f_j*, (𝓔^Δ − φ)(N)⟩_φ| ≤ λ^Δ ‖f_j‖₂ ‖N‖₂.
- **ζ step:** the Kadison–Schwarz inequality 𝓔(a*a) ≥ 𝓔(a)*𝓔(a) replaces "‖h_j(y)‖_op ≤ 1 then integrate y
  out". The two continuations from the shared vertex become 𝓔^{Δ*} ∘ 𝓔^{Δ}, the noncommutative squared graph,
  whose gap is λ^{2Δ} under KMS symmetry.
- **Reduction:** for 𝒜 = L^∞(X) with 𝓔 the walk it is JMR exactly.
- **Open part:** whether the KMS-symmetry hypothesis can be weakened to the "unit-circle approximation" of AKMPSV
  (the directed case), which would give a directed or non-reversible version.

**Payoff, if true.**
- Tensor-network (MPS / finitely correlated) states with a gapped transfer channel ε-fool symmetric multi-point
  observables with JMR's combinatorics η(S) rather than ⌊k/2⌋.
- The D-quasirandom refinement holds with the trivial-subrepresentation multiplicity of the channel's
  multiplicative domain.

**Relevance to us (honest).** Our layer transport has λ_eff ≈ 1 (§3.2), so the conjecture's *gap* half does not
apply. Its *dimension* half (2.2) does, and it is already the fresh-weight lemma.

## 3. Transfers to our problem

The memory object, in the notation of the region stream §5:
- P_l = diag(Φ_l) W_{l+1} is the gated one-step transport.
- The old third cumulant at layer l is T_l = Σ_{s<l} B_s ×³ P_{s→l}, with birth tensors
  B_s = Σ_r w2_{s,r} sym(y_r ⊗ y_r ⊗ e_r).
- D21_l = (T_l)_{aab}.

T_l obeys the exact one-step recursion

T_{l+1} = T_l ×³ P_l + B_{l+1},

i.e. the lifted system (I − 𝒮) T = B with 𝒮 nilpotent. **FC computes the exact Neumann series.**

### T1. The orthogonal hybrid argument (Lindeberg / INW / region §3)

**Objects in each role.**

| role | object |
|---|---|
| pieces | the 16 per-layer closures: replace the true law of a_l by the closure's law, given the true law of z_l |
| down step | Lindeberg swap at layer l |
| up step | regenerate the rest of the network with the closure |
| angle | orthogonality of the hybrid increments over the weight filtration F_l: increment l is a martingale difference in W_{l+1}, W_{l+2}, …, which are fresh |
| global | final MSE = Σ_l K(l) · ‖local_l‖² |

**Not the per-neuron drawing.** The pieces are whole-layer *law replacements*. The angle comes from the fresh
weights' symmetry, not from any neuron graph.

**What it predicts.**
- PRG theory pays L (triangle inequality) for L hybrid steps. Ours pays √L, because the increments are orthogonal
  in expectation over the weights.
- This is the region stream's closing of the books to 2–6 %. It explains region §3.
- It also predicts where the books *fail*: the cross terms E⟨incr_l, incr_m⟩ are non-zero only through components
  that are *not* averaged by fresh weights, i.e. the coherent, low-isotypic channels (§2.2: the mean or Perron
  direction and the trace vectors).
- So MLP 1's 46 % over-prediction (region §3) should be carried by the coherent channel.

**Cheapest decisive test.**
- Project each layer's local error ΔC_l onto the coherent subspace (the span of μμᵀ, the column means, and the
  diagonal), price the cross terms between layers, and check that they account for MLP 1's discrepancy.
- Cost: diagnostic only, on existing atlases.
- **Kill:** if the cross terms sit in the incoherent part, then the orthogonality reading of the hybrid argument
  is wrong at n = 1024.

### T2. The JMR squared graph in the costate, and why the gap is absent (DERIVED)

**Objects in each role.**

| role | object |
|---|---|
| pieces | the coincident index a in D21_ab = T(a, a, b), playing JMR's "shared vertex y with two continuations" |
| down step | integrate a out against the downstream sensitivity |
| angle | the fluctuation of the squared operator P D Pᵀ around its O(n)-average |

**The derivation, under mean-field gates.**
- Linearise the covariance arrow with J_ab ≈ Φ_a Φ_b. Then the backward sensitivity of the final variance of output
  c to the layer-l covariance is rank one: Gᶜ_l = vᶜ_l vᶜ_lᵀ, with vᶜ_l = Φ_l ∘ (W_{l+1} vᶜ_{l+1}).
- With ΔC_l = ½ (e ∘ D21 ∘ Φ + transpose), e = E δ(z), the first-order effect of all third-order content on
  output c is

  s_c = Σ_s Σ_r w2_{s,r} (y_{s,r}ᵀ 𝓜ᶜ_s y_{s,r}) (Z_s uᶜ_s)_r + (the second Wick term, same structure),

  with uᶜ_s = Φ_s ∘ vᶜ_s and the backward Lyapunov recursion

  𝓜ᶜ_s = P_s [diag(e_{s+1} ∘ vᶜ_{s+1}) + 𝓜ᶜ_{s+1}] P_sᵀ.

- **The sum over targets l collapsed.** The pair of continuations from the shared index a became a sandwich, as
  in JMR's ζ step.

**What it predicts and explains.**
- **Costate C3, as a theorem-shaped statement:**
  - forward costs L² n³, because each (source, target) pair needs its own Hadamard square;
  - backward costs L n⁴, because 𝓜ᶜ is one n × n state per output;
  - every mixed order pays min(L² n³, L n⁴).
- **Every age matters (N6, F6.12, F9.2).**
  - JMR's decay λ^{η(S)} requires the squared operator to contract, once its O(n)-average is removed.
  - Here P D Pᵀ = diag(Φ) W D Wᵀ diag(Φ). With W square, W D Wᵀ is a Wishart matrix at aspect ratio 1:
    Marchenko–Pastur on [0, 4 × mean] for flat D.
  - So its fluctuation equals its mean in Frobenius norm and is ≈ 3 × the mean in operator norm.
  - **λ_eff ≈ 1: no contraction with age.** This is the "no spectral gap" fact (F9.2), seen through JMR.
- **Coherent channels are cheap and amplified (N7).** The mean part (2/n) tr(D) diag(Φ²) is the slice channel. It
  is diagonal, so it costs O(n²).

**Test, as specified in v0.**
- Inside FC on bench MLP 0, split each aged source's D21 contribution into its O(n)-mean part (constant in a) and
  its Wishart fluctuation.
- Prediction: the fluctuation share stays ≥ 50 % of the old D21 energy at every age ≥ 2 (λ_eff ≈ 1).
- Kill: a fluctuation share decaying geometrically with age.

**Result (MEASURED, `tests/t2_split.py`, `tests/t2_split_mlp0.json`).**
- **No geometric decay with age: the gap is absent, as predicted.** The mean share of a source is 37–58 % at ages
  1–4 and 50–76 % at ages ≥ 6, at every target. It rises slowly with age.
- **The v0 prediction about the share is wrong.** The fluctuation is *not* ≥ 50 % of the merged old content:
  - The constant (mean) parts of different ages add *coherently*. The rms cosine between the contributions of
    different ages is 0.26–0.41 (max 0.89), and the sum of per-source energies is only 0.23–0.37 of the total energy
    at targets 8–15.
  - So the merged old D21 is **79–81 % constant in a** at targets 8, 12 and 15 (62 % at target 4).
- The old share of D21 energy is 16 %, 49 %, 66 % and 77 % at targets 4, 8, 12 and 15. This is consistent with
  region N6.

**What the constant part is.**
- (1/n) Σ_a κ3(z_a, z_a, z_b) = κ(N_t, z_b), with N_t = (1/n) Σ_a (z_a − m_a)².
- So it is the covariance of each neuron with the layer's **norm process**: the dilation sector of BRIEF item 8, and
  the collective coordinate of F2.4 and of the (2,2) spike.
- Region N7's κ4 column means are the same object one order up, κ(N_t, z_b, z_b). So the trivial-isotypic channel
  is the norm process at every order.

**Which part the final means need (`tests/t2_ablate.py`, raw on MLP 0, FC = 3.24e-8).**

| D21 used by the chain | raw |
|---|---|
| FC (all exact) | 3.24e-8 |
| old (age > 2) constant part removed | 6.8e-7 |
| old fluctuation removed | 1.0e-6 |
| all old removed | 1.77e-6 |
| old constant part from the annealed vector recursion, fluctuation exact | **3.67e-8** |
| same, with old = age > 4 | 3.59e-8 |
| young ages ≤ a exact + annealed constant channel only: a = 2 / 4 / 6 | 1.0e-6 / 4.0e-7 / 1.8e-7 |

**The annealed recursion (JMR's mean term, made computational).**
- For a fresh layer, E_W ‖x diag(P) W‖² = Σ_i x_i² P_i² c_i ≈ ‖x‖² ⟨P² c⟩, where c_i = ‖W_{i·}‖². The same holds
  for x·y.
- So a source's row scalars ‖Y_r‖², Y_r·Z_r, ‖Z_r‖² and Z_r·T_r are frozen when it turns old (exact, from the young
  tier) and multiplied by one scalar gain per layer.
- Its constant part is then (w ∘ q)ᵀ Z_s(t) + 2 (w ∘ q′)ᵀ Y_s(t) + (slice terms), divided by n: a vector times the
  transported legs. Merged over sources, that is a single vector r with r(t+1) = γ_t r(t) P_t + (entrants).
- Cost: O(n²) per layer.
- Its error against the exact constant part is 0–2 % at layers 3–7, 4–5 % at 8–11, and 7–9.5 % at 12–15. It
  grows because the frozen scalars ignore the gate–leg correlations that build up.

**Second network (MLP 1): replicates.**
- FC 1.81e-8 (region's figure reproduced).
- Old constant part removed: 8.7e-7.
- Annealed constant channel: **1.84e-8 (+2 %)**.
- Hybrid at a = 6: 1.85e-7.

**What this settles.**
- The old content separates exactly as JMR's proof does:
  - a mean term, which is cheap and carries ≈ 80 % of the energy and ≈ 40 % of the penalty;
  - a λ term (the fluctuation), with λ_eff ≈ 1, which carries ≈ 20 % of the energy and ≈ 60 % of the penalty.
- The λ term is the true cost object. It is the traceless part of the merged old tensor read through the last fresh
  layer, Σ_{ijk} T°_ijk G_ia G_ja G_kb, and it needs T° in full.
- **Cost model for the competition.** Young ages ≤ a exact at ≈ 4–7 units per source per layer, plus the constant
  channel (≈ 0 units):
  - a = 2: ≈ 0.1–0.15 B, raw 1.0e-6;
  - a = 6: ≈ 0.4 B, raw 1.8e-7.

  Neither competes. The constant channel is a free 1.8–2.1× on any age-truncated chain, and a free ingredient of any
  design that merges old content: only the traceless remainder needs carrying.

### T3. Richardson / weighted-PRG error reduction on the memory recursion (DERIVED)

**Objects in each role.**

| role | object |
|---|---|
| exact one-step operator | 𝒯_l(T) = T ×³ P_l + B_{l+1}; on a state of R CP atoms it costs ≈ 2R/n units of transport plus 2R/n units of D21 readout, cheap for R ≤ n |
| lifted system | (I − 𝒮) T = B |
| approximate inverse | any cheap chain M (a shared-subspace old tier, an age-truncated chain) |
| angle | ε₀ = ‖I − M(I − 𝒮)‖, in the readout norm priced by the fresh-weight lemma |

**Not the per-neuron drawing.** The pieces are lifted layers acting on the whole memory tensor.

**Results.**
1. **Projective compressions gain nothing (Proposition, DERIVED).**
   - Let Ã[s, l] = Π P_s Π ⋯ Π P_{l−1} Π be a Galerkin or shared-subspace compression.
   - The first-order identity Σ_i Ã[s, i−1] (Π A_i Π) Ã[i+1, l] − Σ_i Ã[s, i] Ã[i+1, l] returns Ã itself, because
     sandwiching the exact step between projections re-projects it.
   - Only the leakage (I − Π) A_i, carried *exactly* through all later steps, corrects the error. That leakage is a
     full-rank source born at every layer: carrying it costs what FC costs.
2. **Age truncation is not a preconditioner.**
   - With M = Σ_{k≤a} 𝒮^k one has I − M(I − 𝒮) = 𝒮^{a+1}, and Richardson reproduces the exact Neumann series.
   - The precision is free, but there is no saving.
3. **Where Richardson does apply, it rebuilds the public chain.**
   - Level 0 is a shared-subspace old tier. Its measured ε₀ is ≈ 0.23 of D21 at q = n/4 (region §5 ablation:
     2.3e-7 against 3.5e-8 with ε² law 3.8e-6 ε²), so α < 1 holds in the readout norm.
   - Level 1 would track the leakage exactly for a′ layers and drop it afterwards.
   - This is precisely 504aldo's architecture: young sources dense, old sources in a shared basis. Its old tier costs
     107 units (504aldo F86).

   **Prediction:** a 2-level scheme reaches ε ≈ ε₀ · (share of leakage older than a′). It cannot cost less than the
   young tier plus the q-subspace transport (≈ L q/n units per layer at depth).
   **Kill:** if any 2-level configuration gives raw ≤ 5e-8 at ≤ 10 units per layer inside FC, Richardson is a live
   line. Untested.

**Why the AKMPSV payoff does not transfer.**
- AKMPSV's cheap approximations exist because w × w products are small. Their problem is *precision*, and
  Richardson converts precision cheaply.
- Ours is *dimension*: we need only ≈ 7 % precision (N2), but at full rank. The bulk of A[s, l] is Haar-like
  (participation ratio ≈ n/(2·age), F9.1), and a Haar-like matrix has no cheap operator-norm approximant with
  ε₀ < 1. Its best rank-q approximation has relative Frobenius error √(1 − q/n).

### T3b. The Richardson question, answered: a range no-go and its escape hatches (v1.1, DERIVED)

**The question (coordinator, 21:25 UTC).** Is there an exactly applicable operator against which a cheap
approximate memory can be corrected?

**Yes, there is an exact operator.** The one-step memory map 𝒯_l(T) = T ×³ P_l + B_{l+1} applies exactly to any
CP state with R ≤ n atoms, at ≈ 2R/n units of transport plus 2R/n units of D21 readout. It applies equally in
coordinator note §2's static form: a static tensor read through the frame Γ_t. So the residual
r = B − (I − 𝒮) T̃ of any cheap chain T̃ can be formed exactly.
- For example, region's CP merge at R = n: there the residual is the per-layer merge error, a CP sum of ≤ 3R + n
  atoms.

**But the correction cannot leave the range of the approximate inverse.**
- **Lemma (range).** Let M be any linear approximate inverse (the cheap chain run on given births), and
  P_m = Σ_{i≤m} (I − M L)^i M. Then P_m B ∈ range(M) for every m. So after any number of rounds,

  ‖T − P_m B‖ ≥ dist(T, range(M)).

  Proof: each term (I − ML)^i M B is M applied to something.
- Richardson converts *precision* inside a span. It cannot create *span*.
- In AKMPSV the approximate inverse is a full-range nk × nk matrix (an LU of the lifted Laplacian), so only
  precision is missing. Our missing ingredient is span, so the AKMPSV payoff cannot transfer.

**Consequences for every cheap carrier measured so far.**

| cheap M | range | Richardson limit (≥ error of the best element of the range) |
|---|---|---|
| shared subspace q per leg (504aldo old tier, Oseledets) | legs in a q-dim subspace | the old-tier error at that q (q = n/4: 2.3e-7, region §5) |
| annealed / O(n)-averaged channel (T2 here, costate, F, D) | the trace (norm-process) sector | the traceless remainder: 1.0e-6 if dropped (T2) |
| age truncation | ages ≤ a | the old content: 8.5e-7 at a = 4 (region N6) |
| CP merge at R atoms | rank-R symmetric tensors | rank-limited: 2.0e-7 at R = n (region §8) |
| Perron / spike carriers (team D) | spike + rank 8 | 7.0e-7 (team D) |

For the trace sector the obstruction is sharper than the lemma. The traceless remainder T° is mapped by the
annealed M to ≈ 0, because by Schur orthogonality (§2.2) traceless content feeds the trace channel only
incoherently. So T° lies in the kernel of M in expectation, and the Richardson iterates never see it.

**The no-go, stated for a class.**
- Class: any estimator whose memory is a cheap linear chain M, i.e. (i) a projection, (ii) an equivariant (annealed)
  average, or (iii) a bounded-rank CP or Tucker state, or a sum of these, corrected by any number of exact one-step
  applications.
- Its error is ≥ the distance from the true memory to range(M).
- At n = 1024 the remainder T° that must be in that range is ≈ 20 % of old D21 energy, costs ≈ 1.0e-6 if missed,
  and has effective rank ∝ n per age (atom-complete; team D estimates ≈ n²/12 in total).
- A range that contains it to the needed 5–10 % costs ≥ the 480-unit bilinear bound (region §6).

**So the leaders cannot be running cheap-chain-plus-correction on this representation.**

**Escape hatches the no-go does not cover.**

E1. **An implicit full-range M.**
- The memory is never stored. It is the exact sum Σ_s B_s ×³ P_{s→t}, whose operators (the weights) are known.
  Only its readout is approximated, by an unbiased randomised estimator, and the estimator's noise is corrected by
  averaging.
- The range is full. The obstruction is variance.
- Measured closures:
  - sampling D21 directly needs ≥ 15–60 % noise per layer at budget (F6.10);
  - parity designs give 1.0× (T4);
  - random sketches of an atom-complete sum have relative noise ≈ M/√k for M atoms and k sketches (§3, v0
    reasoning).
- Open only if some structured design kills the *even-chaos* noise. None is known.

E2. **A nonlinear readout of a few transported matrices.**
- The b-mode / Tucker no-go (T6) covers readouts that are *linear* in the carried matrices.
- A readout that is *bilinear* covers rank up to n²: for example D21° ≈ Σ_{k≤K} (X_k ∘ Y_k) with X_k, Y_k
  transported by sandwiches (2 units each), or Hadamard products against the exact covariance.
- It is not closed under transport, since (WᵀXW) ∘ (WᵀYW) ≠ Wᵀ(X ∘ Y)W, so it would need a per-layer refit. It is
  also not ruled out by any measurement here.
- **This is the one class compatible with "the leaders transport ≈ 7 symmetric matrices per layer" that the
  campaign has not tested.**
- Test: fit D21°(t) (the a-traceless old fluctuation, available in `tests/t2_ablate.py`'s stash) by Σ_{k≤K}
  X_k ∘ Y_k, with X_k, Y_k drawn from the span of {C(z_t), C(a_t), Γ_tᵀΓ_t, their one-step images}. Measure the
  readout error against K.
- **Kill:** > 50 % residual at K = 8.

E3. **A different carried object whose exact one-step operator is an n²-state linear map.**
- An example would be the non-Gaussian covariance correction ΔC itself, with old content entering only through a
  closed linear response. Coordinator note §3, D21 as a cross-covariance inside the second chaos, is the natural
  route.
- No such identity is known. Region §3's fact that ΔC is 98 % incoherent and needs D21 at full rank is evidence
  against a closed n² recursion. It is not a proof.

**One-line answer.** The exact operator exists and is cheap. Richardson against it is sound, but it is confined to
the range of the cheap chain, and every cheap range measured excludes the quenched traceless remainder. What the
leaders do must lie in E1 (implicit state, variance-corrected), E2 (a nonlinear readout of a few transported
matrices, the untested one), or E3 (a different closed object).

### T4. Structured designs that fool odd chaos exactly (JMR Thm 1.7, the Z₂ case), MEASURED

**Objects in each role.**

| role | object |
|---|---|
| pieces | input samples |
| structure | the Z₂ group {x, −x}, a Cayley design that fools every odd function of x exactly (JMR: "every odd function is perfectly fooled by such a structured expander") |
| target | the sampled D21. Its signal needs at least one even factor, since E[odd³] = 0. Its leading noise term O_a² O_b is pure odd and is cancelled exactly by the pair |

**Not the per-neuron drawing.** The design acts on the input law, i.e. the Wiener-chaos parity of every layer at
once.

**Test** (`tests/d21_antithetic.py`).
- Bench w1024_d16, MLP 0. 32,768 forward passes per design, in two independent halves.
- Relative Frobenius noise of the off-diagonal D21(z_l), plain against antithetic at an equal number of forward
  passes, and the even share η of the centred variance.

| layer z_l | 1 | 3 | 5 | 8 | 10 | 12 | 15 | 16 |
|---|---|---|---|---|---|---|---|---|
| plain rel. noise | (no signal) | 1.19 | 0.68 | 0.40 | 0.29 | 0.25 | 0.19 | 0.18 |
| antithetic rel. noise | (no signal) | 1.17 | 0.68 | 0.40 | 0.29 | 0.25 | 0.19 | 0.18 |
| variance gain | – | 1.04 | 1.00 | 1.00 | 0.99 | 0.98 | 0.98 | 0.98 |
| η (even share) | 0.000 | 0.33 | 0.38 | 0.42 | 0.43 | 0.44 | 0.45 | 0.45 |

**Readings.**
- z₁ = xW₁ is exactly Gaussian, so D21 = 0 there and the antithetic estimate is exactly 0 up to the centring.
- From the first ReLU on, a third of the centred variance is already even and its share rises to 45 %. The odd
  noise is no longer dominant, and parity designs buy nothing.
- The plain noise reproduces F6.10 (0.28 at layer 10, 0.20 at layer 15).

**Prediction checked:** gain ≈ 1/(1 − (1 − η)³ × (odd-noise share)) ≈ 1 at η ≈ 0.4. **Killed**, as a carrier of
memory by sampling.
- MLP 1 (v1) agrees within noise: gain 1.04, 1.02, 1.02, 1.01, 1.00, 0.98, 0.99 at layers 3–16, with η = 0.33–0.45.
- Higher designs (OU replica designs isolating chaos levels, Kerdock-type spherical designs) face the same
  obstacle: the even chaos at depth is spread over degrees ≥ 4 (Oishi's chaos spectrum: effective degree 8–14).

### T5. INW's "short seed": old content = a few collective coordinates + a fooled remainder? (negative, from existing measurements)

**Objects in each role.**

| role | object |
|---|---|
| seed | the mean or Perron direction (F2.2) and the trace vectors of §2.2 |
| remainder | everything else in T_l |
| claim to test | the remainder's effect on the final means is fooled by the fresh future |

**Verdict.**
- The fresh-weight lemma prices the remainder at K_off ‖δC‖²_F. The measured ΔC is 98 % incoherent, and dropping
  ages > 4 costs 25× the bar (region §3, N6).
- A shared Oseledets subspace of q = n/4 directions, which contains the seed, is 6× worse (region §5).
- **The effective width of our computation is the Frobenius rank of the memory, ∝ n (F9.3), not a short seed.**
- INW's mechanism needs width ≪ 2^{seed}. Our "width" is the state the future must see, and that is n² … n³.
- So INW does not transfer. What does transfer is its *structural* lesson, already in T1: the future is fresh, so the
  past matters only through the Frobenius norm of its state.

### T6. Where a payoff could still come from (SYNTHESIS, flagged for the final version)

Every pseudorandomness theorem above wins through a **gap** (λ < 1), a **small width** (INW), or **small facets
in log** (Bentkus–Nazarov). Our transport has none of these for the incoherent half of the memory. The JMR lens
says the leaders' "zero-cost memory" cannot be a pseudorandom carrier; it must be an **exact algebraic
linearisation**.

The natural candidate is the **lift** that AKMPSV and JMR both use (tensor-lift by the cycle, ρ ⊗ ρ):
- The Hadamard square that breaks linearity, (Y A)²_{ra}, is the diagonal of Aᵀ (y_r y_rᵀ) A.
- That is linear transport on Sym²(ℝⁿ), the space of symmetric matrices.
- The "covariance-response modes" hypothesis for the leaders (504aldo §6, F77: "transport symmetric n × n
  matrices, ≈ 7 sandwiches per layer") is the statement that the Sym²-lift of the old content has small *matrix*
  rank in the b-mode: T(·, ·, b) ≈ Σ_{m ≤ M} Ξ_m g_m(b) with M ≈ 5–10 symmetric Ξ_m transported by sandwiches.

The test is the b-mode (matrix-valued) rank of the merged old tensor at n = 1024, to 7 % in the readout norm.
- If M ≲ 10, this is the leaders' carrier.
- If M ∝ n (as "atom-complete" suggests for the legs), the carrier is something else.

**v1: closed by existing measurements, no new run needed.**
- A b-mode truncation to M modes gives D21 ↦ D21 Π_M, so its readout error is at least the best rank-M error of
  D21 itself.
- Truncating D21 to rank 64 costs 2.5e-7, and rank 256 costs 8.2e-8 (F6.11, est). The part of D21 that the current
  layer does not produce needs rank 128–256 for 60–98 % of its energy.
- The merged old tensor is symmetric, so every mode has the same multilinear rank. Symmetric Tucker of the merged
  old tier at R = 512 is 14× worse (504aldo F46), and the Tucker core costs r³/n² units (F88).
- **So the Sym² lift with M ≈ 5–10 modes is ruled out in either orientation.**

What survives of T6 is item 6 of the summary: the leaders' "≈ 7 symmetric sandwiches per layer" can be the coherent
(norm-process) channel of the old content at O(n²) plus a covariance-shaped κ4 core. The incoherent fluctuation
(≈ 1.0e-6 if dropped) would still have to be carried somehow.

## 4. Tests run

| test | file | network(s) | result |
|---|---|---|---|
| T4 antithetic vs plain sampled D21 noise, n = 1024 | `tests/d21_antithetic.py` | bench w1024_d16 MLPs 0, 1 (N = 32,768) | gain 0.98–1.04× at layers 3–16; η = 0.33–0.45 |
| T2 mean / fluctuation split of old D21 by age, cross-age coherence | `tests/t2_split.py` (on `tests/fc_hooked.py`, region fc.py with three hooks) | MLP 0 | no age decay; merged old content 79–81 % constant in a; ages add coherently |
| T2 ablations: which part the final means need | `tests/t2_ablate.py` | MLP 0 | constant part 6.8e-7, fluctuation 1.0e-6 if dropped (FC 3.24e-8); MLP 1: constant part 8.7e-7 (FC 1.81e-8) |
| T2 annealed (JMR mean-term) vector recursion for the constant channel | `tests/t2_annealed.py` | MLPs 0, 1 | 3.67e-8 (+13 %), 1.84e-8 (+2 %) at O(n²); hybrids a = 2/4/6: 1.0e-6 / 4.0e-7 / 1.8e-7 (MLP 1, a = 6: 1.85e-7) |

All results are in `tests/t2_ablate_results.jsonl`.

## 5. Honest assessment

| kind | items |
|---|---|
| THEOREM (read in source) | JMR Thm 3.1, Thm 1.7–1.9, Prop 5.5; AKMPSV Lemma 6.2, Thm 1.3, the LU / squaring structure |
| THEOREM (standard, from memory) | MOO, HKM, OST, Nazarov, INW, Gowers; arXiv ids given where known |
| DERIVED | the Schur-orthogonality reading of the fresh-weight lemma (§2.2); the costate collapse with the backward Lyapunov recursion (T2, mean-field gates only); the projective fixed-point statement (T3) |
| MEASURED | T4 at n = 1024 on two networks; T2 split on MLP 0; ablations and annealed channel on MLPs 0 and 1 |
| CONJECTURE | the noncommutative ignore-first-step theorem (§2.3) |
| SYNTHESIS | "the payoff must be an exact linearisation (Sym² lift)" (T6) |

**Risks.**
- The costate derivation assumes J_ab ≈ Φ_aΦ_b. The true bivariate kernel derivative adds correlation terms, which
  change the constant but not the scaling.
- The T3 proposition is about first-order error reduction with projective segments. Non-projective cheap
  approximations with ε₀ < 1 are not excluded by proof, only by the Haar-bulk argument.
- The main negative (λ_eff ≈ 1, width ∝ n) is consistent with every established fact, so this report mostly
  *explains*; it does not *overturn*.
- The T2 measurements are inside FC, whose own D21 has 2.5–10 % error (region §5). The split is of FC's memory, not
  of the exact law. The constant-part share is measured on one network, the penalties and the annealed channel on two.
- "Ages add coherently" (cosines up to 0.89) seems to contradict F8.3 (old content nearly orthogonal to other ages).
  It probably does not: the coherence is plausibly all in the constant (norm-process) channel, and F8.3 was measured
  on content beyond it. This was not checked here.

**What the next session should do.**
1. Replicate T2 on MLPs 2–5.
2. Fix the drift of the annealed channel at layers 12–15 with one exact rank-one correction per layer: the
   gate–leg correlation of the Perron direction.
3. Attack only the traceless remainder T° of the merged old tensor, now known to be ≈ 20 % of old energy and
   ≈ 1.0e-6 of penalty. Any merging scheme (region §6 CP merge) should be fitted to T°, not to T. Its reduced norm
   loosens the merge tolerance by ≈ √(1/0.2) ≈ 2.2×.
