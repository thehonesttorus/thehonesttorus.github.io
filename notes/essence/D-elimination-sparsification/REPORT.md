# Team D: elimination, sparsification and scaling

*Status: v0 (1 Oct 2026, ≈ 21:30 UTC). Theory for §1–2, three transfers measured at n = 1024 on the bench (§4), and
three more analysed but not yet run. Every number is at n = 1024, depth 16, on `w1024_d16`, inside the region stream's
FC estimator (raw 3.24e-8 on MLP 0, 1.81e-8 on MLP 1, reproduced here). Code: `fcs.py` (FC plus unbiased atom
resampling, atom-Gram diagnostics and oracle projections), `run_d.py`, `probe_coh.py`. Results: `results/*.json`.*

## 0. The answer in five lines

1. **Every theorem in this domain controls an operator norm.** Kyng–Sachdeva, Spielman–Srivastava, BSS, MSS/Kadison–Singer
   and matrix concentration all bound the worst test direction (Loewner order, operator norm). Their computational payoff
   is the gap between the operator norm and the Frobenius norm of a spread-out error.
2. **Our problem prices errors in the tracial L² norm.** The fresh-weight lemma (region §1) makes the next layer a
   consumer that tests along *random* directions, i.e. the normalised trace τ on M_n. In L²(τ) the gap above is zero,
   so sparsification gains only what the atoms' own coherence γ = ‖Σ A_r‖²_F / Σ ‖A_r‖²_F provides (Proposition D1–D2).
3. **The variance law, derived before running and then measured.** Unbiased Poisson resampling of the old memory atoms
   (Kyng–Sachdeva's move, the sample becomes the state and is transported) has relative variance Σ(1/q_r − 1)‖A_r‖² / ‖S‖²,
   ≈ (1/f − 1)/γ for even norms. Predicted extra raw at f = 0.75: 2.0e-8; measured 2.1e-8 (MLP 0). At f = 0.5:
   predicted 1.3e-7, measured 2.4e-7. Reaching the bar needs f ≳ 0.85, which saves ≤ 15 % of FC's cost. **Killed.**
4. **New structure found on the way.** The old memory is *coherent*, not atom-incoherent: γ grows with depth from 1.2
   to 8.8; adjacent old ages have D21 cosines 0.70–0.88; at t = 14, 75 % of the old D21's energy lies along one Perron
   direction of the propagator in its second leg. But the coherent part is not sufficient: an oracle keeping only the
   Perron component of the old memory scores 1.5e-6, and the best rank-8 second-leg oracle 7.0e-7. The memory is a
   **spike plus an atom-complete bulk**, and the bulk (≈ 8–10 % of the old energy) costs ≈ 20× the bar.
5. **What this means for the programme.** The leaders' "zero-cost memory" cannot be a sampled, sparsified or spiked
   representation of the FC atoms. Whatever they do changes the *object*, not its representation. §3.6 gives the
   one sparsification-flavoured transfer left standing.

## 1. Instances

### 1.1 Kyng–Sachdeva approximate Gaussian elimination (arXiv 1605.02353)

**Statement.** For a Laplacian L of a graph with m edges, a random lower-triangular U with O(m log³ n) nonzeros can be
computed in O(m log³ n) time with L ≈_{1/2} UUᵀ (spectral, w.h.p.). Preconditioned CG then solves Lx = b to accuracy ε
in O(m log³ n log 1/ε) time. Extended to connection Laplacians and SDD in the same paper. Kyng's thesis and later work
(Gao–Kyng–Spielman, arXiv 2303.00709, "robust and practical") give practical constants; the method is the basis of
the Laplacians.jl `approxchol` solver.

**The pieces.** Vertices. Each edge is first split into ρ = Θ(log² n / ε²) parallel copies, so that every multi-edge
has leverage w_e R_e ≤ 1/ρ.

**Down step.** Eliminate vertex v. The Schur complement L − L_{:v} L_{v:} / L_vv is exactly the precision matrix of the
Gaussian free field marginalised over v. It removes the star at v and adds a weighted clique on N(v),
C_v = Σ_{i<j} (w_vi w_vj / w_v) b_ij b_ijᵀ. That is d² edges from d.

**Up step.** Replace the clique by d random edges. For each multi-edge (v, i), pick j with probability w_vj / w_v and
add edge (i, j) with weight w_vi w_vj / (w_vi + w_vj). The sample is unbiased, E[sample] = C_v, and it becomes part of the
matrix that later eliminations act on.

**Certificate (the angle).** Leverage. The normalised increments X_k = L^{-1/2}(sample_k − C_v)L^{-1/2} form a matrix
martingale. Each has ‖X_k‖ ≤ 1/ρ, because every sampled multi-edge keeps leverage ≤ 1/ρ in the *current* Schur
complement. This is the key lemma: leverage scores never increase under elimination. Their predictable quadratic
variation is ≤ (1/ρ)·I per round.

**Global conclusion.** Matrix Freedman (Tropp, arXiv 1101.3039) gives ‖Σ X_k‖ ≤ ε w.h.p. The errors of all
eliminations add as martingale increments in the operator norm, not as a sum of norms.

**Where the payoff comes from, exactly.** Two facts together:
- (a) the sample has d edges, not d², so fill-in stays O(m log² n);
- (b) the error metric is relative spectral.

A sampled clique has relative *Frobenius* error of order 1, since d edges stand in for d²/2. Its relative *spectral*
error is only 1/√ρ, because every sampled edge has small leverage. Without (b) the method fails, as §2 shows.

### 1.2 Spectral sparsification: Spielman–Srivastava, BSS, MSS and Kadison–Singer

- **Spielman–Srivastava** (arXiv 0803.0929). Sample edges with probability ∝ w_e R_e and reweight. With O(n log n / ε²)
  samples, (1 − ε)L ≤ L̃ ≤ (1 + ε)L. The proof is Rudelson / Ahlswede–Winter matrix concentration for
  Σ_e (w_e R_e)^{-1}·L^{-1/2} b_e b_eᵀ L^{-1/2}, a sum of isotropic rank-ones.
- **Batson–Spielman–Srivastava** (arXiv 0808.0163). The deterministic version, with dn rank-one terms and condition
  number (√d + 1)² / (√d − 1)².
  - **Up step.** Add one rank-one term s·v vᵀ and move the barriers u ← u + δ_U, ℓ ← ℓ + δ_L.
  - **Certificate.** The potentials Φ^u(A) = tr(uI − A)^{-1} and Φ_ℓ(A) = tr(A − ℓI)^{-1} do not increase. A good vector
    exists by averaging: Σ_v of the upper barrier's requirement is ≤ Σ_v of the lower barrier's allowance. This is a
    "local check implies a global spectral bound" argument run one rank-one step at a time.
- **Marcus–Spielman–Srivastava** (arXiv 1306.3969; Ramanujan graphs arXiv 1304.4132).
  - Interlacing families: the expected characteristic polynomial of a sum of independent random rank-ones is real-rooted
    (the mixed characteristic polynomial), and some outcome has largest root ≤ the expected polynomial's largest root.
  - The multivariate barrier bounds that root by (1 + √ε)².
  - Consequence (Weaver KS₂): vectors with Σ v vᵀ = I and ‖v‖² ≤ ε split into two halves, each with norm
    ≤ (1/√2 + √ε)².
- **Kadison–Singer as a C*-statement.** Every pure state of the diagonal masa ℓ^∞ ⊂ B(ℓ²) has a unique state
  extension to B(ℓ²). Equivalently (Anderson paving), for every ε there is an r such that every zero-diagonal T can be
  paved: diagonal projections Q_1 … Q_r with Σ Q_j = I and ‖Q_j T Q_j‖ ≤ ε ‖T‖.
  - The down step is the conditional expectation E: B(ℓ²) → ℓ^∞ (forget the off-diagonal piece).
  - The theorem says the forgotten piece is locally invisible on some coarse partition. It is a pure operator-norm
    statement.
  - **In L²(τ) paving is trivial**: a random r-partition gives Σ_j ‖Q_j T Q_j‖²_F ≈ ‖T‖²_F / r.

### 1.3 Operator scaling: Gurvits, GGOW, Kwok–Lau–Ramachandran

**Statement.** Take a completely positive map T(X) = Σ A_i X A_i† on M_n. Scaling alternates two normalisations:
- A_i ← T(I)^{-1/2} A_i;
- A_i ← A_i T*(I)^{-1/2}.

These are the two "forget" steps onto the affine sets {T(I) = I} and {T*(I) = I}. Sinkhorn is the diagonal (commutative)
case.

**Certificate.** The capacity cap(T) = inf_{X≻0} det T(X) / det X. It is > 0 iff T is ε-scalable for all ε, iff the
symbolic matrix Σ x_i A_i has full noncommutative rank (Gurvits 2004).

**Global conclusion and payoff.** Garg–Gurvits–Oliveira–Wigderson (arXiv 1511.03730) show poly(n, b, 1/ε) iterations,
using a lower bound on the capacity from invariant theory (Derksen–Makam degree bounds). This gives a deterministic
polynomial-time algorithm for noncommutative rank (Edmonds' problem over the free skew field). The commutative analogue
(PIT for symbolic determinants) is open and tied to circuit lower bounds.

**The angle.** Kwok–Lau–Ramachandran ("Spectral analysis of matrix scaling and operator scaling", FOCS 2019; arXiv id
1904.03213, from memory): if T is ε-close to doubly balanced and its bipartite operator has spectral gap λ (second
singular value ≤ 1 − λ), scaling converges *linearly*, in O(log(1/ε)/λ) iterations, and the scaling solution is close to
the input. The angle is exactly the gap: the cosine between the two normalisation subspaces.

### 1.4 Shorter: intrinsic freeness and expander hierarchies

- **Bandeira–Boedihardjo–van Handel** (arXiv 2108.06312). Write X = A₀ + Σ g_i A_i and σ(X)² = ‖E(X − A₀)²‖, and let
  v(X)² = ‖Cov(X)‖ (the covariance as an operator on n² entries). Then the spectrum of X is close to that of the free
  model X_free = A₀ ⊗ 1 + Σ A_i ⊗ s_i, with Hausdorff distance ≲ v^{1/2} σ^{1/2} (log n)^{3/4}. In particular
  ‖X‖ ≤ ‖X_free‖ + C v^{1/2} σ^{1/2} (log n)^{3/4}. A 2025 survey is arXiv 2510.01021.
  - This is again a statement about the spectrum (operator norm and spectral distribution) of a quenched matrix. It gives
    deterministic equivalents for *laws*, not for the matrix entries.
- **Expander hierarchies and almost-linear max-flow** (Chen–Kyng–Liu–Peng–Probst Gutenberg–Sachdeva, arXiv 2203.00671;
  expander decompositions, Saranurak–Wang arXiv 1812.08958).
  - The global problem is solved by interior-point steps.
  - Each step is a min-ratio cycle problem on a dynamically maintained, recursively sparsified hierarchy: expanders at
    each level are "forgotten" into a core graph with spectral or ℓ₁ embedding guarantees.
  - The payoff is m^{1+o(1)} total time. It rests on (i) a metric in which compressing an expander costs only a
    subpolynomial distortion, and (ii) a potential argument that tolerates the approximation.

## 2. The noncommutative core: L^∞ versus the tracial L²

**What is known.** The NC form of the payoff is matrix (noncommutative) martingale concentration:
- Ahlswede–Winter; Tropp's matrix Freedman (arXiv 1101.3039);
- Pisier–Xu NC Burkholder–Gundy; Junge–Zeng;
- the MSS root bounds;
- BBvH's intrinsic freeness for the sharp constant.

All of these bound ‖·‖_{L^∞(M_n)} = operator norm. Kadison–Singer is the extreme case: an operator-norm statement about
the conditional expectation onto a masa.

**Proposition D1 (the price is tracial; proved here from the fresh-weight lemma).** Let δ be an error in a layer-l
symmetric object that is independent of W_{l+1}, with w_c the fresh columns. Then
E_W (w_cᵀ δ w_c − (2/n) tr δ)² = 2 (2/n)² ‖δ‖²_F. Write τ = tr/n for the normalised trace. Up to the transfer
coefficient, the final-MSE price is therefore ∝ ‖δ‖²_{L²(M_n, τ)}. A fresh Gaussian consumer is the tracial state: it
averages over test directions instead of maximising over them.

*Consequence.* A (1 ± ε) spectral guarantee on a PSD sum S gives Frobenius error ≤ ε ‖S‖_F (write δ = S^{1/2} E S^{1/2}
with ‖E‖ ≤ ε). It gives no better than that, because sampling errors are spread out, with stable rank ≈ n. Spectral
sparsifiers therefore need Θ(n / ε²) atoms for relative Frobenius error ε. At our tolerance ε ≈ 5–7 % (region N2)
that is ≈ 200–400 n atoms. The memory has only N = t·n ≤ 15 n atoms. **Spectral sparsification (SS, BSS, MSS, KS) cannot
pay at this tolerance under this pricing: a theorem, given D1.**

**Proposition D2 (the variance law of unbiased atom sampling; exact).** Let S = Σ_{r=1}^N A_r, with independent
inclusions q_r and reweighting 1/q_r (Poisson sampling). Then:
- E Ŝ = S;
- E ‖Ŝ − S‖²_F = Σ_r (1/q_r − 1) ‖A_r‖²_F;
- the optimum at expected size k is q_r = min(1, c ‖A_r‖_F);
- without capping, the relative variance is (N_eff / k − 1)/γ, where N_eff = (Σ ‖A_r‖)² / Σ ‖A_r‖² and
  γ = ‖S‖² / Σ ‖A_r‖²;
- for even norms this is (1/f − 1)/γ with f = k/N.

Matrix concentration adds nothing to the second moment: under L²(τ) pricing the second moment *is* the price.

*Deterministic version (subset selection).* Let G_rr′ = ⟨A_r, A_r′⟩_F be the atom Gram. Any reweighting supported on k
atoms has error ‖A(1 − c)‖² ≥ λ_min(G)(N − k), since at least N − k entries of 1 − c equal 1. BSS-type deterministic
selection therefore cannot beat the spectrum of the atom Gram.

**Proposition D3 (positivity is the only route to N-independence).** Suppose the atoms are PSD rank-ones a_r a_rᵀ and
q_r ∝ ‖a_r‖². Then Σ ‖A_r‖_F = tr S, and the relative variance is ≤ PR(S)/k with PR(S) = (tr S)² / ‖S‖²_F (the
participation ratio). This bound is independent of N.

Our atoms are not PSD. A_r = w2_r[(y_r ∘ y_r) z_rᵀ + 2 (y_r ∘ z_r) y_rᵀ] + slice terms:
- the y ⊗ y legs are PSD and w2 > 0;
- but the third leg z_r is sign-indefinite.

Only the left leg u_r = y_r ∘ y_r is entrywise positive, which makes the left-leg Gram (u_r·u_r′) positive. The role of
PR is then played by N_eff / γ, which must be measured (§4.1).

**The precise next statement (conjecture D4, "tracial obstruction").** Suppose an estimator of the final means
realises the memory functional as a linear image of a random or deterministic sub-collection of the FC atoms, possibly
re-weighted and transported. Then its excess raw MSE is
≥ c · Σ_t K(t) · min_{|R| ≤ k_t} ‖Σ_{r ∉ R} P_R^⊥ A_r(t)‖²_F, where P_R^⊥ projects onto the orthogonal complement of the
span of the kept atoms. There is no logarithmic or concentration gain, because the downstream consumer is the tracial
state.

The NC content is this: matrix-martingale methods transfer to our problem only after a change of state from τ to a
vector state concentrated on few directions. That would require a consumer whose test directions are quenched and low-rank.
At depth K_top < K_off (region §2), so the top direction is *cheaper*, not dearer. The one place an L^∞-type state
appears is the bethe (2,2) spike along the Perron/mean direction. It is a single direction, carried at O(n²) already
(region N7).

**Synthesis (not a theorem): operator scaling as the memory clock.** KLR's linear convergence needs a spectral gap of the
bipartite operator. Our alternating pair is the gate algebra (diagonal projections) and the weight algebra (fresh
Gaussian matrices), whose product is the gated propagator. The gated propagator has *no* spectral gap (F9.2) and a
participation ratio ≈ n/(2·age) (F9.1), so the "forget" iteration contracts only polynomially. That is exactly why
content of every age matters (N6, F8.2): the absence of a gap is the memory. The Perron outlier (gain 3–8× the bulk, F2.2)
is the one gapped direction, and it is the coherent part measured in §4.2.

## 3. Transfers

For each: the objects in each role, why it is not the per-neuron drawing, prediction, test, cost, kill.

### 3.1 Kyng–Sachdeva resampling of the memory (measured, killed)

- **Pieces.** The FC memory atoms (s, r): one per (birth layer, birth neuron). Each atom is the rank-≤4 n × n matrix
  A_r(t) built from row r of Y_s(t) = Cov(g_s, z_t), Z_s(t) and the slice leg ΔZ.
  - The atoms are elements of the *integrated-out* field's chaos: the "fill-in" created by marginalising layer s, as in a
    Schur complement. They are not neurons of the drawn graph: a neuron index r labels an atom only through the Wick
    pairing (the hafnian edge) that created it.
- **Down step.** Integrating out the earlier layer creates dense fill (D21 = Σ over all ages and atoms).
- **Up step.** At age 3, keep atom r with probability q_r ∝ ‖A_r‖_F and reweight by 1/q_r. Thereafter transport only
  the kept rows: Y[R], Z[R] and T[R] all transport by right-multiplication with diag(P) W.
  - Cost per (source, target) pair scales by f.
  - The sample becomes the state, exactly as in KS.
- **Angle.** γ and N_eff (Proposition D2).
- **Prediction.** Δraw ≈ 3.8e-6 × Σ_t ŵ_t · relvar_t(D21) (the ε² law, with K_off weights ŵ).
- **Test.** §4.1.
- **Cost at f.** ≈ 840 · (28 + 92 f)/120 units (≈ 28 young pairs exact, ≈ 92 old pairs sampled). f = 0.85 gives
  ≈ 750 u.
- **Kill.** Raw > 1.5 × FC at the f that halves the cost. Killed: f = 0.5 costs 8× in raw.
- **Explains.** Why deterministic atom pruning fails worse than random (region §5: pruning to 50 % gives 1.3e-6;
  unbiased 50 % gives 2.7e-7). The dropped atoms add coherently (γ > 1), so their bias costs γ × their energy, while the
  variance costs only (1/f − 1) × their energy.

### 3.2 BSS / interlacing selection of O(n) representative atoms (theorem-level kill)

- **Pieces.** As in 3.1.
- **Up step.** Deterministic rank-one additions under moving barriers.
- **Kill.** By Propositions D1 and D2: with L²(τ) pricing at ε ≈ 5 %, any spectral guarantee needs ≥ n/ε² atoms
  ≫ N = t n, and any subset selection is bounded below by λ_min of the atom Gram.
- **No run needed.** The §4.1 numbers bound what any selection rule can achieve. The per-target optimal Poisson
  inclusion (with hindsight) is only 1.3–1.4× better than selection at age 3.

### 3.3 Spike plus bulk: carry the coherent Perron part and drop or sample the rest (oracle measured, killed)

- **Objects.** The second leg of the old memory, projected on the Perron right singular direction v_t of the oldest
  propagator: D_old ≈ h vᵀ + R.
  - If h could be carried cheaply, this would be the "one gapped direction" of the synthesis in §2.
  - h vᵀ is cheap: h = diag(Σ_s Y_sᵀ diag(w2 ∘ Z_s v) Y_s) + three analogues. Because the left singular vectors of
    Z_s(t) converge (Oseledets), Z_s(t) v_t ≈ λ(t) c_s with a common scalar. Then Σ_s Y_sᵀ diag(β_s) Y_s is a single
    n × n quadratic form, transported by one sandwich per layer: ≈ 4 sandwiches ≈ 8 u per layer for *all* ages.
- **Not the neuron drawing.** The carrier is a quadratic form on the layer space indexed by the Oseledets direction, not
  by neurons.
- **Prediction.** If the Perron part sufficed, raw ≈ FC at ≈ 1/8 of the cost.
- **Oracle test.** Replace D_old by D_old V Vᵀ inside FC, with V the top-q directions (§4.2).
- **Result.** q = 1 Perron gives 1.5e-6. The best q = 8 (the singular vectors of D_old itself, a strict upper bound for any
  q = 8 b-leg carrier) gives 7.0e-7. **Killed for q ≤ 8.** The bulk is atom-complete, consistent with "shared q = n/4 is
  6× worse" (region §5).

### 3.4 Operator scaling: is a moment closure an alternating-normalisation fixed point? (analysed, not run)

- **Objects.**
  - First normalisation: the forward push-forward through a layer (exact).
  - Second normalisation: the M-projection onto the closure family (match moments; forget higher cumulants).
  - Iterated, this is assumed-density filtering. An EP / Sinkhorn-type two-sided fixed point would add a *backward*
    normalisation from the readout.
- **Analysis.** Our problem is a pure forward marginal: there is no evidence at the output. The backward message is the
  constant function, so the second normalisation is the identity and EP reduces to ADF. There is no fixed-point gain.
- **The one non-trivial backward object** is the costate (sensitivity) of the final mean. It could reweight where
  accuracy is spent. But K_coh ≈ K_off and K_top ≤ K_off (region §2): the costate weight is flat in direction, so the
  reweighting gains nothing in Frobenius.
- **Kill (stated in advance).** If an EP-style two-sided closure on 64 × 16 does not beat ADF by > 10 % at equal cost, drop
  it. Not run: the argument above predicts no gain.

### 3.5 Intrinsic freeness for the non-quenched part (analysed, not run)

- **Objects.** X = the incoherent (j ≠ k) part of the transport Σ W_ai W_bj W_bk κ3(a)_ijk, an order-3 Gaussian chaos in
  the fresh weights, conditionally on the past.
- **What a free model gives.** A deterministic equivalent for its *spectrum and norm*. The memory's ensemble mean is
  zero (F8.5), so the deterministic equivalent of the memory itself is 0 and its omission costs its full energy. That is
  the N6 failure, about 25× the bar.
- **Positive use: pricing only.**
  - The transfer coefficients K(l) are functionals of free multiplicative convolutions of the arcsine-gated Gaussian
    products, so they can be predicted without the region's finite-difference runs.
  - So can the growth of ‖D_old‖_F.
- **Cheapest test.** Compare a free-model prediction of K_off(l) with the region's table (MLP 0). Kill if it is off by
  > 30 %. Not run in v0.

### 3.6 Expander hierarchy: coarsen the oldest ages (proposed; the transfer left standing)

- **Objects.** The sources as nodes of an age hierarchy.
- **Measured.** Old sources are nearly collinear at the D21 level among the *oldest*: at t = 14 the cosines between
  sources of birth layers 0–1, 1–2 and 2–3 are 0.86, 0.88 and 0.85, falling to 0.2–0.4 for ages 3–6. Their legs have
  collapsed onto the common Oseledets directions.
- **Proposal.** As in an expander hierarchy, merge sources whose D21 contributions have cosine > 0.85 into one carrier,
  by CP re-fit in the shared leg subspace: the region's §6 "CP merge" experiment, now with a structural reason to expect
  success at the oldest ages.
- **Cost.** Merging ages ≥ 6 into one carrier cuts the pairs from 120 to ≈ 75, i.e. ≈ 500 u. **On its own this is not
  enough for the bar.** The mid ages 3–6, which are least collinear, carry the K_off-weighted price.
- **Kill.** Raw > 5e-8 with ages ≥ 6 merged into R = n atoms.

## 4. Tests run (n = 1024, bench `w1024_d16`, FC with slices = 2 and κ4 mean-field; truth noise subtracted)

### 4.1 Atom coherence and the variance law (MLP 0; old = age ≥ 3; `results/diag_m0.json`)

| target t | old atoms N | γ_old | γ per source (min / median / max) | N_eff/N | ‖D_old‖/‖D‖ | relvar of D21 at f = 0.5 | at f = 0.25 |
|---|---|---|---|---|---|---|---|
| 3 | 1024 | 1.20 | 1.20 | 0.98 | 0.29 | 0.066 | 0.20 |
| 6 | 4096 | 3.73 | 1.37 / 1.94 / 2.08 | 0.74 | 0.60 | 0.061 | 0.22 |
| 9 | 7168 | 6.49 | 1.53 / 2.27 / 3.64 | 0.57 | 0.75 | 0.033 | 0.15 |
| 12 | 10240 | 7.13 | 1.52 / 2.37 / 5.76 | 0.46 | 0.82 | 0.019 | 0.12 |
| 15 | 13312 | 8.81 | 1.43 / 1.97 / 9.24 | 0.40 | 0.88 | 0.010 | 0.087 |

K_off-weighted relvar(D21), and the ε²-law prediction Δraw = 3.8e-6 × relvar, against measurement:

| f | weighted relvar (age-3 selection, MLP 0) | predicted Δraw | MLP 0 (FC 3.24e-8) | MLP 1 (FC 1.81e-8) | MLP 2 (FC 3.03e-8) |
|---|---|---|---|---|---|
| 0.9 | — | — | — | 1.96e-8 (Δ 0.15e-8) | — |
| 0.75 | 0.0052 | 2.0e-8 | **5.36e-8** (Δ 2.1e-8) | 6.09e-8 (Δ 4.3e-8) | 6.11e-8 (Δ 3.1e-8) |
| 0.5 | 0.033 | 1.3e-7 | **2.72e-7** (Δ 2.4e-7) | 2.86e-7 (Δ 2.7e-7) | 2.15e-7 (Δ 1.8e-7) |
| 0.5, uniform q (no importance) | — | — | 7.14e-7 | — | — |
| 0.25 | 0.147 | 5.6e-7 | — | — | — |

Readings:
- The variance law predicts the f = 0.75 penalty to within 10 % on MLP 0 and within 2× elsewhere.
- The excess at f = 0.5 is consistent with compounding: sampled D21 errors also enter the next layer's Edgeworth slice
  corrections and the κ4 recursion.
- Staying within 1e-8 of FC needs weighted relvar ≤ 0.0026, i.e. f ≳ 0.85; f = 0.9 costs 8 % on MLP 1.
- Norm-proportional importance is worth 2.6× over uniform inclusion at f = 0.5, as D2 predicts (N_eff/N = 0.4–0.7 at depth).

### 4.2 Where the coherence lives (MLP 0; `results/coh_m0.json`, `results/proj*_m0.json`)

| t | sources | mean cosine between old sources | adjacent-age cosines (oldest first) | cos(old, young) | second leg on Perron v | first leg on v | per-source γ |
|---|---|---|---|---|---|---|---|
| 6 | 4 | 0.38 | 0.40, 0.43, 0.38 | 0.40 | 24 % | 0.1 % | 1.4–2.1 |
| 10 | 8 | 0.46 | 0.70, 0.74, 0.68, 0.53, 0.42, 0.31, 0.24 | 0.38 | 63 % | 0.04 % | 1.6–4.3 |
| 14 | 12 | 0.45 | 0.86, 0.88, 0.85, 0.70, 0.55, 0.41, 0.42, 0.29, 0.21, 0.23, 0.18 | 0.38 | 75 % | 0.1 % | 1.4–8.6 |

Oracles inside FC (old memory replaced by its projection, young sources exact):

| old memory kept as | share of ‖D_old‖² kept (t = 6 / 10 / 15) | raw, MLP 0 |
|---|---|---|
| everything (FC) | 1 | 3.24e-8 |
| D_old v vᵀ, v the Perron direction of the oldest propagator | 0.24 / 0.63 / 0.75 | 1.52e-6 |
| D_old V Vᵀ, V the top-8 right singular vectors of D_old (oracle) | 0.83 / 0.90 / 0.93 | 6.98e-7 |
| nothing older than age 2 (region N6) | 0 | 1.9e-6 |

Readings:
- cos(old, young) = 0.38, i.e. R² = 0.14, reproduces F8.1.
- **New:** old content is *not* orthogonal across ages. The oldest sources are nearly collinear, and the coherence is in
  the second leg, along the Perron direction.
- The first leg (the y ∘ y, positive leg) carries no Perron component.

## 5. Honest assessment

- **Theorem.** Propositions D1–D3 are elementary and proved here (Gaussian second moments; Poisson sampling algebra; the
  PSD trace bound). The statements about KS, SS, BSS, MSS, GGOW and BBvH are the published theorems, with arXiv ids
  1605.02353, 0803.0929, 0808.0163, 1306.3969, 1511.03730, 2108.06312 and 2203.00671. KLR's id (1904.03213) is from
  memory and not yet checked.
- **Measurement.** §4 is on MLPs 0–2 with one sampling seed each. The f = 0.75 agreement (Δ 2.1e-8 measured against
  2.0e-8 predicted) is on MLP 0 only; MLP 1's Δ is 2× the prediction.
- **Synthesis.**
  - "The tracial state is why this whole domain does not pay" (D1 plus conjecture D4).
  - "The absence of a spectral gap is the memory" (KLR reading of F9.1–F9.2).
- **Speculation.** 3.6 (the age hierarchy) and the sandwich-transport of the Perron part (3.3) as a cost device, if a
  future object makes the Perron part sufficient.
- **Risks.**
  - The Perron vector used is the top right singular vector of the oldest propagator, not the mean direction μ_z.
  - The oracle for q = 8 uses D_old's own singular vectors, so it is an upper bound for any rank-8 second-leg carrier,
    not an achievable design.
  - Conjecture D4 is stated for linear images of atom sub-collections. It does not cover estimators that change the object
    (for example by re-expanding in another chaos basis), which is where §0.5 points.
