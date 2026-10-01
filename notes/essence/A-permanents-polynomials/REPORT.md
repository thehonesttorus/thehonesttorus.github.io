# Team A: permanents, hafnians and polynomial methods

*Local-to-global essence programme, round 3 (brief: `../BRIEF.md`; mandate: `../INPUT-2026-10-01-local-to-global.md`).
v0, 1 Oct 2026. Labels: **[thm]** a published theorem (source given); **[prop]** proved here; **[meas]** measured here
at n = 1024 on the bench; **[syn]** synthesis; **[spec]** speculation. Code and results are in this directory:
`fc_hook.py` (FC with a hook), `t_atoms.py`, `t_lorentz.py`, `results/`.*

## 0. Summary

1. **The permanent world's speedups all come from one of three structures. The binding object, our memory, has none
   of them.** The three are:
   - nonnegativity, which makes Bethe, belief propagation and JSV work;
   - multiplicativity over rows, which makes the Godsil–Gutman, Clifford and quaternion variance reductions compound;
   - a zero-free disc around a solvable reference, of bounded degree or with expansion, which drives
     Barvinok, Patel–Regts and the cluster expansions.

   The memory (old third-order content) is signed and has zero ensemble mean (F8.5). It is an **additive** sum of
   T = (age)·n nearly orthogonal atoms (F8.3). Around any rank-k reference it needs k ≈ n/4 (F9.5, re-measured
   below). This is a measured classification, not a slogan: §3 gives a number for each item.
2. **Strongest instance: the diagonal-pairing (Bethe) split of the memory [meas].** Write every squared leg
   Y_ra² = (Σ_paths)² as a sum over pairs of paths (p, p′) from the birth neuron r to the target neuron a.
   - The **diagonal pairs p = p′** form the Bethe/annealed value. It is nonnegative and costs one Hadamard-squared
     chain.
   - The **off-diagonal pairs** are signed, and they are the quenched memory.
   - At n = 1024 the Bethe value has cosine **0.71 / 0.88 / 0.89** with the true old D21 at layers 4 / 10 / 14. It
     leaves a relative Frobenius error of **0.70 / 0.48 / 0.47**.
   - A rank-one mean-field of the squared legs, the trace-channel shape, does as well: **0.70 / 0.44 / 0.40**.
   - So the nonnegative sector carries about 80 % of the energy of the memory at depth. The signed remainder is
     **≈ 5–10× too large** for the 5–10 % tolerance (region §4).

   This is the precise sense in which old content is "the off-diagonal pairings that every Bethe, cover or annealed
   method sets to zero". It explains three earlier findings: the interpolation stream's covers (48× worse), its
   trace channel (it *is* the Bethe sector), and F8.3/F8.5.
3. **Sign migration, NC statement [prop].** For an *additive* sum of T rank-one third-order atoms:
   - any unbiased sketch with random signs drawn from a d-dimensional unitary representation (scalars, Clifford,
     quaternions, Haar U(d)) has variance per unit cost equal to that of scalar signs, up to a constant;
   - the Chien–Rasmussen–Sinclair advantage exists only for *products*, where per-factor reductions compound.

   **[meas]** The scalar cube-root-of-unity sketch of our memory has single-sample relative variance
   **0.9e6 / 3.8e6 / 1.1e7** at layers 4 / 10 / 14 (T = 2k / 8k / 12k atoms). It would need
   **3.6e8 – 4.4e9 samples** for 5 % error, about 10⁵–10⁶ times the budget. Dead, with the mechanism identified.
4. **Gate law, Lorentzian check [meas].** Is the homogenised gate-pattern polynomial f(t, t₀) = E ∏(gᵢtᵢ + (1−gᵢ)t₀)
   Lorentzian (log-concave)? At every layer the answer is **no**.
   - Its Hessian at 1 has **51–141 positive eigenvalues** (Lorentzian allows 1).
   - The cause is bulk, not a collective mode: the semicircle edge of the off-diagonal gate covariance exceeds the
     diagonal shift p² by a factor of 1.9 (layer 1) to 4.8 (layer 16).
   - This reproduces F4.2's spectral independence of 2–8.
   - So the Anari–Liu–Oveis Gharan–Vinzant local-to-global theorem (drop i = ∂ᵢ) does not apply to gate patterns.
5. **Barvinok around a solvable reference.**
   - Take the Perron/top-k reference. Its degree is a rank, and rank k = 256 = n/4 still leaves 0.09–0.30 error
     (rank 1: 0.91–1.0) **[meas]**, as F9.5 found.
   - Around the Bethe reference (T2), the test was: Bethe/mean-field reference at full rank, plus rank-k legs for the
     off-diagonal fluctuation only. **[meas, MLPs 0 and 1]** It gives 0.24 / 0.19–0.20 at k = 128 and 0.10–0.11 /
     0.08–0.085 at k = 256 (t = 10 / 14). Rank alone gives 0.28–0.29 / 0.22–0.24 and 0.11–0.12 / 0.086–0.094.
   - So removing the nonnegative sector buys only 10–15 %. **The signed remainder is as high-rank as the whole.**
   - T2 is killed, and F9.4/F9.5 stand.

**Verdict (v0).** The domain's machinery *explains* why the memory is expensive and measures each obstruction. It
produces no carrier.

The memory is the off-diagonal (signed) path-pairing sector. It is:
- 20 % of the energy at depth, and ≈ 45 % of the Frobenius norm;
- full-rank, with no Bethe, cover, Clifford or Lorentzian handle.

Any O(n³)-per-layer carrier must therefore exploit something *outside* this domain's three structures. The candidates
are cancellation in the readout (region §6(b)) and adjoint aggregation.

---

## 1. Instances (the process, to the level of the proofs)

### I1. Jerrum–Sinclair–Vigoda and its constants (Newman–Vardi, arXiv 2012.03367)

**The chain.**
- *Pieces:* edges of a perfect matching of K_{n,n}, with activity 1 on edges of G and λ off G.
- *Down step:* delete an edge, giving a near-perfect matching with holes (u, v).
- *Up step:* slide a hole (swap one matched edge) or re-add an edge.
- *Acceptance:* Metropolis acceptance with weight w(u,v)·λ(M).
- The state space has |Ω| = (n²+1)·n! states: the perfect matchings plus n² hole classes.

**The hole weights are cavity ratios.** The ideal weights are w*(u,v) = λ(P)/λ(N(u,v)) = per/per-with-holes. They put
equal stationary mass 1/(n²+1) on the perfect matchings and on each hole class.

**The angle, i.e. the certificate.**
- Canonical paths unwind the cycles of I ⊕ F.
- The injective encoding η_t(I,F) ≈ I ⊕ F ⊕ M is a matching with at most two holes.
- The hole weights make every hole class as heavy as P, so a hole costs nothing in the congestion bound.
- Weights within a factor 2 cost only a constant. *This is why they can be learned*: anneal λ from 1 (K_{n,n}, where
  w = n is exact) down to 1/n!, by factors 2^{−1/(2i)}. Each stage's weights are then within a factor 2 at the next
  stage, and samples re-estimate them [thm, JSV JACM 2004; BSVV, via N&V §3.3–3.5].

**Constants: where they come from (N&V eqs. 13–25, 55–59) [thm/meas there].**

| factor | bound | origin | N&V verdict |
|---|---|---|---|
| mixing per sample | τ ≤ 7ℓρ/π(P)·ln(1/δ) = 336(n⁴+n²) ln(1/δ), with ℓ ≤ n, ρ ≤ 12n, π(P) ≥ 1/(4(n²+1)) | path length × congestion × one over the mass of perfect matchings | believed loose: resampling time was relaxed to **1 step** (a factor of up to 3.4e7) with 0/640 failures at n ≤ 10; but the chain is smaller than the run there (steps/|Ω| < 1 only at n = 16), so this is untestable |
| warm-up | ≈ 336 n⁵ ln n | ln(1/π(x)) ≈ ln n! | same |
| stages | 0.4 n ln²n ≤ ℓ ≤ 7.22 n ln²n | cooling 2^{−1/(2i)} keeps the weights within a factor 2 | **structural, never relaxed** |
| samples per stage (weights) | 475(n²+1) ln(24ℓ(n²+1)) | **n² cavity ratios**, each to a factor 2^{1/4}, with a Chernoff union bound | partly loose (they fix the union-bound over-count); relaxing by 256× **fails** at n = 8, 10 |
| samples (count) | ≈ (10800/ln 2) n ln²n / ε² | variance (1+9/S)^ℓ − 1 across ℓ telescoping ratios | intrinsic to the telescoping product |
| total | (1,276,800/ln 2) n⁷ ln⁴ n | | crossover with Ryser at n ≈ 68 (1.3e22 steps); n = 100 takes 2.6e23 steps |

**[syn] What is intrinsic.** Two factors are intrinsic: (# defect classes to calibrate = n²) × (# annealing stages
≈ n ln²n). Together they make n³ ln²n cavity-ratio estimates, each at constant precision. The n⁴ mixing bound is
where the slack is believed to be, and nobody has measured it at a size where it matters.

**Transfer of the lesson.**
- In JSV the expensive object is the *calibration of every cavity*. The chain itself is cheap per step.
- In our problem the cavity analogue, the n² pairwise D21 entries per layer, is cheap to store (FC carries it). What is
  expensive is *evaluating* it, because each entry is a sum over T = age·n atoms.
- So our bottleneck is the analogue of JSV's mixing factor (atoms per evaluation), not of its calibration factor.

### I2. Barvinok interpolation, Eldar–Mehraban, and cluster expansions

**Barvinok [thm].**
- Zero-free region: per Z ≠ 0 when |z_ij − 1| ≤ 0.5. The bound cannot exceed √2/2 (1601.07518 Thm 1.3; 1405.1303 had
  0.195). The hafnian has the same 0.5 (1601.07518 Thm 2.3).
- Truncation error: if g(z) = per(J + z(A − J)) has no zeros in |z| ≤ β, the degree-m Taylor polynomial of ln g errs by
  ≤ n/((m+1)β^m(β−1)). So m = O(ln n − ln ε).
- Coefficient cost: g^{(k)}(0) = (n−k)! Σ over k-tuples of ∏(a − 1), which is n^{O(k)}. The log coefficients follow by
  Newton's identities.
- *Pieces* = entries' deviations from the reference J. *Down step* = Taylor coefficient k, k deviations at a time.
  *Angle* = zero-free radius β.

**Patel–Regts [thm, 1607.01167 Thm 3.1].**
- The k-th inverse power sum is Σ_H a_{H,k}·ind(H, G).
- Multiplicativity forces a_{H,k} = 0 unless H is **connected**.
- There are ≤ (eΔ)^k connected sets at a vertex.
- So m = C ln(n/ε) coefficients cost (n/ε)^{O(1)} at bounded degree. *Local-to-global:* only connected local
  structures enter the log.

**Cluster expansions on expanders [thm, JKP 1807.04804; HPR 1806.11548].**
- Polymers are small connected deviations from a ground state.
- Expansion gives the polymer weight ≤ e^{−βα|γ|}(q−1)^{|γ|} (Potts), or ≤ (1+λ)^{−α|γ|} (hard-core).
- Kotecký–Preiss then gives convergence and zero-freeness. Runtime is (n/ε)^{O(log Δ)}.

**Eldar–Mehraban [thm, 1711.09457 Thm 18].**
- Setting: random A with i.i.d. entries of mean μ and variance 1, with μ ≥ (ln ln ln n / ln ln n)^{1/7}, roughly
  1/polyloglog n.
- Result: quasi-polynomial (1 ± 1/poly) approximation for 1 − o(1) of A.
- *Local-to-global, on average:*
  - E|g_A(re^{iθ})|²/(n!)² ≤ e^{r²}. This is a second moment, i.e. local.
  - Jensen's formula then gives E N_r ≤ 4r² roots in the disc.
  - A tube around one of 32/ε⁵ disjoint paths is root-free for most A.
  - Analytic continuation along that tube uses polylog-many derivatives.
- *Mean zero is the boson-sampling regime, believed hard.* The algorithm needs the mean μ to dominate.

**[syn] For us.**
- Our weights have mean zero: the hard regime.
- The only "mean" is the Perron outlier of the gated propagators (F2.2). Interpolating from that rank-one reference
  is Barvinok's path with J → Perron.
- Its degree is a rank. The measurement says it is ≈ n/4 (§0 item 5).

### I3. The Bethe permanent, belief propagation and covers

**Vontobel (1107.4196, Thm 39) [thm].**
- per_B(A) = limsup_M ⟨per(A^{↑M})⟩^{1/M}, averaging over all M-covers.
- The Bethe free energy is convex and belief propagation converges to it.

**Gurvits / Schrijver [thm].** per(A) ≥ per_B(A). Csikvári gave a proof via 2-lifts, which is literally local to
global: lifting cannot increase the permanent, and the limit of lifts is the tree.

**Anari–Rezaei (1811.02933 Thm 4) [thm].** per(A) ≤ √2^n per_B(A), tight for I_{n/2} ⊗ J_2. Also
per(J)/per_B(J) = √(2πn/e)(1 + o(1)).

**All of this needs A ≥ 0.** For signed or complex A, no Bethe guarantee is known. Anari–Rezaei note that even the sign
of per is hard.

**For us [meas + syn].**
- The interpolation stream already ran the cover experiment (`fresh-slate/breakthrough/interpolation/REPORT.md` item
  3). Covers that lift the input converge to the tree and score 2.05e-4, **48× worse** than the Gaussian closure.
  Vontobel's theorem identifies that tree limit as the Bethe value.
- Our §3 measurement locates what the Bethe value misses: the off-diagonal path pairs. These are signed, so they lie
  outside every Bethe theorem.

### I4. Godsil–Gutman estimators, Clifford signs, and sign migration

**Godsil–Gutman [thm].**
- Estimator: per A = E det(A ∘ ε)². The critical ratio is 3^{n/2} for ±1 signs.
- Karmarkar et al. (KKLLL): 2^{n/2} with complex unit or cube-root signs.
- Chien–Rasmussen–Sinclair (CRS): (1 + O(2^{−k/2}))^{n/2} with Clifford Cl_k. The quaternion case is computable in
  polynomial time with ratio (3/2)^{n/2}.
- Moore–Russell (0906.1702 Thms 1–2): with d×d matrices, the unsymmetrised estimator has ratio (1 + O(1/d))ⁿ. But the
  computable symmetrised determinant has ratio Ω(2ⁿ/n^d) at constant d, so d ≈ n² is needed.
- Chien–Harsha–Sinclair–Srinivasan (CHSS, 1101.1169):
  - the determinant over M₂(F) is as hard as the permanent;
  - over a fixed finite-dimensional algebra over F_p (p odd), it is easy if A/rad A is commutative and Mod_pP-hard
    otherwise.

**Mechanism [syn].**
- det² = Σ over pairs (σ, τ) of permutations. The signs kill the off-diagonal pairs σ ≠ τ in expectation and leave
  per = the diagonal.
- The variance is the off-diagonal mass. It factorises over the n rows, so a per-row reduction by a factor of order
  1/d compounds to an exponential gain.
- *Sign migration:* the algebra that cancels the cross terms makes the determinant itself hard.

**For us: the reversal [syn].**
- In Godsil–Gutman the *diagonal* pairing is the target, and the full sum (det²) is the cheap object.
- In our memory the full sum (the exact squared leg) is the target. The *diagonal* pairing (Bethe) is the cheap
  approximation, and the off-diagonal pairs are the content.
- So Godsil–Gutman signs would kill exactly what we need. Used instead as a sketch of the atom sum, they face an
  additive structure (§2 Prop. 1).

### I5. Stable and Lorentzian polynomials (drop i = ∂ᵢ)

**[thm]**
- Brändén–Huh (1902.03719): a polynomial is Lorentzian iff it has nonnegative coefficients and M-convex support, and
  every ∂^α f of degree 2 has ≤ 1 positive eigenvalue. Lorentzian implies 2(1−1/d)-Rayleigh, a weak negative
  dependence.
- Anari–Liu–Oveis Gharan–Vinzant (ALOV, 1811.01816 Thm 2.16): local log-concavity of all quadratic derivatives
  (links), plus indecomposability, implies global log-concavity. This is Oppenheim's trickle-down. The down–up walk
  then mixes in r log(n^r/ε).
- No result found treats Gaussian orthant or gate patterns as stable or Lorentzian.
- Pitt's theorem gives *positive* association for nonnegatively correlated Gaussians, which is the opposite direction.

---

## 2. The noncommutative generalisation

**Known.**
1. *NC randomisation of the permanent* (CRS; Moore–Russell; CHSS, above). Algebra-valued signs reduce variance
   multiplicatively over rows. A dichotomy (commutative semisimple quotient: easy; otherwise hard) says the
   determinant pays the bill.
2. *Operator scaling as noncommutative Sinkhorn / Bethe.*
   - Gurvits's capacity cap(T) = inf det T(X)/det X bounds the mixed discriminant: (n!/nⁿ)·cap ≤ D ≤ cap.
   - Alternating normalisation computes it (Gurvits 2004; Garg–Gurvits–Oliveira–Wigderson, GGOW). This is the
     two-conditional-expectation down–up process of the brief, with an NC van der Waerden bound as its "Bethe"
     certificate.
   - Like the commutative Bethe bounds it lives in the **positive** cone (completely positive maps).
3. *Wick versus free pairings.* For classical Gaussians, moments are hafnians (all pairings). For free semicirculars,
   only the noncrossing pairings survive. Evaluating a network on free inputs would compute the noncrossing part of
   the Wick expansion. That is a biased object, not an estimator.

**Proposition 1 (sign migration for additive sums) [prop].**
- *Setting.* Let D = Σ_{r=1}^T w_r (y_r ∘ y_r) z_rᵀ, with T atoms. Let ξ_r be i.i.d. random elements of a compact
  group, represented by d×d unitaries, with E ξ = 0 and E[ξ ⊗ ξ ⊗ ξ̄²]-type moments chosen so that the estimator
  X = (1/d) Re tr[(Σ ξ_r y_r) ∘ (Σ ξ_r y_r) (Σ w_r ξ_r^{−2} z_r)ᵀ] is unbiased.
- *Variance.* Var X ≥ c_G · Σ_{(r,r′,r″) not all equal} |w_r″ y_r y_r′ z_r″|² / d².
- *Cost.* X costs d² n-vector transports per layer.
- *Conclusion.* The variance per unit cost equals that of d² independent scalar sketches, up to the constant c_G.
- *Why the permanent differs.* For det-type estimators the cross-term mass is a *product* of n per-row factors
  (1 + O(1/d)). Each factor is reduced separately, so the gain is exponential in n at polynomial cost.
- *Proof sketch.* Expand X. The diagonal r = r′ = r″ gives tr(I)/d = 1. Each cross term carries
  (1/d) tr of a word in independent unitaries that is not identically I. Its second moment is Θ(1/d²) by Haar or
  2-design orthogonality (E|tr U|² = 1). Different index triples are uncorrelated. ∎
- *Corollary.* No algebra-valued random sketch makes the memory cheaper than its scalar counterpart. The cheap scalar
  version is measured to need ≥ 3.6e8 samples (§4).

**The gauge-group form (the cleanest NC statement this team found) [syn, with measured angles].**

*The groups and their averages.*
- The weight law is invariant under a gauge group at every hidden layer t: W_{t+1} ↦ g·W_{t+1} with g in
  - the diagonal sign group (Z₂)ⁿ,
  - the hyperoctahedral group B_n = (Z₂)ⁿ ⋊ S_n,
  - or the full O(n),

  each acting on the row (hidden-neuron) index.
- The network *function* is not invariant, so the scored means are not.
- Averaging a quenched functional F(W) over the gauge group of layer t is a conditional expectation E_t. It is the
  Haar average over the group, i.e. the projection onto the gauge-invariant part.
- Different layers' gauge groups commute and are independent, so {E_S : S ⊆ layers} is a **lattice of commuting
  conditional expectations**, a commuting-square lattice in Popa's sense.

*Forget one piece, re-randomise.* Re-randomising the sign of one weight row is the down–up step of the brief, in its
purest form.

*How the existing objects sit in this lattice.*
- **Godsil–Gutman signs are exactly the (Z₂)ⁿ gauge average.** Averaging over sign flips of the intermediate index
  kills the off-diagonal path pairs. So the **Bethe value = E over the (Z₂)ⁿ gauges of all layers**.
- **The region stream's fresh-weight lemma** ("Frobenius norms are the currency") is the O(n) average of the next
  layer.
- **The trace channel / mean-field** is the O(n)-invariant (trivial-isotypic) part of the squared leg in Sym²(ℝⁿ).

*Where the memory sits.*
- The memory is the **non-invariant isotypic part**: degree-2 characters ε_m ε_m′ for the sign group, the traceless
  symmetric component of Sym² for O(n).
- Measured angles at n = 1024 (relative Frobenius norm of D_old outside the invariant part):
  - **0.47–0.48** for all-layer (Z₂)ⁿ;
  - **0.40–0.45** for the O(n) trace part;
  - **0.46–0.48** for (Z₂)ⁿ applied at the *last* layer only.
- So averaging only the last fresh layer loses as much as averaging all of them (§4 cut curve).

*The "next theorem" (Conjecture 1 below):*
- the invariant/non-invariant angle of quenched third-order transport stays open as n → ∞;
- the non-invariant part is regenerated at every fresh layer, as the degree-2 component of the most recent gauge
  group, with coefficients equal to the full old leg.

*Link to team B.* This is the abelian or orthogonal case of Fourier analysis on the gauge group, which is
Jeronimo–Mittal–Roy's tool. The memory is exactly what lies off the trivial representation.

**Conjecture 1 (the Bethe defect of quenched third-order transport is O(1) in width) [conj].**
- *Setting.* Let J_a = G_1⋯G_a with G_i = diag(P_i)W_i, gated He-Gaussian layers with fresh weights. Let B_a be the
  diagonal-pairing value obtained by replacing (J_a)_{ra}² with ((G_1∘G_1)⋯(G_a∘G_a))_{ra}.
- *Claim.* For fixed age a ≥ 2, the relative Frobenius defect ε_a(n) = ‖D − D_B‖_F/‖D‖_F of the (2,1)-slice
  contraction converges to a positive constant ε_a(∞) as n → ∞, and ε_a(∞) increases to ε_∞ ≈ 0.45.
- *Equivalently.* The angle between the quenched leg algebra and its diagonal (commutative, annealed) subalgebra stays
  open: cos θ_a → c_a < 1.
- *Evidence.* 0.70 / 0.48 / 0.47 at n = 1024 for the all-old aggregate at t = 4 / 10 / 14. The width test is §4 T1.
- *Significance.* With fresh weights the Bethe sector is a commuting square: E_diag ∘ E_next = E_diag. The defect is
  the part of the memory that is not in any commuting square. This is the obstruction team F's NC trickle-down would
  need to bound.

**[syn] The general NC statement our domain supports.** The noncommutative versions of the permanent's local-to-global
tools fall into two kinds:
- those that buy *positivity or convexity* (operator scaling, capacity, Bethe for completely positive maps). These are
  efficient.
- those that buy *cancellation* (Clifford signs). These migrate the hardness.

The memory needs cancellation of signed, additive cross terms, so it sits on the wrong side of this line.

---

## 3. Transfers to our problem

Each transfer below states which object plays each role. None of them is the per-neuron or per-layer drawing: the
pieces are **pairs of paths** (Wick pairings of legs), **atoms of the old third cumulant**, or **gate patterns as
points of a polynomial's support**.

### T1. Bethe/diagonal pairing of the squared legs (strongest; measured)

**Roles.**
- *Pieces:* path pairs (p, p′) from birth neuron r (layer s) to target a (layer t). Y_ra² = Σ_{p,p′} w(p)w(p′).
- *Down step:* forget which intermediate neuron each path used, i.e. condition on the diagonal p = p′. That is the
  conditional expectation onto the commutative (annealed) algebra.
- *Up step:* re-pair.
- *Angle:* the cosine between the true old D21 and its Bethe value.
- *Global quantity:* old D21(t) to ≤ 5–10 %.

**Why it is not the per-neuron drawing.** The object is a pairing of paths, a hafnian-type structure. The same neuron
appears in both diagonal and off-diagonal pairs.

**Prediction and measurement [meas, MLP 0].**
- cos = 0.71 / 0.88 / 0.89 and relative error 0.70 / 0.48 / 0.47 at t = 4 / 10 / 14.
- The rank-one mean-field gives 0.70 / 0.44 / 0.40.
- The Bethe sector is ≈ 80 % of the memory's energy at depth.

**Cost.**
- Bethe value: one Hadamard chain per source, n³ per source-layer.
- Mean-field: O(n²) per source-layer, i.e. the trace channel. The interpolation stream closed that recursion with
  R² = 0.99.

**Kill criterion.** Defect > 10 %. **Killed as a stand-alone carrier.**

**What it explains.**
- F8.5: the ensemble mean of the off-diagonal pairs is zero, and they are the quenched part.
- F8.3: they are orthogonal to the young content.
- The interpolation stream's covers (Bethe = cover limit, 48× worse) and its trace channel (= the Bethe sector).

**Width test (pending).** Conjecture 1: measure ε at n = 256 and 512.

### T2. Barvinok around the Bethe reference: low-rank loop correction (measured; killed)

**Roles.**
- *Solvable reference:* the Bethe/mean-field part at full rank, O(n²) per source.
- *Interpolation:* D(z) = D_B + z·(off-diagonal pairs).
- *Degree-1 coefficient:* the fluctuation term, with legs truncated to rank k.
- *Zero-free analogue:* the off-diagonal part is *linear* in z for the (2,1) slice. The question is not degree but the
  **rank** needed for the loop correction once the coherent (nonnegative) sector is removed.

**Prediction.** If the Bethe reference absorbs the Perron/coherent energy, the fluctuation might compress below the
n/4 that F9.5 measured for the full legs.

**Cost model.**
- Bethe/mean-field: O(L² n²) in total, negligible.
- Rank-k fluctuation legs: ≈ 2k/n units per source-layer pair (k×n legs transported by n×n gates), so
  ≈ 120·2k/n·(4 products) ≈ 1 u per 1 % of n.
- k = n/8 costs ≈ 120 u. That is borderline. k ≤ n/16 is needed to fit next to the young tier.

**Test.** `t_atoms.py` (mf_plus_rank), MLPs 0 and 1, t = 10, 14, k = 64, 128, 256.

**Kill criterion.** Relative error > 10 % at k = 128.

**Result [meas].**

| | MLP 0, t = 10 | MLP 0, t = 14 | MLP 1, t = 10 | MLP 1, t = 14 |
|---|---|---|---|---|
| Bethe + rank-128 fluctuation | 0.240 | 0.186 | 0.239 | 0.204 |
| rank 128 alone | 0.283 | 0.216 | 0.294 | 0.242 |
| Bethe + rank-256 fluctuation | 0.106 | 0.079 | 0.101 | 0.085 |
| rank 256 alone | 0.114 | 0.086 | 0.115 | 0.094 |

**Killed.** The off-diagonal sector is not more compressible than the whole.

**What it would overturn.** F9.4 ("better bases save constant factors"), if a nonlinear (Hadamard-aware) basis beats
the linear one.

### T3. Godsil–Gutman / Clifford sketch of the memory (measured; dead with mechanism)

**Roles.**
- *Pieces:* atoms r of the old cumulant, T = age·n.
- *Down step:* random signs ξ_r cancel cross terms in expectation.
- *Up step:* averaging over K sketches.
- *Angle:* the critical ratio E|X|²/|EX|².

**Why it is not the per-neuron drawing.** The sketch vectors u = Σ ξ_r y_r are superpositions of all atoms.

**Measured [meas, MLP 0].**
- Single-sample relative variance 8.9e5 / 3.8e6 / 1.1e7, i.e. ≈ 0.06–0.2 T².
- K for 5 %: 3.6e8 / 1.5e9 / 4.4e9.
- At ≈ K·n² FLOPs per layer, the cost is ≥ 10⁵ times the budget.

**NC version.** Proposition 1: no gain per cost.

**What it explains.**
- F6.10: sampling is a Godsil–Gutman estimator whose signs are the inputs x; 1.4e7 samples for 1 %.
- Why "learned" or random compressions of the memory fail.

**Kill criterion.** Met.

### T4. Lorentzian structure of the gate law (measured; fails)

**Roles.**
- *Pieces:* gates.
- *Down step:* ∂ᵢ, conditioning gate i open.
- *Local check:* Hessian signature of the degree-2 derivatives.
- *Global quantity:* log-concavity, giving down–up mixing and capacity-based counting.

**Measured [meas, MLP 0, closure moments].**
- The Hessian of f(t, t₀) at 1 has 107, 141, 123, 104, 89, 75, 61, 51 positive eigenvalues at layers
  1, 2, 4, 6, 8, 10, 13, 16 (Lorentzian requires 1).
- The off-diagonal gate-covariance edge over p² is 1.9 → 4.8.

**Prediction.** Rescaling all gate correlations by c < 0.2–0.5 would restore the property. That is not available.

**What it explains.** F4.2 (spectral independence 2–8). It closes the "local log-concavity" route for gate patterns.

### T5. JSV annealing with cavity ratios: multi-stage calibration over depth (proposal)

**Roles.**
- *Stages:* layers, or a schedule that turns on the off-diagonal pairing strength z from 0 (Bethe) to 1.
- *Cavity ratios:* per-target ratios D21_old(a,b)/D21_Bethe(a,b).
- *JSV lesson:* only factor-2 accuracy is needed per stage if a cheap refinement exists.

**The obstacle.** Here the refinement step is sampling, which F12.2 kills (realisable variance reduction 1.1–1.5×).

**Status.** No test planned unless T2 leaves a residual whose ratio to the Bethe value is smooth across targets.

---

## 4. Tests run (n = 1024, bench `w1024_d16`)

| test | network | result |
|---|---|---|
| Bethe (diagonal-pairing) value of old D21 (age ≥ 3) | MLP 0, t = 4/10/14 | relative error 0.70/0.48/0.47; cos 0.71/0.88/0.89; old share of FC D21 0.41/0.76/0.84 |
| rank-one mean-field of squared legs | MLP 0 | 0.70/0.44/0.40 |
| Perron/top-k leg truncation (Barvinok around a rank-k reference) | MLP 0 | k=1: 1.00/0.96/0.91; k=16: 0.95/0.77/0.69; k=128: 0.58/0.28/0.22; k=256: 0.30/0.11/0.09 |
| GG sketch, cube-root signs, single-sample relative variance | MLP 0 | 8.9e5/3.8e6/1.1e7 (T = 2048/8192/12288 atoms) → K(5 %) = 3.6e8–4.4e9 |
| gate polynomial Lorentzian check | MLP 0, layers 1–16 | 51–141 positive eigenvalues; ratio of edge to shift 1.9–4.8 |
| Bethe + rank-k loop correction | MLPs 0, 1, t = 10/14 | k=128: 0.24/0.19 (MLP 0), 0.24/0.20 (MLP 1); k=256: 0.11/0.08, 0.10/0.085 |
| replication of rows 1–4 | MLP 1, t = 10/14 | Bethe 0.48/0.48 (cos 0.88/0.88); mean-field 0.45/0.43; GG relvar 5.8e6/9.8e6 |

These test objects inside FC's first-order representation (`fc_hook.py` = `fresh-slate/breakthrough/region/fc.py` +
a hook). They are not end-to-end MSEs. The tolerance used is region §4 (D21 to 5–10 % at layers 6–12).

## 5. Honest assessment (v0)

- **Theorem:** I1–I5 as cited. Proposition 1 is a short calculation (sketched).
- **Measurement:** §4, on two networks (t = 10, 14) and one network (t = 4). All are first-order (FC) objects.
- **Synthesis:**
  - the three-structure classification;
  - "old content = off-diagonal path pairings";
  - the reversal relative to Godsil–Gutman;
  - the positivity/cancellation line for NC tools.
- **Speculation:** Conjecture 1. T5.
- **Risks:**
  - The Bethe value I compute (diagonal pairing through every intermediate layer) is one choice. A partial Bethe
    (exact for the last k layers) interpolates to the exact value, and its curve in k would be more informative.
  - Proposition 1 bounds sketches with i.i.d. group-valued signs. Structured, non-i.i.d. or adaptive sketches are not
    covered.

---

## 6. Final (v1): no-go for three carrier classes, and the escape hatches left

This section answers the coordinator's question (21:25 UTC): either a carrier for the signed remainder at O(n³) per
layer, or a no-go for a class broad enough to say what the leaders must be doing. **No carrier was found. Below is a
no-go for three classes that together cover every structure this domain offers. Each comes with an end-to-end or D21
number at n = 1024, and the escape hatches are named.**

### 6.1 New measurements (final)

**Width test of Conjecture 1** (`t_width.py`, relative Frobenius error of the old (age ≥ 3) D21 of FC; *k = 3, 5 cuts
apply only to sources older than k, so they are lower bounds*):

WIDTH_TABLE

Reading:
- **The O(n)-invariant (mean-field/trace) defect is flat in width.** At t = 10 it is 0.44–0.54 at n = 256,
  0.45–0.47 at n = 512 and 0.44 at n = 1024. At t = 14 it is 0.32–0.52, 0.39–0.42 and 0.40.
- The (Z₂)ⁿ diagonal-chain defect falls from ≈ 0.62 to 0.48 and converges to it from above. That part is only the
  self-averaging of the diagonal.
- **Conjecture 1 is supported in its refined form**: the angle between the memory and the gauge-invariant subalgebra
  stays open, at ≈ 0.40–0.45 of the norm, with no visible decay over 256 → 1024.

**The cut curve: averaging only the last fresh layer equals averaging all of them.**
- Cut at k = 1: 0.477 / 0.460 at n = 1024, against full Bethe 0.482 / 0.466. The same holds at 256 and 512.
- Cut at k = 2 gives the same.
- *Interpretation:* the signed sector is **regenerated at every fresh layer**. It is the degree-2 (ε_m ε_m′)
  component of the newest gauge group, with coefficients equal to the full old legs. A diagonal pairing at layer t−1
  makes all earlier pairings self-average.
- So "local defect at the last layer ≈ global defect": the defects do not accumulate. Local-to-global holds here, but
  in the wrong direction for us: it does not help, because the one local defect is already the whole defect.

**Proposition 1, numerical check** (`t_prop1.py`, synthetic nearly orthogonal atoms, n = 48, T = 192):
- One Haar-U(d) sample has relative variance **0.92 / 0.86 / 0.83×** that of d² scalar cube-root samples at equal
  transport cost (d = 2 / 3 / 4).
- It is unbiased (the mean of 400 samples matches var/400).
- The algebra's readout costs d³ n² against d² n², so there is no net gain.

**End-to-end** (`t_e2e.py`, bench w1024_d16, FC with slices and κ4-mf; old = age ≥ 3; raw = final MSE − truth noise):

E2E_TABLE

The Bethe (gauge-invariant) sector carries ≈ 80 % of the memory's energy. **It buys only ≈ 1.7× over dropping the
memory, and leaves the estimator ≈ 30× above FC.** The signed 20 % of energy is what the score is made of.

The interpolation stream's trace-channel closure TCX scores 1.05e-6, the same number as the Bethe-sector FC on MLP 0
(1.05e-6). That is consistent with "trace channel = Bethe sector".

### 6.2 The no-go (measured; each class with its number)

A carrier of old content at ≤ 7 u per layer must not belong to any of the following classes.

**(N-a) Gauge-invariant carriers.**
- *Definition:* any carrier whose old-content state at a cut is invariant under the sign gauge (Z₂)ⁿ (a fortiori
  under O(n) or B_n) of the next fresh layer. This covers Bethe, belief propagation and graph covers; trace channels
  and mean-field or annealed transport of quadratic objects; the dilation charge (costate C6/C7); and any function of
  row norms, Hadamard squares or gate probabilities.
- *Ceiling:* old D21 error ≥ 0.40–0.48 at n = 1024, flat in width. End-to-end ≥ 1.05e-6 raw on MLP 0 (FC + Bethe
  memory), against 3.2e-8 for FC.
- *Why:* the memory is the non-trivial isotypic component, and a gauge-invariant state cannot represent it.

**(N-b) Rank- or degree-limited carriers.**
- *Definition:* Barvinok or Taylor expansion around any rank-k reference (Perron, Oseledets or Bethe), CP merge to R
  atoms, and Tucker.
- *Ceiling:* rank n/4 gives 0.09–0.11; rank n/8 gives 0.19–0.28 even after the Bethe sector is removed (T2). Region
  CP R = n gives 2.0e-7 end-to-end.
- *Why:* the signed sector is as high-rank as the whole (T2).

**(N-c) Unbiased randomised carriers.**
- *Definition:* Godsil–Gutman, Clifford, quaternion and Haar sketches of the atom sum, and input sampling (the same
  estimator with x as the signs).
- *Ceiling:* 1.5e9–4e9 samples for 5 %.
- *Why:* additive structure, so the algebra gives no gain per cost (Proposition 1, proved and checked).

**[syn] Together.** The memory is (i) gauge-non-invariant, (ii) of rank ≳ n/4 within each age, and (iii) needed
deterministically. Any carrier at the leaders' cost must therefore hold, at each cut, a **sign-sensitive object of
dimension ≳ n²/4 that aggregates all ages**. That is an n×n-sized signed state, not L of them. And its update must be
a contraction with the actual fresh weights (to keep the ε_m ε_m′ components), at O(n³).

### 6.3 Escape hatches (what the leaders could be doing; ranked by this team's credence)

**1. A single aggregated signed n×n state per cut, approximate rather than exact.**
- C2/C3 forbid an *exact* o(L n²) state. The cut curve shows the signed sector is regenerated each layer from the
  full old legs, and that it decays: cuts at k = 3, 5 lose only 0.39 / 0.24, a lower bound.
- So an approximate one-state recursion may exist if the signed contribution through a cut is dominated by its
  *projection onto a fixed n²-dimensional quadratic statistic of the past*: for example the matrix
  Σ_r w_r y′_r y′_rᵀ ⊙ (something one-dimensional in the third leg).
- *Test:* regress the cut-1 off-diagonal term on such aggregates and measure R². A small, cheap next experiment.

**2. Error budgeting across layers (region §6(b)).**
- The D21 tolerance is 5–10 % only at layers 6–12 (K_off large).
- A carrier accurate at those layers and crude elsewhere, plus inter-layer cancellation, is not excluded by any of
  N-a to N-c.

**3. Non-first-order representations.**
- All of N-a to N-c are measured inside FC's first-order chaos frame.
- A leader working in a different frame (for example per-neuron marginal laws with a K = 3 chain and regenerated
  sources, as in the public chain) carries "old content" in a different basis, where the age split may not be the
  natural cut.
- *Risk:* our classes are basis-dependent in their proof but not in their measurements. The e2e numbers are
  frame-independent statements about which *information* a carrier holds.

**What is *not* an escape hatch:**
- more positivity: Lorentzian fails at every layer (T4);
- more symmetry: the dilation sector is worth ≈ 2× (costate);
- noncommutative signs (Proposition 1).

### 6.4 Final honest assessment

- **Theorem:** I1–I5 as cited, and Proposition 1 (short proof; numerically confirmed).
- **Measured at n = 1024 on bench networks:**
  - Bethe/mean-field split (2 networks);
  - width law (256/512/1024, 2 networks each);
  - cut curve;
  - T2 kill;
  - Lorentzian failure;
  - end-to-end Bethe-sector FC (E2E_NETS networks).
- **Synthesis:**
  - the gauge-group reading (Godsil–Gutman = (Z₂)ⁿ average; the fresh-weight lemma = O(n) average; memory = the
    non-invariant isotypic part);
  - "regenerated at the newest fresh layer";
  - the three-structure classification.
- **Conjecture:** Conjecture 1 (refined: the invariant-angle defect is O(1) in width). The evidence is flat over
  256–1024 for the O(n)-invariant part.
- **Risks:**
  - all D21 measurements live in FC's first-order frame;
  - the mean-field used is an oracle from the exact legs (the interpolation stream's recursion reproduces it to
    R² ≈ 0.99);
  - k = 3, 5 cut values exclude young-old sources;
  - Proposition 1 covers i.i.d. group-valued signs only.
