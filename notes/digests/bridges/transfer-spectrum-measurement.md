# Transfer-spectrum measurement: does transported old third-cumulant content concentrate on a few modes?

### A test of the expander / transfer-operator / Gibbs bridge on the Phase 2 problem (width 128, depth 16), with the width scaling to 1024

*Prong-2 measurement for the bridges stream, 2026-10-01. It tests one candidate bridge on the competition problem: [competition-plan.md](../../competition-plan.md) §3.1 and the `old-content/` stream of [streams/README.md](../../streams/README.md). Inputs from other prongs: the vocabulary of Note 1 ([conditional-arrow-algebra.md](../../conditional-arrow-algebra.md): Gibbs measures on histories, the transfer operator §5.3, KMS states); the four-ingredient table of [local-to-global-unlocks.md](../../local-to-global-unlocks.md) §4; the dictionary discipline of [mlp-bridge.md](../../mlp-bridge.md); the theorem lists of the neighbouring digests [expanders.md](expanders.md), [hdx-spectral-independence.md](hdx-spectral-independence.md) and [arxiv-2609.38007.md](arxiv-2609.38007.md), which are used by reference and not re-derived. Nothing below changes the objects of prong 1.*

**Labels.** **Measured**: a number from the runs in §8. **Fact**: a theorem from a source read in this session (quoted or transcribed). **Fact (memory)**: a standard result I did not open in this session. **Derivation**: proved here, with the numerical check named. **Interpretation**, **Conjecture**: my reading. Bridges in §5 carry **THEOREM / KNOWN-LINK / ANALOGY / SPECULATION**.

**Script:** [`streams/theory/transfer_spectrum.py`](../../streams/theory/transfer_spectrum.py) (numpy only). **Raw outputs:** [`transfer-spectrum-results/`](transfer-spectrum-results/).

---

## 0. Verdict in brief

1. **The prediction is half right.** Old content does concentrate on the top singular directions of its propagator, and far more than chance: a random k-subspace keeps nothing until k ≈ n, while the top-k subspace of the propagator keeps 95 % (age 4) to 99.99 % (ages ≥ 13) of the energy of the D21 contribution at k = n/4 (§3.2). But the concentration is **polynomial in age and linear in width**. The propagator has no spectral gap in the expander sense. Its participation ratio follows the free-probability law PR ≈ n / (1 + Σ(r_l − 1)) ≈ n/(2·age) to within 1–4 % for ages ≤ 5 (Measured, §3.1). Its consecutive Lyapunov gaps are 0.013–0.018 per layer, which is 3–5 × 1/(2n) and not O(1).
2. **Modes needed for a 2 % D21 error.** Two error measures are used. ε_own is the error relative to the content's own D21 contribution. ε_tot is the error relative to the whole D21(z_t), which is what the chain needs.
   - Width 128, ε_own, target-side propagator basis: k = 96 (age 2), 64 (age 4), 40 (age 8), 26–28 (ages 13–15). As a fraction of n that is k/n = 0.75 → 0.2.
   - Width 128, ε_tot, merged "old tier" (everything older than w layers, one tensor per target layer), targets t ≥ 10: k = 48–56 (w = 4), 28–48 (w = 6), 12–28 (w = 8).
   - The data-adapted HOSVD basis of the tensor itself, an oracle that a chain does not have, needs about half as many modes: 24–32, 12–24 and 4–16 for the same three thresholds.
   - The second MLP agrees: own 96 / 72 / 40 / 24–28; tiers 48–56, 28–40, 14–24 (§3.6).
3. **The mode count is a fixed fraction of n, not a fixed number.** The Gaussian-ensemble formula of §2.4 needs only the propagator. It reproduces the measured ε_own and the measured k_2 % exactly at 12 of 15 ages (MLP 0) and 11 of 15 (MLP 1), and is within one grid step at the rest (§3.2, §3.6). Run on fresh He-initialised MLPs at widths 128, 256, 512 and 1024, it gives the same k_2 %/n at every width to ±0.02 (§3.4). At n = 1024 a 2 % (own) carrier therefore needs about 0.19–0.5 n ≈ 200–510 modes for ages 15 down to 4.
   - **Cost.** A Tucker core of rank k costs k³/n² units per layer. For the w = 4 tier with the propagator basis (≈ 400–450 modes at n = 1024) that is 55–85 units per layer, against the ≈ 3–4 units per layer that [costmodel](../../streams/costmodel/REPORT.md) leaves for old and κ₄ content.
   - Only the oldest band is affordable: content older than about 10 layers needs about 0.1 n modes for 2 % of D21, i.e. about 1 unit per layer.
   - This agrees with 504aldo's closed arithmetic ("Tucker core … wins only for r ≲ 150").
4. **One genuine "Perron" mode exists.** A single outlier direction emerges with depth: the mean direction μ_z(t). Its overlap with the top left singular vector is 0.6–0.94 at ages 8–15 (Measured, width 128, both MLPs) and 0.85 at n = 1024. The ratio s₁/s₂ reaches 1.5–2.0, and the outlier separates from the bulk more clearly as n grows (§3.1, §3.4). It is width-independent and shared by all sources. It carries 0–20 % of the old content's energy up to age 9 and 4–67 % at ages 10–15, where the old content is only 1–6 % of D21 (§3.1, §3.6). This is the one place where the expander / Perron–Frobenius intuition is literally true, and it is a rank-one effect.
5. **The expander is the ensemble average, not the layer.** The ensemble-averaged covariance map E_W[Wᵀ S W] = (2/n) Tr(S)·I is a conditional expectation onto the scalars: a perfect expander, all non-trivial eigenvalues 0. The averaged third-order map is identically 0 (Derivation, §2.5). A single layer's map S ↦ B S Bᵀ is a completely positive map with **one** Kraus operator, the opposite end from a quantum expander. What does mix in one layer is the *orientation* of the source relative to the future propagator (§2.4). That is why a formula using only the propagator predicts the error.
6. **On the faces themselves, an expander-type certificate holds unpinned.** The gate law of each layer is a law on {0,1}ⁿ (the face law of dictionary v1). Its unpinned spectral-independence constant η₀ = λ_max(Ψ) is 1.8–4.9 at width 128 and 2–8 at width 1024 (Measured, §3.5). It is bounded and does not grow like n, which is the hypothesis of the Anari–Liu–Oveis Gharan local-to-global theorem for the single-gate (Glauber, down-up) walk on faces. Pinnings were not tested.

So the data support "old content concentrates on its propagator's top modes". They do not support "few modes" at the competition width. The cheap-carrier hope survives only for the oldest band, for the single mean mode, and for a structure-adapted basis that is better than the propagator's (the HOSVD gap). The transported covariance J C(a_s) Jᵀ, which the chain can compute, is such a basis: it comes within 2–8 modes of the oracle when the source keeps its slices (§3.6). But it too needs a fixed fraction of n.

---

## 1. The question, the objects, and the dictionary

### 1.1 Setting (as in the atlases)

- Network and input: random ReLU MLP, width n = 128, depth L = 16, no biases, weights W_l with iid N(0, 2/n) entries (Measured: n·Var = 1.95–2.03, kurtosis 2.92–3.07), input x ~ N(0, I_n).
- Conventions: x@W, i.e. z_l = a_{l−1} W_l with a_{−1} = x, and a_l = relu(z_l). In column form A_l = W_lᵀ.
- Gates: Φ_l = P(z_l > 0) (atlas `gate_p`), D_l = diag(Φ_l).
- Atlases: two independent sample atlases of the same MLP (`seed3`, `seed4`; N = 5·10⁵ each). Two MLPs were measured (`mlp_00000`, `mlp_00001`).

### 1.2 Transfer maps, propagators, error measures (definitions)

- **Layer maps.** On symmetric 3-tensors, T_l(X) = A_{l+1}^{⊗3}(Φ_l^{⊗3} ∘ X), the "old content" map for third cumulants with births ignored. On symmetric matrices, 𝒯_l(S) = B_l S B_lᵀ with B_l = A_{l+1} D_l, the covariance-like map; in x@W form this is S ↦ Wᵀ diag(Φ) S diag(Φ) W.
- **Propagator** from a_s to z_t (s < t), with t − s weight matrices and t − s − 1 gates:
  J_{s→t} = A_t D_{t−1} A_{t−1} ⋯ D_{s+1} A_{s+1}, so that T_{t−1}∘⋯∘T_{s}(X) = J_{s→t}^{⊗3} X when the first map is ungated. The gated version G = D_t J (a_s → linearised a_t) is what `old-content/propproj.py` calls M_{s→t}; its PR is printed for cross-checking.
- **Source.** K_s = AD(κ3(a_s)), the all-distinct part of the post-activation third cumulant, from central moments of `post_M3` exactly as in `experiments/oracle_k3.py` (`central3`, `all_distinct`).
- **Transported content.** X_{s,t} = J_{s→t}^{⊗3} K_s, with D21 contribution d_{s,t}[a,b] = X_{s,t}[a,a,b] at layer t. It is compared with the atlas's D21(z_t)[a,b] = κ3(z_a, z_a, z_b).
  - Because K_s contains everything present at layer s, X_{t−w,t} *is* the merged "old tier" at threshold w: all content older than w layers, linearly transported.
- **Projection.** P_k = U_k U_kᵀ, where U_k holds the top-k left singular vectors of J_{s→t}. Then X̂ = P_k^{⊗3} X. Since U_kU_kᵀ J = J V_kV_kᵀ, this is the same as projecting the source on the top-k right singular vectors at birth, so it costs a k³ core and n×k legs.
- **Error measures.**
  ε_own(k) = ‖D21(X̂) − d‖ / ‖d‖,  ε_tot(k) = ‖D21(X̂) − d‖ / ‖D21(z_t)‖,  share = ‖d‖ / ‖D21(z_t)‖,
  and k_2 % = the smallest k on the grid with ε ≤ 0.02.
- **Noise correction.** Every ε is noise-corrected with the independent atlas: ε² = ⟨e_A, e_B⟩ / ⟨d_A, d_B⟩. The projection is linear and the two atlases' Monte Carlo noises are independent, so the cross product removes the noise energy from both numerator and denominator.

### 1.3 Naive points of this dictionary, stated before the results (rule P2.1)

- **N1 — births and renormalisation ignored.** The task's transport is the Φ³ pass-through only. The oracle ladder (§3.1 of the plan) shows that old content also reaches the all-distinct κ3(a) through the [D21(z) ⊗ C] and [K22 ⊗ C] diagrams with coefficients ≈ 3, so the linear transport is not all of what old content does. Measured consequence: even at age 2 the linearly transported content leaves 84 % of D21(z_t) unexplained (§3.2).
- **N2 — no re-masking.** The source is AD-masked at birth and then transported linearly ("noAD"). The alternative bookkeeping that re-zeroes repeated indices at every layer is a Hadamard mask. It does not commute with projection, and `old-content/propproj.py --ad` already showed it destroys the concentration (keep plateaus at 0.86 for k = 64). Both decompositions are exact identities, with different births; the linear-transfer picture applies only to noAD.
- **N3 — Monte Carlo noise.** One atlas's D21 has 1–6 % noise. The transported contribution of old sources is noisier: 2 % at age 1, 34 % at age 15 for the tiny oldest contributions. All ε are noise-corrected, but the oldest ages carry large error bars.
- **N4 — sample size.** One width and two MLPs on the atlas side. The width scaling uses fresh random MLPs and the ensemble formula, not atlases.
- **N5 — mean-field gates.** The gates are expected gates (Φ ∈ [0, 1]), not per-sample masks. At depth 25–30 % of them are frozen at 0 or 1, i.e. dead or always-on neurons on this input law (Measured).

---

## 2. Theory that predicts the measurement

### 2.1 The Sym^k functor: higher-order transfer has no spectrum of its own (THEOREM, elementary)

**Statement.** Let B ∈ ℝ^{n×n} have singular values σ_1 ≥ … ≥ σ_n and eigenvalues λ_1, …, λ_n. The map B^{⊗k} restricted to Sym^k(ℝⁿ), with the Frobenius inner product, has
- singular values {σ_{i_1}⋯σ_{i_k} : i_1 ≤ … ≤ i_k}, and
- eigenvalues {λ_{i_1}⋯λ_{i_k}}.

Moreover (B₂B₁)^{⊗k} = B₂^{⊗k}B₁^{⊗k}.

*Proof.* With B = UΣVᵀ, U^{⊗k} and V^{⊗k} restrict to orthogonal maps of Sym^k. Σ^{⊗k} is diagonal on the orthonormal basis of normalised symmetrised monomials. For the eigenvalues, Schur-triangularise B = QRQ*; then R^{⊗k} is triangular on monomials in lexicographic order.

*Check.* `selftest`, k = 2, n = 9: singular values agree to 7·10⁻¹⁵ and |eigenvalues| to 2.5·10⁻¹⁴.

**Consequence.** For both the covariance-like map (k = 2) and the third-cumulant map (k = 3), the whole transfer spectrum is the Lyapunov / singular spectrum of the matrix product J, raised to the k-th tensor power. A tensor power creates no gap: at order k the ratio of the top two singular values is still σ₁/σ₂.

### 2.2 Free probability: the participation-ratio law (Fact + Derivation)

**Fact** (Pennington–Schoenholz–Ganguli, arXiv:1711.04735, eq. (11) and Supplement Result 1, read). The S-transform turns free multiplicative convolution into a product. For the input–output Jacobian of a deep network, S_{JJᵀ} = ∏_l S_{W_lW_lᵀ} S_{D_l²}.

**Fact** (same source, §2.4–2.5, read).
- Gaussian (Wishart) factors: S_{WWᵀ}(z) = σ_w⁻²(1+z)⁻¹. A 0/1 gate of density p: S_{D²}(z) = (z+1)/(z+p).
- Moments of JJᵀ for L ReLU layers with Gaussian weights: m₁ = (σ_w²p)^L and m₂ = (σ_w²p)^{2L}(L+p)/p. At criticality the variance is L/p.
- Linear Gaussian networks: λ_max = L^{−L}(L+1)^{L+1} (Fuss–Catalan) and variance L.

**Derivation** (the form used here). Write φ = Tr/n and r(x) = φ(x²)/φ(x)².
- For free positive a, b, the free moment formula φ(abab) = φ(a²)φ(b)² + φ(a)²φ(b²) − φ(a)²φ(b)² gives **r(a⊠b) − 1 = (r(a) − 1) + (r(b) − 1)**.
- Since PR(J) = (Σσ²)²/Σσ⁴ = n / r(JᵀJ),
  PR(J_{s→t}) ≈ n / (1 + Σ_{l=s+1}^{t}(r(A_lᵀA_l) − 1) + Σ_{l=s+1}^{t−1}(r(D_l²) − 1)).
- Reference values: r(AᵀA) = 2 for a square Ginibre matrix (Marchenko–Pastur, ratio 1), and r(D²) = 1/p for a 0/1 mask. For ReLU with Gaussian weights this gives r − 1 = L/p, which is PSG's variance.
- *Check.* `selftest`: a product of 6 Gaussian 400×400 matrices with 5 Bernoulli(½) masks gives PR 31.6 against the predicted 33.3. The Ginibre-only value would be 57.1.

**Measured inputs.** r(AᵀA) = 1.99–2.05 per layer. r(D_l²) = 1.0 (layer 0) rising to 1.5–2.35. The expected gates behave like 0/1 masks of density about ½, because of N5.

### 2.3 Products of many matrices: Lyapunov regime, free regime, and the polymer picture (Fact)

From Hanin–Nica, arXiv:1812.05994, read:
- **Thm 1.** Take X^{(i)} = (p n_{i−1})^{−1/2} D^{(i)} W^{(i)}, with Bernoulli(p) diagonal D^{(i)} and iid entries of W satisfying (i) mean 0, variance 1; (ii) symmetry; (iii) finite moments; (iv) no atoms. Let β = (3/p − 1)Σ_i 1/n_i + (μ₄ − 3)/(p n₁)·‖u‖₄⁴. Then (n₀/n_d)‖M^{(d)}u‖² ≈ exp(N(−β/2, β)), in Kolmogorov–Smirnov distance O((Σ n_i⁻²)^{1/5}) and in moments up to O(Σ n_i⁻²).
- **Prop. 2.** For p = ½, the singular values of M^{(d)} equal in law those of the ReLU input–output Jacobian.
- **Corollary 3.** For ReLU, β = 5Σ1/n_i + ….
- **Reading.** The controlling parameter is depth/width. Here β ≈ 5·15/128 ≈ 0.59 at width 128 and 0.07 at width 1024, so at width 1024 the products sit much closer to the free (Fuss–Catalan-type) regime.

**Fact (as cited by Hanin–Nica §1.2, not opened).**
- Furstenberg–Kesten: the top Lyapunov exponent λ_max = lim (1/d) log‖M^{(d)}‖ exists.
- Oseledets: the multiplicative ergodic theorem, giving deterministic exponents and filtrations.
- Isopi–Newman: when d → ∞ first and then n → ∞, the density of normalised Lyapunov exponents is the triangle law h(λ) = 2λ on (0, 1). Tucci obtains the same global law in the other order of limits, but the local statistics depend on the order (Akemann–Burda–Kieburg).

**Fact (Hanin–Nica §1.4, read).** ‖M^{(d)}u‖² is a line-to-line partition function of a directed polymer on the complete multipartite graph. The D's act as "{0, 1}-valued spins on the vertices … restricting the allowed paths", and mean-zero weights produce "significant cancellation", so that Z_d does not grow exponentially when n grows with d. This is literally a sum over histories through open neurons, i.e. through faces, with signed weights. It is used in §5 (B3).

### 2.4 Orientation scrambling and a propagator-only error formula (Derivation)

**Claim.** J_{s→t} = R·A_{s+1}, where A_{s+1} is Gaussian, independent of the source K_s, and right-orthogonally invariant in law. Up to the dependence of the downstream gates on A_{s+1} (see §2.6), the orientation of the source relative to J's right singular vectors is therefore Haar-random. For any projection built from J, the expected D21 error then depends on the source only through O(n)-invariant quadratic data.
- Sym³(ℝⁿ) splits under O(n) into the harmonic (traceless) part and a trace part ≅ ℝⁿ (the vector v_k = Σ_i K_iik).
- By Schur's lemma, the expected error depends on K only through ‖K_harm‖² and ‖K_tr‖².
- An AD source has K_iik = 0, so it is **purely harmonic**: for AD sources the error is source-independent.

**Formula.** Let K = Sym(Z) with Z iid N(0, 1). Then E⟨K, X⟩⟨K, Y⟩ = ⟨X, Sym Y⟩, and ‖Sym(p⊗p⊗q)‖² = ⅓‖p‖⁴‖q‖² + ⅔‖p‖²(p·q)². With G = JJᵀ and G_k = P_k G P_k = U_kΣ_k²U_kᵀ, the cross terms collapse because P_kJJᵀ = G_k, giving
  E‖D21 error‖² / E‖D21‖² = 1 − f(k)/f(n),
  f(k) = ⅓(Σ_a g_a²)(Σ_a g_a) + ⅔Σ_a g_a h_a,  g_a = Σ_{p≤k} σ_p²U_ap²,  h_a = Σ_{p≤k} σ_p⁴U_ap².
All k are computed at once in O(n²) after one SVD.

*Checks.*
- Monte Carlo, n = 32, 300 symmetric Gaussian tensors: predicted / Monte Carlo ε = 0.865/0.869, 0.697/0.694, 0.421/0.422, 0.205/0.206 and 0.074/0.074 at k = 2, 4, 8, 12, 16.
- For AD-masked tensors the errors are 1–3 % larger at n = 32, a 1/n effect of the dropped trace part.
- On the atlas the formula reproduces the measured k_2 % at 12 of 15 ages (§3.2).

### 2.5 The annealed layer map is a twirl (Derivation) and is not the quenched one

For W with iid N(0, 2/n) entries independent of S: E[WᵀSW]_{ab} = Σ_ij S_ij E[W_ia W_jb] = (2/n) Tr(S) δ_ab. All odd moments vanish, so E[W^{⊗3}] = 0. Therefore:
- The **ensemble-averaged** covariance transfer is S ↦ (2/n)Tr(S)·I: up to scale, the trace-preserving conditional expectation onto ℂ·I. All its non-trivial eigenvalues are 0, a perfect expander.
- The ensemble-averaged third-order transfer is **zero**. Old κ₃ content survives only because one particular W is fixed (quenched).
- The **quenched** map S ↦ BSBᵀ is completely positive with a single Kraus operator. **Fact** (Hastings 2007, as recorded in [expanders.md](expanders.md) §9.10): random-unitary channels with D Kraus operators have |λ₂| → 2√(D−1)/D. A single conjugation (D = 1) is the degenerate end: no averaging, so no trivial eigenvector is isolated.
- The Pauli twirl in Chen–Rouzé's proof (ρ − ρ_{−A} = 2^{−2|A|−1}Σ_S [S,[S,ρ]], see [arxiv-2609.38007.md](arxiv-2609.38007.md) §8.7) is the analogue of the *annealed* map: a group average that equals a conditional expectation.

### 2.6 Rank-one outliers (Fact) and the mean mode

**Fact** (Benaych-Georges–Nadakuditi, arXiv:0910.2120, Thms 2.6–2.8, read).
- Setting: X_n is non-negative definite with limiting spectral law μ_X on [a, b], and X̃ = X_n(I + P) with P = θuu*, where u is in generic position.
- Threshold: the top eigenvalue leaves b **iff** θ > 1/T_μ(b⁺), with T_μ(z) = ∫ t/(z−t) dμ(t). It then converges to T_μ⁻¹(1/θ), and |⟨ũ, u⟩|² → −1/(θ²ρT′(ρ) + θ).
- Below the threshold the overlap tends to 0. Square-root edge decay makes the threshold finite (Prop. 2.9).

**Interpretation (mine).** The gates are not free from the weights in one direction. Φ_{l+1} = Φ(A_{l+1}μ_{a,l}/σ), so D_{l+1}A_{l+1}μ_{a,l} is close to a rectification of the mean: the gates are open exactly where the mean pushes the pre-activation up. A generic direction v, independent of A_{l+1}, has mean-square gain E‖D A v‖²/‖v‖² = 2E[Φ²]. That is 0.5 at layer 0 and 0.58–0.95 at layers 1–15 (Measured on both MLPs), i.e. below 1 whenever gates are uncertain. The mean direction is transported with a larger gain: the measured ratio ‖Jμ_a‖²/(‖μ_a‖² · mean s²) is 1.2–3.3 (MLP 0) and 1.2–5.9 (MLP 1) at n = 128, and 4–8 at n ≥ 256 after 8–15 layers. The result is a coherent rank-one direction amplified relative to the bulk at every layer: a growing multiplicative spike, the BGN mechanism iterated. The measured overlap with μ_z(t) and the growth of the outlier with width (§3.1, §3.4) are what the theorem class predicts once the spike passes threshold. This identification is not proved here: the spike is generated by the coupling, not given as an independent P.

### 2.7 Width scaling (Derivation sketch + Conjecture)

At fixed depth all objects have width-independent spectral laws as n → ∞:
- the Wishart factors;
- the empirical law of Φ (mean-field);
- their free products (§2.2).

The error functional f of §2.4 is a normalised trace-type functional of these products. Hence **k_ε/n converges to a constant** for every ε, and the number of modes needed for a fixed D21 accuracy grows linearly with width. The mean spike of §2.6 is the only width-independent structure, and it is rank one. Measured support: §3.4. *Conjecture:* the same holds for structured (non-generic) sources, such as the Wick part of κ3(a_s) built from C(z_s), because their legs are themselves free products of the same factors. Only the HOSVD-type gain of §3.2, a constant factor, would carry over. Not tested at n > 128.

---

## 3. Measurements

### 3.1 (a) Spectra of the propagators J_{s→t} (MLP 0; MLP 1 in §3.6)

Mean over the 16 − age pairs of each age:

| age | PR measured | PR free (§2.2) | r90 | r99 | r999 | PR of D_tJ (propproj convention) | s₁/s₂ | s₁/s₈ | s₁/s₃₂ |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 63.8 | 63.8 | 65.7 | 98.9 | 114.7 | 43.4 | 1.03 | 1.13 | 1.54 |
| 2 | 32.2 | 32.5 | 37.1 | 62.8 | 79.4 | 25.8 | 1.07 | 1.34 | 2.61 |
| 3 | 21.2 | 21.6 | 26.2 | 47.8 | 63.5 | 17.8 | 1.10 | 1.54 | 4.24 |
| 4 | 15.6 | 16.2 | 20.4 | 38.7 | 53.2 | 13.1 | 1.17 | 1.78 | 6.71 |
| 6 | 9.5 | 10.7 | 13.9 | 27.9 | 40.1 | 8.3 | 1.32 | 2.37 | 15.7 |
| 8 | 6.2 | 8.0 | 10.6 | 21.9 | 32.0 | 5.6 | 1.46 | 3.24 | 36.0 |
| 10 | 4.4 | 6.4 | 8.0 | 18.2 | 26.8 | 4.2 | 1.60 | 4.08 | 78.8 |
| 12 | 3.8 | 5.4 | 7.0 | 15.5 | 23.5 | 3.5 | 1.64 | 4.83 | 143 |
| 15 | 3.4 | 4.4 | 6.0 | 13.0 | 20.0 | 3.2 | 1.53 | 5.69 | 339 |

- **Free-probability law.** Free probability holds to 1–4 % up to age 5. Beyond that the measured PR falls below the free value, by 25 % at age 15; the outlier below accounts for this.
- **Cross-check with the old-content stream** (rule C3). Its PR of D_tJ by age on its atlas A1 (a different MLP) was 43.3, 25.4, 17.8, 13.3, 10.5, 8.8, 7.6, 6.5, …. Ours on MLP 0 is 43.4, 25.8, 17.8, 13.1, 10.3, 8.3, 6.8, 5.6: the same law.
- **Lyapunov spectrum.** Finite-time exponents of J_{0→15} (log sᵢ/15): top eight +0.059, +0.030, +0.007, −0.012, −0.026, −0.037, −0.043, −0.057; sᵢ for i = 16, 32, 64: −0.149, −0.330, −0.966. The mean consecutive gap is 0.0134 for i ≤ 16 and 0.0180 for i = 17–64, against 1/(2n) = 0.0039. There is **no O(1) gap anywhere** except at the top in the mean direction.
- **Mean-direction outlier** (Measured; overlaps are squared cosines):

| age | s₁²/mean s² | s₁/s₂ | ⟨u₁(J), μ_z(t)⟩² | ‖Jμ_a(s)‖² / (‖μ_a‖² · mean s²) | energy of the old content in the top mode, 1 − ε_own(1)² |
|---|---|---|---|---|---|
| 2 | 8.5 | 1.07 | 0.09 | 1.17 | 0.00 |
| 4 | 17.9 | 1.17 | 0.42 | 1.70 | 0.01 |
| 6 | 29.4 | 1.32 | 0.58 | 2.35 | 0.01 |
| 8 | 42.3 | 1.46 | 0.63 | 2.83 | 0.09 |
| 10 | 52.8 | 1.60 | 0.72 | 3.27 | 0.08 |
| 12 | 58.3 | 1.64 | 0.87 | 3.01 | 0.11 |
| 14 | 62.8 | 1.68 | 0.92 | 1.52 | 0.52 |
| 15 | 61.9 | 1.53 | 0.93 | 1.18 | 0.12 |

- **Alignment of top-k subspaces** (‖Q_kᵀQ′_k‖²_F/k; random = k/n = 0.06 / 0.12 / 0.25 for k = 8 / 16 / 32):
  - *Target side* (the U_k of all sources s for a fixed target t = 15, against s = 0): 0.34–0.91 at k = 8, 0.52–0.94 at k = 16 and 0.71–1.00 at k = 32. Old sources share a basis at the target, as the backward Oseledets picture predicts.
  - *Source side* (V_k of J_{s→t} against J_{s→15}, for s = 1 and 5): rises from 0.2–0.5 at t = s + 1 to 1. The forward filtration stabilises only over 5–10 layers, so a basis fixed at birth costs 7–25 % more modes (§3.2).

### 3.2 (b) The transported all-distinct third cumulant (MLP 0, noise-corrected)

Medians by age. The columns are ε_own and ε_tot at k = 4, 8, 16, 32, 64, then k_2 % as median (max) over sources.

| age | share | unexpl. | ε_own k=4 / 8 / 16 / 32 / 64 | ε_tot k=4 / 8 / 16 / 32 / 64 | k_2 % own | k_2 % tot |
|---|---|---|---|---|---|---|
| 1 | 0.563 | 0.831 | .999 / .996 / .979 / .868 / .527 | .562 / .565 / .565 / .505 / .310 | 128 (128) | 128 (128) |
| 2 | 0.479 | 0.840 | .990 / .958 / .853 / .531 / .135 | .471 / .438 / .388 / .245 / .060 | 96 (112) | 80 (96) |
| 3 | 0.421 | 0.905 | .959 / .870 / .697 / .342 / .032 | .410 / .294 / .238 / .130 / .013 | 80 (96) | 64 (80) |
| 4 | 0.357 | 0.946 | .919 / .804 / .563 / .226 / .010 | .267 / .234 / .147 / .060 / .004 | 64 (80) | 48 (56) |
| 5 | 0.282 | 0.978 | .875 / .758 / .468 / .146 / .003 | .199 / .151 / .104 / .033 / .001 | 56 (80) | 40 (48) |
| 6 | 0.230 | 0.989 | .895 / .744 / .391 / .093 / .001 | .118 / .090 / .066 / .016 / .000 | 48 (56) | 34 (48) |
| 8 | 0.119 | 0.960 | .760 / .608 / .270 / .046 / .000 | .067 / .053 / .030 / .004 / .000 | 40 (48) | 20 (28) |
| 10 | 0.050 | 0.991 | .617 / .432 / .187 / .021 / .000 | .033 / .026 / .010 / .001 / .000 | 36 (40) | 11 (24) |
| 12 | 0.026 | 0.990 | .669 / .433 / .152 / .011 / .000 | .019 / .011 / .004 / .000 / .000 | 30 (32) | 3.5 (10) |
| 15 | 0.003 | 0.999 | .645 / .446 / .130 / .005 / .000 | .002 / .001 / .000 / .000 / .000 | 28 (28) | 1 (1) |

"unexpl." is ‖D21(z_t) − d‖/‖D21(z_t)‖. At age 1 the remainder is the slice part of κ3(a_s), which the AD source omits by definition. At older ages it is births plus the non-pass-through diagrams (N1). The linearly transported content and the rest are nearly orthogonal.

**Bases compared** (k_2 % of ε_own, median over sources):

| age | U_k(J_{s→t}) | V_k(J_{s→L−1}) fixed at birth | eig C(z_t) | eig J C(a_s) Jᵀ | HOSVD of X (oracle) | random |
|---|---|---|---|---|---|---|
| 2 | 96 | 120 | 96 | 80 | 56 | ≈128 |
| 4 | 64 | 80 | 80 | 56 | 40 | ≈128 |
| 6 | 48 | 56 | 80 | 40 | 32 | ≈128 |
| 8 | 40 | 48 | 80 | 36 | 28 | ≈128 |
| 10 | 36 | 40 | 72 | 32 | 24 | ≈128 |
| 12 | 30 | 32 | 72 | 32 | 24 | ≈128 |
| 15 | 28 | 32 | 64 | 32 | 20 | ≈128 |

- Random: ε_own > 0.89 at k ≤ 64 at every age.
- The current covariance C(z_t) is a poor basis; it is dominated by young content.
- The transported covariance J C(a_s) Jᵀ is better than J at young ages: it closes a third to a half of the J-to-HOSVD gap at ages 2–8, because it encodes where the source's energy sits.
- The oracle HOSVD needs 1.25–1.7 × fewer modes. The actual sources are more concentrated than the generic ones of §2.4, because their legs come from the earlier propagation, but the gain is a constant factor.

**Ensemble formula (§2.4) against measurement.** ε_own at k = 4, 8, 16, 32, 64 (measured | ensemble), and k_2 % (measured / ensemble):

| age | measured | ensemble | k_2 % |
|---|---|---|---|
| 2 | .990 .958 .853 .531 .135 | .988 .956 .850 .575 .135 | 96 / 96 |
| 4 | .919 .804 .563 .226 .010 | .938 .832 .591 .236 .011 | 64 / 64 |
| 6 | .895 .744 .391 .093 .001 | .849 .678 .393 .096 .001 | 48 / 48 |
| 8 | .760 .608 .270 .046 .000 | .715 .522 .246 .040 .000 | 40 / 40 |
| 10 | .617 .432 .187 .021 .000 | .589 .398 .159 .017 .000 | 36 / 32 |
| 12 | .669 .433 .152 .011 .000 | .524 .326 .109 .009 .000 | 30 / 28 |
| 15 | .645 .446 .130 .005 .000 | .461 .262 .069 .003 .000 | 28 / 24 |

**Merged old tier.** For each target t and threshold w the tier is the single tensor X_{t−w,t}. The table gives k_2 % of ε_tot as: propagator basis U_k(J) / birth-fixed V_k / eig C(z_t) / HOSVD oracle.

| t | w = 3 | w = 4 | w = 5 | w = 6 | w = 7 | w = 8 |
|---|---|---|---|---|---|---|
| 8 | 56/80/64/40 | 48/56/48/32 | 32/40/32/20 | 20/32/20/12 | 10/16/12/4 | 1/4/4/4 |
| 10 | 80/96/64/40 | 48/64/56/32 | 40/48/40/20 | 28/40/32/16 | 20/32/24/8 | 12/16/12/4 |
| 12 | 64/80/80/40 | 48/64/64/32 | 48/56/48/24 | 40/48/40/20 | 32/40/32/16 | 24/32/24/12 |
| 14 | 64/80/80/40 | 56/64/80/32 | 48/56/64/32 | 48/48/48/24 | 32/40/40/20 | 28/32/32/16 |
| 15 | 56/56/80/40 | 56/56/64/32 | 48/48/64/32 | 40/40/56/24 | 40/40/48/20 | 28/32/40/16 |

### 3.3 (c) The symmetric-matrix transfer map S ↦ B_l S B_lᵀ

By §2.1 the singular values are s_is_j and the eigenvalues l_il_j (i ≤ j), where s and l are those of B_l (dim Sym² = 8256).

Per layer:
- operator norm s₁² = 1.9 (l = 0) rising to 3.5–5.4;
- spectral radius |l₁|² = 0.51–1.19;
- PR of the Sym² spectrum = 2024 (l = 0), then 750–1270;
- r90 = 970–2780 and r99 = 1760–5100.

The large gap between operator norm and spectral radius is the non-normality of B_l.

Products z_s → z_t, by age (mean):

| age | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 15 |
|---|---|---|---|---|---|---|---|---|---|
| PR(Sym²) | 1031 | 358 | 173 | 94 | 37 | 15 | 8.0 | 6.6 | 5.9 |
| r90 | 1481 | 656 | 365 | 227 | 106 | 55 | 33 | 26 | 18 |
| r99 | 2891 | 1557 | 959 | 639 | 332 | 194 | 128 | 100 | 71 |

The actual transported covariance c = J C(a_s) Jᵀ is 30–100 % of C(z_t), falling with age. Projected on the top-k left singular vectors of J in both indices, its relative error is:
- k = 16: 0.28 (age 4), 0.087 (age 8), 0.032 (age 12);
- k = 32: 0.12, 0.015, 0.002.

Covariance-like content concentrates much faster than third-order content. A Sym² object needs (Σ_{p≤k}σ_p²)² of the energy, not the cube.

### 3.4 (d) Width scaling (fresh He MLPs, depth 16, gates from 2·10⁴ Monte Carlo inputs, ensemble formula)

k_2 %/n for ε_own, mean over pairs:

| age | n = 128 | 256 | 512 | 1024 |
|---|---|---|---|---|
| 1 | 0.947 | 0.946 | 0.946 | 0.946 |
| 2 | 0.707 | 0.708 | 0.704 | 0.709 |
| 4 | 0.488 | 0.503 | 0.495 | 0.499 |
| 6 | 0.370 | 0.387 | 0.381 | 0.386 |
| 8 | 0.299 | 0.316 | 0.310 | 0.315 |
| 10 | 0.247 | 0.264 | 0.261 | 0.265 |
| 12 | 0.213 | 0.229 | 0.227 | 0.229 |
| 15 | 0.172 | 0.191 | 0.188 | 0.192 |

- **Absolute k at n = 1024** (ε_own = 20, 10, 5, 2, 1 %): age 8 needs 162, 214, 263, 322, 363 modes; age 15 needs 91, 124, 156, 197, 227.
- **PR·age/n** = 0.47–0.51 for ages ≤ 7 at all widths, i.e. PR ≈ n/(2·age). It falls to 0.39 at age 15 for n = 1024.
- **Mean-mode outlier at n = 1024**: s₁²/mean = 61 (age 10) and 128 (age 15); s₁/s₂ = 1.15 and 1.45; ⟨u₁, μ_z⟩² = 0.66 and 0.85. The outlier is sharper at larger n, as a BGN spike above threshold should be.
- The atlas MLP at n = 128 gives k_2 %/n of 0.31 (age 8) and 0.22 (age 15), against 0.30 and 0.17 for the fresh n = 128 MLPs.

### 3.5 (e) The face law as a Markov random field: unpinned spectral independence

**Definition.** The gate pattern g = 1[z_l > 0] ∈ {0,1}ⁿ has a law on {0,1}ⁿ: the face law of dictionary v1. Restricted to the uncertain gates (0.02 < p < 0.98), Ψ(i,j) = P(g_j | g_i) − P(g_j | ḡ_i) = Cov_ij/Var_i, and η₀ = λ_max(Ψ) = λ_max(Cor) − 1 ([hdx-spectral-independence.md](hdx-spectral-independence.md) §3.3: ALO Def. 1.1; Chen–Eldan Fact 23). The input is the atlas's `gate_GG` and `gate_p`.

**Width 128 (atlas, MLP 0).**

| layer | 0 | 2 | 4 | 6 | 8 | 10 | 12 | 15 |
|---|---|---|---|---|---|---|---|---|
| uncertain gates m | 128 | 126 | 105 | 78 | 76 | 54 | 69 | 53 |
| η₀ | 1.78 | 3.25 | 4.46 | 3.49 | 4.05 | 3.78 | 4.73 | 4.90 |

**Fresh MLPs, n = 1024.** m = 1024 → 472; η₀ = 2.1, 3.8, 5.1, 6.3, 6.9, 7.3, 7.0, 7.6 (same layers). η₀/m → 0.

**Gaussian surrogate check.** A zero-mean Gaussian surrogate with the same correlations (Sheppard's arcsine law, all p = ½) reproduces layer 0 exactly (1.78), which validates the estimator. At depth it gives η₀ = 15–29. The mean shift that freezes half the gates is what keeps the uncertain gates nearly independent.

### 3.6 Second MLP and the full-source contrast

**Second MLP** (`mlp_00001`, the same two sample seeds; `atlas_mlp1.txt`). Every qualitative statement above repeats.

(a) Propagator spectra:

| age | 1 | 2 | 4 | 6 | 8 | 10 | 12 | 15 |
|---|---|---|---|---|---|---|---|---|
| PR | 63.7 | 32.4 | 16.3 | 10.6 | 7.3 | 4.8 | 3.4 | 2.6 |
| free prediction | 63.7 | 32.6 | 16.4 | 10.9 | 8.2 | 6.6 | 5.5 | 4.4 |

Mean outlier: ⟨u₁, μ_z⟩² = 0.72 (age 8), 0.81 (age 10), 0.89 (age 12), 0.94 (ages 14–15); s₁/s₂ = 1.26 → 1.98; mean-direction gain ratio 3.5 → 5.9.

(b) Transported third cumulant:

| age | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 15 |
|---|---|---|---|---|---|---|---|---|
| k_2 % own, propagator basis | 96 | 72 | 52 | 40 | 32 | 28 | 26 | 28 |
| k_2 % own, HOSVD oracle | 80 | 52 | 40 | 32 | 24 | 20 | 20 | 24 |
| share | 0.395 | 0.235 | 0.140 | 0.093 | 0.064 | 0.036 | 0.015 | 0.003 |

- The ensemble formula gives exactly the measured k_2 % at 11 of 15 ages and is within one grid step at the other four.
- Merged old tier, targets t ≥ 10, k_2 % of ε_tot (propagator / HOSVD): w = 4: 48–56 / 32; w = 6: 28–40 / 20; w = 8: 14–24 / 8–16.
- Energy of the old content in the top (mean) mode: ≤ 4 % at ages ≤ 9; 14, 49, 60, 67 and 35 % at ages 10–14.

(e) Face law: η₀ = 1.69 (layer 0, equal to the Sheppard value), then 2.5–5.6. The Sheppard surrogate gives 2.8–28.6.

**Full source** (κ3(a_s) with slices, so its partial trace is non-zero; MLP 0; `atlas_mlp0_fullsource.txt`).

*Code validation.* At age 1 the transport identity D21(z_{s+1}) = D21(A_{s+1}^{⊗3}κ3(a_s)) must hold at the sample level. Measured: share 1.000, unexplained 0.000.

| age | 2 | 4 | 6 | 8 | 10 | 12 | 15 |
|---|---|---|---|---|---|---|---|
| share | 0.91 | 0.72 | 0.57 | 0.47 | 0.33 | 0.16 | 0.04 |
| unexplained | 0.20 | 0.37 | 0.48 | 0.57 | 0.70 | 0.85 | 0.97 |
| k_2 % own: propagator U_k(J) | 96 | 56 | 48 | 32 | 28 | 22 | 20 |
| k_2 % own: transported covariance eig J C(a_s) Jᵀ | 52 | 32 | 24 | 20 | 18 | 18 | 20 |
| k_2 % own: HOSVD oracle | 40 | 24 | 20 | 16 | 12 | 14 | 16 |
| k_2 % own: ensemble prediction | 96 | 64 | 48 | 40 | 32 | 28 | 24 |

Linear transport of the whole cumulant, slices included, explains far more of D21(z_t) than the AD part alone (unexplained 0.2–0.6 against 0.84–0.96 at ages 2–8). Projecting this source on the propagator basis needs 10–30 % fewer modes than the AD source, and fewer than the ensemble prediction at ages ≥ 7. That is the predicted signature of a non-zero trace part (§2.4), which concentrates on two legs.

The **transported covariance basis** comes within 2–8 modes of the HOSVD oracle at every age ≥ 4. For content of covariance-response type, a basis the chain can compute (one n×n sandwich) is nearly optimal.

Because the share is larger, the ε_tot requirement is stricter. Merged tiers at t ≥ 10 (propagator / HOSVD): w = 4: 56 / 20–24; w = 6: 32–48 / 12–16; w = 8: 16–40 / 8–12. At width 1024 these remain fixed fractions of n by §2.7.

---

## 4. Interpretation

**Does the data support "old content concentrates on few modes"?**
- **Concentration: yes.** Against a random basis the propagator basis is enormously better, and it improves with age: k_2 % (own) falls from 96 at age 2 to about 28 at age 15. Old sources share their target-side subspace.
- **"Few modes": no, at the competition width.**
  1. The concentration law is PR ≈ n/(2·age), the free product of Marchenko–Pastur factors and half-density gates. It has no gap. A 2 % D21 error needs a fixed *fraction* of n, about 0.2–0.5 n for ages 15–4, at every width tested (128–1024), and the ensemble formula that shows this reproduces the atlas measurement.
  2. The only width-independent structure is one outlier, the mean direction. It carries 0–20 % of the old content's energy up to age 9; at ages 10–14 in MLP 1 it carries more, where the old content is only 1–6 % of D21.
  3. The data-adapted oracle basis saves a constant factor of 1.25–1.7, not a power of n.

**How many modes for a 2 % D21 error?** The criterion that matters for the chain is ε_tot, i.e. 2 % of the whole D21(z_t), spent on one tier.

| tier (targets t ≥ 10) | n = 128, propagator basis | n = 128, HOSVD oracle | n = 1024, propagator basis (× 8, §3.4) | Tucker core k³/n² at n = 1024 (units per layer) |
|---|---|---|---|---|
| everything older than 4 layers | 48–56 | 24–32 | ≈ 380–450 | 54–84 |
| older than 6 | 28–48 | 12–24 | ≈ 220–380 | 11–54 |
| older than 8 | 12–28 | 4–16 | ≈ 100–220 | 0.8–11 |
| single source of age ≥ 10 | 1–11 (median), ≤ 24 (max) | – | ≈ 10–190 | ≤ 7 |

**Caveats.**
- The n = 1024 column assumes the old content's *share* of D21 is width-independent. The share is a mean-field quantity (gate statistics), but it was not measured at n = 1024. The oracle1024 stream can check it.
- Leg transport adds k/n units per layer and the birth projection is extra.
- The ≈ 3–4 units per layer available for old content is the costmodel stream's number ("about 45–60 u for old-source and κ₄ content").

**What the cheap-carrier hope can still use (Interpretation).**
1. **The oldest band** (ages ≥ 8–10): about 0.1 n shared modes, under 1 unit per layer.
2. **The mean mode**: one direction, free to track (it is μ_z). It carries up to 20 % of the content of ages ≤ 9 and up to 67 % at ages 10–14 (MLP 1).
3. **A basis better than the propagator's**, from structure the chain already has. The HOSVD gap says the sources are not generic: their legs come from C(z_s). J C(a_s) Jᵀ already closes a third to a half of that gap at ages 2–8. A sharper candidate is the leg structure of the Wick term, the factors J D_Φ C(z_s), which is O(n²) per source.
4. **Covariance-type content is cheap.** Sym² content needs the square, not the cube, of the spectral energy fraction (§3.3). For sources with slices, the transported-covariance basis is within 2–8 modes of the oracle (§3.6). The trace part of a third-order source is carried by a vector and a Gram matrix (§2.4), so old content in "covariance-response" form (504aldo's reading of the leaders) is the right place to look. That agrees with the plan's revision log.

**Rule P2.1 (charging the negative to the dictionary).** The "few modes" negative is most plausibly charged to:
- N1: the Φ³ pass-through is not the whole of old content's route to D21. The diagrams [D21(z) ⊗ C] and [K22 ⊗ C] route it through objects of covariance type.
- N2: the noAD choice; under AD the result is worse, not better.

It is *not* charged to the theory. The next dictionary version should carry old content through the identified first-order diagrams, in their covariance-response form, and redo §3.2 for that carrier.

---

## 5. The user's intuition, tested bridge by bridge

**Verdict.** Gibbs / Markov-random-field structure and expander theory meet in this problem in exactly two places, both measured here:
- the rank-one Perron mode of the propagator;
- the bounded spectral independence of the face law.

They do *not* meet in the bulk transport of old content, which is a free-probability object with no gap. The quenched layer is a one-Kraus channel; only its ensemble average is a perfect expander.

**B1 — THEOREM (§2.1).** The transfer spectra at orders 2 and 3 are the tensor powers of the matrix-product spectrum. No higher-order gap exists that is not already a first-order Lyapunov gap.

**B2 — KNOWN-LINK + Derivation (§2.2–2.3).** Products of random gated matrices obey free multiplicative convolution (Fuss–Catalan-type laws). The participation ratio is n/(1 + Σ(r_l − 1)). The Lyapunov regime (n fixed, d → ∞) and the free regime (d fixed, n → ∞) do not commute (Hanin–Nica §1.1–1.2). The competition shape (d/n = 1/64) is deep in the free regime: concentration ∝ 1/age, not e^{−gap·age}.

**B3 — KNOWN-LINK (Hanin–Nica §1.4) + ANALOGY to Note 1.**
- The old-content transport ‖Jv‖² is a line-to-line polymer partition function over histories through open gates (faces), with signed weights.
- Note 1 §5.2 identifies KMS states of the arrow algebra with positive Gibbs measures on histories. Their transfer operator (§5.3) is a positive operator, so it has a Perron–Frobenius eigenvector and, under primitivity, a gap (Birkhoff contraction of the Hilbert metric, tanh(Δ/4): Fact (memory)).
- Our transfer operator is the same sum over histories with *signed* weights. Sign cancellation is what removes the Perron gap in the bulk (Hanin–Nica: "mean zero … significant cancellation").
- In the mean direction the gates align the signs (rectification). The weights become effectively positive along it and a Perron mode reappears: the measured outlier, B6.
- So the Gibbs-measure picture of Note 1 is the *positive (mean) sector* of the old-content transport, and the bulk is its sign-problem sector. *Interpretation; the identification of the mean sector with a positive cone is not proved.*

**B4 — Derivation (§2.5) + KNOWN-LINK (Hastings, via expanders.md).** The annealed covariance map is a twirl (conditional expectation onto scalars, λ₂ = 0). The quenched map has one Kraus operator. Every expander-type statement about the layer maps is true of the ensemble average, and the average kills third-order content identically.

**B5 — Derivation (§2.4).** Orientation mixing is perfect in one layer: right-orthogonal invariance of A_{s+1}. This is the precise sense in which "the walk could be anywhere" (local-to-global §1, the expander mixing principle) holds here. It holds for the *orientation* of the source, not for its *amplitude*. Its use: any projection's error is predicted from the propagator alone (verified, §3.2), so carrier design reduces to spectral data of J.

**B6 — KNOWN-LINK (BGN Thms 2.6–2.8) + Measured + Interpretation.** A rank-one multiplicative spike above the free-probability threshold gives an isolated top singular value whose singular vector overlaps the spike direction. Measured:
- overlap 0.63–0.93 with μ_z(t) at ages 8–15 (n = 128) and 0.85 at n = 1024;
- s₁/s₂ ≈ 1.5;
- the gain of the mean direction is 1.2–3.3 × the bulk mean at n = 128 (MLP 0), up to 5.9 (MLP 1), and 4–8 at n ≥ 256.
This is the single "Perron–Frobenius / spectral gap" of the problem.

**B7 — KNOWN-LINK (ALO; CLV; hdx digest §3.3) + Measured (§3.5).**
- The face law of each layer is an MRF on {0,1}ⁿ with unpinned η₀ = 1.8–4.9 (n = 128) and 2–8 (n = 1024), bounded and not growing like n. That is the unpinned half of the hypothesis under which the single-gate down-up (Glauber) walk on faces mixes polynomially. ALO Thm 1.3 needs all pinnings; CLV's O(n log n) needs bounded degree, which this dense law does not have.
- **What must be true** for the full theorem: η_i bounded under every pinning of i gates. Pinning is restriction to a sub-frame (Note 2 §2), so this is a statement about restrictions, U6 of local-to-global §5.
- **Test:** pinned η from `gate_GG` conditional on one or two gates. It needs third-order gate statistics; the atlases do not store them.

**B8 — ANALOGY (strong) → Guard: the Chen–Rouzé "time-averaged, detailed-balanced Lindbladian with single-Pauli jumps on A".** [CR]'s construction ([arxiv-2609.38007.md](arxiv-2609.38007.md) §8) has four ingredients. Each maps to something here; two survive.

| CR ingredient | role in [CR] | analogue here | survives? |
|---|---|---|---|
| detailed balance (KMS-symmetric generator) | self-adjointness in ⟨·,·⟩_ρ, a stationary state | the layer map B^{⊗3} is not self-adjoint in any natural inner product. The Φ-weighted inner product would need B = D⁻¹BᵀD, false for random W. The quiver is acyclic (local-to-global §6, "Acyclicity") | **no** |
| stationarity, ℛ[ρ] = ρ | the recovery map fixes the Gibbs state | no stationary old content: the layers are i.i.d. random maps. The ensemble-stationary object is the twirl's fixed point (multiples of I) | **only annealed** |
| time averaging, ℛ_t = (1/t)∫₀ᵗ e^{sℒ}ds, ℰ(ℛ_t†O) ≤ 2/t | a gap-free small Dirichlet form | the old pool at layer t is Σ_s J_{s→t}^{⊗3}B_s, a *discounted* average over ages of the transfer cocycle applied to births. Its recursion P_{t+1} = T_{t+1}P_t + (birth) is a random affine recursion. Fact (memory): Kesten 1973 / Bougerol–Picard 1992, stationary iff the top Lyapunov exponent is negative. The discount factor is the generic gain² ≈ 2EΦ² < 1, so this is a resolvent, not a Cesàro mean, and the gap-free 2/t mechanism has no counterpart | **reshaped** |
| single-Pauli jumps on A generate the algebra on A | local generators; Leibniz walks words back to single sites at cost 2^w | single-gate flips on a block of neurons generate all gate patterns on the block: the Glauber walk on faces of B7. Its certificate (bounded η) is measured unpinned | **yes, on faces** |

The analogue that survives is a time-averaged Glauber dynamics on the face law of a block A of neurons. It is reversible with respect to the face law, which is classical detailed balance. Its Cesàro average would be an approximate recovery map for the gates on A from their boundary; in the dense mean-field case the "boundary" is all other uncertain gates. **SPECULATION:**
- *What it would buy:* a quantitative Markov property for faces. The conditional law of the gates on A given the rest is recovered by a local dynamics.
- *What it would need:* conditional mutual information of the face law across a "cut" to decay. In a dense, mean-field law there is no distance, so the CR / [Y] distance-decay statements have nothing to act on. The right object is the influence matrix (B7), not CMI against distance.
- *Prong-1 guard:* faces are the framework's objects, so this is stated at the level of the theory (a reversible walk on faces with a KMS-type state), not at the level of the old-content carrier.

**B9 — SPECULATION (NCG).** The algebraic shape of the result:
- the annealed map is a conditional expectation (a twirl);
- the quenched map is a single-Kraus conjugation;
- the difference carries all old κ₃ content.

In Connes's terms the old content is the part of the transfer cocycle that the conditional expectation does not see: a fluctuation of a random cocycle around its average. The Perron (mean) mode is the one direction where the cocycle has a fixed sign, i.e. where it is cohomologous to a positive one. Note 1 §6 already says that "is the cocycle a coboundary" is the sufficiency question. *What must be true for this to become a statement:* a definition of the signed-weight transfer cocycle on the history groupoid of a network, and a theorem relating its sign structure to the outlier spectrum. Nothing here proves that.

---

## 6. Known theorems connecting this area to the others

| link | theorem (hypotheses → conclusion) | source | status |
|---|---|---|---|
| transfer spectrum ↔ tensor powers | B^{⊗k}|_{Sym^k} has singular values ∏σ_{i_j}, eigenvalues ∏λ_{i_j} | §2.1, `selftest` | THEOREM (elementary) |
| products ↔ free probability | S_{AB} = S_A S_B for free A, B; S_{JJᵀ} = ∏S_{W_lW_lᵀ}S_{D_l²} | PSG eq. (11), Result 1; BGN §2.5.2 | Fact |
| products ↔ free probability | ReLU, Gaussian weights, criticality: Var(JJᵀ) = L/p; linear Gaussian: λ_max = L^{−L}(L+1)^{L+1} | PSG eq. (19), §2.4 | Fact |
| products ↔ free probability | r(a⊠b) − 1 = (r(a) − 1) + (r(b) − 1); PR = n/r | §2.2 | Derivation |
| products ↔ ergodic theory | Furstenberg–Kesten top exponent; Oseledets MET; Isopi–Newman triangle law; non-commuting limits | as cited in Hanin–Nica §1.2 | Fact (as cited) |
| products ↔ ReLU nets | log-normality of ‖Mu‖² with β = (3/p − 1)Σ1/n_i + …; p = ½ gives the ReLU Jacobian | Hanin–Nica Thm 1, Prop. 2, Cor. 3 | Fact |
| products ↔ Gibbs measures | ‖Mu‖² = line-to-line polymer partition function over open paths, with spins D | Hanin–Nica §1.4 | Fact (identity) + ANALOGY (B3) |
| Gibbs ↔ transfer operators | positive operators contract the Hilbert projective metric by tanh(Δ/4); Ruelle–Perron–Frobenius: for Hölder potentials on mixing subshifts of finite type the transfer operator has a simple top eigenvalue with a gap, the Gibbs state is its eigenmeasure, correlations decay exponentially | Birkhoff 1957; Ruelle, Bowen LNM 470 | Fact (memory) |
| spikes ↔ free probability | rank-one multiplicative spike: outlier iff θ > 1/T_μ(b⁺); eigenvector overlap formula | BGN Thms 2.6–2.8 | Fact |
| expanders ↔ quantum channels | random-unitary channels: |λ₂| → 2√(D−1)/D; quantum Alon–Boppana | Hastings, via expanders.md §9.10 | Fact (there) |
| expanders ↔ twirls ↔ NCG | E_W[WᵀSW] = (2/n)Tr(S)I, E_W[W^{⊗3}] = 0 | §2.5 | Derivation |
| Markov semigroups ↔ Gibbs (quantum) | time-averaged KMS Lindbladian with single-Pauli jumps on A is a quasi-local approximate recovery map; local Markov at all temperatures; gap ⇒ global Markov | Chen–Rouzé Thm III.1, Cors III.1–III.2, App. B (via arxiv-2609.38007.md §8) | Fact (there) |
| random affine recursions | X_n = A_nX_{n−1} + B_n with i.i.d. (A_n, B_n): unique stationary perpetuity iff top Lyapunov exponent < 0 (plus E log⁺‖B‖ < ∞ and non-degeneracy) | Kesten 1973; Bougerol–Picard 1992 | Fact (memory) |
| MRF ↔ HDX ↔ mixing | spectral independence of a law on {0,1}ⁿ and all pinnings ⇔ local spectral expansion of its complex ⇒ Glauber gap (product formula) | ALO Thms 1.3, 1.5, 3.1; CLV; via hdx-spectral-independence.md | Fact (there) |
| MRF ↔ expanders | Dobrushin ‖R‖ < 1 in any matrix norm ⇒ uniqueness and O(n log n) mixing | DGJ, via expanders.md §9.2 | Fact (there) |
| MRF ↔ Gibbs | Hammersley–Clifford: positive Markov random fields = Gibbs laws with clique potentials | via expanders.md §9 | Fact (memory there) |
| local-to-global ↔ this measurement | one layer scrambles orientation (Haar), so the error of any projection is a function of the propagator spectrum; old sources share the target-side top subspace (backward Oseledets) | §2.4, §3.1 | Derivation + Measured |

---

## 7. Messages to the other prongs and streams

**To `old-content/`.**
- The propagator basis needs k ≈ 0.2–0.5 n for 2 % (own), and the fraction is width-independent (§3.4). Rank-k carriers built on the linear Φ³ transport therefore do not scale to n = 1024 except for ages ≳ 8–10.
- The sources-with-slices convention (noAD births B_s, as in `propproj.py`) reports higher keep at equal k than the AD source here: 0.89 against 0.73 at age 8, k = 16, though the conventions also differ in where D21 is read. The harmonic/trace split of §2.4 predicts that the trace part concentrates on two legs. The full-source run (§3.6) confirms the direction: 10–30 % fewer modes than the AD source, and the transported-covariance basis comes within 2–8 modes of the HOSVD oracle.
- Worth testing next: a basis built from the Wick legs J D_Φ C(z_s), which the chain has. The target is to close the 1.25–1.7× HOSVD gap.

**To `oracle1024/`.** Two cheap measurements at n = 1024 would turn this into a closed number:
- the share of old content in D21(z_t) by age;
- the mean-mode energy fraction.

The width scaling of k/n is already settled by the ensemble formula.

**To `costmodel/`.** Use the rows of the §4 table as candidate designs. Only "older than 8 layers, about 0.1–0.2 n modes" plausibly fits the old-content budget.

**To prong 1.**
1. The transfer operator on histories with *signed* weights has a positive (mean) sector with a Perron mode and a sign-cancelling bulk (B3, B6). Is there a cocycle formulation in which the positive sector is the KMS / Gibbs part of Note 1 §5 and the bulk is a fluctuation around a conditional expectation (B9)?
2. Faces carry a reversible walk (single-gate flips) whose certificate (bounded spectral independence) is measured unpinned (B7). This is U6 of local-to-global §5 with a first number attached.

**To prong 3.** Add one guard to local-to-global §6: *a single quenched layer is a one-Kraus channel*. Expander statements about layer maps hold for the ensemble average (a twirl) and fail for the realisation. The realisation's concentration is the free-probability ∝ n/age law, with one BBP-type outlier.

---

## 8. Reproduction, results files, sources

```
source /root/whest/bin/activate; export OPENBLAS_NUM_THREADS=2
cd notes/streams/theory
python transfer_spectrum.py selftest
python transfer_spectrum.py atlas <atlas128>/seed3/mlp_00000.npz <atlas128>/seed4/mlp_00000.npz --save ../../digests/bridges/transfer-spectrum-results/atlas_mlp0.npz
python transfer_spectrum.py atlas <atlas128>/seed3/mlp_00001.npz <atlas128>/seed4/mlp_00001.npz --save .../atlas_mlp1.npz
python transfer_spectrum.py atlas <...mlp_00000 pair> --source full --save .../atlas_mlp0_fullsource.npz
python transfer_spectrum.py scaling --widths 128,256,512,1024 --depth 16 --mc 20000 --save .../scaling.npz
```

Run times on the shared 4-core box: atlas 95–305 s, scaling 125–195 s, selftest 10 s. Outputs: `transfer-spectrum-results/{selftest,atlas_mlp0,atlas_mlp1,atlas_mlp0_fullsource,scaling}.txt` and the `.npz` arrays.

### 8.1 Sources

**Read in this session:**
- J. Pennington, S. Schoenholz, S. Ganguli, *Resurrecting the sigmoid in deep learning through dynamical isometry*, arXiv:1711.04735: §2.3–2.5, eqs. (11)–(19), Supplement Result 1, Example 1.
- B. Hanin, M. Nica, *Products of many large random matrices and gradients in deep neural networks*, arXiv:1812.05994: Thm 1, Prop. 2, Cor. 3, §1.1–1.5.
- F. Benaych-Georges, R. R. Nadakuditi, *The eigenvalues and eigenvectors of finite, low rank perturbations of large random matrices*, arXiv:0910.2120: Thms 2.1–2.3, 2.6–2.8, Props 2.4, 2.9, §2.5, §4.

**Used through the neighbouring digests (read there by parallel runs):** Hastings 2007; Anari–Liu–Oveis Gharan; Chen–Liu–Vigoda; Chen–Eldan; Chen–Rouzé arXiv:2504.02208; Yang arXiv:2609.38007; Dyer–Goldberg–Jerrum.

**From memory, not opened:** Furstenberg–Kesten 1960 and Oseledets 1968 (as cited by Hanin–Nica); Birkhoff 1957 contraction; Ruelle–Perron–Frobenius (Bowen, LNM 470); Kesten 1973 and Bougerol–Picard 1992 on random affine recursions; Voiculescu's asymptotic freeness of independent orthogonally invariant matrices.

**Not retrieved:**
- The ChatGPT share page for arXiv:2609.38007, which is blocked by the egress proxy; see [chat-2609.38007-retrieval-status.md](chat-2609.38007-retrieval-status.md).
- `C:\Users\User\Downloads\expander_survey.pdf`, which is on the user's machine and not in the container.
