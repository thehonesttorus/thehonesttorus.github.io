# Foundations, part 4: measured facts about the object

### What the campaign has measured about deep He-initialised ReLU networks (their laws, their structure, their scaling) and what each fact means for any estimator

*Fresh-slate foundations note, Part 4 (1 Oct 2026). Contract: [BRIEF.md](BRIEF.md). This note compiles the measured facts about the random ReLU networks of this problem that the campaign's streams and digests contain. Each is restated as a statement about the object, not about the estimator that happened to measure it. Each fact carries its conditions (width n, depth L, number of networks, Monte Carlo sample size N, layers), its source, a flag where it was measured only at small width, at the Phase 1 shape or on an estimator's internal state, and what it implies for any estimator: an information constraint, a cost constraint, an opportunity, or nothing. Nothing was re-measured for this note: numbers are quoted from the sources. Section 15 collects the places where two sources disagree, and the places where a statement in circulation (BRIEF §3 included) goes beyond its evidence. Siblings: [foundations-problem.md](foundations-problem.md) (Part 1, the problem from first principles, with its own measurements; cited here as "Part 1", and agreements are recorded inline) and [foundations-unlocks.md](foundations-unlocks.md) (Part 2).*

**Sources, read in full.**
- The plan: [competition-plan.md](../competition-plan.md) §§3.1, 6b, 7.
- Stream reports: [oracle1024](../streams/oracle1024/REPORT.md) (with `results/summary_prelim.md`), [chain128](../streams/chain128/REPORT.md), [old-content](../streams/old-content/REPORT.md), [coef-ensemble](../streams/coef-ensemble/REPORT.md), [submission](../streams/submission/REPORT.md), [costmodel](../streams/costmodel/REPORT.md), [theory](../streams/theory/REPORT.md).
- Digests: [phase2-intel-2026-10-01.md](../digests/phase2-intel-2026-10-01.md) (team EscAI's oracle), [oishi1029-phase1-whitened-mc.md](../digests/oishi1029-phase1-whitened-mc.md), [504aldo-k3-chain.md](../digests/504aldo-k3-chain.md), [bridges/transfer-spectrum-measurement.md](../digests/bridges/transfer-spectrum-measurement.md).

A few definitions were checked in the primary sources behind the digests: EscAI's `research/tenfold_push/README.md` and `research/spectral_push/README.md`, and the Oishi probe notes `analytic-per-layer.md`, `mlmc-surrogates.md`, `jacobian-reuse.md` and `RESEARCH.md`. The bench calibration in [bench/RESULTS.md](bench/RESULTS.md) is quoted where the BRIEF relies on it.

**Citation keys.** plan, oracle1024, chain128, old-content, coef-ens, submission, costmodel, theory, intel (the EscAI material in the phase2-intel digest), Oishi, 504aldo (with its finding numbers F..), transfer (the transfer-spectrum digest), bench, Part 1.

---

## 0. Summary

### 0.1 Five load-bearing statements

1. **What suffices, and where it comes from.**
   - The final readout needs only each final pre-activation's own law to fourth order. Fed the measured (m, s², κ3, κ4), the first-order readout formula has bias 1.2e-10 (n = 1024, 8 networks). The Edgeworth series converges: 6.5e-11, 5.0e-12 and 3.1e-13 at CLT weights 2, 3 and 4 (F3.1).
   - Supplying the per-neuron pre-activation variance, κ3 and κ4 at every layer to the best public representation, whose own joint state stays approximate, cuts its raw MSE 20×, from 2.45e-8 to 1.17e-9 (F3.2).
   - Those marginals are manufactured by pairwise structure. Through each ReLU, the covariance needs the pairwise (2,1), (3,1) and (2,2) cumulant slices, and marginal cumulants do not help it (F6.1). The next (2,1) slice needs the all-distinct third cumulant (about half of it, F6.2) and the (2,1,1) slice of the fourth cumulant (F7.1).
   - The marginal information therefore does not regenerate itself.
2. **Generation at each ReLU is local, tree-shaped and parameter-free at n = 1024.**
   - The new joint third-order content equals the first-order gate diagrams with leg-partition (counting) coefficients. Their error is 0.8–1.0 % of the next (2,1) slice at n = 1024 (one network, preliminary), and it falls as ≈ n^-0.8 at fixed depth (F6.3).
   - The coefficients fitted at n ≤ 160 drift from the counting values. That drift is a finite-width, shape-determined effect, and it is absent at 1024 (F6.4).
   - Every one of these diagrams except one Gaussian triangle is a tree in neuron-index space (F6.7).
   - The (2,1,1) fourth-cumulant slice is indispensable: 3 % error without it, 0.9 % with it. Its rank-one, covariance-shaped regeneration captures a shrinking part of it as width grows: 55 %, 33 %, 20 % at n = 128, 256, 1024 (F7.1).
3. **Memory is large, orthogonal to the present, and of rank proportional to n.**
   - At n = 1024, 38–42 % of the next (2,1) slice at every layer is third-order content born at earlier ReLUs. Its projection on everything built from the current layer has R² ≤ 0.14 (F8.1). It cannot be absorbed into renormalised births (F8.3), and every layer's birth matters at the readout (F6.12).
   - The gated propagators that carry this content have participation ratio ≈ n/(2·age) and no spectral gap (F9.1, F9.2).
   - A 2 % carrier of it needs 0.19–0.5 n modes for ages 15 down to 4, at every width from 128 to 1024 (F9.3; beyond n = 256 this comes from an ensemble formula on fresh networks).
4. **The low-rank structures that remain at every width are one direction and one scalar.**
   - The direction: the mean direction is a Perron outlier of the gated propagators. Its squared overlap with μ_z is 0.85 at n = 1024, its gain is 3–8× the bulk, and it sharpens with width (F2.2). At n = 128 it carries 51–90 % of the slice-containing old content at ages ≥ 6 (F2.3).
   - The scalar: the norm process. At the Phase 1 shape the last layer's ‖y‖² explains 85.6 % of the per-sample output variance (F2.4). At n = 1024 the per-neuron κ3 is largely a common scale mixture (Part 1).
   - The final covariance has stable rank 2.6 in energy, but per-neuron accuracy needs it at ≥ 0.5 n rank (F2.5).
5. **Sensitivity and scaling fix the accuracy and the evidence base.**
   - Final MSE ≈ 4.2e-6·ε² for a relative error ε of the (2,1) slice at every layer. That is the Gaussian-closure MSE times ε², so a 10 % margin over raw 2.1e-8 needs ε ≤ 2.2 %. At n = 128, errors that are random with respect to the true structure cost 4× more per ε² (F10.1).
   - Width scaling is steep but pre-asymptotic below n ≈ 256. The Gaussian-closure distance falls as n^-0.82 from 64 to 128 but as n^-2.0 from 128 to 1024, so a small-width fit overstates the 1024 value 12× (F5.4, Section 13).
   - At the Phase 1 shape (L/n = 1/8) the cumulant expansion is asymptotic, while at the competition shape the readout series converges (F3.5, F10.4).
   - Facts measured only at small width or at the Phase 1 shape are therefore flagged throughout.

### 0.2 What any estimator must respect

| kind | statement | size | facts |
|---|---|---|---|
| information | the final readout needs (m, s², κ3, κ4) of each final pre-activation, nothing more | bias 1.2e-10 | F3.1 |
| information | per-neuron pre-activation variance, κ3, κ4 at every layer carry what the best public representation misses | 20× (2.45e-8 → 1.17e-9) | F3.2 |
| information | non-Gaussian corrections to the final variance | +3.75 %, 27× the 0.14 % tolerance | F3.3 |
| information | pairwise slices, not marginals, set the next covariance and the next (2,1) slice; the all-distinct κ3 is about half of that slice | 44–47 % of D21 (1024) | F6.1, F6.2 |
| information | the transported effect of the (2,1,1) slice of κ4(z) must be captured to ≥ 30 % (≥ 60 % for margin) | closure error 3 % without it, 0.9 % with it | F7.1 |
| accuracy | the (2,1) slice at every layer to ≈ 2–5 % (relative) | 4.2e-6·ε² | F10.1 |
| memory | ≈ 40 % of the (2,1) slice is old, orthogonal to the present, of rank ∝ n | 0.38–0.42; R² ≤ 0.14; 0.2–0.5 n modes | F8.1, F9.3 |
| cost | joint structure shaped like a triangle (a cycle) in index space is unaffordable; trees cost a few products | ≈ 1024 u vs ≈ 1–7 u per term | F6.7 |
| cost | third-order joint structure is invisible to sampling within the budget | D21 to 1 % needs ≥ 1.4e7 samples (B buys 65,536) | F6.10 |
| opportunity | one mean direction (Perron mode) and one norm scalar per layer | overlap 0.85 (1024); 85.6 % (P1) | F2.2–F2.4 |
| opportunity | about half of the last layer's gates are frozen: exactly linear or exactly zero | 54 % outside (0.02, 0.98) | F4.1 |
| opportunity | the gate law is weakly coupled in aggregate | η₀ = 2–8 for up to 1024 gates | F4.2 |
| opportunity | structural constants are properties of the shape (n, L), not of the network | ≤ 0.002 in ε | F11.2 |
| caution | small widths and the Phase 1 shape are pre-asymptotic | up to 12× in projected error | F5.4, §13 |

---

## 1. How to read this note

**Format.** Each fact has an identifier F⟨section⟩.⟨k⟩, a statement in bold, its numbers, and then:
- *Where*: shape (n × L), networks, samples, layers, and the source;
- *Flags*: present only when the fact was not measured at the competition shape on exact laws;
- *Implies*: tagged **[I]** for an information constraint (what an estimator must know, and how accurately), **[C]** for a cost constraint (what it cannot afford or must pay), **[O]** for an opportunity (structure it may exploit), or **[–]** when nothing follows for estimation.

**Flags.**
- **≤256**: measured only at widths ≤ 256, at depth 16 unless stated. Width scaling is pre-asymptotic below n ≈ 256 (F5.4, Section 13).
- **P1**: measured only at the Phase 1 shape, n = 256 and L = 32. There L/n = 1/8, eight times the competition's 1/64; the grader's smoke test also runs this shape.
- **1 net**: one network.
- **prelim**: low sample size, or the production run was not reported.
- **est**: measured on an estimator's internal state (a chain's covariance or cumulant core), not on the exact law.
- **cond**: an information statement obtained by injecting exact quantities into an approximate propagation. It holds conditional on that approximation.
- **ens**: computed from a Gaussian-source ensemble formula on fresh networks, not from exact moments.

**Notation and indexing.**
- z_l = a_{l−1}W_l is a pre-activation, a_l = relu(z_l) and a_0 = x ~ N(0, I_n).
- C is the covariance of the pre-activations and Φ_i = P(z_i > 0) the gate probability.
- Cumulant slices of the pre-activations:
  - D3_a = κ3(z_a, z_a, z_a);
  - D21_ab = κ3(z_a, z_a, z_b), the (2,1) slice;
  - K4_a = κ4(z_a, z_a, z_a, z_a);
  - K22_ab = κ4(z_a, z_a, z_b, z_b);
  - K31_ab = κ4(z_a, z_a, z_a, z_b);
  - K211_abc = κ4(z_a, z_a, z_b, z_c), the (2,1,1) slice.
  - "All-distinct" means the entries with pairwise distinct indices.
- t = m/s is the standardized bias of a pre-activation.
- "Old content" is joint structure born at earlier ReLUs and carried by the linear maps. A "birth" is joint structure created at the current ReLU.
- ε(D21) is the relative Frobenius error of a model of the next layer's D21, noise-corrected unless stated.
- 1 unit = one dense 1024³ product = 2^31 FLOPs, and B = 2^41 = 1024 units.
- Layer indices are quoted as each source gives them. The atlas-based streams (plan §3.1, oracle1024, chain128, coef-ens, old-content, theory, transfer) call the exactly Gaussian first pre-activation xW_1 "layer 0", and their layer l feeds layer l+1 through one ReLU. 504aldo counts L00–L15, Part 1 uses z_1 … z_16, and Oishi counts 1 … 32.

**Exact background used to read the measurements.** These are not measurements; the derivations are in Part 1 §§2–3 and transfer §2.
- **E1. Positive homogeneity.** Every a_l is positively homogeneous of degree 1 in x, so E a_L = E‖x‖ times the spherical average, with E‖x‖ = 31.992 at n = 1024. The input radius alone puts ≈ μ_lμ_lᵀ/2n into every covariance.
- **E2. The ReLU identity.** relu(z) = ½z + ½|z|, so μ_l = ½μ_{l−1}W_l + ½E|z_l|: half of every mean propagates linearly.
- **E3. Quenched and annealed.**
  - W_{l+1} is independent of the law of a_l.
  - Over the weight ensemble, E_W[WᵀSW] = (2/n) tr(S)·I and E_W[W^{⊗3}] = 0. The ensemble-averaged transfer of third-order content is therefore zero, and all transported third-order content belongs to the fixed (quenched) network.
- **E4. Tensor powers add no spectral gap.** B^{⊗k} restricted to Sym^k has singular values ∏ σ_{i_j}. Higher-order transfer therefore has no gap that the matrix product lacks.
- **E5. Every readout sensitivity is a face statistic.** ∂E relu(z)/∂κ_k = E[relu^{(k)}(z)]/k!: the gate probability, half the density at the wall, and the slope and curvature of that density.
- **E6. Order counting.** With ε = n^-1/2, off-diagonal covariances are O(ε). A slice of κ_m is O(ε^{m−2}) when every index multiplicity is even and O(ε^{m−1}) otherwise (checked in F5.2).

---

## 2. The global mode: scale, the mean direction, the norm process

**F2.1 The final layer is mostly a deterministic mean vector.**
- Over neurons, the final-layer means have mean square 0.9095 (rms ≈ 0.95), and the per-neuron variance averages 0.071–0.075.
- So ≈ 92 % of E‖a_16‖² lies in ‖μ_16‖². The infinite-width value is 0.93 (Part 1 §2.5).
- *Where:* 1024 × 16.
  - Mean square: public mini split, 100 networks, truth N = 1e9 (Part 1).
  - Per-neuron variance, four estimates:
    - 0.0748 on the public mini split (EscAI, from the stored mean final variance);
    - 0.0736 in 504aldo F42;
    - 0.0728 on the submission dev set (from its truth-noise floor 3.64e-8 at N = 2e6, 6 networks);
    - 0.0712 on bench w1024_d16 (plain-sampling raw MSE 1.086e-5 × 6,554 samples, 6 networks).
  - 504aldo F62 reports Var(h_L) = 0.079–0.096 on three networks (statistic as reported).
  - At P1 (256 × 32) the mean square is 0.9093 (Oishi).
- *Implies:*
  - [I] The target, an rms per-neuron error of 1e-4, is a relative accuracy of ≈ 1e-4 on O(1) numbers and ≈ 4e-4 of the per-neuron spread.
  - [C] Plain sampling has MSE ≈ 0.0748/N (F12.1).

**F2.2 The mean direction is the one Perron (outlier) mode of the gated propagators.**
- Setup: the gated products are J_{s→t} = W_t D_{t−1} ⋯ D_{s+1} W_{s+1}, where D = diag of the expected gates.
- Their top left singular vector aligns with the target mean direction μ_z(t). Squared overlap:
  - n = 128, first network: 0.09, 0.42, 0.58, 0.63, 0.72, 0.87 and 0.92–0.93 at ages 2, 4, 6, 8, 10, 12 and 14–15;
  - n = 128, second network: 0.72–0.94 at ages 8–15;
  - n = 1024: 0.66 at age 10 and 0.85 at age 15.
- Separation from the bulk:
  - s₁/s₂ rises to 1.5–2.0 at n = 128 and stays near 1.45 from n = 256 on;
  - s₁²/mean s² at age 15 is 62 (atlas network, n = 128), 46 (fresh network, n = 128) and 128 (n = 1024).
- Gain: the mean direction is transported with 1.2–3.3× the bulk gain on one n = 128 network, up to 5.9× on the other, and 3.1–7.8× at n ≥ 256 (ages 8–15).
- Lyapunov spectrum: the finite-time exponents of J_{0→15} at n = 128 are +0.059, +0.030, +0.007, and negative from the fourth on. The top one is the mean direction.
- *Where:* transfer §3.1, §3.4, §3.6. At n = 128: two networks, gates from atlases with N = 5e5. At n = 256–1024: one fresh network per width, gates from 2e4 Monte Carlo inputs.
- *Flags:* ≤256 for the atlas numbers; at 1024, 1 net with Monte Carlo gates.
- *Implies:*
  - [O] One direction per layer, close to μ_z and so cheap to track at O(n²), carries a coherent and amplified part of everything transported.
  - [I] Errors aligned with it are amplified, not damped. Part 1 §3.5 finds the input-averaged propagator's top output direction is the final mean direction (overlap 0.6–0.99 at n = 64–128).

**F2.3 The mean mode carries most of the slice-containing old third-order content at depth, and little of its all-distinct part.**
- Measure: the rank-one projection of transported third-order content onto the top propagator direction u₁(J), counted as a share of that content's own D21 energy.
- For the full source (κ3(a_s) including its index-coincident slices) the share is:
  - 17 % at age 4, 33 % at 5, 51 % at 6, 58 % at 7;
  - 65 % at 8, 75 % at 9, 78 % at 10;
  - 85–90 % at ages 11–15.
- For the all-distinct part alone it is 0–22 % up to age 9, and 3–67 % at ages 10–15, where that part is only 1–6 % of D21.
- A different measure: the old-content stream finds that the leading tensor mode of the aged pool (ages > 4, slices included) lies along the mean direction (cos ≈ 0.93). It holds 60–80 % of the tensor energy at depth but carries "almost no D21".
- *Where:*
  - transfer §3.1 and §3.6: n = 128, one network, noise-corrected with two atlases of N = 5e5;
  - old-content: n = 128, atlases A1, A2, B1 with N = 6e5.
- *Flags:* ≤256; 1 net for the transfer figures. The two measures read differently (Section 15, T1).
- *Implies:*
  - [O] Possibly the strongest single structure available for old content: one direction per layer.
  - [I] How much D21 it carries at n = 1024 is not measured (Section 16).

**F2.4 At the Phase 1 shape, the per-sample output fluctuation is dominated by one scalar, the layer norm [P1].**
- The last layer's ‖y‖² alone explains 85.6 % of the per-sample variance of the output. The first layer's norm explains 21 %.
- The input-averaged Jacobian E[J] = E[f xᵀ] has 79 % of its Frobenius energy in one direction (203× the isotropic share) and 97 % in eight.
- *Where:* Oishi (mlmc-surrogates and active-subspace probes); 256 × 32.
- *Flags:* P1.
- *Competition-shape counterpart (Part 1 §3.3, one network at n = 1024):* the per-neuron third cumulant is largely a common scale mixture.
  - λ3,i ≈ b_l t_i, with R² = 0.52, 0.78, 0.87 and 0.89 at z2, z4, z8 and z16.
  - λ4 does not depend on t.
  - The common-scale part is 30–46 % of the λ3 readout term and 58–80 % of the λ4 term.
- *Implies:*
  - [O] One scalar per layer, the norm process, accounts for most of the per-sample fluctuation and for a large coherent part of the per-neuron non-Gaussianity. The input radius is its exactly known first instance (E1). Carrying its law costs O(1) per layer.
  - [I] It does not carry the neuron-specific remainder.

**F2.5 The final covariance is concentrated in energy, but per-neuron accuracy needs it at near full rank.**
- At n = 1024, the final pre-activation covariance as carried by a V29-derived chain (8 networks) has:
  - stable rank Σλ²/λ₁² = 2.57;
  - entropy effective rank 16.7;
  - participation (Σλ²)²/Σλ⁴ = 5.9;
  - 25 eigen-directions for 90 % of Σλ² and 79 for 99 %, with 6.7 % of Σλ² beyond rank 32.
- At the Phase 1 shape (exact covariances), the variance ranks r90/r99/r99.9 of the activation covariance are:
  - 171/244/255 at layer 0;
  - 58/162/206 at layer 7;
  - 31/115/168 at layer 15;
  - 14/66/126 at layer 23;
  - 8/47/111 at layer 31.
- Yet truncating the covariance of layers ≥ 16 to rank 48, 64 or 128 leaves a final bias 211×, 128× or 9.6× the 4.2e-7 noise floor; only rank 192 gets below it. Energy rank and accuracy rank differ by ≈ 4×.
- *Where:*
  - intel §3.4 (EscAI `spectral_push`). Its README says the figures describe the estimator's covariance, not an exact reference. The chain's covariance sketch error at layer 14 is ≈ 0.7 % (intel §3.2).
  - Oishi §§7.2, 8 (P1).
- *Flags:* est (1024); P1 (ranks).
- *Agreement:* Part 1 §4 measures the same split at n = 1024. A 0.1 % truncation needs rank ≥ 0.93 n at a2, falling to 0.47 n at a15.
- *Implies:*
  - [I][C] The per-neuron quadratic forms w_iᵀΣw_i need a second-order object at ≥ 0.5 n rank; truncation is not available.
  - [O] The energy concentration (mean spike, F2.2) is available for representation, not for truncation.

---

## 3. Per-neuron laws and the readout

**F3.1 The final readout needs each final neuron's own law to fourth order, and the Edgeworth series converges at n = 1024.**
- Inputs fed to the readout E relu(z) for each final pre-activation, and the bias MSE that remains:
  - measured mean and variance only (Gaussian readout): 2.03e-7;
  - adding the measured κ3 and κ4 in the first-order formula: 1.17e-10 (95 % CI 1.08e-10 to 1.31e-10);
  - the complete Edgeworth expansion to CLT weight 2: 6.5e-11;
  - to weight 3, with the true κ5: 5.0e-12;
  - to weight 4, with the true κ6: 3.1e-13.
- *Where:* intel §3.2 (EscAI `marginal_oracle`); 1024 × 16; 8 networks; two independent N = 1e8 replicas per network; cross-replica bias estimate.
- *Agreement:* Part 1 §3.1 (one network, 1.2e6 samples) gives 1.7e-7 Gaussian, ≤ 2.9e-10 first order and ≤ 1.2e-10 second order.
- *Implies:*
  - [I] The final layer requires (m, s², κ3, κ4) of every final pre-activation and nothing else. κ3 and κ4 enter at first order with loose tolerances: ≈ 12 % of the λ3 term and ≈ 40 % of the λ4 term (Part 1 §3.2).
  - [–] The readout formula is not a bottleneck for any target above ≈ 1e-10; the difficulty lies entirely in its inputs.

**F3.2 The per-neuron (marginal) pre-activation laws at every layer carry what the best public representation misses [cond].**
- Experiment: at every layer of a V29-derived K = 3 chain, selected per-neuron pre-activation quantities are replaced by reference values while the means propagate.
- Raw MSE, starting from the chain's 2.45e-8 on networks 0–2:
  - variance replaced: 1.47e-8 (−40 %);
  - variance and κ3: 8.35e-9;
  - variance and κ4: 8.92e-9;
  - post-activation variance only: 2.53e-8, no gain;
  - variance, κ3 and κ4, with the exact Gaussian layer 0 kept, on networks 0–7: 1.165e-9. That is a ratio of 0.0497 (95 % CI 0.046–0.055), and every network improves.
- A linear "deferred-mean" correction of the readout, propagated through Wᵀ with the trajectory held fixed, reaches the same 20×: 2.344e-8 → 1.173e-9 on 8 networks.
- *Where:* intel §3.2 (EscAI `tenfold_push`); 1024 × 16; public mini networks 0–2 for the partial rows and 0–7 for the full row. References come from the public supplemental bake (N = 1e8 per network, so they are noisy).
- *What was replaced* (EscAI README):
  - "κ3" replaces the per-neuron third cumulant, not the joint tensor.
  - "Joint κ4" replaces the per-neuron fourth cumulant and the two chain inputs derived from it: the input of the chain's κ4→κ3 births, and its radial (2,2) model. It does not replace the joint fourth-order tensor.
- *Flags:* cond. The partial rows use 3 networks, and in them the exact layer-0 variance was also replaced by a noisy reference.
- *Implies:*
  - [I] Given the per-neuron variance, κ3 and κ4 of every pre-activation (n numbers each per layer), even a representation whose own joint state is approximate reaches ≈ 1e-9, an order of magnitude below the bar. Producing those marginals is the problem (F3.3, F6.1, F6.3).
  - [I] Correcting post-activation second moments, without the pre-activation shape, carries nothing.

**F3.3 Second moments propagated as Gaussian miss the true final variance by 3.75 %, and correcting the variance alone makes the means worse.**
- Against Monte Carlo, a Gaussian (covariance) propagation underestimates the per-neuron pre-activation variance by 0.04 % at layer 1 and by 3.75 % at layer 16.
- The final pre-activations have excess kurtosis +6.6 % and "std-skew ~1e-3" (as reported).
  - Reading: the skew figure is presumably the signed average over neurons, which is small because each neuron's skewness takes the sign of its t (Part 1 §3.3).
  - The rms standardized skewness is ≈ 0.10 at n = 1024 (Part 1).
- Forcing the true variance into that propagation raises the final MSE from 4.79e-6 to 6.92e-6.
- *Where:* 504aldo F60; 1024 × 16; public mini networks; Monte Carlo with 3e5 samples.
- *P1 counterpart* (Oishi `diag.py`, a different Gaussian-propagation implementation from F3.6's; 256 × 32):
  - the Gaussian propagation is 0.83 % off in the mean at layer 32;
  - it is 11.4 % off in σ at layer 32, against 0.09 % at layer 1;
  - the Gaussian readout fed the true marginal (μ, σ) still errs by 4.06e-6, against 5.96e-5 for the full Gaussian propagation.
- *Implies:*
  - [I] The non-Gaussian contribution to the final variance (3.75 %) is ≈ 27× Part 1's per-neuron tolerance of 0.14 %. The covariance must therefore be propagated with its non-Gaussian corrections, which Part 1 §3.4 traces to the (2,1), (3,1) and (2,2) slices.
  - [I] In a second-moment description, variance and shape errors compensate each other. Variance information helps only together with κ3 and κ4 (consistent with F3.2).

**F3.4 The non-Gaussianity of the marginals grows linearly in depth, and it is small at the competition shape.**
- n = 1024: excess kurtosis +6.6 % at layer 16 (504aldo F60). Part 1 §5 (one network) gives rms λ3 = 0.10 and λ4 = 0.056 at z16, with rms λ3 ≈ (5–7) l/n and λ4 ≈ (3.6–6) l/n from n = 128 to 1024.
- n = 128: kurtosis 0.47–0.64 at layers 10–14 (theory, 2 networks).
- n = 32, depth 6: kurtosis up to 1.0, with correlations 0.21–0.36 (theory §4).
- P1, 256 × 32: z32 excess kurtosis +0.30 to +0.44 over 5 networks, with skew ≈ 0 on average (Oishi).
- *Implies:*
  - [I] At the competition shape the corrections are small (λ ≤ 0.1), so first order in 1/n suffices at the readout.
  - [I] At n ≤ 64 and at the smoke-test shape λ reaches 0.3–1, and an expansion in 1/n is no longer small there (F3.5, F10.4).

**F3.5 At L/n = 1/8, a marginal's cumulant series is asymptotic, while a max-entropy density with 6–8 moments converges [P1].**
- Readout of the final layer from its exact marginal moments.
- Gram–Charlier/Edgeworth MSE by order (minimum at order 4):

  | order | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 |
  |---|---|---|---|---|---|---|---|---|
  | MSE | 1.11e-6 | 1.61e-7 | 2.30e-8 | 3.07e-8 | 2.73e-8 | 3.75e-8 | 4.9e-7 | 5.1e-6 |

- Maximum-entropy densities with M moments: 5.26e-8 at M = 4, 5.08e-9 at M = 6, 8.15e-10 at M = 8, and 1.70e-10 at M = 10 (the noise floor).
- The sharp worst case over all laws sharing the same M moments (Markov–Krein) is 3.43e-4 at M = 2 down to 1.29e-5 at M = 12, decaying as M^-1.65. The actual laws sit ≈ 2,500× inside it.
- Wasserstein distance to a Gaussian: summed over neurons, the squared W₂ distance between each final activation's law and relu of a Gaussian with the same (μ, σ) is 0.86–1.15 % of the total per-sample variance.
- *Where:* Oishi, deterministic-closure probe (32 networks, 8,192 neurons) and mlmc-surrogates probe (5 networks, 65,536 samples); 256 × 32.
- *Flags:* P1.
- *Implies:*
  - [I] At the smoke-test shape and at narrow widths, cumulant truncation of a marginal diverges past fourth order. At the competition shape it converges (F3.1).
  - [O] An exponential-family (max-entropy) description of each marginal converges where the cumulant series does not.
  - [I] The per-neuron laws are smooth and near-Gaussian, far from the worst case that their moments allow.

**F3.6 Under Gaussian propagation the mean error saturates near 1 % with depth while the σ error grows, and the spread of t = μ/σ grows steadily [P1].**
- Gaussian closure against exact means, per layer (40 networks):

  | layer | 1 | 2 | 4 | 8 | 16 | 32 |
  |---|---|---|---|---|---|---|
  | relative rms error of the mean | 0.0047 % (= reference noise) | 0.227 % | 0.446 % | 0.652 % | 0.792 % | 0.827 % |
  | relative error of σ | 0 % | 0.21 % | 1.0 % | 3.3 % | 9.0 % | 16.1 % |
  | rms μ/σ of the pre-activations | 0 | 0.68 | 1.23 | 2.00 | 3.11 | 4.28 |

- The one-dimensional Gaussian readout has an irreducible error of 0.223 % (layer 2), 0.259 % (layer 4) and 0.112 % (layer 32), against 0.0105 % needed.
- *Where:* Oishi `analytic-per-layer` (definitions checked); 256 × 32; 40 networks.
- *Flags:* P1.
- *The spread of t is a mean-field quantity* (Part 1 §5, one network at n = 1024): std(t) = 0.69, 1.24, 2.10, 3.02 and 3.70 at z2, z4, z8, z12 and z16. The infinite-width values are 0.68, 1.24, 2.06, 2.78 and 3.46. At n = 128 the values are only 0.64–2.05.
- *Implies:*
  - [I] A description by second moments saturates at ≈ 1 % relative error of the means at depth, two orders of magnitude above the target.
  - [O] At n = 1024 the spread of standardized biases, and with it the fraction of frozen gates (F4.1), follows mean-field theory to ≈ 7 %.

**F3.7 Readout sensitivities [P1]; their competition-shape values are in Part 1.**
- Final MSE ≈ 0.048·(δμ/σ)² for an error δμ in the pre-activation mean.
- A relative perturbation of σ by 1e-3 costs 4.26e-9; by 1e-2, 8.10e-8.
- Perturbing the standardized cumulants c3–c6 by 1e-2 costs 3.65–3.68e-9.
- Projecting relu onto polynomials of degree ≤ 6 keeps 99.2 % of its variance.
- Competition shape: Part 1 §3.2 gives the tolerances for MSE 2.5e-9 at n = 1024 (one network): 5.6e-5 relative on the mean, 0.14 % on the variance, ≈ 12 % on the λ3 term and ≈ 40 % on the λ4 term.
- *Where:* Oishi; 256 × 32.
- *Flags:* P1.
- *Implies:* [I] The readout is most sensitive to the mean (through the gate probability), then to the variance (through the wall density). The third and fourth cumulants enter weakly.

---

## 4. Gates and faces

**F4.1 About half of the last layer's gates are frozen, and a frozen neuron is exactly linear or exactly zero.**
- 1024 × 16:
  - One fresh network per width, gates from 2e4 Monte Carlo inputs: the number of uncertain gates (0.02 < P(z > 0) < 0.98) falls from 1024 at layer 0 to 472 at layer 15. So 54 % of the last layer's gates are frozen (transfer §3.5).
  - Part 1 §2.5 (bench network, 1.2e6 samples): 3.5 %, 15 % and 26 % of the neurons of z8, z12 and z16 never changed gate.
- n = 128 (two networks, N = 5e5): at layers ≥ 8, 18–44 % of gates have P < 1e-3 or P > 1 − 1e-3, and 32–61 % lie outside (0.02, 0.98) (transfer N5).
- P1 (256 × 32), from Oishi:
  - At the last layer, 25–29 % of neurons are always on in the sample and 27–32 % never fire in 1,024–32,768 samples.
  - The fraction of neurons that fire at least once in a block of 1,024 samples falls from 1.00 at layer 0 to 0.72 at layer 30.
  - Averaged over neurons, exactly half of all activations are non-zero at every layer.
  - The always-on neurons (about a quarter) carry 53–73 % of a sampling estimator's error; the never-firing ones carry none.
- *Implies:*
  - [O] A frozen-on neuron is exactly linear (a = z) and a frozen-off neuron is exactly zero, up to O(sφ(t)/t²).
    - That remainder is below 1e-7 for |t| > 4.75, which holds for 22 % of z16 at n = 1024 (Part 1).
    - The 54 % outside (0.02, 0.98) are nearly so.
    - A large part of the last layer is therefore a known linear map of the previous layer's mean, or zero.
  - [I] Frozen-on neurons pass the previous layer's mean errors with gain 1, and they carry the largest per-neuron variance.
  - [C] On average there is no sparsity to exploit (density 0.5).

**F4.2 The face (gate) law is weakly dependent in aggregate: its unpinned spectral independence is 2–8, and boundedness is not established.**
- Definition: η₀ = λ_max(Cor) − 1, where Cor is the correlation matrix of the uncertain gates (0.02 < p < 0.98).
- n = 128 (two atlas networks): η₀ = 1.69–1.78 at layer 0, equal to the exact Sheppard value, then 2.5–5.6 at layers 2–15.
- n = 1024 (fresh networks, 2e4 samples): η₀ = 2.1, 3.8, 5.1, 6.3, 6.9, 7.3, 7.0 and 7.6 at layers 0, 2, 4, 6, 8, 10, 12 and 15. The Monte Carlo bias is ≈ +0.2.
- Growth with width at depth:
  - layer 9: 4.2, 6.1, 7.5 and 7.5 at n = 128, 256, 512 and 1024;
  - layer 13: 3.7, 4.5, 5.6 and 7.9;
  - it levels off at layers 9–12 and is still rising at layers 13–15.
- Among the uncertain gates the mean |Cor_ij| is 0.08–0.10 at depth (n = 128). A zero-mean Gaussian surrogate with the same correlations gives η₀ = 15–29 at depth: the mean shift that freezes gates decorrelates the remaining ones.
- Pinned values (conditional on some gates) were not measured; they need third-order gate statistics.
- *Where:* transfer §3.5.
- *Flags:* ≤256 for the atlases; ens/Monte Carlo gates at 512–1024.
- *Implies:*
  - [O] At depth the gate law of a layer behaves as a weakly coupled field (η₀ much smaller than the number of uncertain gates), which favours local or mean-field treatments of the gates.
  - [I] The hypothesis of the local-to-global mixing theorems (bounded η under all pinnings) is unverified (Section 15, C1).

**F4.3 With depth, the gates of independent inputs agree, and the face partition is far finer than any covering [P1].**
- Activation-pattern agreement between two inputs with correlation ρ:
  - for independent inputs (ρ = 0): 0.498 (layer 1), 0.805 (layer 8), 0.870 (layer 16), 0.921 (layer 32), and 0.840 over all layers;
  - at layer 1: 0.850 at ρ = 0.9, 0.951 at ρ = 0.99, 0.986 at ρ = 0.999.
- Gate flips across all 8,192 gates under a relative input perturbation t: 0.3 at t = 1e-4, 26.6 at t = 1e-2, 923 at t = 1.
- Linearising the network at a hub input x₀ gives relative errors 1.632, 0.557, 0.156 and 0.037 at ρ = 0, 0.9, 0.99 and 0.999. Var(f)/E[f²] = 0.041, so the linearisation beats predicting zero (relative error < 0.20) only within ρ > 0.995, i.e. 5.7°.
- Covering the sphere by caps: in d = 256 this needs ≈ 1e131 caps of 18.2°, 1e218 of 8.1° and 1e346 of 2.6°.
- *Where:* Oishi `jacobian-reuse`; 256 × 32; 6–8 networks.
- *Flags:* P1.
- *Implies:*
  - [C] Any enumeration, or piecewise-linear covering, of the faces is out.
  - [O] At depth most gates agree across independent inputs (87–92 % at layers 16–32). The input-dependent face structure sits on a nearly input-independent skeleton (consistent with F4.1).

**F4.4 Gates open where the mean pushes up: generic directions are damped and the mean direction is amplified.**
- The mean-square gain of a generic direction through D·A, 2E[Φ²], is 0.5 at layer 0 and 0.58–0.95 at layers 1–15, with one exception of 1.07 (n = 128, two networks).
- Because Φ_{l+1} = Φ(A_{l+1}μ_{a,l}/σ), the gates rectify the mean (F2.2).
- *Where:* transfer §2.6.
- *Flags:* ≤256.
- *Implies:* [I] Errors in generic directions shrink at every layer by 0.58–0.95 in power, while errors along the mean grow (F2.2, F10.3).

**F4.5 Exact per-neuron gate probabilities add no information to a K = 3 description [cond].**
- Supplying the measured gate probabilities to the gate transport of a V29-derived chain gives 2.4516e-8, against 2.4501e-8 without them (3 networks, 1024 × 16).
- *Where:* EscAI `tenfold_push` README (the primary source behind intel §3.2).
- *Flags:* cond.
- *Implies:* [I] Per-neuron gate probabilities are not where the missing information is. At n = 1024 their Gaussian values are accurate to ≈ 1e-3 at the last layer (Part 1 §3.1). The pairwise and joint gate statistics are a different matter (F6.8).

---

## 5. Pairwise (second-order) structure

**F5.1 Correlations between distinct neurons are O(n^-1/2) and grow about threefold with depth.**
- At layer 1, the rms off-diagonal correlation is 0.211 at n = 32 and 0.109 at n = 128 (ratio 1.94, against the predicted 2).
- At n = 128 it is 0.11 at layer 1, 0.25 at layer 7 and 0.29–0.38 at layers 10–14 (two networks, N = 5e5).
- At n = 32 and depth 6 it is 0.21–0.36.
- Part 1 §5 reports √n × rms correlation = 0.74 (a1), 0.87 (a2), 1.1–1.3 (a4) and 1.4 (a8) at n = 32–128. The plan §3.1 estimates ρ ≈ 0.03 at n = 1024.
- *Where:* theory §1.3, §4; Part 1 §5.
- *Flags:* ≤256. The n^-1/2 law was checked at layer 1 only.
- *Implies:* [I] The small parameter of every pairwise expansion is ρ ≈ (1 to 3)·n^-1/2, larger at depth. At n = 1024 it is a few per cent, so first order is accurate, but coherent sums over n terms are not small (Part 1 §3.3).

**F5.2 The order in n^-1/2 of each joint cumulant slice is as counted.**
- Measured at layer 1 (σ-normalised rms), at n = 32 and n = 128:

  | object | n = 32 | n = 128 | ratio | predicted |
  |---|---|---|---|---|
  | ρ (off-diagonal correlation) | 0.211 | 0.109 | 1.94 | 2 |
  | all-distinct κ3 | 0.084 | 0.020 | 4.2 | 4 |
  | D21 | 0.138 | 0.033 | 4.2 | 4 |
  | D3 | 0.288 | 0.071 | 4.1 | 4 |
  | K22 | 0.135 | 0.027 | 5.0 | 4 |
  | K4 | 0.354 | 0.077 | 4.6 | 4 |
  | K211 | 0.050 | 0.006 | 8.3 | 8 |
  | K31 | 0.102 | 0.013 | 7.8 | 8 |

- *Where:* theory §1.3 (claim C4, upheld by two independent verifiers).
- *Flags:* ≤256, layer 1 only. Deeper layers were not tested.
- *Implies:*
  - [I] The target D21 sits at order ε² = 1/n. The κ3 slices (D21, D3) and the even κ4 slices (K22, K4) enter at that order.
  - [I] The (2,1,1) and (3,1) slices of κ4 enter at ε³, i.e. at relative order ε ≈ 3 % of D21 at n = 1024. That is above the 2.2 % bar, which is why the (2,1,1) slice decides the interface (F7.1, F10.1). κ5 enters at ε⁴.
  - [C] A first-order description needs n² objects (the slices) plus the n³ slice K211.

**F5.3 Covariance-type content from earlier layers concentrates on few propagator modes; the current covariance is mostly young.**
- The symmetric-matrix transfer of products z_s → z_t (n = 128, one network) has participation ratio on Sym² of:

  | age | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 15 |
  |---|---|---|---|---|---|---|---|---|---|
  | PR(Sym²) | 1031 | 358 | 173 | 94 | 37 | 15 | 8.0 | 6.6 | 5.9 |

- The transported covariance J C(a_s) Jᵀ is 30–100 % of C(z_t), less at older ages.
- Projected on the top-k left singular vectors of J, its relative error is:
  - at k = 16: 0.28 (age 4), 0.087 (age 8), 0.032 (age 12);
  - at k = 32: 0.12, 0.015 and 0.002.
- One layer's map S ↦ BSBᵀ is strongly non-normal: its operator norm s₁² is 1.9 at layer 0 rising to 3.5–5.4, while its spectral radius |l₁|² is only 0.51–1.19.
- *Where:* transfer §3.3.
- *Flags:* ≤256, 1 net.
- *Implies:*
  - [O] Second-order (Sym²) content from earlier layers is represented well by a few propagator modes. It needs the square of the spectral energy fraction, not its cube as third-order content does.
  - [I] A layer's own covariance is a poor basis for its old third-order content (F9.4).

**F5.4 At n = 1024, the non-Gaussian corrections accumulated in the pairwise law are ≈ 400× the target in MSE, and their width law changes regime between n = 128 and 256.**
- The distance between the true final means and a Gaussian (covariance-closure) propagation is:
  - n = 1024: raw MSE 4.30e-6 ± 0.36e-6 (6 bench networks, truth N = 2e6, noise subtracted), and 4.39e-6 (range 2.5–5.9e-6) on 8 public mini networks with truth N = 1e9 (504aldo F43);
  - n = 128: 2.89e-4 (8 networks); n = 64: 5.1e-4 (8 networks).
- That is n^-0.82 between 64 and 128 but ≈ n^-2.0 between 128 and 1024. A two-width fit at 64–128 predicts 5.2e-5 at 1024, 12× too high.
- *Where:* bench (cited in BRIEF §§1, 4); 504aldo F43.
- *Implies:*
  - [I] A factor ≈ 400 in MSE (20 in rms) separates the Gaussian description from the target at n = 1024.
  - [I] Widths ≤ 128 are in a different regime, so projections to 1024 must be anchored at n ≥ 256.

---

## 6. Third-order joint structure

**F6.1 Within the K = 3 description, the next ReLU sees third-order structure only through D3 and D21.**
- An otherwise memoryless K = 3 propagation is fed the exact D3 and D21 of every pre-activation. It then reproduces the full K = 3 propagation's final MSE exactly.
- Final MSE with partial information (the full propagation gives 4.32e-8):
  - exact D21 alone: 4.4e-7;
  - exact D3 alone: 4.7e-7;
  - D3 set to 0: 3× worse than a wrong D3.
- Part 1 §3.4 finds the same from the covariance side. The covariance through the gate needs the (2,1) slice, and the (3,1) and (2,2) slices of κ4. Marginal cumulants do not improve it at all.
- *Where:* 504aldo F65; 1024 × 16; 8 public mini networks with truth N = 1e9. Part 1: widths 16–128.
- *Implies:* [I] Per layer, the next nonlinearity needs n + n² third-order numbers (D3 and D21). But these are contractions of the full n³ third cumulant of the previous post-activation, which therefore has to be generated or carried implicitly (F6.2).

**F6.2 The all-distinct third cumulant of the post-activation makes up about half of the next (2,1) slice, at every width.**
- Keeping only the index-coincident slices of κ3(a_l) and dropping its all-distinct part leaves ε(D21(l+1)) =
  - 24–60 % at n = 32 (layers 0–3);
  - 43–70 % at n = 128 (depth 16);
  - 44–47 % at n = 1024 (layers 3–14).
- The width trend is flat: the rms over layers 1–6 is 0.561, 0.517 and 0.491 at n = 32, 64 and 128 (ε ∝ n^-0.10).
- *Where:*
  - plan §3.1: n = 32 with one network and N = 4e5; n = 128 with two networks and two atlases of N = 5e5 each;
  - oracle1024: one network, N = 32,768 × 2;
  - chain128 §1: one network, pair atlases.
- *Flags:* 1 net at 1024.
- *Implies:*
  - [I] No truncation by index pattern works: the n³ all-distinct entries have to be represented.
  - [C] They cannot be represented explicitly: one contraction pass costs ≥ 1024 units (Part 1 §4). They must be produced from O(n²) objects at O(n³) cost.

**F6.3 At n = 1024, the new joint third-order structure at each ReLU is the first-order gate diagrams with counting coefficients, to ≈ 1 % of D21, with no fitted constant.**
- To a few per cent, the all-distinct κ3(a) is the Gaussian ρ² (Wick) term plus these first-order gate diagrams:
  - B0: Φ³κ3(z), the pass-through of the current pre-activation's all-distinct third cumulant;
  - B1, B2: the D21(z) hyperedge with one covariance edge, in two orientations;
  - B3: the D3 dressing;
  - B4: the Gaussian ρ³ term;
  - B5: the K22 hyperedge with a covariance edge;
  - B6: the (2,1,1) κ4(z) hyperedge.
- Their coefficients come from leg-partition counting.
- ε(D21(l+1)) of this closure at n = 1024:
  - 0.0–0.9 % at layers 3–5, 0.8–1.0 % at 6–9 and 0.84–0.95 % at 10–14;
  - fitting the coefficients per layer gains ≈ 0.1 point (0.7–0.85 % at 10–14);
  - the fitted coefficients equal the counting values at every layer: B1 2.90–3.13, B2 2.89–2.99, B5 1.25–1.50 and B6 1.40–1.49, against 3, 3, 1.5 and 1.5.
- Width trend at layers 7–14:

  | width | first-order closure | Wick term alone |
  |---|---|---|
  | 128 | 3.7–7.1 % | 7.0–10.3 % |
  | 256 | 2.0–2.9 % | 6.6–8.4 % |
  | 1024 | 0.81–0.98 % | 4.0–4.5 % |

  The closure falls as ≈ n^-0.8 (deep layers 4.7 → 2.6 → 0.9 %) and the Wick term alone as ≈ n^-0.25.
- Adding the Gaussian ρ³ and ρ⁴ Hermite terms changes ≤ 1 point at every width: what the Wick term misses is non-Gaussian.
- *Where:*
  - oracle1024: at 1024 one network with N = 32,768 × 2 replicas (replica cross-product estimate); at 128 and 256 two networks with N = 5e5 and 1e6;
  - plan §3.1: widths 32 and 128, dense atlases;
  - chain128 §1: widths 32–128.
- *Flags:* 1 net and prelim at 1024; the production runs at N = 3.5e6 were not reported. The width exponent rests on one or two networks per width.
- *Implies:*
  - [O] At the competition width, the law that generates joint third-order content is local, short and parameter-free. It is a fixed list of index-tree diagrams on (C, D21(z), D3(z), K22(z), κ4(z)_(2,1,1), κ3(z)), weighted by per-neuron gate and wall weights.
  - [I] It needs the (2,1,1) slice of κ4(z) (F7.1) and the pass-through of old content (F8.1).
  - [–] Fitted per-layer coefficient tables are unnecessary at n = 1024.

**F6.4 At widths ≤ 160 the coefficients are renormalised, but the renormalisation is an ensemble property of the shape and decays with width.**
- Leg-partition (counting) closure, and the same closure with an ensemble table of fitted coefficients. ε_rep by layer band (layers 1–3 / 4–9 / 10–14), noise-corrected:

  | width (networks) | counting coefficients | ensemble table |
  |---|---|---|
  | 64 (6) | 0.045 / 0.074 / 0.093 | 0.038 / 0.046 / 0.031 |
  | 96 (6) | 0.036 / 0.049 / 0.068 | 0.034 / 0.030 / 0.029 |
  | 128 (8) | 0.031 / 0.038 / 0.055 | 0.030 / 0.025 / 0.018 |
  | 160 (3) | 0.023 / 0.029 / 0.023 | 0.022 / 0.018 / 0.012 |

- The fitted coefficients drift with depth:
  - B2 from 3 to 2.3–2.6;
  - B6 from 1.5 to 1.0–1.2;
  - B3 from ≈ 0.4 at layer 1 to ≈ 0 at depth (counting value 1, corrected to 0.5 in F6.5).
- Only B2, B3 and B6 matter in ablation.
- The drift is a property of the shape:
  - a held-out ensemble table equals the per-network fit to ≤ 0.002 in ε;
  - a table fitted at one width loses ≤ 0.002 at any other width in 64–160;
  - a quadratic-in-depth table does as well as a per-layer one.
- The drift decays with width as n^-q, with q ≈ 0.4–1.0 for B1, B2, B4, B5 and B6, and q ≈ 0.1 for B3. At 1024 it is absent (F6.3).
- Second-order diagrams explain part of the drift (theory, n = 128, two networks). With them fixed:
  - B6 returns to 1.45–1.53 at every layer;
  - B3 returns to 0.50–0.82, and B4 to 0.65–1.08;
  - B1 and B2 overshoot (3.35–4.70 and 3.46–3.89).
- The series therefore behaves as an alternating asymptotic one at n = 128. A free fit over ≈ 35 diagrams reaches 1.29–1.74 %.
- *Where:* coef-ens (6, 6, 8 and 3 networks at n = 64, 96, 128 and 160; N = 5e5 per atlas; pairs for noise); theory §3.
- *Flags:* ≤256.
- *Implies:*
  - [I] At the competition width the first-order description stands on its own. A table fitted at small width encodes finite-width effects that vanish by 1024.
  - [O] Anything that is renormalised is a property of (n, L) and can be computed offline.

**F6.5 Two corrections to the first-order list (n = 128).**
- **The D3 dressing B3 has coefficient 0.5, not 1.** The term is symmetric under swapping j and k (|Aut| = 2). Fixing this alone lowers the fixed first-order error at layers 10–14 from 4.15 to 3.08 % on one network and from 4.48 to 3.35 % on the other. At n = 1024 the fitted B3 is 0.17–0.46 (oracle1024), which is consistent.
- **One ε³ diagram is missing.** B7 = 3 sym3[w2_i w2_j Φ_k κ3_ijk C_ij]. It is 2.6–3.3 % of ‖D21‖ at layers 1–5 and 0.6–1.3 % at layers 8–14. Adding it lowers the layer 1–3 error by ≈ 0.55 points (4.40 → 3.84 % and 4.27 → 3.72 %).
- *Where:* theory, claims C2 and C3, each upheld by two independent verifiers; n = 128, two networks.
- *Flags:* ≤256. The 1024 numbers of F6.3 predate both corrections.
- *Implies:* [I] The complete first-order list has eight terms beyond the Wick term. Their coefficients are B0 1, B1 3, B2 3, B3 0.5, B4 1, B5 1.5, B6 1.5 and B7 3.

**F6.6 At small width, second-order and even κ5 terms are several per cent of D21.**
- Shares of ‖D21‖ by diagram class (n = 128, layers 8–14, one network):
  - B0, the pass-through: 53–64 %;
  - the leaf class (all diagrams with a vertex hanging on a single covariance edge): 13–16 %;
  - B6: 6.4–10 %; B2: 6.3–8.5 %;
  - S5c (κ4_iiij with a covariance edge): 2.0–3.7 %;
  - S4d (D21_ki D21_kj): 1.6–2.6 %.
- Adding S4d alone to the corrected first order lowers the error from 2.29–3.80 % to 1.71–2.53 %. The full fixed second order (2.3–3.75 % at layers 10–14) does worse than either, because the leaf terms overshoot.
- At n = 32 and depth 6, the κ5 slices reach 9.6 % (the (3,1,1) slice) and 5.0 % (the (2,2,1) slice) of D21, while κ6 (2,2,2) stays ≤ 0.6 %. There no fixed closure beats a fitted one at layers 3–4: the expansion has stopped converging.
- The transported Monte Carlo noise of each second-order term is ≤ 0.14 % of D21.
- *Where:* theory §§3–4; n = 128 with two networks and N = 5e5 × 2; n = 32 with N = 4e6 × 2.
- *Flags:* ≤256.
- *Implies:* [I] Dropping second-order structure at n = 1024 is justified by the n-scaling (F6.3, F5.2), not by the small-width sizes. At the smoke-test shape and at n ≤ 64, the cumulant expansion of the joint law is not usable.

**F6.7 Every leading diagram except one is a tree in neuron-index space.**
- Every identified first-order diagram except the Gaussian ρ³ triangle ρ_ij ρ_jk ρ_ki is a star or path: one centre and two legs.
- Transporting a star to D21(l+1) costs:
  - one product per distinct leg type on the a-side;
  - one product per leg type on the b-side;
  - one shared final product.
- The identified closure has four leg types, and costs 7.17 units per middle layer at Strassen level 5 (12.8 dense). A path term needs 7 n × n products with no n³ object (theory, as corrected by its verifiers); the C-leaf class shares a single arm.
- The triangle has no star form. It costs ≈ 2n⁴ FLOPs ≈ 1024 units per transition (theory: ≈ 3n units over its roles), and it buys ≤ 1 point of ε (F6.3).
- Terms that need an n³ slice explicitly cost ≈ 2n units and 8 GB each.
- *Where:* costmodel §§1, 3 (metered at n = 1024); theory §5 (exactness verified to 1e-15 against the dense n⁴ transport).
- *Implies:*
  - [C] At n = 1024, joint structure shaped like a tree in index space costs a few dense products per layer, while a cycle (triangle) costs B per layer.
  - [O] The leading-order structure is tree-like, so tree-shaped (Bethe-type) factorisations are natural at this order.

**F6.8 Pairwise gate statistics resum whole diagram classes exactly.**
- **Leaf identity.** All diagrams in which some vertex hangs on a single covariance edge sum to 6 sym3[Φ_k C_ik Γ_ij] − Wick[p(0)]. Here Γ_ij = Cov(1[z_i > 0], relu(z_j)) is the n × n gate–activation cross-covariance, and p_i(0) is the per-neuron density at the wall.
- **T-class identity.** All diagrams with one κ3_ijk hyperedge and everything else on one pair sum to κ3_ijk(Φ³ + Σ_pairs Φ_k Cov(g_i, g_j)), with the pair gate covariances.
- **Dressing.** Degree-1 vertices are dressed to all orders by the true gate probability, and degree-2 vertices by the true wall density.
- **Checks.**
  - Coverage: at network order ≤ 4 there are 187 connected diagrams, all covered; 81 of them have a C-leaf.
  - On a toy model the identities hold: the residuals scale as λ^(N+1) for N = 1–4.
  - At n = 128 the fitted LEAF coefficient is 5.8–7.1 (theory 6), and TWOLEAF is −0.93 to −1.31 (theory −1).
  - The reading that the T-class "absorbs B7's renormalisation" was refuted by the verifiers.
- *Where:* theory §1.6, §2, §3, claim C7.
- *Flags:* ≤256 for the fits; the identities are exact.
- *Implies:* [O] Pairwise gate statistics (P(z_i > 0, z_j > 0) and Cov(1[z_i > 0], a_j), both n × n) and the per-neuron wall densities sum infinite classes of diagrams exactly. They are natural carriers of non-Gaussian joint structure.

**F6.9 The (2,1) slice does not transform like a matrix under congruence.**
- Born (2,1) content of layer l − 1, transported correctly by one dense product family, matches the exact age-1 S21 at layer l with cosine 0.980–0.995.
- Transported congruently, as B ↦ WᵀBW, it reaches |cos| ≤ 0.008.
- The correct transport is a one-sided hub contraction, T′_iij = Σ_a W_ai h_a (CW)_ai (CW)_aj.
- *Where:* intel §3.1 (EscAI transport probe); 1024 × 16.
- *Implies:* [I][C] No symmetric-matrix (covariance-response) representation carries the (2,1) slice through a linear layer. Its transport costs about one dense product family per layer.

**F6.10 The third-order joint structure is invisible to sampling within the budget.**
- At n = 1024, a sample estimate of D21(l+1) from N = 32,768 samples has relative noise 0.88 (layer 1), 0.49 (layer 5), 0.28 (layer 10) and 0.20 (layer 15).
- At depth the noise is ≈ 0.2 √(32,768/N). A 1 % estimate needs N ≈ 1.4e7 at layer 14, 4.5e7 at layer 7 and 2.5e8 at layer 1.
- At n = 128 with N = 5e5, one atlas's D21 has 1–6 % noise, and the transported contribution of old sources has 2 % (age 1) to 34 % (age 15).
- *Where:* oracle1024 (one network); chain128; transfer N3.
- *Flags:* 1 net at 1024; ≤256 for the other figures.
- *Implies:* [C] The whole budget buys 65,536 forward passes (6,554 at the 0.1 floor). At n = 1024 a sampled D21 would carry ≈ 15–60 % relative noise, so third-order joint structure has to be computed, not sampled.

**F6.11 D21 is concentrated in energy, but the part that matters for accuracy is high-rank [est].**
- D21 as computed by an exact factorised K = 3 propagation at n = 1024, used as a proxy for the true slice:
  - its energy ranks are r90 = 7–20 and r99 = 150–540;
  - truncating it at every layer to rank 256, 128, 64 or 8 gives final MSE 8.2e-8, 1.45e-7, 2.5e-7 or 5.9e-7, against ≈ 4e-8 untruncated.
- The part of D21 that the current layer's slices do not produce is 18–23 % of its energy at every layer. It is high-rank: rank 32, 128 and 256 hold 24–78 %, 60–94 % and 83–98 % of it.
- The propagator leg of each source is born as the identity, so it has full rank. Compounding symmetric-Tucker truncation of the sources at rank 512, 256 or 128 costs 14×, 31× or 48×.
- *Where:* 504aldo F46, F64, F65; 1024 × 16; 8 public mini networks.
- *Flags:* est (the K = 3 proxy, not the exact slice).
- *Implies:* [I] For D21, energy rank is not accuracy rank. Truncating it even at rank n/4 costs ≈ 2×.

**F6.12 Every layer's birth matters at the readout, the first ReLU's included.**
- In a K = 3 description at n = 1024 (V16b class, raw 3.7e-8 on 8 networks), skipping the third-order birth at single layers gives:
  - layer 0 (the first ReLU, 15 layers before the readout): 1.37e-7, ≈ 4×;
  - layer 1: 1.68e-7;
  - layers {0, 1}: 3.96e-7;
  - layers {0, 1, 2}: 9.0e-7;
  - layers {13, 14}: 1.32e-7.
- At n = 128 (4 networks, κ4 teacher-forced), forcing the full κ3(a_l) of one layer to the truth multiplies the final MSE by:
  - 0.98 at l = 1, 0.95 at l = 4, 0.87 at l = 8, 0.76 at l = 12 and 0.83 at l = 14;
  - forcing only its all-distinct part: 0.95, 0.93 and 0.97 at l = 4, 8 and 12.
- No single layer dominates.
- *Where:* 504aldo F66; chain128 §5.
- *Flags:* cond; the 128 numbers are ≤256.
- *Implies:*
  - [I] Third-order structure created at the first ReLU still matters at the readout fifteen layers later, so a representation cannot forget old births.
  - [I] Removing the birth of the first ReLU, or of the last two, costs about the same (≈ 4×). Correcting the state to the truth at a single layer gains most at the late layers.

---

## 7. Fourth-order joint structure

**F7.1 The (2,1,1) slice of κ4(z) decides the third-order interface below ≈ 3 %, and its covariance-shaped regeneration captures less as width grows.**
- n = 1024, layers 7–14:
  - closure without the slice: 2.85–3.17 %;
  - with the rank-one regeneration u_i C_jk: 2.46–2.70 %, still above 2.2 %;
  - with the exact slice: 0.81–0.98 %.
- Of the gap between "no slice" and "exact slice", the regeneration closes ≈ 55 % at n = 128, 33 % at 256 and 20 % at 1024.
- At n = 64–160 the regeneration costs 2–3 points of ε at the shallow layers. Its error at layers 1–3 falls only as n^-0.46 (extrapolated 2.4 % at 1024), and refitting the coefficients around it buys nothing there.
- Step level (layers 1–6, n = 32 / 64 / 128), first-order diagram engine:
  - exact slice: 0.035 / 0.037 / 0.022;
  - regenerated: 0.074 / 0.098 / 0.068, ∝ n^-0.05;
  - zeroed: 0.113 / 0.115 / 0.085, ∝ n^-0.20.
- End to end at n = 128 (κ4 otherwise teacher-forced, geometric mean of 4 networks), final MSE:
  - exact slice 2.2e-6;
  - best rank-4 family 3.7e-6;
  - rank 1 4.8e-6;
  - u_i C_jk 6.1e-6;
  - zeroed 1.7e-5.
- *Where:* oracle1024; coef-ens; chain128 §§1, 3.
- *Flags:* 1 net and prelim at 1024; the rest ≤256.
- *Implies:* [I] The (2,1,1) slice is a genuine n³ object.
  - Its transported effect on the next (2,1) slice is ≈ 2.9 % of D21 at 1024. Assuming errors add in quadrature, a model must capture ≥ 30 % of that effect to reach 2.2 %, and ≥ 60 % to reach 1.5 %.
  - The rank-one covariance-shaped regeneration captures ≈ 15–20 %.

**F7.2 At depth the slice's transported effect is low-rank (n = 128); at shallow layers its energy is not.**
- View the slice as a family of n symmetric matrices indexed by the doubled neuron, M^(i)_jk = κ4(z_i, z_i, z_j, z_k).
  - At layers 1–3, rank 32 of 128 captures only 69–81 % of its energy.
  - At depth, 86–91 % (one network) and 56–82 % (the other) is already in the first mode.
- The effect on D21(l+1) is easier to capture at depth. At layer 10 of one network: true slice 0.031, rank 1 0.036, a C_off regeneration 0.038, no slice 0.077.
- At layers 1–5, nothing below rank 32 comes within 1.5× of the true slice.
- End to end, the effect saturates by rank ≈ 4–16. At order 2, ranks 16 and 32 give 6.5e-6 and 6.3e-6, against 5.3e-6 exact.
- *Where:* plan §3.1 (vii) (n = 128, two networks); chain128 §3.
- *Flags:* ≤256.
- *Implies:*
  - [O] At depth a few covariance-response modes may carry the slice's effect (untested at 1024).
  - [I] At shallow layers it is high-rank.

**F7.3 Pairwise fourth-order content at depth is not covariance-shaped [est].**
- Regressing the off-diagonal (2,2) slice of the post-activation κ4 on the pre-activation covariance, through the origin, gives R² = 0.998 at layer 1, falling monotonically to 0.045 at layer 15. The coefficient λ falls from 0.32 to 1.6e-3 (AndreasHad04; 8 networks; that team's chain state).
- A different object follows the covariance closely: the rank-one matrix core of the augmented K = 3 chain, built from the (2,1,1) hub contraction and the symmetrised K31, follows λ C_pre off the diagonal. Its R² is 0.48–0.52 at layer 0, 0.83–0.86 at 1, 0.92–0.93 at 2, ≥ 0.95 from 3 on and 0.985–0.99 at 13–14 (8 networks; 504aldo F105).
- The physical C-aligned fourth-order content has λ ≈ 2.1–2.6e-3, flat in depth. The chain's fitted 7–8e-3 is ≈ 3× that, and its gain came from cancelling the chain's own truncation bias (EscAI).
- *Where:* intel §§2.4, 3.3; 504aldo F105; 1024 × 16.
- *Flags:* est.
- *Implies:*
  - [I] At depth the post-activation (2,2) slice carries real content that is not shaped like the covariance.
  - [I] A covariance-shaped account of fourth-order content is physical only in its small C-aligned part.

**F7.4 Index-coincident fourth-order content carries more of the remaining error at larger width [≤256].**
- Keeping all index-coincident κ4 entries ((4), (3,1), (2,2), (2,1,1)) in a compressed fourth-order state cuts raw MSE by 35 % at width 128 (3 networks) and by 64 % at width 192 (1 network).
- The diagonal alone cuts 24.5 %. The further index-coincident entries therefore add ≈ 14 % relative; the digest attributes this increment to the (2,1,1) slice.
- *Where:* intel §§3.4, 6.4 (EscAI `observable_push`, their reference algorithm).
- *Flags:* ≤256, est.
- *Implies:* [I] The share of the error that sits in fourth-order joint content grows with width, at least from 128 to 192, consistent with F7.1.

**F7.5 At depth the per-neuron fourth cumulant is made mostly from the off-diagonal fourth-order structure of the previous layer (n = 128).**
- A closure that transports only the (4), (3,1) and (2,2) slices of κ4(a) gets K4 = κ4(z)_aaaa wrong by:
  - 16 % at layer 1, 47 % at 4, 65 % at 8;
  - 78–91 % at layers 13–15, which is worse than setting κ4 = 0 at depth.
- Carrying the full κ4 with every same-order diagram leaves 9 %, 12 %, 34 % and 33–46 % at the same layers.
- Part 1 §3.3 (n = 24–32): at a15 the (2,2), (2,1,1) and (1,1,1,1) patterns of κ4(a) each exceed the diagonal pattern by factors 3–10 in their contribution to λ4(z).
- *Where:* chain128 §6 (one network, order-2 Edgeworth step); Part 1 §3.3.
- *Flags:* ≤256; est for the closure numbers.
- *Implies:*
  - [I] The per-neuron κ4 cannot be propagated neuron by neuron: at depth it is mostly the contraction of off-diagonal fourth-order structure.
  - [I] At n = 128, even all first-order diagrams leave a third of K4 unexplained at depth. Fourth-order generation is less well described at first order than third-order generation. This is untested at 1024.

**F7.6 The neuron-averaged fourth cumulant is one coherent scalar per layer.**
- E_w κ4(z_{l+1,i}) = 3σ⁴[Var‖a_l − μ_l‖² − 2‖Σ_l‖_F²]. This is the thin-shell excess of the centred activation norm.
- It is exact, and Part 1 checked it to 0.1–2.2 % at five layers of one n = 1024 network.
- Outside the assigned sources: the bethe design stream's width-64 oracle reports a rank-one spike in κ(z_a, z_a, z_b, z_b), seeded by the input radius and amplified with depth ([CONVERGENCE.md](CONVERGENCE.md)).
- *Where:* Part 1 §3.3.
- *Flags:* 1 net for the check.
- *Implies:* [O] The coherent part of the fourth-order structure is one scalar per layer: the law of the activation norm (F2.4). It is cheap to carry.

---

## 8. Memory across depth (old content)

**F8.1 At n = 1024, about 40 % of the next (2,1) slice is old content, nearly orthogonal to everything built from the current layer.**
- Measure: the pass-through Φ³κ3(z_l) of the current pre-activation's all-distinct third cumulant (B0 of F6.3). It is the all-distinct part of the third-order structure born at earlier ReLUs; the index-coincident part reaches D21(l+1) through other diagrams (B1, B2).
- At n = 1024 it is 0.38–0.42 of ‖D21(l+1)‖ at every layer 1–14 (noise-free replica cross-product).
  - Dropping it leaves ε = 38–42 %.
  - Its held-out projection on all of the current layer's birth terms (slices, Wick, Hermite terms, B1–B6) has R² = 0.01–0.14.
- At n = 128 and 256 the share rises with depth:
  - n = 128: from 0.36–0.41 at layer 1 to 0.58–0.59 at layer 14, with projection R² 0.04–0.69;
  - n = 256: from 0.38 at layer 1 to 0.51–0.53 at layer 14, with R² 0.02–0.46.
- *Where:* oracle1024 `results/summary_prelim.md`. At 1024: one network, N = 32,768 × 2 replicas. At 128 and 256: two networks each, N = 5e5 and 1e6.
- *Flags:* 1 net and prelim at 1024.
- *Implies:* [I] About two fifths of the third-order interface is memory. No function of the current layer's state reconstructs it, so it must be carried.

**F8.2 Content of every age up to about ten layers carries several per cent of D21.**
- Windowed shares (old-content stream, sources with their slices). Share of ‖D21‖ at layers 10–15 held by content older than w layers:

  | window | n = 64 | n = 128 (two networks) | n = 256 |
  |---|---|---|---|
  | w = 1 | 0.71–0.97 | 0.88–0.99 | 0.91–0.98 |
  | w = 2 | 0.52–0.74 | 0.65–0.95 | 0.76–0.93 |
  | w = 4 | 0.26–0.48 | 0.50–0.88 | 0.49–0.77 |

  The share grows from 64 to 128 and is flat from 128 to 256.
- All-distinct convention at n = 128: the age-1 source is 0.75–0.93 of ‖D21‖, ages ≥ 2 together 0.37–0.63, and ages ≥ 5 together 0.06–0.46.
- Per-age shares from the transfer digest (n = 128, one network):

  | age | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 | 15 |
  |---|---|---|---|---|---|---|---|---|---|---|
  | all-distinct source | 0.56 | 0.48 | 0.42 | 0.36 | 0.28 | 0.23 | 0.12 | 0.05 | 0.026 | 0.003 |
  | full source (with slices) | – | 0.91 | – | 0.72 | – | 0.57 | 0.47 | 0.33 | 0.16 | 0.04 |

- *Where:* old-content (atlases A1, A2 and B1 at n = 128, N = 6e5; K64 and K256); transfer §§3.2, 3.6. The three decompositions use different conventions for births and slices, so their numbers are not comparable across the tables.
- *Flags:* ≤256.
- *Implies:* [I] A window shorter than ≈ 8–10 layers drops several per cent of D21 against a 2.2 % bar. At 1024 the share by age is not measured (Section 16).

**F8.3 Old content is orthogonal to young content and is not absorbed into renormalised births.**
- Orthogonality:
  - the D21 contributions of old sources are nearly orthogonal to each other and to the young ones;
  - linearly transported content and the rest of D21 are nearly orthogonal (transfer);
  - even at age 2, linear transport leaves 84 % of D21(z_t) unexplained.
- Absorption, by per-layer least squares of the old pool's D21 (window 4, layers 7–15):
  - on the memoryless birth diagrams, it removes at most ≈ 10–15 % of it, even in-sample;
  - adding the young sources removes ≈ 30–40 % in-sample;
  - on a held-out network the residual stays at 0.35–0.79, against a share of 0.39–0.88, i.e. at the drop-it level.
- The residual is 10–35× the 2.2 % target at every width from 64 to 256, and it does not shrink with n.
- *Where:* old-content (n = 64, 128, 256; an independent-sample replicate reproduces the numbers to ±0.004); transfer §3.2.
- *Flags:* ≤256.
- *Implies:* [I] Old content is not a function of the current layer's state, nor of recent births. Any estimator must carry it explicitly, in some form, across depth.

**F8.4 In a K = 3 description at n = 1024, cutting memory short costs 5× to 60×.**
- Keeping only the last W third-order sources (4 networks), raw MSE against 3.88e-8 with all of them:

  | W | 12 | 8 | 6 | 4 | 3 | 2 |
  |---|---|---|---|---|---|---|
  | raw MSE | 4.64e-8 | 1.86e-7 | 4.35e-7 | 1.02e-6 | 1.54e-6 | 2.24e-6 |

- In the regenerated chain, keeping 4–5 young sources gives 1.35–2.9e-7, against 2.13e-8 with the old tier.
- *Where:* 504aldo F61 and §4; 1024 × 16; public mini networks.
- *Flags:* cond (within a K = 3 description).
- *Implies:* [I] The final means depend on third-order content born up to ≈ 12 layers back. A window of 8 costs 5×.

**F8.5 Old third-order content is purely a property of the fixed network.**
- The ensemble-averaged transfer of third-order content is identically zero: E_W[W^{⊗3}] = 0 (E3). Old content survives only because the weights are quenched.
- *Where:* transfer §2.5 (exact).
- *Implies:* [I] No annealed or self-averaging shortcut exists for old content. It must be computed from the actual weights, network by network.

---

## 9. Gated propagators

**F9.1 The participation ratio of gated products follows the free-probability law PR ≈ n/(2·age).**
- The free-probability prediction is PR(J_{s→t}) ≈ n/(1 + Σ_l (r(A_lᵀA_l) − 1) + Σ_l (r(D_l²) − 1)).
- Measured against it (n = 128, one network):

  | age | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 15 |
  |---|---|---|---|---|---|---|---|---|---|
  | PR measured | 63.8 | 32.2 | 21.2 | 15.6 | 9.5 | 6.2 | 4.4 | 3.8 | 3.4 |
  | free law | 63.8 | 32.5 | 21.6 | 16.2 | 10.7 | 8.0 | 6.4 | 5.4 | 4.4 |

  The agreement is within 4 % up to age 4. Beyond that the measured PR falls below the law, by up to 31 % at ages 9–13; the source attributes this to the mean outlier (F2.2), as an interpretation.
- PR·age/n = 0.47–0.51 for ages ≤ 7 at every width from 128 to 1024 (fresh networks).
- The masked propagator at age 5 has PR 6.0, 10.5 and 23.0 at n = 64, 128 and 256 (old-content).
- Inputs to the law: r(AᵀA) = 1.99–2.05 per layer, and r(D_l²) = 1.0 at layer 0 rising to 1.5–2.35.
- *Where:* transfer §§2.2, 3.1, 3.4; old-content.
- *Flags:* ≤256 for the atlas networks; ens and Monte Carlo gates at 512–1024.
- *Implies:* [C][I] The effective rank of transported content is a fixed fraction of n, decaying only as 1/age.

**F9.2 The propagators have no spectral gap; the competition shape is deep in the free regime.**
- Finite-time Lyapunov exponents of J_{0→15} at n = 128:
  - top eight: +0.059, +0.030, +0.007, −0.012, −0.026, −0.037, −0.043, −0.057;
  - the 16th, 32nd and 64th: −0.149, −0.330, −0.966.
- The mean consecutive gap is 0.013 per layer for i ≤ 16 and 0.018 for i = 17–64. That is 3–5 × 1/(2n), and nowhere of order 1 except at the top (the mean direction).
- Hanin–Nica's log-normality parameter for these products is β ≈ 5 Σ 1/n_i: ≈ 0.59 at n = 128 and ≈ 0.07 at n = 1024. At n = 1024 a propagated squared norm therefore fluctuates by a factor ≈ e^(±0.26) (theorem applied, not measured).
- *Where:* transfer §§2.3, 3.1.
- *Flags:* ≤256 for the exponents.
- *Implies:* [I] There is no exponential forgetting of old content in the bulk: concentration is polynomial in age. Only generic directions are discounted, by the per-layer gain 2E[Φ²] ≈ 0.58–0.95 (F4.4).

**F9.3 The number of modes old content needs is a fixed fraction of n, the same at every width from 128 to 1024 [ens].**
- k_2 % is the number of propagator modes needed for 2 % error on a source's own D21. Its fraction k_2 %/n (fresh He networks, depth 16, gates from 2e4 Monte Carlo inputs):

  | age | n = 128 | n = 256 | n = 512 | n = 1024 |
  |---|---|---|---|---|
  | 1 | 0.947 | 0.946 | 0.946 | 0.946 |
  | 2 | 0.707 | 0.708 | 0.704 | 0.709 |
  | 4 | 0.488 | 0.503 | 0.495 | 0.499 |
  | 6 | 0.370 | 0.387 | 0.381 | 0.386 |
  | 8 | 0.299 | 0.316 | 0.310 | 0.315 |
  | 10 | 0.247 | 0.264 | 0.261 | 0.265 |
  | 12 | 0.213 | 0.229 | 0.227 | 0.229 |
  | 15 | 0.172 | 0.191 | 0.188 | 0.192 |

- Absolute numbers at n = 1024, for 20 %, 10 %, 5 %, 2 % and 1 % error: age 8 needs 162, 214, 263, 322 and 363 modes; age 15 needs 91, 124, 156, 197 and 227.
- The ensemble formula reproduces the exact-atlas k_2 % at 12 of 15 ages on one n = 128 network and 11 of 15 on the other. Its per-k curves deviate by up to 41 % at ages 6–14, and by up to 1.9× at age 15.
- For the merged old tier (targets t ≥ 10, 2 % of the whole D21) at n = 128:

  | tier | propagator basis | HOSVD oracle |
  |---|---|---|
  | born ≥ 4 layers back | 48–56 modes | 24–32 |
  | born ≥ 6 layers back | 28–48 | 12–24 |
  | born ≥ 8 layers back | 12–28 | 4–16 |

  Scaled to n = 1024 the propagator figures become ≈ 384–448, 224–384 and 96–224.
- A Tucker core of rank k costs k³/n² units per layer.
- *Where:* transfer §§3.2, 3.4, 4.
- *Flags:* ens at 256–1024; the n = 1024 tier sizes assume the old-content share is independent of width, which is not measured there.
- *Implies:* [C] Carrying old third-order content in a basis derived from the propagators needs O(n) modes at 1024, ≈ 200–500. Only content older than ≈ 10 layers fits in ≈ 0.1 n modes, i.e. under 1 unit per layer.

**F9.4 Better bases save constant factors, not the linear growth in n.**
- At n = 128:
  - the oracle HOSVD basis of the content itself needs 1.2–1.7× fewer modes than the propagator basis;
  - the transported covariance J C(a_s) Jᵀ closes a third to a half of that gap for all-distinct sources at ages 2–8, and comes within 2–8 modes of the oracle for full sources;
  - the current covariance C(z_t) is a poor basis, because it is dominated by young content;
  - a basis fixed at birth costs 7–25 % more modes, because the forward filtration stabilises only over 5–10 layers.
- Old sources share their target-side subspace: alignment 0.34–0.91 at k = 8, against 0.06 for random subspaces.
- *Where:* transfer §§3.1, 3.2, 3.6.
- *Flags:* ≤256.
- *Implies:* [O] Structure the computation already holds (transported covariances) is a nearly optimal basis for covariance-type content. [C] The mode count stays proportional to n.

**F9.5 At n = 1024 the propagator legs need 0.25–0.4 n directions [est].**
- Propagator legs of age ≥ 5 have ≥ 92 % of their energy in rank 128 and ≥ 99 % in rank 256. Age-1 legs have 53 % and 81 %.
- A shared basis for all sources of age ≥ a, with r columns, raw MSE relative to the untruncated control:
  - a = 4: +1.0 % at r = 384, +6.7 % at 320, +22 % at 256;
  - a = 3: +96 % at r = 256;
  - a = 6: +0.8 % at r = 256;
  - nested at age ≥ 7: fine at r = 224, a cliff of +15 % at r = 128.
- *Where:* 504aldo F66, F72, F73; 1024 × 16; public mini networks.
- *Flags:* est (in a K = 3 chain).
- *Implies:* [C] Consistent with F9.3: at 1024, carriers of old content need ≈ 0.2–0.4 n shared directions, and the accuracy falls off a cliff below that.

**F9.6 One layer scrambles the orientation of a source, so a propagator's spectrum predicts the error of any projection built from it.**
- J_{s→t} = R·A_{s+1}, where A_{s+1} is Gaussian, independent of the source, and right-orthogonally invariant. The source's orientation relative to J is therefore Haar-random, up to the dependence of the downstream gates on A_{s+1}.
- The expected projection error then depends only on the spectrum of J. It is source-independent for all-distinct (harmonic) sources.
- Checks: Monte Carlo at n = 32 agrees to < 1 % in ε, and the formula reproduces the atlas k_2 % at 11–12 of 15 ages.
- *Where:* transfer §§2.4, 3.2.
- *Flags:* ≤256.
- *Implies:* [O] To that accuracy, designing a carrier for old content reduces to spectral data of the gated propagators, which is cheap to compute.

---

## 10. How errors reach the output

**F10.1 Final MSE grows as ε² in the relative error of the (2,1) slice; the constant equals the Gaussian-closure MSE at n = 1024.**
- For a relative rms error ε of D21 at every layer, extra MSE ≈ 4.2e-6·ε² (n = 1024, public mini networks; measured independently by EscAI).
  - ε = 1 recovers the Gaussian closure (4.4e-6).
  - A margin of 10 % over raw 2.1e-8 needs ε ≤ 2.2 %, i.e. 99.95 % of D21's energy reproduced.
- At n = 128 (4 networks, κ4 teacher-forced):
  - over structured (closure) errors, final MSE = a + kε² with k = 0.7–2.0 × that network's Gaussian-closure MSE (correlation 0.86–0.995);
  - injected independent Gaussian noise costs ΔMSE/ε² = 2.1–2.3e-3 (order-2 step, quadratic to ±10 % over ε = 0.02–0.2), i.e. ≈ 4× more per ε² than structured errors.
- The interface ladder at n = 128 (geometric means of final MSE):

  | κ3 description | final MSE | layer-rms ε(D21) |
  |---|---|---|
  | Gaussian closure | 2.13e-4 | 1 |
  | slices only | 9.8e-5 | 0.53 |
  | Wick term | 2.6e-5 | 0.16 |
  | counting closure | 4.3e-6 | 0.10 |
  | all first-order diagrams | 2.2e-6 | 0.036 |

- *Where:* 504aldo F71(f); intel §3.4; chain128 §§2, 4.
- *Flags:* cond (measured inside chains); the injected-noise law is ≤256.
- *Implies:* [I] Whatever represents pairwise third-order content must reproduce D21 to ≈ 2–5 % relative at every layer to reach raw ≈ 1e-8. Errors that are random with respect to the true structure are about four times more harmful than errors that follow it.

**F10.2 Per sample, the final activation is a nearly linear function of the penultimate layer and only weakly of early layers.**
- At n = 1024 (3 networks, N = 32,768), the best linear predictor of the final activations from layer l explains:

  | layer l | h_2 | h_4 | h_8 | h_12 | h_14 | h_15 |
  |---|---|---|---|---|---|---|
  | R² | 0.30 | 0.42 | 0.64 | 0.84 | 0.92 | 0.96 |

- P1 (256 × 32, 64 networks): anchoring a sample ensemble to the exact means at layer L reduces the variance of the final estimate by:

  | L | 1 | 2 | 4 | 8 | 16 | 24 | 30 | 31 |
  |---|---|---|---|---|---|---|---|---|
  | factor | 1.04 | 1.15 | 1.33 | 1.79 | 3.00 | 6.31 | 25.6 | 52.1 |

- At P1, given the exact mean and full covariance of layer 31, only 1/310–1/394 of the per-sample variance remains (16 networks). The mean alone leaves 1/61, and the mean with the diagonal 1/118–1/124.
- *Where:* 504aldo F62; Oishi oracle-anchoring probe.
- *Flags:* the R² row is at 1024 (3 networks); the anchoring rows are P1.
- *Implies:*
  - [I] The per-sample fluctuation of the output is made in the last few layers.
  - [I] At P1, the mean and covariance of the penultimate layer determine all but ≈ 0.3 % of it. The remainder is higher-order structure, which is where the target lives.

**F10.3 Mean errors are damped with depth in generic directions, but not near the output and not along the mean.**
- Measured with exact reverse-mode averages of the propagator (Part 1 §3.5). The cumulative power gain of a mean error from a_l to the output, ‖E J_{l→L}‖_F²/n, at n = 64 / 128:

  | from | input | a1 | a4 | a7 | a10 | a15 |
  |---|---|---|---|---|---|---|
  | gain | 0.034 / 0.0092 | 0.074 / 0.018 | 0.21 / 0.050 | 0.37 / 0.12 | 0.50 / 0.24 | 0.89 / 0.84 |

- Summed over source layers, independent mean errors at every layer reach the output amplified 2.2–2.7×.
- At P1 a relative error ε in the layer-31 means gives final MSE 0.945 ε², so ε ≤ 1.03e-4 is needed for 1e-8 (Oishi).
- Errors along the mean direction are amplified, not damped (F2.2). Errors in the scale (dilation) mode pass with gain exactly 1 (Part 1, E1).
- *Where:* Part 1 §3.5 (n = 64, 128); Oishi.
- *Flags:* ≤256, P1.
- *Implies:* [I] Every layer's means must be right to ≈ 4e-5 rms at n = 1024 (Part 1 §3.6), and coherent or mean-aligned errors are not forgiven.

**F10.4 The expansion is asymptotic at small width, and consistent truncation matters more than order.**
- At n = 128, adding the κ3²/2 Edgeworth term without the matching κ5, κ6 and κ3·κ4 terms makes the final means 1.5–3.4× worse (2.2e-6 → 5.3e-6 with κ4 teacher-forced). It also destabilises the deep layers of a self-consistent chain (one network: 8.3e-5 at order 2 against 6.9e-6 at order 1).
- Edgeworth corrections are not positivity-preserving: at n = 64, order 2, a chain reached negative variances.
- The second-order diagram series overcorrects at n = 128 (F6.4). At the Phase 1 shape the cumulant expansion is asymptotic, with optimal order K = 3; K = 4 is 10–20× worse even if free (504aldo F27–F34).
- At n = 1024, consistent second-order corrections to the marginal moments change raw MSE by < 1 % (3 networks): κ3² added gives ratio 1.0011, κ3·κ4 added 0.9969 (EscAI).
- *Where:* chain128 TL;DR item 4 and §2; theory §3; 504aldo §3.4; EscAI `tenfold_push` README.
- *Flags:* mostly ≤256 and P1.
- *Implies:* [I] At the competition width first order in 1/n is adequate. At narrow widths, and at the smoke-test shape, an inconsistently truncated expansion is worse than a lower-order consistent one.

**F10.5 What a K = 3 description misses is not a function of cheap local statistics.**
- 504aldo:
  - the final-layer residual of a K = 3 chain has |corr| < 0.09 with every local feature tried (F46);
  - an online per-layer mean correction with 13 features gains 1.45× on one base and 1.0× on the regenerated base: "K3 content is fundamentally NONLOCAL" (F47, F68);
  - slice fits add 0.02 to R² and nothing to the MSE (F65);
  - a ridge with 144 constants gains 1.07× (F82).
- EscAI: a 40-feature span explains 8.77 % of the baseline error, and its moment fit gains 3.2 % on held-out networks (n = 1024). At widths 128–192, weight-spectrum predictors of per-network difficulty flip their rank correlation from +0.71 to −0.43 between two tests.
- *Where:* 504aldo F46, F47, F65, F68, F82; intel §3.2; EscAI `spectral_push`.
- *Flags:* cond.
- *Implies:* [I] The remaining information is structural joint content. Calibration cannot stand in for it, and no fitted correction on per-layer summaries recovers it.

---

## 11. Variation across networks

**F11.1 The scored means are quenched.**
- Per-neuron means fluctuate around their ensemble values by O(n^-1/2), far above the ≈ 1.2e-4 rms accuracy the bar needs (BRIEF §3).
- Old third-order content has ensemble mean zero (F8.5).
- *Implies:* [I] Per-network structure must be resolved from the actual weights.

**F11.2 The structural constants are properties of the shape (n, L), not of the network.**
- Renormalised closure coefficients: a held-out ensemble table equals the per-network fit to ≤ 0.002 in ε at every layer and width from 64 to 160, for the mean and the worst network alike (F6.4).
- The three coefficients that matter (B2, B3, B6) vary across networks with sd 0.00–0.14, against a per-fit Monte Carlo sd of 0.00–0.02, so that spread is real. The other four vary more (sd 0.3–1.0 at depth), but in flat directions of the error. The held-out table is still as good as the per-network fit.
- 504aldo's fourth-order λ table varies 2–4 % across networks at layers 1–10 and 6–10 % at layers 11–14. Transferring a frozen table from the other networks costs 1–4 % (F68, §3.3).
- The free-probability laws of the propagators (F9.1, F9.3) and the mean-field spread of t (F3.6) are network-independent to leading order.
- *Where:* coef-ens; 504aldo; transfer.
- *Flags:* ≤256 for the coefficient tables; the λ table is est (1024).
- *Implies:* [O] Constants that describe the structure (diagram coefficients, scaling exponents, mode fractions) can be computed offline per shape and shipped as data, gated on (n, L). Per-network information enters only through the weights.

**F11.3 At n = 1024, network-to-network spread is ±20–40 % for a second-order description and ±12–20 % for a K = 3 description.**
- Gaussian-closure distance:
  - 2.5–5.9e-6 over 8 public mini networks (504aldo F43);
  - 4.30e-6 ± 0.36e-6 (standard error) over 6 bench networks, an sd of ≈ 20 %.
- K = 3 chain per network: 1.88–2.39e-8 (V19, 8 networks) and 1.96–2.61e-8 (V25).
- EscAI's baseline: ≈ 2.0–2.8e-8 over 8 networks.
- Plain sampling: 1.086e-5 ± 0.085e-5 (standard error, 6 networks).
- 504aldo's 8 local networks run 6.5 % above the 100-network public board, at every version.
- *Where:* 504aldo §3.1, F43; bench; EscAI `state_attribution`.
- *Implies:* [I] For measurement, a difference below ≈ 10–15 % is not resolved by 6–8 networks unless the comparison is paired on the same networks.

---

## 12. What sampling sees

**F12.1 The per-sample variance fixes the price of plain sampling.**
- The per-neuron final variance is 0.071–0.075 (F2.1). Plain sampling therefore has MSE ≈ 0.0748/N.
- At the 0.1 B floor (6,554 samples) raw MSE is 1.09e-5 (bench, 6 networks), or 1.2–1.5e-5 on three networks (504aldo F62).
- Raw 1e-8 needs N ≈ 7.5e6 samples ≈ 114 B (Part 1 §4).
- At the Phase 1 shape the constant is C = 5.12e-2 (Oishi). The BRIEF's "≈ 0.05 per sample" is that P1 value.
- *Where:* bench; 504aldo F42, F62; Oishi.
- *Implies:* [C] Sampling alone is ≈ 1,100× short at the floor. Any sampling component must come with a surrogate that removes ≥ 99.9 % of the per-sample variance.

**F12.2 Variance reductions that can actually be realised are small at n = 1024.**
- At n = 1024 (3 networks, N = 32,768), variance-reduction factors:

  | device | factor |
  |---|---|
  | control variate from x | 1.18–1.20× |
  | control variate from h_1 | 1.31–1.35× |
  | antithetic | 1.10× |
  | Sobol | 1.37–1.48× |
  | Sobol with a control variate | 1.2–1.3× |
  | oracle linear control variate from h_15 | only 25–27× |

- Tying an adjusted 3.67e-8 (a K = 3 description) needs 33–40×; tying an adjusted 4.9e-9 needs 245–300×.
- *Where:* 504aldo F62.
- *Implies:* [C] No realisable input-side device closes the gap. The per-sample fluctuation is generated late (F10.2), where no exactly integrable surrogate exists.

**F12.3 The per-sample output has heavy high-degree chaos content [P1].**
- Hermite-chaos spectrum, as fractions of the per-sample variance (120 networks):
  - degree 1: 25.7 %; degree 2: 22.9 %; degree 3: 9.5 %; degree 4: 4.4 %; degree 5: 2.9 %;
  - degree 6: 6.4 %; degree 7: 5.0 %; degrees 8–15: 0–1.4 % each;
  - degrees 16 and 17: 5.4 % and 4.3 %; degrees above 60: 4.8 %;
  - odd degrees together 54.3 %, even 45.7 %.
- Variance left after a rule is exact to degree D:

  | D | 1 | 2 | 3 | 5 | 7 | 9 | 15 | 32 | 64 | 128 | 256 |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | variance left | 74 % | 51 % | 42 % | 33 % | 27 % | 21 % | 20 % | 9.4 % | 4.1 % | 2.2 % | 0.68 % |

  So exactness to degree 5 buys 3.0× and to degree 32 buys 10.7×.
- The degree-1 chaos has 79 % of its energy in one direction (the averaged Jacobian, F2.4). The residual of degree ≥ 4 is isotropic, with effective dimension ≳ 224 of 256.
- The mean chaos order is ≈ 25 (504aldo F21).
- *Where:* Oishi chaos-spectrum and active-subspace probes; 256 × 32.
- *Flags:* P1.
- *Implies:* [C] No low-degree cubature or quasi-Monte Carlo rule is competitive: the kinks spread the variance over high degrees and over all directions.

**F12.4 Per face, the linear map is extremely low-rank in energy but needs about 32 directions for accuracy [P1].**
- The per-sample input–output Jacobian has participation rank 3.95 ± 0.82, and rank 5.6, 11.9 and 18.4 for 90, 99 and 99.9 % of its energy; σ₆₄/σ₁ = 2e-6.
- Truncating it to rank r gives a bias² (against a baseline MSE of 5.73e-6) of:

  | r | 4 | 8 | 18 | 32 | 64 |
  |---|---|---|---|---|---|
  | bias² | 4.6e-2 | 3.2e-3 | 9.6e-6 | 1.7e-8 | 9e-15 |

- Two random hubs share their output-side top-8 subspace (cos² 0.327, 10× chance) far more than their input side (0.076, 2.4× chance).
- *Where:* Oishi `jacobian-reuse`; 8 networks × 384 exact Jacobians; 256 × 32.
- *Flags:* P1.
- *Implies:*
  - [I] Energy rank and accuracy rank differ by ≈ 8×.
  - [O] The shared structure across faces sits on the output side, consistent with the shared target-side subspace of old content (F9.4) and the output mean direction (F2.2).

---

## 13. Scaling with width and depth, collected

| quantity | measured law (values) | widths (networks) | source | reading |
|---|---|---|---|---|
| Gaussian-closure distance of the final means (raw MSE) | n^-0.82 from 64 to 128, then ≈ n^-2.0: 5.1e-4, 2.9e-4, 4.3e-6 | 64 (8), 128 (8), 1024 (6) | bench | regime change between 128 and 256; small-width fits overstate the 1024 value 12× |
| first-order closure error ε(D21), deep layers | ≈ n^-0.8: 4.7 → 2.6 → 0.9 % | 128 (2), 256 (2), 1024 (1) | oracle1024 | the generation law sharpens quickly with width |
| Wick term alone, ε(D21) | ≈ n^-0.25: 7–10 → 4–4.5 % | 128–1024 | oracle1024 | Gaussian structure is not the leading correction |
| index-coincident slices only, ε(D21) | n^-0.10: 0.56 → 0.49 | 32–128 (1) | chain128 | the all-distinct part does not shrink relative to D21 |
| all first-order diagrams with the exact (2,1,1) slice, ε | n^-0.3 to n^-0.6 | 32–128 (1) | chain128 | |
| same, with the slice regenerated or dropped | n^-0.05 to n^-0.2 | 32–128 (1) | chain128 | the (2,1,1) term becomes relatively more important with width |
| share of the κ4 gap that the rank-one regeneration closes | 55 % → 33 % → 20 % | 128, 256, 1024 | oracle1024 | gets worse with width |
| renormalisation drift of the coefficients | n^-0.4 to n^-1.0 (B3: n^-0.1) | 64–160 (3–8) | coef-ens | absent at 1024 |
| old pool share of D21 (window 4, layers 10–15) | 0.26–0.48 → 0.50–0.88 → 0.49–0.77 | 64, 128, 256 | old-content | flat from 128 to 256 |
| pass-through (B0) share of D21 | 0.36–0.59 (128), 0.38–0.53 (256), 0.38–0.42 (1024) | 2, 2, 1 | oracle1024 | about 40 % at 1024, flat in depth there |
| modes for a 2 % old-content carrier, k/n | constant (e.g. age 8: 0.30–0.32) | 128–1024 (ens) | transfer | mode count ∝ n |
| PR of the masked propagator at age 5 | 6.0, 10.5, 23.0 | 64, 128, 256 | old-content | ∝ n |
| mean-mode outlier at age 15: s₁²/mean s², squared overlap with μ | 46 → 128; 0.49 → 0.85 | 128 → 1024 (fresh networks) | transfer | sharpens with width |
| face-law η₀ at depth | × 1.2–2.1 from 128 to 1024; levels off at layers 9–12 | 128–1024 | transfer | slow growth, not shown bounded |
| joint slice sizes at layer 1, n × 4 | ρ ÷ 1.94; κ3 slices ÷ 4.1–4.2; K22, K4 ÷ 4.6–5.0; K211, K31 ÷ 7.8–8.3 | 32 vs 128 | theory | as counted (E6) |
| repeated-index κ4 content's share of the remaining error | −35 % (128) → −64 % (192) | 128 (3), 192 (1) | intel | grows with width |
| marginal non-Gaussianity | rms λ3 ≈ (5–7) l/n, λ4 ≈ (3.6–6) l/n | 128–1024 | Part 1 | linear in depth |
| spread of t = μ/σ at layer 16 | 2.05 (n = 128), 3.11 (n = 256, layer 16 of a 32-layer network), 3.70 (n = 1024); infinite width 3.46 | | Part 1, Oishi | mean field reached by 1024 |
| cumulant expansion of the propagation | asymptotic, best K = 3 at L/n = 1/8 | P1 | 504aldo | |
| readout Edgeworth series | convergent, ×13–16 per CLT weight from weight 2 on, at L/n = 1/64 | 1024 (8) | intel | |

**Depth.** Five things grow with depth:
- the non-Gaussianity of each marginal (F3.4);
- the spread of t, and with it the frozen-gate fraction (F3.6, F4.1);
- the correlations between neurons (F5.1);
- the mean-direction outlier (F2.2);
- the share of old content in the pairwise third-order structure at n ≤ 256 (F8.1).

At n = 1024 that share is flat at ≈ 40 %.

The participation ratio of the propagators falls as 1/age. The damping of mean errors weakens towards the output (F10.3).

---

## 14. The meter and the grader (facts about cost, not about the object)

These are not properties of the networks, but they convert every information requirement above into a price. Sources: costmodel (metered with flopscope 0.12.1 at n = 1024), submission, intel §1, 504aldo §§2.2, 6, Part 1 §4.

**M1. Units.**
- 1 unit = one dense 1024³ product = 2^31 FLOPs, and B = 1024 units.
- One forward pass costs 3.357e7 FLOPs (≈ 2^25), so B buys ≈ 65,500 passes and the 0.1 floor ≈ 6,550.
- Implies: [C] at the floor, a layer has 6.4 units ≈ 6 dense products.

**M2. Products.**
- A dense f32 product costs 0.9995 u. Any f64 operand doubles the bill (1.999 u).
- Batching changes the call count, not the price.
- A same-object Gram costs 0.5002 u. A product (n,n)@(n,r) costs r/n units.
- Implies: [C] carrying one n × n object costs 1–2 units per layer.

**M3. Strassen–Winograd** written as flopscope ops is permitted, with no recursion limit.
- Price per product at levels 1–5: 0.877, 0.772, 0.683, 0.609 and 0.5555 u.
- A family costs ≈ 113 calls at level 5 (67 at level 3), whatever its batch size.
- The f32 relative error is 3.6e-6 at level 5.
- Implies: [C] ≈ 11 products per layer at the floor, paid for in calls.

**M4. Symmetric sandwich** WᵀΦSΦW.
- Its cheapest form costs 1.027 u (one Strassen product plus three half blocks, 214 calls).
- The tagged three-operand einsum costs 1.504 u in only 3 calls. But it raises SymmetryError on exactly symmetric, indefinite matrices of O(1) scale (the float32 allclose check), and the weighted aliased Gram einsum fails the same way. Either one can zero a network.
- The safe weighted Gram (split by sign, aliased) costs 0.503 u.
- Implies: [C] symmetric tricks are safe only on covariance-like (SPD, small-scale) inputs.

**M5. Elementwise operations.**
- Ordinary elementwise ops cost 1 FLOP per element (0.0005 u per n × n pass); exp, log and x² cost 16.
- norm.cdf and norm.pdf are billed at float64: 96 and 54 per element of an n × n array (0.047 and 0.026 u).
- Per-layer Gaussian weights on n-vectors cost 0.0001 u.
- Implies: [O] pairwise nonlinear statistics, such as bivariate gate probabilities and wall densities, are cheap.

**M6. Residual (Python-side) time.**
- It costs ≈ 0.022 ms per flopscope call on the grader (0.042 ms on our box). The 0.4 s cap therefore allows ≈ 18,000 calls per network, ≈ 1,100 per layer.
- Pooled `out=` buffers cost 0.013–0.03 ms per call; fresh n × n results 0.09–0.1 ms.
- Residual timing is strongly machine-dependent: public V29 measured 0.25–0.55 s on different machines.
- Implies: [C] calls, not FLOPs, bind a design that uses many small operations or deep Strassen.

**M7. Grader resources.**
- Throughput ≈ 1e10 FLOP/s: a 0.34 B bill takes 76–100 s per network, against the 120 s cap.
- The solution gets 2 vCPU and the backend 14. The backend caps each array at 4 GiB (a dense f32 n³ tensor is 4 GiB).
- Memory is 8 GB. setup() has 5 s, runs 5–15 times per submission, and a failure there fails the whole submission.
- Implies: [C] no explicit n³ objects.

**M8. Robustness.**
- The smoke test runs an MLP deeper than 16 layers, probably 256 × 32 (the Phase 1 shape). Any table indexed by layer that assumes L = 16 fails the whole submission.
- Assigning to `.shape` fails on the grader client.
- Any failure on a scored network gives that network multiplier 1.0, which is catastrophic.
- C > B returns a zero prediction (504aldo F54).
- The public chains fail on non-suite inputs: 1024 × 32 (MemoryError or budget), and weights ×10 or |W| (SymmetryError) (submission).
- Implies: [C] every structural constant must be gated on the shape (n, L), with a shape-generic fallback. The smoke-test shape is where the 1/n expansion is weakest (F3.5, F10.4).

**M9. Data independence.** The FLOP bill of a fixed estimator does not depend on the data (V29: 260.06 u on every network). Implies: [–] cost can be certified once per shape.

**M10. Fair accounting** (16 Sep 2026).
- A score's benefit must come from the estimation method. Packing is banned.
- Gains from stale symmetry tags (flopscope issues #264, #265, #267) are disqualifiable, even during prize review.
- Implies: [C] only deliberately created symmetry tags may be relied on.

---

## 15. Corrections, tensions and caveats found while compiling

**Statements in BRIEF §3 that go beyond their evidence.**

- **C1. Spectral independence.** "The joint gate law has bounded unpinned spectral independence: η ≈ 1.7–5.6 at width 128 and 2–8 at width 1024."
  - The values are right.
  - Boundedness is not established: at layers ≥ 6, η₀ rises 1.2–2.1× from n = 128 to 1024, and it is still rising at layers 13–15.
  - The n = 1024 values carry a Monte Carlo bias of ≈ +0.2, and pinnings were not tested (F4.2).
- **C2. The 40 % of old content.** "About 40 % of the pairwise third-order structure feeding layer l+1 comes from content generated more than one layer earlier."
  - At 1024 this was measured on one network with N = 32,768 × 2, by a noise-free replica estimator.
  - It is the pass-through of the all-distinct part only.
  - At n = 128–256 the share rises with depth to 0.51–0.59 (F8.1).
- **C3. The 0.3 n modes.** "No low-rank representation of it below ≈ 0.3 n modes exists (widths 64–1024)."
  - Atlas measurements cover n = 64–256; n = 512 and 1024 come from the ensemble formula on fresh networks.
  - The required fraction depends on age: from 0.95 n (age 1) to 0.19 n (age 15) for a 2 % error on the content's own D21. Merged tiers at targets t ≥ 10 need 0.1–0.45 n (F9.3).
- **C4. The variance oracle.** "Giving an analytic propagation the true per-neuron pre-activation variance at every layer cuts its error by 40 %; adding the joint fourth cumulant reaches 1.17e-9."
  - The 40 % row uses 3 networks and a noisy layer-0 reference.
  - "Joint κ4" in that ladder is the per-neuron κ4 plus two chain inputs derived from it, not a joint tensor (EscAI README).
  - So the 1.17e-9 row means per-neuron variance, κ3 and κ4 at every layer (F3.2).
- **C5. Sampling.** "Monte Carlo variance ≈ 0.05 per sample; making a rule exact to Hermite degree 5 buys only ≈ 3×, degree 32 ≈ 10×."
  - These are Phase 1 values (256 × 32).
  - At 1024 the per-neuron final variance is 0.071–0.075 (F2.1, F12.1), and the chaos spectrum has not been measured.
- **C6. The n^-0.8 law.** "The residual of that first-order description falls as ≈ n^-0.8 at fixed depth."
  - It rests on deep-layer values at n = 128 (2 networks), 256 (2) and 1024 (1 network, N = 32k).
  - The 1024 closure omits B7 and uses B3 = 1 (F6.5).

**Other caveats.**

- **C7. EscAI's final covariance spectrum** (stable rank 2.57 and the rest) is the estimator's covariance, not an exact one (F2.5).
- **C8. The C-shaped fourth-order laws.**
  - The (2,2)-slice regression with R² from 0.998 down to 0.045 is on AndreasHad04's chain state.
  - 504aldo's λ-law object is the rank-one core of the augmented chain, not the (2,2) slice (F105).
  - The physical C-aligned λ ≈ 2.1–2.6e-3 is EscAI's estimate (F7.3).
- **C9. A spurious direction.** Oishi's "30.7 % of the error energy lies along the truth direction" is not used here. Their own variance decomposition (ICC ≈ 0) shows it is noise structure of a sampling estimator, not a correctable bias. The structural statement it gestures at, a common scale mode, is F2.4.
- **C10. Two Gaussian-closure values** circulate for the six bench networks: 4.30e-6 (bench, linearised cross-covariance) and 4.10e-6 (design streams, [CONVERGENCE.md](CONVERGENCE.md)). They are different implementations of the closure. The ≈ 400× gap in F5.4 holds for either.

**T1. Tension: the mean mode and old content.**
- The old-content stream reports that the leading tensor mode of the aged pool (ages > 4, slices included) lies along μ, with 60–80 % of the tensor energy but "almost no D21".
- The transfer digest reports that the rank-one projection of the transported full source onto u₁(J) holds 51–90 % of that content's own D21 energy at ages ≥ 6.
- The measures differ: the leading mode of a merged pool against the projection of single sources onto a propagator direction, and normalisation by the whole D21 against the content's own. Both are at n = 128.
- The cheapest resolution: measure the D21 carried by the μ-direction rank-one part of the merged old tier, in one convention, at n = 1024 (streaming contractions as in oracle1024).

---

## 16. Not yet measured at the competition shape

Each item matters to every design. The cheapest measurement known is given.

1. **The old-content share and the first-order closure error** at n = 1024 on more networks and at higher N (the oracle1024 production runs at N = 3.5e6 per replica). Today they rest on one network at N = 32,768 (F6.3, F8.1).
2. **The share of old content by age at 1024, and the fraction carried by the mean mode.** These are the two numbers transfer §7 asks oracle1024 for. They decide whether the mean mode is the lever that F2.3 suggests.
3. **The (2,1,1) slice at 1024**, specifically the rank of its transported effect on D21(l+1). The streaming contraction T(B6) = ⅔T(Y0) + ⅓T(Y2) of oracle1024 gives it per layer without an n³ tensor (F7.1, F7.2).
4. **The fourth-order generation law at 1024.** At n = 128 first-order diagrams leave a third of K4 unexplained at depth (F7.5). Whether that shrinks with n as the third-order law does (F6.3) is unknown.
5. **Pinned spectral independence of the face law**, which needs third-order gate statistics (F4.2).
6. **The norm-mode share of the per-sample output variance, and the chaos spectrum, at 1024.** Both are known only at P1 (F2.4, F12.3).
7. **The per-neuron λ3, λ4 and the common-scale mixture** on more than one network at 1024 (Part 1 used one).
8. **Network-to-network spread of the structural shares** (old-content share, k/n, η₀) over many networks. keenanpepper's 1,000-network atlas at N = 1e9 holds pair blocks from which D21 and K22 are derivable (plan §3.1).
9. **Old versus newly born content in the covariance's non-Gaussian correction.** It is +3.75 % of the final variance in total (F3.3), and its split by age has not been measured.

---

## Appendix A. Where each source's facts went

| source | facts |
|---|---|
| plan §3.1, §6b, §7 | F6.2, F6.3, F7.2, F3.2, F7.4, M2–M10 |
| oracle1024 (REPORT and `results/summary_prelim.md`) | F6.2, F6.3, F6.10, F7.1, F8.1, §13 |
| chain128 | F6.2, F6.3, F6.10, F6.12, F7.1, F7.2, F7.5, F10.1, F10.4, §13 |
| old-content | F2.3, F8.2, F8.3, F9.1, §13, T1 |
| coef-ensemble | F6.4, F7.1, F11.2, §13 |
| submission | F2.1, M7, M8 |
| costmodel | F6.7, M1–M6, M9 |
| theory | F3.4, F5.1, F5.2, F6.4–F6.8, F10.4, §13 |
| intel (EscAI, plus the README definitions) | F2.5, F3.1, F3.2, F4.5, F6.9, F7.3, F7.4, F10.1, F10.4, F10.5, F11.3, M7, M8, M10, C4, C7, C8 |
| Oishi (Phase 1 shape) | F2.1, F2.4, F2.5, F3.3, F3.5–F3.7, F4.1, F4.3, F10.2, F10.3, F12.1, F12.3, F12.4, C5, C9 |
| 504aldo | F2.1, F3.3, F6.1, F6.11, F6.12, F7.3, F8.4, F9.5, F10.1, F10.2, F10.4, F10.5, F11.2, F11.3, F12.1, F12.2, M1, M8 |
| transfer | F2.2, F2.3, F4.1, F4.2, F4.4, F5.3, F6.10, F8.2, F8.3, F8.5, F9.1–F9.4, F9.6, §13, C1, C3, T1 |
| bench (calibration) | F2.1, F5.4, F11.3, F12.1, C10 |
