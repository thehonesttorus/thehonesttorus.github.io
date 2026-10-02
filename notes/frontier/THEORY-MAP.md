# The frontier (504aldo V29, MIT) read in our theory: limitations, unlocks, experiments

User direction (2 Oct, 09:40 UTC): take the frontier approaches and understand them in our abstract theory. Identify their limitations and the unlocks (theoretical, computational and algorithmic), realise the unlocks as changes to an estimator, and test them with the adjusted Phase-2 MSE as the reward.

**Base.** `estimator_v29r3.py` is upstream V29 plus a robustness wrapper (504aldo, MIT, `LICENSE-504aldo-MIT`). It is the public chain the ranks 35–50 clones run.

**Its numbers.**
- On our dev set (w1024_d16, 6 MLPs): raw 1.755e-8 (truth noise subtracted) at C/B 0.2526 steady, so dev adjusted ≈ 4.4e-9. MLP 0 alone: raw 1.942e-8.
- Public board: LB 5.40e-9 (raw 2.13e-8).
- Leaders: J2W 1.6e-9 (raw 1.50e-8 at 0.110 B); marius_binner 1.7e-9 (raw 1.14e-8 at 0.151 B).
- Sources: `../digests/504aldo-k3-chain.md`, `../digests/phase2-intel-2026-10-01.md`.

## 1. What V29 carries, in our language

- **Gaussian part.** Mean and full covariance (μ, C), transported by a Strassen sandwich each layer.
- **Third cumulant: a crossed product along the age axis.**
  - One source is born per ReLU. Each source is a sum of hub blocks Sym(X ⊗ P ⊗ Y), which are our atoms with legs (Y, Y, Z) (FC's Wick legs).
  - **Ages 1–4: exact dense legs.** These are fresh layers, which G10 says no grading of neuron space can compress.
  - **Ages 5–7:** confined to one transported shared basis of rank 384, with static per-source factors.
  - **Ages ≥ 8:** a nested sub-basis of rank 224.
  - This is a two-step stand-in for the Dixmier allocation k(a) ∝ n/a, which team G showed is the score optimum. In operator terms it is AF blocks along the age shift, as with the odometer.
- **Readout.** Each source enters D21 through Hadamard products of its own legs, i.e. the masa conditional expectation E_D. Team G's commuting-square failure is why this is costly: "the bill is set by the number of dense n × n legs formed per layer".
- **Slices and feedback.** Thin legs of rank 16 (residual (2,1) leg and D21 feedback).
- **Fourth cumulant.**
  - The diagonal is exact.
  - The off-diagonal core is memoryless, G_off = λ_l C_off. That is the CAP (dilation) sector of κ4 only. λ is fitted per layer and is ≈ 3× the physical value (EscAI), so its gain is a cancellation.

**Cost anatomy (260 u).** Young K3 116 u (44.5 %), old K3 107 u (41 %), thin legs 25 u, covariance 7 u, births 6 u.

## 2. Limitations

| | limitation | theory | measured |
|---|---|---|---|
| L1 | Old-tier cost is linear in the number of old sources (≈ 1.1 u per old source-layer; joins ≈ 3.5 u each) | the old memory's continuation space is small and fixed (rational memory): per-target readout rank 43–116 of 1024 at 99 % (H §7); one cohort basis of k = 64 suffices for ages ≥ 8 (team F S1) | 107 u = 41 % of the bill; the author: "whatever carries [the leaders'] old-source content costs them ~nothing" |
| L2 | κ4 off-diagonal is CAP-only and memoryless | the sector split predicts a FREE (mixing) κ4 sector, geometrically decaying, so cheap to carry for young ages | EscAI: on an exact base λ hurts; the author: "the top methods carry fourth-cumulant content beyond a memoryless core" |
| L3 | Errors sit in the per-layer pre-activation marginals | fresh-weight lemma: Δv_a = w_aᵀ ΔC w_a has a coherent part (2/n) tr ΔC and an incoherent part ≈ 2√2 rms ΔC_off, so the marginals inherit every off-diagonal covariance error | EscAI oracle: true variance, κ3 and κ4 per layer give 20× (2.3e-8 → 1.2e-9); variance alone −40 % |
| L4 | Young tier 116 u (4 sources × ≈ 4 dense products) | G10: fresh layers cannot be compressed in neuron space, only along the age axis | ages 1–2 exact ≈ 0.12 B in FC's form (H §11) |

## 3. Unlocks turned into experiments (reward: dev adjusted = raw × max(0.1, C/B), paired MLPs)

- **E1. Age-axis allocation by knobs** (coordinator, `sweep.py`).
  - Knobs: AGE_OLD ∈ {2, 3, 4}, R_OLD, AGE_OLD2, R_OLD2.
  - These are Dixmier-guided moves: an earlier join with a larger basis, and a smaller nested rank.
- **E2. Finite-state deep memory** (child session).
  - Replace the per-source nested tier (ages ≥ AGE_OLD2) by one cohort state whose cost does not grow with the number of sources.
  - Candidate: a Tucker core in a transported basis, cores added at join without a fit, k ≈ 64–128. The core readout costs 2 n k³.
- **E3. κ4 beyond memoryless and the marginals** (child session).
  - Diagnose which marginal errors dominate on dev, via an oracle against MC truth per layer.
  - Then replace or augment the λ C_off channel with a physically grounded κ4 carrier with memory, e.g. region's mean-field κ4 recursion with (2,2) column means and the exact diagonal, or a FREE-sector κ4 for young ages.

## 4. Where the remaining accuracy lies (E3), in the paper's terms and ours

- **ARC's cumulant-propagation paper** (arXiv 2605.05179; digest `../digests/paper-2605.05179-cumulant-propagation.md`):
  - MSE ≈ c_K (L/n)^K.
  - For odd K, the paired entries of κ_{K+1} are as large as the generic κ_K entries. The *augmented* K = 3 algorithm therefore tracks every harmonic component of κ4 except the **traceless** one.
- **V29 is an augmented K = 3 chain.** Its regenerated κ4 core G = diag(g) + λ C_off (R² 0.90–0.97 against the aug chain's transported core) is the r = 1 harmonic (matrix-core) projection.
- **The dictionary.** The non-traceless harmonic parts are our **CAP sector** (O(n)-invariant: trace and dilation). The traceless part is our **FREE sector** (through-strings).
- **So "κ4 beyond a memoryless core" (the author's diagnosis of the leaders) means the FREE κ4 sector.** Our sector theory predicts it decays geometrically in age, like the FREE K3 sector (≈ 0.65 per age, costate §8). It is then concentrated in the youngest births, so a young-only FREE-κ4 carrier is the cheapest way to move toward K = 4 accuracy. Full factorised K = 4 costs 24 n⁴L², which is impossible at the budget.
- **Cost estimate.**
  - One FREE κ4 source has 4-legged atoms, ≈ 4–6 dense products per layer (≈ 2.2–3.3 u at Strassen L5).
  - Ages 1–2: ≈ 5–7 u per layer, ≈ 75–100 u in total. Affordable only once the K3 bill is cut (F1, F2).
- **EscAI's oracle** points the same way. True per-layer marginals (variance, κ3, joint κ4) give 20× in raw, and variance alone −40 %. By the fresh-weight lemma, per-neuron variance errors are random quadratic forms of the full post-activation covariance error, which the missing κ4 (2,2) and (3,1) content feeds through the bivariate closure.
- **Status:** not attempted today; the reward gain from cost is larger and surer at the floor (see E1, F1, F2).

## 5. The up-down involution on the age tree (Pearson–Bellissard), and what it says about V29 (2 Oct, 10:50 UTC)

Spark from the user (the other chat's Pearson–Bellissard / Lapidus note), transplanted to the age axis:

- **The tree and F.** The Cantor set of histories is the age tree (ages = cylinder depth). F swaps two siblings at each vertex, F² = 1, and [F, f] is the Haar detail of the readout f along the tree, i.e. the change in the estimate when one more age is resolved separately instead of merged with its sibling.
- **Dixmier trace = measure.** The Dixmier trace at the critical exponent is the measure, here the 1/a density on ages (team G's score optimum). So the score-optimal resolution of an age block is proportional to its measure, log(a_end / a_start).
- **Self-similar string.** Fibonacci blocks [F_k, F_{k+1}) all have measure → log φ, so they form a self-similar string of ratio φ, with equal Dixmier weight per block and per-source rank ∝ 1/a. Team E found the memory's age-graded rank class to be Fibonacci/Penrose.

**What this explains in the frontier.** V29's empirically tuned old tiers *are* this string:

| tier | ages | V29 rank | Dixmier c·n/a at c = 2.2, block average |
|---|---|---|---|
| 1 | 5–7 | 384 | ≈ 382 |
| 2 | 8–12 (sources older than 7 transports) | 224 | ≈ 220 |

- The tier boundaries are Fibonacci ages and the ranks are Dixmier ranks with one constant, c ≈ 2.2.
- The young ages 1–4 stay dense. Dixmier would give 751 / 563 at ages 3 / 4, but a per-source frame costs ≈ 7k/n dense products per layer against 4 dense ones (the Hadamard wall), so it does not pay above k ≈ 585.
- **The next block** is ages 13+ at k ≈ 160 (passed to F2). Depth 16 contains only one decade of age, so deeper blocks are worth little here: the Lapidus oscillations (complex dimensions 2πik/log φ) need several periods of log-age to matter.
- **Haar sparsity ("sharing").** [F, f] vanishing across siblings is the statement that a cohort can be merged without loss. F2's cohort Tucker core is that quotient-tree carrier, and it is exact (the core is additive), so its only approximation is the basis.

**E3 probe (the [D21 ⊗ C] hyperedge, the "refinement by one more layer" term of the newborn).**
- *Setup.* V29 carries it at rank 16 with first-order coefficients. Scaling it by s at birth (`estimator_fbs.py`, H_FB_SCALE), MLPs 0–2, zero cost.

  | s | 0.8 | 1.0 | 1.2 | 1.5 | 2.0 |
  |---|---|---|---|---|---|
  | raw (×1e-8) | 1.507 | 1.467 | 1.430 | **1.409 (−4.0 %)** | 1.432 |

- *Reading.* Rank truncation loses energy, and depth renormalises the first-order coefficient (competition-plan atlas: the D21 diagram coefficients drift with depth).
- *Principled form, no fitted constant.* Rescale the captured part by its own truncated energy, (‖D21‖_F / ‖Q_f B_f‖_F)^γ at n² cost. Running as `fbg050` / `fbg100`.
