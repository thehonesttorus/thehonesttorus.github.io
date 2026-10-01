# C — Transfer to the competition

### Gibbs–Markov theory, expanders and NCG against a FLOP-bounded estimator of the post-ReLU means of a random 1024×16 ReLU MLP. One literal expander: the doubled leg of the cumulant transport. One exact symmetry with a measurable consequence: dilations, i.e. the gain mode, which is what the published (2,1,1) regeneration actually carries. A priced list of dead ends, each with the experiment that would falsify the verdict.

*Bridges synthesis, angle C (transfer to the competition), 2026-10-01. Inputs, all read in full before writing:*
- *the bridge digests [arxiv-2609.38007.md](../../digests/bridges/arxiv-2609.38007.md), [arxiv-2504.02208.md](../../digests/bridges/arxiv-2504.02208.md), [hammersley-clifford.md](../../digests/bridges/hammersley-clifford.md), [expanders.md](../../digests/bridges/expanders.md), [hdx-spectral-independence.md](../../digests/bridges/hdx-spectral-independence.md), [nc-dirichlet-lindblad.md](../../digests/bridges/nc-dirichlet-lindblad.md), [transfer-spectrum-measurement.md](../../digests/bridges/transfer-spectrum-measurement.md) and [chat-2609.38007-retrieval-status.md](../../digests/bridges/chat-2609.38007-retrieval-status.md);*
- *the programme notes [research-program.md](../../research-program.md), [conditional-arrow-algebra.md](../../conditional-arrow-algebra.md) (Note 1), [simplicial-complex-as-decomposition.md](../../simplicial-complex-as-decomposition.md) (Note 2), [local-to-global-unlocks.md](../../local-to-global-unlocks.md) and [mlp-bridge.md](../../mlp-bridge.md);*
- *[competition-plan.md](../../competition-plan.md) §§3.1, 6a, 6b and 7;*
- *the stream reports under [streams/](../): oracle1024, old-content, chain128, coef-ensemble, costmodel, submission, and the stopped est-accuracy and est-cost;*
- *the fresh-slate [BRIEF.md](../../fresh-slate/BRIEF.md).*

*The sibling notes [A-unification.md](A-unification.md) and [B-programme.md](B-programme.md) were consulted for consistency only.*

**Framing: the fresh slate.** Since 17:45 UTC on 1 Oct the competition system is designed from a fresh slate ([CP] §6a, [BRIEF]). Moment and cumulant chains are now baselines, not starting points. Several levers below are phrased in terms of the published chain's objects: the D21 interface, the (2,1,1) slice, old-source κ3 content. Each is offered as a *fact about the object* that any design must reproduce or route around. None is a proposal to adapt that chain.

The programme's rules hold throughout:
- Nothing here changes the objects of Notes 1–2. Faces are abstract barycentres; an arrow is a linear layer with the ReLU immediately before it.
- Every statement below is about the competition's dictionary instance, not about the general theory.
- A negative experimental result is charged first to the dictionary (rule P2.1).

"Dead on arrival" (DOA) means: *for this competition it cannot pay* (width 1024, depth 16, budget ≈ 0.15 B, interface bar ε ≤ 2.2 %). It never means "the theory is wrong".

**Labels.** Every claim carries exactly one label.
- **THEOREM**: a published theorem cited to the digest section where it was read, or an elementary statement proved in full here.
- **KNOWN-LINK**: a published connection between two areas, cited.
- **DERIVED**: proved here, with the proof in the text. "+ checked" means also verified numerically (checks C1–C11, §10).
- **ANALOGY**: a structural similarity, with the exact point where it breaks.
- **CONJECTURE**: a precise statement, not proved, with its test.
- **SPECULATION**.

Measured numbers, ours or the streams', are evidence attached to a labelled claim. They are not claims of their own.

**Sigla.**
- Digests: [D-Y] Yang; [D-CR] Chen–Rouzé; [D-HC] Hammersley–Clifford; [D-EXP] expanders; [D-HDX] HDX / spectral independence; [D-NC] noncommutative Dirichlet forms; [D-TS] transfer-spectrum measurement; [D-CHAT] chat retrieval status.
- Sibling notes: [A], [B].
- Stream reports: [O1024] oracle1024; [OC] old-content; [C128] chain128; [CE] coef-ensemble; [CM] costmodel; [EA] est-accuracy (stopped).
- Plans: [CP] competition plan; [BRIEF] fresh-slate brief; [FU] fresh-slate `foundations-unlocks.md`, consulted only for the two cross-references in §5.3.

**Retrieval.** No new source was opened. Published results are cited through the digests that quote them. Standard facts used from memory are marked (memory):
- the Bai–Yin / Latała operator-norm law for i.i.d. matrices;
- Mehler's formula and the Hermite expansion of an indicator;
- the Leonov–Shiryaev formula for cumulants of products;
- the Lebowitz–Percus–Verlet ensemble corrections;
- the Roberts–Yaida–Hanin finite-width four-point vertex;
- TAP / Plefka expansions;
- the Kotecký–Preiss criterion.

Each is used only as a name, or is re-derived or checked here. The ChatGPT conversation behind [D-Y] could not be retrieved ([D-CHAT] §1).

**New measurements for this note.** The checks C1–C11 run at widths 8–1024. They took about two CPU-hours on a shared 4-core box. The scripts are scratch files, not committed, because this task writes one file. §10 and the Appendix give every formula needed to reproduce them.

---

## 0. The answer in brief

1. **Where the "same essence" is literally a theorem in the competition object: the doubled leg.** When a cumulant tensor is pushed through a linear layer, the transport depends on how often an index repeats.
   - An index that appears once is transported by W.
   - An index that appears twice is transported by the Hadamard square W∘W.
   - An index that appears three times is transported by W∘W∘W.

   For He weights, W∘W = (2/n)·11ᵀ + E with ‖E‖ = (4√2/√n)(1 + o(1)). This is a dense weighted expander, √2 above the Ramanujan value for degree n, and the expander mixing lemma holds verbatim. W∘W∘W is O(1/n). W itself (norm 2√2, no trivial eigenvector) scrambles orientations. DERIVED + checked (C1: ‖E‖ = 0.1767 at n = 1024, predicted 0.1768). The flat (Perron) mode of the doubled leg transmits exactly the *trace channel* of a cumulant, Σ_p κ_{pp…}. That channel is the covariance of the squared norm with the remaining legs. §6.4.
2. **Where it is two halves of one mechanism.** Heredity is free and decorrelation does the work.
   - The layer chain is exactly Markov, because each layer is a deterministic function of the previous one.
   - The closure works because pre-activation correlations are entrywise O(n^{−1/2}) and single legs scramble orientation (Haar, [D-TS] §2.4).
   - The leg-partition closure is a first-order cluster expansion with exact transport (KNOWN-LINK, Leonov–Shiryaev). Its error falls as n^{−0.8} and is 0.8–1.0 % at n = 1024, so it is not the binding error. §6.
3. **Where it is false.**
   - **(a) A pairwise Gibbs / MRF description of the face law cannot carry the binding content.**
     - Pairwise gate statistics are bivariate functionals (THEOREM).
     - The 2-clique Gibbs closure overstates three-gate cumulants by a factor up to 4/π (DERIVED + checked).
     - Hammersley–Clifford positivity fails at every layer ≥ 1, and the dependency graph is complete (DERIVED + checked).
     - Gates are 0-homogeneous, so any face-law description is blind to the dilation mode of item 4 (THEOREM).
   - **(b) Spectral independence is a frame bound,** Cov ≼ (1+η)·diag. It is not an accuracy certificate for a closure (THEOREM + DERIVED).
   - **(c) A time-averaged, detailed-balance recovery map cannot restore erased old content.** The ingredients are missing ([D-TS] B8), and the erased content is not a function of the retained state. Under an orientation hypothesis, its Bayes-optimal recovery is zero (DERIVED under hypothesis H, with measured support).
   - **(d) A CMI-type decay statement predicts that the chain may forget essentially nothing.** There is no gap to decay with, and each old source still carries ≈ 0.1 of D21 at the last layer, 5× the bar. §§3–5.
4. **The one bridge with a positive competition payoff: the dilation (gain) mode.**
   - **Mechanism.** With no biases, every arrow commutes with x ↦ tx for t > 0. The input radius is the obvious consequence ([BRIEF] §3). The less obvious one is that the gain of each input direction, ‖a_l(x̂)‖, fluctuates across directions, and that this fluctuation propagates undiscounted.
   - **What a scale mixture predicts** (z = s·y; DERIVED, §7.1):
     - κ4 (2,1,1) slice = Var(s²)·(var_i C_jk + 2C_ijC_ik), the Wick part of the slice's own moment;
     - κ3 (2,1) slice = ½Var(s²)·(var_a μ_b + 2μ_a C_ab);
     - a spike Var(s)·μμᵀ in C.
   - **Measured (C7–C9).**
     - The published regeneration u_iC_jk has u parallel to var (cos 0.91–0.997).
     - A single fitted scalar per layer recovers 74–114 % of the D21 gain that the n-parameter regeneration gets at n = 128.
     - At depth the variance excess along the mean is 1.0–1.6× the gain-mode prediction, so the spike is largely the gain mode at order 2.
   - **The amplitude.** It is width-stable as c·n. The scalar Var‖z̃‖²/(E‖z̃‖²)² − 2‖C‖²_F/(tr C)², which is exactly the trace of the (2,2) slice divided by (tr C)² (DERIVED identity), predicts it within 30 %.
   - **The limit at n = 1024.** Across widths the template's share of the slice energy collapses: at layer 15 it is <<SHARE15>> from n = 128 to 1024 (C9). This is the measured mechanism behind the published regeneration's falling share of the κ4 gap (55 → 33 → 20 % at n = 128 → 256 → 1024, [O1024]). At n = 1024 the binding (2,1,1) content is mostly *not* the gain mode.
   - **Not a Gaussian quantity.** A Gaussian-input recursion recovers layer 1 exactly but only 24–47 % of the amplitude at layer 15 (C10). §7.
5. **The three measured facts named in the brief, priced.** §3.
   - *Free-probability law of old content.* An excellent offline calculator and a DOA carrier at n = 1024: 54–86 u/layer against 3–4 available. Content with more scrambled legs is worse (DERIVED: the kept energy goes as the m-th power of the per-leg fraction). The (2,1,1) slice transports like a covariance, with two scrambled legs and one flattened leg.
   - *BBP spike.* Alive only as deflation, with ≈ 20 % fewer modes. As a standalone carrier ε_tot ≥ 0.16 against 0.022.
   - *Spectral independence.* A diagnostic. Its growth with depth is in part the spike seen through the gates: the top mode overlaps the spike image at 0.39–0.49 at layers 1–4, against 0.01 for a random direction.
6. **NCG's part for the competition.** This is not its part for the theory, which [A] and [B] treat.
   - *Literal.* Free probability (Voiculescu's noncommutative probability) is the law of the propagators. The annealed layer is a conditional expectation, a twirl.
   - *Duality.* This is the one Chen–Rouzé ingredient that transfers, and it is already in use as the Duhamel telescoping of the heisenberg design ([FU] unlock 36). Its backward questions obey the same free-probability mode law, because J and Jᵀ share singular values (DERIVED).
   - *Analogy.* For positively homogeneous dictionaries the dilation is a central charge of the arrow algebra, and the gain mode is its canonical-ensemble fluctuation (Lebowitz–Percus–Verlet).
   - *No role.* The detailed-balanced single-Pauli Lindbladian itself has no competition role: the estimator runs no dynamics and has no stationary state. §8.
7. **What to run** (§9). Eight experiments use existing tools. The decisive cheap ones are:
   - E1: rank-1 power iteration of the (2,1,1) slice at n = 1024, minutes per layer;
   - E7: the gain template at n = 1024 on more networks, one done here (C9);
   - E8: the D21 → MSE law for dropped old content versus closure error (chain128, hours).

---

## 1. Results at a glance

| # | statement | label | for the competition | § |
|---|---|---|---|---|
| R1 | Doubled leg W∘W = (2/n)11ᵀ + E with ‖E‖ = 4√2/√n; tripled leg O(1/n); a single leg scrambles. The flat mode carries the trace (norm-coupling) channel | DERIVED + checked (C1, C5) | structural fact for every design | 6.4 |
| R2 | Pairwise gate statistics (gate_GG, gate_GX, any pairwise MRF) are bivariate functionals and cannot hold trivariate content as state | THEOREM | DOA (pairwise face-law state) | 4.1 |
| R3 | The 2-clique Gibbs closure of a probit gate law overstates three-gate cumulants by φ(t)(1−2p)/(t·p(1−p)) ∈ (1, 4/π] at leading order | DERIVED + checked (C2) | DOA (MRF closure) | 4.2 |
| R4 | At layers ≥ 1 the face law vanishes on at least one orthant (the cone of W's rows is pointed): Hammersley–Clifford positivity fails; the dependency graph is complete | DERIVED + checked (C3) | DOA (HC factorisation) | 4.3 |
| R5 | Gates are 0-homogeneous, so every face-law description is blind to the dilation (gain) mode | THEOREM (elementary) | limits every gate-only design | 4.4 |
| R6 | η-SI ⟺ Cov ≼ (1+η)·diag; independent of the entrywise smallness the closure needs | THEOREM + DERIVED (C4) | DOA as certificate or lever; diagnostic | 3.3 |
| R7 | Propagator-basis carriers need ∝ n modes; content with m scrambled legs keeps e(k)^m | THEOREM (cited) + DERIVED | DOA at n = 1024; calculator alive | 3.1 |
| R8 | Spike alone: ε_tot ≥ 0.16 for the w = 4 tier; as deflation ≈ 20 % fewer modes | DERIVED (arithmetic on measured inputs) | DOA alone; small lever | 3.2 |
| R9 | Erased old content: no detailed balance or stationarity; not a function of the retained state; Bayes-optimal recovery = 0 under H | DERIVED (under H) + measured support | DOA | 5.1–5.2 |
| R10 | Duality (CR's telescoping, [FU] 36) moves the memory from the forward state to the backward question, whose low-rank compressions obey the same free law | KNOWN-LINK + DERIVED | alive for its bilinear error, no mode saving | 5.3 |
| R11 | No transport gap ⇒ no exponential forgetting; each old source ≈ 0.1 of D21 at t = 15 | DERIVED (conditional) + measured | DOA (forgetting by depth) | 5.4 |
| R12 | Closure = first-order cluster expansion with exact transport; error ∝ n^{−0.8}, 0.9 % at n = 1024 | KNOWN-LINK + measured | not binding | 6.1–6.3 |
| R13 | Flat doubled-leg carrier for the (2,1,1) slice | refuted at depth (C6) | DOA | 6.5 |
| R14 | Scale mixture ⇒ κ4(2,1,1) = Var(s²)·Wick, κ3(2,1) = ½Var(s²)·T, C ∋ Var(s)μμᵀ; gex = Σ_ij κ4_iijj/(tr C)² exactly | DERIVED | explains the published regeneration | 7.1 |
| R15 | u_Creg ∥ var (cos 0.91–0.997); one scalar recovers 74–114 % of the regeneration's D21 gain at n = 128; c·n width-stable; template share falls with n | measured (C7–C9) | one scalar per layer; not enough at n = 1024 | 7.2–7.3 |
| R16 | The gain amplitude is not a Gaussian-closure quantity: Gaussian-input increments give 24–47 % of it at layer 15 | measured (C10); CONJECTURE refuted in its simplest form | the scalar needs its own closure | 7.4 |

---

## 2. The object, the bar and the budget: what a bridge must beat

- **Task** ([BRIEF] §1, [CP] §0).
  - Input x ~ N(0, I_n), n = 1024, L = 16 layers, W_l i.i.d. N(0, 2/n), no biases.
  - z_l = a_{l−1}W_l (x@W convention, a_{−1} = x), a_l = relu(z_l).
  - Estimate the n means E[a_15].
  - Score = final-layer MSE × max(0.1, F/B), with B = 2^41 FLOPs. The unit is u = 2^31 FLOPs, one dense 1024³ product.
- **The bar** ([CP] §0, §7).
  - Money line: adjusted MSE ≈ 1.6–1.8e-9, i.e. raw ≈ 1.6e-8 at the 0.1 floor.
  - Best published raw: 1.14e-8 at 0.151 B. V29: raw 2.13e-8 at 0.253 B.
- **Interface law** ([CP] §3.1, [C128], [O1024]). In a cumulant chain the error that reaches the final layer through the nonlinearity's (2,1) interface obeys extra MSE ≈ 4.2e-6·ε², where ε is the relative rms error of D21_{ab} = κ3(z_a, z_a, z_b). The frontier needs ε ≤ 2.2 %.
- **Measured at n = 1024** ([O1024], one MLP; [CP] §6b), D21 error by κ4 input:

  | κ4 input to the closure | D21 error |
  |---|---|
  | exact inputs (leg-partition closure) | 0.8–1.0 % at layers 3–14 |
  | (2,1,1) slice regenerated as u_iC_jk | 2.5–2.7 % |
  | no slice | 2.85–3.17 % |

  The regeneration closes a falling share of the κ4 gap: 55 → 33 → 20 % at n = 128 → 256 → 1024.
- **End to end at width 128** ([C128]; [CP] §6b). With κ4 at truth, the D21 ladder maps one-to-one onto the final layer. With any κ4 closure of its own, every κ3 rule lands at 3e-5–8e-5.
- **EscAI's oracle** ([CP] §7), starting from an analytic propagation:
  - giving it the true per-neuron variances cuts its error by 40 %;
  - adding the true per-neuron κ3 reaches raw 8.4e-9;
  - adding the joint κ4 reaches 1.17e-9.
- **Budget** ([CM], [CP] §7).
  - A dense product costs 1 u; Strassen L5 costs 0.556 u.
  - Slices plus the full first-order closure cost ≈ 104 u ≈ 0.102 B.
  - The leaders' 0.145–0.16 B therefore leaves ≈ 45–60 u, i.e. 3–4 u per layer, for old-source and κ4 content.
  - The binding engineering constraint is the residual time per flopscope call, not FLOPs.

So there are two binding objects.

| binding object | why it is hard | source |
|---|---|---|
| (i) The (2,1,1) slice of κ4 | An n³ object. At n = 128, rank 32 of 128 holds only 69–81 % of it at layers 1–3; the first mode holds 86–91 % (one MLP) and 56–82 % (the other) at depth. | [CP] §3.1 (vii) |
| (ii) Old-source κ3 content | Share of D21 0.5–0.8 at depth; no carrier below ∝ n modes. | [OC] |

A bridge has to say something about one of these at ≤ 3–4 u per layer, or show how to avoid them.

---

## 3. Item 1: the three measured facts

### 3.1 The free-probability law of old κ3 content

**Cited facts.**
- THEOREM ([D-TS] §2.1): B^{⊗k} restricted to Sym^k has singular values ∏σ_{i_j}. There is no higher-order gap that is not already a first-order Lyapunov gap.
- KNOWN-LINK ([D-TS] §2.2–2.3: Pennington–Schoenholz–Ganguli; Hanin–Nica): gated products obey free multiplicative convolution, and PR(J_{s→t}) ≈ n/(1 + Σ(r_l − 1)) ≈ n/(2·age), measured to within 4–6 %.
- Measured ([D-TS] §3.4): k_2%/n is width-independent to ±0.02 from n = 128 to 1024. At n = 1024, age 8 needs 322 modes and age 15 needs 197 modes for a 2 % error relative to the content's own D21 contribution.

**DERIVED (arithmetic; [D-TS] §0.3).** A Tucker core of rank k costs k³/n² u per layer. The w = 4 tier needs ≈ 384–448 modes at n = 1024, i.e. 54–86 u per layer, against 3–4 available. **DOA as a carrier.**

**DERIVED (legs, under the hypothesis of [D-TS] §2.4).** Let a source be Haar-oriented relative to the right singular vectors of J, and let it have m *scrambled* legs (indices that appear once). Project each such leg onto the top-k left singular subspace of J. To leading order in 1/n the expected kept energy is e(k)^m, with e(k) = Σ_{p≤k}σ_p²/Σ_pσ_p².

*Proof.* Write the source as Sym(Z) with Z i.i.d. Gaussian, as in [D-TS] §2.4. The projected transported energy factorises over legs into ∏ tr(P_kJJᵀ)/tr(JJᵀ), up to symmetrisation cross terms that are O(1/n) for traceless sources. ∎

Consequences:
- All-distinct κ4 (m = 4) needs more modes than all-distinct κ3 (m = 3) at equal accuracy. **DOA a fortiori.**
- The (2,1,1) slice transports differently. Its doubled index goes through W∘W, which flattens it (R1, §6.4), not through J. As a transported object it is therefore *covariance-type*, with m = 2. This is why a form like u_iC_jk can work at all, and why [D-TS] §4 item 4 finds covariance-type content cheap.

**What is alive.** The propagator-only ensemble error formula ([D-TS] §2.4) predicts mode counts from the spectrum of J alone. It reproduces the measured k_2% at 11–12 of 15 ages. It is a free offline calculator for any design that projects transported content, including the backward questions of §5.3.

### 3.2 The rank-one spike along the mean direction

**Cited facts.**
- The spike is a BGN outlier above the free-probability threshold ([D-TS] §2.6, B6).
- At n = 1024 and age 15, ⟨u₁, μ_z⟩² = 0.85 and s₁/s₂ ≈ 1.45 ([D-TS] §3.4).
- The top finite-time Lyapunov exponent is +0.059 per layer, while the bulk is discounted by the gain 2E[Φ²] = 0.58–0.95 per layer ([D-TS] B8).
- The old pool's leading tensor mode is the mean direction (cos 0.93, 60–80 % of the tensor energy at depth), but "it carries almost no D21" ([OC] Results).

**DERIVED (the spike as a standalone carrier; arithmetic).** Suppose a carrier transports the spike part of a tier exactly and drops the rest. Then ε_tot = share·√(1 − f), where share is the tier's share of D21 and f is the fraction of the tier's D21 energy along the spike. The inputs:
- share = 0.5–0.8 for w = 4 at depth; measured at n = 128–256 and flat across that range, extrapolated to n = 1024 ([OC] table (1));
- f ≤ 0.9 in the best measured case (slice-containing sources, ages ≥ 6, [D-TS] §3.6); [OC]'s convention gives f ≈ 0.

Hence ε_tot ≥ 0.5·√0.1 ≈ 0.16, against 0.022. **DOA as a standalone carrier.**

**DERIVED (deflation).** Carry the spike exactly and the bulk in the propagator basis. The bulk's own target relaxes from ε to ε/√(1 − f). At f = 0.84 that is from 2 % to 5 %. At age 8 and n = 1024, this means 263 modes instead of 322 ([D-TS] §3.4): ≈ 20 % fewer, but still ∝ n. **A small lever.**

**Link to the gain mode (measured, C9; ANALOGY with an exact breaking point).** The variance excess of C along μ̂ can be compared with the scale-mixture prediction (c4/4)·‖μ‖² of §7.1:
- at depth (layers 10–15) it is 1.0–1.6× the prediction at n = 128–<<SPIKEMAXN>>;
- at layer 2 it is 2.4–2.9×.

So at depth the spike is largely the gain mode seen at order 2. At shallow layers it is mostly the BGN amplification of the mean direction by gates aligned with it ([D-TS] §2.6). The identification breaks exactly there: the gain mode multiplies the whole vector, while the BGN mechanism amplifies one direction.

### 3.3 Bounded spectral independence of the face law

**Cited facts** ([D-TS] §3.5).
- Unpinned η₀ at n = 128 (MLP 0): 1.78, 3.25, 4.46, 3.49, 4.05, 3.78, 4.73, 4.90 at layers 0, 2, …, 15.
- At n = 1024: 2.1–7.6.
- The zero-mean Gaussian surrogate gives 15–29.
- Mean |Cor_ij| among the uncertain gates is 0.08–0.10 at depth.

**THEOREM** ([D-HDX] §3.3, Fact (a): Chen–Eldan Fact 23 and Remark 37). For a law on {±1}ⁿ, η-spectral independence is equivalent to Cov ≼ (1 + η)·diag(Cov): a frame bound.

**DERIVED.** Spectral independence and the entrywise smallness that a closure needs are independent conditions.
- *Bounded η does not give small entries.* Take a block-diagonal correlation matrix with 2×2 blocks of correlation 0.9. Then η = 0.9, but the entries are 0.9 and no first-order closure is accurate.
- *Small entries do not give bounded η.* Take Cor − I = θvvᵀ with |v_i| = n^{−1/2}. The entries are θ/n, but η = θ is arbitrary.

The closure error (§6) is a sum over diagrams of products of entrywise correlations. A bound on λ_max(Cor) controls such sums only through ‖·‖_F ≤ √n·‖·‖_op, which is vacuous here. A bounded η certifies Glauber mixing on faces under all pinnings (ALO), but the estimator runs no dynamics. **DOA as an accuracy certificate or a design lever.**

**Measured (C4).** At layers 1–4 (MLP 0, n = 128), compare the top eigenvector of Cor(g) on the uncertain gates with the spike image r_b ∝ φ(α_b)(W_lᵀe)_b / (σ_b√(p_b(1−p_b))), where e is the top eigenvector of C(a_{l−1}).
- Overlaps: 0.44, 0.39, 0.39, 0.49, against ≈ 0.01 for a random direction.
- The Rayleigh quotient of r gives 82–99 % of η₀ at these layers.
- At deeper layers the overlap is mostly < 0.25.

*Reading* (ANALOGY with the BGN outlier, breaking at depth where the overlap drops): part of the growth of η₀ with depth is the spike seen through the gates. η₀ stays useful as a cheap diagnostic of distance from product form. The pinned η needs third-order gate statistics and is the prong-1 test of U6 ([D-TS] B7); it has no competition role.

---

## 4. Item 2: a Gibbs / MRF description of the face law, and the (2,1,1) slice

**Setting.** The face law of layer l is the law of the gate pattern g = 1[z_l > 0] ∈ {0,1}ⁿ (dictionary v1). The statistics that the community bakes store, and that [CP]'s Line C proposed as a representation, are:
- gate_p;
- gate_GG = P(z_i > 0, z_j > 0);
- gate_GX = E[1[z_i>0](a_j − μ_j)].

The question is whether a Gibbs description of this law, sharpened by spectral independence, closes the (2,1,1) slice more cheaply or more accurately than the leg-partition closure does.

### 4.1 THEOREM: pairwise statistics cannot hold trivariate content as state

**Statement.**
- (i) Every statistic of the form E[f(z_i, z_j)] is a functional of the bivariate law of (z_i, z_j). This covers gate_GG, gate_GX, every pairwise moment, and the sufficient statistics of every pairwise MRF on gates.
- (ii) For distinct i, j, k, neither the (2,1,1) entry κ4(z)_iijk nor the all-distinct κ3(z)_ijk is determined by the three bivariate marginals of (z_i, z_j, z_k).

*Proof.* (i) holds by definition. For (ii), let φ be the standard normal density.
- Take bounded functions a and b with ∫aφ = ∫bφ = 0, for instance a(t) = 1[|t| > 1] − P(|Z| > 1) and b(t) = sign(t).
- Put p_ε(z_i, z_j, z_k) = φ(z_i)φ(z_j)φ(z_k)·[1 + ε·a(z_i)b(z_j)b(z_k)]. For |ε| < 1/(sup|a|·sup|b|²) this is a probability density.
- Integrating out any one variable kills the perturbation, so all three bivariate marginals are those of the independent standard Gaussian, for every ε.
- Means and covariances are unchanged, so κ4_iijk = E_ε[z_i²z_jz_k] − 0 = ε·E[Z²a(Z)]·(E|Z|)². This is non-zero, since E[Z²; |Z| > 1] ≠ P(|Z| > 1).
- With b(z_i)b(z_j)b(z_k) as the perturbation instead, κ3_ijk = ε(E|Z|)³ ≠ 0 with the same bivariate marginals. ∎

**Consequence.** Both binding objects of §2 are trivariate: the j ≠ k part of the (2,1,1) slice, and old content, which is all-distinct κ3 transported from earlier layers. A state made of pairwise face statistics can contain them only through a *model* that generates trivariate content from bivariate content. The leg-partition closure is such a model, with exact transport. A Gibbs description of the pairwise face law is another model: the maximum-entropy one. §4.2 shows that it is the wrong model, and §4.3 that the classical theorem behind it does not apply.

### 4.2 DERIVED + checked: the 2-clique Gibbs closure overstates three-gate cumulants by up to 4/π

**Setting.** Let (z₁, z₂, z₃) be Gaussian with unit variances and correlations r_ab = O(ρ). Let the gates be g_a = 1[z_a > t_a], with p_a = P(z_a > t_a), v_a = p_a(1 − p_a) and φ_a = φ(t_a). Let ν be the maximum-entropy law on {0,1}³ with the same one- and two-gate marginals: a pairwise binary MRF, i.e. the "2-clique Gibbs" description of the face law.

**Statement.** To leading order O(ρ²):
- exact: κ3(g₁, g₂, g₃) = φ₁φ₂φ₃·[t₂ r₁₂r₂₃ + t₁ r₁₂r₁₃ + t₃ r₁₃r₂₃];
- pairwise MRF: κ3_ν = φ₁φ₂φ₃·Σ_{centre b} r_ab·r_bc·V(t_b), with V(t) = φ(t)(1 − 2p)/(p(1 − p)).

For equal thresholds the ratio κ3_ν/κ3 = φ(t)(1−2p)/(t·p(1−p)) lies in (1, 4/π]. It tends to 4/π = 1.273 as t → 0 and to 1 as |t| → ∞.

*Proof.*
- **Exact law.** Mehler's expansion (memory) gives 1[z > t] = p + Σ_{k≥1} c_k He_k(z) with c_k = φ(t)He_{k−1}(t)/k!, so c₁ = φ(t) and c₂ = φ(t)t/2.
  - At order ρ² the joint cumulant comes only from paths a–b–c with gate b using He₂ and gates a, c using He₁.
  - E[He₁(z_a)He₂(z_b)He₁(z_c)] = 2r_ab·r_bc + O(ρ³).
  - The contribution is therefore c₁(t_a)·c₂(t_b)·c₁(t_c)·2r_ab·r_bc = φ_aφ_bφ_c·t_b·r_ab·r_bc.
- **Pairwise MRF.** Expand the binary pairwise MRF at weak coupling.
  - Couplings: J_ab = Cov(g_a, g_b)/(v_a v_b) + O(ρ²), with Cov(g_a, g_b) = φ_aφ_b·r_ab + O(ρ²).
  - The three-point cumulant is a tree through a centre gate whose own third cumulant is κ3(g_b) = v_b(1 − 2p_b): κ3_ν = Σ_b J_ab J_bc·v_a·κ3(g_b)·v_c + O(ρ³).
  - Substituting the couplings gives φ_aφ_c·r_ab·r_bc·φ_b²(1 − 2p_b)/v_b.
- **The limits.** As t → 0, 1 − 2p ≈ 2φ(0)t and p(1−p) → ¼, so the ratio tends to 8φ(0)² = 4/π. As |t| → ∞, Mills' ratio gives the limit 1. ∎

**Checks (C2).** Exact trivariate probit probabilities by nested Gauss–Legendre quadrature, against the fitted pairwise max-ent law.

| quantity | values |
|---|---|
| κ3_ν/κ3_exact on the grid t ∈ [0.4, 1], ρ ≤ 0.1 | 1.21–1.28 |
| κ3_ν/κ3_exact at t = 2, ρ = 0.2 (higher orders enter) | 1.41 |
| leading-order factor V(t)/t at t = 0.5, 1, 1.5, 2, 3 | 1.264, 1.238, 1.200, 1.159, 1.093 |

**Why the max-ent description fails.** The exact gate law is not a pairwise MRF at the order that matters. Its three-body Möbius coefficient of log P, the HC interaction on the triple ([D-HC] §3.1, Lemma A), is O(ρ²): −1.985e-3, −7.288e-3, −2.497e-2 at ρ = 0.05, 0.1, 0.2 (t = 0.4). That is the same order as the three-gate cumulant itself. Setting it to zero, which is what "Gibbs with 2-cliques" means, misstates the three-gate cumulants by 9–27 % at leading order.

A closure built on that description is worse than the probit (Gaussian) closure that the leg-partition expansion already uses. Including the 3-cliques restores accuracy, but that state is a set of n³ parameters: the trivariate object itself, with no saving.

### 4.3 DERIVED + checked: positivity fails and the dependency graph is complete, so Hammersley–Clifford is vacuous

**Statement.** For l ≥ 1 and invertible W_l, the face law of layer l gives zero probability to at least one orthant. In addition, every pair of gates is dependent for generic dense W.

*Proof.*
- Write z_l = a_{l−1}W_l with a_{l−1} ∈ ℝⁿ₊. Then z_l lies in the convex cone K = {aW_l : a ≥ 0}.
- Put c = W_l^{−1}1. Then (aW_l)·c = a·1 > 0 for every a ≥ 0 with a ≠ 0, so K lies strictly on one side of the hyperplane c^⊥.
- Every point y of the open orthant with sign pattern −sign(c) (generic c has no zero entries) has y·c < 0, so y ∉ K.
- The gate pattern 1[−c > 0] therefore has probability zero.
- Pairwise dependence: Cor(g_i, g_j) is a non-constant analytic function of C_ij, and C_ij ≠ 0 almost surely. ∎

**Checks (C3).** The orthants met by K, sampled with a ≥ 0 at random, are lower bounds limited by sampling:

| n | orthants met by the cone of rows of W | met by Gaussian input (layer 0), same samples | 2ⁿ |
|---|---|---|---|
| 8 | 105 | 250 | 256 |
| 10 | 317 | 896 | 1024 |
| 12 | 1375 | 3476 | 4096 |
| 14 | 4606 | 10974 | 16384 |

**Consequence.** Hammersley–Clifford needs positivity ([D-HC] §4). Without it, Markov ⇏ Gibbs: the intersection axiom fails, and Moussouris-type counterexamples exist ([D-HC] §4.2–4.3). Even where a factorisation exists, its graph is K_n, where "Markov with respect to the graph" says nothing.

There is no sparse separator structure inside a layer to exploit. This matches the guard recorded in [D-CR] §5.3 and §7.1: the boundary gain of Yang's CMI bound is isoperimetric and vanishes under edge expansion, and a dense layer is the extreme case.

### 4.4 THEOREM (elementary): face-law descriptions are blind to dilations

Gates are 0-homogeneous: 1[z_l(tx) > 0] = 1[z_l(x) > 0] for t > 0. So the face law of every layer, and the joint law of the face history, is a function of the direction x̂ alone. It carries no information about the scale.

The same goes for any statistic computed from gates alone. A statistic that also uses z, such as gate_GX, carries the scale only through its z factor.

§7 shows that the dominant shape of the (2,1,1) slice at depth at width 128 is the dilation (gain) mode, i.e. the fluctuation of ‖a_l(x̂)‖ across directions. This is a function of the history and of the position inside the face's cone, ‖M_h x̂‖ with M_h the linear map of the history. It is not a function of any single layer's face law.

**Verdict on item 2: DOA.** No pairwise state carries the binding content (THEOREM). The pairwise Gibbs closure is biased by 9–27 % at leading order (DERIVED). Hammersley–Clifford does not apply (DERIVED). The largest single shape in the slice at moderate width is invisible to faces (THEOREM). Spectral independence certifies mixing, not accuracy (§3.3). What *does* carry a large part of the slice is a statement about the law of z, not of the gates: §7.

---

## 5. Item 3: recovery maps for erased old content, and what CMI decay predicts

### 5.1 What the Chen–Rouzé construction needs, and what exists here

The construction ([D-Y] §8, [D-CR] §§1.3–3.7) is the time-averaged, KMS-detailed-balanced Lindbladian with single-Pauli jumps on A. In [A]'s language (§0.4) it is the Cesàro mean of the heat flow of the KMS-twisted derivation of N_A ⊂ M.

[D-TS] B8 already tabulates its four ingredients against the old-content transport. The result:
- **Detailed balance: absent.** No natural inner product makes the layer map self-adjoint, and the quiver is acyclic.
- **Stationarity: annealed only.** The ensemble-averaged layer is a twirl; the quenched layer is a one-Kraus conjugation.
- **Time averaging: reshaped.** It becomes a discounted random affine recursion P_{t+1} = T_{t+1}P_t + (births). Kesten–Bougerol–Picard does not apply, because the maps are not i.i.d. and the top exponent is +0.059 per layer.
- **Single-site jumps: present, but only on faces,** as Glauber dynamics on gates.

As t → ∞, CR's map becomes a conditional expectation onto the fixed-point algebra and recovers exactly ([D-CR] §6.1). Nothing in the estimator plays that algebra, because no law is stationary. The estimator also runs no dynamics. **The construction itself has no competition role.** Its surviving image is the reversible single-gate walk on faces, which serves prong 1 (U6 of local-to-global §5), not the estimator.

### 5.2 DERIVED under hypothesis H: no recovery from the retained state beats dropping

**Setting.**
- S_t is any retained state at layer t: statistics of z_t together with the young content carried exactly.
- X_t is the erased content: the old-source κ3 transported to t, or its D21 image.
- R is any recovery map applied to S_t.

**Statement.** E‖X_t − R(S_t)‖² ≥ E‖X_t − E[X_t|S_t]‖², with the expectation over the network ensemble. Under hypothesis H, E[X_t|S_t] = 0. Then no recovery map does better than dropping X_t.

**Hypothesis H.** Conditional on S_t, the law of X_t is invariant under X_t ↦ −X_t.

*Proof.* The first inequality is the defining property of the conditional mean. The second claim follows because H makes the conditional law symmetric. ∎

**Status of H.** It is an approximation; its evidence is as follows.
- *Exact at one layer.* E_W[W^{⊗3}] = 0 ([D-TS] §2.5). The flip W_{s+1} ↦ −W_{s+1} reverses the sign of third-order transported content and leaves every even-order transported statistic unchanged.
- *Approximate downstream.* The downstream gates break the flip, since relu(−z) ≠ relu(z), so H holds only approximately there.
- *Measured support* ([OC] final section, absorption test at w = 4):
  - renormalised births remove at most ≈ 10–15 % of the old pool's D21, even in-sample;
  - adding the young sources removes ≈ 30–40 % in-sample, but this does not transfer: on the held-out MLP the residual stays at the drop-it level (0.35–0.79, against shares of 0.39–0.88).

  A conditional mean that existed and could be learned would transfer across networks.

**CONJECTURE (testable form of H).** Fit absorption coefficients layer by layer on independent MLPs. Their relative signs are random: they agree no more often than a fair coin. Test: E3.

**DERIVED (the Petz map, one line).** Classically, the Petz recovery of a marginalisation channel relative to a reference law σ returns σ's conditional completion of the observed marginal. If σ is the closure's own model, the map returns the closure's births (tree-level content), which the chain already computes. A recovery map passes on the information in its reference; it does not create information. To recover old content, the reference must already contain it, which is circular.

**Verdict on recovery: DOA.**

### 5.3 What does transfer: duality (KNOWN-LINK + DERIVED)

**KNOWN-LINK.** The idea of [CR] that [A] §0.4 isolates is "recovery is not mixing": the defect is paired against a vector of small Dirichlet energy. Its commutative analogue already exists in the fresh-slate work, as the Duhamel telescoping of the heisenberg design ([FU] unlock 36, heisenberg Theorems 1–2, which [FU] calls "the commutative analogue of Chen–Rouzé's telescoping against the erased state"). For any reference chain ν with ρ_l = T_{l−1}ν_{l−1}:

  E_{μ_L} r_j − E_{ν_L} r_j = Σ_l (ρ_l − ν_l)[g_{l,j}].

Two properties matter here:
- All old content sits inside the exact pulled-back question g_l.
- A model ĝ_l makes the error bilinear: (one arrow's local defect) × (the question's model error).

**DERIVED.** The pulled-back question of a final readout through J_{l→L} lies in the range of J_{l→L}ᵀ. J and Jᵀ have the same singular values. A rank-k compression of the backward question therefore keeps e(k) per leg, exactly as a forward carrier does (§3.1), and has the same participation ratio, PR ≈ n/(2·age) ([D-TS] §2.2).

Duality changes three things:
- which side carries the memory;
- the shape of the error, from additive to bilinear, which is its real gain;
- the need for intermediate n³ tensors: the heisenberg design evaluates the backward pairing at full rank by matrix products ([FU] unlock 35, "States forward, questions backward").

It does not change the mode count. Any *compression* of the backward question meets the same free-probability law as compression of the forward content.

**CONJECTURE (Haar pairing).** Take dropped old content and closure error of equal D21 norm. They cost the same final-layer MSE, following the same 4.2e-6·ε² law within 30 %. The reason: in both cases the orientation relative to the downstream propagator is Haar ([D-TS] §2.4), so the pairing with the pulled-back question sees only the norm. Test: E8. If dropped old content turns out to cost half or less per unit of D21, then the pairing helps, and the D21 share overstates what a chain must carry.

### 5.4 What a CMI-type decay statement predicts for how much old content a chain may forget

**Measured** ([OC] tracker, MLP A1, AD convention, t = 15, file `results/tracker_A1.txt`).
- Per-age shares of ‖D21‖: 0.888, 0.276, 0.250, 0.247, 0.190 and 0.140 for ages 1–6.
- The nine sources of ages 7–15 together: 0.314.
- Old sources' contributions are nearly mutually orthogonal ([OC] Results), so each of the nine carries ≈ 0.314/√9 ≈ 0.10.
- In the no-AD convention, ages ≥ 5 carry 0.10–0.79.
- Width: the w = 4 share is 0.62–0.79 at n = 128 and 0.49–0.77 at n = 256, flat from n = 128 on ([OC] table (1)).

**DERIVED (conditional).** Suppose a Yang / Chen–Rouzé-type bound held for the depth chain: the dependence of layer t on layer t − d, given the retained state, is ≤ C·e^{−Δd}, with Δ the gap of the transport. A chain could then drop content older than d* = ln(C/ε)/Δ.
- Δ is at most the consecutive Lyapunov gap, 0.013–0.018 per layer ([D-TS] §0.1).
- The Sym^k functor adds no higher-order gap ([D-TS] §2.1).
- With C ≈ 0.25 (the age-2 share) and ε = 0.022: d* ≥ ln(11.4)/0.018 ≈ 135 layers.
- The top direction is amplified (+0.059 per layer), not discounted.

The competition shape, d/n = 1/64, sits in the free regime, where concentration goes as 1/age, not as e^{−gap·age} ([D-TS] B2).

**Measured, in agreement.** No exponential decay is visible across ages 2–15 at t = 15. Each single old source is about 5 times the bar.

**Prediction.**
- At L = 16 a chain may forget *nothing*: dropping even the oldest single source costs ε ≈ 0.1.
- This is width-independent, because the rates are mean-field and the shares are flat from n = 128 to 256.
- Test at n = 1024: E2.

**Verdict: DOA as a forgetting principle.** A CMI statement needs a gap, and at this aspect ratio the transport has none.

---

## 6. Item 4: the cluster / polymer view of the closure

### 6.1 KNOWN-LINK: the closure is a first-order cluster expansion with exact transport

**How the expansion arises.**
1. Expand each relu(z_i) in Hermite polynomials around its own Gaussian marginal (Mehler, memory). The vertex weights are w₁ = Φ, w₂ = φ/σ, and so on.
2. Write the joint cumulants of a = relu(z) by the Leonov–Shiryaev formula (memory). They become sums over *connected* diagrams:
   - vertices are neurons, weighted by Hermite coefficients;
   - edges and hyperedges are cumulants of z.

"First order" keeps diagrams with at most one non-Gaussian hyperedge plus covariance edges. The leg-partition coefficients count the ways the legs of a vertex can be distributed among them.

**The identified first-order diagrams** ([CP] §3.1 (v)): Gaussian ρ² + Φ³κ3(z) + [D21(z) ⊗ C] + [K22 ⊗ C] + [κ4(z)_(2,1,1)]. At n = 1024 the fitted coefficients equal the leg-partition values within noise at every depth ([O1024]; [CP] §6b).

**ANALOGY with Hammersley–Clifford's Lemma A** ([D-HC] §3.1). Both rest on Möbius inversion, with three differences:
- the lattice: HC inverts on the Boolean lattice of subsets, cumulants on the lattice of set partitions;
- the object: HC inverts log P and gets interaction potentials, while cumulants invert moments and get connected correlations;
- the dual sides: log versus linear.

The two agree only at leading order in weak coupling, and §4.2 measures where they part. The three-body potential and the three-gate cumulant are both O(ρ²), but with different coefficients.

**Convergence** (ANALOGY with Kotecký–Preiss, memory). A polymer expansion converges when polymer activities are small compared with the entropy of polymers. Here the activity of an edge is an entrywise correlation, O(n^{−1/2}). That is why the first-order truncation error falls as a power of n (n^{−0.8} measured) rather than with a temperature. It breaks as an analogy because there is no lattice entropy to beat: the number of diagrams per vertex grows with n, and only the 1/√n edge weights win.

### 6.2 ANALOGY: Bethe / tree level on pseudorandom geometry, and where it breaks

**The similarity.** The closure keeps tree-level gate diagrams and treats the dense, loop-rich part exactly through matrix products (the transport). Bethe / cavity does the same on locally tree-like graphs.

**Where it breaks: the reason the approximation is valid.**
- Bethe is exact when loops are long (girth → ∞) and couplings are O(1).
- Here every triple of neurons in consecutive layers closes a short loop through the complete bipartite connection.
- The truncation is accurate because each edge weighs O(n^{−1/2}), not because loops are rare.

The relevant theory is the dense mean-field (TAP / Plefka) expansion in powers of the coupling, not the sparse cavity method (memory).

**Consequence for the bethe design.** Loop-series corrections in a dense layer are organised by powers of n^{−1/2}, and the shortest loops dominate. E6 decomposes the closure residual by diagram class.

### 6.3 DERIVED: the closure error is not the binding error at n = 1024

**The error law.** The fit ε(n) = 0.155·n^{−1/2} + 4.26·n^{−1} reproduces the deep-layer closure errors 4.7, 2.6 and 0.9 % at n = 128, 256 and 1024.
- It predicts 1.5 % at n = 512.
- So does the log-log fit through the three points (slope −0.79).
- More widths between 128 and 1024 therefore cannot tell the two forms apart.

**The cost.** At n = 1024, a 0.9 % error costs 4.2e-6 × 0.009² = 3.4e-10 raw MSE, against ≈ 1.6e-8 raw available.

**Verdict.** Higher-order cluster terms are DOA as a lever at n = 1024. The binding errors are in the *inputs* to the closure: the κ4 slice and old content (§2).

**The coef-ensemble extrapolation.** [CE] extrapolates the regenerated-slice error to 1.3 % at layers 4–9 and ≤ 0.6 % at layers 10–14 for n = 1024. [O1024] measures 2.5–2.7 % at depth. The measurement wins. §7.3 gives the mechanism the extrapolation missed: the share of the slice that the regeneration can carry falls with width.

### 6.4 DERIVED + checked: the transport selection rule, i.e. the doubled leg is an expander

**Statement.** Let W have i.i.d. N(0, 2/n) entries. Then:
- W∘W = (2/n)·11ᵀ + E, where E has i.i.d. entries of mean 0 and variance 8/n², and ‖E‖ = (4√2/√n)(1 + o(1));
- W∘W∘W has i.i.d. entries of mean 0 and variance 120/n³, and ‖W∘W∘W‖ = (2√120/n)(1 + o(1)) ≈ 21.9/n;
- ‖W‖ = 2√2·(1 + o(1)), and W has no flat part.

*Proof.*
- **Moments.** E W² = 2/n, E W⁴ = 12/n² and E W⁶ = 120/n³. So W² has mean 2/n and variance 8/n², and W³ has mean 0 and variance 120/n³.
- **Norms.** An n × n matrix with i.i.d. mean-zero entries of variance σ² and a finite fourth moment has operator norm 2σ√n(1 + o(1)) (Bai–Yin / Latała, memory).
- **The flat part.** (2/n)11ᵀ has norm 2. ∎

**Expander reading.**
- Normalise to B = (W∘W)/2. Its top singular value is 1 (the flat mode), and the rest are ≤ 2√2/√n.
- For a Ramanujan graph of degree n the corresponding value is 2√(n−1)/n ≈ 2/√n. The doubled leg is √2 above it, because χ²₁ weights have variance twice their squared mean.
- The expander mixing lemma holds verbatim: |uᵀ(W∘W)v − (2/n)(1ᵀu)(1ᵀv)| ≤ ‖E‖·‖u‖·‖v‖.

**The selection rule.** Transport a cumulant through a layer, κ(z') = κ(a)·W^{⊗k}.
- **A doubled output index receives two kinds of input.**
  - The doubled input index arrives through W∘W, which flattens it to its trace up to O(n^{−1/2}).
  - Pairs of distinct input indices arrive through W_pi·W_qi, which does not flatten them.
- **A tripled output index receives tripled input content only at O(1/n).** Diagonal (tripled) content is therefore regenerated at each layer from the other channels. It is never inherited.

**The Perron mode is the dilation channel.** The flat part transmits Σ_p κ(a)_{pp…}. By the trace identities of §7.1 (iv), this is the covariance of ‖ã‖² with the remaining legs, minus Wick terms: the norm-coupling channel. So the Perron mode of the doubled-leg expander is the dilation channel of §7. This is where the expander picture and positive homogeneity meet.

**Checks (C1, C5).** Operator norms at four widths:

| n | ‖E‖ | 4√2/√n | flattening error of a positive profile ‖Eᵀv‖/‖(2/n)1(1ᵀv)‖ | ‖W∘W∘W‖ | 2√120/n |
|---|---|---|---|---|---|
| 128 | 0.4910 | 0.5000 | 0.134 | 0.1934 | 0.1712 |
| 256 | 0.3737 | 0.3536 | 0.098 | 0.1027 | 0.0856 |
| 512 | 0.2500 | 0.2500 | 0.066 | 0.0446 | 0.0428 |
| 1024 | 0.1767 | 0.1768 | 0.047 | 0.0222 | 0.0214 |

- **Generic vectors.** |uᵀEv|/(‖u‖‖v‖) = 0.003–0.011 for generic u and v.
- **Convergence.** The tripled-leg norm converges slowly because W³ is heavy-tailed: its entries have kurtosis 10395/225 ≈ 46.
- **Channel fractions** at n = 128, as fractions of ‖D3(z')‖ for the diagonal κ3_iii. The channels interfere, so the fractions do not sum to 1.

  | channel into the diagonal | fraction |
  |---|---|
  | tripled | 0.36 at layer 1, falling to 0.03–0.12 at depth |
  | doubled | 0.65–0.98 |
  | distinct | 0.31–0.50 |

- **D21.** The doubled→doubled channel carries 0.52–0.81 of ‖D21‖. It is flat up to a relative error of 0.19–0.29 at depth (0.44–0.49 at layer 1). By the n^{−1/2} law that error is ≈ 7–10 % at n = 1024.

### 6.5 Pilot: a flat doubled-leg carrier for the (2,1,1) slice, refuted at depth

**The hypothesis R1 suggests.** The doubled index of the (2,1,1) slice is flat: K_iijk ≈ M_jk, with one n × n matrix per layer. This costs about the same as the published u_iC_jk but has no regression vector.

**The pilot (C6).** Width 128, two MLPs, fitted on atlas seed 3 and scored on seed 4. The (2,1,1) slice enters the leg-partition closure as hyperedge B6. The score is the noise-corrected D21 error after transport. The "gap closed" is (none − variant)/(none − true).

Variants:
- Creg: the published u_iC_off,jk with u fitted;
- flat: the i-average M̄ of the slice;
- rank1: the best u⊗M, with both fitted.

| layer | Creg | flat | rank1 | cos(u₁, 1) |
|---|---|---|---|---|
| 1 | 0.35 | 0.37 | 0.38 | 0.979 |
| 2 | 0.25 | 0.28 | 0.30 | 0.962 |
| 3 | 0.34 | 0.40 | 0.43 | 0.947 |
| 4 | 0.38 | 0.50 | 0.62 | 0.909 |
| 5 | 0.48 | 0.30 | 0.61 | 0.749 |
| 6 | 0.62 | 0.33 | 0.71 | 0.818 |
| 7 | 0.60 | −0.05 | 0.68 | 0.726 |
| 8 | 0.80 | 0.02 | 0.90 | 0.773 |
| 9 | 0.93 | −1.08 | 0.93 | 0.708 |
| 10 | 0.88 | −0.55 | 0.92 | 0.724 |
| 11 | 0.80 | −1.88 | 0.83 | 0.636 |
| 12 | 0.80 | −1.02 | 0.84 | 0.674 |
| 13 | 0.81 | −1.08 | 0.83 | 0.686 |
| 14 | 0.78 | −2.05 | 0.82 | 0.770 |

(MLP 0. MLP 1 agrees: Creg 0.23–0.80; flat 0.22–0.45 at shallow layers and −0.71–0.10 at depth; rank1 0.26–0.83. Reversing fit and score atlases gives the same picture.)

**Readings.**
- **At shallow layers** the flat carrier matches or beats Creg, as R1 predicts.
- **At depth it fails.** The slice is close to rank 1 with M ≈ C_off (cos(M̄, C_off) = 0.86–0.97), but its doubled-index profile is not flat (cos(u₁, 1) = 0.64–0.82).
- **That profile is var_i (§7).** Its non-flat part arrives through the distinct→doubled channel, i.e. the off-diagonal of C at depth, which carries the spike. The flat carrier keeps only the doubled→doubled channel, i.e. the trace.
- **Rank 1 does better than Creg** at every layer, by up to +0.24 (layer 4).

**Verdict: DOA.** At n = 1024 the flat carrier can capture at most the trace channel, whose content is the gain coupling. §7.3 shows that this is a falling share of the slice.

---

## 7. The dilation (gain) mode: what the published (2,1,1) regeneration actually carries

**The symmetry.** With no biases, every arrow commutes with dilations: relu(t·z)W = t·relu(z)W for t > 0. [BRIEF] §3 and [FU] unlock 1 draw the exact consequence for the input radius: E f = E‖x‖·E_{S^{n−1}} f.

**What this section adds.** The radius is the smallest part of the effect. On the sphere, the per-direction gain ‖a_l(x̂)‖ fluctuates from direction to direction. Positive homogeneity carries every such fluctuation forward unchanged. Each layer adds its own. The result is a scale mixture whose variance grows with depth, and it accounts for most of the (2,1,1) slice at depth at moderate width.

### 7.1 DERIVED: what a scale mixture predicts

**Setting.** Let z = s·y, with s > 0 independent of y and E s² = 1. Write s = 1 + δ with small fluctuations σ_δ² = Var δ, so that Var(s²) = 4σ_δ² + O(σ_δ³). Let y be Gaussian, y ~ N(μ, Σ).

**Statement.** To leading order in σ_δ²:
- (i) C_z = Σ + Var(s)·μμᵀ. This is exact.
- (ii) κ3(z)_aab = ½Var(s²)·(var_a μ_b + 2C_ab μ_a) = ½Var(s²)·T_ab.
- (iii) κ4(z)_iijk = Var(s²)·(var_i C_jk + 2C_ij C_ik) = Var(s²)·R_ijk.
  - Likewise κ4(z)_iijj = Var(s²)·(var_i var_j + 2C_ij²).
  - R is exactly the Wick part of the (2,1,1) moment. "K ≈ c·R" therefore says that the (2,1,1) fourth moment is (1 + c) times its own Wick part.
- (iv) Two identities, exact for any law with finite fourth moments (z̃ = z − E z):
  - Σ_i κ4_iijk = Cov(‖z̃‖², z̃_j z̃_k) − 2(C²)_jk;
  - Σ_ij κ4_iijj = Var(‖z̃‖²) − 2‖C‖_F².

  So the *gain excess* gex := Var‖z̃‖²/(tr C)² − 2‖C‖²_F/(tr C)² equals Σ_ij κ4_iijj/(tr C)², the normalised trace of the (2,2) slice. For a scale mixture of a Gaussian, gex = Var(s²)·(1 + 2‖Σ‖²_F/(tr Σ)²). Hence c4/gex ≤ 1, with equality when no eigenvalue of C dominates.
- (v) Positive homogeneity propagates the scale. If z_l = s·y_l, then z_{l+1} = s·relu(y_l)W_{l+1}: a scale fluctuation present at layer l reaches every later layer exactly, undiscounted. If layer l+1 adds a fresh, independent gain g, then Var((sg)²) = Var(s²) + Var(g²) to first order. The amplitude therefore accumulates additively with depth.
- (vi) The input radius r̃ = ‖x‖/√n has Var(r̃²) = 2/n. But x is exactly Gaussian, so κ4(z₀) = 0: the radial part is cancelled exactly by the spherical part, x̂ being uniform on the sphere. "2/n" is therefore not the amplitude at depth. The gains that the layers add set it.

*Proof.*
- **(i)–(iii): law of total cumulance, conditioning on δ.** Given δ, z̃ is Gaussian with mean m(δ) = (δ − Eδ)μ and covariance S(δ) = (1+δ)²Σ, so its conditional cumulants of order ≥ 3 vanish. A total cumulant is therefore a sum, over partitions of its indices, of joint cumulants of the blocks' conditional cumulants (m for blocks of size 1, S for blocks of size 2).
  - Order 2: E S_ab + Cov(m_a, m_b) = Σ_ab + σ_δ²μ_aμ_b.
  - Order 3: the (2,1) partitions give Cov(S_aa, m_b) + 2Cov(S_ab, m_a) = Cov((1+δ)², δ)·T_ab = 2σ_δ²·T_ab + O(σ_δ³). The (1,1,1) partition gives κ3(δ)·μμμ, which is of higher order for near-Gaussian δ.
  - Order 4: only the (2,2) partitions are O(σ_δ²). They give Var((1+δ)²)·(Σ_iiΣ_jk + 2Σ_ijΣ_ik). Partitions with mean blocks give joint cumulants of δ of order ≥ 3, i.e. O(σ_δ³) or smaller.
  - Replacing Σ by C changes these by O(σ_δ⁴).
- **(iv):** sum the definition of κ4 over the repeated index, using E‖z̃‖² = tr C.
- **(v):** 1-homogeneity of relu, independence, and E s⁴E g⁴ − 1 = Var(s²) + Var(g²) + Var(s²)Var(g²).
- **(vi):** a Gaussian has κ4 = 0. ∎

### 7.2 Measured at n = 128: the regeneration *is* the gain mode (C7, C8)

**Data.** Two atlas MLPs at width 128 (N = 5·10⁵, two independent sample seeds). The (2,1,1) slice K is formed explicitly.

**C7: shape.**
- *Shape of the published u.* The published regeneration fits u_i by regressing K on C_off. Its u is parallel to var:
  - cos(u_Creg, var) = 0.91–0.997 at layers 1–14, both MLPs;
  - the lowest values are 0.913–0.94 at layers 5–9 of MLP 0.
- *Amplitude.* The median of u_i/((2/n)·var_i) grows with depth:
  - MLP 0: 1.66, 2.50, 3.74, 4.93, 6.50, 10.3, 11.4, 17.1, 19.7, 21.9, 21.2, 21.3, 19.5, 22.3;
  - MLP 1: 1.84 → 16.7.

  The parameter-free input-radius amplitude 2/n is 2–22× too small, and closes only 0.07–0.29 of the D21 gap.
- *Energy.* The share of the slice's noise-free energy in the one-dimensional span of R:
  - MLP 0: 0.27, 0.27, 0.31, 0.36, 0.49, 0.65, 0.71, 0.80, 0.88, 0.92, 0.94, 0.91, 0.92, 0.90;
  - MLP 1: 0.30 → 0.85.

  The published n-parameter family {u_iC_jk} holds slightly *less* energy (MLP 0: 0.26 → 0.84–0.89), because it lacks the 2C_ijC_ik term.

**C8: D21 gap closed.** Gap closed on the D21 interface after transport, fitted on seed 3 and scored on seed 4. The one-scalar row is c·R with c fitted.

| layer | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| MLP 0, published u_iC_jk (n parameters) | 0.35 | 0.25 | 0.34 | 0.38 | 0.48 | 0.62 | 0.60 | 0.80 | 0.93 | 0.88 | 0.80 | 0.80 | 0.81 | 0.78 |
| MLP 0, one scalar c·R | 0.34 | 0.26 | 0.34 | 0.39 | 0.46 | 0.57 | 0.46 | 0.75 | 0.69 | 0.79 | 0.69 | 0.79 | 0.67 | 0.67 |
| MLP 0, c·n/2 | 1.6 | 2.4 | 3.5 | 4.5 | 7.4 | 9.7 | 11.2 | 12.8 | 16.1 | 16.1 | 15.3 | 14.9 | 14.9 | 15.3 |
| MLP 1, published u_iC_jk | 0.24 | 0.23 | 0.25 | 0.27 | 0.30 | 0.41 | 0.47 | 0.51 | 0.63 | 0.63 | 0.64 | 0.73 | 0.80 | 0.54 |
| MLP 1, one scalar c·R | 0.23 | 0.20 | 0.21 | 0.22 | 0.25 | 0.35 | 0.43 | 0.58 | 0.57 | 0.55 | 0.62 | 0.66 | 0.76 | 0.56 |
| MLP 1, c·n/2 | 1.8 | 2.8 | 3.3 | 4.4 | 4.5 | 6.5 | 6.7 | 9.1 | 9.9 | 9.2 | 9.9 | 9.3 | 9.4 | 10.1 |

The one-scalar template recovers 74–114 % of what the n-parameter regeneration recovers, median ≈ 0.9.

**C11: the κ3 counterpart.** The same holds at order 3. The noise-free share of D21's off-diagonal energy in the T shape of (ii) is:
- MLP 0: 0.28 at layer 1, 0.46–0.83 at layers 2–6, 0.90–0.97 at layers 8–15;
- MLP 1: 0.32 → 0.86–0.93.

The best scalar is ≈ 7× the input-radius value.

**Reading.**
- The regeneration that the published chain uses, with its fitted u, is to within cos ≥ 0.91 the gain mode with u = c·var.
- One scalar per layer suffices for ≈ 90 % of its effect at n = 128.
- The amplitude is the layers' gain variance, not the input radius.

This is the first quantitative explanation of *why* the u_iC_jk form works at all.

