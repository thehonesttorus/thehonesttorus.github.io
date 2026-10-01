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

