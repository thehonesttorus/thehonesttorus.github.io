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
- **DERIVED**: proved here, with the proof in the text. "+ checked" means also verified numerically (checks C1–C16, §10).
- **ANALOGY**: a structural similarity, with the exact point where it breaks.
- **CONJECTURE**: a precise statement, not proved, with its test.
- **SPECULATION**.

Measured numbers, ours or the streams', are evidence attached to a labelled claim. They are not claims of their own.

**Sigla.**
- Digests: [D-Y] Yang; [D-CR] Chen–Rouzé; [D-HC] Hammersley–Clifford; [D-EXP] expanders; [D-HDX] HDX / spectral independence; [D-NC] noncommutative Dirichlet forms; [D-TS] transfer-spectrum measurement; [D-CHAT] chat retrieval status.
- Sibling notes: [A], [B].
- Stream reports: [O1024] oracle1024; [OC] old-content; [C128] chain128; [CE] coef-ensemble; [CM] costmodel; [EA] est-accuracy (stopped).
- Plans: [CP] competition plan; [BRIEF] fresh-slate brief; [FU] fresh-slate `foundations-unlocks.md`, consulted only for cross-references (unlocks 1, 35 and 36).

**Retrieval.** No new source was opened. Published results are cited through the digests that quote them. Standard facts used from memory are marked (memory):
- the Bai–Yin / Latała operator-norm law for i.i.d. matrices;
- Mehler's formula and the Hermite expansion of an indicator;
- the Leonov–Shiryaev formula for cumulants of products;
- the Lebowitz–Percus–Verlet ensemble corrections;
- the Roberts–Yaida–Hanin finite-width four-point vertex;
- TAP / Plefka expansions;
- the Kotecký–Preiss criterion.

Each is used only as a name, or is re-derived or checked here. The ChatGPT conversation behind [D-Y] could not be retrieved ([D-CHAT] §1).

**New measurements for this note.** The checks C1–C12 run at widths 8–1024. They took about three hours of wall time on a shared 4-core box. The scripts are scratch files, not committed, because this task writes one file. §10 and the Appendix give every formula needed to reproduce them.

**Review.** A mathematical-correctness review (critic 1, same day) re-derived every DERIVED claim, re-ran the cheap checks and added C13–C16 (§10). Its corrections are folded into the text; the list, with what changed and why, is the review log at the end.

---

## 0. The answer in brief

1. **Where the "same essence" is literally a theorem in the competition object: the doubled leg.** When a cumulant tensor is pushed through a linear layer, the transport depends on how often an index repeats.
   - An index that appears once is transported by W.
   - An index that appears twice is transported by the Hadamard square W∘W.
   - An index that appears three times is transported by W∘W∘W.

   For He weights, W∘W = (2/n)·11ᵀ + E with ‖E‖ = (4√2/√n)(1 + o(1)). This is a dense weighted expander, √2 above the Ramanujan value for degree n, and the expander mixing lemma holds verbatim. W∘W∘W has norm O(1/n). That does not make the tripled channel negligible, because its target, the diagonal D3, is itself O(1/n) per unit: at layer 1 the tripled channel carries ≈ 0.6 of ‖D3(z₁)‖ at every width, falling to 0.02–0.10 at depth (C15, §6.4). W itself (norm 2√2, no outlier singular value) scrambles orientations. DERIVED + checked (C1: ‖E‖ = 0.1767 at n = 1024, predicted 0.1768). The flat (Perron) mode of the doubled leg transmits exactly the *trace channel* of a cumulant, Σ_p κ_{pp…}. That channel is the covariance of the squared norm with the remaining legs. §6.4.
2. **Where it is two halves of one mechanism.** Heredity is free and decorrelation does the work.
   - The layer chain is exactly Markov, because each layer is a deterministic function of the previous one.
   - The closure works because pre-activation correlations are entrywise O(n^{−1/2}) and single legs scramble orientation (Haar, [D-TS] §2.4).
   - The leg-partition closure is a first-order cluster expansion with exact transport. The diagram expansion is a KNOWN-LINK (Leonov–Shiryaev); the identification of the closure with its first-order truncation is the programme's ([CP] §3.1 (v)), and one fitted coefficient (B3) is still off: 0.15–0.40 against a counted 1 (§6.1). Its error falls as n^{−0.8} to n^{−0.9} and is 0.74–1.07 % at n = 1024, so it is not the binding error. §6.
3. **Where it is false.**
   - **(a) A pairwise Gibbs / MRF description of the face law cannot carry the binding content.**
     - Pairwise gate statistics are bivariate functionals (THEOREM).
     - The 2-clique Gibbs closure misstates three-gate cumulants. For equal thresholds it overstates them by a factor in (1, 4/π]; for thresholds of mixed sign the error is unbounded and can flip the sign (DERIVED + checked, C2, C13).
     - Positivity fails at every layer ≥ 1 (DERIVED + checked). At small n the support also fails the weaker safe-symbol closure, which is all the Hammersley–Clifford proof uses, and the conditional-dependence graph is complete or nearly so (measured, C14).
     - Gates are 0-homogeneous, so a face-law description carries no information on the input radius (THEOREM). The per-direction gain of item 4 varies inside each face's cone, so a gate-only state cannot hold its law, only its conditional mean given the faces (elementary).
   - **(b) Spectral independence is a frame bound,** Cov ≼ (1+η)·diag. It is not an accuracy certificate for a closure (THEOREM + DERIVED).
   - **(c) A time-averaged, detailed-balance recovery map cannot restore erased old content.** The ingredients are missing ([D-TS] B8). What does restore part of it is a *symmetry*. In the scale-mixture model, positive homogeneity carries the scale's κ3 signature unchanged, so the erased content's component along the gain template is a function of the retained state plus one accumulated scalar (§5.2 (a), §7.5: model + measured; the measured coefficient is discounted by a few per cent per layer, so the model is approximate). The rest has zero Bayes-optimal recovery under a sign hypothesis H′ (DERIVED under H′; the evidence for H′ is weak, §5.2 (b)).
   - **(d) A CMI-type decay statement lets the chain forget only the oldest one or two sources.** Old content does decay, at ≈ 0.15–0.3 per layer of age. The rate is the gates' discount of bulk directions, not a spectral gap, which would be 10–20× slower. The forgetting horizon is therefore ≈ 8–16 layers, comparable to L = 16. Measured at n = 128 (one MLP, t = 15): sources of age ≥ 14 together carry 0.016 of ‖D21‖ (AD convention), below the bar, but those of age ≥ 7 carry 0.32 (AD) and 0.62 (no-AD). §5.4.
4. **The one bridge with a positive competition payoff: the dilation (gain) mode.**
   - **Mechanism.** With no biases, every arrow commutes with x ↦ tx for t > 0. The input radius is the obvious consequence ([BRIEF] §3). The less obvious one is that the gain of each input direction, ‖a_l(x̂)‖, fluctuates across directions, and that this fluctuation propagates undiscounted (exactly in the scale-mixture model of §7.1; the measured discount is a few per cent per layer, §7.5).
   - **What a scale mixture predicts** (z = s·y; DERIVED, §7.1):
     - κ4 (2,1,1) slice = Var(s²)·(var_i C_jk + 2C_ijC_ik), the Wick part of the slice's own moment;
     - κ3 (2,1) slice = ½Var(s²)·(var_a μ_b + 2μ_a C_ab);
     - a spike Var(s)·μμᵀ in C.
   - **Measured (C7–C9).**
     - The published regeneration u_iC_jk has u parallel to var (cos 0.91–0.997).
     - A single fitted scalar per layer recovers 74–114 % of the D21 gain that the n-parameter regeneration gets at n = 128.
     - At layer 15 the variance excess along the mean is 1.05–1.6× the gain-mode prediction, so the spike at depth is largely the gain mode at order 2.
   - **The amplitude.** It is width-stable as c·n. The scalar Var‖z̃‖²/(E‖z̃‖²)² − 2‖C‖²_F/(tr C)² is exactly the trace of the (2,2) slice divided by (tr C)² (DERIVED identity). In the scale-mixture model it bounds the amplitude from above, c4 ≤ gex. Outside the model there is no such inequality (§7.1 (iv)). Measured: c4/gex = 0.77–1.01 at n = 1024 and 0.5–1.0 at n = 128.
   - **The limit at n = 1024.** Across widths the template's share of the slice energy falls: at layers 6–15 it is 0.47–0.91 at n = 128 and 0.25–0.39 at n = 1024 (two networks each, C9). This matches, and plausibly explains, the published regeneration's falling share of the κ4 gap (55 → 33 → 20 % at n = 128 → 256 → 1024, [O1024], measured on a different network). At n = 1024 the binding (2,1,1) content is mostly *not* the gain mode.
   - **Not a Gaussian quantity.** A Gaussian-input recursion recovers layer 1 exactly but only 24–47 % of the amplitude at layer 15 (C10).
   - **Old κ3 content is mostly the gain mode, in the published convention.** In the no-AD convention the old pool's D21 is 49–97 % along the κ3 template T (87–97 % at layers 10–15; 77–94 % on a third MLP). A scalar accumulator, fed once by each source as it ages out of the young window and never discounted, removes most of its energy at zero cost (n = 128, three MLPs, C12 and C16):
     - residual 0.11–0.34 of ‖D21‖, against 0.20–0.88 if the pool is dropped;
     - still 5–15× the bar. At layers 10–15 the part that needs a real carrier shrinks 2.3–3.1× in D21 norm (2.8–5.4× with the best-fitting scalar, which a chain does not have).

     This had not been tested; [OC]'s absorption features were births and young sources. §7.
5. **The three measured facts named in the brief, priced.** §3.
   - *Free-probability law of old content.* An excellent offline calculator and a DOA carrier at n = 1024: 54–86 u/layer against 3–4 available. Content with more scrambled legs is worse (DERIVED: the kept energy goes as the m-th power of the per-leg fraction). The (2,1,1) slice plausibly transports like a covariance, with two scrambled legs (heuristic: its doubled index is flat at shallow layers but carries the profile var_i at depth, §6.5; test E4).
   - *BBP spike.* Alive only as deflation, with ≈ 10–25 % fewer modes (more at greater age). As a standalone carrier ε_tot ≥ 0.16 against 0.022.
   - *Spectral independence.* A diagnostic. Its growth with depth is in part the spike seen through the gates: the top mode overlaps the spike image at 0.39–0.49 at layers 1–4, against 0.01 for a random direction.
6. **NCG's part for the competition.** This is not its part for the theory, which [A] and [B] treat.
   - *Literal.* Free probability (Voiculescu's noncommutative probability) is the law of the propagators. The annealed layer is a conditional expectation, a twirl.
   - *Duality.* This is the one Chen–Rouzé ingredient that transfers. The heisenberg design already uses it as a Duhamel telescoping ([FU] unlock 36); that this is CR's telescoping in commutative form is an ANALOGY, and the telescoping identity itself is elementary. The linear part of its backward questions obeys the same free-probability mode law, because J and Jᵀ share singular values (DERIVED).
   - *Speculation* (§7.7). For positively homogeneous dictionaries the dilation may be read as a central charge of the arrow algebra, and the gain mode as its canonical-ensemble fluctuation. The Lebowitz–Percus–Verlet comparison behind this is an ANALOGY that breaks because the scale is not conserved.
   - *No role.* The detailed-balanced single-Pauli Lindbladian itself has no competition role: the estimator runs no dynamics and has no stationary state. §8.
7. **What to run** (§9). Ten experiments use existing tools. The decisive cheap ones are:
   - E1: rank-1 power iteration of the (2,1,1) slice at n = 1024, minutes per layer. Is the non-gain remainder one mode?
   - E7: the gain template at n = 1024 on more networks, two done here (C9).
   - E9: the old pool minus the template at n = 256, then the old-content carriers on the residual. Does §7.5 bring the old tier into the budget?

---

## 1. Results at a glance

| # | statement | label | for the competition | § |
|---|---|---|---|---|
| R1 | Doubled leg W∘W = (2/n)11ᵀ + E with ‖E‖ = 4√2/√n; tripled leg O(1/n) in norm, yet ≈ 0.6 of the diagonal target at layer 1, because that target is O(1/n) per unit too; a single leg scrambles. The flat mode carries the trace (norm-coupling) channel | DERIVED + checked (C1, C5, C15) | structural fact for every design | 6.4 |
| R2 | Pairwise gate statistics (gate_GG, gate_GX, any pairwise MRF) are bivariate functionals and cannot hold trivariate content as state | THEOREM | DOA (pairwise face-law state) | 4.1 |
| R3 | The 2-clique Gibbs closure of a probit gate law misstates three-gate cumulants: at leading order by the factor φ(t)(1−2p)/(t·p(1−p)) ∈ (1, 4/π] for equal thresholds; without bound, and possibly in sign, for thresholds of mixed sign | DERIVED + checked (C2, C13) | DOA (MRF closure) | 4.2 |
| R4 | At layers ≥ 1 the face law vanishes on at least one orthant (the cone of W's rows is pointed), so positivity fails. At small n the support also has no safe symbol, and the conditional-dependence graph is complete or nearly so | positivity: DERIVED + checked (C3); safe symbol and completeness: measured (C14, n = 6, 8) | DOA (HC factorisation) | 4.3 |
| R5 | Gates are 0-homogeneous, so every face-law description is blind to the input radius. The per-direction gain varies inside each face's cone, so a gate-only state cannot hold its law | THEOREM (elementary) | limits every gate-only design | 4.4 |
| R6 | η-SI ⟺ Cov ≼ (1+η)·diag; independent of the entrywise smallness the closure needs | THEOREM + DERIVED (C4) | DOA as certificate or lever; diagnostic | 3.3 |
| R7 | Propagator-basis carriers need ∝ n modes; content with m scrambled legs keeps e(k)^m | THEOREM (cited) + DERIVED | DOA at n = 1024; calculator alive | 3.1 |
| R8 | Spike alone: ε_tot ≥ 0.16 for the w = 4 tier; as deflation ≈ 10–25 % fewer modes | DERIVED (arithmetic on measured inputs) | DOA alone; small lever | 3.2 |
| R9 | Erased old content: no detailed-balance recovery (ingredients missing). Its gain-template part is recoverable by the dilation symmetry (accumulator); the rest has Bayes-optimal recovery 0 under H′ | ingredients: cited ([D-TS] B8); T-part: scale-mixture model + measured (C12, C16); the rest: DERIVED under H′, whose evidence is weak | DOA for detailed balance; alive for the symmetry | 5.1–5.2 |
| R10 | Duality (the commutative analogue of CR's telescoping, [FU] 36) moves the memory from the forward state to the backward question, whose low-rank compressions obey the same free law (for the question's linear part) | ANALOGY (CR ↔ Duhamel) + THEOREM (elementary telescoping identity) + DERIVED (mode law) | alive for its bilinear error, no mode saving | 5.3 |
| R11 | Old content decays at ≈ 0.15–0.3 per layer of age, set by the gates' discount, not by a spectral gap. The forgetting horizon (≈ 8–16 layers) is comparable to L, so only the 1–2 oldest sources fall below the bar | measured (C16, n = 128, one MLP) + DERIVED (conditional arithmetic) | DOA (forgetting by depth) | 5.4 |
| R12 | Closure = first-order cluster expansion with exact transport; error ∝ n^{−0.8} to n^{−0.9}, 0.74–1.07 % at n = 1024; fitted coefficients within 0.15 of the counted ones except B3 | KNOWN-LINK (Leonov–Shiryaev) + identification ([CP] §3.1) + measured | not binding | 6.1–6.3 |
| R13 | Flat doubled-leg carrier for the (2,1,1) slice | refuted at depth (C6) | DOA | 6.5 |
| R14 | Scale mixture ⇒ (to leading order in Var s²) κ4(2,1,1) = Var(s²)·Wick, κ3(2,1) = ½Var(s²)·T; C ∋ Var(s)μμᵀ exactly. For any law, gex = Σ_ij κ4_iijj/(tr C)² exactly; c4 ≤ gex holds only inside the model | DERIVED | explains the published regeneration | 7.1 |
| R15 | u_Creg ∥ var (cos 0.91–0.997); one scalar recovers 74–114 % of the regeneration's D21 gain at n = 128; c·n width-stable; template share falls with n | measured (C7–C9) | one scalar per layer; not enough at n = 1024 | 7.2–7.3 |
| R16 | The gain amplitude is not a Gaussian-closure quantity: Gaussian-input increments give 24–47 % of it at layer 15 | measured (C10); CONJECTURE refuted in its simplest form | the scalar needs its own closure | 7.4 |
| R17 | In the published (no-AD) convention the old κ3 pool is 49–97 % along T; an undiscounted scalar accumulator leaves 0.11–0.34 of ‖D21‖ (0.20–0.88 if dropped) | measured (C12, C16; n = 128, three MLPs); CONJECTURE at n = 1024 | old content becomes a residual 2.3–3.1× smaller at layers 10–15 | 7.5 |

---

## 2. The object, the bar and the budget: what a bridge must beat

- **Task** ([BRIEF] §1, [CP] §0).
  - Input x ~ N(0, I_n), n = 1024, L = 16 layers, W_l i.i.d. N(0, 2/n), no biases.
  - z_l = a_{l−1}W_l (x@W convention, a_{−1} = x), a_l = relu(z_l).
  - Estimate the n means E[a_15].
  - Score = final-layer MSE × max(0.1, F/B), with B = 2^41 FLOPs. The unit is u = 2^31 FLOPs, one dense 1024³ product.
- **The bar** ([CP] §0; [phase2-intel](../../digests/phase2-intel-2026-10-01.md)).
  - Money line: adjusted MSE ≤ 1.6–1.8e-9, i.e. raw ≤ 1.6e-8 at the 0.1 floor.
  - Best raw on the board: 1.14e-8 at 0.151 B.
  - Public V29: raw 2.13e-8 at 0.253 B. The same chain without its old-source tier runs at 0.150 B but raw 1.3–2.9e-7, ten times worse.
- **Interface law** ([CP] §3.1, [C128], [O1024]). In a cumulant chain the error that reaches the final layer through the nonlinearity's (2,1) interface obeys extra MSE ≈ 4.2e-6·ε², where ε is the relative rms error of D21_{ab} = κ3(z_a, z_a, z_b). The frontier needs ε ≤ 2.2 %.
- **Measured at n = 1024** ([O1024], two MLPs, layers 1–14; [CP] §6b), D21 error by κ4 input:

  | κ4 input to the closure | D21 error |
  |---|---|
  | exact inputs (leg-partition closure) | 0.74–1.07 % (0.74–0.91 % at layers 7–14) |
  | (2,1,1) slice regenerated as u_iC_jk | 2.14–2.69 % |
  | no slice | 2.4–3.2 % |

  The regeneration closes a falling share of the κ4 gap: 55 → 33 → 20 % at n = 128 → 256 → 1024.
- **End to end at width 128** ([C128]; [CP] §6b). With κ4 at truth, the D21 ladder maps one-to-one onto the final layer. With a κ4 closure of its own, the κ3 rule stops mattering: with κ4 = 0 every κ3 rule sits at 2.0–2.4e-5, and with the memoryless κ4 at 1.4–1.7e-4 (8.0–8.4e-5 at Edgeworth order 2). A dense n⁴ κ4 closure comes within ≈ 3× of the teacher-forced level ([C128] TL;DR 5; [CP] §6b summarises the range as 3e-5–8e-5).
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
- KNOWN-LINK ([D-TS] §2.2–2.3: Pennington–Schoenholz–Ganguli; Hanin–Nica): gated products obey free multiplicative convolution, and PR(J_{s→t}) ≈ n/(1 + Σ(r_l − 1)) ≈ n/(2·age). Measured, this holds to within 4 % up to age 4 (6 % at age 5). At greater ages the measured PR falls below the free value, by up to 31 % at ages 9–13 and 23 % at age 15 ([D-TS] §3.1), which the digest attributes to the mean-direction outlier.
- Measured ([D-TS] §3.4): k_2%/n is width-independent to ±0.02 from n = 128 to 1024. At n = 1024, age 8 needs 322 modes and age 15 needs 197 modes for a 2 % error relative to the content's own D21 contribution.

**DERIVED (arithmetic; [D-TS] §0.3).** A Tucker core of rank k costs k³/n² u per layer. The w = 4 tier needs ≈ 384–448 modes at n = 1024, i.e. 54–86 u per layer, against 3–4 available. **DOA as a carrier.**

**DERIVED (legs, under the hypothesis of [D-TS] §2.4).** Let a source be Haar-oriented relative to the right singular vectors of J, and let it have m *scrambled* legs (indices that appear once). Project each such leg onto the top-k left singular subspace of J. To leading order in 1/PR(J) the expected kept tensor energy is e(k)^m, with e(k) = Σ_{p≤k}σ_p²/Σ_pσ_p².

*Proof.* Write the source as Sym(Z) with Z i.i.d. Gaussian, as in [D-TS] §2.4. The projected transported tensor energy factorises over legs into ∏ tr(P_kJJᵀ)/tr(JJᵀ). The symmetrisation cross terms are of relative order 1/PR(J), since tr(A²) ≤ ‖A‖·tr A. ∎ This is a statement about tensor energy. The D21-level version is E4.

Consequences:
- All-distinct κ4 (m = 4) needs more modes than all-distinct κ3 (m = 3) at equal accuracy. **DOA a fortiori.**
- The (2,1,1) slice plausibly transports differently (a heuristic, not derived; test E4). Its doubled index is not a scrambled leg.
  - At shallow layers the doubled→doubled channel W∘W flattens it (R1, §6.4).
  - At depth its profile is var_i, which arrives through the distinct→doubled channel (§6.5). That is a function of the layer's own variances, not a Haar-oriented leg.

  Either way only the legs j and k are scrambled, so as a transported object the slice should be *covariance-type*, with m = 2. This would explain why a form like u_iC_jk can work at all, and why [D-TS] §4 item 4 finds covariance-type content cheap.

**What is alive.** The propagator-only ensemble error formula ([D-TS] §2.4) predicts mode counts from the spectrum of J alone. It reproduces the measured k_2% at 11–12 of 15 ages. It is a free offline calculator for any design that projects transported content, including the backward questions of §5.3.

### 3.2 The rank-one spike along the mean direction

**Cited facts.**
- The spike behaves as a BGN outlier above the free-probability threshold. This is the digest's *interpretation*, not a cited fact ([D-TS] §2.6, B6): the BGN theorem is cited, but the identification is marked unproved, because the spike is generated by the coupling of gates and weights rather than given as an independent perturbation.
- At n = 1024 and age 15, ⟨u₁, μ_z⟩² = 0.85 and s₁/s₂ ≈ 1.45 ([D-TS] §3.4).
- The top finite-time Lyapunov exponent is +0.059 per layer, while the bulk is discounted by the gain 2E[Φ²] = 0.58–0.95 per layer ([D-TS] B8).
- The old pool's leading tensor mode is the mean direction (cos 0.93, 60–80 % of the tensor energy at depth), but "it carries almost no D21" ([OC] Results).

**DERIVED (the spike as a standalone carrier; arithmetic).** Suppose a carrier transports the spike part of a tier exactly and drops the rest. Then ε_tot = share·√(1 − f), where share is the tier's share of D21 and f is the fraction of the tier's D21 energy along the spike. The inputs:
- share = 0.5–0.8 for w = 4 at depth; measured at n = 128–256 and flat across that range, extrapolated to n = 1024 ([OC] table (1));
- f ≤ 0.9 in the best measured case (slice-containing sources, ages ≥ 6, [D-TS] §3.6); [OC]'s convention gives f ≈ 0.

Hence ε_tot ≥ 0.5·√0.1 ≈ 0.16, against 0.022. **DOA as a standalone carrier.**

The rank-one spike is not the κ3 gain template T of §7.5. T has one leg along μ and the other along var or C. In the published convention, T carries most of the old pool's D21 (§7.5); the spike does not.

**DERIVED (deflation).** Carry the spike exactly and the bulk in the propagator basis. The bulk's own target relaxes from ε to ε/√(1 − f). The measured f depends on the source and the age ([D-TS] §3.6, n = 128, MLP 0):
- full source (slices kept): the top direction carries 0.65 of the transported D21 energy at age 8 and 0.85–0.90 at ages 11–15;
- all-distinct source: 0–0.22 at ages ≤ 9.

With f = 0.65 at age 8, the target relaxes from 2 % to 3.4 %: ≈ 290 modes instead of 322 at n = 1024 ([D-TS] §3.4, interpolated), ≈ 10 % fewer. With f = 0.85–0.90 at age 15, it relaxes to 5.2–6.3 %: ≈ 150 modes instead of 197, ≈ 20–25 % fewer. The mode counts are for the all-distinct source, so pairing them with the full source's f is itself optimistic. Either way the count stays ∝ n. **A small lever.**

**Link to the gain mode (measured, C9; ANALOGY with an exact breaking point).** The variance excess of C along μ̂ can be compared with the scale-mixture prediction (c4/4)·‖μ‖² of §7.1:
- at layer 15 it is 1.05–1.6× the prediction (n = 128–1024, eight networks);
- at layer 10 it is 1.4–2.2×;
- at layer 2 it is 2.4–3.3×.

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
- *Small entries do not give bounded η.* Take Cor = I + θ(vvᵀ − I/n) with v = n^{−1/2}·1, i.e. every off-diagonal correlation equal to θ/n. An exchangeable mixture realises it on {±1}ⁿ for θ ≤ n: draw a fair sign σ, then i.i.d. spins with mean σ√(θ/n). The entries are θ/n, but η = λ_max(Cor) − 1 = θ(1 − 1/n) is arbitrary.

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
- gate_GX = E[1[z_i>0]·(a_j − μ_j)].

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

### 4.2 DERIVED + checked: the 2-clique Gibbs closure misstates three-gate cumulants (by up to 4/π for equal thresholds, without bound for mixed ones)

**Setting.** Let (z₁, z₂, z₃) be Gaussian with unit variances and correlations r_ab = O(ρ). Let the gates be g_a = 1[z_a > t_a], with p_a = P(z_a > t_a), v_a = p_a(1 − p_a) and φ_a = φ(t_a). Let ν be the maximum-entropy law on {0,1}³ with the same one- and two-gate marginals: a pairwise binary MRF, i.e. the "2-clique Gibbs" description of the face law.

**Statement.** To leading order O(ρ²):
- exact: κ3(g₁, g₂, g₃) = φ₁φ₂φ₃·[t₂ r₁₂r₂₃ + t₁ r₁₂r₁₃ + t₃ r₁₃r₂₃];
- pairwise MRF: κ3_ν = φ₁φ₂φ₃·Σ_{centre b} r_ab·r_bc·V(t_b), with V(t) = φ(t)(1 − 2p)/(p(1 − p)).

For equal thresholds the ratio κ3_ν/κ3 = φ(t)(1−2p)/(t·p(1−p)) lies in (1, 4/π]. It tends to 4/π = 1.273 as t → 0 and to 1 as |t| → ∞.

For unequal thresholds write V(t) = f(t)·t, with f(t) = φ(t)(1−2p)/(t·p(1−p)) even and in (1, 4/π], and w_b = r_ab·r_bc. The ratio is Σ_b f(t_b)·t_b w_b / Σ_b t_b w_b.
- If all t_b w_b have one sign, it is a weighted mean of the f(t_b) and stays in (1, 4/π].
- If they have mixed signs (thresholds of both signs, or correlations of both signs, as in a random layer), the exact sum can cancel while the pairwise one does not. The ratio is then unbounded and can be negative (C13).

*Proof.*
- **Exact law.** Mehler's expansion (memory) gives 1[z > t] = p + Σ_{k≥1} c_k He_k(z) with c_k = φ(t)He_{k−1}(t)/k!, so c₁ = φ(t) and c₂ = φ(t)t/2.
  - At order ρ² the joint cumulant comes only from paths a–b–c with gate b using He₂ and gates a, c using He₁.
  - E[He₁(z_a)He₂(z_b)He₁(z_c)] = 2r_ab·r_bc + O(ρ³).
  - The contribution is therefore c₁(t_a)·c₂(t_b)·c₁(t_c)·2r_ab·r_bc = φ_aφ_bφ_c·t_b·r_ab·r_bc.
- **Pairwise MRF.** Expand the binary pairwise MRF at weak coupling.
  - Couplings: J_ab = Cov(g_a, g_b)/(v_a v_b) + O(ρ²), with Cov(g_a, g_b) = φ_aφ_b·r_ab + O(ρ²).
  - The three-point cumulant is a tree through a centre gate whose own third cumulant is κ3(g_b) = v_b(1 − 2p_b): κ3_ν = Σ_b J_ab J_bc·v_a·κ3(g_b)·v_c + O(ρ³).
  - Substituting the couplings gives φ_aφ_c·r_ab·r_bc·φ_b²(1 − 2p_b)/v_b.
- **The limits.** As t → 0, 1 − 2p ≈ 2φ(0)t and p(1−p) → ¼, so the ratio tends to 8φ(0)² = 4/π. As |t| → ∞, Mills' ratio gives the limit 1. That f is monotone between the two limits, hence in (1, 4/π], is checked on a grid of 8·10⁴ points in (0, 8], not proved.
- **Unequal thresholds.** Substitute V(t_b) = f(t_b)·t_b in the two leading terms. ∎

**Checks (C2).** Exact trivariate probit probabilities by nested Gauss–Legendre quadrature, against the fitted pairwise max-ent law.

| quantity | values |
|---|---|
| κ3_ν/κ3_exact, five threshold triples in [−0.5, 2] (including mixed signs), ρ = 0.05, 0.1, 0.2, correlations ρ·(1, 0.8, 0.6) | 1.21–1.28 at 14 of 15 points |
| the same at t = (2, 2, 2), ρ = 0.2 (higher orders enter) | 1.41 |
| leading-order factor V(t)/t at t = 0.5, 1, 1.5, 2, 3 | 1.264, 1.238, 1.200, 1.159, 1.093 |
| re-run (C13) at (0.4, 0.4, 0.4), (0.5, 1, 1.5), (−0.5, 0.5, 1), (2, 2, 2), (1, −0.5, 0.25), ρ = 0.05, 0.1, 0.2 | 1.12–1.28 at 14 of 15 points; 1.41 at (2, 2, 2), ρ = 0.2; the Möbius coefficients at t = 0.4 reproduced to four digits |
| mixed thresholds, equal correlations ρ (C13): t = (−1, 0.5, 0.5), where the exact O(ρ²) term vanishes | 11.4, 5.1, 2.95 at ρ = 0.02, 0.05, 0.1 |
| t = (1, −0.48, −0.5): the pairwise law gets the sign wrong | −0.26, −0.71, −3.2 at ρ = 0.02, 0.05, 0.1 |
| t = (−1, 0.6, 0.6) | 1.31–1.36, above 4/π |

**Why the max-ent description fails.** The exact gate law is not a pairwise MRF at the order that matters. Its three-body Möbius coefficient of log P, the HC interaction on the triple ([D-HC] §3.1, Lemma A), is O(ρ²): −1.985e-3, −7.288e-3, −2.497e-2 at ρ = 0.05, 0.1, 0.2 (t = 0.4). That is the same order as the three-gate cumulant itself. Setting it to zero, which is what "Gibbs with 2-cliques" means, misstates the three-gate cumulants by 9–27 % at leading order for equal thresholds with |t| ≤ 3, and arbitrarily, even in sign, for thresholds of mixed sign.

A closure built on that description is worse than the probit (Gaussian) closure that the leg-partition expansion already uses. Including the 3-cliques restores accuracy, but that state is a set of n³ parameters: the trivariate object itself, with no saving.

### 4.3 Positivity fails (DERIVED + checked); no safe symbol and a complete graph (measured), so Hammersley–Clifford is vacuous

**Statement.**
- (i) DERIVED + checked. For l ≥ 1 and invertible W_l, the face law of layer l gives zero probability to at least one orthant.
- (ii) Measured at n = 6 and 8 on the sampled support (C14); CONJECTURE at n = 1024 (test: repeat C14 at n = 10–14). At layers ≥ 1 the support is not closed under switching any set of gates to a common vacuum, for any vacuum: there is no safe symbol. The one exception (n = 6, layer 2) has a dead gate. Every pair of live gates is *conditionally* dependent given all the other gates, with one exception per width where the support is very sparse, so the Markov graph is complete or nearly so.

*Proof.*
- Write z_l = a_{l−1}W_l with a_{l−1} ∈ ℝⁿ₊. Then z_l lies in the convex cone K = {aW_l : a ≥ 0}.
- Put c = W_l^{−1}1. Then (aW_l)·c = a·1 > 0 for every a ≥ 0 with a ≠ 0, so K lies strictly on one side of the hyperplane c^⊥.
- Every point y of the open orthant with sign pattern −sign(c) (generic c has no zero entries) has y·c < 0, so y ∉ K.
- The gate pattern 1[−c > 0] therefore has probability zero. ∎

The draft also argued completeness from Cor(g_i, g_j) ≠ 0. That argument does not work. The Markov graph is defined by *conditional* dependence given all the other gates, and marginal correlation neither implies nor excludes it. (ii) is therefore a measurement, not a theorem.

**Checks (C3).** The orthants met by K, sampled with a ≥ 0 at random, are lower bounds limited by sampling:

| n | orthants met by the cone of rows of W (8·10⁵ draws of a ≥ 0) | met by Gaussian input, layer 0 (4·10⁵ draws) | 2ⁿ |
|---|---|---|---|
| 8 | 105 | 250 | 256 |
| 10 | 317 | 896 | 1024 |
| 12 | 1375 | 3476 | 4096 |
| 14 | 4606 | 10974 | 16384 |

**Checks (C14).** Gaussian input, 4·10⁶ samples per network. The plug-in conditional mutual information I(g_i; g_j | all other gates) is compared with a null that shuffles g_j within each context. The safe-symbol test is exhaustive over vacua and subsets on the sampled support.

| n, layer | support | pairs conditionally independent (CMI < 3 × null maximum) | vacua with safe-symbol closure |
|---|---|---|---|
| 6, layer 0 | 64 of 64 | 0 of 15 | all 64 (the law is positive) |
| 6, layer 1 | 27 of 64 | 0 of 15 | 0 |
| 6, layer 2 | 18 of 64 | 5 of 15, all involving one dead gate (P(on) = 0) | 2 |
| 6, layer 3 | 14 of 64 | 6 of 15: five involving a dead gate, one between live gates | 0 |
| 8, layers 0 / 1 / 2 / 3 | 256 / 116 / 99 / 52 of 256 | 0 / 0 / 0 / 1 of 28 | 256 / 0 / 0 / 0 |

A second n = 8 network gives 0 of 28 conditionally independent pairs at layers 1–3 and no safe symbol at layer 1. In every case the orthant 1[c < 0] predicted by (i) is never sampled.

**Consequence.** The Hammersley–Clifford proof uses less than positivity: only that the support be closed under switching any set of sites to a vacuum (Pos-safe, [D-HC] §3.3). Positivity failure alone therefore does not settle the matter; the measured absence of a safe symbol does, at small n.
- Without a safe symbol, Markov ⇏ Gibbs in general: the intersection axiom fails, and Moussouris-type counterexamples exist ([D-HC] §4.2–4.3).
- The positive results without positivity ([D-HC] §4.5, the pivot property) were not tested.
- Even where a factorisation exists, its graph is complete or nearly so, and "Markov with respect to the graph" says nothing.

There is no sparse separator structure inside a layer to exploit. This matches the guard recorded in [D-CR] §5.3 and §7.1: the boundary gain of Yang's CMI bound is isoperimetric and vanishes under edge expansion, and a dense layer is the extreme case.

### 4.4 THEOREM (elementary): face-law descriptions are blind to the radius and cannot hold the gain

Gates are 0-homogeneous: 1[z_l(tx) > 0] = 1[z_l(x) > 0] for t > 0. So the face law of every layer, and the joint law of the face history, is a function of the direction x̂ alone. It carries no information about the scale.

The same goes for any statistic computed from gates alone. A statistic that also uses z, such as gate_GX, carries the scale only through its z factor.

The gain mode of §7 is a different object, and 0-homogeneity alone does not cover it. It is the fluctuation of ‖a_l(x̂)‖ across directions on the sphere; the input radius is the part that §7.1 (vi) shows cancels exactly.
- The gain equals ‖M_h x̂‖, with M_h the linear map of the face history h. It depends on the history *and* on the position of x̂ inside the cone of h.
- So it is not a function of the faces, and a state made only of gate statistics cannot hold its law (elementary).
- The faces do carry part of it: E[‖M_h x̂‖² | h] = tr(M_hᵀM_h·E[x̂x̂ᵀ | h]) varies with h. A face-only state can at best supply this conditional mean, not the fluctuation inside cones.

§7 shows that the gain fluctuation's shape dominates the (2,1,1) slice at depth at width 128. How much of it the faces' conditional mean carries was not measured.

**Verdict on item 2: DOA.**
- No pairwise state carries the binding content (THEOREM).
- The pairwise Gibbs closure is biased by 9–27 % at leading order for equal thresholds, and arbitrarily, even in sign, for mixed ones (DERIVED + checked).
- Hammersley–Clifford does not apply (positivity: DERIVED; safe symbol and completeness: measured at small n).
- The largest single shape in the slice at moderate width is not a function of the faces (elementary).
- Spectral independence certifies mixing, not accuracy (§3.3). What *does* carry a large part of the slice is a statement about the law of z, not of the gates: §7.

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

### 5.2 What can be recovered from the retained state: one exact symmetry, then nothing

**Setting.**
- S_t is any retained state at layer t: statistics of z_t together with the young content carried exactly.
- X_t is the erased content: the old-source κ3 transported to t, or its D21 image.
- R is any recovery map applied to S_t.

**(a) Scale-mixture model + measured: the gain-template part of the erased content *is* recoverable, by a symmetry.**
- *The mechanism, exact inside the model of §7.1.* If z_l = s·y_l with s > 0 independent of y_l, positive homogeneity carries s, and with it its κ3 signature c·T_l (T_l = var_aμ_b + 2μ_aC_ab), from one layer to the next with c unchanged (§7.1 (v)). The component of X_t along T_t is then a function of the current (μ, var, C), which every design carries, and of one scalar. For the network the model is an approximation: the per-direction gain is not independent of the direction, and the measured coefficient is discounted (below). This is a model statement checked by measurement, not a derivation.
- *Where the scalar comes from.* It is accumulated from the sources as they age out of the young window, while they are still carried.
- *Measured* (C12, C16, §7.5; n = 128, three MLPs, published no-AD convention).
  - The old pool's D21 is 49–97 % along T (77–94 % on the third MLP).
  - The undiscounted accumulator leaves 0.11–0.34 of ‖D21‖, against 0.20–0.88 if the pool is dropped.
  - The true coefficient is discounted, not exactly conserved: the best scalar is 0.74–1.0 times the accumulated one, i.e. a discount of a few per cent per layer of age (≈ 2–5 %, inferred from the overshoot, not measured source by source).

This is a recovery map. It needs no detailed balance and no time average, only an exact symmetry of the arrows. It is the positive answer that item 3 admits.

**(b) DERIVED under hypothesis H′: the rest has zero Bayes-optimal recovery.** Let Y_t = X_t − (its T-component).

*Statement.* E‖Y_t − R(S_t)‖² ≥ E‖Y_t − E[Y_t|S_t]‖², with the expectation over the network ensemble. Under hypothesis H′, E[Y_t|S_t] = 0. Then no recovery map does better than dropping Y_t.

*Hypothesis H′.* Conditional on S_t, including the accumulator, the law of Y_t is invariant under Y_t ↦ −Y_t.

*Proof.* The first inequality is the defining property of the conditional mean. The second claim follows because H′ makes the conditional law symmetric. ∎

**Status of H′.** It is an assumption with weak evidence.
- *What the flip gives at one layer.* E_W[W^{⊗3}] = 0 ([D-TS] §2.5). The flip W_{s+1} ↦ −W_{s+1} reverses the sign of third-order transported content and leaves every even-order statistic unchanged. But it also reverses the mean, μ ↦ −μ, which is part of S_t. So it shows only that odd content has zero mean given the *even-order* statistics. It says nothing about E[Y_t | S_t] once S_t contains μ. The T-part shows the difference: T is odd in μ, so c·T flips together with μ, and conditioning on μ predicts it. That the non-T remainder has no such predictor is exactly what H′ assumes, not something the flip shows.
- *Downstream.* The downstream gates break the flip, since relu(−z) ≠ relu(z). The flipped network is another network, not a symmetry of the conditional law.
- *Why the T-part is excluded.* The coefficient of the gain mode is not sign-symmetric: Var(s²) > 0, and under the flip T and κ3 change sign together, so c = ⟨κ3, T⟩/‖T‖² does not. That is why H′ excludes it, and why (a) works.
- *Measured support.*
  - In the AD convention, which re-attributes the gain-mode slices to young sources, the old pool is nearly orthogonal to T: cos² ≤ 0.33 (C12).
  - [OC]'s absorption test (no-AD, w = 4): births remove at most ≈ 10–15 % of the pool's D21. Young sources remove ≈ 30–40 % in-sample, but this does not transfer; on the held-out MLP the residual stays at the drop-it level (0.35–0.79, against shares of 0.39–0.88).
  - Neither feature set contained T, so the test supports H′ in the directions it tried, not for the whole pool.

**CONJECTURE (testable form of H′).** Fit absorption coefficients for the non-T residual layer by layer on independent MLPs. Their relative signs are random: they agree no more often than a fair coin. Test: E3.

**DERIVED (the Petz map, one line).** Classically, the Petz recovery of a marginalisation channel relative to a reference law σ returns σ's conditional completion of the observed marginal. If σ is the closure's own model, the map returns the closure's births (tree-level content), which the chain already computes. A recovery map passes on the information in its reference; it does not create information.

**Verdict on recovery.**
- Detailed-balance and time-averaged recovery (§5.1): DOA.
- Symmetry-based recovery of the gain-template part: alive, measured (a).
- Recovery of the remainder: DOA under H′.

### 5.3 What does transfer: duality (ANALOGY + elementary identity + DERIVED)

**ANALOGY** (an in-programme connection, not a published one). The idea of [CR] that [A] §0.4 isolates is "recovery is not mixing": the defect is paired against a vector of small Dirichlet energy. Its commutative analogue already exists in the fresh-slate work, as the Duhamel telescoping of the heisenberg design ([FU] unlock 36, heisenberg Theorems 1–2, which [FU] calls "the commutative analogue of Chen–Rouzé's telescoping against the erased state"). The analogy breaks at the smallness. In CR the paired vector is small because a KMS-detailed-balanced generator controls its Dirichlet energy. Here there is no generator and no Dirichlet form; whatever smallness there is comes from the model error of the question.

**THEOREM (elementary).** For any reference chain ν with ν_0 = μ_0 and ρ_l = T_{l−1}ν_{l−1}:

  E_{μ_L} r_j − E_{ν_L} r_j = Σ_l (ρ_l − ν_l)[g_{l,j}].

*Proof.* Let g_{l,j} = T_l^*⋯T_{L−1}^* r_j be the exact pulled-back question, so g_{l−1,j} = T_{l−1}^* g_{l,j}. Then ρ_l[g_{l,j}] = ν_{l−1}[g_{l−1,j}], and the sum telescopes to ν_0[g_{0,j}] − ν_L[g_{L,j}] = E_{μ_L} r_j − E_{ν_L} r_j. ∎

Two properties matter here:
- All old content sits inside the exact pulled-back question g_l.
- A model ĝ_l makes the error bilinear: (one arrow's local defect) × (the question's model error).

**DERIVED (for the linear part of the question).** The gradient of the pulled-back question of a final readout, J_{l→L}(x)ᵀ∇r at each point with J the gated propagator, lies in the range of J_{l→L}ᵀ. J and Jᵀ have the same singular values. A rank-k compression of the backward question therefore keeps e(k) per leg, exactly as a forward carrier does (§3.1), and has the same participation ratio, PR ≈ n/(2·age) ([D-TS] §2.2).

Duality changes three things:
- which side carries the memory;
- the shape of the error, from additive to bilinear, which is its real gain;
- the need for intermediate n³ tensors: the heisenberg design evaluates the backward pairing at full rank by matrix products ([FU] unlock 35, "States forward, questions backward").

It does not change the mode count. Any *compression* of the backward question meets the same free-probability law as compression of the forward content.

**CONJECTURE (Haar pairing).** Take dropped old content and closure error of equal D21 norm. They cost the same final-layer MSE, following the same 4.2e-6·ε² law within 30 %. The reason: in both cases the orientation relative to the downstream propagator is Haar ([D-TS] §2.4), so the pairing with the pulled-back question sees only the norm. Test: E8. If dropped old content turns out to cost half or less per unit of D21, then the pairing helps, and the D21 share overstates what a chain must carry.

### 5.4 What a CMI-type decay statement predicts for how much old content a chain may forget

**Measured** ([OC] tracker, MLP A1, n = 128, t = 15). The shares of ages 1–6 are [OC]'s file `results/tracker_A1.txt`. Everything else is C16: the same tracker re-run on two independent regenerated atlases of A1, which reproduce [OC]'s shares to ±0.001.
- *Per-age shares of ‖D21‖, AD convention:* 0.888, 0.276, 0.250, 0.247, 0.190, 0.140 for ages 1–6; then 0.131, 0.073, 0.092, 0.041, 0.033, 0.030, 0.027, 0.015, 0.004 for ages 7–15.
- *No-AD convention:* 0.171, 0.104, 0.132, 0.159, 0.121, 0.120; then 0.147, 0.054, 0.113, 0.088, 0.067, 0.060, 0.106, 0.070, 0.011.
- *The old sources are not mutually orthogonal.* Their off-diagonal D21 contributions have pairwise cosines of +0.38 on average (AD; range +0.03 to +0.75) and +0.66 (no-AD). The nine sources of ages 7–15 add up to 0.317 (AD) and 0.616 (no-AD), against quadrature sums of 0.19 and 0.26. [OC] Results calls them nearly orthogonal; the direct measurement does not support that, and the draft's "each of the nine carries ≈ 0.314/√9 ≈ 0.10" does not hold.
- *A fixed source shrinks as it ages,* for every birth layer 0–11. In the AD convention its D21 contribution drops ≈ 3× at the first transport (the all-distinct projection), then by a factor 0.7–0.85 per layer: rate 0.15–0.4, typically 0.25. In the no-AD convention the factor is 0.75–0.9 per layer: rate 0.1–0.3, typically 0.2.
- *Width:* the w = 4 share is 0.62–0.79 at n = 128 and 0.49–0.77 at n = 256, flat from n = 128 on ([OC] table (1)).

**DERIVED (conditional arithmetic).** Suppose a Yang / Chen–Rouzé-type bound held for the depth chain: the dependence of layer t on content born at t − d, given the retained state, is ≤ C·e^{−Δd}. A chain could then drop content older than d* = ln(C/ε)/Δ.
- *The rate.* The draft took Δ to be at most the consecutive Lyapunov gap, 0.013–0.018 per layer ([D-TS] §0.1), which gives d* ≥ 135 layers. That is the wrong rate. A gap separates singular values; what shrinks old content relative to the current state is the discount of bulk directions by the gates, 2E[Φ²] = 0.58–0.95 per leg per layer ([D-TS] B8). For three legs that allows Δ between 0.08 and 0.8 per layer, and the measured rate is 0.1–0.4.
- *The horizon.* With C ≈ 0.25 (the AD age-2 share, as in the draft), ε = 0.022 and Δ = 0.15–0.3, d* ≈ 8–16 layers: comparable to L = 16.
- *What does not help.* The Sym^k functor adds no higher-order gap ([D-TS] §2.1), and the top direction is amplified (+0.059 per layer), not discounted. The free-regime law (concentration ∝ 1/age, [D-TS] B2) concerns the propagator's mode count, not the size of the content.

**Measured, in agreement.** Share of ‖D21(z_15)‖ carried by all sources of age ≥ d (C16):

| d | 7 | 10 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|
| AD | 0.317 | 0.120 | 0.066 | 0.041 | 0.016 | 0.004 |
| no-AD | 0.616 | 0.378 | 0.235 | 0.181 | 0.079 | 0.011 |

Only the sources of age ≥ 14 (AD) or of age 15 (no-AD) fall below the bar.

**Prediction (CONJECTURE).** At n = 1024 the per-age profile is the same within network scatter: the rates are set by gate discounts, which are mean-field, and the w = 4 pool share is flat from n = 128 to 256. A chain at L = 16 may then drop the one or two oldest sources and nothing else. Test: E2.

**Verdict: DOA as a forgetting principle at this depth.** Old content does decay, at a rate set by the gates' discount rather than by a gap, but its horizon is comparable to the depth. Over 16 layers forgetting removes one or two of the eleven old sources (age > 4).

---

## 6. Item 4: the cluster / polymer view of the closure

### 6.1 KNOWN-LINK + identification: the closure is a first-order cluster expansion with exact transport

**How the expansion arises.**
1. Expand each relu(z_i) in Hermite polynomials around its own Gaussian marginal (Mehler, memory). The vertex weights are w₁ = Φ, w₂ = φ/σ, and so on.
2. Write the joint cumulants of a = relu(z) by the Leonov–Shiryaev formula (memory). They become sums over *connected* diagrams:
   - vertices are neurons, weighted by Hermite coefficients;
   - edges and hyperedges are cumulants of z.

"First order" keeps diagrams with at most one non-Gaussian hyperedge plus covariance edges. The leg-partition coefficients count the ways the legs of a vertex can be distributed among them.

**The identified first-order diagrams** ([CP] §3.1 (v)): Gaussian ρ² + Φ³κ3(z) + [D21(z) ⊗ C] + [K22 ⊗ C] + [κ4(z)_(2,1,1)]. The diagram expansion is the published part (Leonov–Shiryaev); the identification of the closure's coefficients with leg-partition counts is the programme's. At n = 1024 the fitted coefficients stay within 0.15 of the leg-partition values at every depth, with one exception ([O1024]): B3, the D3(z) hyperedge with two edges, is fitted at 0.15–0.40 against a counted 1, at every width. [O1024] lists it as an open discrepancy in a small term. The identification therefore holds for five of the six diagram classes, not exactly.

**ANALOGY with Hammersley–Clifford's Lemma A** ([D-HC] §3.1). Both rest on Möbius inversion, with three differences:
- the lattice: HC inverts on the Boolean lattice of subsets, cumulants on the lattice of set partitions;
- the object: HC inverts log P and gets interaction potentials, while cumulants invert moments and get connected correlations;
- the dual sides: log versus linear.

The two agree only at leading order in weak coupling, and §4.2 measures where they part. The three-body potential and the three-gate cumulant are both O(ρ²), but with different coefficients.

**Convergence** (ANALOGY with Kotecký–Preiss, memory). A polymer expansion converges when polymer activities are small compared with the entropy of polymers. Here the activity of an edge is an entrywise correlation, O(n^{−1/2}). That is why the first-order truncation error falls as a power of n (n^{−0.8} to n^{−0.9} measured) rather than with a temperature. It breaks as an analogy because there is no lattice entropy to beat: the number of diagrams per vertex grows with n, and only the 1/√n edge weights win.

### 6.2 ANALOGY: Bethe / tree level on pseudorandom geometry, and where it breaks

**The similarity.** The closure keeps tree-level gate diagrams and treats the dense, loop-rich part exactly through matrix products (the transport). Bethe / cavity does the same on locally tree-like graphs.

**Where it breaks: the reason the approximation is valid.**
- Bethe is exact when loops are long (girth → ∞) and couplings are O(1).
- Here every triple of neurons in consecutive layers closes a short loop through the complete bipartite connection.
- The truncation is accurate because each edge weighs O(n^{−1/2}), not because loops are rare.

The relevant theory is the dense mean-field (TAP / Plefka) expansion in powers of the coupling, not the sparse cavity method (memory).

**Consequence for the bethe design.** Loop-series corrections in a dense layer are organised by powers of n^{−1/2}, and the shortest loops dominate. E6 decomposes the closure residual by diagram class.

### 6.3 Measured fit + arithmetic: the closure error is not the binding error at n = 1024

**The error law.** The fit ε(n) = 0.155·n^{−1/2} + 4.26·n^{−1} reproduces the deep-layer closure errors 4.7, 2.6 and 0.9 % at n = 128, 256 and 1024.
- It predicts 1.5 % at n = 512.
- So does the log-log fit through the three points (slope −0.79).
- More widths between 128 and 1024 therefore cannot tell the two forms apart.
- The final [O1024] report gives deep-layer ranges of 3.7–7.1, 2.0–2.9 and 0.74–0.91 % at the three widths (midpoints 5.4, 2.45 and 0.83 %), i.e. about n^{−0.9}. The two-term fit was made on earlier midpoints. The conclusion below does not depend on which set is used.

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
- **A tripled output index receives tripled input content through W∘W∘W, of norm O(1/n).** That does not make the channel negligible, because its target is small too.
  - D3(z') is O(1/n) per unit: measured rms 15–19/n at n = 128–512 (C15).
  - D3(a) is O(1) per unit wherever gates are uncertain. At layer 0, κ3 of relu of an N(0, 2) unit is 2^{3/2}·0.3265 = 0.923.
  - The channel's output is then ≈ √120·0.923/n ≈ 10/n per unit, against ≈ 17/n for the target.
  - Measured (C15): the tripled channel carries 0.60–0.67 of ‖D3(z₁)‖ at layer 1 at all three widths, 0.33–0.35 at layer 2, 0.18–0.22 at layer 3, and 0.02–0.10 at depth, falling with n.

  Diagonal content is therefore inherited at shallow layers and mostly regenerated at depth. The draft's "never inherited" was wrong.

**The Perron mode is the norm-coupling channel.** The flat part transmits Σ_p κ(a)_{pp…}. By the trace identities of §7.1 (iv), this is the covariance of ‖ã‖² with the remaining legs, minus Wick terms: the norm-coupling channel (exact). In the scale-mixture model of §7 that channel is the dilation (gain) channel. This is where the expander picture and positive homogeneity meet (interpretation). What it carries is the gain's *amplitude*. The var_i profile of the (2,1,1) slice at depth arrives through the distinct→doubled channel instead (§6.5).

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
  | tripled, re-run (C15: fresh networks at n = 128, 256, 512; norm ratios, two replicas) | 0.60, 0.67, 0.66 at layer 1; 0.02–0.10 at depth, falling with n |
  | doubled | 0.65–0.98 |
  | distinct | 0.31–0.50 |

  The draft's 0.36 at layer 1 disagrees with the re-run and with the analytic estimate of ≈ 0.6 (≈ 10/n against ≈ 17/n per unit). It equals 0.60², so C5 may have reported an energy fraction; the atlas networks were not re-run.
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

  So the *gain excess* gex := Var‖z̃‖²/(tr C)² − 2‖C‖²_F/(tr C)² equals Σ_ij κ4_iijj/(tr C)², the normalised trace of the (2,2) slice. For a scale mixture of a Gaussian, gex = Var(s²)·(1 + 2‖Σ‖²_F/(tr Σ)²). Hence c4/gex = 1/(1 + 2‖Σ‖²_F/(tr Σ)²) ≤ 1. The ratio is close to 1 unless a few eigenvalues dominate C; a strong spike lowers it.

  This inequality is a property of the model, not an identity. gex includes the diagonal entries κ4_iiii and κ4_iijj, while c4 is fitted on distinct-index entries only. Adding independent per-coordinate noise with negative excess kurtosis (uniform, say) lowers gex and leaves the distinct-index slice itself unchanged, so outside the model c4 > gex is possible.
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

### 7.3 Across widths up to n = 1024 (C9)

**Method.** Fresh He MLPs, one network per (width, seed), N = 2.5·10⁴ inputs per batch. Nothing of size n³ is formed.
- **‖K‖².** Estimated from two independent batches A and B with the Gram identity ⟨M4_A, M4_B⟩_AD = mean_{s,t}[p₂p₁² − 2p₃p₁ − p₂² + 2p₄], where p_m = Σ_i (z_si z_ti)^m. This is Möbius inversion over the partitions of {i, j, k}. Each batch's cumulant uses its own var and C.
- **Templates.** R and C_off come from a third batch; its two halves give unbiased template norms.
- **What is reported.**
  - share_R: the energy of the slice along R, i.e. what one scalar per layer can carry.
  - share_Creg: the energy captured by the published family u_iC_jk with u free.
  - c4: the slice's amplitude on R.
  - c3: the D21 amplitude on T.
  - gex: the gain excess of §7.1 (iv).
  - share3_T: D21's energy along T.
  - The ratio of the measured variance excess along μ̂ to the gain prediction (c4/4)‖μ‖².

| n | seed | share_R at l = 2 / 6 / 10 / 15 | share_Creg (l = 15) | cos(u, var) | c4·n/2 at l = 2 / 6 / 10 / 15 | gex·n/2 | c3·n | share3_T | spike excess / gain prediction |
|---|---|---|---|---|---|---|---|---|---|
| 128 | 21 | 0.30 / 0.51 / 0.59 / 0.82 | 0.81 | 0.94–0.99 | 2.9 / 9.1 / 8.7 / 11.0 | 3.3 / 9.1 / 12.1 / 18.2 | 5.0 / 7.5 / 6.8 / 6.7 | 0.48 / 0.77 / 0.85 / 0.91 | 2.4 / 1.4 / 1.6 / 1.4 |
| 128 | 22 | – / 0.47 / 0.79 / 0.91 | 0.91 | 0.93–0.98 | – / 5.7 / 7.5 / 7.4 | – / 6.3 / 11.2 / 14.7 | – / 5.8 / 5.0 / 4.4 | – / 0.85 / 0.93 / 0.96 | – / 2.0 / 1.5 / 1.6 |
| 256 | 21 | 0.31 / 0.37 / 0.39 / 0.64 | 0.57 | 0.96–0.99 | 2.8 / 5.6 / 8.8 / 13.5 | 3.0 / 6.0 / 9.7 / 14.7 | 5.0 / 6.9 / 7.9 / 8.2 | 0.60 / 0.76 / 0.78 / 0.85 | 2.7 / 2.1 / 1.4 / 1.05 |
| 256 | 22 | – / 0.37 / 0.40 / 0.66 | 0.63 | 0.97 | – / 5.4 / 6.3 / 8.3 | – / 6.0 / 7.6 / 12.4 | – / 5.9 / 6.3 / 6.0 | – / 0.76 / 0.81 / 0.88 | – / 1.8 / 1.8 / 1.6 |
| 512 | 21 | 0.25 / 0.30 / 0.24 / 0.38 | 0.36 | 0.95–0.99 | 2.5 / 4.8 / 6.0 / 10.3 | 2.8 / 5.9 / 7.7 / 11.9 | 4.4 / 6.9 / 8.5 / 9.0 | 0.45 / 0.77 / 0.76 / 0.79 | 2.9 / 2.4 / 2.2 / 1.4 |
| 512 | 22 | – / 0.37 / 0.54 / 0.66 | 0.60 | 0.97–0.98 | – / 5.9 / 10.4 / 13.7 | – / 6.4 / 9.9 / 14.5 | – / 7.7 / 9.1 / 9.7 | – / 0.78 / 0.85 / 0.89 | – / 2.1 / 1.4 / 1.2 |
| 1024 | 21 | 0.27 / 0.26 / 0.25 / 0.36 | 0.34 | 0.96–0.98 | 2.4 / 4.4 / 6.3 / 10.5 | 2.7 / 5.6 / 7.8 / 10.5 | 4.5 / 7.3 / 8.9 / 9.8 | 0.47 / 0.72 / 0.79 / 0.83 | 3.3 / 2.6 / 2.0 / 1.4 |
| 1024 | 22 | – / 0.28 / 0.30 / 0.39 | 0.36 | 0.95–0.98 | – / 4.6 / 6.6 / 9.3 | – / 5.7 / 8.1 / 10.1 | – / 7.9 / 8.9 / 9.3 | – / 0.76 / 0.82 / 0.86 | – / 2.7 / 2.0 / 1.5 |

**Readings.**
1. **The shape is the same at every width.** cos(u, var) = 0.93–0.99 everywhere. share_Creg ≈ share_R, so the free u buys nothing over u = c·var.
2. **c·n is width-stable to within network-to-network scatter** (±40 % at mid-depth). The gain amplitude is O(1/n) per layer with an O(1) coefficient, the scaling of the finite-width four-point vertex (§7.7).
3. **The template's share of the (2,1,1) slice falls with width.**
   - At layers 6–15 it is 0.47–0.91 at n = 128.
   - It is 0.24–0.66 at n = 512, across two networks.
   - It is 0.25–0.39 at n = 1024, across two networks.
   - The fall from 128 to 256 is clean. Beyond 256, network scatter is comparable to the trend.

   The slice energy per neuron falls as n^{−1.2} at layer 2 (‖K‖²/(n·var⁴) = 2.41, 1.06, 0.47, 0.21). At layer 15 it falls faster (1140, 269, 70.5, 29.7). This is consistent with the template part there riding on the spike, whose weight relative to the bulk decays with n (an interpretation, not separately measured).
4. **This is the mechanism behind the published regeneration's measured fall** (55 → 33 → 20 % of the κ4 gap at n = 128 → 256 → 1024, [O1024], a different network). What the regeneration can carry is the gain mode, and the gain mode is a falling share of the slice.

   At n = 1024, 60–75 % of the (2,1,1) energy is *not* the gain mode. That remainder is the binding object of [CP] §6b. E1 decides whether it is low-rank.
5. **The κ3 template does not collapse.** share3_T is 0.72–0.86 at layers 6–15 at n = 1024. At the D21 interface the scale-mixture shape remains the dominant part of the target at the competition width.
6. **c4 versus gex.** c4/gex = 0.77–1.01 at n = 1024 and 0.5–1.0 at n = 128. As the scale-mixture model predicts (iv), c4/gex ≤ 1 within noise, and it is lowest where the spike is strongest (n = 128 at depth). This is a consistency check of the model, not of an identity.

   The single-mixture relation c3 = c4/2 holds only to within a factor of about 2: c3·n/(c4·n/2) = 0.6–1.9.

   At layer 15 the variance excess along μ̂ is 1.05–1.6× the gain prediction (c4/4)‖μ‖², so the spike at depth is largely the gain mode. At layer 2 it is 2.4–3.3×, i.e. mostly BGN amplification of the mean (§3.2).

### 7.4 The amplitude is not a Gaussian quantity (C10)

**DERIVED prediction to test** (§7.1 (v)). If the direction law at each layer were Gaussian, the gain excess would obey gex_{l+1} = gex_l + Δ_{l+1}, with Δ_{l+1} = gex(relu(y_l)W_{l+1}) for y_l ~ N(μ_l, C_l).

**Measured** (n = 128 and 256, N = 2·10⁵ network samples and 10⁵ Gaussian samples per layer):

| layer | 1 | 2 | 4 | 6 | 8 | 10 | 12 | 15 |
|---|---|---|---|---|---|---|---|---|
| gex·n/2, network (n = 128) | 1.57 | 2.52 | 4.33 | 6.31 | 7.36 | 9.63 | 10.07 | 14.22 |
| Σ Δ·n/2, Gaussian-input recursion (n = 128) | 1.56 | 2.28 | 3.24 | 3.86 | 3.40 | 3.88 | 3.73 | 3.57 |
| same, with the scale-induced spike (gex/4)μμᵀ removed from C_l (n = 128) | 1.56 | 2.31 | 3.33 | 4.17 | 4.07 | 5.14 | 5.74 | 6.75 |
| gex·n/2, network (n = 256) | 1.66 | 2.83 | 5.29 | 8.34 | 11.14 | 12.81 | 15.15 | 18.70 |
| Σ Δ·n/2, Gaussian-input recursion (n = 256) | 1.69 | 2.62 | 4.21 | 5.10 | 5.45 | 5.53 | 4.87 | 4.50 |

**Readings.**
- Layer 1 is reproduced exactly, since z₀ is Gaussian.
- Beyond it, the Gaussian-input increments fall to ≈ 0 while the network's gex keeps growing by ≈ 0.9 (n = 128) to 1.2 (n = 256) per layer, in units of 2/n. At layer 15 the recursion gives 24–47 % of the measured amplitude.
- So the per-layer gain variance is fed mainly by the non-Gaussian structure of the direction law: higher cumulants interacting with the gates. It is not fed by the covariance.

**CONJECTURE refuted in its simplest form.** "The gain mode is carried by an O(n²)-per-layer scalar recursion from (μ, C) alone."

**What survives (CONJECTURE).** A scalar recursion for the (2,2)-trace gex = Σ_ij κ4_iijj/(tr C)² does exist in principle. Its transport goes through the doubled-leg flat mode (R1): tr₂₂(z') ≈ 4·tr₂₂(a) + O(n^{−1/2}). Its gate step needs the first-order births of the trace from κ3 and κ4 of z. Test: E10.

### 7.5 Old κ3 content lies along the gain template in the published convention (C12)

**Setting.** The old-content stream's exact per-source tracker ([OC], `tracker.py`, run on the n = 128 atlases) splits D21(z_l) into the contributions of sources of each age.

**Two conventions.**
- *AD*: every slice is re-attributed to the newest birth.
- *no-AD*: the published convention; transported sources keep their slices.

**The test.** Project the old pool (age > 4) onto the one-scalar gain template T_l = var_aμ_b + 2μ_aC_ab, built from the layer's own (μ, var, C).

| layer | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|---|---|---|
| MLP 0, AD: old share of ‖D21‖ | 0.08 | 0.13 | 0.16 | 0.19 | 0.28 | 0.36 | 0.44 | 0.45 | 0.41 | 0.39 |
| MLP 0, AD: cos²(old, T) | 0.00 | 0.01 | 0.05 | 0.18 | 0.27 | 0.33 | 0.26 | 0.09 | 0.11 | 0.17 |
| MLP 0, no-AD: old share | 0.22 | 0.38 | 0.52 | 0.55 | 0.65 | 0.76 | 0.81 | 0.72 | 0.75 | 0.82 |
| MLP 0, no-AD: cos²(old, T) | 0.61 | 0.82 | 0.89 | 0.90 | 0.95 | 0.97 | 0.94 | 0.92 | 0.94 | 0.95 |
| MLP 0, no-AD: ε after the best scalar c·T | 0.13 | 0.16 | 0.17 | 0.17 | 0.15 | 0.14 | 0.20 | 0.21 | 0.19 | 0.18 |
| MLP 0, no-AD: ε after the undiscounted accumulator | 0.13 | 0.17 | 0.19 | 0.22 | 0.24 | 0.26 | 0.31 | 0.31 | 0.28 | 0.28 |
| MLP 1, no-AD: old share | 0.20 | 0.28 | 0.34 | 0.42 | 0.58 | 0.65 | 0.84 | 0.80 | 0.81 | 0.88 |
| MLP 1, no-AD: cos²(old, T) | 0.49 | 0.65 | 0.74 | 0.79 | 0.87 | 0.89 | 0.93 | 0.91 | 0.94 | 0.94 |
| MLP 1, no-AD: ε after the best scalar c·T | 0.15 | 0.16 | 0.18 | 0.19 | 0.21 | 0.22 | 0.22 | 0.24 | 0.21 | 0.22 |
| MLP 1, no-AD: ε after the undiscounted accumulator | 0.15 | 0.17 | 0.19 | 0.22 | 0.24 | 0.24 | 0.27 | 0.32 | 0.29 | 0.28 |

(MLP 1 in the AD convention agrees with MLP 0: cos² ≤ 0.20.)

A third MLP (A1, seed 770100, the old-content stream's network; C16) replicates the no-AD rows at layers 6–15:
- old share 0.24, 0.43, 0.54, 0.58, 0.65, 0.68, 0.62, 0.69, 0.79, 0.75;
- cos²(old, T) 0.77, 0.83, 0.85, 0.86, 0.88, 0.91, 0.93, 0.94, 0.93, 0.89;
- ε after the best scalar 0.11, 0.18, 0.21, 0.22, 0.22, 0.20, 0.16, 0.17, 0.20, 0.25;
- ε after the accumulator 0.11, 0.18, 0.23, 0.24, 0.24, 0.22, 0.23, 0.26, 0.34, 0.33.

**How the accumulator works.** A source is projected once onto T_l, when it ages out of the young window, giving a_l = ⟨D21(source), T_l⟩/‖T_l‖². Then c_acc(l) = Σ_{l'≤l} a_{l'} is carried forward without discount and applied to the current template.
- This needs no fit and no extra transport.
- It overshoots the best scalar by 0–35 % (A1: c_best/c_acc = 0.74–1.0, lowest at layers 12–15). The T-coefficient of the old pool is therefore discounted by a few per cent per layer of age (≈ 2–5 %, inferred from the overshoot), not conserved exactly.

**Readings.**
- In the convention the published chain uses, the old pool is mostly the gain mode. Its D21 energy at layers 6–15 is 49–97 % along T, and 87–97 % at layers 10–15.
- **This was not tested before.** [OC]'s absorption features were birth diagrams and young sources. With them the in-sample residual stays at 0.25–0.62. The template T comes free from (μ, var, C). One scalar on it brings the old pool from 0.20–0.88 of ‖D21‖ down to 0.11–0.25 (best scalar) or 0.11–0.34 (accumulator), over the three MLPs.
- **The AD convention agrees.** That convention moves the gain-mode slices into the young sources, cos²(young, T) = 0.72–0.92. What it calls old is the non-T residual: 0.35–0.45 at depth, nearly orthogonal to T, as [BRIEF] §3 states.
- **What is left is still 5–15 times the bar.** The template does not remove the old-content problem. At layers 10–15 it shrinks the part that needs a real carrier by a factor of 2.3–3.1 in D21 norm with the accumulator, and 2.8–5.4 with the best scalar, which a chain does not have (n = 128, three MLPs).

**CONJECTURE (n = 1024).** Two predictions:
- In the no-AD convention the old pool stays ≥ 70 % (energy) along T at layers ≥ 8. Support: share3_T of the *total* D21 is still 0.72–0.86 at n = 1024 (C9).
- A carrier then needs to reach only ε_own ≈ 0.022/0.34 to 0.022/0.22, i.e. ≈ 6.5–10 %, on the non-T residual that the accumulator leaves at depth (≈ 9–14 % if the best scalar were available), rather than 2–3 % on the whole pool.

Suppose the residual's mode count follows the single-age law of [D-TS] §3.4 (age 8, n = 1024: 162, 214 and 263 modes at ε_own = 20, 10 and 5 %). At ε_own ≈ 6.5–10 % that is ≈ 215–245 modes, and a Tucker core costs ≈ 9–14 u per layer (≈ 6–11 u at 9–14 %). Compare 54–86 u for the w = 4 tier at 2 % without the template, and 3–4 u available. The draft's 4–9 u assumed ε_own = 10–20 %, which needs a residual of 0.11–0.22, below what the accumulator leaves at depth. This is rough: the residual is not a generic Haar-oriented source, and its own law must be measured. Tests: E9 (the residual's own mode count) and E2.

### 7.6 What this means for the competition

**κ4 (2,1,1).**
- The published regeneration u_iC_jk is the gain mode (cos ≥ 0.91), and one scalar per layer reproduces ≈ 90 % of its effect at n = 128.
- At n = 1024 the gain mode holds 0.25–0.39 of the slice's energy (two networks). Removing an energy share s leaves √(1 − s) of the slice's norm, so if the D21 error is linear in the slice error the gain mode can close only ≈ 13–22 % of the κ4 gap. The same conversion gives 27–70 % at n = 128 (s = 0.47–0.91), where the regeneration closes ≈ 55 %. On the [O1024] networks the regeneration closes 20 % at n = 1024.
- The remaining 60–75 % is the binding error at the competition width. Nothing in this note's bridges carries it: not the MRF (§4), not the expander (§6.5), not the gain (§7.3).
- E1 is the decisive next measurement: is the remainder dominated by one mode?

**Old κ3 content.** The gain template plus a scalar accumulator replaces most of the old tier's D21 effect at no cost (§7.5). The open question moves from "carry the old pool" to "carry the non-T residual". At depth that residual is 2.3–3.1× smaller in D21 norm with the accumulator (2.8–5.4× with the best scalar), at n = 128 on three MLPs.

**The amplitude.**
- For κ3, the accumulator supplies it at no extra cost. A Monte Carlo estimate of ⟨D21, T⟩ with the control variate (z̃² − var) needs ≈ 5·10⁴–1.2·10⁵ samples for 1 % at n = 1024 (C11). Pushing N samples through one dense layer costs N/n u, and the per-sample contraction Σ_ab z̃_a²T_ab z̃_b costs as much again. That is ≈ 100–230 u per layer, or 1.5–3.7 B over 16 layers, against ≈ 0.15 B for a whole leading submission. The draft's "≈ 2–4 u" was off by a factor of ≈ 25–60; Monte Carlo of the κ3 amplitude is not affordable.
- For κ4, Monte Carlo is hopeless: 3·10⁶–8·10⁶ samples for 1 % already at n = 128, and ≥ 10⁹ at n = 1024.
- The κ4 amplitude therefore needs its own closure (E10). An offline table of c·n against l is the fallback; it carries ±30–40 % scatter between networks.

### 7.7 Known links, analogies, and the NCG reading

- **KNOWN-LINK** ([D-TS] §2.3, Hanin–Nica Thm 1, Cor. 3). For a fixed input and random weights, ‖Mu‖² is log-normal with β = 5Σ1/n_i for ReLU. This is the annealed version of the gain accumulation.
  - The quenched quantity measured here (fixed weights, random input, excess over the Gaussian value) grows by ≈ 0.9–1.2 per layer in units of 2/n, i.e. ≈ 1.8–2.4/n per layer. The annealed log-variance grows by 5/n.
  - The two are different quantities. The annealed one includes network-to-network variation and the Gaussian part. Only the O(depth/width) scaling transfers.
  - The link breaks exactly at quenched versus annealed, the guard of [D-TS] B4.
- **KNOWN-LINK** (memory: Roberts–Yaida–Hanin, *The Principles of Deep Learning Theory*) for the published statement. The leading finite-width four-point vertex of a ReLU network is O(depth/width). In the ensemble it is the variance of the stochastic metric, i.e. a scale mixture over the previous layer's empirical second moment.
  - **ANALOGY** for the identification with the quenched gain mode measured here. The vertex averages over networks; the gain mode is one network's fluctuation over inputs. It breaks at quenched versus annealed, as in the previous bullet.
  - What transfers is the scaling: c·n is width-stable, as O(depth/width) predicts.
- **ANALOGY** (memory: Lebowitz–Percus–Verlet ensemble corrections). In a canonical ensemble, correlation functions differ from microcanonical ones by terms proportional to the fluctuation of the conserved quantity times derivatives with respect to it. The gain-mode cumulants (R for κ4, T for κ3, μμᵀ for C) are such terms for the scale.
  - It breaks because the scale is not conserved: every layer injects fresh variance (§7.4), and only its *past* fluctuation is conserved, by 1-homogeneity.
  - It breaks again in the measured discount of the old pool's T-coefficient, a few per cent per layer (§7.5).
- **SPECULATION (NCG).** For a positively homogeneous dictionary, dilations act on every face's state and commute with every arrow. The scale is then a central element of the arrow algebra, states disintegrate over it as a direct integral, and the gain mode is the first-order effect of the fluctuation of this central charge.
  - This is a property of the ReLU dictionary, not of the theory. A tiling or a Bratteli diagram has no dilation, so it fails the drag test. It belongs in the dictionary (prong 2), not in Notes 1–2.
  - It does carry one competition lesson. Any face-only description (R5) must be paired with the scale: a "face × gain" state, not a face state.

---

## 8. The user's intuition, tested on the competition object

The intuition: Gibbs states, the approximate local Markov property and MRFs share their "fundamental essence" with expander theory, the time-averaged detailed-balanced single-Pauli Lindbladian hides a bridge, and NCG has a big part to play. On the competition object it splits three ways. [A] §0 and §8 test it on the abstract stage.

**Literally a theorem.**
1. **The doubled leg is an expander** (R1, §6.4). In every cumulant transport the Hadamard square W∘W is a dense weighted graph √2 above Ramanujan, and the expander mixing lemma holds verbatim. Its Perron mode is the trace (norm-coupling) channel, which in the scale-mixture model is the dilation channel of §7. This is the one place in the competition object where "expander" is not a metaphor.
2. **The annealed layer is a perfect expander,** a twirl onto the scalars ([D-TS] §2.5). The score, however, is quenched: one Kraus operator, no gap.
3. **Spectral independence is local spectral expansion of the face law** ([D-HDX] §3.4, ALO). The measured η₀ is moderate (2–8 at n = 1024); that it stays bounded in n is consistent with the data but not established ([D-TS] §3.5). Even where it holds, it is a frame bound (R6) and certifies Glauber mixing, not the estimator.

**Two halves of one mechanism.**
1. **The closure is heredity plus decorrelation.** The layer chain is exactly Markov: the law of z_{l+1} is a push-forward of the law of z_l, and that is free. The closure is accurate because of decorrelation, which has two faces, the two halves of the transport:
   - single legs scramble orientation (Haar);
   - doubled legs average (the expander).

   The local-to-global step "local Markov + decorrelation ⇒ global accuracy" is literally what the leg-partition closure does. The decorrelation comes from randomness (entrywise n^{−1/2}), not from distance or a spectral gap.
2. **The gain mode is heredity of the scale plus the expander's Perron mode** (interpretation). Positive homogeneity carries a scale fluctuation forward exactly, the hereditary half. On the doubled→doubled channel the expander transmits only the norm-coupling (trace) part, the decorrelation half: everything else on that channel is averaged down to O(n^{−1/2}).

   Together they explain why the *amplitude* of the scale-mixture shape survives transport. They do not explain the rest:
   - the shape's var_i profile arrives through the distinct→doubled channel (§6.5);
   - most of the amplitude is fed afresh at each layer by non-Gaussian structure (§7.4);
   - at n = 1024 the shape holds only 25–39 % of the (2,1,1) slice (§7.3).

   This is the single "same essence" that pays in the competition, and it pays mainly at the D21 interface, where the κ3 shape still holds 72–86 % of the target at n = 1024 (§7.3).

**False.**
1. **"The Markov property is an expansion property."** The face law fails positivity (R4, derived). At the small widths where it could be checked it also has no safe symbol, and its Markov graph is complete or nearly so (R4, measured). Its accuracy-relevant content is not pairwise (R2). A pairwise Gibbs closure is biased: by up to 4/π for equal thresholds, without bound and possibly in sign for mixed ones (R3).
2. **"Approximate Markov ⇒ a forgettable past."** Old content does decay, at ≈ 0.15–0.3 per layer of age, by the gates' discount rather than a gap. But its horizon (≈ 8–16 layers) is comparable to the depth: at L = 16 only the one or two oldest sources fall below the bar (R11).
3. **"The time-averaged detailed-balanced Lindbladian gives a recovery map for erased content."** The estimator has no dynamics, no stationary state and no detailed balance. Part of the erased content *is* recoverable, but by a symmetry, not by detailed balance: the dilation carries the gain-template component unchanged (§5.2 (a)). The rest has Bayes-optimal recovery zero under H′ (R9).
4. **"Spectral independence certifies the estimator"** (R6).

**What the Lindbladian's "subtle essence" turns out to be, for the competition: *pair, don't recover*.** CR never needs the generator to mix. The defect is paired against a smooth question ([A] §0.4). Its commutative analogue (ANALOGY, R10) is the Duhamel telescoping that the heisenberg design already uses ([FU] unlock 36). Its compressions meet the same free-probability law (R10). Its gain is the bilinear error.

**NCG's part in the competition.**
- *Literal:* free probability governs the propagators' mode counts, and the annealed layer is a conditional expectation. Both are offline calculators.
- *Speculation:* the dilation as a central charge (§7.7). The Lebowitz–Percus–Verlet reading behind it is an ANALOGY.
- *None:* spectral triples, Connes and Carlen–Maas metrics, and Lindbladian dynamics have no role in an estimator that only pushes forward and pairs. [A] §0.6 and [B] place NCG's real weight in the theory: modular theory, sufficiency, derivations.

---

## 9. Experiments: predictions, falsifiers, tools, costs

| id | experiment | tool | prediction (label) | falsified if | cost |
|---|---|---|---|---|---|
| E1 | Best family-rank-1 u⊗M of the (2,1,1) slice at n = 1024 by streaming power iteration (u-step E[z̃_i²·z̃ᵀMz̃] minus Wick terms; M-step: the weighted Gram Z̃ᵀdiag(Σ_i u_iz̃_i²)Z̃ minus Wick terms), plus c·R and the flat variant, each scored by D21 error in the oracle's closure | `stream_oracle.py` (add two O(Nn²) contractions) | At layers ≥ 6 the best rank 1 closes ≤ 40 % of the κ4 gap (a rank failure); c·R closes ≈ the regeneration's 20 %; u ∥ var at cos ≥ 0.95 (CONJECTURE) | rank 1 closes ≥ 60 %: the non-gain remainder has one dominant mode and a cheap carrier exists | ≈ 10 iterations × O(Nn²) at N = 32k: about a minute per layer |
| E2 | Old-pool D21 share by age at n = 1024: propagate x_{t+1} = W_{t+1}ᵀD_t x_t from each source's centred activations and read the (2,1) slices | `stream_oracle.py` | The per-age profile of n = 128 (A1, t = 15, AD: 0.13 at age 7 falling to 0.004 at age 15; ages ≥ 14 together 0.016), so that only the 1–2 oldest sources can be dropped (CONJECTURE, §5.4) | the pool of ages ≥ 10 falls below 0.022 (forgetting becomes a lever), or the profile differs from n = 128 by more than 2× at ages 7–12 | 1–3 CPU-h at N = 32k |
| E3 | Sign test of absorption coefficients for the non-T residual (old pool minus its T-component) over 6 MLPs at n = 128 | `moment_atlas_np.py --k3`, `absorb.py`, `tracker.py` | Per-layer coefficient signs agree across MLPs at chance level (CONJECTURE from H′) | sign agreement > 80 % on most layers: a learnable conditional mean exists beyond T | 2–3 h |
| E4 | Kept energy for m scrambled legs: e(k)^m | `transfer_spectrum.py` (ensemble formula, order m) | All-distinct κ4 needs more modes than κ3 (DERIVED, §3.1); the (2,1,1) slice behaves as m = 2 (CONJECTURE, §3.1) | m = 4 needs fewer modes than m = 3 (the derivation's hypothesis fails), or the slice needs as many modes as m = 3 | minutes |
| E5 | gate_GG rebuilt from bivariate cumulants (μ, C, κ3_aab, κ3_abb, κ4_aabb, κ4_aaab) by bivariate Edgeworth | atlases | Agreement to the Edgeworth truncation; no trivariate input changes it (THEOREM R2; consistency check for the faces design) | residual far above the truncation estimate | minutes |
| E6 | Second-order closure residual by diagram class at n = 128 and 256 | `closure2.py` (theory stream) | Short loops (triangles) dominate (ANALOGY of §6.2) | long loops dominate | not priced |
| E7 | Gain template at n = 1024 on 3+ more networks; c·R scored at the D21 level inside the oracle | `gainshare.py` (scratch; formulas in the Appendix), `stream_oracle.py` | share_R ≈ 0.25–0.4 at layers 2–15; c·n/2 ≈ 2.4 (l = 2) to ≈ 10–14 (l = 15) within ±40 %; c·R closes ≈ 13–22 % of the κ4 gap (1 − √(1 − share_R), with share_R measured for 2 networks, C9) | share_R ≥ 0.6 at layers ≥ 6: then the gain mode would be the whole story | ≈ 12 min per network on 3 cores |
| E8 | Final-layer cost per unit D21 error: dropped old content against injected closure error of equal ‖·‖ | `chain128` teacher forcing | Same slope, 4.2e-6·ε², within 30 % (CONJECTURE, Haar pairing, §5.3) | old content costs ≤ 0.5× per unit: the pairing helps | a few hours at n = 128 |
| E9 | Old pool minus the gain template at n = 256 (and 1024 by streaming): T-share, then the old-content carriers run on the residual | `tracker.py`, `carriers.py`, `propproj.py` on a `--k3` atlas | no-AD old pool ≥ 70 % along T at layers ≥ 8; the residual needs ε_own ≈ 6.5–10 %, ≈ 9–14 u per layer (CONJECTURE, §7.5) | T-share < 0.5 at n = 256, or the residual needs ∝ n modes at 2 % | n = 256 atlas ≈ 1–2 h; carriers minutes |
| E10 | A scalar closure for the gain: transport of the (2,2)-trace through the doubled-leg flat mode plus first-order gate births with κ3 and κ4 inputs, compared with measured gex | `gainrec.py` (scratch) | Recovers ≥ 80 % of gex at layer 15, where the Gaussian-input recursion gives 24–47 % (CONJECTURE, §7.4) | < 60 % | an hour of derivation, minutes to run |

The cheap decisive pair is E1 + E7 for the κ4 slice and E9 for old content. E1 answers [CP] §6b's open question, a rank failure or a wrong-M failure, at the real width. E9 decides whether §7.5's template moves the old-content problem into the budget; §7.5's arithmetic says it does not by itself (9–14 u per layer against 3–4).

---

## 10. Checks (C1–C16): what was run, and the numbers

All runs are in the scratchpad (`csynth/`), not committed. Width-128 checks use the coordinating session's atlases: two MLPs, each with two independent sample atlases (sample seeds 3 and 4, N = 5·10⁵ each, with `pre_M211`). The other checks use fresh He MLPs generated with numpy. Fit and score always use different samples. The review's re-runs C13–C16 live in a separate scratchpad (`critic1C/`), also not committed; their formulas are in the Appendix (A4, A6–A8).

| check | what | result | § |
|---|---|---|---|
| C1 | ‖E‖ and ‖W∘W∘W‖ against 4√2/√n and 2√120/n, n = 128–1024 | 0.1767 vs 0.1768 (n = 1024); tripled 0.0222 vs 0.0214 | 6.4 |
| C2 | Probit three-gate law (Gauss–Legendre) vs pairwise max-ent | κ3_ν/κ3 = 1.21–1.28; leading factor 1.264–1.093 at t = 0.5–3; three-body Möbius O(ρ²); mixed thresholds in C13 | 4.2 |
| C3 | Orthants met by the row cone of W (n = 8–14) | 105/256, 317/1024, 1375/4096, 4606/16384 | 4.3 |
| C4 | Top Cor(g) mode vs spike image (n = 128, layers 1–4) | overlap 0.44, 0.39, 0.39, 0.49 (random 0.01); 82–99 % of η₀ | 3.3 |
| C5 | Transport channel fractions (n = 128) | tripled → diagonal 0.36 → 0.03–0.12 (the layer-1 value is contradicted by C15); D21 doubled→doubled 0.52–0.81, flat to 0.19–0.29 | 6.4 |
| C6 | Flat (2,1,1) carrier pilot | refuted at layers ≥ 5: worse than the published form, gap closed −2.05 to +0.33 | 6.5 |
| C7 | u_Creg vs var; slice energy in span(R) (atlases) | cos 0.91–0.997; share 0.27 → 0.94 (MLP 0), 0.30 → 0.85 (MLP 1) | 7.2 |
| C8 | D21 gap closed by c·R vs u_iC_jk (atlases) | ratio 0.74–1.14, median ≈ 0.9 | 7.2 |
| C9 | Gain template across widths (Gram trick) | table of §7.3 | 7.3 |
| C10 | Gaussian-input gain recursion | 24–47 % of gex at layer 15 | 7.4 |
| C11 | κ3 template T: D21 energy share (atlases); Monte Carlo cost of its amplitude | 0.90–0.97 at layers 8–15 (n = 128); 0.72–0.86 at n = 1024 (layers 6–15); 1 % amplitude needs 5·10⁴–1.2·10⁵ samples at n = 1024 with the control variate, i.e. ≈ 100–230 u per layer | 7.2, 7.6 |
| C12 | Old pool along T (tracker, AD and no-AD); undiscounted accumulator | no-AD cos² 0.49–0.97; residual 0.13–0.24 (fit), 0.13–0.32 (accumulator); AD cos² ≤ 0.33 | 7.5 |
| C13 | Three-gate law vs pairwise max-ent, re-run and mixed thresholds (same quadrature) | re-run of C2: 1.12–1.28 at 14 of 15 points; t = (−1, 0.5, 0.5), equal ρ: ratio 11.4, 5.1, 2.95 at ρ = 0.02, 0.05, 0.1; t = (1, −0.48, −0.5): ratio −0.26 to −3.2 (wrong sign); f(t) monotone on (0, 8] | 4.2 |
| C14 | Conditional mutual information of gate pairs, and the safe-symbol test (n = 6, 8; 4·10⁶ inputs; within-context shuffle null) | live pairs conditionally dependent except one pair at layer 3 at each width; no safe symbol at layers ≥ 1 except n = 6, layer 2 (2 vacua, a dead gate present); the predicted empty orthant is never sampled | 4.3 |
| C15 | Tripled→diagonal channel ‖(W∘W∘W)ᵀD3(a)‖/‖D3(z')‖ on fresh networks (n = 128, 256, 512; 3·10⁵ samples, two replicas) | 0.60–0.67 at layer 1, 0.33–0.35 at layer 2, 0.18–0.22 at layer 3, 0.02–0.10 at depth, falling with n; rms D3(z)·n = 15–19 | 6.4 |
| C16 | Per-source tracker on MLP A1 (seed 770100; two regenerated atlases, N = 1.2·10⁵ each; noise-free cross products) | reproduces [OC]'s A1 shares to ±0.001; per-age shares at t = 15 (§5.4); old-source cosines +0.38 (AD), +0.66 (no-AD); fixed-source decay 0.1–0.4 per layer; no-AD cos²(old, T) 0.77–0.94; residual 0.11–0.25 (best scalar), 0.11–0.34 (accumulator); c_best/c_acc 0.74–1.0 | 5.4, 7.5 |

**Validation of the estimators.**
- At layer 0 (z₀ Gaussian), gainscan gives c4 = 0.00 and gex·n/2 = −0.00. gainrec gives Δ₁ = gex₁ to 1 %.
- The transport identity of the pilots agrees with the oracle to 2.6e-16, and with the atlas D21 to 3.8e-6.
- The κ4 amplitudes from the Gram trick (C9, n = 128) agree in size with the explicit-tensor values (C8) on different networks.

---

## 11. Messages to the fresh-slate streams and to the theory

- **faces.**
  1. Gates are 0-homogeneous, so a face-law state is blind to the input radius; and the per-direction gain varies inside each cone, so faces alone cannot hold it (R5). Carry the gain with the faces: a "face × gain" state.
  2. Pairwise face statistics are bivariate functionals (R2) and cannot be the state for trivariate content.
  3. Positivity fails at every layer ≥ 1, and at small widths the support has no safe symbol (R4), so the Hammersley–Clifford factorisation is not available. Where a factorisation exists, its graph is complete or nearly so.
  4. A pairwise max-ent completion biases three-gate cumulants: by up to 4/π for equal thresholds, without bound and possibly in sign for mixed ones (R3).
- **bethe.**
  - Dense layers: the expansion parameter is n^{−1/2} per edge (TAP / Plefka), and short loops dominate (§6.2).
  - The doubled-leg flat mode is the "cavity field" of the trace channel.
  - The (2,1,1) slice is plausibly covariance-type in transport (§3.1; a heuristic, tested by E4).
- **signings.**
  - W ↦ −W kills odd transported content in the annealed average (hypothesis H′, §5.2), and the quenched odd content is the sign problem.
  - The gain mode is even, so it survives every gauge or sign average. It is the part of κ3 and κ4 a signing design gets for free.
- **tropical.** The gain ‖M_h x̂‖ is a function on the fan, constant in the radial direction of each cone and varying across cones. Its variance across directions is the amplitude of the dominant (2,1,1) shape at moderate width and of most of the old pool's D21 (§7.5).
- **heisenberg.**
  - Duality is the commutative analogue of Chen–Rouzé's telescoping (an ANALOGY; the identity itself is elementary). Its compressed backward questions obey the same free-probability law, at least in their linear part, because J and Jᵀ share singular values (R10).
  - Do not expect mode savings from duality. Expect the bilinear error.
  - E8 tests whether old content pairs with the final readout like closure error does.
- **markov.**
  - Old content decays along depth at ≈ 0.15–0.3 per layer of age, by the gates' discount, not by a gap. The horizon (≈ 8–16 layers) is comparable to L = 16, so only the one or two oldest sources can be dropped (R11).
  - The gain-template part of erased content is recoverable by symmetry, i.e. the accumulator of §7.5. The rest has zero Bayes gain under H′ (R9); E3 tests H′.
  - In the published (no-AD) convention the old pool is mostly the gain template (§7.5). The CMI that matters is that of the non-T residual.
- **bench / scaffold.**
  - Two Stage-Q diagnostics cost one sampling pass, O(Nn²): gex per layer, and the template shares share_R and share3_T.
  - The amplitude of a known template is a matched filter: its sample count does not grow with the number of entries. Each sample still costs a forward pass, N/n u per layer, so at n = 1024 it is a diagnostic, not an estimator component (§7.6).
- **Prong 1 / prong 3.**
  - The dilation is a central charge only for positively homogeneous dictionaries. Keep it in the dictionary, not in Notes 1–2 (drag test, §7.7).
  - The doubled-leg expander and the trace identity are dictionary facts too. The abstract counterpart, a conditional expectation whose Perron mode is a central element, is a question for [A]'s stage, not answered here.

---

## Appendix: formulas for reproducing the checks

**A1. The (2,1,1) slice without n³ tensors.** z̃ is centred and Co = C with its diagonal zeroed. For distinct (i, j, k):
- *Regression vector on C_off* (over distinct j, k ≠ i): u_i·(‖Co‖² − 2((Co∘Co)1)_i) = E[z̃_i²(z̃ᵀCo z̃ − 2z̃_i(Co z̃)_i)] − var_i(⟨Co,Co⟩ − 2(Co∘Co)1)_i − 2(Co Co Co)_ii. The draft's denominator ‖Co‖² kept the pairs with j = i or k = i, which the numerator excludes; the difference is a relative O(1/n) for each i.
- *Contraction with R:* ⟨M4, R⟩_AD = E[(z̃²·var)·Q − 2Σ_i var_i z̃_i³(Co z̃)_i + 2(z̃²·(Co z̃)² − (z̃²)ᵀ(Co∘Co)z̃²)], with Q = z̃ᵀCo z̃.
- *Template inner product:* ⟨R_A, R_B⟩_AD = Σ_i var^A_i var^B_i(⟨Co^A,Co^B⟩ − 2(Co^A∘Co^B 1)_i) + 2Σ_i var^A_i(Co^B Co^A Co^B)_ii + 2Σ_i var^B_i(Co^A Co^B Co^A)_ii + 4Σ_i[(Σ_j g_ij)² − Σ_j g_ij²], where g = Co^A∘Co^B.
- *Slice norm:* ‖K‖² = ⟨M4_A, M4_B⟩ − ⟨M4_A, R_B⟩ − ⟨R_A, M4_B⟩ + ⟨R_A, R_B⟩, from independent batches, each with its own (var, C). Here ⟨M4_A, M4_B⟩_AD = mean_{s∈A, t∈B}[p₂p₁² − 2p₃p₁ − p₂² + 2p₄] with p_m = Σ_i(z̃_si z̃_ti)^m, which takes four N_A × n × N_B GEMMs.

**A2. The κ3 template.**
- T_ab = var_aμ_b + 2μ_aC_ab for a ≠ b, and D21_ab = E[z̃_a²z̃_b].
- c3 = ⟨D21, T⟩/‖T‖², estimated per sample as q = Σ_ab z̃_a²T_ab z̃_b. The control variate replaces z̃_a² by z̃_a² − var_a.

**A3. The gain excess.** gex = Var‖z̃‖²/(tr C)² − 2‖C‖²_F/(tr C)², computed in one pass from Σ‖z‖², Σ‖z‖⁴, Σ‖z‖²z, Σz and Σzzᵀ.

**A4. The three-gate law.**
- Exact P(g₁, g₂, g₃): nested Gauss–Legendre quadrature of the trivariate normal on the orthants shifted by t.
- Max-ent fit: Newton's method on the pairwise binary exponential family matching all one- and two-gate marginals.
- Möbius coefficient on the triple: Σ_{S⊆{1,2,3}}(−1)^{3−|S|} log P(1_S).
- Mixed thresholds (C13): equal correlations ρ on all three pairs; the same quadrature and fit.

**A5. The cone test.** Draw 8·10⁵ points a ≥ 0: half with |N(0,1)| entries, half also masked to random faces of ℝⁿ₊. Take sign(aW) and count the distinct patterns. Repeat with 4·10⁵ draws of x ~ N(0, I) for the Gaussian layer.

**A6. The tracker test.** Run `tracker.run_sources(W, lay, ad=…)` on a `--k3` atlas.
- The old pool is the sum of the sources with s < l − w.
- T_l comes from the atlas's (μ, var, C).
- The accumulator adds ⟨D21(source of age w+1), T_l⟩/‖T_l‖² once per source and never discounts.
- C16 uses two atlases of the same network from independent samples. Shares, cosines and pool sizes are cross products, e.g. √⟨d_A, d_B⟩/√⟨D_A, D_B⟩, so the sampling noise of each atlas cancels.

**A7. The tripled channel (C15).** For each layer l ≥ 1, D3(a_{l−1}) and D3(z_l) come from per-unit raw moments of two independent sample replicas. The channel is (W_l∘W_l∘W_l)ᵀD3(a_{l−1}) in the x@W convention, normalised by √⟨D3_A(z_l), D3_B(z_l)⟩.

**A8. Conditional independence and the safe symbol (C14).**
- I(g_i; g_j | rest) is the plug-in conditional mutual information over all contexts of the other n − 2 gates. The null permutes g_j within each context.
- A vacuum o passes the safe-symbol test if, for every sampled pattern x and every sub-mask B of x ⊕ o, the pattern x ⊕ B is also in the sampled support.
- The predicted empty orthant is 1[c < 0] with c = W_l^{−1}1 (x@W convention).

---

## Review log (critic 1: mathematical correctness, 2026-10-01)

Every DERIVED proof was re-derived, and every cited number was checked against its digest or stream report. The cheap checks were re-run and C13–C16 added. Corrections, most important first:

1. **Per-source old content and forgetting** (§5.4, R11, §0 3(d), §8, §11, E2). Three draft claims are refuted by the tracker on the very network cited (C16): "each of the nine old sources carries ≈ 0.314/√9 ≈ 0.10", "no exponential decay is visible across ages 2–15", and "dropping even the oldest single source costs ε ≈ 0.1".
   - The single-source shares fall from 0.13 at age 7 to 0.004 at age 15.
   - The old sources are positively correlated (cosines +0.38 AD, +0.66 no-AD), not orthogonal.
   - A fixed source decays at 0.1–0.4 per layer. The Lyapunov gap was the wrong rate; the gates' discount ([D-TS] B8) sets it.

   The verdict, DOA as a forgetting principle, survives in a weaker form: the horizon is ≈ 8–16 layers, comparable to L.
2. **Monte Carlo cost of the κ3 amplitude** (§7.6). 5·10⁴–1.2·10⁵ samples cost N/n u per dense layer, ≈ 100–230 u per layer with the contraction, not 2–4 u.
3. **The 4/π factor** (§4.2, R3, §0, §8, §11). The bound (1, 4/π] holds for equal thresholds, and whenever all t_b·r_ab·r_bc share a sign. For mixed thresholds the pairwise closure's error is unbounded and can flip sign (C13). That f is monotone between its proved limits is checked numerically, not proved.
4. **Completeness of the dependency graph** (§4.3, R4). The draft proved it from marginal correlation, but the Markov graph is defined by conditional dependence. It is now measured (C14) at n = 6 and 8 and labelled accordingly. The Hammersley–Clifford consequence is restated through the safe-symbol condition, since the proof does not need full positivity; the support fails that condition too at small n.
5. **The tripled channel** (§6.4, R1, §0 1). "Never inherited" is false. The O(1/n) norm of W∘W∘W is matched by an O(1/n) target, and the channel carries 0.60–0.67 of ‖D3(z₁)‖ at layer 1 at every width (analytic ≈ 0.6) and 0.02–0.10 at depth (C15).
6. **The residual factor** (§7.5, §7.6, §0 4, R17). "2.6–4.6×" matched neither method: at layers 10–15 it is 2.3–3.1× with the accumulator and 2.8–5.4× with the best scalar.
   - A third MLP (A1) replicates C12. Its accumulator overshoot is slightly larger at layers 12–15 (c_best/c_acc = 0.74–0.78).
   - The n = 1024 cost arithmetic now uses the residual that the accumulator actually leaves: 9–14 u per layer, not 4–9.
7. **Gap closed** (§7.6, E7). An energy share s closes ≈ 1 − √(1 − s) of the gap in norm, so 0.25–0.39 gives ≈ 13–22 %, not 20–35 %.
8. **Deflation** (§3.2, R8, §0 5). f = 0.84 is not the age-8 value: it is 0.65 for the full source, and ≤ 0.22 for the all-distinct source whose mode counts were used. The saving is ≈ 10 % at age 8 and ≈ 20–25 % at age 15.
9. **The θvvᵀ example** (§3.3). Cor − I = θvvᵀ has a non-zero diagonal. Corrected to Cor = I + θ(vvᵀ − I/n), realised by an exchangeable mixture, with η = θ(1 − 1/n).
10. **c4 ≤ gex** (§7.1 (iv), §7.3, §0 4). It is a property of the scale-mixture model, not an identity; a family of counterexamples is given.
11. **Evidence for H′** (§5.2 (b)). The weight flip also flips μ, which is part of the conditioning state, so it does not show E[Y | S] = 0. H′ is now marked as an assumption with weak evidence.
12. **Faces and the gain** (§4.4, R5, §0 3(a), §11). "Blind to the gain mode" conflated the input radius, to which faces are provably blind, with the per-direction gain, which is not a function of the faces but which they partly predict.
13. **Appendix A1.** The regression denominator must exclude the pairs with j = i or k = i, as the numerator does.

**Labels changed (inflations removed).**
- R10, §5.3: KNOWN-LINK → ANALOGY (CR ↔ Duhamel is an in-programme connection) + THEOREM (the elementary telescoping, now proved) + DERIVED (linear part only).
- §5.2 (a), R9: DERIVED → scale-mixture model + measured.
- R4, §4.3: completeness and the safe symbol, DERIVED + checked → measured at small n, CONJECTURE at n = 1024.
- §3.1, E4: "covariance-type, m = 2" → heuristic / CONJECTURE. Its stated mechanism, flattening by W∘W, fails at depth by §6.5.
- §3.2: "the spike is a BGN outlier", cited fact → the digest's interpretation.
- §7.7: the identification of the quenched gain mode with the Roberts–Yaida–Hanin vertex, KNOWN-LINK → ANALOGY.
- §0 6, §8: the dilation as a central charge, "Analogy" → SPECULATION, as §7.7 already had it.
- §6.3: DERIVED → measured fit + arithmetic. §6.1: KNOWN-LINK → KNOWN-LINK for the expansion + identification for the coefficients (B3 is off).
- §8: "spectral independence holds" → η₀ moderate, boundedness not established.

**Numbers updated to the final stream reports.**
- [O1024] now covers two MLPs: 0.74–1.07 %, 2.14–2.69 % and 2.4–3.2 %; about n^{−0.9}; B3 fitted at 0.15–0.40.
- [C128]'s κ4-closure numbers: 2.0–2.4e-5 with κ4 = 0, 1.4–1.7e-4 memoryless.
- [D-TS]'s PR law holds to 4 % only up to age 4–5, and deviates by up to 31 % at the ages where old content lives.

**Checked and confirmed.**
- §4.1: the counterexample densities.
- §4.2: the leading-order formulas and their limits.
- §4.3: the positivity proof; the predicted orthant is never sampled.
- §6.4: the norms and the expander arithmetic, re-run at four widths.
- §7.1 (i)–(iv) and the two trace identities.
- §7.3: the Gram/Möbius identity; A1's numerator and R-contraction.
- The cost arithmetic of §3.1 (54–86 u) and §3.2 (ε_tot ≥ 0.16); the ε(n) fit; the C8, C9 and C12 table readings.
- §5.3: the Duhamel identity.

---

## Sources

- **Bridge digests:** [arxiv-2609.38007.md](../../digests/bridges/arxiv-2609.38007.md) §8 and §12; [arxiv-2504.02208.md](../../digests/bridges/arxiv-2504.02208.md) §§1.3–3.7, 5.3, 6.1, 7.1; [hammersley-clifford.md](../../digests/bridges/hammersley-clifford.md) §§3.1, 4.1–4.3; [expanders.md](../../digests/bridges/expanders.md) §§3, 9.10, 11.3–11.7; [hdx-spectral-independence.md](../../digests/bridges/hdx-spectral-independence.md) §§3.3–3.5; [nc-dirichlet-lindblad.md](../../digests/bridges/nc-dirichlet-lindblad.md) §4.6; [transfer-spectrum-measurement.md](../../digests/bridges/transfer-spectrum-measurement.md) §§0, 2.1–2.7, 3.4–3.5, 4, 5 (B2, B4–B8); [chat-2609.38007-retrieval-status.md](../../digests/bridges/chat-2609.38007-retrieval-status.md) §1.
- **Programme:** [research-program.md](../../research-program.md); [conditional-arrow-algebra.md](../../conditional-arrow-algebra.md); [simplicial-complex-as-decomposition.md](../../simplicial-complex-as-decomposition.md); [local-to-global-unlocks.md](../../local-to-global-unlocks.md); [mlp-bridge.md](../../mlp-bridge.md).
- **Competition:** [competition-plan.md](../../competition-plan.md) §§0, 3.1, 6a, 6b, 7 and the revision log; [phase2-intel-2026-10-01.md](../../digests/phase2-intel-2026-10-01.md) (board numbers); [fresh-slate/BRIEF.md](../../fresh-slate/BRIEF.md) §§1, 3; [fresh-slate/foundations-unlocks.md](../../fresh-slate/foundations-unlocks.md) unlocks 1, 35 and 36 (cross-references only).
- **Stream reports:** [oracle1024](../oracle1024/REPORT.md), [old-content](../old-content/REPORT.md) (and `results/tracker_A1.txt`), [chain128](../chain128/REPORT.md), [coef-ensemble](../coef-ensemble/REPORT.md), [costmodel](../costmodel/REPORT.md), [submission](../submission/REPORT.md), [est-accuracy](../est-accuracy/REPORT.md) (stopped), [est-cost](../est-cost/REPORT.md) (stopped).
- **Tools used for the checks:** [experiments/oracle_k3.py](../../experiments/oracle_k3.py), [experiments/moment_atlas_np.py](../../experiments/moment_atlas_np.py) atlases, [old-content/tracker.py](../old-content/tracker.py), plus scratch scripts whose formulas are in the Appendix.
- **From memory, as marked:** Bai–Yin / Latała; Mehler; Leonov–Shiryaev; Kotecký–Preiss; TAP / Plefka; Lebowitz–Percus–Verlet (1967); Roberts–Yaida–Hanin (2022).
