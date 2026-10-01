# Fresh-slate design brief — WhestBench Phase 2 (1 Oct 2026)

The user's instruction (1 Oct, verbatim): "The system for the competition should be completely start from 0 fresh slate no implicit bias about how it may look, and the profound key theoretical unlocks are the full driving force for its design which we can explore by a mixture of quicker and then more proper competition setting adjusted mse experimentation." Adapting existing algorithms "as an NCG add-on is always going to be heavily restricted and ruin our development of the abstract framework."

## 1. The task, exactly

- Input x ~ N(0, I_n), n = 1024. Network: a_0 = x, z_l = a_{l−1} W_l, a_l = ReLU(z_l), l = 1..L, L = 16. W_l has i.i.d. N(0, 2/n) entries (He initialisation), no biases; the estimator receives the weights as float32 (L, n, n) in the x @ W convention. The layers' weight matrices are drawn independently.
- Output: the per-neuron means E[a_l] for every layer, shape (L, n). The score uses the final layer only: mean over the evaluation MLPs of MSE_final × max(0.1, F/B), where F is the FLOPs the estimator used and B = 2^41. Truth is a 1e9-sample Monte Carlo bake.
- One forward pass costs about 2^25 FLOPs, so B buys about 65,536 forward passes and the 0.1 floor about 6,554. Unit used in our notes: 1 unit = one (1024 × 1024) @ (1024 × 1024) product = 2^31 FLOPs; B = 1024 units.
- Metering (flopscope 0.12.1, every array operation is billed): dense matmul 2mkn − mn; float64 costs 2× float32; exp/log/x**2 16 per element; norm.cdf 48 per element (float32); same-object Gram (X.T@X written as one operand) half price; Strassen–Winograd written as flopscope ops is permitted. A fair-accounting rule (16 Sep) says a score benefit must come from the estimation method; packing values is banned.
- Grader limits per MLP: 120 s wall, 0.4 s Python-side residual time (time outside flopscope ops; ≈ 0.02–0.04 ms per flopscope call), 8 GB RAM, 2 vCPU; setup() ≤ 5 s; the smoke test runs a 256-wide, 32-deep MLP; any failure zeroes the MLP with multiplier 1.0 (catastrophic).
- Where the bar is (1 Oct): best adjusted 1.6e-9 (raw 1.5e-8 at 0.11 B); best raw 1.14e-8 at 0.151 B. At the 0.1 floor, beating 1.6e-9 needs raw < 1.6e-8; raw 5e-9 at 0.1 B would score 5e-10. **Correction (18:55 UTC): the Gaussian (covariance) closure scores raw ≈ 4.3e-6 at n = 1024, measured on the bench set w1024_d16** (an earlier version of this brief said 4e-5, a factor-10 error); the target raw ≈ 1e-8 is therefore ≈ 400× below the Gaussian closure. Plain sampling at the 0.1 floor is ≈ 1.2e-5 raw. Submissions close 17 Oct 23:59 UTC.

## 2. Fresh-slate rules

1. Do not start from an existing estimator: moment/cumulant chains (with "sources", hubs, tiers), covariance propagation, Monte Carlo / QMC / cubature, the public 504aldo code, Phase 1 methods. They may appear only as baselines to compare against.
2. Start from a theoretical principle of the programme (below) and the problem's own mathematics, and let the estimator follow from it: its state, its per-layer operation, its error mechanism, its cost. State the principle precisely and say where the design uses it in an essential way.
3. Keep the abstract framework clean: the competition realisation must not drag the general theory into niche specifics; negative results are charged first to the dictionary (the way the framework is realised), not to the theory.
4. Be quantitative early: a design that cannot plausibly reach raw ≈ 1e-8 at ≈ 0.1 B at n = 1024 should be found out in hours, not days.

## 3. Facts about the object (measured; stated neutrally, not as algorithms)

- Positive homogeneity (no biases): f(tx) = t f(x) for t > 0, so E f = E‖x‖ · E_{S^{n−1}} f exactly, and every linear region of the network is a polyhedral cone through the origin (a "face").
- The layer weights are mutually independent, so W_{l+1} is independent of everything that shapes the law of a_l; but the scored quantity is quenched (a fixed network), and per-neuron means fluctuate around their ensemble values by O(n^{−1/2}), far above the ≈ 1.2e-4 rms accuracy the bar needs. Per-network structure must be resolved.
- What information suffices (team EscAI's oracle on 8 width-1024 networks): giving an analytic propagation the true per-neuron pre-activation variance at every layer cuts its error by 40 %; adding the true per-neuron third cumulants reaches raw 8.4e-9; adding the joint fourth cumulant reaches 1.17e-9; the readout from exact per-neuron moments has bias 1.2e-10. The per-neuron marginal laws of pre-activations, to about fourth order, at every layer, carry essentially all the information; they in turn depend on joint (pairwise and higher) structure of the previous layer.
- Joint structure across depth: about 40 % of the pairwise third-order structure feeding layer l+1 (the (2,1) slice κ3(z_a, z_a, z_b)) comes from content generated more than one layer earlier; it is essentially orthogonal to anything computable from the current layer's state, and no low-rank representation of it below ≈ 0.3 n modes exists (widths 64–1024). Old content concentrates on the top singular directions of the gated propagators prod diag(Φ) W by a free-probability law (participation ratio ≈ n / (2 · age)), plus one spike along the mean direction that emerges with depth.
- The joint gate (activation-pattern, "face") law of a layer has bounded unpinned spectral independence: η ≈ 1.7–5.6 at width 128 and 2–8 at width 1024.
- How non-Gaussian structure is generated: at each ReLU, the new joint structure is given, to leading order, by first-order "gate diagrams" (Gaussian Wick terms plus one cumulant hyperedge with covariance edges) with combinatorial coefficients; the residual of that first-order description falls as ≈ n^{−0.8} at fixed depth.
- Covariance spectra collapse with depth (effective rank falls layer by layer); a fraction of neurons are nearly always on or always off at depth (exactly linear or exactly zero there).
- Sampling is far behind: Monte Carlo variance ≈ 0.05 per sample; making a rule exact to Hermite degree 5 buys only ≈ 3×, degree 32 ≈ 10× (Phase 1 measurements).

## 4. Evaluation protocol (two stages)

- **Width scaling is pre-asymptotic at small widths (measured 18:55 UTC):** the Gaussian closure's error falls as n^-0.82 between widths 64 and 128 but as ≈ n^-2 between 128 and 1024, so a fit on widths 64–128 over-predicts its width-1024 error 12×. Fit on 256/512 and check directly at 1024 (bench set w1024_d16: numpy predictors are feasible at n = 1024 — `eval_q` on that set, or `run_p_inproc.py`).
- Stage Q (quick, hours): a numpy prototype of the estimator; widths 64, 128, 256 (512 if feasible), depth 16, ≥ 4 MLPs per width; truth from `whest dataset bake` (choose N so the truth noise avg_variance/N is well below the errors you measure, or subtract it). Report raw final-layer MSE (and all layers), its width scaling, the projected cost at n = 1024 in units (notes/streams/costmodel/cost.py has measured prices), and the projected adjusted MSE at n = 1024. The shared benchmark under notes/fresh-slate/bench/ (being built) will provide standard seeds, truth and an evaluator; use it when it appears, bake your own until then.
- Stage P (proper, days): a flopscope implementation at 1024 × 16 through `whest run` with the grader caps, scored against high-N truth; the scaffold under notes/fresh-slate/scaffold/ (being built) gives a grader-safe template.

## 5. Where the theory lives

- The programme: notes/research-program.md, notes/conditional-arrow-algebra.md, notes/local-to-global-unlocks.md, notes/simplicial-complex-as-decomposition.md, notes/mlp-bridge.md. Objects are abstract barycentres, not activation vectors; an arrow is a linear layer with the ReLU immediately before it.
- Bridges (Gibbs states, approximate Markov property, expanders, high-dimensional expanders, spectral independence, detailed-balance Lindbladians, noncommutative Dirichlet forms, transfer spectra): notes/digests/bridges/ and, when finished, notes/bridges-gibbs-markov-expanders-ncg.md.
- The user's input of 1 Oct on matchings, permanents, Bethe approximations, Kasteleyn signs, tropicalisation and cluster algebras: notes/fresh-slate/input-2026-10-01-matchings-bethe-cluster.md (unverified; verification in notes/fresh-slate/foundations.md when it appears).
- Measured facts in detail: notes/competition-plan.md §§3.1, 6b, 7; notes/streams/*/REPORT.md. Read them for facts, not for designs.

## 6. Deliverable of a design stream

`notes/fresh-slate/designs/<key>/`: DESIGN.md (principle → derivation → estimator: state, per-layer operations, cost per layer and total at n = 1024 in units, error mechanism and its predicted scaling), prototype code, Stage Q results (tables), a verdict (projected adjusted MSE at 1024 with its uncertainty), and the next experiments. Push early and often.
