# What six principle-driven designs converged on (1 Oct 2026, 19:15 UTC)

Six fresh-slate design streams started from six different principles (faces and barycentres; Bethe/cavity; matchings and signings; tropical geometry and temperature; the Heisenberg picture with Dirichlet forms; Markov networks and CMI). Within three hours they arrived, independently, at the same structure and the same two binding facts. That convergence is a finding about the problem, and it frames the next round. Sources: notes/fresh-slate/designs/*/DESIGN.md and RESULTS.md.

## The common structure

1. **A Gaussian (or Gaussian-copula) reference propagated exactly is the zeroth order of every design.** Tropical: every realisation of the zero-temperature skeleton reduces to it (no small parameter at L = 16; a 10²³-fold P − Q cancellation at n = 1024). Signings: loop diagrams with fresh-weight legs vanish for n ≥ 96, so one-step sums are tree sums around the Gaussian product. Bethe: the product of node beliefs is the tree reference. Heisenberg: the reference chain of the Duhamel identity. The small parameter of every correction is n^{-1/2}.
2. **Binding fact (i): the covariance must be propagated with its non-Gaussian corrections.** Heisenberg's oracle at width 64: true mean and covariance plus true κ3, κ4 → final MSE 2.8e-6 / 4.3e-6, while its own covariance chain with the same injection stays at 4e-5 – 2e-4; after first order, the residual is covariance drift caused by the joint fourth-cumulant slices κ4(pppq), κ4(ppqq). Markov: its first failure was variance drift, fixed only by putting non-Gaussian structure into the covariance propagation, and the two-site κ4 tree diagrams are leading order for κ4 because they add coherently. Markov's remaining error at depth sits in the layer means of nearly-always-on neurons.
3. **Binding fact (ii): third-order content generated many layers earlier matters at the readout.** Heisenberg: ages ≤ 1 give almost nothing, all ages give 18–23× at width 64, 8 ages keep 92 % of the gain. Markov: final MSE falls with the depth window only when all history is carried; the depth-cut CMI is large and the top n/4 transported directions recover only half of it. Faces: the final-layer third cumulant is born roughly uniformly over all 15 earlier layers; the most recent facets contribute about nothing net; old births do not accumulate because later facets fold them away; non-Gaussianity is carried by facets (codimension 1), not faces. The earlier moment-chain measurements say the same at n = 1024: old content is ≈ 40 % of the (2,1) slice, orthogonal to anything born at the current layer, not compressible below ≈ 0.3 n modes, and not absorbable into a renormalised local closure.
4. **Width laws are steep but pre-asymptotic below n ≈ 256.** First-order designs fall as n^{-2} to n^{-3.2} between 128 and 256; the Gaussian closure falls as n^{-0.82} between 64 and 128 but ≈ n^{-2} between 128 and 1024 (measured: raw 4.3e-6 at 1024 on the bench). Projections from small widths are unreliable; direct width-1024 runs on the bench decide.

## Where the designs stand (projections, being replaced by direct width-1024 runs)

| design | best variant | projected adjusted at 1024 | cost | status |
|---|---|---|---|---|
| heisenberg | first-order Heisenberg–Duhamel, all ages (n×n source-layer pull-back: diagonal + star + exact coincident planes) | **measured at 1024: raw 4.23e-7 ± 0.21e-7 (Gaussian closure 4.10e-6, gain 9.7×) → adjusted ≈ 2e-7**; 8 sources raw 4.90e-7 at ≈ 375 u | ≈ 485 u (0.47 B) | gain ≈ 10× at every width from 256 up (not growing; HD n^-1.86 vs reference n^-1.83); second-order deciding experiment NEGATIVE at 64/128 (true κ4 injection makes it worse) because the transported κ3 slices are 22–31 % wrong from the first transport step (decoupled-gate linear response misses the joint gate law and δ-insertions = markov's curvature passage); next: exact κ3 transport through one ReLU |
| markov | MKV full history (site expansion around a latent Gaussian field, all history carried) | **measured at 1024 (6 networks, paired): MKV-2 with the curvature passage raw 2.42e-7 ± 0.13e-7 at ≈ 500 u → adjusted ≈ 1.2e-7 (best raw of all designs so far); full history without it raw 4.18e-7 ± 0.19e-7 at 550 u → 2.2e-7**; window 4: raw 1.26e-6 at 283 u (3.5e-7); window 1: raw 3.60e-6 at 104 u; Gaussian closure 4.10e-6 | 0.10–0.54 B | closed: not competitive (≈ 140× above the bar); measured 256 → 1024 slope n^-1.8; accuracy gains paid back one-for-one in cost because history costs O(L²) products; reusable: two-site κ4 trees are leading order (none for κ3), and old content crossing a ReLU needs its curvature passage (a leading-order 16–20 % κ3 term, carried exactly by conditioning on each source neuron's own scalar latent) |
| signings | copula v1 (non-Gaussian marginals, latent R by Mehler inversion, carried (2,1) slice, first-order hyperedges) | **measured at 1024 (all 6 networks, paired): raw 1.04e-6 ± 0.07e-6 vs Gaussian closure 4.30e-6 (ratio ≈ 0.24 on every network) → adjusted ≈ 2.3e-7 (1.9–2.6e-7)** | 190–260 u | closed; one-step local rule ≈ 5× above target at 1024 (teacher forcing n^-3.8), drift adds ≈ 10×; old content not carried |
| faces | 'w16': facet births on renormalised legs over all depths, transported along face-averaged arrows, (2,1) slice into Cov(a), κ3/κ4 into the readout | **measured at 1024: raw 3.24e-7 ± 0.31e-7 (12.7× below the Gaussian closure's 4.10e-6 on the same 6 networks) → adjusted ≈ 2e-7 as prototyped (0.62 B), ≈ 1.1e-7 with Strassen (≈ 360 u)**; window 6: raw 4.45e-7 | 0.35–0.62 B | best raw so far; gain over Gaussian grows with width (2.2× at 64 → 12.7× at 1024), so its small-width projection was 50× too pessimistic; ablation at 128: (2,1) slice worth 5×, κ4 2.4×; the exact input-radius κ4 correction changes nothing (w512) |
| tropical | TCT-0 + T-extrapolation | ≈ 4e-7 | 0.02–0.04 B | closed: not competitive |
| bethe | v4: exact bivariate-Gaussian edges + (2,1) slice + fresh hubs + age-1 old triples + the global order parameter carried as a (2,2)-slice edge belief (per-pair variance scale mixture), node κ4 from two-site + tree terms | **measured at 1024 (6 networks, paired): raw 3.76e-7 ± 0.28e-7 (10.9× below Gaussian) at ≈ 245 u with Strassen L3 → adjusted ≈ 9e-8; without old triples raw 5.8e-7 at ≈ 153 u → adjusted 8.7e-8 (best adjusted so far)**; edge0 raw 1.04e-6 at 105 u | 0.15–0.35 B | at depth the global (2,2) spike is 95–97 % explained by ONE strongly non-Gaussian collective coordinate, the top principal coordinate g along the Perron/mean direction (κ3(g) ≈ 0.7, κ4(g) ≈ 1); early layers it is a norm/variance modulation; a collective-coordinate split (v5) double-counted as realised |

For comparison: the public moment chain reaches raw 2.1e-8 at 0.25 B (adjusted 5.4e-9); the leaders raw 1.5e-8 at 0.11 B (1.6e-9) and raw 1.14e-8 at 0.15 B; an oracle with exact per-neuron moments to fourth order at every layer reaches raw 1.17e-9. The best fresh design is 30–60× short of the bar.

## A third binding fact (bethe, 19:12 UTC)

The bethe oracle at width 64 finds that, after the pairs are handled, the binding object is the per-neuron (node) fourth cumulant, which is large and coherently positive because of a global order parameter: a rank-one spike in κ(z_a, z_a, z_b, z_b), seeded by the input radius and amplified with depth. Exact background: with no biases a_l(x) = R · b_l(θ) with R = ‖x‖ independent of the direction θ, so every joint moment factorises into a moment of R (chi law, known exactly) times a spherical-input moment; at depth the layer norm plays the same role. A global scalar mode is cheap to carry (one scalar law per layer), so this may be the first structural fact that a representation can exploit at O(n²) per layer.

## After the direct width-1024 runs (19:45 UTC)

Every design that carries old third-order content lands at raw ≈ 2.4–4.2e-7 at n = 1024 (markov MKV-2 with the curvature passage 2.42e-7, faces 3.24e-7, markov without it 4.18e-7, heisenberg 4.23e-7), about 10× below the Gaussian closure, at 0.35–0.55 B; the designs that do not carry it land at ≈ 1e-6 (bethe 1.04e-6 at 0.10 B, signings 1.04e-6). The best adjusted score is ≈ 1.1e-7 (bethe; faces with Strassen), ≈ 70× above the bar. Two streams independently located the same next-order defect: old content crossing a ReLU needs its curvature passage / δ-insertion terms (the joint gate law), a 16–31 % error in the transported κ3 slices that is harmless at first order and fatal at second order.

## The open question for the next round

Every principle so far has been realised as "a Gaussian reference plus carried corrections", and the corrections that matter are (i) non-Gaussian covariance propagation through each ReLU, driven by joint κ3 (2,1) and κ4 (3,1)/(2,2)/(2,1,1) slices, and (ii) third-order content of every age. Carrying (ii) explicitly costs O(L² n³) (pairs of source and target layers). The question is whether a theoretical unlock gives a representation in which (i) and (ii) are carried at O(L n³) with ≲ 10 dense products per layer, or makes part of them unnecessary, at the accuracy the bar needs (per-neuron mean errors ≈ 1e-4 rms at n = 1024).

## Round 2 results (20:45 UTC): accuracy is reached, the cost of memory binds

- **Region FC v2** (`breakthrough/region/REPORT.md`): a self-contained first-order Wiener-chaos estimator (two-time cross-covariances, every age, exact slices, an O(n²) mean-field κ4 recursion) scores **raw 3.03e-8 ± 0.35e-8 on all 6 bench networks** at n = 1024 (exact Gaussian closure 4.10e-6, ≈ 135× worse). That is within 2× of the leaders' raw. It costs ≈ 840 u dense (≈ 0.5 B with Strassen, if wall-feasible), so adjusted ≈ 1.5–2.5e-8, 10–15× above the bar.
- **Necessary conditions** (region §4), from transfer coefficients whose books close to 2–6 %: a dense covariance arrow to ≈ 0.3 % of ‖C_off‖_F at layers 6–12; D21 at every layer to ≲ 5–10 %, content of every age included (dropping ages > 4 costs 25× the bar); the joint κ4 only through two n-vectors per layer (its diagonal and the column means of the (2,2) slice) to ≈ 30 %; diagonal κ3 and κ4 for the readout. Not necessary: n² κ4 slices, the (3,1) slice, second-order Edgeworth on the covariance, coherent directions, early-layer accuracy.
- **The memory cost bound** (region §6; costate C3): in any atom-by-atom realisation the all-age D21 contraction costs ≥ 480 u (≥ 270 u with Strassen). Every cheap form measured at n = 1024 fails: row basis k = 256 (4× worse), shared Oseledets subspace q = n/4 (6×), Wick-atom pruning to n/4 (100×), slice-atom pruning to n/2 (2.5×). Atoms carry comparable energy.
- **Theory judge** (`judge/theory.md`): every full-history fresh design computes the first-order layer of the moment chain, and all share one approximation: transporting old content with mean gates is wrong at O(1) on coincident index patterns (exact first-order response Φ_j[Φ_i − E(a_i)φ_i/σ_i] against Φ_i²Φ_j: −27 % at t = 0, +22 % at t = 1). Through the ε² law that explains their common 3–4e-7. FC avoids it with exact slices. The one structural lead for cheap memory is the dilation/Perron sector (D1), unproven at n = 1024.
- **Quant judge** (`judge/quant.md`): the 120 s wall cap limits Strassen L5 to about 330 products in a 90 s backend budget (273 ms each), so several stream Strassen prices were infeasible.
- **chain128 final** (`notes/streams/chain128/REPORT.md`, width 128, κ4(z)_{aabc} at Monte Carlo truth): the κ3 closure ladder maps one-to-one onto the final MSE (K=2 2.1e-4; slices only 9.8e-5; leading Wick 2.6e-5; oracle first-order closure 4.3e-6; full first-order diagram engine 2.2e-6), with final MSE = a + k·ε(D21)² (corr 0.86–0.995). The (2,1,1) κ4 slice is the lever: exact 2.2e-6, rank-4 covariance-response family 3.7e-6, rank 1 4.8e-6, u_iC_jk regeneration 6.1e-6, zeroed 1.7e-5.
- **Interpolation final**: trace channel plus a scalar κ4 channel scores 1.05e-6 at n = 1024; not competitive; the no-Onsager theorem stands.

**The open problem is now one object: the all-age memory at about the cost of one layer.** The leaders carry it at about zero cost.

## Round 3 (from 20:45 UTC): the local-to-global essence programme

At the user's direction (`notes/essence/INPUT-2026-10-01-local-to-global.md`), six deep-dive teams study the "forget one piece, re-randomise, converge to global" process where it delivers large computational payoff, its noncommutative generalisation, and non-naive transfers to the memory problem: A permanents and polynomials, B pseudorandomness of layered computation, C free probability and criticality, D elimination and sparsification, E tilings and conic geometry, F the noncommutative framework and synthesis (`notes/essence/BRIEF.md`). In parallel: costate measures the dilation-sector share of the memory at n = 1024; region measures a CP-merged old tier and reprices FC under the wall-feasible Strassen mix; bethe tests a localization estimator conditioned on the collective coordinate.

## Round 3, first results (21:20 UTC)

Measured at n = 1024 (MLP 0 unless stated). Raw = final MSE minus truth noise.

| variant inside FC | raw | source |
|---|---|---|
| FC, all ages | 3.24e-8 (6 MLPs: 3.03e-8) | costate `results/fcdil_w1024.jsonl`, region |
| ages ≤ 2 only | 1.77e-6 | costate |
| ages ≤ 2 + dilation sector of older content (oracle amplitudes) | 7.67e-7 | costate |
| ages ≤ 4 only | 7.79e-7 | costate |
| ages ≤ 4 + dilation sector of older content (oracle / law-level) | 3.84e-7 / 3.71–4.35e-7 | costate |
| old content (age > 2) CP-merged to R = n atoms per layer (6 MLPs) | 2.02e-7; merge alone 3,510 u | region §8 |

- **The dilation lead (D1) fails at FC precision.** Carrying the conserved scale sector of the old content recovers about 2× of what dropping it costs, and leaves FC 12× worse at a window of 4. The coherent sector is real but carries only part of the readout-relevant memory.
- **CP merging is closed** on both accuracy (rank-limited: 8 ALS sweeps change nothing) and cost.
- **FC repriced under the wall-feasible Strassen mix:** ≈ 855 products, ≈ 557 u = 0.544 B, adjusted 1.65e-8 (region §8).
- **Essence team v0s agree on one verdict:** the local-to-global machinery explains why memory is expensive and measures each obstruction, but supplies no carrier. Team A: the nonnegative (Bethe, diagonal-pairing) sector holds about 80 % of memory energy at depth but leaves 0.47–0.70 relative error; additive signed sums defeat Clifford/quaternion sketches; the gate law is not Lorentzian. Team B: the hybrid argument holds in strong form (fresh-weight orthogonality) but none of the gaps does (effective λ ≈ 1, width ∝ n). Team C: content splits into a dilation Jordan sector (coherent across ages, odds law s/(1 − s) = 0.31·k) and a free sector (cross-age cosines ≤ 0.03, energy × g³ per step, PR = n/(2(age + 1)) to 1 % up to age 5); free energy beyond age 3/5/7 is 5.6/2.9/1.3 % of final D21 energy. Team E: an annealed conic "wedge calculus" reproduces the transfer coefficients and the old-content age law from 2-D conic geometry, but the binding object is quenched.
- **Interpolation (final):** in the Wiener chaos of the *weights*, 92–93 % of per-neuron κ3 variance at depth is chaos 1, σ²(Wᵀt)_p with t = Cov(|ã|², a), which has an O(n²) recursion; the (2,1) slice's chaos-1 part is the rank-one spike. Conditioning on the top principal coordinate at n = 1024 removes ≈ 78 % of excess kurtosis, but that coordinate is only mildly non-Gaussian there (κ3 0.18, κ4 0.04).
- **Reading.** Every coherent, annealed or symmetric sector (dilation, Bethe, trace, chaos-1-in-weights) is cheap and worth about 2×. The quenched, signed, incoherent free sector of old content is what costs, and it is needed to 5–10 %.

## Round 3, first positive carriers (21:30 UTC)

| carrier (inside FC, n = 1024) | raw, MLPs 0 / 1 / 2 | cost | source |
|---|---|---|---|
| FC, all ages exact | 3.24 / 1.81 / 3.03e-8 | ≈ 840 u dense, ≈ 557 u wall-feasible | region |
| ages 5–8 Tucker-compressed once at the cut in the static frame, R = n/4, then transported | 3.34 / 1.81 / 3.05e-8 (lossless) | 36 + 384 u once, then 1.5 + 16 u per layer | region §9 |
| same, R = n/8 | 4.00 / 2.66 / 3.98e-8 (+30 %) | 18 + 48 u once, then 0.75 + 2 u per layer | region §9 |
| same, rank-R CP instead of Tucker | 6.1e-8 (n/4), 8.8e-8 (n/8), MLP 0 | — | region §9 |
| old (age > 2) content's channel constant in the repeated index carried by one n-vector (≈ 0 u), remainder exact | 3.67 / 1.84 / – e-8 | O(n²) per layer | team B T2 |
| each source of age a read in the top 2n/a right singular directions of its own propagator (oracle) | 3.31 / 1.95 / – e-8 | ≈ 0.51 of FC (O(n³ log L)) | team D §3.7 |

- **Static once-per-bin compression works where per-layer moving frames failed** (6× worse at n/4) and where CP merging failed (6.7× worse at R = n). About three bins at R = n/8 carry all content older than four layers for ≈ 250 u; the young tier (ages 1–4, exact, ≈ 400 u) is now FC's dominant cost.
- **80 % of old energy is one channel.** Old D21 is 79–81 % constant in the repeated index (the dilation/norm channel), carried at O(n²) per layer. The irreducible object is the a-traceless remainder T° (≈ 45 % of the old tensor's norm), which loosens every merge tolerance ≈ 2.2×.
- **Reading (coordinator note 2).** Freezing the polarization per age block makes the memory's commutator effectively finite rank (R ≈ n/4–n/8 per block), although it is not compact for the global gauge polarization; the per-age resolution law 2n/a makes the age grading Dixmier-critical (total ∝ n ln L).
- **Open:** static bins for ages 3–4 and a cheaper core readout (region); T°-only bins (region); causal multiresolution with frozen dyadic frames (team D); a Hadamard-readout carrier of a few sandwich-transported symmetric matrices, the one untested class consistent with the leaders' cost (costate); width universality and the log-L lower bound (team G).

## Round 3: the causal age-multiresolution carrier passes, its constant does not (21:35 UTC)

| carrier inside FC, n = 1024 | raw, MLPs 0–5 | price | source |
|---|---|---|---|
| FC, all ages exact | 3.24, 1.81, 3.03, … ; mean 3.03e-8 | 0.85 B dense; 557 u wall-feasible | region |
| causal QR-transported frames, k = 2n/age (one conversion SVD at age 3) | 3.40, 1.98, 3.28, 2.35, 4.05, 3.97e-8; mean 3.17e-8 | 0.57 B dense; ≈ 0.53 B wall-feasible | team D §3.8 |
| same, k = 1.5n/age | 4.1e-8 (MLP 0) | 0.47 B dense; ≈ 0.43 B wall-feasible | team D §3.8 |
| Bentley–Saxe odometer, shared dyadic frames at 2n/(youngest age in block) | 3.26e-8, 1.84e-8 (MLPs 0, 1) | ≈ as above | team D §3.8 |
| odometer at n/(youngest age) | 6.8e-8 (fails) | — | team D §3.8 |

- **Note 2's carrier is real and causal:** lossless, O(n³ log L), and the dyadic odometer works.
- **Its constant is ≈ 3.5× too high for the leaders' price.** The atom form needs ≈ 9 products per unit of resolution (floor 7: four (n, n, k) contractions plus three full-coordinate materialisations forced by the Hadamard squares Y∘Y, Y∘Z, Z∘Z, Z∘T). Adjusted ≈ 1.7e-8, against FC's 2.6e-8, the public chain's 5.4e-9 and the leaders' 1.6e-9.
- **Where the bill now sits.** Young ages 1–2 (always full rank) are 34 % of it. The binding consumer is FC's slice correction dk21, which needs the full D21 in the neuron basis every layer: the coordinatewise ReLU forces the Hadamard materialisations (note 2 §2(b)).
- **Next lever:** a cheaper D21 readout from a factored source, and the young tier.

## Round 3: the combined static-bin design is closed on cost (21:40 UTC)

Region REPORT §10 (MLPs 0–2; FC alone mean 2.69e-8 on these three):

| design | raw mean | dense cost | wall-feasible |
|---|---|---|---|
| ages ≤ 4 exact + static Tucker bins every 2 layers, R = n/4 | 3.05e-8 (+13 %) | ≈ 2,190 u (cuts 1,230: core formation 192 per cut; core readout and transport ≈ 525) | ≥ 2 B (the Tucker core work is a 3-way einsum, not Strassen-able) |
| same, R = n/8 | 5.9e-8 | ≈ 863 u | ≈ 0.7–0.8 B |
| ages ≤ 3 exact, R = n/4 | 3.5e-8 | — | — |
| ages ≤ 2 exact, R = n/4 / n/8 | 5.6e-8 / 2.2e-7 | — | — |

- Static bins solve the *accuracy* of old content but not its cost; team D's causal per-age carrier (≈ 0.53 B wall-feasible, lossless) is the cheaper realisation of the same law.
- **The exact young tier (ages ≤ 4) costs ≈ 400 u dense, 2.6× the leaders' whole bill (≈ 154 u).** A design at ≤ 0.15 B needs a new representation of ages 1–4 as well as of old content. For comparison, the public chain's young tier costs 116 u at Strassen L5 (4 sources × 2.17 u per layer), about half of FC's per-pair constant.

## Round 3: synthesis and the deep end (22:00 UTC)

- **Team F (final):** every essence transfer factorises one conditional expectation, the Brauer/Weingarten fresh-weight average, with two exact sectors. CAP (O(n)-invariant: trace, dilation, Bethe, Perron; unipotent at He-criticality) must be pinned; FREE (through-strings) decays in energy ×g³ per step with rank n/(2(a + 1)) and tensorises (Theorem F1). The cap share agrees across teams (0.82 C, 0.82 F, ≈ 0.8 A, 0.75 D). **Cohort carrier S1:** one transported basis per age cohort, cores added without fit: k = 64 for all ages ≥ 8 gives 3.97e-7 vs 3.92e-7 inside costate's exact first order (MLP 0) at ≈ 1–2 u per layer, replacing ≈ 196 u of FC's pairs; k = 32 1.09×; ages ≥ 4 need k = 256. "Deep memory is cheap; the bill is ages 1–7."
- **Team G (v0):** quenched memory = −½[F, X]Ω for F = 2E_D − 1 (the masa expectation); a cut is a commuting square iff diag(P)W is a Bratteli matrix. T3: relative loss ≈ A·e^{−4.6c}, width-universal; PR(Z) = n/(2a) at 256 and 512; c for ≤ 10 % of FC raw: 1.24 / 1.39 / ≈ 1.58 at n = 256 / 512 / 1024. T4: a sharp conditional log bound, but energy decays geometrically, so 2n/a over-resolves old ages; leaner profiles under test.
- **Bethe (localization):** pinning the input-linear shadow of the Perron coordinate improves the Gaussian closure 4.10e-6 → 3.07e-6 (1.34×) but v4 only 1.03× at 3× cost; pinning does not remove the need for old content. Oracle inside v4 (MLP 0): true node v + κ3 at layers ≥ 6 → 5.2e-8 (from 4.3e-7). Localization is closed as a memory carrier.
- **Assembled next (region):** young ages 1–2 exact + team D's frozen frames for ages 3–7 + team F's cohort basis for ages ≥ 8, inside FC on all 6 MLPs, priced dense and wall-feasible.

## The keystone theorem (team G, 21:50 UTC)

**Theorem G10 (team G REPORT §5; proof given there).** If k ≥ 2 independent fresh layers form a λ-quantum expander on ℝⁿ (true with high probability for independent Haar layers, λ → √(2k − 1)/k; Hastings 2007, Pisier 2014), then every grading D of neuron space that each layer moves by at most one level has a counting function that grows exponentially until half the space. Hence no finitely summable grading of neuron space is uniformly compatible with the layers; at best θ-summability (the finite-n form of Connes' 1989 hypertrace obstruction).

- **Checked.** Team F verified the proof (the expander bound, the commutator bound and the shell rank) for orthogonal layers; the gated, non-unitary case is still open.
- **What it explains.** Every static carrier in neuron space failed for this reason: shared Oseledets subspaces, per-layer moving frames, fixed-rank bases.
- **Where compression is possible.** The only amenable direction is the age axis: a single shift ℤ, N(λ) = λ, Dixmier-critical. The working carriers (team D's frozen frames, region's frozen Tucker bins, team F's cohort bases, the dyadic odometer) all live there. Their algebra is a crossed product of an AF algebra of age blocks by ℤ, for the odometer the Bunce–Deddens algebra C(ℤ₂) ⋊ ℤ: the user's "shift on sequences with tail equivalence", exactly.
- **Measured companions (team G T3).** Relative loss ≈ A·e^{−4.6c}, width-universal at n = 256 and 512; the energy profile is subcritical (geometric), so the asymptotically optimal total resolution is O(n) in L, with no gain at L = 16.
- **Housekeeping.** The three background workflows launched before 20:00 (judge panel, foundations check, bridges-synthesis critic) stopped at ≈ 20:10 UTC, most likely in a session restart. Their finished outputs are in the repo (judge/quant.md, judge/theory.md, the foundations notes, bridges-synthesis A/B/C); the judge panel's final decision was never produced and has been overtaken by round 3.

## Depth 32 (team E G3, 22:00 UTC)

At 256 × 32 the causal carrier at k = 2n/a is lossless at every layer (±3 %), and its cost grows ×2.72 from depth 16 to 32 (predicted ×2.69 for the critical 1/a profile, ×2.27 for geometric; FC grows ×4.1). The memory's age-graded rank class matches the Fibonacci/Penrose quasicrystal (log at 99 %), not Thue–Morse (team E §6).

## Free products are two-dimensional (team B G1, 22:00 UTC)

Team B proved that for walks on the d-torus the λ-term rank decays as Δ^{−d/2} and is Dixmier-critical exactly at Pólya's recurrence dimension d = 2 (measured rank·Δ = 12,000 ± 6 % over Δ = 4 … 512 at N = 4096, as predicted), with a log₂ T-frame blockwise estimator where truncation needs ≈ N ln(1/ε) steps. Free variance is additive under free multiplicative convolution, so free products are universally "d_eff = 2" (Ginibre PR·(a + 1)/n = 0.96–1.00): the reason the network sits on the wall. In FC at n = 1024, r₉₀·a/n = 0.58–0.63 for ages 4–14.

## Rational memory at width 1024, and the cost floor (coordinator H, 23:00 UTC)

- **Causal merge at n = 1024** (ages 1–4 exact, R = n/2 merged histories, window 2): price-weighted error 0.69 % of D21 energy (8.3 % amplitude), against 0.18 % (4.3 %) at n = 256. The old-content residual grows slowly with merges at 1024 (≈ 1 % → 4.3 %) where it is flat at 256. This is at the edge of region's tolerance (H REPORT §10).
- **Same curve, not a new one.** At ≈ 4.5 n atoms per deep target the merge is a different point on the resolution–loss curve of team G's c·n/a frames (c = 1.5 ≈ 4.7 n, c = 2 ≈ 5.6 n), and it needs a fit the frames do not.
- **The young pairs set a cost floor.** Ages 1–2 exact are 29 pairs at 6–7 products, ≈ 125 u wall-feasible; with the covariance arrow that is ≈ 0.14 B before any older content. At FC's raw that floors adjusted at ≈ 4.2e-9 (public chain 5.4e-9, leaders 1.6e-9).
- **What would move the score** (H §11): a representation of ages 1–2 that avoids the per-source Hadamard materialisations, or a cut of the per-pair constant (7 → 4, slice legs only where their score weight is large); and FC's raw accuracy.
