# Quantitative judge: the fresh-slate designs and breakthrough proposals at n = 1024

*Judge (quantitative), 1 Oct 2026, ≈ 20:15 UTC. Read: BRIEF.md (with the 18:55 calibration correction), CONVERGENCE.md (19:47 version, with MKV-2), every designs/\*/DESIGN.md, RESULTS.md and VERDICT.md, the three breakthrough REPORT.md drafts (costate, interpolation, region), foundations-problem.md, foundations-unlocks.md, foundations-input-verified.md, foundations-facts.md, scaffold/KIT.md and RESIDUAL.md. Every per-network result file that exists was read (Appendix A), and every product count was taken from the prototype code (Appendix B). Nothing was re-run; all new numbers here come from arithmetic on stored results and on KIT.md prices. Units: 1 u = one dense 1024³ float32 product = 2³¹ FLOPs; B = 1024 u. Adjusted = raw × max(0.1, C/B). The bar is 1.6e-9 adjusted (raw 1.5e-8 at 0.11 B).*

---

## 0. Verdict

1. **Nothing is close.** Every design and breakthrough proposal measured at n = 1024 is 60–260× above the bar. The top cluster (faces win6, bethe edge1 and edge0, the costate scale-mode variants, faces w16, signings v1) sits at **0.95–1.24e-7 adjusted (central values)**, i.e. 60–77× the bar. Six networks cannot resolve the order inside that cluster.
2. **Width will not close the gap.** Every design's error ratio to the Gaussian closure is flat from n = 256–512 to 1024: heisenberg 10.1× / 11.9× / 9.7× at 256 / 512 / 1024, faces w16 9.9× / 11.9× / 12.7×, markov ≈ 10× at 256 and 1024, signings 0.24 of the closure from 128 on, bethe n^-2.0 from 64 to 1024. Each design removes a fixed fraction of a first-order (1/n) error. Measured as an amplitude, the residual left is 24–50 % of the closure's error. The bar needs 6 %; the public chain is at 7 %.
3. **Several stream cost claims do not survive the measured prices.** Strassen L5 costs 273 ms of backend time per product (KIT.md, batch of 4), so more than ≈ 330 products cannot all run at L5 inside a 90 s backend budget (≈ 440 with no margin at all for the rest of the 120 s wall cap). This makes the Strassen prices quoted for faces w16, heisenberg, markov and MKV-2 infeasible. MKV-2 also does ≈ 1,124 products, not the ≈ 900 its "≈ 500 u" implies. On the other side, three designs were priced dense and could use Strassen: faces win6, bethe edge1 and signings. Section 3 lists the corrections. Applying them reorders the middle of the table, but not the verdict.
4. **The full-history designs are excluded by cost alone.** At ≥ 0.37 B the bar needs raw ≤ 2.1–4.3e-9. That is 2.6–5.4× below the best raw any estimator has achieved (1.14e-8) and only 1.8–3.7× above EscAI's oracle with exact per-neuron moments to fourth order (1.17e-9). The designs concerned carry old content through O(L²) source-target pair loops: faces w16, heisenberg, markov full, MKV-2, region FC and costate full.
5. **The fresh designs are re-deriving the moment chain's mathematics.** Their state, interface, error law, generation rule and old-content transport are the published chain's (Section 6). They implement it at a lower rung of its own diagram ladder and at a higher cost. The best fresh full-history raw (MKV-2, 2.42e-7 at 0.75 B) lies inside the range of the published chain *with its old-source tier removed* (1.3–2.9e-7 at 0.15 B; foundations-facts F8.4). That six independent principles land on the same objects is a finding about the problem: the information the readout needs forces them (F3.2, F6.3, F7.1, F8.1, F10.1; region N1–N4).

---

## 1. Method

### 1.1 What counts as a measurement

- **Direct:** the shared bench set w1024_d16: 6 networks, truth N = 2e6, truth noise 3.6e-8 subtracted. Ratios between designs are taken paired, on the same networks.
- **Projected:** two proposals have a width-1024 result on MLP 0 only (costate gl variants, region FC). MLP 0 is the hardest network for the closure (5.00e-6, against a 4.10e-6 mean) and the easiest relative to it for first-order designs: its design/closure ratio is 0.68–0.81 of the six-network mean for faces w16, markov full and MKV-2. Projecting by the closure would therefore flatter these proposals. I project instead through a paired full-history anchor on the same network (costate "full" on MLP 0, 3.92e-7) scaled to the six-network first-order value (4.23e-7, heisenberg). The quoted range runs from the closure-based projection to +12 % above the anchor-based one.
- **Uncertainty:** ±1 s.e. across networks (5–10 % of the mean). The truth-noise subtraction is unbiased with s.e. ≈ 1e-9 on a 6-network mean (bench/eval_p.md), which is negligible at the 1e-7 level. Unpaired differences below ≈ 10–15 % are not resolved by 6 networks (F11.3). A factor 1.3 in adjusted MSE is about the current resolution.

### 1.2 Cost: units, wall time and calls

- **Units.** Dense products are counted from the code (Appendix B). Prices are from KIT.md:

  | engine | units per product | backend time per product |
  |---|---|---|
  | dense | 1.00 | 21 ms |
  | Strassen L1 | 0.877 | 23 ms |
  | Strassen L2 | 0.772 | 47 ms |
  | Strassen L3 (batch of 4) | 0.679 | 71 ms |
  | Strassen L4 (single product; no batch figure) | 0.610 | 158 ms |
  | Strassen L5 (batch of 4) | 0.545 | 273 ms |

  Elementwise work is added where it is not negligible: bethe's 40-node bivariate-normal table, and MKV-2's multi-operand einsums.
- **Wall-feasible price.** Start every product dense. Upgrade products one level at a time, taking the largest units saving per second first, until the backend time reaches **90 s**. That leaves ≈ 30 s of the 120 s cap for client/server round trips (1–3 ms per call, RESIDUAL.md), elementwise work and glue. The C/B range is quoted from this mix (low end) to L3-only pricing (≈ 30–60 s, the safe end). Cutting the budget to 60 s moves central values by 0–8 % (faces win6 1.03e-7, faces w16 1.29e-7, MKV-2 2.06e-7, markov full 3.09e-7) and changes no rank.
- **What the wall cap allows.**

  | price ceiling | products at most | backend time | raw needed to beat the bar |
  |---|---|---|---|
  | 0.10 B | 187 | 51 s | 1.6e-8 |
  | 0.15 B | 281 | — | 1.07e-8 |
  | 0.20 B | 365 | — | 8.0e-9 |
  | 0.25 B | 438 | — | 6.4e-9 |

  At most ≈ 330 products fit at L5 inside 90 s.
- **Calls (the 0.4 s residual cap).** Residual is ≈ 0.025 ms per call in-process, ≈ 0.035–0.04 ms with glue (V29: 13,121 calls → 0.46–0.52 s). RESIDUAL.md's 2× margin is ≈ 6–8k calls. A call estimate is (sequential product stages per layer) × (calls per family: 46 / 64 / 79 at L3 / L4 / L5) + elementwise calls. Section 2 gives per-design estimates, assuming the per-source elementwise work is stacked into (S, n, n) arrays. None of the prototypes does this stacking yet.

### 1.3 Score (0–10)

S = 10 · log₁₀(4.1e-7 / A) / log₁₀(4.1e-7 / 1.6e-9), where A is the central wall-feasible adjusted MSE.

- S = 0: the Gaussian closure at the 0.1 floor, which every stream gets for free.
- S = 10: the bar.
- References: the public chain V29 (5.4e-9) scores 7.8. The same chain without its old tier (2.0–4.4e-8) scores 4.0–5.4. Plain MC at the floor scores below 0.
- One score point is a factor 1.73 in adjusted MSE. A difference of 0.5 points (×1.3) is within the uncertainty.

---

## 2. Ranking at n = 1024

"amp" is the residual amplitude, √(raw / closure raw), with closure raw 4.10e-6 paired on the same 6 networks (4.30e-6 for the bench linearised closure, used by signings and interpolation). The bar needs amp ≈ 0.060.

| rank | design (variant) | raw at 1024 | basis | amp | products (code) | wall-feasible units (C/B range) | adjusted, central (range) | × bar | S |
|---|---|---|---|---|---|---|---|---|---|
| 1 | **faces win6** (renormalised facet births, window 6) | 4.45e-7 ± 0.33e-7 | direct, 6 | 0.33 | 385 | 218 u (0.21–0.26) | **9.5e-8** (8.8e-8 – 1.23e-7) | 60 | 2.6 |
| 2 | **bethe edge1** (+ age-1 old triples; 6-product old step of DESIGN §3) | 8.14e-7 ± 0.51e-7 | direct, 6 | 0.45 | 189 (231 as coded) | 103 u + 14–30 u (0.11–0.16) | **9.9e-8** (8.7e-8 – 1.34e-7); as coded 1.18e-7 | 62 | 2.6 |
| 3 | **costate A2gl** (ages ≤ 2 exact + scale mode at law level) | ≈ 6.0e-7 (4.6–6.8e-7) | MLP 0 only (5.60e-7), projected | 0.38 | 309 | 168 u (0.16–0.21) | **9.9e-8** (7.6e-8 – 1.4e-7) | 62 | 2.6 |
| 4 | costate A1gl | ≈ 8.5e-7 (6.4–9.5e-7) | MLP 0 only (7.84e-7), projected | 0.46 | 218 | 119 u (0.12–0.15) | 1.0e-7 (7.5e-8 – 1.4e-7) | 62 | 2.5 |
| 5 | **bethe edge0** (exact Gaussian edges + (2,1) slice + fresh hubs) | 1.04e-6 ± 0.07e-6 | direct, 6 | 0.50 | 102 | 56 u + 14–30 u (floor) | **1.04e-7** (0.98–1.11e-7) | 65 | 2.5 |
| 6 | costate A3gsl (+ residual slice chain) | ≈ 4.7e-7 (3.6–5.3e-7) | MLP 0 only (4.35e-7), projected | 0.34 | 435 | 253 u (0.25–0.29) | 1.17e-7 (0.89–1.54e-7) | 73 | 2.3 |
| 7 | **faces w16** (all depths) | 3.24e-7 ± 0.31e-7 | direct, 6 | 0.28 | 615 | 381 u (0.37–0.41) | 1.21e-7 (1.10–1.46e-7) | 76 | 2.2 |
| 8 | **signings v1** (copula + carried slice + hyperedges, PQ = 1) | 1.04e-6 ± 0.07e-6 | direct, 6 | 0.49 | 217 | 118 u (0.12–0.15) | 1.24e-7 (1.15–1.64e-7) | 77 | 2.2 |
| 9 | **interpolation**: closure + χ factor + self-averaging readout (deg 2) | 1.41e-6 (leave-one-network-out) | direct, 6 | 0.57 | 31 | ≈ 17 u (floor) | 1.41e-7 (1.25–1.6e-7); in-chain with 2 training networks 1.60e-6 → 1.6e-7 | 88 | 1.9 |
| 9b | **interpolation TC** (trace channel t = Cov(\|ã\|², a): one n-vector per layer, no fitted constants) | 1.45e-6 | stated in commit 202246a's message only; no result file stored, network count not stated | 0.58 | 31 + O(n²) | floor | 1.45e-7 | 91 | 1.9 |
| 10 | **region FC**, second chaos only | ≈ 5.6e-7 (4.2–6.3e-7) | MLP 0 only (5.16e-7), projected | 0.37 | 495 | 295 u (0.29–0.33) | 1.6e-7 (1.2–2.1e-7) | 101 | 1.7 |
| 11 | **markov MKV-2** (curvature passage, full history) | 2.42e-7 ± 0.13e-7 | direct, 6 | 0.24 | ≈ 1,124 + 12–50 u einsum | 755 u (0.75–0.79) | 1.86e-7 (1.72–2.03e-7) | 116 | 1.4 |
| 12 | **heisenberg** HD, 8 sources (A = 7) | 4.90e-7 ± 0.21e-7 | direct, 6 | 0.35 | 659 | 413 u (0.41–0.44) | 2.0e-7 (1.9–2.25e-7) | 124 | 1.3 |
| 13 | region FC + first-order non-Gaussian slices | ≈ 3.8e-7 (2.9–4.3e-7) | MLP 0 only (3.51e-7), projected | 0.30 | 855 | 558 u (0.55–0.57) | 2.1e-7 (1.6–2.5e-7) | 130 | 1.2 |
| 14 | heisenberg HD all ages = costate "full" (exact first order) | 4.23e-7 ± 0.21e-7 (6); costate 4.03e-7 (MLPs 0–2) | direct | 0.32 | 855 | 558 u (0.55–0.57) | 2.3e-7 (2.2–2.5e-7) | 145 | 1.0 |
| 15 | faces mem21 (one-step facet tree + slice) | 2.59e-6 ± 0.14e-6 | direct, 6 | 0.79 | 90 | floor | 2.6e-7 | 162 | 0.8 |
| 16 | **markov** MKV full history | 4.18e-7 ± 0.19e-7 | direct, 6 | 0.32 | ≈ 1,034 | 689 u (0.68–0.69) | 2.8e-7 (2.7–3.0e-7) | 177 | 0.7 |
| 17 | markov window 1 / window 4 | 3.60e-6 / 1.26e-6 | direct, 6 | 0.94 / 0.55 | 186 / 505 | floor / 0.30 | 3.6e-7 / 3.7e-7 | ≈ 230 | 0.2 |
| 18 | **tropical** TCT-0 (= Gaussian closure with exact two-wall lift) | 4.10e-6 ± 0.33e-6 | direct, 6 | 1 | 31 | floor | 4.1e-7 | 256 | 0.0 |

**Note on interpolation TC (row 9b).** The only source for 1.45e-6 is the message of commit 202246a; no result file is stored. The trace channel explains 65 % (layer 2) to 93 % (layer 16) of the variance of the per-neuron κ3 at n = 1024 (MLP 0, N = 262k, results/t5). Its 1.45e-6 is ×2.8–3.0 below the closure with no fitted constant, at the closure's cost. That ties the fitted self-averaging readout (row 9) and is the best accuracy-per-cost add-on the fresh slate has produced. Before this number appeared, my projection was 1.8–3e-6.

**Not scored (no measurement at 1024):**
- **bethe v4/v5** (carrying the global order parameter). At w128 v4 is 10–30× worse than edge. The v5 oracle is worse than v5.

**Reference rows (not fresh):**

| reference | raw | cost | adjusted | S |
|---|---|---|---|---|
| public chain V29 | 2.13e-8 (amp 0.072) | 0.253 B | 5.4e-9 | 7.8 |
| same chain without its old-source tier | 1.3–2.9e-7 | 0.150 B | 2.0–4.4e-8 | 4.0–5.4 |
| leaders | 1.5e-8 / 1.14e-8 | 0.11 / 0.151 B | 1.6e-9 / 1.7e-9 | 10 |
| EscAI oracle (exact per-neuron v, κ3, κ4 at every layer) | 1.17e-9 | — | — | — |

**Calls.** Estimates for engineered ports, under the stacking assumption of §1.2:

| design | calls | estimated residual | notes |
|---|---|---|---|
| closure-level (tropical, interpolation) | 1–1.5k | ≤ 0.06 s | |
| costate gl | ≈ 2.1k (its own count) | ≈ 0.07 s | |
| bethe edge0 | 2–3k | 0.06–0.12 s | the prototype's 40-node Φ₂ loop alone is ≈ 400 calls per layer (6k) |
| heisenberg, region FC, costate full | 4–6k | 0.1–0.23 s | |
| faces w16, faces win6 | 5–7k at the L4/L5 mix | 0.13–0.28 s | over the 2× margin unless L3 families are used |
| bethe edge1 at L5 | ≈ 9–10k | 0.25–0.4 s | must drop to L3 for half its stages, which moves it to ≈ 1.2e-7 |
| signings | 8–11.5k as prototyped | 0.2–0.46 s | at risk: the K = 6 Mehler sums are ≈ 300 elementwise calls per layer unless stacked |
| markov | | | the K = 60 Mehler series in relu_cov_offdiag is ≈ 300 calls per layer as written |
| MKV-2 | | | the numpy prototype took 421–432 s per network (3.6× the wall cap) before any Strassen, because of its multi-operand einsums |

---

## 3. Corrections to the streams' and the dossier's cost and score claims

| claim (where) | correction | effect |
|---|---|---|
| faces w16 "≈ 1.1e-7 with Strassen (≈ 360 u)" (DESIGN §9, CONVERGENCE) | 615 products at L5 need ≈ 168 s of backend, over the cap. The cheapest 90 s mix (86 L3 + 529 L4) is 381 u. | 1.21e-7 (L3-only 1.33e-7) |
| faces win6 priced dense, "≈ 410 u, 1.8e-7" | 385 products. The 90 s mix (133 L4 + 252 L5) is 218 u. | **9.5e-8, the best central value of any design** (1.03e-7 at a 60 s budget) |
| heisenberg "485 u = 0.47 B, adjusted ≈ 2e-7"; "16 families × 79 calls" | 855 products, 7 per pair as the code does. At L5 that is 234 s of backend: infeasible (costate already noted this). Feasible: 558 u. A single family per layer is impossible because each pair's products form a dependent chain of ≥ 3–4 stages. | 2.3e-7 (A = 7: 2.0e-7) |
| markov full "550 u → 2.2e-7" | ≈ 1,034 products (cost_mkv.py also runs the two-site trees over l + 1 sources). At L5 that is 283 s: infeasible. Feasible: 689 u. | 2.8e-7 |
| MKV-2 "≈ 500 u → ≈ 1.2e-7" (markov DESIGN, CONVERGENCE) | Per source and layer the code does 3-channel transport (3 products) + link (1) + dCa (3 + 1) + two-site trees (2). That is ≈ 1,124 products, plus 12–50 u of multi-operand einsums. L5 throughout would take 307 s: infeasible. The numpy prototype runs 430 s per network. | 1.86e-7 (1.72–2.03e-7); ≈ 1.5e-7 even at wall-infeasible all-L5 |
| bethe edge0 "105 u → 1.07e-7" | The 40-node Gauss–Legendre Φ₂ on n² pairs plus the rest of the edge table is ≈ 1.4–2 u per layer, 14–30 u in all, and is not counted. It does not matter: the 102 products cost 56 u at L5, so the design stays at the floor. | 1.04e-7 holds |
| bethe edge1 "195 u → 1.5e-7" | Priced dense, and the code's old step is 9 products, not 6. At L5 (52–63 s) it is 103–126 u, plus 14–30 u elementwise. | central 0.99–1.18e-7 (range 0.87–1.58e-7) |
| signings v1 "190–260 u → 2.3e-7" | The PQ = 1 code is ≈ 217 products (15 per layer with the repeated (R0∘R0) @ (W∘W h h) product reused). At L5 (59 s) that is 118 u. | 1.24e-7 |
| CONVERGENCE "every design that carries old content … at 0.35–0.55 B" | Wall-feasible costs are 0.37–0.79 B. | — |
| CONVERGENCE "best adjusted ≈ 1.1e-7 (bethe; faces with Strassen)" | The top cluster is 0.95–1.24e-7 and includes faces win6, bethe edge1 and edge0, costate gl (projected), faces w16 and signings v1. Its internal order is not resolved. | — |

Two other numbers in circulation need context:
- **"MKV-2: best raw so far."** This is true and significant: paired against faces w16 on the same networks, MKV-2 is lower by 8.3e-8 ± 3.0e-8 (2.8σ). But best raw is not best adjusted. MKV-2 ranks 11th.
- **"bethe and signings both measure 1.04e-6."** These are different code paths with the same mean (signings per network: 1.17 / 1.21 / 0.92 / 1.10 / 1.06 / 0.78e-6). Both are a pair-belief state plus a carried (2,1) slice plus one-step-generated triples. Their agreement is a convergence result (Section 6), not a copying error.

---

## 4. Headroom: what each design would need, and whether its own data make that plausible

The requirement in each row is the raw needed at the design's own wall-feasible cost, then at its cheapest variant.

| design | needs | evidence (width law, ablations, oracles) | plausible? |
|---|---|---|---|
| **faces** (win6; w16) | win6: ≤ 7.5e-9 at 0.21 B, **60×** below its raw. A window-3 variant (≈ 0.13 B) needs ≤ 1.2e-8, but at w128 window 3 is 1.6× worse than window 6. w16 needs ≤ 4.3e-9 (76×). | The gain over the closure has saturated (×11.9 at 512, ×12.7 at 1024). The ablated levers are already inside (slice 5×, κ4 2.4× at w128). The one next-order attempt made it worse: (2,2)/(3,1) with a diagonal inner covariance gave 5.6e-5 against 3.9e-5. The curvature passage on κ3 is ≤ 3 %. External estimates of the next terms: MKV-2's curvature passage ×1.7 at 1024; heisenberg R3 ≤ ×7–20 at w64 with *true* κ3; chain128's "full first-order engine" ≈ ×10 above the Wick rung at w128. Best case ×10–20 gives raw 2–4e-8 at ≥ 0.21 B, i.e. ≥ 4e-9 adjusted. | No |
| **bethe** (edge0; edge1) | ≤ 1.6e-8 at the floor (**65×**); edge1 ≤ 1.3e-8 at 0.12 B (62×). | The width law is n^-2.0 from 64 to 1024, with the ratio to the closure fixed at 0.20–0.25. The largest measured lever of any cheap design is the true node (v, κ3, κ4): ×9–17 at w64. Its realisation fails so far (v4 is 10–30× worse at w128; the v5 oracle is worse). Even the full ×17 is short of ×62. All-distinct old content is not carried (signings L3 split: ≈ 70 % of that sector missed). | No, unless node κ4 lands *and* old content is carried at ≤ 0.12 B |
| **signings** v1 | ≤ 1.3e-8 at 0.12 B (77×). | Teacher forcing closes it. One copula + hyperedge step from an exact state gives 3.6e-5 / 2.55e-6 / 3.5e-7 at n = 128 / 256 / 512 (n^-3.8, then n^-2.9), extrapolating to 4–6e-8 at 1024. The local rule alone is 3–5× above what is needed, before ×10 of drift. | No (closed) |
| **costate** gl variants | ≤ 0.97–1.4e-8 at 0.12–0.17 B (**62×**). | The scale (dilation) mode carries most of the old content's effect on the readout at O(n²) per layer. At 1024 (MLP 0), A1 → A1gl gains ×2.7 at an equal product count, and A3gsl is 1.11× the exact first order at 51 % of its products. Ceiling: the exact first order (≈ 4.0e-7 raw), i.e. ≥ 4.8e-8 adjusted at 0.12 B, 30× above the bar, even with perfect old content. Second order is Θ(L³n³) by its own Theorem C5. | No for the bar. Yes as a cost device: near-exact-first-order accuracy (1.1–1.4× of exact) at ≈ 35–50 % of the products, if it holds on 6 networks. |
| **interpolation** (self-averaging readout; in-chain; χ) | ≤ 1.6e-8 at the floor (88×). | It measures ×3.05 on the closure (leave-one-network-out). TC reaches ×2.8–3.0 without fitted constants (commit message only). χ is an exact −18 % on 6 of 6 networks. On a K = 3 base, fitted per-layer corrections gained ×1.0–1.45 (504aldo F47/F68), and EscAI's 40-feature span explains 8.8 % (F10.5). Old third-order content has ensemble mean zero (F8.5), so no self-averaging map can reach it. | No as a design. As an additive tool, expect ≤ ×1.5 on a first-order base; the deciding run is on such a base. |
| **region** FC | Chaos-only: ≤ 5.5e-9 at 0.29 B (101×). With slices: ≤ 2.9e-9 (130×). | Its ledger gives a sufficient target. With true slices, first-order bivariate Edgeworth (D21 + κ4 (2,2)) leaves 4.4e-9 in the covariance channel, and the readout with true diagonal κ3, κ4 is at the noise floor. FC's D21 error is 14–25 % from layer 2 on; ≤ 4.5–6 % is needed (×16–25 in MSE). Cost must also fall 2–3× through compression of old sources. | Only as a chain with second-order-accurate slices at ≤ 0.15 B, which is the published chain's route (Section 6) |
| **heisenberg**, **markov** (full, MKV-2), **costate full**, region FC + slices | ≤ 2.1–3.9e-9 at 0.41–0.77 B (116–177×). | This is the structural exclusion: their cost alone requires raw below every achieved raw (Section 5). Heisenberg's second-order deciding experiment is negative (R9–R10: true κ4 on computed κ3 hurts). MKV-2's curvature passage bought ×1.7. Markov's own verdict: "no completion of this expansion reaches ≈ 1e-8 at ≤ 0.1 B". | No |
| **tropical** | ≤ 1.6e-8 at the floor (256×). | It equals the closure. The temperature route gives at most ×1.6–2.5 at 64–128, selected against the truth. | No (closed) |

**Width laws, all streams.** No design's residual is higher order in 1/n. The design/closure ratio is constant once n ≥ 256 (Section 0, item 2). The improvement the bar needs, from amp 0.24–0.50 to amp 0.06, cannot come from n = 1024 being "large". It must come from the accuracy of the first-order objects that are carried. In D21 terms: the ε² law, extra MSE ≈ 4.2e-6 ε² ≈ 0.95 × the closure MSE (F10.1; region measured 3.78e-6 ε² on MLP 0), says the designs' effective ε is 24–50 %. The published chain's measured 4–7 % matches its amp of 0.072, which supports this mapping. The bar needs ≈ 2–5 %.

---

## 5. The feasible region, and why the pair-loop designs are out

The bar is the curve raw × max(0.1, C/B) = 1.6e-9. At the wall-feasible Strassen mix, cost translates into a product budget per network:

| C/B | products | per layer | raw needed |
|---|---|---|---|
| 0.10 | ≤ 187 | ≈ 12 | ≤ 1.6e-8 |
| 0.15 | ≤ 281 | ≈ 18 | ≤ 1.07e-8 |
| 0.25 | ≤ 438 | ≈ 28 | ≤ 6.4e-9 |
| 0.37 | ≈ 615 (faces w16) | — | ≤ 4.3e-9 |
| 0.55 | 855 (exact first order) | — | ≤ 2.9e-9 |

The leaders sit at (raw 1.5e-8, 0.11 B) and (1.14e-8, 0.151 B), about 200–280 Strassen products.

- **Pair loops cost Θ(L²) products.** Carrying old content exactly costs ≈ 4–10 products per source-target pair over 105–120 pairs. That is 495–1,124 products, 0.29–0.79 B wall-feasible. The costate stream proves this is unavoidable for *exact* first order (Theorem C3: the double edge factors through a cut only via the n²-dimensional second chaos; Tucker truncation saves nothing above rank √(n/2) ≈ 22).
- **At ≥ 0.37 B the requirement is below everything achieved.** The bar then needs raw ≤ 4.3e-9, which is 2.6× below the leaders' best raw and within 3.7× of EscAI's exact-moments oracle. No first-order improvement listed in Section 4 comes within a factor 10 of that. Full-history pair loops are out *whatever their accuracy becomes*, unless their cost falls by ≥ 2.5×.
- **What remains is ≤ 187–281 products at raw ≤ 1.07–1.6e-8.** At comparable product counts, the fresh designs are at raw 8–10e-7: bethe edge1 189, costate A1gl 218, signings 217. That is 60× the leaders' raw at equal cost. This is the gap, measured at constant cost.

---

## 6. The fresh designs re-derive the moment chain's mathematics

The ruling forbids adapting the public chain, and this report does not propose it. The quantitative record nonetheless shows that every fresh principle, once realised at n = 1024, computes the chain's objects:

| object (the published chain's name, from competition-plan §3.1 and foundations-facts) | faces | bethe | signings | heisenberg | markov | region | costate |
|---|---|---|---|---|---|---|---|
| zeroth order: Gaussian covariance closure (K = 2 chain) | 'gauss', face-measure Mehler | Mehler-12 / exact edges | v0 Gaussian copula | reference chain ν | latent Gaussian field G | gclose | reference ν |
| the interface: D21 = κ3(z_a, z_a, z_b) entering Cov(a) through ½ E[δ(z_a) 1(z_b > 0)] D21 | "(2,1) slice into Cov(a), leading bivariate Edgeworth" | ½ K_ab L₋₁ L₀ in the pair map | carried D + hyperedge term | S(p, q) through E[δ(z_p) H(z_q)] | "(2,1)-slice correction of the field" | §3 D21 term | S_{s→k} into the injection |
| error law ≈ 4.2e-6 ε² (504aldo F71) | — | — | — | — | — | 3.78e-6 ε², "504aldo's ε² law, now derived" | — |
| generation: the first-order gate diagrams (Wick + one cumulant hyperedge; F6.3) | facet births on legs | hubs (local second-derivative vertex) | "trees + at most one hyperedge" | diagonal + star + coincident planes | single-site + two-site κ4 trees | second chaos + slices | source atoms |
| old content: sources transported by gated propagators W Φ W ⋯, young and old tiers | face-averaged arrows P_{s→l} | gated two-step propagator (age 1) | CP sources (v2, unstable) | ū_{l→k}, pulled back | P^{(m→l)}, K^{(m→l)}, window | Z_s(l) | U_{s→k}, A-truncation |
| next order: δ-insertion (curvature passage), κ4 (2,2) and (2,1,1) slices | (2,2)/(3,1) attempt failed; cp ≤ 3 % | Q edge belief (v4) | — | R3/R10: needed, not built | MKV-2: ×1.7 | N3: (2,2) necessary | C5: Θ(L³n³) |

The streams say so themselves:
- faces: "getting the folding right is exactly the (2,1)/(2,1,1) cumulant transport of the existing chains";
- heisenberg: "the computation coincides with the adjoint form of a second-order cumulant chain";
- signings: "that is a chain of cumulant sources … makes it no cheaper";
- tropical: "that collapses onto moment/diagram closures";
- region: FC "is the Heisenberg/Duhamel identity written in the chaos basis", with the chain's ε² law.

**Where on the chain's own ladder.**
- At w128, the fresh designs measure 3.1–6.8e-5: faces w16 3.93e-5, heisenberg 3.16e-5, costate full 3.13e-5, MKV-2 4.8e-5, markov 5.2e-5, bethe edge1 5.3e-5, signings 6.8e-5.
- That is the band of chain128's "leading Wick" rung (2.6e-5) and of its self-consistent chain with κ4 = 0 (2.0–2.4e-5).
- The "full first-order diagram engine" rung is 2.2e-6 (the chain128 figures are geometric means over 4 networks, with κ4 teacher-forced, so 2.2e-6 is the best case for that rung).
- At 1024, the best fresh full-history raw (2.42e-7) lies inside the range of the published chain with its old tier removed (1.3–2.9e-7 at 0.15 B), and is 11× the chain with it (2.13e-8 at 0.253 B).

So the fresh slate has re-derived the chain at its Wick rung, with exact but expensive old content, and without the κ4 (2,1,1)/(2,2) slices that the chain's next rungs need (F7.1; region N3).

**Why this is a finding about the problem, not a failure of the streams.**
- At 1024 the readout needs each layer's per-neuron (v, κ3, κ4) (F3.1–F3.2: given them, even an approximate joint state reaches 1.17e-9).
- Those marginals are manufactured through the pair slices (2,1), (3,1)/(2,2) and (2,1,1) (F6.1, F7.1).
- At n = 1024 these slices are generated by a parameter-free first-order gate-diagram law to ≈ 1 % of D21 (F6.3).
- About 40 % of them is old content, orthogonal to the present (F8.1, F8.3), with ensemble mean zero (F8.5).
- Any principle must therefore produce these objects. The principles differ in bookkeeping and in their exact error identities: Duhamel and the bilinear remainder, the E5 telescoping, the cut-CMI price list, the fresh-weight Frobenius lemma, and the dilation conservation law. They do not differ in what they compute.
- The fresh slate added two structural objects that the published chain does not use explicitly:
  - the **dilation/norm/trace sector** as an exactly closed O(n²)-per-layer carrier (costate C6; interpolation TC; the bethe order parameter; already listed as the strongest low-rank structure in F2.2–F2.4);
  - the **self-averaging split** of the error (interpolation).
- Neither has a measured end-to-end gain beyond ≈ ×3 at 1024: TC ×2.8–3.0 and the self-averaging readout ×3.05 on the closure; the scale mode ×2.7 over an A = 1 truncation.

**What the ruling leaves.** The ruling governs provenance; the numbers govern content. A fresh system reaches the bar only if it delivers all of the following at ≤ 0.10–0.15 B (≤ 187–281 Strassen products):
- the chain's interface accuracy: D21 relative error ≲ 2–5 % at every layer from 1 to 14 (F10.1, region N2);
- the (2,2) slice in Cov(a) and the (2,1,1) slice in the next D21 (region N3, F7.1);
- per-neuron κ3 and κ4 to ≲ 7 % (region N4).

The best fresh designs are at ε ≈ 25–30 % and 0.21–0.77 B.

---

## 7. Breakthrough proposals: what each changes quantitatively

- **costate** (theory and first 1024 runs).
  - C1: the co-state is n-dimensional and exactly closed in value.
  - C3 and C4: exact evaluation still costs Θ(L²n³).
  - C5: second order costs Θ(L³n³).
  - C6: the dilation sector is exactly closed (a conserved charge).
  - Measured at 1024: exact first order is 4.03e-7 at 855 products (MLPs 0–2). The scale mode lowers the error of an A = 1 truncation 2.7× at the same product count (MLP 0: 2.09e-6 → 7.84e-7).
  - Quantitative value: it brings first-order accuracy to within 1.1–1.4× of exact at ≈ 300–435 products instead of 855, i.e. from ≈ 0.55 B to 0.16–0.25 B wall-feasible. It does not lower the first-order floor (≈ 4e-7 raw). Projected adjusted ≈ 1.0e-7 (single network; replicate on all 6).
- **interpolation** (draft).
  - χ factor: an exact, free −18 % on the closure, 6 of 6 networks.
  - Self-averaging readout: ×3.05 leave-one-network-out at closure cost, giving 1.41e-7 adjusted. The in-chain version with 2 training networks gives 1.60e-6.
  - Trace channel (TC): an O(n²)-per-layer carrier of 65–93 % of the per-neuron κ3 variance. The end-to-end raw is 1.45e-6 at 1024 according to commit 202246a (no stored result file), i.e. 1.45e-7 adjusted at the floor.
  - Negative and precise: there is no Onsager memory, and covers converge to a tree that is 48× worse than the closure.
  - Quantitative value: cheap additive corrections. The facts say they shrink on a first-order base (F10.5, F8.5).
- **region** (draft). It contributes the most useful quantitative instrument of the campaign:
  - **transfer coefficients K(l)** whose priced local errors close the books on the closure to 6 % (5.29e-6 priced against 5.00e-6 measured, MLP 0);
  - the anatomy of the closure's error: 71 % off-diagonal non-Gaussian covariance, 28 % per-neuron readout, 0.3 % variance;
  - sufficiency of (D21, κ4 (2,2), diagonal κ3, κ4) at ≈ 5 % accuracy;
  - per-layer tolerances, with layers 6–12 carrying 60 % of the price.
  - Its prototype FC is HD in another basis (Appendix B: same 855 products). It replicates HD's accuracy on MLP 0 and has the same cost problem.
  - Every design should be scored against the N1–N5 ledger, layer by layer, before the next end-to-end run.

---

## 8. Measurement hygiene and the next quantitative checks

1. **Replicate every single-network 1024 number on all 6 networks:** costate A1gl, A2gl, A1gsl, A3gsl; region FC (both variants). Store the interpolation TC run (1.45e-6 is stated in a commit message only) together with the TC + χ and TC-on-a-first-order-base variants. Three of the top six central values rest on MLP 0 alone.
2. **Time the Strassen mix under the grader caps** (120 s wall, client/server transport) before quoting a Strassen-priced score for any design above ≈ 300 products. The Strassen prices quoted for four designs from three streams (faces w16, heisenberg, markov full, MKV-2) are wall-infeasible. faces win6 at 133 L4 + 252 L5 products uses ≈ 90 s of backend: verify it on `run_p.py`.
3. **Log products and calls from code** (as costate does with `nprod`) instead of quoting them from design text. Code and text differ by 1.2–1.4× (MKV-2 ≈ 1,124 vs ≈ 900, bethe edge1 231 vs 195, heisenberg 7 vs 5 products per pair in DESIGN v1).
4. **Use the one-step D21 ε at 1024 as the go/no-go per design.** Measure each design's one-step D21 error from the true state at 1024 against F6.3's 0.9 % for the first-order diagram law. This separates generation error from drift; the region ledger then prices drift per layer. This is cheaper than another end-to-end run, and it says which objects to fix.
5. **Truth near the bar.** N = 2e6 truth (noise 3.6e-8, 2.4× the bar's raw) is adequate for ranking designs at ≥ 1e-7 when paired. To certify a design near 1.5e-8, use N ≥ 1e7 (noise ≤ 7.5e-9) and ≥ 20 networks: the 6-network closure mean has ±8 % s.e. and the network-to-network sd is ≈ 20 % (F11.3).
6. **Numerics.** All fresh numbers are float64 numpy. Strassen L4/L5 in float32 has not been checked for the subtractive cancellations inside D21 assemblies (relative errors of 1e-5 would cost ≈ 1e-10 MSE, harmless now but not at the bar).

---

## Appendix A. Per-network raw MSE at n = 1024 (×1e-7; MLPs 0–5, truth noise subtracted)

| estimator | 0 | 1 | 2 | 3 | 4 | 5 | mean ± s.e. | source |
|---|---|---|---|---|---|---|---|---|
| Gaussian closure (Mehler) | 50.0 | 47.5 | 47.1 | 35.3 | 35.1 | 31.0 | 41.0 ± 3.3 | markov/results/bench1024_gauss.json (identical in faces) |
| bench linearised closure | 52.4 | 49.0 | 50.6 | 35.4 | 39.6 | 31.1 | 43.0 ± 3.6 | signings RESULTS |
| faces w16 | 3.32 | 3.84 | 2.37 | 2.84 | 4.36 | 2.73 | 3.24 ± 0.31 | faces/results/q_predict_hyb21_w16_w1024_d16.json |
| faces win6 | 4.59 | 4.92 | 4.37 | 3.28 | 5.63 | 3.92 | 4.45 ± 0.33 | q_predict_win6_w1024_d16.json |
| faces mem21 | 28.8 | 30.4 | 26.6 | 21.5 | 25.1 | 22.6 | 25.9 ± 1.4 | q_predict_mem21_w1024_d16.json |
| markov full | 3.92 | 3.98 | 4.02 | 4.42 | 5.02 | 3.72 | 4.18 ± 0.19 | bench1024_full.json |
| markov MKV-2 | 2.11 | 2.16 | 2.31 | 2.78 | 2.85 | 2.31 | 2.42 ± 0.13 | bench1024_c16.json |
| markov window 4 / window 1 | 13.9 / 43.4 | 14.2 / 38.9 | 13.0 / 39.9 | 10.9 / 32.8 | 13.5 / 32.7 | 10.0 / 28.2 | 12.6 / 36.0 | bench1024_w4.json, bench1024_w1.json |
| signings v1 | 11.7 | 12.1 | 9.2 | 11.0 | 10.6 | 7.8 | 10.4 ± 0.7 | signings RESULTS |
| costate full / A7 / A3 / A1 / A0 | 3.92 / 4.51 / 10.5 / 20.9 / 30.4 | 4.03 / 4.79 / 10.0 / 22.5 / 31.3 | 4.13 / 5.23 / 11.1 / 19.1 / 29.3 | — | — | — | 4.03 (3 networks) | breakthrough/costate/results/w1024_d16.jsonl |
| costate A1gl / A2gl / A1gsl / A3gsl | 7.84 / 5.60 / 6.59 / 4.35 | — | — | — | — | — | MLP 0 only | w1024_d16_gl.jsonl |
| region FC (chaos only / + Gaussian slices / + first-order slices) | 5.16 / 3.93 / 3.51 | — | — | — | — | — | MLP 0 only | region REPORT §5 |
| interpolation self-averaging (leave-one-out, deg 2) | — | — | — | — | — | — | 14.1 (per network 10.8–19.1) | interpolation REPORT T3 |
| interpolation TC | — | — | — | — | — | — | 14.5 | commit 202246a message only (no result file) |
| bethe edge0 / edge1; heisenberg all / A7 | — | — | — | — | — | — | 10.4 ± 0.7 / 8.1 ± 0.5; 4.23 ± 0.21 / 4.90 ± 0.21 | bethe results_q1024.json (no per-network data), heisenberg RESULTS R11 |

Paired differences (6 networks):

| comparison | difference (raw) |
|---|---|
| faces w16 − MKV-2 | +8.3e-8 ± 3.0e-8 |
| faces w16 − markov full | −9.4e-8 ± 2.4e-8 |
| faces win6 − markov full | +2.7e-8 ± 3.0e-8 (tie) |

Design/closure ratios:

| design | range across networks | geometric mean |
|---|---|---|
| faces w16 | 0.050–0.124 | 0.079 |
| MKV-2 | 0.042–0.081 | 0.060 |
| markov full | 0.078–0.143 | 0.103 |
| signings | 0.18–0.31 | 0.24 |

## Appendix B. Dense n³ products per network, read off the code (L = 16)

- **faces** (`fbt.py predict`, mode 'hyb').
  - Per layer l: the covariance sandwich (2), plus propagator updates P_{s→l} (l − 1), plus per source (legs K 1, D21 2, κ4 1).
  - w16: Σ(5l + 1) = 615. The x-space Jacobian product, which is dead code at w16, is excluded.
  - win6: 2 + min(l − 1, 5) + 4·min(l, 6), minus D21 at the last layer, gives ≈ 385.
  - mem21: 6 per layer, 90.
- **bethe** (`bethe.py`).
  - edge0: `contract` 5 + sandwich 2 per layer; the last layer is diagonal-only (3). Total ≈ 102.
  - edge1: + per layer a two-step propagator (1) + `contract` on P (5) + the pattern correction (3), i.e. 9 per layer, giving 231. The DESIGN §3 count is 6, giving 189.
  - Elementwise: Φ₂ by 40-node Gauss–Legendre ≈ 2,900 FLOPs per pair, plus ladders and the Stein table. That is ≈ 1.4–2 u per layer.
- **signings** (`copula1.py`, PQ = 1): Cz 2, X¹ 1, (R0∘R0)·(W²hh) 1 (computed 3× in the code, reusable), K21·W 1, Y 1, K22·W² 1, K31·W 1, Dn 7. That is 15 per layer; the last layer needs 7. Total ≈ 217.
- **heisenberg** (`hdpull.py`) = **costate full** = **region FC + slices**.
  - Per source-target pair: transport 1 (none for a source's first target), R 1, star 2, D2·U 1, plane and diagonal contractions 2.
  - That is 105 + 120 × 6 + 30 = 855. A = 7: 92 pairs, giving 659 (costate's `nprod` agrees).
  - region FC, chaos only: 105 + 120 × 3 + 30 = 495.
- **markov** (`mkv.py`, var21 = 'full').
  - Per source and transition: transport = PV, KV 2; PPU, PKU, KKU 3; the M assembly grouped 2. That is 7 × 105.
  - Two-site trees: 2 × 120 (l + 1 sources).
  - Per layer: field 2, new-source link 1, symmetrisation 1.
  - Total ≈ 1,034.
- **MKV-2** (`mkv2.py`).
  - Per source and transition: 3-channel transport 3 + link 1 + dCa (3 channel contractions + 1). That is 8 × 105.
  - Two-site trees: 2 × 120.
  - Per layer: 44.
  - Total ≈ 1,124. The multi-operand einsums (M3, M4, M3x, … on (3, S, n) arrays) add ≈ 100–700 FLOPs per (source, neuron) pair per layer, ≈ 12–50 u.
- **closure, tropical, interpolation:** 2 per layer + 1, total ≈ 31.
- **costate gl:** the prototype's own `nprod`: 218 / 268 / 309 / 435.

## Appendix C. Pricing and scoring formulas

- Units:
  - U = Σ_products price(level) + E.
  - Prices and per-product times are from KIT.md: dense 0.9995 u / 21 ms; L1 0.877 / 22.8; L2 0.772 / 47.3; L3 0.679 / 71.2 (batch 4); L4 0.610 / 158.5 (single); L5 0.545 / 273.3 (batch 4).
- Wall-feasible mix: start dense, upgrade dense→L1→…→L5 greedily within a 90 s backend budget. The upper end of the range is L3-only.
- Adjusted: A = raw × max(0.1, U/1024).
- Requirement at cost U: raw* = 1.6e-9 / max(0.1, U/1024).
- Score: S = 10 · log₁₀(4.1e-7/A) / log₁₀(4.1e-7/1.6e-9).
- Residual amplitude: amp = √(raw / closure raw), paired.
