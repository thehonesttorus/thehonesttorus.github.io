# Frontier phase: the public chain read in our theory, changes tested on the adjusted score

2 Oct 2026, coordinator H with child sessions F1 (cost toward the floor) and F2 (finite-state memory).

- **Base.** 504aldo V29 (MIT), `estimator_v29r3.py`.
- **Reward.** Dev adjusted = raw × max(0.1, steady C/B) on w1024_d16 (6 MLPs, truth noise subtracted, paired).
- **Base numbers.** Raw 1.7576e-8, C/B 0.2526, dev adjusted 4.44e-9. Its public LB is 5.40e-9; the leaders are at 1.6e-9.
- **Theory map.** `THEORY-MAP.md` §1–5.

## 1. What our theory says the frontier is (summary of THEORY-MAP)

**Mechanism.** V29 is an augmented K = 3 cumulant chain.
- The third cumulant is a crossed product along the age axis: young ages 1–4 are exact dense atoms, older ages are confined to transported shared bases.
- Its κ4 is the CAP (non-traceless) harmonic sector only, memoryless (λ C_off).

**Our results explain its hand-tuned design.**

| feature | explained by |
|---|---|
| dense young tier | G10 (fresh layers incompressible) and the Hadamard wall (a per-source frame of rank k costs ≈ 7k/n dense products against 4) |
| old-tier boundaries at ages 5 and 8 | Fibonacci blocks of the age tree (Pearson–Bellissard F; team E's Fibonacci/Penrose class) |
| old-tier ranks 384 and 224 | the Dixmier allocation c·n/a with c ≈ 2.2 averaged over those blocks (team G's score optimum) |

**Limitations (L1–L4).**
- L1: the old-tier cost is linear in the number of old sources (107 u).
- L2: κ4 is CAP-only and memoryless, with a fitted λ that cancels other errors.
- L3: the errors sit in the per-layer marginals.
- L4: the young tier costs 116 u.

## 2. Experiments and results (paired)

**E1. Knobs** (MLPs 0–2, `results/E1_table.md`).
- V29's own knobs sit at a local optimum.
- The best changes are thin legs 8/8 (−2.1 % adjusted) and the nested tier from age 6 (−0.8 %).
- Every age-axis reallocation tried is worse, by +1 to +28 %.
- Strassen L6 with 16-leaves: −5.4 % C/B at unchanged raw.

**F2. Finite-state deep memory** (`F2/REPORT.md`). Negative at depth 16.
- The cohort Tucker core is exact (checked to float32 rounding).
- But the deep content needs k ≈ 200: six deep sources span 224 dimensions at 99.98 %, and each mover's 99 % rank is 122–155.
- A dense k³ core pays only beyond ≈ 6k³/n² sources, i.e. ≈ 66 at k = 224, against at most 7 at depth 16.
- At k = 64 it is cheap (−4.8 % C/B) but raw is +26 % to 5.6× worse.
- The deep content is long-lived. Dropping ages ≥ 9 costs 7×, ≥ 12 costs 2×, ≥ 14 costs +8 %.
- By-product: the post-transport basis (raw −1.6 %, cost +1.4 %) is an enabler for rank cuts.

**F1. Toward the floor** (`F1/`).
- An incremental join (exact-core truncation of span(Qp) plus a residual sketch): raw −7 % on MLPs 0–2, −2.5 % on all 6, cost −0.8 %.
- Thin and feedback rank 8: −3.7 % C/B.
- Aggressive rank cuts (e.g. R_OLD 256 / R_OLD2 128: C/B 0.209) raise raw by 50–70 %. The V29 cost–raw frontier is steep.

**E3. Renormalised first-order diagrams** (accuracy at zero cost).

- *Hyperedge scale.* The newborn's [D21 ⊗ C] hyperedge (the D21-feedback thin legs) is carried at rank 16 with first-order coefficients. Scaling it by s:

  | s | 0.8 | 1.0 | 1.2 | 1.5 | 2.0 |
  |---|---|---|---|---|---|
  | raw change, MLPs 0–2 | +2.8 % | 0 | −2.5 % | **−4.0 %** | −2.3 % |

- *Held-out check at s = 1.5.* MLPs 3–5 −3.2 %, all three improving; −3.5 % on all 6.
- *Where the gain comes from.* The shallow layers: ×2 below layer 8 gives −3.6 %, ×2 from layer 8 on gives −0.1 %.
- *Energy matching alone* (no fitted constant, (‖D21‖/‖QB‖)^γ) recovers only −1.6 %. Most of the gain is a genuine coefficient renormalisation.
- *κ4 λ against the stronger hyperedge.* λ wants to fall: LAM 0.80 gives −1.3 % more on MLPs 0–2, but is neutral on all 6 (MLP 5 worse). V29's inflated λ (EscAI: ≈ 3× physical) was partly standing in for the under-weighted hyperedge.
- *Control.* Scaling the slice residual, which is an exact decomposition and not a diagram, is very harmful (+44 to +56 %).

## 3. Combined estimator (6 dev MLPs)

| estimator | raw (6-MLP mean) | steady C/B | dev adjusted | Δ vs V29 |
|---|---|---|---|---|
| V29r3 (base) | 1.7576e-8 | 0.2526 | 4.44e-9 | — |
| syn_rfb8: incj + R_FB 8 + hyperedge ×1.5 + LAM 0.80 + Strassen L6 | 1.7086e-8 (−2.8 %) | 0.2302 | 3.93e-9 | **−11.4 %** |
| **pkg/estimator.py**: incj + R_FB 8 + hyperedge ×1.5 + Strassen L6 (LAM 0.95, defaults baked in, no env) | **1.7005e-8 (−3.2 %)** | **0.2302** | **3.92e-9** | **−11.8 %** |
| pkg5/estimator.py: same with upstream Strassen L5 (residual-safe variant) | 1.7016e-8 (−3.2 %) | 0.2435 | 4.14e-9 | −6.7 % |
| **F1/est_F1_final.py**: incj (JS 32, 1 pass) + thin 8/8 + hyperedge ×1.5 + LAM 0.80 + Strassen L6 (F1's 6-MLP finalist) | **1.719e-8 (−2.2 %)** | **0.2213** | **3.80e-9** | **−14.4 %** (−13.6 % in F1's convention, C/B = max over non-first MLPs) |

Per-MLP raw change of `pkg/` against V29 (MLPs 0–5): −4.7 / −50 (MLP 1's raw is near zero, so its relative change is large) / −5.7 / −1.1 / +2.4 / +0.5 %.

**Grader-harness check, and a caveat.**
- `run_p.py` (`whest run --runner subprocess`, graded caps) on this 4-core container cannot reproduce the grader. On it even the base V29r3 fails all 6 MLPs: residual 0.46–0.48 s against the 0.4 s cap, two worker EOFs, walls 88–109 s. The submission stream had measured 0.29 s residual on an idle machine with its grader emulator.
- Paired against base on the same machine, `pkg/` (Strassen L6) raises the residual by ≈ 10–15 % (0.52–0.54 s on the MLPs that completed) and the wall to 105–112 s. Its first predict also timed out once.
- L6 adds many small ops, so it eats into the residual and wall margins. Treat `pkg/` as needing a grader-emulator check (`notes/streams/submission/scripts/grader_emul.py`) on an idle machine before any upload.
- `pkg5/` keeps upstream L5 and so V29's residual profile.

## 4. What this means

- The changes our theory licenses inside the V29 architecture are worth ≈ 10–15 % on the adjusted score. That would move V29's LB 5.40e-9 to roughly 4.7–4.8e-9.
- The rest of the gap to the leaders (≈ 3×) is architectural. The leaders carry old-source content at almost no cost and κ4 beyond a memoryless core.
- Two points from our theory bear on that gap:
  1. Per F2's measurement and the Hadamard wall, no transported-frame or core representation of the old content is cheap at depth 16. The old content is long-lived (1/a weight) and needs about 200 dimensions.
  2. The accuracy lever is the traceless (FREE) κ4 sector, concentrated in the youngest births.
- A young-only FREE-κ4 carrier, paid for by the cost savings above, is the next concrete design. Not attempted today.

## 5. Break-through attempts after the synthesis (2 Oct, afternoon)

- **Radial (homogeneity) factorisation.**
  - *The identity.* Bias-free ReLU nets satisfy F(x) = s·F(u), with s = ‖x‖/√n and u uniform on the √n-sphere, so E_x F = E[s]·E_u F exactly (E[s] = 1 − 2.44e-4 at n = 1024).
  - *Link to costate.* Costate's dilation template (2μ⊗S + diag(S)⊗μ) is exactly the third-cumulant signature of this radial scale mixture (derived in `radial/fc_sphere.py`).
  - *Test (in FC, MLPs 0–2).* The u-problem is seeded with the sphere's κ4, c = −2/(n+2) times the pairings, through FC's mean-field κ4 channel. **Negative** at full memory: raw 2.69e-8 → 4.99e-8. It reduces the dependence on old memory by 22 % (window 4: 7.9e-7 → 6.1e-7). Multiplying FC-x by E[s] again (double counting) gives 1e-7-level errors, confirming that the x-chain already carries the radial factor.
- **Systematic-bias calibration** (`diag/analyse.py`, best estimator, 6 MLPs). Negative. The coherent bias is ≈ 0.5 % of each layer's MSE. A left-one-out affine final-layer calibration changes raw by +0.2 %, and constant-bias removal by −0.4 %. The residual is per-neuron: errors in the per-neuron variance, κ3 and κ4 marginals (EscAI: oracle marginals give 20×).
- **Oracle harness ready.** `diag/estimator_oracle.py` teacher-forces chosen per-layer pre-activation marginals (μ, var, κ3, κ4) into the best estimator. A plumbing test at layer 0 reproduces the baseline: raw 1.8768e-8 against 1.8773e-8.
- **Blocked on data.** Attribution needs per-layer truth far more precise than our dev set's (N = 2e6). The N = 1e9 per-layer moments of the 1,000 public `full` networks (keenanpepper/arc-whestbench-p2-full1000-N1e9, plus the weights from aicrowd/arc-whestbench-public-2026) are on huggingface.co, which this environment's network policy denies.
