# Stream old-content: a cheap carrier for the transported old-source content of κ3

Status (2026-10-01, evening): widths 64 and 128 done (sample truth and noise-free model chain, two MLPs, noise replicate);
width 256 model chain partly done, width-256 sample-truth atlas building. The verdict below is provisional on width 256.

## Question

The published K = 3 chain carries κ3 as per-layer sources, each transported with dense n×n legs; its old-source tier
(age ≥ 5, shared basis 384 → nested 224) costs 107 of 260 units (≈ 6.7 units/layer). Leaders run 0.11–0.16 B in total
(≈ 7–10 units/layer including everything). Is there a carrier of the old content's effect on
D21(l+1)_ab = κ3(z_a, z_a, z_b) with ε ≤ 2–3 % at a small fraction of that cost (n = 1024)?

## Method

- `tracker.py` — exact per-source decomposition of κ3(z_l) on a full third-moment atlas (width ≤ 256, dense n³).
  Two exact conventions: **AD** (this note's: O_l = AD(Φ³κ3(z_l)) is the old content, everything else is a birth) and
  **no-AD** (the published chain's: sources transported with the full Φ³, X_{s,t} = M_{s→t}^{⊗3} B_s exactly, slice
  corrections absorbed into each birth). Validation: Σ_s X_{s,l} reproduces the atlas κ3(z_l) to 1e-6–1e-5 (float32
  accumulation of the atlas) at every layer (`results/tracker_*.txt`, column 'valid'). A noise-free **model chain**
  (`--model`: births = first-order closure from the chain's own κ3(z), atlas slices) runs at any width from a pair atlas.
- `carriers.py`, `propproj.py` — carriers, each measured as ε_l = ‖D21(Ôld_l) − D21(Old_l)‖ / ‖D21(z_l)‖ (whole D21 as
  denominator, i.e. the ε of the chain's error law). "Old" = sources of age > w; w = 4 (the published young tier) unless
  stated. **dyn** = the carrier feeds its own approximation forward (compounding); **static** = re-fitted to the true
  Old_l each layer (oracle bound). Re-truncations inside dynamic carriers are oracle (dense SVD / HOSVD / ALS); the cost
  column prices a randomized factored version.
- Atlases: seeds 770100 (A, two independent N = 6e5 sample sets A1/A2) and 770101 (B1) at width 128; K64 (seed 770100,
  width 64, N = 6e5); pair atlases P64, P64b, P256 for the model chain; K256 (width 256, N = 2e5) building.

## Results

### Where the D21 content lives (tracker, width 128, A1)
- AD convention: age-1 source ≈ 0.75–0.93 of ‖D21‖; ages ≥ 2 together 0.37–0.63; ages ≥ 5 0.06–0.46 (grows with depth).
  No-AD convention: ages ≥ 5 carry 0.10–0.79 (the slice corrections make young and old partly cancel).
- Old sources' D21 contributions are nearly orthogonal to each other and to the young ones (regressions below).
- The dominant tensor mode of the old pool is the mean direction (cos(v₁, μ) ≈ 0.93, 60–80 % of tensor energy at depth)
  but it carries almost no D21: tensor energy is the wrong metric, only D21 error is reported below.
- Noise: A1 vs A2 (independent samples, same MLP) agree to ±0.005 in every ε; B1 (second MLP) within the same ranges.

### Carrier table (width 128 unless stated, no-AD, w = 4; ε range over layers 5–15; per-layer numbers in `results/`)

| carrier | variant | ε (D21, rel.) | est. cost at n = 1024, units/layer | file |
|---|---|---|---|---|
| none (drop old content) | — | 0.10–0.79 | 0 | carriers_A1_noad_w4 |
| (c) regression on O(n²) chain objects | C, C∘C, var/μ outer, slice transports and sandwiches of D21(z_{l−1}), young D21 (in-sample, per layer) | 0.075–0.45 (AD: 0.06–0.43) | ~6–10 | carriers_A1_{noad,ad}_w4 |
| (c) old pool on young sources' D21 | in-sample | 0.09–0.63 | 0 | carriers_A1_noad_w4 |
| (d) geometric tail | one more exact source, ages ≥ 6 = γ_l × it (+ young) | 0.055–0.52 | ~2–4 | carriers_A1_noad_w4 |
| (a) full symmetric mode family, r dense S_r transported as WᵀΦSΦW | r = 8 / 16 / 32, dyn | 0.028–0.074 / 0.014–0.027 / 0.003–0.005 | ≈ 2.5 r (40 at r = 16) | carriers_A1_noad_w4 |
| (e) mode family, directions fixed to top-r eigenvectors of C_l | r = 16 / 32, static | 0.020–0.036 / 0.007–0.015 | ≈ 2.5 r | carriers_A1_noad_w4 |
| (e) S_r fitted on C-features (C, C∘C, μμᵀ, …) | r = 16, static | 0.08–0.23 (worse than dropping at shallow layers) | ~0 | modesS_A1_noad_w4 |
| (e) low-rank mode family, each S_p rank q | (r, q) = (16, n/8) / (16, n/4) / (32, n/8), dyn | 0.017–0.031 / 0.014–0.027 / 0.006–0.012 | ≈ 4 r q / n (8 / 16 / 16) | modesSdyn_A1_noad_w4 |
| (e) **shared-basis family** (Tucker r × q × q, one U for all modes; q = n/4 ↔ 256 at 1024) | (16, n/4) / (32, n/4) / (32, n/2), dyn | 0.016–0.029 / 0.006–0.012 / 0.002–0.005 | ≈ 2 r q²/n² readout + ~3–6 merge (≈ 5–8 at (32, 256)) | modesU_A1_noad_w4 |
| (b/d) CP / one merged hub of R columns, legs transported, ALS merge | R = n/2 / n, dyn | 0.012–0.025 / 0.007–0.014 | 3 R/n transport+readout + ≥ 4.5 per ALS sweep (≈ 15–20 at R = n/2) | cpdyn_A1_noad_w4 |
| propagator projection (top-k left singular vectors of M_{s→t}, per source) | k = 16 / 32 / 64, static, pool | 0.05–0.096 / 0.014–0.027 / ≤ 0.004 | k ≈ 0.3 n ≈ 300 at 1024: the published shared basis (≈ 6.7) | propproj_A1_noad_w4 |
| any of the above with the AD convention | — | 2–5× worse at equal size (AD re-masking scatters content out of the propagator range) | — | carriers_A1_ad_w4, propproj_model64_ad_w4 |
| any of the above with w = 1 | (32, n/8) low-rank family, dyn | 0.05–0.14 | — | modesSdyn_A1_noad_w1 |

### Width scaling (the deciding question)

| carrier | n = 64 sample (K64) | n = 128 sample (A1; B1) | n = 64 model | n = 128 model | n = 256 model |
|---|---|---|---|---|---|
| full family r = 16, static | 0.011–0.018 | 0.005–0.015; 0.011–0.017 | 0.003–0.007 | 0.012–0.024 | 0.059–0.075 |
| full family r = 32, static | 0.001–0.003 | 0.001–0.004 | ≈ 0 | 0.002–0.006 | 0.023–0.031 |
| low-rank family (32, n/8), dyn | 0.006–0.013 | 0.006–0.012; 0.006–0.010 | 0.009–0.025 | 0.009–0.017 | pending |
| shared basis (32, n/4), dyn | 0.004–0.010 | 0.006–0.012 | — | — | pending |
| CP R = n/2, dyn | 0.012–0.023 | 0.012–0.025 | 0.017–0.045 | — | — |
| propagator k for 2 % (pool) | ≈ 16 | ≈ 32–40 | ≈ 20 | ≈ 40 | pending |

Readings. (i) Propagator projection needs k ∝ n (PR of the masked products ∝ n/(t+1), as for products of random
matrices): at n = 1024 it is the published shared basis, not a few vectors. (ii) On the noise-free model chain the mode
count r needed for a given ε grows ∝ n (r = 16: 0.5 % → 2 % → 7 % for n = 64 → 128 → 256). (iii) On sample truth
the n = 64 → 128 step is flat for r = 32 families and CP at R = n/2, i.e. better than the model chain predicts; the
width-256 sample atlas decides which law holds.

### Agreement / disagreement with EscAI (corpaci/ARCwhitebox, DEADENDS "response-mode expansion", killed 28 Sep)
- Agree: the born (2,1) content does not ride a few congruently transported symmetric matrices. In our terms the
  newborn/young sources have mode spectra spread over ~n modes (age-1 source: r = 16 holds only 58 % of its tensor
  energy, r = 64 96 %), and with w = 1 even the best low-rank family is 5–14 %. A young tier of dense transports is
  required (their "~6 dense n³ products/layer" floor).
- Disagree / extend: they never ran modes (A = 0). Seeded from the *aged* pool (ages ≥ 5) in the published source
  convention, symmetric mode families do carry the old content: at width 128, 16–32 modes (or a shared-basis Tucker
  (32, n/4)) reach 0.6–2.9 % with compounding. Their seeds S_a = C(s_a) correspond to our "C-feature" variant, which is
  measured dead here (ε 0.08–0.23).

## Verdict (provisional, pending width 256)

No carrier reaches ε ≤ 2–3 % at a cost clearly below the published old tier (≈ 6.7 units/layer):
- dead: dropping, windowing, regression on chain matrices, geometric tail, C-seeded modes, AD-convention carriers;
- propagator projection works but needs k ≈ 0.3 n, i.e. it *is* the published shared basis;
- dense mode families need r ≥ 16 (≥ 40 units/layer);
- the two candidates that could undercut the old tier are the **shared-basis Tucker family (r ≈ 32, q ≈ n/4)** and the
  **low-rank family (32, n/8)**, ≈ 5–16 units/layer *if r stays ≈ 32 at n = 1024*. Sample truth at 64 → 128 says it
  might; the model chain at 64 → 128 → 256 says r ∝ n, which would put them at > 100 units. Width-256 sample truth
  will be added here.

## Open issues / next steps
- Width-256 sample-truth atlas (K256, N = 2e5) → modesU / modesSdyn / propproj rows.
- Cost column is an estimate of a randomized factored implementation, not a metered build; the merge step (range
  finder + core projection of the incoming hub) dominates and should be prototyped in flopscope before trusting it.
- Model chain and sample truth disagree on width scaling; the model chain lacks the κ4 diagrams and renormalised
  coefficients (its D21 differs from the atlas by 25–50 % at depth), so sample truth is authoritative.
