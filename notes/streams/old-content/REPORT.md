# Stream old-content: a cheap carrier for the transported old-source content of κ3

Status: done (2026-10-01). Widths 64, 128, 256 on sample truth (full third-moment atlases) and on a noise-free model
chain; two MLPs and an independent-sample replicate at width 128; a Monte Carlo noise check.

## Question

The published K = 3 chain carries κ3 as per-layer sources, each transported with dense n×n legs. Its old-source tier
(age ≥ 5, shared basis 384 → nested 224) costs 107 of 260 units (≈ 6.7 units/layer on average, ≈ 10 on the layers where
it is active). Leaders run 0.11–0.16 B in total (≈ 7–10 units/layer for everything). Is there a carrier of the old
content's effect on D21(l+1)_ab = κ3(z_a, z_a, z_b) with ε ≤ 2–3 % at ≤ 10 units per layer (n = 1024), ideally near zero?

## Verdict

**No.** None of the tested carriers, (a)–(e) plus the propagator projection, reaches ε ≤ 2–3 % at n = 1024 for ≤ 10
units/layer. Several reach it at width 64–128. At every width the size each carrier needs grows at least in
proportion to n, and the step from 128 to 256 breaks the flat 64 → 128 trend that looked promising:

| carrier at fixed relative size, dynamic (compounding), sample truth | n = 64 | n = 128 | n = 256 |
|---|---|---|---|
| shared-basis mode family (Tucker r × q × q), r = 32, q = n/4 | 0.004–0.010 | 0.006–0.012 (B1: 0.005–0.012) | **0.025–0.044** |
| low-rank mode family, 32 modes of rank n/8 | 0.006–0.013 | 0.006–0.012 (B1: 0.006–0.010) | **0.024–0.047** |
| low-rank mode family, 16 modes of rank n/8 | 0.007–0.013 | 0.017–0.031 | **0.045–0.129** |
| full symmetric mode family, r = 16, static (oracle bound) | 0.011–0.018 | 0.005–0.015 | **0.043–0.058** |
| full symmetric mode family, r = 32, static | 0.001–0.003 | 0.001–0.004 | **0.012–0.023** |
| propagator projection, k needed for ≈ 2 % on the pool | ≈ 16 | ≈ 32–40 | ≈ 64–80 |
| participation ratio of the masked propagator M_{s→t}, age 5 | 6.0 | 10.5 | 23.0 |

Extrapolated to n = 1024, the mode count needed for 2 % grows by ≈ 2–2.5× per width doubling (r ≈ 12 → 30 from 128 to
256 for the full family). That gives r ≈ 150–250 modes, i.e. ≥ 60 units/layer for any mode form. Propagator
projection needs k ≈ 0.25–0.3 n ≈ 300, which is exactly the published shared basis. The old content of the third
cumulant is not a few-mode object at competition width: its effective rank grows with n, like the participation
ratio of products of random masked matrices (PR ≈ n/(t+1)).

So the leaders' near-free old content (EscAI: "a construction that the entire measured candidate space ... does not
contain") is not any compression of the transported κ3 old content that we could find. It is either (i) not carried
at all, with its D21 effect absorbed in renormalised closure coefficients or a different closure (see the coef-ensemble
stream), or (ii) a representation outside the source/transport picture altogether.

## Method

- `tracker.py`: exact per-source decomposition of κ3(z_l) on a full third-moment atlas (dense n³, n ≤ 256).
  It supports two exact conventions:
  - **AD** (default): the old content is O_l = AD(Φ³ κ3(z_l)); everything else is a birth.
  - **no-AD** (`--noad`, the published chain's convention): sources are transported with the full Φ³, so
    X_{s,t} = M_{s→t}^{⊗3} B_s exactly, and the slice corrections are absorbed into each birth.

  Validation: Σ_s X_{s,l} reproduces the atlas κ3(z_l) at every layer to 1e-6–1e-4. The 'valid' column of
  `results/tracker_*.txt` reports this; the residue is the float32 accumulation of the atlas. `--model` runs a
  noise-free closure chain (births = first-order closure from the chain's own κ3(z), atlas slices) from a pair-only atlas.
- `carriers.py` and `propproj.py` implement the carriers. The error is ε_l = ‖D21(Ôld_l) − D21(Old_l)‖ / ‖D21(z_l)‖,
  with the whole D21 as denominator, i.e. the ε of the chain's error law (extra MSE ≈ 4.2e-6 ε²).
  - Old = sources of age > w, with w = 4 (the published young tier) unless stated.
  - **dyn**: the carrier feeds its own approximation forward, so errors compound.
  - **static**: re-fitted to the true Old_l at each layer (an oracle bound).
  - Re-truncations inside dynamic carriers are oracle (dense SVD / HOSVD / ALS). Costs are estimates for a randomized
    factored implementation, not metered.
- `run_all.sh`: the width-128 batch. Atlases (not committed, regenerable with `notes/experiments/moment_atlas_np.py`):
  - A1, A2: seed 770100, width 128, N = 6e5 each, independent sample seeds.
  - B1: seed 770101, width 128, N = 6e5.
  - K64: seed 770100, width 64, N = 6e5, `--k3`.
  - K256: seed 770100, width 256, N = 2e5, `--k3` with float32 accumulators.
  - P64, P64b, P256: pair-only atlases for the model chain.
  - t1: seed 770100, width 128, N = 8192, used for the noise check.

## Results

### Where the D21 content lives (tracker, width 128, A1)
- AD convention: the age-1 source is 0.75–0.93 of ‖D21‖. Ages ≥ 2 together give 0.37–0.63, and ages ≥ 5 give 0.06–0.46.
- No-AD convention: ages ≥ 5 carry 0.10–0.79. Young and old partly cancel through the slice corrections.
- Old sources' D21 contributions are nearly orthogonal to each other and to the young ones (see the regressions below).
- The leading tensor mode of the old pool is the mean direction: cos(v₁, μ) ≈ 0.93, holding 60–80 % of the tensor
  energy at depth. It carries almost no D21. Tensor energy is the wrong metric, so only D21 error is reported.
- Noise does not drive any carrier number:
  - A1 vs A2 (independent samples, same MLP) agree to ±0.005 in every ε.
  - B1 (second MLP) falls in the same ranges.
  - At width 128, N = 8192 and N = 6e5 give the same carrier ε (`modesSdyn_t1N8192_noad_w4.txt` vs `modesSdyn_A1_noad_w4.txt`).

### Carrier table (width 128, sample truth A1, no-AD, w = 4; ε range over layers 5–15; per-layer numbers in `results/`)

The cost column gives estimated units per layer at n = 1024: first at the size that works at width 128, then at the
size n = 1024 needs according to the width scaling above.

| carrier | variant | ε at width 128 | est. cost at n = 1024 | file |
|---|---|---|---|---|
| none (drop old content) | — | 0.10–0.79 | 0 | carriers_A1_noad_w4 |
| (c) regression on O(n²) chain objects | C, C∘C, var/μ outer products, slice transports and sandwiches of D21(z_{l−1}), young D21 (in-sample, per layer) | 0.075–0.45 (AD: 0.06–0.43) | ~6–10 | carriers_A1_{noad,ad}_w4 |
| (c) old pool on the young sources' D21 | in-sample | 0.09–0.63 | 0 | carriers_A1_noad_w4 |
| (d) geometric tail | one more exact source, ages ≥ 6 = γ_l × it (+ young) | 0.055–0.52 | ~2–4 | carriers_A1_noad_w4 |
| (a) full symmetric mode family, r dense S_r transported as WᵀΦSΦW | r = 8 / 16 / 32, dyn | 0.028–0.074 / 0.014–0.027 / 0.003–0.005 | ≈ 2.5 r: 40 at r = 16; ≥ 400 at the r that n = 1024 needs | carriers_A1_noad_w4 |
| (e) mode family with directions fixed to the top-r eigenvectors of C_l | r = 16 / 32, static | 0.020–0.036 / 0.007–0.015 | ≈ 2.5 r | carriers_A1_noad_w4 |
| (e) S_r fitted on C-features (C, C∘C, μμᵀ, …; ≈ EscAI's C-seeded modes) | r = 16, static | 0.08–0.23 (worse than dropping at shallow layers) | ~0 | modesS_A1_noad_w4 |
| (e) low-rank mode family, each S_p of rank q | (16, n/8) / (16, n/4) / (32, n/8), dyn | 0.017–0.031 / 0.014–0.027 / 0.006–0.012 | ≈ 4 r q / n: 16 at (32, 128); ≥ 60 scaled | modesSdyn_A1_noad_w4 |
| (e) shared-basis family (Tucker r × q × q, one U for all modes; q = n/4 ↔ 256) | (16, n/4) / (32, n/4) / (32, n/2), dyn | 0.016–0.029 / 0.006–0.012 / 0.002–0.005 | transport (r+q)/n, readout ≈ r q²/n², merge ≈ 3 r q²/n² + 2: ≈ 6–10 at (32, 256); ≈ 50+ at the r that n = 1024 needs | modesU_A1_noad_w4 |
| (b/d) CP, i.e. one merged hub of R columns, legs transported, ALS merge | R = n/2 / n, dyn | 0.012–0.025 / 0.007–0.014 | 3 R/n transport+readout + ≥ 4.5 per ALS sweep: ≈ 15–20 at R = n/2 | cpdyn_A1_noad_w4 |
| propagator projection (top-k left singular vectors of M_{s→t}, per source) | k = 16 / 32 / 64, static, pool | 0.05–0.096 / 0.014–0.027 / ≤ 0.004 | k ≈ 300: the published shared basis (≈ 6.7) | propproj_A1_noad_w4 |
| any carrier in the AD convention | — | 2–5× worse at equal size (re-masking scatters content out of the propagator range) | — | carriers_A1_ad_w4, propproj_model64_ad_w4 |
| any carrier with w = 1 | low-rank (32, n/8), dyn | 0.05–0.14 | — | modesSdyn_A1_noad_w1 |

CP at width 64 (sample truth): R = n/2 gives 0.012–0.023 and R = n gives 0.006–0.012, flat from 64 to 128. Width 256 was
not run: ALS on 256³ is too slow here. The mode families suggest it would degrade the same way.

### Model chain vs sample truth
The noise-free model chain shows the same growth, but earlier: full family r = 16 static gives 0.003–0.007 / 0.012–0.024 /
0.059–0.075 at n = 64 / 128 / 256. Shared basis (32, n/4) dyn gives 0.010–0.019 at n = 64 and 0.011–0.020 at n = 128.
Its D21 differs from the atlas by 25–50 % at depth: no κ4 diagrams, no renormalisation. Sample truth is authoritative;
both give the same verdict.

### Agreement with EscAI (corpaci/ARCwhitebox, DEADENDS "response-mode expansion", killed 28 Sep)
- **Agree.** Born and young content does not ride a few congruently transported symmetric matrices:
  - The age-1 source needs r = 64 for 96 % of its tensor energy at width 128.
  - With w = 1 even the best low-rank family is 5–14 %.
  - A dense young tier is required, consistent with their ~6 dense n³ products per layer.
- **Agree.** C-seeded modes, our "C-feature" fit, are dead: ε 0.08–0.23.
- **Extend.** They never ran modes (A = 0). We did, seeded from the aged pool in the published convention:
  - At width 64–128, mode families carry the old content to 0.4–1.2 % with compounding.
  - The required mode count grows with n; at width 256 the same relative size gives 2.4–4.7 %.
  - Their cost estimate (0.13–0.15 B) applies to a few modes. That is the wrong size at n = 1024.

### Answer to the coordinating session's transfer-operator / expander hypothesis
- The mechanism is right in the published convention: old sources concentrate on the top singular vectors of their
  propagator. Per source, k = 32 keeps 95–99 % at ages ≥ 6 (width 128).
- The PR of M_{s→t} at age 5 is 6.0 / 10.5 / 23.0 at n = 64 / 128 / 256, i.e. ∝ n.
- The pool needs k ≈ 0.25–0.3 n per band, so this is the published shared basis, not a near-free carrier.

## Files
- `tracker.py`, `carriers.py`, `propproj.py`, `run_all.sh`
- `results/tracker_*` per-source D21 shares and validation;
  `results/carriers_*`, `modesS*_*`, `modesU_*`, `cp*_*`, `propproj_*` carrier ε per layer.
  Naming: A1/A2/B1 = width 128 sample truth, K64/K256 = width 64/256 sample truth, `model*` and `*model` = model chain.

## Open issues / next steps
- Costs are estimates for a factored implementation (not metered). The verdict does not depend on them: the required
  sizes grow with n.
- Widths 512–1024 are out of reach for dense n³ tensors here. A streaming version of the D21-ε measurement for one
  carrier at width 1024 (oracle1024-style) would confirm the extrapolation.
- The useful positive result is negative information for Line C: the leaders' cheap old content is not a compression
  of the transported old κ3. Test instead whether *dropping* the old tier and refitting the closure coefficients (or
  regenerating the D21 shortfall from O(n²) state with fitted, ensemble-level coefficients) recovers the published
  raw error. Regression (c) here was in-sample per layer on ≤ 12 features and failed. A learned or table-driven
  correction on the full ensemble has not been tested.
