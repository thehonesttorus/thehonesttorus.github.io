# F1: driving V29's cost down (age axis, joins, thin legs)

Child session of the frontier coordinator, 2 Oct 2026. Base: `../estimator_v29r3.py` (504aldo V29, MIT), unedited.
Reward: dev adjusted = mean raw × max(0.1, C/B) on w1024_d16. Pairs are always on the same MLPs.
C/B is the sweep convention (max over the non-first MLPs of a process). The steady third-MLP value is about 0.003 lower.

## Result

**Best: `est_F1_final.py`.** It is a standalone single file with the defaults baked in. The same arithmetic is `est_incj.py` with the env below.

On all 6 MLPs it reaches **dev adjusted 3.86e-9**, against 4.46e-9 for base: **−13.6 %**. Raw 1.719e-8 (−2.2 %) at C/B 0.2245 (steady 0.2213), against 0.254 for base (−11.6 %).

| config (6 MLPs) | raw MLP 0..5 (×1e-8) | mean raw | C/B | dev adjusted | Δ vs base |
|---|---|---|---|---|---|
| base V29 | 1.942 0.384 2.074 1.195 2.531 2.420 | 1.758e-8 | 0.2540 | 4.464e-9 | 0 |
| knob stack, no incj (thin 8/8, FB 1.5, LAM 0.80, L6) | 1.924 0.317 2.079 1.096 2.584 2.465 | 1.744e-8 | 0.2307 | 4.024e-9 | −9.9 % |
| coordinator syn_rfb8 (incj JS64 JPASS2, R_FB 8, FB, LAM, L6) | 1.845 0.151 1.971 1.062 2.599 2.624 | 1.709e-8 | 0.2328 | 3.977e-9 | −10.9 % |
| incj JS32 JPASS1, R_FB 8, R_RES 16, FB, LAM, L6 | 1.858 0.121 1.955 1.081 2.697 2.659 | 1.728e-8 | 0.2284 | 3.947e-9 | −11.6 % |
| **incj JS32 JPASS1, thin 8/8, FB 1.5, LAM 0.80, L6 (= est_F1_final)** | 1.877 0.124 1.960 1.038 2.720 2.594 | **1.719e-8** | **0.2245** | **3.859e-9** | **−13.6 %** |

The final stack combines:
- **F1 code: the incremental exact-core join** (below), JS = 32, JPASS = 1.
- R_FB = R_RES = 8.
- From the coordinator: D21-feedback hyperedge ×1.5, LAM scale 0.80, and Strassen L6 with 16-leaves.

Base MLPs 3–5 come from the coordinator's `sw_base_345.json`.

## The F1 code change: incremental exact-core join (`est_incj.py`)

**What V29 does.** It rebuilds the shared basis at every join with a one-pass rank-r randomized range finder, with no oversampling, on the weighted leg Gram (old + joiner).

**What F1 does instead:**
- Keep span(Qp) exactly: Qp = w1 ⊙ Qc, orthonormalised once with QR.
- Add a JS-column sketch of the joiner's residual directions, orthogonal to Qp.
- Truncate the (r + JS)-dimensional space to the top-r subspace of the **exact weighted core** (old core lifted by R, plus the joiner's exact Gram in that space). This takes one QR pass from the old-basis start, or 2 passes, or eigh.
- The joiner's factors come from V^T·[Qoᵀ A; Qsᵀ A] at r(r+s)n cost. The old-factor rotation is V1ᵀ R.

**What it buys:**
- A more accurate basis: MLPs 0–2 raw −7 % at equal knobs.
- A slightly cheaper join: C/B −0.002 at JS 64, −0.017 at JS 32 with thin 8.

**What it does not buy:**
- On held-out MLPs 3–5 alone it is raw-neutral (+0.8 %).
- The per-MLP spread (±3–5 %) is about the size of the effect. The honest 6-MLP gain of incj at fixed other knobs is about −2 to −3 % raw plus a small cost saving.

Kill-switch: `F1_INCJ=0` restores V29's join exactly.

## Screens on MLPs 0–2 (base 1.467e-8 at 0.2540, adjusted 3.725e-9)

| config | raws (×1e-8) | mean raw | C/B | adjusted | Δ |
|---|---|---|---|---|---|
| a3_o5_192 (AGE_OLD 3, AGE_OLD2 5, r2 192) | 2.315 1.231 2.793 | 2.113e-8 | 0.2339 | 4.94e-9 | +33 % |
| a3_r384_o6_160 | 2.272 1.185 3.262 | 2.240e-8 | 0.2338 | 5.24e-9 | +41 % |
| a2_r512_o4_192 | 4.367 2.182 4.022 | 3.523e-8 | 0.2569 | 9.05e-9 | +143 % |
| thin8 (R_FB = R_RES = 8) | 1.959 0.439 2.110 | 1.503e-8 | 0.2425 | 3.65e-9 | −2.1 % |
| incj (defaults) | 1.877 0.225 1.984 | 1.362e-8 | 0.2521 | 3.43e-9 | −7.8 % |
| incj + T2EIG (exact tier-2 eigh) | 1.842 0.253 2.019 | 1.371e-8 | 0.2529 | 3.47e-9 | −6.9 % |
| incj a3_o5_192 | 2.821 1.244 2.756 | 2.274e-8 | 0.2319 | 5.27e-9 | +42 % |
| incj r256 / r2 128 | 2.669 1.114 3.083 | 2.289e-8 | 0.2093 | 4.79e-9 | +29 % |
| incj R_FB 8 | 1.896 0.269 1.994 | 1.386e-8 | 0.2446 | 3.39e-9 | −9.0 % |
| incj thin8 | 1.899 0.286 2.011 | 1.398e-8 | 0.2407 | 3.37e-9 | −9.6 % |
| incj thin8 JS32 | 1.890 0.272 2.034 | 1.399e-8 | 0.2375 | 3.32e-9 | −10.8 % |
| incj thin8 r320 / r2 160 | 1.939 0.697 2.410 | 1.682e-8 | 0.2162 | 3.64e-9 | −2.4 % |
| incj thin8 r288 / r2 128 | 2.155 1.029 2.796 | 1.993e-8 | 0.2049 | 4.08e-9 | +9.6 % |
| incj thin8 r384 / r2 160 | 1.849 0.397 2.154 | 1.467e-8 | 0.2324 | 3.41e-9 | −8.5 % |
| incj thin8 AGE_OLD2 6, r2 192 | 1.871 0.503 2.281 | 1.552e-8 | 0.2322 | 3.60e-9 | −3.3 % |
| incj thin8 AGE_OLD 5, AGE_OLD2 8 | 1.999 0.299 1.969 | 1.422e-8 | 0.2494 | 3.55e-9 | −4.8 % |
| stack JS32 (L6, FB 1.5, LAM 0.80, thin8) | 1.857 0.149 1.994 | 1.333e-8 | 0.2257 | 3.01e-9 | −19.2 % |
| stack JS16 | 1.881 0.220 2.032 | 1.378e-8 | 0.2242 | 3.09e-9 | −17.1 % |
| stack JS32 JPASS1 | 1.877 0.124 1.960 | 1.321e-8 | 0.2245 | 2.96e-9 | −20.4 % |
| JE = 2 (join every 2 layers, `est_incj2.py`), MLPs 0–1 only | 1.954 0.324 | n/a | 0.2534 | worse | dropped |

Held-out MLPs 3–5 (base 1.195 / 2.531 / 2.420 → 2.049e-8):

| config | raws (×1e-8) | mean raw | Δ raw vs base |
|---|---|---|---|
| incj | 1.120 2.569 2.507 | 2.065e-8 | +0.8 % |
| incj thin8 JS32 | 1.196 2.729 2.519 | 2.148e-8 | +4.9 % |
| incj FB 1.5, LAM 0.80 | 0.982 2.491 2.593 | 2.022e-8 | −1.3 % |
| incj thin8 FB, LAM | 1.019 2.628 2.561 | 2.069e-8 | +1.0 % |
| final stack (L6) | 1.038 2.720 2.594 | 2.117e-8 | +3.3 % |

The 0–2 screens over-state every raw gain: MLP 1's noise-subtracted raw is about 0.2–0.4e-8 and swings by ±0.3e-8. Read only 6-MLP numbers as final.

## What limits further gains (why 0.1 B is far)

**The young tier.** Ages 1–4 are dense, about 2.2 u per source-layer, about 113 u ≈ 0.11 B on their own.
- Every move that takes a young age out costs far more raw than it saves.
- AGE_OLD 3 gives +19 to +45 % raw for only −8 % cost. AGE_OLD 2 is catastrophic, with or without incj.
- This matches G10: fresh layers cannot be compressed in neuron space.

**The old-tier ranks sit at the cliff.**
- R_OLD 320 or r2 160 already cost +15 % raw on 0–2. 256 / 128 cost +60 %.
- Better truncation (incj) moved the cliff only a little.
- The old tier's cost is linear in the number of old sources, because dense legs are formed (Qc·F, n²r) and contracted per source. Lower ranks alone cannot fix that. It needs a cohort state (E2).

**Joins.** These are about 50 u in V29.
- Incj made them a little cheaper.
- Joining every 2 layers (`est_incj2.py`, JE=2) costs more than it saves: the delayed joiners stay dense.

**Fixed items** (C_pre, births, thin legs ≈ 30 u after thin 8) are already small.

**Summary.** The reachable floor with the V29 architecture is ≈ 0.21–0.22 B. Getting to 0.1 B (lode_dockx) needs a different old-tier representation, plus a young tier of at most 2 dense ages, with the raw loss compensated elsewhere.

## Files

- `est_incj.py`: V29 + incremental join, with env knobs.
  - F1_INCJ, F1_JS, F1_JPASS, F1_JEIG, F1_T2EIG.
  - F1_RFB, F1_RRES.
  - H_FB_SCALE.
- `est_incj2.py`: multi-joiner variant (F1_JE). Not recommended.
- `est_thin.py`: thin-rank wrapper of the unedited base.
- **`est_F1_final.py`: the final pick, standalone.**
- `sweepF1.py`, `cfg*.json`, `summ.py`.
- `results/sw_*.json|log` and `results/sweep*.log`.
- Residuals measured in-process here are inflated by 2 concurrent jobs (0.45–0.6 s), so the grader-faithful timing check still has to be run (coordinator).
