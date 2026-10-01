# What the field knows: the public frontier of the ARC White-Box Estimation Challenge

*Companion to [competition-phase2.md](competition-phase2.md) (rules, environment, score). Assembled on 2026-10-01 from the forum (the Phase 1 technique census of 13 Aug and its 18 Aug update; the Phase 2 write-ups of 26 Aug, 5, 6, 10 and 17 Sep), the public repositories named below, and the community datasets on Hugging Face. Every number is quoted from its source; "raw" means final-layer MSE, "adjusted" means raw × max(0.1, C/B). Phase 1 was 256 wide × 32 deep; Phase 2 is 1024 × 16, which changes which family wins.*

---

## 0. The one-paragraph summary

In Phase 1 (narrow and deep) moment/cumulant closures hit a floor near 9e-7 raw and sampling with variance reduction (antithetic pairs, randomised lattices, Kerdock spherical 5-designs, dead/on/kink routing) owned the frontier at ≈1.5e-7 raw. In Phase 2 (wide and shallow) the balance flipped exactly as the organisers intended: plain sampling sits at ≈1.2e-6 adjusted, the kit's covariance propagation at 4e-6 raw, and **factorised third-cumulant propagation** ("K=3 cumulant chains") reaches 2.1e-8 raw at 0.25–0.5 of the budget, with one such chain published under MIT by its author (team 504aldo, rank ≈10). The top of the board (raw 1.1–1.5e-8 at 0.11–0.16 of the budget) is better on both axes, carries fourth-cumulant content and the full "old source" third-cumulant content at roughly the cost of 504aldo's young tier alone, and has published nothing. Two mechanisms recur in every honest write-up: a chained closure is an **error-compensating dynamical system** (corrections must be fitted on the chain's own rolled-forward trajectory, never injected from truth), and **per-MLP adaptivity is a transfer trap** (in-sample gains invert on fresh networks). Community moment atlases for Phase 2 exist at scale (keenanpepper: 1,000 full-split networks with joint moments at N = 1e9; 14,048 seed-regenerable networks with marginal moments; Ω-sketched third cumulants).

---

## 1. Phase 1 verdicts that carry over (the technique census, forum 18157)

| family | Phase 1 verdict | number |
|---|---|---|
| moment/cumulant closures (two-moment readout) | closed as a competitor; floor budget-independent | raw ≈ 8.6–9.3e-7 |
| k=2 closure scale bias | one scalar multiply gives ~3× | optimal factor 0.9921 ± 0.0001 |
| Edgeworth readout at the last layer, correct sign (third-cumulant term is −(g₁/6)aφ(a)) | ~10× on the readout; in-chain injection degrades | 8.85e-8 readout |
| error-compensation dynamics | per-layer errors anticorrelate; coherent corrections amplify ≈16:1 | measured |
| full moment oracle (exact m, C, γ₃, γ₄ at every layer) | the whole moment-based class is exhausted below the podium | 8.16e-9 adjusted (Phase 1) |
| plain MC + VR: float32 2.00×, shifted lattice + tent fold 1.73×, whitening 1.05×, exact first layer 1.055×, antithetic 1.007–1.05× | the honest Phase 1 frontier | ≈1.5e-7 raw |
| Kerdock / MUB-129 spherical 5-design (66,048 antipodal directions, FWHT at layer 0) | angular problem "solved"; proved within 0.023 % of optimal among fixed nonnegative rules at that support | 2.09–2.43e-7 raw |
| dead/on/kink routing (α = μ/σ < −3 dead, > 3 linear, else sample) | the shared execution chassis of the samplers | 2.18e-7 raw |
| control variates | layer-1 controls lose (0.93×); mid-network controls 2–2.4×; a perfect penultimate control would be 110× | measured |
| learned corrections | work as trajectory-fitted system components (≈8.4× over kprop); die as local patches | measured |
| score identity | above the floor a pure sampler's adjusted score is flat in N (α = 1.01 ± 0.03 on the grader) | measured |

Also measured in Phase 1 and still true: 10–16 % of the variance of the final activation sits at Hermite degrees ≥ 7 (the ReLU kink is a dense medium: 99.1 % of probability mass lies within 1e-3 of a switching surface), the answer is "born early" (layers 0–7 carry half the final mean) while the error is made late, and the grader's suite is 50 public + 50 sealed MLPs with per-MLP seeds identical across submissions.

---

## 2. Phase 2: the public ladder

| source | method | raw | C/B | adjusted | date |
|---|---|---|---|---|---|
| kit `03_covariance_propagation` | Gaussian closure, one Hermite term off-diagonal | 3.93–4.05e-6 | 0.0235 | 3.93e-7 | – |
| baghyeonsig, #328274 (forum 18200) | 0.75/0.25 blend of gain-covariance propagation and whitened antithetic MC (branch errors nearly orthogonal, r = 0.038) | – | 0.093 | 3.04e-7 | 26 Aug |
| venture_agent_ai, #329951 (18214) | exact bivariate ReLU covariance (Price + Plackett, 1-D Gauss–Legendre; −7.4 % alone) + per-layer ridge calibration on the recorded rows (feeding it back is 69× worse) | 1.37e-6 | floor | 1.32e-7 | 5 Sep |
| yujigon, #330001 (18216) | zero fitted constants: exact covariance closure + a rank-4,608 CP channel carrying the inherited joint third cumulant in O(rn); the 64 "calibration" constants of the fitted family absorb accumulated drift, not local skew | 1.07–1.13e-6 | 0.099 | 1.07e-7 | 6 Sep |
| 504aldo V29, #330093 (18218, code MIT) | factorised K=3 chain in flopscope ops + memoryless κ₄ regeneration + rank-16 D21 feedback + shared subspace for old sources + Strassen–Winograd | 2.13e-8 | 0.253 | 5.40e-9 | 10 Sep |
| top of the board (unpublished) | carries κ₄ content beyond a memoryless core and old-source content at ≈ the cost of the young tier | 1.14–1.50e-8 | 0.11–0.16 | 1.6–1.8e-9 | 1 Oct |

The reference implementation's own augmented factorised K=3 chain, un-metered, measures **1.80e-8** on the mini-8 validation set (keenanpepper's teacher bake). So the public knowledge reaches "K=3 plus a memoryless fourth" ≈ 2.1e-8; the private knowledge of the top three is worth another 1.5–2× in accuracy at 1.6–2.3× lower cost.

Pure sampling in Phase 2 sits at ≈1.2e-6 adjusted (65k forward passes buy 1/√N on a 0.075 variance); randomised Sobol collapses to the i.i.d. rate at width 1024 (18219). Analytics beat sampling by ≈200× here; sampling is "at most a floor-priced sidecar".

---

## 3. The published Phase 2 chain, in one page (504aldo, forum 18218)

Backbone: factorised K=3 cumulant propagation of Wu et al. (arXiv 2605.05179, §S.4.3) with the "K3-simple" term selection of ARC's reference code. The third cumulant is a sum of **sources**, one born at each ReLU layer, each stored as an n×n propagator leg P (identity at birth, then P ← W d(w₁) P), an n×n birth leg A = P a_b, and a few thin n×16 legs; up to 15 sources are alive at the last layer. The paper's own algorithm would cost 30 n³ L² = 3.7 × the budget at this shape.

Structural facts the author measured:
1. **The nonlinearity reads only two slices** of the third cumulant: D3ᵢ = κ₃,ᵢᵢᵢ and D21ᵢ꜀ = κ₃,ᵢᵢ꜀ (n×n). Feeding the exact D3 and D21 reproduces the exact result (oracle test F65). Fitted rule: extra final MSE ≈ 4.2e-6 × ε², ε the relative rms error of D21, so D21 must be right to ≈2 %.
2. **Exact 7 → 4 units per source-layer** (unit = one n×n×n matmul = 2³¹ FLOPs; B = 1024 units): the four (2,1) contributions collapse into two dense contractions per source.
3. **κ₄ regenerates from the covariance**: the off-diagonal core of the fourth-cumulant channel is λ_ℓ · C_off with one scalar per layer (R² 0.66–0.97 rising with depth); the diagonal must stay exact; the gain flows only through the κ₄ → κ₃ feed. This halved the score at unchanged cost (4.2e-8 → 2.2e-8). Why the core aligns with C_off is open.
4. **D21 feedback into the births is rank 16**; the final layer needs only var, D3 and the κ₄ diagonal (final-layer skew is essential: D3 := 0 there is 10× worse).
5. **Old sources live in a shrinking shared subspace**: rank ≈ 3n/8 at age 4, n/4 at age 6, 7n/32 at age 7; one notch lower is a cliff.
6. λ tracks mean(κ₄ diagonal)/mean(variance) per MLP (log–log slope ≈ 1).
7. Cost engineering: Strassen–Winograd block products as flopscope ops (sponsor confirmed permitted, no recursion limit, 11 Sep); residual time is dominated by result-buffer allocation, fixed by `out=` into pooled buffers (0.331 → 0.170 s).

Where the FLOPs go (V29, 260 units = 0.254 B): 116 units young sources (dense legs, ≤4 alive), 107 units old sources (shared basis), 25 units thin/elementwise, 7 units covariance, 6 units the whole nonlinear closure. "Our bill without the old tier is 153 units = 0.150 B. The leaders' multipliers were 0.149, 0.152, 0.164."

The dead-end table (25 rows) is the most valuable part: dropping old sources (4× per source skipped), memoryless closures on the (2,1)/(3) slices (45× worse: the next layer needs the fully off-diagonal K3, 18–23 % of D21 energy, high rank, not predictable from live features with +0.02 R²), Gaussian scale mixtures (2.3e-6), CP-rank caps (importance is flat), Tucker confinement (the P leg is born full rank), per-source κ₄ legs (≥0.95 B), learned output corrections (features correlate < 0.09 with the error), MC control-variate hybrids (needs 33–300× variance reduction, best 1.5×), marginal Edgeworth/Gram–Charlier fixes (capped at 4–5e-6), accounting tricks (none: billing is tight), float64 (no gain: all error is closure/truncation), any Hadamard product of factored legs (Khatri–Rao rank r²: "THE wall").

---

## 4. The grader's sandbox, from the harness and meter sources

- The participant process on the grader runs **flopscope-client** (numpy-free; pyzmq + msgpack) against a **flopscope-server** holding every array; arrays are `RemoteArray` handles, **immutable** (`__setitem__` and in-place operators raise), `.T` and basic slicing free; `str` paths only; 4 GiB per array; module-level and `setup()`-time handles survive across MLPs.
- Hardened image: no numpy, stdlib stripped to an allowlist (`itertools, functools, collections, time, json, contextlib` for participants; `pathlib`, `pickle`, `ctypes`, `multiprocessing`, `asyncio`, `socket`, `inspect`, `ast` poisoned), `open()` read-only under `/submission/`, no shell binaries, no network; `fft.*` unavailable on the client; `vectorize`/callback ops unavailable; `stats.*` always float64.
- Billing identity `cost = flop_cost × dtype_rate × complex_factor × weight`, checked **before** the op executes; one worker serves the whole suite; `setup()` runs once per worker (5–15 times per submission) with a 5 s cap that includes interpreter start-up and imports; the watchdog finalises when no MLP finished for 8 min or 40 min total.
- Aliasing is by Python object identity: `fnp.einsum("ji,jk->ik", X, X)` and `fnp.inner(X, X)` get the symmetric-output discount (1,074,265,600 for 1024²), `X.T @ X` does not (2,146,435,072).
- Costs at n = 1024 (float32): matmul 2.146e9; matvec 2.1e6; `eigh` 9.66e9 (9n³); `cholesky` 3.58e8 (n³/3); `svd` thin 2.79e10, top-64 2.68e8; `inv` 2.15e9 (symmetric 1.43e9); `exp/log/arccos/x**2` 16 per element; `where` 4 per element; `norm.cdf` 48 per element (float64 → 96), `norm.pdf` 27 (→ 54); `standard_normal` 16 per element at float32; `as_symmetric` 7n² − 1; `fill_diagonal` n; `reshape`/`copy`/`astype` n²; `zeros`, views, `diag` extract free. No `erf`.
- The grader's smoke test probes a **non-suite shape** (depth 32 was seen): every depth- or width-indexed table must be guarded with a fallback path, or the whole submission fails.

---

## 5. Community datasets (Hugging Face, MIT, keenanpepper)

| repository | contents | size |
|---|---|---|
| `arc-whestbench-p2-higher-moments-2026` | mini split (100 MLPs): per-layer pre/post marginal raw moments to order 4, dense `pre_M11, pre_M21, pre_M31, pre_M22, M11` at N = 1e8; plus a closure-training cache (κ₂₁, κ₃₁, κ₂₂ pair features, Mehler-8 residuals) | 100 × 322 MB + 100 × 151 MB |
| `arc-whestbench-p2-full1000-N1e9` | all 1,000 full-split MLPs: the same joint moments at **N = 1e9** (relative MC noise ≈ 3e-5) | 1,000 × 337 MB |
| `arc-whestbench-p2-d8b-corpus-14k` | 14,048 fresh He-init MLPs (weights regenerate from `torch.Generator("cpu").manual_seed(770000 + idx)`): per-layer post-ReLU marginal moments m₁..m₄ and final means at N = 1e8 (noise floor ≈ 9e-10); for the first 2,048 also a **teacher bake** of ARC's augmented factorised K=3 chain state (κ₄ diagonal, rank-8 factors of the κ₂₂ slice, pre-activation mean/var, teacher final mean; 1.8023e-8 on mini-8) | shards |
| `whest-p2-bakev2-{d8b-sketch-g00, bench-supp, mini-supp}` | Ω-sketched third cumulants (every net) and dense `post_K21, post_K22, gate_GG, gate_GX` blocks (validation nets), marginal cumulants to order 6, gate probabilities, noise floors from independent re-bakes | per layer files |

These are the training and validation surfaces for anything learned; the full-split N = 1e9 joint moments are the exact targets for closure diagnostics (true D21 and D3 at every layer, hence oracle tests of the kind 504aldo ran on 8 networks can be run on 1,000).

---

## 6. What this implies

1. **The public frontier is reproducible today**: 504aldo's V25 (no Strassen, 0.367 B) and V29 (0.253 B) are MIT, self-contained single files, and were already re-submitted by 13 other accounts. Reproducing them locally and on the grid is the first step, not an achievement.
2. **The gap to the money line is a representation problem, not a tuning problem.** The leaders carry old-source and κ₄ content at ≈0.15 B; every compression of 504aldo's leg representation hits the Khatri–Rao wall. The published hints point at (a) a different projection of the fourth order (SSC-style response closures: regenerate κ₄ from a few live modes with a fitted coupling), (b) learned or sketched memories of the all-distinct third cumulant (keenanpepper's Ω-sketch targets exist for exactly this), (c) a state that is "angular rather than moment-based" (trim_qewas's oracle: the angular factor alone reaches 0.28 % closure error).
3. **Fitting is allowed and necessary, with two rules**: fit on the chain's own rolled-forward trajectories at scale (thousands of networks, frozen before evaluation), and never adapt per MLP on in-sample statistics. Shipped tables are data; a reviewer may ask what a table absorbs (yujigon's question), so the derivation should be known.
4. **Robustness is where ranks are lost**: one zeroed MLP costs 9e-3 adjusted against a 1.6e-9 target; the smoke test runs a foreign shape; the residual clock is allocation-bound; `setup()` must survive 5 s cold including imports.

---

## Sources

- Forum: 18157 (technique census, 13/18 Aug 2026), 18193 (Phase 1 solution and Phase 2 ideas, 20 Aug), 18200 (26 Aug), 18214 (5 Sep), 18216 (6 Sep), 18218 (10 Sep, 504aldo), 18219 (17 Sep), 18197 (Phase 2 announcement), 18151/18097/18154/18182/18177/18175/18173/18171/18176 (Phase 1 write-ups cited by the census).
- Code: github.com/504aldo/whest-p2-cumulant-k3 (MIT); github.com/jamesrahenry/arc-whitebox-replication; github.com/01-1/arc-wbe; github.com/Oishi1029/arc-whestbench-2026; github.com/agentaiventure-dot/arc-whest-estimator; github.com/AIcrowd/{whestbench, flopscope, whest-starterkit}.
- Datasets: huggingface.co/datasets/keenanpepper/{arc-whestbench-p2-higher-moments-2026, arc-whestbench-p2-full1000-N1e9, arc-whestbench-p2-d8b-corpus-14k, whest-p2-bakev2-*}; huggingface.co/datasets/aicrowd/arc-whestbench-public-2026@v2-phase2.
- Paper: W. Wu, V. Lecomte, M. Winer, G. Robinson, J. Hilton, P. Christiano, arXiv:2605.05179 (Table 1: K=3 factorised 30n³L² + 39n³L; K=2 basic 7n³L).
