# est-accuracy: lowering the raw MSE of the public chain through its fourth-cumulant handling and births

*Status: STOPPED at 17:45 UTC on 1 Oct 2026 on the user's instruction (the competition system is to be designed from a
fresh slate; adapting the public chain is out). No further public-chain work will be pushed from this stream.*

**Stop note.** Measured so far, on ONE width-1024 MLP only (seed 7301001, truth N = 1e5, so the paired-difference truth
noise is ~0.4–1.5e-9 and MLP-to-MLP variation is unmeasured; the 6-MLP dev bake did not finish in time): relative to the
patched V25 (final MSE 7.4995e-7 incl. truth noise), the exact transported κ4 (4)-slice diagonal in the regenerated core
(`EA_XDIAG=1 EA_XD4=1`: (W∘W)²K4v replacing the projected (W∘W)·cA·K4v, +0.00004 B) gave −5.8e-9; adding the exact
(2,2) part as well (+0.015 B) gave −5.9e-9; restoring the B1 D21-feedback leg to its leg-partition coefficient
(`EA_FBX=2`, the port drops half of B1 and B2 relative to the reference chain) −2.1e-9 (×3: −3.7e-9); the Edgeworth
B3 dressing of the newborn hub vertex (w2 += w5 D3/6) −1.0e-9; these three were additive (combo −8.8e-9, i.e. of the
order of 30–40 % of the chain's ~2e-8 error, if it held on more MLPs). Neutral or harmful: κ4 vertex dressing (−0.2e-9),
dressed birth legs (+0.3e-9), B5 with the chain's rank-2 K22 (+1.1e-9), restoring B2 (`EA_FBY=2`, +5.5e-9), EscAI's
path/star κ4-diagonal births added to the core (+2.3e-9), physical λ (0.3× table, no adaptation: +1.1e-8), the one-step
exact (2,1,1) feed of the κ4 diagonal (+3.3e-8, +0.05 B), and exact-transport (2,2)/(3,1) use-side slices (the partial
(2,2) slice diverges: 5.6e-6). Raw results: `/tmp` only (not committed); code: `scripts/est_lab.py` (env-gated levers),
`scripts/lab_run.py`, `scripts/compare.py`, `scripts/linfit.py`. Unverified on the dev set; treat as leads, not results.

## Question

Can the raw final-layer MSE of the public 504aldo chain (patched V25 bundle, C/B 0.3667) be lowered at width 1024
by better fourth-cumulant handling and better births, within budget? Levers: (a) exact κ4 diagonal / slices instead of
regenerated; (b) a better (2,1,1) slice than the regeneration u_i C_jk; (c) the physical κ4 λ; (d) leg-partition closure
births missing from the chain; (e) EscAI's deployable pieces.

## Method

- Code: `scripts/est_lab.py` = the patched V25 bundle (`../submission/bundles/v25/estimator.py`) plus env-gated levers
  (all default off: with no env set it is bit-identical V25 arithmetic). `scripts/lab_run.py` runs one estimator
  instance across the MLPs in-process, keeps every layer's prediction; `scripts/compare.py` makes the paired tables.
- Truth: `whest dataset bake --width 1024 --depth 16 --n-mlps 6 --n-samples 2e6`, seeds 7301001..7301006 (the
  submission stream's dev set, so V25's raw is cross-checked against its 1.756e-8). Raw = MSE − avg_variance/N.
  Variants are compared PAIRED on the same MLPs (difference of MSEs; truth noise largely cancels).

## Results

(pending)

## Open items

(pending)

## Working notes (to be folded into Results)

- Coefficient audit of the port's births against the reference/lean chain (lean_k3_aug `X1, Y1, X2, Y2`) and the
  leg-partition closure: the port's `Y1 = A d(w2) + Yt` drops the reference's `c_w2 r_w2 WK12` term and its shipped
  config shadows the B2 block (`sb2`). Both dropped terms duplicate a kept one, so the port carries the closure
  diagrams B1 (D21 hyperedge + edge, `Xt = 1.5 d(w2) D21`) and B2 (`Yt = 0.5 d(w1) D21ᵀ d(w3)`) at HALF their
  leg-partition coefficient (0.25 Σ₆perm vs 0.5 Σ₆perm). Levers `EA_FBX`, `EA_FBY` scale them (2 = full).
- The oracle's `CLOSURE_COEF` B3 = 1 looks like a double count: B3's two leaves are symmetric (like B6), so the
  per-hub Edgeworth weight (1/6) w5 D3 gives 0.5 × sym3; the width-1024 fit (0.17–0.46) is consistent with 0.5.
- The port has no B3 (D3 hyperedge + 2 edges) and no κ4-diagonal-hub term in its births; both are a column
  scaling of the newborn hub vertex weight (w2 → w2 + w5 D3/6 + w6 κ4/24, the Edgeworth-dressed E[relu'']):
  levers `EA_B3`, `EA_K4V`, free.
