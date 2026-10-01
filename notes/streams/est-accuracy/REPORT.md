# est-accuracy: lowering the raw MSE of the public chain through its fourth-cumulant handling and births

*Status: IN PROGRESS (started 1 Oct 2026, ~17:00 UTC). Numbers below are measured; open items at the end.*

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
