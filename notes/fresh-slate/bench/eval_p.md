# Stage P: running a flopscope estimator at 1024 × 16 with the grader caps

Environment (every container): `python3 -m venv /root/whest && /root/whest/bin/pip install "whestbench==0.16.1" "flopscope==0.12.1" numpy scipy pyarrow psutil`.

## Truth

- Set `w1024_d16` (when present in `sets/`): 6 MLPs, input seeds 7301001–7301006 (the submission stream's dev set),
  N = 2e6, truth noise avg_variance/N ≈ 3.6e-8 per MLP. That noise is larger than the best raw errors (~1e-8), so
  always report raw = MSE − avg_variance/N (s.e. of the 6-MLP mean of raw ≈ 1e-9) and compare designs PAIRED on the
  same MLPs (the truth noise cancels to first order in a paired difference).
- `whest run` needs the parquet dataset. Re-bake it in your container (identical bytes given the seeds):
  `echo '[7301001,7301002,7301003,7301004,7301005,7301006]' > s.json; OMP_NUM_THREADS=4 whest dataset bake --n-mlps 6 --n-samples 2000000 --width 1024 --depth 16 --mlp-seeds s.json --output dev6`
  (≈ 8 min per MLP on 4 idle cores; run detached). Any other set can be re-baked the same way from its `sets/<name>.json`
  (`seeds`, `n_samples`, `width`, `depth`).

## Quick in-process check (no parquet needed)

`/root/whest/bin/python run_p_inproc.py --estimator path/estimator.py --set w1024_d16 [--mlps 0,1]`
— weights are regenerated from the seeds (bit-identical to the bake), the estimator runs under
`BudgetContext(2^41)`; prints per-MLP raw, C/B, residual, wall and the adjusted score raw × max(0.1, C/B).
Residual times are only meaningful on an idle box; in-process memory is not the grader's (arrays live in the
flopscope server there).

## Grader-faithful run

`/root/whest/bin/python run_p.py --estimator path/estimator.py --dataset dev6 --runner subprocess --max-threads 2 --tag NAME --out NAME.json`
(a copy of the submission stream's `measure_run.py`): `whest run --runner subprocess --format json --profile` with the
graded caps (2^41 FLOPs, 120 s wall, 0.4 s residual, 5 s setup), peak RSS of the worker tracked from outside (8 GB is the
grader's limit; the in-process runner keeps every array in the worker, so local RSS overstates the grader's), truth noise
subtracted. For the grader's client/server transport (residual measured as the client's wall − dispatch), see
`notes/streams/submission/scripts/grader_emul.py`.

## Robustness checks before any submission

- Smoke shape 256 × 32 (set `w256_d32`): the grader smoke-tests a 256-wide, 32-deep MLP; any layer-indexed table that
  assumes depth 16 fails the whole submission.
- Other shapes (512 × 16, 1024 × 4, 256 × 8, 1024 × 32) and adversarial weights (×10, |W|): no exception, finite output,
  within budget. A failure zeroes the MLP at multiplier 1.0.
- Static: never assign `x.shape = ...` (fails on the grader); numpy is not installed at evaluation (flopscope + stdlib
  only); `setup()` ≤ 5 s and idempotent; keep residual ≤ 0.2 s locally (2× margin).
