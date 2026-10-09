# flopscope 0.12.1 billing rules (as installed in the grader venv `/opt/wb`)

Source of truth: `flopscope/data/default_weights.json` (weights, dtype rates), `flopscope/_flops.py`,
`_accumulation/`, and the published `docs/reference/cost-model.md` (copy in the scratchpad `sources/github/`).
Measured bills below come from `scratchpad/fnp/bill_probe.py` at n = 1024 (one n^3 = 1.07e9).

## The one equation
    charged = int(flop_cost * dtype_rate * complex_factor * weight)
- `flop_cost`: textbook operation count (FMA = 2), **not** what BLAS does. A sub-cubic algorithm written out
  of smaller `fnp.matmul` calls is billed as the sum of its calls: one Strassen level = 7/8 of the matmul bill
  (plus 18 n^2/4 additions at weight 1), two levels 49/64 (ratio 0.766), three 0.67.
- `dtype_rate`: float32/float16/int32 and narrower = 1.0; **float64/int64 = 2.0**; float128 = 4.0.
  `astype` to float32 costs one copy (numel at the wider rate), then everything after is half price.
- `weight` tiers {0, 1, 4, 16}: 0 = views/metadata/zeros/empty (`transpose`, `diagonal`, `broadcast_to`, `split`,
  `zeros`); 1 = arithmetic, reductions, contractions, copies, `reshape`/`ravel`/`copy` (always billed numel),
  `concatenate`, `ones`, `eye`, `astype`, `sqrt`, `maximum`, `clip`, comparisons, `polyval`, `stats.norm.*`,
  all of `linalg.*`; 4 = `where(cond, x, y)` (4 x numel!), `take`, `sort`, `argsort`, `unique`, `searchsorted`,
  `bincount`, `histogram`, random permutations; 16 = `exp`, `log`, `sin`, `cos`, `tanh`, `power`, `arctan2`,
  `hypot`, `mod`, `floor_divide`, and every random *distribution* sampler (`random.normal` etc.).
- `complex_factor`: 1 for real dtypes (complex packing is priced out).

## Contractions (the whole bill of a chain)
    flop_cost = (2K - 1) * M      K = contracted size, M = output cells
- `matmul (n,n)@(n,n)`: 2n^3 - n^2  -> measured **4.00 n^3 in float64, 2.00 n^3 in float32**.
- `W @ C @ W.T` as two matmuls: 8 n^3 (f64) / 4 n^3 (f32).
- `fnp.einsum('ij,jk,lk->il', W, Csym, W)` with `Csym = flops.as_symmetric(C, symmetry=(0, 1))` and the same
  array object `W` in both outer slots: the engine proves the output symmetric and bills M = n(n+1)/2 for the
  second step: **6.00 n^3 (f64) / 3.00 n^3 (f32)**. Splitting it into two `matmul` calls loses the saving (8/4).
  So a float32 symmetric sandwich costs 3 n^3 = 3.2e9: the Gaussian pair chain is 16 x 3.2e9 = 5.2e10 = 0.024 B.
- Elementwise on n x n: `multiply` 1 per element (half on a symmetric tensor), `exp` 16, `stats.norm.cdf` 96 (f64) /
  96 (f32: stats kernels always compute in float64, 48 ops x rate 2), `stats.norm.pdf` 54, `polyval` deg d:
  2 d n^2 x rate, `where(C>0, C, 0)` 10 n^2 (compare 2 + where 8 in f64), `maximum(C,0)` 2 n^2 (f64).
  None of these matter next to the n^3 products unless called tens of times per layer.
- `outer(v, v)`: n(n+1)/2; `W @ v`: 2n^2; reductions: numel.

## Linear algebra (per matrix; float64 doubles everything)
`cholesky` n^3/3; `qr` 2(2mnk - 2k^3/3); `solve` 2n^3/3 + 2n^2 nrhs; `inv` 2n^3; `det` 2n^3/3; `eigh` 9n^3;
`eigvalsh` 4n^3/3; `eig` 25n^3; `svd` thin 6ab^2 + 20b^3; `svdvals` 2ab^2 + 2b^3; top-k `svd(k=)` min(4mnk, economy).
`eigh` of the 1024 x 1024 covariance is therefore 9 n^3 (f32) = 3 symmetric sandwiches: affordable once, not per layer.

## Wall time
Residual wall = wall - flopscope backend time - flopscope overhead; the Phase-2 gate is 0.4 s residual per predict,
120 s total per predict, 5 s setup, 8 GB. The einsum path search is unbilled but limited at >= 8 operands.

## Budget arithmetic for Phase 2 (B = 2^41 = 2.2e12, score multiplier max(0.1, C/B))
- The floor 0.1 B = 2.2e11 billed FLOPs = 205 n^3 over 16 layers = **12.8 n^3 per layer** = about 6 float32
  (n,n,n) matmuls, or 4 float32 symmetric sandwiches, or 3 float64 matmuls per layer. Below the floor nothing
  is gained, so a chain should spend exactly this.
- A float64 pair chain with plain matmuls (8 n^3/layer) already uses 0.06 B; the kappa_3 slices of a production
  chain (several more sandwiches per layer) push a float64 chain past the floor (504aldo's open chain: 0.25 B).
- Consequence for design: float32 + symmetric einsum + one Strassen level makes a sandwich 2.6 n^3, so about 5
  sandwiches per layer fit under the floor. Every extra n^3 object transported per layer costs ~0.02 B.
- The naive count in `k3chain3` (returned `fl`, FMA = 2, no dtype rate, no symmetry) overstates the float32
  symmetric bill by ~2.7x and understates a float64 plain-matmul bill by 2x: never read naive counts as billed cost.

## Precision
float32 arithmetic is 1e-7 relative per operation; over 16 layers of 1024-term sums the covariance is accurate to
~1e-6 relative, far below the 1e-4 rms error of a 1e-8-MSE estimator. Accumulations that cancel (the Laplacian
defect delta = e - Lap e, differences of near-equal slices) must be arranged to avoid catastrophic cancellation
in float32 or kept in float64 at twice the price.
