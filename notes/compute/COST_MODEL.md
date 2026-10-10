# The competition cost scorer, in full (flopscope 0.12.1 + whestbench 0.16.1, Phase 2)

Everything here was read in the source of the grader venv (`/opt/wb`) and measured with `flops.BudgetContext`
at n = 1024, L = 16. Probe scripts and raw outputs are in the session scratchpad `costmodel/{A,B,C,D}/` (A:
contractions and symmetry; B: pointwise, array and reduction ops; C: budget machinery, accumulation, timing;
D: linalg, FFT, random, stats). This file supersedes `FLOPSCOPE.md`; corrections to that file are marked [fix].
Units: n^3 = 1.0737e9, n^2 = 1.0486e6; f32 = float32, f64 = float64.

## 1. What is scored

| quantity | value | source |
|---|---|---|
| per-MLP score | `final_layer_mse * max(0.1, C_m / B)`; failures score with multiplier 1.0 on zero predictions | `whestbench/scoring.py:653-661, 958-971` |
| primary metric | **final layer only** (predictions must still be shape (16, 1024)) | `scoring.py:945-948` |
| aggregate | arithmetic mean of per-MLP scores | `scoring.py:1081` |
| B | 2^41 = 2.199e12 FLOPs per MLP | `budget.py` `PHASE2_ROUND` |
| C_m | exactly the billed FLOPs (lambda = 0: residual time is gated, not priced) | `budget.py:131-146` |
| failure gates | FLOPs > B (exception, before the op runs); wall > 120 s; **residual > 0.4 s** | `scoring.py:736-911` |
| setup() | unmetered (outside any context: billed to a 1e15 global default), 5 s cap, runs once per worker process (about 5-15 times per submission) | `subprocess_worker.py:180-210`, template |
| shipped data | up to 50 MiB / 50 files, pickle-free `.npz`, `fnp.load` costs 0 FLOPs | `limits.py`, ship-weights doc |
| allowed code | the grader's CPython, `flopscope` / `flopscope.numpy`, pure-Python stdlib. **No numpy, scipy, torch**; no threads, subprocess, FFI, compiled code | allowed-code doc, `cli.py:3843-3857` |
| fair accounting | a benefit that derives from how computation is accounted (rather than from the estimation method) can be invalidated; packing several values into one element is banned | allowed-code doc |

**The floor.** max(0.1, C/B) means everything below 0.1 B = 2.2e11 FLOPs = **205 n^3** is free in score. At the
floor the score is 0.1 x raw MSE, so beating v56 (adjusted 3.1e-9 at C/B 0.2) needs raw MSE below 3.1e-8, i.e.
RMS below 1.76e-4. Spend up to the floor on accuracy and never beyond unless the accuracy gain is proportionally
larger (above the floor, doubling FLOPs must halve raw MSE to break even).

A Monte Carlo sample is billed exactly `2 d w^2 + 17 w + 4 d w` = 3.36e7 (the harness's own formula,
`budget.py mc_flops_per_sample`): the floor buys about 6,500 samples, B about 65,000.

## 2. The billing equation

`charged = int(flop_cost * dtype_rate * complex_factor * weight)`, floored once, checked **before** the op runs
(an overshoot raises `BudgetExhaustedError`, charges nothing) (`_budget.py:1839-1916`).

- **dtype_rate** from `np.result_type` of the operands (NEP 50: Python scalars are weak) plus any `out=`: rate 1 for
  every dtype of 32 bits or fewer (bool, int8-32, f16, **f32**, c64), rate 2 for f64/int64/uint64/c128, 4 for f128.
  Nothing below f32 is cheaper. Promotion traps: `np.float64(c)` scalars, 0-d `fnp.array(2.0)` (f64 default),
  f32 x int32 array -> f64, constructors `ones/zeros/eye/full/array/linspace` default to f64, `arange` to int64.
- **complex_factor**: 1 for real dtypes (complex linalg x4, complex multiply x6).
- **weight** tiers from `data/default_weights.json` (472 entries; unknown names default to 1):
  - 0: views and layout (`transpose`, `.T`, `diagonal`, `broadcast_to`, `expand_dims`, `squeeze`, `split`, `flip`),
    `zeros`, `empty` (57 ops).
  - 1: arithmetic, `sqrt`, `abs`, `square`, `reciprocal`, `maximum`, comparisons, logical, reductions, contractions,
    every `linalg`/`fft`/`stats` op, copies (`copy`, `astype`, `array`, `reshape`, `ravel`, `concatenate`, `stack`).
  - 4: `where`, `take`, fancy indexing (wrapper x4), `sort`, `argsort`, `unique`, `searchsorted`, `bincount`,
    histograms, permutations.
  - 16: `exp`, `log`, trig, `tanh`, **`power` including every `x**k`** [fix], `mod`, `floor_divide`, `hypot`,
    `arctan2`, every random distribution sampler.
- **flop_cost** is the textbook count with FMA = 2 and a free first copy per output:
  `total = (k-1) * M + alpha - #outputs` (`_accumulation/_cost.py:354`). Values of the data never enter, except for
  result-sized ops (boolean mask indexing: mask size + 4 per selected element).

## 3. Contractions and symmetry (the n^3 terms)

| pattern (n = 1024) | f32 bill | notes |
|---|---|---|
| `A @ B`, `dot`, `einsum`, `tensordot` | 1.999 n^3 | (2K-1) M_out; f64 doubles everything |
| `inner(A, A)`, `einsum('ij,kj->ik', A, A)` (same object) | **1.0005 n^3**, tagged symmetric | `A @ A.T` and `X.T @ X` get **no** discount: `A.T` is a different object |
| `Cs @ Cs` (tagged) | 1.0005 n^3 | `Cs @ Cs.T` 2 n^3 |
| sandwich `einsum('ij,jk,lk->il', W, Cs, W)`, Cs tagged, same W object | **2.9995 n^3** | `W@C@W.T`, `multi_dot`, untagged C, `W.T` in a slot: 4 n^3 |
| sandwich via factor: `T = W @ L; inner(T, T)` | 2.9995 n^3, bit-exact symmetric | + Cholesky 0.333 n^3 if the factor must be made |
| diagonal of W C W^T, any formulation | 2.001 n^3 | 1.55 n^3 with a factor and Strassen-2 on `W @ L` |
| W diag(d) W^T: `einsum('ij,j,kj->ik', W, d, W)` | **1.0015 n^3**, tagged | `(W*d) @ W.T` 2 n^3 |
| rank-k update `inner(U, U)` | (2k-1) n(n+1)/2 | k = 64: 63.6 n^2 (`U @ U.T` 127 n^2) |
| thin `A @ B`, B (n, m) | (2n-1) n m | m = 64: 128 n^2 |
| Monte Carlo moments `einsum('bi,bj->ij', X, X)`, B rows | ~B n^2 | `X.T @ X` 2B n^2; `fnp.cov` f64 16 n^3 at B = 4n |
| 3-index T(w_i, w_i, w_i) for all rows, dense / S3-tagged | 2.0 / 1.0 n^4 (1.1e12) | infeasible; also 4 GiB storage |
| factored cubic `Z = W @ U.T; (Z*Z*Z) @ lam`, rank r | 2 n^2 r + O(n r) | r = 64: 128 n^2; the only workable cumulant form |
| diagonal (independent-input) third cumulant `(W*W*W) @ k3` | 4 n^2 | |

Rules (A, C):
1. **Three or more operands are billed along opt_einsum's `auto` path, whatever `optimize=` says** (`_cost.py:533`);
   each step takes the cheapest of three orbit formulas. Write a chain as one einsum to get the best order billed.
2. **Symmetry is seen only when declared (`SymmetricTensor`) or structural (the same Python object in two slots).**
   Untagged symmetric data, `W.T`, `W[:]`, `fnp.asarray(W)` earn nothing.
3. **[new] The f32 symmetric sandwich can raise `SymmetryError` after being charged.** The inferred-symmetric output
   is validated with `allclose(R, R.T, atol=1e-6, rtol=1e-5)` (`_einsum.py:963`, `_symmetric.py:146-162`); f32 round-off
   on outputs of size ~0.5 or larger fails it (measured: fails at median entry 0.48, passes at 0.19). Our covariances
   are O(1), so **use the factor form** (bit-exact) or f64 (6 n^3). Pre-scaling by 2^-24 to pass the check works
   but the source itself calls scale-to-pass a laundering vector: do not use it.
4. Pre-reduction: a summed index present in one operand only is reduced first (`einsum('ij,jk->i')` ~3 n^2).
5. Tags propagate through elementwise ops when **every** non-scalar operand is tagged; cost is then unique elements
   (0.5005 n^2 per op). Tags are dropped by any dense operand (including broadcast vectors `v[:,None]*S`), `where`,
   `astype`, `fnp.array`, `asarray`, `polyval`, `stats.*`, `reshape`, `triu`, slices that change an axis length,
   and writes (`fill_diagonal`, `copyto`). Tag sources at no extra cost: `outer(v, v)`, Gram einsums, `S @ S`,
   square `zeros/ones/full/eye`, `diag(v)`, unary ops on tagged arrays. Explicit tagging: `symmetrize(X, (0, 1),
   mode='canonical-copy')` n^2 (no validation) or `as_symmetric` 7 n^2.
6. No packed or triangular storage exists; the discount is in the bill only.

## 4. Pointwise, special functions, reductions (the n^2 terms)

| op (n x n, f32) | bill per n^2 | on a tagged matrix |
|---|---|---|
| `+ - * /`, `sqrt`, `abs`, `square`, `reciprocal`, `maximum(x, 0)`, comparisons | 1 | 0.5 |
| `exp`, `log`, `tanh`, **`x**2`, `x**k`, `x**0.5`** | 16 | 8 |
| `where(c, x, y)` (+ the comparison) | 4 (+1) | tag lost |
| `clip(a, lo, hi)` | 1 per bound | |
| `stats.norm.cdf` / `.pdf` / `.ppf` | 96 / 54 / 166, **always billed and returned in f64** (+2 to cast back) | no discount |
| **Phi by Abramowitz-Stegun 26.2.17, hand-built** (abs err 3e-7) | **35** | 17.5 |
| **phi = `exp(x*x*(-0.5)) * c`** | **19** | 9.5 |
| Phi and phi sharing one exp | 36 | 18 |
| `sum`, `max`, `cumsum` along an axis | 1 (N - #outputs) [fix] | 1 (no discount along an axis) |
| `mean` / `var` / `std` | 1 / 4 / 4 | |
| `sort`, `argsort` along an axis | 40 | |
| `polyval` deg d | 2d (f64 if coefficients are a Python list) | tag lost |
| hand-written Horner deg d | 2d | d |
| `outer(v, w)` / `outer(v, v)` | 1 / 0.5 | |
| `diag(A)` extract, `diagonal`, `.T`, slicing | 0 | |
| `reshape`, `ravel`, `copy`, `astype` | 1 (billed even when numpy returns a view) [fix] | |
| fancy indexing | 4 per gathered element | |

- **Not available**: `erf`, `erfc`, `ndtr`, `log_ndtr`, Owen's T, `numpy.polynomial` (no `hermval`), `vectorize`,
  scipy. Hermite polynomials must be built by recurrence; Gauss-Hermite nodes shipped as data.
- **Hermite series sum_k rho^k (c_k c_k^T)/k!, K = 10**: power form 172 n^2; repeated multiply 36 n^2; **Horner in a
  tagged rho with `outer(a_k, a_k)` terms: 14 n^2** (`a_k` must be the same object in both slots).
- Casts: `fnp.array(x64, dtype=f32)` 1/elem (half of `.astype`); a ufunc `dtype="float32"` folds the cast in for free
  (`fnp.exp(A64, dtype="float32")` 16 not 32; `fnp.positive(S64, dtype="float32")` 0.5 and keeps the tag).
- Arrays are immutable: `+=`, `a[i] = v` raise `TypeError`; `out=` never saves.

## 5. Linear algebra, FFT, random

| op (n = 1024) | f32 bill |
|---|---|
| `cholesky` | 0.333 n^3 (cheapest covariance factor) |
| `qr` (n, n) / (n, 64) | 2.67 n^3 / 1.64e7 |
| `eigh` / `eigvalsh` / `eig` | 9 / 1.33 / 25 n^3 |
| `svd` thin / `svdvals` | 26 / 4 n^3 |
| `svd(A, k=k)` (top k) | 4 n^2 k (k = 8: 3.4e7) |
| `solve` (n rhs) / `inv` | 2.67 / 2 n^3 (no triangular solve exists) |
| `norm(A, 2)`, `cond`, `matrix_rank` | 4 n^3 each (use `svdvals(A, k=1)`: 4 n^2) |
| `matrix_power(A, p)` | (floor(log2 p) + popcount(p) - 1) matmuls |
| FFT complex length n | 5 n ceil(log2 n), rfft half |
| `standard_normal(dtype=f32)` / f64 / `normal` | 16 / 32 / 32 per element |
| `random(dtype=f32)`, `integers(int32)` | 1 per element |
| `permutation(n)` | 4 n |
| Haar frame m = 1024 by QR of a Gaussian | 2.68 n^3 (ship it as data instead: 0) |
| `multivariate_normal` | **avoid**: 26 d^3 SVD in f64 |

Failed linalg calls are still billed. `hermitian=True` never changes a bill.

## 6. Time: the 0.4 s residual gate and the call budget

- residual = wall - backend - overhead (`_budget.py:1579-1591`). Backend is time inside numpy calls; overhead is
  flopscope's own dispatch (free, 70-270 us per call, counts only toward the 120 s cap).
- **Every counted flopscope call leaks 20-45 us into residual** (wrapper bookkeeping after the wall sample,
  `_budget.py:1210-1234`). About 10^4 calls per predict exhaust the gate on their own. Write per-layer work as a
  few dozen large calls; never loop over neurons or samples in Python. Batched and looped forms bill identical FLOPs.
- Freeing large temporaries is residual too (20 x 16 MB: 3.6 ms).
- The remote (AIcrowd client/server) path was not probed; per-call latency there may differ.

## 7. Policy classification of cost-reducing techniques

| class | techniques |
|---|---|
| intended, use freely | f32 everywhere; one einsum per chain; same-object Grams; symmetric tags; factor-form sandwich; W diag(d) W^T einsum; Horner with `*`/`+`; hand-built Phi/phi; `x*x` instead of `x**2`; `maximum` instead of `where`; data shipped in the submission (frames, tables, nodes); spending up to the 0.1 B floor |
| legitimate algorithm, but confirm with the organisers before relying on it | Strassen recursion written with `fnp.matmul` (genuinely fewer multiplications; v56 uses six levels, (7/8)^6 = 0.45 of a matmul) |
| pricing gaps: do not use (report instead) | `svd(A, k=n-1, hermitian=True)` as a 4 n^3 full eigendecomposition; `fnp.random.symmetric` as a 2 FLOP/element Gaussian sampler; scaling by 2^-24 to pass symmetry validation |
| prohibited | raw numpy or scipy, `np.asarray(fA)` arithmetic, threads, subprocess, FFI, compute in residual time, packing values into wider elements |

## 8. What this means for the Stage 21 system (budget map, f32, n = 1024)

| component (per predict) | cheapest faithful idiom | billed |
|---|---|---|
| Gaussian-part forward sweep, 16 layers: means `W m` | matvec | 32 n^2 |
| covariance transport `C_z = W K_h W^T` | factor form: `L = cholesky(K_h)`, `T = W @ L`, `inner(T, T)` | 3.33 n^3 x 16 = 53 n^3 |
| post-activation covariance and cut kernel (Hermite series, K = 10, both) | Horner in a tagged rho, `outer(a_k, a_k)` terms, hand-built Phi/phi on vectors | ~30 n^2 x 16 = 0.5 n^3 |
| final-layer per-neuron variances (diagonal only) | `einsum('ij,jk,ik->i')` or row norms of `W @ L` | 2 n^3 (already inside the sweep) |
| per-neuron third cumulants of the final preactivation, rank-r factored sources | `Z = W @ U.T; (Z*Z*Z) @ lam` | 2 n^2 r per layer carried (r = 64: 0.125 n^3) |
| backward (Heisenberg) observables, carried vectors | matvecs `W^T diag(Phi) v` | 2 n^2 per vector per layer |
| spectral diagnostics: gaps g_l, flatness along given directions, undecided cut | vector ops | O(n^2) |
| orthogonal probe block of m inputs (collective cumulants, D3) | shipped Haar frame (0 FLOPs) x shipped radii, forward in batch, `einsum('bi,i->b')` projections | 32 n^2 m: m = 1024 is 32 n^3 = 0.016 B |
| restart window (Thm 5.1): exact recomputation of the last l layers | as the sweep, plus factored births | about 3.5 n^3 per layer: l = 10 is 35 n^3 |
| **total** | | **about 120 n^3 = 0.06 B, under the 0.1 B floor with about 85 n^3 to spare** |

Call count for the whole system is a few hundred to a few thousand counted calls (residual well under 0.4 s).

**Measured, not estimated** (`whest/s21_fnp.py`, `scripts/s21_bill_pass1.py`, net 0, warm worker, real
`BudgetContext`): the pass-1 sweep bills **51.3 n^3 = 0.0251 B** (matmul 30.0, inner 16.0, cholesky 5.0 n^3; all
n^2 work together 0.3 n^3), 2,137 counted calls, residual 0.072 s, wall 2.15 s. Its f32 output agrees with the f64
numpy prototype to 1.9e-6 RMS at every layer (final-layer error against truth 2.016e-3 in both), so f32 round-off is
three orders of magnitude below every accuracy effect discussed in notes/stage21.

Consequences for reading experiments:
1. Accuracy measured with a float64 numpy prototype transfers unchanged to an f32 flopscope implementation only if
   f32 round-off is below the effect being measured (it is ~1e-7 relative per op; check cancellation-heavy
   differences in f32 explicitly).
2. A design is never "too expensive" because of the n^2 terms: they are 1-3% of the n^3 terms when written with the
   idioms above. The n^3 count per layer is what must be planned.
3. Monte Carlo and probe blocks are cheap only in batch; anything per-sample in Python fails the residual gate.
