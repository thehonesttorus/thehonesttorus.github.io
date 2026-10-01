# Porting kit: measured prices (flopscope 0.12.1, n = 1024 unless stated)

*Measured 2026-10-01 by `probe_kit.py` (results/kit_costs.json), in-process meter, 2 BLAS threads, box otherwise
quiet. 1 unit = 2^31 FLOPs (one dense 1024³ float32 product); B = 1024 units. "calls" = logged flopscope ops (the
residual currency, see RESIDUAL.md). "resid ms" = Python-side residual of one invocation in-process; "backend ms" =
numpy/BLAS wall (counts against the 120 s wall cap, not the residual). Helpers named here are in `kit.py`.*

## Reading the table: what matters for a design

- **Products.** Dense f32 is 2mkn − mn FLOPs, any float64 operand bills the whole op 2× (and mixed f32@f64 is
  also slow). Batching gives no discount, only fewer calls. Small batched products are cheap in units and calls.
- **Same-object Gram is half price**: `einsum('ji,jk->ik', X, X)` 0.50 u vs `X.T @ X` 1.0 u (a transposed view is
  another object). It costs ~1.5 ms residual (the symmetric-output machinery) and returns a tagged array.
- **The tagged sandwich einsum raised `SymmetryError` here on an SPD C at n = 1024** (unit-scale inputs, float32
  output check `allclose(atol=1e-6)`). Treat `einsum('ij,ia,jb->ab', C_tagged, W, W)` as unusable; the safe
  sandwich is 2 dense products (2.0 u) or Strassen (2 × 0.557 u at L5).
- **Strassen–Winograd** (kit.Strassen, the lower-op engine): price per 1024³ product L1…L5 = 0.877 / 0.772 /
  0.683 / 0.610 / 0.557; calls per family **independent of the batch** (L1 16, L2 31, L3 46, L4 64, L5 79; hub 87).
  It is much slower in backend wall than BLAS (L5 one product 267 ms vs 21 ms dense; a batch of 4 at L5 1.1 s):
  budget the 120 s wall cap if you use many deep families.
- **Gaussian elementwise.** `stats.norm.cdf/pdf` bill float64 (96 / 54 FLOPs per element as billed): negligible
  on (n,) vectors, 0.047 u per (n, n) matrix plus ~2 ms residual; keep gate statistics at vector level when you can.
  ReLU moments of a vector of Gaussians (`kit.relu_moments`) cost 1e-4 u but 18 calls.
- **Factorizations at n = 1024 (f32):** Cholesky 0.167 u, slogdet/det 0.333 u, solve with n right-hand sides
  1.33 u (one rhs 0.334 u), inv 1.0 u, QR 1.33 u (n × 256: 0.115 u), eigvalsh 0.667 u, **eigh 4.5 u**, **svd 13 u**.
  All one call. float64 doubles. Batched small factorizations are nearly free: (4096, 8, 8) slogdet 0.0009 u,
  (1024, 32, 32) eigh 0.14 u, (64, 128, 128) eigh 0.56 u, each one call.
- **No LU and no Pfaffian primitive** in `fnp.linalg`: det/slogdet use LU internally; a Pfaffian must be built from
  a determinant (Pf² = det, sign by continuity or a skew-LU written by hand) or a Householder tridiagonalisation.
- **Views are logged calls**: reshape (also billed 1 FLOP per element), transpose, swapaxes, moveaxis, diagonal.
  Basic slicing is not logged in-process (but is a client/server round trip on the grader: RESIDUAL.md).
- Gathers: `fnp.take` bills 4 FLOPs per gathered element (0.002 u for n rows of n).

Other measured forms (notes/streams/costmodel/REPORT.md §2): weighted Gram with sign-indefinite weights via two
split-sign aliased Grams 0.503 u in 9 calls; Hadamard-then-contract `einsum('ij,ij,j,jk->ik')` 1.0 u in one call;
rank-r families of sandwiches 1.5–2.0 r units; `norm.ppf` 166 FLOPs/element (float64); `astype` 1–2 per element.

## Table

| primitive | units | calls | resid ms | backend ms |
|---|---|---|---|---|
| matmul (n,n)@(n,n) f32 | 0.9995 | 1 | 0.109 | 20.9 |
| matmul out= pooled | 0.9995 | 1 | 0.092 | 21.0 |
| matmul f64 | 1.9990 | 1 | 0.128 | 43.8 |
| matmul f32 @ f64 (promotes) | 1.9990 | 1 | 0.975 | 63.1 |
| matmul (n,n)@(n,64) | 0.0625 | 1 | 0.091 | 2.3 |
| matmul (n,n)@(n,256) | 0.2499 | 1 | 0.090 | 5.2 |
| matvec (n,n)@(n,) | 0.0010 | 1 | 0.068 | 0.2 |
| batched (8,n,n)@(8,n,n) | 7.9961 | 1 | 0.615 | 187.4 |
| batched small (4096,8,8)@(4096,8,8) | 0.0018 | 1 | 0.062 | 1.0 |
| batched small (1024,32,32)@(1024,32,32) | 0.0308 | 1 | 0.120 | 2.3 |
| A.T @ A (view: full price) | 0.9995 | 1 | 0.131 | 17.7 |
| gram_same einsum('ji,jk->ik',A,A) | 0.5002 | 1 | 1.454 | 19.5 |
| sandwich W^T C W (2 dense) | 1.9990 | 2 | 0.161 | 39.9 |
| diag_sandwich | 1.0005 | 3 | 0.133 | 11.6 |
| tagged einsum sandwich (SPD C) | **raises** `SymmetryError: Tensor not symmetric along axes (0, 1): max d…` | | | |
| Strassen L1 one product | 0.8768 | 16 | 0.794 | 22.8 |
| Strassen L2 one product | 0.7715 | 31 | 1.454 | 47.3 |
| Strassen L3 one product | 0.6829 | 46 | 2.806 | 103.9 |
| Strassen L4 one product | 0.6096 | 64 | 3.668 | 158.5 |
| Strassen L5 one product | 0.5567 | 79 | 5.544 | 266.9 |
| Strassen L5 batch 4 (shared left) | 2.1808 | 79 | 5.941 | 1093.3 |
| Strassen L3 batch 4 (shared left) | 2.7168 | 46 | 3.172 | 284.6 |
| Strassen hub L5 sum_4 X_k Y_k^T | 2.1870 | 87 | 5.733 | 968.9 |
| add (n,n) out= | 0.0005 | 1 | 0.035 | 1.2 |
| multiply (n,n) by row vector out= | 0.0005 | 1 | 0.045 | 0.9 |
| add (n,n) fresh result | 0.0005 | 1 | 0.038 | 1.0 |
| exp (n,n) | 0.0078 | 1 | 0.031 | 0.9 |
| sqrt (n,) | 0.0000 | 2 | 0.039 | 0.0 |
| norm.cdf (n,) +astype | 0.0000 | 2 | 0.114 | 0.0 |
| norm.pdf (n,) +astype | 0.0000 | 2 | 0.059 | 0.0 |
| norm.cdf (n,n) | 0.0469 | 1 | 2.071 | 0.0 |
| relu_moments (n,) | 0.0001 | 18 | 0.837 | 0.1 |
| sum axis 0 (n,n) | 0.0005 | 1 | 0.062 | 0.6 |
| einsum 'ij,ij->i' | 0.0010 | 1 | 0.057 | 0.9 |
| take rows (n idx) | 0.0020 | 1 | 0.069 | 0.8 |
| reshape (n,n)->(2,512,2,512) | 0.0005 | 1 | 0.034 | 0.0 |
| transpose view | 0.0000 | 1 | 0.035 | 0.0 |
| diagonal view | 0.0000 | 1 | 0.032 | 0.0 |
| copyto (n,n) | 0.0005 | 1 | 0.071 | 1.0 |
| cholesky n f32 | 0.1667 | 1 | 0.115 | 34.7 |
| cholesky n f64 | 0.3333 | 1 | 0.056 | 32.0 |
| eigh n f32 | 4.5000 | 1 | 0.075 | 215.6 |
| eigvalsh n f32 | 0.6667 | 1 | 0.075 | 123.2 |
| qr n f32 (reduced) | 1.3333 | 1 | 0.074 | 217.4 |
| qr (n,256) f32 | 0.1146 | 1 | 0.072 | 46.4 |
| solve n, n rhs f32 | 1.3333 | 1 | 0.059 | 120.8 |
| solve n, 1 rhs f32 | 0.3343 | 1 | 0.065 | 33.5 |
| slogdet n f32 | 0.3333 | 1 | 0.070 | 33.4 |
| slogdet n f64 | 0.6667 | 1 | 0.093 | 33.4 |
| det n f32 | 0.3333 | 1 | 0.072 | 34.3 |
| inv n f32 | 1.0000 | 1 | 0.058 | 100.2 |
| svd n f32 | 13.0000 | 1 | 0.062 | 614.5 |
| slogdet batched (4096,8,8) | 0.0009 | 1 | 0.060 | 3.7 |
| cholesky batched (4096,8,8) | 0.0003 | 1 | 0.051 | 2.6 |
| solve batched (4096,8,8)@(4096,8,8) | 0.0026 | 1 | 0.060 | 12.0 |
| eigh batched (1024,32,32) | 0.1406 | 1 | 0.063 | 165.7 |
| slogdet batched (1024,32,32) | 0.0107 | 1 | 0.071 | 17.0 |
| cholesky batched (64,128,128) | 0.0208 | 1 | 0.084 | 25.0 |
| eigh batched (64,128,128) | 0.5625 | 1 | 0.073 | 151.7 |