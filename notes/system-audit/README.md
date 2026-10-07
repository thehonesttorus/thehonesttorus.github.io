# System audit: the cost rules, where the bill goes, and every annealed average

Working note XXIX. A full audit of the current best estimator (`est_v33`: the V32 chain plus the quenched kappa4 pair
class, `notes/retained-transport`) on two axes: computation (is anything naively implemented under flopscope's
billing, and which representational choices dominate the bill) and theory (every place the chain replaces the
actual weights by a weight-averaged form, with its leverage and its quenched remainder measured).

## 1. The billing rules that matter here (flopscope 0.12, `docs/reference/cost-model.md`)

`charged = int(flop_cost x dtype_rate x complex_factor x weight)`.

| rule | consequence for the chain |
|---|---|
| matmul (m,k)(k,n) bills 2mkn - mn; einsum bills its contraction plus a `copyto` of numel(out) when `out=` is given | use `matmul(out=)`, not `einsum(out=)`, for every contraction with a destination |
| symmetric output savings only when operands alias (`inner(A, A)`, `outer(v, v)`) or carry an `as_symmetric` tag; `A @ A` bills full | Grams as `einsum("ia,ib->ab", X, X)` with aliased operands bill half |
| every write bills >= 1 per element: `copy`, `reshape`, `ravel`, `astype`, `concatenate`, `copyto`, `ones`, `full`, `fill_diagonal` (min(m,n)); views (`T`, slicing, `swapaxes`) and `zeros`/`empty` are free | buffers are pooled and written `out=`; zero-filled allocation is free, constant fills are not |
| float64 bills 2x; `stats.norm.pdf/cdf` always compute in float64 (27 and 48 FLOPs per element, x2); a Python float scalar as a `copyto` source resolves to float64 | the page-touch `copyto(b, 0.0)` bills at float64 (1.7 units per process); the stats calls are O(n) and immaterial |
| transcendentals weight 16; gathers and `where(c, x, y)` weight 4 | none on any n^2 or n^3 path |
| QR (reduced) 2(2mnk - 2k^3/3); eigh 9n^3; thin SVD 6ab^2 + 20b^3 | the range-finder QRs (1024 x 320) are 0.18 units each, 3.7 units per network |
| Strassen written as ordinary ops bills its analytical count | the chain's dominant saving; depth is set by the 0.4 s residual cap, not by FLOPs |

Phase 2 scores C = F (lambda = 0) with the residual wall time gated at 0.4 s per network: crossing it fails the
network.

## 2. Where the bill goes (`code/prof_chain.py`, `outputs/profile_v33_net0.txt`)

Every billed op was attributed to its source line through the Python stack at billing time; helpers and the Strassen
class roll up to their call site in the layer loop. Network 0, current best, one cold predict (223.9 units of 2 n^3,
0.2186 B):

| component | units | share |
|---|---|---|
| `_dslices`: the D3/D21 contractions of every source (young dense hub, old tier in factor space, thin legs) | 88.1 | 39.4% |
| transport of the young sources' dense legs, W [A, P] | 65.8 | 29.4% |
| old tier: forming dense legs from the rank-320 basis and the rank-192 nested tier | 22.3 | 10.0% |
| covariance transport C_pre (symmetric Strassen family) | 6.7 | 3.0% |
| joins: range finders, QRs, rotations, the two-step metric | about 25 | about 11% |
| pair-term programs (one fused einsum) | 2.5 | 1.1% |
| all elementwise work, copies, reshapes, reductions | under 10 | under 4% |

By op: matmul 89%, the Strassen block additions about 4%, QR 1.7%. The bill is dense contractions of the source
representation; no naive implementation of any size remains (the V26-V29 rounds removed them). The only free waste is
the float64-rate page touch (`copyto(b, 0.0)`, 1.7 units, first predict of a process only).

**The scored regime is cheaper than our harness.** whestbench's worker loads the estimator once and serves every
network of a run through the same instance (`subprocess_worker.main`, `cli` lines 2496-2511), so every network after
the first runs with the pools warm and at the full Strassen depth (`V26_STRASSEN=6`; the first predict is capped at
`STRASSEN_FIRST=4`). Our harness ran one predict per process. Re-measured with a warm-up predict first
(`code/run_v29w.py`, `outputs/scored_regime_v33_16nets.txt`): C/B **0.1986** against 0.2186, raw identical. Every
C/B in notes XIX-XXVIII is the first-call figure; the scored adjusted error is about 9% lower than reported.

**The residual cap is the binding risk.** The same runs show 0.52-0.54 s of residual Python time per measured predict
and 0.66 s on the first, against the 0.4 s gate, with 16 processes sharing the VM. Section 4 measures it on an idle
machine.

## 3. Leverage of every annealed (weight-averaged) construction (networks 0-15, paired, `outputs/abl_*.txt`)

| construction in the kappa4 sector | how it is annealed | removed or zeroed: raw change |
|---|---|---|
| (2,2) pre-activation slice wk4m = (g4_i + g4_j)/6 | a rank-two sum shape from the transported diagonal | **+1539%** |
| lambda s_off^2 stand-in for the (2+1+1) class | one fitted scalar per layer times the off-diagonal variance | +88.7% |
| (3,1) slice lambda C_off | the same scalar times the covariance | +14.7% |
| kappa4 -> kappa3 feed (hub pair X3, Y3, M_t1) | built from the mean-field core dG | +4.3% |
| adaptive lambda rule (BETA) | a ratio of means | +0.3% +- 0.2 |
| matrix-core projection g_prev (cA, cI) | row sums with fresh-He constants | quenched part added on the diagonal by V33 (-2.25%) |

The (2,2) slice carries by far the most leverage and is the crudest shape, so its quenched remainder is measured first
(section 5).
