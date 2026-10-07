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

## 4. Gate saturation: which neurons the source state can drop (`V33_SAT`, emulated; `outputs/sat_*.txt`)

The retained-transport PDF's saturation certificate says a gate deep in its off state passes almost nothing: a
neuron's row of every leg enters the next transport weighted by w1 = Phi(alpha), and its own reads move a mean of
order phi(alpha). If such rows could be dropped, the transports, the D21 contractions and the old-tier formings, about
80% of the bill, would shrink with the active count. Emulated by w1 -> w1 1[alpha > alpha_c] in every transport and
zeroed D21 rows at the read (networks 0-15, against the current best):

| drop alpha <= | rows dropped, layer 3 -> 15 | raw |
|---|---|---|
| -2.5 | 2% -> 23% (about 11% on average) | +0.22% +- 0.33 |
| -2.0 | 5% -> 28% | +13.9% +- 1.1 |
| -1.5 | 10% -> 34% | +125% |
| -1.0 | 21% -> 38% | +688% |

Only the deep tail is free: at alpha <= -2.5 (w1 <= 0.006) the drop is exact within noise, and a quarter of the rows
at depth go with it. Below -2 the damage is steep. Section 6 separates the transport rows from the reads.

## 5. The (2,2) slice against Monte Carlo truth (`code/k4mc2.py`, `code/k4w_an.py`, `outputs/k22_slice_truth_off0.txt`)

4e6 inputs through network 0 (two halves for the noise), the true pre-activation slice kappa(y_i, y_i, y_j, y_j) at
layers 4, 7, 10, 13 against the chain's wk4m and three constructions (off-diagonal entries; relative error, with the
Monte Carlo's own noise at 0.11, 0.08, 0.06, 0.05 of the signal):

| target layer | chain (g4_i + g4_j)/6 | exact pair class, true K22(x) | pair class, chain's K22(x), rank 4 | scale mixture g(v_i v_j + 2 C_ij^2) |
|---|---|---|---|---|
| 4 | corr 0.537, err 0.168 | 0.376, 0.181 | 0.373, 0.183 | 0.544, 0.165 |
| 7 | 0.602, 0.187 | 0.332, 0.214 | 0.332, 0.217 | 0.612, 0.183 |
| 10 | 0.653, 0.224 | 0.312, 0.263 | 0.312, 0.271 | 0.663, 0.217 |
| 13 | 0.715, 0.292 | 0.229, 0.348 | 0.230, 0.362 | 0.738, 0.276 |

- **The quenched pair class is not the (2,2) slice.** Transported exactly through the real weights from the true
  post-activation slices, it is worse than the chain's annealed shape at every layer. Unlike the diagonal (note XXVIII),
  the off-diagonal (2,2) slice is carried by the classes with three and four distinct indices, which the pair class
  leaves out: the mixture's own structure, not per-neuron quenched content.
- **The chain's shape is close to the mixture's.** Its arithmetic-mean shape and the mixture's product shape agree to
  a few percent; the mixture's 2 g C_ij^2 term is the one structured piece the chain lacks. The chain's residual
  correlates with C_ij^2 increasingly with depth (0.07, 0.10, 0.19, 0.42), with three times the mixture's amplitude at
  layer 13: the (2,2) slice carries more C^2 content than a scale mixture of one gain does.
