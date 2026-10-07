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
at depth go with it. Below -2 the damage is steep, and it is mostly in the reads:

| drop | raw | better on |
|---|---|---|
| alpha <= -2.0, transport rows only | +1.81% +- 0.67 | 3/16 |
| alpha <= -2.0, D21 rows only | +7.79% +- 0.74 | 0/16 |
| alpha <= -1.5, transport rows only | +26.3% +- 1.8 | 0/16 |
| alpha <= -1.0, transport rows only | +159% | 0/16 |

A saturated neuron's own (2,1) reads matter more than its forward transport, and both turn on steeply between -2.5 and
-2. The cost lever is therefore the alpha <= -2.5 drop of both, which removes about 11% of rows on average (23% at
layer 15) from the transports, the D21 contractions and the old-tier formings, about 75% of the bill. Rounded to
Strassen-compatible sizes (multiples of 64, leaf side >= 14) that is 5-6% of the bill.

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

**In the chain the truer shapes lose** (`V33_WK4M`, networks 0-15, paired):

| wk4m | raw | better on |
|---|---|---|
| scale mixture sqrt(g_i g_j)(v_i v_j + 2 C_ij^2), g from the chain's own diagonal | +4.80% +- 0.71 | 0/16 |
| chain's (g4_i + g4_j)/6 + 2 sqrt(g_i g_j) C_ij^2 | +12.08% +- 1.05 | 0/16 |
| geometric mean sqrt(g4_i g4_j)/3 | -0.61% +- 0.42 | 11/16 |

The structure that matches the truth better offline makes the chain worse: the fitted fourth-cumulant sector
compensates the dropped gate-covariance and second-order terms (note XXI section 7), and a truer (2,2) shape breaks
that balance. The (2,2) slice is closed; its quenched remainder is not pair-class content, and its mixture refinements
do not transfer.

## 6. Randomized evaluation of the residual operators, measured (`code/hutch_var.py`, `outputs/hutch_variance_off0.txt`)

A further synthesis (`THEORY.md`, "quenched residual queries") proposes keeping the chain and evaluating each
fixed-weight residual operator either exactly at low rank or by sign probes, Z = eps o (A eps), E Z = diag A,
Var Z_i = sum_(j != i) A_ij^2, and states the acceptance rule as max(0.1, C_new/B)(M_1 + V/K) < max(0.1, C_0/B) M_0.
Its decisive quantity, the off-diagonal energy of the operator, was measured on the chain's stored K22 at four
layers (network 0):

| operator (diagonal = the wanted correction) | probes for SNR 1 | units per layer at SNR 1 |
|---|---|---|
| pair class 3 (W o W) K22 (W o W)^T | 1015-1018 | 3.0 |
| deflated by the row means of W o W (the quenched operator) | 414-486 | 1.2-1.4 |
| remainder after the exact rank-4 part (V33) | 132-320 | 0.4-0.9 |

A weight-sandwiched operator has off-diagonal rows of n entries each of the diagonal's own size, so the probe count
for SNR 1 is of order n before deflation and a few hundred after; a useful 10% error costs a hundred times that. The
exact rank-4 contraction that V33 ships costs 0.05 units per layer with no noise. The synthesis's ordering (exact
low-rank first, probes only for a cheap low-variance remainder) is the right one, and on these operators it ends at
the exact branch.

## 7. Exact billing fixes (`V34_OPT` letters; networks 0-15, scored regime, paired)

Every item below computes the same quantity in fewer billed FLOPs or fewer ops; raw MSE moves only at float rounding
(per-network differences up to 0.16%, mean within +-0.03%).

| letters | what changes | raw | scored C/B |
|---|---|---|---|
| a, b, d | pair-term programs as one einsum per output group (no one-hot group contraction); the join Grams `Sj`, `S_s`, `Qp^T Qp` as aliased two-operand einsums (half billed); thin feedback products as `matmul(out=)` (no `einsum(out=)` copy) | -0.005% +- 0.021 | 0.1986 -> 0.1951 |
| + c | under `JOIN_POST=3` the range finder's QR is skipped: only span(W Qn) is used, and that block orthonormalises it | +0.013% +- 0.017 | -> 0.1934 |
| + e | column r+1 of every L leg is identically zero; the M-leg thin contractions run over r+1 columns | rounding | -> 0.1931 |
| + f, g | slot 0 (the layer-0 source) has identically zero feedback legs and feed vectors, so those blocks run over slots 1..k; D3's first term read off `A*A*w2` (already formed for LP), the feed's `R = rowsum(X3*P)` in one pass, `Xt*P` once | +0.024% +- 0.022 | -0.33% (with section 9) |
| + j | the joiner's A and P legs as one batch: range finder = one mm + one hub, factors = one mm, instead of six products | -0.026% +- 0.020 | -0.06%; 39 ms less local residual |

Files: `outputs/v34_*.txt`, `v35_cand_abcde(fg).txt`, `v35_j*.txt`.

**Steady state.** The measured warm call (predict #2 of a process) still pays the one-time allocation of the level-6
pools; predict #3 onward (`W_WARM=2`) costs 0.1899 B with a-e against 0.1931 at predict #2 (`v34_abcd_predict3_memory.txt`,
`profile_v34_predict3_net0.txt`: 194.4 units, of which the young D21 hub leaf 55.8, the young transport 54.2, the
old-tier formings 21.6, the joins about 20).

## 8. Memory and where the arrays live

The in-process harness peaks at VmHWM 8.04 GiB / VmPeak 9.39 GiB (level-6 pools). whestbench's `SubprocessRunner` sets
`RLIMIT_AS` = 8192 MB before loading the estimator (`code/run_worker.py` drives that exact worker). Run against the
in-process flopscope, it fails: with 96 BLAS threads predict #1 already raises `MemoryError` (VmHWM 4.35 GiB, the rest
of the address space is thread reservations); with 4 or 16 threads predict #1 (level 4) passes and predict #2 (level 6,
and level 5 too at 16 threads) fails (`outputs/worker_rlimit_*.txt`).

This is not the grader's situation. The grader runs the estimator against the **flopscope client**: arrays are
`RemoteArray` handles and the data lives in the flopscope backend ("the instance has 64 GB in total, the rest going to
the flopscope backend and the harness", `docs/reference/estimator-contract.md`; the client's `_remote_array.py`: "the
data lives on the remote grading server, not in this process"). The 8 GB cap bounds a process that holds no array
data, and the graded V28 (level 5, shared basis) is consistent with that. Two consequences: a local
`--runner subprocess` check cannot validate level 6, and on the grader `residual = wall - client dispatch`, so an op's
encode, transport and server time (page faults included) is overhead, not residual; only our own Python between ops
counts against the 0.4 s gate.

**First-call items, not shipped.** The page touch (`copyto(b, 0.0)` billed at the float64 rate) and the `_buf2`
reshape cost 2.77 units on predict #1 and 3.29 on predict #2; the level cap of predict #1 (`STRASSEN_FIRST=4`) costs
about 23.8 units. Over a 50-MLP run that averages about 0.6 units per network (0.3%). Lifting the cap could put
network 1 over the residual gate, which zeroes its prediction and dominates the mean, so the cap stays.

## 9. The saturation drop made real (V35: `V33_SAT=-2.5 V35_SAT_ROUND=64 V35_SATC=1`)

The active set of layer l is the top `na` neurons by alpha = mu/sigma, with `na` = the count above -2.5 rounded up to
a multiple of 64 (so the dropped set is a subset of alpha <= -2.5, and every family keeps Strassen-compatible sides).
Two families shrink:

- **the young D21 hub** forms only the active rows: LA/LP rows gathered in alpha order (`take`, O(n^2)), a hub with
  m = na, and the result put back in neuron order with the dropped rows zero;
- **the young transport** at layer l+1 contracts only over layer l's active rows: at the end of layer l the
  transported legs (w1-scaled, so their dropped rows are zero) are gathered on the active rows into the spare
  ping-pong side, and the family runs `W[:, perm[:na]] @ legs[:na]`. The covariance rides in the same family with
  all rows permuted; its dropped rows' part `W[:, perm[na:]] C[perm[na:]]` is added densely.

The Strassen pools are keyed by the full block shape and handed out as sliced views, so the compacted families run
in the buffers the full families own (no new memory); the leaf minimum of these families is 8 (`V35_SAT_MN`), which
keeps level 6 at na = 832..960.

| (networks 0-15) | raw | C/B | adjusted |
|---|---|---|---|
| emulated mask (same rows zeroed, full families) vs no drop | +0.144% +- 0.164 | 0.1931 | +0.15% |
| compact vs emulated | -0.010% +- 0.024 | 0.1931 -> 0.1852 | **-4.12%** |
| compact vs no drop | +0.134% +- 0.168 | 0.1931 -> 0.1852 | -3.98% |

The compact families reproduce the emulation to rounding; the drop itself is neutral within noise. Same-VM residual
cost: +10 ms locally (`v35_pair_vm2_*.txt`). With SAT = -2.25 the bill falls another 0.9% but raw rises
+0.51% +- 0.37 (adjusted -0.4%, within noise; not adopted).

**The full drop is not available.** Also removing a dropped neuron's own D3 entry and its D21 column (so that nothing
at layer l reads its row of any leg, which would let the transports' output rows and the formings shrink too) costs
raw +12.2% +- 1.4 (`v35_sat_full_drop.txt`). A saturated neuron's forward fluctuation can go, but its own mean is
carried by the tail, where the third-cumulant corrections are large relative to it: its marginal reads must stay.

## 10. Leaf minimum 8 for the narrow families (`V35_SB_MN`, `V35_CPRE_MN`, `V35_JN_MN`)

The rank-320/192 shared-basis families (formings, rotations, factor hubs), the C_pre symmetric family (blocks of side
512) and the join's single products stopped one Strassen level short at the leaf minimum 16. At minimum 8 each runs
one level deeper; exact (raw -0.004% +- 0.023):

| (with j, networks 0-15) | C/B | adjusted | local residual |
|---|---|---|---|
| no change | 0.1845 | | |
| formings, rotations, factor hubs, C_pre at minimum 8 | 0.1810 | -1.88% | about +0 ms net of j |
| also the join products and the lift | 0.1799 | -2.49% | +22-24 ms (two VMs agree) |

Without j the same change costs +55 ms; j's batching pays most of it back.

## 11. Schedule: no join at li = L-2 (`V35_SKIP_JOIN_L2`)

The li = 14 join and nest compress a source that is read for two more layers only. Keeping it dense instead: raw
**-0.43% +- 0.14**, better on 13 of 16 networks, C/B unchanged (the audit's predicted -1.5 units did not
materialise once the compaction landed). Adopted as an accuracy gain at equal cost.

## 12. Validation on all 100 networks (scored regime, paired; `outputs/v35_100nets_*.txt`)

Candidate = V33 + `V34_OPT=abcdefgj V33_SAT=-2.5 V35_SAT_ROUND=64 V35_SATC=1 V35_SKIP_JOIN_L2=1 V35_SB_MN=8
V35_CPRE_MN=8 V35_JN_MN=8` (`code/est_v35.py`), against V33 in the same harness (warm-up predict on another network,
then the measured predict):

| | raw (mean) | C/B | adjusted (mean) | first-call C/B | local residual, measured call (mean / max) |
|---|---|---|---|---|---|
| V33 | 2.2939e-08 | 0.1986 | 4.5553e-09 | 0.2186 | 0.497 / 0.623 s |
| candidate | 2.2899e-08 | **0.1802** | **4.1273e-09** | 0.1983 | 0.524 / 0.613 s |

Raw -0.16% +- 0.10 (better on 58 of 100), adjusted **-9.4%**. From predict #3 on the bill is about 0.003 B lower
again (section 7), so a 50-MLP run averages about C/B 0.177. The local residual rises 5%; on the grader only our own
Python between ops counts (section 8), and the graded V28 measured 0.17 s there.

What this round did and did not do: sections 7 and 10 are billing and Strassen-depth engineering; section 9 is the
retained-transport saturation certificate made exact, with its full-drop limit explained by where a saturated
neuron's mean lives; section 11 is a schedule correction. None of it changes what the chain computes beyond rounding,
except the two measured approximations (the alpha <= -2.5 drop, neutral; the L-2 join, -0.43%).
