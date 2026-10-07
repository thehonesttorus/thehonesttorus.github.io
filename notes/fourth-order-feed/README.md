# The missing fourth-order feed, by index class

Working note XXI. A second synthesis was pasted in response to note XIX. It accepts the two-channel ledger as
progress, makes two corrections to its language, proposes a change of coordinates (common gain and channel
disagreement), and changes the proposed next step: before accumulating the three dropped fourth-cumulant slices,
measure the contracted quantity through which they affect g_4, decompose it by index-multiplicity class without
building a fourth-order tensor, and derive its response to the incoming third cumulant with a scalar
first-variation formula. The synthesis did not have access to note XIX's code or outputs, so part of what it
proposes was already in the ledger; this note says which part, does the parts that were not, and records the
corrections.

What is new here:

- **The two corrections are accepted and note XIX is amended.** A recursion for cumulants is not a probability
  specification, so "the chain is a g-measure" is an analogy in moment coordinates with one quantitative anchor
  (the transport of the scale mode), not an identification; the wording in note XIX is changed. Unit row sums of
  the gain transport follow from positive homogeneity only for the scale-mixture perturbation family with the
  least-squares projection onto the two gain coordinates, which is the family note XIX used; the raw and the
  normalised transport are both kept, and the normalisation is not an independent validation.
- **In common/disagreement coordinates the missing information is channel-specific.** The one-step residual of the
  common gain is within 3% of g_3 at every layer on both networks; the residual forcing of the disagreement
  d = g_4 - g_3 is zero through layer 4 and grows to 0.004-0.005 per layer at layers 12-15, where it is the whole
  of what drives d from -0.004 to zero against a negative source forcing and a retention of 0.9. The missing feed
  is a positive input to g_4 alone, of a fifth to a quarter of g_4 per layer at depth, absent at shallow depth.
- **The contracted missing quantity was already measured** (note XIX's "true-h" column against the measured g_4 is
  Q_total minus Q_retained); its decomposition by index class is the new measurement, done in one Monte Carlo pass
  of 4e5 inputs with power sums and no tensor (`code/k4classes.py`): section 3.

- **The missing feed is one index class, and it is the scale mixture's own term on the common-mode variance.** The
  (3+1) and (1+1+1+1) classes are zero within noise at every depth; the (2+1+1) class is the whole gap (5%, 7%,
  15% of g_4 at targets 8, 11, 15), and it equals 6 g sigma_diag^2 sigma_off^2, the mixture's fourth cumulant on the
  cross term between the diagonal and the off-diagonal (common-mode) parts of the next layer's variance, to
  10-25%. The pair-slice ledger drops it because pairs see only sigma_diag^2. It needs no accumulator: it is a
  closed-form function of objects the chain holds at O(n^2). Section 5 separates its transported and generated
  parts and closes the fourth-cumulant channel of the two-gain recursion with it.

## 1. The two corrections

**A cumulant recursion is not a g-function.** van Enter, Le Ny and Paccaut's dichotomy is between conditional
probability specifications in one direction and in two; a g-function is a normalised next-symbol law given a
one-sided history, with a variation sequence measuring its dependence on the remote past. The chain carries means,
covariances, cumulant sources and a regenerated fourth-cumulant sector, and reads the next gate law from them; that
is one-sided and bounded-memory, but it is not a probability specification on a symbolic alphabet with admissible
histories, and the decay rates measured for selected cumulant responses are not the variation sequence of a
conditional law. Note XIX's sentence "the chain is a g-measure, not a Gibbs measure, by construction" is replaced
by "the chain is one-sided by construction: the analogue of a g-measure, not of a Gibbs measure", and the
quantitative content that survives is the one that was measured: the response of the output to content born at
layer b decays geometrically on every mode but the scale, whose transport is 1. The reversal the synthesis adds is
right and useful: for fixed weights the full hidden vector is an exact non-autonomous Markov state, and memory
appears only when it is replaced by retained information; the chain's error is the memory term of that projection.

**What homogeneity does and does not give.** relu(c x) = c relu(x) preserves a deterministic radial change. The
normalised gains g_3 = kappa_3 / (1.5 mu sigma^2) and g_4 = kappa_4 / (3 sigma^4) are invariant under that change,
so homogeneity alone says nothing about how a perturbation of g_3 transfers into g_4. The unit-row-sum statement of
note XIX is for a specified family: y = G y~ with G independent of y~, E G^2 = 1, Var G^2 = g; for that family
relu(G y~) = G relu(y~) makes the next pre-activation G z~, whose third and fourth cumulants are 1.5 g mu sigma^2 and
3 g sigma^4 with the same g to first order, so the family's image under the layer is the family with the same gain
plus the fresh generation, and (1, 1) is an eigenvector of the exact transport with eigenvalue 1. The projection
onto the gain coordinates (least squares against the patterns with the true marginal moments) adds second-order
terms, and the pair-slice restriction with the fourth-order truncation of the integration by parts loses 7% and
17% of the rows at depth. Note XIX reported both the raw matrix and the normalised one and used the normalised one
as a model; this note keeps that discipline and reports the one-step residuals against the raw matrix, so that the
normalisation is never used to generate the residual it is then credited with explaining.

## 2. Common gain and disagreement

With a = r34 / (r33 + r34) and b = r43 / (r43 + r44) from the row-normalised transport, s the one-loop generation
and e the one-step residual of the two-channel model against the measured gains of the next layer (computed from
the measured gains of layer l, so there is no drift), exactly

    g3_{l+1} = g3_l + a d_l + s3_l + e3_l,      d_{l+1} = (1 - a - b) d_l + (s4_l - s3_l) + (e4_l - e3_l),

with d = g_4 - g_3. The ledger of note XIX in these coordinates (`code/gac3d.py`, `outputs/gac3d_off{0,1}.txt`;
network 0 shown, network 1 within 0.001 at every entry):

| l+1 | g_3, d measured | retention 1 - a - b | source s_4 - s_3 | residual e_3 | residual e_4 - e_3 | e_3 / g_3 | (e_4 - e_3) / abs(d) |
|---|---|---|---|---|---|---|---|
| 2 | 0.0095, -0.0039 | 0.46 | -0.0026 | +0.00005 | +0.00003 | +0.5% | 1% |
| 4 | 0.0127, -0.0037 | 0.72 | -0.0016 | -0.00029 | +0.00022 | -2.2% | 6% |
| 6 | 0.0156, -0.0040 | 0.78 | -0.0023 | -0.00010 | +0.00043 | -0.6% | 11% |
| 8 | 0.0176, -0.0030 | 0.89 | -0.0010 | -0.00037 | +0.00197 | -2.1% | 66% |
| 10 | 0.0197, -0.0031 | 0.88 | -0.0020 | +0.00007 | +0.00174 | +0.4% | 56% |
| 12 | 0.0215, -0.0019 | 0.92 | -0.0017 | +0.00034 | +0.00356 | +1.6% | 185% |
| 14 | 0.0223, -0.0009 | 0.91 | -0.0021 | +0.00057 | +0.00356 | +2.5% | 420% |
| 15 | 0.0223, +0.0005 | 0.97 | -0.0014 | +0.00019 | +0.00471 | +0.9% | - |

- **The common channel is closed.** e_3 is within 3% of g_3 at every layer on both networks, with no trend: what the
  pair-slice, first-order programs retain determines the next third-cumulant gain, and the 7% row-sum deficit of
  the raw transport is compensated within the model (the one loop on the marginal state re-counts it, note XIX
  section 2.3), not by a missing feed.
- **The disagreement channel is not.** The retained model predicts that d stays at -0.004 to -0.005 (the fresh
  generation keeps forcing it negative at -0.002 per layer and the retention is 0.8-0.97), and the measurement has
  it rising to zero from layer 8 on. The residual forcing e_4 - e_3 is +0.0002 through layer 4, +0.0005 at layers
  5-6, +0.001 at 7, +0.002 at 8-11 and +0.0036 to +0.0047 at 12-15. On network 1 the same column reads +0.0005,
  +0.0012, +0.0014, +0.0019, +0.0016, +0.0028, +0.0028, +0.0031, +0.0047, +0.0048 at layers 4-15. It enters almost
  entirely through e_4: the missing information feeds the fourth-cumulant gain and leaves the third alone.
- **Its size.** A fifth to a quarter of g_4 per layer at depth. That is larger than the 14% by which the pair-slice
  restriction of the true fourth cumulant of h fell short of the measured g_4 in note XIX, because the residual here
  is against the mixture-shaped slices of unit gain times the measured gains, and the true slices of y carry
  non-mixture content in the (3,1) and (2,2) patterns as well; both numbers say the dropped content is positive and
  of order 0.003-0.005 per layer in g_4 units at depth.

## 3. The fourth cumulant of the next pre-activation by index class

The synthesis is right that the contracted quantity should be measured before its slices are accumulated, and
right that it needs no tensor. For one sample and one output neuron i, with a_j = W_ij X_j and X the centred
post-activation, the power sums p_r = sum_j a_j^r give the five index-multiplicity classes of (sum_j a_j)^4
exactly: (4) p4; (3+1) 4 (p1 p3 - p4); (2+2) 3 (p2^2 - p4); (2+1+1) 6 (p1^2 p2 - 2 p1 p3 - p2^2 + 2 p4);
(1+1+1+1) the rest. Over a batch they are four matrix products X^{o r} (W^{o r})^T. The Gaussian pairings of each
class are subtracted with the sample covariance of h by inclusion-exclusion over coincident indices, so the five
class cumulants sum to kappa_4(z_i) = E z^4 - 3 (E z^2)^2 exactly, and each is projected on 3 sigma_z^4 to give
its share of g_4. Classes (4) and (2+2) are what the pair-slice ledger retains; (3+1), (2+1+1) and (1+1+1+1) are
the dropped slices. `code/k4classes.py` does this in one pass of 4e5 inputs through the network at source layers
7, 10 and 14 (target pre-activations 8, 11 and 15), in double precision; the plug-in cumulant's finite-sample bias
is of order sigma^4 / T and is ignored at this T.

Network 0 (`outputs/k4classes_off0.txt`; the ledger's true-h column of note XIX is the pair column's independent
check: 0.01391, 0.01639, 0.01960 there against 0.01392, 0.01645, 0.01961 here):

| target z_(l+1) | g_4 total | (4) | (3+1) | (2+2) | (2+1+1) | (1+1+1+1) | pair = (4) + (2+2) | pair / total |
|---|---|---|---|---|---|---|---|---|
| 8 | 0.01470 | 0.00032 | 0.00002 | 0.01359 | 0.00081 | -0.00005 | 0.01392 | 0.947 |
| 11 | 0.01766 | 0.00014 | -0.00001 | 0.01631 | 0.00154 | -0.00032 | 0.01645 | 0.932 |
| 15 | 0.02298 | 0.00014 | 0.00002 | 0.01947 | 0.00307 | 0.00027 | 0.01961 | 0.853 |

| target | measured (2+1+1) | mixture 6 g sigma_diag^2 sigma_off^2 at g = g_4 total | measured (1+1+1+1) | mixture 3 g sigma_off^4 |
|---|---|---|---|---|
| 8 | 0.00081 | 0.00076 | -0.00005 | 0.00020 |
| 11 | 0.00154 | 0.00131 | -0.00032 | 0.00041 |
| 15 | 0.00307 | 0.00236 | 0.00027 | 0.00094 |

- **The missing feed is one class, and it is the mixture's own.** The (3+1) class is zero at every depth and the
  (1+1+1+1) class is within noise of zero; the whole of the gap between the pair restriction and the measured g_4
  (5%, 7%, 15% at targets 8, 11, 15) is the (2+1+1) class. For a scale mixture of gain g the full fourth cumulant
  of z_i is 3 g sigma_i^4 with the full variance sigma_i^2 = sigma_diag^2 + sigma_off^2, where sigma_diag^2 =
  sum_j W_ij^2 Var h_j is what the pair slices see and sigma_off^2 = sum_{j != k} W_ij W_ik Cov(h_j, h_k) is the
  common-mode part; expanding, the (2+2) and (4) classes carry 3 g sigma_diag^4, the (2+1+1) class 6 g sigma_diag^2
  sigma_off^2 and the (1+1+1+1) class 3 g sigma_off^4. Evaluated on the measured variances with the measured total
  gain, the mixture's (2+1+1) term is 0.00076, 0.00131, 0.00236 against the measured 0.00081, 0.00154, 0.00307:
  the class is the mixture acting on the off-diagonal part of the next layer's variance, to 10-25%, with the
  remainder of the right sign for the part generated at the layer rather than transported (section 5).
- **Network 1 repeats it** (`outputs/k4classes_off1.txt`): (2+1+1) = 0.00065, 0.00132, 0.00363 at targets 8, 11, 15
  against the mixture's 0.00076, 0.00133, 0.00311; (3+1) zero; (1+1+1+1) within 0.0003 of zero; pair / total
  0.947, 0.905, 0.826.
- **Why it was invisible to the mean off-diagonal share.** The mean over neurons of sigma_off^2 / sigma^2 is
  negative (-1% to -4%: most post-activations are slightly anticorrelated), which is why note XIX dismissed the
  off-diagonal variance as the cause. The projection on 3 sigma^4 weights the large-variance neurons, which sit
  on the mean direction and have a positive common-mode share; on that weighting the share is +6% at target 15,
  and twice that is the 13% the (2+1+1) class carries.
- **So the dropped slices do not need accumulating.** The feed is a closed-form function of objects the chain
  already has at O(n^2): sigma_diag^2 = (W o W) var_h and sigma_off^2 = sigma^2 - sigma_diag^2. This is what the
  synthesis asked for, in the strong form: not the contracted missing quantity but its formula.

## 4. Which slices of the source law feed the dropped classes

The pair-slice ledger is closed under pairs, and that is why it closes the common channel and not the disagreement.
A pair moment E[h_a^p h_c^q] of the post-activation depends only on the pre-activation variables y_a, y_c, so its
first variation (1/6) sum T_ijk E[d_ijk (h_a^p h_c^q)] involves only the entries of kappa_3(y) with i, j, k in
{a, c}: the diagonal and the (2,1) slice. Likewise the classes (4) and (2+2) of kappa_4(z) involve only pair cumulants
of h. So the retained programs form a closed system on pair slices, exact at first order (note XIX: the sum of the
three first-order pieces equals the true-h pair restriction to four decimals), and the dropped classes are fed by
what pairs cannot see:

- the (2+1+1) and (1+1+1+1) classes of kappa_4(z) are contractions of the triple and quadruple cumulants of h over
  distinct indices;
- those are fed, at first order in the source law, by the all-distinct (1,1,1) slice of kappa_3(y) through products
  of three gate derivatives, T_jkl Cov(h_j' h_k' h_l', h_m), and by the (2,1) slice through the second derivative of
  ReLU, which is a delta function at the gate: T_jjk E[delta(y_j) h_k' h_l h_m]_connected = (phi_j / sigma_j) times a
  conditional moment at y_j = 0. The latter are the "second-order gate terms" that note XVIII section 5 listed as the
  remaining class of the chain's error, and they are the terms in which the gate is not a product;
- the (3+1) class is fed by the (2,1) and (3,1) slices of y with one spectator index, and is incoherent in W (odd
  powers) unless the slice has a coherent structure.

The chain carries the all-distinct slice of kappa_3(y): the sources are CP tensors T_abc = sum_j w2_j Sym(A_aj A_bj
P_cj) + ..., and their all-distinct entries are transported with the pair entries. What the chain does not do is
let them feed the fourth cumulant of the next layer through anything but the product gate. The synthesis's scalar
first-variation formula, delta kappa_4(Y) = L[(Y - m)^4 - 6 v (Y - m)^2 - 4 c_3 (Y - m)], evaluated for L the
Gaussian score of a CP source (H_T(y) = (1/6)[T(zeta, zeta, zeta) - 3 sum_ab T_abc (S^-1)_ab zeta_c], zeta = S^-1
(y - mu)), is the route to that feed without a tensor: a Gaussian sample of the carried state, the score polynomial
at O(n rank) per source, and the fourth-order readout of z_i for all i at O(n^2) per sample. Its noise is the
obstacle: the readout's rms is about ten times sigma^4 and the feed is of order 0.003 x 3 sigma^4, so one output
neuron needs 1e6-1e7 samples and the coherent fit over 1024 neurons about 1e5, at a few 1e6 flops per sample, which
is minutes, not hours. It is the next computation if the class decomposition puts the missing feed in the
all-distinct classes; if it puts it in the (3+1) class, the feed comes from the (3,1) slice of y that the chain's
lambda C_off ansatz under-fits (note XVIII section 2) and the route is the one already open there.

## 5. Generation against transport, and the closed channel

The one-loop generation of z_(l+1) was decomposed by the same five classes at every source layer, with y_l Gaussian at
the true marginal state (`code/k4classes_gauss.py`, `outputs/k4classes_gauss_off0.txt`, T = 3e5), so that the dropped
classes of the true fourth cumulant can be split into what the layer generates and what it transports. At the three
targets where the true classes were measured (`outputs/k4classes_off0.txt`), with the pair ledger of note XIX at the
same layers (one loop + variation under the true third-cumulant slices + variation under the true fourth-cumulant
slices, all pair-restricted) and the mixture's closed-form transport of the dropped classes at the measured gain of
the source layer:

| target | g_4 measured | pair classes, true | pair ledger | dropped classes, true | dropped, one loop | mixture transport of dropped | remainder |
|---|---|---|---|---|---|---|---|
| 8 | 0.01457 | 0.01392 | 0.01387 | 0.00078 | -0.00017 | 0.00083 | +0.00012 |
| 11 | 0.01770 | 0.01645 | 0.01638 | 0.00121 | -0.00012 | 0.00158 | -0.00025 |
| 15 | 0.02278 | 0.01961 | 0.01966 | 0.00337 | -0.00013 | 0.00309 | +0.00041 |

- **The fourth-cumulant channel is closed to 1%.** Pair ledger plus mixture transport of the dropped classes plus
  the one loop's dropped classes gives 0.01453, 0.01784, 0.02262 against the measured 0.01457, 0.01770, 0.02278: 0.3%,
  0.8%, 0.7%. The remainder (+0.0001, -0.0003, +0.0004) is at the noise of the (1+1+1+1) class. The one loop
  generates nothing in the dropped classes at any depth (-0.0001 to -0.0002); what the pair restriction misses is
  transport, and it is the mixture's own, in closed form.
- **The feed the synthesis asked for is therefore derived.** For the gain carried as a scale mixture, the dropped
  classes of the next fourth cumulant are 6 g sigma_diag^2 sigma_off^2 + 3 g sigma_off^4 per output neuron. Nothing
  needs accumulating and no tensor is touched: the two variances are O(n^2) from the carried state.

**Why the two-gain recursion still does not close.** `code/gac3c.py` (`outputs/gac3c_off0.txt`) reruns the closed
recursion of note XIX with the fourth row completed: the mixture's dropped-class transport added to r44 (0.82 to
0.93-0.96 at depth) and the one loop taken with all five classes.

| l | measured g_3, g_4, ratio | A: note XIX, pair | B: pair + dropped-class transport, full one loop | C: rows normalised, full one loop |
|---|---|---|---|---|
| 4 | 0.0127 0.0090 0.71 | 0.0132 0.0091 0.69 | 0.0133 0.0094 0.71 | 0.0138 0.0090 0.65 |
| 7 | 0.0170 0.0128 0.75 | 0.0176 0.0117 0.66 | 0.0177 0.0128 0.72 | 0.0205 0.0127 0.62 |
| 10 | 0.0197 0.0166 0.84 | 0.0203 0.0117 0.58 | 0.0205 0.0147 0.72 | 0.0263 0.0159 0.61 |
| 13 | 0.0220 0.0207 0.94 | 0.0210 0.0091 0.44 | 0.0214 0.0146 0.68 | 0.0316 0.0178 0.57 |
| 15 | 0.0223 0.0228 1.02 | 0.0205 0.0069 0.34 | 0.0210 0.0141 0.67 | 0.0342 0.0184 0.54 |

Variant B tracks the measured g_4 to 2% through layer 7 and then falls behind (0.0141 against 0.0228 at layer 15),
even though its per-step transport is now the mixture's. The ledger says where the rest is: the column "d3 true
against r43 g_3" of the attribution is 0.00181 against 0.00090 at target 8, 0.00251 against 0.00143 at 11 and
0.00229 against 0.00036 at 15. The true third-cumulant pair slices of y feed the next fourth cumulant two to six
times more than the mixture-shaped slices at the fitted gain do, while they feed the next third cumulant exactly as
the mixture does (d3 = 0.0204 against r33 g_3 = 0.0206 at layer 14). The third-cumulant content of the law is not
one number: note XIX's record already showed the two patterns of its (2,1) slice, mu_a S_ac and sigma_a^2 mu_c,
carrying different amplitudes at depth (0.0164 and 0.0234 at layer 15, equal for a mixture), and it is the
non-mixture part that drives the fourth cumulant. This is the synthesis's question answered: the additional
information that changes the fourth-order feed while leaving the third-order behaviour intact is the shape of the
third-cumulant slices beyond their mixture projection, and a two-gain state cannot carry it. The chain can and
does: its sources are the full third-cumulant tensor.

**What this says about the chain's fourth-cumulant sector.** The chain regenerates the sector from a transported
diagonal and a fitted per-layer table, with no feed from the sources. The ledger now gives the three pieces the
table stands in for at depth, in g_4 units per layer: generation about zero (-0.0001 to -0.0007), the mixture's
dropped-class transport +0.001 to +0.002, and the sources' feed into the fourth cumulant through the pair programs
+0.002 to +0.003 (the "d3" column, two to six times what a mixture-shaped feed would be). Both of the latter are
computable inside the chain at O(n^2) per layer: the first from (W o W) var_h and sigma^2, the second from the
(2,2)-pattern term program that the chain already has for pair moments (`TERM_SPECS` column (2,2)) applied to the
carried D3 and D21 slices. That is a derived replacement for the fitted table with the physics that K4SM = 3 and 4
lacked (both tied the sector to the mixture shape, which is the part the table was not missing), and it is the one
chain experiment this note proposes. Its expected size is the table's measured misfit, 15% of the sector at layer 14
(note XVIII section 2), and the sector's leverage on the error was measured there at the few-percent level, so the
experiment is worth one round and not more. The network 1 one-loop decomposition is queued and will be appended to
`outputs/` when it completes.

## 6. What transfers

| from the synthesis | what it is here | number or decision |
|---|---|---|
| a cumulant recursion is not a g-function | accepted; note XIX amended | wording changed; the measured transport of the scale survives |
| homogeneity gives unit row sums only for a specified family and projection | accepted; the family and projection are the ones used; raw and normalised both kept | 7% and 17% raw deficits at depth, now attributed |
| common / disagreement coordinates | the common channel is closed, the disagreement forcing is channel-specific and grows with depth | e_3 within 3%; e_4 - e_3 = 0.004-0.005 per layer at 12-15 |
| measure the missing contracted quantity first | it was the ledger's true-h column; its class decomposition is new | (2+1+1) is the whole gap; (3+1), (1+1+1+1) zero |
| the index-multiplicity classes without a tensor | done, one Monte Carlo pass, power sums | pair classes match the ledger to three digits |
| derive the source response of the missing feed | the mixture's closed form for the dropped classes, and the measured d3 feed in the pair classes | channel closed to 1% at targets 8, 11, 15 |
| a two-channel model with a memory operator | the missing information is the non-mixture shape of the third-cumulant slices, which the chain carries | d3 feed 2-6x the mixture's |
| keep the chain and the table frozen | agreed for this round; one derived replacement proposed | expected size: the 15% misfit of the sector at layer 14 |

