# The one-step fourth-cumulant error: representation, closure and the omitted classes

Working note XXXVI. It runs the discriminating computations that notes XXXII-XXXV set up:
- the checkpoint's cell comparison in the chain's own chart (note XXXV section 2);
- the validation of the adjoint error functional of note XXXI against the measured oracles;
- the split of the one-step fourth-cumulant error into the pair program's closure error and the transport of the
  post-activation classes the chain omits.

Monte Carlo truth throughout (1.6e7 inputs, two independent halves). The production configuration of note XXVIII
(raw 2.2525e-8 and 2.1266e-8 on networks 0 and 1).

## 1. The cells: representation repairs do not move the one-step error

`code/cells.py` takes the one-step dump of note XXXII: the chain fed true pre-activation statistics at layer l. It
transports the y-level arrays the chain computed (K4v, K22) to the three slices of z_(l+1) in several declared tensors,
each slice taken from one tensor, and compares them with the Monte Carlo slices of z_(l+1) from the other half of the
sample. The cells (note XXXV section 2):
- **chain:** the chain's own output.
- **c00:** seed chart with the residual core at 2I.
- **c10:** residual core at the actual metric W W^T.
- **c01:** exact transport of the literal retained tensor D(K4v) + B(K22), which is cell 11 in this chart.
- **c6:** c01 plus the transported y (3,1) class from note XXXIV's closed form.
- **c6-noλ:** c6 without the lam core.

Every cell except c6-noλ carries the chain's own lam core. The replication of the chain's slices from the dumped
arrays is exact to 1e-7.

Relative error of the diagonal | (2,2) | (3,1) slices against truth (`outputs/cells_mc2_off{0,1}.txt`):

| net 1, layer | MC noise | chain | c00 | c10 | c01 | c6 | c6-noλ |
|---|---|---|---|---|---|---|---|
| 2 -> 3 | .089 .106 .524 | .142 .162 .838 | .142 .159 .845 | .142 .159 .846 | .142 .159 .843 | .136 .156 .827 | .160 .168 .877 |
| 6 -> 7 | .047 .056 .226 | .161 .178 .744 | .161 .175 .780 | .161 .175 .781 | .161 .175 .780 | .158 .173 .769 | .233 .213 .914 |
| 10 -> 11 | .036 .043 .124 | .250 .246 .750 | .250 .244 .769 | .250 .244 .769 | .250 .244 .769 | .245 .241 .760 | .343 .297 .952 |
| 13 -> 14 | .031 .037 .088 | .292 .299 .736 | .292 .298 .749 | .292 .298 .749 | .292 .298 .749 | .287 .295 .740 | .413 .371 .966 |

Network 0 is the same to the second digit: chain diagonal 0.137-0.292, c6 0.134-0.287, c6-noλ 0.162-0.426.

- **The representation is not the defect.** The coherent seed chart, the actual metric and the literal transport
  differ from the chain by at most 0.003 on the diagonal and (2,2) slices. They are slightly worse on the (3,1) slice,
  where the chain's lam-only slice happens to fit better. THEORY_3's ledger terms (D_metric, D_known, D_unresolved
  for R_W(E)) are each 0.2-1.6% of the truth and uncorrelated with the residual (|corr| <= 0.07). The (2,2) cross term
  2 c_ab^T K22 c_ab, omitted from every full-matrix slice, is 0.6-1.8% of the main term.
- **The closed-form (3,1)_y class helps slightly and consistently.** c6 is 0.003-0.006 better than c01 in every
  slice, at every layer, on both networks. Its transported size is 1.4-2.7% of the diagonal truth.
- **The lam core carries content that nothing else in the state does, and its share grows with depth.** Removing it
  costs 0.02 on the diagonal at layer 2 and 0.12-0.15 at layers 12-13. The fitted lam itself is flat
  (0.0085-0.0105). This is the depth profile of an accumulating mode, not of a per-layer defect.
- **What remains is not representation.** The one-step error lies in the pair program's y-slice closure, or in the
  y (2,1,1) and (1,1,1,1) classes that the lam core stands in for. Section 3 splits the two.

## 2. The adjoint error functional predicts the oracles

Note XXXI's `errbudget.py` (now `code/errbudget.py`) takes the per-layer errors of the three diagonals that the mean
map reads (variance, kappa3, kappa4), as dumped by the chain. It weights them by the mean map's first-order
coefficients and propagates them to the output with the exact Jacobian at Monte Carlo truth. That predicts the change
of the final MSE when each channel is replaced by truth at every layer. Its first runs lost their output to a
filter; rerun:

| channel replaced at every layer | predicted, net 0 (mc2 / mc3) | measured, net 0 | predicted, net 1 | measured, net 1 |
|---|---|---|---|---|
| variance diagonal | -51.7% / -52.8% | -52.0% / -51.9% | -38.1% | -36.8% |
| kappa3 diagonal (D3) | -16.2% / -15.2% | -14.8% | -14.4% | -15.0% |
| kappa4 diagonal (g4) | -27.7% / -28.5% | -30.1% | -13.5% | -15.5% |
| all three | -95.1% / -104.3% | | -88.6% | |

Measured values are the noise-free oracle extrapolations of note XXXI (`notes/birth-address/outputs/oracle_mc*`).

- **The first-order functional is quantitatively right.** Every prediction is within 2 points of its oracle, so it can
  score a proposed repair from a one-step dump, without a full oracle run.
- **To first order the three diagonals carry the whole error.** What the three channels leave unexplained is 11%
  (net 1) and -4% to +5% (net 0).
- **The per-layer weights rise steadily with depth.** The last two layers carry 0.10-0.13 of |e|^2 each in the
  variance channel and 0.04-0.07 in the others.
- **The variance channel is the largest, and it is downstream.** Note XXXIII's F5 showed that the variance errors are
  produced by the kappa3 readouts and the kappa4 sector through the pair program.

## 3. The class split on the post-activation level

Network 1, Monte Carlo mc4 (1.6e7 inputs, seed 402), one-step dump inputs from mc2 (`outputs/yclasses_mc4_off1.txt`).
Slices are diagonal | (2,2) on 12,000 random pairs | (3,1) off-diagonal, each relative to the truth's norm. The chain
columns give the diagonal and (2,2) slices. "lam -> fitted D" and "lam -> g_conn D" are the chain with its lam core
replaced by candidate D at its fitted coefficient and at the trace-sector coefficient:

| l -> l+1 | abs R / t | closure (diag) | repr (diag) | chain | lam -> fitted D | lam -> g_conn D | D explains | D coef x 2/n | g_conn | R3 / t3 (D3, D21) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 -> 2 | .144 .148 .835 | .005 | .038 | .137 .153 | .128 .150 | .127 .150 | .29 .12 .16 | .0040 | .0056 | .35 .49 |
| 3 -> 4 | .161 .157 .866 | .005 | .027 | .131 .146 | .124 .141 | .125 .141 | .45 .27 .35 | .0063 | .0065 | .31 .49 |
| 5 -> 6 | .205 .184 .894 | .006 | .021 | .154 .161 | .140 .150 | .142 .152 | .56 .37 .46 | .0076 | .0063 | .28 .45 |
| 7 -> 8 | .244 .221 .922 | .009 | .018 | .171 .187 | .151 .169 | .165 .178 | .63 .45 .52 | .0082 | .0054 | .28 .45 |
| 9 -> 10 | .284 .262 .939 | .012 | .015 | .206 .220 | .168 .195 | .192 .210 | .66 .48 .55 | .0089 | .0051 | .28 .46 |
| 11 -> 12 | .383 .327 .957 | .011 | .015 | .269 .264 | .190 .212 | .255 .251 | .76 .60 .64 | .0091 | .0044 | .28 .44 |
| 13 -> 14 | .409 .362 .964 | .017 | .014 | .288 .292 | .195 .225 | .281 .281 | .78 .63 .69 | .0092 | .0039 | .26 .42 |
| 14 -> 15 | .479 .405 .966 | | | | | | .81 .68 .72 | .0092 | .0036 | .26 .43 |

R is reproducible. The split-half correlation is 0.83-0.99 from layer 3 on, and its noise share is below 10% from
layer 4 on (77% at layer 0, where R is small).

- **The pair program is essentially exact; its closure is not the defect.**
  - From true z inputs, its y slices are K4v to 0.5-1.3% (corr 1.000) and K22 to 6-35% (corr 0.93-0.998, worst at
    depth). The (3,1)_y closed form of note XXXIV is good to 11-13% (corr 0.99).
  - Transported, all of it is 0.5-1.7% of the next diagonal and (2,2) slice.
  - The cosine between the closure error and the lam core's misfit is -0.04 to +0.12, so lam is not compensating the
    closure.
  - The representation residual (Ritz approximation, missing (3,1)_y and cross term) is 1.4-3.8% on the diagonal.
  - The (3,1) slice's 0.2-0.5 representation residual is the chain's lam-only (3,1) slice, which lacks the
    transported pair part.
- **The one-step error is the omitted classes net of the lam core.** "Classes - lam" is within 0.01 of the chain's
  whole error on the diagonal at every layer.
- **The omitted classes are large, and are most of the (3,1) slice.**
  - |R|/|t| on the diagonal rises from 0.14 to 0.48 with depth. On the (2,2) slice it is 0.15 to 0.41.
  - On the (3,1) slice it is 0.80-0.97. The (3,1) slice of z' is almost entirely the transport of the y (2,1,1) and
    (1,1,1,1) classes.
  - Note XXI saw only 5-15% of the omitted classes because it measured the projection on 3 sigma^4. Most of R is
    per-neuron structure that averages out in that projection.
- **P1 holds: one radial-tangent mode is the leading part of R, coherently in all three slices.**
  - Candidate D, with one coefficient per layer, explains 29/12/16% of R at layer 1 and 81/68/72% at layer 14
    (diagonal, (2,2), (3,1)).
  - The coefficient fitted per slice agrees across the three slices to 1-5%.
  - D beats the chain's own lam-core shape (A: 50/33/32% on average over layers 2-14) and the centred scale mode
    (B: 54/36/42%). The mean terms matter.
  - At layer 0 -> 1, where z_0 is exactly Gaussian, D's fitted coefficient is 0.2 (explaining 0.4%), not 1. The input
    radius contributes nothing net, as section 4 predicted. The coefficient at depth (2.7-4.7) is generated scale.
- **P2 holds at shallow depth only, and P3 fails.**
  - The fitted coefficient, in Var(s^2) units, rises from 0.004 (layer 1) to 0.0092 and saturates from layer 10 on.
  - g_conn, the connected Var|y|^2, matches it at layers 2-5 (0.0063 against 0.0063 at layer 3). It then falls to
    0.0036.
  - So the soft mode that carries R at depth is not the isotropic norm of y. Its amplitude rises while the norm's
    falls: the scale that couples to the off-diagonal (common-mode) variance is a different collective variable from
    |y|^2.
- **The Gaussian second-order class (C) is minor.** It explains 10-14% of R's diagonal at layers 1-5 with coefficient
  1.5-1.7, and nothing at depth. This confirms note XXI's zero one-loop generation.
- **The prize.** Replacing the lam core by D at its fitted coefficient lowers the chain's one-step diagonal error at
  depth from 0.27-0.29 to 0.19-0.21, and the (2,2) error from 0.27-0.29 to 0.21-0.22. With g_conn's coefficient the
  gain is small at depth, because g_conn is too small there.
- **kappa3: the legs carry a third to a half of the readouts.** R3, the transported all-distinct class of y, is
  26-35% of the D3 slice and 42-49% of the D21 slice at every depth. The scale mode explains 2-39% of it, rising with
  depth. The rest is the legs' own all-distinct content, which section 3c tests against the first-order gated
  transport.

`code/mcstats.py ... post` accumulates the same slices for y_l = relu(z_l) on the same samples as z_(l+1). On any one
sample the empirical cumulant tensor of z_(l+1) = W y_l is exactly W# of the empirical cumulant tensor of y_l. So

    t = T_pair(y) + R,   T_pair(y) = W#[D(k4_y) + S22(K22_y) + S31(K31_y)],

with the (2,2) cross term included, holds exactly per sample. R is exactly the transport of the y (2,1,1) and (1,1,1,1)
classes: the omitted classes, measured, with none of the pair-slice sampling noise. The identity and the transport were
checked against the brute-force empirical fourth-cumulant tensor on a width-7 network (agreement 1e-8).
`code/yclasses.py` reports per transition:
- |R|/|t| and its split-half reproducibility;
- the closed-form candidates for R below;
- the pair program's closure error, transported;
- the exact split of the chain's one-step error into closure, omitted classes net of the lam core, and representation.

**Replicate.** An independent Monte Carlo of network 1 (mc5, seed 403, no chain dump; `outputs/yclasses_mc5_off1.txt`)
reproduces every number to the third digit. Averaged over layers 2-14:
- D explains 62.3/45.9/52.5%, against 63.0/46.1/52.7%.
- The fitted coefficients are 2.69 -> 4.70, against 2.72 -> 4.72.
- g_conn is identical.

These are properties of the network, not of the sample.

**The chain's post-activation arrays are literal cumulants** (`code/ysem.py`, `outputs/ysem_off1.txt`). In the one-step dump
(true pre-activation inputs), against the Monte Carlo of y:
- the mean pk1v, the variance K2v and the kappa3 diagonal K3v match to <= 0.4%;
- K11 is the full off-diagonal covariance, matched to <= 1.4%;
- K21[a, b] = kappa(y_a, y_a, y_b), matched to 2%.

Every ingredient of G D is therefore in the chain's state at the end of each layer.

### 3b. The same split for the third cumulant (prediction written before the data)

The kappa3 readouts carry about 40% of the error (note XXXI). The same identity holds for kappa3:

    (D3, D21 slices of z') = T3_pair(y) + R3,
    T3_pair = W#[diagonal and (2,1) classes of kappa3(y)].

R3 is exactly the transport of the all-distinct class kappa(y_a, y_b, y_c). The pair part is local: the pair program
computes it from pair statistics of z. R3 is what only the legs carry. Checked against the brute-force tensor on the
width-7 network (`code/check_k3_split.py`, 1e-6). By note XXXIV's support preservation, R3 has three sources at first
and second order:
- **The legs' own all-distinct z entries**, gated: Phi_a Phi_b Phi_c T_abc.
- **The Gaussian second-order term**, one facet and two covariances:
  rho_a Phi_b Phi_c C_ab C_ac + rho_b Phi_a Phi_c C_ab C_bc + rho_c Phi_a Phi_b C_ac C_bc, with rho = phi/s.
- **The scale mode's kappa3 tangent**, Cov(s, s^2) Sym3(mu, Sigma) with Cov(s, s^2) ~ Var(s^2)/2. On all-distinct
  indices it is mu_a C_bc + 2 perms, first order in C, and carries the large ReLU means.

`yclasses.py` reports |R3|/|t3| and the scale-mode fit on R3. If the scale-mode term explains R3 with
Var(s^2) ~ g_conn, the legs' error is mostly the scale mode's. If not, it is in the legs' genuine all-distinct
content or in the facet term.

### 3c. The third cumulant's transport law is exact; the chain's error is the legs' representation

`code/mclegs.py` re-runs the Monte Carlo on the same samples with the true gates Phi_l. It accumulates the contractions
E[d_a U_i^2], E[U_i^2 U_c] and E[U_i^3] of U = (W o Phi) d. Their all-distinct parts (inclusion-exclusion with the
true D21 and kappa3 diagonal) are the contractions of the true all-distinct third cumulant T of z_l that the theorem
transports at first order. They are checked against the brute-force tensor in `code/check_legs_contractions.py`. The
second-order one-facet term has a closed form, checked in `code/check_facet_transport.py`:

    kappa(y_a, y_b, y_c) |_(C^2) = rho_a Phi_b Phi_c C_ab C_ac + 2 perms,   rho = phi/s.

Results, network 1 (`outputs/mclegs_mc4_off1.txt`; the independent replicate `outputs/mclegs_mc5_off1.txt` agrees to
the third digit):
- **R3 is the first-order gated transport of the legs plus the one-facet term, both at coefficient 1.**
  - Together they explain 99.1-99.4% of R3 in the D3 slice at every layer from 1 to 14. The fitted coefficients are
    0.98-1.00 (legs) and 0.94-1.00 (facet).
  - At layer 0, where z_0 is Gaussian and T = 0, the facet term alone explains 92% at coefficient 1.000.
  - The legs alone explain 64-88%. The facet term, one gate-boundary density with two covariances (note XXXIII's
    mechanism), is a quarter of the all-distinct transport and is not optional.
  - This is note XXXIV's support-preservation theorem plus the Mehler second order, confirmed at the network's
    correlated reference with nothing fitted.
- **The leg-fed (2,1,1) class of kappa4 is real but minor.** Jointly with D, its coefficient is near its theoretical
  value of 1 (0.83-1.43). It adds 5-15 points of explained R at layers 1-5 and nothing at depth. D plus leg-fed
  explains 44-82%. The rest of R (about 20-30% at depth) is in second-order products beyond the Gaussian C x C term.

**Where the chain's kappa3 readout error is** (`code/k3chain.py`, `outputs/k3chain_off1.txt`). The one-step dump with
every oracle on records the chain's own D3 and D21 at l + 1 before replacement. Its y-level K3v and K21 give the pair
part, and the rest is the chain's legs part. Against the exact law:
- **Total error.** It is D3 3.5-4.8% and D21 5.2-7.2% of the truth slice.
- **The pair part is negligible:** 0.3-1.2%.
- **The legs part carries essentially all of it.** That is an 11-19% error relative to R3, growing with depth.
- **The legs error is orthogonal to both terms of the law.** On the facet term its coefficient is about 0 (explaining
  <= 2%); on the first-order legs term, also about 0 (<= 5%). The chain's legs carry the facet term and the right
  transport. Their error is representation: the finite-rank, compressed storage of the sources.

With note XXXI's precision law (Delta MSE / MSE ~ 180 x the readout error energy), a 4-5% rms readout error is about
30-40% of the MSE, which is the D3 + D21 oracle's measured -40%. This locates note XXXIII's conclusion, that the binding
constraint is computing the retained statistics, in one place: the fidelity of the legs' compressed representation.
Section 3e measures which compression stage it comes from.

### 3e. The legs error is not compression

Each compression stage of the legs was relaxed in turn in the all-oracle one-step dump, and the legs-part error measured
against the exact law (`outputs/k3_representation_budget.txt`; D3 | D21, mean over layers):

| variant | D3 legs error | D21 legs error |
|---|---|---|
| production | 0.0396 | 0.0605 |
| residual leg rank 4 -> 32 | 0.0388 | 0.0598 |
| D21 feedback rank 2 -> 16 | 0.0391 | 0.0596 |
| no confinement of old sources | 0.0384 | 0.0591 |
| shared basis 320 -> 640 | 0.0385 | 0.0592 |
| nested tier off | 0.0395 | 0.0604 |
| three range-finder passes | 0.0393 | 0.0603 |
| exact matmul (no Strassen) | 0.0396 | 0.0605 |
| everything relaxed at once | 0.0371 | 0.0575 |

About 94% of the legs error survives every relaxation, and Mehler order 3 in the pair programs changes nothing. The
error is in what the births contain.

**The births implement the second-order law, and the error is its remainder** (`outputs/k3chain_off1.txt`, rerun with the
law's remainder). The birth hub is B1 = Sym(X1 x P x Y1), with P = I at birth, X1 = 3 a_b and Y1 = a_b d(w2),
a_b = d(Phi) C_off. The chain's Wick weights are exactly the gate coefficients: w1 = Phi, w2 = rho = c(1,2),
w3 = c(1,3), checked numerically. So B1 is the one-facet term exactly. The chain's legs error correlates with minus
the law's own remainder, rem = R3 - legs - facet:
- the correlation is +0.71 at layer 1, falling to +0.43 at layer 13;
- the coefficient is 0.73-0.97, close to 1;
- it explains 50% of the legs error at layer 1 and 18-25% at depth.

The rest is the same remainder born at earlier layers and carried in the legs. The chain is the second-order law, and
its kappa3 error is the law's next order, born fresh at each layer and accumulated.

**The next order, in closed form** (`code/k3rem.py`, `outputs/k3rem_mc4_off1.txt`; contractions checked in
`code/check_k3_remainder_terms.py`). Second-order cross terms between two tangent entries carry
omega1 omega2 prod_v c_v(1, k1 + k2): the derivative orders add, because the cumulant expansion is the exponential of
commuting derivative operators.
- **GC1 dominates.** It is the D21 source with a covariance on its doubled index,
  (Gamma_uv/2) c_u(1,3) Phi_v Phi_w C_uw, with c(1,3) = -alpha phi/s^2 the slope of the gate density.
  - At coefficient 1 it explains 10-53% of rem (34-53% from layer 5 on), with fitted coefficient 0.73-0.94.
  - Its norm is 0.24 of rem at layer 1, rising to 0.75-0.87 at depth.
- **GC2 adds 2-9%.** It is the covariance on the single index, (Gamma_uv/2) rho_u rho_v Phi_w C_vw.
- **The third-order Gaussian terms are negligible.** These are the double edge and the triangle
  rho rho rho C C C. The layer-0 remainder is Monte Carlo noise.
- **GC1 + GC2 + DS explain 41-57% of rem at depth.** The rest is presumably the legs' own gate covariance
  T_abc rho_a rho_b C_ab, the Edgeworth correction of the gates, and Gamma x Gamma.

GC1 is note XXXIII's "third-chaos response of the D21 source", the facet mechanism, now with its exact coefficient.

### 3f. The D21 feedback carries the Gamma x C terms at half weight

The chain's V18 feedback (X1 += 1.5 d(w2) D21, Y1 += 0.5 d(w1) D21^T d(w3)) has exactly the GC1 (via Yt) and GC2 (via
Xt) structures, plus one of the four Gamma x Gamma cases. Normalised against the facet part of the same hub, per
unordered pair of leaves:
- the facet enters at 6 times its theoretical weight (rho_j);
- GC1, GC2 and the Gamma x Gamma case enter at 3 times theirs.

So every D21 feedback term is at half the theorem's weight. Regressing the chain's legs error on -GC1 gives coefficient
0.6-0.7 at every layer (`outputs/k3chain_gc_regression.txt`), and the rank of the D21 approximation does not change it.

Switch V40_FB_SY = 2 doubles the Yt weight (GC1 at the theorem's value); V40_FB_SX = 2 doubles Xt (GC2). In the
all-oracle one-step dump (`outputs/k3_feedback_weights.txt`):
- **With Yt doubled, GC1 disappears from the legs error.** Its coefficient goes from 0.66 to 0.01, -0.03 and -0.10.
- **The legs error falls.** D3 goes 0.0396 -> 0.0379 and D21 0.0605 -> 0.0566 (at layers 10-14, 0.0433 -> 0.0412 and
  0.0676 -> 0.0614). That is 8% and 12% of the error energy.
- **The other variants.** Doubling Xt alone is slightly worse; both doubled is the same as Yt alone. Feedback rank 16
  on top gives 0.0376 / 0.0561.
- **What the hub cannot do.** Being bilinear, it cannot set GC1, GC2 and Gamma x Gamma independently; the theorem
  would also need the three Gamma x Gamma cases the hub lacks.

Note XXXI's precision law (Delta MSE / MSE ~ 180 x the readout error energy) turns the measured drop into about -5% of
the MSE, at no cost. Being derived rather than fitted, it is tested on 32 networks, paired.

### 3g. The derived feedback weight in the chain, and what the two in-chain tests share

**Free-running, paired over 32 networks** (`outputs/fb_paired_32nets.txt`, production against V40_FB_SY = 2, the same
cost): the output is worse on 25 of 32, +2.3% +- 0.5% (mean of per-network ratios).

On network 1 the free-running chain behaves as the theory says in its own channel (`outputs/fb_free_off1.txt`):
- **Mid-depth readouts improve.** The kappa3 readouts are better at layers 3-11 (layer 10: D3 0.0421 -> 0.0394, D21
  0.0639 -> 0.0580).
- **The adjoint kappa3 channel improves.** Its output content falls 10% (0.725 -> 0.650, units 1e-8 n). Per source
  layer it falls 10-17% at layers 11-15, even where the L2 error of D3 rose (layers 14-15).
- **The variance channel worsens.** It rises 7% (1.31 -> 1.40), mostly at layers 8-13 and 15, although the variance
  and covariance L2 errors fell slightly.

Section 3d's G D core failed the same way:
- **The targeted channel.** The kappa4 channel's output content was flat while its L2 error fell 30%.
- **The variance channel** rose 10%.

Two theory-correct local repairs, each verified against the exact law in its own statistic, both lose through the
variance channel. The fitted elements (the lam table, its adaptive rule REF_R and its 0.95 scale) were tuned on the
trajectory of the uncorrected chain, so this is note XXI section 7's compensation, now seen per channel.

**It is not the (3,1) slice acting on the covariance program** (`outputs/k31_covariance_attribution.txt`). With every
other input true, the one-step error of the chain's y covariance is the same, 1.4% at layer 1 falling to 0.4%, whether
the z (3,1) slice is the truth, the lam core, the G D core, or absent.

What remains open is the route by which a better kappa3 or kappa4 state raises the variance channel's output content.
The adjoint functional measures each channel's first-order output content to within 2 points (section 2). The
principled next step is therefore to weight corrections by that functional rather than by L2 accuracy: fit the derived
shapes' few amplitudes (the D coefficient per layer, the feedback weights) to the adjoint-predicted output error.
Re-tuning the fitted tables on output would not be a result.

### 3d. In the chain: G D in place of the lam core (V39_KD)

`est_v29.py` switch V39_KD builds G_l D_l at every layer from the chain's own state:
- the post-activation C, mean, variance, kappa3 diagonal and K21 slice from the end of the previous layer;
- the pre-activation covariance, mean, D3 and D21 of the current layer;
- the exact transport of the pair part (six n^3 products per layer).

It replaces the lam core in the three slices, with G_l the network-1 fit (coefficients 2.06 -> 4.72; layer 1 -> 2 uses
0.20). Bits of V39_KD_BITS select the slices. The pair parts and the K4 -> K3 birth feed are unchanged.
- **One step (all oracles on, `outputs/kd_onestep_off1.txt`).** The chain's own slices reproduce the prediction
  exactly. At layers 12-14 the diagonal error falls from 0.27-0.29 to 0.19-0.21, the (2,2) from 0.27-0.29 to
  0.22-0.23, the (3,1) from 0.73-0.75 to 0.59-0.64.
- **Free-running (`outputs/kd_free_off1.txt`).** The kappa4 diagonal error is lower at every layer from 2 on (0.315 ->
  0.201 at layer 14). The (2,2) slice error falls 0.32 -> 0.23 and the kappa3 readouts improve slightly (D3 0.051 ->
  0.048). The variance and off-diagonal covariance errors rise 5-7%.
- **Output.** The output is worse: net 1 2.1266e-8 -> 2.3086e-8 (+8.6%), net 0 2.2530e-8 -> 2.4873e-8 (+10.4%), at +8%
  of the budget. Every slice subset is worse (`outputs/kd_ablation.txt`):
  - diagonal only +6.5% / +1.8%;
  - (2,2) only +2.9% / +10.5%;
  - (3,1) only +5.0% / +1.1%;
  - the output worsens monotonically as the amplitude falls (0: +71%/+84%; 0.5: +24%/+33%).
- **Why, from the adjoint functional (section 2, rerun on the KD dump).**
  - In absolute terms the kappa4 channel's output content is unchanged: |Delta_k4|^2 = 0.760 against 0.759, in units
    of 1e-8 n. A 30% lower L2 error of the kappa4 diagonal has the same projection on the output.
  - The kappa3 channel improves slightly (0.725 -> 0.690).
  - The variance channel worsens by 10% (1.31 -> 1.45): the changed (2,2) and (3,1) slices reach the next covariance
    through the pair program's use-side terms.

This is note XXI section 7's result again, now located:
- **What it is not.** The compensation is not with the pair closure; the closure is exact (section 3).
- **What it is.** The output reads the kappa4 sector only through two projections: the mean map's weights on the
  diagonal, and the covariance use-side terms of the next pair program. L2 accuracy of the slices is the wrong target.
- **The output-relevant part.** It lies in what D does not capture. The remaining third of R is the theory's
  source 1 (the leg-fed class) and second-order products, which section 3c tests directly.
- **What would not be a result.** Tuning the amplitude on the output would only re-fit lam in another shape.

## 4. Theory: the input radius cancels, the gain mode survives, and it lives in the trace sector

**The radius factorises the whole stack.** The official networks are bias-free ReLU MLPs with x ~ N(0, I_n). Positive
homogeneity gives z_l(x) = r z_l(u) at every layer, with r = |x|/sqrt(n) independent of the direction u (uniform on the
sphere of radius sqrt(n)), and y_l = r relu(zeta_l) likewise. By the law of total cumulance,

    kappa4(z_l) = E[r^4] kappa4(zeta_l) + Var(r^2) Sym3(Sigma, Sigma) + Cov(r, r^3) Sym4(mu, kappa3) + O(n^-2),

with Var(r^2) = 2/n and Cov(r, r^3) = 3/(2n) + O(n^-2) (exact gamma-function values in `code/yclasses.py`). The other
cumulants of the conditional moments are O(n^-2). The bracket is the radial tangent of note XXXII section 4 at
Var(s) = 1/(2n). This identity is suggested by Euler's relation for bias-free ReLU networks in the deep-research
synthesis (its finding 10; derived, not verified there). It is elementary.

**The radius is not a source in the Gaussian-input chain.** The direction carries the compensating law: the uniform
law on the sphere has kappa4(u) = -(2/(n+2)) Sym3(I, I). Since z_0 = W_0 x is exactly Gaussian, the two parts cancel
at layer 0. Every layer map commutes with the radial mixture (homogeneity; at first order, note XXXIV section 2b), so
they cancel at every depth. Put as a perturbation statement: the radial tangent at layer l is the network's exact
response to excess scale variance of the input, and a Gaussian input has none. Candidate D (the radial tangent,
transported, minus its pair-supported part) is therefore not a predicted component at coefficient 1.

**Relation to notes XI, XIX and XXI.** What follows re-derives, with the radius made explicit, the gain-mixture
picture those notes built and measured:
- **Note XI.** GAC carries z = G y with Var G^2 accumulating critically (tau = 1) from a weights-only injection.
- **Note XIX.** It measured the scale's transport at eigenvalue one.
- **Note XXI.** Its power-sum class decomposition of the diagonal found:
  - the (3+1) and (1+1+1+1) classes zero within noise;
  - the (2+1+1) class the whole dropped content, equal to the mixture's 6 g sigma_diag^2 sigma_off^2 at g = the next
    layer's normalised g_4, to 10-25%;
  - the one-loop generation of dropped classes about zero.

Its section 7 is the result that governs any use of this: the derived per-neuron shape (fitted amplitude) and the
derived amplitude (g from D3) both made the chain worse on 16 of 16 networks, +9% and +50%. The fitted lam also
absorbs compensating errors of the pair program.

What is new here:
- **The radius is accounted for.** The cancellation says G_l starts at zero at layer 0.
- **The mean terms are in the shape.** The shape is the raw-vector radial tangent, including the
  Cov(s, s^3) Sym4(mu, kappa3) term the centred mixture lacks.
- **The amplitude comes from the trace sector.** It is the connected Var|y|^2, not D3 or g_4.
- **Section 3's split measures what section 7 could not.** It separates the omitted classes from the closure error
  in all three slices, so it shows whether lam's fitted value is compensating closure error. It also checks note
  XXI's finding of zero one-loop generation (candidate C) per neuron and in every slice.

**What survives: generated scale, transported exactly.** Homogeneity holds from any layer onward: z_m = F_m(y_l) with
F_m positively homogeneous of degree one. A fluctuation of the scale of y_l over input directions, generated by the
gates at layer l, is therefore carried to every deeper layer exactly, with the radial-tangent shape. That is
eigenvalue one, the Ward identity at finite amplitude. Scale generated at different depths accumulates. The prediction
for the omitted classes is

    R_l ~ G_l D_l,   D_l = W#[Var(s^2) Sym3 + Cov(s, s^3) Sym4]_(triple, quadruple) per unit Var(s^2),

with one scalar G_l per layer: the accumulated excess scale variance of y_l. G_l must be visible in the trace sector.
For y = s eta with eta Gaussian, the connected part of Var|y|^2 is Var(s^2) (E|y|^2)^2 to first order. So the
parameter-free estimate is

    g_conn(l) = [sum_ab kappa(y_a, y_a, y_b, y_b) + 4 sum_ab m_a kappa(y_a, y_b, y_b)] / (E|y_l|^2)^2,

the connected part of Var|y_l|^2. It is pair-supported (the K4 diagonal and (2,2) slice, the kappa3 diagonal and D21
slice, the means), so it lies in the chain's retained state. Three predictions, tested in section 3:
- **P1.** D's shape fits R with one coefficient per layer.
- **P2.** That coefficient times Var(r^2) equals g_conn(l).
- **P3.** g_conn grows with depth.

If they hold, the omitted classes have a derived form G_l D_l. Whether the chain can use it depends on note XXI section 7's
compensation, which section 3's closure / classes split measures directly. The full form of D at z_(l+1) needs only
mu, Sigma and the kappa3 diagonal and D21 slice of z_(l+1). Its pair part enters before the transport the chain
already does:

    chain' = T_pair(y slices - G Dpair) + G Dfull.

The cost is O(n^2) beyond the existing transport.

**Where this sits in the literature** (deep-research run of this session; sources verified 3-0 unless noted):
- **ARC's cumulant-propagation estimator** (Wu, Lecomte, Winer, Robinson, Hilton, Christiano 2026, arXiv:2605.05179,
  section 4.1). For odd truncation order it found that it must track the full trace of the next cumulant,
  sum_ij kappa[X_i, X_i, X_j, X_j], the connected part of Var|X|^2. Without it the normalised error is O(1). g_conn is
  that trace with the mean term that ReLU means force. The augmented variant keeps all of kappa4 except its traceless
  (harmonic) part, the complement of the omitted classes here and of THEORY_3's harmonic residual.
- **The finite-width theory's collective observable** (Favaro, Hanin, Marinucci, Nourdin, Peccati 2025; Hanin 2024).
  Under random weights the next layer is conditionally Gaussian given the previous one, and kappa(a,b,c,d) =
  Cov(Sigma_ab, Sigma_cd) + 2 perms. The four-point vertex is the variance of the collective covariance, and for ReLU at
  He init it grows linearly in depth (Roberts-Yaida eq. 5.113-5.120: V/(n K^2) = 5(l-1)/n). That theory averages over
  weights, so it has no third cumulants and no mixed-neuron classes. G_l is its fixed-weight, input-averaged analogue,
  with the mean term the weight average removes.
- **The Stein kernel** (Nourdin-Peccati; Bierme-Bonami-Nourdin-Peccati 2012). For centred F,
  kappa4(F) = 3 Cov(F^2, Gamma_1(F)) with Gamma_1(F) = <DF, -D L^(-1) F>. For a pre-activation, DF is its input-Jacobian
  row, so Gamma_1 is a Jacobian-overlap kernel along the Ornstein-Uhlenbeck interpolation. The scale mode is the part
  of Gamma_1 common to all neurons. This bridge is the synthesis's derived finding 9, not verified; it is not used here.
- **Super-sample covariance in cosmology** (Takada-Hu 2013; Barreira-Schmidt 2017). The squeezed trispectrum is
  response x soft-mode variance x response, rank one for the isotropic background. G_l D_l is that structure with the
  gain as the background mode. The extension to a k-dimensional soft mode (the 64-dimensional common mode of note XX)
  is the natural next term if P1 leaves a structured residual.
