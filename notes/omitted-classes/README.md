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

Two theory-correct local repairs, each verified against the exact law in its own statistic, both raise the final error,
and in both the adjoint's variance channel grows.

**Correction (checkpoint J).** An earlier version of this paragraph called this note XXI section 7's compensation through
the fitted tables. That causal reading is not established. The evidence shows an adverse signed terminal effect, and that
has other sufficient explanations, none involving a fitted element:
- **Ordinary signed cancellation.** Removing one of two opposing errors raises the total.
- **Recentering.** A mean change compels a variance change, v = q - m^2.
- **A closure derivative.** It can be wrong along the free-running trajectory while the closure's values on true inputs
  are accurate.

Moreover, the size of a channel's "content" depends on the chart: transforming between (mean, covariance) and
(mean, raw second moment) coordinates moves 2 A m mdot between the mean and variance attributions.

Section 3h runs checkpoint J's discriminating protocol on the paired arrays, starting with the exact raw-moment /
recentering split and a closed signed budget.

**It is not the (3,1) slice acting on the covariance program** (`outputs/k31_covariance_attribution.txt`). With every
other input true, the one-step error of the chain's y covariance is the same, 1.4% at layer 1 falling to 0.4%, whether
the z (3,1) slice is the truth, the lam core, the G D core, or absent.

The route by which a better kappa3 or kappa4 state raises the final error is found in section 3h. (An earlier version
proposed setting the derived amplitudes by the adjoint-predicted output error. Fitted against final truth, that is output
calibration, not a derivation. The adjoint's right use is to locate the adverse route and to test independently derived
coefficients.)

### 3h. Checkpoint J's protocol: the feedback repair loses a cancellation, the G D core misreads its own slice

Free-running dumps, at the same cost, of three chains on networks 0 and 1:
- production;
- y2 = V40_FB_SY = 2, the derived feedback weight of section 3f;
- KD = V39_KD = 1, the G D core of section 3d.

The comparisons are paired: each repair against production on the same network, with no amplitude fitted to anything.

**Steps 1-2: raw second moment, recentering, and a closed signed budget** (`code/sbudget.py`;
`outputs/protocol_steps12_off1.txt`, `outputs/protocol_all_off0.txt`).
- **Recentering identity.** The per-layer identity Delta var = Delta S2 - (mu_1 + mu_0) Delta mu (and q = v + m^2
  post-activation) splits the variance channel into a raw-second-moment part and a recentering part.
- **Allocation.** With d = e_1 - e_0, the allocation <tau_X, e_0 + e_1> of each adjoint channel sums exactly to
  n Delta MSE. The check 2<e_0, d> + |d|^2 = Delta holds to print precision, and the tangent remainder is |rho|/|d| = 0.03.

| % of the baseline n MSE | net 1 y2 | net 1 KD | net 0 y2 | net 0 KD |
|---|---|---|---|---|
| Delta MSE | +9.82 | +8.56 | +0.51 | +10.48 |
| var channel | +6.09 | +4.32 | -1.36 | +2.80 |
| - raw second moment | +8.22 | +4.34 | -1.50 | +6.24 |
| - recentering | -2.13 | -0.02 | +0.14 | -3.44 |
| kappa3 channel | +3.69 | +1.41 | +2.01 | +1.27 |
| kappa4 channel | +0.10 | +2.84 | -0.29 | +6.26 |
| tangent remainder | -0.05 | -0.02 | +0.14 | +0.15 |

Recentering is not the adverse route. Wherever the variance channel is adverse, the raw second moment carries it, and
recentering opposes it.

**Steps 3-4: first entry into the second-moment block, and the closure derivative there.**
- **Off-diagonal entry** (`code/sbudget.py` step 3, `outputs/protocol_step3_off1.txt`). The raw second-moment change of
  z_l enters through the off-diagonal second moments (net 1: +8.05 of the y2 change, +4.34 of the KD change).
- **The unary route is null.** The kappa3 -> q, kappa4 -> q and (mu, var) -> q entries together carry at most 0.14%.
  Their Gaussian coefficients (checkpoint J proposition 3.1) reproduce the chain's Delta q to 1-10%, with correlation
  at least 0.994.
- **The closure derivative is physical** (`code/covroute.py`, `outputs/protocol_covroute_off1.txt`). The chain's change
  of the post-activation covariance is reproduced by the gated coefficients at the true reference (note XXXIV) to
  2-3%, with correlation 1.000 at every layer. The residual's allocation is at most 0.14% of the MSE.
- **The baseline error routes are all adverse** (`code/errroute.py`). Split the same way against Monte Carlo truth, the
  baseline's covariance error routes do not cancel one another: C +30.5, D21 +14.7, mean products +17.1, K31 +5.1
  (net 1).

These steps rule out recentering and a closure-derivative defect. They read the state accumulated at each layer,
however, so they cannot say where a change first enters, nor whether its terminal effect is a cancellation.
- **The one-step D21 projection is not the full response** (`code/gradd21.py`). Projected on the one-step
  output-relevant direction, the y2 change of D21 is adverse on net 1 (+0.75) and favourable on net 0 (-0.92), with
  the same sign free-running and one-step.
- **What that projection omits.** It counts only the first read of D21, not its continuation through the covariance
  recursion.

**The first-entry ledger** (`code/fentry.py`, `code/pairx.py`; `outputs/fentry_off0.txt`, `outputs/fentry_off1.txt`).
- **State and reference maps.** The observable state is X_k = (m_k, C^y_k), with the full covariance and var_y on its
  diagonal, and Z_k = (mu_k, C^z_k). The reference maps at the Monte Carlo truth are:
  - A_k: Z_k <- X_(k-1), exact (mu = W m, C^z = W C^y W^T);
  - B_k: X_k <- Z_k, the gated transport: dm = Phi dmu + (phi/2s) dvar,
    dvar_y = 2m(1 - Phi) dmu + (Phi - 2m phi/2s) dvar,
    dC^y_ab = (Phi_a Phi_b + rho_a rho_b C_ab) dC_ab + (u_a Phi_b + Phi_a u_b) C_ab.
- **First entries.** For a run, Delta = run - truth, with
  N^z_k = Delta Z_k - A_k Delta X_(k-1) and N^y_k = Delta X_k - B_k Delta Z_k.
- **The ledger is exact.** Telescoping gives Delta m_out = sum_k (prop N^z_k + prop N^y_k) with no remainder: N absorbs
  every nonlinearity. Measured: |sum - e|/|e| = 2e-15.
- **Routes.** N^y is split by route with the gated coefficients:
  - the kappa3 and kappa4 diagonals, into m, var_y and (through Phi) C^y;
  - D21, K22 and K31;
  - the residual, the chain's response minus the reference linearisation.

  N^z is the representation step (0.1% of Delta C^z).
- **Forward and adjoint.** Forward images delta_X give e = sum_X delta_X. The full mean-covariance adjoint
  (Lambda_(k-1) = W^T [Lambda_off o K + diag(beta)] W) gives per-layer entries, agreeing with the forward images to
  1e-16. The Monte Carlo halves give noise-free Grams.

**1. The production chain's error enters through its cumulant state, not its closure.** Route shares
<delta_X, e>/|e|^2:

| | kappa3 diag | kappa4 diag | D21 | K22 | K31 | residual |
|---|---|---|---|---|---|---|
| net 1 | 0.254 | 0.243 | 0.306 | 0.016 | 0.166 | 0.012 |
| net 0 | 0.217 | 0.277 | 0.274 | -0.007 | 0.193 | 0.059 |

The residual includes the reference's Monte Carlo noise; the representation share is at most 0.012.
- **Where the error enters.** The error enters mostly at layers 7-15: the kappa3 diagonal and kappa4 diagonal at the
  last three layers, and D21 and K31 from layer 5 on.
- **Overlaps.** The route images overlap weakly and mostly negatively: noise-free Gram off-diagonals are at most 0.06
  of |e|^2. The residual and representation routes are dominated by the reference's noise.

**How the comparisons split.** For a comparison 0 -> r:
- n Delta MSE = sum_X own_X + sum_(X<Y) pair_XY;
- own_X = Delta|delta_X|^2 and pair_XY = 2 Delta<delta_X, delta_Y>;
- each term is <tau_X, delta_Y^0 + delta_Y^r> with tau = delta^r - delta^0.

tau does not depend on the truth, so Monte Carlo noise enters only through delta^0 + delta^r. The noise bar is
|half 0 - half 1| / 2.

**2. y2: the repair improves its own routes and loses a cancellation** (% of the baseline n MSE).

| | net 1 | net 0 |
|---|---|---|
| Delta MSE | +9.82 | +0.51 |
| own: kappa3 diag | -3.23 +- 0.70 | -3.71 +- 1.15 |
| own: D21 | -2.37 +- 0.56 | -1.31 +- 1.03 |
| own: all routes | -5.09 | -4.86 |
| pairs: all | +14.91 | +5.37 |
| pair kappa3-D21 | +6.25 +- 0.77 | +2.17 +- 0.87 |
| pair kappa3-kappa4 | +3.74 +- 0.53 | +3.99 +- 0.15 |

- **The repair works in its own routes.** On both networks the output images of the kappa3-diagonal and D21 errors
  shrink (net 0's D21 only at 1.3 sigma).
- **The final error still rises.** The change in one route lines up with the error left in another, so the pair
  terms grow.
- **The opposition it removes was helping.** In the baseline the kappa3-diagonal image opposes the D21 image (noise-free
  Gram -0.05 of |e|^2 on net 1) and the kappa4 image (-0.03). The repair removes much of that opposition:
  - on net 1 the combined kappa3-state image |delta_k3 + delta_D21|^2 is unchanged within noise
    (-3.23 - 2.37 + 6.25 = +0.7 +- 1.2);
  - the remaining rise is the kappa3-kappa4 pair (+3.7), then kappa4-D21 and D21-K31 (+1.8, +1.2), plus noisier pairs
    with the residual route.
- **The common pair.** On both networks the kappa3-kappa4 pair alone costs +3.7 to +4.0.
- **Where it first enters.** On net 1 the first adverse entry is D21 and the kappa3 diagonal at layers 7-8
  (+0.51, then +2.23 and +1.09). Layers 1-6 are net favourable (cumulative -0.54).
- **Net 0.** The per-layer kappa3 entries alternate in sign from layer 10 (-0.56, +0.78, +1.13, -1.51, +1.73, -1.44).

This is checkpoint J's first alternative, ordinary signed cancellation. It is now established in the first-entry chart,
where every term is located and the budget has no remainder.

**3. KD: the G D core's slices are closer in L2 and worse in what the output reads.**

| | net 1 | net 0 |
|---|---|---|
| Delta MSE | +8.56 | +10.48 |
| own: K31 | +5.04 +- 0.11 | +2.64 +- 0.73 |
| own: kappa4 diag | +1.92 +- 0.54 | +2.26 +- 0.24 |
| own: all routes | +6.17 | +5.82 |
| pairs: all | +2.39 | +4.66 |
| pair kappa3-kappa4 | +4.96 +- 0.07 | +0.46 +- 0.83 |
| pair D21-K31 | +2.57 +- 1.08 | +2.93 +- 1.07 |

- **The slices are closer in L2.** The (3,1) slice and the kappa4 diagonal are closer to truth in L2: one-step (3,1)
  error 0.73 -> 0.59; diagonal error -30%.
- **Their output images grow.** The own terms are +6.2 and +5.8, more than half of the change on both networks. The
  G D shape errs in the direction the output reads.
- **The pairs add to it rather than offset it.** The main pairs are kappa3-kappa4 on net 1 and D21-K31 on both.
- **Where it first enters.** On net 1 the K31 route turns adverse at layer 5 (+1.60), and the cumulative turns
  positive at layer 8.

**What the protocol settles.**
- **Ruled out.** Recentering, and a closure-derivative defect at the entry.
- **y2.** The failure is cross-route cancellation among the cumulant-state routes. The repair is locally right in its
  own routes.
- **KD.** The failure is a local error in the readout direction of the repaired slices.
- **Fitted-table compensation is not needed for either.**
  - Whether the anti-alignments that y2 removes were produced by a fitted element is a separate question.
  - The pair that recurs, kappa3-kappa4 in three of the four comparisons, involves the kappa4 diagonal. Its (2+1+1)
    stand-in is lam s_off^2, one fitted scalar per layer (system-audit note, section 3).
  - The discriminating test is this ledger with that class at a derived amplitude. That means the D transport of
    section 3 at the derived gain g_conn (the table's "lam -> g_conn D" column), not a coefficient fitted to any output.
- **Step 5 is not needed.** The ledger closes exactly at one step and the residual route carries at most 0.2%, so
  checkpoint J's step 5 (the discarded two-step return) adds nothing here.

**For the estimator.**
- **The kappa3 readout is joint.** The output reads the kappa3 state through the sum of its diagonal and D21 images, and
  that sum overlaps the kappa4 image.
- **Repairs must be paired.** A repair of either cumulant sector must be judged on that joint readout, not per slice and
  not in L2. A kappa3 repair has to be paired with the kappa4 error it currently cancels. That cancellation happens in
  the propagated images, not inside one unit's query (section 3i).

### 3i. Checkpoint K's degree-one state, measured: the physical-cumulant queries are adequate; the loss is in transport

**Checkpoint K's reformulation.** For a bias-free ReLU network every future mean is a functional of one object:
the radius-weighted angular law nu^(1) of the hidden state. After ReLU that object is equivalently the whole support
function h_X(w) = E (w^T X)_+.
- **What the chain already computes.** Each next-layer mean m_(k,i) = E relu(z_(k,i)) is exactly such a degree-one
  query, at w = row i of W_k.
- **How it computes it.** It takes the query from the physical cumulants of the marginal through the Edgeworth map of
  `est_v29` (TERM_SPECS[(1,)]):
  Gauss(mu, s) + k3/6 E relu''' + k4/24 E relu'''' + k3^2/72 E relu^(6) + k3 k4/144 E relu^(7).
- **Checkpoint K's question in these terms.** Does that physical-cumulant state lose degree-one information that the
  suffix queries?

**A refinement of the radial family** (checkpoint K, theorem 5.1). Write the hidden gain as s = 1 + eps, with z = s A
for a Gaussian A independent of eps.
- **Gain variance.** Along v = E eps^2, kappa4 = 12 v sigma^4 + O(v^2), and kappa6 is O(v^2) for a symmetric eps. The
  kappa<=4 map is first-order exact along this direction.
- **Gain skewness.** Along tau3 = E eps^3, at fixed gain mean and variance (a radial-family direction), the cumulant
  generating function gains tau3 (sigma^4 t^4 / 2 + sigma^6 t^6 / 6). So kappa4 = 12 tau3 sigma^4 and
  kappa6 = 120 tau3 sigma^6 appear at the same order.
- **What that does at mu = 0.** The kappa4 term of the mean is -tau3 sigma phi / 2 and the kappa6 term is
  +tau3 sigma phi / 2: they cancel, as the invariance requires.
- **Consequence.** The truncated map's response along such a direction is spurious: its remainder is H = -T4, with
  coefficient -1.

The invariance therefore bites only through skewed or heavy-tailed gain. For a smooth gain at width n that is second
order (tau3 ~ v^2 ~ n^-2).

**The measurement** (`code/kquery.py`; `outputs/kquery_off0.txt`, `outputs/kquery_off1.txt`; Monte Carlo halves for
noise-free norms).
- **(0) Calibration.** The map reproduces the chain's own means from its dumped state to 0.4-0.6% of the chain's
  error.
- **(a) The truncation remainder of the true queries.** H_k = E relu(z_k) - map(true mu, var, k3, k4) is what an exact
  kappa<=4 state still misses. Rms per unit:
  - layers 8-15: |H| ~ 1e-5, against |true - Gauss| ~ 5e-4 and |T4| ~ 1.7e-4, so 2% of the non-Gaussian part of the
    query;
  - the cross-half coefficient of H on T4 is -0.01 to -0.02, with cos(H, T4) -0.2 to -0.47 and cos(H, T3) -0.3 to
    -0.7 on both networks.
  - **Interpretation.** Higher orders cancel part of the kappa4 and kappa3 terms, with the radial family's sign, at
    about 1.5% of the kappa4 term.
- **What exact kappa3, kappa4 would leave.** The output image of -H (full mean-covariance propagation) has energy 3.2%
  (net 0) and 3.5% (net 1) of the baseline MSE. Its share of each run's error is at most 1.4%.
- **(b) Inside one query.** Split the chain's query error as N_m = P3 + P4 + R, with P3 = c3 dk3 and P4 = c4 dk4.
  The kappa3 and kappa4 errors are nearly orthogonal across units (cos -0.02 to -0.13 at layers 4-15); within-query
  compensation is 2<P3, P4> = -6 to -7% of |P3|^2 + |P4|^2. The rest R matches |H|, as it should.

| sum over layers of the query errors (noise-free) | net 1 | net 0 |
|---|---|---|
| y2: Delta |P3|^2 | -11% | -10% |
| y2: Delta |P3 + P4|^2 | -3.8% | -3.8% |
| KD: Delta |P4|^2 | +6.5% | +4.7% |
| KD: Delta |P3 + P4|^2 | +4.5% | +3.9% |

**What this settles.**
- **The coordinates are not the bottleneck.** With exact inputs the physical kappa<=4 state is an adequate coordinate
  system for the degree-one queries the network makes: 2% of their non-Gaussian part, and a 3-4% floor of the current
  MSE.
- **The radial-family mechanism is real but second order.** It is present with the predicted sign at about 1.5% of the
  kappa4 term.
- **The loss is in transport.** The chain's error is the inaccuracy of the propagated cumulant state: the 97% of the
  first-entry ledger in section 3h. In checkpoint K's language, degree-one information is lost in how the state is
  transported, not in how it is read.
- **y2 fails in the signed pairing.** y2 improves the queries themselves on both networks, with no within-query
  compensation lost, yet the output worsens by up to +9.8%. The adverse effect is in the signed suffix pairing of query
  errors across units and layers: checkpoint K's telescoping identity (section 9), which is this note's first-entry
  ledger.
- **KD misreads its own query.** KD worsens the kappa4-weighted query error (c4 is largest where |alpha| is far from 1)
  while its L2 error falls 30%.

### 3j. Quotient before truncation: exact null sources, and what they say about the chain's kappa4 sector

**The identity** (supplied THEORY.md; checked here independently). Take a Gaussian reference N(mu, Sigma), a symmetric G,
normalised symmetric products, and the source responses L_k[T]F = E[T : nabla^k F] / k!. Then for every positively
homogeneous F of degree p,

    L4[Sigma.G] F = ((p - 2)/12) L2[G] F - (1/4) L3[mu.G] F.

- **Proof.** g = G : nabla^2 F has degree p - 2, so E[Sigma : nabla^2 g] = E[(Z - mu) . nabla g] = (p - 2) E g - E[mu . nabla g]
  (Gaussian integration by parts, then Euler). THEORY.md states the case p = 1.
- **The two null tuples.**
  - p = 1, every remaining bias-free ReLU suffix, hence every future mean: (dSigma, dk3, dk4) = (G/12, mu.G/4, Sigma.G)
    is invisible.
  - p = 2, the post-activation covariance E[y_a y_b] the chain carries: the invisible tuple is (0, mu.G/4, Sigma.G),
    with no covariance part.
- **Numerical check** (`code/verify_null.py`). For a two-layer ReLU net in 2D at a correlated, non-centred reference,
  using sector-wise quadrature of F times the Hermite densities:
  - degree 1: null combination 7e-18, against individual terms of 0.002-0.10;
  - degree 2: -6e-17;
  - controls: swapping the two tuples gives 8e-3 and 4e-3.

**Consequences for the chain's channels.**
- **The readout respects the quotient.** Per unit, the Edgeworth coefficients satisfy the p = 1 relation identically:
  c4 s^2 = -cv/12 - c3 mu/4. The variance coefficients satisfy the p = 2 relation: q4 s^2 = -q3 mu/4.
- **The intermediate covariance is not quotient-invariant.** A degree-1 null source still changes the post-activation
  covariance, by L2[G]/12 of the degree-2 observable. Exact dynamics cancels this downstream through the induced
  higher post-activation cumulants. A chain whose covariance channel is exact but whose kappa3/kappa4 transport is
  approximate can therefore break the equivalence. The candidate breaking points are in transport: the kappa4 core
  transport (dG with the fitted lam), the K4 -> K3 birth feed, and the leg truncations.
- **The chain's kappa4 sector is a trace core.** In production (V33_K4Q = 3) the regenerated core is
  G = diag(dG) + lam C_off, with slices
  - g4row = 2 dG,
  - wk4m_ab = (dG_a + dG_b)/3,
  - B_ab = lam C_ab.

  Since var ~ 2 under He scaling, these match the slices of Sigma.G:
  - diag var_a G_aa,
  - (2,2) (var_a G_bb + var_b G_aa + 4 C_ab G_ab)/6,
  - (3,1) (var_a G_ab + C_ab G_aa)/2.

  To first order the whole kappa4 sector is therefore a covariance shift -G/12 plus a kappa3 shift -mu.G/4. The
  exceptions are three:
  - the Euclidean metric METRIC_C = 2 in place of Sigma;
  - the diagonal-only quenched term k4corr;
  - the non-trace shape that KD introduces.
- **The gain mode is invisible.** A mean-one scale mode z = S Z0 (E S = 1, Var S = t) has the law tangent
  (t (Sigma + mu mu^T), 6 t mu.Sigma, 12 t Sigma.Sigma). It is the sum of the k = 4 null with G = 12 t Sigma and the
  k = 3 null (t mu mu^T, 3 t mu.Sigma), so it is invisible to every future mean. Per unit,
  cv (s^2 + mu^2) + 6 c3 mu s^2 + 12 c4 s^4 = 0 identically (this is section 3i's first-order exactness along the gain
  variance).
- **What follows for the chain.** The physical covariance carries the network's gain mode exactly. Its effect on the
  means is cancelled only if the chain's kappa3 and kappa4 carry the same gain amplitude t. Any mismatch between the
  amplitudes implied by the two channels is a spurious, visible error, of 6 c3 mu var dt3 + 12 c4 var^2 dt4 per unit.

**Staged tests** (the experiment VMs are stopped; Azure subscription in Warned state).
- **(A) Continuation consistency on the existing arrays** (`code/nullalign.py`). Per layer, for the truth and each chain
  run:
  - which metric makes the kappa4 slices a trace core (covariance metric against the chain's 2I);
  - the gain amplitudes fitted separately on the kappa4 slices (diag, (2,2), (3,1)) and the kappa3 slices (diag, D21);
  - the mean error that the chain's gain-amplitude inconsistency predicts, against its actual query error P3 + P4;
  - how far its lower-order errors line up with the null companions (dG/12, mu.dG/4) of its kappa4 error trace core,
    and how y2 and KD move that.

  The fitter recovers an exact trace core to 3e-16.
- **(B) The null-source audit** Omega(G) = R4(Sigma.G) + R3(mu.G)/4 + R2(G)/12. This injects the full tuple at one
  layer and measures the chain's output response, whose exact value is 0; no truth is needed. It needs injection hooks
  that place mu.G/4 in the kappa3 state (the readouts and the legs), not only in the readouts.

### 3k. The staged tests, run: the error is transport, the truth's kappa4 is a consistent gain mode, and the chain's transport is nearly gauge-covariant except for rough gauge fields

Everything here ran on Modal (`infra/modal/mjob.py`), after the Azure and AWS pools were dropped.
- **Monte Carlo.** mc2 (seeds 301/302) and mc4 (401/402, post-activations) were regenerated in 16 independent chunks per
  pass (`code/mcchunk.py`, `code/mcmerge.py`): 1.6e7 samples per pass, two halves of 8e6.
- **Comparisons.** The output comparisons use the dataset's own all-layer truth.
- **Baseline.** Production on nets 0-31 has mean raw MSE 2.281e-8 (standard error 0.045e-8). Nets 0 and 1 are typical.

**(1) Where the error is created: intervention telescoping** (`code/itele.py`, consistent oracles V41_ORC_CONSIST = 1,
`outputs/itele_off01.txt`). Oracles D3 D21 G4 WK4M K31 VAR COFF MU are applied at layers 0..k, one half of the Monte
Carlo per run. The step effect is <e_(k-1)^A, e_(k-1)^B> - <e_k^A, e_k^B>, noise-free.

| % of the free n MSE | layers 0-6 (cumulative) | each of layers 7-15 | all layers oracled |
|---|---|---|---|
| net 0 | +4.2 | +7 to +15 | -1.2 |
| net 1 | +2.6 | -1 to +20 (9-15: +5 to +20) | -2.8 |

- **Transport, not readout.** With true inputs at every layer the noise-free output error is zero within the
  oracle-noise bias. That bias appears as net 1's -11 +- 3.5 at layer 0, where the chain is exact.
- **Where.** The error is created uniformly over layers 7-15. The suffix kernel's prior (`../ray-compiler`, section 2)
  puts the higher-order weights there.
- **Endpoint.** The consistent and inconsistent all-oracle endpoints agree to print precision, as they must: when every
  layer's readouts are replaced, the legs are never read.

**(2) The readout at true inputs, all channels** (`code/cread.py`, `outputs/cread_off{0,1}.txt`; output images of the
defect, % of the baseline n MSE).

| | mean | var_y | C_off | all, energy | all, share of e |
|---|---|---|---|---|---|
| net 0 | +9.0 (1.7) | +0.5 | +4.9 (1.1) | +4.0 (3.2) | +6.0 +- 3.2 |
| net 1 | -1.3 (1.6) | -0.3 | -5.1 (1.8) | +1.3 (1.6) | +0.9 +- 1.4 |

The second-moment readout carries at most a few percent. kquery on the new Monte Carlo replicates section 3i: |H| ~ 1e-5
at layers 8-15 against |true - Gauss| ~ 5e-4, and the image of -H is 3.2% of the MSE on net 1 (`outputs/kquery_modal_off*`).

**(3) The truth's kappa4 is the gain mode, with one amplitude across channels** (`code/nullalign.py`,
`outputs/nullalign_off{0,1}.txt`). Gain amplitudes are fitted per slice to the gain-tangent templates.

| layer 14 | t4 diag | t4 (2,2) | t4 (3,1) | t3 diag | t3 D21 |
|---|---|---|---|---|---|
| truth, net 1 | 4.99e-3 (0.93) | 5.02e-3 (0.93) | 5.08e-3 (0.56) | 5.00e-3 (0.92) | 5.00e-3 (0.84) |
| chain, net 1 | 4.22e-3 | 4.25e-3 | 2.88e-3 | 4.91e-3 | 4.94e-3 |
| truth, net 0 | 5.28e-3 (0.95) | 5.26e-3 (0.94) | 5.56e-3 (0.63) | 5.53e-3 (0.93) | 5.48e-3 (0.86) |
| chain, net 0 | 4.36e-3 | 4.44e-3 | 2.98e-3 | 5.44e-3 | 5.41e-3 |

(in brackets: the template's explained fraction)

- **The truth.**
  - Its kappa4 diagonal is 93-98% gain shape at every layer, and its (2,2) slice is 76-85% a trace core in either metric.
  - From layer ~11 on, the kappa3 and kappa4 amplitudes coincide (t3 = t4 within 1-5%). At shallow layers t3 > t4
    (layer 1: 1.6e-3 against 0.8e-3).
- **The chain.**
  - Its kappa3 carries 98% of the true amplitude, its kappa4 diagonal and (2,2) 83-85%, and its (3,1) slice
    (lam C_off) 54-57%.
  - Its (3,1) slice is exactly a 2I trace core, while the truth's is only 37-67% trace core.
- **The per-unit check is null.** The per-unit mean error this inconsistency predicts at the same layer is not aligned
  with the chain's query errors (cos -0.08 to +0.05). The inconsistency does not act through the readout of its own
  layer.

**(4) The null-source audit: close to gauge-covariant for the gain mode, not for rough gauge fields** (V43_NULL_INJ,
`code/nullinj.py`, `code/nullinj2.py`, `code/nullinjv.py`; `../ray-compiler/outputs/nullaudit_off01.txt`).
- **The injection.** At layer l, the exact degree-1 null tuple of a diagonal gauge field G = diag(g), with every
  part the chain represents:
  - covariance: dvar = g/12;
  - kappa3: dD3 = mu g/4, dD21_ab = mu_b g_a/12;
  - kappa4: dg4row = var g, dwk4m = (var_a g_b + var_b g_a)/6, dwk431_ac = C_ac g_c/2;
  - the (2,1,1) entries C_bc g_a/6, through the K4 -> K3 feed.

  A diagonal G has no all-distinct entries. The exact response of every mean is zero.
- **Consistency.** The birth M-block must subtract what the legs actually carry, which is the pre-injection D3 and D21
  (V43_NI_CONS = 1, as V41). A first run without this silently removed the injected kappa3 from the post-activation
  representation, the artifact of stage map 10. That run reported gain-shape ratios of 0.38-0.84, and they are
  withdrawn.
- **The protocol.** Two shapes: g = 12 eps var, the diagonal of the gain mode, and g = 12 eps var o (random signs).
  eps = +-1e-3 with central differences, layers 3, 7, 11, 14, nets 0 and 1. Responses are reported as fractions of the
  covariance part's own response.

| full tuple / covariance part | layer 3 | layer 7 | layer 11 | layer 14 |
|---|---|---|---|---|
| gain shape, net 0 / net 1 | 0.07 / 0.05 | 0.13 / 0.17 | 0.18 / 0.17 | 0.14 / 0.19 |
| random shape, net 0 / net 1 | 0.29 / 0.34 | 0.34 / 0.36 | 0.22 / 0.22 | 0.09 / 0.09 |

- **How covariant it is.** The chain cancels 81-95% of a gain-shaped null perturbation and 64-91% of a rough one. At the
  injection layer the response is 0.01-0.05: the per-unit readout identity holds. The residual is created in the
  transport.
- **Where the residual sits** (net 0, layers 7 and 11, gain | random):
  - **The birth M-block's rank-4 (2,1) residual.** R_RES 64 gives 0.10, 0.10 | 0.22, 0.16; the exact R_RES 1024
    gives 0.10, 0.10 | 0.20, 0.13, against production's 0.13, 0.18 | 0.34, 0.22. The post-activation trace core's
    (2,1) slice, (2 C_ij v_i + var_i v_j)/3, has the full-rank part diag(v) C^y, which rank 4 cannot hold. This is the
    prediction of `../ray-compiler` section 3, confirmed.
  - **The lam feed.** Without the feed the ratios are 0.10, 0.11 | 0.20, 0.15. For rough shapes the injected (2,1,1)
    term is half the needed completion.
  - **The lam machinery is part of the covariance.** lam = 0 gives 0.22, 0.27 | 0.62, 0.48, and freezing the adaptive
    lam rule gives 0.23, 0.25 | 0.34, 0.23. The lam terms are the chain's representation of the gain mode's
    off-diagonal kappa4, and the cancellation needs them.
  - **No effect.** The D21 feedback (V18_NO_FB) does not change the anomaly.
- **Sensitivity.** A 0.1% variance perturbation at layer 3 moves the output by 1.5 times the baseline error. The scale
  of any null component of the chain's state error is therefore large, even at a 5-35% anomaly.

**What (1)-(4) say together.**
- **The chain's state error matters only where it is visible.** The chain is exact at the true state (1). Its readout
  carries a few percent (2). Its transport is 64-95% covariant along exact null directions (4).
- **The truth's kappa4 is mostly null content** (3), the gain mode with one amplitude across channels. The chain
  carries that amplitude at 83-85% in the kappa4 diagonal and (2,2), and 54-57% in (3,1).
- **The deficit is visible.** Since the chain is nearly covariant for the gain shape, the deficit acts like a visible
  covariance-plus-kappa3 error, its null companion. Raising the amplitude alone moves the output by that visible
  equivalent, which other errors currently offset. This is why every derived amplitude made the output worse (derived
  3 g var^2 +53%, x1.19 +58%, KD +9%, net 0): the fitted lam is a compensating counterterm, as note E11 said. The first
  audit's reading that it compensates a gain-mode anomaly is withdrawn.
- **For the estimator.** The anomaly that remains has two located sources, the rank-4 Rres and the lam feed's missing
  per-unit shape. Whether removing it lowers the error is an output question (section 3l).

### 3l. The output tests of section 3k's leads (32 networks, paired against production)

All runs used nets 0-31 with the dataset's own truth, paired per network against the production baseline (mean raw
2.281e-8). FLOPs are from flopscope; the residual wall time is machine-noisy and left out of the comparison.

| variant | raw MSE change (per-net mean +- s.e.) | better on | FLOPs (B) | MSE x FLOPs |
|---|---|---|---|---|
| production (R_RES 4) | 0 | - | 0.2131 | 4.86e-9 |
| birth M-block rank R_RES 16 | -1.38% +- 0.41 | 25/32 | 0.2190 | +1.4% |
| R_RES 32 | -1.98% +- 0.52 | 24/32 | 0.2268 | +4.3% |
| R_RES 64 | -2.17% +- 0.59 | 24/32 | 0.2425 | +11% |
| separable (2,1) projection onto the A leg (V45_ESEP) | +0.36% +- 0.34 | 8/32 | 0.2131 | +0.3% |
| physical-metric (3,1) and feed shapes (V44_PMETRIC) | +0.28% +- 0.18 | 4/32 | 0.2131 | +0.2% |
| physical (2,2) slice (V33_WK4M = 1) | +5.48% +- 0.51 | 0/32 | - | - |
| both of the last two | +8.94% +- 1.09 | 1/32 | - | - |

- **The audit predicted a real accuracy lever.** The rank-4 (2,1) residual is a measured accuracy loss: -2.2% at rank
  64, significant at 3.6 sigma. An earlier output tuning had read the output as flat from rank 4 to 32. The thin legs
  carry R_RES + 2 columns through every layer, though, so the cost grows faster than the MSE falls, and the score is
  worse at every rank tried.
- **The structured fix does not reproduce the gain.** V45_ESEP moves the row-scaled part diag(u) C_off diag(w1) of the
  residual onto the A leg for free; it is the post-activation trace core's full-rank part, the audit's anomaly source.
  It changes nothing. The rank's gain is generic accuracy of the residual, not its gauge structure.
- **Physical shapes lose, as before.** The fitted lam core is a counterterm; physical shapes at unchanged mean
  amplitude make the output worse.
- **Reading.** The chain's free-run error lies almost entirely in visible directions: nullalign (3)-(4), and these
  tests. Gauge covariance is not where the remaining error is, so the quotient programme will not close the gap by
  removing an anomaly. What it does offer is a cost argument: the truth's kappa4 is 93-98% null content, so a chain
  that is covariant by construction need not carry it. That needs a representation in which the gain mode is factored
  out of kappa3 as well, which note E11 found is not where the old-source content lives.

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
