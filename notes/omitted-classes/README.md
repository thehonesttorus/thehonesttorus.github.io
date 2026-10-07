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

*(Results pending: section to be completed from `outputs/yclasses_*.txt`.)*

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
