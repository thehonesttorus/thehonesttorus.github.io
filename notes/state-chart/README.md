# The chain's state chart, audited against the coherent-state theory

Working note XXXV. The checkpoint "K4: coherent state, exact transport, and signed terminal response" (NCG-20261007-I)
makes its first executable step a state-chart audit (its section 10.1): are the chain's arrays literal slices, seed
coefficients, standardised cumulants or fitted gains; which traces does K4Q subtract; is the Ritz approximation
zero-diagonalised before its rows are used; does lambda change? The checkpoint could not see the implementation. This
note answers from `est_v29.py` (V33_K4Q = 3, V33_K4Q_RANK = 4, V25_BETA = 1, METRIC_C = 2), with the algebra checked in
`notes/k4-coherent-state/code/check_k4_contract.py`. It then states what the checkpoint's four-cell comparison becomes
in this chart. No experiment.

Notation follows the checkpoint: B(E) is the pure pair lift (B(E)_iijj = E_ij), D(d) the diagonal class, C(M, N) the
normalised product core, J_n the minimum-norm trace lift, P the vector form of G_n on diagonal traces
(P(r) = 6r/(n+4) - 3 (1^T r) 1/((n+2)(n+4))), H = W o W, M = W W^T.

## 1. The audit

**The post-activation arrays are literal cumulant slices.** K4v is the literal kappa4(y_a): the moment-to-cumulant table
gives pk4 - 4 pk1 pk3 - 3 pk2^2 + 12 pk1^2 pk2 - 6 pk1^4. K22 is the literal kappa(y_a, y_a, y_b, y_b), symmetric and
zero-diagonal. There is no (3,1) array at the y level, and no trace array: the off-diagonal trace is never formed, so
the transported y tensor implicitly has U_ij = 0 for i != j. The diagonal trace u = K4v + K22 1 satisfies the
checkpoint's compatibility condition by construction.

**The pre-activation diagonal is in the seed chart, with the residual core at metric 2I.** With Khat the Ritz
approximation of K22 (K4Q_RANK + 8 = 12 pairs) and Khat_0 its zero-diagonal part (the code subtracts diag Khat before
using row sums, so yes, it is zero-diagonalised: rs = Khat 1 - diag Khat), the code's operations give exactly

    g4row = diag[ L_W S + C(2I, N_E) + C(2I, lam W C_y,off W^T) ],   S = D(K4v) + B(Khat_0),
    N_E = W Diag(P(r)) W^T,   r = E 1,   E = K22 - Khat_0.

This is the checkpoint's seed chart (11) with seed S and declared trace U = Diag(u): the residual trace is
U - C_n S = Diag(r), so the residual core N_res is exactly N_E. The physical y tensor that the diagonal transports is
therefore T_(U,S) = D(K4v) + B(Khat_0) + J_n(Diag r), not the literal D(K4v) + B(K22). The two differ by the harmonic
part Q_n B(E) = B(E) - J_n(Diag r), whose literal entries are (checkpoint (13)): diagonal -P(r)_i, pair
E_ij - (P(r)_i + P(r)_j)/6, (3,1) zero. Of the checkpoint's four mechanisms, this diagonal carries the residual-core metric
error C(M - 2I, N_E) and omits the harmonic residual L_W Q_n B(E). Those are THEORY_3's D_metric and
D_known + D_unresolved.

**The pre-activation (2,2) slice is in a different chart: the pure trace core.** wk4m = (dG_i + dG_j)/3 with
dG = H P(u) + lam s_off^2 is exactly the (2,2) slice of C(2I, W Diag(P(u)) W^T + lam W C_y,off W^T): seed zero, the whole
y tensor replaced by its trace core J_n(Diag u). The Ritz content that the diagonal receives (L_W B(Khat_0)) never
reaches this slice, nor does the 4 M_ab N_ab term of the core slice formula.

**The pre-activation (3,1) slice is the lam constituent only.** wk431 = lam C_z,off, the (3,1) slice of
C(2I, lam W C_y W^T) (with C_y where the lam core of the other slices has C_y,off). The trace core's own (3,1) slice,
(W Diag(P(u)) W^T)_ab, and the seed's (3,1) slice, 3 h_a^T Khat_0 c_ab + (W^o3 Diag(K4v) W^T)_ab (checkpoint section 4.4,
h_a = w_a o w_a, c_ab = w_a o w_b), are both absent.

**lam changes every layer.** It is reset to the fitted table value LAM[l] and rescaled once by
clip(mean(dG_0)/mean(var)/ref, 0.5, 2) with ref a fitted per-layer table (REF_R), where dG_0 = H P(u) + lam_table s_off^2 contains the pair core but not the K4Q
correction. A change to the core therefore moves lam, and through it every slice. Note XXXII section 4 recorded that the
rule is not scale-homogeneous.

**Summary of the chart.** Three output slices, three charts: the diagonal is the seed chart (seed D(K4v) + B(Khat_0),
residual core at 2I); the (2,2) slice is the pure core chart (seed zero, core at 2I); the (3,1) slice is the lam core
alone. Every slice also carries the lam core, and lam is adaptive. The y-level inputs are literal. This is the precise
content of "the slices are not the transport of one declared tensor" (note XXXII section 1).

## 2. The four-cell comparison in this chart

The checkpoint's cells are T_ab = L_W(S + b B(E)) + C(M_0 + a dM, N_res - b N_E), with a the metric and b the seed.
Because N_res = N_E here, the residual core vanishes when the seed is completed:

| cell | metric | seed | tensor |
|---|---|---|---|
| 00 | 2I | D(K4v) + B(Khat_0) | L_W S + C(2I, N_E): the chain's diagonal |
| 10 | W W^T | D(K4v) + B(Khat_0) | L_W S + C(M, N_E): residual core transported exactly |
| 01 | (any) | D(K4v) + B(K22) | L_W[D(K4v) + B(K22)]: exact transport of the literal retained tensor |
| 11 | (any) | D(K4v) + B(K22) | the same tensor as cell 01 |

The interaction T_11 - T_10 - T_01 + T_00 = -C(dM, N_E) holds as the checkpoint states, and here it is the whole metric
effect: once the seed is completed, the metric choice has nothing left to act on. Three consequences:

- The chain is not any of the four cells off the diagonal. Its (2,2) and (3,1) slices come from the pure core chart and the
  lam core. A fifth reference is needed: cell 00 with all three slices taken from its tensor. The comparison
  "chain vs cell 00" isolates the slice incoherence; "cell 00 vs cell 01" the residual harmonic content; "cell 00 vs
  cell 10" the metric.
- Cell 01 is the coherent tensor of note XXXII section 2, with slices diag
  Q = W^o4 K4v + 3 diag(H K22 H^T), pair (H K22 H^T)_ab + (H Diag(K4v) H^T)_ab + 2 c_ab^T K22 c_ab, and (3,1)
  ((H(3 K22 + Diag K4v)) o W) W^T, the last two in the checkpoint's own form h_a^T E h_b + 2 c_ab^T E c_ab and 3 h_a^T E c_ab.
- Note XXXIV extends the literal retained tensor beyond what the code computes. The y (3,1) slice is fixed by objects the
  chain holds (a closed form in C_off, D21 and the pre-activation slices), and the theorem of note XXXIV makes it part of
  the same first-order state. A sixth cell, L_W[D(K4v) + B(K22) + S31(B^y)], is the coherent transport of the complete
  pair-supported y tensor. The difference from cell 01 is exactly the transport of a class the chain never formed.

## 3. Where the remaining questions sit

The checkpoint separates representation, missing shape and law closure (its section 10.3). In this chart:

- **Representation.** The difference between the chain and cell 01 (or the sixth cell): slice incoherence, the
  residual harmonic content, the 2I metric. All of it is deterministic in stored arrays.
- **Missing shape.** The (2,1,1) and (1,1,1,1) classes of the y tensor (H_R in the checkpoint's decomposition (4)), with the
  implicit completion U_off = 0. By note XXXIV their first-order sources are the legs' all-distinct third-cumulant
  entries, their own memory (above all the scale mode), and second-order covariance products.
- **Law closure.** The pair program's bivariate closure at a Gaussian reference, and the mean map. The checkpoint's
  theorem 7.1 says a fourth-order state cannot determine the response universally, and its finite-law examples prove it.
  Along the competition trajectory, though, the oracle ceiling of note XXXI bounds this part for the local maps: with
  every pre-activation statistic they read replaced by truth, the residual is 1.8e-9 on network 1 (6-11% of the present
  error, including the reference noise) and at the reference's noise floor on network 0. The Gaussian-reference response
  geometry is adequate along this trajectory for the local maps. Whether it is adequate for the fourth-cumulant state map
  is not measured.

One cross-check ties the two theory threads together. The checkpoint's repeated-ReLU responses (its (33):
D Var = c/3, D kappa3 = v, D kappa4 = c + 4c^3 for the unary cubic tangent at N(0, 1)) are reproduced to rounding by the
independent quadrature of note XXXIV (`notes/gated-transport/code/check_gated_transport.py`).

## 4. The computation this audit makes possible

Everything in section 2 is deterministic in arrays the chain stores at each layer (K4v, K22, the Ritz pairs, W, C_y) and
was dumped in the one-step run of note XXXII. The cells' slices can be formed offline and compared with each other
(tensor identities, no truth needed: the known interaction -C(dM, N_E) is a bookkeeping test) and with the Monte Carlo
slices of z_(l+1) (the representation and missing-shape split). That is the checkpoint's section 10.2, with the
implementation chart now known.
