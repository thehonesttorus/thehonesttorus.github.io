# N4. Noncommutative Dirichlet forms and Gamma-calculus as the algebra of the layer telescoping

Frame N4 of the NCG round on the heat-defect principle. Status labels: **proved** (complete proof here), **sketch**
(argument given, routine details omitted), **checked** (verified numerically, output quoted), **measured** (Monte Carlo
at small width), **conjecture**.

Code and raw outputs are in `loc/n4/`:
- `chk1_pairdefect.py` checks Theorem B by exact finite differences (outputs `chk1_out.txt`, `chk1_out2.txt`);
- `chk2_closed.py` checks the closed-walk visibility of Theorem C (`chk2_out_n16.txt`, `chk2_out_n16L8.txt`);
- `chk3_gamma_op.py` checks the operator versus Schur carre du champ of Proposition 1.1 (`chk3_out.txt`);
- `chk4_edgeworth5.py` measures the size of the first omitted Edgeworth order (`chk4_out_n128.txt`, `chk4_out_n256.txt`);
- `chk5_hub_identities.py` checks the hub forms of the late-layer Gamma pairings (`chk5_out.txt`).

Notation follows F1:
- H = d_Sigma - (1/2) Hess_m;
- D = HE;
- delta = e - y.grad e - Lap e;
- delta(0) = 2 tr D(0, I) for a homogeneous chain.

One unit is 2n^3 FLOPs at n = 1024; B = 1024 units; the production bill is 208 units.

## 0. Summary

**The question.** F1 wrote the layer telescoping D(Phi(X)) = DPhi.D(X) + d(Phi; X) as an identity. Here it is read as
the diffusion (carre du champ) chain rule of a Dirichlet form, and that form is taken in its noncommutative
(Cipriani-Sauvageot) version on the matrix algebra in which the chain's tables live.

1. **The local defect is a Gamma-pairing (Theorem A, proved).**
   - In moment coordinates the local defect of any layer map is d(Phi) = -D^2 Phi[Gamma(M, M)].
   - Here Gamma(M_p, M_q) = (1/2) <J_1 xi_p, J_1 xi_q> is the Gram matrix of the first-chaos projections of the
     observables whose expectations the map reads.
   - So the defect is the Hessian of the closure in the mixture (moment) affine structure, contracted with the
     first-chaos Gram of its inputs.
   - The exact layer is linear in the law, so its defect is zero.

2. **Carried-pair formula (Theorem B, proved; checked to 1e-5 to 1e-2).**
   - Every readout of the chain has cumulant-exponential (Wick/Edgeworth) form Phi_S = exp(sum_{A in S} kappa_A
     d^A / A!) G over a set S of carried multi-indices. Its local defect is
     tr H Phi_S = sum_{A,B} [1_S(A+B) - 1_S(A) 1_S(B)] (d^{A+B} Phi_S / (A! B!)) Gamma(kappa_A, kappa_B).
   - The visibility weight of each omitted open class follows: N_S(k) = #{a : a, k-a in S}. It is 2 and 1 for GC
     (kappa_3 and kappa_4 paths), 3 for K3, and 4 and 3 (kappa_5 and kappa_6) for the production chain.
   - Omitted *tables* are weighted by tr R(eps)/eps instead. For the off-diagonal covariance this is its first-chaos
     share, which explains mean field's slope of about -1, deepening with L.
   - **Redundancy corollary.** The cheap (n^3) part of the analytic defect is exactly the next-order open walks at
     weights N_S. A single-slope defect merge is therefore a weight-mismatched version of carrying those walks at their
     derived Edgeworth weights.

3. **Closed walks: the Laplacian closes vector states into the trace (Theorem C, proved; checked).**
   - The second chaos is a noncommutative probability space (M_n, tau = tr, omega_ab = <L_a, . L_b>). Cumulants are
     open walks (vector-state moments) plus closed walks (trace moments).
   - Gamma maps open walks to open walks (concatenation). The m-Laplacian maps open walks to closed ones:
     Lap omega_ab(X) = 2 tau(H_a X H_b).
   - With homogeneity and the Euler-Stein rung Lap K_2 = -K_2, an omitted closed walk tau_k in a carried table is
     visible at weight w_k = k - (k/2) X_k / tau_k, where X_k is a third-chaos contraction.
   - Checked (n = 16): w_3 = 3.00, 2.70, 2.21, 1.76, 1.61, 1.81, 1.81 and w_4 = 4.00, 3.50, 2.71, 1.99, 1.72, 2.03,
     2.20 at hidden layers 1-7, with per-neuron correlation 0.94-1.00. So deep closed walks are seen at about weight 2,
     the weight of an open kappa_3 class.
   - This settles the open point of SYNTHESIS Prop 6. Closed omissions are visible to the *exact* defect, but only
     through Laplacians, and the Gamma-built local defects contain no Laplacians. So the *analytic* defect sees them only
     through full-rank Schur Gammas, which is the n^4 class.

4. **Factorization (Theorem D, proved; checked).** A Gamma-term costs n^{tw+1}, where tw is the treewidth of its
   contraction graph with the output indices joined.
   - The noncommutative carre du champ Gamma^op(C, C) = (1/2) sum_j (dC_j)^2 is treewidth 2: n^3, about 8 units per
     source-layer.
   - The commutative (Schur) one, Gamma^S(C, C)_ab = (1/2) sum_j (dC_ab,j)^2, is K4: n^4.
   - They are linked by the Cartan conditional expectation, E_D[Gamma^op] = rowsum Gamma^S. Every finite-Schur-rank
     pairing <u v^T, Gamma^S> is a twisted Gamma^op.
   - **So noncommutativity is exactly what makes the cheap part cheap.** The wall is the commutative (pointwise on
     pairs) carre du champ at full Schur rank.
   - For the chain:
     - every per-neuron readout defect (mean, variance, the D3 and g4 readouts) and every cross term
       Gamma(mu_a, C_ab) is n^3;
     - per-pair defects whose two indices both touch two internal vertices (the pair Gram, the (2,2) hierarchy source,
       the Mehler truncation) and every per-neuron closed walk are n^4;
     - the cross term that carried GC's defect is *free* in V56: Gamma(mu_a, C_ab) = D21_ab - Y_ab.

5. **Response principle (Proposition E, proved).**
   - The defect of a table is minus the heat operator of its error, D(X) = -H(eps_X). With homogeneity, a component
     matched at the base point contributes -(1/2) Lap_m eps.
   - Fitted closures and counterterms produce defect without error to the extent that their response differs from
     the law's.
   - For the production lambda-core this is estimated at about 10% of its hierarchy source (sketch): small.

6. **Curvature (Section 7).**
   - The ray-quotient flow has Bakry-Emery curvature -1/2 (CD(-1/2, infinity)), so it gives no contraction and no a
     priori sign.
   - The Gaussian readout's mixture Hessian is indefinite, so there is no Jensen sign.
   - What survives is positivity of the vector and trace states on squares: *even* walks have a sign, odd walks do not.
     Per-neuron sign rules therefore exist only for even classes; the kappa_4 path is the exhausted gain rule.

7. **The costed component and its prediction.**
   - **Cost.** The cheapest exact-to-first-order defect of the last layer is d_L = -2 g_5 Pi5 - (3/2) g_6 Pi6, where:
     - Pi5 = sum_s rowsum(P_s o w2_s o Z_s^2) and Z_s = (Y M^-1) At_s, one product per source: 4-8 units;
     - Pi6 is one more hub product plus a solve: 5-9 units.

     The last 2 layers cost about 20-40 units, the last 4 about 60-100.
   - **Value.** By the redundancy corollary its value is the kappa_5-level open classes. Monte Carlo gives T5/T4 at
     depth of 0.27-0.33 at n = 128 and 0.15-0.17 at n = 256. That is a width exponent of about -0.9, which
     extrapolates to 0.04-0.10 at n = 1024. It puts the kappa_5 class at about 1-3% of the final MSE per late layer
     and 2-6% overall.
   - **Verdict.**
     - The analytic late defect merge (A-late) is predicted negative: kappa_a X* <= 0.1 against a break-even of 0.33.
     - The class-resolved alternative, V57_P5 (the 5-path at its derived weight, layers 12-15, +16-32 units), is
       predicted at -1.5% to -4% raw: neutral to worse in adjusted MSE.
   - **Gate.** One cheap Monte Carlo oracle experiment decides it (section 10).

**Self-verdict: marginal.** Theory: the noncommutative structure is genuine and decides the cost classification, the
visibility weights and the closed-walk question. Estimator: the frame proves that the affordable part of the defect is
the next open walk, and that walk is worth little.

## 1. Two carres du champ on the chain's tables

### 1.1 Setting

- u(m, Sigma) = E F(m + Sigma^{1/2} Z) is the exact functional.
- A *table* is any quantity the chain carries at a layer, as a function of the input law (m, Sigma): means, the dense
  covariance, kappa_3 slices (D3, D21), the kappa_4 diagonal, source legs.
- **Derivation.** For an exact moment X(m, Sigma) = E_{N(m, Sigma)} xi(x), Stein's lemma gives
  d_j X = E[xi (Sigma^-1 (x - m))_j] = E[D_j xi], the first-chaos coefficient of xi. For cumulants this is the tilt
  identity d kappa_k(z_i1, .., z_ik) = kappa_{k+1}(z_i1, .., z_ik, X). Write d: tables -> tables (x) H_1, where
  H_1 = R^n is the input's first chaos, and call it the *first-chaos derivation*.
- **Carre du champ.** Gamma(X, Y) := (1/2) <dX, dY>_Sigma, the H_1-legs contracted with the metric Sigma.
- **Heat operator.** Traced, tr H = tr d_Sigma - (1/2) Lap_m. Since Lap_m(Phi o X) = DPhi.Lap X + D^2 Phi[grad X, grad X],
  the diffusion chain rule reads

      tr H (Phi o X) = DPhi . tr H X - D^2 Phi[Gamma(X, X)].                                                (1.1)

  This is F1 Prop 8 in carre-du-champ form; (1.1) is its only input.
- **Hierarchy.** For exact moments tr H M = 0. For exact cumulants (Hopf-Cole, F1 Thm 4),
  tr H kappa_k = sum over ordered (A, B), A u B = [k], of Gamma(kappa_|A|, kappa_|B|).

### 1.2 The noncommutative and the commutative carre du champ

Matrix-valued tables (C, D21, wk4m, wk431, and the legs) live in A = C_b^inf(R^n_m) (x) M_N(R), with N the number of
neurons. The localization semigroup acts as T_t = P_t (x) id, where P_t is the heat flow in m. Two Dirichlet
structures carry the same derivation:

- **(op)** A as a von Neumann algebra with trace E_m (x) tr. The derivation d: A -> A (x) R^n takes values in the
  Hilbert A-bimodule with left and right *matrix* multiplication. This is the Cipriani-Sauvageot picture (generator =
  d* d; Cipriani-Sauvageot 2003; Wirth arXiv:2207.09247; Vernooij-Wirth arXiv:2303.15949). Its A-valued carre du champ is
  Gamma^op(F, G) = (1/2) sum_j (d_j F)^* (d_j G).
- **(S)** A as the commutative algebra C_b^inf (x) l^inf([N]^2), with the pointwise (Schur) product on pairs. Its carre
  du champ is Gamma^S(F, G)_ab = (1/2) sum_j d_j F_ab d_j G_ab.

**Proposition 1.1 (proved; checked).** For symmetric F:
1. E_D[Gamma^op(F, F)] = diag(Gamma^op(F, F)) has entries sum_c Gamma^S(F, F)_ac. The conditional expectation onto the
   Cartan (diagonal) subalgebra D_N of the noncommutative carre du champ is the row marginal of the commutative one.
2. For K = sum_r u_r v_r^T (Schur rank R), <K, Gamma^S(F, F)> = (1/2) sum_r sum_j tr(D_{u_r} d_j F D_{v_r} d_j F). A
   Schur multiplier of finite rank turns the commutative pairing into a twisted noncommutative one (a D_N-bimodule
   map X -> D_u X D_v).
3. Gamma^op(F, F) >= 0 as a matrix, and Gamma^S(F, F) >= 0 entrywise.

*Proof.* (1) (sum_j d_jF d_jF)_aa = sum_j sum_c (d_j F_ac)^2. (2) tr(D_u X D_v X) = sum_ab u_a X_ab v_b X_ba with X
symmetric. (3) Sums of squares. QED.

**Check** (`chk3_out.txt`). The covariance tangent was CP-sourced, d_j C = sum_m w_m (A_m P_m^T + P_m A_m^T) l_mj:
- Gamma^op computed by the doubled-edge formula without forming the tensor: relative error 3.4e-16;
- item 1: relative error 1.8e-16;
- item 2: exact.

At n = n_births = n_in = 1024 the tensor route costs 1024 units. Gamma^op and a rank-1 Schur pairing cost about 8
units each per source-layer.

### 1.3 The second chaos as a noncommutative probability space

Truncate a layer's pre-activations at the second chaos (note XLI):

    z_a = mu_a + L_a . x + (1/2)(x^T H_a x - tr H_a),   H_a in M_n(R)_sa.

**The space.** Use the algebra M_n(R) with the trace tau = tr and the vector functionals
omega_ab(X) = <L_a, X L_b>. For k >= 2, the cumulant generating function gives

    kappa(z_a1, .., z_ak) = (1/2) sum_{s in S_k} omega_{a_s1 a_sk}(H_{a_s2} ... H_{a_s(k-1)})
                          + (1/(2k)) sum_{s in S_k} tau(H_{a_s1} ... H_{a_sk}).                              (1.2)

The first sum is the *open walks* (vector-state moments of the noncommuting tuple); the second is the *closed walks*
(trace moments). On the diagonal this is note XLI's kappa_k = k! [omega(H^{k-2})/2 + tau(H^k)/(2k)]. The joint
cumulants of the layer are the mixed moments of a noncommutative distribution; that is where noncommutativity lives in
this chain (SYNTHESIS 1.3.2).

**Proposition 1.2 (heat flow on the second chaos; proved).** Under the input shift m, L_a(m) = L_a + H_a m and H_a is
fixed. For m-independent X and Y:
1. **Laplacian closes.** Lap_m omega_ab(X) = 2 tau(H_a X H_b). The Laplacian maps vector-state moments to
   trace-state moments: it closes an open walk through the input.
2. **Gamma concatenates.** 2 Gamma(omega_ab(X), omega_cd(Y)) = omega_bd(X^T H_a H_c Y) + omega_bc(X^T H_a H_d Y^T)
   + omega_ad(X H_b H_c Y) + omega_ac(X H_b H_d Y^T). Gamma maps pairs of open walks to sums of open walks, joining
   them across one new edge. It never produces a trace.
3. **Positivity.** omega_aa(Z^* Z) >= 0 and tau(Z^* Z) >= 0. Every even walk |H^j L_a|^2 and tau(H^{2j}) is
   nonnegative; odd walks are unsigned (H_a is indefinite: 52% of its eigenvalues are negative, SYNTHESIS W4).

*Proof.* (1) d_j d_j of (L + Hm)^T X (L + Hm). (2) d_j omega_ab(X) = (H_a e_j)^T X L_b + L_a^T X H_b e_j; take
products and sum over j. (3) Positivity of the two functionals. QED.

**Reading.** Local defects are built from Gammas only (1.1). So a chain whose tables are open walks has open-walk
local defects. Closed walks reach the defect only where a *Laplacian* acts: in the exact defect of a carried table
(Theorem C), or, inside the telescoping, through Gammas of pair tables whose CP contraction graph closes a cycle
(Theorem D).

## 2. Theorem A: the local defect is a mixture-Hessian against the first-chaos Gram

**Theorem A (proved).** Let a layer map be written in moment coordinates, N = Phi(M). Here M are the moments of the
layer's pre-activation law that the closure reads, and N are the moments of the post-activation law it outputs.

1. **Defect.** If the inputs are exact (tr H M = 0), then

       tr H N^chain = -D^2 Phi[Gamma(M, M)],   Gamma(M_p, M_q) = (1/2) <E[D xi_p], E[D xi_q]> = (1/2) <J_1 xi_p, J_1 xi_q>,

   where M_p = E xi_p. This is the local defect.
2. **Zero defect.** The local defect vanishes for all exact inputs iff Phi is affine along every direction spanned by
   the first-chaos tangents of the attainable laws. The exact layer (E f(z)^k is linear in the law) is one such map.
3. **Cumulant coordinates.** Changing to cumulant coordinates kappa = K(M) gives F1 Prop 8's bracket,
   d = DPhi.R(X) - D^2 Phi[Gamma(X, X)] - R(Phi(X)), because tr H K(M) = -D^2 K[Gamma(M, M)] is exactly the
   Hopf-Cole source R.

*Proof.* (1.1) with tr H M = 0. Item 2: D^2 Phi vanishes on the span of the Gram's ranges iff Phi is affine there.
Item 3: apply (1.1) to Phi o K^-1 and K o Phi. QED.

- **The noncommutative reading.** Per-neuron tables are elements of the Cartan subalgebra D_N. Pair tables are
  elements of M_N viewed as functions on the pair groupoid. The derivation d sends both into (.) (x) H_1, and Gamma
  contracts the H_1-legs. In CP form, a d-image is a source term with one leg replaced by its first-chaos leg: for a
  birth at m this is l_m, the birth's first chaos, which the Schur identity dresses through the current field
  (section 5.3).
- **Why it matters.**
  - Every local defect is a *linear functional of one Gram matrix*, Gamma(M, M). The coefficients are fixed by the
    closure's second derivatives.
  - So the cost question becomes which entries of that Gram a given closure needs, and how each entry contracts. That
    is Theorem D.

## 3. Theorem B: the carried-pair formula for Wick/Edgeworth readouts

**The readouts.** Every readout of the production chain is a truncated cumulant-exponential (Wick/Edgeworth)
functional (`estimator_final_v56.py`, TERM_SPECS).
- The mean readout (1,) is g + (D3/6) g_3 + (g4/24) g_4 + (D3^2/72) g_6 + (D3 g4/144) g_7, where g_k = d_mu^k E relu.
  This is exp(D3 d^3/6 + g4 d^4/24) g truncated at second order without the g4^2 term.
- The moment readouts (2,), (3,), (4,) are linear in D3 and g4.
- The pair programs (1,1), (2,1), (2,2) read c_off, d21, d21T, d3row/col, g4row/col, wk4m and wk431 at Mehler order 2,
  with products of up to two tables.

**Theorem B (carried-pair formula; proved).** The setting:
- Let G(mu, C) be a Gaussian expectation over a set of neurons. It is heat-consistent: d_{C_ab} G = d_a d_b G for
  a != b, and d_{v_a} G = (1/2) d_a^2 G (Price).
- Let S be a set of multi-indices containing every index of order 1 and 2.
- Let Phi_S = exp(sum_{A in S, |A| >= 3} kappa_A d^A / A!) G, with A! = prod_i A_i!.

If the carried cumulants are exact (tr H kappa_A = tr R(kappa_A) for A in S), then

    tr H Phi_S = sum over ordered (A, B), A, B != 0, of [1_S(A + B) - 1_S(A) 1_S(B)] (d^{A+B} Phi_S / (A! B!)) Gamma(kappa_A, kappa_B).   (3.1)

*Proof.*
1. For the exponential form, d_{kappa_A} Phi_S = d^A Phi_S / A! for every A in S. For |A| <= 2 this uses heat
   consistency of G.
2. Hence d_{kappa_A} d_{kappa_B} Phi_S = d^{A+B} Phi_S / (A! B!).
3. By (1.1), tr H Phi_S = sum_{C in S} (d^C Phi_S / C!) tr H kappa_C - sum_{A,B in S} (d^{A+B} Phi_S / (A! B!))
   Gamma(kappa_A, kappa_B).
4. The hierarchy for joint cumulants gives tr H kappa_C = sum over ordered (A, B) with A + B = C of
   binom(C, A) Gamma(kappa_A, kappa_B), and binom(C, A)/C! = 1/(A! B!).
5. Collect the coefficient of each Gamma(kappa_A, kappa_B). QED.

- **Truncated readouts.** The truncated forms differ from the exponential at order kappa_4 Gamma (the missing g4^2
  term of the mean) and kappa_3 Gamma (the missing D3^2 term of the variance readout). Both are at the kappa_6 level.
- **Check** (`chk1_out*.txt`). The exact chi-square model above with n = 6, a full cumulant tower and exact (m, Sigma)
  dependence; tr D by finite differences against (3.1):
  - relative agreement 3e-7 to 4.6e-3 for S = {1,2}, {1,2,3}, {1,2,3,4} and the gapped {1,2,4} (seeds 0, 2, 3 at
    second-chaos scale 0.04; seed 1 at 0.08);
  - the residual is the truncation of the exponential series. At scale 0.08 with seed 0 the Edgeworth series itself
    diverges (err_7 > err_5) and the agreement degrades to 3e-1, as it must.

### 3.1 Corollary B1: visibility weights of omitted open classes

In the second-chaos model, d mu = L, d kappa_a = a! H^{a-1} L and Gamma(kappa_a, kappa_b) = (a! b!/2) omega(H^{a+b-2})
on the diagonal. The first-order error of Phi_S is sum_{k not in S} (d^k Phi_S / k!) kappa_k, whose open part is
err_k^open = (1/2) d^k Phi_S omega(H^{k-2}). (3.1) then gives, for contiguous S = {1, .., K},

    -tr H Phi_S = sum_{k > K} N_S(k) err_k^open,     N_S(k) = #{a in [1, k-1] : a in S and k - a in S}.      (3.2)

In the readout's own local defect, closed walks have weight 0.

- **Check** (`chk1_out2.txt`). -sum_k N_S(k) err_k^open equals the finite-difference tr D to the same precision for
  every contiguous S. For example, seed 2, S = {1,2,3}: +1.81214e-3 against +1.81215e-3.
- **Detector convention.** In F1's convention err = -delta/(2p) with delta = 2 tr D, so a class seen at weight w has
  p = w and slope -1/(2w).

**Multivariate pair maps.** A multi-index class A of order k carries open walks on k vertices. Each such walk is
reached by cutting one of its k - 1 edges (Proposition 1.2.2). So tr R(kappa_A^open) = (k-1) kappa_A^open (sketch: an
edge-counting argument, exact on the diagonal by (3.2)). The weight of an omitted pair class is k - 1 times the
fraction of its edge cuts landing in carried x carried splits.

**Weights of omitted tables.** When a component eps of a *carried table* is missing, the readout reads the table
exactly, and the table's defect is D(X) = -H(eps) (Proposition E). Its weight is tr H(eps)/eps:
- an omitted exact covariance: tr R(C_ab)/C_ab = (L_a . L_b)/C_ab, the first-chaos share;
- an omitted open part of kappa_k: (k - 1) - (k!/2) tau(H^k)/kappa_k^open, so the closed partner lowers it (a
  homogeneity computation like Theorem C's; sketch);
- an omitted closed part: Theorem C.

| chain (F1/F5 measurements) | leading omitted classes and predicted weights w | predicted p = w | measured p (F1) |
|---|---|---|---|
| MF (diagonal covariance) | off-diagonal covariance in the transport: w = first-chaos share, about 0.4-0.7 at L = 8, falling with depth | 0.4-0.7 (L = 8), smaller at L = 16 | 0.46-0.61 (L = 8), 0.32-0.47 (L = 16) |
| GC (dense covariance) | kappa_3 open, readout and (2,1) map: 2; kappa_4 path: 1; gain kappa_3 half: 1; closed parts and gain kappa_4 half: 0 (F5 referee) | between 1 and 2, nearer 1 if the gain dominates | 1.1-1.2 (n <= 256) |
| K3 (S = {1,2,3}) | kappa_4 open: N = 3; kappa_5: 2; omitted GC feedback and kappa_4 feeds: lower | <= 3 | 2.3-2.8 |
| production (S = all orders <= 4) | readouts: kappa_5 level 4, kappa_6 level 3. Tables: closed walks w_k (Theorem C, about 2-3 at depth); compressions about 2; Mehler truncation about 1; closures and counterterms: response mismatch (Proposition E) | 2-3.5 | pending |

- **What the table supports.** The MF row is new quantitatively. F1's kink heuristic gave MF p = 1/2 by a different
  argument and failed at K3 (F1 referee R2.3). Here the weights are derived, and they order all three measured chains
  correctly:
  - the MF L = 16 slopes steepen (-1.06 to -1.54 against -0.82 to -1.08) because the first-chaos share falls with
    depth;
  - K3 sits below its leading weight 3 because it also omits lower-weight classes.

### 3.2 Corollary B2: each production readout, its carried set and its leading defect

| readout | carried set S | leading local-defect pairs (A, B), A + B not in S | Gamma form (bare second chaos; dressed in 5.3) | cost |
|---|---|---|---|---|
| mean (1,), layers 3-15 | {1,2,3,4} | (1,4), (2,3) [order 5]; (2,4), (3,3) [order 6] | -(g_5/24) Gamma~(mu, g4) - (g_5/12) Gamma~(v, D3) = -2 g_5 Pi5; order 6: -(3/2) g_6 Pi6 (+ the D3^2-truncation term) | n^3 per source |
| second moment (2,) | {1,2,3,4} linear | same pairs, with the g^(2) derivatives | reuses Pi5, Pi6 | about 0 extra |
| (3,), (4,) (post-activation D3, g4 reads) | {1,2,3,4} | same | reuses Pi5, Pi6 | about 0 extra |
| pair maps (1,1), (2,1), (2,2) | all multi-indices of order <= 4 (wk4m, wk431 as closures) | order 5: (4,1), (3,2) classes; the (1,1) + (2,1) and (2,1) + (1,2) Gammas are pair x pair | tree terms (per-neuron x pair) at n^3; pair x pair at n^4 | n^3 + n^4 |
| D3 from sources | open kappa_3 | omitted triangle tau(H_a^3) | closed | n^4 |
| D21 from sources | open (2,1) | omitted mixed triangle tau(H_a^2 H_b) | closed | n^4 |
| g4 = transported + lambda core + V56 path + star | open path + third chaos + closure | omitted 4-cycle 3 tau(H^4), kappa_3 x C product, closure response | closed + open | n^4 + n^3 |

Notation: Gamma~ := <d., d.> = 2 Gamma; Pi5 = omega(H^3) = L^T H^3 L; Pi6 = |H^2 L|^2.

- **The cross term is gone.** For the production chain the GC cross term Gamma(mu_a, C_ab) cancels in (3.1), because
  (1,0) + (1,1) = (2,1) is carried (D21).
- **So the production readouts' local defects start one order up**, at total cumulant order 5 with weight 4. That is
  the localization-visible part of the next order (F1 9), now with exact weights.

### 3.3 Corollary B3: redundancy of the analytic defect

**Corollary B3 (proved in the class model).** For a chain whose carried tables are open walks with law-correct
response, the n^3 part of the analytic defect is

    delta^an = -2 sum_classes N_S(c) err_c^open + (pair-Gram terms).

1. **What the merge becomes.** The defect merge e + delta^an/(2w) with class-resolved w adds sum_c err_c^open. That is
   exactly *carrying the next-order open walks at their derived Edgeworth weights*. The one-coefficient merge
   e - c delta is a weight-mismatched version of the same thing.
2. **What the analytic defect cannot add.** It brings no information beyond the next open walk, except the automatic
   composition through the chain (the S_l transport, which an in-place injection reproduces at first order).
3. **Where the rest lives.** Everything else a positive X* run might reveal is closed walks, compressions or
   calibration: the n^4 or finite-difference-only part.

*Proof.* (3.2) per readout plus (1.1) through the telescoping. QED.

## 4. Theorem C: closed walks become visible through the Laplacian

**Theorem C (proved).** Let z be a pre-activation of the homogeneous network, with chaos kernels K_2 = E D^2 z and
K_3 = E D^3 z at the base law, and let tau_k = tr K_2^k. As a function of the input law (m, Sigma) through the
y-kernel Sigma^{1/2} E_{N(m, Sigma)}[D^2 z] Sigma^{1/2}, it satisfies, at (0, I),

    tr H tau_k = k tau_k - (k/2) X_k,     X_k = sum_j sum_{p+q=k-2} tr(K_{3,j} K_2^p K_{3,j} K_2^q).        (4.1)

*Proof.*
1. **Homogeneity.** Under (m, Sigma) -> (cm, c^2 Sigma) the y-kernel scales by c, since D^2 z is homogeneous of
   degree -1. So tau_k has degree k, and at m = 0 Euler gives tr d_Sigma tau_k = (k/2) tau_k.
2. **Laplacian.** d_j K_2 = K_{3,j} and sum_j d_j^2 K_2 = tr_{34} c_4 = -K_2 (Euler-Stein ladder, rung 2, F1 Thm 3).
   Expanding Lap tr K_2^k by Leibniz gives Lap tau_k = -k tau_k + k X_k.
3. Combine: tr H = tr d_Sigma - (1/2) Lap. QED.

**Corollary C1 (visibility of closed omissions).** Suppose a carried table equals its law value minus a closed walk
c tau_k. Then D(table) = -c tr H tau_k, the error is +c tau_k, and the defect sees the omission at weight
w_k = k - (k/2) X_k/tau_k. At the first hidden layer the layer-1 means vanish, so K_3 = 0 there and w_k = k exactly.

**Check** (`chk2_out_n16.txt`; n = 16, L = 4, 4e5 antithetic samples; K_2 and K_3 by Stein on the exact input gradient,
K_2 = E[J (x) x] and K_3 = E[J (x) He_2(x)]):

| hidden layer | w_3, L = 4 net (per-neuron corr) | w_4, L = 4 (corr) | w_3, L = 8 net (corr) | w_4, L = 8 (corr) |
|---|---|---|---|---|
| 1 | 3.00 (1.00) | 4.00 (1.00) | 3.00 (1.00) | 4.00 (1.00) |
| 2 | 2.49 (1.00) | 3.18 (1.00) | 2.70 (1.00) | 3.50 (1.00) |
| 3 | 2.18 (0.97) | 2.69 (0.95) | 2.21 (1.00) | 2.71 (1.00) |
| 4 | | | 1.76 (0.98) | 1.99 (0.94) |
| 5 | | | 1.61 (0.99) | 1.72 (0.99) |
| 6 | | | 1.81 (1.00) | 2.03 (0.99) |
| 7 | | | 1.81 (0.99) | 2.20 (0.99) |

- **The identity itself.** Checked by exact finite differences on the 2-layer net (closed-form K_2(m, Sigma)):
  tr d_Sigma tau/tau = 1.500000 and 2.000000, and Lap tau/tau = -2.999999 and -3.999999, for k = 3 and 4.
- **The depth trend.** The third-chaos correction X_k grows over the first 3-4 hidden layers and then saturates.
  w_3 settles at about 1.6-1.8 and w_4 at about 1.7-2.2. The per-neuron shape stays (corr 0.94-1.00).
  Deep closed walks are therefore a weight-2 class, like the open kappa_3 class. Their width dependence is not
  measured here (n = 16 only).

**What Theorem C settles.**
- **The open point.** SYNTHESIS Prop 6 and Corollary 7 left open whether closed-walk omissions are visible in a
  composed chain; F5's referee flagged the inconsistency between F5 Prop 5.4 and F5 5.5c.
- **The answer.** Closed omissions are visible to the *exact* defect, at weight k minus a third-chaos term. The term
  grows with depth and saturates: the triangle goes 3, 2.7, 2.2, then about 1.7. The per-neuron shape is preserved
  (corr >= 0.94), so a single-slope detector catches them.
- **Why they are not visible cheaply.**
  - In the telescoping (1.1) only Gammas appear, and by Proposition 1.2.2 Gammas of open tables are open.
  - The closed content must therefore come from Gammas of tables whose CP contraction graph closes a cycle: the pair
    Gram of the transport and birth maps, the gate covariance and the triangle at birth.
  - Those are exactly the n^4 class (Theorem D).
- **Consequence for the gating run.** The exact defect from finite differences contains the closed classes; the
  affordable analytic defect does not. If X* is high because of closed walks, kappa_a of any n^3 analytic delta is low.

## 5. Theorem D: factorization, and where noncommutativity makes things cheap

### 5.1 The cost rule

**Theorem D (proved, given the CP source representation).** Write each tangent of a source-built table in CP form:
legs A_s, P_s (n x n) and first-chaos legs l_s (n_b x n_in, dressed in 5.3). Build a graph whose vertices are the
neuron, birth and input indices and whose edges are the leg factors. Join the output indices of the Gamma-term into a
clique. Then:
1. **Cost.** The Gamma-term costs O(n^{tw+1}) per source-layer pair, with tw the treewidth of that graph after merging
   parallel edges. Merging is free: Hadamard products of n x n matrices.
2. **Treewidth 2 (n^3).**
   - per-neuron open walks, omega_aa(H_a^k) = L_a^T H_a^k L_a: a fan of triangles sharing the output vertex a (for
     k = 2, a triangle with doubled edges, diag(X G X^T));
   - per-pair walks in which one of a, b touches a single internal vertex: (V G V^T)_ab, and the cross terms
     Gamma(mu_a, C_ab), Gamma(v_a, C_ab), Gamma(D3_a, C_ab);
   - Gamma^op(C, C) and its diagonal = rowsum Gamma^S;
   - every finite-Schur-rank pairing <u v^T, Gamma^S>.
3. **K4 (n^4).**
   - per-neuron closed walks tau(H_a^k) for k >= 3;
   - the per-pair Schur Gram Gamma^S(C, C)_ab;
   - per-pair walks with both a and b on two internal vertices, e.g. L_b^T H_a^2 L_b;
   - the (2,2) hierarchy source R(kappa_22) and the Mehler-truncation defect, which contain Gamma^S(C, C).

*Proof.* Item 1 is the standard variable-elimination bound for tensor-network contraction. It is the cost of the best
contraction order. It is not a lower bound against algebraic shortcuts such as low-rank structure of G or fast matrix
multiplication. The cases:
- **Open walk.** omega_aa(H_a^2) = sum_{m,m'} X_am G_mm' X_am' with X = P o At o w. The vertices are a, m, m' with
  edges am (doubled), am' (doubled) and mm': a triangle, eliminated as diag(X G X^T).
- **Closed walk.** tau(H_a^3) has edges am, am', am'', mm', m'm'', m''m: that is K4.
- **Gamma^op.** The doubled-edge formula of `chk3`: m and m' are joined by two parallel edges (via the contracted
  neuron c and the input j), so the graph reduces to a path.
- **Schur pair Gram.** Edges am, am', bm, bm', mm' plus the output clique ab: K4.
- **Finite Schur rank.** The pairing sums over a and b, which deletes the output edge ab and leaves the diamond K4 - e
  (treewidth 2). QED.

**The noncommutative statement.** By Proposition 1.1, the n^3 objects are exactly the functionals of the operator-
valued (noncommutative) carre du champ and of its Cartan marginal. The n^4 wall is the commutative carre du champ on
pairs at full Schur rank. The matrix (convolution) product of the pair groupoid is what makes the cheap part cheap. The
pointwise Schur product on pairs is the wall.

This is the precise form of SYNTHESIS 1.3.2: open walks are vector-state moments, reached at n^3; closed walks are
trace moments, n^4. The aggregate (Cartan-projected, annealed) part of every closed object is cheap; only its flat,
per-pair or per-neuron resolution is not.

### 5.2 What this says about the chain's readouts

- **The final layer.** It emits only means, so all of its local Gammas are per-neuron (Table 3.2): **n^3**. Its table
  defects (D3_L and g4_L: the triangle and the 4-cycle) are **n^4**.
- **Layer L-1.**
  - The covariance map's tree terms are n^3.
  - Its pair x pair terms are n^4. The layer-L variance reads them through diag(W^T d^C W), which for each output
    neuron is a Schur-rank-1 pairing (n^3 per neuron, n^4 over all neurons).
  - Aggregated over output neurons with fixed weights they are n^3. That allows *per-layer scalar* (calibration-type)
    use, but not per-neuron use.
- **Truth-free calibration is affordable only in aggregate.**
  - A per-layer scalar such as a closure amplitude can be fitted to zero the Cartan-projected (row-summed) hierarchy
    violation at n^3 per source-layer, or by Hutchinson probes at negligible variance, since it averages over n^2
    pairs.
  - The quenched per-neuron resolution, which is where DATA_FINDINGS puts the remaining error (signed mean error
    < 1e-5 per layer), stays n^4.
  - The production lambda is already adaptive per network (V25, 1-2% unexplained spread), and the counterterm ceiling
    is -3.8%. **No calibration gain is predicted.**

### 5.3 Dressing: the Schur hub as the field-localization carre du champ, and the free cross term

**Proposition 5.1 (proved).** Let At = L l^T (arms), Y = sum_s (P_s o At_s o w2_s) At_s^T and G = L L^T (the
first-chaos Gram, assumed invertible). In the second-chaos model:
1. Y_ij = L_j . v_i with v_i = H_i L_i.
2. Pi5_i = omega_ii(H_i^3) = sum_s sum_m (P_s o w2_s)_im (Z_s)_im^2, with Z_s = Y G^-1 At_s.
3. Pi6_i = |H_i^2 L_i|^2 = [Yt G^-1 Yt^T]_ii, with Yt = sum_s (P_s o w2_s o Z_s) At_s^T.
4. **Cross term.** Gamma~(mu_a, C_ab) = Y_ab + Q_ab with Q = sum_s (At_s o At_s o w2_s) P_s^T. Since the open (2,1)
   slice is 2Y + Q, **the cross term equals D21 - Y**: free in a chain that forms both hub halves (V56).
5. **Mixed walks** for the (3,1)-closure defect: v_a^T H_b L_a = [(Z o At) d(w2) P^T]_ab and
   v_a^T H_a L_b = [(P o Z) d(w2) At^T]_ab.

*Proof.* The births' first-chaos vectors are l_m = L^T G^-1 At_.m (L invertible); substitute. QED.

**Check** (`chk5_out.txt`): all five identities hold to machine precision.

- **Dressing.** The chain replaces G by its own covariance M = C + eps mean(var) and the first-chaos arms by its full
  arms (note XLI).
- **Gamma reading.** G^-1 replaced by C^-1 is the carre du champ of the *whitened field localization* of the layer law:
  a linear-tilt localization along z_l with control C C^T = C_l^-1 dt, under which F1's Theorem 5 holds for any law.
  - At layer 1 this is the isotropic input localization exactly: z = W^T x, so a tilt along z with metric C^-1 =
    (W^T W)^-1 is the input tilt.
  - At depth it projects on span(z_l) instead of on the input's first chaos.
- **Why note XLI found the dressed hub at weight 1 and the bare one 4-8 times too small.** The dressed hub is the
  carre du champ of a localization that sees the layer's own field. The first-chaos version sees only the 21.5% linear
  share (DATA_FINDINGS).
- **Commutative Petz metrics.** In the commutative case all Petz metrics coincide, and C^-1 is stage-8 Theorem C's
  capacity (regression) pairing. No modular structure enters here; F4 already showed the KMS-mean gauge is not a lever.

## 6. Proposition E: the defect sees response, not value

**Proposition E (proved).**
1. **Defect = minus the heat operator of the error.** Let X^chain be a table and eps_X = kappa_true - X^chain, with
   R computed from exact lower tables. Then D(X^chain) = H X^chain - R(kappa) = -H(eps_X), exactly.
2. **With homogeneity of degree d**, tr H eps = (d/2) eps - (1/2) Lap_m eps at the base. A component that is exact at
   the base point (eps(0, I) = 0), such as a closure fitted there, still contributes -(1/2) Lap_m eps(0, I). Its
   integrated contribution along the localization vanishes, because
   eps(0, I) = E int (1+t)^-2 tr H eps(m_t, Sigma_t) dt = 0 (F1 Theorem 1 for eps, with eps = 0 at point masses).

So the base-point detector reports an error that is not there: the closure's profile b(s) is non-monotone, with
b(1) = 0 and b'(1) != 0.

*Proof.* Item 1: H kappa_true = R(kappa) by the hierarchy. Item 2: the Euler relation at m = 0, then F1 Theorem 1
applied to eps. QED.

**Size for the production closures (sketch).**
- **The lambda core**, X = lambda_l C_ab, with lambda adaptive through the aggregate ratio mean(dG)/mean(var) (V25,
  beta = 1).
  - Homogeneity fixes the Sigma-part of the response (degree 4) for any homogeneous closure. Only the m-Laplacian
    mismatch remains.
  - With Lap C_ab = 2(C_ab - G_ab) (the radial identity), and Lap(aggregate degree-d cumulant)/aggregate =
    d - 2 tr R/aggregate, one gets tr H(lambda C) of about 2.75 lambda C, against tr R(kappa_4-slice) of about
    3 kappa^open.
  - The false defect is therefore about 10% of the slice's hierarchy source.
- **Counterterms** (1 + d)X add d tr H X, with |d| of 0.01-0.17: under 0.5 of the rms error.

**Prediction.** The response mismatch lowers the gating run's X* by a factor of only 0.9-0.97. It is not a dominant
effect. It does explain the F1 referee's finding (R2.4) that one counterterm moves the slope: the counterterm's
response is law-like while its correction is not.

## 7. Curvature (frame question iv)

1. **The m-heat flow.** Its Bakry-Emery operator is Gamma_2(f) = (1/4)|Hess f|^2_HS, so it is CD(0, n).
2. **The ray-quotient flow A = (1/2)(Lap + y.grad).** Gamma_2^A(f) = (1/4)(|Hess f|^2 - |grad f|^2), so it is
   CD(-1/2, infinity): repulsive and negatively curved. No gradient contraction is available.
   - The truth is a (A - 1/2)-eigenfunction (F1 Cor 2.3).
   - Curvature bounds therefore give neither a sign nor a size for delta beyond Cauchy-Schwarz on Gamma.
3. **No Jensen sign.** In mixture coordinates the Gaussian readout's Hessian is the quadratic form
   g_vv c^2 + 2 g_muv a c, where a and c are the tangents of mu and v. Its discriminant g_muv^2 is >= 0, so the form
   is indefinite. A defect that is a Hessian against a PSD Gram therefore has no sign. The same holds for every
   Edgeworth readout, since g_k changes sign with He_{k-2}(alpha).
4. **What survives: positivity on squares (Proposition 1.2.3).** A class has a per-neuron sign iff its walk is even
   and its Edgeworth coefficient has a known sign:
   - the kappa_4 path class (|HL|^2 >= 0) has sign(1 - alpha^2). This is the gain/kappa_4 sign rule, measured exhausted
     on the production residual (SYNTHESIS W4);
   - the kappa_6 path Pi6 (>= 0) has sign He_4(alpha) = alpha^4 - 6 alpha^2 + 3;
   - odd classes (D3, Pi5, the triangle) are unsigned, because H_a is indefinite.
5. **The only useful "curvature" is the per-layer contraction of the defect transport.** The mean-gate Lipschitz
   constant per layer is 2 E[Phi(alpha)^2]. DATA_FINDINGS' 9x loss between layers 3 and 15 corresponds to 0.83 per
   layer, i.e. E[Phi(alpha)^2] of about 0.41. This sets which layers' local defects matter: the last 4 carry 65%,
   consistent with the injection profile.

**Answer to (iv).** No CD(K, N)-type bound gives the sign or the size of the per-layer defect. Sign rules exist
exactly for even walks, by positivity of the two states, and the even classes the production chain still omits
(Pi6, the 4-cycle) are small (section 8).

## 8. The runtime component: the cheapest exact-to-first-order late defect, costed

### 8.1 The final layer (K = 1)

The final layer emits only the mean. By Table 3.2 its local defect at orders 5 and 6 is

    d_L = -2 g_5 Pi5 - (3/2) g_6 Pi6 + (the kappa_6-level terms from truncating the exponential at second order),

with Pi5 and Pi6 in the dressed hub form (Proposition 5.1, G -> M). The flopscope cost at the last layer:
- The V56 hub product Y and the solve M^-1 Y^T are formed at layer 15 already (V56_LAST = 2), so B := Y M^-1 is the
  transpose of V56's solve output: free. At most 1 unit if the production path keeps only column norms.
- Z_s = B At_s for the young dense sources (age <= 4, about 5): 1 unit each, about 0.5 with the chain's Strassen
  levels.
- The old tiers in factor space:
  - tier 1 (r = 320, about 4 sources): (B Qc) once at 0.31 units, then (B Qc) FAt_s at 0.31 units each;
  - tier 2 (r = 192, about 7 sources): B QU once, then 0.19 units each.
- Pi5 = sum_s rowsum(P_s o w2_s o Z_s^2): about 16 n^2 elementwise operations, about 0.03 units.
- **Pi5 total: 4-8 units.**
- Pi6: Yt = sum_s (P_s o w2_s o Z_s) At_s^T is the same product count again, plus one solve on the existing
  factorization and a row norm. **Pi6: +5-9 units.**

### 8.2 The last 2-4 layers

- **Each earlier layer l in {12, 13, 14}.**
  - Its mean and variance readouts need Pi5 and Pi6 at that layer (Y and the factor exist for l >= 3 under V56):
    4-8 units, plus 5-9 for Pi6.
  - The covariance map's tree terms Gamma~(v_a, C_ab), Gamma~(D3_a, C_ab) and the mixed walks of Proposition 5.1(5)
    cost two products per source: about 8-16 units.
  - The pair x pair terms are n^4 and excluded.
- **Transport to the output.** The in-place injection (SYNTHESIS A-late) needs no separate forward pass.
- **Totals.** K = 1: 4-17 units; K = 2: 20-40; K = 4: 60-100. That is +2% to +50% of the bill.

### 8.3 What it is worth: the redundancy corollary and the Monte Carlo size of the kappa_5 class

- **What the cheap defect contains.** By Corollary B3, the n^3 late defect is the kappa_5- and kappa_6-level open
  classes at weights 4 and 3. Its maximal value is what carrying those classes exactly would remove.
- **Their size** (`chk4_out_n128.txt`, `chk4_out_n256.txt`; random MLP, L = 12, two halves of 2^18 samples, signal
  rms by split-half covariance). T_k = g_k kappa_k / k! is the k-th Edgeworth term of the mean readout.

| layer | n = 128: rms T4 | T5/T4 | T6/T4 | r | n = 256: rms T4 | T5/T4 | T6/T4 | r |
|---|---|---|---|---|---|---|---|---|
| 2 | 1.57e-3 | 0.12 | 0.02 | 0.071 | 6.84e-4 | 0.065 | 0.035 | 0.033 |
| 5 | 2.05e-3 | 0.24 | 0.15 | 0.141 | 8.21e-4 | 0.111 | 0.055 | 0.067 |
| 8 | 1.76e-3 | 0.27 | 0.13 | 0.187 | 8.13e-4 | 0.157 | 0.064 | 0.088 |
| 11 | 1.81e-3 | 0.33 | 0.30 | 0.221 | 9.23e-4 | 0.168 | 0.100 | 0.114 |

(r is the median kappa_4/(2 var^2).)

- **Width scaling, measured.**
  - r and rms T4 halve per doubling, as SYNTHESIS 1.3.3 says.
  - T5/T4 falls by 0.5-0.6 per doubling at depth, a width exponent of about -0.9. That is faster than the naive
    open-walk scaling n^{-1/2}, so the open 5-path is probably not the dominant part of kappa_5.
  - T6/T4 falls by about 0.4-0.5 per doubling.
- **Extrapolated to n = 1024, layers 12-15** (depth raises the ratio by about 1.2 from layer 11 to 15): T5/T4 is
  0.04-0.10 and T6/T4 is 0.01-0.03.
- **From size to MSE.**
  - T4 scales like r sigma: about 9e-4 at n = 256, hence about 2.3e-4 at n = 1024. That is about 1.8 times the
    production's final rms error, 1.25e-4. This agrees with the independent estimate from the g4 oracle (15-30% of
    MSE at a 30% table error).
  - Hence rms T5 is 1e-5 to 2.3e-5, and (T5/err)^2 is 0.6-3.4% of the final MSE per late layer's readout. Summed over
    the late layers with transport, the kappa_5-level classes hold 2-6% of the MSE.
  - The open part (what Pi5 carries) is only a fraction of that: the closed 12 tau(H^5) and third-chaos parts are the
    rest.

**Predictions (costed).**

| component | added cost | predicted raw | predicted adjusted (from 3.14e-9) |
|---|---|---|---|
| A-late defect merge, analytic open part, K = 1-4, single slope | +4 to +100 units | -0.5% to -4% (kappa_a X* about 0.02-0.1; weight mismatch costs a further 10-20% of that) | 3.1e-9 to 4.6e-9: neutral to worse |
| **V57_P5**: the 5-path at its derived weight, (g_5/2) Pi5 in the mean and variance readouts of layers 12-15 (class-resolved; no slope) | +16 to +32 units (+8% to +15%) | -1.5% to -4% | 3.2e-9 to 3.5e-9: worse |
| V57_P5 + Pi6, layers 12-15 | +36 to +68 units | -2% to -4.5% | worse |
| V57_P5 at the final layer only | +4 to +8 units | -0.5% to -2% | 3.12e-9 to 3.24e-9: neutral |

Break-even for a +8% to +15% cost is -7.4% to -13% raw. **The frame's runtime component is predicted not to pay**
unless the oracle experiment below finds the kappa_5 class at least 2-3 times as large as extrapolated.

## 9. Predictions for the pending X* run (README section 4, registered X* = 0.45 [0.25, 0.70], slope -0.15 to -0.25)

Each is stated before the run.

1. **Slope.** The visible classes have weights about 1 (Mehler and covariance-type table errors), 2-3 (closed walks
   at depth, compressions) and 3-4 (readout truncations). Hence p_eff = sum w^2 E / sum w E in [2, 3.5] and a slope
   of -0.14 to -0.25.
2. **X*.** The weight spread alone costs about 10% (Corollary 7 of SYNTHESIS with weights {1, 2, 3, 4}: 0.83). The
   response mismatch costs a factor of 0.9-0.97. Invisible components (third-chaos-cancelled closed walks,
   calibration) cost the rest. **X* = 0.35 [0.2, 0.6].**
3. **rms delta / rms err is about 2 p_eff sqrt(X*): 2-5** (OLS: |a| = sqrt(X*) |err|/|delta|, with a = -1/(2 p_eff)).
   F1's K3 value was 4.2-4.9.
4. **Per-layer decomposition (from the exact 2049-run defect on network 0).**
   - The final layer's own local defect, the d_L of section 8.1, captures at most 15% of delta's variance at layer 15.
   - The closed-walk part (triangle and 4-cycle of the late tables) captures 20-50%.
   - **Falsifier.** If the n^3 analytic late defect computed from the chain's own tables captures >= 40% of the exact
     defect (kappa_a >= 0.4), Corollary B3's premise (the cheap part is only the next open walk) is wrong for the
     production chain, and A-late is back on the table.
5. **The closed-walk weight declines with depth and saturates.** At n = 16 it settles at w_3 of about 1.6-1.8 and w_4
   of about 1.7-2.2 by hidden layer 4. For the triangle at hidden layers 12-15 at n = 1024, predict w_3 in [1.2, 2.2].
   This is checkable offline from the MC caches' power sums plus a Stein estimate of K_3 restricted to the top
   collective modes, or directly in the 2049-run decomposition.

## 10. The minimal decisive experiment (E-N4; AWS, unbilled)

It decides V57_P5, and with it the frame's runtime content. It needs no 2049-run job.

**Inputs.**
- W_off{0..3}, truth_off{0..3}, chaindump_{0..3}.npz (out, mu, var, g4row per layer).
- One extra chain run per network with V30_DUMP_LEGS=1 and V29_DUMP_LAYERS=12,13,14,15. It dumps Y, the factor of
  M = C + eps mean(var), the young A_s, P_s, w2_s, t_s and the old-tier factors (Qc, FA, FP; QU, FA2, FP2).

**Steps.**
1. **Monte Carlo** (extend `mccache.py`): per-neuron power sums s1..s6 of the pre-activations at layers 11-15, in two
   independent halves of 2^27 samples each.
   - Per-neuron noise of kappa_5 is about 2e-3 sigma^5, giving T5 noise of about 1.6e-5 per half. That is comparable
     to the extrapolated signal (1-2.3e-5), so every aggregate must use split-half cross products. For example, the
     oracle MSE is <err - T5^(1), err - T5^(2)>, which is noise-free in expectation. The per-neuron fidelity in step 3
     is corrected by the split-half reliability.
   - Cost: 2^28 x 16 x 2 n^2 = 9e15 FLOPs per network in float32 GEMMs, about 40 minutes on one c7i.24xlarge. Four
     networks: about 3 instance-hours.
2. **Oracle value.** T5_l = g_5(mu_l, var_l) kappa5_l^MC / 120 and T6_l = g_6 kappa6_l^MC / 720, with (mu, var) from
   the chain dump.
   - At layer 15: report the OLS slope and explained share of err_15 = truth - out on T5, and the MSE ratio of
     err_15 - T5 and of err_15 - T5 - T6.
   - At layers 12-14: regress the injected error inj_l (DATA_FINDINGS decomposition) on T5_l. Push the explained part
     to the output with the mean gates, as in `ana_loc2.py`.
3. **Fidelity of the dressed 5-path.** Compute Pi5 from the dumped legs (Proposition 5.1, dressed). Report
   corr(60 Pi5, kappa5^MC) and the regression slope per layer, as note XLI's D2 did for the path class.

**Outputs.**
- the oracle MSE reduction from kappa_5 (and kappa_6) at the final readout and at layers 12-15, transported;
- the fidelity (corr, slope) of the dressed Pi5;
- the measured T5/T4 at n = 1024, against the 0.04-0.10 extrapolated here.

**Decision rule.**
- **Build V57_P5** (layers 12-15, derived weight 1, cold screen on 16 networks, then scored) iff both hold:
  - the transported kappa_5 oracle removes >= 12% of the final MSE;
  - the dressed Pi5 has corr >= 0.8 and slope in [0.5, 2] against the open-equivalent kappa_5.
- **If the oracle removes < 8%:** close N4's runtime branch. By Corollary B3 this also closes the analytic A-late
  merge: its affordable part is this class.
- **Prediction:** the oracle removes 2-6% (centre 4%), T5/T4 at layers 12-15 is 0.04-0.10, and the open-equivalent
  60 Pi5 explains less than half of kappa_5 (width exponent -0.9 against -0.5 for the open path).

**An addition to the registered X* run (no new runs).** Compute the n^3 analytic late defect d_L (section 8.1) from
the same leg dumps. Report its share of the exact defect's layer-15 local part (prediction 9.4: <= 15%). This measures
kappa_a for the cheap form directly.

## 11. What noncommutativity buys here, stated plainly

1. **It is genuine, and it decides costs.**
   - The chain's joint cumulants are mixed moments of the noncommuting second-chaos tuple (H_a) under two states: the
     vector states omega_ab and the trace.
   - The localization's heat flow acts on them by a closing map (Lap omega_ab = 2 tau(H_a . H_b)).
   - The local defect is built from the operator-valued carre du champ, which concatenates open walks.
   - The n^3 / n^4 boundary of the whole programme (notes XL-XLI, F1 Proposition 9, F5 6.2) is exactly the boundary
     between functionals of the noncommutative (Cipriani-Sauvageot) carre du champ, including its Cartan marginal and
     finite-Schur-rank pairings, and the commutative Schur carre du champ at full rank.
   - That is the frame's main theorem, and it is not a relabelling: it gives costs.
2. **It decides visibility.**
   - Readout truncations are seen at weight N_S(k). Omitted tables are seen at weight tr R(eps)/eps; for the
     covariance that is the first-chaos share, which explains MF.
   - Closed omissions are seen at weight k - (k/2) X_k/tau_k, through the Laplacian only, so only by the exact
     (finite-difference) defect.
   - Fitted closures produce a small false defect (response, not value).
3. **It does not make the defect cheaper.**
   - The affordable analytic defect is the next-order open walk (Corollary B3).
   - The closed part, which Theorem C shows the exact defect does see, sits behind the commutative wall.
   - Per-neuron sign rules exist only for even walks (positivity of the states).
4. **What to do.**
   - Run E-N4 (about 3 instance-hours) before any A-late build.
   - Add the kappa_a measurement of the cheap form to the registered X* run.
   - Expect both to come out negative for a runtime component. The theory stands either way.

## Appendix A. Numerical outputs

### A.1 Theorem B (carried-pair formula), `chk1_out2.txt` (n = 6, second-chaos scale 0.04)

```
seed 2: S=[1,2]     trD(fd) +1.349868e-03  formula +1.349855e-03  rel 9.9e-06 | -sum N_S err_open +1.349855e-03
        S=[1,2,3]   trD(fd) +1.812152e-03  formula +1.812140e-03  rel 6.4e-06 | -sum N_S err_open +1.812140e-03
        S=[1,2,3,4] trD(fd) -1.646574e-05  formula -1.648198e-05  rel 9.9e-04 | -sum N_S err_open -1.648198e-05
        S=[1,2,4]   trD(fd) -4.760948e-04  formula -4.761074e-04  rel 2.6e-05
seed 3: S=[1,2]     trD(fd) -2.263335e-05  formula -2.264894e-05  rel 6.9e-04   (classes k=3 +3.77e-04 w=2, k=4 -7.32e-04 w=1: cancel)
        S=[1,2,3]   trD(fd) +2.178178e-03  formula +2.178174e-03  rel 1.9e-06
        S=[1,2,3,4] trD(fd) +5.497447e-06  formula +5.472128e-06  rel 4.6e-03
        S=[1,2,4]   trD(fd) -2.258420e-03  formula -2.258436e-03  rel 7.1e-06
```

The seed-3 GC line shows Corollary 7's cancellation in miniature. The two classes' errors are of opposite sign, and
their weights (2 and 1) differ, so the defect almost vanishes while the error does not.

### A.2 Theorem C (closed walks), `chk2_out_n16.txt`

```
(a) n=16 L=4 samples=400000
  layer 1: rms tau3 2.392e-01 rms X3 2.058e-04 aggregate w3 +3.00 (corr +1.00) | rms tau4 1.878e-01 rms X4 1.932e-04 aggregate w4 +4.00 (corr +1.00)
  layer 2: rms tau3 2.293e-01 rms X3 8.241e-02 aggregate w3 +2.49 (corr +1.00) | rms tau4 1.524e-01 rms X4 6.467e-02 aggregate w4 +3.18 (corr +1.00)
  layer 3: rms tau3 2.377e-01 rms X3 1.511e-01 aggregate w3 +2.18 (corr +0.97) | rms tau4 1.539e-01 rms X4 1.163e-01 aggregate w4 +2.69 (corr +0.95)
(b) 2-layer, k=3: tr d_Sigma tau/tau = 1.500000, Lap tau/tau = -2.999999, tr H tau/(k tau) = 1.000000
    2-layer, k=4: tr d_Sigma tau/tau = 2.000000, Lap tau/tau = -3.999999, tr H tau/(k tau) = 1.000000
```

L = 8 (`chk2_out_n16L8.txt`):

```
(a) n=16 L=8 samples=200000
  layer 1: rms tau3 1.338e-01  rms X3 1.892e-04  aggregate w3 +3.00 (corr +1.00) | rms tau4 1.093e-01 rms X4 1.719e-04 aggregate w4 +4.00 (corr +1.00)
  layer 2: rms tau3 4.911e-01  rms X3 1.025e-01  aggregate w3 +2.70 (corr +1.00) | rms tau4 4.920e-01 rms X4 1.238e-01 aggregate w4 +3.50 (corr +1.00)
  layer 3: rms tau3 1.785e-01  rms X3 9.468e-02  aggregate w3 +2.21 (corr +1.00) | rms tau4 9.771e-02 rms X4 6.349e-02 aggregate w4 +2.71 (corr +1.00)
  layer 4: rms tau3 7.317e-02  rms X3 6.277e-02  aggregate w3 +1.76 (corr +0.98) | rms tau4 3.244e-02 rms X4 3.370e-02 aggregate w4 +1.99 (corr +0.94)
  layer 5: rms tau3 1.634e-01  rms X3 1.533e-01  aggregate w3 +1.61 (corr +0.99) | rms tau4 8.691e-02 rms X4 9.942e-02 aggregate w4 +1.72 (corr +0.99)
  layer 6: rms tau3 8.748e-02  rms X3 6.976e-02  aggregate w3 +1.81 (corr +1.00) | rms tau4 3.698e-02 rms X4 3.681e-02 aggregate w4 +2.03 (corr +0.99)
  layer 7: rms tau3 9.875e-02  rms X3 8.019e-02  aggregate w3 +1.81 (corr +0.99) | rms tau4 4.588e-02 rms X4 4.242e-02 aggregate w4 +2.20 (corr +0.99)
(b) 2-layer, neuron 0, k=3: tr d_Sigma tau / tau = 1.500000 (expect 1.5),  Lap tau / tau = -2.999999 (expect -3 + X-term),  tr H tau / (k tau) = 1.000000
     X_k at m=0: 0.000e+00 (zero: layer-1 means vanish at the base law)
(b) 2-layer, neuron 0, k=4: tr d_Sigma tau / tau = 2.000000 (expect 2.0),  Lap tau / tau = -3.999999 (expect -4 + X-term),  tr H tau / (k tau) = 1.000000
     X_k at m=0: 0.000e+00 (zero: layer-1 means vanish at the base law)
```

### A.3 Edgeworth order 5 against order 4, `chk4_out_n128.txt` and `chk4_out_n256.txt`

```
n=128 L=12 samples/half=2^18
 layer |  rms T3      rms T4      rms T5      rms T6  | T5/T4  T6/T4 | r = k4/(2 var^2) median
    0  | 4.619e-07  0.000e+00  1.621e-07  0.000e+00 | 162130607411951575236608.000  0.000 | -0.0002
    1  | 3.187e-03  1.162e-03  3.183e-05  4.720e-05 | 0.027  0.041 | 0.0413
    2  | 4.049e-03  1.567e-03  1.835e-04  3.489e-05 | 0.117  0.022 | 0.0714
    3  | 4.184e-03  1.547e-03  2.164e-04  8.920e-05 | 0.140  0.058 | 0.0912
    4  | 4.731e-03  1.874e-03  3.274e-04  1.419e-04 | 0.175  0.076 | 0.1206
    5  | 4.806e-03  2.050e-03  4.958e-04  3.050e-04 | 0.242  0.149 | 0.1408
    6  | 4.436e-03  1.852e-03  4.524e-04  2.477e-04 | 0.244  0.134 | 0.1649
    7  | 4.357e-03  1.696e-03  4.031e-04  2.417e-04 | 0.238  0.142 | 0.1648
    8  | 4.462e-03  1.760e-03  4.711e-04  2.305e-04 | 0.268  0.131 | 0.1866
    9  | 4.762e-03  1.580e-03  5.690e-04  2.169e-04 | 0.360  0.137 | 0.1907
   10  | 4.894e-03  1.756e-03  5.197e-04  3.161e-04 | 0.296  0.180 | 0.1879
   11  | 4.187e-03  1.807e-03  5.950e-04  5.413e-04 | 0.329  0.300 | 0.2212

n=256 L=12 samples/half=2^18
 layer |  rms T3      rms T4      rms T5      rms T6  | T5/T4  T6/T4 | r = k4/(2 var^2) median
    0  | 0.000e+00  4.802e-05  8.037e-08  0.000e+00 | 0.002  0.000 | -0.0002
    1  | 1.513e-03  5.555e-04  0.000e+00  7.862e-06 | 0.000  0.014 | 0.0203
    2  | 1.886e-03  6.838e-04  4.443e-05  2.374e-05 | 0.065  0.035 | 0.0332
    3  | 2.304e-03  7.266e-04  6.075e-05  2.141e-05 | 0.084  0.029 | 0.0451
    4  | 2.292e-03  8.747e-04  7.229e-05  2.271e-05 | 0.083  0.026 | 0.0591
    5  | 2.210e-03  8.207e-04  9.116e-05  4.526e-05 | 0.111  0.055 | 0.0667
    6  | 2.301e-03  8.879e-04  1.137e-04  5.951e-05 | 0.128  0.067 | 0.0780
    7  | 2.387e-03  7.869e-04  1.104e-04  5.010e-05 | 0.140  0.064 | 0.0850
    8  | 2.383e-03  8.129e-04  1.279e-04  5.173e-05 | 0.157  0.064 | 0.0883
    9  | 2.184e-03  8.838e-04  1.334e-04  8.442e-05 | 0.151  0.096 | 0.0964
   10  | 2.331e-03  8.767e-04  1.449e-04  7.660e-05 | 0.165  0.087 | 0.1052
   11  | 2.129e-03  9.225e-04  1.554e-04  9.197e-05 | 0.168  0.100 | 0.1143
```

### A.4 Proposition 1.1 and Proposition 5.1

```
(1) Gamma_op doubled-edge vs explicit: rel err 3.4e-16
(3) diag(Gamma_op) vs rowsum(Gamma_S): rel err 1.8e-16
(4) Schur-rank-1 pairing via operator products: 2.435967e+03 vs explicit 2.435967e+03 (rel 0.0e+00)
cost at n=nb=nin=1024: tensor route 1024 units; Gamma_op (doubled edge) 8 units; Schur-rank-1 pairing 8 units per rank
(a) Y_ij = L_j.v_i: True   (b) Pi5 hub form: True   (c) Pi6 hub form: True
(d) cross = Y + Q: True ; D21(open) = 2Y + Q: True ; cross = D21 - Y: True   (e) mixed walks: True True
```

## Appendix B. Relation to the literature

- **Noncommutative Dirichlet forms.** The derivation and bimodule language is Cipriani-Sauvageot's (J. Funct. Anal.
  2003), with the GNS and KMS extensions of Wirth (arXiv:2207.09247) and Vernooij-Wirth (arXiv:2303.15949).
- **Curvature.** Noncommutative curvature-dimension conditions are in Wirth-Zhang (arXiv:2105.08303) and in
  Carlen-Maas and Junge-Mei-Parcet. Schur multipliers as bimodule maps on matrix algebras are surveyed in
  arXiv:2510.17732.
- **Classical ingredients.** Bakry-Emery Gamma-calculus, Price's theorem, Stein's lemma and the Edgeworth expansion are
  classical. So are the chi-square cumulant formula (1.2) and the variable-elimination bound of Theorem D.
- **New here, to our knowledge:**
  - the carried-pair formula (3.1) and its weights N_S(k);
  - the table weight tr R(eps)/eps and its explanation of the measured MF, GC and K3 slopes;
  - Theorem C (closed-walk visibility through the rung-2 ladder);
  - the identification of the chain's n^3 / n^4 boundary with the noncommutative versus commutative carre du champ,
    via E_D[Gamma^op] = rowsum Gamma^S and finite Schur rank;
  - the dressed-hub forms of Pi5, Pi6 and the free cross term;
  - Proposition E's response principle for fitted closures.
- **Search.** An alphaXiv search on Dirichlet forms, Schur multipliers and carre du champ found no treatment of moment
  closures. Priority is not claimed.

## Referee report

Adversarial referee for frame N4. My checks are independent of `n4/` and live in `loc/ref_n4/`:
- `r1_thmB_exp.py` tests Theorem B on a readout whose exponential series converges exactly, f = exp(lam z).
- `r2_mf_share.py` measures the MF table weight, which the frame predicted but never computed.

I also re-derived every proof in sections 1-7 and read the V56 hub code (`estimator_final_v56.py` l. 1816-1850,
2725-2760) for the cost claims.

### R1. What is correct (verified)

1. **Theorem B (carried-pair formula) is an exact identity.** The author's check used relu with a truncated
   Edgeworth series, which diverges, so the agreement degraded to 3e-1. Mine uses f = exp(lam z), so
   Phi_S = exp(sum_{k in S} kappa_k lam^k / k!) is summable and d^A Phi = lam^|A| Phi. The model is the exact
   second-chaos one at general (m, Sigma). Formula (3.1) matches finite differences to FD precision:

   | seed | K | finite differences | (3.1) | relative difference |
   |---|---|---|---|---|
   | 0 | 2 | -5.20450e-4 | -5.20444e-4 | 1.2e-5 |
   | 0 | 3 | | | 9.6e-6 |
   | 0 | 5 | | | 3.4e-4 |
   | 1 | 2-5 | | | 1.5e-6 to 4.7e-5 |

   The truth E e^{lam z} has trH = -7e-9 (the control). The proof is the chain rule (1.1) plus the multi-index
   Hopf-Cole hierarchy, and the bookkeeping binom(C, A)/C! = 1/(A! B!) is right. **Corollary B1's** N_S(k) and the
   explicit last-layer form d_L = -2 g_5 Pi5 - (3/2) g_6 Pi6 also re-derive: the (1,4)+(2,3) pairs give 2 g_5 Pi5,
   and the (2,4)+(3,3) pairs give (3/2) g_6 Pi6.

2. **The table weight tr R(eps)/eps predicts the MF slope quantitatively. This is the frame's best result.**
   - **Setup.** The omitted object is o_a = sum_{b != c} W_ba W_ca C^y_bc. Its first-chaos part is f_a, with
     L^y = E[y x^T] estimated in split halves. I measured w = <f, o>/<o, o> on two fresh n = 64 nets, 2e6 samples.
   - **L = 8.** w = 0.99, 0.80, 0.69, 0.57, 0.51, 0.48, 0.45 (post-layers 1-7). The late layers dominate the
     transported error, which gives about 0.45-0.57 against the measured MF p of 0.46-0.61.
   - **L = 16.** w falls to 0.51 to 0.39 over layers 9-15, against the measured 0.32-0.47.
   - **Per neuron.** corr(f, o) is 0.97-0.99.
   - So the MF row of the table is now checked, not just asserted. The weights also order GC (1 and 2) and K3 (at
     most 3; measured 2.3-2.8) correctly. This gives the first derived, truth-free theory of the detector's slope.
3. **Theorem C's proof is correct**: Euler (degree k) plus rung 2 (Lap K_2 = -K_2) plus Leibniz.
4. **These also hold:**
   - Proposition 1.1, Proposition 1.2 and the doubled-edge formula;
   - Proposition 5.1's hub identities. Z_im = v_i . l_m needs L invertible, which the author states.
   - The curvature constants: Gamma_2 = (1/4)|Hess|^2, so CD(0, n); the ray flow is CD(-1/2, inf).
   - The mean-gate contraction 2 E Phi(alpha)^2 = 0.83, so E Phi^2 = 0.42.
   - The OLS relation |delta|/|err| = 2 p sqrt(X*).
   - Proposition E(2), which is F1's Theorem 1 applied to eps; it is valid because eps vanishes at point masses.

### R2. Mathematical errors and overclaims

1. **Theorem C is "checked" in name only for its new part.**
   - Part (b) of `chk2` checks the identity only at the first hidden layer. There K_3 = 0, and the code prints
     X_k = 0.
   - Part (a) does not test the identity at all. It *evaluates* w_k = k - (k/2) X_k/tau_k from Stein estimates of
     K_2 and K_3. The w_3 of about 1.6-1.8 at depth are therefore values of the formula, not a verification of it.
   - The proof is sound, so this is a labelling error. But the third-chaos term X_k, which carries the depth trend,
     has never been compared against a finite difference.
   - All numbers are at n = 16 and L <= 8, and the width dependence is unmeasured. The prediction "w_3 in [1.2, 2.2]
     at n = 1024" has no basis beyond n = 16.
2. **The X* prediction contradicts the frame's own Theorem C.**
   - SYNTHESIS Corollary 7 put closed walks among the invisible classes and got a visible share of 0.35-0.7.
     Theorem C moves them to *visible* at weight about 2, which pushes the visibility bound up.
   - Section 9.2 instead goes from 0.83 (the weight spread) times 0.9-0.97 (response mismatch) to X* = 0.35, by
     charging "the rest" to "third-chaos-cancelled closed walks, calibration".
   - Nothing in the frame cancels closed walks: X_k lowers w_3 from 3 to 1.7, not to 0. To get 0.35 from about 0.78,
     55% of the energy would have to be invisible, and no class is identified for it.
   - **Derived from N4's own results, X* is about 0.5-0.75.** Predictions 9.2 and 9.4 ("<= 15%", "20-50%") are
     guesses, not consequences of the theorems.
3. **Corollary B3 is labelled "proved" but holds only in the exact-table class model, and the decision rests on
   it.**
   - The production chain's late local defects also contain closure and compression maps: the lambda-core g4
     closure, tier-2 compression, the counterterms (1 + d)X, and the de-meaned, amplitude-fitted path term (V56
     subtracts its active mean and scales by V56_A).
   - Those local defects are per-neuron n^3 objects, and they are not next-order open walks. They are response
     mismatches (Proposition E), part real and part false defect.
   - So "the n^3 analytic defect *is exactly* the next open walk" is false for the production chain.
   - Hence "kappa_5 oracle < 8% closes A-late" does not follow. It closes the open-walk part of A-late, not the
     closure-defect part.
   - Separately, "V57_P5 dominates A-late" is an overclaim: V57_P5 drops the kappa_6 and pair-tree terms that A-late
     contains.
4. **"Gamma never produces a trace" is a second-chaos artifact.**
   - In the homogeneous network, d_j kappa_2 picks up tr(K_{3,j} K_2), so Gamma(kappa_2, kappa_2) contains
     sum_j tr(K_{3,j} K_2)^2, a closed contraction.
   - Theorem C's own X_k shows that third chaos mixes the open and closed classes.
   - At n = 1024 deep pre-activations are 78% higher chaos (DATA_FINDINGS), so the clean open/closed dichotomy, and
     with it B3's split, is a model statement.
5. **Smaller points.**
   - **Theorem A(2).** Its "only if" is too strong. D^2 Phi[Gamma] = 0 for every attainable Gram does not force
     D^2 Phi to vanish on the span when D^2 Phi is indefinite and the attainable Grams are not the full PSD cone on
     it. Only "if" is proved.
   - **The sign convention.** g_4 is proportional to (alpha^2 - 1) phi, so the sign stated as sign(1 - alpha^2) needs
     its error-sign convention spelled out.
   - **Proposition E's "about 10%"** (the 2.75 against 3) leaves out the C Lap lambda and 2 grad lambda . grad C
     terms of Lap(lambda C). It is unsupported.
   - **The T5/T4 width exponent of -0.9** rests on one doubling (128 -> 256) with a clipped-to-zero shallow layer.
     Extrapolating it two more doublings carries at least a 2x uncertainty, so the 2-6% band should read about
     1-12%.

### R3. Is the noncommutative content genuine?

Mostly not. It changes no computation.

1. **Section 1.2 is an amplification, not a noncommutative Dirichlet form.**
   - The "Cipriani-Sauvageot" structure is T_t = P_t (x) id_{M_N}: a commutative heat semigroup tensored with
     matrices.
   - Its A-valued carre du champ is the matrix Gram sum_j (d_j F)^T (d_j F).
   - No theorem of noncommutative Dirichlet-form theory is used: no complete Markovianity, closability, intrinsic
     metric, or NC curvature-dimension condition.
   - Proposition 1.1 is two index identities: (sum_j X_j X_j)_aa = sum_{j,c} X_ac^2, and
     tr(D_u X D_v X) = sum u_a v_b X_ab^2.
2. **The slogan "noncommutativity makes the cheap part cheap; the Schur product is the wall" is contradicted by the
   frame's own formula.**
   - The n^3 doubled-edge formula is Gamma^op = A (Wg o (P^T A)) P^T, which is a **Schur product** on birth pairs.
     This is visible in `chk3`.
   - The real criterion is Theorem D's treewidth of the contraction graph: commutative tensor-network elimination
     (variable elimination; Markov-Shi).
   - So the n^3/n^4 boundary is combinatorial, and the "op versus Schur" names attach to it after the fact.
3. **Gamma^op is not what any defect consumes.**
   - The pair maps are Schur functions of the tables, so their Hessians are Schur multipliers h_bc. Through the
     Mehler powers rho^k these have no finite Schur rank.
   - The defect therefore needs h o Gamma^S per pair, which is n^4, as section 5.2 concedes.
   - The only n^3 noncommutative objects are aggregates (Cartan marginals), and the author predicts no calibration
     gain from them. That is correct: the counterterms are truth-fitted on the same 100 scored networks, so
     truth-free aggregate calibration cannot beat them.
4. **Section 1.3 is a vocabulary change.**
   - (M_n, tr, omega_ab) as a "noncommutative probability space" restates the classical Gaussian quadratic-form
     cumulant formula (traces of products; Magnus, Mathai-Provost).
   - No noncommutative-probability machinery is applied: no freeness, no operator-valued subordination, no
     R-transform.
   - Lap omega_ab = 2 tau(H_a . H_b) is two derivatives of a quadratic. It is also already in F1 Proposition 7,
     which has Lap(L^T H^{k-2} L) = 2 tr H^k.
   - A genuine use would be asymptotic freeness of the (H_a), to get annealed closed walks from free cumulants. It
     is absent, and the quenched per-neuron residual (DATA_FINDINGS) suggests it would not pay.
5. **What is real.** Positivity of the two functionals on squares gives the even-walk sign rule. That is
   positivity of a Gram, not noncommutative geometry, and its main instance (the gain rule) is already exhausted.

### R4. Cost accounting

1. **B = Y M^-1 is free.**
   - Production runs V56_GAL = 0, the full solve `_Z56 = solve(M56, Y^T)`, so B = Z56^T is already formed.
   - If V57_CHOL (SYNTHESIS B) is adopted, the chain keeps only column norms of L^-1 Y^T. B then costs one more
     triangular pair, about 0.7-1.3 units.
   - V56's path term is de-meaned and amplitude-fitted. A Pi5 "at derived weight 1, no slope" is inconsistent with
     how the chain uses the sibling Pi4 term, and would likely need the same treatment.
2. **The old-tier cost is undercounted.**
   - The Hadamard P_s o w2_s o Z_s^2 needs P_s in neuron coordinates, and Hadamard products do not commute with the
     basis change.
   - For tier-factored sources this costs a second product (Qc FP_s) of r/n units: 0.31 for tier 1, 0.19 for
     tier 2.
   - That makes Pi5 about 6-11 units, not 4-8, and V57_P5 at layers 12-15 about +22-40 units (+11-19%), with break-
     even at -10% to -16% raw.
   - This strengthens the author's negative verdict.
3. **E-N4's compute is correct.** 2^28 x 16 x 2 n^2 = 9.0e15 FLOPs, about 40 minutes at about 4 TFLOP/s. The noise
   estimate is also right: Var k_5 is about 120 sigma^10 / N, so T5 noise is about 1e-5, comparable to the signal.
   Split-half cross products are therefore mandatory, as specified.

### R5. Relevance at the 1e-8 level

- **Runtime.** The frame proposes a component, V57_P5, that its own (corrected) numbers price as adjusted-negative.
  It argues that the affordable analytic defect cannot exceed it.
- **Limits of that argument.** By R2.3 the argument covers the open-walk content only. The closure and compression
  local defects of layers 12-15 are n^3 and uncovered.
- **Interpretation of the X* run.** Theorem C and the MF check matter here. They make a derived per-class slope
  model possible, and the slope is the quantity the registered rule tests for stability.
- **The value of E-N4.** It is real, cheap (about 3 instance-hours) and correctly noise-controlled. But it decides
  "carry kappa_5", which needs no Dirichlet-form theory, and it does not by itself close A-late (R2.3).

### R6. Novelty

- **Search.** An alphaXiv search on carre du champ / moment closure / Edgeworth / heat consistency found no prior
  statement of (3.1) or of the weights N_S(k).
- **Not new.**
  - Theorem D is standard treewidth contraction.
  - Proposition 1.2 is in F1 Proposition 7.
  - Section 1.3 is the classical quadratic-form cumulant formula.
  - Section 7 is textbook Bakry-Emery.
- **Probably new, though modest:**
  - (3.1) itself;
  - the table-weight rule, now checked on MF;
  - Theorem C's third-chaos correction (unverified numerically);
  - Proposition E's "a fitted closure has defect without error" (a direct corollary of F1 Theorem 1).

### R7. Verdict

- **What survives** is commutative:
  - an exact carried-pair identity;
  - the derived detector weights N_S(k) and tr R(eps)/eps, which quantitatively explain MF, GC and K3 (MF
    confirmed here);
  - the closed-walk visibility theorem.
- **What does not.** The noncommutative framing (Cipriani-Sauvageot, "vector state versus trace", "op versus Schur
  is the n^3/n^4 wall") is relabelling. Its one computational nugget, the doubled-edge formula, is a Schur
  identity, and the frame itself admits that no defect consumes Gamma^op.
- **Errors to fix:**
  - Theorem C's "checked" status;
  - the X* prediction, which is inconsistent with Theorem C;
  - B3's "proved" scope;
  - the old-tier cost;
  - Proposition E's 10%.
- **Status.** Survives as a theory contribution to reading the X* run, not as a noncommutative-geometry result and
  not as a runtime component.

**Recommended additions to the X* run** (no new chains):
1. Fit the slope per class using the derived weights: about 0.5 for table/covariance-type errors, 2 for kappa_3
   readout classes and closed walks, and 4 and 3 for the readout kappa_5 and kappa_6.
2. Test whether a two- or three-weight model beats the single slope on held-out networks.
3. Re-register X* at 0.5-0.75 if closed walks are visible as Theorem C says, so that the run actually tests
   Theorem C.
4. Validate X_k by a finite-difference check at a deep layer, using importance-weighted Stein estimates of
   K_2(m) = E[w(m, xi) z(xi)((xi - m)(xi - m)^T - I)] with w = exp(m.xi - |m|^2/2). These are analytic in m at
   fixed samples, at n <= 8.
