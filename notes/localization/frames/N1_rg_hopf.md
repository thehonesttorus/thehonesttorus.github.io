# N1. Stochastic localization as the exact renormalization group: the chain as a truncated effective action, the noncommutative walk algebra of its defect, and Hopf-algebraic (per-network) renormalization of its counterterms

Frame N1 of the NCG round. Status labels: **proved** (complete proof here or in a cited frame, re-derived),
**sketch** (argument given, routine details omitted), **measured** (numerical, outputs quoted, scripts named),
**conjecture**.
- Code (`loc/n1code/`): `n1_ct.py` GC chain with per-layer counterterms, measuring e, R = de/db, delta = e - Lap e and
  d delta/db; `n1_k3lib.py`, `n1_k3ct.py` the same on F1's dense kappa_3 chain K3 (D3 counterterms); `n1_ana.py`,
  `n1_k3ana.py` the renormalization conditions; `n1_final.py` the production-like pipeline (common counterterms, then
  the projected defect); `n1_sub.py` cheap subspaces; `n1_gap.py` the refit-gap estimator check; `n1_nonlin.py`
  linear-response check in the real chain; `n1_walks.py` walk-algebra identities; `n1_gen.py` 12 extra n = 32
  networks with Monte Carlo truth (`loc/n1data/`).
- Outputs (`loc/n1out/`): `ana_all.txt` (27 GC networks; `ana_L8.txt` first pass), `k3ana_6nets.txt`,
  `final_m3_all.txt`, `sub_L8.txt`, `gap_L8.txt`, `gap_n32.txt`, `nonlin.txt`, `walks.txt`, and the raw
  `ct_*.npz`, `k3ct_*.npz`.
- Networks and truth: the F1/F3 tiny networks (`loc/data`, `loc/f1data`) plus `loc/n1data`. Production facts come
  from `DATA_FINDINGS.md`, `SYNTHESIS.md`, notes XXXIX/XLI and `notes/ray-compiler` section 8.

## 0. Summary and verdict

**Verdict: marginal as a runtime component (the bottleneck is the in-budget defect, which this frame prices but cannot
remove); positive as theory, with three operational results and one zero-compute experiment that decides the
per-network question.**

1. **Localization is the exact RG, and the chain is a truncated effective action (proved; mostly known).** The
   conditional cumulant generating function K(lambda; m, Sigma) of the output is the Wilsonian effective action with
   background field m and remaining UV covariance Sigma; it solves the Polchinski equation, Eldan's process is its
   particle (Foellmer) realization, and F1's cumulant heat hierarchy is its lambda-expansion. Homogeneity makes the
   output a scaling operator of RG eigenvalue 1/2 in scaling variables. The new statement is a **two-flow curvature
   theorem** (Thm 2.3): exact layer maps intertwine the RG flows at successive depths; a closure fails to, its local
   defect is that commutator, and the error is the integral of the commutator over the depth x scale rectangle.

2. **Where noncommutativity genuinely enters: the walk algebra (proved, checked).** In the second-chaos sector every
   joint cumulant table is a sum of *vector-state* moments omega_ij(W) = L_i^T W L_j (open walks) and *trace* moments
   tau(H_i1 ... H_ik) (closed walks) of the noncommuting family {H_i} in M_n. The Polchinski flow **splits exactly**
   (Thm 3.3): the vector-state part solves the quadratic (carre du champ) flow by itself, and the trace part is driven
   by the Laplacian of the vector part, Delta_m omega_ij(W) = 2 tau(H_i W H_j). Consequences:
   - the RG residual of *any* open-walk truncation (every chain in this repo) is exactly a trace-state functional, so
     the part of the heat defect that sees the omitted closed walks **is** the n^4 wall (a "no free defect" theorem);
   - **visibility weights**: a closed walk omitted from a *carried* kappa_k table enters the defect with scaling
     dimension p = k (delta = 2k eps; checked: the ratio tends to 3.0 for the triangle), more visible than the open
     classes (p = 2 for the kappa_3 path, p = 1 for the kappa_4 path); a closed walk omitted from a readout that never
     carried kappa_k is invisible (p = 0). This corrects the visibility table of the synthesis and predicts that the
     production chain's finite-difference defect is dominated by wall terms that no in-budget analytic defect contains:
     in-budget capture kappa_a^2 = 0.45-0.8 (central 0.6).

3. **Hopf algebra and counterterms (proved at first order; measured).** Counterterm amplitudes compose through depth
   as a non-abelian group (characters of the shuffle Hopf algebra of layer-ordered insertions; the linear refit is its
   abelianization). The localization acts on error classes through the **grading** by scaling dimension, and the defect
   is the grading derivation applied to the error, delta = 2 Y eps, the analogue of Connes-Kreimer's beta = Y Res.
   Two operational consequences, measured on 27 tiny networks with the Gaussian-closure chain (n = 32-64, L = 8, 16)
   and 6 with the dense kappa_3 chain K3:
   - **the naive renormalization condition "zero the projected defect" fails** (MSE ratio 0.5-1.8 for production-like
     families, 2-9 with all 35 directions, 2.7 on K3), because the counterterm shapes have different scaling
     dimensions from the error they absorb (measured p_j = 0.3-1.0 against p_err = 1.1 for GC; 0.6-1.0 against 2.5
     for K3's D3 amplitudes). The **grading-corrected condition** (RC-Y: fit the counterterms to -(1/2) Y^{-1} delta,
     one calibrated scalar) works: it comes within 0.04-0.15 of the per-network truth oracle and beats the truth-fitted
     common amplitudes by 1.5-7x. On K3, the most production-like case, the common D3 amplitudes gain nothing
     held-out (ratio 1.11) while the truth-free per-network ones reach 0.62 (oracle 0.55); at depth 16 (GC) the
     common amplitudes again gain nothing (1.11) and RC-Y reaches 0.16;
   - **local (per-layer) renormalization must be BPHZ-ordered**: sequential layer-by-layer fixing with preparation
     reproduces the joint fit (0.12-0.31 against 0.09-0.19), independent per-layer subtraction overcounts
     catastrophically (6-69x worse than doing nothing), the forest-formula phenomenon.

4. **Per-network adaptivity is worth far more than the fitted constants (estimated from the repo's refit data).**
   The 103 counterterms fitted as constants buy about 4-6% of raw. The train/held-out gap of that fit measures the
   between-network spread of the optimal amplitudes: tr(G Sigma_b) = (K/2) gap = 74% of the MSE (K = 50, gap 2.97
   points), so a per-network choice of the same 103 amplitudes would remove **about 45-79% of each network's error,
   central 60%** (upper value from the gap formula, lower after the 1.1-1.5x overestimate the formula shows at small
   n). The span of the counterterm responses covers most of the network's collective error subspace (effective
   dimension about 130-170, consistent with the measured participation ratio); the fitted constants average a
   quenched (network-specific) projection to almost nothing.

5. **The practical target.** A truth-free per-network condition captures pi X of the MSE (pi = share of the error in
   the counterterm span, X = defect visibility): verified at small n to +-0.05 (GC; conservative by up to 0.15 on K3).
   It must be realized **in output space** (e + a P_U delta_hat, a projected defect merge), not as in-chain
   amplitudes: per-network amplitudes leave the linear regime (measured: the chain diverges at the 35-direction
   optimum), and the response matrix R costs 103 JVPs. U can be the top modes of the chain's own output covariance
   (free; holds 52-59% of the residual error after common counterterms at k/n = 0.1 at depth 8 and 85% at depth 16).
   The projection is also a denoiser: it suppresses incoherent error in delta_hat by k/N, so a sketched (cheap, noisy)
   analytic defect becomes usable (measured, depth 16: at noise equal to the defect itself, full merge 1.13 against
   projected 0.37).

**Costed component and prediction (section 6).** "Renormalized projected late defect" (N1-RPD): the open-walk
(in-budget) local defects of layers 12-15, partly free from V56's hub products, the rest by one hub-type product per
late source-layer (dense, about +55 units) or with sketched tangents (about +22 to +32 units; the covariance-map
cross term does not sketch), transported by mean gates and projected on the top-128 modes of the chain's own
last-layer covariance. Predicted raw x (0.77-0.90); adjusted 2.7-3.3e-9 sketched (central 3.0e-9, -5%), 3.1-3.6e-9
dense, against 3.14e-9. It pays only if the in-budget defect keeps X_ib >= 0.3 in the collective subspace and the
sketch noise stays below half the projected defect; a variant without the cross term (+7 to +17 units) reaches
2.5-3.1e-9 if E2 shows the cross term is small for this chain.

**Decisive experiment (section 7).** E0 needs no chain run: the stored 103-direction responses, the V56 errors and
the chain dumps give, per network, the oracle span share pi_R, the free-subspace share pi_U, and the exact
between-network amplitude spread. Rule: pi_U(128) >= 0.4 and pi_R >= 0.4 keep (iii) alive; then E1 on the pending X*
run's delta_hat (networks 0-3) measures pi X directly, and E2 on network 0 splits delta into in-budget and wall parts.

## 1. Setting

As in F1 section 1: X ~ N(0, I_n); F(x) = x_L with x_{l+1} = relu(x_l W_l), positively 1-homogeneous; u(m, Sigma) = E
F(m + Sigma^{1/2} Z); the chain is an estimator E(m, Sigma) exact on point masses and homogeneous; e(y) = E(y, I).
Heat operator H = d_Sigma - (1/2) Hess_m; defect D = H E; reduced defect delta(y) = e - y.grad e - Lap e, delta(0) = 2
tr D(0, I). Eldan path: Sigma_t = I/(1+t), dm_t = Sigma_t dW_t; s = 1/(1+t) is the posterior variance. eps := e -
truth (the chain's excess), err := truth - e = -eps. Units: 1 unit = 2 n^3 FLOPs, B = 1024 units. Production: V56 +
counterterms, raw 1.55e-8, C/B 0.203, adjusted 3.14e-9.

## 2. Stochastic localization as the exact renormalization group

### 2.1 The Wilson-Polchinski form

**Theorem 2.1 (proved; known in substance).** Let W(lambda; m, Sigma) = E exp<lambda, F(m + Sigma^{1/2} Z)> and
K = log W (finite for all lambda since F is Lipschitz). Then
1. (Wilson) d_Sigma W = (1/2) Hess_m W: the Boltzmann weight of the source term solves the linear heat equation;
2. (Polchinski) d_Sigma K = (1/2)(Hess_m K + grad_m K grad_m K^T), i.e. V = -K solves
   d_Sigma V = (1/2)(Hess V - grad V grad V^T);
3. along Eldan's path, K(lambda; m_t, Sigma_t) is the conditional CGF of F(X) given F_t; at t = 0 it is the full
   effective action (the output's cumulant generating function), as t -> infinity it tends to the bare action
   <lambda, F(x)>;
4. expanding in lambda gives F1's cumulant heat hierarchy H kappa_k = (1/2) sum_{A u B = [k]} grad kappa_A (x)
   grad kappa_B.

*Proof.* (1) Price's theorem for each fixed lambda. (2) Hopf-Cole: W = e^K. (3) Gaussian posterior (F1 Thm 1).
(4) Taylor coefficients of K in lambda are the joint cumulants. QED.

- **Reading.** The localization is the Polchinski flow read from the IR end: m_t is the field already integrated
  (observed), Sigma_t the UV covariance not yet integrated. The process m_t is the Polchinski particle of Shi-Tian-Zhang
  (arXiv 2510.04460, their Thm 4) and the Foellmer drift of Bauerschmidt-Bodineau-Dagallier; the elicit survey
  (`elicit/rg_hopf`) has the conventions. Bruned-Laubie-Minguella (arXiv 2608.31049, Aug 2026) prove that a decorated
  Feynman-graph ansatz with BPHZ (Connes-Kreimer) subtractions solves the Polchinski equation: the two languages used in
  this frame are the same renormalization.
- **The source is the network.** Unlike QFT, the "interaction" is the source term lambda.F(x) itself; the Gaussian
  measure is free. So the "vertices" of the effective action are the cumulants of F under the remaining noise, and a
  cumulant chain is a vertex-expansion truncation evaluated at one background point.

**Proposition 2.2 (scaling variables; proved, F1 Cor 2 and Thm 3 restated).** With homogeneity
K(lambda; c m, c^2 Sigma) = K(c lambda; m, Sigma). In scaling variables y = m / sqrt(s), T = -log s, the RG flow of
the mean is generated by A = (1/2)(Lap + y.grad) and the truth is an eigen-operator: A v = v/2 (dimension 1/2 in the
variance unit). The Euler-Stein ladder is the eigen-equation projected on the chaos (Taylor) components at y = 0, and
the reduced defect delta = (1 - 2A) e is the estimator's anomalous dimension operator applied to it. QED (F1).

### 2.2 The chain as a truncated effective action, and the two-flow curvature theorem

The chain carries tables T_l (mean, covariance, sources, closures) for l = 0..L, with T_0 the exact input law,
T_l = Phi_l(T_{l-1}) (closure readout composed with the linear transport), output e = pi(T_L). For any table-valued
function T(m, Sigma) let its **RG residual** be D(T) := H T - R(T), with R the right side of the hierarchy (built from
gradients of lower tables only). Exact cumulant tables have D = 0 at every depth (Thm 2.1 applied to F_l).

**Definition (curvature).** For a layer map Phi and an input table T, d(Phi; T) := D(Phi(T)) - DPhi . D(T)
= DPhi . R(T) - (1/2) D^2 Phi[grad T, grad T] - R(Phi(T)). It depends only on (T, grad T, Phi) (F1 Prop 8).

**Theorem 2.3 (two-flow curvature; proved from F1 Thm 1 and Prop 8).**
1. The exact layer map Phi_l^ex intertwines the RG flows at depths l-1 and l: d(Phi_l^ex; T) = 0 whenever T is an exact
   table. A closure Phi_l does not; d_l := d(Phi_l; T_{l-1}) is the commutator of the closure with the Polchinski
   flow, evaluated on the chain's own trajectory.
2. At every (m, Sigma): tr D(e) = sum_l S_l . tr d_l, S_l = Dpi DPhi_L ... DPhi_{l+1}.
3. truth - e(0, I) = - E int_0^infinity sum_{l=1}^{L} < S_l d_l (m_t, Sigma_t), Q_t > dt, Q_t = Sigma_t^2.

*Proof.* (1) The conditional law of x_l given F_t is the pushforward of that of x_{l-1}; conditional cumulants solve
the hierarchy (Thm 2.1(4) for F_l); hence Phi^ex maps D = 0 tables to D = 0 tables, and by the chain rule identity
D(Phi(T)) = DPhi D(T) + d(Phi;T) the local term vanishes. (2) Iterate the chain rule (F1 Prop 8) at each (m, Sigma).
(3) Insert (2) into F1 Thm 1. QED.

- **Reading.** Two flows act on the space of truncated effective actions: depth (Phi_l) and scale (the localization).
  Exactness is their commutation; the defect is the curvature; the error is the integrated curvature over the
  rectangle [0, L] x [0, infinity) (a discrete non-abelian Stokes formula: going "depth first, then localize" versus
  "localize first, then depth", where the latter is exact because the chain is exact on point masses).
- **Per-layer, truth-free attribution.** err_l := -E int <S_l d_l, Q> dt sums to the error and has a truth-free
  integrand. At the base point and with a per-layer profile exponent p_l, err_l ~ -S_l tr d_l(0, I)/p_l: the object
  A-late of the synthesis injects. DATA_FINDINGS' injection shares (65% from layers 12-15) are the truth-side image of
  this decomposition.
- **What is new.** Items 1-2 are F1 Prop 8; item 3 is F1 Thm 1. The content added here is the identification of d_l
  as the commutator of two flows, which is what licenses the renormalization conditions of section 5 to be imposed
  layer by layer (section 4.3).

### 2.3 The grading: scaling dimensions of error classes and of counterterms

For a function G(m, Sigma) define the **localization log-derivative**
L(G) := 2 d/ds|_{s=1} E_{m ~ N(0,(1-s) I)} G(m, s I) = 2 E tr D_G(0, I) (= delta_G(0) for homogeneous G).

**Theorem 2.4 (grading; proved).**
1. L is linear, L(u) = 0, L(e) = delta(0); hence delta(0) = L(eps) exactly.
2. If eps = sum_c eps_c with averaged-path profiles E_m eps_c(m, sI) = s^{p_c} eps_c(0, I), then
   delta(0) = 2 sum_c p_c eps_c(0, I) =: 2 Y eps, where Y is the grading by scaling dimension; on such a decomposition
   eps = (1/2) Y^{-1} delta.
3. If a counterterm b_j enters the estimator with response R_j(m, Sigma) = dE/db_j of profile exponent p_j, then
   d delta(0)/d b_j = 2 p_j R_j(0, I).

*Proof.* (1) Linearity of expectation and derivative; u is a martingale along the path, so its averaged profile is
constant; F1 Thm 1 deterministic form gives the identification with tr D. (2) Differentiate s^{p_c} at s = 1. (3) Apply
(1)-(2) to R_j = d eps/d b_j. QED.

- **Content.** The theorem is elementary; its content is the hypothesis that classes have *neuron-independent*
  exponents (F1's quenched-amplitude / annealed-exponent conjecture), which section 5 measures for counterterm shapes:
  corr(d delta/d b_j, R_j) = 0.86-0.98 for four of the five GC shapes (0.61-0.72 for the Mehler-1 covariance shape)
  and 0.71-0.91 for K3's D3 shapes, with class exponents that differ from the error's.
- **Connes-Kreimer.** In CK the renormalization-group generator is the grading derivation Y of the Hopf algebra, and
  the beta function is Y applied to the residue of the counterterm, beta = Y Res gamma_- (Connes-Kreimer II,
  hep-th/0003188); the residue is recovered by Y^{-1}. Item 2 is the same algebra: the localization acts on error
  classes by s^Y, the defect is 2Y eps (the "beta function" of the truncation), and the error is (1/2) Y^{-1} of it.
  The single-coefficient merge e - delta/(2p) is the approximation Y ~ p.

## 3. Where genuine noncommutativity enters: the walk algebra of the second chaos

### 3.1 The joint law of the second-chaos field is a noncommutative distribution

Truncated at the second Wiener chaos (note XLI), a layer's pre-activations are
z_i = mu_i + L_i.x + (1/2)(x^T H_i x - tr H_i), with first-chaos vectors L_i in R^n and symmetric second-chaos kernels
H_i in M_n(R) (H_i = sum_b sum_m P_im w2_m l_m l_m^T from the sources). The H_i do not commute. Let A be the algebra
they generate, tau = Tr, and for a word W in A let omega_ij(W) = L_i^T W L_j (vector functionals).

**Proposition 3.1 (joint cumulants; proved, checked).** For x ~ N(m, sI) put L_i' = L_i + H_i m. For k >= 2,

    kappa(z_i1, ..., z_ik) = (s^k / 2k) sum_{sigma in S_k} tau(H_isigma1 ... H_isigmak)
                           + (s^{k-1} / 2) sum_{sigma in S_k} L'_isigma1^T H_isigma2 ... H_isigma(k-1) L'_isigmak,

and kappa_1(z_i) = mu_i + L_i.m + (1/2) m^T H_i m + (s - 1) tr H_i / 2. So kappa_3(i,j,k) =
s^2 [L_i'^T H_j L_k' + L_j'^T H_k L_i' + L_k'^T H_i L_j'] + s^3 tau(H_i H_j H_k).

*Proof.* x = m + sqrt(s) g makes z_i a quadratic form in g with linear coefficient sqrt(s) L_i' and matrix s H_i;
its joint CGF is -(1/2) log det(I - s H_t) - (s/2) tr H_t + (s/2) L_t'^T (I - s H_t)^{-1} L_t' with H_t = sum t_i H_i,
L_t' = sum t_i L_i'. Extract the multilinear coefficient. QED.
*Check* (`n1_walks.py`, n = 6, three neurons, 1e7 samples): kappa_3 Monte Carlo -0.0625 +- 0.0023 against the formula
-0.0613 (open -0.0824, closed +0.0211); mean and covariance to 3 digits.

- **Reading.** The joint law of a layer field in the second-chaos sector is the *-distribution of the noncommuting
  family (H_i) under two kinds of functionals: the trace (closed walks, cycles) and the vector functionals
  (open walks, paths). Single-neuron (diagonal) tables see only one operator H_i (commutative: spectral moments
  tr H_i^k and the spectral measure of H_i in the vector L_i); the joint tables (D21, K22, K31, the gate covariance)
  see words in several H_i. This is the precise place the synthesis located the cost wall. What follows shows that it
  is also where the heat defect lives.

### 3.2 The localization is a derivation on the walk algebra

Localization moves the vector data and fixes the operators: L_i -> L_i + H_i m (Prop 3.1). Write delta_e for the
derivative in the direction e of the input mean.

**Lemma 3.2 (proved).** For words W, V and neurons i, j, k, l:
1. delta_e L_i' = H_i e, delta_e H = 0; hence Lap_m omega'_ij(W) = 2 tau(H_i W H_j)  (the Laplacian closes the open
   walk) and grad_m tau(.) = 0;
2. <grad omega'_ij(W), grad omega'_kl(V)> = omega'_jl(W^T H_i H_k V) + omega'_jk(W^T H_i H_l V^T)
   + omega'_il(W H_j H_k V) + omega'_ik(W H_j H_l V^T)  (the carre du champ concatenates two open walks through H H).

*Proof.* grad_m (L_i'^T W L_j') = H_i W L_j' + H_j W^T L_i'; differentiate again and sum over an orthonormal basis
(sum_e (H_i e)^T W (H_j e) = tau(H_i W H_j)); for 2 expand the inner product of the two gradients. QED.

So on the walk algebra the heat semigroup acts through one derivation (the endpoint move L -> H e); its Laplacian maps
vector functionals to the trace and its carre du champ Gamma is concatenation. This is a noncommutative Dirichlet-form
structure in the Cipriani-Sauvageot sense (a derivation into the bimodule spanned by the open walks), in finite
dimension and with an explicit symbol.

### 3.3 The flow splits; no free defect; visibility weights

Split every exact cumulant table into its open part T^o (degree 2 in L') and closed part T^c (degree 0 in L').

**Theorem 3.3 (proved; checked).** In the second-chaos sector:
1. (Exact split of the Polchinski hierarchy) d_s T^o = R(T^o) and d_s T^c = (1/2) Lap_m T^o. The vector-state
   hierarchy is closed under the RG by itself; the trace hierarchy is driven by the Laplacian of the vector one.
2. (Residual of the open-walk truncation) A chain that carries the open parts exactly and omits the closed parts has
   RG residual D(T^o) = -(1/2) Lap T^o = -d_s T^c: a pure trace-state functional.
3. (Visibility) For a readout h(kappa_1, ..., kappa_K) that reads a carried kappa_k table, the omitted closed walk
   C_k = s^k tau(...) enters the mean defect as h_k D_k = -h_k k s^{k-1} tau(...), against an error h_k C_k: scaling
   dimension **p = k** (delta = 2k eps) at leading order in the non-Gaussianity. For a readout that never carried
   kappa_k, the omitted closed walk does not enter the first-order defect (p = 0), while its omitted open walks enter
   through the gradient Grams with p = 2 (kappa_3 path) and p = 1 (kappa_4 path) (F5 Prop 5.4 in these units).
4. (No free defect) The trace part of the defect of any carried table of order k >= 3 (the D3 diagonal, the (2,1)
   and (3,1) slices, the gate covariances) is a trace moment tau(H_i1 ... H_ik). Open walks reach the current layer
   through its covariance (the Schur identity G_birth = At^T C^{-1} At of note XLI turns L^T H^q L into hub products,
   n^3 per source-layer); a trace does not: tr H_i^3 = tr (At d(c_i) At^T C^{-1})^3 needs an n x n product per neuron,
   n^4 per source-layer, and the joint traces are the K_4-contraction class of note XL. No in-budget form of the defect
   contains this part.

*Proof.* (1) tau-terms have zero m-gradient, so R(T) = R(T^o); Lap T^c = 0. The exact table satisfies d_s T = (1/2)
Lap T + R(T). Fix (m, s) and regard both sides as polynomials in L (equivalently in L' = L + H m, an affine
bijection): open parts and R(T^o) (a sum of products of two gradients, each of degree 1 in L') are homogeneous of
degree 2, closed parts and Lap T^o (Lemma 3.2: each derivative replaces one L' by H e) of degree 0. The identity holds
for all L, so its degree-2 and degree-0 components hold separately. (2) D(T^o) = d_s T^o - (1/2) Lap T^o - R(T^o) =
-(1/2) Lap T^o by (1). (3) Chain rule for H on h(kappa) (F1 Prop 7): the omitted C_k has zero gradient, so it enters
only through h_k D_k with D_k = -d_s C_k = -k C_k / s; the error is h_k C_k at first order. The correction terms come
from the dependence of h_k on (kappa_1, kappa_2), which are exact here; they are of relative order of the non-Gaussian
cumulants. For a readout without kappa_k, C_k appears nowhere in the first-order defect. (4) Lemma 3.2(1); the trace
identity by substituting the Schur form of the birth Gram; note XL Thm 1(iii) for the joint classes. QED.

*Checks* (`n1_walks.py`, `n1out/walks.txt`):
- Polchinski residual of the exact joint kappa_3 by finite differences: 6.2e-8 (against R_3 = -0.2355).
- Residual of the open-walk truncation: -0.09057, against -3 s^2 tau(H_0 H_1 H_2) = -0.09057.
- Visibility of the triangle omitted from a carried kappa_3 (Edgeworth readout of relu, single neuron): tr D / eps =
  4.33, 3.74, 3.39 at H scale 1, 0.5, 0.25, converging linearly to 3.0, i.e. delta/eps -> 6 (p = 3).

- **Operator-valued reading.** The open parts are matrix elements of the M_N-valued map W -> L'^T W L' (N neurons),
  an operator-valued vector state on the algebra of the H_i; the closed parts are its trace. Lemma 3.2 says the heat
  generator maps the first onto the second. Free probability (asymptotic freeness of the H_i, notes F4 and
  `notes/ncg-probability`) predicts the *annealed* trace moments of words in distinct H_i, and for the joint triangle
  that prediction is zero (W1: annealed mean zero). So the trace part of the defect is purely quenched: no counterterm
  constant can absorb it and no free-probability closure can supply it.

**Corollary 3.4 (visibility bound, corrected; proved given the class model).** With uncorrelated per-neuron class
components of energies E_c and dimensions p_c, the explained share of the full defect is
X*_full = (sum p_c E_c)^2 / (sum E_c sum p_c^2 E_c), and that of an in-budget defect that omits the trace terms is
X*_ib = (sum_{open} p_c E_c)^2 / (sum_all E_c . sum_{open} p_c^2 E_c).
- **Production.** The oracle attribution (note XXXI) and W1 give: kappa_3 readouts about 40% of the MSE, of which the
  triangle / joint-gate (closed) part is the W1 pair-gate class (at most a third of the MSE); the kappa_4 diagonal
  15-30% after V56, whose residual includes the 4-cycle (24-30% of the path term, note XLI); about 30% truncation and
  calibration (p about 0). Scanning closed shares 0.15-0.35 (p = 3-4), open shares 0.35-0.6 (p = 1-2), invisible
  0.2-0.3:
  - X*_full = 0.45-0.68;
  - X*_ib = 0.30-0.55;
  - the wall terms hold 35-70% of the *defect's* energy even when they are a minority of the error (they carry
    weights 3-4 against 1-2), so the analytic in-budget capture is kappa_a^2 = X*_ib / X*_full = 0.45-0.8.
- **Correction to the synthesis.** Its Corollary 7 put the readout's closed walks at weight 0. That holds for a
  readout that never carried the class; for the production chain, whose D3 and g4 tables are carried and omit the
  triangle and the 4-cycle, the omitted closed walks are visible at weights 3 and 4. The finite-difference X* run will
  therefore look better than any affordable analytic defect can be (prediction P3 in section 7).

## 4. Hopf-algebraic structure of the counterterms

### 4.1 The counterterm group is non-abelian, and the linear refit is its abelianization

Let the alphabet be the counterterm directions a = (l, X) (layer, statistic), and theta_l(b_l) the rescaling of the
statistics read at layer l. The renormalized chain is e(b) = pi o Phi_L o theta_L(b_L) o ... o Phi_1 o
theta_1(b_1)(T_0).

**Proposition 4.1 (sketch).** The Taylor expansion of e(b) is a sum over depth-ordered words in the alphabet:
e(b) - e(0) = sum_k sum_{a_1 <= ... <= a_k (in depth)} c(a_1 ... a_k) b_{a_1} ... b_{a_k}. Configurations compose by
the deconcatenation (depth-ordered) coproduct; the maps b -> e(b) for varying layer maps are characters of the shuffle
Hopf algebra on the alphabet (Chen series of a product integral). First-order coefficients c(a) = R_a are the
infinitesimal characters (the stored responses of `notes/ray-compiler` section 8); the second-order cross-layer
coefficients c(a a') with a before a' are the pre-Lie insertions "the counterterm at a changes the tables on which
the shape at a' is evaluated". The linear output-metric refit (`code/refit.py`) is the abelianization (Lie-algebra
level).
- *Repo evidence (conjecture for the mechanism).* The first-order refit has twice been pessimistic by 0.6-0.7 points
  (predicted -3.07% / -3.16%, measured -3.76% / -3.83%): a favourable second-order (pre-Lie) term of about a fifth
  of the first-order gain.

### 4.2 The defect is the grading derivation of the counterterm character

By Thm 2.4(3), at first order delta(b) - delta(0) = 2 sum_a p_a R_a b_a = 2 Y (R b): the defect response of a
counterterm configuration is the grading applied to its infinitesimal character. Renormalization conditions written
on the defect must therefore undo the grading before they are compared with the error (section 5.1).

### 4.3 Birkhoff / BPHZ: renormalization conditions are local in depth and must be prepared

**Theorem 4.2 (locality at first order; sketch).** d delta / d b_l = 2 S_l d tr d_l / d b_l + O(|b| |d|).
*Proof sketch.* delta = 2 sum_l' S_l' tr d_l' (Thm 2.3). A counterterm at layer l changes the tables at every l' >= l;
the induced change of d_l' (l' > l) and of S_l' multiplies a closure defect, which is small uniformly near the
trajectory (d vanishes for exact layer maps). Only the counterterm's own local defect contributes at first order. QED.

So each layer's counterterm can be fixed by a condition on that layer's **prepared** local defect (all earlier
counterterms in place), which is the BPHZ recursion phi_-(Gamma) = -R[phi(Gamma) + sum phi_-(gamma) phi(Gamma/gamma)]
with the depth order in place of the subgraph order. What it excludes: fixing each layer independently against the
unprepared defect. The response spans of different layers overlap (they all excite the collective output modes), so
independent subtraction counts the overlap once per layer, exactly the overcounting the forest formula removes.

*Measured* (`n1_ana.py`, GC, 27 networks, RC-Y of section 5 with all directions (35 at L = 8, 75 at L = 16), linear
model):

| group | joint fit | sequential (forward, prepared) | block-Jacobi (independent per layer) |
|---|---|---|---|
| n32 L8 (16 nets) | 0.113 | 0.156 | 13.1 |
| n48 L8 (4) | 0.194 | 0.314 | 6.37 |
| n64 L8 (4) | 0.090 | 0.131 | 18.3 |
| n64 L16 (3) | 0.107 | 0.121 | 68.6 |

### 4.4 What Birkhoff says about the fitted constants: universality fails for this problem

In QFT the counterterms are local (independent of the IR state) because divergences are UV objects. Here the analogous
property would be that the optimal amplitudes are annealed (network-independent). The repo's own refit data say
otherwise (section 5.3): the between-network spread of the optimal amplitudes is an order of magnitude larger than
their mean. The counterterm span is a good *basis* for each network's error, and the constants are its quenched
average. Renormalization here must be per network, i.e. by conditions, not by universal constants.

## 5. Truth-free per-network renormalization conditions

### 5.1 Three candidate conditions

With R = de/db (per network) and Dd = d delta/db:
- **RC-G** (naive: zero the projected defect): b = argmin |delta + Dd b|^2.
- **RC-Y** (grading-corrected): b = argmin |a delta - R b|^2 with one class scalar a ~ -1/(2 p_err), i.e. fit the
  counterterms to -(1/2) Y^{-1} delta in the single-degree approximation. Equivalently, RC-G with each defect
  response rescaled by p_err / p_j (Thm 2.4(3)).
- **RC-L** (local, prepared): RC-Y solved layer by layer in depth order (section 4.3).

By Thm 2.4, RC-G returns b_j = (p_err / p_j) b_j^opt on the span (in the diagonal case), so it is unbiased only if
every counterterm shape has the error's scaling dimension.

**Proposition 5.1 (value of a projected condition; proved under the stated orthogonality).** Let a.delta = X eps + w
with w orthogonal to eps, |w|^2 = X(1-X)|eps|^2 (the OLS form), P the projector on the counterterm span (dimension k),
pi_e and pi_w the shares of eps and w in its range, and nu an incoherent error of the defect estimate with energy N
sigma^2. Then the MSE removed by RC-Y is (2X - X^2) pi_e - X(1-X) pi_w - a^2 k sigma^2 / |eps|^2, which is pi X when
pi_w = pi_e = pi, against X - a^2 N sigma^2 / |eps|^2 for the full merge. The projection wins as soon as the defect
estimate's incoherent noise energy exceeds (1 - pi) X |eps|^2 / (1 - k/N). *Proof.* Expand |eps - P(a delta + a
nu)|^2, with nu independent of everything and isotropic, and the assumption that P w is uncorrelated with P eps (w is
orthogonal to eps by construction; that its projection stays orthogonal to the projected error is the assumption).
QED. The formula is checked against measurement in section 5.2.

### 5.2 Small-n measurement (GC chain with per-layer counterterms)

**Setup** (`n1_ct.py`). GC chain (dense covariance, Mehler order 14) with five per-layer counterterm shapes, the
analogues of the production amplitudes: var:l, cL:l (Mehler order 1 of the off-diagonal covariance), cN:l (orders >= 2),
m3:l (unit kappa_3 readout shape A3/6, the D3 analogue), m4:l (unit kappa_4 shape A4/24, the g4 analogue). For each
network: e, MC truth, R (central differences in b), delta = e - Lap e (2n+1 runs, h = 0.05; h-check against 2h: 3e-4
to 1e-3 of rms delta), Dd (Laplacians at b = +-0.01). Linear-response evaluation; COMMON = truth-fitted amplitudes
pooled over the other networks of the same depth with inner-LOO ridge (the V47_CAL analogue); ORACLE = per-network
truth fit; the merge slope a is leave-one-out across the pool (a = -0.44 to -0.47). Families: m3 (7 directions,
k/n = 0.11-0.22, the production-like ratio k/n ~ 0.1), mean (m3 + m4, 14), all (35).

Networks: n32 L8 (4 from `loc/data` plus 12 new ones in `loc/n1data`, `n1_gen.py`, 2e6 MC samples), n48 L8 (4),
n64 L8 (4), n64 L16 (3). Geometric-mean MSE ratios to the uncorrected chain (`n1out/ana_all.txt`; `ana_L8.txt` is
the first pass on the 12 original L8 networks):

| group | X | merge | m3: oracle / common / RC-G / **RC-Y** | mean: oracle / common / RC-G / **RC-Y** |
|---|---|---|---|---|
| n32 L8 (16) | 0.84 | 0.113 | 0.21 / 0.77 / 0.86 / **0.29** | 0.05 / 0.58 / 0.78 / **0.16** |
| n48 L8 (4) | 0.81 | 0.189 | 0.37 / 0.74 / 1.07 / **0.48** | 0.19 / 0.57 / 1.36 / **0.34** |
| n64 L8 (4) | 0.92 | 0.077 | 0.31 / 0.63 / 0.78 / **0.35** | 0.17 / 0.53 / 0.50 / **0.23** |
| n64 L16 (3) | 0.83 | 0.107 | 0.08 / 1.11 / 1.48 / **0.16** | 0.02 / 0.91 / 1.84 / **0.12** |

(k/n of the m3 family: 0.22, 0.15, 0.11, 0.23; of the mean family twice that.)

- **RC-Y comes within 0.04-0.15 of the per-network truth oracle and beats the truth-fitted common amplitudes by
  1.5-7x.** It is truth-free except for one class scalar. At depth 16 the common amplitudes are useless held out
  (1.11, 0.91), as in production, while RC-Y removes 84-88%.
- **The value formula holds**: on the L8 groups 1 - pi X predicts 0.31, 0.49, 0.37 (m3) and 0.18, 0.34, 0.24 (mean)
  against measured 0.29, 0.48, 0.35 and 0.16, 0.34, 0.23.
- **RC-G fails** (0.5-1.8 for these families; 2-9 with all 35 directions). The anomalous dimensions explain it:

| group | p_err | var | cL | cN | m3 | m4 |
|---|---|---|---|---|---|---|
| n32 L8 | 1.12 | 0.38 (corr 0.94) | 0.42 (0.67) | 0.91 (0.86) | 0.55 (0.93) | 0.88 (0.95) |
| n48 L8 | 1.12 | 0.41 (0.96) | 0.40 (0.66) | 1.04 (0.86) | 0.54 (0.93) | 0.90 (0.95) |
| n64 L8 | 1.12 | 0.40 (0.98) | 0.43 (0.72) | 0.99 (0.89) | 0.59 (0.96) | 0.90 (0.97) |
| n64 L16 | 1.06 | 0.26 (0.92) | 0.29 (0.61) | 0.74 (0.88) | 0.40 (0.88) | 0.72 (0.91) |

  (p_j = <Dd_j, R_j>/(2|R_j|^2), median over directions; corr(Dd_j, R_j) in brackets; p_err = -1/(2a).) The defect
  response of each shape is its own response times 2 p_j to correlation 0.86-0.98 (cL 0.61-0.72), as Thm 2.4(3)
  requires, with p_j between 0.26 and 1.04 while the error has 1.06-1.12. For m3, p_err/p_m3 ~ 1.9 at depth 8 and 2.7
  at depth 16: RC-G overshoots the m3 amplitudes two- to threefold; 1 - pi + (1 - 1.9)^2 pi = 0.84 predicts its
  measured 0.78 (n64 L8).
- **After the common counterterms** (C+Y: COMMON first, then RC-Y on the corrected chain's own defect delta + Dd b_c
  with a refitted slope; pooled over all 27 GC networks), the defect stays visible (X = 0.72-0.83) and the truth-free
  per-network step removes a further 61-72%:

| family | common | common + RC-Y | common + full merge |
|---|---|---|---|
| m3 | 0.776 | 0.297 | 0.098 |
| mean | 0.599 | 0.180 | 0.108 |
| late (layers >= L/2) | 0.459 | 0.177 | 0.129 |
| all | 0.357 | 0.101 | 0.099 |

- **Noise** (synthetic incoherent noise of s times rms(delta) per neuron, 20 draws):

| s | n64 L8: full merge / RC-Y 35 / RC-Y mean (14) | n64 L16: full merge / RC-Y mean (30) |
|---|---|---|
| 0 | 0.08 / 0.10 / 0.26 | 0.16 / 0.18 |
| 0.5 | 0.38 / 0.26 / 0.32 | 0.37 / 0.27 |
| 1 | 1.26 / 0.76 / 0.51 | 1.01 / 0.57 |
| 2 | 4.82 / 2.66 / 1.31 | 3.73 / 1.82 |

  The crossover sits where Prop 5.1 puts it: the projection pays once the defect estimate's noise is about half the
  defect.
- **Locality (section 4.3)** at depth 16: joint 0.107, sequential 0.121, block-Jacobi 68.6.
- **Nonlinearity** (`n1_nonlin.py`, real chain at the fitted amplitudes, n48/n64 L8): m3 amplitudes (|b| <= 1.3)
  follow the linear model to 0.01-0.04 (0.640/0.655, 0.203/0.210, 0.478/0.467); mean amplitudes up to |b| ~ 3 follow
  it; at |b| ~ 8 it breaks (0.342 -> 0.482, 0.557 -> 1.038); the 35-direction optimum (|b| 7-56) makes the chain
  diverge (NaN or ratio 4-8). **Per-network amplitudes must not be applied inside the chain**; the output-space form
  e + P(a delta) has no such limit.

### 5.2b The next-order chain K3 (the closest small-n analogue of production)

K3 (F1's dense third-cumulant chain: first-order Edgeworth mean, kappa_3 in the covariance, Gaussian kappa_3 of
relus, first-order kappa_3 transport; its error is one order up, MSE 10x below GC) with **D3 counterterms** D3:l
(scaling the mean readout's kappa_3 term, l = 3..7; the production D3 amplitudes), n32 L8, 6 networks
(`n1_k3ct.py`, `n1_k3ana.py`, `n1out/k3ana_6nets.txt`):

| net | X | merge | oracle | common | RC-G | **RC-Y** | projected, C_out top 10% (pi_U) | D3 shapes: p (corr) | noise s = 1: merge / projected |
|---|---|---|---|---|---|---|---|---|---|
| s0 | 0.81 | 0.144 | 0.522 | 1.106 | 1.43 | **0.560** | 0.486 (0.59) | 0.97 (0.81) | 0.77 / 0.54 |
| s1 | 0.78 | 0.193 | 0.332 | 0.818 | 3.07 | **0.375** | 0.548 (0.48) | 0.62 (0.81) | 0.97 / 0.63 |
| s10 | 0.89 | 0.150 | 0.509 | 0.815 | 12.6 | **0.580** | 0.638 (0.45) | 0.66 (0.76) | 1.63 / 0.76 |
| s11 | 0.80 | 0.245 | 0.883 | 1.057 | 2.48 | **0.947** | 0.720 (0.31) | 0.68 (0.85) | 0.80 / 0.77 |
| s2 | 0.41 | 0.442 | 0.498 | 0.839 | 1.66 | **0.639** | 0.596 (0.44) | 0.66 (0.71) | 1.01 / 0.65 |
| s3 | 0.68 | 0.324 | 0.730 | 2.914 | 1.58 | **0.744** | 0.958 (0.09) | 0.81 (0.91) | 0.94 / 1.04 |
| geo-mean | | 0.230 | 0.552 | 1.113 | 2.67 | **0.616** | 0.641 | | |

- **This is production's situation in miniature.** Truth-fitted common D3 amplitudes are useless held out (1.11): the
  optimal amplitudes are quenched. The per-network optimum of the same five amplitudes removes 45%; the truth-free
  RC-Y removes 38%, within 0.06 of it.
- **The grading mismatch is larger one order up.** The error's dimension is p_err = 2.4-2.6 (F1: 2.3-2.8), the D3
  shapes' 0.6-1.0: RC-G overshoots three- to fourfold and makes things 1.4-12.6x worse.
- **The free subspace** (top 10% of the K3 chain's own output covariance) holds 0.09-0.59 of the error and the
  projected merge removes 36% on average; under defect noise equal to the defect it still removes up to 46% (5 of 6
  networks gain; one loses 4%) where the full merge is break-even or harmful.

### 5.3 Production: the between-network spread of optimal amplitudes, from the refit data

**Estimator.** Per network k let b_k be its own optimal amplitudes (b_k = G_k^{-1} R_k^T eps_k), b_bar their mean,
Sigma_b their covariance, and assume G_k ~ G. A common fit on K training networks has expected in-sample gain
b_bar' G b_bar + tr(G Sigma_b)/K and held-out gain b_bar' G b_bar - tr(G Sigma_b)/K, while the per-network oracle has
b_bar' G b_bar + tr(G Sigma_b). Hence

    oracle per-network gain = held-out gain + gap/2 + (K/2) gap,      gap = in-sample - held-out.

*Check at small n* (`n1_gap.py`, random half splits): with 12 L8 networks of three widths (K = 6), predicted vs
measured oracle 0.82/0.68 (m3), 0.90/0.67 (m4), 0.88/0.58 (cN), 0.94/0.83 (var); with 16 networks of one width (n32,
K = 8, `gap_n32.txt`) 0.98/0.75, 0.87/0.73, 1.33/0.74, 1.00/0.85. The estimator is an overestimate, by 1.1-1.5x when
the held-out gain is comparable to the gap (m3, m4, var) and badly when the held-out gain is negative (cN). Production
is in the first regime (held-out 3.07, gap 2.97), so a factor of about 1.3 is the right correction there.

**Production numbers** (`notes/ray-compiler/outputs/calfit_103dirs.txt`, ridge 0, K = 50): in-sample -6.04%,
held-out -3.07% (measured -3.76%), gap 2.97 points:
- tr(G Sigma_b) = 25 x 2.97% = **74% of the MSE**; b_bar' G b_bar ~ 4.6%;
- per-network oracle with the same 103 amplitudes ~ **79%** (gap formula), **about 60% (45-70%)** after the small-n
  bias factor;
- the effective dimension this implies for each network's error is 103/(0.6-0.79) ~ 130-170 (if span(R) is generic
  within it), consistent with the participation ratio (53 at layer 15) and "top-128 modes carry 89% of the variance";
- by family (ridge > 0 for most, so these are lower bounds on the spread): coff 18%, D21 16%, D3 6.5%, g4 5.8%, K22
  6%, K31 5.3%; by layer band 0-4: 14%, 5-9: 31%, 10-15: 24%. The early layers' responses fit the per-network error
  as well as the late ones': the span works as a basis of the collective output subspace, not as a set of physically
  local corrections.
- **Coefficient space** (`notes/leaderboard-system/outputs/counterterms_*_{train0-49,all100}.json`; b_all - b_train =
  (b_B - b_A)/2 has covariance Sigma_b/100): per-statistic spread sd(b_k) ~ 10 x rms(b_all - b_train): var 1e-3
  (mean |b| 1e-4), coff 4e-3 (4e-4), D3 0.10 (0.011), D21 0.06 (0.021), g4 0.52 (0.030), K22 0.19 (0.016), K31 0.64
  (0.18). The per-network optimum moves the amplitudes by 3-17 times their fitted values (inflated by collinearity in
  coefficient space; the G-metric figure above is the one that matters). g4 and K31 per-network amplitudes of 1 +- 0.5
  are outside the linear regime, as at small n.
- The fitted var and coff amplitudes sit at 1 +- 1e-3 because the Gaussian sector is protected: the Gaussian closure is
  exactly RG-covariant along Gaussian directions (F1 sec. 6), so these directions carry little of the error and their
  responses have low dimension (p ~ 0.4 at small n).

### 5.4 Where the per-network error lives, and a free subspace

The output-space realization needs a subspace U. Candidates, small n (`n1_sub.py`, k = 7 = the m3 family size):
- share of the error (no counterterms) in the top-k modes of the chain's own output covariance C_out (free in the
  chain), mean over the 24 L8 networks (`sub_L8.txt`): 0.55 (k/n = 0.05), **0.66 (0.10)**, 0.81 (0.20), 0.88 (0.30),
  against k/n at random;
- per network at k = 7 (geometric mean, range): C_out 0.79 (0.46-0.99), the output first-chaos Gram J J^T 0.83
  (0.55-1.00), span(R_m3) 0.70 (0.41-0.96), random 0.16 (0.03-0.43);
- filtered merge e + a P_U delta at noise s = 0 / 1 / 2, geometric means: full 0.11 / 0.98 / 3.46; U = C_out top-7
  0.26 / 0.45 / 0.91; U = span(R_m3) 0.33 / 0.51 / 0.98.
- **The production-like pipeline** (`n1_final.py`, `n1out/final_m3_all.txt`): common truth-fitted m3 amplitudes first
  (the V47_CAL analogue), then the corrected chain's own defect delta_c = delta + Dd b_c projected on the top
  k = 0.1 n modes of C_out, one LOO scalar. Geometric-mean MSE ratios to the chain without counterterms:

| group | common | X after common | pi_C (k = 0.1 n) | s = 0: full / C_out / span(R) | s = 1 | s = 2 |
|---|---|---|---|---|---|---|
| n32 L8 (16) | 0.773 | 0.84 | 0.52 | 0.089 / 0.371 / 0.271 | 0.686 / 0.442 / 0.420 | 2.42 / 0.610 / 0.797 |
| n48 L8 (4) | 0.742 | 0.81 | 0.59 | 0.108 / 0.341 / 0.428 | 0.665 / 0.399 / 0.507 | 2.32 / 0.556 / 0.741 |
| n64 L8 (4) | 0.630 | 0.82 | 0.56 | 0.092 / 0.310 / 0.392 | 0.848 / 0.375 / 0.480 | 3.12 / 0.551 / 0.711 |
| n64 L16 (3) | 1.114 | 0.80 | 0.85 | 0.157 / 0.277 / 0.213 | 1.126 / 0.371 / 0.425 | 3.78 / 0.653 / 1.129 |

  The free subspace does as well as the counterterm span, the projected step halves the common-counterterm residual
  (common x (1 - pi_C X) = 0.630 x 0.54 = 0.34 predicted, 0.310 measured at n64 L8), and it survives defect noise
  that makes the full merge harmful. At depth 16 the residual error is concentrated in the collective modes
  (pi_C = 0.85 at k = 6 of 64): the free projection removes 72% of the error where the common counterterms remove
  none, and 63% at defect noise equal to the defect.
- **Production analogue** (DATA_FINDINGS): the chain's pre-activation error projected on the MC covariance's top
  16/64/128 modes holds 8-15/27-43/46-63% (random 1.6/6.2/12.5%). The free subspace should keep pi_U(128) ~ 0.45-0.65,
  somewhat less than the 0.5-0.79 of span(R) (section 5.3), at no cost and with no response runs.

**Answer to the brief's question ("is the projected local defect computable in budget, since it needs only
projections?").** Not more cheaply than the defect itself, for two reasons, and yes in a weaker sense for a third.
(a) The projection onto a counterterm shape at layer l is a neuron sum of the layer's local defect weighted by a
co-state, so the per-neuron local defect must still be formed at every neuron; projection saves only the transport,
which is cheap anyway. (b) Zeroing the projected defect is the wrong condition (RC-G, Thm 2.4(3)); the right one
(RC-Y) needs the error estimate -(1/2) Y^{-1} delta, then a fit in the span, which needs the per-network responses R
(103 JVPs, out of budget) unless the span is replaced by a free subspace (C_out top modes, as above). (c) What the
projection does buy is tolerance: incoherent errors of the defect estimate are cut by k/N, so a sketched analytic
defect (section 5.5) can replace the dense one.

### 5.5 Cost of the defect the condition needs

Every condition above needs the network's own delta, or its projection. Projection does not make the local defects
cheaper (they are per-neuron sums that must be formed at every neuron of the layer); it makes them **tolerant**: noise
that is incoherent across neurons is cut by k/N ~ 1/8, so sketched tangents are admissible. The in-budget defect of the
production chain, by Thm 3.3 and F1 Prop 7-8, consists of:
1. **GC-type diagonal terms** <grad mu_i, grad v_i> and |grad v_i|^2: these are 2 L^T H L and 4|H L|^2 in dressed
   form, i.e. the diagonal of V56's hub product Y and its Schur hub diag(Y M^{-1} Y^T). **0 extra units** (they are
   cancelled inside the Edgeworth readout by the carried D3 and path class at leading order, so what remains are the
   next items).
2. **Next-order open walks** in the Edgeworth readout's local defect (<grad v, grad D3> ~ L^T H^3 L, |grad D3|^2,
   <grad mu, grad g4>, all open): via the Schur identity, one product Z_s = (Y M^{-1}) At_s per late source-layer
   (1 unit each), or by a sketch: K random input directions transported as first-chaos legs (K n^2 per layer) and
   contracted with the source matrices (K n^2 per source-layer), K/n units per source-layer.
3. **Covariance-map local defects**, which reach the mean only through later layers' variance readouts (so only
   layers 12-14 matter in the late window):
   - diag-type: the first-chaos Gram L_a.L_b (needs the first-chaos map L_l, 1 unit per layer to transport from the
     input, 16 units, or K n^2 per layer sketched) and (H_a L_a).(H_b L_b), the off-diagonal Schur hub (one product
     after V56's solve, 1 unit per layer);
   - cross term sum_m A_am At_am w2_m P_bm (the first-chaos-projected (2,1) slice; At is already the dressed
     first-chaos cross-covariance): one hub-type product per source-layer plus its transport into the next variance,
     about 1 unit per source-layer, and it does not sketch (it needs no input directions).
4. **Wall terms** (triangle in D3, 4-cycle in g4, joint-gate closed walks in D21/K31): excluded (Thm 3.3(4)).

Late window (layers 12-15, 65% of the injected error), with AGE_OLD = 4, i.e. 5 young source-layers per layer plus the
factor-form old tier:

| item | dense | sketched (K = 64-256) |
|---|---|---|
| 1 GC-type diagonal | 0 | 0 |
| 2 next-order open walks (layers 12-15) | 20 | 2-9 (legs and contractions) |
| 3 diag-type (layers 12-14) | 16 + 3 | 1-4 + 3 |
| 3 cross term (layers 12-14) | 15 | 15 |
| projection on top-128 of C_15 (randomized subspace iteration), mean-gate transport | 1 | 1 |
| **total** | **about +55** | **about +22 to +32** |

E2 (section 7) reports the share of each item in the projected exact defect, so that the cheapest sufficient subset is
built; the cross term is the expensive item and, for GC, the dominant one at depth (F1 sec. 8e).

## 6. The estimator component and its prediction

**Component N1-RPD ("renormalized projected late defect").**

    e_RG = e + a_U . U U^T dinj_hat,     dinj_hat = sum_{l=12}^{15} J_{l->L} tr d_l^{open}(0, I),

with U the top-128 eigenvectors of the chain's last-layer covariance, J_{l->L} mean-gate transport, tr d_l^{open} the
in-budget local defect of section 5.5 (sketched or dense), and a_U one scalar fitted offline (truth-based, once; it
absorbs the single-degree approximation of Y^{-1}). Counterterms stay as they are (common constants), so the
component is the C+Y step of section 5.2 realized in output space.

**Prediction** (gain = s_late x pi_U x X_ib, minus sketch noise after projection):
- s_late = 0.65 (layers 12-15), pi_U = 0.5-0.65, X_ib = 0.30-0.55 (Cor. 3.4): gain 0.10-0.23;
- raw 1.55e-8 x (0.77-0.90) = 1.19-1.40e-8;
- sketched: C/B 0.203 + (22-32)/1024 = 0.224-0.234 (x 1.11-1.15): adjusted **2.7-3.3e-9, central 3.0e-9 (-5%)**,
  if the sketch noise after projection stays below half the projected defect;
- dense: C/B 0.257 (x 1.27): adjusted 3.1-3.6e-9, worse than today on most of the range;
- without the cross term (if E2 finds it small for the production chain): +7 to +17 units, adjusted 2.5-3.1e-9;
- window 9-15 (88% of the injection), sketched: gain 0.13-0.31 at about +45 to +60 units (x 1.22-1.29): adjusted
  2.6-3.5e-9.
- **Caveat on the sketch.** The k/N suppression of Prop 5.1 holds for noise incoherent across neurons. A sketch that
  uses the same K input directions for every neuron makes errors that are correlated along the collective modes of
  the tangents, and those survive the projection; they fall only as 1/sqrt(K) times |M|_F / |tr M| for the projected
  Gram M = sum_i lambda_i a_i b_i^T. K must therefore be set by E2 (K = 64-256 costs 2-9 units, so this is a fidelity
  question, not a budget one).

**Break-even**: s_late pi_U X_ib >= (C'/C - 1) / (C'/C) ~ 0.10-0.13 (sketched), 0.21 (dense), 0.03-0.08 (sketched,
no cross term). The sketched form pays on most of the predicted range, but by little; the dense form needs the top of
it; the variant without the cross term pays everywhere if its X_ib holds.

**What would change the verdict.** If E0 finds pi_U(128) >= 0.7, E1 finds X (projected) >= 0.6 with the
finite-difference defect, and E2 finds the open part carries >= 70% of the projected defect with a small cross term,
the late-window component goes to -15% to -25% adjusted and the 9-15 window to -20% to -30%. If E2 finds the wall
terms carry most of the projected defect (as Cor. 3.4 allows), it is under -5% and is not worth building.

## 7. The decisive experiments (AWS; the lead runs them)

**E0 (zero chain runs; minutes on one instance).**
- *Inputs.* The stored responses out_c{j}_{n}.npy (103 directions x 100 networks), out_b_{n}.npy, the V56 outputs
  with counterterms, truth_off{n}.npz, chaindump_{n}.npz (C at the last layer).
- *Outputs per network.* (a) pi_R = |P_{R_n} eps_n|^2/|eps_n|^2 for the 103-span and per family and layer band (exact
  per-network oracle, linear model); (b) pi_U(r) for U = top-r eigenvectors of the chain's last-layer covariance,
  r = 32, 64, 128, 256; (c) principal angles between span(R_n) and U; (d) the per-network optimal amplitudes and their
  spread (sd/|mean| per statistic) and the G-metric tr(G Sigma_b), to check the gap estimate of section 5.3.
- *Decision.* pi_R >= 0.4 and pi_U(128) >= 0.4: per-network renormalization has a ceiling of >= 40% of raw; go to
  E1. pi_R < 0.25: the gap estimate was wrong, close (iii) here.
- *Prediction (P1).* pi_R = 0.5-0.79 (central 0.62); pi_U(128) = 0.45-0.65; amplitude spread sd/|mean| >= 3 for
  D3, g4, K22, K31.

**E1 (uses the pending X* run: networks 0-3, 128 directions, float64, mask-frozen chain).**
- *Outputs.* RC-Y in output space: e + a P_U delta_hat for U = span(R_n), U = top-128 of C_15, U = all (the full
  merge), a fitted on networks 0-1 and applied to 2-3; MSE ratios. Repeat with delta_hat from 16, 32, 64 of the 128
  directions to measure the noise curve (Prop 5.1).
- *Decision.* Held-out projected ratio <= 0.85 with U = top-128 of C_15 at the full direction count: build E2.
- *Prediction (P2).* X (projected) = 0.35-0.65; projected ratio 0.55-0.80 with U = span(R), 0.6-0.85 with
  U = C_15 top-128; the projected ratio degrades less than the full merge as directions are dropped (at 16 directions:
  full merge >= 1.0, projected <= 0.95).

**E2 (network 0; the exact 2049-run defect with per-layer dumps, as F1 specifies).**
- *Outputs.* The in-budget analytic open-walk defect of section 5.5 at layers 12-15 against the exact per-layer
  defect, projected on U, item by item (next-order open walks, diag-type, cross term; dense and sketched at K = 64,
  128, 256), so that the cheapest sufficient subset is chosen; and, as a direct test of Thm 3.3(3)-(4), the triangle
  tr H_i^3 of one late layer computed offline at n^4 (about 1e12 FLOPs per source-layer) and its predicted defect
  contribution -(w3/2) tr H_i^3.
- *Decision.* Build N1-RPD iff the open part explains >= 50% of the projected exact defect.
- *Prediction (P3).* The open (in-budget) part explains 45-80% of the projected exact defect; the triangle and 4-cycle
  contributions, though a minority of the error, explain 20-55% of the defect.

## 8. Claims ledger

| claim | status | evidence |
|---|---|---|
| 2.1 localization = Polchinski flow; the chain = truncated effective action | proved (known) | Shi-Tian-Zhang; BBD; Bruned-Laubie-Minguella |
| 2.2 output = RG eigen-operator of eigenvalue 1/2; defect = anomalous dimension operator | proved (F1) | |
| 2.3 two-flow curvature: local defect = commutator of closure and RG; error = integrated curvature | proved | F1 Thm 1 + Prop 8 |
| 2.4 grading: delta = L(eps) = 2 Y eps; d delta/d b_j = 2 p_j R_j | proved; class-constancy measured | corr 0.84-0.98 |
| 3.1 joint cumulants = trace + vector-state moments of noncommuting (H_i) | proved, checked | MC 1e7 |
| 3.2 localization is a derivation on the walk algebra; Lap closes walks; Gamma concatenates | proved | |
| 3.3 exact split of the Polchinski hierarchy; residual of open-walk truncation = trace functional; visibility p = k for carried closed walks; no free defect | proved, checked | residual to 5 digits; ratio -> 3.0 |
| 3.4 corrected visibility bound; production X*_full 0.45-0.68, X*_ib 0.30-0.55 | proved given the class model; numbers conjectural | |
| 4.1 counterterm group non-abelian; linear refit = abelianization; pre-Lie second order explains the repo's 0.65-point refit pessimism | sketch / conjecture | |
| 4.2 locality at first order; prepared (sequential) conditions = joint, independent = overcounting | sketch; measured | 0.12-0.31 vs 0.09-0.19 vs 6-69 |
| 5.1 projected gain = pi X; noise crossover | proved; measured | +-0.05 |
| RC-G fails, RC-Y works (truth-free per-network amplitudes within 0.04-0.15 of the oracle; common amplitudes useless at depth 16 and on K3) | measured (GC n <= 64, L = 8, 16; K3 n = 32) | `ana_all.txt`, `k3ana_6nets.txt` |
| output-space projection on the chain's own top covariance modes: free, as good as the counterterm span, noise-robust | measured | `final_m3_all.txt`, `sub_L8.txt` |
| production per-network spread tr(G Sigma_b) ~ 74% of MSE; oracle about 60% (45-79%) | estimated from refit data (estimator checked at small n, overestimates 1.1-1.5x) | E0 decides |
| per-network amplitudes break linearity inside the chain | measured | `nonlin.txt` |

**Novelty, honestly.** Theorem 2.1 is known; 2.2-2.4 are F1 with an RG reading (the Connes-Kreimer beta = Y Res
analogy is exact as algebra, and is what makes RC-Y the right condition). The walk-algebra split (Thm 3.3) is, to my
knowledge, new as stated; it is a short computation, and its value is that it locates the heat defect's expensive part
precisely on the trace-state moments and gives the visibility weights. The per-network finding (section 5.3) is a
re-reading of existing refit numbers and needs E0 to be trusted. The renormalization conditions are classical
Galerkin/defect-correction ideas; the measured failure of the naive condition and its grading fix are the new part.

**What the frame does not give.** A cheaper defect. The heat defect of this chain, by Thm 3.3, contains the
closed-walk wall wherever the chain carries a table that omits closed walks; the in-budget part is the open-walk part,
and projection only buys noise tolerance. The per-network ceiling (about 60%) is large, but reaching it needs a
truth-free detector of the collective error coordinates, and the only one on the table is the defect, capped by
X_ib.

## Referee report

Adversarial referee for frame N1. My checks are independent of `n1code/` and live in `loc/ref_n1/`:
- `split_check.py`: Prop 3.1, Lemma 3.2 and Thm 3.3(1)-(2) by finite differences, for commuting and for
  noncommuting kernels;
- `homog_triangle.py`: the triangle's visibility in a *1-homogeneous* relu net, with output in
  `homog_triangle_out.txt`;
- `cor34_scan.py`: a scan of Corollary 3.4's parameter ranges.

I also read the refit data the frame's section 5.3 rests on (`notes/ray-compiler/outputs/calfit_103dirs.txt`), note
XLI sections 1-2, and the F1 referee report. I confirmed that Bruned-Laubie-Minguella (arXiv 2608.31049, 31 Aug 2026,
"On the equivalence between the Polchinski flow and the Connes-Kreimer approaches to perturbative renormalisation")
exists.

### Q1. What is correct (verified)

- **Thm 2.1** is correct and known. It is Price plus Hopf-Cole. The Polchinski reading is standard
  (Bauerschmidt-Bodineau-Dagallier, Prob. Surveys 2024; Shi-Tian-Zhang 2510.04460).
- **Prop 3.1** is correct. I re-derived the multilinear coefficients of -(1/2) log det(I - sH_t) + (s/2) L'^T (I - sH_t)^{-1} L'.
  It is the classical joint-cumulant formula for Gaussian quadratic forms (Magnus; Mathai-Provost), written as words.
- **Lemma 3.2** is correct. I checked both items by hand.
- **Thm 3.3(1)-(2)** is correct. With n = 6, three neurons, m != 0 and s = 0.8:
  - the full Polchinski residual is about 4e-7 (finite-difference level);
  - the open-truncation residual is 0.751606295, against -3 s^2 tau = 0.751606295.
- **Thm 3.3(3), visibility p = k, survives in the homogeneous setting.** This is the most useful thing I can add. The
  author's check uses a non-homogeneous second-chaos toy, in which H is fixed along the path. The production network
  is 1-homogeneous, so its second-chaos kernel Hbar(y) = E grad^2 F(y+g) *moves* with the mean.
  - Rung 2 of the Euler-Stein ladder in matrix form (tr_(3,4) c_4 = -c_2) gives Lap_y Hbar(0) = -Hbar(0).
  - grad Hbar(0) = c_3, the third chaos.
  - For tau = tr Hbar^3 this gives Lap tau = -3 tau + 6 Gamma, with Gamma = sum_e tr(c3[e] c3[e] c2). Homogeneity
    gives tr d_Sigma tau = (3/2) tau at the base point.
  - Hence tr H tau = 3 tau - 3 Gamma, i.e. **p_eff = 3(1 - Gamma/tau)**. The second-chaos value 3 is exact when the
    third chaos vanishes, and the correction is a third-chaos-dressed closed walk.
  - Measured (Hermite Monte Carlo, 2.4e7 samples):

    | net | p_eff | |c3|/|c2| |
    |---|---|---|
    | one hidden layer (n = 5) | 3.00 for every neuron | at the MC noise floor (third chaos absent) |
    | three hidden layers (n = 5) | 2.66-2.97 | 0.2-0.37 |

    The weight-3 claim is therefore sound in the homogeneous problem, with a depth-growing correction that the frame
    should state.
- **Thm 2.4** is correct but elementary: L(G) is linear, and d/ds s^p = p at s = 1. **Prop 5.1** is correct. I
  re-derived (2X - X^2) pi_e - X(1 - X) pi_w. One fix is needed: the crossover condition is on a^2 N sigma^2, not on
  N sigma^2.
- **The refit-gap estimator is algebraically correct for equal G_k.** Pooled OLS is then the mean of the per-network
  projections. I re-derived oracle = held-out + gap/2 + (K/2) gap and reproduced 74% / 79% from the file (train -6.04%,
  held-out -3.07% +- 0.38, K = 50).
- **The arithmetic of section 6 checks out** given its inputs: gain 0.10-0.23, adjusted 2.7-3.3e-9 sketched, break-even
  0.10-0.13.

### Q2. Noncommutative content: real, but it changes no computation in this frame

1. **Every theorem of section 3 holds verbatim for commuting H_i.** `split_check.py` with simultaneously diagonal
   kernels gives:
   - an open-truncation residual equal to -3 s^2 tau to 9 digits;
   - a full residual of 7e-8.

   The split is by *degree in L'*, not by any noncommutative structure. The visibility weights depend only on the
   scaling of a closed walk. Noncommutativity enters only the *cost*: noncommuting H_i cannot be simultaneously
   diagonalized, so the traces of words do not reduce to n^2 contractions of a spectrum. That is the n^4 wall that notes
   XL and XLI already priced. The "*-distribution under vector states and the trace" language is accurate terminology
   for a classical quadratic-form computation.
2. **"Noncommutative Dirichlet form in the Cipriani-Sauvageot sense" is a mislabel.**
   - The semigroup is the classical Gaussian heat semigroup on functions of m in R^n. It acts on polynomials with
     matrix coefficients.
   - The derivation moves the vectors (L' -> L' + He) and never acts on the algebra generated by the H_i.
   - The trace appears because sum_e (H e)^T W (H e) sums over an orthonormal basis, not because of a trace on a von
     Neumann algebra that the semigroup preserves.
   - No completely positive Markov semigroup on A = alg(H_i) is defined or used.
3. **The Hopf / BPHZ / Chen-series layer is relabelling.**
   - The parameter space of the counterterms is additive R^103. The "non-abelian group" is the Taylor expansion of a
     composition of maps (Chen-Fliess), and "abelianization" just means linearization. No coproduct, antipode or
     Birkhoff factorization is computed.
   - The grading "Y" is multiplication by an empirical profile exponent on an ad hoc class decomposition, not a Hopf
     grading. "beta = Y Res" is an analogy of notation.
   - "BPHZ-prepared sequential conditions versus block-Jacobi overcounting" is the textbook behaviour of Gauss-Seidel /
     forward stagewise regression against Jacobi on a non-diagonally-dominant normal matrix with overlapping regressor
     spans. The 6-69x blow-up is expected for Jacobi and says nothing about forests.
   - The proof sketch of Thm 4.2 is also unsound. The later local defects d_l' change at the same order as d_l. Both
     are (closure non-Gaussianity) x b, and a closure's d(Phi; T) does not vanish identically in T. "Locality at first
     order" is therefore not established; only the measurement stands.
4. **RC-Y is not a new renormalization condition.** b = argmin |a delta - R b|^2 gives R b = P_R(a delta). RC-Y is
   exactly F1's merge projected onto span(R). The section 5 tables confirm it:
   - noiseless, the full merge beats RC-Y everywhere (common + full merge 0.098 against common + RC-Y 0.297 for m3);
   - the projection pays only as a denoiser.

   So "RC-Y within 0.04-0.15 of the oracle" is the identity pi X ~ pi at X ~ 0.8. "RC-G fails" is the correct and
   useful negative: do not zero the defect, regress the merge.

### Q3. Mathematical errors and overclaims

1. **The exact split holds only in the second-chaos truncation, and the production chain is not near it.**
   - DATA_FINDINGS gives the first-chaos share as 0.215 at layer 15, so 78% of the variance is chaos >= 2, of unknown
     third-chaos content.
   - Beyond the second chaos, grad tau = c_3-terms != 0. The residual of an open-walk truncation is then not a pure
     trace functional; it also carries third-chaos-dressed terms (Gamma above).
   - "Proved" should read "proved in the second-chaos sector; corrections O(c_3^2) measured 1-11% at depth 4, n = 5".
2. **"No in-budget form of the defect contains this part" is too strong.**
   - Exactly, the claim is right: per-neuron words need n^4.
   - A randomized trace is still possible. With K shared Gaussian probes, H_i g = Lambda(X_i o Lambda^T g) is one n x n
     product per probe and source-layer, and the second application goes through the Schur Gram At^T C^{-1} At. So
     tr H_i^3 for all i costs about 2-4 K units per source-layer.
   - Restricting to the top-r eigenspace of the birth Gram costs about r^2/n units per source-layer.
   - Both are unbiased or controlled approximations whose variance or truncation error is unmeasured. The correct
     statement is "not affordable at a variance that is known to be useful".
   - Note XLI's warning that the hub's quadratic form lives in the covariance's bottom directions makes the low-rank
     route doubtful, but it is not excluded.
3. **Corollary 3.4 is mis-summarized, and its main inference runs the wrong way.**
   - With a single open-class weight, X*_ib = (o p_o)^2 / (o p_o^2) equals the open share o exactly.
   - Scanning the frame's own ranges (`cor34_scan.py`, 126 mixes) gives:

     | quantity | range | median |
     |---|---|---|
     | X*_full | 0.45-0.77 | 0.63 |
     | X*_ib | 0.35-0.60 | 0.50 |
     | X*_ib / X*_full | **0.52-1.25** | 0.78 |
     | wall share of the defect energy | 0.36-0.94 | 0.75 |

   - A ratio above 1 means that dropping the high-p wall terms can *raise* the single-slope explained share: weight
     3-4 components distort the one-coefficient merge.
   - So "kappa_a^2 = capture <= 1, the wall is a loss" is wrong as stated. For a one-slope merge, the in-budget defect
     is not dominated by the full one. The real consequence of the weights is different: **the finite-difference X* of
     the pending run is not an upper bound for the analytic merge, and a two-slope merge (open and closed fitted
     separately) is the right comparison.**
4. **The per-network ceiling (section 5.3) is softer than presented.**
   - The gap is one 50/50 split with held-out se 0.38, so tr(G Sigma_b) = 74 +- about 12%.
   - G_k = G is assumed. Heterogeneous G_k inflates the gap without any amplitude spread.
   - The oracle includes the chance-level projection of incoherent error on a 103-dim span (about 103/1024 = 10% for
     isotropic error), which no truth-free detector can reach.
   - The linear-response oracle is unattainable inside the chain (the frame's own nonlinearity result), and in output
     space R costs 103 JVPs.
   - "Per-network adaptivity is worth 60%" is therefore a diagnostic of how concentrated the error is, not a reachable
     target. E0 measures the honest version of it.
5. **Prop 4.1's explanation of the refit's 0.6-0.7 point pessimism by "pre-Lie insertions" is untested.** Second-order
   response and the shrinkage formula's own conservatism are simpler alternatives.

### Q4. Cost-accounting errors

1. **Transport.**
   - The covariance-map local defects of layers 12-14 reach the mean through later variance readouts, which requires
     a tangent-linear pass through the covariance channel. The F1 referee (R3.1) priced this at about 2-3 units per
     layer.
   - The 1-unit "projection plus mean-gate transport" line undercounts by about +8 to +12 units.
   - Mean-gate-only transport also has unmeasured fidelity (F1 R2.5).
2. **Item 2.**
   - |grad D3|^2 is the open walk L^T H^4 L. It needs a second hub product (q At^T, then the solve) per source-layer
     beyond (Y M^{-1}) At_s, so the dense figure is about 40 units, not 20.
   - The old factor-form tier is not costed.
3. **The in-budget defect is defined for the second-chaos model, not for V56's actual readouts.** V56's readouts are
   Edgeworth readouts with D3, g4, D21, K22 and K31 plus counterterms, a regularized (C + eps)^{-1} hub, and saturation
   drops. F1's "one level up" risk is inherited unresolved.
4. **Corrected costs.**
   - Sketched: +32 to +50 units, C/B 0.234-0.252, cost factor 1.15-1.24. Break-even s pi X >= 0.13-0.19.
   - At the frame's central gain of 0.165: raw 1.29e-8 x 0.243 = **3.15e-9, i.e. break-even with today's 3.14e-9**.
   - Range: 2.8-3.5e-9, before any sketch-noise penalty.
   - The dense variant loses, as the author says.

### Q5. Relevance and experiment design

1. **pi_U is borrowed from the wrong object.** The 46-63% in the top-128 modes in DATA_FINDINGS is for the
   *transported pre-activation* error W_l e_(l-1). It excludes the layer's own readout injection (21% of the final MSE
   is injected at layer 15). The output error of the late injections, which is the component's target, may be less
   collective.
   - E0 should report pi_U for (a) the total output error and (b) the late-injected part, using the 8-network injection
     decomposition. It should also use the output (post-relu) covariance as well as C_15.
2. **E1 and E2 are exposed to the known late-layer kink.** The smoothness pilot found a ~1e-7 kink from layer ~9-10,
   suspected to come from the nested tier-2 compression of old sources, which dominates second differences. The
   frame's experiments use finite-difference defects at layers 12-15, exactly where the kink lives.
   - Freeze the tier-2 compression bases and ranks along with the masks.
   - Re-run the smoothness pilot at layers 12-15 before E1/E2.
   - Otherwise the projected delta_hat measures the kink.
3. **E1 is underpowered.** a is fitted on 2 networks and tested on 2, against a 0.85 threshold. Use leave-one-out
   over the 4 networks, report each network, and bootstrap over directions.
4. **E2 should fit a two-slope merge** (open and closed parts with separate coefficients), as Q3.3 explains, and
   report Gamma/tau for one late layer to test the homogeneous correction to p = 3.
5. **E0 is the best item in the frame:** zero chain runs, cheap, and decisive for pi_U. Run it first.

### Q6. Novelty

- Thm 2.1 is known.
- The RG reading of F1 is known in substance: Polchinski equals localization, and BLM 2608.31049 connects it to CK.
- Prop 3.1 is classical.
- The degree split of Thm 3.3 is a short new observation. It is commutative mathematics.
- The visibility weight p = k is a correct, useful refinement of the synthesis's weight-0 table. I have extended it
  above to the homogeneous network via rung 2.
- The renormalization-condition material is defect correction plus Galerkin projection plus Gauss-Seidel.
- New and useful empirically:
  - RC-G (zeroing the defect) fails because counterterm shapes have different profile exponents from the error;
  - projecting the merge onto the chain's own top covariance modes is a strong denoiser at small n.

### Q7. Verdict

- **Mathematics:** sound where stated in the second-chaos sector, with one overclaim each in Thm 3.3(4) and Thm 4.2,
  and a mis-summarized Cor 3.4.
- **Noncommutative / Hopf content:** vocabulary. Noncommutativity explains a cost already known and changes no formula
  here. The Hopf layer changes no computation.
- **Estimator N1-RPD:** after the cost corrections, its central prediction is break-even (about 3.15e-9 against 3.14e-9).
- **Strongest surviving idea (corrected).** A truth-free late-window merge e + a P_U delta_open:
  - U is the top modes of the chain's own output covariance, used as a noise filter for a sketched open-walk defect;
  - the omitted closed walks are handled by their separate visibility weight p ~ 3(1 - Gamma/tau), not folded into
    one slope;
  - it is gated by E0 (pi_U of the late-injected output error >= 0.4) and by a kink-free, mask- and basis-frozen E1.

Priority: moderate-low. Run E0, because it is free.
