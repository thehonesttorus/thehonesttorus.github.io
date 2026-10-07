# The chain's fourth-cumulant state as one declared tensor

Working note XXXII. A construction, not an experiment. It answers the three obligations set after the K4 theory reset
(THEORY_3, "coherent tensor transport, metric dependence, and non-overlapping residuals", and the synthesis that
followed it), in the coordinates of the production code `est_v29.py` (V33_K4Q = 3, the regenerated core with the
adaptive lambda rule). The obligations:

1. specify one coherent corrected fourth-order tensor and derive every slice the chain uses from it;
2. separate representation error (misuse of retained information) from law-closure error (information the retained
   state does not determine);
3. define the complete terminal comparison, with the signed interactions kept and the radial null respected.

What is new relative to THEORY_3: the code correspondence it could not check is checked here (section 1, exact algebra
in `code/check_k4_contract.py`); the slices of the exactly transported retained tensor are given in closed form with
their costs (section 2); the representation / law-closure split is made entry by entry for what the chain carries
(section 3); and the radial null is written in the chain's own state variables, with a proof that the chain's mean map
already satisfies it, so that any violation is located in the state map (section 4, `code/check_radial.py`). No network
experiment was run for this note. Section 5 records measurements that already existed when it was written; they are not
tests of the constructions here.

## 1. What the chain retains and what it declares (the code contract, checked)

Indices a, b run over the post-activation y of layer l; i, j over the next pre-activation z = W y. H = W o W,
r = H 1 = diag(W W^T), P(v) = cA v + cI (1^T v) 1 with cA = 6/(n+4), cI = -3/((n+2)(n+4)) (`st["cA"]`, `st["cI"]`),
METRIC_C = 2.

**Retained at the y level** (the Wick program's outputs; post-activation assembly and regenerated block of
`notes/birth-address/code/est_v37_audit.py`): the fourth-cumulant diagonal
d = K4v and the (2,2) slice K = sym(K22), zero diagonal. Nothing else of the y fourth cumulant is computed: the pair
program has outputs (1,1), (2,1), (2,2) only, so the y (3,1) slice and the (2,1,1), (1,1,1,1) classes are absent. As a
tensor the retained y state is

    T_y = D(d) + S22(K),     polynomial  T_y[v^4] = (v o v)^T (3K + diag d)(v o v).

**Declared at the z level** (the regenerated block with K4Q = 3), with u = d + K 1, Khat the Ritz approximation of K (K4Q_RANK + 8 = 12
pairs, not 4), Khat0 its zero-diagonal part, E = K - Khat0, s_off^2 = var - H var_y, and lam the adaptive per-layer
scalar:

    g4row = Q(K, d) - R_W(E) + 2 lam s_off^2,      Q = (W o W o W o W) d + 3 diag(H K H^T),
                                                     R_W(E) = 3 diag(H E H^T) - 2 H P(E 1),
    wk4m  = (dG_i + dG_j)/3,                         dG = H P(u) + lam s_off^2,
    wk431 = lam C_off(z),                            wk431[a, c] = kappa(z_a, z_c, z_c, z_c).

The first line is the exact identity THEORY_3 assumed (checked to 4e-15 against the code's own operation sequence), so
its three-way split R_W(E) = D_metric + D_known + D_unresolved applies verbatim. The other two lines have a precise
reading. Let N = W diag(P(u)) W^T + lam W C_y,off W^T. Then:

- wk4m is exactly the (2,2) slice of the product core J(2I, N), whose slices are K_ij = (M_ii N_jj + M_jj N_ii +
  4 M_ij N_ij)/6 at M = 2I;
- the diagonal of the same core, 2 diag N, is what g4row would be without the K4Q correction: g4row adds to it the
  exact-pair content of Khat0 and of d (the Ritz part of the pair class, with the realised metric) that the other two
  slices never receive;
- that core's (3,1) slice is N_ij; the code's wk431 = lam (W C_y W^T)_ij keeps the lam constituent (with C_y in place of
  C_y,off) and drops the pair core's part (W diag(P(u)) W^T)_ij entirely.

So the three declared slices come from three different constituents: the diagonal from the exact pair class (up to
R_W(E)) plus lam, the (2,2) slice from the Euclidean trace core at metric 2I plus lam, the (3,1) slice from lam alone.
They are non-overlapping classes and do form a symmetric tensor; what fails is that this tensor is not the transport of
any declared y-level tensor. A coherent update of all three is therefore a new tensor model, not a more accurate
evaluation of the present one.

## 2. Obligation 1: one coherent tensor

**Declaration.** The coherent corrected z-level tensor is the exact transport of the retained y tensor plus one
explicitly declared constituent for the y classes the chain does not retain:

    T_z = W#[ D(d) + S22(K) ] + W#Lambda_y,

with every slice the next stage reads (diagonal, (2,2), (3,1)) derived from T_z and the realised W (no 2I metric). For the
present lam term the declared constituent is Lambda_y = J(I, lam C_y,off), whose transport J(W W^T, lam W C_y,off W^T) at
metric 2I reproduces the code's lam contributions to g4row and wk4m exactly. Which Lambda_y is right is a law-closure
question (section 3); the coherence requirement is only that it be one declared tensor.

**Slices of the transported retained tensor** (`code/check_k4_contract.py`, against an explicit 4-tensor), with
A = 3K + diag(d) and u_ij = W_i o W_j (rows i, j of W):

    diag_i      = Q_i = (W^o4 d)_i + 3 diag(H K H^T)_i
    (2,2)_ij    = (H K H^T)_ij + (H diag(d) H^T)_ij + 2 u_ij^T K u_ij
    (3,1)_ij    = kappa(i,i,i,j) = [((H A) o W) W^T]_ij

All three come from one product HK: the diagonal is rowsum((HK) o H) + W^o4 d, the (2,2) main term (HK) H^T, the (3,1)
slice ((3HK + H diag d) o W) W^T. Exact cost: three dense n^3 products per layer (about 1.5% of B each over 15 layers, 4.4% of B
together, before Strassen), against one for the diagonal alone (note XXVIII's mode 1). The cross term 2 u_ij^T K u_ij has zero mean
over an independent He row when K has zero diagonal and second moment of order (2/n)^4 ||K||_F^2; against the
coherent main term (2/n)^2 1^T K 1 it is of relative size 2 sqrt(2) ||K||_F / 1^T K 1, about 3/n for a mixture-shaped
K. Its omission must be declared, as a tensor (the (2,2)-class entries 2 u_ij^T K u_ij subtracted), not left implicit.

**The three-way decomposition as the cheap tier.** With K = Khat0 + A_g + E_0 (THEORY_3 section 4: A_g the
minimum-norm table with the row sums s = E 1, E_0 1 = 0), the diagonal splits as THEORY_3 states, and the accessible
part R_W(A_g) = 6 (Hg) o r - 6 W^o4 g - 2 H p is O(n^2). The other slices do not inherit that price: the A_g part of
the (3,1) slice is 3 [(Hg)_i (W W^T)_ij + r_i (W D_g W^T)_ij - 2 (W^o3 D_g W^T)_ij] and the Khat0 part is
3 sum_k lam_k (HV)_ik (W D_(v_k) W^T)_ij up to its diagonal correction, both dense products. So a coherent tensor is cheap only on the diagonal; on
all three slices the dense product HK is the cheaper exact route, and the three-way split is an analysis of the
diagonal, not a cost saving for the full state. Two consequences:

- D_unresolved is not unknown information in this code: K is a stored n x n matrix when it is transported. It is the
  part whose transport costs n^3; THEORY_3's "unresolved" is relative to the row sums, not to the state.
- A diagonal-only repair (adding D_known + D_metric) moves the diagonal further from the other two slices. Whether that
  helps is a terminal question (section 4), and the theory predicts no sign for it.

**Why a product-core state in the covariance metric is not a shortcut here.** THEORY_3's moving-metric decomposition
with M = C (the covariance, transported for free since C_z = W C_y W^T) would carry the trace core covariantly and give
all slices at O(n^2). But its core needs the full trace U_ij = sum_a T_ijaa, and the off-diagonal trace is
U_ij = B_ij + B_ji + sum_(a not in {i,j}) T_aaij (THEORY_3 section 7): the (3,1) and (2,1,1) classes, which the pair state
does not retain. The pair-retained state determines only diag U = u, which is exactly why the code's core is
J(I, diag P(u)). A covariance-metric core therefore needs law-closure input (for the scale mixture, U_ij is proportional
to C_ij); it is a section 3 choice, not a representation repair.

## 3. Obligation 2: representation error and law-closure error

For one transition l -> l+1, with the z_l state taken as given:

| quantity the next stage reads | determined by the retained state? | in the code now | kind of error |
|---|---|---|---|
| y pair classes d, K | yes, through the bivariate closure of the Wick program (Hermite order MEHLER = 2) | computed | law closure of the pair program (its truncation) |
| transport of D(d) + S22(K) to the z slices | yes, deterministic algebra in K, d, W | diagonal up to R_W(E); (2,2) by the 2I trace core; (3,1) not at all | representation |
| y (3,1) slice and its transport | yes: it is a pair statistic of y, fixed by the z pair state under the same bivariate closure | not computed (no (3,1) output) | representation (closure-determined, uncomputed) |
| y (2,1,1) and (1,1,1,1) classes | no: they need triple and quadruple joint information; the triples are partly in the kappa_3 legs (all-distinct entries), the quadruples nowhere | the lam C_off constituent, lam fitted per layer | law closure |
| kappa_5 and higher in the mean map and pair program | no | dropped | law closure |
| errors already in the z_l state | n/a | propagated | input |

The synthesis's identifiability argument locates the law-closure rows exactly. Take as retained features at the y
level the pair-restricted polynomials of degree at most four (and, through the legs, the cubic monomials), and as F any
quantity the next stage reads, for instance kappa(z_i, z_i, z_j, z_k) contracted with W. If the residual r_F of F after
L2 projection on the retained features is nonzero, the laws (1 + t r_F) nu and (1 - t r_F) nu agree on everything
retained and disagree on F. The pair program resolves that freedom by a specific choice: the bivariate Edgeworth law of
the retained pair cumulants, with the non-pair classes set by the lam constituent. Representation error is
everything in the table that the retained state does fix; law-closure error is what that choice gets wrong.

Two consequences for the design:

- Restoring the transport of the retained tensor (section 2) and adding the y (3,1) output to the pair program are pure
  representation repairs: no new observable, no new assumption.
- The (2,1,1) and (1,1,1,1) classes need either new observables or a declared law class. The theory documents give the
  candidates: the scale mixture (closed form 6 g s_diag^2 s_off^2 + 3 g s_off^4 on the diagonal, note XXI), and the
  facet response of the incoming kappa_3 sources (Delta_T, the third-chaos matrix element of the cumulant theory),
  which the legs could feed. Each defines a Lambda_y for section 2.

## 4. Obligation 3: the terminal comparison and its null conditions

**The comparison.** With e the baseline chain's terminal error vector and Delta_X the terminal change produced by
model X propagated through the complete chain (not a sum of per-mechanism changes),

    MSE(X) - MSE(0) = (2 <e, Delta_X> + ||Delta_X||^2) / n,

and for two mechanisms the interaction MSE(a+b) - MSE(a) - MSE(b) + MSE(0) = 2 <Delta_a, Delta_b>/n plus the nonlinear
propagation terms. The comparison is therefore a factorial over the mechanisms, not a list of single changes:
D_metric, D_known, D_unresolved on the diagonal; the (2,2) and (3,1) slices of the transported tensor; the y (3,1)
output. Three controls fix what is being compared.

- **Freeze lam first.** The adaptive rule rescales lam each layer by the ratio of mean(dG) to mean(var) against a fitted
  table (REF_R). A change of the declared diagonal moves that mean and so changes lam, which changes all three slices.
  The first comparison keeps the baseline lam trajectory (recorded per layer) fixed; releasing it is a second,
  separate comparison.
- **Ceilings.** The oracles of section 5 give the terminal effect of the true diagonal, the true sector, and the true
  sector with true kappa_3 readouts. A representation repair cannot exceed the sector ceiling unless it also corrects
  law closure.
- **A first-order prediction before each run.** The adjoint form e_out = sum_l v_(l+1)^T delta_l (the cumulant theory,
  section 6; `notes/birth-address/code/errbudget.py`) predicts <e, Delta_X> from per-layer defects. It has to be
  validated on the existing D3 and g4 oracles before it is used to rank anything.

**The radial null, in the chain's variables** (`code/check_radial.py`). Let S = 1 + eps be independent of the input
direction, with E eps = 0 and E eps^2 = v. Positive homogeneity gives E F(S A) = E F(A) for the terminal mean: the
mean-preserving radial deformation changes every cumulant and no terminal mean. For a state (mu, Sigma, kappa3, kappa4)
its tangent per unit v is

    d mu = 0,   d Sigma = Sigma + mu mu^T,   d kappa3 = 3 kappa3 + 2 Sym3(mu, Sigma),
    d kappa4 = 6 kappa4 + 4 Sym3(Sigma, Sigma) + 3 Sym4(mu, kappa3)

(checked by exact enumeration on a non-Gaussian finite law; the residual is O(v) and halves with v). In the slices
the chain carries:

    diag:   6 d_i + 12 var_i^2 + 12 mu_i D3_i
    (2,2):  6 K_ij + 4 (var_i var_j + 2 C_ij^2) + 6 (mu_i D21_ji + mu_j D21_ij)
    (3,1):  6 B_ij + 12 var_i C_ij + 9 mu_i D21_ij + 3 mu_j D3_i

and the exact layer map sends the tangent of the z_l state to the tangent of the z_(l+1) state, because ReLU and the
linear layer commute with an independent positive scale. The 4 Sym3(Sigma, Sigma) term is the scale mixture; the
3 Sym4(mu, kappa3) term is THEORY_3's signed radial coupling, with coefficient m4 - m1 m3 = 3v.

**The mean map already satisfies it.** At a Gaussian point the Gram-Charlier mean
mu Phi + sigma phi - kappa3 alpha phi/(6 sigma^2) + kappa4 (alpha^2 - 1) phi/(24 sigma^3) has first-order response along
the tangent

    phi [ (sigma^2 + mu^2)/(2 sigma) - mu alpha + sigma (alpha^2 - 1)/2 ] = sigma phi [ (1 + alpha^2)/2 - alpha^2 + (alpha^2 - 1)/2 ] = 0

identically in alpha (to 3e-17 at four test points). The variance response and the kappa3 and kappa4 responses cancel
exactly. So the chain's handling of the radial direction is decided entirely by its state map: whether the
regenerated sector carries 4 Sym3(Sigma, Sigma) + 3 Sym4(mu, kappa3) + 6 kappa4 forward as the exact map does. Two
tests follow, neither needing Monte Carlo:

- **Terminal null.** Inject the input tangent (at the first layer mu = 0 and kappa3 = 0, so it is d Sigma = Sigma and
  d kappa4 = 4 Sym3(Sigma, Sigma)), propagate through the chain, and require zero first-order change of the output. The
  covariance part alone moves every output; the size of the chain's residual against that is its homogeneity defect.
- **Layerwise Ward test.** Push the tangent of the chain's state at z_l through one chain layer and compare with the
  tangent of its state at z_(l+1), slice by slice. This locates the defect, and it applies to every candidate
  T_z of section 2 before any terminal run. Truncation at fourth order allows a residual at kappa_5 order, so the test
  separates large violations from truncation-level ones.

The adaptive lam rule is not homogeneous: it compares mean(dG), of degree four, with mean(var), of degree two,
against a fixed table, so it fixes an absolute scale of the fourth-cumulant sector. That is harmless for the fixed
N(0, I) input but it is a structural reason to expect a Ward defect in exactly the place the oracles point to.

A third, necessary condition applies to any declared state: per pair (i, j), the 3 x 3 Gram matrix of the quadratic
features z_i^2, z_i z_j, z_j^2 after subtracting their best linear predictors (THEORY_3 section 8, with the pair
kappa_3 and kappa_4 slices) must be positive semidefinite. It costs O(n^2) and is a realisability check of the bivariate
laws the pair program implicitly uses, not a sufficiency condition.

## 5. Evidence already in hand (measured before this note)

These come from the oracle runs of note XXXI (`notes/birth-address/outputs/oracle_mc2_off1.txt`,
`oracle_cmp_mc2_off1.txt`, `onestep_mc2_off1.txt`, `oracle_mc3_off0.txt`; Monte Carlo of 1.6e7 inputs, noise-free
extrapolation from two halves). They are not tests of any construction above, but they bear on all three obligations.
Network 1 has the full set; network 0 has an independent replica of the combined oracles.

The chain's declared slices against Monte Carlo truth, per layer (MC noise in brackets):

| layer | diagonal: rel err, corr | (2,2) slice | (3,1) slice | C_off | var |
|---|---|---|---|---|---|
| 3 | 0.125, 0.74 (0.063) | 0.144, 0.57 (0.075) | 0.81, 0.63 (0.40) | 0.005 | < 0.001 |
| 7 | 0.164, 0.79 (0.034) | 0.178, 0.65 (0.040) | 0.74, 0.69 (0.16) | 0.005 | 0.001 |
| 11 | 0.265, 0.74 (0.025) | 0.257, 0.62 (0.031) | 0.76, 0.69 (0.09) | 0.007 | 0.001 |
| 14 | 0.318, 0.76 (0.022) | 0.325, 0.67 (0.026) | 0.75, 0.72 (0.06) | 0.008 | 0.001 |

**The defect is regenerated each step, not inherited.** In a run with every pre-activation statistic replaced by truth
at every layer, the chain's own next-layer slices, computed from true inputs before the replacement (inputs from one
Monte Carlo half, judged against the other), are as wrong as in the free-running chain: diagonal 14-33% (corr
0.64-0.80), (2,2) slice 16-32% (corr 0.47-0.70), (3,1) slice 74-90%, at layers 2-15. Input propagation contributes
almost nothing; the state map produces the error afresh at every layer. This is the representation-plus-closure
quantity that section 6, step 1 splits.

Terminal effect of replacing statistics by truth at every layer (raw 2.1266e-8 for the chain on network 1):

| replaced | change |
|---|---|
| (2,2) slice alone | **+12.7%** (worse) |
| (3,1) slice alone | -6% |
| (2,2) + (3,1) | -2.5% |
| diagonal alone (earlier replica) | -15.5% |
| diagonal + (2,2) + (3,1) | **-33%** |
| variance alone | -37% (noise term at full N 14%) |
| off-diagonal covariance alone | -21% |
| variance + off-diagonal covariance | -37% |
| kappa_3 readouts + diagonal | -56% (earlier replica -63%) |
| kappa_3 readouts + diagonal + (2,2) + (3,1) | **-91%** (raw 1.8e-9) |
| all of the above + variance + off-diagonal covariance | -86% (noisier) |

On network 0 the full oracle set (`oracle_mc2_off0.txt`, raw 2.2525e-8) repeats every entry of the network-1 table with
larger effects: (2,2) slice alone **+24.7%** (worse); (3,1) slice alone -22%; (2,2) + (3,1) -6%; diagonal + (2,2) + (3,1)
**-41%**; kappa_3 readouts + diagonal -72%; with the (2,2) and (3,1) slices also true **-96%** (raw 9.2e-10); variance alone
-52%, covariance alone -34%, both -52%; everything -96%. An independent Monte Carlo replica (`oracle_mc3_off0.txt`) gives
-78% and -98% for the two combined oracles.

- The coherence prediction holds in its sharpest form: one true slice in an otherwise incoherent state makes the
  output worse; the whole sector true is worth twice the diagonal alone; with the kappa_3 readouts true, the sector
  takes the error from -56% to -91%.
- The variance and covariance are accurate to 0.1% and 0.5-0.8%, and true variance alone is still worth 37%: the
  readout's precision demand on the variance is extreme. But they add nothing once the kappa_3 readouts and the sector
  are true (-86% against -91%, within the extrapolation noise), so their errors are downstream of those statistics:
  D21, the (2,2) and the (3,1) slices all enter the next covariance through the pair program.
- With every retained pre-activation statistic true at every layer the residual is 1.8e-9 to 3e-9, of which 6.4e-10
  is the reference's own noise. Law closure of the pair program and the mean map, given a correct retained state, is
  therefore at most about 6-11% of the present error on this network (one network; the extrapolation carries a few
  percent of noise). This bounds the Wick truncation and the mean-map
  truncation. It does not bound the (2,1,1) and (1,1,1,1) law closure, whose effect enters the z slices that the
  oracle replaced; that one has to be separated as in section 6.

## 6. The discriminating computations, in order (proposed, not run)

1. **Representation against law closure, offline on an existing dump.** The one-step run with every oracle on
   (`ostp1`, network 1) stored, per layer, the y state K, d, var_y computed from true z inputs, the chain's own z slices
   before replacement, and the Monte Carlo truth. Evaluating T_z = W#[D(d) + S22(K)] exactly in all three slices from that
   dump splits the chain's one-step defect into representation (chain minus T_z) and law closure (T_z plus a declared
   Lambda_y minus truth), slice by slice, and gives THEORY_3's three diagonal terms their sizes. No new chain run.
2. **Ward tests of the present map and of T_z**, layerwise then terminal (section 4). They need no Monte Carlo and
   decide whether the regenerated sector's defect is a homogeneity defect.
3. **Validate the adjoint prediction** against the existing D3 and g4 oracles.
4. **Only then the terminal factorial** of section 4, with lam frozen, on paired networks, against the sector ceiling.

## Code

- `code/check_k4_contract.py`: the code's diagonal against Q - R_W(E) + 2 lam s_off^2; the three-way split; the slices
  of the transported retained tensor against an explicit 4-tensor; wk4m as the (2,2) slice of J(2I, N). All at
  rounding (4e-15).
- `code/check_radial.py`: the radial tangent on an arbitrary finite law; the first-order cancellation of the
  Gram-Charlier mean map along it.
