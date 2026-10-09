# XLIII. One skeleton under two papers, and where the chain's variance error lives

Status: 9 October 2026. Sections 1-5 record what was read, derived and measured; section 6 is the pre-registration of
the next experiments, committed before their data.

## 1. The request

The user pointed to two papers and suggested that their theory and methods are noncommutative geometry in disguise,
that they combine with this programme, and that the combination could go further:
- arXiv 2510.04929, *Efficient Quantum Hermite Transform* (Jain, Iyer, Somma, Bao, Jordan);
- arXiv 2606.09728, *Quantum Cut Sparsifiers* (Basu, Brakensiek, Kothari, Putterman).

The work was done without workflow fan-out:
- both papers read in full;
- the theory developed directly;
- the measurements run on cached Monte Carlo data (16M inputs, two independent halves; networks 0-1);
- the chain's own internal slices dumped once per network (`code/` here and `workbench/k3work/var_*.py`, `chaindump2.py`).

## 2. What the papers do, and the noncommutative geometry under them

### 2a. Quantum Hermite transform

**The factorization is an identity in SL(2, R).** The paper fast-forwards the oscillator H = (x^2 + p^2)/2 through

    e^(-iHt) = e^(-i tan(t/2) p^2/2) e^(-i sin(t) x^2/2) e^(-i tan(t/2) p^2/2).

It holds because x^2, p^2 and {x, p} close into sp(2, R) = sl(2, R) (their Claim 69). It is the three-shear
decomposition of a rotation, so it holds in every representation of the metaplectic group. In imaginary time it is the
Mehler kernel of the Ornstein-Uhlenbeck semigroup, written as heat flow, then Gaussian tilt, then heat flow.

**The discretization is a rational noncommutative torus, controlled on a spectral truncation.**
- x is a diagonal grid and p is its Fourier conjugate: the clock and shift matrices generate M_M(C), the
  noncommutative torus at theta = 1/M.
- The canonical relation [x, p] = i fails globally. The paper shows it holds to exp(-gamma N) on the span of the first
  N ~ M/log M eigenvectors of H (Theorem 6).
- In noncommutative-geometry terms, this is a spectral truncation in the sense of Connes and van Suijlekom: the
  compressed operator system approximates the Moyal plane's.
- The proof splits nested commutators into few and many nestings:
  - a polynomial of degree t moves the Hermite level by at most t (finite propagation in the level grading);
  - the long tail is controlled by 1/t!.

**Hermite sampling measures the number operator.** It samples |f_hat(v)|^2, the spectral measure of N in the state f,
which is the Wiener chaos decomposition of f. Gaussian Goldreich-Levin finds the heavy chaos coefficients from prefix
weights.

### 2b. Quantum cut sparsifiers

**The Hamiltonian is a group-algebra element.** L_G = sum_e w_e (1 - SWAP_e) is the image of
grad_G = sum_e w_e (id - (e)) in R[S_n] under the permutation representation on (C^2)^(x n).
- The level-k Kikuchi graph is the k-particle sector.
- The harmonic subspaces E_k (their Lemma 3.5) are the Schur-Weyl isotypic components: the Specht module (n-k, k)
  tensored with spin (n-2k)/2.
- The up and down operators satisfy [D, U] = n - 2k on level k (their Lemma 3.3). This is su(2), the compact twin of
  the oscillator sl(2, R) in section 2a.

**Every inequality is proved in the group algebra, hence in all representations at once.** Three are used:
- Alon-Kozma: an expander dominates a scaled complete graph, grad_G >= lambda grad_(K_n).
- The octopus inequality of Caputo, Liggett and Richthammer.
- Chen's moving-particle lemma: the effective resistance R_G(u, v) times L_G dominates L_uv at every level.

Positivity in C*(S_n) is what makes "all levels simultaneously" free.

**The sparsification argument.**
- In each isotypic block the complete-graph operator is the scalar Casimir k(n + 1 - k). Matrix Chernoff then works
  block by block.
- Leverage decays like n/i at level i for expanders.
- It does not decay for important edges (their Remark 8.5: a bridge of a tree keeps its weight at every level).
  Expander decomposition separates the two kinds.

### 2c. The common skeleton, with stage 10

Both papers have the same shape:
- a small algebra acts on an exponentially large space with a level grading:
  - in the Hermite transform, sl(2, R) and the Hermite levels;
  - in the sparsifier paper, R[S_n] with its su(2) commutant and the Kikuchi levels;
- identities and inequalities are proved in the small algebra, hence hold in every representation, hence at every
  level at once;
- finite approximations are controlled on spectral truncations (Hermite transform) or by a universal domination
  (sparsifiers).

The two dual pairs are Howe's (O(n), sl(2, R)) on L^2(R^n) and Schur-Weyl's (S_n, SU(2)) on qubits.

Stage 10 makes the same move for Iwahori-Hecke algebras:
- the gap holds in every *-representation, by Bruhat-monotone coupling and the q-Tits cone lemma;
- its Temperley-Lieb/ASEP corollary is the q-deformation of the path case of the Aldous property, which the sparsifier
  paper uses at q = 1 through the octopus inequality.

## 3. The chain in the same skeleton

### 3a. Dictionary

| chain | skeleton |
|---|---|
| neurons | modes |
| cumulant tower (note XL) | bosonic Fock space; level = token number |
| linear layer W | second quantization: W^(x k) on level k (a GL(n) action that preserves levels) |
| ReLU layer | site-local ladder operators: the Hermite vertex at each neuron |

The Gaussian input is Howe's setting:
- the Gaussian closure lives in the oscillator representation;
- the chain's levels are the token levels of section 2b.

### 3b. Theorem 1 (Mehler-ladder form of the two-site covariance)

Let (h_a, h_b) have means mu, standard deviations s, correlation rho and small joint cumulants k_pq (p copies of a,
q copies of b). Let psi_a = relu(alpha_a + .) in L^2(gamma) and c_m(alpha) = <psi_a, He_m>. Then

    E_G[relu^(p)(h_a) relu^(q)(h_b)] = s_a^(1-p) s_b^(1-q) sum_m c_(m+p)(alpha_a) c_(m+q)(alpha_b) rho^m / m!
                                    = s_a^(1-p) s_b^(1-q) < d^p psi_a , rho^N d^q psi_b >_gamma,

where rho^N = e^(-tN) at e^(-t) = rho is the Mehler (Ornstein-Uhlenbeck) operator and d is the annihilation operator of
the Hermite basis. The first-order Edgeworth expansion of the covariance is therefore

    Cov(relu(h_a), relu(h_b)) = K_G + sum_(p+q>=3) k_pq / (p! q!) * [ladder-inserted Mehler elements].

Here K_G, the Gaussian closure, is a matrix element of the imaginary-time metaplectic rotation of section 2a.

*Proof.* Stein's identity E[F'(z) He_m(z)] = E[F(z) He_(m+1)(z)], applied p and q times, together with Mehler's
formula E[He_j(z_a) He_k(z_b)] = delta_jk k! rho^k.

*Check* (`code/edge_test.py`). On synthetic non-Gaussian laws the corrections reduce the covariance error from
3e-2 to the Monte Carlo floor (4e-4 to 1e-3).

### 3c. Theorem 2 (two-site sufficiency)

Var(h_(t,i)) = sum_ab W_ia W_ib Cov(relu(h_a), relu(h_b)), so a layer reads the previous layer's law only through its
two-site marginals.

In the slice harmonic analysis of section 2b, these are the components of each level that depend on at most two
sites. Section 4c measures how few levels are needed: at width 1024, levels 3 and 4 suffice to 1e-4.

### 3d. Theorem 3 (the gate covariance: what is omitted, and its collective part)

**The expansion.** Expand relu(h_c) in Hermite polynomials of the standardized pre-activation. For centered X,
Leonov-Shiryaev gives

    k(X_c^2, X_d^2, X_e) = 4 C_cd k3(c,d,e) + 2 C_ce k3(c,d,d) + 2 C_de k3(c,c,d) + k5(c,c,d,d,e),

verified by Monte Carlo (`code/ls_check.py`: 0.0994 against 0.0981, standard error 0.001).

**What each term is.** Put w2 = E relu'' = phi(alpha)/s. The (2,2,1) class of the next layer's slice k(z_a, z_a, z_b)
contains:
- the gate covariance sum W_ac W_ad W_be w2_c w2_d Phi_e C_cd k3(c,d,e), which joins two transported tokens by a
  covariance line;
- the GC1/GC2 terms: two-site slices moved along one covariance line. They are carried and folded (note XXXIX).

Configurations in which two of the sites c, d, e coincide are two-site objects, which the chain's pair programs
compute. What it omits is the all-distinct gate covariance: the Khatri-Rao wall of notes XV, XXXV and XL.

**Its collective part is CP.** Write C = sum_(j<=k) lambda_j U_j U_j^T + C_bulk. The collective part of the omitted
class is

    sum_j lambda_j k(zeta_j, zeta_j, y),   zeta_j = W diag(w2 o U_j) (h - mu),   y = W diag(Phi) (h - mu),

minus its coincident-site part. These are third cumulants of linear transforms of h. On a CP source
sum_r u_r x v_r x w_r they are sum_(j,r) lambda_j (W(w2 o U_j o u_r)) x (W(w2 o U_j o v_r)) x (W(Phi o w_r)): CP rank
k per term, at a cost of O(k R n^2) per layer.

**The sparsifier analogy.** This is expander decomposition again:
- the collective modes are the important edges, kept exactly;
- the bulk off-diagonal of C, with entries O(n^(-1/2)), is the expander-like part, whose transported weight is the
  open question;
- the coincident-site part is the mean field, and the chain already carries it.

**What it costs.** At the production CP rank (R ~ 1e3) and k = 16 the collective part costs about 16 units per layer.
It pays only if it carries a large share of the error at few layers. Section 6 measures the ceiling first.

### 3e. A closure of the fourth cumulant by the third (heuristic)

For a field linear-plus-quadratic in a Gaussian latent (second Wiener chaos, the metaplectic case of section 2a),
every cumulant is a walk sum of the same matrices:
- open walks L^T H ... H L;
- closed walks tr(H^k), the closed walks of note XLII's frame N3.

So k4 is determined by k3 and the latent geometry.

**The single-factor closure (derived).** Take h_a = L_a.xi + (1/2) xi^T H_a xi with small H and xi ~ N(0, I). At
leading order:
- kappa3(a) = 3 L_a^T H_a L_a;
- k(a,a,b) = 2 L_a^T H_a L_b + L_a^T H_b L_a;
- k(a,a,a,b) = 6 L_a^T H_a H_b L_a + 6 L_a^T H_a^2 L_b (counting the twelve open walks on {a,a,a,b}).

When neuron a's quadratic part is aligned with its own linear direction (H_a L_a = h_a L_a), eliminating h_a and the
H_b term gives the parameter-free closure

    k(a,a,a,b) = 2 kappa3(a) k(a,a,b) / v_a - (2/3) kappa3(a)^2 C_ab / v_a^2.

**How it fits the truth** (`code/k31_sc.py`, `code/k31_mix.py`):
- **From true inputs it has the right absolute scale with nothing fitted.** The entrywise amplitude is 1.02-1.10, the
  correlation with the true (3,1) slice rises from 0.43 (layer 5) to 0.86 (layers 13-14), and the variance-metric
  amplitude is 1.0-1.1.
- **It is complementary to the chain's scale-mixture (gain) closure.** That closure, 3 tau v_a C_ab, correlates at
  0.70-0.76.
- **The combination removes more than either.** In the variance metric, the chain's closure alone removes 32-53% of
  the true (3,1) term's energy. Adding the second-chaos closure, computed entirely from the chain's own kappa3, D21,
  v and C, at unit coefficient removes 46-64%. The best per-layer mixes sit near (1.1, 0.75) on both networks.
- **A weaker local form has a stable coefficient.** Its exchange part, kappa3(a) k(a,a,b) / v_a alone, explains 9-35%
  of the chain's (3,1) error with a coefficient of about 1.5 at every layer (`outputs/k31_fix_net0.txt`).

## 4. Measurements

### 4a. Oracle ceilings (production counterterms on)

Noise-extrapolated a = 2 MSE(full) - mean(MSE(h0), MSE(h1)), against base 1.478e-8 (network 0) and 1.574e-8
(network 1); Monte Carlo statistics replace the chain's at the listed layers:

| oracle | layers | net 0 | net 1 |
|---|---|---|---|
| VAR (pre-activation variances) | 12-15 | -32% | -37% |
| VAR | all | -56% | -48% (halves differ 2x: noisy) |
| VAR + COFF + D3 + D21 | 12-15 | -52% | -56% |
| D21 | 12-15 | -9% | -13% |
| D3 | 12-15 | -2.5% | -13% |
| COFF | 12-15 | +4.5% | -28% (noisy) |
| G4 (kappa4 diagonal) | 12-15 | +16% | +34% |
| WK4M + K31 | 12-15 | +17% | +15% |

The all-layer rows of section 6a, run after the pre-registration (`outputs/oracle_all_layers.txt`):

| oracle, every layer | wiring | net 0 | net 1 |
|---|---|---|---|
| D21 | default | -26% | -13.5% |
| D21 | consistent (V41_ORC_CONSIST=1) | -20% | -5% |
| K31 | (unaffected by the wiring) | -22% | -14% |
| D21 + K31 | default | -54% | -44% |
| D21 + K31 | consistent | -50% | -33% |
| D21 at 12-15 | consistent | -1% | -10% |

About the wiring:
- The default oracle also feeds the oracle D21 into the birth step's gated subtraction, so newborn sources absorb the
  difference between the oracle and the chain's own value.
- The consistent wiring keeps that subtraction on the chain's own values. This is the experiment as registered: the
  legs keep carrying their own content.
- The halves sit 10-20% above the full-sample values in the D21 runs, so the noise extrapolation carries most of
  those numbers.

The late error lives in the second-order statistics of the late layers, above all the per-neuron variances.

### 4b. The variance error is propagated injection (`outputs/var_decomp_nets01.txt`)

**The split.** With KG(state) the pure Gaussian closure, every per-neuron variance error splits exactly as

    dv = p + (x - g):
- p (propagated): KG of the chain's state minus KG of the true state, read through the next layer;
- x (the chain's correction): what the chain adds beyond KG;
- g (the true defect): what the true law adds beyond KG of the true covariance.

Relative to the true variance, network 0:

| layer | dv mean | dv rms | true defect g (noise-free rms) | injected x - g (noise-free rms) | propagated share of dv energy |
|---|---|---|---|---|---|
| 4 | -1.7e-4 | 4.3e-4 | 4.3e-3 | 1.7e-4 | 0.42 |
| 9 | -3.1e-4 | 9.1e-4 | 7.7e-3 | 4.7e-4 | 0.66 |
| 12 | -5.3e-4 | 1.39e-3 | 9.8e-3 | 6.4e-4 | 0.73 |
| 15 | -1.03e-3 | 2.03e-3 | 1.07e-2 | 7.5e-4 | 0.84 |

Network 1: layer 15 has dv -1.05e-3 / 2.07e-3, injected 7.7e-4, propagated share 0.81.

**What it says.**
- The pure Gaussian closure would be wrong by about 1% per neuron. The chain removes 99.5% of that energy.
- What it leaves is injected afresh at every layer, growing from 1.7e-4 to 7.5e-4.
- Late in the network 80% of the error is inherited.

**The propagation runs through the chain's own variance errors.** Splitting p by its source:
- the chain's variance errors carry the mean, -9.5e-4 at layer 15;
- the chain's means carry -7.5e-6;
- the chain's off-diagonal covariance carries +1.0e-5.

**The uniform part is not a neutral mode.** For V = v f(alpha) the diagonal channel multiplies a uniform relative
variance error by

    lambda = sum (W^TW)_aa v_a dV_a/dv / sum (W^TW)_aa V_a,   v dV/dv / V = 1 - alpha f'(alpha) / (2 f(alpha)).

It measures 0.93-0.97, which damps. The observed 1.0-1.28 growth per layer of the mean error therefore comes from where
the errors sit: on the dead side, where Var(relu) is most sensitive to v.

### 4c. The true law is perturbative to fourth order (`outputs/var_ladder_net{0,1}.txt`)

E3 and E4 are the first-order Edgeworth corrections built from the TRUE two-site cumulants, read through the next
layer's rows:
- E3 from k(a,a,b) and kappa3;
- E4 from k(a,a,b,b), k(a,a,a,b) and kappa4.

| net 0, layer | g | g - E3 | g - E3 - E4 | x - g (injected) | share of x - g carried by x - E3 - E4 |
|---|---|---|---|---|---|
| 4 | 4.3e-3 | 5.4e-4 | 1.6e-5 | 1.7e-4 | 0.99 |
| 9 | 7.7e-3 | 6.7e-4 | 4.9e-5 | 4.7e-4 | 0.99 |
| 15 | 1.07e-2 | 7.0e-4 | 1.1e-4 | 7.5e-4 | 0.98 |

Network 1 is the same to two digits (layer 15: 1.05e-2, 6.5e-4, 1.1e-4, 7.7e-4, 0.98).

**The meaning.** At width 1024 the exact two-site law is first-order Edgeworth to order 4, to 1e-4 at every depth.
No non-perturbative (mixture) closure is needed. The chain's injection is entirely its deviation from that Edgeworth,
that is, errors in its own third- and fourth-order two-site cumulants.

**The size of the residual.** It is about 10x the naive bulk counting (order 5 ~ n^(-3/2)), as expected from the
collective modes.

### 4d. Which slice (`outputs/var_attrib_net{0,1}.txt`)

**The decomposition is exact.** The chain's own Edgeworth terms, built from its dumped slices at its own state,
reproduce its correction x to 1e-4. The dumps are what its Wick stage reads, and the stage is first-order Edgeworth.

**The per-layer injection by term** (share of the noise-free injection energy):

| layer | D21 = k(a,a,b) | K31 = k(a,a,a,b) | K22 | kappa3 | kappa4 | left |
|---|---|---|---|---|---|---|
| net 0, 6 | 0.70 | 0.20 | 0.09 | 0.02 | 0.00 | 4.6e-5 of 2.9e-4 |
| net 0, 12 | 0.80 | 0.13 | 0.05 | 0.02 | 0.00 | 1.1e-4 of 6.4e-4 |
| net 0, 15 | 0.78 | 0.20 | 0.03 | 0.01 | -0.01 | 1.2e-4 of 7.5e-4 |
| net 1, 12 | 0.74 | 0.23 | 0.04 | 0.01 | -0.01 | 1.1e-4 of 6.0e-4 |
| net 1, 15 | 0.76 | 0.21 | 0.03 | 0.02 | -0.01 | 1.2e-4 of 7.7e-4 |

**The two slices that matter.**
- **The (3,1) slice is carried at half amplitude.** The chain's term is about half the truth's at every layer from 6
  to 15 on both networks (2.4e-4 against 5.1e-4 at layer 15 on network 0; 1.1e-4 against 2.1e-4 at layer 6 on
  network 1).
- **D21 is right to 93% in the variance metric, but entrywise its correlation with the truth falls from 0.96
  (layer 5) to 0.87 (layers 13-14)**, and the chain's |D21| is 3.0e-3 against 3.6e-3. The entrywise difference
  (`outputs/d21_probe_net0.txt`):
  - correlates +0.54 with the true D21 and +0.57 with the regression form k3(a) C_ab / v_a;
  - correlates -0.43 with C;
  - holds 20% of its energy in column means: a coupling of neuron b to the layer's total energy.

### 4e. Cheap completions do not reach it (`outputs/d21_fix_net{0,1}.txt`)

The probe fits n^2 completions of the two slices per layer, in sample, against the noise-free variance-metric error:
- a D21 amplitude;
- the regression form;
- C;
- the chain's K31;
- k3 D21 / v;
- 3 v C.

**The ceiling is low.**
- The amplitude alone explains 10-23%: the truth wants D21 2-4% larger. Production's output-fitted V47
  calibration moves it 1-2% smaller, so the counterterms are compensating something else.
- All six terms together explain 14-29%.

The remaining three quarters is three-site structure: the transported kappa3 tensor, which no two-site completion
sees.

## 5. Synthesis

**Where the variance error lives.**
- The late per-neuron variance error, worth a third to a half of the final MSE, is accumulated injection.
- Each layer's injection is the chain's error in two two-site slices:
  - kappa(a,a,b), 70-80%, mostly three-site transport error;
  - kappa(a,a,a,b), 15-25%, carried at half amplitude.
- The exact law needs nothing beyond fourth-order Edgeworth.

**One-step shares do not predict single-slice oracle values (section 6a).**
- Fixing kappa(a,a,b) alone at every layer is worth -5% to -26%. Fixing kappa(a,a,a,b) alone is worth -14% to -22%.
- Together they are worth -33% to -54%, the size of the whole variance ceiling.
- The two slices' one-step errors are only mildly anti-correlated in the variance metric: -0.07 to -0.20
  (`outputs/dk_corr_nets01.txt`). The super-additivity therefore comes from the full chain: counterterms fitted around
  the chain's own slices, propagation, and the Monte Carlo noise of all-layer substitution.
- For the system, the oracle values are the operative facts. The (3,1) slice is the cheaper lever: it is carried at
  half amplitude, and the second-chaos closure of section 3e recovers part of it at n^2 cost.

**Why the heat defect could not see it.** This is consistent with note XLII's Theorem A5 (first-order localization is
blind to the trace-state, closed-walk content). The injection is a quadratic form in the next layer's rows,
w^T (Delta_chain - Delta_true) w, built from two-site slices of the third and fourth cumulants. The identification is a
conjecture, not a theorem.

**What the papers contribute.**
- **The Hermite transform paper** supplies the Mehler-ladder form (Theorem 1) and the metaplectic closure heuristic
  (3e).
- **The sparsifier paper** supplies the design principle for the wall: keep the collective modes exact and treat the
  bulk by its mean field. The chain already does the mean-field part through its coincident-site programs.

**The open theoretical questions.** Both are the sparsifier paper's leverage-decay question in this setting.
- Does the transported weight of the bulk part of the all-distinct gate covariance decay with age, as leverage decays
  with level for expanders?
- Does the network's tower transport have an Aldous property: does every level contract at least as fast as the mean
  level?

## 6. Pre-registration (committed before the data)

### 6a. Ceiling of the two-slice lever (oracle, all layers)

**Setup.** Production estimator with its counterterms, networks 0-1, Monte Carlo files full / h0 / h1, noise-extrapolated
as in 4a. The oracle replaces the chain's readout slices at every layer; its legs keep carrying its own content.

**Predictions:**
- `V37_ORACLE=D21`: -30% to -45% on both networks.
- `V37_ORACLE=K31`: -3% to -12%. This is uncertain, because the k31 counterterms (up to +-0.17) were fitted around
  the chain's own slice.
- `V37_ORACLE=D21,K31`: -35% to -55%.

**Decision.**
- If D21-all is at most -20%, the two-slice lever is not worth an n^3-per-layer component. The note records that and
  the gate-covariance route is closed.
- If D21-all is at least -30%, run experiment G (6b).

**Result** (table in 4a; `outputs/oracle_all_layers.txt`).
- D21 alone: -20% / -5% consistent (-26% / -13.5% with the default wiring). The prediction (-30% to -45%) fails on
  both networks and both wirings.
- K31 alone: -22% / -14%. The prediction (-3% to -12%) fails on the high side.
- D21 + K31: -50% / -33% consistent (-54% / -44% default). The prediction (-35% to -55%) holds on network 0 and on the
  default wiring, and is 2 points short on network 1 consistent.
- Decision: D21 alone stays at or below -20% in the registered (consistent) wiring, so the gate-covariance route is
  closed for now and experiment G is not run. The lever moves to the (3,1) slice (6c).

### 6b. Experiment G: the gate-covariance budget (run only if 6a passes)

**Setup.** Monte Carlo on networks 0-1 at source layers 9 and 13. The fields are y = W diag(Phi)(h - mu),
q = W diag(w2/2)((h - mu)^2 - v) and zeta_j (j <= 16). Third cumulants of type (a, a, b) are formed, and the carried
GC1/GC2 and coincident-site parts are subtracted in closed form from the true two-site slices. Everything is read in the
variance metric two layers on.

**Predictions:**
- (G1) The one-step all-distinct gate covariance has noise-free rms at least 0.3 of the chain's D21-term error, and
  correlates at least 0.4 with it.
- (G2) The top-16 collective part carries at least 50% of the gate covariance's energy (top 4: at least 30%).

**Decision.**
- G1 and G2 both hold: price the collective component on the young sources at the last layers, then a cold screen
  and a scored run.
- G1 holds but G2 fails: the wall stands. Record it; the next question is importance-sampling the bulk.
- G1 fails: the D21 error comes from transport of older content, and the next measurement is by source age.

### 6c. The second-chaos (3,1) closure in the production chain (cold screen)

**The change.** `V58_K31SC=beta` adds, at every layer where the chain has a (3,1) slice,

    wk431[a, c] += beta (2 kappa3(c) D21[c, a] / v_c - (2/3) kappa3(c)^2 C[c, a] / v_c^2).

It uses the chain's own D3, D21, var and C_off after the V47 calibration (wk431[a, c] = kappa(z_a, z_c, z_c, z_c)).
The cost is a few n^2 elementwise operations per layer. Production is unchanged when the switch is off.

**Harness.** Networks 0-15, each variant paired against production rerun in the same batch.

**Predictions:**
- `V58_K31SC=1` (parameter-free): raw -1% to -6% (mean paired change), better on at least 10/16, FLOPs +0.1% or less.
  The basis: the oracle says a perfect (3,1) slice is worth -14% to -22%, and in the variance metric the closure
  removes a further 20-30% of the slice's remaining error energy. The k31 counterterms were fitted around the chain's
  own slice, which can eat part of the gain.
- `V58_K31SC=0.5`: -0.5% to -4%, better on at least 9/16.

**Decision.** If beta = 1 meets -1% or better on at least 10/16, run the scored test: 100 networks, adjusted MSE, then
a held-out refit of the k31 calibration. Otherwise record the result and stop.

**Result** (`outputs/k31sc_cold_16nets/`; base mean raw 1.577e-8; network 0 also ran once locally as a crash check
before the pre-registration was committed, giving +2.3% at beta = 1, the same as the batch):

| variant | mean paired raw change | median | better on | FLOPs |
|---|---|---|---|---|
| `V58_K31SC=1` | +1.28% +- 0.77 (se) | +2.5% | 6/16 | +0.02% |
| `V58_K31SC=0.5` | -0.76% +- 0.37 | -0.3% | 10/16 | +0.02% |

**The parameter-free closure fails its prediction.** The half-strength version meets its own (-0.5% to -4%,
>= 9/16), but the gain is small. By the registered rule there is no scored run.

**The reason.** Built from the chain's own D21, the closure correlates with the true (3,1) slice at 0.53-0.56; built
from the true D21 it correlates at 0.86 (`outputs/k31_sc_net0.txt`). The closure therefore imports the chain's D21
error into the (3,1) slice.

**What the closure says about the two slices.** In the second-chaos regime the true (3,1) slice is largely made of
the true (2,1) slice:

    k(a,a,a,b) ~ 2 kappa3(a) k(a,a,b) / v_a - (2/3) kappa3(a)^2 C_ab / v_a^2.

This structural tie is the natural reading of the oracles' super-additivity in section 6a, although the one-step
error correlation is weak. The two slices are one object: a K31 closure is only as good as the D21 it reads. The
(3,1) lever is therefore not independent of the D21 wall after all.
