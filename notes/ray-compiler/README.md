# The ray compiler and the suffix-risk kernel: what the output can see, checked, and mapped onto the chain

This note takes two supplied developments, the context-typed response compiler and the finite-width He suffix-risk
kernel with its noncentral ray reduction. It checks their claims in exact arithmetic or low-dimensional quadrature
(`code/verify_compiler.py`, `outputs/verify_compiler.txt`) and reads them against the production chain (`est_v29.py`).
The experiments it calls for run on Modal (`infra/modal/mjob.py`).

## 1. What is proved and checked

**(1) The terminal parity law.** relu(t) = t/2 + |t|/2. Hence, for any law of the pre-activation y,

    E relu(y) = E[y]/2 + E|y|/2,

and E|y| depends only on the even part of the law under y -> -y. So the terminal mean reads the law of y only through:
- its mean;
- its even raw moments E y^(2k).

The odd raw moments E y^3, E y^5, ... beyond the mean are invisible. The third cumulant is not invisible: E y^4
contains 4 mu kappa3. This is why the Edgeworth coefficient of kappa3, E relu''' / 6, is proportional to phi'(mu/s) and
vanishes at mu = 0.

On the sphere S^(n-1) the same fact is the supplied spectrum of the single-unit readout. The Funk-Hecke multipliers
mu_l = E[relu(t) P_l(t)] satisfy:
- mu_1 = 1/(2n);
- mu_l = 0 for every odd l >= 3;
- for the even ones, exactly, mu_(2k+2) / mu_(2k) = -(2k - 1)/(n + 2k + 1).

For example mu_4/mu_2 = -1/1027 at n = 1024. All of these are checked in exact rational arithmetic for n = 3, 8, 64 and
1024, and k = 1..5.

One caveat on reading the multiplier. The energy a unit puts in degree 4 is (mu_4^2 N_4)/(mu_2^2 N_2) = 1/12 of its
degree-2 energy at every n. That is the Hermite ratio of |x|. The 1/(n+3) is the harmonic normalisation; it does not
suppress fourth-order directional structure.

**(2) The ray reduction and the two-scalar compiler.** Write a perturbation of a Gaussian reference N(mu, Sigma) as a
density source dp, and read it with a degree-d homogeneous F along rays z = r u. Then

    E_dp F = integral over the sphere of F(u) q(u),   q(u) = integral_0^inf r^(n-1+d) dp(r u) dr.

So a source is invisible to every degree-d readout exactly when q = 0. Along a ray the Gaussian factor is
r^(n-1+d) exp(-A r^2/2 + B r), with A = u^T Sigma^-1 u and B = u^T Sigma^-1 mu. Its radial moments obey the
integration-by-parts recurrence

    m_(k+2) = beta m_(k+1) + (n + d + k) m_k,   beta = B / sqrt(A)   (A normalised to 1).

The null sources are exactly the radial divergences, D_d b = (n + d) b + r db/dr - A r^2 b + B r b. Every polynomial
source therefore reduces, modulo the null space, to c0(u) + c1(u) r: two scalars per direction. The recurrence and the
two-scalar reduction are checked to 1e-14 for n = 5, 40, 200. m_1/m_0 is within 0.04 of sqrt(n + 1/2) + beta/2 there.

**(3) Degree-p null tuples and the context-typed coefficient.** For every positively homogeneous F of degree p,

    (dSigma, dk3, dk4) = (-(p - 2) G/12, mu.G/4, Sigma.G)   is null   (normalised symmetric products),

which is THEORY.md's identity extended to every degree.

- **Check.** Ray integrals in n = 3 at a correlated, non-centred reference give q = 0 to 1e-15 for p = 1, 2, 3. The
  controls are of order one.
- **Gain mode.** (t (Sigma + mu mu^T), 6 t mu.Sigma, 12 t Sigma.Sigma) is null for p = 1 (1e-15) and not for p = 2.
- **Other orders.** The same Gaussian-integration-by-parts and Euler argument gives the k = 3 family
  (dmu, dSigma, dk3) = (-((p - 1)/6) v, (1/3) mu.v, Sigma.v) and the k = 2 radial scaling.
- **The polynomial null space.** Within kappa <= 4 the null space is spanned by these trace-type families,
  Sigma.(lower order) with its mean companions. A kappa4 source is equivalent to kappa <= 3 sources exactly when it is
  a trace core Sigma.G; its Sigma-traceless part is genuine fourth-order content.
- **Context typing.** A core inside a context of derivative order m reads at (p - m - 2)/12. In the kappa3 x kappa4
  term of the Edgeworth mean (m = 3, p = 1) the coefficient is -1/3, not -1/12. This is checked exactly, by
  Gauss-Hermite on the chain's own coefficient formulas. With the outer -1/12 it is off by 62%.

**(4) The minimum-score gauge.** Take the Hermite (L^2(phi)) norm of a source, |S2|^2/2 + |S3|^2/6 + |S4|^2/24, in
whitened coordinates with a = Sigma^(-1/2) mu. Along its null orbit (S2 - G/12, S3 - a.G/4, S4 - I.G) the norm is
stationary exactly when

    c G + 4 tr(G) I + 2 (G a a^T + a a^T G) = 24 D,   c = 4n + 18 + 2|a|^2,   D = S2 + a.S3 + tr12 S4,

which is the supplied equation; the residual is 3e-15. This also identifies the score and D, which the supplied text
left implicit.

The equation is solved in closed form in O(n^2). Split along a and its complement: the operator acts by
c + 4|a|^2, c + 2|a|^2 and c on the parallel-parallel, mixed and perpendicular blocks, and the trace is solved
self-consistently. This agrees with least squares to 7e-15.

**(5) The finite-width transition.** Q_n f(c) = E[a f(c')], with a = (2/n) sqrt(S_A S_B).
- At c = 1, a = 2 S_A/n and E a = 1 exactly.
- For c < 1, E a < 1 strictly (Cauchy-Schwarz), so Q_n is sub-Markov.
- Monte Carlo at n = 16, 64, 256 confirms both: E[a c'] approaches rho(c) as 1/n.

## 2. The suffix-risk kernel as an a priori weight (`code/suffix_kernel.py`, `outputs/suffix_kernel.txt`)

Over random He suffix weights, the expected output energy of a law perturbation at post-activation layer k is the double
integral of k_D, D = 15 - k, against it. Independent pairs of samples sit at the typical cosine c_k:
- c_0 = rho(0) = 0.318, then 0.494, 0.605, ..., 0.923 at layer 14;
- rho(c) = (sqrt(1 - c^2) + (pi - arccos c) c)/pi.

In infinite width, K_D = rho^(oD); finite width moves it by O(D/n), about 1.5% here. Its Taylor coefficients at c_k
weigh the parts of a perturbation by their order in (c - c_k).

| layer k | K' (mean errors) | K''/2 | K'''/6 | K''''/24 |
|---|---|---|---|---|
| 0 | 0.027 | 0.026 | 0.025 | 0.03 |
| 4 | 0.127 | 0.208 | 0.392 | 0.84 |
| 8 | 0.327 | 0.592 | 1.460 | 4.65 |
| 11 | 0.559 | 0.802 | 2.011 | 7.83 |
| 14 | 0.874 | 0.414 | 0.860 | 4.25 |

- **Mean errors.** An error in the mean is damped by the angular contraction prod rho'(c_j) along the suffix: 30 times
  less weight at layer 0 than at layer 14.
- **Higher orders.** These weights peak at layers 10-12. The cone is narrow there, close to the (1 - c)^(3/2) branch
  point of rho, and enough layers remain to amplify. That is where the first-entry ledger found the error entering
  (layers 7-15). This prediction comes from the network's geometry alone, with no fit to the chain.
- **Caveat.** Mapping the order in (c - c_k) onto a cumulant order is a leading-order identification: the
  radius-weighted directional moments of a kappa_j perturbation also have lower-order parts. The table is a prior, not
  a measurement. The measurement is the chain's own suffix response to injected perturbations (section 4).

## 3. What the theory says about the production chain

**(a) The kappa4 sector is a trace core in the wrong metric.**
- **What the chain does.** It regenerates G = diag(dG) + lam C_off and writes its slices with METRIC_C = 2:
  g4row = 2 dG, wk4m = (dG_a + dG_b)/3, wk431 = lam C_ab. That is the trace core (2I).G, exactly.
- **What the identities need.** They hold for the physical trace core Sigma.G. At depth Sigma is not near 2I: the
  pre-activation variance is about 0.16 at layer 15 (2 Var y, with Var y = 0.079 from the Monte Carlo). Note XXXVI
  section 3j's "var ~ 2 under He scaling" is true only at layer 0 and should be read as "var roughly uniform across
  units".
- **Per unit.** The two forms match through the amplitude G^phys_aa = 2 dG_a / var_a. The readout identity
  c4 s^2 = -cv/12 - c3 mu/4 is per unit, so the mean map is unaffected.
- **Across slices.** The (2,2) and (3,1) slices and the K4 -> K3 feed are trace cores of the physical metric only to
  the extent that var_a is constant and the C_ab G_aa terms vanish.
- **The feed.** Its (2,1,1) entries are (2I.G)_aabc = lam C_bc/3. The physical entries are
  (1/6)(var_a G_bc + C_bc G_aa + 4 G_ab C_ac): the shape rho_a dG_a (Phi C Phi)_bc of the newborn's all-distinct kappa3
  is absent, and a scalar lam cannot absorb a per-unit dG_a. Each birth therefore writes kappa3 and kappa4 gain
  amplitudes that the physical covariance does not tie together. By section 1(3) that mismatch is visible, as
  6 c3 mu var dt3 + 12 c4 var^2 dt4 per unit.

**(b) The kappa4 sector can be read as a degree-graded covariance.** By section 1(3), a trace core Sigma.G acts, to
first order, on a readout of degree p as a kappa3 shift -mu.G/4 plus a covariance shift ((p - 2)/12) G. The chain's
readouts have these degrees and weights:
- means: degree 1, weight -1/12;
- post-activation covariance: degree 2, weight 0;
- post-activation kappa3 slices: degree 3, weight +1/12;
- kappa4: degree 4, weight +1/6;
- inside the kappa3 x kappa4 Edgeworth term: -1/3.

So the kappa4 sector need not be carried as slices at all. It can be carried as one symmetric G, the gauge field, read
with degree-dependent weights, while the kappa3 state receives -mu.G/4. This is the representation in which
"quotient before truncation" holds by construction.

**(c) What the terminal law adds.** By itself it adds no free accuracy: the Edgeworth coefficients already encode it.
It does fix what the last layer cannot see, and the risk kernel says how little the early layers' mean errors matter.
Together they bound where work pays.

## 4. Experiments (Modal; 32-network paired evaluation with the dataset's own all-layer truth)

- **Baseline.** Production v29 on nets 0-31 (`base32`): mean raw MSE 2.281e-8, s.d. 0.26e-8, standard error 0.045e-8,
  range 1.72-2.89e-8. Nets 0 and 1 (2.25e-8, 2.13e-8) are typical.
- **Running.** The tests of note XXXVI section 3j:
  - nullalign: the metric of the truth's trace core, and the kappa3 and kappa4 gain amplitudes;
  - cread: the readout at true inputs;
  - kquery: a statistical replicate of section 3i;
  - intervention telescoping: where the error is created, with consistent oracles.

  Monte Carlo mc2 and mc4 were regenerated on Modal in 16 independent chunks each (`k3work/mcchunk.py`,
  `k3work/mcmerge.py`). The shift is fixed by the dataset's own means, so chunks merge exactly.
- **What they returned** (note XXXVI section 3k).
  - Telescoping: the error is created uniformly over layers 7-15, and with true inputs everywhere the output error
    vanishes. The readout at true inputs carries at most a few percent (cread).
  - nullalign: the truth's kappa4 is the gain mode, with t3 = t4 from layer ~11 on. The chain's kappa4 carries 83-85% of
    that amplitude (54-57% in the (3,1) slice), and its kappa3 98%.
  - The Omega audit (V43_NULL_INJ, a diagonal gauge field with every part the chain represents, including the (2,1,1)
    feed term of section 3(a); the M-block subtracts what the legs carry). The chain's output responds to an exact
    null injection with 0.05-0.19 of the covariance part's own response for the gain shape, and 0.09-0.36 for a rough
    shape (`outputs/nullaudit_off01.txt`; the exact value is 0).
  - Where it fails. The residual is created in transport. Its located sources are the birth M-block's rank-4 (2,1)
    residual, as section 3 predicts (the post-activation trace core's (2,1) slice holds the full-rank diag(v) C^y; rank
    64 or the exact rank removes 30-45% of the anomaly), and the lam feed's missing per-unit shape.

## 5. The Gaussian-transport estimator at real magnitudes: exact, 10^4 times quieter, and still far too noisy

Another chat's development proposes a sampled estimator for a cumulant source's effect on the output. One Gaussian
integration by parts turns the score estimator E[h F] into E[U . grad F]: U is the source's Stein field (quadratic for
kappa3, cubic for kappa4), and the expectation is evaluated by one tangent channel through the actual network with the
actual gates. Its use here would be a sampled replacement for the chain's old-source tier, which is about 45% of the
FLOPs. The aged-out sources' first-order effect on the output would be added by sampling instead of carried as
compressed legs; to first order the decomposition is consistent.

**The probe** (`code/transport_probe.py`, `code/transport_probe2.py`; network 0, Monte Carlo truth for the reference and
the sources, 8192 samples with antithetic pairs, `outputs/transport_probe_off0.txt`).
- **The reference at depth is nearly singular.** Dead units make the post-activation covariance rank-deficient, so
  W C^y W^T is too: rank 843-950 of 1024 above 1e-6 lambda_max. The reference is restricted to its support (P the
  pseudo-inverse, U projected onto the range), which is exact because the true kappa3 lies in that range.
- **A source diagonal in unit coordinates is catastrophic.** Take T_aaa = kappa3_a, which is not the physical tensor.
  Its dual coordinates v = P (z - mu) blow up along small-variance directions: the per-sample sd is about 25 per
  neuron, against a signal of order 1e-4. Both estimators are pure noise there.
- **The physical old source.** This is the diagonal class of a source born at the post-activation of layer b0 = b - age,
  carried to z_b by the first-order gated transport: T = sum_i s_i l_i^(x3), with l_i the columns of
  W_b diag(Phi) ... diag(Phi_b0) and s = kappa3(y_b0). Its legs are aligned with Sigma (median l^T P l / |l|^2 of 1.3-3.9,
  against 1/var of 3.8-5.0).

| source (net 0) | per-sample sd per neuron (transport) | score / transport variance | signal rms per neuron | samples for 10% |
|---|---|---|---|---|
| cut 12, age 4 | 1.5e-2 | 2.0e4 | ~2e-4 | ~2e7 |
| cut 9, age 4 | 1.3e-2 | 9.8e3 | ~1.7e-4 | ~5e7 |
| cut 12, age 8 | 3.3e-3 | 6.9e3 | ~5e-5 | ~5e6 |
| whitened-diagonal control (benign) | 7.8e-2 | 2.4e4 | ~1.4e-3 | ~2e6 |

- **The expected-gate control variate.** The same tangent through the fixed gates Phi has an exactly zero mean, since
  E[U] = 0. It removes only a further 1.5-2.3x, at an optimal coefficient of 1.
- **Why.** The signal is a third-derivative, gate-boundary effect, E[T : grad^3 F]/6. The variance is the per-sample
  fluctuation of the actual gates, U . (grad F - E grad F), which no fixed linear control removes.
- **Verdict.** The integration by parts is real (10^4 over the score estimator, consistent with the other chat's 650x
  on its own source). But at the physical source's magnitude a sample is still noise of about 100 times the per-neuron
  signal. 10% accuracy needs 10^6-10^7 samples, about 1e13 FLOPs, against the 1e11 the old tier costs. A sampled old tier
  is not viable at Phase 2 precision. Any sampled component has to target quantities with a far larger signal-to-noise
  ratio per sample than a single source's first-order effect on the output.

## 6. Where the chain's kappa3/kappa4 error lives, and the effective dimension of the state

**The error is visible, not gain** (`code/gainsplit.py`, `outputs/gainsplit_off01.txt`). Each layer's kappa3 diagonal
is split along the gain template 6 mu var, and the kappa4 diagonal along 12 var^2. The chain's error is readout-weighted
by the mean's Edgeworth coefficients c3 and c4, which is what the means read.
- **The truth is mostly gain-shaped.** In norm, the truth's non-gain remainder falls with depth for kappa3 (0.53 at
  layer 2 to 0.26 at layer 15) and grows for kappa4 (0.12 to 0.27). About three quarters of the slices at depth is
  gain-shaped (93% in energy).
- **The chain's error is not.** In readout-weighted norm, the gain component of the chain's error is 1e-5 to 2e-4,
  against 0.7-1.3e-3 for the remainder (the cross-half noise floor is 1-2e-4). The error is 95-99% in the visible part
  on both networks.
- **Consequence.** A gain-quotient representation, carrying the complete gain mode as one exact scalar per layer and
  only the remainder in the legs, cannot repair accuracy. It is a cost question only: does the remainder compress
  better than the whole?

**The state's effective dimension contracts with a power-law spectrum** (`code/effdim.py`,
`outputs/effdim_off0.txt`). These are spectra of the pre-activation covariance C^z_l, network 0.

| layer | participation ratio (tr C)^2/tr C^2 | directions for 90 / 99 / 99.9% of tr C | decay p, lambda_k ~ k^-p (ranks 10-300) |
|---|---|---|---|
| 0 | 511 | 521 / 789 / 915 | 0.33 |
| 3 | 240 | 331 / 654 / 839 | 0.85 |
| 7 | 129 | 229 / 548 / 743 | 1.23 |
| 11 | 81 | 179 / 493 / 700 | 1.41 |
| 15 | 51 | 140 / 428 / 632 | 1.54 |

- **The decay.** The participation ratio halves roughly every three layers early on and every five later. The decay
  exponent rises steadily.
- **No sharp cutoff.** This is the measurable form of the fractal/Ahlfors intuition: there is no low-rank cutoff, but
  strong and increasing concentration. The post-activation spectra are the same within a few percent. At depth the
  means dominate the second moment: mean^2/var of y grows from 0.47 to 12.
- **The fixed-rank tiers sit across this profile.** The old tier's rank 320 lies at the measured output knee
  (288/320/352: +7.5/+2.2/+0.1%), and that knee sits near k99 at depth. The young tier is dense (rank 1024) at every
  layer, while 99% of the state's variance needs 428-548 directions from layer 7 on.
- **What the profile prices.** It suggests where a fixed-rank design spends FLOPs that the state does not use. It
  cannot give the price itself: the legs need Hadamard (Khatri-Rao) products for their readouts, so the leg rank and the
  covariance rank are not the same object.

## 7. Operator theory of the legs: free spectra, the readout law, and the rank a tier needs

The chain truncates the legs of aged sources to a shared basis. Its ladder came from frontier scans (notes XVII, XXIV):
dense to age 4, rank 320 for ages 5-8, and a nested rank-192 sub-basis for ages 9 and older. This section derives the
ranks from the operator theory of the legs, states predictions before measuring, and then measures.

**7.1 The legs are free products** (`code/legspec.py`, `outputs/legspec_off0.txt`).
- **The object.** A source born at the post-activation of layer b reaches the pre-activation of layer l through the
  gated product J(b->l) = W_l D_(l-1) ... W_(b+1) D_b, with D = diag Phi(alpha) the first-order gate. Age a = l - b.
- **The free model.** Treat the He weights and the gates as freely independent. W^T W tends to 2 MP(1), whose
  normalised second moment rho = m2/m1^2 is 2. The gates enter as D^2, with rho_Phi = E Phi^4 / (E Phi^2)^2.
- **The law.** For free positive elements the normalised second moments add, rho(x [x] y) = rho(x) + rho(y) - 1 (the
  first two coefficients of the S-transform). So the leg's participation ratio is

      PR(b->l) = (tr J^T J)^2 / tr (J^T J)^2 = n / (1 + a + sum_(k=b)^(l-1) (rho_Phi(k) - 1)).

  Without gates this is the Fuss-Catalan value n/(a + 1).
- **The check.** On network 0 the law holds to within 2% up to age 8, and within 11% at age 12:

| (birth, age) | (2,1) | (2,4) | (5,2) | (8,4) | (8,6) | (11,4) | (2,12) |
|---|---|---|---|---|---|---|---|
| exact PR | 367.1 | 119.0 | 210.1 | 114.5 | 80.8 | 115.9 | 37.7 |
| free-probability prediction | 367.3 | 119.8 | 210.9 | 115.8 | 80.7 | 115.1 | 42.2 |

- **What it says.** From layer 3 on the gate spread is rho_Phi = 1.9-2.1, so each layer of age costs about two
  Fuss-Catalan factors: PR ~ n/(1 + 2a).
- **The whole spectrum is a free product** (`legspec.py --tails`, `outputs/legtails_off{0,1}.txt`). A surrogate with
  the same gates and fresh He Gaussian weights reproduces the tails eps(a, r), at every rank 64-512 on both networks:
  - to within 1.7% up to age 7;
  - to within 4% at ages 8-10;
  - to within 7% at ages 11-12 (median 1-4%).

  So the rank a leg needs is set by its age and the gate statistics, not by the particular weights.
- **Where the real legs depart.** Their participation ratio falls below the free value at old ages (37.7 against
  42.2 at age 12). The cause is a few heavier top directions, not the tail.
- **Observability.** Weight each row by its observability o_l, the diagonal Gramian. Its terms are the layer's own
  reads, K'_l (c3^2 + rho^2 E Phi^2) for D3 and D21, plus transmission Phi_l^2 o (W_(l+1)^2)^T o_(l+1). The
  sqrt(o)-weighted spectra are only a little more concentrated: PR 280-288 against 336-367 at age 1, and equal at the
  ages the tiers hold.
- **Transmission dominates.** The reads carry 3-14% of o_l up to layer 9, then 18%, 22%, 30%, 44% and 85% at layers
  10-14. The production join metric diag(Phi^2 + 0.1) (`V32_JOIN_POST=3`) is the one-step, uniform-o form of this
  Gramian.

**7.2 What a tier truncates.**
- **The tails.** The tail eps(a, r) is the fraction of |J|_F^2 outside the top r directions. Network 0, mean over
  birth layers 2, 5, 8, 11 (network 1 agrees within 10%):

| age | r = 128 | 192 | 256 | 320 | 384 | 512 |
|---|---|---|---|---|---|---|
| 1 | 0.49 | 0.34 | 0.22 | 0.14 | 0.086 | 0.025 |
| 2 | 0.31 | 0.17 | 0.088 | 0.043 | 0.020 | 3.1e-3 |
| 3 | 0.20 | 0.087 | 0.036 | 0.014 | 4.9e-3 | 4.0e-4 |
| 4 | 0.13 | 0.046 | 0.016 | 4.8e-3 | 1.3e-3 | 5.4e-5 |
| 5 | 0.087 | 0.027 | 7.6e-3 | **1.9e-3** | 4.1e-4 | 1.0e-5 |
| 6 | 0.059 | 0.015 | 3.4e-3 | 6.8e-4 | 1.1e-4 | |
| 8 | 0.028 | 5.0e-3 | 8.1e-4 | 1.1e-4 | 1.1e-5 | |
| 9 | 0.020 | **3.0e-3** | 3.9e-4 | 4.1e-5 | | |
| 10 | 0.014 | 1.7e-3 | 1.8e-4 | 1.5e-5 | | |
| 12 | 7.2e-3 | 6.5e-4 | 5.0e-5 | 2.9e-6 | | |

- **One tail for both tiers.** The production ladder leaves about the same tail, 2-3e-3, at the youngest member of
  each tier (age 5 at 320, age 9 at 192): the frontier scans found equal precision at both boundaries. The spectra
  have a hard lower edge, so the tails fall ever faster: the local exponent p = -dlog eps/dlog r is 6.7-7.8 at
  (age 5, rank 320) and 6.5 at (age 9, rank 192).
- **Why the basis can be shared.** It is the forward Lyapunov subspace: J(b->l) = J(b'->l) J(b->b') for b < b', so
  the range of every older leg lies inside a younger one's. A tier's rank is set by its youngest member.

**7.3 Two laws, and what they predict.**
- **Readout law.** The readouts are Hadamard contractions of leg entries (D3_i = sum_j A_ij^2 P_ij w2_j and the like).
  A tail eps of incoherent entry error is therefore a readout error of relative energy about eps. To first order the
  MSE grows linearly in the tail, Delta MSE = sum_a kappa_a eps_a(r_a), with kappa_a the output weight of age a's
  readouts.
- **The old per-source truncation table.** Frontier tests, exact SVD per source, all ages >= a_0 at rank r, network 0.
  It fits the law with kappa falling by about two per age between ages 4 and 8 (section 7.4 gives the full fit).
- **Error-share law.** Suppose a tier's error falls as r^-p and its rank-proportional cost is a share s of the bill.
  Stationarity of S = MSE x C in the tier's rank gives

      p Delta MSE_t / MSE = s_t.

  - **How flat the optimum is.** Near it, delta log S ~ (1/2)(p + 1) s (delta log r)^2. A +-10% rank change costs under
    1%, which is why the 304/320/336 scan of note XXIV was flat to within its noise.
- **Per-age form.** Stationarity in each age's own rank gives eps_a = gamma N_a r_a MSE / (C kappa_a p_a). Here N_a is
  the number of layers age a spends in the tier and gamma the cost per unit rank per member-layer. Since kappa_a falls
  with age, the optimal tail grows with it: at age 9 it is 3-7 times the age-5 tail, where the ladder holds them
  equal.


**7.4 Predictions, then measurements** (32 networks, paired against production; `outputs/confinement_32nets.txt`).

Stated before the runs:
- **P1, no confinement** (`V21_NO_CONFINE=1`, every source dense). Raw 4-7% lower. That is about 6% if the readout
  error falls with rank at the spectral tails' exponent (p ~ 5-7), and 3.6% if at the steeper r^-9 that the other
  reading of note XXIV's rank knee would imply.
- **P2, one tier** (`V24_AGE_OLD2=0`, every old age at 320). Raw 0.35-0.9% lower if kappa halves per age (0.9% if it
  falls by 0.6 per age), at 3-6% more FLOPs.

| variant | raw, paired | better on | F / B (cold) |
|---|---|---|---|
| production (320, nested 192) | 2.2806e-8 | | 0.2131 |
| P1: no confinement | -4.96% (per net -4.87 +- 0.47) | 32/32 | 0.2918 (+37%) |
| P2: one tier at 320 | -0.36% (per net -0.36 +- 0.14) | 22/32 | 0.2198 (+3.1%) |

Both fall inside the predicted ranges.

**What the two measurements fix.**
- **The exponent.** The shared basis costs X = 5.2% of the dense chain's MSE. Against note XXIV's knee (rank 320 to
  288: +5.3% raw; 320 to 352: -2.1%) this puts the error's local rank exponent at p ~ 5.6-6.7. That is the
  spectrum's own exponent at the age-5 boundary (6.7-7.8), so the readout law holds with the spectral p, not with a
  steeper effective one.
- **The confinement error is an independent, additive error** (`code/confdelta.py`, `outputs/confdelta_off0-7.txt`,
  networks 0-7).
  - What independence predicts. Let delta be the output difference, dense minus production. If e_prod = e_dense -
    delta with delta orthogonal to e_dense, then |delta|^2 = -Delta MSE, the cross term is -2|delta|^2 and
    cos(e, delta) = -sqrt(|delta|^2/MSE).
  - What is measured: |delta|^2 = 5.3% of the MSE, Delta MSE = -4.6%, cross term -10.0% (-10.7% predicted) and
    cos -0.215 (-0.23 predicted). The nested tier alone: |delta|^2 = 0.57%, Delta MSE = -0.56%.
  - Where it sits by layer. The confinement's share of the per-layer error is 2.5% at layer 5, where the tier starts.
    It is a stationary 5-6% from layer 9 on.
  - What this licenses. The linear, additive form of the readout law is the right one.
- **The calibration.** The readout law with kappa_a = kappa_5 q^(a-5), calibrated on P1 and P2, gives kappa_5 = 22.2
  and q = 0.46: kappa halves every 0.9 ages (`code/ladder_opt.py`, `outputs/ladder_opt.txt`). Age 5 carries 4.0 of
  the 5.2 points, age 6 0.7 and age 9 0.3; every other age is under 0.1.
- **Tier 1 is at its optimum.**
  - The cost side. The rank-proportional cost is fitted from P1's and P2's FLOP differences and the profile of note
    XXIX: forming plus the factor-space D21 hub is 2.3e-3 units of 2 n^3 per rank per member-layer, and the join is
    8.8e-3 per rank per layer.
  - The comparison. That makes s_1 = 0.25 of the bill, against p X_1 ~ 7 x 4.8% = 0.33: within the flatness of the
    optimum.
- **Tier 2 is over-provisioned, by a margin worth nothing.** p X_2 ~ 0.02 against s_2 ~ 0.05: the law wants a lower r2.
  But at that flat optimum the score moves by about 0.1%.

**The ladder the theory selects** (`outputs/ladder_opt.txt`; score = raw x C with C the cold FLOPs plus the residual
term, relative to production):

| ladder | predicted confinement error | predicted score |
|---|---|---|
| production: 5-8 at 320, >= 9 at 192 | 5.22% | 0 |
| one tier: >= 5 at 320 (P2, measured +2.7%) | 4.84% | +2.5% |
| best two-tier: 5-8 at 336, >= 9 at 160 | 4.29% | -0.6% |
| best three-tier: 5-6 at 352, 7-9 at 256, >= 10 at 128 | 3.99% | -0.7% |
| age 4 confined too (AGE_OLD 3), best rank 416 | 3.50% | +0.6% |

- **The production ladder is the theory's optimum to within 0.7%.** No tier structure that the law allows is worth
  more. The dense boundary at age 4/5 is right: confining age 4 costs more accuracy than its FLOPs buy at every rank.
  The per-age refinement that section 7.3 predicts (lower ranks at older ages) is real but is worth 0.1 point beyond
  the best two-tier ladder, which needs no new code.

**The one test** (stated before the run). The best two-tier ladder needs only the existing switches:
`V21_R_OLD=336 V24_R_OLD2=160`. Predicted on the 32 networks:
- raw -0.9%, as the confinement error falls from 5.22% to 4.29%: age 5 at 336 gains 1.3 points, ages 9-11 at 160 lose
  0.6;
- FLOPs +0.3% (0.2138 B);
- score about -0.6%.

A result outside raw -0.5% to -1.3% would mean the per-age split (kappa_5, q) is wrong, not just imprecise.

**The test, interim (19 of 32 networks; the Modal workspace was paused mid-run).**

| | predicted | measured |
|---|---|---|
| raw, paired against production | -0.9% (range -0.5% to -1.3%) | -0.12% (per net -0.07 +- 0.36), better on 11/19 |
| F / B | 0.2138 | 0.2137 |

- **The cost model is exact; the error model is not.** By the criterion stated in advance, the per-age split
  (kappa_5, q) is wrong, at about 2 standard errors.
- **Which half is wrong cannot be told from this run.** The prediction is the sum of two parts: -1.3 points for age 5
  at 336 and +0.6 for ages 9-11 at 160.
  - The likelier failure is the second. Rank 160 lies below the measured point (192). The nested sub-basis is a
    two-pass range finder in factor space, not the optimal SVD of the tails, and note XVII's tier-2 scans saw a cliff
    at 128.
  - Splitting the two halves takes two runs: rank 336 alone, and nested rank 160 alone.
- **What stands.** The production ladder stays. The measured optimum (P1, P2 and the error-share law for tier 1) is
  unaffected. What failed is only the extrapolation of the per-age split below the measured tier-2 rank.
