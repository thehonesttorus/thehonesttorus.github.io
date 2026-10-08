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

