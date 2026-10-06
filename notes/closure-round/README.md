# The closure round

Working note XVIII. Note XVII measured that the representation round of the v29 chain is converged (2.7% of its
raw error) and that everything else, including the whole 1.7x to the leaders, is closure: the model by which the
chain replaces the law of the pre-activations by a few carried cumulants. This note opens that round the way Prodan
opens a finite-volume round: name the approximation, write what it drops, size the dropped part, and test derived
replacements end to end on the official truth. It also answers the second reading of the Smale-space paper that
was pasted into the session (section 1), because that reading proposes a representation, and the representation
question is closed.

## 1. The symbolic reader is the moment chain

The pasted reading asks for "a weight-accessible, two-sided continuation coding with controlled local overlap,
correct conditional weights, and a small operator representation for filling or integrating gaps", and shows on
the golden-mean shift that a block transfer matrix integrates non-commuting readouts over all fillings of a hole
without enumerating them. Translated to the network, every piece of that is already in the chain or is excluded by
a number we have:

- **The symbols are gate patterns and the fillings are the 2^n patterns of a layer.** A transfer matrix over symbols
  is exponential; the only small compositional representation of E[prod W_l D_l] under the Gaussian input law is the
  moment representation, in which the conditional expectation over patterns given the carried moments is the gated
  propagator W diag(Phi) at first order and the kappa_3 sources at the next. The "small operator representation for
  integrating gaps" is the chain. Whether "the readout closes in an accessible small compositional representation" is
  the closure question of this note.
- **The Parry law is the annealed law.** Uniform weighting of legal fillings is the maximal-entropy measure; the
  true measure on patterns is the Gaussian pushforward, whose density against the product law is the conditional
  gate covariance of note XIV (the Rota-Baxter defect). Section 7 of the pasted reading asks whether that tilt is
  bounded above and below; it is not (the cell measure is multifractal, note XIV section 5), so Ahlfors regularity
  of the symbolic model does not transfer to the state.
- **The stratification bound kills cell quadrature in our dimension.** The pasted reading's own estimate is
  MSE = O(N^(-1-2 alpha/s)) for a conditional sample per cell. The Gaussian input law on R^1024 has s = 1024 and
  the readout is Lipschitz, alpha = 1, so the exponent is 1 + 1/512: Monte Carlo with a stratification that costs
  more than Monte Carlo. The dimension argument of that section is correct and it is fatal to the idea here.
- **The heat operator exists and the chain already truncates it.** The spectral calculus the pasted reading
  reaches through the Deaconu-Renault groupoid of the shift is, for the Gaussian input space, the Ornstein-Uhlenbeck
  semigroup: the pair kernel E[sigma(h_a) sigma(h_c)] = sum_k c_k(a) c_k(c) rho_ac^k / k! is the Hermite spectral
  expansion of the heat operator with the correlation as time. The chain's pair programs keep k <= 2 (section 2 of
  the `est_v29.py` term table: `('c_off','ones2')` and `('c_off','c_off')`). So the chain is a spectral truncation of
  the heat operator, and the truncation order is a closure parameter we can move.
- **Two-sided compatibility, retained history, local charts.** The sources fixed at birth and transported by
  pullback (note XVI) are the retained history; the shared basis and its joins are the local charts; both rounds are
  measured converged. The paper's bracket has no cheap analogue for the mean because the output boundary is the
  whole final vector, not a conditioning event.

What survives of the pasted reading is its standard of proof: that a representation must transport geometry,
state and readout together. The chain's representation does; its closure is where the state is approximated.

## 2. The closure, named

At every layer the chain holds, for the pre-activation h = W y: the mean mu, the variance var, the off-diagonal
covariance C_off, the third-cumulant diagonal D3 and (2,1) slice D21 (from the transported sources), and a
fourth-cumulant sector: the diagonal g4row, the (2,2) slice wk4m and the (3,1) slice wk431. The post-activation
moments are then computed by Gaussian integration by parts around the Gaussian with (mu, var, C_off), to first
order in the cumulant slices (the `TERM_SPECS` table), with the Gaussian part expanded in Mehler order <= 2 in
C_off. The fourth-cumulant sector is not transported content: it is regenerated each layer as
g4row = 2 dG, wk4m_ac = (dG_a + dG_c)/3, wk431 = lambda_l C_off, with dG a transported diagonal and lambda_l a
per-layer table fitted on public dumps (`LAM`, times an adaptive rule `BETA`). E11 found the fitted sector to be
0.48-0.84 of 3 g sigma^4 and called it a counterterm.

**The fitted sector is the scale mixture's, at 80%.** From the layer dump of network 0 (`code/k4sector.py`), with
g read from the chain's own D3 by projection on 1.5 mu var (g = 0.013, 0.019, 0.022 at layers 5, 10, 14), the
derived slices of the scale mixture z = G y, Var G^2 = g (E10), are in closed form

    k4_a = 3 g var_a^2,   K22_ac = g var_a var_c (+ 2 g C_ac^2),   K31_ac = 3 g var_a C_ac

(checked against `pairvar.sm_slices` to 0.1-5%), and the chain's fitted sector has these shapes exactly:

| layer | g4row / k4_sm (norm, cosine) | wk4m / K22_sm (norm, cosine) | derived K31 vs C_off, vs d(var) C_off (cosine) |
|---|---|---|---|
| 5 | 0.768, 0.997 | 0.773, 0.998 | 0.995, 1.000 |
| 10 | 0.808, 0.996 | 0.815, 0.997 | 0.990, 1.000 |
| 14 | 0.815, 0.992 | 0.830, 0.994 | 0.979, 1.000 |

So the chain's lambda table and its adaptive rule are a fit of one number per layer, 0.77-0.83, in front of the
derived scale-mixture fourth-cumulant sector, and its (3,1) ansatz lambda C_off differs from the derived
3 g d(var) C_off only by the per-neuron variance factor.

**The truth of the fourth-cumulant sector** (`code/k4truth.py`, `outputs/k4truth_off{0,1}.txt`): the Monte Carlo
accumulators of E10 give the pre-activation's kappa_4 slices directly. Projected on the scale-mixture shapes, the
diagonal and the (2,2) slice carry one amplitude g_4 per layer (R^2 0.90-0.94 from layer 4 on), and the (3,1) slice
projected on 3 d(var) C_off carries the same amplitude within noise (its R^2 is low because fourth-order
off-diagonal moments are noisy at T = 4e5, but the implied g tracks the diagonal's at every layer on both networks).
So the sector's shape is the scale mixture's. Its amplitude is not the kappa_3 amplitude: g_4 / g_3 (g_3 from E10's
kappa_3 channel) rises from 0.72-0.75 at layer 5 through 0.84-0.85 at layer 10 to 0.96-1.05 at layers 14-15, on both
networks. The chain's fitted diagonal, 0.77-0.83 of 3 g_3 var^2, is therefore the truth at layers 5-10 and about 15%
low at layer 14: the "counterterm" of E11 is the real amplitude of a mixture whose fourth cumulant grows into its
third with depth. The chain's (3,1) slice, lambda C_off with lambda = 0.0069-0.0103, is 0.55-0.80 of the truth's
3 g_4 d(var) C_off at depth, the one slice it has under-fitted.

| layer (net 0 / net 1) | g_4 from k4 diag | g_4 from K22 | g_4 from K31 | g_3 (E10) | g_4 / g_3 |
|---|---|---|---|---|---|
| 5 | 0.0105 / 0.0100 | 0.0106 / 0.0100 | 0.0089 / 0.0084 | 0.0140 / 0.0138 | 0.75 / 0.72 |
| 10 | 0.0166 / 0.0159 | 0.0171 / 0.0163 | 0.0150 / 0.0139 | 0.0197 / 0.0188 | 0.84 / 0.85 |
| 14 | 0.0215 / 0.0201 | 0.0220 / 0.0210 | 0.0226 / 0.0210 | 0.0223 / 0.0201 | 0.96 / 1.00 |

For a scale mixture z = G y with G^2 = 1 + eps, the exact fourth cumulant is 3 sigma^4 E eps^2 + 1.5 sigma^2 mu^2
E eps^3 + ..., and the third is 1.5 sigma^2 mu (E eps^2 - E eps^3 / 4) + mu^3 E eps^3 / 8 + ...: the two
amplitudes differ by the skewness of the activity fluctuation, and their measured ratio is the statement that the
fluctuation is skewed and that the skewness falls with depth. That is the next derivation; the chain's lambda
table is its fitted shadow.

**The Mehler truncation drops 1e-4 of the covariance.** The off-diagonal correlations of the pre-activations have
rms 0.07, 0.10, 0.12 at layers 5, 10, 14 (99.9th percentile 0.23-0.42). The dropped Mehler terms of the (1,1) and
(2,1) programs, relative to the first-order term in Frobenius norm, are k = 3: 0.8-1.5e-4, k = 4: 0.2-1.4e-4;
summed coherently over entries (the part a transport along the mean direction keeps) k = 4 is 2-4e-4 of k = 1
because even powers of the correlation add. That is the size of the chain's final relative error, so the
truncation is at the precision edge and is free to raise (two Hadamard powers per layer, 2 n^2 flops).

**The Edgeworth tail of the scale mixture.** ReLU is homogeneous, so for the scale mixture E relu(G y) = E[G] m
exactly, with E[G] = 1 - g/8 at first order: the mean error of the mixture is -(g/8) m = -(g/8) sigma (alpha Phi +
phi). The chain represents it through Edgeworth terms: kappa_3 = 1.5 g mu sigma^2 gives -(g/4) alpha^2 phi sigma
and kappa_4 = 3 g sigma^4 gives (g/8)(alpha^2 - 1) phi sigma, whose sum is -(g/8) sigma phi (1 + alpha^2). The
difference is the tail of the even cumulants, -(g/8) mu (Phi - alpha phi):

| alpha | exact / (g sigma) | kappa_3 | kappa_4 | tail | share kept by kappa_3 + kappa_4 |
|---|---|---|---|---|---|
| 0 | -0.050 | 0 | -0.050 | 0 | 1.00 |
| 0.5 | -0.087 | -0.022 | -0.033 | -0.032 | 0.63 |
| 1.0 | -0.135 | -0.061 | 0 | -0.075 | 0.45 |
| 1.5 | -0.191 | -0.073 | +0.020 | -0.139 | 0.28 |
| 2.0 | -0.251 | -0.054 | +0.020 | -0.217 | 0.13 |

The chain's alpha are wide (median 0.4-0.5, a 30-39% tail above 1, 90th percentile 2-4 at layers 5-14), so if
the mixture's mean error reached the output only through these two terms, half of it would be missing and the
chain's error would be 2e-3 per neuron, not 1.5e-4. It is not missing: the arrow-mixture note measured that the
sources and second-order terms carry 70% of E[G_L] - 1 and the kappa_4 channel 30%, and that conditioning the
chain on the mixture on top double-counts (2.19e-8 to 1.84e-6). The per-layer Edgeworth accounting and the
end-to-end accounting are two descriptions of one object; the tail term is a test of which one the chain's
fitted 0.8 is paying for.

## 3. Tests (network 0, base = warm join + feedback, residual rank 4, Strassen level 6: raw 2.1785e-8)

(results pending: the ladder is running; this section is filled in the next commit)
