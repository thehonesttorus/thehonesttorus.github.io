# XLII. Localization: the chain's error as its drift along stochastic localization

Status: in progress (9 October 2026). Sections 1-3 record what is established; section 4 is the pre-registration of
the decisive experiment, committed before its data.

## 1. The question and how it was attacked

The user's direction: put localization processes (Eldan's stochastic localization, the linear-tilt schemes of
Chen-Eldan and Anari-Koehler-Vuong, pinning in simplicial complexes, Kikuchi/token graphs) inside Connes'
noncommutative integration (transverse measures on groupoid algebras), find the natural home, and read off what it
says about a genuinely new Phase 2 system. Six frames were developed in parallel, each with an adversarial referee
(F1 heat-flow defect and iterated expectations; F2 activation-pattern pinning and the transverse measure of linear
regions; F3 collective-mode disintegration; F4 modular tilts and KMS/Petz covariances; F5 convex-order sandwiches;
F6 localization couplings for hybrid estimators). Four Elicit Research Agent sessions surveyed the literature
(deterministic uses of localization, Connes' integration and inference, deterministic estimation of network outputs,
pinning/Kikuchi). Full developments and referee reports: `frames/`; Elicit deliverables: `elicit/`.

## 2. Measured structure of the adopted system's error (100 official networks)

From new caches of the chain's per-layer state (`code/chaindump.py`) and of Monte Carlo statistics at 2^20 inputs per
network (`code/mccache.py`), analysed by `code/ana_loc.py`, `code/ana_loc2.py` (outputs in `outputs/`):

- **The error grows about linearly with depth** (~1e-9 of MSE per layer, 1.55e-8 at layer 15) and is per-neuron
  scatter: the neuron-averaged signed error is below 1e-5 at every layer.
- **It is injected late.** Writing e_l = diag(Phi(alpha_l)) W_l e_(l-1) + inj_l and transporting every injection to
  the last layer with the chain's mean gates (exact telescoping), the share of the final squared error injected at
  layers 12-15 is 0.11 / 0.15 / 0.18 / 0.21, at layers 0-8 together 0.13. Injections at different layers are
  uncorrelated, mean-gate transport is contractive (an injection at layer 3 loses ~9x in squared norm by layer 15),
  and the injected size grows with depth.
- **Deep pre-activations are mostly non-linear in the input**: the first-chaos share of their variance falls from
  0.74 (layer 1) to 0.22 (layer 15).
- **The fluctuations collapse onto collective modes**: the participation ratio of the pre-activation covariance falls
  from 0.39 n to 0.052 n (about 53 modes); the top 128 modes carry 89% of the variance at layer 15.
- **The error is enriched in those modes but not confined to them**: the pre-activation mean error's top-16 / 64 / 128
  shares are 8-15% / 27-43% / 46-63% (random: 1.6 / 6.2 / 12.5%).
- **No sign rule and no gain shape** survive in the production residual (F5's prediction for a chain that carries the
  gain mode): frac(truth < chain) = 0.47-0.55 at every layer from 4 on; |corr| with (1 + alpha^2) phi(alpha) sigma and
  (alpha^2 - 1) phi(alpha) sigma is below 0.03 at layers 7-15.

## 3. What the frames established (refereed)

- **The natural home is commutative, and it is a heat equation.** Any Gaussian-input estimator E(m, Sigma) (a chain run
  on the input law N(m, Sigma)) that is exact on point masses has error equal to minus the integral of its heat
  defect D = d_Sigma E - (1/2) Hess_m E along any Gaussian localization path (Dynkin; F1 Theorem 1, F5 independently;
  Elicit's survey derived the same identity). With homogeneity the leaf space is the ray quotient: e(y) = E(y, I),
  delta = e - y.grad e - Lap e, and truth - e(0) = -(1/2) int (1 + tau)^(-3/2) E_{y ~ N(0, tau I)} delta(y) dtau.
  The exact functional is the Doob martingale; every chain's error is its integrated iterated-estimation gap.
  The 'transverse measure' is this Green (occupation) measure; the groupoids of the problem (rays, rotations, gate
  patterns) are type I or abelian (F2, F6), so no genuinely noncommutative structure carries leverage here.
- **The base-point defect is a truth-free per-neuron error detector at small width.** At n = 48-256 the defect
  delta(0) = e(0) - Lap_m e(0) explains 59-98% (Gaussian closure) and 74-88% (dense kappa_3 chain) of the per-neuron
  error, with err ~ -delta/(2p); a held-out merge e + a delta cuts MSE by 76-81% (F1), and the parameter-free midpoint
  (e + Lap e)/2 cuts a Gaussian closure's MSE 8.7x (F5; reproduced by its referee). Corrected class weights (F5
  referee): first-order localization sees the kappa_4 path at weight 1, the kappa_3 path at 2, the gain's kappa_3
  half at 1, and the gain's kappa_4 half and closed walks at 0.
- **Negatives, each with a quantified wall.** Pinning gates (F2): the omitted pair-gate class is about a third of the
  MSE but no carrier costs less than ~100 units. Convex order (F5): partition-function sandwiches are blind to the
  mean; moment sandwiches are 50-1500x too wide. Sampling hybrids (F6): for any positive-weight rule with nodes
  independent of W, N x MSE >= 0.54 sigma^2, so any unbiased add-on improves the chain by at most ~0.24% at 0.1 B.

## 4. Pre-registration: does the heat defect explain the production chain's error at n = 1024?

**Quantity.** For networks 0-3, the adopted system (`estimator_final_v56.py`, counterterms on) run in float64 with the
saturation drop off (`V60_F64=1 V33_SAT=nan`; the only remaining discontinuities are continuous clips and maxima),
and the directional heat defect d_r = [E(0, I + h v v^T) - E(0, I - h v v^T)]/(2h) - [E(hv, I) + E(-hv, I) -
2E(0, I)]/(2h^2) over K = 128 random unit directions per network (`code/heatdef.py`), h fixed by the smoothness pilot
below. delta_hat = 2 (n/K) sum_r d_r per neuron and layer; the two halves of the directions give the noise correction
(`code/heatdef_ana.py`).

**Smoothness pilot (gate, run first).** Network 0, two directions, h = 0.5 / 0.25 / 0.125, with and without
counterterms. At n = 256 the defect is h-stable and about 1% of the covariance part at layers 1-12 but jumps to ~1e-2
and becomes h-dependent from layer 13 on (a non-smooth path switching on at depth, calibrated for n = 1024). The main
run uses the largest h at which the per-layer defect agrees across h within 20% at every layer; if no h passes at the
last layers, the experiment is reported for the smooth layers only.

**Predictions (registered before the data).**
- X* (noise-free share of the per-neuron error explained by delta at the final layer): F1's prediction 0.6-0.85;
  its referee's prior 0.45 (0.3-0.6). Registered: 0.45, interval [0.25, 0.70].
- Slope a in err ~ a delta: -0.15 to -0.25 (p = 2-3.5).
- Single-direction noise rho (std of n d_r over |mean|): 3-6.
- By layer: X* is not lower at layers 12-15 than at 6-11.

**Decision rule.** X* >= 0.5 at the final layer, with the held-out slope (fitted on networks 0-1, applied to 2-3)
stable within 20%: derive the defect's per-layer local form for the last 3-4 layers in the source representation and
build a runtime delta only if its analytic form captures >= 70% of this measured delta at <= +70% of the bill; then a
cold 16-network screen and a scored run. X* < 0.3: no runtime component; the identity stands as the explanation of
the chain's error. In between: record the by-layer profile and stop.
