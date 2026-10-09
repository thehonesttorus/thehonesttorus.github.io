# Measured structure of the adopted system's error (100 official networks unless stated)

Sources: chain dumps (`workbench/k3work/chaindump.py`, job cdump1), Monte Carlo cache at 2^20 inputs per network
(`mccache.py`, job mcc1), analyses `ana_loc.py` / `ana_loc2.py` (jobs analoc1/analoc2). The ground truth's own noise
is ~1e-10 at the last layer (negligible); MC-cache statistics have relative noise sqrt(2/2^20) = 1.4e-3 on variances.

**Per-layer MSE of the adopted system** grows about linearly with depth, ~1e-9 per layer: 6.5e-10 (layer 1),
1.46e-9 (3), 3.6e-9 (6), 6.9e-9 (9), 1.09e-8 (12), 1.55e-8 (15). The neuron-averaged signed error is below 1e-5 at
every layer: the error is per-neuron scatter, not a common shift.

**Linear (first-chaos) share of the pre-activation variance**, |E[X h_i]|^2 / Var(h_i), median over neurons:
0.735 (layer 1), 0.51 (3), 0.36 (6), 0.29 (9), 0.24 (12), 0.215 (15). Deep pre-activations are 78% higher chaos.

**Pre-activation covariance spectrum (MC)**: participation ratio / n = 0.39 (layer 1), 0.19 (4), 0.12 (7), 0.087
(10), 0.063 (13), 0.052 (15) [~53 effective modes]; top-16 modes carry 8% (1) -> 42% (15) of the variance, top-128
carry 47% -> 89%. The fluctuations collapse onto a few collective modes with depth.

**Where the final-layer error was injected.** Decompose e_l = diag(Phi(alpha_l)) W_l e_(l-1) + inj_l (mean-gate
transport of the previous layer's mean error plus a newly injected part) and push every inj_l to the last layer with
the chain's mean gates (exact telescoping; 8 networks). Share of the final squared error by injection layer:
layers 0-5: 3%; 6-8: 10%; 9-11: 23%; 12: 11%; 13: 15%; 14: 18%; 15: 21%. Injections from different layers are
uncorrelated (shares of |C_l|^2 sum to 0.996). Mean-gate transport is contractive (an injection at layer 3 loses
~9x in squared norm by layer 15), and injection sizes grow with depth (8e-10 at layer 3, 3.2e-9 at layer 15).
=> The last 4 layers' injected error is ~2/3 of the final MSE.

**Share of each layer's own squared error that is newly injected** (rest transported): 0.32 (1), 0.55 (3), 0.44
(5), 0.37 (7), 0.32 (9), 0.27 (11), 0.24 (13), 0.21 (15).

**Error by gate margin** |alpha| = |W out_(l-1)| / sigma (share of the layer's squared mean error):
layer 1: [0,0.5) 0.53, [0.5,1) 0.30, [1,1.5) 0.13; layer 15: [0,0.5) 0.14, [0.5,1) 0.13, [1,1.5) 0.12, [1.5,2) 0.11,
[2,3) 0.17, [3,inf) 0.34. At depth a third of the error sits in near-linear (large-margin) neurons, i.e. transported.

**Collective-mode content of the error.** The chain's pre-activation mean error d_l = W_l e_(l-1), projected on the top
k eigenvectors of the MC covariance Cov(h_l): top-16 / 64 / 128 carry 8-15% / 27-43% / 46-63% (random vector:
1.6% / 6.2% / 12.5%), the enrichment rising with depth. The error is 4-9x concentrated in the collective modes but
37-54% of it lies outside the top 128.

**Earlier oracle attribution (note XXXI, birth-address/README.md section 3, networks 0-1):** replacing the chain's
third-cumulant readouts (D3, D21) by truth at every layer removes ~40% of the final MSE; the kappa4 diagonal 15-30%;
both ~63-72%. The chain's variance has relative error 4e-4 (layer 1) rising to 1.5e-3 (layer 15).
