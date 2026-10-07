# One-sided memory, Gibbs, and the scale as a conserved charge

Working note XIX. The question asked was whether the recent work on Markov partitions, locally Markovian (one-sided)
descriptions and Gibbs measures has subtle links to the chain, beyond the Smale-space reading of note XIV. Ten papers
were read (section 8). The standard is the one of notes XVII and XVIII: a link counts when it changes a number or a
decision; a relabelling is recorded as a relabelling. The development that the reading forced is in section 2: the
open derivation of note XVIII, why the fourth-cumulant amplitude of the pre-activations is half the third-cumulant
amplitude at layer 1 and equal to it at layer 15, is done on the Monte Carlo truth with a two-channel ledger, and it
turns on the one structural fact that the Gibbs literature names and the network has: a conserved mode of
non-summable memory, which here is the scale.

What the reading changed:

- **The chain is a g-measure, not a Gibbs measure, by construction** (section 1). Depth is time, the gate pattern of a
  layer given the carried state is the g-function, and the exact object would also condition on the output, which has
  no cheap conditioning event (note XVIII section 1). The variation sequence of this g-function is summable on every
  mode but one: the scale, whose transport is exactly 1 by positive homogeneity. The uniqueness theory of g-measures
  says bounded-memory closures cannot be exact along such a mode, and the chain's history confirms it: the scale is
  the one mode that had to be given a global variable (GAC's gain, then the coherent third-cumulant sources).
- **The scale is a conserved charge, and its two amplitudes are reconciled** (section 2). For a scale mixture y = G y~
  the map relu(G y) = G relu(y) is exact, so the mixture is an invariant family of the depth dynamics and its transport
  is 1 in every cumulant channel. The fresh non-Gaussianity that ReLU creates at a Gaussian state is not a mixture: its
  fourth-cumulant amplitude is 0.48 of its third at layer 1 on both networks and falls to 0 at depth. The measured
  ratio g_4 / g_3 = 0.50 at layer 1 rising to 1.0 at layer 15 is the relaxation of that skew towards the invariant
  mixture direction, with a derived rate that goes from 0.46 per layer at layer 2 to 0.97 at layer 15 as the gates
  saturate and cumulant orders decouple. The closed two-channel recursion reproduces the measured fourth-cumulant
  amplitude to 1-4% through layer 11 once homogeneity is imposed on the transport rows, and the third-cumulant
  amplitude saturates for a reason the recursion names (section 2.4).
- **Quenched against annealed, with a number** (sections 3 and 6). The finite-width hierarchies (Hanin, Yaida, Celli,
  Kawase-Ota) are the annealed theory: third cumulants vanish by the W -> -W symmetry, and the fourth grows as 5 L / n.
  The chain is the quenched theory of one network under one input law: its content is the third cumulant the annealed
  theory does not have, and its fourth-cumulant growth is 0.023 at layer 15 against the annealed 0.073, because the
  weight fluctuation is frozen in one network. The one place the annealed/quenched distinction costs the chain a
  number is its fourth-cumulant table, fitted once for both networks: the derived quenched replacement moved the two
  networks in opposite directions (-2.5% and +1.1%, note XVIII), and the networks' own g_4 differ by 7% at layer 14.
- **Everything else is a relabelling** (sections 4, 5): the KMS/groupoid characterisation of Gibbs states is, for a
  Gaussian pushforward, the Ornstein-Uhlenbeck modular structure whose truncation the chain already is; the local
  product structure of Gibbs measures is what the pushforward lacks (a common scale across all directions is the
  opposite of a product), and we had measured that; explicit Gibbs constants have no uniform analogue because the only
  contraction rates we can name are depth-dependent and vanish at depth.

## 1. One direction or two

van Enter, Le Ny and Paccaut (2011.14664) separate two notions that are usually conflated. A measure on sequences is
a **g-measure** (locally Markovian, "almost Markov in one direction") when the conditional law of the present symbol
given the whole past is a continuous function of the past; it is **Gibbs** (DLR, "in two directions") when the
conditional law of a finite window given everything outside it is continuous. Neither class contains the other:
long-range Dyson-type Gibbs measures fail the one-sided continuity, and there are g-measures whose backward conditional
law (the present given the future) is discontinuous, so they are not Gibbs. The distinction is exactly the chain's.

**The chain is one-sided.** Depth is time. The symbols are the gate patterns of the layers (note XVIII section 1), and
the chain carries, at layer l, a finite state (mean, variance, off-diagonal covariance, the transported third-cumulant
sources, the regenerated fourth-cumulant sector) from which the law of the next pattern is read: at first order the
gate probabilities Phi(alpha_i), then the cumulant corrections. That reading is a g-function: the present given a
summary of the past, with the summary of bounded size. The exact object of the challenge, E[z_16], is a one-sided
functional of the pattern process, so the one-sided description is the natural estimator. The two-sided object would
condition the pattern of layer l on the output as well; that is the bracket of the Smale-space reading, and note XVIII
section 1 showed it has no cheap conditioning event because the output boundary is the whole final vector. Its exact
form is the n^4 contraction of note XVIII section 5. So the chain's class is fixed: it is a bounded-memory g-measure
approximation of a process that is not, as far as we can tell, Gibbs in any useful sense (section 5).

**The variation sequence has one non-summable mode.** The g-measure literature controls the dependence of the present
on the remote past by the variation sequence var_k(g), how much the conditional law changes when the past is altered
beyond k steps back; summable variation (and, by Johansson and Oberg, square-summable) gives uniqueness of the
g-measure and the exponential mixing that makes bounded-memory truncations converge. The chain's analogue is measured:
the response of the final mean to third-cumulant content born at layer b decays as rho^age with rho = 0.80-0.92 on the
bulk of the channel (note XVII, E10), but with rho = 1 on the mean direction. Section 2 of this note measures the same
thing for the scale mode of the fourth-cumulant channel: the transport of a unit-gain scale mixture has row sums
1.005 and 1.04 at layer 1, and positive homogeneity says exactly 1 for the full cumulants. So the variation sequence
of the network's g-function is geometric on every mode but one, and on that one it is constant: the memory of the
scale is infinite, neither decaying (summable) nor growing (unstable). That is the He-critical network in the
language of the Dyson model: at the border of the uniqueness regime, and it is why no truncation of the carried
state at finite memory can be exact along the scale. The chain's record is the same statement from the other side:
the gain was the one quantity that had to be made a global variable (GAC), and it is now carried by the coherent
sources whose transport is critical.

**What the one-sided reading rules out.** A one-sided finite-memory description is complete when the g-function is
a function of the carried state. It is not: the gate-covariance term of note XVIII section 5 (1.1% of the transported
third cumulant, incoherent, a third of the error) is the dependence of the present pattern on past content that the
carried moments do not summarise, and it costs n^4 per source and layer to carry. In g-measure terms the chain's
truncation error is a non-vanishing tail of the variation sequence in the pattern coordinates, not in the moment
coordinates; the moment coordinates are the ones in which the tail is smallest, which is the representation statement
of note XVIII again.

## 2. The scale as a conserved charge: the two-channel ledger

### 2.1 Homogeneity

Let y be the pre-activation of a layer and suppose it is a scale mixture, y = G y~ with y~ Gaussian with the carried
mean and covariance, G independent of y~, E G^2 = 1 and Var G^2 = g. Positive homogeneity of ReLU gives
relu(G y~) = G relu(y~) exactly, and the next pre-activation is z = W relu(y) = G W relu(y~) = G z~. The scale passes
through the layer unchanged; the layer only adds, inside z~, the fresh non-Gaussianity of relu at a Gaussian state.
So the family of scale mixtures is invariant under the depth dynamics, with G a conserved charge, and if the fresh
non-Gaussianity were itself a scale mixture with a new independent factor G', the product G G' would carry a gain
g + g' at first order: the additive, critical accumulation that Hanin's log-normal product law (2204.01058, eq. 4.7)
gives for the annealed activity norms, and that the recursion kappa_4^(l+1) = 5 K^2 / n + kappa_4^(l) of the same
paper writes for ReLU.

The measured amplitudes of note XVIII contradict the naive version of this: with g_3 the coefficient of
kappa_3(z_k) = 1.5 g mu_k sigma_k^2 and g_4 the coefficient of kappa_4(z_k) = 3 g sigma_k^4, a scale mixture has
g_3 = g_4 to first order (the skewness of G^2 enters g_3 at relative order 0.5 sqrt(g) = 2%, not the measured 50%),
yet g_4 / g_3 is 0.50 at layer 1 and reaches 1.0 only at layer 15. The ledger below resolves it: the fresh
non-Gaussianity of ReLU is not a scale mixture, it is a skew, and the depth dynamics relaxes the skew towards the
invariant mixture direction at a rate that slows with depth.

### 2.2 The ledger

`code/gac3.py` works on the Monte Carlo pair cumulants of the official networks (`mc_cum_off{0,1}.npz`, T = 4e5; the
same accumulators as E10 and note XVIII). For each layer l it takes the true marginal state (mu, S) of y_l, forms the
Gaussian reference of h_l = relu(y_l) and its first variations under the third- and fourth-cumulant slices of y_l by
exact Gaussian integration by parts (`pairvar.Layer.var3`, `var4`: the chain's own programs), pushes the resulting
pair slices of h_l through W_{l+1} with the coherent restrictions

    kappa_3(z_k) ~ sum_a W_ka^3 k3_a + 3 sum_{a != c} W_ka^2 W_kc K21_ac,
    kappa_4(z_k) ~ sum_a W_ka^4 k4_a + 3 sum_{a != c} W_ka^2 W_kc^2 K22_ac,

and projects on the two patterns 1.5 mu_z sigma_z^2 and 3 sigma_z^4. Each channel is evaluated three times: with the
true slices of y_l (a ledger whose sum must reproduce the measured amplitude of z_{l+1}), with the Gaussian reference
alone (the one-loop generation phi_3, phi_4), and with the scale-mixture slices of unit gain (`pairvar.sm_slices`),
once the third-cumulant slices and once the fourth, which gives the 2 x 2 transport matrix M_l of a two-component gain
state. The ledger on network 0 (`outputs/gac3_off0.txt`; network 1, `outputs/gac3_off1.txt`, gives the same picture with
every transport entry within 0.05 and the amplitudes 5-10% lower at depth):

| l | g_3(l+1) measured / true-h / one loop phi_3 | g_4(l+1) measured / true-h / one loop phi_4 | phi_4 / phi_3 | M_l = [r33 r34; r43 r44] | row sums |
|---|---|---|---|---|---|
| 0 | 0.00642 / 0.00659 / 0.00661 | 0.00319 / 0.00317 / 0.00319 | 0.48 | [0 1.08; 0 1.09] (mu = 0) | - |
| 1 | 0.00952 / 0.00952 / 0.00441 | 0.00564 / 0.00558 / 0.00184 | 0.42 | [0.574 0.432; 0.119 0.922] | 1.006, 1.041 |
| 2 | 0.01132 / 0.01117 / 0.00309 | 0.00742 / 0.00725 / 0.00120 | 0.39 | [0.772 0.200; 0.140 0.879] | 0.972, 1.019 |
| 4 | 0.01399 / 0.01381 / 0.00244 | 0.01048 / 0.01020 / 0.00086 | 0.35 | [0.858 0.089; 0.125 0.855] | 0.947, 0.980 |
| 7 | 0.01755 / 0.01734 / 0.00204 | 0.01457 / 0.01391 / 0.00107 | 0.52 | [0.896 0.048; 0.053 0.860] | 0.944, 0.913 |
| 10 | 0.02076 / 0.02036 / 0.00225 | 0.01770 / 0.01639 / 0.00010 | 0.04 | [0.894 0.027; 0.073 0.819] | 0.921, 0.892 |
| 12 | 0.02198 / 0.02140 / 0.00197 | 0.02074 / 0.01832 / -0.00039 | -0.20 | [0.901 0.019; 0.063 0.799] | 0.920, 0.862 |
| 14 | 0.02232 / 0.02190 / 0.00131 | 0.02278 / 0.01960 / -0.00003 | -0.02 | [0.924 0.009; 0.016 0.817] | 0.933, 0.833 |

Four facts are in this table.

- **The first-order programs close on the pair slices.** In both channels the sum one-loop + variation under the true
  third-cumulant slices + variation under the true fourth-cumulant slices equals the true-h column to four decimals at
  every layer: the chain's integration-by-parts programs are exact at first order for what they carry, as note XVIII's
  mixture check found on the chain's states.
- **The pair-slice restriction loses nothing on the third cumulant and up to 14% on the fourth.** The true-h column
  is within 3% of the measured g_3 at every layer (2% at layer 15). For g_4 it is within 3% through layer 7, then
  5% low at layer 8, 10% at layer 12 and 14% at layers 14-15: the three slices of the fourth cumulant of h that the
  restriction drops (the (3,1), (2,1,1) and (1,1,1,1) patterns) carry a seventh of the coherent fourth cumulant at
  depth. The off-diagonal share of the pre-activation variance is -1% to -4% (the post-activations are slightly
  anticorrelated on average), so this is not the variance's off-diagonal part; it is the growth of the off-diagonal
  correlations themselves (rms 0.07 to 0.12 over depth, note XVIII) feeding the mixed patterns.
- **The fresh non-Gaussianity of ReLU is a skew, not a scale.** At layer 0 the pre-activation is exactly Gaussian with
  zero mean, so the amplitudes of layer 1 are pure one-loop generation: phi_4 / phi_3 = 0.483 on network 0 and 0.487
  on network 1, which is the measured ratio 0.498 / 0.457 at layer 1. The generation then falls with depth in both
  channels, faster in the fourth: phi_3 goes from 0.0066 to 0.0013 per layer and phi_4 from 0.0032 to zero at layer 9
  and slightly negative from layer 12 (the pre-activations' fourth cumulant created by ReLU at a nearly saturated gate
  has the opposite sign). The reason is the gate: as the mean direction saturates the gates, ReLU acts linearly on most
  of the mass, and a linear map creates no cumulants.
- **The transport is the mixture's, and it mixes the channels only while the gates are half open.** The row sums of
  M_l are 1.006 and 1.041 at layer 1, where homogeneity predicts 1 for the full cumulants, and they fall to 0.93 and
  0.83 at depth, the same deficits as the true-h restriction (2% and 14%) plus the fourth-order truncation of the
  integration by parts, which for a mixture whose fifth and sixth cumulants are also first order in g costs a few
  percent. The off-diagonal entries are the mixing of the channels by the nonlinearity: a fourth-cumulant gain feeds
  the third at 0.43 per layer at layer 1 and at 0.01 at layer 14; a third-cumulant gain feeds the fourth at 0.12 and
  0.02. Normalising the rows to 1 (homogeneity), the second eigenvalue of M_l, which is the survival per layer of the
  non-mixture component of the state (the skew), is 0.46 at layer 2, 0.72 at layer 4, 0.82 at layer 7, 0.88 at layer
  10 and 0.97 at layer 15. The relaxation towards the invariant mixture is fast while the gates are half open and
  nearly frozen once they are saturated.

### 2.3 The profile

With these pieces the depth profile of g_4 / g_3 is a two-component linear recursion, (g_3, g_4)_{l+1} =
M_l (g_3, g_4)_l + (phi_3, phi_4)_l, with the measured M_l or with its rows normalised to 1 (`code/gac3h.py`,
`outputs/gac3h_off0.txt`):

| l | measured g_3, g_4, ratio | measured M: g_3, g_4, ratio | rows normalised: g_3, g_4, ratio | one-gain (r33 + r34) g | critical sums phi_3, phi_4 |
|---|---|---|---|---|---|
| 1 | 0.00642 0.00319 0.498 | 0.00661 0.00319 0.483 | 0.00661 0.00319 0.483 | 0.00661 | 0.00661 0.00319 |
| 2 | 0.00952 0.00564 0.592 | 0.00957 0.00557 0.582 | 0.00955 0.00542 0.568 | 0.01105 | 0.01102 0.00503 |
| 4 | 0.01272 0.00899 0.707 | 0.01324 0.00911 0.688 | 0.01380 0.00892 0.646 | 0.01600 | 0.01682 0.00736 |
| 6 | 0.01559 0.01162 0.746 | 0.01621 0.01123 0.693 | 0.01819 0.01171 0.644 | 0.01934 | 0.02213 0.00879 |
| 8 | 0.01755 0.01457 0.830 | 0.01840 0.01205 0.655 | 0.02219 0.01445 0.652 | 0.02171 | 0.02688 0.01030 |
| 10 | 0.01967 0.01656 0.842 | 0.02028 0.01170 0.577 | 0.02635 0.01635 0.620 | 0.02366 | 0.03168 0.01053 |
| 12 | 0.02147 0.01956 0.911 | 0.02090 0.01028 0.492 | 0.02994 0.01804 0.603 | 0.02426 | 0.03587 0.01086 |
| 14 | 0.02231 0.02147 0.962 | 0.02069 0.00807 0.390 | 0.03309 0.01899 0.574 | 0.02387 | 0.03947 0.01000 |
| 15 | 0.02232 0.02278 1.020 | 0.02051 0.00690 0.337 | 0.03427 0.01923 0.561 | 0.02359 | 0.04078 0.00997 |

- **The shallow profile is derived.** Through layer 6 the recursion with the measured M_l reproduces the measured
  ratio to 0.05 (0.48, 0.58, 0.64, 0.69, 0.71, 0.69 against 0.50, 0.59, 0.66, 0.71, 0.75, 0.75) and both amplitudes to
  4%. The ratio starts at the generation ratio 0.48 and climbs because the transport mixes the skew into the
  mixture direction at the second-eigenvalue rate.
- **The deep fourth-cumulant amplitude is derived once homogeneity is imposed.** With the measured M_l the derived
  g_4 collapses after layer 8 (to 0.0069 at layer 15 against the measured 0.0228), because the 14-17% per-layer
  deficit of the restricted fourth-cumulant row compounds while the generation has stopped. With the rows normalised
  to 1, which is what homogeneity says of the full cumulants, g_4 is reproduced to 1-4% through layer 11 (0.01635
  against 0.01656 at layer 10) and is 8-16% low at layers 12-15, where the measured g_4 continues to grow at the rate
  of g_3 while the derived generation phi_4 is zero: at depth the fourth-cumulant amplitude grows only by being fed
  from the third (r43 g_3), and the restricted r43 (0.02-0.07) is the entry most affected by the dropped slices.
- **The third-cumulant amplitude saturates, and the recursion says why it cannot derive that.** The measured g_3 is
  0.02198, 0.02231, 0.02232 at layers 13-15: stationary. The measured-M recursion gets this right (0.0205 at layer
  15, 8% low) because its row sum 0.93 is a 7% per-layer loss that balances the one-loop generation 0.0013-0.0020; the
  homogeneity-normalised recursion overshoots by 54% because it forbids the loss and keeps accumulating the one loop,
  and the pure critical sum of one loops overshoots by 83%. Homogeneity is exact, so a true mixture cannot lose 7%
  per layer; what is lost is not the mixture's gain but the part of the one-loop generation that is not fresh. The
  gain inflates the pre-activation covariance along the mean direction by (g / 4) mu mu^T, which at depth is of
  order 5-20% of the off-diagonal covariance (g / 4 = 0.0055 times alpha_i alpha_j against correlations of 0.1); the Gaussian one loop computed on the
  marginal state re-counts that inflation as new third-cumulant generation. The earlier conditional-state check
  (`ncgprob/gain_defs.py`, note XVIII's record) measured the effect: the one-loop injection at layer 14 is 0.00074 on
  the marginal state and 0.00051 on the state with the gain's rank-one part removed, a third less. So the saturation
  of g_3 is the statement that at depth the genuinely fresh generation (two thirds of phi_3, about 0.0010) is cancelled
  by nothing, and the chain's measured 7% loss per layer is the re-counted third. A derivation that carries the gain
  on the conditional state, as GAC and `gac2.py` do for one channel, is required for the deep third-cumulant
  amplitude; section 2.4 reports it.

### 2.4 The two-channel recursion on the conditional state

`code/gac2b.py` is `gac2.py` of note XVIII's record with the fourth-cumulant channel added: the two gains are carried
on the closure's own conditional state (the mean divided by E[G], the covariance with the gain's rank-one part
removed, as GAC does), the generation and the transport are computed on that state, and the output mean is
E[G](g) times the closure mean, with E[G] from the Gamma law of note XVIII. No fitted constant. Four variants on
network 0: the measured transport or its rows normalised to 1, and E[G] fed by g_3 or by g_4
(`outputs/gac2b_off0_*.txt`):

| l | measured g_3, g_4, ratio | conditional state, measured M | conditional state, rows normalised |
|---|---|---|---|
| 1 | 0.0064 0.0032 0.50 | 0.0066 0.0032 0.48 | 0.0066 0.0032 0.48 |
| 2 | 0.0095 0.0056 0.59 | 0.0094 0.0056 0.59 | 0.0093 0.0054 0.58 |
| 4 | 0.0127 0.0090 0.71 | 0.0122 0.0089 0.73 | 0.0126 0.0087 0.69 |
| 6 | 0.0156 0.0116 0.75 | 0.0141 0.0107 0.76 | 0.0156 0.0111 0.71 |
| 8 | 0.0176 0.0146 0.83 | 0.0153 0.0111 0.73 | 0.0181 0.0133 0.74 |
| 10 | 0.0197 0.0166 0.84 | 0.0161 0.0102 0.63 | 0.0204 0.0144 0.71 |
| 12 | 0.0215 0.0196 0.91 | 0.0166 0.0088 0.53 | 0.0224 0.0154 0.69 |
| 14 | 0.0223 0.0215 0.96 | 0.0164 0.0064 0.39 | 0.0238 0.0151 0.63 |
| 15 | 0.0223 0.0228 1.02 | 0.0163 0.0054 0.33 | 0.0243 0.0150 0.62 |

| output of network 0, last-layer mean | MSE | residual scale (x 8) |
|---|---|---|
| closure (Gaussian, no gain) | 4.057e-6 | +0.0130 |
| GAC (one fitted-law gain) | 1.559e-6 | +0.0020 |
| one-channel derived gain (`gac2.py`, note XVIII record) | 1.541e-6 | -0.0011 |
| two channels, measured M, E[G] from g_3 | 1.525e-6 | +0.0008 |
| two channels, rows normalised, E[G] from g_3 | 1.525e-6 | -0.0012 |
| two channels, rows normalised, E[G] from g_4 | 1.670e-6 | +0.0034 |
| two channels, measured M, E[G] from g_4 | 1.972e-6 | +0.0052 |

- **The conditional state removes the double count.** With the rows normalised to 1, the critical accumulation of
  the third-cumulant channel on the conditional state gives g_3 = 0.0243 at layer 15 against the measured 0.0223,
  9% high, where the same recursion on the marginal state was 54% high. The one-loop generation on the conditional
  state at layer 14 is 0.00057 against 0.00131 on the marginal state; the difference is the gain-inflated covariance
  re-counted as generation, as section 2.3 diagnosed. The homogeneous accumulation is then the right law for the
  third cumulant to within the restriction's 2-3%.
- **The fourth-cumulant channel is the one that is still short at depth.** On the conditional state g_4 is derived
  to 0.03 in the ratio through layer 6 (measured M) and then falls behind on both states: 0.0150 against 0.0215 at
  layer 14 with rows normalised. The feed from the third channel, r43 = 0.02-0.08 at depth, is the entry the dropped
  (3,1), (2,1,1) and (1,1,1,1) slices of h affect most, and it is what carries the deep growth of g_4 to the level
  of g_3. So the approach of the ratio to 1 at depth is derived in its mechanism (rows summing to 1 and a dead
  generation make (1, 1) the stationary direction, reached at the second-eigenvalue rate) and not in its rate.
- **The derived scale is right to 1e-3 of the mean.** Both g_3-fed variants give a last-layer mean error of 1.525e-6,
  2% better than GAC's fitted gain and 1% better than the one-channel derivation, with a residual scale component
  of the mean of +0.0008 and -0.0012 (GAC +0.0020). The g_4-fed variants are worse because the derived g_4 is low.
  This is the closure-level estimator (no transported sources), so it is a test of the gain derivation, not of the
  chain; it says the mean's scale is a derived quantity to a part in a thousand, and that the quantity which sets it
  is the third-cumulant channel's gain on the conditional state.

### 2.5 What this says about the chain

The chain's fourth-cumulant sector is regenerated each layer as a transported diagonal plus a per-layer source
(g4row = 2 dG, dG = (W o W)(g_prev - var_prev lambda) + var lambda, note XVIII section 2): it has the form of the
fourth row of the recursion above, with the fitted table lambda_l playing the generation phi_4 + r43 g_3 and the
Hadamard transport playing r44. Note XVIII measured its amplitude at 0.77-0.83 of 3 g_3 var^2 at layers 5-14 and
called the fitted table the shadow of the g_4 / g_3 profile; this section derives the profile, and the table is
indeed its average over the layers that matter for the output. Replacing the table by the measured profile itself
(K4SM = 4, note XVIII) made the chain 12.5% worse, and by a flat 0.8 (K4SM = 3) 2.5% better on one network and 1.1%
worse on the other: the chain's sector is not the truth's slices alone, it also absorbs the dropped (3,1) amplitude
and the second-order gate terms, so the derived profile is a statement about the truth and not a knob. No chain run
was made for this section.

## 3. Non-autonomous and random thermodynamics

Two of the papers treat exactly the structure the depth dynamics has: a sequence of different maps or potentials
rather than one iterated map. Ben Ovadia (2609.12068) builds the thermodynamic formalism for potentials that vary
along an orbit of a driving Gibbs process: a semi-Ruelle operator, conformal measures, Gibbs estimates that hold
on average along the driving rather than uniformly, effective expansion on average, and a Lasota-Yorke inequality
giving quasi-compactness. Yayama (2605.28957) gives the existence of Gibbs measures for a sequence of continuous
functions, with the condition that the subshift be balanced.

**The dictionary.** The driving process is the weight sequence W_1, ..., W_16, i.i.d. He-Gaussian: the simplest
Gibbs process. A quantity computed for one realisation of the weights is quenched; averaged over the weights it is
annealed. The per-layer transfer operator of the pattern process is the sample-wise (semi-Ruelle) operator, and the
chain is its moment truncation; the finite-width hierarchies of section 6 are its annealed average. "Gibbs on
average along the driving" is the statement that constants hold for typical weight sequences, and we have the
number for how typical: the chain's raw error is 2.18e-8 on network 0 and 2.08e-8 on network 1, a 5% spread; the
ledger's generation ratio at layer 1 agrees across the networks to 1% (0.483, 0.487) and its transport entries to
0.05, while the deep amplitudes differ by 5-10% (g_3 at layer 15: 0.0223, 0.0200). At width 1024 the quenched
quantities self-average to that precision and no further.

**Where the chain is annealed and should not be.** The one object in the chain fitted once for both networks is the
fourth-cumulant table lambda_l. It is an annealed object standing in for a quenched one, and the quenched
fluctuation of what it stands for is measurable: g_4 at layer 14 is 0.0215 on network 0 and 0.0201 on network 1
(7%), and the derived replacement of note XVIII (K4SM = 3) moved the two networks' errors in opposite directions
(-2.5%, +1.1%). The per-network spread of the sector's amplitude is the size of what deriving it would gain: that
is the quantitative content of Ben Ovadia's "on average", and it closes the question of fitting the table better.

**Effective expansion on average.** Ben Ovadia's hypothesis is that the driven maps expand on average along the
driving. The network does not: the pair propagator is marginal in the bulk (He criticality; the transport of the
third-cumulant channel is 0.80-0.92 and that of the scale mode is 1), expanding only on the mean direction. The
quasi-compactness he obtains from Lasota-Yorke is what would give a uniform spectral gap of the transfer operator;
the gaps we can name are 0.08-0.20 per layer on the third-cumulant sources (summable memory) and 0 on the scale
(section 1). Yayama's balanced condition, uniform control of cylinder measures along the sequence, is what the
multifractal cell measure of note XIV section 5 denies: the Gibbs constants of the pattern process, if they exist,
degrade with depth. Neither paper has a finite-depth statement or a constant; both place the network outside their
hypotheses in the same way the Smale-space paper did.

## 4. The groupoid reading and the NCG state

Bissacot, Fukushima-Kimura, Pereira Lima and Raszeja (2308.16641) prove, for subshifts on multidimensional lattices,
that the Gibbs measures in the sense of the Gibbs inequality, the DLR measures, and the KMS states of the groupoid
C*-algebra of the Gibbs relation for the cocycle c_f built from the potential, coincide. This is the result that
links the Gibbs literature to the NCG language of notes XIV-XVII: a Gibbs state is a KMS state, the modular
automorphism group is the flow generated by the potential, and the DLR conditional expectations are the groupoid's
conditional expectations onto the unit space.

For the network the translation is complete and short. The measure is the Gaussian pushforward on the pattern
process; its "potential" against the product law is the conditional gate covariance (note XIV's Rota-Baxter defect).
The KMS condition for a Gaussian state is Gaussian integration by parts, E[y_i F(y)] = mu_i E F + sum_j S_ij E d_j F,
the modular flow is the Ornstein-Uhlenbeck semigroup, and the Mehler expansion in the correlation is its spectral
decomposition, which note XVIII section 1 already identified as the heat operator the chain truncates. The chain's
term programs are the KMS identities at first order in the cumulant slices and second order in the correlation;
section 2 of this note shows them exact at first order on the pair slices. So the groupoid characterisation is a
relabelling of what the chain is, and its one theorem beyond the relabelling, that the three notions of Gibbs state
coincide, is automatic for a measure with a smooth density. Nothing in it moves a number.

## 5. Product structure and constants

Climenhaga (2310.17495) proves that a measure with the Gibbs property for a Holder potential on a hyperbolic set has
local product structure: locally it is equivalent to a product of its conditional measures on stable and unstable
leaves. Thiam (2604.17528) gives five equivalent characterisations of Gibbs measures on subshifts of finite type with
explicit constants, through the Birkhoff cone contraction of the transfer operator.

**The pushforward has no product structure, and the scale is why.** The E4 measurements of note XVIII's record
(`mc_fibres0.log`) are the test. For the final pre-activation, the squared norm correlates with the squared
projection on a random direction at 0.109 +- 0.008 (a scale mixture gives the same value for every direction), on
the mean direction at 0.212 and on the top covariance direction at 0.214; conditioning on the norm removes the
excess kurtosis of random directions entirely (0.066 to -0.001), conditioning on the mean projection removes most of
it (to 0.008-0.019). A common scale across all directions is the opposite of a product of leaf measures, and on top
of it the expanding direction is coupled to the scale twice as strongly as the rest. So the pushforward is not a
Gibbs measure for the expanding/marginal splitting in Climenhaga's sense, which was already the conclusion of note
XIV section 5 from the multifractal cell measure; this is the same fact seen from the measure's side.

**No uniform constants.** Thiam's constants come from a cone contraction that is uniform over the shift. The
contractions we can name are the per-layer transport 0.80-0.92 of the third-cumulant sources (bulk), 1 on the mean
direction and the scale, and the second eigenvalue 0.46 -> 0.97 of the two-channel gain transport (section 2): all
depth-dependent, and the last tends to 1. The only uniform constant in the problem is the row sum of the
scale-mixture transport, which is 1 by homogeneity, and that is a conservation law, not a contraction.

## 6. The finite-width hierarchies are the annealed theory

Hanin (2204.01058) solves random fully connected networks as a perturbative hierarchy in 1 / n: the cumulants of
the output over the weight ensemble at fixed input obey recursions in which the fourth cumulant of the
pre-activations grows additively, kappa_4^(l+1) = 5 K^2 / n + kappa_4^(l) for ReLU, so that it is proportional to
sum_l 1 / n_l, and the activity norms obey a log-normal product law. Yaida (1910.00019) writes the same
recursions (R1-R3) for the four-point vertex, V^(l) / K^2 = 5 sum_{l' < l} 1 / n_l' for a single input. Celli
(2605.24072) gives Edgeworth expansions of the output law to order 4m - 1 with total-variation error n^-m, using
the conditionally Gaussian structure of the output given the previous layer. Kawase and Ota (2604.15742) write an
effective theory for the empirical kernel G with a stochastic recursion, a fourth-cumulant transport coefficient
chi and a tadpole source K_1, and find that a closure in G alone breaks down with depth as non-Gaussianity
accumulates, proposing to extend the state.

**Annealed and quenched, with numbers.** All four average over the weights. Under that average the third cumulants of
the pre-activations vanish by the symmetry W -> -W of the ensemble, and the output is conditionally Gaussian given
the previous activations (Celli's structure). The chain lives in the other ensemble: one network, random input. There
the third cumulant does not vanish, it is the chain's main content (kappa_3 = 1.5 g mu sigma^2, the Edgeworth mean
correction -(g / 8) sigma phi (1 + alpha^2) of note XVIII), and there is no conditional Gaussianity, only the
homogeneity scale of section 2. The annealed fourth-cumulant growth 5 L / n is 0.073 at layer 15; the quenched total
norm fluctuation Var |h|^2 / E^2 is 0.033 and its excess over the Gaussian part is 0.023 (E2): the annealed
fluctuation is two to three times the quenched one per network, the difference being the weight fluctuation that one
network has frozen. The chain's κ_4 ledger of section 2 is the single-sample version of Yaida's R1-R3 with the third
row added; its generation phi_4 falls to zero with depth where the annealed 5 / n is constant, because a fixed input
law saturates the gates along the mean direction and the ensemble average does not.

**The same breakdown, the same cure.** Kawase and Ota's finding that the kernel-only closure fails at depth as
non-Gaussianity accumulates is the quenched record of this programme in annealed form: the closure's mean error at
the last layer is 4.1e-6, GAC (kernel plus one global scale) 1.6e-6, the chain (kernel plus transported third
cumulants) 2.2e-8. Their transport coefficient chi is the r44 of section 2 and their tadpole K_1 is the mean-coupled
third-cumulant channel; their proposed extension of the state is the step from GAC to the chain. Their paper does
not have the quenched third-cumulant channel, and it does not have the conserved scale, because the ensemble average
hides both.

## 7. What transfers

| from the reading | what it is here | number or decision |
|---|---|---|
| g-measure against Gibbs (one direction, two directions) | the chain is one-sided by construction; the two-sided object is the n^4 bracket | the estimator's class is fixed; no two-sided representation to look for |
| summable variation and uniqueness | summable on every mode but the scale; the scale has constant variation (He criticality) | row sums 1.005 / 1.04 at layer 1; rho = 1 on the mean direction; the scale must be a global variable |
| Gibbs on average along a driving process | self-averaging of the quenched chain at width 1024 | 5% spread of the error between networks; ledger entries agree to 0.05, deep amplitudes to 5-10% |
| random potentials, non-autonomous Gibbs measures | the fitted lambda table is annealed; the quenched sector fluctuates by the size of the gain from deriving it | 7% spread of g_4; -2.5% / +1.1% on K4SM = 3; stop fitting the table |
| KMS states on the Gibbs groupoid | Gaussian integration by parts, Ornstein-Uhlenbeck, Mehler; the chain's term programs | relabelling; first-order exactness on pair slices verified to four decimals |
| local product structure of Gibbs measures | absent: a common scale across directions, coupled twice as strongly to the mean direction | corr(norm^2, proj^2) 0.11 random, 0.21 mean direction |
| explicit Gibbs constants (cone contraction) | no uniform contraction; the gain transport's second eigenvalue tends to 1 | 0.46 -> 0.97 over depth |
| annealed hierarchies (kappa_4 ~ 5 L / n, log-normal product law, conditional Gaussianity) | the quenched ledger: a third-cumulant row the annealed theory lacks, a conserved scale, generation that dies with depth | 0.073 annealed against 0.023 quenched at layer 15; phi_4 / phi_3 = 0.48 at layer 1, 0 at depth |
| kernel-only closure breaks down at depth | closure -> GAC -> chain | 4.1e-6, 1.6e-6, 2.2e-8 |

## 8. Sources

- A. van Enter, A. Le Ny, F. Paccaut, Markov and almost Markov properties in one, two or more directions, arXiv:2011.14664.
- S. Ben Ovadia, Thermodynamic formalism out of equilibrium, part II: semi-Ruelle operator, conformal measures, and effective expansion and quasi-compactness, arXiv:2609.12068.
- Y. Yayama, Existence of Gibbs measures for sequences of continuous functions, arXiv:2605.28957.
- R. Bissacot, B. H. Fukushima-Kimura, R. Pereira Lima, T. Raszeja, Gibbs measures on multidimensional spaces: equivalences and a groupoid approach, arXiv:2308.16641.
- V. Climenhaga, Gibbs measures have local product structure, arXiv:2310.17495.
- A. Thiam, Gibbs measures on subshifts of finite type: five equivalent characterizations with explicit constants, arXiv:2604.17528.
- B. Hanin, Random fully connected neural networks as perturbatively solvable hierarchies, arXiv:2204.01058.
- S. Yaida, Non-Gaussian processes and neural networks at finite widths, arXiv:1910.00019.
- L. Celli, Optimal non-asymptotic Edgeworth expansions for multivariate neural network outputs, arXiv:2605.24072.
- H. Kawase, T. Ota, Collective kernel EFT for pre-activation ResNets, arXiv:2604.15742.

Code: `code/gac3.py` (the two-channel ledger on the Monte Carlo state), `code/gac3h.py` (the closed recursions),
`code/gac2b.py` (the two-channel recursion on the closure's conditional state), with `code/pairvar.py` as shared
machinery. Outputs in `outputs/`.
