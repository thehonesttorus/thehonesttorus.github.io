# The birth address of the third-cumulant sources, and the precision the readout demands

Working note XXXI. Two theory-led measurements after note XXX closed the conditional-contraction cut: an
observable-sufficiency (quotient) test of the sources' internal address, and the oracle attribution of the chain's
error that note XXVII section 6 says the theory owes.

## 1. The birth address as a quotient

Each third-cumulant source is carried by legs whose rows are the current layer's neurons and whose columns are the
neurons j of its birth layer. The birth index is internal: transport multiplies the legs from the left, every
readout sums over j, and the D21 hub contracts over j, so a permutation or truncation of a source's columns at
birth commutes with everything the chain does to it. The development's Appendix B (an exact quotient exists when a
component is invisible to every later readout) therefore asks: which birth columns are invisible?

Every coefficient that weights birth column j in a readout carries neuron j's gate tail. The Wick coefficients w2_j
(the (2,1) birth term) and e_j, the source diagonal s_j, and the columns of M_b = diag(S3c) + 3 S21^T (the
post-activation cumulants kappa(y_j, y_j, .)) all vanish as alpha_j -> -infinity. The A and P columns themselves do
not (A = W a_b with a_b = diag(w1) C_off, P = W at birth). Prediction: a column's readout weight falls like a tail
factor of its birth neuron, so the columns of saturated birth neurons are an approximate quotient, removable at
birth by one permutation per source with no per-layer cost. Unlike the read-side drop of note XXIX (a saturated
neuron's own reads, which must stay), this removes nothing a saturated neuron reads.

**Kernel measurement** (`code/birthcols.py`, `outputs/birthcols_*.txt`; dumped legs, networks 0 and 1, read layers
5, 9, 13). Coefficient magnitudes against the birth alpha, relative to |alpha| < 0.25 (network 0, read layer 13):

| birth alpha | abs(w2) | abs(e) | abs(s) | column energy sum_i t_ij^2 | phi(alpha)^2 / phi(0)^2 |
|---|---|---|---|---|---|
| <= -2.5 | 0.014 | 0.0015 | 0.00015 | 7.6e-4 | 5e-6 |
| -2.5 to -2.0 | 0.11 | 0.019 | 0.0032 | 0.016 | 0.0068 |
| -2.0 to -1.5 | 0.27 | 0.069 | 0.019 | 0.077 | 0.049 |
| -1.5 to -1.0 | 0.52 | 0.22 | 0.099 | 0.25 | 0.22 |
| -1.0 to -0.5 | 0.80 | 0.49 | 0.33 | 0.53 | 0.58 |

The prediction holds: the column energy of the D3 readout tracks phi(alpha_j)^2 across the whole moderate tail, the
same on both networks and at every read layer. The dropped part of the readouts is small and incoherent
(correlation with the kept part -0.12 to +0.11 down to tau = -1.0):

| drop birth columns with alpha <= | columns dropped (layer 9 / 13) | D3: relative energy dropped | D21: relative energy dropped |
|---|---|---|---|
| -2.5 | 3-4% / 6-7% | 1e-6 - 7e-6 | 3e-6 - 1.3e-5 |
| -2.0 | 6-7% / 9-11% | 2e-5 - 6e-5 | 4e-5 - 9e-5 |
| -1.5 | 12-13% / 15-17% | 3e-4 - 7e-4 | 8e-4 - 2.4e-3 |
| -1.0 | 20-21% / 24% | 3e-3 - 5e-3 | 1e-2 - 2e-2 |

**In the chain** (`V36_BIRTH_SAT`, networks 0-15, paired against the V35 candidate; `outputs/ladder_*.txt`):

| drop birth columns with alpha <= | -2.5 | -2.0 | -1.75 | -1.5 | -1.25 | -1.0 | -0.5 |
|---|---|---|---|---|---|---|---|
| raw, legs masked at birth | -0.01% +- 0.10 | +1.85% | | +20.4% | | +158% | +755% |
| raw, masked inside the readouts only | | +1.75% | +7.1% | +20.1% | +60% | +158% | +754% |

Masking only inside the readouts (legs and joins untouched) does the same damage as masking the legs, so the joins'
basis selection is not involved: the readouts themselves cannot lose this much.

## 2. What the ladder measures: the readout's precision demand

The quotient is real (the weights fall as predicted), and still only alpha <= -2.5 is free, because the final MSE
is extraordinarily sensitive to the third-cumulant readouts. Pairing the two tables (incoherent dropped parts, all
layers at once):

| relative error energy of the kappa_3 readouts | about 1e-4 | about 1.5e-3 | about 1.5e-2 |
|---|---|---|---|
| change of the final MSE | +1.8% | +20% | +158% |

so, at small errors, **Delta MSE / MSE ~ 180 x (relative error energy of the kappa_3 readouts)**: a 1% rms error of
D3 and D21 at every layer costs about 2% of the final MSE, and a 7% rms error would account for all of it. This is
the quantitative form of note XXVII's "an incoherent 10% error of the bulk third-cumulant tensor accounts for the
whole residual", now measured with a controlled perturbation rather than inferred. It also says what any cheaper
representation of the sources must deliver: relative readout error energy below about 1e-5 to stay within 0.2% of
the MSE. The birth-address drop at alpha <= -2.5 meets that and removes 3-7% of birth columns (a few percent of the
young-tier bill); nothing more aggressive does.

## 3. Oracle attribution: where the chain's error lives, by statistic

Note XXVII section 6 owed an attribution of the chain's error to the statistics it carries. It is done here by
substitution: a Monte Carlo pass of 1.6e7 inputs through the network (`code/mcstats.py`, float32 with a per-layer
shift, two independent halves) gives, at every layer, the true pre-activation statistics the chain carries, and the
chain is rerun with one or several of them replaced by truth at every layer (`code/est_v37_audit.py`, switch
`V37_ORACLE`; `code/oracle_one.py`). The replaced value is noisy, so each oracle is run on each half and on the full
sample, and the noise-free effect is the extrapolation a = 2 MSE(N) - mean MSE(N/2) of MSE(N) = a + b/N
(`code/oracle_extrap.py`). Networks 0 and 1, the production configuration of note XXVIII (raw 2.2525e-8 and
2.1266e-8).

**The third-cumulant readouts and the fourth-cumulant diagonal** (`outputs/oracle_mc1_off{0,1}.txt`):

| replaced by truth at every layer | network 0 | network 1 | noise term at full N |
|---|---|---|---|
| D3 = kappa_3(z_i) | -14.8% | -15.0% | +2% |
| D21 = kappa(z_i, z_i, z_c) | -19.4% | -18.0% | +7-8% |
| D3 + D21 | **-40.0%** | **-38.6%** | +9% |
| g4 = kappa_4(z_i) | **-30.1%** | **-15.5%** | +1% |
| D3 + D21 + g4 | **-71.8%** | **-62.6%** | +10% |

The chain's readouts against truth (`outputs/oracle_cmp_mc1_off{0,1}.txt`): D3 has relative error 3.2-5.1% and
correlation 0.999 at every layer from 3 on (the Monte Carlo noise is 0.5-3% there and the chain is exact at layer 1);
D21 5-7%; the variance 4e-4 at layer 1 (noise) rising to 1.5e-3 at layer 15.

- **The kappa_3 readouts carry 40% of the error, and the measured precision demand predicts it.** A 4.3% rms
  relative error of D3 and D21 at every layer is an error energy of about 1.8e-3, which section 2's law
  (Delta MSE / MSE ~ 180 x energy) turns into about 33%; the oracle measures 39-40%.
- **The fourth-cumulant diagonal carries 15-30%, uniformly in depth.** Per layer, the g4 oracle lowers the MSE by
  25-30% (network 0) and 20-30% (network 1) at every layer from 3 to 15, so its defect is at full relative strength
  by layer 3 rather than accumulated. Note XXI section 7 had put the sector's leverage at the few-percent level by
  perturbing its fitted table; the oracle shows that the table's per-neuron residual (30% rms at depth, note XXVIII
  section 4) costs a fifth to a third of the error, the size the fourth Hermite coefficient
  (kappa_4/24)(alpha^2 - 1) phi / sigma^3 predicts for it.
- **The two are roughly additive and slightly synergistic.** Jointly they remove 72% and 63%, against 70% and 55%
  for the sum of the separate oracles: once the kappa_3 readouts are true, the fourth-cumulant error costs more
  (network 1: -24% marginal against -16% alone), which is the compensation of note XXI section 7 seen from the other
  side.
