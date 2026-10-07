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
