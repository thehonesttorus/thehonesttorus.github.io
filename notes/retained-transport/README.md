# Retained transport: two-address-leg memory, the He readout metric, and the quartic pair class

Working note XXVIII. Three documents arrived from the companion synthesis: a 19-page report on the retained
third-to-fourth feed (`he_gaussian_retained_transport.pdf`), a theory note on conditional innovation memory
(`THEORY.md`: a shared address with one quadratic response matrix per neuron, and an exact fresh-He readout metric),
and a follow-up on shared quadratic actions (a canonical split of that memory and an SVD compression of its traceless
remainder). None of them had the chain, its source states or official weights; each ends by naming the measurement on
the chain's actual per-layer states that would decide it. This note does those measurements on official network 0,
on the Azure VMs, with the chain's own source legs and pair evaluator.

## 1. What each document claims, placed on est_v29

| claim | where it lands in the chain | status |
|---|---|---|
| exact retained feed K = 2 D_b S D_a + 2 D_a S^T D_b + diag, at an independent reference | the (2,2) pair program already carries the first-order d21T term at the correlated reference (Mehler order 2), and `pk22` is symmetrised, so both orientations of S enter K22 | already carried |
| the row-sum summaries S a, S^T b, diag S (a pooled scalar misses signed redistribution) | `g_prev = cA (K4 + K22 1) + cI (totals)`: the chain's fourth-cumulant core is fed by exactly these row sums | already carried |
| the full quartic weight dependence 3 (w^2)^T K w^2 - 2 diag(K)^T w^4, and its selected-gate O(mnk) evaluation | the chain replaces sum_b W_ib^2 K22_ab by a mean-field row sum carried by (W o W): the quenched per-neuron part is dropped | measured, section 4 |
| gate commutator [R, A]: the antisymmetric part of S is visible only through gate-ratio differences | an exact identity about which part of D21 the feed reads; the chain computes the feed from the full D21, so it is a diagnostic here, not a saving | recorded |
| two-address-leg family S_>=2(U) = Sym^3(U) + Sym^2(U) (x) U_perp with exact transport and merging | the chain's old tier is Sym^3(Qc) at rank 320 (every leg on the shared basis) | measured, section 2 |
| the He readout metric alpha \|E\|^2 + beta \|ctr E\|^2 and its closed-form optimal correction | applies to any projection whose family contains Sym(P (x) v); including the chain's own tier plus one n-vector | measured, sections 2-3 |
| shared quadratic actions: exact core and trace, SVD of the traceless remainder D^o | decides whether S_>=2 can be carried at the cost of a few actions instead of r(r+1)/2 coordinates | measured, section 3 |

## 2. The two-address-leg family on the chain's old tier (`code/s2fit.py`, `outputs/s2_L9.txt`, `outputs/s2_L11.txt`)

The old sources (age > 4, the chain's shared-basis tier; 5 sources at layer 9, 7 at layer 11) were dumped dense
(`V21_NO_CONFINE=1`) with the gate slope of every layer, written as a bank of stars
T = sum_h [3 w2_h J(A_h, P_h) + J(P_h, M_h)] (checked against the chain's D3 to 1e-16), projected at L0 in a frame U
(top-r left singular vectors of the stacked legs; a star-weighted frame gives the same numbers to 1%), and transported
exactly through the real gated weights G = W_(l+1) diag(w1_l). Scored on the reads the chain uses: relative error of
the D3 diagonal / the D21 slice. The projection, transport and read identities pass at 1e-13 and the He correction is
a verified minimum of its risk (selftest lines in the outputs).

Network 0, L0 = 9 (layers 9 and 14 shown; every intermediate layer lies between):

| r | Sym^3(U), the chain's kind | S_>=2 Frobenius | S_>=2 He metric |
|---|---|---|---|
| 32 | 0.509/0.563 -> 0.414/0.439 | 0.379/0.399 -> 0.327/0.334 | 0.184/0.279 -> 0.126/0.205 |
| 48 | 0.391/0.457 -> 0.300/0.334 | 0.256/0.272 -> 0.220/0.226 | 0.119/0.188 -> 0.081/0.133 |
| 64 | 0.300/0.375 -> 0.222/0.258 | 0.176/0.187 -> 0.151/0.155 | 0.079/0.127 -> 0.056/0.088 |
| 96 | 0.204/0.259 -> 0.136/0.170 | 0.086/0.091 -> 0.073/0.075 | 0.039/0.062 -> 0.027/0.041 |
| 128 | 0.145/0.187 -> 0.094/0.121 | 0.044/0.046 -> 0.038/0.039 | 0.020/0.031 -> 0.013/0.021 |
| 320 (production) | 0.019/0.025 -> 0.013/0.016 | | |

L0 = 11 (7 sources) repeats every column to within 0.01: at r = 64, Sym^3 0.225/0.282, S_>=2 Frobenius 0.134/0.140,
He 0.061/0.091; at r = 96, 0.153/0.197, 0.061/0.064, 0.028/0.042.

- **The He correction is real on the chain's content, and it is the larger of the two effects.** At every rank it
  halves the read error of the Frobenius projection (D3 by 2.2-2.4x, D21 by 1.4-1.8x), at every layer of the
  transport, with no fitted constant. The fresh-He readout law was derived for an independent next matrix; the
  actual quenched, gated transport respects it. The contraction vector ctr E is a first-class error channel of the
  third-cumulant memory, as THEORY.md section 6 says.
- **The free third leg is worth a factor of 3-4 at equal r.** S_>=2 at r = 128 (Frobenius) is where Sym^3 is at
  r = 320; with the He correction, S_>=2 at r = 128 equals the production tier on D3 and is 1.2-1.3x worse on D21.
- **Merging is exact.** All old sources live in one state (U, B_i); their birth labels are not needed. This is the
  structural difference from the chain, whose old-tier cost scales with the number of old sources.
- **Cost.** The full family costs about q/n units of 2 n^3 per layer to transport and as much to read, with
  q = r(r+1)/2: 8 units each at r = 128, plus 32 units per graduating source for the projection. That is more than the
  whole old tier (about 7-10 units per layer). The family is the right object; its r(r+1)/2 coordinates per neuron
  are the price, which is what the shared-actions note addresses (section 3).

## 3. Shared quadratic actions, and the He repair of the chain's own tier (`code/s2fit2.py`, `outputs/s2b_L9.txt`, `outputs/s2b_L11.txt`)

**The traceless remainder has no small action rank.** The He-corrected two-address-leg state was split canonically,
T = U^(x3) H + Phi_U(D) with U^T D = 0, D_i = d_i I / sqrt(r) + D_i^o, keeping the core H and the trace channel d exact.
The traceless remainder carries a large share of the tensor ((1/3)|D^o|^2 = 61 against a core of at most 167 at r = 64;
21 against 205 at r = 128), and its singular spectrum decays slowly and identically at every r and at both L0:

| s (shared actions) | 4 | 8 | 16 | 32 | 64 | 128 | 256 |
|---|---|---|---|---|---|---|---|
| fraction of \|D^o\|^2 left out | 0.90 | 0.82 | 0.67-0.69 | 0.46-0.49 | 0.23-0.25 | 0.06 | 0.003 |

On the reads, at r = 128 (L0 = 9): s = 16 gives 0.083/0.131, s = 64 gives 0.054/0.082, s = 128 gives 0.032/0.049,
against 0.020/0.031 for the full state and 0.019/0.025 for production. Matching production needs s of order 256, where
the source write (the p N r^2 term of the note's own cost formula, with N = 2n hub stars per source) is about 8 units
per graduating source: no saving against the tier it would replace. The algebra is right (core and trace exact, the
SVD tail is the exact local He risk) and the representation is honest; the chain's sources are close to the note's own
flat-spectrum counterexample, not to a small-tail case. Closed for cost.

**The He repair rides on the production tier, but there is little left to win there.** The family Sym^3(U) plus one
scalar-return vector {Phi_U(v I): v perp U} contains the He correction's K(v), so THEORY.md's formula applies to the
chain's own tier at the price of one n-vector (transport n^2, reads about r/(2n) = 0.16 units per layer):

| r = 320, reads D3/D21 | L9 | L11 | L14 | L15 |
|---|---|---|---|---|
| Sym^3 (production) | 0.019/0.025 | 0.015/0.020 (L0 = 11) | 0.013/0.016 | 0.010/0.014 |
| + He repair, U part only | unchanged | unchanged | unchanged | unchanged |
| + scalar return (Frobenius or He, equal here) | 0.014/0.022 | 0.011/0.017 | 0.009/0.014 | 0.008/0.012 |

A quarter off the D3 error and an eighth off D21, nearly free. But the whole representation defect of the shared-basis
tier is 1.4% of the raw error (note XIX, finite-volume round), so this is worth at most a few tenths of a percent and is
not implemented. At r = 128 the same channel takes Sym^3 from 0.145/0.187 to 0.098/0.156: useful only for a much
cheaper tier, which the action spectrum above rules out.

## 4. The quartic weight dependence of the kappa4 pair class (`code/k4mc.py`, `code/k4q_an.py`, `outputs/k4q_off0.txt`)

The chain's regenerated fourth-cumulant core transports its diagonal as t_g = (W o W) g_prev with
g_prev = cA (K4 + K22 1) + cI (totals): the row-sum (mean-field) form of the exact pair class

    Q_i = sum_a W_ia^4 K4_a + 3 sum_(a != b) W_ia^2 W_ib^2 K22_ab = 3 (w_i^2)^T K w_i^2 - 2 diag(K)^T w_i^4,

the PDF's full quartic readout. A Monte Carlo of 6e6 inputs (two halves for the noise) through network 0 gave, at
source layers 3, 6, 9, 12, the true per-neuron pair class of kappa4(z_(l+1)) by note XXI's power sums, the true total
kappa4, and the true post-activation slices K22 and K4; the production chain (current best, raw 2.271e-8 on this
network) was dumped at the same layers. The power-sum class equals the exact contraction of the true K22 to 1e-5, and
the chain's lambda part rebuilds at correlation 1.0000, so both sides are the quantities named.

| target layer | pair class (truth) | quenched part Q - Q_mf (truth) | chain's own Q - Q_mf vs truth | corr(chain exact Q, truth) vs mean-field | chain kappa4 residual rms (noise) | corr(residual, quenched) | LS coefficient | d(mean) per layer of the quenched part |
|---|---|---|---|---|---|---|---|---|
| 4 | 0.0111 | 0.00049 | corr 1.000 | 0.999 vs 0.939 | 0.00166 (0.00086) | +0.32 | 1.09 | 1.2e-5 |
| 7 | 0.0059 | 0.00026 | 1.000 | 0.999 vs 0.947 | 0.00107 (0.00032) | +0.29 | 1.18 | 1.4e-5 |
| 10 | 0.0032 | 0.00014 | 1.000 | 0.997 vs 0.944 | 0.00075 (0.00015) | +0.24 | 1.37 | 1.3e-5 |
| 13 | 0.0023 | 0.00010 | 1.000 | 0.994 vs 0.947 | 0.00083 (0.00009) | +0.14 | 1.13 | 1.4e-5 |

- **The mean-field row sum is right on average and drops a per-neuron part of 4-5% of the pair class.** The two means
  agree to three digits at every layer; the quenched part is purely incoherent, which is the signature of the chain's
  residual (closure round, section 4).
- **The chain already holds everything needed to compute it.** Its own slices give the quenched part at correlation
  1.000 with the truth: the part is fixed by the known weights, not by the slices' errors.
- **It is genuinely missing from the chain's kappa4 diagonal, at the derived coefficient.** It correlates +0.14 to
  +0.32 with the chain's per-neuron residual, with least-squares coefficient 1.0-1.4 against the derived 1. It removes
  1-5% of that residual's rms; the rest is the dropped (2+1+1) class and the regenerated lambda part (note XXI).
- **Its leverage on the mean is 1.2-1.5e-5 rms per layer** through the fourth Hermite coefficient of the gate,
  (kappa4/24) (alpha^2 - 1) phi / sigma^3. Fifteen incoherent layers of that size are about 2.5e-9 of output MSE if
  nothing compensates it, a tenth of the raw error. The chain test (section 5) is the arbiter.

Its cost is one n^3 product per layer, (W o W) K22, on the Strassen family: the PDF's selected-gate O(mnk) evaluation is
the cheaper form if the dense one pays.

## 5. The quartic pair class in the chain (`code/v33_k4q.diff`, switch `V33_K4Q`; Azure, networks 0-15, paired)

Implemented in `est_v29.py`: the previous layer's post-activation (2,2) slice and kappa4 diagonal are kept, and at the
next layer the exact pair class is formed with one n^3 product, (W o W) K22, on the Strassen family. Mode 1 replaces
the mean-field t_g everywhere (fourth-cumulant core, its (2,2) slice and the adaptive lambda rule); mode 2 adds only
the per-neuron quenched part Q/2 - t_g to the kappa4 diagonal and leaves the coherent sector, whose fitted amplitude
compensates other closure defects (note XXI section 7), untouched. Against the current best run on the same code:

| variant | raw | better on | C/B | adjusted |
|---|---|---|---|---|
| current best | 2.4134e-8 (mean of 16) | | 0.2177 | |
| mode 1: exact pair class everywhere | -1.92% +- 0.52 | 14/16 | 0.2266 | +2.09% +- 0.52 |
| mode 2: quenched part on the diagonal only | **-2.20% +- 0.31** | **16/16** | 0.2266 | +1.80% +- 0.31 |

The derived term lowers the raw error on every network at coefficient one, with nothing fitted: the full quartic
weight dependence of the PDF is a real, missing, correctly signed piece of the chain's closure, of the size section 4
predicted (a tenth of the error at most, a fifth of that realised). The dense product costs 0.65 units per layer (+4.1%
of the bill), more than the 2.2% it buys.

**Most of it is cheap.** Writing W o W = (row mean) + xi, the mean-field core contains every term of the exact pair
class except 3 diag(xi K22 xi^T) and the (4)-class weighting (W o W)^2 K4, which is n^2. The first is a quadratic form
whose size is set by the spectrum of K22, so a rank-k eigen-approximation of K22 carries it in proportion to its
spectral energy, and a randomized range finder (the chain's own recipe) costs about 4(k + 8)/n units per layer.
Measured on the dumped K22 (`code/k4lr.py`, `outputs/k4lr_off0.txt`):

| source layer | K22 spectral energy in top 1 / 16 / 64 | quenched part captured, k = 1 (corr, rel err) | k = 16 | k = 64 |
|---|---|---|---|---|
| 3 | 0.15 / 0.44 / 0.74 | 0.975, 0.22 | 0.983, 0.18 | 0.993, 0.12 |
| 6 | 0.32 / 0.64 / 0.84 | 0.988, 0.15 | 0.994, 0.11 | 0.998, 0.07 |
| 9 | 0.53 / 0.78 / 0.89 | 0.995, 0.10 | 0.998, 0.07 | 0.999, 0.05 |
| 12 | 0.63 / 0.84 / 0.92 | 0.996, 0.08 | 0.998, 0.06 | 0.999, 0.04 |

Rank one already reproduces the quenched part at correlation 0.975-0.996, even at layer 3 where it holds 15% of the
spectral energy. Mode 3 (`V33_K4Q=3`, rank `V33_K4Q_RANK`, 0 = the (4)-class term alone) adds this cheap form to the
kappa4 diagonal as mode 2 does: k + 8 eigenpairs from a sketch of the current weight, one power iteration, about
10 n^2 (k + 8) flops per layer.

| mode 3, against the current best | networks | raw | better on | C/B | adjusted |
|---|---|---|---|---|---|
| rank 0: the (4)-class term alone | 0-15 | -0.15% +- 0.05 | 13/16 | (gated with the dense product in this run) | |
| rank 2 | 0-15 | -2.30% +- 0.30 | 16/16 | 0.2185 vs 0.2177 | -1.94% +- 0.30 |
| rank 8 | 0-15 | -2.29% +- 0.29 | 16/16 | 0.2189 | -1.76% +- 0.29 |
| **rank 4** | **0-99** | **-2.25% +- 0.19** | **89/100** | **0.2186 vs 0.2177** | **-1.85% +- 0.19** |

- **The cheap form keeps the whole gain.** Ranks 2-16 give the dense product's -2.2% to -2.3% raw at a cost of
  0.0008-0.0012 of the budget instead of 0.0089; the effect is flat in rank from 2 on.
- **It lives in the pair part, not the diagonal class.** The (4)-class reweighting alone is worth 0.15%: the offline
  correlation at rank one is carried by K22's leading mode, whose quenched contraction ((W o W) v)^2 is the per-neuron
  content the mean-field row sum averages away.
- **New best on the 100 official networks:** mean adjusted 5.0143e-9 (was 5.1091e-9), mean raw 2.2938e-8, C/B 0.2186,
  with the current best flags plus `V33_K4Q=3 V33_K4Q_RANK=4` (`code/est_v33.py`, `outputs/v33_q3r4_100nets.txt`). The
  default path is bit-identical to the previous best (checked on networks 0-15).

So the retained-transport development delivered one concrete, derived improvement: the full quartic weight
dependence of the kappa4 pair class, which the chain's mean-field (W o W) row-sum transport had averaged away. Its
coefficient is one, nothing is fitted, it is better on 89 of 100 networks, and it costs 0.4% of the bill.
