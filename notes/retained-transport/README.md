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
