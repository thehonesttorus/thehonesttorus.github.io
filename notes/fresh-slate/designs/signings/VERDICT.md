# signings stream — verdict (1 Oct 2026)

## The principle and what it turned out to do

Wick/Isserlis (Gaussian moments are hafnians), its Hermite–Mehler form (expectations of products of functions of a
Gaussian are loopless hafnians of half-edges weighted by correlations), and the signings/determinant toolkit
(Kasteleyn, Godsil–Gutman, covers, Bethe). Applied in the problem's own terms (DESIGN.md §1):

1. **Fresh-weight matchings fix the order of every sector** of the next layer's per-neuron cumulants. The layer weights
   are independent of the state, so each unpaired weight leg costs n^{-1/2}. Every loop (cycle) among weighted vertices
   is suppressed by an extra n^{-1/2}. Measured: the loop residual falls ≥ n^{-1.5} (T1). **The step sums are tree
   sums.** Determinantal resummation and signed averages have nothing to resum there.
2. **Signed averages are structurally blind to what is needed.** A Z/2 gauge average annihilates odd-degree content, and
   everything needed beyond the covariance is odd (third cumulants). A Z/3 gauge (cube roots of unity) sketches cubic
   CP content unbiasedly, but needs 1.6e3 → 3.5e8 samples for 10 % accuracy as the number of sources grows from 16 to
   1024 (measured). The infinite graph cover (Bethe) of the layered network is the independent-input model. Its "loop
   corrections" are the covariance and the old content, i.e. everything that matters.
3. **What the principle does give:** an exact closed-form generator of all joint cumulants from a Gaussian-copula state
   (T, R), with Mehler resumming all edge multiplicities (T0 exact). It also gives closed-form one-hyperedge (Edgeworth)
   corrections in the same basis, via the ReLU's δ/indicator profiles. Estimator v1 = copula + carried pairwise
   (2,1) slice + hyperedge corrections.

## Numbers

| | value |
|---|---|
| v1 raw final MSE, widths 64/128/256/512 | 1.75e-4 / 6.8e-5 / 1.41e-5 / 5.1e-6 |
| width slope | n^{-1.76 ± 0.12} |
| projected raw at 1024 (fit on 64–512) | 1.4e-6 (1σ 1.1–1.8e-6) |
| **measured raw at 1024** (bench `w1024_d16`, MLPs 0–2, N = 2e6, noise subtracted) | **1.10e-6** (1.17 / 1.21 / 0.92 e-6); paired Gaussian closure 5.06e-6 (ratio 0.22). Numpy wall ≈ 45 s per MLP |
| cost at 1024 (single-edge form, float32) | ≈ 18 units/layer (≈ 2 covariance, ≈ 6 cumulant contractions, ≈ 8 slice propagation, 2 coincidence corrections) → ≈ 190–260 units = 0.19–0.25 B; elementwise Mehler work < 0.2 units/layer |
| **adjusted MSE at 1024** | **≈ 2.4e-7** (measured raw × 0.19–0.25 cost factor: 2.1e-7 – 2.8e-7) |
| bar (1 Oct) | 1.6e-9 adjusted, so ≈ 150× short (raw needs ≈ 1e-8; v1 is 110× above). **Not competitive as built.** |

## Why, and what would make it competitive

- The single step is good and gets better fast: from an exactly specified state (marginals, covariance, (2,1) slice),
  one copula + hyperedge step has error 3.6e-5 at n = 128 and 2.55e-6 at n = 256 (≈ n^{-3.8}). If that holds, it is
  ≈ 1e-8 at 1024. The **local rule** is not the bottleneck.
- The bottleneck is **drift of the carried state over depth**. The state v1 carries does not contain the third-order
  content born at earlier ReLUs (the all-distinct slice: ≈ 70 % missed by the copula's tree paths at layer 3), nor
  second-order pairwise structure. Its marginals (means, variances, shapes) drift layer by layer. Carrying diagonal
  births as CP sources (v2) recovered only ≈ 20 % of the missing all-distinct sector and was unstable at width 64.
- To be competitive, the state must carry hyperedge content across depth: sources born at every ReLU (diagonal and
  star-shaped/path-born), transported by gate-linear chains. Their contractions with fresh weights are matrix products,
  and the cost is L²/2 products unless old ages are compressed. That is a chain of cumulant sources. The matching
  principle derives and organises it, but makes it no cheaper. Signings do not touch it.

## The deciding experiment

Teacher forcing at n = 512 and 1024 with high-N per-layer statistics, run two ways:
(a) one step from the exact state, to confirm that the ≈ n^{-3.8} decay holds;
(b) a single injection at layer k, then the measured drift per subsequent layer as a function of n.
If (a) reaches ≈ 1e-8 at 1024 and (b) shows the drift is carried by a low-dimensional set of old sources, a copula +
hyperedge step paired with a compressed source bank is worth a Stage P. If the drift needs ≳ 0.3 n modes (as the brief's
old-content fact suggests), it is not, at ≈ 0.1 B.

## Failures charged to the realisation, in order

1. The v0 copula as the whole state: right for one step, wrong from the second. The (2,1) defect alone gave 3 % next-layer
   variance errors at width 64.
2. The first v1 propagation of D missed the coincident-index tree terms (13 % slice error at layer 2). Found by an
   exact-copula check and fixed.
3. v2 sources with all ages and no compression: unstable (indefinite Cov(a)) at width 64.
4. The determinantal/signed machinery found no place in the estimator. This is a negative result about where loops live
   in this problem (in the state, across depth), not a defect of the theory.
