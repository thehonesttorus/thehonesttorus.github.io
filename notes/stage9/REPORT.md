# Stage 9 session report: localization theory applied to the WhestBench mean-estimation problem

Branch `claude/determined-fermat-9hk45i`. Note: `notes/stage9/ncg_mlp_stage9.pdf` (16 pp). Code: `whest/`, `scripts/`, `infra/`.

## 1. Theory (the stage-8 principle made exact for estimators)

1. **Defect–path theorem** (Thm 1.2 in the note). Any estimator that is exact on point masses has error equal to minus the
   expected total drift of the estimator along any complete localization of the input law; the drift is the second
   variation of the estimator against the scheme's driving covariance. The other chat's "heat-flow defect" is the
   Eldan-path corollary; the theorem is scheme-independent and says which errors a scheme sees.
2. **Capacity of input tilts and the midpoint** (Prop. 2.1, Sec. 2.3). The merge coefficient is one half because the
   defect decays along Eldan's path like the kink density (1+tau)^(-1/2) (measured 0.96/0.73/0.46/0.28 vs 0.89/0.71/0.45/0.24
   at tau = 0.25/1/4/16): the parameter-free estimator is (e + Lap e)/2. The Euler–Stein defect delta(0) = e(0) - Lap e(0) is a 2% residual of an
   O(1) identity; the error is a tenth of that. Random probes need ~1e9 samples; k-direction tilts have signal-to-noise
   sqrt(k)*1e-2; only the full coordinate sum (2n+1 chain runs) resolves it. Input-side localization is a diagnostic,
   never an estimator component.
3. **Index-shift lemma and the Gram theorem** (Lemma 3.1, Thm 3.2). d_k' = d_{k+1}: a mean tilt raises the Wiener-chaos
   index; births, arms, slices and gate covariances are the first and second derivatives of the Gaussian moment map.
   Exactly: delta(y) = -<Hess G, Gamma(y)> with Gamma the Gram matrix of the first-layer moment gradients (closed
   forms). Verified to 7e-7 on a tiny network. The organisers' Wick coefficients are these index shifts, and their
   36-term table is the first-order tilt response.
4. **Scale-mixture lemma** (Lemma 4.1). A fluctuating collective scale produces a (2,2)-type kappa4 = 4v Sym(S x S);
   the Gaussian closure errs by -(v/2) sigma (alpha^2-1) phi and the kappa4 mean readout cancels it exactly at O(v).
   The organisers' chain keeps precisely this (one radial scalar per layer) and nothing else of kappa4.
5. **Sensitivity law** (Prop. 4.2): a relative variance error eps moves a mean by eps*sigma*phi/2; 1e-3 on the
   diagonal costs 7x in MSE; the bulk of the covariance may be off by percents; the collective mode is the central
   variable of the stage-8 split.
6. **Cocycle and the Hadamard obstruction** (Sec. 6). Old kappa3 sources are the images of their births under the
   gated transport cocycle, whose effective rank decays 164, 92, 65, 47, 29, 23 with age; but the slices are read
   through Hadamard products of the legs, so Frobenius truncation of the legs costs 2*eps in the slices; the memory
   is compressible only once it has become unimportant.

## 2. Measurements (official network 0, n=1024, L=16, truth from 1e9 samples)

| estimator | final MSE | cost |
|---|---|---|
| Gaussian closure | 4.06e-6 | 0.024 B (f32, symmetric einsum) |
| own kappa3 chain (k3v3) | 5.74e-7 | ~0.5 B |
| Gaussian closure + exact Euler–Stein defect, e - a delta (a = 0.535) | 4.35e-7 | +2049 chains (diagnostic) |
| Gaussian closure + parameter-free midpoint (e + Lap e)/2 | 4.50e-7 (net 1: 4.72e-6 -> 3.45e-7, a = 0.500) | same |
| organisers' K=3 simple chain, numpy port (`whest/kprop3.py`), radial kappa4 off | 5.63e-7 | |
| same, radial kappa4 on (one scalar per layer) | **3.53e-8** | 3.4e12 raw FLOPs = 1.5 B f32 (dense legs, rank 2n per layer) |
| organisers' chain + exact Euler-Stein defect, e - a delta (a = 0.255; corr 0.90, 81% explained) | **6.64e-9** | +2049 chain runs (diagnostic; symbolic version is the open task) |
| open-source chain (504aldo) | 2.13e-8 | 0.25 B |
| leaders (raw) | 1.1–1.5e-8 | 0.10–0.14 B |

- The port matches the torch reference to 1e-15 (per-layer means, covariances, (2,1) and (3) slices, radial core).
- Radial kappa4 scale: c*0.9 -> 8.9e-8, c*1.1 -> 2.7e-7, c*1.25 -> 4.4e-6 (c*1.0 = 3.5e-8 is the sharp optimum: no calibration gain there, and the chain is 2.5–8x sensitive to +-10%).
- Two-tier compression (young window w, shared basis k): (4,384 CP) **3.69e-8** (lossless to 5%); (4,256 CP) 5.2e-8; (2,384 CP) 7.9e-8; (4,128) 2.0e-7; (2,256) 3.2e-7;
  (2,128) 9.3e-7; (2,128)+collective directions 8.4e-7; (1,128) 2.0e-6; (0,128) 3.7e-6. Diagnostics (old sources kept dense):
  dropped (w=2) 3.7e-6; diagonal only 1.0e-6; scaled by 0.9: 5.1e-8; dropped (w=4) 7.4e-7; old slices replaced by their
  rank-1/4/16/64 truncation (w=2): 2.2e-7/2.0e-7/1.6e-7/7.8e-8 — the coherent (separable, O(n^2)) part carries most of the
  old memory; an incoherent tail remains at ages 2–5. The old sources' effect is smooth
  in amplitude and brittle in shape: Hadamard readouts defeat Frobenius truncation (note, Sec. 6).
- Euler-Stein of the window-1 chain (2049 runs, relay job lapx5c; the first runs lapx5/6 measured the randomized Tucker tier by a configuration slip and were discarded): raw 2.50e-6, corr 0.889, 81% explained, a = 0.482 (midpoint 1/2: 4.84e-7) -> **4.82e-7**. a returns to 1/2 for a first-order deficiency: third confirmation of a(d) = 1/(d+1). The correction of a short-memory chain is a factor 5, not the factor 50 of the memory itself: the dropped memory's effect is only 81% a defect.
- Source sparsification (keep the q*n most important hub columns of every birth block, no window): q = 0.5 -> 4.1e-7,
  0.25 -> 1.5e-6, 0.1 -> 2.6e-6; least-squares rescaling of the kept columns changes nothing. The importance profile
  is flat (every hub neuron contributes a comparable, orthogonal rank-one term), the non-sparsifiable regime.
- Billed cost under flopscope (fnp port `whest/kprop3f.py`, float32, MSE identical to float64): (4,128) 0.862 B; (2,128) 0.592 B;
  (1,128) 0.439 B; (0,128) 0.275 B; float64 doubles these. The accurate configuration (4,384) is ~1 B as written, ~0.45 B with the
  identity-leg structure of the birth blocks (3 transported matrices and 4 readout products per source-age instead of 6 and 6) and Strassen.
- Exact Euler–Stein defect of the ported chain on net 0 (2049 chain runs on Modal, 41 s each on 4 cores): corr 0.902, 81% explained, a = 0.255, MSE 3.53e-8 -> **6.64e-9**. The midpoint a = 1/2 does not apply (3.3e-8): a first-order chain's defect is second order in the births and decays like (1+tau)^-1, giving a = 1/3 in theory. Cross-network check on net 1 (job lapx4): corr 0.914, 84% explained, a = 0.257, 4.38e-8 -> **7.15e-9**; the constants cross over without loss (a_0 on net 1: 7.152e-9; a_1 on net 0: 6.646e-9). One constant a = 0.256 fixed offline corrects the chain 5-6x on every network tried.

## 3. Compute

- Modal app `whest9` driven through the AWS relay `w9-1` (`infra/modal_relay.py`: deploy / batch / fetch / run / log).
  The sandbox proxy cannot carry gRPC; the relay can. Up to ~90 concurrent containers observed at cpu=1–4.
  The relay stopped itself twice from idleness while a batch was running; the batch wrapper now keeps it alive.
- flopscope rules, measured (`notes/compute/FLOPSCOPE.md`): float64 bills 2x; `einsum('ij,jk,lk->il', W, C_sym, W)`
  with `as_symmetric` bills 3 n^3 (float32) for W C W^T; `stats.norm.cdf` returns float64 (cast back or the whole
  chain is billed at the float64 rate — this happened in the first billing run); the floor 0.1 B is 12.8 n^3 per layer.

## 4. What this says about a competitive system

- The accuracy of the exact first-order chain (3.5e-8) comes from two things: the exact first-order term table (36
  diagrams, all index shifts) and the one radial kappa4 scalar per layer (the central variable). Our own chain had the
  first in a different form and lacked the second, which alone was worth 16x. Everything else of kappa4 is noise.
- The cost of that chain is the memory of kappa3: 12 n^3 per source-age with the identity-leg structure (24 as the
  reference writes it), 0.5-0.7 B over the network in float32. The floor is 12.8 n^3 per layer in total. A window
  of four ages with a 384-dimensional shared basis is lossless (3.7e-8) and still costs about 0.5 B; two ages with
  that basis give 7.9e-8; anything cheaper loses an order of magnitude. Frobenius compression of young sources does
  not work because the readouts are Hadamard products (note, Sec. 6).
- Therefore a 0.1 B entry at the leaders' 1.1-1.5e-8 is not a compressed first-order chain. The measured smoothness
  of the old sources' effect in amplitude (10% scale -> 1.5x) and the open-source chain's memoryless regeneration of
  kappa4 suggested the same for kappa3, but the shape diagnostic (`scripts/diag_oldshape.py`) rules out the naive
  version: the old sources' slices correlate only 0.3-0.5 with the young ones (residual 0.9 of their norm) and grow
  to 3x their size at depth. The memory of kappa3 is a genuinely new shape at every layer; whatever the leaders do
  at 0.1 B, it is neither a compressed nor a regenerated first-order memory.
- The Euler-Stein merge is the one estimator-agnostic correction the theory provides: a factor 9-14 for the
  Gaussian closure and a factor 5.3 for the exact first-order chain (6.6e-9, twice better than the leaders' raw MSE),
  with one fitted constant whose value the theory predicts from the order of the chain's defect. Making it
  affordable needs the symbolic second-order response (note, Thm 3.2; task left open): that is the system to build.

### 2b. Four-network scan (nets 0-3, 240 runs on the relay; all rows agree to +-20%, net 0 quoted)
- Window ladder (old sources dropped beyond w ages): w=1: 2.5e-6, 2: 1.7e-6, 3: 1.1e-6, 4: 7.5e-7, 6: 3.3e-7, 8: 1.4e-7, 10: 6.3e-8, full 3.5e-8.
- Drop one age a from the full chain: a=1: 5.6e-7, 2: 4.0e-7, 3: 2.6e-7, 4: 1.8e-7, 5: 1.3e-7, 6: 1.0e-7, 7: 7.2e-8, 8: 6.3e-8, 9: 4.8e-8, 10: 4.2e-8, 11: 3.8e-8, >=12: 3.6-3.8e-8 (free). Marginal value decays x0.75 per age; eleven ages matter.
- No cancellation between ages (diag_agegram): at layer 15 the (2,1) slice is 15 contributions of norms 0.19..0.02 of the total, cosines 0.1-0.5, cumulative norm linear.
- Rank-r truncation of the summed old slices: window 4: r=1/4/16/64 -> 1.1/1.0/0.74/0.47e-7; window 2: 2.2/2.0/1.6/0.78e-7; window 1: 4.3/3.8/3.1/1.5e-7.
- Spectral birth truncation (star from top-k eigenpairs, residual exact): window 4, k=1/4/16/64: 1.3/1.2/1.0/1.0e-7 (saturates at 16); residual diagonal only: 5.4/5.1/4.0/2.1e-7; residual rank-k too: 5.2/4.5/2.3e-7. The missing factor 3 (w=4) / 6 (w=2) is the hub-diagonal pairing of the bulk, not its spectrum; the old sum is low rank because the cocycle is, not the births.
- Collective mode (diag_collective, 4 nets): top eigenvector of the pre-activation covariance: 0.9% of the energy at layer 0, 37-48% at layer 15; cos with the mean 0.97-0.98 from layer 12. Transported collective directions of sources born at layer >= 4 stay collective (cos 0.78-0.98); sources born at layers 0-3 have none. Family {T p_l'} at l=15: singular values (1, 0.32, 0.29, 0.23).
- Feature-learnability (diag_features): the exact chain's per-neuron error is not predictable from 34 local chain features (cross-network ridge makes it worse: 1.0e-6 vs 4.4e-8; in-sample -5%); collective share of its error 0.2-5%. The window-1 chain's error is 53-64% collective.
- Coherent + cross decomposition of an old source's readout (note, Prop. 5.2; verified exact to 1e-16 when all parts are kept, `oldmode=cross`):
  old sources (age >= w) read through: coherent only: w=2 1.6e-6 / w=4 7.3e-7 (= dropping them: 1.65e-6 / 7.5e-7);
  coherent + cross: 1.6e-6 / 7.3e-7; coherent + diagonal memory: 1.2e-6 / 5.4e-7; all three: 1.2e-6 / 5.3e-7;
  + star bulk: 7.3e-7 / 3.2e-7; + residual bulk: 3.1e-7 / 1.3e-7; dense: 3.5e-8. Same on nets 1-3 (+-25%).
  **The memory is the bulk hub pairing** (incoherent parts of the birth covariance and residual); the collective/mean-field parts,
  the only ones that accumulate cheaply, are worth less than a third of a decade. The mean-field tier (task 14) is refuted.

## 5. The eight papers and the essence of the obstacle (note, Sec. 8)

The user pointed to eight papers (hierarchic flows / Lempereur–Mallat; Hamiltonian sparsification and seminorm
sparsifiers / Basu–Brakensiek–Putterman et al.; samplizer / Wang–Zhang; adaptive phase estimation / Linden–de Wolf;
quantum Hermite transform / Jain et al.; LTF learning / Krivcenko–Nguyen–de Wolf; path counting via exterior algebra /
Panolan et al.). Read for the memory obstacle they give four laws and three leads:

- (L1) The birth blocks are sums of n rank-one hub terms of comparable, orthogonal weight: the non-sparsifiable
  class of Basu–Brakensiek–Putterman (Thm 1.5), confirmed by the sparsification scan (q=0.5: 12x worse; reweighting factor 1.000).
- (L2) No stochastic estimator reaches the 1e-2 slice precision: sample access costs eps^-2 (samplizer's quadratic gap,
  the O(k^2/eps^2) trials of the exterior-algebra estimator) = 1e4-1e6 probes vs 2n chain runs. Trilinear readouts
  are worse: Gaussian sketches of the hub index have zero mean on trilinear forms.
- (L3) Reorganising the computation (sequential/adaptive, promise of Gaussian weights) buys a constant factor at most
  (Linden–de Wolf: <= 2, with a Farkas dual certificate); our constant is the identity-leg/Strassen 24 -> 10.5 n^3.
- (L4) Crossover law: a source held in its cocycle range of rank r costs O(n^2 r) to transport but Theta(n r^3) to read
  out (the hub-diagonal readout is a Khatri–Rao core); subspace beats dense only for r < n^(2/3) ~ 101. The 99%-energy
  rank of the cocycle is > 256 for ages 1-3 and < 100 from age ~5 (diag_cocycle), so the window must be ~4 dense ages:
  exactly the lossless (4,384) configuration, 4 x 12 n^3 per layer = 3-4x the floor. An exact first-order chain cannot
  reach 0.1 B by any representation of its memory.
- (P1) Eldan's path is an OU semigroup, diagonal in Hermite degree (fast-forwarding): a defect of kink order d decays like
  (1+tau)^(-d/2), giving the merge constant a(d) = 1/(d+1): 1/2 for the Gaussian closure (measured 0.535, 0.500), 1/3 for
  the first-order chain (measured 0.255, 0.257; reduced by direction rotation). The constant is a property of the chain's
  order, fixable offline.
- (P2) The memory is a sum over hub-rooted stars (two input-to-hub paths, one hub-to-output path); the first leg
  accumulates into one n x n matrix (U_{l+1} = W D U_l + D_Phi T_{l+1<-0}); what cannot collapse is the hub-diagonal
  pairing of the identity leg with the readout leg of the same birth. Its mean-field part (variance profile of T o T)
  is O(n^2) and carries ~80% (rank-1 experiments); the rest is the realisation-specific fluctuation of a Gaussian
  matrix product.
- (P3) kappa3 memory, the second-order tilt response Gamma(y) and the Euler–Stein defect are one cocycle-transported
  identity-leg object (Gram theorem): the symbolic correction costs what the memory costs and must share its legs.

Conclusion: a floor-level system is a short-memory chain + O(n^2) derived corrections (mean-field hub pairing for old
ages, radial scalar, Euler–Stein merge with the derived constant). The deciding measurement: is the defect of a
window-one chain as explainable by its own coordinate Laplacian (Gaussian: 89-93%, exact chain: 81-84%)? That is
2049 runs on Modal (next).

## 6. Compute and hand-offs

- Modal: works through the AWS relay only (the sandbox proxy has no gRPC). `python infra/modal_relay.py batch JOB jobs/JOB.tsv --cpu N`;
  results in `s3://claude-whest-9292-97a992/results9/JOB/`, fetched with `python infra/fleet.py get JOB`. If the relay
  stopped itself, `python infra/fleet.py start w9-1` and `python infra/modal_relay.py fetch JOB` recovers results from
  the Modal volume. Concurrency reached 90 containers.
- AWS: spot and on-demand quotas are exhausted/limited (vCPU 200); only the relay runs. Quota increases must be
  requested in the console or the IAM user given `servicequotas:*` (see `notes/compute/SETUP.md`).
- Google Cloud: not set up; the browser-agent prompt is in `notes/compute/GCP_PROMPT.md`; hand the key back as the
  environment secret `GOOGLE_APPLICATION_CREDENTIALS_JSON`.
- Elicit: the three literature reports (localization schemes and estimators; Gaussian-conditioning MLP theory; the leaders' methods) are in `notes/stage9/elicit/`.
- Organisers' reference code: cloned read-only at `/home/user/alignment-research-center/mlp_cumulant_propagation`
  (public); torch 2.14 CPU installed locally for the verification only. Nothing from it is in the repo except the
  term table listing.

## 7. Stage 10: the selected Dirac operator (theory, `notes/stage10/ncg_mlp_stage10.pdf`, 20 pp)

The He-ReLU network with Gaussian input is placed inside random spectral geometry (Dirac ensembles), Sahasrabudhe's
two-level lift and Klartag's contact structure. Proved (with an identity-check script, `scripts/verify_stage10.py`):
- The weights are one random chiral Dirac chain D on the layered neuron set (quadratic Barrett–Glaser action; He =
  unit gain of the gated transfer). The input selects a projection P_x (its gate pattern); the network is the
  compressed resolvent z_L = [(1 - D_down P_x)^-1]_{L0} x. The activation tiling is the partition by the selected Dirac
  operator P_x D P_x; crossing a wall is a rank-one contact update.
- Two levels: walls of each layer are an isotropic random hyperplane process; codegrees on both sides are arc-cosine
  kernels; Crofton counts wall crossings (n x angular length / pi); at He, lengths are preserved and chords contract like
  3 pi / l (the selected geometry folds).
- Kink current (exact up to dead tiles, i.e. inputs where a whole layer is off, probability ~2^-n):
  u = sum over walls of f_{l,a}(0) E[|grad z|^2 (dh_L/dh_l) e_a | z_{l,a} = 0] = half the transported local time at the kinks
  of a Brownian path. Densities are easy (Gaussian right to 1e-3); the on-wall conditional transports are the hard part.
- The third-cumulant memory is a record of contacts: hub = kink, identity leg = pinned coordinate (Price). Frame bound:
  the memory is exactly as compressible as its transport leg TC. The neuron algebra and the modular frame of the layer
  state are complementary (||C|| ~ 2 sqrt(2/n)): spectral tiers are blind to hub diagonals. The memory readout is odd in
  the last weight matrix (symmetric-phase surrogates see none of it).
- Euler–Stein correction = average of the defect over Gaussian balls of uniform radius (OU time ~ Exp(1)); a defect of
  Hermite contact order d gives a = 1/(d+2) (d = 0: 1/2, Gaussian closure and window-1 chain; d = 2: 1/4, exact chain;
  measured 0.50-0.535, 0.482, 0.255-0.257). d is the order of the concentration-function coefficient at the kink that the
  estimator gets wrong.
- Mehler factors are imaginary modular times of the Gibbs state of the Hermite number operator; the ReLU heat trace is
  critical at linear order and its kink term beta^(3/2) fixes the folding exponent (beta_l ~ 9 pi^2 / (2 l^2)).
- Bootstrap: the mean is the unique value of a positivity + Stein moment problem; its first loop equation is the kink
  current; cumulant chains are its truncations without positivity (single-neuron level-one island = Scarf's interval).
The workflow agents' raw findings (readers of the three earlier papers, five lenses) are in `notes/stage10/workflow_findings/`.
