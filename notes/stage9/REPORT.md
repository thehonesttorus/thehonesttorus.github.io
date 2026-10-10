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

## 8. The token picture (stage 10, Section 9)

Input slots are vertices, weights are edges, tokens are Hermite quanta: Gaussian L^2 = bosonic Fock space over the
vertices, the k-th chaos = the k-token (bosonic Kikuchi) level, quadratic forms = elements of the pair-groupoid algebra,
contractions = groupoid convolution.
- Each layer-1 neuron is a one-mode token tower along its edge vector; its quadratic form is w w^T / (2 sqrt(2 pi) |w|).
- Euler–Fock ladder: for any 1-homogeneous F, (N + sum_i a_i a_i) F = F, so the mean of every neuron = 2 tr(its quadratic form).
- Layer 2: Q_c = W1^T D_c W1 (two-step walks i -> a -> j through the hidden layer); cumulants of the linear+quadratic part are
  closed and open walk sums on the hidden codegree graph G = W1 W1^T (log char. function = log det + resolvent). Exactly, distinct
  neurons contribute diagram sums (hub paths + triangles: 0.059 / 0.0115 vs leading diagrams 0.064 / 0.0133 at widths 24 / 96),
  while contacts (repeated neurons) carry most of a single pre-activation's kappa_3 and must be exact.
- Mean token number of layer-l pre-activations = 1/(1 - cos theta_{l-1}) (pair temperature): 1, 1.47, 1.98, ..., 13.0 at layer 16
  (width-256 net: 10.6). The tower does not truncate; Gaussianity is controlled by token contractions (codegree), not token number.
- Two tiers: E_{N(m,S)} F = <m, E grad F> + <S, E Hess F> (one token + two tokens); the heat equation = "a token pair equals two
  displacements"; a chain's heat defect is its violation; localization cools the environment (contact and exchange switch off).
- Two-token walk: Gamma -> W (K o Gamma) W^T with K = Phi Phi^T (independent) + diag(Phi - Phi^2) (contact) + exchange series.
  Independent walkers are O(1) wrong (two layers: 1.9e-1 vs Gaussian closure 1.7e-3) by exactly the contact term
  f_c sum_a W_ca^2 mu_a^2 (corr 0.999); removing it lands within 7.8e-4 of the Gaussian closure.
  => Gaussian closure = walk with exact contacts; kappa_3 chain = + first-order exchange; Euler–Stein = cooling correction.
- Sparsification: two-token objects and hub sums with hub coefficients known at birth accumulate (n^3 per layer; layer-2 stars
  cost one product W2 G). From layer 3 the hub coefficient is itself a walk (three transported tokens meet at the hub): that is
  the memory, and no operation on the input (which contracts input vertices, not hidden hubs) removes it.

## 9. Quantum-cut (Kikuchi) sparsifiers and a two-tier dynamics (stage 10, Section 10)

The paper (Basu–Brakensiek–Kothari–Putterman, arXiv:2606.09728) sparsifies L_G = sum_e w_e (I - swap_e), the Kikuchi Laplacian at
all levels at once. Manoeuvres: up/down operators and harmonic levels E_k; leverage of an edge decays with the level and beats the
dimension n^{k+2} in matrix Chernoff; inequalities proved once in the group algebra C[S_n] (Alon–Kozma via the octopus inequality)
descend to every level through the lift g(sigma) = f({sigma(1..k)}); expander decomposition, Chen's resistance lift, edge moving.
Bridges (tree edges) keep full importance at every level (their Remark 8.5).
- Essence: Schur–Weyl duality. (C^2)^{(x)n} = sum_k V_{n/2-k} (x) S^{(n-k,k)}; levels are su(2) weight spaces, L_G acts as the
  image of one group-algebra element in each two-row irrep, so one sampler serves all levels; k(n+1-k) = C(n,2) - content sum.
- The network has the bosonic counterpart: Howe duality (O(n), sl_2) on Gaussian space; E^- = (1/2) sum a_i a_i lowers token
  pairs, and the Euler–Fock ladder is E^- descent; the mean is the O(n)-invariant sector (radial foliation: leafwise sl_2,
  transverse O(n)). And the hard-core counterpart: gate patterns are points of 2^[n], g^2 = g is exclusion, the contact term is
  the hard-core correction, hidden relabelling S_n is an exact symmetry broken by the realization.
- Memory redundancy rho_eff = |sum P^(c)|^2 / sum |P^(c)|^2 = 1.00: hub sampling has relative error (1-q)/q at any width; every
  hub is a bridge. The paper's mechanism (importance decaying with level) has no counterpart in the readout.
- Higher tier with random dynamics: localization = stochastic flow on first-layer weights (contracting weights, Brownian bias).
  Gate probabilities are martingales driven by wall currents J = E_ball[delta(z) grad z]; the gate covariance (contact and exchange,
  within and across layers) equals the integrated Gram of wall currents (checked to 6 digits at layer 1 = Price's formula).
  Localization freezes deep gates first: disagreement of conditionally independent copies at layer l = theta_{l-1}(t)/pi
  (checked at width 512 to layer 12).
- The two tiers (depth tower of tile algebras; localization filtration) form a non-commuting square; the angle operator is the noise
  stability of the depth-l tiling, whose 2-point marginals are time-correlated Kikuchi weights running from Phi Phi^T (t=0) to K (t=inf).

## 10. Stage 11: tilings, Bratteli diagrams, singular foliations (`notes/stage11/ncg_mlp_stage11.pdf`, 17 pp)

Read in full: Bellissard–Julien–Savinien (tiling groupoids and Bratteli diagrams, arXiv:0911.0080), Julien–Savinien II
(arXiv:1005.2965), Putnam–Treviño (bi-infinite Bratteli diagrams, translation flows, arXiv:2205.01537), Giannakis–Montgomery
(measure-free Koopman–von Neumann, arXiv:2608.11591), Paul (semiclassics beyond Ehrenfest, hal-00617372), Francis (transverse
order k foliations, arXiv:2311.03940), Fischer–Laurent-Gengoux (neighbourhoods of leaves, arXiv:2401.05966), Louis (Nash blowup,
Helffer–Nourrigat cone, arXiv:2509.01133). Checks: `scripts/verify_stage11.py` -> `notes/stage11/verify_stage11.txt` (all pass).
- Cocycle: at fixed input the network is a weighted Bratteli diagram; activations = left-harmonic state, sensitivities =
  right-harmonic state, PT's conserved sum_v nu_r nu_s = <h_l, delta_l> = f at every layer (backprop); neuron = rectangle
  (width h, height delta), ReLU rescaling = diagonal (Teichmüller-type) flow; wall crossing = rank-one move. Rauzy–Veech induction
  is a gated linear network whose activation tiling is the Rauzy cylinders (checked: integer Jacobian det 1, 64/64 cylinders).
- Depth is parabolic: F(theta) = theta - theta^2/3pi - theta^3/18pi^2 + ..., Fatou coordinate 3pi/theta + (3/2) log theta,
  iterative residue 3/2; theta_l = 3pi/(l + 1.5 log l + 4.39 + o(1)) from orthogonal inputs. theta_16 = 0.378 exact, 0.384 Fatou,
  0.589 by the stage-10 law 3pi/l. Correction to stage 10: two independent inputs differ on 3.66 n gates over 16 layers, not 8.32 n.
  Farey (parabolic), not Gauss (hyperbolic); at infinite width sum theta_l = inf so gate sequences are never tail-equivalent;
  conjecture: finite width adds a linear contraction ~1/n (crossover depth ~ n).
- Two tiers = one flow: unresolved gate variance at (depth l, localization time t) depends only on the Fatou time
  tau = Phi(theta_0(t)) + l - 1 — PT's Teichmüller flow (state deformation + diagram shift) realized; one layer = sqrt2/3pi = 0.150
  of sqrt(1+t). Stage-10 width-512 data collapse on tau except the three last-layer points (finite-width signal).
- Walls = stratified singular foliation (leaves = faces). Same-layer crossings: normal crossings, abelian isotropy, Debord, HN cone =
  log cotangent. Deep wall crossing a wall that feeds it: bent by the gated transport T; formal transverse model = three-line
  arrangement y1 y2 (y2 + T y1) = 0, isotropy aff(1) ([E, theta] = theta), linear isotropy 1-dim (checked). Bent crossings = Price
  births of the third cumulant: kappa_3(z_2a) = c3 sum_b T_ab^3 |w_b|^3 for orthogonal rows, c3 = (pi+2)/(2pi)^{3/2}. Gaussian closure is
  exact for normal-crossing wall systems. Faces contractible, walls globally defined => no holonomy (FLG, Francis/Scott): the
  noncommutativity is isotropic and sits on the memory.
- Upper tier = semiclassical deformation: F(m,hbar) = Husimi symbol; exact boundary layer sqrt(hbar) kappa psi(d/sqrt hbar) (checked
  1e-8); target hbar=1 has n/4 unresolved gates per layer, so infinite width is a kinetic (van Hove) limit, not semiclassical.
  GM: moment chains = Fock truncations; depth -> commutative quotient (hbar_l ~ theta_{l-1}^2/2); u = kernel mean embedding, estimator
  error = MMD in the depth-L feature space.
- Mean = border functional on codim-1 faces = supertrace pairing of the wall flip with the normal derivative (checked 8 digits).

## 11. Stage 12: the network as its own two-tier system (`notes/stage12/ncg_mlp_stage12.pdf`, 20 pp)

Theory, with quick identity checks only: `scripts/verify_stage12.py` -> `notes/stage12/verify_stage12.txt` (all pass; one He
network, n=1024, L=16, 2e5 inputs, plus quadrature and small Monte Carlo identities).
- One layer is a two-tier system. Upper tier = quasi-free projection (mu_l, C_l) of z_l; lower tier = gate gas; up map = gating
  then W; closure = the second-order moment chain; memory = non-quasi-free remainder. Exact Euler split at every layer:
  E f = <E J_l, mu_l> + sum_a Cov(J_la, z_la); at l=0 the mean is pure lower tier (stage-11 border formula).
- Self-localization (proved at infinite width): through every row of layer l+1 the law of h_l is a ball whose standardized centre
  is N(0, t_l) across rows, t_l = |m_l|^2/Tr Sigma_l -> cos theta_l/(1 - cos theta_l) = Eldan time with input overlap theta_l.
  Rows = localization paths. Clock recursion t' = E m1^2 / E s1^2 = folding map. Lemma: E_{N(0,t)} Phi(r)Phi(-r) = arccos(t/(1+t))/2pi.
  Measured at l=16: t 11.84 vs 13.18, gate variance 0.0614 vs 0.0629, P(|r|<1) 0.226 vs 0.227.
- Depth = localization time: in the Fatou chart the clock advances by exactly 1 per layer; Szekeres' regular iteration gives
  fractional depth; Eldan time t = fractional depth Phi_F(arccos(t/(1+t))) - Phi_F(pi/2). Stage-11 collapse seen from inside.
- Cooling: gamma_{t,y} is the Gibbs state of |x|^2/2 at beta = 1+t; the network anneals itself on T_l = 1 - cos theta_l ~ 9pi^2/2tau^2.
  At finite width T_l ~ (1 - cos theta_l) + V_l/4 (angular + radial; V = across-input Var log|h|^2). l=16: 0.0779 measured vs
  0.0705 + 0.0080 -> 0.0779. Crossover tau* ~ 2.0 sqrt(n) (= 64 at n=1024; the competition nets are at tau_16 = 23.8).
- Hyperbolic upper tier: Fisher-Rao plane (curvature -1/2); field r = signed distance to the kink geodesic {mu=0},
  r = sqrt2 sinh(d/sqrt2); homogeneity = translation along it; rectification = skew product, base map R = m1/s1 > 0.
  Across-input log-norm increments (4/n)(3/2 - J2(theta)/2pi) ~ 4 theta^2/n are summable: n V_l -> 44.6 (30.8 at l=16), not 2+5l.
  Collective mode = this diffusion (trace share ~ tV/4: 0.093 vs 0.086 at l=16); stage 9's 47% is its Frobenius share.
- NEW (measured): the dominant source of marginal non-Gaussianity of deep pre-activations at width 1024 is the collective mode,
  not the one-body birth at walls (which is O(1/n)). Skewness gamma1 = eta * 1.5 V_{l-1} * r (R^2 0.86-0.91, eta = 0.72, 0.72, 0.71
  at l = 8, 12, 16); excess kurtosis flat in r at ~0.6 * 3V. The same eta fixes the field-compression prediction (10.09 vs 10.14).
  Since ReLU commutes with the radial scale, this part of the memory is a one-dimensional mixture, not a many-body effect.
- Gate gas: dichotomized Gaussian; Ising couplings J_ab = c_ab phi_a phi_b / (p_a q_a p_b q_b); Mattis component with pattern
  r phi(r)/(pq) ~ r|r| (the field itself). Wall-density law: any wall-localized one-body Q averages to (int Q) sin(theta/2)/sqrt(pi)
  ~ 2.66 (int Q)/tau -> harmonic in Fatou time, total ~ log tau (memory spread log-uniformly over depth). Hermite cancellation:
  neuron averages of marginal memory are suppressed by powers of T, mean squares are not (MSE sees it).
- Memory = Edgeworth series stratified by codimension (proved): each term is supported on chain strata; the cumulant on a stratum is
  the joint cumulant of the fields whose walls meet there (Kikuchi level = codimension). Normal crossings carry nothing; pair memory
  lives only on bent crossings (stage 11's aff(1) points). Third-order formula checked (softplus, difference 0.0018 +- 0.0028 and
  -0.0022 +- 0.0039); strata cancel strongly, so all three must be kept together. Birth of third cumulants ~ 0.718/tau.
- Modular: the Gaussian state is KMS on the fibre relation of h_l; extremal decomposition = law of h_l; conditional expectation onto
  the centre = posterior; centres L^inf(h_l) form a decreasing tower. Forward state localizes, backward state delocalizes as
  prod(1 - theta_{k-1}/pi) ~ (tau_l/tau_L)^3 (measured 0.0108 vs 0.0133 at the input, 0.5571 vs 0.5589 at l=12).
- Open: derive eta ~ 0.71; size of the bent-crossing pair term vs tau; radial-angular factorization at depth (integrate the radial
  mixture exactly, run the closure on the angular part); a modular (negentropy) bound on the closure error.

## 12. Stage 13: v56 read through convex geometry (`notes/stage13/ncg_mlp_stage13.pdf`, 21 pp)

Read in full: Klartag-Lehec (slicing, arXiv:2412.15044), Bizeul-Klartag-Lehec (KLS, arXiv:2610.05474), Klartag-Ordentlich
(SDPI under heat flow, arXiv:2406.03427), Klartag (stochastically evolving ellipsoid, arXiv:2504.05042v2), Sahasrabudhe
(exponentially small scales, arXiv:2512.15077), and the caustic-collar billiard paper. Checks: `scripts/verify_stage13.py`
-> `notes/stage13/verify_stage13.txt` (all pass).
- v56 anatomy: public 504aldo v29 chain (Gaussian closure, CP-leg kappa3 sources in young/old tiers, hub slices,
  memoryless kappa4, Wick table, Strassen billing) + our cost engineering + five theory terms: quenched kappa4 pair class,
  scale-mixture pair, the fold (-10% held-out), 103 counterterms (-3.8%), chaos-graded kappa4 diagonal / Schur hub (-9.8%).
  The two big wins are exact derived terms whose refitted amplitude stays at its derived value (0.95, 1.007).
- Eldan's cumulant hierarchy (BKL Lemma 4.3) = infinitesimal law of total cumulance (proved): drift -(m kappa_m + L_m),
  L_3 = 0, L_4 = three kappa3 A^-1 kappa3 pairings (checked 1e-7; 20% of the kappa4 drift). v56's Schur hub
  12[Y C^-1 Y^T]_ii is this between-tier fourth cumulant: explains C^-1, the need for quenched (source-resolved) kappa3,
  and weight 1. Predicts the (2,2) and (3,1) pairings from the same solve: L4(iijj) = <D_i,A^-1 D_j> +
  2<k3(ij.),A^-1 k3(ij.)>, L4(iiij) = 3<D_i, A^-1 k3(ij.)> (test on K31, 74% off).
- Critical gain = conserved dilation charge (relu commutes with scale mixing; transport factor exactly 1). Size law from
  stage 12: increments (4/n)(3/2 - J2/2pi) ~ 4 theta^2/n, saturating (V/4 = 0.0058-0.0072 at l=8-14 vs measured
  0.0070-0.0100; input radius 1/(2n)). Prediction: increments uncorrelated (saturation, not mean reversion).
- First-order cancellation (proved, checked O(eps^2)): the Edgeworth-corrected closure is exact on a scale mixture when it
  carries gamma1 = 6 eps r, gamma2 = 12 eps. So the charge moves the means only through inconsistency of the carried
  variance/kappa3/kappa4 -> build: a dilation-consistent (scale-mixture package) closure, O(n^2)/layer. Truth is ~70% package
  (kurtosis/skew-slope = 1.6-1.8 vs 2).
- Spherical geometry: f = h_{P+} - h_{P-} (tropical), E f = (V_1(P+) - V_1(P-))/sqrt(2pi) (mean widths; stage-11 border
  formula = edge formula of V_1); each layer evaluates the zonoid support function E<w,X>_+; closure = Gaussian
  (ellipsoidal) zonoid; zonoid curvature = wall density (all checked). Klartag dictionary: contacts = frozen gates, free
  directions = unresolved gates; principle: keep the exact constraints you touch (dilation package, Euler-Stein ladder).
- Billiards: support-function envelopes, g + g'' > 0 (= wall density), lambda = t^2 (= our T ~ theta^2/2), rigidity twin
  (closure exact on a full tilt family => Gaussian); Newton-on-high-modes needs family invariants we lack (flex theorem from
  the other side). No algorithm by itself.
- Slicing/heat flow: D(mu||gamma) = (1/2) int (n/(1+t) - E Tr A_t) dt (checked 0.5%) as a memory budget; the localized
  covariance diffuses with kappa3 as diffusion coefficient (quenched content).

## 13. Stage 14: lattices, leaves and ellipsoids (`notes/stage14/ncg_mlp_stage14.pdf`, 14 pp)

Theory with quick checks (`scripts/verify_stage14.py` -> `notes/stage14/verify_stage14.txt`, all pass).
- Klartag's process is spectral: an L-free ellipsoid = a flat torus R^n/L* with spectral gap >= 4pi^2; contacts = gap
  eigenfunctions (first spectral truncation; its operator system is spanned by the contact difference set); the process is a
  frozen-gap Dyson motion of the metric; it stops at a Voronoi-perfect metric (the truncation determines the metric); local
  optima are eutactic = John. Mahler compactness = the thick part. (A2 perfect + eutactic, Z^2 not: checked.)
- Leaves are sliding lattice points (proved): a wall {w.x = b} imposes Klartag's constraint <(x-m)(x)(x-m), A> >= 1 at every
  one of its points; the tangency point binds, giving w^T A^-1 w <= (w.m - b)^2 (a convex matrix-fractional "curved Ryshkov
  body"); at contact the constraint is <y(x)y, dA> = 0 with y the tangency direction. With the centre free and bias-free walls:
  the cone <w(x)w, mm^T - rho^2 Sigma> >= 0, rows = lattice vectors; deeper: a bundle over cells (gated Jacobian rows).
  The layer map z = Wx is Klartag's T (cells -> orthants, ball -> state ellipsoid): the chain's frame.
- The dilation is free for every contact configuration of bias-free walls (proved; dg - 2g = 0 checked), never frozen ->
  the lattice-side reason for the critical gain; biases break it (prediction: biased nets have a damped gain mode).
- Kink chaos law (proved, checked): E_k(r) = He_{k-2}(r)^2 phi(r)^2/k!; ~5-7% below the turning point r^2/4, none below
  r^2/8; beyond ~ phi(r)/(pi k^2 sqrt(4k - r^2)); tail phi(r) K^-3/2/(3 pi) (1.339e-6 vs 1.337e-6). So chaos truncation at K =
  the resolution ellipsoid of Mahalanobis radius 2 sqrt(K); next order removes 1-(K/(K+1))^1.5 (35% at 3->4) of the
  neglected chaos energy; at K=3 radius 3.46 (matches stage 12's "marginal memory vanishes for |r| >= 3").
- Leaf space: Gaussian transverse measures on every stratum; mean = codim-1 pairing with the kink cocycle; tile groupoid
  gap labels sum_C Z gamma(C) (cell masses). NC strata (bent crossings): isotropy aff(1) = the upper tier's location-scale
  algebra (E <-> scale); Aff_+(1) is non-unimodular (Delta = 1/sigma), so the Plancherel weight's modular flow is the dilation
  (Takesaki/Haagerup): the critical gain is the modular Hamiltonian of the NC strata. Radial partition function
  Z_l(beta) = int |h_l(x^)|^beta dsigma (its second cumulant = the dilation charge).
- Sticky localization with walls as obstacles: Ito accounting checked (18.700 +- 0.050 vs 18.683); contacts accumulate to
  perfect (6 in Sym_3). Maximal resolved probe = polar of the Loewner ellipsoid of dual wall points; precision
  A = sum lambda_a x_a x_a^T (polar eutaxy) -> whitened contact normals form a tight frame (1e-14): the local code.
  Conic Steiner formula holds in a layer-1 cell (statistical dimension 3.19 vs n/2 = 5). Siegel analogue = stage 12 field law;
  theta series n(1+st)^-1/2.
- For the chain: carry the dilation (cannot be pinned); effort follows the wall density (resolution r^2/4); positivity on the
  truncated operator system as a truth-free sticky constraint (diagnostic: PSD failures of per-neuron moment matrices).

## 14. Stage 15: deformed products, quasi-free lifts and the error cocycle (`notes/stage15/ncg_mlp_stage15.pdf`, 17 pp)

Read in full: Ebrahimi-Fard–Patras–Tapia–Zambotti (Hopf-algebraic deformations of products and Wick polynomials, arXiv:1710.00735)
and Zhou (non-commutative random surface growth with a reflecting wall, arXiv:2203.15920). Checks: `scripts/verify_stage15.py` ->
`notes/stage15/verify_stage15.txt` (official net 0 with its 1e9-sample truth at every layer; under a minute).
- Stance: start from the problem, not from the chain. Ask what any estimator's error is made of, which parts are predictable, and
  how the instance should shape the computation.
- Common mechanism of both papers: dynamics = convolution with a semigroup of functionals through a coproduct
  (phi_lambda = lambda * id; P_t = (id (x) <.>_t) Delta); classical process = commutative (Gelfand–Tsetlin) shadow of a free NC walk.
- Wall-jet theorem (proved, checked): under the TRUE law, the Wick (Appell) coefficients of relu(z - c) are P(Z > c), p(c), -p'(c), ...
  (the jet of the transverse measure at the wall); joint cumulants of relu outputs = diagrams with jet vertices and joint-cumulant
  hyperedges (no hyperedge inside one vertex). Connected law checked to 2e-13; Cov(relu, relu) on a non-Gaussian pair:
  Gaussian closure -1.2e-2, consistent expansion +4e-4 (30x), true jets with only k11 -2.7e-2 (worse): partial renormalization hurts.
- Rows are random gates (proved): Gaussian rows diagonalize Hermite degree (Mehler) = classical Pauli-path orthogonality; the radial
  (dilation) mode is the only obstruction.
- Error cocycle (leading order; measured on net 0): delta_L = sum_l T_L..T_{l+1} beta_l with T = Phi(r) (.) W (first jet).
  First jet carries 95% (Gaussian closure) / 79% (exact first-order chain) of the final MSE. Sum rule sum birth x gain vs MSE:
  Gaussian 2.59/4.06e-6, collective mode removed 1.50/1.57e-6; exact chain 3.20/3.53e-8. Collective share of inherited error 48-68%
  (Gaussian) vs 0-3% (exact chain). Transport factors 0.65-1.04 (near-critical). 81% of the exact chain's MSE is born in layers
  11-16 (42% in the last two). Last-layer channels (squared-jet fit): covariance 59%, kappa3 31%, kappa4 11%.
- No free lunch (proved): given the state, anisotropic births have zero conditional mean; only the isotropic part is predictable
  from per-neuron features (measured: 19% / 7% of the last-layer remainder). Explains stage 9's ridge failure and the -3.8% ceiling
  of v56's counterterms; on an exact first-order chain counterterms can gain at most ~2.5%.
- Lift: relu = (id + |.|)/2; the fold is the type-B reflecting wall (Weyl group B_n), births come only from it; Hermite even/odd =
  Laguerre(-1/2)/(+1/2) (Zhou's a = +-1/2 parities). Gate gas = arcsine pairs (Grothendieck–Krivine rounding of the quasi-free
  correlation) + star cumulants -(4/pi^2) sum_i prod_j rho_ij (checked: ratio 1.013/1.038/1.074) = the hub stars of stage 9.
- Complexity: Gaussian backbone = Clifford/matchgate; folds = magic, M_inf = 0.188, per neuron M_inf/sqrt(2 pi (1+t)) (checked
  1-9%). Generic quantum mean estimation for all n outputs ~1.1 B at MSE 1.56e-8 (MC 77 B) vs closures ~0.2 B: the resource is structure.
- Theory-native estimator: S1 modular sector exact (the one coherent mode); S2 linear cocycle exact (matvecs + adjoint gains);
  S3 births with true-law jets AND hyperedges at one order, allocated per instance by a fractional knapsack (additivity from the
  cocycle), Heisenberg (adjoint) exact reads of old kappa3 families at the late layers on near-wall rows (~100 n^3), short forward
  windows early. Predictions: exact modular sector takes the Gaussian closure from 4.1e-6 towards 1.6e-6; counterterms <= ~2.5% on an
  exact chain; ~60% of the last-layer birth is row-projected covariance error (test: one MC of Cov(h_15)); early-layer precision can be
  relaxed ~3x at <= 5% cost; biased nets shrink the coherent part.

## 15. Stage 16: the commutant tower, operator growth, transport memory (`notes/stage16/ncg_mlp_stage16.pdf`, 13 pp)

Built on the user's two companion notes (Branching, angular geometry and joint dynamics; Wick geometry, covariance holonomy and
Kikuchi closure) and on the experimental session's stage-13 tests. Checks: `scripts/verify_stage16.py` -> `notes/stage16/verify_stage16.txt`
(~6 min; `SKIP_SLOW=1` skips the network-0 Jacobian scan).
- Audit (all hold): spin branching p+- = (j+1)/(2j+1), j/(2j+1); three-spin logical qubit S12 = -Z, S23 = Z/2 + (sqrt3/2)X,
  [S12,S23] = -i sqrt3 Y; Zhou v1 Thm 4.3 must be read nested (E[Z_a^2 Z_b^2] = ab + 2a^2 + a: MC 3.021 vs 3.010, product 1.33);
  G_1 = diag(1/2, 1/2n, ...); squeeze holonomy of a contact loop -> a^2/4 (ratio 0.985 at a = 0.02); hard-core normalizer; mixture
  kappa4 = Cov(Q); Bernoulli Wick relation -2p(B-p).
- Commutant tower (proved, checked): O(n) -> Brauer (pairings), B_n -> even partitions, S_n -> partitions; dims 3 < 4 < 15 (2k=4),
  15 < 31 < 203 (2k=6). Annealed pre-activation cumulants are Brauer (odd ones vanish); upper-geometry mixtures give only pairing
  kappa4 (the dilation package); folds create the other blocks; three-block content of pre-activations is purely quenched, so no
  counterterm or mean-field tier can reach v56's remaining "bulk three-site structure". Annealed theory = noiseless subsystem of the
  O(n) twirl (multiplicity spaces); the instance lives in the spatial factor the twirl destroys.
- Operator growth (measured, net 0, 49152 inputs, Hutchinson JVPs): first-chaos share of neurons 0.73 -> 0.19 (layers 1 -> 16);
  mean degree of the higher chaos 2.76 -> 17.14, about +1 per layer. Fixed-basis (Pauli-path/input-Hermite) expansions are
  hopeless; per-layer Wick re-centring (frames that co-evolve with the instance's weights) is forced.
- Where the correction lives (net 0): 60-66% of truth - Gaussian closure is an isotropic function of each row's field; the exact
  first-order chain captures 99.1-99.7% of the correction; the residual is incoherent and late (matches the experimental session's
  0.4-3.7% along the mean).
- Transport memory (proved): the covariance chain leaves the frame rotation free; the normalizer theorem forbids diagonal-only
  (slice) closures under any transport outside B_n = Aut(Z^n) cap O(n) -> slices must be regenerated from transported legs (exact
  reason for stage 9's Hadamard/Frobenius obstruction). Lattice reading: consecutive coordinate lattices in general position.
- Ruled out: (i) the leaf/wall-sum estimator with global Dirac energies (carre du champ): exact at layer 1, MSE 3.5e-2 at layer 2
  (Gaussian closure 1.5e-7), 8.3 at layer 16; Dirac energy / variance of z = 1.00, 1.47, 2.54, 13.7 at layers 1, 2, 4, 16 (Poincare
  deficit = the chaos ladder); the exact identity relies on cancellation between layers, so transverse (leaf-conditional) energies would be needed.
  (ii) Minimal-sensitivity Wick frames: no gain over the centred (BPHZ) frame.
- Free and open: row-exchangeable self-calibration (theorem: an estimator can learn its own per-row correction from a random q-row
  subset of its own weights; added MSE (1 - q/n)(1 - R^2) E beta^2 x gain); non-Gaussian frame families; diagram selection by annealed
  two-copy variances; late-layer allocation. Suggested next measurement for the experimental session: R^2 of per-row second-order
  covariance births on their first-order proxy at layers 12-16 (decides whether self-calibration pays ~8x on the 59% channel).
