# Conditional modular flows and the matrix-trace carrier

Working note XXVII. Assessment of a pasted synthesis (its `THEORY.md` and code archive were links into another
sandbox and did not arrive): a hierarchy of conditional states and recovery maps, realised on a tower of algebras by
relative modular amplitudes, with the proposal that the unit of computation should change from transporting all
visible consequences to composing a small family of conditional state changes, and a concrete carrier for the
chain's third-order state. Section 1 checks the mathematics; section 2 says what it would change in the chain and
what has to be true for that; section 3 is the measurement that decides it.

## 1. The construction, checked

On nested matrix algebras A_0 in ... in A_J with a common trace and trace-preserving expectations E_j, a density h
restricts to h_j = E_j(h). The refinement amplitude v_j = h_j^{1/2} h_{j-1}^{-1/2} gives the query map
T_j(a) = E_{j-1}(v_j^* a v_j) and the state map R_j(b) = v_j b v_j^*. Checked:

- T_j is completely positive and unital (T_j(1) = h_{j-1}^{-1/2} E_{j-1}(h_j) h_{j-1}^{-1/2} = 1), and
  phi_{j-1}(T_j a) = phi_j(a) by the trace property; R_j(h_{j-1}) = h_j. This is the finite Accardi-Cecchini
  generalised conditional expectation and its Petz-type adjoint.
- The amplitude is Connes' Radon-Nikodym cocycle u_j(t) = h_j^{it} h_{j-1}^{-it} at t = -i/2 (positive definite
  densities, so the complex powers exist), with the twisted law u(t + s) = u(t) sigma_t^{j-1}(u(s)).
- The amplitudes telescope (v_J ... v_1 = h_J^{1/2} h_0^{-1/2}) and the query maps compose to
  h_0^{-1/2} E_0(h_J^{1/2} a h_J^{1/2}) h_0^{-1/2} without commuting anything; the relative entropy chain
  D(h_J || h_0) = sum_j D(h_j || h_{j-1}) holds because log h_{j-1} lies in A_{j-1}.
- The limitation is real and the example is right: for h = [[1, c], [c, 1]] over the diagonal,
  T(Z) = sqrt(1 - c^2) Z (h^{1/2} = [[alpha, beta], [beta, alpha]] with alpha^2 - beta^2 = sqrt(1 - c^2)); exact
  recovery of the reference state does not make T the identity on retained observables, so the acceptance test for
  any compressed carrier is the source-readout pairing, not state recovery. This is the same criterion as the
  source-query bound of note XXVI and the reading metric of note XXIV.
- The finitely correlated (CP-memory) recursion for additive readouts and the telescoping bound
  |phi_0(T_1...T_J a) - phi_0(T~_1...T~_J a)| <= sum_j ||(T_j - T~_j) a_j|| (contractivity of UCP maps) are standard
  and correct.

## 2. The carrier, and what it would change

**The object.** Everything the chain carries at third order is one symmetric tensor Theta_l (the sum over all live
sources, plus thin parts), evolving linearly between births, Theta_{l+1} = G_l^{(x)3} Theta_l + N_{l+1} with
G_l = W_{l+1} diag(w1_l) the product-gate response, and read only through D3 = diag(Theta) and the repeated-index
slices D21_ic = Theta_iic. The chain stores Theta as CP legs, n hubs per source, so transport and reads cost about n^3
per source per layer; that, not closure, is the bill (note XV: the cost law is the CP rank of the summed cores).

**The proposal.** Store Theta_ijk = Sym tau(a_i a_j b_k) with a_i, b_i real symmetric r x r (section 8 of the
synthesis). Three exact properties:

| operation | matrix-trace carrier | chain (CP legs) |
|---|---|---|
| transport G^{(x)3} | a_i -> sum_k G_ik a_k, b likewise: 2 n^2 r^2, for all sources at once | about 2 n^3 per source (young) or n^2 r_basis per source (old) |
| D3_i = Theta_iii | tau(a_i^2 b_i): n r^3 | Hadamard triple products, n^2 per source |
| D21_ic = Theta_iic | (1/3)[tau(a_i^2 b_c) + tau((b_i a_i + a_i b_i) a_c)]: n r^3 + 2 n^2 r^2 | hub contractions, n^3 per young source |
| product gate | a_i -> Phi_i a_i: n r^2 | row scaling of every leg |

(Symmetric matrices make tau(a_i a_j a_k)-type forms fully S_3-symmetric: cyclicity plus transposition.) At r = 32
transport plus reads is about 2 units per layer (one unit = 2 n^3), against roughly 14 for the chain at its current
bill, and the cost does not grow with the number of live sources. If the whole tensor could be carried this way at
small r with the chain's accuracy, the bill would fall towards the 0.1 B floor, where the score is raw x 0.1.

**What has to be true.** Three things, in order of difficulty:

1. *Representability on the reads.* A carrier of bond r has n r (r + 1) parameters, and its mode unfolding has rank
   at most r (r + 1). The chain's tensor has unfolding rank n (the birth core's C is full rank), so r >= 32 is
   necessary at n = 1024; the tensor itself is generic in that space (in the shared basis the old content is a sum
   of about ten sources times n hubs, CP rank far above 320^{3/2}), so the carrier can only succeed if the reads the
   chain makes, now and after future transports, see much less than the tensor holds. This is measurable on dumped
   tensors (section 3).
2. *Births.* Each layer adds a star core of CP rank n (legs C w1, I, C w1 w2 and the thin residual), which a bond-r
   carrier cannot absorb by direct sum without growing r. A rounding step (refit of carrier + newborn to bond r on the
   reads) must cost a few units per layer at most. The cubic module identity of note XXVI says what truncating the
   auxiliary space loses: compressing a_i by an isometry V changes tau(a a b) by the K^* D K terms, so the rounding
   must be done on the reads, not by projecting the matrices.
3. *Placement.* The cheapest first use does not need births at all: replace only the shared-basis tier (sources older
   than four transports, absorbed one per layer at the join) by a carrier, keeping the young sources as exact CP legs.
   That trades the old tier's per-source cost for one carrier, and the join becomes the rounding step.

## 3. Measurement: representability on the reads (Azure VM, `code/mt_fit.py`, `outputs/`)

Network 0, the chain run with every source dense so that the exact CP legs are available (`run_v29legs.py`, with the
gate w1 now dumped), legs at layers 10-13. The target is the summed tensor of either the shared-basis tier (the six
sources older than four transports at layer 10) or all ten sources, through the chain's own (3,) and (2,1) programs
(the rebuild matches the chain's D3 with correlation 0.998). It is transported exactly to layers 11-13 by
G = W_{l+1} diag(w1_l) with no births, so the comparison is the carried state's own future. The carrier
Theta_ijk = Sym tau(a_i a_j b_k) is fitted (Adam, float32) to the reads D3 and D21 (off-diagonal) of layers 10 and
11, transported exactly with the same G, and scored on layers 12 and 13. Yardstick: the chain's kind of compression,
every selected leg projected on one shared rank-D basis at layer 10 and transported exactly.

Relative errors of the reads, D3 / D21:

| carrier | fitted layers 10, 11 | held-out layers 12, 13 |
|---|---|---|
| matrix-trace r = 32, old tier | 0.000 / 0.055-0.060 | **0.22-0.20 / 0.35-0.32** |
| matrix-trace r = 64, old tier | 0.000 / 0.003-0.017 | **0.30-0.28 / 0.47-0.44** |
| matrix-trace r = 32, all ten sources | 0.000 / 0.074-0.079 | 0.28-0.25 / 0.44-0.40 |
| shared basis rank 320, old tier | 0.019 / 0.023, 0.017 / 0.021 | 0.015 / 0.019, 0.013 / 0.017 |
| shared basis rank 128, old tier | 0.13 / 0.16 | 0.11-0.09 / 0.13-0.12 |
| shared basis rank 320, all ten sources | 0.15 / 0.20 | 0.12-0.11 / 0.15-0.14 |

The carrier reproduces the reads it is fitted to (at r = 64 almost exactly) and not the state: transported two
layers on, its reads are 20-47% wrong, and more parameters make the held-out error worse. That is the source-query
criterion of the synthesis failing in its own terms: matching the present queries does not represent the content
that the future queries read.

The count says why. To 2% the old tier lives in the rank-320 shared basis, where its core is a sum of about six
sources times n hubs, CP rank about 6000 in 320 dimensions, far above the 320^{3/2} that a generic symmetric
tensor of that size needs, so it is generic: about 320^3/6 = 5.5 million degrees of freedom. A bond-r carrier
restricted to that subspace has about 320 r^2 (the a_i are n x r^2 but only their 320-dimensional span matters),
which matches only at r near 130, where its transport (4 n^2 r^2, about 32 units per layer) costs more than the
chain. The future reads do probe the whole core: D21 at one layer is n^2 numbers, and the rows of the propagators at
successive layers are fresh directions in the Lyapunov subspace, so a few layers of reads span it.

The over-determined repeat (dump at layers 9-14, fitted on 9-12, scored on 13-14) settles it: with four layers of
reads the carrier generalises to 9% / 13-14% (D3 / D21), the accuracy of a shared basis of rank 128 (10% / 12-13%),
against 1.3% / 1.7% for the chain's rank 320. And the chain's old tier is already near the information optimum: six
sources times three legs times 320 x 1024 is 5.9 million numbers against the 5.5 million degrees of freedom of a generic
core. The third-order content has no small algebra; the carrier is closed.

## 4. Integrated sampling, derived

The synthesis's central move is to change the unit of computation to conditional refinements that return an
integrated response: sample (or quadrature over) a retained address z, and integrate everything else analytically,
E[F] = E_z[T F (z)]. Feeding chosen values into a middle layer is the case where z is (part of) the state at layer l
and T is the downstream map. The question of whether this can beat the chain is exact, and it has an answer in terms
of four measured numbers.

**The precision ratio.** The output error the chain makes is epsilon ~ 1.5e-4 rms per neuron (white across neurons,
uniform in alpha, closure-round section 4); a single output neuron's spread over the Gaussian input is
sigma_F ~ 0.4 (plain Monte Carlo reaches MSE 1.2e-5 at the 10% floor with about 1.3e4 forward passes). Any estimator
whose per-neuron information comes from samples needs (sigma_F / epsilon)^2 ~ 7 x 10^6 effective samples, and the
budget allows about 10^4 forward passes in total. So sampling must be made 700 times more efficient per neuron
before it competes with the chain, and 2000 times before it beats it.

**Three ways to buy the factor, and what each costs.**

1. *Conditioning (Rao-Blackwell).* The variance left is Var_z E[F | z]. For every cheap address tried (1-8 input
   directions, a 256-direction first-layer frame, the radius, the antithetic pair) E Var(F | z) is at least 84-99% of
   Var F (closure-round section 6), so the factor bought is at most 1.2 and usually 1.01.
2. *Control variates.* A per-sample surrogate S with an exactly known mean reduces the variance to Var(F - S); a
   factor of 700 needs S to reproduce each output to about 4% of its spread on every sample. The surrogates with
   known means are low Wiener chaos (the first chaos of a deep ReLU network carries a fraction of the variance that
   decays with depth) or couplings to the chain's own closure law, and a coupling y ~ q_l of the true state x_l is
   only as close as the closure law is to the true law, which is the error being corrected: the construction is
   circular.
3. *Integration (each sample carries an analytic cloud).* An integrated sample costs a closure run downstream of the
   address. If that run is the chain, the budget holds one such sample. If it is a cheap conditional closure (mean
   field, O(n^2 L) per sample, so 10^3-10^4 samples are affordable), the address has to satisfy two conditions at
   once: E[F | z] may vary by at most sqrt(N) epsilon, i.e. about 3.5% of sigma_F (the address explains at most
   0.1% of the output variance), and given z the neurons must be close enough to independent that the cheap closure
   is accurate to epsilon (the address carries the correlations that make the product-gate closure fail).

**The obstruction.** Condition 3 asks for an address that carries correlation but not variance. In a
Gaussian-reference expansion that is impossible: the gate correlations rho_ab that the closure gets wrong are the
off-diagonal entries of the same covariance whose leading modes carry the output variance; the common-mode subspace
that removes 96% of the defect (note XX, rank 64) is the high-variance subspace. In operator language: the address
algebra would have to be in commuting-square position with the gate algebra (conditional independence of the gates)
while being nearly orthogonal to the readout (no variance), and for this network those two requirements pick the
same subspace. So, under the measured constants, no sampling scheme, adaptive or mid-layer, closes the factor; this is
the integrated-sampling version of the dichotomy in closure-round section 6, with the escape condition stated
explicitly.

## 5. The closure defect is a K_4 contraction

The term that integrated sampling would have had to supply, the all-distinct gate covariance of the third-order
transport (a third of the error, closure-round section 5), is fully determined by what the chain already carries (C
and the source legs); it is not missing information, it is an expensive contraction. Per source, at the next layer,

X_i = sum_{a, b, h} W'_ia W'_ib R'_ab Y_ah X_bh zeta_hi,   R' = D_phi R D_phi,   zeta_hi = sum_c W'_ic Phi_c Z_ch,

with h the source's hub index. Its index graph is the complete graph K_4 on {i, a, b, h} (every pair of indices
shares a dense factor), whose treewidth is 3, so every pairwise contraction order costs n^{tw + 1} = n^4 per source;
algebraic speedups give n^{omega(2, 1, 1)} (about n^3.25, impractical) or about n^3.8 with Strassen blocks. The
spectral split R' = rank-K common mode + bulk turns it into K extra leg transports per source (the cost note XX
priced), and the bulk part is incoherent, so it has no self-averaging (TAP-like) piece that a cheaper equation could
capture: the coherent parts are what the chain's existing terms already hold, which is why its residual is white and
unpredictable from local features.

## 6. Where an unlock has to come from

Taken together: the representation side is at its information optimum (section 3), sampling cannot buy the factor
(section 4), and the dominant identified defect is an n^4 contraction of information the chain already has
(section 5). The leaders' 1.6x better raw at lower cost therefore cannot be this architecture with a cheaper carrier
or a sampling correction. It has to come from an expansion that does not generate the defect, i.e. a reference state
and a set of carried quantities in which the gate algebra and the transported correlations stand in commuting
position, so that the product-gate transport is exact to the order carried. That is the precise form the
conditional-modular question takes here: find the address (the retained algebra at each layer) relative to which the
gates are conditionally independent and the conditional moments are cheap to carry, at a cost that does not scale
with the number of hubs.

Before building toward that, one measurement is owed by the theory itself, because the attribution it rests on is
partly by elimination: the closure round measured the first-order gate covariance at one layer (a third of the
error) and assigned the remaining two thirds to the same class (all-distinct fourth cumulant, second-order gate
terms) without computing them. Computing each dropped class exactly offline on one network (n^4 per source-layer is
affordable once on the VMs) and injecting it into the chain would show whether the classes the theory names account
for the error. If they do, the target above is the whole programme; if they do not, the unattributed part is the next
theoretical object.
