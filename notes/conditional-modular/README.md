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

## 3. Measurement: representability on the reads (Azure VM, `outputs/`)

(running)
