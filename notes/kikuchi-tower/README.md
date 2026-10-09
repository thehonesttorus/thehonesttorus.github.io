# The Kikuchi tower: the MLP as an interacting token walk, and what each move costs

Working note XL. It takes up the proposal of this round: read the MLP as a walk, a random process, on a Kikuchi graph,
generalised to the noncommutative setting of the attached study (the symmetric-group AF algebra, its Young graph and
traces, the translation crossed product, coefficient localizers and their Schur reduction, KMS transfer, strong
minimality). The aim is a cost theory first: which parts of the cumulant tower can ride the chain's existing legs, which
cannot at any price below n^4, and what that licenses in a FLOP-aware system. The predictions in section 5 are stated
before their runs.

## 1. The tower

**Tokens.** The k-th cumulant of a layer's pre-activations is a symmetric function on [n]^k: a vector in the k-token
sector Sym^k(R^n) of the bosonic Fock space over the neurons. A basis vector is a multiset of sites. Its occupation
pattern is a partition lambda of k (the multiplicities, sorted):
- the hard-core pattern lambda = (1^k) is the Kikuchi level k (k-subsets): the bulk;
- the patterns with a part >= 2 are the chain's slices, (2,1), (3), (2,2), (3,1), (2,1,1).

**The Young structure.** The permutations of the sites act on every pattern sector. For a (2,1) slice (an n x n table
off the diagonal) the module splits exactly into:
- the trivial irrep (the constant S0);
- two standard irreps (the row and column effects, S1 and A1);
- the higher irreps (zero row and column sums, S2 and A2).

This is the five-component split of note XXXIX section 5, now recognised as the isotypic decomposition. In the
attached study's terms:
- K0 of the symmetric-group AF algebra is Lambda/(p1 - 1), with the Schur classes s_lambda as the pattern labels;
- its traces are the S_n-invariant (annealed) statistics, the trivial-irrep component;
- the standard irrep is the additive (rank-2) part.

**The layer as a token process.**
- z = W y moves every token independently: W^(x k) on the k-token sector, the second quantisation of W. It preserves
  the CP form of a source: each leg is transported on its own.
- y = f(z) acts by the cumulant (Faa di Bruno / Wick) expansion. Each output token sits at a vertex v, which absorbs p
  internal lines from input hyperedges (cumulants of z on the output sites) with weight E[f^(p)(z_v)]. This is a
  branching, interacting process on the tower.

Its elementary moves, by how many tokens of one CP term they touch:
- **V0, vertex dressing.** A vertex also absorbs a local hyperedge of its own site (D3_v, kappa4_v), a tadpole. The
  weight c(p) becomes c(p) + kappa3 c(p+3)/6 + kappa4 c(p+4)/24.
- **V1, one-token moves.** A vertex joins one token of a term to a covariance or slice hyperedge. Two kinds:
  - the newborn's star is created this way (Phi C w2 C Phi);
  - GC1/GC2 move a token out of a doubly occupied site of the total (2,1) slice, along a covariance edge: a Kikuchi
    move on the covariance graph.
- **V2, two-token edges.** A covariance edge joins two vertices that carry tokens of the same CP term, for instance
  the gate covariance w2_a w2_b Phi_c C_ab T_abc of note XXXV section 5.

## 2. The cost classification

**Theorem 1 (moves and CP rank).** Let a source be T = sum_r u_r x v_r x x_r with independently transported legs.
Then in one layer:
- (i) **V0 rescales legs.** CP rank unchanged; O(n) per layer.
- (ii) **V1 displaces one leg per term.** The CP rank is unchanged and the move costs O(n^2), when the move acts on a
  term whose touched token is localized (the newborn's centre, P = I) or reads a slice of the total tensor (the GC
  class, folded into the newborn's arm, note XXXIX section 11).
- (iii) **V2 is a wall.** The pair state of the two joined tokens becomes d(u) C d(v), which is full rank whenever C
  is. Its exact transport needs W d(u) for every term: n^3 per term and n^4 per source, in any representation that
  transports the pair jointly. This is the Khatri-Rao wall of notes XV and XXXV, now as a statement about the move
  class rather than about one contraction.

*Proof sketch.*
- (i) A tadpole multiplies the vertex weight by a per-site scalar.
- (ii) For a localized token the joined line is a vector, so the term keeps its product form with one leg replaced
  by an n-vector computable in O(n). For the total-slice reading, the moved token's new leg is a row of an n x n table
  the chain already holds, so the displacement is O(n^2) per source (note XXXIX section 11b).
- (iii) The Schur product (u v^T) o C = d(u) C d(v) has the rank of C for legs without zeros. Under the next linear
  layer it becomes W d(u) C d(v) W^T, and the only shared factor across terms is C itself.

**The inventory of the kappa3 bulk at first order, by class.**
- Transport Phi_a Phi_b Phi_c T (free tokens): carried.
- The star (V1, localized centre): carried.
- GC1 and GC2 (V1, total slice): carried, folded.
- The triangle C_ab C_bc C_ca (V2 at birth): measured negligible (note XXXVI section 3e).
- **The gate covariance (V2): the wall.**
- **The tadpole dressing of every vertex that transports or creates bulk (V0): free, and missing.**
  - The pair programs dress the covariance ('d3row' / 'g4row' x 'c_off').
  - The legs, the births and the replicated slices (D21_w = w1^2 D21 w1, D3_w) use the bare Gaussian w1 = Phi.
  - Under the scale (gain) mixture the dressed first vertex is P(z > 0), which the mixture leaves invariant; the bare
    Phi(mu/sigma) is not. The omission is therefore a gauge inconsistency, not just a missing term.
- The (2,1,1) fourth-cumulant feed: a closure (the lam core), flex-limited (notes XXXVII-XXXVIII).

## 3. The noncommutative dictionary

- **The Young AF algebra and its traces.** The token patterns and their annealed statistics, as in section 1.
- **The translation crossed product A x| Z.** The layer shift. A source's birth index is the Z-grading and its age
  the winding.
  - "No new even trace values, a whole new odd group": transport never changes a source's annealed totals but
    carries its orientation, its history.
  - "40% of the (3,1) slice is history, not state" (note XXXVIII) is the measured face of this.
- **Schur elimination and the polynomial inverse (the study's Theorems 4.7 and 4.8).**
  - The birth M block is an exact elimination of the slice sector into the newborn, up to its rank-4 compression.
    The fold's lesson (P9: +140% without the residual correction) is the theorem's hypothesis made concrete: the
    elimination must be exact on the gap-sensitive part, which here is the low-irrep part.
  - The young window is the finite-propagation polynomial approximation of the kappa3 transport's resolvent. Its
    certified error is q^(k+1), and note XXXVI measured q = 0.46 per age.
- **KMS transfer and cohomology.**
  - The readout-weighted (adjoint) measure on token histories is the conformal measure.
  - Cohomologous potentials have diagonally similar transfer matrices. In the chain this is the ReLU network's
    positive-diagonal gauge and the radial (gain) mode.
  - V53 below is the requirement that the vertex weights be computed in gauge-invariant variables.
- **Strong minimality (Bray-Seguin).** The token dynamics has no invariant left ideals: no token subspace is exactly
  preserved.
  - Theorem A(4): finitely many layers' transports of any positive token state dominate the identity.
  - Hence no static sparse support exists, which is why the pattern windows of note XII failed. The only admissible
    reductions are dynamical: transported bases, as in the old tier.
- **Positive extension.** The fourth-cumulant closures as local data with no global tensor, the flex kernel of note
  XXXVIII.

## 4. What the tower says about the bill

The bill is a sum over (token level, age) of the price of exactness. Measured on the adopted system (note XXXIX
section 11g):
- the young tier, about 125 units, pays for the recent births' full-irrep content;
- the old tier, about 60 units (formings about 20, factor-space hubs about 16, joins about 23), pays for old content;
- the covariance, 6.6 units.

**The Young decomposition prices the slice by irrep.**
- An old source's additive (trivial + standard) (2,1) readout needs only its row and column sums.
  - Column sums: (1^T (A o P)) A^T, with 1^T (A o P) = sum_p FA_p o FP_p because Q^T Q = I. This is O(n r).
  - Row sums: diag(Q G Q^T) with G = FA d(w2 s) FP^T. This is O(n r^2).
  - Together that is about 0.05 units per old source-layer, against about 0.6 for the dense forming plus the
    factor-space hub.
- The flat part needs the full Hadamard readout.

**The truth's D21 is 80% additive at depth.** At layers 10-14, S1 + A1 holds 0.80 of the energy and S2 + A2 0.19
(`../leaderboard-system/outputs/d21audit_off0.txt`). The chain's error is all in the flat part.

**Hypothesis H (the irrep split of the old tier).** After five or more transports, a source's flat (higher-irrep)
(2,1) readout is worth little at the output. Strong mixing randomises its flat content, while its additive content,
which couples to the collective modes, keeps acting coherently. If H holds, old content is nearly free in its
additive form, which is the signature of the leaders that note XXXIX's public post could not explain: their bill
equals ours minus the old tier while old content is carried.

## 5. Predictions (cold harness, networks 0-15, paired against the adopted fold system V35 + `V52_FB_FOLD=1` in the same batch)

| switch | what it removes or adds | prediction |
|---|---|---|
| `V54_OLD_D21=2` | the old tier's flat (2,1) readout (additive part kept exactly) | **P16 (H):** raw +0% to +4% |
| `V54_OLD_D21=0` | the old tier's whole (2,1) readout | P17: raw +30% or more |
| `V54_YNG_D21=2` | the young hub's flat (2,1) readout | P18: raw +20% or more (where the fold's value lives) |
| `V54_OLD_D3=0` | the old sources' diagonal readout | P19: raw +5% or more |
| `V53_TADPOLE=1` | tadpole-dressed vertices for legs, births and replicated slices | P20: raw -0.3% to -2%, FLOPs unchanged |

**Decision.**
- If P16 holds, the next step is the collective old tier: old sources read out through the additive part only, in
  factor space, so the factor-space hub goes.
- If P19 fails as well (old D3 small), the dense forming goes too.
- That tier is then priced with the profiler and measured in adjusted MSE as in note XXXIX.
