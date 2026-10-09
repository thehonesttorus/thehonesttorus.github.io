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

## 6. P16-P20, measured (`outputs/irrep_cold_16nets/`)

Cold harness, networks 0-15, every variant paired against V35 + `V52_FB_FOLD=1` rerun in the same batch (mean raw
2.216e-8, 0.1954 B):

| switch | raw against the fold system | better on | FLOPs | prediction |
|---|---|---|---|---|
| `V54_OLD_D21=2` (old tier: additive part only) | +383% +- 18 | 0/16 | +0.01% | P16 (H), 0 to +4%: **fails** |
| `V54_OLD_D21=0` (no old (2,1) readout) | +847% +- 45 | 0/16 | -0.91% | P17, >= +30%: holds |
| `V54_YNG_D21=2` (young hub: additive part only) | +2027% +- 81 | 0/16 | +0.01% | P18, >= +20%: holds |
| `V54_OLD_D3=0` | (invalid run) | | | P19: repeated in section 7 |
| `V53_TADPOLE=1` | +22.5% +- 14.2 | 3/16 | 0.00% | P20, -0.3% to -2%: **fails** |

**H fails, and by a wide margin.**
- The additive part of the old readout recovers about half of what dropping the readout costs (+383% against +847%).
- The flat (higher-irrep) content is where an old source's value lies, as it is for a young one (P18).
- The projection was unweighted. It also counted the saturated rows that production masks after the readout (their
  pre-activation slice is not small), so the test is conservative. A gate-weighted projection could recover more, but
  not the 380 points between H and the measurement.
- Even in a weighted metric, the factor-space additive readout costs about what the hub it would replace costs:
  O(n r^2) per term type per source-layer, against the hub's O(n^2 r) per slab, with r = 320. So the system
  conclusion does not depend on the weighting.

**What this says about the tower.**
- The S_n-isotypic split is the annealed decomposition. The network's quenched weights break the symmetry, and the
  value sits in the quenched (flat) part.
- This matches strong minimality (section 3): no static reduction exists, only dynamical ones. The old tier's
  transported basis is such a reduction; the irrep split is not.
- Task #48 (the collective old tier) is withdrawn.

**The tadpole.**
- Neutral to slightly adverse on 13 networks (-3% to +7%), and unstable on three (+85%, +219%, +35%).
- The dressed w3 and w12 at large |alpha| are Edgeworth tails (He_4(3) = 30): there the truncated series is not a
  reliable vertex weight.
- Section 7 separates the kappa3 tadpole from the kappa4 one, whose input is the closure's diagonal.

**The invalid run.** The `V54_OLD_D3` switch read the young sources' D3 after their left factor LP had been completed.
That added the P*P*s and 2 M*P terms to the young part, so the +6808% measures a corrupted young D3, not the old one.
The switch now reads the young part where production does, and the run is repeated in section 7.

## 7. The level-4 star (stated before its runs)

**The term.** The c(3) vertex with three covariance arms creates, at every layer, the fourth-cumulant bulk

    sum_m c3_m Sym(e_m x a_m x a_m x a_m),   a_m = d(Phi) C_off e_m,   c3 = E f'''(z),

the K_{1,3} tree of the Wick expansion. It shares two things with the kappa3 star sum_m c2_m Sym(e_m x a_m x a_m):
- the same arms;
- the same localized centre.

Its external lines take the same first-order gate, c(1) = Phi per line, as the kappa3 legs do. So, by Theorem 1
(i)-(ii) one level up, the level-4 star rides the kappa3 source's A and P legs at no transport cost. Its diagonal at
every later layer is

    kappa4(z_i)  +=  4 sum_sources sum_m c3_m A_im^3 P_im,

O(n^2) per source-layer, about 0.3% of the bill (`V55_K4STAR=1`).

**Why it is new.**
- The chain's fourth-cumulant state (note XXXII) declares only the diagonal and the (2,2) slice at the y level.
- The (1,1,1,1) class, where the star lives at birth, is absent; only the closure's mixture terms stand in for it.
- The star is quenched: it carries the network's own arms, the kind of content section 6 found valuable.

**What is left out.**
- The star's (2,2) and (3,1) slices: one hub-like product per source-layer.
- The path P4: two c(2) vertices joined by a covariance line, a V2 move. Its diagonal is diag(X C X^T) with
  X = (A o P) d(c2), one n^3 product per source-layer.
- The kappa3 x C class: likewise one n^3 product per source-layer.

**How it enters.**
- The star's diagonal joins the kappa4 diagonal after the closure's mixture gains and its (2,2) block have read the
  closure's own diagonal.
- From there it reaches the mean, the covariance and the y-level (2,2) slice through the exact tadpole terms of the
  Wick tables (for example ('g4row', 'c_off', (1,2), (5,2), 1/12) in the (2,2) program).
- The y-level diagonal's onward transport is O(1/n) (the harmonic projection's c_A = 6/(n+4)), so the star is not
  counted twice through the closure's state.

**Predictions** (cold harness, networks 0-15, paired against the fold system rerun in the same batch):

| switch | prediction |
|---|---|
| `V55_K4STAR=1` | **P21:** raw -1% to -6%, better on at least 11/16; FLOPs +0.2% to +0.6% |
| `V55_K4STAR_LOG=1` (same runs) | **P21a:** the star's rms diagonal is at least 5% of the closure's at layers 8-15 |
| `V54_OLD_D3=0` (corrected) | **P19':** raw +50% or more |
| `V53_TADPOLE=2` (kappa3 tadpole only) | **P20':** no network worse than +10%, mean within +-2% |

If P21 holds, the star goes to the scored regime: all 100 networks, counterterms refitted on 0-49, held-out 50-99,
adjusted MSE, as in note XXXIX section 11.
