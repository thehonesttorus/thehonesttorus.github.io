# Theory judge: what the fresh-slate principles contributed, and what they re-derived

*Judge note, 1 Oct 2026, written against the state of the streams at ≈ 20:30 UTC. It scores the six design streams and the three breakthrough drafts on theoretical depth and remaining promise. Read for it: BRIEF.md (with the 18:55 correction: the Gaussian closure scores raw ≈ 4.1–4.3e-6 at n = 1024), CONVERGENCE.md, every `designs/*/DESIGN.md`, `RESULTS.md` and `VERDICT.md`, the three `breakthrough/*/REPORT.md` drafts (and their result files where cited), `foundations-problem.md`, `foundations-unlocks.md`, `foundations-input-verified.md`, `foundations-facts.md`, `scaffold/KIT.md` and `scaffold/RESIDUAL.md`. The public moment chain appears only as a baseline, as BRIEF §2 rule 1 allows. Its numbers come from `notes/digests/504aldo-k3-chain.md`, cited by finding number (F..). Nothing here proposes adapting it.*

*Labels. **Exact**: an identity proved in the cited note. **Measured**: a number from the cited run; the network set is stated when it is not the six-network bench `w1024_d16`. **Derived here**: computed for this note; the check is in Appendix A. **Judgement**: the judge's reading of the evidence.*

---

## 0. Verdict

1. **Every fresh estimator that propagates non-Gaussian structure computes one object: the first-order layer of the moment chain.** That covers five of the six designs and the region and costate drafts; tropical stayed at the Gaussian closure, and interpolation corrects the closure instead. The object has four parts:
   - a Gaussian (or copula) reference;
   - third-cumulant content born at every ReLU as a hub ("star") diagram, whose vertex coefficient is E relu″ = φ(t)/σ at the birth neuron, with two legs carrying the transported cross-covariance C_s diag(Φ_s) P_{s→l} and one leg carrying the mean-gate propagator P_{s→l} = W_{s+1} Φ_{s+1} ⋯ W_l;
   - injection of the (2,1) slice into Cov(a) by first-order bivariate Edgeworth;
   - a per-neuron Edgeworth readout.

   Faces calls these "facet births on renormalised legs", heisenberg "star diagrams pulled back to the source layer", markov "unary sites on a latent field", region "second-chaos births" and costate "source atoms". Term by term this is the young-source hub of the public K3 chain; in faces' case the coefficients are identical (§1.1). Bethe edge0 and signings v1 are the same chain with history collapsed into the carried (2,1) slice, which is the public repo's memoryless slice closure.
2. **The designs also land where the chain's own ablations land** (§1.2). For example, bethe edge0 measures 1.04e-6 at 0.10 B against F65 keep = 0 at 9.08e-7 and 0.09 B; faces' window 6 measures 4.45e-7 against F61 W = 6 at 4.35e-7. The full-history designs measure 3.2–4.2e-7 against 4.47e-7 for the chain with its residual leg removed (F63). On the adjusted metric every rung of the chain's documented negative frontier beats every fresh design: F65 runs from 9.7e-8 at worst to 2.2e-8 at best, at the chain's pre-Strassen prices, against the fresh best of 1.07e-7.
3. **One approximation, shared by all of them, explains the shared ceiling** (§2, derived here). They all transport old content with mean gates:
   - On the all-distinct pattern this is exact at leading order.
   - On coincident patterns it is wrong at O(1). The exact first-order response is ∂κ(a_i,a_i,a_j)/∂κ(z_i,z_i,z_j) = Φ_j[Φ_i − E(a_i)φ_i/σ_i]; mean gates give Φ_i²Φ_j instead, which is −27 % at t = 0 and +22 % at t = 1. For the diagonal pattern the exact response at t = 0 is 0.341; mean gates give 0.125.
   - Coincident patterns are exactly what the next layer's paired fresh weights contract coherently, so this error does not average out.
   - This is heisenberg's 22–31 % slice error "from the first transport step" (R10), region's 14–25 % D21 error, and the dilation-mode leak of costate C6.
   - Through the ε² law (≈ 4e-6 ε²) it accounts for where every full-history design landed: 3–4e-7.
   - The public chain's per-layer residual leg contains the repair: it re-projects the exact slice response at every layer and births the difference as a new source. The chain's own ablation without that leg is 4.47e-7; with it, 3.9e-8 (F63).
4. **This convergence is a finding about the problem, not a failure of six streams.** Three properties of the problem make this chain its unique first-order perturbation theory:
   - Fresh independent Gaussian weights: Isserlis pairings, loops of fresh-weight legs suppressed, no Onsager memory.
   - n^{-1/2} is the only small parameter: there is no zero-temperature regime and no spectral gap in the signed bulk.
   - The ReLU's derivative profile is a sequence of face measures: the gate, the wall density, and its slope.

   Every principle that was carried to quantitative accuracy reproduced the chain. This campaign's own pre-fresh-slate measurements had already catalogued the complete local law (foundations-facts F6.3: first-order gate diagrams B0–B7 with counting coefficients reproduce the next D21 to ≈ 1 % at n = 1024; one network, preliminary). The fresh designs rebuilt its Wick part and its pass-through term B0, then named the missing members as their "next step": the curvature passage, the joint gate law and the (2,1,1) κ4 slice. The ruling can exclude the chain's code; it cannot exclude this mathematics. The honest test of a fresh design is therefore what its principle adds that the chain lacks.
5. **What the principles added is real, and most of it is framework-grade** (§4). It includes:
   - the Duhamel telescoping and its bilinear remainder (heisenberg);
   - E f = E Δf, facets rather than faces as the carriers of non-Gaussianity, exact telescoping along face-averaged arrows, and the uniform age profile (faces);
   - the double-edge obstruction and the conserved dilation charge (costate);
   - two-site κ4 trees, and the curvature passage obtained by conditioning on each site's own latent variable (markov);
   - loop negligibility of fresh-weight legs and the parity blindness of signings (signings);
   - the global order parameter (bethe);
   - the 10²³ tropical cancellation and the high-temperature phase diagram (tropical);
   - the no-Onsager theorem and the trace channel (interpolation);
   - an error calculus whose books close at n = 1024 (region).
6. **One structure was found that the chain does not use as a derived object: the dilation/Perron sector.**
   - Positive homogeneity makes the scale-mixing law an exactly conserved charge of every arrow (costate C6).
   - Its signatures are measured in five places:
     - the rank-one spike of the old (2,1) slice;
     - the coherent node κ4 (bethe; the thin-shell scalar of foundations);
     - λ3 ∝ t with R² up to 0.89 at n = 1024;
     - Var τ_l growing linearly in depth;
     - one depth-growing mode holding 64 % of the energy of the deviations of the public chain's fitted λ-table from its reference values (F75).
   - Carrying it costs O(n²) per layer.
   - The chain contains it only implicitly or as fitted constants.
   - It is the one open lead for carrying old content cheaply (§5, D1). It is unproven at n = 1024: costate's law-level variant uses half the products of exact first order for ≈ 10 % more error, on one network.
7. **The open question is now narrow.**
   - The local law is known.
   - The binding object is memory: the ≈ 40 % of each D21 that is old, orthogonal to the present and of rank ∝ n.
   - Costate C3 shows that exact first-order evaluation of memory costs Θ(L² n³) in every natural contraction order. At the per-pair constants measured so far (3–7 products for each of ≈ 120 (source, target) pairs) that is ≈ 200–480 units even at Strassen prices.
   - The leaders' raw 1.5e-8 at 0.11 B (≈ 7 units per layer) shows that a cheaper representation reaching the bar exists. Either it exploits structure in the memory, or it brings the per-pair constant to about one product. The public chain's author judged it is not a compressed version of his representation.
   - The dilation/Perron sector is the only structural candidate this round produced.

**Scores** (0–10; theoretical depth and remaining promise, judged together):

| rank | stream | score | depth | promise | estimator verdict |
|---|---|---|---|---|---|
| 1 | costate (breakthrough draft) | **8.0** | structure theorems C1–C6 of the converged object; one conserved charge | highest: the dilation sector at law level, the readout-metric content of memory | first order = the chain; the law-level scale charge is new |
| 2 | region (breakthrough draft) | **7.0** | fresh-weight Frobenius lemma; local errors priced by K(l) close the books to 6 % at 1024; necessary conditions | high, as the evaluator and target-setter for any next design | FC = the chain in the chaos basis |
| 3 | heisenberg | **6.5** | Duhamel identity and bilinear remainder, exact and general; the cleanest localisation of the shared defect (R10) | medium: the identity can host a dilation-invariant reference | collapsed: the adjoint form of the chain (its own words) |
| 4 | faces | **6.0** | E f = E Δf; facets, not faces; exact telescoping; age attribution; the chaos ceiling | low beyond the identities | collapsed: the chain's hub with identical coefficients |
| 5 | bethe | **5.5** | dense power counting; exact edge table; the global order parameter | medium (v3/v4 = the dilation sector) | collapsed: the memoryless slice closure, same number |
| 6 | markov | **5.5** | exact single-site formulas; two-site κ4 trees; the curvature passage from a per-site separator; the CMI price list | low–medium (MKV-2 is O(L² n³)) | collapsed: the chain's sources and windows |
| 7 | interpolation (breakthrough draft) | **5.5** | no Onsager memory; chaos-1 exactness of the closure; covers closed; self-averaging split; trace channel | medium–low | not a chain: corrections to the closure |
| 8 | signings | **4.5** | loop suppression of fresh-weight legs; the Z/2 parity theorem; cover limit | low (closed by its own analysis) | collapsed: the memoryless slice closure at ≈ 2× its cost |
| 9 | tropical | **4.0** | I2, I4; the 10²³ cancellation; the phase diagram; the uncentred split | low (closed) | closed: TCT-0 = the Gaussian closure |

---

## 1. The common object

### 1.1 Term map

Public-chain names are V16/V25 (digest §1.3). In the public chain, A is the transported birth leg a_b = diag(w1) C_off and P the propagator, born as W and updated P ← W diag(w1) P. Its weights are w1 = Φ(t) and w2 = φ(t)/σ = E relu″, and e = 2 Cov(a, gate).

| component | public chain | faces 'w16' | heisenberg `hdpull` | markov MKV | region FC; costate | bethe edge0; signings v1 |
|---|---|---|---|---|---|---|
| reference | covariance propagation + term program | Mehler covariance (R = 6) | Gaussian ν_l with exact (m, C) | latent Gaussian field, Cov_φ | exact bivariate Gaussian closure | exact edge table; copula (T, R) by Mehler inversion |
| birth vertex | w2 | c = φ/(2s) = w2/2 | α₂ (Hermite) | E X̃″ = E φ″ | w2 | L₋₁ (density at 0) |
| cross-covariance leg | A | K = C_s diag(Φ_s) P_{s→l} | R = ρ diag(α₁) U | K = C_m D_m P | Y = S diag(Φ) Z | c_ab L₀ (one step) |
| propagator leg | P | P_{s→l} = W D_β ⋯ W, β = Φ | U ← U Φ W | P^{(m→l)} | Z_s(l) | none (one step) |
| hub contraction | D3 = 3 Σ w2 A²P; D21 ∋ 2(A∘P∘w2)Aᵀ + (A∘A∘w2)Pᵀ | κ3 += 6 Σ P c K²; D21 += 4 MᵀK + 2 (K∘K)ᵀ(P∘c): **same coefficients** | S = 2((α₂U∘R)ᵀR) + (α₂R∘R)ᵀU | 3 p E[X̃″] k²: same coefficient | Σ w2 [Y²Z + 2 YZY] | T^a_abc = Σ c_ab c_bc L₀L₋₁L₀ |
| coincident and diagonal births | residual leg: S3c, S_sep = diag(e) C_off diag(w1) | inside K's diagonal | exact coincident planes D2, diagonal D1 | 3 p² E[(X̃²)′] k, p³ κ3(X) | "exact Gaussian slices of each birth" | carried (2,1) slice |
| exact slice response of transported content (re-birth) | residual leg Rres (rank 16), D21 feedback | absent (κ3-only curvature term: ≤ 3 %) | absent (named, R10) | MKV-2 (1.4× at width 256) | FC: "first-order non-Gaussian slices" (D21 still 14–25 % off) | the slice is re-projected each layer, memoryless |
| (2,1) slice into Cov(a) | term program | ½ (φ/s) Φ D21 + transpose | Stein injection | ½ (E φ″ E φ′ D21 + transpose) | ½ (E[δ(z_a) 1(z_b > 0)] D21 + transpose) | hyperedge terms |
| κ4 | regenerated core diag(dG) + λ C_off, fitted λ table | two-birth tree 48 Mᵀ C M | none (oracle only) | two-site trees | none; costate: law-level scale charge | bethe v3/v4: spike edge belief |
| history | all ages, tiers at rank 384 / 224 | window or all ages | ages ≤ A, or all | window or all | all pairs; A + slice; A + gl | none / age 1 |

Two qualifications:
- Faces' legs keep diag(C_s), which the public chain books in its residual leg.
- Markov splits the Gaussian birth terms between the field and the site, so its coincident coefficient is e − 2Φσφ rather than e.

Neither changes the object.

### 1.2 Same configuration, same numbers

The fresh designs were measured on bench `w1024_d16` (6 networks unless stated). The public chain was measured on its P2-mini networks (4–8). The two sets agree at the Gaussian level: 4.10–4.30e-6 on the bench against 4.39e-6 for covariance propagation on the mini networks (F43).

| configuration | fresh design, raw at 1024 | public chain, raw at 1024 |
|---|---|---|
| (2,1) slice + fresh hubs, history collapsed | bethe edge0 1.04e-6 at 0.10 B; signings v1 1.04e-6 at 0.19–0.25 B | F65 keep = 0: 9.08e-7 at 0.09 B |
| + one age exact | bethe edge1 8.1e-7 at 0.19 B | F65 keep = 1: 6.72e-7 at 0.145 B |
| ages ≤ 3 exact, older into the slice | costate A3slice 6.87e-7 (3 nets; 435 products) | F65 keep = 3: 2.89e-7 at 0.24 B |
| window ≈ 4, older dropped or diagonalised | markov window 4: 1.26e-6; costate A = 3: 1.05e-6 | F61 W = 4: 1.02e-6 |
| window 6 | faces win6 4.45e-7 | F61 W = 6: 4.35e-7 |
| all ages; Gaussian births, mean-gate transport, no residual leg | faces 3.24e-7; costate 4.03e-7 (3 nets); markov 4.18e-7; heisenberg 4.23e-7; region FC 3.51e-7 (1 net) | F63 "S_sep alone": 4.47e-7 |
| + per-layer residual leg (rank 16 suffices) | — | F63: 3.88e-8; K3 reference 4.32e-8 (F43) |
| adjoint (pull-back) evaluation order | heisenberg, costate | F88, dead end #20: no cost drop |
| Gaussian scale mixture | costate gl (law level, retired old content only) | F78: memoryless scale-mixture closure 2.33e-6, dead end #4 (a different use) |

---

## 2. The shared approximation, and why it sets the ceiling

**The responses of one ReLU to a cumulant of its input** (derived here, Appendix A). Take a pair or triple of pre-activations at zero mutual correlation, standardised so that t = μ/σ. The exact first-order response of the post-activation cumulants is:

| input pattern | exact response | mean-gate rule (all designs) |
|---|---|---|
| all distinct, κ(z_i,z_j,z_k) | Φ_iΦ_jΦ_k | Φ_iΦ_jΦ_k (exact) |
| (2,1), κ(z_i,z_i,z_j) | Φ_j[Φ_i − E(a_i)φ_i/σ_i] | Φ_i²Φ_j |
| (3), κ3(z_i) | Φ + ½ tφ ẽ₂ − ẽ₁φ − ẽ₁² tφ, with ẽ₁ = E a/σ, ẽ₂ = E a²/σ² | Φ³ |

Numbers (Φ_j factored out of the (2,1) row):

| t_i | −1 | 0 | +1 | +2 |
|---|---|---|---|---|
| (2,1): exact / mean gate | 0.139 / 0.025 | 0.341 / 0.250 (−27 %) | 0.579 / 0.708 (+22 %) | 0.869 / 0.955 (+10 %) |
| (3): exact / mean gate | 0.131 / 0.004 | 0.341 / 0.125 | 0.528 / 0.596 | — |

At non-zero correlation there are O(ρ) corrections. foundations-facts F6.8 resums them exactly through the pair gate covariances (T-class and leaf identities).

**Why this error survives averaging.** The coincident patterns of κ3(a_l) are the ones that the next layer's fresh weights contract in pairs (Isserlis), and pairs are the coherent, unaveraged channel of the error table in foundations-problem §3.3. Interpolation's T5 (result files, one network) measures the size of that channel at n = 1024. The chaos-1 part of κ3(z_{l+1,p}) in the fresh column w_p, namely 3σ²(t_l·w_p) with t_l = E[‖ã_l‖² ã_l], explains:
- 65 % of the per-neuron κ3 at layer 2;
- 93 % at layers 8–16.

The mean-gate error on this channel has a fixed sign for hot units: −27 % at t = 0. It is the same defect four streams localised:
- heisenberg R10: the transported κ3 slices are 22–31 % wrong from the first transport step, while the sources are exact;
- region FC: D21 is 14–25 % off from layer 2 on;
- costate C6: Euler's relu′(z) z = relu(z) makes the scale direction an eigenvector of eigenvalue 1, and mean gates replace E a by Φ m, so they leak it;
- markov: its curvature passage (the ½ψ″F² term of MKV-2) is the derived fix, and faces' "later facets fold old skew away" describes the same mechanism.

**Arithmetic.** The final MSE grows with the relative error ε of the (2,1) slice at every layer as 3.8e-6 ε² (region) or 4.2e-6 ε² (public F71; foundations-facts F10.1). At ε = 0.15–0.25 the covariance channel alone costs 0.9–2.6e-7. The diagonal κ3 used by the readout carries the same 20–30 % error. Region prices that channel at ≈ 1.5e-6 for ε = 1, so it adds ≈ 1e-7 (estimate). The total, 3–4e-7, is where all five full-history designs measured.

**Two halves of one object.** The memoryless designs, bethe edge0 and signings v1, apply the exact coincident response: their pair maps expand E[relu²(z_a) relu(z_b)] to first order in the carried slice, and the result is Φ_b[Φ_a − E(a_a)φ_a/σ_a], not Φ_a²Φ_b. But they carry no memory. The full-history designs carry memory but transport it with mean gates. Each has half of the chain's K3 layer, and each lands on the chain's corresponding ablation: ≈ 1e-6 (F65 keep = 0) and ≈ 4e-7 (F63 without the residual leg). No stream has yet measured both halves together at n = 1024; markov's MKV-2, still running, is the one attempt.

**The repair.**
- At the slice level the repair is elementwise: n² work per layer.
- For transported content the coincidence happens at intermediate layers, inside the Hadamard square of a propagator (costate's double edge). There the repair is a re-birth at every layer, which is exactly the chain's residual leg.
- The bar needs ε ≲ 4.5–6.3 % (region N2). For a 10 % margin over raw 2.1e-8, foundations-facts F10.1 needs ε ≤ 2.2 %.
- Region's draft shows that a re-birth is necessary but not automatically sufficient: FC re-births first-order slices and still sits at 3.5e-7 on one network. Locating why is its open item.

---

## 3. Design by design

Each entry gives (a) whether the principle drove the estimator or collapsed into the chain, with the coinciding terms named; (b) its exact identities and structural facts; (c) what remains open from that principle.

### 3.1 costate (breakthrough draft): 8.0

**(a)** The Heisenberg principle, taken to its end, drove the stream essentially. Its first-order estimator ("all pairs", 4.03e-7 on 3 networks, 855 products) is the heisenberg estimator, hence the chain. Its contribution is a structure theory of that object and one new sector.

**(b)**
- **C1.** The first-order co-state is n-dimensional and exactly closed by the backward recursion G_l = T_l* G_{l+1} + J_{l+1}.
- **C2–C3 (the double edge).** (UM)∘(VM) = (U⊙V)(M⊙M). Consequences:
  - no exact n²-per-layer forward state exists for old content;
  - the (2,1) slice chain is the double edge averaged over the newest weights: it keeps the paired a = b part and drops the unpaired part, which is of the same order;
  - first order costs Θ(L² n³) in every contraction order considered.
- **C4.** A rank-d propagator has a rank-d² double edge, so Tucker truncation saves nothing below d ≈ √(n/2) ≈ 22.
- **C5.** Second order is Θ(L³ n³).
- **C6.** The dilation charge is exactly conserved: F#(P ⋆ ν) = P ⋆ F#ν.
- C3 and C4 explain from first principles a cluster of the public repo's dead ends, which were found by trial: #6, #18, #20, #21 and #25.
- **The oracle filter at width 128.** The readout reads old content mainly through one spike plus a few modes. Exact diagonal plus a rank-1 off-diagonal recovers 84 % of the old-content gain (from 1.24e-4 to 4.65e-5, against 3.13e-5 with everything).
- **Caveat (Judgement).** C3 is a cut lemma (generic position) plus an analysis of the natural contraction orders. It is not an algorithmic lower bound: randomised or approximate evaluations are outside it.

**(c)**
- The law-level scale charge ("gl"/"gsl") needs its six-network run at 1024. So far, on MLP 0:
  - A3gsl measures 4.35e-7 at 435 products, against 3.9e-7 for all pairs at 855: half the products for ≈ 10 % more error, ≈ 1e-7 adjusted at Strassen prices;
  - at width 128, gl beats all pairs at a quarter to half of the products.
- The principle-level version of the chain's residual leg is C1's recursion with the exact coincident responses of §2 in the linear-response arrow, plus C6 for the scale direction. Whether that can be evaluated below Θ(L² n³) is the open question of the round.

### 3.2 region (breakthrough draft): 7.0

**(a)** The principle, feasibility by necessary conditions, drove the analysis essentially. Its estimator FC is the chain in the chaos basis (Y = S diag(Φ) Z; D21 = Σ w2 [Y²Z + 2 YZY]), with a slice re-birth. It measures 3.51e-7 on one network.

**(b)**
- **The fresh-weight lemma.** Errors enter the next layer as a Gaussian chaos in the fresh columns, priced by rotation invariants, so Frobenius norms are the currency: K_coh ≈ K_off to ±30 %.
- **Transfer coefficients K(l)** measured at 1024.
- **The books close.** The priced local errors sum to 5.29e-6 against a measured 5.00e-6.
- **Anatomy of the Gaussian closure's error:** 71 % non-Gaussian covariance correction, 28 % mean readout, 0.3 % variance.
- **Local sufficiency.**
  - The D21 term carries 92–100 % of ‖ΔC‖².
  - Adding the κ4 (2,2) term leaves 4.4e-9.
  - The (3,1) slice is unnecessary (3e-10).
  - Second-order Edgeworth on the covariance is not needed at 1024.
- **Necessary conditions N1–N5**, including the derivation of the ε² law.

This is the Lipschitz-of-questions calculus (unlock 39) made quantitative, and it corrects two width-64 diagnoses (§6).

**(c)** It is the right evaluator for any next design, at hours rather than days. Its FC defect search should find the mean-gate coincident response of §2. Its pricing of old content by age at 1024 is item 2 of foundations-facts §16 and decides D1. One network so far.

### 3.3 heisenberg: 6.5

**(a)** The Duhamel principle drove the error accounting essentially. The estimator collapsed, as the stream itself says ("coincides with the adjoint form of a second-order cumulant chain"):
- `hdpull.eval_DS` is the chain's hub contraction plus the Gaussian birth part of its residual leg (coincident planes D2, diagonal D1);
- transport is by mean gates.

The pull-back is not a cost lever:
- the public chain already contracts factored legs per (source, target) pair without n³ tensors;
- its adjoint ordering is public dead end #20 (F88);
- costate C3 shows that no natural contraction order changes the cost order.

**(b)**
- **Theorem 1** (Exact): any reference chain differs from the truth by Σ_l (ρ_l − ν_l)[g_l], with the exact pulled-back questions g_l. It is the commutative form of Chen–Rouzé telescoping (unlock 36).
- **Theorem 2** (Exact): the remainder is bilinear, local defect × question-model error.
- R0 checks Theorem 1 to Monte Carlo precision.
- R10 is the round's cleanest localisation of the shared defect: sources exact, transport 22–31 % wrong from step one.
- R11 is an honest correction at 1024: the earlier 30× claim was against a mis-calibrated reference, and the gain is ≈ 10×.

**(c)**
- Its step 1, exact κ3 transport through one ReLU, is §2's response table plus the F6.3/F6.8 diagram list: cheap at the slice level, a re-birth for transported content.
- Its step 2, old joint κ4, should be narrowed by the n = 1024 facts. The (2,2) slice is needed in the covariance (region). The (2,1,1) slice is needed in D21 generation (foundations-facts F7.1: 2.9–3.2 % without it, 0.8–1.0 % with it). The (3,1) slice is not needed.
- The identity's best remaining use is with a dilation-invariant reference family, P ⋆ N(m, C) (costate C6), so that the local defects carry no scale mode.

### 3.4 faces: 6.0

**(a)** The principle drove the diagnosis essentially and the estimator not at all:
- 'w16' is the chain's hub with identical coefficients: 6c = 3 w2, 2c = w2, 4c = 2 w2, with K the chain's A leg (§1.1);
- the stream itself calls the bookkeeping standard;
- win6 lands on the chain's W = 6 window.

**(b)**
- **E2** (Exact): E f = E Δf. The mean is the Gaussian mass of the walls; barycentres are facet integrals. Three streams derived it independently: faces, tropical I2 and foundations-problem 2.1(iv).
- **T1** (Measured): codimension-0 face data carry ≈ 1 % of the skew at layer 2 (rms 0.003 against 0.266). κ(g_i,g_i,g_k) = (1 − 2p_i)Γ_ik vanishes at p = ½. So the face algebra D₀ is the wrong carrier; non-Gaussianity lives on facets. This is a programme-level negative for the "objects are barycentres" dictionary at codimension 0.
- **E5** (Exact): telescoping z̃_l = x P_{x→l} + Σ_s ν_s P_{s→l} along face-averaged arrows. It is the same decomposition as markov's (★) and region's chaos births.
- **T3** (Exact attribution, Measured): the final skew is born uniformly over depth, and the newest facets contribute ≈ 0 net as the first leg.
- **The chaos ceiling** (Measured): chaos 1 and chaos 2 in x each hold 18–23 % of the variance at depth, and ≈ 60 % sits in chaos ≥ 3, at widths 64–512. Input-coordinate expansions are therefore capped, and per-layer renormalisation is forced.
- **T4** (Measured): decoupling gates from face-wise gradients is 5–10× worse than closing the pair law.

**(c)** The proposed facet-conditional covariance Cov(z_l | z_{s,k} = 0) is the curvature passage of §2. Its deciding oracle experiment is largely answered by region §3. The promise lies in the identities, not in the estimator.

### 3.5 bethe: 5.5

**(a)** Power counting drove the design essentially. Edge0 is the chain's memoryless (2,1)-slice closure plus the newborn hub (L₋₁ = w2, L₀ = Φ). It reproduces that closure's number: 1.04e-6 at 0.10 B against F65 keep = 0, 9.08e-7 at 0.09 B; edge1 against keep = 1 gives 8.1e-7 against 6.72e-7. Its best-in-round adjusted score (1.07e-7) is the memoryless closure's (9.1e-8).

**(b)**
- **The lift-limit estimate fails.** On the dense layered graph, node beliefs alone are wrong at the leading quenched order and do not improve with width (2.5e-3 at width 64, 3.3e-3 at 128).
- **Edge beliefs and generated triples** enter at the same order (with signings; unlock 13).
- **An exact bivariate edge table**, built by a Stein recursion, accurate to 1e-16.
- **The global order parameter.** Cavity fields are Gaussian conditionally on q_l(x) = |a_l(x)|²/n. The rank-one spike in κ(z_a,z_a,z_b,z_b) is seeded by the input radius (a third of it at layer 1) and amplified by depth: ε rises from 0.06 to 0.74, against 2/n = 0.03 for the radius alone.

**(c)** v3/v4, carrying Q as a scale mixture, is D1 below; costate C6 makes it exact. Caveat (§6): the width-64 reading "node beliefs bind, pairs do not" fails at 1024.

### 3.6 markov: 5.5

**(a)** The Hammersley–Clifford ansatz drove the state, and the state is the chain's sources:
- sites = births, P^{(m→l)} = the propagator, K^{(m→l)} = C_m D_m P = the A leg;
- the single-site κ3 terms are the hub (3 p E[X̃″] k², same coefficient), the coincident term and the diagonal;
- the window is the chain's windowing (window 4 gives 1.26e-6 against W = 4 at 1.02e-6);
- retiring old sources to a per-neuron projection is the memoryless collapse.

**(b)**
- **Exact single-site Stein formulas** (checked to 1e-16).
- **Two-site κ4 trees.** For κ4 they are leading order, because their weights p_c²p_d² are positive and add coherently (the layer-2 κ4 error falls from 26 % to 4 %). For κ3 they vanish, because E X′ = 0. So the minimal Markov network has pairwise non-Gaussian cliques at fourth order. Faces' 48 Mᵀ C M is the same tree.
- **The curvature passage** is Gaussian conditioning on each site's own latent variable: a junction tree with a per-site separator. It is the principled derivation of the coincident terms of §2 (layer-3 κ3 slope 0.842 → 0.957).
- **CMI as the price of separators** (P3), with the separator test: the top n/4 transported directions recover half the log gap.
- **The history axis.** Every variant lands within 2× in adjusted score, so the accuracy bought with history is paid back in cost.

**(c)** MKV-2 at 1024 was still running (it gained 1.4× over full-history MKV at width 256). If it closes the §2 gap, it will have re-derived the chain's residual leg from a junction-tree principle, at O(L² n³).

### 3.7 interpolation (breakthrough draft): 5.5

**(a)** Not a chain re-derivation. It builds a different kind of object, disorder-conditional corrections, on top of the Gaussian closure.

**(b)**
- **Theorem A**: the smart path in the weights is the OU semigroup, and by Stroock's formula the quenched content is Hermite chaos in W.
- **Chaos-1 exactness of the closure readout**: J₁(F − G) = 0, by the symmetry w ↦ −w.
- **B1**: there is no Onsager memory, because each W_l acts once, on a law independent of it. Old content is quenched transport, not AMP memory.
- **C1–C2 (covers)**: input-shared lifts are function-identical; input-lifted towers tend to the tree, which is the diagonal closure and 48× worse. With signings, this closes the Bethe-on-covers route.
- **The radial factor**: −18 % on the closure, on 6 of 6 networks.
- **The norm process**: Var τ_l grows linearly in depth, with layer-to-layer correlation 0.97. The mixture readout is not better, because the fluctuation is anisotropic.
- **The self-averaging split**: 58 % of the closure's error is a smooth function of (α_j, s_j) (leave-one-network-out).
- **T5** (result files): the trace channel of §2.

**(c)**
- (Judgement) The self-averaging part is plausibly the dilation/trace sector's imprint on the readout: the basis s φ(α) α^k spans the Edgeworth readout of λ3 ∝ t with λ4 constant. If so it is derivable, and the fit is unnecessary.
- It has not been reconciled with public F47, where covariance propagation plus an online correction gave 1.0×.
- T5 is a structural fact. As a cost lever it is small, since the exact contraction costs one product per layer, and its 7 % remainder exceeds the λ3 readout tolerance.

### 3.8 signings: 4.5

**(a)** The principle drove the negative structure essentially. v1 is the memoryless slice closure with a copula marginal (1.04e-6), the same as bethe edge0 and F65 keep = 0, at about twice their cost.

**(b)**
- **Loops of fresh-weight legs** are suppressed by n^{-1/2} per cycle (T1: the loop residual falls at least as n^{-1.5}). The step sums are trees, so determinants have nothing to resum.
- **The Z/2 parity theorem.** Signed averages annihilate odd content, and everything needed beyond the covariance is odd. This closes Godsil–Gutman-type estimators structurally.
- **Z/3 sketches** are unbiased, but need 1.6e3 samples at 16 sources and 3.5e8 at 1024.
- **The cover limit** is the independent-input model.
- **An exact Mehler/hafnian generator** of a copula's joint cumulants.
- **The bottleneck is drift.** The local step error falls as n^{-3.8}, then n^{-2.9}, while drift adds a further 10× at width 512.

**(c)** Closed by its own analysis. The one item to carry forward is its Isserlis sector table, which is where the coherent (paired-leg) channel sits.

### 3.9 tropical: 4.0

**(a)** The principle drove the stream essentially, to its own closure:
- TCT-0 is the Gaussian closure (4.10e-6);
- the temperature extrapolation is the one genuinely tropical estimator: at most 2.5× and only on optimistic selection at widths 64–128, untested at 1024;
- the sampled bridge gains nothing.

**(b)**
- **I2**: the mean is the multiplicity-weighted Gaussian mass of the tropical hypersurface; on the sphere, the total wall mass is (n − 1) times the spherical mean.
- **I4**: E softplus_T − E relu = (π²/6) T² p(0) + O(T⁴).
- **The P − Q cancellation**: 1.2e23 at layer 16, n = 1024, growing by √(4n/π) per layer.
- **The phase diagram**: about a quarter of the gates are hot at every depth, and the path sum is in weak disorder (β_eff ≈ 1.1 against √(2 ln n) ≈ 3.7).
- **The uncentred split is non-perturbative**: correlation −1.00, with E|∇z|²/s² = 17 at depth.
- **Temperature**: the closure error falls 400× by T ≈ 0.8, but none of that accuracy transfers back to T = 0.

**(c)** Closed. The only repair is a centred calculus. At the level of the identity that is unlock 37 (Malliavin–Stein), not tropical geometry, and its conjecture is untested.

---

## 4. Ledger: identities and structural facts worth keeping, regardless of the competition

| # | identity or fact | status | found by | meaning for the framework |
|---|---|---|---|---|
| 1 | Duhamel telescoping of any reference chain against the exact pulled-back questions; bilinear remainder | Exact | heisenberg Thms 1–2 | the error of any closure is a sum of local defects, each paired with an exact question; doubly robust |
| 2 | E f = E Δf; Δa_L = Σ δ(z)‖∇z‖² ∂a_L/∂a_l; barycentres are facet integrals | Exact | faces E2, tropical I2, foundations-problem 2.1 | the mean is carried by codimension-1 faces |
| 3 | codimension-0 face data carry ≈ 1 % of the skew | Measured | faces T1 | the face algebra D₀ is insufficient; arrows act nonlinearly on the boundary |
| 4 | telescoping along averaged arrows: z̃_l = x P_{x→l} + Σ ν_s P_{s→l} | Exact | faces E5 (≡ markov (★) ≡ region FC) | a composite of arrows is the composite of averaged arrows plus transported births |
| 5 | exact age attribution; the final skew is born uniformly over depth | Exact + Measured | faces T3 | no memoryless rule; every birth matters |
| 6 | chaos 1, 2 and ≥ 3 in x hold ≈ 20 %, 20 % and 60 % at depth, independent of width | Measured | faces | input-coordinate expansions are capped; renormalised legs are forced |
| 7 | loops of fresh-weight legs are suppressed; step sums are trees | Derived + Measured | signings T1, bethe §2, unlock 13 | tree gluing across width; loops live only in the state |
| 8 | Z/2 sign averages annihilate odd content; the scored information is sign-odd | Exact | signings, unlocks 5(d)/22 | determinantal and signing tools have no role in the step |
| 9 | covers: shared-input lifts are function-identical; input-lifted towers tend to the diagonal closure | Exact + Measured | interpolation C1–C2, signings | the lift end of a Bethe expansion is the wrong end here |
| 10 | no Onsager term along depth | Exact | interpolation B1 | memory is quenched transport |
| 11 | the closure's readout error starts at chaos 2 in the last weights | Exact | interpolation | quenched corrections organised by weight chaos |
| 12 | the trace channel 3σ²(t_l·w_p) carries 65–93 % of per-neuron κ3 at 1024 | Exact decomposition + Measured (1 net) | interpolation T5 | the coherent (paired) channel dominates at large n |
| 13 | thin-shell identity E_w κ4 = 3σ⁴[Var‖ã‖² − 2‖Σ‖_F²] | Exact (checked to 0.1–2.2 % at 1024) | foundations-problem 3.3 | the coherent κ4 is one scalar per layer |
| 14 | the dilation charge is conserved; the scale direction has eigenvalue 1; mean gates leak it | Exact | costate C6 (with the radial factorisation, foundations 2.1) | a symmetry with a sufficient invariant sector |
| 15 | the global order parameter: a rank-one κ4 spike, seeded by the radius and amplified by depth | Measured (width 64) | bethe; interpolation T1; public F75 | the coherent loop correction is a fluctuating scalar |
| 16 | two-site κ4 trees are leading order (coherent); there are none for κ3 | Derived + Measured | markov §7.3; faces | the minimal Markov network has pair cliques at fourth order |
| 17 | the curvature passage is conditioning on each site's own latent | Derived + Measured | markov MKV-2 | per-site separators |
| 18 | the coincident-pattern ReLU responses (§2) | Derived here | this note | the anatomy of the shared defect |
| 19 | the double edge: no n²-per-layer state; first order Θ(L² n³); Tucker rank d²; second order Θ(L³ n³) | Exact (generic weights) | costate C2–C5 | the cost law belongs to the object, not to the bookkeeping |
| 20 | fresh-weight Frobenius lemma; K(l); the books close to 6 % at 1024; the ε² law | Exact + Measured | region | an error calculus for any design |
| 21 | P − Q cancellation 10²³; hot gates and weak path disorder at L = 16 | Exact + Measured | tropical | the max-plus skeleton carries no information here |
| 22 | single-site Stein formulas; CMI = price of a separator | Exact | markov P3, §2.3 | memory is a property of the resolution, priced exactly |

---

## 5. Open directions at the level of principles, ranked

**D1. The dilation/Perron sector as the carrier of old content.** (Most promising; unproven at 1024.)
- **Principle.**
  - Positive homogeneity gives a gauge action that commutes with every arrow; its invariant sector, the scale-mixing law, is exactly conserved (C6).
  - The mean direction is the one Perron outlier of the gated propagators (foundations-facts F2.2: squared overlap 0.85 at n = 1024, gain 3–8× the bulk).
  - Isserlis makes the paired channels coherent.
- **Evidence.**
  - In the readout, an exact diagonal plus a rank-1 off-diagonal recovers 84 % of the old-content gain (costate, width 128).
  - The mean mode holds 51–90 % of the slice-containing old content's own D21 energy at ages ≥ 6 (F2.3, width 128). This is in tension with the old-content stream's "almost no D21" (T1).
  - The coherent node κ4 (bethe); λ3 ∝ t (foundations); one depth-growing λ mode in the public chain (F75).
- **Bound.**
  - The sector's third-order content lies inside what the (2,1)-slice chain carries. In the chain's own measurements, slice-level memory plus exact young sources bottoms out at adjusted 2.6e-8 (F65 keep = 8; K3 level, pre-Strassen prices), against 2.2e-8 for the full chain, i.e. 15–25× above the bar. On that evidence the sector alone cannot reach the bar. It also carries κ4 content at the level of the law, which F65 does not test.
  - Its value is to carry, at O(n²), the part of memory that the readout reads most. The unpaired, high-rank remainder must still reach ≲ 10 % of itself (ε² law).
- **Deciding measurements.**
  - At n = 1024 on six networks:
    - costate's gl/gsl against all pairs;
    - the D21 carried by the μ-direction rank-one part of the merged old tier, by age, in one convention (foundations-facts §15 T1, §16 item 2);
    - the fraction of the exact-minus-mean-gate residual slice explained by sym(m ⊗ C) and by σ²·1 ⊗ (Wᵀt).
  - If the sector carries ≳ 90 % of old D21 at ages ≥ 4, then exact young pairs plus the conserved charge is the representation to build. If it carries ≲ 60 %, this direction reduces to a cost trim.

**D2. Collective-coordinate (covariance-mixture) states**, generalising the dilation from one coordinate to k (Hubbard–Stratonovich / localisation, unlock 20).
- Third-order content of the form Σ_k sym(b_k ⊗ X_k) is carried by k + 1 sandwiches per layer, which is the target form named in CONVERGENCE.
- Conjecture, with a caution. foundations-facts F6.9 shows that the (2,1) slice does not transform congruently (|cos| ≤ 0.008), and F9.3 that old content needs 0.19–0.5 n modes for 2 %. So only content that the mixture itself generates can be carried this way.
- Test: the readout-weighted share of old (2,1) energy captured by sym(b ⊗ X) for k = 1, 2, 4 and 8 at n = 1024.

**D3. Goal-oriented accuracy allocation.** The bilinear remainder says the question model needs ≈ 6 % relative accuracy, not 2e-4 absolute. Region's K(l) and faces' T3 price each (age, layer). Allocating exact pairs by K(l) × the age share can cut cost but cannot change its order (C3).

**D4. The centred wall identity** (unlock 37). E relu(Z) = μ P(Z > 0) + E[δ(Z) Γ_Z], where Γ is a two-replica overlap of gradients.
- It is the only per-unit object in the ledger that is not a moment.
- Its conjecture, that the centred wall correction is O(n^{-1/2}), is untested.
- Framework value is high; near-term competition value is low.

**Closed by this round, with reasons:**
- the tropical skeleton and temperature (tropical);
- determinants, signings and the Z/3 sketch (signings);
- covers and the Bethe lift limit (signings, interpolation);
- AMP/Onsager along depth (interpolation B1);
- low-rank or Tucker compression of the double edge (costate C4, with public F46/F88);
- the slice chain alone (costate C3, F65);
- windows as the memory strategy (F61 and every design);
- codimension-0 face closures (faces T1);
- gate/gradient decoupling (faces T4).

**For the user.** On the evidence, a fresh design that reaches the bar will contain most of the first-order local law of foundations-facts F6.3, and that is chain mathematics whatever principle derives it. At n = 1024, with exact inputs, the Wick births alone leave 4.0–4.5 % of D21 per step, and the full list including the (2,1,1) slice leaves ≈ 1 %. The principles' remaining job is the representation of memory, where the chain is expensive (its old tier costs 107 units, F86) and where D1 is the only lead. Whether re-derived chain mathematics is admissible under the 1 Oct ruling is the user's decision. Judgement: on the present evidence it is the only reading under which a fresh-slate design has a visible route to 1.6e-9.

---

## 6. Corrections and stale statements found while judging

1. **Stale calibration in the unlocks note.** Unlock 3 still quotes Gaussian closure "raw 4e-5 at n = 1024" and "the next order must be right to about 2 %"; unlock 36 says "a question model with ≈ 2–3 % relative error meets the bar". With the corrected 4.1–4.3e-6, the rms closure error is ≈ 0.066 n^{-1/2} and the bar ≈ 0.004 n^{-1/2}, so the first-order correction must be right to ≈ 6 %, consistent with region's ε ≤ 6.3 %. Heisenberg §4 and tropical §5 quote the same 4e-5, which their own later verdicts superseded.
2. **Heisenberg R2–R3 (width 64).** "The residual after first order is covariance drift from κ4 (pppq) and (ppqq)." At 1024 the (3,1) slice is unnecessary (3e-10) and second-order covariance Edgeworth is not needed. The (2,2) slice into the covariance and the (2,1,1) slice into D21 are needed. The width-64 picture should not steer a second-order programme.
3. **Bethe §9 (width 64).** "Pair parts are not binding; node beliefs are." At 1024 the covariance (pair) channel is 71 % of the closure's error (region), and the node κ4 is needed only to ≈ 40 % (foundations-problem 3.2). The order parameter is real at 1024, but CONVERGENCE's "third binding fact" overstates its priority at that width.
4. **Heisenberg §2, "the Heisenberg advantage".** Old content at full rank, n³ per pair and no forward tensors is the public chain's existing contraction pattern. The adjoint ordering is its dead end #20, and costate C3 shows that no natural ordering changes the cost order.
5. **BRIEF §3, "no low-rank representation below ≈ 0.3 n modes", against costate's "rank-1 spike".** Both hold. The spike carries most of the old slice's energy at depth (width 128); the 0.3 n figure is for 2 % fidelity on the content's own D21 (0.19–0.95 n by age, foundations-facts C3). The dossier should state both.
6. **Interpolation's 58 % self-averaging against public F47** (covariance propagation plus an online mean correction: 1.0×). Unreconciled. §3.7(c) gives the likely explanation.
7. **The n^-0.8 law and the 1024 first-order closure (F6.3)** predate the B3 = ½ and B7 corrections (foundations-facts F6.5 and §15 item C6). Quote them as preliminary, from one network.

---

## Appendix A. Check of §2 (numpy, runs in seconds)

```python
# First-order response of kappa(a_i,a_i,a_j), a = relu(z), to kappa(z_i,z_i,z_j) at zero correlation.
# Perturbation z_j = t_j + g2 + c (g1^2 - 1): kappa(z_i,z_i,z_j) = 2c + O(c^2); every other cumulant moves at O(c^2).
import numpy as np
from math import erfc, sqrt, pi
g = np.linspace(-12, 12, 400001); pg = np.exp(-g*g/2)/sqrt(2*pi)*(g[1]-g[0])
Phi = np.vectorize(lambda t: 0.5*erfc(-t/sqrt(2))); phi = lambda t: np.exp(-t*t/2)/sqrt(2*pi)
def k21(ti, tj, c):
    ai = np.maximum(ti + g, 0.0); u = tj + c*(g*g - 1); h = u*Phi(u) + phi(u)   # z_j direction integrated exactly
    mi, mj = (pg*ai).sum(), (pg*h).sum(); return (pg*(ai - mi)**2*(h - mj)).sum()
for ti in (-1.0, 0.0, 1.0, 2.0):
    resp = (k21(ti, .3, 1e-4) - k21(ti, .3, -1e-4))/2e-4/2
    Ea = ti*Phi(ti) + phi(ti)
    print(ti, resp, Phi(.3)*(Phi(ti) - Ea*phi(ti)), Phi(.3)*Phi(ti)**2)   # numeric = formula != mean-gate
# Output: -1: 0.0856 0.0856 0.0156 | 0: 0.2106 0.2106 0.1545 | 1: 0.3579 0.3579 0.4374 | 2: 0.5368 0.5368 0.5901
# The diagonal check (z = t + g + c (g^2 - 1), kappa3(z) = 6c) gives 0.131 / 0.341 / 0.528 at t = -1 / 0 / 1,
# matching the formula of §2, against Phi^3 = 0.004 / 0.125 / 0.596.
```
