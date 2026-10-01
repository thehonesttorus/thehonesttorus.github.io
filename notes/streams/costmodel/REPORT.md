# Cost model: closure-based K=3 chains under the flopscope 0.12.1 meter (n = 1024, L = 16)

*2026-10-01. Everything below was metered with the in-process flopscope 0.12.1 engine (whestbench 0.16.1 venv) on this container: 4 cores shared with two atlas jobs, `OPENBLAS_NUM_THREADS=1`. 1 unit = one dense 1024³ matmul = 2³¹ FLOPs; B = 2⁴¹ = 1024 units. Files in this folder:*

| file | what |
|---|---|
| `cost.py` | the calculator: per-layer and total units, C/B, flopscope calls, residual and memory estimates for designs (i)–(v); `--check` re-verifies every formula against the metered JSONs |
| `probe_ops.py` → `ops_measured.json` | 119 metered forms: matmuls, batches, Grams, sandwiches, Strassen L0–L5, Hadamard-then-contract, rank-r families, Gaussian weights, elementwise |
| `probe_designs.py` → `designs_measured.json` | metered skeletons of one middle transition of each design, plus the symmetry-check hazard and the split-sign Gram |
| `audit_v29.py` → `v29_ledger.json` | 504aldo's namespace-tagged V29 (`estimator_v29_ns.py`, MIT) run under the meter on a random He MLP, ledger by layer × family |
| `strassen_v29.py` | V29's `_Strassen` class copied verbatim (MIT), used by both probes |

## 1. Headline

1. **The V29 ledger reproduces exactly under the meter.** Steady-state call: **260.06 units = 0.2540 B, 13,121 flopscope calls**; the first call of a process, which runs the staged L4 Strassen, came to 273.47 units. The published numbers are 260.1 units, 273.5 units and 13,121 ops. The bill does not depend on the data (a random He MLP was used). The calculator's structural replay of V29 gives **256.6 units (−1.3 %)**, and every layer is within 0.42 units of the metered one.
2. **Only the n³ products matter, and each product type has a measured price.** Dense f32: 0.9995. f64: 1.999 (any f64 operand doubles the bill). Batching changes the call count but not the price. Aliased Gram: 0.5002. Strassen-Winograd at L1–L5 through the V29 kernel: 0.877 / 0.772 / 0.683 / 0.609 / **0.5555** per product. A batched family costs **113 calls at L5, whatever the batch size** (67 at L3, 1 dense). My closed-form Strassen recursion reproduces the meter to 4 digits for single products, batches and the hub kernel.
3. **The cheapest symmetric sandwich Wᵀ diag(p) S diag(p) W is 1.027 u**: one Strassen product S X plus the 3 independent half-size output blocks (214 calls). The other forms:
   - factor-Gram, when S = UUᵀ is available: 1.056 u (114 calls);
   - two Strassen products: 1.111 u;
   - tagged 3-operand einsum: **1.504 u in only 3 calls**;
   - dense: 2.0 u.

   The tagged einsum carries a hazard. Its output is checked with `np.allclose(atol=1e-6, rtol=1e-5)` on the float32 result. On an exactly symmetric, *indefinite*, O(1)-scale S it raises `SymmetryError` (a rounding asymmetry of 7e-6 on an output of magnitude 7.3). It passes on SPD inputs and at entry scale ≤ 0.1. **Use it only for covariances (V25 shipped it there), never for modes.**
4. **Closure transport collapses to a small number of products.** Every identified closure diagram (Wick, B1, B2, B3, B5, the r = 1 κ₄ core) is a star with one centre and two legs. The transport of a star to D21(l+1) needs:
   - one a-leg product G = Eᵀ diag(x) W per distinct leg type (E, x);
   - one b-leg product E·Lt per leg type, where Lt sums every term that shares the leg;
   - **one single final product Lᵀ W**, into which every b-index contribution merges.

   The identified closure has **4 leg types**: (C,Φ), (C,w2), (D21zᵀ,w2), (D21z,Φ). With the two slice products and the covariance, a middle layer costs **7.17 units at Strassen L5 (12.8 dense)**. The calculator agrees with the metered skeleton to 0.02 units at L0, L3, L4 and L5.
5. **The design table at L5 (§4):**

   | design | units | C/B | calls | memory |
   |---|---|---|---|---|
   | (ii) closure + slices | 104 | 0.102 | 8,300 | 2.1 GB |
   | (ii) leading Wick only | 59 | 0.057 | | |
   | (iii)/(iv) with r = k = 8 dense modes | 280 | 0.274 | 10,900 | 5.6 GB |
   | (iii)/(iv) with the 8 modes in a shared q = 256 basis | 137 | 0.134 | | |
   | (v) K=2 | 19 (26 dense) | | | |

   Dense symmetric modes cost **1.57 units per mode-layer**: a 1.01 u sandwich plus a 0.503 u birth Gram. With eight of them the chain (280 u) already exceeds V29's whole bill. **Within the leaders' 0.145–0.16 B, a closure chain can afford the full first-order closure (104 u) plus about 45–60 u for old-source and κ₄ content. That fits modes kept in a shared basis, not dense modes.**
6. **Residual time is the binding engineering constraint, not FLOPs.** Residual scales with the number of flopscope calls:
   - on this box: 0.042 ms per call. V29 itself measured **0.546 s here (over the 0.4 s cap)**, against 0.25–0.32 s on its author's box;
   - on the grader, scaled by V29: about 0.022 ms per call.

   Design (ii) at L5 needs about 8,300 calls (≈ 0.18 s at grader scale, 0.35 s here). At L3 it needs 6,000 calls and 126 u; dense it needs 2,100 calls and 179 u. **Every Strassen family costs about 113 calls regardless of batch, so the lever is fewer, fatter families.**

## 2. Measured op table (one call each unless stated; units; residual on this box)

| form | units | flopscope calls | residual ms | note |
|---|---|---|---|---|
| (n,n)@(n,n) f32 | 0.9995 | 1 | 0.04–0.08 | = n²(2n−1) |
| (n,n)@(n,n) f64, or f32 @ f64 | 1.9990 | 1 | 0.06 | one f64 operand promotes the whole op |
| matmul(out=pooled) | 0.9995 | 1 | 0.04 | out= costs nothing extra on matmul |
| W.T @ G (strided view) | 0.9995 | 1 | 0.09 | |
| (n,n)@(n,r), r = 16/32/64/128/256/384 | 0.0156/0.0312/0.0625/0.125/0.250/0.375 | 1 | 0.03–0.04 | linear in r |
| (n,n)@(n,) | 0.0010 | 1 | 0.03 | |
| batched (k,n,n)@(k,n,n), k = 2/4/8/16 | 2.0/4.0/8.0/16.0 | 1 | 0.1–0.6 | no batching discount; calls 1 instead of k |
| broadcast (n,n)@(k,n,n), or (n,n)@(n,kn) | same as batched | 1 | | |
| einsum('ji,jk->ik', X, X), inner(X, X) (same object) | **0.5002** | 1 | 0.09 | symmetric output, tagged |
| X.T @ X (view = another object) | 0.9995 | 1 | | no discount |
| einsum('ji,j,jk->ik', G, d, G) (weighted alias) | **error** | | | `SymmetryError` (float32 check), unshippable |
| weighted Gram Gᵀ diag(d) G, plain | 1.000 | 2 | | |
| weighted Gram, **split-sign aliased** (gather rows with d > 0 and d ≤ 0, scale by √\|d\|, two aliased Grams) | **0.503** | 9 | 1.5 | the cheapest weighted Gram for sign-indefinite d |
| as_symmetric (n,n) | 0.0034 | 1 | 0.09 | 7n² − 1 |
| symmetrize canonical-copy | 0.0005 | 1 | | |
| **sandwich**: dense X.T @ (S @ X) | 1.9995 | 3 | 0.2 | |
| as_symmetric + einsum('ji,jk,kl->il', X, S, X) | **1.5037** | 3 | 0.3 | SymmetryError hazard (§1.3); S already tagged: 1.5002 in 2 calls |
| einsum('ij,ia,jb->ab', S_tag, X, X) (V25 form) | 1.5002 | 2 | | same hazard |
| untagged 3-operand einsum | 1.9995 | 2 | | |
| Cholesky + Uᵀ X + aliased Gram | 1.6669 | 4 | | chol = n³/3 |
| factor given: Uᵀ X + aliased Gram | 1.5002 | 3 | | |
| two Strassen products, L3/L4/L5 | 1.366/1.218/**1.111** | 134/182/226 | 7–10 | |
| **sym3**: S X Strassen + 3 of 4 output blocks, L4/L5 | 1.080/**1.027** | 192/214 | 9 | V29's C_pre form; S X can ride as a slot of another family |
| factor given: Uᵀ X Strassen L5 + aliased Gram | 1.056 | 114 | 5 | |
| diagonal only: sum(X * (S @ X), 0) | 1.0010 | 4 | | the last layer's var |
| sandwich in f64 | 3.0073 | 3 | | ×2 |
| Strassen L0…L5, one (n,n)@(n,n) | 0.9995 / 0.8768 / 0.7715 / 0.6829 / 0.6090 / **0.5555** | 1 / 23 / 45 / 67 / 91 / **113** | 1.7 / 0.8 / 1.4 / 3.6 / 3.4 / 4.6 | f32 relative error vs f64: 3.4e-7 / 5.8e-7 / 1.2e-6 / 1.7e-6 / 2.4e-6 / 3.6e-6 |
| Strassen batch of 4, L3/L5 | 2.732/2.222 | 67/113 | 3.5/6.2 | calls independent of the batch |
| Strassen hub Σ_k X_k Y_kᵀ, k = 4, L3/L5 | 2.723/2.186 | 69/127 | 3.2/5.6 | |
| Hadamard-then-contract ((H∘G)·w) @ W, any order | 1.0005 | 3 | 0.15–0.19 | the Hadamards cost n² each |
| einsum('ij,ij,j,jk->ik', H, G, w, W) | 1.0005 | **1** | 0.05 | same FLOPs, one call |
| type A (G1∘G2∘z)ᵀ @ W; einsum('ka,ka,k,kb->ab') | 1.0005 | 3; 1 | | |
| type A with Strassen L5 | 0.5565 | 116 | 4.6 | |
| G∘G (same object), square(G) | 0.0005 | 1 | | no alias discount for pointwise ops |
| rank-r family, loop of tagged sandwiches, r = 1/4/8/16/32 (SPD modes at scale 0.03) | 1.50 r | 2r + 1 | 0.34–13.5 | hazard on indefinite O(1) modes |
| rank-r family, batched M @ X then Xᵀ @ T | 2.0 r | 3 | 0.2–2.4 | |
| rank-r family, batched einsum('ji,mjk,kl->mil') | 2.0 r | 2 | 0.1–1.0 | a (1,2)-tagged batch gives `ValueError: operands could not be broadcast` (flopscope) |
| rank-r family in a shared q-column basis (transport U, r diagonal cores), q = 64/256 | 0.063/0.25 for any r ≤ 32 | 3 | 0.1 | |
| norm.cdf / norm.pdf on (n,) (billed f64) | 1.0e-4 / 5.6e-5 | 1 | 0.05 | per-layer Gaussian weights Φ, φ, w2…w5, mean: **0.0001 u, 28 calls** |
| norm.cdf / norm.pdf on (n,n) | 0.0469 / 0.0264 | 1 | | 96 / 54 FLOPs per element (f64) |
| exp, arcsin on (n,n) f32 | 0.0078 | 1 | | 16 per element |
| multiply, add, sqrt, copyto, sum(axis) on (n,n) | 0.0005 | 1 | 0.03–0.1 | 1 per element |
| astype f32 → f64 (n,n) | 0.0010 | 1 | | |
| QR (n,64) / (n,384) reduced | 0.0076 / 0.2461 | 1 | | |

Per-call residual: about 0.03–0.05 ms per op on this box for pooled (out=) ops. Fresh-result ops and Strassen internals cost more: the fused-leaf families run 0.04 ms per call.

## 3. The closure terms and their cheapest forms

Notation: W is the layer weight in the stored (in, out) layout, rows i are post-activation neurons of layer s, and columns a are pre-activations of s+1. The target is D21(s+1)_ab = Σ W_ia W_ja W_kb κ3(a_s)_ijk. Under the sym3 role assignments, a star diagram (centre c with weight x_c, legs E_cp to leaf p with weight x_p and E_cq to leaf q with weight x_q) transports in three ways:

- **Centre on b.** Σ_c W_cb x_c Gp_ca Gq_ca, with Gp = E diag(x_p) W. This is (x∘Gp∘Gq)ᵀ W and goes into **L**.
- **Centre on a, leaf q on b.** Σ_c (x_c W_ca Gp_ca) Σ_q E_cq x_q W_qb. This is [x_q ∘ (Eᵀ Lt)]ᵀ W, with Lt = x∘W∘Gp, and also goes into **L**.
- **Centre on a, leaf p on b.** The same with p and q swapped.

Everything lands in one matrix L, and **D21(s+1) = Lᵀ W is one product**. Per distinct leg type (E, x) the transport pays one a-leg product G[E,x] and one b-leg product E·Lt. All terms that share a leg sum their Lt before the product. D3(s+1) = diag(D21(s+1)) costs nothing extra.

| closure term (oracle_k3 basis, coefficient) | what it needs | n³ products | cheapest metered form, units per transition (L5 Strassen / dense) |
|---|---|---|---|
| covariance C_pre(s+1) = Wᵀ C_a W | symmetric sandwich | 1 + ¾ | **sym3 1.03 / tagged einsum 1.50** (3 calls; acceptable on SPD C_a). Last layer: diagonal only, 0.556 / 1.0 |
| exact slices S21 = κ3(a)_(2,1), D3 | S21 W and S21ᵀ(W∘W) into L; D3∘W∘W into L at n² | 2 | 1.11 / 2.0 (two slots of the A-family) |
| leading Wick Σ_cyc Φ_iΦ_j w2_k C_ik C_jk | leg (C_off, Φ) | 2 | 1.11 / 2.0 |
| B1 sym[w2 Φ w2 D21z_ik C_jk] (×3) | leg (D21zᵀ, w2), plus (C,Φ) | 2 | 1.11 / 2.0 |
| B2 sym[w3 Φ Φ D21z_ik C_ij] (×3) | leg (D21z, Φ), plus (C,Φ) | 2 | 1.11 / 2.0 |
| B3 sym[w5 D3z_i Φ Φ C_ij C_ik] (×1) | rides on (C,Φ): centre weight only | 0 | n² |
| B5 sym[w3 w2 Φ K22_ij C_ik] (×1.5), K22 = λ C_off regenerated | leg (C_off, w2) | 2 | 1.11 / 2.0 |
| B6 κ4(z)_(2,1,1) as the published r = 1 core u_i C_jk (×1.5) | v_a (W^TΦCΦW)_ab: v folds into the (C,Φ) Lt as a column scaling; the diagonal term is n² | **0** | n² |
| B6 as a rank-r family Σ_m u^m_i M^m_jk | r sandwiches Wᵀ Φ Mᵐ Φ W. The same product transports the family and gives its D21 term (v^m_a S^m_ab + v^m_b S^m_aa at n²). r births (weighted Grams G_Cᵀ diag(d_m) G_C, which land directly in the pre-activation space of s+1) | r × (1 + ¾) + r Grams | dense modes: r × (1.01 + 0.503). Shared q = 256 basis: 1.4–3.2 per layer for r = 1–16 |
| B0 Φ³κ3(z) all-distinct (old-source memory) | k carried covariance-response modes (x^m, N^m): N^m → Wᵀ Φ N^m Φ W, x^m → Wᵀ(Φx^m). D21 term at n². Births as above | k × (1 + ¾) + k Grams | as B6: 1.57 per mode-layer dense; 1.4–3.2 per layer in a q = 256 basis (k = 1–16) |
| B4 Gaussian ρ³ | path part: leg (C∘C, x). **Triangle ρ_ij ρ_jk ρ_ki: no star form, ≈ 2n⁴ FLOPs ≈ 1024 u per transition** | 2 (paths); triangle unaffordable | 1.11 / 2.0 for the paths. The oracle says the whole ρ³ term buys ≤ 1 point of ε: drop it |
| final Lᵀ W | | 1 | 0.556 / 1.0 |
| Gaussian weights Φ, φ, w2…w5 (n-vectors), mean | `stats.norm` (float64 billing) | 0 | 0.0001, 28 calls |
| gate statistics on n×n (norm.cdf, arcsin, exp) | elementwise | 0 | 0.008–0.047 each |

Grouping that the skeletons and the calculator use:

- **A-family** (one Strassen family): the a-legs, the C_a W slot, the mode slots and the two slice products;
- **B-family**: the b-legs;
- the final product;
- **3-block family**: the symmetric output blocks of the covariance and the modes, one level shallower.

That gives four families per transition, about 113 calls each at L5. The metered skeletons match the calculator to 0.02 u. The difference is exactly the skeleton's stand-in for the nonlinearity:

| one middle transition | calculator | metered skeleton | calls (calc + 40 stand-in) | residual here |
|---|---|---|---|---|
| (ii) 4 leg types, dense (L0) | 12.776 | 12.797 | 73 + 40 = 113 / 118 | 3.9 ms |
| (ii) L3 | 8.809 | 8.830 | 381 + 40 / 426 | 16.8 ms |
| (ii) L4 | 7.859 | 7.881 | 499 + 40 / 544 | 23.2 ms |
| (ii) L5 | 7.167 | 7.188 | 611 + 40 / 656 | 28.6 ms |
| (ii) Wick leg only, L5 | 3.819 | 3.838 | 581 + 40 / 620 | 24.1 ms |
| (v) cov sym3 L5 / tagged einsum | 1.028 / 1.503 | 1.048 / 1.520 | 214 / 3 (+40) | 7.9 / 1.1 ms |

(The rows for designs iii and iv are in `designs_measured.json` and appear in `python cost.py --check`.)

## 4. Design table

Defaults: the identified closure (4 leg types), covariance by sym3, Strassen L5, r = k = 8, one birth Gram per mode per layer. Residual is estimated two ways:

- "here": calls × 0.042 ms, as measured for V29 on this box;
- "grader": calls × 0.022 ms, which is V29's published 0.25–0.32 s for its 13,121 calls.

Memory counts Strassen pools shared by geometry (V29 style) plus state. The real V29 peaks at 5.5 GB; the cap is 8 GB.

**Totals at Strassen L5** (`python cost.py --q 256 --layers`):

| design | units | C/B | calls | resid here | resid grader | mem GB |
|---|---|---|---|---|---|---|
| (i) V29, structural replay (metered: 260.06 u, 13,121 calls, 0.546 s here) | 256.6 | 0.251 | 13,121 | 551 ms | 289 ms | 5.5 |
| (ii) memoryless slices + closure births, D21/D3 transported | **104.1** | **0.102** | 8,304 | 349 ms | 183 ms | 2.1 |
| (ii) leading Wick leg only | 58.9 | 0.057 | 7,896 | 332 ms | 174 ms | |
| (ii) 3 leg types (without B5) / 5 (with B4 paths) | 89.0 / 119.1 | 0.087 / 0.116 | | | | |
| (iii) (ii) + rank-8 (2,1,1) family, dense modes | 280.5 | 0.274 | 10,912 | 458 ms | 240 ms | 5.6 |
| (iii-sb) the same 8 modes in a shared q = 256 basis | **137.5** | **0.134** | 9,386 | 394 ms | 206 ms | 2.1 |
| (iv) (ii) + 8 covariance-response modes, dense | 280.5 | 0.274 | 10,912 | 458 ms | 240 ms | 5.6 |
| (iv-sb) the same in a q = 256 basis | 137.5 | 0.134 | 9,386 | 394 ms | 206 ms | 2.1 |
| (v) K=2 covariance chain | 18.7 | 0.018 | 4,332 | 182 ms | 95 ms | 0.5 |

Designs (iii) and (iv) carry the same machinery (symmetric modes with a sandwich transport and Gram births), so their costs coincide at equal mode counts. They differ in what the modes hold and in how many birth terms are injected: `--births 3` covers the Wick, B1 and B2 cross terms, which give 393 u dense or 151.5 u in a basis at 8 modes.

**Per layer, L5** (L00 = input Gram; L15 trimmed to diagonals):

| design | L00 | L01 | L02–L14 (each) | L15 | total |
|---|---|---|---|---|---|
| (i) V29 replay | 0.6 | 3.3 | 5.7, 8.1, 10.5, 15.2, 17.4, 18.9, 20.6, 21.7, 22.7, 23.7, 24.7, 25.7, 26.7 | 11.4 | 256.6 (metered 260.06) |
| (ii) | 0.6 | 4.0 | 7.38 | 3.6 | 104.1 |
| (iii)/(iv) dense, 8 modes | 0.6 | 16 | 19.7 | 8.0 | 280.5 |
| (iii-sb)/(iv-sb), 8 modes, q = 256 | 0.6 | 6.4 | 9.7 | 3.8 | 137.5 |
| (v) | 0.6 | 1.2 | 1.24 | 0.8 | 18.7 |

**The Strassen-level trade** (FLOPs against calls, i.e. residual):

| design | dense (L0, cov by tagged einsum) | L3 | L4 | L5 |
|---|---|---|---|---|
| (ii) units / calls | 179.1 / 2,116 | 126.2 / 5,990 | 113.0 / 7,358 | 104.1 / 8,304 |
| (ii) with the covariance by tagged einsum (3 calls) | 179.1 / 2,116 | 130.3 / 4,954 | | 110.7 / 6,932 |
| (iii-sb) units / calls | 212.5 / 3,198 | 159.6 / 7,072 | 146.5 / 8,440 | 137.5 / 9,386 |
| (iii) dense modes units / calls | 440.2 / 4,738 | 323.5 / 8,598 | 295.9 / 9,966 | 280.5 / 10,912 |
| (v) units / calls | 25.8 / 1,266 (kit: 24.08 u) | 21.3 / 3,306 | 19.5 / 4,002 | 18.7 / 4,332 (25.3 / 1,378 with the tagged-einsum covariance) |

**Mode-count sweep at L5** (dense modes / q = 256 basis):

| r | units | C/B |
|---|---|---|
| 1 | 126.1 / 125.2 | 0.123 / 0.122 |
| 4 | 192.3 / 130.5 | 0.188 / 0.127 |
| 8 | 280.5 / 137.5 | 0.274 / 0.134 |
| 16 | 456.8 / 151.6 | 0.446 / 0.148 |

Without births, 8 dense modes cost 224.1 u.

### V29 check

The metered ledger (`v29_ledger.json`) by family has young_transport 60.5, hub 54.8, shared 27.8, old_legs 27.8, j_rf 20.8, fb 10.1, j_rot 9.1, j_proj 7.5, cpre 7.1 and j_tier2 6.7 units, the same as the published §2.1 ledger.

The replay prices each family from its shapes:

- young transport: W times 2·min(l−1,3)+2 slots, shared-left Strassen L5. Matches 4.34 u per layer exactly from L05 on;
- hub: 2·min(l,4) terms through the L5 hub kernel; matches 4.36;
- C_pre: 3 half-blocks at L4, 0.471;
- joins: 4 (n,n)@(n,384) + QR(n,384) + projections + core + Qc + rotations;
- tier-2 moves;
- old legs and shared contractions: (n,384) at L3, (n,224) at L2;
- the thin families: per-source constants read off the ledger (0.167 u per source-layer).

What the replay misses: 0.08–0.4 units at L01–L04 in young_transport, about 0.3 units per layer in the tier-2 era (L08–L14) and 0.4 units at L15, about 3.5 u in all. Every term that dominates the bill is reproduced from shapes and measured prices.

## 5. What this means for a chain at the leaders' bill (0.145–0.16 B = 148–164 u)

- **The full identified first-order closure is cheap: 104 u at L5, 126 u at L3.** Each leg type costs 1.11 u per layer at L5 (≈ 15.6 u over the chain). The leading-Wick-only chain is 59 u. The published chain spends 115.6 u on its young tier alone for what is, per the oracle ladder, essentially the same first-order content.
- **Old-source content (design iv) and κ₄ beyond r = 1 (design iii) cannot be carried as dense symmetric modes.** At 1.57 u per mode-layer, 8 dense modes cost 176 u on top of (ii). In a shared q = 256 basis the marginal cost is 1.4 / 1.8 / 2.2 / 3.2 u per layer for r = 1 / 4 / 8 / 16 (33 u for 8 modes over the chain), which lands (iii-sb)/(iv-sb) at **0.134 B**, under the leaders' bill. Whether a shared q = 256 basis carries the content is the accuracy question for the oracle surface. The cost side says it is the only representation that fits.
- **The r = 1 κ₄ core (B6) and B3 are free.** They ride on the (C,Φ) leg. The Gaussian ρ³ triangle is unaffordable at any form (≈ 1024 u per transition) and must stay dropped.
- **Calls, not FLOPs, set the Strassen level.** (ii) at L5 needs 8,300 calls, 0.35 s on this box. Reaching the plan's 2× residual margin (≤ 0.2 s on the grader, i.e. ≤ about 9,000 calls at 0.022 ms; ≤ about 4,800 calls at this box's 0.042 ms) needs one of three things: L3–L4 (+9–22 u), merging families, or the tagged-einsum covariance (3 calls instead of 214; (ii) at L3 then needs 4,954 calls and 130 u). Recommended starting point: L3 with the einsum covariance.
- **Memory.** A Strassen L5 slot costs about 0.23 GB of pooled scratch (0.12 GB at L3/L4). Every slot added to the largest family adds that, so dense-mode designs at L5 approach V29's 5.5 GB peak. The skeleton for (iv) k = 8 at L5, whose families do not share pools, was OOM-killed at 7.0 GB RSS.

## 6. Hazards found (report-worthy to flopscope)

1. **The result-symmetry check of a symmetric einsum is a float32 `np.allclose(atol=1e-6, rtol=1e-5)`.** A mathematically symmetric sandwich `einsum('ji,jk,kl->il', X, S_tagged, X)` raises `SymmetryError("max deviation = inf")` whenever float32 rounding leaves more than about 1e-6 of absolute asymmetry on near-zero outputs. Measured on an exactly symmetric indefinite N(0,1) matrix: output max 7.3, asymmetry 7.0e-6. It passes at scale 0.1 and on SPD inputs, so it depends on the data. In a submission that is a zeroed MLP with multiplier 1.0. The weighted alias `einsum('ji,j,jk->ik', G, d, G)` fails the same way. The reported deviation "inf" is a placeholder in `_pointwise._validate_result_symmetry`, not the measured value.
2. **A batched symmetric operand tagged on axes (1,2) breaks the 3-operand einsum**: `ValueError: operands could not be broadcast together with shapes (r,n,n) (n,r,n)`. So batched tagged sandwiches are not available, and an untagged batched einsum costs the full 2 u per mode.
3. **`stats.norm.*` always bills float64.** This is negligible here (1e-4 u per layer on vectors).

## 7. Reproduce

```
source /root/whest/bin/activate; cd notes/streams/costmodel
OPENBLAS_NUM_THREADS=1 python probe_ops.py ops_measured.json         # ~3 min
OPENBLAS_NUM_THREADS=1 python probe_designs.py designs_measured.json # ~6 min, 7 GB peak at the larger skeletons
OPENBLAS_NUM_THREADS=1 python audit_v29.py v29_ledger_full.json      # ~4 min (two predicts); needs the 504aldo clone (path in the script);
                                                                    # v29_ledger.json here is its layer x family part
python cost.py --check            # formulas vs every metered number
python cost.py --q 256 --layers   # the design table; --lev, --r, --k, --legs, --cov, --births, --q
```
