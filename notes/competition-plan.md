# Competition plan: sixteen days, three lines of attack, one grid

*Written 2026-10-01 against [competition-phase2.md](competition-phase2.md) (rules), [competition-landscape.md](competition-landscape.md) (public frontier) and [infra/azure](../infra/azure/README.md) (the grid). Submissions close 17 Oct 23:59 UTC; team freeze 2 Oct 23:59 UTC; write-up 24 Oct. The plan is revised as results come in; the revision log is at the end.*

---

## 0. The target and the arithmetic

| what | number |
|---|---|
| money line (public board, 1 Oct) | adjusted ≤ 1.6–1.8e-9 |
| equivalent raw at the floor (C/B ≤ 0.10) | ≤ 1.6e-8 |
| equivalent raw at C/B = 0.15 | ≤ 1.1e-8 |
| best public method (504aldo V29) | raw 2.13e-8 at 0.253 B → 5.4e-9 |
| same chain without the old-source tier | 0.150 B but raw 1.3–2.9e-7 |
| reference augmented K=3, un-metered | raw 1.80e-8 |
| ground-truth floor | 7.5e-11 |

So the problem is: deliver the accuracy of an augmented K=3 chain (≈1.5e-8) **or better**, at **≤ 0.12 B**, with zero failures on 100 fresh MLPs. Three lines, in order of expected value per day.

---

## 1. Line A: reproduce, then compress (days 1–4)

1. Run 504aldo's V25 and V29 through our runner on `mini` (100) and `full` (1,000) and on fresh seeds; record raw, C/B, residual time, memory. Expected: raw 2.1–2.3e-8, C/B 0.253/0.367, residual 0.17–0.34 s. This calibrates the grid against the public board (504aldo reports local 8-dump numbers 6.5 % above the board).
2. Re-run the oracle tests on **1,000 networks instead of 8** using keenanpepper's N = 1e9 joint moments (true D21, D3, κ₄ diagonal at every layer): the D21 ε² law, the value of each old source by age, the λ law for κ₄. Cheap on the grid, and it tells us where the 1.5–2× to the leaders can come from with a thousand-network statistical power nobody has used yet.
3. Try the compressions the author did **not** try, each measured against the ε ≤ 2 % rule: (a) a shared basis learned offline (fixed, data-independent subspaces for old sources at each age, fitted on the 14k corpus, so that no range finder is paid online); (b) a **learned memory** for the old tier: predict the old sources' D21 contribution from the young tier's state and the Ω-sketches (the author's +0.02 R² was on 8 networks with cheap features; the sketch targets exist for thousands); (c) a response closure for κ₄ beyond one scalar (SSC-style: a few covariance-response modes with a fitted coupling), evaluated on the same 1,000-network oracle surface.

Line A's deliverable by day 4 is a submission at or below 504aldo's score with a measured margin on fresh seeds, plus the oracle tables.

## 2. Line B: trajectory-fitted corrections at scale (days 3–10)

The census's one positive learned result is a per-layer corrector fitted DAgger-style on the chain's own rolled-forward trajectories; yujigon showed the 64 constants of the cheap family absorb accumulated drift. Nobody has fitted such a corrector on top of a K=3 chain with thousands of networks. On the grid:

1. Roll the chain forward on the d8b corpus (14,048 networks with N = 1e8 marginal targets; 2,048 with the teacher bake), recording the live state at every layer.
2. Fit per-layer corrections sequentially (layer ℓ fitted on states already corrected at layers < ℓ), with the state features that are cheap online (μ/σ, variances, D3, the κ₄ diagonal, the λ observable, the shared-basis energies), frozen before any evaluation; validate on `mini`/`full` and on fresh seeds; report the public-split number last.
3. Keep every fitted table small and shape-guarded (the smoke-test trap), and keep the derivation of what each table absorbs.

Expected value: the Phase 1 analogue gave 8.4× over plain kprop; on a chain already at its truncation floor the honest expectation is 1.2–2×, which is the whole gap.

## 3. Line C: the structural line (days 1–16, in the background)

This is the prong-2 realization of the programme in [research-program.md](research-program.md), kept as a separate line so it cannot be dragged by Lines A and B and they cannot be dragged by it. The public forensics say the Phase 1 front was "deterministic, non-moment, angular"; the Euler–Stein facet identity writes E[h] as a sum over activation cones; our framework's objects are exactly those cones (faces) and the arrows between them. Concretely:

1. Measure on the Phase 2 shape what the earlier notes measured on toy networks: the face process (activation patterns per layer) on He-init 1024×16 networks under the Gaussian state, its memory, the coboundary residual of the derived kernel, the overlap-continuity of arrows. The Azure atlas jobs can record pattern statistics alongside moments at no extra cost.
2. Ask the question the field has not: is the final-layer mean better represented by a **gate decomposition** E[h_L] = Σ over gate configurations of the last two layers than by moments? Gate probabilities and gate–activation cross moments (`gate_GG`, `gate_GX` in the community bakes) are the natural state; a closure in that state is a different projection, and the census's open problem O3 (the joint law of near-threshold neurons) is exactly where the moment chains lose accuracy.
3. Anything from Line C that produces a number competitive with Line A is promoted into the submission path; anything that does not stays as a write-up contribution (the algorithmic-contribution prize is discretionary and judged on the write-up).

### 3.1 The oracle ladder (what the all-distinct third cumulant is made of)

The published chain's interface to each nonlinearity is the (2,1) slice D21(l+1)_{ab} = κ3(z_a, z_a, z_b) of the next pre-activation, with the error law extra MSE ≈ 4.2e-6 ε² (ε = relative rms error of D21; the frontier needs ε ≤ 2.2 %). D21(l+1) is the transport W_{l+1}ᵀ⊗3 of the full post-activation cumulant κ3(a_l), whose all-distinct part (i ≠ j ≠ k) no chain can store at n = 1024. So the structural question is: *what is the all-distinct part of κ3(a) made of, and how little of it does D21(l+1) need?* The tool is [experiments/oracle_k3.py](experiments/oracle_k3.py) on full third-moment atlases built by [experiments/moment_atlas_np.py](experiments/moment_atlas_np.py) `--k3` (exact at small width; the transport identity holds at the sample level, so every ε below is a pure representation error plus Monte Carlo noise in the target).

Width 32, depth 5, N = 4e5, one MLP (ε of D21(l+1), all entries; the off-diagonal-only numbers are within 0.01 of these; atlas `atlas_smoke32k4`):

| layer | slices only ("memoryless") | + leading Wick term Σ_cyc Φ_iΦ_j w2_k C_ik C_jk + Φ³κ3(z) | + Gaussian ρ³, ρ⁴ terms | + fitted first-order cumulant diagrams (κ3 and κ4 hyperedges) | + top n/4 hub columns of the Wick residual |
|---|---|---|---|---|---|
| 0 (Gaussian input) | 0.242 | 0.020 | 0.011 | 0.008 | 0.014 |
| 1 | 0.479 | 0.109 | 0.111 | 0.035 | 0.049 |
| 2 | 0.560 | 0.114 | 0.116 | 0.044 | 0.060 |
| 3 | 0.595 | 0.164 | 0.176 | 0.062 | 0.059 |

Readings. (i) The slices alone lose 24–60 % of D21(l+1): the all-distinct part is not optional (this is 504aldo's dead end 2, now quantified at the tensor level). (ii) The leading Wick term, which every factorised K = 3 chain carries implicitly in O(n³) (G = Wᵀ diag(Φ) C, then Σ_k W_kb w2_k G_ak²), recovers most of it: 2 % at the Gaussian input layer, 11–16 % deeper. (iii) The deeper residual is **not** higher-order Gaussian structure: adding the ρ³ and ρ⁴ Hermite terms (checked against an exact quadrature reference in `--selftest`) changes nothing, so the residual is non-Gaussian, i.e. it is carried by the cumulants of z entering through diagrams the leading term does not have. (iv) A least-squares fit of the residual on seven candidate diagram tensors identifies it: the ρ³ Gaussian term enters with coefficient 1.0 (as it must), the "D21(z) hyperedge + one C edge" term sym[w2_i Φ_j w2_k D21(z)_ik C_jk] with coefficient 2–3 and the "(2,1,1) fourth-cumulant hyperedge" term sym[w2_i Φ_j Φ_k κ4(z)_iijk] with coefficient 1.1–1.4 at every layer; the fit explains 82–86 % of the residual energy and leaves 3.5–6 % of D21(l+1), with the Monte Carlo noise floor of this atlas still to be subtracted. So to a few per cent the all-distinct κ3(a) is the sum of the first-order-in-cumulant gate diagrams: Gaussian ρ² Wick + Φ³κ3(z) + [D21(z) ⊗ C] + [κ4(z)_(2,1,1)], and nothing else of comparable size. (v) After the Wick term the residual is compressible (hub n/4 of the residual beats hub n/4 of the raw all-distinct part by 2–3×), which the raw tensor is not — the compressible object is the *residual after the closure*, not the tensor.

**Width 128, depth 16, N = 5e5 per atlas, two MLPs, two independent atlases each** (so the Monte Carlo noise of the target is measured and subtracted; `--pair` mode). ε of D21(l+1), noise-corrected, evaluated on the atlas the model was not built from; ranges over the two MLPs:

| layers | noise of one atlas's D21 | slices only | + leading Wick | + Gaussian ρ³, ρ⁴ | derived closure without the κ4 term | fitted first-order diagrams (no κ4 term yet) | top n/4 hub columns |
|---|---|---|---|---|---|---|---|
| 0 (Gaussian input) | 0.06 | 0.27–0.28 | 0.04 | 0.04 | 0.04 | 0.04 | 0.18–0.19 |
| 1–3 | 0.025–0.045 | 0.43–0.50 | 0.08–0.09 | 0.08–0.09 | 0.07–0.08 | 0.07–0.08 | 0.18–0.25 |
| 4–9 | 0.01–0.02 | 0.50–0.63 | 0.06–0.10 | 0.06–0.10 | 0.08–0.10 | 0.06–0.08 | 0.08–0.18 |
| 10–14 | ≈ 0.01 | 0.55–0.70 | 0.05–0.08 | 0.04–0.08 | 0.08–0.09 | 0.03–0.07 | 0.06–0.09 |

Readings at width 128 (first batch, no fourth-cumulant slice). (i) The all-distinct part of κ3(a) carries 50–70 % of D21(l+1) at every depth and this does not shrink from width 32 to 128: the "memoryless" chain is wrong by a factor two on its own interface, so any chain that reaches 2e-8 must carry the all-distinct part somehow (the published one does, in factorised form through its "sources"). (ii) The leading Wick term leaves 5–10 % at depth, shrinking only slowly with width (layers 1–3: 11–16 % at width 32, 8–9 % at width 128, roughly n^(−1/4)). (iii) The Gaussian ρ³, ρ⁴ terms buy at most one point; the residual is non-Gaussian. (iv) Hub columns of the raw all-distinct tensor are useless at depth 0–3 (0.18–0.25) and only reach 6–9 % at depth 10–14; the compressible object is the residual after the closure, not the tensor.

**Second batch: atlases with the (2,1,1) fourth-cumulant slice κ4(z)_iijk** (`--k4`; seeds 3 and 4; raw tables in [experiments/results/oracle-k3/](experiments/results/oracle-k3/)). ε of D21(l+1), noise-corrected, cross-evaluated, both MLPs:

| layers | + leading Wick | derived first-order closure (leg-partition coefficients, κ4 slice included) | same closure with the slice regenerated as u_i C_jk (the published chain's r = 1 core) | fitted first-order closure (same seven diagrams, coefficients fitted per layer) |
|---|---|---|---|---|
| 1–3 | 0.08–0.09 | 0.037–0.052 | 0.06–0.07 | 0.036–0.050 |
| 4–9 | 0.06–0.10 | 0.035–0.048 | 0.045–0.075 | 0.018–0.037 |
| 10–14 | 0.045–0.085 | 0.039–0.058 | 0.04–0.06 | 0.014–0.023 |

Readings. (v) With the (2,1,1) slice present, every fitted coefficient agrees with the leg-partition counting at the shallow layers: D21(z)-hyperedge-plus-edge 2.97–3.06 (theory 3), its w3 partner 2.9–3.05 (3), the ρ³ Gaussian term 0.9–1.0 (1), the (2,2)-slice hyperedge 1.2–1.3 (1.5), the (2,1,1) hyperedge 1.4–1.45 (1.5). The all-distinct third cumulant of the post-activation is therefore identified: it is the sum of the first-order gate diagrams, Gaussian ρ² + Φ³κ3(z) + [D21(z) ⊗ C] + [K22 ⊗ C] + [κ4(z)_(2,1,1)], and the earlier negative (2,2) coefficient was the missing κ4 slice being absorbed. (vi) With depth the fitted coefficients drift smoothly away from the first-order values (ρ³ term 1.0 → 0.45, (2,2) term 1.3 → 0.5, (2,1,1) term 1.45 → 1.15, the D21 terms 3 → 2.5–3.3): second-order diagrams renormalise the first-order ones, and the fitted closure reaches 1.4–2.3 % at layers 10–14 while the un-renormalised one stays at 4–6 %. The frontier needs 2.2 %. Whether the renormalised coefficients are an ensemble property of (n, L) (then a 7 × 16 table fitted offline is a legitimate data file) or per-MLP is the next measurement: the second MLP's fit, and then the 1,000-network atlas at width 1024. (vii) The (2,1,1) slice is a genuine n³ object: as a family of n symmetric "covariance-response" matrices indexed by the doubled neuron, rank 32 of 128 captures only 69–81 % of it at layers 1–3, though 86–91 % (one MLP) / 56–82 % (the other) at depth is already in the first mode. Its *effect* on D21(l+1) is easier: at depth the r = 1 family and even the C_off regeneration recover most of it (e.g. layer 10, one MLP: true slice 0.031, rank 1 0.036, C_off 0.038, no slice 0.077), while at layers 1–5 nothing below rank 32 gets within 1.5× of the true slice. (viii) The published chain's reported floor is reproduced by this ladder: its closure is first-order K = 3 plus the regenerated κ4 core, i.e. the "u_i C_jk" column, 4–7 % of D21, and 4.2e-6 × (0.05–0.07)² ≈ 1–2e-8 is exactly where its raw error sits.

What this changes in the plan: the interface error of a K = 3 chain is now a quantified ladder with a known floor, and the two levers below the floor are (a) renormalised first-order coefficients (a data table, if they are an ensemble property) and (b) a better-than-rank-1 carrier for the (2,1,1) fourth-cumulant slice at the shallow layers, where its transport as a rank-r family costs r sandwiches Wᵀ Φ Mʳ Φ W (2r units per layer) instead of the n⁴ of the exact slice. Both are oracle-testable at width 1024 on keenanpepper's N = 1e9 joint-moment atlas (D21 and K22 are derivable from its pair blocks; the (2,1,1) slice is not stored there and needs the GPU atlas jobs). The two identified κ3-diagrams need D21(z) (an n² object the chain already has) and the (2,1,1) slice of κ4(z) (an n³ object it does not); their transports to D21(l+1) cost O(n³) and O(n⁴) respectively. Whether the κ4 term can be dropped or regenerated (504aldo's λ-law regenerates the (2,2) slice from C; the (2,1,1) slice is the analogous question) at width 1024, where ρ ~ 0.03 instead of 0.18, is the width-scaling measurement now running at width 128 (two independent atlases per MLP so the Monte Carlo noise floor is measured and subtracted). At width 1024 the test needs N ≥ 1e8 joint moments, which is keenanpepper's N = 1e9 atlas on the grid.

---

## 4. The grid, as used by the three lines

| job type | pool | what | scale |
|---|---|---|---|
| estimator evaluation on parquet shards | CPU pools, 4 vCPU per task | `whest run` with graded caps, wall cap relaxed | 1,000 MLPs in ≈15 min at 100 tasks |
| chain roll-outs with state recording | CPU pools | the chain in plain numpy (not metered), dumping per-layer state | 14k networks |
| fresh-seed bakes | GPU pool | `whest dataset bake --torch`, N = 1e8 (floor 7.5e-10) or 1e9 | 1,000 networks ≈ 125 GPU-hours at 1e8 |
| moment atlases | GPU pool | marginal moments to order 6, gates, dense pair blocks | as needed |
| fits | one large CPU node or GPU | ridge / small nets on recorded trajectories | hours |

Staging first: the public dataset (77 GB), keenanpepper's full-split joint moments (337 GB), the d8b corpus and teacher bake, the Ω-sketch bakes. All Hugging Face → blob, once, in-region.

---

## 5. Submission discipline

- Every candidate runs on `mini`, on `full`, and on ≥ 200 fresh-seed MLPs with N ≥ 1e8 truth before a submission slot is spent; 10 slots per UTC day, failures count.
- Every candidate passes: `whest validate`, the subprocess runner at the graded caps with `--max-threads 1`, a depth-32 and a width-4/depth-2 shape probe, a cold-start `setup()` timing under 2.5 s, peak memory under 6 GB, residual under 0.2 s on every MLP of `mini` measured on an otherwise idle machine. (Measured here on 2026-10-01: the public V25 chain, whose author reports 0.28 s residual, read 0.396 s and 0.410 s on a 4-core box that was also running two bakes, and the second MLP was zeroed. The residual cap is a Python-speed cap; a 1.4× slower interpreter turns a rank-10 estimator into a failed one, so the margin must be at least 2×.)
- Randomness only from `mlp.seed` / `ctx.seed`; no per-MLP adaptive constants fitted on in-sample statistics; every shipped table documented with what it absorbs.
- The designated final submission is chosen on fresh-seed MSE, not on the public board.

---

## 6. Timeline

| days | milestone |
|---|---|
| 1–2 (Oct 1–2) | grid live (needs network access or a machine with `az`); datasets staged; 504aldo V25/V29 reproduced locally and on `mini`; **team freeze by 2 Oct 23:59 UTC** |
| 3–4 | 1,000-network oracle tables; first submission (reproduction with margin) |
| 5–8 | Line A compressions and Line B fits in parallel on the grid; daily submissions of the best fresh-seed candidate |
| 9–12 | integrate; robustness campaign; Line C measurements written up |
| 13–16 (Oct 13–16) | freeze, designate, keep one fallback submission at a known score |
| 17 Oct | submissions close 23:59 UTC |
| by 24 Oct | algorithmic-contribution write-up |

---

## 6b. Campaign status, 1 Oct 17:00 UTC (stream reports in [streams/](streams/))

- **The D21 interface is solved at the real shape, given exact inputs** (oracle1024, preliminary: one MLP, N = 32k × 2 replicas, replica cross-product with jackknife errors; production N = 3.5e6 running). At width 1024 the leg-partition first-order closure, with no fitted coefficient, gives D21(l+1) errors of 0.8–1.0 % at layers 3–14, below the 2.2 % bar; the closure error falls as ≈ n^(−0.8) (deep layers 4.7 → 2.6 → 0.9 % at n = 128 → 256 → 1024) against n^(−0.25) for the leading Wick term. The fitted coefficients equal the leg-partition values within noise at every depth: the drift seen at width 128 is a finite-width effect, so no coefficient table is needed.
- **The binding object is the (2,1,1) fourth-cumulant slice.** With the exact slice the closure reaches 0.9 %; with the published chain's r = 1 regeneration u_i C_jk it stays at 2.5–2.7 % (above the bar), and the regeneration gets relatively worse with width. Independently, the dense width-128 chain (chain128) finds end to end that with κ4 held at truth the D21 ladder translates one-to-one into the final layer (slices only 1e-4 → Wick 8.6e-6 → closure 5.4e-6 → diagram engine 1.9e-6, a 50× range), while with any κ4 closure of its own every κ3 rule lands at 3e-5–8e-5: the fourth-cumulant closure, not the κ3 interface, is the binding error. EscAI's oracle points the same way (joint fourth cumulant worth ≈ 20×). The coefficient stream is redirected to carriers of the (2,1,1) slice.
- **Old-source κ3 content has no carrier clearly cheaper than the published old tier** (old-content, provisional pending width 256): propagator projection needs k ≈ 0.3 n directions (the published shared basis); dense symmetric mode families need r ≥ 16 (≥ 40 units per layer); the only candidates that could undercut the published ≈ 6.7 units per layer are a shared-basis Tucker family (r ≈ 32, q ≈ n/4) and a low-rank family (32, n/8), and only if r stays near 32 at n = 1024.
- **Fallback submission:** the patched V25 bundle runs clean on six competition-shape MLPs on an idle box (0 failures, C/B 0.3667, residual 0.20–0.22 s per MLP, 3 GB peak, MSE minus truth noise 1.76e-8, i.e. adjusted ≈ 6.4e-9). The original V29 and V25 fail the robustness shapes (MemoryError at 1024×32, SymmetryError on adversarial weights); the patch routes non-suite shapes and failures to a float64 covariance fallback. V29 still fails the 0.4 s residual cap in local emulation (0.46–0.48 s idle); a client/server emulation of the grader is running to decide whether that is real.

## 7. Outside evidence and prices (1 Oct, from the parallel campaign)

- **Where the error lives (team EscAI's public oracle, 8 width-1024 networks):** replacing the pre-activation per-neuron variance at every layer cuts raw MSE by 40 %; adding the per-neuron third cumulant reaches 8.4e-9 and the joint fourth cumulant 1.17e-9, while the K3/K4 readout is not the bottleneck (bias 1.2e-10). Keeping the repeated-index κ4 entries cut raw by 35–64 % at widths 128–192. Their covariance-response-mode design was stopped at its first gate on a structural control (congruent transport of the (2,1) slice: |cos| ≤ 0.008); no mode was run. Digest: [digests/phase2-intel-2026-10-01.md](digests/phase2-intel-2026-10-01.md) §3.
- **Prices (measured with flopscope at n = 1024; [streams/costmodel/REPORT.md](streams/costmodel/REPORT.md)):** a dense product costs 1.0 unit, a same-object Gram 0.5, Strassen–Winograd L5 0.556 in 113 calls per family. Every identified first-order closure term except the Gaussian ρ³ triangle transports to D21(l+1) by factorised products: ~7.2 units per middle layer at L5, about 104 units = 0.102 B for slices plus the full first-order closure over the chain. That price does not include carrying the old-source content (the Φ³κ3(z) term), which is the open question of the old-content stream; the only mode form that fits the leaders' 0.145–0.16 B is a shared 256-column basis (137.5 units for 8 modes). The binding constraint is residual time per flopscope call (~0.022 ms on grader hardware, 0.042 ms here), not FLOPs.
- **Rules and grader (reviewed):** a fair-accounting rule since 16 Sep (benefit must come from the estimation method; packing banned; re-scoring possible during prize review); Strassen confirmed permitted (11 Sep); the smoke test runs a 256-wide, 32-deep MLP and any layer-indexed table that assumes depth 16 fails the whole submission; assigning `x.shape` fails on the grader. The Phase 2 algorithmic-contribution prize is $20,000.
- **Azure grid:** 24 defects found offline and fixed (retired node image, autoscale formula, data-plane auth, container task user, GPU image, grid.py crashes, job limits, shared login state); [../infra/azure/CHECKS.md](../infra/azure/CHECKS.md) holds the first-run checklist.

---

## Revision log

- 2026-10-01 (evening), after reading the published chain's code, findings log and ledger ([digests/504aldo-k3-chain.md](digests/504aldo-k3-chain.md)):
  - **Line A, corrected.** Every leg-based compression of the third-cumulant sources is closed by the author's arithmetic (Tucker, CP caps, rank ladders, shared bases below the age law, hub merges, adjoints; "no leg-based representation of this expansion reaches 0.15 × B"). The author's own ceiling: the K3 + memoryless-κ₄ closure floors at raw ≈1.94e-8 on the board scale; at his multiplier 0.253 that is rank ≈9, and rank 5 would need raw 1.5e-8, *below* the floor. So "reproduce then compress" cannot reach the money line; what remains of Line A is the reproduction (a known score with margin) and the 1,000-network oracle tables, whose purpose is now to find which **representation** the leaders use, not to tune this one. The leaders' bill equals his chain minus the old tier (0.15 B), and his reading is that they transport symmetric n×n objects (covariance-response modes) rather than per-source legs.
  - **Line B, corrected.** Fitted corrections on this chain are measured dead at every level (output ridge |corr| < 0.09; online per-layer correction 1.45× on an older base but 1.01× on the regenerated one; slice fits +0.02 R²; a 144-constant ridge 1.07×). Only trajectory fits on a *different* representation could be alive; Line B is demoted to a check, not a line.
  - **What is open, by the author's own gate**: carry fully off-diagonal third-cumulant content below ≈1 unit per old source-layer without forming dense legs, with D21 accurate to ≤2.2 %; and fourth-cumulant content beyond the memoryless core (the exact (2,1,1) core is worth −9 % raw at 0.5–1.0 B; the only untested lever on his chain is a 9-mode regeneration, ≤1 %). The engineering facts (residual clock allocation-bound, pooled `out=` buffers, Strassen permitted, smoke test at depth 32, 5.5 GB peak) transfer verbatim.
  - **Line C, promoted.** The structural question ("not a per-source hub-tensor expansion at all") is now the main line: a representation of the non-Gaussian content in terms of gate statistics (P(z_i>0, z_j>0), E[1[z_i>0](a_j − μ_j)], which the community bakes already store) or of covariance-response modes, evaluated first as an oracle on the 1,000-network joint-moment atlas before any chain is written.
- 2026-10-01 (16:40): section 7 added from the parallel campaign's support workflow (EscAI oracle, measured prices, rules, Azure fixes); the parallel streams are indexed in [streams/README.md](streams/README.md).
- 2026-10-01 (evening): second oracle batch with the (2,1,1) fourth-cumulant slice; the all-distinct third cumulant identified as the first-order gate diagrams with coefficients matching leg-partition counting; the chain floor reproduced as the regenerated-κ4 column of the ladder; raw tables saved under experiments/results/oracle-k3/.
- 2026-10-01 (later): V25 reproduced locally on a 2-MLP, N = 1e6 dev bake: final-layer MSE 1.00e-7 against a truth floor of 7.5e-8 (0.0748/1e6), i.e. chain error ≈ 2.5e-8 at 0.3667 B, consistent with its author's 2.13e-8 on the public mini set. Section 3.1 added with the first oracle ladder (width 32); the Oishi Phase 1 digest saved (its chaos-spectrum numbers are the quantitative reason sampling is out of the Phase 2 race: exactness to Hermite degree 5 buys 3×, degree 32 buys 10.7×, and the frontier needs ~300×).
- 2026-10-01: first version, written before any grid job has run; numbers are the field's, not ours.
