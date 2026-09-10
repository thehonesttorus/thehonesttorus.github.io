# Completeness critique of the literature study

The adversarial critic's report on the synthesis above. It overturns two of the synthesis's own claims
and surfaces the material the study missed.

---

I have what I need. Writing the critique.

# Completeness critique — ARC WhestBench Phase 2 literature study

**Bottom line: the per-paper NCG rulings are sound and I could not overturn one. But the synthesis's headline structural argument is arithmetically wrong, its top reading-list item is disconfirmed by our own data, two of its five nominated "where the work actually lives" directions are killed by our own coherence number, and the study never retrieved the one document that matters most — arXiv:2605.05179, the benchmark's own source paper, which contains a closed form for the exact channel obstruction 5 identifies as the 3.1× lever, plus a public reference implementation.**

---

## 1. The synthesis's central claim is wrong, and it wrote off the right direction

§1 and §5: *"the published oracle for perfect marginals at every layer is 7.07e-7, so a free, exact, level-infinity version of every per-layer state approximation in this corpus lands at 1.23x... the leader is 36.6x beyond that oracle. Whatever the leader is doing, it is not a better marginal closure. That single fact invalidates the corpus at a stroke."*

`results.txt` §9 says both halves of this and they are not the same oracle:

- 7.07e-7 is the chained **second-order** (Gaussian-marginal) oracle. Its one-step version is the `S_2` line, 1.8e-7.
- The one-step `S_3` oracle is 2.2e-8 and the one-step `S_4` oracle is **1.9e-9 — below the leader's 1.93e-8**. `results.txt` states this in terms: *"The marginal operator system S_4 determines the target to 1.9e-09, below the leaderboard frontier."*

So 7.07e-7 caps **order-2** closures, not "every per-layer state approximation." Worse, it does not cap **us**: our estimator carries the full covariance `C` and the two-point tensor `K21`, both joint objects, so we were never inside the family the oracle bounds. The synthesis used a ceiling on a family we do not belong to, to argue that the direction we are already in is closed. That inference is what licensed dismissing the whole corpus "at a stroke" *and* licensed §5's conclusion that the remaining 45× "lives in the joint, in a sampling design, or in exploiting the He ensemble" — i.e. anywhere but where the evidence points.

The evidence points at a better cumulant closure. See §5 below.

---

## 2. Papers accepted too readily — the connection is verbal

**Abanov et al., arXiv:2507.21135, Appendix C. Reading-list item #1 should be struck: our data disconfirms it, it does not predict it.**

The credited claim is that random matrices are "deep quantum," their commutator-Laplacian spectrum is power-law rather than sqrt, so there is no IR/UV separation and low modes carry no privileged share — offered as "the only mechanistic in-corpus prediction of obstruction 1." Run the arithmetic that the synthesis did not:

- `results.txt` §5: the `h_8 → z_15` transport needs **175 singular directions for 99% of its energy**. Rank-128 truncation therefore discards roughly 1–5% of the transport's energy.
- `results.txt` §1: rank-128 truncation costs **5.544e-7 against 2.190e-7 exact — a 153% degradation.**

A "no privileged modes / no spectral separation" mechanism predicts damage proportional to discarded energy. It predicts rank-128 should be nearly free. It is off by roughly **thirty-fold in the wrong direction**. The Abanov mechanism is a *resolution* story; obstruction 1 is measured to be a *conditioning* story (coherence 0.0887–0.1106, sign agreement 0.506). These are different mechanisms and the data separates them. The synthesis correctly withdrew Freeman et al.'s abelian no-go for exactly this failure mode — mispredicting where we measured — and then put a paper with the same defect at the top of the reading list. Withdraw it, or demote it to the same "credits to withdraw" list.

Same test kills the credit given to **Steinacker 1911.03162** ("most of the configuration space is basically 'white noise'") and confirms the synthesis's own caveat on **Schneiderbauer–Steinacker 1601.08007** §7.1: a spectral-gap criterion is blind to the mechanism that binds.

**Flora et al. — "independent reproduction of the signature of our site_pow flatness."** Their plateau is in *k* (how many moments); ours is in the *exponent of an importance weight* (which ranking). Different axes. "Same signature" is a word, not a measurement. Downgrade.

**§3.0's provenance is itself incomplete in the way the study warns against.** The proposed reweighted sparse re-summation is credited to Carathéodory/Steinitz and column-subset selection. Its actual established form is **DEIM** (Chaturantabut & Sorensen, *SIAM J. Sci. Comput.* 32(5):2737–2764, 2010) — an oblique projection `U(PᵀU)⁻¹Pᵀ` that selects `k` indices *and re-solves the weights of survivors*, exactly the construction proposed — and **Q-DEIM** (Drmač & Gugercin, arXiv:1505.00370), whose whole point is choosing the selection by pivoted QR so as to minimise the amplification constant `‖(PᵀU)⁻¹‖`. That constant is the DEIM-side analogue of our 11×. The study proposed an experiment whose error theory and selection algorithm already exist and did not cite them.

---

## 3. Papers dismissed too fast

Honestly: **none of the 42, on the merits.** I tried to break the Latrémolière, Rieffel Thm 8.2, Aguilar–Kaad–Kyed and Connes–van Suijlekom rulings and each holds. Two *arguments* are broken even though the verdicts survive:

**d'Alessandro–Roch i Carceller–Tavakoli (§3.2).** The kill leans on *"we under-estimate κ₃ by 64% at layer 16, i.e. we are too Gaussian, i.e. strictly interior."* The cited measurement does not say that. `results.txt` §8's `|tail|/|k3_depth8| = 0.644` is the share of the depth-8 cumulant lying in ages 4–8 — not our residual error, and it is measured *with the kernel off*. §8 also records `|F|/|tail| = 2.17–2.47`: the deployed correction is **2.2–2.5× larger in norm** than the sector it stands in for. "Strictly interior to the cone" is not established; the deployed state's κ₃ could be over- or under-shot per neuron. Also, §7's 2% positivity measurement is about the **covariance** eigenvalue, not about the (μ, C, κ₃) triple being a valid moment sequence — a strictly stronger condition nobody measured. The verdict survives on the compute argument (a level-2 moment matrix has side 5.25e5) and on the fact that we have no unknown operators. Drop the cone-direction argument; it will mislead the next reader.

**Becker–Li / Freeman et al. purification (§3.1).** The "capped at 2%" claim is right for positivity. But the ruling then asserts the third-cumulant path sum "does not embed in a density-operator picture at all," which forecloses a construction that is real: carrying `S` with `C = SSᵀ` changes the *metric in which any subsequent low-rank step measures error* from `‖ΔC‖` to `‖ΔS‖`. That is not worth chasing here (obstruction 1 kills the low-rank step regardless), but the stated reason is stronger than the evidence.

---

## 4. Retrieval gaps

- **Connes 1988** (never retrieved) and the **~7 uncovered papers**: immaterial. The corpus fails on structural axes that a seventh cluster would not repair.
- **The real retrieval gap is that the study never read arXiv:2605.05179.** The provenance note correctly identifies it (Wu, Lecomte, Winer, Robinson, Hilton, Christiano, ARC) and correctly says no NCG paper critiques it — and then treats it purely as a provenance fact. `finite-resolution-question-algebra.md` §7 admits: *"Only the main text of the ARC paper was read... ARC's K ≥ 2 variants were not re-implemented."* I retrieved the full text including the technical supplement. It is decisive. See §5 and §6.
- **The reference implementation is public and unmentioned anywhere in the study or the notes:** `https://github.com/alignment-research-center/mlp_kprop` (basic / augmented / factorized / factorized-augmented, K = 1..4).

---

## 5. Missing literature, named — and two of the study's own nominations killed

### 5a. What the source paper contains that the study never saw

**arXiv:2605.05179, §4.1 "power cumulants" + footnote 17 (p. 32) + §S.3.3.** The paper's own worked example is *literally our dominant channel*:

```
κ₃[φ(Z)]_{i,i,j}  =  Σ_{k≥1} (1/k!) (ĉ^{φ²Z}_k)_i (b̂^{φZ}_k)_j · Cov[Z]^k_{i,j}
                     − 2 (b̂^{φZ}_0)_i Σ_{k≥1} (b̂^{φZ}_k)_i (b̂^{φZ}_k)_j · Cov[Z]^k_{i,j}
                  =  κ₂[φ(Z)², φ(Z)]_{i,j} − 2 κ₁[φ(Z)] κ₂[φ(Z)]_{i,j}
```

`κ₃[φ(Z)]_{i,i,j}` **is** our `K21_ab = κ(z_a, z_a, z_b)`. The right-hand side is a Hadamard power series in the correlation matrix — computationally identical in shape to the loop already in `/home/user/thehonesttorus.github.io/notes/whestbench/estimator.py` lines ~262–275, which computes `C = Σ_k ρ^k · outer(d_k, d_k)/k! · so`. The only new ingredient is the Hermite coefficients of **ReLU²**, given in closed form in ARC's §S.3.4. `estimator.py`'s `hermite_coeffs()` returns only the coefficients of ReLU itself — **the estimator implements no power cumulants at all.**

Cost of the closed form at `HERMITE_ORDER = 8`, n = 1024, 16 layers: `8 × n² × 16 ≈ 1.3e8` FLOPs = **0.006% of 2⁴¹**. The synthesis §5 declares the binding engineering question to be *"can the deep two-point channel be delivered for under 37% of 2⁴¹ instead of 49%?"* The source paper writes the freshly-generated part of that channel in closed form for a rounding error.

Two more facts from the same paper the study should have had:

- **ARC §6.4 / Fig. 3(b), the power-cumulant ablation:** dropping power cumulants and the `(K+1)`-trace degrades MSE from `O(1/n^K)` to **`O(1)`** — flat in width. These two adjustments are the entire difference between a competitive and a non-competitive estimator, and we implement neither properly (our `g4` is a single scalar per layer; ARC's augmented variant carries the whole 4th-order tensor minus its traceless part).
- **ARC §4.2 / §S.4.2, the augmented algorithm:** for `K ≥ 3` it "improves the constant factor of the leading-order term in MSE, while keeping the leading-order term in FLOPs the same." A free constant-factor improvement at K=3 that we do not have.

### 5b. The leader is almost certainly the reference algorithm, and the study never checked

ARC conjectures `MSE ≲ c_K (L/n)^K` (variance-normalised) and factorized runtime `≲ c''_K L² n^K` for `K ≥ 3`. At n=1024, L=16, neuron variance 7.8e-2:

| | var-normalised MSE | vs `(L/n)^K` | factorized cost `L²n^K` |
|---|---|---|---|
| ours (8.675e-7 raw) | 1.11e-5 | `c₃ ≈ 2.9` | — |
| leader (1.93e-8 raw) | 2.47e-7 | `c₃ ≈ 0.065` / `c₄ ≈ 4.2` | — |
| K=3 factorized | — | — | `256 × 1.07e9 = 2.75e11` = **12.5% of 2⁴¹** |
| K=4 factorized | — | — | `2.8e14` = 128× budget |

**The leader sits at 14.7% compute. A factorized K=3 costs 12.5% × c''₃.** That is not a coincidence worth ignoring. A properly-implemented factorized (augmented) K=3 lands at exactly the leader's operating point, and our 45× gap is the gap between a full K=3 and our heavily approximated one (one exact source age, a rank-one `K21` slice along `muh`, and fitted scalar kernels for everything else). This back-of-envelope costs five minutes and reframes the entire program; the study did not do it.

### 5c. Two of the study's own nominated directions are killed by its own number

**Randomized NLA / leverage-score sampling (§5 item 4) is dead, by two lines of arithmetic the study skipped.** Sample `m` of the `n` per-source-unit terms uniformly, rescale by `n/m`. Variance ≈ `(Σ‖T‖)²/m`, so relative sd of the aggregate = `(Σ‖T‖)/(√m ‖ΣT‖) = 1/(0.09√m) = 11.1/√m`. For 10% relative accuracy, `m ≈ 12,300 > n = 1024`. **Any unbiased term-sampling estimator needs more term evaluations than computing all of them exactly.** Leverage scores give `(1±ε)` relative guarantees on `‖Ax‖` — a *norm*. For a cancelling sum the guarantee degrades by exactly the condition number `1/0.09`. This is the same `⟨s⟩⁻²` law as the sign problem; it is a no-go, not a home.

**The sign-problem direction (§5 item 2) is over-promised.** I searched it (arXiv:2604.24290 review; 2311.13002, 2103.08948 contour optimisation; Langfeld–Lucini–Rago LLR). Contour deformation / Lefschetz thimbles need a *holomorphic parametrisation of the integrand* so the contour can be moved; our terms are indexed by a discrete unit label `β` with no complex structure to deform. Complex Langevin and constrained-path likewise need a stochastic process over a continuum. Only the **density-of-states / LLR** idea survives structurally — compute a smooth density first, do the oscillatory sum in one dimension afterwards — and that, as the synthesis itself notes, is what the rank-one kernel already does. Report the class honestly: one transferable idea and one confirmed no-go, not "the class the corpus is missing entirely."

### 5d. Output-weighted model order reduction is *already measured dead* — say so

The prompt asks for balanced/Hankel-error MOR. Worth stating as a settled negative rather than an open direction:

- For a single input–output map, **balanced truncation ranks by exactly the singular values of that map**. Our transport truncation *is* balanced truncation of `h_ℓ → z_15`. Measured: obstruction 1(a).
- The `site_pow` sweep (`|p_β| · ‖U_β‖^s`, `s = 0/1/3` → 6.844 / 6.791 / 6.772e-7) is a crude controllability×observability sweep spanning the balanced ranking. **Flat to 1%.**
- Every classical MOR error bound (`‖G−G_r‖_∞ ≤ 2Σσ_i`) is a tail-sum, i.e. triangle-inequality, i.e. 11× loose before it starts.
- **Goal-oriented / DWR** (Becker & Rannacher, *Acta Numerica* 10:1–102, 2001; and arXiv:2305.15285, 2207.11233 for the modern discrete treatment) is the correct name for the synthesis's axis (B) — "a seminorm induced by future observables." It is exactly the adjoint weight. And the DWR literature's own known pathology is ours: when the cell-wise contributions are signed and cancelling, taking absolute values to build the estimator inflates the effectivity index. That is our 11×, in the field that owns it.

The one MOR family *not* foreclosed is **interpolatory/moment-matching** (IRKA, Loewner) — match the aggregate at selected points instead of bounding a norm. In our setting it has no frequency variable and collapses into resummation in the source-depth index; see §6 runner-up.

### 5e. Genuinely missing, with identifiers

| Direction | Papers | Why |
|---|---|---|
| **Finite-width perturbation theory — the theory that owns our object** | Yaida, arXiv:1910.00019; Roberts–Yaida–Hanin, arXiv:2106.10165 | Derives the `L/n` recursions for connected 2- and 4-point functions. Our "universal depth schedule" `a_l` (§12: across-network sd 0.00314 below the across-seed spread) is an ensemble-universal quantity these give in closed form. **The 15 fitted `(a, g4, c2p)` triples in `estimator.py` are candidates for derivation rather than Monte-Carlo fitting** — which would also remove their 8.8% seed noise. This literature is absent from all 42 papers and is the same family as 2605.05179 itself. |
| **Renormalized memory kernels** | Price & Stinis, arXiv:1707.01955 (*Renormalized Reduced Order Models with Memory for Long Time Prediction*, MMS 17(1):68–91) | Truncated Mori–Zwanzig memory terms whose coefficients are **renormalised by fitting to short-time data**, and which are found to follow a **universal power law in time**. That is `results.txt` §8 + §12 verbatim, in the field that owns it, with a functional form for the schedule. The synthesis says "no memory kernel anywhere in the corpus" and then names no MZ paper. Also arXiv:2101.05873 (Lin–Lu, data-driven MZ) and arXiv:2601.07101 (history-enriched linear MOR — closure error from instantaneous-only dependence). |
| **Optimal blending of a biased analytic model with an unbiased sampler under a budget** | Peherstorfer–Willcox–Gunzburger, *SIAM J. Sci. Comput.* 38(5):A3163–A3194 (2016); Gorodetsky–Geraci–Eldred–Jakeman, arXiv:1811.04988 | `estimator.py` has `W2P = 0.75` (swept, §19) and `MC_MAX = 2400` (guessed). MFMC/ACV give the optimal weight and the optimal budget split in closed form from correlations and cost ratios — the FLOP meter is exactly their objective. The synthesis calls this "textbook control variates"; it also **mislabels the mathematics** — with a *deterministic* analytic model there is no control variate, only bias–variance shrinkage, and the optimal weight is `Var_MC/(Var_MC + Bias²)`, not `Cov/Var`. |
| **Deterministic low-rank + stochastic residual, with a proven optimal split** | Meyer–Musco–Musco–Woodruff, arXiv:2010.09649 (Hutch++); Epperly–Tropp–Webber, arXiv:2301.07825 (XTrace) | Our estimator *is* this architecture (analytic source + MC residual). Hutch++ proves the optimal budget split between the deterministic sketch and the stochastic remainder. We have never optimised ours. |
| **Reweighted sparse re-summation (for §3.0)** | Chaturantabut–Sorensen, SISC 32(5) 2010; Drmač–Gugercin, arXiv:1505.00370 | Selection *plus* weight re-solve, with an explicit amplification constant to minimise. |
| **Conditioning of summation** | Higham, *Accuracy and Stability of Numerical Algorithms*, ch. 4 | `Σ‖T‖/‖ΣT‖ ≈ 11` is the condition number of a summation. Correctly named by the synthesis; I confirm it is the right frame and that the classical remedies (compensated summation, higher precision) do not apply — our error is truncation, not roundoff. |

---

## 6. The single highest-value experiment the study did not propose

**Port ARC's own reference implementation — `mlp_kprop`, factorized augmented, K = 3 — to the WhestBench harness and profile it with flopscope at n = 1024, L = 16 on the public mini split.**

Why this and nothing else:

1. **It is the only candidate with a 45× ceiling.** Every other item in the study's §3 is capped by measurement at ≤ 2%, ≤ 3.3%, or 1.23×. The conjectured `c₃(L/n)³` scaling and the `L²n³ = 12.5%` cost land on the leader's exact operating point (1.93e-8 at 14.7%).
2. **It is not a guess about what to build — it is a run of code that already exists**, at a URL the study never surfaced.
3. **It settles the strategic question either way.** If it scores near 2e-8 at ~13%, the gap was implementation completeness and we adopt it. If it does not, then a full K=3 provably does not reach the frontier at depth 16 and every remaining hypothesis about the leader is falsified in one measurement — which is worth more than another cluster of spectral triples.
4. **It subsumes the two upgrades we are provably missing**: power cumulants (ablation: `O(1)` vs `O(1/n^K)`) and the augmented `(K+1)`-tensor (free at leading order for `K ≥ 3`), against our single scalar `g4` per layer.

**Do this first, today, as the falsifiable one-hour version:** replace the estimator's leading-order two-point source construction (`A`, `X1`, `X2`, `Nt` in `estimator.py`) with ARC's closed-form power-cumulant series
`K21_ij = Σ_k (1/k!) ĉ^{ReLU²}_k(t_i) · d_k(t_j) · ρ_ij^k · (scale)` (footnote 17, p. 32; ReLU² Hermite coefficients in §S.3.4),
computed in the Hadamard-power loop the estimator already runs for the covariance. Cost ≈ **0.006% of budget** against the 49% the synthesis names as the binding constraint. Score it against the exact depth-8 two-point channel (`results.txt` §2, 2.316e-7). If it matches, the entire "deliver the deep two-point channel under 37%" problem dissolves and the freed budget buys a fourth-cumulant channel.

**Runner-up, if only a zero-FLOP diagnostic is affordable:** the per-layer **oracle-substitution error budget** — propagate with the model to layer ℓ, substitute the true state there, run the model forward, and plot final MSE against ℓ, separately for μ, C, and κ₃. `results.txt` §9 varied the *readout* order under an oracle state; nobody has ever varied the *state* under the deployed readout. The study asserts "the binding constraint is PROPAGATION" seven times without ever measuring **propagation of what**. Forty-two papers were searched for tools without knowing which quantity to apply them to. This is the DWR error budget computed by brute force, it costs no metered FLOPs, and it would have redirected the entire study.

**Two free micro-fixes found in the code while checking, worth one line each:** the shrinkage weight `w` at `estimator.py:233–236` is a *single global scalar per layer* built from `noise = 6σ⁶/N` — the Gaussian-only term of the third-moment variance, ignoring `κ₆ + 9κ₄κ₂ + 9κ₃²`, and summed as if the per-neuron residuals were independent when they share the same 2400 samples; and `rng.standard_normal` at line 149 is plain i.i.d. sampling with no radial Rao-Blackwellisation, despite `E[h] = E‖x‖ · E[h(x̂)]` holding exactly by positive homogeneity (worth only ~0.5% of MC variance at n=1024 — small, but free and exact).