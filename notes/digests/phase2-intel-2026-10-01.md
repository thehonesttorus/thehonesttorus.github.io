# Phase 2 intelligence sweep, 10 Sep to 1 Oct 2026

*Swept on 2026-10-01 across the AIcrowd forum (read through Exa, since the proxy blocks aicrowd.com), the public leaderboard and submission pages (also through Exa), GitHub (repository search, issue trackers, and `git clone` of every relevant public repository), Hugging Face (MCP search over datasets, models, spaces and collections), PyPI, and alphaXiv. Already digested and not repeated here: [competition-landscape.md](../competition-landscape.md), [504aldo-k3-chain.md](504aldo-k3-chain.md) (through F105), the Phase 1 census, and Oishi1029. Every claim below gives its source and date. Quotes are verbatim. Numbers from other teams are their own measurements and have not been reproduced here.*

---

## 0. Summary

1. **Nobody in the top 5 has disclosed a method.** J2W, marius_binner, Luna, suliman_tadros and mliston have published nothing (forum, GitHub, Hugging Face, blogs). What is new is indirect: the leaders' trajectory, a top-7 team (AndreasHad04) describing its representation, and a third team's (EscAI) cost anatomy of the leaders.
2. **The largest public research drop since 504aldo** is team EscAI's MIT repository `corpaci/ARCwhitebox` (29 Sep; its README cites submission 332675 at rank 42, adjusted 5.0e-9; on the 1 Oct board EscAI is 43=, adjusted 4.8e-9, raw 2.19e-8, util 0.2209, last submitted 30 Sep). It contains three findings that bear directly on our plan:
   - **Covariance-response modes specified and killed.** They wrote out "covariance-response modes" with symmetric-matrix transport `B ← WᵀBW` (the exact idea in our Line C) and killed the expansion on 28 Sep at a preregistered gate (G0). The gate ran on a skeleton with no modes (A = 0) and failed on in-place consumption of the born (2,1) slice; the mode-seed gate (G1) was never run. The evidence that bears on the modes themselves is a follow-up probe: congruent (`WᵀBW`-type) transport of the born content gives |cos| ≤ 0.008 against the true slice.
   - **The error lives in the pre-activation marginals.** In an oracle that replaces pre-activation variance, κ₃ and joint κ₄ at every layer, raw MSE drops about 20× (2.3e-8 to 1.2e-9). Variance alone gives −40%.
   - **504aldo's κ₄ gain is cancellation.** His fitted λ is about 3× the physical value; on an exact base the same channel hurts.
3. **Rules: one substantive clarification and no other changes.**
   - **Fair accounting (16 Sep).** The starter kit and flopscope docs now say metering does not establish eligibility, and the Sponsor may "invalidate, re-score, or disqualify" any submission whose benefit "derives from how computation is accounted", including during prize review. Packing is banned explicitly.
   - **Strassen confirmed permitted (11 Sep).** Strassen–Winograd as flopscope ops is allowed with no recursion limit.
   - **No version changes.** flopscope 0.12.1 and whestbench 0.16.1 are still the latest on PyPI; there has been no library release since 29 Aug (the evaluator reported whestbench 0.16.0 on 29 Sep).
   - **Dates unchanged.** Team freeze is 2 Oct, close 17 Oct, write-ups 24 Oct. **Winners are announced 15 Nov** (stated by a participant on the forum, who cites Rules §6 for the 14-day document window; the organizer's reply does not contradict it).
4. **Grader traps confirmed on 29 Sep:**
   - The smoke test still runs an MLP deeper than 16 layers.
   - `x.shape = ...` works locally but raises on the grader's flopscope-client.
   - On the grader, a 0.34 B chain takes 76–100 s of wall time per MLP, against the 120 s cap.
5. **No new community datasets since 13 Sep, and no relevant new papers.** EscAI independently verified that keenanpepper's `bakev2-mini-supp` weight hashes and `gt_mean` arrays match the public `v2-phase2` mini split bit for bit. They also found that 504aldo's local dumps come from a **different HF snapshot** (`ce926b7d…`) than the current `@v2-phase2` tag (`aa99830f…`).

---

## 1. Rule and grader changes

### 1.1 Fair accounting and packing (16 Sep 2026)

Commits by SP Mohanty on 16 Sep 2026: [whest-starterkit `5eb9aa1`](https://github.com/AIcrowd/whest-starterkit/commit/5eb9aa1455fcb3216af55994bdf25dc242b95797) "docs: align Phase 2 packing guidance with challenge rules", and [flopscope `80e72a92`](https://github.com/AIcrowd/flopscope/commit/80e72a92df67dfa0bd49310fd32c0ff1d357d7f4) (PR [#261](https://github.com/AIcrowd/flopscope/pull/261)). These are the only upstream commits since 10 Sep in either repository. whestbench has none since 29 Aug. Our [starter-kit-operational-facts.md](starter-kit-operational-facts.md) already records the packing ban; the broader wording is what matters:

> "Under the official Rules, a submission's score benefit must derive from its estimation method. The Sponsor may invalidate, re-score, or disqualify submissions whose benefit instead derives from how computation is accounted, including during or after grading and during prize review." (`docs/concepts/allowed-code.md`, starter kit, 16 Sep)

> "**Do not pack several independent values into one machine element to obtain an accounting advantage.** This includes packing booleans into a wider integer for bitwise operations that process those booleans together." (same)

> "Earlier FlopScope documentation described packing as "not banned and not judged" and small sub-32-bit packing gains as "in-bounds". That wording does not describe Phase 2 eligibility." (same)

> "**Metering does not establish competition eligibility.**" … "These are accounting safeguards, not a guarantee that every possible accounting advantage is prevented or permitted" (flopscope `docs/reference/cost-model.md`, 16 Sep)

> "`ctx.summary()` reports the FLOPs metered on your run; it does not determine competition eligibility." (starter kit `docs/reference/flopscope-primer.md`, 16 Sep; previously "`ctx.summary()` on your own run is always ground truth.")

### 1.2 Strassen confirmed permitted (11 Sep 2026)

504aldo, forum topic 18218 post #3, [14 Sep 2026](https://discourse.aicrowd.com/t/everything-we-tried-a-factorized-k-3-cumulant-propagation-estimator-at-0-25-x-b-where-its-flops-go-and-25-measured-dead-ends-team-504aldo-rank-10/18218/3):

> "Small update since the post: the sponsor confirmed (2026-09-11) that Strassen-Winograd as flopscope ops is permitted with no recursion-level limit."

This was already in our landscape note and is restated because §1.1 postdates it. Strassen bills its genuine analytical count, so it is "benefit from the estimation method" in the sense of §1.1. Exploiting billing quirks is not (§1.4).

### 1.3 Smoke test and client parity (29 Sep 2026)

[whestbench #149](https://github.com/AIcrowd/whestbench/issues/149), AndreasHad04, 29 Sep, open:

> "Submission 332734 was rejected before grading with `Error : Smoke test failed`, which is the whole message the submissions API returns. Its only defect was a 16-entry per-layer table indexed by the layer number. … So the evaluator's smoke MLP appears to be deeper than 16 layers (Phase 1 was 256 x 32)."

> "`whest validate --estimator ...` predicts on `sample_mlp(width=4, depth=2)` … so it passes. `whest run --split mini` uses the graded 1024 x 16 shape, so it passes."

Our plan already probes depth 32 and width-4/depth-2; this confirms the trap is still live on 29 Sep with whestbench 0.16.0. Corroboration (ours, from metadata only): the organizers' Hugging Face dataset [`aicrowd/whestbench-smoke-mlp`](https://huggingface.co/datasets/aicrowd/whestbench-smoke-mlp) (last updated 27 Jul, before Phase 2) holds one row in an 8.5 MB parquet file. That matches one 256×32 float32 MLP (256²·32·4 B = 8.4 MB), not 1024×16 (67 MB), so the smoke MLP is probably still the Phase 1 shape. Rows were not read.

[flopscope #267](https://github.com/AIcrowd/flopscope/issues/267), AndreasHad04, 29 Sep:

> "In-process, arrays inherit numpy's writable `.shape`, so `x.shape = new_shape` reshapes in place. The evaluator's `RemoteArray` does not allow it. … **It is not billed in-process.** … **It raises on flopscope-client** (`AttributeError: can't set attribute 'shape'`). Code that runs fine locally fails on the evaluator at that line. Two of our submissions (332205, 332206) failed this way."

### 1.4 Billing-quirk reports filed since 10 Sep (all open)

| issue | date | author | title |
|---|---|---|---|
| [flopscope #263](https://github.com/AIcrowd/flopscope/issues/263) | 17 Sep | srdey2002 | Masked sum retains invalid input symmetry and rejects valid NumPy results |
| [#262](https://github.com/AIcrowd/flopscope/issues/262) | 17 Sep | srdey2002 | `linalg.pinv(..., rtol=None)` silently ignores NumPy's Array API tolerance |
| [#264](https://github.com/AIcrowd/flopscope/issues/264) | 28 Sep | thibaudlepan77-svg | Module-level `random.shuffle` leaves a stale symmetry tag, contractions on the shuffled array bill at the orbit rate |
| [#265](https://github.com/AIcrowd/flopscope/issues/265) | 28 Sep | thibaudlepan77-svg | `choose` with a positional `out` writes into a symmetric-tagged buffer and the tag survives |
| [#266](https://github.com/AIcrowd/flopscope/issues/266) | 28 Sep | thibaudlepan77-svg | `random.permutation` and `choice(replace=False)` bill the pool length while returning a full copy, a 1e6-element copy costs 4 FLOPs |
| [#267](https://github.com/AIcrowd/flopscope/issues/267) | 29 Sep | AndreasHad04 | `.shape` assignment is free in-process, keeps a stale symmetry tag, and raises on flopscope-client |
| [#268](https://github.com/AIcrowd/flopscope/issues/268) | 30 Sep | NathanZane | einsum silently ignores explicit dtype, producing a different numerical result |
| [#269](https://github.com/AIcrowd/flopscope/issues/269) | 30 Sep | NathanZane | Complete QR on tall matrices seems to undercount Q formation |

Three of these (#264, #265, #267) are stale symmetry tags that make a later contraction bill at the cheaper orbit rate. Under §1.1, gaining from them is disqualifiable even after grading. **We must not depend on any symmetry tag we did not create deliberately through `as_symmetric` or the documented Gram/`einsum` aliasing.**

### 1.5 Timeline, prizes and logistics

- **Dates are unchanged.** The latest statement remains the 18197 announcement (23 Aug): team freeze **2 Oct 23:59 UTC**, close **17 Oct 23:59 UTC**, write-ups **24 Oct 23:59 UTC**.
- **Winners are announced on 15 November.** Participant oleksandr_barskyi, [18197 #3, 26 Sep 2026](https://discourse.aicrowd.com/t/phase-2-of-the-arc-white-box-estimation-challenge-is-live/18197/3): "Winners are announced on 15 November and Section 6 of the Official Rules gives 14 days to sign and return the Prize Winner Documents." Only the 14 days is attributed to §6; the date is his statement. The organizer's reply (#4) repeats the 14 days and does not contradict the date, but we could not read the Rules page to confirm it. This date is not in our notes.
- **Payment and eligibility.** mohanty, [18197 #4, 1 Oct 2026, 00:21 UTC](https://discourse.aicrowd.com/t/phase-2-of-the-arc-white-box-estimation-challenge-is-live/18197/4):
  > "Payment: international bank transfer (SWIFT) in USD to an account in the winner's own name."
  > "Winners should expect to provide photo ID, proof of address and bank details, and to return the documents within 14 days."
  > "Sanctions: under Section 3 of the Rules, residents of areas under comprehensive US, EU, UN or UK sanctions can't receive prizes, though they may still take part and keep their rank. Eligibility depends on where you live, not your citizenship."
- **Phase 2 algorithmic-contribution prize: $20,000** (Phase 1: $10,000), at ARC's discretion. Source: [topic 18041](https://discourse.aicrowd.com/t/algorithmic-contribution-prize-guidelines-how-arc-judges-these-prizes-discretion-technical-writeups-llm-usage/18041) (23 Jun; mechanics edited 7 Sep). This is older but missing from our notes. "A PDF technical write-up + exactly one submission ID. Both are mandatory … the submission you intend to reference must already be successfully graded before that phase's submission deadline." ARC is "less interested in methods that rely heavily on sampling, fine-tuned constants, careful performance optimization, and opaque LLM-optimized code."
- **Phase 1 results** have still not been announced as of 1 Oct. The 18197 announcement says "The Phase 1 private re-evaluations are in progress"; no later post exists.
- **Information leak, fixed.** On [9 Sep](https://discourse.aicrowd.com/t/errors-are-leaking-information/18217) marius_binner reported that the global submission feed leaked other teams' error messages ("I can now tell that this guy is attempting to use ridge regression and gradient boosting"). mohanty: "the global submission message is currently leaking some message, and we will fix it right away." Per-submission pages are RBAC-protected.

### 1.6 What the grader UI says about the final score

From [submission #330153](https://www.aicrowd.com/challenges/arc-white-box-estimation-challenge-2026/submissions/330153) (tristan_miano, graded about 27 Sep):

> "Final score · Full test set Graded on all 100 MLPs — the number that ranks you. The private split is sealed until results release." … "These 50 MLPs are graded once and revealed at Phase 2 close."

The same page calls the 50 gated MLPs "The 50 holdout MLPs — a different partition from the public split — that decide the final rank." That wording conflicts with our notes: prize ranking comes only from a fresh private re-run of each team's single designated submission (challenge overview page; 18118, 31 Jul). The UI text is probably loose and describes only the public board (50 public plus 50 gated MLPs, the gated half revealed at close), but the conflict is unresolved. Ask the organizers rather than assume either reading. The same page gives the grader's speed: **mean effective compute 7.57e11 FLOPs (34.4% of B), per-MLP wall 76.3–100.3 s** (median about 78 s), which is about 1e10 FLOP/s. A chain much above 0.4 B is at risk of the 120 s wall cap. At the leaders' 0.11–0.15 B it is not.

---

## 2. The leaderboard and what can be inferred about the leaders

### 2.1 Board on 1 Oct 2026

Source: [leaderboard](https://www.aicrowd.com/challenges/arc-white-box-estimation-challenge-2026/leaderboards), read through Exa on 1 Oct. Top 20 shown; the full 50 rows are in [competition-phase2.md](../competition-phase2.md).

| # | team | adjusted | raw | util | entries | last submission |
|---|---|---|---|---|---|---|
| 1 | J2W | 1.6e-9 | 1.50e-8 | 0.1101 | 364 | 29 Sep 18:17 |
| 2 | marius_binner | 1.7e-9 | **1.14e-8** | 0.1507 | 187 | 30 Sep 21:55 |
| 3 | Luna | 1.8e-9 | 1.41e-8 | 0.1306 | **16** | 30 Sep 15:13 |
| 4 | suliman_tadros | 2.0e-9 | 1.58e-8 | 0.1298 | 80 | 28 Sep |
| 5 | mliston | 2.1e-9 | 1.70e-8 | 0.1258 | 64 | 28 Sep |
| 6 | MeatProxy | 2.4e-9 | 1.52e-8 | 0.1592 | **13** | 30 Sep |
| 7= | a_s6 / emanuel_ruzak / jlacombe / AndreasHad04 / thegamers | 2.5e-9 | 1.55–2.00e-8 | 0.126–0.163 | 12–108 | 21 Sep – 1 Oct |
| 12 | oqaris | 2.6e-9 | 1.98e-8 | 0.1314 | 102 | 28 Sep |
| 13= | Puffi / dipdecir / bryant_le | 2.7e-9 | 1.76–2.00e-8 | 0.133–0.156 | 16–80 | 25–30 Sep |
| 16= | reds / carlos_rodriguez / **lode_dockx** / shiv_m | 2.8e-9 | 1.61–**2.76e-8** | **0.1015**–0.171 | 10–163 | 21–30 Sep |

Readings:

- **J2W is at the 0.1 floor**, so it gains from now on only through raw.
- **marius_binner has the best raw on the board (1.14e-8)** at 0.15 B; raw about 1.1e-8 is demonstrably reachable.
- **Luna and MeatProxy reached the top 6 with 16 and 13 entries**: late or very efficient entrants.
- **lode_dockx is a pure floor rider**: raw 2.76e-8 at 0.1015 B. That is about 25% worse raw than the public chain (about 2.2e-8 at 0.25 B), at about 40% of its cost.
- **Ranks 35–50 sit at raw 2.06–2.23e-8 and util 0.20–0.24**, the signature of 504aldo V29 clones. jamesrahenry, [18218 #2, 13 Sep](https://discourse.aicrowd.com/t/everything-we-tried-a-factorized-k-3-cumulant-propagation-estimator-at-0-25-x-b-where-its-flops-go-and-25-measured-dead-ends-team-504aldo-rank-10/18218/2): "What does it feel like to have 13 different exact submissions under other people's accounts?"

### 2.2 J2W's trajectory

| date | adjusted | raw | util | source |
|---|---|---|---|---|
| 1 Sep | ≈5e-9 | 2.27e-8 | 0.22 | 504aldo findings log F59 (our digest) |
| 5 Sep | — | — | 0.213 → 0.186 | F77 |
| 6 Sep | — | ≈2.0e-8 (1.97–2.02e-8 for the 0.15–0.16 B cluster that includes J2W) | 0.164 | F86 |
| ≈27 Sep | 1.8e-9 | 1.66e-8 | 0.1089 | EscAI `research/tenfold_push/README.md` |
| 28 Sep 19:47 UTC | 1.7e-9 | 1.62e-8 | 0.1082 | EscAI `research/contraction_push/leaderboard_update.json` |
| 1 Oct | 1.6e-9 | 1.50e-8 | 0.1101 | live board |

Cost fell 2× in September and then stopped at the floor. Raw improved over the whole month: from 2.27e-8 (1 Sep) to about 2.0e-8 (6 Sep), then 1.66e-8 (about 27 Sep), then 1.50e-8 (1 Oct). The last step is about 10% in four days.

### 2.3 J2W's identity, Phase 1 signal (circumstantial)

- **Team and account.** The team page lists one member, organizer **joe_wanza** ([team page](https://www.aicrowd.com/challenges/arc-white-box-estimation-challenge-2026/teams/J2W)).
- **Phase 1 rank.** Our Oishi digest records joe_wanza at **#2 in Phase 1** on 7 Aug (1.00e-9).
- **Probable GitHub account: `jwanza`.** On 17 Jul jwanza filed [flopscope #147](https://github.com/AIcrowd/flopscope/issues/147) (and its duplicate #149), proposing a fused `gather_segment_sum` op for "Lookup-table-based estimators [that] compute `out[r, o] = sum_s table[flat_index[r, s], o]`", noting "The common quantized lookup-table case is an `int16` table accumulated in `int32`". Mohanty closed it on 15 Aug as out of scope: "The surface we are holding still for now is the set of ops that already exist."
- **What this suggests.** If jwanza is joe_wanza, J2W's Phase 1 method used quantized lookup tables with gather/segment sums.
- **Since then:**
  - Flopscope prices computed-index gathers at weight 4.0 ("No free-gather discount", cost-model doc).
  - The 16 Sep rule bans packing, though it says "Choosing a lower-precision dtype or quantizing values is not, by itself, this packing technique", so int16 tables are not banned as such.
  - J2W's Phase 2 cost history (0.22 → 0.11 B) looks like a chain being compressed, not a lookup-table trick.
- **Nothing public ties J2W's Phase 2 method to lookup tables.**

### 2.4 AndreasHad04 (rank 7) disclosed its representation (22 Sep 2026)

[18218 #4](https://discourse.aicrowd.com/t/everything-we-tried-a-factorized-k-3-cumulant-propagation-estimator-at-0-25-x-b-where-its-flops-go-and-25-measured-dead-ends-team-504aldo-rank-10/18218/4):

> "We had a kappa4 to kappa3 FEED arm sitting unpriced for a week, and F68 gave us something specific enough to test before building anything, so we tested it first."

> "Our chain already materialises the (2,2) slice of kappa4 of the post-activation, which is your `G`, and both candidate covariances."

> "8 networks, our shipped configuration, cap 3776."

> "Your F86 says 95% of your bill is the K3 source machinery and ours is a different factorisation, so the law may hold exactly in your chain and we would not see it here."

Their measurement: the post-activation (2,2) κ₄ slice is regressed on C_pre off the diagonal, through the origin. R² is 0.998 at layer 1, falling monotonically to 0.045 at layer 15, with λ falling from 3.2e-1 to 1.6e-3. Permuted control R² ≤ 4.7e-7.

504aldo's reply ([#5, 22 Sep](https://discourse.aicrowd.com/t/everything-we-tried-a-factorized-k-3-cumulant-propagation-estimator-at-0-25-x-b-where-its-flops-go-and-25-measured-dead-ends-team-504aldo-rank-10/18218/5)) and commit [`1a4083f`](https://github.com/504aldo/whest-p2-cumulant-k3/commit/1a4083f89184d6b1cc55d2e55a5a561d2ae50a8d) (F105):

> "F68's G is not the (2,2) slice of kappa4. It is the r=1 matrix core of the augmented chain. Its off-diagonal comes from the (2,1,1) hub contraction and Sym(K31); the (2,2) slice enters G only through its row sums, on the diagonal"

> "From V22 on, the code uses the pre-activation law, `G_pre_off = lambda * C_pre_off`. On 8 MLPs its R^2 is 0.48-0.52 at layer 0, >= 0.95 from layer 3 and 0.985-0.99 at depth."

> "replacing G_off by lambda_l C_off in the chain gives 2.23e-8, vs 2.04e-8 with the dense channel and 4.20e-8 without it."

What this tells us:

- **A second top-10 representation exists.** It is "a different factorisation" from 504aldo's per-source legs and materialises the post-activation κ₄(2,2) slice (an n×n object) and both covariances. Its κ₄→κ₃ feed arm was "sitting unpriced" when they posted, so it is a candidate, not a confirmed part of the shipped chain. The measurement "decided a build" for them, without saying which way. It reaches raw 1.55e-8 at 0.163 B, better raw than any leg chain (504aldo's float64 ceiling is about 2.07e-8). "Cap 3776" is an unexplained configuration parameter.
- **The (2,2) slice is not C-shaped at depth**, so whatever AndreasHad04 carries there is real off-C content.

### 2.5 Cost anatomy of the leaders, from EscAI (29 Sep 2026)

`corpaci/ARCwhitebox`, [RUNGS.md](https://github.com/corpaci/ARCwhitebox/blob/main/RUNGS.md):

> "**The leaders' cost signature** (0.1496B) equals our young-tier-only bill to 0.3% — they carry old-source content at near-zero marginal cost by a construction that the entire measured candidate space of this repo does not contain."

`research/astra_oldtier.md` §4:

> "For a young-only leader at 0.131-0.15B to reach 1.8e-9 needs raw ≤ 1.8e-9/0.131 = **1.37e-8**, i.e. -34% vs the converged public raw (2.1e-8) … **Conclusion: the 0.1496B-cluster leaders carry old content at near-zero marginal bill.** … (iii) the #1 entry (1.8e-9) could instead be a 0.1-floor rider with raw ~1.8e-8 (= 1.8e-9/0.1), only 14% better than public raw at ≤ 0.1B — not excludable from the anatomy match"

504aldo, [18218 #3, 14 Sep](https://discourse.aicrowd.com/t/everything-we-tried-a-factorized-k-3-cumulant-propagation-estimator-at-0-25-x-b-where-its-flops-go-and-25-measured-dead-ends-team-504aldo-rank-10/18218/3):

> "From the leaderboard I see someone trimmed 7-22% of the chain's cost (4.2e-9 to 5.2e-9) without losing accuracy; cost anatomy (fig 2 in the post) did not find that for me."

**Inference (ours).** J2W now sits at **0.110 B, below the 0.131–0.15 B young-tier-only bill** of the public representation. So J2W does not run the public chain as published with a cheap old tier. Either its young tier is engineered differently (EscAI's measured young-block floor is 0.055–0.07 B, so a leaner young tier fits under 0.11 B), or it runs a different representation altogether. marius_binner (0.151 B, raw 1.14e-8) sits on the 0.15 B cluster line, but with raw 46% better than the public chain.

---

## 3. New code: team EscAI, `corpaci/ARCwhitebox` (MIT, pushed 29 Sep 2026)

Source: <https://github.com/corpaci/ARCwhitebox>, two commits, both 29 Sep (`c0a1c6c` "Public research corpus: WHEST Phase 2 cumulant estimator + falsification ledger"; `144f2d6`). Author Luiza Corpaci, team EscAI. Graded **submission 332675**: raw 2.19e-8, util 0.2308, adjusted 5.0e-9, zero failures, rank 42 (45 tied on 28 Sep). The 1 Oct board shows a later EscAI entry (30 Sep 07:22) at 43=, adjusted 4.8e-9, util 0.2209, same raw. 502 files.

**Shipped estimator.** 504aldo V29 is embedded verbatim with attribution, plus their own Strassen join routing and a "rung-7 squeeze" (leaf-16, stage-off, old-source rank 352): 0.2430 → 0.2304 B at raw 2.27e-8.

### 3.1 Covariance-response modes and symmetric-matrix transport: specified, killed at the first gate (28 Sep)

[`research/astra_responsemode_spec.md`](https://github.com/corpaci/ARCwhitebox/blob/main/research/astra_responsemode_spec.md) defines the state:

> "B_a, a = 1..A | (n,n) sym, **tagged** | covariance-response modes: B_a(l) = T_{l:s_a} S_a T_{l:s_a}ᵀ, the two-sided transport of a symmetric seed S_a frozen at seed layer s_a. A = 4: seeds at s_a = 1, 5, 9, 13, S_a = C(s_a)"

with update `B_a <- Wᵀ B_a W` ("tagged sandwich 1.5 u each") and `B_a <- (w1 w1ᵀ) o B_a` at the ReLU. Its own up-front verdict:

> "**The bill target is reachable; the raw target is not — under every closure my results support.** … the full expansion below prices at **0.13–0.15B** … Raw: … its omission projects via the eps^2 law … to **raw ~0.8–1.7e-7**"

> "What the leaders' tier implies is NOT this expansion with derived closures: either fitted closures + a weaker error coupling (= G0 passing for them), or constructions absent from all public evidence."

The kill ([`research/DEADENDS.md`](https://github.com/corpaci/ARCwhitebox/blob/main/research/DEADENDS.md), "Response-mode expansion (rebuild thread, killed at G0 2026-09-28)"). The harness `port_rm.py` is a skeleton with "A=0 modes, analytic born slices consumed in place". The B_a modes were therefore never run, and the mode-seed gate G1 was never reached. What was measured is the born-slice consumption, plus a post-freeze transport probe (last quote):

> "P2 as written FAILS at object level, not tolerance level: born D21_hat vs the lean chain's window-1 (age-1) S21 content has cosine 0.005-0.018"

> "Mechanism: born D21/D3 are POST-activation cumulants of the current ReLU; the chain consumes PRE-activation cumulants, i.e. the born content of the PREVIOUS layer after one linear transport. The (2,1) slice does not ride C congruently through W (one-sided hub contraction, T'_iij = sum_a W_ai h_a (CW)_ai (CW)_aj): transporting it needs one dense product family per layer"

> "Measured damage (nets 60-65, f64): born-off skeleton 3.66e-6 (covprop class, sane); born ON 2.42e-4 (66x WORSE)"

> "G0 … extra-MSE = A*delta + K*delta^2 with A ~ 5.1-5.5e-5 … K ~ 5.95-6.11e-5 … coupling is ~14x STRONGER than the lean chain's 4.2e-6 law."

> "born tensor of layer l-1 with ONE correct dense transport vs exact age-1 S21 at l: cos 0.980-0.995 … congruent (census-compatible) transport control: |cos| <= 0.008 … Born content is right; the expansion cannot transport it; the transport that fixes it (~6 dense n^3-class products/layer) is the measured young-block floor. The circle closes."

### 3.2 Where the remaining error lives (oracle interventions, 27–28 Sep)

[`research/tenfold_push/README.md`](https://github.com/corpaci/ARCwhitebox/blob/main/research/tenfold_push/README.md). Reference moments come from keenanpepper's `whest-p2-bakev2-mini-supp` (N = 1e8); hashes verified against public mini weights. Base: their frozen V29-derived estimator, NumPy backend.

Readout truncation is not the bottleneck:

| scalar readout formula fed with reference moments | final cross-replica bias MSE |
|---|---|
| Gaussian | 2.03e-7 |
| existing K3/K4 mean formula | **1.17e-10** |
| Edgeworth to CLT weight 2 / 3 / 4 | 6.5e-11 / 5.0e-12 / 3.1e-13 |

Interventions at every layer (the mean propagates naturally; selected fields are replaced by reference values):

| intervention | nets | mean raw MSE |
|---|---|---|
| baseline | 0–2 | 2.450e-8 |
| pre-activation variance | 0–2 | **1.471e-8** (−40%) |
| pre-activation variance + κ₃ | 0–2 | 8.35e-9 |
| pre-activation variance + joint κ₄ | 0–2 | 8.92e-9 |
| post-activation variance only | 0–2 | 2.527e-8 (no gain) |
| pre-act variance + κ₃ + joint κ₄, exact Gaussian layer 0 | 0–7 | **1.165e-9** (ratio 0.0497, 95% CI [0.046, 0.055], every net improves) |

> "correcting three pre-activation moment channels throughout the network reduces eight-network raw MSE by about 20×. The tested cheap approximations do not reproduce that reduction."

The same 20× is reached through a **linear deferred-mean path**, which leaves the trajectory fixed and only propagates the readout's local correction through Wᵀ ([`research/observable_push/README.md`](https://github.com/corpaci/ARCwhitebox/blob/main/research/observable_push/README.md)): "Eight-network raw MSE falls from **2.34422e-8 to 1.17282e-9**, about **20× lower**." Deployable versions fail. On held-out nets 4–7, their existing moment fit applied through deferred transport gains 3.2% (ratio 0.9684), and the 40-feature propagated-response ridge gains 0.7% (0.9928); the 40-feature span explains only 8.77% of baseline error.

Covariance error structure ([`research/contraction_push/README.md`](https://github.com/corpaci/ARCwhitebox/blob/main/research/contraction_push/README.md)), layer 14:

- The covariance sketch has relative RMS error ≈ 0.70%.
- Post-K22 error is 26.8%.
- An oracle low-rank (4–32) symmetric fit of the covariance residual leaves 60–67% of the held-out sketch error and about 80% of next-layer variance error.

> "some late covariance error is coherent and generalizes across sketch directions. Even reference-assisted low-rank recovery does not remove most of the measured variance error."

### 3.3 The κ₄ λ-law is cancellation, and a dataset-snapshot mismatch

[`research/astra_k4_overshoot.md`](https://github.com/corpaci/ARCwhitebox/blob/main/research/astra_k4_overshoot.md), 25 Sep:

> "The chain's kappa4 core carries **~3x the physical C-aligned content and mostly non-physical structure at depth** — a truncation artifact of the least-converged tracked object"

> "the fitted lambda (7-8e-3) is fitted to the chain artifact; the physical value is 2.1-2.6e-3 ~= **0.3x the fitted table**"

> "504aldo's win living in cancellation: their table (peak 1.1e-2) is 3-5x physical; an oversized C-shaped correction acts as an error-canceller against their pruned base's C-shaped truncation bias"

RUNGS rung 3 records their attempt to transfer it: the λ refit on their own data, holdout raw 4.036 → 3.484e-8 (−13.7%), "**FAIL — shipped dormant.** … the dense augmented chain is *worse* than the plain chain here".

On data provenance:

> "504aldo's dumps are built from snapshot `ce926b7d1ced65faa668c01b94a43f8c4f43783d` … Our mini_p2.npz was extracted from the currently-pinned @v2-phase2 snapshot `aa99830fdc09fad15407b10e8e3459d3e18bba0a` … Different snapshots => plausibly different weight draws (a seed-protocol version bump changes every MLP) and possibly different GT quality."

**Corroboration (ours).** A committed HF cache in [`barnobarno666/ARC-White-Box-Estimation-2026`](https://github.com/barnobarno666/ARC-White-Box-Estimation-2026) (pushed 5 Sep) holds `refs/v2-phase2 = aa99830fdc09fad15407b10e8e3459d3e18bba0a`. The tag pointed at `aa99830f` by 5 Sep at the latest (the ref file was committed in `0887a59`, 5 Sep 16:33 UTC), so 504aldo's `ce926b7d` is a different snapshot, presumably the pre-30-Aug state (the dataset's last-modified date is 30 Aug). Whether the MLPs differ is unresolved. EscAI's suggested check is to compare `mlp_seed` of row 0 against `6319981554997072999`.

### 3.4 Other measured negatives and positives worth keeping

- **eps² law, measured independently** (RUNGS rungs 8–9): "extra MSE ≈ 4.2e-6·eps² … every compressive construction we measured … lands at eps 0.15–8 against a bar of 0.03–0.05". Cup riding `Sym(v⊗C)` gave 6.8× worse; frozen-factor truncation +550…3000%; column sketches worse than dropping.
- **Young-block cost floor 0.055–0.07 B** (`astra_youngblock.md`), from "no-free-sandwich" and "materialization" lemmas. The audit corrects one claim: "computing the symmetrization of `X Y^T` requires both dot products per off-diagonal output. Symmetry alone does not halve the leading FLOPs of one GEMM." (`top3_push/README.md`)
- **Gaussian cross-moment terms to Hermite order 4 (Mehler-4):** 0.4% gain on 3 nets, not promoted.
- **Gaussian scale-mixture κ₄ closure:** about 10% worse. **Empirical state regression:** 17.5% worse in closed loop. **Learned transported memory:** flat (1.0016). **Spectral rank allocation and SVD residual predictors:** worse.
- **Keeping all repeated κ₄ entries ((4), (3,1), (2,2), (2,1,1)) in a compressed fourth-order state is positive:** −35% raw at width 128 (3 nets) and −64% at width 192 (1 net). This uses the reference algorithm at small width, not the competition chain (`observable_push`).
- **An exact O(n³) shortcut** for the transported diagonal of the covariance path/star births (`nonlinear_push/direct_diagonal.py`): "Official flopscope measures **4,303,356,928 float32 FLOPs at n=1024** for the complete component; applying it on all 15 transitions would use **2.935% of B**." Not integrated.
- **Local resource caution:** with their frozen baseline under the official subprocess runner on a Mac, "All failures exceeded the **0.4-second residual-time cap**" (2.3 / 3.1 / 0.8 s), while the online entry has zero failures. Residual timing is very machine-dependent, which agrees with our plan's 2× margin rule.
- **Final pre-activation covariance spectrum of the estimator (8 nets):** stable rank 2.57, participation rank 5.9, 25 directions for 90% energy, 79 for 99% (`spectral_push`).

---

## 4. Other public code since 10 Sep (low value)

| repo | pushed | what |
|---|---|---|
| [anhminhzui-dev/arc-whest](https://github.com/anhminhzui-dev/arc-whest) | 10–11 Sep | Starter-kit fork; full-covariance Gaussian propagation (submission 330458, raw ≈4e-6). Notes the ARC torch K=3 reference reaches 1.18e-7 but is not grader-ready. |
| [barnobarno666/ARC-White-Box-Estimation-2026](https://github.com/barnobarno666/ARC-White-Box-Estimation-2026) | 5–14 Sep | Abandoned ("Ran out of tokens and time"); useful only for the HF ref in §3.3. |
| [MurtuzaShaikh26/ARC-White-Box-Challenge-26](https://github.com/MurtuzaShaikh26/ARC-White-Box-Challenge-26) | 8 Sep | "Scaling law for the bias constants: width holds, depth does not". Calibration-constant study at mean-propagation level. |
| [504aldo/whest-p2-cumulant-k3](https://github.com/504aldo/whest-p2-cumulant-k3) | 22 Sep | Only the F105 correction (§2.4); already in our digest. |

GitHub searches for repositories by leaderboard names (J2W/jwanza, marius_binner, Luna/dogus_ozel, suliman_tadros, mliston, MeatProxy, AndreasHad04, thegamers, oqaris, Puffi) found nothing new. jwanza's public repositories are unrelated (Tsetlin machines, gem5).

---

## 5. Datasets and literature

- **Hugging Face:** no new whest-related datasets, models, spaces or collections after 13 Sep 2026 (searches: whest, whestbench, arc whitebox, white-box estimation, relu mlp cumulant, mlp moments; keenanpepper namespace listed). The newest remain keenanpepper's `whest-p2-bakev2-{d8b-sketch-g00, bench-supp, mini-supp}` (13 Sep) and `arc-whestbench-p2-full1000-N1e9` (updated 13 Sep), all in our landscape note.
- **Verification from EscAI** (`tenfold_push/mini_data_verification.json`): "All eight stored weight SHA-256 values match the actual local float32 weights; all stored `gt_mean` arrays equal the local 1e9-sample truth bit-for-bit" for mini networks 0–7 in `bakev2-mini-supp`. Regenerating d8b weights from the documented torch recipe on an ARM Mac **did not** reproduce the stored hash (`top3_push/README.md`).
- **alphaXiv, papers after May 2026: nothing that changes the method.**

  | paper | date | content | why it does not apply |
  |---|---|---|---|
  | [2608.20483](https://www.alphaxiv.org/abs/2608.20483) | 20 Aug | Leaky-ReLU local linearization for uncertainty propagation | exact only within one activation pattern |
  | [2609.27244](https://www.alphaxiv.org/abs/2609.27244) | 23 Sep | full-covariance Gaussian moment smoothing for Bayesian neural networks | Gaussian closure only |
  | [2605.24072](https://www.alphaxiv.org/abs/2605.24072) | 22 May | Edgeworth expansions of finite-width network outputs over the weight distribution | annealed; by 504aldo's triage rule annealed results have κ₃ = 0 for our quenched problem |
  | [2609.32695](https://www.alphaxiv.org/abs/2609.32695) | 26 Sep | Kac–Rice affine geometry of Gaussian ReLU networks | not an estimation method |

- **PyPI:** flopscope 0.12.1 and whestbench 0.16.1 are the latest (checked 1 Oct). Neither library has had a release since 29 Aug.
- **Not reachable:** Discord (`discord.gg/4gyQvzWPJ`); the Rules page (rendered client-side, empty through Exa). Topic IDs 18220–18227 return nothing through Exa; the category listing on 1 Oct shows no new topics after 18219 (17 Sep).

---

## 6. What this changes for us

1. **Line C cannot be "covariance-response modes" as written.**
   - **Why.** Our plan's revision log names "symmetric n×n objects (covariance-response modes)" as the leaders' probable representation. EscAI wrote exactly that state (`B_a ← WᵀB_aW`, four C-snapshot seeds), priced it at 0.13–0.15 B, and killed the expansion on 28 Sep before any mode was run (the gate harness had A = 0). The born (2,1) slice is a post-activation object and does not transport congruently through W. In-place consumption is 66× worse, and congruent transport of the born content gives |cos| ≤ 0.008 against the true slice. The error coupling is 14× stronger than the chain's. The case against modes is therefore structural, from the transport probe, not a measured mode run.
   - **What would revive it.** A transport rule for the (2,1) slice that uses only elementwise and congruent operations. Nobody has one, and the correct transport costs about 6 dense products per layer, which is the young-block floor.
   - **Action.** Do not rebuild that spec. Keep the gate-statistics branch of Line C (P(z_i>0, z_j>0), gate–activation cross-moments), which nobody has tested publicly. Screen anything new with their `research/eps_prescreen.py` and `port_rm.py` (MIT) before building it.
2. **Point the 1,000-network oracle (Line A step 2) at the pre-activation marginals, not only at D21.** EscAI's 8-network oracle says:
   - Most of the remaining raw error is in pre-activation variance, κ₃ and κ₄ marginals: variance alone gives −40%, all three give 20×.
   - The same 20× comes through a cheap linear deferred-mean path.
   - The readout is not the bottleneck (formula bias 1.2e-10).

   **Action:** replicate their intervention table on 1,000 networks with keenanpepper's N = 1e9 full-split moments. Their targets were noisy (1e8 samples) and covered 8 nets; the grid gives about 125× the networks. Then decompose the variance error into its sources (κ₃ feedback into C, κ₄, the Mehler covariance map). The variance-only oracle (−40%, 3 nets, with layer-0 variance replaced by a noisy reference) is the largest single-channel effect measured. EscAI caution that "noisy/interacting interventions are not a unique causal decomposition" and that replacing post-activation variance alone gains nothing. Nobody has a deployable way to compute that variance, and its cost is unknown.
3. **Treat every λ-table as base-specific.**
   - The public κ₄ "regeneration" gain is cancellation against the pruned base's bias (λ ≈ 3× physical; physical λ ≈ 2.1–2.6e-3, flat in depth).
   - It was fitted on HF snapshot `ce926b7d`, not the current `@v2-phase2` (`aa99830f`).
   - **Action:** pin our staged dataset by commit hash and record `mlp_seed` of row 0. If we change the base chain, refit λ (or drop the channel) on that base, and test on held-out nets before trusting the 2× from 504aldo's V17.
4. **Lever (b) of our §3.1 oracle ladder (a better-than-rank-1 carrier for the (2,1,1) κ₄ slice) now has outside support.**
   - EscAI's small-width runs (reference algorithm, width 128, 3 nets) found that keeping all repeated κ₄ entries cuts raw by 35%. The diagonal alone gives 24.5% of that, so the (2,1,1)-specific increment is about 14% relative. At width 192 the cut is 64%, on 1 net.
   - AndreasHad04 (rank 7, raw 1.55e-8) materialises the post-activation κ₄(2,2) slice in a non-leg factorization (its κ₄→κ₃ feed was unbuilt on 22 Sep) and finds that slice far from C-shaped at depth.

   Two independent teams found real off-C fourth-order content. Neither isolated the (2,1,1) slice at competition width, so this is support, not confirmation. Our diagram identification (first-order gate diagrams including [κ₄(2,1,1)]) is still unpublished and is a candidate write-up contribution; cite EscAI and AndreasHad04.
5. **The money line moved and the leader is at the cost floor.**
   - J2W: raw 1.50e-8 at 0.110 B (adjusted 1.6e-9). Its raw improved 34% over September (2.27e-8 on 1 Sep): about 10% since about 27 Sep (1.66e-8) and 7% since 28 Sep (1.62e-8).
   - Beating J2W on the private re-run likely needs **raw ≲ 1.3–1.4e-8 at ≤ 0.10 B** by 17 Oct, or **≲ 0.9–1.0e-8 at 0.15 B**.
   - Since marius has 1.14e-8 at 0.15 B, a cost cut alone would put him first. J2W does not run the public chain's young tier as published (0.11 B is below its 0.131–0.15 B bill).
   - Our plan's §0 targets (raw ≤ 1.6e-8 at ≤ 0.10 B) are therefore about 10–20% too loose. **Revise them to raw ≤ 1.3e-8 at ≤ 0.10 B.**
6. **Rules and grader hygiene to add to §5 (submission discipline):**
   - (a) Use no symmetry tag we did not create through `as_symmetric` or the documented Gram/`einsum` aliasing. Stale tags after `shuffle`, `choose`, `permutation` or `.shape` (flopscope #264–#267) are disqualifiable under the 16 Sep fair-accounting text, "including … during prize review".
   - (b) Never assign to `.shape`; it fails on the grader.
   - (c) Keep the depth-32 smoke probe; it is still live on 29 Sep.
   - (d) Grader wall time is about 1e10 FLOP/s (0.34 B ≈ 78–100 s per MLP), so stay ≤ 0.35 B for wall-cap margin. This is irrelevant at the target cost.
   - (e) Strassen–Winograd is explicitly permitted.
7. **Dates and admin:**
   - Team freeze **2 Oct 23:59 UTC** (tomorrow).
   - Designate the final submission before **17 Oct 23:59 UTC**; it must already be graded for the write-up to count.
   - Write-up PDF and one submission ID by **24 Oct**.
   - Winners **15 Nov**; documents within 14 days (photo ID, proof of address, bank details; SWIFT, USD, account in the winner's own name).
8. **The algorithmic-prize field ($20k) is crowded with public, MIT, falsification-ledger corpora** (504aldo, 10 Sep; EscAI, 29 Sep), both LLM-assisted with disclosure.
   - Our write-up must cite both and claim only what is beyond them: the 1,000-network oracle tables, the diagram identification of the all-distinct κ₃, and any working representation.
   - ARC states it is "less interested in … fine-tuned constants", so a fitted λ-table or ridge correction will not carry a prize write-up.

---

## Sources (all accessed 2026-10-01)

| source | date of content |
|---|---|
| Forum category listing, <https://discourse.aicrowd.com/c/white-box-estimation-challenge-2026/2991/l/latest> | state on 1 Oct |
| Phase 2 announcement and replies, <https://discourse.aicrowd.com/t/phase-2-of-the-arc-white-box-estimation-challenge-is-live/18197> (#3, #4) | 23 Aug; 26 Sep; 1 Oct |
| 504aldo thread replies, <https://discourse.aicrowd.com/t/everything-we-tried-a-factorized-k-3-cumulant-propagation-estimator-at-0-25-x-b-where-its-flops-go-and-25-measured-dead-ends-team-504aldo-rank-10/18218> (#2–#5) | 13, 14, 22, 22 Sep |
| "Errors are leaking information", <https://discourse.aicrowd.com/t/errors-are-leaking-information/18217> | 9 Sep |
| Algorithmic prize guidelines, <https://discourse.aicrowd.com/t/algorithmic-contribution-prize-guidelines-how-arc-judges-these-prizes-discretion-technical-writeups-llm-usage/18041> | 23 Jun; edited 7 Sep |
| Leaderboard, <https://www.aicrowd.com/challenges/arc-white-box-estimation-challenge-2026/leaderboards> | 1 Oct |
| Submission #330153, <https://www.aicrowd.com/challenges/arc-white-box-estimation-challenge-2026/submissions/330153> | ≈27 Sep |
| Team pages J2W, Luna, `…/teams/J2W`, `…/teams/Luna` | 1 Oct |
| whest-starterkit commit `5eb9aa1`; flopscope commit `80e72a92` / PR #261 | 16 Sep |
| flopscope issues #262–#269, <https://github.com/AIcrowd/flopscope/issues> | 17–30 Sep |
| flopscope issues #147, #149 (jwanza) | 17 Jul; closed 15–16 Aug |
| whestbench issue #149, <https://github.com/AIcrowd/whestbench/issues/149> | 29 Sep |
| `corpaci/ARCwhitebox` (EscAI), <https://github.com/corpaci/ARCwhitebox>: README, RUNGS.md, PHASE2_REVIEW.md, research/DEADENDS.md, astra_responsemode_spec.md, astra_k4_overshoot.md, astra_oldtier.md, research/{top3,tenfold,contraction,observable,memory,nonlinear,spectral}_push/README.md | 24–29 Sep (pushed 29 Sep) |
| 504aldo commit `1a4083f` (F105) | 22 Sep |
| `barnobarno666/ARC-White-Box-Estimation-2026` (HF ref file) | 5 Sep |
| `anhminhzui-dev/arc-whest`; `MurtuzaShaikh26/ARC-White-Box-Challenge-26` | 8–11 Sep |
| Hugging Face MCP searches; keenanpepper dataset listing | 1 Oct (newest item 13 Sep) |
| `aicrowd/whestbench-smoke-mlp` metadata (1 row, 8.5 MB parquet), <https://huggingface.co/datasets/aicrowd/whestbench-smoke-mlp> | updated 27 Jul |
| alphaXiv `discover_papers` (two queries, after 15 May 2026) | 1 Oct |
| PyPI `pip index versions flopscope / whestbench` | 1 Oct |
