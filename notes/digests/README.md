# Digests

Machine-generated extractions (made in this session by reader agents from the primary sources named in each file), kept because the sources are large and the facts are needed repeatedly. They are *not* primary sources: every number should be checked against the original before it is relied on in a submission or a write-up.

| file | source | what it holds |
|---|---|---|
| `paper-2605.05179-cumulant-propagation.md` | Wu, Lecomte, Winer, Robinson, Hilton, Christiano, arXiv:2605.05179 (full text) | MLP definition, Hermite conventions, the diagram summation formula, the basic / augmented / factorised algorithms, FLOP polynomials (Table 1), error theorems, experiments, baselines |
| `starter-kit-operational-facts.md` | AIcrowd/whest-starterkit (README, CHANGELOG, how-to, troubleshooting, reference, examples) | round constants, CLI flags, packaging and submission, example estimators and their FLOP ledgers, every cost figure in the docs, performance tips, score-report fields, common errors, algorithm ideas |
| `harness-and-flopscope-internals.md` | AIcrowd/whestbench and AIcrowd/flopscope sources (0.16.1 / 0.12.1) | scoring pipeline verbatim, runners and the subprocess protocol, dataset loader forms, bake internals, estimator contract, the flopscope cost model (weights, dtype rates, formulas, measured costs at n = 1024), client/server architecture and the hardened sandbox |
| `504aldo-k3-chain.md` | github.com/504aldo/whest-p2-cumulant-k3 (MIT; estimators, docs, findings log F1–F105, lean reference chain) | the published Phase 2 K=3 chain layer by layer with code, the FLOP ledger, the ablation ladder, every finding with numbers, the dead-end catalogue, the ground-truth harness, engineering facts |
| `forum-phase1-census-and-solution-threads.txt` | forum topics 18193 and 18157 (plain text) | the Phase 1 technique census with its 18 Aug update, and the Phase 1 solution / Phase 2 ideas thread |
| `oishi1029-phase1-whitened-mc.md` | github.com/Oishi1029/arc-whestbench-2026 (entrant bin_yong_bong, Phase 1: README, WRITEUP, RESEARCH, LITERATURE, probes, estimators) | the whitened/antithetic Monte Carlo estimator ladder with every measured gain, the chaos-spectrum and cubature probes (what exactness in degree D buys), the oracle-anchoring ladder by layer, the closed sampling routes, the Phase 1 leaderboard reading |
