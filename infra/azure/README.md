# Azure experiment grid for the WhestBench Phase 2 work

Purpose: run thousands of estimator evaluations and offline computations in parallel, with
the official harness, so that an experiment on the full public split (1,000 MLPs) or on
freshly baked MLPs takes minutes, not days, and so that offline precomputation (shipped as
data files, which the rules permit) is limited by ideas rather than by machines.

Design principles

- **Cost is not a constraint; quota is.** Batch-service-mode pools are limited by each Batch
  account's core quota per VM family (not by the subscription's regional vCPU quota), and new
  accounts often start at 0 for some families. The grid therefore shards across many regions
  (one Batch account and one autoscaling pool per region) and the driver round-robins tasks
  across them. `00_quotas.sh` prints the headroom; request increases (portal > Batch account >
  Quotas) where the table is small.
- **Embarrassingly parallel by construction.** One task = (one estimator bundle, one
  dataset shard of ~16 MLPs). FLOP accounting is hardware-independent, so results are
  identical to the grader's for `flops_used` and MSE; only wall-clock figures differ, and
  the runner records residual time separately so the 0.4 s cap can be checked on the
  grader's 2-vCPU footprint by a final `whest run --runner docker` before submitting.
- **Data lives in blob, in-region.** The 70 GB `full` split and 7 GB `mini` split are staged
  once from Hugging Face into the `dataset` container (`03_stage_dataset.sh`, runs inside
  Azure because this container's egress policy blocks huggingface.co); tasks download only
  their shard. Our own bakes (fresh seeds, `whest dataset bake`) go to the same container.
- **Everything is a file.** Estimator bundles are tarballs in the `submissions` container;
  results are one JSON per task in `results`; `grid.py collect` turns them into a table.
- **Official code path.** Tasks run `whest run --runner subprocess` on a local parquet
  shard with the graded defaults (budget 2^41, 0.4 s residual cap, 5 s setup cap) and a
  relaxed wall-clock cap, inside the pinned harness image. Nothing is re-implemented.

Layout

| file | role |
|---|---|
| `env.sh` | names, regions, VM size, slots per node; source it first |
| `00_quotas.sh` | vCPU quota table per region/family, Batch account quotas |
| `01_bootstrap.sh` | resource group, storage + containers, container registry, Batch accounts (one per region) |
| `02_build_image.sh` | cloud build of `runner/Dockerfile` into the registry |
| `03_stage_dataset.sh` | one Batch task that copies the HF dataset into blob |
| `04_pools.sh` | autoscaling container pools, one per region, scale-to-zero when idle |
| `05_gpu_pool.sh` | one autoscaling GPU pool (A100 by default) for bakes and moment atlases |
| `grid.py` | `submit` a grid (estimators × shards), `status`, `collect` results into a table; `bake`, `atlas` on the GPU pool |
| `runner/run_shard.py` | in-container: fetch shard + bundle, run the harness, upload JSON |
| `runner/stage_dataset.py` | in-container: HF → blob copy |
| `runner/bake_shard.py`, `runner/bake_moments.py` | in-container (GPU): one bake slice; one moment atlas |
| `tests/run_offline.sh` | offline dry run of everything above against a fake `az`, checked against the real CLI's help and the Batch SDK models |
| `CHECKS.md` | what was validated offline, the fixes, and the first-run checklist once Azure is reachable |

Order of operations (first time): `az login` → `01_bootstrap.sh` → `00_quotas.sh` (request Batch
quota, edit `env.sh`) → `02_build_image.sh` → `04_pools.sh` (and `05_gpu_pool.sh`) →
`03_stage_dataset.sh` → `grid.py submit ...`. The step-by-step first run, with what to check
after each step, is in [CHECKS.md](CHECKS.md).

Status: written on 2026-10-01 against the documented `az batch` / harness interfaces. The
in-container part (`runner/run_shard.py` on a parquet shard with a bundled estimator, through
`whest run --runner subprocess`, producing the JSON the collector reads) is **tested locally**;
the Azure provisioning and `grid.py` are **not yet executed**, because this session's network
policy denies management.azure.com. They were validated offline on 2026-10-01 (azure-cli 2.90.0:
every az command and flag, the pool/task JSON against the Batch SDK models, the task command
lines against the runner argparse, the autoscale formula against the documented language;
`tests/run_offline.sh`, findings and fixes in [CHECKS.md](CHECKS.md)). It is run end to
end as soon as that host (and login.microsoftonline.com, graph.microsoft.com,
*.batch.azure.com, *.blob.core.windows.net, *.azurecr.io) is allowed, or from any machine where
`az login` works.

## GPU bakes and moment atlases (added 2026-10-01, afternoon)

- `05_gpu_pool.sh` creates an autoscaling GPU pool (A100 by default) from the HPC image; `runner/Dockerfile.gpu`
  is the bake image (torch 2.4.1 + CUDA 12.4, the official v2-phase2 bake stack), built by `02_build_image.sh`.
  Batch exposes the GPUs to container tasks itself, so tasks do not pass `--gpus`.
- `grid.py bake --name fresh-A --n-mlps 1000 --n-samples 100000000 --slices 16` bakes 1,000 fresh-seed MLPs with
  the official recipe on 16 GPU tasks (≈ 450 s per MLP at N = 1e8 on an A100-class card; N = 1e9 is ≈ 4,500 s on
  H200 per the community atlases) and leaves the slices under `bakes/<name>/`; merge with `whest dataset merge`.
- `grid.py atlas --name d8b --start 0 --count 2048 --pairs` computes per-layer marginal moments to order 6, gate
  probabilities and dense pair blocks for seed-regenerable MLPs (`runner/bake_moments.py`, same conventions as the
  community atlases `keenanpepper/arc-whestbench-p2-*`).
- `runner/stage_dataset.py --all-files --prefix <name>` stages a community dataset repository into blob
  (`HF_DATASET=<repo> HF_REVISION=main STAGE_ARGS="--all-files --prefix <name>" ./03_stage_dataset.sh`).

Community resources worth staging first (Hugging Face, MIT): `keenanpepper/arc-whestbench-p2-higher-moments-2026`
(mini split, N = 1e8 joint moments), `keenanpepper/arc-whestbench-p2-full1000-N1e9` (337 GB, all 1,000 full-split MLPs,
joint moments at N = 1e9), `keenanpepper/arc-whestbench-p2-d8b-corpus-14k` (14,048 seed-regenerable MLPs, marginal
moments + final means at N = 1e8, plus a teacher bake of the augmented factored K=3 reference chain),
`keenanpepper/whest-p2-bakev2-*` (Ω-sketched third cumulants).
