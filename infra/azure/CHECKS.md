# Offline validation of the Azure Batch grid

*2026-10-01. Everything below was checked without reaching Azure (management.azure.com is blocked
from this container). The end of the document lists what can only be verified with the network
open, in the order to run it.*

## How it was checked

| item | version / source |
|---|---|
| Azure CLI | 2.90.0 (`/root/azcli`, core 2.90.0; batch module on `azure-batch` 15.0.0b1, `azure-mgmt-batch` 17.3.0) |
| quota extension | not installable here (aka.ms is blocked); its argument definitions were read from `Azure/azure-cli-extensions` `src/quota` (API 2023-02-01) |
| shellcheck | 0.11.0 (`pip install shellcheck-py`) |
| Python | 3.11 (`py_compile` on every `.py`), whestbench 0.16.1 / flopscope 0.12.1 for the runner argparse |
| Azure documentation (via Exa) | *Container workloads on Azure Batch* and *Autoscale compute nodes in a Batch pool* (both ms.date 2026-05-19), *Batch service quotas and limits*, *Role-based access control for Azure Batch*, Azure/Batch retirement tracker (issue #140), ACR Tasks free-trial Q&A |

**Re-run it all:** `REAL_AZ=/root/azcli/bin/az RUNNER_PYTHON=/root/whest/bin/python infra/azure/tests/run_offline.sh`
(about 2 minutes, mostly fetching `az ... --help` once per command). It runs every numbered script and every
`grid.py` subcommand against a fake `az` (`tests/bin/az` → `tests/fake_az.py`), which logs each call and
captures each `--json-file` payload. It then checks:

1. **every az call** against the help of the real installed CLI: the command group and subcommand
   exist, and every flag is an argument of that command (`tests/check_az_calls.py`);
2. **every pool and task JSON** against the Batch SDK models the CLI deserializes them into
   (`BatchPoolCreateContent`, `BatchTaskCreateContent`). It checks property names recursively, enum
   values and documented-required fields (`tests/validate_batch_json.py`). This matters because these
   models keep unknown keys silently, so a typo would never error offline;
3. **every task command line** (`tests/check_payloads.py`). It must be `/bin/bash -c '<cmd>'`. The
   `<cmd>` is shell-split as bash would split it, then fed to the argparse of the runner script it
   calls. Task ids must be valid and unique, no `--gpus`, admin user. Every **autoscale formula** is
   linted against the documented operators, functions, variables and sample methods;
4. negative cases: bad `PREFIX`, `HOME_REGION` outside `REGIONS`, bad `--name`, duplicate bundle
   tags, bad `--runner` must all be rejected.

**Result on the fixed tree:** all steps pass. 330 az calls (34 distinct command+flag sets), 0
failing. 30 JSON objects validate with 0 problems. 18 task command lines and 12 autoscale formulas
pass with 0 problems. shellcheck is clean, and py_compile is clean. The 17 distinct az commands
quoted in the first-run checklist below were checked the same way: 0 failing (flags only; the output
keys their `--query` expressions read were checked against `azure.cli.core.util.todict` of the SDK models).

**Negative control:** the same harness run on the original files (commit 995baf8, with `FAKE_AZ_ACR_PASS` set to a password without quotes: the default one contains a `"`, which breaks the original heredoc pool JSON, F21, before the formula lint is reached) flags the `%`
operator in all 11 CPU pool formulas. It also flags `--gpus` in all 7 GPU tasks, 4 invalid task
ids, the stage task's unwrapped command line, non-admin container tasks, the `collect` crash and the
`--shards 7` crash. Its az calls and JSON property names were all valid: the original's defects
were semantic, not syntactic.

## Checks and results

| # | check | result |
|---|---|---|
| C1 | every `az` group/subcommand/flag in `env.sh`, `00`–`05`, `grid.py` exists in 2.90.0 | pass, both before and after the fixes (list in `az_calls.txt` of a harness run) |
| C2 | `az storage blob list --num-results "*"` | valid ("Provide \"*\" to return all") |
| C3 | `--expiry` format `%Y-%m-%dT%H:%MZ` for `generate-sas` | valid (help: `Y-m-d'T'H:M'Z'`) |
| C4 | `az batch task create --json-file` with a JSON **array** | valid: the CLI bulk-adds and chunks by 100 itself. It raises on client errors other than `TaskExists` and returns per-task rows `{status, taskId, error}` with exit 0 (see F11) |
| C5 | `az batch pool create --json-file` shape | REST camelCase JSON deserialized into `BatchPoolCreateContent`; all keys valid, including `nodeAgentSKUId`, `containerConfiguration.type=dockerCompatible`, `taskSchedulingPolicy.nodeFillType=pack`, `osDisk.diskSizeGB` |
| C6 | container task JSON | `containerSettings{imageName, containerRunOptions}`, `userIdentity.autoUser{scope, elevationLevel}`, `constraints{maxWallClockTime, maxTaskRetryCount}`: valid; the registry may be omitted on the task because it is given on the pool |
| C7 | node image for container pools | **original was retired** → F1 |
| C8 | autoscale formula syntax | **original invalid** → F2, F3 |
| C9 | auth for Batch data plane and blob data plane | **original would get 403s** → F4, F5, F6 |
| C10 | `az vm list-usage` output fields | `currentValue`/`limit` (int), `name.value`: as the script reads them |
| C11 | `az batch account show` output keys | az prints `dedicatedCoreQuotaPerVmFamily` (camel-cased SDK attribute), **not** the REST `...PerVMFamily` the original read → F18 |
| C12 | `az batch job task-counts show` output | `{taskCounts: {active, running, completed, succeeded, failed}, taskSlotCounts: ...}`: as `status` reads it |
| C13 | `az quota create` (comment in 00) | the extension takes `--resource-name --scope --limit-object value=N --resource-type dedicated`; the original comment had a stray `name=...` token and no `--resource-type` → F18 |
| C14 | resource names | storage `whestsa<8 hex>` (15), ACR `whestacr<8 hex>` (16), Batch `whestb<region≤10><4 hex>` (≤ 20): lowercase alnum, within limits; distinct across the 10 regions |
| C15 | `taskSlotsPerNode` 16 on 64 vCPU | ≤ min(4 × vCPU, 256): valid; `autoScaleEvaluationInterval` PT5M = documented minimum |
| C16 | shellcheck 0.11.0 on all scripts | original: info-level only (SC2016 intended, SC2086 on `/tmp/pool_$r.json`); fixed: clean |
| C17 | `py_compile` on all Python | pass |
| C18 | `grid.py` argument parsing and every code path with dummy values | submit (2 bundles × shards 0–3, `--extra`, single shard, TaskExists), status `--failures`, collect `--csv` on a real `run_shard.py` result, bake, atlas `--pairs`: pass; original crashed on `--shards 7` and on `collect` with an empty bundle → F14, F15 |
| C19 | task command lines → runner argparse | all accepted by `run_shard.py`, `stage_dataset.py`, `bake_shard.py`, `bake_moments.py` |
| C20 | `whest run` / `whest dataset bake` flags used by the runners | all exist in whestbench 0.16.1 (`--estimator --dataset --split --runner --format json --wall-time-limit --max-threads --profile --n-mlps`; `--n-mlps --n-samples --width --depth --split --mlp-seeds --output --slice K/N --chunk-size --torch --device`) |
| C21 | Dockerfile pins | `whestbench>=0.16.1,<0.17` and `flopscope>=0.12.1,<0.13` resolve on PyPI (Requires-Python ≥ 3.10; images ship 3.11); `pytorch/pytorch:2.4.1-cuda12.4-cudnn9-runtime` exists on Docker Hub; `whestbench[gpu]` extra exists |

## Bugs found and fixed

Ordered by impact. "Would have" describes the first run with the network open.

| # | file | defect | would have | fix |
|---|---|---|---|---|
| F1 | 04, 05 (`env.sh`) | Node image `microsoft-azure-batch/ubuntu-server-container/20-04-lts` with `batch.node.ubuntu 20.04`. Ubuntu 20.04 support in Batch was retired on **2025-04-23**, and Microsoft now marks the `microsoft-azure-batch` container images deprecated | every CPU pool creation rejected | `POOL_IMAGE=microsoft-dsvm:ubuntu-hpc:2204:latest`, `NODE_AGENT_SKU="batch.node.ubuntu 22.04"`: the documented container image (Docker/Moby and NVIDIA runtime preinstalled, Gen2; D64ds_v5, F64s_v2 and NC A100 v4 all support Gen2). `ubuntu-hpc` also has `2404` SKUs in use; switch both variables when 22.04 nears end of life |
| F2 | 04 | Autoscale formula used `pending % slots`: `%` is not in the autoscale operator table | every pool create fails with an invalid-formula error (`autoscale enable` too) | `ceil(pending / slots)` (`ceil` is a documented built-in) |
| F3 | 04, 05 | `max($PendingTasks.GetSample(TimeInterval_Minute * 2, 0))` with no guard for missing samples (new pool, sample lag) | formula evaluation errors / no scaling while samples are missing | Documented pattern: `GetSamplePercent` guard, keep the current target when fewer than 25 % of samples exist, `max(0, min(MAX, need))` (one `autoscale_formula` helper in `env.sh` for both pools) |
| F4 | 03, 04, 05, grid.py | `az batch account login` without `--shared-key-auth` uses Entra ID. Owner/Contributor carry no Batch *DataActions*, and Microsoft's tutorial states the "Azure Batch Data Contributor" role "is required to create pools, jobs, and tasks" | 403 on every pool/job/task call | `batch_login` helper with `--shared-key-auth` (needs only listKeys on the account, which Owner/Contributor have) |
| F24 | env.sh, grid.py | (found in review) `az batch account login` stores the current account in `~/.azure/config`, which every az process shares. `grid.py` and the scripts switch accounts region by region, so a second az process (`grid.py status` while `submit` runs, or `04_pools.sh`) moves the first one to another region's account mid-run. `az batch job create` accepts a pool id that does not exist in that account, so the tasks sit pending with no error | jobs and tasks silently created in the wrong region's account, never scheduled | `batch_login` (both copies) takes the credentials from `login --shared-key-auth --show` and exports them as `AZURE_BATCH_ACCOUNT`/`AZURE_BATCH_ENDPOINT`/`AZURE_BATCH_ACCESS_KEY`/`AZURE_BATCH_AUTH_MODE`, which az reads before the config file (knack `CLIConfig.get`; the env names are in `az batch job set --help`) and which stay private to the process. The fake az now refuses data-plane calls made without them |
| F5 | 01 | `az storage container create --auth-mode login` needs a *Storage Blob Data* role (Owner lacks data actions), and `\|\| true` swallowed the failure | containers never created; every upload/SAS step fails later with a confusing error | account key (`storage_key` helper), no `\|\| true` (create is idempotent) |
| F6 | 03 | User-delegation SAS (`--auth-mode login --as-user`) needs `generateUserDelegationKey` (a Storage Blob Data role) | staging never starts | account-key SAS, like grid.py |
| F7 | 03 (+ stage_dataset.py) | Stage task ran as the default non-admin user. Batch maps the task user into the container, and `stage_dataset.py` did `makedirs("/data/stage")` | `PermissionError` at start | `autoUser {scope: pool, elevationLevel: admin}` (root in the container, as the docs advise); scratch and `HF_HOME` (incl. the hf_xet chunk cache) moved to `$AZ_BATCH_TASK_WORKING_DIR` |
| F8 | 03 | `commandLine` not wrapped in a shell, with `'...'` quoting around the SAS URL. Batch does not run command lines under a shell | quotes and `&` passed literally / mis-split | task JSON built in Python: `/bin/bash -c` + `shlex` quoting (as grid.py) |
| F9 | grid.py (bake, atlas) | `containerRunOptions: "--gpus all"`. The docs: Batch enables the GPUs for container tasks on GPU pools, "you shouldn't include the `--gpus` argument" | container create failure or duplicated device request on every GPU task | removed (`--shm-size 8g` kept) |
| F10 | 02 | Only the CPU image was built; `whest-bake:latest` (05 prefetch, bake/atlas tasks) was never built | GPU pool image prefetch fails, nodes unusable | 02 builds `Dockerfile.gpu` → `$GPU_IMAGE` as well (`SKIP_GPU_IMAGE=1` to skip). If ACR Tasks are refused (`TasksOperationsNotAllowed`, free-trial subscriptions), it falls back to local `docker build` + `az acr login` + `docker push` |
| F11 | grid.py, 03 | Bulk-add results not inspected. A client error aborted grid.py mid-way **before the manifest was written**, and `TaskExists` rows were silent | partial submissions that `status`/`collect` cannot see | per-task status rows checked (`TaskExists` reported as "left as is"); manifest written in a `finally` |
| F12 | grid.py | Round-robin over all `REGIONS` although 04 skips regions without a Batch account | submit crashes at the first such region, after uploading and submitting the earlier ones | `usable_regions`: login + pool check first, skip and report the rest |
| F13 | grid.py | Task ids: only `_` replaced (needlessly), `.` kept (`est_v1.2` → invalid id → bulk add rejected); tags truncated to 30 chars with no collision check; `--tag` passed unquoted | rejected batches; two bundles sharing a basename silently overwrite each other's results | ids sanitized to `[A-Za-z0-9_-]{≤64}`; duplicate ids or bundle tags rejected; whole command `shlex.join`-quoted |
| F14 | grid.py | `collect` printed `None` with `:>6`, which raised `TypeError` when any bundle had no results yet | `collect` unusable while a grid is running | formatter handles `None` (verified on a real `run_shard.py` result) |
| F15 | grid.py | `--shards 7` (single index) raised `ValueError` | crash | accepts `A-B` or `N` |
| F16 | grid.py, 03 | Jobs never completed. Each submission left one active job per region, and active jobs are capped at 100–300 per Batch account | job creation fails after ~100 experiments | after adding tasks: `az batch job set --on-all-tasks-complete terminatejob` (set after the tasks, as the REST docs require for jobs without a job-manager task); `--keep-job-open` to opt out |
| F17 | runner | Shard/bake scratch in the container's `/tmp`. On `ubuntu-hpc`, Docker's data root (container layers) is on the OS disk, shared by 16 concurrent tasks (≈ 1.1 GB parquet each) | OS-disk exhaustion → nodes go *Unusable* (a known `ubuntu-hpc` failure mode) | scratch under `$AZ_BATCH_TASK_WORKING_DIR` (resource disk, 2.4 TB on D64ds_v5); `osDisk.diskSizeGB=128` (`OS_DISK_GB`, 0 = image default); `HF_HOME=/tmp/hf` instead of `/data/hf` in both images |
| F18 | 00 | Printed subscription vCPU quotas, which do not limit Batch-service-mode pools; read `dedicatedCoreQuotaPerVMFamily` (az prints `...PerVmFamily`, so the families column was always empty); wrong `az quota create` example | misleading quota table | Batch account quotas first (per family of `VM_SIZE`/`GPU_VM_SIZE`, resolved with `az batch location list-skus`), accounts-per-region, then subscription quotas marked informational; corrected `az quota create` |
| F19 | 01 | Not idempotent for ACR/storage re-runs; no resource-provider registration; every Batch account linked the HOME_REGION storage account as auto-storage although nothing uses auto-storage | first run on a fresh subscription may fail with `MissingSubscriptionRegistration`; a needless cross-region link | `az provider register --wait` for Batch/Storage/ContainerRegistry; show-before-create for storage and ACR; `--storage-account` dropped |
| F20 | env.sh | Without `az login`, `SUB_ID` was empty and names silently became `whestsa`/`whestacr` (non-unique); `PREFIX` unchecked | resources created under colliding names or opaque failures | fail fast with a message; `PREFIX` must be 1–8 lowercase alnum |
| F21 | 04, 05 | Pool JSON with the ACR password written to world-readable `/tmp/pool_*.json` via a heredoc (a `"` in any value breaks the JSON) | credential left on disk | JSON built by Python (`pool_json` helper), `umask 077`, `mktemp`, deleted on exit |
| F22 | 03 | GNU-only `date -d '+2 days'` | fails on macOS | `utc_in_hours` (Python) |
| F23 | grid.py | SAS lifetime fixed at 72 h; `env()` hid `env.sh` errors | tasks queued behind quota for > 3 days fail on download; opaque `KeyError`s | `--sas-hours` (default 168); `env.sh` failures reported |

Also changed: `grid.py` defaults bake/atlas images to `GPU_IMAGE` from `env.sh` and checks that the
GPU pool exists. It also validates `--runner` choices, skips `__pycache__` when packing a bundle,
deletes the temporary tarball, and passes `--no-progress` to `download-batch`. `03_stage_dataset.sh`
takes `HF_DATASET`/`HF_REVISION`/`STAGE_ARGS` overrides (community atlases) and refuses to run when
`HOME_REGION` has no pool.

## What can only be verified with the network open (run in this order)

Allow: `management.azure.com`, `login.microsoftonline.com`, `graph.microsoft.com`, `*.batch.azure.com`,
`*.blob.core.windows.net`, `*.azurecr.io` (plus `aka.ms` and `azcliextensionsync.blob.core.windows.net`
only for `az extension add --name quota`). Then, from `infra/azure`, run the numbered scripts as
they are. Commands below that use the `env.sh` helpers (`storage_key`, `batch_login`,
`batch_account_name`, `autoscale_formula`) assume a scratch shell prepared with `bash`, then
`source ./env.sh; set +eu`. `env.sh` turns on `set -euo pipefail`, which would otherwise close the
shell on the first failing command.

1. **Login.** `az login --use-device-code`, then `az account show -o table` (and `az account set -s <id>`).
   In the scratch shell, `echo $STORAGE $ACR $(batch_account_name eastus)` should print names with the
   8-hex subscription suffix. Note the subscription type: a free trial blocks ACR Tasks (step 5) and has
   tiny Batch quotas.
2. **`./01_bootstrap.sh`.** Expect `bootstrap done`. Check:
   - `az batch account list -g whest-rg -o table`: one account per region. Drop from `REGIONS` any
     region that printed "(failed in …)".
   - `az storage container list --account-name $STORAGE --account-key "$(storage_key)" -o table`:
     `dataset`, `submissions`, `results`.
   - `az storage account show -n $STORAGE -g whest-rg --query allowSharedKeyAccess -o tsv` must not
     be `false`. An Azure Policy can disable shared keys; the grid relies on account-key SAS.
   - `az batch account show -n "$(batch_account_name eastus)" -g whest-rg --query allowedAuthenticationModes -o tsv`
     must include `SharedKey`.
3. **`./00_quotas.sh`.** In the Batch account table, `standardDDSv5Family` (and the GPU family) must
   be ≥ 64 × the nodes you want per region. Request increases in the portal (Batch account > Quotas)
   where it is 0 or small, and set `MAX_NODES_PER_POOL` to at most quota / 64. Every region must
   resolve the family name; a `?` means `VM_SIZE` is not offered by Batch there.
4. **Image still supported.** After `batch_login eastus`:
   `az batch pool supported-images list --query "[?imageReference.offer=='ubuntu-hpc'].{sku:imageReference.sku,agent:nodeAgentSkuId,caps:capabilities,eol:batchSupportEndOfLife}" -o table`.
   (az prints the snake_case SDK attribute camel-cased, `nodeAgentSkuId`; the REST spelling
   `nodeAgentSKUId` would come back empty.) Expect `2204` / `batch.node.ubuntu 22.04` with
   `DockerCompatible`. If it shows an end-of-life date before 17 Oct 2026 (the Phase 2 submission
   deadline) or is missing, switch `POOL_IMAGE`/`NODE_AGENT_SKU` to the listed `2404` /
   `batch.node.ubuntu 24.04` pair.
5. **`./02_build_image.sh`.** Both images build. `az acr repository list -n $ACR -o table` shows
   `whest-runner` and `whest-bake`. On `TasksOperationsNotAllowed` (free-trial subscription), either
   upgrade to pay-as-you-go or run it where `docker` is installed (the script falls back to a local
   build + push). Both base images come from Docker Hub; if a cloud build fails on `toomanyrequests`
   (Docker Hub pull limit on the shared build agents), import them once with
   `az acr import -n $ACR --source docker.io/library/python:3.11-slim` (and the `pytorch/pytorch` tag)
   and point the `FROM` lines at `$ACR.azurecr.io/...`.
6. **`./04_pools.sh`** (and `./05_gpu_pool.sh` once GPU quota exists). Not checkable offline: the
   service accepting `osDisk.diskSizeGB` (if rejected, `OS_DISK_GB=0 ./04_pools.sh`) and the formula
   evaluating at runtime. Check both on one pool:
   - `az batch pool autoscale evaluate --pool-id whest-pool-eastus --auto-scale-formula "$(autoscale_formula '$TargetDedicatedNodes' 16 20)"`:
     the result string should contain `$TargetDedicatedNodes=0` and no `error`.
   - `az batch pool show --pool-id whest-pool-eastus --query "{state:state,alloc:allocationState,resizeErrors:resizeErrors,autoScaleRun:autoScaleRun}" -o json`:
     `autoScaleRun.error` must be null.
7. **Stage the small split first.** `STAGE_ARGS="--splits mini" ./03_stage_dataset.sh`. Within ≤ 5
   min plus node boot (≈ 5–10 min incl. image pull) the pool grows to one node:
   `az batch node list --pool-id whest-pool-eastus --query "[].{id:id,state:state}" -o table`.
   This verifies at runtime the container-task user mapping, HF reachability from Azure and the SAS
   upload. On failure, read `az batch task file download --job-id <job> --task-id stage --file-path stderr.txt --destination ./stage.err`.
   Expect `v2-phase2/data/mini-0000x-of-00007.parquet` blobs.
8. **Smoke grid.** `./grid.py submit --name smoke --bundle <estimator dir> --split mini --shards 0 --n-mlps 2`,
   then `./grid.py status --name smoke --failures` until done, then `./grid.py collect --name smoke`.
   Compare with a local `runner/run_shard.py` run of the same bundle and shard: `flops_used` and MSE
   must match exactly (only wall-clock differs). The job should reach *completed* by itself (F16).
9. **Full split.** `./03_stage_dataset.sh` (full + mini, 77 GB; a few hours in one task), then the
   real grid. Watch for `AccountVMSeriesCoreQuotaReached` in `resizeErrors` (step 3) and for the pool
   shrinking to 0 after the tasks finish.
10. **GPU.** After 05: `./grid.py atlas --name smoke --count 1 --n-samples 1000000`. Read the task's
    `stdout.txt`: `bake_moments.py` silently falls back to CPU when CUDA is not visible, so check the
    throughput is GPU-class. Then `./grid.py bake --name fresh-A --n-mlps 16 --slices 1 --n-samples 1000000`
    and `whest dataset merge` on the downloaded slice.

Cleanup when done: `az batch pool delete --pool-id <pool> --yes` (pools scale to zero by themselves,
but they keep counting toward the pool quota), `az batch job delete --job-id <job> --yes`.
