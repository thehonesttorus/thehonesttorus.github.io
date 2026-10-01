#!/usr/bin/env bash
# shellcheck disable=SC2317  # "exit" after "return" is reached when this file is executed instead of sourced
# Shared configuration for the Azure experiment grid.  Source this file first.
# Everything is parameterised so the same scripts work from this container (once
# management.azure.com is allowed by the environment's network policy) or from a laptop.
# Note: sourcing turns on `set -euo pipefail` in the calling shell (meant for the numbered scripts).

set -euo pipefail

_env_fail() { echo "env.sh: $*" >&2; }

export PREFIX="${PREFIX:-whest}"                 # 1-8 lowercase letters/digits: used in storage/registry/Batch names
[[ "$PREFIX" =~ ^[a-z0-9]{1,8}$ ]] || { _env_fail "PREFIX must be 1-8 lowercase letters or digits (got '$PREFIX')"; return 1 2>/dev/null || exit 1; }
export SUB_ID="${SUB_ID:-$(az account show --query id -o tsv 2>/dev/null || true)}"
[ -n "$SUB_ID" ] || { _env_fail "no Azure subscription: run 'az login' (and 'az account set -s <id>') or export SUB_ID"; return 1 2>/dev/null || exit 1; }
export RG="${RG:-${PREFIX}-rg}"
export HOME_REGION="${HOME_REGION:-eastus}"       # storage and registry live here; must also be in REGIONS (03 stages from its pool)

# Regions to shard Batch pools across.  Batch-service-mode pools draw on per-Batch-account core
# quotas (one account per region), so more regions = more cores.  Order = preference.
# Edit after 00_quotas.sh prints the table.
export REGIONS="${REGIONS:-eastus eastus2 westus2 westus3 northeurope westeurope uksouth centralus southcentralus australiaeast}"

# Pool shape.  D-series v5 are plain x86 CPU boxes with 4 GB per vCPU.  The grader gives a
# solution 2 vCPU + 14 backend vCPU + 8 GB; for throughput we run one evaluation per 4 vCPU
# and relax the wall-clock cap, since FLOP accounting is hardware-independent.
export VM_SIZE="${VM_SIZE:-Standard_D64ds_v5}"    # 64 vCPU, 256 GB; fallback: Standard_D32ds_v5, Standard_F64s_v2
export VCPU_PER_NODE="${VCPU_PER_NODE:-64}"
export VCPU_PER_TASK="${VCPU_PER_TASK:-4}"
export MAX_NODES_PER_POOL="${MAX_NODES_PER_POOL:-20}"   # autoscale ceiling per region; raise after quota increases
export USE_SPOT="${USE_SPOT:-0}"                  # 1 = low-priority/spot nodes (cheaper, preemptible); we do not care about cost, so 0

# Node OS image for container pools, publisher:offer:sku:version.  The previous
# microsoft-azure-batch:ubuntu-server-container:20-04-lts / batch.node.ubuntu 20.04 pair was
# retired on 2025-04-23 (Batch rejects new pools with it).  The documented replacement for
# container workloads is microsoft-dsvm:ubuntu-hpc:2204 with batch.node.ubuntu 22.04 (Docker
# and the NVIDIA container runtime preinstalled, Gen2 only).  Check with:
#   az batch pool supported-images list --query "[?imageReference.offer=='ubuntu-hpc']" -o table
export POOL_IMAGE="${POOL_IMAGE:-microsoft-dsvm:ubuntu-hpc:2204:latest}"
export NODE_AGENT_SKU="${NODE_AGENT_SKU:-batch.node.ubuntu 22.04}"
export OS_DISK_GB="${OS_DISK_GB:-128}"            # docker's data root (images + container /tmp) is on the OS disk of ubuntu-hpc; 0 = image default

# Storage layout
export STORAGE="${STORAGE:-${PREFIX}sa$(echo "$SUB_ID" | tr -d '-' | cut -c1-8)}"   # globally unique, 3-24 lowercase alnum
export CONTAINER_DATASET="dataset"      # parquet shards of arc-whestbench-public-2026@v2-phase2 (+ our own bakes)
export CONTAINER_SUBMISSIONS="submissions"  # estimator bundles (tar.gz of a submission folder)
export CONTAINER_RESULTS="results"      # one JSON per task
export ACR="${ACR:-${PREFIX}acr$(echo "$SUB_ID" | tr -d '-' | cut -c1-8)}"
export IMAGE="${IMAGE:-whest-runner:latest}"      # CPU runner (runner/Dockerfile)
export GPU_IMAGE="${GPU_IMAGE:-whest-bake:latest}" # GPU bake image (runner/Dockerfile.gpu)

# Dataset
export HF_DATASET="${HF_DATASET:-aicrowd/arc-whestbench-public-2026}"
export HF_REVISION="${HF_REVISION:-v2-phase2}"

batch_account_name() { echo "${PREFIX}b$(echo "$1" | tr -d '-' | cut -c1-10)$(echo "$SUB_ID" | tr -d '-' | cut -c1-4)"; }
pool_name() { echo "${PREFIX}-pool-$1"; }

# Shared-key login: the Owner/Contributor roles carry no Batch data-plane permissions, so the
# default Entra ID login gets 403 on pool/job/task calls unless "Azure Batch Data Contributor"
# is assigned.  Shared key needs only listKeys on the account, which Owner/Contributor have.
batch_login() { az batch account login -n "$(batch_account_name "$1")" -g "$RG" --shared-key-auth -o none; }

# Storage account key (data-plane calls with --auth-mode login would need a Storage Blob Data role).
storage_key() { az storage account keys list -n "$STORAGE" -g "$RG" --query '[0].value' -o tsv; }

# UTC time N hours from now in the format az storage expects (portable: no GNU date -d).
utc_in_hours() { python3 -c 'import datetime,sys; print((datetime.datetime.now(datetime.timezone.utc)+datetime.timedelta(hours=float(sys.argv[1]))).strftime("%Y-%m-%dT%H:%MZ"))' "$1"; }

# Autoscale formula: enough nodes for the pending (active + running) tasks, capped, scaled to zero
# when idle.  The autoscale language has no % operator; ceil() is a built-in.  Falls back to the
# current target when fewer than 25% of the last 3 minutes of samples are available (new pool).
#   $1 target variable ($TargetDedicatedNodes or $TargetLowPriorityNodes), $2 task slots per node, $3 max nodes
autoscale_formula() {
  local target="$1" slots="$2" maxn="$3"
  cat <<FORMULA
pct = \$PendingTasks.GetSamplePercent(TimeInterval_Minute * 3);
pending = pct < 25 ? -1 : max(\$PendingTasks.GetSample(TimeInterval_Minute * 3));
need = pending < 0 ? ${target} : ceil(pending / ${slots});
${target} = max(0, min(${maxn}, need));
\$NodeDeallocationOption = taskcompletion;
FORMULA
}

# Pool definition (Batch REST JSON, as accepted by `az batch pool create --json-file`).
#   $1 pool id, $2 VM size, $3 task slots per node, $4 container image to prefetch, $5 autoscale formula
#   ACR_USER / ACR_PASS must be in the environment (registry credentials for image pulls).
pool_json() {
  POOL_ID="$1" POOL_VM_SIZE="$2" POOL_SLOTS="$3" POOL_CONTAINER="$4" POOL_FORMULA="$5" python3 - <<'PY'
import json, os
e = os.environ
pub, offer, sku, ver = e["POOL_IMAGE"].split(":")
registry = e["ACR"] + ".azurecr.io"
pool = {
    "id": e["POOL_ID"],
    "vmSize": e["POOL_VM_SIZE"],
    "virtualMachineConfiguration": {
        "imageReference": {"publisher": pub, "offer": offer, "sku": sku, "version": ver},
        "nodeAgentSKUId": e["NODE_AGENT_SKU"],
        "containerConfiguration": {
            "type": "dockerCompatible",
            "containerImageNames": [e["POOL_CONTAINER"]],
            "containerRegistries": [{"registryServer": registry, "username": e["ACR_USER"], "password": e["ACR_PASS"]}],
        },
    },
    "taskSlotsPerNode": int(e["POOL_SLOTS"]),
    "taskSchedulingPolicy": {"nodeFillType": "pack"},
    "enableAutoScale": True,
    "autoScaleFormula": e["POOL_FORMULA"],
    "autoScaleEvaluationInterval": "PT5M",
}
if int(e.get("OS_DISK_GB") or 0) > 0:
    pool["virtualMachineConfiguration"]["osDisk"] = {"diskSizeGB": int(e["OS_DISK_GB"])}
print(json.dumps(pool, indent=1))
PY
}

# Check a `az batch task create --json-file` result (list of per-task results).  The CLI already
# fails on client errors other than TaskExists; TaskExists rows (re-submission) are reported, not fatal.
check_task_results() {
  python3 -c '
import json, sys
res = json.load(open(sys.argv[1])) or []
exists = [r for r in res if ((r.get("error") or {}).get("code") == "TaskExists")]
bad = [r for r in res if str(r.get("status", "")).lower() != "success" and r not in exists]
if exists:
    print(len(exists), "task(s) already existed and were left as they are", file=sys.stderr)
for r in bad:
    print("task", r.get("taskId"), "not added:", r.get("status"), json.dumps(r.get("error")), file=sys.stderr)
sys.exit(1 if bad else 0)
' "$1"
}
