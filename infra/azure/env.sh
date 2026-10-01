#!/usr/bin/env bash
# Shared configuration for the Azure experiment grid.  Source this file first.
# Everything is parameterised so the same scripts work from this container (once
# management.azure.com is allowed by the environment's network policy) or from a laptop.

set -euo pipefail

export PREFIX="${PREFIX:-whest}"                 # short, lowercase, <= 8 chars: used in resource names
export SUB_ID="${SUB_ID:-$(az account show --query id -o tsv 2>/dev/null || true)}"
export RG="${RG:-${PREFIX}-rg}"
export HOME_REGION="${HOME_REGION:-eastus}"       # storage, registry, key vault live here

# Regions to shard Batch pools across.  Quotas are per region per VM family, so more
# regions = more cores.  Order = preference.  Edit after 00_quotas.sh prints the table.
export REGIONS="${REGIONS:-eastus eastus2 westus2 westus3 northeurope westeurope uksouth centralus southcentralus australiaeast}"

# Pool shape.  D-series v5 are plain x86 CPU boxes with 4 GB per vCPU.  The grader gives a
# solution 2 vCPU + 14 backend vCPU + 8 GB; for throughput we run one evaluation per 4 vCPU
# and relax the wall-clock cap, since FLOP accounting is hardware-independent.
export VM_SIZE="${VM_SIZE:-Standard_D64ds_v5}"    # 64 vCPU, 256 GB; fallback: Standard_D32ds_v5, Standard_F64s_v2
export VCPU_PER_NODE="${VCPU_PER_NODE:-64}"
export VCPU_PER_TASK="${VCPU_PER_TASK:-4}"
export MAX_NODES_PER_POOL="${MAX_NODES_PER_POOL:-20}"   # autoscale ceiling per region; raise after quota increases
export USE_SPOT="${USE_SPOT:-0}"                  # 1 = low-priority/spot nodes (cheaper, preemptible); we do not care about cost, so 0

# Storage layout
export STORAGE="${STORAGE:-${PREFIX}sa$(echo "$SUB_ID" | tr -d '-' | cut -c1-8)}"   # globally unique, <= 24 chars
export CONTAINER_DATASET="dataset"      # parquet shards of arc-whestbench-public-2026@v2-phase2 (+ our own bakes)
export CONTAINER_SUBMISSIONS="submissions"  # estimator bundles (tar.gz of a submission folder)
export CONTAINER_RESULTS="results"      # one JSON per task
export ACR="${ACR:-${PREFIX}acr$(echo "$SUB_ID" | tr -d '-' | cut -c1-8)}"
export IMAGE="${IMAGE:-whest-runner:latest}"

# Dataset
export HF_DATASET="aicrowd/arc-whestbench-public-2026"
export HF_REVISION="v2-phase2"

batch_account_name() { echo "${PREFIX}b$(echo "$1" | tr -d '-' | cut -c1-10)$(echo "$SUB_ID" | tr -d '-' | cut -c1-4)"; }
pool_name() { echo "${PREFIX}-pool-$1"; }
