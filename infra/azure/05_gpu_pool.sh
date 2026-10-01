#!/usr/bin/env bash
# One autoscaling GPU pool (home region by default) for Monte-Carlo bakes and moment atlases.
# GPU core quota (Batch account quota for the NCADS_A100_v4 / NCads_H100_v5 family) is usually 0
# on a new Batch account: request it first (00_quotas.sh prints it).  The ubuntu-hpc node image
# carries the NVIDIA driver, Docker and the NVIDIA container runtime; Batch exposes the GPUs to
# container tasks by itself (tasks must not pass --gpus).  Build $GPU_IMAGE with 02_build_image.sh.
source "$(dirname "$0")/env.sh"
export GPU_REGION="${GPU_REGION:-$HOME_REGION}"
export GPU_VM_SIZE="${GPU_VM_SIZE:-Standard_NC24ads_A100_v4}"   # 1x A100 80GB; alt: Standard_NC40ads_H100_v5, Standard_NC4as_T4_v3
export GPU_MAX_NODES="${GPU_MAX_NODES:-8}"
batch_login "$GPU_REGION"
ACR_USER=$(az acr credential show -n "$ACR" --query username -o tsv)
ACR_PASS=$(az acr credential show -n "$ACR" --query 'passwords[0].value' -o tsv)
export ACR_USER ACR_PASS
pool="${PREFIX}-gpu-${GPU_REGION}"
# shellcheck disable=SC2016  # literal autoscale variable name
formula=$(autoscale_formula '$TargetDedicatedNodes' 1 "$GPU_MAX_NODES")
umask 077
tmpd=$(mktemp -d); trap 'rm -rf "$tmpd"' EXIT
if az batch pool show --pool-id "$pool" -o none 2>/dev/null; then
  az batch pool autoscale enable --pool-id "$pool" --auto-scale-formula "$formula" --auto-scale-evaluation-interval PT5M -o none
  echo "updated $pool"
else
  pool_json "$pool" "$GPU_VM_SIZE" 1 "$ACR.azurecr.io/$GPU_IMAGE" "$formula" > "$tmpd/gpu_pool.json"
  az batch pool create --json-file "$tmpd/gpu_pool.json" -o none
  echo "created $pool ($GPU_VM_SIZE, max $GPU_MAX_NODES nodes)"
fi
