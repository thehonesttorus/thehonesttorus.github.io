#!/usr/bin/env bash
# One autoscaling GPU pool (home region by default) for Monte-Carlo bakes and moment atlases.
# GPU quotas (NCads A100 v4 / NC H100 v5 families) are usually 0 on a fresh subscription:
# request them first (00_quotas.sh prints the families). The ubuntu-hpc image carries the
# NVIDIA driver and docker; tasks run the CUDA image with --gpus all.
source "$(dirname "$0")/env.sh"
export GPU_REGION="${GPU_REGION:-$HOME_REGION}"
export GPU_VM_SIZE="${GPU_VM_SIZE:-Standard_NC24ads_A100_v4}"   # 1x A100 80GB; alt: Standard_NC40ads_H100_v5, Standard_NC4as_T4_v3
export GPU_MAX_NODES="${GPU_MAX_NODES:-8}"
export GPU_IMAGE="${GPU_IMAGE:-whest-bake:latest}"
acct=$(batch_account_name "$GPU_REGION")
az batch account login -n "$acct" -g "$RG"
acr_user=$(az acr credential show -n "$ACR" --query username -o tsv)
acr_pass=$(az acr credential show -n "$ACR" --query 'passwords[0].value' -o tsv)
pool="${PREFIX}-gpu-${GPU_REGION}"
formula="pending = max(\$PendingTasks.GetSample(TimeInterval_Minute * 2, 0));
\$TargetDedicatedNodes = min($GPU_MAX_NODES, pending);
\$NodeDeallocationOption = taskcompletion;"
cat > /tmp/gpu_pool.json <<JSON
{
  "id": "$pool",
  "vmSize": "$GPU_VM_SIZE",
  "virtualMachineConfiguration": {
    "imageReference": {"publisher": "microsoft-dsvm", "offer": "ubuntu-hpc", "sku": "2204", "version": "latest"},
    "nodeAgentSKUId": "batch.node.ubuntu 22.04",
    "containerConfiguration": {
      "type": "dockerCompatible",
      "containerImageNames": ["$ACR.azurecr.io/$GPU_IMAGE"],
      "containerRegistries": [{"registryServer": "$ACR.azurecr.io", "username": "$acr_user", "password": "$acr_pass"}]
    }
  },
  "taskSlotsPerNode": 1,
  "enableAutoScale": true,
  "autoScaleFormula": $(python3 -c 'import json,sys; print(json.dumps(sys.argv[1]))' "$formula"),
  "autoScaleEvaluationInterval": "PT5M"
}
JSON
if az batch pool show --pool-id "$pool" -o none 2>/dev/null; then
  az batch pool autoscale enable --pool-id "$pool" --auto-scale-formula "$formula" --auto-scale-evaluation-interval PT5M -o none; echo "updated $pool"
else
  az batch pool create --json-file /tmp/gpu_pool.json -o none && echo "created $pool ($GPU_VM_SIZE, max $GPU_MAX_NODES nodes)"
fi
