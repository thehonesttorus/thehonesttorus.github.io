#!/usr/bin/env bash
# Create (or update) one autoscaling container pool per region.  Autoscale formula: enough
# nodes for the pending tasks, capped at MAX_NODES_PER_POOL, scaled to zero when idle.
source "$(dirname "$0")/env.sh"
acr_user=$(az acr credential show -n "$ACR" --query username -o tsv)
acr_pass=$(az acr credential show -n "$ACR" --query 'passwords[0].value' -o tsv)
slots=$(( VCPU_PER_NODE / VCPU_PER_TASK ))
for r in $REGIONS; do
  acct=$(batch_account_name "$r")
  az batch account login -n "$acct" -g "$RG" 2>/dev/null || { echo "no batch account in $r, skipping"; continue; }
  pool=$(pool_name "$r")
  if [ "$USE_SPOT" = "1" ]; then target='$TargetLowPriorityNodes'; else target='$TargetDedicatedNodes'; fi
  formula="pending = max(\$PendingTasks.GetSample(TimeInterval_Minute * 2, 0));
need = pending / $slots + (pending % $slots > 0 ? 1 : 0);
$target = min($MAX_NODES_PER_POOL, need);
\$NodeDeallocationOption = taskcompletion;"
  cat > /tmp/pool_$r.json <<JSON
{
  "id": "$pool",
  "vmSize": "$VM_SIZE",
  "virtualMachineConfiguration": {
    "imageReference": {"publisher": "microsoft-azure-batch", "offer": "ubuntu-server-container", "sku": "20-04-lts", "version": "latest"},
    "nodeAgentSKUId": "batch.node.ubuntu 20.04",
    "containerConfiguration": {
      "type": "dockerCompatible",
      "containerImageNames": ["$ACR.azurecr.io/$IMAGE"],
      "containerRegistries": [{"registryServer": "$ACR.azurecr.io", "username": "$acr_user", "password": "$acr_pass"}]
    }
  },
  "taskSlotsPerNode": $slots,
  "taskSchedulingPolicy": {"nodeFillType": "pack"},
  "enableAutoScale": true,
  "autoScaleFormula": $(python3 -c 'import json,sys; print(json.dumps(sys.argv[1]))' "$formula"),
  "autoScaleEvaluationInterval": "PT5M"
}
JSON
  if az batch pool show --pool-id "$pool" -o none 2>/dev/null; then
    az batch pool autoscale enable --pool-id "$pool" --auto-scale-formula "$formula" --auto-scale-evaluation-interval PT5M -o none
    echo "updated autoscale for $pool"
  else
    az batch pool create --json-file /tmp/pool_$r.json -o none && echo "created $pool ($VM_SIZE, $slots slots/node, max $MAX_NODES_PER_POOL nodes)"
  fi
done
