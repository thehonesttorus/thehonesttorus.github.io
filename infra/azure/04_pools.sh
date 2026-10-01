#!/usr/bin/env bash
# Create (or update the autoscale formula of) one autoscaling container pool per region.
# Autoscale: enough nodes for the pending tasks, capped at MAX_NODES_PER_POOL, scaled to zero
# when idle (formula in env.sh).  Changing VM_SIZE / POOL_IMAGE of an existing pool needs
# `az batch pool delete --pool-id <pool>` first (pools cannot be re-imaged in place).
source "$(dirname "$0")/env.sh"
ACR_USER=$(az acr credential show -n "$ACR" --query username -o tsv)
ACR_PASS=$(az acr credential show -n "$ACR" --query 'passwords[0].value' -o tsv)
export ACR_USER ACR_PASS
slots=$(( VCPU_PER_NODE / VCPU_PER_TASK ))
# shellcheck disable=SC2016  # literal autoscale variable names, not shell expansions
if [ "$USE_SPOT" = "1" ]; then target='$TargetLowPriorityNodes'; else target='$TargetDedicatedNodes'; fi
formula=$(autoscale_formula "$target" "$slots" "$MAX_NODES_PER_POOL")
umask 077                                    # the pool JSON carries the registry password
tmpd=$(mktemp -d); trap 'rm -rf "$tmpd"' EXIT
for r in $REGIONS; do
  batch_login "$r" 2>/dev/null || { echo "no batch account in $r, skipping"; continue; }
  pool=$(pool_name "$r")
  if az batch pool show --pool-id "$pool" -o none 2>/dev/null; then
    az batch pool autoscale enable --pool-id "$pool" --auto-scale-formula "$formula" --auto-scale-evaluation-interval PT5M -o none
    echo "updated autoscale for $pool"
  else
    pool_json "$pool" "$VM_SIZE" "$slots" "$ACR.azurecr.io/$IMAGE" "$formula" > "$tmpd/pool_$r.json"
    az batch pool create --json-file "$tmpd/pool_$r.json" -o none
    echo "created $pool ($VM_SIZE, $slots slots/node, max $MAX_NODES_PER_POOL nodes)"
  fi
done
