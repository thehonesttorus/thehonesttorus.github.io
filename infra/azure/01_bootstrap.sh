#!/usr/bin/env bash
# One-time: resource group, storage account + containers, container registry, one Batch
# account per region (Batch-service pool allocation, auto-storage linked).  Idempotent.
source "$(dirname "$0")/env.sh"
az group create -n "$RG" -l "$HOME_REGION" -o none
az storage account create -n "$STORAGE" -g "$RG" -l "$HOME_REGION" --sku Standard_LRS --kind StorageV2 --allow-blob-public-access false -o none
for c in "$CONTAINER_DATASET" "$CONTAINER_SUBMISSIONS" "$CONTAINER_RESULTS"; do
  az storage container create -n "$c" --account-name "$STORAGE" --auth-mode login -o none || true
done
az acr create -n "$ACR" -g "$RG" --sku Standard --admin-enabled true -l "$HOME_REGION" -o none
for r in $REGIONS; do
  acct=$(batch_account_name "$r")
  if ! az batch account show -n "$acct" -g "$RG" -o none 2>/dev/null; then
    echo "creating batch account $acct in $r"
    az batch account create -n "$acct" -g "$RG" -l "$r" --storage-account "$STORAGE" -o none || echo "  (failed in $r; drop it from REGIONS or request Batch access)"
  fi
done
echo "bootstrap done: rg=$RG storage=$STORAGE acr=$ACR"
