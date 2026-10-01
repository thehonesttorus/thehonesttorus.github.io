#!/usr/bin/env bash
# One-time: resource providers, resource group, storage account + containers, container
# registry, one Batch account per region (Batch-service pool allocation).  Idempotent.
source "$(dirname "$0")/env.sh"
for ns in Microsoft.Batch Microsoft.Storage Microsoft.ContainerRegistry; do
  az provider register -n "$ns" --wait -o none     # new subscriptions may not have these registered
done
az group create -n "$RG" -l "$HOME_REGION" -o none
if ! az storage account show -n "$STORAGE" -g "$RG" -o none 2>/dev/null; then
  az storage account create -n "$STORAGE" -g "$RG" -l "$HOME_REGION" --sku Standard_LRS --kind StorageV2 \
    --allow-blob-public-access false --min-tls-version TLS1_2 -o none
fi
key=$(storage_key)
for c in "$CONTAINER_DATASET" "$CONTAINER_SUBMISSIONS" "$CONTAINER_RESULTS"; do
  # account key, not --auth-mode login: login mode needs a Storage Blob Data role, which Owner lacks
  az storage container create -n "$c" --account-name "$STORAGE" --account-key "$key" -o none
done
if ! az acr show -n "$ACR" -g "$RG" -o none 2>/dev/null; then
  az acr create -n "$ACR" -g "$RG" --sku Standard --admin-enabled true -l "$HOME_REGION" -o none
fi
for r in $REGIONS; do
  acct=$(batch_account_name "$r")
  if ! az batch account show -n "$acct" -g "$RG" -o none 2>/dev/null; then
    echo "creating batch account $acct in $r"
    # no --storage-account: tasks use SAS URLs, nothing needs auto-storage, and the storage
    # account lives in HOME_REGION only
    az batch account create -n "$acct" -g "$RG" -l "$r" -o none || echo "  (failed in $r; drop it from REGIONS or request Batch access)"
  fi
done
echo "bootstrap done: rg=$RG storage=$STORAGE acr=$ACR"
