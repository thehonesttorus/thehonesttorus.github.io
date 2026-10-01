#!/usr/bin/env bash
# Stage the public dataset into blob storage once, from inside Azure (Hugging Face is reachable
# there).  Runs as a single Batch task in the home region; afterwards every grid task reads its
# shard from blob in-region.  Optional: HF_TOKEN env var is not needed for this public dataset.
source "$(dirname "$0")/env.sh"
acct=$(batch_account_name "$HOME_REGION")
az batch account login -n "$acct" -g "$RG"
pool=$(pool_name "$HOME_REGION")
job="stage-dataset-$(date +%Y%m%d%H%M%S)"
az batch job create --id "$job" --pool-id "$pool" -o none
sas=$(az storage container generate-sas -n "$CONTAINER_DATASET" --account-name "$STORAGE" --permissions racwl --expiry "$(date -u -d '+2 days' +%Y-%m-%dT%H:%MZ)" --auth-mode login --as-user -o tsv)
cat > /tmp/stage_task.json <<JSON
[{
  "id": "stage",
  "commandLine": "python /app/stage_dataset.py --dataset $HF_DATASET --revision $HF_REVISION --dest 'https://$STORAGE.blob.core.windows.net/$CONTAINER_DATASET?$sas'",
  "containerSettings": {"imageName": "$ACR.azurecr.io/$IMAGE"},
  "constraints": {"maxWallClockTime": "PT12H", "maxTaskRetryCount": 2}
}]
JSON
az batch task create --job-id "$job" --json-file /tmp/stage_task.json -o none
echo "staging job $job submitted in pool $pool; watch with: az batch task show --job-id $job --task-id stage"
