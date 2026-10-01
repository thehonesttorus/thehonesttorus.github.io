#!/usr/bin/env bash
# Stage the public dataset into blob storage once, from inside Azure (Hugging Face is reachable
# there).  Runs as a single Batch task on the HOME_REGION CPU pool (04_pools.sh first);
# afterwards every grid task reads its shard from blob.  The dataset is public: no HF token.
# Overrides: HF_DATASET / HF_REVISION (another repo), STAGE_ARGS (extra stage_dataset.py flags,
# e.g. STAGE_ARGS="--splits mini" or "--all-files --prefix d8b" for a community atlas).
source "$(dirname "$0")/env.sh"
case " $REGIONS " in *" $HOME_REGION "*) ;; *) echo "HOME_REGION=$HOME_REGION must be listed in REGIONS (it needs a pool)" >&2; exit 1;; esac
batch_login "$HOME_REGION"
pool=$(pool_name "$HOME_REGION")
az batch pool show --pool-id "$pool" -o none || { echo "pool $pool not found: run 04_pools.sh first" >&2; exit 1; }
job="stage-dataset-$(date -u +%Y%m%d%H%M%S)"
key=$(storage_key)
sas=$(az storage container generate-sas -n "$CONTAINER_DATASET" --account-name "$STORAGE" --account-key "$key" \
      --permissions racwl --expiry "$(utc_in_hours 48)" -o tsv)
umask 077
task_json=$(mktemp); trap 'rm -f "$task_json" "$task_json.out"' EXIT
# Built in Python so the SAS (& % =) and the quoting survive: Batch does not run the command
# line under a shell, so it is wrapped in /bin/bash -c '...'.  Admin (root in the container)
# because Batch maps the task user into the container and a non-admin user cannot write outside
# the task directory.
SAS="$sas" IMG="$ACR.azurecr.io/$IMAGE" python3 - "$task_json" <<'PY'
import json, os, shlex, sys
e = os.environ
dest = f"https://{e['STORAGE']}.blob.core.windows.net/{e['CONTAINER_DATASET']}?{e['SAS']}"
cmd = ["python", "/app/stage_dataset.py", "--dataset", e["HF_DATASET"], "--revision", e["HF_REVISION"], "--dest", dest]
cmd += shlex.split(e.get("STAGE_ARGS", ""))
task = {"id": "stage",
        "commandLine": "/bin/bash -c " + shlex.quote(shlex.join(cmd)),
        "containerSettings": {"imageName": e["IMG"], "containerRunOptions": "--rm"},
        "userIdentity": {"autoUser": {"scope": "pool", "elevationLevel": "admin"}},
        "constraints": {"maxWallClockTime": "PT12H", "maxTaskRetryCount": 2}}
json.dump([task], open(sys.argv[1], "w"), indent=1)
PY
az batch job create --id "$job" --pool-id "$pool" -o none
az batch task create --job-id "$job" --json-file "$task_json" -o json > "$task_json.out"
check_task_results "$task_json.out"
# terminate the job once its task completes, so it stops counting against the active-job quota
az batch job set --job-id "$job" --on-all-tasks-complete terminatejob -o none
echo "staging job $job submitted in pool $pool; watch with: az batch task show --job-id $job --task-id stage --query '{state:state,exit:executionInfo.exitCode}'"
