#!/usr/bin/env bash
# Offline dry run of the whole grid against a fake `az` (tests/bin/az -> fake_az.py): runs every
# numbered script and every grid.py subcommand with dummy values, then checks
#   - every az call against the REAL installed CLI's help (command exists, every flag exists),
#   - every pool/task JSON against the Batch SDK models the CLI deserializes into,
#   - every task command line against the runner scripts' argparse, and the autoscale formulas.
# Needs: a real az (REAL_AZ, default: the first az on PATH), python3; AZ_PYTHON = az's own
# python (for azure.batch.models; default: the python next to REAL_AZ); RUNNER_PYTHON = a python
# with numpy (default python3).  Exit code != 0 if any check fails.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd); top=$(dirname "$here")
REAL_AZ="${REAL_AZ:-$(command -v az)}"
REAL_AZ=$(readlink -f "$REAL_AZ")
AZ_PYTHON="${AZ_PYTHON:-$(dirname "$REAL_AZ")/python}"
[ -x "$AZ_PYTHON" ] || AZ_PYTHON=$(head -1 "$REAL_AZ" | sed -n 's/^#!//p' | awk '{print $1}')
RUNNER_PYTHON="${RUNNER_PYTHON:-python3}"
W="${WORKDIR:-$(mktemp -d)}"; mkdir -p "$W"
export FAKE_AZ_LOG="$W/calls.log" FAKE_AZ_CAPTURE="$W/capture" FAKE_AZ_PYTHON="$AZ_PYTHON"
export PATH="$here/bin:$PATH" OPENBLAS_NUM_THREADS=1
rm -rf "$FAKE_AZ_CAPTURE" "$FAKE_AZ_LOG"; : > "$FAKE_AZ_LOG"
log="$W/run.log"; : > "$log"
fails=0
step() {  # name, command...
  local name="$1"; shift
  if ( "$@" ) >>"$log" 2>&1; then echo "ok   $name"; else echo "FAIL $name (see $log)"; fails=$((fails+1)); fi
}
expect_fail() {
  local name="$1"; shift
  if ( "$@" ) >>"$log" 2>&1; then echo "FAIL $name (should have been rejected)"; fails=$((fails+1)); else echo "ok   $name (rejected as expected)"; fi
}

# bundles and a results tree for collect
mkdir -p "$W/bundles/est_v1.2" "$W/bundles/est-b" "$W/res"
echo "class Estimator: pass" > "$W/bundles/est_v1.2/estimator.py"
echo "class Estimator: pass" > "$W/bundles/est-b/estimator.py"
tar -czf "$W/bundles/packed.tar.gz" -C "$W/bundles/est-b" estimator.py
if [ -n "${RESULTS_SAMPLE:-}" ] && [ -d "$RESULTS_SAMPLE" ]; then
  for f in "$RESULTS_SAMPLE"/*.json; do mkdir -p "$W/res/exp1/est_v1.2" "$W/res/exp1/packed"; cp "$f" "$W/res/exp1/est_v1.2/"; cp "$f" "$W/res/exp1/packed/"; done
else
  "$RUNNER_PYTHON" - "$W/res" <<'PY'
import json, os, sys
for tag in ("est_v1.2", "packed"):
    d = os.path.join(sys.argv[1], "exp1", tag); os.makedirs(d, exist_ok=True)
    for s in range(2):
        per = [{"adjusted_final_layer_score": 1e-6 * (s + 1), "final_layer_mse": 2e-6, "all_layers_mse": 3e-6,
                "flops_used": 2**40, "residual_wall_time_s": 0.1, "budget_exhausted": False, "traceback": None}] * 16
        json.dump({"task": {"harness_seconds": 10}, "report": {"results": {"per_mlp": per}}}, open(os.path.join(d, f"full-{s:05d}-of-00063.json"), "w"))
    json.dump({"task": {"harness_seconds": 1}, "report": None}, open(os.path.join(d, "full-00002-of-00063.json"), "w"))
PY
fi
export FAKE_AZ_RESULTS_SRC="$W/res"

cd "$top"
step "01_bootstrap (fresh)"           env FAKE_AZ_EXISTS=0 ./01_bootstrap.sh
step "01_bootstrap (re-run)"          env FAKE_AZ_EXISTS=1 ./01_bootstrap.sh
step "00_quotas"                      env FAKE_AZ_EXISTS=1 ./00_quotas.sh
step "02_build_image"                 ./02_build_image.sh
step "04_pools (create)"              env FAKE_AZ_POOL_EXISTS=0 ./04_pools.sh
step "04_pools (create, spot)"        env FAKE_AZ_POOL_EXISTS=0 USE_SPOT=1 REGIONS=westus2 ./04_pools.sh
step "04_pools (update)"              env FAKE_AZ_POOL_EXISTS=1 ./04_pools.sh
step "05_gpu_pool (create)"           env FAKE_AZ_POOL_EXISTS=0 ./05_gpu_pool.sh
step "05_gpu_pool (update)"           env FAKE_AZ_POOL_EXISTS=1 ./05_gpu_pool.sh
step "03_stage_dataset"               ./03_stage_dataset.sh
step "03_stage_dataset (community)"   env HF_DATASET=keenanpepper/arc-whestbench-p2-d8b-corpus-14k HF_REVISION=main STAGE_ARGS="--all-files --prefix d8b" ./03_stage_dataset.sh
expect_fail "03_stage_dataset (HOME_REGION not in REGIONS)" env REGIONS=westus2 ./03_stage_dataset.sh
expect_fail "env.sh (bad PREFIX)"     env PREFIX=Bad-Prefix ./00_quotas.sh
cd "$W"
step "grid submit (2 bundles x shards 0-3, extra)" "$RUNNER_PYTHON" "$top/grid.py" submit --name exp1 --bundle "$W/bundles/est_v1.2" --bundle "$W/bundles/packed.tar.gz" --split full --shards 0-3 --n-mlps 4 --extra="--fail-fast --debug"
step "grid submit (single shard, mini, subprocess, TaskExists)" env FAKE_AZ_TASK_EXISTS=1 "$RUNNER_PYTHON" "$top/grid.py" submit --name exp2 --bundle "$W/bundles/est-b" --split mini --shards 5 --runner subprocess --regions "eastus westus2"
step "grid status --failures"          "$RUNNER_PYTHON" "$top/grid.py" status --name exp1 --failures
step "grid collect --csv"              "$RUNNER_PYTHON" "$top/grid.py" collect --name exp1 --csv "$W/exp1.csv"
step "grid bake"                       "$RUNNER_PYTHON" "$top/grid.py" bake --name fresh-A --n-mlps 32 --n-samples 1000000 --slices 4
step "grid atlas --pairs"              "$RUNNER_PYTHON" "$top/grid.py" atlas --name d8b --start 10 --count 3 --pairs
expect_fail "grid submit (bad --name)" "$RUNNER_PYTHON" "$top/grid.py" submit --name "bad name!" --bundle "$W/bundles/est-b"
expect_fail "grid submit (duplicate bundle tags)" "$RUNNER_PYTHON" "$top/grid.py" submit --name exp3 --bundle "$W/bundles/est-b" --bundle "$W/bundles/est-b/"
expect_fail "grid submit (bad --runner)" "$RUNNER_PYTHON" "$top/grid.py" submit --name exp4 --bundle "$W/bundles/est-b" --runner docker
cat "$W/exp1.csv" >> "$log" 2>/dev/null || true

echo "--- az calls vs. real CLI help ($REAL_AZ)"
"$RUNNER_PYTHON" "$here/check_az_calls.py" "$FAKE_AZ_LOG" --az "$REAL_AZ" --cache "${AZHELP_CACHE:-$W/azhelp}" > "$W/az_calls.txt" || fails=$((fails+1))
grep -E '^FAIL|calls,' "$W/az_calls.txt"
echo "--- pool/task JSON vs. Batch SDK models ($AZ_PYTHON)"
"$AZ_PYTHON" "$here/validate_batch_json.py" "$FAKE_AZ_CAPTURE" || fails=$((fails+1))
echo "--- task command lines vs. runner argparse; autoscale formulas"
"$RUNNER_PYTHON" "$here/check_payloads.py" "$FAKE_AZ_CAPTURE" "$top/runner" || fails=$((fails+1))
echo "work dir: $W  (run.log, calls.log, az_calls.txt, capture/)"
if [ "$fails" -eq 0 ]; then echo "ALL OFFLINE CHECKS PASSED"; else echo "$fails check group(s) failed"; exit 1; fi
