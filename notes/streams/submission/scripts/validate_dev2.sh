#!/bin/bash
# Idle-machine validation, pass 2 (after finding that V29's in-process RSS breaks the subprocess
# runner's 8 GB RLIMIT_AS): V29 via the local runner (no AS limit; arrays live in the flopscope
# server on the grader), V25 via subprocess + local. Run with NOTHING else on the box.
OUT=$1; mkdir -p $OUT
export PATH=/root/whest/bin:$PATH
HERE=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-0/sub; DS=$S/dev6
M="/root/whest/bin/python $HERE/scripts/measure_run.py --dataset $DS --max-threads 2"
for rep in 2 3; do
  $M --estimator $HERE/bundles/v25/estimator.py --runner subprocess --tag v25_sub_r$rep --out $OUT/v25_sub_r$rep.json
done
for rep in 1 2 3; do
  $M --estimator $HERE/bundles/v29/estimator.py --runner local --tag v29_local_r$rep --out $OUT/v29_local_r$rep.json
done
$M --estimator $S/orig/v29/estimator.py --runner local --tag v29orig_local_r1 --out $OUT/v29orig_local_r1.json
$M --estimator $S/orig/v25/estimator.py --runner subprocess --tag v25orig_sub_r1 --out $OUT/v25orig_sub_r1.json
$M --estimator $HERE/bundles/v25/estimator.py --runner local --tag v25_local_r1 --out $OUT/v25_local_r1.json
for v in v29 v25; do
  /root/whest/bin/python $HERE/scripts/setup_time.py $HERE/bundles/$v/estimator.py 5 --out $OUT/${v}_setup_time.json
done
