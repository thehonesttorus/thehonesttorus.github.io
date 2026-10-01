#!/bin/bash
# Idle-machine validation on the dev set (1024 x 16, 6 MLPs, N = 2e6). Run with NOTHING else on the box.
# usage: validate_dev.sh OUTDIR
OUT=$1; mkdir -p $OUT
export PATH=/root/whest/bin:$PATH
HERE=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-0/sub; DS=$S/dev6
M="/root/whest/bin/python $HERE/scripts/measure_run.py --dataset $DS --max-threads 2"
for rep in 1 2 3; do
  for v in v29 v25; do
    $M --estimator $HERE/bundles/$v/estimator.py --runner subprocess --tag ${v}_sub_r$rep --out $OUT/${v}_sub_r$rep.json
  done
done
for v in v29 v25; do
  $M --estimator $S/orig/$v/estimator.py --runner subprocess --tag ${v}orig_sub_r1 --out $OUT/${v}orig_sub_r1.json
  $M --estimator $HERE/bundles/$v/estimator.py --runner local --tag ${v}_local_r1 --out $OUT/${v}_local_r1.json
done
for v in v29 v25; do
  OMP_NUM_THREADS=2 /root/whest/bin/python $HERE/scripts/setup_time.py $HERE/bundles/$v/estimator.py 5 --out $OUT/${v}_setup_time.json
done
