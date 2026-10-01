#!/bin/bash
# Robustness sweep: estimator x {off-suite shapes, adversarial weights}, subprocess runner (RLIMIT_AS 8 GB).
# usage: robustness.sh OUTDIR THREADS WALL LABEL:ESTIMATOR_PY [LABEL:ESTIMATOR_PY ...]
OUT=$1; TH=$2; WALL=$3; shift 3
export PATH=/root/whest/bin:$PATH
S=/tmp/claude-0/sub
HERE=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p $OUT
for spec in "$@"; do
  lab=${spec%%:*}; est=${spec#*:}
  for d in ${SETS:-rob_w1024_d32 rob_w512_d16 rob_w1024_d4 rob_w256_d8 rob_adv_x10 rob_adv_abs rob_w1024_d16}; do
    /root/whest/bin/python $HERE/scripts/measure_run.py --estimator $est --dataset $S/$d \
      --runner subprocess --max-threads $TH --wall-time-limit $WALL --tag ${lab}_${d} --out $OUT/${lab}_${d}.json
  done
done
