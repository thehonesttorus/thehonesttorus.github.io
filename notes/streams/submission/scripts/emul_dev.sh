#!/bin/bash
# Client/server (grader-transport) emulation on the dev set; run on an idle box.
# usage: emul_dev.sh OUTDIR REPS LABEL:ESTIMATOR_PY [...]
OUT=$1; REPS=$2; shift 2; mkdir -p $OUT
HERE=$(cd "$(dirname "$0")/.." && pwd)
for rep in $(seq 1 $REPS); do
  for spec in "$@"; do
    lab=${spec%%:*}; est=${spec#*:}
    /root/whest/bin/python $HERE/scripts/grader_emul.py run --estimator $est --dataset /tmp/claude-0/sub/dev6 \
      --server-threads ${SERVER_THREADS:-3} --tag ${lab}_emul_r$rep --out $OUT/${lab}_emul_r$rep.json
  done
done
