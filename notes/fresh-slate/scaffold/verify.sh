#!/bin/bash
# End-to-end grader-cap check of an estimator: bakes (once) a small 1024x16 set and the 256x32 smoke
# shape, runs `whest run` in an isolated subprocess with the graded caps, prints one line per set.
# usage: verify.sh ESTIMATOR_PY [N_MLPS=2] [DATA_DIR=/tmp/claude-0/scaffold-data]
# Truth here is N = 20000 samples (noise floor ~4e-6): this checks failures, caps, cost and residual,
# not accuracy. For accuracy use high-N truth (the shared bench under notes/fresh-slate/bench/).
set -e
EST=$(readlink -f "$1"); NM=${2:-2}; D=${3:-/tmp/claude-0/scaffold-data}
export PATH=/root/whest/bin:$PATH
mkdir -p $D
bake() { [ -d $D/$1 ] || whest dataset bake --n-mlps $2 --n-samples $3 --width $4 --depth $5 --output $D/$1 > $D/$1.bake.log 2>&1; }
bake suite1024x16 4 20000 1024 16
bake smoke256x32 1 2000 256 32
bake w512x16 1 2000 512 16
HERE=$(dirname "$(readlink -f "$0")")
for s in smoke256x32:1 w512x16:1 suite1024x16:$NM; do
  set_=${s%%:*}; n=${s#*:}
  out=$D/run_$(basename $(dirname $EST))_$set_.json
  whest run --estimator $EST --runner subprocess --dataset $D/$set_ --n-mlps $n --format json \
    --max-threads 2 --wall-time-limit 120 --residual-wall-time-limit 0.4 --setup-timeout 5 > $out 2> $out.err || true
  python3 $HERE/run_summary.py $out $set_
done
