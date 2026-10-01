#!/bin/sh
# pair ladders at width 1024 on the production replicas (two per MLP), one MLP per process, 2 BLAS threads each
# usage: sh run_prod_ladders.sh DIR   (DIR holds s<seed>_r{1,2}.npz from `stream_oracle.py stats`)
D=$1; PY=${PY:-/root/whest/bin/python}; export OPENBLAS_NUM_THREADS=2
for s in 770000 770001; do
  ( $PY stream_oracle.py ladder $D/s${s}_r1.npz $D/s${s}_r2.npz --gauss-cache $D/gauss_$s.npz \
      --out results/width1024_mlp${s}_pair.json > $D/ladder_$s.log 2>&1
    grep -v "layer" $D/ladder_$s.log | sed "s#$D/##g" > results/width1024_mlp${s}_pair.txt ) &
done
wait
