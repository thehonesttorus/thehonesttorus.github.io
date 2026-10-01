#!/bin/sh
# streaming pair ladders at small widths on the same seed MLPs (two replicas each), for the width trend
# usage: sh run_widths.sh OUTDIR   (results tables written to results/)
O=$1; PY=${PY:-/root/whest/bin/python}; export OPENBLAS_NUM_THREADS=1
mkdir -p $O
for w in 128 256; do
  N=500000; [ $w = 256 ] && N=1000000
  for s in 770000 770001; do
    for r in 1 2; do
      [ -f $O/w${w}_s${s}_r$r.npz ] || $PY stream_oracle.py stats --seed $s --width $w --n-samples $N --chunk 4096 --sample-seed $((3000 + r)) --out $O/w${w}_s${s}_r$r.npz
    done
    $PY stream_oracle.py ladder $O/w${w}_s${s}_r1.npz $O/w${w}_s${s}_r2.npz --out results/width${w}_mlp${s}_pair.json | grep -v "layer" | sed "s#$O/##g" > results/width${w}_mlp${s}_pair.txt
  done
done
