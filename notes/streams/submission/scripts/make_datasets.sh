#!/bin/bash
# Rebuild every dataset this stream used (not committed: 385 MB + ~0.7 GB).
#   dev6: 6 MLPs, 1024 x 16, N = 2e6 (35 min on 3 threads); dev1/dev2/dev3 = its first 1/2/3 rows
#   rob_*: 1 MLP each at off-suite shapes (seed 4242, N = 2000) + adversarial edits of rob_w1024_d16
set -e
S=${1:-/tmp/claude-0/sub}; mkdir -p $S; cd $S
export PATH=/root/whest/bin:$PATH WHEST_SKIP_HARDWARE_FALLBACK_PROBES=1
HERE=$(cd "$(dirname "$0")" && pwd)
echo '[7301001,7301002,7301003,7301004,7301005,7301006]' > s6.json; echo '[4242]' > s1b.json
OMP_NUM_THREADS=3 whest dataset bake --n-mlps 6 --n-samples 2000000 --width 1024 --depth 16 --mlp-seeds s6.json --output dev6
for k in 1 2 3; do /root/whest/bin/python - $k <<'PY'
import sys, pyarrow.parquet as pq, glob, os, json
k = int(sys.argv[1]); src = 'dev6'; dst = f'dev{k}'
os.makedirs(dst + '/data', exist_ok=True); m = json.load(open(src + '/metadata.json')); m['n_mlps'] = k
json.dump(m, open(dst + '/metadata.json', 'w'))
f = sorted(glob.glob(src + '/data/*.parquet'))[0]; pq.write_table(pq.read_table(f).slice(0, k), dst + '/data/' + os.path.basename(f))
PY
done
for s in "1024 32" "512 16" "1024 4" "256 8" "1024 16"; do set -- $s
  whest dataset bake --n-mlps 1 --n-samples 2000 --width $1 --depth $2 --mlp-seeds s1b.json --output rob_w$1_d$2; done
for m in x10 abs; do /root/whest/bin/python $HERE/make_adversarial.py rob_w1024_d16 rob_adv_$m $m; done
