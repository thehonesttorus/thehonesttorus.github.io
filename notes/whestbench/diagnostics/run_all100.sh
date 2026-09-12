set -e
SP=$(ls -d /tmp/claude-0/-home-user-thehonesttorus-github-io/*/scratchpad | head -1)
cd $SP/est
for R in 16-44 44-72 72-100; do
  echo "=== building $R ==="
  python3 build_split.py $R
  echo "=== running $R ==="
  whest run --estimator estimator.py --dataset ../ds_${R/-/_} --split mini --runner local --json \
      > run_${R}.json 2> run_${R}.err || echo "run $R failed"
  rm -rf ../ds_${R/-/_}
  python3 - <<PY
import json
d=json.load(open("run_${R}.json"))
if 'results' in d:
    r=d['results']
    print("  $R : n=%d  adjusted %.6e  mse %.6e  util %.5f  failed %d" % (
        d['run_config']['n_mlps'], r['adjusted_final_layer_score'], r['final_layer_mse'],
        r['mean_compute_utilization'], r['n_failed_mlps']))
else: print("  $R : ", d.get('error'))
PY
done
echo ALLDONE
