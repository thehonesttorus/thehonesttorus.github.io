"""Knob sweep of the public V29 chain (504aldo, MIT) on the dev set w1024_d16, paired on the same MLPs.
Each config = env overrides; runs run_p_inproc.py in a subprocess; a queue with W concurrent workers.
Reward: adjusted = mean raw x max(0.1, steady C/B) (steady = C/B of the non-first MLPs in a process).
usage: python3 sweep.py CONFIGS.json MLPS WORKERS THREADS"""
import json, os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.join(HERE, '../../fresh-slate/bench')
cfgs = json.load(open(sys.argv[1])); mlps = sys.argv[2]; W = int(sys.argv[3]); th = sys.argv[4]
os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
queue = list(cfgs.items()); running = []
def launch(name, cfg):
    env = dict(os.environ, OMP_NUM_THREADS=th, OPENBLAS_NUM_THREADS=th, MKL_NUM_THREADS=th)
    env.update({k: str(v) for k, v in cfg.get('env', {}).items()})
    est = os.path.join(HERE, '..', cfg.get('estimator', 'estimator_v29r3.py'))
    out = os.path.join(HERE, 'results', f'sw_{name}.json')
    log = open(os.path.join(HERE, 'results', f'sw_{name}.log'), 'w')
    p = subprocess.Popen(['/root/whest/bin/python', 'run_p_inproc.py', '--estimator', est, '--set', 'w1024_d16',
                          '--mlps', mlps, '--json', out], cwd=BENCH, env=env, stdout=log, stderr=subprocess.STDOUT)
    return (name, p, out, time.time())
while queue or running:
    while queue and len(running) < W:
        n_, c_ = queue.pop(0); running.append(launch(n_, c_))
    time.sleep(5)
    for r in list(running):
        name, p, out, t0 = r
        if p.poll() is not None:
            running.remove(r)
            try:
                rows = json.load(open(out))
                rows = rows['per_mlp'] if isinstance(rows, dict) else rows
                ok = [x for x in rows if x.get('ok')]
                raw = sum(x['raw'] for x in ok) / max(1, len(ok))
                cbs = [x['cb'] for x in ok]
                cb = max(cbs[1:]) if len(cbs) > 1 else cbs[0]
                print(json.dumps(dict(cfg=name, n_ok=len(ok), raw=raw, raws=[round(x['raw'] * 1e8, 4) for x in ok],
                                      cb=round(cb, 4), adjusted=raw * max(0.1, cb), min=round((time.time() - t0) / 60, 1))), flush=True)
            except Exception as e:
                print(json.dumps(dict(cfg=name, error=str(e)[:200])), flush=True)
print('SWEEP_DONE', flush=True)
