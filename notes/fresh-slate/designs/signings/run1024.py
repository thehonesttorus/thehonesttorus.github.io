import sys, time, json, numpy as np
sys.path.insert(0, '../../bench'); import eval_q, preds
mlps = [int(x) for x in sys.argv[1].split(',')]
out = {}
for f in ['v1pq1', 'gauss']:
    res = eval_q.eval_set(getattr(preds, f), 'w1024_d16', mlps=mlps, verbose=True)
    out[f] = res
json.dump(out, open('/root/sg/w1024_%s.json' % sys.argv[1].replace(',', '_'), 'w'), default=str)
