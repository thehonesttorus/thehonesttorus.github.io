"""Summarise every sw_<cfg>.json in results/ (robust to the sweep driver dying): per-MLP raw, mean raw, steady C/B
(max over the non-first MLPs of the run), adjusted = mean raw x max(0.1, steady C/B).  usage: python3 summarize.py [prefix]"""
import glob, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
pre = sys.argv[1] if len(sys.argv) > 1 else ''
for f in sorted(glob.glob(os.path.join(HERE, 'results', f'sw_{pre}*.json'))):
    try:
        rows = json.load(open(f))['per_mlp']
    except Exception:
        continue
    ok = [r for r in rows if r.get('ok')]
    if not ok:
        print(os.path.basename(f), 'no ok rows'); continue
    raw = sum(r['raw'] for r in ok) / len(ok)
    cbs = [r['cb'] for r in ok]
    cb = max(cbs[1:]) if len(cbs) > 1 else cbs[0]
    name = os.path.basename(f)[3:-5]
    mlps = ','.join(str(r['mlp']) for r in ok)
    per = ' '.join('%.3f' % (r['raw'] * 1e8) for r in ok)
    print('%-24s mlps %-12s raw %.4e [%s] cb %.4f adj %.4e' % (name, mlps, raw, per, cb, raw * max(0.1, cb)))
