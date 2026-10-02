import json,sys,os
for f in sys.argv[1:]:
    d=json.load(open(f)); rows=[r for r in d['per_mlp'] if r.get('ok')]
    raw=sum(r['raw'] for r in rows)/len(rows); cbs=[r['cb'] for r in rows]; cb=max(cbs[1:]) if len(cbs)>1 else cbs[0]
    print(f"{os.path.basename(f)[3:-5]:24s} n={len(rows)} raws={[round(r['raw']*1e8,4) for r in rows]} raw={raw:.4e} cb={cb:.4f} adj={raw*max(.1,cb):.4e}")
