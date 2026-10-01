import numpy as np, sys, chain as ch
from run_variants import load_truth
t=load_truth('/root/work/truth128')[int(sys.argv[1])]; W=np.asarray(t['weights'],dtype=np.float64); g=np.asarray(t['all_layer_means'])
for kw in eval(sys.argv[2]):
    m=ch.Chain(W,record=False,**kw).run()['means']
    print(kw, ' '.join(f'{x:.1e}' for x in np.mean((m-g)**2,1)), flush=True)
