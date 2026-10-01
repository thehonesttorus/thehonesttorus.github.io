import sys; sys.path.insert(0,'../../bench'); import bench, bethe, numpy as np, warnings, time
warnings.filterwarnings('ignore')
S=bench.load_set('w1024_d16'); mode=sys.argv[1]; idx=[int(x) for x in sys.argv[2].split(',')]
for i in idx:
    W=bench.weights(S,i).astype(np.float64); T=S['means'][i]; nz=S['noise'][i]
    t0=time.time()
    if mode=='gauss':
        es=[bethe.estimate_gauss(W), bethe.localized(W, bethe.estimate_gauss, K=7)]
    else:
        es=[bethe.localized(W, bethe.estimate_v4, K=int(sys.argv[3]))]
    print(mode, i, ' '.join(f'{((e[-1]-T[-1])**2).mean()-nz:.3e}' for e in es), '%.0fs'%(time.time()-t0), flush=True)
