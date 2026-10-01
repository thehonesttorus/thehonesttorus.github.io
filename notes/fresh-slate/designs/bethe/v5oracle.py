import sys; sys.path.insert(0,'../../bench'); import bench, bethe, numpy as np, warnings
warnings.filterwarnings('ignore')
S=bench.load_set(sys.argv[1]); N=int(float(sys.argv[3]))
for i in [int(x) for x in sys.argv[2].split(',')]:
    W=bench.weights(S,i).astype(np.float64); T=S['means'][i]; L,n,_=W.shape
    rng=np.random.default_rng(11); a=rng.standard_normal((N,n)); glaw=[]
    for l in range(L):
        z=a@W[l]; a=np.maximum(z,0); zc=z-z.mean(0); C=zc.T@zc/N
        u=np.linalg.eigh(C)[1][:,-1]; u=u if u.sum()>=0 else -u; g=zc@u; g/=g.std()
        M=np.array([np.mean(g**r) for r in range(13)])
        try: glaw.append(bethe.quad_from_moments(M,6))
        except Exception as ex: glaw.append(None)
    e0=bethe.estimate_v5(W); e1=bethe.estimate_v5(W,glaw=glaw); e2=bethe.estimate_edge(W)
    print(i, ' '.join(f'{((e[-1]-T[-1])**2).mean():.2e}' for e in (e2,e0,e1)), 'k3/k4 g at 15: %.2f %.2f'%(np.sum(glaw[15][1]*glaw[15][0]**3), np.sum(glaw[15][1]*glaw[15][0]**4)-3), flush=True)
