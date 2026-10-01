import sys; sys.path.insert(0,'../../bench'); import bench, bethe, numpy as np, warnings
warnings.filterwarnings('ignore')
S=bench.load_set(sys.argv[1]); idx=[int(x) for x in sys.argv[2].split(',')]
cfg={'edge':bethe.estimate_edge,'v4':bethe.estimate_v4,'v4half':lambda W: bethe.estimate_v4(W,qscale=0.5),'v4nogen':lambda W: bethe.estimate_v4(W,qgen=False)}
print('mlp', ' '.join(cfg))
res={k:[] for k in cfg}
for i in idx:
  W=bench.weights(S,i).astype(np.float64); T=S['means'][i]
  for k,f in cfg.items(): res[k].append(((f(W)[-1]-T[-1])**2).mean())
  print(i, ' '.join(f'{res[k][-1]:.2e}' for k in cfg), flush=True)
print('mean', ' '.join(f'{np.mean(v):.2e}' for v in res.values()))
