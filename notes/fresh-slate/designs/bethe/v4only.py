import sys; sys.path.insert(0,'../../bench'); import bench, bethe, numpy as np, warnings, json
warnings.filterwarnings('ignore')
S=bench.load_set(sys.argv[1])
for i in [int(x) for x in sys.argv[2].split(',')]:
    W=bench.weights(S,i).astype(np.float64); T=S['means'][i]
    e=bethe.estimate_v4(W, old=int(sys.argv[3]) if len(sys.argv)>3 else 1)
    lay=((e-T)**2).mean(1)
    print(i, f'{lay[-1]:.4e}', 'noise', f"{S['noise'][i]:.2e}", 'per-layer', ' '.join(f'{x:.1e}' for x in lay), flush=True)
