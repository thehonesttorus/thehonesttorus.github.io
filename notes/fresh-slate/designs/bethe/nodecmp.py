import sys, numpy as np
sys.path.insert(0, '../../bench'); sys.path.insert(0, '.')
import bench, bethe, oracle
S = bench.load_set(sys.argv[1]); i = int(sys.argv[2]); N = float(sys.argv[3])
W = bench.weights(S, i).astype(np.float64); Tm = S['means'][i]
tr = oracle.mc_state(W, int(N))
rel = lambda a, b: np.sqrt(np.mean((a - b) ** 2)) / np.sqrt(np.mean(b ** 2))
for name, f in (('edge', bethe.estimate_edge_info), ('v4', lambda W: bethe.estimate_v4(W, verbose=True))):
    e, info = f(W)
    print(name, 'layer: rel err m v k3 k4 | mse(mu)')
    for l in range(len(info)):
        t = tr[l]; d = info[l]
        print(f"  {l:2d}: {rel(d['m'],t['m']):.4f} {rel(d['v'],np.diag(t['C'])):.4f} {rel(d['k3'],np.diag(t['K'])):.3f} {rel(d['k4'],t['k4']):.3f} | {((e[l]-Tm[l])**2).mean():.2e}", flush=True)
