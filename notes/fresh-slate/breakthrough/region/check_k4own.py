import sys, numpy as np
sys.path.insert(0, '../../bench')
import bench, fc, mc_analyse as ma
i = int(sys.argv[1]) if len(sys.argv) > 1 else 0
S = bench.load_set('w1024_d16'); W = bench.weights(S, i)
D = np.load(f'data/mc1024_mlp{i}_N262144.npz')
k4f = {}
tr = []; p = fc.run(W, trace=tr, slices=2, k4f=k4f, k4own=(sys.argv[2] == 'own'), k4mf=(sys.argv[2] == 'mf'))
print('raw', ((p[-1] - S['means'][i][-1]) ** 2).mean() - S['noise'][i])
for l in range(1, 16):
    cA = ma.central(D, 'A', l); cB = ma.central(D, 'B', l)
    K22 = 0.5 * (cA['K22'] + cB['K22']); np.fill_diagonal(K22, 0)
    cm = K22.mean(0); k4 = 0.5 * (cA['k4'] + cB['k4'])
    c22m, _, k4m = k4f[l]
    c22m = c22m[0]
    r1 = np.linalg.norm(c22m - cm) / np.linalg.norm(cm); r2 = np.linalg.norm(k4m - k4) / np.linalg.norm(k4)
    print(l, f"colmean K22: rel err {r1:.3f} (|MC| {np.linalg.norm(cm):.2e}, ratio {np.dot(c22m, cm) / np.dot(cm, cm):.2f});  k4 diag: rel err {r2:.3f} ratio {np.dot(k4m, k4) / np.dot(k4, k4):.2f}", flush=True)
