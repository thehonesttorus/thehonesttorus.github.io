import sys, numpy as np, copula1
w, seed = int(sys.argv[1]), int(sys.argv[2])
W = np.load(f"/root/sg/truth/W_w{w}_s{seed}.npy").astype(float)
d = np.load(f"/root/sg/truth/w{w}_s{seed}.npz"); N = int(d['n']); mean = d['S'] / N; noise = (d['Q'] / N - mean ** 2) / N
O = dict(np.load(f"/root/sg/layerstats_w{w}_s{seed}.npz"))
import ast
for use in (ast.literal_eval(sys.argv[3]) if len(sys.argv) > 3 else [(), ('D',), ('marg',), ('marg', 'D'), ('cov',), ('cov', 'D'), ('marg', 'cov'), ('marg', 'cov', 'D')]):
    cfg = {'oracle': dict(O, use=use)} if use else {}
    est = copula1.estimate(W, cfg)
    err = ((est - mean) ** 2).mean(1) - noise.mean(1)
    print(f"{'+'.join(use) or 'v1':14s} final {err[-1]:.2e}  L4 {err[3]:.2e} L8 {err[7]:.2e} L12 {err[11]:.2e}", flush=True)
