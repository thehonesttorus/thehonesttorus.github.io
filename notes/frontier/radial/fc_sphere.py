"""Radial (homogeneity) factorisation test in region's FC (no fitted constants).
Bias-free ReLU MLPs are positively homogeneous: F(x) = s F(u), x = r v, s = r/sqrt(n), u = sqrt(n) v uniform on the
sqrt(n)-sphere, r independent of v.  Hence E_x F = E[s] E_u F exactly.  The u-problem has no radial (dilation) mode;
its input is non-Gaussian only through the sphere's fourth cumulant kappa4(u)_ijkl = c (d_ij d_kl + d_ik d_jl + d_il d_jk),
c = -2/(n+2), so the layer-0 pre-activations z = W^T u carry kappa4(z_a,z_a,z_b,z_b) = c (S_aa S_bb + 2 S_ab^2),
kappa4(z_a) = 3 c S_aa^2, kappa4(z_a,z_a,z_a,z_b) = 3 c S_aa S_ab (S = W^T W).  FC runs the u-problem with that seed in
its mean-field kappa4 channel, and the outputs are multiplied by E[s].
usage: python3 fc_sphere.py MLPS [window]"""
import json, os, sys, time
import numpy as np
from scipy.special import gammaln
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '../../fresh-slate/bench'))
sys.path.insert(0, os.path.join(HERE, '../../fresh-slate/breakthrough/region'))
import bench, fc

def es(n):
    return float(np.exp(0.5 * np.log(2.0 / n) + gammaln((n + 1) / 2.0) - gammaln(n / 2.0)))

def seed(W0):
    n = W0.shape[0]
    S = W0.T.astype(np.float64) @ W0.astype(np.float64)
    c = -2.0 / (n + 2)
    d = np.diag(S)
    K22 = c * (np.outer(d, d) + 2 * S * S)
    c22 = (K22.sum(0) - np.diag(K22)) / (n - 1)            # column means over a != b (mean-field form)
    K31 = 3 * c * d[:, None] * S; np.fill_diagonal(K31, 0)
    k4 = 3 * c * d * d
    return {0: (np.broadcast_to(c22[None, :], (n, n)).copy(), K31, k4)}

if __name__ == '__main__':
    S_ = bench.load_set('w1024_d16')
    mlps = [int(a) for a in sys.argv[1].split(',')]
    window = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] != 'all' else None
    for i in mlps:
        W = bench.weights(S_, i); T = S_['means'][i]; nz = float(S_['noise'][i]); n = W.shape[1]
        out = dict(mlp=i, window=window)
        t0 = time.time()
        px = fc.run(W, slices=2, k4mf=True, window=window)
        pu = fc.run(W, slices=2, k4mf=True, window=window, k4f=seed(W[0])) * es(n)
        pu0 = fc.run(W, slices=2, k4mf=True, window=window) * es(n)     # radial factor alone, no sphere seed
        for k, p in (('x', px), ('u', pu), ('u_noseed', pu0)):
            out[k] = float(((p[-1] - T[-1]) ** 2).mean() - nz)
        out['sec'] = round(time.time() - t0, 1)
        print(json.dumps(out), flush=True)
