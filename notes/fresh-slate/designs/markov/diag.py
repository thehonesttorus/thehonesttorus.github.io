import sys, numpy as np
sys.path.insert(0, '.')
from mkv import mkv
from bake import weights
TR = '/tmp/claude-0/-home-user-thehonesttorus-github-io/fc3b2401-9ebc-5116-8952-39e8e61f85e7/scratchpad/truth'
def truth_cums(t):
    m1, m2, m3, m4 = [t['zmom'][:, k] for k in range(4)]
    v = m2 - m1**2; k3 = m3 - 3*m2*m1 + 2*m1**3
    k4 = m4 - 4*m3*m1 - 3*m2**2 + 12*m2*m1**2 - 6*m1**4
    return m1, v, k3, k4
def report(width, seed, **kw):
    t = np.load(f'{TR}/mom_w{width}_s{seed}.npz'); W = weights(width, 16, seed)
    mu, v, k3, k4 = truth_cums(t)
    _, dg = mkv(W, return_all=True, **kw)
    print("layer  relerr_mu  relerr_var  err_k3/rms  err_k4/rms   (rms k3, rms k4)")
    for l in range(16):
        d = dg[l]
        r = lambda a, b: np.sqrt(np.mean((a-b)**2))/np.sqrt(np.mean(b**2))
        print(f"{l+1:3d} {r(d['mu'],mu[l]):.2e} {r(d['v'],v[l]):.2e} {r(d['k3'],k3[l]):.2e} {r(d['k4'],k4[l]):.2e}   ({np.sqrt(np.mean(k3[l]**2)):.2e}, {np.sqrt(np.mean(k4[l]**2)):.2e})")
if __name__ == '__main__':
    report(int(sys.argv[1]), int(sys.argv[2]), w=int(sys.argv[3]))
