"""layer-3 decomposition: K21 of a_2 (MC vs v1/v0), and kappa3(z_3) sector split."""
import numpy as np, stageq, copula1 as c1, copula as c0
W, mean, noise, _ = stageq.load(64, 0); n = 64; off = ~np.eye(n, dtype=bool)
rng = np.random.default_rng(9); N = 0; Sa = np.zeros(n); E2 = np.zeros((n,n)); E21 = np.zeros((n,n)); E3 = np.zeros(n); Sz=np.zeros((4,n))
mu2 = mean[1]
while N < 3e6:
    x = rng.standard_normal((1<<17, n)); a = np.maximum(np.maximum(x@W[0],0)@W[1],0); ac = a - mu2
    z = a@W[2] - mu2@W[2]
    E2 += ac.T@ac; E21 += (ac**2).T@ac; E3 += (ac**3).sum(0)
    for k in range(4): Sz[k] += (z**(k+1)).sum(0)
    N += 1<<17
E2/=N; E21/=N; E3/=N
w = W[2]
r = lambda x: np.sqrt(np.mean(x**2))
# MC sector values of kappa3(z_3)
s3 = E3@w**3; s21 = 3*np.sum(w**2*((E21*off)@w),0)
m = Sz/N; k3 = m[2]-3*m[0]*m[1]+2*m[0]**3
s111 = k3 - s3 - s21
for name, mod in [('v0', c0), ('v1', c1)]:
    est, dg = mod.estimate(W, {}, return_diag=True)
    # recompute this layer's K21 from the state entering layer 2->3: re-run step to grab internals
    print(name, "k3(z3) err", r(dg[1]['k3']-k3), "of", r(k3))
print("MC sectors rms: {3}", r(s3), "{2,1}", r(s21), "{1,1,1} (residual)", r(s111))
# v1 sector values
est, dg = c1.estimate(W, {'src':1}, return_diag=True); dg[1]['s111'] = dg[1]['k3'] - dg[1]['s3'] - dg[1]['s21']
d = dg[1]
print("v1 sector errs: {3}", r(d['s3']-s3), "{2,1}", r(d['s21']-s21), "{1,1,1}", r(d['s111']-s111))
print("K21(a_2) offdiag err", r((d['K21']-E21)[off]), "of", r(E21[off]))
