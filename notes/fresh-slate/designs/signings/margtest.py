"""how well does a 4-cumulant (Fleishman) marginal reproduce E a, Var a, k3 a of a = ReLU(z), given exact z cumulants?"""
import numpy as np, copula, copula1
w, seed, t = 128, 0, 14
W = np.load(f"/root/sg/truth/W_w{w}_s{seed}.npy").astype(float); O = np.load(f"/root/sg/layerstats_w{w}_s{seed}.npz")
rng = np.random.default_rng(3); S = np.zeros((3, w)); N = 0
while N < 2e6:
    a = rng.standard_normal((1 << 15, w))
    for l in range(t + 1): a = np.maximum(a @ W[l], 0)
    S[0] += a.sum(0); S[1] += (a ** 2).sum(0); S[2] += (a ** 3).sum(0); N += 1 << 15
m = S / N; mu = m[0]; var = m[1] - mu ** 2; k3 = m[2] - 3 * mu * m[1] + 2 * mu ** 3
s = np.sqrt(O['var'][t]); g1 = O['k3'][t] / s ** 3; g2 = O['k4'][t] / s ** 4
a = copula.fleishman(np.clip(g1, -1.5, 1.5), np.clip(g2, -1, 4))
pr = copula1.ext_profiles(O['mz'][t], s, a, dmax=0)
h, F2, F3 = pr['h'][:, 0], pr['F2'][:, 0], pr['F3'][:, 0]
pvar = F2 - h ** 2; pk3 = F3 - 3 * h * F2 + 2 * h ** 3
r = lambda x: np.sqrt(np.mean(x ** 2))
print("std. skew rms", r(g1), "exkurt rms", r(g2))
print("mean a: rms err", r(h - mu), " var a rel", r(pvar / var - 1), " k3 a", r(pk3 - k3), "/", r(k3))
# Gaussian-marginal baseline
pr0 = copula1.ext_profiles(O['mz'][t], s, np.tile([1., 0, 0], (w, 1)), dmax=0)
print("gaussian marginal: mean err", r(pr0['h'][:, 0] - mu))
# covariance of a_14 by MC vs the v1 model (forced state)
rng = np.random.default_rng(4); C = np.zeros((w, w)); Sa = np.zeros(w); N = 0
while N < 2e6:
    a = rng.standard_normal((1 << 15, w))
    for l in range(t + 1): a = np.maximum(a @ W[l], 0)
    Sa += a.sum(0); C += a.T @ a; N += 1 << 15
ma = Sa / N; Ca = C / N - np.outer(ma, ma)
O2 = dict(O)
est, dg = copula1.estimate(W, {'oracle': dict(O2, use=('marg', 'cov', 'D')), 'keepCa': 1}, return_diag=True)
Cm = dg[t]['Ca']; off = ~np.eye(w, dtype=bool); Wn = W[t + 1]
print("Cov(a) offdiag err", r((Cm - Ca)[off]), "/", r(Ca[off]), " diag rel", r((np.diag(Cm) - np.diag(Ca))[np.diag(Ca) > 1e-6] / np.diag(Ca)[np.diag(Ca) > 1e-6]))
vd = np.einsum('aj,a->j', Wn ** 2, np.diag(Cm) - np.diag(Ca)); vo = np.einsum('aj,ab,bj->j', Wn, (Cm - Ca) * off, Wn)
print("var(z) err from diag", r(vd), " from offdiag", r(vo), " var z rms", r(O['var'][t + 1]))
# where is the offdiag error: copula part vs hyperedge part
print("offdiag err without hyperedge term", r((dg[t]['Ca0'] - Ca)[off]))
