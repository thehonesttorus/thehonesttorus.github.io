# Structure of the network gain over inputs: log-gain increments across layers, cylinder share,
# and whether the Gaussian closure equals the gain-tilted ("beta = 2") mean.
import numpy as np, sys
from closure import closure
from cylinders import frame
n, L, s, T = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(float(sys.argv[4]))
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz")
oc = closure(Ws); Wf = [W.T.astype(np.float32) for W in Ws]
P = frame(Ws, "jac", 6).astype(np.float32)
rng = np.random.default_rng(3); B = 10000
logg = []; cyl = []; A = np.zeros((L, n)); Bv = np.zeros((L, n)); Nrm2 = np.zeros(L); Nrm2u = np.zeros(L); Hm = np.zeros((L, n)); Hu = np.zeros((L, n))
for _ in range(T//B):
    x = rng.standard_normal((B, n)).astype(np.float32); r2 = (x.astype(np.float64)**2).sum(1)
    cyl.append(((x @ P.T) > 0) @ (2**np.arange(6)))
    h = x; lg = []
    for l in range(L):
        h = np.maximum(h @ Wf[l], 0); hd = h.astype(np.float64); nn = (hd*hd).sum(1)
        lg.append(np.log(nn/r2))
        A[l] += (hd*np.sqrt(nn)[:, None]).sum(0); Nrm2[l] += nn.sum()
        Bv[l] += (hd*np.sqrt(nn/r2)[:, None]).sum(0); Nrm2u[l] += (nn/r2).sum()
        Hm[l] += hd.sum(0)
    logg.append(np.array(lg).T)
logg = np.concatenate(logg); cyl = np.concatenate(cyl)
A /= T; Bv /= T; Nrm2 /= T; Nrm2u /= T; Hm /= T
print(f"n={n} L={L} s={s}, T={T}")
print("layer  oracle_c   Var(log g)  1-E[sqrt g]/sqrt(E g)   cyl-share(Var log g)   |  KMS test: rel.dist(closure, tilted)  rel.dist(closure, truth)")
for l in [1, 3, 7, 11, 15]:
    if l >= L: continue
    m = tr["m"][l]; mc = oc[l]["m"]; c = mc @ (mc-m)/(mc @ mc)
    g = np.exp(logg[:, l]); vb = 1 - np.mean(np.sqrt(g))/np.sqrt(np.mean(g))
    tot = logg[:, l].var(); means = np.bincount(cyl, weights=logg[:, l], minlength=64)/np.maximum(np.bincount(cyl, minlength=64), 1)
    share = np.var(means[cyl])/tot
    tiltA = A[l]/np.sqrt(Nrm2[l])          # E[h ||h||]/sqrt(E||h||^2)  (x-space)
    print(f"{l+1:4d}   {c:+.4f}    {tot:.4f}     {vb:+.4f}                {share:.3f}                |  A: {np.linalg.norm(mc-tiltA)/np.linalg.norm(mc):.2e}   truth: {np.linalg.norm(mc-m)/np.linalg.norm(mc):.2e}")
d = np.diff(np.concatenate([np.zeros((len(logg), 1)), logg], 1), axis=1)
Cd = np.cov(d.T); sd = np.sqrt(np.diag(Cd)); Rd = Cd/np.outer(sd, sd)
print("log-gain increments: std per layer", np.round(sd, 4))
print("lag-1 autocorrelation", np.round([Rd[i, i+1] for i in range(L-1)], 3))
print("lag-2 autocorrelation", np.round([Rd[i, i+2] for i in range(L-2)], 3))
print("Var(log g_L) = %.4f ; sum of increment variances = %.4f" % (logg[:, -1].var(), np.sum(sd**2)))
