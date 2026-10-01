"""E2: where does the Gaussian-closure error live in the tropical-curvature split (I2)?

E[a_lj] = beta_lj + T_lj,  beta_lj = E[delta(z_lj) |grad_x z_lj|^2]  (own-wall birth),
                           T_lj    = E[g_lj (Delta a_{l-1} W_l)_j]     (inherited curvature through the gate).
Measured by MC with exact per-sample input gradients; delta by a Gaussian kernel of width 0.03 s_lj (bias O(1e-3) rel).
Compared with their factorised (zeroth-order TCT) forms using the TRUE per-neuron mu, s:
beta_G = s phi(mu/s),  T_G = mu Phi(mu/s).  m_G = beta_G + T_G is the one-layer Gaussian readout.
Usage: python e2_wall_split.py seed n L N
"""
import numpy as np, sys
from scipy.stats import norm
seed, n, L, N = [int(float(a)) for a in (sys.argv[1:] + ['0', '64', '8', '200000'][len(sys.argv) - 1:])][:4]
rng = np.random.default_rng(seed)
W = [rng.standard_normal((n, n)) * np.sqrt(2 / n) for _ in range(L)]
# pass 1: per-neuron mean/std of pre-activations
S1 = np.zeros((L, n)); S2 = np.zeros((L, n)); M1 = np.zeros((L, n))
CH = 4000
for c in range(0, N, CH):
    a = rng.standard_normal((CH, n))
    for l in range(L):
        z = a @ W[l]; S1[l] += z.sum(0); S2[l] += (z * z).sum(0); a = np.maximum(z, 0); M1[l] += a.sum(0)
mu = S1 / N; s = np.sqrt(S2 / N - mu ** 2); m = M1 / N
bw = 0.03 * s
# pass 2: beta by kernel, with exact input gradients
B = np.zeros((L, n)); G2 = np.zeros((L, n)); Nb = N // 2
for c in range(0, Nb, 1000):
    X = rng.standard_normal((1000, n)); a = X; J = np.broadcast_to(np.eye(n), (1000, n, n)).copy()
    for l in range(L):
        z = a @ W[l]; G = J @ W[l]                      # G[b, i, j] = d z_lj / d x_i
        g2 = (G * G).sum(1)
        K = np.exp(-0.5 * (z / bw[l]) ** 2) / (np.sqrt(2 * np.pi) * bw[l])
        B[l] += (K * g2).sum(0); G2[l] += g2.sum(0)
        gate = z > 0; a = np.where(gate, z, 0); J = G * gate[:, None, :]
beta = B / Nb; G2 = G2 / Nb
Tt = m - beta
h = mu / s; betaG = s * norm.pdf(h); TG = mu * norm.cdf(h); mG = betaG + TG
rms = lambda v: np.sqrt((v ** 2).mean())
print(f"n={n} L={L} N={N}  (rms over neurons; MC noise of m ~ {np.sqrt(0.5/N):.1e})")
print("layer  rms(m)   err_m=m-mG  err_beta  err_T    corr(eb,eT)  E|grad z|^2/s^2  E|grad z|^2/(mu^2+s^2)  err(beta vs pG(0)E|grad|^2)")
for l in range(L):
    eb = beta[l] - betaG[l]; et = Tt[l] - TG[l]; em = m[l] - mG[l]
    print(f"{l+1:4d}  {rms(m[l]):.3f}   {rms(em):.2e}   {rms(eb):.2e}  {rms(et):.2e}  {np.corrcoef(eb, et)[0,1]:+.2f}      {np.median(G2[l]/s[l]**2):7.2f}          {np.median(G2[l]/(mu[l]**2+s[l]**2)):5.2f}              {rms(beta[l]-norm.pdf(h[l])/s[l]*G2[l]):.2e}")
