"""E3: temperature as a deformation parameter (softplus_T network, ReLU at T = 0).

Question: does the Gaussian closure become accurate fast enough with T that a Richardson / polynomial-in-T^2
extrapolation of closure values G(T), T >= T0, to T = 0 beats the closure at T = 0?
F(T): Monte Carlo truth of the softplus_T network (common random numbers across T).
G(T): Gaussian closure at T (per-neuron Hermite coefficients by Gauss-Hermite quadrature, covariance by Mehler
      series to order K; at T = 0 it reproduces TCT-0).
Usage: python e3_temperature.py set mlp_index N"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../bench"))
import numpy as np, bench
from math import factorial

def act(z, T):
    return np.maximum(z, 0) if T == 0 else T * np.logaddexp(0, z / T)

xg, wg = np.polynomial.hermite_e.hermegauss(160); wg = wg / wg.sum()
K = 40
He = np.zeros((K + 1, len(xg))); He[0] = 1; He[1] = xg
for k in range(1, K): He[k + 1] = xg * He[k] - k * He[k - 1]
kfact = np.array([float(factorial(k)) for k in range(K + 1)])

def closure(Ws, T):
    n = Ws[0].shape[0]; mu = np.zeros(n); C = Ws[0].T @ Ws[0]; out = []
    for l, W in enumerate(Ws):
        if l > 0:
            mu = m @ W; C = W.T @ Cov @ W
        s = np.sqrt(np.diag(C))
        F = act(mu[:, None] + s[:, None] * xg[None, :], T)          # (n, q)
        c = (F * wg) @ He.T                                          # c_k = E[f He_k]
        m = c[:, 0]; out.append(m)
        rho = C / np.outer(s, s)
        Cov = np.zeros((n, n)); rk = np.ones((n, n))
        for k in range(1, K + 1):
            rk = rk * rho; Cov += np.outer(c[:, k], c[:, k]) * rk / kfact[k]
        np.fill_diagonal(Cov, (F * F) @ wg - m * m)
    return np.array(out)

if __name__ == "__main__":
    setname, i, N = sys.argv[1], int(sys.argv[2]), int(float(sys.argv[3]))
    S = bench.load_set(setname); Ws = [w.astype(np.float64) for w in bench.weights(S, i)]
    n = Ws[0].shape[0]; L = len(Ws)
    Ts = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8]
    rng = np.random.default_rng(1000 + i); acc = np.zeros((len(Ts), n)); B = 100000
    for b in range(N // B):
        X = rng.standard_normal((B, n))
        for ti, T in enumerate(Ts):
            a = X
            for W in Ws: a = act(a @ W, T)
            acc[ti] += a.sum(0)
    Fm = acc / (B * (N // B))
    Gm = np.array([closure(Ws, T)[-1] for T in Ts])
    truth0 = S["means"][i][-1]
    print(f"{setname} mlp{i} n={n}: MC N={N} (noise MSE ~{0.05/N:.0e})")
    if truth0 is not None: print("  check MC(T=0) vs bench truth: mse %.1e" % ((Fm[0] - truth0) ** 2).mean())
    print("  T      mse(G-F)   rms(F(T)-F(0))")
    for ti, T in enumerate(Ts):
        print(f"  {T:4.2f}   {((Gm[ti]-Fm[ti])**2).mean():.2e}   {np.sqrt(((Fm[ti]-Fm[0])**2).mean()):.3e}")
    # extrapolate closure values from T >= T0 to T = 0 by polynomial in T^2 (and in T)
    for T0 in (0.05, 0.1, 0.2):
        idx = [k for k, T in enumerate(Ts) if T >= T0]
        for deg in (1, 2, 3):
            if len(idx) <= deg: continue
            t2 = np.array([Ts[k] for k in idx]) ** 2
            V = np.vander(t2, deg + 1)
            coef = np.linalg.lstsq(V, Gm[idx], rcond=None)[0]
            ext = coef[-1]
            print(f"  extrapolate G from T>={T0} deg {deg} in T^2: mse vs F(0) {((ext-Fm[0])**2).mean():.2e}   (closure at T=0: {((Gm[0]-Fm[0])**2).mean():.2e})")
