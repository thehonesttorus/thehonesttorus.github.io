"""Is the gate-pattern generating polynomial of a layer Lorentzian (degree-2 check at the all-ones point)?
f(t, t0) = E prod_i (g_i t_i + (1 - g_i) t0), g_i = 1[z_i > 0], z ~ N(mu, S) from the exact Gaussian closure.
Lorentzian => log-concave on the positive orthant => the Hessian at 1 has exactly one positive eigenvalue.
Pair orthant probabilities by the tetrachoric series (order 12). Usage: python t_lorentz.py MLP"""
import sys, json
import numpy as np
from scipy.special import ndtr
from numpy.polynomial.hermite_e import hermeval
sys.path.insert(0, "../../fresh-slate/bench"); sys.path.insert(0, "../../fresh-slate/breakthrough/region")
import bench, gclose as g
from math import factorial

mlp = int(sys.argv[1]) if len(sys.argv) > 1 else 0
S_ = bench.load_set("w1024_d16"); W = bench.weights(S_, mlp).astype(np.float64)
_, states = g.run(W, "exact")
L, n, _ = W.shape
res = []
for l in range(L):
    if l == 0:
        mu = np.zeros(n); S = W[0].T @ W[0]
    else:
        m, C = states[l - 1]; mu = m @ W[l]; S = W[l].T @ C @ W[l]
    s = np.sqrt(np.diag(S)); a = mu / s; R = S / np.outer(s, s); np.fill_diagonal(R, 0)
    p = ndtr(a); ph = np.exp(-a * a / 2) / np.sqrt(2 * np.pi)
    q = np.outer(p, p)
    for k in range(1, 13):
        c = np.zeros(k); c[-1] = 1
        hk = hermeval(a, c)  # He_{k-1}(a)
        q += R ** k / factorial(k) * np.outer(ph * hk, ph * hk)
    np.fill_diagonal(q, 0)
    cov = q - np.outer(p, p); np.fill_diagonal(cov, 0)
    H = np.zeros((n + 1, n + 1)); H[:n, :n] = q
    Egc = (p * (n - 1) - q.sum(1))
    H[:n, n] = H[n, :n] = Egc
    one = np.ones(n); Q1 = (np.outer(one, one) - np.outer(p, one) - np.outer(one, p) + q); np.fill_diagonal(Q1, 0)
    H[n, n] = Q1.sum()
    ev = np.linalg.eigvalsh(H)
    evc = np.linalg.eigvalsh(cov)
    r = dict(layer=l + 1, max_abs_rho=float(np.abs(R).max()), n_pos=int((ev > 1e-9 * ev[-1]).sum()),
             top_ev=[float(x) for x in ev[-4:][::-1]], cov_top=[float(x) for x in evc[-3:][::-1]],
             cov_bottom=float(evc[0]), mean_p=float(p.mean()), diag_shift=float(-(p * p).mean()))
    print(json.dumps(r), flush=True); res.append(r)
json.dump(res, open(f"results/lorentz_mlp{mlp}.json", "w"), indent=1)
