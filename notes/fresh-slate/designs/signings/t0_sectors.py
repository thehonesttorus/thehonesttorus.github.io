"""T0b: per-sector check of the Mehler generators against Monte Carlo joint cumulants of a_1 (exact copula)."""
import sys, numpy as np, copula
n = int(sys.argv[1]); N = int(float(sys.argv[2]))
rng = np.random.default_rng(1)
W = rng.standard_normal((2, n, n)) * np.sqrt(2 / n)
C1 = W[0].T @ W[0]; s = np.sqrt(np.diag(C1)); R = C1 / np.outer(s, s)
mu, st, dg = copula.layer_step(np.zeros(n), s, np.tile([1., 0, 0], (n, 1)), R, W[1])
w = W[1]
# exact means of a_1: s/sqrt(2pi)
A = []
E2 = np.zeros((n, n)); E22 = np.zeros((n, n)); E21 = np.zeros((n, n)); E31 = np.zeros((n, n))
E211 = np.zeros((n, n, n)); E111 = np.zeros((n, n, n)); M = np.zeros(4 * n).reshape(4, n)
done = 0; m0 = s / np.sqrt(2 * np.pi)
while done < N:
    x = rng.standard_normal((1 << 17, n)); a = np.maximum(x @ W[0], 0) - m0
    E2 += a.T @ a; E22 += (a**2).T @ (a**2); E21 += (a**2).T @ a; E31 += (a**3).T @ a
    E211 += np.einsum('ia,ib,ic->abc', a**2, a, a); E111 += np.einsum('ia,ib,ic->abc', a, a, a)
    done += 1 << 17
E2/=done; E22/=done; E21/=done; E31/=done; E211/=done; E111/=done
v = np.diag(E2)
K22 = E22 - np.outer(v, v) - 2 * E2**2
K21 = E21  # centred: kappa(a,a,b) = E[A^2 B]
K31 = E31 - 3 * v[:, None] * E2
K211 = E211 - v[:, None, None] * E2[None] - 2 * E2[:, :, None] * E2[:, None, :]
off = ~np.eye(n, dtype=bool)
def sec22(K): return 3 * np.einsum('aj,bj,ab->j', w**2, w**2, K * off)
def sec21(K): return 3 * np.einsum('aj,bj,ab->j', w**2, w, K * off)
def sec31(K): return 4 * np.einsum('aj,bj,ab->j', w**3, w, K * off)
dist = np.ones((n, n, n), bool)
for i in range(n): dist[i, i, :] = dist[i, :, i] = dist[:, i, i] = False
def sec211(K): return 6 * np.einsum('aj,bj,cj,abc->j', w**2, w, w, K * dist)
def sec111(K): return np.einsum('aj,bj,cj,abc->j', w, w, w, K * dist)
mc = {'3:{2,1}': sec21(K21), '3:{1,1,1}path': sec111(E111), '4:{2,2}': sec22(K22), '4:{3,1}': sec31(K31), '4:{2,1,1}': sec211(K211)}
for k, v_ in mc.items():
    e = dg['terms'][k]
    print(f"{k:16s} est rms {np.sqrt(np.mean(e**2)):.4f}  MC rms {np.sqrt(np.mean(v_**2)):.4f}  diff rms {np.sqrt(np.mean((e-v_)**2)):.4f}")
