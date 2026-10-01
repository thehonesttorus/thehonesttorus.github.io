"""Mean-field recursion for the two fourth-order objects the means need (REPORT N7): the diagonal kappa4(z_l) and the
column means c22(z_l)_k = (1/n) sum_i kappa4(z_i, z_i, z_k, z_k). Memoryless in the all-distinct part, which the
fresh-weight lemma makes incoherent in every index sum; coherent sums close on themselves because (1/n) sum_i W_ri W_r'i
-> (2/n) delta_rr' in a coherent sum.
"""
import numpy as np
from math import comb, factorial
from scipy.special import ndtr

SQ2PI = np.sqrt(2 * np.pi)


def _He(k, x):
    h0, h1 = np.ones_like(x), x
    if k == 0:
        return h0
    for j in range(1, k):
        h0, h1 = h1, x * h1 - j * h0
    return h1


def gauss_derivs(mu, v, jmax=4, dmax=6):
    """E_G[(relu^j)^{(d)}(z)], z ~ N(mu, v), for j = 1..jmax, d = 0..dmax (distributional derivatives)."""
    s = np.sqrt(v); al = mu / s
    P = ndtr(al); p = np.exp(-0.5 * al * al) / SQ2PI
    I = [P, p, P - al * p, (al * al + 2) * p, -(al ** 3 + 3 * al) * p + 3 * P]
    MG = [P] + [sum(comb(k, j) * mu ** (k - j) * s ** j * I[j] for j in range(k + 1)) for k in range(1, jmax + 1)]
    out = {}
    for j in range(1, jmax + 1):
        for d in range(dmax + 1):
            if d < j:
                out[j, d] = factorial(j) / factorial(j - d) * MG[j - d]
            elif d == j:
                out[j, d] = factorial(j) * P
            else:
                k = d - j - 1
                out[j, d] = factorial(j) * _He(k, -al) * p / s ** (k + 1)
    return out, MG, P, p, s, al


def post_k4(mu, v, k3, k4):
    """Per-neuron kappa4 of relu(z) with Edgeworth corrections (k3/6 F''' + k4/24 F'''' + k3^2/72 F^(6))."""
    E, MG, P, p, s, al = gauss_derivs(mu, v)
    M = {j: E[j, 0] + k3 / 6 * E[j, 3] + k4 / 24 * E[j, 4] + k3 * k3 / 72 * E[j, 6] for j in range(1, 5)}
    m = M[1]
    mu2 = M[2] - m * m
    mu4 = M[4] - 4 * m * M[3] + 6 * m * m * M[2] - 3 * m ** 4
    return mu4 - 3 * mu2 * mu2


def step(mu, S, D21, k4z, c22z, Wn):
    """(mu, S, D21, k4, c22) of z_l -> (k4, c22) of z_{l+1} = a_l Wn."""
    n = len(mu)
    v = np.diag(S).copy(); s = np.sqrt(v); al = mu / s
    P = ndtr(al); p = np.exp(-0.5 * al * al) / SQ2PI; w2 = p / s
    m = mu * P + s * p
    D3 = np.diag(D21).copy() if D21 is not None else np.zeros(n)
    k4a = post_k4(mu, v, D3, k4z)
    eG1 = 2 * m * (1 - P); eG2 = 2 * (P - m * w2)
    C = S.copy(); np.fill_diagonal(C, 0)
    Dz = np.zeros((n, n)) if D21 is None else D21.copy()
    np.fill_diagonal(Dz, 0)
    K22z = 0.5 * (c22z[:, None] + c22z[None, :])
    T = (C * np.outer(eG1, eG1) + 0.5 * C * C * np.outer(eG2, eG2) - 2 * (C * np.outer(P, P)) ** 2
         + 0.25 * K22z * np.outer(eG2, eG2) + 0.5 * (Dz * np.outer(eG2, eG1) + Dz.T * np.outer(eG1, eG2)))
    np.fill_diagonal(T, 0)
    cK = T.mean(0)                                   # column means of K22(a)
    W2 = Wn.astype(np.float64) ** 2
    k4n = k4a @ (W2 * W2) + 3 * (cK @ W2) * W2.sum(0)
    c22n = (2.0 / n) * ((k4a + n * cK) @ W2)
    return k4n, c22n
