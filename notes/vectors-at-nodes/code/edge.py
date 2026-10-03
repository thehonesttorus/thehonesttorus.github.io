# Edge/neuron formula for the mean of a bias-free ReLU net with Gaussian input.
#   E F_j = sum_{(l,i)} int_{z_li=0} delta_li |grad z_li| dgamma_{d-1}
#         = sum_{(l,i)} lim (1/2eps) E[delta_li |grad z_li|^2 ; |z_li|<eps]
# Window bias is linear in eps (kinked conditional density), so Richardson: 2 S(eps/2) - S(eps).
import numpy as np, sys

def net(d, n, L, rng):
    Ws = [rng.standard_normal((n, d)) * np.sqrt(2/d)]
    Ws += [rng.standard_normal((n, n)) * np.sqrt(2/n) for _ in range(L-1)]
    return Ws

def chunk_stats(Ws, X, js, epss):
    L = len(Ws); N, d = X.shape
    H = X; Jp = np.broadcast_to(np.eye(d), (N, d, d))
    Zs, G2 = [], []
    for W in Ws:
        Z = H @ W.T
        Jz = np.einsum('ab,nbc->nac', W, Jp)
        Zs.append(Z); G2.append((Jz**2).sum(2))
        H = np.maximum(Z, 0); Jp = (Z > 0)[:, :, None] * Jz
    out = {}
    for j in js:
        F = H[:, j]
        gradF = Jp[:, j, :]                                # grad F_j (N, d)
        dl = [None]*L
        g = np.zeros((N, Ws[-1].shape[0])); g[:, j] = 1.0
        dl[L-1] = g
        for l in range(L-1, 0, -1):
            dl[l-1] = (dl[l] * (Zs[l] > 0)) @ Ws[l]
        S = np.zeros((len(epss), L)); A = np.zeros((len(epss), L))
        for a, e in enumerate(epss):
            for l in range(L):
                w = (np.abs(Zs[l]) < e) * G2[l] / (2*e)
                S[a, l] = (dl[l] * w).sum(1).sum()
                A[a, l] = (np.abs(dl[l]) * w).sum(1).sum()
        # gradient spread E|gradF(X)-gradF(Y)| using a shifted copy as Y
        spread = np.linalg.norm(gradF - np.roll(gradF, 1, axis=0), axis=1).sum()
        zL = Zs[-1][:, j]
        out[j] = (F.sum(), (F**2).sum(), S, A, spread, zL.sum(), (zL**2).sum())
    return out

def run(d, n, L, N, seed, js, epss=(0.02, 0.01), chunk=100_000):
    rng = np.random.default_rng(seed)
    Ws = net(d, n, L, rng)
    acc = None
    for c in range(N // chunk):
        X = rng.standard_normal((chunk, d))
        o = chunk_stats(Ws, X, js, epss)
        if acc is None: acc = {j: [np.array(v, dtype=object) if False else v for v in o[j]] for j in js}
        else:
            for j in js:
                acc[j] = [a + b for a, b in zip(acc[j], o[j])]
    res = {}
    for j in js:
        s1, s2, S, A, sp, _, _ = acc[j]
        mc = s1/N; se = np.sqrt(s2/N - mc**2)/np.sqrt(N)
        S /= N; A /= N
        Sx = 2*S[1] - S[0]; Ax = 2*A[1] - A[0]           # Richardson eps -> 0
        res[j] = dict(mc=mc, se=se, S=Sx, A=Ax, Sraw=S, spread=sp/N)
    return res

if __name__ == '__main__':
    d, n, L, N = 8, 32, 4, 2_000_000
    res = run(d, n, L, N, seed=1, js=[0, 1, 2, 3, 4, 5])
    print(f"d={d} n={n} L={L} N={N}")
    for j, r in res.items():
        tot = r['S'].sum(); tv = r['A'].sum()
        print(f"j={j}: MC E F={r['mc']:.4f}+-{r['se']:.4f} | edge sum={tot:.4f} (raw eps .02/.01: {r['Sraw'][0].sum():.4f}/{r['Sraw'][1].sum():.4f})"
              f" | per layer {np.round(r['S'],4)} | ||Lap F||={tv:.4f} rho={r['mc']/tv:.3f} | spread={r['spread']:.4f} <= {np.sqrt(np.pi/2)*tv:.4f}")
