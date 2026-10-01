"""T4: in-chain disorder-conditional (self-averaging) correction.  Gaussian closure (Hermite order K) where, at every
layer l, the predicted mean is corrected by e_l(f_j) = basis(f_j) @ c_l, f_j = (alpha_j, s_j).  c_l are fitted
sequentially (layer by layer, with earlier corrections applied) on TRAINING networks (fresh seeds, MC truth), then the
chain is evaluated on the bench networks.  Usage: t4_chaincorr.py [deg] [K] [ridge]"""
import sys, glob, json
import numpy as np
from common import bench, phi, Phi, chi_mean_ratio
from t2_closures import hk

def basis(a, s, deg):
    cols = [s * a ** k for k in range(deg + 1)] + [s * phi(a) * a ** k for k in range(deg + 1)] + [a ** k for k in range(deg + 1)]
    return np.stack(cols, 1)

class Chain:
    def __init__(self, W, K):
        self.W = W.astype(np.float64); self.K = K; self.l = 0
        self.mu = None; self.C = None
    def pre(self):
        W = self.W[self.l]
        if self.l == 0: m = np.zeros(W.shape[0]); S = W.T @ W
        else: m = self.mu @ W; S = W.T @ self.C @ W
        v = np.diag(S); s = np.sqrt(v); a = m / s
        self.m, self.S, self.s, self.a = m, S, s, a
        return a, s
    def post(self, corr):
        m, S, s, a = self.m, self.S, self.s, self.a
        P = Phi(a); p = phi(a)
        mu = m * P + s * p + corr; sec = (m * m + s * s) * P + m * s * p
        R = S / np.outer(s, s); H = hk(a, self.K); C = np.zeros_like(S); Rk = np.ones_like(S); f = 1.0
        for k in range(1, self.K + 1):
            Rk = Rk * R; f *= k; C += np.outer(s * H[k - 1], s * H[k - 1]) * Rk / f
        np.fill_diagonal(C, np.maximum(sec - mu * mu, 1e-12))
        self.mu, self.C = mu, C; self.l += 1
        return mu

def fit(train, deg, K, ridge, L):
    chains = [Chain(W, K) for W, _ in train]; coefs = []
    for l in range(L):
        X = []; Y = []
        for ch, (_, T) in zip(chains, train):
            a, s = ch.pre(); X.append(basis(a, s, deg)); Y.append(T[l] - (ch.m * Phi(a) + s * phi(a)))
        X = np.concatenate(X); Y = np.concatenate(Y)
        c = np.linalg.solve(X.T @ X + ridge * np.eye(X.shape[1]), X.T @ Y) if l > 0 else np.zeros(X.shape[1])
        coefs.append(c)
        for ch in chains:
            ch.post(basis(ch.a, ch.s, deg) @ c)
    return coefs

def predict(W, coefs, deg, K):
    ch = Chain(W, K); out = []
    for c in coefs:
        a, s = ch.pre(); out.append(ch.post(basis(a, s, deg) @ c))
    return np.stack(out)

if __name__ == "__main__":
    deg = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    ridge = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6
    files = sorted(glob.glob("results/train/w1024_d16_s*.npz"))
    train = [(bench.weights_from_seed(int(f.split("_s")[-1][:-4]), 1024, 16), np.load(f)["means"]) for f in files]
    print("training networks:", len(train), flush=True)
    coefs = fit(train, deg, K, ridge, 16)
    S = bench.load_set("w1024_d16"); rows = []
    zero = [np.zeros_like(c) for c in coefs]
    for i in range(len(S["seeds"])):
        W = bench.weights(S, i); T = S["means"][i]; nz = S["noise"][i]
        p0 = predict(W, zero, deg, K); p1 = predict(W, coefs, deg, K)
        r0 = ((p0[-1] - T[-1]) ** 2).mean() - nz; r1 = ((p1[-1] - T[-1]) ** 2).mean() - nz
        r0c = ((chi_mean_ratio(1024) * p0[-1] - T[-1]) ** 2).mean() - nz
        al = ((p1 - T) ** 2).mean(1)
        rows.append((r0, r0c, r1)); print(i, f"closure {r0:.3e}  +chi {r0c:.3e}  in-chain corrected {r1:.3e}  bias {np.mean(p1[-1]-T[-1]):+.1e}", flush=True)
    rows = np.array(rows); print("mean:", " ".join(f"{x:.3e}" for x in rows.mean(0)))
    json.dump(dict(deg=deg, K=K, ridge=ridge, ntrain=len(train), rows=rows.tolist(), coefs=[c.tolist() for c in coefs]),
              open(f"results/t4_deg{deg}_K{K}_n{len(train)}.json", "w"), indent=1)
