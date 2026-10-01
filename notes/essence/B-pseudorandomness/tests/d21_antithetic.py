"""Test T4 (B-pseudorandomness): does a Z2-structured sample design (antithetic pairs {x, -x}),
which fools every odd function exactly, cut the sampling noise of the (2,1) slice D21(z_l) at n = 1024?

For each pre-activation layer z_l (l = 1..16, z_1 = x W_1) we estimate D21_ab = kappa3(z_a, z_a, z_b) from
N forward passes, in two independent halves, with
  plain : N i.i.d. Gaussian inputs
  anti  : N/2 i.i.d. Gaussian inputs and their mirrors (same number of forward passes)
noise^2 of the full estimate = ||D1 - D2||_F^2 / 4 ; signal^2 = <D1, D2> (unbiased for ||D21||^2).
Also reports eta_l = fraction of the centred variance of z_l in its even part (z(x)+z(-x))/2.
Usage: python d21_antithetic.py MLP N
"""
import sys, json, time
import numpy as np
sys.path.insert(0, "/home/user/thehonesttorus.github.io/notes/fresh-slate/bench")
import bench

LAYERS = [1, 3, 5, 8, 10, 12, 15, 16]

def forward_collect(W, X):
    """returns dict l -> z_l (float32) for l in LAYERS"""
    out = {}
    a = X
    for l in range(1, W.shape[0] + 1):
        z = a @ W[l - 1]
        if l in LAYERS:
            out[l] = z
        a = np.maximum(z, 0)
    return out

class Acc:
    def __init__(self, n, m0):
        self.m0 = m0.astype(np.float64)
        self.N = 0
        self.s1 = np.zeros(n); self.s2 = np.zeros((n, n)); self.s21 = np.zeros((n, n)); self.s3 = np.zeros(n)
    def add(self, z):
        c = z.astype(np.float64) - self.m0
        self.N += c.shape[0]
        self.s1 += c.sum(0)
        self.s2 += c.T @ c
        c2 = c * c
        self.s21 += c2.T @ c
        self.s3 += (c2 * c).sum(0)
    def d21(self):
        N = self.N
        m = self.s1 / N; E2 = self.s2 / N; E21 = self.s21 / N
        d = np.diag(E2)
        # kappa3(a,a,b) = E[c_a^2 c_b] - 2 m_a E[c_a c_b] - m_b E[c_a^2] + 2 m_a^2 m_b
        return E21 - 2 * m[:, None] * E2 - d[:, None] * m[None, :] + 2 * (m * m)[:, None] * m[None, :]

def run(mlp, N, chunk=2048, seed=12345):
    S = bench.load_set("w1024_d16")
    W = bench.weights(S, mlp)
    n = W.shape[1]
    rng = np.random.default_rng(seed + mlp)
    # pilot means (independent sample) for numerically safe centring
    pil = forward_collect(W, rng.standard_normal((4096, n), dtype=np.float32))
    m0 = {l: pil[l].mean(0) for l in LAYERS}
    res = {}
    acc = {mode: {h: {l: Acc(n, m0[l]) for l in LAYERS} for h in (0, 1)} for mode in ("plain", "anti")}
    evar = {l: 0.0 for l in LAYERS}; tvar = {l: 0.0 for l in LAYERS}
    t0 = time.time()
    for h in (0, 1):
        for start in range(0, N // 2, chunk):
            X = rng.standard_normal((chunk, n), dtype=np.float32)
            zs = forward_collect(W, X)
            for l in LAYERS: acc["plain"][h][l].add(zs[l])
            Xh = X[: chunk // 2]
            zp = forward_collect(W, Xh); zm = forward_collect(W, -Xh)
            for l in LAYERS:
                acc["anti"][h][l].add(np.concatenate([zp[l], zm[l]], 0))
                c_p = zp[l] - m0[l]; c_m = zm[l] - m0[l]
                ev = 0.5 * (c_p + c_m); od = 0.5 * (c_p - c_m)
                evar[l] += float(((ev - 0) ** 2).sum()); tvar[l] += float((ev ** 2).sum() + (od ** 2).sum())
        print(f"half {h} done {time.time()-t0:.0f}s", flush=True)
    for l in LAYERS:
        row = {}
        for mode in ("plain", "anti"):
            D1 = acc[mode][0][l].d21(); D2 = acc[mode][1][l].d21()
            off = ~np.eye(n, dtype=bool)
            noise2 = ((D1 - D2) ** 2).sum() / 4; sig2 = (D1 * D2).sum()
            noise2o = ((D1 - D2)[off] ** 2).sum() / 4; sig2o = (D1[off] * D2[off]).sum()
            row[mode] = dict(rel_noise=float(np.sqrt(noise2 / sig2)) if sig2 > 0 else None,
                             rel_noise_off=float(np.sqrt(noise2o / sig2o)) if sig2o > 0 else None)
        pn, an = row["plain"]["rel_noise_off"], row["anti"]["rel_noise_off"]
        row["var_gain_off"] = (pn / an) ** 2 if (pn and an) else None  # layer 1 is exactly Gaussian: no signal
        row["eta_even"] = evar[l] / tvar[l]
        res[l] = row
        print(l, json.dumps(row), flush=True)
    return res

if __name__ == "__main__":
    mlp = int(sys.argv[1]); N = int(sys.argv[2])
    r = run(mlp, N)
    json.dump(r, open(f"/home/user/thehonesttorus.github.io/notes/essence/B-pseudorandomness/tests/d21_antithetic_mlp{mlp}_N{N}.json", "w"), indent=1)
