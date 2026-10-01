"""T12: where is TC's remaining error?  Oracle substitutions from the exact-centring MC (t9, MLP 0, N=131072):
per-layer z variance v (diag S), per-neuron kappa3(z), z mean m."""
import glob, numpy as np
from common import bench, phi, Phi, chi_mean_ratio
from t2_closures import hk
import tc
S_ = bench.load_set("w1024_d16"); W = bench.weights(S_, 0).astype(np.float64); T = S_["means"][0]; nz = S_["noise"][0]
d = dict(np.load(sorted(glob.glob("results/t9_w1024_d16_0_N*.npz"))[-1])); L, n = 16, 1024; sig2 = 2.0 / n
import os
if os.path.exists("results/t13_w1024_d16_0.npz"):
    h = np.load("results/t13_w1024_d16_0.npz"); Nh = float(h["N"]); dl = h["s1"] / Nh
    d["m2"] = h["s2"] / Nh - dl ** 2; d["m3"] = h["s3"] / Nh - 3 * dl * h["s2"] / Nh + 2 * dl ** 3
    d["mz"] = h["mz"].astype(np.float64) + dl; print("using t13, N =", Nh)
def run(var=False, k3o=False, mo=False, tmc=False, X4=None):
    out = []; mu = None; C = None; t = np.zeros(n)
    _, ts = tc.predict(W, ret_state=True)
    if tmc: ts = d['t']
    for l in range(L):
        Wl = W[l]
        if l == 0: m = np.zeros(n); S = Wl.T @ Wl; y = np.zeros(n)
        else: m = mu @ Wl; S = Wl.T @ C @ Wl; y = Wl.T @ ts[l - 1]
        if mo: m = d["mz"][l]
        if var:
            s0 = np.sqrt(np.diag(S)); s1 = np.sqrt(d["m2"][l]); S = S * np.outer(s1 / s0, s1 / s0)
        v = np.diag(S); s = np.sqrt(v); a = m / s; P = Phi(a); p = phi(a)
        k3 = d["m3"][l] if k3o else 3 * sig2 * y
        mu_n = m * P + s * p - k3 * a * p / (6 * v); sec = (m*m+v)*P + m*s*p + k3 * p / (3 * s)
        if X4 is not None and l > 0:
            k4 = 3 * sig2 ** 2 * X4[l - 1]; mu_n = mu_n + k4 * p * (a * a - 1) / (24 * v * s); sec = sec - k4 * a * p / (12 * v)
        R = S / np.outer(s, s); H = hk(a, 2)
        Cn = np.outer(s*H[0], s*H[0]) * R + 0.5*np.outer(s*H[1], s*H[1]) * R * R
        u1 = p / s; Cn += 0.5 * sig2 * (np.outer(u1, P * y) + np.outer(P * y, u1))
        np.fill_diagonal(Cn, sec - mu_n ** 2); mu, C = mu_n, Cn; out.append(mu)
    return np.stack(out) * chi_mean_ratio(n)
import sys
Xmc = np.load('/tmp/claude-0/s/t15_w1024_0.npz')['X']
for kw in [dict(var=True, k3o=True, mo=True, X4=Xmc), dict(var=True, k3o=True, X4=Xmc), dict(X4=Xmc)]:
    p = run(**kw); print({k: (v if k != 'X4' else 'X') for k, v in kw.items()}, f"final raw {((p[-1]-T[-1])**2).mean()-nz:.3e}", flush=True)
for kw in []:
    p = run(**kw); print(kw, f"final raw {((p[-1]-T[-1])**2).mean()-nz:.3e}", " per-layer:", " ".join(f"{x:.1e}" for x in ((p-T)**2).mean(1)[[3,7,11,15]]), flush=True)
