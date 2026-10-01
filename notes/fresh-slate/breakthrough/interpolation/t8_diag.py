"""T8: TC diagnostics on bench MLP 0 at 1024: (i) self-consistent t vs MC t per layer; (ii) TC with oracle t (MC);
(iii) per-layer MSE of closure vs TC."""
import numpy as np, glob
from common import bench
import tc
S = bench.load_set("w1024_d16"); W = bench.weights(S, 0); T = S["means"][0]; nz = S["noise"][0]
p, ts = tc.predict(W, ret_state=True)
tmc = np.load(sorted(glob.glob("results/t5_w1024_d16_0_N*.npz"))[-1])["t"]
for l in range(16):
    c = np.corrcoef(ts[l], tmc[l])[0, 1]; rel = np.sqrt(((ts[l]-tmc[l])**2).mean()/(tmc[l]**2).mean())
    print(f"layer {l+1}: corr(t_tc, t_mc) {c:.4f}  rel rms err {rel:.3f}")
p0 = tc.predict(W, inject=False)
print("per-layer MSE closure:", " ".join(f"{x:.1e}" for x in ((p0 - T) ** 2).mean(1)))
print("per-layer MSE TC     :", " ".join(f"{x:.1e}" for x in ((p - T) ** 2).mean(1)))
# oracle t: monkeypatch by running predict with t replaced each layer
import types
def predict_oracle_t(W):
    Wd = W.astype(np.float64); n = W.shape[1]; sig2 = 2.0 / n
    from common import phi, Phi
    from t2_closures import hk
    out = []; mu = None; C = None
    for l in range(16):
        Wl = Wd[l]
        if l == 0: m = np.zeros(n); Sz = Wl.T @ Wl; y = np.zeros(n)
        else: m = mu @ Wl; Sz = Wl.T @ C @ Wl; y = Wl.T @ tmc[l - 1]
        v = np.diag(Sz); s = np.sqrt(v); a = m / s; P = Phi(a); pp = phi(a)
        k3 = 3 * sig2 * y
        mu_n = m * P + s * pp - k3 * a * pp / (6 * v); sec = (m*m+v)*P + m*s*pp + k3 * pp / (3 * s)
        R = Sz / np.outer(s, s); H = hk(a, 2)
        Cn = np.outer(s*H[0], s*H[0]) * R + 0.5*np.outer(s*H[1], s*H[1]) * R * R
        u1 = pp / s; Cn += 0.5 * sig2 * (np.outer(u1, P * y) + np.outer(P * y, u1))
        np.fill_diagonal(Cn, sec - mu_n ** 2); mu, C = mu_n, Cn; out.append(mu)
    return np.stack(out)
po = predict_oracle_t(W)
print("oracle-t TC final raw:", ((po[-1]-T[-1])**2).mean()-nz, " TC:", ((p[-1]-T[-1])**2).mean()-nz)
print("per-layer MSE oracle-t:", " ".join(f"{x:.1e}" for x in ((po - T) ** 2).mean(1)))
