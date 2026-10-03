# Residue from a short Monte Carlo run, via the coherent fourth-order excess (all moments from the samples):
#   gamma_hat = avg_{k!=l} [E z_k^2 z_l^2 - E z_k^2 E z_l^2 - 2 (E z_k z_l)^2 + 2 mu_k^2 mu_l^2] / (E z_k^2 E z_l^2)
# computed with O(T n^2) work through the identity sum_{k,l} a_k a_l X_kl = a^T X a.
import numpy as np, sys
from closure import closure
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ts = [int(t) for t in sys.argv[4].split(",")]; reps = int(sys.argv[5])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); m = tr["m"][-1]
mc = closure(Ws)[-1]["m"]; c_or = mc @ (mc-m)/(mc @ mc)
print(f"n={n} L={L} s={s}: closure {np.mean((mc-m)**2):.2e}, oracle scale {np.mean((mc*(1-c_or)-m)**2):.2e}, 8c = {8*c_or:.4f}", flush=True)
Wf = [W.T.astype(np.float32) for W in Ws]; rng = np.random.default_rng(200 + s)
for T in Ts:
    gs, mses = [], []
    for r in range(reps):
        h = rng.standard_normal((T, n)).astype(np.float32)
        for W in Wf[:-1]: h = np.maximum(h @ W, 0)
        z = (h @ Wf[-1]).astype(np.float64)
        mu = z.mean(0); Ez2 = (z*z).mean(0); a = 1/Ez2
        y = (z*z)*a                                  # z_k^2 / E z_k^2
        sq = (y.sum(1)**2 - (y*y).sum(1)).mean()      # avg over samples of sum_{k!=l} y_k y_l
        M2 = (z.T @ z)/T                              # second-moment matrix
        Ga = a[:, None]*M2*a[None, :]
        gauss = 2*(np.sum(Ga*M2) - np.sum(np.diag(Ga)*np.diag(M2)))   # 2 sum_{k!=l} M2_kl^2 a_k a_l
        mm = (mu**2)*a
        corr = 2*(mm.sum()**2 - (mm*mm).sum())
        base = (n*n - n)                              # sum_{k!=l} E y_k E y_l = n(n-1)
        gam = (sq - base - gauss + corr)/(n*(n-1))
        gs.append(gam); mses.append(np.mean((mc*(1-gam/8)-m)**2))
    print(f"   T={T:5d} ({T/65374*100:.2f}% budget): gamma_hat {np.mean(gs):.4f} +- {np.std(gs):.4f} -> MSE {np.mean(mses):.2e} (worst {np.max(mses):.2e})", flush=True)
