# Calibrate the residue with a SECOND-moment statistic: the variance of the normalized squared norm
#   g(x) = (1/n) sum_k z_k(x)^2 / q_k   (final-layer pre-activations, q_k = E z_k^2 from the closure),
# minus its Gaussian-closure value. gamma_hat = Var_hat(g) - Var_Gauss(g);  m <- (1 - gamma_hat/8) m_closure.
import numpy as np, sys
from closure import closure
from residue import closure_with_residue
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ts = [int(t) for t in sys.argv[4].split(",")]; reps = int(sys.argv[5])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); m = tr["m"][-1]
oc = closure(Ws, keep=True); d = oc[-1]
mc = d["m"]; mu = d["mu"]; S = d["R"]*np.outer(d["sig"], d["sig"]); q = mu**2 + np.diag(S)
VG = np.sum((2*S*S + 4*np.outer(mu, mu)*S)/np.outer(q, q))/n**2
c_or = mc @ (mc-m)/(mc @ mc)
w95 = closure_with_residue(Ws, tau=0.95)[-1]["m"]
print(f"n={n} L={L} s={s}: MC@B {tr['v'][-1].mean()/65536:.2e} | closure {np.mean((mc-m)**2):.2e} | weights-only residue(.95) {np.mean((w95-m)**2):.2e} | oracle scale {np.mean((mc*(1-c_or)-m)**2):.2e} (c={c_or:.4f}, 8c={8*c_or:.4f})", flush=True)
Wf = [W.T.astype(np.float32) for W in Ws]; rng = np.random.default_rng(100 + s)
for T in Ts:
    mses, gs = [], []
    for r in range(reps):
        h = rng.standard_normal((T, n)).astype(np.float32)
        for W in Wf[:-1]: h = np.maximum(h @ W, 0)
        z = (h @ Wf[-1]).astype(np.float64)
        g = ((z*z)/q).mean(1); gam = g.var(ddof=1) - VG
        gs.append(gam); mses.append(np.mean((mc*(1-gam/8)-m)**2))
    print(f"   T={T:5d} ({T/65374*100:.2f}% budget): gamma_hat {np.mean(gs):.4f} +- {np.std(gs):.4f} -> MSE {np.mean(mses):.2e} (worst {np.max(mses):.2e})", flush=True)
