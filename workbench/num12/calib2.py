import numpy as np, sys
from closure import closure, phi
from scipy.special import ndtr
n, L, s, Tc = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz")
m, v = tr["m"][-1], tr["v"][-1]
o = closure(Ws)[-1]; a = o["mu"]/o["sig"]; sg = o["sig"]; mc = o["m"]
feats = {"m": [mc], "1,m": [np.ones(n), mc], "m,sig": [mc, sg], "m,sig,sig*phi": [mc, sg, sg*phi(a)],
         "m,sig,sig*phi,1": [mc, sg, sg*phi(a), np.ones(n)]}
rng = np.random.default_rng(3); base = np.mean((mc-m)**2)
print(f"n={n} L={L}: closure {base:.2e}, MC@B {v.mean()/65536:.2e}, Tc={Tc}")
for k, F in feats.items():
    X = np.array(F).T; out = []
    for rep in range(200):
        r = mc - (m + rng.standard_normal(n)*np.sqrt(v/Tc))
        b, *_ = np.linalg.lstsq(X, r, rcond=None); out.append(np.mean((mc - X@b - m)**2))
    b, *_ = np.linalg.lstsq(X, mc-m, rcond=None)
    print(f"   {k:18s}: {np.mean(out):.2e} (oracle {np.mean((mc-X@b-m)**2):.2e}) coef {np.round(b,4)}")
# scale ratio sum(m_true)/sum(m_closure)
print("   global scale ratio true/closure:", (m.sum()/mc.sum()))
