# Generator-evaluator suite at the competition shape (official network, width 1024, depth 16).
import numpy as np, time, sys
from scipy.special import ndtr
i = int(sys.argv[1]) if len(sys.argv) > 1 else 0
Wcol = np.load(f"../official/W_off{i}.npy").astype(np.float64)        # column convention: z = W h
mt = np.load(f"../official/truth_off{i}.npz")["m"][-1].astype(np.float64)
Wf = [np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol]       # row form: h @ Wf
n = 1024; rng = np.random.default_rng(11); N = 8192
def tail(h1):                                   # layers 2..L on a batch of layer-1 activations
    h = h1.astype(np.float32)
    for W in Wf[1:]: h = np.maximum(h @ W, 0)
    return h.astype(np.float64)
X = rng.standard_normal((N, n)).astype(np.float64)
Z1 = X @ Wcol[0].T                              # layer-1 pre-activations
F = tail(np.maximum(Z1, 0)); V = F.var(0).mean(); mu_mc = F.mean(0)
print(f"plain MC: avg per-neuron Var {V:.4f}  -> MSE at the 10% floor (6554 samples) {V/6554:.2e}")
# E1 radial: E F(X) = E[R] E F(Theta); Var of E[R] F(Theta) sample vs F(X)
r = np.linalg.norm(X, axis=1); ER = np.sqrt(2)*np.exp(np.math.lgamma((n+1)/2) - np.math.lgamma(n/2)) if hasattr(np, "math") else None
from math import lgamma, exp, sqrt
ER = sqrt(2)*exp(lgamma((n+1)/2) - lgamma(n/2)); Fr = F/r[:, None]*ER
print(f"E1 radial: Var ratio {Fr.var(0).mean()/V:.4f}")
# E2 optimal importance sampling bound for the vector residual R = F - mu (best constant analytic part)
Rn = np.linalg.norm(F - mt, axis=1)
print(f"E2 optimal IS (q* ~ phi*||R||): best possible variance factor {np.mean(Rn)**2/np.mean(Rn**2):.4f}")
# E3 corrected restriction hierarchy in the W1 right-singular frame, exact psi correction at layer 1
U, s, Vt = np.linalg.svd(Wcol[0]); Xc = X @ Vt.T          # coordinates in the right-singular frame
def psi(a, sd):
    out = np.maximum(a, 0).copy(); m = sd > 1e-9
    t = a[:, m]/sd[m]; out[:, m] = a[:, m]*ndtr(t) + sd[m]*np.exp(-0.5*t*t)/np.sqrt(2*np.pi); return out
levels = [0, 16, 64, 256, 1024]; Fl = []
for k in levels:
    Wk = Wcol[0] @ Vt[:k].T                        # rows of W1 restricted to the revealed frame
    a = Xc[:, :k] @ Wk.T if k > 0 else np.zeros((N, n))
    sd = np.sqrt(np.maximum((Wcol[0]**2).sum(1) - (Wk**2).sum(1), 0))
    Fl.append(tail(psi(a, sd)) if k < n else F)
F0 = Fl[0][0]
print(f"E3 level-0 evaluator (deterministic): MSE vs truth {np.mean((F0-mt)**2):.2e}")
tot = 0.0
for j in range(1, len(levels)):
    D = Fl[j] - Fl[j-1]; Vj = D.var(0).mean(); cj = 2.0 if j < len(levels)-1 else 1.0 + 15/16
    tot += np.sqrt(Vj*cj); print(f"   level {levels[j-1]:4d} -> {levels[j]:4d}: Var(D_j)/Var(F) {Vj/V:.4f}, mean |E D_j| {np.sqrt(np.mean(D.mean(0)**2)):.2e}")
print(f"E3 MLMC work / plain MC work (optimal allocation): {tot**2/V:.3f}   (<1 means the hierarchy wins)")
# E4 paired source decomposition around x0 = revealed part (k=256), v = complement: D = (F(x0+v)+F(x0-v))/2 - F(x0)
k = 256; x0 = Xc[:, :k] @ Vt[:k]; v = X - x0
Fp = tail(np.maximum((x0 + v) @ Wcol[0].T, 0)); Fm = tail(np.maximum((x0 - v) @ Wcol[0].T, 0)); F00 = tail(np.maximum(x0 @ Wcol[0].T, 0))
Dp = (Fp + Fm)/2 - F00
print(f"E4 paired correction (k={k}): Var(D)/Var(F) {Dp.var(0).mean()/V:.4f}; Var(F(x0))/Var(F) {F00.var(0).mean()/V:.4f}")
print(f"E5 frontier scale: the open-source chain's final MSE is ~2.1e-8; to match it by sampling at the 10% floor needs Var/{V/6554/2.1e-8:.0f}")
