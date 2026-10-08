# Input-chaos view: the second-chaos (Hessian) content born at layer b is sum_e rho_{b,e} (J_b)_e (J_b)_e^T in INPUT space,
# with J_b the mean (closure-gated) Jacobian input -> pre-activations of layer b and rho = phi(a)/sigma the kink density.
# Is the content of OLD births (b <= l-5) low-rank in input space?
import numpy as np, sys
sys.path.insert(0, "../num12")
from closure import closure
from scipy.special import ndtr
Ws = list(np.load("../official/W_off0.npy").astype(np.float64)); oc = closure(Ws)
J = Ws[0].copy(); Gs = []
for b in range(16):
    if b > 0: J = Ws[b] @ (ndtr(oc[b-1]["mu"]/oc[b-1]["sig"])[:, None]*J)
    rho = np.exp(-0.5*(oc[b]["mu"]/oc[b]["sig"])**2)/np.sqrt(2*np.pi)/oc[b]["sig"]
    Gs.append(J.T @ (rho[:, None]*J))
for l in [9, 12, 15]:
    G = sum(Gs[b] for b in range(l-4)); ev = np.linalg.eigvalsh(G)[::-1]; ev = np.maximum(ev, 0)
    cs = np.cumsum(ev)/ev.sum()
    print(f"layer {l+1}: old births b<= {l-4}: input-space energy in top r: " + " ".join(f"r={r}:{cs[r-1]:.3f}" for r in [1, 8, 32, 128, 384]) + f" | participation ratio {ev.sum()**2/(ev**2).sum():.1f}")
    Gy = Gs[l-5]; ev = np.maximum(np.linalg.eigvalsh(Gy)[::-1], 0); cs = np.cumsum(ev)/ev.sum()
    print(f"          single birth b={l-4}: " + " ".join(f"r={r}:{cs[r-1]:.3f}" for r in [1, 8, 32, 128, 384]) + f" | PR {ev.sum()**2/(ev**2).sum():.1f}")
