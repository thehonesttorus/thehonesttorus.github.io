# Scale-mixture (gain) theory for the per-neuron fourth cumulant of z16:  z = G y, y_i ~ N(mu_i, s_i^2), G^2 ~ Gamma(1/g, g)
import numpy as np, sys
from math import lgamma, exp
sys.path.insert(0, "../num12")
from gac import gac
f = np.load("final_stats_off0.npz"); k4mc, var_mc, k3mc = f["k4_mc"], f["var_mc"], f["k3_mc"]; g4c = f["g4_c"]; mu_true = f["mu_true"]
Ws = list(np.load("../official/W_off0.npy").astype(np.float64)); last = gac(Ws)[-1]; gam = last["gamma"]
def EGk(k, g):  # E[G^k] for G^2 ~ Gamma(shape 1/g, scale g)
    a = 1/g; return exp(lgamma(a + k/2) - lgamma(a) - (k/2)*np.log(a))
def mix_k4(mu, s2, g):
    E = [EGk(k, g) for k in range(5)]
    y1 = mu; y2 = mu**2 + s2; y3 = mu**3 + 3*mu*s2; y4 = mu**4 + 6*mu**2*s2 + 3*s2**2
    z1, z2, z3, z4 = E[1]*y1, E[2]*y2, E[3]*y3, E[4]*y4
    return z4 - 4*z3*z1 - 3*z2**2 + 12*z2*z1**2 - 6*z1**4, z2 - z1**2
sig4 = var_mc**2
def report(name, k4):
    sl = np.sum(k4*k4mc)/np.sum(k4mc*k4mc); rms = np.sqrt(np.mean(((k4 - k4mc)/sig4)**2))
    print(f"{name:52s} slope vs MC {sl:.3f}  mean diff/sig^4 {np.mean((k4-k4mc)/sig4):+.4f}  rms {rms:.4f}  (MC noise ~0.0135)")
report("chain fitted regen (g4row)", g4c)
for gname, g in [("GAC weights-only gamma", gam), ("gamma x 1.2", 1.2*gam)]:
    # conditional state: match the TRUE mean and the MC variance of z (marginal), solve for (mu_c, s2_c) given g
    g1, g2 = EGk(1, g), EGk(2, g)
    mu_c = mu_true/g1; s2_c = (var_mc + mu_true**2)/g2 - mu_c**2
    k4, v = mix_k4(mu_c, s2_c, g); report(f"scale mixture, {gname} = {g:.4f}", k4)
print(f"MC median k4/sig^4 {np.median(k4mc/sig4):.4f};  3*gamma = {3*gam:.4f}")
# regression of MC k4 on [chain g4, mixture k4] to see what each carries
k4mix, _ = mix_k4(mu_true/EGk(1, gam), (var_mc + mu_true**2)/EGk(2, gam) - (mu_true/EGk(1, gam))**2, gam)
X = np.stack([g4c, k4mix - 0], 1); c, *_ = np.linalg.lstsq(X/sig4[:, None], k4mc/sig4, rcond=None)
res = k4mc/sig4 - (X/sig4[:, None]) @ c
print(f"joint fit k4_MC ~ {c[0]:.3f} chain + {c[1]:.3f} mixture: residual rms {np.sqrt(np.mean(res**2)):.4f}")
