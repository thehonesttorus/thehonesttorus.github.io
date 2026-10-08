# Operator theory of the kappa3 legs (ray-compiler note section 7): free-probability effective rank of gated Jacobian
# products, against their exact spectra, and the observability-weighted (balanced) spectra that decide the rank a tier needs.
#   python legspec.py NET MCPREFIX [--tails]
# A source born at the post-activation of layer b is carried to the pre-activation of layer l by
#   J(b->l) = W_l diag(Phi_(l-1)) ... W_(b+1) diag(Phi_b)          (first-order product gate, Phi = Phi(mu/sigma), truth)
# Free probability (W_k^T W_k -> 2 MP(1), rho = 2; gates rho_Phi = E Phi^4 / (E Phi^2)^2; rho(a [x] b) = rho(a) + rho(b) - 1):
#   PR(b->l) = n / (1 + a + sum_(k=b)^(l-1) (rho_Phi(k) - 1)),   a = l - b.
# Observability of row i at layer l (diagonal Gramian under the product gate), readout weights for the two channels the
# legs feed: D3 through the mean, r3 = c3^2 (c3 = E relu''' / 6), and D21 through the covariance, r21 = rho^2 E Phi^2
# (rho = phi(alpha)/sigma); both scaled by the suffix weight K'(l) of the layer:
#   o_l = s_l (r3_l + r21_l) + Phi_l^2 o (W_(l+1)^2)^T o_(l+1)
# Reports, per (b, age): exact PR vs the prediction, the energy fraction captured by the top r singular directions
# (r = 64 ... 512), and the same for the observability-weighted legs diag(o_l)^(1/2) J(b->l).
# --tails: instead, every age 1..12 at the four birth layers, the tail fraction eps(a, r) = 1 - captured energy at a fine
# rank grid, its local log-slope p = -dlog eps / dlog r, and the same tails for a FREE SURROGATE: the same gates, fresh He
# Gaussian weights (seed 7). If the legs are free products, the surrogate's spectrum is the network's: the rank a given
# tail needs is then a function of the gate statistics alone.
import sys, math, numpy as np
from math import erf
net, pre = int(sys.argv[1]), sys.argv[2]
W = np.load(f"../official/W_off{net}.npy").astype(np.float64)
L, n, _ = W.shape
F = np.load(f"{pre}_full.npz")
mu, var = F["mu"].astype(np.float64), F["var"].astype(np.float64)
sd = np.sqrt(var); al = mu / sd
Phi = 0.5 * (1 + np.vectorize(erf)(al / math.sqrt(2))); phi = np.exp(-al * al / 2) / math.sqrt(2 * math.pi)
rhoPhi = np.mean(Phi ** 4, axis=1) / np.mean(Phi ** 2, axis=1) ** 2
if "--tails" in sys.argv:
    ranks = (64, 96, 128, 160, 192, 224, 256, 288, 320, 352, 384, 448, 512)
    rng = np.random.default_rng(7)
    G = rng.standard_normal((L, n, n)) * math.sqrt(2.0 / n)          # free surrogate weights (same He law)
    print(f"net {net}: tails eps(a, r) = 1 - (energy of J(b->l) in its top r directions); 'free' = same gates, fresh He weights")
    print("  gain per transport 2 E Phi^2 by layer " + " ".join(f"{2 * float(np.mean(Phi[k] ** 2)):.3f}" for k in range(L)))
    acc = {}
    for b in (2, 5, 8, 11):
        J = np.eye(n); Jf = np.eye(n)
        for a in range(1, L - b):
            l = b + a
            J = W[l] @ (Phi[l - 1][:, None] * J); Jf = G[l] @ (Phi[l - 1][:, None] * Jf)
            if a > 12:
                break
            s = np.linalg.svd(J, compute_uv=False) ** 2; sf = np.linalg.svd(Jf, compute_uv=False) ** 2
            e = 1 - np.cumsum(s) / s.sum(); ef = 1 - np.cumsum(sf) / sf.sum()
            pr, prf = s.sum() ** 2 / np.sum(s * s), sf.sum() ** 2 / np.sum(sf * sf)
            pred = n / (1 + a + np.sum(rhoPhi[b:l] - 1))
            acc.setdefault(a, []).append((e, ef, pr, prf, pred, float(s.sum() / n)))
            print(f"   b={b:2d} a={a:2d}: PR {pr:6.1f} free-surrogate {prf:6.1f} predicted {pred:6.1f} | |J|_F^2/n {s.sum() / n:.3e} | "
                  "eps " + " ".join(f"{e[r - 1]:.1e}" for r in ranks) + " | free " + " ".join(f"{ef[r - 1]:.1e}" for r in ranks), flush=True)
    print("  mean over birth layers: age | PR exact / free / predicted | eps(r) for r = " + " ".join(str(r) for r in ranks))
    for a in sorted(acc):
        e = np.mean([x[0] for x in acc[a]], axis=0); ef = np.mean([x[1] for x in acc[a]], axis=0)
        print(f"   a={a:2d} | {np.mean([x[2] for x in acc[a]]):6.1f} / {np.mean([x[3] for x in acc[a]]):6.1f} / "
              f"{np.mean([x[4] for x in acc[a]]):6.1f} | " + " ".join(f"{e[r - 1]:.2e}" for r in ranks)
              + " | free " + " ".join(f"{ef[r - 1]:.2e}" for r in ranks))
    print("  local log-slope p = -dlog eps/dlog r between consecutive ranks (exact legs):")
    for a in sorted(acc):
        e = np.mean([x[0] for x in acc[a]], axis=0)
        ps = [-(math.log(max(e[r2 - 1], 1e-12)) - math.log(max(e[r1 - 1], 1e-12))) / math.log(r2 / r1) for r1, r2 in zip(ranks[:-1], ranks[1:])]
        print(f"   a={a:2d}: " + " ".join(f"{p:5.1f}" for p in ps))
    sys.exit()
# suffix weights K'(l) of a mean error at the post-activation of layer l (infinite-width kernel, as suffix_kernel.py)
c = [1 / math.pi]
for _ in range(L - 1):
    x = c[-1]; c.append((math.sqrt(1 - x * x) + (math.pi - math.acos(x)) * x) / math.pi)
Kp = []
for l in range(L):
    k = 1.0
    for j in range(l, L - 1):
        k *= 0.5 + math.asin(c[j]) / math.pi
    Kp.append(k)
c3 = (-al) * phi / (6 * var)                       # He1(-alpha) phi / (6 sigma^2)
r21 = (phi / sd) ** 2 * np.mean(Phi ** 2, axis=1)[:, None]
o = np.zeros((L, n))
o[L - 1] = Kp[L - 1] * c3[L - 1] ** 2
for l in range(L - 2, -1, -1):
    o[l] = Kp[l] * (c3[l] ** 2 + r21[l]) + Phi[l] ** 2 * ((W[l + 1] ** 2).T @ o[l + 1])
ranks = (64, 128, 192, 256, 320, 384, 512)
print(f"net {net}: gate spread rho_Phi by layer " + " ".join(f"{x:.2f}" for x in rhoPhi))
print("  observability: share of o_l carried by the readout term (rest = transmission) "
      + " ".join(f"{float(np.sum(Kp[l] * (c3[l] ** 2 + r21[l])) / np.sum(o[l])):.2f}" for l in range(L - 1)))
print("  observability concentration: participation n_eff = (sum o)^2 / sum o^2 by layer "
      + " ".join(f"{float(o[l].sum() ** 2 / np.sum(o[l] ** 2)):.0f}" for l in range(L)))
print("  (b, age): PR exact | predicted | energy captured by top r = " + " ".join(str(r) for r in ranks)
      + " || observable PR | observable energy captured by top r")
for b in (2, 5, 8, 11):
    J = np.eye(n)
    for a in range(1, L - b):
        l = b + a
        J = W[l] @ (Phi[l - 1][:, None] * J)
        if a not in (1, 2, 4, 6, 8, 10, 12) or l >= L:
            continue
        sv2 = np.linalg.svd(J, compute_uv=False) ** 2
        pr = sv2.sum() ** 2 / np.sum(sv2 ** 2)
        pred = n / (1 + a + np.sum(rhoPhi[b:l] - 1))
        cap = np.cumsum(sv2) / sv2.sum()
        Jo = np.sqrt(o[l])[:, None] * J
        so2 = np.linalg.svd(Jo, compute_uv=False) ** 2
        pro = so2.sum() ** 2 / np.sum(so2 ** 2); capo = np.cumsum(so2) / so2.sum()
        print(f"   ({b:2d},{a:2d}) l={l:2d}: {pr:6.1f} | {pred:6.1f} | " + " ".join(f"{cap[r - 1]:.3f}" for r in ranks)
              + f" || {pro:6.1f} | " + " ".join(f"{capo[r - 1]:.3f}" for r in ranks), flush=True)
