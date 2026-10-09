# Localization diagnostics, part 2 (the chain dump's 'mu' is the post-activation mean; the pre-activation mean is
# W_l out_(l-1), exactly, since the mean is linear):
#   python ana_loc2.py DIR NET [NET ...]
# E1 share of the squared mean error carried by neurons binned by the chain's gate margin |alpha| = |mu_pre| / sigma
# E5 injected error per layer: e_l - diag(Phi(alpha_l)) W_l e_(l-1), as a share of |e_l|^2
# E7 collective-mode content of the error: the chain's pre-activation mean error d_l = W_l e_(l-1) projected on the top
#    k eigenvectors of the Monte Carlo pre-activation covariance Cov(h_l) (k = 16, 64, 128), against the k/n of a
#    random vector; and the same for the variance-weighted error d_l / sigma_l
import sys, numpy as np
from scipy.special import ndtr
D = sys.argv[1]; nets = [int(x) for x in sys.argv[2:]]
bins = [0, 0.5, 1, 1.5, 2, 3, 99]; ks = (16, 64, 128)
R = {"e1": [], "e5": [], "e7": [], "e7w": []}
for net in nets:
    cd = np.load(f"{D}/chaindump_{net}.npz"); mc = np.load(f"{D}/mccache_{net}.npz")
    mt = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
    W = np.load(f"{D}/W_off{net}.npy").astype(np.float64)
    S = float(mc["S"]); L, n = mt.shape
    out = cd["out"]; err = out - mt; var = cd["var"]; m1 = mc["s1"] / S
    e1 = np.zeros((L, len(bins) - 1)); e5 = np.zeros(L); e7 = np.zeros((L, len(ks))); e7w = np.zeros((L, len(ks)))
    for l in range(L):
        mu_pre = np.zeros(n) if l == 0 else W[l] @ out[l - 1]
        sig = np.sqrt(np.maximum(var[l], 1e-12)); al = mu_pre / sig
        tot = np.sum(err[l] ** 2)
        for b in range(len(bins) - 1):
            sel = (np.abs(al) >= bins[b]) & (np.abs(al) < bins[b + 1])
            e1[l, b] = np.sum(err[l][sel] ** 2) / tot
        if l > 0:
            inj = err[l] - ndtr(al) * (W[l] @ err[l - 1])
            e5[l] = np.sum(inj ** 2) / tot
            d = W[l] @ err[l - 1]
            Hc = mc["Hc"][l].astype(np.float64) - np.outer(m1[l], m1[l])
            ev, V = np.linalg.eigh(Hc); V = V[:, ::-1]
            pr = (V.T @ d) ** 2; prw = (V.T @ (d / sig)) ** 2
            e7[l] = [pr[:k].sum() / pr.sum() for k in ks]; e7w[l] = [prw[:k].sum() / prw.sum() for k in ks]
    for k, v in (("e1", e1), ("e5", e5), ("e7", e7), ("e7w", e7w)):
        R[k].append(v)
    print(f"net {net} done", flush=True)
A = {k: np.mean(v, 0) for k, v in R.items()}
print("layer | sq. error share by |alpha| " + " ".join(f"[{bins[b]},{bins[b+1]})" for b in range(len(bins) - 1))
      + " | injected/total | pre-act. mean error in top-16/64/128 Cov modes (random: "
      + "/".join(f"{k / 1024:.3f}" for k in ks) + ") | same, sigma-weighted")
for l in range(len(A["e5"])):
    print(f"{l:2d} | " + " ".join(f"{x:6.3f}" for x in A["e1"][l]) + f" | {A['e5'][l]:.3f} | "
          + " ".join(f"{x:.3f}" for x in A["e7"][l]) + " | " + " ".join(f"{x:.3f}" for x in A["e7w"][l]))
