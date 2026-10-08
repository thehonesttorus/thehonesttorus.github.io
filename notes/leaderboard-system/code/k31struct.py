# What does the true (3,1) slice look like? (note XXXIX). K31[a, c] = kappa(z_a, z_a, z_a, z_c), zero diagonal, the
# chain's wk431 convention (the Monte Carlo table transposed). Explained shares of its noise-corrected energy by
# candidate shapes X, each with one free coefficient per row (an upper bound for "row-scaled X" closures) and with the
# coefficient the theory names. Candidates come from the gain mode (kappa_3 tangent 6 t mu.Sigma, kappa_4 tangent
# 12 t Sigma.Sigma: K31 = 12 t var_a C_ac, D21 = 2 t (2 mu_a C_ac + mu_c var_a)) and from the collective (B-independent)
# model, in which K31 and D21 share their column space.
#   python k31struct.py NET MCFULL MCH0 MCH1 [CHAIN_DUMP]
# Shares use the truth's own C, var, mu, D21 (structure of the truth); with CHAIN_DUMP also the chain's D21 and C_off as
# X (what the chain could read from its own state). Plain metric and the read metric c(1,3)_a Phi_c of the covariance
# program (c(1,3) = -alpha phi / sigma^2).
import sys, numpy as np
from math import erf

net, ff, f0, f1 = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
ch = np.load(sys.argv[5]) if len(sys.argv) > 5 else None
F, H0, H1 = np.load(ff), np.load(f0), np.load(f1)
n = F["mu"].shape[1]
Phi = np.vectorize(lambda a: 0.5 * (1 + erf(a / np.sqrt(2))))


def zd(x):
    x = np.array(x, dtype=np.float64); np.fill_diagonal(x, 0.0); return x


def flat(D):
    S = (D + D.T) / 2; A = (D - D.T) / 2; off = ~np.eye(n, dtype=bool)
    c = S[off].mean(); al = (S.sum(1) - (n - 1) * c) / (n - 2); be = A.sum(1) / n
    return zd(S - c - al[:, None] - al[None, :] + A - (be[:, None] - be[None, :]))


def shares(K, noise, X, w=None):
    # per-row free coefficient (upper bound), one global coefficient, energy explained over noise-corrected energy
    Kw, Xw = (K, X) if w is None else (K * w, X * w)
    E = np.sum(Kw ** 2) - noise
    num = np.sum(Kw * Xw, axis=1); den = np.maximum(np.sum(Xw ** 2, axis=1), 1e-300)
    rowfit = np.sum(num ** 2 / den) / E
    g = np.sum(num) / np.sum(den); glob = (2 * g * np.sum(num) - g * g * np.sum(den)) / E
    return rowfit, glob, g


print(f"net {net}: explained share of the noise-corrected K31 energy (per-row coefficient | one coefficient, its value)")
for l in range(1, 16):
    K = zd(F["K31"][l].T); dn = zd(H0["K31"][l].T) - zd(H1["K31"][l].T)
    var = F["var"][l].astype(np.float64); mu = F["mu"][l].astype(np.float64); sd = np.sqrt(var); al = mu / sd
    C = zd(F["cov"][l]); D = zd(F["D21"][l])
    w = ((-al * np.exp(-0.5 * al * al) / np.sqrt(2 * np.pi)) / var)[:, None] * Phi(al)[None, :]
    noise, noisew = np.sum(dn ** 2) / 4, np.sum((dn * w) ** 2) / 4
    cand = {"C": C, "var_a C": var[:, None] * C, "C var_c": C * var[None, :], "C2": zd(C @ C),
            "D21": D, "D21^T": D.T, "flat D21": flat(D), "gain D21 (3var/mu)": (3 * var / np.where(np.abs(mu) > 1e-12, mu, 1e-12))[:, None] * flat(D)}
    if ch is not None and f"D21_{l}" in ch:
        cand["chain flat D21"] = flat(zd(ch[f"D21_{l}"])); cand["chain C_off"] = zd(ch[f"C_off_{l}"])
    out = []
    for name, X in cand.items():
        r, gl, g = shares(K, noise, X)
        rw, glw, _ = shares(K, noisew, X, w)
        out.append(f"{name}: {r:.3f}|{gl:.3f} ({g:.3g}) w {rw:.3f}|{glw:.3f}")
    print(f"{l:2d} noise {noise / np.sum(K ** 2):.3f} | " + " ; ".join(out), flush=True)
