# Note XLVII: the crossed-product (Wiener-chaos) state of the network, built with exact first-jet transport, and its
# third-cumulant slices against the Monte Carlo truth and against the production chain.
#   python chaos_state.py NET [JETS]      JETS = true (Monte Carlo mean/variance at every fold, default) | chain
# The input x ~ N(0, I) is the free (semicircular) system. Every pre-activation is a functional of x:
#   z^l = mu^l + L^l x + (1/2) sum_s V_is J2_s ((l_s . x)^2 - |l_s|^2) + ...
# - L^l = E[dz^l/dx] (first chaos), L^0 = W_0, L^(m+1) = W_(m+1) diag(Phi_m) L^m;
# - a fold at (m, a) is a rank-one second-chaos source: birth direction l_s = row a of L^m, jet J2_s = p_(z_a)(0) =
#   phi(alpha)/S (the wall density, renormalised at the law's own mean and variance), leg V_s = column a of W_(m+1),
#   then conjugated forward by the transports, V <- W diag(Phi) V (sources carried as histories, never compressed);
# Exact Wick contractions (second chaos, path terms) give, with U = L^l Lambda (U_is = Cov(z_i, l_s . x)):
#   kappa3(z_i, z_j, z_k) = sum_s J2_s (V_is U_js U_ks + V_js U_is U_ks + V_ks U_is U_js)
#   D21_ij = kappa(i, i, j) = sum_s J2_s (2 V_is U_is U_js + V_js U_is^2),  D3_i = 3 sum_s J2_s V_is U_is^2
# Compared on neurons with true alpha > -2.5 (active by active, off-diagonal), energies noise-free across the halves.
import os, sys, numpy as np
from scipy.special import ndtr
net = int(sys.argv[1]); jets = sys.argv[2] if len(sys.argv) > 2 else "true"
OFF = os.environ.get("OFFDIR", "../official")
W = np.load(f"{OFF}/W_off{net}.npy").astype(np.float64)
T = {h: np.load(f"mc2_off{net}_{h}.npz") for h in ("full", "h0", "h1")}
c2 = np.load(f"chaindump2_{net}.npz")
n = W.shape[1]
off = ~np.eye(n, dtype=bool)
f64 = lambda a: np.asarray(a, dtype=np.float64)
SQ2PI = np.sqrt(2 * np.pi)
EVAL = [int(x) for x in os.environ.get("EVAL", "2,3,4,6,8,10,12,13,14,15").split(",")]


def state_jets(m):
    if jets == "true" or m == 0 or m > 14:
        mu, v = f64(T["full"]["mu"][m]), f64(T["full"]["var"][m])
    else:
        mu, v = W[m] @ f64(c2["out"][m - 1]), f64(c2[f"var_{m}"])
    S = np.sqrt(v); al = mu / S
    return ndtr(al), np.exp(-0.5 * al * al) / SQ2PI / S


L = W[0].copy()
Lam = np.zeros((n, 0)); J2s = np.zeros(0); V = np.zeros((n, 0))
print(f"=== network {net}, jets = {jets}: chaos-state slices vs truth and vs chain (active set; energies noise-free)")
print("layer | N sources | D3: R2 chaos, R2 chain, best scale chaos | D21: R2 chaos, R2 chain, best scale chaos | "
      "D21: share of the chain's error explained by (chaos - chain) at coefficient c | cov: R2 of L L^T for the truth's "
      "off-diagonal covariance", flush=True)
for m in range(0, 15):
    Phi, J2 = state_jets(m)
    Wn = W[m + 1]
    V = Wn @ (Phi[:, None] * V)                                # conjugate the old sources by this transport
    Lam = np.concatenate([Lam, L.T], axis=1)                   # births: l_s = row a of L^m
    J2s = np.concatenate([J2s, J2])
    V = np.concatenate([V, Wn], axis=1)                        # legs of the newborn: columns of W_(m+1)
    L = Wn @ (Phi[:, None] * L)
    l = m + 1
    if l not in EVAL:
        continue
    U = L @ Lam
    VU = V * U
    D3c = 3.0 * (VU * U) @ J2s
    D21c = 2.0 * (VU * J2s[None, :]) @ U.T + ((U * U) * J2s[None, :]) @ V.T
    np.fill_diagonal(D21c, 0.0)
    mu1, v1 = f64(T["full"]["mu"][l]), f64(T["full"]["var"][l])
    act = (mu1 / np.sqrt(v1)) > -2.5
    mk = off & act[:, None] & act[None, :]
    out = [f"  {l:2d} | {Lam.shape[1]:5d}"]
    for key, X, tk, msk in (("D3", D3c, "k3", act), ("D21", D21c, "D21", mk)):
        ip = lambda A, B: float(np.sum((A * B)[msk]))
        Tf, T0, T1 = (f64(T[h][tk][l]) for h in ("full", "h0", "h1"))
        if key == "D21":
            Tf, T0, T1 = (t.copy() for t in (Tf, T0, T1))
            for t in (Tf, T0, T1):
                np.fill_diagonal(t, 0.0)
        ET = ip(T0, T1)
        r2 = lambda Y: 1.0 - ip(T0 - Y, T1 - Y) / ET
        sc = ip(X, Tf) / ip(X, X)
        cs = "   n/a"
        if l <= 14:
            Ch = f64(c2[f"D3_{l}"]) if key == "D3" else f64(c2[f"D21_{l}"]).copy()
            if key == "D21":
                np.fill_diagonal(Ch, 0.0)
            cs = f"{100 * r2(Ch):6.1f}"
            if key == "D21":
                dX = X - Ch; e0, e1, ef = T0 - Ch, T1 - Ch, Tf - Ch
                c = ip(ef, dX) / ip(dX, dX)
                sh = 1.0 - ip(e0 - c * dX, e1 - c * dX) / ip(e0, e1)
                extra = f" | chain-err share {100 * sh:5.1f}% (c {c:+.2f})"
        out.append(f"{key}: chaos {100 * r2(X):6.1f} chain {cs} scale {sc:+.2f}")
        if key == "D21" and l <= 14:
            out[-1] += extra
    C1 = f64(T["full"]["cov"][l]); C1 = 0.5 * (C1 + C1.T)
    G = L @ L.T
    ipc = lambda A, B: float(np.sum((A * B)[mk]))
    out.append(f"cov R2 {100 * (1 - ipc(C1 - G, C1 - G) / ipc(C1, C1)):5.1f}")
    print(" | ".join(out), flush=True)
