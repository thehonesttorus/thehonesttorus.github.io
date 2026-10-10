# Note XLIV section 9 (active-set version: neurons with alpha = mu / S > -2.5 only, the set the chain's saturation drop keeps) (stage-13 note, Theorem 5.4 / Cor. 5.5 / Remark 5.6, against our data).
#   python -I dil_ledger.py LOCDIR MC2DIR CD2DIR NET
# (A) Apparent dilation charge of each slice. For a common-scale mixture z = (1 + delta)(mu + y), Var delta = eps, the
#     first-order cumulants are (Ctilde = non-dilation covariance, approximated by the covariance itself):
#       C: eps (Ctilde + mu mu^T)    k3_a: 6 eps mu_a Ct_aa          D21_ac: 2 eps (2 mu_a Ct_ac + mu_c Ct_aa)
#       K22_ab: eps (4 Ct_aa Ct_bb + 8 Ct_ab^2)    K31_ab: 12 eps Ct_aa Ct_ab        k4_a: 12 eps Ct_aa^2.
#     eps_S = <f_S, S> / <f_S, f_S> is the least-squares charge each slice carries, for the truth (Monte Carlo) and for
#     the chain's own slices; for C it is the Rayleigh quotient along the mean, mu^T C mu / |mu|^4.
# (B) One-step mean-error ledger. At the true state the first-order Edgeworth closure moves E relu(z_a) by
#       dm = phi/(2S) dvar  -  R phi/(6 S^2) dk3  +  (R^2 - 1) phi/(24 S^3) dk4  +  Phi dmu,
#     R = mu/S; the three slice terms are evaluated with the chain's error (chain minus truth) at each layer. All second
#     moments of errors are noise-free (product of the chain-minus-truth differences against the two independent halves).
import sys, numpy as np
from scipy.special import ndtr
loc, mcd, cdd, net = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
W = np.load(f"{loc}/W_off{net}.npy")
c2 = np.load(f"{cdd}/chaindump2_{net}.npz")
mt = np.load(f"{loc}/truth_off{net}.npz")["m"].astype(np.float64)
out = c2["out"].astype(np.float64)
F = {}
for h in ("full", "h0", "h1"):
    z = np.load(f"{mcd}/mc2_off{net}_{h}.npz")
    F[h] = {k: z[k] for k in ("mu", "var", "k3", "k4", "cov", "D21", "K22", "K31")}
n = W.shape[1]
off = ~np.eye(n, dtype=bool)
sym = lambda A: 0.5 * (A + A.T)
f64 = lambda a: np.asarray(a, dtype=np.float64)
SQ2PI = np.sqrt(2 * np.pi)


ACT = None


def ip(A, B):
    m = off & ACT[:, None] & ACT[None, :]
    return float(np.sum(A[m] * B[m]))


def eps_of(f, S, offdiag):
    return (ip(f, S) / ip(f, f)) if offdiag else float(f[ACT] @ S[ACT] / (f[ACT] @ f[ACT]))


def forms(mu, Ct, v):
    s3 = 6.0 * mu * v
    f21 = 2.0 * (2.0 * mu[:, None] * Ct + mu[None, :] * v[:, None])
    f22 = 4.0 * np.outer(v, v) + 8.0 * Ct * Ct
    f31 = 12.0 * v[:, None] * Ct
    f4 = 12.0 * v * v
    return s3, f21, f22, f31, f4


print(f"=== network {net}", flush=True)
print("(A, active set) apparent charge eps_S x 1e3 (truth | chain), layer: Rayleigh | k3 | D21 | K22 | K31 | k4; active fraction", flush=True)
rows = []
for s in range(3, 15):
    mu_T = f64(F["full"]["mu"][s]); v_T = f64(F["full"]["var"][s]); C_T = sym(f64(F["full"]["cov"][s]))
    ACT = (mu_T / np.sqrt(v_T)) > -2.5
    mA = mu_T * ACT
    rho_T = float(mA @ C_T @ mA / (mA @ mA) ** 2)
    f3, f21, f22, f31, f4 = forms(mu_T, C_T, v_T)
    mu_C = f64(W[s]) @ out[s - 1]
    C_C = sym(f64(c2[f"C_off_{s}"])); np.fill_diagonal(C_C, f64(c2[f"var_{s}"]))
    mCA = mu_C * ACT
    rho_C = float(mCA @ C_C @ mCA / (mCA @ mCA) ** 2)
    sl_C = dict(k3=f64(c2[f"D3_{s}"]), D21=f64(c2[f"D21_{s}"]), K22=sym(f64(c2[f"wk4m_{s}"])), K31=f64(c2[f"wk431_{s}"]).T, k4=f64(c2[f"g4row_{s}"]))
    fm = dict(k3=(f3, False), D21=(f21, True), K22=(f22, True), K31=(f31, True), k4=(f4, False))
    et, ec = {}, {}
    for k, (f, od) in fm.items():
        a = [eps_of(f, f64(F[h][k][s]), od) for h in ("h0", "h1")]
        et[k] = 0.5 * (a[0] + a[1])
        ec[k] = eps_of(f, sl_C[k], od)
    rows.append((s, rho_T, rho_C, et, ec))
    print(f"  {s:2d} truth {1e3*rho_T:6.2f} | " + " | ".join(f"{1e3*et[k]:6.2f}" for k in fm)
          + f"     chain {1e3*rho_C:6.2f} | " + " | ".join(f"{1e3*ec[k]:6.2f}" for k in fm) + f"     active {ACT.mean():.3f}", flush=True)

