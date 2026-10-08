# The five-component audit of the D21 table, ported to the chain's own tables (note XXXIX; the capture-audit and
# depth-profile checkpoints of the parallel line measured the physical table only).
#   python d21audit.py NET CHAIN_DUMP MCFILE [MC_H0 MC_H1]      (SLICE=K22 or SLICE=K31 audits the fourth-cumulant slices:
#   the chain's wk4m against the truth's K22, the chain's wk431 against the truth's K31 transposed, as the V37 oracles use them)
# With the two independent Monte Carlo halves, every error energy is also reported noise-corrected: the noise energy of the
# full-sample truth is |comp(D_h0 - D_h1)|^2 / 4, subtracted component by component ("nc" columns).
# D21[i, c] = kappa(z_i, z_i, z_c), zero diagonal. Off the diagonal every table splits exactly and orthogonally into
#   S0  the mean of the symmetric part;
#   S1  alpha_i + alpha_j (the additive symmetric part, from the row sums);
#   S2  the symmetric remainder (zero row sums);
#   A1  beta_i - beta_j (the additive antisymmetric part);
#   A2  the antisymmetric remainder (zero row sums: a circulation).
# S0 + S1 + A1 is rank <= 3 (u 1^T + 1 v^T + c 1 1^T). Reported per layer, for the truth (Monte Carlo), the chain and the
# chain's error: the components' shares of energy, plain and in the read metric w_ic = rho_i Phi_c (rho = phi(alpha)/sigma,
# the weight with which the covariance program reads row i), the spectral capture of S2 and A2, and what the production
# feedback's rank-2 compression keeps of each component (best rank 2, an upper bound for the range finder).
import sys, numpy as np
from math import erf

net, chain_f, mc_f = int(sys.argv[1]), sys.argv[2], sys.argv[3]
ch = np.load(chain_f); mc = np.load(mc_f)
halves = [np.load(f) for f in sys.argv[4:6]] if len(sys.argv) >= 6 else None
import os
SL = os.environ.get("SLICE", "D21")
CKEY = {"D21": "D21", "K22": "wk4m", "K31": "wk431"}[SL]
TR = (lambda x: x.T) if SL == "K31" else (lambda x: x)
n = mc["mu"].shape[1]
RANKS = (2, 8, 16, 32, 64, 128)


def zd(x):
    x = np.array(x, dtype=np.float64); np.fill_diagonal(x, 0.0); return x


def split(D):
    S = (D + D.T) / 2; A = (D - D.T) / 2
    off = ~np.eye(n, dtype=bool)
    c = S[off].mean()
    rs = S.sum(1)                                   # off-diagonal row sums (diagonal is zero)
    al = (rs - (n - 1) * c) / (n - 2)
    S0 = zd(np.full((n, n), c)); S1 = zd(al[:, None] + al[None, :]); S2 = zd(S - S0 - S1)
    be = A.sum(1) / n
    A1 = zd(be[:, None] - be[None, :]); A2 = zd(A - A1)
    return dict(S0=S0, S1=S1, S2=S2, A1=A1, A2=A2)


def nrm(x, w=None):
    return float(np.sum((x if w is None else x * w) ** 2))


def cap(X, ranks=RANKS):
    s = np.linalg.svd(X, compute_uv=False) ** 2
    t = s.sum(); c = np.cumsum(s) / t
    return [float(c[min(r, len(c)) - 1]) for r in ranks]


def phi(a):
    return np.exp(-0.5 * a * a) / np.sqrt(2 * np.pi)


Phi = np.vectorize(lambda a: 0.5 * (1 + erf(a / np.sqrt(2))))
print(f"net {net} slice {SL}: chain {chain_f}, truth {mc_f}")
print("layer | energy shares S0 S1 S2 A1 A2: truth / chain / chain error (error energy over truth energy) | read-weighted error"
      " shares | rank-2 keeps of the chain table (S1+A1, S2+A2) | S2, A2 truth capture at ranks " + "/".join(map(str, RANKS)))
for l in range(1, 16):
    key = f"{CKEY}_{l}"
    if key not in ch:
        continue
    Dt = zd(TR(mc[SL][l])); Dc = zd(ch[key])
    sd = np.sqrt(np.maximum(mc["var"][l].astype(np.float64), 1e-30)); al = mc["mu"][l].astype(np.float64) / sd
    w = (phi(al) / sd)[:, None] * Phi(al)[None, :]
    ct, cc = split(Dt), split(Dc)
    ce = {k: cc[k] - ct[k] for k in ct}
    Et, Ec, Ee = nrm(Dt), nrm(Dc), nrm(Dc - Dt)
    Etw, Eew = nrm(Dt, w), nrm(Dc - Dt, w)
    sh = lambda comp, E, ww=None: " ".join(f"{nrm(comp[k], ww) / E:.3f}" for k in ("S0", "S1", "S2", "A1", "A2"))
    U, s, Vt = np.linalg.svd(Dc)
    P2 = U[:, :2] @ (U[:, :2].T @ Dc @ Vt[:2].T) @ Vt[:2]          # best rank 2 of the chain table
    pc = split(P2)
    add_keep = (nrm(pc["S1"] + pc["A1"] + pc["S0"]) / max(nrm(cc["S1"] + cc["A1"] + cc["S0"]), 1e-300))
    flat_keep = nrm(pc["S2"] + pc["A2"]) / max(nrm(cc["S2"] + cc["A2"]), 1e-300)
    if halves is not None:
        dn = split(zd(TR(halves[0][SL][l])) - zd(TR(halves[1][SL][l])))
        ncs = " ".join(f"{(nrm(ce[k]) - nrm(dn[k]) / 4) / Et:+.4f}" for k in ("S0", "S1", "S2", "A1", "A2"))
        ncw = " ".join(f"{(nrm(ce[k], w) - nrm(dn[k], w) / 4) / Etw:+.4f}" for k in ("S0", "S1", "S2", "A1", "A2"))
        print(f"{l:2d} nc | chain error minus noise, plain: {ncs} | read-weighted: {ncw} | noise/error plain "
              f"{sum(nrm(dn[k]) for k in dn) / 4 / Ee:.3f}", flush=True)
    print(f"{l:2d} | {sh(ct, Et)} / {sh(cc, Ec)} / {sh(ce, Et)} (err {Ee / Et:.4f}) | {sh(ce, Etw, w)} (err {Eew / Etw:.4f}) | "
          f"{add_keep:.3f} {flat_keep:.3f} | " + " ".join(f"{x:.3f}" for x in cap(ct["S2"])) + " ; "
          + " ".join(f"{x:.3f}" for x in cap(ct["A2"])), flush=True)
