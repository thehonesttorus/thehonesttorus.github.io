# Two-address-leg source memory S_>=2(U) = Sym^3(U) + Sym^2(U) (x) U_perp (THEORY_2 of the conditional-innovation
# synthesis) measured on the chain's actual old-tier third-cumulant content, against the chain's own compression
# (every leg projected on a shared basis: Sym^3(U)).
#
#   T = sum_s sum_h [ 3 w2_h J(A_h, P_h) + J(P_h, M_h) ],  J(a, b) = Sym(a a b) / 3 normalised as in THEORY_2
#   Phi_U(B)_ijk = (u_i^T B_k u_j + u_i^T B_j u_k + u_j^T B_i u_k) / 3,   transport  G^(x3) Phi_U(B) = Phi_GU(G B)
#   projection   B = (3I - 2P) S,  S_k = sum_ij U_ia U_jb T_ijk = (1/3) sum_h [b_kh c c^T + a_kh (c d^T + d c^T)]
#   He metric    alpha |E|_F^2 + beta |ctr E|^2 (alpha = 2 + 4/n, beta = 1 + 8/n):  B* = B0 + [c_U (P e)_i + c_Q (Q e)_i] I_r
#
# Reads are the chain's: D3 (diagonal) and D21 (repeated-index slice, zero diagonal) at L0 .. L0+NL-1 after exact
# gated transport G = W_{l+1} diag(w1_l).  python s2fit.py NET L0 NL "32,48,64,96" [FRAMES]   (FRAMES: apm,ab)
import sys, pickle, time, numpy as np

def selftest():
    # Phi_U(B) from the projection equals Pi_2 T on an explicit tensor; transport identity; He correction optimality
    rng = np.random.default_rng(0); n, r, H = 9, 3, 5
    a, b = rng.standard_normal((n, H)), rng.standard_normal((n, H))
    T = np.zeros((n, n, n))
    for h in range(H):
        x, y = a[:, h], b[:, h]
        T += (np.einsum("i,j,k->ijk", x, x, y) + np.einsum("i,j,k->ijk", x, y, x) + np.einsum("i,j,k->ijk", y, x, x)) / 3
    U = np.linalg.qr(rng.standard_normal((n, r)))[0]; P = U @ U.T; I = np.eye(n)
    Pi = (np.einsum("ia,jb,kc->ijkabc", P, P, I) + np.einsum("ia,jb,kc->ijkabc", P, I, P) + np.einsum("ia,jb,kc->ijkabc", I, P, P)
          - 2 * np.einsum("ia,jb,kc->ijkabc", P, P, P))
    T0 = np.einsum("ijkabc,abc->ijk", Pi, T)
    B = project(U, [(a, b)])
    print("selftest projection:", np.abs(phi_full(U, B) - T0).max(), flush=True)
    G = rng.standard_normal((n, n))
    lhs = np.einsum("ia,jb,kc,abc->ijk", G, G, G, phi_full(U, B)); rhs = phi_full(G @ U, np.einsum("ij,jab->iab", G, B))
    print("selftest transport:", np.abs(lhs - rhs).max(), flush=True)
    d3, d21 = reads_phi(U, B); print("selftest reads:", np.abs(d3 - np.einsum("iii->i", T0)).max(),
                                      np.abs(d21 - (np.einsum("iij->ij", T0) - np.diag(np.einsum("iii->i", T0)))).max(), flush=True)
    al, be = 2 + 4 / n, 1 + 8 / n
    Bs = he_correct(U, B, [(a, b)], al, be)
    risk = lambda X: al * np.sum((T - X) ** 2) + be * np.sum(np.einsum("ijj->i", T - X) ** 2)
    base = risk(phi_full(U, Bs)); worst = 0.0
    for _ in range(20):
        D = rng.standard_normal(Bs.shape); D = D + D.transpose(0, 2, 1)
        worst = min(worst, risk(phi_full(U, Bs + 1e-3 * D)) - base)
    print(f"selftest He optimum: risk {base:.6f}, frobenius projection {risk(phi_full(U, B)):.6f}, min perturbation change {worst:.2e} (>= 0 expected)", flush=True)

def phi_full(U, B):
    t = np.einsum("ia,kab,jb->ijk", U, B, U)
    return (t + t.transpose(0, 2, 1) + t.transpose(2, 0, 1)) / 3     # u_i B_k u_j + u_i B_j u_k + u_j B_i u_k

def project(U, bank):
    # bank: list of (a, b) leg arrays (n, H) for J(a_h, b_h); returns B (n, r, r) of Pi_2 T in the frame U (orthonormal)
    n, r = U.shape; S = np.zeros((n, r * r))
    for a, b in bank:
        c, d = U.T @ a, U.T @ b                                          # (r, H)
        cc = np.einsum("ah,bh->hab", c, c).reshape(c.shape[1], r * r)
        cd = np.einsum("ah,bh->hab", c, d); cd = (cd + cd.transpose(0, 2, 1)).reshape(c.shape[1], r * r)
        S += (b @ cc + a @ cd) / 3.0
    S = S.reshape(n, r, r)
    return 3.0 * S - 2.0 * np.einsum("ka,abc->kbc", U, np.einsum("la,lbc->abc", U, S))   # (3I - 2P) on the free index

def ctr_bank(bank):
    # ctr(T)_i = sum_j T_ijj for T = sum_h J(a_h, b_h):  (2 a_i (a.b) + b_i |a|^2) / 3
    out = 0.0
    for a, b in bank:
        out = out + (2.0 * a @ np.einsum("ih,ih->h", a, b) + b @ np.einsum("ih,ih->h", a, a)) / 3.0
    return out

def ctr_phi(U, B):
    g = np.einsum("jab,jb->a", B, U)                                    # sum_j B_j u_j
    M = U.T @ U
    return (2.0 * U @ g + np.einsum("iab,ba->i", B, M)) / 3.0

def he_correct(U, B0, bank, al, be):
    n, r = U.shape
    e = ctr_bank(bank) - ctr_phi(U, B0); eU = U @ (U.T @ e); eQ = e - eU
    cU, cQ = 3 * be / (be * (r + 2) + 3 * al), 3 * be / (be * r + 3 * al)
    return B0 + (cU * eU + cQ * eQ)[:, None, None] * np.eye(r)[None]

def reads_phi(U, B):
    n, r = U.shape
    d3 = np.einsum("ia,iab,ib->i", U, B, U)
    Z = np.einsum("ia,ib->iab", U, U).reshape(n, r * r)
    q = np.einsum("iab,ib->ia", B, U)
    d21 = (Z @ B.reshape(n, r * r).T + 2.0 * q @ U.T) / 3.0
    np.fill_diagonal(d21, 0.0)
    return d3, d21

if __name__ == "__main__":
    selftest()
    net, L0, NL = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rs = [int(v) for v in sys.argv[4].split(",")]
    frames = (sys.argv[5] if len(sys.argv) > 5 else "apm,ab").split(",")
    D = pickle.load(open(f"v29legs_off{net}.pkl", "rb")); legs = {g["layer"]: g for g in D["legs"]}
    Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64)
    g0 = legs[L0]; k = g0["A"].shape[0]; n = g0["A"].shape[1]
    sel = np.arange(min(k, L0 - 4))                                     # the chain's shared-basis tier (age > 4)
    def m_leg(g, j): return g["P"][j] * g["s"][j][None, :] + 3.0 * g["A"][j] * g["e"][j][None, :] + g["Z"][j] @ g["L"][j].T
    src = [dict(A=g0["A"][j].astype(np.float64), P=g0["P"][j].astype(np.float64), M=m_leg(g0, j).astype(np.float64),
                w2=g0["w2b"][j].astype(np.float64)) for j in sel]
    def reads(srcs):
        D3 = np.zeros(n); D21 = np.zeros((n, n))
        for s in srcs:
            A, P, M, w2 = s["A"], s["P"], s["M"], s["w2"]
            D3 += 3.0 * ((A * A * P) @ w2) + (M * P * P).sum(1)
            D21 += ((A * A * w2) @ P.T + 2.0 * (A * P * w2) @ A.T + (2.0 / 3.0) * (M * P) @ P.T + (1.0 / 3.0) * (P * P) @ M.T)
        np.fill_diagonal(D21, 0.0)
        return D3, D21
    layers = [L0 + i for i in range(NL)]
    G, targets, cur = {}, {}, src
    for m in layers:
        targets[m] = reads(cur)
        if m + 1 in layers:
            G[m] = Wcol[m + 1] * legs[m]["w1"].astype(np.float64)[None, :]
            cur = [dict(A=G[m] @ s["A"], P=G[m] @ s["P"], M=G[m] @ s["M"], w2=s["w2"]) for s in cur]
    # the bank in J form: J(A, 3 w2 P) and J(P, M)
    bank = []
    for s in src:
        bank.append((s["A"], 3.0 * s["P"] * s["w2"][None, :]))
        bank.append((s["P"], s["M"]))
    # consistency of the J form with the chain's reads at L0
    Tc = np.zeros(n)
    for a, b in bank: Tc += (a * a * b).sum(1)
    print(f"net {net} L0 {L0}: {len(sel)} old sources; J-form D3 vs chain-form D3 max rel diff "
          f"{np.abs(Tc - targets[L0][0]).max() / np.abs(targets[L0][0]).max():.2e}", flush=True)
    def errs(get):
        out = []
        for m in layers:
            d3, d21 = get(m)
            out.append(f"L{m}:{np.linalg.norm(d3 - targets[m][0]) / np.linalg.norm(targets[m][0]):.3f}/"
                       f"{np.linalg.norm(d21 - targets[m][1]) / np.linalg.norm(targets[m][1]):.3f}")
        return " ".join(out)
    stack = {"apm": np.concatenate([np.concatenate([s["A"], s["P"], s["M"]], 1) for s in src], 1),
             "ab": np.concatenate([a * np.sqrt(np.linalg.norm(b, axis=0))[None, :] for a, b in bank], 1)}
    al, be = 2 + 4 / n, 1 + 8 / n
    for fr in frames:
        Usv = np.linalg.svd(stack[fr], full_matrices=False)[0]
        for r in rs + ([320] if fr == "apm" else []):
            U = Usv[:, :r]; t0 = time.time()
            # chain-style yardstick: every leg on U (Sym^3(U))
            cur = [dict(A=U @ (U.T @ s["A"]), P=U @ (U.T @ s["P"]), M=U @ (U.T @ s["M"]), w2=s["w2"]) for s in src]
            ys = {}
            for m in layers:
                ys[m] = reads(cur)
                if m in G: cur = [dict(A=G[m] @ s["A"], P=G[m] @ s["P"], M=G[m] @ s["M"], w2=s["w2"]) for s in cur]
            print(f"  frame {fr} r {r:3d}  Sym3(U)  [chain kind]   " + errs(lambda m: ys[m]), flush=True)
            if r > 160: continue
            B0 = project(U, bank)
            for tag, B in (("S>=2 frobenius", B0), ("S>=2 He-metric", he_correct(U, B0, bank, al, be))):
                Uc, Bc, rd = U, B, {}
                for m in layers:
                    rd[m] = reads_phi(Uc, Bc)
                    if m in G: Uc, Bc = G[m] @ Uc, np.einsum("ij,jab->iab", G[m], Bc, optimize=True)
                print(f"  frame {fr} r {r:3d}  {tag:15s} " + errs(lambda m: rd[m]), flush=True)
            q = r * (r + 1) // 2
            print(f"      cost per layer in units of 2n^3: transport {q / n:.2f}, D21 read {q / n:.2f}, one graduating source "
                  f"{4 * q / n:.2f}; Sym3 at r: ~{r / n * 3:.2f} per source  ({time.time() - t0:.0f}s)", flush=True)
