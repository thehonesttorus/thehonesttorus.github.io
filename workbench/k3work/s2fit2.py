# Follow-ups to s2fit.py on the same dumped old-tier content:
#  (1) the He-metric repair on the chain's OWN family. F = Sym^3(U) (+) {Phi_U(v I): v perp U} ("scalar return": one
#      n-vector beyond the chain's shared-basis tier). Frobenius projection onto F keeps P^(x3) T and the trace part
#      (tr D_k / r) I of the two-address-leg remainder; the He correction of THEORY_2 section 7 applies verbatim because
#      K(v) lies in F. Also the U-part-only repair (stays inside Sym^3(U)).
#  (2) the canonical split of the He-corrected two-address-leg state, T = U^(x3) H + Phi_U(D), U^T D = 0, D_i = d_i I/sqrt(r)
#      + D_i^o (shared-quadratic-actions note): singular spectrum of the traceless remainder and the reads after keeping
#      H and d exact and D^o at rank s.
#   python s2fit2.py NET L0 NL "128,320" "64,96,128" "8,16,32,64"
import sys, time, pickle, numpy as np
from s2fit import project, he_correct, reads_phi, ctr_bank, ctr_phi

net, L0, NL = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
r_halo = [int(v) for v in sys.argv[4].split(",")]
r_act = [int(v) for v in sys.argv[5].split(",")]
s_act = [int(v) for v in sys.argv[6].split(",")]
D = pickle.load(open(f"v29legs_off{net}.pkl", "rb")); legs = {g["layer"]: g for g in D["legs"]}
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64)
g0 = legs[L0]; k = g0["A"].shape[0]; n = g0["A"].shape[1]
sel = np.arange(min(k, L0 - 4))
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
bank = []
for s in src:
    bank.append((s["A"], 3.0 * s["P"] * s["w2"][None, :])); bank.append((s["P"], s["M"]))
def errs(get):
    out = []
    for m in layers:
        d3, d21 = get(m)
        out.append(f"L{m}:{np.linalg.norm(d3 - targets[m][0]) / np.linalg.norm(targets[m][0]):.3f}/"
                   f"{np.linalg.norm(d21 - targets[m][1]) / np.linalg.norm(targets[m][1]):.3f}")
    return " ".join(out)
def run_phi(U, B):
    Uc, Bc, rd = U, B, {}
    for m in layers:
        rd[m] = reads_phi(Uc, Bc)
        if m in G: Uc, Bc = G[m] @ Uc, np.einsum("ij,jab->iab", G[m], Bc, optimize=True)
    return errs(lambda m: rd[m])
al, be = 2 + 4 / n, 1 + 8 / n
stack = np.concatenate([np.concatenate([s["A"], s["P"], s["M"]], 1) for s in src], 1)
Usv = np.linalg.svd(stack, full_matrices=False)[0]
print(f"net {net} L0 {L0}: {len(sel)} old sources", flush=True)

# (1) the chain's family plus the scalar-return channel
for r in r_halo:
    U = Usv[:, :r]; P = U @ U.T; t0 = time.time()
    B0 = project(U, bank)                                   # two-address-leg Frobenius projection
    BP = np.einsum("ik,kab->iab", P, B0, optimize=True)     # neuron index in U: Phi_U(BP) = P^(x3) T (the chain's Sym^3)
    DQ = B0 - BP
    trq = np.einsum("iaa->i", DQ) / r
    I = np.eye(r)[None]
    variants = [("Sym3 (chain)", BP)]
    # Sym3 + He repair of the U part only (stays in Sym^3(U))
    eS = ctr_bank(bank) - ctr_phi(U, BP); eSU = U @ (U.T @ eS)
    variants.append(("Sym3 + He U-part", BP + (3 * be / (be * (r + 2) + 3 * al) * eSU)[:, None, None] * I))
    Bh = BP + trq[:, None, None] * I
    variants.append(("Sym3 + scalar return (Frob)", Bh))
    eh = ctr_bank(bank) - ctr_phi(U, Bh); ehU = U @ (U.T @ eh); ehQ = eh - ehU
    variants.append(("Sym3 + scalar return (He)", Bh + (3 * be / (be * (r + 2) + 3 * al) * ehU + 3 * be / (be * r + 3 * al) * ehQ)[:, None, None] * I))
    for tag, B in variants:
        print(f"  r {r:3d} {tag:28s} " + run_phi(U, B), flush=True)
    print(f"      ({time.time() - t0:.0f}s; scalar channel: one n-vector, transport n^2, reads ~ r/(2n) = {r / 2048:.2f} units per layer)", flush=True)
    del B0, BP, DQ, variants

# (2) canonical split and the traceless action spectrum of the He-corrected two-address-leg state
for r in r_act:
    U = Usv[:, :r]; P = U @ U.T; t0 = time.time()
    B = he_correct(U, project(U, bank), bank, al, be)
    BP = np.einsum("ik,kab->iab", P, B, optimize=True); Dm = B - BP
    d = np.einsum("iaa->i", Dm) / np.sqrt(r)
    Do = Dm - (d / np.sqrt(r))[:, None, None] * np.eye(r)[None]
    iu = np.triu_indices(r); wgt = np.where(iu[0] == iu[1], 1.0, np.sqrt(2.0))
    X = Do[:, iu[0], iu[1]] * wgt[None, :]                  # svec rows (Frobenius-isometric)
    Uu, lam, Vt = np.linalg.svd(X, full_matrices=False)
    E = lam ** 2; tot = E.sum()
    # tensor-norm ledger: |T|^2 = |H|^2 + (1/3) sum |D_i|^2 ; |Phi(BP)|^2 from the canonical identity
    Hn = np.sum(np.einsum("ia,ibc->abc", U, BP) ** 2)       # |L|^2 >= |Sym L|^2 (upper bound on the core norm)
    print(f"  r {r:3d}: traceless remainder energy (1/3)|D^o|^2 = {tot / 3:.3e}; trace channel (1/3)|d|^2 = {np.sum(d * d) / 3:.3e}; "
          f"core <= {Hn:.3e}", flush=True)
    print("      tail fraction of |D^o|^2 beyond s: " + "  ".join(f"s={s}:{E[s:].sum() / tot:.3f}" for s in (4, 8, 16, 32, 64, 128, 256, 512) if s < len(E)), flush=True)
    for s in s_act:
        Xs = (Uu[:, :s] * lam[:s]) @ Vt[:s]
        Dos = np.zeros_like(Do); Dos[:, iu[0], iu[1]] = Xs / wgt[None, :]; Dos = Dos + Dos.transpose(0, 2, 1)
        Dos[:, np.arange(r), np.arange(r)] *= 0.5
        Bs = BP + (d / np.sqrt(r))[:, None, None] * np.eye(r)[None] + Dos
        p = s + 1
        print(f"  r {r:3d} s {s:3d}  He S>=2, core+trace exact, D^o rank s  " + run_phi(U, Bs) +
              f"   [~units/layer: transport+read {2 * (r + p) / 1024:.2f}, write {(1024 * r + r ** 3 + p * r * r + 1024 * p) / 2 / 1024 ** 2:.2f}]", flush=True)
    print(f"  r {r:3d} full He S>=2 (s = all)                   " + run_phi(U, B) + f"  ({time.time() - t0:.0f}s)", flush=True)
