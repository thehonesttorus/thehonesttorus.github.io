# Representability of the chain's summed third-cumulant tensor by a matrix-trace carrier
#     Theta_ijk = Sym tau(a_i a_j b_k),   a_i, b_i real symmetric r x r,
# measured on the reads the chain uses (D3 = diagonal, D21 = repeated-index slices), at the dump layer and after
# exact transport to later layers (a_i -> sum_k G_ik a_k, G = W_{l+1} diag(w1_l), the product-gate response).
#
#   python mt_fit.py NET L0 R [SUBSET] [ITERS]        SUBSET: old (age > 4, the shared-basis tier) | all
# needs v29legs_off{NET}.pkl from run_v29legs.py NET "L0,L0+1,L0+2,L0+3" (with the w1 dump).
import sys, pickle, time, numpy as np, torch

net, L0, R = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
subset = sys.argv[4] if len(sys.argv) > 4 else "old"
iters = int(sys.argv[5]) if len(sys.argv) > 5 else 3000
torch.set_num_threads(int(__import__("os").environ.get("MT_THREADS", "8")))
D = pickle.load(open(f"v29legs_off{net}.pkl", "rb"))
legs = {g["layer"]: g for g in D["legs"]}
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64)    # chain's W at layer li is Wcol[li]
g0 = legs[L0]
k = g0["A"].shape[0]
# the chain confines min(k, L0 - AGE_OLD) leading slots (oldest first) at layer L0: that is the shared-basis tier
sel = np.arange(min(k, L0 - 4)) if subset == "old" else np.arange(k)
n = g0["A"].shape[1]
rres = g0["Z"].shape[2] - 2

def m_leg(g, j):
    # M = P diag(s) + 3 A diag(e) + Z L^T (thin residual; the factor 3 is folded into Z at birth)
    return g["P"][j] * g["s"][j][None, :] + 3.0 * g["A"][j] * g["e"][j][None, :] + g["Z"][j] @ g["L"][j].T

src = [dict(A=g0["A"][j].astype(np.float64), P=g0["P"][j].astype(np.float64), M=m_leg(g0, j).astype(np.float64),
            w2=g0["w2b"][j].astype(np.float64)) for j in sel]

def reads(srcs):
    # Theta = sum_s sum_j [w2_j Sym(A_j, A_j, P_j) + (1/3) Sym(M_j, P_j, P_j)]  (the chain's (3,) and (2,1) programs)
    D3 = np.zeros(n); D21 = np.zeros((n, n))
    for s in srcs:
        A, P, M, w2 = s["A"], s["P"], s["M"], s["w2"]
        D3 += 3.0 * ((A * A * P) @ w2) + (M * P * P).sum(1)
        D21 += ((A * A * w2) @ P.T + 2.0 * (A * P * w2) @ A.T + (2.0 / 3.0) * (M * P) @ P.T + (1.0 / 3.0) * (P * P) @ M.T)
    np.fill_diagonal(D21, 0.0)
    return D3, D21

# exact CP transport of the selected content to L0+1..L0+3 (no newborns: this is the carried state's own future)
NL = int(__import__("os").environ.get("MT_NL", "4")); NTR = int(__import__("os").environ.get("MT_NTRAIN", "2"))
layers = [L0 + i for i in range(NL)]
G = {}
targets = {}
cur = src
for m in layers:
    targets[m] = reads(cur)
    if m + 1 in layers:
        w1 = legs[m]["w1"].astype(np.float64)
        Gm = Wcol[m + 1] * w1[None, :]
        G[m] = Gm
        cur = [dict(A=Gm @ s["A"], P=Gm @ s["P"], M=Gm @ s["M"], w2=s["w2"]) for s in cur]

# sanity: the chain's own dumped D3 at L0 (all sources, plus thin feedback / feed terms) against the CP rebuild
D3all, _ = reads([dict(A=g0["A"][j].astype(np.float64), P=g0["P"][j].astype(np.float64), M=m_leg(g0, j).astype(np.float64),
                       w2=g0["w2b"][j].astype(np.float64)) for j in range(k)])
cc = np.corrcoef(D3all, g0["D3"].astype(np.float64))[0, 1]
print(f"net {net} L0 {L0} k {k} subset {subset} ({len(sel)} sources) r {R}: CP rebuild vs chain D3 corr {cc:.3f}", flush=True)

T = {m: (torch.tensor(targets[m][0], dtype=torch.float32), torch.tensor(targets[m][1], dtype=torch.float32)) for m in layers}
Gt = {m: torch.tensor(G[m], dtype=torch.float32) for m in G}
scale = {m: (T[m][0].pow(2).sum(), T[m][1].pow(2).sum()) for m in layers}
train = layers[:NTR]
test = layers[NTR:]

torch.manual_seed(0)
s0 = float(np.abs(targets[L0][0]).mean() ** (1.0 / 3.0)) / R ** 0.5
INIT = __import__("os").environ.get("MT_INIT", "rand")
if INIT == "svd":
    # a_i, b_i in the leading subspace of the content: a_i = sum_m Q_im E_m with E_m an orthonormal basis of Sym_r,
    # Q = top r(r+1)/2 left singular vectors of the stacked legs at L0 (scaled to the D3 magnitude)
    Dm = min(R * (R + 1) // 2, n)        # at most n independent directions in the content
    st = np.concatenate([np.concatenate([s["A"], s["P"], s["M"]], 1) for s in src], 1)
    Q = np.linalg.svd(st, full_matrices=False)[0][:, :Dm]
    E = np.zeros((Dm, R, R)); t = 0
    for p_ in range(R):
        for q_ in range(p_, R):
            if t == Dm: break
            if p_ == q_: E[t, p_, p_] = 1.0
            else: E[t, p_, q_] = E[t, q_, p_] = 2 ** -0.5
            t += 1
    rng = np.random.default_rng(0)
    Ra, Rb = (np.linalg.qr(rng.standard_normal((Dm, Dm)))[0] for _ in range(2))
    a0 = np.einsum("im,mpq->ipq", Q @ Ra, E); b0 = np.einsum("im,mpq->ipq", Q @ Rb, E)
    Ua = torch.nn.Parameter(torch.tensor(a0, dtype=torch.float32))
    Ub = torch.nn.Parameter(torch.tensor(b0, dtype=torch.float32))
else:
    Ua = torch.nn.Parameter(s0 * torch.randn(n, R, R))
    Ub = torch.nn.Parameter(s0 * torch.randn(n, R, R))

def mt_reads(a, b):
    X = a @ a
    Y = b @ a + a @ b
    D3 = torch.einsum("ipq,iqp->i", X, b)
    D21 = (X.reshape(n, -1) @ b.reshape(n, -1).T + Y.reshape(n, -1) @ a.reshape(n, -1).T) / 3.0
    D21 = D21 - torch.diag(torch.diagonal(D21))
    return D3, D21

if INIT == "svd":
    with torch.no_grad():
        d3 = mt_reads(0.5 * (Ua + Ua.transpose(1, 2)), 0.5 * (Ub + Ub.transpose(1, 2)))[0]
        lam = (T[L0][0].abs().mean() / d3.abs().mean().clamp_min(1e-30)) ** (1.0 / 3.0)
        Ua.mul_(lam); Ub.mul_(lam)
        s0 = float(Ua.abs().mean())

def all_reads():
    a = 0.5 * (Ua + Ua.transpose(1, 2)); b = 0.5 * (Ub + Ub.transpose(1, 2))
    out = {}
    for m in layers:
        out[m] = mt_reads(a, b)
        if m in Gt:
            a = torch.einsum("ik,kpq->ipq", Gt[m], a); b = torch.einsum("ik,kpq->ipq", Gt[m], b)
    return out

def rel(out, m):
    e3 = ((out[m][0] - T[m][0]).pow(2).sum() / scale[m][0]).item()
    e21 = ((out[m][1] - T[m][1]).pow(2).sum() / scale[m][1]).item()
    return e3 ** 0.5, e21 ** 0.5

opt = torch.optim.Adam([Ua, Ub], lr=0.03 * s0)
sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, iters + 1)
t0 = time.time()
for it in range(iters + 1):
    opt.zero_grad()
    out = all_reads()
    loss = sum((out[m][0] - T[m][0]).pow(2).sum() / scale[m][0] + (out[m][1] - T[m][1]).pow(2).sum() / scale[m][1] for m in train)
    loss.backward()
    opt.step()
    sched.step()
    if it % 500 == 0 or it == iters:
        with torch.no_grad():
            out = all_reads()
            msg = " ".join(f"L{m}:{rel(out, m)[0]:.3f}/{rel(out, m)[1]:.3f}" for m in layers)
        print(f"  it {it:5d} loss {loss.item():.4f}  rel err D3/D21  {msg}  ({time.time() - t0:.0f}s)", flush=True)

# yardstick: the chain's kind of compression, all selected legs projected on one shared rank-Dk basis at L0
# (top left singular vectors of the stacked legs), then transported exactly
stack = np.concatenate([np.concatenate([s["A"], s["P"], s["M"]], 1) for s in src], 1)
Usv, _, _ = np.linalg.svd(stack, full_matrices=False)
for Dk in (64, 128, 320):
    Q = Usv[:, :Dk]
    cur = [dict(A=Q @ (Q.T @ s["A"]), P=Q @ (Q.T @ s["P"]), M=Q @ (Q.T @ s["M"]), w2=s["w2"]) for s in src]
    errs = []
    for m in layers:
        d3, d21 = reads(cur)
        errs.append(f"L{m}:{np.linalg.norm(d3 - targets[m][0]) / np.linalg.norm(targets[m][0]):.3f}/"
                    f"{np.linalg.norm(d21 - targets[m][1]) / np.linalg.norm(targets[m][1]):.3f}")
        if m in G:
            cur = [dict(A=G[m] @ s["A"], P=G[m] @ s["P"], M=G[m] @ s["M"], w2=s["w2"]) for s in cur]
    print(f"  yardstick shared basis rank {Dk}: " + " ".join(errs), flush=True)
print(f"params: matrix-trace {2 * n * R * (R + 1) // 2}, CP legs {3 * n * n * len(src)}; per-layer transport+reads ~ "
      f"{(4 * n * n * R * R + 3 * n * R ** 3) / 2.15e9:.2f} units (one unit = 2 n^3)", flush=True)
