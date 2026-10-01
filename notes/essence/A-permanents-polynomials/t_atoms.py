"""Team A tests on the memory atoms of FC at n = 1024 (bench w1024_d16).
At target layers t, the old sources (age >= 3) define D_old = sum_s sum_r w_r [Y_ra^2 Z_rb + 2 Y_ra Z_ra Y_rb].
We measure the share of ||D_old|| captured by:
  diag  : the diagonal path pairing (Bethe / annealed-in-intermediate-indices value): Y o Y -> (SP o SP)(W o W)(G o G)...
  mf    : the rank-one mean-field of the squared legs (row sums x column sums; trace-channel-like)
  rank k: leg propagator truncated to its top-k singular directions (Perron / low-degree-around-reference expansion)
and the single-sample relative variance of the unbiased Godsil-Gutman-type sketch with cube-root-of-unity signs.
Usage: python t_atoms.py MLP [targets]"""
import sys, json, time
import numpy as np
sys.path.insert(0, "../../fresh-slate/bench")
import bench
import fc_hook as F

mlp = int(sys.argv[1]) if len(sys.argv) > 1 else 0
targets = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [6, 10, 14]
S = bench.load_set("w1024_d16"); W = bench.weights(S, mlp)
F.DIAG = True
rng = np.random.default_rng(1)
out = {}
NSK = int(sys.argv[3]) if len(sys.argv) > 3 else 16

def rel(A, B):
    return float(np.linalg.norm(A - B) / np.linalg.norm(B))

def hook(t, sources, D21):
    if t not in targets:
        return False
    t0 = time.time()
    n = D21.shape[0]
    old = [s for s in sources if t - s["s"] >= 3]
    D = np.zeros((n, n)); Dd = np.zeros((n, n)); Dmf = np.zeros((n, n)); Dk = {k: np.zeros((n, n)) for k in (1, 16, 128, 256)}; Dmk = {k: np.zeros((n, n)) for k in (64, 128, 256)}
    U = np.zeros(n, complex); V = np.zeros(n, complex)
    Ys, Zs, ws = [], [], []
    for s in old:
        w = s["w2"].astype(np.float32); Z = s["Z"]; SP = s["SP"].astype(np.float32)
        Y = SP @ Z
        D += (((Y * Y) * w[:, None]).T @ Z + 2 * (((Y * Z) * w[:, None]).T @ Y)).astype(np.float64)
        Y2d = (SP * SP) @ s["Z2"]; YZd = np.diag(SP)[:, None] * s["Z2"]
        Dd += (((Y2d * w[:, None]).T @ Z) + 2 * ((YZd * w[:, None]).T @ Y)).astype(np.float64)
        Y2 = (Y * Y).astype(np.float64); YZ = (Y * Z).astype(np.float64)
        Y2m = np.outer(Y2.sum(1), Y2.sum(0)) / Y2.sum()
        YZm = np.outer(YZ.sum(1), YZ.sum(0)) / YZ.sum()
        Dmf += ((Y2m * w[:, None]).T @ Z + 2 * ((YZm * w[:, None]).T @ Y))
        Y2 = (Y * Y).astype(np.float64)
        u_, sv, vt = np.linalg.svd(Z.astype(np.float64), full_matrices=False)
        for k in Dk:
            Zk = ((u_[:, :k] * sv[:k]) @ vt[:k]).astype(np.float32); Yk = SP @ Zk
            Dk[k] += (((Yk * Yk) * w[:, None]).T @ Zk + 2 * (((Yk * Zk) * w[:, None]).T @ Yk)).astype(np.float64)
            if k in Dmk:   # Bethe/mean-field reference at full rank + rank-k legs for the fluctuation only
                Y2k = (Yk * Yk).astype(np.float64); YZk = (Yk * Zk).astype(np.float64)
                F2 = Y2k - np.outer(Y2k.sum(1), Y2k.sum(0)) / Y2k.sum(); FZ = YZk - np.outer(YZk.sum(1), YZk.sum(0)) / YZk.sum()
                Dmk[k] += (((Y2m * w[:, None]).T @ Z + 2 * ((YZm * w[:, None]).T @ Y))
                           + ((F2 * w[:, None]).T @ Zk + 2 * ((FZ * w[:, None]).T @ Yk)))
        for k in (16,):
            pass
        Ys.append(Y); Zs.append(Z); ws.append(w)
    Dall = D21
    # Godsil-Gutman sketch: xi_r cube roots of unity; u = sum xi y, v = sum w conj(xi)^2 z
    om = np.exp(2j * np.pi / 3)
    e1 = []
    for k in range(NSK):
        u = np.zeros(n, complex); v = np.zeros(n, complex)
        for Y, Z, w in zip(Ys, Zs, ws):
            xi = om ** rng.integers(0, 3, size=n)
            u += Y.T.astype(np.float64) @ xi; v += Z.T.astype(np.float64) @ (w * np.conj(xi) ** 2)
        Dh = np.real(np.outer(u * u, v) + 2 * np.outer(u * v, u))
        e1.append(float(np.sum((Dh - D) ** 2) / np.sum(D ** 2)))
    r = dict(t=t, n_old_sources=len(old), atoms=len(old) * n,
             share_old=float(np.linalg.norm(D) / np.linalg.norm(Dall)),
             diag_pairing_relerr=rel(Dd, D), diag_cos=float(np.sum(Dd * D) / np.linalg.norm(Dd) / np.linalg.norm(D)),
             meanfield_relerr=rel(Dmf, D),
             rank_relerr={k: rel(Dk[k], D) for k in Dk}, mf_plus_rank_relerr={k: rel(Dmk[k], D) for k in Dmk},
             sketch_single_relvar=float(np.mean(e1)), sketch_K_for_5pct=float(np.mean(e1) / 0.05 ** 2),
             secs=time.time() - t0)
    print(json.dumps(r), flush=True)
    out[t] = r
    return t >= max(targets)

F.HOOK = hook
t0 = time.time()
F.run(W, k4mf=True, slices=False, readout=True)
json.dump(out, open(f"results/atoms_mlp{mlp}.json", "w"), indent=1)
print("total", time.time() - t0)
