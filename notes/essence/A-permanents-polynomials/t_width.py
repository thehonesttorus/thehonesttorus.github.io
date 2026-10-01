"""Width test of Conjecture 1 and the single-cut Bethe curve.
For old sources (age >= 3) at target t: relative Frobenius error of
  bethe : diagonal path pairing at every intermediate layer
  mf    : rank-one mean-field of the squared legs
  cut_k : exact legs except a forced diagonal pairing at ONE layer, k layers before the target
          (Y o Y -> (Y' o Y') (J o J), Y' = legs at t-k, J = exact last-k gated product)
  rank f: legs truncated to rank f*n
Usage: python t_width.py SET MLP targets"""
import sys, json, time
import numpy as np
sys.path.insert(0, "../../fresh-slate/bench")
import bench
import fc_hook as F

name, mlp = sys.argv[1], int(sys.argv[2]); targets = [int(x) for x in sys.argv[3].split(",")]
KS = (1, 2, 3, 5)
S = bench.load_set(name); W = bench.weights(S, mlp)
F.DIAG = True; F.GS.clear()
snap = {}; out = {}

def rel(A, B):
    return float(np.linalg.norm(A - B) / np.linalg.norm(B))

def D_of(Y, Z, w, Y2=None, YZ=None):
    Y2 = Y * Y if Y2 is None else Y2; YZ = Y * Z if YZ is None else YZ
    return (((Y2 * w[:, None]).T @ Z) + 2 * ((YZ * w[:, None]).T @ Y)).astype(np.float64)

def hook(t, sources, D21):
    for tt in targets:
        if tt - t in KS:
            snap[(tt, tt - t)] = {s["s"]: s["Z"].copy() for s in sources}
    if t not in targets:
        return False
    n = D21.shape[0]
    old = [s for s in sources if t - s["s"] >= 3]
    D = 0; Db = 0; Dm = 0; Dc = {k: 0 for k in KS}; Dr = {f: 0 for f in (8, 4)}; chk = []
    for k in KS:
        J = F.GS[t - k]
        for i in range(t - k + 1, t):
            J = J @ F.GS[i]
        Dc[k] = (J, 0)
    for s in old:
        w = s["w2"].astype(np.float32); Z = s["Z"]; SP = s["SP"].astype(np.float32)
        Y = SP @ Z
        D = D + D_of(Y, Z, w)
        Db = Db + D_of(Y, Z, w, (SP * SP) @ s["Z2"], np.diag(SP)[:, None] * s["Z2"])
        Y2 = (Y * Y).astype(np.float64); YZ = (Y * Z).astype(np.float64)
        Dm = Dm + D_of(Y, Z, w, np.outer(Y2.sum(1), Y2.sum(0)) / Y2.sum(), np.outer(YZ.sum(1), YZ.sum(0)) / YZ.sum())
        for k in KS:
            J, acc = Dc[k]
            if t - s["s"] - k >= 1 and s["s"] in snap.get((t, k), {}):
                Zp = snap[(t, k)][s["s"]]; Yp = SP @ Zp; JJ = J * J
                if len(chk) < 1: chk.append(rel(Zp @ J, Z))
                acc = acc + D_of(Y, Z, w, (Yp * Yp) @ JJ, (Yp * Zp) @ JJ)
            else:
                acc = acc + D_of(Y, Z, w)
            Dc[k] = (J, acc)
        u_, sv, vt = np.linalg.svd(Z.astype(np.float64), full_matrices=False)
        for f in Dr:
            k = n // f; Zk = ((u_[:, :k] * sv[:k]) @ vt[:k]).astype(np.float32)
            Dr[f] = Dr[f] + D_of(SP @ Zk, Zk, w)
    r = dict(set=name, mlp=mlp, t=t, n=n, n_old=len(old), share_old=float(np.linalg.norm(D) / np.linalg.norm(D21)),
             bethe=rel(Db, D), mf=rel(Dm, D), cut={k: rel(Dc[k][1], D) for k in KS}, cut1_vs_bethe=rel(Dc[1][1], Db), cut1_vs_mf=rel(Dc[1][1], Dm), bethe_vs_mf=rel(Db, Dm),
             rank={f"n/{f}": rel(Dr[f], D) for f in Dr}, transport_check=chk)
    print(json.dumps(r), flush=True); out[t] = r
    return t >= max(targets)

F.HOOK = hook
F.run(W, k4mf=True)
json.dump(out, open(f"results/width_{name}_mlp{mlp}.json", "w"), indent=1)
