"""G4: is each loop-length class (= source age) of the memory blockwise low rank in a FROZEN frame?
For every FC source s, at the freeze time m = s + FA (default FA = 3, when it turns old) build one basis per leg mode
(Y-mode, Z-mode; top-R eigenvectors of the energy-weighted leg Gram in layer-m neuron space), project the legs once,
and from then on only transport them exactly (static tensor read through the moving frame). At every later target t
compare the source's D21 contribution
  total D_s(t)                     vs  D_s^R(t)
  loop part L_s = D_s - D_s^Bethe  vs  L_s^R = D_s^R - Bethe(projected legs)   (Bethe = rank-one mean-field of squared legs)
and the Bethe part alone. Output per (age at target, R). usage: python t_loops.py SET MLP [FA]"""
import sys, json
import numpy as np
sys.path.insert(0, "../../fresh-slate/bench")
import bench, fc_hook as F
name, mlp = sys.argv[1], int(sys.argv[2]); FA = int(sys.argv[3]) if len(sys.argv) > 3 else 3
S = bench.load_set(name); W = bench.weights(S, mlp); n = W.shape[1]
RS = [n // 16, n // 8, n // 4, 3 * n // 8]
F.DIAG = True; F.GS.clear()
rows = []

def mf(A):
    A = A.astype(np.float64); return np.outer(A.sum(1), A.sum(0)) / A.sum()

def Dof(Y, Z, w, bethe=False):
    Y = Y.astype(np.float64); Z = Z.astype(np.float64)
    Y2 = Y * Y; YZ = Y * Z
    if bethe: Y2, YZ = mf(Y2), mf(YZ)
    return (Y2 * w[:, None]).T @ Z + 2 * ((YZ * w[:, None]).T @ Y)

def rel(A, B): return float(np.linalg.norm(A - B) / np.linalg.norm(B))

def hook(t, sources, D21):
    G = F.GS[-1]
    for s in sources:
        w = s["w2"].astype(np.float64); SP = s["SP"].astype(np.float32)
        age = t - s["s"]
        if "fz" in s:                      # transport frozen projected legs one layer
            for R in RS:
                s["fz"][R] = (s["fz"][R][0] @ G, s["fz"][R][1] @ G)
        if age == FA:                      # freeze: per-mode bases in the current (layer-t) frame
            Y = (SP @ s["Z"]).astype(np.float64); Z = s["Z"].astype(np.float64)
            ey = (w ** 2) * (Z ** 2).sum(1); ez = (w ** 2) * ((Y ** 2).sum(1) ** 2)
            _, Uy = np.linalg.eigh((Y * ey[:, None]).T @ Y); _, Uz = np.linalg.eigh((Z * ez[:, None]).T @ Z)
            Uy, Uz = Uy[:, ::-1], Uz[:, ::-1]
            s["fz"] = {}
            for R in RS:
                Py = Uy[:, :R] @ Uy[:, :R].T; Pz = Uz[:, :R] @ Uz[:, :R].T
                s["fz"][R] = ((Y @ Py).astype(np.float32), (Z @ Pz).astype(np.float32))   # Y-mode legs, Z-mode legs
        if "fz" in s and age > FA:
            Y = SP @ s["Z"]; Z = s["Z"]
            D = Dof(Y, Z, w); Db = Dof(Y, Z, w, True); L = D - Db
            r = dict(t=t, s=s["s"], age=age, loop_share=float(np.linalg.norm(L) / np.linalg.norm(D)),
                     share_of_D21=float(np.linalg.norm(D) / np.linalg.norm(D21)))
            for R in RS:
                Yr, Zr = s["fz"][R]
                # Tucker: the Y legs that pair with Z use the Y-mode basis, Z legs the Z-mode basis
                DR = Dof(Yr, Zr, w); DbR = Dof(Yr, Zr, w, True)
                r[f"tot_{R}"] = rel(DR, D); r[f"loop_{R}"] = rel(DR - DbR, L); r[f"bethe_{R}"] = rel(DbR, Db)
            rows.append(r)
    return False

F.HOOK = hook
F.run(W, k4mf=True)
json.dump(rows, open(f"results/loops_{name}_mlp{mlp}_fa{FA}.json", "w"), indent=1)
import collections
agg = collections.defaultdict(list)
for r in rows: agg[r["age"]].append(r)
for a in sorted(agg):
    rs = agg[a]; m = lambda k: np.mean([x[k] for x in rs])
    print(f"age {a:2d} (n={len(rs)}): loop share {m('loop_share'):.2f}  " + "  ".join(
        f"R={R}: tot {m(f'tot_{R}'):.3f} loop {m(f'loop_{R}'):.3f}" for R in RS), flush=True)
