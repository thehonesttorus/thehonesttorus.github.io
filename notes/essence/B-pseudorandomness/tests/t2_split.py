"""Test T2 (B-pseudorandomness): inside FC (region stream, slices=2, k4mf=True) on bench w1024_d16, split each source's
D21 contribution at target layer t into its O(n)-average over the last fresh layer (the part constant in the
repeated index a: the 'trace / trivial-isotypic' channel, JMR's mean term) and the Wishart fluctuation (JMR's
lambda term). Prediction (lambda_eff ~ 1): the fluctuation share does not decay with age.
Also: the Gram (cosines) of the per-age contributions (orthogonality across ages, F8.3).
usage: python t2_split.py MLP"""
import sys, json
import numpy as np
sys.path.insert(0, '/home/user/thehonesttorus.github.io/notes/fresh-slate/breakthrough/region')
sys.path.insert(0, '/home/user/thehonesttorus.github.io/notes/fresh-slate/bench')
import bench
import fc_hooked as F

TARGETS = [4, 8, 12, 15]
out = {}

def contrib(src):
    w = src["w2"].astype(np.float64)
    Zf = src["Z"].astype(np.float64)
    Y = src["Yp"].astype(np.float64) if "Yp" in src else src["SP"].astype(np.float64) @ Zf
    D = ((Y * Y) * w[:, None]).T @ Zf + 2 * (((Y * Zf) * w[:, None]).T @ Y)
    if "Delta" in src:
        T = src["Delta"].astype(np.float64) @ Zf
        D += (Zf * Zf).T @ T + 2 * ((Zf * T).T @ Zf)
    return D

def hook(t, sources):
    if t not in TARGETS:
        return
    rows = []; Ds = []
    for src in sources:
        age = t - src["s"]
        D = contrib(src)
        Dm = np.broadcast_to(D.mean(0, keepdims=True), D.shape)      # constant in a
        e = (D * D).sum(); em = (Dm * Dm).sum()
        rows.append(dict(age=int(age), energy=float(e), mean_share=float(em / e)))
        Ds.append(D)
    tot = sum(Ds)
    old = sum(D for D, r in zip(Ds, rows) if r["age"] > 2)
    G = np.array([[float((a * b).sum()) for b in Ds] for a in Ds])
    d = np.sqrt(np.diag(G)); cos = G / np.outer(d, d)
    off = cos[~np.eye(len(Ds), dtype=bool)]
    Om = np.broadcast_to(old.mean(0, keepdims=True), old.shape)
    out[t] = dict(rows=rows, total_energy=float((tot * tot).sum()),
                  old_share_of_total=float((old * old).sum() / (tot * tot).sum()),
                  old_mean_share=float((Om * Om).sum() / (old * old).sum()),
                  sum_energies_over_total=float(sum(r["energy"] for r in rows) / (tot * tot).sum()),
                  cos_rms_offdiag=float(np.sqrt((off ** 2).mean())), cos_max_abs=float(np.abs(off).max()))
    print(t, json.dumps({k: v for k, v in out[t].items() if k != "rows"}), flush=True)
    for r in rows:
        print("   age %2d  energy/total %.3f  mean(a-const) share %.3f" % (r["age"], r["energy"] / out[t]["total_energy"], r["mean_share"]), flush=True)

if __name__ == "__main__":
  F.HOOK = hook
  mlp = int(sys.argv[1])
  S = bench.load_set('w1024_d16')
  W = bench.weights(S, mlp)
  F.run(W, slices=2, k4mf=True)
  json.dump(out, open(f't2_split_mlp{mlp}.json', 'w'), indent=1)
