"""Per-source tracker of the third cumulant at small width (numpy only, dense n^3 tensors).

Input: an atlas npz from notes/experiments/moment_atlas_np.py --k3 (full raw third moments pre_M3, post_M3).

Exact sample-level decomposition.  For every layer l let
    O_l = AD( Phi_l (x) Phi_l (x) Phi_l . kappa3(z_l) )          (the "old content": all-distinct entries, gate-weighted)
    B_l = kappa3(a_l) - O_l                                       (the "birth" of layer l: everything else, incl. all slices)
with Phi_l = P(z_l > 0) (atlas gate_p) and AD = zero every entry with a repeated index.  Since kappa3 of a linear map is
the map applied three times (exactly, also for sample moments), kappa3(z_{l+1}) = W_{l+1}^{(x)3} (O_l + B_l), so
    kappa3(z_l) = sum_{s < l} X_{s,l},   X_{s,s+1} = W_{s+1}^{(x)3} B_s,   X_{s,l+1} = W_{l+1}^{(x)3} AD(Phi_l^3 . X_{s,l})
plus the input's sample kappa3(z_0) (pure Monte Carlo noise, true value 0), tracked as source -1.  The age of X_{s,l} is
l - s (age 1 = the newborn of layer l-1 after one transport).  Every source's contribution to the (2,1) slice
D21(z_l)_ab = kappa3(z_a, z_a, z_b) is therefore known exactly; the sum over sources must reproduce the atlas's
kappa3(z_l) (checked: 'valid' column, relative Frobenius error of the whole tensor and of D21).

    python tracker.py ATLAS.npz [--out results/tracker_X.txt] [--save d21.npz]

Printed per layer l (target interface D21(z_l), l = 1..L-1):
  valid      ||sum_s X_{s,l} - kappa3(z_l)|| / ||kappa3(z_l)||  and the same for D21
  |src-1|    D21 share of the input-noise source
  age a      ||D21 contribution of the source of age a|| / ||D21(z_l)||   (a = 1..6, then 'old7+' = sum of ages >= 7)
  old>w      ||D21 of all sources of age > w|| / ||D21(z_l)||  for w = 1, 4
  spec       energy fraction of the old (age > 1) tensor captured by the top r modes of its (n, n^2) unfolding,
             i.e. by r "covariance-response" pairs v_r (x) S_r with S_r symmetric
"""
import sys, os, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "experiments"))
from oracle_k3 import central3, all_distinct, slices, hermite_model, residual_basis, CLOSURE_COEF  # noqa: E402


def T3(K, W):
    """W^{(x)3} K: (W^T)_{ai} (W^T)_{bj} (W^T)_{ck} K_ijk."""
    n = K.shape[0]
    T = (K.reshape(n * n, n) @ W).reshape(n, n, n)           # k -> c
    T = np.einsum("ijc,jb->ibc", T, W, optimize=True)         # j -> b
    return np.einsum("ibc,ia->abc", T, W, optimize=True)      # i -> a


def phi3(K, Phi):
    return K * Phi[:, None, None] * Phi[None, :, None] * Phi[None, None, :]


def d21(K):
    return slices(K)[1]


def load_atlas(path):
    z = np.load(path)
    W = z["weights"].astype(np.float64)
    L, n, _ = W.shape
    pre_s, post_s = z["pre_s"], z["post_s"]
    lay = []
    for l in range(L):
        mu = pre_s[0, l]; var = pre_s[1, l] - mu ** 2
        M11 = z["pre_M11"][l].astype(np.float64)
        C = M11 - np.outer(mu, mu)
        K3z = central3(z["pre_M3"][l], M11, mu)
        mua = post_s[0, l]
        K3a = central3(z["post_M3"][l], z["post_M11"][l].astype(np.float64), mua)
        Ca = z["post_M11"][l].astype(np.float64) - np.outer(mua, mua)
        lay.append(dict(mu=mu, var=var, C=C, Ca=Ca, K3z=K3z, K3a=K3a, Phi=z["gate_p"][l].astype(np.float64)))
    return W, lay, int(z["n_samples"])


def merge_atlases(paths, out):
    """average of atlases of the same MLP = one atlas with the summed sample count (all fields are sample means)."""
    zs = [np.load(p) for p in paths]
    assert all(np.array_equal(zs[0]["weights"], q["weights"]) for q in zs[1:])
    Ns = np.array([float(q["n_samples"]) for q in zs]); w = Ns / Ns.sum()
    res = {}
    for k in zs[0].files:
        if k in ("weights", "name", "mlp_seed", "sample_seed"):
            res[k] = zs[0][k]
        elif k == "n_samples":
            res[k] = int(Ns.sum())
        else:
            res[k] = sum(wi * q[k].astype(np.float64) for wi, q in zip(w, zs))
    np.savez(out, **res)


def run_sources(W, lay, keep_tensors=False):
    """exact per-source recursion. Returns
       D[s+1, l] = D21 of X_{s,l} (s = -1..L-2; row 0 = input-noise source), shape (L, L, n, n)
       X[l] = dict s -> X_{s,l} if keep_tensors (only the current layer is kept otherwise)
       valid[l] = (tensor rel err, D21 rel err) of sum_s X_{s,l} vs atlas kappa3(z_l)"""
    L, n, _ = W.shape
    D = np.zeros((L, L, n, n))
    cur = {-1: lay[0]["K3z"].copy()}
    D[0, 0] = d21(cur[-1])
    valid = [(0.0, 0.0)]
    hist = [dict(cur)] if keep_tensors else None
    for l in range(L - 1):
        Phi = lay[l]["Phi"]
        Ol = all_distinct(phi3(lay[l]["K3z"], Phi))
        Bl = lay[l]["K3a"] - Ol
        nxt = {s: T3(all_distinct(phi3(X, Phi)), W[l + 1]) for s, X in cur.items()}
        nxt[l] = T3(Bl, W[l + 1])
        cur = nxt
        tot = sum(cur.values())
        Kt = lay[l + 1]["K3z"]
        valid.append((float(np.linalg.norm(tot - Kt) / np.linalg.norm(Kt)),
                      float(np.linalg.norm(d21(tot) - d21(Kt)) / np.linalg.norm(d21(Kt)))))
        for s, X in cur.items():
            D[s + 1, l + 1] = d21(X)
        if keep_tensors:
            hist.append(dict(cur))
    return D, valid, hist


def model_birth(o, K3z_chain, K3a_atlas):
    """noise-free birth of the 'model' chain: all-distinct part = first-order closure without the old-content term
    (Gaussian rho^2 + rho^3 Wick, the D21(z)-hyperedge diagrams B1, B2 and the D3 diagram B3 with the leg-partition
    coefficients, all built from the chain's own kappa3(z)); slice entries = the atlas's slices of kappa3(a) minus the
    slices of the chain's Phi^3 kappa3(z), so the model's post-activation slices equal the atlas's."""
    oo = dict(o); oo["K3z"] = K3z_chain
    basis = residual_basis(oo)                      # B0..B4 (no K22 / K211 keys)
    ad = hermite_model(o["C"], o["mu"], o["var"], 4) + sum(CLOSURE_COEF[b] * basis[b] for b in (1, 2, 3, 4))
    old_full = phi3(K3z_chain, o["Phi"])
    sl = (K3a_atlas - all_distinct(K3a_atlas)) - (old_full - all_distinct(old_full))
    return all_distinct(ad) + sl


def model_states(W, lay):
    """kappa3(z_l) of the noise-free model chain (list over l) and its births (list over l)."""
    L, n, _ = W.shape
    K = np.zeros((n, n, n)); Ks, Bs = [K], []
    for l in range(L - 1):
        o = dict(lay[l])
        B = model_birth(o, K, lay[l]["K3a"])
        Bs.append(B)
        K = T3(all_distinct(phi3(K, o["Phi"])) + B, W[l + 1])
        Ks.append(K)
    return Ks, Bs


def model_lay(W, lay):
    """a copy of the atlas layer list whose K3z / K3a are the model chain's (so old_pool runs on the model)."""
    Ks, Bs = model_states(W, lay)
    out = []
    for l in range(len(lay)):
        d = dict(lay[l]); d["K3z"] = Ks[l]
        if l < len(Bs):
            d["K3a"] = all_distinct(phi3(Ks[l], d["Phi"])) + Bs[l]
        out.append(d)
    return out


def old_pool(W, lay, w, include_noise=False):
    """generator over layers l = 1..L-1 of (Old_l, In_l): Old_l = sum of sources of age > w at layer l (true tensors),
    In_l = the source of age exactly w at layer l (it joins the old pool at l+1).  Recursion used by the carriers:
        Old_{l+1} = W_{l+1}^{(x)3} AD( Phi_l^3 . (Old_l + In_l) )."""
    L, n, _ = W.shape
    cur = {-1: lay[0]["K3z"].copy()} if include_noise else {}
    for l in range(L - 1):
        Phi = lay[l]["Phi"]
        Ol = all_distinct(phi3(lay[l]["K3z"], Phi))
        Bl = lay[l]["K3a"] - Ol
        cur = {s: T3(all_distinct(phi3(X, Phi)), W[l + 1]) for s, X in cur.items()}
        cur[l] = T3(Bl, W[l + 1])
        lp = l + 1
        old = sum((X for s, X in cur.items() if lp - s > w), np.zeros((n, n, n)))
        inc = sum((X for s, X in cur.items() if lp - s == w), np.zeros((n, n, n)))
        yield lp, old, inc


def unfold_energy(K, ranks):
    n = K.shape[0]
    M = K.reshape(n, n * n)
    ev = np.linalg.eigvalsh(M @ M.T)[::-1]
    ev = np.maximum(ev, 0); e = np.cumsum(ev) / ev.sum()
    return [float(e[r - 1]) for r in ranks]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("atlas")
    ap.add_argument("--save", default=None)
    ap.add_argument("--model", action="store_true", help="track the noise-free closure chain driven by the atlas instead")
    args = ap.parse_args()
    W, lay, N = load_atlas(args.atlas)
    L, n, _ = W.shape
    print(f"{args.atlas}: width {n}, depth {L}, N = {N}")
    if args.model:
        mlay = model_lay(W, lay)
        print("MODEL chain (first-order closure births, atlas slices): rel. error of its D21(z_l) vs the atlas's, per layer")
        print(" ".join(f"{l}:{np.linalg.norm(d21(mlay[l]['K3z']) - d21(lay[l]['K3z'])) / np.linalg.norm(d21(lay[l]['K3z'])):.3f}"
                       for l in range(1, L)), flush=True)
        lay = mlay
    D, valid, hist = run_sources(W, lay, keep_tensors=False)
    ranks = (1, 2, 4, 8, 16, 32, 64)
    print(f"{'l':>2} {'valid K3':>8} {'D21':>8} | {'src-1':>6} | " + " ".join(f"age{a:<2d}" for a in range(1, 7))
          + f" {'old7+':>6} | {'old>1':>6} {'old>4':>6} | spec(old>1) r=" + ",".join(str(r) for r in ranks))
    for l in range(1, L):
        Dt = d21(lay[l]["K3z"]); nt = np.linalg.norm(Dt)
        ages = {}
        for s in range(-1, l):
            ages[l - s] = D[s + 1, l]
        noise = np.linalg.norm(D[0, l]) / nt
        share = [np.linalg.norm(ages[a]) / nt if (a in ages and l - a >= 0) else float("nan") for a in range(1, 7)]
        old7 = sum((D[s + 1, l] for s in range(0, l) if l - s >= 7), np.zeros((n, n)))
        o1 = sum((D[s + 1, l] for s in range(0, l) if l - s > 1), np.zeros((n, n)))
        o4 = sum((D[s + 1, l] for s in range(0, l) if l - s > 4), np.zeros((n, n)))
        print(f"{l:>2} {valid[l][0]:8.1e} {valid[l][1]:8.1e} | {noise:6.3f} | " + " ".join(f"{v:5.3f}" for v in share)
              + f" {np.linalg.norm(old7)/nt:6.3f} | {np.linalg.norm(o1)/nt:6.3f} {np.linalg.norm(o4)/nt:6.3f} |", end="", flush=True)
        print("", flush=True)
    # spectra of the old pool (age > 1) and of single sources by age, from the recursion
    print("\nunfolding spectra (energy captured by top r modes):  old = age>1 pool, inc = the age-1 source, ages 2/4/8 single sources")
    print(f"{'l':>2} {'what':>5} | " + " ".join(f"r={r:<3d}" for r in ranks))
    gen = old_pool(W, lay, 1)
    for l, old, inc in gen:
        rows = [("old", old), ("inc", inc)]
        for nm, K in rows:
            if np.linalg.norm(K) == 0:
                continue
            print(f"{l:>2} {nm:>5} | " + " ".join(f"{v:5.3f}" for v in unfold_energy(K, ranks)), flush=True)
    if args.save:
        np.savez_compressed(args.save, D=D.astype(np.float32), valid=np.array(valid))


if __name__ == "__main__":
    main()
