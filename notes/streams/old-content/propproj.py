"""Transfer-operator test: does transported old content concentrate on the top singular directions of its propagator?

For source s (born at the post-activation of layer s) at layer t (age t - s), the masked propagator is
    M_{s->t} = diag(Phi_t) W_t^T diag(Phi_{t-1}) W_{t-1}^T ... diag(Phi_{s+1}) W_{s+1}^T     (post_s -> gated z_t)
and the source enters D21(t+1) through Y_{s,t} = AD(Phi_t^3 . X_{s,t}),  X_{s,t+1} = W_{t+1}^{(x)3} Y_{s,t}.
We project Y_{s,t} in all three indices on the top-k left singular vectors U_k of M_{s->t} and report
    keep_k = 1 - ||D21(W^{(x)3} Y_hat) - D21(X_{s,t+1})|| / ||D21(X_{s,t+1})||      (per source, by age)
and, for the old pool (ages > w, every source with its own basis), eps_k = || sum error || / ||D21(z_{t+1})||.
Also the participation ratio PR = (sum s^2)^2 / sum s^4 and r90 (rank holding 90 % of sum s^2) of M_{s->t}.
Cost of the carrier at n = 1024: per age band k transported vectors (k/n units) plus a k^3 core.

    python propproj.py ATLAS.npz [--model] [--w 1] [--ks 4,8,16,32,64]
"""
import sys, os, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "experiments"))
from tracker import load_atlas, load_pairs_atlas, model_lay, T3, phi3, d21  # noqa: E402
from oracle_k3 import all_distinct  # noqa: E402


def proj3(Y, U):
    n, k = U.shape
    T = (Y.reshape(n * n, n) @ U).reshape(n, n, k)
    T = np.einsum("ijc,jb->ibc", T, U, optimize=True)
    T = np.einsum("ibc,ia->abc", T, U, optimize=True)          # core (k, k, k)
    T = (T.reshape(k * k, k) @ U.T).reshape(k, k, n)
    T = np.einsum("ijc,bj->ibc", T, U, optimize=True)
    return np.einsum("ibc,ai->abc", T, U, optimize=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("atlas"); ap.add_argument("--model", action="store_true")
    ap.add_argument("--noad", action="store_true", help="published-chain convention: sources transported without the AD mask")
    ap.add_argument("--w", type=int, default=1); ap.add_argument("--ks", default="4,8,16,32,64")
    args = ap.parse_args()
    if "pre_M3" in np.load(args.atlas).files:
        W, lay, N = load_atlas(args.atlas)
    else:
        W, lay, N = load_pairs_atlas(args.atlas); args.model = True
    if args.model:
        lay = model_lay(W, lay, cache=args.atlas.replace('.npz', '_model.npz'))
    L, n, _ = W.shape
    ks = [int(x) for x in args.ks.split(",")]
    print(f"{args.atlas}{' [MODEL chain]' if args.model else ''}: width {n}, depth {L}, N = {N}, old = age > {args.w} at layer t{', no-AD sources' if args.noad else ''}")
    print("rows: t (the layer of X_{s,t}), age, PR and r90 of M_{s->t}, share = ||D21 contrib at t+1|| / ||D21(z_{t+1})||, keep_k per k")
    print("pool rows: eps_k of the old pool (each source projected on its own basis) relative to ||D21(z_{t+1})||, and 'none' = pool share")
    cur, Ms = {}, {}
    agg = {}   # age -> list of keep_k rows
    for l in range(L - 1):
        Phi = lay[l]["Phi"]
        mask = (lambda K: K) if args.noad else all_distinct
        Ol = mask(phi3(lay[l]["K3z"], Phi))
        Bl = lay[l]["K3a"] - Ol
        # sources at layer l are X_{s,l}; build their contribution to D21(l+1) and its projections
        Dn = d21(lay[l + 1]["K3z"]); dn = np.linalg.norm(Dn)
        pool_true = np.zeros((n, n)); pool_hat = {k: np.zeros((n, n)) for k in ks}
        nxt = {}
        for s, X in cur.items():
            age = l - s
            M = Phi[:, None] * Ms[s]                        # diag(Phi_l) W_l^T ... : post_s -> gated z_l
            U, sv, _ = np.linalg.svd(M)
            e = sv ** 2; pr = e.sum() ** 2 / (e ** 2).sum(); r90 = int(np.searchsorted(np.cumsum(e) / e.sum(), 0.9) + 1)
            Y = mask(phi3(X, Phi))
            Xn = T3(Y, W[l + 1]); nxt[s] = Xn
            Dt = d21(Xn); nt = np.linalg.norm(Dt)
            keeps = []
            for k in ks:
                Dh = d21(T3(proj3(Y, U[:, :k]), W[l + 1]))
                keeps.append(1 - np.linalg.norm(Dh - Dt) / nt)
                if age + 1 > args.w:
                    pool_hat[k] += Dh
            if age + 1 > args.w:
                pool_true += Dt
            agg.setdefault(age, []).append([pr, r90] + keeps)
            print(f"t={l:>2} age={age:>2} PR={pr:6.1f} r90={r90:4d} share={nt/dn:5.3f} | " + " ".join(f"k{k}:{v:+.3f}" for k, v in zip(ks, keeps)), flush=True)
        if np.linalg.norm(pool_true) > 0:
            print(f"POOL t+1={l+1:>2} none={np.linalg.norm(pool_true)/dn:5.3f} | "
                  + " ".join(f"k{k}:{np.linalg.norm(pool_hat[k]-pool_true)/dn:5.3f}" for k in ks), flush=True)
        for s in list(Ms):
            Ms[s] = W[l + 1].T @ (Phi[:, None] * Ms[s])
        nxt[l] = T3(Bl, W[l + 1]); Ms[l] = W[l + 1].T.copy()
        cur = nxt
    print("\nby age (mean over sources): PR, r90, keep_k")
    for age in sorted(agg):
        a = np.mean(np.array(agg[age]), 0)
        print(f"age={age:>2} n_src={len(agg[age]):>2} PR={a[0]:6.1f} r90={a[1]:6.1f} | " + " ".join(f"k{k}:{v:+.3f}" for k, v in zip(ks, a[2:])))


if __name__ == "__main__":
    main()
