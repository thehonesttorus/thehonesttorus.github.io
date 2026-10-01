"""Layer-by-layer validation of the ReLU step of chain.py against a MC atlas (moment_atlas_np.py --k3 --k4).

For every layer l the step is fed the atlas's TRUE pre-activation cumulants of z_l (mu, C, kappa3, kappa4_{aabc})
and its outputs are compared with the atlas's true post-activation cumulants of a_l and with the next layer's true
pre-activation objects.  Relative rms errors (||pred - true|| / ||true||):
  mu, var, Coff        mean, variance, off-diagonal covariance of a_l
  D3a, D21a            (3,) and (2,1) slices of kappa3(a_l)
  D21+ (mode)          D21(z_{l+1}) from the transported full kappa3(a_l) built by the step (true slices replaced by
                       predicted ones) with the all-distinct rule 'mode'
  X+                   kappa4(z_{l+1})_{aabc} from the memoryless kappa4 closure vs the atlas (all entries / (2,1,1) part)
With --pair B.npz (an independent atlas of the same MLP) the noise of each target is printed (rel(A,B)/sqrt 2).

    python validate.py atlas.npz [--pair atlasB.npz] [--order 2]
"""
import sys, argparse
import numpy as np
import chain as ch
import oracle_k3 as ok


def rel(a, b):
    return float(np.sqrt(np.sum((a - b) ** 2) / np.sum(b ** 2)))


def offd(M):
    M = M.copy(); np.fill_diagonal(M, 0.0); return M


def truths(z, l, L):
    mu_a, C_a, K3a = ch.atlas_post(z, l)
    idx = np.arange(len(mu_a))
    D3a = K3a[idx, idx, idx]; D21a = offd(K3a[idx, idx, :])
    out = dict(mu=mu_a, var=np.diag(C_a).copy(), Coff=offd(C_a), D3a=D3a, D21a=D21a, K3a=K3a)
    if l + 1 < L:
        st1 = ch.atlas_state(z, l + 1)
        out["D21n"] = offd(st1.K3[idx, idx, :]); out["Xn"] = st1.X
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("atlas"); ap.add_argument("--pair"); ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--layers", type=str, default=None)
    a = ap.parse_args()
    z = np.load(a.atlas)
    zb = np.load(a.pair) if a.pair else None
    W = z["weights"].astype(np.float64)
    L, n, _ = W.shape
    layers = range(L) if a.layers is None else [int(x) for x in a.layers.split(",")]
    print(f"{a.atlas}: width {n} depth {L} N {int(z['n_samples'])} order {a.order}")
    hdr = f"{'l':>2} {'mu':>7} {'var':>7} {'Coff':>7} {'D3a':>7} {'D21a':>7} | D21+: {'none':>6} {'wick':>6} {'clos':>6} {'true':>6} | {'X+':>6} {'X211+':>6} {'u.C':>6}"
    if zb is not None:
        hdr += " | noise: mu var Coff D3a D21a D21+"
    print(hdr)
    for l in layers:
        st = ch.atlas_state(z, l)
        T = truths(z, l, L)
        mu_a, C_a, K3a, k4s, info = ch.relu_step(st, "closure", order=a.order, want_k4=(l + 1 < L))
        idx = np.arange(n)
        e = [rel(mu_a, T["mu"]), rel(np.diag(C_a), T["var"]), rel(offd(C_a), T["Coff"]),
             rel(K3a[idx, idx, idx], T["D3a"]), rel(offd(K3a[idx, idx, :]), T["D21a"])]
        s = f"{l:>2} " + " ".join(f"{v:7.4f}" for v in e)
        if l + 1 < L:
            Wn = W[l + 1]
            res = []
            for mode in ("none", "wick", "closure"):
                _, _, K3m, _, _ = ch.relu_step(st, mode, order=a.order, want_k4=False)
                res.append(rel(offd(ok.transport_d21(K3m, Wn)), T["D21n"]))
            res.append(rel(offd(ok.transport_d21(T["K3a"], Wn)), T["D21n"]))
            Xp = ch.transport_k4_slices(*k4s, Wn)
            Xt = T["Xn"]
            r211 = rel(ok.all_distinct(Xp), ok.all_distinct(Xt))
            Co = offd(ch.atlas_state(z, l + 1, with_k3=False, with_k4=False).C)
            K211 = ok.all_distinct(Xt)
            u = np.einsum("ijk,jk->i", K211, Co) / float(np.sum(Co * Co))
            ruc = rel(ok.all_distinct(np.einsum("i,jk->ijk", u, Co)), K211)
            s += " | " + " ".join(f"{v:6.3f}" for v in res) + f" | {rel(Xp, Xt):6.3f} {r211:6.3f} {ruc:6.3f}"
        if zb is not None:
            Tb = truths(zb, l, L)
            nz = [rel(T[k], Tb[k]) / np.sqrt(2) for k in ("mu", "var", "Coff", "D3a", "D21a")]
            if l + 1 < L:
                nz.append(rel(T["D21n"], Tb["D21n"]) / np.sqrt(2))
            s += " | " + " ".join(f"{v:6.4f}" for v in nz)
        print(s, flush=True)


if __name__ == "__main__":
    main()
