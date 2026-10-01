"""Cheap carriers of the transported old-source content of the third cumulant, tested against the exact per-source
tracker (tracker.py) at small width.

Target.  With a young window w (sources of age <= w are carried exactly, as dense/young tiers do), the old pool
Old_l = sum_{age > w} X_{s,l} obeys  Old_{l+1} = W^{(x)3} AD(Phi_l^3 . (Old_l + In_l)),  In_l = the source of age w.
A carrier replaces Old_l by a cheap representation; its error is
    eps_l = || D21(Old_hat_l) - D21(Old_l) || / || D21(z_l) ||        (D21 of the whole pre-activation as denominator)
which is the D21 epsilon of the published chain's error law (extra MSE ~ 4.2e-6 eps^2; frontier eps <= 2.2 %).
'Dynamic' carriers feed their own approximation forward (compounding); 'static' ones are re-fitted to the true
Old_l at every layer (an oracle lower bound for that representation class).

Carriers (cost in units per layer at n = 1024, 1 unit = one n^3 matmul = 2 n^3 flops):
  modes-static r   (a) best rank-r family of covariance-response pairs Old ~ sum_r v_r (x) S_r (SVD of the (n, n^2)
                   unfolding), D21 read from the symmetrised family: D21_ab = sum_r (2 v_a S_ab + v_b S_aa)/3
  modes-dyn r      (a) the same family transported (v -> W^T Phi v, S -> W^T Phi S Phi W: one sandwich per mode),
                   the incoming source added and the sum re-truncated to r modes every layer (oracle re-truncation)
  modes-C r        (e) family with FIXED tilt directions: v-space = top-r eigenvectors of the pre-activation C_l
                   (a chain has C), S_r = Old(v_r, ., .)  (static projection, tells whether old modes align with C)
  cp-static R      (b) nonsymmetric CP (hub columns x_r (x) y_r (x) z_r) of Old_l by ALS, D21 = (X*Y) Z^T symmetrised
  cp-dyn R         (b/d) CP factors transported (each leg x -> W^T Phi x), incoming source merged by warm-started ALS
  regress          (c) least-squares fit (per layer, in-sample) of D21(Old_l) on O(n^2) matrices a chain has
"""
import sys, os, argparse, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "experiments"))
from tracker import load_atlas, load_pairs_atlas, old_pool, run_sources, model_lay, T3, phi3, d21  # noqa: E402
from oracle_k3 import all_distinct, sym3  # noqa: E402


def nrm(x):
    return float(np.linalg.norm(x))


# ---------------------------------------------------------------- (a) covariance-response mode families
def modes_fit(K, r):
    """best rank-r unfolding approximation K ~ sum_r v_r (x) S_r (v orthonormal, S_r = K(v_r, ., .))."""
    n = K.shape[0]
    M = K.reshape(n, n * n)
    ev, U = np.linalg.eigh(M @ M.T)
    V = U[:, ::-1][:, :r]
    S = (V.T @ M).reshape(r, n, n)
    S = 0.5 * (S + S.transpose(0, 2, 1))
    return V, S


def modes_dense(V, S):
    n = V.shape[0]
    return sym3((V @ S.reshape(V.shape[1], n * n)).reshape(n, n, n))


def modes_d21(V, S):
    """D21 of sym3(sum_r v_r (x) S_r) without forming the tensor: (2 v_a S_ab + v_b S_aa)/3."""
    diagS = np.einsum("raa->ra", S)
    return (2.0 * np.einsum("ar,rab->ab", V, S) + (diagS.T @ V.T)) / 3.0


# ---------------------------------------------------------------- (b) CP / hub columns by ALS
def khatri_rao(B, C):
    n, R = B.shape
    return (B[:, None, :] * C[None, :, :]).reshape(n * C.shape[0], R)


def cp_als(K, R, iters=30, init=None, seed=0, tol=1e-6):
    n = K.shape[0]
    rng = np.random.default_rng(seed)
    if init is None:
        # HOSVD-style init: leading mode-1 singular vectors, padded with random columns
        A = rng.standard_normal((n, R)); B = rng.standard_normal((n, R)); C = rng.standard_normal((n, R))
    else:
        A, B, C = [x.copy() for x in init]
    K1 = K.reshape(n, n * n)                    # mode 1: i | (j k)
    K2 = K.transpose(1, 0, 2).reshape(n, n * n)  # mode 2: j | (i k)
    K3 = K.transpose(2, 0, 1).reshape(n, n * n)  # mode 3: k | (i j)
    nK = nrm(K); prev = None
    for it in range(iters):
        A = K1 @ khatri_rao(B, C) @ np.linalg.pinv((B.T @ B) * (C.T @ C))
        B = K2 @ khatri_rao(A, C) @ np.linalg.pinv((A.T @ A) * (C.T @ C))
        C = K3 @ khatri_rao(A, B) @ np.linalg.pinv((A.T @ A) * (B.T @ B))
        if it % 5 == 4 or it == iters - 1:
            fit = nrm(K - np.einsum("ir,jr,kr->ijk", A, B, C, optimize=True)) / nK
            if prev is not None and abs(prev - fit) < tol:
                break
            prev = fit
    return A, B, C


def cp_dense(A, B, C):
    return np.einsum("ir,jr,kr->ijk", A, B, C, optimize=True)


def cp_d21_sym(A, B, C):
    """D21 of the symmetrised CP tensor."""
    T = cp_dense(A, B, C)
    return d21(sym3(T))


# ---------------------------------------------------------------- (c) regression features
def slice_transport_features(M, Phi, W):
    """D21(l+1) of W^{(x)3} Phi^3 . (slice tensor with (i,i,k) entries M_ik), its two distinct placements."""
    G = (Phi ** 2)[:, None] * M * Phi[None, :]
    f_a = (W * W).T @ G @ W                       # (i,i,k) -> (a,a,b)
    f_b = (W * (G @ W)).T @ W                     # (i,k,i) -> (a,b,a)/(b,a,a) mix
    return f_a, f_b


def regress_eps(target, feats, denom):
    X = np.stack([f.ravel() for f in feats], 1)
    coef, *_ = np.linalg.lstsq(X, target.ravel(), rcond=None)
    return nrm(target.ravel() - X @ coef) / denom, coef


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("atlas")
    ap.add_argument("--w", type=int, default=1, help="young window: sources of age <= w are exact")
    ap.add_argument("--which", default="modes,regress", help="comma list: modes, modesdyn, modesC, cp, cpdyn, regress")
    ap.add_argument("--ranks", default="1,2,4,8,16,32")
    ap.add_argument("--cpranks", default="32,64,128,256")
    ap.add_argument("--cplayers", default="", help="layers for cp-static (default all)")
    ap.add_argument("--iters", type=int, default=30)
    ap.add_argument("--noad", action="store_true", help="published-chain source convention (no AD mask on transported sources)")
    ap.add_argument("--model", action="store_true", help="run on the noise-free closure chain driven by the atlas (tracker.model_lay)")
    args = ap.parse_args()
    if "pre_M3" in np.load(args.atlas).files:
        W, lay, N = load_atlas(args.atlas)
    else:
        W, lay, N = load_pairs_atlas(args.atlas); args.model = True
    if args.model:
        lay = model_lay(W, lay, cache=args.atlas.replace('.npz', '_model.npz'))
    L, n, _ = W.shape
    which = set(args.which.split(","))
    ranks = [int(x) for x in args.ranks.split(",")]
    cpranks = [int(x) for x in args.cpranks.split(",")]
    cplayers = set(int(x) for x in args.cplayers.split(",")) if args.cplayers else None
    print(f"{args.atlas}{' [MODEL chain]' if args.model else ''}: width {n}, depth {L}, N = {N}, young window w = {args.w}{', no-AD sources' if args.noad else ''}", flush=True)
    global MASK
    MASK = (lambda K: K) if args.noad else all_distinct
    pool = list(old_pool(W, lay, args.w, ad=not args.noad))      # (l, Old_l, In_l), l = 1..L-1
    D21z = {l: d21(lay[l]["K3z"]) for l in range(L)}

    print(f"\nold share: ||D21(Old_l)|| / ||D21(z_l)||")
    print(" ".join(f"{l}:{nrm(d21(O))/nrm(D21z[l]):.3f}" for l, O, _ in pool), flush=True)

    if "modes" in which:
        print(f"\n(a) modes-static: eps (rel. ||D21(z_l)||) of the symmetrised rank-r family; tensor energy residual in []")
        print(f"{'l':>2} | " + " ".join(f"r={r:<11d}" for r in ranks))
        for l, O, _ in pool:
            if l <= args.w:
                continue
            den = nrm(D21z[l]); row = []
            V, S = modes_fit(O, max(ranks))
            for r in ranks:
                Kr = (V[:, :r] @ S[:r].reshape(r, n * n)).reshape(n, n, n)
                e = nrm(modes_d21(V[:, :r], S[:r]) - d21(O)) / den
                res = nrm(Kr - O) / nrm(O)
                row.append(f"{e:5.3f} [{res:4.2f}]")
            print(f"{l:>2} | " + " ".join(row), flush=True)

    if "modesS" in which:
        # (e) structure of the mode matrices S_p of the optimal rank-r family: (i) each S_p truncated to rank q
        # (transport then costs 2 q/n units per mode instead of one sandwich); (ii) each S_p replaced by its least-squares
        # fit on matrices a chain already has at layer l (C, C*C, mu mu^T, var/mu outer products, I)
        rr = [r for r in (4, 8, 16) if r <= n]; qs = [q for q in (4, 8, 16, 32, 64) if q < n]
        print(f"\n(e) modes with structured S_p: eps for r in {rr}; columns: exact S | rank-q S for q in {qs} | S fitted on C-features")
        for l, O, _ in pool:
            if l <= args.w:
                continue
            den = nrm(D21z[l]); C = lay[l]["C"]; mu = lay[l]["mu"]; var = np.diag(C); one = np.ones(n)
            F = [C, C * C, np.outer(mu, mu), np.outer(var, one) + np.outer(one, var), np.outer(mu, one) + np.outer(one, mu), np.eye(n),
                 np.outer(np.sqrt(var), np.sqrt(var))]
            X = np.stack([f.ravel() for f in F], 1)
            out = []
            V, S = modes_fit(O, max(rr))
            for r in rr:
                Vr, Sr = V[:, :r], S[:r]
                row = [nrm(modes_d21(Vr, Sr) - d21(O)) / den]
                for q in qs:
                    Sq = np.empty_like(Sr)
                    for p_ in range(r):
                        ev, U = np.linalg.eigh(Sr[p_]); idx = np.argsort(-np.abs(ev))[:q]
                        Sq[p_] = (U[:, idx] * ev[idx]) @ U[:, idx].T
                    row.append(nrm(modes_d21(Vr, Sq) - d21(O)) / den)
                Sf = np.stack([(X @ np.linalg.lstsq(X, Sr[p_].ravel(), rcond=None)[0]).reshape(n, n) for p_ in range(r)])
                row.append(nrm(modes_d21(Vr, Sf) - d21(O)) / den)
                out.append(f"r={r}: " + " ".join(f"{v:5.3f}" for v in row))
            print(f"{l:>2} | " + " | ".join(out), flush=True)

    if "modesSdyn" in which:
        # dynamic version of the low-rank-S family: state = sym3(sum_p v_p (x) U_p L_p U_p^T); transported exactly
        # (v -> W^T Phi v, U_p -> W^T Phi U_p; plus the AD slice correction when not --noad), the incoming source added,
        # then re-truncated (oracle: unfolding SVD to r modes, each S_p to rank q)
        pairs = [(4, n // 8), (8, n // 8), (8, n // 4), (16, n // 16), (16, n // 8), (16, n // 4), (32, n // 8)]
        print(f"\n(e) modesS-dyn: eps per layer for (r, q) in {pairs}; cost at n = 1024 ~ 2 r q / n units transport+readout")
        state = {pq: None for pq in pairs}
        for l, O, In in pool:
            den = nrm(D21z[l]); row = []
            for (r, q) in pairs:
                if state[(r, q)] is None:
                    Kh = O.copy()
                else:
                    Kprev, Inprev, Phi = state[(r, q)]
                    Kh = T3(MASK(phi3(Kprev + Inprev, Phi)), W[l])
                if nrm(Kh) == 0:
                    row.append(0.0); state[(r, q)] = (Kh, In, lay[l]["Phi"]); continue
                V, S = modes_fit(Kh, r)
                for p_ in range(r):
                    ev, U = np.linalg.eigh(S[p_]); idx = np.argsort(-np.abs(ev))[:q]
                    S[p_] = (U[:, idx] * ev[idx]) @ U[:, idx].T
                Khat = modes_dense(V, S)
                row.append(nrm(d21(Khat) - d21(O)) / den)
                state[(r, q)] = (Khat, In, lay[l]["Phi"])
            print(f"{l:>2} | " + " ".join(f"{v:5.3f}" for v in row), flush=True)

    if "modesU" in which:
        # shared-basis mode family: Old ~ sym3( sum_p v_p (x) U L_p U^T ), V (n x r) and ONE shared U (n x q) = a
        # Tucker (r, q, q) form (HOSVD of the old tensor); transport = r + q vectors (q/n units) + r q x q cores, readout
        # D21 = (2 v_a S_ab + v_b S_aa)/3 with S_p = U L_p U^T.  Static (re-fitted to Old_l) and dynamic (transported V, U
        # and cores, incoming source added, HOSVD re-truncation) versions.
        pairs = [(8, n // 8), (16, n // 8), (8, n // 4), (16, n // 4), (32, n // 4), (16, n // 2), (32, n // 2)]
        def tucker_fit(K, r, q):
            M1 = K.reshape(n, n * n)
            V = np.linalg.eigh(M1 @ M1.T)[1][:, ::-1][:, :r]
            M2 = K.transpose(1, 0, 2).reshape(n, n * n)
            U = np.linalg.eigh(M2 @ M2.T)[1][:, ::-1][:, :q]
            S = (V.T @ M1).reshape(r, n, n)
            core = np.einsum("pij,ia,jb->pab", S, U, U, optimize=True)
            Sq = np.einsum("ia,pab,jb->pij", U, core, U, optimize=True)
            return V, 0.5 * (Sq + Sq.transpose(0, 2, 1))
        print(f"\n(e) modesU (shared basis U, Tucker (r,q,q)): eps for (r, q) in {pairs}; static | dynamic")
        state = {pq: None for pq in pairs}
        for l, O, In in pool:
            den = nrm(D21z[l]); st, dy = [], []
            for (r, q) in pairs:
                if nrm(O) == 0:
                    st.append(0.0); dy.append(0.0); state[(r, q)] = None; continue
                V, S = tucker_fit(O, r, q)
                st.append(nrm(modes_d21(V, S) - d21(O)) / den)
                if state[(r, q)] is None:
                    Kh = O.copy()
                else:
                    Kprev, Inprev, Phi = state[(r, q)]
                    Kh = T3(MASK(phi3(Kprev + Inprev, Phi)), W[l])
                V, S = tucker_fit(Kh, r, q)
                Khat = modes_dense(V, S)
                dy.append(nrm(d21(Khat) - d21(O)) / den)
                state[(r, q)] = (Khat, In, lay[l]["Phi"])
            print(f"{l:>2} | " + " ".join(f"{v:5.3f}" for v in st) + " | " + " ".join(f"{v:5.3f}" for v in dy), flush=True)

    if "modesC" in which:
        print(f"\n(e) modes-C: v-space fixed to the top-r eigenvectors of C_l (and of Phi C Phi of layer l-1 transported);"
              f" S_r = Old(v_r,.,.); eps")
        print(f"{'l':>2} | " + " ".join(f"r={r:<5d}" for r in ranks))
        for l, O, _ in pool:
            if l <= args.w:
                continue
            den = nrm(D21z[l])
            ev, U = np.linalg.eigh(lay[l]["C"]); U = U[:, ::-1]
            row = []
            for r in ranks:
                V = U[:, :r]
                S = (V.T @ O.reshape(n, n * n)).reshape(r, n, n)
                row.append(f"{nrm(modes_d21(V, S) - d21(O)) / den:5.3f}")
            print(f"{l:>2} | " + " ".join(row), flush=True)

    if "modesdyn" in which:
        print(f"\n(a) modes-dyn: transported family, incoming source added, re-truncated to r modes each layer; eps")
        print(f"{'l':>2} | " + " ".join(f"r={r:<5d}" for r in ranks))
        state = {r: None for r in ranks}
        rows = {}
        for l, O, In in pool:
            den = nrm(D21z[l])
            row = []
            for r in ranks:
                if state[r] is None:
                    Kh = O.copy()                         # first layer with a non-empty pool: start exact
                else:
                    Kprev, Inprev, Phi = state[r]
                    Kh = T3(MASK(phi3(Kprev + Inprev, Phi)), W[l])
                V, S = modes_fit(Kh, r)
                Khat = modes_dense(V, S)
                row.append(f"{nrm(d21(Khat) - d21(O)) / den:5.3f}")
                state[r] = (Khat, In, lay[l]["Phi"])
            print(f"{l:>2} | " + " ".join(row), flush=True)

    if "cp" in which:
        print(f"\n(b) cp-static: nonsymmetric CP rank R of Old_l (ALS, {args.iters} it), eps of D21 of the symmetrised CP; tensor residual in []")
        print(f"{'l':>2} | " + " ".join(f"R={R:<11d}" for R in cpranks))
        for l, O, _ in pool:
            if l <= args.w or (cplayers and l not in cplayers):
                continue
            den = nrm(D21z[l]); row = []
            for R in cpranks:
                t0 = time.time()
                A, B, C = cp_als(O, R, iters=args.iters)
                T = sym3(cp_dense(A, B, C))
                row.append(f"{nrm(d21(T) - d21(O)) / den:5.3f} [{nrm(T - O)/nrm(O):4.2f}]")
            print(f"{l:>2} | " + " ".join(row), flush=True)

    if "cpdyn" in which:
        print(f"\n(b/d) cp-dyn: CP factors transported leg by leg, incoming source merged by warm-started ALS ({args.iters} it); eps")
        print(f"{'l':>2} | " + " ".join(f"R={R:<5d}" for R in cpranks))
        state = {R: None for R in cpranks}
        for l, O, In in pool:
            den = nrm(D21z[l]); row = []
            for R in cpranks:
                if state[R] is None:
                    tgt, init = O, None
                else:
                    (A, B, C), Inprev, Phi = state[R]
                    # exact transport of the previous carrier (+ incoming) with the AD projection, then refit;
                    # warm start = the carrier's own transported legs
                    tgt = T3(MASK(phi3(cp_dense(A, B, C) + Inprev, Phi)), W[l])
                    M = W[l].T * Phi[None, :]
                    init = (M @ A, M @ B, M @ C)
                if np.linalg.norm(tgt) == 0:
                    row.append("  -  "); continue
                A, B, C = cp_als(tgt, R, iters=args.iters, init=init)
                T = sym3(cp_dense(A, B, C))
                row.append(f"{nrm(d21(T) - d21(O)) / den:5.3f}")
                state[R] = ((A, B, C), In, lay[l]["Phi"])
            print(f"{l:>2} | " + " ".join(row), flush=True)

    if "ageregress" in which:
        # (c/d) the old pool's D21 regressed on the D21 of the exactly carried young sources (ages 1..w), and the
        # 'geometric tail' (d): the chain carries one more source (age w+1) and the rest (ages >= w+2) is gamma_l x it
        D, _, _ = run_sources(W, lay, ad=not args.noad)
        print(f"\n(c/d) old pool (age > w) D21 vs the young sources' D21 (per-layer in-sample least squares); eps rel ||D21(z_l)||")
        print(f"{'l':>2} | {'none':>5} {'young':>6} {'young+T':>7} | tail: {'none':>5} {'gamma':>6} {'g,y':>6} | gamma")
        for l in range(1, L):
            den = nrm(D21z[l])
            age = {l - s: D[s + 1, l] for s in range(0, l)}
            tgt = sum((age[a] for a in age if a > args.w), np.zeros((n, n)))
            if nrm(tgt) == 0:
                continue
            yf = [age[a] for a in range(1, args.w + 1) if a in age]
            e_y, _ = regress_eps(tgt, yf, den)
            e_yt, _ = regress_eps(tgt, yf + [y.T for y in yf], den)
            tail = sum((age[a] for a in age if a > args.w + 1), np.zeros((n, n)))
            if (args.w + 1) in age and nrm(tail) > 0:
                f1 = age[args.w + 1]
                e_g, cg = regress_eps(tail, [f1], den)
                e_gy, _ = regress_eps(tail, [f1] + yf, den)
                print(f"{l:>2} | {nrm(tgt)/den:5.3f} {e_y:6.3f} {e_yt:7.3f} | tail: {nrm(tail)/den:5.3f} {e_g:6.3f} {e_gy:6.3f} | {cg[0]:+.3f}", flush=True)
            else:
                print(f"{l:>2} | {nrm(tgt)/den:5.3f} {e_y:6.3f} {e_yt:7.3f} |", flush=True)

    if "regress" in which:
        print(f"\n(c) regression of D21(Old_l) on O(n^2) chain objects (in-sample, per layer); eps rel ||D21(z_l)||")
        print(f"{'l':>2} | {'none':>5} {'C-set':>6} {'+sliceT':>7} {'+sandw':>7} {'+young':>7} {'all':>6} | coefficients(all)")
        prev = None
        for l, O, In in pool:
            den = nrm(D21z[l]); tgt = d21(O)
            C = lay[l]["C"]; var = np.diag(C)
            Co = C - np.diag(var)
            young = D21z[l] - tgt
            f_C = [Co, Co * Co, np.outer(var, np.ones(n)) * Co, Co * np.outer(np.ones(n), var)]
            if prev is not None:
                Phi_p, Dold_p, Din_p, D21_p = prev
                st = list(slice_transport_features(Dold_p + Din_p, Phi_p, W[l]))
                st += list(slice_transport_features(D21_p, Phi_p, W[l]))
                Mw = W[l].T * Phi_p[None, :]
                sw = [Mw @ (Dold_p + Din_p) @ Mw.T, Mw @ D21_p @ Mw.T]
            else:
                st, sw = [], []
            yf = [young, young.T]
            e_none = nrm(tgt) / den
            e_C, _ = regress_eps(tgt, f_C, den)
            e_st, _ = regress_eps(tgt, f_C + st, den) if st else (float("nan"), None)
            e_sw, _ = regress_eps(tgt, f_C + st + sw, den) if st else (float("nan"), None)
            e_y, _ = regress_eps(tgt, f_C + yf, den)
            e_all, coef = regress_eps(tgt, f_C + st + sw + yf, den)
            print(f"{l:>2} | {e_none:5.3f} {e_C:6.3f} {e_st:7.3f} {e_sw:7.3f} {e_y:7.3f} {e_all:6.3f} | "
                  + " ".join(f"{c:+.2g}" for c in coef), flush=True)
            prev = (lay[l]["Phi"], tgt, d21(In), D21z[l])


if __name__ == "__main__":
    main()
