# Merge mclatent2.py chunks: the latent third and fourth cumulants at layers l and l + 1, and how the latent fourth
# cumulant of layer l + 1 splits into collective transport, idiosyncratic transmission and births.
#   python latent_an2.py NET "CHUNK_GLOB"
# For each K-vector v (t, t1, s1, a1; see mclatent2.py): covariance S, kappa_3 and kappa_4 as full K-tensors. The gain
# shapes in latent coordinates are kappa_4 = g (S S + S S + S S) and kappa_3 = (g/2) sym(m (x) S) with m = U^T mu (the
# latent coordinates of the pre-activation mean, which d has subtracted). Reported, with g fitted on the truth t1:
#   the non-gain share of kappa_4(t1); the share of kappa_4(t1) and of its non-gain part explained by kappa_4(s1)
#   (everything that is transmitted linearly) and by kappa_4(a1) = kappa_4(t)[A, A, A, A] (the collective part only);
#   the same for kappa_3; and the consistency check kappa_4(a1) = kappa_4(t)[A^4].
import sys, glob, re, numpy as np
net, pat = int(sys.argv[1]), sys.argv[2]
files = sorted(glob.glob(pat))
T = np.load(f"mc2_off{net}_full.npz")


def load(fs):
    tot = None
    for f in fs:
        z = np.load(f)
        if tot is None:
            tot = {k: np.array(z[k], np.float64) for k in z.files if k not in ("src", "K")}
            tot["src"] = z["src"]; tot["K"] = int(z["K"])
        else:
            for k in z.files:
                if k not in ("src", "K") and not k.startswith("A_"):
                    tot[k] += z[k]
    return tot


def tensors(a, l, v, K):
    c = a["c"]; iu = np.triu_indices(K); npair = len(iu[0])
    m = a[f"m1_{v}_{l}"] / c; M2 = a[f"m2_{v}_{l}"] / c
    m3 = a[f"m3_{v}_{l}"] / c; m4 = a[f"m4_{v}_{l}"] / c
    pidx = np.zeros((K, K), int); pidx[iu] = np.arange(npair); pidx[(iu[1], iu[0])] = np.arange(npair)
    M3 = m3[pidx]                                   # K x K x K: E[v_p v_q v_r]
    M4 = m4[pidx][:, :, pidx]                       # K x K x K x K: E[v_p v_q v_r v_s]
    # central moments (the means are small: d was centred on the true mean)
    S = M2 - np.outer(m, m)
    C3 = (M3 - np.einsum("p,qr->pqr", m, M2) - np.einsum("q,pr->pqr", m, M2) - np.einsum("r,pq->pqr", m, M2)
          + 2 * np.einsum("p,q,r->pqr", m, m, m))
    E = lambda *ix: None
    C4 = (M4 - np.einsum("p,qrs->pqrs", m, M3) - np.einsum("q,prs->pqrs", m, M3) - np.einsum("r,pqs->pqrs", m, M3)
          - np.einsum("s,pqr->pqrs", m, M3)
          + np.einsum("p,q,rs->pqrs", m, m, M2) + np.einsum("p,r,qs->pqrs", m, m, M2) + np.einsum("p,s,qr->pqrs", m, m, M2)
          + np.einsum("q,r,ps->pqrs", m, m, M2) + np.einsum("q,s,pr->pqrs", m, m, M2) + np.einsum("r,s,pq->pqrs", m, m, M2)
          - 3 * np.einsum("p,q,r,s->pqrs", m, m, m, m))
    k4 = C4 - (np.einsum("pq,rs->pqrs", S, S) + np.einsum("pr,qs->pqrs", S, S) + np.einsum("ps,qr->pqrs", S, S))
    return S, C3, k4


def gshape4(S):
    return np.einsum("pq,rs->pqrs", S, S) + np.einsum("pr,qs->pqrs", S, S) + np.einsum("ps,qr->pqrs", S, S)


def gshape3(m, S):
    return 0.5 * (np.einsum("p,qr->pqr", m, S) + np.einsum("q,pr->pqr", m, S) + np.einsum("r,pq->pqr", m, S))


en = lambda X: float(np.sum(X * X))
ex = lambda X, P: 1 - en(X - P) / en(X)
full = load(files); K = full["K"]
idx = {int(re.search(r"_(\d+)\.npz$", f).group(1)): f for f in files}
h0 = load([f for i, f in idx.items() if i % 2 == 0]); h1 = load([f for i, f in idx.items() if i % 2 == 1])
z0 = np.load(files[0])
print(f"net {net}: {len(files)} chunks, {int(full['c'])} samples, K = {K}")
for l in [int(x) for x in full["src"]]:
    A = np.asarray(z0[f"A_{l}"], np.float64)
    R = {v: tensors(full, l, v, K) for v in ("t", "t1", "s1", "a1")}
    n0 = tensors(h0, l, "t1", K)[2]; n1 = tensors(h1, l, "t1", K)[2]
    S1, C31, k41 = R["t1"]
    noise = en(n0 - n1) / 4 / en(k41)
    g = float(np.sum(k41 * gshape4(S1))) / en(gshape4(S1))
    ng = lambda v: R[v][2] - g * gshape4(R[v][0])
    Rt1 = ng("t1")
    k4A = np.einsum("pqrs,pa,qb,rc,sd->abcd", R["t"][2], A.T, A.T, A.T, A.T, optimize=True)   # kappa_4(t)[A^4]
    sv = np.linalg.svd(A, compute_uv=False)
    lam_t, lam_t1 = np.linalg.eigvalsh(R["t"][0])[::-1], np.diag(S1)
    print(f"source layer {l:2d} -> {l + 1:2d}: noise {noise:.3f}; gain g {g:.4f}, non-gain share of kappa_4(t1) {en(Rt1) / en(k41):.3f}"
          f"; A singular values {sv[0]:.2f} .. {sv[-1]:.2f} (median {np.median(sv):.2f})"
          f"\n   kappa_4(t1) explained by kappa_4(s1) {ex(k41, R['s1'][2]):+.3f}, by kappa_4(a1) {ex(k41, R['a1'][2]):+.3f}"
          f" | non-gain part explained by non-gain of s1 {ex(Rt1, ng('s1')):+.3f}, of a1 {ex(Rt1, ng('a1')):+.3f}"
          f"; norms relative to it: s1 {en(ng('s1')) / en(Rt1):.3f}, a1 {en(ng('a1')) / en(Rt1):.3f}"
          f"\n   check kappa_4(a1) = kappa_4(t)[A^4]: rel diff {np.sqrt(en(R['a1'][2] - k4A) / en(R['a1'][2])):.1e};"
          f" non-gain share of kappa_4(t) at layer l {en(R['t'][2] - float(np.sum(R['t'][2] * gshape4(R['t'][0]))) / en(gshape4(R['t'][0])) * gshape4(R['t'][0])) / en(R['t'][2]):.3f}", flush=True)
    # third cumulant, with the mean-coupled gain shape in latent coordinates
    U1 = None
    mu1 = np.asarray(T["mu"][l + 1], np.float64)
    ev, V = np.linalg.eigh(np.asarray(T["cov"][l + 1], np.float64)); U1 = V[:, ::-1][:, :K]
    m1 = U1.T @ mu1
    g3s = gshape3(m1, S1); g3 = float(np.sum(C31 * g3s)) / en(g3s); R3 = C31 - g3 * g3s
    print(f"   kappa_3(t1): gain g {g3:.4f} explains {ex(C31, g3 * g3s):+.3f}; explained by kappa_3(s1) {ex(C31, R['s1'][1]):+.3f},"
          f" by kappa_3(a1) {ex(C31, R['a1'][1]):+.3f}; non-gain part by s1 {ex(R3, R['s1'][1] - g3 * gshape3(m1, R['s1'][0])):+.3f}", flush=True)
