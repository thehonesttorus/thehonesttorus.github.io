"""Backend-agnostic (numpy or flopscope.numpy) two-tier K=3 chain, for measuring the billed cost under the grader's
rules. Same algorithm as whest/kprop3c.py (dense young sources + Tucker tier), written without in-place updates.
  from whest.kprop3f import kprop3f_chain; out = kprop3f_chain(W, dict(window=2, k=128), backend="fnp", dtype="float32")
"""
import math
import numpy as np


class Backend:
    def __init__(self, name="numpy", dtype="float64"):
        self.name = name; self.dt = np.float32 if dtype == "float32" else np.float64
        if name == "fnp":
            import flopscope as fs, flopscope.numpy as fnp
            self.xp = fnp; self.fs = fs
            self.cdf = lambda x: fs.stats.norm.cdf(x); self.pdf = lambda x: fs.stats.norm.pdf(x)
        else:
            from scipy.special import ndtr
            self.xp = np; self.fs = None
            self.cdf = ndtr; self.pdf = lambda x: np.exp(-0.5 * x * x) / math.sqrt(2 * math.pi)

    def arr(self, a):
        return self.xp.array(np.asarray(a, dtype=self.dt))

    def sandwich(self, W, C):
        """W C W^T with the symmetric-output discount when billed."""
        if self.name == "fnp":
            Cs = self.fs.as_symmetric(C, symmetry=(0, 1)); return self.xp.einsum("ij,jk,lk->il", W, Cs, W)
        return W @ C @ W.T

    def zero_diag(self, X):
        xp = self.xp; return X - xp.diag(xp.diag(X))


def wick_all(be, m, var):
    """Wick coefficients E[d^k relu^p], k = 0..4, p = 1..4, on the backend (n-vectors; cost negligible)."""
    xp = be.xp; sig = xp.sqrt(var); alpha = m / sig; P = be.cdf(alpha); p_ = be.pdf(alpha)
    T = [P, p_, P - alpha * p_, (alpha ** 2 + 2) * p_, 3 * P - (alpha ** 3 + 3 * alpha) * p_]   # E[X^j 1[X>-alpha]]
    He = [1.0, alpha, alpha ** 2 - 1, alpha ** 3 - 3 * alpha]
    w = {}
    for p in (1, 2, 3, 4):
        for k in range(5):
            if k >= p:
                kk = k - p + 1; f = math.factorial(p)
                if kk == 0: w[(k, p)] = f * (sig * p_ + m * P)
                elif kk == 1: w[(k, p)] = f * P
                else: w[(k, p)] = f * (-1) ** (kk - 2) * sig ** (-(kk - 1)) * He[kk - 2] * p_
            else:
                q = p - k; fall = math.prod(range(q + 1, p + 1))
                w[(k, p)] = fall * sum(math.comb(q, j) * m ** (q - j) * sig ** j * T[j] for j in range(q + 1))
    return w


def slices_from_legs(be, legs):
    A, B, C = legs; xp = be.xp
    S21 = ((A * B) @ C.T + (A * C) @ B.T + (B * C) @ A.T) / 3.0
    S3 = xp.diag(S21); return be.zero_diag(S21), S3


def tucker_slices(be, Q, S):
    xp = be.xp
    M = xp.einsum("ip,pqs,iq->is", Q, S, Q)       # sum_pq Q_ip Q_iq S_pqs
    S21 = M @ Q.T; S3 = xp.sum(M * Q, axis=1); return be.zero_diag(S21), S3


def sym3(be, S):
    xp = be.xp
    return (S + xp.transpose(S, (0, 2, 1)) + xp.transpose(S, (1, 0, 2)) + xp.transpose(S, (1, 2, 0)) + xp.transpose(S, (2, 0, 1)) + xp.transpose(S, (2, 1, 0))) / 6.0


def rotate_core(be, S, R):
    return be.xp.einsum("pqs,ap,bq,cs->abc", S, R, R, R)


def merge_into_tucker(be, Q, S, legs, k, rng, fixed=None):
    xp = be.xp; A, B, C = legs; n = A.shape[0]
    if Q is None:
        G = xp.concatenate([A, B, C], axis=1)
    else:
        nu = xp.sqrt(xp.einsum("pqs,pqs->p", S, S)); G = xp.concatenate([Q * nu[None, :], A, B, C], axis=1)
    Om = be.arr(rng.standard_normal((G.shape[1], k + 8))); Y = G @ Om; Y = G @ (G.T @ Y)
    if fixed is not None:
        F = xp.linalg.qr(fixed)[0]; Y = Y - F @ (F.T @ Y)
        Qn = xp.concatenate([F, xp.linalg.qr(Y)[0][:, :k - F.shape[1]]], axis=1)
    else:
        Qn = xp.linalg.qr(Y)[0][:, :k]
    Sn = None
    if Q is not None:
        P = Qn.T @ Q; Sn = rotate_core(be, S, P)
    a, b, c = Qn.T @ A, Qn.T @ B, Qn.T @ C
    blk = xp.einsum("ar,br,cr->abc", a, b, c)
    Sn = blk if Sn is None else Sn + blk
    return Qn, sym3(be, Sn)


def kprop3f_chain(W, opts=None, backend="numpy", dtype="float64", record=None):
    o = dict(window=2, k=128, c4scale=1.0, seed=0, collective=0, radial=1); o.update(opts or {})
    be = Backend(backend, dtype); xp = be.xp; L, n, n_in = W.shape; rng = np.random.default_rng(o["seed"])
    eye = be.arr(np.eye(n)); ones = be.arr(np.ones(n))
    mu = be.arr(np.zeros(n_in)); C = be.arr(np.eye(n_in)); young = []; Q = S = None; c4 = 0.0; M = be.arr(np.eye(n_in))
    a22 = 3.0 / (n * (n + 2)); out = []
    for l in range(L):
        Wl = be.arr(W[l])
        # ---- linear step
        m = Wl @ mu; Cz = be.sandwich(Wl, C); young = [(tuple(Wl @ Lg for Lg in legs), age) for legs, age in young]
        if Q is not None:
            Q, R = xp.linalg.qr(Wl @ Q); S = rotate_core(be, S, R)
        Mz = be.sandwich(Wl, M)
        # ---- nonlinear step
        var = xp.clip(xp.diag(Cz), 1e-30, None); Soff = be.zero_diag(Cz); w = wick_all(be, m, var)
        K3_21 = xp.zeros((n, n), dtype=be.dt) if backend == "numpy" else be.arr(np.zeros((n, n))); K3_3 = be.arr(np.zeros(n))
        for legs, age in young:
            s21, s3 = slices_from_legs(be, legs); K3_21 = K3_21 + s21; K3_3 = K3_3 + s3
        if Q is not None:
            s21, s3 = tucker_slices(be, Q, S); K3_21 = K3_21 + s21; K3_3 = K3_3 + s3
        Md = xp.diag(Mz); K4_22 = be.zero_diag(c4 * (xp.outer(Md, Md) + 2 * Mz * Mz) / 3.0); K4_4 = c4 * Md ** 2
        pK = {}
        for p in (1, 2, 3, 4):
            pK[(p,)] = w[(0, p)] + K3_3 * w[(3, p)] / 6.0 + K4_4 * w[(4, p)] / 24.0
        def two_index(pi, pj):
            wi1, wj1 = w[(1, pi)], w[(1, pj)]; wi2, wj2 = w[(2, pi)], w[(2, pj)]
            return (Soff * xp.outer(wi1, wj1) + 0.5 * (K3_21.T * xp.outer(wi1, wj2) + K3_21 * xp.outer(wi2, wj1))
                    + 0.5 * Soff * Soff * xp.outer(wi2, wj2) + 0.25 * K4_22 * xp.outer(wi2, wj2))
        p1, p2, p3, p4 = pK[(1,)], pK[(2,)], pK[(3,)], pK[(4,)]
        Choff = be.zero_diag(two_index(1, 1)); Ch = Choff + xp.diag(p2 - p1 ** 2)
        K3h_3 = p3 - 3 * p2 * p1 + 2 * p1 ** 3
        K3h_21 = be.zero_diag(two_index(2, 1) - 2 * p1[:, None] * Choff)
        if o["radial"]:
            K4h_4 = p4 - 4 * p3 * p1 - 3 * p2 ** 2 + 12 * p2 * p1 ** 2 - 6 * p1 ** 4
            K4h_22 = be.zero_diag(two_index(2, 2) - 2 * p1[:, None] * K3h_21.T - 2 * p1[None, :] * K3h_21 - 2 * Choff ** 2 - 4 * xp.outer(p1, p1) * Choff)
            c4 = float(a22 * (xp.sum(K4h_22) + xp.sum(K4h_4))) * o["c4scale"]
        Phi = w[(1, 1)]; w2 = w[(2, 1)]
        young = [(tuple(Lg * Phi[:, None] for Lg in legs), age + 1) for legs, age in young]
        if Q is not None:
            Q, R = xp.linalg.qr(Phi[:, None] * Q); S = rotate_core(be, S, R)
        star = (Phi[:, None] * Soff, 3.0 * eye, (w2[:, None] * Soff * Phi[None, :]).T)
        own21 = be.arr(np.zeros((n, n))); own3 = be.arr(np.zeros(n))
        for legs, age in young:
            s21, s3 = slices_from_legs(be, legs); own21 = own21 + s21; own3 = own3 + s3
        if Q is not None:
            s21, s3 = tucker_slices(be, Q, S); own21 = own21 + s21; own3 = own3 + s3
        s21, s3 = slices_from_legs(be, star); own21 = own21 + s21; own3 = own3 + s3
        R21 = K3h_21 - own21; R3 = K3h_3 - own3
        newborn = tuple(xp.concatenate([F, G], axis=1) for F, G in zip(star, (xp.diag(R3) + 3.0 * R21.T, eye, eye)))
        young.append((newborn, 0))
        keep = []; fixed = None
        if o["collective"]:
            v = p1
            for _ in range(4): v = Ch @ v; v = v / xp.sqrt(xp.sum(v * v))
            fixed = xp.stack([p1, xp.diag(Ch), ones, v], axis=1)
        for legs, age in young:
            if age >= o["window"]: Q, S = merge_into_tucker(be, Q, S, legs, o["k"], rng, fixed)
            else: keep.append((legs, age))
        young = keep; mu, C, M = p1, Ch, eye
        out.append(np.asarray(p1, dtype=np.float64))
        if record is not None: record[l] = dict(c4=c4, young=len(young))
    return np.array(out)
