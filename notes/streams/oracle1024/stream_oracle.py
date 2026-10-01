"""Streaming oracle for the K = 3 closure ladder at the real shape (width 1024), without n^3 tensors.

Every column of notes/experiments/oracle_k3.py is a relative error of the next layer's (2,1) slice
D21(l+1)_ab = kappa3(z_a, z_a, z_b) = sum_ijk W_ia W_ja W_kb kappa3(a_l)_ijk.  Every model of the
all-distinct part of kappa3(a_l) used there is a sum of
  (a) "triangle-form" tensors  X_ijk = f_i g_j h_k A_ij B_jk D_ik  (A, B, D n x n or absent), built
      from n x n statistics (C, D21(z), K22, gates, Gaussian weights, regenerated u_i C_jk), and
  (b) two genuinely n^3 objects, Phi^3 kappa3(z) and the (2,1,1) slice of kappa4(z), which only ever
      enter through their transport and are therefore accumulated per sample as contractions:
        Phi^3 kappa3(z) :  E[x_a^2 x_b],              x = (Phi o y) @ W_{l+1},  y = z - mu
        (2,1,1) slice   :  E[s_a x_a x_b], E[x_a^2 s_b], s = (w2 o y^2) @ W_{l+1}
      minus their Gaussian parts (triangle forms) and minus their coincident-index parts (n x n slices).
Transport of a triangle form costs O(n^3) when one of A, B, D is absent (all closure diagrams) and
O(n^4) (one GEMM per output column) for the Gaussian rho^3 / rho^4 triangles.  The all-distinct
projection is exact: T(all_distinct X) = T(X) - [S_A + S_B + S_C - 2 S_diag], where S_* are the
transports of the three two-index restrictions X_iik, X_kii, X_iki and the diagonal X_iii.

Subcommands
  stats   two-pass streaming statistics of one MLP / sample seed (same samples in both passes, so every
          central moment is the exact sample central moment, as in moment_atlas_np + oracle_k3)
            python stream_oracle.py stats --seed 770000 --width 1024 --n-samples 2000000 --sample-seed 1 --out S.npz
  verify  check every streamed term and ladder column against the dense oracle on an atlas built with
          moment_atlas_np.py --k3 --k4 from the same MLP seed and sample seed
            python stream_oracle.py verify ATLAS.npz STATS.npz
  ladder  single-replica ladder (python stream_oracle.py ladder S.npz) or pair ladder with noise floor
          (python stream_oracle.py ladder A.npz B.npz) for two replicas of the same MLP

Conventions (identical to oracle_k3): Phi = empirical gate P(z > 0); w2, w3, w5 Gaussian derivative
weights from (mu, var); the Hermite K3z term uses the empirical Phi (oracle_k3 uses the Gaussian
Phi(alpha) there only; the difference is reported by `verify`).
"""
import argparse, os, sys, time
import numpy as np
from math import erf, sqrt, pi, factorial

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "experiments"))
CLOSURE_COEF = [1.0, 3.0, 3.0, 1.0, 1.0, 1.5, 1.5]          # oracle_k3.CLOSURE_COEF (B0..B6)


# ----------------------------------------------------------------------------------------- weights
def mlp_weights(seed, width, depth):
    rng = np.random.default_rng(seed)
    return (rng.standard_normal((depth, width, width)) * np.sqrt(2.0 / width)).astype(np.float32)


# ----------------------------------------------------------------------------------------- streaming
def stream_stats(W, n_samples, chunk, sample_seed, log=None):
    """pass 1: mu, var, gates, C of z_l and mu of a_l; pass 2 (same samples): central statistics and the
    per-sample contractions that carry the n^3 objects into D21(l+1)."""
    L, n, _ = W.shape
    f64 = np.float64
    S1 = np.zeros((L, n)); S2 = np.zeros((L, n)); G = np.zeros((L, n)); SA = np.zeros((L, n)); M11 = np.zeros((L, n, n))
    t0 = time.time()

    def batches():
        rng = np.random.default_rng(sample_seed)
        done = 0
        while done < n_samples:
            m = min(chunk, n_samples - done)
            yield rng.standard_normal((m, n), dtype=np.float32)
            done += m

    for x in batches():
        a_prev = x
        for l in range(L):
            z = a_prev @ W[l]
            a = np.maximum(z, 0.0)
            S1[l] += z.sum(0, dtype=f64); S2[l] += (z * z).sum(0, dtype=f64); G[l] += (z > 0).sum(0)
            SA[l] += a.sum(0, dtype=f64); M11[l] += z.T @ z
            a_prev = a
    N = float(n_samples)
    mu = S1 / N; var = S2 / N - mu ** 2; Phi = G / N; mua = SA / N
    C = M11 / N - np.einsum("li,lj->lij", mu, mu)
    sig = np.sqrt(var); alpha = mu / sig
    w2 = np.exp(-0.5 * alpha ** 2) / np.sqrt(2 * pi) / sig
    Co = C.copy()
    for l in range(L):
        np.fill_diagonal(Co[l], 0.0)
    if log:
        log(f"pass 1 done {time.time()-t0:.0f}s")
    mu32, mua32 = mu.astype(np.float32), mua.astype(np.float32)
    Phi32, w232, Co32 = Phi.astype(np.float32), w2.astype(np.float32), Co.astype(np.float32)
    D3z = np.zeros((L, n)); M4z = np.zeros((L, n)); D21z = np.zeros((L, n, n)); M22 = np.zeros((L, n, n)); M31 = np.zeros((L, n, n))
    D3a = np.zeros((L, n)); D21a = np.zeros((L, n, n))
    TB0 = np.zeros((L - 1, n, n)); TS1 = np.zeros((L - 1, n, n)); TS2 = np.zeros((L - 1, n, n)); Uq = np.zeros((L - 1, n))
    done = 0
    for x in batches():
        a_prev = x
        for l in range(L):
            z = a_prev @ W[l]
            a = np.maximum(z, 0.0)
            y = z - mu32[l]; ac = a - mua32[l]
            y2 = y * y; y3 = y2 * y
            D3z[l] += y3.sum(0, dtype=f64); M4z[l] += (y2 * y2).sum(0, dtype=f64)
            D21z[l] += y2.T @ y; M22[l] += y2.T @ y2; M31[l] += y3.T @ y
            ac2 = ac * ac
            D3a[l] += (ac2 * ac).sum(0, dtype=f64); D21a[l] += ac2.T @ ac
            if l < L - 1:
                Wn = W[l + 1]
                xx = (y * Phi32[l]) @ Wn; s = (y2 * w232[l]) @ Wn
                x2 = xx * xx
                TB0[l] += x2.T @ xx; TS1[l] += (s * xx).T @ xx; TS2[l] += x2.T @ s
                q = ((y @ Co32[l]) * y).sum(1)
                Uq[l] += y2.T @ q
            a_prev = a
        done += x.shape[0]
        if log and (done // chunk) % 50 == 0:
            log(f"pass 2 {done}/{n_samples} {time.time()-t0:.0f}s")
    D3z /= N; M4z /= N; D21z /= N; M22 /= N; M31 /= N; D3a /= N; D21a /= N; TB0 /= N; TS1 /= N; TS2 /= N; Uq /= N
    G22 = M22 - np.einsum("li,lj->lij", var, var) - 2 * C ** 2          # E[y_i^2 y_j^2] - var var - 2 C^2 (diag: kappa4 marginal)
    F211 = M31 - 3 * var[:, :, None] * C                               # F_iik = E[y_i^3 y_k] - 3 var_i C_ik
    return dict(n_samples=n_samples, sample_seed=sample_seed, mu=mu, var=var, Phi=Phi, mua=mua, C=C,
                D3z=D3z, M4z=M4z, D21z=D21z, G22=G22, F211=F211, D3a=D3a, D21a=D21a, TB0=TB0, TS1=TS1, TS2=TS2, Uq=Uq)


# ----------------------------------------------------------------------------------------- triangle forms
class Tri:
    """X_ijk = coef f_i g_j h_k A_ij B_jk D_ik ; vecs[s] on slot s, mats[(s,t)] on slot pair s<t (None = ones)."""
    def __init__(self, vecs, mats, coef=1.0):
        self.v = list(vecs); self.m = {k: v for k, v in mats.items() if v is not None}; self.c = coef

    def permute(self, p):
        inv = np.argsort(p)
        v = [None] * 3
        for s in range(3):
            v[inv[s]] = self.v[s]
        m = {}
        for (s, t), M in self.m.items():
            u, w = inv[s], inv[t]
            m[(min(u, w), max(u, w))] = M if u < w else M.T
        return Tri(v, m, self.c)

    def dense(self):
        n = len(next(x for x in self.v if x is not None))
        one = np.ones(n)
        f, g, h = [x if x is not None else one for x in self.v]
        T = self.c * np.einsum("i,j,k->ijk", f, g, h)
        for (s, t), M in self.m.items():
            T = T * {(0, 1): M[:, :, None], (1, 2): M[None, :, :], (0, 2): M[:, None, :]}[(s, t)]
        return T


PERMS = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]


def sym3(terms):
    return [Tri(t.permute(p).v, t.permute(p).m, t.c / 6.0) for t in terms for p in PERMS]


def _vec(x, n):
    return np.ones(n) if x is None else x


def _mat(M, n):
    return np.ones((n, n)) if M is None else M


def tri_full(t, W, block=32):
    """sum_ijk W_ia W_ja W_kb X_ijk."""
    n = W.shape[0]
    f, g, h = (_vec(x, n) for x in t.v)
    A, B, D = t.m.get((0, 1)), t.m.get((1, 2)), t.m.get((0, 2))
    if A is None:
        P = W.T @ (f[:, None] * _mat(D, n)); Q = W.T @ (g[:, None] * _mat(B, n))      # (a, k)
        return t.c * ((P * Q) @ (h[:, None] * W))
    if B is None:
        P = A @ (g[:, None] * W); Q = _mat(D, n) @ (h[:, None] * W)                     # (i, a), (i, b)
        return t.c * ((W * f[:, None] * P).T @ Q)
    if D is None:
        P = A.T @ (f[:, None] * W); Q = B @ (h[:, None] * W)                            # (j, a), (j, b)
        return t.c * ((W * g[:, None] * P).T @ Q)
    # full triangle: Y[a,k] = sum_ij (W_ia f_i) (W_ja g_j) A_ij B_jk D_ik, one GEMM per a (fp32)
    A32, B32, D32 = A.astype(np.float32), B.astype(np.float32), D.astype(np.float32)
    U = (W * f[:, None]).astype(np.float32); V = (W * g[:, None]).astype(np.float32)
    Y = np.zeros((n, n))
    for a0 in range(0, n, block):
        a1 = min(n, a0 + block)
        Mb = np.matmul(A32[None, :, :] * V.T[a0:a1, None, :], B32[None])               # (b, i, k)
        Y[a0:a1] = np.einsum("bik,ib,ik->bk", Mb, U[:, a0:a1], D32, optimize=True)
    return t.c * (Y @ (h[:, None] * W))


def restrictions(t):
    """X_iik, X_kii, X_iki as (i, k) matrices and X_iii."""
    n = len(next(x for x in t.v if x is not None)) if any(x is not None for x in t.v) else next(iter(t.m.values())).shape[0]
    f, g, h = (_vec(x, n) for x in t.v)
    A, B, D = (_mat(t.m.get(k), n) for k in ((0, 1), (1, 2), (0, 2)))
    dA, dB, dD = np.diag(A), np.diag(B), np.diag(D)
    XA = t.c * (f * g * dA)[:, None] * h[None, :] * B * D                 # X_iik
    XB = t.c * f[None, :] * (g * h * dB)[:, None] * A.T * D.T             # X_kii : f_k g_i h_i A_ki B_ii D_ki
    XC = t.c * (f * h * dD)[:, None] * g[None, :] * A * B.T               # X_iki : f_i g_k h_i A_ik B_ki D_ii
    xd = t.c * f * g * h * dA * dB * dD
    return XA, XB, XC, xd


def slice_transport(XA, XB, XC, xd, W):
    """transport of the non-all-distinct part of a tensor given its restrictions (inclusion-exclusion)."""
    W2 = W * W
    S = W2.T @ XA @ W
    S += (W * (XB @ W)).T @ W
    S += (W * (XC @ W)).T @ W
    S -= 2.0 * (W2 * xd[:, None]).T @ W
    return S


def T_ad(terms, W):
    """transport of all_distinct(sum of triangle terms)."""
    out = 0.0
    for t in terms:
        out = out + tri_full(t, W) - slice_transport(*restrictions(t), W)
    return out


# ----------------------------------------------------------------------------------------- per-layer terms
def hermite_coeffs(mu, var, pmax):
    from numpy.polynomial.hermite_e import hermeval
    sig = np.sqrt(var); alpha = mu / sig
    phi = np.exp(-0.5 * alpha ** 2) / np.sqrt(2 * pi)
    PhiG = 0.5 * (1.0 + np.vectorize(erf)(alpha / sqrt(2)))
    c = {1: sig * PhiG}
    for p in range(2, pmax + 1):
        c[p] = sig * hermeval(-alpha, [0] * (p - 2) + [1]) * phi / factorial(p)
    return c


def hermite_terms(C, mu, var, degree):
    """triangle terms of the Gaussian Hermite expansion of total degree exactly `degree` (oracle_k3.hermite_model)."""
    sig = np.sqrt(var)
    rho = C / np.outer(sig, sig); np.fill_diagonal(rho, 0.0)
    c = hermite_coeffs(mu, var, degree - 2)
    out = []
    for p in range(1, degree - 1):
        for q in range(1, degree - 1):
            r = degree - p - q
            if r < 1:
                continue
            al, be, ga = (p + q - r) // 2, (q + r - p) // 2, (p + r - q) // 2
            if min(al, be, ga) < 0 or (p + q - r) % 2:
                continue
            coef = factorial(p) * factorial(q) * factorial(r) / (factorial(al) * factorial(be) * factorial(ga))
            out.append(Tri([c[p], c[q], c[r]], {(0, 1): rho ** al if al else None, (1, 2): rho ** be if be else None,
                                                (0, 2): rho ** ga if ga else None}, coef))
    return out


def layer_terms(st, l, W, with_h8=True, log=None):
    """all candidate transported terms (n x n, D21(l+1) space) of layer l."""
    Wn = W[l + 1].astype(np.float64)
    n = Wn.shape[0]
    mu, var, Phi, C = st["mu"][l], st["var"][l], st["Phi"][l], st["C"][l]
    sig = np.sqrt(var); alpha = mu / sig
    phi = np.exp(-0.5 * alpha ** 2) / np.sqrt(2 * pi)
    w2 = phi / sig; w3 = -alpha * phi / sig ** 2; w5 = (alpha ** 3 - 3 * alpha) * (-1) * phi / sig ** 4
    D21z, D3z = st["D21z"][l], st["D3z"][l]
    Co = C.copy(); np.fill_diagonal(Co, 0.0)
    K22 = st["G22"][l].copy(); np.fill_diagonal(K22, 0.0)
    out = {"D21": st["D21z"][l + 1].copy()}
    t0 = time.time()
    # slices of kappa3(a): restrictions all = D21a, diagonal D3a
    D21a = st["D21a"][l]
    out["slices"] = slice_transport(D21a, D21a, D21a, st["D3a"][l], Wn)
    # B0: Phi^3 kappa3(z), all-distinct
    P2 = Phi * Phi
    XA = P2[:, None] * Phi[None, :] * D21z
    out["B0"] = st["TB0"][l] - slice_transport(XA, XA, XA, Phi ** 3 * D3z, Wn)
    # leading Wick (gaussian rho^2 with empirical Phi) + B0
    wick = [Tri([Phi, Phi, w2], {(0, 2): C, (1, 2): C}), Tri([Phi, w2, Phi], {(0, 1): C, (1, 2): C.T}),
            Tri([w2, Phi, Phi], {(0, 1): C.T, (0, 2): C.T})]
    out["wick"] = T_ad(wick, Wn) + out["B0"]
    out["H4"] = T_ad(hermite_terms(C, mu, var, 4), Wn)
    out["H6"] = T_ad(hermite_terms(C, mu, var, 6), Wn)
    if with_h8:
        out["H8"] = T_ad(hermite_terms(C, mu, var, 8), Wn)
    if log:
        log(f"  layer {l}: gaussian terms {time.time()-t0:.0f}s")
    out["B1"] = T_ad(sym3([Tri([w2, Phi, w2], {(0, 2): D21z, (1, 2): Co})]), Wn)
    out["B2"] = T_ad(sym3([Tri([w3, Phi, Phi], {(0, 2): D21z, (0, 1): Co})]), Wn)
    out["B3"] = T_ad(sym3([Tri([w5 * D3z, Phi, Phi], {(0, 1): Co, (0, 2): Co})]), Wn)
    out["B4"] = out["H6"]
    out["B5"] = T_ad(sym3([Tri([w3, w2, Phi], {(0, 1): K22, (0, 2): Co})]), Wn)
    # B6: sym3(w2_i Phi_j Phi_k kappa4(z)_iijk) all-distinct = 2/3 T_ad(Y0) + 1/3 T_ad(Y2)
    F_iik = st["F211"][l]; F_kii = st["G22"][l]; F_iii = np.diag(F_kii)       # F_kii[i,k] = F_{k i i}; G22 symmetric
    wv = w2 * var
    gauss0 = [Tri([wv, Phi, Phi], {(1, 2): C}), Tri([w2, Phi, Phi], {(0, 1): C, (0, 2): C}, 2.0)]
    gauss2 = [Tri([Phi, Phi, wv], {(0, 1): C}), Tri([Phi, Phi, w2], {(0, 2): C, (1, 2): C}, 2.0)]
    full0 = st["TS1"][l] - sum(tri_full(t, Wn) for t in gauss0)
    full2 = st["TS2"][l] - sum(tri_full(t, Wn) for t in gauss2)
    Y0 = slice_transport((w2 * Phi)[:, None] * Phi[None, :] * F_iik, w2[None, :] * P2[:, None] * F_kii.T,
                         (w2 * Phi)[:, None] * Phi[None, :] * F_iik, w2 * P2 * F_iii, Wn)
    Y2 = slice_transport(w2[None, :] * P2[:, None] * F_kii.T, (w2 * Phi)[:, None] * Phi[None, :] * F_iik,
                         (w2 * Phi)[:, None] * Phi[None, :] * F_iik, w2 * P2 * F_iii, Wn)
    out["B6"] = (2.0 / 3.0) * (full0 - Y0) + (1.0 / 3.0) * (full2 - Y2)
    # regenerated (2,1,1) slice u_i C_jk (one least-squares u per doubled index, on C_off)
    cc = float(np.sum(Co * Co))
    full_u = st["Uq"][l] - var * float(np.sum(C * Co)) - 2.0 * np.einsum("ij,jk,ki->i", C, Co, C)
    u = (full_u - 2.0 * np.sum(F_iik * Co, 1)) / cc
    out["u"] = u
    out["B6reg"] = T_ad(sym3([Tri([w2 * u, Phi, Phi], {(1, 2): Co})]), Wn)
    if log:
        log(f"  layer {l}: all terms {time.time()-t0:.0f}s")
    return out


# ----------------------------------------------------------------------------------------- ladder
BASIS = ["B0", "B1", "B2", "B3", "B4", "B5", "B6"]


def models(t, coef=None, coef_reg=None):
    """the ladder's D21(l+1) predictions from a layer's term dict; coef: D21-space fitted coefficients."""
    m = {}
    m["memless"] = t["slices"]
    m["wick"] = t["slices"] + t["wick"]
    m["herm2"] = t["slices"] + t["H4"] + t["B0"]
    m["herm3"] = m["herm2"] + t["H6"]
    if "H8" in t:
        m["herm4"] = m["herm3"] + t["H8"]
    base = t["slices"] + t["H4"]
    m["closure_noK4"] = base + sum(c * t[b] for c, b in zip(CLOSURE_COEF[:6], BASIS[:6]))
    m["closure"] = m["closure_noK4"] + CLOSURE_COEF[6] * t["B6"]
    m["closure_reg"] = m["closure_noK4"] + CLOSURE_COEF[6] * t["B6reg"]
    if coef is not None:
        m["fit"] = m["wick"] + sum(c * t[b] for c, b in zip(coef, BASIS))
    if coef_reg is not None:
        m["fit_reg"] = m["wick"] + sum(c * t[b] for c, b in zip(coef_reg, BASIS[:6] + ["B6reg"]))
    return m


def fit_coef(t, reg=False):
    """least squares of the D21(l+1) residual after the leading Wick model on the transported basis."""
    names = BASIS[:6] + (["B6reg"] if reg else ["B6"])
    X = np.stack([t[b].ravel() for b in names], 1); y = (t["D21"] - t["slices"] - t["wick"]).ravel()
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    return coef


def rel(a, b, off=False):
    if off:
        a = a.copy(); b = b.copy(); np.fill_diagonal(a, 0.0); np.fill_diagonal(b, 0.0)
    return float(np.sqrt(np.sum((a - b) ** 2) / np.sum(b ** 2)))


COLS = ["memless", "wick", "herm2", "herm3", "herm4", "closure_noK4", "closure", "closure_reg", "fit", "fit_reg"]


def ladder_single(st, W, cache=None, log=None, with_h8=True):
    L = W.shape[0]
    rows = []
    for l in range(L - 1):
        t = layer_terms(st, l, W, with_h8=with_h8, log=log)
        if cache is not None:
            cache[l] = t
        coef, coef_r = fit_coef(t), fit_coef(t, reg=True)
        m = models(t, coef, coef_r)
        rows.append(dict(l=l, **{k: rel(v, t["D21"]) for k, v in m.items()}, coef=coef, coef_reg=coef_r))
    return rows


def ladder_pair(tA, tB):
    """per layer: noise of D21, oracle-style cross error (A model vs B target, target-noise-subtracted) and the
    replica cross-product estimate eps^2 = <M_A - D_A, M_B - D_B> / <D_A, D_B> (noise-free in expectation)."""
    rows = []
    for l in sorted(tA):
        a, b = tA[l], tB[l]
        DA, DB = a["D21"], b["D21"]
        e_noise = rel(DA, DB) / sqrt(2)
        cA, cB = fit_coef(a), fit_coef(b)
        crA, crB = fit_coef(a, True), fit_coef(b, True)
        mA, mB = models(a, cA, crA), models(b, cB, crB)
        den = float(np.sum(DA * DB))
        r = dict(l=l, noise=e_noise, coefA=cA, coefB=cB, coefrA=crA)
        for k in mA:
            ex = rel(mA[k], DB)
            r[k + "_x"] = ex
            r[k + "_xc"] = sqrt(max(ex ** 2 - e_noise ** 2, 0.0))
            cp = float(np.sum((mA[k] - DA) * (mB[k] - DB))) / den
            r[k + "_cp"] = np.sign(cp) * sqrt(abs(cp))
        # noise of the transported Phi^3 kappa3(z) and kappa4 terms relative to ||D21||
        for k in ("B0", "B6"):
            r[k + "_noise"] = float(np.sqrt(np.sum((a[k] - b[k]) ** 2) / 2.0 / np.sum(DB ** 2)))
        rows.append(r)
    return rows


# ----------------------------------------------------------------------------------------- verification
def verify(atlas_path, stats_path):
    import oracle_k3 as ok
    z = np.load(atlas_path); st = dict(np.load(stats_path))
    W = z["weights"].astype(np.float64)
    L, n, _ = W.shape
    print(f"verify {atlas_path} vs {stats_path}: width {n}, depth {L}, N = {int(z['n_samples'])} / {int(st['n_samples'])}")
    print("relative differences, streamed vs dense (per term, in D21(l+1) space; columns: eps dense, eps streamed)")
    for l in range(L - 1):
        o = ok.layer_objects(z, l)
        Wn = W[l + 1]
        t = layer_terms(st, l, z["weights"])
        K3m = ok.slices_only(o["K3a"])
        dense = {"D21": o["D21"], "slices": ok.transport_d21(K3m, Wn)}
        Kw = ok.wick_model(o["C"], o["Phi"], o["w2"], o["K3z"])
        dense["wick"] = ok.transport_d21(Kw, Wn)
        B = ok.residual_basis(o)
        for i, b in enumerate(B):
            dense[f"B{i}"] = ok.transport_d21(b, Wn)
        dense["H4"] = ok.transport_d21(ok.hermite_model(o["C"], o["mu"], o["var"], 4), Wn)
        dense["H6"] = ok.transport_d21(ok.hermite_model(o["C"], o["mu"], o["var"], 6), Wn) - dense["H4"]
        dense["H8"] = ok.transport_d21(ok.hermite_model(o["C"], o["mu"], o["var"], 8), Wn) - dense["H4"] - dense["H6"]
        K211 = o["K211"]; Co = ok.offdiag(o["C"])
        u = np.einsum("ijk,jk->i", K211, Co) / float(np.sum(Co * Co))
        o_reg = dict(o); o_reg["K211"] = ok.all_distinct(np.einsum("i,jk->ijk", u, Co))
        dense["B6reg"] = ok.transport_d21(ok.residual_basis(o_reg)[6], Wn)
        diffs = " ".join(f"{k}:{rel(t[k], dense[k]):.1e}" for k in dense) + f" u:{rel(t['u'], u):.1e}"
        # ladder columns: dense oracle definitions vs streamed
        e_d = {"memless": rel(dense["slices"], o["D21"]), "wick": rel(dense["slices"] + dense["wick"], o["D21"]),
               "herm4": rel(ok.transport_d21(K3m + ok.hermite_model(o["C"], o["mu"], o["var"], 8, o["K3z"]), Wn), o["D21"]),
               "closure": rel(ok.transport_d21(K3m + ok.closure_model(o), Wn), o["D21"]),
               "closure_reg": rel(ok.transport_d21(K3m + ok.closure_model(o_reg), Wn), o["D21"])}
        coef, fit, _ = ok.fit_residual(ok.all_distinct(o["K3a"]) - Kw, B)
        e_d["fit(tensor)"] = rel(ok.transport_d21(K3m + Kw + fit, Wn), o["D21"])
        m = models(t, fit_coef(t), fit_coef(t, True))
        e_s = {k: rel(m[k], t["D21"]) for k in ("memless", "wick", "herm4", "closure", "closure_reg")}
        e_s["fit(D21)"] = rel(m["fit"], t["D21"])
        print(f"{l:>2} {diffs}")
        print("   " + "  ".join(f"{k} {e_d[k]:.4f}/{e_s.get(k, e_s.get('fit(D21)')):.4f}" for k in e_d))


# ----------------------------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd")
    s = sub.add_parser("stats")
    s.add_argument("--seed", type=int, required=True); s.add_argument("--width", type=int, default=1024)
    s.add_argument("--depth", type=int, default=16); s.add_argument("--n-samples", type=int, required=True)
    s.add_argument("--chunk", type=int, default=1024); s.add_argument("--sample-seed", type=int, required=True)
    s.add_argument("--out", required=True)
    v = sub.add_parser("verify"); v.add_argument("atlas"); v.add_argument("stats")
    g = sub.add_parser("ladder"); g.add_argument("stats", nargs="+"); g.add_argument("--no-h8", action="store_true")
    g.add_argument("--out", default=None, help="write a json of the rows")
    args = ap.parse_args()
    log = lambda m: print(time.strftime("%H:%M:%S"), m, flush=True)
    if args.cmd == "stats":
        W = mlp_weights(args.seed, args.width, args.depth)
        st = stream_stats(W, args.n_samples, args.chunk, args.sample_seed, log=log)
        st.update(mlp_seed=args.seed, width=args.width, depth=args.depth)
        np.savez(args.out, **{k: (v.astype(np.float32) if isinstance(v, np.ndarray) and v.ndim == 3 else v) for k, v in st.items()})
        log(f"saved {args.out}")
    elif args.cmd == "verify":
        verify(args.atlas, args.stats)
    elif args.cmd == "ladder":
        run_ladder(args.stats, with_h8=not args.no_h8, out=args.out, log=log)


def load_stats(p):
    z = np.load(p)
    st = {k: (z[k].astype(np.float64) if z[k].dtype == np.float32 else z[k]) for k in z.files}
    W = mlp_weights(int(st["mlp_seed"]), int(st["width"]), int(st["depth"]))
    return st, W


def run_ladder(paths, with_h8=True, out=None, log=None):
    import json
    stA, W = load_stats(paths[0])
    cA = {}
    rows = ladder_single(stA, W, cache=cA, with_h8=with_h8, log=log)
    print(f"\n{paths[0]}: width {W.shape[1]}, N = {int(stA['n_samples'])}  (single replica, within-sample eps)")
    cols = [c for c in COLS if c in rows[0]]
    print(f"{'l':>2} " + " ".join(f"{c:>12}" for c in cols) + "  fit coefficients B0..B6")
    for r in rows:
        print(f"{r['l']:>2} " + " ".join(f"{r[c]:12.4f}" for c in cols) + "  " + " ".join(f"{c:+.2f}" for c in r["coef"]), flush=True)
    res = dict(single=[{k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in r.items()} for r in rows])
    if len(paths) > 1:
        stB, WB = load_stats(paths[1])
        assert np.array_equal(W, WB)
        cB = {}
        ladder_single(stB, W, cache=cB, with_h8=with_h8, log=log)
        prow = ladder_pair(cA, cB)
        print(f"\npair {paths[0]} | {paths[1]}: noise-corrected eps of D21(l+1)")
        print("  _xc = model from A vs target B, target noise subtracted (oracle_k3 --pair convention)")
        print("  _cp = replica cross-product <M_A-D_A, M_B-D_B>/<D_A,D_B> (model and target noise removed)")
        print(f"{'l':>2} {'noise':>6} {'B0n':>6} {'B6n':>6} | " + " ".join(f"{c[:11]:>11}" for c in cols) + " | xc: " + " ".join(f"{c[:11]:>11}" for c in cols))
        for r in prow:
            print(f"{r['l']:>2} {r['noise']:6.4f} {r['B0_noise']:6.4f} {r['B6_noise']:6.4f} | " + " ".join(f"{r[c + '_cp']:11.4f}" for c in cols)
                  + " | xc: " + " ".join(f"{r[c + '_xc']:11.4f}" for c in cols)
                  + "  coefA " + " ".join(f"{c:+.2f}" for c in r["coefA"]), flush=True)
        res["pair"] = [{k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in r.items()} for r in prow]
    if out:
        with open(out, "w") as f:
            json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
