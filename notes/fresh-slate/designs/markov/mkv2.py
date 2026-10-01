"""MKV-2: every site s carries its whole downstream influence as a function of its own scalar latent
xi_s (standardised g_s):  F_{s->j}(xi) = sum_i a_i[s, j] e_i^{(s)}(xi), e^{(s)} an orthonormal basis of
span{X, xi X, X^2} with 1 and xi projected out (X = ReLU residual).  Crossing a ReLU (layer l -> l+1) is the
Gaussian conditional expectation E[phi(z_b) | xi_s] to second order in Delta = c xi + F, minus the pure
c^2 xi^2 term (counted by the fresh site at b), Galerkin-projected back onto e^{(s)}; its He1 part is
folded into the link c (and thereby into the field covariance).  See DESIGN.md §7, item 6.
"""
import numpy as np
from scipy.special import ndtr
from mkv import relu_moments, relu_cov_offdiag, pdf

_GLx, _GLw = np.polynomial.legendre.leggauss(120)
NB = 3


def site_tables(mu, v):
    """Per-site orthonormal basis on quadrature nodes and the moment/projection tables."""
    n = len(mu); s = np.sqrt(v); t = mu / s
    m1 = mu * ndtr(t) + s * pdf(t); p = ndtr(t)
    u0 = np.clip(-t, -11.9, 11.9); lo, hi = -12.0, 12.0
    nodes, wts = [], []
    for a, b in [(np.full(n, lo), u0), (u0, np.full(n, hi))]:
        h = (b - a) / 2; m = (b + a) / 2
        nodes.append(m[:, None] + h[:, None] * _GLx[None]); wts.append(h[:, None] * _GLw[None])
    xi = np.concatenate(nodes, 1); w = np.concatenate(wts, 1) * pdf(xi)
    w = w / w.sum(1, keepdims=True)
    X = np.maximum(mu[:, None] + s[:, None] * xi, 0) - m1[:, None] - p[:, None] * s[:, None] * xi
    raw = [X, xi * X, X * X]
    E = lambda f: (w * f).sum(1)
    basis = []
    for f in raw:
        f = f - E(f)[:, None] - E(f * xi)[:, None] * xi
        for e in basis:
            f = f - E(f * e)[:, None] * e
        nr = np.sqrt(np.maximum(E(f * f), 0))
        ok = nr > 1e-7 * (s + 1e-30) ** (1 if len(basis) == 0 else 2)
        f = np.where(ok[:, None], f / np.where(ok, nr, 1)[:, None], 0.0)
        basis.append(f)
    e = np.stack(basis)                      # (NB, n, Q)
    normX = np.sqrt(E(X * X))
    h2 = xi * xi - 1; h3 = xi ** 3 - 3 * xi
    T = dict(normX=normX)
    T["M3"] = np.einsum('nq,inq,jnq,knq->nijk', w, e, e, e)
    T["M4"] = np.einsum('nq,inq,jnq,knq,lnq->nijkl', w, e, e, e, e)
    T["M2x"] = np.einsum('nq,inq,jnq->nij', w * xi, e, e)
    T["M3x"] = np.einsum('nq,inq,jnq,knq->nijk', w * xi, e, e, e)
    T["Mh2"] = np.einsum('nq,inq->ni', w * h2, e)
    T["M2h2"] = np.einsum('nq,inq,jnq->nij', w * h2, e, e)
    T["Mh3"] = np.einsum('nq,inq->ni', w * h3, e)
    T["Pxi"] = np.einsum('nq,inq,knq->nik', w * xi, e, e)     # proj of xi e_k on e_i
    T["L1"] = np.einsum('nq,knq->nk', w * xi * xi, e)          # He1 coef of xi e_k
    T["P2"] = T["M3"]                                          # proj of e_j e_k on e_i
    T["L2"] = T["M2x"]                                         # He1 coef of e_j e_k
    return T


def cumulants(a, c, T, R=None):
    """Per-target delta k3, k4 summed over sites. a: (NB, S, J), c: (S, J) std links, T: site tables."""
    EF2 = (a * a).sum(0)
    F3 = np.einsum('sijk,isn,jsn,ksn->sn', T["M3"], a, a, a)
    F2x = np.einsum('sij,isn,jsn->sn', T["M2x"], a, a)
    Fh2 = np.einsum('si,isn->sn', T["Mh2"], a)
    F4 = np.einsum('sijkl,isn,jsn,ksn,lsn->sn', T["M4"], a, a, a, a)
    F3x = np.einsum('sijk,isn,jsn,ksn->sn', T["M3x"], a, a, a)
    F2h2 = np.einsum('sij,isn,jsn->sn', T["M2h2"], a, a)
    Fh3 = np.einsum('si,isn->sn', T["Mh3"], a)
    d3 = (F3 + 3 * c * F2x + 3 * c * c * Fh2).sum(0)
    d4 = (F4 - 3 * EF2 ** 2 + 4 * c * F3x + 6 * c * c * F2h2 + 4 * c ** 3 * Fh3).sum(0)
    if R is not None:   # two-site kappa4 trees within the source layer, latent correlation R (zero diag)
        Q1 = F2x; Q2 = c * Fh2
        RQ1, RQ2 = R @ Q1, R @ Q2
        d4 = d4 + 3 * (Q1 * RQ1).sum(0) + 12 * (Q2 * RQ1).sum(0) + 12 * (Q2 * RQ2).sum(0)
    return d3, d4


def mkv2(W, w=16, K=60, return_all=False, two_site=True, curvature=True, curv_age=99):
    W = np.asarray(W, dtype=np.float64)
    L, n, _ = W.shape
    mu = np.zeros(n); C = W[0].T @ W[0]; k3 = np.zeros(n); k4 = np.zeros(n)
    k3old = np.zeros(n); k4old = np.zeros(n)
    sources = []
    means = np.zeros((L, n)); diag = []
    for l in range(L):
        v = np.diag(C).copy(); sd = np.sqrt(v)
        m1, m2, p = relu_moments(mu, v, k3, k4)
        means[l] = m1
        diag.append(dict(mu=mu.copy(), v=v, k3=k3.copy(), k4=k4.copy()))
        if l == L - 1:
            break
        Wn = W[l + 1]
        dens = pdf(mu / sd) / sd
        U = Wn * dens[:, None]; V = Wn * p[:, None]
        # Gaussian field
        Ca = relu_cov_offdiag(mu, C, K); np.fill_diagonal(Ca, m2 - m1 ** 2)
        dCa = np.zeros((n, n)); tdiag = np.zeros(n)
        newsrcs = []
        for sd_ in sources:
            a, c, T = sd_["a"], sd_["c"], sd_["T"]
            # (2,1)-slice of z_l from this source: kappa(z_b, z_b, z_d), factorised in d
            Lk = np.einsum('sijk,isb,jsb->ksb', T["M3"], a, a) + 2 * c[None] * np.einsum('sik,isb->ksb', T["M2x"], a) \
                + (c * c)[None] * T["Mh2"].T[:, :, None]
            Lc = np.einsum('sij,isb,jsb->sb', T["M2x"], a, a) + 2 * c * np.einsum('si,isb->sb', T["Mh2"], a)
            dCa += np.einsum('ksb,ksd->bd', Lk, a) + Lc.T @ c
            tdiag += (Lk * a).sum((0, 1)) + (Lc * c).sum(0)
            # transport across the ReLU at layer l
            lin = a * p[None, None, :]
            if curvature and sd_["age"] <= curv_age:
                quad = c[None] * np.einsum('sik,ksb->isb', T["Pxi"], a) \
                    + 0.5 * np.einsum('sijk,jsb,ksb->isb', T["P2"], a, a)
                he1 = c * np.einsum('sk,ksb->sb', T["L1"], a) + 0.5 * np.einsum('sjk,jsb,ksb->sb', T["L2"], a, a)
                pre = lin + quad * dens[None, None, :]
                cpre = c * p[None, :] + he1 * dens[None, :]
            else:
                pre = lin; cpre = c * p[None, :]
            sd_["a"] = np.einsum('isb,bj->isj', pre, Wn); sd_["c"] = cpre @ Wn; sd_["age"] += 1
            newsrcs.append(sd_)
        # non-Gaussian (2,1) correction of Cov(a_l), off-diagonal
        dCa = 0.5 * (dCa * dens[:, None] * p[None, :]); dCa = dCa + dCa.T
        np.fill_diagonal(dCa, 0.0)
        Cn = Wn.T @ (Ca + dCa) @ Wn
        mun = m1 @ Wn
        # fresh sites at layer l
        T = site_tables(mu, v)
        a0 = np.zeros((NB, n, n)); a0[0] = T["normX"][:, None] * Wn
        R = C / np.outer(sd, sd); np.fill_diagonal(R, 0.0)
        newsrcs.append(dict(a=a0, c=(C * p[None, :]) @ Wn / sd[:, None], T=T, age=1, R=R if two_site else None))
        Wp = Wn * p[:, None]
        k3old = (Wp ** 3).T @ k3old; k4old = (Wp ** 4).T @ k4old
        k3 = k3old.copy(); k4 = k4old.copy(); keep = []
        for sd_ in newsrcs:
            d3, d4 = cumulants(sd_["a"], sd_["c"], sd_["T"], sd_["R"])
            k3 += d3; k4 += d4
            if sd_["age"] >= w:
                k3old += d3; k4old += d4
            else:
                keep.append(sd_)
        sources = keep
        mu, C = mun, Cn
    return (means, diag) if return_all else means
