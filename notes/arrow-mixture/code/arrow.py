# Arrow-mixture closure ("bundle closure"): the law of each layer is carried as a discrete base measure
# along k directions (the cloud of component means) times a shared Gaussian fiber.  Components born at a
# layer persist through w-1 further ReLU steps (the path memory / border depth) and are then merged by
# moment matching, their spread returned to the fiber.  w = 1 is the ordinary Gaussian closure.
#   base: 'eig'    top-k eigenvectors of the fiber covariance (value quadrature, q Gauss-Hermite nodes)
#         'mean'   the mean direction (gain mode), k = 1
#         'neuron' top-k neurons by single-wall kappa3, split by SIGN (two truncated-Gaussian children)
import numpy as np, sys, itertools
from math import factorial, sqrt, pi
from scipy.stats import norm
from numpy.polynomial.hermite_e import hermegauss


def relu_stats(mu, sig, K):
    """A[m] = sig^m E relu^(m)(z), z ~ N(mu, sig^2) (A[0] = mean); also var and third central moment."""
    a = mu / sig; ph = norm.pdf(a); Ph = norm.cdf(a)
    A = np.empty((K + 1,) + np.shape(mu))
    A[0] = mu * Ph + sig * ph
    A[1] = sig * Ph
    He = [np.ones_like(a), a]
    for j in range(2, K):
        He.append(a * He[-1] - (j - 1) * He[-2])
    for m in range(2, K + 1):
        A[m] = sig * ph * ((-1) ** m) * He[m - 2]
    E2 = (mu * mu + sig * sig) * Ph + mu * sig * ph
    E3 = (mu ** 3 + 3 * mu * sig * sig) * Ph + (mu * mu + 2 * sig * sig) * sig * ph
    var = E2 - A[0] ** 2
    k3 = E3 - 3 * A[0] * E2 + 2 * A[0] ** 3
    return A, var, k3


def trunc_moments(c):
    """standard normal truncated to (-inf, c) and (c, inf): mass, mean, variance (arrays over c)."""
    ph = norm.pdf(c); Ph = norm.cdf(c)
    pL = np.clip(Ph, 1e-300, None); pR = np.clip(1 - Ph, 1e-300, None)
    mL = -ph / pL; vL = 1 - c * ph / pL - mL ** 2
    mR = ph / pR; vR = 1 + c * ph / pR - mR ** 2
    return (pL, mL, np.clip(vL, 1e-12, None)), (pR, mR, np.clip(vR, 1e-12, None))


def arrow(Ws, k=0, q=3, w=2, base="eig", K=10, maxM=4096, verbose=False, gain=False, lev=None):
    """gain: carry the GAC scale mixture on top (the base then has a radial coordinate with infinite memory);
    lev: total-level sparse grid for the k-direction product quadrature (keep children with sum |j - centre| <= lev)."""
    if gain:
        from gac import inject, EG
    n = Ws[0].shape[0]; L = len(Ws)
    wts = np.ones(1); Mh = np.zeros((1, n)); Cf = np.eye(n); labels = [()]
    gam = 0.0; g1 = 1.0; prev = None
    gh_x, gh_w = hermegauss(q); gh_w = gh_w / gh_w.sum()
    out = []
    for l, W in enumerate(Ws):
        mu = Mh @ W.T                      # (M, n) component pre-activation means
        S = W @ Cf @ W.T                   # shared fiber
        S = 0.5 * (S + S.T)
        if gain:
            mbar = wts @ Mh; mu0 = W @ mbar
            S0 = S + W @ ((np.sqrt(wts)[:, None] * (Mh - mbar)).T @ (np.sqrt(wts)[:, None] * (Mh - mbar))) @ W.T
            if l > 0:
                gam = gam + inject(W, mu0, S0, prev); g1_new = EG(gam)
                # GAC bookkeeping: marginal moments are fixed (E z = g1 mu, E zz^T = S + mu mu^T), so the new gain
                # variance rescales the conditional means by rho = g1_old/g1_new and returns (1-rho^2) E[mu mu^T] to the fiber
                rho = g1 / g1_new
                S = S + (1 - rho * rho) * ((np.sqrt(wts)[:, None] * mu).T @ (np.sqrt(wts)[:, None] * mu))
                mu = rho * mu; g1 = g1_new
                mu0 = rho * mu0; S0 = S + W @ ((np.sqrt(wts)[:, None] * (Mh - mbar)).T @ (np.sqrt(wts)[:, None] * (Mh - mbar))) @ W.T * (rho * rho)
            sig0 = np.sqrt(np.diag(S0)); prev = (mu0, sig0, S0 / np.outer(sig0, sig0))
        # ---- split along k directions (sequential conditioning of the fiber) ----
        if k > 0:
            if base == "eig":
                ev, U = np.linalg.eigh(S); dirs = [U[:, -1 - d] for d in range(k)]
            elif base == "mean":
                mbar = wts @ mu
                if np.linalg.norm(mbar) < 1e-8:
                    ev, U = np.linalg.eigh(S); dirs = [U[:, -1]]
                else:
                    dirs = [mbar / np.linalg.norm(mbar)]
            elif base == "rand":
                rng = np.random.default_rng(1000 + l); Q, _ = np.linalg.qr(rng.standard_normal((n, k))); dirs = [Q[:, d] for d in range(k)]
            elif base == "neuron":
                mbar = wts @ mu; sig0 = np.sqrt(np.diag(S))
                _, _, k3 = relu_stats(mbar, sig0, 2)
                idx = np.argsort(-np.abs(k3))[:k]
                dirs = [np.eye(n)[i] for i in idx]
            elif base == "sign":
                mbar = wts @ mu; sig0 = np.sqrt(np.diag(S))
                _, _, k3 = relu_stats(mbar, sig0, 2)
                idx = np.argsort(-np.abs(k3))[:k]
                dirs = [np.eye(n)[i] for i in idx]
            else:
                raise ValueError(base)
            # sequential conditioning: (lam_d, r_d) with the fiber reduced after each direction
            lams, rs = [], []
            for u in dirs:
                lam = float(u @ S @ u); r = (S @ u) / lam
                lams.append(lam); rs.append(r)
                if base == "sign":
                    c = -(mu @ u) / sqrt(lam)
                    (pL, mL, vL), (pR, mR, vR) = trunc_moments(c)
                    new_w = np.concatenate([wts * pL, wts * pR])
                    new_mu = np.concatenate([mu + (mL * sqrt(lam))[:, None] * r[None, :],
                                             mu + (mR * sqrt(lam))[:, None] * r[None, :]])
                    vbar = float(np.concatenate([wts * pL * vL, wts * pR * vR]).sum() / new_w.sum())
                    labels = [lab + (0,) for lab in labels] + [lab + (1,) for lab in labels]
                    wts, mu = new_w, new_mu
                    S = S - lam * (1 - vbar) * np.outer(r, r)
                else:
                    S = S - lam * np.outer(r, r)
            if base != "sign":
                # product Gauss-Hermite grid over the k directions; optional total-level sparse grid with the
                # per-direction variance restored by stretching the kept nodes (second moments stay exact)
                grid = np.array(list(itertools.product(range(q), repeat=len(dirs))))          # (G, k)
                gw = np.prod(gh_w[grid], axis=1); gx = gh_x[grid]
                if lev is not None:
                    keep = np.abs(grid - (q - 1) // 2).sum(1) <= lev
                    grid, gw, gx = grid[keep], gw[keep], gx[keep]
                    gw = gw / gw.sum()
                    m2 = (gw[:, None] * gx * gx).sum(0)
                    gx = gx / np.sqrt(m2)[None, :]
                off = gx @ (np.sqrt(np.array(lams))[:, None] * np.array(rs))                  # (G, n)
                wts = (wts[:, None] * gw[None, :]).ravel()
                mu = (mu[:, None, :] + off[None, :, :]).reshape(-1, n)
                labels = [lab + tuple(g) for lab in labels for g in grid]
            # fold this layer's k split indices into one label entry
            kk = len(dirs)
            labels = [lab[:-kk] + (lab[-kk:],) for lab in labels]
        else:
            labels = [lab + ((),) for lab in labels]
        # ---- ReLU per component, shared reduced fiber ----
        sig = np.sqrt(np.clip(np.diag(S), 1e-12, None)); R = S / np.outer(sig, sig)
        A, var, _ = relu_stats(mu, sig[None, :], K)      # (K+1, M, n)
        post = A[0]                                      # (M, n) component post means
        sw = np.sqrt(wts)[:, None]
        Cpost = np.zeros((n, n)); Rk = np.ones((n, n))
        for m in range(1, K + 1):
            Rk = Rk * R
            G = (sw * A[m]).T @ (sw * A[m])
            Cpost += G * Rk / factorial(m)
        np.fill_diagonal(Cpost, wts @ var)
        out.append(g1 * (wts @ post))
        # ---- merge: drop the oldest split beyond memory w, moment-match groups ----
        if len(labels[0]) > w:
            labels = [lab[1:] for lab in labels]
        groups = {}
        for j, lab in enumerate(labels):
            groups.setdefault(lab, []).append(j)
        if len(groups) < len(labels):
            gw = np.array([wts[g].sum() for g in groups.values()])
            gm = np.array([(wts[g] @ post[g]) / wts[g].sum() for g in groups.values()])
            D = np.concatenate([np.sqrt(wts[g])[:, None] * (post[g] - gm[i][None, :]) for i, g in enumerate(groups.values())])
            Cpost = Cpost + D.T @ D
            wts, Mh, labels = gw, gm, list(groups.keys())
        else:
            Mh = post
        if len(wts) > maxM:
            raise RuntimeError(f"{len(wts)} components > maxM")
        Cf = Cpost
        if verbose:
            print(f"layer {l}: {len(wts)} components", flush=True)
    return np.array(out)


if __name__ == "__main__":
    # usage: arrow.py <W.npy> <truth.npz> k q w base
    Ws = list(np.load(sys.argv[1])); mt = np.load(sys.argv[2])["m"]
    k, q, w, base = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    out = arrow(Ws, k=k, q=q, w=w, base=base)
    ref = arrow(Ws, k=0)
    mse = lambda a, b: np.mean((a - b) ** 2)
    print(f"k={k} q={q} w={w} base={base}: final MSE {mse(out[-1], mt[-1]):.3e}  (closure {mse(ref[-1], mt[-1]):.3e})")
