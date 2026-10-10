"""Exact or symbolic checks of the finite claims in "Approximation bridges" (relayed 10 Oct 2026), and of Stage 25's
ungated-recursion identity. No network experiment; every check is a statement with a proof in the note.

  python notes/stage25/theory25/check_bridges.py"""
import itertools, math, numpy as np, sympy as sp
from fractions import Fraction as Fr
from scipy import integrate, special

rng = np.random.default_rng(25)
ok = lambda b: "PASS" if b else "FAIL"


def b1_split():
    # occupation dual of a Beta(1/2,1/2) energy split: P(k of 2 particles stay at i) = C(2,k)(1/2)_k(1/2)_(2-k)/(1)_2
    poch = lambda a, k: math.prod(Fr(a) + j for j in range(k))
    P = [math.comb(2, k) * poch(Fr(1, 2), k) * poch(Fr(1, 2), 2 - k) / poch(Fr(1), 2) for k in (2, 1, 0)]
    phi = sp.symbols("phi"); E = lambda e: sp.integrate(e, (phi, 0, 2 * sp.pi)) / (2 * sp.pi)
    mom = [E(sp.cos(phi) ** 4), E(sp.cos(phi) ** 2 * sp.sin(phi) ** 2)]
    return f"[B1a] double-occupancy split ii/ij/jj = {P[0]}/{P[1]}/{P[2]}; E cos^4 = {mom[0]}, E cos^2 sin^2 = {mom[1]}: {ok(P == [Fr(3, 8), Fr(1, 4), Fr(3, 8)] and mom == [sp.Rational(3, 8), sp.Rational(1, 8)])}"


def b1_quartic():
    out = []
    for d in (3, 4, 6):
        eta = sp.symbols(f"e0:{d}", positive=True); S = sum(e ** 2 for e in eta)
        HS = 0
        for i, j in itertools.combinations(range(d), 2):
            QS = S - eta[i] ** 2 - eta[j] ** 2 + sp.Rational(3, 4) * (eta[i] + eta[j]) ** 2
            HS += S - QS
        lhs = sp.expand(HS.subs(eta[-1], 1 - sum(eta[:-1])))
        rhs = sp.expand((sp.Rational(d + 2, 4) * (S - sp.Rational(3, d + 2))).subs(eta[-1], 1 - sum(eta[:-1])))
        out.append(sp.simplify(lhs - rhs) == 0)
    return f"[B1b] H_(K_d) (S - 3/(d+2)) = (d+2)/4 (S - 3/(d+2)) on the simplex, d = 3, 4, 6: {ok(all(out))}"


def b3_ito():
    # finite state space: m(J, h) = E_p F, p ~ exp(phi'J phi/2 + h'phi). Ito drift for dh = C dB, dJ = -C^2 dt is
    # -grad_J m : C^2 + (1/2) Hess_h m : C^2; claim: equals -V' C^2 a with V = Cov(phi, F), a = E phi.
    S, D, q = 9, 4, 3; phi = rng.standard_normal((S, D)); F = rng.standard_normal((S, q))
    J = 0.1 * np.eye(D); h = 0.3 * rng.standard_normal(D); B = rng.standard_normal((D, D)); C = 0.5 * (B + B.T)
    def m(J, h):
        w = np.exp(0.5 * np.einsum("sd,de,se->s", phi, J, phi) + phi @ h); p = w / w.sum(); return p @ F, p
    _, p = m(J, h); a = p @ phi; V = (phi - a).T @ (p[:, None] * (F - p @ F))
    eps = 1e-4; C2 = C @ C
    gJ = (m(J + eps * C2, h)[0] - m(J - eps * C2, h)[0]) / (2 * eps)                       # grad_J m : C^2 (directional)
    Hh = np.zeros(q)
    for k in range(D):                                                                        # Hess_h m : C^2 = sum_k (C e_k)' Hess (C e_k)
        u = C[:, k]; Hh += (m(J, h + eps * u)[0] - 2 * m(J, h)[0] + m(J, h - eps * u)[0]) / eps ** 2
    drift = -gJ + 0.5 * Hh; claim = -V.T @ C2 @ a
    # zero-drift preservation: choose C with V'C = 0 (D = 4 > q = 3 leaves one direction)
    _, _, Vt = np.linalg.svd(V.T); c = Vt[-1]; C0 = np.outer(c, c)
    d0 = -(m(J + eps * C0 @ C0, h)[0] - m(J - eps * C0 @ C0, h)[0]) / (2 * eps) + 0.5 * (m(J, h + eps * c)[0] - 2 * m(J, h)[0] + m(J, h - eps * c)[0]) / eps ** 2
    return f"[B3] Ito drift of the final means = -V'C^2 a: max dev {np.max(np.abs(drift - claim)):.1e}; with V'C = 0 drift {np.max(np.abs(d0)):.1e}, martingale part V'C = {np.max(np.abs(V.T @ C0)):.1e}: {ok(np.max(np.abs(drift - claim)) < 1e-5 and np.max(np.abs(d0)) < 1e-5)}"


def b4_quadrature():
    M = 12; h = math.pi / math.sqrt(M); js = np.arange(-M, M + 1)
    def T(EZ, q, sig): return EZ / 2 + h * sig / math.pi * np.sum(np.exp(js * h) * q)
    laws, worst = [], 0.0
    for mu, s in ((0.0, 1.0), (0.7, 0.4), (-1.3, 0.5), (2.5, 0.3), (0.05, 1e-3)):              # Gaussians N(mu, s^2)
        sig = math.sqrt(mu * mu + s * s); t = sig * np.exp(js * h)
        q = np.array([integrate.quad(lambda z: z * z / (z * z + tt * tt) * math.exp(-0.5 * ((z - mu) / s) ** 2) / (s * math.sqrt(2 * math.pi)),
                                     mu - 12 * s, mu + 12 * s, points=[0.0] if abs(mu) < 12 * s else None, limit=400, epsabs=1e-14, epsrel=1e-12)[0] for tt in t])
        exact = s * (special.erfc(-mu / (s * math.sqrt(2))) / 2 * mu / s + math.exp(-mu * mu / (2 * s * s)) / math.sqrt(2 * math.pi))
        err = abs(T(mu, q, sig) - exact) / sig; worst = max(worst, err); laws.append(f"N({mu},{s}^2) {err:.1e}")
    for pts, w in (([0.0, 1.0], [0.5, 0.5]), ([-2.0, 0.0, 3.0], [0.2, 0.5, 0.3]), ([1e-6, 1.0], [0.9, 0.1])):   # atoms, incl. at 0
        pts, w = np.array(pts), np.array(w); sig = math.sqrt(w @ pts ** 2); t = sig * np.exp(js * h)
        q = np.array([w @ (pts ** 2 / (pts ** 2 + tt ** 2)) for tt in t]); exact = w @ np.maximum(pts, 0)
        err = abs(T(w @ pts, q, sig) - exact) / sig; worst = max(worst, err); laws.append(f"atoms {err:.1e}")
    nu = 3.0; sig = math.sqrt(nu / (nu - 2)); t = sig * np.exp(js * h)                         # Student t_3 (finite variance, heavy tails)
    pdf = lambda z: special.gamma(2) / (math.sqrt(nu * math.pi) * special.gamma(1.5)) * (1 + z * z / nu) ** -2
    q = np.array([2 * integrate.quad(lambda z: z * z / (z * z + tt * tt) * pdf(z), 0, np.inf, limit=400, epsabs=1e-14)[0] for tt in t])
    exact = integrate.quad(lambda z: z * pdf(z), 0, np.inf, epsabs=1e-14)[0]; err = abs(T(0.0, q, sig) - exact) / sig; worst = max(worst, err)
    bound = (2 / (1 - math.exp(-math.pi * math.sqrt(M))) + 2 / math.pi) * math.exp(-math.pi * math.sqrt(M))
    return (f"[B4] 25-query resolvent quadrature (M=12): worst |error|/sigma over 9 laws {worst:.2e} <= stated bound {bound:.2e}"
            f" (< 5e-5): {ok(worst <= bound and bound < 5e-5)}; per law: " + ", ".join(laws) + f", t_3 {err:.1e}")


def b4_bracket():
    fails = 0
    for _ in range(500):
        k = int(rng.integers(2, 12)); G = rng.standard_normal((k, k)); H = G @ G.T * rng.random()
        v = rng.standard_normal(k); v /= np.linalg.norm(v); s = math.exp(rng.uniform(-3, 2)); A = H + s * s * np.eye(k)
        q = v @ H @ np.linalg.solve(A, v); x = np.linalg.solve(A, v) + 0.3 * rng.standard_normal(k) * rng.random()
        J = 2 * v @ x - x @ A @ x; r = v - A @ x; lo, hi = 1 - s * s * J - r @ r, 1 - s * s * J
        fails += not (lo - 1e-12 <= q <= hi + 1e-12 and hi - lo <= r @ r + 1e-12)
    return f"[B4b] variational bracket q_s in [1 - s^2 J - |r|^2, 1 - s^2 J], width <= |r|^2, 500 random cases: {ok(fails == 0)}"


def b5_leakage():
    fails = 0
    for _ in range(300):
        N, p = int(rng.integers(6, 16)), int(rng.integers(1, 4)); G = rng.standard_normal((N, N)); A = G @ G.T / N
        P = np.zeros((N, N)); P[:p, :p] = np.eye(p)
        a = np.linalg.eigvalsh(A[:p, :p])[0]; eta = np.linalg.norm(A[p:, :p], 2); tau = a * rng.random()
        w, U = np.linalg.eigh(A); delta = 0.1 * rng.random()
        e = np.where(w < tau, rng.uniform(-1, 1, N), rng.uniform(-delta, delta, N)); E = U @ np.diag(e) @ U.T
        lhs = np.linalg.norm(E[:p, :p], 2); rhs = delta + (1 - delta) * min(1.0, eta ** 2 / (a - tau) ** 2)
        fails += lhs > rhs + 1e-12
    return f"[B5a] leakage bound |P e(A) P| <= delta + (1-delta) min(1, eta^2/(a-tau)^2), 300 random cases: {ok(fails == 0)}"


def b5_clock():
    # reversible chain on R + E where every hidden state has the same escape mass: A 1_E = a 1_E; one step p = 1/a is exact
    fails = 0
    for _ in range(200):
        nE, nR = int(rng.integers(2, 6)), int(rng.integers(1, 4)); a = rng.uniform(0.1, 0.9)
        B = rng.random((nE, nE)); B = 0.5 * (B + B.T); np.fill_diagonal(B, 0)
        B *= (1 - a) / B.sum(1).max() * rng.uniform(0.3, 1)                                       # hidden-to-hidden, rows <= 1-a
        piE = np.full(nE, 1.0 / (nE + nR)); KEE = B + np.diag((1 - a) - B.sum(1))                 # rows of K_EE sum to 1 - a
        KER = np.full((nE, nR), a / nR); piR = np.full(nR, 1.0 / (nE + nR))
        KRE = (piE[:, None] * KER).T / piR[:, None]; KRR = np.diag(1 - KRE.sum(1))                # detailed balance, stochastic
        if (KRR < 0).any(): continue
        f = rng.standard_normal((nE + nR, 3)); pi = np.concatenate([piR, piE]); mu = pi @ f
        lam = piR / piR.sum(); fR, fE = f[:nR], f[nR:]
        muh = (lam @ (fR + KRE @ fE / a)) / (lam @ (1 + KRE @ np.ones(nE) / a))
        fails += not np.allclose(muh, mu)
    return f"[B5b] constant escape profile A 1_E = a 1_E makes p = 1/a exact for every reward: {ok(fails == 0)}"


def ungated():
    # Stage 25 Theorem 2.1 and the exact error identity m_hat - m = sum_l Gamma_l (A_hat_l - A_l), on a 2-D input net
    # integrated exactly by angular quadrature (the means are E R * average over the circle; F is piecewise linear in angle)
    d, n, L = 2, 5, 3; Ws = [rng.standard_normal((n, d if l == 0 else n)) * math.sqrt(2 / (d if l == 0 else n)) for l in range(L)]
    th = (np.arange(2_000_000) + 0.5) * 2 * math.pi / 2_000_000; X = np.stack([np.cos(th), np.sin(th)]) * math.sqrt(math.pi / 2)  # E R = sqrt(pi/2)
    zs, h = [], X
    for W in Ws: z = W @ h; zs.append(z); h = np.maximum(z, 0)
    m = [np.zeros(d)] + [np.maximum(z, 0).mean(1) for z in zs]; A = [np.abs(z).mean(1) for z in zs]
    rec = [np.zeros(d)]
    for l in range(L): rec.append(0.5 * (Ws[l] @ rec[-1] + A[l]))
    from whest_closure import closure                                                               # defined below
    mh, Ah = closure(Ws)
    def G(l):
        P = np.eye(n)
        for k in range(L - 1, l, -1): P = P @ Ws[k]
        return 0.5 ** (L - l) * P
    ident = sum(G(l) @ (Ah[l] - A[l]) for l in range(L))
    return (f"[U] ungated recursion m_l = (W_l m_(l-1) + E|z_l|)/2: max dev {max(np.max(np.abs(rec[l + 1] - m[l + 1])) for l in range(L)):.1e};"
            f" closure error identity m_hat_L - m_L = sum_l Gamma_l (A_hat_l - A_l): max dev {np.max(np.abs((mh - m[-1]) - ident)):.1e}"
            f" (closure error itself {np.max(np.abs(mh - m[-1])):.1e}): "
            f"{ok(max(np.max(np.abs(rec[l + 1] - m[l + 1])) for l in range(L)) < 1e-12 and np.max(np.abs((mh - m[-1]) - ident)) < 1e-12)}")


if __name__ == "__main__":
    import sys, types
    mod = types.ModuleType("whest_closure")
    def closure(Ws):
        """Gaussian closure: means and E|z| per layer under the moment-matched Gaussian (exact relu moments, Hermite series)."""
        from scipy.special import ndtr
        mu_h = np.zeros(Ws[0].shape[1]); C_h = np.eye(Ws[0].shape[1]); Ah = []
        for W in Ws:
            mu = W @ mu_h; C = W @ C_h @ W.T; sd = np.sqrt(np.diag(C)); a = mu / sd
            ph = np.exp(-a * a / 2) / math.sqrt(2 * math.pi); Ph = ndtr(a)
            m = sd * (ph + a * Ph); Ah.append(sd * (2 * ph + a * (2 * Ph - 1)))
            # covariance of relu(z) by the Hermite series in the correlation (exact diagonal)
            R = C / np.outer(sd, sd); He = [np.ones_like(a), -a]
            for k in range(2, 40): He.append(-a * He[-1] - (k - 1) * He[-2])
            dk = [m / sd, Ph] + [ph * He[k - 2] for k in range(2, 41)]
            K = sum(np.outer(sd * dk[k], sd * dk[k]) * R ** k / math.factorial(k) for k in range(1, 41))
            np.fill_diagonal(K, sd * sd * ((a * a + 1) * Ph + a * ph) - m * m)
            mu_h, C_h = m, K
        return mu_h, Ah
    mod.closure = closure; sys.modules["whest_closure"] = mod
    for fn in (b1_split, b1_quartic, b3_ito, b4_quadrature, b4_bracket, b5_leakage, b5_clock, ungated): print(fn(), flush=True)


def gamma_norms(path):
    """Operator norms of the ungated transports Gamma_l = 2^-(L-l+1) W_L...W_(l+1) for the given weights, against the
    Fuss-Catalan edge for iid N(0, 2/n) weights: |Gamma| -> (1/2) 2^(-k/2) ((k+1)^(k+1) / k^k)^(1/2), k = L - l."""
    W = [w.astype(np.float64) for w in np.load(path)]; L = len(W); out = []
    for l in range(L):                                   # 0-based: Gamma for z_l
        k = L - 1 - l; v = rng.standard_normal(W[0].shape[0]); v /= np.linalg.norm(v)
        for _ in range(60):                               # power iteration on P^T P, P = W_(L-1) ... W_(l+1)
            u = v.copy()
            for j in range(l + 1, L): u = W[j] @ u
            for j in range(L - 1, l, -1): u = W[j].T @ u
            v = u / np.linalg.norm(u)
        u = v.copy()
        for j in range(l + 1, L): u = W[j] @ u
        out.append((k, 0.5 ** (k + 1) * np.linalg.norm(u), 0.5 * 2 ** (-k / 2) * math.sqrt((k + 1) ** (k + 1) / k ** k if k else 1.0)))
    return out


if __name__ == "__main__" and len(__import__("sys").argv) > 1:
    for k, g, fc in gamma_norms(__import__("sys").argv[1])[::-1]:
        print(f"[G] k = {k:2d}: |Gamma| = {g:.4f}  (Fuss-Catalan edge {fc:.4f})")
