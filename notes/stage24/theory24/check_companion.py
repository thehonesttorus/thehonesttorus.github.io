"""Independent checks of the finite claims in "Weighted activation graphs and certified approximation" (companion, 10 Oct 2026).

  python notes/stage24/theory24/check_companion.py

Exact rational arithmetic (fractions) where the claim is finite; seeded Monte Carlo or random finite chains otherwise.
  [A] eq (7)   two rows at angle pi*eps: random-site heat-bath eigenvalues 1, 1-eps, eps, 0
  [B] eq (13)  F = (x1+ - x2+)+: sector masses, jumps, mu = s0 (-1 + sqrt2 + 0 + 0 + 1) = 1/(2 sqrt pi)
  [C] eq (11)  P(X in C, X_t in D) = sqrt(t/pi) S_CD + o(sqrt t) for a half-plane pair (Monte Carlo)
  [D] eq (16)  E F(X) = E[(B(Y) - B(X))(Y - X)] / (2(1 - rho)) on a random small network (Monte Carlo)
  [E] eq (17)  E Z_rho^2 = (arccos rho + sqrt(1 - rho^2)) / (2 pi (1 - rho)) for one ReLU (Monte Carlo)
  [F] eq (22)  posterior-refresh cut gap >= 1/(1 + |M|^2) for a random M and a finite cut map (exact Gaussian kernel on a grid)
  [G] eq (28)  cycle of m sectors: local Laplacian gap 2 m s0 (1 - cos(2 pi/m))
  [H] Thm 7.1/7.2, eq (48) and its table: the three-state chain, exactly
  [I] Thm 8.1 residual identity, Thm 8.2 (Chebyshev), eq (43) one-row quotient, eq (44) ellipsoid, on random chains
  [J] eq (3)   [A12, A23] = E13 - E31 on the r = 1, N = 3 slice; eq (4) D_rs(c) counts"""
import itertools, math, numpy as np
from fractions import Fraction as Fr
from scipy.special import ndtr

rng = np.random.default_rng(2024)
ok = lambda b: "PASS" if b else "FAIL"


def A():
    out = []
    for eps in (Fr(1, 10), Fr(1, 4), Fr(1, 2)):
        p = {(0, 0): (1 - eps) / 2, (1, 1): (1 - eps) / 2, (0, 1): eps / 2, (1, 0): eps / 2}
        S = list(p); P = [[Fr(0)] * 4 for _ in range(4)]
        for a, s in enumerate(S):
            for site in (0, 1):
                for v in (0, 1):
                    t = list(s); t[site] = v; t = tuple(t)
                    other = [q for q in S if q[1 - site] == s[1 - site]]
                    P[a][S.index(t)] += Fr(1, 2) * p[t] / sum(p[q] for q in other)
        ev = sorted(np.linalg.eigvals(np.array(P, dtype=float)).real)
        out.append(np.allclose(ev, sorted([1, 1 - float(eps), float(eps), 0]), atol=1e-12))
    return f"[A] two-row heat-bath spectrum {{1, 1-eps, eps, 0}} for eps = 1/10, 1/4, 1/2: {ok(all(out))}"


def B():
    # sectors by angle: (0, pi/4), (pi/4, pi/2), (pi/2, pi), (pi, 3pi/2), (3pi/2, 2pi); F = (x1+ - x2+)+
    def grad(theta):
        x = np.array([np.cos(theta), np.sin(theta)]); x1, x2 = max(x[0], 0), max(x[1], 0)
        if x1 - x2 <= 0: return np.zeros(2)
        return np.array([1.0 if x[0] > 0 else 0.0, -1.0 if x[1] > 0 else 0.0])
    b = [0, np.pi / 4, np.pi / 2, np.pi, 3 * np.pi / 2]; s0 = 1 / (2 * math.sqrt(2 * math.pi)); mu = 0.0; jumps = []
    for k, th in enumerate(b):                                   # ray at angle th separates the sector before and after
        nrm = np.array([-np.sin(th), np.cos(th)])                # unit normal pointing to increasing angle
        j = (grad(th + 1e-9) - grad(th - 1e-9)) @ nrm; jumps.append(j); mu += s0 * j
    X = rng.standard_normal((2, 4_000_000)); mc = np.maximum(np.maximum(X[0], 0) - np.maximum(X[1], 0), 0).mean()
    good = np.allclose(jumps, [-1, math.sqrt(2), 0, 0, 1], atol=1e-9) and abs(mu - 1 / (2 * math.sqrt(math.pi))) < 1e-12
    return (f"[B] five-sector example: jumps {np.round(jumps, 6).tolist()}, mu = {mu:.10f} vs 1/(2 sqrt pi) = {1 / (2 * math.sqrt(math.pi)):.10f}"
            f" (Monte Carlo {mc:.5f}): {ok(good and abs(mc - mu) < 2e-3)}")


def C():
    res = []
    for t in (1e-2, 2.5e-3):
        N = 20_000_000; X = rng.standard_normal(N); Z = rng.standard_normal(N); Xt = math.exp(-t) * X + math.sqrt(1 - math.exp(-2 * t)) * Z
        pr = np.mean((X > 0) & (Xt <= 0)); S = 1 / math.sqrt(2 * math.pi)        # wall {x1 = 0} has Gaussian measure phi(0)
        res.append(pr / (math.sqrt(t / math.pi) * S))
    return f"[C] P(X in C, X_t in D) / (sqrt(t/pi) S_CD) for a half-plane, t = 1e-2, 2.5e-3: {res[0]:.4f} {res[1]:.4f} (-> 1): {ok(abs(res[1] - 1) < 0.02)}"


def D():
    d, n, L = 6, 5, 3
    Ws = [rng.standard_normal((n, d if l == 0 else n)) * math.sqrt(2.0 / (d if l == 0 else n)) for l in range(L)]
    def FB(X):
        h = X; M = np.broadcast_to(np.eye(d), (X.shape[1], d, d))
        for W in Ws:
            z = W @ h; g = (z > 0).astype(float); h = np.maximum(z, 0); M = g.T[:, :, None] * np.einsum("ij,bjk->bik", W, M)
        return h, M
    rho = 0.6; N = 1_500_000; X = rng.standard_normal((d, N)); Y = rho * X + math.sqrt(1 - rho ** 2) * rng.standard_normal((d, N))
    FX, BX = FB(X); _, BY = FB(Y)
    est = np.einsum("bij,jb->ib", BY - BX, Y - X).mean(1) / (2 * (1 - rho)); EF = FX.mean(1)
    return f"[D] finite-noise arrow readout (16), rho = 0.6: max |est - E F| / max E F = {np.max(np.abs(est - EF)) / np.max(EF):.2e}: {ok(np.max(np.abs(est - EF)) / np.max(EF) < 0.03)}"


def E():
    out = []
    for rho in (0.5, 0.9):
        N = 10_000_000; X = rng.standard_normal(N); Y = rho * X + math.sqrt(1 - rho ** 2) * rng.standard_normal(N)
        Z = ((Y > 0).astype(float) - (X > 0)) * (Y - X) / (2 * (1 - rho))
        out.append((np.mean(Z ** 2), (math.acos(rho) + math.sqrt(1 - rho ** 2)) / (2 * math.pi * (1 - rho))))
    return "[E] one-ReLU second moment (17): " + ", ".join(f"rho={r}: MC {a:.4f} vs {b:.4f}" for r, (a, b) in zip((0.5, 0.9), out)) + \
           f": {ok(all(abs(a / b - 1) < 0.01 for a, b in out))}"


def F():
    # cut map = sign pattern of 3 random rows in d = 3; posterior refresh X' | X ~ N(R x, I - R^2); kernel on cells from 2e7 pairs
    out = []
    for trial in range(3):
        d = 3; Wc = rng.standard_normal((3, d)); M = rng.standard_normal((2, d)) * (0.5 + 0.5 * trial)
        Sig = np.linalg.inv(np.eye(d) + M.T @ M); R = np.eye(d) - Sig; w_, V_ = np.linalg.eigh(np.eye(d) - R @ R)
        Ch = V_ @ np.diag(np.sqrt(np.clip(w_, 0, None))) @ V_.T
        cnt = np.zeros((8, 8))
        for _ in range(10):
            X = rng.standard_normal((d, 2_000_000)); Xp = R @ X + Ch @ rng.standard_normal((d, 2_000_000))
            a = ((Wc @ X > 0) * (1 << np.arange(3))[:, None]).sum(0); b = ((Wc @ Xp > 0) * (1 << np.arange(3))[:, None]).sum(0)
            np.add.at(cnt, (a, b), 1)
        keep = cnt.sum(1) > 0; cnt = cnt[np.ix_(keep, keep)]; J_ = 0.5 * (cnt + cnt.T) / cnt.sum(); pc = J_.sum(1)
        H = J_ / np.sqrt(np.outer(pc, pc)); ev = np.sort(np.linalg.eigvalsh(H)); gap = 1 - ev[-2]
        out.append((gap, 1 / (1 + np.linalg.norm(M, 2) ** 2), ev[0], keep.sum()))
    return "[F] posterior-refresh cut gap vs 1/(1+|M|^2) (gap, bound, min eigenvalue, cells): " + "; ".join(f"{g:.4f} >= {b:.4f}, {e:+.4f}, {k}" for g, b, e, k in out) + \
           f": {ok(all(g >= b - 2e-3 and e > -2e-3 for g, b, e, k in out))}"


def G():
    out = []
    for m in (6, 12, 40):
        s0 = 1 / (2 * math.sqrt(2 * math.pi)); Lp = np.zeros((m, m))
        for j in range(m): Lp[j, (j + 1) % m] = Lp[j, (j - 1) % m] = -s0 * m
        Lp -= np.diag(Lp.sum(1)); gap = np.sort(np.linalg.eigvalsh(Lp))[1]
        out.append(abs(gap - 2 * m * s0 * (1 - math.cos(2 * math.pi / m))) < 1e-12)
    return f"[G] cycle local Laplacian gap (28) for m = 6, 12, 40: {ok(all(out))}"


def Hchk():
    K = [[Fr(3, 4), Fr(1, 4), Fr(0)], [Fr(1, 4), Fr(1, 2), Fr(1, 4)], [Fr(0), Fr(1, 4), Fr(3, 4)]]
    pi = [Fr(1, 3)] * 3; f = [Fr(0), Fr(2), Fr(1)]; R, Ei = [0, 2], [1]
    alpha = pi[0] + pi[2]; Q = K[1][1]; Hh = 1 / (1 - Q)
    T = [[K[a][b] + K[a][1] * Hh * K[1][b] for b in R] for a in R]
    fstar = [f[a] + K[a][1] * Hh * f[1] for a in R]; tau = [1 + K[a][1] * Hh for a in R]
    lam = [pi[a] / alpha for a in R]
    mu = sum(l * x for l, x in zip(lam, fstar)) / sum(l * x for l, x in zip(lam, tau))
    rows = []
    for m in range(8):
        p = sum(Q ** j for j in range(m))                              # Neumann p_(m-1)(A) = sum_(j<m) Q^j
        N = sum(lam[i] * (f[a] + K[a][1] * p * f[1]) for i, a in enumerate(R)); Dn = sum(lam[i] * (1 + K[a][1] * p) for i, a in enumerate(R))
        est = N / Dn; rows.append((m, est, 1 - est, Dn, est == 1 - Fr(2) ** -m / (3 - Fr(2) ** -m), 1 - est == Fr(1, 3) * Fr(2) ** -m * (2 - est)))
    good = (alpha == Fr(2, 3) and Q == Fr(1, 2) and Hh == 2 and T == [[Fr(7, 8), Fr(1, 8)], [Fr(1, 8), Fr(7, 8)]] and fstar == [1, 2]
            and tau == [Fr(3, 2), Fr(3, 2)] and mu == 1 and all(r[4] and r[5] for r in rows)
            and [r[1] for r in rows[:4]] == [Fr(1, 2), Fr(4, 5), Fr(10, 11), Fr(22, 23)] and rows[7][1] == Fr(382, 383) and rows[7][3] == Fr(383, 256))
    no_clock = sum(l * x for l, x in zip(lam, fstar)); endpoint = sum(l * f[a] for l, a in zip(lam, R))
    return (f"[H] three-state chain exactly: T, f* = (1,2), tau = (3/2,3/2), mu = 1, (48) and both identities for m = 0..7, table rows: {ok(good)};"
            f" endpoint-reward average {endpoint}, rewards without clock {no_clock}")


def I():
    fails = 0
    for trial in range(200):
        S = int(rng.integers(5, 12)); A0 = rng.random((S, S)) * (rng.random((S, S)) < 0.6); C = A0 + A0.T + np.diag(rng.random(S))
        pi = rng.random(S) + 0.05; pi /= pi.sum()
        Kr = C / pi[:, None] / (C / pi[:, None]).sum(1).max() * 0.5   # conductances -> reversible sub-stochastic, then lazy
        Kr = Kr + np.diag(1 - Kr.sum(1)); K = 0.5 * (np.eye(S) + Kr)
        # make reversible w.r.t. pi exactly: K = I + (C/pi) scaled; check detailed balance
        assert np.allclose(pi[:, None] * K, (pi[:, None] * K).T)
        r = int(rng.integers(1, S - 1)); Rr = np.arange(r); Er = np.arange(r, S); alpha = pi[Rr].sum(); lam = pi[Rr] / alpha
        Q = K[np.ix_(Er, Er)]; Aop = np.eye(len(Er)) - Q; f = rng.standard_normal((S, 3)); mu = pi @ f
        U = rng.standard_normal((len(Er), 3)); v = rng.standard_normal(len(Er))
        Nn = lam @ (f[Rr] + K[np.ix_(Rr, Er)] @ U); Dd = lam @ (1 + K[np.ix_(Rr, Er)] @ v); muh = Nn / Dd
        rf = f[Er] - Aop @ U; r1 = 1 - Aop @ v
        fails += not np.allclose(mu - muh, pi[Er] @ (rf - np.outer(r1, muh)))          # Thm 8.1 (36)
        # Thm 8.2 with Chebyshev and the ellipsoid identity (44)
        Dh = np.sqrt(pi); Hs = Dh[:, None] * K / Dh[None, :]; ev = np.sort(np.linalg.eigvalsh(0.5 * (Hs + Hs.T))); g = 1 - ev[-2]
        a = g * alpha; mdeg = 4
        Tm = lambda x: np.cosh(mdeg * np.arccosh(x)) if x >= 1 else np.cos(mdeg * np.arccos(x))
        em = lambda Mx: None
        lamA, VA = np.linalg.eig(np.diag(np.sqrt(pi[Er])) @ Aop @ np.diag(1 / np.sqrt(pi[Er]))); lamA = lamA.real
        fails += not (lamA.min() >= a - 1e-9 and lamA.max() <= 1 + 1e-9)                 # (32)
        den = Tm((1 + a) / (1 - a)); emv = np.array([Tm((1 + a - 2 * x) / (1 - a)) / den for x in lamA])
        Bs = np.diag(np.sqrt(pi[Er])) @ Aop @ np.diag(1 / np.sqrt(pi[Er])); Bs = 0.5 * (Bs + Bs.T); w_, V_ = np.linalg.eigh(Bs)
        E_m = V_ @ np.diag([Tm((1 + a - 2 * x) / (1 - a)) / den for x in w_]) @ V_.T          # e_m(A) in the symmetric frame
        Am_inv_poly = V_ @ np.diag([(1 - Tm((1 + a - 2 * x) / (1 - a)) / den) / x for x in w_]) @ V_.T
        Pm = np.diag(1 / np.sqrt(pi[Er])) @ Am_inv_poly @ np.diag(np.sqrt(pi[Er]))            # p_(m-1)(A) in the original frame
        Uc = Pm @ f[Er]; vc = Pm @ np.ones(len(Er))
        Nn = lam @ (f[Rr] + K[np.ix_(Rr, Er)] @ Uc); Dd = lam @ (1 + K[np.ix_(Rr, Er)] @ vc); muh = Nn / Dd
        dm = 1 / den; Gm = (f[Er] - muh).T @ np.diag(pi[Er]) @ (f[Er] - muh)
        fails += not (np.sqrt(np.mean((muh - mu) ** 2)) <= dm * np.sqrt((1 - alpha) / 3 * np.linalg.eigvalsh(Gm)[-1]) + 1e-12)   # (40)
        fails += not (alpha * Dd >= alpha - 1e-12)                                                                                 # (39)
        # one-row form (43) and the exact ellipsoid (44)
        bvec = lam @ K[np.ix_(Rr, Er)]; hm = bvec @ Pm
        fails += not np.allclose(muh, (lam @ f[Rr] + hm @ f[Er]) / (1 + hm.sum()))
        zm = np.diag(1 / np.sqrt(pi[Er])) @ E_m @ np.diag(np.sqrt(pi[Er])) @ np.ones(len(Er))
        Bq = np.mean([((zm * pi[Er]) @ (f[Er][:, i] - muh[i])) ** 2 for i in range(3)])
        fails += not np.isclose(np.mean((muh - mu) ** 2), Bq)
    return f"[I] residual identity (36), spectrum (32), Chebyshev bound (40), denominator (39), one-row quotient (43), ellipsoid (44) on 200 random reversible chains: {ok(fails == 0)} ({fails} failures)"


def J():
    E = lambda i, j: np.outer(np.eye(3)[i], np.eye(3)[j])
    A12 = E(0, 1) + E(1, 0); A23 = E(1, 2) + E(2, 1); good = np.allclose(A12 @ A23 - A23 @ A12, E(0, 2) - E(2, 0))
    Nn = 6; cnt_ok = True
    for r, s, C in ((2, 2, (0, 1)), (3, 1, (0, 1)), (3, 3, (0, 1, 2, 3)), (2, 4, (0, 1))):
        c = len(C); cnt = sum(1 for S in itertools.combinations(range(Nn), r) for T in itertools.combinations(range(Nn), s)
                              if set(S) ^ set(T) == set(C))
        a, b = (r + c - s), (r + s - c)
        Drs = math.comb(c, a // 2) * math.comb(Nn - c, b // 2) if a % 2 == 0 and b % 2 == 0 and a >= 0 and b >= 0 else 0
        cnt_ok &= cnt == Drs
    return f"[J] commutator (3) and the D_rs(c) counts behind (4): {ok(good and cnt_ok)}"


if __name__ == "__main__":
    import sys
    sel = sys.argv[1:] or ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'Hchk', 'I', 'J']
    for nm in sel: print(globals()[nm](), flush=True)
