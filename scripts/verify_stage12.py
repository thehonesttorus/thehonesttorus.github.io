"""Checks for the stage-12 note (the network as its own two-tier system).

  1. Self-localization on a random He network (n=1024, L=16, 2e5 inputs): energy fraction of h_l in its mean versus
     cos(theta_l); internal clock t_l = |m_l|^2/Tr Cov(h_l); gate variance versus theta_{l-1}/2pi; log-norm variance versus
     the correlated-increment prediction; field law of the next layer (variance, compression by the collective mode, fraction
     with |r|<1); collective mode (trace share along the mean versus t V/4, and the Frobenius share used in stage 9);
     backward clock |E delta_l|^2/E|delta_l|^2 versus prod (1-theta_{k-1}/pi); the marginal memory of each neuron
     (E relu(z_a) minus its Gaussian closure) against the Gram-Charlier wall terms, binned by |r|; the collective-mode law
     of the marginal cumulants (skewness linear in the field with slope 1.5 V, flat excess kurtosis 3 V).
  2. Lemma: E_{r~N(0,t)} Phi(r)Phi(-r) = arccos(t/(1+t))/2pi; the clock recursion t' = E m1^2/E s1^2 equals the angle map.
  3. Fisher-Rao distance of N(mu, sigma^2) to the kink geodesic {mu=0} is sqrt2 asinh(|r|/sqrt2).
  4. Stratification of the third derivative of a two-layer suffix (softplus, exact chain rule) against Hermite integration
     by parts: deep-wall layer + bent crossings + hub layer = E[d^3 S], with Monte Carlo standard errors.
  5. Birth law of third cumulants: E_{N(0,t)} kappa3(relu(r+xi))^2 ~ (int kappa3^2) / sqrt(2 pi t).
  6. Wall-density law, Hermite cancellation, saturation of the across-input log-norm variance, Ising couplings of the
     gate gas, and the power law of the backward clock.
"""
import numpy as np
from math import pi, acos, cos, sin, sqrt, asinh, acosh, erf, log
from scipy.stats import norm, multivariate_normal
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

phi, Phi = norm.pdf, norm.cdf
F = lambda th: acos((sin(th) + (pi - th) * cos(th)) / pi)
J2 = lambda t: 3 * sin(t) * cos(t) + (pi - t) * (1 + 2 * cos(t) ** 2)
m1 = lambda r: phi(r) + r * Phi(r)
m2 = lambda r: (1 + r * r) * Phi(r) + r * phi(r)
m3 = lambda r: (r ** 3 + 3 * r) * Phi(r) + (r * r + 2) * phi(r)

print("== 1. Self-localization of a He network (n=1024, L=16)")
rng = np.random.default_rng(5)
n, L, N, B = 1024, 16, 200000, 10000
Ws = [rng.normal(0, sqrt(2 / n), size=(n, n)).astype(np.float32) for _ in range(L)]
c_out = rng.normal(size=n).astype(np.float32)
th = [pi / 2]
for _ in range(L):
    th.append(F(th[-1]))
Sh = [np.zeros(n) for _ in range(L + 1)]
Sh2 = np.zeros(L + 1)
Sz = [np.zeros(n) for _ in range(L + 1)]
Sg = [np.zeros(n) for _ in range(L + 1)]
keep = (4, 8, 12, 16)
Szz = {l: np.zeros((n, n)) for l in keep}
Sz3 = {l: np.zeros(n) for l in keep}
Sz4 = {l: np.zeros(n) for l in keep}
Shv2 = {l: np.zeros(n) for l in keep}
Sd = [np.zeros(n) for _ in range(L + 1)]
Sd2 = np.zeros(L + 1)
logn = [[] for _ in range(L + 1)]
for it in range(N // B):
    h = rng.normal(size=(B, n)).astype(np.float32)
    Sh[0] += h.sum(0)
    Sh2[0] += (h.astype(np.float64) ** 2).sum()
    logn[0].append(np.log((h.astype(np.float64) ** 2).sum(1)))
    gates = np.empty((L + 1, B, n), dtype=bool)
    for l in range(1, L + 1):
        z = h @ Ws[l - 1].T
        gates[l] = z > 0
        Sz[l] += z.sum(0)
        Sg[l] += gates[l].sum(0)
        if l in keep:
            z64 = z.astype(np.float64)
            Szz[l] += z64.T @ z64
            Sz3[l] += (z64 ** 3).sum(0)
            Sz4[l] += (z64 ** 4).sum(0)
        h = np.maximum(z, 0)
        Sh[l] += h.sum(0)
        Sh2[l] += (h.astype(np.float64) ** 2).sum()
        if l in keep:
            Shv2[l] += (h.astype(np.float64) ** 2).sum(0)
        logn[l].append(np.log((h.astype(np.float64) ** 2).sum(1)))
    d = np.broadcast_to(c_out, (B, n)).astype(np.float32)
    for l in range(L, -1, -1):
        Sd[l] += d.sum(0)
        Sd2[l] += (d.astype(np.float64) ** 2).sum()
        if l >= 1:
            d = (d * gates[l]) @ Ws[l - 1]
Vpred = 2.0 / n
print("   l  rho meas/pred   t_l meas/pred    gate var meas/pred   Var log|h|^2 meas/pred")
V = {}
tmeas = {}
for l in range(L + 1):
    m = Sh[l] / N
    rho = m @ m / (Sh2[l] / N)
    tp = cos(th[l]) / (1 - cos(th[l])) if l > 0 else 0.0
    tmeas[l] = rho / (1 - rho)
    if l >= 1:
        Vpred += 4 / n * (1.5 - J2(th[l - 1]) / (2 * pi))
    V[l] = np.concatenate(logn[l]).var()
    gv = ""
    if l >= 1:
        p = Sg[l] / N
        gv = "%.4f / %.4f" % ((p * (1 - p)).mean(), th[l - 1] / (2 * pi))
    if l in (1, 2, 4, 8, 12, 16):
        print("  %2d  %.4f / %.4f   %6.2f / %6.2f   %s   %.4f / %.4f" % (l, rho, cos(th[l]), tmeas[l], tp, gv, V[l], Vpred))
print("   l  Var(r) meas | t_{l-1} meas | compressed pred | t_{l-1} inf.width   P(|r|<1) meas/pred"
      "   coll. share meas / t V/4   Frobenius share")
ZZ = np.random.default_rng(11).normal(size=2_000_000)
shares, fvar = {}, {}
for l in keep:
    mu = Sz[l] / N
    C = Szz[l] / N - np.outer(mu, mu)
    r = mu / np.sqrt(np.diag(C))
    tp = cos(th[l - 1]) / (1 - cos(th[l - 1]))
    w = mu / np.linalg.norm(mu)
    share = (w @ C @ w) / np.trace(C)
    shares[l], fvar[l] = share, r.var()
    tb, eps = tmeas[l - 1] / (1 - share), share / tmeas[l - 1]
    rho_ = sqrt(tb) * ZZ
    comp = np.mean(rho_ ** 2 / (1 + eps * rho_ ** 2))
    Coff = C - np.diag(np.diag(C))
    e = np.linalg.eigvalsh(Coff)
    e = e[np.argsort(-np.abs(e))]
    print("  %2d   %6.2f | %6.2f | %6.2f | %6.2f     %.3f / %.3f       %.3f / %.3f          %.3f"
          % (l, r.var(), tmeas[l - 1], comp, tp, (np.abs(r) < 1).mean(), erf(1 / sqrt(2 * tp)), share,
             tmeas[l - 1] * V[l - 1] / 4, e[0] ** 2 / np.sum(e ** 2)))
print("   backward clock |E delta_l|^2 / E|delta_l|^2 (readout c ~ N(0,I)) against prod_{k>l}(1-theta_{k-1}/pi)")
for l in (0, 2, 4, 8, 12, 15, 16):
    md = Sd[l] / N
    pred = np.prod([1 - th[k - 1] / pi for k in range(l + 1, L + 1)])
    print("  l=%2d  %.4f / %.4f" % (l, md @ md / (Sd2[l] / N), pred))
print("   marginal memory D_a = (E relu(z_a) - s_a m1(r_a))/s_a against Gram-Charlier"
      " -g1 He1 phi/6 + g2 He2 phi/24 + g1^2 He4 phi/72")
for l in (8, 12, 16):
    mu = Sz[l] / N
    s2 = np.diag(Szz[l]) / N - mu ** 2
    s = np.sqrt(s2)
    r = mu / s
    k3 = Sz3[l] / N - 3 * mu * (Szz[l].diagonal() / N) + 2 * mu ** 3
    m4c = Sz4[l] / N - 4 * mu * Sz3[l] / N + 6 * mu ** 2 * Szz[l].diagonal() / N - 3 * mu ** 4
    g1, g2 = k3 / s ** 3, m4c / s ** 4 - 3
    D = (Sh[l] / N - s * m1(r)) / s
    GC = (-g1 * r / 6 + g2 * (r * r - 1) / 24 + g1 ** 2 * (r ** 4 - 6 * r * r + 3) / 72) * phi(r)
    floor = np.sqrt((Shv2[l] / N - (Sh[l] / N) ** 2) / N) / s
    print("  l=%2d: mean |g1| = %.4f, mean |g2| = %.4f; corr(D, GC) = %.3f" % (l, np.abs(g1).mean(), np.abs(g2).mean(),
                                                                         np.corrcoef(D, GC)[0, 1]))
    A1 = np.vstack([r, np.ones(n)]).T
    (b1, a1), *_ = np.linalg.lstsq(A1, g1, rcond=None)
    R2 = 1 - np.sum((g1 - A1 @ [b1, a1]) ** 2) / np.sum((g1 - g1.mean()) ** 2)
    A2 = np.vstack([r * r, np.ones(n)]).T
    (b2, a2), *_ = np.linalg.lstsq(A2, g2, rcond=None)
    R22 = 1 - np.sum((g2 - A2 @ [b2, a2]) ** 2) / np.sum((g2 - g2.mean()) ** 2)
    print("     collective-mode law: gamma1 = %.4f r %+.4f (R2 %.3f), model 1.5 V_{l-1} = %.4f (ratio %.2f);"
          " gamma2 = %.4f r^2 %+.4f (R2 %.3f), model 3 V_{l-1} = %.4f (ratio %.2f)"
          % (b1, a1, R2, 1.5 * V[l - 1], b1 / (1.5 * V[l - 1]), b2, a2, R22, 3 * V[l - 1], a2 / (3 * V[l - 1])))
    eta = b1 / (1.5 * V[l - 1])
    tb, eps = tmeas[l - 1] / (1 - shares[l]), eta * shares[l] / tmeas[l - 1]
    rho_ = sqrt(tb) * ZZ
    print("     field compression with the measured efficiency eta=%.2f: %.2f (measured Var r %.2f)"
          % (eta, np.mean(rho_ ** 2 / (1 + eps * rho_ ** 2)), fvar[l]))
    for lo, hi in [(0, 1), (1, 2), (2, 3), (3, 99)]:
        sel = (np.abs(r) >= lo) & (np.abs(r) < hi)
        if sel.sum() < 5:
            continue
        print("     |r| in [%d,%2d): %4d neurons  rms D %.5f  rms GC %.5f  rms(D-GC) %.5f  MC floor %.5f"
              % (lo, hi, sel.sum(), np.sqrt(np.mean(D[sel] ** 2)), np.sqrt(np.mean(GC[sel] ** 2)),
                 np.sqrt(np.mean((D - GC)[sel] ** 2)), np.sqrt(np.mean(floor[sel] ** 2))))

print("== 2. Field-law lemma and the clock recursion")
for t in (0.5, 2.0, 13.0):
    lhs = quad(lambda r: Phi(r) * Phi(-r) * norm.pdf(r, 0, sqrt(t)), -np.inf, np.inf)[0]
    num = quad(lambda r: m1(r) ** 2 * norm.pdf(r, 0, sqrt(t)), -np.inf, np.inf)[0]
    den = quad(lambda r: (m2(r) - m1(r) ** 2) * norm.pdf(r, 0, sqrt(t)), -np.inf, np.inf)[0]
    c = cos(F(acos(t / (1 + t))))
    print("  t=%5.1f: E Phi(r)Phi(-r) = %.8f, arccos(t/(1+t))/2pi = %.8f ; t' closure = %.6f, angle map = %.6f"
          % (t, lhs, acos(t / (1 + t)) / (2 * pi), num / den, c / (1 - c)))
print("  single-neuron SNR map at r=0: %.4f" % (m1(0) / sqrt(m2(0) - m1(0) ** 2)))

print("== 3. Fisher-Rao distance to the kink geodesic")
dFR = lambda a, b: sqrt(2) * acosh(1 + ((a[0] - b[0]) ** 2 / 2 + (a[1] - b[1]) ** 2) / (2 * a[1] * b[1]))
for mu_, sg in [(0.5, 1.0), (3.0, 0.7), (-2.0, 0.3)]:
    res = minimize_scalar(lambda ls: dFR((mu_, sg), (0, np.exp(ls))), bounds=(-10, 10), method="bounded")
    print("  r=%6.3f: min distance %.6f ; sqrt2 asinh(|r|/sqrt2) = %.6f" % (mu_ / sg, res.fun, sqrt(2) * asinh(abs(mu_ / sg) / sqrt(2))))

print("== 4. Stratified third derivative of a two-layer suffix (softplus, exact chain rule)")
rng = np.random.default_rng(7)
mu3 = np.array([0.3, -0.2, 0.5])
A = rng.normal(size=(3, 3))
Sig = A @ A.T / 3 + 0.5 * np.eye(3)
Wv = np.array([0.9, -0.7, 0.6])
P = np.linalg.inv(Sig)
Lc = np.linalg.cholesky(Sig)
kap = rng.normal(size=(3, 3, 3))
kap = sum(kap.transpose(p) for p in [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]) / 6
sp = lambda u, b: np.logaddexp(0, b * u) / b
s1 = lambda u, b: 1 / (1 + np.exp(-b * u))
s2 = lambda u, b: b * s1(u, b) * (1 - s1(u, b))
s3 = lambda u, b: b * b * s1(u, b) * (1 - s1(u, b)) * (1 - 2 * s1(u, b))
for beta in (5.0, 20.0):
    Nn, Bt = 20_000_000, 2_000_000
    sums = np.zeros(5)
    sq = np.zeros(2)
    for it in range(Nn // Bt):
        z = mu3 + rng.normal(size=(Bt, 3)) @ Lc.T
        y = sp(z, beta) @ Wv
        S = sp(y, beta)
        v = (z - mu3) @ P
        kv = np.einsum("bcd,ib,ic,id->i", kap, v, v, v) - 3 * np.einsum("bcd,bc,id->i", kap, P, v)
        g = Wv * s1(z, beta)
        a = s3(y, beta) * np.einsum("bcd,ib,ic,id->i", kap, g, g, g)
        b_ = 3 * s2(y, beta) * np.einsum("bcc,ib,ic->i", kap, g, Wv * s2(z, beta))
        c_ = s1(y, beta) * np.einsum("bbb,ib->i", kap, Wv * s3(z, beta))
        sums += [(S * kv).sum(), a.sum(), b_.sum(), c_.sum(), (S * kv - a - b_ - c_).sum()]
        sq += [((S * kv) ** 2).sum(), ((S * kv - a - b_ - c_) ** 2).sum()]
    T, Aa, Bb, Cc, Dd = sums / Nn
    seT, seD = np.sqrt(sq / Nn - np.array([T, Dd]) ** 2) / sqrt(Nn)
    print("  beta=%4.0f: Hermite IBP %.4f (se %.4f) | deep-wall %.4f + crossings %.4f + hub %.4f = %.4f ; difference %.4f (se %.4f)"
          % (beta, T, seT, Aa, Bb, Cc, Aa + Bb + Cc, Dd, seD))

print("== 5. Birth law of third cumulants")
k3 = lambda r: m3(r) - 3 * m1(r) * m2(r) + 2 * m1(r) ** 3
I3 = quad(lambda r: k3(r) ** 2, -np.inf, np.inf)[0]
print("  c3 = %.6f = (pi+2)/(2pi)^1.5 = %.6f ; int kappa3^2 dr = %.5f" % (k3(0), (pi + 2) / (2 * pi) ** 1.5, I3))
for t in (1, 4, 13, 50):
    Bt_ = quad(lambda r: k3(r) ** 2 * norm.pdf(r, 0, sqrt(t)), -np.inf, np.inf)[0]
    print("  t=%3d: E kappa3^2 = %.5f ; asymptote %.5f" % (t, Bt_, I3 / sqrt(2 * pi * t)))

print("== 6. Wall-density law, Hermite cancellation, saturation, Ising couplings, backward power law")
I1 = quad(lambda r: Phi(r) * Phi(-r), -np.inf, np.inf)[0]
print("  int Phi(r)Phi(-r) dr = %.10f ; 1/sqrt(pi) = %.10f" % (I1, 1 / sqrt(pi)))
for t in (1.0, 13.0, 100.0):
    tq = acos(t / (1 + t))
    Ephi = quad(lambda r: phi(r) * norm.pdf(r, 0, sqrt(t)), -np.inf, np.inf)[0]
    EH2 = quad(lambda r: (r * r - 1) * phi(r) * norm.pdf(r, 0, sqrt(t)), -np.inf, np.inf)[0]
    print("  t=%5.1f: E phi = %.8f, sin(theta/2)/sqrt(pi) = %.8f ; gate var %.6f vs I1 sin(theta/2)/sqrt(pi) = %.6f ;"
          " E He2 phi = %.3e vs -(1+t)^-1.5/sqrt(2pi) = %.3e"
          % (t, Ephi, sin(tq / 2) / sqrt(pi), tq / (2 * pi), I1 * sin(tq / 2) / sqrt(pi), EH2, -(1 + t) ** -1.5 / sqrt(2 * pi)))
thl = [pi / 2]
for _ in range(100000):
    thl.append(F(thl[-1]))
inc = np.array([4 * (1.5 - J2(x) / (2 * pi)) for x in thl[:-1]])
print("  (3/2 - J2/2pi)/theta^2 at theta=0.1, 0.01: %.5f, %.6f" % ((1.5 - J2(0.1) / (2 * pi)) / 0.01, (1.5 - J2(0.01) / (2 * pi)) / 1e-4))
print("  n*V_l (across inputs): l=16 %.2f, l=100 %.2f, l=1e5 %.2f ; Hanin-Nica over (W,x) at l=16: %.0f" %
      (2 + inc[:16].sum(), 2 + inc[:100].sum(), 2 + inc.sum(), 2 + 5 * 16))
Vinf = 2 + inc.sum()
print("  crossover t V/4 = 1: tau* = 3pi sqrt(2n/Vinf) = %.3f sqrt(n)" % (3 * pi * sqrt(2 / Vinf)))
for ra, rb, cc in [(0.3, -0.8, 0.02), (1.5, 0.4, 0.01), (-1.0, 2.0, 0.03)]:
    P11 = multivariate_normal(mean=[0, 0], cov=[[1, cc], [cc, 1]]).cdf([ra, rb])
    pa, pb = Phi(ra), Phi(rb)
    J = log(P11 * (1 - pa - pb + P11) / ((pa - P11) * (pb - P11)))
    print("  Ising coupling r=(%.1f,%.1f), c=%.2f: log odds ratio %.5f ; c phi phi/(p q p q) %.5f"
          % (ra, rb, cc, J, cc * phi(ra) * phi(rb) / (pa * (1 - pa) * pb * (1 - pb))))
Phi0 = lambda x: 3 * pi / x + 1.5 * log(x)
for l0, l1 in [(100, 1000), (1000, 10000)]:
    prod = np.prod([1 - thl[k - 1] / pi for k in range(l0 + 1, l1 + 1)])
    print("  backward product l=%d..%d: %.5f ; (Phi0(theta_l0)/Phi0(theta_l1))^3 = %.5f" % (l0, l1, prod, (Phi0(thl[l0]) / Phi0(thl[l1])) ** 3))
