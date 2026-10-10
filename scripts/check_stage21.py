"""Exact (and small-n statistical) checks of Stage 21's statements, independent of the library's series where possible.
python scripts/check_stage21.py"""
import sys, os, itertools, numpy as np
from scipy import integrate, stats
from scipy.special import ndtr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest import gate_channel as gc
from whest.relu_gauss import relu_moments

rng = np.random.default_rng(21); FAIL = []


def ok(name, err, tol):
    good = bool(np.all(np.asarray(err) <= tol)); print(f"{'PASS' if good else 'FAIL'}  {name}: {np.max(err):.2e}")
    if not good: FAIL.append(name)


def eu2(m1, m2, s1, s2, r):
    """E[relu(z1) relu(z2)] by adaptive quadrature over the positive quadrant (independent of any series)."""
    cov = np.array([[s1 * s1, r * s1 * s2], [r * s1 * s2, s2 * s2]]); iv = np.linalg.inv(cov); dt = np.linalg.det(cov)
    f = lambda y, x: x * y * np.exp(-0.5 * np.array([x - m1, y - m2]) @ iv @ np.array([x - m1, y - m2])) / (2 * np.pi * np.sqrt(dt))
    return integrate.dblquad(f, 0, m1 + 12 * s1, 0, m2 + 12 * s2, epsabs=1e-13, epsrel=1e-12)[0]


# ---- Thm 2.1 (Price channel): d E[u_i u_j] / d Sigma_ij = Pr(z_i > 0, z_j > 0), d E[u_i^2]/d Sigma_ii = Pr(z_i > 0)
e = []
for _ in range(4):
    m1, m2 = rng.normal(size=2) * 0.7; s1, s2 = 0.6 + rng.random(2); r = rng.uniform(-0.6, 0.6); h = 1e-4
    c12 = r * s1 * s2
    up = eu2(m1, m2, s1, s2, (c12 + h) / (s1 * s2)); dn = eu2(m1, m2, s1, s2, (c12 - h) / (s1 * s2))
    orth = stats.multivariate_normal(mean=[0, 0], cov=[[s1 * s1, c12], [c12, s2 * s2]], abseps=1e-12, releps=1e-12).cdf([m1, m2])   # Pr(z1 > 0, z2 > 0)
    e.append(abs((up - dn) / (2 * h) - orth))
    g = lambda v: integrate.quad(lambda x: x * x * np.exp(-(x - m1) ** 2 / (2 * v)) / np.sqrt(2 * np.pi * v), 0, np.inf, epsabs=1e-14)[0]
    e.append(abs((g(s1 * s1 + h) - g(s1 * s1 - h)) / (2 * h) - ndtr(m1 / s1)))
ok("Thm 2.1: Price derivative = co-activation kernel (off-diagonal and diagonal), by quadrature and Genz", e, 2e-6)
# the same identity for the library's series kernels at a random n = 6 Gaussian state, as a matrix statement dS^u = K o dS^z
n = 6; A = rng.standard_normal((n, n)); C = A @ A.T / n + 0.3 * np.eye(n); mu = 0.5 * rng.standard_normal(n)
dX = rng.standard_normal((n, n)); dX = (dX + dX.T) / 2; h = 1e-5
Sp, _ = gc.second_moment(mu, C + h * dX, 30); Sm, _ = gc.second_moment(mu, C - h * dX, 30)
Kc = gc.cut_kernel(mu, C, 30); dS = (Sp - Sm) / (2 * h)
dX0 = dX - np.diag(np.diag(dX)); S0p, _ = gc.second_moment(mu, C + h * dX0, 30); S0m, _ = gc.second_moment(mu, C - h * dX0, 30)
ok("Thm 2.1, off-diagonal directions: dS^u = K o dS^z exactly", np.abs((S0p - S0m) / (2 * h) - Kc * dX0).max(), 1e-6)
print(f"      Thm 2.1, general direction: |dS^u - K o dS^z| = {np.abs(dS - Kc * dX).max():.2e} (the omitted wall coupling)")
ok("Thm 2.1 corrected: dS^u = K o dS^z + (diag(dS) psi + psi^T diag(dS))/2, psi_ij = p_i(0) E[u_j | z_i = 0]",
   np.abs(dS - gc.full_tangent(mu, C, Kc, dX)).max(), 1e-6)
mp, Kp, _, _ = relu_moments(mu, C + h * dX, 30); mm, Km, _, _ = relu_moments(mu, C - h * dX, 30); dK = (Kp - Km) / (2 * h)
off = ~np.eye(n, dtype=bool)
print(f"      centred covariance: off-diagonal |dKh - K o dSigma| = {np.abs(dK - Kc * dX)[off].max():.2e} vs raw S^u {np.abs(dS - Kc * dX)[off].max():.2e}"
      f" (the wall coupling minus d(m m^T) leaves only correlation-order terms)")

# ---- Prop 2.2: Kraus form over the 2^n cuts, Stinespring dual, level-2 marginal; level-3 marginal via Edgeworth/Stein
n = 4; A = rng.standard_normal((n, n)); C = A @ A.T / n + 0.2 * np.eye(n); mu = 0.4 * rng.standard_normal(n); W = rng.standard_normal((n, n)) * np.sqrt(2 / n)
cuts = list(itertools.product((0, 1), repeat=n)); pS = []
for S in cuts:                                         # Pr(sign pattern S) for z ~ N(mu, C): flip coordinates not in S
    sg = np.array([1 if b else -1 for b in S]); Cs = C * np.outer(sg, sg)
    pS.append(stats.multivariate_normal(mean=np.zeros(n), cov=Cs, maxpts=4_000_000, abseps=1e-11, releps=1e-11).cdf(sg * mu))
pS = np.array(pS); X = rng.standard_normal((n, n)); X = X @ X.T
Kraus = sum(p * np.diag(S) @ W @ X @ W.T @ np.diag(S) for p, S in zip(pS, cuts))
K2 = sum(p * np.outer(S, S) for p, S in zip(pS, cuts))
V = np.concatenate([np.sqrt(p) * np.diag(S) @ W for p, S in zip(pS, cuts)], 0)   # Stinespring: R^n -> R^n (x) l^2(cuts)
ok("Prop 2.2: sum_S p(S) T_S X T_S^T = K o (W X W^T), K = level-2 marginal; V^*V = F", [abs(pS.sum() - 1), np.abs(Kraus - K2 * (W @ X @ W.T)).max(),
   np.abs(V.T @ V - gc.channel_dual(W, np.diag(K2))).max()], 1e-6)
ok("Prop 2.2: tetrachoric series kernel = Genz orthant probabilities", np.abs(gc.cut_kernel(mu, C, 40) - K2).max(), 2e-6)
# level 3 at the Gaussian point (distinct indices): E[u1 u2 u3 H_123(z)] = Pr(z1, z2, z3 > 0), H the multivariate Hermite
n = 3; A = rng.standard_normal((n, n)); C = A @ A.T / n + 0.5 * np.eye(n); mu = 0.3 * rng.standard_normal(n); P = np.linalg.inv(C)
Z = rng.multivariate_normal(mu, C, size=8_000_000); Y = (Z - mu) @ P
H3 = Y[:, 0] * Y[:, 1] * Y[:, 2] - (Y[:, 0] * P[1, 2] + Y[:, 1] * P[0, 2] + Y[:, 2] * P[0, 1])
val = np.maximum(Z, 0).prod(1) * H3; orth3 = np.mean((Z > 0).all(1))
ok("Prop 2.2: level-3 tangent dE[u1u2u3]/dkappa_123 = Pr(1,2,3 in S) (z-score, 8e6 samples)", abs(val.mean() - orth3) / (val.std() / np.sqrt(len(val))), 5)

# ---- Prop 2.3: trace preservation is flatness; mirrored orthogonal layers are exactly trace preserving
n = 8; W = rng.standard_normal((n, n)) * np.sqrt(2 / n); Pa = rng.random(n); X = rng.standard_normal((n, n)); X = X @ X.T
Kc = np.outer(Pa, Pa); np.fill_diagonal(Kc, Pa)
ok("Prop 2.3: Tr Phi(X) = Tr(F X)", abs(np.trace(gc.channel(W, Kc, X)) - np.trace(gc.channel_dual(W, Pa) @ X)), 1e-12)
O = np.linalg.qr(rng.standard_normal((n, n)))[0]; Wm = np.concatenate([O, -O], 0); mu_h = rng.random(n)
al = (Wm @ mu_h) / 0.7
ok("Prop 2.3(b): mirrored orthogonal layer, F = I exactly", np.abs(gc.channel_dual(Wm, ndtr(al)) - np.eye(n)).max(), 1e-12)
# He criticality: E F on mu_hat^perp = Pi_perp (statistical, n = 400, 20 draws)
n = 400; vals = []
for _ in range(20):
    W = rng.standard_normal((n, n)) * np.sqrt(2 / n); mh = np.abs(rng.standard_normal(n)) + 0.2; s = np.sqrt(2 / n * 0.6 * n)
    Pa = ndtr((W @ mh) / s); vals.append(np.mean(gc.flatness(gc.channel_dual(W, Pa), mh)))
ok("Prop 2.3: E Tr_perp(F - I)/(n-1) = 0 at He (20 draws, n = 400; |mean| below 4 standard errors)", abs(np.mean(vals)) / (np.std(vals) / np.sqrt(20)) / 4, 1)

# ---- Thm 3.1: closed form, the trajectory table, neutrality of the scale mode
rs = np.linspace(0.01, 0.97, 40)
gq = [2 * integrate.quad(lambda x: ndtr(np.sqrt(r / (1 - r)) * x) ** 2 * np.exp(-x * x / 2) / np.sqrt(2 * np.pi), -np.inf, np.inf, epsabs=1e-13)[0] for r in rs]
ok("Thm 3.1(a): 2 E Phi(alpha)^2 = 1/2 + arcsin(rho)/pi = f'(rho) (adaptive quadrature; numerical derivative)",
   [np.abs(np.array(gq) - gc.gap_closed(rs)).max(),
    np.abs((gc.relu_corr_map(rs + 1e-6) - gc.relu_corr_map(rs - 1e-6)) / 2e-6 - gc.gap_closed(rs)).max()], 2e-8)
rho, g = gc.trajectory(16)
tab = [.500, .603, .664, .707, .738, .763, .783, .800, .814, .826, .837, .846, .854, .862, .868, .874]
ok("Thm 3.1: trajectory table g_1..g_16 (to the printed 3 decimals)", np.abs(g - np.array(tab)).max(), 5e-4)
loc = {(1, .1): None, (2, .1): 12, (2, .05): 14, (2, .02): 16, (2, .01): None, (3, .1): 9, (3, .05): 11, (3, .02): 13, (3, .01): 14,
       (4, .1): 8, (4, .05): 9, (4, .02): 11, (4, .01): 12}
ok("Cor 3.2: localisation-radius table", [0 if gc.radii(g, k, eta) == v else 1 for (k, eta), v in loc.items()], 0)
n = 50; A = rng.standard_normal((n, n)); C = A @ A.T / n + 0.2 * np.eye(n); mu = rng.standard_normal(n) * 0.5; lam = 1.7
m1, K1, _, _ = relu_moments(mu, C); m2, K2_, _, _ = relu_moments(lam * mu, lam * lam * C)
ok("Thm 3.1(c): scale mode is neutral (homogeneity m(l mu, l^2 C) = l m, Kh -> l^2 Kh)", [np.abs(m2 - lam * m1).max(), np.abs(K2_ - lam * lam * K1).max()], 1e-12)
# Thm 3.1(b) statistical check at n = 1024 with Gaussian weights and trajectory alphas (the note's own MC, larger n)
n = 1024; r0 = 0.8; gg = gc.gap_closed(r0); out = {1: [], 2: [], 3: []}
for _ in range(40):
    W = rng.standard_normal((n, n)) * np.sqrt(2 / n); Pa = ndtr(rng.normal(0, np.sqrt(r0 / (1 - r0)), n)); v = rng.standard_normal((n, 4))
    for k in out: out[k] += list(gc.level_gain(W, Pa, v, k))
print("      Thm 3.1(b), n=1024, rho=0.8: " + ", ".join(f"k={k}: {np.mean(v):.4f} vs g^k {gg ** k:.4f}" for k, v in out.items()))
ok("Thm 3.1(b): level-k contraction at n = 1024 within 2.5% of g^k (k <= 3)", [abs(np.mean(v) / gg ** k - 1) for k, v in out.items()], 0.025)

# ---- Prop 4.1: Szegedy phases
t = 0.01; Rm = np.diag(np.exp(-np.arange(4) * t / 2)); Pm = Rm @ Rm.T
# two-projection geometry: projections onto span{e_i} and onto span{cos th e_i + sin th f_i} meeting at angle arccos(sigma)
th = np.arccos(np.diag(Rm)); dim = 8; Pp = np.zeros((dim, dim)); Qp = np.zeros((dim, dim))
for i, a in enumerate(th):
    e1 = np.zeros(dim); e1[2 * i] = 1; f1 = np.zeros(dim); f1[2 * i] = np.cos(a); f1[2 * i + 1] = np.sin(a)
    Pp += np.outer(e1, e1); Qp += np.outer(f1, f1)
U = (2 * Pp - np.eye(dim)) @ (2 * Qp - np.eye(dim)); ph = np.sort(np.abs(np.angle(np.linalg.eigvals(U))))
ok("Prop 4.1: reflection-walk phases +-2 arccos(e^(-rt/2)); level-1 gap 2 arccos e^(-t/2) ~ 2 sqrt t",
   [np.abs(ph - np.sort(np.repeat(gc.szegedy_phases(np.arange(4), t), 2))).max(), abs(gc.szegedy_phases(1, t) / (2 * np.sqrt(t)) - 1) - 2 * t], 1e-9)

# ---- Prop 4.2: lifting a spectral approximation to Sym^r; per-neuron bound
d = 4; A = rng.standard_normal((d, d)); Cm = A @ A.T + 0.3 * np.eye(d); eps = 0.05
B = rng.standard_normal((d, d)); B = (B + B.T) / 2; Cr = np.linalg.cholesky(Cm); Ev = np.linalg.eigvalsh(Cr.T @ np.linalg.inv(Cm) @ Cr)
Mid = Cr @ (np.eye(d) + eps * B / np.abs(np.linalg.eigvalsh(B)).max()) @ Cr.T          # (1-eps) C <= Ct <= (1+eps) C
e = []
for r in (2, 3):
    Cr_ = Cm; Ct_ = Mid
    for _ in range(r - 1): Cr_ = np.kron(Cr_, Cm); Ct_ = np.kron(Ct_, Mid)
    lo = np.linalg.eigvalsh(Ct_ - (1 - eps) ** r * Cr_).min(); hi = np.linalg.eigvalsh((1 + eps) ** r * Cr_ - Ct_).min()
    e += [max(0, -lo), max(0, -hi)]
ok("Prop 4.2(a): (1-eps)^r C^(x r) <= Ct^(x r) <= (1+eps)^r C^(x r)", e, 1e-9)
from scripts.kik_merge import gauss_pos
mk = rng.normal(size=200); wk = rng.standard_normal((200, d)); q0 = np.einsum("ki,ij,kj->k", wk, Cm, wk); q1 = np.einsum("ki,ij,kj->k", wk, Mid, wk)
sk = np.sqrt(q0); bound = 0.5 * eps * sk * np.exp(-0.5 * (mk / sk) ** 2) / np.sqrt(2 * np.pi)
G = lambda m_, q_: gauss_pos(m_, q_ + m_ * m_)
ok("Prop 4.2(b): |G(m, w^T Ct w) - G(m, w^T C w)| <= eps sigma phi(alpha)/2 (1 + O(eps))", np.max(np.abs(G(mk, q1) - G(mk, q0)) / bound) - 1, 2 * eps)
# (G(m, tau^2) = E relu(N(m, tau^2)); the bound uses d G / d tau^2 = phi(alpha) / (2 tau))

# ---- Thm 4.3: orthogonal-probe trace estimation
n = 30; A = rng.standard_normal((n, n)); A = (A + A.T) / 2; m = 5
est = np.array([gc.orth_probe_trace(A, m, rng) for _ in range(40000)])
ok("Thm 4.3: unbiased, exact variance formula (z-scores), exact at m = n",
   [abs(est.mean() - np.trace(A)) / (est.std() / 200), abs(est.var() / gc.orth_probe_var(A, m) - 1) / (np.sqrt(2 / 40000) * 3) / 5,
    abs(gc.orth_probe_trace(A, n, rng) - np.trace(A)) * 1e9], 5)

# ---- Prop 5.1: Wick-orthogonality test on a synthetic i.i.d.-row readout y_k = G(w_k.mu, w_k^T C w_k)
n = 64; A = rng.standard_normal((n, n)); C = A @ A.T / n; mu = np.abs(rng.normal(size=n)); mh = mu / np.linalg.norm(mu)
# Ybar(a) = E[y | a]: a degree-8 polynomial regression fitted on an independent large sample of rows
Wb = rng.standard_normal((400000, n)) * np.sqrt(2 / n); ab = Wb @ mh; yb = G(Wb @ mu, np.einsum("ki,ij,kj->k", Wb, C, Wb))
Vb = np.vander(ab / ab.std(), 9); coef = np.linalg.lstsq(Vb, yb, rcond=None)[0]
Tst = []
for _ in range(200):
    Wf = rng.standard_normal((1024, n)) * np.sqrt(2 / n); a = Wf @ mh; y = G(Wf @ mu, np.einsum("ki,ij,kj->k", Wf, C, Wf))
    eta = y - np.vander(a / ab.std(), 9) @ coef; Tst.append([np.mean(eta), np.mean(eta * a / ab.std()), np.mean(eta * (a * a / ab.var() - 1))])
Tst = np.array(Tst); z = np.abs(Tst.mean(0)) / (Tst.std(0) / np.sqrt(len(Tst)))
ok("Prop 5.1: T_psi has mean zero for psi = 1, He1, He2 (z-scores over 200 draws of 1024 i.i.d. rows)", z, 4)

print(f"\n{'ALL PASS' if not FAIL else 'FAILED: ' + ', '.join(FAIL)}"); sys.exit(1 if FAIL else 0)
