"""Checks of whest/kikuchi.py against every finite statement and worked example of the Kikuchi note.
python scripts/check_kikuchi.py      (prints one line per check; exits non-zero on any failure)"""
import sys, os, itertools, numpy as np
from scipy import integrate
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest import kikuchi as kk

rng = np.random.default_rng(20261010); FAIL = []


def ok(name, err, tol=1e-9):
    good = bool(np.all(np.asarray(err) <= tol)); print(f"{'PASS' if good else 'FAIL'}  {name}: {np.max(err):.2e}")
    if not good: FAIL.append(name)


def law(natoms=7):
    p = rng.random(natoms); return p / p.sum(), rng.standard_normal(natoms)


# ------------------------------------------------------------------ Sec 2: hard-core up/down, Kikuchi Laplacian, modular flow
N = 5; dim = 1 << N
U = np.zeros((dim, dim))
for S in range(dim):
    for i in range(N):
        if not (S >> i) & 1: U[S | 1 << i, S] = 1
D = U.T; comm = D @ U - U @ D; errs = []
for S in range(dim):
    errs.append(abs(comm[S, S] - (N - 2 * bin(S).count("1")))); errs.append(np.abs(np.delete(comm[S], S)).max())
ok("[D,U] = (N - 2k) I on each sector (eq. 1)", errs)
wed = {(i, j): rng.random() for i in range(N) for j in range(i + 1, N)}; LG = kk.kikuchi_laplacian(N, wed)
psi = rng.standard_normal(dim); form = 0.0
for (i, j), w in wed.items():
    for Zs in range(dim):
        if (Zs >> i) & 1 or (Zs >> j) & 1: continue
        form += w * (psi[Zs | 1 << i] - psi[Zs | 1 << j]) ** 2
ok("L_G >= 0 and its quadratic form (eq. 2)", [max(0, -np.linalg.eigvalsh(LG).min()), abs(psi @ LG @ psi - form)])
pS = rng.random(8); pS /= pS.sum(); a_, b_ = rng.standard_normal((2, 8, 8)); tt = 0.37
sig = lambda X: np.diag(pS ** (1j * tt)) @ X @ np.diag(pS ** (-1j * tt))
ok("modular flow: multiplicative and state-invariant (Sec. 2.3)",
   [np.abs(sig(a_ @ b_) - sig(a_) @ sig(b_)).max(), abs(np.trace(np.diag(pS) @ sig(a_)) - np.trace(np.diag(pS) @ a_))])

# ------------------------------------------------------------------ Sec 3: conditional gate tower
Ng = 4; cfgs = list(itertools.product((0, 1), repeat=Ng)); pr = rng.random(len(cfgs)); pr /= pr.sum()
prob = dict(zip(cfgs, pr)); Fv = dict(zip(cfgs, rng.standard_normal(len(cfgs))))
for k in (0, 1, 2):
    r = kk.gate_tower(prob, Fv, k)
    ok(f"Thm 3.1 at k={k}: D = U*, D f_(k+1) = f_k, norm identity", [r["adjoint_err"], r["tower_err"], r["norm_identity_err"]], 1e-12)
# Example 3.2: signs with 3/8 on equal, 1/8 on unequal
st = [(1, 1), (1, -1), (-1, 1), (-1, -1)]; P2 = np.array([3, 1, 1, 3]) / 8
def cond(i, f):
    out = np.zeros(4)
    for x in range(4):
        same = [y for y in range(4) if st[y][i] == st[x][i]]; out[x] = P2[same] @ f[same] / P2[same].sum()
    return out
X = np.array([s[0] for s in st], float); Y = np.array([s[1] for s in st], float)
ok("Example 3.2: E1 E2 X = X/4, E2 E1 X = Y/2", [np.abs(cond(0, cond(1, X)) - X / 4).max(), np.abs(cond(1, cond(0, X)) - Y / 2).max()], 1e-15)
ok("Example 3.2: Gaussian halfspaces at angle pi/4 give P(equal signs) = 3/4", abs((1 - (np.pi / 4) / np.pi) - 0.75), 1e-15)

# ------------------------------------------------------------------ Sec 4: moment interval, cut relaxation, gap identity
e = []
for _ in range(200):
    p, z = law(); z += rng.normal() * rng.random(); t, q = p @ z, p @ z ** 2; m = p @ np.maximum(z, 0); s = np.sqrt(q)
    lo, hi, D_ = kk.moment_interval(t, q)
    e += [max(0, lo - m), max(0, m - hi), abs(kk.cut_relaxation_value(t, q) - hi), abs(hi - m - p @ (np.abs(z) - s) ** 2 / (4 * s))]
ok("Thm 4.1: t_+ <= m <= (t+s)/2, rank-two relaxation value, gap identity (7)", e, 1e-12)
t, q = 0.3, 1.7; s = np.sqrt(q); pp = (1 + t / s) / 2
ok("Thm 4.1 sharpness: +-s law attains the upper end, {0, q/t} law the lower end",
   [abs(pp * s - (t + s) / 2), abs((t * t / q) * (q / t) - t)], 1e-15)

# ------------------------------------------------------------------ Sec 5: bounded-effect restart and filtered state
e1, e2, e3 = [], [], []
for sgn in (1, -1):
    for _ in range(100):
        p, z = law(9); z = z + sgn * abs(rng.normal()); t, q = p @ z, p @ z ** 2; m = p @ np.maximum(z, 0)
        R = kk.restart(t, q); nu = p * (z - R["c"]) ** 2 / R["v"]; g = R["g"](z)
        e1.append(abs(m - (R["lo"] + R["Delta"] * nu @ g))); e2 += [max(0, -g.min()), max(0, g.max() - 1)]
        e3.append(abs(nu.sum() - 1))
ok("Thm 5.1: m = t_+ + Delta E_nu g_c, both anchor signs (eq. 9)", e1, 1e-12)
ok("Thm 5.1: 0 <= g_c <= 1 and nu is a probability law", e2 + e3, 1e-12)
# filtered moment update (Sec. 5.1) on a finite law of a 3-slot activation vector
p = rng.random(6); p /= p.sum(); A = rng.random((6, 3)); w = rng.standard_normal(3); Z = A @ w; t, q = p @ Z, p @ Z ** 2
c = np.sign(t) * np.sqrt(q); v = p @ (Z - c) ** 2; nuA = (p * (Z - c) ** 2) @ A / v
E2 = np.einsum("x,xi,xa->ia", p, A, A); E3 = np.einsum("x,xi,xa,xb->iab", p, A, A, A)
upd = (np.einsum("a,b,iab->i", w, w, E3) - 2 * c * E2 @ w + c * c * (p @ A)) / v
uvec = np.concatenate([np.ones((6, 1)), A], 1); Gm = np.einsum("x,xi,xj->ij", p, uvec, uvec); bb = np.concatenate([[-c], w])
f = rng.random(6); Gf = np.einsum("x,xi,xj->ij", p * f, uvec, uvec)
gg = rng.random(); Pg = np.array([[gg, np.sqrt(gg * (1 - gg))], [np.sqrt(gg * (1 - gg)), 1 - gg]])
ok("Sec 5.1: filtered moment update, coherent slot state, G_f <= G, ancilla projection",
   [np.abs(nuA - upd).max(), abs(bb @ Gf @ bb / (bb @ Gm @ bb) - (p * (Z - c) ** 2 / v) @ f),
    max(0, -np.linalg.eigvalsh(Gm - Gf).min()), np.abs(Pg @ Pg - Pg).max()], 1e-12)

# ------------------------------------------------------------------ Sec 6: squared-shell state
e = []
for _ in range(200):
    p, z = law(8); t, q, r = p @ z, p @ z ** 2, p @ z ** 4; s = np.sqrt(q); m = p @ np.maximum(z, 0)
    Qd = kk.quartic(q, r); nu4 = p * (z * z - q) ** 2 / Qd["d4"]; gap = (t + s) / 2 - m
    R = kk.restart(t, q); nu = p * (z - R["c"]) ** 2 / R["v"]
    e += [abs(gap - Qd["d4"] / (4 * s ** 3) * nu4 @ Qd["effect"](z)), max(0, gap - Qd["gap_bound"]),
          abs(R["Delta"] * nu @ R["g"](z) + Qd["d4"] / (4 * s ** 3) * nu4 @ Qd["effect"](z) - R["Delta"])]
ok("Sec 6: gap = d4/(4 s^3) E_nu4 s^2/(|Z|+s)^2 <= d4/(4 s^3), and the quadratic/quartic consistency constraint", e, 1e-12)
roots = np.array([-1.3, 0.2, 0.9]); pr3 = np.array([0.2, 0.5, 0.3]); hcoef = np.polyfit(roots, np.maximum(roots, 0), 2)
ok("Sec 6: moment-null mechanism (support on the roots of p: E Z_+ = E h(Z))", abs(pr3 @ np.maximum(roots, 0) - pr3 @ np.polyval(hcoef, roots)), 1e-12)

# ------------------------------------------------------------------ Sec 7-8: spectral hierarchy
pi = np.pi; q = 1 - 1 / pi; r = 4.5 - 8 / pi; s6 = 37.5 - 88 / pi; s8 = 472.5 - 1280 / pi; M = [1, 0, q, 0, r, 0, s6, 0, s8]
o1, o2 = kk.spectral(0.0, M, 1), kk.spectral(0.0, M, 2)
table = dict(L1=0.2013452858, L2=0.2331052396, U2=0.3531900565, U1=0.4128226356, R1=0.6613046299, R2=0.7097568528)
got = dict(L1=o1["L"], L2=o2["L"], U2=o2["U"], U1=o1["U"], R1=o1["R"], R2=o2["R"])
ok("Sec 9.1 table: L1, L2, U2, U1, R1, R2 (rounded to 1e-10)", [abs(got[k] - v) for k, v in table.items()], 6e-11)
ok("R1 = (r - q^2)/(4 q^{3/2})", abs(o1["R"] - (r - q * q) / (4 * q ** 1.5)), 1e-12)
e, ex = [], []
for _ in range(40):
    p, z = law(6); z += 0.3 * rng.normal(); m = p @ np.maximum(z, 0); t = p @ z; Mz = [p @ z ** j for j in range(25)]
    prevU, prevL = np.inf, -np.inf
    for k in range(1, 4):
        o = kk.spectral(t, Mz, k, dps=80)
        e += [max(0, m - o["U"]), max(0, o["L"] - m), max(0, o["U"] - prevU - 1e-12), max(0, prevL - o["L"] - 1e-12), max(0, o["U"] - m - o["R"])]
        prevU, prevL = o["U"], o["L"]
    import mpmath as mp; mp.mp.dps = 80                                            # exact moments to degree 24
    Mx = [mp.fsum(mp.mpf(float(pa)) * mp.mpf(float(za)) ** j for pa, za in zip(p, z)) for j in range(25)]
    o6 = kk.spectral(t, Mx, 6, dps=80); ex += [abs(o6["U"] - m), abs(o6["L"] - m)]     # 6 atoms: cyclic space stabilises
ok("Thm 7.1/7.2/8.1 on finite laws: L_k <= m <= U_k, monotone in k, U_k - m <= R_k", e, 1e-10)
ok("exact termination once the cyclic space stabilises (6 atoms, k = 6)", ex, 1e-12)
z3 = np.array([-1.0, 0.0, 2.0]); p3 = np.array([0.3, 0.3, 0.4]); M3 = [p3 @ z3 ** j for j in range(9)]     # the zero atom drops out of (lambda/q) tau
ok("three-atom law with an atom at 0: odd two-dimensional lower space already exact", abs(kk.spectral(p3 @ z3, M3, 2)["L"] - p3 @ np.maximum(z3, 0)), 1e-12)

# ------------------------------------------------------------------ Sec 9: pair overlap, Rayleigh closure, pair chain, three lengths
e = []
for _ in range(100):
    p = rng.random(7); p /= p.sum(); P = rng.random(7) * (rng.random(7) > 0.2); Qv = rng.random(7) * (rng.random(7) > 0.2)
    o = kk.pair_overlap(p, P, Qv); e += [abs(o["m"] - (o["EP"] - o["kappa"] * o["Enu_inv"])), max(0, o["Enu_inv2"] - 1 / o["kappa"])]
ok("Prop 9.1: m = E P - kappa E_nu 1/max(P,Q), E_nu max^-2 <= 1/kappa", e, 1e-12)
Einv = 2 * integrate.quad(lambda x: np.exp(-x * x / 2) - np.exp(-x * x), 0, np.inf, epsabs=1e-14)[0]
mR = 1 / np.sqrt(2 * pi) - (1 / (2 * pi)) * Einv
ok("Sec 9.1: E 1/max(R1,R2) = sqrt(pi)(sqrt2 - 1) and E F = 1/(2 sqrt pi) (eq. 19)", [abs(Einv - np.sqrt(pi) * (np.sqrt(2) - 1)), abs(mR - 1 / (2 * np.sqrt(pi)))], 1e-12)
ok("Sec 9.1: baseline 0.3989423, correction 0.1168475, gate occupation 3/8",
   [abs(1 / np.sqrt(2 * pi) - 0.3989423), abs(1 / np.sqrt(2 * pi) - mR - 0.1168475), abs(0.25 + 0.125 - 3 / 8)], 6e-8)
Xg = rng.standard_normal((2, 4_000_000)); Fm = np.maximum(np.maximum(Xg[0], 0) - np.maximum(Xg[1], 0), 0).mean()
ok("Sec 9.1: Monte Carlo of E(X_+ - Y_+)_+ agrees with 1/(2 sqrt pi) (4e6 samples, 5 sigma)", abs(Fm - mR) / (0.33 / 2000), 5)
p = rng.random(5); p /= p.sum(); A = rng.random((5, 4)) * (rng.random((5, 4)) > 0.3); pw, rw = np.array([1.0, 0, 0.7, 0]), np.array([0, 0.5, 0, 1.2])
pi_, K, pairs, nu, cnd = kk.pair_chain(p, A, pw, rw); fpair = rng.standard_normal(len(pairs))
Sf = cnd.T @ fpair; condvar = nu @ ((cnd.T @ fpair ** 2) - Sf ** 2)
ok("Sec 9: pair chain K = S*S is a reversible Markov contraction with eq. (18)",
   [np.abs(K.sum(1) - 1).max(), np.abs(pi_[:, None] * K - (pi_[:, None] * K).T).max(), abs(fpair @ (pi_ * fpair) - fpair @ (pi_ * (K @ fpair)) - condvar)], 1e-12)
A2 = np.array([[1, 1, 0, 0], [0, 0, 1, 1.0]]); _, K2, pr2, _, _ = kk.pair_chain(np.array([.5, .5]), A2, np.array([1, 0, 1, 0.]), np.array([0, 1, 0, 1.]))
ok("Sec 9: (1,1,0,0)/(0,0,1,1) example gives two disconnected pair states", abs(np.sort(np.linalg.eigvals(K2).real)[-2] - 1), 1e-12)
u, v = rng.standard_normal(6), rng.standard_normal(6); Xs = rng.standard_normal((2_000_000, 6))
val = np.maximum(np.maximum(Xs @ u, 0) - np.maximum(Xs @ v, 0), 0); ang = np.arccos(u @ v / np.linalg.norm(u) / np.linalg.norm(v))
csgn = np.mean(np.sign(Xs @ u) * np.sign(Xs @ v))
ok("eq. (20) three lengths, by Monte Carlo (z-score)", abs(val.mean() - kk.three_lengths(u, v)) / (val.std() / np.sqrt(len(val))), 5)
ok("Sec 9.2: sign correlation 1 - 2 alpha/pi (z-score)", abs(csgn - (1 - 2 * ang / pi)) / (1 / np.sqrt(len(Xs))), 5)
cc = 0.3; ok("Sec 9.2: heat-bath (E1+E2)/2 eigenvalues 1, (1+c)/2, (1-c)/2, 0", np.abs(kk.two_sign_heat_bath(cc) - np.array([1, (1 + cc) / 2, (1 - cc) / 2, 0])), 1e-12)

# ------------------------------------------------------------------ Sec 10: conditional brackets, injection
e, e22 = [], []
for _ in range(50):
    na, nf = 32, 5; p = rng.random(na); p /= p.sum(); Zv = rng.standard_normal((na, nf)) + 0.3 * rng.standard_normal(nf)
    m = p @ np.maximum(Zv, 0); labs = [np.zeros(na, int)]
    for lev in range(5): labs.append(labs[-1] * 2 + rng.integers(0, 2, na))
    los, his = zip(*[kk.conditional_brackets(p, Zv, lb) for lb in labs])
    for a_, b_ in zip(los, los[1:]): e.append(max(0, (a_ - b_).max()))
    for a_, b_ in zip(his, his[1:]): e.append(max(0, (b_ - a_).max()))
    for lo_, hi_ in zip(los, his):
        e += [max(0, (lo_ - m).max()), max(0, (m - hi_).max()), max(0, np.sqrt(np.mean(((lo_ + hi_) / 2 - m) ** 2)) - kk.rms_bound(lo_, hi_))]
    G0 = (his[2] - los[2]) / 2; dC = (G0 - (his[3] - los[3]) / 2)
    e22.append(abs((np.sum(G0 ** 2) - np.sum((G0 - dC) ** 2)) - (2 * G0 @ dC - dC @ dC)))
ok("Thm 10.1: l_k up, u_k down, l <= m <= u, midpoint RMS bound (21)", e, 1e-12)
ok("eq. (22): squared certificate drop identity", e22, 1e-12)
aa, bb_ = rng.random(2); ok("lower residual: E[Z_+|C] - (E[Z|C])_+ = min(E Z_+, E Z_-)", abs(aa - max(aa - bb_, 0) - min(aa, bb_)), 1e-15)
zz, hh = rng.standard_normal(1000), rng.standard_normal(1000)
ok("eq. (23) injection identity", np.abs(0.5 * (np.maximum(zz + hh, 0) + np.maximum(zz - hh, 0)) - np.maximum(zz, 0) - kk.injection(zz, hh)).max(), 1e-14)
inj = integrate.quad(lambda x: x / 2 * np.exp(-x * x / 2) * 2 / np.sqrt(2 * pi), 0, np.inf)[0]
ok("Sec 10.1 example: E U/2 = 1/sqrt(2 pi)", abs(inj - 1 / np.sqrt(2 * pi)), 1e-12)

# ------------------------------------------------------------------ Sec 11: sparsifier congruence, law-dependent process
N = 6; wed = {(i, j): 1.0 for i in range(N) for j in range(i + 1, N)}; Lm = kk.kikuchi_laplacian(N, wed)
Ls = kk.kikuchi_laplacian(N, kk.sparsify(wed, 60, rng)); lo_, hi_ = kk.relative_spectrum(Lm, Ls)
Cmap = rng.standard_normal((1 << N, 9)); lo2, hi2 = kk.relative_spectrum(Cmap.T @ Lm @ Cmap, Cmap.T @ Ls @ Cmap)
print(f"      sampled sparsifier (60 of 15 edges with replacement): relative spectrum [{lo_:.3f}, {hi_:.3f}]")
ok("Prop 11.1: congruence keeps the relative spectrum inside [1-eps, 1+eps]", [max(0, lo_ - lo2 - 1e-9), max(0, hi2 - hi_ - 1e-9)], 1e-9)
Lg = -rng.random((5, 5)); Lg = (Lg + Lg.T) / 2; np.fill_diagonal(Lg, 0); np.fill_diagonal(Lg, -Lg.sum(1)); pv = rng.random(5); pv /= pv.sum()
Qp = np.diag(1 / pv) @ Lg; Hp = np.diag(pv ** -.5) @ Lg @ np.diag(pv ** -.5)
ok("Sec 11: Q_p 1 = 0, detailed balance in L^2(p), ker H_p = sqrt p",
   [np.abs(Qp @ np.ones(5)).max(), np.abs(np.diag(pv) @ Qp - (np.diag(pv) @ Qp).T).max(), np.abs(Hp @ np.sqrt(pv)).max()], 1e-12)

# ------------------------------------------------------------------ Sec 12: Schur refinement
e = []
for _ in range(30):
    nd = 9; Bm = rng.standard_normal((nd, nd)); H = Bm @ Bm.T + 0.5 * np.eye(nd); gam = np.linalg.eigvalsh(H).min()
    bv = rng.standard_normal(nd); Cm = rng.standard_normal((3, nd)); Pset = [0, 1, 2]; J = [3, 4]; Qset = list(range(3, nd))
    sc = kk.schur(H, bv, Cm, Pset, Qset, gam); en = kk.schur_enrich(H, bv, Cm, Pset, J, gam)
    sc2 = kk.schur(H, bv, Cm, Pset + J, [i for i in range(nd) if i not in Pset + J], gam)
    yPJ = Cm[:, Pset + J] @ np.linalg.solve(H[np.ix_(Pset + J, Pset + J)], bv[Pset + J])
    e += [sc["identity_err"], max(0, sc["err2"] - sc["cert_spec"]), max(0, sc["cert_spec"] - sc["cert_trace"]), max(0, -sc["EP"]),
          max(0, -np.linalg.eigvalsh(sc["UP"]).min()), np.abs(en["y"] - yPJ).max(), abs(en["E"] - sc2["EP"]), np.abs(en["U"] - sc2["UP"]).max(),
          abs(sc["cert_trace"] - sc2["cert_trace"] - en["trace_decrease"]), max(0, sc2["cert_trace"] - sc["cert_trace"])]
ok("Thm 12.1: identity (24), certificate (25), enrichment updates, trace decrease, monotonicity", e, 1e-9)
H2 = np.eye(2); b2 = np.ones(2); C2 = np.array([[1.0, -1.0]])
ok("Sec 12.1 example: empty trial 0, first coordinate 1, exact 0", [abs(C2 @ np.linalg.solve(H2, b2)).max(), abs((C2[:, [0]] @ b2[[0]])[0] - 1)], 1e-15)
Lh, dd = 4, 3; Ts = [np.linalg.qr(rng.standard_normal((dd, dd)))[0] * 0.9 for _ in range(Lh)]
Mm = np.eye((Lh + 1) * dd)
for l, T in enumerate(Ts): Mm[(l + 1) * dd:(l + 2) * dd, l * dd:(l + 1) * dd] = -T
ok("Sec 12.2: H = M*M >= 2/((L+1)(L+2)) for contractive T", max(0, 2 / ((Lh + 1) * (Lh + 2)) - np.linalg.eigvalsh(Mm.T @ Mm).min()), 1e-12)
Pm = rng.standard_normal((6, 6)); Pm = (Pm + Pm.T) / 2; Pm *= 0.7 / np.abs(np.linalg.eigvalsh(Pm)).max(); Sm = np.eye(6) - Pm; rv = rng.standard_normal(6)
uu, rr = np.zeros(6), rv.copy(); e = []
for _ in range(25):
    uu, rr = uu + rr, Pm @ rr; e += [np.abs(rr - (rv - Sm @ uu)).max(), max(0, np.linalg.norm(np.linalg.solve(Sm, rv) - uu) - np.linalg.norm(rr) / 0.3)]
ok("Sec 12.2: residual recycling r_t = r - S u_t and the gap bound", e, 1e-12)

# ------------------------------------------------------------------ Sec 13: observable code
e = []
for _ in range(30):
    dm = 6; Vq = np.linalg.qr(rng.standard_normal((dm, dm)))[0]; lam = np.concatenate([[0, 0], 0.5 + rng.random(dm - 2)]); Kop = Vq @ np.diag(lam) @ Vq.T
    Pm = Vq[:, :2] @ Vq[:, :2].T; Os = []
    for _ in range(4):
        B = rng.standard_normal((dm, dm)); B = (B + B.T) / 2; ci = rng.normal(); O = B - Pm @ B @ Pm + ci * Pm; Os.append(O)
    Rr = rng.standard_normal((dm, dm)); rho = Rr @ Rr.T; rho = (Pm @ rho @ Pm + 0.02 * rho); rho /= np.trace(rho)
    ob = kk.observable_code_bound(Kop, Os, rho); e += [max(0, ob["actual"] - ob["bound"]), max(0, ob["bound"] - ob["bound_eps"] - 1e-12)]
ok("Thm 13.1: RMS <= 2b sqrt(q(1-q)) + d q <= 2b sqrt(eps/gamma) + d eps/gamma", e, 1e-12)

print(f"\n{'ALL PASS' if not FAIL else 'FAILED: ' + ', '.join(FAIL)}")
sys.exit(1 if FAIL else 0)
