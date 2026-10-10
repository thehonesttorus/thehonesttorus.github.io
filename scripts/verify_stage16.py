"""Stage 16 checks: audit of the two companion notes (branching angular geometry; Gaussian holonomy), the partition-algebra
hierarchy (Brauer < even partitions < partitions), the carre-du-champ/replica covariance identity, the chaos structure of the
official network's neurons as functions of the input, and the variational choice of the Wick reference.

    python scripts/verify_stage16.py [DATA_DIR]
Transcript: notes/stage16/verify_stage16.txt
"""
import sys, os, itertools, numpy as np
from math import factorial, comb
from scipy.special import ndtr
from scipy import integrate
from scipy.linalg import expm

SQ2PI = np.sqrt(2 * np.pi)
phi = lambda x: np.exp(-0.5 * np.asarray(x, float) ** 2) / SQ2PI
rng = np.random.default_rng(16)


def hdr(s):
    print("\n" + "=" * 100 + "\n" + s + "\n" + "=" * 100, flush=True)


# ==================================================================================================================
hdr("1. Audit of the companion notes")
# 1a. spin-1/2 branching from spin j: p_+ = (j+1)/(2j+1), p_- = j/(2j+1), from K = c dim V_nu / (d dim V_lambda)
for j in (0.5, 1.0, 3.5):
    pp, pm = (2 * j + 2) / (2 * (2 * j + 1)), (2 * j) / (2 * (2 * j + 1))
    print(f"j={j}: p+={pp:.4f} (=(j+1)/(2j+1) {(j+1)/(2*j+1):.4f}), p-={pm:.4f}, sum {pp+pm:.4f}")
# 1b. the three-spin logical qubit and the exchange operators
k = lambda s: np.eye(8)[int(s, 2)]
L0 = (k('010') - k('100')) / np.sqrt(2); L1 = (2 * k('001') - k('010') - k('100')) / np.sqrt(6)
def swap(a, b):
    P = np.zeros((8, 8))
    for i in range(8):
        bits = list(format(i, '03b')); bits[a], bits[b] = bits[b], bits[a]; P[int(''.join(bits), 2), i] = 1
    return P
Bl = np.stack([L0, L1], 1)
S12, S23 = Bl.T @ swap(0, 1) @ Bl, Bl.T @ swap(1, 2) @ Bl
Zm, Xm, Ym = np.diag([1., -1.]), np.array([[0., 1.], [1., 0.]]), np.array([[0, -1j], [1j, 0]])
print(f"S12 = -Z: {np.allclose(S12, -Zm)}; S23 = Z/2 + (sqrt3/2) X: {np.allclose(S23, Zm/2 + np.sqrt(3)/2*Xm)};"
      f" [S12,S23] = -i sqrt3 Y: {np.allclose(S12@S23 - S23@S12, -1j*np.sqrt(3)*Ym)}")
# 1c. the SO(2) (signed Poisson) check of the multi-time display: E[Z_a^2 Z_b^2] = ab + 2a^2 + a
a_, b_ = 0.7, 1.9; M = 2_000_000
Za = rng.poisson(a_ / 2, M) - rng.poisson(a_ / 2, M); Zb = Za + rng.poisson((b_ - a_) / 2, M) - rng.poisson((b_ - a_) / 2, M)
print(f"E[Z_a^2 Z_b^2] MC {np.mean(Za**2*Zb**2):.4f} vs ab+2a^2+a = {a_*b_+2*a_**2+a_:.4f}; product of separate means {a_*b_:.4f}")
# 1d. the isotropic angular reference G_1 = diag(1/2, 1/2n, ...)
n5 = 5; th = rng.standard_normal((400000, n5)); th /= np.linalg.norm(th, axis=1, keepdims=True)
psi = np.hstack([np.ones((len(th), 1)), th]) / np.sqrt(2); G1 = psi.T @ psi / len(th)
print(f"G_1 diagonal {np.round(np.diag(G1), 4)} (pred 0.5, {1/(2*n5):.4f}); max off-diagonal {np.max(np.abs(G1 - np.diag(np.diag(G1)))):.1e}")
# 1e. holonomy of the symmetric-squeeze connection around a contact-preserving square loop (Q11 = 1 fixed)
def loop_rotation(a, steps=4000):
    path = [(0, 0), (a, 0), (a, a), (0, a), (0, 0)]; T = np.eye(2)
    for (u0, v0), (u1, v1) in zip(path[:-1], path[1:]):
        for s in range(steps):
            t0, t1 = s / steps, (s + 1) / steps; tm = (t0 + t1) / 2
            Q = lambda tt: np.array([[1, v0 + tt * (v1 - v0)], [v0 + tt * (v1 - v0), 1 + u0 + tt * (u1 - u0)]])
            dQ = Q(t1) - Q(t0); T = expm(0.5 * dQ @ np.linalg.inv(Q(tm))) @ T
    return np.arctan2(T[1, 0], T[0, 0]), np.max(np.abs(T @ T.T - np.eye(2)))
for a in (0.02, 0.05, 0.1):
    ang, orth = loop_rotation(a)
    print(f"square loop side {a}: rotation {ang:+.6e}, a^2/4 = {a*a/4:.6e}, ratio {ang/(a*a/4):+.4f}; |T T^t - I| {orth:.1e}")
# 1f. hard-core normalizer: only signed permutations keep span{z_i z_j, i<j}
def leaks(R):
    return max(abs(R[j, i] * R[k2, i]) for j in range(3) for k2 in range(3) if j != k2 for i in range(3))
Rp = np.array([[0, -1, 0], [0, 0, 1], [1, 0, 0.]]); Rr = np.linalg.qr(rng.standard_normal((3, 3)))[0]
print(f"coefficient of z_i^2 in (Rz)_j (Rz)_k, j != k: signed permutation {leaks(Rp):.2e}, random rotation {leaks(Rr):.3f}")
# 1g. conditionally Gaussian mixture: kappa4(X) = Cov(Q_ij,Q_kl) + Cov(Q_ik,Q_jl) + Cov(Q_il,Q_jk)
nn, MM = 3, 2_000_000
A = rng.standard_normal((MM, nn, 2)) * 0.6; Q = np.eye(nn) + np.einsum('mia,mja->mij', A, A)
X = np.einsum('mij,mj->mi', np.linalg.cholesky(Q), rng.standard_normal((MM, nn)))
def k4(i, j, k_, l):
    E = lambda *idx: np.mean(np.prod([X[:, t] for t in idx], 0))
    return E(i, j, k_, l) - E(i, j) * E(k_, l) - E(i, k_) * E(j, l) - E(i, l) * E(j, k_)
CQ = lambda i, j, k_, l: np.mean(Q[:, i, j] * Q[:, k_, l]) - np.mean(Q[:, i, j]) * np.mean(Q[:, k_, l])
for idx in [(0, 0, 0, 0), (0, 0, 1, 1), (0, 1, 2, 2)]:
    i, j, k_, l = idx
    print(f"kappa4{idx}: MC {k4(*idx):+.4f}; Cov(Q) formula {CQ(i,j,k_,l)+CQ(i,k_,j,l)+CQ(i,l,j,k_):+.4f}")
# 1h. a Wick transform does not respect a realization relation: B ~ Bern(p), b^2 - b = 0 but :b^2: - :b: = -2p(B-p)
p_ = 0.3; Bv = np.array([0., 1.])
A1 = lambda b: b - p_; A2 = lambda b: b * b - 2 * p_ * b + 2 * p_ ** 2 - p_
print(f"Bernoulli({p_}): :b^2:-:b: at B=0,1 -> {A2(Bv) - A1(Bv)}; -2p(B-p) -> {-2*p_*(Bv-p_)}")


# ==================================================================================================================
hdr("2. The commutant hierarchy on (R^n)^{(x)2k}: Brauer (O(n)) < even partitions (B_n) < partitions (S_n)")
def inv_dim_finite(n, power, signed):
    tot, cnt = 0.0, 0
    for perm in itertools.permutations(range(n)):
        fixed = [i for i in range(n) if perm[i] == i]
        if not signed:
            tot += len(fixed) ** power; cnt += 1
        else:
            nf = len(fixed)   # sign patterns: chi = sum over fixed points of signs
            for s in range(nf + 1):           # s = number of minus signs among fixed points
                tot += comb(nf, s) * (nf - 2 * s) ** power * 2 ** (n - nf); cnt += comb(nf, s) * 2 ** (n - nf)
    return tot / cnt
def haar_moment(n, power, M=200000):
    tr = np.empty(M)
    for b in range(0, M, 5000):
        G = rng.standard_normal((5000, n, n)); Qm, Rm = np.linalg.qr(G)
        Qm = Qm * np.sign(np.diagonal(Rm, axis1=1, axis2=2))[:, None, :]; tr[b:b+5000] = np.trace(Qm, axis1=1, axis2=2)
    return np.mean(tr ** power), np.std(tr ** power) / np.sqrt(M)
for n_, p2 in ((4, 4), (6, 6)):
    o, ose = haar_moment(n_ + 2, p2)
    print(f"2k={p2}, n={n_}: S_n {inv_dim_finite(n_, p2, False):.1f} (Bell: {15 if p2==4 else 203}); "
          f"B_n {inv_dim_finite(n_, p2, True):.1f} (even partitions: {4 if p2==4 else 31}); "
          f"O(n+2) Haar {o:.2f} +- {ose:.2f} (pairings: {3 if p2==4 else 15})")


# ==================================================================================================================
hdr("3. Carre du champ: Cov(f,g) = int_0^1 E[grad f(X) . grad g(X_rho)] d rho, and the chaos structure of the neurons")
# 3a. the replica (accumulated product-defect) identity on a small He network
ns, Ls = 32, 3
Ws = [rng.standard_normal((ns, ns)) * np.sqrt(2 / ns) for _ in range(Ls)]
def fwd_jac(x):
    h = x; J = np.broadcast_to(np.eye(ns), (len(x), ns, ns)).copy()
    for W in Ws:
        z = h @ W.T; g = (z > 0).astype(float); h = z * g; J = g[:, :, None] * np.einsum('ij,mjk->mik', W, J)
    return h, J
Xa = rng.standard_normal((300000, ns)); ha, _ = fwd_jac(Xa[:, :])
Cmc = np.cov(ha.T)
nodes, wts = np.polynomial.legendre.leggauss(10); nodes = (nodes + 1) / 2; wts = wts / 2
Crep = np.zeros((ns, ns))
for rho, wq in zip(nodes, wts):
    X1 = rng.standard_normal((60000, ns)); X2 = rho * X1 + np.sqrt(1 - rho ** 2) * rng.standard_normal((60000, ns))
    _, J1 = fwd_jac(X1); _, J2 = fwd_jac(X2); Crep += wq * np.einsum('mik,mjk->ij', J1, J2) / len(X1)
print(f"n={ns}, L={Ls}: |Cov_MC - replica integral|_F / |Cov|_F = {np.linalg.norm(Cmc-Crep)/np.linalg.norm(Cmc):.4f}"
      f" (MC noise ~ {1/np.sqrt(60000):.3f}); diagonal ratio mean {np.mean(np.diag(Crep)/np.diag(Cmc)):.4f}")

# 3b. chaos structure of the official network's neurons as functions of the input (Hutchinson forward-mode JVPs)
D_DIR = sys.argv[1] if len(sys.argv) > 1 else "/tmp/claude-0/-home-user-thehonesttorus-github-io/9523ef6a-ff73-5ab2-8be3-275df2b5fcc6/scratchpad/data"
HAVE = os.path.exists(os.path.join(D_DIR, "W_off0.npy"))
if HAVE and not os.environ.get("SKIP_SLOW"):
    W0 = [np.ascontiguousarray(w.T).astype(np.float32) for w in np.load(os.path.join(D_DIR, "W_off0.npy"))]
    Lr, nr = len(W0), W0[0].shape[0]
    S = {k: np.zeros((Lr, nr)) for k in ("h", "hh", "g2", "pair")}; Npair = 0; Bt = 2048
    for it in range(12):
        x1 = rng.standard_normal((Bt, nr), dtype=np.float32); x2 = rng.standard_normal((Bt, nr), dtype=np.float32)
        v = rng.standard_normal((Bt, nr), dtype=np.float32)
        h1, h2, d1, d2 = x1, x2, v, v
        for l in range(Lr):
            z1 = h1 @ W0[l]; z2 = h2 @ W0[l]; g1 = z1 > 0; g2 = z2 > 0
            d1 = (d1 @ W0[l]) * g1; d2 = (d2 @ W0[l]) * g2; h1 = z1 * g1; h2 = z2 * g2
            for h, d in ((h1, d1), (h2, d2)):
                S["h"][l] += h.sum(0, dtype=np.float64); S["hh"][l] += (h.astype(np.float64) ** 2).sum(0); S["g2"][l] += (d.astype(np.float64) ** 2).sum(0)
            S["pair"][l] += (d1.astype(np.float64) * d2).sum(0)
        Npair += Bt
    N2 = 2 * Npair
    var = S["hh"] / N2 - (S["h"] / N2) ** 2; G = S["g2"] / N2; F = S["pair"] / Npair
    print(f"official network 0, {N2} inputs: per layer, totals over neurons (first chaos F, variance V, gradient energy G)")
    print(" layer  F/V(first-chaos share)  1-F/V(higher chaos)  <k>_{>=2} = 1+(G-V)/(V-F)")
    for l in range(Lr):
        Fv, Vv, Gv = F[l].sum(), var[l].sum(), G[l].sum()
        print(f"  {l+1:2d}       {Fv/Vv:.4f}               {1-Fv/Vv:.4f}               {1+(Gv-Vv)/max(Vv-Fv,1e-12):.2f}")


# ==================================================================================================================
hdr("4. What any estimator must add to the Gaussian closure: isotropic (field) and anisotropic parts, network 0")
if HAVE:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from whest.relu_gauss import gauss_chain
    Wn = np.load(os.path.join(D_DIR, "W_off0.npy")).astype(np.float64); mt = np.load(os.path.join(D_DIR, "truth_off0.npz"))["m"]
    rec = {}; outG = gauss_chain(Wn, K=10, record=rec)
    KP = os.path.join(os.path.dirname(D_DIR.rstrip("/")), "kprop3", "kprop3_0.npz")
    outX = np.load(KP)["out"] if os.path.exists(KP) else None
    print(" layer | truth - Gauss: field-explained share | exact chain - Gauss: field-explained share | chain captures of (truth-Gauss)")
    for l in (3, 7, 11, 15):
        r, s_ = rec[l]["alpha"], rec[l]["sigma"]
        Bf = np.stack([s_ * phi(r) * np.polynomial.hermite_e.hermeval(r, np.eye(8)[k]) for k in range(7)] + [s_ * ndtr(r), s_ * r * ndtr(r)], 1)
        def fshare(y):
            c, *_ = np.linalg.lstsq(Bf, y, rcond=None); return 1 - np.mean((y - Bf @ c) ** 2) / np.mean(y ** 2)
        dT = mt[l] - outG[l]
        line = f"  {l+1:2d}   |              {fshare(dT):.3f}                    |"
        if outX is not None:
            dX = outX[l] - outG[l]
            line += f"               {fshare(dX):.3f}                       |   {1 - np.mean((dT - dX)**2)/np.mean(dT**2):.4f}"
        print(line)


# ==================================================================================================================
hdr("5. Is the moment-matched (centred, BPHZ) Wick frame optimal? Minimal-sensitivity frames in one dimension")
def gjet(k_, m, v, c):
    s_ = np.sqrt(v); al = (m - c) / s_
    if k_ == 0: return s_ * (al * ndtr(al) + phi(al))
    if k_ == 1: return ndtr(al)
    return phi(al) * np.polynomial.hermite_e.hermeval(-al, np.eye(k_ - 1)[k_ - 2]) / s_ ** (k_ - 1)
def frame_series(kap, v, c, K):
    lc = np.zeros(K + 1); lc[2] = (kap[2] - v) / 2
    for j in range(3, min(len(kap), K + 1)): lc[j] = kap[j] / factorial(j)
    nu = np.zeros(K + 1); nu[0] = 1
    for m_ in range(1, K + 1): nu[m_] = sum(kk * lc[kk] * nu[m_ - kk] for kk in range(1, m_ + 1)) / m_
    return sum(nu[k_] * gjet(k_, kap[1], v, c) for k_ in range(K + 1))
gx, gw = np.polynomial.hermite_e.hermegauss(80); gw = gw / gw.sum()
for a, c in ((0.15, 0.3), (0.15, -0.8), (0.3, 0.5)):
    Z = gx + a * (gx * gx - 1); mom = [np.sum(gw * Z ** k_) for k_ in range(7)]; kap = [0.0] * 7
    for m_ in range(1, 7): kap[m_] = mom[m_] - sum(comb(m_ - 1, j - 1) * kap[j] * mom[m_ - j] for j in range(1, m_))
    gst = (-1 + np.sqrt(1 + 4 * a * (c + a))) / (2 * a)
    ex = integrate.quad(lambda g: (g + a * (g * g - 1) - c) * phi(g), gst, np.inf, epsabs=1e-14)[0]
    for K in (3, 4):
        vs = np.linspace(0.5 * kap[2], 1.5 * kap[2], 401); vals = np.array([frame_series(kap[:5], v, c, K) for v in vs])
        d = np.gradient(vals, vs); idx = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
        pms = ", ".join(f"{vs[i]/kap[2]:.2f}:{vals[i]-ex:+.1e}" for i in idx) or "none"
        print(f"a={a}, wall {c}, order {K}: centred frame error {frame_series(kap[:5], kap[2], c, K)-ex:+.2e};"
              f" stationary frames (v/k2: error) {pms}")


# ==================================================================================================================
hdr("6. The wall (leaf) decomposition E f = sum_walls p(0) E[|grad z|^2 df/dh | wall] with global Dirac energies (network 0)")
if HAVE:
    from whest.relu_gauss import hermite_relu
    def gate_m2(mu, C, K=12):
        s_ = np.sqrt(np.diag(C)); al = mu / s_; d = hermite_relu(al, K + 1)
        rho = C / np.outer(s_, s_); np.fill_diagonal(rho, 0.0); P = np.outer(d[1], d[1]); rk = np.ones_like(rho); fact = 1.0
        for k_ in range(1, K + 1):
            rk = rk * rho; fact *= k_; P += np.outer(d[k_ + 1], d[k_ + 1]) * rk / fact
        np.fill_diagonal(P, ndtr(al)); return P
    Lw, nw = Wn.shape[0], Wn.shape[1]; Kj = np.eye(nw); acc = np.zeros(nw); wall = []
    for l in range(Lw):
        WK = Wn[l] @ Kj @ Wn[l].T; Edir = np.diag(WK).copy(); Kj = gate_m2(rec[l]["mu"], rec[l]["C"]) * WK
        r, s_ = rec[l]["alpha"], rec[l]["sigma"]
        acc = (ndtr(r) * (Wn[l] @ acc) if l > 0 else 0) + phi(r) / s_ * Edir; wall.append(acc.copy())
        if l in (0, 1, 3, 15):
            print(f"  layer {l+1:2d}: wall-sum MSE {np.mean((acc-mt[l])**2):.2e} vs Gaussian closure {np.mean((outG[l]-mt[l])**2):.2e};"
                  f" mean Dirac energy / variance of z: {np.mean(Edir/s_**2):.3f}")
