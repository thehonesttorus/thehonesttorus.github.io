"""Checks accompanying the study of 'An isomorphism of the free group factors' (OpenAI, 23 Sep 2026).

    python scripts/verify_ffac.py          transcript: notes/stage17/verify_ffac.txt

1. Fox cocycle on F_3 and the prefix identity D(w_m) = T_m h (Lemma 4.2), free prefixes on a ball.
2. Prefix sums of free Haar unitaries: spectrum of T_m T_m^* / m -> free Poisson(1) (Lemma 4.3), random-matrix model.
3. Cutoff inversion (Lemma 4.1): the trade-off ||h||^2 = tau(R) <= 1/(dm), ||T h - delta_e||^2 = mu_m([0,d]).
4. The free-sum norm bound (Lemma 2.5) for conjugates of S = arg(C), against the free-CLT edge.
5. The mechanism of Proposition 3.1 in a matrix model: letter velocity ||s(h)|| vs word deviation ||s(D_w) - S||.
6. The amplification arithmetic of Section 6.
7. Deductions from the companion study: the steering (Pareto) law ||T_m h - delta_e|| ||h|| >= (1+o(1))/(2 sqrt m) for the
   coefficient problem (Tikhonov attains 1/2, the paper's cutoff 2/pi), the resulting cost of the Fox-prefix route against the
   universal bound (moved-letter count > 2/eps - 1), and ||C - A_w|| = 2.
"""
import numpy as np
from scipy.stats import unitary_group

rng = np.random.default_rng(17)


def hdr(s):
    print("\n" + "=" * 100 + "\n" + s + "\n" + "=" * 100, flush=True)


# ------------------------------------------------------------------------------------------------------------------
hdr("1. Fox cocycle on F_3: D_{gb} = D_g + lambda(g) D_b, and D(w_m) = T_m h for the prefix words")
def red(w):
    out = []
    for x in w:
        if out and out[-1] == -x: out.pop()
        else: out.append(x)
    return tuple(out)
def mul(u, v): return red(u + v)
def inv(u): return tuple(-x for x in reversed(u))
def shift(g, k):                       # (lambda(g) k)(b) = k(g^{-1} b): keys move to g*key
    out = {}
    for key, c in k.items():
        kk = mul(g, key); out[kk] = out.get(kk, 0.0) + c
    return out
def add(k1, k2, s=1.0):
    out = dict(k1)
    for key, c in k2.items(): out[key] = out.get(key, 0.0) + s * c
    return {a: b for a, b in out.items() if abs(b) > 1e-14}
def D(word, h):                        # h: dict letter j (1..3) -> coefficient vector
    tot, prefix = {}, ()
    for x in word:
        j = abs(x)
        if x > 0: term = shift(prefix, h.get(j, {}))
        else: term = shift(prefix, {kk: -c for kk, c in shift((x,), h.get(j, {})).items()})   # D_{x^-1} = -lambda(x^-1) h
        tot = add(tot, term); prefix = mul(prefix, (x,))
    return tot
hvec = {(): 0.7, (2,): -0.3, (1, 3): 0.2, (-2, 1): 0.5}
for m in (1, 3, 5):
    b = [None] + [tuple([2] * j + [3] + [-2] * j) for j in range(1, m + 1)]
    p = [()]
    for j in range(1, m + 1): p.append(mul(p[-1], (1,) + b[j]))
    w = mul(p[m], (1,))
    lhs = D(w, {1: hvec})
    rhs = {}
    for j in range(m + 1): rhs = add(rhs, shift(p[j], hvec))
    diff = add(lhs, rhs, -1.0)
    print(f"m={m}: |w| = {len(w)}, |D(w) - T_m h|_1 = {sum(abs(c) for c in diff.values()):.1e} over {len(lhs)} support points")
# cocycle identity on random words, and freeness of the prefixes on a ball
words = [red(tuple(rng.choice([1, -1, 2, -2, 3, -3], size=rng.integers(1, 9)))) for _ in range(200)]
hh = {1: {(): 0.4, (3,): 1.1}, 2: {(1,): -0.6}, 3: {(): 0.2, (-1, 2): 0.3}}
bad = 0
for _ in range(200):
    u, v = words[rng.integers(200)], words[rng.integers(200)]
    bad += sum(abs(c) for c in add(D(mul(u, v), hh), add(D(u, hh), shift(u, D(v, hh))), -1.0).values()) > 1e-12
print(f"cocycle identity D(uv) = D(u) + lambda(u) D(v): {200-bad}/200 random pairs exact")
m = 3; b = [None] + [tuple([2] * j + [3] + [-2] * j) for j in range(1, m + 1)]; p = [()]
for j in range(1, m + 1): p.append(mul(p[-1], (1,) + b[j]))
import itertools
seen, collide, trivial = {}, 0, 0
letters = [(j, s) for j in range(1, m + 1) for s in (1, -1)]
for L in range(1, 5):
    for seq in itertools.product(letters, repeat=L):
        if any(seq[i][0] == seq[i + 1][0] and seq[i][1] == -seq[i + 1][1] for i in range(L - 1)): continue
        g = ()
        for j, s in seq: g = mul(g, p[j] if s > 0 else inv(p[j]))
        trivial += (g == ()); collide += (g in seen); seen[g] = seq
print(f"reduced words of length <= 4 in p_1..p_3: {len(seen)} distinct images, {trivial} trivial, {collide} collisions (free basis)")

# ------------------------------------------------------------------------------------------------------------------
hdr("2. Prefix sums of free Haar unitaries: spectrum of Q_m/m = T_m T_m^*/m -> free Poisson(1)")
N = 1200
def nu_cdf(d):   # free Poisson(1): nu([0,d]) = (1/pi) [ (u/2) sqrt(4-u^2) + 2 arcsin(u/2) ] at u = sqrt(d)
    u = np.sqrt(np.atleast_1d(d)); return (u / 2 * np.sqrt(4 - u * u) + 2 * np.arcsin(u / 2)) / np.pi
cat = [1, 1, 2, 5, 14, 42]
spectra = {}
for m in (8, 32, 128):
    T = np.eye(N, dtype=complex)
    for _ in range(m): T += unitary_group.rvs(N, random_state=rng)
    ev = np.linalg.eigvalsh(T @ T.conj().T) / m; spectra[m] = ev
    mom = [np.mean(ev ** r) for r in range(1, 5)]
    print(f"m={m:3d}: moments of Q/m {', '.join(f'{x:.3f}' for x in mom)} (Catalan 1, 2, 5, 14);"
          f" mu([0,0.05]) {np.mean(ev <= 0.05):.4f} vs nu {nu_cdf(0.05)[0]:.4f}; max eig {ev.max():.2f} (edge 4)")

# ------------------------------------------------------------------------------------------------------------------
hdr("3. Cutoff inversion: ||h||^2 = tau(R_m) <= 1/(dm) and ||T_m h - delta_e||^2 = mu_m([0,d])")
for m in (32, 128):
    ev = spectra[m] * m
    for d in (0.02, 0.1):
        keep = ev > d * m
        print(f"m={m:3d}, d={d}: ||h||^2 = {np.mean(np.where(keep, 1 / np.where(keep, ev, 1), 0)):.4f} (bound 1/(dm) {1/(d*m):.4f});"
              f" error^2 = {np.mean(~keep):.4f} (nu([0,d]) {nu_cdf(d)[0]:.4f})")
print("asymptotics: error^2 ~ (2/pi) sqrt(d) and ||h||^2 = tau(R) ~ 2/(pi m sqrt(d)); both <= eta^2 needs d ~ (pi eta^2/2)^2, m ~ 4/(pi^2 eta^4):")
for eta in (0.3, 0.1, 0.03):
    d = (np.pi * eta ** 2 / 2) ** 2; mm = 4 / (np.pi ** 2 * eta ** 4)
    print(f"  eta={eta}: d ~ {d:.1e}, m ~ {mm:.1e}, word length |w_m| ~ m^2 ~ {mm**2:.1e}")

# ------------------------------------------------------------------------------------------------------------------
hdr("4. Free-sum norm: ||sum_g k(g) A_g S A_g^*|| <= 3 pi ||k||_2  (free CLT edge: 2 (pi/sqrt3) ||k||_2 for spread k)")
N = 600; theta = rng.uniform(-np.pi, np.pi, N); Sd = np.diag(theta)
for G, kind in ((40, "spread"), (40, "one big + small"), (8, "few")):
    k = rng.standard_normal(G)
    if kind == "one big + small": k[0] = 6.0
    Ssum = np.zeros((N, N), complex)
    for g in range(G):
        U = unitary_group.rvs(N, random_state=rng); Ssum += k[g] * (U @ Sd @ U.conj().T)
    nrm = np.max(np.abs(np.linalg.eigvalsh(Ssum)))
    print(f"G={G:2d} ({kind}): ||s(k)|| = {nrm:.3f}; 3 pi ||k|| = {3*np.pi*np.linalg.norm(k):.3f};"
          f" 2||k|| ||S||_2 + max|k| ||S|| = {2*np.linalg.norm(k)*np.pi/np.sqrt(3) + np.max(np.abs(k))*np.pi:.3f};"
          f" free-CLT edge {2*np.pi/np.sqrt(3)*np.linalg.norm(k):.3f}")

# ------------------------------------------------------------------------------------------------------------------
hdr("5. The mechanism of Proposition 3.1 in a matrix model: letters barely move while the word rotates by S")
N = 160
A = [unitary_group.rvs(N, random_state=rng) for _ in range(3)]; C = unitary_group.rvs(N, random_state=rng)
w_, V_ = np.linalg.eig(C); S = (V_ * np.angle(w_)) @ np.linalg.inv(V_); S = (S + S.conj().T) / 2
S = S - np.trace(S).real / N * np.eye(N)          # tau(S) = 0 exactly; the conjugates of S span a copy of l^2(Gamma)
ctr = lambda Y: Y - np.trace(Y) / N * np.eye(N)   # stay off the trivial (identity) direction of the conjugation action
def prefix_words(m):
    P = [np.eye(N, dtype=complex)]; A2p = np.eye(N, dtype=complex); A2 = A[1]; Acur = np.eye(N, dtype=complex)
    for j in range(1, m + 1):
        A2p = A2p @ A2; Bj = A2p @ A[2] @ A2p.conj().T; Acur = Acur @ A[0] @ Bj; P.append(Acur.copy())
    return P
for m, K in ((16, 24), (64, 96)):
    P = prefix_words(m)
    LT = lambda Y: ctr(sum(Pj @ Y @ Pj.conj().T for Pj in P))       # conjugation action of T = sum_j p_j
    LTs = lambda Y: ctr(sum(Pj.conj().T @ Y @ Pj for Pj in P))      # of T^*
    # spectrum bound X >= ||Q|| from the left-regular model T = sum_j P_j (as an operator)
    Tm = sum(P); X = 1.05 * np.linalg.norm(Tm @ Tm.conj().T, 2)
    # g(Q) = (1 - (1 - Q/X)^K)/Q = (1/X) sum_{i<K} (1 - Q/X)^i ; s(h) = L_{T^*} g(L_T L_{T^*}) (S)
    Y = S.copy(); acc = S / X
    for i in range(1, K):
        Y = Y - LT(LTs(Y)) / X; acc = acc + Y / X
    s_h = LTs(acc); s_Dw = LT(s_h)
    # coefficient-side norms from the left-regular model: h = T^* g(Q) delta_e, D_w - delta_e = -(1 - Q/X)^K delta_e
    ev = np.linalg.eigvalsh(Tm @ Tm.conj().T)
    g = (1 - (1 - ev / X) ** K) / np.maximum(ev, 1e-300)
    h2 = np.mean(g ** 2 * ev); e2 = np.mean((1 - ev / X) ** (2 * K))
    print(f"m={m}, K={K}: ||h||_2 = {np.sqrt(h2):.3f} -> letter velocity ||s(h)|| = {np.linalg.norm(s_h, 2):.3f} (<= 3pi||h|| = {3*np.pi*np.sqrt(h2):.2f});"
          f" ||D_w - delta_e||_2 = {np.sqrt(e2):.3f} -> word deviation ||s(D_w) - S|| = {np.linalg.norm(s_Dw - S, 2):.3f}; ||S|| = {np.linalg.norm(S, 2):.2f}")

# ------------------------------------------------------------------------------------------------------------------
hdr("6. Amplification arithmetic: N_s^t = N_{1 + (s-1) t^-2}")
amp = lambda s, t: 1 + (s - 1) / t ** 2
print(f"N_3^sqrt2 = N_{amp(3, np.sqrt(2)):.0f}, N_5^sqrt2 = N_{amp(5, np.sqrt(2)):.0f}: so N_3 ~ N_5 gives N_2 ~ N_3")

# ------------------------------------------------------------------------------------------------------------------
hdr("7. Deductions: steering law sqrt(m) ||T_m h - delta_e|| ||h|| >= 1/2 + o(1); cost of the route; ||C - A_w|| = 2")
from scipy.integrate import quad
# free Poisson(1) in the variable u = sqrt(x): dnu = sqrt(4 - u^2)/pi du on [0, 2], no singularity at 0
I = lambda f, lo=0.0, pts=None: quad(lambda u: f(u * u) * np.sqrt(4 - u * u) / np.pi, lo, 2, limit=500, points=pts)[0]
print("limit law; Tikhonov h = T^*(Q + md)^{-1} delta_e (the Pareto frontier) vs the paper's cutoff h = T^* Q^{-1} 1{Q > md} delta_e:")
for d in (1e-2, 1e-3, 1e-4, 1e-5):
    r = np.sqrt(d)
    e2t = I(lambda x: (d / (x + d)) ** 2, pts=[r]); at = I(lambda x: x / (x + d) ** 2, pts=[r])   # at = m ||h||^2
    e2c = nu_cdf(d)[0]; ac = I(lambda x: 1 / x, lo=r)
    print(f"  d={d:.0e}: Tikhonov err^2 {e2t:.2e}, m||h||^2 {at:.2e}, sqrt(m) err ||h|| = {np.sqrt(e2t * at):.4f} (-> 1/2);"
          f"  cutoff {np.sqrt(e2c * ac):.4f} (-> 2/pi = {2 / np.pi:.4f})")
print("finite m, random-matrix spectra of section 2 (trace in place of the delta_e spectral measure):")
for m in (32, 128):
    ev = spectra[m]
    for d in (0.03, 0.01):
        e2t = np.mean((d / (ev + d)) ** 2); at = np.mean(ev / (ev + d) ** 2)
        keep = ev > d; e2c = np.mean(~keep); ac = np.mean(np.where(keep, 1 / np.where(keep, ev, 1), 0))
        print(f"  m={m:3d}, d={d}: sqrt(m) err ||h||: Tikhonov {np.sqrt(e2t * at):.4f}, cutoff {np.sqrt(e2c * ac):.4f}")
print("cost of absorbing at accuracy eps (letters and word both within eps): moved-letter occurrences l_1(w) = m + 1")
for eps in (0.3, 0.1, 0.03):
    print(f"  eps={eps}: universal l_1(w) > 2/eps - 1 = {2 / eps - 1:.1f};  Fox-prefix route: sufficient with the 3pi bounds"
          f" m = 81 pi^4/(4 eps^4) = {81 * np.pi ** 4 / (4 * eps ** 4):.1e}, first-order necessary m ~ pi^4/(36 eps^4) = {np.pi ** 4 / (36 * eps ** 4):.1e}")
N = 400
A = [unitary_group.rvs(N, random_state=rng) for _ in range(3)]; C = unitary_group.rvs(N, random_state=rng)
Aw = A[0] @ A[1] @ A[2] @ A[1].conj().T @ A[0]                      # w = a b_1 a
z = np.linalg.eigvals(C.conj().T @ Aw)
print(f"||C - A_w|| = {np.linalg.norm(C - Aw, 2):.4f} (free: 2), ||C - A_w||_2 = {np.linalg.norm(C - Aw) / np.sqrt(N):.4f} (free: sqrt2 = {np.sqrt(2):.4f});"
      f" eigenvalues of C^* A_w: max gap on the circle {np.max(np.diff(np.sort(np.angle(z)))):.3f} (Haar: -> 0)")
