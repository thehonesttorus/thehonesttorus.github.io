import numpy as np, itertools
from numpy.polynomial.hermite_e import hermegauss
# (a) Digit interleaving: {0,1}^{2n} -> (Z/2^n)^2, x from odd positions, y from even positions.
#     Check that every Gabber-Galil map preserves EVERY depth-m cylinder partition (m = 1..2n) of the binary tree.
n = 5; M = 2**n
def to_xy(bits):   # bits[0] is the first (coarsest) digit
    x = sum(bits[2*k] << k for k in range(n)); y = sum(bits[2*k+1] << k for k in range(n)); return x, y
def to_bits(x, y):
    b = []
    for k in range(n): b += [(x >> k) & 1, (y >> k) & 1]
    return tuple(b)
maps = [lambda x,y: (x, (y+2*x)%M), lambda x,y: (x, (y-2*x)%M), lambda x,y: (x, (y+2*x+1)%M), lambda x,y: (x, (y-2*x-1)%M),
        lambda x,y: ((x+2*y)%M, y), lambda x,y: ((x-2*y)%M, y), lambda x,y: ((x+2*y+1)%M, y), lambda x,y: ((x-2*y-1)%M, y)]
words = list(itertools.product([0,1], repeat=2*n))
ok = True
for f in maps:
    img = {w: to_bits(*f(*to_xy(w))) for w in words}
    for m in range(1, 2*n+1):
        pref = {}
        for w, v in img.items(): pref.setdefault(w[:m], set()).add(v[:m])
        if any(len(s) != 1 for s in pref.values()): ok = False
print(f"(a) all 8 Gabber-Galil maps send depth-m cylinders to depth-m cylinders for every m=1..{2*n}: {ok}")
# Sanov: A=[[1,2],[0,1]], B=[[1,0],[2,1]] generate a free group: check no short relation among reduced words up to length 8
A = np.array([[1,2],[0,1]]); B = np.array([[1,0],[2,1]]); Ai = np.array([[1,-2],[0,1]]); Bi = np.array([[1,0],[-2,1]])
gens = {'a':A,'A':Ai,'b':B,'B':Bi}; inv = {'a':'A','A':'a','b':'B','B':'b'}
rel = 0
for L in range(1, 9):
    for w in itertools.product('aAbB', repeat=L):
        if any(w[i+1] == inv[w[i]] for i in range(L-1)): continue
        P = np.eye(2, dtype=int)
        for c in w: P = P @ gens[c]
        if np.array_equal(P, np.eye(2, dtype=int)): rel += 1
print(f"    Sanov check: reduced words of length <= 8 equal to the identity: {rel}")
# (b) Galton-Watson duality for tanh (ordered / chaotic) vs ReLU (critical): low-degree leakage P(1<=Z_D<=t)
z, w = hermegauss(160); w = w / w.sum()
def he(k, x):
    h0, h1 = np.ones_like(x), x
    if k == 0: return h0
    for j in range(1, k): h0, h1 = h1, x*h1 - j*h0
    return h1
from math import factorial
def tanh_pgf(sw, sb, R=60):
    q = 1.0
    for _ in range(2000): q = sw**2 * np.sum(w * np.tanh(np.sqrt(q)*z)**2) + sb**2
    coef = np.array([(sw**2/q) * np.sum(w*np.tanh(np.sqrt(q)*z)*he(k, z))**2 / factorial(k) for k in range(R)])
    coef[0] += sb**2 / q
    return coef / coef.sum(), q
def relu_pgf(R=60):
    from math import comb, pi
    from fractions import Fraction
    a = np.zeros(R)
    def bh(k):
        r = Fraction(1)
        for i in range(k): r *= Fraction(1,2) - i; r /= i+1
        return r
    for k in range(R//2):
        a[2*k] += float(bh(k)) * (-1)**k
        if k >= 1: j = k-1; a[2*k] += comb(2*j, j) / (4**j * (2*j+1))
    a /= pi; a[1] += 0.5; return a / a.sum()
def compose(p, q, R):
    out = np.zeros(R); out[0] = p[-1]
    for c in p[::-1][1:]:
        out = np.convolve(out, q)[:R]; out[0] += c
    return out
def leak(p, D, t, R=60):
    b = p.copy(); res = []
    for d in range(1, D+1):
        if d > 1: b = compose(p, b, R)
        res.append(b[1:t+1].sum())
    return np.array(res)
for name, (p, q) in {"tanh ordered (sw=0.8, sb=0.3)": tanh_pgf(0.8, 0.3), "tanh chaotic (sw=2.5, sb=0)": tanh_pgf(2.5, 0.0),
                     "ReLU He (critical)": (relu_pgf(), None)}.items():
    m = np.sum(np.arange(len(p)) * p)
    L = leak(p, 32, 3)
    # extinction probability q_ext: smallest fixed point of the pgf
    s = 0.0
    for _ in range(20000): s = np.polyval(p[::-1], s)
    fq = np.sum(np.arange(1, len(p)) * p[1:] * s**np.arange(len(p)-1))
    print(f"(b) {name}: mean offspring chi_1={m:.3f}, extinction q={s:.3f}, f'(q)={fq:.3f}; "
          f"P(1<=Z_D<=3) at D=4,8,16,32: {L[3]:.2e} {L[7]:.2e} {L[15]:.2e} {L[31]:.2e}")
