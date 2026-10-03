# Li-Boyle finite Fibonacci codes (their App. E): cyclic seeds F0 with k0 zeros, k1 ones; F^(n) = n-fold
# inflation 1->10, 0->1. Compare the logical dimension of LB's code with the entropy and Betti bounds:
#   k <= N h_nu(l) <= N (b1(G_l) - 1) max_v nu(v),  l = correctable window f_n + 1.
import numpy as np, itertools
from collections import Counter
def infl(w, n):
    for _ in range(n): w = ''.join('10' if c == '1' else '1' for c in w)
    return w
def canon(w): return min(w[i:] + w[:i] for i in range(len(w)))
def law(w, l):
    ww = w + w[:l]; return Counter(ww[i:i+l] for i in range(len(w)))
def H(c):
    p = np.array(list(c.values()), float); p /= p.sum(); return -(p*np.log2(p)).sum()
fib = [1, 1]
for _ in range(30): fib.append(fib[-1] + fib[-2])
for (k0, k1, n) in [(3, 4, 4), (4, 5, 4), (5, 6, 5), (5, 7, 6), (6, 7, 6)]:
    m = k0 + k1
    seeds = sorted({canon(''.join('0' if i in z else '1' for i in range(m)))
                    for z in map(set, itertools.combinations(range(m), k0))})
    seeds = [s for s in seeds if len({s[i:] + s[:i] for i in range(m)}) == m]
    def nbrs(s):
        out = set()
        for i in range(m):
            j = (i + 1) % m
            if s[i] != s[j]:
                t = list(s); t[i], t[j] = t[j], t[i]; out.add(canon(''.join(t)))
        return out
    chosen = []
    for s in seeds:
        if all(s not in nbrs(c) for c in chosen): chosen.append(s)
    F = [infl(s, n) for s in seeds]
    N = len(F[0]); lw = fib[n + 1]; lc = fib[n] + 1
    same = all(law(f, lw) == law(F[0], lw) for f in F)
    nl, nl1 = law(F[0], lc), law(F[0], lc - 1)
    h = H(nl) - H(nl1)
    b1 = len(nl) - len(nl1) + 1
    maxnu = max(nl1.values())/N
    k = np.log2(len(chosen))
    print(f"(k0,k1,n)=({k0},{k1},{n}) N={N}: window laws identical (len {lw}): {same} | l={lc} | LB code k={k:.2f} bits "
          f"| entropy bound N h={N*h:.2f} | Betti bound N(b1-1)max nu={N*(b1-1)*maxnu:.2f} (b1={b1}) | k l/N={k*lc/N:.3f}", flush=True)
