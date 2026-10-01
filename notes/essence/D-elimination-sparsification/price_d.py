"""T2 pricing (units = dense 1024^3 f32 products) of the age-multiresolution carriers at L = 16, n = 1024.
Per (source, target): dense young source = 7 products (FC). Factored source of resolution x = k/n: 3 materialisations
(Z, Y, T = A U^T) + 4 (n,n,k) contractions + 1 final (n,k,n) = 8x, frame transport G^T U = x, QR of n x k ~ 1.33 x^2
(kit: QR n x n 1.33 u), A <- A R^T for 3 factors = 3 x^2. Odometer: transport, QR and the final U^T once per block."""
import math
L, c = 16, 2.0
def k(a): return min(1.0, c / a)
def ages(t): return range(1, t + 1)
def causal(dense_frac=0.75, c_=c):
    tot = 0; res = 0
    for t in range(1, L):
        for a in ages(t):
            x = min(1.0, c_ / a) if a >= 3 else 1.0
            if a < 3 or x >= dense_frac:
                tot += 7; res += 1
            else:
                tot += 8 * x + x + 1.33 * x * x + 3 * x * x; res += x
    return tot, res
def odometer(c_=c):
    tot = 0; res = 0
    for t in range(1, L):
        old = [a for a in ages(t) if a >= 3]          # ages 3..t, births in order; Bentley-Saxe blocks by count
        tot += 7 * min(t, 2); res += min(t, 2)
        m = len(old); blocks = []; start = 0
        for j in reversed(range(m.bit_length())):     # binary decomposition, oldest (largest) blocks first
            if m >> j & 1:
                blocks.append(old[::-1][start:start + 2 ** j]); start += 2 ** j
        for b in blocks:
            x = min(1.0, c_ / min(b))
            tot += len(b) * (7 * x + 3 * x * x) + x + 1.33 * x * x + x; res += len(b) * x
    return tot, res
fc = 7 * sum(range(1, L))
for name, (tot, res) in [('FC dense', (fc, sum(range(1, L)))), ('causal c=2', causal()), ('causal c=1.5', causal(c_=1.5)),
                         ('odometer c=2', odometer()), ('odometer c=1.5', odometer(1.5))]:
    print(f"{name:16s} memory products {tot:6.1f} u  resolution {res:5.1f} n  ({tot/res:.2f} per unit resolution)  + arrow 30 u -> {tot+30:6.1f} u = {(tot+30)/1024:.3f} B")
