# Tightness of k*l <= N log2 N: uniform l-type over a binary alphabet (every l-word exactly m times, N = m 2^l).
# Exact BEST count: #cyclic words with marked origin = N t_w(G) prod_v (d_v-1)! / prod_e m_e!  (G = m x de Bruijn(2, l-1)),
# so #necklaces >= t_w(G) prod_v (d_v-1)! / prod_e m_e! / ... ; greedy pruning of window-local conflicts costs
# log2(N^2 2^l + 1) bits; periodic/self-conflicting words are negligible (<= N^2 2^{N/2+l}).
import numpy as np
from math import lgamma, log
def log2fact(x): return lgamma(x + 1)/log(2)
for l in range(4, 15):
    V = 2**(l - 1)
    for m in [l*l, 4*l*l]:
        N = m*2**l
        # matrix-tree for the de Bruijn multigraph: out-degree 2m, edges v -> (2v+a) mod V with multiplicity m
        if V <= 4096:
            L = np.zeros((V, V))
            for v in range(V):
                L[v, v] += 2*m
                for a in (0, 1): L[v, (2*v + a) % V] -= m
            sign, ld = np.linalg.slogdet(L[1:, 1:]); logt = ld/log(2)
        else:
            logt = float('nan')
        logS = logt + V*log2fact(2*m - 1) - 2**l*log2fact(m)          # = log2(#cyclic words with origin / N)
        k = logS - np.log2(N**2 * 2**l + 1) - 1
        print(f"l={l:2d} m={m:4d} N={N:8d}: log2 t_w={logt:9.1f} | log2 #type-class/N={logS:10.1f} | k >= {k:10.1f} = {k/N:.4f} N"
              f" | k*l/(N log2 N) >= {k*l/(N*np.log2(N)):.3f} | upper bound N h = {N:.0f} (h=1)")
