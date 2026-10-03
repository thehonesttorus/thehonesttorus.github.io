# Fibonacci (Sturmian) statistics: entanglement spectrum of the Li-Boyle code = factor frequencies.
# (a) at most 3 distinct eigenvalues (three-distance), each of the form |p + q*alpha| (gap labels Z + alpha Z);
# (b) S(l) = log2 l + O(1);  (c) conditional entropy h(l) = freq(right-special factor) -> k <= N h(l) ~ c N / l.
import numpy as np
phi = (1 + 5**0.5)/2
alpha = 1/phi**2                          # frequency of the letter 0 in the Fibonacci word 0100101001001...
M = 2_000_000
n = np.arange(M)
s = (np.floor((n + 2)*alpha) - np.floor((n + 1)*alpha)).astype(np.int8)    # Sturmian coding of rotation by alpha
def factor_freqs(l):
    v = np.zeros(M - l + 1, dtype=np.int64)
    for j in range(l): v = (v << 1) | s[j:M - l + 1 + j]
    u, c = np.unique(v, return_counts=True)
    return u, c/c.sum()
def three_distance(l):
    # arcs of the circle cut by {-j alpha mod 1 : j = 0..l}: exact factor frequencies of length-l factors
    pts = np.sort(np.mod(-np.arange(l + 1)*alpha, 1.0))
    gaps = np.diff(np.concatenate([pts, [pts[0] + 1]]))
    return np.sort(gaps)
def in_gap_group(x, Q=400):
    # find integers p, q with |x - (p + q alpha)| < 1e-9 (gap-labelling group Z + alpha Z)
    for q in range(-Q, Q + 1):
        p = round(x - q*alpha)
        if abs(x - p - q*alpha) < 1e-9: return (p, q)
    return None
print(" l | #factors | distinct freqs | all in Z+aZ | S(l) - log2(l+1) | N h(l) * l / N")
for l in [2, 3, 5, 8, 13, 21, 34, 55, 89, 144]:
    if l <= 34:
        u, f = factor_freqs(l)
    td = three_distance(l)
    vals = np.unique(np.round(td, 10))
    labels = [in_gap_group(x) for x in vals]
    S = -(td*np.log2(td)).sum()
    # conditional entropy: H(l) - H(l-1) from exact three-distance spectra
    td1 = three_distance(l - 1); S1 = -(td1*np.log2(td1)).sum()
    h = S - S1
    emp = (f"(empirical {len(u)} factors, max |freq-exact| {np.abs(np.sort(f) - td).max():.1e})" if l <= 34 else "")
    print(f"{l:3d} | {len(td):5d} | {len(vals)} {np.round(vals,5)} | {all(x is not None for x in labels)} {labels} | {S - np.log2(l+1):+.4f} | {h*l:.4f} {emp}")
