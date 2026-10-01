"""Note 3's tilings row: memory rank law of three substitution sequences by spectral type.
Memory at window A = Toeplitz autocorrelation matrix T_A = [C(i - j)]_{i,j < A} of the centred +-1 sequence (the age-graded
correlation operator; its eigenvalues discretise the spectral measure at resolution 1/A).  k_q(A) = number of eigenvalues
holding a fraction q of tr T_A.  Fibonacci: pure point; Thue-Morse: singular continuous; Rudin-Shapiro: Lebesgue."""
import numpy as np
N = 1 << 21
def fib(N):
    s = 'a'
    while len(s) < N:
        s = ''.join('ab' if c == 'a' else 'a' for c in s)
    return np.array([1.0 if c == 'a' else -1.0 for c in s[:N]])
k = np.arange(N)
tm = 1.0 - 2.0 * (np.vectorize(lambda x: bin(x).count('1') & 1)(k))
rs = 1.0 - 2.0 * (np.vectorize(lambda x: bin(x & (x >> 1)).count('1') & 1)(k))
seqs = {'Fibonacci (pure point)': fib(N), 'Thue-Morse (sing. cont.)': tm, 'Rudin-Shapiro (Lebesgue)': rs}
As = [32, 64, 128, 256, 512, 1024]
for name, u in seqs.items():
    u = u - u.mean()
    F = np.fft.rfft(u, 2 * N); C = np.fft.irfft(F * np.conj(F))[:max(As) + 1] / N
    out = []
    for A in As:
        idx = np.abs(np.arange(A)[:, None] - np.arange(A)[None, :]); ev = np.sort(np.linalg.eigvalsh(C[idx]))[::-1]
        ev = np.maximum(ev, 0); c = np.cumsum(ev) / ev.sum()
        out.append((A, int(np.searchsorted(c, 0.9) + 1), int(np.searchsorted(c, 0.99) + 1)))
    k90 = np.array([o[1] for o in out]); k99 = np.array([o[2] for o in out])
    p90 = np.polyfit(np.log(As[2:]), np.log(k90[2:]), 1)[0]; p99 = np.polyfit(np.log(As[2:]), np.log(k99[2:]), 1)[0]
    print(f'{name:28s} k90(A) {list(k90)}  k99(A) {list(k99)}  growth exponents (A >= 128): {p90:.2f} / {p99:.2f}')
