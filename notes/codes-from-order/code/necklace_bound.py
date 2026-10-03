# Orbit (Li-Boyle type) codes on the cycle Z_N: codewords |Psi_w> = N^{-1/2} sum_t |T^t w>.
# Theorem: if the code corrects erasure of an interval of length l, all codewords share the l-window law nu,
# and k = log2 dim <= N h_nu(l) <= N H_nu(l)/l <= N log2 p_nu(l)/l <= N log2 N / l.
# Here: exhaustive check over all binary necklaces of length N (largest type classes vs. the bound),
# greedy construction of actual codes, and explicit Knill-Laflamme verification for small N.
import numpy as np, sys
from collections import defaultdict

def rotations(N):
    x = np.arange(2**N, dtype=np.int64)
    R = np.empty((N, 2**N), dtype=np.int64)
    for t in range(N):
        R[t] = ((x << t) | (x >> (N - t))) & ((1 << N) - 1)
    return R

def windows(words, N, l):
    # cyclic l-windows of each word (rows), as integers; bit i of word = symbol at site i
    W = np.empty((len(words), N), dtype=np.int64)
    ext = words | (words << N)                      # unroll the cycle
    for t in range(N):
        W[:, t] = (ext >> t) & ((1 << l) - 1)
    return W

def H(p):
    p = p[p > 0]; return -(p*np.log2(p)).sum()

def type_stats(counts):
    # counts: vector over 2^l window values (cyclic counts, total N); nu = counts/N
    N = counts.sum(); nu = counts/N
    l = int(np.log2(len(counts)))
    # (l-1)-marginal: window value v -> prefix (low l-1 bits) ... symbol order: bit j = site t+j
    pre = np.zeros(2**(l-1)); np.add.at(pre, np.arange(2**l) & ((1 << (l-1)) - 1), nu)
    h = H(nu) - H(pre)                                # conditional entropy of last symbol given first l-1
    return h, H(nu), (counts > 0).sum()

def main(N, ls):
    R = rotations(N)
    canon = R.min(0)
    ndist = np.array([len(set(R[:, w])) for w in range(0)])  # placeholder
    reps = np.unique(canon)
    # aperiodic representatives: all N rotations distinct
    Rr = R[:, reps]
    aper = np.array([len(np.unique(Rr[:, i])) == N for i in range(len(reps))])
    reps = reps[aper]
    print(f"N={N}: {len(reps)} aperiodic binary necklaces")
    for l in ls:
        Wn = windows(reps, N, l)
        C = np.zeros((len(reps), 2**l), dtype=np.int32)
        for t in range(N):
            np.add.at(C, (np.arange(len(reps)), Wn[:, t]), 1)
        keys, inv, cnt = np.unique(C, axis=0, return_inverse=True, return_counts=True)
        # for each type class: log2(size) vs N*h_nu(l)
        worst = -1e9; best = None
        for c in np.argsort(-cnt)[:2000]:
            h, Hl, p = type_stats(keys[c].astype(float))
            slack = N*h - np.log2(cnt[c])
            worst = max(worst, -slack)
            if best is None or np.log2(cnt[c]) > best[0]: best = (np.log2(cnt[c]), N*h, Hl, p)
        print(f"  l={l}: largest class log2 size {best[0]:.2f} | its N*h_nu {best[1]:.2f} | N*H_nu/l {N*best[2]/l:.2f} | N log2 p/l {N*np.log2(best[3])/l:.2f}"
              f" | N log2 N/l {N*np.log2(N)/l:.2f} | max over classes of [log2 size - N h_nu] = {worst:.2f} (theorem: <= 0)")
    return reps

def build_code(N, l, reps):
    # pick the largest l-type class, then greedily choose necklaces pairwise satisfying recoverability:
    # no rotations of w, w' (w != w', or w = w' with t != t') agree outside an interval of length l.
    Wn = windows(reps, N, l)
    C = np.zeros((len(reps), 2**l), dtype=np.int32)
    for t in range(N): np.add.at(C, (np.arange(len(reps)), Wn[:, t]), 1)
    keys, inv, cnt = np.unique(C, axis=0, return_inverse=True, return_counts=True)
    cls = reps[inv == np.argmax(cnt)]
    full = (1 << N) - 1
    def rots(w): return [((w << t) | (w >> (N-t))) & full for t in range(N)]
    def outside_masks():
        ms = []
        for s in range(N):
            m = 0
            for j in range(l): m |= 1 << ((s + j) % N)
            ms.append(full & ~m)
        return ms
    OM = outside_masks()
    def clash(a, b, same):
        ra, rb = rots(a), rots(b)
        for i, x in enumerate(ra):
            for j, y in enumerate(rb):
                if same and i == j: continue
                d = x ^ y
                if any((d & m) == 0 for m in OM): return True
        return False
    good = [w for w in cls if not clash(w, w, True)]
    code = []
    for w in good:
        if all(not clash(w, v, False) for v in code): code.append(w)
    return cls, good, code

def kl_check(N, l, code):
    # explicit states; K = sites 0..l-1 (low bits). Check Tr_{K^c}|Psi_a><Psi_b| = delta_ab rho_K.
    full = (1 << N) - 1
    dimK = 2**l
    def psi(w):
        v = np.zeros(2**N)
        for t in range(N): v[((w << t) | (w >> (N-t))) & full] += 1
        return v/np.linalg.norm(v)
    # reshape: index = high bits (K^c) * 2^l + low bits (K)
    M = [psi(w).reshape(2**(N-l), dimK) for w in code]
    rho = [[M[a].T @ M[b] for b in range(len(code))] for a in range(len(code))]   # (dimK x dimK) = Tr_{K^c}|Psi_b><Psi_a| (real)
    ref = rho[0][0]
    err = max(np.abs(rho[a][b] - (ref if a == b else 0)).max() for a in range(len(code)) for b in range(len(code)))
    ev = np.linalg.eigvalsh(ref); ev = ev[ev > 1e-12]
    return err, -(ev*np.log2(ev)).sum()

if __name__ == '__main__':
    for N, ls in [(16, [3, 4, 5, 6, 8]), (20, [3, 4, 5, 6, 8, 10])]:
        reps = main(N, ls)
    for N, l in [(12, 3), (12, 4), (14, 4), (16, 4)]:
        R = rotations(N); canon = R.min(0); reps = np.unique(canon)
        reps = np.array([w for w in reps if len(set(((w << t) | (w >> (N-t))) & ((1 << N)-1) for t in range(N))) == N])
        cls, good, code = build_code(N, l, reps)
        msg = f"N={N} l={l}: type class {len(cls)}, self-recoverable {len(good)}, greedy code size {len(code)} (k={np.log2(max(len(code),1)):.2f})"
        if N <= 14 and len(code) >= 2:
            err, S = kl_check(N, l, code[:8])
            msg += f" | KL max deviation (first {min(8,len(code))} codewords) {err:.2e}; S(rho_K)={S:.3f}, bound k <= N S/l = {N*S/l:.2f}"
        print(msg, flush=True)
