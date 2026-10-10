"""Single-pass contraction of every statistic the Kikuchi pipeline consumes (one forward pass per sample batch).

  python scripts/kik_contract.py DATA NET TASK NBATCH        writes $OUT/kc_{NET}_{TASK}.npz
  python scripts/kik_contract.py --closure DATA NET           writes $OUT/kcg_{NET}.npz (Gaussian-closure moments)

Per batch of 4096 inputs X ~ N(0, I), one float32 forward pass gives every layer's pre-activation z_l. Accumulated:
  pw   (L, n, P)  power sums  sum z^p, p = 1..P (P = 8 below the last layer, 16 at the last), float64
  pos  (L, n)     sum z_+
  cells           final layer, nested reveal of K = 8 penultimate gates per final neuron (order fixed by a seeded
                  pilot, greedy single-gate width decrease, eq. 22): counts, sum Z, sum Z^2, sum Z_+ per cell
  M0, M1, M2      single-gate candidates: sum 1{A_i>0}, sum 1{A_i>0} Z_j, sum 1{A_i>0} Z_j^2   (n, n)
  pair            P = p_j.A, Q = r_j.A (positive / negative parts of the final row): sums of P, Q, P^2, Q^2, PQ,
                  min(P, Q), min/max
Power sums travel as float64 (the order-16 Hankel forms need it); the cell and candidate sums as float32 (relative
rounding 6e-8, far below their sampling error); counts as int32.
"""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

BS, K, PLOW, PTOP, PILOT = 4096, 8, 8, 16, 4


def load(D, net):
    return [np.ascontiguousarray(w, dtype=np.float32) for w in np.load(f"{D}/W_off{net}.npy")]


def forward(W, X):
    zs, h = [], X
    for Wl in W:
        z = Wl @ h; zs.append(z); h = np.maximum(z, 0)
    return zs


def gate_order(W, net):
    """Seeded pilot (same for every task of a net): per final neuron, the K penultimate gates with the largest
    single-gate decrease of the conditional width g = (sqrt(p q) - |t|)/2 summed over the two cells."""
    n = W[0].shape[0]; rng = np.random.default_rng(777_000 + net)
    M0 = np.zeros(n); M1 = np.zeros((n, n)); M2 = np.zeros((n, n)); S1 = np.zeros(n); S2 = np.zeros(n); N = 0
    for _ in range(PILOT):
        zs = forward(W, rng.standard_normal((n, BS), dtype=np.float32))
        G = (zs[-2] > 0).astype(np.float32); Z = zs[-1]
        M0 += G.sum(1); M1 += G @ Z.T; M2 += G @ (Z * Z).T; S1 += Z.sum(1); S2 += (Z * Z).sum(1); N += BS
    p1 = M0[:, None] / N; t1 = M1 / N; q1 = M2 / N; t = S1[None, :] / N; q = S2[None, :] / N
    g = lambda p, tt, qq: 0.5 * (np.sqrt(np.maximum(p * qq, 0)) - np.abs(tt))
    dec = g(1.0, t, q) - g(p1, t1, q1) - g(1 - p1, t - t1, q - q1)          # (gate i, neuron j)
    return np.argsort(-dec, axis=0)[:K].T.copy()                           # (n, K)


def closure(D, net):
    from whest.relu_gauss import relu_moments
    W = [np.asarray(w, dtype=np.float64) for w in np.load(f"{D}/W_off{net}.npy")]
    mu = np.zeros(W[0].shape[0]); C = W[0] @ W[0].T; T, Q, Mh = [], [], []
    for l in range(len(W)):
        T.append(mu.copy()); Q.append(np.diag(C) + mu * mu)
        m, Kh, _, _ = relu_moments(mu, C, 10); Mh.append(m)
        if l + 1 < len(W):
            mu = W[l + 1] @ m; C = W[l + 1] @ Kh @ W[l + 1].T
    np.savez(f"{os.environ.get('OUT', '.')}/kcg_{net}.npz", t=np.array(T), q=np.array(Q), m=np.array(Mh))
    print(f"net {net}: Gaussian closure moments written", flush=True)


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "--closure":
        closure(a[1], int(a[2])); sys.exit()
    D, net, task, nb = a[0], int(a[1]), int(a[2]), int(a[3])
    W = load(D, net); L, n = len(W), W[0].shape[0]; t0 = time.time()
    order = gate_order(W, net); t1 = time.time()
    Wp = np.maximum(W[-1], 0); Wr = np.maximum(-W[-1], 0)
    pw = [np.zeros((n, PTOP if l == L - 1 else PLOW)) for l in range(L)]; pos = np.zeros((L, n))
    nc = n << K; cnt = np.zeros(nc, np.int64); cz = np.zeros(nc); cz2 = np.zeros(nc); czp = np.zeros(nc)
    M0 = np.zeros(n); M1 = np.zeros((n, n)); M2 = np.zeros((n, n))
    pair = np.zeros((7, n)); base = (np.arange(n, dtype=np.int64) << K)[:, None]
    rng = np.random.default_rng(1_000_003 * (net + 1) + task)
    for b in range(nb):
        zs = forward(W, rng.standard_normal((n, BS), dtype=np.float32))
        for l in range(L):
            z = zs[l]; P = pw[l].shape[1]; x = z.copy()
            for p in range(P):
                pw[l][:, p] += x.sum(1, dtype=np.float64)
                if p + 1 < P: x *= z
            pos[l] += np.maximum(z, 0).sum(1, dtype=np.float64)
        A = np.maximum(zs[-2], 0); G = A > 0; Gf = G.astype(np.float32); Z = zs[-1]; Zp = np.maximum(Z, 0); Z2 = Z * Z
        code = np.zeros((n, BS), np.int64)
        for k in range(K):
            code |= G[order[:, k]].astype(np.int64) << k
        flat = (base + code).ravel()
        cnt += np.bincount(flat, minlength=nc); cz += np.bincount(flat, Z.ravel(), nc)
        cz2 += np.bincount(flat, Z2.ravel(), nc); czp += np.bincount(flat, Zp.ravel(), nc)
        M0 += Gf.sum(1, dtype=np.float64); M1 += Gf @ Z.T; M2 += Gf @ Z2.T
        Pp = Wp @ A; Qq = Wr @ A; mn = np.minimum(Pp, Qq); mx = np.maximum(Pp, Qq)
        for i, v in enumerate((Pp, Qq, Pp * Pp, Qq * Qq, Pp * Qq, mn, mn / np.where(mx > 0, mx, 1))):
            pair[i] += v.sum(1, dtype=np.float64)
        if b == 0:
            print(f"net {net} task {task}: pilot {t1 - t0:.1f}s, first batch {time.time() - t1:.1f}s", flush=True)
    f32 = lambda x: x.astype(np.float32)
    np.savez(f"{os.environ.get('OUT', '.')}/kc_{net}_{task}.npz", N=nb * BS, order=order.astype(np.int16),
             pw_low=np.array(pw[:-1]), pw_top=pw[-1], pos=pos,
             cnt=cnt.astype(np.int32), cz=f32(cz), cz2=f32(cz2), czp=f32(czp),
             M0=M0, M1=f32(M1), M2=f32(M2), pair=pair)
    print(f"net {net} task {task}: {nb * BS} samples in {time.time() - t0:.0f}s", flush=True)
