"""The full-history wall graph (companion Section 3): one-particle gap bound from all n L gates at once.

  python scripts/s24_fullgates.py DATA OUT NET [NSAMP]

Accumulates the covariance of the stacked gate vector g = (g_1, ..., g_L) in {0,1}^(nL) and every gate's flip rate under
the Mehler rotation (t = T_ROT), then reports the top eigenvalues of the rate-normalised covariance
Gamma^(-1/2) Cov(g) Gamma^(-1/2) (Gamma = one-way flux = rate / 2): lambda_lin = 1 / lambda_max bounds the gap of the
full-history wall chain from above (Proposition 3.1 of Stage 24). Also: lambda_max of the gate correlation of the whole
stack, the share of each layer in the top eigenvector, and the same bound for windows of consecutive layers."""
import sys, os, math, time, numpy as np
from scipy.sparse.linalg import eigsh

BS, T_ROT = 2048, 0.01


def run(D, OUT, net, N):
    W = [np.ascontiguousarray(w, dtype=np.float32) for w in np.load(f"{D}/W_off{net}.npy")]; L, n = len(W), W[0].shape[0]
    d = W[0].shape[1]; NL = n * L
    S1 = np.zeros(NL); S2 = np.zeros((NL, NL)); flips = np.zeros(NL); t0 = time.time()
    for b in range(N // BS):
        r = np.random.default_rng([51_000 + net, b])
        X = r.standard_normal((d, BS), dtype=np.float32); Y = r.standard_normal((d, BS), dtype=np.float32)
        Xt = X * np.float32(math.cos(T_ROT)) + Y * np.float32(math.sin(T_ROT))
        G = np.empty((NL, BS), np.float32); h, ht = X, Xt
        for l, Wl in enumerate(W):
            z = Wl @ h; zt = Wl @ ht; G[l * n:(l + 1) * n] = (z > 0); flips[l * n:(l + 1) * n] += np.sum((z > 0) != (zt > 0), 1)
            h = np.maximum(z, 0); ht = np.maximum(zt, 0)
        S1 += G.sum(1, dtype=np.float64); S2 += G @ G.T
    m = S1 / N; S2 /= N; S2 -= np.outer(m, m); C = S2; rate = flips / (N * T_ROT); v = np.diag(C).copy()   # in place (memory)
    live = (v > 1e-6) & (rate > 0); idx = np.flatnonzero(live)
    Cl = C[np.ix_(idx, idx)].astype(np.float32); del C, S2; sd = np.sqrt(v[idx]).astype(np.float32); Gm = (rate[idx] / 2).astype(np.float32)
    R = Cl / sd[:, None]; R /= sd[None, :]; Mr = Cl; Mr /= np.sqrt(Gm)[:, None]; Mr /= np.sqrt(Gm)[None, :]
    er, ur = eigsh(R, k=6, which="LA"); em, um = eigsh(Mr, k=6, which="LA")
    er, ur, em, um = er[::-1], ur[:, ::-1], em[::-1], um[:, ::-1]
    lay = idx // n; share = np.array([np.sum(um[lay == l, 0] ** 2) for l in range(L)])
    lines = [f"net {net}: N = {N}, {len(idx)} live gates of {NL}, {time.time() - t0:.0f}s",
             f"  gate correlation of the whole stack: top eigenvalues " + " ".join(f"{x:.1f}" for x in er) + f" (MP noise edge {(1 + math.sqrt(len(idx) / N)) ** 2:.2f})",
             f"  rate-normalised covariance: top eigenvalues " + " ".join(f"{x:.1f}" for x in em) + f" -> full-history wall-chain gap <= 1/lambda_max = {1 / em[0]:.4f} per unit angle",
             "  layer shares of the top rate-normalised eigenvector: " + " ".join(f"{s_:.3f}" for s_ in share)]
    for w in (1, 2, 4, 8, 16):
        best = []
        for l0 in range(0, L - w + 1, max(1, w // 2)):
            sel = (lay >= l0) & (lay < l0 + w); sub = Mr[np.ix_(sel, sel)]
            best.append((l0, 1 / eigsh(sub, k=1, which="LA")[0][0]))
        lines.append(f"  window of {w:2d} layers: gap bound per start layer " + " ".join(f"{l0 + 1}:{g:.3f}" for l0, g in best))
    txt = "\n".join(lines); print(txt, flush=True)
    with open(f"{OUT}/s24_fullgates_net{net}.txt", "w") as f: f.write(txt + "\n")


if __name__ == "__main__":
    D, OUT, net = sys.argv[1], sys.argv[2], int(sys.argv[3]); N = int(sys.argv[4]) if len(sys.argv) > 4 else 1 << 16
    run(D, OUT, net, N)
