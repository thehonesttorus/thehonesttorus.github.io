# Monte Carlo pre-activation statistics of an official network, for the oracle attribution of note XXXI.
#   python mcstats.py NET NSAMP SEED OUTPREFIX
# Accumulates, for every layer l (pre-activation z_l = W_l y_(l-1), y = relu), in two independent halves of the
# samples: the mean, variance, third and fourth cumulant diagonals and the (2,1) slice E[(z_i - m_i)^2 (z_c - m_c)].
# Writes OUTPREFIX_h0.npz, OUTPREFIX_h1.npz (each half) and OUTPREFIX_full.npz (all samples).
import sys, time, numpy as np
net, N, seed, outp = int(sys.argv[1]), int(float(sys.argv[2])), int(sys.argv[3]), sys.argv[4]
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float32)   # (L, n_out, n_in): z_l = Wcol[l] @ y_(l-1)
L, n, _ = Wcol.shape
B = 32768
rng = np.random.default_rng(seed)
acc = [dict(c=0, s1=np.zeros((L, n)), s2=np.zeros((L, n)), s3=np.zeros((L, n)), s4=np.zeros((L, n)),
            m11=np.zeros((L, n, n)), m21=np.zeros((L, n, n)), yfin=np.zeros(n)) for _ in range(2)]
m0 = None   # per-layer shift (first batch's mean): moments are accumulated for d = z - m0, so nothing cancels
t0 = time.time(); done = 0; nb = 0
while done < N:
    h = 0 if done < N // 2 else 1
    b = min(B, N - done, (N // 2 - done) if h == 0 else N - done)
    Y = rng.standard_normal((n, b), dtype=np.float32)
    a = acc[h]; a["c"] += b
    if m0 is None:
        m0 = np.zeros((L, n), np.float32); Yt = Y
        for l in range(L):
            Zt = Wcol[l] @ Yt; m0[l] = Zt.mean(1); Yt = np.maximum(Zt, 0.0)
    for l in range(L):
        Z = Wcol[l] @ Y
        D = Z - m0[l][:, None]
        D2 = D * D
        a["s1"][l] += D.sum(1, dtype=np.float64); a["s2"][l] += D2.sum(1, dtype=np.float64)
        a["s3"][l] += (D2 * D).sum(1, dtype=np.float64); a["s4"][l] += (D2 * D2).sum(1, dtype=np.float64)
        a["m11"][l] += D @ D.T
        a["m21"][l] += D2 @ D.T
        Y = np.maximum(Z, 0.0, out=Z)
    a["yfin"] += Y.sum(1, dtype=np.float64)
    done += b; nb += 1
    if nb % 32 == 0:
        print(f"{done} samples, {time.time() - t0:.0f}s", flush=True)


def finish(a):
    c = a["c"]; md = a["s1"] / c; e2 = a["s2"] / c; e3 = a["s3"] / c; e4 = a["s4"] / c   # raw moments of d
    var = e2 - md * md
    k3 = e3 - 3 * md * e2 + 2 * md ** 3
    k4 = e4 - 4 * md * e3 + 6 * md * md * e2 - 3 * md ** 4 - 3 * var * var
    E11 = a["m11"] / c; E21 = a["m21"] / c                     # E[d_i d_c], E[d_i^2 d_c]
    # E[(d_i - md_i)^2 (d_c - md_c)] = E[d_i^2 d_c] - 2 md_i E[d_i d_c] - md_c E[d_i^2] + 2 md_i^2 md_c
    D21 = E21 - 2 * md[:, :, None] * E11 - md[:, None, :] * e2[:, :, None] + 2 * (md * md)[:, :, None] * md[:, None, :]
    return dict(n=c, mu=(m0 + md).astype(np.float32), var=var.astype(np.float32), k3=k3.astype(np.float32),
                k4=k4.astype(np.float32), D21=D21.astype(np.float32), yfin=(a["yfin"] / c).astype(np.float32))


f0, f1 = finish(acc[0]), finish(acc[1])
full = {k: acc[0][k] + acc[1][k] for k in ("c", "s1", "s2", "s3", "s4", "m11", "m21", "yfin")}
ff = finish(full)
for tag, d in (("h0", f0), ("h1", f1), ("full", ff)):
    np.savez(f"{outp}_{tag}.npz", **d)
print(f"done {N} samples in {time.time() - t0:.0f}s", flush=True)
