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
t0 = time.time(); done = 0
while done < N:
    h = 0 if done < N // 2 else 1
    b = min(B, N - done, (N // 2 - done) if h == 0 else N - done)
    Y = rng.standard_normal((n, b), dtype=np.float32)
    a = acc[h]; a["c"] += b
    for l in range(L):
        Z = Wcol[l] @ Y
        Z64 = Z.astype(np.float64)
        a["s1"][l] += Z64.sum(1); Z2 = Z64 * Z64
        a["s2"][l] += Z2.sum(1); a["s3"][l] += (Z2 * Z64).sum(1); a["s4"][l] += (Z2 * Z2).sum(1)
        a["m11"][l] += (Z @ Z.T).astype(np.float64)
        a["m21"][l] += ((Z * Z) @ Z.T).astype(np.float64)
        Y = np.maximum(Z, 0.0)
    a["yfin"] += Y.astype(np.float64).sum(1)
    done += b
    if done % (B * 64) == 0:
        print(f"{done} samples, {time.time() - t0:.0f}s", flush=True)


def finish(a):
    c = a["c"]; m = a["s1"] / c; e2 = a["s2"] / c; e3 = a["s3"] / c; e4 = a["s4"] / c
    var = e2 - m * m
    k3 = e3 - 3 * m * e2 + 2 * m ** 3
    k4 = e4 - 4 * m * e3 + 6 * m * m * e2 - 3 * m ** 4 - 3 * var * var
    E11 = a["m11"] / c; E21 = a["m21"] / c                     # E[z_i z_c], E[z_i^2 z_c]
    # E[(z_i - m_i)^2 (z_c - m_c)] = E[z_i^2 z_c] - 2 m_i E[z_i z_c] - m_c E[z_i^2] + 2 m_i^2 m_c
    D21 = E21 - 2 * m[:, :, None] * E11 - m[:, None, :] * e2[:, :, None] + 2 * (m * m)[:, :, None] * m[:, None, :]
    return dict(n=c, mu=m.astype(np.float32), var=var.astype(np.float32), k3=k3.astype(np.float32),
                k4=k4.astype(np.float32), D21=D21.astype(np.float32), yfin=(a["yfin"] / c).astype(np.float32))


f0, f1 = finish(acc[0]), finish(acc[1])
full = {k: acc[0][k] + acc[1][k] for k in ("c", "s1", "s2", "s3", "s4", "m11", "m21", "yfin")}
ff = finish(full)
for tag, d in (("h0", f0), ("h1", f1), ("full", ff)):
    np.savez(f"{outp}_{tag}.npz", **d)
print(f"done {N} samples in {time.time() - t0:.0f}s", flush=True)
