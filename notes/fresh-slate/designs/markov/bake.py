"""Own truth bake (until notes/fresh-slate/bench/ exists).
He MLP, x ~ N(0, I), antithetic pairs (x, -x). Stores weights seed, per-layer means of a_l,
per-layer raw moments of z_l up to order 4 (for diagnostics), per-neuron standard error of the
final-layer mean.  Usage: python bake.py WIDTH DEPTH SEED NSAMP OUT.npz"""
import sys, time
import numpy as np, os
MOM = os.environ.get('MOM', '0') == '1'

def weights(width, depth, seed):
    rng = np.random.default_rng(seed)
    return np.stack([(rng.standard_normal((width, width)) * np.sqrt(2.0 / width)).astype(np.float32)
                     for _ in range(depth)])

def bake(width, depth, seed, nsamp, chunk=32768):
    W = weights(width, depth, seed)
    rng = np.random.default_rng(10_000 + seed)
    S = np.zeros((depth, width)); Z = np.zeros((depth, 4, width))
    fs = np.zeros(width); fs2 = np.zeros(width); done = 0
    while done < nsamp:
        x = rng.standard_normal((chunk, width), dtype=np.float32)
        x = np.concatenate([x, -x])
        a = x
        for l in range(depth):
            z = a @ W[l]
            a = np.maximum(z, 0)
            S[l] += a.sum(0, dtype=np.float64)
            if MOM:
                z64 = z.astype(np.float64); p = z64.copy()
                for k in range(4):
                    Z[l, k] += p.sum(0); p *= z64
        pair = (a[:chunk].astype(np.float64) + a[chunk:]) / 2
        fs += pair.sum(0); fs2 += (pair ** 2).sum(0)
        done += 2 * chunk
    npair = done / 2
    se2 = (fs2 / npair - (fs / npair) ** 2) / npair
    return dict(W_seed=seed, means=S / done, zmom=Z / done, se2=se2, n=done)

if __name__ == "__main__":
    w, d, seed, n, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(float(sys.argv[4])), sys.argv[5]
    t = time.time(); r = bake(w, d, seed, n)
    np.savez(out, **r); print(out, "done", time.time() - t, "s, truth noise", r["se2"].mean())
