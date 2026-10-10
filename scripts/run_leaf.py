"""Compile exact great-circle leaves of an official network.

  python scripts/run_leaf.py DATA NET SEED [--planes K] [--check]
Writes $OUT/leaf_{NET}_{SEED}.npz with g (K, n), events (K, L), wall residuals, closure flags, seconds.
--check: first run a small random network (n=24, L=3) and compare the leaf mean with brute-force quadrature.
"""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.leaf import compile_leaf, random_frame, c_radial


def brute(W, B, Q=200000):
    t = (np.arange(Q) + 0.5) * 2 * np.pi / Q
    X = B @ np.vstack([np.cos(t), np.sin(t)])
    H = X
    for Wl in W: H = np.maximum(Wl @ H, 0)
    return H.mean(axis=1)


if __name__ == "__main__":
    a = sys.argv[1:]
    if "--check" in a:
        rng = np.random.default_rng(0); n, L = 24, 3
        W = rng.standard_normal((L, n, n)) * np.sqrt(2 / n); B = random_frame(n, rng)
        r = compile_leaf(W, B); gb = brute(W, B)
        print(f"check n={n} L={L}: events {r['events'].tolist()}, max|g - brute| = {np.abs(r['g'] - gb).max():.2e}, "
              f"wall residual {r['wall_residual']:.2e}, closed {r['closed']}")
        a.remove("--check")
        if not a: sys.exit()
    D, net, seed = a[0], int(a[1]), int(a[2]); K = int(a[a.index("--planes") + 1]) if "--planes" in a else 1
    W = np.load(f"{D}/W_off{net}.npy").astype(np.float64); n, L = W.shape[1], W.shape[0]
    rng = np.random.default_rng(1000 * seed + net); G = []; E = []; R = []; C = []; T = []
    for k in range(K):
        B = random_frame(n, rng); t0 = time.time(); r = compile_leaf(W, B); dt = time.time() - t0
        G.append(r["g"]); E.append(r["events"]); R.append(r["wall_residual"]); C.append(r["closed"]); T.append(dt)
        print(f"net {net} seed {seed} plane {k}: events {int(r['events'].sum())} (per layer {r['events'].tolist()}), "
              f"wall residual {r['wall_residual']:.1e}, closed {r['closed']}, {dt:.0f}s, c_n g rms {np.sqrt(np.mean((c_radial(n) * r['g'])**2)):.4f}", flush=True)
    out = os.environ.get("OUT", ".")
    np.savez(f"{out}/leaf_{net}_{seed}.npz", g=np.array(G), events=np.array(E), wall=np.array(R), closed=np.array(C), secs=np.array(T))
