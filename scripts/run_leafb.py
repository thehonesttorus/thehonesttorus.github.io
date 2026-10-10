"""Exact great-circle leaves with the batched compiler: python scripts/run_leafb.py DATA NET SEED0 K
Leaf k uses the Haar frame drawn from default_rng(1000 * (SEED0 + k) + NET) (the same frames as scripts/run_leaf.py).
Writes $OUT/leafb_{NET}_{SEED0}.npz: g (K, L, n) leaf means of every layer's activations, walls (K, L), secs (K,)."""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.leaf import random_frame
from whest.leafb import compile_leaf_b

D, net, s0, K = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
W = [np.ascontiguousarray(w, dtype=np.float64) for w in np.load(f"{D}/W_off{net}.npy")]
G, Wl, T = [], [], []
for k in range(K):
    B = random_frame(W[0].shape[1], np.random.default_rng(1000 * (s0 + k) + net))
    t0 = time.time(); r = compile_leaf_b(W, B); T.append(time.time() - t0); G.append(r["g"]); Wl.append(r["walls"])
    print(f"net {net} leaf {s0 + k}: walls {int(r['walls'].sum())}, {T[-1]:.1f}s", flush=True)
np.savez(f"{os.environ.get('OUT', '.')}/leafb_{net}_{s0}.npz", g=np.array(G), walls=np.array(Wl), secs=np.array(T))
