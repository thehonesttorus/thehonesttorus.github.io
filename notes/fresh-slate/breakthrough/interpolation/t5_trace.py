"""T5: Wiener-chaos split (in the column w_p) of the per-neuron third cumulant kappa3(z_{l+1,p}) = kappa3(a_l)[w_p^{x3}]:
chaos-1 part = 3 sigma^2 (t_l . w_p), t_l = E[|a~_l|^2 a~_l] (the 'trace channel', an n-vector); the rest is chaos 3.
MC at n=1024 on a bench network; centering with bake means (exact to noise)."""
import sys, json
import numpy as np
from common import bench
name, i, N = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
S = bench.load_set(name); W = bench.weights(S, i); T = S["means"][i]; L, n, _ = W.shape
sig2 = 2.0 / n
rng = np.random.default_rng(7); chunk = 8192
t = np.zeros((L, n)); m2 = np.zeros((L, n)); m3 = np.zeros((L, n)); m2a = np.zeros((L, n))
done = 0
while done < N:
    h = rng.standard_normal((chunk, n)).astype(np.float32)
    for l in range(L):
        mz = (T[l - 1] @ W[l]) if l > 0 else np.zeros(n)
        z = h @ W[l]
        y = (z - mz.astype(np.float32)).astype(np.float64)
        m2[l] += (y ** 2).sum(0); m3[l] += (y ** 3).sum(0)
        h = np.maximum(z, 0.0)
        at = h.astype(np.float64) - T[l]
        t[l] += ((at ** 2).sum(1)[:, None] * at).sum(0)
        m2a[l] += (at ** 2).sum(0)
    done += chunk
t /= done; m2 /= done; m3 /= done; m2a /= done
out = {}
for l in range(1, L):
    k3 = m3[l]; v = m2[l]
    c1 = 3 * sig2 * (t[l - 1] @ W[l].astype(np.float64))
    # noise of k3 estimate ~ sqrt(15 v^3 / N)
    noise = np.sqrt(15 * v ** 3 / done)
    r = dict(layer=l + 1, rms_k3=float(np.sqrt((k3 ** 2).mean())), rms_noise=float(np.sqrt((noise ** 2).mean())),
             rms_c1=float(np.sqrt((c1 ** 2).mean())), corr=float(np.corrcoef(k3, c1)[0, 1]),
             frac_explained=float(1 - ((k3 - c1) ** 2).mean() / (k3 ** 2).mean()),
             mean_k3=float(k3.mean()), mean_c1=float(c1.mean()))
    out[l + 1] = r
    print({k: (round(x, 4) if isinstance(x, float) else x) for k, x in r.items()}, flush=True)
np.savez(f"results/t5_{name}_{i}_N{done}.npz", t=t, m2=m2, m3=m3, m2a=m2a)
json.dump(out, open(f"results/t5_{name}_{i}_N{done}.json", "w"), indent=1)
