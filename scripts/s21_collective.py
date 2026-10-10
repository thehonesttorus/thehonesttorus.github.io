"""The collective tier of Stage 21, measured in one forward pass per batch.

  python scripts/s21_collective.py DATA NET TASK NBATCH        writes $OUT/s21c_{NET}_{TASK}.npz

Rows of W_l are i.i.d. N(0, s^2 I) (s^2 = 2/n) and independent of h_(l-1). Split a row as w = a mu_hat + w_perp with
mu_hat = E h_(l-1) / |E h_(l-1)| (the exact means of the truth file). For fixed x, w_perp . h(x) ~ N(0, s^2 r(x)^2),
r^2 = |h|^2 - c^2, c = mu_hat . h. Hence the coherent function of Prop 5.1 is exactly

    Ybar_l(a) = E[ m_(l, j) | a_j = a ] = E_x G(a c(x), s^2 r(x)^2),   G(m, v) = E relu(N(m, v)),

a functional of the joint law of the two collective scalars (c, r^2). Per selected layer l the task accumulates
sum_x G(a_j c(x), s^2 r(x)^2) for the n realised a_j = w_j . mu_hat, and for every layer the mixed power sums of
(c, r^2) up to total order 6 (i + 2j <= 6)."""
import sys, os, time, numpy as np
from scipy.special import ndtr

BS = 4096; YB_LAYERS = (1, 2, 4, 8, 12, 15)          # layers (0-based) whose means get the coherent function
MIX = [(i, j) for i in range(7) for j in range(4) if 0 < i + 2 * j <= 6]


def G(m, v):
    sd = np.sqrt(v); a = m / sd
    return m * ndtr(a) + sd * np.exp(-0.5 * a * a) / np.sqrt(2 * np.pi)


if __name__ == "__main__":
    D, net, task, nb = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    W = [np.ascontiguousarray(w, dtype=np.float32) for w in np.load(f"{D}/W_off{net}.npy")]; L, n = len(W), W[0].shape[0]
    truth = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64); s2 = 2.0 / n
    mh = [truth[l] / np.linalg.norm(truth[l]) for l in range(L)]                       # post-activation mean directions
    a = {l: W[l].astype(np.float64) @ mh[l - 1] for l in YB_LAYERS}                      # a_j for layer l's rows
    yb = {l: np.zeros(n) for l in YB_LAYERS}; mix = np.zeros((L, len(MIX)))
    rng = np.random.default_rng(7_700_003 * (net + 1) + task); t0 = time.time()
    for b in range(nb):
        h = rng.standard_normal((n, BS), dtype=np.float32); cs = []
        for l, Wl in enumerate(W):
            h = np.maximum(Wl @ h, 0); hd = h.astype(np.float64)
            c = mh[l] @ hd; r2 = np.maximum(np.einsum("ib,ib->b", hd, hd) - c * c, 0); cs.append((c, r2))
            for k, (i, j) in enumerate(MIX): mix[l, k] += np.sum(c ** i * r2 ** j)
        for l in YB_LAYERS:
            c, r2 = cs[l - 1]
            yb[l] += G(a[l][:, None] * c[None, :], s2 * r2[None, :]).sum(1)
    np.savez(f"{os.environ.get('OUT', '.')}/s21c_{net}_{task}.npz", N=nb * BS, mix=mix, mixidx=np.array(MIX),
             **{f"yb_{l}": yb[l] for l in YB_LAYERS}, **{f"a_{l}": a[l] for l in YB_LAYERS})
    print(f"net {net} task {task}: {nb * BS} samples in {time.time() - t0:.0f}s", flush=True)
