"""The terminal ellipsoid of the companion's eq. (44) on the real networks: its visible spectrum and where errors sit.

  python scripts/s24_ellipsoid.py DATA OUT NET [NET ...]        (V56=dir with nXXX/out_X.npy, optional)

For full resolution (hidden states = inputs) the ellipsoid B = (1/n) sum_i (f_i - mu_i) (x) (f_i - mu_i) has the nonzero
spectrum of the output covariance S / n (S = Cov F(X), X ~ N(0, I)); an error vector e = mu - b is resolved in the
eigenbasis of S. Reported: participation ratio of S, share of tr S in the top 1/4/16/64/256 eigenvectors and in the bottom
half, the share of |e|^2 in eigenvector bands for b = Gaussian closure, CC1, v56, the number of top eigenvectors needed to
hold 90/99% of |e|^2, and the Rayleigh ratio e^T S e / (|e|^2 tr S / n)."""
import sys, os, math, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest import gate_channel as gc
from whest.cc1 import CC1

BS, N = 4096, 1 << 16
BANDS = ((0, 16), (16, 64), (64, 256), (256, 512), (512, 1024))


def run(D, OUT, net, V56):
    W32 = [np.ascontiguousarray(w, dtype=np.float32) for w in np.load(f"{D}/W_off{net}.npy")]; W64 = [w.astype(np.float64) for w in W32]
    n = W32[-1].shape[0]; d = W32[0].shape[1]; mu = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)[-1]
    S1 = np.zeros(n); S2 = np.zeros((n, n))
    for b in range(N // BS):
        h = np.random.default_rng([41_000 + net, b]).standard_normal((d, BS), dtype=np.float32)
        for Wl in W32: h = np.maximum(Wl @ h, 0)
        S1 += h.sum(1, dtype=np.float64); S2 += h @ h.T
    m = S1 / N; S = S2 / N - np.outer(m, m); lam, U = np.linalg.eigh(S); lam, U = lam[::-1], U[:, ::-1]
    bs = {"closure": gc.closure_states(W64, keep_C=False)[-1]["m"], "CC1": CC1(Q=7).run(W64)[0][-1]}
    f56 = f"{V56}/n{net:03d}/out_{net}.npy" if V56 else ""
    if f56 and os.path.exists(f56): bs["v56"] = np.load(f56)[-1]
    pr = lam.sum() ** 2 / np.sum(lam ** 2)
    line = [f"net {net}: tr S/n {lam.sum() / n:.4f}, PR {pr:.1f}; trace share top 1/4/16/64/256: " +
            "/".join(f"{lam[:k].sum() / lam.sum():.3f}" for k in (1, 4, 16, 64, 256)) + f", bottom half {lam[n // 2:].sum() / lam.sum():.4f}"]
    for bn, bv in bs.items():
        e = mu - bv; c = (U.T @ e) ** 2; cs = np.cumsum(c) / c.sum()
        line.append(f"   b={bn:8s} rms {math.sqrt(np.mean(e ** 2)):.3e}: error share in bands " + " ".join(f"{c[a:b_].sum() / c.sum():.3f}" for a, b_ in BANDS) +
                    f" | top-k for 90%/99% of |e|^2: {int(np.searchsorted(cs, 0.9)) + 1}/{int(np.searchsorted(cs, 0.99)) + 1}"
                    f" | Rayleigh ratio {e @ S @ e / (e @ e) / (lam.sum() / n):.2f}")
    line.append("   variance share in the same bands: " + " ".join(f"{lam[a:b_].sum() / lam.sum():.3f}" for a, b_ in BANDS))
    txt = "\n".join(line); print(txt, flush=True)
    with open(f"{OUT}/s24_ellipsoid_net{net}.txt", "w") as f: f.write(txt + "\n")


if __name__ == "__main__":
    D, OUT = sys.argv[1], sys.argv[2]; V56 = os.environ.get("V56", "")
    for net in map(int, sys.argv[3:]): run(D, OUT, net, V56)
