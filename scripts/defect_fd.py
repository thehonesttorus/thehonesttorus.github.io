"""Downstream heat-equation defects of the Gaussian closure, by exact finite differences.

Feed the Gaussian N(mu_k, C_k) of the closure at layer k (k = 0: the input N(t, I) at t = 0) and view the downstream
closure's output e as a function of the fed mean mu. The exact downstream mean Phi(mu, C) = E_{N(mu,C)} F_k satisfies,
by homogeneity and Stein,   Phi = mu . grad_mu Phi + tr(C hess_mu Phi).   The downstream defect
    delta_k = e - mu . grad_mu e - tr(C hess_mu e)
vanishes for exact downstream computation; it tests only errors born after layer k.

  python scripts/defect_fd.py DATA NET K CHUNK NCHUNKS     (partial sums over directions; pool with --merge)
  python scripts/defect_fd.py --merge DIR DATA NET K
Directions: eigenvectors u_i of C_k scaled by sqrt(lambda_i) (k = 0: the coordinate axes); second differences with step
h = 1e-2 of each direction's scale. The chunk with index 0 also evaluates the radial term mu . grad_mu e.
"""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.relu_gauss import relu_moments

H = 1e-2


def downstream(W, k, mu, C, K=10):
    """Means of the activations of layers k..L-1 (0-based) when N(mu, C) is fed as the pre-activation of layer k."""
    out = []
    for l in range(k, len(W)):
        m, Kh, _, _ = relu_moments(mu, C, K); out.append(m)
        if l + 1 < len(W):
            mu = W[l + 1] @ m; C = W[l + 1] @ Kh @ W[l + 1].T
    return np.array(out)


def state(W, k):
    """Gaussian-closure pre-activation mean and covariance at layer k (k = 0: layer 0's pre-activation N(0, W0 W0^T))."""
    mu = np.zeros(W[0].shape[0]); C = W[0] @ W[0].T
    for l in range(k):
        m, Kh, _, _ = relu_moments(mu, C, 10); mu = W[l + 1] @ m; C = W[l + 1] @ Kh @ W[l + 1].T
    return mu, C


def directions(W, k, C):
    if k == 0:                       # input tilt t (metric I) acts on layer 0's mean as W0 t: directions = columns of W0
        return W[0]
    lam, U = np.linalg.eigh(C); return U * np.sqrt(np.clip(lam, 0, None))[None, :]


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "--merge":
        R, D, net, k = a[1], a[2], int(a[3]), int(a[4])
        parts = sorted(glob.glob(f"{R}/dfd_{net}_{k}_*.npz")); z0 = np.load(parts[0])
        lap = sum(np.load(p)["lap"] for p in parts); ndir = sum(int(np.load(p)["ndir"]) for p in parts)
        e = z0["e"]; rad = next(np.load(p)["rad"] for p in parts if np.load(p)["has_rad"])
        delta = e - rad - lap
        truth = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)[k:]
        err = e - truth
        print(f"net {net} k={k}: {len(parts)} chunks, {ndir} directions")
        for j in (0, len(e) // 2, len(e) - 1):
            el, dl = err[j], delta[j]; c = np.corrcoef(el, dl)[0, 1]; ast = (el @ dl) / (dl @ dl)
            f = lambda aa: np.mean((el - aa * dl) ** 2)
            print(f"  layer {k + j:2d}: MSE {np.mean(el**2):.3e} | defect rms {np.sqrt(np.mean(dl**2)):.2e} corr {c:+.3f} a* {ast:.3f} "
                  f"MSE(a*) {f(ast):.3e} MSE(1/2) {f(.5):.3e} MSE(1/3) {f(1/3):.3e}")
        np.savez(f"{R}/defect_{net}_{k}.npz", e=e, delta=delta, lap=lap, rad=rad)
        sys.exit()
    D, net, k, chunk, nch = a[0], int(a[1]), int(a[2]), int(a[3]), int(a[4])
    W = [np.ascontiguousarray(w, dtype=np.float64) for w in np.load(f"{D}/W_off{net}.npy")]
    mu, C = state(W, k); U = directions(W, k, C)
    e = downstream(W, k, mu, C)
    idx = np.array_split(np.arange(U.shape[1]), nch)[chunk]
    lap = np.zeros_like(e)
    for i in idx:
        d = H * U[:, i]
        lap += (downstream(W, k, mu + d, C) + downstream(W, k, mu - d, C) - 2 * e) / H ** 2
    rad = np.zeros_like(e); has = chunk == 0 and k > 0
    if has:
        rad = (downstream(W, k, mu * (1 + H), C) - downstream(W, k, mu * (1 - H), C)) / (2 * H)
    np.savez(f"{os.environ.get('OUT', '.')}/dfd_{net}_{k}_{chunk}.npz", e=e, lap=lap, rad=rad, has_rad=(has or k == 0), ndir=len(idx))
    print(f"net {net} k={k} chunk {chunk}/{nch}: {len(idx)} directions done", flush=True)
