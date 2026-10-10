"""The diagonal-cumulant chain: the exact mean recursion t_l = W_l m_(l-1) with the four-cumulant mean map at every
layer, against the Gaussian closure.  python scripts/kik_chain.py RESULTS DATA NET [NET ...]

State per layer: pre-activation mean mu and covariance C (the closure's Gaussian bivariate Hermite series for the
off-diagonal post-activation covariance, K = 10), plus each neuron's standardised kappa_3, kappa_4.  The mean and the
diagonal second moment of relu(z) use the Gram-Charlier density phi (1 + k3/6 He3 + k4/24 He4 + k3^2/72 He6):
  E z_+   = sigma (d0 + k3/6 d3 + k4/24 d4 + k3^2/72 d6),  d_k = E[(U + a)_+ He_k] = He_{k-2}(-a) phi(a)  (k >= 2)
  E z_+^2 = sigma^2 (e0 + k3/6 e3 + k4/24 e4 + k3^2/72 e6), e_k = 2 He_{k-3}(-a) phi(a)               (k >= 3)
Here the kappas are ORACLE inputs (the pooled Monte Carlo cumulants of each layer's pre-activations): the run
measures what a compiler that delivers only the n diagonal kappas per layer would achieve, and how the remaining
error splits between the mean map and the covariance."""
import sys, os, numpy as np
from scipy.special import ndtr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scripts.kik_merge import pool_net

SQ = np.sqrt(2 * np.pi)


def hermite_relu(al, K):
    from whest.relu_gauss import hermite_relu as h
    return h(al, K)


def step(mu, C, k3, k4, K=10, mode="h4"):
    var = np.clip(np.diag(C), 1e-30, None); sd = np.sqrt(var); a = mu / sd; ph = np.exp(-a * a / 2) / SQ; P = ndtr(a)
    d = hermite_relu(a, K)
    rho = C / np.outer(sd, sd); np.fill_diagonal(rho, 0.0)
    Kh = np.zeros_like(C); rk = np.ones_like(rho); fact = 1.0
    for k in range(1, K + 1):
        rk = rk * rho; fact *= k; Kh += np.outer(d[k], d[k]) * rk / fact
    Kh *= np.outer(sd, sd)
    m = sd * d[0]; s2 = var * ((1 + a * a) * P + a * ph)
    if mode in ("h3", "h4"):
        c3 = k3; c4 = k4 if mode == "h4" else 0 * k4; c33 = c3 * c3 if mode == "h4" else 0 * c3
        d3, d4, d6 = -a * ph, (a * a - 1) * ph, (a ** 4 - 6 * a * a + 3) * ph
        e3, e4, e6 = 2 * ph, -2 * a * ph, 2 * (3 * a - a ** 3) * ph
        m = sd * (d[0] + c3 / 6 * d3 + c4 / 24 * d4 + c33 / 72 * d6)
        s2 = var * (((1 + a * a) * P + a * ph) + c3 / 6 * e3 + c4 / 24 * e4 + c33 / 72 * e6)
    np.fill_diagonal(Kh, s2 - m * m)
    return m, Kh


def kappas(a, l, N):
    L = a["pos"].shape[0]; P = (a["pw_low"][l] if l < L - 1 else a["pw_top"]) / N
    t, q, m3, m4 = P[:, 0], P[:, 1], P[:, 2], P[:, 3]; sd = np.sqrt(q - t * t)
    return (m3 - 3 * t * q + 2 * t ** 3) / sd ** 3, (m4 - 4 * t * m3 + 6 * t * t * q - 3 * t ** 4) / sd ** 4 - 3, t, sd * sd


def run(W, a, N, mode):
    mu = np.zeros(W[0].shape[0]); C = W[0] @ W[0].T; out, terr, verr = [], [], []
    for l in range(len(W)):
        k3, k4, tm, vm = kappas(a, l, N)
        terr.append(np.sqrt(np.mean((mu - tm) ** 2))); verr.append(np.sqrt(np.mean((np.diag(C) - vm) ** 2)))
        m, Kh = step(mu, C, k3, k4, mode=mode); out.append(m)
        if l + 1 < len(W):
            mu = W[l + 1] @ m; C = W[l + 1] @ Kh @ W[l + 1].T
    return np.array(out), np.array(terr), np.array(verr)


if __name__ == "__main__":
    R, D = sys.argv[1], sys.argv[2]; rms = lambda x: np.sqrt(np.mean(x ** 2))
    for net in map(int, sys.argv[3:]):
        a = pool_net(R, net); N = a["N"]; truth = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
        W = [np.asarray(w, dtype=np.float64) for w in np.load(f"{D}/W_off{net}.npy")]
        line = [f"net {net}:"]
        for mode in ("gauss", "h3", "h4"):
            o, te, ve = run(W, a, N, mode)
            e = [rms(o[l] - truth[l]) for l in range(len(W))]
            line.append(f"{mode}: final {e[-1]:.2e} (layers 1/4/8/12: {e[1]:.1e} {e[4]:.1e} {e[8]:.1e} {e[12]:.1e}; final t err {te[-1]:.1e}, var err {ve[-1]:.1e})")
            np.save(f"{R}/kchain_{net}_{mode}.npy", o)
        print("  ".join(line[:1]) + "\n  " + "\n  ".join(line[1:]), flush=True)
