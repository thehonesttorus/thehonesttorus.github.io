"""Exact checks on small bias-free ReLU networks (x ~ N(0, I_d)) of the arrow identities of Stage 24.

(1) Kink-sum formula. F(x) = a_L(x) is piecewise linear and positively homogeneous, so grad F = M(x) and, by Stein's
    lemma twice, E F = E[Delta F] with Delta F = sum_(l,i) delta(z_(l,i)) |grad z_(l,i)|^2 A_l(x) e_i, A_l = downstream
    gated product. Equivalently mu = sum_(l,i) int_{z_(l,i) = 0} |grad z_(l,i)| A_l e_i dsigma_gamma.
    Checked by a window estimator delta_eps(z) = 1{|z| < eps}/(2 eps) against direct Monte Carlo of E F.
(2) The closure's mean recursion unrolls to the mean-field kink sum: m_L = sum_l T_(L<-l) (sd_l * phi(alpha_l)).
(3) Arrow weights at layer 1: the stationary flux of the Mehler rotation x_t = x cos t + y sin t through the facet
    between patterns S and S + {i} equals (1/2pi) P(sign pattern of z_(-i) | z_i = 0) (Schur complement), and the
    per-neuron flip rate is exactly 1/pi; checked by counting flips over a small rotation step.
  python notes/stage24/theory24/check_kink_sum.py"""
import sys, os, itertools, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))
from whest import gate_channel as gc

rng = np.random.default_rng(24)


def forward(Ws, X):
    """X (d, N). Returns pre-activations z_l (n, N) and post-activations."""
    zs, h = [], X
    for W in Ws: z = W @ h; zs.append(z); h = np.maximum(z, 0)
    return zs, h


def kink_sum(Ws, N=2_000_000, eps=2e-3, B=50_000, t=1e-3):
    """Two estimates of sum_(l,i) int_{z_(l,i)=0} |grad z| A_l e_i dsigma_gamma per layer l, and direct E F.
    window: E[delta_eps(z) |grad z|^2 A_l e_i]; flips: Rice's formula for the Mehler rotation x_t = x cos t + y sin t,
    E[sum over gate flips in [0, t] of |grad z| A_l e_i] = t sqrt(2/pi) int_{z=0} |grad z| A_l e_i dsigma_gamma."""
    L = len(Ws); n = Ws[-1].shape[0]; d = Ws[0].shape[1]
    kw = np.zeros((L, n)); kf = np.zeros((L, n)); EF = np.zeros(n); dead = 0
    for _ in range(N // B):
        X = rng.standard_normal((d, B)); Y = rng.standard_normal((d, B)); Xt = X * np.cos(t) + Y * np.sin(t)
        zs, h = forward(Ws, X); zt, _ = forward(Ws, Xt); EF += h.sum(1)
        gates = [(z > 0).astype(float) for z in zs]; flips = [(a > 0) != (b > 0) for a, b in zip(zs, zt)]
        dead += sum(int(np.sum(g.sum(0) == 0)) for g in gates[:-1])
        G = np.broadcast_to(Ws[0][None], (B,) + Ws[0].shape).copy()
        for l in range(L):
            if l > 0: G = np.einsum("ij,bjk->bik", Ws[l], gates[l - 1].T[:, :, None] * G)
            g2 = np.sum(G * G, 2).T                                                 # (n, B) |grad z_(l,i)|^2
            vw = (np.abs(zs[l]) < eps) / (2 * eps) * g2
            vf = flips[l] * np.sqrt(g2) / (t * np.sqrt(2 / np.pi))
            for k in range(l + 1, L): vw = gates[k] * (Ws[k] @ vw); vf = gates[k] * (Ws[k] @ vf)
            kw[l] += vw.sum(1); kf[l] += vf.sum(1)
    return kw / N, kf / N, EF / N, dead


if __name__ == "__main__":
    for (d, n, L) in ((24, 24, 3), (32, 32, 4)):
        Ws = [rng.standard_normal((n, d if l == 0 else n)) * np.sqrt(2.0 / (d if l == 0 else n)) for l in range(L)]
        kw, kf, EF, dead = kink_sum(Ws)
        rel = lambda a: np.max(np.abs(a - EF)) / np.max(np.abs(EF))
        print(f"[1] d=n={n}, L={L} (samples with a fully dead hidden layer: {dead}): rms E F = {np.sqrt(np.mean(EF ** 2)):.4f}; "
              f"kink sum, window estimator: max rel. dev {rel(kw.sum(0)):.2e}; flip estimator: max rel. dev {rel(kf.sum(0)):.2e}")
        print(f"    per-layer kink contributions |kappa_l| (flip estimator) = " + " ".join(f"{np.linalg.norm(k):.4f}" for k in kf) +
              f"; |E F| = {np.linalg.norm(EF):.4f}")
        # [2] closure mean = mean-field kink sum
        st = gc.closure_states(Ws, keep_C=False); m = st[-1]["m"]
        T = np.eye(n); acc = np.zeros(n)
        for l in range(L - 1, -1, -1):
            acc += T @ (st[l]["sd"] * np.exp(-0.5 * st[l]["alpha"] ** 2) / np.sqrt(2 * np.pi))
            T = (T * st[l]["Pa"][None, :]) @ Ws[l] if l > 0 else T
        print(f"[2] closure m_L = {np.array2string(m, precision=6)}\n    mean-field kink sum = {np.array2string(acc, precision=6)}; max |diff| {np.max(np.abs(m - acc)):.1e}")
    # [3] layer-1 arrow weights: flux through facets and flip rates under the Mehler rotation
    d = n = 5; W = rng.standard_normal((n, d)); C = W @ W.T
    t, N = 2e-3, 4_000_000
    X = rng.standard_normal((d, N)); Y = rng.standard_normal((d, N)); Xt = X * np.cos(t) + Y * np.sin(t)
    s0, s1 = (W @ X > 0), (W @ Xt > 0); flips = s0 != s1
    print(f"[3] layer-1 flip rate per neuron (theory 1/pi = {1 / np.pi:.4f}): " + " ".join(f"{r:.4f}" for r in flips.mean(1) / t))
    # facet flux for one (S, i): count transitions S -> S + {i} with exactly one flip, versus (1/2pi) P(pattern_{-i} | z_i = 0)
    i = 0; code0 = (s0 * (1 << np.arange(n))[:, None]).sum(0); one = flips.sum(0) == 1
    from scipy.stats import multivariate_normal
    Sc = C[1:, 1:] - np.outer(C[1:, 0], C[0, 1:]) / C[0, 0]
    out = []
    for rest in itertools.product([0, 1], repeat=n - 1):
        sg = np.array([1 if r else -1 for r in rest], float); Cs = Sc * np.outer(sg, sg)
        orth = multivariate_normal(mean=np.zeros(n - 1), cov=Cs, abseps=1e-7, releps=1e-6).cdf(np.zeros(n - 1) + 1e-12, lower_limit=np.full(n - 1, -np.inf)) if False else None
        # P(sg * z_rest > 0 | z_i = 0) = P(-sg * z_rest < 0) -> cdf of N(0, Cs) at 0 for the negated vector
        orth = multivariate_normal(mean=np.zeros(n - 1), cov=Cs).cdf(np.zeros(n - 1))
        S = sum((1 << (k + 1)) for k, r in enumerate(rest) if r)                    # pattern without i
        up = np.mean(one & flips[i] & (code0 == S)) / t                             # S -> S + {i}
        dn = np.mean(one & flips[i] & (code0 == S + 1)) / t                         # S + {i} -> S
        out.append((orth / (2 * np.pi), up, dn))
    out = np.array(out)
    print(f"    facet fluxes for neuron 0 over all {len(out)} rest-patterns: max |MC up - theory| {np.max(np.abs(out[:, 1] - out[:, 0])):.1e}, "
          f"max |MC down - theory| {np.max(np.abs(out[:, 2] - out[:, 0])):.1e} (theory values range {out[:, 0].min():.4f}..{out[:, 0].max():.4f}); "
          f"sum of theory fluxes {out[:, 0].sum():.5f} = 1/(2 pi) = {1 / (2 * np.pi):.5f}")
