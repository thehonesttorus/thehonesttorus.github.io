import numpy as np, math
Phi = lambda x: 0.5 * (1 + np.vectorize(math.erf)(x / math.sqrt(2))); phi = lambda x: np.exp(-x * x / 2) / math.sqrt(2 * math.pi)
def gate(mu, var):
    s = np.sqrt(var); a = mu / s; Ph = Phi(a); ph = phi(a); m = mu * Ph + s * ph
    return Ph, ph / s, m
def kappa_formula(mu, var, C, a, b, c):
    Ph, rho, m = gate(mu, var); e2 = 2 * Ph - 2 * m * rho - 2 * Ph * Ph; f2 = 2 * m * (1 - Ph)
    return (e2[a] * Ph[b] * Ph[c] * C[a, b] * C[a, c] + f2[a] * (rho[b] * Ph[c] * C[a, b] * C[b, c] + Ph[b] * rho[c] * C[a, c] * C[b, c]))
# (a) per-entry check with common random numbers: kappa(y0,y0,y1,y2) at +-eps
mu = np.array([0.3, -0.5, 0.8]); var = np.array([1.2, 0.7, 1.5]); C = np.array([[0, 0.6, -0.4], [0.6, 0, 0.5], [-0.4, 0.5, 0]])
rng = np.random.default_rng(0); U = rng.standard_normal((3, 40_000_000))
def kap(eps):
    S = np.diag(var) + eps * C; Lc = np.linalg.cholesky(S); Z = Lc @ U + mu[:, None]; Y = np.maximum(Z, 0)
    Y -= Y.mean(1, keepdims=True); y0, y1, y2 = Y
    return (y0 * y0 * y1 * y2).mean() - (y0 * y0).mean() * (y1 * y2).mean() - 2 * (y0 * y1).mean() * (y0 * y2).mean()
for eps in (0.1, 0.2):
    k2 = (kap(eps) + kap(-eps) - 2 * kap(0.0)) / (2 * eps * eps)
    print(f"eps {eps}: MC second-order coefficient {k2:.5f}  formula {kappa_formula(mu, var, C, 0, 1, 2):.5f}")
# (b) transport check: code's Cd vs brute force
n = 7; rng = np.random.default_rng(3); W = rng.standard_normal((n, n)) * np.sqrt(2 / n); H = W * W
mu = rng.standard_normal(n) * 0.5; var = 0.5 + rng.random(n); A = rng.standard_normal((n, n)) * 0.3; Cz = A + A.T; np.fill_diagonal(Cz, 0)
src = open("yclasses.py").read(); frag = src[src.index("    Cz = sym0(F[\"cov\"][l])"):src.index("    def fit(")]
frag = "\n".join(x[4:] for x in frag.splitlines()).replace('Cz = sym0(F["cov"][l]); ', '')
g = dict(zip(("Ph", "dl", "m"), gate(mu, var)))
exec(frag)
bf = np.zeros(n)
for i in range(n):
    for a in range(n):
        for b in range(n):
            for c in range(n):
                if len({a, b, c}) == 3:
                    bf[i] += 6 * H[i, a] * W[i, b] * W[i, c] * kappa_formula(mu, var, Cz, a, b, c)
print("transport: max rel diff", np.abs(Cd - bf).max() / np.abs(bf).max())
