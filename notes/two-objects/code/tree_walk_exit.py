# Object 1 -> object 2: nearest-neighbour walk on the binary code tree with backtracking probability q(depth)
# and forward probabilities (1-q) * mu(vb)/mu(v). Claim: its exit (harmonic) measure on the Cantor set is mu,
# because the probability of escaping from a child back to its parent depends only on depth.
# Truncated at depth D: exit = first hitting of depth D (absorbing), which is the exit measure of the depth-D cylinders.
import numpy as np
rng = np.random.default_rng(7)
D = 6
cond = {(): None}
mu = {(): 1.0}
for k in range(D):
    for v in [p for p in mu if len(p) == k]:
        a = rng.uniform(0.1, 0.9); mu[v + (0,)] = mu[v]*a; mu[v + (1,)] = mu[v]*(1 - a)
q = lambda depth: 0.0 if depth == 0 else 0.3 + 0.1*np.sin(depth)        # depth-dependent backtracking, < 1/2
R = 400000
counts = {}
for r in range(R):
    v = ()
    while len(v) < D:
        if rng.random() < q(len(v)):
            v = v[:-1]
        else:
            p0 = mu[v + (0,)]/mu[v]
            v = v + ((0,) if rng.random() < p0 else (1,))
    counts[v] = counts.get(v, 0) + 1
leaves = [p for p in mu if len(p) == D]
emp = np.array([counts.get(p, 0)/R for p in leaves]); th = np.array([mu[p] for p in leaves])
z = (emp - th)/np.sqrt(th*(1 - th)/R)
print(f"{len(leaves)} cylinders at depth {D}: max |empirical - mu| = {np.abs(emp - th).max():.2e}, max |z-score| = {np.abs(z).max():.2f} (chi2/dof = {np.mean(z**2):.2f})")
