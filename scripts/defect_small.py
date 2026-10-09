"""Heat-flow (Euler-Stein) defect on a tiny network: e(y) = chain at input law N(y, I); delta(0) = e(0) - Lap_y e(0)
(exact answers satisfy e = Lap e at y = 0 by homogeneity + Stein), per neuron, against the true error; and the
distribution of the chain's excess second variation over input directions.
  python scripts/defect_small.py n=256 L=4 N=10000000 h=0.02 [chain opts]"""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.k3chain3 import k3_chain3
from whest.relu_gauss import relu_moments
kw = dict(n=256, L=4, N=10_000_000, seed=0, h=0.02); opts = {}
for a in sys.argv[1:]:
    k, v = a.split("=")
    if k in kw: kw[k] = type(kw[k])(v)
    else: opts[k] = v if k in ("k4", "k22gate") else (float(v) if "." in v else int(v))
n, L, N, seed, h = kw["n"], kw["L"], kw["N"], kw["seed"], kw["h"]
rng = np.random.default_rng(seed); W = rng.standard_normal((L, n, n)) * np.sqrt(2.0 / n)
truth = np.zeros((L, n)); B = 100_000; t0 = time.time()
for _ in range(N // B):
    hh = rng.standard_normal((B, n))
    for l in range(L):
        hh = np.maximum(hh @ W[l].T, 0); truth[l] += hh.sum(0)
truth /= N; print(f"MC truth: {N} samples in {time.time()-t0:.0f}s, noise ~{np.sqrt(0.34/N):.1e}", flush=True)
def gauss_chain_m(y):
    mu = y.copy(); C = np.eye(n)
    for l in range(L):
        mz = W[l] @ mu; Cz = W[l] @ C @ W[l].T; mu, C, _, _ = relu_moments(mz, Cz, 8)
    return mu
def chain(y, kind):
    if kind == "gauss": return gauss_chain_m(y)
    out, _ = k3_chain3(W, dict(opts, m0=y)); return out[-1]
v_coh = W[0].T @ np.ones(n); v_coh /= np.linalg.norm(v_coh)
v_rnd = rng.standard_normal(n); v_rnd /= np.linalg.norm(v_rnd)
for kind in ("gauss", "k3v3"):
    t0 = time.time(); e0 = chain(np.zeros(n), kind); err = e0 - truth[-1]
    d2 = np.zeros((n, n))
    for i in range(n):
        ei = np.zeros(n); ei[i] = h
        d2[i] = (chain(ei, kind) + chain(-ei, kind) - 2 * e0) / h ** 2
    lap = d2.sum(0); delta = e0 - lap
    a = np.dot(delta, err) / np.dot(delta, delta); r2 = 1 - np.sum((err - a * delta) ** 2) / np.sum(err ** 2)
    print(f"[{kind}] ({time.time()-t0:.0f}s) MSE {np.mean(err**2):.3e} |err| {np.sqrt(np.mean(err**2)):.2e} |delta| {np.sqrt(np.mean(delta**2)):.2e} "
          f"corr(err,delta) {np.corrcoef(err, delta)[0,1]:+.3f} a {a:+.3f} MSE explained {r2:.3f} merge MSE {np.mean((err - a*delta)**2):.3e} "
          f"midpoint (e+lap)/2 MSE {np.mean(((e0+lap)/2 - truth[-1])**2):.3e}", flush=True)
    ex = d2 - lap[None, :] / n
    u, sv, vt = np.linalg.svd(ex, full_matrices=False)
    print(f"        excess second variation over directions: top sv " + " ".join(f"{x:.2e}" for x in sv[:5]) + f"; top-1 share {sv[0]**2/np.sum(sv**2):.3f}, top-4 {np.sum(sv[:4]**2)/np.sum(sv**2):.3f}; "
          f"|top dir . coherent W1^T 1| {abs(np.dot(u[:,0], v_coh)):.3f}")
    for name, v in (("coherent", v_coh), ("random", v_rnd), ("top-excess", u[:, 0])):
        d2v = (chain(h * v, kind) + chain(-h * v, kind) - 2 * e0) / h ** 2; exv = d2v - lap / n
        print(f"        dir {name:11s}: |d2_v| {np.sqrt(np.mean(d2v**2)):.2e} |excess| {np.sqrt(np.mean(exv**2)):.2e} corr(excess, err) {np.corrcoef(exv, err)[0,1]:+.3f} corr(e0 - n d2_v, err) {np.corrcoef(e0 - n*d2v, err)[0,1]:+.3f}", flush=True)
