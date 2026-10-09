"""Heat-defect (Laplacian) probes of a chain estimator on an official network.
e(y) := chain output (final-layer means) when the input law is N(y, I).  Hutchinson probes of the mean-Laplacian:
d2_s = (e(h z_s) + e(-h z_s) - 2 e(0)) / h^2 with z_s ~ N(0, I), so that mean_s d2_s -> Lap_y e(0) and the defect
delta(0) = e(0) - Lap_y e(0) is the quantity the heat-flow theorem pairs with the chain's error.
  python scripts/lap_probe.py net=0 kind=gauss|k3v3 seeds=0-15 h=0.02 [chain opts]
writes $OUT/lap_<kind>_<net>_s<a>-<b>.npz with e0, d2 (S, n), seeds, h and the raw e(+-h z)."""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.k3chain3 import k3_chain3
from whest.relu_gauss import relu_moments
kw = dict(net=0, kind="gauss", seeds="0-15", h=0.02, K=8); opts = {}
for a in sys.argv[1:]:
    k, v = a.split("=")
    if k in kw: kw[k] = type(kw[k])(v)
    else: opts[k] = v if k in ("k4", "k22gate") else (float(v) if "." in v else int(v))
D = os.environ.get("DATA", "data"); OUT = os.environ.get("OUT", "."); net, kind, h = kw["net"], kw["kind"], kw["h"]
s0, s1 = [int(x) for x in kw["seeds"].split("-")]; seeds = list(range(s0, s1 + 1))
W = np.load(f"{D}/W_off{net}.npy").astype(np.float64); L, n, _ = W.shape
truth = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)[-1]

def gauss_chain_m(y):
    mu = y.copy(); C = np.eye(n)
    for l in range(L):
        mz = W[l] @ mu; Cz = W[l] @ C @ W[l].T; mu, C, _, _ = relu_moments(mz, Cz, kw["K"])
    return mu

def chain(y):
    if kind == "gauss": return gauss_chain_m(y)
    out, _ = k3_chain3(W, dict(opts, m0=y)); return out[-1]

t0 = time.time(); e0 = chain(np.zeros(n)); err = e0 - truth
print(f"net {net} {kind} {opts}: e(0) MSE {np.mean(err**2):.4e} ({time.time()-t0:.0f}s/run)", flush=True)
d2 = np.zeros((len(seeds), n)); ep = np.zeros_like(d2); em = np.zeros_like(d2)
for i, s in enumerate(seeds):
    z = np.random.default_rng(1000 + s).standard_normal(n)
    ep[i] = chain(h * z); em[i] = chain(-h * z); d2[i] = (ep[i] + em[i] - 2 * e0) / h ** 2
    lap = d2[: i + 1].mean(0); delta = e0 - lap; a = np.dot(delta, err) / np.dot(delta, delta)
    print(f"  seed {s} ({time.time()-t0:.0f}s): {i+1} probes  corr(err,delta) {np.corrcoef(err, delta)[0,1]:+.3f}  a {a:+.3f}  "
          f"merge MSE {np.mean((err - a*delta)**2):.4e}  (raw {np.mean(err**2):.4e})", flush=True)
np.savez(f"{OUT}/lap_{kind}_{net}_s{s0}-{s1}.npz", e0=e0, d2=d2, ep=ep, em=em, seeds=np.array(seeds), h=h, truth=truth, opts=str(opts))
