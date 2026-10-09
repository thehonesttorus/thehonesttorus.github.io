"""Exact (coordinate) Laplacian of a chain estimator on an official network, split by coordinate blocks.
e(y) := chain output when the input law is N(y, I);  d2_i = (e(h e_i) + e(-h e_i) - 2 e(0)) / h^2.
  python scripts/lap_exact.py net=0 kind=gauss|k3v3 coords=0-31 h=0.02 [chain opts]
writes $OUT/lapx_<kind>_<net>_c<a>-<b>.npz with e0, d2 (block, n), coords, h."""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.k3chain3 import k3_chain3
from whest.relu_gauss import relu_moments
from whest.kprop3 import kprop3_chain
from whest.kprop3c import kprop3c_chain
kw = dict(net=0, kind="gauss", coords="0-31", h=0.02, K=8); opts = {}
for a in sys.argv[1:]:
    k, v = a.split("=")
    if k in kw: kw[k] = type(kw[k])(v)
    else: opts[k] = v if k in ("k4", "k22gate", "oldmode", "tier", "tag") else (float(v) if "." in v else int(v))
D = os.environ.get("DATA", "data"); OUT = os.environ.get("OUT", "."); net, kind, h = kw["net"], kw["kind"], kw["h"]
c0, c1 = [int(x) for x in kw["coords"].split("-")]; coords = list(range(c0, c1 + 1))
W = np.load(f"{D}/W_off{net}.npy").astype(np.float64); L, n, _ = W.shape
truth = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)[-1]

def gauss_chain_m(y):
    mu = y.copy(); C = np.eye(n)
    for l in range(L):
        mz = W[l] @ mu; Cz = W[l] @ C @ W[l].T; mu, C, _, _ = relu_moments(mz, Cz, kw["K"])
    return mu

def chain(y):
    if kind == "gauss": return gauss_chain_m(y)
    if kind == "kprop3": return kprop3_chain(W, radial=bool(opts.get("radial", 1)), m0=y)[-1]
    if kind == "kprop3c": return kprop3c_chain(W, {k: v for k, v in opts.items() if k != "tag"}, m0=y)[-1]
    out, _ = k3_chain3(W, dict(opts, m0=y)); return out[-1]

t0 = time.time(); e0 = chain(np.zeros(n)); print(f"net {net} {kind}: e(0) MSE {np.mean((e0-truth)**2):.4e} ({time.time()-t0:.0f}s/run)", flush=True)
d2 = np.zeros((len(coords), n))
for i, c in enumerate(coords):
    y = np.zeros(n); y[c] = h
    d2[i] = (chain(y) + chain(-y) - 2 * e0) / h ** 2
    if i % 8 == 7: print(f"  coord {c} ({time.time()-t0:.0f}s)", flush=True)
np.savez(f"{OUT}/lapx_{kind}{opts.get('tag', '')}_{net}_c{c0}-{c1}.npz", e0=e0, d2=d2, coords=np.array(coords), h=h, truth=truth, opts=str(opts))
