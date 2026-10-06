# Which per-neuron term is the last-step closure missing? Regress the residual (truth - chain output) on derived
# per-neuron features of the last layer's pre-activation state (chain's own Edgeworth terms, the scale-mixture tail,
# higher Hermite terms), with and without the chain's own correction as a regressor.
import numpy as np, pickle, sys
from scipy.special import ndtr
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
P = pickle.load(open(f"v29err_off{net}.pkl", "rb")); out, mt = P["out"][-1], P["truth"][-1]; D = {d["layer"]: d for d in P["dumps"]}; d = D[15]
mu, var, D3, g4, pk1v, W_all = d["mu"], d["var"], d["D3"], d["g4row"], d["pk1v"], d["W_all"]
sig = np.sqrt(var); a = mu / sig; Ph = ndtr(a); ph = np.exp(-a*a/2)/np.sqrt(2*np.pi); m = sig*(a*Ph+ph)
r = mt - out; need = mt - m; corr_chain = out - m
WP = [(0,1),(0,2),(0,3),(0,4),(1,1),(1,2),(2,1),(2,2),(3,1),(3,2),(3,3),(3,4),(4,1),(4,2),(4,3),(4,4),(5,1),(5,2),(6,1),(6,2),(7,1)]
col = lambda p: W_all[:, WP.index(p)]
v = 1.5*mu*var; g = float(v@D3/(v@v)); print(f"net {net}: last layer g(D3) = {g:.4f}; residual rms {np.sqrt(np.mean(r**2)):.3e}; needed-correction rms {np.sqrt(np.mean(need**2)):.3e}; chain correction rms {np.sqrt(np.mean(corr_chain**2)):.3e}")
# Hermite coefficients of relu at (mu, sigma): E f^(k) sigma^k = c_k: c3 = -a ph, c4 = (a^2-1) ph, c5 = -(a^3-3a) ph, c6 = (a^4-6a^2+3) ph (He_{k-2}(a) (-1)^k phi), as sigma^(1-k) multiples
He = {0: np.ones_like(a), 1: a, 2: a*a-1, 3: a**3-3*a, 4: a**4-6*a*a+3, 5: a**5-10*a**3+15*a}
ck = lambda k: ((-1)**k) * He[k-2] * ph / sig**(k-1)   # E[f^(k)(h)] for k >= 2
feats = {
  "chain k3 term D3 c3/6": D3 * ck(3) / 6,
  "chain k4 term g4 c4/24": g4 * ck(4) / 24,
  "D3^2 c6/72": D3*D3*ck(6)/72 if 6 in He else None,
  "sm tail -(g/8) mu (Phi - a phi)": -(g/8)*mu*(Ph - a*ph),
  "sm k5 term: 1.5 g mu sig^2 * sig^2 c5 (shape)": g*mu*var*ck(5),
  "sm k6 term: g sig^6 c6 (shape)": g*var**3*ck(6),
  "m": m, "sigma phi": sig*ph, "mu Phi": mu*Ph, "D3 (raw)": D3, "D3 c5 sig^2": D3*ck(5)*var,
}
feats = {k: v for k, v in feats.items() if v is not None}
names = list(feats); X = np.stack([feats[k] for k in names], 1)
def ols(y, X, names):
    Xc = np.column_stack([X, np.ones(len(y))]); beta, *_ = np.linalg.lstsq(Xc, y, rcond=None); pred = Xc @ beta
    return 1 - np.sum((y-pred)**2)/np.sum(y**2), beta, pred
R2, beta, pred = ols(r, X, names)
print(f"residual explained by all features jointly: R^2 {R2:.3f}")
for k, b in zip(names, beta): 
    f = feats[k]; print(f"   {k:46s} coef {b:+.3f}  single-feature R^2 {ols(r, f[:,None], [k])[0]:.3f}  (feature rms {np.sqrt(np.mean(f**2)):.2e})")
print(f"needed correction explained by the chain's own two terms: R^2 {ols(need, np.stack([feats['chain k3 term D3 c3/6'], feats['chain k4 term g4 c4/24']],1), [])[0]:.3f}; by the chain's full correction: {ols(need, corr_chain[:,None], [])[0]:.3f}")
# residual against the error in the state: compare the chain's mu_pre against what the truth implies? not available; instead the Gaussian-closure-only error vs alpha bins is above.
