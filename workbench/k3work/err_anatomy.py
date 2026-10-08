# Where does the chain's output error live? Per-neuron error of the final mean against features of the last layer's
# pre-activation state (alpha = mu/sigma, sigma, the carried kappa3 diagonal, the gate saturation), from v29err_off{net}.pkl.
import numpy as np, pickle, sys
from scipy.special import ndtr
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
P = pickle.load(open(f"v29err_off{net}.pkl", "rb")); out, mt = P["out"], P["truth"]; D = {d["layer"]: d for d in P["dumps"]}
err = out[-1] - mt[-1]; n = len(err); mse = np.mean(err**2)
print(f"net {net}: raw MSE {mse:.4e}; rms err {np.sqrt(mse):.3e}; |truth| rms {np.sqrt(np.mean(mt[-1]**2)):.3f}; mean err {err.mean():+.3e} (coherent part of MSE {err.mean()**2/mse:.1%})")
d = D[max(D)]; print(f"state from layer {max(D)}"); mu, var, D3 = d["mu"], d["var"], d["D3"]; sig = np.sqrt(var); a = mu / sig; Ph = ndtr(a); ph = np.exp(-a*a/2)/np.sqrt(2*np.pi); m = sig*(a*Ph+ph)
print(f"last-layer pre-activation: alpha quantiles 1/10/50/90/99% {np.quantile(a,[.01,.1,.5,.9,.99]).round(2)}")
def bins(x, name, edges):
    print(f"  by {name}:")
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (x >= lo) & (x < hi)
        if s.sum() == 0: continue
        print(f"    [{lo:+6.2f},{hi:+6.2f}) n={s.sum():4d}  share of MSE {np.sum(err[s]**2)/np.sum(err**2):6.1%}  rms err {np.sqrt(np.mean(err[s]**2)):.2e}  mean err {err[s].mean():+.2e}  rms truth {np.sqrt(np.mean(mt[-1][s]**2)):.3f}")
bins(a, "alpha", [-9, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 99])
bins(sig / sig.mean(), "sigma / mean sigma", [0, 0.5, 0.75, 1, 1.25, 1.5, 2, 9])
rel = err / np.maximum(np.abs(mt[-1]), 1e-3)
print(f"  relative error rms {np.sqrt(np.mean(rel**2)):.2e}; corr(err, truth) {np.corrcoef(err, mt[-1])[0,1]:+.3f}; corr(err, m_gauss - truth) {np.corrcoef(err, m - mt[-1])[0,1]:+.3f}; corr(err, D3) {np.corrcoef(err, D3)[0,1]:+.3f}; corr(err, alpha) {np.corrcoef(err, a)[0,1]:+.3f}")
# the Gaussian closure's own error at the last step (from the chain's pre-activation state) vs the chain's error
e_g = m - mt[-1]; print(f"  last-step Gaussian mean from the chain's (mu, sigma): MSE {np.mean(e_g**2):.3e}; the chain's correction (out - m): rms {np.sqrt(np.mean((out[-1]-m)**2)):.3e}; corr(correction, needed) {np.corrcoef(out[-1]-m, mt[-1]-m)[0,1]:+.3f}; slope {np.sum((out[-1]-m)*(mt[-1]-m))/np.sum((mt[-1]-m)**2):.3f}")
# projection of the error on the mean direction and on the top eigenvector of the last covariance
C = d["C_off"] + np.diag(var); ev, V = np.linalg.eigh(C); u1 = V[:, -1]; mh = mt[-1] / np.linalg.norm(mt[-1])
print(f"  error along the output mean direction: {(err @ mh)**2 / np.sum(err**2):.1%} of MSE; along the top eigenvector of the last pre-activation covariance: {(err @ u1)**2/np.sum(err**2):.1%}; top-16 eigenvectors: {np.sum((V[:, -16:].T @ err)**2)/np.sum(err**2):.1%}")
