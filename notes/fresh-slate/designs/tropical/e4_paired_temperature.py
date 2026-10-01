"""E4: temperature bridge with a sampled residual (hybrid; labelled as such).
F(0) = G(T) + E[a_L(x;0) - a_L(x;T)] + e(T): the high-T closure G(T) plus a paired-sample estimate of the
tropical-to-finite difference.  Measures, per T: e(T) = mse(G(T) - F(T)) and Vd(T) = mean per-neuron variance of the
paired difference, and the projected final MSE with Ns paired samples: e(T) + Vd(T)/Ns (vs plain MC V0/(2 Ns)).
Usage: set i N"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../bench"))
import numpy as np, bench
from e3_temperature import closure, act
setname, i, N = sys.argv[1], int(sys.argv[2]), int(float(sys.argv[3]))
S = bench.load_set(setname); Ws = [w.astype(np.float64) for w in bench.weights(S, i)]; n = Ws[0].shape[0]
t = S["means"][i][-1]
Ts = [0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.2]
rng = np.random.default_rng(7 + i); B = 50000
s1 = np.zeros((len(Ts), n)); s2 = np.zeros((len(Ts), n)); f1 = np.zeros((len(Ts), n)); v0 = np.zeros(n); m0 = np.zeros(n)
for b in range(N // B):
    X = rng.standard_normal((B, n)); a0 = X
    for W in Ws: a0 = np.maximum(a0 @ W, 0)
    m0 += a0.sum(0); v0 += (a0 ** 2).sum(0)
    for k, T in enumerate(Ts):
        a = X
        for W in Ws: a = act(a @ W, T)
        d = a0 - a; s1[k] += d.sum(0); s2[k] += (d * d).sum(0); f1[k] += a.sum(0)
M = B * (N // B); V0 = (v0 / M - (m0 / M) ** 2).mean()
print(f"{setname} mlp{i}: plain per-sample variance V0 = {V0:.3e}")
for Ns in (3277, 6554):
    print(f"  Ns={Ns}: plain MC (2Ns samples) {V0/(2*Ns):.2e}")
for k, T in enumerate(Ts):
    G = closure(Ws, T)[-1]; FT = f1[k] / M; Vd = (s2[k] / M - (s1[k] / M) ** 2).mean()
    e = ((G - FT) ** 2).mean()
    est = G + s1[k] / M   # with M samples (diagnostic): should approach truth
    print(f"  T={T:4.2f}  e(T)={e:.2e}  Vd={Vd:.2e}  (Vd/V0={Vd/V0:.3f})  proj mse @3277 pairs {e+Vd/3277:.2e}   check G+mean diff vs truth {((est-t)**2).mean():.2e}")
