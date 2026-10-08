# Does the first-order (integration-by-parts) closure reproduce the pair moments of a pure scale mixture z = G y?
# Exact: E r(z_a) r(z_c) = E[G^2] M(mu, S)_ac = M(mu, S)_ac (homogeneity). The chain sees the mixture's own moments
# (mu_z, S_z, kappa3, kappa4 slices) and applies first-order IBP. The residual is the second-order defect, compared with
# the size of each first-order slice term (especially the (3,1) slice).
import numpy as np, pickle, sys
sys.path.insert(0, "../ncgprob"); sys.path.insert(0, "../num12")
from pairvar import Layer, sm_slices, gauss_raw, cum_from_raw, mehler
from closure import relu_coeffs, relu2_coeffs
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
D = pickle.load(open(f"v29dump_off{net}.pkl", "rb"))
rng = np.random.default_rng(0); idx = rng.choice(1024, 256, replace=False)   # a 256-neuron sub-block (Mehler sums are O(n^2 K))
for d in D:
    mu, var, C_off, D3 = d["mu"][idx], d["var"][idx], d["C_off"][np.ix_(idx, idx)], d["D3"][idx]
    v = 1.5 * d["mu"] * d["var"]; g = float(v @ d["D3"] / (v @ v))
    S = C_off + np.diag(var); sig = np.sqrt(var)
    # exact pair raw moment of the mixture = Gaussian pair raw moment at (mu, S)
    A = relu_coeffs(mu, sig, 20); B = relu2_coeffs(mu, sig, 20)
    R = S / np.outer(sig, sig); exact11 = mehler(A, A, R, 0, 0, 14); exact21 = mehler(B, A, R, 0, 0, 14) * (1 + 3 * g / 8)
    # the mixture's own moments to first order in g (what the chain carries): raw moments scale by E G^r
    cr = {1: -1/8, 2: 0.0, 3: 3/8, 4: 1.0}; raw = gauss_raw(mu, S)
    s1 = raw["s1"] * (1 + cr[1] * g); S2 = raw["S2"]; M21 = raw["M21"] * (1 + cr[3] * g); M31 = raw["M31"] * (1 + g); M22 = raw["M22"] * (1 + g)
    s3 = raw["s3"] * (1 + cr[3] * g); s4 = raw["s4"] * (1 + g)
    cum = cum_from_raw(s1, S2, M21, M31, M22, s3, s4)
    mu_z, S_z = cum["mu"], cum["S"]; sig_z = np.sqrt(np.diag(S_z))
    Az = relu_coeffs(mu_z, sig_z, 20); Bz = relu2_coeffs(mu_z, sig_z, 20)
    Lz = Layer(mu_z, S_z, Az, Bz, K=14)
    # first-order variations of the raw pair moments (Layer._perturbed returns cumulant deltas; we need raw e11/e21 deltas) -> recompute directly
    k3, K21, k4, K31, K22 = cum["k3"], cum["K21"], cum["k4"], cum["K31"], cum["K22"]
    sigz = sig_z; Rz = Lz.R
    def dE3(P, Q):
        k3s = k3 / sigz**3; K21s = K21 / np.outer(sigz**2, sigz); K12s = K21s.T
        return (k3s[:, None] * mehler(P, Q, Rz, 3, 0, 14) + 3 * K21s * mehler(P, Q, Rz, 2, 1, 14) + 3 * K12s * mehler(P, Q, Rz, 1, 2, 14) + k3s[None, :] * mehler(P, Q, Rz, 0, 3, 14)) / 6
    def dE4(P, Q, use31=True, use22=True, use4=True):
        k4s = k4 / sigz**4; K31s = K31 / np.outer(sigz**3, sigz); K22s = K22 / np.outer(sigz**2, sigz**2); K13s = K31s.T
        out = np.zeros_like(Rz)
        if use4: out += k4s[:, None] * mehler(P, Q, Rz, 4, 0, 14) + k4s[None, :] * mehler(P, Q, Rz, 0, 4, 14)
        if use31: out += 4 * K31s * mehler(P, Q, Rz, 3, 1, 14) + 4 * K13s * mehler(P, Q, Rz, 1, 3, 14)
        if use22: out += 6 * K22s * mehler(P, Q, Rz, 2, 2, 14)
        return out / 24
    off = ~np.eye(len(idx), dtype=bool)
    for name, P, Q, ex in (("(1,1) raw E r r'", Az, Az, exact11), ("(2,1) raw E r^2 r'", Bz, Az, exact21)):
        G0 = mehler(P, Q, Rz, 0, 0, 14); t3 = dE3(P, Q); t4 = dE4(P, Q); t31 = dE4(P, Q, use22=False, use4=False); t22 = dE4(P, Q, use31=False, use4=False); t44 = dE4(P, Q, use31=False, use22=False)
        first = G0 + t3 + t4; res = ex - first
        nz = np.linalg.norm
        print(f"layer {d['layer']} g {g:.4f} {name}: |exact-Gaussian(mu_z,S_z)| {nz((ex-G0)[off]):.3e} = |k3 term| {nz(t3[off]):.3e} + |k4 term| {nz(t4[off]):.3e} [(3,1) {nz(t31[off]):.3e}, (2,2) {nz(t22[off]):.3e}, diag {nz(t44[off]):.3e}]; residual after first order {nz(res[off]):.3e} ({nz(res[off])/nz((ex-G0)[off]):.2%} of the first-order correction; vs (3,1) term {nz(res[off])/nz(t31[off]):.2f}); cosine(residual, (3,1) term) {np.sum(res[off]*t31[off])/(nz(res[off])*nz(t31[off])):+.3f}, cosine(residual, k3 term) {np.sum(res[off]*t3[off])/(nz(res[off])*nz(t3[off])):+.3f}")
