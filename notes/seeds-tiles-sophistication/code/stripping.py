# The network strips input randomness: (a) activation-pattern sign disagreement between independent inputs per layer,
# (b) quenched 'randomness fraction' Var_theta F_j / E_theta F_j^2 of the final layer vs P(Z_L>0) = 1 - b_0^(L).
import numpy as np
spec = np.load("spec.npy")
rho = lambda c: c/2 + (np.sqrt(np.maximum(1-c*c, 0)) + c*np.arcsin(c))/np.pi
rng = np.random.default_rng(3)
n, L, M = 512, 16, 4000
Ws = [rng.standard_normal((n, n)) * np.sqrt(2.0/n) for _ in range(L)]
D = rng.standard_normal((M, n)); D /= np.linalg.norm(D, axis=1, keepdims=True)
H = D; c = (D[:M//2] * D[M//2:]).sum(1).mean()
print(" layer | sign disagreement (indep. inputs): observed  predicted arccos(rho^(l-1)(0))/pi | 'structure' bits (constant over all inputs)")
ct = 0.0
for l, W in enumerate(Ws, 1):
    Z = H @ W.T
    S = Z > 0
    dis = np.mean(S[:M//2] != S[M//2:])
    frozen = np.mean(np.all(S, axis=0) | np.all(~S, axis=0))
    if l in (1, 2, 4, 8, 12, 16):
        print(f"  {l:3d}  |   {dis:.4f}      {np.arccos(ct)/np.pi:.4f}      | {frozen:.3f}")
    ct = rho(ct)
    H = np.maximum(Z, 0)
F = H                                    # final post-ReLU outputs on unit directions
frac = np.array([F.var(0).sum() / (F**2).mean(0).sum()])
print(f"final layer: sum_j Var_theta F_j / sum_j E F_j^2 = {frac.mean():.4f}   (annealed infinite-width prediction P(Z_L>0) = {1-spec[L-1][0]:.4f})")
