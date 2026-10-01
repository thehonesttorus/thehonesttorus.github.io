"""T4: the degree-0 (face-wise constant) route to second moments.

For degree-1 homogeneous f, g (Euler + Stein):  E[f g] = E[grad f . grad g] + E[g Lap f].
With J_l = d z_l / d x (constant on every face), G_l = J_l^T J_l obeys exactly G_{l+1} = W^T D_l G_l D_l W, D_l = diag(g_l(x)).
Face decoupling (gates independent of the gradient Gram):  E G_{l+1} ~= W^T (E[g g^T] o E G_l) W.
Measured: (a) the split E[z z^T] = E G + E[z Lap z^T] (sizes); (b) one-step relative error of the decoupled E G_{l+1};
(c) for comparison, one-step relative error of the Gaussian (Mehler) closure for E[z z^T] at l+1 given the true (mu, C) at l.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bench"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench import weights_from_seed
import fbt

def run(n, N, seed=11, L=16, chunk=2000):
    W = weights_from_seed(seed, n, L).astype(np.float64)
    rng = np.random.default_rng(17 + seed)
    acc = [dict(G=0, gGg=0, gg=0, zz=0, z=0, a=0, aa=0) for _ in range(L)]
    done = 0
    while done < N:
        x = rng.standard_normal((chunk, n)); h = x; J = np.broadcast_to(np.eye(n), (chunk, n, n))
        for l in range(L):
            z = h @ W[l]; J = J @ W[l]                    # J[b] = dz_l/dx, shape (n_in, n)
            g = (z > 0).astype(float); a = z * g
            G = np.einsum("bij,bik->bjk", J, J, optimize=True)
            A = acc[l]
            A["G"] = A["G"] + G.sum(0); A["gGg"] = A["gGg"] + np.einsum("bj,bjk,bk->jk", g, G, g, optimize=True)
            A["gg"] = A["gg"] + g.T @ g; A["zz"] = A["zz"] + z.T @ z; A["z"] = A["z"] + z.sum(0)
            A["a"] = A["a"] + a.sum(0); A["aa"] = A["aa"] + a.T @ a
            h = a; J = J * g[:, None, :]
        done += chunk
    acc = [{k: v / N for k, v in A.items()} for A in acc]
    rel = lambda a, b: float(np.linalg.norm(a - b) / np.linalg.norm(b))
    reld = lambda a, b: float(np.linalg.norm(np.diag(a) - np.diag(b)) / np.linalg.norm(np.diag(b)))
    print(f"# n={n} N={N} seed={seed}")
    print("z-layer | |E G|/|E zz^T| | |E zLapz|/|E zz^T| | decoupled E G_{l+1}: rel err (all / diag) | Gaussian-Mehler E zz^T_{l+1}: rel err (all / diag) | exact gGg transport check")
    for l in range(L - 1):
        A, B = acc[l], acc[l + 1]
        EG, EZZ = A["G"], A["zz"]
        lap = EZZ - EG
        dec = W[l + 1].T @ (A["gg"] * EG) @ W[l + 1]
        exact = W[l + 1].T @ A["gGg"] @ W[l + 1]
        mu = A["z"]; C = EZZ - np.outer(mu, mu)
        Ca = fbt.mehler_cov(mu, C, 6)
        var = np.diag(C); s = np.sqrt(var); t = mu / s
        Ea = fbt.readout(mu, var, 0 * mu, 0 * mu)
        Ea2 = (mu * mu + var) * fbt.ndtr(t) + mu * s * fbt._phi(t)
        np.fill_diagonal(Ca, Ea2 - Ea * Ea)
        Saa = Ca + np.outer(Ea, Ea)
        gz = W[l + 1].T @ Saa @ W[l + 1]
        print(f"{l+1:7d} | {np.linalg.norm(EG)/np.linalg.norm(EZZ):.3f} | {np.linalg.norm(lap)/np.linalg.norm(EZZ):.3f} | "
              f"{rel(dec, B['G']):.4f} / {reld(dec, B['G']):.4f} | {rel(gz, B['zz']):.4f} / {reld(gz, B['zz']):.4f} | {rel(exact, B['G']):.1e}", flush=True)

if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]))
