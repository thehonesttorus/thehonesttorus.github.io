"""Proposition 1 check: unbiased sketches of D = sum_r w_r [(y_r o y_r) z_r^T + 2 (y_r o z_r) y_r^T] with
scalar cube-root signs (averaged over d^2 samples = equal cost) vs Haar U(d) signs (one sample, X = (1/d) Re tr).
Synthetic Gaussian atoms (nearly orthogonal, as the memory, F8.3)."""
import numpy as np
from scipy.stats import unitary_group
rng = np.random.default_rng(0)
n, T, reps = 48, 192, 200
Y = rng.standard_normal((T, n)) / np.sqrt(n); Z = rng.standard_normal((T, n)) / np.sqrt(n); w = np.ones(T)
D = ((Y * Y) * w[:, None]).T @ Z + 2 * (((Y * Z) * w[:, None]).T @ Y)
nD = np.sum(D ** 2)
om = np.exp(2j * np.pi / 3)
def scalar():
    xi = om ** rng.integers(0, 3, T)
    u = Y.T @ xi; v = Z.T @ (w * np.conj(xi) ** 2)
    return np.real(np.outer(u * u, v) + 2 * np.outer(u * v, u))
def mat(d):
    U = unitary_group.rvs(d, size=T, random_state=rng) if d > 1 else None
    Ui = np.linalg.inv(U); U2i = Ui @ Ui
    u = np.einsum('rij,ra->aij', U, Y); v = np.einsum('rij,r,rb->bij', U2i, w, Z)
    # term1: (1/d) tr(u_a u_a v_b); term2: (1/d) tr(u_a v_a u_b)
    uu = np.einsum('aij,ajk->aik', u, u); uv = np.einsum('aij,ajk->aik', u, v)
    t1 = np.einsum('aik,bki->ab', uu, v); t2 = np.einsum('aik,bki->ab', uv, u)
    return np.real(t1 + 2 * t2) / d
for d in (2, 3, 4):
    es = [np.sum((np.mean([scalar() for _ in range(d * d)], 0) - D) ** 2) / nD for _ in range(reps)]
    em = [np.sum((mat(d) - D) ** 2) / nD for _ in range(reps)]
    bias = np.sum((np.mean([mat(d) for _ in range(400)], 0) - D) ** 2) / nD
    print(f"d={d}: scalar x{d*d} relvar {np.mean(es):.1f}  |  U({d}) one sample relvar {np.mean(em):.1f}  ratio {np.mean(em)/np.mean(es):.2f}  (bias check, mean of 400 U-samples: {bias:.2f} vs expected var/400 = {np.mean(em)/400:.2f})")
