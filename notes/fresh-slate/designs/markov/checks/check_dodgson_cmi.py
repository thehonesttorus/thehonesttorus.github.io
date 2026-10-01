"""P2: Dodgson condensation and CMI = -1/2 log(1 - r^2) for a Gaussian on a ⊔ B ⊔ c (a, c singletons)."""
import numpy as np
rng = np.random.default_rng(0)
for trial in range(5):
    k = 5; A = rng.standard_normal((k, k + 3)); M = A @ A.T
    a, B, c = [0], [1, 2, 3], [4]
    d = lambda r, s=None: np.linalg.det(M[np.ix_(r, r if s is None else s)])
    aB, Bc = a + B, B + c
    lhs = d(a + B + c) * d(B); rhs = d(aB) * d(Bc) - d(aB, Bc) ** 2
    r2 = d(aB, Bc) ** 2 / (d(aB) * d(Bc))
    cmi_dodgson = -0.5 * np.log(1 - r2)
    cmi_direct = 0.5 * np.log(d(aB) * d(Bc) / (d(B) * d(a + B + c)))
    print(f"exchange rel. resid {abs(lhs-rhs)/abs(lhs):.1e}  CMI {cmi_dodgson:.6f} vs {cmi_direct:.6f}")
# Markov case: precision zero on (a,c) => r = 0
P = np.linalg.inv(M); P[0, 4] = P[4, 0] = 0; P += np.eye(5) * (1e-9 - min(0, np.linalg.eigvalsh(P).min()) * 1.1)
M = np.linalg.inv(P); d = lambda r, s=None: np.linalg.det(M[np.ix_(r, r if s is None else s)])
print("Markov: cross minor", d([0, 1, 2, 3], [1, 2, 3, 4]))
