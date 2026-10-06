# Absolute correlated-gating error of the bulk kappa3 transport: exact trivariate orthant P vs Phi_a Phi_b Phi_c,
# as an rms ratio over random distinct triples; plus the first-order Mehler pair term and its coherent top-1 share.
import pickle, numpy as np
from scipy.stats import norm, multivariate_normal
D = {d["layer"]: d for d in pickle.load(open("dump_oracleall_off0.pkl", "rb"))}
rng = np.random.default_rng(0)
for l in [3, 6, 9, 12, 14]:
    d = D[l]; s = np.sqrt(d["var"]); a = d["mu"] / s; ph, Ph = norm.pdf(a), norm.cdf(a)
    R = d["C_pre"] / np.outer(s, s); np.fill_diagonal(R, 1.0)
    ev, U = np.linalg.eigh(R); u1 = U[:, -1]; l1 = ev[-1]
    trip = rng.choice(1024, size=(400, 3))
    trip = trip[(trip[:, 0] != trip[:, 1]) & (trip[:, 1] != trip[:, 2]) & (trip[:, 0] != trip[:, 2])]
    P, F, M1, Mc = [], [], [], []
    for (i, j, k) in trip:
        idx = [i, j, k]; S3 = R[np.ix_(idx, idx)]
        P.append(multivariate_normal(mean=np.zeros(3), cov=S3).cdf(a[idx]))
        F.append(Ph[i]*Ph[j]*Ph[k])
        M1.append(ph[i]*ph[j]*R[i, j]*Ph[k] + ph[i]*ph[k]*R[i, k]*Ph[j] + ph[j]*ph[k]*R[j, k]*Ph[i])
        rc = lambda x, y: l1*u1[x]*u1[y]
        Mc.append(ph[i]*ph[j]*rc(i, j)*Ph[k] + ph[i]*ph[k]*rc(i, k)*Ph[j] + ph[j]*ph[k]*rc(j, k)*Ph[i])
    P, F, M1, Mc = map(np.array, (P, F, M1, Mc))
    rms = lambda x: np.sqrt(np.mean(x**2))
    print(f"layer {l:2d}: rms(P - PhiPhiPhi)/rms(PhiPhiPhi) {rms(P-F)/rms(F):.3f} | first-order pair term {rms(M1)/rms(F):.3f}, "
          f"its top-1 coherent part {rms(Mc)/rms(F):.3f} | resid after first order {rms(P-F-M1)/rms(F):.3f} | mean(P-F)/rms F {np.mean(P-F)/rms(F):+.3f}")
