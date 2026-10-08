# (E3) Mixed free cumulants of the gate projections: exact trivariate orthant P(a,b,c) vs the product gate, split into the
# first-order (odd in R: incoherent) and second-order (even: has an annealed part) Mehler terms, over random triples.
import pickle, numpy as np
from scipy.stats import norm, multivariate_normal
D = {d["layer"]: d for d in pickle.load(open("../k3work/dump_oracleall_off0.pkl", "rb"))}
rng = np.random.default_rng(0)
print("layer | rms(P-F)/rms F | 1st order rms | 2nd order (Mehler) rms, mean | exact residual after 1st: rms, mean | after 2nd: rms, mean   (all relative to rms F)")
for l in [3, 6, 9, 12, 14]:
    d = D[l]; s = np.sqrt(d["var"]); a = d["mu"] / s; ph, Ph = norm.pdf(a), norm.cdf(a)
    R = d["C_pre"] / np.outer(s, s); np.fill_diagonal(R, 1.0)
    trip = rng.choice(1024, size=(1500, 3)); trip = trip[(trip[:, 0] != trip[:, 1]) & (trip[:, 1] != trip[:, 2]) & (trip[:, 0] != trip[:, 2])][:600]
    P, F, M1, M2 = [], [], [], []
    for (i, j, k) in trip:
        idx = [i, j, k]; S3 = R[np.ix_(idx, idx)]
        P.append(multivariate_normal(mean=np.zeros(3), cov=S3).cdf(a[idx])); F.append(Ph[i] * Ph[j] * Ph[k])
        pairs = [(i, j, k), (i, k, j), (j, k, i)]
        M1.append(sum(ph[x] * ph[y] * R[x, y] * Ph[z] for x, y, z in pairs))
        # second order: same-pair rho^2 terms (Hermite He_1 = a) and the two-different-pairs product terms
        m2 = sum(0.5 * R[x, y]**2 * (a[x] * ph[x]) * (a[y] * ph[y]) * Ph[z] for x, y, z in pairs)
        m2 += ph[i] * ph[j] * ph[k] * (R[i, j] * R[i, k] + R[i, j] * R[j, k] + R[i, k] * R[j, k])
        M2.append(m2)
    P, F, M1, M2 = map(np.array, (P, F, M1, M2)); rms = lambda x: np.sqrt(np.mean(x**2)); nF = rms(F)
    print(f"  {l:2d} | {rms(P-F)/nF:.4f} | {rms(M1)/nF:.4f} | {rms(M2)/nF:.5f}, {M2.mean()/nF:+.5f} | {rms(P-F-M1)/nF:.5f}, {(P-F-M1).mean()/nF:+.5f} | {rms(P-F-M1-M2)/nF:.5f}, {(P-F-M1-M2).mean()/nF:+.5f}", flush=True)
