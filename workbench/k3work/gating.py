# Correlated gating of the bulk kappa3: exact first-order gate P(z_a>0,z_b>0,z_c>0) vs the factorized Phi_a Phi_b Phi_c.
# Relative correction (Mehler, first order in rho): sum over pairs (phi_a phi_b / Phi_a Phi_b) rho_ab.
# Split rho into the coherent top-k eigen part and the bulk remainder; check against exact trivariate orthant probabilities.
import pickle, numpy as np
from scipy.stats import norm, multivariate_normal
D = {d["layer"]: d for d in pickle.load(open("dump_oracleall_off0.pkl", "rb"))}
rng = np.random.default_rng(0)
for l in [3, 6, 9, 12, 14]:
    d = D[l]; s = np.sqrt(d["var"]); a = d["mu"] / s; ph, Ph = norm.pdf(a), norm.cdf(a)
    R = d["C_pre"] / np.outer(s, s); np.fill_diagonal(R, 1.0)
    ev, U = np.linalg.eigh(R); ev, U = ev[::-1], U[:, ::-1]
    q = ph / Ph
    trip = rng.choice(1024, size=(4000, 3))
    trip = trip[(trip[:, 0] != trip[:, 1]) & (trip[:, 1] != trip[:, 2]) & (trip[:, 0] != trip[:, 2])]
    def corr(Rm):
        i, j, k = trip.T
        return q[i]*q[j]*Rm[i, j] + q[i]*q[k]*Rm[i, k] + q[j]*q[k]*Rm[j, k]
    Roff = R - np.eye(1024)
    out = [np.sqrt(np.mean(corr(Roff)**2))]
    for k in [1, 4]:
        Rc = (U[:, :k] * ev[:k]) @ U[:, :k].T; np.fill_diagonal(Rc, 0)
        out.append(np.sqrt(np.mean(corr(Rc)**2)))
    # exact orthant check on 200 triples
    ex = []
    for (i, j, k) in trip[:200]:
        S3 = R[np.ix_([i, j, k], [i, j, k])]
        p = multivariate_normal(mean=np.zeros(3), cov=S3).cdf(a[[i, j, k]])
        ex.append(p / (Ph[i]*Ph[j]*Ph[k]) - 1)
    ex = np.array(ex)
    print(f"layer {l:2d}: rms rel. gate correction {out[0]:.3f} (top-1 part {out[1]:.3f}, top-4 {out[2]:.3f}); exact orthant rms {np.sqrt(np.mean(ex**2)):.3f}, mean {ex.mean():+.3f}; lambda1 {ev[0]:.0f}")
