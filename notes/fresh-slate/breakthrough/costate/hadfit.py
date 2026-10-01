"""Team B's E2 test: fit the a-traceless old fluctuation D21°(t) = D_old - 1 (1^T D_old / n) (age > 2, inside FC) by
Hadamard products of a few transported symmetric matrices. Unconstrained LS over all pairs B_i o B_j lower-bounds the
residual of any sum_{k<=K} X_k o Y_k with X_k, Y_k in span(B); rank-K fits by alternating LS on G = alpha^T beta."""
import sys, os, json, numpy as np
from scipy.special import ndtr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../bench")); sys.path.insert(0, os.path.join(HERE, "../region"))
import bench, fc_dil, gclose
mlp = int(sys.argv[1]); layers = range(6, 14)
S_ = bench.load_set("w1024_d16"); W = bench.weights(S_, mlp).astype(np.float64)
cap = {}
fc_dil.run(bench.weights(S_, mlp), slices=2, k4mf=True, capture=cap, capA=2)
L, n, _ = W.shape
P = {}; Gam = {}
for t in range(L):
    mu, S = cap[t]["mu"], cap[t]["S"]; P[t] = ndtr(mu / np.sqrt(np.diag(S)))
Gam[0] = W[0].copy()
for t in range(1, L):
    Gam[t] = Gam[t - 1] @ (P[t - 1][:, None] * W[t])
out = []
for t in layers:
    mu, S, Do = cap[t]["mu"], cap[t]["S"], cap[t]["Dold"]
    Dc = Do - Do.mean(0, keepdims=True)                     # a-traceless fluctuation (columns centred over a)
    s2 = np.diag(S).copy(); one = np.ones(n)
    Ca = gclose.cov_exact(mu, S, P[t])
    G = Gam[t].T @ Gam[t]
    Sp = cap[t - 1]["S"]; Gp = Gam[t - 1].T @ Gam[t - 1]
    B = dict(Cz=S, Ca=Ca, GG=G, img_Cz=W[t].T @ Sp @ W[t], img_GG=W[t].T @ Gp @ W[t], diagC=np.diag(s2), I=np.eye(n),
             mm=np.outer(mu, mu), ss=np.outer(s2, s2), sm=np.outer(s2, mu), ms=np.outer(mu, s2),
             om=np.outer(one, mu), mo=np.outer(mu, one), os_=np.outer(one, s2), so=np.outer(s2, one), oo=np.outer(one, one))
    names = list(B); M = len(names)
    feats, fn = [], []
    for i in range(M):
        for j in range(i, M):
            F = B[names[i]] * B[names[j]]
            F = F - F.mean(0, keepdims=True)
            nf = np.linalg.norm(F)
            if nf < 1e-12 * max(1.0, np.linalg.norm(Dc)): continue
            feats.append((F / nf).astype(np.float32).ravel()); fn.append((i, j))
    X = np.stack(feats, 1)                                  # (n^2, p)
    y = Dc.ravel().astype(np.float32)
    coef, *_ = np.linalg.lstsq(X, y, rcond=1e-6)
    res = y - X @ coef
    frac_all = float(res @ res / (y @ y))
    # single-factor baselines: linear in the basis (no Hadamard) and the C6 template
    Xl = np.stack([(lambda F: (F - F.mean(0, keepdims=True)).ravel())(B[k]).astype(np.float32) for k in names], 1)
    cl, *_ = np.linalg.lstsq(Xl, y, rcond=1e-6); rl = y - Xl @ cl
    K = 2 * mu[:, None] * S + mu[None, :] * s2[:, None]; Kc = (K - K.mean(0, keepdims=True)).ravel()
    rk = y - (Kc @ y / (Kc @ Kc)) * Kc
    rec = dict(mlp=mlp, t=t, n_feat=len(fn), resid_all_pairs=frac_all, resid_linear=float(rl @ rl / (y @ y)),
               resid_C6=float(rk @ rk / (y @ y)), old_traceless_share=float(np.sum(Dc * Dc) / np.sum(Do * Do)))
    out.append(rec); print(json.dumps(rec), flush=True)
json.dump(out, open(os.path.join(HERE, "results_live", f"hadfit_mlp{mlp}.json"), "w"), indent=1)
