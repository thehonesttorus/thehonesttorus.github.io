# E0 (note XLII, frame N1): how much of each network's final-layer error do per-network counterterm amplitudes reach?
#   python e0_pernet.py RESPDIR(with calresp/, lamresp/, caldirs.json) OUTDIR(ofin56_{n}.npy) TRUTHDIR [LOCDIR]
# R_n = the 103 stored linear responses (1024 x 103, last layer), e_n = the adopted system's final-layer error.
#   pi_R      = |P_R e|^2 / |e|^2 (chance level 103/1024 = 0.10)
#   honest per-network gain: amplitudes fitted on a random half of the neurons, MSE measured on the other half
#     (10 random splits; ridge chosen on the fitting half), against the common amplitudes already in the system (b = 0)
#     and against chance (the same fit on a random vector with e's norm)
#   the between-network spread of the per-network optimal amplitudes, and how well they are predicted by the common fit
# With LOCDIR (chaindump/mccache/W/truth for some networks): pi_U(r) for the top-r eigenvectors of the gated MC
# pre-activation covariance of layer 15 (D Cov(h_15) D, D = diag Phi(alpha)), r = 16/32/64/128/256 (chance r/1024).
import json, os, sys, numpy as np
from scipy.special import ndtr
rd, od, td = sys.argv[1], sys.argv[2], sys.argv[3]; ld = sys.argv[4] if len(sys.argv) > 4 else None
dirs = json.load(open(f"{rd}/caldirs.json")); K = len(dirs)
rng = np.random.default_rng(5)
rows = []
for n in range(100):
    fb = f"{rd}/lamresp/out_b_{n}.npy"
    if not os.path.exists(fb) or not all(os.path.exists(f"{rd}/calresp/out_c{j}_{n}.npy") for j in range(K)):
        continue
    b = np.load(fb)[-1]
    R = np.stack([(np.load(f"{rd}/calresp/out_c{j}_{n}.npy")[-1] - b) / dirs[j][2] for j in range(K)], 1)
    e = np.load(f"{od}/ofin56_{n}.npy")[-1] - np.load(f"{td}/truth_off{n}.npz")["m"].astype(np.float64)[-1]
    Q, _ = np.linalg.qr(R)
    piR = np.sum((Q.T @ e) ** 2) / np.sum(e ** 2)
    gains = []; gch = []
    for s in range(10):
        p = rng.permutation(1024); a, t = p[:512], p[512:]
        for vec, store in ((e, gains), (rng.standard_normal(1024) * np.linalg.norm(e) / 32.0, gch)):
            best = None
            for lam in (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0):
                G = R[a].T @ R[a]; pen = lam * np.trace(G) / K
                bb = -np.linalg.solve(G + pen * np.eye(K), R[a].T @ vec[a])
                # inner validation inside the fitting half to choose the ridge
                a1, a2 = a[:256], a[256:]
                G1 = R[a1].T @ R[a1]; b1 = -np.linalg.solve(G1 + lam * np.trace(G1) / K * np.eye(K), R[a1].T @ vec[a1])
                v = np.mean((vec[a2] + R[a2] @ b1) ** 2)
                if best is None or v < best[0]:
                    best = (v, bb)
            store.append(np.mean((vec[t] + R[t] @ best[1]) ** 2) / np.mean(vec[t] ** 2))
    bopt = -np.linalg.lstsq(R, e, rcond=None)[0]
    rows.append(dict(n=n, mse=np.mean(e ** 2), piR=piR, cv=np.mean(gains), cv_chance=np.mean(gch), bopt=bopt))
    if len(rows) % 20 == 0:
        print(f"{len(rows)} networks", flush=True)
print(f"\nE0 over {len(rows)} networks (adopted system with its common counterterms; final layer)")
piR = np.array([r["piR"] for r in rows]); cv = np.array([r["cv"] for r in rows]); ch = np.array([r["cv_chance"] for r in rows])
print(f"pi_R (share of the error in the 103-response span): mean {piR.mean():.3f}, median {np.median(piR):.3f}, "
      f"10-90% [{np.quantile(piR, .1):.3f}, {np.quantile(piR, .9):.3f}]   (chance 0.101)")
print(f"neuron-split cross-validated per-network refit: held-out MSE ratio mean {cv.mean():.3f}, median {np.median(cv):.3f} "
      f"(chance-vector control {ch.mean():.3f}); MSE-weighted ratio "
      f"{np.sum(cv * [r['mse'] for r in rows]) / np.sum([r['mse'] for r in rows]):.3f}")
B = np.array([r["bopt"] for r in rows])
print(f"per-network optimal amplitudes: |mean over networks| / std over networks, median over the 103 directions: "
      f"{np.median(np.abs(B.mean(0)) / B.std(0)):.3f}  (small = per-network spread dominates any common shift)")
if ld:
    out = []
    for n in range(100):
        f = f"{ld}/mccache_{n}.npz"
        if not os.path.exists(f) or not os.path.exists(f"{ld}/chaindump_{n}.npz"):
            continue
        mc = np.load(f); cd = np.load(f"{ld}/chaindump_{n}.npz")
        W = np.load(f"{ld}/W_off{n}.npy").astype(np.float64)
        S = float(mc["S"]); m1 = mc["s1"][15] / S
        Cv = mc["Hc"][15].astype(np.float64) - np.outer(m1, m1)
        mu = W[15] @ cd["out"][14]; Dg = ndtr(mu / np.sqrt(cd["var"][15]))
        ev, V = np.linalg.eigh(Dg[:, None] * Cv * Dg[None, :]); V = V[:, ::-1]
        e = np.load(f"{od}/ofin56_{n}.npy")[-1] - np.load(f"{td}/truth_off{n}.npz")["m"].astype(np.float64)[-1]
        pr = (V.T @ e) ** 2
        out.append([pr[:r].sum() / pr.sum() for r in (16, 32, 64, 128, 256)])
    if out:
        print(f"pi_U(r) of the final error in the top-r modes of the gated layer-15 covariance ({len(out)} networks), "
              f"r = 16/32/64/128/256: " + " ".join(f"{x:.3f}" for x in np.mean(out, 0)) + "   (chance 0.016/0.031/0.062/0.125/0.25)")
