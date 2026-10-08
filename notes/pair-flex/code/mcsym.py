# One Monte Carlo chunk for test S of note XXXVIII (notes/pair-flex): the collective symbols of the post-activation
# fourth cumulant.
#   python mcsym.py NET NSAMP CHUNK K LAYERS OUTFILE [MCPREFIX]      e.g. python mcsym.py 1 1e6 0 64 2,5,8,11,14 sym1_0.npz mc4_off1
# At each layer l in LAYERS: y = relu(z_l), e = y - mu_y (true mean from MCPREFIX_full.npz), U = top-K eigenvectors of the
# true Cov(y) (same file), t = U^T e, P = the pair products t_k t_m (k <= m), S_km = sum_c U_ck U_cm e_c^2. Accumulates the
# raw sums from which symtest.py forms, exactly, N_a^off = U^T Psi_a U (unit a's (2,1,1) block projected on the collective
# space, every index coincidence removed) and T4^off (the all-distinct class projected):
#   s1..s4 = sum e^p (n);  C = sum e e^T, D = sum e^2 e^T (n x n);  et, e2t, e3t = sum e^p t^T (n x K);
#   etp, e2tp = sum e^p P^T (n x q);  e2S = sum e^2 S^T (n x q);
#   t1 = sum t (K);  tp = sum P (q);  tpt = sum P t^T (q x K);  tptp = sum P P^T (q x q);  S1 = sum S (q);  SS = sum S S^T.
# Products in float32 per batch, sums in float64. Seed default_rng([4713, NET, CHUNK]).
import sys, time, numpy as np
net, N, chunk, K = int(sys.argv[1]), int(float(sys.argv[2])), int(sys.argv[3]), int(sys.argv[4])
LAY = [int(x) for x in sys.argv[5].split(",")]; outf = sys.argv[6]
pre = sys.argv[7] if len(sys.argv) > 7 else f"mc4_off{net}"
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float32)
L, n, _ = Wcol.shape
F = np.load(f"{pre}_full.npz")
iu = np.triu_indices(K); q = len(iu[0])
U, MY, UU = {}, {}, {}
for l in LAY:
    Cy = np.asarray(F["cov_y"][l], np.float64); Cy = 0.5 * (Cy + Cy.T)
    ev, V = np.linalg.eigh(Cy)
    U[l] = np.ascontiguousarray(V[:, ::-1][:, :K], np.float32)
    MY[l] = np.asarray(F["mu_y"][l], np.float32)
    UU[l] = np.ascontiguousarray((U[l][:, iu[0]] * U[l][:, iu[1]]).T, np.float32)      # q x n
del F
keys = dict(s1=(n,), s2=(n,), s3=(n,), s4=(n,), C=(n, n), D=(n, n), et=(n, K), e2t=(n, K), e3t=(n, K), etp=(n, q), e2tp=(n, q),
            e2S=(n, q), t1=(K,), tp=(q,), tpt=(q, K), tptp=(q, q), S1=(q,), SS=(q, q))
acc = {l: {k: np.zeros(s) for k, s in keys.items()} for l in LAY}
rng = np.random.default_rng([4713, net, chunk])
B = 8192; done = 0; t0 = time.time(); lmax = max(LAY)
while done < N:
    b = min(B, N - done)
    Y = rng.standard_normal((n, b), dtype=np.float32)
    for l in range(lmax + 1):
        Z = Wcol[l] @ Y
        Y = np.maximum(Z, 0.0, out=Z)
        if l in acc:
            a = acc[l]
            e = Y - MY[l][:, None]; E2 = e * e; E3 = E2 * e
            t = U[l].T @ e; P = t[iu[0]] * t[iu[1]]; S = UU[l] @ E2
            a["s1"] += e.sum(1, dtype=np.float64); a["s2"] += E2.sum(1, dtype=np.float64)
            a["s3"] += E3.sum(1, dtype=np.float64); a["s4"] += (E2 * E2).sum(1, dtype=np.float64)
            a["C"] += e @ e.T; a["D"] += E2 @ e.T
            a["et"] += e @ t.T; a["e2t"] += E2 @ t.T; a["e3t"] += E3 @ t.T
            a["etp"] += e @ P.T; a["e2tp"] += E2 @ P.T; a["e2S"] += E2 @ S.T
            a["t1"] += t.sum(1, dtype=np.float64); a["tp"] += P.sum(1, dtype=np.float64); a["tpt"] += P @ t.T
            a["tptp"] += P @ P.T; a["S1"] += S.sum(1, dtype=np.float64); a["SS"] += S @ S.T
    done += b
    if (done // B) % 16 == 0:
        print(f"{done} samples, {time.time() - t0:.0f}s", flush=True)
out = {"c": done, "K": K, "layers": np.array(LAY)}
for l in LAY:
    out[f"U_{l}"] = U[l]; out[f"my_{l}"] = MY[l]
    for k in keys:
        out[f"{k}_{l}"] = acc[l][k]
np.savez(outf, **out)
print(f"done {done} samples in {time.time() - t0:.0f}s", flush=True)
