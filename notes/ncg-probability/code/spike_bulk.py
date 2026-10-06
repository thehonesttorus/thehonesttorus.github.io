# Test of the "spike + annealed bulk" covariance: Sigma_r = top-r eigenpart + (trace of the rest/(n-r)) (I - P_r).
# The next layer's per-neuron variance var_i = w_i' Sigma w_i is what the gate reads.  Quenched full Sigma vs Sigma_r:
# relative rms error of var_i and the implied final-mean MSE (delta m = phi(a) sigma/2 * delta var/var), official net 0.
import pickle, numpy as np
from scipy.stats import norm
D = {d["layer"]: d for d in pickle.load(open("../k3work/dump_oracleall_off0.pkl", "rb"))}
Wc = np.load("../official/W_off0.npy").astype(np.float64)
print("layer | PR of Sigma | bulk share of trace | r: rel rms err of var_i -> implied mean-MSE   (frontier chain final MSE 2.2e-8)")
for l in [4, 8, 12, 14]:
    Sig = D[l]["C_pre"].copy(); np.fill_diagonal(Sig, D[l]["var"]); W = Wc[l + 1]
    ev, U = np.linalg.eigh(Sig); ev, U = ev[::-1], U[:, ::-1]
    PR = ev.sum()**2 / (ev**2).sum()
    var_full = np.einsum("ij,ij->i", W @ Sig, W)
    mu_next = W @ D[l]["mu"]  # proxy for the next-layer mean (chain pre-activation mean transported); only the scale matters
    out = []
    for r in [0, 4, 16, 64, 128, 256, 512, 768]:
        P = U[:, :r]; spike = (P * ev[:r]) @ P.T
        sb = ev[r:].sum() / (1024 - r)
        Sr = spike + sb * (np.eye(1024) - P @ P.T)
        var_r = np.einsum("ij,ij->i", W @ Sr, W)
        rel = (var_r - var_full) / var_full
        s = np.sqrt(var_full); a = mu_next / s
        dm = norm.pdf(a) * s / 2 * rel
        out.append(f"{r}: {np.sqrt(np.mean(rel**2)):.4f} -> {np.mean(dm**2):.1e}")
    print(f"  {l:2d} | {PR:6.1f} | {ev[64:].sum()/ev.sum():.2f} (beyond r=64) | " + "  ".join(out))
