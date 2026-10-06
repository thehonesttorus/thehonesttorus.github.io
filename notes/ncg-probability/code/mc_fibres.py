# Monte Carlo tests of the noncommutative-probability reading, official network NET (1e6 inputs):
#  (E6) scale ladder: free-limit constant mean vs quenched spread of true means vs closure / GAC / chain errors
#  (E2) the gain is the trace fluctuation: Var(|h_l|^2)/E^2 minus the Gaussian part, per layer, vs GAC's gamma_l and L/n
#  (E5) law of G^2: empirical vs Gamma vs lognormal; E[G] from each
#  (E4) quasi-free fibres: excess kurtosis / skewness of linear functionals of z_L inside radial bins vs inside directional bins
import numpy as np, sys, time
sys.path.insert(0, "../num12")
from gac import gac, EG
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
T = int(float(sys.argv[2])) if len(sys.argv) > 2 else 1_000_000
Wc = np.load(f"../official/W_off{net}.npy"); mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
n, L = Wc.shape[1], Wc.shape[0]
Wf = [np.ascontiguousarray(W.T).astype(np.float32) for W in Wc]   # row convention: h @ Wf
rng = np.random.default_rng(12345)
# directions for E4: top covariance eigenvectors of the chain's last-layer pre-activation (from GAC state), plus random
g = gac(list(Wc.astype(np.float64)))
S_last = g[-1]["S"]; mu_last = g[-1]["mu"]; gam = np.array([d["gamma"] for d in g])
ev, U = np.linalg.eigh(S_last); Utop = U[:, -4:][:, ::-1]
Q, _ = np.linalg.qr(rng.standard_normal((n, 32))); dirs = np.concatenate([Utop, mu_last[:, None] / np.linalg.norm(mu_last), Q], axis=1).astype(np.float32)  # 4 top + mean + 32 random
nd = dirs.shape[1]
B = 8192; done = 0; t0 = time.time()
norm_s1 = np.zeros(L); norm_s2 = np.zeros(L); norm_s3 = np.zeros(L); norm_s4 = np.zeros(L)
proj = np.empty((T, nd), dtype=np.float32); radius = np.empty(T, dtype=np.float32); rad_prev = np.empty(T, dtype=np.float32)
while done < T:
    b = min(B, T - done)
    h = rng.standard_normal((b, n), dtype=np.float32)
    for l in range(L):
        z = h @ Wf[l]
        if l == L - 1:
            proj[done:done + b] = z @ dirs; radius[done:done + b] = np.einsum("ij,ij->i", z, z); rad_prev[done:done + b] = np.einsum("ij,ij->i", h, h)
        h = np.maximum(z, 0)
        q = np.einsum("ij,ij->i", h, h).astype(np.float64)
        norm_s1[l] += q.sum(); norm_s2[l] += (q**2).sum(); norm_s3[l] += (q**3).sum(); norm_s4[l] += (q**4).sum()
    done += b
print(f"net {net}: {T} samples in {time.time()-t0:.0f}s", flush=True)
# ---- E6: scale ladder
m = mt[-1]; print(f"[E6] last layer: mean of true means {m.mean():.4f}, spread across neurons (quenched 1/sqrt(n) fluctuation) rms {m.std():.4f}; free-limit constant-vector MSE {np.mean((m - m.mean())**2):.3e}; closure 4.06e-6, GAC 1.56e-6, K3 chain 2.19e-8, leaders ~1.2e-8 (raw)")
# ---- E2/E5: norm fluctuation per layer vs gamma_l
print("[E2] layer | Var(|h|^2)/E^2 (total) | Gaussian part 2tr(S^2)/tr(S)^2 (closure) | excess = gain | GAC gamma_l | l/n")
for l in range(L):
    E1, E2 = norm_s1[l] / T, norm_s2[l] / T; tot = E2 / E1**2 - 1
    # Gaussian part from the post covariance of the closure chain state: use GAC's conditional (M, C) at layer l
    Mv = g[l]["M"]; Cv = None
    line = f"  {l:2d} | {tot:.5f} |"
    print(line + f" (gamma_l {gam[l]:.5f}, l/n {(l+1)/n:.5f})")
# ---- E5: law of G^2 at the last layer (post-activation norm)
q = radius.astype(np.float64)  # pre-activation |z_L|^2
qn = q / q.mean(); Vq = qn.var(); skew = ((qn - 1)**3).mean() / Vq**1.5; kurt = ((qn - 1)**4).mean() / Vq**2 - 3
A = 1 / Vq
print(f"[E5] |z_L|^2/E: Var {Vq:.5f} (gamma_L from GAC {gam[-1]:.5f}); skew {skew:.3f} (Gamma law 2/sqrt(A)={2/np.sqrt(A):.3f}, lognormal {(np.exp(np.log1p(Vq))+2)*np.sqrt(np.expm1(np.log1p(Vq))):.3f}); excess kurt {kurt:.3f} (Gamma 6/A={6/A:.3f})")
print(f"     E[sqrt(q/E q)] empirical {np.sqrt(qn).mean():.6f}; Gamma law {EG(Vq, 'gamma'):.6f}; lognormal {EG(Vq, 'lognormal'):.6f}; first order 1-V/8 {1 - Vq/8:.6f}")
# ---- E4: quasi-free fibres. Standardise projections; excess kurtosis and skewness of each projection overall, within radial bins, within directional bins
P = proj.astype(np.float64); P = (P - P.mean(0)) / P.std(0)
def ex_kurt(x): x = x - x.mean(); v = x.var(); return (x**4).mean() / v**2 - 3
def skw(x): x = x - x.mean(); v = x.var(); return (x**3).mean() / v**1.5
nb = 16
rbins = np.quantile(q, np.linspace(0, 1, nb + 1)); rid = np.clip(np.searchsorted(rbins, q, side="right") - 1, 0, nb - 1)
def within(stat, cond_id):
    out = []
    for j in range(nd):
        vals = []; w = []
        for bidx in range(nb):
            sel = cond_id == bidx
            if sel.sum() > 1000: vals.append(stat(P[sel, j])); w.append(sel.sum())
        out.append(np.average(vals, weights=w))
    return np.array(out)
labels = ["eig1", "eig2", "eig3", "eig4", "mean"] + ["rand"] * 32
K_all = np.array([ex_kurt(P[:, j]) for j in range(nd)]); S_all = np.array([skw(P[:, j]) for j in range(nd)])
K_rad = within(ex_kurt, rid); S_rad = within(skw, rid)
# directional bins: condition on the top eigen-projection (eig1) and on the mean projection
for cname, cj in [("eig1", 0), ("mean", 4)]:
    cb = np.quantile(P[:, cj], np.linspace(0, 1, nb + 1)); cid = np.clip(np.searchsorted(cb, P[:, cj], side="right") - 1, 0, nb - 1)
    K_dir = within(ex_kurt, cid); S_dir = within(skw, cid)
    print(f"[E4] conditioning on {cname}-projection bins vs radial bins (|z_L|^2):")
    print("      direction   excess-kurt: all / radial-bins / dir-bins   |  skew: all / radial-bins / dir-bins")
    for j in list(range(5)) + [5]:
        tag = labels[j] if j < 5 else "rand(avg of 32)"
        if j < 5:
            print(f"      {tag:9s}   {K_all[j]:+.4f} / {K_rad[j]:+.4f} / {K_dir[j]:+.4f}   |  {S_all[j]:+.4f} / {S_rad[j]:+.4f} / {S_dir[j]:+.4f}")
        else:
            print(f"      {tag:9s}   {K_all[5:].mean():+.4f} / {K_rad[5:].mean():+.4f} / {K_dir[5:].mean():+.4f}   |  {S_all[5:].mean():+.4f} / {S_rad[5:].mean():+.4f} / {S_dir[5:].mean():+.4f}")
# independence of radius and direction: corr(|z|^2, proj^2) for each direction, scale-mixture prediction = Var(G^2)-driven positive value equal for all directions
c = np.array([np.corrcoef(q, P[:, j]**2)[0, 1] for j in range(nd)])
print(f"[E4] corr(|z_L|^2, proj^2): eig1 {c[0]:+.4f} eig2 {c[1]:+.4f} eig3 {c[2]:+.4f} mean {c[4]:+.4f} random mean {c[5:].mean():+.4f} +- {c[5:].std():.4f}  (scale mixture: all equal; shift mixture: direction-specific)")
np.savez(f"mc_fibres_off{net}.npz", norm_s1=norm_s1, norm_s2=norm_s2, gam=gam, K_all=K_all, K_rad=K_rad, S_all=S_all, S_rad=S_rad, c=c, Vq=Vq)
