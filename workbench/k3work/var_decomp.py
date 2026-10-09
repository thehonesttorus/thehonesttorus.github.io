# Where does the chain's late pre-activation variance error come from? (note XLIII section 3)
#   python var_decomp.py LOCDIR MC2DIR NET [NET ...]
# For every source layer s and target t = s + 1, with KG(mu, C) the pure Gaussian-closure post-activation covariance
# (bivariate normal ReLU covariance; Hermite series in the correlation, 40 terms (max |rho| 0.6), exact diagonal), the identity
#   dv_i = var_chain_t,i - var_true_t,i = p_i + x_i - g_i
# splits the per-neuron variance error into
#   p (propagated):        diag W_t [KG(chain state_s) - KG(true state_s)] W_t^T
#                          further split into the parts from the chain's means, its variances and its off-diagonal C
#   x (chain correction):  var_chain_t - diag W_t KG(chain state_s) W_t^T  (everything the chain adds beyond the
#                          Gaussian closure: kappa3/kappa4 terms, mixture gains, counterterms)
#   g (true defect):       var_true_t - diag W_t KG(true state_s) W_t^T  (what the exact law adds beyond the Gaussian
#                          closure of the true covariance)
# True states come from the 16M-input Monte Carlo (full and two independent halves A, B); the noise-free energy of a
# component uses the product of its A and B versions. LOCDIR: chaindump_{n}.npz, W_off{n}.npy. MC2DIR: mc2_off{n}_*.npz.
import sys, numpy as np
from scipy.special import ndtr

loc, mcd = sys.argv[1], sys.argv[2]; nets = [int(v) for v in sys.argv[3:]]
M = 40
SQ2PI = np.sqrt(2.0 * np.pi)


def gcoef(al):
    """g[m] = c_m(alpha) / sqrt(m!) for m = 1..M, where relu(alpha + z) = sum_m c_m He_m(z) / m!."""
    phi = np.exp(-0.5 * al * al) / SQ2PI
    g = np.zeros((M + 1, al.size))
    g[1] = ndtr(al)
    h_prev, h = np.zeros_like(al), np.ones_like(al)          # h_k = He_k(alpha) / sqrt(k!), k = -1, 0
    for m in range(2, M + 1):
        k = m - 2                                           # need h_{m-2}
        if k >= 1:
            h_prev, h = h, (al * h - np.sqrt(k - 1.0) * h_prev) / np.sqrt(float(k))
        g[m] = (-1.0) ** m * phi * h / np.sqrt(m * (m - 1.0))
    return g


def relu_var(mu, v):
    s = np.sqrt(v); al = mu / s; Phi = ndtr(al); phi = np.exp(-0.5 * al * al) / SQ2PI
    m1 = s * (al * Phi + phi); m2 = v * ((al * al + 1.0) * Phi + al * phi)
    return m2 - m1 * m1


def KG(mu, C):
    v = np.diag(C).copy(); s = np.sqrt(v); al = mu / s
    R = C / np.outer(s, s); np.fill_diagonal(R, 0.0)
    g = gcoef(al) * s[None, :]
    K = np.zeros_like(C); P = np.ones_like(C)
    for m in range(1, M + 1):
        P *= R
        K += np.outer(g[m], g[m]) * P
    K[np.diag_indices_from(K)] = relu_var(mu, v)
    return K


def qdiag(W, K):
    return np.einsum("ij,ij->i", W @ K, W)


def stats(name, comp, ref, A=None, B=None):
    """mean and rms of comp/ref; noise-free rms from halves; share of the variance-error energy it carries."""
    r = comp / ref
    out = f"{name}: mean {r.mean():+.2e} rms {np.sqrt(np.mean(r * r)):.2e}"
    if A is not None:
        ra, rb = A / ref, B / ref
        out += f" (noise-free rms {np.sqrt(max(np.mean(ra * rb), 0.0)):.2e})"
    return out


for net in nets:
    cd = np.load(f"{loc}/chaindump_{net}.npz")
    W = np.load(f"{loc}/W_off{net}.npy").astype(np.float64)
    T = {h: np.load(f"{mcd}/mc2_off{net}_{h}.npz") for h in ("full", "h0", "h1")}
    out = cd["out"]
    print(f"=== network {net} (all per-neuron quantities divided by the true variance at the target layer)", flush=True)
    rows = []
    for s in range(0, 15):
        t = s + 1
        mu_c = np.zeros(1024) if s == 0 else W[s] @ out[s - 1]
        Cc = cd["C"][s].astype(np.float64); Cc = 0.5 * (Cc + Cc.T); np.fill_diagonal(Cc, cd["var"][s])
        if s == 0:
            Cc = W[0] @ W[0].T
        Kc = KG(mu_c, Cc)
        res = {}
        for h in ("full", "h0", "h1"):
            mu_t = T[h]["mu"][s].astype(np.float64)
            Ct = T[h]["cov"][s].astype(np.float64); Ct = 0.5 * (Ct + Ct.T)
            Kt = KG(mu_t, Ct)
            vt_t = T[h]["var"][t].astype(np.float64)
            g = vt_t - qdiag(W[t], Kt)
            p = qdiag(W[t], Kc - Kt)
            res[h] = dict(g=g, p=p, vt=vt_t)
            if h == "full":
                # split of p: chain means only / chain variances only / chain off-diagonal only (others true)
                Cv = Ct.copy(); np.fill_diagonal(Cv, np.diag(Cc))
                Co = Cc.copy(); np.fill_diagonal(Co, np.diag(Ct))
                res[h]["p_mu"] = qdiag(W[t], KG(mu_c, Ct) - Kt)
                res[h]["p_var"] = qdiag(W[t], KG(mu_t, Cv) - Kt)
                res[h]["p_off"] = qdiag(W[t], KG(mu_t, Co) - Kt)
        vc_t = cd["var"][t]
        x = vc_t - qdiag(W[t], Kc)
        F, A, B = res["full"], res["h0"], res["h1"]
        ref = F["vt"]
        dv = vc_t - ref
        print(f" layer {t:2d}: dv mean {np.mean(dv / ref):+.2e} rms {np.sqrt(np.mean((dv / ref) ** 2)):.2e} | "
              f"check p+x-g-dv rms {np.sqrt(np.mean(((F['p'] + x - F['g'] - dv) / ref) ** 2)):.1e}", flush=True)
        print("    " + stats("g (true defect)", F["g"], ref, A["g"], B["g"]), flush=True)
        print("    " + stats("x (chain corr.)", x, ref), flush=True)
        print("    " + stats("x - g (injected)", x - F["g"], ref, x - A["g"], x - B["g"]), flush=True)
        print("    " + stats("p (propagated)", F["p"], ref, A["p"], B["p"]) + " = "
              + ", ".join(stats(k, F[k], ref).split(":")[1].split(" rms")[0].replace(" mean", k + " mean")
                          for k in ("p_mu", "p_var", "p_off")), flush=True)
        # energy shares of dv carried by p and by (x - g), noise-corrected through the halves
        da, db = vc_t - A["vt"], vc_t - B["vt"]
        e_dv = np.mean((da / ref) * (db / ref))
        sh_p = np.mean((A["p"] / ref) * (db / ref) + (B["p"] / ref) * (da / ref)) / (2 * e_dv)
        sh_i = np.mean(((x - A["g"]) / ref) * (db / ref) + ((x - B["g"]) / ref) * (da / ref)) / (2 * e_dv)
        print(f"    projection shares of the noise-free dv energy: propagated {sh_p:+.2f}, injected {sh_i:+.2f}", flush=True)
        rows.append([t, np.mean(dv / ref), np.mean(F["g"] / ref), np.mean(x / ref), np.mean(F["p"] / ref),
                     np.mean(F["p_mu"] / ref), np.mean(F["p_var"] / ref), np.mean(F["p_off"] / ref), sh_p, sh_i])
    np.save(f"var_decomp_{net}.npy", np.array(rows))
    print(flush=True)
