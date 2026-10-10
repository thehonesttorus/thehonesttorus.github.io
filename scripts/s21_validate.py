"""Stage 21 at n = 1024 on the real networks, from the Gaussian-part trajectory (pass 1) and the pooled kik1 contractions.

  python scripts/s21_validate.py DATA RESULTS NET [NET ...]

Per layer: realised rho_(l-1) = |mu|^2 / Tr S^u against the arc-cosine trajectory; Var_j alpha_j against rho/(1-rho);
realised gap g_l = 2 mean Phi(alpha)^2 against f'(rho_realised) and the trajectory value; level-k contraction of
random incoherent directions through the real (W_l, Phi(alpha_l)) for k = 1, 2, 3 (product kernel), and for k = 2 with
the real pair kernel K_l and with the exact tangent (wall coupling included); flatness spectrum of Pi_perp (F_l - I) Pi_perp.
Per net: localisation radii from the realised gaps; Wick-orthogonality split (coherent function of a_k = w_k . mu_hat
versus incoherent residual) of the true final means, of the closure's errors and of the final kappa_3, kappa_4."""
import sys, os, numpy as np
from scipy.special import ndtr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest import gate_channel as gc

rms = lambda x: np.sqrt(np.mean(np.square(x)))


def coherent_split(y, a, deg=8):
    """Fit Ybar(a) by polynomial regression in a/sd(a) (the n output neurons are n i.i.d. probes); return the
    coherent share Var(Ybar)/Var(y), the residual, and the fit."""
    x = a / a.std(); V = np.vander(x, deg + 1); c = np.linalg.lstsq(V, y, rcond=None)[0]; fit = V @ c
    return 1 - np.var(y - fit) / np.var(y), y - fit, fit


def analyse(D, R, net, rng):
    W = [np.asarray(w, dtype=np.float64) for w in np.load(f"{D}/W_off{net}.npy")]; L, n = len(W), W[0].shape[0]
    truth = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
    st = gc.closure_states(W); rho_t, g_t = gc.trajectory(L)
    print(f"\n=== net {net} ===")
    print(" l | rho_prev real / traj | Var(alpha) vs rho/(1-rho) | g real  f'(rho real)  traj | level-k/g^k k=1 2 3 | pair-kernel k=2, exact tangent | flatness mean sd [min, max]")
    greal = []
    for l in range(L):
        s = st[l]; r = s["rho_prev"]; a = s["alpha"]; Pa = s["Pa"]; g = 2 * np.mean(Pa ** 2); greal.append(g)
        mh = s["mu_prev"] / max(np.linalg.norm(s["mu_prev"]), 1e-300) if l > 0 else None
        V = rng.standard_normal((n, 16))
        if mh is not None: V -= np.outer(mh, mh @ V)
        lk = [np.mean(gc.level_gain(W[l], Pa, V, k)) / g ** k for k in (1, 2, 3)]
        # level 2 with the real pair kernel and with the exact tangent, X = v v^T - (|v|^2/n) I projected (incoherent, traceless)
        k2, k2f = [], []
        if l > 0:
            Kc = gc.cut_kernel(s["mu"], s["C"])
            for j in range(4):
                v = V[:, j]; X = np.outer(v, v); dS = W[l] @ X @ W[l].T
                k2.append(np.sum((Kc * dS) ** 2) / np.sum(X * X)); k2f.append(np.sum(gc.full_tangent(s["mu"], s["C"], Kc, dS) ** 2) / np.sum(X * X))
        fl = gc.flatness(gc.channel_dual(W[l], Pa), s["mu_prev"]) if l > 0 else np.zeros(1)
        print(f"{l + 1:2d} | {r:.3f} / {rho_t[l]:.3f} | {np.var(a):7.3f} vs {r / (1 - r) if r < 1 else np.inf:7.3f} | {g:.3f}  {gc.gap_closed(r):.3f}  {g_t[l]:.3f} |"
              f" {lk[0]:.3f} {lk[1]:.3f} {lk[2]:.3f} | " + (f"{np.mean(k2) / g ** 2:.3f}, {np.mean(k2f) / g ** 2:.3f}" if k2 else "  -  ,   -  ") +
              f" | {fl.mean():+.4f} {fl.std():.3f} [{fl.min():+.2f}, {fl.max():+.2f}]")
    greal = np.array(greal)
    print("  localisation radii from the realised gaps (eta = 0.1 / 0.05 / 0.02 / 0.01): " +
          "; ".join(f"k={k}: " + "/".join(str(gc.radii(greal, k, e) or ">16") for e in (0.1, 0.05, 0.02, 0.01)) for k in (1, 2, 3, 4)) +
          f"; level-1 amplitude from layer 2: {np.sqrt(np.prod(greal[2:])):.3f}")
    # Wick-orthogonality split at the final layer: a_k = w_k . mu_hat with mu_hat the true penultimate mean direction
    mh = truth[-2] / np.linalg.norm(truth[-2]); a = W[-1] @ mh; y = truth[-1]
    sh, eta, fit = coherent_split(y, a)
    ec = st[-1]["m"] - y; she, _, _ = coherent_split(ec, a)
    T = lambda r_: [np.mean(r_ * p) / (np.std(r_ * p) / np.sqrt(n)) for p in (np.ones(n), a / a.std(), (a / a.std()) ** 2 - 1)]
    line = f"  Wick split, final layer: coherent share of the true means {sh:.4f} (residual rms {rms(eta):.2e}); of the closure error {she:.3f} (closure rms err {rms(ec):.2e})"
    kf = f"{R}/kikm_{net}.npz"
    if os.path.exists(kf):
        z = np.load(kf); s3, _, _ = coherent_split(z["c3"], a); s4, _, _ = coherent_split(z["c4"], a)
        line += f"; of kappa3 {s3:.3f}, of kappa4 {s4:.3f}"
    print(line)
    # T_psi on the closure's residual against the truth's own coherent fit (detects coherent leakage)
    print("  T_psi z-scores (psi = 1, He1, He2) of closure - Ybar_truth: " + " ".join(f"{t:+.1f}" for t in T(st[-1]["m"] - fit)) +
          "; of the truth residual itself: " + " ".join(f"{t:+.1f}" for t in T(eta)))
    np.savez(f"{R}/s21v_{net}.npz", g=greal, a=a, fit=fit)


if __name__ == "__main__":
    D, R = sys.argv[1], sys.argv[2]; rng = np.random.default_rng(2121)
    for net in map(int, sys.argv[3:]): analyse(D, R, net, rng)
