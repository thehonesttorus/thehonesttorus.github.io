# Structure of the adopted chain's pre-activation variance error (note XLII section 7).
#   python ana_var.py LOCDIR MC2DIR NET [NET ...]
# LOCDIR holds chaindump_{n}.npz and W_off{n}.npy; MC2DIR holds mc2_off{n}_{full,h0,h1}.npz (16M-input Monte Carlo,
# two independent halves). Per layer: the relative variance error r_i = var_chain/var_true - 1, its rms against the
# Monte Carlo noise (from the halves), and the share of its noise-free energy explained by candidate per-neuron shapes:
#   gain:    mu_i^2 / var_i            (the radial/gain mode's exact contribution to Var(h_i) is mu_i^2 (n/E|X|^2 - 1))
#   alpha:   mu_i / sigma_i
#   kurt:    true kappa4_i / var_i^2   (excess kurtosis)    skew: true kappa3_i / sigma_i^3
#   coll:    share of var_i carried by the top-16 eigenmodes of the true covariance C_l
#   and the joint fit of all five (noise-corrected R^2 = explained / noise-free energy).
import sys, numpy as np
loc, mcd = sys.argv[1], sys.argv[2]; nets = [int(x) for x in sys.argv[3:]]
for net in nets:
    cd = np.load(f"{loc}/chaindump_{net}.npz")
    F, A, B = (np.load(f"{mcd}/mc2_off{net}_{h}.npz") for h in ("full", "h0", "h1"))
    W = np.load(f"{loc}/W_off{net}.npy").astype(np.float64)
    print(f"=== network {net}: layer | rms rel. var error (chain) | MC noise rms | noise-free rms | "
          f"share explained: gain alpha kurt skew coll | joint | mean rel. error")
    for l in range(1, 16):
        vt = F["var"][l].astype(np.float64); va = A["var"][l].astype(np.float64); vb = B["var"][l].astype(np.float64)
        vc = cd["var"][l]
        mu = F["mu"][l].astype(np.float64)
        r = vc / vt - 1.0
        ra = vc / va - 1.0; rb = vc / vb - 1.0
        noise = 0.5 * np.mean((ra - rb) ** 2)                       # ~ var of the full-sample noise in r
        sig_e = max(np.mean(np.mean(ra * rb) - np.mean(ra) * np.mean(rb)), 1e-30)
        feats = {
            "gain": mu * mu / vt,
            "alpha": mu / np.sqrt(vt),
            "kurt": F["k4"][l].astype(np.float64) / vt ** 2,
            "skew": F["k3"][l].astype(np.float64) / vt ** 1.5,
        }
        C = F["cov"][l].astype(np.float64); C = 0.5 * (C + C.T); np.fill_diagonal(C, vt)
        ev, V = np.linalg.eigh(C); V = V[:, ::-1]; ev = ev[::-1]
        feats["coll"] = (V[:, :16] ** 2 @ ev[:16]) / vt
        shares = []
        rc = rb - rb.mean()       # explained share uses one half as response and the other half's features? keep simple:
        for k, f in feats.items():
            fc = f - f.mean()
            # noise-corrected share: cov(ra, f) cov(rb, f) / (var(f) * signal energy)
            ca = np.mean((ra - ra.mean()) * fc); cb = np.mean((rb - rb.mean()) * fc)
            shares.append(ca * cb / (np.mean(fc * fc) * sig_e))
        X = np.stack([f - f.mean() for f in feats.values()], 1)
        ba = np.linalg.lstsq(X, ra - ra.mean(), rcond=None)[0]; bb = np.linalg.lstsq(X, rb - rb.mean(), rcond=None)[0]
        joint = np.mean((X @ ba) * (X @ bb)) / sig_e
        print(f"  {l:2d} | {np.sqrt(np.mean(r ** 2)):.2e} | {np.sqrt(noise):.2e} | {np.sqrt(sig_e):.2e} | "
              + " ".join(f"{s:+.2f}" for s in shares) + f" | {joint:+.2f} | {np.mean(r):+.2e}", flush=True)
