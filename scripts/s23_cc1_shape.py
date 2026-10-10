"""Shape versus moments at the final layer (python scripts/s23_cc1_shape.py KIK1 DATA NET ...).
CC1's final pre-activation mixture z_k | node ~ N(t_k + Wb_k dc_q, rho_q vt_k). Compared:
  (a) CC1 as computed;  (b) the mixture readout re-centred on the exact Monte Carlo (t_k, var_k) [vt_k := var_k - Wb_k^2 Vc];
  (c) the Gaussian readout at exact (t_k, var_k);  (d) the Hermite readout at exact (t, var, k3, k4) (same-sample floor).
Also CC1's own t and var errors, and its mixture skewness against the Monte Carlo skewness."""
import sys, os, numpy as np
from scipy.special import ndtr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.cc1 import CC1, trunc_moments
from scripts.kik_merge import pool_net
from scripts.kik_ladder import hermite_means
rms = lambda x: np.sqrt(np.mean(np.square(x)))
K1, D = sys.argv[1], sys.argv[2]
for net in map(int, sys.argv[3:]):
    W = [np.asarray(w, dtype=np.float64) for w in np.load(f"{D}/W_off{net}.npy")]
    T = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
    a = pool_net(K1, net); N = a["N"]; P = a["pw_top"] / N; t, q, m3, m4 = P[:, 0], P[:, 1], P[:, 2], P[:, 3]; var = q - t * t
    m = CC1(Q=7); out, _ = m.run(W); pre = m.last["pre"]; nodes = pre["nodes"]
    Vc_prev = sum(w * dc * dc for w, dc, r in nodes)
    def mix_readout(tz, vt):
        y = 0
        for w, dc, r in nodes:
            sd = np.sqrt(r * vt); al = (tz + pre["Wb"] * dc) / sd; y = y + w * sd * trunc_moments(al, 1)[1]
        return y
    t_cc, v_cc = pre["mu_z"], pre["vt"] + pre["Wb"] ** 2 * Vc_prev
    k3_cc = sum(w * (3 * (pre["Wb"] * dc) * r * pre["vt"] + (pre["Wb"] * dc) ** 3) for w, dc, r in nodes) / v_cc ** 1.5
    k3_mc = (m3 - 3 * t * q + 2 * t ** 3) / var ** 1.5
    ya = out[-1]; yb = mix_readout(t, var - pre["Wb"] ** 2 * Vc_prev)
    g, h3, h34, *_ = hermite_means(t, q, m3, m4)
    mc = a["pos"][-1] / N
    print(f"net {net}: (a) CC1 {rms(ya - T[-1]):.2e} | (b) mixture at exact (t, var) {rms(yb - T[-1]):.2e} [same-sample {rms(yb - mc):.2e}] |"
          f" (c) Gaussian at exact {rms(g - T[-1]):.2e} [same-sample {rms(g - mc):.2e}] | (d) Hermite k3+k4 [same-sample {rms(h34 - mc):.2e}]")
    print(f"   CC1 moment errors: t {rms(t_cc - t):.2e}, var {rms(v_cc - var):.2e} (rel {rms((v_cc - var) / var):.2e});"
          f" mixture skewness rms {rms(k3_cc):.3f} vs MC {rms(k3_mc):.3f}, corr {np.corrcoef(k3_cc, k3_mc)[0, 1]:.3f}, residual rms {rms(k3_mc - k3_cc):.3f}", flush=True)
