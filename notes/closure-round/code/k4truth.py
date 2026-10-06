# Fourth-cumulant slices of the PRE-activation y_l from the Monte Carlo accumulators, projected on the scale-mixture shapes:
# k4 on 3 var^2, K22 (off) on var var^T, K31 (off) on 3 d(var) C_off and on C_off; coefficient = the implied g per slice.
import numpy as np, sys
from pairvar import y_cumulants
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
D = np.load(f"mc_cum_off{net}.npz"); L = D["s1y"].shape[0]
try: gm = np.load(f"k3channel_off{net}.npz")["gm"]
except Exception: gm = None
LAM = [4.9895e-03, 8.0876e-03, 9.8291e-03, 1.0549e-02, 1.0851e-02, 1.0828e-02, 1.0589e-02, 1.0048e-02, 9.6483e-03, 9.1720e-03, 8.7730e-03, 8.3938e-03, 8.0555e-03, 7.6770e-03, 7.2588e-03, 7.2588e-03]
proj = lambda A, B: (float(np.sum(A * B) / np.sum(B * B)), float(1 - np.sum((A - np.sum(A * B) / np.sum(B * B) * B)**2) / np.sum(A * A)))
print(f"net {net}: implied g per kappa4 slice of the pre-activation (MC, T = {int(D['T'])}); last column: chain table LAM[l]*0.95 for the (3,1) slice, and E10's g_l")
print(" l | g from k4 diag (R^2) | g from K22 off (R^2) | g from K31 off vs 3 d(var) C (R^2) | K31 off vs C: coef (R^2) | mean var | 3 g var | LAM*0.95 | g_l E10")
for l in range(1, L):
    Y = y_cumulants(D, l); var = np.diag(Y["S"]); C_off = Y["S"] - np.diag(var); n = len(var)
    k4, K31, K22 = Y["k4"], Y["K31"].copy(), Y["K22"].copy(); np.fill_diagonal(K31, 0); np.fill_diagonal(K22, 0)
    g4, r4 = proj(k4, 3 * var * var); g22, r22 = proj(K22, np.outer(var, var) - np.diag(var * var)); g31, r31 = proj(K31, 3 * var[:, None] * C_off); c31, rc = proj(K31, C_off)
    gl = gm[l] if gm is not None else float("nan")
    print(f" {l:2d} | {g4:+.4f} ({r4:.2f}) | {g22:+.4f} ({r22:.2f}) | {g31:+.4f} ({r31:.2f}) | {c31:+.4f} ({rc:.2f}) | {var.mean():.3f} | {3*g31*var.mean():.4f} | {LAM[l]*0.95:.4f} | {gl:.4f}", flush=True)
