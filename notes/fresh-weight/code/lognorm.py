# Note XLIV section 9 (P1): is the radial gain a saturating walk, as the stage-13 note predicts, or mean-reverting?
#   python -I lognorm.py LOCDIR NET N
# Forward N Gaussian inputs through the network (h_l = relu(W_l h_(l-1)), width n) and record, for the post-activation
# squared norm rho_l = |h_l|^2 across inputs:
#   V_l = Var_x log rho_l, the increments Delta_l = log rho_l - log rho_(l-1), Var Delta_l, corr(Delta_l, Delta_(l+1)),
#   the AR(1) slope b_l = Cov(Delta_(l+1), log rho_l) / Var(log rho_l), and the leading law
#   Var Delta_l ~ 4 theta_(l-1)^2 / n with theta from the arc-cosine recursion cos theta_(l+1) = k1(cos theta_l),
#   theta_0 = pi/2 (independent inputs), k1(r) = (sqrt(1 - r^2) + (pi - arccos r) r) / pi.
# Also the empirical angle between independent inputs' representations, cos theta_l = <h_l(x), h_l(x')> / (|h||h'|).
import sys, numpy as np
loc, net, N = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
W = np.load(f"{loc}/W_off{net}.npy")
n = W.shape[1]
rng = np.random.default_rng(12345)
H = rng.standard_normal((N, n), dtype=np.float32)
L = W.shape[0]
logr = np.zeros((L + 1, N)); cosang = np.zeros(L + 1)
logr[0] = np.log(np.sum(H.astype(np.float64) ** 2, axis=1))
cosang[0] = 0.0
half = N // 2
for l in range(L):
    H = H @ W[l].T
    np.maximum(H, 0.0, out=H)
    logr[l + 1] = np.log(np.sum(H.astype(np.float64) ** 2, axis=1))
    A, B = H[:half].astype(np.float64), H[half:2 * half].astype(np.float64)
    cosang[l + 1] = float(np.mean(np.sum(A * B, axis=1) / np.sqrt(np.sum(A * A, axis=1) * np.sum(B * B, axis=1))))
V = logr.var(axis=1)
D = np.diff(logr, axis=0)                               # D[l] = log rho_(l+1) - log rho_l, l = 0..L-1
th = [np.pi / 2]
k1 = lambda r: (np.sqrt(1 - r * r) + (np.pi - np.arccos(r)) * r) / np.pi
c = 0.0
for l in range(L):
    c = k1(c); th.append(np.arccos(min(c, 1.0)))
print(f"=== network {net}, N = {N}, n = {n}")
print("layer  V_l/4     Var Delta   4 theta^2/n  ratio   corr(Delta_l, Delta_l+1)  AR1 slope b   cos(theta) measured / recursion")
for l in range(1, L + 1):
    vd = D[l - 1].var()
    law = 4.0 * th[l - 1] ** 2 / n
    cc = np.corrcoef(D[l - 1], D[l])[0, 1] if l < L else float("nan")
    b = np.cov(D[l], logr[l])[0, 1] / logr[l].var() if l < L else float("nan")
    print(f"  {l:2d}  {V[l] / 4:8.5f}  {vd:9.3e}  {law:9.3e}  {vd / law:5.2f}   {cc:+7.3f}                 {b:+7.3f}      "
          f"{cosang[l]:.4f} / {np.cos(th[l]):.4f}")
