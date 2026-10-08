# Signed budget of a paired chain comparison (checkpoint J, protocol steps 1-2; note XXXVI section 3h).
#   python sbudget.py NET MCPREFIX DUMP0 DUMP1
# DUMP0 / DUMP1: free-running dumps (oracle_one.py NET dump none, KEEP_EXTRA=mu,pk1v,K2v) of the baseline and repaired chain.
# Step 1 (exact bookkeeping of the arrays, per layer):
#   pre-activation: Delta var = Delta S2 - (mu_1 + mu_0) Delta mu,      S2 = var + mu^2
#   post-activation: Delta v = Delta q - (m_1 + m_0) Delta m,            q = v + m^2
#   and for each run its error against Monte Carlo truth split the same way: e_var = e_S2 - (mu_hat + mu*) e_mu.
# Step 2 (signed output budget): e_r = out_r - truth at the last layer, d = e_1 - e_0;
#   |e_1|^2 - |e_0|^2 = 2 <e_0, d> + |d|^2 = <d, e_0 + e_1>   (checked).
#   Tangent ledger with the Jacobian at Monte Carlo truth (as errbudget.py): d = sum_X tau_X + rho, tau_X the propagated
#   change of channel X's mean-map injections (var, kappa3, kappa4 diagonals), rho the remainder; the allocation
#   a_X = <tau_X, e_0 + e_1> sums with a_rho to the exact change of n x MSE. The var channel is split per layer into its
#   raw-second-moment part (coef_var Delta S2) and its recentering part (-coef_var (mu_1 + mu_0) Delta mu).
import sys, math, numpy as np
net, pre, d0p, d1p = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64); L, n, _ = Wcol.shape
mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
F = np.load(f"{pre}_full.npz"); A0, A1 = np.load(d0p), np.load(d1p)
mu, var = F["mu"].astype(np.float64), F["var"].astype(np.float64)
sig = np.sqrt(var); al = mu / sig; ph = np.exp(-al * al / 2) / math.sqrt(2 * math.pi); Ph = 0.5 * (1 + np.vectorize(math.erf)(al / math.sqrt(2)))
coef = {"var": ph / (2 * sig), "k3": -al * ph / (6 * var), "k4": (al * al - 1) * ph / (24 * var * sig)}
cname = {"var": "var", "k3": "D3", "k4": "g4row"}
g = lambda A, k, l: A[f"{k}_{l}"].astype(np.float64) if f"{k}_{l}" in A.files else None


def zmean(A, l):
    # the dump's "mu" is the post-activation mean except at the last layer; the pre-activation mean is W_l m_(l-1), exactly
    m = g(A, "pk1v", l - 1)
    return None if m is None else Wcol[l] @ m


def prop(l, v):
    for k in range(l + 1, L):
        v = Ph[k] * (Wcol[k] @ v)
    return v


out0, out1 = g(A0, "pk1v", L - 1), g(A1, "pk1v", L - 1)
e0, e1 = out0 - mt[-1], out1 - mt[-1]; d = e1 - e0; s = e0 + e1
dM = float(e1 @ e1 - e0 @ e0)
print(f"net {net}: {d0p} -> {d1p}")
print(f"  MSE {e0 @ e0 / n:.5e} -> {e1 @ e1 / n:.5e} ({100 * dM / (e0 @ e0):+.2f}%); check 2<e0,d>+|d|^2 = {2 * e0 @ d + d @ d:.4e} vs {dM:.4e}; "
      f"<e0,d>/|e0||d| {e0 @ d / np.linalg.norm(e0) / np.linalg.norm(d):+.3f}, |d|/|e0| {np.linalg.norm(d) / np.linalg.norm(e0):.3f}")
# ---- step 1: exact bookkeeping per layer ----
print("  step 1, per layer (relative norms; pre-activation | post-activation):")
print("    layer: |dvar|/|var*|  |dS2|  |rec|  corr(dS2,rec) | err0 var: S2-part rec-part | err1 var: S2-part rec-part || |dv|/|v*| |dq| |rec| | err0 v: q rec | err1 v: q rec")
for l in range(1, L):
    m0z, m1z, v0z, v1z = zmean(A0, l), zmean(A1, l), g(A0, "var", l), g(A1, "var", l)
    if any(x is None for x in (m0z, m1z, v0z, v1z)):
        continue
    vt, mut = var[l], mu[l]
    dvar = v1z - v0z; dmu = m1z - m0z; dS2 = dvar + (m1z + m0z) * dmu; rec = -(m1z + m0z) * dmu
    r = lambda x, ref: np.linalg.norm(x) / np.linalg.norm(ref)
    ev = lambda vz, mz: (vz - vt + (mz * mz - mut * mut), -(mz + mut) * (mz - mut))   # (e_S2, recentering) with e_var = sum
    e0S, e0R = ev(v0z, m0z); e1S, e1R = ev(v1z, m1z)
    row = (f"    {l:2d}: {r(dvar, vt):.2e} {r(dS2, vt):.2e} {r(rec, vt):.2e} {np.corrcoef(dS2, rec)[0, 1]:+.2f} | "
           f"{r(e0S, vt):.2e} {r(e0R, vt):.2e} | {r(e1S, vt):.2e} {r(e1R, vt):.2e}")
    my0, my1, vy0, vy1 = g(A0, "pk1v", l - 1), g(A1, "pk1v", l - 1), g(A0, "K2v", l - 1), g(A1, "K2v", l - 1)
    if l - 1 >= 0 and my0 is not None and vy0 is not None and my1 is not None and vy1 is not None:
        vyt, myt = F["var_y"][l - 1].astype(np.float64), F["mu_y"][l - 1].astype(np.float64)
        dv = vy1 - vy0; dm = my1 - my0; dq = dv + (my1 + my0) * dm; recy = -(my1 + my0) * dm
        eq0, er0 = vy0 - vyt + (my0 ** 2 - myt ** 2), -(my0 + myt) * (my0 - myt)
        eq1, er1 = vy1 - vyt + (my1 ** 2 - myt ** 2), -(my1 + myt) * (my1 - myt)
        row += (f" || y{l - 1}: {r(dv, vyt):.2e} {r(dq, vyt):.2e} {r(recy, vyt):.2e} | {r(eq0, vyt):.2e} {r(er0, vyt):.2e} | "
                f"{r(eq1, vyt):.2e} {r(er1, vyt):.2e}")
    print(row, flush=True)
# ---- step 2: signed output budget, tangent ledger at Monte Carlo truth ----
tau = {X: np.zeros(n) for X in ("var", "k3", "k4")}; per = {}
tauS2 = np.zeros(n); tauRec = np.zeros(n)
for l in range(1, L):
    for X in ("var", "k3", "k4"):
        a, b = g(A0, cname[X], l), g(A1, cname[X], l)
        if a is None or b is None:
            continue
        p = prop(l, coef[X][l] * (b - a)); tau[X] += p; per[(X, l)] = p
    m0z, m1z = zmean(A0, l), zmean(A1, l)
    if m0z is not None and m1z is not None and ("var", l) in per:
        pr = prop(l, -coef["var"][l] * (m1z + m0z) * (m1z - m0z)); tauRec += pr; per[("rec", l)] = pr
        per[("S2", l)] = per[("var", l)] - pr; tauS2 += per[("S2", l)]
rho = d - sum(tau.values())
alloc = {X: float(tau[X] @ s) for X in tau}; alloc["remainder"] = float(rho @ s)
print(f"  step 2: signed allocation of n x Delta MSE = {dM:.4e} (sum check {sum(alloc.values()):.4e}); |rho|/|d| {np.linalg.norm(rho) / np.linalg.norm(d):.3f}")
for k, v in alloc.items():
    print(f"    {k:9s}: {v:+.4e}  ({100 * v / (e0 @ e0):+.2f}% of baseline n x MSE)")
print(f"    var split: raw second moment {float(tauS2 @ s):+.4e}, recentering {float(tauRec @ s):+.4e}")
print("  per source layer, allocation <prop(l, .), e0 + e1> as % of baseline (var = S2 + rec | k3 | k4):")
for l in range(1, L):
    f = lambda key: 100 * float(per[key] @ s) / float(e0 @ e0) if key in per else float("nan")
    print(f"    layer {l:2d}: var {f(('var', l)):+6.2f} = S2 {f(('S2', l)):+6.2f} + rec {f(('rec', l)):+6.2f} | k3 {f(('k3', l)):+6.2f} | k4 {f(('k4', l)):+6.2f}")
# ---- step 3: where the raw second moment change of z_l enters (first entry into the second-moment block) ----
# S2_z,l = sum_a H_la q_(l-1),a + sum_(a != b) W_la W_lb (C_ab + m_a m_b);   the unary route, to first order at the truth:
#   dq = 2 m dmu + Phi dvar + (phi / 3s) dkappa3 - (alpha phi / 12 s^2) dkappa4     (checkpoint J, proposition 3.1)
# Pieces of Delta S2_z,l: k3 -> q (direct cubic entry), k4 -> q (direct quartic entry), (mu, var) -> q (inside the block),
# q residual (the chain's actual Delta q minus the three, i.e. its closure's departure from the Gaussian coefficients),
# and the off-diagonal part (Delta S2 - H Delta q). Each is allocated by <prop(l, coef_var * piece), e0 + e1>.
print("  step 3: entries into Delta S2 of z_l, allocation as % of baseline n x MSE; and the unary closure check")
print("    layer: k3->q  k4->q  (mu,var)->q  q-resid  offdiag || Delta q check: |pred-actual|/|actual| corr")
tot = {}
for l in range(2, L):
    mz0, mz1 = zmean(A0, l - 1), zmean(A1, l - 1)
    need = [g(A0, k, l - 1) for k in ("var", "D3", "g4row", "pk1v", "K2v")] + [g(A1, k, l - 1) for k in ("var", "D3", "g4row", "pk1v", "K2v")]
    v0z, v1z = g(A0, "var", l), g(A1, "var", l)
    if any(x is None for x in need) or mz0 is None or v0z is None or v1z is None:
        continue
    va0, D30, g40, m0, vy0, va1, D31, g41, m1, vy1 = need
    mu_, s_ = mu[l - 1], sig[l - 1]; a_ = mu_ / s_; c_ = ph[l - 1]; P_ = Ph[l - 1]
    my_t = F["mu_y"][l - 1].astype(np.float64)
    dq_act = (vy1 + m1 * m1) - (vy0 + m0 * m0)
    p_k3 = c_ / (3 * s_) * (D31 - D30)
    p_k4 = -a_ * c_ / (12 * s_ * s_) * (g41 - g40)
    p_mv = 2 * my_t * (mz1 - mz0) + P_ * (va1 - va0)
    p_res = dq_act - p_k3 - p_k4 - p_mv
    H = Wcol[l] * Wcol[l]
    m0z_l, m1z_l = zmean(A0, l), zmean(A1, l)
    dS2 = (v1z - v0z) + (m1z_l + m0z_l) * (m1z_l - m0z_l)
    p_off = dS2 - H @ dq_act
    cv = coef["var"][l]
    al_ = {k: 100 * float(prop(l, cv * x) @ s) / float(e0 @ e0) for k, x in
           (("k3", H @ p_k3), ("k4", H @ p_k4), ("mv", H @ p_mv), ("res", H @ p_res), ("off", p_off))}
    for k, v in al_.items():
        tot[k] = tot.get(k, 0.0) + v
    pred = p_k3 + p_k4 + p_mv
    print(f"    {l:2d}: {al_['k3']:+6.2f} {al_['k4']:+6.2f} {al_['mv']:+6.2f} {al_['res']:+6.2f} {al_['off']:+6.2f} || "
          f"{np.linalg.norm(pred - dq_act) / max(np.linalg.norm(dq_act), 1e-300):.3f} {np.corrcoef(pred, dq_act)[0, 1]:+.3f}", flush=True)
print("    total: " + " ".join(f"{k} {v:+.2f}" for k, v in tot.items()))
