# Note XLIX section 9: the second-chaos four-loop of the (3,1) slice, kappa(Q_i,Q_i,Q_i,Q_j) = 48 tr(A_i^3 A_j),
# A_i = (1/2) Lambda D_i Lambda^T, D_i = diag(J2_s V_is), so the term is 3 tr(B_i^3 B_j), B_i = Lambda D_i Lambda^T.
# Evaluated for NROWS sampled active rows i and all j, then combined with the ladder, cross-age pairs and secondary
# stars on that subset.   python chaos_loop.py NET L MCDIR WDIR LADDER_NPZ [NROWS]
import sys, os, re, math, numpy as np
net, L, MC, WD, LAD = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
NROWS = int(sys.argv[6]) if len(sys.argv) > 6 else 64
# reuse chaos_qq.py up to the terms (it defines blocks, Ub, K_cross, K_sec, mask, w, t0, t1, tf, cum)
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "chaos_qq.py")).read()
src = src[:src.index("yc, ys = K_cross[mask]")]
sys.argv = ["chaos_qq.py", str(net), str(L), MC, WD, LAD]
exec(src)
Lam = np.concatenate([Lb for (Lb, _, _) in blocks], axis=1)          # n_in x S
J2all = np.concatenate([Jb for (_, Jb, _) in blocks])
Vall = np.concatenate([Vb for (_, _, Vb) in blocks], axis=1)          # n x S
Aall = Vall * J2all[None, :]                                          # a_is = J2_s V_is
rows_act = np.nonzero(mask.any(axis=1))[0]
rng = np.random.default_rng(7 + L)
rows = np.sort(rng.choice(rows_act, size=min(NROWS, len(rows_act)), replace=False))
K_loop = np.zeros((len(rows), n))
for r, i in enumerate(rows):
    Bi = (Lam * Aall[i][None, :]) @ Lam.T                              # n_in x n_in
    B3 = Bi @ Bi @ Bi
    dv = np.einsum("sv,sv->v", Lam, B3 @ Lam)                          # diag(Lam^T B_i^3 Lam), length S
    K_loop[r] = 3.0 * (Aall @ dv)                                      # tr(B_i^3 B_j) = sum_v a_jv dv_v
# subset of the mask on the sampled rows, in the mask's own ordering
idx_i, idx_j = np.nonzero(mask)
pos = {i: r for r, i in enumerate(rows)}
sub = np.array([k for k, i in enumerate(idx_i) if i in pos])
yl = np.array([K_loop[pos[idx_i[k]], idx_j[k]] for k in sub])
ws, a0, a1, af = w[sub], t0[sub], t1[sub], tf[sub]
def r2s(y):
    ip = lambda A, B: float(np.sum(A * B * ws * ws))
    return 1.0 - ip(a0 - y, a1 - y) / ip(a0, a1), ip(y, af) / ip(y, y)
yc, ysec, lad = K_cross[mask][sub], K_sec[mask][sub], cum[sub]
print(f"=== network {net}, (3,1) slice of layer {L + 1}: four-loop on {len(rows)} sampled rows ({len(sub)} entries)")
for name, y in (("ladder", lad), ("ladder + cross pairs + sec. stars", lad + yc + ysec), ("four-loop alone", yl),
                ("ladder + pairs + stars + four-loop", lad + yc + ysec + yl)):
    R, s = r2s(y)
    print(f"  {name:36s} R2 {100 * R:6.1f}%  scale {s:+.3f}")
