# Two-channel transport of the exact local defects: critical total power q and decaying mean fraction c.
#   output scale contribution of layer l = 1/2 dq_l/q_l + 1/2 (dc_l / c_L) prod_{i>l} kappa'(c_i)
import numpy as np, sys
from ledger import gstep, run_from
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); hm = np.load(f"hmom_n{n}_L{L}_s{s}.npz")
mstar = tr["m"]; Cstar = [hm["S2"][l] - np.outer(hm["S1"][l], hm["S1"][l]) for l in range(L)]
sstar = [(np.zeros(n), Ws[0] @ Ws[0].T)] + [(Ws[l] @ mstar[l-1], Ws[l] @ Cstar[l-1] @ Ws[l].T) for l in range(1, L)]
mL = mstar[-1]; scL = lambda v: (mL @ v)/(mL @ mL)
cs, rows = [], []
for l in range(L):
    mloc, Cloc = gstep(*sstar[l]); ms = mstar[l]
    A = ms @ ms; Q = A + np.trace(Cstar[l]); c = A/Q; cs.append(c)
    dA = mloc @ mloc - A; dQ = (mloc @ mloc + np.trace(Cloc)) - Q
    dc = ((1-c)*dA - c*(dQ - dA))/Q
    mmC = ms @ (Cloc - Cstar[l]) @ ms/(A*Q)          # covariance defect along the mean direction (relative)
    rows.append((dQ/Q, dc, mmC))
cs = np.array(cs); kp = 1 - np.arccos(np.clip(cs, -1, 1))/np.pi
R = [run_from(Ws, l, *sstar[l]) for l in range(L)] + [mstar[-1]]
print(f"c_L = {cs[-1]:.3f}")
print("layer   dq/q       dc        m-aligned cov defect | predicted out (q-chan + c-chan) | ledger carried")
tot = 0
for l in range(L):
    dq, dc, mm = rows[l]
    pred = 0.5*dq + 0.5*dc/cs[-1]*np.prod(kp[l+1:])
    tot += pred
    print(f"{l+1:4d}  {dq:+.5f}  {dc:+.5f}   {mm:+.5f}               | {pred:+.5f}                        | {scL(R[l]-R[l+1]):+.5f}")
print(f"sum predicted {tot:+.5f}   vs ledger total {scL(R[0]-mL):+.5f}")
