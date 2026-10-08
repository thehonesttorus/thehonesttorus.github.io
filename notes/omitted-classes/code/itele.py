# Intervention telescoping: where the chain's output error is created (note XXXVI section 3k).
#   python itele.py NET FREE.npy PREFIX_A PREFIX_B [GROUPPREFIX ...]
# FREE.npy: the free-running chain's saved means (oracle_one.py SAVE_OUT). PREFIX_A/B + "_k.npy", k = 0..15: runs with every
# readout oracle (D3 D21 G4 WK4M K31 VAR COFF MU) at layers 0..k, from the two independent Monte Carlo halves A and B.
# Telescoping (exact):  |e_free|^2 - |e_15|^2 = sum_k (|e_(k-1)|^2 - |e_k|^2),  e_(-1) = e_free,
# and the k-th term is the output effect of the error the chain creates at step k, carried forward by the chain's own
# dynamics. The oracle values carry Monte Carlo noise; with independent halves the cross product <e_k^A, e_k^B> is the
# noise-free error energy (the noise enters linearly and independently), so the noise-free created-error effect at step k
# is <e_(k-1)^A, e_(k-1)^B> - <e_k^A, e_k^B>; |A - B| / 2 is the noise bar.
# Optional GROUPPREFIX_A/B_k.npy: runs that also replace, at layer k only, a single group (e.g. kappa3 readouts only) on top
# of the all-oracle state at layers < k; their difference from (k-1) splits the created error by stage.
import sys, numpy as np
net, freep, pa, pb = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
groups = sys.argv[5:]
mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
L, n = mt.shape
out = lambda p: np.load(p).astype(np.float64)[-1] - mt[-1]
ef = out(freep); E0 = float(ef @ ef)
A = {-1: ef}; B = {-1: ef}
for k in range(L):
    try:
        A[k] = out(f"{pa}_{k}.npy"); B[k] = out(f"{pb}_{k}.npy")
    except FileNotFoundError:
        break
K = max(A)
cross = {k: float(A[k] @ B[k]) for k in A}
print(f"net {net}: free-running n MSE {E0:.4e}; all-oracle through layer {K}: noise-free {100 * cross[K] / E0:.1f}% of free "
      f"(half A {100 * float(A[K] @ A[K]) / E0:.1f}%, half B {100 * float(B[K] @ B[K]) / E0:.1f}%)")
print("  step k: created-error effect = <e_(k-1)^A, e_(k-1)^B> - <e_k^A, e_k^B>, % of the free n MSE (+- half difference); cumulative")
cum = 0.0
for k in range(K + 1):
    d = cross[k - 1] - cross[k]
    dA = float(A[k - 1] @ A[k - 1] - A[k] @ A[k]); dB = float(B[k - 1] @ B[k - 1] - B[k] @ B[k])
    cum += d
    print(f"    {k:2d}: {100 * d / E0:+7.2f} (+- {100 * abs(dA - dB) / 2 / E0:.2f})   cum {100 * cum / E0:+7.2f}")
for gp in groups:
    print(f"  group runs {gp}: step k effect of replacing only this group at layer k (on top of all oracles at < k)")
    for k in range(K + 1):
        try:
            gA, gB = out(f"{gp}_A_{k}.npy"), out(f"{gp}_B_{k}.npy")
        except FileNotFoundError:
            continue
        d = cross[k - 1] - float(gA @ gB)
        print(f"    {k:2d}: {100 * d / E0:+7.2f} of the full step {100 * (cross[k - 1] - cross[k]) / E0:+7.2f}")
