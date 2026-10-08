# Analysis of the V43 null-source audit (note ray-compiler, section 4).
#   python nullinj.py NET DIR      (DIR holds ni_off<NET>_L<layer>_<gtype>_P<parts>_<+|->.npy from oracle_one.py SAVE_OUT)
# Response of the output to the injection, by central difference: r = (m(+eps) - m(-eps)) / 2. The exact response of the
# full tuple (PARTS = 15) is zero; the chain's is its gauge anomaly. Reported against the response to the covariance part
# alone (PARTS = 1), which is what the null tuple's other parts must cancel, and against the baseline error |e|.
import sys, glob, os, re, numpy as np
net, d = int(sys.argv[1]), sys.argv[2]
mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
base = np.load(f"/data/out/base32/free_{net}.npy").astype(np.float64)
e = base[-1] - mt[-1]
R = {}
for f in glob.glob(f"{d}/ni_off{net}_L*_P*_+.npy"):
    m = re.search(r"_L(\d+)_(\w+?)_P(\d+)_\+\.npy$", f)
    L_, g, P = int(m.group(1)), m.group(2), int(m.group(3))
    fm = f.replace("_+.npy", "_-.npy")
    if os.path.exists(fm):
        R[(L_, g, P)] = (np.load(f).astype(np.float64) - np.load(fm).astype(np.float64)) / 2
print(f"net {net}: |e| (baseline output error) {np.linalg.norm(e):.3e}")
print("  layer gtype : |r_cov| (P=1)  |r_tuple-no-feed|/|r_cov| (P=7)  |r_full|/|r_cov| (P=15)  |r_full|/|e|   cos(r_full, e)")
for L_, g in sorted({(k[0], k[1]) for k in R}):
    rc, r7, r15 = (R.get((L_, g, P)) for P in (1, 7, 15))
    if rc is None:
        continue
    c = np.linalg.norm(rc[-1])
    f7 = np.linalg.norm(r7[-1]) / c if r7 is not None else float("nan")
    f15 = np.linalg.norm(r15[-1]) / c if r15 is not None else float("nan")
    fe = np.linalg.norm(r15[-1]) / np.linalg.norm(e) if r15 is not None else float("nan")
    ce = float(r15[-1] @ e / (np.linalg.norm(r15[-1]) * np.linalg.norm(e))) if r15 is not None else float("nan")
    print(f"   {L_:2d}  {g:4s} : {c:.3e}       {f7:.3f}                          {f15:.3f}                {fe:.3e}   {ce:+.3f}")
    # per-layer profile of the full-tuple response relative to the covariance-only response
    if r15 is not None:
        prof = [np.linalg.norm(r15[k]) / max(np.linalg.norm(rc[k]), 1e-30) for k in range(L_, r15.shape[0])]
        print("        |r_full|/|r_cov| by layer from the injection: " + " ".join(f"{x:.2f}" for x in prof))
