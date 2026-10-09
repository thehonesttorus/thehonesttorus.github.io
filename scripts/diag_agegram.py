"""Per-age structure of the third-cumulant memory: at chosen layers, the (2,1)-slice contribution of each source age,
its norm, the cosine matrix across ages, and the norm of the cumulative sum over ages (cancellation between ages).
  python scripts/diag_agegram.py DATA NET [layers=6,10,15]"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest import kprop3c
from whest.kprop3c import kprop3c_chain, slices_from_legs
D, net = sys.argv[1], int(sys.argv[2]); layers = [int(x) for x in (sys.argv[3].split("=")[1] if len(sys.argv) > 3 else "6,10,15").split(",")]
W = np.load(f"{D}/W_off{net}.npy").astype(np.float64)
per_age = {}
orig = kprop3c.nonlin_step
def hooked(st, o, rng, record=None):
    l = len(per_age); per_age[l] = {e[1]: slices_from_legs(e[0]) for e in st.young} if l in layers else None
    return orig(st, o, rng, record)
kprop3c.nonlin_step = hooked
kprop3c_chain(W, dict(window=99, k=8))
np.set_printoptions(precision=2, suppress=True, linewidth=200)
for l in layers:
    d = per_age[l]; ages = sorted(d); S = [d[a][0] for a in ages]; S3 = [d[a][1] for a in ages]
    tot = sum(S); tot3 = sum(S3)
    nrm = np.array([np.linalg.norm(x) for x in S]); nrm3 = np.array([np.linalg.norm(x) for x in S3])
    print(f"layer {l}: ages {ages}")
    print("  |s_a| / |total| (2,1):", nrm / np.linalg.norm(tot))
    print("  |s_a| / |total| (3,)  :", nrm3 / np.linalg.norm(tot3))
    cum = [np.linalg.norm(sum(S[:k + 1])) / np.linalg.norm(tot) for k in range(len(S))]
    print("  |sum_{a<=A} s_a| / |total| (2,1):", np.array(cum))
    cum3 = [np.linalg.norm(sum(S3[:k + 1])) / np.linalg.norm(tot3) for k in range(len(S3))]
    print("  |sum_{a<=A} s_a| / |total| (3,)  :", np.array(cum3))
    G = np.array([[np.sum(x * y) / (np.linalg.norm(x) * np.linalg.norm(y)) for y in S] for x in S])
    print("  cosine matrix across ages (2,1):"); print(G[:8, :8])
    cosT = np.array([np.sum(x * tot) / (np.linalg.norm(x) * np.linalg.norm(tot)) for x in S])
    print("  cosine of each age with the total:", cosT)
