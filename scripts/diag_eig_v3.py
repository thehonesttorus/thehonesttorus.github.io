"""Chain vs Monte Carlo cache: eigenvalues of the pre-activation covariance (top modes) and the chain's error in the
top-k eigen-subspace of the MC covariance vs the bulk.   python scripts/diag_eig_mc.py DATA NET [opt=val]"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.k3chain3 import k3_chain3 as k3_chain2
D, net = sys.argv[1], int(sys.argv[2]); opts = {}
for a in sys.argv[3:]:
    k, v = a.split("="); opts[k] = v if k in ("k4", "k22gate") else (float(v) if "." in v else int(v))
W = np.load(f"{D}/W_off{net}.npy"); mc = np.load(f"{D}/mccache_{net}.npz"); S = float(mc["S"])
rec = {}; out, fl = k3_chain2(W, opts, record=rec); L, n, _ = W.shape
z1 = mc["s1"] / S
print("layer | top-4 eigenvalues: MC | chain/MC - 1 | rel err of C in top-4 subspace (||Q^T dC Q||_F/||Q^T C Q||_F) | in bulk (complement, off-diag) | trace rel err | mu-direction variance rel err")
for l in range(L):
    Cm = mc["Hc"][l].astype(np.float64) - np.outer(z1[l], z1[l]); C = rec[l]["C"]
    w, Q = np.linalg.eigh(Cm); idx = np.argsort(w)[::-1][:4]; w = w[idx]; Q = Q[:, idx]
    wc = np.einsum("ik,ij,jk->k", Q, C, Q)
    dC = C - Cm; top = np.linalg.norm(Q.T @ dC @ Q) / np.linalg.norm(Q.T @ Cm @ Q)
    P = np.eye(n) - Q @ Q.T; dB = P @ dC @ P; CB = P @ Cm @ P
    offb = np.linalg.norm(dB - np.diag(np.diag(dB))) / np.linalg.norm(CB - np.diag(np.diag(CB)))
    u = z1[l] / np.linalg.norm(z1[l]); vm = u @ Cm @ u; vc = u @ C @ u
    print(f"{l:5d} | " + " ".join(f"{x:7.3f}" for x in w) + " | " + " ".join(f"{x:+7.4f}" for x in wc / w - 1) + f" | {top:8.4f} | {offb:8.4f} | {np.trace(C)/np.trace(Cm)-1:+8.4f} | {vc/vm-1:+8.4f}")
