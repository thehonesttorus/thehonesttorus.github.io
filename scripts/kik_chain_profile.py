import sys, numpy as np
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__file__), ".."))
from scripts.kik_chain import step, kappas
from scripts.kik_merge import pool_net
R, D, net = sys.argv[1], sys.argv[2], int(sys.argv[3])
a = pool_net(R, net); N = a["N"]; truth = np.load(f"{D}/truth_off{net}.npz")["m"]
W = [np.asarray(w, dtype=np.float64) for w in np.load(f"{D}/W_off{net}.npy")]
rms = lambda x: np.sqrt(np.mean(x ** 2))
for mode, fixv in (("h4", False), ("h4", True)):
    mu = np.zeros(1024); C = W[0] @ W[0].T; rows = []
    for l in range(16):
        k3, k4, tm, vm = kappas(a, l, N); ve = rms(np.diag(C) - vm) / np.mean(vm); te = rms(mu - tm)
        if fixv:
            dd = np.sqrt(vm / np.diag(C)); C = C * np.outer(dd, dd)
        m, Kh = step(mu, C, k3, k4, mode=mode)
        rows.append(f"{l}:t{te:.1e} v{ve:.1e} m{rms(m - truth[l]):.1e}")
        if l < 15: mu = W[l + 1] @ m; C = W[l + 1] @ Kh @ W[l + 1].T
    print(("oracle-var " if fixv else "") + mode + ": " + " ".join(rows))
