import numpy as np, sys
from closure import closure
from edgeworth import cumulants_next, edgeworth_shift
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = np.load(f"W_n{n}_L{L}_s{s}.npy"); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); m = tr["m"]
out = closure(list(Ws), keep=True)
o1, o2 = out[0], out[1]
e0 = o2["m"]-m[1]
k3, k4, parts = cumulants_next(Ws[1], o1["mu"], o1["sig"], o1["R"])
for name, sel in [("k3 only", ("D3","P3","T3")), ("k4 only", ("D4","PP4","PD4")), ("all", ("D3","P3","T3","D4","PP4","PD4")),
                  ("D3+D4", ("D3","D4")), ("no T3", ("D3","P3","D4","PP4","PD4")), ("no PD4", ("D3","P3","T3","D4","PP4"))]:
    K3 = sum(parts[t] for t in sel if t.endswith("3")) if any(t.endswith("3") for t in sel) else 0*k3
    K4 = sum(parts[t] for t in sel if t.endswith("4")) if any(t.endswith("4") for t in sel) else 0*k4
    d = edgeworth_shift(o2["mu"], o2["sig"], K3, K4)
    e = e0+d
    print(f"{name:8s}: MSE {np.mean(e**2):.3e} (closure {np.mean(e0**2):.3e}, noise {tr['v'][1].mean()/tr['T']:.1e}) mean err {e.mean():+.2e}")
for t in parts: print(t, "rms %.3e" % np.sqrt(np.mean(parts[t]**2)), "mean %+.3e" % parts[t].mean())
