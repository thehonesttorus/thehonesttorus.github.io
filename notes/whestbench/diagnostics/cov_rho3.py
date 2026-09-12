# Choose the coefficient on held-out generalisation, not on the fit set, and compare the two
# parametrisations (per-layer g4 schedule against a flat constant).
import numpy as np
import cov_rho as CR
FIT=[0,2,4,6]; HOLD=list(range(8,16))
def sc(nets,**kw): return float(np.mean([CR.predict(k,**kw) for k in nets]))
b=sc(FIT); bh=sc(HOLD)
print(f"shipped (A4 only):  fit {b:.4e}   held {bh:.4e}\n")
print(f"{'spec':22s}  {'fit':>11} {'fit ratio':>10} {'held':>11} {'held ratio':>11}")
for v in (0.12,0.17,0.22,0.27,0.33):
    f=sc(FIT,S3=v); h=sc(HOLD,S3=v)
    print(f"{'S3 (x g4_l) = %.2f'%v:22s}  {f:11.4e} {f/b:10.4f} {h:11.4e} {h/bh:11.4f}", flush=True)
for v in (0.005,0.007,0.010,0.013):
    f=sc(FIT,S4=v); h=sc(HOLD,S4=v)
    print(f"{'S4 (constant) = %.3f'%v:22s}  {f:11.4e} {f/b:10.4f} {h:11.4e} {h/bh:11.4f}", flush=True)
