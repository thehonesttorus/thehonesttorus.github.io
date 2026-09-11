# The theory-signed calibration hurts monotonically, so the closure is calibrated AROUND its own
# sigma deficit: some other error of the opposite sign is cancelling against it.  Fit the calibration
# end to end with a free sign and shape instead, which is the method that produced both shipped
# kernels.  v_l = A (l/15)^p, applied as Sz <- Sz (1 + v_l).
import numpy as np, sys, itertools
import varcal as V
FIT=[0,2,4,6]; HOLD=list(range(8,16))
def sc(v,nets): return float(np.mean([V.predict(k,v) for k in nets]))
L=15; ls=(np.arange(1,L+1)/L)
base_f=sc(np.zeros(L),FIT)
print(f"deployed  fit {base_f:.4e}\n")
print(f"{'A':>10} {'p':>4}  {'fit(0,2,4,6)':>13}  {'vs deployed':>12}")
best=(base_f,0.0,0.0)
for A in (-0.006,-0.003,-0.0015,-0.0007,0.0007,0.0015,0.003):
    for p in (0.0,1.0,2.0):
        v=A*ls**p; f=sc(v,FIT)
        if f<best[0]: best=(f,A,p)
        print(f"{A:10.4f} {p:4.1f}  {f:13.4e}  {f/base_f:12.4f}", flush=True)
f,A,p=best
print(f"\nbest: A={A}, p={p},  fit {f:.4e}")
if A!=0.0:
    v=A*ls**p
    print(f"held-out(8-15):  deployed {sc(np.zeros(L),HOLD):.4e}   calibrated {sc(v,HOLD):.4e}")
