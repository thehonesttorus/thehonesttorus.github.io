# How much does the final-layer mean amplify a relative error in the propagated sigma?  Perturb
# Sz <- Sz(1+v) at one layer at a time, and uniformly, and read the excess MSE.  The direct one-layer
# effect is phi(t) sigma epsilon with epsilon = v/2, so the ratio measures the amplification the
# remaining layers apply to it.
import numpy as np
import varcal as V
n=1024; L=16
NETS=[0,1,2,3]
base=np.array([V.predict(k,np.zeros(15)) for k in NETS])
print(f"deployed mean MSE {base.mean():.4e}\n")
EPS=7e-4          # relative sigma perturbation
v=2*EPS
print(f"relative sigma perturbation epsilon = {EPS:.1e}\n")
print(f"{'perturbed at':>14} {'MSE':>11} {'excess':>11} {'rms excess':>11} {'amplification':>14}")
direct=None
for l0 in list(range(1,16))+['all']:
    vv=np.zeros(15)
    if l0=='all': vv[:]=v
    else: vv[l0-1]=v
    m=np.array([V.predict(k,vv) for k in NETS])
    ex=max(m.mean()-base.mean(),0.0); rms=np.sqrt(ex)
    # direct single-layer effect: rms over neurons of phi(t) sigma epsilon, taken from layer 15 scale
    if direct is None: direct=rms   # layer 1 as the reference for "one layer, fully propagated"
    print(f"{str(l0):>14} {m.mean():11.4e} {ex:11.4e} {rms:11.4e} {rms/EPS:14.2f}", flush=True)
