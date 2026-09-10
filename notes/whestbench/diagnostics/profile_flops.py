import numpy as np, flopscope as flops, flopscope.numpy as fnp, warnings, time, sys
warnings.simplefilter('ignore')
import estimator as E, twin2 as T
from types import SimpleNamespace
seeds,Y=T.load_nets(); W=T.regen(seeds[0]); mlp=SimpleNamespace(width=1024,depth=16,weights=W,seed=123,name='n')
B=2**41
# monkeypatch matmul/einsum to attribute cost: simpler -> measure total with pieces disabled
def run(**kw):
    est=E.FiniteResolutionCumulantEstimator()
    for k,v in kw.items(): setattr(est,k,v)
    with flops.BudgetContext(flop_budget=B) as ctx:
        t0=time.time(); P=est.predict(mlp,B); dt=time.time()-t0
    mse=np.mean((np.asarray(P).astype(np.float64)-Y[0])**2,axis=1)
    return ctx.flops_used/B, dt, mse[-1]
for kw in [dict(),dict(USE_MC=False),dict(USE_SOURCE=False),dict(USE_MC=False,USE_SOURCE=False),dict(MC_MAX=2000),dict(T_ACT=2.5)]:
    f,dt,m=run(**kw); print(f"{kw}: budget fraction {f:.4f}  time {dt:.1f}s  final MSE {m:.2e}", flush=True)
