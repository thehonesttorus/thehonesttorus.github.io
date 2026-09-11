import numpy as np
import cov4 as CV, ablate as AB
G_S=AB.G_S
FIT=[0,2,4,6]; HOLD=list(range(8,16))
def sc(nets,**kw): return float(np.mean([CV.predict(k,**kw) for k in nets]))
b=sc(FIT); bh=sc(HOLD)
print(f"deployed fit {b:.4e}  held {bh:.4e}\n")
print(f"{'spec':38s} {'fit':>11} {'held':>11} {'held ratio':>11}")
c=dict(A4=0.0020,form=0)
print(f"{'constant 0.0020':38s} {sc(FIT,**c):11.4e} {sc(HOLD,**c):11.4e} {sc(HOLD,**c)/bh:11.4f}",flush=True)
SCH=G_S/4.0/0.002      # normalised so gamma=1 means A4_l = g4_l/4
for g in (0.0008,0.0014,0.002,0.003):
    kw=dict(A4=g,form=0,sched=SCH)
    f=sc(FIT,**kw); h=sc(HOLD,**kw)
    print(f"{'A4_l = %.4f * (g4_l/4)/0.002'%g:38s} {f:11.4e} {h:11.4e} {h/bh:11.4f}",flush=True)
# a milder depth trend
ls=np.arange(1,16)/15.0
for q in (0.3,0.6):
    for g in (0.0016,0.0022,0.003):
        kw=dict(A4=g,form=0,sched=ls**q/np.mean(ls**q))
        f=sc(FIT,**kw); h=sc(HOLD,**kw)
        print(f"{'A4=%.4f * (l/15)^%.1f norm'%(g,q):38s} {f:11.4e} {h:11.4e} {h/bh:11.4f}",flush=True)
