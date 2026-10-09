# Parameter-free second-chaos closure of the (3,1) slice: K31[a,b] = 2 k3(a) D21(a,b)/v_a - (2/3) k3(a)^2 C_ab / v_a^2
# (single-factor quadratic latent model, leading order). Against Monte Carlo truth, from true inputs and from the chain's.
import sys, numpy as np
NET_ARG=int(sys.argv[1])
sys.argv=[sys.argv[0]]+[".",".","0","1"]
exec(open("var_ladder.py").read().split("cd = np.load")[0])
net=NET_ARG  # NET_FIX: var_ladder's header overwrites net
T={h:np.load(f"../data/mc2/mc2_off{net}_{h}.npz") for h in ("full","h0","h1")}
c2=np.load(f"chaindump2_{net}.npz"); W=np.load(f"../data/loc/W_off{net}.npy").astype(np.float64)
def t31(mu,C,X):
    v=np.diag(C).copy(); s=np.sqrt(v); al=mu/s; R=C/np.outer(s,s); np.fill_diagonal(R,0)
    c=ccoef(al,M+5); Y=(X/6.0)*s[:,None]**-2*series(c,c,R,3,1,0); Y=Y+Y.T; np.fill_diagonal(Y,0); return Y
msk=~np.eye(1024,dtype=bool)
for s_ in (5,8,11,13,14):
    t=s_+1
    F=T["full"]; mu=F["mu"][s_].astype(np.float64); C=F["cov"][s_].astype(np.float64); C=.5*(C+C.T); v=np.diag(C).copy()
    Coff=C.copy(); np.fill_diagonal(Coff,0); ref=F["var"][t].astype(np.float64)
    sc=lambda k3,D21: 2*k3[:,None]*D21/v[:,None]-(2/3)*(k3**2)[:,None]*Coff/(v**2)[:,None]
    Dt=F["D21"][s_].astype(np.float64).copy(); np.fill_diagonal(Dt,0)
    Kt=F["K31"][s_].astype(np.float64).copy(); np.fill_diagonal(Kt,0)
    Ktrue_sc=sc(F["k3"][s_].astype(np.float64),Dt)
    k3c=c2[f"D3_{s_}"].astype(np.float64); Dc=c2[f"D21_{s_}"].astype(np.float64).copy(); np.fill_diagonal(Dc,0)
    Kchain_sc=sc(k3c,Dc)
    Kc=c2[f"wk431_{s_}"].astype(np.float64).T.copy(); np.fill_diagonal(Kc,0)
    def ent(X): 
        a=np.sum(Kt[msk]*X[msk])/np.sum(X[msk]**2); return a,np.corrcoef(Kt[msk],X[msk])[0,1]
    q=lambda X: qdiag(W[t],t31(mu,C,X))/ref
    yA=q(T["h0"]["K31"][s_].astype(np.float64)*msk); yB=q(T["h1"]["K31"][s_].astype(np.float64)*msk); ey=np.mean(yA*yB)
    def vm(X):
        x=q(X); amp=np.mean(x*(yA+yB))/(2*np.mean(x*x)); 
        # residual energy if X replaces the truth with unit amplitude
        r1=np.mean((yA-x)*(yB-x))/ey; return amp, 1-r1
    out=[]
    for name,X in (("chain closure",Kc),("2nd-chaos (true inputs)",Ktrue_sc),("2nd-chaos (chain inputs)",Kchain_sc),("chain + 2nd-chaos/2",Kc+0.5*Kchain_sc)):
        a,r=ent(X); amp,red=vm(X)
        out.append(f"{name}: entry amp {a:.2f} corr {r:.3f}; var-metric amp {amp:.2f}, energy removed at unit amplitude {red:+.2f}")
    print(f"layer {s_}->{t} (true K31 term nf {np.sqrt(ey):.2e}):\n   "+"\n   ".join(out),flush=True)
