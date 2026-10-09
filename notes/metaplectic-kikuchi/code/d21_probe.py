import sys, numpy as np
sys.argv=[sys.argv[0]]+[".",".","0","1"]
exec(open("var_ladder.py").read().split("cd = np.load")[0])
net=int(sys.argv[-1]) if False else 0
T={h:np.load(f"../data/mc2/mc2_off{net}_{h}.npz") for h in ("full","h0","h1")}
c2=np.load(f"chaindump2_{net}.npz"); W=np.load(f"../data/loc/W_off{net}.npy").astype(np.float64)
cd=np.load(f"../data/loc/chaindump_{net}.npz")
def t21(mu,C,D21):
    v=np.diag(C).copy(); s=np.sqrt(v); al=mu/s; R=C/np.outer(s,s); np.fill_diagonal(R,0)
    c=ccoef(al,M+5); X=0.5*D21/s[:,None]*series(c,c,R,2,1,0); X=X+X.T; np.fill_diagonal(X,0); return X
msk=~np.eye(1024,dtype=bool)
for s_ in (5,8,11,13,14):
    t=s_+1
    Dc=c2[f"D21_{s_}"].astype(np.float64).copy(); np.fill_diagonal(Dc,0)
    F=T["full"]; mu=F["mu"][s_].astype(np.float64); C=F["cov"][s_].astype(np.float64); C=.5*(C+C.T); v=np.diag(C).copy()
    k3=F["k3"][s_].astype(np.float64)
    Coff=C.copy(); np.fill_diagonal(Coff,0)
    Dt={h:T[h]["D21"][s_].astype(np.float64).copy() for h in ("full","h0","h1")}
    for h in Dt: np.fill_diagonal(Dt[h],0)
    dA,dB=Dt["h0"]-Dc,Dt["h1"]-Dc
    print(f"layer {s_}: |D21 true| {np.sqrt(np.mean(Dt['full'][msk]**2)):.3e} |chain| {np.sqrt(np.mean(Dc[msk]**2)):.3e} "
          f"entry corr(true,chain) {np.corrcoef(Dt['full'][msk],Dc[msk])[0,1]:.4f}; diff nf rms {np.sqrt(max(np.mean(dA[msk]*dB[msk]),0)):.3e} "
          f"(MC noise rms {np.sqrt(np.mean(((Dt['h0']-Dt['h1'])[msk]/2)**2)):.3e})")
    # row/column structure of the difference: row means and column means (rank-1 additive parts)
    dd=0.5*(dA+dB)
    rm=dd.mean(1); cm=dd.mean(0)
    print(f"   diff energy in row means {np.sum(rm**2)*1024/np.sum(dd[msk]**2):.3f}, column means {np.sum(cm**2)*1024/np.sum(dd[msk]**2):.3f}")
    cands={"true":Dt["full"],"k3C/v":k3[:,None]*Coff/v[:,None],"C":Coff,"Ct":Coff*0+ (k3[None,:]*Coff/v[None,:])}
    for k,X in cands.items():
        xa=X[msk]; r=np.mean(dA[msk]*xa)/np.sqrt(np.mean(dA[msk]*dB[msk])*np.mean(xa*xa)) if np.mean(dA[msk]*dB[msk])>0 else np.nan
        print(f"   entry corr(diff, {k}) (noise-free) {r:+.3f}")
    # variance-metric: which neurons: correlation of the D21-term error with target features
    ref=F["var"][t].astype(np.float64)
    eA=qdiag(W[t],t21(mu,C,Dt["h0"]))/ref-qdiag(W[t],t21(mu,C,Dc))/ref; eB=qdiag(W[t],t21(mu,C,Dt["h1"]))/ref-qdiag(W[t],t21(mu,C,Dc))/ref
    print(f"   var-metric D21-term error nf rms {np.sqrt(max(np.mean(eA*eB),0)):.2e}, mean {np.mean(0.5*(eA+eB)):+.2e}")
