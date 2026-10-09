# variance-metric value of cheap n^2 completions of the chain's D21 and (3,1) slices (fit per layer, cross halves)
import sys, numpy as np
net=int(sys.argv[1])
sys.argv=[sys.argv[0]]+[".",".","0","1"]
exec(open("var_ladder.py").read().split("cd = np.load")[0])
T={h:np.load(f"../data/mc2/mc2_off{net}_{h}.npz") for h in ("full","h0","h1")}
c2=np.load(f"chaindump2_{net}.npz"); W=np.load(f"../data/loc/W_off{net}.npy").astype(np.float64)
def tt(mu,C,X,p,q,div):
    v=np.diag(C).copy(); s=np.sqrt(v); al=mu/s; R=C/np.outer(s,s); np.fill_diagonal(R,0)
    c=ccoef(al,M+5); Y=div(s)*X*series(c,c,R,p,q,0); Y=Y+Y.T; np.fill_diagonal(Y,0); return Y
t21=lambda mu,C,X: tt(mu,C,0.5*X,2,1,lambda s: 1/s[:,None])
t31=lambda mu,C,X: tt(mu,C,X/6.0,3,1,lambda s: s[:,None]**-2)
res=[]
for s_ in range(3,15):
    t=s_+1
    F=T["full"]; mu=F["mu"][s_].astype(np.float64); C=F["cov"][s_].astype(np.float64); C=.5*(C+C.T); v=np.diag(C).copy()
    Coff=C.copy(); np.fill_diagonal(Coff,0)
    k3c=c2[f"D3_{s_}"].astype(np.float64)                     # the chain's own kappa3 diagonal
    Dc=c2[f"D21_{s_}"].astype(np.float64).copy(); np.fill_diagonal(Dc,0)
    Kc=c2[f"wk431_{s_}"].astype(np.float64).T.copy(); np.fill_diagonal(Kc,0)
    ref=F["var"][t].astype(np.float64)
    q=lambda X: qdiag(W[t],X)/ref
    # chain-computable candidate directions (all from the chain's own D3, D21 and the covariance)
    reg=k3c[:,None]*Coff/v[:,None]
    colm=np.ones((1024,1))@Dc.mean(0,keepdims=True); np.fill_diagonal(colm,0)
    feats={"D21c":q(t21(mu,C,Dc)),"reg":q(t21(mu,C,reg)),"C":q(t21(mu,C,Coff*np.sqrt(v)[:,None])),
           "K31c":q(t31(mu,C,Kc)),"k3D21/v":q(t31(mu,C,k3c[:,None]*Dc/v[:,None])),"svC":q(t31(mu,C,3*v[:,None]*Coff))}
    # targets: true D21 term minus chain D21 term, true K31 term minus chain K31 term (halves)
    y={}
    for h in ("h0","h1"):
        G=T[h]; Dt=G["D21"][s_].astype(np.float64).copy(); np.fill_diagonal(Dt,0)
        Kt=G["K31"][s_].astype(np.float64).copy(); np.fill_diagonal(Kt,0)
        y[h]=q(t21(mu,C,Dt))-feats["D21c"] + q(t31(mu,C,Kt))-feats["K31c"]
    yA,yB=y["h0"],y["h1"]; ey=np.mean(yA*yB)
    def fit(keys):
        X=np.stack([feats[k] for k in keys],1); G=X.T@X; b=np.linalg.solve(G,0.5*X.T@(yA+yB))
        rA,rB=yA-X@b,yB-X@b; return b,1-np.mean(rA*rB)/ey
    out=[t,np.sqrt(ey)]
    for keys in (["D21c"],["D21c","reg"],["D21c","reg","C"],["K31c","k3D21/v"],["D21c","reg","C","K31c","k3D21/v","svC"]):
        b,sh=fit(keys); out.append((keys,b,sh))
    res.append(out)
    print(f"layer {t}: D21+K31 term error nf {np.sqrt(ey):.2e}")
    for keys,b,sh in out[2:]:
        print(f"   {'+'.join(keys):40s} share {sh:+.3f}  coef " + " ".join(f"{x:+.3f}" for x in b))
