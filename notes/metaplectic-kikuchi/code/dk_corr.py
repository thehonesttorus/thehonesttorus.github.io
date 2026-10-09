# Are the chain's D21-term and K31-term errors (variance metric) anti-correlated?
import sys, numpy as np
NET_ARG=int(sys.argv[1])
sys.argv=[sys.argv[0]]+[".",".","0","1"]
exec(open("var_ladder.py").read().split("cd = np.load")[0])
net=NET_ARG
T={h:np.load(f"../data/mc2/mc2_off{net}_{h}.npz") for h in ("full","h0","h1")}
c2=np.load(f"chaindump2_{net}.npz"); W=np.load(f"../data/loc/W_off{net}.npy").astype(np.float64)
def tt(mu,C,X,p,q,div):
    v=np.diag(C).copy(); s=np.sqrt(v); al=mu/s; R=C/np.outer(s,s); np.fill_diagonal(R,0)
    c=ccoef(al,M+5); Y=div(s)*X*series(c,c,R,p,q,0); Y=Y+Y.T; np.fill_diagonal(Y,0); return Y
t21=lambda mu,C,X: tt(mu,C,0.5*X,2,1,lambda s: 1/s[:,None])
t31=lambda mu,C,X: tt(mu,C,X/6.0,3,1,lambda s: s[:,None]**-2)
msk=~np.eye(1024,dtype=bool)
for s_ in (5,8,11,13,14):
    t=s_+1
    F=T["full"]; mu=F["mu"][s_].astype(np.float64); C=F["cov"][s_].astype(np.float64); C=.5*(C+C.T); ref=F["var"][t].astype(np.float64)
    q=lambda X: qdiag(W[t],X)/ref
    Dc=c2[f"D21_{s_}"].astype(np.float64)*msk; Kc=c2[f"wk431_{s_}"].astype(np.float64).T*msk
    d={};k={}
    for h in ("h0","h1"):
        d[h]=q(t21(mu,C,Dc))-q(t21(mu,C,T[h]["D21"][s_].astype(np.float64)*msk))
        k[h]=q(t31(mu,C,Kc))-q(t31(mu,C,T[h]["K31"][s_].astype(np.float64)*msk))
    dd=np.mean(d["h0"]*d["h1"]); kk=np.mean(k["h0"]*k["h1"]); dk=0.5*np.mean(d["h0"]*k["h1"]+d["h1"]*k["h0"])
    print(f"layer {t}: |d| {np.sqrt(dd):.2e} |k| {np.sqrt(kk):.2e} corr(d,k) {dk/np.sqrt(dd*kk):+.3f} | |d+k|^2/(|d|^2+|k|^2) {(dd+kk+2*dk)/(dd+kk):.3f}",flush=True)
