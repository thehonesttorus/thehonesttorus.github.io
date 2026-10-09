# Can the chain's own D3 and D21 repair its (3,1) slice? Variance-metric fits of the K31-term error only.
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
rows=[]
for s_ in range(3,15):
    t=s_+1
    F=T["full"]; mu=F["mu"][s_].astype(np.float64); C=F["cov"][s_].astype(np.float64); C=.5*(C+C.T); v=np.diag(C).copy()
    Coff=C.copy(); np.fill_diagonal(Coff,0); ref=F["var"][t].astype(np.float64)
    k3c=c2[f"D3_{s_}"].astype(np.float64); Dc=c2[f"D21_{s_}"].astype(np.float64).copy(); np.fill_diagonal(Dc,0)
    Kc=c2[f"wk431_{s_}"].astype(np.float64).T.copy(); np.fill_diagonal(Kc,0)
    q=lambda X: qdiag(W[t],t31(mu,C,X))/ref
    fc=q(Kc); fe=q(k3c[:,None]*Dc/v[:,None]); fs=q(3*v[:,None]*Coff)
    y={}
    for h in ("h0","h1"):
        Kt=T[h]["K31"][s_].astype(np.float64).copy(); np.fill_diagonal(Kt,0); y[h]=q(Kt)-fc
    yA,yB=y["h0"],y["h1"]; ey=np.mean(yA*yB)
    out=[]
    for name,X in (("amp",np.stack([fc],1)),("exch",np.stack([fe],1)),("amp+exch",np.stack([fc,fe],1)),("amp+exch+svC",np.stack([fc,fe,fs],1))):
        G=X.T@X; b=np.linalg.solve(G,0.5*X.T@(yA+yB)); rA,rB=yA-X@b,yB-X@b; out.append((name,b,1-np.mean(rA*rB)/ey))
    print(f"layer {t}: K31-term error nf {np.sqrt(ey):.2e} | "+" | ".join(f"{n} share {sh:+.2f} coef "+",".join(f"{x:+.2f}" for x in b) for n,b,sh in out),flush=True)
