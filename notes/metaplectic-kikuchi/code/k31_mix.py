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
for s_ in range(3,15):
    t=s_+1
    F=T["full"]; mu=F["mu"][s_].astype(np.float64); C=F["cov"][s_].astype(np.float64); C=.5*(C+C.T); v=np.diag(C).copy()
    Coff=C.copy(); np.fill_diagonal(Coff,0); ref=F["var"][t].astype(np.float64)
    k3c=c2[f"D3_{s_}"].astype(np.float64); Dc=c2[f"D21_{s_}"].astype(np.float64).copy(); np.fill_diagonal(Dc,0)
    vc=c2[f"var_{s_}"].astype(np.float64)
    Cc=c2[f"C_off_{s_}"].astype(np.float64); Cc=.5*(Cc+Cc.T); np.fill_diagonal(Cc,0)
    sc=2*k3c[:,None]*Dc/vc[:,None]-(2/3)*(k3c**2)[:,None]*Cc/(vc**2)[:,None]   # entirely from the chain's own state
    Kc=c2[f"wk431_{s_}"].astype(np.float64).T.copy(); np.fill_diagonal(Kc,0)
    q=lambda X: qdiag(W[t],t31(mu,C,X))/ref
    yA=q(T["h0"]["K31"][s_].astype(np.float64)*msk); yB=q(T["h1"]["K31"][s_].astype(np.float64)*msk); ey=np.mean(yA*yB)
    X=np.stack([q(Kc),q(sc)],1); G=X.T@X; b=np.linalg.solve(G,0.5*X.T@(yA+yB)); r=1-np.mean((yA-X@b)*(yB-X@b))/ey
    def red(w): x=X@np.array(w); return 1-np.mean((yA-x)*(yB-x))/ey
    print(f"layer {t}: best (chain, 2nd-chaos) = ({b[0]:.2f}, {b[1]:.2f}) removes {r:+.2f} | fixed (1,0) {red([1,0]):+.2f} (1,0.5) {red([1,.5]):+.2f} (1,1) {red([1,1]):+.2f} (1.3,0.7) {red([1.3,.7]):+.2f}",flush=True)
