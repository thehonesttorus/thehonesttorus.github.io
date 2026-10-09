import sys, numpy as np
sys.argv=[sys.argv[0]]+[".",".","0","1"]
exec(open("var_ladder.py").read().split("cd = np.load")[0])
net=0
T={h:np.load(f"../data/mc2/mc2_off{net}_{h}.npz") for h in ("full","h0","h1")}
c2=np.load(f"chaindump2_{net}.npz"); W=np.load(f"../data/loc/W_off{net}.npy").astype(np.float64)
def t31(mu,C,K31):
    v=np.diag(C).copy(); s=np.sqrt(v); al=mu/s; R=C/np.outer(s,s); np.fill_diagonal(R,0)
    c=ccoef(al,M+5); X=(K31/6.0)*s[:,None]**-2*series(c,c,R,3,1,0); X=X+X.T; np.fill_diagonal(X,0); return X
for s_ in (5,8,11,13,14):
    t=s_+1
    out={}
    for h in ("full","h0","h1"):
        F=T[h]; mu=F["mu"][s_].astype(np.float64); C=F["cov"][s_].astype(np.float64); C=.5*(C+C.T); v=np.diag(C).copy()
        K31=F["K31"][s_].astype(np.float64).copy(); np.fill_diagonal(K31,0)
        Coff=C.copy(); np.fill_diagonal(Coff,0)
        ev,U=np.linalg.eigh(C); U=U[:,::-1]; ev=ev[::-1]
        cand={"true":K31,"svC":3*v[:,None]*Coff}
        for k in (1,4,16):
            Ck=(U[:,:k]*ev[:k])@U[:,:k].T; Ckoff=Ck.copy(); np.fill_diagonal(Ckoff,0)
            cand[f"coll{k}"]=3*np.diag(Ck)[:,None]*Ckoff
        D21=F["D21"][s_].astype(np.float64).copy(); np.fill_diagonal(D21,0)
        k3=F["k3"][s_].astype(np.float64)
        cand["k3D21/v"]=k3[:,None]*D21/v[:,None]
        out[h]={k:qdiag(W[t],t31(mu,C,X)) for k,X in cand.items()}
        out[h]["vt"]=F["var"][t].astype(np.float64)
        if h=="full":
            # entrywise fits
            msk=~np.eye(1024,dtype=bool)
            for k,X in cand.items():
                if k=="true": continue
                a=np.sum(K31[msk]*X[msk])/np.sum(X[msk]**2); r=np.corrcoef(K31[msk],X[msk])[0,1]
                print(f" layer {s_}->{t} entrywise K31 vs {k}: amp {a:.3f} corr {r:.3f}")
    ref=out["full"]["vt"]; A,B=out["h0"],out["h1"]
    yA,yB=A["true"]/ref,B["true"]/ref; ey=np.mean(yA*yB)
    for k in A:
        if k in ("true","vt"): continue
        xA,xB=A[k]/ref,B[k]/ref
        amp=np.mean(xA*yB+xB*yA)/np.mean(2*xA*xB); sh=np.mean(xA*yB)*np.mean(xB*yA)/(np.mean(xA*xB)*ey)
        print(f"   var-metric: true T31 on {k}: amp {amp:+.3f}, explained share {sh:.3f}")
    # the chain's own
    cK=c2[f"wk431_{s_}"].astype(np.float64).T.copy(); np.fill_diagonal(cK,0)
    F=T["full"]; mu=F["mu"][s_].astype(np.float64); C=F["cov"][s_].astype(np.float64); C=.5*(C+C.T)
    xc=qdiag(W[t],t31(mu,C,cK))/ref
    amp=np.mean(xc*yB+xc*yA)/np.mean(2*xc*xc); sh=np.mean(xc*yB)*np.mean(xc*yA)/(np.mean(xc*xc)*ey)
    print(f"   var-metric: true T31 on chain's own K31 (at true state): amp {amp:+.3f}, explained share {sh:.3f}; nf rms true {np.sqrt(ey):.2e}")
