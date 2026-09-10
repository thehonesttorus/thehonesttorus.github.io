# Item 10: does the fitted low-parameter correction act as a self-energy for the omitted sector?
# At each layer compare  tail = (exact sources, depth 8) - (exact sources, depth 3)  against
# F = the 2-parameter aligned correction fitted to the Monte-Carlo residual of the depth-3 model.
import numpy as np
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=T.load_nets()
def analyse(k, N_MC=2400, Kret=3, Kfull=8, seed=0):
    W=T.regen(seeds[k]); Y=Yall[k]
    rng=np.random.default_rng(seed); X=rng.standard_normal((N_MC,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    HS=np.load('../pass2_stats.npz')['HS'] if k==0 else None
    sources={}; rows=[]
    # oracle-input single-layer study: use true means and sampled covariances so the comparison is clean
    Sig=None
    for l in range(L):
        mu=(Y[l-1].astype(np.float64)@W[l].astype(np.float64)) if l>0 else np.zeros(n)
        if HS is not None: Sz=W[l].T.astype(np.float64)@HS[l]@W[l].astype(np.float64)
        else:
            zl=Zs[l].astype(np.float64); Sz=np.cov(zl.T,bias=True)
        sg=np.sqrt(np.diag(Sz)); t=mu/sg; ph=pdf(t); Ph=ndtr(t)
        m1=mu*Ph+sg*ph; m2=(mu*mu+sg*sg)*Ph+mu*sg*ph; m3=mu*m2+2*sg*sg*m1
        M=Ph[:,None]*Sz.copy(); np.fill_diagonal(M,0)
        sources[l]=(M, ph/sg, 2*m1*(1-Ph), m3-3*m2*m1+2*m1**3, Ph)
    for l in range(4,L):
        Wl=W[l].astype(np.float64)
        mu=Y[l-1].astype(np.float64)@Wl
        Sz=W[l].T.astype(np.float64)@HS[l]@Wl if HS is not None else np.cov(Zs[l].astype(np.float64).T,bias=True)
        sg=np.sqrt(np.diag(Sz)); t=mu/sg; s3=sg**3
        def k3_depth(K):
            out=np.zeros(n); U=Wl.copy(); s_lo=max(l-K,0)
            for src in range(l-1,s_lo-1,-1):
                Ms,ps,qs,k3hs,d1s=sources[src]; Nt=Ms.T@U; UN=U*Nt
                out+=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
                if src>s_lo: U=W[src]@(d1s[:,None]*U)
            return out
        A=k3_depth(Kret); B=k3_depth(Kfull); tail=B-A
        zl=Zs[l].astype(np.float64); a=zl.mean(0); zc=zl-a; k3s=(zc**3).mean(0)
        r=(k3s-A)/s3; Xb=np.column_stack([t,np.ones(n)]); cf=np.linalg.lstsq(Xb,r,rcond=None)[0]; F=(Xb@cf)*s3
        c=np.corrcoef(F,tail)[0,1]
        # how much of the tail does the 2-parameter aligned model explain, at oracle coefficients?
        cft=np.linalg.lstsq(Xb,tail/s3,rcond=None)[0]; Fo=(Xb@cft)*s3
        r2=1-np.mean((tail-Fo)**2)/np.mean(tail**2)
        rows.append((l+1,c,np.linalg.norm(F)/np.linalg.norm(tail),r2,np.linalg.norm(tail)/np.linalg.norm(B)))
    return rows
print("Does the fitted 2-parameter correction stand in for the omitted source layers?")
print("corr(F,tail) = correlation of the MC-fitted correction with the exact depth-3->8 tail")
print("|F|/|tail|   = magnitude ratio")
print("R2 aligned   = share of the tail an oracle 2-parameter aligned model could explain")
print("|tail|/|k3|  = size of the omitted sector relative to the full depth-8 cumulant\n")
print("layer | corr(F,tail) | |F|/|tail| | R2 aligned | |tail|/|k3|")
for row in analyse(0):
    print(f"{row[0]:5d} |   {row[1]:+.3f}     |   {row[2]:.2f}     |   {row[3]:.3f}    |   {row[4]:.3f}")
