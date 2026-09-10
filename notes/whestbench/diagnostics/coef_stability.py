# The reduced memory kernel has rank 2 in the basis {sigma^3, t sigma^3}.  Are its two coefficients
# stable across networks?  If so they can be shipped as constants and the Monte-Carlo pass dropped.
import numpy as np
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=T.load_nets()
def coefs(k, N_MC=2400, K=1, seed=0):
    W=T.regen(seeds[k])
    rng=np.random.default_rng(seed); X=rng.standard_normal((N_MC,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    out=[]; sources={}
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); M=G*0.5; M[np.diag_indices(n)]=0
    sources[0]=(M,INV/s,m,s2*s*(2*INV-1.5*INV+2*INV**3),np.full(n,0.5)); Sig=C
    for l in range(1,L):
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sig)@Wl; sz2=np.diag(Sz).copy(); sz=np.sqrt(sz2); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        k3=np.zeros(n); U=Wl.astype(np.float64); s_lo=max(l-K,0)
        for src in range(l-1,s_lo-1,-1):
            Ms,ps,qs,k3hs,d1s=sources[src]; Nt=Ms.T@U; UN=U*Nt
            k3+=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
            if src>s_lo: U=W[src]@(d1s[:,None]*U)
        z=Zs[l].astype(np.float64); a=z.mean(0); zc=z-a; v=(zc*zc).mean(0)
        k3s=(zc**3).mean(0); k4s=(zc**4).mean(0)-3*v*v
        r=(k3s-k3)/s3; Xb=np.column_stack([t,np.ones(n)]); cf=np.linalg.lstsq(Xb,r,rcond=None)[0]
        g4=np.mean(k4s/s4)
        out.append((cf[0],cf[1],g4))
        k3=k3+(Xb@cf)*s3; resid=k3s-k3
        noise=6*np.mean(s3*s3)/N_MC; sig=max(np.mean(resid**2)-noise,0); k3=k3+(sig/(sig+noise))*resid
        k4=g4*s4
        t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)+k4*he2*ph/(24*s3)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*s4)-k4*t*ph/(12*sz2)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rho=Sz/so; ds=T.hermite_d(t,ph,Ph,8); rk=rho.copy(); C=np.zeros((n,n)); fact=1.0
        for kk in range(1,9):
            fact*=kk; C+=rk*np.outer(ds[kk-1],ds[kk-1])/fact
            if kk<8: rk*=rho
        C*=so; C[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
        Msrc=Sz*Ph[:,None]; Msrc[np.diag_indices(n)]=0
        m3=mu*m2+2*sz2*Psi; sources[l]=(Msrc,ph/sz,2*Psi*(1-Ph),m3-3*m2*Psi+2*Psi**3,Ph)
        m=m_next; Sig=C
    return np.array(out)
A=np.array([coefs(k) for k in range(8)])   # (net, layer, 3)
print("Reduced memory kernel fitted per network, basis {t*sigma^3, sigma^3}, plus the kurtosis constant.")
print("Fitted on a 2400-row Monte-Carlo sample, source depth K=1.\n")
print("layer |    a = coef of t*sigma^3    |    b = coef of sigma^3     |   g4 = mean excess kurtosis")
print("      |  mean     sd     sd/|mean|  |  mean     sd    sd/|mean| |  mean     sd    sd/|mean|")
for l in range(L-1):
    row=""
    for j in range(3):
        v=A[:,l,j]; mu_=v.mean(); sd=v.std()
        row+=f"| {mu_:+.5f} {sd:.5f}  {sd/max(abs(mu_),1e-12):6.3f} "
    print(f"{l+2:5d} "+row)
print("\nNoise floor check: same network, three different Monte-Carlo seeds.")
B=np.array([coefs(0,seed=s) for s in [0,1,2]])
for l in [4,9,14]:
    print(f"  layer {l+2:2d}: a = "+", ".join(f"{x:+.5f}" for x in B[:,l,0])+f"   (across-network sd {A[:,l,0].std():.5f})")
