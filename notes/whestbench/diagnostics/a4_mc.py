# The fourth-cumulant coefficient is NOT universal (per-network optima span 0.0008 to >0.0034), so
# estimate it from the Monte-Carlo sample we already pay for.  For centred variables
#     Cov(z_a^2, z_b^2) = kappa_aabb + 2 Cov(z_a,z_b)^2,
# so with s_i = sum_a v_a zc_{i,a}^2 and v = sigma phi,
#     Var(s) = sum_ab v_a v_b kappa_aabb + 2 v^T (Sz . Sz) v,
# and under the ansatz kappa_aabb = g4 sigma_a^2 sigma_b^2 that first term is g4 (v . sigma^2)^2.  So
#     g4_est = [ Var(s) - 2 v^T (Sz.Sz) v ] / (v . sigma^2)^2,     A4 = g4_est/4
# at 2 N n + 3 n^2 flops per layer, i.e. 0.006% of budget.  Blend with the shipped constant for noise.
import numpy as np
from scipy.special import ndtr
import ablate as AB, twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=AB.seeds,AB.Yall
A_S,G_S,C_S=AB.A_S,AB.G_S,AB.C_S
def predict(k,A4=0.0020,w=0.0,lo=0.0,hi=0.02,N=2400,seed=0,report=False):
    W=AB.getW(k); Y=Yall[k]
    rng=np.random.default_rng(seed); X=rng.standard_normal((N,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    sources={}; est=[]
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so0=np.outer(s,s)
    rho=np.clip(G/so0,-1,1); C=so0*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); C=C.astype(np.float32)
    M=(G*0.5).astype(np.float32); M[np.diag_indices(n)]=0
    sources[0]=(M,(INV/s).astype(np.float32),m.astype(np.float32),
                (s2*s*(2*INV-1.5*INV+2*INV**3)).astype(np.float32),np.full(n,0.5,np.float32))
    Sig=C; m=m.astype(np.float32)
    for l in range(1,L):
        Wl=W[l]; mu=(m@Wl).astype(np.float64); Sz=((Wl.T@Sig)@Wl).astype(np.float64)
        sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        U=Wl; Ms,ps,qs,k3hs,d1s=sources[l-1]; Nt=Ms.T@U; UN=U*Nt
        k3=((ps[:,None]*np.float32(3))*UN*Nt+(qs[:,None]*np.float32(3))*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
        X1=(UN*2)*ps[:,None]+(U*U)*qs[:,None]
        X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(U*U)*k3hs[:,None]
        K21=(X1.T@Nt+X2.T@U).astype(np.float64); np.fill_diagonal(K21,0.0)
        k3=k3.astype(np.float64)+A_S[l-1]*t*s3; g4=G_S[l-1]
        z=Zs[l]; a=z.mean(0); zc=z-a; zc2=zc*zc; k3s=(zc2*zc).mean(0).astype(np.float64)
        rr=k3s-k3; nz=6*np.mean(s3*s3)/N; sg=max(np.mean(rr**2)-nz,0); k3=k3+(sg/(sg+nz))*rr
        muh=mu/np.linalg.norm(mu); Ap=zc@muh.astype(np.float32)
        u=0.75*(C_S[l-1]*s3)+0.25*(((zc2).T@Ap).astype(np.float64)/N-K21@muh)
        K21=K21+np.outer(u,muh); np.fill_diagonal(K21,0.0)
        v4=sz*ph
        a4=A4
        if w>0.0:
            sv=(zc2@v4.astype(np.float32)).astype(np.float64)
            q=float(v4@((Sz*Sz)@v4)); den=float(v4@sz2)**2
            g4e=(float(sv.var()) - 2.0*q)/den
            a4e=min(max(g4e/4.0,lo),hi)
            a4=w*a4e+(1.0-w)*A4
            est.append(a4e)
        k4v=g4*s4; t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)+k4v*he2*ph/(24*s3)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*s4)-k4v*t*ph/(12*sz2)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rhoc=Sz/so; ds=T.hermite_d(t,ph,Ph,8); rk=rhoc.copy()
        Cn=np.zeros((n,n),np.float32); fact=1.0
        for kk in range(1,9):
            fact*=kk; Cn=Cn+(rk*np.outer(ds[kk-1],ds[kk-1])/fact).astype(np.float32)
            if kk<8: rk*=rhoc
        Cn=Cn*so.astype(np.float32); F2=ph/sz
        Cn=Cn+(0.5*(K21*np.outer(F2,Ph)+K21.T*np.outer(Ph,F2))).astype(np.float32)
        Cn=Cn+(a4*np.outer(v4,v4)).astype(np.float32)
        Cn[np.diag_indices(n)]=((m2+dE2)-m_next*m_next).astype(np.float32)
        Msrc=Sz*Ph[:,None]; np.fill_diagonal(Msrc,0.0); m3=mu*m2+2*sz2*Psi
        sources[l]=(Msrc.astype(np.float32),(ph/sz).astype(np.float32),(2*Psi*(1-Ph)).astype(np.float32),
                    (m3-3*m2*Psi+2*Psi**3).astype(np.float32),Ph.astype(np.float32))
        Sig=Cn.astype(np.float32); m=m_next.astype(np.float32)
    if report: return float(np.mean((m.astype(np.float64)-Y[-1])**2)), est
    return float(np.mean((m.astype(np.float64)-Y[-1])**2))
if __name__=="__main__":
    import sys
    _,e=predict(0,w=1.0,report=True)
    print("net 0, raw per-layer A4 estimates:", " ".join(f"{x:.4f}" for x in e),"\n")
    _,e1=predict(5,w=1.0,report=True)
    print("net 5 (wanted a large A4):        ", " ".join(f"{x:.4f}" for x in e1),"\n")
    FIT=[0,2,4,6]; HOLD=list(range(8,16))
    def sc(nets,**kw): return float(np.mean([predict(k,**kw) for k in nets]))
    b=sc(FIT,w=0.0); bh=sc(HOLD,w=0.0)
    print(f"shipped constant 0.0020:  fit {b:.4e}  held {bh:.4e}\n")
    print(f"{'blend w':>8} {'clip hi':>8}  {'fit':>11} {'held':>11} {'held ratio':>11}")
    for hi in (0.006,0.012):
        for w in (0.25,0.5,0.75,1.0):
            f=sc(FIT,w=w,hi=hi); h=sc(HOLD,w=w,hi=hi)
            print(f"{w:8.2f} {hi:8.3f}  {f:11.4e} {h:11.4e} {h/bh:11.4f}", flush=True)
