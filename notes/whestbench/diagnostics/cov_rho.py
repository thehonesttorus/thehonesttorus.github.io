# Section 35 says the whole error is the propagated sigma, injected at about 3.3e-4 per layer by the
# OFF-DIAGONAL covariance closure.  So complete that closure to the next order.  For Cov(h_a,h_b) only
# the joint law of the PAIR matters, and the cumulant expansion over index patterns in {a,b} gives,
# after subtracting E[h_a]E[h_b] (which kills every piece surviving at rho = 0):
#
#  (a,a,a) mult 1, third order.  (1/6) kappa_aaa E[relu'''(z_a) relu(z_b)], relu''' = delta'.
#     E[delta'(z_a) relu(z_b)] = -[ p_a'(0) Psi_b + p_a(0) Phi_b rho sigma_b/sigma_a ] and the first
#     piece cancels against E[h_a]E[h_b], leaving   -(1/6) kappa_aaa phi_a Phi_b Sz_ab / sigma_a^3,
#     i.e.  S1 = -(1/6) Sz . [ (lambda_3 phi) (x) Phi  +  Phi (x) (lambda_3 phi) ].
#  (a,a,b) mult 3.  We already keep the leading (1/2) kappa_aab F2_a Phi_b.  Its next order in rho
#     comes from Phi(m(0)/s) = Phi(t_b - rho t_a) ~ Phi_b - rho t_a phi_b:
#     S2 = -(1/2) K21 . Sz . [ (t phi/sigma^2) (x) (phi/sigma) ] + transpose.
#  (a,a,a,b) mult 4, fourth order.  (1/6) kappa_aaab E[relu'''(z_a) relu'(z_b)] with kappa_aaab
#     vanishing at rho = 0; the ansatz kappa_aaab = g4 sigma_a^2 Sz_ab gives
#     S3 = -(g4/6) Sz . [ (t phi) (x) Phi + Phi (x) (t phi) ].
# Every one of these is elementwise: n^2 flops.  Coefficients are PREDICTED, then fitted.
import numpy as np
from scipy.special import ndtr
import ablate as AB, twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=AB.seeds,AB.Yall
A_S,G_S,C_S=AB.A_S,AB.G_S,AB.C_S
A4=0.0020
def predict(k,S1=0.0,S2=0.0,S3=0.0,S4=0.0,a4=A4,N=2400,seed=0):
    W=AB.getW(k); Y=Yall[k]
    rng=np.random.default_rng(seed); X=rng.standard_normal((N,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    sources={}
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
        u=0.75*(C_S[l-1]*s3)+0.25*((zc2.T@Ap).astype(np.float64)/N-K21@muh)
        K21=K21+np.outer(u,muh); np.fill_diagonal(K21,0.0)
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
        v4=sz*ph
        Cn=Cn+(a4*np.outer(v4,v4)).astype(np.float32)
        if S1:
            g=(k3/s3)*ph
            Cn=Cn+(S1*Sz*(np.outer(g,Ph)+np.outer(Ph,g))).astype(np.float32)
        if S2:
            w1=t*ph/sz2; w2=ph/sz
            Cn=Cn+(S2*Sz*(K21*np.outer(w1,w2)+K21.T*np.outer(w2,w1))).astype(np.float32)
        if S3:
            tp=t*ph
            Cn=Cn+(S3*g4*Sz*(np.outer(tp,Ph)+np.outer(Ph,tp))).astype(np.float32)
        if S4:
            tp=t*ph
            Cn=Cn+(S4*Sz*(np.outer(tp,Ph)+np.outer(Ph,tp))).astype(np.float32)
        Cn[np.diag_indices(n)]=((m2+dE2)-m_next*m_next).astype(np.float32)
        Msrc=Sz*Ph[:,None]; np.fill_diagonal(Msrc,0.0); m3=mu*m2+2*sz2*Psi
        sources[l]=(Msrc.astype(np.float32),(ph/sz).astype(np.float32),(2*Psi*(1-Ph)).astype(np.float32),
                    (m3-3*m2*Psi+2*Psi**3).astype(np.float32),Ph.astype(np.float32))
        Sig=Cn.astype(np.float32); m=m_next.astype(np.float32)
    return float(np.mean((m.astype(np.float64)-Y[-1])**2))
if __name__=="__main__":
    FIT=[0,2,4,6]; HOLD=list(range(8,16))
    def sc(nets,**kw): return float(np.mean([predict(k,**kw) for k in nets]))
    b=sc(FIT); bh=sc(HOLD)
    print(f"shipped (A4 only):  fit {b:.4e}   held {bh:.4e}\n")
    print(f"{'spec':28s}  {'fit':>11} {'ratio':>8} {'held':>11} {'held ratio':>11}")
    best=(b,{})
    for name,gs in (('S1',(0.5,1.0,1.5,2.5)),
                    ('S3',(0.12,0.17,0.25,0.4)),
                    ('S4',(0.002,0.005,0.01,0.02,0.04))):
        for v in gs:
            kw={name:v}; f=sc(FIT,**kw)
            if f<best[0]: best=(f,kw)
            print(f"{name+'='+format(v,'.4g'):28s}  {f:11.4e} {f/b:8.4f}", flush=True)
    print(f"\nbest single term {best[1]}  fit {best[0]:.4e}  held {sc(HOLD,**best[1]):.4e}")
    # combine the best single term with the other family
    for extra,gs in (('S4',(0.004,0.008)),('S1',(0.4,0.8)),('S3',(0.08,0.15))):
        if extra in best[1]: continue
        for v in gs:
            kw=dict(best[1]); kw[extra]=v; f=sc(FIT,**kw)
            print(f"{str(kw):28s}  {f:11.4e} {f/b:8.4f} {sc(HOLD,**kw):11.4e}", flush=True)
