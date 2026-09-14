# The Mac audit measures a percent-level Price defect in a Gaussian closure.  Mine agrees once
# normalised the same way.  The question that decides whether it matters for the SCORE is where the
# defect lives: our error budget (section 35) says each layer contributes about the same 6e-5 rms to
# the final mean, so a defect concentrated elsewhere is not the carrier.
# For each layer l, take the deployed trajectory's (m_l, Sigma_l), and measure the Price defect of the
# map from that state to the FINAL mean:   P_v = D_Sigma[v v^T] M - (1/2) D^2_m[v,v] M.
# Absolute norm, not just the ratio -- the ratio hides how much error it can actually inject.
import numpy as np
from scipy.special import ndtr
import ablate as AB, twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
A_S,G_S,C_S=AB.A_S,AB.G_S,AB.C_S
A4=0.0020; S3=0.22
W=[w for w in AB.getW(0)]; Y=AB.Yall[0]

def suffix(m0,S0,start):
    """Deployed closure (including the two shipped covariance terms) from `start` to the end."""
    m=m0.astype(np.float64).copy(); Sg=S0.astype(np.float64).copy(); src=None
    for l in range(start,L):
        Wl=W[l].astype(np.float64); mu=m@Wl; Sz=(Wl.T@Sg)@Wl
        sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        k3=np.zeros(n); K21=np.zeros((n,n))
        if src is not None:
            Ms,ps,qs,k3hs=src; Nt=Ms.T@Wl; UN=Wl*Nt
            k3=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*Wl+k3hs[:,None]*(Wl*Wl*Wl)).sum(0)
            X1=(UN*2)*ps[:,None]+(Wl*Wl)*qs[:,None]
            X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(Wl*Wl)*k3hs[:,None]
            K21=X1.T@Nt+X2.T@Wl; np.fill_diagonal(K21,0.0)
        i=min(l-1,14); g4=G_S[i]
        k3=k3+A_S[i]*t*s3
        nm=np.linalg.norm(mu); muh=mu/nm if nm>0 else mu
        K21=K21+np.outer(C_S[i]*s3,muh); np.fill_diagonal(K21,0.0)
        k4v=g4*s4; t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)+k4v*he2*ph/(24*s3)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*s4)-k4v*t*ph/(12*sz2)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rho=Sz/so; ds=T.hermite_d(t,ph,Ph,8); rk=rho.copy(); C=np.zeros((n,n)); fact=1.0
        for kk in range(1,9):
            fact*=kk; C+=rk*np.outer(ds[kk-1],ds[kk-1])/fact
            if kk<8: rk*=rho
        C*=so; F2=ph/sz
        C=C+0.5*(K21*np.outer(F2,Ph)+K21.T*np.outer(Ph,F2))
        v4=sz*ph; tp=t*ph
        C=C+A4*np.outer(v4,v4)
        C=C+S3*g4*Sz*(np.outer(tp,Ph)+np.outer(Ph,tp))
        C[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
        Msrc=Sz*Ph[:,None]; np.fill_diagonal(Msrc,0.0); m3=mu*m2+2*sz2*Psi
        src=(Msrc,ph/sz,2*Psi*(1-Ph),m3-3*m2*Psi+2*Psi**3)
        m=m_next; Sg=C
    return m

# deployed trajectory states
def trajectory():
    st={}
    W0=W[0].astype(np.float64); G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1)
    C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi))
    st[0]=(m.copy(),C.copy())
    for l in range(1,L):
        # advance one layer using the same code path
        mm,SS=st[l-1]
        # reuse suffix() for one step by calling it with start=l-1 on a truncated horizon
        st[l]=None
    return st
# simpler: run the full closure once, recording (m, Sigma) before each layer
def states():
    out={}
    W0=W[0].astype(np.float64); G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1)
    C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi))
    Msrc0=(G*0.5); np.fill_diagonal(Msrc0,0.0)
    src=(Msrc0,INV/s,m,s2*s*(2*INV-1.5*INV+2*INV**3)); Sg=C
    out[1]=(m.copy(),Sg.copy())
    for l in range(1,L-1):
        Wl=W[l].astype(np.float64); mu=m@Wl; Sz=(Wl.T@Sg)@Wl
        sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        Ms,ps,qs,k3hs=src; Nt=Ms.T@Wl; UN=Wl*Nt
        k3=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*Wl+k3hs[:,None]*(Wl*Wl*Wl)).sum(0)
        X1=(UN*2)*ps[:,None]+(Wl*Wl)*qs[:,None]
        X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(Wl*Wl)*k3hs[:,None]
        K21=X1.T@Nt+X2.T@Wl; np.fill_diagonal(K21,0.0)
        i=l-1; g4=G_S[i]; k3=k3+A_S[i]*t*s3
        nm=np.linalg.norm(mu); muh=mu/nm
        K21=K21+np.outer(C_S[i]*s3,muh); np.fill_diagonal(K21,0.0)
        k4v=g4*s4; t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)+k4v*he2*ph/(24*s3)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*s4)-k4v*t*ph/(12*sz2)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rho=Sz/so; ds=T.hermite_d(t,ph,Ph,8); rk=rho.copy(); C=np.zeros((n,n)); fact=1.0
        for kk in range(1,9):
            fact*=kk; C+=rk*np.outer(ds[kk-1],ds[kk-1])/fact
            if kk<8: rk*=rho
        C*=so; F2=ph/sz
        C=C+0.5*(K21*np.outer(F2,Ph)+K21.T*np.outer(Ph,F2))
        v4=sz*ph; tp=t*ph
        C=C+A4*np.outer(v4,v4)+S3*g4*Sz*(np.outer(tp,Ph)+np.outer(Ph,tp))
        C[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
        Msrc=Sz*Ph[:,None]; np.fill_diagonal(Msrc,0.0); m3=mu*m2+2*sz2*Psi
        src=(Msrc,ph/sz,2*Psi*(1-Ph),m3-3*m2*Psi+2*Psi**3)
        m=m_next; Sg=C
        out[l+1]=(m.copy(),Sg.copy())
    return out

ST=states()
print("Price defect of the deployed closure, state at layer l -> FINAL mean.  Network 0.")
print("The 'inject' column is ||P_v|| scaled by the perturbation size actually present: the per-layer")
print("relative sigma error measured in section 35, so it is comparable to the 6e-5 rms per-layer")
print("contribution that makes up our whole error.\n")
SIG_ERR={1:0.0,2:2.6e-4,3:7.3e-4,4:1.01e-3,5:1.37e-3,6:1.71e-3,7:1.99e-3,8:2.34e-3,9:2.67e-3,
         10:2.99e-3,11:3.19e-3,12:3.56e-3,13:3.95e-3,14:4.40e-3,15:4.75e-3}
print(f"{'layer':>5} {'||dS||':>10} {'||d2/2||':>10} {'||P_v||':>10} {'P/resp':>8} {'P/||M||':>10}")
for l in range(1,15):
    m0,S0=ST[l]
    rs=np.random.default_rng(100+l); v=rs.standard_normal(n); v/=np.linalg.norm(v)
    s=t=0.05; G=np.outer(v,v)
    base=suffix(m0,S0,l)
    dS=(suffix(m0,S0+t*G,l)-suffix(m0,S0-t*G,l))/(2*t)
    d2=(suffix(m0+s*v,S0,l)-2*base+suffix(m0-s*v,S0,l))/(s*s)
    P=dS-0.5*d2; rp=max(np.linalg.norm(dS),np.linalg.norm(0.5*d2))
    print(f"{l:5d} {np.linalg.norm(dS):10.3e} {np.linalg.norm(0.5*d2):10.3e} {np.linalg.norm(P):10.3e}"
          f" {np.linalg.norm(P)/rp:8.2%} {np.linalg.norm(P)/np.linalg.norm(base):10.3e}", flush=True)
