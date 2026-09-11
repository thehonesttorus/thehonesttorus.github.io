# The (a,a,b,b) pattern paid.  The same expansion has more free terms:
#   (a,a,a,b), multiplicity 4!/(3!1!) = 4:
#       dCov = (1/6) kappa_aaab E[relu'''(z_a)] E[relu'(z_b)] = -(1/6) kappa_aaab (t phi/s^2)_a Phi_b
#       with kappa_aaab ~ g4 s_a^3 s_b  ->  -(g4/6) (s t phi)_a (s Phi)_b, symmetrised.
#   second order in kappa_3, the cross term kappa_aab kappa_abb:
#       dCov ~ c (K21 . K21^T) (x) (t phi/s^2)_a (t phi/s^2)_b   -- K21 is already computed.
# Every one of these is n^2 flops.  Greedy forward selection with the shipped term in place.
import numpy as np
from scipy.special import ndtr
import ablate as AB, twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=AB.seeds,AB.Yall
A_S,G_S,C_S=AB.A_S,AB.G_S,AB.C_S
def predict(k,co=None,N=2400,seed=0):
    co=co or {}
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
        z=Zs[l]; a=z.mean(0); zc=z-a; k3s=(zc**3).mean(0).astype(np.float64)
        rr=k3s-k3; nz=6*np.mean(s3*s3)/N; sg=max(np.mean(rr**2)-nz,0); k3=k3+(sg/(sg+nz))*rr
        muh=mu/np.linalg.norm(mu); Ap=zc@muh.astype(np.float32)
        u=0.75*(C_S[l-1]*s3)+0.25*(((zc*zc).T@Ap).astype(np.float64)/N-K21@muh)
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
        v1=sz*ph; v2=sz*Ph; v3=sz*t*ph; g=t*ph/sz2
        add=np.zeros((n,n))
        if co.get('T1'): add=add+co['T1']*np.outer(v1,v1)
        if co.get('T2'): add=add+co['T2']*(np.outer(v3,v2)+np.outer(v2,v3))
        if co.get('T3'): add=add+co['T3']*np.outer(v3,v3)
        if co.get('T4'): add=add+co['T4']*(np.outer(v1,v2)+np.outer(v2,v1))
        if co.get('T5'): add=add+co['T5']*(K21*K21.T)*np.outer(g,g)
        if co.get('T6'): add=add+co['T6']*(Sz*Sz)*np.outer(F2,F2)
        if co.get('T7'): add=add+co['T7']*Sz*np.outer(v1,v1)
        if add.any(): Cn=Cn+add.astype(np.float32)
        Cn[np.diag_indices(n)]=((m2+dE2)-m_next*m_next).astype(np.float32)
        Msrc=Sz*Ph[:,None]; np.fill_diagonal(Msrc,0.0); m3=mu*m2+2*sz2*Psi
        sources[l]=(Msrc.astype(np.float32),(ph/sz).astype(np.float32),(2*Psi*(1-Ph)).astype(np.float32),
                    (m3-3*m2*Psi+2*Psi**3).astype(np.float32),Ph.astype(np.float32))
        Sig=Cn.astype(np.float32); m=m_next.astype(np.float32)
    return float(np.mean((m.astype(np.float64)-Y[-1])**2))
if __name__=="__main__":
    FIT=[0,2,4,6]; HOLD=list(range(8,16))
    def sc(nets,co): return float(np.mean([predict(k,co) for k in nets]))
    cur={'T1':0.0020}; b=sc(FIT,cur)
    print(f"shipped (T1=0.0020):  fit {b:.4e}   held {sc(HOLD,cur):.4e}\n")
    GRID={'T2':(-0.004,-0.0012,0.0012,0.004),'T3':(-0.006,-0.002,0.002,0.006),
          'T4':(-0.0012,-0.0004,0.0004,0.0012),'T5':(-3.0,-1.0,1.0,3.0),
          'T6':(-0.2,0.2,0.6),'T7':(-0.02,-0.006,0.006,0.02)}
    for rnd in (1,2):
        print(f"--- round {rnd} ---   base fit {b:.4e}")
        best=(b,None)
        for name,gs in GRID.items():
            if name in cur: continue
            for v in gs:
                c=dict(cur); c[name]=v; f=sc(FIT,c)
                mark=" *" if f<best[0] else ""
                if f<best[0]: best=(f,(name,v))
                print(f"   {name}={v:<9.4g} fit {f:.4e}  {f/b:7.4f}{mark}", flush=True)
        if best[1] is None or best[0]>=b*0.999:
            print("   no further term improves the fit"); break
        cur[best[1][0]]=best[1][1]; b=best[0]
        print(f"   keep {best[1][0]}={best[1][1]}   fit {b:.4e}   held {sc(HOLD,cur):.4e}", flush=True)
    print(f"\nfinal {cur}")
    print(f"fit {sc(FIT,cur):.4e}   held {sc(HOLD,cur):.4e}")
