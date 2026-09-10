# Is the two-point kernel coefficient universal in depth, like the diagonal one?
# Fitted from the deployed pipeline's own Monte-Carlo residual along the mean direction:
#   u_resid = u_MC - K21_model @ muhat,   regress u_resid on sigma^3  ->  c_l
import numpy as np
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=T.load_nets()
SCHED=np.array([[0.00094,0.00715],[0.00783,0.01340],[0.01134,0.01947],[0.01406,0.02302],[0.01658,0.02679],
 [0.01853,0.03065],[0.02040,0.03478],[0.02217,0.03824],[0.02331,0.04162],[0.02464,0.04615],
 [0.02555,0.04940],[0.02684,0.05096],[0.02781,0.05280],[0.02838,0.05462],[0.02940,0.05873]])
def coefs(k, N_MC=2400, seed=0):
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
        Ms,ps,qs,k3hs,d1s=sources[l-1]; U=Wl.astype(np.float64); Nt=Ms.T@U; UN=U*Nt
        k3=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
        X1=(UN*2)*ps[:,None]+(U*U)*qs[:,None]; X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(U*U)*k3hs[:,None]
        K21m=X1.T@Nt+X2.T@U; np.fill_diagonal(K21m,0.0)
        k3=k3+SCHED[l-1,0]*t*s3; g4=SCHED[l-1,1]
        z=Zs[l].astype(np.float64); a=z.mean(0); zc=z-a; v=(zc*zc).mean(0)
        k3s=(zc**3).mean(0)
        resid=k3s-k3; noise=6*np.mean(s3*s3)/N_MC; sig=max(np.mean(resid**2)-noise,0); k3=k3+(sig/(sig+noise))*resid
        muh=mu/np.linalg.norm(mu); Ap=zc@muh
        u_mc=(zc*zc).T@Ap/N_MC; u_res=u_mc-K21m@muh
        c=float(np.dot(u_res,s3)/np.dot(s3,s3)); r2=1-np.sum((u_res-c*s3)**2)/np.sum(u_res**2)
        out.append((c,r2))
        K21=K21m+np.outer(u_res,muh); np.fill_diagonal(K21,0.0)
        F2=ph/sz
        k4v=g4*s4; t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)+k4v*he2*ph/(24*s3)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*s4)-k4v*t*ph/(12*sz2)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rhoc=Sz/so; ds=T.hermite_d(t,ph,Ph,8); rk=rhoc.copy(); C=np.zeros((n,n)); fact=1.0
        for kk in range(1,9):
            fact*=kk; C+=rk*np.outer(ds[kk-1],ds[kk-1])/fact
            if kk<8: rk*=rhoc
        C*=so; C=C+0.5*(K21*np.outer(F2,Ph)+K21.T*np.outer(Ph,F2))
        C[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
        Msrc=Sz*Ph[:,None]; Msrc[np.diag_indices(n)]=0
        m3=mu*m2+2*sz2*Psi; sources[l]=(Msrc,ph/sz,2*Psi*(1-Ph),m3-3*m2*Psi+2*Psi**3,Ph)
        m=m_next; Sig=C
    return np.array(out)
A=np.array([coefs(k) for k in range(8)])
B=np.array([coefs(0,seed=s) for s in [0,1,2,3]])
print("Two-point kernel coefficient c_l:  u_resid ~ c_l * sigma^3, fitted from the pipeline's own MC residual.\n")
print("layer |  c across 8 networks: mean     sd    sd/|mean| | R2 of the sigma^3 fit | c across 4 seeds on net 0: sd")
for l in range(L-1):
    v=A[:,l,0]; r=A[:,l,1]; sv=B[:,l,0]
    print(f"{l+2:5d} |  {v.mean():+.5f}  {v.std():.5f}   {v.std()/max(abs(v.mean()),1e-12):6.3f}  |        {r.mean():.3f}          |  {sv.std():.5f}")
print("\nshipped schedule (mean over the 8 fitting networks):")
print("    C2P = (" + ", ".join(f"{x:.5f}" for x in A[:,:,0].mean(0)) + ")")
