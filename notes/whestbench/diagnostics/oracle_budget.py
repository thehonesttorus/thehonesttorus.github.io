# The oracle-substitution error budget.  We have asserted "the binding constraint is propagation" many
# times without ever measuring propagation of WHAT.  Run the closure to layer l, substitute the TRUE
# value of one component there, continue, and record the final-layer MSE as a function of l.
#   mean      : true layer means from the 1e9-sample ground truth
#   covariance: high-sample covariance of the post-activations
#   both      : the full second-order state
# Costs no metered FLOPs.
import numpy as np
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=T.load_nets()
HS=np.load('../pass2_stats.npz')['HS']
SCHED=np.array([[0.00094,0.00715],[0.00783,0.01340],[0.01134,0.01947],[0.01406,0.02302],[0.01658,0.02679],
 [0.01853,0.03065],[0.02040,0.03478],[0.02217,0.03824],[0.02331,0.04162],[0.02464,0.04615],
 [0.02555,0.04940],[0.02684,0.05096],[0.02781,0.05280],[0.02838,0.05462],[0.02940,0.05873]])
C2P=np.array([0.00396,0.07819,0.14538,0.22329,0.30105,0.38174,0.45643,0.55252,0.63162,0.70963,
              0.80605,0.92641,1.00444,1.09326,1.18779])
def predict(k, sub_layer=None, sub_what=None, N_MC=2400, seed=0):
    W=T.regen(seeds[k]); Y=Yall[k].astype(np.float64)
    rng=np.random.default_rng(seed); X=rng.standard_normal((N_MC,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    outs=[]; src=None
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); M=G*0.5; M[np.diag_indices(n)]=0
    src=(M,INV/s,m,s2*s*(2*INV-1.5*INV+2*INV**3)); outs.append(m); Sig=C
    for l in range(1,L):
        if sub_layer==l:
            if sub_what in ('mean','both'): m=Y[l-1].copy()
            if sub_what in ('cov','both'):  Sig=HS[l].copy()
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sig)@Wl; sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        Ms,ps,qs,k3hs=src
        U=Wl.astype(np.float64); Nt=Ms.T@U; UN=U*Nt
        k3=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
        X1=(UN*2)*ps[:,None]+(U*U)*qs[:,None]
        X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(U*U)*k3hs[:,None]
        K21=X1.T@Nt+X2.T@U; np.fill_diagonal(K21,0.0)
        k3=k3+SCHED[l-1,0]*t*s3; g4=SCHED[l-1,1]
        z=Zs[l].astype(np.float64); a=z.mean(0); zc=z-a; k3s=(zc**3).mean(0)
        resid=k3s-k3; noise=6*np.mean(s3*s3)/N_MC; sig=max(np.mean(resid**2)-noise,0); k3=k3+(sig/(sig+noise))*resid
        muh=mu/np.linalg.norm(mu); Ap=zc@muh
        u_mc=(zc*zc).T@Ap/N_MC-K21@muh
        u=0.75*(C2P[l-1]*s3)+0.25*u_mc
        K21=K21+np.outer(u,muh); np.fill_diagonal(K21,0.0)
        k4v=g4*s4; t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)+k4v*he2*ph/(24*s3)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*s4)-k4v*t*ph/(12*sz2)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rhoc=Sz/so; ds=T.hermite_d(t,ph,Ph,8); rk=rhoc.copy(); Cn=np.zeros((n,n)); fact=1.0
        for kk in range(1,9):
            fact*=kk; Cn+=rk*np.outer(ds[kk-1],ds[kk-1])/fact
            if kk<8: rk*=rhoc
        Cn*=so; F2=ph/sz
        Cn=Cn+0.5*(K21*np.outer(F2,Ph)+K21.T*np.outer(Ph,F2))
        Cn[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
        Msrc=Sz*Ph[:,None]; np.fill_diagonal(Msrc,0.0)
        m3=mu*m2+2*sz2*Psi
        src=(Msrc,ph/sz,2*Psi*(1-Ph),m3-3*m2*Psi+2*Psi**3)
        m=m_next; Sig=Cn; outs.append(m)
    return np.mean((np.array(outs)-Y)**2,axis=1)[-1]
base=predict(0)
print(f"Oracle-substitution error budget, network 0.  Baseline final-layer MSE {base:.4e}\n")
print("substitute at layer |   mean only   |   cov only    |    both       | (fraction of baseline)")
for l in range(1,L):
    a=predict(0,l,'mean'); b=predict(0,l,'cov'); c=predict(0,l,'both')
    print(f"        {l+1:2d}          | {a:.3e} ({a/base:.2f}) | {b:.3e} ({b/base:.2f}) | {c:.3e} ({c/base:.2f})", flush=True)
