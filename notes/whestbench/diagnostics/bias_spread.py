# Is the sigma deficit a coherent per-layer offset, or per-neuron scatter?  And does an INCOHERENT
# perturbation amplify like a coherent one?  Both questions have to be answered before the budget in
# bias_signs.py means anything.
import numpy as np
import varcal as V, ablate as AB
import bias_signs as B          # reuses its deployed trace
n=1024; L=16; NB=1000000; NETS=[0,1]
SD={k:np.load(f"mom_{k}_{NB}.npz")['sd'] for k in NETS}
print("Per-neuron distribution of  sigma_model/sigma_true - 1  along the deployed trajectory.\n")
print(f"{'layer':>5} {'mean':>12} {'rms':>12} {'sd':>12} {'coherent share':>15}")
for i in range(L-1):
    ms=[];rs=[];ss=[]
    for k in NETS:
        r=B.R[k][i][2]/SD[k][i+1]-1.0
        ms.append(r.mean()); rs.append(np.sqrt(np.mean(r*r))); ss.append(r.std())
    m,rr,s=np.mean(ms),np.mean(rs),np.mean(ss)
    print(f"{i+1:5d} {m:12.3e} {rr:12.3e} {s:12.3e} {(m/rr)**2:15.4f}")
# amplification of an INCOHERENT per-neuron sigma perturbation of the same rms
print("\nAmplification of a per-neuron RANDOM sigma perturbation (rms epsilon = 7e-4) against the")
print("uniform one measured earlier.  Sz <- Sz * (1+v_b) with v_b = 2 epsilon xi_b, xi ~ N(0,1).\n")
def predict_rand(k,eps,layers,rng_seed=11):
    import numpy as np
    from scipy.special import ndtr
    import twin2 as T
    INV=T.INV; pdf=T.pdf; A_S,G_S,C_S=AB.A_S,AB.G_S,AB.C_S
    W=AB.getW(k); Y=AB.Yall[k]; N=2400
    rg=np.random.default_rng(0); X=rg.standard_normal((N,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    pr=np.random.default_rng(rng_seed)
    sources={}
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); C=C.astype(np.float32)
    M=(G*0.5).astype(np.float32); M[np.diag_indices(n)]=0
    sources[0]=(M,(INV/s).astype(np.float32),m.astype(np.float32),
                (s2*s*(2*INV-1.5*INV+2*INV**3)).astype(np.float32),np.full(n,0.5,np.float32))
    Sig=C; m=m.astype(np.float32)
    for l in range(1,L):
        Wl=W[l]; mu=(m@Wl).astype(np.float64); Sz=((Wl.T@Sig)@Wl).astype(np.float64)
        if l in layers:
            f=1.0+2*eps*pr.standard_normal(n)
            Sz=Sz*np.sqrt(np.outer(f,f))
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
        Cn[np.diag_indices(n)]=((m2+dE2)-m_next*m_next).astype(np.float32)
        Msrc=Sz*Ph[:,None]; np.fill_diagonal(Msrc,0.0); m3=mu*m2+2*sz2*Psi
        sources[l]=(Msrc.astype(np.float32),(ph/sz).astype(np.float32),(2*Psi*(1-Ph)).astype(np.float32),
                    (m3-3*m2*Psi+2*Psi**3).astype(np.float32),Ph.astype(np.float32))
        Sig=Cn.astype(np.float32); m=m_next.astype(np.float32)
    return float(np.mean((m.astype(np.float64)-Y[-1])**2))
NN=[0,1,2,3]
base=np.mean([V.predict(k,np.zeros(15)) for k in NN])
eps=7e-4
print(f"{'perturbed at':>14} {'uniform amp':>12} {'random amp':>12}")
UNI={1:0.83,5:0.45,10:0.23,15:0.05,'all':3.57}
for l0 in [1,5,10,15,'all']:
    ly=set(range(1,16)) if l0=='all' else {l0}
    m=np.mean([predict_rand(k,eps,ly) for k in NN])
    print(f"{str(l0):>14} {UNI[l0]:12.2f} {np.sqrt(max(m-base,0))/eps:12.2f}", flush=True)
