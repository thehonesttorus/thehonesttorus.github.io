# Oracle test of two-point cumulant corrections to the off-diagonal covariance closure (network 0)
import numpy as np, math
from scipy.special import ndtr
W=np.load('data/W0.npy').astype(np.float64); Y=np.load('data/Y0.npy').astype(np.float64); n=1024
S=np.load('pass2_stats.npz'); HS=S['HS']; P3=np.load('pass3_stats.npz'); S1=P3['S1']; S2=P3['S2']; S11=P3['S11']; S21=P3['S21']; S22=P3['S22']; S31=P3['S31']
pdf=lambda t: np.exp(-0.5*t*t)/np.sqrt(2*np.pi)
def He(k,x):
    a,b=np.ones_like(x),x
    if k==0: return a
    for m in range(1,k): a,b=b,x*b-m*a
    return b
print("layer | induced next-layer sigma2 rel err (mean, rms): Gauss | +k3 two-point | +k3 mixed | +k4 two-point | aligned-fit residual rms (G / full) | max|rho|")
for l in range(1,15):
    mu=Y[l-1]@W[l]; Szz=W[l].T@HS[l]@W[l]; s=np.sqrt(np.diag(Szz)); t=mu/s; ph=pdf(t); Ph=ndtr(t); rho=Szz/np.outer(s,s); np.fill_diagonal(rho,0)
    # Gaussian closure covariance (Hermite series, exact diagonal)
    ds=[Ph]+[ph*He(k-2,-t) for k in range(2,13)]
    C=np.zeros((n,n)); rk=np.ones((n,n))
    for k in range(1,13): rk=rk*rho; C+=rk*np.outer(ds[k-1],ds[k-1])/math.factorial(k)
    C*=np.outer(s,s)
    m1=mu*Ph+s*ph; m2=(mu*mu+s*s)*Ph+mu*s*ph; np.fill_diagonal(C,m2-m1*m1)
    # two-point cumulants from pass3 (centered with pass3 moments)
    a=S1[l]; Cz=S11[l]-np.outer(a,a); va=np.diag(Cz)
    K21=S21[l]-2*a[:,None]*S11[l]-a[None,:]*S2[l][:,None]+2*(a*a)[:,None]*a[None,:]          # kappa_aab = E[zt_a^2 zt_b]
    # E[zt_a^2 zt_b^2] = S22 - 2 a_b S21_ab - 2 a_a S21_ba + a_b^2 S2_a + a_a^2 S2_b + 4 a_a a_b S11_ab - 3 a_a^2 a_b^2
    M22=S22[l]-2*a[None,:]*S21[l]-2*a[:,None]*S21[l].T+(a*a)[None,:]*S2[l][:,None]+(a*a)[:,None]*S2[l][None,:]+4*np.outer(a,a)*S11[l]-3*np.outer(a*a,a*a)
    K22=M22-np.outer(va,va)-2*Cz*Cz
    # E[zt_a^3 zt_b] = S31 - 3 a_a S21_ab - a_b S3_a... need S3 (per-neuron third raw): use S21 diagonal: S21_aa = E[z_a^3]
    S3=np.diag(S21[l]); M31=S31[l]-3*a[:,None]*S21[l]-a[None,:]*S3[:,None]+3*(a*a)[:,None]*S11[l]+3*(a*a)[:,None]*a[None,:]*a[:,None]*0  # fix below
    # exact: E[(z_a-a_a)^3 (z_b-a_b)] = S31 - a_b S3_a - 3 a_a S21_ab + 3 a_a a_b S2_a + 3 a_a^2 S11_ab - 3 a_a^2 a_b a_a ... derive cleanly:
    # (z_a-a)^3 (z_b-b) = (z_a^3 -3a z_a^2 +3a^2 z_a - a^3)(z_b - b)
    M31=(S31[l]-3*a[:,None]*S21[l]+3*(a*a)[:,None]*S11[l]-(a**3)[:,None]*a[None,:]) - a[None,:]*(S3[:,None]-3*a[:,None]*S2[l][:,None]+3*(a*a)[:,None]*a[:,None]-(a**3)[:,None])
    K31=M31-3*va[:,None]*Cz
    k3d=np.diag(K21).copy(); k4d=np.diag(K22)/1.0  # diag of K22 = E[zt^4]-3v^2-2v^2?? careful: K22_aa = E[zt^4]-v^2-2v^2 = k4 ; fine
    k4d=np.diag(M22)-3*va*va
    F1=Ph; F2=ph/s; F3=-t*ph/(s*s); F4=(t*t-1)*ph/s**3
    d3=0.5*(K21*np.outer(F2,F1)+K21.T*np.outer(F1,F2))
    # mixed terms: (k3_a/6) Cov_G(F_a''',G_b) ~ (k3_a/6) rho_ab c1(F''')_a c1(G)_b with c1(F''')=F4*s_a, c1(G)=s_b*Phi_b
    dmix=(k3d/6)[:,None]*rho*np.outer(F4*s,s*Ph)+(k3d/6)[None,:]*rho*np.outer(s*Ph,F4*s)
    d4=0.25*K22*np.outer(F2,F2)+(1/6)*(K31*np.outer(F3,F1)+K31.T*np.outer(F1,F3))
    for D in (d3,dmix,d4): np.fill_diagonal(D,0)
    Wn=W[l+1]; Strue=HS[l+1]; Coff_true=Strue.copy(); np.fill_diagonal(Coff_true,0)
    base=np.einsum('ij,ik,kj->j',Wn,Strue,Wn); tn=(Y[l]@Wn)/np.sqrt(base)
    def induced(Cm):
        Coff=Cm.copy(); np.fill_diagonal(Coff,0); return (np.einsum('ij,ik,kj->j',Wn,Coff,Wn)-np.einsum('ij,ik,kj->j',Wn,Coff_true,Wn))/base
    res=[]; errs=[]
    for Cm in (C, C+d3, C+d3+dmix, C+d3+dmix+d4):
        e=induced(Cm); res.append((e.mean(),np.sqrt(np.mean(e**2)))); errs.append(e)
    Xb=np.column_stack([np.ones(n),tn,tn*tn]); 
    def alres(e): cf=np.linalg.lstsq(Xb,e,rcond=None)[0]; return np.sqrt(np.mean((e-Xb@cf)**2))
    print(f"{l+1:5d} | "+" | ".join(f"{m:+.1e} {r:.1e}" for m,r in res)+f" | {alres(errs[0]):.1e} / {alres(errs[3]):.1e} | {np.abs(rho).max():.2f}", flush=True)
