# A coherent (single positive-extension) parameterisation makes a SHARP prediction about the functional
# form of the kernels we fitted empirically.  Binary extension X|T ~ N(m+Tu, Sigma-uu^T+TB) gives
#   kappa_3 = R(u,B):   diagonal  3 u_a B_aa ,   two-point  2 u_a B_ab + u_b B_aa
# Taking B aligned with the pre-activation covariance, B = beta*Sz, and u = alpha*muhat (since
# t*sigma = mu, and our diagonal kernel's fitted form a*t*sigma^3 forces u proportional to mu):
#   diagonal   -> 3 alpha beta muhat_a sz2_a  ~  t sigma^3          MATCHES our fitted form
#   two-point  -> alpha beta [ 2 muhat_a Sz_ab + muhat_b sz2_a ]
# So coherence predicts the two-point kernel should go as sigma^2 (not the sigma^3 we fitted), and
# should carry an extra term 2 muhat_a Sz_ab that we do not have at all.  Test both.
import numpy as np, math
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; pdf=T.pdf
seeds,Yall=T.load_nets(); W=T.regen(seeds[0]); Y=Yall[0].astype(np.float64)
HS=np.load('../pass2_stats.npz')['HS']
srcs={}
for l in range(L):
    mu=(Y[l-1]@W[l].astype(np.float64)) if l>0 else np.zeros(n)
    Sz=W[l].T.astype(np.float64)@HS[l]@W[l].astype(np.float64); sg=np.sqrt(np.diag(Sz)); t=mu/sg
    ph=pdf(t); Ph=ndtr(t); m1=mu*Ph+sg*ph; m2=(mu*mu+sg*sg)*Ph+mu*sg*ph; m3=mu*m2+2*sg*sg*m1
    M=Ph[:,None]*Sz.copy(); np.fill_diagonal(M,0)
    srcs[l]=(M, ph/sg, 2*m1*(1-Ph), m3-3*m2*m1+2*m1**3, Ph, mu, Sz, sg, t)
print("Which functional form does the deep two-point tail's mean-direction coefficient actually take?")
print("u_far = K21(ages 2..8) @ muhat, regressed on sigma^p and on the coherent two-term basis.\n")
print("layer |  R2 on s^2  |  R2 on s^3  |  R2 on s^4 | best fitted power p | R2 coherent 2-term | R2 {s^2,s^3}")
for l in range(5,L):
    Wl=W[l].astype(np.float64); M_,p_,q_,k3h_,d1_,mu,Sz,sg,t=srcs[l]
    muh=mu/np.linalg.norm(mu)
    U=Wl.copy(); Kfar=np.zeros((n,n))
    for age,s in enumerate(range(l-1,max(l-9,-1),-1)):
        Ms,ps,qs,k3hs,d1s,_,_,_,_=srcs[s]; Nt=Ms.T@U; UN=U*Nt
        X1=(UN*2)*ps[:,None]+(U*U)*qs[:,None]
        X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(U*U)*k3hs[:,None]
        Kc=X1.T@Nt+X2.T@U; np.fill_diagonal(Kc,0.0)
        if age>=1: Kfar+=Kc
        U=W[s]@(d1s[:,None]*U)
    u=Kfar@muh
    def r2(B):
        c=np.linalg.lstsq(B,u,rcond=None)[0]; return 1-np.sum((u-B@c)**2)/np.sum(u**2)
    rr={}
    for p in (2,3,4): rr[p]=r2(np.column_stack([sg**p]))
    ps_=np.linspace(1.0,5.0,81); best=max(ps_,key=lambda p:r2(np.column_stack([sg**p])))
    coh=np.column_stack([2*(Sz@muh)*muh*0+2*muh*(Sz@muh), muh*0+sg**2*float(muh@muh)])
    coh=np.column_stack([2*muh*(Sz@muh), sg**2])          # the two coherent terms
    r_coh=r2(coh); r_23=r2(np.column_stack([sg**2,sg**3]))
    print(f"{l+1:5d} |   {rr[2]:.4f}    |   {rr[3]:.4f}    |   {rr[4]:.4f}   |        {best:.2f}         |      {r_coh:.4f}        |    {r_23:.4f}")
