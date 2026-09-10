# The diagonal channel's deep tail turned out to be rank one in {sigma^3, t sigma^3} with a universal
# depth schedule, which made it shippable at zero cost.  The two-point channel at depth is worth 3.7x
# but costs ~49% of budget.  Ask the same question of it: is the AGGREGATE deep two-point contribution
# low-dimensional and analytic?
#   K21_deep = sum over source ages 2..8 of the two-point tensor kappa(z_a, z_a, z_b)
#   its effect enters as  dC = 0.5*(K21 (x) (F2, Phi) + transpose),  F2 = phi(t)/sigma
# Measured with oracle inputs (1e9-sample truth means, high-sample covariances) so the target is clean.
import numpy as np, sys, time
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; pdf=T.pdf
seeds,Yall=T.load_nets()
def analyse(netfile, k, HS, AGE_NEAR=1, AGE_FAR=8):
    W=T.regen(seeds[k]); Y=Yall[k].astype(np.float64)
    srcs={}
    for l in range(L):
        mu=(Y[l-1]@W[l].astype(np.float64)) if l>0 else np.zeros(n)
        Sz=W[l].T.astype(np.float64)@HS[l]@W[l].astype(np.float64); sg=np.sqrt(np.diag(Sz)); t=mu/sg
        ph=pdf(t); Ph=ndtr(t); m1=mu*Ph+sg*ph; m2=(mu*mu+sg*sg)*Ph+mu*sg*ph; m3=mu*m2+2*sg*sg*m1
        M=Ph[:,None]*Sz.copy(); np.fill_diagonal(M,0)
        srcs[l]=(M, ph/sg, 2*m1*(1-Ph), m3-3*m2*m1+2*m1**3, Ph)
    rows=[]
    for l in range(5,L):
        Wl=W[l].astype(np.float64); mu=Y[l-1]@Wl
        Sz=W[l].T.astype(np.float64)@HS[l]@Wl; sg=np.sqrt(np.diag(Sz)); t=mu/sg; s3=sg**3
        muh=mu/np.linalg.norm(mu)
        U=Wl.copy(); K_near=np.zeros((n,n)); K_far=np.zeros((n,n))
        for age,s in enumerate(range(l-1,max(l-1-AGE_FAR,-1),-1)):
            Ms,ps,qs,k3hs,d1s=srcs[s]; Nt=Ms.T@U; UN=U*Nt
            X1=(UN*2)*ps[:,None]+(U*U)*qs[:,None]
            X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(U*U)*k3hs[:,None]
            Kc=X1.T@Nt+X2.T@U; np.fill_diagonal(Kc,0.0)
            if age<AGE_NEAR: K_near+=Kc
            else: K_far+=Kc
            U=W[s]@(d1s[:,None]*U)
        # how much of the far tail's ACTION is the rank-one mean-direction form?
        u_far=K_far@muh; R1=np.outer(u_far,muh); np.fill_diagonal(R1,0.0)
        F2=ph_=pdf(t)/sg; Ph=ndtr(t)
        def act(K):
            D=0.5*(K*np.outer(F2,Ph)+K.T*np.outer(Ph,F2)); np.fill_diagonal(D,0.0)
            Wn=W[l+1].astype(np.float64) if l+1<L else W[l].astype(np.float64)
            return np.einsum('ij,ik,kj->j',Wn,D,Wn)
        a_far=act(K_far); a_r1=act(R1)
        r1_share=1-np.sum((a_far-a_r1)**2)/np.sum(a_far**2)
        # is u_far analytic?  regress on sigma^3 * {1, t, t^2, t^3}
        def r2(B):
            c=np.linalg.lstsq(B,u_far,rcond=None)[0]; return 1-np.sum((u_far-B@c)**2)/np.sum(u_far**2), c
        b1,_=r2(np.column_stack([s3]))
        b2,c2=r2(np.column_stack([s3,t*s3]))
        b3,c3=r2(np.column_stack([s3,t*s3,t*t*s3]))
        b4,c4=r2(np.column_stack([s3,t*s3,t*t*s3,t**3*s3]))
        rows.append((l+1, np.linalg.norm(K_far)/np.linalg.norm(K_near+K_far), r1_share, b1,b2,b3,b4, c3[0],c3[1],c3[2]))
    return rows
HS0=np.load('../pass2_stats.npz')['HS']
print("Two-point channel: is the deep tail's aggregate low-dimensional and analytic?  (network 0, oracle inputs)")
print("far share = ||K21 from ages 2..8|| / ||K21 from ages 1..8||")
print("rank1     = share of the far tail's induced next-layer variance captured by the mean-direction rank-one form")
print("R2 of u_far = K21_far @ muhat regressed on sigma^3 * {1, t, t^2, t^3}\n")
print("layer | far share | rank1 | R2:{s3} {s3,ts3} {+t2s3} {+t3s3} | fitted c0 c1 c2")
for r in analyse(None,0,HS0):
    print(f"{r[0]:5d} |   {r[1]:.3f}   | {r[2]:+.3f} | {r[3]:.3f} {r[4]:.3f} {r[5]:.3f} {r[6]:.3f} | {r[7]:+.4f} {r[8]:+.4f} {r[9]:+.4f}")
