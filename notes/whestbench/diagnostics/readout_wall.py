# Where does the remaining 45x live -- in our model of the cumulants, or in the Edgeworth readout?
#
# (A) READOUT WALL.  E[relu(z)] = Psi(mu,sigma) + sigma * sum_{k>=3} c_k phi(t) He_{k-2}(-t),  c_k=E[He_k(x)]/k!
#     exactly (this identity is derived, not fitted: int_{-t}^inf (x+t)phi He_k = phi(t)He_{k-2}(-t)).
#     Estimate c_k empirically at the SAMPLE moments so the O(1/sqrt N) error in mu cancels between the
#     two sides, and watch the truncation converge.  Our estimator stops at the k=4 / k3^2 (k=6) order.
# (B) CUMULANT ORACLE.  Substitute the MC-measured k3 (and k4) at every layer of the deployed closure
#     and read off the final MSE.
import numpy as np, sys, math, os
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=T.load_nets()
SCHED=np.array([[0.00094,0.00715],[0.00783,0.01340],[0.01134,0.01947],[0.01406,0.02302],[0.01658,0.02679],
 [0.01853,0.03065],[0.02040,0.03478],[0.02217,0.03824],[0.02331,0.04162],[0.02464,0.04615],
 [0.02555,0.04940],[0.02684,0.05096],[0.02781,0.05280],[0.02838,0.05462],[0.02940,0.05873]])
C2P=np.array([0.00396,0.07819,0.14538,0.22329,0.30105,0.38174,0.45643,0.55252,0.63162,0.70963,
              0.80605,0.92641,1.00444,1.09326,1.18779])
NB=int(sys.argv[1]) if len(sys.argv)>1 else 200000
KMAX=12

def he_all(x,K):
    """He_0..He_K evaluated at x (probabilists')."""
    H=[np.ones_like(x),x]
    for k in range(2,K+1): H.append(x*H[k-1]-(k-1)*H[k-2])
    return H

def mc_stats(W,NB,seed):
    """Sample means of He_k(x_std) per layer plus mean(relu(z)), accumulated in chunks."""
    rb=np.random.default_rng(seed); ch=4000; done=0
    S1=np.zeros((L,n)); S2=np.zeros((L,n)); SR=np.zeros((L,n)); PW=np.zeros((L,KMAX+1,n))
    Zc=[]  # we need centred/standardised moments, so do two passes: pass 1 for mu,sigma
    while done<NB:
        b=min(ch,NB-done); Xb=rb.standard_normal((b,n)).astype(np.float32); hb=Xb
        for l in range(L):
            z=(hb@W[l]).astype(np.float64); S1[l]+=z.sum(0); S2[l]+=(z*z).sum(0)
            SR[l]+=np.maximum(z,0).sum(0); hb=np.maximum(z,0).astype(np.float32)
        done+=b
    mu=S1/NB; v=S2/NB-mu*mu; sd=np.sqrt(np.maximum(v,1e-30)); relu=SR/NB
    rb=np.random.default_rng(seed); done=0
    while done<NB:
        b=min(ch,NB-done); Xb=rb.standard_normal((b,n)).astype(np.float32); hb=Xb
        for l in range(L):
            z=(hb@W[l]).astype(np.float64); x=(z-mu[l])/sd[l]; H=he_all(x,KMAX)
            for k in range(3,KMAX+1): PW[l,k]+=H[k].sum(0)
            hb=np.maximum(z,0).astype(np.float32)
        done+=b
    return mu,sd,relu,PW/NB

print(f"N = {NB} samples per network.\n")
print("(A) READOUT: does the exact Hermite/Gram-Charlier series close the gap?")
print("    residual r = mean(relu(z)) - Psi(mu_hat,sigma_hat), both at the SAME sample moments.")
print("    'rel' = || r - series_K || / || r ||;  'mse' = mean((Psi+series_K - Y_true)^2) at that layer.")
hdr=f"{'net':>3} {'layer':>5} {'||r||':>10} " + " ".join(f"K={k:<2d}".rjust(9) for k in [3,4,6,8,10,12])
print(hdr)
NETS=[0,1]
CACHE={}
for k in NETS:
    fn=f"rw_{k}_{NB}.npz"
    if os.path.exists(fn):
        d=np.load(fn); mu,sd,relu,Hm=d['mu'],d['sd'],d['relu'],d['Hm']
    else:
        W=T.regen(seeds[k]); mu,sd,relu,Hm=mc_stats(W,NB,7000+k)
        np.savez_compressed(fn,mu=mu,sd=sd,relu=relu,Hm=Hm)
    CACHE[k]=(mu,sd,relu,Hm)
    Y=Yall[k]
    for l in [4,8,12,15]:
        t=mu[l]/sd[l]; ph=pdf(t); Ph=ndtr(t); Psi=mu[l]*Ph+sd[l]*ph
        r=relu[l]-Psi; Hm2=he_all(-t,KMAX)   # He_{k-2}(-t)
        row=[]; ser=np.zeros(n); fact=1.0
        out={}
        for kk in range(3,KMAX+1):
            fact=float(math.factorial(kk))
            ser=ser+sd[l]*(Hm[l,kk]/fact)*ph*Hm2[kk-2]
            if kk in (3,4,6,8,10,12): out[kk]=(np.linalg.norm(r-ser)/np.linalg.norm(r), np.mean((Psi+ser-Y[l])**2))
        print(f"{k:3d} {l:5d} {np.linalg.norm(r):10.3e} " + " ".join(f"{out[kk][0]:9.4f}" for kk in [3,4,6,8,10,12]))
        print(f"{'':3s} {'mse':>5s} {np.mean((Psi-Y[l])**2):10.3e} " + " ".join(f"{out[kk][1]:9.3e}" for kk in [3,4,6,8,10,12]))
