# At layer 1 EVERYTHING except the cumulants is exact: mu_1 = Y_0 W_1 and Sigma_0 is the closed-form
# arcsine covariance of relu(W_0 x).  So the residual there is a pure readout question.  Build the
# Gram-Charlier series with coefficients measured at the EXACT (mu, sigma) and watch it converge on
# the 1e9-sample ground truth.
import numpy as np, sys, os, math
from scipy.special import ndtr
import ablate as AB
n=1024; L=16; INV=1/np.sqrt(2*np.pi); pdf=lambda t: np.exp(-0.5*t*t)*INV
NB=int(sys.argv[1]) if len(sys.argv)>1 else 1000000
KM=10
def he(x,K):
    H=[np.ones_like(x),x]
    for k in range(2,K+1): H.append(x*H[k-1]-(k-1)*H[k-2])
    return H
for k in [0,1]:
    W=AB.getW(k); Y=AB.Yall[k]
    G=(W[0].T@W[0]).astype(np.float64); s2=np.diag(G); s=np.sqrt(s2); so=np.outer(s,s)
    rho=np.clip(G/so,-1,1)
    C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(s*INV,s*INV)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi))
    W1=W[1].astype(np.float64)
    mu=Y[0]@W1                                 # exact
    sz2=np.diag((W1.T@C)@W1); sz=np.sqrt(sz2)  # exact
    t=mu/sz; ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph
    fn=f"l1_{k}_{NB}.npz"
    if os.path.exists(fn): Hm=np.load(fn)['Hm']
    else:
        rb=np.random.default_rng(4242+k); ch=5000; done=0; Hm=np.zeros((KM+1,n))
        while done<NB:
            b=min(ch,NB-done); X=rb.standard_normal((b,n)).astype(np.float32)
            z=(np.maximum(X@W[0],0)@W[1]).astype(np.float64)
            x=(z-mu)/sz; H=he(x,KM)
            for kk in range(3,KM+1): Hm[kk]+=H[kk].sum(0)
            done+=b
        Hm/=NB; np.savez_compressed(fn,Hm=Hm)
    H2=he(-t,KM)
    print(f"\nnet {k}: layer 1, N={NB}.  Y-noise floor about {np.mean((2.6e-5)**2):.1e} MSE.")
    print(f"  Gaussian readout            MSE {np.mean((Psi-Y[1])**2):.4e}")
    ser=np.zeros(n)
    for kk in range(3,KM+1):
        ser=ser+sz*(Hm[kk]/float(math.factorial(kk)))*ph*H2[kk-2]
        lab=f"  + Hermite order {kk:2d}"
        print(f"{lab:28s} MSE {np.mean((Psi+ser-Y[1])**2):.4e}   "
              f"lambda_{kk} = {Hm[kk].std():.4f} (rms over neurons)")
    # what our estimator actually uses: kappa3, kappa4 and the kappa3^2 piece of order 6
    l3=Hm[3]; l4=Hm[4]
    d=sz*ph*(l3/6*H2[1]+l4/24*H2[2]+(l3*l3)/72*H2[4])
    print(f"  deployed readout (k3,k4,k3^2)  MSE {np.mean((Psi+d-Y[1])**2):.4e}")
