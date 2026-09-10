# Defect A: the homogeneity defect of our deterministic closure.
#   E_p = D_m[m] Mhat + 2 D_Sigma[Sigma] Mhat - p Mhat,  p = 1
# Equivalently, is the closure exactly scale-equivariant: (m, Sigma) -> (c m, c^2 Sigma) giving c * Mhat?
import numpy as np
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; pdf=T.pdf
seeds,Yall=T.load_nets(); W=[w.astype(np.float64) for w in T.regen(seeds[0])]
HS=np.load('../pass2_stats.npz')['HS']; Hm=np.load('../pass2_stats.npz')['Hm']
SCHED=np.array([[0.00094,0.00715],[0.00783,0.01340],[0.01134,0.01947],[0.01406,0.02302],[0.01658,0.02679],
 [0.01853,0.03065],[0.02040,0.03478],[0.02217,0.03824],[0.02331,0.04162],[0.02464,0.04615],
 [0.02555,0.04940],[0.02684,0.05096],[0.02781,0.05280],[0.02838,0.05462],[0.02940,0.05873]])
C2P=np.array([0.00396,0.07819,0.14538,0.22329,0.30105,0.38174,0.45643,0.55252,0.63162,0.70963,
              0.80605,0.92641,1.00444,1.09326,1.18779])
def closure(m0, S0, start, kernels=True):
    m=m0.copy(); Sg=S0.copy()
    for l in range(start,L):
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sg)@Wl; sz2=np.diag(Sz).copy(); sz=np.sqrt(sz2); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        Ms=Sz*Ph[:,None]; np.fill_diagonal(Ms,0.0)
        ps=ph/sz; qs=2*Psi*(1-Ph); m3=mu*m2+2*sz2*Psi; k3hs=m3-3*m2*Psi+2*Psi**3
        k3=np.zeros(n); K21=np.zeros((n,n))
        if l>start:
            U=Wl; Nt=Msp.T@U; UN=U*Nt
            k3=((psp[:,None]*3)*UN*Nt+(qsp[:,None]*3)*UN*U+k3hp[:,None]*(U*U*U)).sum(0)
            X1=(UN*2)*psp[:,None]+(U*U)*qsp[:,None]
            X2=(Nt*Nt)*psp[:,None]+(UN*2)*qsp[:,None]+(U*U)*k3hp[:,None]
            K21=X1.T@Nt+X2.T@U; np.fill_diagonal(K21,0.0)
        if kernels:
            i=l-start if l-start<15 else 14
            k3=k3+SCHED[i,0]*t*s3; g4=SCHED[i,1]; K21=K21+np.outer(C2P[i]*s3, mu/np.linalg.norm(mu))
            np.fill_diagonal(K21,0.0)
        else: g4=0.0
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
        C[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
        Msp,psp,qsp,k3hp=Ms,ps,qs,k3hs
        m=m_next; Sg=C
    return m
START=8
m0=Hm[START].copy(); S0=HS[START].copy()+1e-3*np.mean(np.diag(HS[START]))*np.eye(n)
base=closure(m0,S0,START)
print(f"Defect A: homogeneity of the deterministic closure, layers {START}..15, p = 1\n")
print(f"{'c':>8} | ||Mhat(c m, c^2 Sigma)|| / (c ||Mhat(m,Sigma)||) - 1")
for c in (0.5,0.8,1.25,2.0,4.0):
    out=closure(c*m0,c*c*S0,START)
    print(f"{c:8.2f} | {np.linalg.norm(out)/(c*np.linalg.norm(base))-1:+.3e}   max elementwise rel dev "
          f"{np.abs(out/(c*base)-1).max():.3e}")
print("\n(kernels on; the shipped schedules are degree-correct by construction: a_l*t*sigma^3 and c_l*sigma^3 both scale as c^3)")
