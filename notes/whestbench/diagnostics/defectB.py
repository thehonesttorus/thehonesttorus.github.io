# Defect B, the Gaussian differentiation (Price) relation.  For the TRUE functional Price's theorem
# gives, exactly,   D_Sigma[G] M = (1/2) G : D^2_m M.
# With G = v v^T this is   D_Sigma[v v^T] M = (1/2) D^2_m[v,v] M.
# It is the relation linking a DEGREE-TWO perturbation to a product of two DEGREE-ONE derivatives,
# i.e. the Leibniz/multiplicativity relation across the grading.  A separately-fitted moment closure
# need not satisfy it.  Measured here for our closure by central finite differences, and for the true
# network by Monte Carlo with common random numbers as a control.
import numpy as np, time
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
def closure(m0,S0,start,kern=True):
    """deterministic closure from layer `start` to the end; returns the final-layer mean"""
    m=m0.copy(); Sg=S0.copy(); Msp=psp=qsp=k3hp=None
    for l in range(start,L):
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sg)@Wl; sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        k3=np.zeros(n); K21=np.zeros((n,n))
        if Msp is not None:
            U=Wl; Nt=Msp.T@U; UN=U*Nt
            k3=((psp[:,None]*3)*UN*Nt+(qsp[:,None]*3)*UN*U+k3hp[:,None]*(U*U*U)).sum(0)
            X1=(UN*2)*psp[:,None]+(U*U)*qsp[:,None]
            X2=(Nt*Nt)*psp[:,None]+(UN*2)*qsp[:,None]+(U*U)*k3hp[:,None]
            K21=X1.T@Nt+X2.T@U; np.fill_diagonal(K21,0.0)
        i=min(l-start,14); g4=0.0
        if kern:
            k3=k3+SCHED[i,0]*t*s3; g4=SCHED[i,1]
            nm=np.linalg.norm(mu); muh=mu/nm if nm>0 else mu
            K21=K21+np.outer(C2P[i]*s3,muh); np.fill_diagonal(K21,0.0)
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
        Msrc=Sz*Ph[:,None]; np.fill_diagonal(Msrc,0.0)
        m3=mu*m2+2*sz2*Psi
        Msp,psp,qsp,k3hp=Msrc,ph/sz,2*Psi*(1-Ph),m3-3*m2*Psi+2*Psi**3
        m=m_next; Sg=C
    return m
def psqrt(S):
    ev,V=np.linalg.eigh(S); ev=np.clip(ev,0,None); return (V*np.sqrt(ev))@V.T
NS=200_000; CH=20_000
rng=np.random.default_rng(31415); XI=[rng.standard_normal((CH,n)) for _ in range(NS//CH)]
def truth(mm,SS,start):
    A=psqrt(SS); acc=np.zeros(n)
    for xi in XI:
        H=mm[None,:]+xi@A
        for k in range(start,L): H=np.maximum(H@W[k],0.0)
        acc+=H.mean(0)
    return acc/len(XI)
print("Defect B: the Price relation  D_Sigma[v v^T] M = (1/2) D^2_m[v,v] M.")
print("Exactly zero for the true functional.  Reported relative to ||M||.\n")
print(" start | closure: ||P_v||/||M||  | truth (control): ||P_v||/||M|| | closure ||M-truth||/||M||")
for START in (4,8,11,13):
    m0=Hm[START].copy(); S0=HS[START].copy()+1e-3*np.mean(np.diag(HS[START]))*np.eye(n)
    rs=np.random.default_rng(START); v=rs.standard_normal(n); v/=np.linalg.norm(v)
    G=np.outer(v,v); s=t=0.05
    base=closure(m0,S0,START)
    dS=(closure(m0,S0+t*G,START)-closure(m0,S0-t*G,START))/(2*t)
    d2=(closure(m0+s*v,S0,START)-2*base+closure(m0-s*v,S0,START))/(s*s)
    Pc=dS-0.5*d2
    tb=truth(m0,S0,START)
    tdS=(truth(m0,S0+t*G,START)-truth(m0,S0-t*G,START))/(2*t)
    td2=(truth(m0+s*v,S0,START)-2*tb+truth(m0-s*v,S0,START))/(s*s)
    Pt=tdS-0.5*td2
    print(f" {START:5d} |      {np.linalg.norm(Pc)/np.linalg.norm(base):.3e}       |        {np.linalg.norm(Pt)/np.linalg.norm(tb):.3e}          |      {np.linalg.norm(base-tb)/np.linalg.norm(tb):.3e}", flush=True)
