# Our source construction is the leading-order (one power of rho per edge) diagram truncation.
# ARC's power-cumulant identity gives the SAME object exactly, as a Hadamard power series:
#   K21_ab = kappa(h_a, h_a, h_b) = Cov(h_a^2, h_b) - 2 E[h_a] Cov(h_a, h_b),  h = relu(z), z Gaussian
#   Cov(f(z_a), g(z_b)) = sum_{k>=1} rho^k bf_k(t_a) bg_k(t_b) / k! * (sigma scalings)
# Compare both against a large Monte Carlo on a REAL layer.
import numpy as np, math
from numpy.polynomial.hermite_e import hermeval
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; SQ2PI=np.sqrt(2*np.pi); pdf=lambda x: np.exp(-0.5*x*x)/SQ2PI
seeds,Yall=T.load_nets(); W=T.regen(seeds[0])
def relu_coeffs(t,kmax):
    ph=pdf(t); Ph=ndtr(t); out=np.empty((kmax+1,t.size)); out[0]=t*Ph+ph; out[1]=Ph
    x=-t; hs=[np.ones_like(x),x]
    for m in range(2,kmax): hs.append(x*hs[-1]-(m-1)*hs[-2])
    for k in range(2,kmax+1): out[k]=ph*hs[k-2]
    return out
def relu2_coeffs(t,kmax,NQ=6001):
    g=np.linspace(-14,14,NQ); w=pdf(g)*(g[1]-g[0])
    f=np.maximum(t[:,None]+g[None,:],0.0)**2
    out=np.empty((kmax+1,t.size))
    for k in range(kmax+1):
        c=np.zeros(k+1); c[k]=1.0; out[k]=(f*(hermeval(g,c)*w)[None,:]).sum(1)
    return out
NS=600_000; ch=10_000; rng=np.random.default_rng(11)
for LAYER in (3,7,11):
    # accumulate the true two-point cumulant of h_{LAYER+1} and the moments of z_LAYER
    m1=np.zeros(n); S=np.zeros((n,n)); S2=np.zeros((n,n)); zm=np.zeros(n); zS=np.zeros((n,n)); N=0
    for _ in range(NS//ch):
        h=rng.standard_normal((ch,n)).astype(np.float32)
        for l in range(LAYER): h=np.maximum(h@W[l],0)
        z=(h@W[LAYER]).astype(np.float64); hp=np.maximum(z,0)
        zm+=z.sum(0); zS+=z.T@z; m1+=hp.sum(0); S+=hp.T@hp; S2+=(hp*hp).T@hp; N+=ch
    zm/=N; zS=zS/N-np.outer(zm,zm); m1/=N; Cov=S/N-np.outer(m1,m1); M2=S2/N-np.outer((S/N).diagonal()*0+ (S2*0).diagonal(),np.zeros(n))
    A2=S2/N - np.outer((np.diag(S)/N), m1)          # Cov(h_a^2, h_b)
    K21_true=A2-2*m1[:,None]*Cov
    np.fill_diagonal(K21_true,0.0)
    sg=np.sqrt(np.diag(zS)); t=zm/sg; rho=zS/np.outer(sg,sg); np.fill_diagonal(rho,0.0)
    ph=pdf(t); Ph=ndtr(t); Psi=zm*Ph+sg*ph
    # (a) our leading-order diagram source, evaluated in its own layer:  K21_ab = q_a M_ba + p_b M_ab^2
    M=Ph[:,None]*zS.copy(); np.fill_diagonal(M,0.0)
    p=ph/sg; q=2*Psi*(1-Ph)
    K21_diag=q[:,None]*M.T+p[None,:]*(M*M); np.fill_diagonal(K21_diag,0.0)
    # (b) ARC power-cumulant closed form
    KM=10; b=relu_coeffs(t,KM); c=relu2_coeffs(t,KM)
    covAB=np.zeros((n,n)); covA2B=np.zeros((n,n)); rk=rho.copy()
    for k in range(1,KM+1):
        f=math.factorial(k)
        covAB+=rk*np.outer(b[k],b[k])/f
        covA2B+=rk*np.outer(c[k],b[k])/f
        if k<KM: rk=rk*rho
    covAB*=np.outer(sg,sg); covA2B*=np.outer(sg*sg,sg)
    K21_exact=covA2B-2*Psi[:,None]*covAB; np.fill_diagonal(K21_exact,0.0)
    def rel(X): return np.linalg.norm(X-K21_true)/np.linalg.norm(K21_true)
    print(f"layer {LAYER+1:2d}: ||K21_true||={np.linalg.norm(K21_true):.4f} | "
          f"leading-order diagrams rel err {rel(K21_diag):.4f} | power-cumulant closed form rel err {rel(K21_exact):.4f} | "
          f"corr(diag,true)={np.corrcoef(K21_diag.ravel(),K21_true.ravel())[0,1]:.4f} corr(exact,true)={np.corrcoef(K21_exact.ravel(),K21_true.ravel())[0,1]:.4f}", flush=True)
