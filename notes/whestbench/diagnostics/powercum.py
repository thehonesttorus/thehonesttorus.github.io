# Verify ARC's power-cumulant identity for the two-point third cumulant, the channel that dominates
# our error, and price it.
#   kappa_3[phi(Z)]_{iij} = Cov(phi(Z_i)^2, phi(Z_j)) - 2 E[phi(Z_i)] Cov(phi(Z_i), phi(Z_j))
# Both terms are Hadamard power series in the correlation rho_ij:
#   Cov(f(Z_i), g(Z_j)) = sum_{k>=1} rho^k * bf_k(t_i) * bg_k(t_j) / k!    (scaled by sigma powers)
# with bf_k the Hermite coefficients of f in the standardised variable.
import numpy as np
from numpy.polynomial.hermite_e import hermeval
from scipy.special import ndtr
SQ2PI=np.sqrt(2*np.pi); pdf=lambda x: np.exp(-0.5*x*x)/SQ2PI
NQ=4001
def herm_coeffs(t, kmax, power):
    """b_k(t) = E[ ((t+zeta)_+)^power He_k(zeta) ] by Gauss-Hermite-style quadrature on a fine grid."""
    g=np.linspace(-12,12,NQ); w=pdf(g)*(g[1]-g[0])
    f=np.maximum(t[:,None]+g[None,:],0.0)**power                       # (n, NQ)
    out=np.empty((kmax+1,t.size))
    for k in range(kmax+1):
        c=np.zeros(k+1); c[k]=1.0
        out[k]=(f*(hermeval(g,c)*w)[None,:]).sum(1)
    return out
def analytic_relu_coeffs(t,kmax):
    """closed form for ReLU: b_0=t*Phi+phi, b_1=Phi, b_k=phi*He_{k-2}(-t) for k>=2"""
    ph=pdf(t); Ph=ndtr(t); out=np.empty((kmax+1,t.size)); out[0]=t*Ph+ph; out[1]=Ph
    x=-t; hp=np.ones_like(x); h=x
    for k in range(2,kmax+1):
        out[k]=ph*(hp if k==2 else h)
        if k>2: hp,h=h,x*h-(k-2)*hp
    # recompute cleanly: He_{k-2}(-t)
    hs=[np.ones_like(x),x]
    for m in range(2,kmax): hs.append(x*hs[-1]-(m-1)*hs[-2])
    for k in range(2,kmax+1): out[k]=ph*hs[k-2]
    return out
KMAX=12
tt=np.array([-1.5,-0.5,0.0,0.4,1.0,2.0])
bq=herm_coeffs(tt,KMAX,1); ba=analytic_relu_coeffs(tt,KMAX)
print("ReLU Hermite coefficients: quadrature vs closed form, max abs diff =", np.abs(bq-ba).max())
cq=herm_coeffs(tt,KMAX,2)
print("ReLU^2 coefficients (quadrature), k=0..4 at t=0.4:", np.round(cq[:5,3],5))
# --- validate the series and the identity against Monte Carlo on bivariate Gaussians ---
rng=np.random.default_rng(0); N=4_000_000
print("\n ti     tj    rho | Cov(A,B): series/MC | Cov(A^2,B): series/MC | K21: series/MC")
for (ti,tj,rho) in [(0.4,-0.3,0.30),(1.0,0.5,0.50),(-0.5,0.8,0.20),(0.2,0.2,0.60),(1.5,-1.0,0.45)]:
    z1=rng.standard_normal(N); z2=rho*z1+np.sqrt(1-rho*rho)*rng.standard_normal(N)
    A=np.maximum(ti+z1,0); B=np.maximum(tj+z2,0)
    t2=np.array([ti,tj]); b=analytic_relu_coeffs(t2,KMAX); c=herm_coeffs(t2,KMAX,2)
    ks=np.arange(1,KMAX+1); fac=np.array([np.math.factorial(int(k)) if hasattr(np,'math') else float(__import__('math').factorial(int(k))) for k in ks])
    covAB=np.sum(rho**ks*b[1:,0]*b[1:,1]/fac)
    covA2B=np.sum(rho**ks*c[1:,0]*b[1:,1]/fac)
    EA=b[0,0]
    k21=covA2B-2*EA*covAB
    mc_ab=np.cov(A,B,bias=True)[0,1]; mc_a2b=np.cov(A*A,B,bias=True)[0,1]
    mc_k21=np.mean((A-A.mean())**2*(B-B.mean()))
    print(f"{ti:+.1f} {tj:+.1f} {rho:.2f} | {covAB:+.5f}/{mc_ab:+.5f} | {covA2B:+.5f}/{mc_a2b:+.5f} | {k21:+.6f}/{mc_k21:+.6f}")
