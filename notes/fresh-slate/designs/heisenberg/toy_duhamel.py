"""Toy exactness check of the Heisenberg-Duhamel identity and of the kappa_3 Stein pairing.
term_l = E_{rho_l}[g_l] - E_{nu_l}[g_l], rho_l = T nu_{l-1}; sum_l term_l = truth - closure."""
import numpy as np, sys
sys.path.insert(0,'/home/user/thehonesttorus.github.io/notes/fresh-slate/designs/heisenberg')
from hd import *
n=int(sys.argv[1]); L=int(sys.argv[2]); N=int(float(sys.argv[3])); seed=int(sys.argv[4]) if len(sys.argv)>4 else 0
rng=np.random.default_rng(seed); Ws=rng.normal(size=(L,n,n))*np.sqrt(2/n)
# reference chain (m_l, C_l) of z_l, l=0..L-1 (z_0 = x W_1)
ms=[np.zeros(n)]; Cs=[Ws[0].T@Ws[0]]
for l in range(L-1):
    G=Gauss(ms[-1],Cs[-1],K=2); Ca=G.cov_a(); ms.append(G.Ea@Ws[l+1]); Cs.append(Ws[l+1].T@Ca@Ws[l+1])
Gf=Gauss(ms[-1],Cs[-1],K=2); closure_out=Gf.Ea
def run_from(z,l):   # z = samples of z_l, push to final a
    a=np.maximum(z,0)
    for k in range(l+1,L): a=np.maximum(a@Ws[k],0)
    return a
def E_from_gauss(l, push_one, rng, chunk=1_000_000):
    """push_one=False: E_{nu_l}[g_l];  True: E_{T nu_{l-1}}[g_l] (sample z_{l-1} ~ nu_{l-1})."""
    ll = l-1 if push_one else l
    Lc=np.linalg.cholesky(Cs[ll]+1e-12*np.eye(n)); S=np.zeros(n); done=0
    while done<N:
        b=min(chunk,N-done); z=ms[ll]+rng.standard_normal((b,n))@Lc.T
        S+=run_from(z,ll).sum(0); done+=b
    return S/N
truth=E_from_gauss(0,False,np.random.default_rng(seed+7))
terms=[]
for l in range(1,L):
    r=np.random.default_rng(1000+l)   # common random numbers are not used: independent estimates
    terms.append(E_from_gauss(l,True,np.random.default_rng(2000+l))-E_from_gauss(l,False,np.random.default_rng(3000+l)))
terms=np.array(terms)
print("rms(truth-closure) %.4e  rms(truth-closure-sum terms) %.2e   MC se ~ %.1e"%(np.sqrt(((truth-closure_out)**2).mean()),
      np.sqrt(((truth-closure_out-terms.sum(0))**2).mean()), 0.6/np.sqrt(N)))
cl=closure(Ws)[-1]
print("closure check", np.abs(cl-closure_out).max())
for l in range(1,L):
    for second in [False,True]:
        h=hd(Ws,A=None,only_src=l-1,second=second)[-1]-cl
        e=terms[l-1]
        print("term l=%d  rms exact %.3e  rms HD(k3%s) %.3e  rms resid %.3e  corr %.3f"%(l,np.sqrt((e**2).mean()),"+2nd" if second else "",np.sqrt((h**2).mean()),np.sqrt(((e-h)**2).mean()),np.corrcoef(e,h)[0,1]))
