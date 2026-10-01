import sys; sys.path.insert(0,'../../bench'); sys.path.insert(0,'.'); import bench, bethe, oracle, numpy as np, warnings
warnings.filterwarnings('ignore')
S=bench.load_set(sys.argv[1]); i=int(sys.argv[2]); W=bench.weights(S,i).astype(np.float64); n=W.shape[1]
tr=oracle.mc_state(W,int(float(sys.argv[3])))
def chain(mode):
    m=np.zeros(n); C=W[0].T@W[0]; K=np.zeros((n,n)); k4=np.zeros(n); Qz=np.zeros((n,n)); prev=None; Cs=[]
    for l in range(16):
        Cs.append(C)
        if mode=='edge': mu,Ca,Ka,k3a,k4a,c,L0,Lm1=bethe.relu_map_edges(m,C,K,k4); Q=None
        else:
            mu,Ca,Ka,k3a,k4a,c,L,Q,R=bethe.relu_map_v4(m,C,K,k4,Qz,need4=l<15); L0,Lm1=L[0],L[-1]
        if l==15: break
        Wn=W[l+1]; W2=Wn*Wn
        M=c*c*(L0**2)[:,None]*Lm1[None,:]; spec=(k3a,Ka-M,c,L0,Lm1)
        Kn=bethe.contract(*spec,Wn,Wn)
        if prev is not None:
            pspec,Wl=prev; P=Wl@(L0[:,None]*Wn); Kz=K.copy(); k3z=np.diag(K).copy(); np.fill_diagonal(Kz,0)
            Kn=Kn+bethe.contract(*pspec,P,P)-bethe.contract(k3z*L0**3,Kz*(L0**2)[:,None]*L0[None,:],None,L0,Lm1,Wn,Wn)
        prev=(spec,Wn)
        if mode!='edge':
            Qn=W2.T@((Q+np.diag(k4a))@W2); X=L0[:,None]*Wn; U=c@X; s=(c*c)@(X*X); h4=2*(L0*(1-L0)-mu*Lm1); H=(W2*h4[:,None]).T@(U*U-s)
            Qn=Qn+H+H.T; np.fill_diagonal(Qn,0); k4=bethe.k4_next(Wn,k4a,Q,R,c,L,mu); Qz=Qn
        else: k4=(Wn**4).T@k4a
        m=mu@Wn; C=Wn.T@Ca@Wn; K=Kn
    return Cs
off=~np.eye(n,dtype=bool)
for mode in ('edge','v4'):
    Cs=chain(mode)
    print(mode)
    for l in (1,3,6,9,12,15):
        Ct=tr[l]['C']; Co=Cs[l]; d=Co-Ct
        et=np.linalg.eigvalsh(Ct)[-1]; eo=np.linalg.eigvalsh(Co)[-1]
        u=np.linalg.eigh(Ct)[1][:,-1]
        print(f"  l{l}: rel err v {np.sqrt(np.mean(np.diag(d)**2)/np.mean(np.diag(Ct)**2)):.4f} offdiag {np.sqrt(np.mean(d[off]**2)/np.mean(Ct[off]**2)):.4f} | top eig own {eo:.3f} true {et:.3f} | err along top dir {u@d@u:.4f}",flush=True)
