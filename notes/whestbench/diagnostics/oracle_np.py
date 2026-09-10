# Decompose pipeline error: oracle covariance (sample cov at every layer) vs oracle mean (truth) with propagated cov
import sys, time, numpy as np
sys.path.insert(0,'.')
import estimator as E
from types import SimpleNamespace
W=np.load('../data/W0.npy'); Y=np.load('../data/Y0.npy').astype(np.float64)
S=np.load('../pass2_stats.npz'); HS=S['HS']
mlp=SimpleNamespace(width=1024,depth=16,weights=[W[l] for l in range(16)],seed=0,name='net0')
class Twin(E.CumulantClosureEstimator):
    oracle_cov=False; oracle_mean=False; cov_scale=0.0
    def predict(self, mlp, budget):
        # copy of the loop with hooks: swap Sig/m by oracle values after each layer
        xp=self.xp; n,L=mlp.width,mlp.depth; F32=np.float32
        Ws=[xp.asarray(w) for w in mlp.weights]; K=self.K_SOURCES; outs=[]; sources={}
        W0=Ws[0]; G=W0.T@W0; s2=xp.diag(G); s=xp.sqrt(s2); m=s*F32(E.INV_SQRT2PI); so=xp.outer(s,s)
        rho=xp.clip(G/so,-1,1); arc=xp.arccos(rho)
        C=so*((xp.sqrt(1-rho*rho)+rho*(np.pi-arc))*F32(1/(2*np.pi)))-xp.outer(m,m)
        var=s2*F32(0.5-1/(2*np.pi)); C=C-xp.diag(xp.diag(C))+xp.diag(var)
        d1=xp.zeros(n,dtype=xp.float32)+F32(0.5); M=G*F32(0.5); M=M-xp.diag(xp.diag(M))
        sources[0]=(M,F32(E.INV_SQRT2PI)/s,m,s2*s*F32(2*E.INV_SQRT2PI-1.5*E.INV_SQRT2PI+2*E.INV_SQRT2PI**3),d1)
        outs.append(m); Sig=C
        for l in range(1,L):
            if self.oracle_cov: Sig=HS[l].astype(np.float32)
            if self.oracle_mean: m=Y[l-1].astype(np.float32)
            Wl=Ws[l]; mu=m@Wl; Sz=self.sandwich(Wl,Sig)*(1+self.cov_scale); sz2=xp.diag(Sz); sz=xp.sqrt(sz2); t=mu/sz
            ph=E._phi(xp,t); Ph=E._Phi(xp,t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph
            k3=None; act=np.nonzero(np.abs(t)<self.T_ACT)[0]; s_lo=max(l-K,0)
            if act.size>0 and K>0:
                U=Wl[:,act]; k3a=np.zeros(act.size,dtype=np.float32)
                for src in range(l-1,s_lo-1,-1):
                    Ms,ps,qs,k3hs,d1s=sources[src]; Nt=Ms.T@U; UN=U*Nt
                    k3a+=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
                    if src>s_lo: U=Ws[src]@(d1s[:,None]*U)
                k3=np.zeros(n,dtype=np.float32); k3[act]=k3a
            m_next=Psi
            if k3 is not None:
                tph=t*ph; m_next=m_next-k3*tph/(sz2*6)
                if self.USE_K3SQ: t2=t*t; m_next=m_next+(k3*k3)*((t2-6)*t2+3)*ph/(sz2*sz2*sz*72)
            so=xp.outer(sz,sz); rho=Sz/so; ds=self.hermite_coeffs(t,ph,Ph,self.HERMITE_ORDER); rk=rho; C=None; fact=1.0
            for k in range(1,self.HERMITE_ORDER+1):
                fact*=k; dk=ds[k-1]; term=rk*xp.outer(dk,dk*F32(1/fact)); C=term if C is None else C+term
                if k<self.HERMITE_ORDER: rk=rk*rho
            C=C*so; var=m2-Psi*Psi
            if k3 is not None and self.VAR_K3: var=var+k3*ph/(sz*3)-2*Psi*(m_next-Psi)
            C=C-xp.diag(xp.diag(C))+xp.diag(var)
            if K>0:
                M=Sz*Ph[:,None]; M=M-xp.diag(xp.diag(M)); sources[l]=(M,ph/sz,2*Psi*(1-Ph),self._k3_relu(xp,mu,sz,sz2,Ph,ph,Psi,m2),Ph)
            m=m_next; Sig=C; outs.append(m)
        return xp.stack(outs)
def run(**kw):
    est=Twin(xp=np, symm=lambda C:C)
    for k,v in kw.items(): setattr(est,k,v)
    P=est.predict(mlp,0); return np.mean((P.astype(np.float64)-Y)**2,axis=1)
for kw in [dict(K_SOURCES=0),dict(K_SOURCES=15),dict(K_SOURCES=0,oracle_cov=True),dict(K_SOURCES=15,oracle_cov=True),dict(K_SOURCES=0,oracle_mean=True),dict(K_SOURCES=15,oracle_mean=True),
           dict(K_SOURCES=15,cov_scale=0.002),dict(K_SOURCES=15,cov_scale=0.004),dict(K_SOURCES=15,cov_scale=0.006),dict(K_SOURCES=0,VAR_K3=False),dict(K_SOURCES=15,VAR_K3=False),dict(K_SOURCES=15,USE_K3SQ=False)]:
    mse=run(**kw); print(kw, f"final={mse[-1]:.2e}", " ".join(f"{x:.1e}" for x in mse[1:]), flush=True)
