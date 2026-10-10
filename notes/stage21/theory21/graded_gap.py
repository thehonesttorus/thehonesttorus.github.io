import numpy as np
f=lambda r:(np.sqrt(1-r*r)+(np.pi-np.arccos(r))*r)/np.pi
rho=[0.0]
for l in range(16): rho.append(f(rho[-1]))
g=[0.5+np.arcsin(rho[l-1])/np.pi for l in range(1,17)]   # g_l = f'(rho_{l-1}), layers 1..16
print("layer  rho_{l-1}   g_l=f'(rho)")
for l in range(1,17): print(f"{l:3d}   {rho[l-1]:.3f}     {g[l-1]:.3f}")
# amplitude suppression of an incoherent level-k perturbation born at layer b, read at layer L=16:
# variance factor prod_{b<l<=16} g_l^k  (level k), amplitude = sqrt
L=16
print("\nlocalisation: number of most recent layers needed so that older level-k contributions are suppressed (amplitude) below eta")
for k in [1,2,3,4]:
    row=[]
    for eta in [0.1,0.05,0.02,0.01]:
        amp=1.0; ell=0
        for l in range(L,0,-1):
            amp*=np.sqrt(g[l-1]**k); ell+=1
            if amp<eta: break
        row.append(ell if amp<eta else '>16')
    print(f" k={k}: eta=0.1:{row[0]}  0.05:{row[1]}  0.02:{row[2]}  0.01:{row[3]}")
# Monte Carlo check of level-k contraction E||(A W^T)^{(x)k} v^{(x)k}||^2 = (E||A W^T v||^{2k}) vs g^k, n=512, rho=0.8
n=512; rng=np.random.default_rng(0); r0=0.8
for k in [1,2,3]:
    vals=[]
    for t in range(200):
        W=rng.normal(0,np.sqrt(2/n),(n,n)); mu=rng.normal(size=n); mu/=np.linalg.norm(mu)
        m=W.T@mu*np.sqrt(r0/(1-r0)) ; sig=np.sqrt(2/n*n)          # alpha_j = m_j/sigma with Var(alpha)=rho/(1-rho)
        from scipy.stats import norm
        a=norm.cdf(m/np.sqrt(2/n*n)*np.sqrt(n/2)) if False else norm.cdf(rng.normal(0,np.sqrt(r0/(1-r0)),n))
        v=rng.normal(size=n); v-=v@mu*mu; v/=np.linalg.norm(v)
        vals.append(np.linalg.norm(a*(W.T@v))**(2*k))
    gg=0.5+np.arcsin(r0)/np.pi
    print(f"MC n={n}, rho=0.8, k={k}:  E||AW^T v||^(2k) = {np.mean(vals):.4f}   g^k = {gg**k:.4f}")
