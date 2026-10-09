import numpy as np
rng=np.random.default_rng(3); N=20_000_000
A=np.array([[1,.3,.2],[.1,1,.4],[.2,-.3,1.]])
for eps in (0.05,0.025):
    Z=rng.standard_normal((N,3))
    X=Z@A.T + eps*np.stack([Z[:,0]*Z[:,1], Z[:,1]*Z[:,2]-0.3*(Z[:,0]**2-1), Z[:,2]*Z[:,0]+0.5*(Z[:,1]**2-1)],1)
    X=(X-X.mean(0))/X.std(0)                      # standardized (h tilde)
    R=np.corrcoef(X.T)
    k3=lambda i,j,k: np.mean(X[:,i]*X[:,j]*X[:,k])
    c,d,e=0,1,2
    H2c,H2d=X[:,c]**2-1,X[:,d]**2-1
    lhs=np.mean((H2c-H2c.mean())*(H2d-H2d.mean())*X[:,e])   # third cumulant (all centered)
    rhs=4*R[c,d]*k3(c,d,e)+2*R[c,e]*k3(c,d,d)+2*R[d,e]*k3(c,c,d)
    print(f"eps {eps}: kappa(He2_c,He2_d,h_e) MC {lhs:+.5f}  formula {rhs:+.5f}  (gate-cov term {4*R[c,d]*k3(c,d,e):+.5f}); MC se ~{np.std((H2c)*(H2d)*X[:,e])/np.sqrt(N):.5f}")
