import numpy as np
rng=np.random.default_rng(3)
for n in [512,1024]:
    L=16; N=2000
    Ws=[rng.normal(0,np.sqrt(2/n),(n,n)) for _ in range(L-1)]
    X=rng.normal(size=(N,n))
    U=X
    for W in Ws: U=np.maximum(U@W,0)          # u_{L-1} = u_15, rows = inputs
    mu=U.mean(0); muh=mu/np.linalg.norm(mu)
    R=np.linalg.norm(U,axis=1); P=U@muh; Up=U-np.outer(P,muh); Q=np.linalg.norm(Up,axis=1)
    V=Up/Q[:,None]; s2b=(Q/R)**2
    G=V@V.T; iu=np.triu_indices(N,1); rp=G[iu]
    Uh=U/R[:,None]; ct=(Uh@Uh.T)[iu]
    C=np.cov(Up.T); ev=np.linalg.eigvalsh(C); neff=ev.sum()**2/(ev**2).sum()
    print(f"n={n}: sin^2(beta) mean {s2b.mean():.4f} sd {s2b.std():.4f} | cos(theta) pairs mean {ct.mean():.4f} sd {ct.std():.4f} max {ct.max():.4f}")
    print(f"   transverse overlap rho_perp: rms {np.sqrt((rp**2).mean()):.4f}  (1/sqrt(n)={1/np.sqrt(n):.4f})  max|.| {np.abs(rp).max():.4f}  E|rho|^4 {np.mean(rp**4):.2e}  E|rho|^6 {np.mean(rp**6):.2e}")
    print(f"   n_eff from pairs 1/E rho^2 = {1/np.mean(rp**2):.1f};  participation ratio of transverse cov (N={N} samples) = {neff:.1f}")
    print(f"   Var log R^2 over inputs = {np.var(np.log(R**2)):.4f}")
