import numpy as np
lam = lambda m: (17 * 9.0**(m-1) - 1) / 4
print("Seed length M(t,eps): smallest M with sum_{m>M} c_m e^{-2 t lam_m} <= eps^2")
for kind, c in [("classical (2^{m-1} characters)", lambda m: 2.0**(m-1)), ("noncommutative (3*4^{m-1} Weyl labels)", lambda m: 3*4.0**(m-1))]:
    print(" ", kind)
    for t in [1e-4, 1e-3, 1e-2, 1e-1]:
        row = []
        for eps in [1e-3, 1e-6, 1e-12, 1e-24]:
            M = 0
            while sum(c(m) * np.exp(-2*t*lam(m)) for m in range(M+1, M+60)) > eps**2: M += 1
            row.append(M)
        print(f"    t={t:g}: M for eps=1e-3,1e-6,1e-12,1e-24 -> {row}")
spec = np.load("spec.npy")
print("Galton-Watson survival P(Z_D>0) = 1 - b_0^(D) vs Slack asymptote 9 pi^2/(2 D^2):")
for D in [1, 2, 4, 8, 16]:
    print(f"   D={D:2d}: {1-spec[D-1][0]:.4f}   asymptote {9*np.pi**2/(2*D**2):.4f}")
# conditional mean of Z_D given survival
for D in [4, 8, 16]:
    p = spec[D-1]; r = np.arange(len(p))
    print(f"   D={D:2d}: E[Z_D]={np.sum(r*p):.4f} (critical => 1), E[Z_D | Z_D>0]={np.sum(r*p)/(1-p[0]):.2f}, median order given survival = {np.searchsorted(np.cumsum(p[1:])/(1-p[0]), 0.5)+1}")
