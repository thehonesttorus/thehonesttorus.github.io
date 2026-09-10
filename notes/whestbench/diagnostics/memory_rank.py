# Can the many cumulant source channels be compressed into a small memory state?
# Measure the contribution c_s of each source layer s to kappa_3 at the FINAL layer, then ask how
# low-rank the family {c_s} is.  This is the input-output (Hankel) question, not the transport-spectrum
# question that failed earlier: it asks whether the sources' EFFECT compresses, not whether the
# transport does.
import numpy as np
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; pdf=T.pdf
seeds,Yall=T.load_nets(); HS=np.load('../pass2_stats.npz')['HS']
W=T.regen(seeds[0]); Y=Yall[0].astype(np.float64)
srcs={}
for l in range(L):
    mu=(Y[l-1]@W[l].astype(np.float64)) if l>0 else np.zeros(n)
    Sz=W[l].T.astype(np.float64)@HS[l]@W[l].astype(np.float64); sg=np.sqrt(np.diag(Sz)); t=mu/sg
    ph=pdf(t); Ph=ndtr(t); m1=mu*Ph+sg*ph; m2=(mu*mu+sg*sg)*Ph+mu*sg*ph; m3=mu*m2+2*sg*sg*m1
    M=Ph[:,None]*Sz.copy(); np.fill_diagonal(M,0)
    srcs[l]=(M, ph/sg, 2*m1*(1-Ph), m3-3*m2*m1+2*m1**3, Ph)
l=15; Wl=W[l].astype(np.float64)
mu=Y[l-1]@Wl; Sz=W[l].T.astype(np.float64)@HS[l]@Wl; sg=np.sqrt(np.diag(Sz)); t=mu/sg
U=Wl.copy(); C=[]
for s in range(l-1,l-13,-1):
    M,p,q,k3h,d1=srcs[s]; Nt=M.T@U
    C.append(((p[:,None]*3)*(U*Nt)*Nt+(q[:,None]*3)*(U*Nt)*U+k3h[:,None]*(U*U*U)).sum(0))
    U=W[s]@(d1[:,None]*U)
C=np.array(C)                      # (source age, output neuron)
Cn=C/np.linalg.norm(C,axis=1,keepdims=True)
print("Per-source contributions to kappa_3(z_16), network 0.  Source age 1 = previous layer.\n")
print("pairwise correlation between per-source contribution vectors:")
print("      " + " ".join(f"a{j+1:<5d}" for j in range(8)))
for i in range(8):
    print(f"  a{i+1:<3d}" + " ".join(f"{np.dot(Cn[i],Cn[j]):+.3f}" for j in range(8)))
sv=np.linalg.svd(C,compute_uv=False); print(f"\nsingular values of the 12 x {n} source-contribution stack, as share of total energy:")
print("  " + " ".join(f"{x:.3f}" for x in (sv**2/np.sum(sv**2))))
print(f"  cumulative: " + " ".join(f"{x:.3f}" for x in np.cumsum(sv**2)/np.sum(sv**2)))
tail=C[3:].sum(0); full=C.sum(0)
print(f"\ntail (sources of age >= 4) as share of full: {np.linalg.norm(tail)/np.linalg.norm(full):.3f}")
s3=sg**3
for name,B in [("{sigma^3}",np.column_stack([s3])),
               ("{sigma^3, t sigma^3}",np.column_stack([s3,t*s3])),
               ("{sigma^3, t sigma^3, t^2 sigma^3}",np.column_stack([s3,t*s3,t*t*s3])),
               ("age-1 contribution direction",np.column_stack([C[0]])),
               ("age-1,2,3 directions",C[:3].T),
               ("{sigma^3,t sigma^3} + age-1,2,3",np.column_stack([s3,t*s3,C[:3].T]))]:
    cf=np.linalg.lstsq(B,tail,rcond=None)[0]; r2=1-np.sum((tail-B@cf)**2)/np.sum(tail**2)
    print(f"  R2 of the tail from {name:36s} (rank {B.shape[1]}): {r2:.3f}")
