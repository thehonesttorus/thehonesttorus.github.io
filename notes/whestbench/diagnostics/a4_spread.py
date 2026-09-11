# Is the fourth-cumulant coefficient universal, or does it want to be estimated per network?
# The Monte-Carlo sample we already pay for can estimate it for free:
#   s_i = sum_a v_a zc_{i,a}^2   with v = sigma phi
#   Var(s) = sum_ab v_a v_b ( kappa_aabb + 2 Cov_ab^2 )   ->   g4_est = [Var(s) - 2 v^T (Sz.Sz) v]/(v.sigma^2)^2
# costing 2 N n + 3 n^2 flops per layer.  First check whether the per-network optimum even moves.
import numpy as np
import cov_terms as CT
GRID=(0.0008,0.0014,0.0020,0.0026,0.0034)
print(f"{'net':>4} " + " ".join(f"A4={g:<7.4f}" for g in GRID) + "   argmin")
for k in (0,1,2,3,4,5,6,7,9,11,13,15):
    v=[CT.predict(k,{'T1':g}) for g in GRID]
    i=int(np.argmin(v))
    print(f"{k:4d} " + " ".join(f"{x:10.4e}" for x in v) + f"   {GRID[i]:.4f}", flush=True)
