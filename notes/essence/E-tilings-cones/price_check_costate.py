"""Check the price list (wedge K_off x region atlas ||dC||^2, normalised to 3.76e-6) against costate's measured
truncations inside FC at n = 1024 (CONVERGENCE.md 'Round 3, first results', MLP 0; FC all ages 3.24e-8).
Predicted added MSE of dropping content older than A (C's age = k - s - 1, age 0 = newest source):
   sum_k P_k * err_k / Etot_k,  err_k = |sum_{age>A} gamma|^2 KK_k (dilation, coherent) + sum_{age>A} E(1-share) (free),
or the free part only when the dilation sector is carried exactly."""
import json, numpy as np
exec(open('wedge_check_any.py').read().split("reg = ")[0])
n, L = 1024, 16
Km, Ko = wedge(n, L)
odd = {1: 0.71, 3: 0.82, 5: 0.78, 7: 0.74, 9: 0.69, 11: 0.50, 13: 0.51, 15: 0.42}
dC = np.interp(np.arange(L), sorted(odd), [odd[k] for k in sorted(odd)])
lay = np.arange(1, L - 1)
P = Ko[lay] * dC[lay]; P *= 3.76e-6 / P.sum()
d = json.load(open('../C-free-probability-criticality/results/agespec_w1024_d16_0.json'))
lays = {x['layer']: x for x in d['layers']}
def added(A, dil):
    tot = 0.0
    for j, k in enumerate(lay):
        rows = [r for r in d['rows'] if r['k'] == k and r['age'] > A]
        free = sum(r['E'] * (1 - r['share']) for r in rows)
        coh = (sum(r['gamma'] for r in rows)) ** 2 * lays[k]['KK']
        tot += P[j] * (free + (0 if dil else coh)) / lays[k]['Etot']
    return tot
meas = {('2', False): 1.77e-6, ('2', True): 7.67e-7, ('4', False): 7.79e-7, ('4', True): 3.84e-7}
base = 3.24e-8
for lab in ('2', '4'):
    for dil in (False, True):
        for A in (int(lab) - 1, int(lab)):
            print(f"costate 'ages <= {lab}'{' + dilation' if dil else ''}: measured {meas[(lab, dil)]:.2e}; "
                  f"predicted (C age > {A} dropped) {base + added(A, dil):.2e}")
