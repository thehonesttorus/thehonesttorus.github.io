"""Layer-by-layer tolerance that a carrier of the QUENCHED remainder of old content must meet at n = 1024 x 16.

Price per layer (final MSE per unit squared relative D21 error at layer k):
    P_k = c * K_off(k) * ||dC_k||_F^2,  K_off from the wedge calculus (Proposition E1),  ||dC_k||^2 from region's MC atlas
    (MLP 0, odd layers, interpolated),  c fixed so that sum_k P_k = 3.76e-6 (region: dropping D21 entirely costs that).
Quenched remainder at layer k after a cheap carrier (ages <= A explicit + dilation sector of all ages exact):
    r_k(A) = sum_{age > A} E_{s->k} (1 - share_{s->k}) / Etot_k     (team C's agespec, 6 networks, residuals ~orthogonal)
Budget B_q for the remainder; equal-price allocation eps_k^2 = B_q / (#layers * P_k); the carrier must reproduce the
remainder to relative accuracy delta_k = eps_k / sqrt(r_k) (>= 1 means the layer's remainder can be dropped)."""
import json, numpy as np, sys
exec(open('wedge_check_any.py').read().split("reg = ")[0])
n, L = 1024, 16
Km, Ko = wedge(n, L)
odd = {1: 0.71, 3: 0.82, 5: 0.78, 7: 0.74, 9: 0.69, 11: 0.50, 13: 0.51, 15: 0.42}
dC = np.interp(np.arange(L), sorted(odd), [odd[k] for k in sorted(odd)])
lay = np.arange(1, L - 1)                       # covariance channel: layers 1..14 (K_off(15) = 0)
P = Ko[lay] * dC[lay]; P *= 3.76e-6 / P.sum()
C = '../C-free-probability-criticality/results/agespec_w1024_d16_%d.json'
R = {A: np.zeros(len(lay)) for A in (0, 2, 4, 6)}
for i in range(6):
    d = json.load(open(C % i)); Et = {x['layer']: x['Etot'] for x in d['layers']}
    for A in R:
        for j, k in enumerate(lay):
            rem = sum(r['E'] * (1 - r['share']) for r in d['rows'] if r['k'] == k and r['age'] > A)
            R[A][j] += rem / Et[k] / 6
for Bq in (5e-9, 1.5e-8):
    eps = np.sqrt(Bq / (len(lay) * P))
    print(f'\nbudget for the quenched remainder B_q = {Bq:.1e} (bar = 1.5e-8); uniform eps would be {np.sqrt(Bq / P.sum()):.3f}')
    print(' k  price-share  eps_k(D21) |  r_k(A): remainder share of D21 energy, A=0/2/4/6  |  required delta_k = eps/sqrt(r), A=0/2/4/6')
    for j, k in enumerate(lay):
        rr = [R[A][j] for A in (0, 2, 4, 6)]
        dd = ['  -  ' if r <= 0 else f'{min(eps[j] / np.sqrt(r), 9.99):5.2f}' for r in rr]
        print(f'{k:2d}   {P[j] / P.sum():.3f}       {eps[j]:.3f}     |  ' + ' '.join(f'{r:.3f}' for r in rr) + '  |  ' + ' '.join(dd))
