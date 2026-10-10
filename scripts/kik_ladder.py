"""Shape ladder at every layer from the pooled contractions (python scripts/kik_ladder.py RESULTS DATA NET ...).

For each layer's pre-activation z_i: Gaussian at (t, q), and Hermite (Gram-Charlier) shapes with the exact kappa_3, and
kappa_3 + kappa_4 (+ the kappa_3^2 He_6 term), each compared with the same-sample mean E z_+ (sampling noise cancels).
Also the first-order sensitivities of the final mean to each contraction, which set the accuracy a compiler must reach:
dm/dt = P(Z > 0) (exact, Gaussian value used), dm/dq, dm/dkappa3, dm/dkappa4 from the Hermite formula."""
import sys, os, glob, numpy as np
from scipy.special import ndtr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scripts.kik_merge import pool_net

def hermite_means(t, q, m3, m4):
    sd = np.sqrt(q - t * t); al = t / sd; ph = np.exp(-al * al / 2) / np.sqrt(2 * np.pi)
    c3 = (m3 - 3 * t * q + 2 * t ** 3) / sd ** 3; c4 = (m4 - 4 * t * m3 + 6 * t * t * q - 3 * t ** 4) / sd ** 4 - 3
    d0 = al * ndtr(al) + ph; d3 = -al * ph; d4 = (al * al - 1) * ph; d6 = (al ** 4 - 6 * al * al + 3) * ph
    return sd * d0, sd * (d0 + c3 / 6 * d3), sd * (d0 + c3 / 6 * d3 + c4 / 24 * d4 + c3 ** 2 / 72 * d6), c3, c4, sd, al, ph

rms = lambda x: np.sqrt(np.mean(x ** 2))
R, D = sys.argv[1], sys.argv[2]
for net in map(int, sys.argv[3:]):
    a = pool_net(R, net); N = a["N"]; mc = a["pos"] / N; L = mc.shape[0]
    print(f"net {net}: same-sample shape errors per layer (RMS over neurons):  gauss | +k3 | +k3+k4   [rms k3, mean k4]")
    for l in range(L):
        P = (a["pw_low"][l] if l < L - 1 else a["pw_top"]) / N
        g, h3, h34, c3, c4, sd, al, ph = hermite_means(P[:, 0], P[:, 1], P[:, 2], P[:, 3])
        if l in (0, 1, 2, 4, 8, 12, 14, 15):
            print(f"  layer {l:2d}: {rms(g - mc[l]):.2e} | {rms(h3 - mc[l]):.2e} | {rms(h34 - mc[l]):.2e}   [{rms(c3):.3f}, {c4.mean():+.3f}]")
    # sensitivities at the last layer (Hermite-Gaussian part): dm/dt, dm/dq at fixed central shape, dm/dk3, dm/dk4
    dt = ndtr(al); dq = ph / (2 * sd); dk3 = sd * (-al * ph) / 6; dk4 = sd * (al * al - 1) * ph / 24
    tgt = 6e-5
    print(f"  final-layer accuracy each contraction needs for {tgt:.0e} RMS alone: t {tgt / rms(dt):.1e}, q {tgt / rms(dq):.1e}, "
          f"kappa3 {tgt / rms(dk3):.1e} (rms kappa3 {rms(c3):.3f}), kappa4 {tgt / rms(dk4):.1e} (mean kappa4 {c4.mean():.3f})")
