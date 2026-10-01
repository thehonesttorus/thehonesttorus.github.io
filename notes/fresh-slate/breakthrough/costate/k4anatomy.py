"""Top modes of the true joint kappa_4 slices (MC atlas) vs scale-mode templates."""
import sys, numpy as np
f = sys.argv[1]; st = dict(np.load(f))
cos = lambda a, b: abs(a @ b) / np.linalg.norm(a) / np.linalg.norm(b)
for l in (3, 7, 11, 15):
    c2 = st["c2"][l]; v = np.diag(c2); m = st["m"][l] if "m" in st else None
    K22 = st["s22"][l] - np.outer(v, v) - 2 * c2 * c2; K31 = st["s31"][l] - 3 * v[:, None] * c2
    for nm, K in (("K22", K22), ("K31", K31)):
        O = K - np.diag(np.diag(K)); U, s, Vt = np.linalg.svd(O); e = s ** 2 / np.sum(s ** 2)
        cands = dict(s2=v, s=np.sqrt(v), one=np.ones_like(v))
        if m is not None: cands.update(m=m, m2=m * m, s2m2=v + m * m)
        al = " ".join(f"{k}:{cos(U[:, 0], x):.2f}/{cos(Vt[0], x):.2f}" for k, x in cands.items())
        print(f"L{l:2d} {nm} off-diag energy top1 {e[0]:.2f} top8 {e[:8].sum():.2f} top32 {e[:32].sum():.2f} | u1/v1: {al}")
