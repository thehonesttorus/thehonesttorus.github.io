# Merge mclatent.py chunks and split each true slice into collective, residual and cross parts.
#   python latent_an.py NET "CHUNK_GLOB"
# For v in (d, c, r) (full fluctuation, collective part U U^T d, residual): diagonal cumulants var, k3, k4, the
# covariance, D21_ij = kappa(v_i, v_i, v_j) and K31_ij = kappa(v_i, v_i, v_i, v_j) (mcstats.py's formulas). Since
# c and r are uncorrelated, var_d = var_c + var_r and C_d = C_c + C_r, so the gain shapes split too:
# S31_d = 3 var_d C_d = S31_c + S31_r + 3 (var_c C_r + var_r C_c). With one gain g fitted on d, the non-gain part
# R_d = K31_d - g S31_d splits exactly into R_c = K31_c - g S31_c, R_r = K31_r - g S31_r and the cross remainder.
# The latent model (gain times collective law plus Gaussian bulk) predicts R_d = R_c: the non-gain content collective.
# Chunks with even index form half 0 and odd index half 1; the half-difference gives the noise floor.
import sys, glob, re, numpy as np
from math import sqrt, pi, erf
net, pat = int(sys.argv[1]), sys.argv[2]
files = sorted(glob.glob(pat))
phi = lambda a: np.exp(-0.5 * a * a) / sqrt(2 * pi)
Phi = np.vectorize(lambda a: 0.5 * (1 + erf(a / sqrt(2))))
T = np.load(f"mc2_off{net}_full.npz")


def load(fs):
    tot = None
    for f in fs:
        z = np.load(f)
        if tot is None:
            tot = {k: np.array(z[k], np.float64) for k in z.files if k not in ("layers", "K")}
            tot["layers"] = z["layers"]; tot["K"] = int(z["K"])
        else:
            for k in z.files:
                if k not in ("layers", "K"):
                    tot[k] += z[k]
    return tot


def finish(a, l, v):
    c = a["c"]; g = lambda k: a[f"{k}_{v}_{l}"] / c
    md, e2, e3, e4 = g("s1"), g("s2"), g("s3"), g("s4")
    var = e2 - md * md; k3 = e3 - 3 * md * e2 + 2 * md ** 3
    k4 = e4 - 4 * md * e3 + 6 * md * md * e2 - 3 * md ** 4 - 3 * var * var
    E11, E21, E31 = g("m11"), g("m21"), g("m31")
    A_, B_ = md[:, None], md[None, :]; e2i, e3i = e2[:, None], e3[:, None]
    cov = E11 - A_ * B_
    D21 = E21 - 2 * A_ * E11 - B_ * e2i + 2 * A_ * A_ * B_
    Ex3y = E31 - B_ * e3i - 3 * A_ * E21 + 3 * A_ * B_ * e2i + 3 * A_ * A_ * E11 - 3 * A_ ** 3 * B_
    K31 = Ex3y - 3 * var[:, None] * cov
    for M in (D21, K31):
        np.fill_diagonal(M, 0.0)
    return dict(var=var, k3=k3, k4=k4, cov=cov, D21=D21, K31=K31)


chunks = {int(re.search(r"_(\d+)\.npz$", f).group(1)): f for f in files}
full = load(files); h0 = load([f for i, f in chunks.items() if i % 2 == 0]); h1 = load([f for i, f in chunks.items() if i % 2 == 1])
print(f"net {net}: {len(files)} chunks, {int(full['c'])} samples, K = {full['K']}")
en = lambda X, w=1.0: float(np.sum((X * w) ** 2))
for l in [int(x) for x in full["layers"]]:
    F = {v: finish(full, l, v) for v in ("d", "c", "r")}
    H0, H1 = finish(h0, l, "d"), finish(h1, l, "d")
    mu = np.asarray(T["mu"][l], np.float64); var = F["d"]["var"]; s = np.sqrt(var); al = mu / s
    c1, c2, c3, c4 = Phi(al), phi(al) / s, -al * phi(al) / var, (al * al - 1) * phi(al) / (var * s)
    out = [f"layer {l:2d}"]
    for nm, w in (("K31", c3[:, None] * c1[None, :]), ("D21", c2[:, None] * c1[None, :]), ("k4", c4), ("k3", c3)):
        X = {v: F[v][nm] for v in F}
        if nm in ("K31", "k4"):
            Sh = {v: (3 * F[v]["var"][:, None] * F[v]["cov"]) if nm == "K31" else 3 * F[v]["var"] ** 2 for v in F}
        else:
            # gain shapes of the third-order slices need the mean, which only d carries: 1.5 mu var, (2 mu_i C_ij + mu_j var_i)/2
            Sh = {"d": (0.5 * (2 * mu[:, None] * F["d"]["cov"] + F["d"]["var"][:, None] * mu[None, :])) if nm == "D21" else 1.5 * mu * F["d"]["var"]}
        if Sh["d"].ndim == 2:
            np.fill_diagonal(Sh["d"], 0.0)
            for v in ("c", "r"):
                if v in Sh:
                    np.fill_diagonal(Sh[v], 0.0)
        for tag, ww in (("plain", 1.0), ("read", w)):
            g = float(np.sum(X["d"] * Sh["d"] * ww * ww)) / en(Sh["d"], ww)
            Rd = X["d"] - g * Sh["d"]
            noise = en(H0[nm] - H1[nm], ww) / 4 / en(X["d"], ww)
            line = (f"  {nm:3s} {tag:5s}: g {g:+.4f}, non-gain share {en(Rd, ww) / en(X['d'], ww):.3f} (noise {noise:.3f})"
                    f" | energy shares of the slice: c {en(X['c'], ww) / en(X['d'], ww):.3f}, r {en(X['r'], ww) / en(X['d'], ww):.3f}")
            if "c" in Sh:
                Rc = X["c"] - g * Sh["c"]; Rr = X["r"] - g * Sh["r"]; Rx = Rd - Rc - Rr
                ex = lambda P: 1 - en(Rd - P, ww) / en(Rd, ww)
                line += (f" | non-gain part: explained by collective {ex(Rc):+.3f}, by residual {ex(Rr):+.3f}, by c+r {ex(Rc + Rr):+.3f};"
                         f" norms c {en(Rc, ww) / en(Rd, ww):.3f} r {en(Rr, ww) / en(Rd, ww):.3f} cross {en(Rx, ww) / en(Rd, ww):.3f}")
            else:
                # third order: no gain shape for c, r separately; report the collective share of the non-gain part directly
                ex = lambda P: 1 - en(Rd - P, ww) / en(Rd, ww)
                line += f" | slice explained by collective {1 - en(X['d'] - X['c'], ww) / en(X['d'], ww):+.3f}"
            out.append(line)
    print("\n".join(out), flush=True)
