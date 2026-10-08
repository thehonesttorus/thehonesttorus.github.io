# Merge mcchunk.py chunks into the files mcstats.py writes: OUTPREFIX_h0.npz, OUTPREFIX_h1.npz, OUTPREFIX_full.npz.
#   python mcmerge.py OUTPREFIX "H0_GLOB" "H1_GLOB"
# The two halves are disjoint sets of chunks (independent seeds); finish() is mcstats.py's, unchanged.
import sys, glob, numpy as np
outp, g0, g1 = sys.argv[1], sys.argv[2], sys.argv[3]
KEYS = ("s1", "s2", "s3", "s4", "m11", "m21", "m22", "m31")


def load(files):
    acc, accy, m0, m0y = None, None, None, None
    for f in files:
        z = np.load(f)
        if acc is None:
            m0, m0y = z["m0"], z["m0y"]
            acc = dict(c=0, yfin=np.zeros_like(z["yfin"]), **{k: np.zeros_like(z[k]) for k in KEYS})
            if "cy" in z:
                accy = dict(c=0, **{k: np.zeros_like(z[k + "_y"]) for k in KEYS})
        assert np.array_equal(m0, z["m0"]), "chunks from different passes"
        acc["c"] += int(z["c"]); acc["yfin"] += z["yfin"]
        for k in KEYS:
            acc[k] += z[k]
        if accy is not None:
            accy["c"] += int(z["cy"])
            for k in KEYS:
                accy[k] += z[k + "_y"]
    return acc, accy, m0, m0y


def finish(a, shift):
    c = a["c"]; md = a["s1"] / c; e2 = a["s2"] / c; e3 = a["s3"] / c; e4 = a["s4"] / c
    var = e2 - md * md
    k3 = e3 - 3 * md * e2 + 2 * md ** 3
    k4 = e4 - 4 * md * e3 + 6 * md * md * e2 - 3 * md ** 4 - 3 * var * var
    E11 = a["m11"] / c; E21 = a["m21"] / c
    D21 = E21 - 2 * md[:, :, None] * E11 - md[:, None, :] * e2[:, :, None] + 2 * (md * md)[:, :, None] * md[:, None, :]
    E22 = a["m22"] / c
    A_, B_ = md[:, :, None], md[:, None, :]
    E21T = np.swapaxes(E21, 1, 2)
    e2i, e2j = e2[:, :, None], e2[:, None, :]
    Ex2y2 = (E22 - 2 * B_ * E21 - 2 * A_ * E21T + B_ * B_ * e2i + A_ * A_ * e2j + 4 * A_ * B_ * E11
             - 3 * A_ * A_ * B_ * B_)
    cov = E11 - A_ * B_
    K22 = Ex2y2 - var[:, :, None] * var[:, None, :] - 2 * cov * cov
    E31 = a["m31"] / c; e3i = e3[:, :, None]
    Ex3y = E31 - B_ * e3i - 3 * A_ * E21 + 3 * A_ * B_ * e2i + 3 * A_ * A_ * E11 - 3 * A_ ** 3 * B_
    K31 = Ex3y - 3 * var[:, :, None] * cov
    out = dict(n=c, mu=(shift + md).astype(np.float32), var=var.astype(np.float32), k3=k3.astype(np.float32),
               k4=k4.astype(np.float32), D21=D21.astype(np.float32), cov=cov.astype(np.float32),
               K22=K22.astype(np.float32), K31=K31.astype(np.float32))
    if "yfin" in a:
        out["yfin"] = (a["yfin"] / c).astype(np.float32)
    return out


f0s, f1s = sorted(glob.glob(g0)), sorted(glob.glob(g1))
assert f0s and f1s and not set(f0s) & set(f1s)
h = [load(f0s), load(f1s)]
res = []
for acc, accy, m0, m0y in h:
    d = finish(acc, m0)
    if accy is not None:
        d.update({f"{k}_y": v for k, v in finish(accy, m0y).items() if k not in ("n", "yfin")})
    res.append(d)
acc = {k: h[0][0][k] + h[1][0][k] for k in ("c", "yfin") + KEYS}
full = finish(acc, h[0][2])
if h[0][1] is not None:
    accy = {k: h[0][1][k] + h[1][1][k] for k in ("c",) + KEYS}
    full.update({f"{k}_y": v for k, v in finish(accy, h[0][3]).items() if k not in ("n", "yfin")})
for tag, d in (("h0", res[0]), ("h1", res[1]), ("full", full)):
    np.savez(f"{outp}_{tag}.npz", **d)
print(f"merged {len(f0s)} + {len(f1s)} chunks: {res[0]['n']} + {res[1]['n']} samples -> {outp}_{{h0,h1,full}}.npz", flush=True)
