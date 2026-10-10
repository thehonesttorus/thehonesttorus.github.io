# Note XLIV section 9 (active-set version, alpha = mu / S > -2.5 only): how much of each slice's energy is gain-shaped? cos^2 between the slice and its package form,
# for the truth (Monte Carlo halves averaged) and the chain's own slices, at every layer; plus the same for the chain's error
# (chain minus truth) projected on the package form (noise-free through the two halves).
#   python -I gain_share.py LOCDIR MC2DIR CD2DIR NET
import sys, numpy as np
loc, mcd, cdd, net = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
c2 = np.load(f"{cdd}/chaindump2_{net}.npz")
F = {}
for h in ("h0", "h1"):
    z = np.load(f"{mcd}/mc2_off{net}_{h}.npz")
    F[h] = {k: z[k] for k in ("mu", "var", "k3", "k4", "cov", "D21", "K22", "K31")}
n = F["h0"]["mu"].shape[1]
off = ~np.eye(n, dtype=bool)
f64 = lambda a: np.asarray(a, dtype=np.float64)
print(f"=== network {net}: energy share of the gain-shaped package, truth | chain (cos^2 x 100) and share of the chain error energy")
print("layer   k3: truth|chain|err   D21: truth|chain|err   K22: truth|chain|err   K31: truth|chain|err   k4: truth|chain|err")
for s in range(3, 15):
    mu = 0.5 * (f64(F["h0"]["mu"][s]) + f64(F["h1"]["mu"][s])); v = 0.5 * (f64(F["h0"]["var"][s]) + f64(F["h1"]["var"][s]))
    C = 0.5 * (f64(F["h0"]["cov"][s]) + f64(F["h1"]["cov"][s])); C = 0.5 * (C + C.T)
    act = (mu / np.sqrt(v)) > -2.5
    mk = off & act[:, None] & act[None, :]
    forms = dict(k3=6 * mu * v * act, k4=12 * v * v * act,
                 D21=np.where(mk, 2 * (2 * mu[:, None] * C + mu[None, :] * v[:, None]), 0.0),
                 K22=np.where(mk, 4 * np.outer(v, v) + 8 * C * C, 0.0), K31=np.where(mk, 12 * v[:, None] * C, 0.0))
    sl = dict(k3=f64(c2[f"D3_{s}"]), k4=f64(c2[f"g4row_{s}"]), D21=f64(c2[f"D21_{s}"]),
              K22=0.5 * (f64(c2[f"wk4m_{s}"]) + f64(c2[f"wk4m_{s}"]).T), K31=f64(c2[f"wk431_{s}"]).T)
    dot = lambda a, b: float(np.sum(a * b * (mk if a.ndim == 2 else act)))
    row = []
    for k in ("k3", "D21", "K22", "K31", "k4"):
        f = forms[k]
        T = [f64(F[h][k][s]) for h in ("h0", "h1")]
        Tm = 0.5 * (T[0] + T[1])
        cos2 = lambda x: dot(f, x) ** 2 / (dot(f, f) * dot(x, x))
        ch = sl[k]
        e0, e1 = ch - T[0], ch - T[1]
        shr = 0.5 * (dot(f, e0) * dot(f, e1) + dot(f, e1) * dot(f, e0)) / (dot(f, f) * 0.5 * (dot(e0, e1) + dot(e1, e0)))
        row.append(f"{100 * cos2(Tm):5.1f}|{100 * cos2(ch):5.1f}|{100 * shr:5.1f}")
    print(f"  {s:2d}   " + "   ".join(row), flush=True)
