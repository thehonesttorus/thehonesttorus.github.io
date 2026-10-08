# Renormalized comparison of closures (note XXXIX). The chain's fitted elements are counterterms: they absorb part of the
# dynamics the closure omits. Adding a derived term T to the chain while keeping the counterterms at the values fitted
# without it counts that physics twice. The fair comparison renormalizes both: the counterterm amplitudes (the 103
# per-layer amplitudes of the statistics the Wick stage reads, V47_CAL) are refitted in the output metric for the base
# and for base + T alike, on the training networks, and the two are compared on held-out networks.
# The refit uses the base's linear responses R_j (first order in the amplitudes, as in calfit.py):
#   out(X + cal_b) ~ out_X + R b,   b = argmin sum_train |e_X + R b|^2 + ridge (chosen inside the training set).
#   python refit.py RESPDIR RBASEDIR DIRS_JSON BASEPAT NAME=PAT [NAME=PAT ...] [--train 0-49] [--val 50-99]
# PAT is a path pattern with {n}, e.g. out/o35_{n}.npy (every layer's means, as SAVE_OUT writes them).
# Prints, on the held-out networks: raw change of X against the base (no refit), the refitted base, the refitted X, the
# renormalized change (refitted X against refitted base, paired), and the free amplitude a of T when the direction
# d_T = out_X - out_base is added to the counterterm family. Writes refit_coef_NAME.json (fit on the training networks)
# and refit_coef_NAME_all.json (fit on every network, for submission).
import sys, os, json, numpy as np


def nets(s):
    a, _, b = s.partition("-"); return list(range(int(a), int(b) + 1)) if b else [int(a)]


args = [a for a in sys.argv[1:]]
tr = nets(args[args.index("--train") + 1]) if "--train" in args else list(range(50))
va = nets(args[args.index("--val") + 1]) if "--val" in args else list(range(50, 100))
for flag in ("--train", "--val"):
    if flag in args:
        i = args.index(flag); del args[i:i + 2]
rd, rbd, dj, bpat = args[:4]
cands = [a.split("=", 1) for a in args[4:]]
dirs = json.load(open(dj)); K = len(dirs)
off = os.environ.get("OFFICIAL", "../official")
RIDGES = (0.0, 1e-4, 1e-3, 1e-2, 1e-1, 1.0)


def have(n, pat):
    return os.path.exists(pat.format(n=n))


resp, truth = {}, {}
for n in tr + va:
    if not all(os.path.exists(f"{rd}/out_c{j}_{n}.npy") for j in range(K)) or not os.path.exists(f"{rbd}/out_b_{n}.npy"):
        continue
    b = np.load(f"{rbd}/out_b_{n}.npy")[-1]
    resp[n] = np.stack([(np.load(f"{rd}/out_c{j}_{n}.npy")[-1] - b) / dirs[j][2] for j in range(K)], axis=1)
    truth[n] = np.load(f"{off}/truth_off{n}.npz")["m"].astype(np.float64)[-1]
base = {n: np.load(bpat.format(n=n))[-1] for n in resp if have(n, bpat)}


def fit(err, cols, ns):
    # err: n -> (1024,), cols: n -> (1024, k). Ridge chosen on the second half of ns from a fit on the first half.
    k = next(iter(cols.values())).shape[1]
    def gram(sub):
        G = sum(cols[m].T @ cols[m] for m in sub); g = sum(cols[m].T @ err[m] for m in sub); return G, g
    def solve(G, g, rg):
        return -np.linalg.solve(G + rg * np.trace(G) / k * np.eye(k), g)
    h1, h2 = ns[:len(ns) // 2], ns[len(ns) // 2:]
    G1, g1 = gram(h1)
    rg = min(RIDGES, key=lambda r: np.mean([np.mean((err[m] + cols[m] @ solve(G1, g1, r)) ** 2) for m in h2]))
    G, g = gram(ns)
    return solve(G, g, rg), rg


def mse(err, cols, b, ns):
    return np.array([np.mean((err[m] + (0 if b is None else cols[m] @ b)) ** 2) for m in ns])


def pct(x1, x0):
    r = x1 / x0 - 1
    return 100 * (x1.mean() / x0.mean() - 1), 100 * r.std() / np.sqrt(len(r)), int((x1 < x0).sum())


eb = {n: base[n] - truth[n] for n in base}
trb = [n for n in tr if n in eb]; vab = [n for n in va if n in eb]
bb, rgb = fit(eb, resp, trb)
mb0, mb1 = mse(eb, resp, None, vab), mse(eb, resp, bb, vab)
print(f"{K} directions; base {bpat}: train {len(trb)}, held-out {len(vab)}")
print(f"  base: held-out raw {mb0.mean():.5e}; refitted (ridge {rgb:g}) {mb1.mean():.5e} "
      "({:+.2f}% +- {:.2f}, better on {})".format(*pct(mb1, mb0)))
json.dump({"dirs": dirs, "coef": [float(x) for x in bb], "ridge": rgb}, open("refit_coef_base.json", "w"))
ba, _ = fit(eb, resp, trb + vab)
json.dump({"dirs": dirs, "coef": [float(x) for x in ba]}, open("refit_coef_base_all.json", "w"))
for name, pat in cands:
    ec = {n: np.load(pat.format(n=n))[-1] - truth[n] for n in eb if have(n, pat)}
    trc = [n for n in trb if n in ec]; vac = [n for n in vab if n in ec]
    if not trc or not vac:
        print(f"  {name}: missing outputs ({len(trc)} train, {len(vac)} held-out)"); continue
    bc, rgc = fit(ec, resp, trc)
    mc0, mc1 = mse(ec, resp, None, vac), mse(ec, resp, bc, vac)
    mb0v, mb1v = mse(eb, resp, None, vac), mse(eb, resp, bb, vac)
    # free amplitude of T: the direction d_T = out_X - out_base appended to the counterterm family, fitted with it
    ext = {n: np.concatenate([resp[n], (ec[n] - eb[n])[:, None]], axis=1) for n in ec}
    bf, rgf = fit(eb, ext, trc)
    mf1 = mse(eb, ext, bf, vac)
    print(f"  {name}: held-out raw vs base {'{:+.2f}% +- {:.2f} (better on {})'.format(*pct(mc0, mb0v))}; "
          f"refitted {mc1.mean():.5e} (ridge {rgc:g}); renormalized change vs refitted base "
          + "{:+.2f}% +- {:.2f} (better on {}/".format(*pct(mc1, mb1v)) + f"{len(vac)})"
          + f"; free amplitude a = {bf[-1]:.3f}, refitted with a free: " + "{:+.2f}% +- {:.2f}".format(*pct(mf1, mb1v)[:2]))
    json.dump({"dirs": dirs, "coef": [float(x) for x in bc], "ridge": rgc}, open(f"refit_coef_{name}.json", "w"))
    ca, _ = fit(ec, resp, trc + vac)
    json.dump({"dirs": dirs, "coef": [float(x) for x in ca]}, open(f"refit_coef_{name}_all.json", "w"))
