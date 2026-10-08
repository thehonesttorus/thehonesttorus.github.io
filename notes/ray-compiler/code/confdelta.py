# Output perturbation of the shared-basis confinement: delta = out(variant) - out(production), against the chain's own
# error e = out(production) - truth.  MSE(variant) - MSE(prod) = 2 <e, delta>/n + |delta|^2/n.
#   python confdelta.py DIR NETS    (DIR holds out_prod_N.npy, out_nc_N.npy, out_1t_N.npy from oracle_one.py SAVE_OUT)
import sys, numpy as np
d, nets = sys.argv[1], [int(x) for x in sys.argv[2].split(",")]
off = "../official"   # directory holding truth_off{n}.npz
rows = {"nc": [], "1t": []}
for n in nets:
    m = np.load(f"{off}/truth_off{n}.npz")["m"].astype(np.float64)
    p = np.load(f"{d}/out_prod_{n}.npy")
    e = p[-1] - m[-1]; E = float(np.mean(e * e))
    for v in ("nc", "1t"):
        try:
            q = np.load(f"{d}/out_{v}_{n}.npy")
        except FileNotFoundError:
            continue
        dl = q[-1] - p[-1]
        D2 = float(np.mean(dl * dl)); X = float(2 * np.mean(e * dl))
        cos = float(e @ dl / (np.linalg.norm(e) * np.linalg.norm(dl) + 1e-300))
        prof = [float(np.mean((q[l] - p[l]) ** 2) / max(np.mean((p[l] - m[l]) ** 2), 1e-300)) for l in range(q.shape[0])]
        rows[v].append((n, E, D2, X, cos, prof))
        print(f"net {n} {v}: MSE prod {E:.4e} | |delta|^2/n {D2:.3e} = {D2 / E:.4f} MSE | cross 2<e,delta>/n {X / E:+.4f} MSE | "
              f"cos(e, delta) {cos:+.3f} | net change {(D2 + X) / E:+.4f}")
for v, rs in rows.items():
    if not rs:
        continue
    E = np.array([r[1] for r in rs]); D2 = np.array([r[2] for r in rs]); X = np.array([r[3] for r in rs])
    P = np.mean([r[5] for r in rs], axis=0)
    print(f"{v}: mean |delta|^2 / MSE {np.mean(D2 / E):.4f} | mean cross {np.mean(X / E):+.4f} | mean net {np.mean((D2 + X) / E):+.4f} | "
          f"mean cos {np.mean([r[4] for r in rs]):+.3f}")
    print(f"   |delta_l|^2 / |e_l|^2 by layer: " + " ".join(f"{x:.3f}" for x in P))
