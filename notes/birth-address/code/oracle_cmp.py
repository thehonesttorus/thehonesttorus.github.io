# The chain's third-cumulant readouts against Monte Carlo truth, per layer (note XXXI).   python oracle_cmp.py NET PREFIX
import sys, numpy as np
net, pre = int(sys.argv[1]), sys.argv[2]
ch = np.load(f"chain_off{net}.npz"); F = np.load(f"{pre}_full.npz"); H0 = np.load(f"{pre}_h0.npz"); H1 = np.load(f"{pre}_h1.npz")
print(f"net {net}: {int(F['n'])} samples")
off = ~np.eye(1024, dtype=bool)
for l in range(1, 16):
    if f"D3_{l}" not in ch.files:
        continue
    c3, t3 = ch[f"D3_{l}"].astype(np.float64), F["k3"][l].astype(np.float64)
    nz3 = np.linalg.norm(H0["k3"][l] - H1["k3"][l]) / 2 / np.linalg.norm(t3)        # MC noise of the full estimate
    cv, tv = ch[f"var_{l}"].astype(np.float64), F["var"][l].astype(np.float64)
    line = (f"layer {l:2d}: var rel err {np.linalg.norm(cv - tv) / np.linalg.norm(tv):.2e} | D3 corr {np.corrcoef(c3, t3)[0, 1]:.4f} "
            f"rel err {np.linalg.norm(c3 - t3) / np.linalg.norm(t3):.3f} (MC noise {nz3:.3f}) skew rms {np.sqrt(np.mean(t3 ** 2 / tv ** 3)):.3f}")
    if f"D21_{l}" in ch.files:
        c21, t21 = ch[f"D21_{l}"].astype(np.float64)[off], F["D21"][l].astype(np.float64)[off]
        nz21 = np.linalg.norm((H0["D21"][l] - H1["D21"][l])[off]) / 2 / np.linalg.norm(t21)
        line += (f" | D21 corr {np.corrcoef(c21, t21)[0, 1]:.4f} rel err {np.linalg.norm(c21 - t21) / np.linalg.norm(t21):.3f} "
                 f"(MC noise {nz21:.3f})")
    print(line, flush=True)
