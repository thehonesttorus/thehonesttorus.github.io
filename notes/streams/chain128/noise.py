"""Per-layer Monte Carlo noise of one atlas's pre-activation objects, from two independent atlases of the same MLP:
noise_l = rel(A_l, B_l) / sqrt(2) for var, Coff, D3, D21, K4, K22.  Writes results/noise_mlp<i>.json.
    python noise.py A.npz B.npz <mlp index>
"""
import json, sys
import numpy as np
import chain as ch
import oracle_k3 as ok


def objs(z, l):
    st = ch.atlas_state(z, l)
    D3, D21, K4, K31, K22 = ch.slices_from_state(st)
    return dict(var=st.var, Coff=ok.offdiag(st.C), D3=D3, D21=D21, K4=K4, K22=K22)


def main():
    A, B = np.load(sys.argv[1]), np.load(sys.argv[2])
    L = A["weights"].shape[0]
    out = {k: [] for k in ("var", "Coff", "D3", "D21", "K4", "K22")}
    for l in range(L):
        a, b = objs(A, l), objs(B, l)
        for k in out:
            den = np.sum(b[k] ** 2)
            out[k].append(float(np.sqrt(np.sum((a[k] - b[k]) ** 2) / den / 2)) if den > 0 else 0.0)
    json.dump(out, open(f"results/noise_mlp{sys.argv[3]}.json", "w"))
    for k, v in out.items():
        print(k, " ".join(f"{x:.3f}" for x in v))


if __name__ == "__main__":
    main()
