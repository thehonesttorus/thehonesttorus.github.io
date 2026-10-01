"""Compare a chain's per-layer pre-activation state with an atlas's truth (relative rms errors per layer).
    python diag_state.py atlas.npz "dict(k3mode='closure')" ["dict(...)" ...]
columns: mu, var, Coff, D3, D21 of z_l; K4 (diagonal kappa4), K22 (off-diagonal (2,2) slice), X211 (all-distinct
(2,1,1) slice, 'mem' closure only); final column: MSE of the post-activation means vs the atlas's means.
"""
import sys
import numpy as np
import chain as ch
import oracle_k3 as ok


def rel(a, b):
    return float(np.sqrt(np.sum((a - b) ** 2) / max(np.sum(b ** 2), 1e-300)))


def main():
    z = np.load(sys.argv[1])
    W = z["weights"].astype(np.float64)
    L, n, _ = W.shape
    idx = np.arange(n)
    truth = []
    for l in range(L):
        st = ch.atlas_state(z, l)
        D3, D21, K4, K31, K22 = ch.slices_from_state(st)
        truth.append(dict(mu=st.mu, var=st.var, Coff=ok.offdiag(st.C), D3=D3, D21=D21, K4=K4, K22=K22,
                          X211=ok.all_distinct(st.X)))
    for spec in sys.argv[2:]:
        kw = eval(spec)
        if "atlas" in kw.get("k4mode", "") or "force" in kw:
            kw["atlas"] = z
        chn = ch.Chain(W, record=True, **kw)
        out = chn.run()
        print(f"\n{spec}")
        print(f"{'l':>2} {'mu':>7} {'var':>7} {'Coff':>7} {'D3':>7} {'D21':>7} {'K4':>7} {'K22':>7} | {'mse(a)':>8}")
        for l, r in enumerate(out["rec"]):
            T = truth[l]
            e = [rel(r["mu"], T["mu"]), rel(r["var"], T["var"]), rel(ok.offdiag(r["C"]), T["Coff"])]
            e += [rel(r["D3"], T["D3"]), rel(r["D21"], T["D21"]), rel(r["K4"], T["K4"]), rel(r["K22"], T["K22"])] if l else [0] * 4
            mse = float(np.mean((out["means"][l] - z["post_s"][0, l]) ** 2))
            print(f"{l:>2} " + " ".join(f"{v:7.4f}" for v in e) + f" | {mse:8.2e}", flush=True)


if __name__ == "__main__":
    main()
