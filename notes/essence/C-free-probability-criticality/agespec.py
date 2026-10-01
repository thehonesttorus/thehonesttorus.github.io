"""Age spectrum of the content-transfer operator (team C, essence programme).

For every (source s, kink layer k) pair of the exact first-order co-state (costate.py, all pairs, no closure)
records, at width n:
  E      = |S_{s->k}|_F^2                     (D21 energy of the pair's (2,1) slice)
  gamma  = <S, K_k> / <K_k, K_k>             (dilation / scale-mixture amplitude, K = 2 m (x) C + diag(C) (x) m)
  share  = gamma^2 |K|^2 / E                  (share of the pair's slice in the dilation sector)
  tr     = tr(U^T U)/n, PR = tr(U^T U)^2 / tr((U^T U)^2)   (free multiplicative convolution: mean and participation)
and per layer: g_l = tr(M^T M)/n of the one-step mean-gate map M = diag(Phi_l) W_{l+1}, and the Gram matrix of
the dilation-deflated residuals R_{s->k} = S - gamma K across ages (cross-age cosines; coherent vs incoherent sum).

Usage: python agespec.py <set> <mlp index> [out.json]
"""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
FS = os.path.join(HERE, "../../fresh-slate")
sys.path.insert(0, os.path.join(FS, "bench"))
sys.path.insert(0, os.path.join(FS, "breakthrough/costate"))
sys.path.insert(0, os.path.join(FS, "designs/heisenberg"))
import bench
from costate import source_atoms, pair_slice
from hd import Gauss, inject


def run(Ws, K=8, coinc=False):
    Ws = np.asarray(Ws, dtype=np.float64)
    L, n, _ = Ws.shape
    m = np.zeros(n); C = Ws[0].T @ Ws[0]
    live = []
    rows, layers = [], []
    for l in range(L):
        G = Gauss(m, C, K=K)
        Kg = 2 * G.m[:, None] * G.C + (G.m[None, :] * np.diag(G.C)[:, None])
        KK = float(np.sum(Kg * Kg))
        S = np.zeros((n, n)); res = []; amps = []
        for f in live:
            Sf = pair_slice(f["a2"], f["U"], f["V"], f["X"], f["dl"], coinc)
            S += Sf
            E = float(np.sum(Sf * Sf)); Eoff = E - float(np.sum(np.diag(Sf) ** 2))
            gam = float(np.sum(Kg * Sf)) / KK if KK > 0 else 0.0
            R = Sf - gam * Kg
            res.append(R); amps.append(gam)
            UtU = f["U"].T @ f["U"]
            tr = float(np.trace(UtU)) / n; tr2 = float(np.sum(UtU * UtU)) / n
            rows.append(dict(s=f["s"], k=l, age=l - f["s"] - 1, E=E, Eoff=Eoff, gamma=gam,
                             share=gam * gam * KK / E if E > 0 else 0.0, tr=tr, PR=n * tr * tr / tr2))
        lay = dict(layer=l, KK=KK, Etot=float(np.sum(S * S)), Phi2=float(np.mean(G.Phi ** 2)),
                   meanPhi=float(np.mean(G.Phi)), t_rms=float(np.sqrt(np.mean(G.t ** 2))))
        if res:
            Rm = np.array([r.ravel() for r in res])
            Gr = Rm @ Rm.T
            lay["resgram"] = Gr.tolist(); lay["ages"] = [l - f["s"] - 1 for f in live]; lay["amps"] = amps
            Ssum_res = Rm.sum(0)
            lay["Eres_sum"] = float(Ssum_res @ Ssum_res); lay["Eres_inc"] = float(np.trace(Gr))
            ga = sum(amps); lay["Escale_sum"] = ga * ga * KK
            lay["Escale_inc"] = sum(a * a for a in amps) * KK
            del Rm
        layers.append(lay)
        if live:
            D = np.diag(S).copy(); dEa, dC = inject(G, D, S)
        else:
            dEa, dC = np.zeros(n), np.zeros((n, n))
        Ea = G.Ea + dEa
        if l + 1 == L:
            break
        W = Ws[l + 1]
        M = G.Phi[:, None] * W
        lay["g"] = float(np.sum(M * M)) / n
        for f in live:
            for key in ("U", "V", "X"):
                if f[key] is not None:
                    f[key] = f[key] @ M
        a2, R_, Dl, dl = source_atoms(G, K)
        live.append(dict(s=l, a2=a2, dl=dl, U=W.copy(), V=R_ @ W, X=(Dl @ W) if coinc else None))
        m = Ea @ W
        C = W.T @ (G.cov_a() + dC) @ W
        print(f"layer {l} done {time.time()-T0:.0f}s live={len(live)}", flush=True)
    return rows, layers


if __name__ == "__main__":
    T0 = time.time()
    name, i = sys.argv[1], int(sys.argv[2])
    out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(HERE, f"results/agespec_{name}_{i}.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    Sset = bench.load_set(name)
    W = bench.weights(Sset, i)
    rows, layers = run(W)
    json.dump(dict(set=name, mlp=i, rows=rows, layers=layers), open(out, "w"))
    print("wrote", out, f"{time.time()-T0:.0f}s")
