"""G1 measurement (team B): blockwise frozen-frame rank of the memory's lambda term, by dyadic age block, at n = 1024.
Inside FC (region fc.py via ../tests/fc_hooked.py; slices=2, k4mf=True) on bench w1024_d16.
For a source s, its legs at target t are L(t) = [Y_s(t); Z_s(t); T_s(t)] (3n x n rows; columns = target neurons), and
L(t+1) = L(t) G_t with G_t = diag(P_t) W_{t+1}. Its D21 contribution D_s(t) is the readout.
(1) rank law: energy ranks r90/r99 of L(t) (gram of the column space) against age a = t - s.
(2) frozen frame per dyadic block: at block heads a0 in {1,2,4,8} freeze Q = top-r right singular vectors of L(t0),
    r = ceil(c n / a0); afterwards carry only the coefficients L(t0) Q (3n x r, fixed) and the frame F = Q^T G...G (r x n,
    transported at n^2 r per layer). Error of the D21 readout of the source, relative to its exact D_s(t), at every t in
    the block, and the leg error.  Plain truncation (drop the source) has relative error 1 by definition.
usage: python g1_blocks.py MLP c"""
import sys, json, math
import numpy as np
sys.path.insert(0, '/home/user/thehonesttorus.github.io/notes/essence/B-pseudorandomness/tests')
sys.path.insert(0, '/home/user/thehonesttorus.github.io/notes/fresh-slate/breakthrough/region')
sys.path.insert(0, '/home/user/thehonesttorus.github.io/notes/fresh-slate/bench')
import bench
import fc_hooked as F

mlp = int(sys.argv[1]); c = float(sys.argv[2])
S = bench.load_set('w1024_d16'); W = bench.weights(S, mlp); W64 = W.astype(np.float64)
n = W.shape[1]
F.GATES = []
SRC = {1, 2, 3, 5}           # birth layers followed
state = {}                   # s -> dict(head, coeffs (Y,Z,T in Q coords), frame F)
rows = []

def legs(src):
    Zf = src["Z"].astype(np.float64); Y = src["SP"].astype(np.float64) @ Zf
    T = src["Delta"].astype(np.float64) @ Zf if "Delta" in src else np.zeros_like(Zf)
    return Y, Zf, T

def readout(Y, Z, T, w):
    return ((Y * Y) * w[:, None]).T @ Z + 2 * (((Y * Z) * w[:, None]).T @ Y) + (Z * Z).T @ T + 2 * ((Z * T).T @ Z)

def hook(t, sources):
    P, l = F.GATES[-1]
    G = P[:, None] * W64[l + 1]            # transport applied at this iteration
    for s, st in state.items():
        st["F"] = st["F"] @ G
    for src in sources:
        s = src["s"]
        if s not in SRC: continue
        a = t - s
        Y, Z, T = legs(src); w = src["w2"].astype(np.float64)
        Lg = np.vstack([Y, Z, T])
        ev = np.linalg.eigvalsh(Lg.T @ Lg)[::-1]; ev = np.maximum(ev, 0); cs = np.cumsum(ev) / ev.sum()
        r90 = int(np.searchsorted(cs, 0.90) + 1); r99 = int(np.searchsorted(cs, 0.99) + 1)
        Dex = readout(Y, Z, T, w)
        row = dict(s=s, t=t, age=a, r90=r90, r99=r99, pr=float(ev.sum() ** 2 / (ev ** 2).sum()))
        if a in (1, 2, 4, 8) or s not in state:
            r = min(n, math.ceil(c * n / a))
            _, V = np.linalg.eigh(Lg.T @ Lg); Q = V[:, ::-1][:, :r]
            state[s] = dict(head=a, r=r, CY=Y @ Q, CZ=Z @ Q, CT=T @ Q, F=Q.T.copy())
        st = state[s]
        Ya, Za, Ta = st["CY"] @ st["F"], st["CZ"] @ st["F"], st["CT"] @ st["F"]
        Dap = readout(Ya, Za, Ta, w)
        row.update(head=st["head"], r=st["r"],
                   leg_err=float(np.linalg.norm(np.vstack([Ya, Za, Ta]) - Lg) / np.linalg.norm(Lg)),
                   d21_err=float(np.linalg.norm(Dap - Dex) / np.linalg.norm(Dex)))
        rows.append(row)
        print(json.dumps(row), flush=True)

F.HOOK = hook
F.run(W, slices=2, k4mf=True)
json.dump(rows, open(f'g1_blocks_mlp{mlp}_c{c:g}.json', 'w'), indent=1)
