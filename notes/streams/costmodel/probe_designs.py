"""Metered skeletons of one middle-layer transition of each chain design, plus the extra forms the
cost model needs (split-sign weighted Gram, the SymmetryError hazard of the tagged sandwich).

A skeleton runs the exact flopscope op stream of one transition l -> l+1 of a design on random
data of the right shapes (FLOP counts are data-independent), grouped the way an estimator would
group them (one batched family per role, pooled out= buffers), and reports units, flopscope call
count and residual wall time of the second (warm) run. cost.py uses these to validate its
per-layer prices and to set its calls/residual columns.

Run:  OPENBLAS_NUM_THREADS=1 python probe_designs.py designs_measured.json
"""
import json
import sys
import time
import warnings

import numpy as np

warnings.simplefilter("ignore")
import flopscope as flops  # noqa: E402
import flopscope.numpy as fnp  # noqa: E402

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from strassen_v29 import _Strassen  # noqa: E402

N = 1024
H = N // 2
UNIT = float(2 ** 31)
f32 = np.float32
R = np.random.default_rng(1)
A = fnp.asarray
OUT = []


def spd(scale=1.0):
    X = R.standard_normal((N, N)) / np.sqrt(N)
    S = X @ X.T + np.eye(N)
    return (((S + S.T) / 2) * scale).astype(f32)


def gen(scale=1.0):
    return (R.standard_normal((N, N)) * scale).astype(f32)


def vec(lo=0.1, hi=0.9):
    return R.uniform(lo, hi, N).astype(f32)


def measure(name, prep, fn, reps=2, **meta):
    rec = dict(name=name, **meta)
    try:
        with flops.BudgetContext(flop_budget=int(1e18), quiet=True) as ctx:
            st = prep()
            for _ in range(reps):
                f0, n0 = ctx.flops_used, len(ctx.op_log)
                r0 = ctx.residual_wall_time_s
                b0 = ctx.flopscope_backend_time_s
                t0 = time.perf_counter()
                fn(st)
                t1 = time.perf_counter()
            ops = {}
            for r in ctx.op_log[n0:]:
                k = r.op_name
                c, m = ops.get(k, (0, 0))
                ops[k] = (c + int(r.flop_cost), m + 1)
            rec.update(units=(ctx.flops_used - f0) / UNIT, calls=len(ctx.op_log) - n0,
                       residual_ms=1e3 * (ctx.residual_wall_time_s - r0),
                       backend_ms=1e3 * (ctx.flopscope_backend_time_s - b0), wall_ms=1e3 * (t1 - t0),
                       ops={k: [v[0] / UNIT, v[1]] for k, v in ops.items()})
    except Exception as exc:  # noqa: BLE001
        rec["error"] = f"{type(exc).__name__}: {str(exc)[:200]}"
    OUT.append(rec)
    json.dump(dict(n=N, results=OUT), open(sys.argv[1] if len(sys.argv) > 1 else "designs_measured.json", "w"), indent=1)
    msg = rec.get("error") or f"{rec['units']:8.3f} u  calls {rec['calls']:5d}  resid {rec['residual_ms']:7.2f} ms  wall {rec['wall_ms']:6.0f} ms"
    print(f"{name:72s} {msg}", flush=True)
    return rec


# ------------------------------------------------------------------ hazard: tagged sandwich ---------
def hazard(scale, kind):
    def prep():
        X0 = R.standard_normal((N, N))
        S = ((X0 + X0.T) / 2 * scale).astype(f32) if kind == "indef" else spd(scale)
        return dict(W=A((R.standard_normal((N, N)) * np.sqrt(2 / N)).astype(f32)), S=A(S))

    def fn(st):
        St = flops.as_symmetric(st["S"], symmetry=(0, 1))
        fnp.einsum("ji,jk,kl->il", st["W"], St, st["W"])
    return prep, fn


for kind in ("indef", "spd"):
    for scale in (1.0, 0.1, 0.01, 1e-3):
        p, f = hazard(scale, kind)
        measure(f"hazard tagged sandwich, {kind} symmetric S, entry scale {scale}", p, f, hazard=True)


# ------------------------------------------------------------------ split-sign weighted Gram -------
def wgram_prep():
    G = gen(0.03)
    d = R.standard_normal(N).astype(f32)
    return dict(G=A(G), d=A(d), pos=A(np.nonzero(d > 0)[0]), neg=A(np.nonzero(d <= 0)[0]))


def wgram_split(st):
    """G^T diag(d) G for sign-indefinite d: gather the rows of each sign (4 flops/elt), scale by
    sqrt|d|, two aliased Grams over disjoint contraction ranges (together 0.5 u)."""
    G, d = st["G"], st["d"]
    sq = fnp.sqrt(fnp.abs(d))
    Gp = fnp.take(G * sq[:, None], st["pos"], axis=0)
    Gn = fnp.take(G * sq[:, None], st["neg"], axis=0)
    return fnp.einsum("ji,jk->ik", Gp, Gp) - fnp.einsum("ji,jk->ik", Gn, Gn)


def wgram_plain(st):
    G, d = st["G"], st["d"]
    return (G * d[:, None]).T @ G


measure("weighted Gram G^T diag(d) G, plain product", wgram_prep, wgram_plain)
measure("weighted Gram, split-sign aliased (take rows by sign)", wgram_prep, wgram_split)


# ------------------------------------------------------------------ design skeletons ---------------
class Fam:
    """a batched family of b independent (n,n)@(n,n) products: stacked operands in pooled buffers,
    dense (one batched matmul call) or Strassen at `lev` through the copied V29 kernel."""

    def __init__(self, b, lev, shape_l=(N, N), shape_r=(N, N)):
        self.b, self.lev = b, lev
        self.X = fnp.empty((b, 1) + shape_l, dtype=f32); fnp.copyto(self.X, 0.0)
        self.Y = fnp.empty((b, 1) + shape_r, dtype=f32); fnp.copyto(self.Y, 0.0)
        self.O = fnp.empty((b, 1, shape_l[0], shape_r[1]), dtype=f32); fnp.copyto(self.O, 0.0)
        self.sm = _Strassen(f32, b) if lev > 0 else None
        self.shape = (shape_l[0], shape_l[1], shape_r[1])

    def run(self):
        if self.sm is None:
            fnp.matmul(self.X, self.Y, out=self.O)
        else:
            self.sm.mm(self.X, self.Y, self.O, self.sm.level(*self.shape, self.lev))
        return self.O


def load(dst, src):
    fnp.copyto(dst, src)


def skeleton(design, lev, r=0, k=0, cov="sym3", wick_only=False):
    """one middle transition l -> l+1.
    design: 'v' K=2 covariance only; 'ii' closure births, D21/D3 transported; 'iii' ii + rank-r (2,1,1)
    family (r sandwiches); 'iv' ii + k covariance-response modes (k sandwiches + k birth Grams).
    cov: 'einsum' (tagged 3-operand einsum, 1.5 u) or 'sym3' (S W product in a family + 3 output blocks)."""
    n_aleg = 1 if wick_only else 4
    n_bleg = 1 if wick_only else 4
    nsand = r + k            # mode sandwiches (iii, iv)

    def prep():
        st = dict(W=A((R.standard_normal((N, N)) * np.sqrt(2 / N)).astype(f32)), Ca=A(spd()), Coff=A(gen(0.03)),
                  D21=A(gen(0.01)), S21=A(gen(0.01)), Phi=A(vec()), w2=A(vec(0.2, 0.6)), w3=A(vec(-0.3, 0.3)),
                  w5=A(vec(-0.3, 0.3)), D3=A(vec(-0.1, 0.1)), v=A(vec(-0.1, 0.1)),
                  M=[A(spd(0.03)) for _ in range(nsand)], dmode=[A(vec(-1, 1)) for _ in range(k)])
        st["Cat"] = flops.as_symmetric(st["Ca"], symmetry=(0, 1))
        st["WW"] = fnp.multiply(st["W"], st["W"])
        st["XP"] = fnp.empty((N, N), dtype=f32); st["XW2"] = fnp.empty((N, N), dtype=f32)
        st["L"] = fnp.empty((N, N), dtype=f32); st["T"] = fnp.empty((N, N), dtype=f32)
        st["XPT"] = fnp.empty((N, N), dtype=f32)
        # families: a-legs (+ C_a W for the covariance, + modes M (Phi W)), b-legs, slices, final, mode 3-block outputs
        n_a = (n_aleg if design != "v" else 0) + (1 if cov == "sym3" else 0) + nsand
        st["FA"] = Fam(n_a, lev) if n_a else None
        st["FB"] = Fam(n_bleg, lev) if design != "v" else None
        st["FS"] = Fam(2, lev) if design != "v" else None
        st["FF"] = Fam(1, lev) if design != "v" else None
        n3 = (1 if cov == "sym3" else 0) + nsand
        st["F3"] = Fam(3 * n3, max(lev - 1, 0), (H, N), (N, H)) if n3 else None
        st["FK"] = Fam(k, lev) if k else None       # birth Grams of the modes: (d_m * G_C)^T G_C
        return st

    def fn(st):
        W, Phi, w2 = st["W"], st["Phi"], st["w2"]
        fnp.multiply(W, Phi[:, None], out=st["XP"])
        fnp.multiply(W, w2[:, None], out=st["XW2"])
        i = 0
        FA = st["FA"]
        if design != "v":
            load(FA.X[0, 0], st["Coff"]); load(FA.Y[0, 0], st["XP"]); i = 1
            if not wick_only:
                load(FA.X[1, 0], st["Coff"]); load(FA.Y[1, 0], st["XW2"])
                load(FA.X[2, 0], st["D21"].T); load(FA.Y[2, 0], st["XW2"])
                load(FA.X[3, 0], st["D21"]); load(FA.Y[3, 0], st["XP"]); i = 4
        if cov == "sym3":
            load(FA.X[i, 0], st["Ca"]); load(FA.Y[i, 0], W); icov = i; i += 1
        for m in range(nsand):
            load(FA.X[i + m, 0], st["M"][m]); load(FA.Y[i + m, 0], st["XP"])
        if FA is not None:
            FA.run()
        if cov == "einsum":
            fnp.einsum("ji,jk,kl->il", W, st["Cat"], W)
        if design != "v":
            G = FA.O
            GC = G[0, 0]
            FB, L = st["FB"], st["L"]
            # Lt for the b-leg products: sums of vertex-weighted W * G1 (n^2 ops), 2-3 terms per b-leg
            for j in range(n_bleg):
                Lt = FB.Y[j, 0]
                fnp.multiply(st["XW2"], G[j, 0], out=Lt)
                fnp.multiply(W, G[min(j + 1, n_aleg - 1), 0], out=st["T"])
                fnp.multiply(st["T"], st["w3"][:, None], out=st["T"])
                fnp.add(Lt, st["T"], out=Lt)
            fnp.multiply(st["XP"], st["v"][None, :], out=st["T"]); fnp.add(FB.Y[0, 0], st["T"], out=FB.Y[0, 0])  # B6 r=1 core
            load(FB.X[0, 0], st["Coff"])
            if not wick_only:
                load(FB.X[1, 0], st["Coff"]); load(FB.X[2, 0], st["D21"]); load(FB.X[3, 0], st["D21"].T)
            FB.run()
            FS = st["FS"]
            load(FS.X[0, 0], st["S21"]); load(FS.Y[0, 0], W)
            load(FS.X[1, 0], st["S21"].T); load(FS.Y[1, 0], st["WW"])
            FS.run()
            # type-A and b-leg outputs, slices, D3 -> L (n^2 ops)
            fnp.multiply(GC, GC, out=L); fnp.multiply(L, w2[:, None], out=L)
            for j in range(1, n_aleg):
                fnp.multiply(GC, G[j, 0], out=st["T"]); fnp.multiply(st["T"], st["w3"][:, None], out=st["T"]); fnp.add(L, st["T"], out=L)
            for j in range(n_bleg):
                fnp.multiply(FB.O[j, 0], Phi[:, None], out=st["T"]); fnp.add(L, st["T"], out=L)
            fnp.add(L, FS.O[1, 0], out=L)
            fnp.multiply(W, FS.O[0, 0], out=st["T"]); fnp.multiply(st["T"], 2.0, out=st["T"]); fnp.add(L, st["T"], out=L)
            fnp.multiply(st["WW"], st["D3"][:, None], out=st["T"]); fnp.add(L, st["T"], out=L)
            FF = st["FF"]
            load(FF.X[0, 0], L.T); load(FF.Y[0, 0], W)
            FF.run()
            fnp.diagonal(FF.O[0, 0])
        F3 = st["F3"]
        if F3 is not None:
            # symmetric outputs X^T T from 3 of the 4 blocks: (X^T)[:h] T[:, :h], (X^T)[:h] T[:, h:], (X^T)[h:] T[:, h:]
            j = 0
            srcs = ([(W, icov)] if cov == "sym3" else []) + [(st["XP"], i + m) for m in range(nsand)]
            for X, slot in srcs:
                fnp.copyto(st["XPT"], X.T)
                T = FA.O[slot, 0]
                load(F3.X[j, 0], st["XPT"][:H]); load(F3.Y[j, 0], T[:, :H])
                load(F3.X[j + 1, 0], st["XPT"][:H]); load(F3.Y[j + 1, 0], T[:, H:])
                load(F3.X[j + 2, 0], st["XPT"][H:]); load(F3.Y[j + 2, 0], T[:, H:])
                j += 3
            F3.run()
        if k:
            FK = st["FK"]
            for m in range(k):
                fnp.multiply(FA.O[0, 0], st["dmode"][m][:, None], out=st["T"])
                load(FK.X[m, 0], st["T"].T); load(FK.Y[m, 0], FA.O[0, 0])
            FK.run()
        # nonlinearity on n x n (term program stand-in: ~40 elementwise ops) and Wick weights
        for _ in range(40):
            fnp.multiply(st["Coff"], st["Coff"], out=st["T"])
    return prep, fn


for lev in (0, 3, 4, 5):
    for cov in ("einsum", "sym3"):
        measure(f"design v (K=2 covariance) cov={cov} L{lev}", *skeleton("v", lev, cov=cov), design="v", lev=lev, cov=cov)
for lev in (0, 3, 4, 5):
    measure(f"design ii (closure births, 4 leg types) L{lev}", *skeleton("ii", lev), design="ii", lev=lev, cov="sym3")
    measure(f"design ii-wick (leading Wick leg only) L{lev}", *skeleton("ii", lev, wick_only=True), design="ii-wick", lev=lev, cov="sym3")
for r in (1, 4, 8):
    for lev in ((0, 5) if r < 8 else (0,)):      # r = 8 at L5 was OOM-killed at 5.6 GB RSS (unshared pools)
        measure(f"design iii r={r} L{lev}", *skeleton("iii", lev, r=r), design="iii", lev=lev, r=r, cov="sym3")
for k in (4, 8):
    for lev in ((0, 5) if k == 4 else (0,)):      # k = 8 at L5 exceeded this box's free memory next to other jobs
        measure(f"design iv k={k} L{lev}", *skeleton("iv", lev, k=k), design="iv", lev=lev, k=k, cov="sym3")

json.dump(dict(n=N, results=OUT), open(sys.argv[1] if len(sys.argv) > 1 else "designs_measured.json", "w"), indent=1)
print("done", flush=True)
