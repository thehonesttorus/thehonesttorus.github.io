"""Measured flopscope 0.12.1 prices of the operations a closure-based K=3 chain needs at n = 1024.

Every form is run twice inside one in-process BudgetContext (first run = warm-up of the numpy kernels
and flopscope caches, second run measured): billed FLOPs (in units of one 1024^3 matmul = 2^31 FLOPs),
number of flopscope calls (op_log records), backend / overhead / residual wall time of the measured run,
and the per-op breakdown. Output: JSON (argv[1]) consumed by cost.py.

Run:  OPENBLAS_NUM_THREADS=1 python probe_ops.py ops_measured.json
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
UNIT = float(2 ** 31)
BUDGET = int(1e18)
f32, f64 = np.float32, np.float64
R = np.random.default_rng(0)


def he(n=N, dt=f32):
    return (R.standard_normal((n, n)) * np.sqrt(2.0 / n)).astype(dt)


def spd(n=N, dt=f32):
    X = R.standard_normal((n, n)) / np.sqrt(n)
    S = X @ X.T + np.eye(n)
    return ((S + S.T) / 2).astype(dt)


def symm(n=N, dt=f32):
    X = R.standard_normal((n, n))
    return ((X + X.T) / 2).astype(dt)


Wn = he()
Sn = spd()
Gn = R.standard_normal((N, N)).astype(f32)
Hn = R.standard_normal((N, N)).astype(f32)
pn = R.uniform(0.2, 0.9, N).astype(f32)
wn = R.standard_normal(N).astype(f32)
RESULTS = []


def measure(name, prep, fn, note="", reps=2):
    """prep() -> tuple of flopscope arrays (built inside the context, not timed);
    fn(*args) is run `reps` times; the last run is the measurement."""
    rec = dict(name=name, note=note)
    try:
        with flops.BudgetContext(flop_budget=BUDGET, quiet=True) as ctx:
            args = prep()
            for rep in range(reps):
                f0, n0 = ctx.flops_used, len(ctx.op_log)
                b0, o0, r0 = ctx.flopscope_backend_time_s, ctx.flopscope_overhead_time_s, ctx.residual_wall_time_s
                t0 = time.perf_counter()
                out = fn(*args)
                t1 = time.perf_counter()
                f1, n1 = ctx.flops_used, len(ctx.op_log)
                b1, o1, r1 = ctx.flopscope_backend_time_s, ctx.flopscope_overhead_time_s, ctx.residual_wall_time_s
            ops = {}
            for r in ctx.op_log[n0:n1]:
                k = f"{r.op_name}[{r.resolved_dtype}]"
                c, m = ops.get(k, (0, 0))
                ops[k] = (c + int(r.flop_cost), m + 1)
            rec.update(flops=int(f1 - f0), units=(f1 - f0) / UNIT, calls=int(n1 - n0),
                       wall_ms=1e3 * (t1 - t0), backend_ms=1e3 * (b1 - b0), overhead_ms=1e3 * (o1 - o0),
                       residual_ms=1e3 * (r1 - r0), ops={k: [v[0] / UNIT, v[1]] for k, v in ops.items()})
            if isinstance(out, dict):
                rec.update(out)
    except Exception as exc:  # noqa: BLE001
        rec.update(error=f"{type(exc).__name__}: {str(exc)[:200]}")
    RESULTS.append(rec)
    msg = rec.get("error") or f"{rec['units']:.4f} u  calls {rec['calls']}  resid {rec['residual_ms']:.2f} ms  wall {rec['wall_ms']:.0f} ms"
    print(f"{name:70s} {msg}", flush=True)
    return rec


A = fnp.asarray

# ---------------------------------------------------------------- 1. plain products -----------------
measure("mm_f32: (n,n)@(n,n) float32", lambda: (A(Wn), A(Gn)), lambda a, b: a @ b)
measure("mm_f64: (n,n)@(n,n) float64", lambda: (A(Wn.astype(f64)), A(Gn.astype(f64))), lambda a, b: a @ b)
measure("mm_mixed: f32 @ f64 (promotes, billed f64)", lambda: (A(Wn), A(Gn.astype(f64))), lambda a, b: a @ b)


def _mm_out(a, b, o):
    fnp.matmul(a, b, out=o)


measure("mm_out_pooled: matmul(out=preallocated)", lambda: (A(Wn), A(Gn), fnp.empty((N, N), dtype=f32)), _mm_out)
measure("mm_transposed_view: W.T @ G (strided left operand)", lambda: (A(Wn), A(Gn)), lambda a, b: a.T @ b)
for r in (16, 32, 64, 128, 256, 384):
    measure(f"thin_mm r={r}: (n,n)@(n,{r})", lambda r=r: (A(Wn), A(Gn[:, :r].copy())), lambda a, b: a @ b)
measure("matvec: (n,n)@(n,)", lambda: (A(Wn), A(wn)), lambda a, b: a @ b)

# ---------------------------------------------------------------- 2. batched products ---------------
for k in (2, 4, 8, 16):
    measure(f"batched k={k}: (k,n,n)@(k,n,n)",
            lambda k=k: (A(np.stack([Wn] * k)), A(np.stack([Gn] * k))), lambda a, b: a @ b)
    measure(f"bcast k={k}: (n,n)@(k,n,n)",
            lambda k=k: (A(Wn), A(np.stack([Gn] * k))), lambda a, b: a @ b)
    measure(f"concat k={k}: (n,n)@(n,k*n) (one wide product)",
            lambda k=k: (A(Wn), A(np.concatenate([Gn] * k, 1))), lambda a, b: a @ b)
measure("loop k=8: 8 separate (n,n)@(n,n) calls", lambda: (A(Wn), A(Gn)),
        lambda a, b: [a @ b for _ in range(8)])

# ---------------------------------------------------------------- 3. Gram / symmetric ---------------
measure("gram_alias: einsum('ji,jk->ik', X, X) same object", lambda: (A(Gn),), lambda x: fnp.einsum("ji,jk->ik", x, x))
measure("gram_inner: inner(X, X) = X X^T same object", lambda: (A(Gn),), lambda x: fnp.inner(x, x))
measure("gram_matmul_T: X.T @ X (view is another object)", lambda: (A(Gn),), lambda x: x.T @ x)
measure("as_symmetric (n,n) f32", lambda: (A(Sn),), lambda s: flops.as_symmetric(s, symmetry=(0, 1)))
measure("symmetrize canonical-copy", lambda: (A(Sn),), lambda s: flops.symmetrize(s, mode="canonical-copy") if hasattr(flops, "symmetrize") else fnp.symmetrize(s, mode="canonical-copy"))
measure("weighted alias einsum('ji,j,jk->ik', G, d, G)", lambda: (A(Gn), A(wn)), lambda g, d: fnp.einsum("ji,j,jk->ik", g, d, g))


# ---------------------------------------------------------------- 4. sandwiches W^T diag(p) S diag(p) W
def sw_plain(W, p, S):
    X = W * p[:, None]
    T = S @ X
    return X.T @ T


def sw_einsum3(W, p, S):
    X = W * p[:, None]
    St = flops.as_symmetric(S, symmetry=(0, 1))
    return fnp.einsum("ji,jk,kl->il", X, St, X)


def sw_einsum3_pretag(W, p, St):
    X = W * p[:, None]
    return fnp.einsum("ji,jk,kl->il", X, St, X)


def sw_v25(W, p, St):
    X = W * p[:, None]
    return fnp.einsum("ij,ia,jb->ab", St, X, X)


def sw_untagged(W, p, S):
    X = W * p[:, None]
    return fnp.einsum("ji,jk,kl->il", X, S, X)


def sw_chol(W, p, S):
    X = W * p[:, None]
    L = fnp.linalg.cholesky(S)
    Y = L.T @ X
    return fnp.einsum("ji,jk->ik", Y, Y)


def sw_factor(W, p, U):
    """S given as a factor U (S = U U^T, e.g. a Gram the chain already holds)"""
    X = W * p[:, None]
    Y = U.T @ X
    return fnp.einsum("ji,jk->ik", Y, Y)


def sw_diag_only(W, p, S):
    X = W * p[:, None]
    return fnp.sum(X * (S @ X), axis=0)


measure("sandwich_plain: X=diag(p)W; X.T @ (S @ X)", lambda: (A(Wn), A(pn), A(Sn)), sw_plain)
measure("sandwich_einsum3: as_symmetric(S) + einsum('ji,jk,kl->il',X,S,X)", lambda: (A(Wn), A(pn), A(Sn)), sw_einsum3)
measure("sandwich_einsum3_pretagged (S already tagged)", lambda: (A(Wn), A(pn), flops.as_symmetric(A(Sn), symmetry=(0, 1))), sw_einsum3_pretag)
measure("sandwich_v25: einsum('ij,ia,jb->ab', S_tag, X, X)", lambda: (A(Wn), A(pn), flops.as_symmetric(A(Sn), symmetry=(0, 1))), sw_v25)
measure("sandwich_untagged_einsum3 (S not tagged)", lambda: (A(Wn), A(pn), A(Sn)), sw_untagged)
measure("sandwich_cholesky: L=chol(S); Y=L^T X; Y^T Y aliased", lambda: (A(Wn), A(pn), A(Sn)), sw_chol)
measure("sandwich_factor: S=UU^T given; Y=U^T X; Y^T Y aliased", lambda: (A(Wn), A(pn), A(np.linalg.cholesky(Sn.astype(f64)).astype(f32))), sw_factor)
measure("sandwich_diag_only: sum(X*(S@X), 0)", lambda: (A(Wn), A(pn), A(Sn)), sw_diag_only)
measure("sandwich_einsum3_f64", lambda: (A(Wn.astype(f64)), A(pn.astype(f64)), A(Sn.astype(f64))), sw_einsum3)


# ---------------------------------------------------------------- 5. Strassen-Winograd ---------------
def _strassen_prod(lev, k=1, hub=False, fresh=False):
    """One product (or a batch of k) through the copied V29 kernel at `lev` levels.
    Returns prep, fn. The pools are allocated in prep (first-touch billed there, as in V29's
    persistent pools) and the measured second run reuses them."""
    def prep():
        sm = _Strassen(f32, k)
        X = A(np.stack([Wn] * k)[:, None])        # (k,1,n,n)
        Y = A(np.stack([Gn] * k)[:, None])
        if hub:
            out = fnp.empty((1, N, N), dtype=f32)
        else:
            out = fnp.empty((k, 1, N, N), dtype=f32)
        return sm, X, Y, out

    def fn(sm, X, Y, out):
        if hub:
            sm.hub(X, Y, out, sm.level(N, N, N, lev))
        else:
            sm.mm(X, Y, out, sm.level(N, N, N, lev))
        return {}
    return prep, fn


ref64 = Wn.astype(f64) @ Gn.astype(f64)
for lev in (0, 1, 2, 3, 4, 5):
    prep, fn = _strassen_prod(lev)
    rec = measure(f"strassen L{lev}: one (n,n)@(n,n)", prep, fn)
    # accuracy of the level (separate un-metered check in plain numpy semantics through flopscope)
    try:
        with flops.BudgetContext(flop_budget=BUDGET, quiet=True):
            sm, X, Y, out = prep()
            fn(sm, X, Y, out)
            o = np.asarray(out)[0, 0].astype(f64)
        rec["rel_err_vs_f64"] = float(np.linalg.norm(o - ref64) / np.linalg.norm(ref64))
        print(f"    rel err L{lev}: {rec['rel_err_vs_f64']:.2e}")
    except Exception as exc:  # noqa: BLE001
        rec["rel_err_vs_f64"] = f"ERR {exc}"
for lev in (3, 5):
    prep, fn = _strassen_prod(lev, k=4)
    measure(f"strassen L{lev}: batch of 4 products", prep, fn)
    prep, fn = _strassen_prod(lev, k=4, hub=True)
    measure(f"strassen hub L{lev}: sum_k X_k Y_k^T, k=4", prep, fn)


def strassen_sandwich(lev):
    def prep():
        sm = _Strassen(f32, 1)
        X = fnp.multiply(A(Wn), A(pn)[:, None])
        XT = fnp.empty((N, N), dtype=f32)
        fnp.copyto(XT, X.T)
        return sm, A(Sn), X, XT, fnp.empty((1, 1, N, N), dtype=f32), fnp.empty((1, 1, N, N), dtype=f32)

    def fn(sm, S, X, XT, T, out):
        l = sm.level(N, N, N, lev)
        sm.mm(S[None, None], X[None, None], T, l)       # T = S X
        sm.mm(XT[None, None], T, out, l)                # X^T T
        return {}
    return prep, fn


for lev in (3, 4, 5):
    prep, fn = strassen_sandwich(lev)
    measure(f"sandwich_strassen L{lev}: two Strassen products", prep, fn)


def strassen_sym3(lev):
    """V29 _sym_product idea: T = S X (Strassen), then only blocks 11, 12, 22 of X^T T
    (a batch-3 family of (n/2, n) @ (n, n/2) products), 21 = 12^T."""
    h = N // 2

    def prep():
        sm = _Strassen(f32, 3)
        X = fnp.multiply(A(Wn), A(pn)[:, None])
        XT = fnp.empty((N, N), dtype=f32)
        fnp.copyto(XT, X.T)
        return sm, A(Sn), X, XT, fnp.empty((1, 1, N, N), dtype=f32), fnp.empty((3, 1, h, h), dtype=f32), fnp.empty((3, 1, h, N), dtype=f32), fnp.empty((3, 1, N, h), dtype=f32), fnp.empty((N, N), dtype=f32)

    def fn(sm, S, X, XT, T, out3, Lb, Rb, out):
        l = sm.level(N, N, N, lev)
        sm.mm(S[None, None], X[None, None], T, l)
        T2 = T[0, 0]
        fnp.copyto(Lb[0, 0], XT[:h]); fnp.copyto(Lb[1, 0], XT[:h]); fnp.copyto(Lb[2, 0], XT[h:])
        fnp.copyto(Rb[0, 0], T2[:, :h]); fnp.copyto(Rb[1, 0], T2[:, h:]); fnp.copyto(Rb[2, 0], T2[:, h:])
        sm.mm(Lb, Rb, out3, sm.level(h, N, h, lev))
        fnp.copyto(out[:h, :h], out3[0, 0]); fnp.copyto(out[:h, h:], out3[1, 0])
        fnp.copyto(out[h:, :h], out3[1, 0].T); fnp.copyto(out[h:, h:], out3[2, 0])
        return {}
    return prep, fn


for lev in (4, 5):
    prep, fn = strassen_sym3(lev)
    measure(f"sandwich_strassen_sym3 L{lev}: S X then 3 of 4 output blocks", prep, fn)


def strassen_factor_gram(lev):
    """S = U U^T given: Y = U^T X by Strassen, then aliased Gram Y^T Y (0.5 u)."""
    def prep():
        sm = _Strassen(f32, 1)
        X = fnp.multiply(A(Wn), A(pn)[:, None])
        U = np.linalg.cholesky(Sn.astype(f64)).astype(f32)
        UT = A(np.ascontiguousarray(U.T))
        return sm, UT, X, fnp.empty((1, 1, N, N), dtype=f32)

    def fn(sm, UT, X, Y):
        sm.mm(UT[None, None], X[None, None], Y, sm.level(N, N, N, lev))
        Y2 = Y[0, 0]
        return {"_": fnp.einsum("ji,jk->ik", Y2, Y2)}
    return prep, fn


for lev in (5,):
    prep, fn = strassen_factor_gram(lev)
    measure(f"sandwich_factor_strassen L{lev}: Y=U^T X Strassen + aliased Gram", prep, fn)


# ---------------------------------------------------------------- 6. Hadamard-then-contract ---------
measure("hc_a: ((H*G)*w[None,:]) @ W", lambda: (A(Hn), A(Gn), A(wn), A(Wn)), lambda h, g, w, W: ((h * g) * w[None, :]) @ W)
measure("hc_b: (H*G) @ (w[:,None]*W)", lambda: (A(Hn), A(Gn), A(wn), A(Wn)), lambda h, g, w, W: (h * g) @ (w[:, None] * W))
measure("hc_einsum4: einsum('ij,ij,j,jk->ik', H, G, w, W)", lambda: (A(Hn), A(Gn), A(wn), A(Wn)), lambda h, g, w, W: fnp.einsum("ij,ij,j,jk->ik", h, g, w, W))
measure("hc_einsum3: einsum('ij,ij,jk->ik', H, G, wW)", lambda: (A(Hn), A(Gn), A(Wn)), lambda h, g, W: fnp.einsum("ij,ij,jk->ik", h, g, W))
measure("hc_typeA: (G1*G2*z[:,None]).T @ W  (center on b-index)", lambda: (A(Hn), A(Gn), A(wn), A(Wn)), lambda h, g, z, W: (h * g * z[:, None]).T @ W)
measure("hc_typeA_einsum: einsum('ka,ka,k,kb->ab', G1, G2, z, W)", lambda: (A(Hn), A(Gn), A(wn), A(Wn)), lambda h, g, z, W: fnp.einsum("ka,ka,k,kb->ab", h, g, z, W))
measure("self_hadamard G*G (same object)", lambda: (A(Gn),), lambda g: g * g)
measure("square(G)", lambda: (A(Gn),), lambda g: fnp.square(g))


def hc_strassen(lev):
    def prep():
        sm = _Strassen(f32, 1)
        return sm, A(Hn), A(Gn), A(wn), A(Wn), fnp.empty((N, N), dtype=f32), fnp.empty((1, 1, N, N), dtype=f32)

    def fn(sm, h, g, z, W, L, out):
        fnp.multiply(h, g, out=L)
        fnp.multiply(L, z[:, None], out=L)
        sm.mm(fnp.swapaxes(L, 0, 1)[None, None], W[None, None], out, sm.level(N, N, N, lev))
        return {}
    return prep, fn


prep, fn = hc_strassen(5)
measure("hc_typeA_strassen L5: L=(G1*G2)*z (out=), L^T @ W Strassen", prep, fn)


# ---------------------------------------------------------------- 7. rank-r symmetric family transport
def fam_loop_tagged(W, p, M):
    X = W * p[:, None]
    out = []
    for m in range(M.shape[0]):
        St = flops.as_symmetric(M[m], symmetry=(0, 1))
        out.append(fnp.einsum("ji,jk,kl->il", X, St, X))
    return {}


def fam_batched_matmul(W, p, M):
    X = W * p[:, None]
    T = fnp.matmul(M, X)            # (r,n,n)
    return {"_": fnp.matmul(X.T, T)}


def fam_batched_einsum_tagged(W, p, Mt):
    X = W * p[:, None]
    return {"_": fnp.einsum("ji,mjk,kl->mil", X, Mt, X)}


def fam_shared_basis(W, p, U, c):
    """M_m = U diag(c_m) U^T with a shared q-column basis U: transport U once, the r cores stay q x q."""
    X = W * p[:, None]
    Ut = X.T @ U              # n x q   (W^T diag(p) U)
    d = fnp.einsum("aq,mq,aq->ma", Ut, c, Ut)   # diag of each transported mode (n per mode)
    return {"_": d}


for r in (1, 4, 8, 16, 32):
    Mn = np.stack([symm() for _ in range(r)])
    measure(f"family_loop_tagged r={r}: r x (as_symmetric + einsum3)", lambda Mn=Mn: (A(Wn), A(pn), A(Mn)), fam_loop_tagged)
    measure(f"family_batched_matmul r={r}: M@X then X^T@T", lambda Mn=Mn: (A(Wn), A(pn), A(Mn)), fam_batched_matmul)
    measure(f"family_batched_einsum_tagged r={r}: einsum('ji,mjk,kl->mil') M tagged (1,2)",
            lambda Mn=Mn: (A(Wn), A(pn), flops.as_symmetric(A(Mn), symmetry=(1, 2))), fam_batched_einsum_tagged)
    for q in (64, 256):
        measure(f"family_shared_basis r={r} q={q}: transport U (n x q) + r diag cores",
                lambda q=q, r=r: (A(Wn), A(pn), A(R.standard_normal((N, q)).astype(f32)), A(R.standard_normal((r, q)).astype(f32))), fam_shared_basis)

# ---------------------------------------------------------------- 8. Gaussian weights / elementwise --
v = R.standard_normal(N).astype(f32)
measure("norm.cdf (n,) f32 in", lambda: (A(v),), lambda x: flops.stats.norm.cdf(x))
measure("norm.pdf (n,) f32 in", lambda: (A(v),), lambda x: flops.stats.norm.pdf(x))
measure("norm.cdf (n,) + astype f32", lambda: (A(v),), lambda x: flops.stats.norm.cdf(x).astype(f32))
measure("exp (n,) f32", lambda: (A(v),), lambda x: fnp.exp(x))
measure("pdf via exp: exp(-x*x/2)*c (n,) f32", lambda: (A(v),), lambda x: fnp.exp(x * x * (-0.5)) * 0.3989422804014327)
measure("norm.cdf (n,n)", lambda: (A(Gn),), lambda x: flops.stats.norm.cdf(x))
measure("norm.pdf (n,n)", lambda: (A(Gn),), lambda x: flops.stats.norm.pdf(x))
measure("exp (n,n) f32", lambda: (A(Gn),), lambda x: fnp.exp(x))
measure("arcsin (n,n) f32", lambda: (A(np.clip(Gn, -0.99, 0.99)),), lambda x: fnp.arcsin(x))
measure("sqrt (n,n) f32", lambda: (A(np.abs(Gn)),), lambda x: fnp.sqrt(x))
measure("multiply (n,n)", lambda: (A(Gn), A(Hn)), lambda a, b: a * b)
measure("multiply out= (n,n)", lambda: (A(Gn), A(Hn), fnp.empty((N, N), dtype=f32)), lambda a, b, o: fnp.multiply(a, b, out=o))
measure("row-scale W*p[:,None]", lambda: (A(Wn), A(pn)), lambda W, p: W * p[:, None])
measure("fill_diagonal(X, 0)", lambda: (fnp.empty((N, N), dtype=f32),), lambda x: fnp.fill_diagonal(x, 0.0))
measure("sum(X, axis=0)", lambda: (A(Gn),), lambda x: fnp.sum(x, axis=0))
measure("copyto (n,n)", lambda: (fnp.empty((N, N), dtype=f32), A(Gn)), lambda o, x: fnp.copyto(o, x))
measure("astype f32->f64 (n,n)", lambda: (A(Gn),), lambda x: x.astype(f64))
measure("qr (n, 64) reduced", lambda: (A(Gn[:, :64].copy()),), lambda x: fnp.linalg.qr(x))
measure("qr (n, 384) reduced", lambda: (A(Gn[:, :384].copy()),), lambda x: fnp.linalg.qr(x))


def wick_weights(mu, var):
    """per-layer Gaussian ReLU weights a K=3 closure needs: Phi, phi, w2..w5 (n-vectors)."""
    sig = fnp.sqrt(var)
    al = mu / sig
    Phi = flops.stats.norm.cdf(al).astype(f32)
    phi = flops.stats.norm.pdf(al).astype(f32)
    isg = 1.0 / sig
    w2 = phi * isg
    w3 = -(al * w2) * isg
    w4 = (al * al - 1.0) * w2 * isg * isg
    w5 = -(al * al * al - 3.0 * al) * w2 * isg * isg * isg
    m1 = mu * Phi + sig * phi
    return {}


measure("wick_weights per layer (Phi, phi, w2..w5, mean)", lambda: (A(v), A(np.abs(v) + 0.5)), wick_weights)

json.dump(dict(n=N, unit_flops=UNIT, flopscope=flops.__version__, results=RESULTS),
          open(sys.argv[1] if len(sys.argv) > 1 else "ops_measured.json", "w"), indent=1)
print("done", flush=True)
