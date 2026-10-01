"""Measure the porting kit's primitives under the in-process flopscope 0.12.1 meter: units (2^31 FLOPs),
flopscope calls, residual ms and backend ms per call of the primitive (3 warm repetitions, last kept).
usage: OPENBLAS_NUM_THREADS=2 python probe_kit.py OUT.json [filter-substring]"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import flopscope as flops
import flopscope.numpy as fnp
import kit

U = float(2 ** 31)
f32, f64 = fnp.float32, fnp.float64
n = 1024
rng = np.random.default_rng(0)


def arr(shape, dtype=np.float32, spd=False, scale=1.0):
    x = rng.standard_normal(shape) * scale
    if spd:
        m = shape[-1]
        x = x @ np.swapaxes(x, -1, -2) / m + np.eye(m)
    return fnp.asarray(x.astype(dtype))


with flops.BudgetContext(flop_budget=10 ** 16, quiet=True):
    A = arr((n, n)); B = arr((n, n)); A64 = arr((n, n), np.float64)
    SPD = arr((n, n), spd=True); SPD64 = arr((n, n), np.float64, spd=True)
    R64c = arr((n, 64)); R256c = arr((n, 256)); V = arr((n,)); Vn = arr((n, 1))
    B8 = arr((8, n, n)); C8 = arr((8, n, n))
    S8 = arr((4096, 8, 8)); T8 = arr((4096, 8, 8)); SPD8 = arr((4096, 8, 8), spd=True)
    S32 = arr((1024, 32, 32)); T32 = arr((1024, 32, 32)); SPD32 = arr((1024, 32, 32), spd=True)
    SPD128 = arr((64, 128, 128), spd=True)
    OUT = fnp.empty((n, n), dtype=f32); fnp.copyto(OUT, 0.0)
    OUT2 = fnp.empty((n, n), dtype=f32); fnp.copyto(OUT2, 0.0)
    OUTV = fnp.empty((n,), dtype=f32); fnp.copyto(OUTV, 0.0)
    IDX = fnp.asarray(rng.integers(0, n, n).astype(np.int64))
    SW = kit.Strassen(f32, 8)
    X1 = A[None, None]; Y1 = B[None, None]
    OS = fnp.empty((1, 1, n, n), dtype=f32); fnp.copyto(OS, 0.0)
    Y4 = arr((4, 1, n, n)); OS4 = fnp.empty((4, 1, n, n), dtype=f32); fnp.copyto(OS4, 0.0)
    GX = SW.gridc(A, X1); GY = SW.gridc(B, Y1); GO = SW.gridc(OS, OS)
    GY4 = SW.gridc(Y4, Y4); GO4 = SW.gridc(OS4, OS4)
    GH4 = SW.gridc(Y4, Y4, hubT=True)
    OH = fnp.empty((1, n, n), dtype=f32); fnp.copyto(OH, 0.0); GOH = SW.gridc(OH, OH[None])

CASES = [
    ("matmul (n,n)@(n,n) f32", lambda: A @ B),
    ("matmul out= pooled", lambda: fnp.matmul(A, B, out=OUT)),
    ("matmul f64", lambda: A64 @ A64),
    ("matmul f32 @ f64 (promotes)", lambda: A @ A64),
    ("matmul (n,n)@(n,64)", lambda: A @ R64c),
    ("matmul (n,n)@(n,256)", lambda: A @ R256c),
    ("matvec (n,n)@(n,)", lambda: A @ V),
    ("batched (8,n,n)@(8,n,n)", lambda: fnp.matmul(B8, C8)),
    ("batched small (4096,8,8)@(4096,8,8)", lambda: fnp.matmul(S8, T8)),
    ("batched small (1024,32,32)@(1024,32,32)", lambda: fnp.matmul(S32, T32)),
    ("A.T @ A (view: full price)", lambda: A.T @ A),
    ("gram_same einsum('ji,jk->ik',A,A)", lambda: kit.gram_same(A)),
    ("sandwich W^T C W (2 dense)", lambda: kit.sandwich(A, B, OUT, OUT2)),
    ("diag_sandwich", lambda: kit.diag_sandwich(A, B, OUT, OUTV)),
    ("tagged einsum sandwich (SPD C)", lambda: fnp.einsum("ij,ia,jb->ab", flops.as_symmetric(SPD, symmetry=(0, 1)), A, A)),
    ("Strassen L1 one product", lambda: SW.mm(X1, Y1, OS, 1, XG=GX, YG=GY, OG=GO)),
    ("Strassen L2 one product", lambda: SW.mm(X1, Y1, OS, 2, XG=GX, YG=GY, OG=GO)),
    ("Strassen L3 one product", lambda: SW.mm(X1, Y1, OS, 3, XG=GX, YG=GY, OG=GO)),
    ("Strassen L4 one product", lambda: SW.mm(X1, Y1, OS, 4, XG=GX, YG=GY, OG=GO)),
    ("Strassen L5 one product", lambda: SW.mm(X1, Y1, OS, 5, XG=GX, YG=GY, OG=GO)),
    ("Strassen L5 batch 4 (shared left)", lambda: SW.mm(X1, Y4, OS4, 5, XG=GX, YG=GY4, OG=GO4)),
    ("Strassen L3 batch 4 (shared left)", lambda: SW.mm(X1, Y4, OS4, 3, XG=GX, YG=GY4, OG=GO4)),
    ("Strassen hub L5 sum_4 X_k Y_k^T", lambda: SW.hub(Y4, Y4, OH, 5, XG=GY4, YG=GH4, OG=GOH)),
    ("add (n,n) out=", lambda: fnp.add(A, B, out=OUT)),
    ("multiply (n,n) by row vector out=", lambda: fnp.multiply(A, V[None, :], out=OUT)),
    ("add (n,n) fresh result", lambda: A + B),
    ("exp (n,n)", lambda: fnp.exp(A)),
    ("sqrt (n,)", lambda: fnp.sqrt(fnp.abs(V))),
    ("norm.cdf (n,) +astype", lambda: flops.stats.norm.cdf(V).astype(f32)),
    ("norm.pdf (n,) +astype", lambda: flops.stats.norm.pdf(V).astype(f32)),
    ("norm.cdf (n,n)", lambda: flops.stats.norm.cdf(A)),
    ("relu_moments (n,)", lambda: kit.relu_moments(V, fnp.abs(V) + 1.0)),
    ("sum axis 0 (n,n)", lambda: fnp.sum(A, axis=0)),
    ("einsum 'ij,ij->i'", lambda: fnp.einsum("ij,ij->i", A, B)),
    ("take rows (n idx)", lambda: kit.gather_rows(A, IDX)),
    ("reshape (n,n)->(2,512,2,512)", lambda: fnp.reshape(A, (2, 512, 2, 512))),
    ("transpose view", lambda: fnp.transpose(A)),
    ("diagonal view", lambda: fnp.diagonal(A)),
    ("copyto (n,n)", lambda: fnp.copyto(OUT, A)),
    ("cholesky n f32", lambda: fnp.linalg.cholesky(SPD)),
    ("cholesky n f64", lambda: fnp.linalg.cholesky(SPD64)),
    ("eigh n f32", lambda: fnp.linalg.eigh(SPD)),
    ("eigvalsh n f32", lambda: fnp.linalg.eigvalsh(SPD)),
    ("qr n f32 (reduced)", lambda: fnp.linalg.qr(A)),
    ("qr (n,256) f32", lambda: fnp.linalg.qr(R256c)),
    ("solve n, n rhs f32", lambda: fnp.linalg.solve(SPD, A)),
    ("solve n, 1 rhs f32", lambda: fnp.linalg.solve(SPD, Vn)),
    ("slogdet n f32", lambda: fnp.linalg.slogdet(SPD)),
    ("slogdet n f64", lambda: fnp.linalg.slogdet(SPD64)),
    ("det n f32", lambda: fnp.linalg.det(SPD)),
    ("inv n f32", lambda: fnp.linalg.inv(SPD)),
    ("svd n f32", lambda: fnp.linalg.svd(A)),
    ("slogdet batched (4096,8,8)", lambda: fnp.linalg.slogdet(SPD8)),
    ("cholesky batched (4096,8,8)", lambda: fnp.linalg.cholesky(SPD8)),
    ("solve batched (4096,8,8)@(4096,8,8)", lambda: fnp.linalg.solve(SPD8, T8)),
    ("eigh batched (1024,32,32)", lambda: fnp.linalg.eigh(SPD32)),
    ("slogdet batched (1024,32,32)", lambda: fnp.linalg.slogdet(SPD32)),
    ("cholesky batched (64,128,128)", lambda: fnp.linalg.cholesky(SPD128)),
    ("eigh batched (64,128,128)", lambda: fnp.linalg.eigh(SPD128)),
]

flt = sys.argv[2] if len(sys.argv) > 2 else ""
out = {}
for name, fn in CASES:
    if flt and flt not in name:
        continue
    rec = None
    try:
        for rep in range(3):
            with flops.BudgetContext(flop_budget=10 ** 16, quiet=True) as ctx:
                fn()
            rec = dict(units=ctx.flops_used / U, calls=len(ctx.op_log), resid_ms=1e3 * ctx.residual_wall_time_s,
                       backend_ms=1e3 * ctx.flopscope_backend_time_s, overhead_ms=1e3 * ctx.flopscope_overhead_time_s,
                       ops=sorted({r.op_name for r in ctx.op_log}))
    except Exception as e:  # noqa: BLE001
        rec = dict(error=f"{type(e).__name__}: {str(e)[:200]}")
    out[name] = rec
    if "error" in rec:
        print(f"{name:42s} ERROR {rec['error']}", flush=True)
    else:
        print(f"{name:42s} {rec['units']:9.4f} u {rec['calls']:4d} calls  resid {rec['resid_ms']:7.3f} ms  "
              f"backend {rec['backend_ms']:8.2f} ms", flush=True)
json.dump(out, open(sys.argv[1], "w"), indent=1)
