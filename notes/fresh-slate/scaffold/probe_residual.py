"""Residual seconds per flopscope call, by call kind, in-process (full flopscope) or client/server
(flopscope-client against a flopscope-server process, the grader's transport).

  in-process:     OPENBLAS_NUM_THREADS=2 /root/whest/bin/python probe_residual.py OUT.json
  client/server:  /root/whest/bin/python probe_residual.py OUT.json --emulate
                  (starts `python -m flopscope_server` from /root/whest with 2 threads and re-runs this
                   script as the client under /root/fcli/bin/python)
Each kind runs N calls in a loop inside one BudgetContext; residual / N is reported (best of 3 runs).
"""
import json, os, subprocess, sys, time

N = 2000


def kinds(fnp, flops):
    f32 = fnp.float32

    def mk(shape):
        b = fnp.empty(shape, dtype=f32)
        fnp.copyto(b, 1.0)
        return b
    v = mk((1024,)); v2 = mk((1024,)); vo = mk((1024,))
    m = mk((1024, 1024)); m2 = mk((1024, 1024)); mo = mk((1024, 1024))
    s = mk((64, 64)); so = mk((64, 64))
    st = mk((8, 7, 64, 64)); sto = mk((8, 7, 64, 64))
    a0, a1, o2 = st[:, 0], st[:, 1], sto[:, 2]
    ks = {
        "vector add out= (n,)": lambda: fnp.add(v, v2, out=vo),
        "vector add fresh (n,)": lambda: v + v2,
        "small matmul out= (64,64)": lambda: fnp.matmul(s, s, out=so),
        "(n,n) add out=": lambda: fnp.add(m, m2, out=mo),
        "(n,n) add fresh": lambda: m + m2,
        "strided-slab add out=, prebuilt views": lambda: fnp.add(a0, a1, out=o2),
        "strided-slab add out=, slicing in the call": lambda: fnp.add(st[:, 0], st[:, 1], out=sto[:, 2]),
        "basic slice only (no op)": lambda: st[:, 3],
        "reshape (logged view)": lambda: fnp.reshape(s, (2, 32, 2, 32)),
        "vector norm.cdf + astype": lambda: flops.stats.norm.cdf(v).astype(f32),
    }
    return ks


def run(emulated):
    import flopscope as flops
    import flopscope.numpy as fnp
    out = {}
    with flops.BudgetContext(flop_budget=10 ** 16, quiet=True):
        ks = kinds(fnp, flops)
    for name, fn in ks.items():
        best = None
        for rep in range(3):
            with flops.BudgetContext(flop_budget=10 ** 16, quiet=True) as ctx:
                t = time.perf_counter()
                for _ in range(N):
                    fn()
                w = time.perf_counter() - t
            r = dict(resid_ms_per_call=1e3 * ctx.residual_wall_time_s / N, wall_ms_per_call=1e3 * w / N,
                     logged_ops_per_call=len(ctx.op_log) / N if hasattr(ctx, "op_log") else None)
            if best is None or r["resid_ms_per_call"] < best["resid_ms_per_call"]:
                best = r
        out[name] = best
        print(f"{name:45s} resid {best['resid_ms_per_call']:.4f} ms/call  wall {best['wall_ms_per_call']:.4f} ms/call",
              flush=True)
    return out


if __name__ == "__main__":
    dst = sys.argv[1]
    if "--client" in sys.argv:
        json.dump(run(True), open(dst, "w"), indent=1)
    elif "--emulate" in sys.argv:
        url = "ipc:///tmp/probe-residual.sock"
        senv = dict(os.environ, OPENBLAS_NUM_THREADS="2", OMP_NUM_THREADS="2", MKL_NUM_THREADS="2")
        srv = subprocess.Popen([sys.executable, "-m", "flopscope_server", "--url", url, "--timeout", "900"], env=senv,
                               stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        time.sleep(3)
        cenv = dict(os.environ, FLOPSCOPE_SERVER_URL=url, OPENBLAS_NUM_THREADS="1")
        try:
            subprocess.run(["/root/fcli/bin/python", os.path.abspath(__file__), dst, "--client"], env=cenv, check=True)
        finally:
            srv.terminate()
    else:
        json.dump(run(False), open(dst, "w"), indent=1)
