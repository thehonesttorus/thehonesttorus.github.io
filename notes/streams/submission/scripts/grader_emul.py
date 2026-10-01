"""Grader-transport emulation: the estimator runs in a client process against flopscope-client
0.12.1 (every op is a ZMQ round-trip), while a flopscope-server 0.12.1 process holds the arrays
and does the numpy work -- the architecture of the AIcrowd grader (2 vCPU solution process,
14 vCPU backend). Residual is the client's own `wall - dispatch` (flopscope-client
_budget._decompose_timing), i.e. exactly what the 0.4 s cap is applied to on the grader.

Orchestrator (run with the full-flopscope venv, which also has flopscope-server installed):
    /root/whest/bin/python grader_emul.py run --estimator EST.py --dataset DIR --out OUT.json \
        [--n-mlps N] [--server-threads 3] [--client-python /root/fcli/bin/python]
It starts the server, runs the client driver (this file, `client` mode, in the client venv),
samples peak RSS of both processes from outside, and writes one JSON.
"""
import argparse, glob, json, os, subprocess, sys, threading, time

URL_DEFAULT = "ipc:///tmp/fs-emul.sock"


def client_main(a):
    t_proc = float(a.t0)
    import importlib.util
    import numpy as np
    import pyarrow.parquet as pq
    import flopscope as flops
    import flopscope.numpy as fnp
    from whestbench import MLP, SetupContext

    meta = json.load(open(os.path.join(a.dataset, "metadata.json")))
    nsamp = float(meta["n_samples"])
    files = sorted(glob.glob(os.path.join(a.dataset, "data", "*.parquet")))
    tab = pq.read_table(files[0])
    nrows = tab.num_rows if not a.n_mlps else min(a.n_mlps, tab.num_rows)

    spec = importlib.util.spec_from_file_location("estimator", a.estimator)
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, os.path.dirname(os.path.abspath(a.estimator)))
    spec.loader.exec_module(mod)
    est = mod.Estimator()
    width, depth = int(meta["width"]), int(meta["depth"])
    t = time.time()
    with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):   # the grader's setup session
        est.setup(SetupContext(width=width, depth=depth, flop_budget=2 ** 41, api_version="1.0",
                               submission_dir=os.path.dirname(os.path.abspath(a.estimator)), seed=0))
    setup_only = time.time() - t
    setup_window = time.time() - t_proc   # interpreter start + imports + load + setup
    rows = []
    for i in range(nrows):
        r = tab.slice(i, 1).to_pylist()[0]
        w = np.asarray(r["weights"], dtype=np.float32)
        truth = np.asarray(r["final_means"], dtype=np.float32)
        with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):   # harness-side upload session
            weights = [fnp.asarray(w[l]) for l in range(w.shape[0])]
        mlp = MLP(width=w.shape[1], depth=w.shape[0], weights=weights, seed=int(r["mlp_seed"]) % (2 ** 63))
        rec = {"mlp_index": i, "mlp_seed": int(r["mlp_seed"]), "noise": float(r["avg_variance"]) / nsamp}
        ctx = flops.BudgetContext(flop_budget=2 ** 41, quiet=True)
        try:
            with ctx:
                if a.profile_out and i == nrows - 1:   # profile the last (steady-state) MLP only
                    import cProfile
                    prof = cProfile.Profile()
                    pred = prof.runcall(est.predict, mlp, 2 ** 41)
                    prof.dump_stats(a.profile_out)
                else:
                    pred = est.predict(mlp, 2 ** 41)
                arr = fnp.asarray(pred, dtype=fnp.float32)
                del pred
                fl = ctx.flops_used
            with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):   # harness-side fetch session
                p = np.asarray(arr.tolist(), dtype=np.float32)
            rec.update(ok=True, flops_used=int(fl), finite=bool(np.isfinite(p).all()),
                       final_layer_mse=float(np.mean((p[-1] - truth) ** 2)))
        except Exception as e:  # noqa
            rec.update(ok=False, error=f"{type(e).__name__}: {e}"[:500])
        rec.update(wall_time_s=ctx.wall_time_s, backend_s=ctx.flopscope_backend_time_s,
                   overhead_s=ctx.flopscope_overhead_time_s, residual_wall_time_s=ctx.residual_wall_time_s)
        rows.append(rec)
        print(json.dumps(rec), file=sys.stderr, flush=True)
        del weights, mlp
    print(json.dumps({"setup_only_s": setup_only, "setup_window_s": setup_window, "per_mlp": rows}))


def run_main(a):
    import psutil
    url = a.url
    sock = url.replace("ipc://", "")
    if os.path.exists(sock):
        os.remove(sock)
    senv = dict(os.environ)
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        senv[k] = str(a.server_threads)
    # the grader caps a single array at 4 GiB (starter kit); the server's own default is 100 MB
    senv.setdefault("FLOPSCOPE_MAX_ARRAY_BYTES", str(4 * 1024 ** 3))
    srv = subprocess.Popen([sys.executable, "-m", "flopscope_server", "--url", url, "--timeout", "900"],
                           env=senv, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    for _ in range(100):
        if os.path.exists(sock):
            break
        time.sleep(0.1)
    cenv = dict(os.environ)
    cenv["FLOPSCOPE_SERVER_URL"] = url
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        cenv[k] = "1"
    t0 = time.time()
    cli = subprocess.Popen([a.client_python, os.path.abspath(__file__), "client", "--estimator", os.path.abspath(a.estimator),
                            "--dataset", os.path.abspath(a.dataset), "--t0", repr(t0)] + (["--profile-out", a.profile_out] if a.profile_out else []) + (["--n-mlps", str(a.n_mlps)] if a.n_mlps else []),
                           env=cenv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    peak = {"server": 0, "client": 0}
    stop = threading.Event()

    def poll():
        ps, pc = psutil.Process(srv.pid), psutil.Process(cli.pid)
        while not stop.is_set():
            for key, p in (("server", ps), ("client", pc)):
                try:
                    peak[key] = max(peak[key], p.memory_info().rss)
                except psutil.Error:
                    pass
            time.sleep(0.05)

    th = threading.Thread(target=poll, daemon=True)
    th.start()
    out, err = cli.communicate()
    stop.set()
    th.join()
    srv.terminate()
    try:
        srv.wait(10)
    except Exception:
        srv.kill()
    try:
        res = json.loads(out.strip().splitlines()[-1])
    except Exception:
        res = {"per_mlp": [], "client_stdout_tail": out[-2000:]}
    per = res.get("per_mlp", [])
    ok = [r for r in per if r.get("ok")]
    summ = {"n_mlps": len(per), "n_errors": len(per) - len(ok)}
    if ok:
        summ.update(
            C_over_B_max=max(r["flops_used"] for r in ok) / 2 ** 41,
            C_over_B_mean=sum(r["flops_used"] for r in ok) / len(ok) / 2 ** 41,
            final_layer_mse=sum(r["final_layer_mse"] for r in ok) / len(ok),
            truth_noise_floor=sum(r["noise"] for r in ok) / len(ok),
            residual_max_s=max(r["residual_wall_time_s"] for r in ok),
            residual_list_s=[round(r["residual_wall_time_s"], 4) for r in per],
            wall_max_s=max(r["wall_time_s"] for r in ok),
            n_over_residual_cap=sum(1 for r in ok if r["residual_wall_time_s"] > 0.4),
            n_over_wall_cap=sum(1 for r in ok if r["wall_time_s"] > 120.0),
            n_nonfinite=sum(1 for r in ok if not r["finite"]),
        )
        summ["final_mse_minus_noise"] = summ["final_layer_mse"] - summ["truth_noise_floor"]
    rec = {"tag": a.tag, "mode": "client-server emulation", "server_threads": a.server_threads, "client_threads": 1,
           "setup_window_s": res.get("setup_window_s"), "setup_only_s": res.get("setup_only_s"),
           "peak_rss_mb_server": round(peak["server"] / 2 ** 20), "peak_rss_mb_client": round(peak["client"] / 2 ** 20),
           "summary": summ, "per_mlp": per, "client_stderr_tail": err[-3000:], "server_stderr_tail": (srv.stderr.read() or "")[-1500:]}
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    json.dump(rec, open(a.out, "w"), indent=1)
    print(json.dumps({"tag": a.tag, "setup_window_s": rec["setup_window_s"], "rss_server": rec["peak_rss_mb_server"],
                      "rss_client": rec["peak_rss_mb_client"], **{k: v for k, v in summ.items()}}))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["run", "client"])
    ap.add_argument("--estimator", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--n-mlps", type=int, default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--tag", default="emul")
    ap.add_argument("--t0", default="0")
    ap.add_argument("--url", default=URL_DEFAULT)
    ap.add_argument("--server-threads", type=int, default=3)
    ap.add_argument("--client-python", default="/root/fcli/bin/python")
    ap.add_argument("--profile-out", default=None, help="client mode: cProfile the last MLP's predict() into this file")
    a = ap.parse_args()
    if a.mode == "client":
        client_main(a)
    else:
        run_main(a)


if __name__ == "__main__":
    main()
