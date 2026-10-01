"""Run `whest run` exactly as the grader would (subprocess or local runner, graded caps) and
record what the report does not: peak RSS of the worker process tree and total wall time.

    python measure_run.py --estimator bundles/v29/estimator.py --dataset DIR --runner subprocess \
        --max-threads 2 --tag v29_dev6_sub_r1 --out results/v29_dev6_sub_r1.json [--n-mlps N] [--extra "..."]

Writes {"tag", "cmd", "harness_wall_s", "peak_rss_mb", "rc", "summary", "per_mlp", "report"}.
The summary subtracts the truth noise floor avg_variance/N (read from the baked parquet) from the
final-layer MSE when the dataset directory is given.
"""
import argparse, glob, json, os, subprocess, sys, threading, time

import psutil


def tree_rss(p):
    tot = 0
    try:
        procs = [p] + p.children(recursive=True)
    except psutil.Error:
        return 0
    for q in procs:
        try:
            tot += q.memory_info().rss
        except psutil.Error:
            pass
    return tot


def worker_rss(p):
    """RSS of the largest single descendant (the subprocess worker) or of p itself (local)."""
    best = 0
    try:
        procs = [p] + p.children(recursive=True)
    except psutil.Error:
        return 0
    for q in procs:
        try:
            best = max(best, q.memory_info().rss)
        except psutil.Error:
            pass
    return best


def truth_noise(dataset):
    try:
        import pyarrow.parquet as pq
        meta = json.load(open(os.path.join(dataset, "metadata.json")))
        n = float(meta.get("n_samples"))
        files = sorted(glob.glob(os.path.join(dataset, "data", "*.parquet")))
        t = pq.read_table(files[0], columns=["mlp_id", "avg_variance"])
        return {int(i): float(v) / n for i, v in zip(t.column("mlp_id").to_pylist(), t.column("avg_variance").to_pylist())}, n
    except Exception as e:  # noqa
        print("noise floor unavailable:", e, file=sys.stderr)
        return {}, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--estimator", required=True)
    ap.add_argument("--dataset", default=None)
    ap.add_argument("--runner", default="subprocess")
    ap.add_argument("--max-threads", type=int, default=2)
    ap.add_argument("--n-mlps", type=int, default=None)
    ap.add_argument("--wall-time-limit", type=float, default=120.0)
    ap.add_argument("--extra", default="")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cmd = ["whest", "run", "--estimator", os.path.abspath(a.estimator), "--runner", a.runner,
           "--format", "json", "--profile", "--max-threads", str(a.max_threads),
           "--wall-time-limit", str(a.wall_time_limit)]
    if a.dataset:
        cmd += ["--dataset", os.path.abspath(a.dataset)]
    if a.n_mlps:
        cmd += ["--n-mlps", str(a.n_mlps)]
    if a.extra:
        cmd += a.extra.split()
    env = dict(os.environ)
    env.setdefault("WHEST_SKIP_HARDWARE_FALLBACK_PROBES", "1")
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env[k] = str(a.max_threads)
    t0 = time.time()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env,
                            cwd=os.path.dirname(os.path.abspath(a.estimator)))
    ps = psutil.Process(proc.pid)
    peak = {"tree": 0, "worker": 0}
    stop = threading.Event()

    def poll():
        while not stop.is_set():
            peak["tree"] = max(peak["tree"], tree_rss(ps))
            peak["worker"] = max(peak["worker"], worker_rss(ps))
            time.sleep(0.05)

    th = threading.Thread(target=poll, daemon=True)
    th.start()
    out, err = proc.communicate()
    stop.set()
    th.join()
    wall = time.time() - t0
    report = None
    try:
        report = json.loads(out)
    except Exception:
        i = out.find("{")
        if i >= 0:
            try:
                report = json.loads(out[i:])
            except Exception:
                report = None
    res = (report or {}).get("results", report or {})
    per = res.get("per_mlp", []) if isinstance(res, dict) else []
    noise, nsamp = truth_noise(a.dataset) if a.dataset else ({}, None)
    rows = []
    for i, m in enumerate(per):
        r = {k: m.get(k) for k in ("mlp_id", "mlp_index", "final_layer_mse", "all_layers_mse", "flops_used",
                                   "compute_utilization", "wall_time_s", "residual_wall_time_s",
                                   "flopscope_backend_time_s", "flopscope_overhead_time_s",
                                   "budget_exhausted", "time_exhausted", "residual_wall_time_exhausted",
                                   "error_code", "adjusted_final_layer_score", "score_multiplier") if k in m}
        rows.append(r)
    nz = list(noise.values())
    summ = {}
    if isinstance(res, dict):
        summ = {k: res.get(k) for k in ("adjusted_final_layer_score", "final_layer_mse", "all_layers_mse",
                                         "mean_effective_compute", "mean_compute_utilization", "n_failed_mlps",
                                         "mean_score_multiplier") if k in res}
    if rows:
        fl = [r.get("flops_used") for r in rows if r.get("flops_used") is not None]
        rs = [r.get("residual_wall_time_s") for r in rows if r.get("residual_wall_time_s") is not None]
        wt = [r.get("wall_time_s") for r in rows if r.get("wall_time_s") is not None]
        if fl:
            summ["C_over_B_max"] = max(fl) / 2 ** 41
            summ["C_over_B_mean"] = sum(fl) / len(fl) / 2 ** 41
        if rs:
            summ["residual_max_s"] = max(rs)
            summ["residual_mean_s"] = sum(rs) / len(rs)
            summ["residual_list_s"] = rs
        if wt:
            summ["wall_max_s"] = max(wt)
            summ["wall_mean_s"] = sum(wt) / len(wt)
        failed = [r for r in rows if r.get("error_code") or r.get("budget_exhausted") or r.get("time_exhausted")
                  or r.get("residual_wall_time_exhausted")]
        summ["n_failed_rows"] = len(failed)
    if nz and summ.get("final_layer_mse") is not None:
        k = len(rows) if rows else len(nz)
        nf = sum(nz[:k]) / k
        summ["truth_noise_floor"] = nf
        summ["final_mse_minus_noise"] = summ["final_layer_mse"] - nf
        summ["truth_n_samples"] = nsamp
    rec = {"tag": a.tag, "cmd": cmd, "harness_wall_s": round(wall, 2), "peak_rss_mb_tree": round(peak["tree"] / 2 ** 20),
           "peak_rss_mb_worker": round(peak["worker"] / 2 ** 20), "rc": proc.returncode, "nproc": os.cpu_count(),
           "loadavg_at_end": os.getloadavg(), "summary": summ, "per_mlp": rows, "report": report,
           "stderr_tail": err[-3000:]}
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out, "w") as f:
        json.dump(rec, f, indent=1, default=str)
    print(json.dumps({"tag": a.tag, "rc": proc.returncode, "wall": rec["harness_wall_s"],
                      "peak_rss_mb_worker": rec["peak_rss_mb_worker"], **{k: v for k, v in summ.items() if k != "residual_list_s"}},
                     default=str))
    return 0 if report is not None else 2


if __name__ == "__main__":
    sys.exit(main())
