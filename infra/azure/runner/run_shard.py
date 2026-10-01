"""One grid task: evaluate one estimator bundle on one dataset shard with the official
harness and write the JSON report.  Works identically on a laptop (local paths) and inside
an Azure Batch container (blob URLs with SAS).

    python run_shard.py --bundle <dir | .tar.gz | https://...blob...tar.gz?sas> \
        --shard <baked dataset dir | parquet file | https://...parquet?sas> --split full \
        --out <file.json | https://...blob/.../name.json?sas> [--tag NAME] \
        [--wall-time-limit 600] [--max-threads 4] [--n-mlps N] [--width 1024 --depth 16]

A parquet shard (one file of the public dataset, ~16 MLPs) is wrapped into a baked-dataset
directory (metadata.json + data/<split>-00000-of-00001.parquet), which is what `whest run
--dataset <dir>` accepts.  The graded caps (2**41 FLOPs, 0.4 s residual, 5 s setup) are kept;
only the wall-clock cap is relaxed, because FLOP accounting does not depend on the machine.
"""
import argparse, datetime, json, os, platform, shutil, subprocess, sys, tarfile, tempfile, time, urllib.parse


def is_url(s):
    return s.startswith("http://") or s.startswith("https://")


def fetch(src, dest):
    """Copy a local path or download a blob URL (with SAS) to dest."""
    if is_url(src):
        from azure.storage.blob import BlobClient
        bc = BlobClient.from_blob_url(src)
        with open(dest, "wb") as f:
            bc.download_blob(max_concurrency=8).readinto(f)
    else:
        shutil.copyfile(src, dest)
    return dest


def put(src_path, dest):
    if is_url(dest):
        from azure.storage.blob import BlobClient
        bc = BlobClient.from_blob_url(dest)
        with open(src_path, "rb") as f:
            bc.upload_blob(f, overwrite=True)
    else:
        os.makedirs(os.path.dirname(os.path.abspath(dest)) or ".", exist_ok=True)
        shutil.copyfile(src_path, dest)


def prepare_bundle(bundle, work):
    bdir = os.path.join(work, "bundle")
    os.makedirs(bdir, exist_ok=True)
    if os.path.isdir(bundle):
        shutil.copytree(bundle, bdir, dirs_exist_ok=True)
    else:
        local = fetch(bundle, os.path.join(work, "bundle.tar.gz"))
        with tarfile.open(local) as tf:
            tf.extractall(bdir)
        # tolerate a single top-level folder inside the tarball
        entries = os.listdir(bdir)
        if "estimator.py" not in entries and len(entries) == 1 and os.path.isdir(os.path.join(bdir, entries[0])):
            bdir = os.path.join(bdir, entries[0])
    est = os.path.join(bdir, "estimator.py")
    if not os.path.exists(est):
        sys.exit(f"no estimator.py in bundle {bundle}")
    return est


def prepare_dataset(shard, split, width, depth, work):
    if os.path.isdir(shard) and os.path.exists(os.path.join(shard, "metadata.json")):
        return shard
    ddir = os.path.join(work, "dataset")
    os.makedirs(os.path.join(ddir, "data"), exist_ok=True)
    local = fetch(shard, os.path.join(ddir, "data", f"{split}-00000-of-00001.parquet"))
    import pyarrow.parquet as pq
    n = pq.ParquetFile(local).metadata.num_rows
    meta = {
        "schema_version": "3.0", "format": "hf-datasets-parquet", "backend": "shard",
        "seed_protocol": {"name": "whestbench_explicit_per_mlp_seeds", "version": "3.0"},
        "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "split": split, "config": "default", "n_mlps": int(n), "n_samples": 1_000_000_000,
        "width": width, "depth": depth, "source_shard": os.path.basename(urllib.parse.urlparse(shard).path),
    }
    with open(os.path.join(ddir, "metadata.json"), "w") as f:
        json.dump(meta, f, indent=1)
    return ddir


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--shard", required=True)
    ap.add_argument("--split", default="full")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", default="")
    ap.add_argument("--wall-time-limit", type=float, default=600.0)
    ap.add_argument("--max-threads", type=int, default=int(os.environ.get("OMP_NUM_THREADS", "4")))
    ap.add_argument("--n-mlps", type=int, default=None)
    ap.add_argument("--width", type=int, default=1024)
    ap.add_argument("--depth", type=int, default=16)
    ap.add_argument("--runner", default="local", help="local: fast, in-process; subprocess: the grader's transport (adds ~30 s per run)")
    ap.add_argument("--extra", default="", help="extra args passed verbatim to `whest run`")
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    work = tempfile.mkdtemp(prefix="shard-")
    try:
        est = prepare_bundle(args.bundle, work)
        ddir = prepare_dataset(args.shard, args.split, args.width, args.depth, work)
        cmd = ["whest", "run", "--estimator", est, "--dataset", ddir, "--split", args.split,
               "--runner", args.runner, "--format", "json", "--wall-time-limit", str(args.wall_time_limit),
               "--max-threads", str(args.max_threads), "--profile"]
        if args.n_mlps:
            cmd += ["--n-mlps", str(args.n_mlps)]
        if args.extra:
            cmd += args.extra.split()
        print("running:", " ".join(cmd), flush=True)
        t1 = time.time()
        proc = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(est))
        report = None
        try:
            report = json.loads(proc.stdout)
        except Exception:
            # the CLI may print non-JSON before the report; find the last JSON object
            idx = proc.stdout.rfind("\n{")
            if idx >= 0:
                try:
                    report = json.loads(proc.stdout[idx + 1:])
                except Exception:
                    report = None
        out = {
            "task": {"tag": args.tag, "bundle": args.bundle.split("?")[0], "shard": args.shard.split("?")[0],
                      "split": args.split, "host": platform.node(), "cpu_count": os.cpu_count(),
                      "max_threads": args.max_threads, "wall_time_limit": args.wall_time_limit,
                      "started_utc": datetime.datetime.fromtimestamp(t0, datetime.timezone.utc).isoformat(),
                      "harness_seconds": round(time.time() - t1, 1), "returncode": proc.returncode},
            "report": report,
            "stderr_tail": proc.stderr[-4000:],
            "stdout_tail": None if report is not None else proc.stdout[-4000:],
        }
        local_out = os.path.join(work, "result.json")
        with open(local_out, "w") as f:
            json.dump(out, f)
        put(local_out, args.out)
        if report is None:
            print("harness produced no JSON report; rc", proc.returncode, file=sys.stderr)
            print(proc.stderr[-2000:], file=sys.stderr)
            return 2
        res = report.get("results", report)
        keys = ("adjusted_final_layer_score", "final_layer_mse", "all_layers_mse", "mean_effective_compute",
                "mean_compute_utilization", "n_failed_mlps")
        print("result:", {k: res.get(k) for k in keys}, "n_mlps:", len(res.get("per_mlp", [])), flush=True)
        return 0
    finally:
        if not args.keep:
            shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
