#!/usr/bin/env python3
"""Submit, watch and collect an experiment grid on the Azure Batch pools created by 04_pools.sh.

    ./grid.py submit  --name exp1 --bundle path/to/estimator_dir [--bundle other.tar.gz ...] \
                      --split full|mini [--shards 0-62 | --shards 7] [--wall-time-limit 600] [--n-mlps N]
    ./grid.py status  --name exp1 [--failures]
    ./grid.py collect --name exp1 [--csv exp1.csv]
    ./grid.py bake    --name fresh-A --n-mlps 1000 [--n-samples 100000000] [--slices 16]
    ./grid.py atlas   --name d8b --count 2048 [--start 0] [--pairs]

One task = (bundle, shard).  Tasks are spread round-robin over the regions in env.sh that have
a Batch account and a pool; each region gets a job named <name>-<region>, which terminates by
itself once all its tasks have completed (so it stops counting against the active-job quota;
use a new --name for a new submission).  Results land in the `results` container under
<name>/<bundle-tag>/<shard>.json and are aggregated by `collect` (mean of per-MLP adjusted
scores across all shards, i.e. exactly the leaderboard quantity over the shards run).

Uses the `az` CLI only (no SDK), so it runs wherever `az login` works.  Batch calls use
shared-key login (see env.sh); storage calls use the account key.
"""
import argparse, collections, datetime, json, os, re, shlex, subprocess, sys, tarfile, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ID_RE = re.compile(r"^[A-Za-z0-9_-]{1,40}$")      # --name: becomes part of job ids and blob paths
TASK_USER = {"autoUser": {"scope": "pool", "elevationLevel": "admin"}}  # root inside the container


def env():
    p = subprocess.run(["bash", "-c", f"source {shlex.quote(os.path.join(HERE, 'env.sh'))} >/dev/null && env -0"],
                       capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"env.sh failed: {p.stderr.strip()[-800:]}")
    return dict(item.split("=", 1) for item in p.stdout.split("\0") if "=" in item)


def az(*args, json_out=True):
    cmd = ["az", *args]
    if json_out:
        cmd += ["-o", "json"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"az {' '.join(args[:3])} failed: {p.stderr.strip()[:800]}")
    return json.loads(p.stdout) if json_out and p.stdout.strip() else None


def account_key(e):
    return az("storage", "account", "keys", "list", "-n", e["STORAGE"], "-g", e["RG"])[0]["value"]


def container_sas(e, key, container, perms, hours):
    expiry = (datetime.datetime.now(datetime.timezone.utc) + datetime.timedelta(hours=hours)).strftime("%Y-%m-%dT%H:%MZ")
    return az("storage", "container", "generate-sas", "-n", container, "--account-name", e["STORAGE"],
              "--account-key", key, "--permissions", perms, "--expiry", expiry)


def blob_url(e, container, name, sas):
    return f"https://{e['STORAGE']}.blob.core.windows.net/{container}/{name}?{sas}"


def container_url(e, container, sas):
    return f"https://{e['STORAGE']}.blob.core.windows.net/{container}?{sas}"


def batch_account_name(e, region):
    sub = e["SUB_ID"].replace("-", "")
    return f"{e['PREFIX']}b{region.replace('-', '')[:10]}{sub[:4]}"


def batch_login(e, region):
    """Point this process's az calls at one region's Batch account (shared key).  `az batch account
    login` alone keeps the account in ~/.azure/config, which every az process shares, so a concurrent
    run (`status` while `submit` runs, or 04_pools.sh) would switch the account under this one and
    its jobs would land, silently pending, in another region's account.  The AZURE_BATCH_* variables
    take precedence over that file and are private to this process and its az children."""
    acct = batch_account_name(e, region)
    c = az("batch", "account", "login", "-n", acct, "-g", e["RG"], "--shared-key-auth", "--show")
    os.environ.update(AZURE_BATCH_ACCOUNT=c["account"], AZURE_BATCH_ENDPOINT=c["endpoint"],
                      AZURE_BATCH_ACCESS_KEY=c["primaryKey"], AZURE_BATCH_AUTH_MODE="shared_key")
    return acct


def pool_exists(pool):
    try:
        az("batch", "pool", "show", "--pool-id", pool, "--select", "id")
        return True
    except RuntimeError:
        return False


def ensure_job(job, pool):
    try:
        az("batch", "job", "create", "--id", job, "--pool-id", pool, json_out=False)
    except RuntimeError as ex:
        if "JobExists" not in str(ex) and "already exists" not in str(ex):
            raise
        print(f"job {job} exists; adding to it", file=sys.stderr)


def add_tasks(job, tasks):
    """Bulk-add (100 per request); returns the number of tasks not added.  The CLI fails on client
    errors other than TaskExists and reports per-task status rows, which are checked here."""
    n_bad = 0
    for k in range(0, len(tasks), 100):
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
            json.dump(tasks[k:k + 100], f)
        try:
            res = az("batch", "task", "create", "--job-id", job, "--json-file", f.name) or []
        finally:
            os.unlink(f.name)
        for r in res if isinstance(res, list) else [res]:
            status = str(r.get("status", "success")).lower()
            code = (r.get("error") or {}).get("code")
            if status != "success":
                if code == "TaskExists":
                    print(f"  task {r.get('taskId')} already exists, left as is", file=sys.stderr)
                else:
                    n_bad += 1
                    print(f"  task {r.get('taskId')} NOT added: {status} {code}", file=sys.stderr)
    return n_bad


def finish_job(job, keep_open):
    if not keep_open:
        az("batch", "job", "set", "--job-id", job, "--on-all-tasks-complete", "terminatejob", json_out=False)


def list_shards(e, key, split):
    blobs = az("storage", "blob", "list", "-c", e["CONTAINER_DATASET"], "--account-name", e["STORAGE"],
               "--account-key", key, "--prefix", f"{e['HF_REVISION']}/data/{split}-", "--num-results", "*")
    return sorted(b["name"] for b in blobs if b["name"].endswith(".parquet"))


def bundle_tag(bundle):
    base = os.path.basename(os.path.normpath(bundle))
    for ext in (".tar.gz", ".tgz"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return re.sub(r"[^A-Za-z0-9._-]", "-", base)


def upload_bundle(e, key, bundle, name):
    tag = bundle_tag(bundle)
    tmp = None
    if os.path.isdir(bundle):
        tmp = tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=False).name
        with tarfile.open(tmp, "w:gz") as tf:
            for root, dirs, files in os.walk(bundle):
                dirs[:] = [d for d in dirs if d != "__pycache__"]
                for f in files:
                    full = os.path.join(root, f)
                    tf.add(full, arcname=os.path.relpath(full, bundle))
        path = tmp
    else:
        path = bundle
    blob = f"{name}/{tag}.tar.gz"
    try:
        az("storage", "blob", "upload", "-c", e["CONTAINER_SUBMISSIONS"], "--account-name", e["STORAGE"],
           "--account-key", key, "-f", path, "-n", blob, "--overwrite", json_out=False)
    finally:
        if tmp:
            os.unlink(tmp)
    return tag, blob


def task_id(*parts):
    return re.sub(r"[^A-Za-z0-9_-]", "-", "-".join(parts))[:64]


def usable_regions(e, regions, pool_of):
    """Regions whose Batch account accepts a login and that have the pool; others are skipped."""
    ok = []
    for r in regions:
        try:
            batch_login(e, r)
        except RuntimeError as ex:
            print(f"skip {r}: no Batch login ({str(ex)[:200]})", file=sys.stderr)
            continue
        if not pool_exists(pool_of(r)):
            print(f"skip {r}: pool {pool_of(r)} not found (04_pools.sh)", file=sys.stderr)
            continue
        ok.append(r)
    return ok


def parse_shards(spec, shards):
    if not spec:
        return shards
    if "-" in spec:
        lo, hi = (int(x) for x in spec.split("-", 1))
    else:
        lo = hi = int(spec)
    return shards[lo:hi + 1]


def check_name(name):
    if not ID_RE.match(name):
        sys.exit(f"--name {name!r}: use 1-40 letters, digits, '-' or '_' (it becomes part of Batch job ids)")


def build_submit_tasks(e, args, bundles, shards, sas_data, sas_sub, sas_res):
    image = f"{e['ACR']}.azurecr.io/{e['IMAGE']}"
    tasks = []
    for tag, bblob in bundles:
        for s in shards:
            sname = os.path.basename(s).replace(".parquet", "")
            out = blob_url(e, e["CONTAINER_RESULTS"], f"{args.name}/{tag}/{sname}.json", sas_res)
            cmd = ["python", "/app/run_shard.py",
                   "--bundle", blob_url(e, e["CONTAINER_SUBMISSIONS"], bblob, sas_sub),
                   "--shard", blob_url(e, e["CONTAINER_DATASET"], s, sas_data), "--split", args.split,
                   "--out", out, "--tag", tag, "--wall-time-limit", str(args.wall_time_limit),
                   "--max-threads", e["VCPU_PER_TASK"], "--runner", args.runner]
            if args.n_mlps:
                cmd += ["--n-mlps", str(args.n_mlps)]
            if args.extra:
                cmd += [f"--extra={args.extra}"]
            tasks.append({"id": task_id(tag[:40], sname),
                          "commandLine": "/bin/bash -c " + shlex.quote(shlex.join(cmd)),
                          "containerSettings": {"imageName": image, "containerRunOptions": "--rm"},
                          "constraints": {"maxWallClockTime": f"PT{args.task_hours}H", "maxTaskRetryCount": 1},
                          "userIdentity": TASK_USER})
    ids = [t["id"] for t in tasks]
    dup = [i for i, c in collections.Counter(ids).items() if c > 1]
    if dup:
        sys.exit(f"task id collision (bundle tags too similar or duplicated): {dup[:5]}")
    return tasks


def cmd_submit(args):
    check_name(args.name)
    tags = [bundle_tag(b) for b in args.bundle]
    if len(set(tags)) != len(tags):
        sys.exit(f"bundle tags must be distinct (they name the result folders): {tags}")
    e = env()
    key = account_key(e)
    shards = parse_shards(args.shards, list_shards(e, key, args.split))
    if not shards:
        sys.exit("no shards found; run 03_stage_dataset.sh first")
    regions = usable_regions(e, (args.regions or e["REGIONS"]).split(), lambda r: f"{e['PREFIX']}-pool-{r}")
    if not regions:
        sys.exit("no region has a usable Batch account + pool")
    sas_data = container_sas(e, key, e["CONTAINER_DATASET"], "rl", args.sas_hours)
    sas_sub = container_sas(e, key, e["CONTAINER_SUBMISSIONS"], "rl", args.sas_hours)
    sas_res = container_sas(e, key, e["CONTAINER_RESULTS"], "rcwl", args.sas_hours)
    bundles = [upload_bundle(e, key, b, args.name) for b in args.bundle]
    tasks = build_submit_tasks(e, args, bundles, shards, sas_data, sas_sub, sas_res)
    manifest = {"name": args.name, "split": args.split, "bundles": [t for t, _ in bundles], "shards": shards,
                "regions": regions, "n_tasks": len(tasks), "jobs": {},
                "submitted_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    per_region = collections.defaultdict(list)
    for i, t in enumerate(tasks):
        per_region[regions[i % len(regions)]].append(t)
    mpath = f"{args.name}.manifest.json"
    n_bad = 0
    try:
        for r, ts in per_region.items():
            batch_login(e, r)
            job = f"{args.name}-{r}"
            ensure_job(job, f"{e['PREFIX']}-pool-{r}")
            manifest["jobs"][r] = {"job": job, "n_tasks": len(ts)}
            bad = add_tasks(job, ts)
            n_bad += bad
            finish_job(job, args.keep_job_open)
            print(f"{r}: {len(ts) - bad} tasks in job {job}")
    finally:
        with open(mpath, "w") as f:      # written even if a region fails half-way, so status/collect work
            json.dump(manifest, f, indent=1)
    print(f"submitted {len(tasks) - n_bad} tasks ({len(bundles)} bundles x {len(shards)} shards) across "
          f"{len(per_region)} regions; manifest {mpath}")
    if n_bad:
        sys.exit(f"{n_bad} tasks were not added")


def cmd_status(args):
    e = env()
    with open(f"{args.name}.manifest.json") as f:
        man = json.load(f)
    total = collections.Counter()
    for r, j in man["jobs"].items():
        batch_login(e, r)
        counts = az("batch", "job", "task-counts", "show", "--job-id", j["job"])
        tc = counts.get("taskCounts", counts)
        row = {k: tc.get(k, 0) for k in ("active", "running", "completed", "succeeded", "failed")}
        total.update(row)
        print(f"{r:<16} {row}")
        if args.failures and row["failed"]:
            failed = az("batch", "task", "list", "--job-id", j["job"], "--filter", "executionInfo/result eq 'failure'",
                        "--query", "[].{id:id, exit:executionInfo.exitCode, error:executionInfo.failureInfo.code}")
            for t in failed or []:
                print(f"    failed {t}")
    print("total", dict(total))


def summarize(tag, files):
    import statistics as st
    per_mlp, shards_done = [], 0
    for fn in files:
        with open(fn) as f:
            res = json.load(f)
        rep = res.get("report")
        if not rep:
            continue
        shards_done += 1
        per_mlp.extend(rep.get("results", rep).get("per_mlp", []))
    if not per_mlp:
        return (tag, shards_done, 0, None, None, None, None, None, None)
    adj = [m["adjusted_final_layer_score"] for m in per_mlp if m.get("adjusted_final_layer_score") is not None]
    mse = [m["final_layer_mse"] for m in per_mlp if m.get("final_layer_mse") is not None]
    allm = [m["all_layers_mse"] for m in per_mlp if m.get("all_layers_mse") is not None]
    fl = [m.get("flops_used") or m.get("effective_compute") or 0 for m in per_mlp]
    res_t = [m.get("residual_wall_time_s") or 0 for m in per_mlp]
    failed = sum(1 for m in per_mlp if m.get("budget_exhausted") or m.get("time_exhausted") or m.get("residual_wall_time_exhausted")
                 or m.get("combined_budget_exhausted") or m.get("error") or m.get("traceback"))
    return (tag, shards_done, len(per_mlp), st.mean(adj) if adj else None, st.mean(mse) if mse else None,
            st.mean(allm) if allm else None, st.mean(fl) / 2**41 if fl else None, max(res_t) if res_t else None, failed)


def print_table(rows, csv=None):
    def f(x, w):
        return f"{x:{w}.3e}" if isinstance(x, float) else f"{str(x):>{w}}"
    print(f"{'bundle':<32} {'shards':>6} {'mlps':>5} {'adjusted':>11} {'final_mse':>11} {'all_mse':>11} {'util':>9} {'max_resid':>9} {'failed':>6}")
    for r in rows:
        print(f"{r[0]:<32} {r[1]:>6} {r[2]:>5} {f(r[3], 11)} {f(r[4], 11)} {f(r[5], 11)} {f(r[6], 9)} {f(r[7], 9)} {f(r[8], 6)}")
    if csv:
        with open(csv, "w") as fcsv:
            fcsv.write("bundle,shards,mlps,adjusted,final_mse,all_mse,util,max_residual_s,failed\n")
            for r in rows:
                fcsv.write(",".join("" if x is None else str(x) for x in r) + "\n")


def cmd_collect(args):
    e = env()
    key = account_key(e)
    with open(f"{args.name}.manifest.json") as f:
        man = json.load(f)
    tmp = tempfile.mkdtemp()
    az("storage", "blob", "download-batch", "-s", e["CONTAINER_RESULTS"], "-d", tmp, "--pattern", f"{args.name}/*",
       "--account-name", e["STORAGE"], "--account-key", key, "--no-progress", json_out=False)
    rows = []
    for tag in man["bundles"]:
        d = os.path.join(tmp, args.name, tag)
        files = [os.path.join(d, fn) for fn in sorted(os.listdir(d)) if fn.endswith(".json")] if os.path.isdir(d) else []
        rows.append(summarize(tag, files))
    print_table(rows, args.csv)


def gpu_job(e, args, prefix):
    region = args.region or e["HOME_REGION"]
    batch_login(e, region)
    pool = f"{e['PREFIX']}-gpu-{region}"
    if not pool_exists(pool):
        sys.exit(f"GPU pool {pool} not found: run 05_gpu_pool.sh (GPU_REGION={region})")
    job = f"{prefix}-{args.name}"
    ensure_job(job, pool)
    return job


def cmd_bake(args):
    """Shard a fresh-seed bake over the GPU pool: M slices of one seeds.json, then `whest dataset merge`
    offline once all slices are in blob under bakes/<name>/."""
    check_name(args.name)
    import secrets
    if args.seeds_file:
        with open(args.seeds_file) as f:
            seeds = json.load(f)
    else:
        seeds = [secrets.randbits(63) for _ in range(args.n_mlps)]
    if len(seeds) != args.n_mlps:
        sys.exit(f"seeds file has {len(seeds)} seeds, --n-mlps is {args.n_mlps}")
    e = env()
    key = account_key(e)
    sf = f"{args.name}.seeds.json"
    with open(sf, "w") as f:
        json.dump(seeds, f)
    az("storage", "blob", "upload", "-c", e["CONTAINER_DATASET"], "--account-name", e["STORAGE"], "--account-key", key,
       "-f", sf, "-n", f"bakes/{args.name}/seeds.json", "--overwrite", json_out=False)
    sas_r = container_sas(e, key, e["CONTAINER_DATASET"], "rl", args.sas_hours)
    sas_w = container_sas(e, key, e["CONTAINER_DATASET"], "rcwl", args.sas_hours)
    job = gpu_job(e, args, "bake")
    image = f"{e['ACR']}.azurecr.io/{args.image or e['GPU_IMAGE']}"
    tasks = []
    for k in range(args.slices):
        cmd = ["python", "/app/bake_shard.py", "--name", args.name,
               "--seeds", blob_url(e, e["CONTAINER_DATASET"], f"bakes/{args.name}/seeds.json", sas_r),
               "--n-mlps", str(args.n_mlps), "--n-samples", str(args.n_samples), "--slice", f"{k}/{args.slices}",
               "--split", args.split, "--dest", container_url(e, e["CONTAINER_DATASET"], sas_w)]
        # no --gpus: Batch enables the GPUs for container tasks on GPU pools itself
        tasks.append({"id": f"slice-{k:03d}", "commandLine": "/bin/bash -c " + shlex.quote(shlex.join(cmd)),
                      "containerSettings": {"imageName": image, "containerRunOptions": "--rm --shm-size 8g"},
                      "constraints": {"maxWallClockTime": "PT24H", "maxTaskRetryCount": 1},
                      "userIdentity": TASK_USER})
    n_bad = add_tasks(job, tasks)
    finish_job(job, args.keep_job_open)
    print(f"submitted {len(tasks) - n_bad} bake slices for {args.n_mlps} MLPs x {args.n_samples} samples in job {job}; seeds in {sf}")
    if n_bad:
        sys.exit(f"{n_bad} tasks were not added")


def cmd_atlas(args):
    """Moment atlases (bake_moments.py) for a range of seed-regenerable MLPs on the GPU pool."""
    check_name(args.name)
    e = env()
    key = account_key(e)
    sas_w = container_sas(e, key, e["CONTAINER_DATASET"], "rcwl", args.sas_hours)
    job = gpu_job(e, args, "atlas")
    image = f"{e['ACR']}.azurecr.io/{args.image or e['GPU_IMAGE']}"
    dest = container_url(e, e["CONTAINER_DATASET"], sas_w)
    tasks = []
    for i in range(args.start, args.start + args.count):
        cmd = ["python", "/app/bake_moments.py", "--weights-from", "seeds", "--seed0", str(args.seed0), "--idx", str(i),
               "--n-samples", str(args.n_samples)] + (["--pairs"] if args.pairs else []) + \
              ["--out", dest, "--tag", f"atlas/{args.name}/mlp_{i:05d}"]
        tasks.append({"id": f"mlp-{i:05d}", "commandLine": "/bin/bash -c " + shlex.quote(shlex.join(cmd)),
                      "containerSettings": {"imageName": image, "containerRunOptions": "--rm"},
                      "constraints": {"maxWallClockTime": "PT12H", "maxTaskRetryCount": 1},
                      "userIdentity": TASK_USER})
    n_bad = add_tasks(job, tasks)
    finish_job(job, args.keep_job_open)
    print(f"submitted {len(tasks) - n_bad} atlas tasks in job {job}")
    if n_bad:
        sys.exit(f"{n_bad} tasks were not added")


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--name", required=True)
        p.add_argument("--sas-hours", type=float, default=168.0, help="lifetime of the SAS URLs in the task command lines")
        p.add_argument("--keep-job-open", action="store_true", help="do not auto-terminate the job when its tasks complete")

    s = sub.add_parser("submit"); common(s)
    s.add_argument("--bundle", action="append", required=True)
    s.add_argument("--split", default="full"); s.add_argument("--shards", help="A-B (inclusive) or a single index")
    s.add_argument("--regions"); s.add_argument("--wall-time-limit", type=float, default=600.0)
    s.add_argument("--n-mlps", type=int); s.add_argument("--extra", default="")
    s.add_argument("--task-hours", type=int, default=6)
    s.add_argument("--runner", default="local", choices=["local", "subprocess"],
                   help="local (fast; 8 GB limit advisory) or subprocess (grader transport)")
    b = sub.add_parser("bake"); common(b)
    b.add_argument("--n-mlps", type=int, required=True)
    b.add_argument("--n-samples", type=int, default=100_000_000); b.add_argument("--slices", type=int, default=16)
    b.add_argument("--split", default="dev"); b.add_argument("--seeds-file"); b.add_argument("--region")
    b.add_argument("--image", help="repository:tag in the registry (default: GPU_IMAGE from env.sh)")
    a = sub.add_parser("atlas"); common(a)
    a.add_argument("--start", type=int, default=0); a.add_argument("--count", type=int, required=True)
    a.add_argument("--seed0", type=int, default=770000); a.add_argument("--n-samples", type=int, default=100_000_000)
    a.add_argument("--pairs", action="store_true"); a.add_argument("--region")
    a.add_argument("--image", help="repository:tag in the registry (default: GPU_IMAGE from env.sh)")
    st_ = sub.add_parser("status"); st_.add_argument("--name", required=True)
    st_.add_argument("--failures", action="store_true", help="list failed tasks with exit codes")
    c = sub.add_parser("collect"); c.add_argument("--name", required=True); c.add_argument("--csv")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    {"submit": cmd_submit, "status": cmd_status, "collect": cmd_collect, "bake": cmd_bake, "atlas": cmd_atlas}[args.cmd](args)


if __name__ == "__main__":
    main()
