#!/usr/bin/env python3
"""Submit, watch and collect an experiment grid on the Azure Batch pools created by 04_pools.sh.

    ./grid.py submit  --name exp1 --bundle path/to/estimator_dir [--bundle other.tar.gz ...] \
                      --split full|mini [--shards 0-62] [--wall-time-limit 600] [--n-mlps N]
    ./grid.py status  --name exp1
    ./grid.py collect --name exp1 [--csv exp1.csv]

One task = (bundle, shard).  Tasks are spread round-robin over the regions in env.sh; each
region gets a job named <name>-<region>.  Results land in the `results` container under
<name>/<bundle-tag>/<shard>.json and are aggregated by `collect` (mean of per-MLP adjusted
scores across all shards, i.e. exactly the leaderboard quantity over the shards run).

Uses the `az` CLI only (no SDK), so it runs wherever `az login` works.
"""
import argparse, datetime, json, os, shlex, subprocess, sys, tarfile, tempfile, collections

HERE = os.path.dirname(os.path.abspath(__file__))


def env():
    out = subprocess.run(["bash", "-c", f"source {HERE}/env.sh >/dev/null 2>&1; env"], capture_output=True, text=True).stdout
    e = dict(line.split("=", 1) for line in out.splitlines() if "=" in line)
    return e


def az(*args, json_out=True):
    cmd = ["az", *args]
    if json_out:
        cmd += ["-o", "json"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"az {' '.join(args[:3])} failed: {p.stderr.strip()[:500]}")
    return json.loads(p.stdout) if json_out and p.stdout.strip() else None


def account_key(e):
    keys = az("storage", "account", "keys", "list", "-n", e["STORAGE"], "-g", e["RG"])
    return keys[0]["value"]


def container_sas(e, key, container, perms, hours=72):
    expiry = (datetime.datetime.now(datetime.timezone.utc) + datetime.timedelta(hours=hours)).strftime("%Y-%m-%dT%H:%MZ")
    return az("storage", "container", "generate-sas", "-n", container, "--account-name", e["STORAGE"],
              "--account-key", key, "--permissions", perms, "--expiry", expiry)


def blob_url(e, container, name, sas):
    return f"https://{e['STORAGE']}.blob.core.windows.net/{container}/{name}?{sas}"


def batch_login(e, region):
    acct = subprocess.run(["bash", "-c", f"source {HERE}/env.sh >/dev/null 2>&1; batch_account_name {region}"],
                          capture_output=True, text=True).stdout.strip()
    az("batch", "account", "login", "-n", acct, "-g", e["RG"], json_out=False)
    return acct


def list_shards(e, key, split):
    blobs = az("storage", "blob", "list", "-c", e["CONTAINER_DATASET"], "--account-name", e["STORAGE"],
               "--account-key", key, "--prefix", f"{e['HF_REVISION']}/data/{split}-", "--num-results", "*")
    return sorted(b["name"] for b in blobs if b["name"].endswith(".parquet"))


def upload_bundle(e, key, bundle, name):
    tag = os.path.basename(os.path.normpath(bundle)).replace(".tar.gz", "")
    if os.path.isdir(bundle):
        tmp = tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=False).name
        with tarfile.open(tmp, "w:gz") as tf:
            for root, _, files in os.walk(bundle):
                for f in files:
                    full = os.path.join(root, f)
                    tf.add(full, arcname=os.path.relpath(full, bundle))
        path = tmp
    else:
        path = bundle
    blob = f"{name}/{tag}.tar.gz"
    az("storage", "blob", "upload", "-c", e["CONTAINER_SUBMISSIONS"], "--account-name", e["STORAGE"],
       "--account-key", key, "-f", path, "-n", blob, "--overwrite", json_out=False)
    return tag, blob


def cmd_submit(args):
    e = env()
    key = account_key(e)
    regions = (args.regions or e["REGIONS"]).split()
    shards = list_shards(e, key, args.split)
    if args.shards:
        lo, hi = (int(x) for x in args.shards.split("-"))
        shards = shards[lo:hi + 1]
    if not shards:
        sys.exit("no shards found; run 03_stage_dataset.sh first")
    sas_data = container_sas(e, key, e["CONTAINER_DATASET"], "rl")
    sas_sub = container_sas(e, key, e["CONTAINER_SUBMISSIONS"], "rl")
    sas_res = container_sas(e, key, e["CONTAINER_RESULTS"], "rcwl")
    bundles = [upload_bundle(e, key, b, args.name) for b in args.bundle]
    image = f"{e['ACR']}.azurecr.io/{e['IMAGE']}"
    tasks = []
    for tag, bblob in bundles:
        for s in shards:
            sname = os.path.basename(s).replace(".parquet", "")
            out = blob_url(e, e["CONTAINER_RESULTS"], f"{args.name}/{tag}/{sname}.json", sas_res)
            cmdline = (f"python /app/run_shard.py --bundle {shlex.quote(blob_url(e, e['CONTAINER_SUBMISSIONS'], bblob, sas_sub))} "
                       f"--shard {shlex.quote(blob_url(e, e['CONTAINER_DATASET'], s, sas_data))} --split {args.split} "
                       f"--out {shlex.quote(out)} --tag {tag} --wall-time-limit {args.wall_time_limit} "
                       f"--max-threads {e['VCPU_PER_TASK']}" + (f" --n-mlps {args.n_mlps}" if args.n_mlps else "")
                       + (f" --extra {shlex.quote(args.extra)}" if args.extra else ""))
            tasks.append({"id": f"{tag[:30]}-{sname}"[:64].replace("_", "-"),
                          "commandLine": f"/bin/bash -c {shlex.quote(cmdline)}",
                          "containerSettings": {"imageName": image, "containerRunOptions": "--rm"},
                          "constraints": {"maxWallClockTime": f"PT{args.task_hours}H", "maxTaskRetryCount": 1},
                          "userIdentity": {"autoUser": {"scope": "task", "elevationLevel": "nonadmin"}}})
    manifest = {"name": args.name, "split": args.split, "bundles": [t for t, _ in bundles], "shards": shards,
                "regions": regions, "n_tasks": len(tasks), "jobs": {}, "submitted_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    # round-robin over regions, chunked bulk-add of 100
    per_region = collections.defaultdict(list)
    for i, t in enumerate(tasks):
        per_region[regions[i % len(regions)]].append(t)
    for r, ts in per_region.items():
        batch_login(e, r)
        job = f"{args.name}-{r}"
        pool = f"{e['PREFIX']}-pool-{r}"
        try:
            az("batch", "job", "create", "--id", job, "--pool-id", pool, json_out=False)
        except RuntimeError as ex:
            if "JobExists" not in str(ex):
                raise
        for k in range(0, len(ts), 100):
            with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
                json.dump(ts[k:k + 100], f)
            az("batch", "task", "create", "--job-id", job, "--json-file", f.name, json_out=False)
        manifest["jobs"][r] = {"job": job, "n_tasks": len(ts)}
        print(f"{r}: {len(ts)} tasks in job {job}")
    with open(f"{args.name}.manifest.json", "w") as f:
        json.dump(manifest, f, indent=1)
    print(f"submitted {len(tasks)} tasks ({len(bundles)} bundles x {len(shards)} shards) across {len(per_region)} regions")


def cmd_status(args):
    e = env()
    man = json.load(open(f"{args.name}.manifest.json"))
    total = collections.Counter()
    for r, j in man["jobs"].items():
        batch_login(e, r)
        counts = az("batch", "job", "task-counts", "show", "--job-id", j["job"])
        tc = counts.get("taskCounts", counts)
        row = {k: tc.get(k, 0) for k in ("active", "running", "completed", "succeeded", "failed")}
        total.update(row)
        print(f"{r:<16} {row}")
    print("total", dict(total))


def cmd_collect(args):
    e = env()
    key = account_key(e)
    man = json.load(open(f"{args.name}.manifest.json"))
    tmp = tempfile.mkdtemp()
    az("storage", "blob", "download-batch", "-s", e["CONTAINER_RESULTS"], "-d", tmp, "--pattern", f"{args.name}/*",
       "--account-name", e["STORAGE"], "--account-key", key, json_out=False)
    rows = []
    for tag in man["bundles"]:
        d = os.path.join(tmp, args.name, tag)
        per_mlp = []
        shards_done = 0
        harness_s = 0.0
        for fn in sorted(os.listdir(d)) if os.path.isdir(d) else []:
            res = json.load(open(os.path.join(d, fn)))
            rep = res.get("report")
            harness_s += res["task"].get("harness_seconds", 0)
            if not rep:
                continue
            shards_done += 1
            per_mlp.extend(rep.get("results", rep).get("per_mlp", []))
        if not per_mlp:
            rows.append((tag, shards_done, 0, None, None, None, None, None, None))
            continue
        import statistics as st
        adj = [m.get("adjusted_final_layer_score") for m in per_mlp if m.get("adjusted_final_layer_score") is not None]
        mse = [m.get("final_layer_mse") for m in per_mlp if m.get("final_layer_mse") is not None]
        allm = [m.get("all_layers_mse") for m in per_mlp if m.get("all_layers_mse") is not None]
        fl = [m.get("flops_used") or m.get("effective_compute") or 0 for m in per_mlp]
        res_t = [m.get("residual_wall_time_s") or 0 for m in per_mlp]
        failed = sum(1 for m in per_mlp if m.get("budget_exhausted") or m.get("time_exhausted") or m.get("residual_wall_time_exhausted") or m.get("combined_budget_exhausted") or m.get("error") or m.get("traceback"))
        rows.append((tag, shards_done, len(per_mlp), st.mean(adj) if adj else None, st.mean(mse) if mse else None,
                     st.mean(allm) if allm else None, st.mean(fl) / 2**41 if fl else None, max(res_t) if res_t else None, failed))
    print(f"{'bundle':<32} {'shards':>6} {'mlps':>5} {'adjusted':>11} {'final_mse':>11} {'all_mse':>11} {'util':>7} {'max_resid':>9} {'failed':>6}")
    for r in rows:
        f = lambda x, w: (f"{x:{w}.3e}" if isinstance(x, float) else f"{str(x):>{w}}")
        print(f"{r[0]:<32} {r[1]:>6} {r[2]:>5} {f(r[3],11)} {f(r[4],11)} {f(r[5],11)} {f(r[6],7)} {f(r[7],9)} {r[8]:>6}")
    if args.csv:
        with open(args.csv, "w") as fcsv:
            fcsv.write("bundle,shards,mlps,adjusted,final_mse,all_mse,util,max_residual_s,failed\n")
            for r in rows:
                fcsv.write(",".join("" if x is None else str(x) for x in r) + "\n")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("submit"); s.add_argument("--name", required=True); s.add_argument("--bundle", action="append", required=True)
    s.add_argument("--split", default="full"); s.add_argument("--shards"); s.add_argument("--regions")
    s.add_argument("--wall-time-limit", type=float, default=600.0); s.add_argument("--n-mlps", type=int)
    s.add_argument("--extra", default=""); s.add_argument("--task-hours", type=int, default=6)
    st_ = sub.add_parser("status"); st_.add_argument("--name", required=True)
    c = sub.add_parser("collect"); c.add_argument("--name", required=True); c.add_argument("--csv")
    args = ap.parse_args()
    {"submit": cmd_submit, "status": cmd_status, "collect": cmd_collect}[args.cmd](args)


if __name__ == "__main__":
    main()
