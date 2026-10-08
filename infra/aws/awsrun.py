#!/usr/bin/env python3
"""AWS fleet runner: many 4-core jobs per large instance, logs streamed through S3 (beside awsjob.py).

Credentials come from the environment (AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY, AWS_DEFAULT_REGION); nothing here
prints them. Every instance is tagged Project=claude-whest and Fleet=whest; nothing else is touched. The bucket, the
security group and the instance role are awsjob.py's (scratchpad/aws/state.json).

Design:
- Data: s3://BUCKET/data/{official,k3work}/ is synced to /opt/data on every instance at boot. The official networks
  and truth files are there; Monte Carlo files go to data/k3work when a job needs them.
- Code: every batch carries a content-addressed bundle of workbench/{k3work,official,num12,ncgprob}/*.py
  (s3://BUCKET/bundle/<hash>.tgz). Each task runs in its own directory, with its own copy of the code and symlinks
  to the data, so concurrent tasks and concurrent code versions never collide.
- Scheduling: tasks are spread over the running fleet in proportion to its cores. Each instance runs its share with
  xargs -P (cores / threads) slots. Each task uploads log_<tag>.txt and its new files to s3://BUCKET/results/JOB/ the
  moment it finishes, and `watch` prints them as they arrive.
- Cost: instances stop themselves after IDLE minutes with no task running (a systemd watchdog; stop, not
  terminate, so a restart is quick).

  python3 infra/aws/awsrun.py launch N TYPE [spot]     launch N fleet instances (e.g. 2 c7a.48xlarge)
  python3 infra/aws/awsrun.py list                     fleet instances, state, cores, readiness
  python3 infra/aws/awsrun.py start|stop all|NAME
  python3 infra/aws/awsrun.py terminate NAME
  python3 infra/aws/awsrun.py data LOCALFILE... [--k3work]   upload data files (default data/official/)
  python3 infra/aws/awsrun.py batch JOB FILE [--threads T] [--only NAME,...] [--slots K]
                                                             tasks "tag<TAB>command" (one per line); returns at once
  python3 infra/aws/awsrun.py watch JOB                      stream results as they arrive; summary.txt at the end
  python3 infra/aws/awsrun.py get JOB [PATTERN]              download results/JOB/ to scratchpad/aws/results/JOB
"""
import base64, fnmatch, hashlib, io, json, os, sys, tarfile, time
import boto3

TAG = {"Key": "Project", "Value": "claude-whest"}
FLEET = {"Key": "Fleet", "Value": "whest"}
HERE = os.path.dirname(os.path.abspath(__file__))
WB = os.path.join(HERE, "..", "..", "workbench")
CODE_DIRS = ("k3work", "official", "num12", "ncgprob")
SCRATCH = os.environ.get("CLAUDE_AWS_DIR", "/tmp/claude-0/-home-user-thehonesttorus-github-io/"
                         "2cab30f0-c454-5769-8d9a-8195cd80a5f3/scratchpad/aws")
ROLE = "claude-ec2-runner"
IDLE_MIN = 30
S = boto3.session.Session()
REGION = S.region_name
ec2, s3, ssm = S.client("ec2"), S.client("s3"), S.client("ssm")


def state():
    return json.load(open(os.path.join(SCRATCH, "state.json")))


def boot_script(bucket):
    return f"""#!/bin/bash
set -eux
export DEBIAN_FRONTEND=noninteractive
apt-get update -y
apt-get install -y python3-venv python3-pip unzip
python3 -m venv /opt/venv
/opt/venv/bin/pip install -q --upgrade pip
/opt/venv/bin/pip install -q numpy scipy pyarrow flopscope==0.12.1 whestbench==0.16.1
curl -sSL https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip -o /tmp/awscli.zip
unzip -q /tmp/awscli.zip -d /tmp && /tmp/aws/install
mkdir -p /opt/data/official /opt/data/k3work /opt/run && chmod -R 777 /opt/data /opt/run
aws s3 sync s3://{bucket}/data/ /opt/data/ --only-show-errors --region {REGION}
cat > /usr/local/bin/idlewatch.sh <<'EOW'
#!/bin/bash
touch /opt/run/.last
while true; do
  sleep 60
  if pgrep -f "/opt/run/.*/task.sh" >/dev/null; then touch /opt/run/.last; fi
  if [ $(( $(date +%s) - $(stat -c %Y /opt/run/.last) )) -gt {IDLE_MIN * 60} ]; then shutdown -h now; fi
done
EOW
chmod +x /usr/local/bin/idlewatch.sh
cat > /etc/systemd/system/idlewatch.service <<'EOS'
[Unit]
Description=stop the instance when no task has run for a while
[Service]
ExecStartPre=/bin/touch /opt/run/.last
ExecStart=/usr/local/bin/idlewatch.sh
Restart=always
[Install]
WantedBy=multi-user.target
EOS
systemctl daemon-reload && systemctl enable --now idlewatch.service
nproc > /opt/data/nproc.txt
touch /opt/data/READY
"""


def fleet(states=("pending", "running", "stopping", "stopped")):
    f = [{"Name": "tag:Project", "Values": [TAG["Value"]]}, {"Name": "tag:Fleet", "Values": [FLEET["Value"]]},
         {"Name": "instance-state-name", "Values": list(states)}]
    out = []
    for r in ec2.describe_instances(Filters=f)["Reservations"]:
        out += r["Instances"]
    return sorted(out, key=lambda x: _name(x))


def _name(x):
    return next((t["Value"] for t in x.get("Tags", []) if t["Key"] == "Name"), "?")


def _vcpus(itype):
    return ec2.describe_instance_types(InstanceTypes=[itype])["InstanceTypes"][0]["VCpuInfo"]["DefaultVCpus"]


def launch(n, itype, spot=False):
    st = state()
    ami = ssm.get_parameter(Name="/aws/service/canonical/ubuntu/server/24.04/stable/current/amd64/hvm/ebs-gp3/ami-id")["Parameter"]["Value"]
    have = {_name(x) for x in fleet()}
    for _ in range(n):
        i = 1
        while f"fleet{i}" in have:
            i += 1
        name = f"fleet{i}"; have.add(name)
        kw = dict(ImageId=ami, InstanceType=itype, MinCount=1, MaxCount=1, UserData=boot_script(st["bucket"]),
                  IamInstanceProfile={"Name": ROLE}, SecurityGroupIds=[st["sg"]],
                  InstanceInitiatedShutdownBehavior="stop",
                  BlockDeviceMappings=[{"DeviceName": "/dev/sda1", "Ebs": {"VolumeSize": 64, "VolumeType": "gp3",
                                                                            "Throughput": 500, "Iops": 6000,
                                                                            "DeleteOnTermination": True}}],
                  MetadataOptions={"HttpTokens": "required"},
                  TagSpecifications=[{"ResourceType": t, "Tags": [TAG, FLEET, {"Key": "Name", "Value": name}]}
                                     for t in ("instance", "volume")])
        if spot:
            kw["InstanceMarketOptions"] = {"MarketType": "spot", "SpotOptions": {"SpotInstanceType": "persistent",
                                                                                 "InstanceInterruptionBehavior": "stop"}}
        try:
            r = ec2.run_instances(**kw)["Instances"][0]
        except Exception as e:
            print(f"{name}: launch refused: {str(e)[:400]}"); return
        print(f"launched {name}: {r['InstanceId']} {itype} {'spot' if spot else 'on-demand'}")


def _ssm(iids, script, timeout=600, wait=True):
    c = ssm.send_command(InstanceIds=iids, DocumentName="AWS-RunShellScript",
                         Parameters={"commands": [script], "executionTimeout": [str(timeout)]})["Command"]["CommandId"]
    if not wait:
        return c, {}
    res = {}
    for _ in range(timeout // 3 + 20):
        time.sleep(3)
        for iid in iids:
            if iid in res:
                continue
            try:
                inv = ssm.get_command_invocation(CommandId=c, InstanceId=iid)
            except ssm.exceptions.InvocationDoesNotExist:
                continue
            if inv["Status"] not in ("Pending", "InProgress", "Delayed"):
                res[iid] = (inv["Status"], inv.get("StandardOutputContent", ""), inv.get("StandardErrorContent", ""))
        if len(res) == len(iids):
            break
    return c, res


def show():
    xs = fleet()
    online = {}
    run = [x["InstanceId"] for x in xs if x["State"]["Name"] == "running"]
    if run:
        for i in ssm.describe_instance_information(Filters=[{"Key": "InstanceIds", "Values": run}])["InstanceInformationList"]:
            online[i["InstanceId"]] = i.get("PingStatus")
        _, res = _ssm([i for i in run if online.get(i) == "Online"] or run[:0],
                      "test -f /opt/data/READY && echo READY || echo booting; pgrep -fc '/opt/run/.*/task.sh' || true",
                      timeout=60) if any(online.get(i) == "Online" for i in run) else (None, {})
    else:
        res = {}
    for x in xs:
        iid = x["InstanceId"]; r = res.get(iid)
        extra = " ".join(r[1].split()) if r else ""
        print(f"{_name(x):8s} {iid} {x['InstanceType']:14s} {x['State']['Name']:9s} ssm={online.get(iid, '-'):7s} {extra}")


def power(action, which):
    xs = fleet()
    sel = xs if which == "all" else [x for x in xs if _name(x) == which]
    if not sel:
        sys.exit(f"no fleet instance {which}")
    ids = [x["InstanceId"] for x in sel]
    if action == "terminate" and which == "all":
        sys.exit("terminate one instance at a time")
    {"start": ec2.start_instances, "stop": ec2.stop_instances, "terminate": ec2.terminate_instances}[action](InstanceIds=ids)
    print(f"{action}: " + " ".join(_name(x) for x in sel))


def data(files, sub="official"):
    b = state()["bucket"]
    for f in files:
        s3.upload_file(f, b, f"data/{sub}/{os.path.basename(f)}")
        print(f"uploaded {os.path.basename(f)} -> data/{sub}/")


def bundle():
    buf = io.BytesIO(); h = hashlib.sha256(); members = []
    for d in CODE_DIRS:
        p = os.path.join(WB, d)
        for f in sorted(os.listdir(p)):
            if f.endswith(".py"):
                b = open(os.path.join(p, f), "rb").read(); h.update(d.encode()); h.update(f.encode()); h.update(b)
                members.append((f"{d}/{f}", b))
    code = h.hexdigest()[:20]
    key = f"bundle/{code}.tgz"; bkt = state()["bucket"]
    try:
        s3.head_object(Bucket=bkt, Key=key)
    except Exception:
        with tarfile.open(fileobj=buf, mode="w:gz") as tf:
            for name, b in members:
                ti = tarfile.TarInfo(name); ti.size = len(b); tf.addfile(ti, io.BytesIO(b))
        s3.put_object(Bucket=bkt, Key=key, Body=buf.getvalue())
        print(f"uploaded code bundle {code} ({len(members)} files)")
    return code


TASK_SH = r"""#!/bin/bash
# one task: $1 = tag. Own directory, own code copy, symlinked data; uploads its log and new files when done.
set -u
J=__JDIR__; T="$1"; W="$J/w/$T"; B=__BUCKET__; R=__REGION__; TH=__TH__; JOB=__JOB__
mkdir -p "$W"
for d in k3work official num12 ncgprob; do mkdir -p "$W/$d"; cp "$J/code/$d/"*.py "$W/$d/" 2>/dev/null || true; done
for sub in official k3work; do for f in /opt/data/$sub/*; do [ -e "$f" ] && ln -sf "$f" "$W/$sub/" ; done; done
cd "$W/k3work"; ls -A > "$W/.before"; mkdir -p "$W/out"
export PATH=/opt/venv/bin:$PATH OUT="$W/out" PYTHONPATH="$W/num12" PYTHONWARNINGS=ignore \
       OMP_NUM_THREADS=$TH OPENBLAS_NUM_THREADS=$TH MKL_NUM_THREADS=$TH
t0=$(date +%s)
bash "$J/tasks/$T.cmd" > "$W/so.txt" 2> "$W/se.txt"; rc=$?
{ cat "$W/so.txt"; if [ -s "$W/se.txt" ]; then echo; echo "[stderr]"; tail -c 6000 "$W/se.txt"; fi
  echo; echo "[exit $rc after $(( $(date +%s) - t0 ))s on $(hostname), $(nproc) cores, threads $TH]"; } > "$W/log_$T.txt"
for f in $(ls -A); do
  if [ -f "$f" ] && [ ! -L "$f" ] && [[ "$f" != *.py ]] && ! grep -qxF "$f" "$W/.before"; then
    aws s3 cp "$f" "s3://$B/results/$JOB/$f" --quiet --region $R; fi
done
for f in $(ls -A "$W/out" 2>/dev/null); do aws s3 cp "$W/out/$f" "s3://$B/results/$JOB/$f" --quiet --region $R; done
aws s3 cp "$W/log_$T.txt" "s3://$B/results/$JOB/log_$T.txt" --quiet --region $R
rm -rf "$W"
"""


def batch(job, path, threads=4, only=None, nslots=None):
    st = state(); bkt = st["bucket"]
    tasks = [l.rstrip("\n").split("\t", 1) for l in open(path) if l.strip() and not l.startswith("#")]
    run = [x for x in fleet(("running",)) if only is None or _name(x) in only]
    if not run:
        sys.exit("no running fleet instance (awsrun.py start all, or launch)")
    code = bundle()
    cores = {x["InstanceId"]: _vcpus(x["InstanceType"]) for x in run}
    # weighted round robin: instance i gets tasks in proportion to its slots (--slots K caps them, e.g. for timing)
    slots = {i: max(1, c // threads) if nslots is None else nslots for i, c in cores.items()}
    order = sorted(slots, key=lambda i: -slots[i]); share = {i: [] for i in order}; load = {i: 0.0 for i in order}
    for t in tasks:
        i = min(order, key=lambda k: (load[k] + 1) / slots[k]); share[i].append(t); load[i] += 1
    for iid, mine in share.items():
        if not mine:
            continue
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:gz") as tf:
            for tag, cmd in mine:
                b = (cmd + "\n").encode(); ti = tarfile.TarInfo(f"tasks/{tag}.cmd"); ti.size = len(b); tf.addfile(ti, io.BytesIO(b))
            sh = (TASK_SH.replace("__JDIR__", f"/opt/run/{job}").replace("__BUCKET__", bkt).replace("__REGION__", REGION)
                  .replace("__TH__", str(threads)).replace("__JOB__", job)).encode()
            ti = tarfile.TarInfo("task.sh"); ti.size = len(sh); ti.mode = 0o755; tf.addfile(ti, io.BytesIO(sh))
            tags = ("\n".join(t for t, _ in mine) + "\n").encode()
            ti = tarfile.TarInfo("tags.txt"); ti.size = len(tags); tf.addfile(ti, io.BytesIO(tags))
        s3.put_object(Bucket=bkt, Key=f"jobs/{job}/{iid}.tgz", Body=buf.getvalue())
        script = f"""set -u
for i in $(seq 90); do [ -f /opt/data/READY ] && break; sleep 10; done
[ -f /opt/data/READY ] || {{ echo "bootstrap not finished"; exit 1; }}
J=/opt/run/{job}; mkdir -p $J/code && cd $J
aws s3 cp s3://{bkt}/bundle/{code}.tgz code.tgz --quiet --region {REGION} && tar xzf code.tgz -C code
aws s3 cp s3://{bkt}/jobs/{job}/{iid}.tgz job.tgz --quiet --region {REGION} && tar xzf job.tgz
aws s3 sync s3://{bkt}/data/ /opt/data/ --only-show-errors --region {REGION}
touch /opt/run/.last
nohup bash -c "xargs -P {slots[iid]} -I{{}} $J/task.sh {{}} < $J/tags.txt" > $J/xargs.log 2>&1 &
echo "{job}: {len(mine)} tasks, {slots[iid]} slots x {threads} threads on $(hostname) ($(nproc) cores)"
"""
        _ssm([iid], script, timeout=1200, wait=False)
        print(f"{job}: {len(mine)} tasks -> {iid} ({slots[iid]} slots)")
    os.makedirs(os.path.join(SCRATCH, "results", job), exist_ok=True)
    json.dump([t for t, _ in tasks], open(os.path.join(SCRATCH, "results", job, "tags.json"), "w"))


def _line(tag, log):
    body = log.split("\n[stderr]\n")[0]
    rc = next((l.split()[1] for l in log.splitlines() if l.startswith("[exit ")), "?")
    tail = [l for l in body.splitlines() if l.strip() and not l.startswith("[exit ")][-2:]
    if rc != "0":
        tail = [l for l in log.splitlines() if l.strip()][-4:]
    return f"[{tag} exit {rc}] " + " | ".join(tail)


def watch(job, timeout=7200):
    bkt = state()["bucket"]; d = os.path.join(SCRATCH, "results", job); os.makedirs(d, exist_ok=True)
    tags = json.load(open(os.path.join(d, "tags.json")))
    seen = {}; t0 = time.time()
    while len(seen) < len(tags) and time.time() - t0 < timeout:
        pag = s3.get_paginator("list_objects_v2")
        for page in pag.paginate(Bucket=bkt, Prefix=f"results/{job}/log_"):
            for o in page.get("Contents", []):
                tag = o["Key"].rsplit("/log_", 1)[1][:-4]
                if tag in seen or tag not in tags:
                    continue
                log = s3.get_object(Bucket=bkt, Key=o["Key"])["Body"].read().decode()
                open(os.path.join(d, f"log_{tag}.txt"), "w").write(log)
                seen[tag] = _line(tag, log); print(seen[tag], f"(+{time.time() - t0:.0f}s)", flush=True)
        if len(seen) < len(tags):
            time.sleep(10)
    open(os.path.join(d, "summary.txt"), "w").write("\n".join(seen[t] for t in tags if t in seen) + "\n")
    print(f"{job}: {len(seen)}/{len(tags)} results in {time.time() - t0:.0f}s", flush=True)


def get(job, pattern="*"):
    bkt = state()["bucket"]; d = os.path.join(SCRATCH, "results", job); os.makedirs(d, exist_ok=True); n = 0
    for page in s3.get_paginator("list_objects_v2").paginate(Bucket=bkt, Prefix=f"results/{job}/"):
        for o in page.get("Contents", []):
            name = o["Key"].split("/", 2)[2]
            if fnmatch.fnmatch(name, pattern):
                s3.download_file(bkt, o["Key"], os.path.join(d, name)); n += 1
    print(f"{job}: {n} files -> {d}")


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a:
        sys.exit(__doc__)
    c = a[0]
    if c == "launch": launch(int(a[1]), a[2], spot=("spot" in a[3:]))
    elif c == "list": show()
    elif c in ("start", "stop", "terminate"): power(c, a[1])
    elif c == "data": data([x for x in a[1:] if not x.startswith("--")], "k3work" if "--k3work" in a else "official")
    elif c == "batch": batch(a[1], a[2], int(a[a.index("--threads") + 1]) if "--threads" in a else 4,
                             a[a.index("--only") + 1].split(",") if "--only" in a else None,
                             int(a[a.index("--slots") + 1]) if "--slots" in a else None)
    elif c == "watch": watch(a[1])
    elif c == "get": get(a[1], a[2] if len(a) > 2 else "*")
    else: sys.exit(__doc__)
