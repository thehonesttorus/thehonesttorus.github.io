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
- Regions: EC2 quotas are per region, so the fleet spans REGIONS (env CLAUDE_AWS_REGIONS=r1,r2,... overrides).
  The bucket stays in the home region (AWS_DEFAULT_REGION, the bucket's); every instance syncs from and uploads to
  it with that --region. Another region gets its own claude-whest-sg (no inbound) in its default VPC, created on
  first launch there (state.json "sg_by_region"); the instance profile is IAM, hence global.

  python3 infra/aws/awsrun.py launch N TYPE [spot] [--region R]   launch N fleet instances (e.g. 2 c7a.48xlarge)
  python3 infra/aws/awsrun.py list                     fleet instances in all regions, state, readiness
  python3 infra/aws/awsrun.py start|stop all|NAME
  python3 infra/aws/awsrun.py terminate NAME
  python3 infra/aws/awsrun.py data LOCALFILE... [--k3work]   upload data files (default data/official/)
  python3 infra/aws/awsrun.py batch JOB FILE [--threads T] [--only NAME,...] [--slots K] [--mem GB]
                                                             tasks "tag<TAB>command" (one per line); returns at once.
                                                             Slots = min(cores / T, 0.94 memory / GB), GB default 10
  python3 infra/aws/awsrun.py watch JOB                      stream results as they arrive; summary.txt at the end
  python3 infra/aws/awsrun.py get JOB [PATTERN]              download results/JOB/ to scratchpad/aws/results/JOB
"""
import base64, fnmatch, functools, hashlib, io, json, os, sys, tarfile, time
from concurrent.futures import ThreadPoolExecutor
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
AMI_PARAM = "/aws/service/canonical/ubuntu/server/24.04/stable/current/amd64/hvm/ebs-gp3/ami-id"
S = boto3.session.Session()
REGION = S.region_name  # home region: the bucket and state.json's "sg" live here
REGIONS = os.environ.get("CLAUDE_AWS_REGIONS", "eu-north-1,eu-west-1,eu-central-1,eu-west-2,"
                         "us-east-1,us-east-2,us-west-2").split(",")
_CLIENTS = {}


def cl(service, region=REGION):
    """boto3 client for service in region, cached (create them in the main thread: sessions are not thread-safe)"""
    if (service, region) not in _CLIENTS:
        _CLIENTS[service, region] = S.client(service, region_name=region)
    return _CLIENTS[service, region]


s3 = cl("s3")


def state():
    return json.load(open(os.path.join(SCRATCH, "state.json")))


def save(st):
    p = os.path.join(SCRATCH, "state.json")
    json.dump(st, open(p + ".tmp", "w"), indent=1); os.replace(p + ".tmp", p)


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
    """fleet instances in all REGIONS (queried in parallel); each carries its region as x["Region"]"""
    f = [{"Name": "tag:Project", "Values": [TAG["Value"]]}, {"Name": "tag:Fleet", "Values": [FLEET["Value"]]},
         {"Name": "instance-state-name", "Values": list(states)}]
    clients = [(r, cl("ec2", r)) for r in REGIONS]
    out = []
    with ThreadPoolExecutor(len(clients)) as pool:
        for reg, rs in pool.map(lambda rc: (rc[0], rc[1].describe_instances(Filters=f)["Reservations"]), clients):
            for r in rs:
                out += [dict(x, Region=reg) for x in r["Instances"]]
    return sorted(out, key=lambda x: (len(_name(x)), _name(x)))  # fleet2 before fleet10


def _name(x):
    return next((t["Value"] for t in x.get("Tags", []) if t["Key"] == "Name"), "?")


@functools.lru_cache(maxsize=None)
def _itype(itype, region=REGION):
    return cl("ec2", region).describe_instance_types(InstanceTypes=[itype])["InstanceTypes"][0]


def _vcpus(itype, region=REGION):
    return _itype(itype, region)["VCpuInfo"]["DefaultVCpus"]


def _mem_gb(itype, region=REGION):
    return _itype(itype, region)["MemoryInfo"]["SizeInMiB"] / 1024


def _sg(region):
    """the security group for region: state.json's "sg" at home, else claude-whest-sg in that region's default VPC
    (no inbound rules; SSM needs only the default all-egress rule), found or created once, kept in "sg_by_region"."""
    st = state()
    if region == REGION:
        return st["sg"]
    if region not in st.get("sg_by_region", {}):
        e = cl("ec2", region)
        vpc = e.describe_vpcs(Filters=[{"Name": "isDefault", "Values": ["true"]}])["Vpcs"][0]["VpcId"]
        old = e.describe_security_groups(Filters=[{"Name": "group-name", "Values": ["claude-whest-sg"]},
                                                  {"Name": "vpc-id", "Values": [vpc]},
                                                  {"Name": "tag:Project", "Values": [TAG["Value"]]}])["SecurityGroups"]
        sg = old[0]["GroupId"] if old else e.create_security_group(
            GroupName="claude-whest-sg", Description="claude-whest runner: no inbound", VpcId=vpc,
            TagSpecifications=[{"ResourceType": "security-group", "Tags": [TAG]}])["GroupId"]
        st = state(); st.setdefault("sg_by_region", {})[region] = sg; save(st)
        print(f"security group {sg} in {region} ({vpc})")
    return st["sg_by_region"][region]


def launch(n, itype, spot=False, region=None):
    region = region or REGION
    if region not in REGIONS:
        sys.exit(f"region {region} is not in REGIONS {REGIONS} (set CLAUDE_AWS_REGIONS)")
    st = state(); sg = _sg(region)
    ami = cl("ssm", region).get_parameter(Name=AMI_PARAM)["Parameter"]["Value"]
    have = {_name(x) for x in fleet()}
    for _ in range(n):
        i = 1
        while f"fleet{i}" in have:
            i += 1
        name = f"fleet{i}"; have.add(name)
        kw = dict(ImageId=ami, InstanceType=itype, MinCount=1, MaxCount=1, UserData=boot_script(st["bucket"]),
                  IamInstanceProfile={"Name": ROLE}, SecurityGroupIds=[sg],
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
            r = cl("ec2", region).run_instances(**kw)["Instances"][0]
        except Exception as e:
            print(f"{name}: launch refused in {region}: {str(e)[:400]}"); return
        print(f"launched {name}: {r['InstanceId']} {itype} {'spot' if spot else 'on-demand'} in {region}")


def _ssm(xs, script, timeout=600, wait=True):
    """run script on fleet instances xs (one send_command per region); returns {region: command id}, {iid: result}"""
    cmds = {reg: cl("ssm", reg).send_command(
                InstanceIds=[x["InstanceId"] for x in xs if x["Region"] == reg], DocumentName="AWS-RunShellScript",
                Parameters={"commands": [script], "executionTimeout": [str(timeout)]})["Command"]["CommandId"]
            for reg in sorted({x["Region"] for x in xs})}
    if not wait:
        return cmds, {}
    res = {}
    for _ in range(timeout // 3 + 20):
        time.sleep(3)
        for x in xs:
            iid, c = x["InstanceId"], cl("ssm", x["Region"])
            if iid in res:
                continue
            try:
                inv = c.get_command_invocation(CommandId=cmds[x["Region"]], InstanceId=iid)
            except c.exceptions.InvocationDoesNotExist:
                continue
            if inv["Status"] not in ("Pending", "InProgress", "Delayed"):
                res[iid] = (inv["Status"], inv.get("StandardOutputContent", ""), inv.get("StandardErrorContent", ""))
        if len(res) == len(xs):
            break
    return cmds, res


def show():
    xs = fleet()
    online = {}
    run = [x for x in xs if x["State"]["Name"] == "running"]
    for reg in sorted({x["Region"] for x in run}):
        ids = [x["InstanceId"] for x in run if x["Region"] == reg]
        for i in cl("ssm", reg).describe_instance_information(
                Filters=[{"Key": "InstanceIds", "Values": ids}])["InstanceInformationList"]:
            online[i["InstanceId"]] = i.get("PingStatus")
    up = [x for x in run if online.get(x["InstanceId"]) == "Online"]
    _, res = _ssm(up, "test -f /opt/data/READY && echo READY || echo booting; pgrep -fc '/opt/run/.*/task.sh' || true",
                  timeout=60) if up else (None, {})
    for x in xs:
        iid = x["InstanceId"]; r = res.get(iid)
        extra = " ".join(r[1].split()) if r else ""
        print(f"{_name(x):8s} {x['Region']:14s} {iid} {x['InstanceType']:14s} {x['State']['Name']:9s} "
              f"ssm={online.get(iid, '-'):7s} {extra}")


def power(action, which):
    xs = fleet()
    sel = xs if which == "all" else [x for x in xs if _name(x) == which]
    if not sel:
        sys.exit(f"no fleet instance {which}")
    if action == "terminate" and which == "all":
        sys.exit("terminate one instance at a time")
    for reg in sorted({x["Region"] for x in sel}):  # one call per region; a refusal (e.g. quota) skips only that region
        e = cl("ec2", reg); xr = [x for x in sel if x["Region"] == reg]
        try:
            {"start": e.start_instances, "stop": e.stop_instances, "terminate": e.terminate_instances}[action](
                InstanceIds=[x["InstanceId"] for x in xr])
            print(f"{action}: " + " ".join(_name(x) for x in xr) + f" ({reg})")
        except Exception as ex:
            print(f"{action} refused in {reg} for {' '.join(_name(x) for x in xr)}: {str(ex)[:400]}")


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


def batch(job, path, threads=4, only=None, nslots=None, mem=10.0):
    st = state(); bkt = st["bucket"]
    tasks = [l.rstrip("\n").split("\t", 1) for l in open(path) if l.strip() and not l.startswith("#")]
    run = [x for x in fleet(("running",)) if only is None or _name(x) in only]
    if not run:
        sys.exit("no running fleet instance (awsrun.py start all, or launch)")
    code = bundle()
    inst = {x["InstanceId"]: x for x in run}
    cores = {x["InstanceId"]: _vcpus(x["InstanceType"], x["Region"]) for x in run}
    gbs = {x["InstanceId"]: _mem_gb(x["InstanceType"], x["Region"]) for x in run}
    # slots: cores / threads, capped by memory (--mem GB per task, 6% kept for the system); --slots K overrides.
    # Weighted round robin: instance i gets tasks in proportion to its slots.
    slots = {i: max(1, min(c // threads, int(0.94 * gbs[i] // mem))) if nslots is None else nslots
             for i, c in cores.items()}
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
        _ssm([inst[iid]], script, timeout=1200, wait=False)
        print(f"{job}: {len(mine)} tasks -> {_name(inst[iid])} {iid} {inst[iid]['Region']} ({slots[iid]} slots)")
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
    if c == "launch": launch(int(a[1]), a[2], spot=("spot" in a[3:]),
                             region=a[a.index("--region") + 1] if "--region" in a else None)
    elif c == "list": show()
    elif c in ("start", "stop", "terminate"): power(c, a[1])
    elif c == "data": data([x for x in a[1:] if not x.startswith("--")], "k3work" if "--k3work" in a else "official")
    elif c == "batch": batch(a[1], a[2], int(a[a.index("--threads") + 1]) if "--threads" in a else 4,
                             a[a.index("--only") + 1].split(",") if "--only" in a else None,
                             int(a[a.index("--slots") + 1]) if "--slots" in a else None,
                             float(a[a.index("--mem") + 1]) if "--mem" in a else 10.0)
    elif c == "watch": watch(a[1])
    elif c == "get": get(a[1], a[2] if len(a) > 2 else "*")
    else: sys.exit(__doc__)
