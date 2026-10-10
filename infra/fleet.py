#!/usr/bin/env python3
"""Fleet runner for this branch (derived from the earlier infra/aws/awsrun.py design; same bucket, security group and
instance role, but its own Fleet tag so it never collides with other sessions' instances).

  python3 infra/fleet.py launch N TYPE [spot] [--region R]   e.g. launch 1 c7i.48xlarge spot
  python3 infra/fleet.py list
  python3 infra/fleet.py start|stop all|NAME ; terminate NAME
  python3 infra/fleet.py batch JOB FILE [--threads T] [--slots K] [--mem GB] [--only NAME,..]
        FILE lines: "tag<TAB>command"; the command runs in a private work dir with this repo's whest/ and scripts/
        copied in, DATA=/opt/data/official (W_off*.npy, truth_off*.npz), OUT=$W/out; new files and the log go to
        s3://BUCKET/results9/JOB/.
  python3 infra/fleet.py watch JOB            stream logs and output files as each task finishes (3 s polls)
  python3 infra/fleet.py ready NAME           wait for running + SSM online + data bootstrap
  python3 infra/fleet.py get JOB [PATTERN]    download results to SCRATCH/results/JOB
  python3 infra/fleet.py ssm NAME 'CMD'       run a shell command on one instance (diagnostics)
"""
import base64, fnmatch, functools, hashlib, io, json, os, sys, tarfile, time
from concurrent.futures import ThreadPoolExecutor
import boto3

TAG = {"Key": "Project", "Value": "claude-whest"}
FLEET = {"Key": "Fleet", "Value": "whest9"}
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
CODE_DIRS = ("whest", "scripts")
SCRATCH = os.environ.get("CLAUDE_FLEET_DIR", os.path.expanduser("~/.whest9"))
BUCKET = "claude-whest-9292-97a992"; SG_HOME = "sg-0eb5e94c0923d2de3"; ROLE = "claude-ec2-runner"
IDLE_MIN = 45
AMI_PARAM = "/aws/service/canonical/ubuntu/server/24.04/stable/current/amd64/hvm/ebs-gp3/ami-id"
S = boto3.session.Session(); REGION = S.region_name or "eu-north-1"
REGIONS = os.environ.get("CLAUDE_AWS_REGIONS", "eu-north-1,eu-west-1,eu-central-1,us-east-1,us-east-2,us-west-2").split(",")
_C = {}
def cl(service, region=REGION):
    if (service, region) not in _C:
        _C[service, region] = S.client(service, region_name=region)
    return _C[service, region]
s3 = cl("s3")
def _state_path(): return os.path.join(SCRATCH, "state.json")
def state():
    p = _state_path(); return json.load(open(p)) if os.path.exists(p) else {"sg_by_region": {REGION: SG_HOME}}
def save(st):
    os.makedirs(SCRATCH, exist_ok=True); json.dump(st, open(_state_path(), "w"), indent=1)

def boot_script():
    return f"""#!/bin/bash
set -eux
export DEBIAN_FRONTEND=noninteractive
apt-get update -y
apt-get install -y python3-venv python3-pip unzip
python3 -m venv /opt/venv
/opt/venv/bin/pip install -q --upgrade pip
/opt/venv/bin/pip install -q numpy scipy pyarrow boto3 flopscope==0.12.1 whestbench==0.16.1
curl -sSL https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip -o /tmp/awscli.zip
unzip -q /tmp/awscli.zip -d /tmp && /tmp/aws/install
mkdir -p /opt/data/official /opt/run && chmod -R 777 /opt/data /opt/run
aws s3 sync s3://{BUCKET}/data/official/ /opt/data/official/ --only-show-errors --region {REGION}
cat > /usr/local/bin/idlewatch.sh <<'EOW'
#!/bin/bash
touch /opt/run/.last
while true; do
  sleep 60
  if pgrep -f "/opt/run/.*/task.sh" >/dev/null; then touch /opt/run/.last; fi
  pgrep -f "modal_app.py batch" > /dev/null && touch /opt/run/.last
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
    clients = [(r, cl("ec2", r)) for r in REGIONS]; out = []
    with ThreadPoolExecutor(len(clients)) as pool:
        for reg, rs in pool.map(lambda rc: (rc[0], rc[1].describe_instances(Filters=f)["Reservations"]), clients):
            for r in rs:
                out += [dict(x, Region=reg) for x in r["Instances"]]
    return sorted(out, key=lambda x: (len(_name(x)), _name(x)))
def _name(x): return next((t["Value"] for t in x.get("Tags", []) if t["Key"] == "Name"), "?")
@functools.lru_cache(maxsize=None)
def _itype(itype, region=REGION):
    return cl("ec2", region).describe_instance_types(InstanceTypes=[itype])["InstanceTypes"][0]
def _vcpus(itype, region=REGION): return _itype(itype, region)["VCpuInfo"]["DefaultVCpus"]
def _mem_gb(itype, region=REGION): return _itype(itype, region)["MemoryInfo"]["SizeInMiB"] / 1024
def _sg(region):
    st = state()
    if region not in st.get("sg_by_region", {}):
        e = cl("ec2", region)
        vpc = e.describe_vpcs(Filters=[{"Name": "isDefault", "Values": ["true"]}])["Vpcs"][0]["VpcId"]
        old = e.describe_security_groups(Filters=[{"Name": "group-name", "Values": ["claude-whest-sg"]}, {"Name": "vpc-id", "Values": [vpc]}])["SecurityGroups"]
        sg = old[0]["GroupId"] if old else e.create_security_group(GroupName="claude-whest-sg", Description="claude-whest runner: no inbound", VpcId=vpc,
                                                                  TagSpecifications=[{"ResourceType": "security-group", "Tags": [TAG]}])["GroupId"]
        st.setdefault("sg_by_region", {})[region] = sg; save(st)
    return st["sg_by_region"][region]

def launch(n, itype, spot=False, region=None):
    region = region or REGION; sg = _sg(region)
    ami = cl("ssm", region).get_parameter(Name=AMI_PARAM)["Parameter"]["Value"]
    have = {_name(x) for x in fleet()}
    for _ in range(n):
        i = 1
        while f"w9-{i}" in have: i += 1
        name = f"w9-{i}"; have.add(name)
        kw = dict(ImageId=ami, InstanceType=itype, MinCount=1, MaxCount=1, UserData=boot_script(),
                  IamInstanceProfile={"Name": ROLE}, SecurityGroupIds=[sg], InstanceInitiatedShutdownBehavior="stop",
                  BlockDeviceMappings=[{"DeviceName": "/dev/sda1", "Ebs": {"VolumeSize": 64, "VolumeType": "gp3", "Throughput": 500, "Iops": 6000, "DeleteOnTermination": True}}],
                  MetadataOptions={"HttpTokens": "required"},
                  TagSpecifications=[{"ResourceType": t, "Tags": [TAG, FLEET, {"Key": "Name", "Value": name}]} for t in ("instance", "volume")])
        if spot:
            kw["InstanceMarketOptions"] = {"MarketType": "spot", "SpotOptions": {"SpotInstanceType": "persistent", "InstanceInterruptionBehavior": "stop"}}
        try:
            r = cl("ec2", region).run_instances(**kw)["Instances"][0]
        except Exception as e:
            print(f"{name}: launch refused in {region}: {str(e)[:400]}"); return
        print(f"launched {name}: {r['InstanceId']} {itype} {'spot' if spot else 'on-demand'} in {region}")

def _ssm(xs, script, timeout=600, wait=True):
    cmds = {reg: cl("ssm", reg).send_command(InstanceIds=[x["InstanceId"] for x in xs if x["Region"] == reg], DocumentName="AWS-RunShellScript",
                                             Parameters={"commands": [script], "executionTimeout": [str(timeout)]})["Command"]["CommandId"]
            for reg in sorted({x["Region"] for x in xs})}
    if not wait: return cmds, {}
    res = {}
    for _ in range(timeout // 3 + 20):
        time.sleep(3)
        for x in xs:
            iid, c = x["InstanceId"], cl("ssm", x["Region"])
            if iid in res: continue
            try: inv = c.get_command_invocation(CommandId=cmds[x["Region"]], InstanceId=iid)
            except c.exceptions.InvocationDoesNotExist: continue
            if inv["Status"] not in ("Pending", "InProgress", "Delayed"):
                res[iid] = (inv["Status"], inv.get("StandardOutputContent", ""), inv.get("StandardErrorContent", ""))
        if len(res) == len(xs): break
    return cmds, res

def show():
    xs = fleet(); online = {}
    run = [x for x in xs if x["State"]["Name"] == "running"]
    for reg in sorted({x["Region"] for x in run}):
        ids = [x["InstanceId"] for x in run if x["Region"] == reg]
        for i in cl("ssm", reg).describe_instance_information(Filters=[{"Key": "InstanceIds", "Values": ids}])["InstanceInformationList"]:
            online[i["InstanceId"]] = i.get("PingStatus")
    up = [x for x in run if online.get(x["InstanceId"]) == "Online"]
    _, res = _ssm(up, "test -f /opt/data/READY && echo READY || echo booting; pgrep -fc '/opt/run/.*/task.sh' || true; uptime | sed 's/.*load/load/'", timeout=60) if up else (None, {})
    for x in xs:
        iid = x["InstanceId"]; r = res.get(iid); extra = " ".join(r[1].split()) if r else ""
        print(f"{_name(x):6s} {x['Region']:12s} {iid} {x['InstanceType']:14s} {x['State']['Name']:9s} ssm={online.get(iid, '-'):7s} {extra}")

def power(action, which):
    xs = fleet(); sel = xs if which == "all" else [x for x in xs if _name(x) == which]
    if not sel: sys.exit(f"no fleet instance {which}")
    if action == "terminate" and which == "all": sys.exit("terminate one instance at a time")
    for reg in sorted({x["Region"] for x in sel}):
        e = cl("ec2", reg); xr = [x for x in sel if x["Region"] == reg]
        try:
            {"start": e.start_instances, "stop": e.stop_instances, "terminate": e.terminate_instances}[action](InstanceIds=[x["InstanceId"] for x in xr])
            print(f"{action}: " + " ".join(_name(x) for x in xr) + f" ({reg})")
        except Exception as ex:
            print(f"{action} refused in {reg}: {str(ex)[:300]}")

def bundle():
    buf = io.BytesIO(); h = hashlib.sha256(); members = []
    for d in CODE_DIRS:
        p = os.path.join(ROOT, d)
        for f in sorted(os.listdir(p)):
            if f.endswith(".py"):
                b = open(os.path.join(p, f), "rb").read(); h.update(d.encode()); h.update(f.encode()); h.update(b); members.append((f"{d}/{f}", b))
    code = h.hexdigest()[:20]; key = f"bundle9/{code}.tgz"
    try: s3.head_object(Bucket=BUCKET, Key=key)
    except Exception:
        with tarfile.open(fileobj=buf, mode="w:gz") as tf:
            for name, b in members:
                ti = tarfile.TarInfo(name); ti.size = len(b); tf.addfile(ti, io.BytesIO(b))
        s3.put_object(Bucket=BUCKET, Key=key, Body=buf.getvalue()); print(f"uploaded code bundle {code} ({len(members)} files)")
    return code

TASK_SH = r"""#!/bin/bash
set -u
J=__JDIR__; T="$1"; W="$J/w/$T"; B=__BUCKET__; R=__REGION__; TH=__TH__; JOB=__JOB__
mkdir -p "$W/out"
for d in whest scripts; do mkdir -p "$W/$d"; cp "$J/code/$d/"*.py "$W/$d/" 2>/dev/null || true; done
cd "$W"
export PATH=/opt/venv/bin:$PATH OUT="$W/out" DATA=/opt/data/official PYTHONPATH="$W" PYTHONWARNINGS=ignore \
       OMP_NUM_THREADS=$TH OPENBLAS_NUM_THREADS=$TH MKL_NUM_THREADS=$TH
t0=$(date +%s)
bash "$J/tasks/$T.cmd" > "$W/so.txt" 2> "$W/se.txt"; rc=$?
{ cat "$W/so.txt"; if [ -s "$W/se.txt" ]; then echo; echo "[stderr]"; tail -c 6000 "$W/se.txt"; fi
  echo; echo "[exit $rc after $(( $(date +%s) - t0 ))s on $(hostname), $(nproc) cores, threads $TH]"; } > "$W/log_$T.txt"
for f in $(ls -A "$W/out" 2>/dev/null); do aws s3 cp "$W/out/$f" "s3://$B/results9/$JOB/$f" --quiet --region $R; done
aws s3 cp "$W/log_$T.txt" "s3://$B/results9/$JOB/log_$T.txt" --quiet --region $R
rm -rf "$W"
"""

def batch(job, path, threads=8, only=None, nslots=None, mem=6.0):
    tasks = [l.rstrip("\n").split("\t", 1) for l in open(path) if l.strip() and not l.startswith("#")]
    run = [x for x in fleet(("running",)) if only is None or _name(x) in only]
    if not run: sys.exit("no running fleet instance")
    code = bundle(); inst = {x["InstanceId"]: x for x in run}
    cores = {i: _vcpus(x["InstanceType"], x["Region"]) for i, x in inst.items()}
    gbs = {i: _mem_gb(x["InstanceType"], x["Region"]) for i, x in inst.items()}
    slots = {i: max(1, min(c // threads, int(0.94 * gbs[i] // mem))) if nslots is None else nslots for i, c in cores.items()}
    order = sorted(slots, key=lambda i: -slots[i]); share = {i: [] for i in order}; load = {i: 0.0 for i in order}
    for t in tasks:
        i = min(order, key=lambda k: (load[k] + 1) / slots[k]); share[i].append(t); load[i] += 1
    for iid, mine in share.items():
        if not mine: continue
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:gz") as tf:
            for tag, cmd in mine:
                b = (cmd + "\n").encode(); ti = tarfile.TarInfo(f"tasks/{tag}.cmd"); ti.size = len(b); tf.addfile(ti, io.BytesIO(b))
            sh = (TASK_SH.replace("__JDIR__", f"/opt/run/{job}").replace("__BUCKET__", BUCKET).replace("__REGION__", REGION)
                  .replace("__TH__", str(threads)).replace("__JOB__", job)).encode()
            ti = tarfile.TarInfo("task.sh"); ti.size = len(sh); ti.mode = 0o755; tf.addfile(ti, io.BytesIO(sh))
            tags = ("\n".join(t for t, _ in mine) + "\n").encode(); ti = tarfile.TarInfo("tags.txt"); ti.size = len(tags); tf.addfile(ti, io.BytesIO(tags))
        s3.put_object(Bucket=BUCKET, Key=f"jobs9/{job}/{iid}.tgz", Body=buf.getvalue())
        script = f"""set -u
for i in $(seq 90); do [ -f /opt/data/READY ] && break; sleep 10; done
[ -f /opt/data/READY ] || {{ echo "bootstrap not finished"; exit 1; }}
J=/opt/run/{job}; mkdir -p $J/code && cd $J
aws s3 cp s3://{BUCKET}/bundle9/{code}.tgz code.tgz --quiet --region {REGION} && tar xzf code.tgz -C code
aws s3 cp s3://{BUCKET}/jobs9/{job}/{iid}.tgz job.tgz --quiet --region {REGION} && tar xzf job.tgz
touch /opt/run/.last
nohup bash -c "xargs -P {slots[iid]} -I{{}} $J/task.sh {{}} < $J/tags.txt" > $J/xargs.log 2>&1 &
echo "{job}: {len(mine)} tasks, {slots[iid]} slots x {threads} threads on $(hostname) ($(nproc) cores)"
"""
        _ssm([inst[iid]], script, timeout=1200, wait=False)
        print(f"{job}: {len(mine)} tasks -> {_name(inst[iid])} {iid} {inst[iid]['Region']} ({slots[iid]} slots x {threads} threads)")
    d = os.path.join(SCRATCH, "results", job); os.makedirs(d, exist_ok=True)
    json.dump([t for t, _ in tasks], open(os.path.join(d, "tags.json"), "w"))

def _line(tag, log):
    body = log.split("\n[stderr]\n")[0]
    rc = next((l.split()[1] for l in log.splitlines() if l.startswith("[exit ")), "?")
    tail = [l for l in body.splitlines() if l.strip() and not l.startswith("[exit ")][-2:]
    if rc != "0": tail = [l for l in log.splitlines() if l.strip()][-4:]
    return f"[{tag} exit {rc}] " + " | ".join(tail)

def watch(job, timeout=7200, every=3.0):
    """Stream a job: every finished task's log AND output files are downloaded the moment they land in S3
    (tasks upload outputs before their log, so a log's arrival means its outputs are complete)."""
    d = os.path.join(SCRATCH, "results", job); os.makedirs(d, exist_ok=True)
    tags = json.load(open(os.path.join(d, "tags.json"))); seen = {}; have = set(os.listdir(d)); t0 = time.time()
    while len(seen) < len(tags) and time.time() - t0 < timeout:
        logs = []
        for page in s3.get_paginator("list_objects_v2").paginate(Bucket=BUCKET, Prefix=f"results9/{job}/"):
            for o in page.get("Contents", []):
                name = o["Key"].split("/", 2)[2]
                if name.startswith("log_"):
                    logs.append((name, o["Key"]))
                elif name not in have:
                    s3.download_file(BUCKET, o["Key"], os.path.join(d, name)); have.add(name)
        for name, key in logs:
            tag = name[4:-4]
            if tag in seen or tag not in tags: continue
            log = s3.get_object(Bucket=BUCKET, Key=key)["Body"].read().decode()
            open(os.path.join(d, name), "w").write(log)
            seen[tag] = _line(tag, log); print(seen[tag], f"(+{time.time() - t0:.0f}s)", flush=True)
        if len(seen) < len(tags): time.sleep(every)
    open(os.path.join(d, "summary.txt"), "w").write("\n".join(seen[t] for t in tags if t in seen) + "\n")
    print(f"{job}: {len(seen)}/{len(tags)} results in {time.time() - t0:.0f}s", flush=True)

def ready(which, timeout=1500):
    """Block until the named instance is running, its SSM agent is online and its bootstrap wrote /opt/data/READY."""
    t0 = time.time()
    while time.time() - t0 < timeout:
        xs = [x for x in fleet(("pending", "running")) if _name(x) == which]
        if xs and xs[0]["State"]["Name"] == "running":
            x = xs[0]
            info = cl("ssm", x["Region"]).describe_instance_information(Filters=[{"Key": "InstanceIds", "Values": [x["InstanceId"]]}])["InstanceInformationList"]
            if info and info[0].get("PingStatus") == "Online":
                _, res = _ssm([x], "test -f /opt/data/READY && echo READY $(ls /opt/data/official | wc -l) files $(nproc) cores || echo booting", timeout=60)
                out = " ".join(" ".join(v[1].split()) for v in res.values()) if res else ""
                if "READY" in out:
                    print(f"{which} ready after {time.time() - t0:.0f}s: {out}", flush=True); return True
        time.sleep(10)
    sys.exit(f"{which} not ready after {timeout}s")

def get(job, pattern="*"):
    d = os.path.join(SCRATCH, "results", job); os.makedirs(d, exist_ok=True); n = 0
    for page in s3.get_paginator("list_objects_v2").paginate(Bucket=BUCKET, Prefix=f"results9/{job}/"):
        for o in page.get("Contents", []):
            name = o["Key"].split("/", 2)[2]
            if fnmatch.fnmatch(name, pattern): s3.download_file(BUCKET, o["Key"], os.path.join(d, name)); n += 1
    print(f"{job}: {n} files -> {d}")

if __name__ == "__main__":
    a = sys.argv[1:]
    if not a: sys.exit(__doc__)
    c = a[0]
    if c == "launch": launch(int(a[1]), a[2], spot=("spot" in a[3:]), region=a[a.index("--region") + 1] if "--region" in a else None)
    elif c == "list": show()
    elif c == "ready": ready(a[1])
    elif c in ("start", "stop", "terminate"): power(c, a[1])
    elif c == "batch": batch(a[1], a[2], int(a[a.index("--threads") + 1]) if "--threads" in a else 8,
                             a[a.index("--only") + 1].split(",") if "--only" in a else None,
                             int(a[a.index("--slots") + 1]) if "--slots" in a else None, float(a[a.index("--mem") + 1]) if "--mem" in a else 6.0)
    elif c == "watch": watch(a[1])
    elif c == "get": get(a[1], a[2] if len(a) > 2 else "*")
    elif c == "ssm":
        xs = [x for x in fleet(("running",)) if _name(x) == a[1]]; _, res = _ssm(xs, a[2], timeout=300)
        for k, v in res.items(): print(k, v[0]); print(v[1]); print(v[2])
    else: sys.exit(__doc__)
