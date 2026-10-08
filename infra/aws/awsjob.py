#!/usr/bin/env python3
"""AWS experiment runner (replaces infra/azure after the Azure subscription was lost).

Credentials come from the environment (AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY, AWS_DEFAULT_REGION); nothing here
prints them. Every resource it creates is tagged Project=claude-whest, and it refuses to stop, start or terminate an
instance without that tag. Jobs run through SSM Run Command (no SSH): the instance downloads the code bundle from S3,
runs the command in the background, and uploads out/<job>/ (log.txt, DONE with the exit code, anything else written
there) to s3://<bucket>/results/<job>/.

  python3 infra/aws/awsjob.py setup                         create the S3 bucket and the security group (idempotent)
  python3 infra/aws/awsjob.py launch NAME TYPE [GB] [spot]  launch an Ubuntu 24.04 instance with the Python stack
  python3 infra/aws/awsjob.py list                          our instances and their state
  python3 infra/aws/awsjob.py stop|start|terminate NAME
  python3 infra/aws/awsjob.py ready NAME                    has the bootstrap finished (SSM online and READY marker)
  python3 infra/aws/awsjob.py bundle DIR [DIR ...]          tar the given scratch dirs (code only) and upload as the bundle
  python3 infra/aws/awsjob.py submit JOB NAME 'CMD'         run CMD in /opt/work/src on instance NAME (returns at once)
  python3 infra/aws/awsjob.py results JOB [--wait SECONDS]  download results/<job>/ and print log.txt
"""
import io, json, os, secrets, sys, tarfile, time
import boto3

TAG = {"Key": "Project", "Value": "claude-whest"}
SCRATCH = os.environ.get("CLAUDE_AWS_DIR", "/tmp/claude-0/-home-user-thehonesttorus-github-io/"
                         "2cab30f0-c454-5769-8d9a-8195cd80a5f3/scratchpad/aws")
STATE = os.path.join(SCRATCH, "state.json")
ROLE = "claude-ec2-runner"
S = boto3.session.Session()
REGION = S.region_name
ec2, s3, ssm = S.client("ec2"), S.client("s3"), S.client("ssm")

BOOT = """#!/bin/bash
set -eux
export DEBIAN_FRONTEND=noninteractive
apt-get update -y
apt-get install -y python3-venv python3-pip unzip numactl htop
python3 -m venv /opt/venv
/opt/venv/bin/pip install -q --upgrade pip
/opt/venv/bin/pip install -q numpy scipy pyarrow flopscope==0.12.1 whestbench==0.16.1
curl -sSL https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip -o /tmp/awscli.zip
unzip -q /tmp/awscli.zip -d /tmp && /tmp/aws/install
mkdir -p /opt/work && chmod 777 /opt/work
/opt/venv/bin/python -c "import numpy; numpy.show_config()" > /opt/work/numpy_config.txt 2>&1
nproc > /opt/work/nproc.txt
touch /opt/work/READY
"""


def state():
    return json.load(open(STATE)) if os.path.exists(STATE) else {}


def save(st):
    os.makedirs(SCRATCH, exist_ok=True); json.dump(st, open(STATE, "w"), indent=1)


def ours(name=None, states=("pending", "running", "stopping", "stopped")):
    f = [{"Name": "tag:Project", "Values": [TAG["Value"]]}, {"Name": "instance-state-name", "Values": list(states)}]
    if name:
        f.append({"Name": "tag:Name", "Values": [name]})
    out = []
    for r in ec2.describe_instances(Filters=f)["Reservations"]:
        out += r["Instances"]
    return out


def one(name):
    xs = ours(name)
    if len(xs) != 1:
        sys.exit(f"expected exactly one tagged instance named {name}, found {len(xs)}")
    return xs[0]


def setup():
    st = state()
    if "bucket" not in st:
        acct = S.client("sts").get_caller_identity()["Account"]
        b = f"claude-whest-{acct[-4:]}-{secrets.token_hex(3)}"
        s3.create_bucket(Bucket=b, CreateBucketConfiguration={"LocationConstraint": REGION})
        s3.put_bucket_tagging(Bucket=b, Tagging={"TagSet": [TAG]})
        s3.put_public_access_block(Bucket=b, PublicAccessBlockConfiguration={k: True for k in (
            "BlockPublicAcls", "IgnorePublicAcls", "BlockPublicPolicy", "RestrictPublicBuckets")})
        st["bucket"] = b
    if "sg" not in st:
        vpc = ec2.describe_vpcs(Filters=[{"Name": "isDefault", "Values": ["true"]}])["Vpcs"][0]["VpcId"]
        sg = ec2.create_security_group(GroupName="claude-whest-sg", Description="claude-whest runner: no inbound",
                                       VpcId=vpc, TagSpecifications=[{"ResourceType": "security-group", "Tags": [TAG]}])
        st["sg"] = sg["GroupId"]
    save(st); print(f"bucket {st['bucket']} | security group {st['sg']} | region {REGION}")


def launch(name, itype, gb=256, spot=False):
    st = state()
    if ours(name):
        sys.exit(f"an instance named {name} already exists")
    ami = ssm.get_parameter(Name="/aws/service/canonical/ubuntu/server/24.04/stable/current/amd64/hvm/ebs-gp3/ami-id")["Parameter"]["Value"]
    kw = dict(ImageId=ami, InstanceType=itype, MinCount=1, MaxCount=1, UserData=BOOT,
              IamInstanceProfile={"Name": ROLE}, SecurityGroupIds=[st["sg"]],
              BlockDeviceMappings=[{"DeviceName": "/dev/sda1", "Ebs": {"VolumeSize": int(gb), "VolumeType": "gp3", "DeleteOnTermination": True}}],
              MetadataOptions={"HttpTokens": "required"},
              TagSpecifications=[{"ResourceType": t, "Tags": [TAG, {"Key": "Name", "Value": name}]} for t in ("instance", "volume")])
    if spot:
        kw["InstanceMarketOptions"] = {"MarketType": "spot", "SpotOptions": {"SpotInstanceType": "persistent",
                                                                             "InstanceInterruptionBehavior": "stop"}}
    r = ec2.run_instances(**kw)["Instances"][0]
    print(f"launched {name}: {r['InstanceId']} {itype} {'spot' if spot else 'on-demand'} (ami {ami})")


def show():
    for x in ours():
        nm = next((t["Value"] for t in x.get("Tags", []) if t["Key"] == "Name"), "?")
        print(f"{nm:10s} {x['InstanceId']} {x['InstanceType']:14s} {x['State']['Name']:9s} {x.get('LaunchTime')}")


def power(action, name):
    iid = one(name)["InstanceId"]
    {"stop": ec2.stop_instances, "start": ec2.start_instances, "terminate": ec2.terminate_instances}[action](InstanceIds=[iid])
    print(f"{action} {name} ({iid})")


def run_ssm(iid, script, timeout=600):
    c = ssm.send_command(InstanceIds=[iid], DocumentName="AWS-RunShellScript",
                         Parameters={"commands": [script], "executionTimeout": [str(timeout)]})["Command"]["CommandId"]
    for _ in range(timeout // 3 + 20):
        time.sleep(3)
        try:
            inv = ssm.get_command_invocation(CommandId=c, InstanceId=iid)
        except ssm.exceptions.InvocationDoesNotExist:
            continue
        if inv["Status"] not in ("Pending", "InProgress", "Delayed"):
            return inv["Status"], inv.get("StandardOutputContent", ""), inv.get("StandardErrorContent", "")
    return "Timeout", "", ""


def ready(name):
    iid = one(name)["InstanceId"]
    info = ssm.describe_instance_information(Filters=[{"Key": "InstanceIds", "Values": [iid]}])["InstanceInformationList"]
    if not info or info[0].get("PingStatus") != "Online":
        print(f"{name}: SSM agent not online yet"); return False
    stt, out, err = run_ssm(iid, "test -f /opt/work/READY && echo READY $(cat /opt/work/nproc.txt) cores || "
                                 "(echo NOT-READY; tail -3 /var/log/cloud-init-output.log)", timeout=60)
    print(f"{name}: {out.strip()} {err.strip()[:200]}"); return "READY" in out


def bundle(dirs):
    st = state(); buf = io.BytesIO()
    root = os.path.dirname(SCRATCH)
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        for d in dirs:
            p = os.path.join(root, d)
            for f in sorted(os.listdir(p)):
                if f.endswith((".py", ".sh")):
                    tf.add(os.path.join(p, f), arcname=f"{d}/{f}")
    buf.seek(0); s3.put_object(Bucket=st["bucket"], Key="bundle/codebundle.tgz", Body=buf.getvalue())
    print(f"uploaded bundle ({len(buf.getvalue()) // 1024} KB) from {dirs}")


def submit(job, name, cmd):
    st = state(); iid = one(name)["InstanceId"]; b = st["bucket"]
    script = f"""set -u
for i in $(seq 120); do [ -f /opt/work/READY ] && break; sleep 10; done
[ -f /opt/work/READY ] || {{ echo "bootstrap not finished after 20 min"; exit 1; }}
mkdir -p /opt/work/src && cd /opt/work/src
aws s3 cp s3://{b}/bundle/codebundle.tgz /opt/work/codebundle.tgz --quiet --region {REGION} && tar xzf /opt/work/codebundle.tgz -C /opt/work/src
mkdir -p out/{job}
cat > /opt/work/cmd_{job}.sh <<'CMDEOF'
{cmd}
CMDEOF
cat > /opt/work/run_{job}.sh <<RUNEOF
#!/bin/bash
cd /opt/work/src
export PATH=/opt/venv/bin:\\$PATH OUT=out/{job}
bash /opt/work/cmd_{job}.sh > out/{job}/log.txt 2>&1
echo \\$? > out/{job}/DONE
aws s3 cp out/{job} s3://{b}/results/{job}/ --recursive --quiet --region {REGION}
RUNEOF
chmod +x /opt/work/run_{job}.sh
nohup /opt/work/run_{job}.sh >/dev/null 2>&1 &
echo "launched {job} on $(hostname), $(nproc) cores"
"""
    stt, out, err = run_ssm(iid, script, timeout=1500)
    print(f"{job} {name} {stt}: {out.strip()[-300:]} {err.strip()[-300:]}")


def results(job, wait=0):
    st = state(); b = st["bucket"]; t0 = time.time(); pre = f"results/{job}/"
    while True:
        keys = [o["Key"] for o in s3.list_objects_v2(Bucket=b, Prefix=pre).get("Contents", [])]
        if any(k.endswith("/DONE") for k in keys) or time.time() - t0 >= wait:
            break
        time.sleep(20)
    d = os.path.join(SCRATCH, "results", job); os.makedirs(d, exist_ok=True)
    for k in keys:
        p = os.path.join(d, k[len(pre):]); os.makedirs(os.path.dirname(p), exist_ok=True)
        s3.download_file(b, k, p)
    done = os.path.exists(os.path.join(d, "DONE"))
    print(f"job {job}: {'DONE exit ' + open(os.path.join(d, 'DONE')).read().strip() if done else 'not finished'} | {len(keys)} files in {d}")
    lp = os.path.join(d, "log.txt")
    if os.path.exists(lp):
        print(open(lp).read()[-4000:])


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a:
        sys.exit(__doc__)
    c = a[0]
    if c == "setup": setup()
    elif c == "launch": launch(a[1], a[2], *(a[3:4] or [256]), spot=("spot" in a[4:]))
    elif c == "list": show()
    elif c in ("stop", "start", "terminate"): power(c, a[1])
    elif c == "ready": ready(a[1])
    elif c == "bundle": bundle(a[1:])
    elif c == "submit": submit(a[1], a[2], a[3])
    elif c == "results": results(a[1], int(a[a.index("--wait") + 1]) if "--wait" in a else 0)
    else: sys.exit(__doc__)
