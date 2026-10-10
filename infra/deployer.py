#!/usr/bin/env python3
"""Deploy infra/modal_gw.py from a small AWS instance: the only step that needs Modal's gRPC client.

  python infra/deployer.py deploy      start (or first launch) the t3.small "d10", run `modal deploy`, record the
                                       gateway URL and token in ~/.whest10/gw.json; the instance stops itself 30 min later
  python infra/deployer.py status | stop | terminate

Secrets never appear in a command line or log: the Modal token (from ~/.modal.toml, profile research-4885) and the
gateway token travel as one S3 object that the instance reads into its environment and deletes at once (and the
sandbox deletes it again afterwards). The instance has its own Fleet tag, so infra/fleet.py never schedules work on it.
"""
import io, json, os, re, secrets, sys, tarfile, time, tomllib, uuid
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fleet

NAME, ITYPE, PROFILE = "d10", "t3.small", os.environ.get("MODAL_PROFILE", "research-4885")
FLEET_TAG = {"Key": "Fleet", "Value": "whest10-deploy"}
REGION, BUCKET = fleet.REGION, fleet.BUCKET
HOME = os.path.expanduser("~/.whest10"); CFG = os.path.join(HOME, "gw.json")
ec2, ssm, s3 = fleet.cl("ec2"), fleet.cl("ssm"), fleet.s3

BOOT = r"""#!/bin/bash
set -eux
export DEBIAN_FRONTEND=noninteractive
apt-get update -y && apt-get install -y python3-venv unzip
python3 -m venv /opt/mvenv && /opt/mvenv/bin/pip install -q --upgrade pip && /opt/mvenv/bin/pip install -q modal==1.6.1 boto3
curl -sSL https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip -o /tmp/awscli.zip && unzip -q /tmp/awscli.zip -d /tmp && /tmp/aws/install
mkdir -p /opt/run && touch /opt/run/.last
cat > /usr/local/bin/idlestop.sh <<'EOW'
#!/bin/bash
while true; do sleep 60; if [ $(( $(date +%s) - $(stat -c %Y /opt/run/.last) )) -gt 1800 ]; then shutdown -h now; fi; done
EOW
chmod +x /usr/local/bin/idlestop.sh
cat > /etc/systemd/system/idlestop.service <<'EOS'
[Unit]
Description=stop the deploy instance when idle
[Service]
ExecStartPre=/bin/touch /opt/run/.last
ExecStart=/usr/local/bin/idlestop.sh
Restart=always
[Install]
WantedBy=multi-user.target
EOS
systemctl daemon-reload && systemctl enable --now idlestop.service
touch /opt/run/READY
"""


def instance():
    f = [{"Name": "tag:Project", "Values": [fleet.TAG["Value"]]}, {"Name": "tag:Fleet", "Values": [FLEET_TAG["Value"]]},
         {"Name": "tag:Name", "Values": [NAME]}, {"Name": "instance-state-name", "Values": ["pending", "running", "stopping", "stopped"]}]
    xs = [x for r in ec2.describe_instances(Filters=f)["Reservations"] for x in r["Instances"]]
    return xs[0] if xs else None


def up():
    x = instance()
    if x is None:
        ami = fleet.cl("ssm").get_parameter(Name=fleet.AMI_PARAM)["Parameter"]["Value"]
        x = ec2.run_instances(ImageId=ami, InstanceType=ITYPE, MinCount=1, MaxCount=1, UserData=BOOT,
                              IamInstanceProfile={"Name": fleet.ROLE}, SecurityGroupIds=[fleet._sg(REGION)],
                              InstanceInitiatedShutdownBehavior="stop", MetadataOptions={"HttpTokens": "required"},
                              BlockDeviceMappings=[{"DeviceName": "/dev/sda1", "Ebs": {"VolumeSize": 16, "VolumeType": "gp3", "DeleteOnTermination": True}}],
                              TagSpecifications=[{"ResourceType": t, "Tags": [fleet.TAG, FLEET_TAG, {"Key": "Name", "Value": NAME}]} for t in ("instance", "volume")])["Instances"][0]
        print(f"launched {NAME} {x['InstanceId']} ({ITYPE})", flush=True)
    iid = x["InstanceId"]
    st = x["State"]["Name"]
    if st == "stopping":
        ec2.get_waiter("instance_stopped").wait(InstanceIds=[iid]); st = "stopped"
    if st == "stopped":
        ec2.start_instances(InstanceIds=[iid]); print(f"starting {NAME}", flush=True)
    ec2.get_waiter("instance_running").wait(InstanceIds=[iid], WaiterConfig={"Delay": 5, "MaxAttempts": 120})
    for _ in range(120):                       # SSM agent online
        info = ssm.describe_instance_information(Filters=[{"Key": "InstanceIds", "Values": [iid]}])["InstanceInformationList"]
        if info and info[0].get("PingStatus") == "Online":
            break
        time.sleep(5)
    for _ in range(80):                        # boot script finished (first launch only takes a while)
        if "READY" in run(iid, "test -f /opt/run/READY && echo READY || echo booting", 60)[1]:
            return iid
        time.sleep(10)
    raise RuntimeError("deploy instance did not become ready")


def run(iid, script, timeout):
    cid = ssm.send_command(InstanceIds=[iid], DocumentName="AWS-RunShellScript",
                           Parameters={"commands": [script], "executionTimeout": [str(timeout)]})["Command"]["CommandId"]
    for _ in range(timeout // 2 + 30):
        time.sleep(2)
        try:
            inv = ssm.get_command_invocation(CommandId=cid, InstanceId=iid)
        except ssm.exceptions.InvocationDoesNotExist:
            continue
        if inv["Status"] not in ("Pending", "InProgress", "Delayed"):
            return inv["Status"], inv.get("StandardOutputContent", ""), inv.get("StandardErrorContent", "")
    return "Timeout", "", ""


def deploy():
    os.makedirs(HOME, exist_ok=True); os.chmod(HOME, 0o700)
    cfg = json.load(open(CFG)) if os.path.exists(CFG) else {}
    gw_token = cfg.get("token") or secrets.token_urlsafe(32)
    prof = tomllib.load(open(os.path.expanduser("~/.modal.toml"), "rb"))[PROFILE]
    iid = up()
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        tf.add(os.path.join(fleet.ROOT, "infra/modal_gw.py"), arcname="infra/modal_gw.py")
    code_key = f"bundle10/deploy_{uuid.uuid4().hex[:12]}.tgz"; s3.put_object(Bucket=BUCKET, Key=code_key, Body=buf.getvalue())
    sec_key = f"tmp10/{uuid.uuid4().hex}.json"
    s3.put_object(Bucket=BUCKET, Key=sec_key, ServerSideEncryption="AES256",
                  Body=json.dumps({"MODAL_TOKEN_ID": prof["token_id"], "MODAL_TOKEN_SECRET": prof["token_secret"], "GW_TOKEN": gw_token}).encode())
    script = f"""set -eu
touch /opt/run/.last; cd /opt/run && rm -rf dep && mkdir dep && cd dep
aws s3 cp s3://{BUCKET}/{code_key} c.tgz --quiet --region {REGION} && tar xzf c.tgz
eval "$(/opt/mvenv/bin/python - <<'PY'
import json, shlex, boto3
s3 = boto3.client("s3", region_name="{REGION}")
o = json.loads(s3.get_object(Bucket="{BUCKET}", Key="{sec_key}")["Body"].read()); s3.delete_object(Bucket="{BUCKET}", Key="{sec_key}")
print(" ".join(f"export {{k}}={{shlex.quote(v)}};" for k, v in o.items()))
PY
)"
/opt/mvenv/bin/modal deploy infra/modal_gw.py 2>&1 | tail -40
touch /opt/run/.last
"""
    try:
        t0 = time.time(); st, out, err = run(iid, script, 1500)
    finally:
        s3.delete_object(Bucket=BUCKET, Key=sec_key); s3.delete_object(Bucket=BUCKET, Key=code_key)
    print(out[-4000:], err[-2000:], sep="\n")
    m = re.search(r"https://[A-Za-z0-9.-]+--whest10-gw\.modal\.run", out)
    if st != "Success" or not m:
        sys.exit(f"deploy failed ({st}) after {time.time() - t0:.0f}s")
    fd = os.open(CFG, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    os.write(fd, json.dumps({"url": m.group(0), "token": gw_token}).encode()); os.close(fd)
    print(f"gateway {m.group(0)} deployed in {time.time() - t0:.0f}s; config in {CFG}")


if __name__ == "__main__":
    c = sys.argv[1] if len(sys.argv) > 1 else ""
    if c == "deploy":
        deploy()
    elif c == "status":
        x = instance(); print(f"{NAME}: {x['InstanceId']} {x['InstanceType']} {x['State']['Name']}" if x else f"{NAME}: none")
    elif c in ("stop", "terminate"):
        x = instance()
        if x:
            (ec2.stop_instances if c == "stop" else ec2.terminate_instances)(InstanceIds=[x["InstanceId"]]); print(f"{c} {NAME}")
    else:
        sys.exit(__doc__)
