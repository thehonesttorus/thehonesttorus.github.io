#!/usr/bin/env python3
"""Drive infra/modal_app.py from the AWS relay instance (its egress carries gRPC; the sandbox's does not).
  python infra/modal_relay.py deploy [NAME]                 bundle whest/ scripts/ infra/modal_app.py -> S3 -> relay -> modal deploy
  python infra/modal_relay.py batch JOB FILE [--cpu N ...]  run a batch on Modal through the relay; results land in
                                                            s3://BUCKET/results9/JOB/ (watch with infra/fleet.py watch JOB)
  python infra/modal_relay.py mc JOB NET N SEED0 NSEEDS [--gpu L4] [--cov 1,7]
  python infra/modal_relay.py log JOB                       print the relay-side log of a job
  python infra/modal_relay.py fetch JOB                     copy a finished job's results from the Modal volume to S3
  python infra/modal_relay.py run ARGS...                   run "python infra/modal_app.py ARGS" on the relay (e.g. run ls official)
"""
import hashlib, io, json, os, sys, tarfile, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fleet
RELAY = os.environ.get("WHEST_RELAY", "w9-1"); BUCKET = fleet.BUCKET; REGION = fleet.REGION; ROOT = fleet.ROOT


def bundle(extra=()):
    buf = io.BytesIO(); h = hashlib.sha256()
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        for rel in [f"{d}/{f}" for d in ("whest", "scripts") for f in sorted(os.listdir(os.path.join(ROOT, d))) if f.endswith(".py")] + ["infra/modal_app.py"] + list(extra):
            b = open(os.path.join(ROOT, rel), "rb").read(); h.update(rel.encode()); h.update(b)
            ti = tarfile.TarInfo(rel); ti.size = len(b); tf.addfile(ti, io.BytesIO(b))
    code = h.hexdigest()[:20]; fleet.s3.put_object(Bucket=BUCKET, Key=f"bundle9/relay_{code}.tgz", Body=buf.getvalue()); return code


def relay(script, wait=True, timeout=1200):
    xs = [x for x in fleet.fleet(("running",)) if fleet._name(x) == RELAY]
    if not xs: sys.exit(f"relay {RELAY} is not running")
    _, res = fleet._ssm(xs, script, timeout=timeout, wait=wait)
    for k, v in res.items(): print(v[1][-6000:]); print(v[2][-2000:] if v[2].strip() else "", end="")


PRE = f"""set -u; export PATH=/opt/venv/bin:$PATH; mkdir -p /opt/modalrun && cd /opt/modalrun
aws s3 cp s3://{BUCKET}/bundle9/relay_CODE.tgz code_CODE.tgz --quiet --region {REGION} && rm -rf code_CODE && mkdir code_CODE && tar xzf code_CODE.tgz -C code_CODE && cd code_CODE
"""

if __name__ == "__main__":
    a = sys.argv[1:]
    if not a: sys.exit(__doc__)
    c = a[0]
    if c == "deploy":
        code = bundle(); relay(PRE.replace("CODE", code) + "modal deploy --strategy recreate infra/modal_app.py 2>&1 | tail -5", timeout=1500)
    elif c in ("batch", "mc"):
        job = a[1]
        if c == "batch":
            tasks = [l.rstrip("\n").split("\t", 1) for l in open(a[2]) if l.strip() and not l.startswith("#")]
            os.makedirs(os.path.join(ROOT, "jobs"), exist_ok=True); open(os.path.join(ROOT, f"jobs/{job}.tsv"), "w").write("\n".join("\t".join(t) for t in tasks) + "\n")
            d = os.path.join(fleet.SCRATCH, "results", job); os.makedirs(d, exist_ok=True); json.dump([t for t, _ in tasks], open(os.path.join(d, "tags.json"), "w"))
            code = bundle(extra=[f"jobs/{job}.tsv"]); opts = " ".join(a[3:])
            cmd = f"python infra/modal_app.py batch {job} jobs/{job}.tsv {opts}"
        else:
            code = bundle(); opts = " ".join(a[2:]); cmd = f"python infra/modal_app.py mc {job} {opts}"
        # the batch client keeps the relay alive (its idle watch looks at /opt/run/.last) and ships results to S3 at the end
        script = PRE.replace("CODE", code) + f"""nohup bash -c "(while true; do touch /opt/run/.last; sleep 60; done) & KA=\\$!; {cmd} > /opt/modalrun/{job}.log 2>&1; python infra/modal_app.py get {job} /opt/modalrun/out/{job} >> /opt/modalrun/{job}.log 2>&1; aws s3 cp /opt/modalrun/out/{job} s3://{BUCKET}/results9/{job}/ --recursive --quiet --region {REGION}; aws s3 cp /opt/modalrun/{job}.log s3://{BUCKET}/results9/{job}/relay_log.txt --quiet --region {REGION}; kill \\$KA" > /dev/null 2>&1 &
echo "relay: {job} started (code {code})"
"""
        relay(script, timeout=120)
    elif c == "run":                      # any client subcommand of infra/modal_app.py, synchronously (e.g. run "ls official")
        code = bundle(); relay(PRE.replace("CODE", code) + f"python infra/modal_app.py {' '.join(a[1:])} 2>&1 | tail -60", timeout=900)
    elif c == "fetch":                    # fetch a job's results from the Modal volume to S3 (e.g. after the relay was interrupted)
        job = a[1]; code = bundle()
        relay(PRE.replace("CODE", code) + f"python infra/modal_app.py get {job} /opt/modalrun/out/{job} 2>&1 | tail -2; aws s3 cp /opt/modalrun/out/{job} s3://{BUCKET}/results9/{job}/ --recursive --quiet --region {REGION}; ls /opt/modalrun/out/{job} | wc -l", timeout=900)
    elif c == "log":
        relay(f"tail -40 /opt/modalrun/{a[1]}.log 2>/dev/null || echo 'no log yet'", timeout=60)
    else:
        sys.exit(__doc__)
