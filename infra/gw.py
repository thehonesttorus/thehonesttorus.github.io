#!/usr/bin/env python3
"""Client for the whest10 Modal gateway (infra/modal_gw.py), over plain HTTPS.

  python infra/gw.py ping
  python infra/gw.py run JOB TASKS [--cpu 8] [--mem 16384] [--timeout 3600] [--env K=V ...]
        TASKS: lines "tag<TAB>command" (or a .json list of task dicts). Uploads the code tree if it changed, spawns every
        task at once, then streams: each finished task is written to RESULTS/JOB/<tag>/ the moment it is reported, and one
        line per task is appended to RESULTS/JOB/stream.txt (tag, exit, seconds, last stdout line).
  python infra/gw.py submit JOB TASKS [...]        the same without watching
  python infra/gw.py watch JOB                     resume streaming a submitted job
  python infra/gw.py cancel JOB
  python infra/gw.py seed                          copy the official networks from S3 into the Modal data volume

From Python: `for tag, res in gw.map(job, tasks): ...` yields results in completion order.
Config ~/.whest10/gw.json {"url", "token"} is written by infra/deployer.py.
"""
import base64, hashlib, io, json, os, sys, tarfile, time
import requests

HOME = os.path.expanduser("~/.whest10"); CFG = os.path.join(HOME, "gw.json")
RESULTS = os.environ.get("GW_RESULTS", os.path.join(HOME, "results"))
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
CODE_DIRS = ("whest", "scripts", "systems", "infra")
CODE_EXT = (".py", ".txt", ".json", ".cfg", ".tsv")
_S = None


def _cfg():
    return json.load(open(CFG))


def _sess():
    global _S
    if _S is None:
        _S = requests.Session()
        _S.headers["Authorization"] = f"Bearer {_cfg()['token']}"
        _S.verify = os.environ.get("REQUESTS_CA_BUNDLE") or (os.path.exists("/root/.ccr/ca-bundle.crt") and "/root/.ccr/ca-bundle.crt") or True
    return _S


def _req(method, path, tries=6, **kw):
    """HTTPS with backoff on transport errors and 5xx; 4xx are final (401 = bad token, 404 = unknown job)."""
    url = _cfg()["url"].rstrip("/") + path; tmo = kw.pop("timeout", 120)
    for k in range(tries):
        try:
            r = _sess().request(method, url, timeout=tmo, **kw)
            if r.status_code < 500:
                if r.status_code >= 400:
                    raise RuntimeError(f"{method} {path}: HTTP {r.status_code} {r.text[:500]}")
                return r
            err = f"HTTP {r.status_code} {r.text[:300]}"
        except (requests.ConnectionError, requests.Timeout) as e:
            err = f"{type(e).__name__}: {str(e)[:300]}"
        if k == tries - 1:
            raise RuntimeError(f"{method} {path}: {err}")
        time.sleep(min(30, 2 ** k))


def bundle():
    """Deterministic tar.gz of the code tree; the hash names the snapshot so unchanged code is never re-uploaded."""
    files = []
    for d in CODE_DIRS:
        for base, dirs, fs in os.walk(os.path.join(ROOT, d)):
            dirs[:] = sorted(x for x in dirs if not x.startswith((".", "__")))
            for f in sorted(fs):
                p = os.path.join(base, f)
                if f.endswith(CODE_EXT) and os.path.getsize(p) < 8 << 20:
                    files.append(os.path.relpath(p, ROOT))
    buf = io.BytesIO(); h = hashlib.sha256()
    with tarfile.open(fileobj=buf, mode="w:gz", compresslevel=6) as tf:
        for rel in files:
            b = open(os.path.join(ROOT, rel), "rb").read(); h.update(rel.encode() + b"\0" + b)
            ti = tarfile.TarInfo(rel); ti.size = len(b); ti.mtime = 0; ti.mode = 0o644
            tf.addfile(ti, io.BytesIO(b))
    return h.hexdigest()[:24], buf.getvalue()


def push():
    sha, blob = bundle()
    if not _req("GET", f"/code/{sha}").json()["have"]:
        _req("POST", f"/code/{sha}", data=blob, headers={"Content-Type": "application/gzip"}, timeout=300)
        print(f"code {sha}: uploaded {len(blob) / 1e6:.2f} MB", file=sys.stderr)
    return sha


def submit(job, tasks, cpu=8, mem=16384, timeout=3600, env=None, code=None):
    code = code or push()
    tasks = [t if isinstance(t, dict) else {"tag": t[0], "cmd": t[1]} for t in tasks]
    tags = [t["tag"] for t in tasks]
    if len(set(tags)) != len(tags):
        raise ValueError("duplicate task tags")
    out = []
    for i in range(0, len(tasks), 400):          # keep each request small; all chunks are spawned within seconds
        out.append(_req("POST", "/submit", json={"job": job, "code": code, "tasks": tasks[i:i + 400], "cpu": cpu,
                                                 "mem": mem, "timeout": timeout, "env": env or {}}, timeout=300).json())
    return {"job": job, "code": code, "spawned": sum(o["spawned"] for o in out)}


def _save(job, tag, res, failed=False):
    d = os.path.join(RESULTS, job, tag); os.makedirs(d, exist_ok=True)
    if failed:
        open(os.path.join(d, "error.txt"), "w").write(res["error"]); line = f"{tag}\tFAILED\t-\t{res['error'][:200]}"
    else:
        open(os.path.join(d, "stdout.txt"), "w").write(res["stdout"])
        if res["stderr"].strip():
            open(os.path.join(d, "stderr.txt"), "w").write(res["stderr"])
        blob = base64.b64decode(res["files"]) if res.get("files") else b""
        if res.get("stored"):
            blob = _req("GET", "/file", params={"path": res["stored"]}, timeout=600).content
        if blob:
            with tarfile.open(fileobj=io.BytesIO(blob)) as tf:
                tf.extractall(d, filter="data")
        json.dump({k: res[k] for k in ("exit", "secs", "cpus", "stored")}, open(os.path.join(d, "meta.json"), "w"))
        last = next((l for l in reversed(res["stdout"].splitlines()) if l.strip()), "")
        if res["exit"] != 0:
            last = (last + " | " + res["stderr"].strip().splitlines()[-1] if res["stderr"].strip() else last)[:400]
        line = f"{tag}\t{res['exit']}\t{res['secs']}\t{last}"
    with open(os.path.join(RESULTS, job, "stream.txt"), "a") as f:
        f.write(line + "\n")
    return line


def map(job, tasks=None, every=2.0, quiet=False, **kw):
    """Submit (if tasks are given) and yield (tag, result) in completion order; a failed task yields {"error": ...}."""
    if tasks is not None:
        submit(job, tasks, **kw)
    have = set(); total = None; t0 = time.time(); last_note = ""
    sp = os.path.join(RESULTS, job, "stream.txt")
    if os.path.exists(sp):                       # resuming: tasks already saved locally are not fetched again
        have = {l.split("\t", 1)[0] for l in open(sp) if l.strip()}
    while total is None or len(have) < total:
        r = _req("POST", "/poll", json={"job": job, "have": sorted(have)}, timeout=180).json()
        total = r["total"]
        for tag, res in list(r["done"].items()) + list(r["failed"].items()):
            failed = tag in r["failed"]
            line = _save(job, tag, res, failed=failed); have.add(tag)
            if not quiet:
                print(line, flush=True)
            yield tag, res
        note = f"[{job}] {len(have)}/{total} done, {r['running']} running, {r['pending'] - r['running']} queued, {time.time() - t0:.0f}s"
        if not quiet and note.split(",")[0] != last_note.split(",")[0]:
            print(note, file=sys.stderr, flush=True); last_note = note
        if len(have) < total:
            time.sleep(every)


def _tasks(path):
    if path.endswith(".json"):
        return json.load(open(path))
    return [l.rstrip("\n").split("\t", 1) for l in open(path) if l.strip() and not l.startswith("#")]


def seed():
    """Copy s3://BUCKET/data/official/* into the Modal volume via presigned URLs (the data never transits the sandbox)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import fleet
    keys = []
    for page in fleet.s3.get_paginator("list_objects_v2").paginate(Bucket=fleet.BUCKET, Prefix="data/official/"):
        keys += [o["Key"] for o in page.get("Contents", [])]
    files = {k.rsplit("/", 1)[1]: fleet.s3.generate_presigned_url("get_object", Params={"Bucket": fleet.BUCKET, "Key": k}, ExpiresIn=6 * 3600)
             for k in keys if not k.endswith("/")}
    cid = _req("POST", "/seed", json={"files": files}).json()["call"]
    print(f"seeding {len(files)} files (call {cid})", flush=True)
    while True:
        r = _req("GET", f"/call/{cid}").json()
        if r["state"] != "pending":
            print(r); return r
        time.sleep(10)


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a:
        sys.exit(__doc__)
    opt = lambda k, d: type(d)(a[a.index(k) + 1]) if k in a else d
    env = dict(a[i + 1].split("=", 1) for i, x in enumerate(a) if x == "--env")
    kw = dict(cpu=opt("--cpu", 8.0), mem=opt("--mem", 16384), timeout=opt("--timeout", 3600), env=env)
    c = a[0]
    if c == "ping":
        print(_req("GET", "/ping").json())
    elif c == "push":
        print(push())
    elif c == "submit":
        print(submit(a[1], _tasks(a[2]), **kw))
    elif c == "run":
        for _ in map(a[1], _tasks(a[2]), **kw):
            pass
    elif c == "watch":
        for _ in map(a[1]):
            pass
    elif c == "cancel":
        print(_req("POST", "/cancel", json={"job": a[1]}).json())
    elif c == "seed":
        seed()
    else:
        sys.exit(__doc__)
