"""Modal app "whest10": an HTTPS gateway in front of a pool of runners.

The sandbox's proxy cannot carry gRPC, so the Modal client cannot run there. Only `modal deploy` of this file needs the
client, and infra/deployer.py does that once from a small AWS instance. Everything else, uploading code, fanning out
tasks and collecting results, is plain HTTPS to the gateway (client: infra/gw.py).

  POST /code/{sha}   body = tar.gz of the code tree; stored once on the volume whest10-store
  GET  /code/{sha}   {"have": bool}
  POST /submit       {"job", "code", "tasks": [{"tag", "cmd", "cpu"?, "mem"?, "timeout"?, "env"?}], "cpu", "mem", "timeout", "env"}
                     every task is spawned at once; the call returns immediately
  POST /poll         {"job", "have": [tags already received]} -> every task finished since, with its log and output files
  POST /cancel       {"job"}
  GET  /ping         data and runner sanity

Each task runs `bash -c cmd` in a fresh copy of the code tree with DATA=/data/official (volume whest-data) and
OUT=./out; whatever it writes to $OUT comes back with its log. Requests carry `Authorization: Bearer <GW_TOKEN>`.
"""
import asyncio, base64, hmac, io, os, shutil, subprocess, tarfile, time
import modal

APP = "whest10"
app = modal.App(APP)
data_vol = modal.Volume.from_name("whest-data", create_if_missing=True)
store = modal.Volume.from_name("whest10-store", create_if_missing=True)
jobs = modal.Dict.from_name("whest10-jobs", create_if_missing=True)
gw_secret = modal.Secret.from_dict({"GW_TOKEN": os.environ.get("GW_TOKEN", "")})

run_image = (modal.Image.debian_slim(python_version="3.12")
             .pip_install("numpy", "scipy", "pyarrow", "flopscope==0.12.1", "whestbench==0.16.1"))
gw_image = modal.Image.debian_slim(python_version="3.12").pip_install("fastapi[standard]")
INLINE_MAX = 48 << 20          # outputs above this stay on the store volume and are fetched with GET /file


@app.function(image=run_image, volumes={"/data": data_vol, "/store": store}, cpu=8.0, memory=16384, timeout=4 * 3600)
def runner(job: str, tag: str, cmd: str, code: str, env: dict, budget: int):
    t0 = time.time()
    try:
        jobs.put(f"start/{job}/{tag}", t0)
    except Exception:
        pass
    src = f"/store/code/{code}.tgz"
    if not os.path.exists(src):
        store.reload()
    work = f"/tmp/w_{tag}_{int(t0 * 1e3) % 10**9}"
    os.makedirs(work); tarfile.open(src).extractall(work, filter="data")
    os.makedirs(f"{work}/out", exist_ok=True)
    th = str((env or {}).get("NTHREADS", max(1, os.cpu_count() or 1)))
    e = dict(os.environ, DATA="/data/official", OUT=f"{work}/out", PYTHONPATH=work, PYTHONWARNINGS="ignore",
             OMP_NUM_THREADS=th, OPENBLAS_NUM_THREADS=th, MKL_NUM_THREADS=th)
    e.update({k: str(v) for k, v in (env or {}).items()})
    try:
        p = subprocess.run(["bash", "-c", cmd], cwd=work, env=e, capture_output=True, text=True, timeout=max(60, budget - 60))
        rc, so, se = p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired as x:
        rc, so, se = -9, (x.stdout or b"").decode() if isinstance(x.stdout, bytes) else (x.stdout or ""), f"[timeout after {budget - 60}s]"
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        tf.add(f"{work}/out", arcname=".")
    blob = buf.getvalue(); stored = None
    if len(blob) > INLINE_MAX:
        stored = f"out/{job}/{tag}.tgz"; os.makedirs(os.path.dirname(f"/store/{stored}"), exist_ok=True)
        open(f"/store/{stored}", "wb").write(blob); store.commit(); blob = b""
    shutil.rmtree(work, ignore_errors=True)
    return {"exit": rc, "stdout": so[-400_000:], "stderr": se[-40_000:], "secs": round(time.time() - t0, 2),
            "cpus": os.cpu_count(), "files": blob, "stored": stored}


@app.function(image=gw_image, secrets=[gw_secret], volumes={"/store": store, "/data": data_vol}, timeout=900,
              scaledown_window=1200)
@modal.concurrent(max_inputs=64)
@modal.asgi_app(label="whest10-gw")
def gateway():
    from fastapi import FastAPI, HTTPException, Request
    from fastapi.responses import Response
    web = FastAPI()
    token = os.environ.get("GW_TOKEN", "")
    variants = {}

    def auth(req):
        got = req.headers.get("authorization", "")
        if not token or not hmac.compare_digest(got, f"Bearer {token}"):
            raise HTTPException(status_code=401)

    def variant(cpu, mem, timeout):
        key = (float(cpu), int(mem), int(timeout))
        if key not in variants:
            variants[key] = runner.with_options(cpu=key[0], memory=key[1], timeout=key[2])
        return variants[key]

    @web.get("/ping")
    async def ping(req: Request):
        auth(req)
        off = sorted(os.listdir("/data/official")) if os.path.isdir("/data/official") else []
        return {"ok": True, "official_files": len(off), "first": off[:3], "store_code": len(os.listdir("/store/code")) if os.path.isdir("/store/code") else 0}

    @web.get("/code/{sha}")
    async def have_code(sha: str, req: Request):
        auth(req)
        if not os.path.exists(f"/store/code/{sha}.tgz"):
            await store.reload.aio()
        return {"have": os.path.exists(f"/store/code/{sha}.tgz")}

    @web.post("/code/{sha}")
    async def put_code(sha: str, req: Request):
        auth(req)
        body = await req.body()
        os.makedirs("/store/code", exist_ok=True)
        tmp = f"/store/code/.{sha}.part"; open(tmp, "wb").write(body); os.replace(tmp, f"/store/code/{sha}.tgz")
        await store.commit.aio()
        return {"stored": sha, "bytes": len(body)}

    @web.post("/submit")
    async def submit(req: Request):
        auth(req)
        b = await req.json()
        job, code, tasks = b["job"], b["code"], b["tasks"]
        if not os.path.exists(f"/store/code/{code}.tgz"):
            await store.reload.aio()
            if not os.path.exists(f"/store/code/{code}.tgz"):
                raise HTTPException(status_code=400, detail=f"code {code} not uploaded")
        sem = asyncio.Semaphore(128)

        async def one(t):
            cpu, mem, tmo = t.get("cpu", b.get("cpu", 8)), t.get("mem", b.get("mem", 16384)), t.get("timeout", b.get("timeout", 3600))
            async with sem:
                env = {"NTHREADS": str(max(1, int(round(float(cpu))))), **b.get("env", {}), **t.get("env", {})}
                fc = await variant(cpu, mem, tmo).spawn.aio(job, t["tag"], t["cmd"], code, env, int(tmo))
            return t["tag"], fc.object_id
        pairs = await asyncio.gather(*[one(t) for t in tasks])
        rec = await jobs.get.aio(f"job/{job}", None) or {"t": time.time(), "tasks": {}}
        rec["tasks"].update(dict(pairs)); rec["code"] = code
        await jobs.put.aio(f"job/{job}", rec)
        return {"job": job, "spawned": len(pairs), "total": len(rec["tasks"])}

    @web.post("/poll")
    async def poll(req: Request):
        auth(req)
        b = await req.json()
        rec = await jobs.get.aio(f"job/{b['job']}", None)
        if rec is None:
            raise HTTPException(status_code=404, detail="unknown job")
        have = set(b.get("have", [])); todo = [(t, c) for t, c in rec["tasks"].items() if t not in have]
        sem = asyncio.Semaphore(96)

        async def chk(tag, cid):
            async with sem:
                try:
                    r = await modal.FunctionCall.from_id(cid).get.aio(timeout=0)
                    r = dict(r); r["files"] = base64.b64encode(r["files"]).decode()
                    return tag, "done", r
                except TimeoutError:
                    return tag, "pending", None
                except Exception as e:
                    return tag, "failed", {"error": f"{type(e).__name__}: {str(e)[:3000]}"}
        res = await asyncio.gather(*[chk(t, c) for t, c in todo])
        pend = [t for t, s, _ in res if s == "pending"]
        started = 0
        if pend:
            started = sum([1 for x in await asyncio.gather(*[jobs.contains.aio(f"start/{b['job']}/{t}") for t in pend]) if x])
        return {"done": {t: r for t, s, r in res if s == "done"}, "failed": {t: r for t, s, r in res if s == "failed"},
                "pending": len(pend), "running": started, "total": len(rec["tasks"])}

    @web.get("/file")
    async def get_file(path: str, req: Request):
        auth(req)
        await store.reload.aio()
        p = os.path.normpath(f"/store/{path}")
        if not p.startswith("/store/out/") or not os.path.exists(p):
            raise HTTPException(status_code=404)
        return Response(open(p, "rb").read(), media_type="application/gzip")

    @web.post("/cancel")
    async def cancel(req: Request):
        auth(req)
        b = await req.json(); rec = await jobs.get.aio(f"job/{b['job']}", None) or {"tasks": {}}
        n = 0
        for t, c in rec["tasks"].items():
            if t in set(b.get("tags", rec["tasks"])):
                try:
                    await modal.FunctionCall.from_id(c).cancel.aio(); n += 1
                except Exception:
                    pass
        return {"cancelled": n}

    @web.post("/seed")
    async def seed(req: Request):
        """Fill /data/official from presigned HTTPS URLs ({"files": {name: url}}), in parallel; returns a call id."""
        auth(req)
        b = await req.json()
        fc = await fetch_data.spawn.aio(b["files"])
        return {"call": fc.object_id}

    @web.get("/call/{cid}")
    async def call(cid: str, req: Request):
        auth(req)
        try:
            return {"state": "done", "result": await modal.FunctionCall.from_id(cid).get.aio(timeout=0)}
        except TimeoutError:
            return {"state": "pending"}
        except Exception as e:
            return {"state": "failed", "error": f"{type(e).__name__}: {str(e)[:3000]}"}

    return web


@app.function(image=run_image, volumes={"/data": data_vol}, cpu=4.0, memory=8192, timeout=3600)
def fetch_data(files: dict):
    import concurrent.futures, urllib.request
    os.makedirs("/data/official", exist_ok=True)

    def get(item):
        name, url = item
        dst = f"/data/official/{name}"
        if os.path.exists(dst):
            return 0
        with urllib.request.urlopen(url, timeout=600) as r, open(dst + ".part", "wb") as f:
            shutil.copyfileobj(r, f, 1 << 22)
        os.replace(dst + ".part", dst); return os.path.getsize(dst)
    with concurrent.futures.ThreadPoolExecutor(32) as ex:
        sizes = list(ex.map(get, files.items()))
    data_vol.commit()
    return {"fetched": sum(1 for s in sizes if s), "bytes": sum(sizes), "total": len(os.listdir("/data/official"))}
