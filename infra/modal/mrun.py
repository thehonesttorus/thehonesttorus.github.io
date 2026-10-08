#!/usr/bin/env python3
"""Modal experiment runner, second design (app "claude-wb2"; supersedes mjob.py).

What changed against mjob.py, and why (Modal docs: guide/scale, guide/function-invocation-methods,
guide/managing-deployments, reference/modal.Volume):
- Code is not baked into the image. Each batch carries an immutable, content-addressed snapshot of
  workbench/{k3work,official,num12,ncgprob}/*.py on the volume (/data/code/<hash>/), uploaded once by the client.
  Editing code therefore never redeploys. Batches with different code run side by side, and no warm container can
  serve stale code. mjob.py's per-version function variants (one ephemeral secret each, 20+ variants) are gone.
- The image holds only the dependencies and changes only when this file's remote functions change. Deploys use
  `--strategy recreate`: under the default rolling strategy an old container keeps taking new inputs until its
  replacement is warm, which is where mjob.py's stale-code results came from.
- Containers are sized to the work: cpu=4 by default (the chain's measured use is 3.8 cores of 8). The workspace's
  container allowance (100 on Starter) is then the parallelism, at a quarter of the reserved cores.
- Every call is spawned at once and results print as they finish, not in submission order. A summary goes to
  scratchpad/modal/JOB/summary.txt.

Credentials come from the Modal client's own config (~/.modal.toml, mode 600, or MODAL_TOKEN_ID / MODAL_TOKEN_SECRET);
nothing here reads or prints them. Data lives on the volume whest-data:
  /data/official  W_off<id>.npy, truth_off<id>.npz
  /data/k3work    Monte Carlo files, chain dumps and other inputs that later jobs read
  /data/code/H    code snapshots (H = content hash)
  /data/out/JOB   log_<tag>.txt and every new file a job wrote

  python3 infra/modal/mrun.py deploy                      deploy now (only needed when this file changes)
  python3 infra/modal/mrun.py probe N [CPU] [SECONDS]     N sleeping calls: the peak number running at once
  python3 infra/modal/mrun.py run JOB 'CMD' [--cpu N] [--mem GB] [--timeout S] [--keep]
  python3 infra/modal/mrun.py fan JOB 'CMD with {net}' NETS [opts]     one container per net (NETS: 0-15 or 0,3,7)
  python3 infra/modal/mrun.py batch JOB FILE [opts]                    one container per line "tag<TAB>command"
  python3 infra/modal/mrun.py get JOB [PATTERN]                        copy /data/out/JOB to scratchpad/modal/JOB
  python3 infra/modal/mrun.py ls [DIR]                                 list a volume directory
"""
import concurrent.futures as cf, fcntl, fnmatch, hashlib, os, shutil, subprocess, sys, time
import modal

HERE = os.path.dirname(os.path.abspath(__file__))
WB = os.path.join(HERE, "..", "..", "workbench")
CODE_DIRS = ("k3work", "official", "num12", "ncgprob")
APP = "claude-wb2"
SCRATCH = os.environ.get("CLAUDE_MODAL_DIR", "/tmp/claude-0/-home-user-thehonesttorus-github-io/"
                         "2cab30f0-c454-5769-8d9a-8195cd80a5f3/scratchpad/modal")

app = modal.App(APP)
vol = modal.Volume.from_name("whest-data", create_if_missing=True)
image = (modal.Image.debian_slim(python_version="3.12")
         .apt_install("curl")
         .pip_install("numpy", "scipy", "pyarrow", "flopscope==0.12.1", "whestbench==0.16.1"))


def _workdir(code):
    vol.reload()
    snap = f"/data/code/{code}"
    if not os.path.isdir(snap):
        raise RuntimeError(f"code snapshot {code} is not on the volume")
    for d in CODE_DIRS:
        w = f"/work/{d}"; os.makedirs(w, exist_ok=True)
        if os.path.isdir(f"{snap}/{d}"):
            for f in os.listdir(f"{snap}/{d}"):
                shutil.copy2(f"{snap}/{d}/{f}", w)
        src = f"/data/{d}"
        if os.path.isdir(src):
            for f in os.listdir(src):
                if not f.endswith(".py") and not os.path.exists(f"{w}/{f}"):
                    os.symlink(f"{src}/{f}", f"{w}/{f}")
    return set(os.listdir("/work/k3work"))


def _threads(cpu):
    n = str(max(1, int(round(cpu))))
    return {k: n for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}


@app.function(image=image, volumes={"/data": vol}, cpu=4.0, memory=16384, timeout=3600)
def runcmd(job, tag, cmd, code, cpu=4.0, keep=False):
    before = _workdir(code)
    out = f"/data/out/{job}"; os.makedirs(out, exist_ok=True)
    env = dict(os.environ, OUT=out, PYTHONPATH="/work/num12", PYTHONWARNINGS="ignore", **_threads(cpu))
    t0 = time.time()
    p = subprocess.run(["bash", "-c", cmd], cwd="/work/k3work", env=env, capture_output=True, text=True)
    log = p.stdout + ("\n[stderr]\n" + p.stderr[-6000:] if p.stderr.strip() else "")
    log += f"\n[exit {p.returncode} after {time.time() - t0:.0f}s on {os.cpu_count()} visible cores, cpu={cpu}, code {code[:10]}]\n"
    for f in sorted(set(os.listdir("/work/k3work")) - before):
        src = f"/work/k3work/{f}"
        if os.path.isfile(src) and not os.path.islink(src) and not f.endswith(".py"):
            if keep:
                shutil.copy2(src, f"/data/k3work/{f}")
            shutil.move(src, f"{out}/{f}")
    open(f"{out}/log_{tag}.txt", "w").write(log)
    vol.commit()
    return tag, p.returncode, log


@app.function(image=image, cpu=1.0, memory=1024, timeout=900)
def hold(i, seconds):
    t0 = time.time()
    time.sleep(seconds)
    return i, os.environ.get("MODAL_TASK_ID", ""), t0, time.time()


@app.function(image=image, volumes={"/data": vol}, timeout=600)
def listdir(d):
    vol.reload()
    p = f"/data/{d}".rstrip("/")
    if not os.path.isdir(p):
        return f"{p}: missing"
    return "\n".join(f"{os.path.getsize(os.path.join(p, f)):>12d}  {f}" for f in sorted(os.listdir(p)))


# ------------------------------------------------------------------ client side
def _lock():
    os.makedirs(SCRATCH, exist_ok=True)
    lk = open(os.path.join(SCRATCH, "mrun.lock"), "w")
    fcntl.flock(lk, fcntl.LOCK_EX)                       # one deploy or upload at a time across parallel callers
    return lk


def ensure_deployed(force=False):
    stamp = os.path.join(SCRATCH, "mrun_deployed.txt")
    hc = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()
    with _lock():
        if not force and os.path.exists(stamp) and open(stamp).read().strip() == hc:
            return
        r = subprocess.run(["modal", "deploy", "--strategy", "recreate", os.path.abspath(__file__)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit("modal deploy failed:\n" + r.stdout[-2000:] + r.stderr[-2000:])
        open(stamp, "w").write(hc)
        print(f"deployed {APP} (runner {hc[:10]})", file=sys.stderr)


def ensure_code():
    """Upload the content-addressed code snapshot once; return its hash."""
    h, files = hashlib.sha256(), []
    for d in CODE_DIRS:
        p = os.path.join(WB, d)
        for f in sorted(os.listdir(p)):
            if f.endswith(".py"):
                h.update(d.encode()); h.update(f.encode()); h.update(open(os.path.join(p, f), "rb").read())
                files.append((d, f))
    code = h.hexdigest()[:20]
    stamp = os.path.join(SCRATCH, "mrun_code.txt")
    with _lock():
        if os.path.exists(stamp) and code in open(stamp).read().split():
            return code
        with vol.batch_upload(force=True) as b:
            for d, f in files:
                b.put_file(os.path.join(WB, d, f), f"/code/{code}/{d}/{f}")
        open(stamp, "a").write(code + "\n")
        print(f"uploaded code snapshot {code} ({len(files)} files)", file=sys.stderr)
    return code


def fn(name, **opts):
    for i in range(6):
        try:
            f = modal.Function.from_name(APP, name)
            return f.with_options(**opts) if opts else f
        except modal.exception.ResourceExhaustedError:
            time.sleep(2 ** i)
    raise SystemExit("Modal lookup kept hitting rate limits")


def _nets(s):
    out = []
    for part in s.split(","):
        a, _, b = part.partition("-")
        out += list(range(int(a), int(b) + 1)) if b else [int(a)]
    return out


def _opts(a):
    o = {"cpu": 4.0, "mem": 16, "timeout": 3600, "keep": "--keep" in a}
    for k in ("cpu", "mem", "timeout"):
        if f"--{k}" in a:
            o[k] = float(a[a.index(f"--{k}") + 1])
    return o


def _function(o):
    if (o["cpu"], o["mem"], o["timeout"]) == (4.0, 16, 3600):
        return fn("runcmd")
    return fn("runcmd", cpu=o["cpu"], memory=int(o["mem"] * 1024), timeout=int(o["timeout"]))


def _line(t, rc, log):
    body = log.split("\n[stderr]\n")[0]
    tail = [l for l in body.splitlines() if l.strip()][-2:]
    if rc != 0:
        tail = [l for l in log.splitlines() if l.strip()][-4:]
    return f"[{t} exit {rc}] " + " | ".join(tail)


def run_jobs(job, jobs, o):
    """Spawn every (tag, cmd) at once; print each result as it finishes; write a summary."""
    code = ensure_code()
    f = _function(o)
    t0 = time.time()
    calls = [(t, f.spawn(job, t, cmd, code, o["cpu"], o["keep"])) for t, cmd in jobs]
    print(f"{job}: spawned {len(calls)} calls in {time.time() - t0:.1f}s (code {code[:10]}, cpu {o['cpu']})", flush=True)
    lines = {}
    with cf.ThreadPoolExecutor(max_workers=min(len(calls), 128)) as ex:
        futs = {ex.submit(fc.get): t for t, fc in calls}
        for fu in cf.as_completed(futs):
            t = futs[fu]
            try:
                s = _line(*fu.result())
            except Exception as e:                    # a container failure must not hide the other results
                s = f"[{t} FAILED] {type(e).__name__}: {str(e)[:300]}"
            lines[t] = s
            print(s, f"(+{time.time() - t0:.0f}s)", flush=True)
    d = os.path.join(SCRATCH, job); os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "summary.txt"), "w").write("\n".join(lines[t] for t, _ in jobs if t in lines) + "\n")
    print(f"{job}: {len(lines)} results in {time.time() - t0:.0f}s", flush=True)


def probe(n, cpu, secs):
    f = fn("hold", cpu=cpu, memory=int(max(1024, 1024 * cpu)))
    t0 = time.time()
    calls = [f.spawn(i, secs) for i in range(n)]
    res = [c.get() for c in calls]
    ev = sorted([(r[2], 1) for r in res] + [(r[3], -1) for r in res])
    cur = peak = 0
    for _, dl in ev:
        cur += dl; peak = max(peak, cur)
    st = sorted(r[2] - t0 for r in res)
    print(f"probe: {n} calls of {secs:.0f}s at cpu={cpu}: peak {peak} running at once, {len({r[1] for r in res})} containers; "
          f"start offsets min {st[0]:.0f}s median {st[len(st) // 2]:.0f}s max {st[-1]:.0f}s; wall {time.time() - t0:.0f}s")


def get(job, pattern="*"):
    d = os.path.join(SCRATCH, job); os.makedirs(d, exist_ok=True)
    names = []
    for e in vol.listdir(f"/out/{job}"):
        name = os.path.basename(e.path)
        if fnmatch.fnmatch(name, pattern):
            with open(os.path.join(d, name), "wb") as fh:
                for chunk in vol.read_file(e.path):
                    fh.write(chunk)
            names.append(name)
    print(f"{job}: {len(names)} files -> {d}")


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a:
        sys.exit(__doc__)
    c = a[0]
    if c == "get":
        get(a[1], a[2] if len(a) > 2 and not a[2].startswith("--") else "*"); sys.exit()
    ensure_deployed(force=(c == "deploy"))
    if c == "deploy":
        pass
    elif c == "probe":
        probe(int(a[1]), float(a[2]) if len(a) > 2 else 1.0, float(a[3]) if len(a) > 3 else 60.0)
    elif c == "ls":
        print(fn("listdir").remote(a[1] if len(a) > 1 else ""))
    elif c == "run":
        o = _opts(a); run_jobs(a[1], [("0", a[2])], o)
        print(open(os.path.join(SCRATCH, a[1], "summary.txt")).read())
    elif c == "fan":
        run_jobs(a[1], [(str(n), a[2].replace("{net}", str(n))) for n in _nets(a[3])], _opts(a))
    elif c == "batch":
        run_jobs(a[1], [l.rstrip("\n").split("\t", 1) for l in open(a[2]) if l.strip() and not l.startswith("#")], _opts(a))
    else:
        sys.exit(__doc__)
