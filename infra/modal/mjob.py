#!/usr/bin/env python3
"""Modal experiment runner: a third compute pool beside infra/aws.

Credentials come from the Modal client's own config (~/.modal.toml, mode 600, or MODAL_TOKEN_ID / MODAL_TOKEN_SECRET);
nothing here reads or prints them. The app is deployed once as "claude-whest" and called by name, so parallel callers
create no ephemeral apps (Modal rate-limits app creation). The code (workbench/{k3work,official,num12,ncgprob}/*.py)
is baked in at deploy time; every call first checks a hash of that code and redeploys when it changed.
Data lives on the volume whest-data:
  /data/official  W_off<id>.npy, truth_off<id>.npz (the dataset's own all-layer means, so output MSE needs no Monte Carlo)
  /data/k3work    Monte Carlo files, chain dumps and other inputs that later jobs read
  /data/out/JOB   log_<tag>.txt and every new file a job wrote
Each call builds /work/{k3work,official,num12,ncgprob} from the code plus symlinks to the volume's data, runs the command
in /work/k3work with OUT=/data/out/JOB, and moves new non-code files to /data/out/JOB (and, with --keep, also to
/data/k3work so that later jobs can read them).

  python3 infra/modal/mjob.py deploy                               (re)deploy now
  python3 infra/modal/mjob.py bench                                cores, BLAS, a 1024^3 matmul timing
  python3 infra/modal/mjob.py prep 0 1 2 3 4 5 6                   download official shards into the volume (one container each)
  python3 infra/modal/mjob.py run JOB 'CMD' [--cpu N] [--mem GB] [--timeout S] [--keep]
  python3 infra/modal/mjob.py fan JOB 'CMD with {net}' NETS [opts] one container per net (NETS: 0-15 or 0,3,7)
  python3 infra/modal/mjob.py batch JOB FILE [opts]                one container per line "tag<TAB>command" of FILE
  python3 infra/modal/mjob.py get JOB [PATTERN]                    copy /data/out/JOB to scratchpad/modal/JOB, print the logs
  python3 infra/modal/mjob.py ls [DIR]                             list a volume directory
"""
import fcntl, fnmatch, hashlib, os, shutil, subprocess, sys, time
import modal

HERE = os.path.dirname(os.path.abspath(__file__))
WB = os.path.join(HERE, "..", "..", "workbench")
CODE_DIRS = ("k3work", "official", "num12", "ncgprob")
APP = "claude-whest"
SCRATCH = os.environ.get("CLAUDE_MODAL_DIR", "/tmp/claude-0/-home-user-thehonesttorus-github-io/"
                         "2cab30f0-c454-5769-8d9a-8195cd80a5f3/scratchpad/modal")
HF = "https://huggingface.co/datasets/aicrowd/arc-whestbench-public-2026/resolve/v2-phase2/data/mini-0000{k}-of-00007.parquet"

app = modal.App(APP)
vol = modal.Volume.from_name("whest-data", create_if_missing=True)
image = (modal.Image.debian_slim(python_version="3.12")
         .apt_install("curl")
         .pip_install("numpy", "scipy", "pyarrow", "flopscope==0.12.1", "whestbench==0.16.1"))
for _d in CODE_DIRS:
    image = image.add_local_dir(os.path.join(WB, _d), f"/code/{_d}", ignore=lambda p: not str(p).endswith(".py"))


def _workdir():
    vol.reload()
    for d in CODE_DIRS:
        w = f"/work/{d}"; os.makedirs(w, exist_ok=True)
        for f in os.listdir(f"/code/{d}"):
            shutil.copy2(f"/code/{d}/{f}", w)
        src = f"/data/{d}"
        if os.path.isdir(src):
            for f in os.listdir(src):
                if not f.endswith(".py") and not os.path.exists(f"{w}/{f}"):
                    os.symlink(f"{src}/{f}", f"{w}/{f}")
    return set(os.listdir("/work/k3work"))


def _threads(cpu):
    n = str(max(1, int(cpu)))
    return {k: n for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}


def _wb_hash(root):
    h = hashlib.sha256()
    for d in CODE_DIRS:
        p = os.path.join(root, d)
        for f in sorted(os.listdir(p)):
            if f.endswith(".py"):
                h.update(f.encode()); h.update(open(os.path.join(p, f), "rb").read())
    return h.hexdigest()


@app.function(image=image, volumes={"/data": vol}, cpu=16.0, memory=32768, timeout=3600)
def runcmd(job, tag, cmd, cpu=16.0, keep=False, code=None):
    if code is not None and _wb_hash("/code") != code:      # a warm container of an older deployment: refuse
        return tag, 97, f"[stale code in container: expected {code[:10]}, have {_wb_hash('/code')[:10]}; rerun]\n"
    before = _workdir()
    out = f"/data/out/{job}"; os.makedirs(out, exist_ok=True)
    env = dict(os.environ, OUT=out, PYTHONPATH="/work/num12", PYTHONWARNINGS="ignore", **_threads(cpu))
    t0 = time.time()
    p = subprocess.run(["bash", "-c", cmd], cwd="/work/k3work", env=env, capture_output=True, text=True)
    log = p.stdout + ("\n[stderr]\n" + p.stderr[-6000:] if p.stderr.strip() else "")
    log += f"\n[exit {p.returncode} after {time.time() - t0:.0f}s on {os.cpu_count()} visible cores, cpu={cpu}]\n"
    for f in sorted(set(os.listdir("/work/k3work")) - before):
        src = f"/work/k3work/{f}"
        if os.path.isfile(src) and not os.path.islink(src) and not f.endswith(".py"):
            if keep:
                shutil.copy2(src, f"/data/k3work/{f}")
            shutil.move(src, f"{out}/{f}")
    open(f"{out}/log_{tag}.txt", "w").write(log)
    vol.commit()
    return tag, p.returncode, log


@app.function(image=image, cpu=8.0, memory=8192, timeout=600)
def bench():
    import numpy as np, io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        np.show_config()
    a = np.random.default_rng(0).standard_normal((1024, 1024)).astype(np.float32)
    t0 = time.time()
    for _ in range(20):
        a @ a
    dt = (time.time() - t0) / 20
    return f"cores visible {os.cpu_count()}; 1024^3 float32 matmul {dt * 1e3:.1f} ms ({2 * 1024 ** 3 / dt / 1e9:.0f} GFLOP/s)"


@app.function(image=image, volumes={"/data": vol}, cpu=4.0, memory=16384, timeout=3600)
def prep(k):
    _workdir()
    os.makedirs("/data/official", exist_ok=True); os.makedirs("/data/k3work", exist_ok=True)
    pq = f"/tmp/mini{k}.parquet"
    subprocess.run(["curl", "-sSL", "--retry", "4", "-o", pq, HF.format(k=k)], check=True)
    p = subprocess.run(["python", "/work/official/extract.py", pq], cwd="/data/official", capture_output=True, text=True)
    os.remove(pq); vol.commit()
    return k, p.returncode, p.stdout.strip(), p.stderr[-2000:]


@app.function(image=image, volumes={"/data": vol}, timeout=600)
def listdir(d):
    vol.reload()
    p = f"/data/{d}".rstrip("/")
    if not os.path.isdir(p):
        return f"{p}: missing"
    return "\n".join(f"{os.path.getsize(os.path.join(p, f)):>12d}  {f}" for f in sorted(os.listdir(p)))


# ------------------------------------------------------------------ client side
def _code_hash():
    h = hashlib.sha256(open(os.path.abspath(__file__), "rb").read())
    for d in CODE_DIRS:
        p = os.path.join(WB, d)
        for f in sorted(os.listdir(p)):
            if f.endswith(".py"):
                h.update(f.encode()); h.update(open(os.path.join(p, f), "rb").read())
    return h.hexdigest()


def ensure_deployed(force=False):
    os.makedirs(SCRATCH, exist_ok=True)
    stamp = os.path.join(SCRATCH, "deployed_hash.txt")
    with open(os.path.join(SCRATCH, "deploy.lock"), "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)                  # one deploy at a time across parallel callers
        hc = _code_hash()
        if not force and os.path.exists(stamp) and open(stamp).read().strip() == hc:
            return
        r = subprocess.run(["modal", "deploy", os.path.abspath(__file__)], capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit("modal deploy failed:\n" + r.stdout[-2000:] + r.stderr[-2000:])
        open(stamp, "w").write(hc)
        print(f"deployed {APP} (code {hc[:10]})", file=sys.stderr)


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
    o = {"cpu": 16.0, "mem": 32, "timeout": 3600, "keep": "--keep" in a}
    for k in ("cpu", "mem", "timeout"):
        if f"--{k}" in a:
            o[k] = float(a[a.index(f"--{k}") + 1])
    return o


def _run(o):
    # the code hash in env gives every code version its own container pool: no warm container of an older
    # deployment can serve the call (Modal pools containers per option set)
    return fn("runcmd", cpu=o["cpu"], memory=int(o["mem"] * 1024), timeout=int(o["timeout"]), env={"WB_CODE": _wb_hash(WB)})


def _report(t, rc, log):
    body = log.split("\n[stderr]\n")[0]
    tail = [l for l in body.splitlines() if l.strip()][-2:] + [log.strip().splitlines()[-1]]
    if rc != 0:
        tail = [l for l in log.splitlines() if l.strip()][-4:]
    print(f"[{t} exit {rc}] " + " | ".join(tail), flush=True)


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
    for n in sorted(names):
        if n.startswith("log_"):
            print(f"--- {n}\n" + open(os.path.join(d, n)).read()[-3000:])


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
    elif c == "bench":
        print(fn("bench").remote())
    elif c == "prep":
        for k, rc, so, se in fn("prep").map([int(x) for x in a[1:]]):
            print(f"shard {k}: exit {rc} nets {so} {se.strip()[-300:]}")
    elif c == "ls":
        print(fn("listdir").remote(a[1] if len(a) > 1 else ""))
    elif c == "run":
        o = _opts(a); tag, rc, log = _run(o).remote(a[1], "0", a[2], o["cpu"], o["keep"], _wb_hash(WB))
        print(log[-4000:])
    elif c in ("batch", "fan"):
        o = _opts(a); f = _run(o)
        if c == "batch":
            jobs = [l.rstrip("\n").split("\t", 1) for l in open(a[2]) if l.strip() and not l.startswith("#")]
        else:
            jobs = [(str(n), a[2].replace("{net}", str(n))) for n in _nets(a[3])]
        calls = [(t, f.spawn(a[1], t, cmd, o["cpu"], o["keep"], _wb_hash(WB))) for t, cmd in jobs]
        for t, fc in calls:
            try:
                _report(*fc.get())
            except Exception as ex:                      # a container failure must not hide the other results
                print(f"[{t} FAILED] {type(ex).__name__}: {str(ex)[:300]}", flush=True)
    else:
        sys.exit(__doc__)
