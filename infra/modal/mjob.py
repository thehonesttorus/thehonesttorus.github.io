#!/usr/bin/env python3
"""Modal experiment runner: a third compute pool beside infra/aws.

Credentials come from the Modal client's own config (~/.modal.toml, mode 600, or MODAL_TOKEN_ID / MODAL_TOKEN_SECRET);
nothing here reads or prints them. Code ships with every call (workbench/{k3work,official,num12,ncgprob}/*.py). Data
lives on the volume whest-data:
  /data/official  W_off<id>.npy, truth_off<id>.npz (the dataset's own all-layer means, so output MSE needs no Monte Carlo)
  /data/k3work    Monte Carlo files, chain dumps and other inputs that later jobs read
  /data/out/JOB   log_<tag>.txt and every new file a job wrote
Each call builds /work/{k3work,official,num12,ncgprob} from the shipped code plus symlinks to the volume's data, runs the
command in /work/k3work with OUT=/data/out/JOB, and moves new non-code files to /data/out/JOB (and, with --keep, also
to /data/k3work so that later jobs can read them).

  python3 infra/modal/mjob.py bench                                cores, BLAS, a 1024^3 matmul timing
  python3 infra/modal/mjob.py prep 0 1 2 3 4 5 6                    download official shards into the volume (one container each)
  python3 infra/modal/mjob.py run JOB 'CMD' [--cpu N] [--mem GB] [--timeout S] [--keep]
  python3 infra/modal/mjob.py fan JOB 'CMD with {net}' NETS [opts]  one container per net (NETS: 0-15 or 0,3,7)
  python3 infra/modal/mjob.py batch JOB FILE [opts]                one container per line "tag<TAB>command" of FILE
  python3 infra/modal/mjob.py get JOB [PATTERN]                    copy /data/out/JOB to scratchpad/modal/JOB and print the logs
  python3 infra/modal/mjob.py ls [DIR]                             list a volume directory
"""
import fnmatch, os, shutil, subprocess, sys, time
import modal

WB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "workbench")
CODE_DIRS = ("k3work", "official", "num12", "ncgprob")
SCRATCH = os.environ.get("CLAUDE_MODAL_DIR", "/tmp/claude-0/-home-user-thehonesttorus-github-io/"
                         "2cab30f0-c454-5769-8d9a-8195cd80a5f3/scratchpad/modal")
HF = "https://huggingface.co/datasets/aicrowd/arc-whestbench-public-2026/resolve/v2-phase2/data/mini-0000{k}-of-00007.parquet"

app = modal.App("claude-whest")
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


@app.function(image=image, volumes={"/data": vol}, cpu=16.0, memory=32768, timeout=3600)
def runcmd(job, tag, cmd, cpu=16.0, keep=False):
    before = _workdir()
    out = f"/data/out/{job}"; os.makedirs(out, exist_ok=True)
    env = dict(os.environ, OUT=out, PYTHONPATH="/work/num12", **_threads(cpu))
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
    return f"cores visible {os.cpu_count()}; 1024^3 float32 matmul {dt * 1e3:.1f} ms ({2 * 1024 ** 3 / dt / 1e9:.0f} GFLOP/s)\n" + \
        "\n".join(l for l in buf.getvalue().splitlines() if "openblas" in l.lower() or "name" in l.lower())[:800]


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


def _fn(o):
    return runcmd.with_options(cpu=o["cpu"], memory=int(o["mem"] * 1024), timeout=int(o["timeout"]))


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
    with modal.enable_output() if "--verbose" in a else open(os.devnull) as _, app.run():
        if c == "bench":
            print(bench.remote())
        elif c == "prep":
            for k, rc, so, se in prep.map([int(x) for x in a[1:]]):
                print(f"shard {k}: exit {rc} nets {so} {se.strip()[-300:]}")
        elif c == "ls":
            print(listdir.remote(a[1] if len(a) > 1 else ""))
        elif c == "run":
            o = _opts(a); tag, rc, log = _fn(o).remote(a[1], "0", a[2], o["cpu"], o["keep"])
            print(log[-4000:])
        elif c == "batch":                     # FILE: one "tag<TAB>command" per line, each in its own container
            o = _opts(a); f = _fn(o)
            jobs = [l.rstrip("\n").split("\t", 1) for l in open(a[2]) if l.strip() and not l.startswith("#")]
            calls = [(t, f.spawn(a[1], t, cmd, o["cpu"], o["keep"])) for t, cmd in jobs]
            for t, fc in calls:
                tag, rc, log = fc.get()
                body = log.split("\n[stderr]\n")[0]
                tail = [l for l in body.splitlines() if l.strip()][-2:] + [log.strip().splitlines()[-1]]
                print(f"[{t} exit {rc}] " + " | ".join(tail), flush=True)
        elif c == "fan":
            o = _opts(a); nets = _nets(a[3]); f = _fn(o)
            calls = [(n, f.spawn(a[1], str(n), a[2].replace("{net}", str(n)), o["cpu"], o["keep"])) for n in nets]
            for n, fc in calls:
                tag, rc, log = fc.get()
                body = log.split("\n[stderr]\n")[0]
                tail = [l for l in body.splitlines() if l.strip()][-2:] + [log.strip().splitlines()[-1]]
                print(f"[net {n} exit {rc}] " + " | ".join(tail), flush=True)
        else:
            sys.exit(__doc__)
