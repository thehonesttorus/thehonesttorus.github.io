#!/usr/bin/env python3
"""Modal app "whest9" for this branch. The sandbox's proxy cannot carry gRPC, so the client side of this file runs on
the AWS relay instance (infra/modal_relay.py ships the code there and drives it over SSM). Volume whest-data (shared
with the earlier sessions): /data/official (W_off*.npy, truth_off*.npz), /data/code9/<hash>/ (code snapshots),
/data/out9/JOB/ (results: log_<tag>.txt and every file a task wrote to $OUT).
  python infra/modal_app.py deploy
  python infra/modal_app.py batch JOB FILE [--cpu N] [--mem GB] [--timeout S]   lines "tag<TAB>command"; the command
        runs in /work with whest/ and scripts/ copied in, DATA=/data/official, OUT=/work/out, PYTHONPATH=/work.
  python infra/modal_app.py mc JOB NET N SEED0 NSEEDS [--gpu L4|A10G|A100] [--cov l1,l2]   GPU Monte Carlo moments
  python infra/modal_app.py get JOB DEST
"""
import hashlib, io, os, shutil, subprocess, sys, tarfile, time
import modal
APP = "whest9"; HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, ".."); CODE_DIRS = ("whest", "scripts")
app = modal.App(APP)
vol = modal.Volume.from_name("whest-data")
image = modal.Image.debian_slim(python_version="3.12").pip_install("numpy", "scipy", "pyarrow", "flopscope==0.12.1", "whestbench==0.16.1")
gpu_image = modal.Image.debian_slim(python_version="3.12").pip_install("numpy", "torch")


def _workdir(code):
    vol.reload(); snap = f"/data/code9/{code}"
    if not os.path.isdir(snap):
        raise RuntimeError(f"code snapshot {code} missing on the volume")
    if os.path.isdir("/work"): shutil.rmtree("/work")
    shutil.copytree(snap, "/work"); os.makedirs("/work/out", exist_ok=True)


@app.function(image=image, volumes={"/data": vol}, cpu=8.0, memory=16384, timeout=3600)
def runcmd(job, tag, cmd, code, cpu=8.0):
    _workdir(code); out = f"/data/out9/{job}"; os.makedirs(out, exist_ok=True)
    th = str(max(1, int(round(cpu))))
    env = dict(os.environ, OUT="/work/out", DATA="/data/official", PYTHONPATH="/work", PYTHONWARNINGS="ignore",
               OMP_NUM_THREADS=th, OPENBLAS_NUM_THREADS=th, MKL_NUM_THREADS=th)
    t0 = time.time()
    p = subprocess.run(["bash", "-c", cmd], cwd="/work", env=env, capture_output=True, text=True)
    log = p.stdout + ("\n[stderr]\n" + p.stderr[-6000:] if p.stderr.strip() else "")
    log += f"\n[exit {p.returncode} after {time.time() - t0:.0f}s on modal cpu={cpu} ({os.cpu_count()} visible), code {code[:10]}]\n"
    for f in os.listdir("/work/out"):
        shutil.move(f"/work/out/{f}", f"{out}/{f}")
    open(f"{out}/log_{tag}.txt", "w").write(log); vol.commit()
    return tag, p.returncode, log


@app.function(image=image, volumes={"/data": vol}, timeout=600)
def put_code(code, files):
    snap = f"/data/code9/{code}"
    for rel, data in files.items():
        p = f"{snap}/{rel}"; os.makedirs(os.path.dirname(p), exist_ok=True); open(p, "wb").write(data)
    vol.commit(); return code


@app.function(image=gpu_image, volumes={"/data": vol}, gpu="L4", timeout=3600)
def mc_gpu(job, net, N, seed, cov_layers=()):
    """Monte Carlo of pre-activation moments on a GPU: per-neuron s1..s4, gate counts, post-activation sums h1, h2,
    and E[z z^T] sums for cov_layers. float32 forward passes in batches, float64 accumulation."""
    import numpy as np, torch
    dev = "cuda"; W = np.load(f"/data/official/W_off{net}.npy")          # (L, n_out, n_in), z = W h
    Wt = [torch.from_numpy(np.ascontiguousarray(w.T)).to(dev) for w in W]  # h @ W^T as (B, n_in) @ (n_in, n_out)
    L, n = len(Wt), Wt[0].shape[0]; B = 32768
    S = {k: torch.zeros((L, n), dtype=torch.float64, device=dev) for k in ("s1", "s2", "s3", "s4", "pos", "h1", "h2")}
    C = {l: torch.zeros((n, n), dtype=torch.float64, device=dev) for l in cov_layers}
    g = torch.Generator(device=dev); g.manual_seed(int(seed)); t0 = time.time(); done = 0
    while done < N:
        h = torch.randn((B, n), generator=g, device=dev)
        for l in range(L):
            z = h @ Wt[l]; h = torch.relu(z); z64 = z.double(); z2 = z64 * z64
            S["s1"][l] += z64.sum(0); S["s2"][l] += z2.sum(0); S["s3"][l] += (z2 * z64).sum(0); S["s4"][l] += (z2 * z2).sum(0)
            S["pos"][l] += (z > 0).sum(0); h64 = h.double(); S["h1"][l] += h64.sum(0); S["h2"][l] += (h64 * h64).sum(0)
            if l in C: C[l] += z64.T @ z64
        done += B
    out = f"/data/out9/{job}"; os.makedirs(out, exist_ok=True)
    np.savez(f"{out}/mcc_{net}_{seed}.npz", N=done, **{k: v.cpu().numpy() for k, v in S.items()}, **{f"cov{l}": C[l].cpu().numpy() for l in C})
    vol.commit(); return net, seed, done, time.time() - t0


@app.function(image=image, volumes={"/data": vol}, timeout=600)
def listdir(d):
    vol.reload(); p = f"/data/{d}".rstrip("/")
    return "\n".join(f"{os.path.getsize(os.path.join(p, f)):>12d}  {f}" for f in sorted(os.listdir(p))) if os.path.isdir(p) else f"{p}: missing"


# ---------------------------------------------------------------- client side (runs on the relay)
def _code_files():
    files = {}; h = hashlib.sha256()
    for d in CODE_DIRS:
        for f in sorted(os.listdir(os.path.join(ROOT, d))):
            if f.endswith(".py"):
                b = open(os.path.join(ROOT, d, f), "rb").read(); files[f"{d}/{f}"] = b; h.update(f"{d}/{f}".encode()); h.update(b)
    return h.hexdigest()[:20], files


def _fn(name, **opts):
    f = modal.Function.from_name(APP, name); return f.with_options(**opts) if opts else f


def _arg(a, k, d, conv=float):
    return conv(a[a.index(f"--{k}") + 1]) if f"--{k}" in a else d


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a: sys.exit(__doc__)
    c = a[0]
    if c == "deploy":
        subprocess.run(["modal", "deploy", "--strategy", "recreate", os.path.abspath(__file__)], check=True)
    elif c == "batch":
        job, path = a[1], a[2]; cpu = _arg(a, "cpu", 8.0); mem = _arg(a, "mem", 16.0); to = _arg(a, "timeout", 3600, int)
        code, files = _code_files(); _fn("put_code").remote(code, files)
        tasks = [l.rstrip("\n").split("\t", 1) for l in open(path) if l.strip() and not l.startswith("#")]
        f = _fn("runcmd", cpu=cpu, memory=int(mem * 1024), timeout=to)
        calls = [(t, f.spawn(job, t, cmd, code, cpu)) for t, cmd in tasks]; t0 = time.time()
        print(f"{job}: {len(calls)} tasks spawned (cpu {cpu}, code {code})", flush=True)
        for t, fc in calls:
            try:
                tag, rc, log = fc.get(); tail = [l for l in log.splitlines() if l.strip()][-2:]
                print(f"[{tag} exit {rc}] " + " | ".join(tail) + f" (+{time.time()-t0:.0f}s)", flush=True)
            except Exception as ex:
                print(f"[{t} FAILED] {type(ex).__name__}: {str(ex)[:300]}", flush=True)
    elif c == "mc":
        job, net, N, seed0, nseeds = a[1], int(a[2]), int(float(a[3])), int(a[4]), int(a[5])
        gpu = _arg(a, "gpu", "L4", str); cov = [int(x) for x in _arg(a, "cov", "", str).split(",") if x]
        f = _fn("mc_gpu", gpu=gpu)
        calls = [f.spawn(job, net, N, seed0 + s, cov) for s in range(nseeds)]
        for fc in calls:
            try: print("done", fc.get(), flush=True)
            except Exception as ex: print("FAILED", str(ex)[:300], flush=True)
    elif c == "get":
        job, dest = a[1], a[2]; os.makedirs(dest, exist_ok=True)
        for e in vol.listdir(f"/out9/{job}"):
            name = os.path.basename(e.path)
            with open(os.path.join(dest, name), "wb") as fh:
                for chunk in vol.read_file(e.path): fh.write(chunk)
        print(f"{job}: fetched to {dest}")
    elif c == "ls":
        print(_fn("listdir").remote(a[1] if len(a) > 1 else ""))
    else:
        sys.exit(__doc__)
