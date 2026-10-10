#!/usr/bin/env python3
"""Score an estimator on the official networks, all in parallel on Modal, with results streamed as they finish.

  python systems/score.py JOB EST_FILE [--nets 0-99] [--cpu 8] [--mem 16384] [--save-out] [--env K=V ...] [--base JOB0]

Each network runs the scored regime (systems/v56/run_scored.py: a warm-up predict on another network, then the measured
predict). Every finished network prints its line at once, followed by the running summary. At the end
RESULTS/JOB/summary.json holds the per-network rows. With --base, the summary is paired against an earlier job.
"""
import json, os, re, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "infra"))
import gw

PAT = re.compile(r"net (\d+) \S+: raw (\S+)\s+C/B (\S+)\s+adjusted (\S+)\s+wall (\S+)s\s+residual (\S+)s")


def nets(spec):
    out = []
    for part in spec.split(","):
        a, _, b = part.partition("-"); out += list(range(int(a), int(b) + 1)) if b else [int(a)]
    return out


def summary(rows, base=None):
    k = sorted(rows); raw = np.array([rows[i]["raw"] for i in k]); cb = np.array([rows[i]["cb"] for i in k])
    adj = raw * np.maximum(0.1, cb)
    s = (f"n={len(k)}  raw {raw.mean():.4e}  C/B {cb.mean():.4f}  adjusted {adj.mean():.4e} (mean raw x mean C/B "
         f"{raw.mean() * max(0.1, cb.mean()):.4e})  wall max {max(rows[i]['wall'] for i in k):.0f}s")
    if base:
        common = [i for i in k if i in base]
        if common:
            d = np.array([np.log(rows[i]["raw"] * max(.1, rows[i]["cb"]) / (base[i]["raw"] * max(.1, base[i]["cb"]))) for i in common])
            s += f"  | vs base on {len(common)}: adjusted {100 * (np.exp(d.mean()) - 1):+.2f}% +- {100 * d.std(ddof=1) / np.sqrt(len(d)) if len(d) > 1 else 0:.2f} (geo), better on {(d < 0).sum()}/{len(d)}"
    return s


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) < 2:
        sys.exit(__doc__)
    job, est = a[0], a[1]
    opt = lambda k, d: type(d)(a[a.index(k) + 1]) if k in a else d
    env = dict(a[i + 1].split("=", 1) for i, x in enumerate(a) if x == "--env")
    base = None
    if "--base" in a:
        base = {int(k): v for k, v in json.load(open(os.path.join(gw.RESULTS, opt("--base", ""), "summary.json")))["rows"].items()}
    tasks = []
    for i in nets(opt("--nets", "0-99")):
        save = f"SAVE_OUT=$OUT/out_{i}.npy " if "--save-out" in a else ""
        tasks.append({"tag": f"n{i:03d}", "cmd": f"{save}python systems/v56/run_scored.py {i} {job} {est}"})
    rows = {}
    for tag, res in gw.map(job, tasks, cpu=opt("--cpu", 8.0), mem=opt("--mem", 16384), timeout=opt("--timeout", 3600), env=env, quiet=True):
        m = PAT.search(res.get("stdout", "")) if "error" not in res else None
        if not m:
            print(f"{tag}: FAILED {res.get('error') or res.get('stderr', '')[-600:]}", flush=True); continue
        i = int(m.group(1)); rows[i] = dict(raw=float(m.group(2)), cb=float(m.group(3)), wall=float(m.group(5)), resid=float(m.group(6)), secs=res["secs"])
        print(f"{tag}: raw {rows[i]['raw']:.4e} C/B {rows[i]['cb']:.4f} wall {rows[i]['wall']:.0f}s task {res['secs']:.0f}s   || {summary(rows, base)}", flush=True)
    if rows:
        json.dump({"est": est, "env": env, "rows": rows, "summary": summary(rows, base)},
                  open(os.path.join(gw.RESULTS, job, "summary.json"), "w"), indent=1)
        print("FINAL", summary(rows, base), flush=True)
