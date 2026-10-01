#!/usr/bin/env python3
"""Check every logged az call (one JSON argv per line, from fake_az.py) against the help of the
REAL installed az CLI: the command group/subcommand must exist and every --flag / -x must be one
of its arguments.  Help pages are fetched once per command (3 in parallel) and cached.

    python3 check_az_calls.py <calls.log> [--az /path/to/real/az] [--cache DIR]
"""
import argparse, concurrent.futures as cf, json, os, re, subprocess, sys

ap = argparse.ArgumentParser()
ap.add_argument("log")
ap.add_argument("--az", default=os.environ.get("REAL_AZ", "az"))
ap.add_argument("--cache", default=os.path.join(os.environ.get("TMPDIR", "/tmp"), "azhelp-cache"))
args = ap.parse_args()
os.makedirs(args.cache, exist_ok=True)

calls = [json.loads(l) for l in open(args.log) if l.strip()]
def split(argv):
    path = []
    for t in argv:
        if t.startswith("-"):
            break
        path.append(t)
    flags = [t.split("=", 1)[0] for t in argv[len(path):] if t.startswith("-") and not re.match(r"^-\d", t)]
    return " ".join(path), flags

def help_for(cmd):
    fn = os.path.join(args.cache, cmd.replace(" ", "_") + ".txt")
    if not os.path.exists(fn):
        p = subprocess.run([args.az, *cmd.split(), "--help"], capture_output=True, text=True,
                           env={**os.environ, "OPENBLAS_NUM_THREADS": "1"})
        with open(fn, "w") as f:
            f.write(f"rc={p.returncode}\n" + p.stdout + p.stderr)
    text = open(fn).read()
    rc = int(text.split("\n", 1)[0][3:])
    opts = set()
    for line in text.splitlines():
        m = re.match(r"^\s{4}(-[^:]*?)\s*(?::|$)", line)
        if m:
            opts.update(t for t in m.group(1).split() if t.startswith("-"))
    ok = rc == 0 and "Command" in text and "is misspelled or not recognized" not in text
    return ok, opts

cmds = sorted({split(c)[0] for c in calls})
with cf.ThreadPoolExecutor(3) as ex:
    helps = dict(zip(cmds, ex.map(help_for, cmds)))

bad = 0
seen = set()
for argv in calls:
    cmd, flags = split(argv)
    ok, opts = helps[cmd]
    key = (cmd, tuple(sorted(set(flags))))
    if key in seen:
        continue
    seen.add(key)
    unknown = [f for f in flags if f not in opts]
    if not ok or unknown:
        bad += 1
        print(f"FAIL az {cmd}: " + ("command not found; " if not ok else "") + (f"unknown flags {unknown}" if unknown else ""))
    else:
        print(f"ok   az {cmd} {' '.join(sorted(set(flags)))}")
print(f"{len(calls)} calls, {len(seen)} distinct command+flag sets, {bad} failing")
sys.exit(1 if bad else 0)
