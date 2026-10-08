# Bake a configuration into a single estimator file (note XXXIX): every os.environ.get("NAME", "default") switch named in
# the configuration gets the configured value as its default, so the file runs with no environment at all.
#   python bake.py SRC OUT CONFIG_FILE
# CONFIG_FILE: whitespace-separated NAME=VALUE pairs (values may not contain whitespace; V47_CAL strings do not).
import re, sys

src, out, cfg = sys.argv[1], sys.argv[2], sys.argv[3]
pairs = dict(kv.split("=", 1) for kv in open(cfg).read().split())
code = open(src).read()
done = set()
for name, val in pairs.items():
    pat = re.compile(r'_os\.environ\.get\(\s*"' + re.escape(name) + r'"\s*,\s*"[^"]*"\s*\)')
    code, k = pat.subn(f'_os.environ.get("{name}", "{val}")', code)
    if k == 0:
        sys.exit(f"switch {name} not found as _os.environ.get(\"{name}\", \"...\") in {src}")
    done.add(name)
open(out, "w").write(code)
print(f"baked {len(done)} switches into {out}")
