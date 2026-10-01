#!/usr/bin/env python3
"""Offline stand-in for the `az` CLI, used by run_offline.sh.  Logs every call (argv as one JSON
line) to $FAKE_AZ_LOG, copies every --json-file payload to $FAKE_AZ_CAPTURE, and answers with
canned, realistically shaped output (camelCase keys as az prints them; --query applied with
jmespath when available; -o json|tsv|none).  Nothing here talks to Azure.

Knobs: FAKE_AZ_EXISTS=1 (resources already exist), FAKE_AZ_POOL_EXISTS=0|1, FAKE_AZ_TASK_EXISTS=1
(report the first task of each bulk add as TaskExists), FAKE_AZ_RESULTS_SRC=<dir> (tree copied
by `storage blob download-batch`), FAKE_AZ_ACR_PASS=<registry password>.
"""
import json, os, re, shutil, sys

SUB = "0123abcd-4567-89ef-0123-456789abcdef"
argv = sys.argv[1:]
env = os.environ
with open(env.get("FAKE_AZ_LOG", os.devnull), "a") as f:
    f.write(json.dumps(argv) + "\n")

path = []
for t in argv:
    if t.startswith("-"):
        break
    path.append(t)
cmd = " ".join(path)


def opt(*names, default=None):
    for i, t in enumerate(argv):
        for n in names:
            if t == n and i + 1 < len(argv):
                return argv[i + 1]
            if t.startswith(n + "="):
                return t.split("=", 1)[1]
    return default


def capture(kind):
    src = opt("--json-file")
    d = env.get("FAKE_AZ_CAPTURE")
    if src and d:
        os.makedirs(d, exist_ok=True)
        n = len(os.listdir(d))
        shutil.copy(src, os.path.join(d, f"{n:04d}_{kind}.json"))
    return json.load(open(src)) if src else None


def not_found(what):
    sys.stderr.write(f"(ResourceNotFound) The {what} was not found.\n")
    sys.exit(3)


exists = env.get("FAKE_AZ_EXISTS") == "1"
out = None
if cmd == "account show":
    out = {"id": SUB, "name": "fake-sub", "state": "Enabled"}
elif cmd in ("provider register", "group create", "storage account create", "acr create", "batch account create",
             "batch account login", "batch job create", "batch job set", "batch pool autoscale enable", "acr login",
             "storage blob upload"):
    out = {}
elif cmd == "storage container create":
    out = {"created": not exists}
elif cmd in ("storage account show", "acr show"):
    if not exists:
        not_found(cmd.split()[0])
    out = {"name": opt("-n", "--name"), "provisioningState": "Succeeded"}
elif cmd == "batch account show":
    if not exists:
        not_found("Batch account")
    out = {"name": opt("-n", "--name"), "dedicatedCoreQuota": 500, "lowPriorityCoreQuota": 100, "poolQuota": 100,
           "activeJobAndJobScheduleQuota": 300, "dedicatedCoreQuotaPerVmFamilyEnforced": True,
           "dedicatedCoreQuotaPerVmFamily": [{"name": "standardDDSv5Family", "coreQuota": 500},
                                             {"name": "standardNCADSA100v4Family", "coreQuota": 0}]}
elif cmd == "storage account keys list":
    out = [{"keyName": "key1", "value": "ZmFrZWtleQ==", "permissions": "FULL"}]
elif cmd == "storage container generate-sas":
    out = f"se=2026-10-08T00%3A00Z&sp={opt('--permissions')}&spr=https&sv=2022-11-02&sr=c&sig=a%2Bb%2Fc%3D"
elif cmd == "storage blob list":
    prefix = opt("--prefix", default="")
    n = 7 if "/mini-" in prefix else 63
    out = [{"name": f"{prefix}{i:05d}-of-{n:05d}.parquet", "properties": {"contentLength": 1}} for i in range(n)]
    out.append({"name": prefix.rsplit("/", 1)[0] + "/README.md", "properties": {}})
elif cmd == "storage blob download-batch":
    src, dest = env.get("FAKE_AZ_RESULTS_SRC"), opt("-d", "--destination")
    if src:
        shutil.copytree(src, dest, dirs_exist_ok=True)
    out = []
elif cmd == "acr credential show":
    # default stresses quoting; FAKE_AZ_ACR_PASS gives a realistic (base64-like) one
    out = {"username": "fakeacr", "passwords": [{"name": "password", "value": env.get("FAKE_AZ_ACR_PASS", "p\"a'ss+/=word")}]}
elif cmd == "acr build":
    out = {"status": "Succeeded"}
elif cmd == "batch pool show":
    if env.get("FAKE_AZ_POOL_EXISTS", "1") != "1":
        sys.stderr.write("(PoolNotFound) The specified pool does not exist.\n")
        sys.exit(1)
    out = {"id": opt("--pool-id")}
elif cmd == "batch pool create":
    capture("pool")
    out = None
elif cmd == "batch task create":
    tasks = capture("tasks")
    tasks = tasks if isinstance(tasks, list) else [tasks]
    out = [{"status": "success", "taskId": t.get("id"), "eTag": "0x1", "error": None} for t in tasks]
    if env.get("FAKE_AZ_TASK_EXISTS") == "1" and out:
        out[0] = {"status": "clientError", "taskId": out[0]["taskId"], "error": {"code": "TaskExists"}}
elif cmd == "batch job task-counts show":
    out = {"taskCounts": {"active": 3, "running": 16, "completed": 5, "succeeded": 4, "failed": 1},
           "taskSlotCounts": {"active": 3, "running": 16, "completed": 5, "succeeded": 4, "failed": 1}}
elif cmd == "batch task list":
    out = [{"id": "bundle-full-00003-of-00063", "executionInfo": {"exitCode": 2, "failureInfo": {"code": "FailureExitCode"}}}]
elif cmd == "batch location list-skus":
    out = [{"name": "Standard_D64ds_v5", "familyName": "standardDDSv5Family"},
           {"name": "Standard_NC24ads_A100_v4", "familyName": "standardNCADSA100v4Family"}]
elif cmd == "batch location quotas show":
    out = {"accountQuota": 3}
elif cmd == "vm list-usage":
    out = [{"currentValue": 0, "limit": 10, "name": {"value": "cores", "localizedValue": "Total Regional vCPUs"}, "unit": "Count"},
           {"currentValue": 0, "limit": 100, "name": {"value": "standardDDSv5Family", "localizedValue": "x"}, "unit": "Count"}]
else:
    sys.stderr.write(f"fake_az: no canned answer for '{cmd}'\n")
    sys.exit(2)

q = opt("--query")
if q and out is not None:
    try:
        import jmespath
        out = jmespath.search(q, out)
    except ImportError:
        pass
fmt = opt("-o", "--output", default="json")
if out is None or fmt == "none":
    sys.exit(0)
if fmt == "tsv":
    rows = out if isinstance(out, list) else [out]
    for r in rows:
        if isinstance(r, dict):
            print("\t".join("" if v is None else str(v) for v in r.values()))
        elif r is not None:
            print(r)
else:
    print(json.dumps(out, indent=2))
