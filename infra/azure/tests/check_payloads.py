#!/usr/bin/env python3
"""Offline checks of what the scripts would send to Batch (captured by fake_az.py):

1. every task commandLine is `/bin/bash -c '<cmd>'`, and <cmd> (shell-split exactly as bash
   would) is accepted by the argparse of the runner script it calls (run_shard.py,
   stage_dataset.py, bake_shard.py, bake_moments.py) -- the runner's main() is entered with that
   argv and stopped right after parse_args;
2. task ids are unique per job file and <= 64 chars of [A-Za-z0-9_-];
3. container tasks do not pass --gpus (Batch adds GPUs itself on GPU pools) and run as admin
   (Batch maps the task user into the container);
4. every autoscale formula only uses the documented operators/functions/variables (no %),
   ends statements with ';', and assigns a target variable.

    python3 check_payloads.py <capture-dir> <runner-dir>
"""
import argparse, glob, json, os, re, runpy, shlex, sys, types

cap, runner = sys.argv[1], sys.argv[2]
errs, n_tasks, n_pools = [], 0, 0

class Parsed(Exception):
    pass

def parse_with_runner(script, argv):
    orig = argparse.ArgumentParser.parse_args
    def fake_parse(self, args=None, namespace=None):
        raise Parsed(orig(self, args, namespace))
    argparse.ArgumentParser.parse_args = fake_parse
    saved_argv, stubs = sys.argv, []
    for stub in ("torch", "huggingface_hub", "azure.storage.blob"):   # heavy/absent imports, unused before parsing
        if stub not in sys.modules:
            sys.modules[stub] = types.ModuleType(stub)
            stubs.append(stub)
    for attr in ("HfApi", "hf_hub_download"):
        if not hasattr(sys.modules["huggingface_hub"], attr):
            setattr(sys.modules["huggingface_hub"], attr, None)
    try:
        sys.argv = [script] + argv
        runpy.run_path(os.path.join(runner, os.path.basename(script)), run_name="__main__")
        return None, "main() returned without parsing"
    except Parsed as p:
        return p.args[0], None
    except SystemExit as e:
        return None, f"argparse rejected: exit {e.code}"
    finally:
        argparse.ArgumentParser.parse_args = orig
        sys.argv = saved_argv
        for k in stubs:          # only the stubs; real modules (numpy) cannot be re-imported
            del sys.modules[k]

FUNCS = {"avg", "ceil", "floor", "len", "lg", "ln", "log", "max", "min", "norm", "percentile", "rand", "range",
         "round", "std", "stop", "sum", "time", "val"}
METHODS = {"GetSample", "GetSamplePeriod", "Count", "HistoryBeginTime", "GetSamplePercent"}
SVARS = {"$TargetDedicatedNodes", "$TargetLowPriorityNodes", "$NodeDeallocationOption", "$TargetDedicated",
         "$TargetLowPriority", "$CPUPercent", "$ActiveTasks", "$RunningTasks", "$PendingTasks", "$SucceededTasks",
         "$FailedTasks", "$TaskSlotsPerNode", "$CurrentDedicatedNodes", "$CurrentLowPriorityNodes", "$UsableNodeCount",
         "$PreemptedNodeCount"}
DEALLOC = {"requeue", "terminate", "taskcompletion", "retaineddata"}

def lint_formula(f, where):
    if len(f) > 8000:
        errs.append(f"{where}: formula longer than 8 KB")
    if "%" in f:
        errs.append(f"{where}: '%' is not an autoscale operator")
    stmts = [s.strip() for s in f.split(";")]
    if stmts[-1]:
        errs.append(f"{where}: last statement not terminated by ';'")
    for v in re.findall(r"\$\w+", f):
        if v not in SVARS:
            errs.append(f"{where}: unknown service variable {v}")
    for fn in re.findall(r"(?<![.\w$])([A-Za-z_]\w*)\s*\(", f):
        if fn not in FUNCS:
            errs.append(f"{where}: {fn}() is not an autoscale function")
    for meth in re.findall(r"\.(\w+)\s*\(", f):
        if meth not in METHODS:
            errs.append(f"{where}: .{meth}() is not a sample method")
    for opnd in re.findall(r"\$NodeDeallocationOption\s*=\s*(\w+)", f):
        if opnd not in DEALLOC:
            errs.append(f"{where}: bad $NodeDeallocationOption {opnd}")
    if not re.search(r"\$Target(Dedicated|LowPriority)(Nodes)?\s*=", f):
        errs.append(f"{where}: no target node variable assigned")
    if f.count("(") != f.count(")"):
        errs.append(f"{where}: unbalanced parentheses")

for fn in sorted(glob.glob(os.path.join(cap, "*.json"))):
    data = json.load(open(fn))
    name = os.path.basename(fn)
    if name.endswith("_pool.json"):
        n_pools += 1
        lint_formula(data.get("autoScaleFormula", ""), name)
        continue
    tasks = data if isinstance(data, list) else [data]
    ids = [t.get("id", "") for t in tasks]
    if len(set(ids)) != len(ids):
        errs.append(f"{name}: duplicate task ids")
    for t in tasks:
        n_tasks += 1
        tid = t.get("id", "")
        if not re.fullmatch(r"[A-Za-z0-9_-]{1,64}", tid):
            errs.append(f"{name}: invalid task id {tid!r}")
        cs = t.get("containerSettings") or {}
        if "--gpus" in cs.get("containerRunOptions", ""):
            errs.append(f"{name}/{tid}: containerRunOptions passes --gpus")
        au = (t.get("userIdentity") or {}).get("autoUser") or {}
        if au.get("elevationLevel") != "admin":
            errs.append(f"{name}/{tid}: container task not admin (cannot write outside the task dir)")
        cl = t.get("commandLine", "")
        outer = shlex.split(cl)
        if outer[:2] != ["/bin/bash", "-c"] or len(outer) != 3:
            errs.append(f"{name}/{tid}: commandLine is not /bin/bash -c '<cmd>': {cl[:80]}")
            continue
        inner = shlex.split(outer[2])
        if inner[0] != "python" or not inner[1].startswith("/app/"):
            errs.append(f"{name}/{tid}: unexpected command {inner[:2]}")
            continue
        ns, err = parse_with_runner(inner[1], inner[2:])
        if err:
            errs.append(f"{name}/{tid}: {inner[1]} {err}")
        elif os.path.basename(inner[1]) == "run_shard.py" and not (ns.out.startswith("https://") and "sig=" in ns.out):
            errs.append(f"{name}/{tid}: --out is not a SAS URL")

for e in errs:
    print("FAIL", e)
print(f"checked {n_tasks} tasks and {n_pools} pool formulas, {len(errs)} problems")
sys.exit(1 if errs else 0)
