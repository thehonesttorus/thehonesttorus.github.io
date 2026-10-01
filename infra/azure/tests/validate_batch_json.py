#!/usr/bin/env python3
"""Validate captured `--json-file` payloads against the Batch data-plane models that the installed
az CLI deserializes them into (azure-batch 15.x: BatchPoolCreateContent for pool create,
BatchTaskCreateContent for each task of a bulk add).  Those models silently keep unknown keys, so
a misspelled property would only surface as wrong behaviour once the service sees it; this walks
every key against the REST field names, checks enum values and the fields documented as Required.

    /root/azcli/bin/python validate_batch_json.py <capture-dir>      (needs az's own python)
"""
import enum, glob, json, os, re, sys
import azure.batch.models as m
from azure.batch._model_base import Model

def required_fields(cls):
    doc = cls.__doc__ or ""
    req = set()
    for name, body in re.findall(r":ivar (\w+):(.*?)(?=:vartype|\Z)", doc, re.S):
        if "Required." in body:
            req.add(name)
    return req

def enum_types(cls, attr):
    ann = str(cls.__annotations__.get(attr, ""))
    out = []
    for n in re.findall(r"_models\.(\w+)", ann):
        t = getattr(m, n, None)
        if isinstance(t, type) and issubclass(t, enum.Enum):
            out.append(t)
    return out

def check(d, cls, path, errs):
    if not isinstance(d, dict):
        errs.append(f"{path}: expected object for {cls.__name__}, got {type(d).__name__}")
        return
    inst = cls(d)
    fields = {rf._rest_name: attr for attr, rf in cls._attr_to_rest_field.items()}
    for k in d:
        if k not in fields:
            near = [f for f in fields if f.lower() == k.lower()]
            errs.append(f"{path}.{k}: not a property of {cls.__name__}" + (f" (did you mean {near[0]}?)" if near else ""))
    for attr in required_fields(cls):
        rest = cls._attr_to_rest_field[attr]._rest_name
        if d.get(rest) is None:
            errs.append(f"{path}: missing required {rest} ({cls.__name__})")
    for rest, attr in fields.items():
        if rest not in d:
            continue
        val = getattr(inst, attr)
        for et in enum_types(cls, attr):
            vals = {str(x.value).lower() for x in et}
            for v in (val if isinstance(val, list) else [val]):
                if isinstance(v, str) and v.lower() not in vals:
                    errs.append(f"{path}.{rest}: {v!r} not in {et.__name__} {sorted(vals)}")
        items = val if isinstance(val, list) else [val]
        raw = d[rest] if isinstance(d[rest], list) else [d[rest]]
        for i, (v, r) in enumerate(zip(items, raw)):
            if isinstance(v, Model):
                check(r, type(v), f"{path}.{rest}" + (f"[{i}]" if isinstance(val, list) else ""), errs)

errs, n = [], 0
for fn in sorted(glob.glob(os.path.join(sys.argv[1], "*.json"))):
    data = json.load(open(fn))
    if fn.endswith("_pool.json"):
        check(data, m.BatchPoolCreateContent, os.path.basename(fn), errs); n += 1
    else:
        for i, t in enumerate(data if isinstance(data, list) else [data]):
            check(t, m.BatchTaskCreateContent, f"{os.path.basename(fn)}[{i}]", errs); n += 1
for e in errs:
    print("FAIL", e)
print(f"validated {n} objects, {len(errs)} problems")
sys.exit(1 if errs else 0)
