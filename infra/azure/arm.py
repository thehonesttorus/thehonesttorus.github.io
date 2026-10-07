#!/usr/bin/env python3
"""Minimal Azure Resource Manager client for experiment infrastructure.

Reads AZURE_TENANT_ID, AZURE_CLIENT_ID, AZURE_CLIENT_SECRET, AZURE_SUBSCRIPTION_ID from the environment
(never from arguments or files) and never prints the token.

  python3 infra/azure/arm.py get    <path>                 e.g. /resourcegroups?api-version=2021-04-01
  python3 infra/azure/arm.py put    <path> <json|@file>
  python3 infra/azure/arm.py post   <path> [json|@file]
  python3 infra/azure/arm.py delete <path>
  python3 infra/azure/arm.py wait   <azure-asyncoperation-url>   poll until Succeeded/Failed

A path starting with "/" and not "/subscriptions" is prefixed with /subscriptions/$AZURE_SUBSCRIPTION_ID.
Mutating calls (put/post/delete) are refused unless the path is inside a resource group named rg-claude-*,
or is a Microsoft.Quota request, or is a VM run-command/power action inside rg-claude-*.
"""
import json, os, sys, time, urllib.request, urllib.parse, urllib.error

ARM = "https://management.azure.com"


def token():
    body = urllib.parse.urlencode({
        "client_id": os.environ["AZURE_CLIENT_ID"], "client_secret": os.environ["AZURE_CLIENT_SECRET"],
        "scope": "https://management.azure.com/.default", "grant_type": "client_credentials"}).encode()
    url = f"https://login.microsoftonline.com/{os.environ['AZURE_TENANT_ID']}/oauth2/v2.0/token"
    with urllib.request.urlopen(urllib.request.Request(url, data=body)) as r:
        return json.load(r)["access_token"]


def full(path):
    if path.startswith("https://"):
        return path
    if not path.startswith("/subscriptions"):
        path = f"/subscriptions/{os.environ['AZURE_SUBSCRIPTION_ID']}{path}"
    return ARM + path


def allowed_mutation(url):
    p = urllib.parse.urlparse(url).path.lower()
    if "/providers/microsoft.quota/" in p:
        return True
    parts = p.split("/")
    if "resourcegroups" in parts:
        i = parts.index("resourcegroups")
        return i + 1 < len(parts) and parts[i + 1].startswith("rg-claude-")
    return False


def call(method, path, data=None):
    url = full(path)
    if method != "GET" and not allowed_mutation(url):
        sys.exit(f"refused: {method} outside rg-claude-* resource groups: {url}")
    if isinstance(data, str) and data.startswith("@"):
        data = open(data[1:]).read()
    req = urllib.request.Request(url, method=method, data=(data.encode() if isinstance(data, str) else None),
                                 headers={"Authorization": f"Bearer {token()}", "Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req) as r:
            raw = r.read().decode(); hdr = dict(r.headers); code = r.status
    except urllib.error.HTTPError as e:
        raw = e.read().decode(); hdr = dict(e.headers); code = e.code
    out = {"http": code}
    for h in ("Azure-AsyncOperation", "Location"):
        if h in hdr:
            out[h] = hdr[h]
    try:
        out["body"] = json.loads(raw) if raw else None
    except json.JSONDecodeError:
        out["body"] = raw[:2000]
    return out


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    verb, path = sys.argv[1].lower(), sys.argv[2]
    data = sys.argv[3] if len(sys.argv) > 3 else None
    if verb == "wait":
        while True:
            r = call("GET", path); st = (r.get("body") or {}).get("status")
            if st in ("Succeeded", "Failed", "Canceled"):
                print(json.dumps(r, indent=1)); return
            time.sleep(10)
    print(json.dumps(call(verb.upper(), path, data), indent=1))


if __name__ == "__main__":
    main()
