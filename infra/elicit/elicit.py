#!/usr/bin/env python3
"""Elicit API v2 client (Research Agent sessions, reports, paper search, library full text).

The key is read from ~/.elicit_hdr (a 600-mode file holding the Authorization header line) and never printed.
Requests go through curl, which the agent proxy is configured for (Python's default client is refused at the edge).
Polling stays far below the shared limit of 100 requests per minute per IP.

  python3 infra/elicit/elicit.py usage
  python3 infra/elicit/elicit.py agent "QUERY" [--name NAME]      start a Research Agent session; prints its id
  python3 infra/elicit/elicit.py say SESSION "MESSAGE"            follow-up turn in a session
  python3 infra/elicit/elicit.py status [SESSION...]              status of named / all saved sessions
  python3 infra/elicit/elicit.py events SESSION [--all]           new events since the last call (text entries only)
  python3 infra/elicit/elicit.py dump SESSION                     full reduced history + sources + artifacts -> OUT/
  python3 infra/elicit/elicit.py wait SESSION... [--every 60]     block until every session is idle, then dump them
  python3 infra/elicit/elicit.py report "QUESTION" [--papers N] [--extract K]
  python3 infra/elicit/elicit.py search "QUERY" [--n N] [--since YEAR]
  python3 infra/elicit/elicit.py fulltext SOURCE_ID
"""
import json, os, subprocess, sys, time

API = "https://elicit.com/api/v2"
HDR = os.path.expanduser("~/.elicit_hdr")
OUT = os.environ.get("ELICIT_OUT", "/tmp/claude-0/-home-user-thehonesttorus-github-io/"
                     "2cab30f0-c454-5769-8d9a-8195cd80a5f3/scratchpad/elicit")
STATE = os.path.join(OUT, "sessions.json")


def call(method, path, body=None, timeout=60):
    cmd = ["curl", "-sS", "-m", str(timeout), "-X", method, "-H", f"@{HDR}", "-H", "Content-Type: application/json",
           "-w", "\n%{http_code}", API + path]
    if body is not None:
        cmd += ["--data-binary", "@-"]
    r = subprocess.run(cmd, input=(json.dumps(body) if body is not None else None), capture_output=True, text=True)
    txt, _, code = r.stdout.rpartition("\n")
    if not code.isdigit() or int(code) >= 400:
        raise RuntimeError(f"{method} {path}: HTTP {code} {txt[:400]} {r.stderr[:200]}")
    return json.loads(txt) if txt.strip() else {}


def _state():
    os.makedirs(OUT, exist_ok=True)
    return json.load(open(STATE)) if os.path.exists(STATE) else {}


def _save(st):
    json.dump(st, open(STATE, "w"), indent=1)


def _text(ev):
    """Flatten one reduced event to a line of text (keeps message / thinking / tool entries readable)."""
    t = ev.get("type", "?")
    for k in ("text", "content", "message", "summary", "title", "question"):
        v = ev.get(k)
        if isinstance(v, str) and v.strip():
            return f"[{t}] {v}"
        if isinstance(v, list):
            s = " ".join(x.get("text", "") if isinstance(x, dict) else str(x) for x in v)
            if s.strip():
                return f"[{t}] {s}"
    return f"[{t}] " + json.dumps({k: v for k, v in ev.items() if k not in ("eventId",)})[:600]


def events(sid, all_=False):
    st = _state(); cur = None if all_ else st.get(sid, {}).get("cursor")
    d = call("GET", f"/sessions/agents/{sid}/events" + (f"?cursor={cur}" if cur else ""))
    st.setdefault(sid, {})["cursor"] = d.get("cursor"); _save(st)
    return d


def dump(sid):
    d = call("GET", f"/sessions/agents/{sid}/events")
    name = _state().get(sid, {}).get("name", sid[:8])
    base = os.path.join(OUT, name); os.makedirs(base, exist_ok=True)
    json.dump(d, open(os.path.join(base, "events.json"), "w"), indent=1)
    lines = [_text(e) for e in d.get("events", [])]
    open(os.path.join(base, "transcript.md"), "w").write("\n\n".join(lines) + "\n")
    for kind in ("sources", "artifacts"):
        try:
            x = call("GET", f"/sessions/agents/{sid}/{kind}")
            json.dump(x, open(os.path.join(base, f"{kind}.json"), "w"), indent=1)
        except Exception as e:
            print(f"{kind}: {e}")
    try:
        arts = json.load(open(os.path.join(base, "artifacts.json")))
        for a in arts.get("deliveredOutputs", []) or []:
            aid = a.get("artifactId") or a.get("id")
            if aid:
                c = call("GET", f"/sessions/agents/{sid}/artifacts/{aid}/content")
                json.dump(c, open(os.path.join(base, f"output_{aid}.json"), "w"), indent=1)
        for a in arts.get("artifacts", []) or []:
            aid = a.get("artifactId") or a.get("id")
            if aid:
                u = call("GET", f"/sessions/agents/{sid}/artifacts/{aid}/download").get("url")
                fn = os.path.basename(a.get("name") or a.get("filename") or aid)
                if u:
                    subprocess.run(["curl", "-sS", "-m", "120", "-o", os.path.join(base, fn), u], check=False)
    except Exception as e:
        print(f"artifacts: {e}")
    print(f"{name}: {len(lines)} events -> {base}")


def status(sids):
    st = _state(); sids = sids or list(st)
    for s in sids:
        d = call("GET", f"/sessions/agents/{s}")
        print(f"{st.get(s, {}).get('name', '?'):14s} {s} {d.get('status')}  {d.get('title', '')}")


def main(a):
    if not a:
        sys.exit(__doc__)
    c = a[0]; opt = lambda k, dflt=None: a[a.index(k) + 1] if k in a else dflt
    if c == "usage":
        print(call("GET", "/usage"))
    elif c == "agent":
        d = call("POST", "/sessions/agents", {"query": a[1]})
        st = _state(); st[d["sessionId"]] = {"name": opt("--name", d["sessionId"][:8]), "query": a[1][:300]}; _save(st)
        print(d["sessionId"], d.get("status"), d.get("url"))
    elif c == "say":
        print(call("POST", f"/sessions/agents/{a[1]}/messages", {"message": a[2]}))
    elif c == "status":
        status(a[1:])
    elif c == "events":
        d = events(a[1], "--all" in a)
        print("status:", d.get("status"))
        for e in d.get("events", []):
            print(_text(e)[:4000])
    elif c == "dump":
        dump(a[1])
    elif c == "wait":
        sids = [x for x in a[1:] if not x.startswith("--") and x != opt("--every")]
        every = int(opt("--every", "60")); left = set(sids)
        while left:
            for s in sorted(left):
                d = call("GET", f"/sessions/agents/{s}")
                if d.get("status") != "processing":
                    print(f"{s}: {d.get('status')}", flush=True); dump(s); left.discard(s)
            if left:
                time.sleep(every)
    elif c == "report":
        d = call("POST", "/sessions/reports", {"researchQuestion": a[1], "maxSearchPapers": int(opt("--papers", "300")),
                                                "maxExtractPapers": int(opt("--extract", "40"))})
        print(json.dumps(d, indent=1))
    elif c == "search":
        body = {"query": a[1], "maxResults": int(opt("--n", "50"))}
        if opt("--since"):
            body["filters"] = {"minYear": int(opt("--since"))}
        d = call("POST", "/search/papers", body)
        for p in d.get("papers", d.get("results", [])):
            print(f"- {p.get('title')} ({p.get('year')}) {p.get('doi') or p.get('url') or ''}")
    elif c == "fulltext":
        print(call("GET", f"/library/sources/{a[1]}/full-text").get("markdown", "")[:200000])
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
