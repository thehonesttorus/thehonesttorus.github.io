"""Run kaggle/my_agent.py through the *official* ARC-AGI-3-Agents framework loop locally.

This is the deployment-faithful path (same Agent.main() loop Kaggle's gateway drives),
as opposed to tools/harness.py which talks to the engine directly and is faster.

Requires a clone of https://github.com/arcprize/ARC-AGI-3-Agents at ./vendor/ARC-AGI-3-Agents
(or set ARC_AGENTS_DIR). Game sources are downloaded on first use into ./environment_files
via the arc-agi package (anonymous key is fine) and cached.

Usage: python tools/play_local_framework.py [--games ls20,vc33] [--max-actions 500]
"""
from __future__ import annotations

import argparse
import importlib.util
import logging
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VENDOR = Path(os.environ.get("ARC_AGENTS_DIR", ROOT / "vendor" / "ARC-AGI-3-Agents"))
if not VENDOR.exists():
    raise SystemExit(f"framework not found at {VENDOR}; clone ARC-AGI-3-Agents there or set ARC_AGENTS_DIR")
sys.path.insert(0, str(VENDOR))
sys.path.insert(0, str(ROOT))
os.environ.setdefault("NCE_PKG_DIR", str(ROOT))

import arc_agi  # noqa: E402
from arc_agi import OperationMode  # noqa: E402


def load_my_agent():
    spec = importlib.util.spec_from_file_location("nce_kaggle_agent", ROOT / "kaggle" / "my_agent.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.MyAgent


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--games", default=None)
    p.add_argument("--max-actions", type=int, default=500)
    a = p.parse_args()
    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    arc = arc_agi.Arcade(operation_mode=OperationMode.NORMAL)
    envs = arc.get_environments()
    gids = [e.game_id.split("-")[0] for e in envs]
    if a.games:
        gids = [g for g in a.games.split(",") if g in gids]
    MyAgent = load_my_agent()
    MyAgent.MAX_ACTIONS = a.max_actions
    rows = []
    for gid in gids:
        env = arc.make(gid)
        agent = MyAgent(card_id="local", game_id=gid, agent_name=f"nce.{gid}", ROOT_URL="http://localhost", record=False, arc_env=env, tags=["local"])
        agent.main()
        final = agent.frames[-1]
        rows.append((gid, final.levels_completed, agent.action_counter, final.state.name))
        print(f"{gid}: levels={final.levels_completed} actions={agent.action_counter} state={final.state.name}", flush=True)
    sc = arc.get_scorecard()
    print("scorecard:", getattr(sc, "score", sc))


if __name__ == "__main__":
    main()
