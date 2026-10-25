"""Build the Kaggle submission notebook for the ARC Prize 2026 - ARC-AGI-3 competition.

Mirrors the official starter kit's notebook (arcprize/ARC-AGI-3-Kaggle-Starter,
scripts/build_notebook.py), with one change: the whole ``nce`` package is shipped,
not a single agent file. The package is zipped, base64-embedded in a notebook cell,
and extracted to /kaggle/working/nce_pkg at run time (no internet is needed).

Usage:  python tools/build_kaggle_notebook.py [--accelerator cpu|t4|p100|rtx6000] [--out notebooks/submission.ipynb]
"""
from __future__ import annotations

import argparse
import base64
import io
import json
import zipfile
from pathlib import Path
from textwrap import dedent

ROOT = Path(__file__).resolve().parents[1]
PKG_DIR = ROOT / "nce"
AGENT_SRC = ROOT / "kaggle" / "my_agent.py"

_ACCELERATORS = {
    "cpu": {"name": "none", "gpu": False},
    "t4": {"name": "nvidiaTeslaT4", "gpu": True},
    "p100": {"name": "nvidiaTeslaP100", "gpu": True},
    "rtx6000": {"name": "nvidiaRtx6000", "gpu": True},
}


def code_cell(source: str) -> dict:
    return {"cell_type": "code", "metadata": {"trusted": True}, "outputs": [], "execution_count": None, "source": source}


def markdown_cell(source: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": source}


def zip_package_b64() -> str:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(PKG_DIR.rglob("*.py")):
            zf.write(p, arcname=str(Path("nce") / p.relative_to(PKG_DIR)))
    return base64.b64encode(buf.getvalue()).decode("ascii")


def build(accelerator: str) -> dict:
    if accelerator not in _ACCELERATORS:
        raise SystemExit(f"unknown accelerator {accelerator!r}")
    accel = _ACCELERATORS[accelerator]
    agent_body = AGENT_SRC.read_text()
    pkg_b64 = zip_package_b64()

    install_cell = code_cell(
        "!pip install --no-index --find-links \\\n"
        "    /kaggle/input/competitions/arc-prize-2026-arc-agi-3/arc_agi_3_wheels \\\n"
        "    arc-agi python-dotenv"
    )
    unpack_cell = code_cell(dedent(f'''\
        import base64, io, os, time, zipfile
        # A global deadline for every game (the notebook's 9-hour limit, minus a margin for submission writing):
        # games may run in batches, so a per-game wall clock alone cannot guarantee the notebook finishes.
        os.environ.setdefault("NCE_DEADLINE", str(time.time() + 8.25 * 3600))
        # 48 vCPUs for the hidden games: several seconds of compute per action are available, so every engine budget
        # (synthesis, planning, Kato) is raised together
        os.environ.setdefault("NCE_COMPUTE_SCALE", "4")
        _B64 = """{pkg_b64}"""
        os.makedirs("/kaggle/working/nce_pkg", exist_ok=True)
        with zipfile.ZipFile(io.BytesIO(base64.b64decode(_B64))) as zf:
            zf.extractall("/kaggle/working/nce_pkg")
        os.environ["NCE_PKG_DIR"] = "/kaggle/working/nce_pkg"
        print("nce package files:", sum(len(f) for _, _, f in os.walk("/kaggle/working/nce_pkg")))
        '''))
    write_agent_cell = code_cell("%%writefile /tmp/my_agent.py\n" + agent_body)
    run_cell = code_cell(dedent('''\
        import os

        if os.getenv('KAGGLE_IS_COMPETITION_RERUN'):
            !curl --fail --retry 999 --retry-all-errors --retry-delay 5 \\
                  --retry-max-time 600 http://gateway:8001/api/games
            !cp -r /kaggle/input/competitions/arc-prize-2026-arc-agi-3/ARC-AGI-3-Agents \\
                   /kaggle/working/ARC-AGI-3-Agents
            !cp /tmp/my_agent.py \\
                /kaggle/working/ARC-AGI-3-Agents/agents/templates/my_agent.py
            with open('/kaggle/working/ARC-AGI-3-Agents/agents/__init__.py', 'w') as f:
                f.write("""from typing import Type
        from dotenv import load_dotenv
        from .agent import Agent, Playback
        from .swarm import Swarm
        from .templates.random_agent import Random
        from .templates.my_agent import MyAgent

        load_dotenv()

        AVAILABLE_AGENTS: dict[str, Type[Agent]] = {
            'random': Random,
            'myagent': MyAgent,
        }
        """)
            with open('/kaggle/working/ARC-AGI-3-Agents/.env', 'w') as f:
                f.write("""SCHEME=http
        HOST=gateway
        PORT=8001
        ARC_API_KEY=test-key-123
        ARC_BASE_URL=http://gateway:8001/
        OPERATION_MODE=online
        ENVIRONMENTS_DIR=
        RECORDINGS_DIR=/kaggle/working/server_recording
        """)
            !cd /kaggle/working/ARC-AGI-3-Agents && \\
                NCE_PKG_DIR=/kaggle/working/nce_pkg MPLBACKEND=agg \\
                python main.py --agent myagent
        '''))
    dummy_cell = code_cell(dedent('''\
        import os
        if not os.getenv('KAGGLE_IS_COMPETITION_RERUN'):
            import pandas as pd
            submission = pd.DataFrame(data=[['1_0', '1', True, 1]],
                                      columns=['row_id', 'game_id', 'end_of_game', 'score'])
            submission.to_parquet('/kaggle/working/submission.parquet', index=False)
            submission.head()
        '''))
    return {
        "metadata": {
            "kernelspec": {"language": "python", "display_name": "Python 3", "name": "python3"},
            "language_info": {"name": "python", "mimetype": "text/x-python", "file_extension": ".py", "pygments_lexer": "ipython3"},
            "kaggle": {"accelerator": accel["name"], "isInternetEnabled": False, "isGpuEnabled": accel["gpu"], "language": "python", "sourceType": "notebook"},
        },
        "nbformat_minor": 4,
        "nbformat": 4,
        "cells": [
            markdown_cell("# ARC Prize 2026 — ARC-AGI-3 — Noncommutative Correspondence Engine\n\nBuilt by `tools/build_kaggle_notebook.py`; edit the `nce` package and `kaggle/my_agent.py`, not this notebook."),
            install_cell, unpack_cell, write_agent_cell, run_cell, dummy_cell,
        ],
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--accelerator", default="cpu", choices=sorted(_ACCELERATORS))
    p.add_argument("--out", default=str(ROOT / "notebooks" / "submission.ipynb"))
    a = p.parse_args()
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(build(a.accelerator), indent=1))
    print(f"wrote {out} (accelerator={a.accelerator}, package bytes b64={len(zip_package_b64())})")


if __name__ == "__main__":
    main()
