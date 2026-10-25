"""Kaggle / ARC-AGI-3-Agents adapter for the Noncommutative Correspondence Engine.

Contract enforced by the ARC-AGI-3-Agents framework (used verbatim by the Kaggle
competition rerun):
  - subclass ``agents.agent.Agent``,
  - class must be named ``MyAgent``,
  - implement ``is_done(frames, latest_frame) -> bool`` and
    ``choose_action(frames, latest_frame) -> GameAction``.

The framework calls ``choose_action`` with the latest observation and then
``take_action``; we observe the resulting frame on the *next* ``choose_action``
call, so this adapter keeps the last action and feeds (obs, action) pairs to the
engine's ``observe`` hook before asking it to act again.

Engine protocol (see nce/engine.py and tools/harness.py):
    engine = Engine(game_id, win_levels, available_actions)
    action_id, data = engine.act(obs)              # data = {"x":..,"y":..} for ACTION6 else None
    engine.observe(obs, action_id, data)           # optional
    engine.on_game_over(obs)                       # optional
    engine.want_stop(obs) -> bool                  # optional: engine's own budget policy
where obs = dict(frames=[np.int8 (64,64), ...], state=str, levels_completed=int,
                 win_levels=int, available_actions=list[int], full_reset=bool).
"""
from __future__ import annotations

import os
import sys
import time
from typing import Any, Optional

import numpy as np
from arcengine import FrameData, GameAction, GameState

from agents.agent import Agent

# The engine package is shipped next to this file (locally) or extracted to a
# known location by the submission notebook (on Kaggle); both are tried.
for _cand in (
    os.environ.get("NCE_PKG_DIR"),
    os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."),
    "/kaggle/working/nce_pkg",
    "/tmp/nce_pkg",
):
    if _cand and os.path.isdir(os.path.join(_cand, "nce")) and _cand not in sys.path:
        sys.path.insert(0, _cand)

from nce.proxy import EngineProxy as Engine  # noqa: E402  (process-per-game; NCE_USE_PROCESS=0 for in-process)


def _obs_from_frame(frame: FrameData) -> dict[str, Any]:
    return dict(
        frames=[np.asarray(f, dtype=np.int8) for f in frame.frame],
        state=frame.state.name if hasattr(frame.state, "name") else str(frame.state),
        levels_completed=int(frame.levels_completed),
        win_levels=int(frame.win_levels),
        available_actions=list(frame.available_actions or []),
        full_reset=bool(frame.full_reset),
    )


class MyAgent(Agent):
    """Framework adapter: delegates every decision to the correspondence engine."""

    # The framework stops the loop once action_counter exceeds MAX_ACTIONS; the
    # engine owns the real budget policy (RHAE-aware), so this is only a backstop.
    MAX_ACTIONS = int(os.environ.get("NCE_MAX_ACTIONS", "100000"))
    # Wall-clock backstop per game, seconds (Kaggle runs all hidden games in threads).
    WALL_CLOCK_LIMIT = float(os.environ.get("NCE_WALL_CLOCK", str(4 * 3600)))

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self._engine: Optional[Engine] = None
        self._last_action: Optional[tuple[int, Optional[dict[str, Any]]]] = None
        self._t0 = time.time()

    @property
    def name(self) -> str:
        return f"{super().name}.nce"

    def _ensure_engine(self, obs: dict[str, Any]) -> Engine:
        if self._engine is None:
            self._engine = Engine(self.game_id, obs["win_levels"], obs["available_actions"])
        return self._engine

    def cleanup(self, scorecard=None) -> None:  # type: ignore[override]
        try:
            if self._engine is not None and hasattr(self._engine, "close"):
                self._engine.close()
        finally:
            super().cleanup(scorecard)

    def is_done(self, frames: list[FrameData], latest_frame: FrameData) -> bool:
        if latest_frame.state is GameState.WIN:
            return True
        if time.time() - self._t0 > self.WALL_CLOCK_LIMIT:
            return True
        if self._engine is not None and hasattr(self._engine, "want_stop"):
            try:
                return bool(self._engine.want_stop(_obs_from_frame(latest_frame)))
            except Exception:
                return False
        return False

    def choose_action(self, frames: list[FrameData], latest_frame: FrameData) -> GameAction:
        obs = _obs_from_frame(latest_frame)
        engine = self._ensure_engine(obs)

        # Deliver the outcome of our previous action to the engine.
        if self._last_action is not None and hasattr(engine, "observe"):
            aid, data = self._last_action
            try:
                engine.observe(obs, aid, data)
            except Exception:
                pass

        if obs["state"] in ("NOT_PLAYED", "GAME_OVER"):
            if obs["state"] == "GAME_OVER" and hasattr(engine, "on_game_over"):
                try:
                    engine.on_game_over(obs)
                except Exception:
                    pass
            aid, data = 0, None
        else:
            try:
                aid, data = engine.act(obs)
            except Exception as e:  # never crash the swarm: fall back to RESET-free safe action
                avail = [a for a in obs["available_actions"] if a != 0] or [1]
                aid, data = avail[0], ({"x": 0, "y": 0} if avail[0] == 6 else None)
                print(f"[nce] engine error on {self.game_id}: {e!r}", file=sys.stderr)

        action = GameAction.from_id(int(aid))
        if action.is_complex():
            d = data or {"x": 0, "y": 0}
            action.set_data({"x": int(d["x"]), "y": int(d["y"])})
        action.reasoning = {"nce": getattr(engine, "last_reason", "")}
        self._last_action = (int(aid), data)
        return action
