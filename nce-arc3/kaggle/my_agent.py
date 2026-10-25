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
    # Absolute epoch deadline shared by all games (set by the submission notebook); 0 disables it.
    DEADLINE = float(os.environ.get("NCE_DEADLINE", "0"))
    # Wall-clock backstop per game, seconds; only used when no global deadline is set (a per-game limit would stop
    # a progressing game while notebook time remains).
    WALL_CLOCK_LIMIT = float(os.environ.get("NCE_WALL_CLOCK", "0" if DEADLINE else str(4 * 3600)))

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self._engine: Optional[Engine] = None
        self._last_action: Optional[tuple[int, Optional[dict[str, Any]]]] = None
        self._last_reason = ""
        self._step_failed = False
        self._t0 = time.time()

    @property
    def name(self) -> str:
        return f"{super().name}.nce"

    def _ensure_engine(self, obs: dict[str, Any]) -> Optional[Engine]:
        """Create the engine from the first real frame (never from the framework's placeholder, which carries
        win_levels = 0 and no available actions)."""
        if self._engine is None and obs["state"] != "NOT_PLAYED" and obs["frames"] and obs["available_actions"]:
            self._engine = Engine(self.game_id, obs["win_levels"], obs["available_actions"])
        return self._engine

    def take_action(self, action: GameAction) -> Optional[FrameData]:  # type: ignore[override]
        """Submit our own (action id, coordinates). GameAction members are process-wide singletons shared by every
        game thread, so their action_data / reasoning fields are never written or read here. A failed or invalid
        step returns None (the game continues) and the next choose_action resynchronises the engine."""
        aid, data = self._last_action if self._last_action is not None else (int(action.value), None)
        act = GameAction.from_id(int(aid))
        payload = None
        if act.is_complex():
            d = data or {"x": 0, "y": 0}
            payload = {"x": int(min(63, max(0, int(d["x"])))), "y": int(min(63, max(0, int(d["y"]))))}
        try:
            raw = self.arc_env.step(act, data=payload, reasoning={"nce": self._last_reason[:200]})
            if raw is None:
                raise ValueError("no frame returned")
            return FrameData.model_validate(self._convert_raw_frame_data(raw))
        except Exception as e:
            print(f"[nce] step failed on {self.game_id}: {e!r}", file=sys.stderr)
            self._step_failed = True
            return None

    def cleanup(self, scorecard=None) -> None:  # type: ignore[override]
        try:
            if self._engine is not None and hasattr(self._engine, "close"):
                self._engine.close()
        finally:
            super().cleanup(scorecard)

    def is_done(self, frames: list[FrameData], latest_frame: FrameData) -> bool:
        if latest_frame.state is GameState.WIN:
            return True
        if self.WALL_CLOCK_LIMIT and time.time() - self._t0 > self.WALL_CLOCK_LIMIT:
            return True
        if self.DEADLINE and time.time() > self.DEADLINE:
            return True
        if self._engine is not None and hasattr(self._engine, "want_stop"):
            try:
                return bool(self._engine.want_stop(_obs_from_frame(latest_frame)))
            except Exception:
                return False
        return False

    def choose_action(self, frames: list[FrameData], latest_frame: FrameData) -> GameAction:
        if len(frames) > 4:
            del frames[1:-2]              # the framework keeps every frame; only the latest few are ever read
        obs = _obs_from_frame(latest_frame)
        engine = self._ensure_engine(obs)

        # Deliver the outcome of our previous action to the engine; after a failed step the outcome is unknown,
        # so the engine adopts the current state without storing a demonstration.
        if engine is not None and self._last_action is not None:
            aid, data = self._last_action
            try:
                if self._step_failed:
                    engine.resync(obs)
                else:
                    engine.observe(obs, aid, data)
            except Exception:
                pass
        self._step_failed = False

        if engine is None or obs["state"] in ("NOT_PLAYED", "GAME_OVER"):
            if engine is not None and obs["state"] == "GAME_OVER" and hasattr(engine, "on_game_over"):
                try:
                    engine.on_game_over(obs)
                except Exception:
                    pass
            aid, data = 0, None
        else:
            try:
                aid, data = engine.act(obs)
            except Exception as e:  # never crash the swarm
                avail = [a for a in obs["available_actions"] if a != 0] or [1]
                aid, data = avail[0], ({"x": 0, "y": 0} if avail[0] == 6 else None)
                print(f"[nce] engine error on {self.game_id}: {e!r}", file=sys.stderr)

        if engine is not None and obs["available_actions"] and aid != 0 and int(aid) not in obs["available_actions"]:
            aid = [a for a in obs["available_actions"] if a != 0][0] if [a for a in obs["available_actions"] if a != 0] else 0
            data = {"x": 0, "y": 0} if aid == 6 else None
        self._last_reason = str(getattr(engine, "last_reason", "")) if engine is not None else "reset"
        self._last_action = (int(aid), data)
        return GameAction.from_id(int(aid))     # take_action reads self._last_action, never the shared member
