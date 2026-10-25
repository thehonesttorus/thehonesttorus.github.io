"""Observation unpacking: the only interface-level facts the engine relies on.

An observation is a dict (see tools/harness.py and kaggle/my_agent.py):
    frames: list of int8 arrays of shape (64, 64), values 0..15; several when the
            environment animates in response to one action. The last one is the
            state the next action acts on.
    state: "NOT_PLAYED" | "NOT_FINISHED" | "WIN" | "GAME_OVER"
    levels_completed, win_levels: ints
    available_actions: list of action ids (0 = RESET, 1..5 simple, 6 takes x,y, 7 simple)
    full_reset: bool

Nothing here interprets the content of a frame.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

GRID = 64
COLORS = 16
RESET, A1, A2, A3, A4, A5, CLICK, A7 = 0, 1, 2, 3, 4, 5, 6, 7


@dataclass(frozen=True)
class Obs:
    frames: tuple[np.ndarray, ...]      # each (64, 64) int8
    state: str
    levels_completed: int
    win_levels: int
    available_actions: tuple[int, ...]
    full_reset: bool

    @property
    def frame(self) -> np.ndarray:
        """The terminal frame of the observation (the state acted on next)."""
        return self.frames[-1]

    @property
    def animated(self) -> bool:
        return len(self.frames) > 1

    @property
    def terminal(self) -> bool:
        return self.state in ("WIN", "GAME_OVER")

    @property
    def simple_actions(self) -> tuple[int, ...]:
        return tuple(a for a in self.available_actions if a not in (RESET, CLICK))

    @property
    def can_click(self) -> bool:
        return CLICK in self.available_actions


def parse(o: dict[str, Any]) -> Obs:
    frames = tuple(np.ascontiguousarray(np.asarray(f, dtype=np.int8)) for f in o["frames"])
    if not frames:
        frames = (np.zeros((GRID, GRID), dtype=np.int8),)
    return Obs(
        frames=frames,
        state=str(o["state"]),
        levels_completed=int(o["levels_completed"]),
        win_levels=int(o.get("win_levels", 0)),
        available_actions=tuple(int(a) for a in o.get("available_actions", [])),
        full_reset=bool(o.get("full_reset", False)),
    )


def one_hot(frame: np.ndarray) -> np.ndarray:
    """Frame as a measure on cells x colors: (16, 64, 64) indicator tensor."""
    f = np.asarray(frame).astype(np.int64)
    out = np.zeros((COLORS, GRID, GRID), dtype=np.int8)
    np.put_along_axis(out, f[None, :, :], 1, axis=0)
    return out


def frame_key(frame: np.ndarray) -> bytes:
    """Exact identity of a frame (for exact verification, never for novelty heuristics)."""
    return np.ascontiguousarray(frame, dtype=np.int8).tobytes()
