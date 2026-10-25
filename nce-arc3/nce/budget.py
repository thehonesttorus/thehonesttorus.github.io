"""Action accounting under RHAE.

Scoring facts (technical report, Section 4, and arc_agi/scorecard.py):
  S_l = min(1.15, (h_l / a_l)^2) per completed level, a_l = actions since the previous
  level completion (every submitted action counts, RESETs included);
  E = sum_l l*S_l / sum_l l, capped at (sum over completed l) / sum_l l;
  the official evaluation stops a level after 5 h_l actions.

Human baselines are unknown for hidden games, so the engine's budget policy is
baseline-free: it tracks actions per level and wall clock, and exposes the RHAE
estimate only when baselines are supplied (local scoreboard).
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional, Sequence


@dataclass
class LevelLedger:
    actions_per_completed_level: list[int] = field(default_factory=list)
    since_level: int = 0
    total: int = 0
    resets: int = 0
    game_overs: int = 0
    started: float = field(default_factory=time.time)

    def record(self, action_id: int, levels_completed_after: int, state_after: str) -> None:
        self.total += 1
        self.since_level += 1
        if action_id == 0:
            self.resets += 1
        if state_after == "GAME_OVER":
            self.game_overs += 1
        completed = len(self.actions_per_completed_level)
        while completed < levels_completed_after:
            self.actions_per_completed_level.append(self.since_level)
            self.since_level = 0
            completed += 1

    @property
    def elapsed(self) -> float:
        return time.time() - self.started

    def rhae(self, baselines: Optional[Sequence[int]], win_levels: int) -> Optional[float]:
        if not baselines:
            return None
        n = win_levels or len(baselines)
        w = list(range(1, n + 1))
        tot = sum(w)
        s = 0.0
        for l, a in enumerate(self.actions_per_completed_level, start=1):
            if l - 1 < len(baselines) and a > 0:
                s += w[l - 1] * min(1.15, (baselines[l - 1] / a) ** 2)
        cap = sum(w[: len(self.actions_per_completed_level)]) / tot
        return min(s / tot, cap)
