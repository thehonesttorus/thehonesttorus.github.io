"""NCEAgent: the in-episode correspondence solver (Section 10, harness protocol).

    agent = NCEAgent(game_id, win_levels, available_actions)
    action_id, data = agent.act(obs)
    agent.observe(obs_after, action_id, data)
    agent.on_game_over(obs)

obs = dict(frames=[np.int8 (64,64), ...], state, levels_completed, win_levels, available_actions, full_reset)
"""
from __future__ import annotations

import sys
import time
from typing import Optional

import numpy as np

from . import algebra, kato, synth
from .library import Library
from .store import Store
from .words import ClauseSet, Hypothesis, predicted_frame


class NCEAgent:
    SOFT_STEP_S = 0.2
    SYNTH_BUDGET_S = 0.25

    def __init__(self, game_id: str, win_levels: int, available_actions: list[int]):
        self.game_id = game_id
        self.win_levels = int(win_levels)
        self.lib = Library(win_levels=max(2, self.win_levels))
        self.store = Store()
        self.H: dict[int, list[Hypothesis]] = {}
        self.partials: dict[int, list[ClauseSet]] = {}
        self.prev: Optional[algebra.BaseAlgebra] = None
        self.level = 0
        self.actions_since_level = 0
        self.total_actions = 0
        self.last_reason = ""
        self.visited_level: set[bytes] = set()
        self.diag = dict(synth_calls=0, deaths=0, agreements=0, kato_steps=0, fallthrough=0, try_once=0, step_time_max=0.0, step_time_sum=0.0)

    # ------------------------------------------------------------------ protocol
    def act(self, obs: dict) -> tuple[int, Optional[dict]]:
        t0 = time.time()
        try:
            A = self.prev if self.prev is not None else algebra.build(obs["frames"][-1])
            self.prev = A
            avail = [a for a in obs.get("available_actions", []) if a != 0] or [1]
            choice = kato.kato_select(A, avail, self.H, self.partials, self.store, self.actions_since_level)
            if choice is not None:
                if not self.store.tried(choice[0]):
                    self.diag["try_once"] += 1
                    self.last_reason = f"try-once a{choice[0]}"
                else:
                    self.diag["kato_steps"] += 1
                    self.last_reason = f"kato a{choice[0]}"
            else:
                choice = self._fallthrough(A, avail)
                self.diag["fallthrough"] += 1
        except Exception as e:  # floor policy: never leave a game idle
            print(f"[nce] act error on {self.game_id}: {e!r}", file=sys.stderr)
            avail = [a for a in obs.get("available_actions", []) if a != 0] or [1]
            a = avail[self.total_actions % len(avail)]
            choice = (a, {"x": 0, "y": 0} if a == 6 else None)
            self.last_reason = "floor"
        dt = time.time() - t0
        self.diag["step_time_max"] = max(self.diag["step_time_max"], dt)
        self.diag["step_time_sum"] += dt
        return int(choice[0]), choice[1]

    def observe(self, obs: dict, aid: int, data: Optional[dict]) -> None:
        try:
            self._observe(obs, int(aid), data)
        except Exception as e:
            print(f"[nce] observe error on {self.game_id}: {e!r}", file=sys.stderr)
            self.prev = algebra.build(obs["frames"][-1])

    def on_game_over(self, obs: dict) -> None:
        self.diag["deaths"] += 1

    def want_stop(self, obs: dict) -> bool:
        return False

    def diagnostics(self) -> dict:
        """Agent-internal counts only (blind-scoreboard discipline): never words, never frames."""
        n = max(self.total_actions, 1)
        return dict(
            actions=self.total_actions, level=self.level,
            tried={a: len(v) for a, v in sorted(self.store.records.items())},
            H_sizes={a: len(v) for a, v in sorted(self.H.items())},
            partials={a: len(v) for a, v in sorted(self.partials.items()) if v},
            nondeterministic=sorted(self.store.nondeterministic),
            synth_calls=self.diag["synth_calls"], try_once=self.diag["try_once"], kato=self.diag["kato_steps"],
            fallthrough=self.diag["fallthrough"], deaths=self.diag["deaths"],
            step_ms_avg=round(1000 * self.diag["step_time_sum"] / n, 1), step_ms_max=round(1000 * self.diag["step_time_max"], 1),
            unexplained_classes=self.diag.get("unexplained", 0), types=len({c.key for recs in self.store.records.values() for r in recs for c in r.A0.comps}),
        )

    # ------------------------------------------------------------------ learning
    def _observe(self, obs: dict, aid: int, data: Optional[dict]) -> None:
        frames = obs["frames"]
        A1 = algebra.build(frames[-1])
        flag1 = str(obs["state"])
        self.total_actions += 1
        self.actions_since_level += 1
        A0 = self.prev
        if obs["levels_completed"] > self.level:
            self._on_level_end(A1, obs["levels_completed"])
            self.prev = A1
            return
        if A0 is None or obs.get("full_reset"):
            self.prev = A1
            return
        if flag1 == "GAME_OVER":
            self.store.append(aid, A0, data, A1, "GAME_OVER", self.level)
            self.prev = A1
            return
        rec = self.store.append(aid, A0, data, A1, flag1, self.level)
        self.visited_level.add(A1.key())
        # animation sub-steps as serial compositions under the null action (stored only, v1)
        if len(frames) > 1:
            prev = A0
            for f in frames:
                Af = algebra.build(f)
                self.store.append_null(aid, prev, Af, self.level)
                prev = Af
        # residual test: hypotheses that fail the new demonstration die (T2, T6)
        hyps = self.H.get(aid, [])
        survivors = [h for h in hyps if synth.exact(h.cs, rec.A0, rec.data, rec.A1, rec.flag1)]
        for h in survivors:
            h.n_verified += 1
        if hyps and not survivors:
            self.diag["deaths"] += 0  # placeholder for a finer diagnostic
        if not survivors:
            self.diag["synth_calls"] += 1
            r = synth.synthesize(aid, self.store, self.lib, budget_s=self.SYNTH_BUDGET_S)
            survivors = r.hypotheses
            self.partials[aid] = r.partials
            self.diag["unexplained"] = self.diag.get("unexplained", 0) + r.unexplained_classes
            if survivors:
                self.last_reason += f" | a{aid}: {survivors[0].cs.describe()}"
        self.H[aid] = survivors
        self.prev = A1

    def _on_level_end(self, A_new: algebra.BaseAlgebra, new_level: int) -> None:
        self.level = new_level
        self.actions_since_level = 0
        self.visited_level = set()
        # carry H and the store across levels (Section 8.5); refresh the alphabet snapshot (T11)
        n_types = len({c.key for recs in self.store.records.values() for r in recs for c in r.A0.comps})
        self.lib.snapshot_alphabet(n_types)

    # ------------------------------------------------------------------ fallthrough (C5 third clause, v1)
    def _fallthrough(self, A: algebra.BaseAlgebra, avail: list[int]) -> tuple[int, Optional[dict]]:
        """All hypotheses agree and nothing is newborn: move to a predicted state not yet visited this level."""
        best = None
        for (a, data) in kato.candidates(A, avail):
            hyps = self.H.get(a) or []
            if not hyps:
                continue
            fr = predicted_frame(hyps[0].cs, A, data)
            if fr is None:
                continue
            key = np.ascontiguousarray(fr, dtype=np.int8).tobytes()
            if key == A.key():
                continue  # predicted no-op: wasted action (C8)
            if key not in self.visited_level:
                self.last_reason = f"unvisited-state a{a}"
                return (a, data)
            if best is None:
                best = (a, data)
        if best is None:
            a = avail[self.total_actions % len(avail)]
            best = (a, kato.representative_data(a, A))
        self.last_reason = f"fallthrough a{best[0]}"
        return best
