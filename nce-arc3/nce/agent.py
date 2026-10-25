"""NCEAgent: the in-episode correspondence solver (Section 10, harness protocol).

    agent = NCEAgent(game_id, win_levels, available_actions)
    action_id, data = agent.act(obs)
    agent.observe(obs_after, action_id, data)
    agent.on_game_over(obs)

obs = dict(frames=[np.int8 (64,64), ...], state, levels_completed, win_levels, available_actions, full_reset)

Decision order each turn (Sections 5-6, C5): try every available action once; else, if the
goal coherent state has a reachable splitting target under the consensus model, take the
plan's first action (the step is itself an experiment where hypotheses disagree: T7); else
the Kato rule; else the goal-observable newborn fallthrough.
"""
from __future__ import annotations

import sys
import time
from typing import Optional

import numpy as np

from . import algebra, kato, plan, synth
from .goals import GoalSet
from .library import Library
from .store import Store
from .words import ClauseSet, Hypothesis, predicted_frame


class NCEAgent:
    SYNTH_BUDGET_S = 0.25
    PLAN_NODES = 60
    PLAN_TIME_S = 0.08

    def __init__(self, game_id: str, win_levels: int, available_actions: list[int]):
        self.game_id = game_id
        self.win_levels = int(win_levels)
        self.lib = Library(win_levels=max(2, self.win_levels))
        self.store = Store()
        self.goals = GoalSet(self.lib)
        self.H: dict[int, list[Hypothesis]] = {}
        self.partials: dict[int, list[ClauseSet]] = {}
        self.prev: Optional[algebra.BaseAlgebra] = None
        self.level = 0
        self.actions_since_level = 0
        self.total_actions = 0
        self.last_reason = ""
        self.visited_level: set[bytes] = set()
        self.plan_path: list = []
        self.diag = dict(synth_calls=0, deaths=0, kato_steps=0, plan_steps=0, fallthrough=0, try_once=0,
                         step_time_max=0.0, step_time_sum=0.0, unexplained=0, goal_probes=0,
                         synth_reshaped=0, synth_fill_inconsistent=0, synth_creation_failed=0, synth_combo_fail=0,
                         synth_no_classes_fail=0, synth_mismatch_hist={})

    # ------------------------------------------------------------------ protocol
    def act(self, obs: dict) -> tuple[int, Optional[dict]]:
        t0 = time.time()
        try:
            A = self.prev if self.prev is not None else algebra.build_cached(obs["frames"][-1])
            self.prev = A
            avail = [a for a in obs.get("available_actions", []) if a != 0] or [1]
            choice = self._decide(A, avail)
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
            self.prev = algebra.build_cached(obs["frames"][-1])

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
            plan=self.diag["plan_steps"], fallthrough=self.diag["fallthrough"], deaths=self.diag["deaths"],
            goals_alive=len(self.goals.G), goals_killed=self.goals.killed, goals_confirmed=self.goals.confirmed,
            goals_pending=len(self.goals.pending), goals_regen=self.goals.regenerated,
            step_ms_avg=round(1000 * self.diag["step_time_sum"] / n, 1), step_ms_max=round(1000 * self.diag["step_time_max"], 1),
            unexplained_classes=self.diag["unexplained"], types=self._n_types(),
            synth_reshaped=self.diag["synth_reshaped"], synth_fill_inconsistent=self.diag["synth_fill_inconsistent"],
            synth_creation_failed=self.diag["synth_creation_failed"], synth_combo_fail=self.diag["synth_combo_fail"],
            synth_mismatch_hist=self.diag["synth_mismatch_hist"],
        )

    # ------------------------------------------------------------------ decision
    def _decide(self, A: algebra.BaseAlgebra, avail: list[int]) -> tuple[int, Optional[dict]]:
        # 1. try once
        for a in avail:
            if not self.store.tried(a):
                self.diag["try_once"] += 1
                self.last_reason = f"try-once a{a}"
                return (a, kato.representative_data(a, A))
        # 2. goal-splitting over the consensus model (Section 6.4), with decision-relevant interleaving (T7)
        self.goals.ensure(A)
        self.goals.resolve_pending(self.H, predicted_frame)
        Hc = plan.consensus(self.H)
        restrict = None
        if Hc and self.goals.G:
            reach = plan.plan_reachable(A, avail, Hc, self.PLAN_NODES, self.PLAN_TIME_S)
            path = self.goals.split_target(reach, self.store)
            if path:
                first = path[0]
                restrict = {(first[0], None if first[1] is None else (first[1]["x"], first[1]["y"]))}
                choice = kato.kato_select(A, avail, self.H, self.partials, self.store, self.actions_since_level, restrict=restrict)
                self.diag["plan_steps"] += 1
                self.last_reason = f"plan a{first[0]} len{len(path)}" + (" +kato" if choice else "")
                return choice if choice else first
        # 3. Kato over all candidates
        choice = kato.kato_select(A, avail, self.H, self.partials, self.store, self.actions_since_level)
        if choice is not None:
            self.diag["kato_steps"] += 1
            self.last_reason = f"kato a{choice[0]}"
            return choice
        # 4. fallthrough (C5 third clause, v1: a predicted state not yet visited this level)
        self.diag["fallthrough"] += 1
        return self._fallthrough(A, avail)

    def _fallthrough(self, A: algebra.BaseAlgebra, avail: list[int]) -> tuple[int, Optional[dict]]:
        cands = kato.candidates(A, avail)
        if not cands:
            return (avail[0], None)
        start = self.total_actions % len(cands)
        order = cands[start:] + cands[:start]
        fallback = None
        for (a, data) in order:
            hyps = self.H.get(a) or []
            if not hyps:
                continue
            fr = predicted_frame(hyps[0].cs, A, data)
            if fr is None:
                continue
            key = np.ascontiguousarray(fr, dtype=np.int8).tobytes()
            if key == A.key():
                continue  # predicted no-op (C8)
            if key not in self.visited_level:
                self.last_reason = f"unvisited a{a}"
                return (a, data)
            if fallback is None:
                fallback = (a, data)
        if fallback is None:
            fallback = order[0]
        self.last_reason = f"fallthrough a{fallback[0]}"
        return fallback

    # ------------------------------------------------------------------ learning
    def _observe(self, obs: dict, aid: int, data: Optional[dict]) -> None:
        frames = obs["frames"]
        A1 = algebra.build_cached(frames[-1])
        flag1 = str(obs["state"])
        self.total_actions += 1
        self.actions_since_level += 1
        A0 = self.prev
        if obs["levels_completed"] > self.level:
            if A0 is not None:
                self.goals.positive(A0, aid, data)      # the completing action: a demonstration of the win correspondence
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
        self.goals.negative(A1)                          # a non-winning post-state is a certified negative
        if len(frames) > 1:                              # animation sub-steps as serial compositions (stored, v1)
            prev = A0
            for f in frames:
                Af = algebra.build_cached(f)
                self.store.append_null(aid, prev, Af, self.level)
                prev = Af
        hyps = self.H.get(aid, [])
        survivors = [h for h in hyps if synth.exact(h.cs, rec.A0, rec.data, rec.A1, rec.flag1)]
        for h in survivors:
            h.n_verified += 1
        if not survivors:
            self.diag["synth_calls"] += 1
            r = synth.synthesize(aid, self.store, self.lib, budget_s=self.SYNTH_BUDGET_S)
            survivors = r.hypotheses
            self.partials[aid] = r.partials
            self.diag["unexplained"] += r.unexplained_classes
            self.diag["synth_reshaped"] += r.reshaped_records > 0
            self.diag["synth_fill_inconsistent"] += r.fill_inconsistent
            self.diag["synth_creation_failed"] += r.creation_failed
            if not survivors and r.combos_tried:
                self.diag["synth_combo_fail"] += 1
                b = "0" if r.min_mismatch_cells == 0 else ("1-4" if r.min_mismatch_cells <= 4 else ("5-32" if r.min_mismatch_cells <= 32 else ("33-256" if r.min_mismatch_cells <= 256 else ">256")))
                self.diag["synth_mismatch_hist"][b] = self.diag["synth_mismatch_hist"].get(b, 0) + 1
        self.H[aid] = survivors
        self.prev = A1

    def _on_level_end(self, A_new: algebra.BaseAlgebra, new_level: int) -> None:
        self.level = new_level
        self.actions_since_level = 0
        self.visited_level = set()
        self.goals.new_level()
        self.lib.snapshot_alphabet(self._n_types())     # T11: refresh the alphabet snapshot

    def _n_types(self) -> int:
        keys = set()
        for recs in self.store.records.values():
            for r in recs:
                keys.update(c.key for c in r.A0.frame_vectors())
        return len(keys)
