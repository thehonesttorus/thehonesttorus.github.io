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
import traceback
from typing import Optional

import numpy as np

from . import algebra, kato, maternity, plan, synth
from .goals import GoalSet
from .split import SplitRules
from .library import Library
from .store import Store
from .words import ClauseSet, Hypothesis, clause_set_log_weight, predicted, predicted_frame, word_context


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
        self.split = SplitRules()
        self.H: dict[int, list[Hypothesis]] = {}
        self.partials: dict[int, list[ClauseSet]] = {}
        self.prev: Optional[algebra.BaseAlgebra] = None
        self.level = 0
        self.actions_since_level = 0
        self.total_actions = 0
        self.last_reason = ""
        self.visited_level: set[bytes] = set()
        self.plan_path: list = []
        self.memory = np.full((64, 64), -1, dtype=np.int16)   # remembered canvas: last colour seen at a cell not covered by a mover type
        self.mover_types: set[tuple[int, bytes]] = set()      # types seen to move, vanish or appear
        self.clock: frozenset = frozenset()                   # the level's autonomous colours (Store.clock_colours)
        self._pcache: dict = {}
        self.diag = dict(synth_calls=0, deaths=0, kato_steps=0, plan_steps=0, fallthrough=0, try_once=0,
                         step_time_max=0.0, step_time_sum=0.0, unexplained=0, goal_probes=0,
                         synth_reshaped=0, synth_fill_inconsistent=0, synth_creation_failed=0, synth_combo_fail=0,
                         synth_no_classes_fail=0, synth_mismatch_hist={}, promotions=0, level_gated=0, death_gated=0,
                         death_unexplained=0)

    # ------------------------------------------------------------------ protocol
    def act(self, obs: dict) -> tuple[int, Optional[dict]]:
        t0 = time.time()
        try:
            if self.prev is None:
                A = self.perceive(obs["frames"][-1])
                self._seed_memory(A)
            else:
                A = self.prev
            self.prev = A
            avail = [a for a in obs.get("available_actions", []) if a != 0] or [1]
            choice = self._decide(A, avail)
        except Exception as e:  # floor policy: never leave a game idle
            print(f"[nce] act error on {self.game_id}: {e!r}\n{traceback.format_exc()}", file=sys.stderr)
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
            print(f"[nce] observe error on {self.game_id}: {e!r}\n{traceback.format_exc()}", file=sys.stderr)
            self.resync(obs)

    def resync(self, obs: dict) -> None:
        """Adopt the observed state without storing a demonstration (a lost response, or recovery after an error).
        Never raises."""
        try:
            frames = obs.get("frames") or []
            if not frames:
                return
            A = self.perceive(frames[-1])
            lc = int(obs.get("levels_completed", self.level))
            if lc > self.level:
                self._on_level_end(A, lc)
            elif lc < self.level:
                self._level_restart(A, lc)
            self.prev = A
        except Exception as e:  # pragma: no cover
            print(f"[nce] resync error on {self.game_id}: {e!r}", file=sys.stderr)
            self.prev = None

    def _level_restart(self, A: algebra.BaseAlgebra, level: int) -> None:
        """The game went back to an earlier level (a full reset): the level's per-visit state restarts; demonstrations
        keep their level labels, so earlier words stay available behind LevelGate."""
        self.level = level
        self.actions_since_level = 0
        self.visited_level = set()
        self.memory[:] = -1
        self._seed_memory(A)
        self.goals.new_level()

    def on_game_over(self, obs: dict) -> None:
        self.diag["deaths"] += 1

    def want_stop(self, obs: dict) -> bool:
        return False

    def perceive(self, frame: np.ndarray) -> algebra.BaseAlgebra:
        """Base algebra through the learned split rules (Section 2.5), memoised per frame and rule set."""
        key = (np.ascontiguousarray(frame, dtype=np.int8).tobytes(), self.split.learned)
        A = self._pcache.get(key)
        if A is None:
            A = self.split.apply(algebra.build_cached(frame))
            if len(self._pcache) >= 1024:
                self._pcache.clear()
            self._pcache[key] = A
        return A

    def _reperceive_store(self) -> None:
        """Split rules changed: rebuild every stored demonstration's algebras and decompositions; words must re-verify."""
        from .match import decompose
        for recs in self.store.records.values():
            for r in recs:
                r.A0 = self.perceive(r.A0.frame)
                r.A1 = self.perceive(r.A1.frame)
                r.tr = decompose(r.A0, r.A1)
                for k in ("_gates", "_tr_short"):            # caches keyed by component ids are stale
                    r.__dict__.pop(k, None)
        self.store._clock_cache = None
        self.H = {}
        self.partials = {}
        if self.prev is not None:
            self.prev = self.perceive(self.prev.frame)

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
            synth_mismatch_hist=self.diag["synth_mismatch_hist"], splits=self.diag.get("splits", 0), split_rules=sum(len(v) for v in self.split.remainders.values()),
            promotions=self.diag["promotions"], library_letters=len(self.lib.dims), beta_c=round(self.lib.beta_c, 4),
            clock_colours=len(self.store.clock_colours()), memorised_words=self.diag.get("memorised_words", 0),
            level_gated=self.diag["level_gated"], death_gated=self.diag["death_gated"], death_unexplained=self.diag["death_unexplained"],
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
        self.goals.resolve_pending(self.H, lambda cs, A, d: predicted_frame(cs, A, d, self.memory))
        Hc = plan.consensus(self.H)
        restrict = None
        if Hc and self.goals.G:
            reach = plan.plan_reachable(A, avail, Hc, self.PLAN_NODES, self.PLAN_TIME_S, perceive=self.perceive, memory=self.memory)
            path = self.goals.split_target(reach, self.store)
            if path:
                first = path[0]
                restrict = {(first[0], None if first[1] is None else (first[1]["x"], first[1]["y"]))}
                choice = kato.kato_select(A, avail, self.H, self.partials, self.store, self.actions_since_level, restrict=restrict, memory=self.memory)
                self.diag["plan_steps"] += 1
                self.last_reason = f"plan a{first[0]} len{len(path)}" + (" +kato" if choice else "")
                return choice if choice else first
        # 3. Kato over all candidates
        choice = kato.kato_select(A, avail, self.H, self.partials, self.store, self.actions_since_level, memory=self.memory)
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
        order = [(a, d) for (a, d) in order if self.store.seen_post(a, A, d) not in (A.key(), b"GAME_OVER")] or order
        for (a, data) in order:
            hyps = self.H.get(a) or []
            if not hyps:
                continue
            fr, flag = predicted(hyps[0].cs, A, data, self.memory)
            if fr is None or flag == "GAME_OVER":
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
        frames = obs.get("frames") or []
        if not frames:
            return
        A1 = self.perceive(frames[-1])
        flag1 = str(obs["state"])
        self.total_actions += 1
        self.actions_since_level += 1
        A0 = self.prev
        if int(obs["levels_completed"]) < self.level:
            self._level_restart(A1, int(obs["levels_completed"]))
            self.prev = A1
            return
        if obs["levels_completed"] > self.level:
            if A0 is not None:
                self.goals.positive(A0, aid, data)      # the completing action: a demonstration of the win correspondence
            self._on_level_end(A1, obs["levels_completed"])
            self.prev = A1
            return
        if A0 is None or obs.get("full_reset"):
            self.prev = A1
            self.memory[:] = -1
            self._seed_memory(A1)
            self.actions_since_level = 0
            return
        if aid == 0:
            self.actions_since_level = 0                 # back at the level's initial state: replay length restarts
        if flag1 == "GAME_OVER":
            self.store.append(aid, A0, data, A1, "GAME_OVER", self.level)
            hyps = [h for h in self.H.get(aid, []) if synth.exact(h.cs, A0, data, A1, "GAME_OVER")]
            for h in hyps:
                h.contexts.add(word_context(h.cs, A0, data))
            if not hyps and self.H.get(aid):
                self._synthesize(aid)             # the word must now carry a death gate
            else:
                self.H[aid] = hyps
            self.prev = A1
            return
        if self.split.learn(A0, A1):                      # Section 2.5: a mover type earned its own stratum
            self.diag["splits"] = self.diag.get("splits", 0) + 1
            A1 = self.perceive(A1.frame)
            self._reperceive_store()
            A0 = self.prev
        rec = self.store.append(aid, A0, data, A1, flag1, self.level, memory=self.memory)
        self._update_memory(rec)
        self.visited_level.add(A1.key())
        self.goals.negative(A1)                          # a non-winning post-state is a certified negative
        if len(frames) > 1:                              # animation sub-steps as serial compositions (stored, v1)
            prev = A0.frame
            for f in frames:
                self.store.append_null(aid, prev, f, self.level)
                prev = f
        clock = self.store.clock_colours()
        if clock != self.clock:
            self.clock = clock
            self.H = {}                            # the flow changed: every word is re-derived in the new quotient
            self.partials = {}
        hyps = self.H.get(aid, [])
        survivors = [h for h in hyps if synth.exact(h.cs, rec.A0, rec.data, rec.A1, rec.flag1, rec.memory, clock)]
        for h in survivors:
            h.n_verified += 1
            h.contexts.add(word_context(h.cs, rec.A0, rec.data))
        if not survivors:
            survivors = self._synthesize(aid)
        self.H[aid] = survivors
        self.prev = A1

    def _synthesize(self, aid: int) -> list[Hypothesis]:
        self.diag["synth_calls"] += 1
        r = synth.synthesize(aid, self.store, self.lib, budget_s=self.SYNTH_BUDGET_S)
        self.partials[aid] = r.partials
        self.diag["unexplained"] += r.unexplained_classes
        self.diag["synth_reshaped"] += r.reshaped_records > 0
        self.diag["synth_fill_inconsistent"] += r.fill_inconsistent
        self.diag["synth_creation_failed"] += r.creation_failed
        self.diag["level_gated"] += r.level_gated
        self.diag["death_gated"] += r.death_gated
        self.diag["death_unexplained"] += r.death_unexplained
        self.diag["memorised_words"] = self.diag.get("memorised_words", 0) + (r.memorised > 0)
        if not r.hypotheses and r.combos_tried:
            self.diag["synth_combo_fail"] += 1
            b = "0" if r.min_mismatch_cells == 0 else ("1-4" if r.min_mismatch_cells <= 4 else ("5-32" if r.min_mismatch_cells <= 32 else ("33-256" if r.min_mismatch_cells <= 256 else ">256")))
            self.diag["synth_mismatch_hist"][b] = self.diag["synth_mismatch_hist"].get(b, 0) + 1
        self.H[aid] = r.hypotheses
        return r.hypotheses

    def _update_memory(self, rec) -> None:
        """Canvas memory: cells not covered by a mover-type component take the observed colour (strata persist
        under occlusion). When a new mover type is identified the whole history is re-read (Store.recompute_memories)."""
        tr = rec.tr
        before = len(self.mover_types)
        for i, j, g in tr.moved:
            self.mover_types.add(rec.A0.comps[i].type); self.mover_types.add(rec.A1.comps[j].type)
        for i in tr.annihilated:
            self.mover_types.add(rec.A0.comps[i].type)
        for j in tr.created:
            self.mover_types.add(rec.A1.comps[j].type)
        if len(self.mover_types) != before:
            self.memory = self.store.recompute_memories(self.mover_types)
            self.H = {}                            # fills were fitted against the old memories
            return
        covered = np.zeros(64 * 64, dtype=bool)
        for c in rec.A1.comps:
            if c.type in self.mover_types:
                covered[c.cells] = True
        mem = self.memory.ravel()
        frame = rec.A1.frame.ravel()
        mem[~covered] = frame[~covered]

    def _seed_memory(self, A: algebra.BaseAlgebra) -> None:
        covered = np.zeros(64 * 64, dtype=bool)
        for c in A.comps:
            if c.type in self.mover_types:
                covered[c.cells] = True
        self.memory.ravel()[~covered] = A.frame.ravel()[~covered]

    def _on_level_end(self, A_new: algebra.BaseAlgebra, new_level: int) -> None:
        self.level = new_level
        self.actions_since_level = 0
        self.visited_level = set()
        self.memory[:] = -1
        self._seed_memory(A_new)
        self.goals.new_level()
        self.lib.snapshot_alphabet(self._n_types())     # T11: refresh the alphabet snapshot
        # the maternity ward (T5, T10): exact promotion against the corpus of verified words
        lib2, names = maternity.promote(self.lib, maternity.corpus_of(self.H))
        if names:
            self.lib = lib2
            self.goals.lib = lib2
            self.diag["promotions"] += len(names)
            for hs in self.H.values():
                for h in hs:
                    h.logw = clause_set_log_weight(h.cs, lib2)

    def _n_types(self) -> int:
        keys = set()
        for recs in self.store.records.values():
            for r in recs:
                keys.update(c.key for c in r.A0.frame_vectors())
        return len(keys)
