"""Stored transitions (demonstrations, T1) and the ignorance table (the compact part, T12).

Spec Sections 4.1 and 4.7. Nothing unexplained is ever dropped. The ignorance table never
selects an action; it is the non-determinism detector: two transitions of the same action
with identical (pre-frame, datum) but different post-frames admit no word with constant
parameters, so the action is marked non-deterministic on observables.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .algebra import BaseAlgebra
from .match import Transition, decompose


def _mask_bytes(rec: "Record", c: int) -> tuple[bytes, bytes]:
    """Packed masks of colour c before and after the transition (cached on the record)."""
    cache = rec.__dict__.setdefault("_masks", {})
    if c not in cache:
        cache[c] = (np.packbits(rec.A0.frame == c).tobytes(), np.packbits(rec.A1.frame == c).tobytes())
    return cache[c]


def _pointed_key(rec: "Record"):
    if not rec.data:
        return None
    cid = int(rec.A0.cell2comp[int(rec.data["y"]), int(rec.data["x"])])
    return rec.A0.comps[cid].key if cid >= 0 else None


def _own_frame_change(rec: "Record", c: int) -> list[tuple]:
    """The change of colour c in its own frame: lost and gained cells as offsets from the pre-mask's bbox corner,
    read from the min corner and from the max corner (two descriptors)."""
    m0 = rec.A0.frame == c
    m1 = rec.A1.frame == c
    ys, xs = np.nonzero(m0)
    if ys.size == 0:
        return []
    lost = np.argwhere(m0 & ~m1)
    gained = np.argwhere(m1 & ~m0)
    lo = np.array([ys.min(), xs.min()])
    hi = np.array([ys.max(), xs.max()])
    def rel(pts, o):
        return tuple(sorted(map(tuple, (pts - o).tolist())))
    return [("min", rel(lost, lo), rel(gained, lo)), ("max", rel(lost, hi), rel(gained, hi))]


def _change_descriptor(rec: "Record", c: int) -> tuple:
    """What happened to colour c, up to position: cells lost and gained, and the colours they went to / came from."""
    lost = (rec.A0.frame == c) & (rec.A1.frame != c)
    gained = (rec.A0.frame != c) & (rec.A1.frame == c)
    return (int(lost.sum()), int(gained.sum()),
            tuple(sorted(int(v) for v in np.unique(rec.A1.frame[lost]))),
            tuple(sorted(int(v) for v in np.unique(rec.A0.frame[gained]))))


@dataclass
class Record:
    A0: BaseAlgebra
    data: Optional[dict]
    A1: BaseAlgebra
    flag1: str
    level: int
    tr: Transition
    memory: Optional[np.ndarray] = None     # remembered canvas (64x64, -1 unknown) before the transition

    @property
    def key(self) -> tuple[bytes, Optional[tuple[int, int]]]:
        d = (int(self.data["x"]), int(self.data["y"])) if self.data else None
        return (self.A0.key(), d)


@dataclass
class Store:
    records: dict[int, list[Record]] = field(default_factory=dict)        # action id -> transitions
    null_records: dict[int, list[Record]] = field(default_factory=dict)   # animation sub-steps per action
    ignorance: dict[int, dict] = field(default_factory=dict)              # action -> {(pre_key, datum): post_key}
    nondeterministic: set[int] = field(default_factory=set)
    order: list[Record] = field(default_factory=list)                     # all records in time order
    _clock_cache: Optional[tuple[int, frozenset]] = None

    def append(self, a: int, A0: BaseAlgebra, data: Optional[dict], A1: BaseAlgebra, flag1: str, level: int,
               memory: Optional[np.ndarray] = None) -> Record:
        rec = Record(A0, data, A1, flag1, level, decompose(A0, A1), None if memory is None else memory.copy())
        self.records.setdefault(a, []).append(rec)
        self.order.append(rec)
        table = self.ignorance.setdefault(a, {})
        k = rec.key
        post = A1.key() if flag1 != "GAME_OVER" else b"GAME_OVER"
        if k in table and table[k] != post:
            self.nondeterministic.add(a)
        table[k] = post
        return rec

    def append_null(self, a: int, frame0: np.ndarray, frame1: np.ndarray, level: int) -> None:
        """Animation sub-steps are kept as data (serial compositions of the null action, N7): the frame pairs. No rule
        reads their decomposition yet, so it is neither perceived nor computed."""
        self.null_records.setdefault(a, []).append((np.array(frame0, dtype=np.int8), np.array(frame1, dtype=np.int8), level))

    def transitions(self, a: int) -> list[Record]:
        return self.records.get(a, [])

    def tried(self, a: int) -> bool:
        return bool(self.records.get(a))

    def count(self) -> int:
        return sum(len(v) for v in self.records.values())

    def seen_post(self, a: int, A0: BaseAlgebra, data: Optional[dict]) -> Optional[bytes]:
        d = (int(data["x"]), int(data["y"])) if data else None
        return self.ignorance.get(a, {}).get((A0.key(), d))

    def recompute_memories(self, mover_types: set, first_frame: Optional[np.ndarray] = None) -> np.ndarray:
        """Canvas memory per record from the stored history: a cell remembers the colour of the most recent
        frame in which no mover-type component covered it. Returns the live memory after the last record."""
        mem = np.full(64 * 64, -1, dtype=np.int16)

        def write(A: BaseAlgebra) -> None:
            covered = np.zeros(64 * 64, dtype=bool)
            for c in A.comps:
                if c.type in mover_types:
                    covered[c.cells] = True
            mem[~covered] = A.frame.ravel()[~covered]

        prev_level = None
        for rec in self.order:
            if prev_level is None or rec.level != prev_level:
                mem[:] = -1                       # a new level is a new canvas
                write(rec.A0)
            rec.memory = mem.copy().reshape(64, 64)
            if rec.flag1 != "GAME_OVER":
                write(rec.A1)
            prev_level = rec.level
        return mem.reshape(64, 64)

    # ---- the action-independent flow (T2, T12) -------------------------------------------------------------
    def clock_colours(self, min_records: int = 5) -> frozenset[int]:
        """Colours that evolve autonomously of the action (the level's modular flow sigma_t, T2/T12). Over the current
        level's non-terminal records (null steps included), a colour c is a clock colour iff
          (a) its post-state is a function of its own pre-state (no pre-mask of c ever had two post-masks),
          (b) it changes in every such record, at least `min_records` times (a change no action withholds),
          (c) it is not the canvas (the top stratum only receives what others vacate; rendering convention), and
          (d) action-independence is witnessed: the same change in the colour's own frame (lost and gained cells
              relative to its pre-bbox, from either corner) occurs under two different action contexts, a context
              being the action id and the pointed component's shape type.
        Verification is then taken modulo its cells, and it is never a word of any action."""
        if self._clock_cache is not None and self._clock_cache[0] == len(self.order):
            return self._clock_cache[1]
        level = self.order[-1].level if self.order else 0
        items = [(a, r) for a, recs in self.records.items() for r in recs if r.level == level and r.flag1 != "GAME_OVER"]
        out: set[int] = set()
        if len(items) >= min_records and all(r.tr.changed_cells > 0 for _, r in items):
            cands = None
            for _, r in items:
                d = r.A0.frame != r.A1.frame
                ch = set(np.unique(r.A0.frame[d]).tolist()) | set(np.unique(r.A1.frame[d]).tolist())
                cands = ch if cands is None else (cands & ch)
                if not cands:
                    break
            for c in sorted(cands or ()):
                c = int(c)
                if any(len(r.A0.strata) > 1 and c in r.A0.strata[0] for _, r in items):
                    continue
                posts: dict[bytes, bytes] = {}
                contexts: dict[tuple, set] = {}
                ok = True
                for a, r in items:
                    pre, post = _mask_bytes(r, c)
                    if pre == post or posts.setdefault(pre, post) != post:
                        ok = False
                        break
                    ctx = (a, _pointed_key(r))
                    for desc in _own_frame_change(r, c):
                        contexts.setdefault(desc, set()).add(ctx)
                if ok and any(len(v) >= 2 for v in contexts.values()):
                    out.add(c)
        res = frozenset(out)
        self._clock_cache = (len(self.order), res)
        return res

    def retire_level(self, a: int, level: int) -> None:
        self.records[a] = [r for r in self.records.get(a, []) if r.level >= level]
