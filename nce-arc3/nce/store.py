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


@dataclass
class Record:
    A0: BaseAlgebra
    data: Optional[dict]
    A1: BaseAlgebra
    flag1: str
    level: int
    tr: Transition

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

    def append(self, a: int, A0: BaseAlgebra, data: Optional[dict], A1: BaseAlgebra, flag1: str, level: int) -> Record:
        rec = Record(A0, data, A1, flag1, level, decompose(A0, A1))
        self.records.setdefault(a, []).append(rec)
        table = self.ignorance.setdefault(a, {})
        k = rec.key
        post = A1.key() if flag1 != "GAME_OVER" else b"GAME_OVER"
        if k in table and table[k] != post:
            self.nondeterministic.add(a)
        table[k] = post
        return rec

    def append_null(self, a: int, A0: BaseAlgebra, A1: BaseAlgebra, level: int) -> None:
        self.null_records.setdefault(a, []).append(Record(A0, None, A1, "NOT_FINISHED", level, decompose(A0, A1)))

    def transitions(self, a: int) -> list[Record]:
        return self.records.get(a, [])

    def tried(self, a: int) -> bool:
        return bool(self.records.get(a))

    def count(self) -> int:
        return sum(len(v) for v in self.records.values())

    def seen_post(self, a: int, A0: BaseAlgebra, data: Optional[dict]) -> Optional[bytes]:
        d = (int(data["x"]), int(data["y"])) if data else None
        return self.ignorance.get(a, {}).get((A0.key(), d))

    def retire_level(self, a: int, level: int) -> None:
        self.records[a] = [r for r in self.records.get(a, []) if r.level >= level]
