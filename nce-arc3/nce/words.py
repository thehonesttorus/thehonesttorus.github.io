"""Words, clauses, clause sets: prediction, rendering, fitting and exact verification.

Spec: final-architecture.md Sections 3.3, 3.5, 4.3-4.5.

A hypothesis for an action id is a ClauseSet: a parallel composition of clauses acting on
orthogonal summands of the base algebra (components), plus creation clauses and an optional
Over (GAME_OVER) flag. A clause is

    selector (F->C) . [Bind . Gate]? . body (C->C letters with fitted units) . terminator

with the `Fill[c]` letter attached whenever the clause vacates cells. Components selected by
no clause are acted on by the identity (the empty word, the vacuum).

Letters here are the families of Section 3.3 with their fitted parameters:
  selectors:  ("all",) ("pointed",) ("colour", c) ("rank", r) ("type", key)
  bind/gate:  ("contact", type_or_None) ("coincide", type_or_None) ("adjacent1", type_or_None)
              ("samecolour", None) ("sameshape", None)    -- apply body iff the predicate holds
  body:       ("step", (dy, dx)) ("move", (dy, dx)) ("rot", h) ("recolour", c) ("copy", (dy, dx))
  terminator: ("stay",) ("annihilate",)
  fill:       colour c of vacated cells (None when nothing is vacated)

The exact test (T6) is cell-for-cell equality of render(predict(h, Psi_t)) with Psi_{t+1}.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Optional

import numpy as np

from .algebra import BaseAlgebra, Component, GRID, d4_apply, d4_inverse, transport_cells

Param = Any
Letter = tuple  # (kind, param)

# D4 composition table: compose(a, b)(x) = d4_apply(d4_apply(x, b), a)
_D4_COMPOSE: Optional[np.ndarray] = None


def d4_compose(a: int, b: int) -> int:
    """Element c with d4_apply(x, c) == d4_apply(d4_apply(x, b), a) for every array x."""
    global _D4_COMPOSE
    if _D4_COMPOSE is None:
        probe = np.arange(6).reshape(2, 3)
        table = np.zeros((8, 8), dtype=np.int64)
        for i in range(8):
            for j in range(8):
                target = d4_apply(d4_apply(probe, j), i)
                for c in range(8):
                    if d4_apply(probe, c).shape == target.shape and np.array_equal(d4_apply(probe, c), target):
                        table[i, j] = c
                        break
        _D4_COMPOSE = table
    return int(_D4_COMPOSE[a, b])


def d4_relative(g_before: int, g_after: int) -> int:
    """h with g_after = compose(h, g_before)."""
    for h in range(8):
        if d4_compose(h, g_before) == g_after:
            return h
    raise RuntimeError("D4 relative element not found")  # pragma: no cover


# ------------------------------------------------------------------------------ data model

@dataclass(frozen=True)
class Clause:
    selector: Letter
    bind: Optional[Letter]
    body: tuple[Letter, ...]
    term: Letter
    fill: Optional[int]

    def letters(self) -> list[str]:
        """Letter family names along the Bratteli path (for KMS weights)."""
        names = [{"all": "SelAll", "pointed": "SelPointed", "colour": "SelColour", "rank": "SelRank", "type": "SelType"}[self.selector[0]]]
        if self.bind is not None:
            names += ["Bind", "Gate"]
        for kind, _ in self.body:
            names.append({"step": "Step", "move": "Move", "rot": "Rot", "recolour": "Recolour", "copy": "Copy"}[kind])
        if self.fill is not None:
            names.append("Fill")
        if self.term[0] == "annihilate":
            names.append("Annihilate")
        return names

    def vacates(self) -> bool:
        return self.term[0] == "annihilate" or any(k in ("step", "move", "rot") for k, _ in self.body)

    def describe(self) -> str:
        s = f"{self.selector[0]}" + (f"[{self.selector[1]}]" if len(self.selector) > 1 and self.selector[1] is not None and self.selector[0] != 'type' else "")
        if self.selector[0] == "type":
            s += "[type]"
        if self.bind:
            s += f"·gate({self.bind[0]})"
        for k, p in self.body:
            s += f"·{k}{p if k != 'rot' else '[' + str(p) + ']'}"
        if self.fill is not None:
            s += f"·fill{self.fill}"
        if self.term[0] == "annihilate":
            s += "·annihilate"
        return s


@dataclass(frozen=True)
class Creation:
    key: bytes            # canonical shape key (must be in the alphabet)
    color: int
    g: int
    offset: tuple[int, int]      # absolute bbox corner, or offset relative to the pointed component when relative
    relative: bool = False       # Copy/placement relative to the pointed component (SelPointed . Copy . units)

    def letters(self) -> list[str]:
        return ["SelPointed", "Copy"] if self.relative else ["Create"]

    def cells(self, A: BaseAlgebra, data: Optional[dict]) -> Optional[np.ndarray]:
        from .algebra import canonical_image_from_key, cells_from_canonical
        canon = canonical_image_from_key(self.key)
        if self.relative:
            if not data:
                return None
            cid = int(A.cell2comp[int(data["y"]), int(data["x"])])
            if cid < 0:
                return None
            base = A.comps[cid].offset
            off = (base[0] + self.offset[0], base[1] + self.offset[1])
        else:
            off = self.offset
        return cells_from_canonical(canon, self.g, off)


@dataclass(frozen=True)
class ClauseSet:
    clauses: tuple[Clause, ...] = ()
    creations: tuple[Creation, ...] = ()
    over: bool = False

    def letters(self) -> list[list[str]]:
        out = [c.letters() for c in self.clauses] + [c.letters() for c in self.creations]
        if self.over:
            out.append(["Over"])
        return out

    def describe(self) -> str:
        parts = [c.describe() for c in self.clauses] + [f"create(c{c.color})" for c in self.creations]
        if self.over:
            parts.append("over")
        return " ⊕ ".join(parts) if parts else "identity"


@dataclass
class Hypothesis:
    cs: ClauseSet
    logw: float
    n_verified: int = 0
    beyond_ceiling: bool = False

    def describe(self) -> str:
        return f"{self.cs.describe()} (logw {self.logw:.2f}, verified {self.n_verified})"


# ------------------------------------------------------------------------------ selection

def select(selector: Letter, A: BaseAlgebra, data: Optional[dict]) -> list[Component]:
    kind = selector[0]
    if kind == "all":
        return list(A.comps)
    if kind == "pointed":
        if not data:
            return []
        x, y = int(data["x"]), int(data["y"])
        cid = int(A.cell2comp[y, x])
        return [A.comps[cid]] if cid >= 0 else []
    if kind == "colour":
        return [c for c in A.comps if c.color == selector[1]]
    if kind == "rank":
        return [c for c in A.comps if c.stratum == selector[1]]
    if kind == "type":
        return [c for c in A.comps if c.key == selector[1]]
    raise ValueError(kind)


def _pair_predicate(bind: Letter, comp: Component, A: BaseAlgebra) -> bool:
    kind, other_key = bind
    others = [o for o in A.comps if o.id != comp.id and (other_key is None or o.key == other_key)]
    if kind == "contact":
        nb = set(A.neighbours(comp.id))
        return any(o.id in nb for o in others)
    if kind == "coincide":
        return any(o.bbox == comp.bbox and o.g == comp.g for o in others)
    if kind == "adjacent1":
        return any(abs(o.bbox[0] - comp.bbox[0]) + abs(o.bbox[1] - comp.bbox[1]) == 1 and o.g == comp.g for o in others)
    if kind == "samecolour":
        return any(o.color == comp.color for o in others)
    if kind == "sameshape":
        return any(o.key == comp.key for o in others)
    raise ValueError(kind)


# ------------------------------------------------------------------------------ prediction

@dataclass
class Fate:
    comp: Component
    kind: str                    # "identity" | "move" | "recolour" | "annihilate"
    cells: Optional[np.ndarray]  # new cells (None if annihilated or off-grid)
    color: int
    fill: Optional[int]
    off_grid: bool = False


def apply_clause(cl: Clause, comp: Component, A: BaseAlgebra) -> Fate:
    if cl.bind is not None and not _pair_predicate(cl.bind, comp, A):
        return Fate(comp, "identity", comp.cells, comp.color, None)
    g, (oy, ox), color = comp.g, comp.offset, comp.color
    moved = False
    for kind, param in cl.body:
        if kind in ("step", "move"):
            oy, ox = oy + param[0], ox + param[1]
            moved = True
        elif kind == "rot":
            g = d4_compose(param, g)
            moved = True
        elif kind == "recolour":
            color = int(param)
        elif kind == "copy":
            pass  # handled as a creation at the clause-set level (not in v1)
    if cl.term[0] == "annihilate":
        return Fate(comp, "annihilate", None, color, cl.fill)
    if moved:
        cells = transport_cells(comp, g, (oy, ox))
        if cells is None:
            return Fate(comp, "move", None, color, cl.fill, off_grid=True)
        return Fate(comp, "move", cells, color, cl.fill)
    if color != comp.color:
        return Fate(comp, "recolour", comp.cells, color, None)
    return Fate(comp, "identity", comp.cells, comp.color, None)


def predict(cs: ClauseSet, A: BaseAlgebra, data: Optional[dict]) -> tuple[list[Fate], list[Creation], str]:
    """Apply the clause set to the base algebra. Components claimed by several clauses: the first wins."""
    fates: dict[int, Fate] = {}
    for cl in cs.clauses:
        for comp in select(cl.selector, A, data):
            if comp.id in fates:
                continue
            fates[comp.id] = apply_clause(cl, comp, A)
    for comp in A.comps:
        if comp.id not in fates:
            fates[comp.id] = Fate(comp, "identity", comp.cells, comp.color, None)
    flag = "GAME_OVER" if cs.over else "NOT_FINISHED"
    return [fates[c.id] for c in A.comps], list(cs.creations), flag


def render(A: BaseAlgebra, fates: list[Fate], creations: list[Creation], data: Optional[dict] = None) -> Optional[np.ndarray]:
    """Vacate, then paint survivors in HN flag order (largest mass first, rarest last). None if any move left the grid."""
    canvas = A.frame.copy().ravel()
    for f in fates:
        if f.off_grid:
            return None
        if f.kind in ("move", "annihilate"):
            if f.fill is None:
                return None  # a vacating clause without a fitted fill cannot be rendered
            canvas[f.comp.cells] = f.fill
    survivors = [(A.hist[f.color] if f.color < 16 else 0, f.color, f.cells) for f in fates if f.kind != "annihilate" and f.cells is not None]
    for c in creations:
        cells = c.cells(A, data)
        if cells is None:
            return None
        survivors.append((A.hist[c.color], c.color, cells))
    survivors.sort(key=lambda t: -t[0])
    for _, color, cells in survivors:
        canvas[cells] = color
    return canvas.reshape(GRID, GRID)


def exact(cs: ClauseSet, A0: BaseAlgebra, data: Optional[dict], A1: BaseAlgebra, flag1: str) -> bool:
    fates, creations, flag = predict(cs, A0, data)
    if (flag == "GAME_OVER") != (flag1 == "GAME_OVER"):
        return False
    if flag1 == "GAME_OVER" and cs.over:
        return True  # the terminal transition's frame is not predicted beyond the flag
    out = render(A0, fates, creations, data)
    return out is not None and np.array_equal(out, A1.frame)


def predicted_frame(cs: ClauseSet, A: BaseAlgebra, data: Optional[dict]) -> Optional[np.ndarray]:
    fates, creations, _ = predict(cs, A, data)
    return render(A, fates, creations, data)


# ------------------------------------------------------------------------------ KMS weight

def clause_set_log_weight(cs: ClauseSet, lib) -> float:
    """log P of a clause set = sum over its clauses of the conditional Parry kernel of each path (Prop. 1.3 factorization)."""
    total = 0.0
    for letters in cs.letters():
        total += lib.log_P(letters)
    return total
