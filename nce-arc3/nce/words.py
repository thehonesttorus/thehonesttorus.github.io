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

from .algebra import BaseAlgebra, Component, GRID, _make_component, d4_apply, d4_inverse, transport_cells

VIRTUAL_ID = -1   # the pointed cell as a virtual component: the finest member of the SelPointed family

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

REVEAL = 16   # Fill family member: the remembered canvas content under the vacated cells (d = 17 with the palette)


MASS_GATES = ("massge", "massle")
SEL_NAME = {"all": "SelAll", "pointed": "SelPointed", "pointedcell": "SelPointed", "colour": "SelColour", "rank": "SelRank", "type": "SelType"}
BODY_NAME = {"step": "Step", "move": "Move", "rot": "Rot", "recolour": "Recolour", "copy": "Copy", "grow": "Grow", "erode": "Erode"}


@dataclass(frozen=True)
class Clause:
    selector: Letter
    bind: Optional[Letter]
    body: tuple[Letter, ...]
    term: Letter
    fill: Optional[int]
    top: bool = False      # Layer family (d = 2): paint the moved components above every stratum instead of in HN order
    tag: Optional[tuple] = field(default=None, compare=False)   # the behaviour class this clause was fitted to

    def letters(self) -> list[str]:
        """Letter family names along the Bratteli path (for KMS weights)."""
        names = [SEL_NAME[self.selector[0]]]
        if self.bind is not None:
            names += ["Meter" if self.bind[0] in MASS_GATES else "Bind", "Gate"]
        for kind, _ in self.body:
            names.append(BODY_NAME[kind])
        if self.fill is not None:
            names.append("Fill")
        if self.top:
            names.append("Layer")
        if self.term[0] == "annihilate":
            names.append("Annihilate")
        return names

    def lettered(self) -> list[tuple]:
        """The same path with each letter's fitted unit: (family name, parameter)."""
        out: list[tuple] = [(SEL_NAME[self.selector[0]], tuple(self.selector))]
        if self.bind is not None:
            out += [("Meter" if self.bind[0] in MASS_GATES else "Bind", tuple(self.bind)), ("Gate", None)]
        for kind, p in self.body:
            out.append((BODY_NAME[kind], (kind, p)))
        if self.fill is not None:
            out.append(("Fill", self.fill))
        if self.top:
            out.append(("Layer", True))
        if self.term[0] == "annihilate":
            out.append(("Annihilate", None))
        return out

    def vacates(self) -> bool:
        return self.term[0] == "annihilate" or any(k in ("step", "move", "rot", "erode") for k, _ in self.body)

    def describe(self) -> str:
        s = f"{self.selector[0]}" + (f"[{self.selector[1]}]" if len(self.selector) > 1 and self.selector[1] is not None and self.selector[0] != 'type' else "")
        if self.selector[0] == "type":
            s += "[type]"
        if self.bind:
            s += f"·gate({self.bind[0]})"
        for k, p in self.body:
            s += f"·{k}{p if k != 'rot' else '[' + str(p) + ']'}"
        if self.fill is not None:
            s += "·reveal" if self.fill == REVEAL else f"·fill{self.fill}"
        if self.top:
            s += "·top"
        if self.term[0] == "annihilate":
            s += "·annihilate"
        return s


@dataclass(frozen=True)
class Creation:
    key: bytes            # canonical shape key (must be in the alphabet)
    color: int
    g: int
    offset: tuple[int, int]      # absolute bbox corner, or offset relative to the anchor component when relative
    relative: bool = False       # Copy/placement relative to the pointed component (SelPointed . Copy . units)
    anchor: Optional[tuple] = None   # ("mover", selector) : placed relative to each component the clause selector picks
    gate: Optional[tuple] = None     # (predicate): the creation happens iff the pointed component satisfies it

    def letters(self) -> list[str]:
        if self.anchor is not None:
            return [SEL_NAME[self.anchor[1][0]], "Copy"]
        if self.gate is not None:
            # typed paths: F -SelPointed-> C -Bind-> P -Gate-> C -Copy-> C, or F -Create-> C -Bind-> P -Gate-> C
            return ["SelPointed", "Bind", "Gate", "Copy"] if self.relative else ["Create", "Bind", "Gate"]
        return ["SelPointed", "Copy"] if self.relative else ["Create"]

    def lettered(self) -> list[tuple]:
        names = self.letters()
        param = (self.key, self.color, self.g, tuple(self.offset))
        return [(n, None) for n in names[:-1]] + [(names[-1], param)]

    def cells_list(self, A: BaseAlgebra, data: Optional[dict]) -> Optional[list[np.ndarray]]:
        from .algebra import canonical_image_from_key, cells_from_canonical
        if self.gate is not None:
            if not data:
                return []
            cid = int(A.cell2comp[int(data["y"]), int(data["x"])])
            if cid < 0 or not gate_holds(self.gate, (), A.comps[cid], A):
                return []
        canon = canonical_image_from_key(self.key)
        if self.anchor is not None:
            bases = [c.offset for c in select(self.anchor[1], A, data)]
            if not bases:
                return []
        elif self.relative:
            if not data:
                return None
            cid = int(A.cell2comp[int(data["y"]), int(data["x"])])
            if cid < 0:
                return None
            bases = [A.comps[cid].offset]
        else:
            bases = [(0, 0)]
            canon_cells = cells_from_canonical(canon, self.g, self.offset)
            return None if canon_cells is None else [canon_cells]
        out = []
        for base in bases:
            cells = cells_from_canonical(canon, self.g, (base[0] + self.offset[0], base[1] + self.offset[1]))
            if cells is None:
                return None
            out.append(cells)
        return out


@dataclass(frozen=True)
class ClauseSet:
    clauses: tuple[Clause, ...] = ()
    creations: tuple[Creation, ...] = ()
    over: bool = False
    over_gate: Optional[tuple] = None    # (selector, gate, body): the action kills iff a selected component satisfies the gate
    level_gated: bool = False            # LevelGate family: the word is asserted for the current level only
    # The memorised part (T12): (pre-frame key, datum, post-frame bytes or b"GAME_OVER", residual cells) for the
    # demonstrations the generalising part does not reproduce. Admitted only under the code-length test in synth.
    exceptions: tuple = ()

    def letters(self) -> list[list[str]]:
        out = [c.letters() for c in self.clauses] + [c.letters() for c in self.creations]
        if self.over:
            out.append(["Over"])
        if self.over_gate is not None:
            out.append([SEL_NAME[self.over_gate[0][0]], "Meter" if self.over_gate[1][0] in MASS_GATES else "Bind", "Gate"])
            out.append(["Over"])
        if self.level_gated:
            out.append(["LevelGate"])
        return out

    def lettered(self) -> list[list[tuple]]:
        out = [c.lettered() for c in self.clauses] + [c.lettered() for c in self.creations]
        if self.over:
            out.append([("Over", None)])
        if self.over_gate is not None:
            out.append([(SEL_NAME[self.over_gate[0][0]], tuple(self.over_gate[0])),
                        ("Meter" if self.over_gate[1][0] in MASS_GATES else "Bind", tuple(self.over_gate[1])), ("Gate", None)])
            out.append([("Over", None)])
        if self.level_gated:
            out.append([("LevelGate", None)])
        return out

    def describe(self) -> str:
        parts = [c.describe() for c in self.clauses] + [f"create(c{c.color})" for c in self.creations]
        if self.over:
            parts.append("over")
        if self.over_gate is not None:
            parts.append(f"over if {self.over_gate[0][0]}:{self.over_gate[1][0]}")
        if self.level_gated:
            parts.append("levelgate")
        if self.exceptions:
            parts.append(f"+{len(self.exceptions)} memorised")
        return " ⊕ ".join(parts) if parts else "identity"


@dataclass
class Hypothesis:
    cs: ClauseSet
    logw: float
    n_verified: int = 0
    beyond_ceiling: bool = False
    contexts: set = field(default_factory=set)   # word contexts (see word_context) on which the word was verified

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
    if kind == "pointedcell":
        c = pointed_cell(A, data) if data else None
        return [c] if c is not None else []
    if kind == "colour":
        return [c for c in A.comps if c.color == selector[1]]
    if kind == "rank":
        return [c for c in A.comps if c.stratum == selector[1]]
    if kind == "type":
        return [c for c in A.comps if c.key == selector[1]]
    raise ValueError(kind)


def pointed_cell(A: BaseAlgebra, data: dict) -> Optional[Component]:
    """The pointed cell as a one-cell virtual component of its parent's colour and stratum (id VIRTUAL_ID)."""
    x, y = int(data["x"]), int(data["y"])
    cid = int(A.cell2comp[y, x])
    if cid < 0:
        return None
    parent = A.comps[cid]
    return _make_component(VIRTUAL_ID, parent.color, np.array([y * GRID + x], dtype=np.int64), parent.stratum)


def _pair_predicate(bind: Letter, comp: Component, A: BaseAlgebra) -> bool:
    kind, other_key = bind
    if kind == "massge":                 # the moment map's colour mass (an occupation number), a frame-level test
        return int(A.hist[other_key[0]]) >= other_key[1]
    if kind == "massle":
        return int(A.hist[other_key[0]]) <= other_key[1]
    if kind == "iscolour":
        return comp.color == other_key
    if kind == "istype":
        return comp.key == other_key
    if kind == "isrank":
        return comp.stratum == other_key
    others = [o for o in A.comps if o.id != comp.id and (other_key is None or o.key == other_key)]
    if kind == "contact":
        if comp.id < 0:
            return False
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


def gate_holds(gate: Letter, body: tuple, comp: Component, A: BaseAlgebra, data: Optional[dict] = None) -> bool:
    """Gates: a pair or unary predicate on the pre-state; a predicate on the transported component ("collide" with
    a component of a type, or any non-top component when the type is None; "offgrid"; colour collisions); a meter
    on a colour mass; or ("ptr", inner): the unary predicate `inner` applied to the pointed component (SelPointed
    composed into the gate: the datum decides, not the clause's own component)."""
    kind, key = gate
    if kind == "cofire":
        return True          # a clause-level condition, applied by predict(); every selected component passes here
    cache = A.__dict__.setdefault("_gatec", {})
    ck = (gate, body, comp.id if comp.id >= 0 else ("v", int(comp.cells[0])), _datum(data) if kind == "ptr" else None)
    v = cache.get(ck)
    if v is None:
        v = _gate_holds(gate, body, comp, A, data)
        cache[ck] = v
    return v


def _gate_holds(gate: Letter, body: tuple, comp: Component, A: BaseAlgebra, data: Optional[dict] = None) -> bool:
    kind, key = gate
    if kind == "ptr":
        if not data:
            return False
        cid = int(A.cell2comp[int(data["y"]), int(data["x"])])
        return cid >= 0 and _pair_predicate(key, A.comps[cid], A)
    if kind in ("notcollide", "notoffgrid", "notcollidecolour", "notinside"):
        return not gate_holds((kind[3:], key), body, comp, A, data)
    if kind not in ("collide", "offgrid", "collidecolour", "inside"):
        return _pair_predicate(gate, comp, A)
    g, (oy, ox) = comp.g, comp.offset
    for k, param in body:
        if k in ("step", "move"):
            oy, ox = oy + param[0], ox + param[1]
        elif k == "rot":
            g = d4_compose(param, g)
    cache = A.__dict__.setdefault("_transport", {})
    tkey = (comp.id if comp.id >= 0 else ("v", int(comp.cells[0])), g, oy, ox)
    if tkey in cache:
        cells = cache[tkey]
    else:
        cells = transport_cells(comp, g, (oy, ox))
        cache[tkey] = cells
    if cells is None:
        return kind == "offgrid"
    if kind == "offgrid":
        return False
    if kind in ("collidecolour", "inside"):
        new = np.setdiff1d(cells, comp.cells, assume_unique=True)
        under = A.frame.ravel()[new]
        return bool(np.any(under == key)) if kind == "collidecolour" else bool(np.all(under == key))
    occ = A.cell2comp.ravel()[cells]
    ids = set(occ[occ >= 0].tolist()) - {comp.id}
    top = len(A.strata) > 1
    for oid in ids:
        o = A.comps[oid]
        if key is None:
            if not (top and o.stratum == 0):
                return True
        elif o.key == key:
            return True
    return False


def grown_cells(cells: np.ndarray, d: tuple[int, int], c: int, frame: np.ndarray) -> np.ndarray:
    """Cells one step along d from `cells`, inside the grid, not in `cells`, showing colour c in `frame`."""
    ys, xs = np.divmod(cells, GRID)
    ys2, xs2 = ys + d[0], xs + d[1]
    ok = (ys2 >= 0) & (ys2 < GRID) & (xs2 >= 0) & (xs2 < GRID)
    sh = ys2[ok] * GRID + xs2[ok]
    sh = sh[frame.ravel()[sh] == c]
    return np.setdiff1d(sh, cells, assume_unique=False)


def eroded_cells(cells: np.ndarray, d: tuple[int, int]) -> np.ndarray:
    """P ^ T_d P: the cells of `cells` whose predecessor along d is also in `cells`."""
    ys, xs = np.divmod(cells, GRID)
    ys2, xs2 = ys - d[0], xs - d[1]
    ok = (ys2 >= 0) & (ys2 < GRID) & (xs2 >= 0) & (xs2 < GRID)
    prev = np.full(cells.shape, -1, dtype=np.int64)
    prev[ok] = ys2[ok] * GRID + xs2[ok]
    keep = np.isin(prev, cells)
    return cells[keep]


# ------------------------------------------------------------------------------ prediction

@dataclass
class Fate:
    comp: Component
    kind: str                    # "identity" | "move" | "recolour" | "annihilate"
    cells: Optional[np.ndarray]  # new cells (None if annihilated or off-grid)
    color: int
    fill: Optional[int]
    off_grid: bool = False
    top: bool = False


def apply_clause(cl: Clause, comp: Component, A: BaseAlgebra, data: Optional[dict] = None) -> Optional[Fate]:
    """The clause's fate for a component it selects; None if the gate excludes the component."""
    if cl.bind is not None and not gate_holds(cl.bind, cl.body, comp, A, data):
        return None
    g, (oy, ox), color = comp.g, comp.offset, comp.color
    moved = False
    cur: Optional[np.ndarray] = None      # cells under boundary letters (Grow, Erode), applied in sequence
    grew = eroded = False
    for kind, param in cl.body:
        if kind in ("step", "move"):
            oy, ox = oy + param[0], ox + param[1]
            moved = True
        elif kind == "rot":
            g = d4_compose(param, g)
            moved = True
        elif kind == "recolour":
            color = int(param)
        elif kind == "grow":
            # Grow[d, c]: P -> P v (T_d P ^ Q_c), the component advances one cell along d into cells of colour c
            (dy, dx), c = param
            base = comp.cells if cur is None else cur
            cur = np.union1d(base, grown_cells(base, (dy, dx), c, A.frame))
            grew = True
        elif kind == "erode":
            # Erode[d]: P -> P ^ T_d P, the component loses its trailing edge along d (the two-projection flow)
            cur = eroded_cells(comp.cells if cur is None else cur, param)
            eroded = True
        elif kind == "copy":
            pass  # handled as a creation at the clause-set level (not in v1)
    if cl.term[0] == "annihilate":
        return Fate(comp, "annihilate", None, color, cl.fill)
    if cur is not None:
        if cur.size == 0:
            return Fate(comp, "annihilate", None, color, cl.fill)
        if eroded:
            return Fate(comp, "move", cur, color, cl.fill, top=cl.top)
        return Fate(comp, "grow", cur, color, None, top=cl.top)
    if moved:
        cells = transport_cells(comp, g, (oy, ox))
        if cells is None:
            return Fate(comp, "move", None, color, cl.fill, off_grid=True)
        return Fate(comp, "move", cells, color, cl.fill, top=cl.top)
    if color != comp.color:
        return Fate(comp, "recolour", comp.cells, color, None)
    return Fate(comp, "identity", comp.cells, comp.color, None)


def predict(cs: ClauseSet, A: BaseAlgebra, data: Optional[dict]) -> tuple[list[Fate], list[Creation], str]:
    """Apply the clause set to the base algebra. Components claimed by several clauses: the first wins.
    Coupled clauses (gate ("cofire", class)) run in a second pass and fire iff the clause fitted to that class
    changed something in the first pass: a correspondence on a pair, through the library's pair vertex P."""
    fates: dict[int, Fate] = {}
    fired: set = set()
    coupled = [cl for cl in cs.clauses if cl.bind is not None and cl.bind[0] == "cofire"]
    for cl in cs.clauses:
        if cl.bind is not None and cl.bind[0] == "cofire":
            continue
        for comp in select(cl.selector, A, data):
            if comp.id in fates:
                continue
            f = apply_clause(cl, comp, A, data)
            if f is not None:
                fates[comp.id] = f
                if f.kind != "identity" and cl.tag is not None:
                    fired.add(cl.tag)
    for cl in coupled:
        if cl.bind[1] not in fired:
            continue
        for comp in select(cl.selector, A, data):
            if comp.id in fates:
                continue
            f = apply_clause(cl, comp, A, data)
            if f is not None:
                fates[comp.id] = f
    for comp in A.comps:
        if comp.id not in fates:
            fates[comp.id] = Fate(comp, "identity", comp.cells, comp.color, None)
    dead = cs.over
    if not dead and cs.over_gate is not None:
        sel, gate, body = cs.over_gate
        dead = any(gate_holds(gate, body, comp, A, data) for comp in select(sel, A, data))
    memo = memorised(cs, A, data)
    if memo is not None:
        dead = memo == b"GAME_OVER"
    flag = "GAME_OVER" if dead else "NOT_FINISHED"
    return [fates[c.id] for c in A.comps] + [f for k, f in fates.items() if k < 0], list(cs.creations), flag


def render(A: BaseAlgebra, fates: list[Fate], creations: list[Creation], data: Optional[dict] = None,
           memory: Optional[np.ndarray] = None) -> Optional[np.ndarray]:
    """Vacate (constant fill or remembered canvas), then paint survivors in HN flag order (largest mass
    first, rarest last), movers flagged `top` last. None if any move left the grid or a fill is missing."""
    canvas = A.frame.astype(np.int16).ravel()

    def vacate(f: Fate) -> bool:
        if f.fill is None:
            return False  # a vacating clause without a fitted fill cannot be rendered
        if f.fill == REVEAL:
            if memory is None:
                canvas[f.comp.cells] = -1          # unknown canvas: wildcard cells
            else:
                canvas[f.comp.cells] = memory.ravel()[f.comp.cells]   # -1 where never observed
        else:
            canvas[f.comp.cells] = f.fill
        return True

    late: list[Fate] = []     # virtual (sub-component) fates: painted after their parents
    for f in fates:
        if f.off_grid:
            return None
        if f.comp.id < 0:
            late.append(f)
            continue
        if f.kind in ("move", "annihilate") and not vacate(f):
            return None
    survivors = [(1 if f.top else 0, A.hist[f.color] if f.color < 16 else 0, f.color, f.cells) for f in fates
                 if f.comp.id >= 0 and f.kind != "annihilate" and f.cells is not None]
    for c in creations:
        cells_list = c.cells_list(A, data)
        if cells_list is None:
            return None
        for cells in cells_list:
            survivors.append((0, A.hist[c.color], c.color, cells))
    survivors.sort(key=lambda t: (t[0], -t[1]))
    for _, _, color, cells in survivors:
        canvas[cells] = color
    for f in late:
        if f.kind in ("move", "annihilate") and not vacate(f):
            return None
    for f in late:
        if f.kind != "annihilate" and f.cells is not None:
            canvas[f.cells] = f.color
    return canvas.reshape(GRID, GRID)


def frame_matches(pred: Optional[np.ndarray], frame: np.ndarray, ignore: Optional[frozenset] = None) -> bool:
    """Cell-for-cell equality, with -1 cells of the prediction (unknown canvas under a reveal) as wildcards.
    Cells showing a colour in `ignore` (the level's clock colours) before or after are outside the comparison:
    the exact test is taken in the quotient by the autonomous flow."""
    if pred is None:
        return False
    known = pred >= 0
    if ignore:
        for c in ignore:
            known &= (pred != c) & (frame != c)
    return bool(np.array_equal(pred[known], frame[known].astype(np.int16)))


def concretize(pred: Optional[np.ndarray], A: BaseAlgebra) -> Optional[np.ndarray]:
    """Replace wildcard cells by the top-stratum colour for Kato/planning renderings (the best available guess)."""
    if pred is None:
        return None
    if np.any(pred < 0):
        pred = pred.copy()
        pred[pred < 0] = A.strata[0][0] if A.strata else 0
    return pred.astype(np.int8)


MEMO_CELL = float(np.log(16 * GRID * GRID))   # code length of one memorised cell: a colour and a lattice site


def _datum(data: Optional[dict]) -> Optional[tuple]:
    return (int(data["x"]), int(data["y"])) if data else None


def memorised(cs: ClauseSet, A: BaseAlgebra, data: Optional[dict]) -> Optional[bytes]:
    """The memorised post-state for exactly this (pre-frame, datum), if the clause set holds one."""
    if not cs.exceptions:
        return None
    key, d = A.key(), _datum(data)
    for k, dd, post, _ in cs.exceptions:
        if k == key and dd == d:
            return post
    return None


def exact(cs: ClauseSet, A0: BaseAlgebra, data: Optional[dict], A1: BaseAlgebra, flag1: str, memory: Optional[np.ndarray] = None,
          ignore: Optional[frozenset] = None) -> bool:
    memo = memorised(cs, A0, data)
    if memo is not None:
        if memo == b"GAME_OVER" or flag1 == "GAME_OVER":
            return (memo == b"GAME_OVER") == (flag1 == "GAME_OVER")
        return memo == np.ascontiguousarray(A1.frame, dtype=np.int8).tobytes()
    fates, creations, flag = predict(cs, A0, data)
    if (flag == "GAME_OVER") != (flag1 == "GAME_OVER"):
        return False
    if flag1 == "GAME_OVER":
        return True  # the terminal transition's frame is not predicted beyond the flag
    out = render(A0, fates, creations, data, memory)
    return frame_matches(out, A1.frame, ignore)


def predicted_frame(cs: ClauseSet, A: BaseAlgebra, data: Optional[dict], memory: Optional[np.ndarray] = None) -> Optional[np.ndarray]:
    memo = memorised(cs, A, data)
    if memo is not None:
        return None if memo == b"GAME_OVER" else np.frombuffer(memo, dtype=np.int8).reshape(GRID, GRID).copy()
    fates, creations, _ = predict(cs, A, data)
    return concretize(render(A, fates, creations, data, memory), A)


def word_context(cs: ClauseSet, A: BaseAlgebra, data: Optional[dict]) -> tuple:
    """What the word sees in a situation: per clause the multiset of (component type, gate outcome) over its
    selection, and the same for the death gate. The newborn rule (T2, T7) is relative to the word: a situation
    is a new experiment for a word iff the word's own context there is one it has never been verified on."""
    out = []
    for cl in cs.clauses:
        sel = select(cl.selector, A, data)
        out.append(tuple(sorted((c.type, True if cl.bind is None else gate_holds(cl.bind, cl.body, c, A, data)) for c in sel)))
    if cs.over_gate is not None:
        sel_, gate, body = cs.over_gate
        out.append(tuple(sorted((c.type, gate_holds(gate, body, c, A, data)) for c in select(sel_, A, data))))
    return tuple(out)


def predicted(cs: ClauseSet, A: BaseAlgebra, data: Optional[dict], memory: Optional[np.ndarray] = None) -> tuple[Optional[np.ndarray], str]:
    """(predicted frame or None, predicted flag). A GAME_OVER prediction carries no frame."""
    memo = memorised(cs, A, data)
    if memo is not None:
        if memo == b"GAME_OVER":
            return None, "GAME_OVER"
        return np.frombuffer(memo, dtype=np.int8).reshape(GRID, GRID).copy(), "NOT_FINISHED"
    fates, creations, flag = predict(cs, A, data)
    if flag == "GAME_OVER":
        return None, flag
    return concretize(render(A, fates, creations, data, memory), A), flag


# ------------------------------------------------------------------------------ KMS weight

def clause_set_log_weight(cs: ClauseSet, lib) -> float:
    """log P of a clause set = sum over its clauses of the conditional Parry kernel of each path (Prop. 1.3 factorization)."""
    total = 0.0
    for letters in cs.letters():
        total += lib.log_P(letters)
    # the memorised part pays its code length in the same units: beta_c log(16 * 4096) per residual cell
    total -= lib.beta_c * MEMO_CELL * float(sum(e[3] for e in cs.exceptions))
    return total
