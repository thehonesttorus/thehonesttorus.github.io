"""Synthesis of clause sets in Parry-KMS order with exact verification (Sections 4.5-4.6).

The enumeration is residual-guided (kms-library graft): the stored transitions of an action
are decomposed into behaviour classes (components that moved by the same gauge element,
recoloured to the same colour, were annihilated), and for each class the candidate clauses
are those whose selector (and optional pair gate) separates exactly the components showing
that behaviour from the ones that do not, in every stored transition, with one parameter
vector shared across all of them (C11: a linear pass over the components' invariants, never
a search over frames). Candidate clause sets are ranked by the conditional Parry kernel at
beta_c and accepted only on cell-for-cell reproduction of every stored transition (T6).
Every exact clause set found is kept (T7). Partially explaining clauses are returned as
partial words: they may drive the Kato score (C10) but never stand as the action's law.
"""
from __future__ import annotations

import itertools
import time
from dataclasses import dataclass, field, replace
from typing import Optional

import numpy as np

from .algebra import GRID, BaseAlgebra, Component
from .library import Library
from .store import Record, Store
from .words import REVEAL, VIRTUAL_ID, Clause, ClauseSet, Creation, Hypothesis, clause_set_log_weight, d4_relative, eroded_cells, exact, gate_holds, grown_cells, pointed_cell, predicted_frame, select, word_context, _pair_predicate

K_MAX = 32
MAX_STEP_WORD = 6          # Step^k beats Move for k <= 6 under the Parry order (Section 0, item 1)
MAX_COMBOS = 256


@dataclass
class SynthResult:
    hypotheses: list[Hypothesis] = field(default_factory=list)
    partials: list[ClauseSet] = field(default_factory=list)     # exact for their class, not for the whole transition
    classes: int = 0
    unexplained_classes: int = 0
    exhausted: bool = False
    # count-only diagnostics
    reshaped_records: int = 0        # records with a reshaped frame vector (same colour, new key)
    fill_inconsistent: int = 0       # classes killed by an inconsistent vacated-cell colour
    creation_failed: bool = False    # creations neither constant nor pointed-relative
    combos_tried: int = 0
    min_mismatch_cells: int = -1     # best combination's cell mismatch (0 when exact)
    level_gated: bool = False
    death_gated: bool = False
    clock: frozenset = frozenset()   # colours verified modulo (the level's autonomous flow)
    death_unexplained: bool = False  # terminal records exist but no gate separates them from the live ones


# ------------------------------------------------------------------ behaviour classes

def _signature(rec: Record, cid: int, use_rot: bool = False) -> Optional[tuple]:
    """Outcome of component cid in record rec: None for identity/reshape, else a hashable signature.
    With use_rot, a reshape that is also congruent under a D4 change reads as an in-place Rot."""
    tr = rec.tr
    for i, j, g in tr.moves:
        if i == cid:
            if g.is_identity:
                return None
            h = d4_relative(g.g_before, g.g_after) if g.rotated else 0
            return ("move", g.dy, g.dx, h)
    if use_rot:
        for i, j, g in tr.inplace_rot:
            if i == cid:
                return ("move", 0, 0, d4_relative(g.g_before, g.g_after))
    for i, j in tr.recolours:
        if i == cid:
            return ("recolour", rec.A1.comps[j].color)
    if cid in tr.annihilated:
        return ("annihilate",)
    for i, j in tr.reshaped:
        if i == cid:
            a, b = rec.A0.comps[i], rec.A1.comps[j]
            if a.stratum == 0 and len(rec.A0.strata) > 1:
                return None
            if b.n < a.n and np.all(np.isin(b.cells, a.cells, assume_unique=True)):
                for d in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    if np.array_equal(eroded_cells(a.cells, d), b.cells):
                        return ("erode", d)
                return None
            if b.n <= a.n or not np.all(np.isin(a.cells, b.cells, assume_unique=True)):
                return None
            new = np.setdiff1d(b.cells, a.cells, assume_unique=True)
            cols = np.unique(rec.A0.frame.ravel()[new])
            if cols.size != 1:
                return None
            c = int(cols[0])
            for d in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                if np.array_equal(np.sort(grown_cells(a.cells, d, c, rec.A0.frame)), new):
                    return ("grow", d, c)
            return None
    return None


def _cell_signature(rec: Record) -> Optional[tuple]:
    """Outcome of the pointed cell when its parent component has no outcome of its own: a sub-component recolour."""
    if not rec.data:
        return None
    x, y = int(rec.data["x"]), int(rec.data["y"])
    cid = int(rec.A0.cell2comp[y, x])
    if cid < 0 or _signature(rec, cid, True) is not None:
        return None
    c0, c1 = int(rec.A0.frame[y, x]), int(rec.A1.frame[y, x])
    if c0 == c1:
        return None
    return ("cellrecolour", c1)


def _body_for(sig: tuple) -> tuple[tuple, ...]:
    if sig[0] == "move":
        _, dy, dx, h = sig
        body: list[tuple] = []
        if h:
            body.append(("rot", h))
        l1 = abs(dy) + abs(dx)
        if l1 > 0:
            if l1 <= MAX_STEP_WORD:
                body += [("step", (int(np.sign(dy)), 0))] * abs(dy) + [("step", (0, int(np.sign(dx))))] * abs(dx)
            else:
                body.append(("move", (dy, dx)))
        return tuple(body)
    if sig[0] in ("recolour", "cellrecolour"):
        return (("recolour", sig[1]),)
    if sig[0] == "grow":
        return (("grow", (sig[1], sig[2])),)
    if sig[0] == "erode":
        return (("erode", sig[1]),)
    return ()


def _fill_admissible(rec: Record, cids: list[int], sig: Optional[tuple] = None) -> Optional[set[int]]:
    """Admissible Fill members for the cells vacated by components cids in rec: a constant palette colour
    if the vacated cells show one colour; REVEAL if they show the remembered canvas wherever it is known
    (an unknown canvas cannot falsify a reveal). None when nothing was vacated (unconstrained)."""
    A0, A1, tr = rec.A0, rec.A1, rec.tr
    new_cover = np.zeros(64 * 64, dtype=bool)
    for i, j, g in tr.moves:
        if not g.is_identity:
            new_cover[A1.comps[j].cells] = True
    for i, j, g in tr.inplace_rot:
        new_cover[A1.comps[j].cells] = True
    vac_all = []
    for cid in cids:
        cells = A0.comps[cid].cells
        if sig is not None and sig[0] == "erode":
            cells = np.setdiff1d(cells, eroded_cells(cells, sig[1]), assume_unique=True)
        vac = cells[~new_cover[cells]]
        if vac.size:
            vac_all.append(vac)
    if not vac_all:
        return None
    vac = np.concatenate(vac_all)
    shown = A1.frame.ravel()[vac].astype(np.int16)
    adm: set[int] = set()
    colours = np.unique(shown)
    if colours.size == 1:
        adm.add(int(colours[0]))
    if rec.memory is None:
        adm.add(REVEAL)
    else:
        mem = rec.memory.ravel()[vac]
        known = mem >= 0
        if np.array_equal(mem[known], shown[known]):
            adm.add(REVEAL)
    return adm


# ------------------------------------------------------------------ selector fitting

def _selector_candidates(recs: list[Record]) -> list[tuple]:
    cands: list[tuple] = [("all",), ("pointed",)]
    colours, ranks, keys = set(), set(), set()
    for r in recs:
        for c in r.A0.comps:
            colours.add(c.color); ranks.add(c.stratum); keys.add(c.key)
    if any(r.data for r in recs):
        cands.append(("pointedcell",))
    cands += [("colour", c) for c in sorted(colours)]
    cands += [("rank", k) for k in sorted(ranks)]
    cands += [("type", k) for k in sorted(keys)]
    return cands


def _predicate_candidates(recs: list[Record]) -> list[tuple]:
    keys = sorted({c.key for r in recs for c in r.A0.comps})
    preds: list[tuple] = []
    for kind in ("contact", "coincide", "adjacent1"):
        preds.append((kind, None))
        preds += [(kind, k) for k in keys]
    preds += [("samecolour", None), ("sameshape", None)]
    preds += [("iscolour", c) for c in sorted({c.color for r in recs for c in r.A0.comps})]
    preds += [("istype", k) for k in keys]
    preds += [("isrank", s) for s in sorted({c.stratum for r in recs for c in r.A0.comps})]
    return preds


class _Selections:
    """Memo of selector applications per record (the fit loops apply the same selector many times)."""

    def __init__(self, recs: list[Record]):
        self.recs = recs
        self.cache: dict[tuple, list] = {}

    def get(self, selector: tuple, ri: int) -> list:
        key = (selector, ri)
        v = self.cache.get(key)
        if v is None:
            r = self.recs[ri]
            v = select(selector, r.A0, r.data)
            self.cache[key] = v
        return v


def _fits_exactly(selector: tuple, recs: list[Record], targets: list[set[int]], memo: Optional[_Selections] = None) -> bool:
    memo = memo or _Selections(recs)
    for ri, B in enumerate(targets):
        if {c.id for c in memo.get(selector, ri)} != B:
            return False
    return True


def _is_superset(selector: tuple, recs: list[Record], targets: list[set[int]], memo: _Selections) -> bool:
    for ri, B in enumerate(targets):
        if B and not B <= {c.id for c in memo.get(selector, ri)}:
            return False
    return True


class _GateTable(dict):
    """gate_holds(pred, body, comp, A0) per component id, computed on first access."""

    def __init__(self, rec: Record, pred: tuple, body: tuple):
        super().__init__()
        self.rec, self.pred, self.body = rec, pred, body

    def __missing__(self, cid: int) -> bool:
        A = self.rec.A0
        comp = A.comps[cid] if cid >= 0 else pointed_cell(A, self.rec.data)
        v = bool(gate_holds(self.pred, self.body, comp, A)) if comp is not None else False
        self[cid] = v
        return v


def _gate_table(rec: Record, pred: tuple, body: tuple) -> _GateTable:
    """Lazy per-record table of a gate, shared across classes for body-independent predicates (pair, unary and
    mass tests) and per body for the collision tests."""
    dep = pred[0] in ("collide", "offgrid", "collidecolour", "inside", "notcollide", "notoffgrid", "notcollidecolour", "notinside")
    key = (pred, body if dep else ())
    cache = rec.__dict__.setdefault("_gates", {})
    t = cache.get(key)
    if t is None:
        t = _GateTable(rec, pred, body)
        cache[key] = t
    return t


def _mass_gate_candidates(pos: list[Record], neg: list[Record]) -> list[tuple]:
    """Meter gates separating two sets of records exactly by one colour mass (an occupation number of the moment
    map): mass_c >= t on every positive and < t on every negative, or the mirror. t is the positives' extreme
    value, i.e. the gate claims nothing outside the masses it was demonstrated on (the support of the data)."""
    out: list[tuple] = []
    if not pos or not neg:
        return out
    for c in range(16):
        p = [int(r.A0.hist[c]) for r in pos]
        n = [int(r.A0.hist[c]) for r in neg]
        if min(p) > max(n):
            out.append(("massge", (c, min(p))))
        elif max(p) < min(n):
            out.append(("massle", (c, max(p))))
    return out


def _fits_with_gate(selector: tuple, pred: tuple, recs: list[Record], targets: list[set[int]], body: tuple = (),
                    memo: Optional[_Selections] = None) -> bool:
    memo = memo or _Selections(recs)
    for ri, B in enumerate(targets):
        r = recs[ri]
        selected = memo.get(selector, ri)
        if not B <= {c.id for c in selected}:
            return False
        t = _gate_table(r, pred, body)
        for comp in selected:
            if t[comp.id] != (comp.id in B):
                return False
    return True


def _holds_on_positives(pred: tuple, body: tuple, recs: list[Record], targets: list[set[int]]) -> bool:
    for ri, B in enumerate(targets):
        if not B:
            continue
        t = _gate_table(recs[ri], pred, body)
        if not all(t[cid] for cid in B):
            return False
    return True


# ------------------------------------------------------------------ main entry

def synthesize(a: int, store: Store, lib: Library, budget_s: float = 0.25, k_max: int = K_MAX) -> SynthResult:
    t0 = time.time()
    res = SynthResult()
    recs_all = store.transitions(a)
    recs = [r for r in recs_all if r.flag1 != "GAME_OVER"]
    terminal = [r for r in recs_all if r.flag1 == "GAME_OVER"]
    res.clock = store.clock_colours()
    if not recs:
        if terminal:
            cs = ClauseSet(over=True)
            res.hypotheses.append(Hypothesis(cs, clause_set_log_weight(cs, lib), len(terminal)))
        res.exhausted = True
        return res

    # 1. behaviour classes across all non-terminal records. Readings are tried in Parry order: first every
    #    congruent in-place D4 change as an identity with a changed rendering (cost 0), then as a Rot (cost log 7).
    readings = _readings(recs, lib)
    for use_rot, short in readings:
        _synthesize_reading(a, recs, terminal, lib, t0, budget_s, k_max, res, use_rot, short)
        if res.hypotheses:
            break
    # LevelGate (T9, T13): when no word covers every level, the current level's records alone, at the
    # cost of a LevelGate letter (d = win_levels). Earlier levels are kept as data, never discarded.
    if not res.hypotheses and len({r.level for r in recs}) > 1:
        cur = max(r.level for r in recs)
        recs = [r for r in recs if r.level == cur]
        terminal = [r for r in terminal if r.level == cur]
        for use_rot, short in readings:
            _synthesize_reading(a, recs, terminal, lib, t0, budget_s, k_max, res, use_rot, short)
            if res.hypotheses:
                break
        for h in res.hypotheses:
            h.cs = replace(h.cs, level_gated=True)
            h.logw = clause_set_log_weight(h.cs, lib)
        res.level_gated = bool(res.hypotheses)
    # death gate (Over . Gate . Bind . Sel): the action kills iff a selected component satisfies the gate
    if res.hypotheses and terminal:
        _fit_death_gate(recs, terminal, res, lib, t0, budget_s)
    for h in res.hypotheses:
        h.contexts = {word_context(h.cs, r.A0, r.data) for r in recs + terminal}
    return res


def _fit_death_gate(live: list[Record], terminal: list[Record], res: SynthResult, lib: Library, t0: float, budget_s: float) -> None:
    allr = live + terminal
    selectors = _selector_candidates(allr)
    preds = _predicate_candidates(allr)
    keys = sorted({c.key for r in allr for c in r.A0.comps})
    bodies = sorted({cl.body for h in res.hypotheses for cl in h.cs.clauses if cl.body}, key=len)
    gates: list[tuple[tuple, tuple]] = [(p, ()) for p in preds]
    gates += [(g_, ()) for g_ in _mass_gate_candidates(terminal, live)]     # death when a colour mass crosses a value
    for b in bodies:
        gates.append((("offgrid", None), b))
        gates.append((("collide", None), b))
        gates += [(("collide", k), b) for k in keys]
    for sel in selectors:
        for gate, body in gates:
            if time.time() - t0 > budget_s:
                return
            if all(any(gate_holds(gate, body, c, r.A0) for c in select(sel, r.A0, r.data)) for r in terminal) and \
               not any(any(gate_holds(gate, body, c, r.A0) for c in select(sel, r.A0, r.data)) for r in live):
                for h in res.hypotheses:
                    h.cs = replace(h.cs, over_gate=(sel, gate, body))
                    h.logw = clause_set_log_weight(h.cs, lib)
                res.death_gated = True
                return
    res.death_unexplained = True


def _short_shift(lib: Library) -> int:
    """The translation length beyond which the Parry order prefers annihilation plus a pointed copy to a move:
    4 l > d_Annihilate * d_SelPointed * d_Copy (log weights), so l* = that product / 4 (no free constant)."""
    return max(1, int(lib.dims["Annihilate"] * lib.dims["SelPointed"] * lib.dims["Copy"]) // 4)


def _readings(recs: list[Record], lib: Library) -> list[tuple[bool, bool]]:
    """Readings of the demonstrations in Parry order: (identity, Rot) x (long moves as moves, as annihilate + copy).
    The second reading of each pair exists only when some record makes it differ."""
    rot = any(r.tr.inplace_rot for r in recs)
    short = any(r.data for r in recs) and any(g.l1 > _short_shift(lib) for r in recs for _, _, g in r.tr.moved)
    out = [(False, False)]
    if rot:
        out.append((True, False))
    if short:
        out.append((False, True))
        if rot:
            out.append((True, True))
    return out


def _with_reading(recs: list[Record], lib: Library, short: bool):
    """Context manager swapping each record's decomposition for its short-shift reading (cached on the record)."""
    class _Ctx:
        def __enter__(self_):
            self_.saved = []
            if short:
                from .match import decompose
                ms = _short_shift(lib)
                for r in recs:
                    alt = r.__dict__.get("_tr_short")
                    if alt is None:
                        alt = decompose(r.A0, r.A1, max_shift=ms)
                        r.__dict__["_tr_short"] = alt
                    self_.saved.append((r, r.tr))
                    r.tr = alt
            return self_

        def __exit__(self_, *exc):
            for r, tr in self_.saved:
                r.tr = tr
            return False
    return _Ctx()


def _synthesize_reading(a: int, recs: list[Record], terminal: list[Record], lib: Library, t0: float, budget_s: float,
                        k_max: int, res: SynthResult, use_rot: bool, short: bool = False) -> None:
    with _with_reading(recs, lib, short):
        _synthesize_reading_inner(a, recs, terminal, lib, t0, budget_s, k_max, res, use_rot)


def _synthesize_reading_inner(a: int, recs: list[Record], terminal: list[Record], lib: Library, t0: float, budget_s: float,
                              k_max: int, res: SynthResult, use_rot: bool) -> None:
    classes: dict[tuple, list[set[int]]] = {}
    for ri, r in enumerate(recs):
        for c in r.A0.comps:
            if c.color in res.clock:
                continue                      # the autonomous flow is not a word of this action
            sig = _signature(r, c.id, use_rot)
            if sig is None:
                continue
            classes.setdefault(sig, [set() for _ in recs])[ri].add(c.id)
        sig = _cell_signature(r)
        if sig is not None:
            classes.setdefault(sig, [set() for _ in recs])[ri].add(VIRTUAL_ID)
    res.classes = len(classes)
    res.unexplained_classes = 0
    res.fill_inconsistent = 0
    res.partials = []
    res.reshaped_records = sum(1 for r in recs if any(r.A0.comps[i].stratum > 0 for i, _ in r.tr.reshaped))

    # creations: constant (Create[tau,g,t]) or relative to the pointed component (SelPointed . Copy . units).
    # A "created" component that lies entirely inside cells vacated this transition and shows the remembered
    # canvas is not a creation: it is revealed by the fill (occlusion), and is excluded here.
    def _unexplained_created(r: Record) -> list[int]:
        out = []
        vacated = np.zeros(64 * 64, dtype=bool)
        for i, j, g in r.tr.moves:
            if not g.is_identity:
                vacated[r.A0.comps[i].cells] = True
        for i in r.tr.annihilated:
            vacated[r.A0.comps[i].cells] = True
        for j in r.tr.created:
            cells = r.A1.comps[j].cells
            if r.A1.comps[j].color in res.clock:
                continue                      # the autonomous flow is not a creation of this action
            if r.memory is not None and np.all(vacated[cells]) and np.all(r.memory.ravel()[cells] == r.A1.frame.ravel()[cells]):
                continue
            out.append(j)
        return out
    creations: list[Creation] = []
    created_lists = [_unexplained_created(r) for r in recs]
    created_sets = [tuple(sorted((r.A1.comps[j].key, r.A1.comps[j].color, r.A1.comps[j].g, r.A1.comps[j].offset) for j in cl)) for r, cl in zip(recs, created_lists)]
    creation_ok = all(cs == created_sets[0] for cs in created_sets)
    if creation_ok:
        creations = [Creation(k, col, g, off) for (k, col, g, off) in created_sets[0]]
    elif all(r.data for r in recs):
        rel_sets = []
        for r, cl in zip(recs, created_lists):
            cid = int(r.A0.cell2comp[int(r.data["y"]), int(r.data["x"])])
            if cid < 0:
                rel_sets = None; break
            base = r.A0.comps[cid].offset
            rel_sets.append(tuple(sorted((r.A1.comps[j].key, r.A1.comps[j].color, r.A1.comps[j].g, (r.A1.comps[j].offset[0] - base[0], r.A1.comps[j].offset[1] - base[1])) for j in cl)))
        if rel_sets and all(cs == rel_sets[0] for cs in rel_sets):
            creation_ok = True
            creations = [Creation(k, col, g, off, relative=True) for (k, col, g, off) in rel_sets[0]]
        elif rel_sets:
            # gated creation: the same pointed-relative (or constant) set in every record where something was created,
            # nothing elsewhere, and a predicate on the pointed component separating the two groups exactly
            nonempty = [i for i, cs in enumerate(rel_sets) if cs]
            if nonempty and all(rel_sets[i] == rel_sets[nonempty[0]] for i in nonempty):
                pointed = [r.A0.comps[int(r.A0.cell2comp[int(r.data["y"]), int(r.data["x"])])] for r in recs]
                want = [bool(cs) for cs in rel_sets]
                for pred in _predicate_candidates(recs):
                    if all(gate_holds(pred, (), p, r.A0) == w for p, r, w in zip(pointed, recs, want)):
                        creation_ok = True
                        spec = rel_sets[nonempty[0]]
                        const = all(created_sets[i] == created_sets[nonempty[0]] for i in nonempty)
                        if const:
                            creations = [Creation(k, col, g, off, gate=pred) for (k, col, g, off) in created_sets[nonempty[0]]]
                        else:
                            creations = [Creation(k, col, g, off, relative=True, gate=pred) for (k, col, g, off) in spec]
                        break
    if not creation_ok and any(created_lists):
        # a copy left relative to each moved component's old place (SelX . Copy . units): one creation per mover
        mover_sets = []
        for r, cl in zip(recs, created_lists):
            movers = [(i, r.A0.comps[i]) for i, j, g in r.tr.moved]
            created = [r.A1.comps[j] for j in cl]
            if not movers or len(created) != len(movers):
                mover_sets = None; break
            assigned = []
            used = set()
            for i, m in movers:
                best = None
                for cj in created:
                    if cj.id in used:
                        continue
                    rel = (cj.offset[0] - m.offset[0], cj.offset[1] - m.offset[1])
                    if best is None or abs(rel[0]) + abs(rel[1]) < abs(best[1][0]) + abs(best[1][1]):
                        best = (cj, rel)
                if best is None:
                    mover_sets = None; break
                used.add(best[0].id)
                assigned.append((best[0].key, best[0].color, best[0].g, best[1], m.key, m.color))
            if mover_sets is None:
                break
            mover_sets.append(tuple(sorted(assigned)))
        if mover_sets and all(ms == mover_sets[0] for ms in mover_sets):
            spec = {(k, col, g, rel) for (k, col, g, rel, mk, mc) in mover_sets[0]}
            if len(spec) == 1:
                k, col, g, rel = next(iter(spec))
                mover_colours = {mc for (_, _, _, _, mk, mc) in mover_sets[0]}
                anchor_sel = ("colour", next(iter(mover_colours))) if len(mover_colours) == 1 else ("all",)
                creation_ok = True
                creations = [Creation(k, col, g, rel, anchor=("mover", anchor_sel))]

    selectors = _selector_candidates(recs)
    predicates = _predicate_candidates(recs)
    memo = _Selections(recs)

    # 2. per class: clauses whose selector/gate separates exactly, with constant fill
    alternatives: list[list[Clause]] = []
    for sig, targets in classes.items():
        body = _body_for(sig)
        term = ("annihilate",) if sig[0] == "annihilate" else ("stay",)
        vacates = sig[0] not in ("recolour", "cellrecolour", "grow")
        fills_ok: list[int] = []
        if vacates:
            common: Optional[set[int]] = None
            for r, B in zip(recs, targets):
                if not B:
                    continue
                adm = _fill_admissible(r, sorted(B), sig)
                if adm is None:
                    continue
                common = adm if common is None else (common & adm)
            if common is None:
                fills_ok = [int(recs[0].A0.strata[0][0])]        # nothing vacated anywhere: fill never rendered
            elif not common:
                res.unexplained_classes += 1
                res.fill_inconsistent += 1
                continue
            else:
                fills_ok = sorted(common, key=lambda c: (c != REVEAL, c))   # reveal first: it is the canvas, not a guess
        fill: Optional[int] = fills_ok[0] if vacates else None
        found: list[Clause] = []
        for selector in selectors:
            if time.time() - t0 > budget_s:
                break
            if _fits_exactly(selector, recs, targets, memo):
                found.append(Clause(selector, None, body, term, fill))
        if not found:
            preds = predicates
            if body and any(k in ("step", "move", "rot") for k, _ in body):
                # a move that sometimes does not happen: gated on what the transported component would meet
                keys = sorted({c.key for r in recs for c in r.A0.comps})
                cols = sorted({c.color for r in recs for c in r.A0.comps})
                preds = preds + [("notcollide", None), ("notoffgrid", None)] + [("notcollide", k) for k in keys] \
                        + [("inside", c) for c in cols] + [("notcollidecolour", c) for c in cols] \
                        + [("collide", None), ("offgrid", None)] + [("collide", k) for k in keys] \
                        + [("collidecolour", c) for c in cols] + [("notinside", c) for c in cols]
            supers = [s_ for s_ in selectors if _is_superset(s_, recs, targets, memo)]
            if supers:
                kept = []
                for p_ in preds:                      # a gate must admit every positive
                    if time.time() - t0 > budget_s:
                        break
                    if _holds_on_positives(p_, body, recs, targets):
                        kept.append(p_)
                preds = kept
            for selector in supers:
                # frame-level meter gates: records where the class acts versus records where a selected component does not
                pos = [r for ri, r in enumerate(recs) if targets[ri]]
                neg = [r for ri, r in enumerate(recs) if {c.id for c in memo.get(selector, ri)} - targets[ri]]
                if pos and neg and not (set(map(id, pos)) & set(map(id, neg))):
                    for g_ in _mass_gate_candidates(pos, neg):
                        if _fits_with_gate(selector, g_, recs, targets, body, memo):
                            found.append(Clause(selector, g_, body, term, fill))
                for pred in preds:
                    if time.time() - t0 > budget_s:
                        break
                    if _fits_with_gate(selector, pred, recs, targets, body, memo):
                        found.append(Clause(selector, pred, body, term, fill))
        if found and vacates and len(fills_ok) > 1:
            found = found + [Clause(cl.selector, cl.bind, cl.body, cl.term, f) for cl in found for f in fills_ok[1:]]
        if found and sig[0] in ("move", "grow", "erode"):
            # Layer family: the same clauses with the mover painted above every stratum (fitted by the exact test)
            found = found + [Clause(cl.selector, cl.bind, cl.body, cl.term, cl.fill, top=True) for cl in found]
        if not found:
            res.unexplained_classes += 1
        else:
            found.sort(key=lambda cl: -lib.log_P(cl.letters()))
            alternatives.append(found[:16])
            res.partials.append(ClauseSet((found[0],)))

    # 3. assemble clause sets (one clause per explainable class) and verify exactly on every record.
    #    The exact test is the arbiter: an unexplained class may be a decomposition artifact (occlusion,
    #    reveal) that identity rendering reproduces, so combinations are tried regardless.
    res.creation_failed = not creation_ok
    if True:   # the exact test is the arbiter: even with unexplained classes or creations the identity combination may be exact in the quotient
        combos = itertools.product(*alternatives) if alternatives else [()]
        for n, combo in enumerate(combos):
            if n >= MAX_COMBOS or time.time() - t0 > budget_s:
                break
            cs = ClauseSet(tuple(combo), tuple(creations))
            res.combos_tried += 1
            if all(exact(cs, r.A0, r.data, r.A1, r.flag1, r.memory, res.clock) for r in recs):
                res.hypotheses.append(Hypothesis(cs, clause_set_log_weight(cs, lib), len(recs)))
                res.min_mismatch_cells = 0
            elif res.min_mismatch_cells != 0:
                mm = 0
                for r in recs:
                    fr = predicted_frame(cs, r.A0, r.data, r.memory)
                    if fr is None:
                        mm += 4096
                    else:
                        d = fr != r.A1.frame
                        for c in res.clock:
                            d &= (fr != c) & (r.A1.frame != c)
                        mm += int(np.count_nonzero(d))
                if res.min_mismatch_cells < 0 or mm < res.min_mismatch_cells:
                    res.min_mismatch_cells = mm
        res.hypotheses.sort(key=lambda h: -h.logw)
        del res.hypotheses[k_max:]
    res.exhausted = time.time() - t0 <= budget_s
