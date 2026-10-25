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
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .algebra import BaseAlgebra, Component
from .library import Library
from .store import Record, Store
from .words import Clause, ClauseSet, Creation, Hypothesis, clause_set_log_weight, d4_relative, exact, predicted_frame, select, _pair_predicate

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


# ------------------------------------------------------------------ behaviour classes

def _signature(rec: Record, cid: int) -> Optional[tuple]:
    """Outcome of component cid in record rec: None for identity/reshape, else a hashable signature."""
    tr = rec.tr
    for i, j, g in tr.moves:
        if i == cid:
            if g.is_identity:
                return None
            h = d4_relative(g.g_before, g.g_after) if g.rotated else 0
            return ("move", g.dy, g.dx, h)
    for i, j in tr.recolours:
        if i == cid:
            return ("recolour", rec.A1.comps[j].color)
    if cid in tr.annihilated:
        return ("annihilate",)
    return None


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
    if sig[0] == "recolour":
        return (("recolour", sig[1]),)
    return ()


def _fill_colour(rec: Record, cids: list[int]) -> Optional[int]:
    """Colour of the cells vacated by the components cids in rec (must be a single colour), None if inconsistent."""
    A0, A1, tr = rec.A0, rec.A1, rec.tr
    new_cover = np.zeros(64 * 64, dtype=bool)
    for i, j, g in tr.moves:
        if not g.is_identity:
            new_cover[A1.comps[j].cells] = True
    for i, j in tr.recolours:
        new_cover[A1.comps[j].cells] = True
    for j in tr.created:
        new_cover[A1.comps[j].cells] = True
    colours = set()
    for cid in cids:
        cells = A0.comps[cid].cells
        vac = cells[~new_cover[cells]]
        if vac.size:
            colours.update(np.unique(A1.frame.ravel()[vac]).tolist())
    if len(colours) > 1:
        return None
    return colours.pop() if colours else -1   # -1: nothing vacated (fill unconstrained)


# ------------------------------------------------------------------ selector fitting

def _selector_candidates(recs: list[Record]) -> list[tuple]:
    cands: list[tuple] = [("all",), ("pointed",)]
    colours, ranks, keys = set(), set(), set()
    for r in recs:
        for c in r.A0.comps:
            colours.add(c.color); ranks.add(c.stratum); keys.add(c.key)
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
    return preds


def _fits_exactly(selector: tuple, recs: list[Record], targets: list[set[int]]) -> bool:
    for r, B in zip(recs, targets):
        sel = {c.id for c in select(selector, r.A0, r.data)}
        if sel != B:
            return False
    return True


def _fits_with_gate(selector: tuple, pred: tuple, recs: list[Record], targets: list[set[int]]) -> bool:
    for r, B in zip(recs, targets):
        sel = {c.id for c in select(selector, r.A0, r.data)}
        if not B <= sel:
            return False
        for cid in sel:
            if _pair_predicate(pred, r.A0.comps[cid], r.A0) != (cid in B):
                return False
    return True


# ------------------------------------------------------------------ main entry

def synthesize(a: int, store: Store, lib: Library, budget_s: float = 0.25, k_max: int = K_MAX) -> SynthResult:
    t0 = time.time()
    res = SynthResult()
    recs_all = store.transitions(a)
    recs = [r for r in recs_all if r.flag1 != "GAME_OVER"]
    terminal = [r for r in recs_all if r.flag1 == "GAME_OVER"]
    if not recs:
        if terminal:
            cs = ClauseSet(over=True)
            res.hypotheses.append(Hypothesis(cs, clause_set_log_weight(cs, lib), len(terminal)))
        res.exhausted = True
        return res

    # 1. behaviour classes across all non-terminal records
    classes: dict[tuple, list[set[int]]] = {}
    for ri, r in enumerate(recs):
        for c in r.A0.comps:
            sig = _signature(r, c.id)
            if sig is None:
                continue
            classes.setdefault(sig, [set() for _ in recs])[ri].add(c.id)
    res.classes = len(classes)
    res.reshaped_records = sum(1 for r in recs if any(r.A0.comps[i].stratum > 0 for i, _ in r.tr.reshaped))

    # creations: constant (Create[tau,g,t]) or relative to the pointed component (SelPointed . Copy . units)
    creations: list[Creation] = []
    created_sets = [tuple(sorted((r.A1.comps[j].key, r.A1.comps[j].color, r.A1.comps[j].g, r.A1.comps[j].offset) for j in r.tr.created)) for r in recs]
    creation_ok = all(cs == created_sets[0] for cs in created_sets)
    if creation_ok:
        creations = [Creation(k, col, g, off) for (k, col, g, off) in created_sets[0]]
    elif all(r.data for r in recs):
        rel_sets = []
        for r in recs:
            cid = int(r.A0.cell2comp[int(r.data["y"]), int(r.data["x"])])
            if cid < 0:
                rel_sets = None; break
            base = r.A0.comps[cid].offset
            rel_sets.append(tuple(sorted((r.A1.comps[j].key, r.A1.comps[j].color, r.A1.comps[j].g, (r.A1.comps[j].offset[0] - base[0], r.A1.comps[j].offset[1] - base[1])) for j in r.tr.created)))
        if rel_sets and all(cs == rel_sets[0] for cs in rel_sets):
            creation_ok = True
            creations = [Creation(k, col, g, off, relative=True) for (k, col, g, off) in rel_sets[0]]

    selectors = _selector_candidates(recs)
    predicates = _predicate_candidates(recs)

    # 2. per class: clauses whose selector/gate separates exactly, with constant fill
    alternatives: list[list[Clause]] = []
    for sig, targets in classes.items():
        body = _body_for(sig)
        term = ("annihilate",) if sig[0] == "annihilate" else ("stay",)
        vacates = sig[0] != "recolour"
        fill: Optional[int] = None
        if vacates:
            fills = {_fill_colour(r, sorted(B)) for r, B in zip(recs, targets) if B}
            fills.discard(-1)
            if len(fills) > 1 or None in fills:
                res.unexplained_classes += 1
                res.fill_inconsistent += 1
                continue
            fill = fills.pop() if fills else int(recs[0].A0.strata[0][0])   # nothing vacated: fill never rendered
        found: list[Clause] = []
        for selector in selectors:
            if time.time() - t0 > budget_s:
                break
            if _fits_exactly(selector, recs, targets):
                found.append(Clause(selector, None, body, term, fill))
        if not found:
            for selector in selectors:
                for pred in predicates:
                    if time.time() - t0 > budget_s:
                        break
                    if _fits_with_gate(selector, pred, recs, targets):
                        found.append(Clause(selector, pred, body, term, fill))
        if not found:
            res.unexplained_classes += 1
        else:
            found.sort(key=lambda cl: -lib.log_P(cl.letters()))
            alternatives.append(found[:8])
            res.partials.append(ClauseSet((found[0],)))

    # 3. assemble clause sets (one clause per class) and verify exactly on every record
    res.creation_failed = not creation_ok
    if res.unexplained_classes == 0 and creation_ok:
        combos = itertools.product(*alternatives) if alternatives else [()]
        for n, combo in enumerate(combos):
            if n >= MAX_COMBOS or time.time() - t0 > budget_s:
                break
            cs = ClauseSet(tuple(combo), tuple(creations))
            res.combos_tried += 1
            if all(exact(cs, r.A0, r.data, r.A1, r.flag1) for r in recs):
                res.hypotheses.append(Hypothesis(cs, clause_set_log_weight(cs, lib), len(recs)))
                res.min_mismatch_cells = 0
            elif res.min_mismatch_cells != 0:
                mm = 0
                for r in recs:
                    fr = predicted_frame(cs, r.A0, r.data)
                    mm += 4096 if fr is None else int(np.count_nonzero(fr != r.A1.frame))
                if res.min_mismatch_cells < 0 or mm < res.min_mismatch_cells:
                    res.min_mismatch_cells = mm
        res.hypotheses.sort(key=lambda h: -h.logw)
        del res.hypotheses[k_max:]
    res.exhausted = time.time() - t0 <= budget_s
    return res
