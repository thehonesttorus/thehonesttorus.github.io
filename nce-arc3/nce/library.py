"""The initial graded library: letters as families with statistical dimensions, and its critical state.

Spec: final-architecture.md Sections 3.1-3.5. Vertices V = {F, C, P, TOP}. Every letter is a
family of correspondences; its dimension is the number of summands (Theorem 1.5(iii),
Theorem 2.1(i)), fixed by interface facts (palette 16, lattice 4096, |D4| = 8, win_levels)
or by a snapshot |Xi_0| = max(2, |Xi|) of the shape alphabet taken at level start and at
every library change (T11). Units (specific gauge elements) are fitted parameters, never
edges (Section 3.2); the d = 1 subgraph {SelPointed, SelAll: F->C; Over: F->TOP} is acyclic,
so beta_c exists (Theorem 2.1(iv)) and this is asserted at every refresh. (SelPointed now has
two members, component and cell, so the d = 1 subgraph is {SelAll: F->C; Over: F->TOP}.)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from . import kms

F, C, P, TOP = "F", "C", "P", "TOP"
GRID_CELLS = 64 * 64
PALETTE = 16
D4_ORDER = 8


@dataclass
class LetterSpec:
    name: str
    src: str
    dst: str
    dim: float


def letter_specs(xi0: int, win_levels: int, promoted: Optional[list[LetterSpec]] = None) -> list[LetterSpec]:
    xi0 = max(2, int(xi0))
    specs = [
        LetterSpec("SelPointed", F, C, 2),          # the pointed component, or the pointed cell (finest projection)
        LetterSpec("SelAll", F, C, 1),
        LetterSpec("SelColour", F, C, PALETTE),
        LetterSpec("SelRank", F, C, PALETTE),
        LetterSpec("SelType", F, C, xi0),
        LetterSpec("Step", C, C, 4),
        LetterSpec("Move", C, C, 127 * 127 - 1),   # translations v in [-63,63]^2 \ {0}: 16128 summands
        LetterSpec("Rot", C, C, D4_ORDER - 1),
        LetterSpec("Recolour", C, C, PALETTE - 1),
        LetterSpec("Fill", C, C, PALETTE + 1),      # palette colours + the remembered canvas (reveal)
        LetterSpec("Layer", C, C, 2),              # mover painted in HN order, or above every stratum
        LetterSpec("Grow", C, C, 4 * PALETTE),      # P -> P v (T_d P ^ Q_c): advance one cell along d into colour c
        LetterSpec("Erode", C, C, 4),              # P -> P ^ T_d P: lose the trailing edge along d
        LetterSpec("Copy", C, C, 2),
        LetterSpec("Annihilate", C, F, PALETTE),
        LetterSpec("Create", F, C, xi0 * D4_ORDER * GRID_CELLS),
        LetterSpec("Bind", C, P, 6 * (xi0 + 1) + 2 * (PALETTE + xi0 + 8)),   # pair predicates; unary colour/type/rank tests on the component or on the pointed one
        LetterSpec("Meter", C, P, 2 * PALETTE * (GRID_CELLS + 1)),   # a colour mass compared with a threshold, either side
        LetterSpec("Gate", P, C, 2),
        LetterSpec("Const", F, F, 2),
        LetterSpec("LevelGate", F, F, max(2, int(win_levels))),
        LetterSpec("Over", F, TOP, 1),
    ]
    if promoted:
        specs.extend(promoted)
    return specs


@dataclass
class Library:
    """The typed library with its critical KMS state. Immutable except through `refresh`."""
    win_levels: int
    xi0: int = 2
    promoted: list[LetterSpec] = field(default_factory=list)
    composites: dict[str, tuple[str, ...]] = field(default_factory=dict)   # promoted name -> the letter path it abbreviates
    gauge_subgroup: tuple[int, ...] = tuple(range(8))      # H <= D4 (Section 8.4)
    beta_c: float = 0.0
    log_v: dict[str, float] = field(default_factory=dict)
    u: dict[str, float] = field(default_factory=dict)
    dims: dict[str, float] = field(default_factory=dict)
    _graph: kms.Library = field(default_factory=kms.Library, repr=False)

    def __post_init__(self) -> None:
        self.refresh()

    # ---- construction / refresh (T11: recompute only when the library changes) ----------
    def refresh(self) -> None:
        g = kms.Library()
        for s in letter_specs(self.xi0, self.win_levels, self.promoted):
            g.add(kms.Primitive(s.name, s.src, s.dst, float(s.dim)))
        self._graph = g
        self.dims = {s.name: float(s.dim) for s in letter_specs(self.xi0, self.win_levels, self.promoted)}
        self.assert_beta_c_exists()
        st = kms.CriticalState.of(g)
        if st is None:  # pragma: no cover (guarded by the assertion)
            raise RuntimeError("no critical state")
        self.beta_c = st.beta_c
        self.log_v = {v: float(np.log(st.v[st.index[v]])) for v in g.vertices}
        self.u = {v: float(st.u[st.index[v]]) for v in g.vertices}

    def assert_beta_c_exists(self) -> None:
        ones = [p for p in self._graph.primitives.values() if p.dim <= 1.0 + 1e-12]
        sub = kms.Library({p.name: p for p in ones}, list(self._graph.vertices))
        if ones and kms.spectral_radius(sub.matrix(0.0)) >= 1.0 - 1e-12:
            raise AssertionError("d = 1 subgraph has a cycle: beta_c would not exist")
        if kms.spectral_radius(self._graph.matrix(0.0)) <= 1.0:
            raise AssertionError("r(A_0) <= 1: library too sparse for a critical state")

    def snapshot_alphabet(self, n_types: int) -> None:
        new = max(2, int(n_types))
        if new != self.xi0:
            self.xi0 = new
            self.refresh()

    def with_promoted(self, spec: LetterSpec, composite: tuple[str, ...]) -> "Library":
        comps = dict(self.composites)
        comps[spec.name] = tuple(composite)
        return Library(self.win_levels, self.xi0, list(self.promoted) + [spec], comps, self.gauge_subgroup)

    def abbreviate(self, letters: list[str]) -> list[str]:
        """Rewrite a letter path with the promoted composites (longest first): the shorter coding (T5, T10)."""
        out = tuple(letters)
        for name, comp in sorted(self.composites.items(), key=lambda kv: -len(kv[1])):
            out = kms.substitute(out, comp, name)
        return list(out)

    # ---- the conditional Parry kernel (E3 reading) -------------------------------------
    def log_P(self, letters: list[str]) -> float:
        """log P(w | start) = log v_r - log v_s - beta_c * sum log d_e for a typed path of letter names."""
        if not letters:
            return 0.0
        if self.composites:
            letters = self.abbreviate(letters)
        prims = [self._graph.primitives[n] for n in letters]
        for a, b in zip(prims, prims[1:]):
            if a.target != b.source:
                raise ValueError(f"not a typed path: {a.name}({a.target}) -> {b.name}({b.source})")
        return self.log_v[prims[-1].target] - self.log_v[prims[0].source] - self.beta_c * float(sum(np.log(p.dim) for p in prims))

    def log_weight_cost(self, letters: list[str]) -> float:
        """-beta_c * sum log d_e, the path's energy term alone (used for ranking clauses of equal endpoints)."""
        if self.composites:
            letters = self.abbreviate(letters)
        return -self.beta_c * float(sum(np.log(self.dims[n]) for n in letters))

    @property
    def graph(self) -> kms.Library:
        return self._graph
