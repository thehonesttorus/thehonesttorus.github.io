"""Critical KMS (Parry) state of a typed correspondence library, and library growth.

Theory: Theorem 2.1(iv) (Laca-Neshveyev / Parry), Proposition 9.3 with amendment E6,
and the derivation in docs/promotion-inequality.md.

A library is a finite directed multigraph G = (V, E, s, r). Each edge e is a typed
primitive correspondence with statistical dimension d_e >= 1 (Longo dimension; the
product of dimensions along a path is the dimension of the composite, Theorem 1.5(iii)).

    A_beta(u, v) = sum_{e: u -> v} d_e^{-beta}

beta_c is the unique root of r(A_beta) = 1 (r = spectral radius) when it exists; the
critical KMS state gives a word w = e_1...e_k from s to r the conditional weight

    omega(w | start at s) = (v_r / v_s) * prod_j d_{e_j}^{-beta_c}

where v is the Perron right eigenvector. These weights define a Markov chain on
vertices (they sum to 1 over the edges leaving s), which is the Parry measure.

Edges of dimension exactly 1 (gauge unitaries) have weight 1 at every beta; a cycle of
such edges makes r(A_beta) >= 1 for all beta and beta_c undefined. Gauge elements are
therefore not library edges: they act on vertex types (the base algebra is the crossed
product by the gauge group), see docs/promotion-inequality.md.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence

import numpy as np


@dataclass(frozen=True)
class Primitive:
    name: str
    source: str      # vertex type
    target: str      # vertex type
    dim: float       # statistical dimension d_e >= 1

    def __post_init__(self) -> None:
        if self.dim < 1.0:
            raise ValueError(f"statistical dimension must be >= 1, got {self.dim} for {self.name}")


@dataclass
class Library:
    """Typed library graph of primitive correspondences."""
    primitives: dict[str, Primitive] = field(default_factory=dict)
    vertices: list[str] = field(default_factory=list)

    def add(self, p: Primitive) -> None:
        if p.name in self.primitives:
            raise KeyError(f"primitive {p.name} already in library")
        self.primitives[p.name] = p
        for v in (p.source, p.target):
            if v not in self.vertices:
                self.vertices.append(v)

    def remove(self, name: str) -> None:
        del self.primitives[name]

    def copy(self) -> "Library":
        return Library(dict(self.primitives), list(self.vertices))

    def edges_from(self, vertex: str) -> list[Primitive]:
        return [p for p in self.primitives.values() if p.source == vertex]

    def matrix(self, beta: float) -> np.ndarray:
        idx = {v: i for i, v in enumerate(self.vertices)}
        A = np.zeros((len(self.vertices), len(self.vertices)))
        for p in self.primitives.values():
            A[idx[p.source], idx[p.target]] += p.dim ** (-beta)
        return A


def spectral_radius(A: np.ndarray) -> float:
    if A.size == 0:
        return 0.0
    return float(np.max(np.abs(np.linalg.eigvals(A))))


def perron_vectors(A: np.ndarray, iters: int = 10_000, tol: float = 1e-13) -> tuple[np.ndarray, np.ndarray, float]:
    """Right (v, sum 1) and left (u, u.v = 1) Perron vectors and the Perron root.

    Power iteration on (A + I) keeps the iteration aperiodic for periodic graphs.
    """
    n = A.shape[0]
    B = A + np.eye(n)
    v = np.ones(n) / n
    u = np.ones(n) / n
    lam = 1.0
    for _ in range(iters):
        v_new = B @ v
        lam = float(v_new.sum())
        v_new /= lam
        u_new = B.T @ u
        u_new /= u_new.sum()
        done = np.abs(v_new - v).max() < tol and np.abs(u_new - u).max() < tol
        v, u = v_new, u_new
        if done:
            break
    rho = lam - 1.0
    u = u / float(u @ v)
    return v, u, rho


def critical_beta(lib: Library, lo: float = 0.0, hi: float = 64.0, tol: float = 1e-10) -> Optional[float]:
    """Solve r(A_beta) = 1 by bisection. r is decreasing in beta (every d_e >= 1).

    Returns None when no root exists in [lo, hi]: r(A_lo) < 1 means the library is too
    sparse for a critical state (every word weight would exceed the Parry normalisation);
    r(A_hi) > 1 means a dimension-1 cycle (gauge elements left as edges) pins r >= 1.
    """
    # Dimension-1 edges keep weight 1 at every beta; if they alone carry a cycle, r(A_beta)
    # stays >= 1 for all beta and the critical point is at infinity (no KMS state).
    gauge_only = Library({k: p for k, p in lib.primitives.items() if p.dim <= 1.0 + 1e-12}, list(lib.vertices))
    if gauge_only.primitives and spectral_radius(gauge_only.matrix(0.0)) >= 1.0 - 1e-12:
        return None
    r_lo = spectral_radius(lib.matrix(lo))
    r_hi = spectral_radius(lib.matrix(hi))
    if r_lo < 1.0 - 1e-12 or r_hi > 1.0 + 1e-12:
        return None
    if abs(r_lo - 1.0) <= 1e-12:
        return lo
    a, b = lo, hi
    while b - a > tol:
        m = 0.5 * (a + b)
        if spectral_radius(lib.matrix(m)) >= 1.0:
            a = m
        else:
            b = m
    return 0.5 * (a + b)


@dataclass
class CriticalState:
    """The critical KMS state of a library: beta_c and the Perron data at beta_c."""
    beta_c: float
    v: np.ndarray            # right Perron vector, indexed like lib.vertices, sum 1
    u: np.ndarray            # left Perron vector, u.v = 1
    index: dict[str, int]

    @classmethod
    def of(cls, lib: Library) -> Optional["CriticalState"]:
        beta = critical_beta(lib)
        if beta is None:
            return None
        v, u, _ = perron_vectors(lib.matrix(beta))
        return cls(beta, v, u, {x: i for i, x in enumerate(lib.vertices)})

    def log_weight(self, lib: Library, word: Sequence[str]) -> float:
        """log omega(word | start) = log v_r - log v_s - beta_c * sum log d_e. Requires a typed path."""
        if not word:
            return 0.0
        ps = [lib.primitives[n] for n in word]
        for a, b in zip(ps, ps[1:]):
            if a.target != b.source:
                raise ValueError(f"word is not a path: {a.name} -> {b.name}")
        s, r = ps[0].source, ps[-1].target
        return float(np.log(self.v[self.index[r]]) - np.log(self.v[self.index[s]]) - self.beta_c * sum(np.log(p.dim) for p in ps))

    def weight(self, lib: Library, word: Sequence[str]) -> float:
        return float(np.exp(self.log_weight(lib, word)))

    def transition_probabilities(self, lib: Library, vertex: str) -> dict[str, float]:
        """Parry transition probabilities over the edges leaving `vertex`; they sum to 1."""
        return {p.name: self.weight(lib, [p.name]) for p in lib.edges_from(vertex)}


def rank_words(lib: Library, state: CriticalState, words: Iterable[Sequence[str]]) -> list[tuple[float, tuple[str, ...]]]:
    """Words in decreasing critical-state weight: the canonical search order (Proposition 9.3(i))."""
    scored = [(state.log_weight(lib, w), tuple(w)) for w in words]
    scored.sort(key=lambda t: -t[0])
    return scored


def corpus_log_weight(lib: Library, state: CriticalState, corpus: Iterable[Sequence[str]]) -> float:
    """S_L = sum over explanation words of log omega_L(w)  (code length of the corpus under L)."""
    return float(sum(state.log_weight(lib, w) for w in corpus))


def substitute(word: Sequence[str], composite: Sequence[str], new_name: str) -> tuple[str, ...]:
    """Replace every non-overlapping occurrence of `composite` in `word` by `new_name`."""
    word = tuple(word)
    comp = tuple(composite)
    k = len(comp)
    out: list[str] = []
    i = 0
    while i < len(word):
        if k and word[i:i + k] == comp:
            out.append(new_name)
            i += k
        else:
            out.append(word[i])
            i += 1
    return tuple(out)


@dataclass
class PromotionDecision:
    promote: bool
    gain: float               # S_{L'} - S_L
    beta_c_before: float
    beta_c_after: Optional[float]
    occurrences: int


def promotion_decision(lib: Library, corpus: Sequence[Sequence[str]], composite: Sequence[str], new_name: str,
                       d_new: Optional[float] = None, kappa: float = 0.0) -> PromotionDecision:
    """Library growth ("maternity ward") with amendment E6 applied.

    Adjoining the composite as a new primitive raises beta_c (entropy rate goes up), which
    lowers every other word's weight; the composite's occurrences in the corpus of
    explanations get the shorter, lower-dimension coding. Promote iff the corpus's total
    critical-state log-weight increases by more than kappa:

        S_{L'} - S_L > kappa,   S_L = sum_{w in corpus} log omega_L(w).

    d_new defaults to the product of the composite's dimensions (Theorem 1.5(iii)); a
    smaller value encodes a Wold-newborn compression of the composite. The before/after
    states are computed exactly (no first-order approximation).
    """
    comp = tuple(composite)
    ps = [lib.primitives[n] for n in comp]
    for a, b in zip(ps, ps[1:]):
        if a.target != b.source:
            raise ValueError("composite is not a typed path")
    prod_dim = float(np.prod([p.dim for p in ps]))
    if d_new is None:
        d_new = prod_dim
    if d_new > prod_dim + 1e-12:
        raise ValueError("d_new cannot exceed the product of the composite's dimensions")
    before = CriticalState.of(lib)
    if before is None:
        raise ValueError("library has no critical state")
    occurrences = sum(1 for w in corpus for i in range(len(w) - len(comp) + 1) if tuple(w[i:i + len(comp)]) == comp)
    lib2 = lib.copy()
    lib2.add(Primitive(new_name, ps[0].source, ps[-1].target, d_new))
    after = CriticalState.of(lib2)
    if after is None:
        return PromotionDecision(False, float("-inf"), before.beta_c, None, occurrences)
    s_before = corpus_log_weight(lib, before, corpus)
    s_after = corpus_log_weight(lib2, after, [substitute(w, comp, new_name) for w in corpus])
    gain = s_after - s_before
    return PromotionDecision(gain > kappa, gain, before.beta_c, after.beta_c, occurrences)
