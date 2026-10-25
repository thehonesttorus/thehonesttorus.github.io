"""The maternity ward: library growth by exact promotion (T5, T10; amendment E6).

At a level end the corpus of explanations is the set of verified words (every member of every
coherent set), each letter carrying its fitted unit. Candidate composites are the typed
sub-paths of length 2..4 that occur in at least two corpus words. A composite's dimension as
a new letter is the number of *distinct* unit tuples with which it occurs in the corpus (the
Wold-newborn compression of Theorem 2.1: the realised summands, never the full product of the
parts' dimensions). It is adjoined iff the corpus's total critical-state log-weight increases
(kms.promotion_decision recomputes the critical state exactly before and after; adjoining
raises beta_c, which every other word pays for). Words are abbreviated by promoted letters
when weighed (Library.abbreviate), so the Parry order among exact hypotheses, the Kato purity
and the goal weights all move. Verification is untouched: words stay exact or not.
"""
from __future__ import annotations

from collections import Counter, defaultdict

from . import kms
from .library import LetterSpec, Library
from .words import Hypothesis

MAX_PROMOTIONS_PER_CALL = 3
MIN_LEN, MAX_LEN = 2, 4

Lettered = list[tuple]          # [(family name, unit), ...]


def corpus_of(H: dict[int, list[Hypothesis]]) -> list[Lettered]:
    return [path for hs in H.values() for h in hs for path in h.cs.lettered() if path]


def _abbreviate(lib: Library, word: Lettered) -> Lettered:
    """Abbreviate a lettered word by the library's composites (units of a composite: the tuple of its parts')."""
    out = list(word)
    for name, comp in sorted(lib.composites.items(), key=lambda kv: -len(kv[1])):
        k = len(comp)
        res: Lettered = []
        i = 0
        while i < len(out):
            if k and tuple(n for n, _ in out[i:i + k]) == tuple(comp):
                res.append((name, tuple(u for _, u in out[i:i + k])))
                i += k
            else:
                res.append(out[i])
                i += 1
        out = res
    return out


def candidates(words: list[Lettered]) -> list[tuple[tuple[str, ...], int]]:
    """Composites (letter-name sub-paths) occurring in >= 2 words, with their realised summand counts."""
    counts: Counter = Counter()
    units: dict[tuple[str, ...], set] = defaultdict(set)
    for w in words:
        names = tuple(n for n, _ in w)
        seen: set[tuple[str, ...]] = set()
        for k in range(MIN_LEN, MAX_LEN + 1):
            for i in range(len(names) - k + 1):
                comp = names[i:i + k]
                seen.add(comp)
                units[comp].add(tuple(u for _, u in w[i:i + k]))
        counts.update(seen)
    out = [(c, len(units[c])) for c, n in counts.items() if n >= 2]
    out.sort(key=lambda cn: (-counts[cn[0]], cn[1], len(cn[0])))
    return out


def promote(lib: Library, corpus: list[Lettered], kappa: float = 0.0) -> tuple[Library, list[str]]:
    """Promote up to MAX_PROMOTIONS_PER_CALL composites, each decided exactly against the current corpus."""
    names: list[str] = []
    if not corpus:
        return lib, names
    words = [_abbreviate(lib, w) for w in corpus]
    for comp, d_new in candidates(words):
        if len(names) >= MAX_PROMOTIONS_PER_CALL:
            break
        name = f"K{len(lib.composites) + 1}"
        plain = [[n for n, _ in w] for w in words]
        try:
            d = kms.promotion_decision(lib.graph, plain, comp, name, d_new=float(d_new), kappa=kappa)
        except (ValueError, KeyError):
            continue
        if not d.promote:
            continue
        ps = [lib.graph.primitives[n] for n in comp]
        try:
            lib = lib.with_promoted(LetterSpec(name, ps[0].source, ps[-1].target, float(d_new)), comp)
        except AssertionError:
            continue       # the composite would break the critical state's existence: not adjoined
        words = [_abbreviate(lib, w) for w in words]
        names.append(name)
    return lib, names
