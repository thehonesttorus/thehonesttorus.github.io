import numpy as np
import pytest

from nce import algebra


def _frame_with_shapes():
    F = np.full((64, 64), 5, dtype=np.int8)            # dominant colour (top HN stratum)
    F[10:13, 10:15] = 2                                # 3x5 rectangle
    F[30:32, 40:41] = 2                                # 2x1 bar, same colour, separate component
    F[50, 5] = 9                                       # single cell
    F[20:23, 20:23] = 7                                # 3x3 block
    F[21, 21] = 11                                     # hole colour inside the block
    L = np.zeros((4, 4), dtype=bool); L[:, 0] = True; L[3, :] = True   # an L shape
    F[40:44, 50:54][L] = 3
    return F


def test_label4_numpy_matches_scipy_counts():
    F = _frame_with_shapes()
    for c in np.unique(F):
        lab_s, n_s = algebra.label4(F == c)
        lab_n, n_n = algebra._label4_numpy(F == c)
        assert n_s == n_n
        # same partition up to relabelling
        assert len(set(zip(lab_s.ravel().tolist(), lab_n.ravel().tolist()))) == n_s + (1 if (F != c).any() else 0)


def test_build_components_strata_contact():
    F = _frame_with_shapes()
    A = algebra.build(F)
    assert A.strata[0] == [5]                                   # maximal slope = dominant colour
    types = {(c.color, c.n) for c in A.comps}
    assert (2, 15) in types and (2, 2) in types and (9, 1) in types and (7, 8) in types and (11, 1) in types and (3, 7) in types
    assert all(A.cell2comp[np.divmod(c.cells, 64)].tolist() == [c.id] * c.n for c in A.comps)
    # the hole (colour 11) is in contact with the block (colour 7) and with nothing else
    hole = next(c for c in A.comps if c.color == 11)
    block = next(c for c in A.comps if c.color == 7)
    nb = A.neighbours(hole.id)
    assert nb == [block.id]
    # the two colour-2 components have different shape keys (3x5 vs 2x1)
    k = [c.key for c in A.comps if c.color == 2]
    assert len(set(k)) == 2
    assert next(c for c in A.comps if c.n == 1).hull_rank == 1
    bar = next(c for c in A.comps if c.n == 2); assert bar.hull_rank == 2
    assert block.hull_rank == 3


def test_canonical_key_is_d4_and_translation_invariant():
    L = np.zeros((4, 3), dtype=bool); L[:, 0] = True; L[3, :] = True
    keys = set()
    for g in range(8):
        img = algebra.d4_apply(L, g)
        key, gg, canon = algebra.canonical_form(img)
        keys.add(key)
        assert np.array_equal(algebra.d4_apply(img, gg), canon)
    assert len(keys) == 1
    # translation invariance through build()
    F1 = np.zeros((64, 64), dtype=np.int8); F1[5:9, 5:8][L] = 4
    F2 = np.zeros((64, 64), dtype=np.int8); F2[40:44, 20:23][L] = 4
    a1 = [c for c in algebra.build(F1).comps if c.color == 4][0]
    a2 = [c for c in algebra.build(F2).comps if c.color == 4][0]
    assert a1.key == a2.key and a1.g == a2.g and a1.offset == (5, 5) and a2.offset == (40, 20)


def test_d4_inverse_and_transport():
    L = np.zeros((4, 3), dtype=bool); L[:, 0] = True; L[3, :] = True
    for g in range(8):
        h = algebra.d4_inverse(g)
        assert np.array_equal(algebra.d4_apply(algebra.d4_apply(L, g), h), L)
    F = np.zeros((64, 64), dtype=np.int8); F[5:9, 5:8][L] = 4
    comp = [c for c in algebra.build(F).comps if c.color == 4][0]
    # pure translation keeps the same gauge element
    cells = algebra.transport_cells(comp, comp.g, (10, 20))
    F2 = algebra.render([(4, cells)], np.zeros((64, 64), dtype=np.int8))
    comp2 = [c for c in algebra.build(F2).comps if c.color == 4][0]
    assert comp2.key == comp.key and comp2.g == comp.g and comp2.offset == (10, 20)
    # rotation by a D4 element changes g but not the key
    g2 = (comp.g + 1) % 4
    cells3 = algebra.transport_cells(comp, g2, (0, 0))
    comp3 = [c for c in algebra.build(algebra.render([(4, cells3)], np.zeros((64, 64), dtype=np.int8))).comps if c.color == 4][0]
    assert comp3.key == comp.key and comp3.g == g2
    assert algebra.transport_cells(comp, comp.g, (62, 62)) is None


def test_build_speed():
    import time
    rng = np.random.default_rng(0)
    F = rng.integers(0, 16, size=(64, 64), dtype=np.int8)      # worst case: ~thousands of tiny components
    t = time.time(); A = algebra.build(F); dt = time.time() - t
    assert sum(c.n for c in A.comps) == 64 * 64                  # uncapped: every cell belongs to a component
    assert dt < 4.0
    F2 = _frame_with_shapes()
    t = time.time()
    for _ in range(20): algebra.build(F2)
    assert (time.time() - t) / 20 < 0.02                        # typical frame: well under 20 ms
