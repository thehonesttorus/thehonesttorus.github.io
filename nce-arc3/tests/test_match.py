import numpy as np

from nce import algebra, match


def _base():
    F = np.full((64, 64), 5, dtype=np.int8)
    F[10:13, 10:15] = 2          # rectangle A
    F[30:33, 40:45] = 2          # rectangle B (same type as A)
    F[50, 5] = 9                 # dot
    F[20:23, 20:23] = 7          # block
    return F


def test_translation_is_matched_by_parry_cost_not_proximity_heuristic():
    F0 = _base()
    F1 = F0.copy()
    F1[10:13, 10:15] = 5; F1[11:14, 12:17] = 2      # rectangle A moves by (1, 2)
    A0, A1 = algebra.build(F0), algebra.build(F1)
    tr = match.decompose(A0, A1)
    moved = tr.moved
    assert len(moved) == 1
    i, j, g = moved[0]
    assert (g.dy, g.dx) == (1, 2) and not g.rotated
    assert A0.comps[i].color == 2 and A0.comps[i].bbox[0] == 10
    assert not tr.created and not tr.annihilated and not tr.recolours
    assert [A0.comps[i].color for i, _ in tr.reshaped] == [5]       # the top stratum's complement shape changed
    assert tr.charge.sum() == 0


def test_rotation_recolour_creation_annihilation():
    F0 = _base()
    F1 = F0.copy()
    # dot annihilated
    F1[50, 5] = 5
    # block recoloured in place
    F1[20:23, 20:23] = 12
    # new component created
    F1[60:62, 60:63] = 14
    # rectangle B rotated 90 degrees (3x5 -> 5x3) and carried to a disjoint position
    F1[30:33, 40:45] = 5; F1[40:45, 50:53] = 2
    A0, A1 = algebra.build(F0), algebra.build(F1)
    tr = match.decompose(A0, A1)
    assert any(A0.comps[i].color == 7 and A1.comps[j].color == 12 for i, j in tr.recolours)
    assert [A0.comps[i].color for i in tr.annihilated] == [9] and not tr.occluded
    assert [A1.comps[j].color for j in tr.created] == [14]
    rot = [g for i, j, g in tr.moved if g.rotated]
    assert len(rot) == 1 and rot[0].l1 > 0
    assert tr.charge[9] == -1 and tr.charge[14] == 1 and tr.charge[7] == -1 and tr.charge[12] == 1


def test_in_place_rotation_of_overlapping_cells_is_a_reshape_not_a_rot_letter():
    # A D4 change whose cells overlap the original is identity + occlusion/reshape (cost 0),
    # never a Rot letter (cost log 7): the cheaper word wins under the Parry kernel.
    F0 = _base()
    F1 = F0.copy()
    F1[30:33, 40:45] = 5; F1[30:35, 40:43] = 2
    A0, A1 = algebra.build(F0), algebra.build(F1)
    tr = match.decompose(A0, A1)
    assert not [g for _, _, g in tr.moved if g.rotated]
    assert any(A0.comps[i].color == 2 and A0.comps[i].bbox[0] == 30 for i, _ in tr.reshaped)
    assert not tr.created and not tr.annihilated


def test_null_transition():
    F0 = _base()
    A0, A1 = algebra.build(F0), algebra.build(F0.copy())
    tr = match.decompose(A0, A1)
    assert tr.is_null and not tr.moved and not tr.created and not tr.annihilated
    assert all(g.is_identity for _, _, g in tr.moves) and len(tr.moves) == len(A0.comps)


def test_greedy_fallback_agrees_on_simple_case():
    F0 = _base(); F1 = F0.copy(); F1[10:13, 10:15] = 5; F1[12:15, 10:15] = 2
    A0, A1 = algebra.build(F0), algebra.build(F1)
    saved = match._lsa
    try:
        match._lsa = None
        tr = match.decompose(A0, A1)
    finally:
        match._lsa = saved
    assert len(tr.moved) == 1 and (tr.moved[0][2].dy, tr.moved[0][2].dx) == (2, 0)


def test_shared_translation_is_matched_as_one_gauge_element():
    """Three identical blocks all shift right by 3: one clause explains all three, so the shared gauge element wins
    over the per-pair reading (one block jumping 9, two staying)."""
    F0 = np.zeros((64, 64), dtype=np.int8)
    F1 = np.zeros((64, 64), dtype=np.int8)
    for x in (10, 13, 16):
        F0[20:22, x:x + 2] = 3
        F1[20:22, x + 3:x + 5] = 3
    A0, A1 = algebra.build(F0), algebra.build(F1)
    tr = match.decompose(A0, A1)
    shifts = sorted((g.dy, g.dx) for _, _, g in tr.moves if A0.comps[_].color == 3) if False else \
        sorted((g.dy, g.dx) for i, j, g in tr.moves if A0.comps[i].color == 3)
    assert shifts == [(0, 3)] * 3, shifts


def _blocks(positions, colour=3):
    F = np.zeros((64, 64), dtype=np.int8)
    for (y, x) in positions:
        F[y:y + 2, x:x + 2] = colour
    return F


def test_history_resolves_a_jump_over_an_identical_block():
    """A mover jumps 8 over a static twin at the midpoint. Alone, 'both moved 4' is the shorter description; with the
    action's history (the mover has jumped 8 before) the identity plus the known jump costs nothing."""
    A0, A1 = algebra.build(_blocks([(20, 10), (20, 14)])), algebra.build(_blocks([(20, 14), (20, 18)]))
    key = A0.comps[[c.color for c in A0.comps].index(3)].type
    tr = match.decompose(A0, A1, prefer=frozenset({(key, 0, 8)}))
    shifts = sorted((g.dy, g.dx) for i, j, g in tr.moves if A0.comps[i].color == 3)
    assert shifts == [(0, 0), (0, 8)], shifts


def test_history_resolves_two_long_moves_that_would_cross():
    A0 = algebra.build(_blocks([(5, 40), (2, 50)]))
    A1 = algebra.build(_blocks([(5, 48), (11, 50)]))
    key = A0.comps[[c.color for c in A0.comps].index(3)].type
    tr = match.decompose(A0, A1, prefer=frozenset({(key, 0, 8), (key, 9, 0)}))
    shifts = sorted((g.dy, g.dx) for i, j, g in tr.moves if A0.comps[i].color == 3)
    assert shifts == [(0, 8), (9, 0)], shifts
