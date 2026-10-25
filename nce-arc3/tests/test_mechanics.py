"""End-to-end toy mechanics: the agent must form exact words for tiny synthetic games.

The toys here are not ARC games and encode no game knowledge: each is one correspondence
family from the library exercised in isolation (a gauge move with a lethal wall, a cell-level
recolour under the pointed selector, an in-place D4 change, boundary growth, a level-gated
rule) so the synthesis path for that family is verified.
"""
import numpy as np
import pytest

from nce import algebra, kms, maternity
from nce.agent import NCEAgent
from nce.library import Library
from nce.store import Store
from nce import synth
from nce.words import predicted


def obs_of(frame, state="NOT_FINISHED", avail=(1, 2, 3, 4), full_reset=False, levels=0):
    return dict(frames=[frame.copy()], state=state, levels_completed=levels, win_levels=1,
                available_actions=list(avail), full_reset=full_reset)


def drive(agent, env, n):
    """The harness loop: act, step, observe; RESET after a death."""
    o = env.reset_obs()
    for _ in range(n):
        if o["state"] == "GAME_OVER":
            agent.on_game_over(o)
            aid, data = 0, None
        else:
            aid, data = agent.act(o)
        o = env.step(aid, data)
        agent.observe(o, aid, data)
    return o


class Mover:
    """A 2x2 block moved one cell by actions 1-4 (up/down/left/right); entering the wall kills."""
    DIRS = {1: (-1, 0), 2: (1, 0), 3: (0, -1), 4: (0, 1)}

    def __init__(self):
        self.reset_obs()

    def frame(self):
        F = np.zeros((64, 64), dtype=np.int8)
        F[0, :] = F[-1, :] = F[:, 0] = F[:, -1] = 1
        F[self.y:self.y + 2, self.x:self.x + 2] = 3
        return F

    def reset_obs(self):
        self.y, self.x = 2, 30
        self.dead = False
        return obs_of(self.frame(), full_reset=True)

    def step(self, aid, data):
        if aid == 0:
            return self.reset_obs()
        dy, dx = self.DIRS[aid]
        y, x = self.y + dy, self.x + dx
        if y < 1 or x < 1 or y + 1 > 62 or x + 1 > 62:
            return obs_of(self.frame(), state="GAME_OVER")
        self.y, self.x = y, x
        return obs_of(self.frame())


def test_mover_learns_gauge_words_and_a_death_gate():
    env = Mover()
    agent = NCEAgent("toy-mover", 1, [1, 2, 3, 4])
    drive(agent, env, 40)
    d = agent.diagnostics()
    assert all(d["H_sizes"].get(a, 0) > 0 for a in (2, 3, 4)), d
    # action 1 (up) kills from the start row: its coherent set must carry a death gate (the toy starts at row 2)
    assert d["deaths"] >= 1
    hs = agent.H.get(1) or []
    assert hs and any(h.cs.over_gate is not None or h.cs.over for h in hs), [h.cs.describe() for h in hs]
    # the top word for "down" predicts the next frame exactly
    A = agent.perceive(env.frame())
    fr, flag = predicted(agent.H[2][0].cs, A, None, agent.memory)
    env2 = Mover(); env2.y, env2.x = env.y, env.x
    assert flag == "NOT_FINISHED" and np.array_equal(fr, env2.step(2, None)["frames"][0])


class CellToggle:
    """Clicking a cell of the big region toggles it between colours 4 and 7."""

    def __init__(self):
        self.F = np.zeros((64, 64), dtype=np.int8)
        self.F[10:40, 10:40] = 4
        self.F[12, 12] = 7
        self.F[20:22, 30:32] = 7

    def reset_obs(self):
        return obs_of(self.F, avail=(6,), full_reset=True)

    def step(self, aid, data):
        if aid == 6 and data is not None:
            x, y = int(data["x"]), int(data["y"])
            if self.F[y, x] == 4:
                self.F[y, x] = 7
            elif self.F[y, x] == 7:
                self.F[y, x] = 4
        return obs_of(self.F, avail=(6,))


def test_cell_toggle_is_a_pointed_cell_word():
    env = CellToggle()
    agent = NCEAgent("toy-toggle", 1, [6])
    # scripted clicks on distinct cells of both colours, then the agent's own clicks
    rng = np.random.default_rng(0)
    o = env.reset_obs()
    agent.act(o)
    for k in range(10):
        y, x = int(rng.integers(10, 40)), int(rng.integers(10, 40))
        if k % 4 == 3:
            y, x = 12, 12
        o = env.step(6, {"x": x, "y": y})
        agent.observe(o, 6, {"x": x, "y": y})
    drive(agent, env, 6)
    hs = agent.H.get(6) or []
    assert hs, agent.diagnostics()
    A = agent.perceive(env.F)
    for (y, x) in ((15, 15), (12, 12), (20, 30), (39, 39)):
        fr, flag = predicted(hs[0].cs, A, {"x": x, "y": y}, agent.memory)
        exp = env.F.copy(); exp[y, x] = 7 if exp[y, x] == 4 else 4
        assert fr is not None and np.array_equal(fr, exp), (y, x, hs[0].cs.describe())


class Rotator:
    """Action 1 applies a D4 quarter turn to an L-tromino about its bbox corner (an in-place change)."""

    def __init__(self):
        self.g = 0
        self.canon = np.array([[1, 0], [1, 1]], dtype=bool)

    def frame(self):
        F = np.zeros((64, 64), dtype=np.int8)
        img = algebra.d4_apply(self.canon, self.g)
        h, w = img.shape
        F[20:20 + h, 20:20 + w][img] = 4
        return F

    def reset_obs(self):
        return obs_of(self.frame(), avail=(1,), full_reset=True)

    def step(self, aid, data):
        from nce.words import d4_compose
        if aid == 1:
            self.g = d4_compose(1, self.g)
        return obs_of(self.frame(), avail=(1,))


def test_in_place_rotation_is_found_on_the_second_reading():
    env = Rotator()
    agent = NCEAgent("toy-rot", 1, [1])
    drive(agent, env, 8)
    hs = agent.H.get(1) or []
    assert hs, agent.diagnostics()
    assert any(k == "rot" for cl in hs[0].cs.clauses for k, _ in cl.body), hs[0].cs.describe()
    A = agent.perceive(env.frame())
    fr, _ = predicted(hs[0].cs, A, None, agent.memory)
    assert np.array_equal(fr, env.step(1, None)["frames"][0])


class GrowBar:
    """A 4-tall bar of colour 14 advances one column to the right into a colour-0 field on every action 1."""

    def __init__(self):
        self.w = 3

    def frame(self):
        F = np.full((64, 64), 9, dtype=np.int8)
        F[5:9, 5:45] = 0
        F[5:9, 5:5 + self.w] = 14
        F[30:40, 30:40] = 2
        return F

    def reset_obs(self):
        return obs_of(self.frame(), avail=(1, 2), full_reset=True)

    def step(self, aid, data):
        if aid == 1:
            self.w = min(self.w + 1, 40)
        return obs_of(self.frame(), avail=(1, 2))


def test_boundary_growth_is_a_grow_word():
    env = GrowBar()
    agent = NCEAgent("toy-grow", 1, [1, 2])
    drive(agent, env, 10)
    hs = agent.H.get(1) or []
    assert hs, agent.diagnostics()
    assert any(k == "grow" for cl in hs[0].cs.clauses for k, _ in cl.body), hs[0].cs.describe()
    A = agent.perceive(env.frame())
    fr, _ = predicted(hs[0].cs, A, None, agent.memory)
    assert np.array_equal(fr, env.step(1, None)["frames"][0])
    assert agent.H.get(2) and all(not cl.body for cl in agent.H[2][0].cs.clauses)   # action 2 is the identity word


def test_level_gate_when_the_rule_changes_between_levels():
    lib = Library(win_levels=3)
    store = Store()

    def frame(x):
        F = np.zeros((64, 64), dtype=np.int8); F[10:12, x:x + 2] = 5; return F
    for x in (10, 12, 14):      # level 0: action 1 moves right
        store.append(1, algebra.build(frame(x)), None, algebra.build(frame(x + 1)), "NOT_FINISHED", 0)
    for x in (30, 29, 28):      # level 1: action 1 moves left
        store.append(1, algebra.build(frame(x)), None, algebra.build(frame(x - 1)), "NOT_FINISHED", 1)
    r = synth.synthesize(1, store, lib, budget_s=2.0)
    assert r.hypotheses and r.level_gated and all(h.cs.level_gated for h in r.hypotheses)
    A = algebra.build(frame(20))
    fr, _ = predicted(r.hypotheses[0].cs, A, None)
    assert np.array_equal(fr, frame(19))


def test_promotion_is_exact_and_abbreviates_words():
    lib = Library(win_levels=2)
    word = [("SelColour", ("colour", 3)), ("Step", ("step", (0, 1))), ("Step", ("step", (0, 1))), ("Fill", 0)]
    corpus = [list(word) for _ in range(12)] + [[("SelAll", ("all",)), ("Recolour", ("recolour", 5))]] * 5
    lib2, names = maternity.promote(lib, corpus)
    assert names, "twelve repeats of one four-letter unit must pay for a composite letter of dimension 1"
    assert lib2.beta_c > lib.beta_c                       # adjoining raises the critical temperature (E6)
    plain = [n for n, _ in word]
    assert len(lib2.abbreviate(plain)) < len(plain)
    assert lib2.log_P(plain) > lib.log_P(plain)           # the promoted word is coded shorter
    other = ["SelColour", "Rot", "Annihilate"]            # a word outside the corpus pays for the higher beta_c
    assert lib2.log_weight_cost(other) < lib.log_weight_cost(other)
    comp = tuple(lib2.composites[names[0]])
    d = kms.promotion_decision(lib.graph, [plain] * 12 + [["SelAll", "Recolour"]] * 5, comp, "Kx", d_new=1.0)
    assert d.promote and d.gain > 0


class MoverWithCountdown(Mover):
    """The mover of `Mover`, plus a 1-cell-thick bar of colour 12 that loses its rightmost cell on every action."""

    def reset_obs(self):
        self.bar = 20
        return super().reset_obs()

    def frame(self):
        F = super().frame()
        F[60, 5:5 + self.bar] = 12
        return F

    def step(self, aid, data):
        if aid != 0:
            self.bar = max(self.bar - 1, 1)
        return super().step(aid, data)


def test_action_independent_countdown_is_an_erode_clause_in_every_word():
    env = MoverWithCountdown()
    agent = NCEAgent("toy-countdown", 1, [1, 2, 3, 4])
    drive(agent, env, 16)
    d = agent.diagnostics()
    assert all(d["H_sizes"].get(a, 0) > 0 for a in (2, 3, 4)), d
    for a in (2, 3, 4):
        assert any(k == "erode" for cl in agent.H[a][0].cs.clauses for k, _ in cl.body), agent.H[a][0].cs.describe()
    A = agent.perceive(env.frame())
    fr, _ = predicted(agent.H[3][0].cs, A, None, agent.memory)
    assert np.array_equal(fr, env.step(3, None)["frames"][0])
