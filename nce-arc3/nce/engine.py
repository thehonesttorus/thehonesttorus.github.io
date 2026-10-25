"""TEMPORARY stub engine (random actions) used only to validate plumbing. Replaced by the real engine."""
import random

class Engine:
    def __init__(self, game_id, win_levels, available_actions):
        self.rng = random.Random(0); self.last_reason = "stub"
    def act(self, obs):
        avail = [a for a in obs["available_actions"] if a != 0] or [1]
        a = self.rng.choice(avail)
        return a, ({"x": self.rng.randint(0, 63), "y": self.rng.randint(0, 63)} if a == 6 else None)
    def observe(self, obs, aid, data): pass
    def on_game_over(self, obs): pass
    def want_stop(self, obs): return False
