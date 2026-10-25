import random
class RandomAgent:
    def __init__(self, game_id, win_levels, available_actions):
        self.rng = random.Random(0); self.avail = available_actions
    def act(self, o):
        avail = [x for x in o["available_actions"] if x != 0] or [1]
        a = self.rng.choice(avail)
        return a, ({"x": self.rng.randint(0,63), "y": self.rng.randint(0,63)} if a == 6 else None)
