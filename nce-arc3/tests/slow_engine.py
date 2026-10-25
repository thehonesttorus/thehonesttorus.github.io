"""A deliberately slow engine for proxy tests: act() sleeps longer than the proxy timeout."""
import time


class SlowEngine:
    def __init__(self, game_id, win_levels, available_actions):
        self.last_reason = "slow"

    def act(self, obs):
        time.sleep(2.0)
        return 2, None

    def observe(self, obs, aid, data):
        pass

    def resync(self, obs):
        pass
