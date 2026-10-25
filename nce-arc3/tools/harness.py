"""Fast local evaluation harness for ARC-AGI-3 agents (engine-level, offline).

Agent protocol (duck-typed):
    agent = AgentCls(game_id, win_levels, available_actions)
    action_id, data = agent.act(obs)      # obs: dict(frames=[np.ndarray 64x64,...], state, levels_completed, available_actions, full_reset)
    agent.observe(obs)                    # optional, called after each step with the new obs

Scores per level use the per-game baseline_actions from metadata and the RHAE formula.
Usage: python harness.py <module:Class> [--games ls20,vc33] [--max-actions 2000] [--procs 4] [--time-per-game 600]
"""
import argparse, importlib, json, logging, os, sys, time, traceback
from multiprocessing import Pool
import numpy as np

def rhae(level_actions, baselines, n_levels):
    # level_actions: list of actions used for each completed level (sequential)
    w = list(range(1, n_levels + 1)); tot = sum(w)
    s = 0.0; k = len(level_actions)
    for l, a in enumerate(level_actions, start=1):
        h = baselines[l - 1] if l - 1 < len(baselines) else None
        if h is None or a <= 0: continue
        s += w[l - 1] * min(1.15, (h / a) ** 2)
    cap = sum(w[:k]) / tot
    return min(s / tot, cap), cap

def play_one(args):
    spec, gid, max_actions, time_limit, verbose = args
    logging.disable(logging.CRITICAL)
    import arc_agi
    from arcengine import GameAction, GameState
    mod, cls = spec.split(":")
    sys.path.insert(0, os.getcwd())
    AgentCls = getattr(importlib.import_module(mod), cls)
    arc = arc_agi.Arcade()
    env = arc.make(gid)
    info = next(e for e in arc.get_environments() if e.game_id.startswith(gid))
    f = env.reset()
    def obs_of(f):
        return dict(frames=[np.asarray(x) for x in f.frame], state=f.state.name, levels_completed=f.levels_completed,
                    win_levels=f.win_levels, available_actions=list(f.available_actions), full_reset=f.full_reset)
    o = obs_of(f)
    agent = AgentCls(gid, f.win_levels, list(f.available_actions))
    t0 = time.time(); n = 0; level_actions = []; since_level = 0; last_levels = 0; err = None
    try:
        while n < max_actions and time.time() - t0 < time_limit:
            if o["state"] == "WIN": break
            if o["state"] in ("GAME_OVER", "NOT_PLAYED"):
                aid, data = 0, None
                if hasattr(agent, "on_game_over"): agent.on_game_over(o)
            else:
                aid, data = agent.act(o)
            a = GameAction.from_id(aid)
            f = env.step(a, data=data if a.is_complex() else None)
            n += 1; since_level += 1
            if f is None: break
            o = obs_of(f)
            if hasattr(agent, "observe"): agent.observe(o, aid, data)
            if o["levels_completed"] > last_levels:
                for _ in range(o["levels_completed"] - last_levels):
                    level_actions.append(since_level); since_level = 0
                last_levels = o["levels_completed"]
                if verbose: print(f"  {gid}: level {last_levels} at action {n}", flush=True)
    except Exception as e:
        err = traceback.format_exc()
    score, cap = rhae(level_actions, info.baseline_actions or [], f.win_levels if f else len(info.baseline_actions))
    return dict(game=gid, levels=last_levels, win_levels=info.baseline_actions and len(info.baseline_actions), actions=n,
                level_actions=level_actions, baselines=info.baseline_actions, score=round(score*100,2), cap=round(cap*100,1),
                seconds=round(time.time()-t0,1), error=err)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("spec"); p.add_argument("--games", default=None); p.add_argument("--max-actions", type=int, default=2000)
    p.add_argument("--procs", type=int, default=4); p.add_argument("--time-per-game", type=float, default=600)
    p.add_argument("--out", default=None); p.add_argument("-v", action="store_true")
    a = p.parse_args()
    logging.disable(logging.CRITICAL)
    import arc_agi
    arc = arc_agi.Arcade()
    gids = [e.game_id.split("-")[0] for e in arc.get_environments()]
    if a.games: gids = [g for g in a.games.split(",") if g in gids] or a.games.split(",")
    jobs = [(a.spec, g, a.max_actions, a.time_per_game, a.v) for g in gids]
    t = time.time()
    with Pool(a.procs) as pool:
        res = pool.map(play_one, jobs)
    tot = np.mean([r["score"] for r in res])
    print(f"{'game':6} {'lv':>5} {'acts':>6} {'score':>7} {'cap':>6} {'sec':>6}  level_actions / baselines")
    for r in res:
        print(f"{r['game']:6} {r['levels']:>2}/{r['win_levels']:<2} {r['actions']:>6} {r['score']:>7} {r['cap']:>6} {r['seconds']:>6}  {r['level_actions']} / {r['baselines']}" + (f"  ERROR: {r['error'].splitlines()[-1]}" if r['error'] else ""))
    print(f"MEAN RHAE over {len(res)} games: {tot:.2f}%   levels: {sum(r['levels'] for r in res)}   wall {time.time()-t:.0f}s")
    if a.out: json.dump(res, open(a.out, "w"), indent=1)
