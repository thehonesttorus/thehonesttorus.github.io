import numpy as np

from nce.proxy import EngineProxy


def _obs():
    return dict(frames=[np.zeros((64, 64), dtype=np.int8)], state="NOT_FINISHED", levels_completed=0,
                win_levels=6, available_actions=[1, 2, 3, 4, 6], full_reset=False)


def test_process_proxy_round_trip_and_close():
    p = EngineProxy("t", 6, [1, 2, 3, 4, 6])
    assert p._conn is not None
    for _ in range(20):
        aid, data = p.act(_obs())
        assert aid in (1, 2, 3, 4, 6)
        assert (data is None) if aid != 6 else (0 <= data["x"] <= 63 and 0 <= data["y"] <= 63)
        p.observe(_obs(), aid, data)
    assert p.want_stop(_obs()) is False
    p.on_game_over(_obs())
    p.close()
    assert p._conn is None


def test_in_process_mode_and_bad_spec_raises():
    p = EngineProxy("t2", 6, [1], use_process=False)
    assert p._conn is None and p.act(_obs())[0] in (1, 2, 3, 4, 6)
    try:
        EngineProxy("t3", 6, [1], engine_spec=("nce.nonexistent", "Engine"))
        raised = False
    except ModuleNotFoundError:
        raised = True
    assert raised
