"""Process-per-game engine proxy.

The ARC-AGI-3-Agents framework (which the Kaggle gateway drives) runs every game in its
own *thread* of one Python process. With 110 hidden games and a pure-Python engine, the
GIL would serialize all games onto one core (9 h / 110 games is ~5 CPU-minutes each).
`EngineProxy` runs the real engine in a child process and exposes the same protocol, so
the swarm's threads only shuttle observations and actions while the 48 vCPUs do the work.

Protocol mirrored (see nce/engine.py):
    act(obs) -> (action_id, data | None)
    observe(obs, action_id, data)
    on_game_over(obs)
    want_stop(obs) -> bool
    last_reason : str

Failure policy: see EngineProxy. A learned engine is never silently replaced after one bad call.
"""
from __future__ import annotations

import multiprocessing as mp
import os
import sys
import traceback
from typing import Any, Optional

import numpy as np


def _pack(obs: dict[str, Any]) -> dict[str, Any]:
    # frames as contiguous int8 arrays pickle cheaply; everything else is plain.
    return dict(
        frames=[np.ascontiguousarray(np.asarray(f, dtype=np.int8)) for f in obs["frames"]],
        state=str(obs["state"]),
        levels_completed=int(obs["levels_completed"]),
        win_levels=int(obs.get("win_levels", 0)),
        available_actions=[int(a) for a in obs.get("available_actions", [])],
        full_reset=bool(obs.get("full_reset", False)),
    )


def _worker(conn, engine_spec: tuple[str, str], init_args: tuple) -> None:  # pragma: no cover (child)
    """Serve one engine. Every message is handled in its own try block: an error is reported for that message and
    the worker keeps its state and keeps serving (a learned engine is never thrown away for one bad call)."""
    try:
        mod, cls = engine_spec
        pkg_dir = os.environ.get("NCE_PKG_DIR")
        if pkg_dir and pkg_dir not in sys.path:
            sys.path.insert(0, pkg_dir)
        import importlib
        Engine = getattr(importlib.import_module(mod), cls)
        engine = Engine(*init_args)
        conn.send((0, "ready", None))
    except Exception:
        try:
            conn.send((0, "error", traceback.format_exc()))
        except Exception:
            pass
        return
    while True:
        try:
            msg = conn.recv()
        except (EOFError, OSError):
            return
        seq, kind = msg[0], msg[1]
        try:
            if kind == "act":
                out = engine.act(msg[2])
                conn.send((seq, "ok", (out, getattr(engine, "last_reason", ""))))
            elif kind == "observe":
                if hasattr(engine, "observe"):
                    engine.observe(msg[2], msg[3], msg[4])
                conn.send((seq, "ok", None))
            elif kind == "resync":
                if hasattr(engine, "resync"):
                    engine.resync(msg[2])
                conn.send((seq, "ok", None))
            elif kind == "game_over":
                if hasattr(engine, "on_game_over"):
                    engine.on_game_over(msg[2])
                conn.send((seq, "ok", None))
            elif kind == "want_stop":
                conn.send((seq, "ok", bool(engine.want_stop(msg[2])) if hasattr(engine, "want_stop") else False))
            elif kind == "close":
                conn.send((seq, "ok", None))
                return
            else:
                conn.send((seq, "error", f"unknown message {kind!r}"))
        except Exception:
            try:
                conn.send((seq, "error", traceback.format_exc()))
            except Exception:
                return


class EngineProxy:
    """Engine in a child process with the in-process protocol. Failure policy:
    - an error inside one call is reported and the worker keeps its learned state;
    - a timeout returns None for that call (the caller uses its floor action) and later stale replies are discarded
      by sequence number, so the pipe never desynchronises;
    - a dead worker is replaced by a fresh worker process (state is lost, parallelism is kept); only if no process
      can be started does the game run in-process."""

    def __init__(self, game_id: str, win_levels: int, available_actions: list[int],
                 engine_spec: tuple[str, str] = ("nce.engine", "Engine"), timeout: Optional[float] = None,
                 use_process: Optional[bool] = None):
        self.game_id = game_id
        self.last_reason = ""
        self._spec = engine_spec
        self._init = (game_id, win_levels, list(available_actions))
        self._timeout = float(os.environ.get("NCE_CALL_TIMEOUT", "300")) if timeout is None else float(timeout)
        self._local = None
        self._conn = None
        self._proc = None
        self._seq = 0
        self.restarts = 0
        if use_process is None:
            use_process = os.environ.get("NCE_USE_PROCESS", "1") != "0"
        self._use_process = use_process
        if use_process:
            self._start_worker()
        if self._conn is None:
            self._local = self._make_local()

    def _start_worker(self) -> None:
        try:
            ctx = mp.get_context("spawn")
            parent, child = ctx.Pipe()
            self._proc = ctx.Process(target=_worker, args=(child, self._spec, self._init), daemon=True)
            self._proc.start()
            self._conn = parent
            startup = float(os.environ.get("NCE_PROCESS_STARTUP_TIMEOUT", "90"))
            if not parent.poll(startup):
                raise TimeoutError("engine worker did not start")
            _, kind, payload = parent.recv()
            if kind != "ready":
                raise RuntimeError(payload)
        except Exception as e:
            print(f"[nce] process engine unavailable for {self.game_id} ({e!r}); running in-process", file=sys.stderr)
            self._shutdown()

    def _make_local(self):
        import importlib
        mod, cls = self._spec
        return getattr(importlib.import_module(mod), cls)(*self._init)

    def _shutdown(self) -> None:
        try:
            if self._conn is not None:
                self._conn.close()
        except Exception:
            pass
        try:
            if self._proc is not None:
                if self._proc.is_alive():
                    self._proc.terminate()
                self._proc.join(timeout=5)
        except Exception:
            pass
        self._conn = None
        self._proc = None

    def _call(self, kind: str, *args: Any) -> tuple[bool, Any]:
        """(ok, payload). ok is False on error or timeout; the worker is kept unless it died."""
        if self._conn is None:
            return False, None
        self._seq += 1
        seq = self._seq
        try:
            self._conn.send((seq, kind) + tuple(args))
            while True:
                if not self._conn.poll(self._timeout):
                    print(f"[nce] worker timeout on {self.game_id} ({kind}); continuing with the floor action", file=sys.stderr)
                    return False, None
                rseq, rkind, payload = self._conn.recv()
                if rseq != seq:
                    continue                      # a stale reply to a call that had timed out
                if rkind == "error":
                    print(f"[nce] worker error on {self.game_id} ({kind}): {payload}", file=sys.stderr)
                    return False, None
                return True, payload
        except (EOFError, OSError, BrokenPipeError) as e:
            print(f"[nce] worker died on {self.game_id} ({e!r}); starting a fresh worker", file=sys.stderr)
            self._shutdown()
            self.restarts += 1
            if self._use_process:
                self._start_worker()
            if self._conn is None and self._local is None:
                self._local = self._make_local()
            return False, None

    @staticmethod
    def _floor(obs: dict[str, Any]) -> tuple[int, Optional[dict[str, Any]]]:
        avail = [int(a) for a in obs.get("available_actions", []) if int(a) != 0] or [1]
        a = avail[0]
        return a, ({"x": 0, "y": 0} if a == 6 else None)

    # ---- protocol -------------------------------------------------------------
    def act(self, obs: dict[str, Any]) -> tuple[int, Optional[dict[str, Any]]]:
        if self._conn is not None:
            ok, out = self._call("act", _pack(obs))
            if ok and out is not None:
                (aid, data), self.last_reason = out
                return int(aid), data
            if self._conn is not None:
                self.last_reason = "floor"
                return self._floor(obs)
        aid, data = self._local.act(obs)
        self.last_reason = getattr(self._local, "last_reason", "")
        return int(aid), data

    def _forward(self, kind: str, method: str, obs: dict[str, Any], *args: Any) -> None:
        if self._conn is not None:
            self._call(kind, _pack(obs), *args)
            if self._conn is not None:
                return
        if self._local is not None and hasattr(self._local, method):
            getattr(self._local, method)(obs, *args)

    def observe(self, obs: dict[str, Any], aid: int, data: Optional[dict[str, Any]]) -> None:
        self._forward("observe", "observe", obs, int(aid), data)

    def resync(self, obs: dict[str, Any]) -> None:
        self._forward("resync", "resync", obs)

    def on_game_over(self, obs: dict[str, Any]) -> None:
        self._forward("game_over", "on_game_over", obs)

    def want_stop(self, obs: dict[str, Any]) -> bool:
        if self._conn is not None:
            ok, out = self._call("want_stop", _pack(obs))
            if self._conn is not None:
                return bool(out) if ok else False
        if self._local is not None and hasattr(self._local, "want_stop"):
            return bool(self._local.want_stop(obs))
        return False

    def close(self) -> None:
        if self._conn is not None:
            try:
                self._call("close")
            finally:
                self._shutdown()
