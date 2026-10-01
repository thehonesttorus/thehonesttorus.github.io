"""Time the grader's 5 s setup window exactly as SubprocessRunner sees it: from spawning the
worker (interpreter start + imports + estimator load + setup()) to the start response.

    python setup_time.py ESTIMATOR.py [repeats] [--out results.json]
"""
import json, os, sys, time
from pathlib import Path

from whestbench.runner import SubprocessRunner, ResourceLimits
from whestbench.sdk import SetupContext
from whestbench.runner import EstimatorEntrypoint


def main():
    est = os.path.abspath(sys.argv[1])
    reps = int(sys.argv[2]) if len(sys.argv) > 2 and not sys.argv[2].startswith("--") else 5
    out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else None
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
    times = []
    for _ in range(reps):
        r = SubprocessRunner()
        ctx = SetupContext(width=1024, depth=16, flop_budget=2 ** 41, api_version="1.0",
                           submission_dir=str(Path(est).parent), seed=0)
        lim = ResourceLimits(setup_timeout_s=60.0, predict_timeout_s=30.0, memory_limit_mb=8192, flop_budget=2 ** 41)
        t = time.perf_counter()
        r.start(EstimatorEntrypoint(file_path=Path(est), class_name="Estimator"), ctx, lim)
        times.append(time.perf_counter() - t)
        r.close()
    rec = {"estimator": est, "setup_window_s": times, "max_s": max(times), "cap_s": 5.0}
    print(json.dumps(rec))
    if out:
        json.dump(rec, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
