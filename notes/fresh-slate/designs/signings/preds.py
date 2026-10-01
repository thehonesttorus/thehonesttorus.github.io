import numpy as np, copula, copula1
def v0(W, cfg=None): return copula.estimate(np.asarray(W, float))
def v1(W, cfg=None): return copula1.estimate(np.asarray(W, float), {})
def v2(W, cfg=None): return copula1.estimate(np.asarray(W, float), {'src': 1})
def v1pq1(W, cfg=None):
    copula1_PQ = copula1.PQ
    import copula as c; c.PQ = 1; copula1.PQ = 1
    try: return copula1.estimate(np.asarray(W, float), {})
    finally: c.PQ = 3; copula1.PQ = copula1_PQ
import sys as _s, os as _o
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '../../bench'))
def gauss(W, cfg=None):
    import eval_q; return eval_q.baseline_gauss(np.asarray(W, np.float32))
