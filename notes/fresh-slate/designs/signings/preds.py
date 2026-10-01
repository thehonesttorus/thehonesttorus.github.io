import numpy as np, copula, copula1
def v0(W): return copula.estimate(np.asarray(W, float))
def v1(W): return copula1.estimate(np.asarray(W, float), {})
def v2(W): return copula1.estimate(np.asarray(W, float), {'src': 1})
def v1pq1(W):
    copula1_PQ = copula1.PQ
    import copula as c; c.PQ = 1; copula1.PQ = 1
    try: return copula1.estimate(np.asarray(W, float), {})
    finally: c.PQ = 3; copula1.PQ = copula1_PQ
