import numpy as np
from hd import hd, closure
def pred_closure(W): return closure(np.asarray(W, np.float64))
def pred_hd(W): return hd(np.asarray(W, np.float64), A=None)
def pred_hd2(W): return hd(np.asarray(W, np.float64), A=None, second=True)
def pred_hd2full(W): return hd(np.asarray(W, np.float64), A=None, second="full")
