import numpy as np
from hd import hd, closure
def pred_closure(W): return closure(np.asarray(W, np.float64))
def pred_hd(W): return hd(np.asarray(W, np.float64), A=None)
def pred_hd2(W): return hd(np.asarray(W, np.float64), A=None, second=True)
def pred_hd2full(W): return hd(np.asarray(W, np.float64), A=None, second="full")
def pred_hd_star(W): return hd(np.asarray(W, np.float64), A=None, diagrams="star")
def pred_hd_startri(W): return hd(np.asarray(W, np.float64), A=None, diagrams="startri")
def pred_hd_A3(W): return hd(np.asarray(W, np.float64), A=3)
def pred_hd_A7(W): return hd(np.asarray(W, np.float64), A=7)
from hd import hd2
def pred_hd2k4(W): return hd2(np.asarray(W, np.float64))
