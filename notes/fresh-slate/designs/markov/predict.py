import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mkv import mkv
def full(W): return mkv(W, w=16, var21='full')
def w2(W): return mkv(W, w=2, var21='full')
def gauss(W): return mkv(W, gauss=True, ng_cov=False, var21=False)
