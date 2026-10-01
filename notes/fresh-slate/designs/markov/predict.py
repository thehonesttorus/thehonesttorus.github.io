import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mkv import mkv
def full(W): return mkv(W, w=16, var21='full')
def w2(W): return mkv(W, w=2, var21='full')
def gauss(W): return mkv(W, gauss=True, ng_cov=False, var21=False)
from mkv2 import mkv2
def w4(W): return mkv(W, w=4, var21='full')
def c16(W): return mkv2(W, w=16)
def c4(W): return mkv2(W, w=4)
