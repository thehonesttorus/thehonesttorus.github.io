"""Shared helpers for the interpolation stream (numpy, float64 closures)."""
import os, sys
import numpy as np
from scipy.special import ndtr, gammaln
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "bench"))
import bench  # noqa

SQ2PI = np.sqrt(2 * np.pi)
def phi(a): return np.exp(-0.5 * a * a) / SQ2PI
def Phi(a): return ndtr(a)

def chi_mean_ratio(n):
    """E r / sqrt(E r^2) for r ~ chi_n."""
    return np.exp(0.5 * np.log(2.0 / n) + gammaln((n + 1) / 2) - gammaln(n / 2))

def mc_layers(W, N, seed=0, chunk=4096, stats=None):
    """Monte Carlo through the net, float32; returns per-layer sample means and calls stats(l, A) per chunk."""
    rng = np.random.default_rng(seed)
    L, n, _ = W.shape
    acc = np.zeros((L, n))
    done = 0
    while done < N:
        m = min(chunk, N - done)
        h = rng.standard_normal((m, n)).astype(np.float32)
        x = h
        for l in range(L):
            h = np.maximum(h @ W[l], 0.0)
            acc[l] += h.sum(0, dtype=np.float64)
            if stats is not None:
                stats(l, h, x)
        done += m
    return acc / N
