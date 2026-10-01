"""P5: size of the sign cancellation in the tropical (Newton-polytope) form a_l = max(P_l, Q_l) - Q_l.

Gaussian mean widths are additive under Minkowski sums, so with w_P = E P, w_Q = E Q (exact, no closure):
  w_Q'(j) = sum_i W+_ij w_Q(i) + |W-_ij| w_P(i),   w_P'(j) = sum_i W+_ij w_P(i) + |W-_ij| w_Q(i)
and the ReLU sets P_l = max(P', Q'), Q_l = Q', so w_P = w_Q' + m_l with m_l = E a_l (O(1), taken from TCT-0).
Layer 0: P = x_i (a segment-free linear form: w = 0), Q = 0.  Report median w_Q / m by layer."""
import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__file__))
from tct0 import tct0
for n in (256, 1024):
    rng = np.random.default_rng(0); L = 16
    Ws = [rng.standard_normal((n, n)) * np.sqrt(2 / n) for _ in range(L)]
    m = tct0(Ws) if n <= 256 else None
    wP = np.zeros(n); wQ = np.zeros(n)  # x_i = h_{point e_i}: mean width 0
    out = []
    for l, W in enumerate(Ws):
        Wp = np.maximum(W, 0); Wm = np.maximum(-W, 0)
        wQn = wP @ Wm + wQ @ Wp
        ml = m[l] if m is not None else np.full(n, 0.55)
        wP, wQ = wQn + ml, wQn
        out.append(np.median(wQ / ml))
    print(f"n={n}: median E[Q_l]/E[a_l] by layer:", " ".join(f"{v:.1e}" for v in out))
