# Stream old-content: a cheap carrier for the transported old-source content of κ3

Status: in progress (tracker and carrier harness written and validated on a small test atlas; width-128 N = 6e5 atlases building).

## Question
The published K = 3 chain carries the third cumulant as per-layer sources, each transported with dense n×n legs; its
old-source tier costs 107 of 260 units. Leaders spend ~7–10 units per layer in total. Is there a carrier of the old
content's effect on D21(l+1) = κ3(z_a, z_a, z_b) reaching ε ≤ 2–3 % at ≤ 10 units per layer (n = 1024)?

## Method
- `tracker.py`: exact sample-level per-source decomposition of κ3(z_l) at width 128 (dense n³ tensors).
- `carriers.py`: carriers (a) covariance-response mode families, (b) CP / hub columns, (c) regression on O(n²) objects,
  (d) merged aggregate source (CP-dyn), (e) structure-suggested families; static (oracle refit) and dynamic (compounding).

Results, verdict: pending.
