"""G3 cost predictions at 256 x 32 vs 1024/256 x 16, in units of 'full-resolution source-steps' (one source carried at
rank n for one layer = 1; FC = all pairs).  Profiles: (A) critical k(a) = min(1, 2/a) over all ages (team D's law),
(B) geometric / fixed window w (team G's subcritical profile), (C) measured: window w*(t) needed for <= 10 % loss
(g3_c2w*.json, layers <= 24) times 2n/a resolution.  Also the wedge price weighting at L = 32."""
import numpy as np, json
exec(open('wedge_check_any.py').read().split("reg = ")[0])
def cost(L, kfun, win=lambda t: 10**9):
    return sum(sum(kfun(a) for a in range(1, min(t, win(t)) + 1)) for t in range(1, L))
full = lambda a: 1.0
crit = lambda a: 1.0 if a <= 2 else min(1.0, 2.0 / a)
for L in (16, 32):
    print(f'L={L}: FC {cost(L, full):.0f}; critical 2n/a {cost(L, crit):.1f}; fixed window 8 x 2n/a {cost(L, crit, lambda t: 8):.1f};'
          f' fixed window 12 {cost(L, crit, lambda t: 12):.1f}; measured window 0.75t {cost(L, crit, lambda t: max(4, int(round(0.75 * t)))):.1f}')
for name, kf, wf in [('FC', full, lambda t: 10**9), ('critical', crit, lambda t: 10**9), ('window 8', crit, lambda t: 8), ('0.75 t', crit, lambda t: max(4, int(round(0.75 * t))))]:
    r = cost(32, kf, wf) / cost(16, kf, wf)
    print(f'  ratio cost(32)/cost(16) {name}: {r:.2f}   (L^2 -> 4.1, L log L -> 2.5, L -> 2.07)')
Km, Ko = wedge(256, 32); w = Ko[:31] / Ko[:31].sum(); c = np.cumsum(w[::-1])[::-1]
first = int(np.argmax(c <= 0.95)); print(f'wedge price at 256x32: layers >= {first} carry 95 % of the covariance price; layers < 8 carry {w[:8].sum():.3f}')
ptw = lambda t: 10**9 if t >= first else 2
print(f'  price-weighted critical carrier (all ages only at targets >= {first}, 2 ages elsewhere): {cost(32, crit, ptw):.1f} vs critical {cost(32, crit):.1f}')
