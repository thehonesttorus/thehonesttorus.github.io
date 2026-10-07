# Pairwise cross-term changes of a comparison from fentry.py's noise-free route Grams (note XXXVI section 3h).
#   python pairx.py FENTRY_LOG
# n Delta MSE = sum_X Delta|delta_X|^2 + sum_(X != Y) Delta<delta_X, delta_Y>; with the Grams printed per run in units of
# that run's |e|^2, the change in units of the baseline |e_0|^2 is G^r |e_r|^2 / |e_0|^2 - G^0. Prints the diagonal (own)
# and the symmetric pair terms 2 Delta<delta_X, delta_Y> for X < Y, largest first.
import sys, re
txt = open(sys.argv[1]).read().splitlines()
R = ("k3", "k4", "D21", "K22", "K31", "resid", "rep")
runs = []
i = 0
while i < len(txt):
    m = re.match(r"run (\d+) \((\S+)\): n MSE ([0-9.e+-]+)", txt[i])
    if m:
        E = float(m.group(3)); G = None
        for j in range(i, min(i + 30, len(txt))):
            if "noise-free Gram" in txt[j]:
                G = [[float(x) for x in txt[j + 1 + a].split(":")[1].split()] for a in range(len(R))]
                break
        runs.append((m.group(2), E, G))
    i += 1
b, E0, G0 = runs[0]
for name, E, G in runs[1:]:
    s = E / E0
    print(f"{b} -> {name}: noise-free changes, % of baseline n MSE")
    own = {R[a]: 100 * (G[a][a] * s - G0[a][a]) for a in range(len(R))}
    pairs = {(R[a], R[c]): 200 * (G[a][c] * s - G0[a][c]) for a in range(len(R)) for c in range(a + 1, len(R))}
    print("  own (diagonal): " + "  ".join(f"{k} {v:+.2f}" for k, v in own.items()))
    print(f"  sum own {sum(own.values()):+.2f}, sum pairs {sum(pairs.values()):+.2f}, total {sum(own.values()) + sum(pairs.values()):+.2f} "
          f"(actual {100 * (s - 1):+.2f})")
    print("  pairs 2 Delta<delta_X, delta_Y>, |.| >= 0.5: " + "  ".join(f"{x}-{y} {v:+.2f}" for (x, y), v in
                                                                 sorted(pairs.items(), key=lambda t: -abs(t[1])) if abs(v) >= 0.5))
