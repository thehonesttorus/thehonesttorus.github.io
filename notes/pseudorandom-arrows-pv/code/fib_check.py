# Golden-mean process: continuation classes, Bratteli incidence, KMS_beta recursion, Vershik commutator defect.
import numpy as np, itertools
p = 0.6                     # P(0 | last=0); after a 1 the next symbol is 0
N = 14                      # depth of truncation
def prob(word):
    pr, last = 1.0, 0
    for b in word:
        if last == 1: pr *= (1.0 if b == 0 else 0.0)
        else: pr *= (p if b == 0 else 1 - p)
        last = b
    return pr
words = {n: [w for w in itertools.product([0,1], repeat=n) if prob(w) > 0] for n in range(N+1)}
H = 5                        # horizon used to compare continuation laws exactly (conditional law of next H symbols)
def cont_law(w):
    pw = prob(w)
    return tuple(round(prob(w + z) / pw, 12) for z in itertools.product([0,1], repeat=H))
classes = {}
for n in range(0, N - H + 1):
    groups = {}
    for w in words[n]: groups.setdefault(cont_law(w), []).append(w)
    classes[n] = list(groups.values())
print("number of continuation classes |P_n|:", [len(classes[n]) for n in classes])
print("class sizes at n=1..8:", [sorted(len(C) for C in classes[n]) for n in range(1, 9)])
# incidence between levels n -> n+1
def cls_index(n, w):
    for i, C in enumerate(classes[n]):
        if w in C: return i
for n in [3, 6]:
    A = np.zeros((len(classes[n+1]), len(classes[n])), int)
    for j, C in enumerate(classes[n]):
        for b in [0, 1]:
            w = C[0] + (b,)
            if prob(w) > 0: A[cls_index(n+1, w), j] += 1
    print(f"incidence level {n}->{n+1} (rows: classes at n+1, ordered by size):\n", A)
# KMS_beta weights: x_C = sum_b nu_C(b)^beta x_{C.b}; stationary => x^{(n)} = r^{-n} v with M(beta) v = r v
for beta in [0.0, 0.5, 1.0, 2.0]:
    M = np.array([[p**beta, (1-p)**beta], [1.0, 0.0]])       # states (A: last 0, B: last 1)
    ev, V = np.linalg.eig(M); k = np.argmax(ev.real); v = np.abs(V[:, k].real)
    print(f"beta={beta}: Perron root r(beta)={ev[k].real:.6f}  (beta=0: golden ratio {(1+5**.5)/2:.6f}; beta=1: 1)  log r = pressure-type term {np.log(ev[k].real):+.4f}")
