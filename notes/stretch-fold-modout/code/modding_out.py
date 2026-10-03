# A Cantor-structured input set x(xi) = sum_k r^(k-1) g_k(xi_k) pushed through a random He ReLU net.
# For pairs whose first differing digit is k, track the chord distance sqrt(2(1-corr)) layer by layer,
# compare with the infinite-width two-point law corr_l = rho^{o l}(c_k), and count PB eps-chain classes.
import numpy as np, itertools
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
rho = lambda c: c/2 + (np.sqrt(np.maximum(1-c*c, 0)) + c*np.arcsin(np.clip(c, -1, 1)))/np.pi
rng = np.random.default_rng(11)
n, L, K, r = 512, 16, 8, 0.55
G = rng.standard_normal((K, 2, n))
codes = np.array(list(itertools.product([0, 1], repeat=K)))
X = sum(r**k * G[k, codes[:, k]] for k in range(K))
X /= np.linalg.norm(X, axis=1, keepdims=True)
first_diff = np.full((len(codes), len(codes)), -1)
for i in range(len(codes)):
    d = codes != codes[i]
    first_diff[i] = np.where(d.any(1), d.argmax(1), -1)
def chord_by_digit(H):
    Hn = H / np.linalg.norm(H, axis=1, keepdims=True)
    C = Hn @ Hn.T
    return np.array([np.sqrt(2*(1 - C[first_diff == k].mean())) for k in range(K)]), np.array([C[first_diff == k].mean() for k in range(K)])
H = X; d0, c0 = chord_by_digit(H)
rows = [("input", d0, d0)]
Ws = [rng.standard_normal((n, n)) * np.sqrt(2.0/n) for _ in range(L)]
pred_c = c0.copy()
for l, W in enumerate(Ws, 1):
    H = np.maximum(H @ W.T, 0); pred_c = rho(pred_c)
    if l in (1, 4, 8, 16):
        d, _ = chord_by_digit(H); rows.append((f"layer {l}", d, np.sqrt(2*(1-pred_c))))
print("chord distance between inputs whose first differing digit is k (k=1 coarsest .. 8 finest)")
print("            " + "  ".join(f"k={k+1:<5d}" for k in range(K)))
for name, d, pdist_ in rows:
    print(f"{name:>9} obs " + "  ".join(f"{v:7.4f}" for v in d))
    if name != "input": print(f"{'':>9} pred" + "  ".join(f"{v:7.4f}" for v in pdist_))
# eps-chain (single-linkage) classes at a fixed scale, per layer
H = X; eps = 0.6 * d0[2]
out = []
for l in range(L+1):
    if l > 0: H = np.maximum(H @ Ws[l-1].T, 0)
    Hn = H / np.linalg.norm(H, axis=1, keepdims=True)
    Z = linkage(pdist(Hn), method='single')
    if l in (0, 1, 2, 4, 8, 12, 16): out.append(f"l={l}: {len(set(fcluster(Z, eps, criterion='distance')))}")
print(f"number of eps-chain classes (eps = {eps:.3f}) among {len(codes)} inputs: " + ", ".join(out))
