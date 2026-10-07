import numpy as np
rng = np.random.default_rng(5); n, L = 7, 4
W = rng.standard_normal((L, n, n)) * np.sqrt(2 / n)
W[0] = rng.standard_normal((n, n)) / np.sqrt(n)
np.save("official/W_off9.npy", W.astype(np.float32))
