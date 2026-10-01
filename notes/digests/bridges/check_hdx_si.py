"""Numerical sanity checks for notes/digests/bridges/hdx-spectral-independence.md.

Pure numpy. Each check corresponds to a labelled statement in the digest:
  C1  ALO Theorem 3.1: spectrum of the 1-skeleton walk P_0 of X_mu equals
      spec(Psi_mu/(n-1)) + (n-1) copies of -1/(n-1) + {1}.
  C2  Markov-projection identity: D(P || Q_Markov) = sum_i I(X_{i+1}; X_{<i} | X_i).
  C3  Commuting square at a layer: E_{<=l} E_{>=l} = E_{l} iff past _||_ future | present.
  C4  Classical recovery: heat-bath restricted to A, run to stationarity, maps
      mu_{A^c} (x) unif_A back to mu exactly; it acts only through sigma_{boundary(A)}
      when mu is Markov; the time-averaged map has Dirichlet form <= c/t.
  C5  Alev-Lau bound and ALO Theorem 1.3 bound versus the true Glauber gap.
Run:  python3 check_hdx_si.py
"""
import itertools
import numpy as np

rng = np.random.default_rng(0)


def subsets(n):
    return [tuple(s) for s in itertools.product([0, 1], repeat=n)]


# ---------------------------------------------------------------- C1
def check_alo_spectrum(n=5, trials=3):
    worst = 0.0
    for _ in range(trials):
        conf = subsets(n)
        w = rng.random(len(conf)) ** 3          # a generic full-support law on {0,1}^n
        mu = dict(zip(conf, w / w.sum()))
        P1 = np.array([sum(p for s, p in mu.items() if s[i]) for i in range(n)])
        # pairwise probabilities Pr[x_i = a, x_j = b]
        def pr(i, a, j, b):
            return sum(p for s, p in mu.items() if s[i] == a and s[j] == b)
        # influence matrix Psi(i,j) = Pr[j|i] - Pr[j|not i]
        Psi = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i != j:
                    Psi[i, j] = pr(i, 1, j, 1) / P1[i] - pr(i, 0, j, 1) / (1 - P1[i])
        # 1-skeleton walk of X_mu on 2n vertices (i, a): weight of edge {(i,a),(j,b)} = Pr[x_i=a, x_j=b]
        V = [(i, a) for i in range(n) for a in (1, 0)]
        W = np.zeros((2 * n, 2 * n))
        for u, (i, a) in enumerate(V):
            for v, (j, b) in enumerate(V):
                if i != j:
                    W[u, v] = pr(i, a, j, b)
        P0 = W / W.sum(axis=1, keepdims=True)
        ev_P = np.sort(np.linalg.eigvals(P0).real)
        ev_pred = np.sort(np.concatenate([np.linalg.eigvals(Psi).real / (n - 1),
                                          -np.ones(n - 1) / (n - 1), [1.0]]))
        worst = max(worst, np.max(np.abs(ev_P - ev_pred)))
        lam2 = np.sort(ev_P)[-2]
        lmax = np.max(np.linalg.eigvals(Psi).real)
        assert abs(lam2 - max(lmax / (n - 1), -1 / (n - 1))) < 1e-9
    return worst


# ---------------------------------------------------------------- C2, C3
def markov_projection(P):
    """P: array of shape (q,)*L joint law. Returns Q(x)=P(x0) prod P(x_{i+1}|x_i)."""
    L = P.ndim
    Q = np.ones_like(P)
    idx = np.indices(P.shape)
    m0 = P.sum(axis=tuple(range(1, L)))
    Q *= m0[idx[0]]
    for i in range(L - 1):
        pair = P.sum(axis=tuple(k for k in range(L) if k not in (i, i + 1)))
        cond = pair / pair.sum(axis=1, keepdims=True)
        Q *= cond[idx[i], idx[i + 1]]
    return Q


def H(p):
    p = p[p > 0]
    return -(p * np.log(p)).sum()


def cmi_past(P, i):
    """I(X_{i+1}; X_{<i} | X_i) in nats."""
    L = P.ndim
    keep = tuple(range(i + 2))
    Pm = P.sum(axis=tuple(range(i + 2, L))) if i + 2 < L else P
    # H(X_{i+1}|X_i) - H(X_{i+1}|X_{<=i})
    p_i_ip1 = Pm.sum(axis=tuple(range(i)))
    p_i = p_i_ip1.sum(axis=1)
    p_le_i = Pm.sum(axis=i + 1)
    return (H(p_i_ip1.ravel()) - H(p_i)) - (H(Pm.ravel()) - H(p_le_i.ravel()))


def check_markov_projection(q=3, L=4):
    P = rng.random((q,) * L) ** 2
    P /= P.sum()
    Q = markov_projection(P)
    D = (P * np.log(P / Q)).sum()
    S = sum(cmi_past(P, i) for i in range(1, L - 1))
    return D, S


def check_commuting_square(q=3):
    # three layers X0 - X1 - X2 ; E_{<=1} = E[. | X0,X1], E_{>=1} = E[. | X1,X2], E_1 = E[. | X1]
    def defect(P):
        f = rng.standard_normal((q, q, q))
        Ege = (P * f).sum(axis=0, keepdims=True) / P.sum(axis=0, keepdims=True)       # E[f|X1,X2]
        Ege = np.broadcast_to(Ege, P.shape)
        Ele_ge = (P * Ege).sum(axis=2, keepdims=True) / P.sum(axis=2, keepdims=True)  # E[E[f|X1X2]|X0X1]
        E1 = (P * f).sum(axis=(0, 2), keepdims=True) / P.sum(axis=(0, 2), keepdims=True)
        return np.max(np.abs(np.broadcast_to(Ele_ge, P.shape) - np.broadcast_to(E1, P.shape)))
    P = rng.random((q, q, q)); P /= P.sum()
    Q = markov_projection(P)
    return defect(P), defect(Q)


# ---------------------------------------------------------------- C4
def ising_chain(n, beta, h):
    conf = list(itertools.product([-1, 1], repeat=n))
    w = np.array([np.exp(beta * sum(s[i] * s[i + 1] for i in range(n - 1)) + sum(h[i] * s[i] for i in range(n)))
                  for s in conf])
    return conf, w / w.sum()


def heat_bath_generator(conf, mu, sites):
    """Generator L = sum_{a in sites} (E_a - I) of heat-bath updates on the given sites (Schrodinger picture on laws)."""
    N = len(conf)
    index = {s: k for k, s in enumerate(conf)}
    L = np.zeros((N, N))
    for k, s in enumerate(conf):
        for a in sites:
            s_flip = list(s); s_flip[a] = -s[a]; kf = index[tuple(s_flip)]
            Z = mu[k] + mu[kf]
            L[k, k] += mu[k] / Z - 1.0          # stay with prob mu(s)/Z
            L[k, kf] += mu[kf] / Z               # move with prob mu(s')/Z
    return L                                     # row-stochastic generator: law' = law @ expm(tL)


def expm(M, t):
    w, V = np.linalg.eig(M)
    return (V @ np.diag(np.exp(t * w)) @ np.linalg.inv(V)).real


def check_recovery(n=6, beta=0.8):
    h = rng.standard_normal(n) * 0.3
    conf, mu = ising_chain(n, beta, h)
    A = [2, 3]
    LA = heat_bath_generator(conf, mu, A)
    # nu = mu_{A^c} (x) unif_A
    index = {s: k for k, s in enumerate(conf)}
    nu = np.zeros(len(conf))
    for k, s in enumerate(conf):
        key = tuple(x for i, x in enumerate(s) if i not in A)
        marg = sum(mu[index[t]] for t in conf if tuple(x for i, x in enumerate(t) if i not in A) == key)
        nu[k] = marg / 2 ** len(A)
    rec = nu @ expm(LA, 200.0)
    err_long = np.abs(rec - mu).sum()
    # time average R_t = (1/t) int_0^t e^{sL} ds, Dirichlet form of R_t f for |f| <= 1
    # (Heisenberg picture: f -> expm(tL) f ; L is mu-reversible)
    ts = [1.0, 4.0, 16.0, 64.0]
    D = np.diag(mu)
    out = []
    for t in ts:
        grid = np.linspace(0, t, 801)
        Rt = sum(expm(LA, s) for s in grid) / len(grid)
        worst = 0.0
        for _ in range(50):
            f = rng.choice([-1.0, 1.0], size=len(conf))
            g = Rt @ f
            dir_form = -(g @ D @ LA @ g)
            worst = max(worst, dir_form)
        out.append((t, worst, worst * t))
    # locality: the limit map is E_A, which depends only on sigma at sites adjacent to A
    EA = expm(LA, 200.0)
    loc_ok = True
    for k, s in enumerate(conf):
        for k2, s2 in enumerate(conf):
            if all(s[i] == s2[i] for i in range(n) if i in (1, 4)) and all(s[i] == s2[i] for i in A):
                # same boundary {1,4} and same A: transition rows to A-configurations must agree
                for a_conf in itertools.product([-1, 1], repeat=len(A)):
                    p1 = sum(EA[k, index[t]] for t in conf if tuple(t[i] for i in A) == a_conf)
                    p2 = sum(EA[k2, index[t]] for t in conf if tuple(t[i] for i in A) == a_conf)
                    if abs(p1 - p2) > 1e-8:
                        loc_ok = False
    return err_long, out, loc_ok


# ---------------------------------------------------------------- C5
def glauber_gap_vs_bounds(n=6, beta=0.4):
    h = rng.standard_normal(n) * 0.2
    conf, mu = ising_chain(n, beta, h)
    LG = heat_bath_generator(conf, mu, list(range(n))) / n + np.eye(len(conf))   # discrete-time Glauber
    ev = np.sort(np.linalg.eigvals(LG).real)
    gap = 1 - ev[-2]
    # spectral independence eta_i: max over pinnings of i sites of lambda_max(Psi) (signed, binary)
    index = {s: k for k, s in enumerate(conf)}
    etas = []
    for i in range(n - 1):
        worst = -np.inf
        for pinned in itertools.combinations(range(n), i):
            for vals in itertools.product([-1, 1], repeat=i):
                sub = [(s, mu[index[s]]) for s in conf if all(s[p] == v for p, v in zip(pinned, vals))]
                free = [j for j in range(n) if j not in pinned]
                Z = sum(p for _, p in sub)
                def pr(cond):
                    return sum(p for s, p in sub if cond(s)) / Z
                m = len(free)
                Psi = np.zeros((m, m))
                for a, u in enumerate(free):
                    pu = pr(lambda s: s[u] == 1)
                    for b, v in enumerate(free):
                        if u != v:
                            Psi[a, b] = (pr(lambda s: s[u] == 1 and s[v] == 1) / pu
                                         - pr(lambda s: s[u] == -1 and s[v] == 1) / (1 - pu))
                worst = max(worst, np.max(np.linalg.eigvals(Psi).real))
        etas.append(worst)
    alo = (1 / n) * np.prod([1 - etas[i] / (n - i - 1) for i in range(n - 1)])
    return gap, alo, etas


if __name__ == "__main__":
    print("C1  ALO Thm 3.1 spectrum identity, max |difference| over 3 random laws on {0,1}^5:",
          f"{check_alo_spectrum():.2e}")
    D, S = check_markov_projection()
    print(f"C2  D(P||Markov projection) = {D:.10f} ;  sum of CMI(past;next|present) = {S:.10f}")
    dP, dQ = check_commuting_square()
    print(f"C3  commuting-square defect: generic law {dP:.3e} ; its Markov projection {dQ:.3e}")
    err, out, loc_ok = check_recovery()
    print(f"C4  heat-bath on A=[2,3], t=200: ||recovered - mu||_1 = {err:.2e} ; limit depends only on sigma_(1,4): {loc_ok}")
    for t, w, wt in out:
        print(f"    time-average t={t:5.1f}: max Dirichlet form of R_t f over 50 sign vectors = {w:.4f}, times t = {wt:.3f}")
    gap, alo, etas = glauber_gap_vs_bounds()
    print(f"C5  Ising chain n=6: Glauber gap {gap:.4f} ; ALO Thm 1.3 lower bound {alo:.4f} ; eta_i = "
          + ", ".join(f"{e:.3f}" for e in etas))
