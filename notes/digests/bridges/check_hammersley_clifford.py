"""Numerical sanity checks for notes/digests/bridges/hammersley-clifford.md.

Pure numpy (no scipy). Each check corresponds to a labelled statement in the digest:
  H1  Moebius inversion of log P at a vacuum: Phi_A vanishes off cliques for a positive
      Gibbs law on the 4-cycle, sums back to log P, and does not vanish for a generic law.
  H2  Approximate HC (digest Sec. 5.3, Mine): |Phi_A| <= 2^{|A|-2} max |log cross-ratio|.
  H3  Moussouris' support: every law on it is pairwise Markov on C4, every edge pattern
      occurs in it (so no edge factor can vanish), its flip graph is an 8-cycle with no
      commuting square; a strictly positive Markov extension exists (Gandolfi-Lenarda
      Ex. 6.2) for the uniform law but not for P* (Gandolfi-Lenarda Lemma 5.2):
      the cross-ratio constraints have non-zero holonomy.
  H4  Brown-Poulin wheel: Gibbs state of a non-commuting clique Hamiltonian on the
      5-vertex wheel is an exact quantum Markov network.
  H5  Leifer-Poulin Ex. 4.8: Heisenberg 3-chain Gibbs state has log rho 2-local (cumulants
      on {A,C} and {A,B,C} vanish) yet I(A:C|B) > 0; the diagonal (Ising) chain has I = 0.
  H6  Lauritzen-Zwiernik two-clique criterion: Tr T(R) <= 1 with
      T = exp(log rho_AB + log rho_BC - log rho_B), equality for a Markov state, and
      Tr rho(log rho - log T) = I(A:C|B) (their Lemma 3.2 with Delta_R = 0; their D for an
      unnormalised second argument carries the extra -1 + Tr T).
  H7  Junge-Renner-Sutter-Wilde-Winter Cor. 4.1: I(A:C|B) >= -2 log F(rho, R(rho_AB))
      for the beta_0-averaged rotated Petz map (log base 2).
  H8  Quantitative intersection (digest Sec. 13, Mine): for h(b,c,d) = P(A=1|b,c,d),
      ||h - E[h|C]|| <= (||h - E[h|C,D]|| + ||h - E[h|B,C]||) / (1 - c_F), with c_F the cosine
      of the Friedrichs angle between L2(C,D) and L2(B,C) in L2(P_BCD) (the norm of the
      two-block Gibbs sampler off the functions of C); c_F = 1 when the support is not
      axis-connected (X_A = X_B = X_D), where intersection fails.
Run:  python3 check_hammersley_clifford.py
"""
import itertools
import numpy as np

rng = np.random.default_rng(1)
LOG2 = np.log(2.0)

# ------------------------------------------------------------------ helpers (classical)
EDGES_C4 = {(0, 1), (1, 2), (2, 3), (0, 3)}


def adjacent(i, j, edges):
    return (min(i, j), max(i, j)) in edges


def is_clique(A, edges):
    return all(adjacent(i, j, edges) for i, j in itertools.combinations(A, 2))


def moebius(g, n, q, vac):
    """Phi_A(x) = sum_{B subset A} (-1)^{|A\\B|} g(x^B), x^B = x on B, vacuum off B."""
    confs = list(itertools.product(range(q), repeat=n))
    phi = {}
    for r in range(n + 1):
        for A in itertools.combinations(range(n), r):
            vals = {}
            for x in confs:
                s = 0.0
                for rb in range(len(A) + 1):
                    for B in itertools.combinations(A, rb):
                        xb = tuple(x[i] if i in B else vac[i] for i in range(n))
                        s += (-1) ** (len(A) - len(B)) * g[xb]
                vals[x] = s
            phi[A] = vals
    return phi, confs


def check_H1_H2():
    n, q = 4, 3
    vac = (0,) * n
    confs = list(itertools.product(range(q), repeat=n))
    # positive Gibbs law with random edge and vertex potentials on C4
    pot_e = {e: rng.normal(size=(q, q)) for e in EDGES_C4}
    pot_v = [rng.normal(size=q) for _ in range(n)]
    logp = {x: sum(pot_e[e][x[e[0]], x[e[1]]] for e in EDGES_C4)
            + sum(pot_v[i][x[i]] for i in range(n)) for x in confs}
    Z = np.log(sum(np.exp(v) for v in logp.values()))
    g = {x: v - Z for x, v in logp.items()}
    phi, _ = moebius(g, n, q, vac)
    off = max(abs(v) for A, d in phi.items() if not is_clique(A, EDGES_C4) for v in d.values())
    recon = max(abs(sum(phi[A][x] for A in phi) - g[x]) for x in confs)
    # generic positive law (not Markov on C4)
    w = rng.random(len(confs)) + 0.1
    g2 = {x: np.log(v) for x, v in zip(confs, w / w.sum())}
    phi2, _ = moebius(g2, n, q, vac)
    off2 = max(abs(v) for A, d in phi2.items() if not is_clique(A, EDGES_C4) for v in d.values())
    print(f"H1 Gibbs on C4: max |Phi_A| off cliques = {off:.2e}; reconstruction error = {recon:.2e}")
    print(f"H1 generic law: max |Phi_A| off cliques = {off2:.3f} (non-zero, as expected)")
    # H2: approximate HC bound, with a perturbed (approximately Markov) law
    eps = 0.05
    g3 = {x: g[x] + eps * rng.normal() for x in confs}
    phi3, _ = moebius(g3, n, q, vac)
    worst_ratio = 0.0
    for A, d in phi3.items():
        if is_clique(A, EDGES_C4):
            continue
        i, j = next((i, j) for i, j in itertools.combinations(A, 2) if not adjacent(i, j, EDGES_C4))
        rest = [k for k in A if k not in (i, j)]
        for x in confs:
            m = 0.0
            for rb in range(len(rest) + 1):
                for B in itertools.combinations(rest, rb):
                    def xs(S):
                        return tuple(x[k] if k in S else vac[k] for k in range(n))
                    L = g3[xs(B + (i, j))] - g3[xs(B + (i,))] - g3[xs(B + (j,))] + g3[xs(B)]
                    m = max(m, abs(L))
            bound = 2 ** (len(A) - 2) * m
            if bound > 0:
                worst_ratio = max(worst_ratio, abs(d[x]) / bound)
    print(f"H2 approximate HC: max |Phi_A| / (2^(|A|-2) max|log cross-ratio|) = {worst_ratio:.3f} (<= 1)")


# ------------------------------------------------------------------ H3 Moussouris
OMEGA0 = [(0, 0, 0, 0), (1, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 0),
          (1, 1, 1, 1), (0, 1, 1, 1), (0, 0, 1, 1), (0, 0, 0, 1)]
OMEGA00 = [(1, 1, 1, 1), (0, 1, 1, 1), (1, 0, 1, 1), (0, 0, 1, 1),
           (1, 1, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 0, 0)]


def cmi_classical(P, A, C, B):
    """I(X_A : X_C | X_B) in nats for a dict law on {0,1}^4."""
    def marg(S):
        m = {}
        for x, p in P.items():
            k = tuple(x[i] for i in S)
            m[k] = m.get(k, 0.0) + p
        return m
    pAB, pBC, pB, pABC = marg(A + B), marg(B + C), marg(B), marg(A + B + C)
    I = 0.0
    for x, p in P.items():
        if p <= 0:
            continue
        a = tuple(x[i] for i in A); b = tuple(x[i] for i in B); c = tuple(x[i] for i in C)
        I += p * np.log(pABC[a + b + c] * pB[b] / (pAB[a + b] * pBC[b + c]))
    return I


def cross_ratio_rows():
    """The 8 relations (5) of Gandolfi-Lenarda: pairs (1,3)|(2,4) and (2,4)|(1,3)."""
    rows = []
    for (x, y) in [(0, 2), (1, 3)]:
        others = [k for k in range(4) if k not in (x, y)]
        for vals in itertools.product([0, 1], repeat=2):
            def c(vx, vy):
                z = [None] * 4
                z[x], z[y] = vx, vy
                z[others[0]], z[others[1]] = vals
                return tuple(z)
            rows.append((c(1, 1), c(0, 0), c(1, 0), c(0, 1)))   # p11 p00 = p10 p01
    return rows


def positive_markov_extension_residual(P):
    """Least-squares residual of: log-ratios fixed on supp P, plus the 8 cross-ratio equations."""
    confs = list(itertools.product([0, 1], repeat=4))
    idx = {x: k for k, x in enumerate(confs)}
    supp = [x for x in confs if P.get(x, 0) > 0]
    rows, rhs = [], []
    x0 = supp[0]
    for x in supp[1:]:
        r = np.zeros(16); r[idx[x]] = 1; r[idx[x0]] = -1
        rows.append(r); rhs.append(np.log(P[x]) - np.log(P[x0]))
    for a, b, c, d in cross_ratio_rows():
        r = np.zeros(16); r[idx[a]] += 1; r[idx[b]] += 1; r[idx[c]] -= 1; r[idx[d]] -= 1
        rows.append(r); rhs.append(0.0)
    M, v = np.array(rows), np.array(rhs)
    sol, *_ = np.linalg.lstsq(M, v, rcond=None)
    return np.linalg.norm(M @ sol - v)


def check_H3():
    # every law on Omega0 is pairwise Markov on C4
    worst = 0.0
    for _ in range(20):
        w = rng.random(8) + 0.05
        P = dict(zip(OMEGA0, w / w.sum()))
        worst = max(worst, abs(cmi_classical(P, [0], [2], [1, 3])), abs(cmi_classical(P, [1], [3], [0, 2])))
    patterns = {e: {(x[e[0]], x[e[1]]) for x in OMEGA0} for e in EDGES_C4}
    all_patterns = all(len(s) == 4 for s in patterns.values())
    # flip graph of the support
    S = set(OMEGA0)
    deg = {x: sum(1 for i in range(4) if tuple(x[k] ^ (k == i) for k in range(4)) in S) for x in S}
    squares = 0
    for x in S:
        for i, j in itertools.combinations(range(4), 2):
            xi = tuple(x[k] ^ (k == i) for k in range(4)); xj = tuple(x[k] ^ (k == j) for k in range(4))
            xij = tuple(x[k] ^ (k in (i, j)) for k in range(4))
            if xi in S and xj in S and xij in S:
                squares += 1
    print(f"H3 Moussouris support: max CMI over random laws on it = {worst:.1e}; "
          f"every edge pattern occurs: {all_patterns}; flip-graph degrees = {sorted(set(deg.values()))}; "
          f"commuting squares inside = {squares}")
    Pm = {x: 1 / 8 for x in OMEGA0}
    Pstar = {x: (2 / 9 if x == (1, 1, 1, 1) else 1 / 9) for x in OMEGA00}
    cm = max(abs(cmi_classical(Pstar, [0], [2], [1, 3])), abs(cmi_classical(Pstar, [1], [3], [0, 2])))
    print(f"H3 P* (Gandolfi-Lenarda Lemma 5.2): pairwise CMIs = {cm:.1e}")
    print(f"H3 positive Markov extension residual: Moussouris uniform = "
          f"{positive_markov_extension_residual(Pm):.2e} (extension exists), "
          f"P* = {positive_markov_extension_residual(Pstar):.3f} (holonomy: none exists)")


# ------------------------------------------------------------------ helpers (quantum)
I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zp = np.array([[1, 0], [0, -1]], dtype=complex)


def op(n, d):
    """Tensor product on n qubits; d maps site -> 2x2 matrix."""
    out = np.array([[1.0 + 0j]])
    for k in range(n):
        out = np.kron(out, d.get(k, I2))
    return out


def herm_fun(M, f):
    w, V = np.linalg.eigh((M + M.conj().T) / 2)
    return (V * f(w)) @ V.conj().T


def ptrace(rho, keep, n):
    rho = rho.reshape([2] * (2 * n))
    tr = [k for k in range(n) if k not in keep]
    for k in sorted(tr, reverse=True):
        nn = rho.ndim // 2
        rho = np.trace(rho, axis1=k, axis2=k + nn)
    m = 2 ** len(keep)
    return rho.reshape(m, m)


def S(rho):
    w = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    w = w[w > 1e-15]
    return float(-(w * np.log(w)).sum())


def qcmi(rho, A, C, B, n):
    return (S(ptrace(rho, sorted(A + B), n)) + S(ptrace(rho, sorted(B + C), n))
            - S(ptrace(rho, sorted(B), n)) - S(ptrace(rho, sorted(A + B + C), n)))


def gibbs(H):
    r = herm_fun(H, lambda w: np.exp(w - w.max()))
    return r / np.trace(r).real


def check_H4():
    n = 5   # sites 0..3 = rim 1..4, site 4 = centre 5
    h5 = op(n, {0: Zp, 1: Zp, 4: Y}); hC = op(n, {1: Zp, 2: Zp, 4: X})
    h4 = op(n, {2: Zp, 3: Zp, 4: Y}); hB = op(n, {3: Zp, 0: Zp, 4: X})
    worst = 0.0
    for beta in [0.3, 1.0, 2.5]:
        rho = gibbs(beta * (h5 + hC + h4 + hB))
        worst = max(worst, abs(qcmi(rho, [0], [2], [1, 3, 4], n)), abs(qcmi(rho, [1], [3], [0, 2, 4], n)))
    comm = np.linalg.norm(h5 @ hC - hC @ h5)
    # generic non-commuting clique Hamiltonian on the same wheel: Markov fails
    paulis = [X, Y, Zp]
    Hgen = np.zeros((32, 32), dtype=complex)
    for tri in [(0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)]:
        for a, b, c in itertools.product(range(3), repeat=3):
            Hgen += rng.normal() * 0.3 * op(n, {tri[0]: paulis[a], tri[1]: paulis[b], tri[2]: paulis[c]})
    rg = gibbs(Hgen)
    gen = qcmi(rg, [0], [2], [1, 3, 4], n)
    print(f"H4 Brown-Poulin wheel: max CMI over shields and beta = {worst:.1e}; "
          f"||[h5,hC]|| = {comm:.1f} (terms do not commute); generic triangle H: I(1:3|245) = {gen:.2e}")


def cumulant(Hm, Xset, n):
    """Brown-Poulin cumulant K_X = prod_{i in X}(id - E_i) prod_{i notin X} E_i (H),
    E_i = normalised partial trace on site i, tensored back with the identity."""
    def E(M, i):
        keep = [k for k in range(n) if k != i]
        return lift(ptrace(M, keep, n) / 2, keep, n)
    M = Hm.copy()
    for i in range(n):
        M = (M - E(M, i)) if i in Xset else E(M, i)
    return M


def check_H5():
    n = 3
    heis = lambda a, b: sum(op(n, {a: P, b: P}) for P in (X, Y, Zp))
    H = heis(0, 1) + heis(1, 2)
    Hz = op(n, {0: Zp, 1: Zp}) + op(n, {1: Zp, 2: Zp})
    res = []
    for beta in [0.5, 1.0, 2.0, 4.0]:
        rho = gibbs(-beta * H)
        res.append(qcmi(rho, [0], [2], [1], n) / LOG2)
    rho = gibbs(-1.0 * H)
    L = herm_fun(rho, np.log)
    k02 = np.linalg.norm(cumulant(L, {0, 2}, n)); k012 = np.linalg.norm(cumulant(L, {0, 1, 2}, n))
    k01 = np.linalg.norm(cumulant(L, {0, 1}, n))
    ising = qcmi(gibbs(-1.0 * Hz), [0], [2], [1], n)
    print("H5 Heisenberg 3-chain: I(A:C|B) in bits at beta=0.5,1,2,4: "
          + ", ".join(f"{v:.4f}" for v in res))
    print(f"H5 cumulants of log rho: |K_AC| = {k02:.1e}, |K_ABC| = {k012:.1e}, |K_AB| = {k01:.2f}; "
          f"Ising chain I = {ising:.1e}")


def rel_ent(r, s):
    return float(np.trace(r @ (herm_fun(r, lambda w: np.log(np.clip(w, 1e-300, None)))
                               - herm_fun(s, np.log))).real)


def lift(M, keep, n):
    """Embed an operator on sites `keep` (sorted) into n qubits, identity elsewhere."""
    k = len(keep)
    perm_in = keep + [i for i in range(n) if i not in keep]
    big = np.kron(M, np.eye(2 ** (n - k)))
    big = big.reshape([2] * (2 * n))
    inv = np.argsort(perm_in)
    big = np.transpose(big, list(inv) + [n + i for i in inv])
    return big.reshape(2 ** n, 2 ** n)


def check_H6():
    n = 3
    worst_tr, worst_id = -np.inf, 0.0
    for _ in range(10):
        G = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
        rho = G @ G.conj().T; rho /= np.trace(rho).real
        rAB, rBC, rB = ptrace(rho, [0, 1], n), ptrace(rho, [1, 2], n), ptrace(rho, [1], n)
        T = herm_fun(lift(herm_fun(rAB, np.log), [0, 1], n) + lift(herm_fun(rBC, np.log), [1, 2], n)
                     - lift(herm_fun(rB, np.log), [1], n), np.exp)
        trT = np.trace(T).real
        worst_tr = max(worst_tr, trT - 1)
        lhs = rel_ent(rho, T)
        worst_id = max(worst_id, abs(lhs - qcmi(rho, [0], [2], [1], n)))
    # Markov state: commuting Gibbs on the chain
    rho = gibbs(-0.8 * (op(n, {0: X, 1: Zp}) + op(n, {1: Zp, 2: Y})))
    rAB, rBC, rB = ptrace(rho, [0, 1], n), ptrace(rho, [1, 2], n), ptrace(rho, [1], n)
    T = herm_fun(lift(herm_fun(rAB, np.log), [0, 1], n) + lift(herm_fun(rBC, np.log), [1, 2], n)
                 - lift(herm_fun(rB, np.log), [1], n), np.exp)
    print(f"H6 Lauritzen-Zwiernik: max(Tr T - 1) over random states = {worst_tr:.2e} (<= 0); "
          f"identity Tr rho(log rho - log T) = I, error = {worst_id:.1e}; Markov state: Tr T = {np.trace(T).real:.10f}, "
          f"||T - rho|| = {np.linalg.norm(T - rho):.1e}")


def fidelity(r, s):
    sr = herm_fun(r, lambda w: np.sqrt(np.clip(w, 0, None)))
    return float(np.sum(np.sqrt(np.clip(np.linalg.eigvalsh(sr @ s @ sr), 0, None))))


def check_H7():
    n = 3
    worst = np.inf
    for _ in range(5):
        G = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
        rho = G @ G.conj().T; rho /= np.trace(rho).real
        rAB, rBC, rB = ptrace(rho, [0, 1], n), ptrace(rho, [1, 2], n), ptrace(rho, [1], n)
        sBC = lift(rBC, [1, 2], n); sB = lift(rB, [1], n)
        x_in = lift(rAB, [0, 1], n)   # rho_AB (x) 1_C as input to the map
        ts = np.linspace(-12, 12, 2401); dt = ts[1] - ts[0]
        beta0 = (np.pi / 2) / (np.cosh(np.pi * ts) + 1)
        out = np.zeros((8, 8), dtype=complex)
        sBC_h = herm_fun(sBC, np.sqrt); sB_mh = herm_fun(sB, lambda w: w ** -0.5)
        for t, wgt in zip(ts, beta0):
            s = t / 2
            U = herm_fun(sB, lambda w: w ** (1j * s)); V = herm_fun(sBC, lambda w: w ** (-1j * s))
            inner = sB_mh @ U @ x_in @ U.conj().T @ sB_mh
            # partial trace over C then re-tensor with 1_C is already implicit: x_in = rho_AB (x) 1
            # but the Petz map uses rho_B^{-1/2} X_AB rho_B^{-1/2} (x) 1_C: x_in already has 1_C
            rec = V @ sBC_h @ inner @ sBC_h @ V.conj().T
            out += wgt * rec * dt
        tr_out = np.trace(out).real
        out = out / tr_out
        I = qcmi(rho, [0], [2], [1], n) / LOG2
        bound = -2 * np.log2(fidelity(rho, out))
        worst = min(worst, I - bound)
        print(f"   H7 trial: Tr R(rho_AB) before normalisation = {tr_out:.4f}, I = {I:.4f} bits, -2log2F = {bound:.4f}")
    print(f"H7 JRSWW universal recovery: min over random states of I(A:C|B) - (-2 log2 F) = {worst:.3f} (>= 0)")


def block_sampler_data(P):
    """P: array over (a,b,c,d) binary. Returns h, weights and the operators on functions of (b,c,d)."""
    pts = list(itertools.product([0, 1], repeat=3))           # (b,c,d)
    w = np.array([P[:, b, c, d].sum() for b, c, d in pts])
    h = np.array([P[1, b, c, d] / P[:, b, c, d].sum() if P[:, b, c, d].sum() > 0 else 0.0 for b, c, d in pts])
    k = {pt: i for i, pt in enumerate(pts)}
    def cond(keep):
        M = np.zeros((8, 8))
        for (b, c, d) in pts:
            grp = [q for q in pts if all(q[t] == (b, c, d)[t] for t in keep)]
            tot = sum(w[k[q]] for q in grp)
            for q in grp:
                M[k[(b, c, d)], k[q]] = w[k[q]] / tot if tot > 0 else 0.0
        return M
    E_notB = cond([1, 2])   # E[. | C, D]
    E_notD = cond([0, 1])   # E[. | B, C]
    E_C = cond([1])
    return h, w, E_notB, E_notD, E_C


def wnorm(f, w):
    return float(np.sqrt(np.sum(w * f * f)))


def check_H8():
    worst = 0.0
    for _ in range(200):
        P = rng.random((2, 2, 2, 2)) ** 2 + 1e-3
        P /= P.sum()
        h, w, EB, ED, EC = block_sampler_data(P)
        T = EB @ ED
        W, Wi = np.diag(np.sqrt(w)), np.diag(1 / np.sqrt(w))
        cF = np.linalg.norm(W @ T @ (np.eye(8) - EC) @ Wi, 2)
        lhs = wnorm(h - EC @ h, w)
        rhs = (wnorm(h - EB @ h, w) + wnorm(h - ED @ h, w)) / (1 - cF)
        worst = max(worst, lhs / rhs)
    # degenerate law: X_A = X_B = X_D uniform, C constant 0
    P = np.zeros((2, 2, 2, 2)); P[0, 0, 0, 0] = P[1, 1, 0, 1] = 0.5
    h, w, EB, ED, EC = block_sampler_data(P)
    m = w > 0
    T = EB @ ED
    # restrict to the support of P_BCD
    idx = np.where(m)[0]
    Ts, ECs, ws = T[np.ix_(idx, idx)], EC[np.ix_(idx, idx)], w[idx]
    Ws, Wsi = np.diag(np.sqrt(ws)), np.diag(1 / np.sqrt(ws))
    cF_deg = np.linalg.norm(Ws @ Ts @ (np.eye(len(idx)) - ECs) @ Wsi, 2)
    hs = h[idx]
    viol = wnorm(hs - ECs @ hs, ws)
    devB = wnorm(hs - EB[np.ix_(idx, idx)] @ hs, ws); devD = wnorm(hs - ED[np.ix_(idx, idx)] @ hs, ws)
    print(f"H8 quantitative intersection: max lhs/rhs over 200 positive laws = {worst:.3f} (<= 1); "
          f"degenerate law X_A=X_B=X_D: c_F = {cF_deg:.3f}, ||h-E[h|B,D-free]|| = {devB:.1e}, {devD:.1e}, "
          f"but ||h - E[h|C]|| = {viol:.3f}")


if __name__ == "__main__":
    check_H1_H2()
    check_H3()
    check_H4()
    check_H5()
    check_H6()
    check_H7()
    check_H8()
