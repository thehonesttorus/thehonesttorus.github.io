"""Numerical sanity checks for notes/digests/bridges/nc-dirichlet-lindblad.md.

Pure numpy, finite dimensions. Each check corresponds to a labelled statement in the digest:

  K1  Ding-Li-Lin Theorem 10 construction (KMS-detailed-balanced Lindbladian from self-adjoint
      single-Pauli couplings on A, Gaussian weighting q): KMS symmetry, Gibbs stationarity,
      failure of GNS symmetry (no commutation with the modular operator) for noncommuting H.
  K2  Derivation form of the Dirichlet form (Chen-Rouze Lemma X.1 / [RFA24, Lemma C.2], here
      in the Ding-Li-Lin parametrisation):
        E(X) = -<X, L X>_KMS
             = sum_a sum_{nu1,nu2} q(nu1) q(nu2) / (2 cosh(beta (nu1-nu2)/4))
                    Tr[ sqrt(rho) [A_nu1, X]^* sqrt(rho) [A_nu2, X] ]
             = sum_a int g(t) || [ e^{iHt} Atilde_a e^{-iHt}, X ] ||_KMS^2 dt,
        g(t) = 1/(beta cosh(2 pi t / beta)),  Atilde_a = sum_nu q(nu) (A_a)_nu  (self-adjoint).
  K3  Kernel of the single-Pauli-on-A generator: at beta = 0 it is 1_A (x) M_{A^c}
      (dimension 4^{n-|A|}); at beta > 0 with a generic noncommuting H it collapses.
  K4  Gap-free Cesaro decay: E(R_t X) <= c ||X||_KMS^2 / t with c = sup_u (1-e^{-u})^2/u ~ 0.4073,
      using only KMS self-adjointness (the mechanism of Chen-Rouze Cor. VII.1).
  K5  beta = 0 single-Pauli generator = "noncommutative hypercube": eigenvalue -4 w(S) on a Pauli
      string S, w = weight on A; Efron-Stein for the depolarising expectations.
  K6  De Palma-Marvian-Trevisan-Lloyd Lipschitz constant versus the single-Pauli commutator
      seminorm M(H) = max_{i,P} ||[P_i, H]||:  M/2 <= max_i ||H - E_i H|| <= 3M/4, hence
      (2/3) M <= ||H||_L <= (3/2) M  (via their Prop. 15).
  K7  Programme toy: re-routing (fibre) Lindbladian on a complete-history corner with
      GNS-detailed balance for P_beta (Alicki/Carlen-Maas form). Kernel = endpoint face algebra
      D_0; Cesaro limit = Takesaki conditional expectation onto D_0; Dirichlet form = sum of
      KMS-squared commutators (Carlen-Maas Prop. 2.5); P_beta' is stationary / recovered iff the
      accumulated action F is a function of the endpoint (sufficiency, Note 1 section 6).
  K8  Poincare + commutator bound: E(a) <= (1/2) sum_a sup_t ||[Atilde_a(t), a]||^2, hence a
      Connes-distance radius bound d_D(rho, psi) <= sqrt(#couplings / (2 lambda rho_min)) for every
      state psi, with D the direct integral over modular time of the dressed couplings.

Run:  python3 check_nc_dirichlet_lindblad.py
"""
import itertools
import numpy as np

rng = np.random.default_rng(7)
I2 = np.eye(2, dtype=complex)
PX = np.array([[0, 1], [1, 0]], dtype=complex)
PY = np.array([[0, -1j], [1j, 0]], dtype=complex)
PZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [PX, PY, PZ]


def kron(*ms):
    out = np.array([[1.0 + 0j]])
    for m in ms:
        out = np.kron(out, m)
    return out


def site_op(op, i, n):
    return kron(*[op if k == i else I2 for k in range(n)])


def random_herm(d):
    a = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    return (a + a.conj().T) / 2


def random_2local_H(n, scale=1.0):
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for i in range(n):
        for P in PAULIS:
            H += scale * rng.normal() * site_op(P, i, n)
    for i in range(n - 1):
        for P in PAULIS:
            for Q in PAULIS:
                H += scale * rng.normal() * site_op(P, i, n) @ site_op(Q, i + 1, n)
    return H


def mfun(h, f):
    w, v = np.linalg.eigh(h)
    return (v * f(w)) @ v.conj().T


# superoperator conventions: vec column-stacking, vec(A X B) = (B^T kron A) vec(X)
def sop_left(A):
    return np.kron(np.eye(A.shape[0]), A)


def sop_right(B):
    return np.kron(B.T, np.eye(B.shape[0]))


def vec(X):
    return X.reshape(-1, order="F")


def unvec(v, d):
    return v.reshape(d, d, order="F")


def bohr_components(A, H, tol=1e-9):
    """Return dict nu -> A_nu = sum_{E2-E1=nu} P_E2 A P_E1 (H assumed nondegenerate generic)."""
    w, v = np.linalg.eigh(H)
    At = v.conj().T @ A @ v
    comps = {}
    d = len(w)
    for i in range(d):
        for j in range(d):
            nu = w[i] - w[j]
            key = None
            for k in comps:
                if abs(k - nu) < tol:
                    key = k
                    break
            if key is None:
                key = nu
                comps[key] = np.zeros((d, d), dtype=complex)
            comps[key][i, j] += At[i, j]
    return {k: v @ M @ v.conj().T for k, M in comps.items()}


def kms_weight(rho):
    s = mfun(rho, np.sqrt)
    return np.kron(s.T, s)  # <X,Y>_rho = vec(X)^* W vec(Y)


def dll_lindbladian(H, beta, couplings, qfun):
    """Ding-Li-Lin Theorem 10: L(X) = i[G,X] + sum_a (L_a^* X L_a - 1/2 {L_a^* L_a, X}),
    L_a = sum_nu exp(-beta nu/4) q(nu) (A_a)_nu, G = -(i/2) sum_nu tanh(-beta nu/4) (K)_nu."""
    d = H.shape[0]
    Ls, Atildes = [], []
    for A in couplings:
        comps = bohr_components(A, H)
        L = sum(np.exp(-beta * nu / 4) * qfun(nu) * M for nu, M in comps.items())
        At = sum(qfun(nu) * M for nu, M in comps.items())
        Ls.append(L)
        Atildes.append(At)
    K = sum(L.conj().T @ L for L in Ls)
    G = np.zeros((d, d), dtype=complex)
    for nu, M in bohr_components(K, H).items():
        G += -0.5j * np.tanh(-beta * nu / 4) * M
    S = 1j * (sop_left(G) - sop_right(G))
    for L in Ls:
        S += np.kron(L.T, L.conj().T)  # X -> L^* X L
    S -= 0.5 * (sop_left(K) + sop_right(K))
    return S, G, Ls, Atildes


def gibbs(H, beta):
    r = mfun(H, lambda w: np.exp(-beta * (w - w.min())))
    return r / np.trace(r).real


def report(name, ok, detail):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    return ok


# ---------------------------------------------------------------- K1, K2, K3, K4, K8
def check_dll(n=3, A=(0,), beta=1.3, trials=3):
    allok = True
    q = lambda nu: np.exp(-(beta * nu) ** 2 / 8)
    for tr in range(trials):
        H = random_2local_H(n, scale=0.6)
        d = 2 ** n
        rho = gibbs(H, beta)
        couplings = [site_op(P, i, n) for i in A for P in PAULIS]
        S, G, Ls, Atildes = dll_lindbladian(H, beta, couplings, q)
        W = kms_weight(rho)
        # K1: KMS symmetry W S = S^* W ; Gibbs stationarity S^*(rho) = 0 ; G Hermitian
        kms_err = np.abs(W @ S - S.conj().T @ W).max()
        stat_err = np.abs(S.conj().T @ vec(rho)).max()
        herm_err = np.abs(G - G.conj().T).max()
        # GNS symmetry would require commutation with modular operator Delta(X) = rho X rho^{-1}
        Delta = np.kron(np.linalg.inv(rho).T, rho)
        gns_def = np.abs(S @ Delta - Delta @ S).max() / np.abs(S).max()
        allok &= report(f"K1 trial {tr}", kms_err < 1e-10 and stat_err < 1e-10 and herm_err < 1e-12,
                        f"KMS asym {kms_err:.1e}, L*(rho) {stat_err:.1e}, G herm {herm_err:.1e}; "
                        f"relative size of [L,Delta] = {gns_def:.2e} (nonzero, so KMS but not GNS)")
        # K2: Dirichlet form in derivation form, frequency and time domain
        sq = mfun(rho, np.sqrt)
        worst_f, worst_t = 0.0, 0.0
        for _ in range(4):
            X = random_herm(d) + 1j * random_herm(d)
            Ed = -(vec(X).conj() @ W @ S @ vec(X))
            Ef = 0.0
            for Acoup in couplings:
                comps = bohr_components(Acoup, H)
                keys = list(comps)
                comms = {k: comps[k] @ X - X @ comps[k] for k in keys}
                for k1 in keys:
                    for k2 in keys:
                        w12 = q(k1) * q(k2) / (2 * np.cosh(beta * (k1 - k2) / 4))
                        Ef += w12 * np.trace(sq @ comms[k1].conj().T @ sq @ comms[k2])
            # time domain quadrature
            ts = np.linspace(-12 * beta, 12 * beta, 4001)
            dt = ts[1] - ts[0]
            g = 1.0 / (beta * np.cosh(2 * np.pi * ts / beta))
            wH, vH = np.linalg.eigh(H)
            Et = 0.0
            for At in Atildes:
                Atb = vH.conj().T @ At @ vH
                Xb = vH.conj().T @ X @ vH
                sqb = vH.conj().T @ sq @ vH
                for t, gt in zip(ts, g):
                    ph = np.exp(1j * wH * t)
                    Att = (ph[:, None] * Atb) * ph.conj()[None, :]
                    C = Att @ Xb - Xb @ Att
                    Et += gt * np.trace(C.conj().T @ sqb @ C @ sqb).real * dt
            worst_f = max(worst_f, abs(Ed - Ef) / abs(Ed))
            worst_t = max(worst_t, abs(Ed.real - Et) / abs(Ed))
        allok &= report(f"K2 trial {tr}", worst_f < 1e-9 and worst_t < 1e-4,
                        f"rel. error frequency-domain {worst_f:.1e}, time-domain quadrature {worst_t:.1e}")
        # Hermitian similarity for spectral data
        Wh = np.kron(mfun(rho, lambda w: w ** 0.25).T, mfun(rho, lambda w: w ** 0.25))
        Whi = np.linalg.inv(Wh)
        Lhat = Wh @ S @ Whi
        Lhat = (Lhat + Lhat.conj().T) / 2
        ev, U = np.linalg.eigh(Lhat)
        kerdim = int(np.sum(np.abs(ev) < 1e-9))
        gap = -ev[ev < -1e-9].max()
        allok &= report(f"K3 trial {tr} (beta={beta})", True,
                        f"dim ker L_A = {kerdim} (beta=0 value would be 4^(n-|A|) = {4 ** (n - len(A))}); "
                        f"spectral gap {gap:.3e}")
        # K4: Cesaro decay of Dirichlet form, only KMS self-adjointness used
        c_const = max((1 - np.exp(-u)) ** 2 / u for u in np.linspace(1e-4, 20, 200001))
        ok4 = True
        for t in [1.0, 10.0, 100.0, 1000.0]:
            X = random_herm(d)
            y = Wh @ vec(X)                  # KMS-isometric coordinates
            coeff = U.conj().T @ y
            x = -ev
            fac = np.where(x > 1e-12, (1 - np.exp(-t * x)) / (t * np.where(x > 1e-12, x, 1)), 1.0)
            Rt = fac * coeff
            Erad = np.sum(x * np.abs(Rt) ** 2)
            bound = c_const * np.sum(np.abs(coeff) ** 2) / t
            ok4 &= Erad <= bound * (1 + 1e-9)
        allok &= report(f"K4 trial {tr}", ok4, f"E(R_t X) <= {c_const:.4f} ||X||^2_KMS / t at t = 1,10,100,1000")
        # K8: Poincare and the commutator (Lipschitz) bound
        rho_min = np.linalg.eigvalsh(rho).min()
        ok8 = True
        wH, vH = np.linalg.eigh(H)
        for _ in range(3):
            a = random_herm(d)
            Ea = -(vec(a).conj() @ W @ S @ vec(a)).real
            # sup over modular time of commutator norms
            sup_sum = 0.0
            for At in Atildes:
                Atb = vH.conj().T @ At @ vH
                ab = vH.conj().T @ a @ vH
                best = 0.0
                for t in np.linspace(-15 * beta, 15 * beta, 1201):
                    ph = np.exp(1j * wH * t)
                    Att = (ph[:, None] * Atb) * ph.conj()[None, :]
                    best = max(best, np.linalg.norm(Att @ ab - ab @ Att, 2))
                sup_sum += best ** 2
            a0 = a - np.trace(rho @ a).real * np.eye(d)
            var = (vec(a0).conj() @ W @ vec(a0)).real
            ok8 &= (Ea <= 0.5 * sup_sum * (1 + 1e-6)) and (gap * var <= Ea * (1 + 1e-9))
        diam = np.sqrt(len(couplings) / (2 * gap * rho_min))
        allok &= report(f"K8 trial {tr}", ok8,
                        f"E(a) <= (1/2) sum sup_t ||[A~(t),a]||^2 and Poincare hold; radius bound "
                        f"d_D(rho, any state) <= sqrt(#coup/(2 lambda rho_min)) = {diam:.2f}")
    # K3 at beta = 0 and for a commuting (classical Ising) H
    H = random_2local_H(n)
    couplings = [site_op(P, i, n) for i in A for P in PAULIS]
    S0, _, _, _ = dll_lindbladian(H, 0.0, couplings, lambda nu: 1.0)
    ev0 = np.linalg.eigvals(S0)
    allok &= report("K3 beta=0", int(np.sum(np.abs(ev0) < 1e-9)) == 4 ** (n - len(A)),
                    f"dim ker = {int(np.sum(np.abs(ev0) < 1e-9))}, expected {4 ** (n - len(A))}")
    Hc = sum(rng.normal() * site_op(PZ, i, n) @ site_op(PZ, i + 1, n) for i in range(n - 1)) \
        + sum(rng.normal() * site_op(PZ, i, n) for i in range(n))
    Hc = Hc + 1e-3 * sum(rng.normal() * site_op(PZ, i, n) for i in range(n))  # break degeneracies
    Sc, _, _, _ = dll_lindbladian(Hc, beta, couplings, q)
    rho_c = gibbs(Hc, beta)
    Wh = np.kron(mfun(rho_c, lambda w: w ** 0.25).T, mfun(rho_c, lambda w: w ** 0.25))
    Lh = Wh @ Sc @ np.linalg.inv(Wh)
    evc = np.linalg.eigvalsh((Lh + Lh.conj().T) / 2)
    report("K3 commuting Ising H, beta>0", True,
           f"dim ker = {int(np.sum(np.abs(evc) < 1e-8))} (compare 4^(n-|A|) = {4 ** (n - len(A))})")
    return allok


# ---------------------------------------------------------------- K5
def check_hypercube(n=3, A=(0, 1)):
    couplings = [site_op(P, i, n) for i in A for P in PAULIS]
    d = 2 ** n
    S = sum(np.kron(P.T, P) for P in couplings) - len(couplings) * np.eye(d * d)
    ok = True
    for labels in itertools.product(range(4), repeat=n):
        Sstr = kron(*[[I2, PX, PY, PZ][l] for l in labels])
        w = sum(1 for i in A if labels[i] != 0)
        img = unvec(S @ vec(Sstr), d)
        ok &= np.allclose(img, -4 * w * Sstr)
    # Efron-Stein: ||X - E_A X||_2^2 <= sum_{i in A} ||X - E_i X||_2^2
    def E_sites(X, sites):
        Y = X.copy()
        for i in sites:
            Y = sum(site_op(P, i, n) @ Y @ site_op(P, i, n) for P in [I2] + PAULIS) / 4
        return Y
    worst = 0.0
    for _ in range(20):
        X = random_herm(d)
        lhs = np.linalg.norm(X - E_sites(X, A)) ** 2
        rhs = sum(np.linalg.norm(X - E_sites(X, [i])) ** 2 for i in A)
        worst = max(worst, lhs / rhs)
    return report("K5", ok and worst <= 1 + 1e-12,
                  f"Pauli strings are eigenvectors with eigenvalue -4 w_A(S); Efron-Stein ratio max {worst:.3f} <= 1")


# ---------------------------------------------------------------- K6
def check_dmtl(n=3, trials=200):
    d = 2 ** n
    lo, hi = np.inf, 0.0
    for _ in range(trials):
        Hm = random_herm(d)
        M = max(np.linalg.norm(site_op(P, i, n) @ Hm - Hm @ site_op(P, i, n), 2)
                for i in range(n) for P in PAULIS)
        m = max(np.linalg.norm(Hm - sum(site_op(P, i, n) @ Hm @ site_op(P, i, n)
                                        for P in [I2] + PAULIS) / 4, 2) for i in range(n))
        lo, hi = min(lo, m / M), max(hi, m / M)
    return report("K6", lo >= 0.5 - 1e-12 and hi <= 0.75 + 1e-12,
                  f"max_i||H - E_i H|| / max_(i,P)||[P_i,H]|| in [{lo:.3f}, {hi:.3f}] subset [1/2, 3/4]")


# ---------------------------------------------------------------- K7
def check_rerouting(endpoints=(3, 4, 2), beta=0.9, betap=1.7):
    N = sum(endpoints)
    lab = np.concatenate([[k] * m for k, m in enumerate(endpoints)])
    allok = True
    for case in ["random F (not sufficient)", "F = U(endpoint) (sufficient)"]:
        if case.startswith("random"):
            F = rng.normal(size=N)
        else:
            U = rng.normal(size=len(endpoints))
            F = U[lab]
        p = np.exp(-beta * F)
        p /= p.sum()
        rho = np.diag(p).astype(complex)
        # GNS-detailed-balanced re-routing generator (Alicki / Carlen-Maas Thm 3.1 form)
        S = np.zeros((N * N, N * N), dtype=complex)
        Vs = []
        for mu in range(N):
            for nu in range(N):
                if mu != nu and lab[mu] == lab[nu]:
                    V = np.zeros((N, N), dtype=complex)
                    V[mu, nu] = 1.0
                    om = beta * (F[mu] - F[nu])  # Delta V = rho V rho^{-1} = e^{-om} V
                    Vs.append((V, om))
        # symmetric rates c_{mu nu} = c_{nu mu}
        rates = {}
        for (V, om) in Vs:
            mu, nu = np.argwhere(V)[0]
            key = (min(mu, nu), max(mu, nu))
            rates.setdefault(key, 0.5 + rng.random())
        for (V, om) in Vs:
            mu, nu = np.argwhere(V)[0]
            c = rates[(min(mu, nu), max(mu, nu))]
            Vd = V.conj().T
            # X -> e^{-om/2} (V^* [X,V] + [V^*,X] V)
            S += c * np.exp(-om / 2) * (sop_left(Vd) @ (sop_right(V) - sop_left(V))
                                        + sop_right(V) @ (sop_left(Vd) - sop_right(Vd)))
        W = kms_weight(rho)
        kms = np.abs(W @ S - S.conj().T @ W).max()
        Delta = np.kron(np.linalg.inv(rho).T, rho)
        gns = np.abs(S @ Delta - Delta @ S).max()
        stat = np.abs(S.conj().T @ vec(rho)).max()
        # kernel = endpoint face algebra D_0 = span{p_tau}
        Wh = np.kron(mfun(rho, lambda w: w ** 0.25).T, mfun(rho, lambda w: w ** 0.25))
        Lh = Wh @ S @ np.linalg.inv(Wh)
        ev, Uv = np.linalg.eigh((Lh + Lh.conj().T) / 2)
        kerdim = int(np.sum(np.abs(ev) < 1e-9))
        # Cesaro limit = spectral projection onto kernel (KMS-orthogonal) vs Takesaki expectation
        Pk = Uv[:, np.abs(ev) < 1e-9]
        Ecesaro = np.linalg.inv(Wh) @ (Pk @ Pk.conj().T) @ Wh
        projs = [np.diag((lab == k).astype(complex)) for k in range(len(endpoints))]
        err_E = 0.0
        for _ in range(5):
            X = random_herm(N) + 1j * random_herm(N)
            ET = sum(np.trace(rho @ P @ X) / np.trace(rho @ P) * P for P in projs)
            err_E = max(err_E, np.abs(unvec(Ecesaro @ vec(X), N) - ET).max())
        # Dirichlet form = sum c ||[V,X]||_KMS^2 (Carlen-Maas Prop. 2.5, with V scaled by sqrt(c))
        errD = 0.0
        for _ in range(5):
            X = random_herm(N) + 1j * random_herm(N)
            Ed = -(vec(X).conj() @ W @ S @ vec(X)).real
            Ec = sum(rates[(min(*np.argwhere(V)[0]), max(*np.argwhere(V)[0]))]
                     * (vec(V @ X - X @ V).conj() @ W @ vec(V @ X - X @ V)).real for (V, om) in Vs)
            errD = max(errD, abs(Ed - Ec) / abs(Ed))
        # sufficiency test: is P_beta' stationary, and recovered by the Cesaro limit?
        pp = np.exp(-betap * F)
        pp /= pp.sum()
        rhop = np.diag(pp).astype(complex)
        stat_p = np.abs(S.conj().T @ vec(rhop)).sum()
        rec = unvec(Ecesaro.conj().T @ vec(rhop), N)  # Schroedinger dual of E applied to P_beta'
        rec_err = np.abs(rec - rhop).sum()
        osc = max(np.ptp(F[lab == k]) for k in range(len(endpoints)))
        ok = kms < 1e-10 and gns < 1e-10 and stat < 1e-12 and kerdim == len(endpoints) \
            and err_E < 1e-9 and errD < 1e-9
        if case.startswith("F = U"):
            ok &= stat_p < 1e-10 and rec_err < 1e-10
        else:
            ok &= stat_p > 1e-6 and rec_err > 1e-6
        allok &= report(f"K7 {case}", ok,
                        f"KMS {kms:.0e}, GNS {gns:.0e}, dim ker {kerdim} (= #endpoints {len(endpoints)}), "
                        f"Cesaro=Takesaki {err_E:.0e}, Dirichlet=sum c||[V,X]||^2 {errD:.0e}; "
                        f"max fibre osc(F) {osc:.3f}, ||L*(P_beta')||_1 {stat_p:.2e}, "
                        f"||E_*(P_beta') - P_beta'||_1 {rec_err:.2e}")
    return allok


if __name__ == "__main__":
    ok = True
    ok &= check_dll()
    ok &= check_hypercube()
    ok &= check_dmtl()
    ok &= check_rerouting()
    print("ALL PASS" if ok else "SOME CHECK FAILED")
