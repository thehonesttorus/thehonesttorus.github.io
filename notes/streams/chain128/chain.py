"""Dense K=3 moment chain for small-width ReLU MLPs (numpy, float64).

Propagates per layer the pre-activation mean mu, covariance C, the FULL third cumulant kappa3(z) (n^3 entries)
and a fourth-cumulant closure (the slices kappa4(z)_{aabc}, which contain the (4), (3,1), (2,2) and (2,1,1)
slices), through the ReLU and the linear map z_{l+1} = a_l @ W_{l+1} (exact multilinear transport).

ReLU step (all from Gaussian expectations of derivatives of relu^p under the Gaussian with matching mean and
covariance, corrected by the Edgeworth operator  E[F(z)] = E_G[exp(sum_r kappa_r/r! d^r) F]  truncated at
first order in kappa3 and kappa4 plus the kappa3^2 / 2 term):
  1-D   E[a_i^p] = F_p(0) + D3/6 F_p(3) + K4/24 F_p(4) + D3^2/72 F_p(6),  F_p(k) = E_G[(relu^p)^{(k)}(Y_i)] (exact,
        closed form through truncated Gaussian moments; = the paper's Hermite coefficients of powers, Prop S.3.8)
  2-D   E[a_i^p a_j^q] = sum_{(r,s)} c_rs E_G[(relu^p)^{(r)}(z_i) (relu^q)^{(s)}(z_j)], with the c_rs from the
        kappa3 slices (D3, D21) and kappa4 slices (K4, K31, K22) of the pair, and every Gaussian expectation
        evaluated to all orders in C_ij by the Mehler series sum_k C_ij^k / k! F_i(r+k) G_j(s+k)
        -> mean, variance, off-diagonal covariance, the (3,) and (2,1) slices of kappa3(a), and the (4), (3,1),
        (2,2) slices of kappa4(a)
  3-D   the all-distinct part of kappa3(a): Gaussian Hermite (Wick) term + first-order cumulant diagrams, built
        with notes/experiments/oracle_k3.py (hermite_model, residual_basis, CLOSURE_COEF); the variant decides
        which terms and which coefficients:
          'none'     slices only (memoryless)
          'wick'     Gaussian rho^2 + Phi^3 kappa3(z)                              (leading Wick births)
          'closure'  + all first-order diagrams with leg-partition coefficients (oracle closure_model)
          'fit'      same diagrams, per-layer coefficient table (fitted on atlases)
Fourth-cumulant closure for z_{l+1}: kappa4(a) restricted to its (4), (3,1), (2,2) slices (memoryless kappa4),
transported exactly to kappa4(z_{l+1})_{aabc}; options to replace it by an atlas's truth (teacher forcing),
to zero it, or to regenerate its (2,1,1) slice as u_i C_jk.

Usage as a library: Chain(W, **opts).run() -> dict with per-layer post-activation means (L, n) and diagnostics.
"""
import os, sys
import numpy as np
from math import factorial, pi
from scipy.special import ndtr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "experiments"))
import oracle_k3 as ok  # noqa: E402

KMAX = 40   # Mehler series length (|rho| <= 0.6 converges to < 1e-9 relative)


# ----------------------------------------------------------------------------------------------------------
# 1-D Gaussian tables
# ----------------------------------------------------------------------------------------------------------

def herme_vals(x, mmax):
    H = np.empty((mmax + 1,) + np.shape(x))
    H[0] = 1.0
    if mmax >= 1:
        H[1] = x
    for m in range(1, mmax):
        H[m + 1] = x * H[m] - m * H[m - 1]
    return H


def relu_pow_table(mu, var, p, kmax):
    """F[k, i] = E[(relu^p)^{(k)}(Y_i)], Y_i ~ N(mu_i, var_i), k = 0..kmax (distributional derivatives)."""
    sig = np.sqrt(var); al = mu / sig
    phi = np.exp(-0.5 * al * al) / np.sqrt(2 * pi); Phi = ndtr(al)
    n = mu.shape[0]
    F = np.zeros((kmax + 1, n))
    He = herme_vals(-al, max(kmax - p - 1, 0))
    for k in range(p + 1, kmax + 1):
        F[k] = factorial(p) * sig ** (-(k - p)) * He[k - p - 1] * phi
    # k <= p: sigma^{p-k} int_{-alpha}^inf (alpha+u)^p He_k(u) phi(u) du via truncated moments I_m
    I = np.zeros((2 * p + 1, n)); I[0] = Phi
    if 2 * p >= 1:
        I[1] = phi
    for m in range(2, 2 * p + 1):
        I[m] = (m - 1) * I[m - 2] + (-al) ** (m - 1) * phi
    from numpy.polynomial.hermite_e import herme2poly
    from math import comb
    for k in range(0, min(p, kmax) + 1):
        hk = herme2poly([0] * k + [1])          # power-basis coefficients of He_k(u)
        acc = np.zeros(n)
        for r in range(p + 1):
            cr = comb(p, r) * al ** (p - r)
            for s, hs in enumerate(hk):
                if hs != 0.0:
                    acc += cr * hs * I[r + s]
        F[k] = sig ** (p - k) * acc
    return F


# ----------------------------------------------------------------------------------------------------------
# Edgeworth operators (dicts {(r, s): coefficient array}) and the 2-D Mehler expectation
# ----------------------------------------------------------------------------------------------------------

def op_mul(A, B):
    out = {}
    for (r1, s1), a in A.items():
        for (r2, s2), b in B.items():
            key = (r1 + r2, s1 + s2)
            out[key] = out.get(key, 0.0) + a * b
    return out


def op_add(*ops):
    out = {}
    for O in ops:
        for k, v in O.items():
            out[k] = out.get(k, 0.0) + v
    return out


def mehler(C, F, G, r, s, kmax=KMAX):
    """M_ij = E_G[f^{(r)}(z_i) g^{(s)}(z_j)] = sum_k C_ij^k / k! F_i(r+k) G_j(s+k)."""
    out = np.zeros_like(C)
    Ck = np.ones_like(C)
    for k in range(kmax + 1):
        out += Ck * np.outer(F[r + k], G[s + k])
        Ck = Ck * C / (k + 1)
    return out


class State:
    """pre-activation cumulants of one layer: mu, C (n x n), K3 (n^3) or None, X = kappa4_{aabc} (n^3) or None."""
    def __init__(self, mu, C, K3=None, X=None):
        self.mu, self.C, self.K3, self.X = mu, C, K3, X

    @property
    def var(self):
        return np.diag(self.C).copy()


def slices_from_state(st):
    n = st.mu.shape[0]
    idx = np.arange(n)
    if st.K3 is not None:
        D3 = st.K3[idx, idx, idx].copy(); D21 = st.K3[idx, idx, :].copy(); np.fill_diagonal(D21, 0.0)
    else:
        D3 = np.zeros(n); D21 = np.zeros((n, n))
    if st.X is not None:
        K4 = st.X[idx, idx, idx].copy()
        K31 = st.X[idx, idx, :].copy(); np.fill_diagonal(K31, 0.0)    # kappa4_{aaab}
        K22 = st.X[:, idx, idx].copy(); np.fill_diagonal(K22, 0.0)    # kappa4_{aabb}
    else:
        K4 = np.zeros(n); K31 = np.zeros((n, n)); K22 = np.zeros((n, n))
    return D3, D21, K4, K31, K22


def relu_moments(st, pmax=4, order=2, k4=True):
    """raw 1-D moments E[a^p] (p=1..pmax) and pair raw moments E[a_i^p a_j^q] (p+q <= pmax, p,q >= 1).
    order: 0 = Gaussian, 1 = first order in kappa3/kappa4, 2 = also kappa3^2/2."""
    mu, C, var = st.mu, st.C, st.var
    n = mu.shape[0]
    D3, D21, K4, K31, K22 = slices_from_state(st)
    if not k4:
        K4 = K4 * 0; K31 = K31 * 0; K22 = K22 * 0
    kneed = KMAX + 7
    F = {p: relu_pow_table(mu, var, p, kneed) for p in range(1, pmax + 1)}
    # 1-D
    m1d = {}
    for p in range(1, pmax + 1):
        v = F[p][0].copy()
        if order >= 1:
            v += D3 / 6 * F[p][3] + K4 / 24 * F[p][4]
        if order >= 2:
            v += D3 ** 2 / 72 * F[p][6]
        m1d[p] = v
    # 2-D operator
    one = {(0, 0): 1.0}
    O3 = {(3, 0): D3[:, None] / 6, (2, 1): D21 / 2, (1, 2): D21.T / 2, (0, 3): D3[None, :] / 6}
    O4 = {(4, 0): K4[:, None] / 24, (3, 1): K31 / 6, (2, 2): K22 / 4, (1, 3): K31.T / 6, (0, 4): K4[None, :] / 24}
    if order == 0:
        Op = one
    elif order == 1:
        Op = op_add(one, O3, O4)
    else:
        O3sq = op_mul(O3, O3)
        Op = op_add(one, O3, O4, {k: 0.5 * v for k, v in O3sq.items()})
    pair = {}
    for p in range(1, pmax):
        for q in range(1, pmax - p + 1):
            acc = np.zeros((n, n))
            for (r, s), c in Op.items():
                acc += c * mehler(C, F[p], F[q], r, s)
            pair[(p, q)] = acc
    gate = F[1][1].copy()
    if order >= 1:
        gate += D3 / 6 * F[1][4] + K4 / 24 * F[1][5]
    if order >= 2:
        gate += D3 ** 2 / 72 * F[1][7]
    return m1d, pair, F, gate


def central_pair(m1d, pair, n, pmax=4):
    """central pair moments mu_pq (i != j entries meaningful) from raw ones."""
    from math import comb
    m = m1d[1]
    R = {}
    for p in range(0, pmax + 1):
        for q in range(0, pmax + 1 - p):
            if p == 0 and q == 0:
                R[(0, 0)] = np.ones((n, n))
            elif q == 0:
                R[(p, 0)] = np.broadcast_to(m1d[p][:, None], (n, n))
            elif p == 0:
                R[(0, q)] = np.broadcast_to(m1d[q][None, :], (n, n))
            else:
                R[(p, q)] = pair[(p, q)]
    mi = m[:, None]; mj = m[None, :]
    cen = {}
    for p in range(0, pmax + 1):
        for q in range(0, pmax + 1 - p):
            acc = np.zeros((n, n))
            for r in range(p + 1):
                for s in range(q + 1):
                    acc = acc + comb(p, r) * comb(q, s) * (-mi) ** (p - r) * (-mj) ** (q - s) * R[(r, s)]
            cen[(p, q)] = acc
    return cen


# ----------------------------------------------------------------------------------------------------------
# transports
# ----------------------------------------------------------------------------------------------------------

def transport3(K, W):
    T = np.tensordot(K, W, axes=([2], [0]))                 # ijc
    T = np.einsum("ijc,jb->ibc", T, W, optimize=True)
    return np.einsum("ibc,ia->abc", T, W, optimize=True)


def transport_k4_slices(G, T, S, W):
    """X_abc = kappa4(z)_{aabc} for z = a @ W when kappa4(a) is supported on its (4), (3,1), (2,2) slices:
    G_i = kappa4_iiii, T_ij = kappa4_iiij (i != j), S_ij = kappa4_iijj (i != j, symmetric)."""
    T = T.copy(); np.fill_diagonal(T, 0.0); S = S.copy(); np.fill_diagonal(S, 0.0)
    W2 = W * W
    X = np.einsum("ia,ib,ic->abc", W2 * G[:, None], W, W, optimize=True)
    TW = T @ W
    X += 2 * np.einsum("ia,ib,ic->abc", TW * W, W, W, optimize=True)
    X += np.einsum("ia,ib,ic->abc", W2, TW, W, optimize=True)
    X += np.einsum("ia,ib,ic->abc", W2, W, TW, optimize=True)
    Q = np.einsum("ij,jb,jc->ibc", S, W, W, optimize=True)         # Q_i,bc = sum_j S_ij W_jb W_jc
    X += np.einsum("ia,ibc->abc", W2, Q, optimize=True)
    X += np.einsum("ia,ib,iac->abc", W, W, Q, optimize=True)
    X += np.einsum("ia,ic,iab->abc", W, W, Q, optimize=True)
    return X


# ----------------------------------------------------------------------------------------------------------
# atlas helpers (truth objects for validation and teacher forcing)
# ----------------------------------------------------------------------------------------------------------

def atlas_state(z, l, with_k3=True, with_k4=True):
    """exact pre-activation cumulants of layer l from an atlas npz (moment_atlas_np.py --k3 [--k4])."""
    pre_s = z["pre_s"]
    mu, m2 = pre_s[0, l], pre_s[1, l]
    M11 = z["pre_M11"][l].astype(np.float64)
    C = M11 - np.outer(mu, mu)
    K3 = ok.central3(z["pre_M3"][l], M11, mu) if (with_k3 and "pre_M3" in z.files) else None
    X = None
    if with_k4 and "pre_M211" in z.files:
        M21 = z["pre_M21"][l].astype(np.float64)
        M3 = z["pre_M3"][l]
        M211 = z["pre_M211"][l]
        var = m2 - mu ** 2
        Eu2uu = (M211 - np.einsum("k,ij->ijk", mu, M21) - np.einsum("j,ik->ijk", mu, M21) + np.einsum("j,k,i->ijk", mu, mu, m2)
                 - 2 * np.einsum("i,ijk->ijk", mu, M3 - np.einsum("k,ij->ijk", mu, M11) - np.einsum("j,ik->ijk", mu, M11)
                                 + np.einsum("i,j,k->ijk", mu, mu, mu))
                 + np.einsum("i,jk->ijk", mu ** 2, M11 - np.outer(mu, mu)))
        X = Eu2uu - np.einsum("i,jk->ijk", var, C) - 2 * np.einsum("ij,ik->ijk", C, C)
    return State(mu.astype(np.float64), C, K3, X)


def atlas_post(z, l):
    """exact post-activation mean, covariance and kappa3 of layer l."""
    mu = z["post_s"][0, l].astype(np.float64)
    M11 = z["post_M11"][l].astype(np.float64)
    C = M11 - np.outer(mu, mu)
    np.fill_diagonal(C, z["post_s"][1, l] - mu ** 2)
    K3 = ok.central3(z["post_M3"][l], M11, mu) if "post_M3" in z.files else None
    return mu, C, K3


# ----------------------------------------------------------------------------------------------------------
# the ReLU step
# ----------------------------------------------------------------------------------------------------------

def relu_step(st, k3mode="closure", coef=None, order=2, k4carry=True, phi_mode="corrected", herm_deg=4,
              want_k4=True):
    """returns (mu_a, C_a, K3_a or None, (G, T, S) kappa4(a) slices or None, info)."""
    n = st.mu.shape[0]
    k2 = st.K3 is None and k3mode == "k2"
    pmax = 4 if want_k4 else 3
    m1d, pair, F, gate = relu_moments(st, pmax=pmax, order=(0 if k2 else order))
    cen = central_pair(m1d, pair, n, pmax=pmax)
    mu_a = m1d[1]
    var_a = m1d[2] - mu_a ** 2
    C_a = cen[(1, 1)].copy(); np.fill_diagonal(C_a, var_a)
    info = dict(gate=gate)
    if k2:
        return mu_a, C_a, None, None, info
    # third cumulant slices
    D3a = m1d[3] - 3 * m1d[2] * mu_a + 2 * mu_a ** 3
    D21a = cen[(2, 1)].copy(); np.fill_diagonal(D21a, 0.0)      # kappa3(a)_{iij}
    # all-distinct part
    mu, var, C = st.mu, st.var, st.C
    if k3mode == "none":
        Kd = np.zeros((n, n, n))
    else:
        D3, D21, K4, K31, K22 = slices_from_state(st)
        Phi = gate if phi_mode == "corrected" else ndtr(mu / np.sqrt(var))
        K3z = st.K3 if st.K3 is not None else np.zeros((n, n, n))
        o = dict(mu=mu, var=var, C=C, K3z=K3z, Phi=Phi, K22=K22)
        if st.X is not None:
            o["K211"] = ok.all_distinct(st.X)
        else:
            o["K211"] = np.zeros((n, n, n))
        H = ok.hermite_model(C, mu, var, herm_deg)
        if k3mode == "wick":
            c = [1.0, 0, 0, 0, 0, 0, 0]
        elif k3mode == "closure":
            c = ok.CLOSURE_COEF
        elif k3mode == "fit":
            c = coef
        else:
            raise ValueError(k3mode)
        B = ok.residual_basis(o) if any(abs(x) > 0 for x in c[1:]) else [ok.all_distinct(np.einsum("i,j,k,ijk->ijk", Phi, Phi, Phi, K3z))]
        Kd = H + sum(ci * b for ci, b in zip(c, B))
    K3a = Kd
    idx = np.arange(n)
    K3a[idx, idx, :] = D21a; K3a[idx, :, idx] = D21a; K3a[:, idx, idx] = D21a.T
    K3a[idx, idx, idx] = D3a
    k4s = None
    if want_k4:
        G = m1d[4] - 4 * m1d[3] * mu_a + 6 * m1d[2] * mu_a ** 2 - 3 * mu_a ** 4 - 3 * var_a ** 2
        T = cen[(3, 1)] - 3 * cen[(2, 0)] * cen[(1, 1)]
        S = cen[(2, 2)] - cen[(2, 0)] * cen[(0, 2)] - 2 * cen[(1, 1)] ** 2
        np.fill_diagonal(T, 0.0); np.fill_diagonal(S, 0.0)
        k4s = (G, T, S)
    return mu_a, C_a, K3a, k4s, info


# ----------------------------------------------------------------------------------------------------------
# the chain
# ----------------------------------------------------------------------------------------------------------

class Chain:
    """opts:
      k3mode   'k2' (A: covariance propagation, Gaussian closure to all orders in C),
               'none' | 'wick' | 'closure' | 'fit'  (dense K=3 with the all-distinct birth rule)
      coefs    (L, 7) table for 'fit'
      k4mode   'mem'   kappa4(z) from the transported (4),(3,1),(2,2) slices of kappa4(a) (memoryless kappa4)
               'zero'  kappa4(z) = 0
               'atlas' kappa4(z)_{aabc} from the atlas at every layer (teacher-forced fourth cumulant)
               'atlas211'   only the (2,1,1) all-distinct part from the atlas, the rest 'mem'
               'reg211'     the (2,1,1) part regenerated as u_i C_jk (u fitted per layer on the atlas slice), rest 'mem'
               'zero211'    'mem' with its (2,1,1) all-distinct part zeroed
      order    Edgeworth order (0, 1, 2)
      force    dict layer -> what to replace from the atlas after that layer's ReLU: 'k3' (kappa3(a)), 'all'
               (mu, C, kappa3 of a), 'k3d' (only the all-distinct part of kappa3(a))
      d21_noise  (eps, seed): multiply the all-distinct kappa3(a) contribution to D21(l+1) ... implemented as a
               relative perturbation of D21(z_{l+1}) for the error law
      atlas    npz (needed for k4mode atlas*, force)
    """
    def __init__(self, W, k3mode="closure", coefs=None, k4mode="mem", order=2, force=None, atlas=None,
                 d21_noise=None, phi_mode="corrected", herm_deg=4, record=True):
        self.W = np.asarray(W, dtype=np.float64)
        self.k3mode, self.coefs, self.k4mode, self.order = k3mode, coefs, k4mode, order
        self.force = force or {}
        self.atlas = atlas
        self.d21_noise = d21_noise
        self.phi_mode, self.herm_deg = phi_mode, herm_deg
        self.record = record

    def run(self):
        W = self.W
        L, n, _ = W.shape
        k2 = self.k3mode == "k2"
        st = State(np.zeros(n), W[0].T @ W[0], None if k2 else np.zeros((n, n, n)), None if k2 else np.zeros((n, n, n)))
        means = np.zeros((L, n))
        rec = []
        rng = np.random.default_rng(self.d21_noise[1]) if self.d21_noise else None
        for l in range(L):
            if self.k4mode.startswith("atlas") or self.k4mode in ("reg211",):
                if self.atlas is not None and l > 0 and st.X is not None:
                    Xt = atlas_state(self.atlas, l, with_k3=False).X
                    if self.k4mode == "atlas":
                        st.X = Xt
                    elif self.k4mode == "atlas211":
                        st.X = (st.X - ok.all_distinct(st.X)) + ok.all_distinct(Xt)
                    elif self.k4mode == "reg211":
                        K211 = ok.all_distinct(Xt)
                        Co = ok.offdiag(st.C)
                        u = np.einsum("ijk,jk->i", K211, Co) / float(np.sum(Co * Co))
                        st.X = (st.X - ok.all_distinct(st.X)) + ok.all_distinct(np.einsum("i,jk->ijk", u, Co))
            if self.k4mode == "zero211" and st.X is not None:
                st.X = st.X - ok.all_distinct(st.X)
            if self.k4mode == "zero" and st.X is not None:
                st.X = np.zeros_like(st.X)
            coef = self.coefs[l] if self.coefs is not None else None
            want_k4 = (not k2) and l < L - 1
            mu_a, C_a, K3a, k4s, info = relu_step(st, self.k3mode, coef, self.order, phi_mode=self.phi_mode,
                                                  herm_deg=self.herm_deg, want_k4=want_k4)
            if self.record:
                D3, D21, K4, K31, K22 = slices_from_state(st)
                rec.append(dict(mu=st.mu.copy(), var=st.var, D21=D21, D3=D3, K4=K4, K22=K22, C=st.C.copy()))
            if l in self.force and self.atlas is not None:
                mu_t, C_t, K3_t = atlas_post(self.atlas, l)
                what = self.force[l]
                if what == "all":
                    mu_a, C_a, K3a = mu_t, C_t, K3_t.copy()
                elif what == "k3":
                    K3a = K3_t.copy()
                elif what == "k3d":
                    K3a = (K3a - ok.all_distinct(K3a)) + ok.all_distinct(K3_t)
            means[l] = mu_a
            if l == L - 1:
                break
            Wn = W[l + 1]
            mu_z = mu_a @ Wn
            C_z = Wn.T @ C_a @ Wn
            if k2:
                st = State(mu_z, C_z, None, None)
                continue
            K3z = transport3(K3a, Wn)
            if self.d21_noise is not None and self.d21_noise[0] > 0:
                eps = self.d21_noise[0]
                idx = np.arange(n)
                D21 = K3z[idx, idx, :].copy(); np.fill_diagonal(D21, 0.0)
                rms = np.sqrt(np.mean(D21 ** 2))
                E = rng.standard_normal((n, n)) * eps * rms; np.fill_diagonal(E, 0.0)
                # perturb kappa3(z)_{aab} and its symmetric images (the (2,1) slice) by E
                K3z[idx, idx, :] += E; K3z[idx, :, idx] += E; K3z[:, idx, idx] += E.T
            X = transport_k4_slices(*k4s, Wn) if self.k4mode != "zero" else np.zeros((n, n, n))
            st = State(mu_z, C_z, K3z, X)
        return dict(means=means, rec=rec)


def paper_k2(W):
    """the paper's Algorithm 2 (covariance propagation, leading-order off-diagonal update), for reference."""
    W = np.asarray(W, dtype=np.float64)
    L, n, _ = W.shape
    mu = np.zeros(n); C = W[0].T @ W[0]
    means = np.zeros((L, n))
    for l in range(L):
        var = np.diag(C).copy()
        F1 = relu_pow_table(mu, var, 1, 1); F2 = relu_pow_table(mu, var, 2, 0)
        a, b, c = F1[0], F2[0], F1[1]
        Cn = C * np.outer(c, c); np.fill_diagonal(Cn, b - a * a)
        means[l] = a
        if l == L - 1:
            break
        mu = a @ W[l + 1]; C = W[l + 1].T @ Cn @ W[l + 1]
    return means
