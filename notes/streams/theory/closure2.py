"""Second-order diagram closure of the all-distinct post-activation third cumulant kappa3(a_i, a_j, a_k), a = relu(z).

Theory (REPORT.md, section 1).  For a near-Gaussian z with mean mu, covariance C and higher cumulants kappa_m, the
joint cumulant of f(z_i), f(z_j), f(z_k) on distinct i, j, k is the sum over CONNECTED hypergraphs on the three
vertices {i, j, k} whose hyperedges are C edges (two distinct vertices) and cumulant hyperedges kappa_m (m >= 3, any
leg multiplicities r_v, at least two distinct vertices; single-vertex hyperedges are the marginal and are absorbed in
the vertex weights).  The weight of one labelled hypergraph is

    prod_h kappa_h / prod_v r_{h,v}!   x   prod_types 1/(multiplicity)!   x   prod_v w_v(deg v),

    w_v(d) = E[f^(d)(z_v)]  (degree-1 vertices: the true gate P(z_v > 0); d >= 2: Gaussian relu value
             He_{d-2}(-alpha) phi(alpha) / sigma^{d-1}, plus the marginal dressing where it enters at the order kept).

Equivalently (oracle_k3's phrasing) each vertex of degree d carries w(d)/d! times the number of leg partitions.  A
term below is written coef * sym3[X] with X one labelled representative; coef = 6/|Aut(X)| x prod 1/r! x prod 1/mult!.

Order counting at width n, eps = n^-1/2 (fixed random He-init MLP, measured in REPORT.md): C_off ~ eps; every
kappa3 slice ~ eps^2; kappa4 slices with all multiplicities even ((2,2), (4)) ~ eps^2, others ((2,1,1), (3,1)) ~ eps^3;
kappa5 ~ eps^4; kappa6 (2,2,2) ~ eps^4.  The target all-distinct kappa3(a) is O(eps^2).

    leading (eps^2):      Wick C C path, Phi^3 kappa3_ijk                                    (oracle_k3.wick_model)
    first order (eps^3):  B1 B2 (D21 + edge), B4 (Gaussian rho^3), B5 (K22 + edge), B6 (K211),  B7 (kappa3_ijk + edge,
                          MISSING from oracle_k3's basis); B3 (D3 dressing of the Wick vertex) is eps^4 and its
                          coefficient is 0.5, not 1.0 as in oracle_k3.CLOSURE_COEF.
    second order (eps^4): S-terms below (30 shapes; three need atlas fields that moment_atlas_np.py does not build).

Usage
    python closure2.py --pair A.npz B.npz [--layers 0-14] [--out results/x.txt]   pair evaluation (moment_atlas_np format)
    python closure2.py --pair56 A.npz B.npz                                       same for atlas_k56.py atlases (+ kappa5/6)
    python closure2.py --selftest                                                 n = 3 identity against the enumerator
"""
import os
import sys
import time
from math import sqrt, pi, erf

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "experiments"))
import oracle_k3 as ok  # noqa: E402

all_distinct, sym3, offdiag, rel, hermite_model = ok.all_distinct, ok.sym3, ok.offdiag, ok.rel, ok.hermite_model


# ---------------------------------------------------------------------------------------------------------------
# vertex weights
# ---------------------------------------------------------------------------------------------------------------
def relu_w(mu, var, dmax=8):
    """Gaussian vertex weights w(d) = E[relu^(d)(z)], z ~ N(mu, var): w0 = E relu, w1 = Phi(alpha),
    w(d) = He_{d-2}(-alpha) phi(alpha) / sigma^{d-1} for d >= 2."""
    from numpy.polynomial.hermite_e import hermeval
    sig = np.sqrt(var); al = mu / sig
    phi = np.exp(-0.5 * al ** 2) / sqrt(2 * pi)
    Phi = 0.5 * (1.0 + np.vectorize(erf)(al / sqrt(2)))
    w = {0: sig * (al * Phi + phi), 1: Phi}
    for d in range(2, dmax + 1):
        w[d] = hermeval(-al, [0] * (d - 2) + [1]) * phi / sig ** (d - 1)
    return w


# ---------------------------------------------------------------------------------------------------------------
# the diagram terms
# ---------------------------------------------------------------------------------------------------------------
def E(spec, *ops):
    return np.einsum(spec, *ops)


# order tags: 'L' leading eps^2, '1' first order eps^3 (oracle basis), '1+' first order missing from the oracle basis,
# '2' second order eps^4 available from moment_atlas_np atlases, '2new' second order needing new atlas fields
TERM_INFO = {
    # name: (order, group, coef, description)
    "B0": ("L", "B", 1.0, "Phi_i Phi_j Phi_k kappa3_ijk"),
    "B1": ("1", "B", 3.0, "D21_ik + C_jk: w2_i Phi_j w2_k"),
    "B2": ("1", "B", 3.0, "D21_ik + C_ij: w3_i Phi_j Phi_k"),
    "B3": ("2", "B", 0.5, "D3_i dressing of the Wick vertex: w5_i Phi_j Phi_k C_ij C_ik (oracle_k3 uses 1.0)"),
    "B4": ("1", "B", 1.0, "Gaussian rho^3 (hermite degree 6 - 4)"),
    "B5": ("1", "B", 1.5, "K22_ij + C_ik: w3_i w2_j Phi_k"),
    "B6": ("1", "B", 1.5, "K211_i;jk: w2_i Phi_j Phi_k"),
    "B7": ("1+", "T+e", 3.0, "kappa3_ijk + C_ij: w2_i w2_j Phi_k"),
    "S1": ("2", "rho4", 1.0, "Gaussian rho^4 (hermite degree 8 - 6)"),
    "S2a": ("2", "dress", 1.0, "Wick with true gates P(z>0) instead of Phi(alpha) at the degree-1 vertices"),
    "S2c": ("2", "dress", 0.125, "kappa4_iiii dressing of the Wick vertex: w6_i Phi_j Phi_k C_ij C_ik"),
    "S3a1": ("2", "k3+2e", 1.5, "kappa3_ijk + C_ij^2: w3_i w3_j Phi_k"),
    "S3a2": ("2", "k3+2e", 3.0, "kappa3_ijk + C_ij C_ik: w3_i w2_j w2_k"),
    "S3b1": ("2", "k3+2e", 1.5, "D21_ik + C_jk^2: w2_i w2_j w3_k"),
    "S3b2": ("2", "k3+2e", 1.5, "D21_ik + C_ij^2: w4_i w2_j Phi_k"),
    "S3b3": ("2", "k3+2e", 3.0, "D21_ik + C_ij C_jk: w3_i w2_j w2_k"),
    "S3b4": ("2", "k3+2e", 3.0, "D21_ik + C_ij C_ik: w4_i Phi_j w2_k"),
    "S3b5": ("2", "k3+2e", 3.0, "D21_ik + C_jk C_ik: w3_i Phi_j w3_k"),
    "S4a": ("2", "k3xk3", 0.5, "kappa3_ijk^2: w2 w2 w2"),
    "S4b": ("2", "k3xk3", 3.0, "kappa3_ijk D21_ik: w3_i Phi_j w2_k"),
    "S4c": ("2", "k3xk3", 0.75, "D21_ik D21_jk: w2 w2 w2"),
    "S4d": ("2", "k3xk3", 0.75, "D21_ki D21_kj: Phi_i Phi_j w4_k"),
    "S4e": ("2", "k3xk3", 1.5, "D21_ik D21_kj: w2_i Phi_j w3_k"),
    "S5a": ("2", "k4odd+e", 1.5, "K211_i;jk + C_jk: w2 w2 w2"),
    "S5b": ("2", "k4odd+e", 3.0, "K211_i;jk + C_ij: w3_i w2_j Phi_k"),
    "S5c": ("2", "k4odd+e", 1.0, "K31_ij (kappa4_iiij) + C_ik: w4_i Phi_j Phi_k"),
    "S5d": ("2", "k4odd+e", 1.0, "K31_ij + C_jk: w3_i w2_j Phi_k"),
    "S6a": ("2", "K22+2e", 0.75, "K22_ij + C_ik^2: w4_i w2_j w2_k"),
    "S6b": ("2", "K22+2e", 0.75, "K22_ij + C_ik C_jk: w3_i w3_j w2_k"),
    "S6c": ("2", "K22+2e", 1.5, "K22_ij + C_ik C_ij: w4_i w3_j Phi_k"),
    "S7a": ("2", "K22xk3", 0.75, "K22_ij kappa3_ijk: w3_i w3_j Phi_k"),
    "S7b": ("2", "K22xk3", 0.75, "K22_ij D21_ik: w4_i w2_j Phi_k"),
    "S7c": ("2", "K22xk3", 0.75, "K22_ij D21_ki: w3_i w2_j w2_k"),
    "S8": ("2", "K22xK22", 0.1875, "K22_ij K22_jk: w2_i w4_j w2_k"),
    "S9a": ("2new", "k5", 0.75, "kappa5 (2,2,1) slice kappa5_iijjk: w2_i w2_j Phi_k"),
    "S9b": ("2new", "k5", 0.5, "kappa5 (3,1,1) slice kappa5_iiijk: w3_i Phi_j Phi_k"),
    "S10": ("2new", "k6", 0.125, "kappa6 (2,2,2) slice kappa6_iijjkk: w2 w2 w2"),
    # leaf-resummed closure (section 5 of REPORT.md): every diagram in which some vertex k is attached by a single
    # C edge to i sums to Phi_k C_ik Gamma_ij, Gamma_ij = Cov(1[z_i>0], a_j) (exact pair object); diagrams with two
    # such leaves on one center are counted twice and are subtracted exactly with the true density p_i(0) = E relu''.
    "LEAF": ("R", "leaf", 6.0, "sym3[Phi_k C_ik Gamma_ij], Gamma_ij = Cov(1[z_i>0], relu(z_j))  (all C-leaf diagrams)"),
    "TWOLEAF": ("R", "leaf", -1.0, "Wick form with the true density p_i(0) at the center (double-counted two-leaf class)"),
    "TWOLEAFe": ("R", "leaf", -1.0, "same with the Edgeworth density w2 + D3/6 w5 + K4/24 w6"),
    "TCL": ("R", "tclass", 1.0, "kappa3_ijk (Phi_i Phi_j Phi_k + sum_pairs Phi_k Cov(1[z_i>0], 1[z_j>0])): every diagram with one "
            "kappa3_ijk hyperedge and the rest on one pair (B0, B7, S3a1, S4b, S7a, ... to all orders); needs gate_GG"),
    "G3t": ("1", "gauss-nonleaf", 1.0, "Gaussian rho^3 triangle C_ij C_jk C_ik: w2 w2 w2"),
    "G4a": ("2", "gauss-nonleaf", 0.75, "Gaussian rho^4 C_ij^2 C_jk^2: w2_i w4_j w2_k"),
    "G4b": ("2", "gauss-nonleaf", 1.5, "Gaussian rho^4 C_ij^2 C_jk C_ik: w3_i w3_j w2_k"),
}
NONLEAF1 = ["B6", "B7", "G3t"]
NONLEAF2 = ["S3a1", "S3a2", "S3b1", "S3b2", "S3b3", "S4a", "S4b", "S4c", "S4d", "S4e", "S5a", "S5b", "S6a", "S6b",
            "S7a", "S7b", "S7c", "S8", "G4a", "G4b"]
FIRST = ["B0", "B1", "B2", "B3", "B4", "B5", "B6"]          # oracle_k3.residual_basis order
ORACLE_COEF = [1.0, 3.0, 3.0, 1.0, 1.0, 1.5, 1.5]            # oracle_k3.CLOSURE_COEF (B3 overcounted by 2)
LEG_COEF = [TERM_INFO[t][2] for t in FIRST]                  # corrected: B3 = 0.5
SECOND_AVAIL = [t for t, v in TERM_INFO.items() if v[0] == "2" and t not in ("B3", "G4a", "G4b")]
SECOND_NEW = [t for t, v in TERM_INFO.items() if v[0] == "2new"]


def terms(o, names=None):
    """dict name -> all-distinct n^3 tensor (UNSCALED: the closure adds TERM_INFO[name][2] * tensor).
    o: mu, var, C (full), Phi (degree-1 vertex weight), w (dict d -> Gaussian w(d), d >= 2; w[1] Gaussian Phi(alpha)),
       T (kappa3 full n^3), K4f (kappa4_{iijk} full n^3), optional P (kappa5_{iijjk}), Q (kappa5_{iiijk}), S (kappa6_{iijjkk})."""
    n = len(o["mu"])
    idx = np.arange(n)
    C, T, w, Phi = o["C"], o["T"], o["w"], o["Phi"]
    Co = offdiag(C)
    Td = all_distinct(T)
    Dm = T[idx, idx, :]                        # D[i, k] = kappa3_iik
    D3 = np.diag(Dm).copy(); Do = offdiag(Dm)
    K4f = o.get("K4f")
    if K4f is not None:
        U = all_distinct(K4f)
        K22 = offdiag(K4f[:, idx, idx])           # kappa4_iijj
        V = offdiag(K4f[idx, idx, :])             # kappa4_iiij  (i tripled)
        K4 = K4f[idx, idx, idx].copy()
    w2, w3, w4, w5, w6 = w[2], w[3], w[4], w[5], w[6]
    want = set(names) if names is not None else None
    out = {}

    def add(name, fn):
        if want is None or name in want:
            out[name] = all_distinct(fn())

    add("B0", lambda: E("i,j,k,ijk->ijk", Phi, Phi, Phi, T))
    add("B1", lambda: sym3(E("i,j,k,ik,jk->ijk", w2, Phi, w2, Do, Co)))
    add("B2", lambda: sym3(E("i,j,k,ik,ij->ijk", w3, Phi, Phi, Do, Co)))
    add("B3", lambda: sym3(E("i,j,k,i,ij,ik->ijk", w5, Phi, Phi, D3, Co, Co)))
    add("B4", lambda: hermite_model(C, o["mu"], o["var"], 6) - hermite_model(C, o["mu"], o["var"], 4))
    add("B7", lambda: sym3(E("i,j,k,ijk,ij->ijk", w2, w2, Phi, Td, Co)))
    add("S1", lambda: hermite_model(C, o["mu"], o["var"], 8) - hermite_model(C, o["mu"], o["var"], 6))
    add("S2a", lambda: ok.wick_model(C, Phi, w2, np.zeros_like(T)) - hermite_model(C, o["mu"], o["var"], 4))
    add("S3a1", lambda: sym3(E("i,j,k,ijk,ij->ijk", w3, w3, Phi, Td, Co * Co)))
    add("S3a2", lambda: sym3(E("i,j,k,ijk,ij,ik->ijk", w3, w2, w2, Td, Co, Co)))
    add("S3b1", lambda: sym3(E("i,j,k,ik,jk->ijk", w2, w2, w3, Do, Co * Co)))
    add("S3b2", lambda: sym3(E("i,j,k,ik,ij->ijk", w4, w2, Phi, Do, Co * Co)))
    add("S3b3", lambda: sym3(E("i,j,k,ik,ij,jk->ijk", w3, w2, w2, Do, Co, Co)))
    add("S3b4", lambda: sym3(E("i,j,k,ik,ij,ik->ijk", w4, Phi, w2, Do, Co, Co)))
    add("S3b5", lambda: sym3(E("i,j,k,ik,jk,ik->ijk", w3, Phi, w3, Do, Co, Co)))
    add("S4a", lambda: E("i,j,k,ijk->ijk", w2, w2, w2, Td * Td))
    add("S4b", lambda: sym3(E("i,j,k,ijk,ik->ijk", w3, Phi, w2, Td, Do)))
    add("S4c", lambda: sym3(E("i,j,k,ik,jk->ijk", w2, w2, w2, Do, Do)))
    add("S4d", lambda: sym3(E("i,j,k,ki,kj->ijk", Phi, Phi, w4, Do, Do)))
    add("S4e", lambda: sym3(E("i,j,k,ik,kj->ijk", w2, Phi, w3, Do, Do)))
    if K4f is not None:
        add("B5", lambda: sym3(E("i,j,k,ij,ik->ijk", w3, w2, Phi, K22, Co)))
        add("B6", lambda: sym3(E("i,j,k,ijk->ijk", w2, Phi, Phi, U)))
        add("S2c", lambda: sym3(E("i,j,k,i,ij,ik->ijk", w6, Phi, Phi, K4, Co, Co)))
        add("S5a", lambda: sym3(E("i,j,k,ijk,jk->ijk", w2, w2, w2, U, Co)))
        add("S5b", lambda: sym3(E("i,j,k,ijk,ij->ijk", w3, w2, Phi, U, Co)))
        add("S5c", lambda: sym3(E("i,j,k,ij,ik->ijk", w4, Phi, Phi, V, Co)))
        add("S5d", lambda: sym3(E("i,j,k,ij,jk->ijk", w3, w2, Phi, V, Co)))
        add("S6a", lambda: sym3(E("i,j,k,ij,ik->ijk", w4, w2, w2, K22, Co * Co)))
        add("S6b", lambda: sym3(E("i,j,k,ij,ik,jk->ijk", w3, w3, w2, K22, Co, Co)))
        add("S6c", lambda: sym3(E("i,j,k,ij,ik,ij->ijk", w4, w3, Phi, K22, Co, Co)))
        add("S7a", lambda: sym3(E("i,j,k,ij,ijk->ijk", w3, w3, Phi, K22, Td)))
        add("S7b", lambda: sym3(E("i,j,k,ij,ik->ijk", w4, w2, Phi, K22, Do)))
        add("S7c", lambda: sym3(E("i,j,k,ij,ki->ijk", w3, w2, w2, K22, Do)))
        add("S8", lambda: sym3(E("i,j,k,ij,jk->ijk", w2, w4, w2, K22, K22)))
    add("G3t", lambda: E("i,j,k,ij,jk,ik->ijk", w2, w2, w2, Co, Co, Co))
    add("G4a", lambda: sym3(E("i,j,k,ij,jk->ijk", w2, w4, w2, Co * Co, Co * Co)))
    add("G4b", lambda: sym3(E("i,j,k,ij,jk,ik->ijk", w3, w3, w2, Co * Co, Co, Co)))
    if "GG" in o:
        Gc = offdiag(o["GG"] - np.outer(Phi, Phi))
        add("TCL", lambda: Td * (E("i,j,k->ijk", Phi, Phi, Phi) + E("k,ij->ijk", Phi, Gc) + E("i,jk->ijk", Phi, Gc)
                                 + E("j,ik->ijk", Phi, Gc)))
    if "Gam" in o:
        add("LEAF", lambda: sym3(E("k,ik,ij->ijk", Phi, Co, offdiag(o["Gam"]))))
        add("TWOLEAF", lambda: ok.wick_model(C, Phi, o["w2hat"], np.zeros_like(T)))
        if K4f is not None:
            add("TWOLEAFe", lambda: ok.wick_model(C, Phi, w2 + D3 / 6 * w5 + K4 / 24 * w6, np.zeros_like(T)))
    if "P" in o:
        add("S9a", lambda: sym3(E("i,j,k,ijk->ijk", w2, w2, Phi, o["P"])))
    if "Q" in o:
        add("S9b", lambda: sym3(E("i,j,k,ijk->ijk", w3, Phi, Phi, o["Q"])))
    if "S" in o:
        add("S10", lambda: E("i,j,k,ijk->ijk", w2, w2, w2, o["S"]))
    return out


def wick_gauss(o):
    """the leading Gaussian rho^2 term with Gaussian gates Phi(alpha) (oracle_k3.closure_model's base)."""
    return hermite_model(o["C"], o["mu"], o["var"], 4)


# ---------------------------------------------------------------------------------------------------------------
# objects from atlases
# ---------------------------------------------------------------------------------------------------------------
def _attach_pair(o, pair, l, h_index=1):
    if pair is not None:
        o["Gam"] = pair["Gam"][l]
        o["w2hat"] = pair["w2hat"][h_index, l]
    return o


class Atlas:
    """moment_atlas_np.py atlas (raw moments) with the big tensors loaded once."""

    def __init__(self, path, pairfile=None):
        self.pair = dict(np.load(pairfile)) if pairfile else None
        z = np.load(path)
        self.path = path
        self.files = z.files
        self.W = z["weights"].astype(np.float64)
        self.N = int(z["n_samples"])
        for k in ("pre_s", "post_s", "gate_p", "pre_M11", "pre_M21", "pre_M22", "post_M11"):
            setattr(self, k, z[k])
        self.gate_GG = z["gate_GG"] if "gate_GG" in z.files else None
        self.pre_M3 = z["pre_M3"]; self.post_M3 = z["post_M3"]
        self.pre_M211 = z["pre_M211"] if "pre_M211" in z.files else None

    def objects(self, l):
        mu, m2 = self.pre_s[0, l], self.pre_s[1, l]
        var = m2 - mu ** 2
        M11 = self.pre_M11[l].astype(np.float64)
        C = M11 - np.outer(mu, mu)
        T = ok.central3(self.pre_M3[l], M11, mu)
        mu_a = self.post_s[0, l]
        K3a = ok.central3(self.post_M3[l], self.post_M11[l].astype(np.float64), mu_a)
        mu1 = self.pre_s[0, l + 1]
        T1 = ok.central3(self.pre_M3[l + 1], self.pre_M11[l + 1].astype(np.float64), mu1)
        n = len(mu); idx = np.arange(n)
        o = dict(mu=mu, var=var, C=C, T=T, K3a=K3a, D21=T1[idx, idx, :].copy(), Phi=self.gate_p[l].astype(np.float64),
                 w=relu_w(mu, var))
        if self.gate_GG is not None:
            o["GG"] = self.gate_GG[l].astype(np.float64)
        if self.pre_M211 is not None:
            M21 = self.pre_M21[l].astype(np.float64); M3 = self.pre_M3[l]; M211 = self.pre_M211[l]
            Eu2uu = (M211 - np.einsum("k,ij->ijk", mu, M21) - np.einsum("j,ik->ijk", mu, M21) + np.einsum("j,k,i->ijk", mu, mu, m2)
                     - 2 * np.einsum("i,ijk->ijk", mu, M3 - np.einsum("k,ij->ijk", mu, M11) - np.einsum("j,ik->ijk", mu, M11)
                                     + np.einsum("i,j,k->ijk", mu, mu, mu))
                     + np.einsum("i,jk->ijk", mu ** 2, C))
            Cf = C.copy(); np.fill_diagonal(Cf, var)
            o["K4f"] = Eu2uu - np.einsum("i,jk->ijk", var, Cf) - 2 * np.einsum("ij,ik->ijk", Cf, Cf)
        return _attach_pair(o, self.pair, l)


class Atlas56:
    """atlas_k56.py atlas (central moments, with the kappa5 / kappa6 slices)."""

    def __init__(self, path, pairfile=None):
        self.pair = dict(np.load(pairfile)) if pairfile else None
        z = np.load(path)
        self.path = path
        self.W = z["weights"].astype(np.float64)
        self.N = int(z["n_samples"])
        self.d = {k: z[k] for k in z.files if k not in ("name",)}

    def objects(self, l):
        d = self.d
        n = d["pre_mean"].shape[1]; idx = np.arange(n)
        mu = d["pre_mean"][l]; C = d["pre_C"][l].copy(); var = np.diag(C).copy()
        cm = CentralMoments(var, C, d["pre_m3"][l], d["pre_m211"][l], d["pre_m221"][l], d["pre_m311"][l],
                            d["pre_m222"][l], d["pre_cm"][4, l])
        o = dict(mu=mu, var=var, C=C, T=d["pre_m3"][l], K4f=cm.cumulant("iijk"),
                 P=all_distinct(cm.cumulant("iijjk")), Q=all_distinct(cm.cumulant("iiijk")), S=all_distinct(cm.cumulant("iijjkk")),
                 Phi=d["gate_p"][l].astype(np.float64), w=relu_w(mu, var), K3a=d["post_m3"][l],
                 D21=d["pre_m3"][l + 1][idx, idx, :].copy())
        return _attach_pair(o, self.pair, l)


class CentralMoments:
    """central moments E[prod u] of multisets of the labels i, j, k as (n, n, n)-broadcast arrays, from the stored
    slices, and joint cumulants by the partition formula without singleton blocks (u is centred)."""

    BASE = {  # sorted count pattern -> (attribute, number of indices)
        (2,): ("var", 1), (1, 1): ("C", 2), (3,): ("D3", 1), (2, 1): ("Dm", 2), (1, 1, 1): ("m3", 3),
        (4,): ("cm4", 1), (3, 1): ("m31", 2), (2, 2): ("m22", 2), (2, 1, 1): ("m211", 3),
        (2, 2, 1): ("m221", 3), (3, 1, 1): ("m311", 3), (2, 2, 2): ("m222", 3)}

    def __init__(self, var, C, m3, m211, m221=None, m311=None, m222=None, cm4=None):
        n = len(var); idx = np.arange(n)
        self.n = n
        self.var, self.C, self.m3 = var, C, m3
        self.Dm = m3[idx, idx, :]; self.D3 = m3[idx, idx, idx]
        self.m211 = m211; self.m31 = m211[idx, idx, :]; self.m22 = m211[:, idx, idx]
        self.cm4 = m211[idx, idx, idx] if cm4 is None else cm4
        self.m221, self.m311, self.m222 = m221, m311, m222

    def moment(self, labels):
        if len(labels) == 1:
            return 0.0
        cnt = {s: labels.count(s) for s in "ijk" if labels.count(s)}
        order = sorted(cnt, key=lambda s: (-cnt[s], s))
        pat = tuple(cnt[s] for s in order)
        name, k = self.BASE[pat]
        A = getattr(self, name)
        # A's index a carries the label order[a]; for symmetric ties the order is irrelevant
        sub = "".join(order)
        outl = "".join(s for s in "ijk" if s in cnt)
        X = np.einsum(f"{sub}->{outl}", A)
        sh = [self.n if s in cnt else 1 for s in "ijk"]
        return X.reshape(sh)

    def cumulant(self, labels):
        import itertools
        pos = list(range(len(labels)))
        tot = 0.0
        for part in _partitions_no_singletons(pos):
            b = len(part)
            term = (-1) ** (b - 1) * _fact(b - 1)
            prod = 1.0
            for blk in part:
                prod = prod * self.moment("".join(labels[p] for p in blk))
            tot = tot + term * prod
        return np.broadcast_to(tot, (self.n,) * 3).copy()


def _fact(m):
    from math import factorial
    return factorial(m)


def _partitions_no_singletons(pos):
    if not pos:
        yield []
        return
    first, rest = pos[0], pos[1:]
    import itertools
    for r in range(1, len(rest) + 1):
        for comb in itertools.combinations(rest, r):
            remaining = [p for p in rest if p not in comb]
            for p in _partitions_no_singletons(remaining):
                yield [[first] + list(comb)] + p


# ---------------------------------------------------------------------------------------------------------------
# transport and evaluation
# ---------------------------------------------------------------------------------------------------------------
def transport(X, W):
    """D21(l+1)_{ab} = sum_{ijk} W_ia W_ja W_kb X_ijk  (BLAS: (n^2, n) @ (n, n), batched W^T @ ., then a diagonal contraction)."""
    n = X.shape[0]
    T1 = (X.reshape(n * n, n) @ W).reshape(n, n, n)          # [i, j, b]
    T2 = np.matmul(W.T[None, :, :], T1)                      # [i, a, b] = sum_j W_ja T1[i, j, b]
    return np.einsum("ia,iab->ab", W, T2)


# ---------------------------------------------------------------------------------------------------------------
# cheap transports of path-shaped terms (cost accounting, section 4 of REPORT.md)
# ---------------------------------------------------------------------------------------------------------------
# every term whose factors live on two of the three vertex pairs is a PATH  X_ijk = A_ik B_jk  (center k; the vertex
# weights are absorbed: A = diag(a) M diag(c), B = diag(b) N).  Its symmetrised all-distinct transport needs only
# n x n matmuls:  T(sym3 X) = (1/3)(R_center + R_end_i + R_end_j) - (i = j coincidence slice), with
#   G = W^T A, H = W^T B                                   2 matmuls
#   R_center = (G o H) W                                   1
#   R_end_i  = (A (W o H^T))^T W   ,  R_end_j likewise      2 + 2
#   coincidence slice Y_ik = A_ik B_ik: (W o W)^T Y W, (W o (Y W))^T W     2 + 2
# = 11 matmuls (11 units at n = 1024), and several path terms sharing one arm B cost the same 11 together.
PATH = {  # name: (A_ik, B_jk) builders on the objects (center k carries its weight inside A)
    "WICK": lambda o, q: (q["Phi"][:, None] * q["Co"] * q["w2"][None, :], q["Phi"][:, None] * q["Co"]),
    "B1": lambda o, q: (q["w2"][:, None] * q["Do"] * q["w2"][None, :], q["Phi"][:, None] * q["Co"]),
    "B2": lambda o, q: (q["Phi"][:, None] * q["Do"].T * q["w3"][None, :], q["Phi"][:, None] * q["Co"]),
    "B3": lambda o, q: (q["Phi"][:, None] * q["Co"] * (q["w5"] * q["D3"])[None, :], q["Phi"][:, None] * q["Co"]),
    "B5": lambda o, q: (q["w2"][:, None] * q["K22"] * q["w3"][None, :], q["Phi"][:, None] * q["Co"]),
    "S5c": lambda o, q: (q["Phi"][:, None] * q["V"].T * q["w4"][None, :], q["Phi"][:, None] * q["Co"]),
    "S4d": lambda o, q: (q["Phi"][:, None] * q["Do"].T * q["w4"][None, :], q["Phi"][:, None] * q["Do"].T),
    # the whole C-leaf class in one path transport: center i, arms Phi_k C_ik and Gamma_ij
    "LEAF": lambda o, q: (q["Phi"][:, None] * q["Co"], offdiag(o["Gam"]).T),
}


def path_objects(o):
    n = len(o["mu"]); idx = np.arange(n)
    Dm = o["T"][idx, idx, :]
    q = dict(Phi=o["Phi"], Co=offdiag(o["C"]), Do=offdiag(Dm), D3=np.diag(Dm).copy())
    for d in range(2, 7):
        q[f"w{d}"] = o["w"][d]
    if "K4f" in o:
        q["K22"] = offdiag(o["K4f"][:, idx, idx]); q["V"] = offdiag(o["K4f"][idx, idx, :])
    return q


def transport_path(A, B, W):
    """T(all_distinct(sym3[A_ik B_jk])) with 11 n x n matmuls (no n^3 object)."""
    G = W.T @ A; H = W.T @ B                       # [a, k]
    Rc = (G * H) @ W
    Ri = (A @ (W * H.T)).T @ W
    Rj = (B @ (W * G.T)).T @ W
    Y = A * B                                      # coincidence slice i = j:  X_iik = A_ik B_ik
    Yc = (W * W).T @ Y @ W
    P = Y @ W
    Ye = (W * P).T @ W                             # b on one of the coinciding ends (both ends give the same)
    return ((Rc + Ri + Rj) - (Yc + 2 * Ye)) / 3.0


def cost_check(path, l=10, pairfile=None):
    """dense n^4 transport of the path terms vs the 11-matmul formula, on one atlas layer."""
    A = Atlas(path, pairfile); o = A.objects(l); q = path_objects(o); W = A.W[l + 1]
    TT = terms(o)
    for name, fn in PATH.items():
        if name == "LEAF" and "Gam" not in o:
            continue
        Am, Bm = fn(o, q)
        t0 = time.time(); fast = transport_path(Am, Bm, W); t1 = time.time()
        if name == "WICK":
            dense_t = ok.wick_model(o["C"], o["Phi"], o["w"][2], np.zeros_like(o["T"]))
            dense = transport(dense_t, W)
            fast = 3 * fast     # the Wick term is the sum of the 3 labelled diagrams = 3 sym3[...]
        else:
            dense = transport(TT[name], W)
            # closure2 tensors are sym3[X] with X written on (i, j, k); PATH writes the same diagram with the
            # center on k, so the two must agree up to the term's symmetry factor
            pass
        t2 = time.time()
        print(f"  {name:6s} rel.diff 11-matmul path transport vs dense n^4 transport {rel(fast, dense):.1e}   (fast {t1 - t0:.3f}s)")


def corrected(e, en):
    return float(np.sqrt(max(e * e - en * en, 0.0)))


def fit_tensor(R, basis):
    X = np.stack([b.ravel() for b in basis], 1); y = R.ravel()
    G = X.T @ X; r = X.T @ y
    coef = np.linalg.solve(G, r)
    return coef, (X @ coef).reshape(R.shape)


def analyse_pair(A, B, layers, out=print, both=True):
    """noise-corrected cross-evaluation (models built from A, target from B; and B -> A when both)."""
    W = A.W; L, n, _ = W.shape
    assert np.array_equal(A.W, B.W)
    names_all = None
    rows = []
    for l in layers:
        t0 = time.time()
        oa, ob = A.objects(l), B.objects(l)
        Wn = W[l + 1]
        res = {}
        for tag, (o1, o2) in {"AB": (oa, ob), "BA": (ob, oa)}.items():
            if tag == "BA" and not both:
                continue
            TT = terms(o1)
            Kd = all_distinct(o1["K3a"])
            K3m = o1["K3a"] - Kd
            D21t = o2["D21"]
            en = rel(o1["D21"], o2["D21"]) / sqrt(2)
            Tr = {k: transport(v, Wn) for k, v in TT.items()}
            base_m = transport(K3m, Wn)
            Tw_gauss = transport(wick_gauss(o1), Wn)
            Kw_gate = ok.wick_model(o1["C"], o1["Phi"], o1["w"][2], o1["T"])     # Wick (gate Phi) + B0, oracle fit base
            Tw_gate = transport(Kw_gate, Wn)
            res[tag] = dict(o1=o1, TT=TT, Kd=Kd, K3m=K3m, D21t=D21t, en=en, Tr=Tr, base_m=base_m,
                            Tw_gauss=Tw_gauss, Kw_gate=Kw_gate, Tw_gate=Tw_gate)
        yield l, res
        del res


def summarise(l, res, Wn, emit):
    """compute the closure ladder for one layer and orientation-average it."""
    out = {}
    for tag, r in res.items():
        Tr, D21t, en = r["Tr"], r["D21t"], r["en"]
        o1 = r["o1"]; n = len(o1["mu"]); idx = np.arange(n); sg = np.sqrt(o1["var"])
        def nrm(X, k):
            return float(np.sqrt(np.mean(X ** 2)))
        sc = {"rho": nrm(offdiag(o1["C"] / np.outer(sg, sg)), 2),
              "k3_ijk": nrm(all_distinct(o1["T"] / np.einsum("i,j,k->ijk", sg, sg, sg)), 3),
              "D21_iik": nrm(offdiag(o1["T"][idx, idx, :] / np.outer(sg ** 2, sg)), 2),
              "D3_iii": nrm(o1["T"][idx, idx, idx] / sg ** 3, 1)}
        if "K4f" in o1:
            K = o1["K4f"] / np.einsum("i,j,k->ijk", sg ** 2, sg, sg)
            sc.update({"K22_iijj": nrm(offdiag(K[:, idx, idx]), 2), "K211_iijk": nrm(all_distinct(K), 3),
                       "K31_iiij": nrm(offdiag(K[idx, idx, :]), 2), "K4_iiii": nrm(K[idx, idx, idx], 1)})
        nD = float(np.sqrt(np.sum(D21t ** 2)))
        base = r["base_m"]

        def ev(T):
            return corrected(rel(T, D21t), en)
        o = {}
        o["noise"] = en
        o["wick"] = ev(base + r["Tw_gate"])
        cl_or = base + r["Tw_gauss"] + sum(c * Tr[t] for c, t in zip(ORACLE_COEF, FIRST) if t in Tr)
        o["cl1_oracle"] = ev(cl_or)
        cl1 = base + r["Tw_gauss"] + sum(c * Tr[t] for c, t in zip(LEG_COEF, FIRST) if t in Tr)
        o["cl1"] = ev(cl1)
        cl1p = cl1 + 3.0 * Tr["B7"]
        o["cl1+B7"] = ev(cl1p)
        sec = [t for t in SECOND_AVAIL if t in Tr]
        cl2 = cl1p + sum(TERM_INFO[t][2] * Tr[t] for t in sec)
        o["cl2"] = ev(cl2)
        new = [t for t in SECOND_NEW if t in Tr]
        if new:
            o["cl2+k56"] = ev(cl2 + sum(TERM_INFO[t][2] * Tr[t] for t in new))
        # group ablations: cl1+B7 plus one group
        groups = {}
        for t in sec + new:
            groups.setdefault(TERM_INFO[t][1], []).append(t)
        for g, ts in groups.items():
            o["+" + g] = ev(cl1p + sum(TERM_INFO[t][2] * Tr[t] for t in ts))
            o["-" + g] = ev(cl2 + (sum(TERM_INFO[t][2] * Tr[t] for t in new) if new else 0) - sum(TERM_INFO[t][2] * Tr[t] for t in ts))
        if "LEAF" in Tr:
            newl = [t for t in SECOND_NEW if t in Tr]
            lr1 = base + 6.0 * Tr["LEAF"] + Tr["B0"] + sum(TERM_INFO[t][2] * Tr[t] for t in NONLEAF1)
            o["LR1"] = ev(lr1 - Tr["TWOLEAF"])
            o["LR1e"] = ev(lr1 - Tr["TWOLEAFe"])
            lr2 = sum(TERM_INFO[t][2] * Tr[t] for t in NONLEAF2)
            o["LR2"] = ev(lr1 - Tr["TWOLEAF"] + lr2)
            o["LR2e"] = ev(lr1 - Tr["TWOLEAFe"] + lr2)
            if newl:
                o["LR2+k56"] = ev(lr1 - Tr["TWOLEAF"] + lr2 + sum(TERM_INFO[t][2] * Tr[t] for t in newl))
            # free refit of the leaf-resummed basis (theory: LEAF 6, TWOLEAF -1, B6 1.5, B7 3, G3t 1)
            TT0 = r["TT"]
            fixed2 = sum(TERM_INFO[t][2] * TT0[t] for t in NONLEAF2) + (sum(TERM_INFO[t][2] * TT0[t] for t in newl) if newl else 0)
            lb = ["LEAF", "TWOLEAF", "B6", "B7", "G3t"]
            cLR, _ = fit_tensor(r["Kd"] - TT0["B0"] - fixed2, [TT0[t] for t in lb])
            o["fitLR"] = ev(base + Tr["B0"] + lr2 + (sum(TERM_INFO[t][2] * Tr[t] for t in newl) if newl else 0)
                            + sum(c * Tr[t] for c, t in zip(cLR, lb)))
            o["coefLR"] = dict(zip(lb, cLR))
            if "TCL" in Tr:
                lrt1 = base + 6.0 * Tr["LEAF"] - Tr["TWOLEAF"] + Tr["TCL"] + 1.5 * Tr["B6"] + Tr["G3t"]
                o["LRT1"] = ev(lrt1)
                rest2 = [t for t in NONLEAF2 if t not in ("S3a1", "S4b", "S7a")]
                o["LRT2"] = ev(lrt1 + sum(TERM_INFO[t][2] * Tr[t] for t in rest2))
                lbt = ["LEAF", "TWOLEAF", "TCL", "B6", "G3t"]
                fixt = sum(TERM_INFO[t][2] * TT0[t] for t in rest2)
                cT, _ = fit_tensor(r["Kd"] - fixt, [TT0[t] for t in lbt])
                o["fitLRT"] = ev(base + sum(TERM_INFO[t][2] * Tr[t] for t in rest2) + sum(c * Tr[t] for c, t in zip(cT, lbt)))
                o["coefLRT"] = dict(zip(lbt, cT))
        # single-term ablations: cl1+B7 plus one term, cl2 minus one term
        o["add1"] = {t: ev(cl1p + TERM_INFO[t][2] * Tr[t]) for t in sec + new}
        o["drop1"] = {t: ev(cl2 + (sum(TERM_INFO[u][2] * Tr[u] for u in new) if new else 0) - TERM_INFO[t][2] * Tr[t]) for t in sec + new}
        # term sizes relative to ||D21||
        o["size"] = {t: TERM_INFO[t][2] * float(np.sqrt(np.sum(Tr[t] ** 2))) / nD for t in Tr}
        # fits in tensor space (built on the model atlas): oracle fit (Wick gate + B0..B6 free) and the refit with the
        # second-order terms and B7 at fixed coefficients
        TT, Kd, Kw = r["TT"], r["Kd"], r["Kw_gate"]
        basis = [TT[t] for t in FIRST]
        coef1, fit1 = fit_tensor(Kd - Kw, basis)
        o["fit1"] = ev(base + r["Tw_gate"] + sum(c * Tr[t] for c, t in zip(coef1, FIRST)))
        o["coef1"] = coef1
        fixed = 3.0 * TT["B7"] + sum(TERM_INFO[t][2] * TT[t] for t in sec if t != "S2a") + \
            (sum(TERM_INFO[t][2] * TT[t] for t in new) if new else 0)
        Tfixed = 3.0 * Tr["B7"] + sum(TERM_INFO[t][2] * Tr[t] for t in sec if t != "S2a") + \
            (sum(TERM_INFO[t][2] * Tr[t] for t in new) if new else 0)
        coef2, _ = fit_tensor(Kd - Kw - fixed, basis)
        o["fit2"] = ev(base + r["Tw_gate"] + Tfixed + sum(c * Tr[t] for c, t in zip(coef2, FIRST)))
        o["coef2"] = coef2
        # fit with B7 free too (8 free first-order coefficients) on top of the fixed second order
        fixed_noB7 = fixed - 3.0 * TT["B7"]
        coef3, _ = fit_tensor(Kd - Kw - fixed_noB7, basis + [TT["B7"]])
        o["coef3"] = coef3
        o["fit3"] = ev(base + r["Tw_gate"] + Tfixed - 3.0 * Tr["B7"] + sum(c * Tr[t] for c, t in zip(coef3, FIRST + ["B7"])))
        # all terms free (first + second order, 8 + len(sec) coefficients) in tensor space: does the span suffice?
        allb = FIRST + ["B7"] + [t for t in sec if t != "S2a"] + new
        coefa, _ = fit_tensor(Kd - Kw, [TT[t] for t in allb])
        o["fitall"] = ev(base + r["Tw_gate"] + sum(c * Tr[t] for c, t in zip(coefa, allb)))
        o["coefall"] = dict(zip(allb, coefa))
        # residual energy explained in tensor space
        def r2(M):
            return 1.0 - float(np.sum((Kd - M) ** 2)) / float(np.sum(Kd ** 2))
        o["tensR2_cl1+B7"] = r2(wick_gauss(r["o1"]) + sum(c * TT[t] for c, t in zip(LEG_COEF, FIRST)) + 3 * TT["B7"])
        o["tensR2_cl2"] = r2(wick_gauss(r["o1"]) + sum(c * TT[t] for c, t in zip(LEG_COEF, FIRST)) + 3 * TT["B7"]
                             + sum(TERM_INFO[t][2] * TT[t] for t in sec + new))
        o["scales"] = sc
        out[tag] = o
    # transported noise of every term (A vs B estimate of the same object)
    if "AB" in res and "BA" in res:
        nD = float(np.sqrt(np.sum(res["AB"]["D21t"] ** 2)))
        out["tnoise"] = {t: TERM_INFO[t][2] * float(np.sqrt(np.sum((res["AB"]["Tr"][t] - res["BA"]["Tr"][t]) ** 2) / 2)) / nD
                         for t in res["AB"]["Tr"]}
    return out


def fmt_layer(l, s):
    keys = [k for k in s["AB"] if not isinstance(s["AB"][k], (dict, np.ndarray))]
    avg = {k: np.mean([s[t][k] for t in ("AB", "BA") if t in s]) for k in keys}
    return avg


def main_pair(pa, pb, layers, outpath, k56=False, pairs=(None, None)):
    A = (Atlas56 if k56 else Atlas)(pa, pairs[0]); B = (Atlas56 if k56 else Atlas)(pb, pairs[1])
    W = A.W
    f = open(outpath, "w") if outpath else None

    def emit(s=""):
        print(s, flush=True)
        if f:
            f.write(s + "\n"); f.flush()
    emit(f"pair {pa} | {pb}: width {W.shape[1]}, depth {W.shape[0]}, N = {A.N} + {B.N}")
    emit("eps = noise-corrected relative rms error of D21(l+1), model from one atlas, target from the other, averaged over both orientations")
    allrows = []
    for l, res in analyse_pair(A, B, layers):
        s = summarise(l, res, W[l + 1], emit)
        res.clear()
        avg = fmt_layer(l, s)
        allrows.append((l, s, avg))
        if outpath:
            import pickle
            with open(outpath + ".pkl", "wb") as fh:
                pickle.dump([(ll, ss, aa) for ll, ss, aa in allrows], fh)
        cols = ["noise", "wick", "cl1_oracle", "cl1", "cl1+B7", "cl2"] + (["cl2+k56"] if "cl2+k56" in avg else []) + ["fit1", "fit2", "fit3", "fitall"]
        emit(f"L{l:02d} " + " ".join(f"{c}={avg[c]:.4f}" for c in cols))
        emit(f"     tensor R2: cl1+B7 {avg['tensR2_cl1+B7']:.4f}  cl2 {avg['tensR2_cl2']:.4f}")
        if "LR1" in avg:
            emit("     leaf-resummed: " + " ".join(f"{c}={avg[c]:.4f}" for c in ("LR1", "LR1e", "LR2", "LR2e", "LR2+k56", "fitLR", "LRT1", "LRT2", "fitLRT") if c in avg)
                 + "   coefLR AB: " + " ".join(f"{t}:{c:+.3f}" for t, c in s["AB"]["coefLR"].items()))
            if "coefLRT" in s["AB"]:
                emit("     coefLRT AB (theory LEAF 6, TWOLEAF -1, TCL 1, B6 1.5, G3t 1): "
                     + " ".join(f"{t}:{c:+.3f}" for t, c in s["AB"]["coefLRT"].items()))
        emit("     groups (+g on cl1+B7 | -g from cl2): " + "  ".join(f"{k[1:]}:{avg[k]:.4f}|{avg['-' + k[1:]]:.4f}" for k in avg if k.startswith("+")))
        for tag in ("AB", "BA"):
            emit(f"     coef {tag} fit1 (B0..B6): " + " ".join(f"{c:+.3f}" for c in s[tag]["coef1"]))
            emit(f"     coef {tag} fit2 (B0..B6 | 2nd order + B7 fixed): " + " ".join(f"{c:+.3f}" for c in s[tag]["coef2"]))
            emit(f"     coef {tag} fit3 (B0..B6, B7 | 2nd order fixed): " + " ".join(f"{c:+.3f}" for c in s[tag]["coef3"]))
        a1 = {t: np.mean([s[g]["add1"][t] for g in ("AB", "BA")]) for t in s["AB"]["add1"]}
        d1 = {t: np.mean([s[g]["drop1"][t] for g in ("AB", "BA")]) for t in s["AB"]["drop1"]}
        emit("     add1 (cl1+B7 + term): " + " ".join(f"{t}:{a1[t]:.4f}" for t in a1))
        emit("     drop1 (cl2 - term):   " + " ".join(f"{t}:{d1[t]:.4f}" for t in d1))
        emit("     coef AB fitall: " + " ".join(f"{t}:{c:+.2f}" for t, c in s["AB"]["coefall"].items()))
        emit("     scales (rms, normalised by sigmas): " + " ".join(f"{k}:{v:.2e}" for k, v in s["AB"]["scales"].items()))
        sz = s["AB"]["size"]
        emit("     size |coef*T(term)|/|D21|: " + " ".join(f"{t}:{sz[t]:.4f}" for t in sz))
        if "tnoise" in s:
            emit("     tnoise of term:            " + " ".join(f"{t}:{s['tnoise'][t]:.4f}" for t in s["tnoise"]))
    if f:
        f.close()
    return allrows


def parse_layers(s, L):
    if s is None:
        return list(range(L - 1))
    out = []
    for part in s.split(","):
        if "-" in part:
            a, b = part.split("-"); out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


if __name__ == "__main__":
    args = sys.argv[1:]
    opt = {}
    pos = []
    i = 0
    while i < len(args):
        if args[i] in ("--layers", "--out", "--pairfiles"):
            opt[args[i]] = args[i + 1]; i += 2
        else:
            pos.append(args[i]); i += 1
    if pos and pos[0] == "--selftest":
        import toy_diagrams
        toy_diagrams.check_closure2_identity()
    elif pos and pos[0] == "--costcheck":
        cost_check(pos[1], int(opt.get("--layers", 10)), opt.get("--pairfiles"))
    elif pos and pos[0] in ("--pair", "--pair56"):
        L = np.load(pos[1])["weights"].shape[0]
        pf = tuple(opt["--pairfiles"].split(",")) if "--pairfiles" in opt else (None, None)
        main_pair(pos[1], pos[2], parse_layers(opt.get("--layers"), L), opt.get("--out"), k56=(pos[0] == "--pair56"), pairs=pf)
