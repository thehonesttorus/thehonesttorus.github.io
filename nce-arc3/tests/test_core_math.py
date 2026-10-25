"""Unit tests for the framework's operators: kms, tyler, halmos, wold."""
import numpy as np
import pytest

from nce import halmos, kms, tyler, wold


# ----------------------------------------------------------------------------- kms

def _lib():
    L = kms.Library()
    L.add(kms.Primitive("a", "X", "X", 2.0))
    L.add(kms.Primitive("b", "X", "Y", 3.0))
    L.add(kms.Primitive("c", "Y", "X", 1.5))
    L.add(kms.Primitive("d", "Y", "Y", 6.0))
    return L


def test_critical_beta_root_and_parry_normalisation():
    L = _lib()
    st = kms.CriticalState.of(L)
    assert st is not None
    assert abs(kms.spectral_radius(L.matrix(st.beta_c)) - 1.0) < 1e-6
    for v in L.vertices:
        probs = st.transition_probabilities(L, v)
        assert abs(sum(probs.values()) - 1.0) < 1e-6          # Parry measure is a Markov chain
    # lower dimension => higher weight for single-edge words from the same vertex
    assert st.weight(L, ["c"]) > 0 and st.weight(L, ["d"]) > 0
    assert st.log_weight(L, ["a", "b"]) == pytest.approx(
        np.log(st.v[st.index["Y"]]) - np.log(st.v[st.index["X"]]) - st.beta_c * (np.log(2.0) + np.log(3.0)))


def test_adjoining_an_edge_raises_beta_c_amendment_E6():
    L = _lib()
    b0 = kms.critical_beta(L)
    L2 = L.copy()
    L2.add(kms.Primitive("e_new", "X", "Y", 4.0))
    b1 = kms.critical_beta(L2)
    assert b1 is not None and b0 is not None and b1 > b0


def test_gauge_cycle_has_no_critical_state():
    L = kms.Library()
    L.add(kms.Primitive("g", "X", "X", 1.0))
    L.add(kms.Primitive("h", "X", "X", 2.0))
    assert kms.critical_beta(L) is None


def test_promotion_is_frequency_gated():
    L = _lib()
    comp = ["a", "b"]
    rare = [("a", "b"), ("d",), ("c", "a"), ("d", "d"), ("c",), ("d", "c", "a")]
    frequent = [("a", "b")] * 12 + [("d",), ("c",)]
    d_rare = kms.promotion_decision(L, rare, comp, "ab", d_new=2.0)
    d_freq = kms.promotion_decision(L, frequent, comp, "ab", d_new=2.0)
    assert d_rare.beta_c_after > d_rare.beta_c_before
    assert d_freq.gain > d_rare.gain
    assert d_freq.promote
    assert kms.substitute(("c", "a", "b", "d"), comp, "ab") == ("c", "ab", "d")


# --------------------------------------------------------------------------- tyler

def test_tyler_recovers_gauge_of_isotropic_data():
    rng = np.random.default_rng(0)
    m, n = 5, 2000
    Z = rng.standard_normal((n, m))
    Z /= np.linalg.norm(Z, axis=1, keepdims=True)
    A = rng.standard_normal((m, m)) + 3 * np.eye(m)       # invertible gauge
    S = Z @ A.T
    res = tyler.tyler_fixed_point(S)
    assert res.converged and res.residual < 1e-6
    W = tyler.whiten(S, res.G)
    cov = W.T @ W / n
    assert np.allclose(cov, np.eye(m) / m, atol=0.02)      # balanced: I/m


def test_hn_filtration_detects_overfilled_subspace():
    rng = np.random.default_rng(1)
    m = 6
    # 70% of the mass in a 2-dimensional subspace violates nu(V) < dim V / m = 1/3
    V = rng.standard_normal((700, 2)) @ np.eye(2, m)
    W = rng.standard_normal((300, m))
    S = np.concatenate([V, W])
    flag = tyler.hn_filtration(S, iters=300)
    assert not flag.stable
    assert flag.blocks[0].shape[1] == 2
    assert flag.slopes[0] > flag.slopes[-1]
    assert abs(flag.masses[0] - 0.7) < 0.02
    assert tyler.is_subspace_stable(W)


# -------------------------------------------------------------------------- halmos

def test_tangent_tax_identity_and_commutator_identity():
    rng = np.random.default_rng(2)
    n = 7
    Pb = rng.standard_normal((n, 3))
    Qb = rng.standard_normal((n, 2))
    P, Q = halmos.projector(Pb), halmos.projector(Qb)
    y = Q @ rng.standard_normal(n)
    f, tax = halmos.min_norm_interpolant(Pb, Qb, y)
    assert np.allclose(Q @ f, y, atol=1e-8)
    th = halmos.principal_angles(Pb, Qb)
    # ||f||^2 = <y, (I + tan^2 Theta) y> computed in the Halmos frame
    QA, QB = halmos.orthonormal_basis(Pb), halmos.orthonormal_basis(Qb)
    U, s, Vt = np.linalg.svd(QA.T @ QB)
    coords = Vt @ (QB.T @ y)
    tan2 = (1 - s ** 2) / s ** 2
    assert f @ f == pytest.approx(float(np.sum((1 + tan2) * coords ** 2)), rel=1e-8)
    assert tax >= -1e-9
    # amendment E10: ||[P,Q]||_HS^2 = 2 Tr(P Q P^perp Q) = 1/2 sum sin^2 2theta
    lhs = halmos.commutator_hs(P, Q) ** 2
    assert lhs == pytest.approx(2 * np.trace(P @ Q @ (np.eye(n) - P) @ Q), rel=1e-8)
    assert lhs == pytest.approx(0.5 * np.sum(np.sin(2 * th) ** 2), rel=1e-8)


def test_kato_speed_limit_attained_by_commutator_direction():
    rng = np.random.default_rng(3)
    n = 6
    P = halmos.projector(rng.standard_normal((n, 2)))
    Q = halmos.projector(rng.standard_normal((n, 3)))
    G = P @ Q - Q @ P                                      # optimal direction
    rate, bound = halmos.kato_speed_limit(P, Q, G)
    assert rate == pytest.approx(bound, rel=1e-8)
    Pdot = G @ P - P @ G
    assert np.allclose(halmos.kato_generator(P, Pdot), G, atol=1e-8)


def test_disagreement_zero_for_identical_predictions():
    v = np.ones((3, 4))
    assert halmos.disagreement(v) == 0.0
    w = np.array([[1, 0, 0, 0], [0, 1, 0, 0]], dtype=float)
    assert halmos.disagreement(w) == 0.0                   # orthogonal rays commute as projections
    u = np.array([[1, 0, 0, 0], [1, 1, 0, 0]], dtype=float)
    assert halmos.disagreement(u) == pytest.approx(0.5)    # theta = 45 deg: 2 cos^2 sin^2 = 1/2


# ---------------------------------------------------------------------------- wold

def test_wold_filtration_newborn_and_absorb():
    F = wold.WoldFiltration(4)
    x = np.array([1.0, 0, 0, 0])
    r, nr = F.newborn(x)
    assert nr == 1.0 and F.absorb(x, grade=1) and F.rank == 1
    assert not F.absorb(2 * x, grade=1)                     # already explained
    y = np.array([1.0, 1.0, 0, 0])
    r, nr = F.newborn(y)
    assert np.allclose(r, [0, 1, 0, 0]) and nr == pytest.approx(1.0)
    assert F.newborn_fraction(y) == pytest.approx(0.5)
    assert F.absorb(y, grade=2) and F.rank == 2
    assert np.allclose(F.grade_projector(2), np.outer([0, 1, 0, 0], [0, 1, 0, 0]))
    tr = wold.transformed_residuals(np.array([1.0, 2.0]))
    assert np.allclose(tr["r2"], [1, 4]) and np.allclose(tr["rr"], [1, 2, 4])
    assert wold.hankel_rank_zero([np.zeros(3), np.zeros((2, 2))])
    assert not wold.hankel_rank_zero([np.zeros(3), np.array([0, 1e-3])])
