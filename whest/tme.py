"""The transport-memory estimator, implemented line by line from "The system, derived" (stage-14 note, Section 6).

State carried from layer to layer: the Gaussian part (m, C) of the pre-activations, the coherent pair-kurtosis
matrix K^z, and third-cumulant history registers. A register (P, Q) represents sum_k sym(P_:k x P_:k x Q_:k).
Each birth layer b contributes a cherry register (T A_b C_b, T diag a2_b) and a coincidence register (T, T Delta_b^T),
which share the transport T = T_{l<-b}; so per birth we carry T, TX = T A_b C_b and TD = T Delta_b^T.

Conventions: W[l] has shape (n_out, n_in) and z_l = W[l] @ h_{l-1}; the note's W^T is our W. Layers are 0-based.
Every coefficient is a closed-form Gaussian integral (Prop. jets, Cor. gauss, Cor. corr); no fitted constant,
no inverse, no amplitude gauge.
"""
import numpy as np
from scipy.special import ndtr

SQ2PI = np.sqrt(2 * np.pi)


def tokens(m, var):
    """Wall jets and exact moments of relu(z), z ~ N(m, var), per neuron (Prop. jets, Cor. gauss)."""
    sigma = np.sqrt(np.clip(var, 1e-300, None)); alpha = m / sigma
    Phi = ndtr(alpha); ph = np.exp(-0.5 * alpha * alpha) / SQ2PI
    G = sigma * (alpha * Phi + ph)
    a1 = Phi; a2 = ph / sigma; a3 = -alpha * ph / sigma ** 2; a4 = (alpha * alpha - 1) * ph / sigma ** 3
    G2 = sigma ** 2 * ((1 + alpha * alpha) * Phi + alpha * ph)
    M3 = sigma ** 3 * ((alpha ** 3 + 3 * alpha) * Phi + (alpha * alpha + 2) * ph)
    M4 = sigma ** 4 * ((alpha ** 4 + 6 * alpha * alpha + 3) * Phi + (alpha ** 3 + 5 * alpha) * ph)
    k3 = M3 - 3 * G2 * G + 2 * G ** 3
    k4 = M4 - 4 * M3 * G - 3 * G2 ** 2 + 12 * G2 * G ** 2 - 6 * G ** 4
    F1 = 2 * G * (1 - Phi); F2 = 2 * (a1 - G * a2)
    return dict(sigma=sigma, alpha=alpha, a1=a1, a2=a2, a3=a3, a4=a4, G=G, G2=G2, k3=k3, k4=k4, F1=F1, F2=F2)


def coincidence_rows(t, C, Kz, M=None, k3z=None):
    """Delta with rows delta_a: the births of Prop. coinc, the (1/4) F2^a a2^b K^z_ab term of (L4), and the exact
    first-order transfer of the coincident class of kappa3(z). The register transport A^{x3} gives a1_a^2 a1_b on an
    (a,a,b) entry, but two derivatives on the same variable give E[1{z_a>0}] = a1_a, so the exact coefficient is
    (1/2) F2^a a1^b = (a1 - G a2)_a a1_b (the 3 kappa_aab d_a^2 d_b / 6 term of the Edgeworth operator on
    Cov(X_a, u_b)); the difference is added here from the (2,1) readout M = kappa3(z)_aab. On the diagonal the exact
    coefficient of kappa3(z_a) in kappa3(u_a) is a1 - G a2 + a3 (G^2 - G2/2) in place of a1^3."""
    a1, a2, a3, F1, F2, k3, G, G2 = t["a1"], t["a2"], t["a3"], t["F1"], t["F2"], t["k3"], t["G"], t["G2"]
    Cd = np.diag(C)
    D = (a1[None, :] * C * (F1 - 2 * a1 * a2 * Cd)[:, None]
         + a2[None, :] * C * C * (0.5 * F2 - a1 * a1)[:, None]
         + 0.25 * F2[:, None] * a2[None, :] * Kz)
    diag = (k3 - 3 * a1 * a1 * a2 * Cd * Cd) / 3.0
    if M is not None:
        D = D + (0.5 * F2 - a1 * a1)[:, None] * a1[None, :] * M
        c_aaa = a1 - G * a2 + a3 * (G * G - 0.5 * G2)
        diag = diag + (c_aaa - a1 ** 3) * k3z / 3.0
    np.fill_diagonal(D, diag)
    return D


def readouts(regs):
    """k3_a = 3 sum_reg sum_k P_ak^2 Q_ak and M_ab = sum_reg [(P o P) Q^T + 2 (P o Q) P^T]_ab."""
    k3 = 0.0; M = 0.0
    for r in regs:
        T, TX, TD, a2b = r["T"], r["TX"], r["TD"], r["a2"]
        Q = T * a2b[None, :]                                           # cherry Q = T diag(a2_b)
        k3 = k3 + 3.0 * (np.einsum("ak,ak,ak->a", TX, TX, Q) + np.einsum("ak,ak,ak->a", T, T, TD))
        M = M + (TX * TX) @ Q.T + 2.0 * (TX * Q) @ TX.T + (T * T) @ TD.T + 2.0 * (T * TD) @ T.T
    return k3, M


def tme_chain(W, record=None, nreg=None, coinc_transfer=True):
    """Returns the activation means of every layer, shape (L, n). `nreg`: keep only the nreg youngest births
    (None = all; the derived system keeps all)."""
    W = [np.ascontiguousarray(Wl, dtype=np.float64) for Wl in W]
    L = len(W); n = W[0].shape[0]
    m = np.zeros(n); C = W[0] @ W[0].T; Kz = np.zeros((n, n)); regs = []
    out = np.empty((L, n))
    for l in range(L):
        t = tokens(m, np.diag(C))
        a1, a2, a3, a4, G, G2 = t["a1"], t["a2"], t["a3"], t["a4"], t["G"], t["G2"]
        if regs:
            k3, M = readouts(regs)
        else:
            k3 = np.zeros(n); M = np.zeros((n, n))
        Kd = np.diag(Kz).copy()
        Eu = G + k3 * a3 / 6.0 + Kd * a4 / 24.0                                       # (L1)
        out[l] = Eu
        if record is not None:
            record[l] = dict(m=m.copy(), var=np.diag(C).copy(), k3=k3.copy(), Kd=Kd, Eu=Eu.copy(), alpha=t["alpha"].copy())
        if l == L - 1:
            break
        # (L2)-(L3): activation covariance
        cov = (np.outer(a1, a1) * C + 0.5 * np.outer(a2, a2) * C * C
               + 0.5 * (M * np.outer(a2, a1) + M.T * np.outer(a1, a2)) + 0.25 * np.outer(a2, a2) * Kz)
        varu = G2 - G * G + k3 * (a2 - G * a3) / 3.0 + Kd * (a3 - G * a4) / 12.0
        np.fill_diagonal(cov, varu)
        # (L5): pair kurtosis of the activations and its coherent transport
        F1, F2 = t["F1"], t["F2"]
        Ku = np.outer(F1, F1) * C + (0.5 * np.outer(F2, F2) - 2.0 * np.outer(a1 * a1, a1 * a1)) * C * C + 0.25 * np.outer(F2, F2) * Kz
        np.fill_diagonal(Ku, t["k4"])
        Wn = W[l + 1]; W2 = Wn * Wn
        Kz_new = W2 @ Ku @ W2.T
        dK = 3.0 * np.diag(Kz_new) - 2.0 * (W2 * W2) @ np.diag(Ku)
        np.fill_diagonal(Kz_new, dK)
        # registers: transport the old ones by W A_l, then add this layer's births
        A = a1
        for r in regs:
            r["T"] = Wn @ (A[:, None] * r["T"]); r["TX"] = Wn @ (A[:, None] * r["TX"]); r["TD"] = Wn @ (A[:, None] * r["TD"])
        Delta = coincidence_rows(t, C, Kz, M if coinc_transfer else None, k3)
        regs.append(dict(T=Wn.copy(), TX=Wn @ (A[:, None] * C), TD=Wn @ Delta.T, a2=a2.copy(), born=l))
        if nreg is not None and len(regs) > nreg:
            regs = regs[-nreg:]
        # Gaussian part
        m = Wn @ Eu; C = Wn @ cov @ Wn.T; Kz = Kz_new
    return out
