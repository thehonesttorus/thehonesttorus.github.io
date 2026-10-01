"""Rational-memory signal experiments (coordinator, 1 Oct 2026).

Model: first-order Wiener-chaos memory with mean gates (the Wick part of FC; coordinator note 1 §2), built on the exact
Gaussian closure. For a birth at layer s, neuron r (weight w2 = E delta(z_{s,r}) = p/s), the legs at target t are
  Z_s(t) = W_{s+1} D_{s+1} W_{s+2} ... D_{t-1} W_t      (rows r: layer-s neurons; cols: layer-t neurons)
  Y_s(t) = S_s D_s Z_s(t)                               (S_s = C(z_s), D = diag(P), P = P(z > 0))
and its contribution to D21(z_t) is w2 [ (y∘y) zᵀ + 2 (y∘z) yᵀ ].  All old atoms at a cut t0 share the future
propagator Q_{t0->t} = D_{t0} W_{t0+1} ... D_{t-1} W_t:  Y_s(t) = Y_s(t0) Q,  Z_s(t) = Z_s(t0) Q.
The readout that matters is the covariance correction E∘D21 with E_ab ≈ (p_a/s_a) P_b, priced by K_off(t).

X1  continuation (Hankel) spectrum of the old atoms' future responses: eigenvalues of the price-weighted response Gram.
X1b how much of the ACTUAL memory trajectory (all coefficients 1) the top-r response directions carry.
X1c per-target matrix rank of the old memory's readout (a lower bound for any atom carrier at that target).
X2  response-norm relabelling: fit R merged atoms at the cut, transported by the true future propagators, to the old
    memory's future readouts (Adam on the exact gradient), against pruning to the R strongest atoms.
"""
import argparse, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "notes", "fresh-slate", "bench"))
sys.path.insert(0, os.path.join(ROOT, "notes", "fresh-slate", "breakthrough", "region"))
import bench            # noqa: E402
import gclose           # noqa: E402

# region's measured K_off(l) (MLP 0, n = 1024, per unit ||dC(a_l)||_F^2); l = 15 does not feed a covariance
KOFF = np.array([2.6e-8, 7.4e-8, 8.6e-8, 1.8e-7, 2.3e-7, 2.7e-7, 4.4e-7, 5.4e-7, 6.0e-7, 6.2e-7, 7.0e-7, 7.6e-7,
                 6.3e-7, 4.7e-7, 2.9e-7, 0.0])


def closure_states(W):
    W = W.astype(np.float64)
    L, n, _ = W.shape
    st = []
    m = C = None
    for l in range(L):
        if l == 0:
            mu = np.zeros(n); S = W[0].T @ W[0]
        else:
            mu = m @ W[l]; S = W[l].T @ C @ W[l]
        m, C, P, p, s = gclose.closure_C(mu, S, "exact")
        st.append(dict(mu=mu, S=S, P=P, w2=p / s, e=p / s))
    return W, st


def legs_at_cut(W, st, t0, smax):
    """Stacked legs at target t0 for sources s = 0..smax: Y (Natoms, n), Z (Natoms, n), weights w (Natoms,), source id."""
    Ys, Zs, ws, sid = [], [], [], []
    for s in range(smax + 1):
        Z = W[s + 1].copy()
        for t in range(s + 1, t0):
            Z = (Z * st[t]["P"][None, :]) @ W[t + 1]
        Y = (st[s]["S"] * st[s]["P"][None, :]) @ Z
        Ys.append(Y); Zs.append(Z); ws.append(st[s]["w2"]); sid.append(np.full(Z.shape[0], s))
    return np.vstack(Ys), np.vstack(Zs), np.concatenate(ws), np.concatenate(sid)


def future_props(W, st, t0, tmax):
    """Q_{t0->t} for t = t0..tmax."""
    n = W.shape[1]
    Q = np.eye(n); out = [Q]
    for t in range(t0, tmax):
        Q = (Q * st[t]["P"][None, :]) @ W[t + 1]
        out.append(Q)
    return out


def readout(A, B, w, e, f):
    """Price-relevant readout E∘D21 = diag(e) [Σ w ((a∘a) bᵀ + 2 (a∘b) aᵀ)] diag(f) for atoms (rows of A, B)."""
    U = (A * A) * w[:, None]; V = (A * B) * w[:, None]
    D = U.T @ B + 2.0 * (V.T @ A)
    return (e[:, None] * D) * f[None, :]


def gram(A, B, w, e, f):
    """⟨R_p, R_q⟩_F for R_p = w_p diag(e)[(a∘a) bᵀ + 2 (a∘b) aᵀ]diag(f)."""
    U = e[None, :] * (A * A); V = e[None, :] * (A * B)
    Bf = B * f[None, :]; Af = A * f[None, :]
    G = (U @ U.T) * (Bf @ Bf.T) + 2 * (U @ V.T) * (Bf @ Af.T) + 2 * (V @ U.T) * (Af @ Bf.T) + 4 * (V @ V.T) * (Af @ Af.T)
    return G * np.outer(w, w)


def x1(W, st, t0, age_min, tmax, weights):
    L, n, _ = W.shape
    smax = t0 - age_min
    Y0, Z0, w, sid = legs_at_cut(W, st, t0, smax)
    Qs = future_props(W, st, t0, tmax)
    P = Y0.shape[0]
    G = np.zeros((P, P)); targets = []; tot = 0.0
    for k, t in enumerate(range(t0, tmax + 1)):
        om = weights[t]
        if om == 0: continue
        A = Y0 @ Qs[k]; B = Z0 @ Qs[k]
        e = st[t]["e"]; f = st[t]["P"]
        G += om * gram(A, B, w, e, f)
        M = readout(A, B, w, e, f)
        sv = np.linalg.svd(M, compute_uv=False); en = np.cumsum(sv ** 2) / np.sum(sv ** 2)
        targets.append(dict(t=t, fro2=float(np.sum(sv ** 2)), r90=int(np.searchsorted(en, 0.90) + 1),
                            r99=int(np.searchsorted(en, 0.99) + 1)))
        tot += om * float(np.sum(sv ** 2))
    lam, V = np.linalg.eigh(G); lam = lam[::-1]; V = V[:, ::-1]
    lam = np.maximum(lam, 0)
    cum = np.cumsum(lam) / lam.sum()
    one = np.ones(P)
    proj = (V.T @ one) ** 2 * lam
    cum1 = np.cumsum(proj) / proj.sum()
    def r_at(c, frac): return int(np.searchsorted(c, frac) + 1)
    out = dict(n=n, t0=t0, age_min=age_min, atoms=P, sources=int(smax + 1),
               hankel_r90=r_at(cum, 0.90), hankel_r99=r_at(cum, 0.99), hankel_r999=r_at(cum, 0.999),
               traj_r90=r_at(cum1, 0.90), traj_r99=r_at(cum1, 0.99), traj_r999=r_at(cum1, 0.999),
               traj_check=float(one @ G @ one / tot), targets=targets)
    return out, (Y0, Z0, w, Qs, G, lam, V)


def x2_fit(W, st, t0, tmax, weights, Y0, Z0, w, Qs, R, iters=1500, lr=0.02, seed=0, log_every=0, fitwin=99):
    """Fit R merged atoms (Yh, Zh at the cut; weights absorbed) to the old memory's priced future readouts."""
    n = W.shape[1]
    ks_all = [k for k, t in enumerate(range(t0, tmax + 1)) if weights[t] > 0]
    ts_all = [t0 + k for k in ks_all]
    targ_all = []; norm_all = 0.0
    for k, t in zip(ks_all, ts_all):
        A = Y0 @ Qs[k]; B = Z0 @ Qs[k]
        M = readout(A, B, w, st[t]["e"], st[t]["P"]); targ_all.append(M); norm_all += weights[t] * np.sum(M * M)
    sel = [i for i, k in enumerate(ks_all) if k <= fitwin]
    ks = [ks_all[i] for i in sel]; ts = [ts_all[i] for i in sel]; targ = [targ_all[i] for i in sel]
    norm = sum(weights[t] * np.sum(M * M) for M, t in zip(targ, ts))
    # init: the R atoms with the largest priced response, weights absorbed as w^(1/3)
    resp = np.zeros(Y0.shape[0])
    for k, t in zip(ks, ts):
        A = Y0 @ Qs[k]; B = Z0 @ Qs[k]
        e = st[t]["e"]; f = st[t]["P"]
        U = e[None, :] * (A * A); V = e[None, :] * (A * B); Bf = B * f[None, :]; Af = A * f[None, :]
        resp += weights[t] * (w * w) * (np.sum(U * U, 1) * np.sum(Bf * Bf, 1) + 4 * np.sum(U * V, 1) * np.sum(Bf * Af, 1)
                                         + 4 * np.sum(V * V, 1) * np.sum(Af * Af, 1))
    idx = np.argsort(resp)[::-1][:R]
    c = np.cbrt(w[idx])[:, None]
    Yh = Y0[idx] * c; Zh = Z0[idx] * c
    one = np.ones(R)
    def loss_grad(Yh, Zh):
        L_ = 0.0; gY = np.zeros_like(Yh); gZ = np.zeros_like(Zh)
        for M, k, t in zip(targ, ks, ts):
            Q = Qs[k]; e = st[t]["e"]; f = st[t]["P"]; om = weights[t]
            A = Yh @ Q; B = Zh @ Q
            E = readout(A, B, one, e, f) - M
            L_ += om * np.sum(E * E)
            GE = 2 * om * E
            Ge = (e[:, None] * GE) * f[None, :]          # pull back through diag(e) . diag(f)
            gA = 2 * A * (B @ Ge.T) + 2 * (A @ Ge.T) * B + 2 * ((A * B) @ Ge)
            gB = (A * A) @ Ge + 2 * (A @ Ge.T) * A
            gY += gA @ Q.T; gZ += gB @ Q.T
        return L_ / norm, gY / norm, gZ / norm
    prune_loss = loss_grad(Yh, Zh)[0]
    mY = np.zeros_like(Yh); vY = np.zeros_like(Yh); mZ = np.zeros_like(Zh); vZ = np.zeros_like(Zh)
    b1, b2, eps = 0.9, 0.999, 1e-12
    scale = np.sqrt(np.mean(Yh * Yh)) + 1e-12
    hist = []
    for it in range(1, iters + 1):
        Lv, gY, gZ = loss_grad(Yh, Zh)
        mY = b1 * mY + (1 - b1) * gY; vY = b2 * vY + (1 - b2) * gY * gY
        mZ = b1 * mZ + (1 - b1) * gZ; vZ = b2 * vZ + (1 - b2) * gZ * gZ
        step = lr * scale * (0.3 if it > 0.7 * iters else 1.0)
        Yh -= step * (mY / (1 - b1 ** it)) / (np.sqrt(vY / (1 - b2 ** it)) + eps)
        Zh -= step * (mZ / (1 - b1 ** it)) / (np.sqrt(vZ / (1 - b2 ** it)) + eps)
        if it % 100 == 0 or it == 1:
            hist.append((it, float(Lv)))
            if log_every: print(f"    it {it} rel loss {Lv:.4f}", flush=True)
    def eval_all(Yh, Zh):
        tot = 0.0; beyond = 0.0; nb = 0.0
        for M, k, t in zip(targ_all, ks_all, ts_all):
            Q = Qs[k]; E = readout(Yh @ Q, Zh @ Q, one, st[t]["e"], st[t]["P"]) - M
            tot += weights[t] * np.sum(E * E)
            if k > fitwin:
                beyond += weights[t] * np.sum(E * E); nb += weights[t] * np.sum(M * M)
        return tot / norm_all, (beyond / nb if nb > 0 else float('nan'))
    ev = eval_all(Yh, Zh)
    return dict(R=R, fitwin=fitwin, prune_rel=float(prune_loss), fit_rel=float(loss_grad(Yh, Zh)[0]),
                all_rel=float(ev[0]), beyond_rel=float(ev[1]), hist=hist)


def gradcheck():
    rng = np.random.default_rng(0)
    n, L = 12, 6
    W = (rng.standard_normal((L, n, n)) * np.sqrt(2.0 / n))
    W, st = closure_states(W)
    t0, tmax = 3, 5
    weights = np.ones(L)
    Y0, Z0, w, sid = legs_at_cut(W, st, t0, 1)
    Qs = future_props(W, st, t0, tmax)
    R = 3
    Yh = rng.standard_normal((R, n)); Zh = rng.standard_normal((R, n))
    # reuse internals by a tiny fit with iters=0 is awkward; replicate loss/grad here
    one = np.ones(R)
    targ = [readout(Y0 @ Qs[k], Z0 @ Qs[k], w, st[t0 + k]["e"], st[t0 + k]["P"]) for k in range(len(Qs))]
    def lossf(Yh, Zh):
        s = 0.0
        for k in range(len(Qs)):
            t = t0 + k
            E = readout(Yh @ Qs[k], Zh @ Qs[k], one, st[t]["e"], st[t]["P"]) - targ[k]
            s += np.sum(E * E)
        return s
    def grad(Yh, Zh):
        gY = np.zeros_like(Yh); gZ = np.zeros_like(Zh)
        for k in range(len(Qs)):
            t = t0 + k; Q = Qs[k]; e = st[t]["e"]; f = st[t]["P"]
            A = Yh @ Q; B = Zh @ Q
            E = readout(A, B, one, e, f) - targ[k]
            Ge = (e[:, None] * (2 * E)) * f[None, :]
            gA = 2 * A * (B @ Ge.T) + 2 * (A @ Ge.T) * B + 2 * ((A * B) @ Ge)
            gB = (A * A) @ Ge + 2 * (A @ Ge.T) * A
            gY += gA @ Q.T; gZ += gB @ Q.T
        return gY, gZ
    gY, gZ = grad(Yh, Zh)
    h = 1e-6; i, j = 1, 4
    Yp = Yh.copy(); Yp[i, j] += h; Ym = Yh.copy(); Ym[i, j] -= h
    num = (lossf(Yp, Zh) - lossf(Ym, Zh)) / (2 * h)
    Zp = Zh.copy(); Zp[i, j] += h; Zm = Zh.copy(); Zm[i, j] -= h
    numz = (lossf(Yh, Zp) - lossf(Yh, Zm)) / (2 * h)
    print("gradcheck Y:", num, gY[i, j], " Z:", numz, gZ[i, j])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="w256_d16"); ap.add_argument("--mlps", default="0")
    ap.add_argument("--t0", type=int, default=8); ap.add_argument("--age", type=int, default=3)
    ap.add_argument("--tmax", type=int, default=14)
    ap.add_argument("--fit", default="", help="comma list of R as fractions of n, e.g. 0.03125,0.0625,0.125")
    ap.add_argument("--iters", type=int, default=1500)
    ap.add_argument("--fitwin", type=int, default=99, help="fit only targets t0..t0+fitwin; evaluate on all")
    ap.add_argument("--gradcheck", action="store_true")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    if a.gradcheck:
        gradcheck(); sys.exit(0)
    S = bench.load_set(a.set)
    res = []
    for i in [int(x) for x in a.mlps.split(",")]:
        t_ = time.time()
        W, st = closure_states(bench.weights(S, i))
        r, (Y0, Z0, w, Qs, G, lam, V) = x1(W, st, a.t0, a.age, a.tmax, KOFF)
        r["mlp"] = i; r["sec_x1"] = round(time.time() - t_, 1)
        print(json.dumps({k: v for k, v in r.items() if k != "targets"}), flush=True)
        print("  per-target readout rank r90/r99:", [(d["t"], d["r90"], d["r99"]) for d in r["targets"]], flush=True)
        fits = []
        for fr in [float(x) for x in a.fit.split(",") if x]:
            R = max(1, int(round(fr * W.shape[1])))
            t_ = time.time()
            fr_ = x2_fit(W, st, a.t0, a.tmax, KOFF, Y0, Z0, w, Qs, R, iters=a.iters, fitwin=a.fitwin)
            fr_["sec"] = round(time.time() - t_, 1)
            print("  X2", json.dumps({k: v for k, v in fr_.items() if k != "hist"}), flush=True)
            fits.append(fr_)
        r["x2"] = fits
        res.append(r)
    if a.out:
        json.dump(res, open(a.out, "w"), indent=1)


def x3_als(W, st, t0, tmax, weights, Y0, Z0, w, Qs, R, sweeps=4, cg_iters=25, gd_iters=40, lr=0.02, fitwin=99,
           init="top"):
    """Cheaper fitter: alternate an exact linear least-squares solve for the third legs Zh (conjugate gradient on the
    normal equations; the readout is linear in Zh for fixed Yh) with a few Adam steps on Yh. Reports the loss on the
    fit window and on all future targets after each sweep, and the work in 'readout passes' (one pass = transport and
    read out R atoms at every fitted target)."""
    ks_all = [k for k, t in enumerate(range(t0, tmax + 1)) if weights[t] > 0]
    ts_all = [t0 + k for k in ks_all]
    targ_all = []; norm_all = 0.0
    for k, t in zip(ks_all, ts_all):
        M = readout(Y0 @ Qs[k], Z0 @ Qs[k], w, st[t]["e"], st[t]["P"]); targ_all.append(M)
        norm_all += weights[t] * np.sum(M * M)
    sel = [i for i, k in enumerate(ks_all) if k <= fitwin]
    ks = [ks_all[i] for i in sel]; ts = [ts_all[i] for i in sel]; targ = [targ_all[i] for i in sel]
    norm = sum(weights[t] * np.sum(M * M) for M, t in zip(targ, ts))
    one = np.ones(R)
    # init
    if init == "top":
        resp = np.zeros(Y0.shape[0])
        for k, t in zip(ks, ts):
            A = Y0 @ Qs[k]; B = Z0 @ Qs[k]; e = st[t]["e"]; f = st[t]["P"]
            U = e[None, :] * (A * A); V = e[None, :] * (A * B); Bf = B * f[None, :]; Af = A * f[None, :]
            resp += weights[t] * (w * w) * (np.sum(U * U, 1) * np.sum(Bf * Bf, 1) + 4 * np.sum(U * V, 1) * np.sum(Bf * Af, 1)
                                             + 4 * np.sum(V * V, 1) * np.sum(Af * Af, 1))
        idx = np.argsort(resp)[::-1][:R]; c = np.cbrt(w[idx])[:, None]
        Yh = Y0[idx] * c; Zh = Z0[idx] * c
    else:  # principal directions of the weighted y-cloud, Zh from the LS solve
        Ys = Y0 * np.sqrt(w)[:, None]
        _, s_, Vt = np.linalg.svd(Ys, full_matrices=False)
        Yh = Vt[:R] * (s_[:R, None] / np.sqrt(Y0.shape[0])) ** (2 / 3)
        Zh = np.zeros((R, Y0.shape[1]))
    passes = 0

    def fwd(Yh, Zh, k, t):
        Q = Qs[k]; return readout(Yh @ Q, Zh @ Q, one, st[t]["e"], st[t]["P"])

    def adjZ(Yh, X, k, t):      # gradient wrt Zh of <X, readout> with X already the residual
        Q = Qs[k]; e = st[t]["e"]; f = st[t]["P"]
        A = Yh @ Q; Ge = (e[:, None] * X) * f[None, :]
        gB = (A * A) @ Ge + 2 * (A @ Ge.T) * A
        return gB @ Q.T

    def loss_on(Yh, Zh, which):
        tot = 0.0
        for M, k, t in which:
            E = fwd(Yh, Zh, k, t) - M; tot += weights[t] * np.sum(E * E)
        return tot
    fit_list = list(zip(targ, ks, ts)); all_list = list(zip(targ_all, ks_all, ts_all))

    def solve_Z(Yh, Zh):
        nonlocal passes
        def Nop(Z):
            nonlocal passes
            out = np.zeros_like(Z)
            for M, k, t in fit_list:
                out += weights[t] * adjZ(Yh, fwd(Yh, Z, k, t), k, t)
            passes += 1
            return out
        rhs = np.zeros_like(Zh)
        for M, k, t in fit_list:
            rhs += weights[t] * adjZ(Yh, M, k, t)
        x = Zh.copy(); r = rhs - Nop(x); p = r.copy(); rs = np.sum(r * r)
        for _ in range(cg_iters):
            Ap = Nop(p); alpha = rs / max(np.sum(p * Ap), 1e-300)
            x += alpha * p; r -= alpha * Ap; rs_new = np.sum(r * r)
            if rs_new < 1e-20 * np.sum(rhs * rhs): break
            p = r + (rs_new / rs) * p; rs = rs_new
        return x

    def grad_Y(Yh, Zh):
        nonlocal passes
        gY = np.zeros_like(Yh)
        for M, k, t in fit_list:
            Q = Qs[k]; e = st[t]["e"]; f = st[t]["P"]
            A = Yh @ Q; B = Zh @ Q
            E = readout(A, B, one, e, f) - M
            Ge = (e[:, None] * (2 * weights[t] * E)) * f[None, :]
            gA = 2 * A * (B @ Ge.T) + 2 * (A @ Ge.T) * B + 2 * ((A * B) @ Ge)
            gY += gA @ Q.T
        passes += 1
        return gY
    hist = []
    mY = np.zeros_like(Yh); vY = np.zeros_like(Yh); b1, b2 = 0.9, 0.999; it_global = 0
    scale = np.sqrt(np.mean(Yh * Yh)) + 1e-12
    for sw in range(sweeps):
        Zh = solve_Z(Yh, Zh)
        hist.append(dict(sweep=sw, stage="Z", fit=loss_on(Yh, Zh, fit_list) / norm, all=loss_on(Yh, Zh, all_list) / norm_all,
                         passes=passes))
        for _ in range(gd_iters):
            it_global += 1
            g = grad_Y(Yh, Zh)
            mY = b1 * mY + (1 - b1) * g; vY = b2 * vY + (1 - b2) * g * g
            Yh -= lr * scale * (mY / (1 - b1 ** it_global)) / (np.sqrt(vY / (1 - b2 ** it_global)) + 1e-30)
        hist.append(dict(sweep=sw, stage="Y", fit=loss_on(Yh, Zh, fit_list) / norm, all=loss_on(Yh, Zh, all_list) / norm_all,
                         passes=passes))
    Zh = solve_Z(Yh, Zh)
    hist.append(dict(sweep=sweeps, stage="Z", fit=loss_on(Yh, Zh, fit_list) / norm, all=loss_on(Yh, Zh, all_list) / norm_all,
                     passes=passes))
    return dict(R=R, fitwin=fitwin, init=init, hist=hist)
