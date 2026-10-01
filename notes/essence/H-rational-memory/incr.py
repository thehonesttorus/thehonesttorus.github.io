"""Causal response-matched odometer (coordinator experiment).

Young sources (ages 1..AY at the current target) are exact. When a source turns AY+1 it is folded into a merged set of
R histories: fit R new histories (ALS: CG least squares for third legs, Adam steps for first legs, warm-started from
the previous merged set) to the priced readouts of (previous merged set + new source) over targets t..t+WIN, using the
true future propagators. Evaluate at every target: the priced readout error of the carried old content against the
exact old content (all atoms of age > AY), and its share of the total priced D21 energy (young + old).
"""
import json, sys, time
import numpy as np
import rmem, bench

def legs_source(W, st, s, t):
    Z = W[s + 1].copy()
    for u in range(s + 1, t):
        Z = (Z * st[u]["P"][None, :]) @ W[u + 1]
    Y = (st[s]["S"] * st[s]["P"][None, :]) @ Z
    return Y, Z

def step(Y, Z, st, t, W):  # transport legs from target t to t+1
    G = (np.ones(1) * st[t]["P"])[None, :]
    return (Y * G) @ W[t + 1], (Z * G) @ W[t + 1]

def run(set_name, mlp, AY=2, R_frac=0.25, WIN=2, sweeps=2, cg_iters=15, gd_iters=20, tstart=None, log=print):
    S = bench.load_set(set_name)
    W, st = rmem.closure_states(bench.weights(S, mlp))
    L, n, _ = W.shape
    R = max(1, int(round(R_frac * n)))
    one = np.ones(R)
    merged = None          # (Yh, Zh) at the current target
    rows = []; work_passes = 0
    for t in range(AY + 1, L - 1):           # targets with a covariance price (t <= 14)
        if rmem.KOFF[t] == 0: continue
        # exact old content at t (evaluation only): all sources with age t - s > AY
        Mold = np.zeros((n, n)); Myoung = np.zeros((n, n))
        for s in range(0, t):
            Y, Z = legs_source(W, st, s, t)
            Mx = rmem.readout(Y, Z, st[s]["w2"], st[t]["e"], st[t]["P"])
            if t - s > AY: Mold += Mx
            else: Myoung += Mx
        # fold in the source that turns AY+1 now: s_new = t - AY - 1
        s_new = t - AY - 1
        Yn, Zn = legs_source(W, st, s_new, t)
        wn = st[s_new]["w2"]
        if merged is None:
            Yin, Zin, win_ = Yn, Zn, wn
        else:
            Yin = np.vstack([merged[0], Yn]); Zin = np.vstack([merged[1], Zn]); win_ = np.concatenate([np.ones(merged[0].shape[0]), wn])
        if Yin.shape[0] <= R:
            c = np.cbrt(win_)[:, None]; merged = (Yin * c, Zin * c)
        else:
            # targets for the fit
            tmax = L - 2
            targets = [u for u in range(t, min(t + WIN, tmax) + 1) if rmem.KOFF[u] > 0]
            Qs = rmem.future_props(W, st, t, max(targets))
            # warm start: previous merged set if it exists, else strongest atoms
            if merged is not None and merged[0].shape[0] == R:
                Yh0, Zh0 = merged
            else:
                c = np.cbrt(win_)[:, None]
                nrm = np.sum((Yin * c) ** 2, 1) * np.sum((Zin * c) ** 2, 1)
                idx = np.argsort(nrm)[::-1][:R]; Yh0, Zh0 = (Yin * c)[idx], (Zin * c)[idx]
            Yh, Zh, passes = als_fit(W, st, t, targets, Qs, Yin, Zin, win_, Yh0.copy(), Zh0.copy(), sweeps, cg_iters, gd_iters)
            work_passes += passes
            merged = (Yh, Zh)
        # evaluate carried old content at t
        Mc = rmem.readout(merged[0], merged[1], np.ones(merged[0].shape[0]), st[t]["e"], st[t]["P"])
        E = Mc - Mold
        o2 = np.sum(Mold ** 2); tot2 = np.sum((Mold + Myoung) ** 2); e2 = np.sum(E ** 2)
        rows.append(dict(t=t, old_rel=float(e2 / o2), tot_rel=float(e2 / tot2), old_share=float(o2 / tot2), K=float(rmem.KOFF[t])))
        log(f"t={t:2d} old-rel {e2/o2:.4f} total-rel {e2/tot2:.4f} (old share {o2/tot2:.2f}) passes so far {work_passes}")
        # transport merged set to t+1 for the next step
        if t + 1 <= L - 1:
            merged = step(merged[0], merged[1], st, t, W)
    priced = sum(r["K"] * r["tot_rel"] for r in rows) / sum(r["K"] for r in rows)
    return dict(set=set_name, mlp=mlp, AY=AY, R=R, WIN=WIN, sweeps=sweeps, rows=rows, priced_total_rel=priced, passes=work_passes)

def als_fit(W, st, t, targets, Qs, Yin, Zin, win_, Yh, Zh, sweeps, cg_iters, gd_iters, lr=0.02):
    R = Yh.shape[0]; one = np.ones(R)
    ks = [u - t for u in targets]
    Ms = [rmem.readout(Yin @ Qs[k], Zin @ Qs[k], win_, st[u]["e"], st[u]["P"]) for k, u in zip(ks, targets)]
    wts = [rmem.KOFF[u] for u in targets]
    passes = 0
    def fwd(Yh, Zh, k, u): Q = Qs[k]; return rmem.readout(Yh @ Q, Zh @ Q, one, st[u]["e"], st[u]["P"])
    def adjZ(Yh, X, k, u):
        Q = Qs[k]; e = st[u]["e"]; f = st[u]["P"]; A = Yh @ Q; Ge = (e[:, None] * X) * f[None, :]
        return ((A * A) @ Ge + 2 * (A @ Ge.T) * A) @ Q.T
    def solveZ(Yh, Zh):
        nonlocal passes
        def Nop(Zx):
            nonlocal passes
            out = np.zeros_like(Zx)
            for M, k, u, om in zip(Ms, ks, targets, wts): out += om * adjZ(Yh, fwd(Yh, Zx, k, u), k, u)
            passes += 1; return out
        rhs = np.zeros_like(Zh)
        for M, k, u, om in zip(Ms, ks, targets, wts): rhs += om * adjZ(Yh, M, k, u)
        x = Zh.copy(); r = rhs - Nop(x); p = r.copy(); rs = np.sum(r * r)
        for _ in range(cg_iters):
            Ap = Nop(p); a = rs / max(np.sum(p * Ap), 1e-300); x += a * p; r -= a * Ap; rn = np.sum(r * r)
            if rn < 1e-20 * np.sum(rhs * rhs): break
            p = r + (rn / rs) * p; rs = rn
        return x
    def gradY(Yh, Zh):
        nonlocal passes
        g = np.zeros_like(Yh)
        for M, k, u, om in zip(Ms, ks, targets, wts):
            Q = Qs[k]; e = st[u]["e"]; f = st[u]["P"]; A = Yh @ Q; B = Zh @ Q
            E = rmem.readout(A, B, one, e, f) - M; Ge = (e[:, None] * (2 * om * E)) * f[None, :]
            gA = 2 * A * (B @ Ge.T) + 2 * (A @ Ge.T) * B + 2 * ((A * B) @ Ge); g += gA @ Q.T
        passes += 1; return g
    mY = np.zeros_like(Yh); vY = np.zeros_like(Yh); b1, b2 = 0.9, 0.999; itg = 0
    scale = np.sqrt(np.mean(Yh * Yh)) + 1e-12
    for sw in range(sweeps):
        Zh = solveZ(Yh, Zh)
        for _ in range(gd_iters):
            itg += 1; g = gradY(Yh, Zh)
            mY = b1 * mY + (1 - b1) * g; vY = b2 * vY + (1 - b2) * g * g
            Yh -= lr * scale * (mY / (1 - b1 ** itg)) / (np.sqrt(vY / (1 - b2 ** itg)) + 1e-30)
    Zh = solveZ(Yh, Zh)
    return Yh, Zh, passes

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="w256_d16"); ap.add_argument("--mlp", type=int, default=0)
    ap.add_argument("--AY", type=int, default=2); ap.add_argument("--R", type=float, default=0.25)
    ap.add_argument("--win", type=int, default=2); ap.add_argument("--sweeps", type=int, default=2)
    ap.add_argument("--cg", type=int, default=15); ap.add_argument("--gd", type=int, default=20)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    t_ = time.time()
    res = run(a.set, a.mlp, a.AY, a.R, a.win, a.sweeps, a.cg, a.gd, log=lambda s: print(s, flush=True))
    res["sec"] = round(time.time() - t_, 1)
    print(json.dumps({k: v for k, v in res.items() if k != "rows"}), flush=True)
    if a.out: json.dump(res, open(a.out, "w"), indent=1)
