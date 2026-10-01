"""Prong 2 of the research programme (see ../research-program.md): targeted experiments
that read an actual ReLU MLP into the abstract barycentre-arrow objects of
../conditional-arrow-algebra.md and ../simplicial-complex-as-decomposition.md and test
the theory's dichotomies against that realization.

Pure numpy.  No torch, no scipy.

Run:
    python3 mlp_bridge.py                       # trained numpy stand-in (2-d inputs, widths 24,24,24)
    python3 mlp_bridge.py --mode random         # untrained stand-in with the same architecture
    python3 mlp_bridge.py --weights net.npz     # your own network (format below)
    python3 mlp_bridge.py --weights net.npz --inputs X.npy --restrict 12 --input-face pattern

Weight file format (.npz): keys W1,b1,W2,b2,...,WL,bL with W_l of shape (n_l, n_{l-1}),
b_l of shape (n_l,), forward pass h_l = relu(W_l h_{l-1} + b_l) for l < L and logits
W_L h_{L-1} + b_L (no ReLU on the last layer).  Transposed weights (n_{l-1}, n_l) are
detected and transposed.  --inputs X.npy is an (N, n_0) array of inputs; without it the
inputs are N(0, I) samples.

THE DICTIONARY (the one modelling decision, stated so that a negative result can be charged
to the dictionary and not to the theory; see research-program.md, rule P2.3):

  object (face) at hidden layer l  := nonzero pattern of the ReLU output of layer l, as a
                                      subset of the units; its dimension d = number of units
  object at layer 0                := one trivial face '*' (or, with --input-face pattern,
                                      the nonzero pattern of the input vector)
  object at the output layer       := the argmax class (a vertex of the output simplex), or
                                      the sign of a single logit
  arrow  sigma -> tau              := an observed transition: the ReLU that produced sigma,
                                      then W_{l+1}, landing in tau at layer l+1
  multiplicity of an arrow         := empirical count
  state on layer 0                 := the input distribution
  cocycle F(gamma)                 := -log P(gamma | s(gamma)), the empirical kernel
                                      (the Doob-normalised walk of Note 1, Section 5.3)
  restriction (--restrict k)       := the face is read on a sub-frame of k units only
                                      (pinching to a subalgebra of the frame, Note 2, Section 2)

What is tested (each experiment targets one exact statement of the theory):

  E1  the quiver: faces, arrows, dimension jumps |d(tau)-d(sigma)| > 1, heights (number of
      distinct histories arriving at a face), coverage.
  E2  sufficiency of the barycentre as conditional independence: I(next ; past | present)
      per layer and per lag, against a permutation null (Note 1, Section 6).
  E3  the cohomology class of F: weighted least-squares projection of F on coboundaries
      (residual = non-product part of log P), the within-endpoint spread of the accumulated
      action F(mu) (Note 1, Section 10, item 2) and the centrality defect (Note 1, Section 6:
      face sufficient <=> F coboundary <=> P central/tail-invariant).
  E4  the composition fork (Note 2, Section 3): is the outgoing kernel of a face close to that
      of an overlapping face?  correlation between the overlap cosine |s&s'|/sqrt(|s||s'|)
      and the fidelity of the two next-face laws, against a null that re-pairs patterns and
      kernels at random; and whether the arrows transport the overlap (corr of G before and
      the expected G after one step).
  E6  Betti numbers b0, b1 of the 2-skeleton of each layer nerve (patterns only, Note 2,
      Section 7), when the width allows.
"""
import argparse, json, sys, time
import numpy as np

# ----------------------------------------------------------------------------- utilities

def relu(x):
    return np.maximum(x, 0.0)


def entropy_from_keys(keys):
    """Plug-in entropy (nats) of the empirical law of an integer-key array."""
    _, c = np.unique(keys, return_counts=True)
    p = c / c.sum()
    return float(-(p * np.log(p)).sum())


def combine(*cols):
    """Mixed-radix combination of several non-negative integer id arrays into one key array."""
    cols = [np.asarray(c, dtype=np.int64) for c in cols]
    key = np.zeros_like(cols[0])
    for c in cols:
        m = int(c.max()) + 1
        if key.max() > (np.iinfo(np.int64).max // max(m, 1)) - 1:
            # overflow: fall back to row-unique relabelling
            key = np.unique(np.stack([key, c], 1), axis=0, return_inverse=True)[1].reshape(-1).astype(np.int64)
        else:
            key = key * m + c
    return key


def cmi(A, B, C):
    """I(A;B|C) in nats, plug-in."""
    return (entropy_from_keys(combine(A, C)) + entropy_from_keys(combine(B, C))
            - entropy_from_keys(combine(A, B, C)) - entropy_from_keys(C))


def permute_within(B, C, rng):
    """Permute B within the strata of C (destroys A-B dependence given C, keeps marginals)."""
    idx = np.argsort(C, kind="stable")
    perm = np.lexsort((rng.random(len(C)), C))
    out = np.empty_like(B)
    out[idx] = B[perm]
    return out


def ids_of_rows(bool_rows):
    """Map boolean rows to dense integer ids; return ids and the distinct rows (as bool matrix)."""
    packed = np.packbits(bool_rows, axis=1)
    uniq, inv = np.unique(packed, axis=0, return_inverse=True)
    inv = inv.reshape(-1)
    # recover the distinct boolean rows
    first = np.zeros(len(uniq), dtype=np.int64)
    first[inv[::-1]] = np.arange(len(inv))[::-1]
    return inv.astype(np.int64), bool_rows[first]


# ----------------------------------------------------------------------------- networks

def load_weights(path):
    z = np.load(path)
    Ws, bs = [], []
    l = 1
    while f"W{l}" in z:
        W = np.asarray(z[f"W{l}"], dtype=np.float64)
        b = np.asarray(z[f"b{l}"], dtype=np.float64).reshape(-1) if f"b{l}" in z else None
        Ws.append(W)
        bs.append(b)
        l += 1
    if not Ws:
        sys.exit("no W1 in weight file")
    # orient: W_l must map n_{l-1} -> n_l
    n_in = Ws[0].shape[1] if bs[0] is None or Ws[0].shape[0] == bs[0].shape[0] else Ws[0].shape[0]
    fixed = []
    prev = n_in
    for W, b in zip(Ws, bs):
        if W.shape[1] != prev and W.shape[0] == prev:
            W = W.T
        if W.shape[1] != prev:
            sys.exit(f"cannot orient weight of shape {W.shape} after width {prev}")
        if b is None:
            b = np.zeros(W.shape[0])
        fixed.append((W, b))
        prev = W.shape[0]
    return [w for w, _ in fixed], [b for _, b in fixed]


def random_net(widths, rng):
    Ws, bs = [], []
    for a, b in zip(widths[:-1], widths[1:]):
        Ws.append(rng.standard_normal((b, a)) * np.sqrt(2.0 / a))
        bs.append(rng.standard_normal(b) * 0.1)
    return Ws, bs


def teacher_label(X):
    f = np.sin(2.2 * X[:, 0]) + np.sin(2.2 * X[:, 1]) + 0.5 * X[:, 0] * X[:, 1]
    return (f > 0).astype(np.float64)


def train_net(widths, rng, steps=4000, batch=256, lr=3e-3, n_train=20000):
    """Minimal numpy MLP trainer (Adam, logistic loss, one logit) on a synthetic 2-d task."""
    assert widths[0] == 2 and widths[-1] == 1
    Ws, bs = random_net(widths, rng)
    X = rng.standard_normal((n_train, 2)) * 1.3
    y = teacher_label(X)
    m = [np.zeros_like(W) for W in Ws]; v = [np.zeros_like(W) for W in Ws]
    mb = [np.zeros_like(b) for b in bs]; vb = [np.zeros_like(b) for b in bs]
    b1, b2, eps = 0.9, 0.999, 1e-8
    for t in range(1, steps + 1):
        idx = rng.integers(0, n_train, batch)
        x, yy = X[idx], y[idx]
        hs, pre = [x], []
        h = x
        for l, (W, b) in enumerate(zip(Ws, bs)):
            z = h @ W.T + b
            pre.append(z)
            h = relu(z) if l < len(Ws) - 1 else z
            hs.append(h)
        logit = hs[-1][:, 0]
        p = 1.0 / (1.0 + np.exp(-logit))
        g = ((p - yy) / batch)[:, None]
        gW, gb = [None] * len(Ws), [None] * len(Ws)
        for l in range(len(Ws) - 1, -1, -1):
            gW[l] = g.T @ hs[l]
            gb[l] = g.sum(0)
            if l > 0:
                g = (g @ Ws[l]) * (pre[l - 1] > 0)
        for l in range(len(Ws)):
            m[l] = b1 * m[l] + (1 - b1) * gW[l]; v[l] = b2 * v[l] + (1 - b2) * gW[l] ** 2
            mb[l] = b1 * mb[l] + (1 - b1) * gb[l]; vb[l] = b2 * vb[l] + (1 - b2) * gb[l] ** 2
            Ws[l] -= lr * (m[l] / (1 - b1 ** t)) / (np.sqrt(v[l] / (1 - b2 ** t)) + eps)
            bs[l] -= lr * (mb[l] / (1 - b1 ** t)) / (np.sqrt(vb[l] / (1 - b2 ** t)) + eps)
    Xt = rng.standard_normal((5000, 2)) * 1.3
    acc = float((((forward(Ws, bs, Xt)[-1][:, 0] > 0) == (teacher_label(Xt) > 0.5))).mean())
    return Ws, bs, acc


def forward(Ws, bs, X):
    """Return [h_0=X, h_1, ..., h_{L-1} (ReLU outputs), logits]."""
    hs = [X]
    h = X
    for l, (W, b) in enumerate(zip(Ws, bs)):
        z = h @ W.T + b
        h = relu(z) if l < len(Ws) - 1 else z
        hs.append(h)
    return hs


# ----------------------------------------------------------------------------- the dictionary

def read_faces(Ws, bs, X, restrict, input_face, rng):
    """Apply the dictionary.  Returns per layer: ids (n,), patterns (bool matrix of distinct faces),
    dims (per sample), and the layer names."""
    hs = forward(Ws, bs, X)
    L = len(Ws)
    layers = []   # list of dicts: name, ids, patterns(bool, distinct), dim(per sample)
    # layer 0
    if input_face == "pattern":
        B0 = X > 0
        ids0, pats0 = ids_of_rows(B0)
    else:
        ids0 = np.zeros(len(X), dtype=np.int64)
        pats0 = np.ones((1, 1), dtype=bool)
    layers.append(dict(name="input", ids=ids0, patterns=pats0, dim=pats0.sum(1)[ids0]))
    # hidden layers
    for l in range(1, L):
        B = hs[l] > 0
        if restrict is not None and restrict < B.shape[1]:
            cols = np.sort(rng.choice(B.shape[1], restrict, replace=False))
            B = B[:, cols]
        ids, pats = ids_of_rows(B)
        layers.append(dict(name=f"hidden{l}", ids=ids, patterns=pats, dim=pats.sum(1)[ids], width=B.shape[1]))
    # output layer: argmax class, or sign of a single logit
    logits = hs[-1]
    if logits.shape[1] == 1:
        cls = (logits[:, 0] > 0).astype(np.int64)
        pats = np.eye(2, dtype=bool)
    else:
        cls = logits.argmax(1).astype(np.int64)
        pats = np.eye(logits.shape[1], dtype=bool)
    layers.append(dict(name="output", ids=cls, patterns=pats, dim=np.ones(len(X), dtype=np.int64)))
    return layers


# ----------------------------------------------------------------------------- experiments

def E1_quiver(layers, out):
    print("\n== E1  the quiver read off the network =======================================")
    n = len(layers[0]["ids"])
    res = []
    for l in range(len(layers) - 1):
        a, b = layers[l], layers[l + 1]
        src, dst = a["ids"], b["ids"]
        arrows = np.unique(combine(src, dst))
        n_faces_b = int(b["ids"].max()) + 1
        dd = b["dim"] - a["dim"]
        jump = float((np.abs(dd) > 1).mean()) if a["name"] != "input" and b["name"] != "output" else float("nan")
        # heights: number of distinct histories (face sequences from layer 1) arriving at each face of layer l+1
        hist = np.unique(np.stack([layers[k]["ids"] for k in range(1, l + 2)], 1), axis=0)
        h_counts = np.bincount(hist[:, -1], minlength=n_faces_b)
        _, cnt = np.unique(dst, return_counts=True)
        coverage = float(cnt[cnt >= 2].sum() / n)
        row = dict(layer=b["name"], faces=n_faces_b, arrows=int(len(arrows)),
                   mean_dim=float(b["dim"].mean()), width=b.get("width", None),
                   frac_dim_jump_gt1=jump,
                   heights_mean=float(h_counts[h_counts > 0].mean()), heights_max=int(h_counts.max()),
                   frac_faces_unique_history=float((h_counts == 1).sum() / (h_counts > 0).sum()),
                   coverage_seen_twice=coverage)
        res.append(row)
        print(f"  {a['name']:>8} -> {b['name']:<8} faces={row['faces']:5d} arrows={row['arrows']:6d} "
              f"mean d={row['mean_dim']:.1f}" + (f"/{row['width']}" if row['width'] else "") +
              f"  |dd|>1: {row['frac_dim_jump_gt1']:.2f}  heights mean/max={row['heights_mean']:.1f}/{row['heights_max']}"
              f"  unique-history faces={row['frac_faces_unique_history']:.2f}  coverage(seen>=2)={coverage:.2f}")
    out["E1"] = res


def E2_memory(layers, out, rng, n_null=3):
    print("\n== E2  memory: I(next ; past | present), per layer and lag  (nats; null = permutation) ==")
    print("      sufficiency of the barycentre <=> all of these vanish (Note 1, Section 6)")
    res = []
    for l in range(1, len(layers) - 1):
        present, nxt = layers[l]["ids"], layers[l + 1]["ids"]
        H_next_given_present = entropy_from_keys(combine(nxt, present)) - entropy_from_keys(present)
        row = dict(present=layers[l]["name"], next=layers[l + 1]["name"],
                   H_next_given_present=H_next_given_present, lags=[])
        line = f"  present={layers[l]['name']:<8} next={layers[l+1]['name']:<8} H(next|present)={H_next_given_present:.3f}  "
        for lag in range(1, l + 1):
            past = layers[l - lag]["ids"]
            if past.max() == 0:
                continue
            val = cmi(nxt, past, present)
            null = float(np.mean([cmi(nxt, permute_within(past, present, rng), present) for _ in range(n_null)]))
            excess = max(val - null, 0.0)
            frac = excess / H_next_given_present if H_next_given_present > 1e-12 else 0.0
            row["lags"].append(dict(lag=lag, cmi=val, null=null, excess=excess, frac_of_H=frac))
            line += f"lag{lag}: {val:.3f}-{null:.3f}={excess:.3f} ({100*frac:.0f}% of H)  "
        # whole past at once
        if l >= 1:
            past_all = combine(*[layers[k]["ids"] for k in range(0, l)])
            if past_all.max() > 0:
                val = cmi(nxt, past_all, present)
                null = float(np.mean([cmi(nxt, permute_within(past_all, present, rng), present) for _ in range(n_null)]))
                row["whole_past"] = dict(cmi=val, null=null, excess=max(val - null, 0.0))
                line += f"| whole past: {val:.3f}-{null:.3f}={max(val-null,0):.3f}"
        print(line)
        res.append(row)
    out["E2"] = res


def coboundary_projection(src, dst, F, w, n_nodes, iters=2000, tol=1e-12):
    """Weighted least squares  min_U sum_g w_g (F_g - (U[dst]-U[src]))^2  by conjugate gradients
    on the weighted Laplacian.  Returns U and the residual F - dU."""
    def L(U):
        d = U[dst] - U[src]
        out = np.zeros(n_nodes)
        np.add.at(out, dst, w * d)
        np.add.at(out, src, -w * d)
        return out
    rhs = np.zeros(n_nodes)
    np.add.at(rhs, dst, w * F)
    np.add.at(rhs, src, -w * F)
    U = np.zeros(n_nodes)
    r = rhs - L(U)
    p = r.copy()
    rr = r @ r
    for _ in range(iters):
        if rr < tol:
            break
        Lp = L(p)
        alpha = rr / max(p @ Lp, 1e-300)
        U += alpha * p
        r -= alpha * Lp
        rr_new = r @ r
        p = r + (rr_new / rr) * p
        rr = rr_new
    return U, F - (U[dst] - U[src])


def E3_cocycle(layers, out):
    print("\n== E3  the class of the cocycle F = -log P(arrow | source) ====================")
    # global face numbering across layers
    offsets, tot = [], 0
    for lay in layers:
        offsets.append(tot)
        tot += int(lay["ids"].max()) + 1
    src_all, dst_all, F_all, w_all, layer_of = [], [], [], [], []
    per_sample_F = np.zeros(len(layers[0]["ids"]))
    per_sample_F_by_layer = []
    for l in range(len(layers) - 1):
        s, d = layers[l]["ids"], layers[l + 1]["ids"]
        key = combine(s, d)
        uk, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
        inv = inv.reshape(-1)
        src_cnt = np.bincount(s)
        us, ud = s[np.searchsorted(key, uk, sorter=np.argsort(key))], None
        # recover the (source,dest) of each unique arrow
        order = np.argsort(key, kind="stable")
        firsts = order[np.searchsorted(key[order], uk)]
        us, ud = s[firsts], d[firsts]
        F = -np.log(cnt / src_cnt[us])
        src_all.append(us + offsets[l]); dst_all.append(ud + offsets[l + 1])
        F_all.append(F); w_all.append(cnt.astype(np.float64)); layer_of.append(np.full(len(uk), l))
        per_sample_F += F[inv]
        per_sample_F_by_layer.append(per_sample_F.copy())
    src = np.concatenate(src_all); dst = np.concatenate(dst_all)
    F = np.concatenate(F_all); w = np.concatenate(w_all); lo = np.concatenate(layer_of)
    U, resid = coboundary_projection(src, dst, F, w, tot)
    Fc = F - (w * F).sum() / w.sum()   # centre: a constant is a coboundary too (up to the layer index)
    rho = float(np.sqrt((w * resid ** 2).sum() / (w * Fc ** 2).sum()))
    print(f"  arrows={len(F)}  faces={tot}  residual fraction  ||F - dU||_w / ||F - mean||_w = {rho:.3f}")
    print("  (0 would mean P(tau|sigma) = g(sigma) h(tau) on the support of the quiver: the face is sufficient)")
    per_layer = []
    for l in range(len(layers) - 1):
        m = lo == l
        Fl, wl = F[m], w[m]
        var_l = (wl * (Fl - (wl * Fl).sum() / wl.sum()) ** 2).sum()
        if var_l < 1e-12 * max(wl.sum(), 1.0):
            per_layer.append(float("nan"))      # constant F on this layer: nothing to explain
            continue
        _, res_l = coboundary_projection(src[m], dst[m], Fl, wl, tot)   # fit on this layer's arrows only
        per_layer.append(float(np.sqrt((wl * res_l ** 2).sum() / var_l)))
    print("  per-layer residual fractions (U fitted per layer): " + "  ".join(f"{layers[l]['name']}->{layers[l+1]['name']}: {r:.3f}" for l, r in enumerate(per_layer)))
    # within-endpoint spread of the accumulated action, and the centrality defect
    rows = []
    for l in range(1, len(layers)):
        Fmu = per_sample_F_by_layer[l - 1]
        end = layers[l]["ids"]
        tot_var = float(Fmu.var())
        # within-endpoint variance
        order = np.argsort(end, kind="stable")
        e_sorted = end[order]; f_sorted = Fmu[order]
        bounds = np.flatnonzero(np.diff(e_sorted)) + 1
        groups = np.split(f_sorted, bounds)
        within = sum(len(g) * g.var() for g in groups) / len(Fmu)
        R2 = 1.0 - within / tot_var if tot_var > 1e-15 else float("nan")
        # centrality defect: sum_tau P(tau) [ log h_tau - H(history | tau) ]
        hist_key = np.unique(np.stack([layers[k]["ids"] for k in range(1, l + 1)], 1), axis=0, return_inverse=True)[1].reshape(-1)
        H_hist_given_end = entropy_from_keys(combine(hist_key, end)) - entropy_from_keys(end)
        h_tau = np.bincount(np.unique(np.stack([hist_key, end], 1), axis=0)[:, 1])
        _, ecnt = np.unique(end, return_counts=True)
        E_log_h = float((ecnt / ecnt.sum() * np.log(h_tau[h_tau > 0])).sum())
        defect = E_log_h - H_hist_given_end
        rows.append(dict(endpoint=layers[l]["name"], sd_total=float(np.sqrt(tot_var)), sd_within=float(np.sqrt(within)),
                         R2_endpoint=R2, centrality_defect=defect, E_log_heights=E_log_h))
        print(f"  endpoint={layers[l]['name']:<8} action F(mu): sd total={np.sqrt(tot_var):.3f} within-endpoint={np.sqrt(within):.3f} "
              f"R^2(endpoint)={R2:.3f} | centrality defect = E log h - H(hist|end) = {E_log_h:.3f}-{H_hist_given_end:.3f} = {defect:.3f} nats")
    out["E3"] = dict(residual_fraction=rho, per_layer_residual=per_layer, endpoints=rows)


def E4_tolerance(layers, out, rng, min_count=20, max_faces=400, n_null=5):
    print("\n== E4  composition fork: overlap cosine of faces vs fidelity of their outgoing kernels ==")
    print("      the groupoid composition is blind to overlaps; the frame composition weights by G.")
    print("      a positive relation means the kernel is G-continuous: the tolerance structure carries information")
    res = []
    for l in range(1, len(layers) - 1):
        ids, pats, nxt = layers[l]["ids"], layers[l]["patterns"], layers[l + 1]["ids"]
        if nxt.max() == 0:
            print(f"  {layers[l]['name']}: next layer has a single face; skipped")
            continue
        cnt = np.bincount(ids)
        nonempty = pats.sum(1) > 0
        big = np.flatnonzero((cnt >= min_count) & nonempty[:len(cnt)])
        if len(big) < 8:
            print(f"  {layers[l]['name']}: only {len(big)} faces with >= {min_count} samples; skipped")
            continue
        if len(big) > max_faces:
            big = rng.choice(big, max_faces, replace=False)
        M = int(nxt.max()) + 1
        K = np.zeros((len(big), M))
        pos = {f: i for i, f in enumerate(big)}
        sel = np.isin(ids, big)
        rows = np.array([pos[f] for f in ids[sel]])
        np.add.at(K, (rows, nxt[sel]), 1.0)
        K /= K.sum(1, keepdims=True)
        sq = np.sqrt(K)
        Fid = sq @ sq.T
        B = pats[big].astype(np.float64)
        d = B.sum(1)
        G = (B @ B.T) / np.sqrt(np.outer(d, d))
        # transported overlap: expected overlap cosine of the targets, E[G(tau,tau') | tau~P(.|s), tau'~P(.|s')]
        Pn = layers[l + 1]["patterns"].astype(np.float64)
        dn = np.maximum(Pn.sum(1), 1.0)          # an empty face has overlap 0 with everything
        Gn = (Pn @ Pn.T) / np.sqrt(np.outer(dn, dn))
        Gout = K @ Gn @ K.T
        iu = np.triu_indices(len(big), 1)
        g, f, go = G[iu], Fid[iu], Gout[iu]
        if g.std() < 1e-12 or f.std() < 1e-12:
            print(f"  {layers[l]['name']}: overlap or fidelity constant over the {len(big)} faces; skipped")
            continue
        r = float(np.corrcoef(g, f)[0, 1])
        r_out = float(np.corrcoef(g, go)[0, 1]) if go.std() > 1e-12 else float("nan")
        nulls = []
        for _ in range(n_null):
            Bp = B[rng.permutation(len(big))]
            dp = Bp.sum(1)
            Gp = (Bp @ Bp.T) / np.sqrt(np.outer(dp, dp))
            nulls.append(float(np.corrcoef(Gp[iu], f)[0, 1]))
        bins = [(0.0, 0.5), (0.5, 0.8), (0.8, 1.0001)]
        binned = [(float(f[(g >= a) & (g < b)].mean()) if ((g >= a) & (g < b)).any() else float("nan"), int(((g >= a) & (g < b)).sum())) for a, b in bins]
        res.append(dict(layer=layers[l]["name"], n_faces=int(len(big)), corr=r, null_corr=float(np.mean(nulls)),
                        null_sd=float(np.std(nulls)), corr_G_transported=r_out, fidelity_by_overlap_bin=binned))
        print(f"  {layers[l]['name']:<8} faces={len(big):4d}  corr(G, Fid)={r:+.3f}  null={np.mean(nulls):+.3f}+-{np.std(nulls):.3f}  "
              f"corr(G_in, E G_out)={r_out:+.3f}  "
              f"mean Fid for G<0.5: {binned[0][0]:.3f} (n={binned[0][1]}), 0.5<=G<0.8: {binned[1][0]:.3f} (n={binned[1][1]}), G>=0.8: {binned[2][0]:.3f} (n={binned[2][1]})")
    out["E4"] = res


def E6_betti(layers, out, max_width=48, max_tri=40000):
    print("\n== E6  Betti numbers b0, b1 of the 2-skeleton of each layer nerve (patterns only) =====")
    res = []
    for l in range(1, len(layers) - 1):
        P = layers[l]["patterns"]
        w = P.shape[1]
        if w > max_width:
            print(f"  {layers[l]['name']}: width {w} > {max_width}; skipped")
            continue
        active = np.flatnonzero(P.any(0))
        P = P[:, active].astype(np.int64)
        co = (P.T @ P) > 0          # pairs co-active in some face
        V = P.shape[1]
        edges = [(i, j) for i in range(V) for j in range(i + 1, V) if co[i, j]]
        eidx = {e: k for k, e in enumerate(edges)}
        tris = []
        for i in range(V):
            for j in range(i + 1, V):
                if not co[i, j]:
                    continue
                both = (P[:, i] & P[:, j]).astype(bool)
                if not both.any():
                    continue
                third = np.flatnonzero(P[both].any(0))
                for k in third:
                    if k > j:
                        tris.append((i, j, k))
        if len(tris) > max_tri:
            print(f"  {layers[l]['name']}: {len(tris)} triangles > {max_tri}; skipped")
            continue
        d1 = np.zeros((len(edges), V))
        for k, (i, j) in enumerate(edges):
            d1[k, i] = -1; d1[k, j] = 1
        r1 = np.linalg.matrix_rank(d1) if len(edges) else 0
        if tris:
            d2 = np.zeros((len(tris), len(edges)))
            for k, (i, j, m) in enumerate(tris):
                d2[k, eidx[(j, m)]] = 1; d2[k, eidx[(i, m)]] = -1; d2[k, eidx[(i, j)]] = 1
            r2 = np.linalg.matrix_rank(d2)
        else:
            r2 = 0
        b0 = V - r1
        b1 = len(edges) - r1 - r2
        res.append(dict(layer=layers[l]["name"], vertices=V, edges=len(edges), triangles=len(tris), b0=int(b0), b1=int(b1)))
        print(f"  {layers[l]['name']:<8} V={V} E={len(edges)} T={len(tris)}  b0={b0}  b1={b1}   (full simplex on V would have E={V*(V-1)//2}, b1=0)")
    out["E6"] = res


# ----------------------------------------------------------------------------- main

def markov_surrogate(layers, rng):
    """Histories sampled from the fitted layer-by-layer kernels: same quiver, same arrows and
    multiplicities in expectation, but no memory beyond the current face.  Used to calibrate E2/E3."""
    n = len(layers[0]["ids"])
    sur = [dict(layers[0])]
    cur = layers[0]["ids"].copy()
    for l in range(1, len(layers)):
        s, d = layers[l - 1]["ids"], layers[l]["ids"]
        M = int(d.max()) + 1
        # sample next face for each sample given its current face, from the empirical kernel
        order = np.argsort(s, kind="stable")
        s_sorted, d_sorted = s[order], d[order]
        bounds = np.r_[0, np.flatnonzero(np.diff(s_sorted)) + 1, len(s)]
        nxt = np.empty(n, dtype=np.int64)
        for i in range(len(bounds) - 1):
            face = s_sorted[bounds[i]]
            pool = d_sorted[bounds[i]:bounds[i + 1]]
            where = np.flatnonzero(cur == face)
            if len(where):
                nxt[where] = rng.choice(pool, len(where), replace=True)
        # faces of the current layer that never occur as sources (cannot happen: every face has a successor)
        sur.append(dict(name=layers[l]["name"], ids=nxt, patterns=layers[l]["patterns"], dim=layers[l]["patterns"].sum(1)[nxt]))
        cur = nxt
    return sur


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["trained", "random", "weights"], default="trained")
    ap.add_argument("--weights", help=".npz with W1,b1,...; implies --mode weights")
    ap.add_argument("--inputs", help=".npy inputs (N, n_0); default N(0,I)")
    ap.add_argument("--widths", default="2,24,24,24,1", help="stand-in architecture, comma separated")
    ap.add_argument("--n-samples", type=int, default=40000)
    ap.add_argument("--restrict", type=int, default=None, help="read faces on a random sub-frame of k units per layer")
    ap.add_argument("--input-face", choices=["trivial", "pattern"], default="trivial")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--train-steps", type=int, default=4000)
    ap.add_argument("--self-check", action="store_true", help="also run E2/E3 on a Markov surrogate of the histories")
    ap.add_argument("--json", help="write the results to this json file")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    if args.weights:
        args.mode = "weights"
    out = dict(args=vars(args))
    if args.mode == "weights":
        Ws, bs = load_weights(args.weights)
        print(f"loaded {args.weights}: widths {[Ws[0].shape[1]] + [W.shape[0] for W in Ws]}")
    else:
        widths = [int(x) for x in args.widths.split(",")]
        if args.mode == "trained":
            Ws, bs, acc = train_net(widths, rng, steps=args.train_steps)
            print(f"trained stand-in {widths} on the synthetic 2-d task: test accuracy {acc:.3f}  ({time.time()-t0:.1f}s)")
            out["accuracy"] = acc
        else:
            Ws, bs = random_net(widths, rng)
            print(f"random stand-in {widths}")
    n0 = Ws[0].shape[1]
    if args.inputs:
        X = np.load(args.inputs).astype(np.float64)
        if args.n_samples < len(X):
            X = X[rng.choice(len(X), args.n_samples, replace=False)]
    else:
        X = rng.standard_normal((args.n_samples, n0)) * (1.3 if n0 == 2 else 1.0)
    print(f"samples: {len(X)}  input dim {n0}  dictionary: faces=ReLU patterns"
          + (f" restricted to {args.restrict} units" if args.restrict else "") + f", input face={args.input_face}, output face=argmax class")
    layers = read_faces(Ws, bs, X, args.restrict, args.input_face, rng)
    E1_quiver(layers, out)
    E2_memory(layers, out, rng)
    E3_cocycle(layers, out)
    E4_tolerance(layers, out, rng)
    E6_betti(layers, out)
    if args.self_check:
        print("\n######## self-check: the same quiver with Markov histories (no memory beyond the face) ########")
        sur = markov_surrogate(layers, rng)
        out_sur = {}
        E2_memory(sur, out_sur, rng)
        E3_cocycle(sur, out_sur)
        out["markov_surrogate"] = out_sur
    print(f"\ndone in {time.time()-t0:.1f}s")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=1, default=float)


if __name__ == "__main__":
    main()
