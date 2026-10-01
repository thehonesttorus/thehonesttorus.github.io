"""Are the renormalised closure coefficients an ensemble property of (n, L)?  (stream coef-ensemble)

Two stages, so the expensive tensor work is done once per atlas:

  python coef_table.py cache ATLAS.npz OUT.npz
      per layer l = 0..L-2 of a k3+k4 atlas (moment_atlas_np.py --k3 --k4), everything the closure ladder needs in
      transported (n x n) form plus the tensor-space normal equations of the residual fit:
        D21    target D21(l+1) of this atlas
        Tm     transport of the (3,) and (2,1) slices of kappa3(a_l)          (carried exactly by the chain)
        Twick  transport of oracle_k3.wick_model (leading Wick, atlas gates, + Phi^3 kappa3(z))
        Th4    transport of the Gaussian rho^2 Wick term (hermite_model degree 4), the base of the closure
        TB     transports of the 8 basis tensors: oracle_k3.residual_basis B0..B6, and B7 = B6 with the (2,1,1)
               slice replaced by its r = 1 regeneration u_i C_jk (the published chain's core)
        G, b, rr   <B_i, B_j>, <B_i, R>, <R, R> with R = all-distinct kappa3(a) - Gaussian rho^2 Wick
                   (tensor-space least squares of the residual on the basis; B7 included so both variants refit)

  python coef_table.py analyse CACHEDIR [--out DIR]
      groups caches by width; MLPs are identified by mlp seed, a second cache of the same MLP with another sample
      seed is its pair (independent Monte Carlo atlas).  Per width and layer:
        coefficient tables (ensemble fit over all MLPs, tensor space and D21 space), across-MLP spread, and the
        pair-to-pair spread of the per-MLP fit (its Monte Carlo noise);
        leave-one-out eps of D21(l+1) for every MLP (model from the MLP's own atlas A, coefficients from the other
        MLPs), for the ladder
          wick   (a) leading Wick only
          leg    (b) closure with the leg-partition coefficients (oracle_k3.CLOSURE_COEF)
          legR   (b') same with the (2,1,1) slice regenerated as u_i C_jk (the published chain's floor)
          own    (c) per-MLP fitted closure (tensor space, fitted on the held-out MLP itself: an oracle)
          ownD   (c') per-MLP fit in D21 space
          ens    (d) ensemble table, tensor space
          ensD   (d') ensemble table, D21 space (regress transported residual on transported basis)
          ensR   (e) ensemble table refitted with the regenerated (2,1,1) slice, tensor space
          ensRD  (e') the same in D21 space
        evaluated within-atlas (raw) and, for MLPs with a pair atlas B, against B's D21 with B's noise energy
        subtracted (eps_corr = sqrt(eps^2 - eps_noise^2), eps_noise = rel(D21_A, D21_B)/sqrt 2, as in
        oracle_k3.analyse_pair).
"""
import argparse, glob, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "experiments"))
import oracle_k3 as ok  # noqa: E402

NAMES = ["B0 k3z", "B1 D21+C", "B2 w3 D21", "B3 D3", "B4 rho3", "B5 K22", "B6 K211"]
LEG = np.array(ok.CLOSURE_COEF)


def transport(K, W):
    """D21_ab = sum_ijk W_ia W_ja W_kb K_ijk (BLAS for the two big contractions)."""
    T = np.tensordot(K, W, axes=([2], [0]))          # i j b
    T = np.tensordot(T, W, axes=([1], [0]))          # i b a
    return np.einsum("iba,ia->ab", T, W)


def cache(path, out):
    z = np.load(path)
    W = z["weights"].astype(np.float64)
    L, n, _ = W.shape
    keys = ["D21", "Tm", "Twick", "Th4", "TB", "G", "b", "rr", "Tkd"]
    acc = {k: [] for k in keys}
    for l in range(L - 1):
        o = ok.layer_objects(z, l)
        Wn = W[l + 1]
        K3m = ok.slices_only(o["K3a"]); Kd = ok.all_distinct(o["K3a"])
        Kw = ok.wick_model(o["C"], o["Phi"], o["w2"], o["K3z"])
        Kh4 = ok.hermite_model(o["C"], o["mu"], o["var"], 4)
        B = ok.residual_basis(o)
        assert len(B) == 7
        K211 = o["K211"]; Co = ok.offdiag(o["C"])
        u = np.einsum("ijk,jk->i", K211, Co) / float(np.sum(Co * Co))
        K211r = ok.all_distinct(np.einsum("i,jk->ijk", u, Co))
        B.append(ok.all_distinct(ok.sym3(np.einsum("i,j,k,ijk->ijk", o["w2"], o["Phi"], o["Phi"], K211r))))
        R = Kd - Kh4
        X = np.stack([x.ravel() for x in B], 1)
        acc["G"].append(X.T @ X); acc["b"].append(X.T @ R.ravel()); acc["rr"].append(float(R.ravel() @ R.ravel()))
        acc["D21"].append(o["D21"]); acc["Tm"].append(transport(K3m, Wn)); acc["Twick"].append(transport(Kw, Wn))
        acc["Th4"].append(transport(Kh4, Wn)); acc["TB"].append(np.stack([transport(x, Wn) for x in B]))
        acc["Tkd"].append(transport(Kd, Wn))
        print(f"{path} layer {l} done", flush=True)
    res = {k: np.array(v) for k, v in acc.items()}
    np.savez(out, width=n, depth=L, mlp_seed=int(z["mlp_seed"]), sample_seed=int(z["sample_seed"]),
             n_samples=int(z["n_samples"]), **res)


# ---------------------------------------------------------------- analysis
def rel(a, b):
    return float(np.sqrt(np.sum((a - b) ** 2) / np.sum(b ** 2)))


def corr(e, en):
    return float(np.sqrt(max(e * e - en * en, 0.0)))


def tensor_fit(Gs, bs, idx):
    G = sum(g[np.ix_(idx, idx)] for g in Gs); b = sum(v[idx] for v in bs)
    return np.linalg.lstsq(G, b, rcond=1e-12)[0]


def d21_fit(caches, l, idx):
    """least squares of the transported residual (Tkd - Th4) on the transported basis TB[idx], stacked over caches."""
    X = np.concatenate([np.stack([c["TB"][l, k].ravel() for k in idx], 1) for c in caches])
    y = np.concatenate([(c["Tkd"][l] - c["Th4"][l]).ravel() for c in caches])
    return np.linalg.lstsq(X, y, rcond=None)[0]


IDX = [0, 1, 2, 3, 4, 5, 6]          # true (2,1,1) slice
IDXR = [0, 1, 2, 3, 4, 5, 7]         # regenerated slice
MODELS = ["wick", "leg", "legR", "own", "ownD", "ens", "ensD", "ensR", "ensRD"]


def predict(c, l, coef, idx):
    return c["Tm"][l] + c["Th4"][l] + np.tensordot(coef, c["TB"][l, idx], axes=1)


def analyse(cdir, outdir):
    files = sorted(glob.glob(os.path.join(cdir, "*.npz")))
    caches = [dict(np.load(f)) for f in files]
    widths = sorted({int(c["width"]) for c in caches})
    os.makedirs(outdir, exist_ok=True)
    summary = {}
    for w in widths:
        cs = [c for c in caches if int(c["width"]) == w]
        seeds = sorted({int(c["mlp_seed"]) for c in cs})
        A, Bp = {}, {}
        for s in seeds:
            mine = sorted([c for c in cs if int(c["mlp_seed"]) == s], key=lambda c: int(c["sample_seed"]))
            A[s] = mine[0]
            if len(mine) > 1:
                Bp[s] = mine[1]
        L1 = A[seeds[0]]["D21"].shape[0]
        lines = [f"# width {w}: {len(seeds)} MLPs (seeds {seeds[0]}..{seeds[-1]}), {len(Bp)} with a pair atlas; "
                 f"N = {int(A[seeds[0]]['n_samples'])} per atlas"]
        # ---- coefficient tables
        tabs = {}
        for name, idx, space in [("tensor", IDX, "T"), ("D21", IDX, "D"), ("tensor_reg", IDXR, "T"), ("D21_reg", IDXR, "D")]:
            ens, sd, pairsd = [], [], []
            for l in range(L1):
                if space == "T":
                    ens.append(tensor_fit([A[s]["G"][l] for s in seeds], [A[s]["b"][l] for s in seeds], idx))
                    per = np.array([tensor_fit([A[s]["G"][l]], [A[s]["b"][l]], idx) for s in seeds])
                    pp = np.array([tensor_fit([Bp[s]["G"][l]], [Bp[s]["b"][l]], idx) - tensor_fit([A[s]["G"][l]], [A[s]["b"][l]], idx) for s in Bp])
                else:
                    ens.append(d21_fit([A[s] for s in seeds], l, idx))
                    per = np.array([d21_fit([A[s]], l, idx) for s in seeds])
                    pp = np.array([d21_fit([Bp[s]], l, idx) - d21_fit([A[s]], l, idx) for s in Bp])
                sd.append(per.std(0, ddof=1) if len(seeds) > 1 else np.zeros(len(idx)))
                pairsd.append(np.sqrt(np.mean(pp ** 2, 0) / 2) if len(pp) else np.full(len(idx), np.nan))
            tabs[name] = dict(ens=np.array(ens), sd=np.array(sd), pairsd=np.array(pairsd))
            nm = NAMES[:6] + (["B6 K211"] if idx == IDX else ["B6 uC"])
            lines.append(f"\n## coefficients, {name} fit, ensemble over all {len(seeds)} MLPs: value (across-MLP sd of per-MLP fits / MC sd of one per-MLP fit from pairs)")
            lines.append(" l | " + " | ".join(f"{x:>21}" for x in nm))
            for l in range(L1):
                lines.append(f"{l:>2} | " + " | ".join(f"{tabs[name]['ens'][l, k]:+6.2f} ({tabs[name]['sd'][l, k]:4.2f}/{tabs[name]['pairsd'][l, k]:4.2f})" for k in range(len(idx))))
        # ---- leave-one-out
        eps_raw = {m: np.zeros((len(seeds), L1)) for m in MODELS}
        eps_x = {m: np.full((len(seeds), L1), np.nan) for m in MODELS}
        eps_rep = {m: np.full((len(seeds), L1), np.nan) for m in MODELS}
        noise = np.full((len(seeds), L1), np.nan)
        for i, s in enumerate(seeds):
            tr = [t for t in seeds if t != s]
            c = A[s]
            for l in range(L1):
                ens = {}
                if tr:
                    ens["ens"] = (tensor_fit([A[t]["G"][l] for t in tr], [A[t]["b"][l] for t in tr], IDX), IDX)
                    ens["ensD"] = (d21_fit([A[t] for t in tr], l, IDX), IDX)
                    ens["ensR"] = (tensor_fit([A[t]["G"][l] for t in tr], [A[t]["b"][l] for t in tr], IDXR), IDXR)
                    ens["ensRD"] = (d21_fit([A[t] for t in tr], l, IDXR), IDXR)

                def preds(c):
                    p = {"wick": c["Tm"][l] + c["Twick"][l],
                         "leg": predict(c, l, LEG, IDX), "legR": predict(c, l, LEG, IDXR),
                         "own": predict(c, l, tensor_fit([c["G"][l]], [c["b"][l]], IDX), IDX),
                         "ownD": predict(c, l, d21_fit([c], l, IDX), IDX)}
                    for m, (cf, ix) in ens.items():
                        p[m] = predict(c, l, cf, ix)
                    return p
                pred = preds(c)
                for m, p in pred.items():
                    eps_raw[m][i, l] = rel(p, c["D21"][l])
                if s in Bp:
                    tb = Bp[s]["D21"][l]
                    en = rel(c["D21"][l], tb) / np.sqrt(2)
                    noise[i, l] = en
                    predB = preds(Bp[s])
                    for m, p in pred.items():
                        eps_x[m][i, l] = corr(rel(p, tb), en)
                        # MC noise of the model's own inputs (C, kappa3(z), kappa4 slice, slices of kappa3(a)), from the
                        # spread of the same model built on the two atlases; subtracted too -> pure representation error
                        mn = rel(p, predB[m]) * np.linalg.norm(predB[m]) / np.linalg.norm(tb) / np.sqrt(2)
                        eps_rep[m][i, l] = corr(eps_x[m][i, l], mn)
        hdr = " l | noise | " + " ".join(f"{m:>6}" for m in MODELS)
        lines.append(f"\n## held-out eps of D21(l+1), cross-evaluated on the pair atlas and noise-corrected; mean over the {len(Bp)} pair MLPs [max over MLPs of own ens ensD ensR ensRD]")
        lines.append(hdr)
        for l in range(L1):
            lines.append(f"{l:>2} | {np.nanmean(noise[:, l]):5.3f} | " + " ".join(f"{np.nanmean(eps_x[m][:, l]):6.3f}" for m in MODELS)
                         + "   [" + " ".join(f"{np.nanmax(eps_x[m][:, l]):.3f}" for m in ("own", "ens", "ensD", "ensR", "ensRD")) + "]")
        lines.append(f"\n## held-out representation error: as above with the model-input MC noise also subtracted (eps_rep); mean over pair MLPs [max]")
        lines.append(hdr)
        for l in range(L1):
            lines.append(f"{l:>2} | {'':5} | " + " ".join(f"{np.nanmean(eps_rep[m][:, l]):6.3f}" for m in MODELS)
                         + "   [" + " ".join(f"{np.nanmax(eps_rep[m][:, l]):.3f}" for m in ("own", "ens", "ensD", "ensR", "ensRD")) + "]")
        lines.append(f"\n## held-out eps of D21(l+1), within-atlas (raw, includes the atlas's own transported noise), mean over all {len(seeds)} MLPs [max over MLPs of own ens ensD ensR ensRD]")
        lines.append(hdr)
        for l in range(L1):
            lines.append(f"{l:>2} | {'':5} | " + " ".join(f"{eps_raw[m][:, l].mean():6.3f}" for m in MODELS)
                         + "   [" + " ".join(f"{eps_raw[m][:, l].max():.3f}" for m in ("own", "ens", "ensD", "ensR", "ensRD")) + "]")
        txt = "\n".join(lines)
        print(txt, flush=True)
        with open(os.path.join(outdir, f"width{w}.txt"), "w") as f:
            f.write(txt + "\n")
        summary[w] = dict(n_mlps=len(seeds), n_pairs=len(Bp), noise=np.nanmean(noise, 0).tolist(),
                          eps_x={m: np.nanmean(eps_x[m], 0).tolist() for m in MODELS},
                          eps_raw={m: eps_raw[m].mean(0).tolist() for m in MODELS},
                          eps_rep={m: np.nanmean(eps_rep[m], 0).tolist() for m in MODELS},
                          eps_rep_all={m: eps_rep[m].tolist() for m in MODELS},
                          eps_x_all={m: eps_x[m].tolist() for m in MODELS}, eps_raw_all={m: eps_raw[m].tolist() for m in MODELS},
                          seeds=seeds, pair_seeds=sorted(Bp),
                          coef={k: v["ens"].tolist() for k, v in tabs.items()},
                          coef_sd={k: v["sd"].tolist() for k, v in tabs.items()},
                          coef_pairsd={k: v["pairsd"].tolist() for k, v in tabs.items()})
    with open(os.path.join(outdir, "summary.json"), "w") as f:
        json.dump(summary, f)


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    a = sp.add_parser("cache"); a.add_argument("atlas"); a.add_argument("out")
    b = sp.add_parser("analyse"); b.add_argument("cachedir"); b.add_argument("--out", default=os.path.join(HERE, "results"))
    args = ap.parse_args()
    if args.cmd == "cache":
        cache(args.atlas, args.out)
    else:
        analyse(args.cachedir, args.out)


if __name__ == "__main__":
    main()
