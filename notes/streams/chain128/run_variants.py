"""Run chain variants on baked ground truth and record per-layer MSE of the post-activation means.

    python run_variants.py --dataset /root/work/truth128 --mlps 0,1 --variants A,B,C --out results/raw
Variants (see REPORT.md):
  A0   paper's Algorithm 2 (K=2, leading-order off-diagonal update)
  A    K=2 covariance propagation, Gaussian closure to all orders in C
  M    dense K=3, slices only (memoryless all-distinct part)
  B    dense K=3, Wick-only births (Gaussian rho^2 + Phi^3 kappa3(z))
  C    dense K=3, first-order closure with leg-partition coefficients
  D    dense K=3, per-layer fitted coefficients (--coefs table.npy)
  C0   C with kappa4(z) = 0
  C211 C with the (2,1,1) part of kappa4 zeroed
  CE-atlas / CE-atlas211 / CE-reg211   C with the fourth cumulant teacher-forced from --atlas (all / (2,1,1) / u C)
  F<l>  C with kappa3(a_l) replaced by the atlas's (teacher forcing at layer l), needs --atlas
  FD<l> same, only the all-distinct part
  N<eps>  C with a relative rms perturbation eps of D21(z) at every layer (error law)
"""
import argparse, glob, json, os, time
import numpy as np
import chain as ch


def load_truth(d):
    import pyarrow.parquet as pq
    files = sorted(glob.glob(os.path.join(d, "data", "*.parquet")))
    out = []
    for f in files:
        t = pq.read_table(f)
        cols = t.column_names
        for r in range(t.num_rows):
            row = {c: t.column(c)[r].as_py() for c in cols if c in ("mlp_name", "mlp_seed", "weights", "all_layer_means", "avg_variance", "final_layer_variance", "all_layer_variances")}
            out.append(row)
    return out


_TRUTH_CACHE = {}


def truth_objs(z):
    key = id(z)
    if key not in _TRUTH_CACHE:
        L = z["weights"].shape[0]
        objs = []
        for l in range(L):
            st = ch.atlas_state(z, l)
            D3, D21, K4, K31, K22 = ch.slices_from_state(st)
            objs.append(dict(var=st.var, D3=D3, D21=D21, K4=K4, K22=K22, Coff=ch.ok.offdiag(st.C)))
        _TRUTH_CACHE.clear(); _TRUTH_CACHE[key] = objs
    return _TRUTH_CACHE[key]


def interface_errors(rec, T):
    """per-layer relative rms errors of the chain's pre-activation objects vs the atlas (atlas noise included)."""
    def rel(a, b):
        return float(np.sqrt(np.sum((a - b) ** 2) / max(np.sum(b ** 2), 1e-300)))
    out = {k: [] for k in ("var", "Coff", "D3", "D21", "K4", "K22")}
    for l, r in enumerate(rec):
        out["var"].append(rel(r["var"], T[l]["var"])); out["Coff"].append(rel(ch.ok.offdiag(r["C"]), T[l]["Coff"]))
        for k in ("D3", "D21", "K4", "K22"):
            out[k].append(rel(r[k], T[l][k]) if l else 0.0)
    return out


def run_variant(v, W, atlas=None, coefs=None):
    if v == "A0":
        return dict(means=ch.paper_k2(W), rec=[])
    if v == "A":
        return ch.Chain(W, k3mode="k2", record=False).run()
    base = dict(record=atlas is not None)
    if "@" in v:                      # perturbation @ base chain:  N<eps>@k3:k4 | F<l>@k3:k4 | FD<l>@k3:k4 | FA<l>@k3:k4
        pert, b = v.split("@")
        k3m, k4m = b.split(":")[:2]
        kw = dict(k3mode=k3m, k4mode=k4m, coefs=coefs, atlas=atlas, **base)
        if pert.startswith("N"):
            eps, seed = (pert[1:].split("s") + ["1234"])[:2]
            kw["d21_noise"] = (float(eps), int(seed))
        elif pert.startswith("FD"):
            kw["force"] = {int(pert[2:]): "k3d"}
        elif pert.startswith("FA"):
            kw["force"] = {int(pert[2:]): "all"}
        elif pert.startswith("F"):
            kw["force"] = {int(pert[1:]): "k3"}
        return ch.Chain(W, **kw).run()
    if ":" in v:                      # generic  k3mode:k4mode[:order]
        parts = v.split(":")
        kw = dict(k3mode=parts[0], k4mode=parts[1], coefs=coefs, atlas=atlas)
        if len(parts) > 2:
            kw["order"] = int(parts[2])
        return ch.Chain(W, **kw, **base).run()
    if v == "M":
        return ch.Chain(W, k3mode="none", **base).run()
    if v == "B":
        return ch.Chain(W, k3mode="wick", **base).run()
    if v == "C":
        return ch.Chain(W, k3mode="closure", **base).run()
    if v == "D":
        return ch.Chain(W, k3mode="fit", coefs=coefs, **base).run()
    if v in ("C0", "B0", "M0", "D0"):
        mode = dict(C0="closure", B0="wick", M0="none", D0="fit")[v]
        return ch.Chain(W, k3mode=mode, coefs=coefs, k4mode="zero", **base).run()
    if v == "C211":
        return ch.Chain(W, k3mode="closure", k4mode="zero211", **base).run()
    if v.startswith("CE-"):
        return ch.Chain(W, k3mode="closure", k4mode=v[3:], atlas=atlas, **base).run()
    if v.startswith("DE-"):
        return ch.Chain(W, k3mode="fit", coefs=coefs, k4mode=v[3:], atlas=atlas, **base).run()
    if v.startswith("FD"):
        return ch.Chain(W, k3mode="closure", force={int(v[2:]): "k3d"}, atlas=atlas, **base).run()
    if v.startswith("FA"):
        return ch.Chain(W, k3mode="closure", force={int(v[2:]): "all"}, atlas=atlas, **base).run()
    if v.startswith("F"):
        return ch.Chain(W, k3mode="closure", force={int(v[1:]): "k3"}, atlas=atlas, **base).run()
    if v.startswith("N"):
        eps = float(v[1:])
        return ch.Chain(W, k3mode="closure", d21_noise=(eps, 1234), **base).run()
    raise ValueError(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True); ap.add_argument("--mlps", default="0")
    ap.add_argument("--variants", default="A,C"); ap.add_argument("--out", required=True)
    ap.add_argument("--atlas", default=None, help="atlas npz pattern with {i} for the MLP index")
    ap.add_argument("--coefs", default=None)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    truth = load_truth(a.dataset)
    coefs = np.load(a.coefs) if a.coefs else None
    for i in [int(x) for x in a.mlps.split(",")]:
        row = truth[i]
        W = np.asarray(row["weights"], dtype=np.float64)
        gt = np.asarray(row["all_layer_means"], dtype=np.float64)
        atlas = np.load(a.atlas.format(i=i)) if a.atlas else None
        for v in a.variants.split(","):
            fn = os.path.join(a.out, f"{v}_mlp{i}.json")
            if os.path.exists(fn):
                continue
            t0 = time.time()
            out = run_variant(v, W, atlas, coefs)
            m = out["means"]
            mse = np.mean((m - gt) ** 2, axis=1)
            res = dict(variant=v, mlp=i, per_layer_mse=mse.tolist(), final_mse=float(mse[-1]),
                       secs=time.time() - t0, means_final=m[-1].tolist())
            if atlas is not None and out.get("rec"):
                res["eps"] = interface_errors(out["rec"], truth_objs(atlas))
            json.dump(res, open(fn, "w"))
            print(f"mlp {i} {v:>10}: final {mse[-1]:.3e}  mean-over-layers {mse.mean():.3e}  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
