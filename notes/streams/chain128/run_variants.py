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


def run_variant(v, W, atlas=None, coefs=None):
    if v == "A0":
        return ch.paper_k2(W)
    if v == "A":
        return ch.Chain(W, k3mode="k2", record=False).run()["means"]
    base = dict(record=False)
    if v == "M":
        return ch.Chain(W, k3mode="none", **base).run()["means"]
    if v == "B":
        return ch.Chain(W, k3mode="wick", **base).run()["means"]
    if v == "C":
        return ch.Chain(W, k3mode="closure", **base).run()["means"]
    if v == "D":
        return ch.Chain(W, k3mode="fit", coefs=coefs, **base).run()["means"]
    if v == "C0":
        return ch.Chain(W, k3mode="closure", k4mode="zero", **base).run()["means"]
    if v == "C211":
        return ch.Chain(W, k3mode="closure", k4mode="zero211", **base).run()["means"]
    if v.startswith("CE-"):
        return ch.Chain(W, k3mode="closure", k4mode=v[3:], atlas=atlas, **base).run()["means"]
    if v.startswith("DE-"):
        return ch.Chain(W, k3mode="fit", coefs=coefs, k4mode=v[3:], atlas=atlas, **base).run()["means"]
    if v.startswith("FD"):
        return ch.Chain(W, k3mode="closure", force={int(v[2:]): "k3d"}, atlas=atlas, **base).run()["means"]
    if v.startswith("FA"):
        return ch.Chain(W, k3mode="closure", force={int(v[2:]): "all"}, atlas=atlas, **base).run()["means"]
    if v.startswith("F"):
        return ch.Chain(W, k3mode="closure", force={int(v[1:]): "k3"}, atlas=atlas, **base).run()["means"]
    if v.startswith("N"):
        eps = float(v[1:])
        return ch.Chain(W, k3mode="closure", d21_noise=(eps, 1234), **base).run()["means"]
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
            m = run_variant(v, W, atlas, coefs)
            mse = np.mean((m - gt) ** 2, axis=1)
            json.dump(dict(variant=v, mlp=i, per_layer_mse=mse.tolist(), final_mse=float(mse[-1]),
                           secs=time.time() - t0, means_final=m[-1].tolist()), open(fn, "w"))
            print(f"mlp {i} {v:>10}: final {mse[-1]:.3e}  mean-over-layers {mse.mean():.3e}  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
