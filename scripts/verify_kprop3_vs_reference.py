"""Layer-by-layer check of whest/kprop3.py against the organisers torch reference (needs torch and the cloned
alignment-research-center/mlp_cumulant_propagation on PYTHONPATH):  PYTHONPATH=.../src:. python scripts/verify_kprop3_vs_reference.py 16 4"""
import logging; logging.disable(logging.WARNING)
import sys, numpy as np, torch
from mlp_kprop.mlp import MLP
from mlp_kprop.kprop_harmonic import mlp_kprop, SIMPLE
from whest.kprop3 import kprop3_chain, linear_step, nonlin_step, K3State
torch.set_default_dtype(torch.float64); torch.set_grad_enabled(False)
n, L = int(sys.argv[1]) if len(sys.argv) > 1 else 16, int(sys.argv[2]) if len(sys.argv) > 2 else 4
rng = np.random.default_rng(0); W = rng.standard_normal((L, n, n)) * np.sqrt(2.0 / n)
mlp = MLP(input_dim=n, hidden_dim=n, output_dim=n, num_layers=L + 1)
with torch.no_grad():
    for l in range(L):
        mlp.Ws[l].weight.copy_(torch.tensor(W[l]))
        if mlp.Ws[l].bias is not None: mlp.Ws[l].bias.zero_()
    mlp.Ws[L].weight.copy_(torch.eye(n))
    if mlp.Ws[L].bias is not None: mlp.Ws[L].bias.zero_()
K_in = {1: torch.zeros(n), 2: torch.eye(n)}
use_avg = len(sys.argv) > 3 and sys.argv[3] == "avg"
ref = mlp_kprop(mlp, K_in, k_max=3, kind=SIMPLE, factor=True, use_avg_metric=use_avg, output_all=True)
rec = {}; out = kprop3_chain(W, record=rec)
for l in range(L):
    pre = ref[f"pre{l}"]; act = ref[f"act{l}"]; r = rec[l]
    m_ref = pre[1].core.numpy(); S_ref = pre[2].core.numpy()
    if 3 in pre: K3 = pre[3]; S21_ref = K3.get_dslice((2, 1)).numpy(); S3_ref = K3.get_dslice((3,)).numpy()
    else: S21_ref = np.zeros((n, n)); S3_ref = np.zeros(n)
    c4_ref = float(pre[4].core) if 4 in pre else float("nan"); metric = pre[4].metric if 4 in pre else None
    mu_ref = act[1].core.numpy(); Ch_ref = act[2].core.numpy(); c4h_ref = float(act[4].core) if 4 in act else float("nan")
    K3h = act[3]; S21h_ref = K3h.get_dslice((2, 1)).numpy(); S3h_ref = K3h.get_dslice((3,)).numpy()
    print(f"layer {l}: pre mean {np.max(np.abs(m_ref - r['m'])):.1e} | pre var {np.max(np.abs(np.diag(S_ref) - r['var'])):.1e} | pre K3_21 {np.max(np.abs(S21_ref - r['K3_21'])):.1e} (scale {np.max(np.abs(S21_ref)):.2e}) | pre K3_3 {np.max(np.abs(S3_ref - r['K3_3'])):.1e} | "
          f"act mean {np.max(np.abs(mu_ref - r['mu_h'])):.1e} | act cov {np.max(np.abs(Ch_ref - r['Ch'])):.1e} | act K3_21 {np.max(np.abs(S21h_ref - r['K3h_21'])):.1e} | act K3_3 {np.max(np.abs(S3h_ref - r['K3h_3'])):.1e} | act c4 ref {c4h_ref:.3e} mine {r['c4']:.3e}", flush=True)
    if l == 0 and metric is not None: print("   metric ndim", metric.ndim, "first entries", metric.flatten()[:3].numpy())
print("final mean max diff:", np.max(np.abs(ref[f"act{L-1}"][1].core.numpy() - out[-1])))
