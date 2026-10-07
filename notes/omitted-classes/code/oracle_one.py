# One chain run for the oracle attribution (note XXXI).   python oracle_one.py NET ORACLE MCFILE TAG
# ORACLE: none | dump | "+"-joined subset of D3 D21 G4 WK4M K31 VAR COFF. "dump" saves the chain's own per-layer D3, D21,
# var, kappa_4 diagonal (g4row), (2,2) and (3,1) slices (wk4m, wk431) and off-diagonal covariance for comparison.
import sys, os, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
net, orc, mcf, tag = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
os.environ.update({"V29_WARM_JOIN": "1", "V29_WARM_FB": "1", "V17_R_RES": "4", "V26_STRASSEN": "6", "V26_STRASSEN_MIN": "16",
                   "V32_JOIN_SMM": "1", "V32_ROT_SMM": "1", "V32_JOIN_POST": "3", "V32_JP_C": "0.1", "V21_R_OLD": "320",
                   "V24_R_OLD2": "192", "V18_R_FB": "2", "V24_AGE_OLD2": "8", "V33_K4Q": "3", "V33_K4Q_RANK": "4",
                   "V34_OPT": "abcde"})
# "dump:ORACLES" dumps a run with oracles on (the one-step closure test reads it), written to chain_off{net}_o.npz
dump = orc.startswith("dump"); orcs = orc.split(":", 1)[1] if orc.startswith("dump:") else ("" if dump else orc)
if dump:
    os.environ["V29_DUMP_LAYERS"] = ",".join(str(x) for x in range(16))
if orcs and orcs != "none":
    os.environ["V37_ORACLE"] = orcs.replace("+", ","); os.environ["V37_ORACLE_FILE"] = mcf
spec = importlib.util.spec_from_file_location("estv29", "est_v29.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{net}.npy"); mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**44, wall_time_limit_s=900.0, quiet=True):
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
pl = np.mean((out - mt) ** 2, axis=1) if mt.shape == out.shape else None
print(f"net {net} {tag}: raw {np.mean((out[-1] - mt[-1]) ** 2):.5e}" +
      ("" if pl is None else "  per-layer " + " ".join(f"{x:.2e}" for x in pl)), flush=True)
if dump:
    keep = ("D3", "D21", "var", "g4row", "wk4m", "wk431", "C_off") + (("K22", "K4v", "K2v", "mu") if orcs else ())
    keep = keep + tuple(x for x in os.environ.get("KEEP_EXTRA", "").split(",") if x)   # e.g. K21,K3v,K11,pk1v
    got = {}
    for d in mod.DUMPS:
        for k in keep:
            if d.get(k) is not None and (k, d["layer"]) not in got:
                got[(k, d["layer"])] = np.asarray(d[k], np.float32)
    got.update({(f"{k}own", l): np.asarray(v, np.float32) for k, lst in mod.OWN.items() for l, v in lst})
    np.savez(os.environ.get("DUMP_OUT", f"chain_off{net}{'_o' if orcs else ''}.npz"), **{f"{k}_{l}": v for (k, l), v in got.items()})
