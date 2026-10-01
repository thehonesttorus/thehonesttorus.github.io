"""One line per `whest run --format json` report: failures, C/B, residual / wall per MLP, MSE."""
import json, sys
try:
    d = json.load(open(sys.argv[1]))
except Exception as e:  # noqa: BLE001
    print(f"{sys.argv[2]:14s} NO REPORT ({e}); see {sys.argv[1]}.err"); sys.exit(0)
r = d["results"]
pm = r.get("per_mlp", [])
res = [p.get("residual_wall_time_s", 0) for p in pm]
wall = [p.get("wall_time_s", 0) for p in pm]
cb = [p.get("flops_used", 0) / 2 ** 41 for p in pm]
print(f"{sys.argv[2]:14s} fails {r['n_failed_mlps']}/{len(pm)}  C/B max {max(cb) if cb else 0:.4f}  "
      f"residual s {', '.join(f'{x:.3f}' for x in res)}  wall max {max(wall) if wall else 0:.1f} s  "
      f"final MSE {r['final_layer_mse']:.3e}  adjusted {r['adjusted_final_layer_score']:.3e}  "
      f"failures {r.get('failure_breakdown')}")
