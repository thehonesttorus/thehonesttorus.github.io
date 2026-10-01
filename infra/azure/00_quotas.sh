#!/usr/bin/env bash
# Print the vCPU quota headroom for the candidate regions and VM families, so REGIONS /
# VM_SIZE / MAX_NODES_PER_POOL in env.sh can be set to what the subscription actually allows.
# Quota increases are requested in the portal (Subscriptions > Usage + quotas) or with
#   az quota create --resource-name standardDDSv5Family --scope /subscriptions/$SUB_ID/providers/Microsoft.Compute/locations/<region> --limit-object value=1000 name=standardDDSv5Family
source "$(dirname "$0")/env.sh"
printf "%-18s %-28s %8s %8s %8s\n" region family current limit free
for r in $REGIONS; do
  az vm list-usage --location "$r" -o json 2>/dev/null | python3 -c '
import json,sys
r=sys.argv[1]
rows=json.load(sys.stdin)
for x in rows:
    name=x["name"]["value"]
    if name in ("cores","standardDDSv5Family","standardDSv5Family","standardFSv2Family","standardDv5Family","standardEDSv5Family","lowPriorityCores","standardNCADSA100v4Family","standardNCASv3_T4Family"):
        cur,lim=int(x["currentValue"]),int(x["limit"])
        print(f"{r:<18} {name:<28} {cur:8d} {lim:8d} {lim-cur:8d}")
' "$r"
done
echo
echo "Batch service quotas (per Batch account, if the account exists):"
for r in $REGIONS; do
  acct=$(batch_account_name "$r")
  az batch account show -n "$acct" -g "$RG" -o json 2>/dev/null | python3 -c '
import json,sys
a=json.load(sys.stdin); q=a.get("dedicatedCoreQuota"); lp=a.get("lowPriorityCoreQuota"); fam=a.get("dedicatedCoreQuotaPerVMFamily") or []
print(sys.argv[1], "dedicated:", q, "lowpri:", lp, "families:", {f["name"]: f["coreQuota"] for f in fam if f.get("coreQuota")})
' "$r" || true
done
