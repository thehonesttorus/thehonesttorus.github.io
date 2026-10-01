#!/usr/bin/env bash
# Print the quota headroom that actually limits the grid, so REGIONS / VM_SIZE /
# MAX_NODES_PER_POOL in env.sh can be set to what the subscription allows.  Run after
# 01_bootstrap.sh (the Batch account quotas only exist once the accounts do).
#
# 1. Batch account core quotas (per region, per VM family).  Batch-service-mode pools (the
#    default, used here) draw ONLY on these, not on the subscription's regional vCPU quotas.
#    New accounts often start at 0 for some families.  Increase: portal > Batch account >
#    Quotas > Request quota increase (or a support request with quotaTicketDetails, type
#    "Dedicated", VMFamily e.g. standardDDSv5Family).
# 2. Batch accounts per region for the subscription (az batch location quotas show).
# 3. Subscription regional vCPU quotas (az vm list-usage): informational only; they apply to
#    user-subscription-mode Batch accounts.  Increase with the quota extension:
#      az extension add --name quota
#      az quota create --resource-name standardDDSv5Family --resource-type dedicated \
#        --scope /subscriptions/$SUB_ID/providers/Microsoft.Compute/locations/<region> --limit-object value=1000
source "$(dirname "$0")/env.sh"
GPU_VM_SIZE="${GPU_VM_SIZE:-Standard_NC24ads_A100_v4}"

echo "== Batch account quotas (these limit the pools); VM_SIZE=$VM_SIZE GPU_VM_SIZE=$GPU_VM_SIZE"
for r in $REGIONS; do
  acct=$(batch_account_name "$r")
  fam=$(az batch location list-skus -l "$r" --query "[?name=='$VM_SIZE'].familyName | [0]" -o tsv 2>/dev/null || true)
  gfam=$(az batch location list-skus -l "$r" --query "[?name=='$GPU_VM_SIZE'].familyName | [0]" -o tsv 2>/dev/null || true)
  [ -n "$fam" ] || echo "$r: $VM_SIZE is not offered by Batch in this region (az batch location list-skus -l $r)"
  if ! acct_json=$(az batch account show -n "$acct" -g "$RG" -o json 2>/dev/null); then
    echo "$r: no Batch account $acct (run 01_bootstrap.sh)"; continue
  fi
  printf '%s' "$acct_json" | python3 -c '
import json, sys
r, fam, gfam = sys.argv[1], sys.argv[2], sys.argv[3]
a = json.load(sys.stdin)
# az prints ...PerVmFamily (camel-cased SDK attribute), the REST API ...PerVMFamily
rows = a.get("dedicatedCoreQuotaPerVmFamily") or a.get("dedicatedCoreQuotaPerVMFamily") or []
fams = {f.get("name"): f.get("coreQuota") for f in rows}
keys = ["dedicatedCoreQuota", "lowPriorityCoreQuota", "poolQuota", "activeJobAndJobScheduleQuota", "dedicatedCoreQuotaPerVmFamilyEnforced"]
cols = [k + "=" + str(a.get(k)) for k in keys] + [(fam or "?") + "=" + str(fams.get(fam)), (gfam or "?") + "=" + str(fams.get(gfam))]
print(r.ljust(16), " ".join(cols))
' "$r" "$fam" "$gfam"
done

echo
echo "== Batch accounts per region (subscription)"
for r in $REGIONS; do
  printf '%-16s %s\n' "$r" "$(az batch location quotas show -l "$r" --query accountQuota -o tsv 2>/dev/null || echo '?')"
done

echo
echo "== Subscription vCPU quotas (informational: only user-subscription-mode Batch accounts use these)"
printf "%-16s %-28s %8s %8s %8s\n" region family current limit free
for r in $REGIONS; do
  az vm list-usage --location "$r" -o json 2>/dev/null | python3 -c '
import json,sys
r=sys.argv[1]
for x in json.load(sys.stdin):
    name=x["name"]["value"]
    if name in ("cores","standardDDSv5Family","standardDSv5Family","standardFSv2Family","standardDv5Family","standardEDSv5Family","lowPriorityCores","standardNCADSA100v4Family","standardNCadsH100v5Family","standardNCASv3_T4Family"):
        cur,lim=int(x["currentValue"]),int(x["limit"])
        print(f"{r:<16} {name:<28} {cur:8d} {lim:8d} {lim-cur:8d}")
' "$r" || true
done
