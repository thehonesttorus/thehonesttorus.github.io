#!/usr/bin/env python3
"""Fetch job results from the experiment storage container.

  python3 infra/azure/results.py <job> [--wait SECONDS]   download work/results/<job>/* into <scratchpad>/azure/results/<job>/
                                                          (waiting up to SECONDS for the DONE marker) and print log.txt
"""
import os, sys, time
from azure.storage.blob import ContainerClient
S = "/tmp/claude-0/-home-user-thehonesttorus-github-io/2cab30f0-c454-5769-8d9a-8195cd80a5f3/scratchpad/azure"
job = sys.argv[1]; wait = int(sys.argv[sys.argv.index("--wait") + 1]) if "--wait" in sys.argv else 0
c = ContainerClient.from_container_url(open(f"{S}/work_sas").read().strip())
t0 = time.time()
while True:
    names = [b.name for b in c.list_blobs(name_starts_with=f"results/{job}/")]
    if any(n.endswith("/DONE") for n in names) or time.time() - t0 >= wait: break
    time.sleep(20)
d = f"{S}/results/{job}"; os.makedirs(d, exist_ok=True)
for n in names:
    p = os.path.join(d, n[len(f"results/{job}/"):]); os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, "wb").write(c.download_blob(n).readall())
done = os.path.exists(f"{d}/DONE")
print(f"job {job}: {'DONE exit ' + open(f'{d}/DONE').read().strip() if done else 'not finished'} | {len(names)} files in {d}")
if os.path.exists(f"{d}/log.txt"): print(open(f"{d}/log.txt").read()[-4000:])
