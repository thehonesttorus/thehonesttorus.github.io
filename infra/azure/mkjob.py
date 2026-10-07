#!/usr/bin/env python3
"""Build the body of an Azure managed run command that runs infra/azure/vmjob.sh on the experiment VM.

  python3 infra/azure/mkjob.py <job-name> '<shell command run in /opt/work/src, outputs to $OUT>'

Writes <scratchpad>/azure/job_<name>.json (mode 600; it contains the container SAS as a protected parameter)
and prints the arm.py command that submits it.
"""
import json, os, sys
S = "/tmp/claude-0/-home-user-thehonesttorus-github-io/2cab30f0-c454-5769-8d9a-8195cd80a5f3/scratchpad/azure"
job, cmd = sys.argv[1], sys.argv[2]
vm = os.environ.get("VM", "claude-vm1")
here = os.path.dirname(os.path.abspath(__file__))
body = {"location": "swedencentral", "properties": {
    "source": {"script": open(os.path.join(here, "vmjob.sh")).read()},
    "parameters": [{"name": "JOB", "value": job}, {"name": "CMD", "value": cmd}],
    "protectedParameters": [{"name": "WORKSAS", "value": open(f"{S}/work_sas").read().strip()}],
    "asyncExecution": False, "timeoutInSeconds": 1800, "treatFailureAsDeploymentFailure": False}}
p = f"{S}/job_{job}.json"
open(p, "w").write(json.dumps(body)); os.chmod(p, 0o600)
print(f"python3 infra/azure/arm.py put \"/resourceGroups/rg-claude-chain-20261007/providers/Microsoft.Compute/virtualMachines/{vm}/runCommands/{job}?api-version=2024-07-01\" @{p}")
