#!/usr/bin/env python3
"""Submit a prepared job (scratchpad/azure/job_<job>.json from mkjob.py) to one of the experiment VMs.

    python3 infra/azure/submit.py <job> [vm]      vm: claude-vm1 (default), claude-vm2, claude-vm3

The job file is VM-agnostic; only the run-command path names the VM. Results land in the same blob
folder (results/<job>/) whichever VM runs it, so results.py is unchanged.
"""
import json, os, subprocess, sys

RG = "rg-claude-chain-20261007"
SCRATCH = os.environ.get("CLAUDE_AZURE_DIR", "/tmp/claude-0/-home-user-thehonesttorus-github-io/"
                         "2cab30f0-c454-5769-8d9a-8195cd80a5f3/scratchpad/azure")


def main():
    job = sys.argv[1]
    vm = sys.argv[2] if len(sys.argv) > 2 else "claude-vm1"
    path = (f"/resourceGroups/{RG}/providers/Microsoft.Compute/virtualMachines/{vm}/runCommands/{job}"
            "?api-version=2024-07-01")
    here = os.path.dirname(os.path.abspath(__file__))
    out = subprocess.run([sys.executable, os.path.join(here, "arm.py"), "put", path, "@" + os.path.join(SCRATCH, f"job_{job}.json")],
                         capture_output=True, text=True)
    try:
        print(job, vm, json.loads(out.stdout)["http"])
    except Exception:
        print(job, vm, "submit failed:", (out.stderr or out.stdout)[-400:])
        sys.exit(1)


if __name__ == "__main__":
    main()
