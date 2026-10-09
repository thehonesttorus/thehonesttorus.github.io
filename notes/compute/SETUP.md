# Compute setup notes (what this session can and cannot reach)

## AWS (live)
- Credentials for IAM user `claude-runner` are in the environment; region eu-north-1; bucket `claude-whest-9292-97a992`
  holds the 100 official networks under `data/official/` (6.7 GB) and all job results.
- `infra/fleet.py` launches instances tagged `Project=claude-whest, Fleet=whest9` (name `w9-N`), runs a batch file
  of `tag<TAB>command` lines in parallel slots, streams logs to S3. Spot c7i.48xlarge (192 vCPU) is ~$1.5/h.
- Instances stop themselves after 45 idle minutes.
- The IAM user cannot read or request service quotas (`servicequotas:*` denied). Quota increases must be requested
  in the AWS console (EC2 > Limits, "Running On-Demand Standard instances" and "All Standard Spot Instance Requests"
  in eu-north-1 and one US region), or the policy must be extended with `servicequotas:GetServiceQuota` and
  `servicequotas:RequestServiceQuotaIncrease`.
- Monte Carlo at 1e9 samples per network is far cheaper on GPUs (g6.xlarge L4 spot ~$0.18/h); the G-instance
  quota is unknown for the same reason.

## Modal (not reachable from this session)
- No Modal token is present in this environment (`~/.modal.toml`, `MODAL_TOKEN_ID`/`MODAL_TOKEN_SECRET` absent).
- To enable: add `MODAL_TOKEN_ID` and `MODAL_TOKEN_SECRET` in the cloud environment's settings (Network secrets /
  environment variables). A new session picks them up; the runner in `infra/modal/` of the other branch then works.

## Google Cloud (not reachable; project not created)
- `gcloud` is installed but has no credentials. A service-account key is the simplest hand-off: set
  `GOOGLE_APPLICATION_CREDENTIALS_JSON` (the key file's contents) and `GCP_PROJECT` as environment secrets.
- A browser-agent prompt to create the project is in `notes/compute/GCP_PROMPT.md`.

## Elicit (live)
- API key works against `https://elicit.com/api/v2` (reports, paper search, research-agent sessions). The monthly
  usage was at 61% before this session; three reports were run (see `notes/stage9/elicit/`).
