# Compute setup notes (what this session can and cannot reach)

## AWS (live)
- Credentials for IAM user `claude-runner` are in the environment; region eu-north-1; bucket `claude-whest-9292-97a992`
  holds the 100 official networks under `data/official/` (6.7 GB) and all job results.
- `infra/fleet.py` launches instances tagged `Project=claude-whest, Fleet=whest9` (name `w9-N`), runs a batch file
  of `tag<TAB>command` lines in parallel slots, streams logs to S3. Spot c7i.48xlarge (192 vCPU) is ~$1.5/h.
- Instances stop themselves after 45 idle minutes.
- The IAM user cannot read or request service quotas (`servicequotas:*` denied). Quota increases must be requested
  in the AWS console (EC2 > Limits: "Running On-Demand Standard instances" and "All Standard Spot Instance Requests"
  in eu-north-1 and one US region), or the IAM policy extended with `servicequotas:GetServiceQuota` and
  `servicequotas:RequestServiceQuotaIncrease`.
- Monte Carlo at 1e9 samples per network is far cheaper on GPUs (g6.xlarge L4 spot ~$0.18/h); the G-instance
  quota is unknown for the same reason.

## Modal (live, through the AWS relay)
- Token for workspace `research-4885` is in `~/.modal.toml` here and on the relay instance `w9-1` (mode 600). The
  sandbox proxy cannot carry gRPC, so the Modal client runs on `w9-1`: `infra/modal_relay.py` bundles `whest/`,
  `scripts/`, `infra/modal_app.py` (+ a jobs file) to S3, runs the client over SSM, and the results come back to
  `s3://BUCKET/results9/JOB/` (`infra/fleet.py watch JOB`, `infra/fleet.py get JOB`).
- App `whest9` (deployed): `runcmd` (CPU tasks, code snapshot on volume `whest-data` under `/data/code9/<hash>`,
  official networks at `/data/official`, results at `/data/out9/JOB`), `mc_gpu` (L4 Monte Carlo moments),
  `put_code`, `listdir`. Batch: `python infra/modal_relay.py batch JOB jobs/JOB.tsv --cpu 8 --mem 16`.
- The relay stops itself after 45 idle minutes: `python infra/fleet.py start w9-1`, wait for `READY`.
- Modal's limits: default 100 concurrent containers per function (raise in the dashboard if sweeps need more).

## Google Cloud (not reachable; project not created)
- `gcloud` is installed but has no credentials. A service-account key is the simplest hand-off: set
  `GOOGLE_APPLICATION_CREDENTIALS_JSON` (the key file's contents) and `GCP_PROJECT` as environment secrets.
- A browser-agent prompt to create the project is in `notes/compute/GCP_PROMPT.md`.

## Elicit (live)
- API key works against `https://elicit.com/api/v2` (reports, paper search, research-agent sessions). The monthly
  usage was at 61% before this session; three reports were run (see `notes/stage9/elicit/`).
