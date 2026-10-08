# Workbench: how to resume the omitted-classes and continuation-equivalence work

On 2026-10-08 the Azure experiment VMs became unavailable, and the working code was moved here from a session scratchpad.
The research record is `notes/omitted-classes/README.md`; sections 3h-3j are the latest.

## Layout

- `k3work/` holds the estimator and every analysis script.
  - **The estimator.** `est_v29.py` is the production chain. Run it through `oracle_one.py`, which sets the production
    environment.
  - **Audit switches in `est_v29.py`.**
    - `V37_ORACLE`: subset of D3 D21 G4 WK4M K31 VAR COFF MU, read from `V37_ORACLE_FILE`.
    - `V37_ORACLE_LAYERS`: comma list of layers; empty means every layer.
    - `V39_KD`: the G D core.
    - `V40_FB_SX` and `V40_FB_SY`: the feedback weights.
  - **`oracle_one.py` environment.**
    - `KEEP_EXTRA`: extra dump keys, e.g. `pk1v,K2v,K11`.
    - `DUMP_OUT`: dump path.
    - `SAVE_OUT`: saves only the predicted means of every layer.
- `num12/` and `ncgprob/` are helper modules that the old job bundle shipped alongside; `official/` scripts use
  `PYTHONPATH=../num12`.
- `official/` holds `extract.py`, the evaluation scripts and `pipeline.sh`. The weights themselves (6.3 GB) are not
  committed.
- `stagemap/` holds read-only maps of the transport stages of `est_v29.py`, for the kappa3 legs and the kappa4 sector.
  They come from a design workflow on 2026-10-08; its other parts had not finished.

## Software

Python 3.12 in a venv, with:
- numpy and scipy (OpenBLAS);
- pyarrow;
- `flopscope==0.12.1` and `whestbench==0.16.1`.

There is no GPU code.

## Data to regenerate (it lived only on the VM disks)

1. **Official networks.**
   - Source: the Hugging Face dataset `aicrowd/arc-whestbench-public-2026`, revision `v2-phase2`, files
     `data/mini-0000{0..6}-of-00007.parquet` (100 networks).
   - Extraction: in `official/`, run `python extract.py <shard>.parquet` for each shard. It writes `W_off{id}.npy`
     (column convention, h_l = relu(W h_(l-1))) and `truth_off{id}.npz`.
   - The scripts expect `../official` relative to `k3work/`.
2. **Monte Carlo truth** (run in `k3work/`). `python mcstats.py NET 1.6e7 SEED PREFIX [post]` writes
   `PREFIX_{full,h0,h1}.npz`. The sets used, and what each is for:
   - `mc2_off0` (seed 301) and `mc2_off1` (seed 302), no `post`: oracle inputs;
   - `mc4_off0` (seed 401, `post`) and `mc4_off1` (seed 402, `post`): the reference truth for fentry, kquery, nullalign
     and cread;
   - `mc5_off1` (seed 403, `post`): a replicate;
   - legs: `python mclegs.py NET 1.6e7 SEED mc4_off<NET>`, with the same seeds (used in section 3c).

   Each 1.6e7-sample pass took about an hour on 96 vCPUs. On small instances, split each pass into chunks by seed and
   merge them. The accumulators are sums, but `mcstats.py` has no merge step yet.
3. **Chain dumps** (in `k3work/`).
   - Free-running:
     `KEEP_EXTRA=pk1v,K2v,K11 DUMP_OUT=chain_off<NET>_fk_<tag>.npz python oracle_one.py NET dump none free`, with these
     tags:
     - `base`, no extra environment;
     - `y2`, with `V40_FB_SY=2`;
     - `kd`, with `V39_KD=1`.
   - All-oracle, one per independent half:
     `KEEP_EXTRA=K11,pk1v,K2v DUMP_OUT=chain_off<NET>_or_h<h>.npz python oracle_one.py NET dump:D3+D21+G4+WK4M+K31+VAR+COFF mc2_off<NET>_h<h>.npz or<h>`,
     for h = 0 and 1.

## Pending, in order

1. **The second-moment readout at true inputs.**
   `python cread.py NET mc4_off<NET> chain_off<NET>_or_h0.npz chain_off<NET>_or_h1.npz chain_off<NET>_fk_base.npz`
2. **The continuation-equivalence tests of section 3j** (gain-amplitude consistency, the metric of the kappa4 trace core,
   null-companion alignment).
   `python nullalign.py NET mc4_off<NET> chain_off<NET>_fk_base.npz chain_off<NET>_fk_y2.npz chain_off<NET>_fk_kd.npz`
3. **Intervention telescoping.** Run with `V37_ORACLE=D3+D21+G4+WK4M+K31+VAR+COFF+MU`, `V37_ORACLE_LAYERS=0,...,k` and
   `SAVE_OUT` for k = -1..15. Consecutive differences give the output effect of the error created at each layer,
   propagated by the chain itself.
4. **The null-source audit** Omega(G) = R4(Sigma.G) + R3(mu.G)/4 + R2(G)/12, whose exact value is 0. It needs injection
   hooks that put mu.G/4 into the kappa3 legs, not only the readouts; see `stagemap/map_linear-legs.json`.

## Compute notes

On the old VMs, a chain run took under 2 minutes at 16 BLAS threads, and an analysis script under 2 minutes at 48
threads. Each job needed roughly 10 GB of memory or less; fentry with three truths needs about 5 GB. The data come to
about 20 GB in total.
