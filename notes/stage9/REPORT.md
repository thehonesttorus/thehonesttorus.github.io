# Stage 9 session report: localization theory applied to the WhestBench mean-estimation problem

Branch `claude/determined-fermat-9hk45i`. Note: `notes/stage9/ncg_mlp_stage9.pdf` (16 pp). Code: `whest/`, `scripts/`, `infra/`.

## 1. Theory (the stage-8 principle made exact for estimators)

1. **Defect–path theorem** (Thm 1.2 in the note). Any estimator that is exact on point masses has error equal to minus the
   expected total drift of the estimator along any complete localization of the input law; the drift is the second
   variation of the estimator against the scheme's driving covariance. The other chat's "heat-flow defect" is the
   Eldan-path corollary; the theorem is scheme-independent and says which errors a scheme sees.
2. **Capacity of input tilts** (Prop. 2.1). The Euler–Stein defect delta(0) = e(0) - Lap e(0) is a 2% residual of an
   O(1) identity; the error is a tenth of that. Random probes need ~1e9 samples; k-direction tilts have signal-to-noise
   sqrt(k)*1e-2; only the full coordinate sum (2n+1 chain runs) resolves it. Input-side localization is a diagnostic,
   never an estimator component.
3. **Index-shift lemma and the Gram theorem** (Lemma 3.1, Thm 3.2). d_k' = d_{k+1}: a mean tilt raises the Wiener-chaos
   index; births, arms, slices and gate covariances are the first and second derivatives of the Gaussian moment map.
   Exactly: delta(y) = -<Hess G, Gamma(y)> with Gamma the Gram matrix of the first-layer moment gradients (closed
   forms). Verified to 7e-7 on a tiny network. The organisers' Wick coefficients are these index shifts, and their
   36-term table is the first-order tilt response.
4. **Scale-mixture lemma** (Lemma 4.1). A fluctuating collective scale produces a (2,2)-type kappa4 = 4v Sym(S x S);
   the Gaussian closure errs by -(v/2) sigma (alpha^2-1) phi and the kappa4 mean readout cancels it exactly at O(v).
   The organisers' chain keeps precisely this (one radial scalar per layer) and nothing else of kappa4.
5. **Sensitivity law** (Prop. 4.2): a relative variance error eps moves a mean by eps*sigma*phi/2; 1e-3 on the
   diagonal costs 7x in MSE; the bulk of the covariance may be off by percents; the collective mode is the central
   variable of the stage-8 split.
6. **Cocycle and the Hadamard obstruction** (Sec. 6). Old kappa3 sources are the images of their births under the
   gated transport cocycle, whose effective rank decays 164, 92, 65, 47, 29, 23 with age; but the slices are read
   through Hadamard products of the legs, so Frobenius truncation of the legs costs 2*eps in the slices; the memory
   is compressible only once it has become unimportant.

## 2. Measurements (official network 0, n=1024, L=16, truth from 1e9 samples)

| estimator | final MSE | cost |
|---|---|---|
| Gaussian closure | 4.06e-6 | 0.024 B (f32, symmetric einsum) |
| own kappa3 chain (k3v3) | 5.74e-7 | ~0.5 B |
| Gaussian closure + exact Euler–Stein defect, e - a delta (a = 0.535) | 4.35e-7 | +2049 chains (diagnostic) |
| Gaussian closure + parameter-free midpoint (e + Lap e)/2 | 4.50e-7 | same |
| organisers' K=3 simple chain, numpy port (`whest/kprop3.py`), radial kappa4 off | 5.63e-7 | |
| same, radial kappa4 on (one scalar per layer) | **3.53e-8** | 3.4e12 raw FLOPs = 1.5 B f32 (dense legs, rank 2n per layer) |
| open-source chain (504aldo) | 2.13e-8 | 0.25 B |
| leaders (raw) | 1.1–1.5e-8 | 0.10–0.14 B |

- The port matches the torch reference to 1e-15 (per-layer means, covariances, (2,1) and (3) slices, radial core).
- Radial kappa4 scale: c*0.9 -> 8.9e-8, c*1.1 -> 2.7e-7, c*1.25 -> 4.4e-6 (c*1.0 = 3.5e-8 is the sharp optimum: no calibration gain there, and the chain is 2.5–8x sensitive to +-10%).
- Two-tier compression (young window w, shared basis k): (4,384 CP) **3.69e-8** (lossless to 5%); (4,256 CP) 5.2e-8; (2,384 CP) 7.9e-8; (4,128) 2.0e-7; (2,256) 3.2e-7;
  (2,128) 9.3e-7; (2,128)+collective directions 8.4e-7; (1,128) 2.0e-6; (0,128) 3.7e-6. Diagnostics (old sources kept dense):
  dropped (w=2) 3.7e-6; diagonal only 1.0e-6; scaled by 0.9: 5.1e-8; dropped (w=4) 7.4e-7. The old sources' effect is smooth
  in amplitude and brittle in shape: Hadamard readouts defeat Frobenius truncation (note, Sec. 6).
- Billed cost under flopscope (fnp port `whest/kprop3f.py`, float32, MSE identical to float64): (4,128) 0.862 B; (2,128) 0.592 B;
  (1,128) 0.439 B; (0,128) 0.275 B; float64 doubles these. The accurate configuration (4,384) is ~1 B as written, ~0.45 B with the
  identity-leg structure of the birth blocks (3 transported matrices and 4 readout products per source-age instead of 6 and 6) and Strassen.
- Exact Euler–Stein defect of the ported chain on net 0 (2049 chain runs on Modal): PENDING (job lapx2).

## 3. Compute

- Modal app `whest9` driven through the AWS relay `w9-1` (`infra/modal_relay.py`: deploy / batch / fetch / run / log).
  The sandbox proxy cannot carry gRPC; the relay can. Up to ~90 concurrent containers observed at cpu=1–4.
  The relay stopped itself twice from idleness while a batch was running; the batch wrapper now keeps it alive.
- flopscope rules, measured (`notes/compute/FLOPSCOPE.md`): float64 bills 2x; `einsum('ij,jk,lk->il', W, C_sym, W)`
  with `as_symmetric` bills 3 n^3 (float32) for W C W^T; `stats.norm.cdf` returns float64 (cast back or the whole
  chain is billed at the float64 rate — this happened in the first billing run); the floor 0.1 B is 12.8 n^3 per layer.

## 4. What this says about a competitive system

- The accuracy of the exact first-order chain (3.5e-8) comes from two things: the exact first-order term table (36
  diagrams, all index shifts) and the one radial kappa4 scalar per layer (the central variable). Our own chain had the
  first in a different form and lacked the second, which alone was worth 16x. Everything else of kappa4 is noise.
- The cost of that chain is the memory of kappa3: 12 n^3 per source-age with the identity-leg structure (24 as the
  reference writes it), 0.5-0.7 B over the network in float32. The floor is 12.8 n^3 per layer in total. A window
  of four ages with a 384-dimensional shared basis is lossless (3.7e-8) and still costs about 0.5 B; two ages with
  that basis give 7.9e-8; anything cheaper loses an order of magnitude. Frobenius compression of young sources does
  not work because the readouts are Hadamard products (note, Sec. 6).
- Therefore a 0.1 B entry at the leaders' 1.1-1.5e-8 is not a compressed first-order chain. The measured smoothness
  of the old sources' effect in amplitude (10% scale -> 1.5x) and the open-source chain's memoryless regeneration of
  kappa4 point to the same idea for kappa3: keep one or two dense ages, regenerate the rest from the current state
  with per-layer scalars fitted on the public networks against the dense chain's slices (see the shape diagnostic,
  `scripts/diag_oldshape.py`).
- The Euler-Stein merge is the one estimator-agnostic correction the theory provides; it is a factor 9 for the
  Gaussian closure at full scale and its value for the first-order chain is measured by the Modal job `lapx2`
  (PENDING). Making it affordable needs the symbolic second-order response (note, Thm 3.2; task left open).

## 5. Compute and hand-offs

- Modal: works through the AWS relay only (the sandbox proxy has no gRPC). `python infra/modal_relay.py batch JOB jobs/JOB.tsv --cpu N`;
  results in `s3://claude-whest-9292-97a992/results9/JOB/`, fetched with `python infra/fleet.py get JOB`. If the relay
  stopped itself, `python infra/fleet.py start w9-1` and `python infra/modal_relay.py fetch JOB` recovers results from
  the Modal volume. Concurrency reached 90 containers.
- AWS: spot and on-demand quotas are exhausted/limited (vCPU 200); only the relay runs. Quota increases must be
  requested in the console or the IAM user given `servicequotas:*` (see `notes/compute/SETUP.md`).
- Google Cloud: not set up; the browser-agent prompt is in `notes/compute/GCP_PROMPT.md`; hand the key back as the
  environment secret `GOOGLE_APPLICATION_CREDENTIALS_JSON`.
- Elicit: three reports were run earlier (`notes/stage9/elicit/` if present, else the scratchpad summaries fed the
  literature in the note).
- Organisers' reference code: cloned read-only at `/home/user/alignment-research-center/mlp_cumulant_propagation`
  (public); torch 2.14 CPU installed locally for the verification only. Nothing from it is in the repo except the
  term table listing.
