# Heisenberg–Duhamel stream — VERDICT (1 Oct 2026)

**Principle.** The Heisenberg picture on Gaussian space. The readout is pulled back through the Koopman operators of the arrows; a Gaussian reference chain (closure) is pushed forward. Exact identity (Theorem 1, verified at toy scale): truth − reference = Σ_l (ρ_l − ν_l)[g_l]. Each term is the local non-Gaussianity created by one arrow acting on a Gaussian, paired with the *exact* pulled-back observable g_l. Any model ĝ of g leaves a **bilinear** remainder Σ_l (ρ_l − ν_l)[g_l − ĝ_l] (Theorem 2). Stein/Edgeworth turns each pairing into cumulant sources contracted with Gaussian-averaged derivatives of g_l.

**Estimator (first-order HD).** Gaussian reference (m, C), plus κ₃ sources at every layer evaluated at the *source* layer on pulled-back direction matrices U = W_{l+1}Φ⋯W_k:
- diagonal + star diagrams + exact coincident planes;
- all ages; no forward n³ tensors and no low-rank assumption on old content;
- Stein injection of the diagonal κ₃ and (2,1) slice into the next layer's mean and covariance.

Cost: 7 products per (source, target) pair, ≈ 485 u (0.47 B) for all ages, ≈ 375 u for 8 sources. `hdpull.py` equals the full-tensor prototype to 1e-15.

**Measured directly at n = 1024** (bench w1024_d16, 6 networks, paired, truth noise subtracted):

| | raw MSE | adjusted |
|---|---|---|
| Gaussian closure (reference) | 4.10e-6 | ≈ 4.1e-7 (0.1 floor) |
| first-order HD, all ages | **4.23e-7 ± 2.1e-8** | **≈ 2.0e-7** at 0.47 B |
| + exact coincident κ₃ transport (re-birth) | 4.02e-7 ± 2.8e-8 | ≈ 1.9e-7 |
| HD, 8 sources | 4.90e-7 | ≈ 1.8e-7 at 0.37 B |

The gain is ≈ 10× in raw MSE, constant from width 256 to 1024 (HD ∝ n^{-1.86}, reference ∝ n^{-1.83}). The earlier "30×" was against a mis-calibrated reference and is withdrawn. **Not competitive** with the bar (adjusted 1.6e-9); it sits on the shared full-history ceiling of 3–4e-7.

**Why it stops there (charged to the realisation, not the principle).**
1. The remaining error is covariance drift from joint κ₄ slices κ₄(pppq) and κ₄(ppqq). In the width-64 oracle, true κ₃ and κ₄ on the own chain gain 7–20× over first order. Age-0 or single-site κ₄ is the wrong object; the needed κ₄ is old, multi-site content.
2. Second order cannot help until κ₃ transport is exact. Adding the *true* κ₄ on top of computed first order makes things worse (gain 0.27–0.83 at widths 64 and 128). The transported κ₃ slices are 20–30 % wrong from the first transport step, while the sources are exact.
3. The exact zero-correlation coincident response removes only a small part of that error (layer-2 slice 22 → 19 %; −5 % MSE at 1024). The rest is the O(ρ) part of the one-step response (joint gate law, δ-insertions with correlation edges); mean |ρ| is still 0.1–0.2 at n = 256–1024.

**What survives.**
- The exact Duhamel accounting and the bilinear remainder: a clean way to say which object each layer's error pairs with.
- The source-layer pull-back: full-rank old content at n³ per (source, target) pair, the backward answer to "no forward low-rank state below 0.3 n".

**Deciding next experiment** (for whoever continues). Implement the resummed one-step response of κ₃ through a ReLU at non-zero correlation (pair gate covariances, foundations-facts F6.8) as a pull-back re-birth. Require transported-slice error ≤ 3–5 % at width 64. Then re-run the true-κ₄ injection: proceed to old κ₄ transport only if it then gains ≥ 5×.

Files: DESIGN.md (derivation, §6 status), RESULTS.md (R0–R12), hd.py (full-tensor prototype), hdpull.py (n × n implementation), oracle*.py, decide.py, k3err*.py.
