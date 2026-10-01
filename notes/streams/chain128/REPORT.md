# chain128 — does a better D21 interface give a lower final-layer MSE?

*Stream report, 2026-10-01. Status: in progress — interim numbers below, error-law and teacher-forcing runs queued.*

## Question
Build a dense K=3 moment chain at width 128 (full κ3(z), n³ = 2.1M entries; plus a fourth-cumulant closure) and test end
to end whether the closures identified by the oracle ladder (notes/competition-plan.md §3.1) translate into lower
final-layer MSE; measure the law MSE vs ε(D21) at width 128 and say what it implies at width 1024.

## Method (code in this directory)
- `chain.py` — dense chain. State per layer: μ, C, full κ3(z), and either X = κ4(z)_{aabc} (slices) or the full κ4(z)
  (n⁴, float32). ReLU step: Edgeworth operator E[F(z)] = E_G[exp(Σ κ_r/r! ∂^r) F] truncated at first order in κ3, κ4
  plus κ3²/2, on exact Gaussian expectations of derivatives of reluᵖ (truncated-moment closed forms = the paper's Hermite
  coefficients of powers, Prop S.3.8) with the Mehler series to all orders in C_ij (40 terms). Gives μ(a), C(a), the
  (3) and (2,1) slices of κ3(a) and the (4),(3,1),(2,2) slices of κ4(a). All-distinct κ3(a): either
  `oracle_k3` (`hermite_model` + `residual_basis` × `CLOSURE_COEF` / fitted table) or a generic **diagram engine**
  (`triple_engine`): all connected diagrams on three distinct vertices with ≤ 1 hyperedge (κ3 or κ4, any index
  pattern) and Gaussian C edges (multiplicity ≤ 6 without / ≤ 2 with a hyperedge). Exact multilinear transport.
- κ4 closures: `zero`; `mem` (κ4(a) restricted to its (4),(3,1),(2,2) slices, transported exactly to κ4(z)_{aabc});
  `dense` (full κ4 carried: all-distinct κ4(a) = Φ⁴κ4(z) + 16 Gaussian trees + 12 κ3-hyperedge-plus-edge diagrams;
  (2,1,1) slice of κ4(a) from the engine with vertex function relu² − 2m·relu; slices as in `mem`; n⁵ transport);
  `dense2` (+ the same-order κ3×κ3 and κ4-hyperedge-plus-edge diagrams); `atlas` (κ4(z)_{aabc} teacher-forced from a MC
  atlas at every layer: isolates the κ3 interface).
- `validate.py` (step fed the atlas's true z cumulants), `diag_state.py`, `fit_coefs.py`, `run_variants.py`, `summarize.py`.

Ground truth: `whest dataset bake --n-mlps 8 --n-samples 1e7 --width 128 --depth 16`, seeds 128000–128007; truth noise
floor avg_variance/N = 1.2e-8 (negligible here). Atlases: `moment_atlas_np.py --k3 --k4`, N = 5e5, MLPs 0–3 (one atlas
each; a second-seed atlas of MLP 0 for the noise floor is queued). Width-32 atlas (depth 8, N = 2e6, two seeds) for
step validation.

Reference implementation: ARC's `mlp_kprop` is not public (`git clone github.com/alignment-research-center/mlp_kprop`
asks for credentials; not on PyPI). The dense chain + diagram engine is the reference here.

## Results so far

### Step validation (width 32, atlas truth fed in; relative rms errors; noise = one-atlas MC noise)
Layer 0 (Gaussian input) reproduces the atlas to its MC noise (μ 1e-4, var 2e-4, C_off 1.3e-3, D3 4e-4, D21 3e-3), so
the Gaussian machinery (truncated moments, Mehler) is exact. Deeper layers have truncation errors above noise (μ 0.1–0.3 %,
var 0.2–1.6 %, C_off 1–3 %, D21(a) 2–7 %) at this very non-Gaussian width; order-1 Edgeworth is slightly better than
order 2 at width 32. D21(l+1) error of the transported κ3(a) with the true slices: memoryless 25–70 %, leading Wick 2–14 %,
oracle closure 0.7–13 %, **diagram engine 0.4–6 %**. The engine's Gaussian part equals `hermite_model` to 1e-16.

### End to end, width 128, depth 16 (final-layer MSE of post-activation means vs N = 1e7 truth; mean over MLPs)
Own fourth-cumulant closure, 8 MLPs:

| variant | κ3 all-distinct rule | κ4 closure | final MSE | layer 4 | layer 8 | layer 12 |
|---|---|---|---|---|---|---|
| A0 paper Alg. 2 | — (K=2) | — | 2.88e-4 | 8.4e-5 | 1.7e-4 | 2.2e-4 |
| A | — (K=2, Mehler to all orders) | — | 2.77e-4 | 7.9e-5 | 1.6e-4 | 2.1e-4 |
| M0 / M | slices only | zero / mem | 9.9e-5 / 8.6e-5 | | | |
| B0 / B | leading Wick | zero / mem | 4.2e-5 / 1.33e-4 | | | |
| C0 / C | oracle closure (leg-partition coefs) | zero / mem | 3.2e-5 / 8.4e-5 | 7.7e-6 / 1.6e-6 | 1.5e-5 / 1.4e-5 | 2.1e-5 / 4.2e-5 |
| D0 / D | fitted per-layer coefs (atlases 0,1) | zero / mem | 2.9e-5 / 8.4e-5 | | | |
| engine | diagram engine | zero / mem | 3.1e-5 / 8.0e-5 | | | |

Teacher-forced fourth cumulant (κ4(z)_{aabc} from the atlas at every layer), MLP 0:

| κ3 all-distinct rule | final MSE | mean over layers | D21 ε by layer (1, 4, 8, 12, 15), atlas noise included |
|---|---|---|---|
| slices only | 1.0e-4 | 4.4e-5 | |
| leading Wick | 8.6e-6 | 2.9e-6 | |
| oracle closure | 5.4e-6 | 1.6e-6 | 0.061 0.039 0.063 0.120 0.172 |
| fitted coefs | 2.0e-6 | 0.9e-6 | |
| diagram engine | **1.9e-6** | 0.7e-6 | 0.061 0.024 0.040 0.060 0.070 |

Reading so far. (i) With the fourth cumulant held at truth, the D21 interface ladder translates one-to-one into the
final layer: slices-only 1e-4 → Wick 8.6e-6 → closure 5.4e-6 → fitted / engine 2e-6 (50× from memoryless to engine).
(ii) Without it, every κ3 rule lands at 3e-5 (κ4 = 0) or 8e-5 (memoryless κ4): **at width 128 the fourth-cumulant
closure, not the κ3 interface, is the binding error.** The memoryless κ4 closure misses 16 % of κ4(z) at layer 1 and
diverges to 70–90 % by layer 10 (its D21 error follows: 6 % → 37 %). (iii) The dense κ4 closure tracks κ4(z) to 7–16 % up
to layer 6 and keeps the chain at the teacher-forced level through layer 8 (layer-8 MSE 9e-7 vs 1e-5 with κ4 = 0), then
degrades (K4 error 30–67 % at layers 7–15; final 8e-6 on MLP 0, 3.6e-5 on MLP 2).

(Tables regenerate from `results/raw`, `results/eps` with `summarize.py`.)
