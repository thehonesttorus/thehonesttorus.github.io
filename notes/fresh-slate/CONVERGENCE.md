# What six principle-driven designs converged on (1 Oct 2026, 19:15 UTC)

Six fresh-slate design streams started from six different principles (faces and barycentres; Bethe/cavity; matchings and signings; tropical geometry and temperature; the Heisenberg picture with Dirichlet forms; Markov networks and CMI). Within three hours they arrived, independently, at the same structure and the same two binding facts. That convergence is a finding about the problem, and it frames the next round. Sources: notes/fresh-slate/designs/*/DESIGN.md and RESULTS.md.

## The common structure

1. **A Gaussian (or Gaussian-copula) reference propagated exactly is the zeroth order of every design.** Tropical: every realisation of the zero-temperature skeleton reduces to it (no small parameter at L = 16; a 10²³-fold P − Q cancellation at n = 1024). Signings: loop diagrams with fresh-weight legs vanish for n ≥ 96, so one-step sums are tree sums around the Gaussian product. Bethe: the product of node beliefs is the tree reference. Heisenberg: the reference chain of the Duhamel identity. The small parameter of every correction is n^{-1/2}.
2. **Binding fact (i): the covariance must be propagated with its non-Gaussian corrections.** Heisenberg's oracle at width 64: true mean and covariance plus true κ3, κ4 → final MSE 2.8e-6 / 4.3e-6, while its own covariance chain with the same injection stays at 4e-5 – 2e-4; after first order, the residual is covariance drift caused by the joint fourth-cumulant slices κ4(pppq), κ4(ppqq). Markov: its first failure was variance drift, fixed only by putting non-Gaussian structure into the covariance propagation, and the two-site κ4 tree diagrams are leading order for κ4 because they add coherently. Markov's remaining error at depth sits in the layer means of nearly-always-on neurons.
3. **Binding fact (ii): third-order content generated many layers earlier matters at the readout.** Heisenberg: ages ≤ 1 give almost nothing, all ages give 18–23× at width 64, 8 ages keep 92 % of the gain. Markov: final MSE falls with the depth window only when all history is carried; the depth-cut CMI is large and the top n/4 transported directions recover only half of it. Faces: the final-layer third cumulant is born roughly uniformly over all 15 earlier layers; the most recent facets contribute about nothing net; old births do not accumulate because later facets fold them away; non-Gaussianity is carried by facets (codimension 1), not faces. The earlier moment-chain measurements say the same at n = 1024: old content is ≈ 40 % of the (2,1) slice, orthogonal to anything born at the current layer, not compressible below ≈ 0.3 n modes, and not absorbable into a renormalised local closure.
4. **Width laws are steep but pre-asymptotic below n ≈ 256.** First-order designs fall as n^{-2} to n^{-3.2} between 128 and 256; the Gaussian closure falls as n^{-0.82} between 64 and 128 but ≈ n^{-2} between 128 and 1024 (measured: raw 4.3e-6 at 1024 on the bench). Projections from small widths are unreliable; direct width-1024 runs on the bench decide.

## Where the designs stand (projections, being replaced by direct width-1024 runs)

| design | best variant | projected adjusted at 1024 | cost | status |
|---|---|---|---|---|
| heisenberg | first-order Heisenberg–Duhamel, all ages | ≈ 1e-7 (raw 4–7e-7) | ≈ 0.26 B | strongest; second order bounded 7–20× by oracle |
| markov | MKV-2 full history | ≈ 1e-7 (range 2e-8 – 5e-7) | ≈ 0.49 B | closed: not competitive |
| signings | copula v1 | ≈ 3e-7 | — | closed |
| faces | 'mem' / 'lin21' | 1.6e-6 / 4e-6 (from small-width fits) | 0.1 / 0.6 B | direct 1024 run in flight; facet-conditional covariance test pending |
| tropical | TCT-0 + T-extrapolation | ≈ 4e-7 | 0.02–0.04 B | closed: not competitive |
| bethe | v1/v2 node + pair beliefs | pending | ≈ 0.10 B (+ old triples) | v2 in progress |

For comparison: the public moment chain reaches raw 2.1e-8 at 0.25 B (adjusted 5.4e-9); the leaders raw 1.5e-8 at 0.11 B (1.6e-9) and raw 1.14e-8 at 0.15 B; an oracle with exact per-neuron moments to fourth order at every layer reaches raw 1.17e-9. The best fresh design is 30–60× short of the bar.

## The open question for the next round

Every principle so far has been realised as "a Gaussian reference plus carried corrections", and the corrections that matter are (i) non-Gaussian covariance propagation through each ReLU, driven by joint κ3 (2,1) and κ4 (3,1)/(2,2)/(2,1,1) slices, and (ii) third-order content of every age. Carrying (ii) explicitly costs O(L² n³) (pairs of source and target layers). The question is whether a theoretical unlock gives a representation in which (i) and (ii) are carried at O(L n³) with ≲ 10 dense products per layer, or makes part of them unnecessary, at the accuracy the bar needs (per-neuron mean errors ≈ 1e-4 rms at n = 1024).
