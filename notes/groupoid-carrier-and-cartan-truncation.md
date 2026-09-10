# The network as a groupoid convolution algebra, and which truncations are admissible

*A neuron is a unit, not an edge. Then matrix multiplication is groupoid convolution, the activation is
a Cartan projection, the 0/1 code is the Cartan projection lattice, and the cumulant expansion is a path
sum with diagonal sources. This last identification rules out the obvious finite-resolution truncation,
and the ruling-out is measured.*

Companion to [`finite-resolution-question-algebra.md`](finite-resolution-question-algebra.md), 10 September 2026.
Measurements are on WhestBench Phase 2 networks (width 1024, depth 16) against the shipped `1e9`-sample
ground truth; the scripts are in [`whestbench/diagnostics/`](whestbench/diagnostics/).

---

## 0. Summary

The proposal under test: *view each neuron as an edge in a groupoid, the activation as a self-loop/unit,
and the weights as a function on the edges of the C\*-algebra.*

1. **The assignment has to be swapped, and the swap is forced by index bookkeeping.** Groupoid convolution
   sums over the *units*. Matrix multiplication sums over *neurons*. So the neuron is the unit and the weight
   entry is the arrow. With neurons as edges (vertices = layers), every composable pair factors uniquely,
   convolution degenerates to a pointwise product, and no contraction happens at all.
2. **Everything the proposal wanted survives the swap, and more sharply.** With neurons as units, the
   diagonal subalgebra `C_0(X)` — the span of the units, i.e. the self-loops — is exactly where the
   activation lives. `(C*(G), C_0(X))` is a Cartan pair, and the 0/1 activation code is its projection
   lattice. That is a second, independent motivation for the code, on a different algebra from the one in
   the companion note.
3. **The network's nonlinearity is the state-dependence of a Cartan projection.** ReLU is `z ↦ P(z) z`
   with `P` a diagonal projection. The algebraic operation is multiplication by a Cartan element; only the
   *dependence of `P` on the point* is nonlinear. Estimation is then: average an alternating product of
   arrows and state-dependent Cartan projections against a Gaussian state.
4. **The cumulant expansion is a groupoid path sum with Cartan-supported sources**, graded by path length.
   That grading is what my estimator computes, so the picture is descriptive of working code, not decorative.
5. **It makes a falsifiable prediction, and the prediction is false.** Connes–van Suijlekom spectral
   truncation says: truncate the transport to its top spectral modes. Measured, that loses more than
   depth gains — rank 128 is 2.4× worse than exact, and rank 16 is no worse than rank 64, so it is a bias,
   not a resolution effect. The groupoid picture explains the failure: a spectral projection of the
   convolution part does not commute with `C_0(X)`, so it does not preserve the Cartan pair that the
   sources live in.
6. **The measurement that does support the framework**: the correction that matters is the *off-diagonal*
   two-point cumulant, not the per-neuron diagonal. Restricting the two-point channel to the nearest
   source costs 3.1× in final-layer MSE; restricting the diagonal channel the same way costs 6%.

---

## 1. Why the neuron must be the unit

For a discrete groupoid `G` with unit space `X = G^(0)`, convolution is

```
(f1 * f2)(γ) = Σ_{αβ = γ}  f1(α) f2(β),
```

a sum over factorisations of `γ`. The index summed over is the intermediate **object**. Now compare a layer:

```
(W_{l+1} W_l)_{ik} = Σ_j  W_l[i,j] W_{l+1}[j,k],
```

summed over the neuron `j`. For convolution to reproduce this, `j` must be an object of the groupoid.

Take the proposal literally instead — vertices are layers `{0,…,L}`, edges are neurons, `e^l_j : l → l+1`.
A composable pair `(e^l_i, e^{l+1}_j)` has exactly one factorisation, so

```
(f1 * f2)(e^l_i e^{l+1}_j) = f1(e^l_i) f2(e^{l+1}_j),
```

a pointwise product on length-two paths. Nothing is contracted, because the intermediate object is the
*layer*, and each layer is a single object. The weight `W_l[i,j]`, indexed by a composable pair, is then a
2-cochain `G^(2) → C` — which is the slot a twist `σ` occupies in `C*(G, σ)`. That is a real structure, but
a twist *multiplies* the convolution; it does not perform it. So the literal reading cannot express a layer.

## 2. The corrected carrier

Let the unit space be the neurons themselves,

```
X  =  ⊔_{l=0}^{L} V_l ,      |V_l| = n ,
```

and `G = X × X` the pair groupoid, so `C*(G) = M_N`, `N = (L+1)n`, with convolution equal to matrix
multiplication and `C_0(X)` the diagonal. The network's weights are the single block-superdiagonal element

```
W  =  Σ_l Σ_{ij}  W_l[i,j] · δ_{(i∈V_l → j∈V_{l+1})}  ∈  C*(G),
```

and the `(V_0, V_L)` block of the convolution power `W^{*L}` is exactly `W_L ⋯ W_1`. Path length is the
block grading. (An étale sub-groupoid generated only by the layer arrows gives the same thing more
economically; the pair groupoid is the shortest correct statement.)

`(M_N, D_N)` is a Cartan pair in Renault's sense — `D_N` is a MASA, its normalisers generate, and
`E(a) = diag(a)` is a faithful conditional expectation — and the groupoid it reconstructs is the pair
groupoid on the neurons. So the user's "activation as self-loop/unit" is right on the nose once the swap is
made: **the self-loops span the diagonal, and the activation is diagonal.**

## 3. The activation, precisely

ReLU is not an element of `C*(G)`; it is not linear. What is in the algebra is the gate. Writing
`z(x)` for the pre-activations at some layer,

```
ReLU(z)  =  P(z) · z ,      P(z) = diag( 1[z_j > 0] )  ∈  projections of C_0(X).
```

So the layer map is *multiplication by a Cartan projection*, and the whole nonlinearity is the dependence
of that projection on the point `x`. The forward pass is an alternating word

```
W · P · W · P · ⋯ · W ,
```

off-diagonal, diagonal, off-diagonal, …, and the challenge target is the Gaussian-state expectation of
such a word. This is the same alternation my estimator runs: the linear-response transport it
back-propagates is literally `W_s D^{(s)} W_{s+1} D^{(s+1)} ⋯`, with `D = E[P] = diag(Φ(t_j))`.

**The 0/1 code, again.** The companion note motivates the code on the *question* algebra (functions of the
input), where it is the smallest subalgebra resolving the sign of every parameter-labelled question. Here it
appears on a different algebra for a different reason: the code is the projection lattice of the Cartan
subalgebra, `P_w ∈ C_0(X)`, and all the noncommutativity of the network sits in the arrows. Two
independent derivations of the same object on two different algebras is a better sign than either alone.

## 4. The cumulant expansion is a path sum with diagonal sources

The non-Gaussian correction my estimator carries has the form

```
κ_3(z_{l,a})  =  Σ_{s<l}  Σ_{β ∈ V_s}  ( source at unit β )  ·  ( transport along paths β → a ),
```

where the source coefficients are

```
p_β = φ(t_β)/σ_β ,     q_β = 2 Ψ_β (1 − Φ_β) ,     κ_{3,β} = third cumulant of ReLU at unit β ,
```

all **functions on the unit space** — Cartan elements — and the transport is convolution. So the expansion
is graded by (source unit, path length), which is exactly the groupoid grading of §2. The tensor itself
transports linearly,

```
T^{(l+1)}  =  (Λ_l ⊗ Λ_l ⊗ Λ_l)[ T^{(l)} ]  +  S^{(l)} ,
```

with `S^{(l)}` the freshly generated Cartan-supported source; the factorised form is closed under transport
precisely *because* the source is diagonal, one set of factors per source layer. That is why the cost is
`O(K)` matmuls per layer for `K` retained source layers, and it is the reason a compression is wanted.

## 5. The prediction, and the measurement that kills it

The companion note's finite-resolution machinery is Connes–van Suijlekom spectral truncation: replace `a`
by `PaP` for a spectral projection `P`. Applied here the prescription is unambiguous — truncate the far
transports to their top spectral modes, keeping the near ones sharp. The transports are strongly
contracting (the map from layer 8 to layer 15 needs 23 singular directions for half its energy, 85 for 90%),
so this looks like it should be nearly free.

Measured on two networks, keeping the newest 2 source layers exact and truncating the rest, final-layer MSE:

| source depth | transport rank | net 0 | net 1 |
|---|---|---|---|
| 1 | exact | 9.52e-7 | 9.35e-7 |
| 8 | exact | 2.19e-7 | 1.56e-7 |
| 8 | 128 | 5.54e-7 | 4.82e-7 |
| 8 | 64 | 5.92e-7 | 6.21e-7 |
| 8 | 32 | 6.23e-7 | 6.85e-7 |
| 8 | 16 (K₀=3) | 4.57e-7 | 5.12e-7 |

Truncation costs 2.4–3.1× against exact, and rank 16 is *not* worse than rank 64 — the damage is a bias
that does not decrease as resolution is added. So this is not a resolution trade-off at all.

**Why, in the framework's own terms.** The sources are supported on the unit space and weighted diagonally.
A spectral projection `P` of the convolution part does not commute with `C_0(X)`, so `a ↦ PaP` is not a
conditional expectation onto a subalgebra containing the diagonal: it mixes the unit index and destroys
exactly the Cartan support the sources live in. The truncation is inadmissible for this algebra, and the
error it makes is structural rather than small.

The admissible finite-resolution moves here are the ones compatible with the Cartan pair — conditional
expectations onto intermediate subalgebras containing `C_0(X)`, corners `1_Y A 1_Y` for a sub-unit-space
`Y ⊂ X`, or coarse-grainings of the units (groupoid quotients). Note this does not contradict the
companion note: **the two truncations live on different algebras.** Chaos truncation on the *question*
algebra retains the vacuum and is exact on the target (companion note §4.1); spectral truncation on the
*neuron* algebra breaks a Cartan pair. Conflating them is what the framework is for avoiding.

## 6. What the framework got right: noncommutativity is where the error is

The third cumulant has a diagonal part `κ(z_a,z_a,z_a)` and an off-diagonal part `K21(a,b) = κ(z_a,z_a,z_b)`.
The framework's central claim is that the obstruction is noncommutativity, so the off-diagonal part should
carry the correction. Measured, mean final-layer MSE over three networks, varying the depth of each channel
independently:

| diagonal depth | two-point depth | mean final MSE |
|---|---|---|
| 4 | 1 | 7.77e-7 |
| 8 | 1 | 7.29e-7 |
| 4 | full | 3.95e-7 |
| 8 | full | 2.32e-7 |

Deepening the *diagonal* channel from 4 to 8 sources buys 6% when the two-point channel is shallow.
Deepening it with the two-point channel at full depth buys 41%. Restricting the two-point channel to the
nearest source, at fixed diagonal depth 8, costs 3.1×. The depth of the path sum matters almost entirely
through its off-diagonal component.

Consistently, `K21` is close to rank one: at layers 13–16 the leading singular component carries 89–93% of
its energy, and its right singular vector correlates 0.98 with the propagated mean direction. At resolution
one, the off-diagonal cumulant is a single common mode.

## 7. Standing measurement and limits

The estimator built on this expansion scores, on the 16-network public shard under the official harness
(`whest run --dataset … --split mini`), zero failures:

| | |
|---|---|
| adjusted final-layer score | 9.61e-8 |
| raw final-layer MSE | 9.61e-7 |
| compute utilisation | 9.79% |

The Phase 2 leader stands at 2.8e-9 adjusted (1.93e-8 raw, 14.7% compute), so this is roughly 34× off the
frontier and the picture in this note has not closed that gap.

What the note does not claim: that the groupoid carrier produces a faster estimator. It does not, so far.
Its contributions are one correction (the unit/edge swap), one structural identification (activation as a
state-dependent Cartan projection, code as the Cartan projection lattice), one falsified prediction with a
principled explanation (spectral truncation is inadmissible on this algebra), and one confirmed one
(the off-diagonal channel dominates). Under the challenge's score, `MSE × max(0.1, C/B)`, the path-sum
depth that the picture organises is bought at a linear price in compute, and the optimum sits at 4–6
source layers — worth about 1.6× over the standing number, not the 34× that would matter.
