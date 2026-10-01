# The MLP bridge: dictionary v1 and the first targeted experiments

### Reading an actual ReLU network into the barycentre-arrow objects, and testing the theory's dichotomies on it

*Prong 2 of the [research programme](research-program.md). Inputs from the other prongs used here: the exact statements being tested, quoted from Note 1 ([conditional arrow algebra](conditional-arrow-algebra.md)) and Note 2 ([simplicial complex as decomposition](simplicial-complex-as-decomposition.md)); from prong 3 ([local-to-global unlocks](local-to-global-unlocks.md)) only the list of quantities it asked to have measured. Script: [`experiments/mlp_bridge.py`](experiments/mlp_bridge.py), pure numpy, 3 seconds per run on the stand-ins. Status labels: **Dictionary** (a modelling choice), **Measured** (a number from a run), **Reading** (interpretation), **Naive** (a known weakness of the dictionary, to be charged first when a result is negative, programme rule P2.1).*

---

## 0. What a bridge is, and what it is not

The theory (Notes 1–2) is conditional on a quiver $\Lambda$ of arrows between faces and a cocycle $F$ on the arrows; it says nothing about how either is read off a network (Note 1 §2.1, "the one open modelling decision"). A bridge is a **dictionary**: an explicit rule that reads a trained network and an input law into $(\Lambda,F,\text{states})$. The experiments then test the theory's *dichotomies* (sufficient / not, coboundary / not, groupoid composition / frame composition, trivial / non-trivial nerve) on the realization produced by the dictionary. A dichotomy coming out one way is a fact about the realization under that dictionary; a result that looks "negative" is first a statement that the dictionary is too naive (rule P2.1), because the theory's statements are theorems about the abstract objects and cannot fail.

Everything below respects the standing constraints: objects are faces (abstract barycentres), never activation vectors; an arrow is the ReLU that produced the source face together with the following linear map; weights are not elements of the arrow algebra; no coordinates, metric or walk is assumed.

---

## 1. Dictionary v1

| abstract object | read as | status |
|---|---|---|
| face $\sigma\in K_\ell$ at hidden layer $\ell$ | the set of units with nonzero ReLU output at layer $\ell$; $d(\sigma)=|\sigma|$ | **Dictionary** |
| face at layer $0$ | one trivial face $*$ (option: the nonzero pattern of the input vector) | **Dictionary**, **Naive N5** |
| face at the output layer | the argmax class, a vertex of the output simplex (or the sign of a single logit) | **Dictionary** |
| arrow $\gamma:\sigma\to\tau$ | an observed transition: the ReLU that produced $\sigma$, then $W_{\ell+1}$, landing in $\tau$ | **Dictionary** |
| multiplicity of $\gamma$ | its empirical count | **Dictionary** |
| state on layer $0$ (origin weights $w_{\sigma_0}$, Note 1 §5.3) | the input law: the training law for the stand-ins, a user file otherwise | **Dictionary** |
| cocycle $F(\gamma)$ | $-\log\hat P(\gamma\mid s(\gamma))$, the empirical Markov kernel, i.e. the Doob-normalised walk of Note 1 §5.3 with $\beta=1$ | **Dictionary**, **Naive N3** |
| restriction to a sub-frame (Note 2 §2) | `--restrict k`: faces are read on a random subset of $k$ units per layer | **Dictionary**, **Naive N6** |

**Naive points, stated in advance.**
- **N1** Faces are full activation patterns, so the resolution of the objects is fixed by the width; in a wide network almost every input has its own face and no statistics are possible. The theory's answer is restriction to a sub-frame (Note 2 §2), but v1 chooses the sub-frame at random (N6).
- **N2** The state is the input law of the task. It enters only as origin weights, as Note 1 §5.3 says it should, but nothing in v1 distinguishes inputs the network was trained on from inputs it was not.
- **N3** $F$ is read from data, not from the weights. The weight-only readings of Note 1 §2.2 (which arrows exist from the restricted map $W[\tau,\sigma]$; a log-gain; a measure on outgoing arrows) are not implemented. Consequently $F$ is a Markov kernel by construction, and the question "is the law on histories of Gibbs form for *some* cocycle" becomes "is the face process Markov" (see §3).
- **N4** The output face is a single vertex; the multi-vertex output faces (sets of classes with positive logit) are not used.
- **N5** Layer $0$ is trivial for inputs that can be negative, so the first arrow carries no ReLU.
- **N6** Random sub-frames. The theory suggests the restriction should be the *minimal sufficient subalgebra* (Note 1 §6, Jenčová–Petz), not a random one.

---

## 2. The experiments and the statement each one targets

| id | what is computed | theory statement targeted |
|---|---|---|
| E1 | faces, arrows, mean $d$, fraction of arrows with $|d(\tau)-d(\sigma)|>1$, heights $h_\tau$ (distinct histories arriving at $\tau$), fraction of faces with a unique ancestor history | the quiver of Note 1 §2.1; "dimension can jump"; heights and the coherence excess of Note 2 §6–7 |
| E2 | $I(\sigma_{\ell+1};\sigma_{\ell-j}\mid\sigma_\ell)$ for each lag $j$, and with the whole past, minus a permutation null; as a fraction of $H(\sigma_{\ell+1}\mid\sigma_\ell)$ | sufficiency of the barycentre, Note 1 §6: all of these vanish iff the current face carries the information of the history |
| E3 | weighted least-squares projection of $F$ on coboundaries $\delta U$, residual fraction $\|F-\delta U\|/\|F-\bar F\|$ (globally and per layer); $R^2$ of the accumulated action $F(\mu)$ explained by the endpoint face; the centrality defect $\mathbb E_\tau[\log h_\tau]-H(\mu\mid\tau)$ | the coboundary criterion of Note 1 §6 and §10 item 2: face sufficient $\iff F=\delta U\iff$ the measure is central (tail-invariant). Exact identity used: $F=\delta U+\text{const}$ on a layer $\iff\hat P(\tau\mid\sigma)=g(\sigma)h(\tau)$ on the support of the quiver |
| E4 | over pairs of faces in one layer: correlation between the overlap cosine $G_{\sigma\sigma'}=|\sigma\cap\sigma'|/\sqrt{|\sigma||\sigma'|}$ and the fidelity $\sum_\tau\sqrt{\hat P(\tau|\sigma)\hat P(\tau|\sigma')}$ of their outgoing kernels, against a null that re-pairs patterns and kernels; and the correlation between $G$ before a step and the expected $G$ of the targets after it | the composition fork of Note 2 §3 and §7: the groupoid composition is blind to overlaps, the frame (operator-system) composition weights by $G$. A positive relation means the kernel is $G$-continuous and the tolerance relation of Connes–van Suijlekom carries information |
| E6 | Betti numbers $b_0,b_1$ of the 2-skeleton of each layer's nerve (vertices $=$ units, faces $=$ co-active sets) | Note 2 §7, item 1: invariants of the layer complexes from patterns alone |
| self-check | E2 and E3 re-run on histories sampled from the fitted Markov kernels (same quiver, no memory) | calibration: which quantities measure memory and which measure the kernel |

E5 (propagation numbers against observed memory depth) is deferred: the memory profile of E2 is its empirical side, the operator-system side needs prong 1's finite-resolution definitions on a realized quiver.

---

## 3. Measured: stand-ins

Architecture $2\to24\to24\to24\to1$, inputs $\mathcal N(0,1.3^2I)$, 40 000 samples. "Trained" is a numpy Adam run (4000 steps) on the synthetic rule $\operatorname{sign}(\sin2.2x_1+\sin2.2x_2+0.5x_1x_2)$, test accuracy 0.98–0.99; "random" is the same architecture at initialisation (same seed, same draw).

### 3.1 The quiver (E1), trained, seed 0

| transition | faces at target | arrows | mean $d$ / width | $|\Delta d|>1$ | heights mean / max | faces with unique history |
|---|---|---|---|---|---|---|
| $*\to$ hidden 1 | 236 | 236 | 13.4 / 24 | – | 1.0 / 1 | 1.00 |
| hidden 1 $\to$ 2 | 357 | 1097 | 11.9 / 24 | 0.62 | 3.1 / 40 | 0.38 |
| hidden 2 $\to$ 3 | 370 | 1342 | 11.1 / 24 | 0.54 | 6.6 / 79 | 0.22 |
| hidden 3 $\to$ class | 2 | 415 | 1 | – | 1278 / 1529 | 0.00 |

**Reading.** Dimension jumps by more than one on more than half of the arrows (0.52–0.82 across all runs): the "projection information" is genuinely not a step of $\pm1$. Heights grow with depth and the fraction of faces with a unique ancestor history falls from 1 to about 0.2–0.4: the coherence excess of Note 2 §6 (faces reachable from a common past but not on a common history) is already substantial at depth 3. The random network has fewer faces (218/210/237) but the same profile.

### 3.2 Memory (E2) and the class of $F$ (E3)

Excess conditional information, in nats and as a fraction of $H(\text{next}\mid\text{present})$; null-subtracted. "Residual" is the per-layer coboundary residual fraction of E3; $R^2$ is the share of the variance of the accumulated action explained by the endpoint face; "defect" is the centrality defect in nats.

| run | $I(\sigma_3;\sigma_1\mid\sigma_2)$ | $I(\text{class};\sigma_2\mid\sigma_3)$ | residual h1$\to$h2 / h2$\to$h3 / h3$\to$class | $R^2$ endpoint h2 / h3 / class | defect h2 / h3 |
|---|---|---|---|---|---|
| trained, seed 0 | 0.377 (30 %) | 0.065 (40 %) | 0.46 / 0.50 / 0.82 | 0.64 / 0.47 / 0.02 | 0.42 / 0.59 |
| trained, seed 1 | 0.345 (25 %) | 0.047 (35 %) | 0.42 / 0.47 / 0.73 | – / 0.51 / – | – / 0.58 |
| trained, seed 2 | 0.354 (26 %) | 0.037 (38 %) | 0.46 / 0.48 / 0.72 | – / 0.50 / – | – / 0.56 |
| random, seed 0 | 0.371 (31 %) | (class constant) | 0.46 / 0.48 / – | 0.71 / 0.60 / 0.00 | 0.45 / 0.62 |
| random, seed 1 | 0.457 (34 %) | 0.035 (39 %) | 0.45 / 0.42 / 0.79 | – / 0.64 / – | – / 0.62 |
| random, seed 2 | 0.346 (34 %) | 0.078 (36 %) | 0.42 / 0.43 / 0.79 | – / 0.62 / – | – / 0.77 |
| trained s0, restrict 12 | 0.433 (37 %) | 0.094 (38 %) | 0.51 / 0.50 / 0.86 | 0.52 / 0.46 / 0.13 | 0.46 / 0.64 |
| trained s0, restrict 6 | 0.372 (32 %) | 0.257 (50 %) | 0.53 / 0.56 / 0.80 | 0.46 / 0.43 / 0.00 | 0.57 / 0.75 |
| **Markov surrogate** of trained s0 | **0.001 (0 %)** | **0.000 (0 %)** | 0.46 / 0.50 / 0.81 | 0.63 / 0.43 / 0.01 | 0.41 / 0.67 |

Deeper stand-in $2\to32^{\times4}\to1$ (seed 1, 60 000 samples): memory 28 % at hidden 2$\to$3; at hidden 3$\to$4 lag 1 gives 31 %, lag 2 gives 34 %, the whole past 0.53 nats of 1.25; residuals 0.43 / 0.41 / 0.50 / 0.65; $R^2$ 0.67 / 0.50 / 0.43 / 0.00.

**Reading.** Two different failures of sufficiency, which the self-check separates cleanly.

1. **The realized law on histories is not Markov in the faces.** About a quarter to a third of the uncertainty about the next face, given the current face, is removed by the previous face, and in the deeper network the lag-2 face adds as much again. The surrogate, which has the same quiver and the same arrows with the same multiplicities but no memory, measures 0.000–0.001: the estimator is null-calibrated and the memory is real. In the theory's terms the empirical law is **outside the Gibbs class** of Note 1 §5.2: it is not $w_{\sigma_0}e^{-F(\mu)}$ for any cocycle $F$, because any such law is Markov. The barycentre is not sufficient, and the reason is not a non-trivial class of a cocycle but history dependence that no cocycle on arrows can express.
2. **Even within the Gibbs class the cocycle is not a coboundary.** Restricting attention to the Markov part (the fitted kernel, which is the law of the surrogate), the per-layer residual of $F$ against coboundaries is 0.42–0.56 on hidden transitions and 0.72–0.86 on the class transition: $\hat P(\tau\mid\sigma)$ is far from product form $g(\sigma)h(\tau)$ on the support of the quiver. The endpoint face explains about 45–65 % of the variance of the accumulated action at depth 2–3 and essentially none of it at the class, and the centrality defect grows with depth. These three numbers are unchanged by the surrogate (0.46/0.50 vs 0.46/0.50; defects 0.41/0.67 vs 0.42/0.59), which confirms they are properties of the kernel, as Note 1 §6 says they should be, not of the memory.

**Training is invisible to dictionary v1.** Trained and random networks of the same architecture give the same memory fraction (25–30 % vs 31–34 %), the same residuals (0.42–0.50) and the same $G$-continuity (§3.3). Three seeds are not enough to call the 5-point difference in memory fraction a training effect. Under rule P2.1 this is read as: *the dictionary does not see what training does*. Faces, arrows and the input law are all the same kind of object before and after training; what changes under training must live in a finer reading (§4).

**Random restriction does not restore sufficiency.** Reading faces on 12 or 6 of the 24 units leaves the memory fraction at 32–50 % and raises the residuals. Pinching to a random sub-frame coarsens the present without coarsening the past, so the past explains more, not less. This is the expected behaviour of a *wrong* restriction and is what N6 predicted; the right restriction is the one the theory names (§4, D2).

### 3.3 The composition fork (E4)

| run | layer | faces used | corr$(G,\text{Fid})$ | null | corr$(G_{\text{in}},\mathbb E\,G_{\text{out}})$ | mean Fid, $G<0.5$ / $0.5$–$0.8$ / $\ge0.8$ |
|---|---|---|---|---|---|---|
| trained s0 | hidden 1 | 134 | +0.40 | 0.00 ± 0.02 | +0.81 | 0.00 / 0.02 / 0.20 |
| trained s0 | hidden 2 | 227 | +0.38 | 0.00 ± 0.00 | +0.67 | 0.00 / 0.01 / 0.18 |
| trained s0 | hidden 3 | 186 | +0.18 | −0.01 ± 0.01 | +0.17 | 0.50 / 0.55 / 0.78 |
| random s0 | hidden 1 | 99 | +0.40 | −0.01 ± 0.01 | +0.83 | 0.00 / 0.00 / 0.19 |
| random s0 | hidden 2 | 116 | +0.42 | 0.00 ± 0.01 | +0.66 | 0.00 / 0.00 / 0.15 |
| trained seeds 1, 2; random seeds 1, 2 | hidden 2 | 196, 239; 108, 106 | +0.35, +0.36; +0.41, +0.44 | ≈ 0 | +0.57, +0.36; +0.75, +0.63 | – |
| deep s1 | hidden 1–4 | 267–316 | +0.32, +0.33, +0.35, +0.39 | ≈ 0 | +0.77, +0.55, +0.63, +0.37 | – |

**Reading.** The outgoing kernel of a face is $G$-continuous: faces sharing most of their vertices have overlapping next-face laws (fidelity 0.15–0.20 at $G\ge0.8$ against 0.00 below $G<0.5$), and the correlation is six to a hundred null standard deviations from zero in every layer of every run. The arrows also *transport* the overlap: the expected overlap of the targets correlates at 0.4–0.8 with the overlap of the sources. Under dictionary v1 the realization therefore behaves like the **frame picture** of Note 2 §3, where arrows between overlapping faces compose with the overlap cosine, rather than the orthogonal groupoid picture, which discards $G$. This is the one dichotomy on which v1 gives a definite answer, and the answer is the same for trained and random networks and at every width tried.

### 3.4 Nerves (E6)

All layer nerves of the 24- and 32-unit stand-ins have $b_0=1$, $b_1=0$, and are within a few edges of the full simplex (275 of 276 pairs of units co-active somewhere at hidden 1). A 6-unit random sub-frame produced one layer with $b_1=1$.

**Reading.** At these widths the topology of the layer complexes is trivial; every pair and almost every triple of units is jointly active for some input. This is the empirical side of Note 1 §10: for a finite network the NCG topology is trivial and the content is in the measure theory (face dimensions, the class of $F$, the states). The simplicial complex matters through its facets and its measure (E1–E4), not through its homology.

---

## 4. Charging the negative results to the dictionary: v2 candidates

Rule P2.1 says: before concluding anything about the theory, name what in the dictionary could have produced the result and change that. The three results that need charging are the non-Markov memory, the invisibility of training, and the failure of random restriction.

- **D1 (weight-only arrows).** Read $\Lambda_\ell$ from the restricted linear maps: $\sigma\to\tau$ exists iff the cone $\{h\ge0\text{ supported on }\sigma\}$ meets the pattern $\tau$ under $h\mapsto W_{\ell+1}h+b_{\ell+1}$, and carry the rank and the singular values of the block $W[\tau,\sigma]$ as the "dimension transported" (Note 1 §4.3, Murray–von Neumann comparison). This removes the dependence of the quiver on the input law (N2, N3) and separates what the weights allow from what the data uses.
- **D2 (restriction chosen by the theory).** Replace random sub-frames by the smallest set of units $S_\ell$ such that the restricted face $\sigma|_{S_\ell}$ makes the next face conditionally independent of the past: the empirical minimal sufficient sub-frame. Note 1 §6 says the minimal sufficient subalgebra is generated by the Connes cocycles; in the realized setting this is a greedy search over units with E2 as the objective, and whether $|S_\ell|$ shrinks under training is the first place training could become visible.
- **D3 (faces at a tolerance).** Read a face as an equivalence class of patterns sharing at least $k$ vertices (Connes–van Suijlekom's tolerance relation, Note 1 §7.2), which the $G$-continuity of §3.3 says is the natural topology on faces. This is the honest version of N1: coarsen by the structure the realization exhibits, not by width.
- **D4 (the state as the task).** Compare the quantities above under the training law and under a shifted or adversarial input law; a dictionary in which training is visible must at least see the difference between data the network was shaped by and data it was not.
- **D5 (multi-vertex output faces).** Use $\{c:\text{logit}_c>0\}$ or a top-$k$ set as the output face so that the last arrow has a non-trivial target complex (N4).

The order of leverage is D2, D1, D3: D2 is the only one that can make the memory vanish, D1 is the only one that reads the weights, D3 is the only one the data already asked for.

---

## 5. Messages to the other prongs

**To prong 1 (theory).** Three facts from the realization, offered as questions and not as definitions. (a) The realized law on histories is not Markov in the faces; the Gibbs/KMS family of Note 1 §5 is the Markov sub-family. Is there a canonical extension of the cocycle formalism to laws with memory, e.g. a cocycle on the path groupoid that is not of the form $F(\mu)-F(\nu)$ with $F$ additive, and what replaces the coboundary criterion there? (b) The arrows are $G$-continuous and transport $G$. Is there a theorem that the frame composition with kernel $G$ is the *only* composition rule compatible with a $G$-continuous transfer operator, which would make §3.3 a derivation rather than an observation? (c) The minimal sufficient sub-frame (D2) is what the Jenčová–Petz theorem constructs; what is its description for the layered quiver, and is it a face of the coherence complex?

**To prong 3 (unlocks).** Of the quantities prong 3 asked for: the quiver's square alternating sums of $F$ (U3) are computable from the output of E3 and will be added next; the links of the layer complexes are nearly complete at the widths tried (E6), so their $\gamma_j$ will be small for the trivial reason that the complexes are nearly simplices, and a meaningful link spectrum needs D3 first; the incidence patterns of consecutive layers do not repeat (face counts 236, 357, 370), so U5 needs an architecture with repeated structure before it can be tested.

---

## 6. Running it on real weights

```
python3 notes/experiments/mlp_bridge.py --weights net.npz --inputs X.npy
python3 notes/experiments/mlp_bridge.py --weights net.npz --inputs X.npy --restrict 16 --input-face pattern --json out.json
```

`net.npz` holds `W1,b1,...,WL,bL` with `W_l` of shape `(n_l, n_{l-1})` (transposed weights are detected), ReLU after every layer except the last; `X.npy` is an `(N, n_0)` array of inputs. If almost every input has its own face (E1 reports the fraction of samples whose face is seen at least twice), the resolution is too fine for statistics and `--restrict` is the v1 remedy; D2–D3 above are the v2 remedies. `--self-check` re-runs E2/E3 on the Markov surrogate so that estimator artefacts are not read as findings (programme rule P2.4).
