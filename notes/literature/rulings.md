# Per-paper rulings

All 49 papers, after read and adversarial challenge. Cluster keys: metric (quantum
metric spaces), trunc (spectral truncation), fuzzy (fuzzy/matrix geometry), connesd (Connes distance),
sdp (moment hierarchies), closure (dynamical closure/Koopman), state (state approximation), found
(foundational).

| paper | cluster | verdict | reason |
|---|---|---|---|
| Frederic Latremoliere, "How to approximate the flat spectral triple of a quantum torus by fuzzy tori: a twiste | metric | not-relevant | I retrieved the full text and checked each claim; the original verdict stands, and the best adversarial card is defeated by our own measurements rather than by the paper's irrelevance alone.

(1) The Fejer card dies on obstruction 1's exact form. Our failure signature is not merely "compression fail |
| F. D'Andrea, F. Lizzi, P. Martinetti, "Spectral geometry with a cut-off: topological and metric aspects" (2013 | trunc | not-relevant | The original verdict survives, and the strongest adversarial reading reinforces it. This is a 2013 math-ph paper on how Connes' spectral distance behaves under a momentum cut-off: Gromov-Hausdorff convergence of truncated state spaces (asymptotic, no rate), a minimal-length result under a bounded or |
| Konrad Aguilar, Jens Kaad, David Kyed, "Polynomial approximation of quantum Lipschitz functions", Documenta Ma | metric | not-relevant | Full text read and independently verified; the other analyst's verdict survives, though one of their stated reasons needed correcting in a way that strengthens the dismissal. The paper resolves a domain technicality — polynomial-domain L_{D_q} versus maximal-domain L^max_{D_q} on the Podles sphere a |
| David Kerr, Hanfeng Li, "On Gromov-Hausdorff convergence for operator metric spaces" (2004) | metric | not-relevant | Verdict confirmed on full text, with the reasoning tightened rather than overturned. The paper contains: an amalgamated-sum construction existing only to give a triangle inequality; a theorem identifying two definitions of a matricial quantum Gromov-Hausdorff distance (dist_op = dist_s, Thm 3.7); Li |
| Francesco Flora, Losel Matos, Tamas Krivachy, Antonio Acin, "Moment Optimization in the Navascues-Pironio-Acin | sdp | not-relevant | The challenge fails on four points, the first two decisive.

(a) WRONG PRIMITIVE — the paper forbids the only move that has ever worked here. Confirmed in the text: the decision variable is strictly binary with fixed Hamming weight (Sec. 3.1 lines 626, 641, 655; the entire Gumbel top-k apparatus, li |
| Yendrembam Chaoba Devi, Alpesh Patil, A. Bose, F. G. Scholtz et al, "Revisiting Connes' finite spectral distan | connesd | not-relevant | Both prongs fail on inspection of the actual text.

PRONG A fails on the hypothesis, not the conclusion. The collapse theorem's premise is a transitive unitary action on trace-norm level sets of Delta-rho. We have no such action. The input distribution N(0,I) is O(n)-invariant, but the WEIGHTS are f |
| A. Abanov, Luca Candelori, H. Steinacker, Kharen Musaelian et al, "Quantum Geometry of Data" (2025) | fuzzy | not-relevant | The challenge is the strongest available and still fails on three counts. (i) It is retrospective, qualitative and strictly weaker than what we have already measured. Both passages are one-paragraph asides with no equations, no experiment, and an explicit "out of the scope of this work"; the paper h |
| Lukas Schneiderbauer, H. Steinacker, "Measuring finite quantum geometries via quasi-coherent states" (2016) | fuzzy | not-relevant | Full text read and the other analyst's claims independently verified — their assessment is accurate in every particular I checked, including the verbatim hierarchy quote (lines 4369-4370), eq. (6.9), the squashed CP^2_N figures, the 8.6e-14 torus zero mode, and the N=231/d=6 scaling claim. The adver |
| A. Connes, Walter D. van Suijlekom, "Spectral Truncations in Noncommutative Geometry and Operator Systems", Co | trunc | not-relevant | Each of the four legs breaks, and they break for reasons already measured rather than for reasons of taste.

(1) is true and is precisely why the paper cannot help: it supplies the REPRESENTATION, and obstruction 4 measures the representation defect at the readout to be ~500x smaller than the total  |
| A. P. Balachandran, G. Bimonte, Elisa Ercolessi, P. Teotonio-Sobrinho et al, "Finite quantum physics and nonco | fuzzy | not-relevant | The challenge lands one real correction and no usable content, so the verdict stands.

The correction is real and should be recorded: Section 6.3 does NOT belong to the failed corner-restriction family. It is a group-averaging conditional expectation onto an invariant subalgebra, which sits on the a |
| Frederic Latremoliere, "Approximation of quantum tori by finite quantum tori for the quantum Gromov-Hausdorff  | metric | not-relevant | Every strand of the challenge breaks on contact with the text or with our own diagnostics.

(1) Lemma 3.7's template is vacuous for us in the precise technical sense that matters. A bound of the form //a - Pa// <= delta*L(a) is UNIFORM over the Lipschitz ball. A uniform-over-the-ball bound can only  |
| G. Fiore, F. Pisacane, "Energy cutoff, effective theories, noncommutativity, fuzzyness: the case of O(D)-covar | trunc | not-relevant | The challenge is the best available reading and it still fails, on four independently sufficient grounds, all checked against the retrieved text rather than the prior summary. (1) THE COUNTERTERM IS NOT A FITTING PROCEDURE. Fiore-Pisacane's f_1, f_2 are fixed by exact so(D+1) representation theory o |
| David Kerr, "Matricial quantum Gromov-Hausdorff distance" (2002) | metric | not-relevant | The adversarial pass strengthens the original verdict rather than weakening it. Every claim in the assessment checks out against the full text. The one route that looked live — matrix states seeing what scalar states cannot, rhyming with "marginals do not determine the joint" — dies on a fact the an |
| A. Licht, "Coulomb Potential of a Point Mass in Theta Noncommutative Geometry" (2006) | connesd | not-relevant | All three hooks collapse under our own measurements, and two of them collapse into things we have already run.

HOOK 1 IS SELF-REFUTING AGAINST OBSTRUCTION 1. A latent-shift Gaussian representation of the source non-Gaussianity is, concretely, a scalar latent zeta times a direction v, giving kappa_3 |
| H. Steinacker, "On the quantum structure of space-time, gravity, and higher spin in matrix models", Classical  | fuzzy | not-relevant | All three adversary cases fail on inspection, and I could not construct a fourth.

(1) dies on decision value and on object mismatch. The "white noise" line is one sentence of prose motivating why the IKKT ACTION concentrates on almost-commutative configurations; it is not a theorem, carries no quan |
| Samuel Buckley-Bonanno, Noah I. Eckstein, Susanne F. Yelin, "Quantum simulation of gauge theories on dynamical | fuzzy | not-relevant | Every steelman closes against the text, and three of them close against our own measured numbers.

(A) fails on structure of the cancellation. BCH/Magnus cancellation is exact because e^A e^B e^{-A} e^{-B} is a Lie-group identity: the terms come in *matched pairs* that annihilate symbolically. Our c |
| D. Karabali, S. Randjbar-Daemi, V. P. Nair, "Fuzzy spaces, the M(atrix) model and the quantum Hall effect" (20 | fuzzy | not-relevant | Each hook dies on contact with our measurements, for independent reasons. (1) The soft-penalty prescription runs the wrong way on compute. It is a well-definedness device for a field theory, not a reduction technique: it ADDS an action term and keeps a full-rank operator, so it saves zero FLOPs. Und |
| M. Rieffel, "Metrics on states from actions of compact groups", Documenta Mathematica (1998) | metric | not-relevant | The verdict survives, though two of the analyst's supporting arguments do not and should not be reused.

Attack (1) corrects their reasoning without moving the answer. Ergodicity is available on M_N via clock-shift, so the theorems can be instantiated -- and instantiating them shows they say nothing |
| Nicola d'Alessandro, C. R. I. Carceller, Armin Tavakoli, "Semidefinite block-matrix relaxations for computing  | sdp | not-relevant | Full text retrieved and read independently; the assessment's summary is accurate in every checkable particular (BMM definition and PSD-by-construction proof; NPA recovered at Theta(X)=tr(X psi); type-(i)/type-(ii) dichotomy; localising BMMs adding constraints without variables; Lambda_i = Omega_{Lam |
| Wei Wu, "Non-Commutative Metrics on Matrix State Spaces" (2004) | metric | not-relevant | Each of the four arguments above collapses on contact with the actual text, and I checked the text myself rather than relying on the summary.

(2) is the strongest and it is an equivocation, refutable by a theorem. Our algebra is COMMUTATIVE — functions of the Gaussian input x. For a commutative dom |
| Shivraj Prajapat, Yendrembam Chaoba Devi, A. Mukhopadhyay, F. G. Scholtz et al, "Connes distance function on f | connesd | not-relevant | The challenge fails on two grounds, the second of which is new and stronger than the first analyst's argument. (1) PROVENANCE OF THE HOOK: the sup-over-a-dual-ball template is Connes (1989) and Rieffel, not this paper. This paper fixes D a priori from the monopole-twisted geometry on S^3/S^2, inheri |
| Zhen-Peng Xu, Rene Schwonnek, A. Winter, "Bounding the Joint Numerical Range of Pauli Strings by Graph Paramet | sdp | not-relevant | All three legs of my own challenge break, and the first two break in ways that make the paper LESS useful than the analyst thought, not more.

LEG 1 FAILS BECAUSE THM. 26 IS CANCELLATION-BLIND BY CONSTRUCTION. I checked the proof: it is Cauchy-Schwarz, /sum_i a_i <A_i>/^2 <= (sum_i a_i^2/w_i)(sum_i  |
| M. Rieffel, "Metrics on state spaces", Documenta Mathematica (1999) | metric | not-relevant | Every strand of the challenge breaks, and two break on evidence the original analyst did not cite.

(1) fails on a detail I checked in the source and which strengthens their verdict: our propagated object is NOT a state. Obstruction 3 records exactly one negative eigenvalue at about -1e-4 of the tra |
| A. Connes, "Compact metric spaces, Fredholm modules, and hyperfiniteness", Ergodic Theory and Dynamical System | found | not-relevant | Every leg breaks, and leg (2) — the one that would have mattered — breaks on an experiment we have already run.

Leg (2) is refuted directly by our own notes, not by argument. notes/groupoid-carrier-and-cartan-truncation.md lines 174-180 records that the Cartan-compatible corner 1_Y A 1_Y "does slig |
| D. Freeman, D. Giannakis, Brian Mintz, J. Slawinska, "Data assimilation in operator algebras", PNAS (2022) | closure | not-relevant | Verdict survives on independently verified full text, and one of the analyst's concessions should be withdrawn as too generous.

What the paper is: sequential Bayesian filtering re-derived inside a non-abelian operator algebra, so that finite-dimensional truncation Pi_L A Pi_L is completely positive |
| A. P. Balachandran, Francisco Calderon, V. P. Nair, S. Vaidya et al, "Uncertainties in quantum measurements: a | state | not-relevant | Each limb of the challenge dies on contact with the actual text or with our own numbers.

(1) The lift-ambiguity limb is a redescription of something we already have in strictly stronger form. finite-resolution-question-algebra.md line 69-73 already derives the positive-extension set for our restric |
| M. Rieffel, "Gromov-Hausdorff Distance for Quantum Metric Spaces" (2000) | metric | not-relevant | The Thm 8.2 challenge is the best case available and it fails on three independent grounds, each checked against the paper's text rather than the summary.

(i) THE HYPOTHESIS DOING ALL THE WORK IS THE ONE WE LACK. delta_n -> 0 only because the truncation index and the Lip-norm are the SAME grading:  |
| D. Giannakis, Michael Montgomery, "Koopman and transfer operator techniques from the perspective of quantum th | closure | not-relevant | Retrieved and read in FULL TEXT (4,835 lines) via alphaXiv, not abstract. I independently confirm the analyst's factual claims: it is a survey chapter of the authors' own prior work [17,18,25,28,30,31,33] with no new theorems (Lemma 6 is "= [30, Lemma 8]", Prop. 9 is "= [31, Theorem 6]") and no nume |
| Simon Becker, Wuchen Li, "Quantum Statistical Learning via Quantum Wasserstein Natural Gradient", Journal of S | state | not-relevant | Every load-bearing piece of the challenge above turns out to be either not this paper's contribution, or to point at the wrong metric, or to be unaffordable. Taking them in order.

ATTRIBUTION. None of the transferable content is Becker-Li's. Metric-aware/natural-gradient preconditioning is Amari (1 |
| C. Rovelli, Simone Speziale, "A semiclassical tetrahedron" (2006) | fuzzy | not-relevant | Full text read from the arXiv source; the original assessment is accurate and I found nothing it missed. The paper constructs one state in a fixed SU(2) intertwiner space of dimension O(min j) so that six pre-chosen geometric operators have classical expectation values in the uniform large-spin limi |
| P. Bertozzini, R. Conti, "Non-Commutative Geometry, Categories and Quantum Physics" (2008) | found | not-relevant | All three attacks die on inspection, and the first dies on our own measured data.

(1) The rank-one reading collapses three ways. (a) The rank-one Fell bundle is rank one *on every arrow (A,B) simultaneously* — a separate line fibre over each of the /O/^2 pairs. That is n^2 independent complex param |
| Mateus Araujo, Andrew J. P. Garner, Miguel Navascues, "Non-commutative optimization problems with differential | sdp | not-relevant | Every strand of the challenge dies on our own numbers, and the bracket strand dies twice.

(1) BRACKET WIDTH, the decisive kill. Target: MSE ~1e-9 means per-neuron absolute error ~3e-5 against a scored-layer neuron variance of 7.8e-2 (sigma ~0.28), i.e. relative accuracy ~1e-4. Their brackets at lev |
| L. Barbieri-Viale, "The Infinitesimal Structure of Quantum Information" (2026) | state | not-relevant | Every strand of the challenge is closed by the paper's own text, and I checked each one against the source rather than the summary. (1) The convolution hook dies at line 2039, where the paper pre-emptively disavows exactly the reading I was constructing: "In classical applications of dual numbers (s |
| D. Anshu, Therese Jekel, Basa Landry, "Quantum Wasserstein distances for quantum permutation groups", Journal  | metric | not-relevant | Every hook above dissolves on contact with the actual text, and I verified each by grep rather than by trusting either summary. (1) Lemma 4.15 is the best hook and it still fails, for three independent reasons. First, "the cost of a product relaxation equals the variance of the mixture" is not a tec |
| J. Schwarz, Bastian Boll, Daniel Gonzalez-Alvarado, Christoph Schnorr et al, "Quantum State Assignment Flows", | state | not-relevant | All four angles die, and three of them die on measurements we already have rather than on judgment.

(1) The metric angle is killed by diagnostic block 1 of results.txt, which the previous analyst did not invoke and which is the decisive datum. We already ran the unit-space truncation with site_pow  |
| Markus Faulhuber, Thomas Strohmer, "Quantum paving: When sphere packings meet Gabor frames" (2024) | found | not-relevant | Each of the three arguments dies on inspection of the actual text, and two of them die on our own measurements.

(1) The Janssen/Poisson step needs an exact discrete abelian group in the summation index to have a Pontryagin dual - that is the entire content of Lemma 3.5 and eq. (12) (the adjoint lat |
| Hua-qing Zhou, Ting Gao, Fengli Yan, "Optimal convex approximation of qubit states and geometry of completely  | state | not-relevant | Each line of attack collapses on inspection of the actual text. (1) The shared object is the shared object at d=2: the convex set is the hull of at most six points in a three-parameter Bloch ball, and the whole technical apparatus is a case enumeration over which of six p_i vanish (Appendix A, cases |
| R. P. Kostecki, "Quantum theory as inductive inference" (2010) | state | not-relevant | The original verdict survives, though the strongest argument against it is better than the prior analyst made it, and fails for a sharper reason than any they gave.

RETRIEVAL AND VERIFICATION. I obtained the full unabridged v4 text and independently re-ran their keyword census: moment 0, cumulant 0 |
| M. Karasev, "Adiabatics using quantum action" (2014) | closure | not-relevant | The verdict survives, though one of its supporting arguments does not and should be retired.

Point (e) is correct and the original reason should not lean on it: cancellation is not an argument against averaging. The real argument is that the &-projection is not a DEFINABLE OPERATION on our object,  |
| J. Barrett, "Matrix geometries and fuzzy spaces as finite spectral triples" (2015) | fuzzy | not-relevant | Every one of the four pushes dies on contact with the measured facts, and two of them die in a way that actually STRENGTHENS the original verdict.

(1) The Frobenius projector kills itself when instantiated. piK requires a *-algebra A acting on H; we have no algebra action on our transported objects |
| Li-qiang Zhang, Deng-hui Yu, Chang-shui Yu, "The best approximation of a given qubit state with the limited pu | state | not-relevant | The challenge is the best case available and it still fails, for reasons that are structural rather than rhetorical.

KILLER 1 — Caratheodory is post-hoc, and the bottleneck is upstream of it. The constructive algorithm is: hold d+2 vectors, find a linear dependence, eliminate one, repeat. It takes  |
| H. Steinacker, "Emergent geometry and gravity from matrix models: an introduction" (2010) | fuzzy | not-relevant | All three cards break, and (A) breaks in the most instructive way.

(A) fails on two counts. First, the aggregate it preserves is the wrong aggregate. Eq. (21) preserves ONE scalar functional, Tr. Our aggregate is the 1024-vector obtained by contracting the source path sum against the propagated mea |
| M. Vancliff, "The Interplay of Algebra and Geometry in the Setting of Regular Algebras" (2015) | found | not-relevant | The determinantal-codimension argument is the best case available and it fails on a decisive test: it has no discriminating power between our failure and our success. The identical codimension count applies to the omitted sector O, which results.txt line 237-239 records as having rank >= 8 while its |
| Bhishan Jacelon, "Metrics on trace spaces", Journal of Functional Analysis (2021) | metric | not-relevant | Full text retrieved and checked against the analyst's claims; every checkable particular (Thm A statement, Lipschitz family Def 4.3/4.6, transport constant Def 2.11, Hall's marriage eigenvalue matching, Pomeau-Manneville CLT in Thm 5.1, Rieffel [38] in the bibliography) is accurate. All three advers |
| D. Oriti, Daniele Pranzetti, J. Ryan, Lorenzo Sindoni, "Generalized quantum gravity condensates for homogeneou | fuzzy | not-relevant | The challenge fails on all three points, and the failure is verifiable rather than a matter of taste.

Point (2) collapses first. The "negative prediction" is a motivational argument in a physics introduction with no metric, no norm, no bound, and no counterexample — the word "probably" is in the se |
| D. Voiculescu, "Hybrid Normed Ideal Perturbations of n-Tuples of Operators II: Weak Wave Operators" (2017) | found | not-relevant | The Kuroda hook is the only real one and it closes on a precise mathematical distinction: a DIRECT SUM is not a SUM. /X (+) X/_Phi < 2/X/_Phi holds because a symmetric norming function is sublinear on the singular-value sequence of orthogonally supported blocks — nothing cancels, the singular values |
| A. Connes, "Trace de Dixmier, modules de Fredholm et geometrie riemannienne" (1988) | found | not-relevant | The verdict survives, but their reasoning does not; the correct kill is different and stronger.

The steelman (Connes' trace theorem as an alternative representation of E_gamma[h_16]) dies four times over. (a) Evaluating mu_k(M_{h_16} /D_gamma/^{-d}) requires the matrix of M_{h_16} in the Hermite ba |
| Quantum Mechanics for Closure of Dynamical Systems (Dartmouth, 2022), arXiv:2208.03390 | closure | not-relevant | Their verdict survives, and I confirm their reading of the paper: full text checked, no numbered theorems, positivity-via-positive-compression (secs. 4.3.1, 9.1) is the entire structure-preservation contribution, the only convergence statement is the iterated limit eq. 6.10 with no finite-L rate, co |
| Second quantization for classical nonlinear dynamics (2025), arXiv:2501.07419 | closure | not-relevant | All three prongs die on our own measurements, and prong 1 dies on the sharpest one.

PRONG 1 fails on obstruction 4, and fatally. The paper's entire m-grading apparatus (Thm. 6, sec. 6, the 4.6M moment coefficients) is a REPRESENTATION theorem for the observable at readout. Obstruction 4 measures ou |
