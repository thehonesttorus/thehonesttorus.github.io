# Measured algebras, detail modules and the return term

Working note XXVI. Assessment of a pasted synthesis (two parts, received as text; the attached `THEORY.md` and
`scattering_measured_algebras.zip` did not arrive in this session): that the more fundamental home of wavelets and
scattering is a measured algebra with its subalgebras and the operations between them, that algebraic grading
explains demodulation by the modulus, that a finite-index conditional expectation bounds scattering energy, that
noncommutative Dirichlet calculus is the source of scales and differences, and that the third-order transport needs
the action of the retained algebra on the detail, not only the detail's energy; plus a noncommutative-torus
version of scattering with a paired modulus. Each identity was checked by hand; the one that makes a measurable
prediction for the chain was implemented and run on the Azure VMs (section 3).

## 1. The claims, checked

| claim | status | note |
|---|---|---|
| module frame a = sum_j u_j E(u_j^* a) gives E(a^* a) = sum_j E(a^* u_j) E(u_j^* a) | correct | multiply by a^*, apply E, use E(x b) = E(x) b for b in B |
| M_2 over the diagonal with u_0 = I, u_1 = sigma_x: a = c_0 + u_1 c_1, c_1 = diag(a_21, a_12) | correct | the retained diagonal acts on the detail by the swap |
| grading: A_chi A_psi in A_{chi psi}, a in A_chi implies a^* a and \|a\| fixed | correct | exact demodulation of a homogeneous component |
| selection rule E(abc) = 0 unless chi psi eta = 1 | correct | averaging over the group is the projection onto the fixed algebra |
| E(h) >= lambda h for h >= 0 implies \|\|E h\|\|_2^2 >= lambda \|\|h\|\|_2^2, hence R_{m+1} <= (1 - lambda) R_m after the first modulus | correct | tau(E(h)^2) = tau(h E(h)) and tau(h^{1/2}(E(h) - lambda h) h^{1/2}) >= 0; the modulus preserves the tracial 2-norm |
| best constant 1/N for averaging over N points, 1/k for pinching into k blocks, none for a non-atomic space | correct | standard |
| A = T_{t_J}, W_j = (T_{2 t_j} - T_{2 t_{j+1}})^{1/2} is a Parseval bank | correct | the sum telescopes to T_0 = I |
| T_t(a^* b) - T_t(a)^* T_t(b) = 2 int_0^t T_s Gamma(T_{t-s} a, T_{t-s} b) ds | correct | differentiate F(s) = T_s[(T_{t-s} a)^* (T_{t-s} b)] |
| [K, f(a)] = sum_ij f^[1](lambda_i, lambda_j) P_i [K, a] P_j; for f = ReLU the divided differences lie in [0, 1], so the commutator contracts in Hilbert-Schmidt norm | correct | Daleckii-Krein; a Schur multiplier with entries in [0, 1] |
| P a^3 P = a_0^3 + a_0 K^* K + K^* K a_0 + K^* D K, and the ordered law P abc P = a_0 b_0 c_0 + Y_a X_b c_0 + a_0 Y_b X_c + Y_a D_b X_c | correct | insert P + Q twice; the example a_t = [[0, 1], [1, t]] has compressed moments 0, 1, t |
| the operator modulus is not 1-Lipschitz in Hilbert-Schmidt norm (constant sqrt 2, Araki-Yamagami); the paired modulus (\|a\|, \|a^*\|)/sqrt 2 is norm preserving and 1-Lipschitz in tracial L^2 | correct | the counterexample A = e_11, B = [[1/2, sqrt3/2], [0, 0]] gives 1 against 3/2; the proof via tau(\|x\|\|y\|) - tau(xy) = 2 tau(x_+ y_- + x_- y_+) >= 0 and the block H(a) is sound |

The mathematics is right throughout and the synthesis is careful about its own scope (it says where the index
bound gives nothing, that Morita equivalence does not preserve a modulus cascade, and that a cascade indexed by
words is not a representation).

## 2. Placement against the measurements

**The index bound gives nothing for the network, for a reason the synthesis states.** The coarse readout of the
neural problem is the mean over a continuous Gaussian input, a conditional expectation onto the scalars of a
non-atomic space, whose Pimsner-Popa constant is zero. Per neuron the same failure is visible: for x = ReLU(h),
E[x]^2 / E[x^2] is 1/pi at alpha = 0 and tends to zero as alpha -> -infinity (a third or more of the neurons sit
below alpha = -1 at depth), so no uniform lambda exists. This is the zero-gap statement of note XXV from the other
side: positive descendants are not forced to be visible to the mean.

**The exact gradings are the ones already used.** A selection rule needs a group acting by automorphisms under which
the carried objects are homogeneous. The network has two exact ones: the parity of the Gaussian input (odd chaos
has mean zero, built into the term programs' Hermite tables) and the dilation x -> c x (positive homogeneity: the
scale mixture and its radial zero mode, notes XIX to XXI). A random dense W is a symmetry only in law, so it gives
traces, not exact selection rules for the per-neuron quantities the score reads.

**The Dirichlet identities are the Gaussian calculus the chain already runs.** For the Ornstein-Uhlenbeck semigroup
Gamma is the carre du champ, the heat-defect identity is the Mehler expansion of the pair programs, and their order
is converged (closure-round). The divided-difference contraction is a stability statement; the chain's residual is
not a stability failure but an incoherent transport defect (closure-round section 5).

**The cubic identity explains the V32 join, and predicts its next step.** At a join the chain keeps the content in
span(Qn) and discards the rest (Q). In the exact dynamics the discarded content is transported by the next W, and
W does not preserve the orthogonal splitting, so part of it returns into the retained range: this is the term
Y D X of the ordered law, with X the retained-to-omitted coupling, D the transport inside the omitted space, Y the
way back. The shipped join drops the return entirely; projecting after W (note XXIV) keeps its first step, and
weighting the rows by Phi(alpha)^2 keeps the gate inside that step. The same law says the return continues through
the next transport W_{l+1} diag(Phi): a two-step reading metric O = c I + D_Phi W_{l+1}^T W_{l+1} D_Phi / 2 should
recover a further share at the price of two n x n x r products per join. That is a measurement (section 3).

**The chain already carries detail modules.** The proposal's target, a hierarchy of retained algebras with cheap
detail modules carrying the source and continuation actions, has instances in the chain: the S21 residual leg
(rank 4: the slice content the star core does not hold), the feedback thin legs (rank 16: the continuation of D21),
the nested tier inside the shared basis, and now the state-weighted return at the joins. The synthesis's three
discriminating questions have measured answers here. Does multiplication close cheaply? Not for the (i, i) slices:
the Hadamard product of legs has Khatri-Rao rank r^2 in any basis (note XV), the wall. Does transport preserve the
representation? Exactly, between joins (the shared basis is transported with the legs). Does the discarded
information matter to the readout? The representation's whole defect is now about 4.4% of raw at rank 320, and the
rest of the error (about 95%) is closure, a third of it the n^4 gate-covariance term.

**The nontracial state is the measured part.** The synthesis insists that a change of representation must keep the
weighting of the state that reads. In the chain that is literal: the join's projection is now taken in the state
diag(Phi(alpha)^2 + c) after W, which is a weighted, nontracial compression of the transported content, and it is
what made the new best (5.2925e-9).

## 3. Measurement: the second return step (`V32_JOIN_POST=4`, Azure VM, `outputs/return2_16nets.txt`)

The two-step reading metric O = c I + D_Phi W_{l+1}^T W_{l+1} D_Phi / 2, with D_Phi the gate of the join layer and
W_{l+1} the next layer's weights (He normalisation, so the second term has the mean of the one-step version), and
the O-orthogonal projection onto span(W Qn) through a Cholesky factor of Y^T O Y. Paired against the current best
(one-step weighted metric, rank 320) on official networks 0-15:

| metric | raw | C/B | adjusted against current best |
|---|---|---|---|
| two-step, c = 0.1 | -1.03% +- 0.24 (14/16 better) | 0.2349 (+2.3%) | +1.25% |
| two-step, c = 0.3 | -0.59% +- 0.23 (12/16) | 0.2349 | +1.70% |

The prediction holds: the omitted content keeps returning through the second transport, and keeping that return
lowers the raw error on 14 of 16 networks. It does not pay: the two extra n x n x r products per join cost 2.3% of
the bill for 1% of raw, and 1% of raw is a quarter of what 32 ranks are worth here (+4.5% raw from rank 320 to 288),
so it cannot buy a rank reduction either. The return series converges fast enough that its first step (the
gateway W and the gate) is where the value is; that step costs nothing beyond the projection and is in the new best.

Two other weightings tried on the way, the local reads weighted by their Edgeworth coefficient (alpha phi(alpha))^2
or by phi(alpha)^2 instead of a constant floor, are ill-posed: both weights vanish together with Phi(alpha)^2 on the
saturated-off neurons, the weighted basis divides by the square root of the weight, and the runs overflow (NaN on
some networks, killed on others). With a floor added they reduce to the measured constant version; on network 0 the
(alpha phi)^2 form was 1.5% worse in raw before any floor.

## 3b. The ladder under the new metric (`outputs/age_gates_16nets.txt`)

With the projection now taken in the reading state, the question was whether the age gates should move (confine
earlier, nest earlier or later) and whether the thin legs' ranks are still right. Paired against the current best
on networks 0-15 (an earlier submission of this sweep dropped the first flag of each variant in the job builder;
these are the corrected runs):

| variant | raw | C/B | adjusted |
|---|---|---|---|
| confine after age 3 (rank 320 / 352 / 384 / 416) | +16.9% / +5.3% / +0.6% / -2.5% | 0.2149 / 0.2251 / 0.2356 / 0.2466 | +9.4% / +3.2% / +3.2% / +4.7% |
| confine after age 2 (rank 384 / 448) | +33.2% / +5.5% | 0.2228 / 0.2464 | +29.3% / +13.2% |
| nested tier after age 6 / 8 (rank 192) | +4.75% / -0.79% | 0.2281 / 0.2308 | +4.07% / -0.27% +- 0.41 |
| confine after 3 and nest after 6 | +19.8% | 0.2134 | +11.4% |
| feedback thin legs rank 16 -> 8 | +1.70% +- 0.21 | 0.2220 | **-1.67% +- 0.21** |

Confining earlier does not pay at any rank: the age-3 and age-2 content needs a rank whose formation and joins
cost about what the dense legs it replaces cost, which is the diffusion-wavelet ladder of note XXV read
quantitatively (the scaling space after three transports is still wider than the knee). The age gates of v29
(4 and 7) stay. The feedback legs carry less than their rank 16: rank 8 loses 1.7% of raw and saves 3.3% of the bill,
consistently on every network. Lower is better down to rank 1-2 (16 networks: rank 12 / 8 / 6 / 4 / 3 / 2 / 1 give
-1.03 / -1.67 / -1.92 / -2.27 / -2.34 / -2.59 / -2.60% adjusted; removing the legs is +3.1%), and on all 100 networks
(`outputs/feedback_rank.txt`) rank 4 gives -2.75% +- 0.10 (mean adjusted 5.1478e-9) and rank 2 with the nested gate at
8 gives **-3.45% +- 0.17, mean adjusted 5.1091e-9 (raw 2.3468e-8 at C/B 0.2177)**, the current best: add
`V18_R_FB=2 V24_AGE_OLD2=8` to the configuration of note XXIV. The D21 feedback carries two or three directions of
content, not sixteen. This closes the parameter work; what follows is the representation question.

## 3c. A second pasted synthesis: Markov geometry, index and source-query bounds

A further synthesis (pasted; its companion PDF was a link to another sandbox and did not arrive) builds the calculus
from a single conditional expectation: the Stinespring identity Phi(a^* b) - Phi(a)^* Phi(b) = d(a)^* d(b) with
d = (I - VV^*) pi(.) V, the composition law Gamma_{Psi Phi} = Psi(Gamma_Phi) + Gamma_Psi(Phi, Phi) (total
covariance), the relative bimodule A (x)_B A with delta_E a = a (x) 1 - 1 (x) a and partial^* partial = I - E, the
hierarchy L = sum c_j (I - E_j) with the martingale details as eigencomponents, the circle example
\|\|[F, M_f]\|\|_HS^2 = 4 sum_k \|k\| \|f^(k)\|^2 with Index T_u = A_- - A_+ and (1/4)\|\|[F, M_u]\|\|^2 = A_- + A_+ for
unitary symbols, the Pimsner-Popa uniform bound \|\|a - Ea\|\|^2 <= Ind_PP(E) \|\|Gamma_E(a, a)\|\|, Fisher-information
loss \|\|h - Eh\|\|^2, Petz sufficiency through Connes cocycles, and, for the network, the signed source-query error
epsilon(c, s) = sum_{j not kept} <Q_j c, Q_j s> with backward readouts c_l = T_l^* c_{l+1} and the exact
accumulation <c_N, e_N> = <c_0, e_0> + sum_l <c_{l+1}, r_l>. I checked the identities (the Fourier count of the
Hilbert-Schmidt norm, the trace identity for the index, the bimodule energy, the projection counterexample
(5/6, 1/3, -1/6)); they are correct.

Its operational conclusion, that the source and the query must be jointly compressible and the query is the
backward readout through the actual transport, is the principle the V32 join implements and measures: the
one-step backward readout (W, then the gate) is the reading metric of the new best, and the two-step readout
(section 3) is real but costs more than it returns. The intrinsic geometry it proposes from the transport,
D_T = [[0, T^*], [T, 0]] with energy 2 sum \|T_ij\|^2 \|f_j - g_i\|^2, is the post-W metric written as a Dirac
operator. So the "missing implication" it names (cheap geometry compatible with the transport, plus regular sources
and queries, gives cheap signed transport) has a measured instance with its limit: the first backward step pays
7.86%, the second does not, and the closure-side residual is not a compression error at all.

## 4. Decision

| from the synthesis | what it is here | decision |
|---|---|---|
| algebra-valued frames, detail modules | the chain's residual and feedback legs, nested tier, and the weighted join | already present; the frame language is accurate |
| grading and selection rules | parity and dilation, both used | no new exact rule for random W |
| finite-index scattering bound | lambda = 0 for the Gaussian mean; no uniform per-neuron constant | consistent with zero gap; no lever |
| Dirichlet calculus, heat defect, divided differences | OU carre du champ = the converged Mehler programs; a Lipschitz statement | relabelling |
| cubic module identity, return term Y D X | explains why the post-W, gate-weighted join works | second step measured: real (-1.03% raw, 14/16) but costs more than it returns (+1.25% adjusted) |
| nontracial (KMS-symmetric) state | the weighted join metric | measured: part of the new best |
| paired modulus on the noncommutative torus | a correct noncommutative scattering construction | theory only; no bearing on the estimator |
