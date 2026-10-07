# Wavelets, scattering, and the structure underneath them

Working note XXV. Asked: go deep into wavelets and scattering transforms and their development in noncommutative
geometry, starting from Lempereur and Mallat, *Hierarchic flows to estimate and sample high-dimensional
probabilities* (arXiv 2405.03468), and Bodmann and Emilsdottir, *A scattering transform for graphs based on heat
semigroups* (arXiv 2208.12773); then, mid-way, whether this is a special case of something more fundamental in NCG
that links the theories of the earlier notes. Read in full or in the relevant sections: those two; Arias, Barbieri
and Hernandez, *Scattering networks on noncommutative finite groups* (2505.20950); Marcolli and Paolucci,
*Cuntz-Krieger algebras and wavelets on fractals* (0908.0596); Farsi, Gillaspy, Julien, Kang and Packer, *Wavelets
and spectral triples for fractal representations of Cuntz algebras* (1603.06979); Gerontogiannis, Goffeng and
Mesland, *Heat operators and isometry groups of Cuntz-Krieger algebras* (2406.07416); Gerontogiannis and Goffeng,
*Five shades of KMS* (2605.31390); Kribs, *Quantum channels, wavelets, dilations and representations of O_n*
(math/0309390); Kim, Kribs, Lozano, Pereira and Plosker, *Quasiorthogonality of commutative algebras, complex
Hadamard matrices, and mutually unbiased measurements* (2504.18741); Bratteli and Jorgensen on Cuntz algebras and
multiresolution (funct-an/9612003); Jaffe, Jiang, Liu, Ren and Wu, *Quantum Fourier analysis* (2002.03477). The
measurements are in note XXIV (the V32 family) and are summarised in section 5.

## 1. What the sources establish

**Scattering (Mallat; Bodmann-Emilsdottir; Arias-Barbieri-Hernandez).** A scattering transform alternates a
Parseval filter bank with a pointwise modulus, U[p]f = |psi_{j_m} * ... |psi_{j_1} * f||, and reads the result with
a low-pass, S[p]f = phi * U[p]f. Three theorems carry over to every version read here: a Pythagoras identity per
layer (the propagated energy splits into what goes on and what is read out), non-expansiveness and Lipschitz
stability (||Sf - Sf'|| <= ||f - f'||), and energy decay under an admissibility condition. Bodmann and Emilsdottir
get the decay from positivity: the heat semigroup T_t = e^{-t Delta/2} is positivity preserving, so
||T|g||| >= ||T g||, so with S_t = (I - e^{-t Delta})^{1/2} the high-pass satisfies ||S|g||| <= ||S g|| (a
Beurling-Deny inequality that removes the nonlinearity), and ||g_k|| <= (1 - e^{-t lambda_max})^{k/2} ||f||. On a
finite group the filters are class functions with matrix Fourier coefficients gamma(pi) I, Parseval iff
sum_j |gamma_j(pi)|^2 = 1, and the decay rate is 1 - min_pi |gamma_0(pi)|^2; for nonnegative signals (what a
modulus produces) a theorem of Kueh, Olson, Rockmore and Tan forces a fixed fraction of the energy onto the
representations that contain the trivial one, so the modulus itself pushes energy to low frequency.

**Hierarchic flows (Lempereur-Mallat).** The law of a field factorises across scales,
p(phi) = w^{-1} p_J(phi_J) prod_j pbar_j(phibar_j | phi_j), with phi_j = G_j phi_{j-1} the low-pass and
phibar_j = Gbar_j phi_{j-1} the wavelet coefficients, renormalised to unit variance coordinate by coordinate
(D_j = diag(sigma^{-1})). Each conditional is an exponential family in a scattering covariance: the coefficients
phi * psi~ and their moduli filtered by a second wavelet transform, |phi * psi~_{j'}|^q * psi~_l, whose
correlations Lambda_2 (power spectrum), Lambda_3 (modulus against coefficient, the phase alignment that carries
skewness) and Lambda_4 (modulus envelope against modulus envelope) are kept only for the same wavelet at the same
position, O(log^3 d) parameters in all. The reason given is structural: wavelet coefficients are nearly
uncorrelated across scales, their moduli are strongly correlated, and pairs built from different wavelets have
nearly disjoint frequency supports, so under stationarity and locality they are negligible. What is proved is the
algebraic form of the hierarchic energies (Theorems 4.1, 5.1) and lower bounds on log-Sobolev constants; that the
conditionals have bounded log-Sobolev constants is stated as a conjecture, supported by relaxation-time
measurements on phi^4.

**Wavelets as noncommutative geometry (Bratteli-Jorgensen; Marcolli-Paolucci; Farsi et al.; Kribs; Gerontogiannis,
Goffeng, Mesland).** A multiresolution analysis of scale N is a representation of the Cuntz algebra O_N: the filters
give isometries S_i f(z) = m_i(z) f(z^N) with S_i^* S_j = delta_ij and sum_i S_i S_i^* = 1, and the latter is the
perfect-reconstruction condition. Kribs shows these representations are classified by a single unital completely
positive map Phi(X) = sum_i A_i X A_i^* on a finite-dimensional space (the compressions of the isometries), and
that every such map dilates to a representation of O_n; wavelets are the case where the Kraus operators are
polyphase filters. On the path space of a Cuntz or Cuntz-Krieger algebra (the Cantor set of a Markov partition),
Farsi et al. prove that the wavelet spaces are exactly the eigenspaces of the Laplace-Beltrami operator of the
Pearson-Bellissard spectral triple, W_n = direct sum over words |gamma| = n of E_gamma, with the Cuntz isometries
carrying E_gamma unitarily onto E_{i gamma}; with uneven weights the wavelets are orthonormal for the inner product
weighted by the Markov measure, not the flat one. Gerontogiannis, Goffeng and Mesland build the log-Laplacian on
the Deaconu-Renault groupoid of the Markov chain, diagonalised by Haar wavelets with eigenvalues that grow linearly
in the word length (logarithmically in scale), heat operator a Riesz potential up to rank one, spectral triples
exhausting K^1(O_A). Gerontogiannis and Goffeng show the eigenvalues are statistically linear in length (for
tau-almost every path, lambda(pi_n x)/n -> lambda_A tau(F_A)) and recover the KMS state of the gauge action, the
Gibbs measure of the Markov chain, from local heat traces at the critical time t_c where the pressure vanishes.

## 2. The structure underneath: Cartan pairs in relative position

Read together, the sources say that a wavelet or scattering transform is three pieces of operator-algebra data,
each of which exists far more generally:

1. **A commutative subalgebra where the nonlinearity acts.** The modulus is pointwise: functional calculus in the
   diagonal (position) algebra Delta. In the Cuntz-Krieger picture this is C(Sigma_A), the Cartan subalgebra of the
   groupoid algebra, with the conditional expectation onto it given by restriction to the unit space.
2. **A second commutative subalgebra where the linear step is diagonal, in a definite relative position to the
   first.** Convolution filters are diagonal in the Fourier algebra F Delta F^*. The relative position of two
   maximal abelian subalgebras is measured by Weiner's quasi-orthogonality Q(A, B) = Tr(T_A T_B), which is 1
   exactly when E_A E_B = E_{C1} (a commuting square), n when they coincide, and for Delta against U Delta U^* equals
   the squared Frobenius norm of the doubly stochastic matrix U o conj(U) (Kim, Kribs et al., Proposition 2.2 and
   Corollary 3.8): quasi-orthogonal if and only if sqrt(n) U is a complex Hadamard matrix. The Fourier matrix is the
   canonical Hadamard matrix, so position and frequency are an exact commuting square, and Jones' spin-model
   construction turns any complex Hadamard matrix into a subfactor whose planar-algebra Fourier transform (the
   quantum Fourier analysis of Jaffe, Liu and coauthors) is the ordinary one for the Fourier matrix.
3. **A perfect-reconstruction decomposition between scales, orthogonal for a specified state.** Cuntz relations, or
   more generally a Pimsner-Popa basis of an inclusion with a conditional expectation; multiresolution is a tower
   of such inclusions (a filtration of conditional expectations: for Haar wavelets, literally the martingale
   differences of the dyadic filtration), and the orthogonality is in the GNS inner product of the state that
   defines the expectation (the Markov/KMS measure in Farsi et al., the trace in Kim-Kribs, where a conditional
   expectation is an orthogonal projection only for that inner product).

So the more profound object is a **Cartan pair together with a second Cartan subalgebra and a tower of state-preserving
conditional expectations**: Renault's correspondence between Cartan pairs and twisted etale groupoids, Popa's
commuting squares and quasi-orthogonality, and the Jones basic construction. Classical wavelet scattering is the
case where the second Cartan is the Fourier one (a group, a Hadamard matrix, an exact commuting square), which is
what makes the cross-wavelet terms vanish and convolution cost n log n.

**The ReLU network is the same three pieces with the structure removed.** The gate acts in the neuron diagonal
Delta_l (item 1, exactly). The weights move the diagonal to W_l Delta W_l^T, a Gaussian random relative position
(item 2). For a Haar-random unitary E sum_ij |U_ij|^4 = 2n/(n + 1) (3n/(n + 2) for a real orthogonal one), so
Q(Delta, U Delta U^*) is about 2 to 3: one or two units of defect out of a possible n - 1, small relative to the
dimension and of order one in absolute terms, re-drawn at every layer (a Gaussian W is the non-orthogonal version of
the same thing). The chain's carried state is a graded tower (item 3): sources indexed by birth layer, i.e. by word
length in the path expansion of the propagators, which is the Fock degree of the per-layer correspondence.

## 3. What the dictionary explains in the earlier notes

| earlier finding | reading in the common structure |
|---|---|
| the chain's closures are exact for traces and fail per neuron by incoherent O(n^{-1/2}) terms (notes XIII, XVIII; annealed against quenched, task 7) | Delta and W Delta W^T are asymptotically free (Voiculescu), i.e. a commuting square in the large-n limit for normalised traces; the per-neuron residual is the finite-n defect of that commuting square, random in sign because the relative position is random |
| the dominant identified residual is the all-distinct gate covariance, about a third of the error (closure-round section 5) | it is the "different position" term that the scattering-spectra reduction drops; Mallat may drop it because the Fourier pair is an exact commuting square with disjoint supports, we may not because the Gaussian pair is not |
| the (W o W) transport of the fourth-cumulant diagonal carries quenched content cheaply (note XV) | W o W (normalised) is the doubly stochastic block-index matrix whose distance from J/n is the quasi-orthogonality defect; the chain carries that defect exactly at first order |
| the hidden separator of dimension about 64 (note XX) | a third subalgebra (the common-mode span) relative to which the gate algebra and the transported content are conditionally in commuting position; carrying it costs 64 transports per source-layer, the price of making the square commute |
| Markov partitions and Gibbs measures (note XIX) | the path-space filtration of a Markov partition is the cylinder-set tower of C(Sigma_A) in O_A; its martingale differences are the spectral triple's eigenspaces (Farsi et al.); its Gibbs measure is the KMS state recovered from heat traces (Gerontogiannis-Goffeng) |
| He criticality and marginal error gain (note XVII) | the scattering energy-decay theorems need a gap (1 - beta < 1); at criticality the mean-direction transport has gain one, so the cascade has zero gap and errors accumulate in quadrature rather than being forgotten |
| the age-tiered shared basis (note XV) | a diffusion-wavelet ladder: the scaling spaces of the propagator products shrink with age (the Lyapunov knee near n/3 after four layers); the per-layer attenuation of transported cumulant content (0.80-0.92) plays the role of the heat semigroup e^{-t} on word length, with eigenvalues statistically linear in length as in Gerontogiannis-Goffeng; unlike a Cuntz tree the source list does not branch, so there is no critical temperature, only a geometric memory of six or seven layers |
| the V32 join (note XXIV) | item 3: the compression must be orthogonal for the state that reads the content. The shipped join projected in the flat metric of the wrong layer |

The last row is the one that paid, and it was predicted by the structure before it was measured: Farsi et al.'s
wavelets are orthonormal for the Markov measure, Kim-Kribs's conditional expectations are orthogonal projections
only for the trace, and the chain's join used neither the metric of the layer that reads its legs nor the weights
with which that layer reads them.

## 4. Mallat's two prescriptions, applied

- **Renormalise before you truncate.** Lempereur and Mallat divide each wavelet coefficient by its standard
  deviation so the conditional problem is well conditioned. Our version renormalises by observability rather than
  variance: a row of the carried legs at layer l reaches later layers through its gate Phi(alpha) and is read
  locally through D3 and D21, so the join projects in diag(Phi(alpha)^2 + c) after W. Measured below: worth 1.4
  points of raw at zero cost on top of the flat post-W projection.
- **Carry long-range interactions at coarse scales.** The coarse scale here is the common-mode subspace of note XX,
  which captures 96% of the gate-covariance term at rank 64 and is priced at 64 transports per source-layer. The
  reason it cannot be cheaper is visible in the formula of closure-round section 5: the correction contracts the
  fine (per-hub) third cumulant with the coarse correlation, D3corr_i ~ sum_e w2_e W_ie Phi_e sum_k lambda_k
  ([W diag(v_k phi) C]_ie)^2, one n^3 product per coarse direction. Mallat's hierarchy avoids this because
  stationarity makes every K_{j,j'} a convolution; a random W gives no such reduction. The hierarchic route stays
  closed at the price note XX found, and the dichotomy of closure-round section 6 (conditioning on any cheap
  statistic leaves at least 99% of the variance) rules out the sampling version.

## 5. Measurements (Azure VMs; full tables in note XXIV, `../geometry-memory/outputs/`)

| lever (from this note's structure) | nets | raw | C/B | adjusted against shipped |
|---|---|---|---|---|
| projection in the reading layer's metric (post-W), rank 320, Strassen join | 100 | +4.30% | 0.2326 | -5.45% (5.432e-9) |
| + Strassen rotations and D21 lift (arithmetic) | 16 | +3.61% | 0.2295 | -7.34% |
| + renormalised reading metric diag(Phi^2 + 0.1) | 100 | +2.98% | 0.2296 | **-7.86% (5.2925e-9)** |

## 6. What does not transfer

- The stability and energy theorems are statements about the transform, not about the closure error of its
  moments; at He criticality they give no contraction (zero gap).
- Bounded log-Sobolev constants for wavelet conditionals are a conjecture even for phi^4.
- The Cuntz-Krieger results live on autonomous dynamics (fixed isometries, eigenvalues depending only on the word);
  our transports are fresh random matrices at every layer, so only the statistical (Oseledets) version applies.
- The scattering covariance collapses to O(log^3 d) coefficients by stationarity; the neuron index has permutation
  symmetry only in law, so the analogous averages are self-averaging traces, not the per-neuron means the score
  needs.

## 7. Decision

| from the sources | what it is here | decision |
|---|---|---|
| scattering cascade, Pythagoras, Beurling-Deny decay | the network's layer is the same cascade with a random filter bank; zero gap at criticality | explains marginal error gain; no lever |
| scattering covariance, same-position reduction | the chain's slice closure is this reduction; its dropped term is the dominant residual | names the residual; the reduction's justification fails without a commuting square |
| hierarchic conditionals over coarse scales | common-mode separator of dimension 64 | priced in note XX; closed |
| renormalisation of coefficients | projection in the gate-renormalised reading metric | measured: part of the new best (-7.86%, 5.2925e-9) |
| wavelets = Cuntz representations = Laplacian eigenspaces orthogonal for the Markov/KMS state | compressions must be orthogonal for the reading state | measured: the V32 family |
| Cartan pairs, commuting squares, quasi-orthogonality, Hadamard matrices, quantum Fourier analysis | the common structure of scattering, the Cuntz-Krieger wavelets, the Markov-Gibbs picture, the separators and the chain; the network is its random-position case | the frame for task 2: a representation must either make the square commute (pay the separator) or carry the defect cheaply, as W o W does at first order |
