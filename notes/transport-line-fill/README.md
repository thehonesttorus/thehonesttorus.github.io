# Source as transport, and the exact line fill, placed

Working note XXIII. A further theory note was uploaded (`THEORY_correlated_boundary_response.md`, developed from
the 46-page manuscript assessed in note XXII). It derives, for an arbitrary correlated Gaussian reference with no
small-factor assumption:

- a pooled coherent form of the omitted fourth-cumulant classes and its first variation, with the moving-mean
  correction integrated inside the same expectation so that no baseline third moments are needed (Proposition 1,
  Theorem 2);
- the incoming third-cumulant source realised as the minimum-energy quadratic transport v_T = (1/3) Sigma grad s_T,
  so that E[s_T f] = E[v_T . grad f] and only first derivatives of ReLU appear (Theorem 3), with the score and
  transport representatives differing by an explicit Gaussian divergence;
- a Ward change of integrand from the full scale source, usable as a zero-mean control without replacing the
  pair source (Theorem 4);
- exact integration of one whitened Gaussian direction through all gate crossings of a single ReLU layer, by
  truncated Gaussian moments of a degree-seven polynomial per interval, compiled with prefix sums to
  O(n^2 + m n + n log n) per line (Theorems 5 and section 7);
- the conditional-integration projection P_a as an operator whose commutator with a readout measures the
  variance removed, Var f - Var P_a f = (1/4) ||[2 P_a - I, M_f] 1||^2 (section 8), and a predictable-adaptivity
  contract with the certificate ||f - P_a f||^2 >= (a . E[xi f])^2 (section 9).

The note is careful about its own scope: it is a local one-layer evaluator, its line fill does not extend to a
deep network restricted to a line, and the number of outer samples needed on the actual states is unknown. Its
small checks (a two-dimensional correlated example, finite laws, constructed width-256 algebra tests) verify the
identities, not anything about the challenge networks.

## What our measurements already say about it

**The target quantity is at the noise floor on these networks.** The evaluator computes the response of the
omitted (3,1), (2,1,1), (1,1,1,1) classes to the incoming pair-slice third-cumulant source. Note XXI section 5
decomposed the true fourth cumulant of the next pre-activation by class on the Monte Carlo law, and subtracted the
mixture's closed-form transport and the one-loop generation of the omitted classes. The remainder, which is what a
third-cumulant response into the omitted classes could be, is +0.0001, -0.0003, +0.0004 at targets 8, 11, 15 (g4
units), within the noise of the (1,1,1,1) class and at most 2% of g4. A faster or exact evaluator of that remainder
cannot move the score. The feed that drives the deep g4/g3 profile is in the retained classes: the true
third-cumulant slices feed the next fourth cumulant 2-6 times more than mixture-shaped ones (note XXI section 5),
and the chain's existing first-order programs already compute that response exactly (the ledger's retained column
matches the true pair classes to four decimals, note XIX).

**The line fill is Rao-Blackwellisation along one Gaussian direction, and its ceiling is measured.** The commutator
identity of section 8 and the certificate of section 9 say the variance removed by filling direction a is at least
the squared first-chaos coefficient of the integrand along a, and in general the variance of the integrand that the
conditional expectation in that direction explains. For the challenge quantity, the final mean over a 1024-
dimensional input, note XVIII section 6 measured every cheap conditioning: radial 0.7%, the best 256-direction frame
16%, antithetic 9%, with 99% of the variance left under any retained description that fits the budget (573x short).
One direction is at most a few parts in a thousand of that. The note's 31.6x variance reduction is in two dimensions,
where one direction is half the space; the note itself states the 1/n gap for random single-coordinate fills in n
dimensions. And for the deep network, the restriction to a line has many more than n breakpoints, which the note
also says, so the exact fill does not apply to the quantity we estimate.

**Where it would apply, it meets the same price.** The one place in the chain where a one-layer correlated-Gaussian
response is the error is the gate covariance of the all-distinct third-cumulant transport (note XVIII section 5,
note XX): about a third of the error, carried by a rank-64 common mode. A sampling evaluator of the full
gate-correlated transport would need, for every output neuron separately (the pooled form gives one coherent
scalar, not n reusable corrections, as the note says), a score-weighted expectation per source at O(n^2) per
sample; the correction is 1% of a quantity whose score-weighted estimator has variance of order its square, so the
sample count to resolve it is far above n, and n samples per source-layer already cost n^3. The rank-64 deterministic
route of note XX is cheaper and was already priced out.

## What transfers

| from the note | what it is here | decision |
|---|---|---|
| pooled omitted-class identity and moving-mean compile | correct; the same quantity as note XXI's class decomposition, summed | no use: the quantity is at the noise floor here |
| source as minimum-energy transport | a correct and useful change of representative (score multiplication to a first derivative); lower variance in the note's example | recorded; no target in the chain needs it |
| exact line fill through all gate crossings | Rao-Blackwellisation along one direction of a one-layer functional | bounded by note XVIII section 6 for the global problem; does not extend to depth |
| commutator measure of the variance a conditional fill removes | the operator form of the conditioning ledger of note XVIII section 6 | same numbers, new notation |
| Ward zero-mean control | a valid control variate for a pair-source integrand | no pair-source integrand is sampled in the chain |
| predictable adaptivity contract | correct and standard | relevant only to a sampling estimator |

The theory is sound as far as it can be checked, and it closes the formal question of how to evaluate the missing
feed at a correlated reference without a factor assumption. On these networks the missing feed it evaluates is
measured to be negligible, and the variance-reduction primitive it builds is capped by a measurement we already
have. Nothing here changes the plan: the open problem remains a cheaper representation of the per-neuron
third-cumulant transport (task 2), and the bounded levers being tested on the Azure VM.
