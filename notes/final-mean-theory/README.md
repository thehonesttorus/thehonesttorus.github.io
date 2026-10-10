# What the theory says about the final-layer mean: a proved design

Working note XLVIII. Theory only: every design choice below is either forced by a proved statement or labelled as a
free choice, with its consequence stated. No experiment is used as an argument. Two kinds of input are used:
- the coupled-observable design study (the bilateral residual identity and the regularised control lemma);
- the stage-15 wall-jet expansion.

Measured numbers that appear elsewhere in our notes are cited only as context, never as premises.

## 0. Setting

- Input x ~ N(0, I_n). Pre-activations z^0 = W_0 x and z^(l+1) = W_(l+1) h^l; post-activations h^l = relu(z^l).
- l = 0..L with L = 15, n = 1024, He weights (σ^2 = 2/n). Target: m := E h^L ∈ R^n.
- Write μ^l = E z^l, s_a^2 = Var z^l_a, r_a = μ_a/s_a, and φ^l = E|z^l| (the fold means).

## 1. The target is an exact linear image of one-dimensional data

**Theorem 1 (fold transport identity).** For every l, m^l = (1/2) W_l m^(l-1) + (1/2) φ^l, with W_0 m^(-1) := E z^0 = 0. Hence

    m^L = Σ_(l=0..L) 2^-(L-l+1) P_(L,l) φ^l,      P_(L,l) = W_L W_(L-1) ... W_(l+1),  P_(L,L) = I.

*Proof.* relu(t) = (t + |t|)/2 and E z^l = W_l m^(l-1). Induct on l. ∎

**Corollary 1.1 (exact error identity).** For any estimator written as m̂ = Σ 2^-(L-l+1) P_(L,l) φ̂^l,

    m̂ - m = Σ_l 2^-(L-l+1) P_(L,l) (φ̂^l - φ^l)

exactly. The transports P are ungated weight products, known exactly and applied as matrix-vector products.

**Corollary 1.2 (marginal sufficiency).** φ^l_a = E|z^l_a| depends only on the law of the one scalar z^l_a. No joint
object enters the target.

## 2. One-dimensional decoding has a bilateral certificate, so degree six suffices

Fix a neuron and a reference N(μ̂, ŝ^2), and write ζ = (z - μ̂)/ŝ. Define:
- h_k = E He_k(ζ), the Hermite moments of the true marginal relative to the reference;
- c_k = E_ref[g(z) He_k(ζ)], for g = |.| (or relu).

Stein's lemma gives c_k = ŝ^k E_ref g^(k):
- for relu: c_0 = ŝ(r̂Φ + φ), c_1 = ŝΦ, and c_k = ŝ He_(k-2)(-r̂) φ(r̂) for k ≥ 2;
- for |.|: twice these for k ≥ 2, and c_1 = ŝ(2Φ - 1).

These are the wall jets of stage 15.

**Theorem 2 (bilateral decoding certificate).** Suppose the density ratio of the true marginal to the reference lies in
L^2(ref). Then E g(z) = Σ_k c_k h_k / k!. A decoder that uses estimates ĥ_k for k ≤ K satisfies

    |ĝ - E g| ≤ |Σ_(k≤K) c_k (ĥ_k - h_k)/k!|  +  (Σ_(k>K) c_k^2/k!)^(1/2) (Σ_(k>K) h_k^2/k!)^(1/2).

*Proof.* Expand g in the orthogonal Hermite basis of L^2(ref) and apply Parseval. Cauchy-Schwarz on the tail gives the
second term. This is the coupled-observable study's identity τ(ρA) - τ(E_N ρ E_N A) = τ((ρ - E_N ρ)(A - E_N A)), with
N the polynomials of degree ≤ K. ∎

**Size of the tail at K = 6.** Take the reference matched to the true mean and variance, so h_1 = h_2 = 0. Then:
- h_3 = γ_3 and h_4 = γ_4 (skewness and excess kurtosis);
- h_5 = γ_5;
- h_6 = γ_6 + 10 γ_3^2;
- h_7 = γ_7 + 35 γ_3 γ_4;
- and so on.

Neuron cumulants scale as γ_j = O(n^-(j-2)/2), so Σ_(k>6) h_k^2/k! = O(n^-2). The readout tail at r̂ = 0 is
φ(0)^2 ŝ^2 Σ_(j≥3) He_(2j)(0)^2/(2j+2)!. Its terms decay like j^(-5/2), and the sum is about 2.5e-3 ŝ^2. The second term
of Theorem 2 is therefore O(ŝ n^-1) times a small constant, a few times 1e-6 per neuron. That is far below the 1e-4 per
output at which the present system errs.

**Corollary 2.1 (what the decoder must know).** Up to a remainder that is negligible at the target precision, the
fold mean of every neuron is determined by six numbers:
- the marginal's mean and variance;
- κ_3 and κ_4;
- the products κ_3^2 and κ_3 κ_4, which consistent Edgeworth order requires;

each multiplied by an explicit jet c_k(r̂). Combined with Corollary 1.1, the final error is exactly the sum of the
transported errors in these low-order marginal moments, plus the certified remainder.

**The pair version.** The same argument in L^2 of a bivariate reference applies to the pair map
K_ab = Cov(relu z_a, relu z_b). Truncating at total Hermite degree six, with the pair cumulants κ(a,a,b), κ(a,a,b,b),
κ(a,a,a,b) and their consistent products, leaves a remainder of the same negligible order. So the post-activation
covariance needs exactly the pair cumulants up to order four.

## 3. Pair and triple structure is necessary information (not a representational habit)

**Theorem 3 (fresh-row reading).** Let w ~ N(0, σ^2 I) be independent of a symmetric matrix Δ and of a symmetric
3-tensor Θ. Then

    E (w^T Δ w - σ^2 Tr Δ)^2 = 2σ^4 ‖Δ‖_F^2,      E <Θ, w^⊗3>^2 = σ^6 (6‖Θ‖_F^2 + 9‖tr Θ‖^2),

where (tr Θ)_i = Σ_j Θ_ijj.

*Proof.* Wick's theorem on six (respectively four) Gaussian factors: count the cross pairings and the pairings with one
internal pair per copy. ∎

The variance and skewness of neuron a at layer l are contacts of the previous post-activation law:
- s_a^2 = w_a^T K w_a;
- κ_3(z_a) = <κ_3(h), w_a^⊗3>.

The rows of W_l are independent of that law, so by Theorem 3 the contacts read the errors of K and κ_3(h) in Frobenius
norm. No subspace of dimension below n can serve all rows. And K^l depends on every entry of C^l through the entrywise
pair map. So:

**Corollary 3.1.** The full pair structure C^l is needed at every layer whose fold means carry weight in Corollary 1.1,
and so are the triple-structure contacts. This is a statement about information. Carrying C^l as a dense matrix is the
cheapest exact carrier, since C^l = W_l K^(l-1) W_l^T costs n^3, but it is not required by the theorem.

**On the final layer specifically.** The final error over the last weights equals an MMD between the true penultimate
law and the estimator's implied law. In the arc-cosine kernel this is exact (note XLVII §5); for the fold means it is
the kernel E|w.h||w.h'| = (2σ^2/π)|h||h'|(√(1-c^2) + c arcsin c). Theorems 2 and 3 are the per-row and the information
forms of that statement.

## 4. Consequence: the structure is forced, the error is a ledger, and the ceiling is the noise floor

Theorems 1-3 force the architecture. Every correct estimator of m must compute, at every layer:
- the mean (exact, by Theorem 1);
- the pair structure;
- κ_3 contacts in Frobenius-faithful form;
- κ_4 contacts.

It must then decode them with consistent degree-six jets. Corollaries 1.1 and 2.1 give the error exactly:

    m̂ - m = Σ_l 2^-(L-l+1) P_(L,l) [ Σ_(k≤6) (c^l_k/k!)(ĥ^l_k - h^l_k) + R_l ],     |R_l| = O(n^-1) per neuron (Theorem 2).

**What is free.** Only three things:
- (i) the representation of the pair and triple structure;
- (ii) the contraction schedule;
- (iii) where approximation error is admitted, weighted by the response of the final means (note XLVII §6).

**What is ruled out.** Fitted counterterms. By the consistency corollary of stage 15 (jets and hyperedges are one
deformation), a coefficient fitted to the output adds a term not in the ledger. The ledger is minimised by making each
ĥ_k equal its true value, not by compensating one ĥ_k against another.

**The ceiling.** If every low-order marginal moment of every layer is computed to its own next-order accuracy
(relative O(1/n)), the error is the certified remainder: per neuron O(1e-6), summed through the transports. This is at
or below the 8e-11 noise floor of our 1e9-sample truth. So the theoretical ceiling of the forced architecture lies two
orders of magnitude below the present 1.5e-8. Its distance from that ceiling is not a representation problem. It is the
inexactness of the present computation:
- fitted counterterms and a fitted κ_4 coefficient λ;
- compressed history;
- fourth-order regeneration that omits coherent classes.

## 5. The design the theorems leave: a consistent, exact order-1/n state

**State per layer** (exact to relative O(1/n)):
- μ^l, through Theorem 1 (ungated transport of fold means);
- C^l = W_l K^(l-1) W_l^T, with K from the degree-six pair map of §2;
- the third-cumulant history as sources: each fold's birth tensor carried forward by first-jet transport, uncompressed
  or compressed with a certified Frobenius loss (Theorem 3 converts that loss into contact error);
- the fourth-cumulant slices regenerated exactly from carried objects. The facet-mechanism theory (note XXXIII §2)
  states that at first order every omitted fourth-order class is generated by the third-cumulant sources, the
  covariance and the dilation charge, so the regeneration is a closed formula, not a fit.

**Decoder:** Corollary 2.1 per neuron, and the pair version for K.

**No fitted parameter appears anywhere.**

**Cost (to be proved in full in §6).** Exact first-jet transport of the history costs one n^3 product per (source
layer, current layer) pair, Σ_l l = 120 of them. Each per-pair readout (the (2,1), (3,1), (2,2) contractions) is a
bounded number of n^3 products. The total is therefore c · 120 · n^3 with c a small constant, the same order as the
present 0.22 B (which compresses at rank 320/192 yet still spends two thirds of its bill there).

## 6. What remains to prove (theory, in this order)

1. **Completeness of the order-1/n diagram list.** For each slice of each layer (D3, D21, κ_4, K22, K31, and the
   covariance corrections): every coherent diagram of the fold map at order ≤ 1/n, written in carried objects. Then a
   proof that the omitted diagrams are O(n^(-3/2)) in the Frobenius norm that Theorem 3 says the contacts read.
2. **The regeneration theorem for the fourth-order slices.** The closed formula, from sources, covariance and the
   dilation charge, with its remainder.
3. **The contraction schedule and its exact FLOP count.** Including the response-weighted compression of old history
   (note XLVII §6: effective rank 26-150 for content born at layer 12 or earlier), via Eckart-Young in the response
   metric, with Theorem 3 converting the discarded singular mass into a certified contact error.
4. **The resulting adjusted-score bound.** MSE ceiling from §4 times the cost from item 3.
