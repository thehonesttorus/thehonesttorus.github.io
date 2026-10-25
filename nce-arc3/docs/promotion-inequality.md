# Cross-check: the corrected library-promotion ("maternity ward") inequality after amendment E6

Setting (Theorem 2.1 (iv), Proposition 9.3). Library graph G = (V, E, s, r); primitive e has statistical dimension d_e > 1 (gauge unitaries have d_e = 1 and must be handled as weight-1 edges; see note at the end). A_β(u, v) = Σ_{e: u→v} d_e^{−β}. β_c solves r(A_{β_c}) = 1 with right eigenvector v and left eigenvector u, u·v = 1. The critical KMS (Parry) weight of a word w = e_1…e_k from s to r, conditional on starting at s, is

  ω(w) = (v_r / v_s) Π_j d_{e_j}^{−β_c}.

Fact (E6). λ(β) := r(A_β) is strictly decreasing in β (every entry decreases), and adjoining any edge with positive weight increases λ(β) at every β. So adjoining a primitive raises β_c: δβ_c = δλ / (−λ'(β_c)) > 0, with

  λ'(β_c) = −Σ_e u_{s(e)} v_{r(e)} d_e^{−β_c} log d_e  (< 0),
  δλ = u_s v_r d_new^{−β_c}   (first-order PF perturbation for a new edge s→r of weight d_new^{−β_c}; the old path is not removed, so no subtraction).

The monograph's δβ_c < 0 is the wrong sign; it also subtracted the unfactored path weight, which has no basis since the path remains available.

Criterion. The quantity the theory lets us compare is the critical-state weight (code length) of the *corpus of explanations* the agent has accumulated: W_obs = multiset of words used to explain the observed transitions (and level transitions). Define the corpus log-weight under library L,

  S_L = Σ_{w ∈ W_obs} log ω_L(w) = Σ_w [ log v_{r(w)} − log v_{s(w)} − β_c Σ_{j} log d_{e_j} ].

Promotion of the composite c = e_1…e_k (s→r) to a primitive e_new of dimension d_new ≤ Π_j d_{e_j} replaces every occurrence of c in W_obs by e_new, changes β_c to β_c' = β_c + δβ_c and (v, u) to first order. Promote iff

  S_{L'} − S_L > κ_lib.

First-order evaluation (n_c = number of occurrences of c in W_obs; C_tot = Σ_{w ∈ W_obs} Σ_j log d_{e_j} the total log-dimension mass of the corpus before replacement; C_tot' the mass after replacement):

  S_{L'} − S_L ≈ n_c β_c ( Σ_j log d_{e_j} − log d_new )        [gain: shorter, lower-dimension explanation of what was actually seen]
                 − δβ_c · C_tot'                                  [cost: the whole corpus is now coded at a higher critical temperature]
                 + Σ_w [ δ log v_{r(w)} − δ log v_{s(w)} ]       [boundary terms from the eigenvector change; zero for closed words, small otherwise]

So a compression is promoted when it is used often (n_c large) and compresses a lot (Σ log d − log d_new large), and is blocked when the alphabet's entropy-rate rise δβ_c, multiplied by how much other explanation mass there is, outweighs it. This is minimum description length in the theory's own currency, with the critical KMS state playing the role of the universal code (Shannon–McMillan–Breiman for the Parry measure). κ_lib ≥ 0 is the only free constant; κ_lib = 0 is the theory's default.

Two consequences for the engine:
1. Promotion is naturally *frequency-gated*: a word seen once never qualifies unless its compression gain alone beats the global cost, which it cannot when the corpus is non-trivial. This prevents library explosion without an ad-hoc cap.
2. Because β_c rises after each promotion, every later search is run at the recomputed β_c (T11), and the agent's prior over long words becomes *flatter in the new primitives and steeper in the old ones*: compositional depth is re-counted in the enlarged alphabet, which is exactly how the T10 depth ceiling log T / log(1/λ) is escaped.

Note on d_e = 1 edges. Gauge unitaries (D4 elements, translations) have statistical dimension 1 and contribute weight 1 at every β, so a library with a d=1 cycle has r(A_β) ≥ 1 for all β and no critical β_c exists unless the gauge group is quotiented out: words are taken modulo the gauge group (the base algebra is already the crossed product by the gauge group, Component 1), so gauge elements are not edges of G but act on the vertices' types. Designers must handle this explicitly; otherwise the Parry prior is undefined.
