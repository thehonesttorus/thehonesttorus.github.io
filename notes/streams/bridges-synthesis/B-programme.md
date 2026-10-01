# Bridges, angle B: what Gibbs–Markov theory and expanders give the barycentre-arrow programme

### Note 1's sufficiency criterion corrected, the Chen–Rouzé construction transplanted to histories, and positivity quantified as a uniform local gap

*Bridges synthesis, angle B (transfer to the programme), 2026-10-01. Written for prong 1 and prong 3 of the [research programme](../../research-program.md); every statement is about the abstract objects (faces = abstract barycentres, arrows = a ReLU followed by its linear layer, the quiver $\Lambda$, a real cocycle $F$ on arrows, states), never about activation vectors, and every one would be stated identically for a Bratteli diagram or a tiling AF filtration (drag test, rule C2). Inputs read in full before writing: the bridge digests [arxiv-2609.38007.md](../../digests/bridges/arxiv-2609.38007.md) (Yang), [arxiv-2504.02208.md](../../digests/bridges/arxiv-2504.02208.md) (Chen–Rouzé), [hammersley-clifford.md](../../digests/bridges/hammersley-clifford.md), [expanders.md](../../digests/bridges/expanders.md) (Tao; Hoory–Linial–Wigderson, the probable `expander_survey.pdf`), [hdx-spectral-independence.md](../../digests/bridges/hdx-spectral-independence.md), [nc-dirichlet-lindblad.md](../../digests/bridges/nc-dirichlet-lindblad.md), [transfer-spectrum-measurement.md](../../digests/bridges/transfer-spectrum-measurement.md), [chat-2609.38007-retrieval-status.md](../../digests/bridges/chat-2609.38007-retrieval-status.md); the programme notes [conditional-arrow-algebra.md](../../conditional-arrow-algebra.md) (Note 1), [simplicial-complex-as-decomposition.md](../../simplicial-complex-as-decomposition.md) (Note 2), [local-to-global-unlocks.md](../../local-to-global-unlocks.md), [mlp-bridge.md](../../mlp-bridge.md); for context [competition-plan.md](../../competition-plan.md) §§3.1, 6b, 7 and the stream reports. Nothing below redefines an object of Notes 1–2; corrections travel as messages (rule C1, §11).*

**Labels.** Every claim carries exactly one label.
- **THEOREM**: a published theorem cited to the digest section where it is quoted (and proved or checked there), or a standard fact whose proof I give.
- **KNOWN-LINK**: a published connection between two areas, cited.
- **DERIVED**: proved here; the proof is in the text. Most DERIVED finite statements were also checked numerically (§12; script in the appendix); those that were not are short identities whose proof is given.
- **ANALOGY**: a structural similarity, with the exact point where it breaks.
- **CONJECTURE**: a precise statement, not proved.
- **SPECULATION**: a direction, not yet a statement.

These map onto the programme's labels as Fact = THEOREM / KNOWN-LINK, Derivation = DERIVED, Conjecture = CONJECTURE, Interpretation = ANALOGY / SPECULATION.

**What could not be read.** The shared chat behind arXiv:2609.38007 was not retrieved (eleven routes failed; [chat-2609.38007-retrieval-status.md](../../digests/bridges/chat-2609.38007-retrieval-status.md) §1). Nothing here describes it. The local file `expander_survey.pdf` is not in the container; [expanders.md](../../digests/bridges/expanders.md) §0 uses Hoory–Linial–Wigderson, served under that file name, as the stand-in.

---

## 1. Results at a glance

| # | statement | label | § |
|---|---|---|---|
| R1 | Note 1 §6's "the pair (input face, output face) is sufficient iff $F$ is a coboundary" is false at depth 1: on $\{a,b\}\to\{c,d\}$ fully joined, every $F$ makes the pair sufficient, and only a 3-dimensional subspace of the 4-dimensional space of cochains consists of coboundaries | DERIVED | 3.1 |
| R2 | Corrected statement: $F$ is a coboundary iff the pair is sufficient **and** the induced function $\Phi(s,r)$ on the reachability graph $R\subseteq K_0\times K_L$ is a coboundary there. Note 1's equivalence holds for every $F$ iff the bigons span the cycle space, $p(B)=Z_1(\Gamma)$; this holds when $R$ is a forest (one input face, as in dictionary v1, or one output face) and when each component of $R$ has a waist face (any complete layered quiver of depth $\ge2$) | DERIVED | 3.2–3.3 |
| R3 | The endpoint statistic depends on how the Gibbs family is normalised. For the globally normalised family Note 1 is right. With fixed origin weights $w$ (Note 1 §5.2's corner KMS family) and two or more input faces, the endpoint face is sufficient iff $F=\delta U$ **and** the law of $U(r(\mu))$ under uniform counting of the histories from $s$ is the same for any two origins sharing an endpoint. "$U=0$ on the input layer" is neither necessary nor sufficient there | DERIVED | 3.4 |
| R4 | "The face is Markov" and "the face is sufficient for the Gibbs family" are both Petz sufficiency of a face subalgebra, for two different families (the past as parameter, the temperature as parameter). They are logically independent: every Gibbs law is Markov, few are sufficient | DERIVED | 4 |
| R5 | For re-routing jumps inside a depth window $W$, reversibility of the jump process with respect to $P$ is exactly the KMS$_\beta$ quasi-invariance of $P$ with Radon–Nikodym cocycle $e^{-\beta c_F}$. The Cesàro mean of the window dynamics converges to the window's conditional expectation, which recovers $P$ from its restriction outside $W$, iff the jumps connect each window class (for square-flip jumps: iff the window fibres are flip-connected). For Gibbs laws this recovery is exact and has zero radius | DERIVED | 5 |
| R6 | **Positivity quantified gives a uniform local gap.** For a positive chain, the projective (cross-ratio) diameter $\Delta_\ell$ of each layer kernel is invariant under every conditioning. The single-site window sampler has absolute spectral gap $\ge(1-2\tanh(\Delta/4))/w$ for every window length $w$ and boundary condition. For every $\Delta<\infty$, block samplers of length $b>e^{\Delta/2}-1$ have gap $\ge(b-(e^{\Delta/2}-1))/(w+b-1)$. The single-site case is Dobrushin's condition on a path run through the Dyer–Goldberg–Jerrum path coupling, with the classical projective-diameter bound on the influences | DERIVED (single-site case: KNOWN-LINK, Dobrushin–DGJ) | 6.2 |
| R7 | Forgetting with a local certificate: the KMS state's restriction to the first $t$ layers depends on the terminal condition at most through $\prod_{\ell\ge t}\tanh(\Delta_\ell/4)$ in total variation (kernel $A_\ell$ joins layers $\ell,\ell+1$). In the stationary case the normalised transfer operator has $\lvert\lambda_2\rvert\le\tanh(\Delta/4)$. Both are classical contraction statements (Birkhoff). This has the shape of U1 of the local-to-global note but is not U1: U1's link-spectral hypothesis plays no role, and its literal conclusion is trivial for Markov laws | DERIVED; reading as U1: ANALOGY | 6.3 |
| R8 | Quantitative sufficiency: $\mathrm{Var}_\beta(F\mid s,r)\le S^2/(2\,\mathrm{gap})$, with $S$ the largest square defect. The relative-entropy sufficiency defect is $\le\frac{(\beta'-\beta)^2}2\max\mathrm{Var}$, and the exact Fisher-information loss is $\mathbb E_\beta\mathrm{Var}_\beta(F\mid s,r)$ | DERIVED | 6.4 |
| R9 | **Two halves of one class.** On a complete layered quiver, $H^1(\Gamma)$ is an extension of the layerwise zigzag part $\bigoplus_\ell H^1(\Gamma_\ell)$ by a mismatch part $M\cong\bigoplus_{0<\ell<L}\mathbb R^{K_\ell}/\mathbb R$. Mixing (the projective diameters) sees only the zigzag image, $\Delta_\ell=\lvert\beta\rvert Z_\ell$. Sufficiency sees the whole class. Vanishing square defects imply vanishing zigzags ($Z_\ell\le2S$), but rank-one transfer (perfect mixing) does not imply sufficiency | DERIVED | 6.5 |
| R10 | U3 in $\ell^2$ form: if $H^1(X_\square)=0$, then $\mathrm{dist}_2(F,B^1)\le\lVert\delta_\square F\rVert_2/\sqrt{\mu_1}$ and the fibre oscillation is $\le2\sqrt L\,\lVert\delta_\square F\rVert_2/\sqrt{\mu_1}$. U3's own hypothesis already forces $p(B)=Z_1(\Gamma)$. Measured: $\mu_1=n^2$ for complete layered quivers of width $n$, at every depth tried | DERIVED; depth-independence CONJECTURE | 7 |
| R11 | Yang's Lemma III.4 holds for any subalgebra with a trace-preserving conditional expectation. In the commutative (Gibbs-on-histories) resolution it is dominated by R8: it loses at least one power of $\lvert\beta'-\beta\rvert$ | THEOREM + DERIVED | 8.1–8.2 |
| R12 | A conjecture and a speculation say where Chen–Rouzé and Yang would carry real content: CMI decay across a depth buffer for coherent (non-diagonal) KMS states, and a noncommutative uniform local gap | CONJECTURE (8.4); SPECULATION (8.5) | 8.3–8.5 |
| R13 | Bounded *unpinned* spectral independence $\eta_0$ of a state on the face algebra gives exactly one thing: every linear observable $\sum_ua_ue_u$, the face dimension included, has variance at most $(1+\eta_0)$ times its product-state value. Mixing needs all pinnings (U6), and ALO's product bound needs in addition a cap $\eta_i\le\theta(n-i-1)$, $\theta<1$, on the deepest pinnings: a uniform $\eta_i\le\eta$ alone does not give $n^{-(1+\eta)}$ (counterexample in §9). With the cap, at $\eta\approx7$ and $n=1024$ the bound is $10^{-21}$–$10^{-34}$: numerically vacuous | DERIVED (Prop 9.1); THEOREM (ALO, with the cap) | 9 |

The verdict on the user's intuition, in the programme's terms, is §10.

---

## 2. Objects

Everything is finite. Paths are written forward, left to right (Note 1 writes them right to left; nothing depends on this).

**D1 (layered quiver, histories, trimming).** Layers $K_0,\dots,K_L$ are finite sets of faces. Arrows $\gamma\in\Lambda_\ell$ go from $s(\gamma)\in K_\ell$ to $r(\gamma)\in K_{\ell+1}$, and parallel arrows are allowed. A *complete history* is $\mu=(\gamma_1,\dots,\gamma_L)$ with $r(\gamma_i)=s(\gamma_{i+1})$, $s(\mu)=s(\gamma_1)\in K_0$ and $r(\mu)=r(\gamma_L)\in K_L$; its face sequence is $\sigma_0,\dots,\sigma_L$. $\Lambda$ is *trimmed* if every face and arrow lies on a complete history. Trimming discards nothing the Gibbs states see. With a root adjoined, this is exactly a finite Bratteli diagram (Note 1 §3.4). The *fibre* of $(s,r)$ is the set $H_{s,r}$ of complete histories from $s$ to $r$.

**D2 (cocycle and Gibbs families).**
- $F:E^1\to\mathbb R$, with $F(\mu)=\sum_iF(\gamma_i)$.
- $\Gamma$ is the underlying undirected multigraph; $C_1(\Gamma)$, $Z_1(\Gamma)$ are its chains and cycles.
- $F$ is a *coboundary* if $F=\delta U$, $(\delta U)(\gamma)=U(r\gamma)-U(s\gamma)$. Over $\mathbb R$ this holds iff $F$ vanishes on $Z_1(\Gamma)$, and $[F]\in H^1(\Gamma;\mathbb R)$.

Two Gibbs families appear in Note 1 §5.2:
- the **fixed-origin-weight family** $P^w_\beta(\mu)=w_{s(\mu)}e^{-\beta F(\mu)}/Z_{s(\mu)}(\beta)$, $w>0$ fixed. These are the KMS$_\beta$ states of the complete-history corner $\bigoplus_{\sigma_0}M_{N_L(\sigma_0)}$;
- the **globally normalised family** $P^{\rm gl}_\beta(\mu)=e^{-\beta F(\mu)}/Z(\beta)$, the unique KMS$_\beta$ state of the algebra of all complete histories ("the pair-groupoid algebra" of Note 1 §5.2). This is the family that `notes/checks/graph_algebra_checks.py` actually tests: its function `sufficient` asks whether $F(\mu)$ is a function of the statistic.

**D3 (bigons, reachability graph, zigzag obstruction).**
- $p:\mathbb R^{\{\text{complete histories}\}}\to C_1(\Gamma)$ sends $e_\mu$ to the chain $\sum_ie_{\gamma_i}$.
- $\partial_h e_\mu=[r(\mu)]-[s(\mu)]$.
- $B=\mathrm{span}\{e_\mu-e_{\mu'}:(s,r)(\mu)=(s,r)(\mu')\}\subseteq\ker\partial_h$, and $p(B)\subseteq Z_1(\Gamma)$ is the *bigon space*.
- The *reachability graph* $R$ is the simple bipartite graph on $K_0\sqcup K_L$ with an edge $\{s,r\}$ iff $H_{s,r}\neq\emptyset$.
- The *zigzag obstruction space* is $O(\Lambda):=Z_1(\Gamma)/p(B)$.

**D4 (squares, flips, the square complex).**
- A *square* $\square$ is a pair of length-2 paths with common ends, $\sigma\to\tau\to\rho$ and $\sigma\to\tau'\to\rho$. Its defect is $\delta F(\square)=F(\sigma\tau)+F(\tau\rho)-F(\sigma\tau')-F(\tau'\rho)$, the alternating sum of U3 in [local-to-global-unlocks.md](../../local-to-global-unlocks.md) §5.
- Two histories in a fibre are joined by a *flip* if they differ at one inner position, i.e. in $(\gamma_i,\sigma_i,\gamma_{i+1})$. A fibre is *flip-connected* if its flip graph is connected.
- The *square complex* $X_\square$ is the 2-complex with 1-skeleton $\Gamma$ and one 2-cell per square. On a trimmed quiver every square lies on a complete history.

**D5 (layer transfer matrices, projective diameter).**
- $A_\ell(\sigma,\tau)=\sum_{\gamma\in\Lambda_\ell:\sigma\to\tau}e^{-\beta F(\gamma)}$ for $\sigma\in K_\ell$, $\tau\in K_{\ell+1}$. This is Note 1 §5.2's weighted matrix, one layer at a time.
- The layer pair is *positive* if $A_\ell>0$ entrywise, i.e. completely joined.
- For a positive matrix $A$,
$$
\Delta(A)=\max_{\sigma,\sigma',\tau,\tau'}\log\frac{A(\sigma,\tau)A(\sigma',\tau')}{A(\sigma,\tau')A(\sigma',\tau)}\in[0,\infty),
$$
its *projective diameter* (the diameter of its image in Hilbert's projective metric; here only the cross-ratio definition is used).

---

## 3. Note 1 §6, corrected

### 3.1 The counterexample

**DERIVED (R1).** Take $L=1$, $K_0=\{a,b\}$, $K_1=\{c,d\}$, one arrow for each of the four pairs, and $F=1$ on $a\to c$, $0$ elsewhere.
- Every fibre $H_{s,r}$ is a single arrow, so the pair $(s,r)$ is sufficient for both families (the conditional law given the pair is a point mass).
- $F$ is not a coboundary, because $F(ac)-F(ad)-F(bc)+F(bd)=1\neq0$, while every coboundary $U(r\gamma)-U(s\gamma)$ gives $0$ on this alternating 4-cycle.

So Note 1 §6's second bullet ("the pair is sufficient iff $F$ is a coboundary with arbitrary $U$") fails. Its forward direction (coboundary ⇒ pair sufficient) is true, and it is the only direction Note 1 §9 checked. The counterexample was first reported in [expanders.md](../../digests/bridges/expanders.md) §11.1 item 5. It is verified here (§12, C1), and the next theorem says exactly what replaces the statement.

### 3.2 The corrected statement

**DERIVED (Lemma 3.1).** If $\Lambda$ is trimmed, then $p(\ker\partial_h)=Z_1(\Gamma)$, and $q:e_\mu\mapsto e_{(s(\mu),r(\mu))}$ induces an isomorphism $\ker\partial_h/B\cong Z_1(R)$.

*Proof.*
- **$\subseteq$.** $\partial_\Gamma p(e_\mu)=[r(\mu)]-[s(\mu)]=\partial_he_\mu$, so $p$ maps $\ker\partial_h$ into cycles.
- **$\supseteq$, setup.** Use trimming to choose, for each face $\rho$, a prefix $\alpha_\rho$ (a path from $K_0$ to $\rho$, empty if $\rho\in K_0$) and a suffix $\beta_\rho$ (from $\rho$ to $K_L$). For an arrow $\gamma:\sigma\to\tau$, both $\mu_\gamma:=\alpha_\sigma\gamma\beta_\tau$ and $\nu_\tau:=\alpha_\tau\beta_\tau$ are complete histories. Their chains satisfy
$$
e_\gamma=p(e_{\mu_\gamma})-p(e_{\nu_\tau})+\bar\alpha_\tau-\bar\alpha_\sigma ,
$$
where $\bar\alpha_\rho$ denotes the chain of $\alpha_\rho$.
- **$\supseteq$, conclusion.** For a cycle $z=\sum n_\gamma e_\gamma$, the terms $\sum_\gamma n_\gamma(\bar\alpha_{r\gamma}-\bar\alpha_{s\gamma})=\sum_\rho\bar\alpha_\rho\,(\partial z)(\rho)$ vanish. So $z=p(y)$ with $y=\sum_\gamma n_\gamma(e_{\mu_\gamma}-e_{\nu_{r\gamma}})$, and $\partial_hy=\partial_\Gamma z=0$.
- **The quotient.** $q$ maps onto $C_1(R)$, and $\ker q=B$. Since $\partial_Rq=\partial_h$, $q(\ker\partial_h)=Z_1(R)$: any preimage of a cycle of $R$ lies in $\ker\partial_h$. $\square$

**DERIVED (Theorem 3.2).** Let $\Lambda$ be trimmed and let $\beta$ take at least two values.
- **(a) Pair sufficiency.** For either family of D2, the pair $(s,r)$ is a sufficient statistic iff $F(\mu)=\Phi(s(\mu),r(\mu))$ for some function $\Phi$ on the edges of $R$, iff $F$ vanishes on $p(B)$.
- **(b) Coboundary.** $F$ is a coboundary iff (a) holds and $\Phi$ is a coboundary on $R$, i.e. $\Phi(s,r)=V(r)-W(s)$.
- **(c) When the two coincide.** "Pair sufficient ⇒ coboundary" holds for every $F$ on $\Lambda$ iff $p(B)=Z_1(\Gamma)$, i.e. iff $O(\Lambda)=0$. Always $\dim O(\Lambda)\le\dim Z_1(R)$.

*Proof.*
- **(a)** For $\mu,\mu'$ in one fibre, $P(\mu)/P(\mu')=e^{-\beta(F(\mu)-F(\mu'))}$ in both families: origin weights and normalisations cancel. So the conditional law on a fibre is the Gibbs law $\propto e^{-\beta F}$ there, and it is the same at two values of $\beta$ iff $F$ is constant on the fibre. In chain language that reads $\langle F,p(e_\mu-e_{\mu'})\rangle=0$.
- **(b)** $F$ is a coboundary iff it vanishes on $Z_1(\Gamma)=p(\ker\partial_h)$ (Lemma 3.1), i.e. iff the functional $y\mapsto\sum_\mu y_\mu F(\mu)$ vanishes on $\ker\partial_h$. That happens iff it vanishes on $B$ (which is (a)) and the induced functional on $\ker\partial_h/B\cong Z_1(R)$, namely $z\mapsto\sum_{(s,r)}z_{(s,r)}\Phi(s,r)$, vanishes. The latter says $\Phi\in B^1(R)$.
- **(c)** The pair-sufficient cochains are the annihilator of $p(B)$, and the coboundaries are the annihilator of $Z_1(\Gamma)$. By Lemma 3.1, $O(\Lambda)=p(\ker\partial_h)/p(B)$ is a quotient of $\ker\partial_h/B\cong Z_1(R)$. $\square$

*Remark.* [expanders.md](../../digests/bridges/expanders.md) §11.1 states the extra condition as "$\Phi(\sigma_0,\tau)-\Phi(\sigma_0',\tau)$ independent of $\tau$". That rectangle condition is (b) when $R$ is complete bipartite; in general the condition is the cycle condition on $R$.

### 3.3 When Note 1's equivalence does hold

**DERIVED (Corollary 3.3).** $O(\Lambda)=0$, so pair sufficiency ⇔ coboundary, in each of the following cases:
- **(i) $R$ is a forest.** In particular there is a single input face (dictionary v1's trivial layer-0 face `*`, [mlp-bridge.md](../../mlp-bridge.md) §1) or a single output face.
- **(ii) Each connected component of $R$ has a waist.** A waist is a face $\tau$ at an inner layer that every origin of the component reaches and that reaches every output of the component.
- **(iii) $\Lambda$ is completely layered and $L\ge2$.** Every inner face is a waist.

In the other direction, at depth 1 without parallel arrows, $p(B)=0$ and $O(\Lambda)=H_1$ of the bipartite graph. Every quiver of depth 1 with a cycle is a counterexample.

*Proof.*
- (i) $Z_1(R)=0$ forces $\ker\partial_h=B$, so $Z_1(\Gamma)=p(B)$.
- (ii) Let $F$ be pair-sufficient, and pick paths $\alpha_s:s\leadsto\tau$ and $\beta_r:\tau\leadsto r$ for the origins and outputs of the component. Then $\alpha_s\beta_r\in H_{s,r}$, so $\Phi(s,r)=F(\alpha_s)+F(\beta_r)$ is separable, and (b) applies.
- (iii) follows from (ii). $\square$

*Checked* (§12, C2–C3), on random trimmed quivers with parallel arrows:
- $p(\ker\partial_h)=Z_1(\Gamma)$ and $\dim\ker\partial_h/B=\dim Z_1(R)$ held in all 243 cases;
- "coboundary ⇔ $\Phi$ coboundary on $R$" held for every pair-sufficient $F$;
- the equivalence failed in 26 of 243 quivers, and in each of them a generic pair-sufficient $F$ was not a coboundary;
- every quiver with a forest $R$, or with a waist (270/270), satisfied it.

### 3.4 The endpoint statistic and the normalisation of the family

**DERIVED (Theorem 3.4).** Let $\Lambda$ be trimmed and let $\beta$ range over an open interval.
- **(a) Global family.** The endpoint face is sufficient for $\{P^{\rm gl}_\beta\}$ iff $F(\mu)=V(r(\mu))$ on complete histories, iff $F=\delta U$ with $U$ constant on the origins of each component of $R$. This is Note 1 §6's first bullet, which is correct for this family.
- **(b) Fixed origin weights.** Fix $w>0$. The endpoint face is sufficient for $\{P^w_\beta\}$ iff $F=\delta U$ for some $U$ and, whenever two origins $s,s'$ share an endpoint, $U\circ r$ has the same law under the uniform counting measure on complete histories from $s$ as under the one from $s'$. The values of $U$ on $K_0$ play no role.

*Proof.*
- **(a)** Sufficiency means $F(\mu)-F(\mu')=0$ for co-terminal $\mu,\mu'$ (the normalisation is global, so it cancels). Then $F$ is pair-sufficient with $\Phi(s,r)=V(r)$ separable, and Theorem 3.2(b) applies. Conversely, $F=\delta U$ with $U\equiv c$ on the origins of a component gives $F(\mu)=U(r\mu)-c$, and each output lies in one component.
- **(b), setup.** Write $Z_s(\beta)=\sum_{s(\mu)=s}e^{-\beta F(\mu)}$ and $\mathbb E_{\beta,s}$ for the Gibbs mean on histories from $s$. For co-terminal $\mu$ (from $s$) and $\mu'$ (from $s'$),
$$
\frac{d}{d\beta}\log\frac{P^w_\beta(\mu)}{P^w_\beta(\mu')}=-(F(\mu)-F(\mu'))+\mathbb E_{\beta,s}F-\mathbb E_{\beta,s'}F .
$$
Sufficiency means this vanishes on the interval.
- **(b), step 1.** With $s=s'$ it gives pair sufficiency, $F(\mu)=\Phi(s,r)$.
- **(b), step 2.** For $s\neq s'$ sharing $\tau$: $\Phi(s,\tau)-\Phi(s',\tau)=m_s(\beta)-m_{s'}(\beta)$ with $m_s=\mathbb E_{\beta,s}F$. The left side does not depend on $\beta$ and the right side does not depend on $\tau$. Fix $\beta_0$ and put $c(s)=-m_s(\beta_0)$. Then $V(\tau):=\Phi(s,\tau)+c(s)$ does not depend on which origin $s$ reaching $\tau$ is used, since $\Phi(s,\tau)-\Phi(s',\tau)=(m_s-m_{s'})(\beta_0)=c(s')-c(s)$. So $\Phi(s,\tau)=V(\tau)-c(s)$, and $F=\delta U$ with $U|_{K_L}=V$ and $U|_{K_0}=c$.
- **(b), step 3.** With this $U$, $\mathbb E_{\beta,s}F=\mathbb E_{\beta,s}[U(r)]-c(s)$, and the condition of step 2 becomes $\mathbb E_{\beta,s}U(r)=\mathbb E_{\beta,s'}U(r)$ on the interval. Equivalently, $\frac d{d\beta}\log\frac{M_s(\beta)}{M_{s'}(\beta)}=0$ with $M_s(\beta)=\sum_{s(\mu)=s}e^{-\beta U(r\mu)}$. By analyticity, $M_s/N_s=M_{s'}/N_{s'}$ for all $\beta$, where $N_s=M_s(0)$ counts the histories from $s$. Equal Laplace transforms of finitely supported laws mean equal laws.
- **(b), converse.** If $F=\delta U$ and the laws agree, the ratio above equals $\frac{w_s}{w_{s'}}\frac{M_{s'}(\beta)}{M_s(\beta)}=\frac{w_sN_{s'}}{w_{s'}N_s}$, a constant. $\square$

*Examples* (§12, C5).
- Depth 1, origins $s,s'$ each joined to $t_1,t_2$, $U(t_1)=0$, $U(t_2)=1$, $U(s)=0$, $U(s')=5$: sufficient for $P^w$ and **not** for $P^{\rm gl}$.
- Origins $s\to\{t_1,t_2\}$ and $s'\to\{t_1\}$, $F=\delta U$, $U|_{K_0}=0$, $U(t_1)\neq U(t_2)$: sufficient for $P^{\rm gl}$ and **not** for $P^w$. (With $U(t_1)=U(t_2)$ it is sufficient for both; re-checked in the critique pass.)
- On 263 random quivers with random coboundaries, the criterion (b) predicted sufficiency correctly in all 263.

### 3.5 Consequences

**DERIVED (Corollary 3.5).**
- (i) Note 1 §6's chain "face sufficient ⇔ $F\sim0$ ⇔ $\mathbb P_\beta$ central ⇔ tracial on $A_L$ ⇔ modular flow trivial" holds verbatim for $P^{\rm gl}$ (with "$F\sim0$" read as $F=\delta U$, $U$ constant on the input faces of each component). It also holds for $P^w$ when there is one input face.
- (ii) For $P^w$ with two or more input faces, endpoint sufficiency does not imply centrality. Under (b), $P^w_\beta(\mu)=(w_s/M_s(\beta))e^{-\beta U(r\mu)}$, so histories into $\tau$ from different origins carry different weights unless $w_s\propto M_s(\beta)$.

*Proof.* (i) is Theorem 3.4(a) together with Note 1's own argument. (ii) is read off the displayed formula. $\square$

What this changes elsewhere (rule C4):
- the realised quantities of [mlp-bridge.md](../../mlp-bridge.md) E3 used a single input face, so its readings are unaffected;
- [expanders.md](../../digests/bridges/expanders.md) §11.1's remark that the first bullet "is unaffected" is true for $P^{\rm gl}$ and false for $P^w$.

---

## 4. The Markov property and the sufficiency of the barycentre: one Petz question, two families

**DERIVED (Proposition 4.1).** Let $P$ be any law on complete histories, Gibbs or not, and $0<\ell<L$. Take any reference law $R_{<\ell}$ of full support on pasts. The following are equivalent:
- (i) $P$ is Markov at $\ell$: past $\perp$ future given $\sigma_\ell$.
- (ii) The statistic $\sigma_\ell$ is sufficient (Fisher–Neyman) for the family $\{P(\sigma_\ell,\text{future}\mid\text{past}=\pi)\}_\pi$, in which the parameter is the past.
- (iii) $I_P(\text{past}:\text{future}\mid\sigma_\ell)=0$.
- (iv) The subalgebra of functions of (past, $\sigma_\ell$) is Petz-sufficient for the pair $\{P,\ R_{<\ell}\otimes P_{\sigma_\ell,\text{future}}\}$. Its sufficiency defect is
$$
D\big(P\,\big\|\,R_{<\ell}\otimes P_{\sigma_\ell,\rm fut}\big)-D\big(P_{\rm past,\sigma_\ell}\,\big\|\,R_{<\ell}\otimes P_{\sigma_\ell}\big)=I_P(\text{past}:\text{future}\mid\sigma_\ell).
$$

*Proof.*
- (i) ⇔ (ii) is the definition of sufficiency, and (i) ⇔ (iii) is standard.
- The identity in (iv) follows by expanding the four relative entropies. The $\log R_{<\ell}$ terms cancel, and what remains is $H(\sigma_\ell,\mathrm{fut})-H(\mathrm{all})+H(\mathrm{past},\sigma_\ell)-H(\sigma_\ell)$. This is the classical form of [arxiv-2504.02208.md](../../digests/bridges/arxiv-2504.02208.md) §6.3, there with $R$ uniform.
- Zero defect is sufficiency (Petz; Jenčová–Petz Thm 1 as quoted in Note 1 §6). $\square$

**DERIVED (Corollary 4.2).** For the Gibbs families, (i) holds at every layer, for every $F$ and every $\beta$ (THEOREM: Note 1 §5.3, where $P_\beta$ is the Doob chain), while endpoint sufficiency for $\{P_\beta\}$ holds only as in Theorem 3.4. So "the face is sufficient for the past" and "the face is sufficient for the scale $\beta$ of the cocycle" are independent. Both are instances of Jenčová–Petz, applied to two different families. They are exactly the "two different failures of sufficiency" of [mlp-bridge.md](../../mlp-bridge.md) §3.2: memory, and a non-coboundary kernel.

**THEOREM ([hdx-spectral-independence.md](../../digests/bridges/hdx-spectral-independence.md) §6.4, U7; checked there, C2; for three variables it is Ibinson–Linden–Winter's identity, [hammersley-clifford.md](../../digests/bridges/hammersley-clifford.md) §5.3(i)).** $D(P\Vert Q_{\rm Markov})=\sum_iI(\sigma_{i+1};\sigma_{<i}\mid\sigma_i)$, and this equals the relative-entropy distance of $P$ from the closure of the KMS family. With Proposition 4.1, measured memory is therefore a sum of Petz sufficiency defects, one per layer.

---

## 5. Window samplers: KMS is detailed balance, and the Chen–Rouzé map on histories

Fix a window $W=\{a+1,\dots,b-1\}$ of inner layers. The *window re-routing relation* $G_W$ consists of the pairs $(\mu',\mu)$ of complete histories that agree everywhere except in the faces of $W$ and the arrows incident to them. It carries the cocycle $c_F(\mu',\mu)=F(\mu')-F(\mu)$. $G_W$ is a finite equivalence relation; its classes are the window fibres, indexed by the data outside $W$.

**DERIVED (Proposition 5.1, R5).** Let a jump process on complete histories have rates $k(\mu\to\mu')=a(\mu,\mu')\,m(\beta c_F(\mu',\mu))$. Here $a=a^{\rm T}\ge0$ is supported on $G_W$, and $m>0$ satisfies $m(x)=e^{-x}m(-x)$ (for instance Metropolis, $m(x)=\min(1,e^{-x})$, or heat-bath, $m(x)=(1+e^x)^{-1}$). A probability $P$ that is positive on each class:
- is reversible for $k$ iff $P(\mu')/P(\mu)=e^{-\beta c_F(\mu',\mu)}$ for every jump with $a(\mu,\mu')>0$;
- if the $a$-graph connects each class of $G_W$, this holds iff $P$ is quasi-invariant on $G_W$ with Radon–Nikodym cocycle $e^{-\beta c_F}$, i.e. iff $P$ restricts to a KMS$_\beta$ state of $(C^*(G_W),\alpha^{c_F})$. The KMS–quasi-invariance correspondence is the THEOREM of Renault and Neshveyev quoted in Note 1 §5.2.

*Proof.* Reversibility is $P(\mu)k(\mu\to\mu')=P(\mu')k(\mu'\to\mu)$, i.e. $P(\mu')/P(\mu)=m(x)/m(-x)=e^{-x}$ with $x=\beta c_F(\mu',\mu)$. The cocycle identity carries this from the jump edges to the whole class. $\square$

**DERIVED (same proof).** If the windows cover the inner layers and every fibre is flip-connected, the measures reversible for every window sampler are exactly the $\psi(s,r)e^{-\beta F(\mu)}$, because the window relations then generate the relation "same origin and same endpoint". These are the KMS$_\beta$ states of the relation "same origin and same endpoint", with the boundary weights free.

**DERIVED (Proposition 5.2).** Let $\mathcal L$ be such a reversible generator and $\mathcal R_t=t^{-1}\int_0^te^{s\mathcal L}ds$.
- (a) $\mathcal R_t\to\Pi$, where $\Pi f(\mu)=\mathbb E_P[f\mid\text{jump-component of }\mu]$.
- (b) $\Pi=\mathbb E_P[\,\cdot\mid\text{outside }W]$ iff every class of $G_W$ is connected by jumps. In that case, for **every** law $Q$ with $Q|_{\text{outside }W}=P|_{\text{outside }W}$, one has $Q\Pi=P$: exact recovery.
- (c) If $P=P_\beta$, then $\mathbb E_P[\,\cdot\mid\text{outside }W]$ depends only on the two faces $\sigma_a,\sigma_b$ adjacent to $W$. The recovery map is local, with zero radius.
- (d) [DERIVED] $\mathcal E(\mathcal R_tf)\le0.41\,t^{-1}\lVert f\rVert^2$ with no gap: with $x\ge0$ the spectral variable of $-\mathcal L$ in $L^2(P)$, $\mathcal E(\mathcal R_tf)=\int x\big(\tfrac{1-e^{-tx}}{tx}\big)^2d\mu_f(x)\le t^{-1}\sup_u\tfrac{(1-e^{-u})^2}{u}\lVert f\rVert^2$. The same derivation is in [expanders.md](../../digests/bridges/expanders.md) §9.9 and is checked in [nc-dirichlet-lindblad.md](../../digests/bridges/nc-dirichlet-lindblad.md) K4; [CR] Cor. VII.1 is the operator-norm version, $\le2/t$. The constant is $\sup_u(1-e^{-u})^2/u\approx0.4073$.

*Proof.*
- (a) Reversibility makes $\mathcal L$ self-adjoint in $L^2(P)$. Then $\mathcal Lf=0$ iff $\mathcal E(f,f)=0$ iff $f$ is constant across every jump. The mean ergodic theorem does the rest.
- (b) The classes of $G_W$ are the level sets of the outside data, and $\Pi(\mu,\mu')=P(\mu'\mid\text{class})$. Hence $(Q\Pi)(\mu')=Q(\text{class})P(\mu'\mid\text{class})=P(\mu')$, because $Q$ and $P$ give each class the same mass. This is the commutative form of the argument of [arxiv-2504.02208.md](../../digests/bridges/arxiv-2504.02208.md) §6.1.
- (c) is the Markov property of $P_\beta$ (Note 1 §5.3). $\square$

**ANALOGY (the Chen–Rouzé construction, transplanted).** The "time-averaged detailed-balanced Lindbladian based on single-Pauli jumps on $A$" ([arxiv-2504.02208.md](../../digests/bridges/arxiv-2504.02208.md) §1.4) corresponds, ingredient by ingredient, to the window sampler. The last column marks where the analogy breaks.

| [CR] ingredient | programme analogue | status | where it breaks |
|---|---|---|---|
| erase $A$: $\rho_{-A}=\mathrm{Tr}_A\rho\otimes\tau_A$ | any $Q$ agreeing with $P$ outside $W$ | DERIVED (5.2b) | — |
| single-site Paulis generate $B(\mathcal H_A)$ | square flips generate $G_W$ iff the window fibres are flip-connected | DERIVED (5.2b) | a quiver whose window fibres are not flip-connected (the two-route example of [expanders.md](../../digests/bridges/expanders.md) §11.1) makes recovery *fail*, not merely slow |
| KMS detailed balance, exact thanks to a coherent term | reversibility ⇔ KMS quasi-invariance on $G_W$ (5.1); the re-routing jumps $e_{\mu\nu}$ are eigenvectors of the modular operator of a state diagonal in the history basis, so Alicki's construction is GNS-, hence KMS-, detailed-balanced ([nc-dirichlet-lindblad.md](../../digests/bridges/nc-dirichlet-lindblad.md) §6.1). Diagonality alone would not give GNS: every state is diagonal in its eigenbasis, and KMS-but-not-GNS generators exist on $M_2$ (Carlen–Maas App. B, ibid. §2.2) | DERIVED | no coherent correction is needed: the modular flow $\alpha^F_{-\beta t}$ preserves every window algebra ([arxiv-2609.38007.md](../../digests/bridges/arxiv-2609.38007.md) B10) |
| Cesàro mean, Dirichlet form $\le2/t$ with no gap | the same (5.2d) | DERIVED | unnecessary under positivity: §6 gives exponential convergence of the heat-bath window dynamics |
| quasi-locality from Lieb–Robinson | exact locality: $\Pi$ sees only $\sigma_a,\sigma_b$ (5.2c) | THEOREM (Markov) | the main difficulty of [CR] is absent |
| $e^{\mu\lvert A\rvert}$ from Pauli word length and Gibbs conjugation | the flip-graph diameter, and the block length $b>e^{\Delta/2}-1$ of Theorem 6.4 | ANALOGY | a different mechanism: a projective-diameter cost (temperature × zigzag class), not a word-length cost |
| a uniform local gap would give global Markov (Cor. B.2) | positivity ⇒ a uniform local gap (Thm 6.4), forgetting (Cor 6.5) and quantitative sufficiency (Thm 6.7); global Markov is free | DERIVED | in the programme the gap is a theorem under positivity, not a hypothesis |

**DERIVED (summary of 5.1–5.2 and §6).** In the diagonal KMS family of Note 1 the construction is **exact and local, and its content reduces to the gap.** **SPECULATION.** Its genuinely noncommutative content (an approximate conditional expectation where no exact one exists) would appear only once the state has coherences between histories (§8).

---

## 6. Positivity quantified: a uniform local gap, forgetting, and quantitative sufficiency

This section is the programme's answer to "a trickle-down theorem for local gaps of KMS samplers is the missing piece". **DERIVED (Theorem 6.4).** In the commutative resolution the piece is not missing. It is supplied by a one-dimensional path-coupling argument whose only certificate is the projective diameter of each layer kernel, which is Hammersley–Clifford's positivity hypothesis made quantitative.

**KNOWN-LINK (the ingredients are classical).** For $b=1$, Theorem 6.4 is Dobrushin's uniqueness condition on a path graph, run through the Dyer–Goldberg–Jerrum path coupling ([expanders.md](../../digests/bridges/expanders.md) §9.2: contraction $1-(1-\lVert R\rVert)/n$ per step; §11.1 there already proposes this reduction for the re-routing chain), with the influences of a site on its two neighbours bounded by Lemma 6.2. Lemma 6.2 is the classical bound of the Dobrushin coefficient by Birkhoff's contraction coefficient $\tanh(\Delta/4)$ (from memory). Lemma 6.3 is the contraction-to-gap theorem that [expanders.md](../../digests/bridges/expanders.md) §11.1 cites from memory (Chen 1998; Levin–Peres–Wilmer Thm 13.1). The proofs are given in full below. What is the programme's own is the reading of the certificate as the projective diameter of a layer kernel, and the use made of it in §§6.3–6.6.

Setting for §§6.1–6.4: a *positive chain* $x_0,\dots,x_m$ with finite state spaces and law
$$P(x)\propto w(x_0)\prod_{\ell=0}^{m-1}A_\ell(x_\ell,x_{\ell+1})\,v(x_m),\qquad A_\ell>0 .$$
Write $\Delta_\ell=\Delta(A_\ell)$, $t_\ell=\tanh(\Delta_\ell/4)$ and $\delta=\max_\ell t_\ell$. The face process of $P_\beta$ on a positive quiver is such a chain, with $A_\ell$ as in D5.

### 6.1 Three lemmas

**DERIVED (Lemma 6.1).** Let $p,q$ be probability vectors with $p_i/q_i\in[c,ce^D]$. Then $\mathrm{TV}(p,q)\le\tanh(D/4)$, and the bound is attained.

*Proof* (also in [nc-dirichlet-lindblad.md](../../digests/bridges/nc-dirichlet-lindblad.md) §6.1(f)).
- Put $\rho=p/q\in[m,M]$ with $\mathbb E_q\rho=1$ and $M=me^D$.
- On $[m,M]$, $(\rho-1)_+\le(M-1)(\rho-m)/(M-m)$ (a chord above a convex function). So $\mathrm{TV}=\mathbb E_q(\rho-1)_+\le(M-1)(1-m)/(M-m)$.
- Over $m\in[e^{-D},1]$ this is largest at $m=e^{-D/2}$, where it equals $(e^{D/2}-1)/(e^{D/2}+1)=\tanh(D/4)$. $\square$

**DERIVED (Lemma 6.2: conditioning does not change the certificate).** For positive $h,f$ and positive $A$, define
- the forward Doob kernel $P(x,y)=A(x,y)h(y)/\sum_{y'}A(x,y')h(y')$,
- the backward kernel $B(y,x)=f(x)A(x,y)/\sum_{x'}f(x')A(x',y)$.

Both have Dobrushin coefficient $\max_{x,x'}\mathrm{TV}(\text{row }x,\text{row }x')\le\tanh(\Delta(A)/4)$.

*Proof.* The ratio of two rows of $P$ is a constant times $A(x,\cdot)/A(x',\cdot)$, whose multiplicative range is at most $e^{\Delta(A)}$ by the definition of $\Delta$. The same holds for $B$. Apply Lemma 6.1. $\square$

**DERIVED (from Lemma 6.2).** Every conditioning of a positive chain (pinning sites, fixing boundary values, changing $w$ or $v$) gives positive chains with the same kernels $A_\ell$ and new boundary weights, whose transition kernels are Doob transforms. So **the certificate is hereditary for free.** This is the quantitative form of "pinnings of Markov fields are Markov fields" ([hammersley-clifford.md](../../digests/bridges/hammersley-clifford.md) §5.4).

**DERIVED (Lemma 6.3: contraction ⇒ gap; path coupling).** Let $Q$ be reversible on a finite set carrying the path metric $d$ of a graph. If $W_1(Q(x,\cdot),Q(y,\cdot))\le1-\kappa$ for adjacent $x,y$, then every eigenvalue $\lambda\neq1$ of $Q$ satisfies $\lvert\lambda\rvert\le1-\kappa$.

*Proof.*
- Along a geodesic, the triangle inequality for $W_1$ gives $W_1(Q(x,\cdot),Q(y,\cdot))\le(1-\kappa)d(x,y)$ for all $x,y$.
- Let $Qf=\lambda f$ with $\lambda\neq1$. Then $f\perp1$ in $L^2(\pi)$, so $f$ is not constant and $\mathrm{Lip}(f)>0$.
- $\lvert\lambda\rvert\,\lvert f(x)-f(y)\rvert=\lvert Qf(x)-Qf(y)\rvert\le\mathrm{Lip}(f)\,W_1\le\mathrm{Lip}(f)(1-\kappa)d(x,y)$. Maximise over $x,y$. $\square$

### 6.2 The uniform local gap

**DERIVED (Theorem 6.4, R6).** Pin $x_0$ and $x_m$ to arbitrary values; the free sites form the window $\{1,\dots,m-1\}$ of size $w=m-1$. For $b\ge1$, the *block heat-bath sampler* $Q_b$ picks $t$ uniformly from $\{2-b,\dots,m-1\}$ ($w+b-1$ choices) and resamples $x_{B_t}$, $B_t=\{t,\dots,t+b-1\}\cap\{1,\dots,m-1\}$, from its conditional law. It is reversible, and every eigenvalue $\lambda\neq1$ satisfies
$$
\lvert\lambda\rvert\ \le\ 1-\kappa_b,\qquad \kappa_b=\frac{b-2\sum_{k=1}^{b}\delta^k}{w+b-1}.
$$
In particular:
- **Single site.** For $b=1$ the bound refines to $\kappa_1=\min_j\big[1-t_{j-1}\mathbf 1_{j\ge2}-t_j\mathbf 1_{j\le w-1}\big]/w$. It is positive whenever $\Delta<2\log3$, i.e. $\tanh(\Delta/4)<\frac12$.
- **Any $\Delta<\infty$.** $\kappa_b>0$ as soon as $b>e^{\Delta/2}-1$, because $\sum_{k\ge1}\delta^k=\delta/(1-\delta)=(e^{\Delta/2}-1)/2$.

The bound is uniform in the window length $w$ (as $c/w$) and in the boundary condition. Internal pinnings only improve it.

*Proof.*
- **Setup.** Each block update is a conditional expectation, hence $P$-self-adjoint, and $Q_b$ is their average. Use the Hamming metric on the free sites. Let $x,y$ differ only at site $j$, and couple the two chains by choosing the same $t$.
- **Blocks containing $j$.** There are exactly $b$ such $t$. These blocks see identical boundary values, so the coupled resamples coincide and the distance drops to $0$.
- **The block just right of $j$** ($t=j+1$). Its conditional law is a bridge chain started from $x_j$ resp. $y_j$ with the same forward Doob kernels. By Lemma 6.2 they contract total variation by $\le\delta$ per step. A step-by-step maximal coupling (identical once equal) makes the $k$-th block site disagree with probability $\le\delta^k$. So the expected number of new disagreements is $\le\sum_{k=1}^b\delta^k$, while site $j$ still differs.
- **The block just left of $j$** ($t=j-b$). The same, using the backward kernels of Lemma 6.2.
- **Every other $t$.** Identical laws; the distance stays $1$.
- **Assembly.** $\mathbb E\,d'\le\big[2(1+\sum_k\delta^k)+(w+b-1-b-2)\big]/(w+b-1)=1-\kappa_b$. For $b=1$, keep the separate $t_{j-1},t_j$. Lemma 6.3 finishes. $\square$

*Checked* (§12, G2–G3): the absolute gap was $\ge\kappa_1$ on all 47 random chains with $\kappa_1>0$ (smallest ratio 1.001), and $\ge\kappa_b$ for blocks at $\Delta\in[2.4,3.3]$, beyond the single-site threshold. Re-checked in the critique pass with independent code: 238/238 random chains (smallest ratio 1.0001), adversarial two-state chains near $\Delta=2\log3$ (smallest ratio 1.025), and 28/28 block cases.

### 6.3 Forgetting: U1 as a theorem

**DERIVED (Corollary 6.5, R7).** Let $P^v$ and $P^{v'}$ differ only in the terminal weights. For every $t$,
$$
\mathrm{TV}\big(P^v|_{x_0,\dots,x_t},\,P^{v'}|_{x_0,\dots,x_t}\big)\ \le\ \prod_{\ell=t}^{m-1}\tanh(\Delta_\ell/4),
$$
and symmetrically the law of the late layers forgets the origin weights.

*Proof.*
- Read backwards, $P^v$ is a Markov chain started from its law at $x_m$, with kernels $B_\ell(x_{\ell+1},x_\ell)\propto f_\ell(x_\ell)A_\ell(x_\ell,x_{\ell+1})$, where $f_\ell$ are the forward weights from $w$. These kernels do not involve $v$.
- Each $B_\ell$ contracts total variation by $\le t_\ell$ (Lemma 6.2).
- The joint law of $x_0,\dots,x_t$ is the law of $x_t$ pushed through further kernels that are common to both chains, so its total variation is at most that of $x_t$. $\square$

*Reading for the programme* (DERIVED for the forgetting statement; ANALOGY for its identification with U1). In Note 1 §5.3 the Doob walk depends on the output condition only through the backward partition function $Z$, and that dependence decays at the rate $\tanh(\Delta_\ell/4)$ per layer, with a certificate that is local (one layer pair at a time). This is the *shape* U1 of [local-to-global-unlocks.md](../../local-to-global-unlocks.md) §5 asks for, but it is not U1. *Where it breaks:* U1's hypothesis (one-sided local spectral expansion of the layer and history complexes) plays no role here, and U1's literal conclusion ("the KMS state is determined up to $\lambda^t$ by the statistics of histories of length $t$") is exact and trivial for KMS states, which are Markov and hence determined by their consecutive-pair statistics. What holds is forgetting of boundary conditions at rate $\lambda=\max_\ell\tanh(\Delta_\ell/4)$, a classical contraction statement. U1 as stated stays open.

**DERIVED (Corollary 6.6, the stationary Perron gap).** If $A_\ell=A>0$ for all $\ell$ (the repeated pattern of U5), then $\lvert\lambda_2(A)\rvert/\rho(A)\le\tanh(\Delta(A)/4)$.

*Proof.*
- Let $P=D_h^{-1}AD_h/\rho(A)$ with $h$ the Perron vector. Then $\delta(P)\le\tanh(\Delta/4)$ (Lemma 6.2).
- For real measures $\nu$ of mass zero, $\lVert\nu P\rVert_1\le\delta(P)\lVert\nu\rVert_1$. Split $\nu=\nu_+-\nu_-$: the two parts have equal mass, and each pair of rows differs by $\le2\delta(P)$ in $\ell^1$.
- On complex zero-mass measures, $\lVert P^n\rVert\le2\delta^n$. Gelfand's formula then bounds the spectral radius on that subspace by $\delta$.
- Left eigenvectors with $\lambda\ne1$ have zero mass, and $\mathrm{spec}(P)=\mathrm{spec}(A)/\rho(A)$. $\square$

This is the total-variation form of the Birkhoff–Hopf contraction, which [transfer-spectrum-measurement.md](../../digests/bridges/transfer-spectrum-measurement.md) §6 lists as a fact from memory (KNOWN-LINK, from memory and not opened: Birkhoff's theorem is the contraction in Hilbert's projective metric with the same coefficient $\tanh(\Delta/4)$, and the eigenvalue bound follows from it). *Drag test, the Penrose diagram* (DERIVED, checked in §12):
- the Fibonacci incidence matrix $A=\begin{pmatrix}1&1\\1&0\end{pmatrix}$ has a zero, but $A^2=\begin{pmatrix}2&1\\1&1\end{pmatrix}>0$ has $\Delta=\log2$;
- so the central measure of Connes's Penrose AF algebra is approached from any terminal condition at rate $\le\tanh(\log2/4)=0.172$ per two levels;
- the exact rate is $\varphi^{-4}=0.146$.

### 6.4 Quantitative sufficiency

Setting: a single-arrow positive quiver, or more generally a positive chain with $A_\ell=e^{-u\,g_\ell}$ and the additive functional $F(x)=\sum_\ell g_\ell(x_\ell,x_{\ell+1})$. Pin both ends, $x_0=a$ and $x_m=c$; the fibre law is $P_u\propto e^{-uF}$ for either family of D2. Let $S$ be the largest change of $F$ under a single-site change (on a quiver: the largest square defect $\lvert\delta F(\square)\rvert$).

**DERIVED (Theorem 6.7, R8).**
- (a) *Exact.* The Fisher information about $u$ lost by keeping only the ends is $I_{\rm full}(u)-I_{\rm ends}(u)=\mathbb E_u\,\mathrm{Var}_u(F\mid\text{ends})$.
- (b) $\mathrm{Var}_u(F\mid a,c)\le S^2/(2\,\mathrm{gap}(Q_1))$, and $\le S^2/(2\kappa_1(u))$ when $\kappa_1(u)>0$, where $\kappa_1(u)$ is computed from $\Delta(e^{-ug_\ell})=\lvert u\rvert Z_\ell$, with $Z_\ell$ the largest alternating 4-cycle (zigzag) sum of $g_\ell$.
- (c) For $\beta<\beta'$, the relative-entropy sufficiency defect satisfies
$$\delta_{\rm ends}(\beta,\beta'):=\mathbb E_{P_\beta}D\big(P_\beta(\cdot\mid\text{ends})\,\Vert\,P_{\beta'}(\cdot\mid\text{ends})\big)\le\tfrac{(\beta'-\beta)^2}{2}\sup_{u\in[\beta,\beta']}\max_{\rm ends}\mathrm{Var}_u(F\mid\text{ends}).$$
- (d) The defects in (a) and (c) vanish iff $F$ is constant on every fibre; when the fibres are flip-connected (for instance when every layer pair is positive) this is $S=0$.

*Proof.*
- (a) In an exponential family the score is $-(F-\mathbb E F)$. Apply the law of total variance.
- (b) Use the Poincaré inequality $\mathrm{Var}\le\mathcal E(F,F)/\mathrm{gap}$ with $\mathcal E(F,F)=\tfrac12\sum\pi(x)Q(x,y)(F(x)-F(y))^2\le S^2/2$, then Theorem 6.4.
- (c) On a fibre, $D(P_\beta\Vert P_{\beta'})=\int_\beta^{\beta'}(\beta'-u)\,\mathrm{Var}_u(F)\,du$. Average over the ends.
- (d) The laws are positive on each fibre, so the variance vanishes iff $F$ is constant there; along a flip path $F$ changes only by square defects. $\square$

*Checked* (§12, G4): $\mathrm{Var}\le S^2/(2\,\mathrm{gap})$ held on every chain tested, with the largest ratio 0.175. The gain over the trivial bound $\mathrm{osc}^2/4\le(wS)^2/4$ is a factor of order $w$: the variance grows linearly in depth, not quadratically, as for a sum along a mixing chain. Near $\kappa_1\to0$ the Poincaré form can be worse than the trivial bound; use the smaller of the two.

### 6.5 Complete layered quivers: one certificate, two halves of one class

**DERIVED (Proposition 6.8, R9).** Let $\Lambda$ be completely layered with single arrows, $n_\ell=\lvert K_\ell\rvert$, $L\ge2$, and let $\Gamma_\ell$ be the complete bipartite graph of layer pair $\ell$.
- **(a)** There is an exact sequence
$$0\to M\to H^1(\Gamma)\to\textstyle\bigoplus_\ell H^1(\Gamma_\ell)\to0,\qquad M\cong\bigoplus_{0<\ell<L}\mathbb R^{K_\ell}/\mathbb R .$$
The middle map is restriction. $M$ is spanned by the *mismatch functions* $m_\ell=a_{\ell-1}-b_\ell$ of layerwise coboundaries $F|_{\Gamma_\ell}=a_\ell(\tau)-b_\ell(\sigma)$. The dimensions are $\sum_\ell(n_\ell-1)(n_{\ell+1}-1)$ and $\sum_{0<\ell<L}(n_\ell-1)$.
- **(b)** $\Delta(A_\ell)=\lvert\beta\rvert Z_\ell$, where $Z_\ell=\max\lvert F(\sigma\tau)+F(\sigma'\tau')-F(\sigma\tau')-F(\sigma'\tau)\rvert$ is a norm on $H^1(\Gamma_\ell)$. The projective diameters, and hence the certificates of Theorem 6.4 and Corollary 6.5, see only the image of $[F]$ in $\bigoplus_\ell H^1(\Gamma_\ell)$. The exact gap does depend on the mismatch part, which enters the fibre law as single-site fields $e^{-\beta m_\ell(\sigma_\ell)}$ (§12, M: with every $\Delta_\ell$ fixed, the fibre gap moved from 0.133 to 0.250 as the fields grew).
- **(c)** $Z_\ell\le2S$, because each zigzag is the difference of two squares through a neighbouring layer pair.
- **(d)** $S=0$ ⇔ $F$ is a coboundary ⇔ pair sufficiency. Endpoint sufficiency for $P^{\rm gl}$ needs in addition $U$ constant on the input layer (Theorem 3.4(a)), so with two or more input faces it is strictly stronger than $S=0$.
- **(e)** $Z\equiv0$ ⇔ every $A_\ell$ has rank one; then the faces of different layers are independent under $P^w_\beta$ and under every pinning: perfect one-step decorrelation, a transfer kernel that is a conditional expectation onto the constants. Perfect decorrelation is **strictly weaker** than (d). If $Z\equiv0$, then $\mathrm{Var}_\beta(F\mid\text{ends})=\sum_{0<\ell<L}\mathrm{Var}_{\propto e^{-\beta m_\ell}}(m_\ell)$, which is non-zero as soon as one mismatch is non-constant.
- **(f)** If $\lvert\beta\rvert S<\log3$, then
$$\mathrm{gap}\ \ge\ \frac{1-2\tanh(\lvert\beta\rvert S/2)}{L-1},\qquad \mathrm{Var}_\beta(F\mid s,r)\ \le\ \frac{(L-1)S^2}{2\big(1-2\tanh(\lvert\beta\rvert\max_\ell Z_\ell/4)\big)} .$$

*Proof.*
- (a) Layer pairs share no arrows, so restriction is onto. A class restricting to $0$ on every $\Gamma_\ell$ is represented by a layerwise coboundary, unique up to one constant per pair (each $\Gamma_\ell$ is connected). Modulo global coboundaries what remains is the mismatches modulo constants. The dimensions add up to $\dim Z_1(\Gamma)=\sum n_\ell n_{\ell+1}-\sum n_\ell+1$.
- (b) The 4-cycles generate $H^1(\Gamma_\ell)$, and the cross-ratios of $e^{-\beta F}$ are exactly $e^{\mp\beta(\text{zigzag})}$.
- (c) Use layer pair $\ell+1$ if $\ell\le L-2$, else $\ell-1$, as in §3.3.
- (d) Squares make complete fibres flip-connected, and Corollary 3.3(iii) applies.
- (e) When zigzags vanish, the square defects are differences $m_\ell(\tau)-m_\ell(\tau')$, the face process is a product, and $F=\text{const}(s,r)+\sum_\ell m_\ell(\sigma_\ell)$.
- (f) Theorem 6.4 with $\Delta\le2\lvert\beta\rvert S$, then Theorem 6.7(b). $\square$

*Checked* (§12, C4): squares span $Z_1$ and $Z_\ell\le2S$ held on complete quivers of shapes $[2,2,2]$ through $[4,3,2,3,2]$.

**This is the precise sense in which sufficiency and expansion are two halves of one mechanism.** One cochain, $\delta F$ evaluated on the 4-cycles of the quiver, carries both:
- its *zigzag* part (cycles alternating inside one layer pair) is a projective diameter, so it governs the certificates for mixing, forgetting and the local gap (the exact gap also feels the mismatch part, as single-site fields, (b));
- the whole class, including the *mismatch* part, governs sufficiency;
- the bound (f) puts the whole class in the numerator and only the zigzag part in the denominator.

### 6.6 Sparse quivers: coarse-graining, and parallel arrows

**DERIVED (Proposition 6.9).** Suppose $r\mid L$ and the quiver is *$r$-primitive*: every product $A_{kr}A_{kr+1}\cdots A_{kr+r-1}$ is positive. Parallel arrows are allowed.
- (i) The coarse face process $x_k=\sigma_{kr}$ is a positive chain with kernels $A^{(r)}_k$, so Theorem 6.4 and Corollaries 6.5–6.6 apply to it.
- (ii) Write $\mathrm{seg}_k$ for the part of the history between layers $(k-1)r$ and $kr$. Then
$$\mathrm{Var}_\beta(F\mid s,r)=\sum_k\mathbb E\,\mathrm{Var}_\beta\big(F(\mathrm{seg}_k)\mid x_{k-1},x_k\big)+\mathrm{Var}_\beta\Big(\sum_k\bar g_k(x_{k-1},x_k)\,\Big|\,s,r\Big),$$
with $\bar g_k(\sigma,\tau)=\mathbb E_\beta[F(\mathrm{seg}_k)\mid\sigma,\tau]$. The first sum is bounded by the oscillation of $F$ over $r$-step bigons; the second by the Poincaré step of Theorem 6.7(b) on the coarse chain (Theorem 6.4 for its gap), with $S$ replaced by the largest single-site change of $\sum_k\bar g_k$. Since $\bar g_k$ is not the coarse chain's energy, the Fisher-information reading of Theorem 6.7(a) does not transfer.
- (iii) With $r=1$ and parallel arrows, the first sum is the length-1 bigon (parallel-arrow) variance. In that case $\Delta(A_\ell)$ is $\lvert\beta\rvert$ times the zigzag norm of the *face-pair free energy* $\bar F_\beta(\sigma,\tau)=-\beta^{-1}\log\sum_{\gamma:\sigma\to\tau}e^{-\beta F(\gamma)}$. As $\beta\to0$, $\Delta(A_\ell)$ tends to the largest $\lvert\log\rvert$ of the cross-ratios of the multiplicities, a purely combinatorial (entropic) zigzag.

*Proof.*
- (i) Given the coarse faces, the segments are independent (Markov property), and the coarse kernel is the product of the layer kernels.
- (ii) is the law of total variance, using that independence.
- (iii) follows from the definitions. $\square$

### 6.7 What §6 does and does not do

- **ANALOGY (trickle-down).** In the Anari–Liu–Oveis Gharan complex $X_\mu$ of the face process (in its multi-spin form, Chen–Liu–Vigoda; parts = layers), links are pinnings ([hdx-spectral-independence.md](../../digests/bridges/hdx-spectral-independence.md) §3.2). By Lemma 6.2, every pinned chain carries the same cross-ratios, so heredity costs nothing. Theorem 6.4 then replaces trickle-down by direct one-dimensional path coupling. The single-site threshold $\tanh(\Delta/4)<\tfrac12$ has the same form as the top-link hypothesis $\gamma\le\tfrac12$ of Oppenheim's trickle-down ([local-to-global-unlocks.md](../../local-to-global-unlocks.md) §2.2); its exact counterpart is classical, Dobrushin's condition $t_{j-1}+t_j<1$ on the total influence on a site (§6.2). *Where it breaks:* trickle-down acts on all links of a high-dimensional complex and needs spectral data on them. Here the depth direction is one-dimensional, which is why the elementary argument suffices. The identification of the top-link eigenvalue with $\tanh(\Delta/4)$ would need Birkhoff–Hopf's bound on maximal correlation, quoted from memory only, so it is not claimed.
- **Scope: positivity.** If $A_\ell$ has a zero (a sparse layer pair), then $\Delta_\ell=\infty$ and the statements of §6 do not apply at that resolution, though Proposition 6.9 may apply at a coarser one. This is Hammersley–Clifford's positivity condition ([hammersley-clifford.md](../../digests/bridges/hammersley-clifford.md) §1.2, §4) in quantitative form: positivity in Besag's sense makes the support a product (a safe symbol, the weaker hypothesis the Hammersley–Clifford proof actually uses, only makes it closed under switching sites to the vacuum, [hammersley-clifford.md](../../digests/bridges/hammersley-clifford.md) §1.2), and a finite projective diameter makes it a product with bounded distortion. **ANALOGY** (the reading of [hammersley-clifford.md](../../digests/bridges/hammersley-clifford.md) §13 B3, whose inequality is that digest's own elementary theorem, not a published link): positivity versus a quantitative gap is the same hypothesis at two resolutions; there it appears as the Friedrichs angle of the two-block sampler, here as the projective diameter. *Where it breaks:* B3 bounds a two-block intersection defect on a general graph by a Friedrichs angle, Theorem 6.4 bounds a many-site gap on a path by projective diameters, and no theorem here relates the two constants.
- **Low temperature. DERIVED (from Theorem 6.4).** $\Delta=\lvert\beta\rvert Z$ grows with $\beta$, so the block length needed grows like $e^{\Delta/2}$. At every finite $\beta$ the block sampler with $b=\lceil e^{\Delta/2}\rceil$ has gap $>1/(w+e^{\Delta/2})$, so what deteriorates is the block length (the cost of one exact block update) and, once $w\lesssim e^{\Delta/2}$, the gap itself; the single-site certificate is lost beyond $\Delta=2\log3$. **ANALOGY:** this matches the expander guard that low temperature is where expansion fails ([expanders.md](../../digests/bridges/expanders.md) §9.7). *Where it breaks:* in one dimension there is no phase coexistence, so the gap never becomes exponentially small in the size, as it does on expanders.

---

## 7. The square complex: U3 in $\ell^2$ form

**DERIVED (Theorem 7.1, R10).** Let $\Lambda$ be trimmed, use unweighted inner products, and let $\delta_\square:C^1\to C^2(X_\square)$ evaluate a cochain on squares. Suppose $H^1(X_\square;\mathbb R)=0$, and let $\mu_1>0$ be the least eigenvalue of $\delta_\square^*\delta_\square$ on $(B^1)^\perp$. Then:
- (a) $\mathrm{dist}_2(F,B^1)\le\lVert\delta_\square F\rVert_2/\sqrt{\mu_1}$.
- (b) For $\mu,\mu'$ in one fibre, $\lvert F(\mu)-F(\mu')\rvert\le2\sqrt L\,\lVert\delta_\square F\rVert_2/\sqrt{\mu_1}$. Hence $\mathrm{Var}(F\mid s,r)\le L\lVert\delta_\square F\rVert_2^2/\mu_1$ under any law on the fibre.
- (c) $H^1(X_\square)=0$ ⇔ squares span $Z_1(\Gamma)$ ⇒ $p(B)=Z_1(\Gamma)$. So U3's coboundary-expansion hypothesis *contains* the hypothesis that restores Note 1 §6 (§3.3).

*Proof.*
- (a) Write $F=P_{B^1}F+G$ with $G\in(B^1)^\perp=\ker\delta^*$. Since squares are cycles, $\delta_\square\delta U=0$, so $\delta_\square F=\delta_\square G$. $H^1(X_\square)=0$ makes $\delta_\square$ injective on $(B^1)^\perp$, so $\lVert\delta_\square G\rVert^2\ge\mu_1\lVert G\rVert^2$.
- (b) In $F(\mu)-F(\mu')$ the part $\delta U$ cancels, and $\lvert G(\mu)\rvert\le\sqrt L\lVert G\rVert_2$ by Cauchy–Schwarz over the $L$ distinct arrows of $\mu$.
- (c) $H^1=\ker\delta_\square/B^1$, and squares $\subseteq p(B)\subseteq Z_1$. $\square$

*Checked* (§12, H): the oscillation bound held with ratio $\le0.453$ over 200 random cochains.

**CONJECTURE (7.2).** For the complete layered quiver of uniform width $n$ and every depth $L\ge2$, $\mu_1(X_\square)=n^2$, with one 2-cell per unordered pair of length-2 paths (D4) and unweighted inner products (with ordered pairs every value doubles).
- Measured: $\mu_1=4.000$ for $n=2$ and $9.000$ for $n=3$, at every $L=2,\dots,8$. Re-measured in the critique pass with independent code: $16.000$ for $n=4$ ($L=2,\dots,5$) and $25.000$ for $n=5$ ($L=2,\dots,4$).
- For nine non-uniform shapes (depth 2–4) $\mu_1$ took values from 6 to 16, depending on the widths; depth-independence was tested only for uniform widths.

If true, the $\ell^2$ coboundary expansion of the square complex is depth-uniform, and Theorem 7.1 gives a sufficiency defect linear in $L$ from square data alone, with no gap hypothesis. This complements Theorem 6.7, which needs positivity but is pointwise in $S$.

**SPECULATION (7.3, a Garland-type bound for square complexes; not yet a precise statement, since "a function of local data" is not specified).** For trimmed layered quivers, $\mu_1(X_\square)$ is bounded below by a function of local data: the spectra of the links of faces in $X_\square$ (the graph of squares through a face) and the flip-connectivity of length-2 fibres. *What must be true:* a localisation identity of Garland type for cube-like 2-complexes ([hdx-spectral-independence.md](../../digests/bridges/hdx-spectral-independence.md) §2.1 gives the simplicial one). This would make U3 local in the sense prong 3 asked for; it is the open half of U3.

---

## 8. The coherent resolution: where Chen–Rouzé and Yang would actually be needed

### 8.1 Yang's localisation lemma for an arbitrary subalgebra

**THEOREM** ([arxiv-2609.38007.md](../../digests/bridges/arxiv-2609.38007.md) §5.4–5.5, Yang Lemma III.4 and Cor. III.5). For faithful $\rho,\sigma$, the tensor subalgebra $B(\mathcal H_R)\otimes1$ and $0<\alpha\le1$,
$$\delta_R\le(\tfrac1\alpha+3)[\mathrm{Tr}\rho^{1+\alpha}\sigma^{-\alpha}]^{1/(1+\alpha)}\mathfrak q_R^{2\alpha/(1+\alpha)},$$
where $\mathfrak q_R=\frac12\int\eta_R(t)\,dt/\lvert\sinh\pi t\rvert$ is the modular-time leakage of the Connes cocycle $u_t=\sigma^{it}\rho^{-it}$.

**DERIVED (extension to subalgebras; numerically checked in [chat-2609.38007-retrieval-status.md](../../digests/bridges/chat-2609.38007-retrieval-status.md) §7.1 and [arxiv-2609.38007.md](../../digests/bridges/arxiv-2609.38007.md) B7).** The same inequality holds for any unital $*$-subalgebra $N\subseteq M_n$, with $\mathsf E_N$ its trace-preserving conditional expectation, $\rho_N=\mathsf E_N\rho$, $\sigma_N=\mathsf E_N\sigma$, $\delta_N=D(\rho\Vert\sigma)-D(\rho_N\Vert\sigma_N)$, and $\eta_N(t)=\lVert u_t-\mathsf E_Nu_t\rVert_\infty$.

*Proof.* Yang's argument uses the tensor structure only through four facts. Each holds for $W_N(Y)=Y\rho_N^{-1/2}\rho^{1/2}$ on $L^2(N,\mathrm{Tr})$:
- (i) $W_N$ is an isometry with $W_N\rho_N^{1/2}=\rho^{1/2}$: $\mathrm{Tr}(\rho\,\rho_N^{-1/2}Y^*Y\rho_N^{-1/2})=\mathrm{Tr}(\mathsf E_N(\rho)\rho_N^{-1/2}Y^*Y\rho_N^{-1/2})=\mathrm{Tr}(Y^*Y)$, by bimodularity and trace preservation;
- (ii) $W_N^\dagger\Delta W_N=\Delta_N$: $\langle W_NY,\Delta W_NY'\rangle=\mathrm{Tr}(\sigma Y'\rho_N^{-1}Y^*)=\mathrm{Tr}(\sigma_NY'\rho_N^{-1}Y^*)$, because $Y'\rho_N^{-1}Y^*\in N$;
- (iii) $O\rho^{1/2}=W_N(O\rho_N^{1/2})$ for $O\in N$, so the trial vector $\mathsf E_N(B_s)\rho^{1/2}$ lies in $\mathrm{Ran}\,W_N$, and $\lVert B_s-\mathsf E_NB_s\rVert\le\mathfrak q_N/s$ because $\mathsf E_N$ is linear and fixes scalars;
- (iv) $\lVert\rho_N^{1/2}\rVert=1$ and $\langle\rho_N^{1/2},\Delta_N\rho_N^{1/2}\rangle=\mathrm{Tr}\,\sigma_N=1$, which is all the Jensen step needs.

The bound $\lVert(\Delta+s)^{1/2}O\rho^{1/2}\rVert^2\le(1+s)\lVert O\rVert^2$ and the small-resolvent step use no subalgebra structure. $\square$

This is a quantitative, one-directional Petz theorem (cocycle almost in $N$ ⇒ sufficiency almost holds) for the inclusions of the programme: $D_0\subset D$, a finite stage of a Bratteli diagram, a tiling patch algebra.

### 8.2 In the commutative resolution it is the wrong tool

**DERIVED (Remark 8.2, R11).** Take $\rho=P_\beta$, $\sigma=P_{\beta'}$ on complete histories and $N$ = functions of the pair $(s,r)$.
- Then $u_t=e^{it(\beta-\beta')F(\mu)}$ times a phase depending only on the origin.
- If $F$ is non-constant on some fibre, $\eta_N(t)\sim\lvert t\rvert\,\lvert\beta-\beta'\rvert\max_\mu\lvert F(\mu)-\bar F_{\rm fibre}\rvert$ as $\beta'\to\beta$, so $\mathfrak q_N\asymp\lvert\beta-\beta'\rvert$ and the bound is of order $\lvert\beta-\beta'\rvert^{2\alpha/(1+\alpha)}$, with exponent $\le1$.
- The exact defect is $\int_\beta^{\beta'}(\beta'-u)\,\mathbb E\mathrm{Var}_u\,du\asymp(\beta'-\beta)^2$ (Theorem 6.7).

So the leakage criterion of [arxiv-2609.38007.md](../../digests/bridges/arxiv-2609.38007.md) B8 and the leakage form of quantitative U3 in [chat-2609.38007-retrieval-status.md](../../digests/bridges/chat-2609.38007-retrieval-status.md) §7.4 are correct but lose at least one power of $\lvert\beta'-\beta\rvert$ against Theorems 6.7 and 7.1. Their exponent $2\alpha/(1+\alpha)$ is the price of noncommutative generality. They become the right tool only for states with coherences.

### 8.3 Coherent states on histories

**DERIVED (Lemma 8.3).** Fix one origin and let $\mathcal H=\ell^2(\text{complete histories from it})$. For $0<a<L$,
$$\mathcal H=\bigoplus_{\sigma\in K_a}\mathcal H^{<}_\sigma\otimes\mathcal H^{>}_\sigma ,$$
where $\mathcal H^{<}_\sigma$ is spanned by prefixes ending at $\sigma$ and $\mathcal H^{>}_\sigma$ by suffixes starting there. A history through $\sigma$ is the same thing as a (prefix, suffix) pair.
- Note 1's modular Hamiltonian $H_F=\mathrm{diag}F(\mu)$ acts on each block as $H^<_\sigma\otimes1+1\otimes H^>_\sigma$, since $F$ is additive over arrows.
- Its Gibbs state is therefore a classical mixture over the cut face of product states. This is the history-space form of "exactly Markov" (Note 1 §5.3; [arxiv-2609.38007.md](../../digests/bridges/arxiv-2609.38007.md) B10). $\square$

**CONJECTURE (8.4, coherent window Markov property).** Let $H=H_F+V$, $V=\sum_XV_X$, where each $V_X$ is self-adjoint and lies in the span of matrix units $e_{\mu\nu}$ with $\mu,\nu$ differing only in the layers of a depth window $X$. Assume $\lvert X\rvert\le R_0$, $\lVert V_X\rVert\le J$, and that each layer meets at most $d$ windows. Let $\rho=e^{-\beta H}/Z$. Cut depth into $A$ = layers $<a$, a buffer $B$ = layers $a..b$, and $C$ = the rest, and assume that $V$ commutes with the face projections at layers $a$ and $b$. Then, with $C_\beta,c_\beta$ depending only on $\beta,J,R_0,d$ and the largest face count,
$$I_\rho(A:C\mid B)\le C_\beta\exp\big(C_\beta\,g-c_\beta(b-a)\big),\qquad g=\sum_{X\cap A\neq\emptyset\neq X\setminus A}\lVert V_X\rVert ,$$
where the conditional mutual information is defined blockwise through Lemma 8.3. That is legitimate because the commutation hypothesis keeps the cut faces classical registers for $\rho$. Without it, a window term that changes a cut face creates coherences between different summands of Lemma 8.3, no canonical algebra of "layers $\le b$" contains them, and the CMI would first have to be defined (the direct-sum guard).

*What must be true:*
- (i) Yang's cut comparison (Lemma III.1) for the reference $e^{-\beta(H-V_{\rm cut})}$, which factorises across the cut conditionally on the cut face by Lemma 8.3;
- (ii) a Lieb–Robinson bound in depth for $H-V_{\rm cut}$. The local dimension is the largest face count, the hard constraints are the direct-sum structure, and arrows are not reversible, so the light cone must be defined from quiver path length ([local-to-global-unlocks.md](../../local-to-global-unlocks.md) §6);
- (iii) the subalgebra form of §8.1 for the non-factor retained algebra $\bigoplus_\sigma B(\mathcal H^<_\sigma)\otimes1$.

*Consistency:* at $V=0$ the CMI is $0$ (Lemma 8.3). This is the setting in which the "memory as noncommutativity" bridge of [arxiv-2609.38007.md](../../digests/bridges/arxiv-2609.38007.md) B13 would become a statement.

**SPECULATION (8.5, noncommutative uniform local gap; not a precise statement until the "layer-transfer completely positive maps of $\rho$" are defined, e.g. for a finitely correlated structure along depth).** For $\rho$ as in 8.4, consider the KMS-symmetric window Lindbladians with square-flip couplings, built by the Ding–Li–Lin / Chen–Kastoryano–Gilyén construction ([nc-dirichlet-lindblad.md](../../digests/bridges/nc-dirichlet-lindblad.md) §2.3). Suppose the layer-transfer completely positive maps of $\rho$ have finite projective diameter in the noncommutative Hilbert metric. Then window generators of length $w$ have gap $\ge c/w$ uniformly over boundary conditions, for blocks of length depending only on that diameter.

*Status:*
- at $V=0$ the state is diagonal. The diagonal (classical) sector of these generators is a single-flip jump process with filtered Metropolis rates, whose gap follows from Theorem 6.4 only after a comparison of Dirichlet forms with the heat-bath sampler, with a constant depending on $\Delta$ and the face counts. The off-diagonal (dephasing) sector is not covered by Theorem 6.4. Neither step is written out here;
- the noncommutative Hilbert-metric contraction theory (Reeb–Kastoryano–Wolf, cited from memory, not read) would replace Lemma 6.2;
- no $k$-level noncommutative local-to-global theorem was found in the literature searched ([hdx-spectral-independence.md](../../digests/bridges/hdx-spectral-independence.md) §4.5: only the two-level approximate tensorization of Bardet–Capel–Rouzé).

*If true:* a Chen–Rouzé Cor. B.2-type global Markov property for 8.4's states, and Theorem 6.7 for coherent families.

---

## 9. The face law and spectral independence

A state $\omega$ on the face algebra $D_0(K)\cong C(K)$ of a layer (Note 2 §1.1) is a law on $\{0,1\}^V$ supported on faces, with vertex projections $e_u$. Its influence matrix under Lüders pinning, $\Psi_\omega(u,v)=\omega^{e_u}(e_v)-\omega^{1-e_u}(e_v)$, is ALO's ([hdx-spectral-independence.md](../../digests/bridges/hdx-spectral-independence.md) §6.2).

**DERIVED (Proposition 9.1, R13).** Restrict to the non-deterministic vertices and put $D=\mathrm{diag}\,\mathrm{Var}_\omega(e_u)$.
- $\Psi_\omega=D^{-1}\mathrm{Cov}_\omega-I$, so $\eta_0(\omega):=\lambda_{\max}(\Psi_\omega)=\lambda_{\max}(\mathrm{Cor}_\omega)-1$.
- For every real $a$, $\mathrm{Var}_\omega(\sum_ua_ue_u)\le(1+\eta_0)\sum_ua_u^2\mathrm{Var}_\omega(e_u)$. In particular the face dimension $d=\sum_ue_u$ has $\mathrm{Var}_\omega(d)\le(1+\eta_0)\sum_up_u(1-p_u)$.

*Proof.* $\omega(e_v\mid e_u=1)-\omega(e_v\mid e_u=0)=\mathrm{Cov}_{uv}/\mathrm{Var}_u$. Moreover $D^{-1}\mathrm{Cov}$ is similar to $\mathrm{Cor}$, and $a^{\rm T}\mathrm{Cov}\,a\le\lambda_{\max}(\mathrm{Cor})\,a^{\rm T}Da$. $\square$

*Applied to the measurement* ([transfer-spectrum-measurement.md](../../digests/bridges/transfer-spectrum-measurement.md) §3.5; I re-derived the definition used there, $\eta_0=\lambda_{\max}(\mathrm{Cor})-1$):
- the measured $\eta_0\approx1.7$–$8$ says every linear statistic of the barycentre's included-vertex pattern, on the uncertain units the measurement used, has variance at most $2.7$–$9$ times its value under the product state with the same marginals;
- that is all *unpinned* spectral independence gives.

**THEOREM** (ALO Thm 1.3, as quoted in [hdx-spectral-independence.md](../../digests/bridges/hdx-spectral-independence.md) §3.5 and §6.2, re-read in the critique pass). If $\omega$ and its iterated Lüders pinnings have $\eta_i\le\eta_i^*$, the single-vertex walk has spectral gap at least $\frac1n\prod_{i=0}^{n-2}\big(1-\frac{\eta_i^*}{n-i-1}\big)$.

**DERIVED (correction).** A uniform bound $\eta_i\le\eta$ alone does **not** give $n^{-(1+\eta)}$: when $n-i-1\le\eta$ the factors may vanish. Counterexample (§12, A): the law on $\{0,1\}^4$ with mass $\propto1$ on $0000$ and $1111$ and $\propto\varepsilon$ elsewhere has $\max_i\eta_i\to3$, so the claimed bound would be $4^{-4}\approx3.9\cdot10^{-3}$, while its Glauber gap is $6\cdot10^{-7}$ at $\varepsilon=10^{-6}$. ALO's Remark 1.7 states the form $1-1/d^{1+\alpha}$ loosely: it follows from their product formula only while every factor stays positive, and then only up to a constant depending on $\alpha$. Their hardcore bound uses a cap on the deepest pinnings, $\eta_i\le\theta(n-i-1)$ with $\theta=\lambda/(1+\lambda)<1$ (ALO Thm 1.8, Rem 1.10). In general, single-site marginals in $[b,1-b]$ under all pinnings give the cap with $\theta=1-2b$, and under a cap $\theta<1$ the product bound is at least a constant depending on $\eta$ and $\theta$ times $n^{-(1+\eta)}$. **DERIVED (arithmetic):** with $\eta_i\le\min(7,\theta(n-i-1))$ and $n=1024$ the product bound is $10^{-21.0}$, $10^{-27.4}$ and $10^{-34.4}$ for $\theta=0.5$, $0.9$, $0.99$: polynomial in form, vacuous in practice.

So even a confirmed pinned bound at the measured level, with the cap, would make U6 a qualitative statement only. The optimal Chen–Liu–Vigoda route needs bounded degree, which the dense face law lacks.

**ANALOGY (two Markov random fields; do not conflate).**
- The history law is a one-dimensional field over *layers*, with huge spin spaces (faces). Under positivity Theorem 6.4 governs it.
- A layer's face law is a dense, mean-field field over *units*. Theorem 6.4 says nothing about it, and spectral independence is the only certificate available.
- *Where the two meet:* the face law at layer $\ell$ is the layer-$\ell$ marginal of the history law. For a Gibbs history law, the marginals at distant layers decorrelate at the rate of Corollary 6.5. That says nothing about the correlations *within* a layer, which are what $\eta_0$ measures.

---

## 10. The user's intuition, tested in the programme

| claim, in programme terms | verdict | label | § |
|---|---|---|---|
| Gibbs states on histories are Markov random fields | literally a theorem, for every $F$ and $\beta$ (Doob chain). Conversely every positive Markov chain on the quiver is a KMS state for some $F$ | THEOREM (Note 1 §5.3; [hdx-spectral-independence.md](../../digests/bridges/hdx-spectral-independence.md) U7) | 4 |
| the Markov property is an expansion property | **false**: it holds with $\Delta=\infty$ (hard constraints) and at every temperature, while expansion degrades with $\beta\times$zigzag | DERIVED | 4, 6.7 |
| Hammersley–Clifford positivity and an expander gap share an essence | **a theorem here**: a finite projective diameter (positivity quantified) gives a uniform local gap and exponential forgetting with explicit constants, and it is inherited by every conditioning | DERIVED (Lemma 6.2–Cor 6.6; classical in substance: Dobrushin–DGJ path coupling and Birkhoff contraction, §6.2); the two-resolutions reading of [hammersley-clifford.md](../../digests/bridges/hammersley-clifford.md) B3 is ANALOGY | 6 |
| the sufficiency of the barycentre (the NCG/Petz side) and expansion are one essence | **two halves of one class**: the zigzag image of $[F]$ controls the expansion certificates (the projective diameters; the exact gap also feels the mismatch part) and the whole class controls sufficiency. On complete quivers sufficiency implies perfect expansion, not conversely | DERIVED (6.8) | 6.5 |
| the time-averaged detailed-balanced single-Pauli Lindbladian hides an analogue | yes: the window re-routing sampler. KMS = detailed balance; Cesàro limit = the window conditional expectation; single Paulis ↔ square flips. In the diagonal KMS family it is exact and local, so it is trivial there; its substance appears only for coherent states | DERIVED (5.1–5.2) + CONJECTURE (8.4) + SPECULATION (8.5) | 5, 8 |
| a trickle-down theorem for local gaps of KMS samplers is the missing piece | in the commutative positive case **not missing**: the one-dimensional path coupling of Theorem 6.4 (Dobrushin's condition in one dimension) supplies it. Missing for coherent states and for dense face laws | DERIVED (classical: Dobrushin–DGJ) + SPECULATION (8.5) | 6.7, 8.5, 9 |
| expansion of the quiver itself would help | beside the point: what matters is the projective diameter of the transfer kernels. The quiver graph's expansion is the "other expander" of [hdx-spectral-independence.md](../../digests/bridges/hdx-spectral-independence.md) §6.5 | ANALOGY (guard) | 6.7 |
| NCG ties it together | yes at the operator-algebra and measure layer: groupoid KMS cocycles, Petz–Jenčová sufficiency, conditional expectations, Dirichlet forms as derivations, and now the cohomology of the arrow cocycle on the quiver's square complex. No spectral-triple or K-theory statement arises in this angle | KNOWN-LINK + DERIVED | 3, 4, 7, 8 |

---

## 11. Messages to the other prongs (rule C1), and the competition context

**To prong 1 (theory).**
1. **Two corrections to Note 1 §6.** Replace the second bullet by Theorem 3.2, recording $O(\Lambda)=Z_1(\Gamma)/p(B)$ and its sufficient conditions (Corollary 3.3). State which family the first bullet and the chain of equivalences refer to: correct for $P^{\rm gl}$ and for one input face; for $P^w$ with several input faces use Theorem 3.4(b) and Corollary 3.5(ii).
2. Record Proposition 4.1: Markov and barycentre sufficiency are both Jenčová–Petz sufficiency, for different families.
3. Record Lemmas 6.1–6.3, Theorem 6.4, Corollaries 6.5–6.6, Theorem 6.7 and Propositions 6.8–6.9 as Derivations, citing Dobrushin–DGJ and Birkhoff alongside Lemmas 6.2–6.3 and the single-site case of Theorem 6.4, which are classical. The new definitions are the projective diameter of a layer kernel, the zigzag and mismatch parts of $[F]$, and $r$-primitivity.
4. Define $X_\square$ and record Theorem 7.1. Open: Conjecture 7.2 and Speculation 7.3.
5. Lemma 8.3 is the structure of the history space across a depth cut; Conjecture 8.4 and Speculation 8.5 are where Chen–Rouzé and Yang would carry real content.

**To prong 2 (bridge).**
1. Dictionary v1's fitted kernels have many zero transitions, so $\Delta_\ell=\infty$ at a single layer. Report the onset depth $r$ at which multi-step kernels become positive, and their projective diameters $\Delta^{(r)}$. These are the coordinate-free certificates of Theorem 6.4 and Corollary 6.5.
2. For the Markov part (the fitted kernel), report $\max_\square\lvert\delta F(\square)\rvert$, the zigzag norms $Z_\ell$, and the pair-fibre variance against $S^2/(2\,\mathrm{gap})$. A violation is charged first to the estimator (P2.4) and then to the dictionary (P2.1).
3. Memory: report per-layer CMI with the whole past. Their sum is the KL distance to the KMS family (§4), with a Markov surrogate calibrating it at $0$.
4. The face law: unpinned $\eta_0$ yields only Proposition 9.1; the next number to measure is pinned $\eta$.
5. Corrections 1–2 above do not affect E3, which used a single input face.

**To prong 3 (unlocks).**
- **U1:** its conclusion type, forgetting of boundary conditions at a geometric rate, holds in the commutative resolution under ($r$-)positivity with the local certificate $\max_\ell\tanh(\Delta_\ell/4)$ (Corollary 6.5, classical). U1 as stated (link spectra ⇒ transfer gap) is not proved here, and its literal conclusion is trivial for Markov laws (§6.3).
- **U3:** what square sums control is pair sufficiency; the coboundary condition additionally needs $O(\Lambda)=0$. The $\ell^2$ form is proved (Theorem 7.1); the local (Garland) half is open (Speculation 7.3).
- **U6:** unchanged ([hdx-spectral-independence.md](../../digests/bridges/hdx-spectral-independence.md) §6.2), with Proposition 9.1 as the only consequence of the present measurement.
- **New reduction**, in the format of the unlocks note: *global object*: the forgetting of boundary conditions, the window recovery, the sufficiency defect; *local certificate*: the projective diameter of one layer kernel; *explicit bound*: Theorem 6.4, Corollary 6.5, Theorem 6.7; *cost*: block length $e^{\Delta/2}$.
- **New guard:** positivity. A transfer with zeros or signs is outside every statement of §6.

**Competition context ([competition-plan.md](../../competition-plan.md) §§3.1, 6b, 7; stream reports), without dragging the theory.**
- The quantities on which the competition turns are moment and cumulant objects: the D21 interface, old-source content, and the (2,1,1) κ4 slice that [oracle1024](../oracle1024/REPORT.md) and [chain128](../chain128/REPORT.md) identify as the binding error. None of them is a KMS state on histories.
- The old-content transport measured in [transfer-spectrum-measurement.md](../../digests/bridges/transfer-spectrum-measurement.md) is a *signed* sum over histories through open gates, so it lies outside the positivity hypothesis of §6: Theorem 6.4 and Corollary 6.6 make no prediction of a gap there. The measured free-probability law with no gap, and the absence of a few-mode carrier in [old-content](../old-content/REPORT.md), are consistent with that and say nothing against §6.
- **ANALOGY:** the one positive sector, the rectified mean direction, is the one place a Perron mode appears ([transfer-spectrum-measurement.md](../../digests/bridges/transfer-spectrum-measurement.md) B3, B6). *Where it breaks:* that sector is not proved to be a positive cone.
- Under rule P2.1 these negatives are charged to the expectation that the signed old-content transport behaves like a positive (KMS) transfer operator, not to the theory.
- Nothing in this note bears on Lines A–B.

---

## 12. Numerical checks

All checks were run in the session scratchpad. The compact script in the appendix reproduces C1, C2, G2/G4 and G5. Pure numpy. Entries marked "critique pass" were run by the correctness critic with independent code.

| id | statement | result |
|---|---|---|
| C1 | depth-1 counterexample (R1) | pair sufficient: True; coboundary: False; $\dim Z_1=1$, $\dim p(B)=0$, $\dim Z_1(R)=1$ |
| C2 | Lemma 3.1 and Theorem 3.2(b) on random trimmed quivers with parallel arrows | $p(\ker\partial_h)=Z_1$ and $\dim\ker\partial_h/B=\dim Z_1(R)$: 243/243; for pair-sufficient $F$, coboundary ⇔ $\Phi\in B^1(R)$: 243/243 |
| C3 | when Note 1's equivalence holds | $Z_1=p(B)$ in 217/243; $R$ a forest ⇒ holds in all cases; where it fails, generic pair-sufficient $F$ is never a coboundary; waist in every component ⇒ holds: 270/270 |
| C4 | complete layered quivers | squares span $Z_1$ ($\dim$ 3, 6, 12, 17 for four shapes); $\max Z_\ell\le2S$ in all four |
| C5 | Theorem 3.4 | the two examples of §3.4 as stated (the second needs $U(t_1)\neq U(t_2)$, critique pass); criterion (b) right on 263/263 random quivers with random coboundaries (critique pass, independent code: 263/263, and criterion (a) 263/263) |
| G1 | Lemma 6.1 | max TV$/\tanh(D/4)=1.0000$ over 20,000 pairs (attained) |
| G2 | Theorem 6.4, $b=1$ | absolute gap $\ge\kappa_1$ on 47/47 chains with $\kappa_1>0$; median gap$/\kappa_1=1.37$ (critique pass: 238/238, smallest ratio 1.0001) |
| G3 | Theorem 6.4, blocks | $\Delta\in[2.43,3.32]$ (beyond $2\log3=2.197$): smallest $b$ with $\kappa_b>0$ is 2–3, and the absolute gap $\ge\kappa_b$ in 10/10 |
| G4 | Theorem 6.7 | Var$\,\le S^2/(2\,$gap$)$ always (max ratio 0.175); KL defect below its bound in 8/8 |
| G5 | Corollary 6.5 | max TV$/\prod\tanh(\Delta/4)=0.62$ over 200 chains |
| H | Theorem 7.1 and Conjecture 7.2 | $\mu_1=4.000$ ($n=2$) and $9.000$ ($n=3$) for $L=2..8$; oscillation bound ratio $\le0.453$ (critique pass: $\mu_1=16.000$ for $n=4$, $L\le5$, and $25.000$ for $n=5$, $L\le4$) |
| P | Corollary 6.6, Penrose | $\tanh(\log2/4)=0.1716\ge\varphi^{-4}=0.1459$ |
| M | Proposition 6.8(b): what the gap sees | complete quiver of width 3 and depth 5, zigzag part fixed (every $\Delta_\ell$ unchanged), mismatch fields scaled by 0, 1, 3, 6: fibre gap 0.133, 0.133, 0.150, 0.250 (critique pass) |
| A | §9, ALO product form | $\{0,1\}^4$, mass on $0000,1111$ and $\varepsilon$ elsewhere: $\max_i\eta_i=2.78$, $3.00$, $3.00$ and Glauber gap $5.9\cdot10^{-3}$, $6.0\cdot10^{-5}$, $6.0\cdot10^{-7}$ at $\varepsilon=10^{-2},10^{-4},10^{-6}$, against the claimed $n^{-(1+\eta)}\approx3.9\cdot10^{-3}$ (critique pass) |

These are checks of algebra and of the stated inequalities on small instances, not experiments on any network.

---

## Appendix: compact reproduction script

```python
import itertools, numpy as np
rng = np.random.default_rng(0); rk = lambda M: np.linalg.matrix_rank(M, 1e-9) if M.size else 0

def quiver(sizes, p, mult=1):                      # random trimmed layered quiver; arrow = (l, i, j, copy)
    E = [(l, i, j, c) for l in range(len(sizes)-1) for i in range(sizes[l]) for j in range(sizes[l+1])
         if rng.random() < p for c in range(rng.integers(1, mult+1))]
    L = len(sizes)-1
    while True:
        f = {(0, i) for i in range(sizes[0])}; b = {(L, i) for i in range(sizes[L])}
        for l in range(L): f |= {(l+1, e[2]) for e in E if e[0] == l and (l, e[1]) in f}
        for l in range(L-1, -1, -1): b |= {(l, e[1]) for e in E if e[0] == l and (l+1, e[2]) in b}
        E2 = [e for e in E if (e[0], e[1]) in f & b and (e[0]+1, e[2]) in f & b]
        if len(E2) == len(E): return sorted(f & b), E
        E = E2

def histories(E, L):
    H = [[e] for e in E if e[0] == 0]
    for l in range(1, L): H = [h+[e] for h in H for e in E if e[0] == l and e[1] == h[-1][2]]
    return H

def structure(V, E, L):                             # d (coboundary), p (history -> chain), bigons, R
    vi = {v: k for k, v in enumerate(V)}; ei = {e: k for k, e in enumerate(E)}; H = histories(E, L)
    d = np.zeros((len(E), len(V)))
    for e in E: d[ei[e], vi[(e[0]+1, e[2])]] += 1; d[ei[e], vi[(e[0], e[1])]] -= 1
    P = np.zeros((len(E), len(H)))
    for k, h in enumerate(H):
        for e in h: P[ei[e], k] += 1
    fib = {}
    for k, h in enumerate(H): fib.setdefault((h[0][1], h[-1][2]), []).append(k)
    Bc = [np.eye(len(H))[k] - np.eye(len(H))[ks[0]] for ks in fib.values() for k in ks[1:]]
    PB = P @ np.array(Bc).T if Bc else np.zeros((len(E), 0))
    S = sorted({a for a, _ in fib}); T = sorted({c for _, c in fib}); M = np.zeros((len(fib), len(S)+len(T)))
    for k, (a, c) in enumerate(fib): M[k, S.index(a)] = -1; M[k, len(S)+T.index(c)] = 1
    return d, P, PB, fib, M

def in_range(A, x): return np.linalg.norm(A @ np.linalg.lstsq(A, x, rcond=None)[0] - x) < 1e-8

print("C1 depth 1, K22: F = 1 on one arrow")
V, E = quiver([2, 2], 1.0); d, P, PB, fib, M = structure(V, E, 1); F = np.eye(len(E))[0]
print("   pair-sufficient", all(np.ptp((P.T @ F)[k]) < 1e-12 for k in fib.values()), " coboundary", in_range(d, F))

print("C2 cycle-space theorem on random quivers")
ok = holds = forest_ok = n = 0
for _ in range(200):
    L = int(rng.integers(1, 5)); V, E = quiver(list(rng.integers(1, 5, L+1)), rng.uniform(.3, .9), 2)
    if not E: continue
    d, P, PB, fib, M = structure(V, E, L); n += 1
    dimZ1 = len(E) - rk(d); dimpB = rk(PB); dimZ1R = M.shape[0] - rk(M)
    X = rng.normal(size=len(E))
    if PB.shape[1]: U = np.linalg.svd(PB, full_matrices=False)[0][:, :rk(PB)]; X -= U @ (U.T @ X)
    Phi = np.array([(P.T @ X)[ks[0]] for ks in fib.values()])
    ok += in_range(d, X) == in_range(M, Phi); holds += dimZ1 == dimpB; forest_ok += (dimZ1R > 0) or (dimZ1 == dimpB)
print(f"   {n} quivers: [coboundary <=> Phi coboundary on R] for pair-sufficient F: {ok}/{n};"
      f" Z1 = p(B): {holds}/{n}; R forest => Z1 = p(B): {forest_ok}/{n}")

def cross(K):
    lK = np.log(K); return max(abs(lK[i, j]+lK[k, l]-lK[i, l]-lK[k, j])
                               for i, k in itertools.combinations(range(K.shape[0]), 2)
                               for j, l in itertools.combinations(range(K.shape[1]), 2))

def heatbath(sizes, K, a, c):                       # single-site heat-bath on window 1..m-1, ends pinned
    m = len(sizes)-1; st = list(itertools.product(*[range(sizes[k]) for k in range(1, m)])); ix = {s: i for i, s in enumerate(st)}
    wt = lambda x: np.prod([K[k][x[k], x[k+1]] for k in range(m)])
    pi = np.array([wt((a,)+s+(c,)) for s in st]); pi /= pi.sum(); Q = np.zeros((len(st), len(st)))
    for s in st:
        for j in range(1, m):
            ts = [tuple(s[:j-1]) + (v,) + tuple(s[j:]) for v in range(sizes[j])]
            ws = np.array([wt((a,)+t+(c,)) for t in ts])
            for t, w in zip(ts, ws / ws.sum()): Q[ix[s], ix[t]] += w / (m-1)
    D = np.sqrt(pi); ev = np.sort(np.linalg.eigvalsh((D[:, None]*Q/D[None, :] + (D[:, None]*Q/D[None, :]).T)/2))
    return st, pi, 1-ev[-2], 1-max(abs(ev[0]), ev[-2])

print("G2/G4 positive chains: abs gap >= kappa = min_j(1 - t_{j-1} - t_j)/w ; Var(F|ends) <= S^2/(2 gap)")
good = tot = 0; ratios = []
for _ in range(40):
    m = int(rng.integers(3, 6)); sizes = list(rng.integers(2, 4, m+1)); beta = rng.uniform(.1, 1.)
    G = [rng.normal(size=(sizes[k], sizes[k+1])) for k in range(m)]; K = [np.exp(-beta*g) for g in G]
    t = [np.tanh(cross(k)/4) for k in K]; w = m-1
    kap = min(1-(t[j-1] if j > 1 else 0)-(t[j] if j < w else 0) for j in range(1, w+1))/w
    st, pi, gap, agap = heatbath(sizes, K, 0, 0)
    Fv = np.array([sum(G[k][x[k], x[k+1]] for k in range(m)) for x in [(0,)+s+(0,) for s in st]])
    S = max(abs(sum(G[k][y[k], y[k+1]] for k in range(m)) - Fv[i]) for i, s in enumerate(st)
            for j in range(1, w+1) for v in range(sizes[j]) for y in [list((0,)+s+(0,))] if not y.__setitem__(j, v))
    var = pi @ Fv**2 - (pi @ Fv)**2; ratios.append(var / (S**2/(2*gap)))
    if kap > 0: tot += 1; good += agap >= kap - 1e-12
print(f"   gap >= kappa in {good}/{tot} chains with kappa > 0; max Var/(S^2/2gap) = {max(ratios):.3f}")

print("G5 forgetting: TV(x_t | v) vs (x_t | v') <= prod_{k>=t} tanh(Delta_k/4)")
worst = 0
for _ in range(200):
    m = int(rng.integers(3, 7)); sizes = list(rng.integers(2, 4, m+1)); K = [np.exp(-rng.uniform(.1, 2)*rng.normal(size=(sizes[k], sizes[k+1]))) for k in range(m)]
    w0 = rng.uniform(.2, 1, sizes[0]); vs = [rng.uniform(.01, 1, sizes[-1]) for _ in range(2)]; t = int(rng.integers(0, m))
    def marg(v):
        f = w0.copy(); h = v.copy()
        for k in range(t): f = f @ K[k]
        for k in range(m-1, t-1, -1): h = K[k] @ h
        return f*h/(f @ h)
    worst = max(worst, .5*np.abs(marg(vs[0])-marg(vs[1])).sum() / np.prod([np.tanh(cross(K[k])/4) for k in range(t, m)]))
print(f"   max ratio {worst:.3f}")
```

Output of this script (seed 0): C1 pair-sufficient True, coboundary False; C2 165/165 quivers satisfy the cycle-space theorem, $Z_1=p(B)$ in 148/165, forest ⇒ equivalence 165/165; G2/G4 gap $\ge\kappa$ in 35/35, max Var$/(S^2/2\,$gap$)=0.175$; G5 max ratio 0.591. The counts differ from the table in §12, which comes from larger scratch runs with other seeds; the conclusions are the same.

---

## Sources

Read for this note: the eight bridge digests and five programme notes listed in the header, all in full; [competition-plan.md](../../competition-plan.md) §§0, 3, 3.1, 6b, 7; the stream reports in [streams/](../); and `notes/checks/graph_algebra_checks.py` (to fix which Gibbs family Note 1 §9 tested).

No primary paper was opened for the first version of this note; the critique pass opened one (below). Every published theorem used is quoted from the digest section cited at its use, where its source was read:
- Yang Lemma III.4: [arxiv-2609.38007.md](../../digests/bridges/arxiv-2609.38007.md) §5;
- Chen–Rouzé: [arxiv-2504.02208.md](../../digests/bridges/arxiv-2504.02208.md);
- Hammersley–Clifford and its positivity analyses: [hammersley-clifford.md](../../digests/bridges/hammersley-clifford.md);
- ALO / Alev–Lau / Oppenheim: [hdx-spectral-independence.md](../../digests/bridges/hdx-spectral-independence.md);
- Kastoryano–Brandão, Mossel–Sly, HLW: [expanders.md](../../digests/bridges/expanders.md);
- Vernooij–Wirth, Carlen–Maas, Alicki: [nc-dirichlet-lindblad.md](../../digests/bridges/nc-dirichlet-lindblad.md);
- Renault/Neshveyev and Jenčová–Petz: Note 1 §§5.2, 6.

Opened in the critique pass: N. Anari, K. Liu, S. Oveis Gharan, arXiv:2001.00303v3 (via alphaXiv): Defs 1.1–1.2, Thm 1.3, Thms 1.5–1.6, Rem 1.7, Thm 1.8, Lemma 1.12.

From memory, flagged where used and not relied on for any DERIVED statement: Birkhoff–Hopf contraction in Hilbert's projective metric (§6.3, §6.7); the bound of the Dobrushin coefficient by Birkhoff's coefficient (§6.2; Lemma 6.2 proves the form used); the contraction-to-gap theorem as cited from memory in [expanders.md](../../digests/bridges/expanders.md) §11.1 (§6.2; Lemma 6.3 proves it); and the noncommutative Hilbert-metric theory of Reeb–Kastoryano–Wolf (Speculation 8.5).

Not retrieved: the ChatGPT conversation behind arXiv:2609.38007; the local `expander_survey.pdf`.

---

## Critique log (correctness pass, 2026-10-01)

Every DERIVED proof was re-derived, and the finite statements were re-checked with independent code: Lemma 3.1, Theorem 3.2, Corollary 3.3 (333 random quivers), Theorem 3.4 (263/263), Proposition 4.1, Lemmas 6.1–6.3, Theorem 6.4 (238 random and adversarial chains, 28 block cases), Corollaries 6.5–6.6, Theorem 6.7(a) for both families, Proposition 6.8, Theorem 7.1, Conjecture 7.2 (now also $n=4,5$), Remark 8.2 (defect $\propto\Delta\beta^2$, bound $\propto\Delta\beta$) and Proposition 9.1. They hold. Corrected:

1. **§9, the ALO bound.** "$\eta_i\le\eta$ for every pinning ⇒ gap $\ge n^{-(1+\eta)}$" is false: the product factors can vanish at the deepest pinnings (counterexample, §12 A). Replaced by ALO's product form plus the cap $\eta_i\le\theta(n-i-1)$, with the arithmetic redone ($10^{-21}$–$10^{-34}$ at $\eta=7$, $n=1024$).
2. **§6.7, positivity.** "A safe symbol makes the support a product" is false; a safe symbol only gives closure under switching to the vacuum ([hammersley-clifford.md](../../digests/bridges/hammersley-clifford.md) §1.2). Besag's positivity gives the product.
3. **§6.5(b), what mixing sees.** The projective diameters see only the zigzag image of $[F]$, but the exact gap also depends on the mismatch part (single-site fields; §12 M). The text now says "certificates".
4. **§5 table, KMS and GNS.** "KMS = GNS because the state is diagonal" is not a valid reason. The reason is that the re-routing jumps are modular eigenvectors.
5. **§3.4, example 2** needs $U(t_1)\neq U(t_2)$. **§6.5(d)**: endpoint sufficiency for $P^{\rm gl}$ is stronger than $S=0$ when there are several input faces. **§6.7**: the low-temperature statement ("$c$ deteriorates") is replaced by what Theorem 6.4 gives. **§6.6(ii)** and **§6.4(b)**: missing hypotheses added. **§8.1**: $0<\alpha\le1$. **§9**: "fluctuates" made "variance".

Labels changed:

- **R7, §6.3, §11 (U1).** "U1 becomes a theorem" was inflated. The forgetting bound is DERIVED (and classical). Its identification with U1 is ANALOGY: U1's hypothesis is unused, and U1's literal conclusion is trivial for Markov laws.
- **§6.2, R6, §10.** Theorem 6.4 for $b=1$ and Lemmas 6.2–6.3 are classical (Dobrushin–DGJ path coupling, Birkhoff's coefficient, contraction ⇒ gap). They are now credited as a KNOWN-LINK; the proofs stay.
- **§6.7, §10 row 3.** The citation of [hammersley-clifford.md](../../digests/bridges/hammersley-clifford.md) B3 changed from KNOWN-LINK to ANALOGY, because B3 is that digest's own derivation, not a published link.
- **§5.2(d).** Changed from THEOREM (cited) to DERIVED, with the one-line proof.
- **7.3 and 8.5.** Changed from CONJECTURE to SPECULATION: neither is a precise statement yet. Speculation 8.5's "at $V=0$ this is Theorem 6.4" is replaced by what is actually covered.
- **Conjecture 8.4.** Now assumes that $V$ commutes with the cut-face projections, so that its CMI is defined, and names what its constants depend on.
