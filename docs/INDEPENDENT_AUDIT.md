# Independent adversarial mathematical audit

Date: 2026-11-16 (Asia/Shanghai). Scope: the actual v1.1 `main.tex`, not the earlier worker's conclusions. This is a bounded independent AI audit, not human refereeing, a formal proof-assistant certificate, or a guarantee of correctness. Nothing was uploaded, published, or sent externally.

## Outcome

No counterexample, circular dependence, or substantive gap in the characteristic-zero proof of Theorem 1.1 was found in the checks below. I retained the main proof and the finite-generator Koszul argument. I did find one genuine missing range qualification in Remark 10.5: with V=k^N fixed, the content-[n] component of Λ^n(V*) has dimension one only when n≤N, not for n>N. This is fixed by stating the range and taking N≥n degree by degree. It does not affect the theorem.

Several other changes prevent overinterpretation rather than repair a false main result:
- The DST Koszul statement is explicitly about the multilinear species, not an unqualified assertion about the ordinary algebra U(k).
- The Segre-product argument now spells out its degree-2 kernel, quadratic-presentation hypothesis, finite-generation hypotheses, and ordinary (not super) tensor multiplication.
- The dictionary now specifies singleton generators, absence of empty-set generators, the formally adjoined arity-zero unit, and products only on disjoint label sets. In particular no geometric L_0 or generator m_0 is assumed.
- “Not contained in LM/DST” became the narrower “not explicitly stated.” No claim that earlier work cannot imply the theorem is made.
- LM/DST's coefficient conventions are stated at the point of use, without asserting that all of LM's geometry required characteristic zero.

## 1. Polarisation (Proposition 3.1)

I checked functoriality in V rather than assuming that relabelling indices alone was enough. The maps on F(vw) are induced by Sym² of the linear map; the cubic exchange relation is trilinear in its inputs, so the substitution maps really preserve the defining ideal. The symmetry of [D(v),E(w)] makes the map from Sym²(V) well-defined. Jacobi and [E,E]=0 kill the cubic relation.

For a multihomogeneous kernel element x of multiplicities m_j, substitute v_i -> Σ_{p∈P_j}t_p v_p and fold each P_j back onto its original v_i. Folding the coefficient of t_1...t_n gives (∏m_j!)x. The coefficient is nonzero in characteristic zero. Each homogeneous component after substitution stays in the kernel. Overlap between the *names* of the old and new indices causes no problem: these are linear maps, not permutations, and the composite on the source support is exactly the required diagonal scalar map. This proves reduction to multilinear components without a representation-theoretic semisimplicity assumption.

## 2. Presentation of k (Proposition 4.1)

I checked both inverse maps. Commuting ad D_i permit rearrangement; the mixed symmetry [D_i,E_j]=[D_j,E_i] exchanges the distinguished E index with any D index. The derivations ∂_i of the free Lie algebra commute on generators and hence everywhere. The ideal of all ∂_S[G_i,G_j] is stable under them. Thus the semidirect product exists before making any injectivity assertion.

The map to g kills every displayed relator because it is ad(D_S)[E_i,E_j]. The reverse map on D_i,E_i satisfies all the presentation relations and sends each iterated D-bracket back to G_T. The composites fix the generating sets. Consequently no unstated relation of the ideal k has been discarded. The size-0 and size-1 derivative relators are exactly [E_i,E_j] and [E_i,F_sj]-[E_j,F_is]. This also checks the presentation of l and the asserted kernel (F^sm∩I)/I_l.

## 3. Excess filtration and Lemma A

For |S|=m, a split with |P|=a has weight max(a−1,0)+max(m−a−1,0). If m≥2, the endpoint terms have weight m−1 and every interior term m−2. For i≠j the endpoint brackets cannot cancel in the free Lie algebra: their unordered pairs of letters are distinct. For i=j the whole relator and its proposed leading expression are zero. Cases m=0,1 are handled separately and exactly.

The filtration is increasing by highest weight; this direction is important. The leading ideal of all of I contains the ideal J of the displayed leading relators, not conversely a priori. Hence k^top surjects onto gr(k), giving dim k≤dim k^top. Each fixed content is finite-dimensional because all generator contents are nonempty and bracket length is bounded by content size. Thus equality of dimensions forces J_mu=(gr I)_mu. Nonnegative weights then give J_0=I_l: a nonzero homogeneous iterated adjoint of a relator can have weight zero only if the relator and every adjoint letter have weight zero. An element in I∩F_0 is its own leading component and therefore lies in I_l once that dimension equality holds.

This checks precisely the direction of Lemma A. It does not assume the injectivity it is meant to prove. There is no use of a low-degree flatness experiment in this argument.

## 4. Lemma B and the corollaries

The exchange relation is [E_j,G_T]=[E_t,G_{T−t+j}]. All ad E commute. For two splittings W=C+T with |T|=s fixed, a multiset element in T'−T is available in C, and exchanging it for an element in T−T' strictly decreases the multiset discrepancy. This includes repeated-index cases, not just set splittings.

For the proposed inverse map G_T -> Y_{T,|T|}, both sides of every exchange relation map to the same generator Y_{T+j,|T|}. This produces a genuine homomorphism into the explicitly constructed free-Lie semidirect product. The reverse map on Y_{W,s} is well-defined by the exchange argument and intertwines the E derivations. The two composites fix generators, which rules out hidden relations among the Y's; it is not merely a surjective normal-form argument.

The excess-zero part selects precisely s=2. A free Lie subalgebra generated by a subset of these free generators is free (inclusion admits the retraction killing the complementary generators). Once θ is injective this proves both the multiset-indexed structure statement about h and freeness of all [D_i,E_j] for i≤j, including diagonal terms.

## 5. PBW reduction: termination and completeness

Fixed finite-length lexicographic order over the well-ordered alphabet D_1<E_1<D_2<E_2<... is well-founded. Degree-lex therefore is a compatible well-founded monoid order even with infinitely many indices. Every right-hand monomial is strictly smaller. An alternative purely finite-content observation is that each reduction preserves the multiset of letters and length, so only finitely many candidate words can occur for any starting polynomial.

No rule has a reducible word properly included in another leading word: all leading words have length two and are pairwise distinct. Two overlapping leading pairs require strictly decreasing indices k>j>i. The forbidden adjacent type DE leaves exactly DDD, EDD, EED, EEE. There are no additional repeated-index overlaps. I re-expanded the printed EDD and EED routes; their normal polynomials agree with the manuscript. Disjoint ambiguities are covered by the standard reduction theorem and need not be listed as critical overlaps.

**Bergman source-access limit.** I looked for earlier copies in the project, workspace, temporary/reference files and previous manuscript attachments; no usable local Bergman original was located. The legitimate public academic-course copy http://jthartwig.net/F20-699/papers/Bergman-1978.pdf was located by search. HTTPS retrieval failed with curl TLS error 35; HTTP retrieval returned 500. No PDF was recovered. Search indexing returned the original theorem's hypotheses and parts (a),(a'),(b), but not the full part (c) in the excerpt. Secondary public formulations agree with the quoted normal-basis consequence. I did read the author's corrections page https://math.berkeley.edu/~gbergman/papers/updates/diamond.html; it does not change the hypotheses/conclusion being used here. I therefore independently checked the *application* of the standard Diamond Lemma but do not claim to have inspected the full original Theorem 1.2 or its proof in this audit. This source-access limitation is not a counterexample to the algebraic argument. PBW itself is used as a standard theorem, not reproved here.

## 6. Counting and all-degree Hilbert-series matching

The chains have strictly decreasing indices and type E...ED...D; exactly one chain occurs for each number of D's. The signed pair involution works also when c or w is empty: case II has c nonempty; case I creates a chain and leaves a normal suffix; case II puts a nonbad pair at the beginning of the normal word. The two cases reverse each other. Hence the inverse-series identity holds in every degree, independently of computation.

The equality Y=H(x)−e^{xy}−e^{xz}+1 and multiplication by y−z give H(−x)=e^{−x(y+z)}(1−Y). Cancelling y−z is valid in Q[y,z][[x]] although *dividing* by it as a unit would not be. The remaining denominators 1−Y and H(−x) are units, so logarithms are legitimate. PBW gives H_g=log H_U(g); the previously proved semidirect presentation subtracts xy to give H_k. Lemma B independently gives the identical H_k^top. Only then is Lemma A used, and then polarisation. I found no circular use of the main theorem in this dimension chain.

For the Eulerian corollary, differentiating −log Q gives y(A(x,y)−1), so the n≥2 coefficient is y A_{n−1}(y). Thus g and k have A(n−1,a−1), NOT A(n,a). By contrast U(k) has A(n,a). This shift agrees with the independent LM dictionary. The n=1 components are correctly excluded from that corollary. The Carlitz counting proof and zero boundary Eulerian coefficients are consistent.

## 7. Koszul Remark 10.5: checked, with the finite-N qualification fixed

I read the downloaded BF text (`../refs/bf.txt`), including §1 conventions, 1.2, 1.14, 1.16, 1.17, Proposition 3(b), and Theorem 4(b) with its short proof. These are results for ordinary, not necessarily commutative, finitely generated connected degree-one-generated algebras with finitely generated homogeneous relation ideals. Exterior algebras are allowed. The finite N restriction supplies the finite hypotheses; the argument is not silently applying BF to infinitely many generators.

Rearranging (V⊗U)^⊗2 gives the relation space Λ²V⊗Sym²U, of dimension 3 binom(N,2). I checked the cross-relator sign explicitly: [D_i,E_j]−[D_j,E_i] becomes (v_iv_j−v_jv_i)⊗(de+ed). Its ordinary tensor-pairing annihilator is Sym²V*⊗(U*)^⊗2 + (V*)^⊗2⊗Λ²U*. This is exactly the kernel in degree 2 of the map to Λ(V*)∘Sym(U*). The map is surjective on generators, and BF Proposition 3(b) ensures that the Segre product is quadratically presented, so checking the degree-2 kernel suffices in all degrees. BF 1.17(b) includes polynomial rings with zero monomial ideal; 1.16(10) gives Koszulness of the exterior factor and duality invariance; Theorem 4(b) gives the Segre closure. Thus the Koszul claim about U(g_N) has a valid non-circular proof. The numerical identity alone would not prove Koszulness, and is not being used to do so.

With N fixed, Λ^n(V*)_[n] is 1-dimensional for n≤N and zero for n>N; the old wording omitted this range. The repaired explanation explicitly takes N≥n. There is no claim of an ordinary-algebra Koszul transfer from DST's species assertion.

## 8. LM/DST identification and limits of prior-work conclusions

Read directly in downloaded sources, not just the revision worker's notes:
- LM (`../refs/math_0001003.txt`): Theorem 2.2, Theorem 2.3 and its point-counting calculation, the beginning of the toric description in §2.4, and §§3.1–3.6 (especially (3.5),(3.6),(3.17),(3.19)).
- DST (`../refs/2112.06015.txt`): §§2.1–2.4, including the Cauchy-product convention, weight grading for Koszulness, Propositions 2.4.2 and 2.4.3 and their printed arguments.
- BF as specified above. The other downloaded Norledge paper is not used by the manuscript and was not independently re-audited.

The content-[n] component of the associative ideal of relators is spanned by u r v of disjoint squarefree contents; positivity of contents prevents a repeated-index relator from contributing there. Expansion of r is literally Σ_{i∈I,j∈J}(G_I G_J−G_J G_I), one copy per ordered two-partition. This agrees with DST. The generators, not the entire homology representations, are symmetric: each m_t is a fundamental class with trivial S_t action. Their degrees 2t−2 are even; multiplication of homology classes is the boundary-inclusion product, with no super signs or Poincaré-duality regrading. D-degree n−p is half the homological degree. The E-degree is p, the number of generators. In particular m_1 must not be omitted.

LM begins its homology algebra at n=1 without a unit; DST uses monoids in species. Adjoining the formal unit at arity zero reconciles the EGF's constant 1. There is no need to invent L_0. For n=2 the relation is commuting singleton generators, giving the two degree pieces of H_*(P¹). The identification is a species identification across disjoint sets; an ordinary algebra identification including repeated labels would require the additional Schur/coinvariant construction. The revised manuscript does not claim it.

LM's EGF yields A_n for U(k); multiplying by e^{xy} and taking the PBW logarithm recovers the manuscript's g dimensions. This alternative route uses the presentation of k and Lemma B, not θ-injectivity. Thus its dependency graph is also non-circular.

The report by Laubie is credited for the ordering, but DST actually prints the ordering and its normal-monomial counting argument. The unpublished report is not required to use or inspect Proposition 2.4.2. I did not read the journal/paywalled version or reprove every geometric theorem/citation of LM, and do not claim such coverage. The paper correctly pins its numbering to arXiv versions. Novelty remains “未检索到”; this audit is not a new comprehensive literature search.

## 9. Computational checks actually done in this audit

Re-ran, successfully: `ambiguities.py`, `pbw_check.py`, `prior_work_check.py`, `egf_check.py`. The exact overlap normal forms agree; cubic ranks are 104=104 for U(g) and 212=212 for its dual; U(k) has the Eulerian rows through n=6; the stated closed-series simplifications vanish. Captured output: `audit/out_rerun_checks.txt`.

I also wrote a small independent diagnostic, `audit/audit_repeated_words.py`, which imports none of the author's verification modules. It checks every possible first reduction against leftmost normal reduction on all words in D_1,E_1,D_2,E_2,D_3,E_3 of length ≤5, including repeated indices (1+6+36+216+1296+7776 words). All passed: `audit/out_repeated_words.txt`. This is only a finite counterexample search, not an all-degree proof.

Read `kproj.py` and its relevant coordinate/bracketing routines in `kmodel.py`. The projection recursion is consistent: bracketing by a small letter preserves the split by presence of a big generator, while bracketing by a big letter always lands in positive excess. The anchored multilinear free-Lie coordinates use an integral left-normed basis; ranks can therefore be compared modulo p and over Q. The modular-to-rational inequality has the correct direction. I did NOT rerun n=8 rational or n=9 modular computations in this bounded audit, and do not independently certify the previously recorded runtime/memory figures. The all-degree proof does not need those runs.

## 10. Files, build and remaining limits

The pre-audit v1.1 files are preserved as `main_v1.1_before_audit.tex`, `paper_v1.1_before_audit.pdf`, `CHANGES_v1.1_before_audit.md`; v1.0 originals remain. `audit/manuscript_audit.patch` records the exact modifications. The changes do not alter any manuscript verification algorithms or their claimed results, so `code.zip` is unchanged; audit-only diagnostics are separately saved under `audit/`.

The final source was compiled twice with Tectonic. Both final runs succeeded with no undefined-reference/citation or overfull/underfull warnings in their output. Build output: `audit/out_build.txt`. This is not a visual examination of every rendered page.

Remaining limits: no human expert review, no machine-checked formal proof, original Bergman PDF not recovered, cited prior-work proofs not recursively re-proved in full, no all-source novelty guarantee, and large previous computational runs not repeated here. No guarantee of infallibility is made.

Final PDF: `/workspace/algebra-open-problems/tao-lie/paper/paper.pdf`
Pages: 17.
SHA-256: `1bc4dcac91f28b3152327e65d9e8c99ec1f4ba3fd4576d273955ba3e536699b7`.
Source SHA-256: `edba196f0309fc9b1e84ea86f7dc7a5e067980e4d730bdbdf6c5d97de662b106`.
Unchanged `code.zip` SHA-256: `46240b46d0cca52ccc8ca5890576fc7250eac2d92a659657497f44c17de3bd86`.
